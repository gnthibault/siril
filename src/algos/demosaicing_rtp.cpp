#include "librtprocess.h"
#include "core/siril.h"
#include "algos/demosaicing.h"
#ifdef __cplusplus
extern "C" {
#endif
#include "core/proto.h"
#include "gui/progress_and_log.h"
#ifdef __cplusplus
}
#endif

/* 0 corresponds to red, 1 corresponds to green channel one,
 * 2 corresponds to blue, and 1 corresponds to green channel two */
static void pattern_to_cfarray(sensor_pattern pattern, unsigned int cfarray[2][2]) {
	switch (pattern) {
		case BAYER_FILTER_RGGB:
			cfarray[0][0] = 0; cfarray[0][1] = 1;
			cfarray[1][0] = 1; cfarray[1][1] = 2;
			break;
		case BAYER_FILTER_BGGR:
			cfarray[0][0] = 2; cfarray[0][1] = 1;
			cfarray[1][0] = 1; cfarray[1][1] = 0;
			break;
		case BAYER_FILTER_GBRG:
			cfarray[0][0] = 1; cfarray[0][1] = 2;
			cfarray[1][0] = 0; cfarray[1][1] = 1;
			break;
		case BAYER_FILTER_GRBG:
			cfarray[0][0] = 1; cfarray[0][1] = 0;
			cfarray[1][0] = 2; cfarray[1][1] = 1;
			break;
		case XTRANS_FILTER:
			// take a deep breath
		default:
			break;
	}
}

/* 0 corresponds to red, 1 corresponds to green channel one,
 * 2 corresponds to blue, and 3 corresponds to green channel two */
static void pattern_to_cfarray2(sensor_pattern pattern, unsigned int cfarray[2][2]) {
	switch (pattern) {
		case BAYER_FILTER_RGGB:
			cfarray[0][0] = 0; cfarray[0][1] = 1;
			cfarray[1][0] = 3; cfarray[1][1] = 2;
			break;
		case BAYER_FILTER_BGGR:
			cfarray[0][0] = 2; cfarray[0][1] = 1;
			cfarray[1][0] = 3; cfarray[1][1] = 0;
			break;
		case BAYER_FILTER_GBRG:
			cfarray[0][0] = 1; cfarray[0][1] = 2;
			cfarray[1][0] = 0; cfarray[1][1] = 3;
			break;
		case BAYER_FILTER_GRBG:
			cfarray[0][0] = 1; cfarray[0][1] = 0;
			cfarray[1][0] = 2; cfarray[1][1] = 3;
			break;
		case XTRANS_FILTER:
			// take a deep breath
		default:
			break;
	}
}

static bool progress(double p) {
	// p is [0, 1] progress of the debayer process
	return true;
}

WORD *debayer_buffer_new_ushort(WORD *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, unsigned int xtrans[6][6]) {

	// super-pixel is handled by siril code, not librtprocess
	if (interpolation == BAYER_SUPER_PIXEL)
		return debayer_buffer_superpixel_ushort(buf, width, height, pattern);

	unsigned int cfarray[2][2];
	float rgb_cam[3][4] = { 1.0f };	// our white balance: we don't care
	int i, rx = *width, ry = *height;
	long j, nbpixels = rx * ry;
	long n = nbpixels * 3;
	// 1. convert input data to float (memory size: 2 times original)
	float **rawdata = (float **)malloc(ry * sizeof(float *));
	if (!rawdata) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	rawdata[0] = (float *)malloc(nbpixels * sizeof(float));
	if (!rawdata[0]) {
		PRINT_ALLOC_ERR;
		free(rawdata);
		return NULL;
	}
	// TODO: vectorize!
	for (j = 0; j < nbpixels; j++)
		rawdata[0][j] = (float)buf[j];

	for (i=1; i<ry; i++)
		rawdata[i] = rawdata[i - 1] + rx;

	// 2. allocate the demosaiced image buffer (memory size: 6 times original)
	float *newdata = (float *)malloc(n * sizeof(float));
	if (!newdata) {
		free(rawdata);
		PRINT_ALLOC_ERR;
		return NULL;
	}

	float **red = (float **)malloc(ry * sizeof(float *));
	red[0] = newdata;
	for (i=1; i<ry; i++)
		red[i] = red[i - 1] + rx;

	float **green = (float **)malloc(ry * sizeof(float *));
	green[0] = red[0] + nbpixels;
	for (i=1; i<ry; i++)
		green[i] = green[i - 1] + rx;

	float **blue = (float **)malloc(ry * sizeof(float *));
	blue[0] = green[0] + nbpixels;
	for (i=1; i<ry; i++)
		blue[i] = blue[i - 1] + rx;

	// 3. process
	siril_debug_print("calling librtprocess ushort (%d)\n", interpolation);
	rpError retval;
	switch (interpolation) {
		case BAYER_VNG:
			pattern_to_cfarray2(com.debayer.bayer_pattern, cfarray);
			retval = vng4_demosaic(rx, ry, rawdata, red, green, blue, cfarray, progress);
			break;
		case BAYER_BILINEAR:
//		case BAYER_NEARESTNEIGHBOR:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			/* bayerfast: This demosaicer is not intended for final
			 * output, only for fast preview. */
			retval = bayerfast_demosaic(rx, ry, rawdata, red, green, blue, cfarray, progress, 1.0);
			break;
		default:
		case BAYER_RCD:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = rcd_demosaic(rx, ry, rawdata, red, green, blue, cfarray, progress);
			break;
		case BAYER_AHD:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = ahd_demosaic(rx, ry, rawdata, red, green, blue, cfarray, rgb_cam, progress);
			break;
		//case BAYER_AMAZE:
			// retval = amaze_demosaic // need documentation about arguments
		case BAYER_DCB:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = dcb_demosaic(rx, ry, rawdata, red, green, blue, cfarray, progress, 1, TRUE);
			break;
		case BAYER_HPHD:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = hphd_demosaic(rx, ry, rawdata, red, green, blue, cfarray, progress);
			break;
		case BAYER_IGV:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = igv_demosaic(rx, ry, rawdata, red, green, blue, cfarray, progress);
			break;
		case BAYER_LMMSE:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = lmmse_demosaic(rx, ry, rawdata, red, green, blue, cfarray, progress, 1);
			// need documentation about last argument, 'iterations'
			break;
		case XTRANS:
			//retval = xtransfast_demosaic(rx, ry, rawdata, red, green, blue, xtrans, progress);
			retval = markesteijn_demosaic(rx, ry, rawdata, red, green, blue, xtrans, rgb_cam, progress, 1, TRUE, 16, FALSE);
			break;
	}
	free(rawdata[0]);
	free(rawdata);

	// 4. get the result in WORD (memory size: 3 times original)
	WORD *newfitdata = (WORD *)malloc(n * sizeof(WORD));
	if (!newfitdata) {
		PRINT_ALLOC_ERR;
		retval = RP_MEMORY_ERROR;
	}
	else {
		for (j = 0; j < n; j++) {
			newfitdata[j] = roundf_to_WORD(newdata[j]);
			// this rounding ^ is required because librtprocess
			// often returns data out of expected range
		}
	}

	free(newdata);
	free(blue);
	free(green);
	free(red);

	return retval == RP_NO_ERROR ? newfitdata : NULL;
}

// Warning: buf may be destroyed in case of failure, to avoid data duplication
// freeing buf is left to the caller
float *debayer_buffer_new_float(float *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, unsigned int xtrans[6][6]) {

	// super-pixel is handled by siril code, not librtprocess
	if (interpolation == BAYER_SUPER_PIXEL)
		return debayer_buffer_superpixel_float(buf, width, height, pattern);

	unsigned int cfarray[2][2];
	float rgb_cam[3][4] = { 1.0f };	// our white balance: we don't care
	int i, rx = *width, ry = *height;
	long j, nbpixels = rx * ry;
	long n = nbpixels * 3;
	// 1. prepare input data for librtprocess
	float **rawdata = (float **)malloc(ry * sizeof(float *));
	if (!rawdata) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	rawdata[0] = buf; // no duplication, input will be overwritten
	// TODO: do the following only for interpolations that need a conversion
	// TODO: vectorize!
#ifdef SIRIL_OUTPUT_DEBUG
	float min = 100000.0f, max = 0.0f;
#endif
	for (j = 0; j < nbpixels; j++) {
		if (buf[j] < 0.0f)
			buf[j] = 0.0f;
		else buf[j] *= 65535.0f;
#ifdef SIRIL_OUTPUT_DEBUG
		if (buf[j] > max) max = buf[j];
		if (buf[j] < min) min = buf[j];
#endif
	}
#ifdef SIRIL_OUTPUT_DEBUG
	fprintf(stdout, "****** before debayer, data is [%f, %f] (should be [0, 65535]) ******\n", min, max);
#endif

	for (i=1; i<ry; i++)
		rawdata[i] = rawdata[i - 1] + rx;

	// 2. allocate the demosaiced image buffer
	float *newdata = (float *)malloc(n * sizeof(float));
	if (!newdata) {
		free(rawdata);
		PRINT_ALLOC_ERR;
		return NULL;
	}

	float **red = (float **)malloc(ry * sizeof(float *));
	red[0] = newdata;
	for (i=1; i<ry; i++)
		red[i] = red[i - 1] + rx;

	float **green = (float **)malloc(ry * sizeof(float *));
	green[0] = red[0] + nbpixels;
	for (i=1; i<ry; i++)
		green[i] = green[i - 1] + rx;

	float **blue = (float **)malloc(ry * sizeof(float *));
	blue[0] = green[0] + nbpixels;
	for (i=1; i<ry; i++)
		blue[i] = blue[i - 1] + rx;

	// 3. process
	siril_debug_print("calling librtprocess float (%d)\n", interpolation);
	rpError retval;
	switch (interpolation) {
		case BAYER_VNG:
			pattern_to_cfarray2(com.debayer.bayer_pattern, cfarray);
			retval = vng4_demosaic(rx, ry, rawdata, red, green, blue, cfarray, progress);
			break;
		case BAYER_BILINEAR:
//		case BAYER_NEARESTNEIGHBOR:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			/* bayerfast: This demosaicer is not intended for final
			 * output, only for fast preview. */
			retval = bayerfast_demosaic(rx, ry, rawdata, red, green, blue, cfarray, progress, 1.0);
			break;
		default:
		case BAYER_RCD:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = rcd_demosaic(rx, ry, rawdata, red, green, blue, cfarray, progress);
			break;
		case BAYER_AHD:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = ahd_demosaic(rx, ry, rawdata, red, green, blue, cfarray, rgb_cam, progress);
			break;
		//case BAYER_AMAZE:
			// retval = amaze_demosaic // need documentation about arguments
		case BAYER_DCB:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = dcb_demosaic(rx, ry, rawdata, red, green, blue, cfarray, progress, 1, TRUE);
			break;
		case BAYER_HPHD:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = hphd_demosaic(rx, ry, rawdata, red, green, blue, cfarray, progress);
			break;
		case BAYER_IGV:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = igv_demosaic(rx, ry, rawdata, red, green, blue, cfarray, progress);
			break;
		case BAYER_LMMSE:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = lmmse_demosaic(rx, ry, rawdata, red, green, blue, cfarray, progress, 1);
			// need documentation about last argument, 'iterations'
			break;
		case XTRANS:
			//retval = xtransfast_demosaic(rx, ry, rawdata, red, green, blue, xtrans, progress);
			retval = markesteijn_demosaic(rx, ry, rawdata, red, green, blue, xtrans, rgb_cam, progress, 1, TRUE, 16, FALSE);
			break;
	}
	free(rawdata);

	// 4. convert back to siril range if needed
	// TODO: do the following only for interpolations that needed a conversion
#ifdef SIRIL_OUTPUT_DEBUG
	min = 100000.0f; max = 0.0f;
#endif
	for (j = 0; j < n; j++) {
#ifdef SIRIL_OUTPUT_DEBUG
		if (newdata[j] > max)
			max = newdata[j];
		if (newdata[j] < min)
			min = newdata[j];
#endif
		newdata[j] = newdata[j] * 1.52590219e-5f; // 1/65535
	}
#ifdef SIRIL_OUTPUT_DEBUG
	fprintf(stdout, "****** after debayer, data is [%f, %f] (should be [0, 65535]) ******\n", min, max);
#endif

	free(blue);
	free(green);
	free(red);

	return retval == RP_NO_ERROR ? newdata : NULL;
}

