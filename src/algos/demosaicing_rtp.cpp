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

/* This function retrieve the xtrans matrix from the FITS header */
static int retrieveXTRANSPattern(char *bayer, unsigned int xtrans[6][6]) {
	int x, y, i = 0;

	if (strlen(bayer) != 36) {
		siril_log_color_message(_("FITS header does not contain a proper XTRANS pattern, demosaicing cannot be done"), "red");
		return 1;
	}

	for (x = 0; x < 6; x++) {
		for (y = 0; y < 6; y++) {
			switch (bayer[i]) {
				default:	// shouldn't default be an error?
				case 'R':
					xtrans[x][y] = 0;
					break;
				case 'G':
					xtrans[x][y] = 1;
					break;
				case 'B':
					xtrans[x][y] = 2;
					break;
			}
			i++;
		}
	}
	return 0;
}

static bool progress(double p) {
	// p is [0, 1] progress of the debayer process
	return true;
}

int debayer(fits* fit, interpolation_method interpolation) {
	rpError retval;
	unsigned int i, cfarray[2][2], xtrans_array[6][6];
	long j, nbpixels = fit->naxes[0] * fit->naxes[1];
	long n = nbpixels * fit->naxes[2];
	float rgb_cam[3][4] = { 1.0f };	// our white balance: we don't care

	/* prepare the buffer */
	fprintf(stdout, "starting librtprocess debayer\n");

	float **rawdata = (float **)malloc(fit->ry * sizeof(float *));
	if (!rawdata) {
		PRINT_ALLOC_ERR;
		return -1;
	}
	if (fit->type == DATA_USHORT) {
		rawdata[0] = (float *)malloc(nbpixels * sizeof(float));
		if (!rawdata[0]) {
			PRINT_ALLOC_ERR;
			free(rawdata);
			return -1;
		}
		// TODO: vectorize!
		for (j = 0; j < n; j++)
			rawdata[0][j] = (float)fit->data[j];
	}
	else if (fit->type == DATA_FLOAT) {
		rawdata[0] = fit->fdata;
		/* some demosaicing functions support [0, 1] range with a factor
		 * given in argument, but not all so for now we'll do it here */
		// TODO: vectorize!
#ifdef SIRIL_OUTPUT_DEBUG
		float min = 100000.0f, max = 0.0f;
#endif
		for (j = 0; j < n; j++) {
			if (fit->fdata[j] < 0.0f)
				fit->fdata[j] = 0.0f;
			else fit->fdata[j] = fit->fdata[j] * 65535.0f;
#ifdef SIRIL_OUTPUT_DEBUG
			if (fit->fdata[j] > max)
				max = fit->fdata[j];
			if (fit->fdata[j] < min)
				min = fit->fdata[j];
#endif
		}
#ifdef SIRIL_OUTPUT_DEBUG
		fprintf(stdout, "****** before debayer, data is [%f, %f] (should be [0, 65535]) ******\n", min, max);
#endif
	} else {
		free(rawdata);
		return -1;
	}

	for (i=1; i<fit->ry; i++)
		rawdata[i] = rawdata[i - 1] + fit->rx;

	// allocate the demosaiced image buffer
	n *= 3;
	float *newdata = (float *)malloc(n * sizeof(float));
	if (!newdata) {
		free(rawdata);
		PRINT_ALLOC_ERR;
		return -1;
	}

	float **red = (float **)malloc(fit->ry * sizeof(float *));
	red[0] = newdata;
	for (i=1; i<fit->ry; i++)
		red[i] = red[i - 1] + fit->rx;

	float **green = (float **)malloc(fit->ry * sizeof(float *));
	green[0] = red[0] + nbpixels;
	for (i=1; i<fit->ry; i++)
		green[i] = green[i - 1] + fit->rx;

	float **blue = (float **)malloc(fit->ry * sizeof(float *));
	blue[0] = green[0] + nbpixels;
	for (i=1; i<fit->ry; i++)
		blue[i] = blue[i - 1] + fit->rx;

	/* process */
	switch (interpolation) {
		case BAYER_VNG:
			pattern_to_cfarray2(com.debayer.bayer_pattern, cfarray);
			retval = vng4_demosaic(fit->rx, fit->ry, rawdata, red, green, blue, cfarray, progress);
			break;
		case BAYER_BILINEAR:
		case BAYER_NEARESTNEIGHBOR:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			/* bayerfast: This demosaicer is not intended for final
			 * output, only for fast preview. */
			retval = bayerfast_demosaic(fit->rx, fit->ry, rawdata, red, green, blue, cfarray, progress, 1.0);
			break;
		default:
		case BAYER_RCD:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = rcd_demosaic(fit->rx, fit->ry, rawdata, red, green, blue, cfarray, progress);
			break;
		case BAYER_AHD:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = ahd_demosaic(fit->rx, fit->ry, rawdata, red, green, blue, cfarray, rgb_cam, progress);
			break;
		//case BAYER_AMAZE:
			// retval = amaze_demosaic // need documentation about arguments
		case BAYER_DCB:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = dcb_demosaic(fit->rx, fit->ry, rawdata, red, green, blue, cfarray, progress, 1, TRUE);
			break;
		case BAYER_HPHD:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = hphd_demosaic(fit->rx, fit->ry, rawdata, red, green, blue, cfarray, progress);
			break;
		case BAYER_IGV:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = igv_demosaic(fit->rx, fit->ry, rawdata, red, green, blue, cfarray, progress);
			break;
		case BAYER_LMMSE:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = lmmse_demosaic(fit->rx, fit->ry, rawdata, red, green, blue, cfarray, progress, 1);
			// need documentation about last argument, 'iterations'
			break;
		//case BAYER_SUPER_PIXEL:
			// is it supported by librtprocess? our old code below:
			// retval = super_pixel(buf, newbuf, *width, *height, pattern);
		case XTRANS:
			retrieveXTRANSPattern(fit->bayer_pattern, xtrans_array);
			//retval = xtransfast_demosaic(fit->rx, fit->ry, rawdata, red, green, blue, xtrans_array, progress);
			retval = markesteijn_demosaic(fit->rx, fit->ry, rawdata, red, green, blue, xtrans_array, rgb_cam, progress, 1, TRUE, 16, FALSE);
			break;
	}

	fprintf(stdout, "saving librtprocess debayer\n");
	/* get the result */
	if (fit->type == DATA_USHORT) {
		WORD *newfitdata = (WORD *)realloc(fit->data, n * sizeof(WORD));
		if (!newfitdata) {
			PRINT_ALLOC_ERR;
			retval = RP_MEMORY_ERROR;
			goto free_and_return;
		}
		for (j = 0; j < n; j++)
			newfitdata[j] = (WORD)newdata[j];	// or is rounding required?
		fit->naxes[2] = 3;
		fit->naxis = 3;
		free(newdata);
		fit_replace_buffer(fit, newfitdata, DATA_USHORT);
	}
	else if (fit->type == DATA_FLOAT) {
		// TODO: vectorize!
#ifdef SIRIL_OUTPUT_DEBUG
		float min = 100000.0f, max = 0.0f;
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
		fit->naxes[2] = 3;
		fit->naxis = 3;
		free(fit->fdata);
		fit_replace_buffer(fit, newdata, DATA_FLOAT);
	}

free_and_return:
	free(blue);
	free(green);
	free(red);
	free(rawdata);

	return retval != RP_NO_ERROR;
}

