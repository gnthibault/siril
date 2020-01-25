#include "librtprocess.h"
#include "core/siril.h"
#include "algos/demosaicing.h"
#ifdef __cplusplus
extern "C" {
#endif
#include "core/proto.h"
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

int debayer(fits* fit, interpolation_method interpolation, gboolean stretch_cfa) {
	rpError retval;
	unsigned int i, cfarray[2][2];
	long j, nbpixels = fit->naxes[0] * fit->naxes[1];
	long n = nbpixels * fit->naxes[2];

	/* prepare the buffer */
	if (fit->type == DATA_USHORT) {
		return 1;
	}
	fprintf(stdout, "starting librtprocess debayer\n");

	float **rawdata = (float **)malloc(fit->ry * sizeof(float *));
	if (!rawdata) {
		PRINT_ALLOC_ERR;
		return -1;
	}
	rawdata[0] = fit->fdata;
	for (i=1; i<fit->ry; i++)
		rawdata[i] = rawdata[i - 1] + fit->rx;
	// TODO: vectorize!
	for (j = 0; j < n; j++)
		fit->fdata[j] = fit->fdata[j] * 65535.0f;

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
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);

			retval = bayerfast_demosaic(fit->rx, fit->ry, rawdata, red, green, blue, cfarray, progress, 1.0);
			break;
		default:
		case BAYER_RCD:
			pattern_to_cfarray(com.debayer.bayer_pattern, cfarray);
			retval = rcd_demosaic(fit->rx, fit->ry, rawdata, red, green, blue, cfarray, progress);
			break;
			/* 
			   ahd_demosaic
			   amaze_demosaic
			   bayerfast_demosaic
			   dcb_demosaic
			   hphd_demosaic
			   igv_demosaic
			   lmmse_demosaic
			   vng4_demosaic
			   markesteijn_demosaic
			   xtransfast_demosaic
			   */
	}

	/* get the result */
	// TODO: vectorize!
	for (j = 0; j < n; j++)
		newdata[j] = newdata[j] * 1.52590219e-5f; // 1/65535

	fprintf(stdout, "saving ibrtprocess debayer\n");
	free(blue);
	free(green);
	free(red);
	free(rawdata);
	fit->naxes[2] = 3;
	fit->naxis = 3;
	free(fit->fdata);
	fit_replace_buffer(fit, newdata, DATA_FLOAT);

	return retval != RP_NO_ERROR;
}

