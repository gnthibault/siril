/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2019 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

/* Management of Siril's internal image format: unsigned 16-bit FITS */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <gsl/gsl_statistics.h>
#ifdef _WIN32
#include <windows.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "io/sequence.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "algos/statistics.h"
#include "io/single_image.h"

static char *MIPSHI[] = {"MIPS-HI", "CWHITE", "DATAMAX", NULL };
static char *MIPSLO[] = {"MIPS-LO", "CBLACK", "DATAMIN", NULL };
static char *PixSizeX[] = { "XPIXSZ", "XPIXELSZ", "PIXSIZE1", NULL };
static char *PixSizeY[] = { "YPIXSZ", "YPIXELSZ", "PIXSIZE2", NULL };
static char *BinX[] = { "XBINNING", "BINX", NULL };
static char *BinY[] = { "YBINNING", "BINY", NULL };
static char *Focal[] = { "FOCAL", "FOCALLEN", NULL };
static char *Exposure[] = { "EXPTIME", "EXPOSURE", NULL };

#define __tryToFindKeywords(fptr, type, keywords, value) \
{ \
	int __iter__ = 0; \
	int __status__; \
	do { \
		__status__ = 0; \
		fits_read_key(fptr, type, keywords[__iter__], value, NULL, &__status__); \
		__iter__++; \
	} while ((keywords[__iter__]) && (__status__ > 0)); \
}

static void read_fits_date_obs_header(fits *fit) {
	int status = 0;
	fits_read_key(fit->fptr, TSTRING, "DATE-OBS", &(fit->date_obs), NULL, &status);

	status = 0;
	char ut_start[FLEN_VALUE];
	/** Case seen in some FITS files. Needed to get date back in SER conversion **/
	fits_read_key(fit->fptr, TSTRING, "UT-START", &ut_start, NULL,
				&status);
	if (ut_start[0] != '\0' && fit->date_obs[2] == G_DIR_SEPARATOR) {
		int year, month, day;
		sscanf(fit->date_obs, "%02d/%02d/%04d", &day, &month, &year);
		g_snprintf(fit->date_obs, sizeof(fit->date_obs), "%04d-%02d-%02dT%s",
				year, month, day, ut_start);
	}
}

static int fit_stats(fits *fit, double *mini, double *maxi) {
	int status = 0;
	int ii, anaxis;
	long npixels = 1;
	long anaxes[3] = {1,1,1}, firstpix[3] = {1,1,1};
	double *pix, sum = 0.;
	double meanval = 0., minval = 1.E33, maxval = -1.E33;

	/* initialize value in case where it does not work */
	*mini = 0;
	*maxi = 0;

	fits_get_img_dim(fit->fptr, &anaxis, &status);
	fits_get_img_size(fit->fptr, 3, anaxes, &status);

    if (status) {
       fits_report_error(stderr, status); /* print error message */
       return(status);
    }

    npixels = anaxes[0];  /* no. of pixels to read in each row */
	pix = (double *) malloc(npixels * sizeof(double)); /* memory for 1 row */
	if (pix == NULL) {
		PRINT_ALLOC_ERR;
		return (1);
	}

	/* loop over all planes of the cube (2D images have 1 plane) */
	for (firstpix[2] = 1; firstpix[2] <= anaxes[2]; firstpix[2]++) {
		/* loop over all rows of the plane */
		for (firstpix[1] = 1; firstpix[1] <= anaxes[1]; firstpix[1]++) {
			/* give starting pixel coordinate and number of pixels to read */
			if (fits_read_pix(fit->fptr, TDOUBLE, firstpix, npixels, NULL, pix,
					NULL, &status))
				break; /* jump out of loop on error */

			for (ii = 0; ii < npixels; ii++) {
				sum += pix[ii]; /* accumulate sum */
				if (pix[ii] < minval)
					minval = pix[ii]; /* find min and  */
				if (pix[ii] > maxval)
					maxval = pix[ii]; /* max values    */
			}
		}
	}    /* end of loop over planes */
	free(pix);

	if (status) {
		fits_report_error(stderr, status); /* print any error message */
	} else {
		if (npixels > 0)
			meanval = sum / npixels;
		siril_debug_print("  sum of pixels = %g\n", sum);
		siril_debug_print("  mean value    = %g\n", meanval);
		siril_debug_print("  minimum value = %g\n", minval);
		siril_debug_print("  maximum value = %g\n", maxval);
		*maxi = maxval;
		*mini = minval;
	}
	return status;
}

/* copy the complete header in a heap-allocated string */
static void fits_read_history(fitsfile *fptr, GSList **history, int *status) {
	int i, hdupos, nkeys;
	char card[FLEN_CARD];
	GSList *list = NULL;

	fits_get_hdu_num(fptr, &hdupos); /*Get the current HDU position */
	for (; !*status; hdupos++) {
		fits_get_hdrspace(fptr, &nkeys, NULL, status);
		for (i = 1; i <= nkeys; i++) {
			if (fits_read_record(fptr, i, card, status))
				break;
			if (!strncmp(card, "HISTORY", 7)) {
				list = g_slist_prepend(list, g_strdup(card + 8));
			}
		}
		fits_movrel_hdu(fptr, 1, NULL, status);
	}
	if (*history)
		g_slist_free_full(*history, free);
	list = g_slist_reverse(list);
	*history = list;
}

/* reading the FITS header to get useful information
 * stored in the fit, requires an opened file descriptor */
static void read_fits_header(fits *fit) {
	/* about the status argument: http://heasarc.gsfc.nasa.gov/fitsio/c/c_user/node28.html */
	int status = 0;
	double zero, mini, maxi;

	fit_stats(fit, &mini, &maxi);

	__tryToFindKeywords(fit->fptr, TUSHORT, MIPSHI, &fit->hi);
	__tryToFindKeywords(fit->fptr, TUSHORT, MIPSLO, &fit->lo);

	if (fit->orig_bitpix == SHORT_IMG) {
		if (fit->lo)
			fit->lo += 32768;
		if (fit->hi)
			fit->hi += 32768;
	}

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "BSCALE", &zero, NULL, &status);
	if (!status && 1.0 != zero) {
		siril_log_message(_("Loaded FITS file has "
				"a BSCALE different than 1 (%f)\n"), zero);
		status = 0;
		/* We reset the scaling factors as we don't use it */
		fits_set_bscale(fit->fptr, 1.0, 0.0, &status);
	}
	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "BZERO", &zero, NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "DATAMAX", &(fit->data_max), NULL, &status);
	if (status == KEY_NO_EXIST) {
		fit->data_max = maxi;
	}

	status = 0;
	fits_read_key(fit->fptr, TUSHORT, "MAXCFA", &(fit->maximum_pixel_value),
			NULL, &status);


	/*******************************************************************
	 * ************* CAMERA AND INSTRUMENT KEYWORDS ********************
	 * ****************************************************************/

	__tryToFindKeywords(fit->fptr, TFLOAT, PixSizeX, &fit->pixel_size_x);
	__tryToFindKeywords(fit->fptr, TFLOAT, PixSizeY, &fit->pixel_size_y);
	__tryToFindKeywords(fit->fptr, TUINT, BinX, &fit->binning_x);
	if (fit->binning_x <= 0) fit->binning_x = 1;
	__tryToFindKeywords(fit->fptr, TUINT, BinY, &fit->binning_y);
	if (fit->binning_y <= 0) fit->binning_y = 1;

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "INSTRUME", &(fit->instrume), NULL,
			&status);

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "TELESCOP", &(fit->telescop), NULL,
			&status);

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "OBSERVER", &(fit->observer), NULL,
			&status);

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "BAYERPAT", &(fit->bayer_pattern), NULL,
			&status);

	status = 0;
	fits_read_key(fit->fptr, TINT, "XBAYROFF", &(fit->bayer_xoffset), NULL,
			&status);

	status = 0;
	fits_read_key(fit->fptr, TINT, "YBAYROFF", &(fit->bayer_yoffset), NULL,
			&status);

	read_fits_date_obs_header(fit);

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "DATE", &(fit->date), NULL,
			&status);

	__tryToFindKeywords(fit->fptr, TDOUBLE, Focal, &fit->focal_length);
	if (!sequence_is_loaded() || com.seq.current == 0)
		fprintf(stdout,
				"Read from FITS header: pix size %gx%g, binning %ux%u, focal %g\n",
				fit->pixel_size_x, fit->pixel_size_y, fit->binning_x,
				fit->binning_y, fit->focal_length);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CCD-TEMP", &(fit->ccd_temp), NULL,
			&status);	// Non-standard keywords used in MaxIm DL

	__tryToFindKeywords(fit->fptr, TDOUBLE, Exposure, &fit->exposure);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "APERTURE", &(fit->aperture), NULL,
			&status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "ISOSPEED", &(fit->iso_speed), NULL,
			&status);	// Non-standard keywords used in MaxIm DL

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CVF", &(fit->cvf), NULL, &status);

	/*******************************************************************
	 * ******************* PLATE SOLVING KEYWORDS **********************
	 * ****************************************************************/
	status = 0;
	fits_read_key(fit->fptr, TUINT, "EQUINOX", &(fit->wcs.equinox), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "OBJCTRA", &(fit->wcs.objctra), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "OBJCTDEC", &(fit->wcs.objctdec), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CRPIX1", &(fit->wcs.crpix1), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CRPIX2", &(fit->wcs.crpix2), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CRVAL1", &(fit->wcs.crval1), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CRVAL2", &(fit->wcs.crval2), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CDELT1", &(fit->wcs.cdelt1), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CDELT2", &(fit->wcs.cdelt2), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CROTA1", &(fit->wcs.crota1), NULL, &status);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "CROTA2", &(fit->wcs.crota2), NULL, &status);


	/*******************************************************************
	 * ************************* DFT KEYWORDS **************************
	 * ****************************************************************/

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "DFTNORM0", &(fit->dft.norm[0]), NULL,
			&status);
	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "DFTNORM1", &(fit->dft.norm[1]), NULL,
			&status);
	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "DFTNORM2", &(fit->dft.norm[2]), NULL,
			&status);

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "DFTORD", &(fit->dft.ord), NULL,
			&status);

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "DFTTYPE", &(fit->dft.type), NULL,
			&status);

	status = 0;
	fits_read_history(fit->fptr, &(fit->history), &status);
}

/* copy the complete header in a heap-allocated string */
static char *copy_header(fits *fit) {
	int status = 0;
	int i, hdupos, nkeys, strsize, strlength;
	char card[FLEN_CARD];
	char *header;

	/* each line in the FITS header is 80 character wide
	 * in our string, we also keep new lines, so that's 81.
	 * initial allocation is 20 times 81 = 1620
	 * reallocations are 7 lines, 567 */
	strsize = 1620;
	header = malloc(strsize);
	if (header == NULL) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	header[0] = '\0';
	strlength = 0;
	fits_get_hdu_num(fit->fptr, &hdupos); /*Get the current HDU position */
	for (; !status; hdupos++) {
		fits_get_hdrspace(fit->fptr, &nkeys, NULL, &status);
		//fprintf(stdout, "Header listing for HDU #%d:\n", hdupos);
		for (i = 1; i <= nkeys; i++) {
			int cardlen;
			char *newstr;
			if (fits_read_record(fit->fptr, i, card, &status))
				break;
			cardlen = strlen(card);
			if (strlength + cardlen + 1 >= strsize) {
				strsize += 567;
				newstr = realloc(header, strsize);
				if (newstr == NULL) {
					PRINT_ALLOC_ERR;
					free(header);
					return NULL;
				}
				header = newstr;
			}
			strcpy(header + strlength, card);
			strlength += cardlen;
			strcpy(header + strlength, "\n");
			strlength++;
		}
		fits_movrel_hdu(fit->fptr, 1, NULL, &status);
	}

	if (header[0] == '\0') {
		free(header);
		header = NULL;
	}
	if (!header)
		return NULL;
	if (!sequence_is_loaded() || com.seq.current == 0)
		siril_debug_print("%s", header);// don't display for all frames of a sequence
	return header;
}

static void report_fits_error(int status) {
	if (status) {
		char errmsg[FLEN_ERRMSG];
		while (fits_read_errmsg(errmsg)) {
			siril_log_message(_("FITS error: %s\n"), errmsg);
		}
	}
}

static void conv_8_to_16(WORD *data, unsigned int nbdata) {
	int i;

	for (i = 0; i < nbdata; i++) {
		double tmp = (double) data[i] / UCHAR_MAX_DOUBLE * USHRT_MAX_DOUBLE;
		data[i] = round_to_WORD(tmp);
	}
}

/* convert FITS data formats to siril native.
 * nbdata is the number of pixels, w * h.
 * from is not freed, to must be allocated and can be the same as from */
static void convert_data(int bitpix, const void *from, WORD *to, unsigned int nbdata, gboolean values_above_1) {
	int i;
	BYTE *data8;
	int16_t *data16;
	double *pixels_double;
	double norm = 1.0;
	long *sdata32;	// TO BE TESTED on 32-bit arch, seems to be a cfitsio bug
	unsigned long *data32;

	switch (bitpix) {
		case BYTE_IMG:
			data8 = (BYTE *)from;
			for (i = 0; i < nbdata; i++)
				to[i] = (WORD)data8[i];
			break;
		case USHORT_IMG:	// siril 0.9 native
			// nothing to do
			break;
		case SHORT_IMG:
			// add 2^15 to the read data to obtain unsigned
			data16 = (int16_t *)from;
			for (i = 0; i < nbdata; i++) {
				int sum = 32768 + (int)data16[i];
				to[i] = (WORD)sum;
			}
			break;
		case ULONG_IMG:		// 32-bit unsigned integer pixels
			data32 = (unsigned long *)from;
			for (i = 0; i < nbdata; i++)
				to[i] = (WORD)(data32[i] >> 16);
			break;
		case LONG_IMG:		// 32-bit signed integer pixels
			sdata32 = (long *)from;
			for (i = 0; i < nbdata; i++)
				to[i] = (WORD)((sdata32[i] >> 16) + 32768);
			break;
		case DOUBLE_IMG:	// 64-bit floating point pixels
		case FLOAT_IMG:		// 32-bit floating point pixels
			pixels_double = (double *)from;
			/* various data values can be found in a float or
			 * double image. Sometimes it's normalized between 0
			 * and 1, but sometimes, DATAMIN and DATAMAX give the range.
			 */
			if (!values_above_1) norm = USHRT_MAX_DOUBLE;
			for (i = 0; i < nbdata; i++) {
				to[i] = round_to_WORD(norm * pixels_double[i]);
			}
			break;
		case LONGLONG_IMG:	// 64-bit integer pixels
		default:
			siril_log_message(_("Unknown FITS data format in internal conversion\n"));
	}
}

/* read buffer from an already open FITS file, fit should have all metadata
 * correct, and convert the buffer to fit->data with the given type, which
 * currently should be TBYTE or TUSHORT because fit doesn't contain other data.
 * filename is for error reporting
 */
static int read_fits_with_convert(fits* fit, const char* filename) {
	int status = 0, zero = 0, datatype;
	BYTE *data8;
	double *pixels_double;
	uint32_t *pixels_long;
	long orig[3] = { 1L, 1L, 1L };
	// orig ^ gives the coordinate in each dimension of the first pixel to be read
	unsigned int nbdata = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];

	fits_movabs_hdu(fit->fptr, 1, 0, &status); // make sure reading primary HDU

	switch (fit->bitpix) {
	case BYTE_IMG:
		data8 = malloc(nbdata * sizeof(BYTE));
		datatype = fit->bitpix == BYTE_IMG ? TBYTE : TSBYTE;
		fits_read_pix(fit->fptr, datatype, orig, nbdata, &zero,
				data8, &zero, &status);
		if (status) break;
		convert_data(fit->bitpix, data8, fit->data, nbdata, FALSE);
		free(data8);
		break;
	case SHORT_IMG:
		fits_read_pix(fit->fptr, TSHORT, orig, nbdata, &zero,
				fit->data, &zero, &status);
		if (status) break;
		convert_data(fit->bitpix, fit->data, fit->data, nbdata, FALSE);
		fit->bitpix = USHORT_IMG;
		break;
	case USHORT_IMG:
		// siril 0.9 native, no conversion required
		fits_read_pix(fit->fptr, TUSHORT, orig, nbdata, &zero,
				fit->data, &zero, &status);
		if (status == NUM_OVERFLOW) {
			// in case there are errors, we try short data
			status = 0;
			fits_read_pix(fit->fptr, TSHORT, orig, nbdata, &zero, fit->data,
					&zero, &status);
			if (status)
				break;
			convert_data(SHORT_IMG, fit->data, fit->data, nbdata, FALSE);
			if (fit->lo)
				fit->lo += 32768;
			if (fit->hi)
				fit->hi += 32768;
			fit->bitpix = USHORT_IMG;
		}
		break;
	case ULONG_IMG:		// 32-bit unsigned integer pixels
	case LONG_IMG:		// 32-bit signed integer pixels
		pixels_long = malloc(nbdata * sizeof(uint32_t));
		status = 0;
		datatype = fit->bitpix == LONG_IMG ? TLONG : TULONG;
		fits_read_pix(fit->fptr, datatype, orig, nbdata, &zero,
				pixels_long, &zero, &status);
		if (status) break;
		convert_data(fit->bitpix, pixels_long, fit->data, nbdata, FALSE);
		free(pixels_long);
		fit->bitpix = USHORT_IMG;
		break;
	case DOUBLE_IMG:	// 64-bit floating point pixels
	case FLOAT_IMG:		// 32-bit floating point pixels
		pixels_double = malloc(nbdata * sizeof(double));
		fits_read_pix(fit->fptr, TDOUBLE, orig, nbdata, &zero,
				pixels_double, &zero, &status);
		if (status) break;
		convert_data(fit->bitpix, pixels_double, fit->data, nbdata, fit->data_max > 1.0);
		free(pixels_double);
		fit->bitpix = USHORT_IMG;
		break;
	case LONGLONG_IMG:	// 64-bit integer pixels
	default:
		siril_log_message(_("FITS image format %d is not supported by Siril.\n"), fit->bitpix);
		return -1;
	}

	if (status) {
		siril_log_message(_("Fitsio error reading data, file: %s.\n"), filename);
		report_fits_error(status);
		return -1;
	}

	return 0;
}

/* This function reads partial data on one layer from the opened FITS and
 * convert it to siril's format (USHORT) */
static int internal_read_partial_fits(fitsfile *fptr, unsigned int ry,
		int bitpix, WORD *dest, gboolean values_above_1, int layer,
		const rectangle *area) {
	int datatype;
	BYTE *data8;
	double *pixels_double;
	long *pixels_long;
	long fpixel[3], lpixel[3], inc[3] = { 1L, 1L, 1L };
	int zero = 0, status = 0;

	/* fpixel is first pixel, lpixel is last pixel, starts with value 1 */
	fpixel[0] = area->x + 1;        // in siril, it starts with 0
	fpixel[1] = ry - area->y - area->h + 1;
	fpixel[2] = layer + 1;
	lpixel[0] = area->x + area->w;  // with w and h at least 1, we're ok
	lpixel[1] = ry - area->y;
	lpixel[2] = layer + 1;

	unsigned int nbdata = area->w * area->h;

	switch (bitpix) {
		case BYTE_IMG:
			data8 = malloc(nbdata * sizeof(BYTE));
			datatype = bitpix == BYTE_IMG ? TBYTE : TSBYTE;
			fits_read_subset(fptr, datatype, fpixel, lpixel, inc, &zero, data8,
					&zero, &status);
			if (status) break;
			convert_data(bitpix, data8, dest, nbdata, FALSE);
			free(data8);
			break;
		case SHORT_IMG:
			fits_read_subset(fptr, TSHORT, fpixel, lpixel, inc, &zero, dest,
					&zero, &status);
			convert_data(bitpix, dest, dest, nbdata, FALSE);
			break;
		case USHORT_IMG:
			fits_read_subset(fptr, TUSHORT, fpixel, lpixel, inc, &zero, dest,
					&zero, &status);
			break;
		case ULONG_IMG:		// 32-bit unsigned integer pixels
		case LONG_IMG:		// 32-bit signed integer pixels
			pixels_long = malloc(nbdata * sizeof(long));
			status = 0;
			datatype = bitpix == LONG_IMG ? TLONG : TULONG;
			fits_read_subset(fptr, datatype, fpixel, lpixel, inc, &zero,
					pixels_long, &zero, &status);
			if (status) break;
			convert_data(bitpix, pixels_long, dest, nbdata, FALSE);
			free(pixels_long);
			break;
		case DOUBLE_IMG:	// 64-bit floating point pixels
		case FLOAT_IMG:		// 32-bit floating point pixels
			pixels_double = malloc(nbdata * sizeof(double));
			fits_read_subset(fptr, TDOUBLE, fpixel, lpixel, inc, &zero, pixels_double,
					&zero, &status);
			if (status) break;
			convert_data(bitpix, pixels_double, dest, nbdata, values_above_1);
			free(pixels_double);
			break;
		case LONGLONG_IMG:	// 64-bit integer pixels
		default:
			siril_log_message(_("FITS image format %d is not supported by Siril.\n"), bitpix);
			return -1;
	}
	return status;
}

static int siril_fits_create_diskfile(fitsfile **fptr, const char *filename, int *status) {
	gchar *localefilename = get_locale_filename(filename);
	fits_create_diskfile(fptr, localefilename, status);
	g_free(localefilename);
	return *status;
}

static void save_wcs_keywords(fits *fit) {
	int status = 0;

	if (fit->wcs.equinox > 0) {
		fits_update_key(fit->fptr, TUINT, "EQUINOX", &(fit->wcs.equinox),
						"Equatorial equinox", &status);
		status = 0;
		fits_update_key(fit->fptr, TSTRING, "CTYPE1", "RA---TA", "Coordinate type for the first axis", &status);
		status = 0;
		fits_update_key(fit->fptr, TSTRING, "CTYPE2", "DEC--TA", "Coordinate type for the second axis", &status);
		status = 0;
		fits_update_key(fit->fptr, TSTRING, "OBJCTRA", &(fit->wcs.objctra),	"Image center R.A. (hms)", &status);
		status = 0;
		fits_update_key(fit->fptr, TSTRING, "OBJCTDEC", &(fit->wcs.objctdec), "Image center declination (dms)", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CRPIX1", &(fit->wcs.crpix1), "Axis1 reference pixel", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CRPIX2", &(fit->wcs.crpix2), "Axis2 reference pixel", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CRVAL1", &(fit->wcs.crval1), "Axis1 reference value", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CRVAL2", &(fit->wcs.crval2), "Axis2 reference value", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CDELT1", &(fit->wcs.cdelt1), "Axis1 scale", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CDELT2", &(fit->wcs.cdelt2), "Axis2 scale", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CROTA1", &(fit->wcs.crota1), "Axis1 rotation angle (deg)", &status);
		status = 0;
		fits_update_key(fit->fptr, TDOUBLE, "CROTA2", &(fit->wcs.crota2), "Axis2 rotation angle (deg)", &status);
	}
}

static void save_fits_header(fits *fit) {
	int i, status = 0;
	double zero;
	char comment[FLEN_COMMENT];

	if (fit->hi) { /* may not be initialized */
		fits_update_key(fit->fptr, TUSHORT, "MIPS-HI", &(fit->hi),
				"Upper visualization cutoff ", &status);
		fits_update_key(fit->fptr, TUSHORT, "MIPS-LO", &(fit->lo),
				"Lower visualization cutoff ", &status);
	}
	status = 0;
	switch (fit->bitpix) {
	case BYTE_IMG:
	case SHORT_IMG:
		zero = 0.0;
		break;
	default:
	case USHORT_IMG:
		zero = 32768.0;
		break;
	}
	fits_update_key(fit->fptr, TDOUBLE, "BZERO", &zero,
			"offset data range to that of unsigned short", &status);

	status = 0;
	zero = 1.0;
	fits_update_key(fit->fptr, TDOUBLE, "BSCALE", &zero, "default scaling factor",
			&status);

	status = 0;
	if (fit->maximum_pixel_value)
		fits_update_key(fit->fptr, TUSHORT, "MAXCFA",
				&(fit->maximum_pixel_value), "maximum raw pixel value", &status);

	/*******************************************************************
	 * ************* CAMERA AND INSTRUMENT KEYWORDS ********************
	 * ******************** AND DATES **********************************
	 * ****************************************************************/

	status = 0;
	if (fit->instrume[0] != '\0')
		fits_update_key(fit->fptr, TSTRING, "INSTRUME", &(fit->instrume),
				"instrument name", &status);
	status = 0;
	if (fit->telescop[0] != '\0')
		fits_update_key(fit->fptr, TSTRING, "TELESCOP", &(fit->telescop),
				"telescope used to acquire this image", &status);
	status = 0;
	if (fit->observer[0] != '\0')
		fits_update_key(fit->fptr, TSTRING, "OBSERVER", &(fit->observer),
				"observer name", &status);
	status = 0;
	int itmp;
	char fit_date[40];
	fits_get_system_time(fit_date, &itmp, &status);
	fits_update_key(fit->fptr, TSTRING, "DATE", fit_date,
			"UTC date that FITS file was created", &status);

	status = 0;
	if (fit->date_obs[0] != '\0')
		fits_update_key(fit->fptr, TSTRING, "DATE-OBS", &(fit->date_obs),
				"YYYY-MM-DDThh:mm:ss observation start, UT", &status);

	/* all keywords below are non-standard */
	status = 0;
	if (fit->pixel_size_x > 0.)
		fits_update_key(fit->fptr, TFLOAT, "XPIXSZ", &(fit->pixel_size_x),
				"X pixel size microns", &status);
	if (fit->pixel_size_y > 0.)
		fits_update_key(fit->fptr, TFLOAT, "YPIXSZ", &(fit->pixel_size_y),
				"Y pixel size microns", &status);

	status = 0;
	if (fit->binning_x)
		fits_update_key(fit->fptr, TUINT, "XBINNING", &(fit->binning_x),
				"Camera binning mode", &status);
	if (fit->binning_y)
		fits_update_key(fit->fptr, TUINT, "YBINNING", &(fit->binning_y),
				"Camera binning mode", &status);

	status = 0;
	if (fit->focal_length > 0.)
		fits_update_key(fit->fptr, TDOUBLE, "FOCALLEN", &(fit->focal_length),
				"Camera focal length", &status);

	status = 0;
	if (fit->ccd_temp)
		fits_update_key(fit->fptr, TDOUBLE, "CCD-TEMP", &(fit->ccd_temp),
				"CCD temp in C", &status);

	status = 0;
	if (fit->exposure > 0.)
		fits_update_key(fit->fptr, TDOUBLE, "EXPTIME", &(fit->exposure),
				"Exposure time [s]", &status);

	status = 0;
	if (fit->aperture > 0.)
		fits_update_key(fit->fptr, TDOUBLE, "APERTURE", &(fit->aperture),
				"Aperture of the instrument", &status);

	status = 0;
	if (fit->iso_speed > 0.)
		fits_update_key(fit->fptr, TDOUBLE, "ISOSPEED", &(fit->iso_speed),
				"ISO camera setting", &status);

	status = 0;
	if (fit->bayer_pattern[0] != '\0') {
		fits_update_key(fit->fptr, TSTRING, "BAYERPAT", &(fit->bayer_pattern),
				"Bayer color pattern", &status);

		status = 0;
		fits_update_key(fit->fptr, TINT, "XBAYROFF", &(fit->bayer_xoffset),
				"X offset of Bayer array", &status);

		status = 0;
		fits_update_key(fit->fptr, TINT, "YBAYROFF", &(fit->bayer_yoffset),
				"Y offset of Bayer array", &status);
	}

	status = 0;
	if (fit->cvf > 0.)
		fits_update_key(fit->fptr, TDOUBLE, "CVF", &(fit->cvf),
				"Conversion factor (e-/adu)", &status);

	/*******************************************************************
	 * ******************* PLATE SOLVING KEYWORDS **********************
	 * ****************************************************************/

	save_wcs_keywords(fit);

	/*******************************************************************
	 * ******************** PROGRAMM KEYWORDS **************************
	 * ****************************************************************/

	status = 0;
	char programm[32];
	sprintf(programm, "%s v%s", PACKAGE, VERSION);
	programm[0] = toupper(programm[0]);			// convert siril to Siril
	fits_update_key(fit->fptr, TSTRING, "PROGRAM", programm,
			"Software that created this HDU", &status);

	/*******************************************************************
	 * ********************* HISTORY KEYWORDS **************************
	 * ****************************************************************/

	status = 0;
	if (fit->history) {
		GSList *list;
		for (list = fit->history; list; list = list->next) {
			fits_write_history(fit->fptr, (char *)list->data, &status);
		}
	}

	status = 0;
	if (com.history) {
		for (i = 0; i < com.hist_display; i++) {
			if (com.history[i].history[0] != '\0')
				fits_write_history(fit->fptr, com.history[i].history, &status);
		}
	}

	/*******************************************************************
	 * ************************* DFT KEYWORDS **************************
	 * ****************************************************************/

	status = 0;
	if (fit->dft.type[0] != '\0') {
		if (fit->dft.type[0] == 'S')
			strcpy(comment, "Module of a Discrete Fourier Transform");
		else if (fit->dft.type[0] == 'P')
			strcpy(comment, "Phase of a Discrete Fourier Transform");
		else
			status = 1;			// should not happen
		fits_update_key(fit->fptr, TSTRING, "DFTTYPE", &(fit->dft.type),
				comment, &status);
	}

	status = 0;
	if (fit->dft.ord[0] != '\0') {
		if (fit->dft.ord[0] == 'C')
			strcpy(comment, "Low spatial freq. are located at image center");
		else if (fit->dft.ord[0] == 'R')
			strcpy(comment, "High spatial freq. are located at image center");
		else
			status = 1;			// should not happen
		fits_update_key(fit->fptr, TSTRING, "DFTORD", &(fit->dft.ord), comment,
				&status);
	}

	for (i = 0; i < fit->naxes[2]; i++) {
		if (fit->dft.norm[i] > 0.) {
			char str1[] = "DFTNORM";
			char str2[] = "Normalisation value for channel #";
			char key_str[FLEN_KEYWORD], comment_str[FLEN_VALUE];
			sprintf(key_str, "%s%d", str1, i);
			sprintf(comment_str, "%s%d", str2, i);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, key_str, &(fit->dft.norm[i]),
					comment_str, &status);
		}
	}
}

/********************** public functions ************************************/

double get_exposure_from_fitsfile(fitsfile *fptr) {
	double exp = 0.0;

	__tryToFindKeywords(fptr, TDOUBLE, Exposure, &exp);
	return exp;
}

int import_metadata_from_fitsfile(fitsfile *fptr, fits *to) {
	fits from = { 0 };

	from.fptr = fptr;

	read_fits_header(&from);
	copy_fits_metadata(&from, to);
	return 0;
}

// return 0 on success, fills realname if not NULL with the opened file's name
int readfits(const char *filename, fits *fit, char *realname) {
	int status, retval;
	char *name = NULL;
	gchar *basename;
	image_type imagetype;
	unsigned int nbdata;
	double offset;

	fit->naxes[2] = 1; //initialization of the axis number before opening : NEED TO BE OPTIMIZED

	if (stat_file(filename, &imagetype, &name)) {
		siril_log_message(_("%s.[any_allowed_extension] not found.\n"),
				filename);
		free(name);
		return 1;
	}
	if (imagetype != TYPEFITS) {
		siril_log_message(
				_("The file %s is not a FITS file or doesn't exists with FITS extensions.\n"),
						filename);
		free(name);
		return 1;
	}

	if (realname)
		strcpy(realname, name);

	status = 0;
	siril_fits_open_diskfile(&(fit->fptr), name, READONLY, &status);
	if (status) {
		report_fits_error(status);
		free(name);
		return status;
	}
	free(name);

	status = 0;
	fits_get_img_param(fit->fptr, 3, &(fit->bitpix), &(fit->naxis), fit->naxes, &status);
	if (status) {
		siril_log_message(
				_("FITSIO error getting image parameters, file %s.\n"), filename);
		report_fits_error(status);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return status;
	}

	/* since USHORT_IMG is a cfitsio trick and doesn't really exist in the
	 * file, when we read it, the bitpix is given as SHORT_IMG. If BZERO is
	 * 2^15, it means that it is USHORT data and we force the bitpix to it
	 * in order to read it properly. Same thing for LONG and ULONG.
	 * https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node23.html
	 */
	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "BZERO", &offset, NULL, &status);
	if (!status) {
		if (fit->bitpix == SHORT_IMG && offset != 0.0) {
			fit->bitpix = USHORT_IMG;
		}
		else if (fit->bitpix == LONG_IMG && offset != 0.0)
			fit->bitpix = ULONG_IMG;
	} else {
		/* but some software just put unsigned 16-bit data in the file
		 * and don't set the BZERO keyword... */
		if (status == KEY_NO_EXIST && fit->bitpix == SHORT_IMG)
			fit->bitpix = USHORT_IMG;
	}
	// and we store the original bitpix to reuse it during later partial
	// reads when we have no access to the header or a fits * struct.
	fit->orig_bitpix = fit->bitpix;

	fit->rx = fit->naxes[0];
	fit->ry = fit->naxes[1];
	nbdata = fit->rx * fit->ry;

	if (fit->naxis == 3 && fit->naxes[2] != 3) {
		siril_log_message(_("Unsupported FITS image with %ld channels.\n"),
				fit->naxes[2]);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return -1;
	}
	// comment the above if and uncomment below for tests with 4 channels or more
	//if (fit->naxis == 3) fit->naxes[2] = 3;

	if (fit->naxis == 2 && fit->naxes[2] == 0) {
		fit->naxes[2] = 1;
		/* naxes[2] is set to 1 because:
		 * - it doesn't matter, since naxis is 2, it's not used
		 * - it's very convenient to use it in multiplications as the number of layers
		 */
	}
	if (fit->bitpix == LONGLONG_IMG) {
		siril_log_message(
				_("FITS images with 64 bits signed integer per pixel.channel are not supported.\n"));
		status = 0;
		fits_close_file(fit->fptr, &status);
		return -1;
	}

	/* realloc fit->data to the image size */
	WORD *olddata = fit->data;
	if ((fit->data = realloc(fit->data, nbdata * fit->naxes[2] * sizeof(WORD)))
			== NULL) {
		PRINT_ALLOC_ERR;
		status = 0;
		fits_close_file(fit->fptr, &status);
		if (olddata)
			free(olddata);
		return -1;
	}

	if (fit->naxis == 3) {
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data + nbdata;
		fit->pdata[BLAYER] = fit->data + nbdata * 2;
	} else {
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data;
		fit->pdata[BLAYER] = fit->data;
	}

	read_fits_header(fit);	// stores useful header data in fit

	retval = read_fits_with_convert(fit, filename);
	fit->top_down = FALSE;

	if (!retval) {
		// copy the entire header
		if (fit->header)
			free(fit->header);
		fit->header = copy_header(fit);

		basename = g_path_get_basename(filename);
		siril_log_message(_("Reading FITS: file %s, %ld layer(s), %ux%u pixels\n"),
				basename, fit->naxes[2], fit->rx, fit->ry);
		g_free(basename);
	}

	status = 0;
	fits_close_file(fit->fptr, &status);
	return retval;
}

int siril_fits_open_diskfile(fitsfile **fptr, const char *filename, int iomode, int *status) {
	gchar *localefilename = get_locale_filename(filename);
	fits_open_diskfile(fptr, localefilename, iomode, status);
	g_free(localefilename);
	return *status;
}

// reset a fit data structure, deallocates everything in it and zero the data
void clearfits(fits *fit) {
	if (fit == NULL)
		return;
	if (fit->data)
		free(fit->data);
	if (fit->header)
		free(fit->header);
	if (fit->history)
		g_slist_free_full(fit->history, free);
	if (fit->stats) {
		int i;
		for (i = 0; i < fit->naxes[2]; i++)
			free_stats(fit->stats[i]);
		free(fit->stats);
	}
	memset(fit, 0, sizeof(fits));
}

/* Read a rectangular section of a FITS image in Siril's format, pointed by its
 * exact filename. Only layer layer is read.
 * Returned fit->data is upside-down. */
int readfits_partial(const char *filename, int layer, fits *fit,
		const rectangle *area, gboolean do_photometry) {
	int status;
	unsigned int nbdata;
	double offset, data_max = 0.0;
	double mini, maxi;

	status = 0;
	if (siril_fits_open_diskfile(&(fit->fptr), filename, READONLY, &status)) {
		report_fits_error(status);
		return status;
	}

	status = 0;
	fits_get_img_param(fit->fptr, 3, &(fit->bitpix), &(fit->naxis), fit->naxes,
			&status);
	if (status) {
		report_fits_error(status);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return status;
	}

	fit_stats(fit, &mini, &maxi);

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "BZERO", &offset, NULL, &status);
	if (!status) {
		if (fit->bitpix == SHORT_IMG && offset != 0.0) {
			fit->bitpix = USHORT_IMG;
		}
		else if (fit->bitpix == LONG_IMG && offset != 0.0)
			fit->bitpix = ULONG_IMG;
	} else {
		/* but some software just put unsigned 16-bit data in the file
		 * and don't set the BZERO keyword... */
		if (status == KEY_NO_EXIST && fit->bitpix == SHORT_IMG)
			fit->bitpix = USHORT_IMG;
	}
	fit->orig_bitpix = fit->bitpix;

	if (do_photometry) {
		read_fits_date_obs_header(fit);
		status = 0;
		__tryToFindKeywords(fit->fptr, TDOUBLE, Exposure, &fit->exposure);
	}

	if (fit->naxis == 2 && fit->naxes[2] == 0)
		fit->naxes[2] = 1;	// see readfits for the explanation
	if (layer > fit->naxes[2] - 1) {
		siril_log_message(_("FITS read partial: there is no layer %d in the image %s\n"),
				layer + 1, filename);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return -1;
	}

	if (fit->naxis == 3 && fit->naxes[2] != 3) {
		siril_log_message(_("Unsupported FITS image format (%ld axes).\n"),
				fit->naxes[2]);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return -1;
	}

	status = 0;
	if (fit->bitpix == FLOAT_IMG)
		fits_read_key(fit->fptr, TDOUBLE, "DATAMAX", &data_max, NULL, &status);
	if (status == KEY_NO_EXIST) {
		data_max = maxi;
	}


	nbdata = area->w * area->h;
	/* realloc fit->data to the image size */
	WORD *olddata = fit->data;
	if ((fit->data = realloc(fit->data, nbdata * sizeof(WORD))) == NULL) {
		PRINT_ALLOC_ERR;
		status = 0;
		fits_close_file(fit->fptr, &status);
		if (olddata)
			free(olddata);
		return -1;
	}
	fit->pdata[RLAYER] = fit->data;
	fit->pdata[GLAYER] = fit->data;
	fit->pdata[BLAYER] = fit->data;

	status = internal_read_partial_fits(fit->fptr, fit->naxes[1],
			fit->bitpix, fit->data, data_max > 1.0, layer, area);

	if (status) {
		report_fits_error(status);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return -1;
	}

	fit->naxes[0] = area->w;
	fit->naxes[1] = area->h;
	fit->rx = fit->naxes[0];
	fit->ry = fit->naxes[1];
	fit->naxes[2] = 1;	
	fit->naxis = 2;

	status = 0;
	fits_close_file(fit->fptr, &status);
	fprintf(stdout, _("Loaded partial FITS file %s\n"), filename);
	return 0;
}

/* read subset of an opened fits file.
 * The rectangle's coordinates x,y start at 0,0 for first pixel in the image.
 * layer and index also start at 0.
 * buffer has to be allocated with enough space to store the area.
 */
int read_opened_fits_partial(sequence *seq, int layer, int index, WORD *buffer,
		const rectangle *area) {
	int status;

	if (!seq || !seq->fptr || !seq->fptr[index]) {
		printf("data initialization error in read fits partial\n");
		return 1;
	}
	if (area->x < 0 || area->y < 0 || area->x >= seq->rx || area->y >= seq->ry
			|| area->w <= 0 || area->h <= 0 || area->x + area->w > seq->rx
			|| area->y + area->h > seq->ry) {
		fprintf(stderr, "partial read from FITS file has been requested outside image bounds or with invalid size\n");
		return 1;
	}

#ifdef _OPENMP
	assert(seq->fd_lock);
	omp_set_lock(&seq->fd_lock[index]);
#endif

	status = internal_read_partial_fits(seq->fptr[index], seq->ry, seq->bitpix, buffer, seq->data_max > 1.0, layer, area);

#ifdef _OPENMP
	omp_unset_lock(&seq->fd_lock[index]);
#endif
	if (status)
		return 1;

	/* reverse the read data, because it's stored upside-down */
	WORD *swap = malloc(area->w * sizeof(WORD));
	int i;
	for (i = 0; i < area->h/2 ; i++) {
		memcpy(swap, buffer + i*area->w, area->w*sizeof(WORD));
		memcpy(buffer + i*area->w, buffer + (area->h - i - 1)*area->w, area->w*sizeof(WORD));
		memcpy(buffer + (area->h - i - 1)*area->w, swap, area->w*sizeof(WORD));
	}
	free(swap);

	return 0;
}

/* creates, saves and closes the file associated to f, overwriting previous  */
int savefits(const char *name, fits *f) {
	int status, i;
	long orig[3] = { 1L, 1L, 1L }, pixel_count;
	char filename[256];
	BYTE *data8;

	f->naxes[0] = f->rx;
	f->naxes[1] = f->ry;

	if (f->naxis == 3 && f->naxes[2] != 3) {
		printf("Trying to save a FITS color file with more than 3 channels?");
		return 1;
	}

	if (!ends_with(name, com.ext)) {
		snprintf(filename, 255, "%s%s", name, com.ext);
	} else {
		snprintf(filename, 255, "%s", name);
	}

	g_unlink(filename); /* Delete old file if it already exists */

	status = 0;
	if (siril_fits_create_diskfile(&(f->fptr), filename, &status)) { /* create new FITS file */
		report_fits_error(status);
		return 1;
	}
	status = 0;
	/* some float cases where it is USHORT saved as float */
	if (f->bitpix != BYTE_IMG && f->data_max > 1.0 && f->data_max <= USHRT_MAX) {
		f->bitpix = USHORT_IMG;
	}
	if (fits_create_img(f->fptr, f->bitpix, f->naxis, f->naxes, &status)) {
		report_fits_error(status);
		return 1;
	}

	pixel_count = f->naxes[0] * f->naxes[1] * f->naxes[2];

	status = 0;
	switch (f->bitpix) {
	case BYTE_IMG:
		data8 = calloc(pixel_count, sizeof(BYTE));
		WORD norm = get_normalized_value(f);
		for (i = 0; i < pixel_count; i++) {
			if (norm == USHRT_MAX)
				data8[i] = conv_to_BYTE(f->data[i]);
			else
				data8[i] = (BYTE) (f->data[i]);
		}
		if (fits_write_pix(f->fptr, TBYTE, orig, pixel_count, data8, &status)) {
			report_fits_error(status);
			if (data8)
				free(data8);
			return 1;
		}
		free(data8);
		break;
	case SHORT_IMG:
		if (f->orig_bitpix == BYTE_IMG) {
			conv_8_to_16(f->data, pixel_count);
		}
		if (fits_write_pix(f->fptr, TSHORT, orig, pixel_count, f->data, &status)) {
			report_fits_error(status);
			return 1;
		}
		break;
	case USHORT_IMG:
		if (f->orig_bitpix == BYTE_IMG) {
			conv_8_to_16(f->data, pixel_count);
		}
		if (fits_write_pix(f->fptr, TUSHORT, orig, pixel_count, f->data, &status)) {
			report_fits_error(status);
			return 1;
		}
		break;
	case LONG_IMG:
	case LONGLONG_IMG:
	case FLOAT_IMG:
	case DOUBLE_IMG:
	default:
		siril_log_message(_("ERROR: trying to save a FITS image "
				"with an unsupported format (%d).\n"), f->bitpix);
		fits_close_file(f->fptr, &status);
		return 1;
	}

	if (!status) {
		save_fits_header(f);
		// copy the entire header
		if (f->header)
			free(f->header);
		f->header = copy_header(f);
	}

	status = 0;
	fits_close_file(f->fptr, &status);
	if (!status) {
		siril_log_message(_("Saving FITS: file %s, %ld layer(s), %ux%u pixels\n"),
				filename, f->naxes[2], f->rx, f->ry);
	}
	return 0;
}

/* Duplicates some of a fits data into another, with various options; the third
 * parameter, oper, indicates with bits what operations will be done:
 *
 * - CP_ALLOC: allocates the to->data pointer to the size of from->data and
 *   sets to->pdata; required if data is not already allocated with the
 *   correct size or at all. No data is copied
 * - CP_INIT: initialize to->data with zeros, same size of the image in from,
 *   but no other data is modified. Ignored if not used with CP_ALLOC.
 * - CP_COPYA: copies the actual data, from->data to to->data on all layers,
 *   but no other information from the source. Should not be used with CP_INIT
 * - CP_FORMAT: copy all metadata and leaves data to null
 * - CP_EXTRACT: same as CP_FORMAT | CP_COPYA, but only for layer number passed
 *   as argument and sets layer number information to 1, without allocating
 *   data to one layer only
 * - CP_EXPAND: forces the destination number of layers to be taken as 3, but
 *   the other operations have no modifications, meaning that if the source
 *   image has one layer, the output image will have only one actual layer
 *   filled, and two filled with random data unless CP_INIT is used to fill it
 *   with zeros.
 *
 * Example: to duplicate a fits from one to an unknown-allocated other, those
 * flags should be used:	CP_ALLOC | CP_COPYA | CP_FORMAT
 *
 */
int copyfits(fits *from, fits *to, unsigned char oper, int layer) {
	int depth, i;
	unsigned int nbdata = from->rx * from->ry;

	if ((oper & CP_EXPAND))
		depth = 3;
	else if ((oper & CP_EXTRACT))
		depth = 1;
	else depth = from->naxes[2];

	if ((oper & CP_FORMAT)) {
		// copying metadata, not data or stats which are kept null
		memcpy(to, from, sizeof(fits));
		to->naxis = depth == 3 ? 3 : from->naxis;
		to->naxes[2] = depth;
		if (depth != from->naxes[2]) {
			to->maxi = -1.0;
		}
		to->stats = NULL;
		to->fptr = NULL;
		to->data = NULL;
		to->pdata[0] = NULL;
		to->pdata[1] = NULL;
		to->pdata[2] = NULL;
		to->header = NULL;
		to->history = NULL;
	}

	if ((oper & CP_ALLOC)) {
		// allocating to->data and assigning to->pdata
		WORD *olddata = to->data;
		if (!(to->data = realloc(to->data, nbdata * depth * sizeof(WORD)))) {
			PRINT_ALLOC_ERR;
			if (olddata)
				free(olddata);
			return -1;
		}
		to->pdata[RLAYER] = to->data;
		if (depth == 3) {
			to->pdata[GLAYER] = to->data + nbdata;
			to->pdata[BLAYER] = to->data + 2 * nbdata;
		} else {
			to->pdata[GLAYER] = to->data;
			to->pdata[BLAYER] = to->data;
		}

		if ((oper & CP_INIT)) {
			// clearing to->data allocated above
			memset(to->data, 0, nbdata * depth * sizeof(WORD));
		}
	}

	if ((oper & CP_COPYA)) {
		// copying data and stats
		memcpy(to->data, from->data, nbdata * depth * sizeof(WORD));
		if (from->stats) {
			for (i = 0; i < from->naxes[2]; i++) {
				if (from->stats[i])
					add_stats_to_fit(to, i, from->stats[i]);
			}
		} else {
			invalidate_stats_from_fit(to);
		}
	}

	if ((oper & CP_EXTRACT)) {
		to->rx = from->rx;
		to->ry = from->ry;
		to->lo = from->lo;
		to->hi = from->hi;
		to->bitpix = from->bitpix;
		to->naxis = 2;
		to->naxes[0] = from->naxes[0];
		to->naxes[1] = from->naxes[1];
		to->naxes[2] = 1;
		if (depth != from->naxes[2]) {
			to->maxi = -1.0;
		}
		memcpy(to->data, from->pdata[layer],
				nbdata * to->naxes[2] * sizeof(WORD));
		if (from->stats && from->stats[layer])
			add_stats_to_fit(to, 0, from->stats[layer]);
	}
	update_used_memory();
	return 0;
}

/* copy non-mandatory keywords from 'from' to 'to' */
int copy_fits_metadata(fits *from, fits *to) {
	to->pixel_size_x = from->pixel_size_x;
	to->pixel_size_y = from->pixel_size_y;
	to->binning_x = from->binning_x;
	to->binning_y = from->binning_y;

	to->maximum_pixel_value = from->maximum_pixel_value;

	strncpy(to->date_obs, from->date_obs, FLEN_VALUE);
	strncpy(to->date, from->date, FLEN_VALUE);
	strncpy(to->instrume, from->instrume, FLEN_VALUE);
	strncpy(to->telescop, from->telescop, FLEN_VALUE);
	strncpy(to->observer, from->observer, FLEN_VALUE);
	strncpy(to->dft.type, from->dft.type, FLEN_VALUE);
	strncpy(to->dft.ord, from->dft.ord, FLEN_VALUE);
	strncpy(to->bayer_pattern, from->bayer_pattern, FLEN_VALUE);

	to->bayer_xoffset = from->bayer_xoffset;
	to->bayer_yoffset = from->bayer_yoffset;
	to->focal_length = from->focal_length;
	to->iso_speed = from->iso_speed;
	to->exposure = from->exposure;
	to->aperture = from->aperture;
	to->ccd_temp = from->ccd_temp;
	to->cvf = from->cvf;
	to->dft.norm[0] = from->dft.norm[0];
	to->dft.norm[1] = from->dft.norm[1];
	to->dft.norm[2] = from->dft.norm[2];

	return 0;
}

int save1fits16(const char *filename, fits *fit, int layer) {
	if (layer != RLAYER) {
		int nbdata = fit->naxes[0] * fit->naxes[1];
		memcpy(fit->data, fit->data + layer * nbdata, nbdata * sizeof(WORD));
	}
	fit->naxis = 2;
	fit->naxes[2] = 1;
	return savefits(filename, fit);
}

/* this method converts 24-bit RGB or BGR data (no padding) to 48-bit FITS data.
 * order is RGB when inverted is FALSE, BGR when inverted is TRUE
 * fit->data has to be already allocated and fit->rx and fit->ry must be correct */
void rgb24bit_to_fits48bit(unsigned char *rgbbuf, fits *fit, gboolean inverted) {
	int i, j, nbdata;
	nbdata = fit->rx * fit->ry;
	WORD *rdata, *gdata, *bdata;
	rdata = fit->pdata[RLAYER] = fit->data;
	gdata = fit->pdata[GLAYER] = fit->data + nbdata;
	bdata = fit->pdata[BLAYER] = fit->data + 2 * nbdata;
	for (i = 0; i < fit->rx; ++i) {
		for (j = 0; j < fit->ry; ++j) {
			if (inverted)
				*bdata++ = (WORD) *rgbbuf++;
			else
				*rdata++ = (WORD) *rgbbuf++;
			*gdata++ = (WORD) *rgbbuf++;
			if (inverted)
				*rdata++ = (WORD) *rgbbuf++;
			else
				*bdata++ = (WORD) *rgbbuf++;
		}
	}
}

/* this method converts 8-bit gray data to 16-bit FITS data.
 * fit->data has to be already allocated and fit->rx and fit->ry must be correct */
void rgb8bit_to_fits16bit(unsigned char *graybuf, fits *fit) {
	WORD *data;
	int nbdata, i;
	fit->pdata[0] = fit->data;
	fit->pdata[1] = fit->data;
	fit->pdata[2] = fit->data;
	data = fit->data;
	nbdata = fit->rx * fit->ry;
	for (i = 0; i < nbdata; ++i) {
		*data++ = *graybuf++;
	}
}

/* this method converts 48-bit RGB or BGR data (no padding) to 48-bit FITS data.
 * order is RGB when inverted is FALSE, BGR when inverted is TRUE
 * the endianness of the data, since we have two byte per value, may not match the endianness
 * of our FITS files, so the change_endian parameter allows to flip the endian.
 * fit->data has to be already allocated and fit->rx and fit->ry must be correct */
void rgb48bit_to_fits48bit(WORD *rgbbuf, fits *fit, gboolean inverted,
		gboolean change_endian) {
	int i, j, nbdata;
	nbdata = fit->rx * fit->ry;
	WORD *rdata, *gdata, *bdata, curval;
	rdata = fit->pdata[0] = fit->data;
	gdata = fit->pdata[1] = fit->data + nbdata;
	bdata = fit->pdata[2] = fit->data + 2 * nbdata;
	for (i = 0; i < fit->rx; ++i) {
		for (j = 0; j < fit->ry; ++j) {
			curval = *rgbbuf++;
			if (change_endian)
				curval = (curval >> 8) | (curval << 8);
			if (inverted)
				*bdata++ = curval;
			else
				*rdata++ = curval;

			curval = *rgbbuf++;
			if (change_endian)
				curval = (curval >> 8) | (curval << 8);
			*gdata++ = curval;

			curval = *rgbbuf++;
			if (change_endian)
				curval = (curval >> 8) | (curval << 8);
			if (inverted)
				*rdata++ = curval;
			else
				*bdata++ = curval;
		}
	}
}

/* this method flips top-bottom of fit data.
 * fit->rx, fit->ry, fit->naxes[2] and fit->pdata[*] are required to be assigned correctly */
void fits_flip_top_to_bottom(fits *fit) {
	int line, axis, line_size;
	WORD *swapline, *src, *dst;

	line_size = fit->rx * sizeof(WORD);
	swapline = malloc(line_size);

	for (axis = 0; axis < fit->naxes[2]; axis++) {
		for (line = 0; line < fit->ry / 2; line++) {
			src = fit->pdata[axis] + line * fit->rx;
			dst = fit->pdata[axis] + (fit->ry - line - 1) * fit->rx;

			memcpy(swapline, src, line_size);
			memcpy(src, dst, line_size);
			memcpy(dst, swapline, line_size);
		}
	}
	free(swapline);
}

/* This function copies an area from the fits 'from' on layer 'layer' into
 * another and initializes all relevant data */
/* the crop function does the same but in place and for all channels without
 * reallocating */
void extract_region_from_fits(fits *from, int layer, fits *to,
		const rectangle *area) {
	int x, y, d, ystart, yend;
	clearfits(to);
	to->data = malloc(area->w * area->h * sizeof(WORD));

	d = 0;
	ystart = from->ry - area->y - area->h;
	yend = from->ry - area->y;
	for (y = ystart; y < yend; y++) {
		for (x = area->x; x < area->x + area->w; x++) {
			to->data[d++] = from->pdata[layer][x + y * from->rx];
		}
	}

	to->rx = area->w;
	to->ry = area->h;
	to->naxes[0] = area->w;
	to->naxes[1] = area->h;
	to->naxes[2] = 1;
	to->naxis = 2;
	to->pdata[0] = to->data;
	to->pdata[1] = to->data;
	to->pdata[2] = to->data;
	to->bitpix = (from->bitpix ? from->bitpix : USHRT_MAX);
}

/* creates a new fit image from scratch (NULL fit) or into a fits * previously
 * allocated, with non-cleared (random) data. */
int new_fit_image(fits **fit, int width, int height, int nblayer) {
	gint npixels;
	WORD *data;
	assert(width > 0);
	assert(height > 0);
	assert(nblayer == 1 || nblayer == 3);

	npixels = width * height;
	data = malloc(npixels * nblayer * sizeof(WORD));
	if (data == NULL) {
		PRINT_ALLOC_ERR;
		return -1;
	}

	if (*fit)
		clearfits(*fit);
	else {
		*fit = calloc(1, sizeof(fits));
		if (!*fit) {
			PRINT_ALLOC_ERR;
			free(data);
			return -1;
		}
	}

	(*fit)->bitpix = USHORT_IMG;
	if (nblayer > 1)
		(*fit)->naxis = 3;
	else (*fit)->naxis = 2;
	(*fit)->rx = width;
	(*fit)->ry = height;
	(*fit)->naxes[0] = width;
	(*fit)->naxes[1] = height;
	(*fit)->naxes[2] = nblayer;
	(*fit)->data = data;
	(*fit)->pdata[RLAYER] = (*fit)->data;
	if (nblayer > 1) {
		(*fit)->pdata[GLAYER] = (*fit)->data + npixels;
		(*fit)->pdata[BLAYER] = (*fit)->data + npixels * 2;
	}
	else {
		(*fit)->pdata[GLAYER] = (*fit)->data;
		(*fit)->pdata[BLAYER] = (*fit)->data;
	}
	return 0;
}

/* In-place conversion to one channel.
 * See copyfits with CP_EXTRACT for the same in a new fits */
void keep_first_channel_from_fits(fits *fit) {
	if (fit->naxis == 1)
		return;
	fit->naxis = 2;
	fit->naxes[2] = 1;
	fit->data = realloc(fit->data, fit->rx * fit->ry * sizeof(WORD));
	fit->pdata[RLAYER] = fit->data;
	fit->pdata[GLAYER] = fit->data;
	fit->pdata[BLAYER] = fit->data;
	if (fit->maxi > 0) {
		if (fit->maxi != fit_get_max(fit, 0))
			fit->maxi = 0;
		if (fit->mini != fit_get_min(fit, 0))
			fit->mini = 0;
	}
}
