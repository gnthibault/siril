/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2017 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"
#include "gui/callbacks.h"

static char *MIPSHI[] = {"MIPS-HI", "CWHITE", NULL };
static char *MIPSLO[] = {"MIPS-LO", "CBLACK", NULL };
static char *PixSizeX[] = { "XPIXSZ", "XPIXELSZ", NULL };
static char *PixSizeY[] = { "YPIXSZ", "YPIXELSZ", NULL };
static char *BinX[] = { "XBINNING", "BINX", NULL };
static char *BinY[] = { "YBINNING", "BINY", NULL };
static char *Focal[] = { "FOCAL", "FOCALLEN", NULL };
static char *Exposure[] = { "EXPTIME", "EXPOSURE", NULL };

#define __tryToFindKeywords(fptr, type, keyword, value) \
{ \
	int __iter__ = 0; \
	int __status__; \
	do { \
		__status__ = 0; \
		fits_read_key(fptr, type, keyword[__iter__], value, NULL, &__status__); \
		__iter__++; \
	} while ((keyword[__iter__]) && (__status__ > 0)); \
}

/* reading the FITS header to get useful information */
static void read_fits_header(fits *fit) {
	/* about the status argument: http://heasarc.gsfc.nasa.gov/fitsio/c/c_user/node28.html */
	int status = 0;
	int zero;

	__tryToFindKeywords(fit->fptr, TUSHORT, MIPSHI, &fit->hi);
	__tryToFindKeywords(fit->fptr, TUSHORT, MIPSLO, &fit->lo);

	status = 0;
	fits_read_key(fit->fptr, TINT, "BSCALE", &zero, NULL, &status);
	if (!status && 1 != zero)
		siril_log_message(
				_("Loaded FITS file has a BSCALE different than 1 (%d)\n"), zero);
	status = 0;
	fits_read_key(fit->fptr, TINT, "BZERO", &zero, NULL, &status);

	if (fit->bitpix == SHORT_IMG && zero == 32768)
		fit->bitpix = USHORT_IMG;

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
	fits_read_key(fit->fptr, TSTRING, "DATE-OBS", &(fit->date_obs), NULL,
			&status);

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

	/*******************************************************************
	 * ************************* DFT KEYWORDS **************************
	 * ****************************************************************/

	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "DFT_NOR0", &(fit->dft_norm[0]), NULL,
			&status);
	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "DFT_NOR1", &(fit->dft_norm[1]), NULL,
			&status);
	status = 0;
	fits_read_key(fit->fptr, TDOUBLE, "DFT_NOR2", &(fit->dft_norm[2]), NULL,
			&status);

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "DFT_ORD", &(fit->dft_ord), NULL,
			&status);

	status = 0;
	fits_read_key(fit->fptr, TSTRING, "DFT_TYPE", &(fit->dft_type), NULL,
			&status);

	status = 0;
	fits_read_key(fit->fptr, TUSHORT, "DFT_RX", &(fit->dft_rx), NULL, &status);
	fits_read_key(fit->fptr, TUSHORT, "DFT_RY", &(fit->dft_ry), NULL, &status);
}

// return 0 on success, fills realname if not NULL with the opened file's name
int readfits(const char *filename, fits *fit, char *realname) {
	int status;
	long orig[3] = { 1L, 1L, 1L };
	// orig ^ gives the coordinate in each dimension of the first pixel to be read
	int zero = 0;
	char name[256], *msg = NULL, *basename;
	image_type imagetype;
	int i;
	unsigned int nbdata;
	double *pixels_double;
	signed long *pixels_long;
	unsigned long offset;
	BYTE *data8;

	fit->naxes[2] = 1; //initialization of the axis numer before opening : NEED TO BE OPTIMIZED

	if (stat_file(filename, &imagetype, name)) {
		msg = siril_log_message(_("%s.[any_allowed_extension] not found.\n"),
				filename);
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		return 1;
	}
	if (imagetype != TYPEFITS) {
		msg = siril_log_message(
				_("The file %s is not a FITS file or doesn't exists with FITS extensions.\n"),
						filename);
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		return 1;
	}

	if (realname)
		strcpy(realname, name);
	status = 0;
	fits_open_diskfile(&(fit->fptr), name, READONLY, &status);
	if (status) {
		report_fits_error(status);
		return status;
	}

	status = 0;
	fits_get_img_param(fit->fptr, 3, &(fit->bitpix), &(fit->naxis), fit->naxes,
			&status);
	if (status) {
		msg = siril_log_message(
				_("FITSIO error getting image parameters, file %s.\n"), filename);
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		report_fits_error(status);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return status;
	}

	fit->rx = fit->naxes[0];
	fit->ry = fit->naxes[1];
	nbdata = fit->rx * fit->ry;

	if (fit->naxis == 3 && fit->naxes[2] != 3) {
		msg = siril_log_message(_("Unknown FITS image format (%ld axes).\n"),
				fit->naxes[2]);
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		status = 0;
		fits_close_file(fit->fptr, &status);
		return -1;
	}
	if (fit->naxis == 2 && fit->naxes[2] == 0) {
		fit->naxes[2] = 1;
		/* naxes[2] is set to 1 because:
		 * - it doesn't matter, since naxis is 2, it's not used
		 * - it's very convenient to use it in multiplications as the number of layers
		 */
	}
	if (fit->bitpix == LONGLONG_IMG) {
		msg = siril_log_message(
				_("FITS images with 64 bits signed integer per pixel.channel are not supported.\n"));
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		status = 0;
		fits_close_file(fit->fptr, &status);
		return -1;
	}

	/* realloc fit->data to the image size */
	WORD *olddata = fit->data;
	if ((fit->data = realloc(fit->data, nbdata * fit->naxes[2] * sizeof(WORD)))
			== NULL) {
		fprintf(stderr, "readfits: error realloc %s %lu\n", filename,
				nbdata * fit->naxes[2]);
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

	read_fits_header(fit);

	status = 0;
	switch (fit->bitpix) {
	case SBYTE_IMG:
	case BYTE_IMG:
		data8 = calloc(fit->rx * fit->ry * fit->naxes[2], sizeof(BYTE));
		fits_read_pix(fit->fptr, TBYTE, orig, nbdata * fit->naxes[2], &zero,
				data8, &zero, &status);
		for (i=0; i < fit->rx * fit->ry * fit->naxes[2]; i++)
				fit->data[i] = (WORD)data8[i];
		free(data8);
		break;
	case USHORT_IMG:
		fits_read_pix(fit->fptr, TUSHORT, orig, nbdata * fit->naxes[2], &zero,
				fit->data, &zero, &status);
		break;
	case SHORT_IMG:
		fits_read_pix(fit->fptr, TSHORT, orig, nbdata * fit->naxes[2], &zero,
				fit->data, &zero, &status);
		break;
	case ULONG_IMG:		// 32-bit unsigned integer pixels
	case LONG_IMG:		// 32-bit signed integer pixels
		fits_read_key(fit->fptr, TULONG, "BZERO", &offset, NULL, &status);
		pixels_long = (signed long*) malloc(
				nbdata * fit->naxes[2] * sizeof(signed long));
		status = 0;
		fits_read_pix(fit->fptr, TLONG, orig, nbdata * fit->naxes[2], &zero,
				pixels_long, &zero, &status);
		long m_long = gsl_stats_long_max(pixels_long, 1, nbdata * fit->naxes[2]);
		for (i = 0; i < nbdata * fit->naxes[2]; i++) {
			double shift = (double) (0x80000000UL - offset) / (double) UINT_MAX;
			if (m_long> USHRT_MAX) {
				double pixel = (double) (pixels_long[i] / (double) UINT_MAX);
				fit->data[i] = round_to_WORD((pixel + shift) * USHRT_MAX_DOUBLE);
			}
			else {
				double pixel = (double) (pixels_long[i]);
				fit->data[i] = round_to_WORD((pixel + shift));
			}
		}
		free(pixels_long);
		fit->bitpix = USHORT_IMG;
		break;
	case DOUBLE_IMG:	// 64-bit floating point pixels
	case FLOAT_IMG:		// 32-bit floating point pixels
		pixels_double = (double *) malloc(
				nbdata * fit->naxes[2] * sizeof(double));
		fits_read_pix(fit->fptr, TDOUBLE, orig, nbdata * fit->naxes[2], &zero,
				pixels_double, &zero, &status);
		/* The BITPIX seems strange for some pictures of FLOAT_IMG.
		 * Indeed, it looks like that some FLOAT_IMG are in fact integers.
		 * In order to read most pictures we test the maximum value to
		 * know if the range is [0, 1] or not */
		double m_double = gsl_stats_max(pixels_double, 1,
				nbdata * fit->naxes[2]);
		for (i = 0; i < nbdata * fit->naxes[2]; i++) {
			if (m_double > 1.0) {
				fit->data[i] = round_to_WORD(pixels_double[i]);
			} else
				fit->data[i] = round_to_WORD(
				USHRT_MAX_DOUBLE * pixels_double[i]);
		}
		free(pixels_double);
		fit->bitpix = USHORT_IMG;		// image is converted in 16-bit unsigned
		break;
	case LONGLONG_IMG:	// 64-bit integer pixels
		break;

	default:
		msg = siril_log_message(_("Unknown FITS image format.\n"));
	}
	if (msg) {
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		status = 0;
		fits_close_file(fit->fptr, &status);
		update_used_memory();
		return -1;
	}
	if (status) {
		msg = siril_log_message(_("Fitsio error reading data, file: %s.\n"),
				filename);
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		report_fits_error(status);
	}

	if (fit->header)
		free(fit->header);
	fit->header = list_header(fit);

	if (gtk_widget_get_visible(lookup_widget("data_dialog")))// update data if already shown
		show_FITS_header(fit);

	status = 0;
	fits_close_file(fit->fptr, &status);
	basename = g_path_get_basename(filename);
	siril_log_message(_("Reading FITS: file %s, %ld layer(s), %ux%u pixels\n"),
			basename, fit->naxes[2], fit->rx, fit->ry);
	g_free(basename);
	return 0;
}

char *list_header(fits *fit) {
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
	if (header == NULL)
		return NULL;
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
		fprintf(stdout, "%s", header);// don't display for all frames of a sequence
	return header;
}

// reset a fit data structure, deallocates everything in it
void clearfits(fits *fit) {
	if (fit == NULL)
		return;
	if (fit->data)
		free(fit->data);
	if (fit->header)
		free(fit->header);
	memset(fit, 0, sizeof(fits));
}

void report_fits_error(int status) {
	if (status) {
		char errmsg[FLEN_ERRMSG];
		while (fits_read_errmsg(errmsg)) {
			siril_log_message(_("FITS error: %s\n"), errmsg);
		}
	}
}

/* Read a rectangular section of a FITS image in Siril's format, pointed by its
 * exact filename. Only layer layer is read. */
int readfits_partial(const char *filename, int layer, fits *fit,
		const rectangle *area) {
	int status;
	unsigned int nbdata;
	long fpixel[3], lpixel[3], inc[3] = { 1L, 1L, 1L };
	int zero = 0;

	status = 0;
	if (fits_open_diskfile(&(fit->fptr), filename, READONLY, &status))
		report_fits_error(status);
	if (status)
		return status;

	status = 0;
	fits_get_img_param(fit->fptr, 3, &(fit->bitpix), &(fit->naxis), fit->naxes,
			&status);
	if (status) {
		report_fits_error(status);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return status;
	}
	fit->rx = fit->naxes[0];	// size of the real image
	fit->ry = fit->naxes[1];
	if (fit->naxis == 2 && fit->naxes[2] == 0)
		fit->naxes[2] = 1;	// see readfits for the explanation
	if (layer > fit->naxes[2] - 1) {
		siril_log_message(
				_("FITS read partial: there is no layer %d in the image %s\n"),
				layer + 1, filename);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return -1;
	}

	/* fpixel is first pixel, lpixel is last pixel, starts with value 1 */
	fpixel[0] = area->x + 1;	// in siril, it starts with 0
	fpixel[1] = fit->ry - area->y - area->h;
	fpixel[2] = layer + 1;
	lpixel[0] = area->x + area->w;	// with w and h at least 1, we're ok
	lpixel[1] = fit->ry - area->y - 1;
	lpixel[2] = layer + 1;
	fit->rx = area->w;		// size of the area to be read
	fit->ry = area->h;
	nbdata = fit->rx * fit->ry;

	if (fit->naxis == 3 && fit->naxes[2] != 3) {
		siril_log_message(_("Unknown FITS image format (%ld axes).\n"),
				fit->naxes[2]);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return -1;
	}
	if (fit->bitpix != SHORT_IMG && fit->bitpix != USHORT_IMG
			&& fit->bitpix != BYTE_IMG) {
		siril_log_message(
				_("Only Siril FITS images can be used with partial image reading.\n"));
		return -1;
	}
	fit->naxes[2] = 1;	// force to 1 layer
	fit->bitpix = USHORT_IMG;

	/* realloc fit->data to the image size */
	WORD *olddata = fit->data;
	if ((fit->data = realloc(fit->data, nbdata * sizeof(WORD))) == NULL) {
		fprintf(stderr, "readfits: error realloc %s %u\n", filename, nbdata);
		status = 0;
		fits_close_file(fit->fptr, &status);
		if (olddata)
			free(olddata);
		return -1;
	}
	fit->pdata[RLAYER] = fit->data;
	fit->pdata[GLAYER] = fit->data;
	fit->pdata[BLAYER] = fit->data;

	status = 0;
	fits_read_subset(fit->fptr, TUSHORT, fpixel, lpixel, inc, &zero, fit->data,
			&zero, &status);
	if (status) {
		report_fits_error(status);
		status = 0;
		fits_close_file(fit->fptr, &status);
		return -1;
	}
	status = 0;
	fits_close_file(fit->fptr, &status);
	siril_log_message(_("Loaded partial FITS file %s\n"), filename);
	return 0;
}

/* read subset of an opened fits file.
 * The rectangle's coordinates x,y start at 0,0 for first pixel in the image.
 * layer and index also start at 0.
 * buffer has to be allocated with enough space to store the area.
 */
int read_opened_fits_partial(sequence *seq, int layer, int index, WORD *buffer,
		const rectangle *area) {
	int status = 0;
	long fpixel[3], lpixel[3], inc[3] = { 1L, 1L, 1L };
	int zero = 0;

	if (!seq || !seq->fptr || !seq->fptr[index]) {
		printf("data initialization error in read fits partial\n");
		return 1;
	}
	if (area->x < 0 || area->y < 0 || area->x >= seq->rx || area->y >= seq->ry
			|| area->w <= 0 || area->h <= 0 || area->x + area->w > seq->rx
			|| area->y + area->h > seq->ry) {
		printf(
				"partial read from FITS file has been requested outside image bounds or with invalid size\n");
		return 1;
	}

	/* first and last pixel, in 1 .. rx and 1 .. ry coordinates */
	fpixel[0] = area->x + 1;	// in siril, it starts with 0
	fpixel[1] = seq->ry - area->y - area->h + 1;
	fpixel[2] = layer + 1;
	lpixel[0] = area->x + area->w;	// with w and h at least 1, we're ok
	lpixel[1] = seq->ry - area->y;
	lpixel[2] = layer + 1;

#ifdef _OPENMP
	assert(seq->fd_lock);
	omp_set_lock(&seq->fd_lock[index]);
#endif
	fits_read_subset(seq->fptr[index], TUSHORT, fpixel, lpixel, inc, &zero,
			buffer, &zero, &status);
#ifdef _OPENMP
	omp_unset_lock(&seq->fd_lock[index]);
#endif
	if (status) {
		report_fits_error(status);
		fprintf(stderr, "fpixel: %ld,%ld,%ld - lpixel: %ld,%ld,%ld\n",
				fpixel[0], fpixel[1], fpixel[2], lpixel[0], lpixel[1],
				lpixel[2]);
		return 1;
	}

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
	char filename[256], *msg;
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

	unlink(filename); /* Delete old file if it already exists */

	status = 0;
	if (fits_create_diskfile(&(f->fptr), filename, &status)) { /* create new FITS file */
		report_fits_error(status);
		return 1;
	}
	if (fits_create_img(f->fptr, f->bitpix, f->naxis, f->naxes, &status)) {
		report_fits_error(status);
		return 1;
	}

	pixel_count = f->naxes[0] * f->naxes[1] * f->naxes[2];

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
		if (fits_write_pix(f->fptr, TSHORT, orig, pixel_count, f->data, &status)) {
			report_fits_error(status);
			return 1;
		}
		break;
	case USHORT_IMG:
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
		msg = siril_log_message(
				_("ERROR: trying to save a FITS image "
				"with an unsupported format (%d).\n"), f->bitpix);
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		fits_close_file(f->fptr, &status);
		return 1;
	}

	if (!status)
		save_fits_header(f);
	fits_close_file(f->fptr, &status);
	if (!status)
		siril_log_message(_("Saving FITS: file %s, %ld layer(s), %ux%u pixels\n"),
				filename, f->naxes[2], f->rx, f->ry);
	return 0;
}

void save_fits_header(fits *fit) {
	int i, status = 0;
	int zero;
	unsigned int offset = 0;
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
		zero = 0;
		break;
	default:
	case USHORT_IMG:
		zero = 32768;
		break;
	}
	fits_update_key(fit->fptr, TUINT, "BZERO", &zero,
			"offset data range to that of unsigned short", &status);

	status = 0;
	zero = 1;
	fits_update_key(fit->fptr, TUINT, "BSCALE", &zero, "default scaling factor",
			&status);

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
		fits_update_key(fit->fptr, TUINT, "XBAYROFF", &(offset),
				"X offset of Bayer array", &status);

		status = 0;
		fits_update_key(fit->fptr, TUINT, "YBAYROFF", &(offset),
				"Y offset of Bayer array", &status);
	}

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
	if (fit->dft_type[0] != '\0') {
		if (fit->dft_type[0] == 'S')
			strcpy(comment, "Module of a Discrete Fourier Transform");
		else if (fit->dft_type[0] == 'P')
			strcpy(comment, "Phase of a Discrete Fourier Transform");
		else
			status = 1;			// should not happen
		fits_update_key(fit->fptr, TSTRING, "DFT_TYPE", &(fit->dft_type),
				comment, &status);
	}

	status = 0;
	if (fit->dft_ord[0] != '\0') {
		if (fit->dft_ord[0] == 'C')
			strcpy(comment, "Low spatial freq. are located at image center");
		else if (fit->dft_ord[0] == 'R')
			strcpy(comment, "High spatial freq. are located at image center");
		else
			status = 1;			// should not happen
		fits_update_key(fit->fptr, TSTRING, "DFT_ORD", &(fit->dft_ord), comment,
				&status);
	}

	for (i = 0; i < fit->naxes[2]; i++) {
		if (fit->dft_norm[i] > 0.) {
			char str1[] = "DFT_NOR";
			char str2[] = "Normalisation value for channel #";
			char key_str[FLEN_KEYWORD], comment_str[FLEN_VALUE];
			sprintf(key_str, "%s%d", str1, i);
			sprintf(comment_str, "%s%d", str2, i);
			status = 0;
			fits_update_key(fit->fptr, TDOUBLE, key_str, &(fit->dft_norm[i]),
					comment_str, &status);
		}
	}

	status = 0;
	if (fit->dft_rx) { /* may not be initialized */
		fits_update_key(fit->fptr, TUSHORT, "DFT_RX", &(fit->dft_rx),
				"Original width size", &status);
		fits_update_key(fit->fptr, TUSHORT, "DFT_RY", &(fit->dft_ry),
				"Original height size", &status);
	}
}

/* Duplicates some of a fits data into another, with various options; the third
 * parameter, oper, indicates with bits what operations will be done:
 *
 * - CP_ALLOC: allocates the to->data pointer to the size of from->data and
 *   sets to->pdata; required if data is not already allocated with the
 *   correct size or at all. No data is copied
 * - CP_COPYA: copies the actual data, from->data to to->data on all layers,
 *   but no other information from the source
 * - CP_INIT: initialize to->data with zeros, same size of the image in from,
 *   but no other data is modified
 * - CP_FORMAT: copy all other information than data and pdata
 * - CP_EXTRACT: same as CP_FORMAT | CP_COPYA, but only for layer number passed
 *   as argument and sets layer number information to 1, without allocating
 *   data to one layer only
 * - CP_EXPAND: forces the destination number of layers to be taken as 3, but
 *   the other operations have no modifications, meaning that if the source
 *   image has one layer, the output image will have only one actual layer
 *   filled, and two filled with random data.
 *
 * Example: to duplicate a fits from one to an unknown-allocated other, those
 * flags should be used:	CP_ALLOC | CP_COPYA | CP_FORMAT
 *
 */
int copyfits(fits *from, fits *to, unsigned char oper, int layer) {
	int depth;
	unsigned int nbdata = from->rx * from->ry;

	if ((oper & CP_EXPAND)) {
		depth = 3;
	} else {
		depth = from->naxes[2];
	}

	if ((oper & CP_ALLOC)) {
		WORD *olddata = to->data;
		if ((to->data = realloc(to->data, nbdata * depth * sizeof(WORD)))
				== NULL) {
			fprintf(stderr, "copyfits: error reallocating data\n");
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
	}
	//	memcpy(to->r,from->r,from->rx * from->ry*sizeof(WORD));

	if ((oper & CP_INIT)) {
		memset(to->data, 0, nbdata * depth * sizeof(WORD));
	}

	if ((oper & CP_COPYA)) {
		memcpy(to->data, from->data, nbdata * depth * sizeof(WORD));
	}

	if ((oper & CP_FORMAT)) {
		to->rx = from->rx;
		to->ry = from->ry;
		to->lo = from->lo;
		to->hi = from->hi;
		to->bitpix = from->bitpix;
		to->naxis = depth == 3 ? 3 : from->naxis;
		to->naxes[0] = from->naxes[0];
		to->naxes[1] = from->naxes[1];
		to->naxes[2] = depth;
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
		memcpy(to->data, from->pdata[layer],
				nbdata * to->naxes[2] * sizeof(WORD));
	}
	update_used_memory();
	return 0;
}

/* copy non-mandatory keywords from 'from' to 'to' */
int copy_header(fits *from, fits *to) {
	to->pixel_size_x = from->pixel_size_x;
	to->pixel_size_y = from->pixel_size_y;
	to->binning_x = from->binning_x;
	to->binning_y = from->binning_y;

	strncpy(to->date_obs, from->date_obs, FLEN_VALUE);
	strncpy(to->date, from->date, FLEN_VALUE);
	strncpy(to->instrume, from->instrume, FLEN_VALUE);
	strncpy(to->dft_type, from->dft_type, FLEN_VALUE);
	strncpy(to->dft_ord, from->dft_ord, FLEN_VALUE);

	to->focal_length = from->focal_length;
	to->iso_speed = from->iso_speed;
	to->exposure = from->exposure;
	to->aperture = from->aperture;
	to->ccd_temp = from->ccd_temp;
	to->dft_norm[0] = from->dft_norm[0];
	to->dft_norm[1] = from->dft_norm[1];
	to->dft_norm[2] = from->dft_norm[2];
	to->dft_rx = from->dft_rx;
	to->dft_ry = from->dft_ry;

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

int new_fit_image(fits *fit, int width, int height, int nblayer) {
	gint npixels;
	WORD *data;
	assert(width > 0);
	assert(height > 0);
	assert(nblayer <= 3);

	npixels = width * height;
	data = calloc(npixels, sizeof(WORD) * nblayer);

	if (data != NULL) {
		clearfits(fit);
		fit->bitpix = USHORT_IMG;
		if (nblayer > 1)
			fit->naxis = nblayer;
		else
			fit->naxis = 2;
		fit->rx = width;
		fit->ry = height;
		fit->naxes[0] = width;
		fit->naxes[1] = height;
		fit->naxes[2] = nblayer;
		fit->data = data;
		fit->pdata[RLAYER] = fit->data;
		if (nblayer > 1) {
			fit->pdata[GLAYER] = fit->data + npixels;
			fit->pdata[BLAYER] = fit->data + npixels * 2;
		}
		else {
			fit->pdata[GLAYER] = fit->data;
			fit->pdata[BLAYER] = fit->data;
		}
		return 0;
	} else
		return 1;
}

void keep_first_channel_from_fits(fits *fit) {
	if (fit->naxis == 1)
		return;
	fit->naxis = 1;
	fit->naxes[2] = 1;
	fit->data = realloc(fit->data, fit->rx * fit->ry * sizeof(WORD));
	fit->pdata[RLAYER] = fit->data;
	fit->pdata[GLAYER] = fit->data;
	fit->pdata[BLAYER] = fit->data;
	// mini and maxi could be wrong now
}
