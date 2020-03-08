/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef HAVE_LIBTIFF
#define uint64 uint64_hack_
#define int64 int64_hack_
#include <tiffio.h>
#undef uint64
#undef int64
#endif
#ifdef HAVE_LIBJPEG
#include <jpeglib.h>
#endif
#ifdef HAVE_LIBPNG
#include <png.h>
#include <setjmp.h>
#endif
#ifdef HAVE_LIBRAW
#include <libraw/libraw.h>
#endif
#ifdef HAVE_LIBHEIF
#include <libheif/heif.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "algos/geometry.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "single_image.h"

/********************* TIFF IMPORT AND EXPORT *********************/

#ifdef HAVE_LIBTIFF

static int readtifstrip(TIFF* tif, uint32 width, uint32 height, uint16 nsamples, WORD **data) {
	uint32 rowsperstrip;
	uint16 config;
	int retval = nsamples;

	TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
	TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);

	const unsigned long npixels = width * height;
	*data = malloc(npixels * sizeof(WORD) * nsamples);
	WORD *gbuf[3] = {*data, *data, *data};
	if (nsamples == 4) {
		siril_log_message(_("Alpha channel is ignored.\n"));
	}
	if ((nsamples == 3) || (nsamples == 4)) {
		gbuf[GLAYER] = *data + npixels;
		gbuf[BLAYER] = *data + npixels * 2;
	}

	const tmsize_t scanline = TIFFScanlineSize(tif);
	WORD *buf = (WORD *)_TIFFmalloc(TIFFStripSize(tif));
	if (!buf) {
		PRINT_ALLOC_ERR;
		return -1;
	}
	for (uint32 row = 0; row < height; row += rowsperstrip){
		uint32 nrow = (row + rowsperstrip > height ? height - row : rowsperstrip);
		switch (config) {
		case PLANARCONFIG_CONTIG:
			if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, row, 0), buf, nrow * scanline) < 0) {
				siril_log_message(_("An unexpected error was encountered while trying to read the file.\n"));
				retval = -1;
				break;
			}
			for (int i = 0; i < width * nrow; i++) {
				*gbuf[RLAYER]++ = buf[i * nsamples + 0];
				if ((nsamples == 3) || (nsamples == 4)) {
					*gbuf[GLAYER]++ = buf[i * nsamples + 1];
					*gbuf[BLAYER]++ = buf[i * nsamples + 2];
				}
			}
			break;
		case PLANARCONFIG_SEPARATE:
			if (nsamples >= 3)		//don't need to read the alpha
				nsamples = 3;
			for (int j = 0; j < nsamples; j++) {	//loop on the layer
				if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, row, j), buf, nrow * scanline) < 0) {
					siril_log_message(_("An unexpected error was encountered while trying to read the file.\n"));
					retval = -1;
					break;
				}
				for (int i = 0; i < width * nrow; i++)
					*gbuf[j]++ = buf[i];
			}
			break;
		default:
			siril_log_message(_("Unknown TIFF file.\n"));
			retval = -1;
		}
	}
	_TIFFfree(buf);
	return retval;
}

static int readtifstrip32(TIFF* tif, uint32 width, uint32 height, uint16 nsamples, float **data) {
	uint32 rowsperstrip;
	uint16 config;
	int retval = nsamples;

	TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
	TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);

	const unsigned long npixels = width * height;
	*data = malloc(npixels * sizeof(float) * nsamples);
	float *gbuf[3] = { *data, *data, *data };
	if (nsamples == 4) {
		siril_log_message(_("Alpha channel is ignored.\n"));
	}
	if ((nsamples == 3) || (nsamples == 4)) {
		gbuf[1] = *data + npixels;
		gbuf[2] = *data + npixels * 2;
	}

	const tmsize_t scanline = TIFFScanlineSize(tif);
	float *buf = (float *)_TIFFmalloc(TIFFStripSize(tif));
	if (!buf) {
		PRINT_ALLOC_ERR;
		return -1;
	}
	for (uint32 row = 0; row < height; row += rowsperstrip) {
		uint32 nrow =
				(row + rowsperstrip > height ? height - row : rowsperstrip);
		switch (config) {
		case PLANARCONFIG_CONTIG:
			if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, row, 0), buf, nrow * scanline) < 0) {
				siril_log_message(_("An unexpected error was encountered while trying to read the file.\n"));
				retval = -1;
				break;
			}
			for (int i = 0; i < width * nrow; i++) {
				*gbuf[RLAYER]++ = buf[i * nsamples + 0];
				if ((nsamples == 3) || (nsamples == 4)) {
					*gbuf[GLAYER]++ = buf[i * nsamples + 1];
					*gbuf[BLAYER]++ = buf[i * nsamples + 2];
				}
			}
			break;
		case PLANARCONFIG_SEPARATE:
			if (nsamples >= 3)		//don't need to read the alpha
				nsamples = 3;
			for (int j = 0; j < nsamples; j++) {	//loop on the layer
				if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, row, j),
						buf, nrow * scanline) < 0) {
					siril_log_message(_("An unexpected error was encountered while trying to read the file.\n"));
					retval = -1;
					break;
				}
				for (int i = 0; i < width * nrow; i++)
					*gbuf[j]++ = buf[i];
			}
			break;
		default:
			siril_log_message(_("Unknown TIFF file.\n"));
			retval = -1;
		}
	}
	_TIFFfree(buf);
	return retval;
}

static int readtif8bits(TIFF* tif, uint32 width, uint32 height, uint16 nsamples, WORD **data) {
	int retval = nsamples;

	const unsigned long npixels = width * height;
	*data = malloc(npixels * sizeof(WORD) * nsamples);
	WORD *gbuf[3] = { *data, *data, *data };
	if (nsamples == 4) {
		siril_log_message(_("Alpha channel is ignored.\n"));
	}
	if ((nsamples == 3) || (nsamples == 4)) {
		gbuf[1] = *data + npixels;
		gbuf[2] = *data + npixels * 2;
	}

	/* get the data */
	uint32 *raster = (uint32*) _TIFFmalloc(npixels * sizeof(uint32));
	if (raster != NULL) {
		if (TIFFReadRGBAImage(tif, width, height, raster, 0)) {
			for (int j = 0; j < height; j++) {
				int istart = j * width;
				for (int i = 0; i < width; i++) {
					*gbuf[RLAYER]++ = (WORD)TIFFGetR(raster[istart + i]);
					if ((nsamples == 3) || (nsamples == 4)) {
						*gbuf[GLAYER]++ = (WORD)TIFFGetG(raster[istart + i]);
						*gbuf[BLAYER]++ = (WORD)TIFFGetB(raster[istart + i]);
					}
				}
			}
		}
		else {
			siril_log_message(_("An unexpected error was encountered while trying to read the file.\n"));
			retval = -1;
		}
		_TIFFfree(raster);
	}
	else retval = -1;
	return retval;
}

static uint16_t get_compression_mode() {
	GtkToggleButton *button;

	button = GTK_TOGGLE_BUTTON(lookup_widget("radiobuttonCompDeflate"));
	if (gtk_toggle_button_get_active(button))
		return (uint16_t) COMPRESSION_ADOBE_DEFLATE;
	else
		return (uint16_t) COMPRESSION_NONE;
}

static TIFF* Siril_TIFFOpen(const char *name, const char *mode) {
#ifdef _WIN32
	wchar_t *wname;

	wname = g_utf8_to_utf16(name, -1, NULL, NULL, NULL);
	if (wname == NULL) {
		return NULL;
	}

	TIFF* tif = TIFFOpenW(wname, mode);
	g_free(wname);
	return tif;
#else
	return(TIFFOpen(name, mode));
#endif
}

/* reads a TIFF file and stores it in the fits argument.
 * If file loading fails, the argument is untouched.
 */
int readtif(const char *name, fits *fit) {
	int retval = 0;
	uint32 height, width;
	uint16 nbits, nsamples, color;
	WORD *data = NULL;
	float *fdata = NULL;
	
	TIFF* tif = Siril_TIFFOpen(name, "r");

	if (!tif) {
		siril_log_message(_("Could not open the TIFF file %s\n"), name);
		return -1;
	}
	
	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
	TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);	
	TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &nbits);
	TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &color);
	TIFFGetField(tif, TIFFTAG_MINSAMPLEVALUE, &(fit->lo));
	TIFFGetField(tif, TIFFTAG_MAXSAMPLEVALUE, &(fit->hi));

	const unsigned long npixels = width * height;

	switch(nbits){
		case 8:
			/* High level functions in readtif8bits: should read every 8-bit TIFF file */
			retval = readtif8bits(tif, width, height, nsamples, &data);
			break;

		case 16:
			retval = readtifstrip(tif, width, height, nsamples, &data);
			break;

		case 32:
			retval = readtifstrip32(tif, width, height, nsamples, &fdata);
			break;

		default :
			siril_log_message(_("Siril cannot read this TIFF format.\n"));
			retval = -1;
	}
	TIFFClose(tif);
	if (retval < 0) {
		free(data);
		free(fdata);
		return -1;
	}
	clearfits(fit);
	fit->rx = width;
	fit->ry = height;
	fit->naxes[0] = width;
	fit->naxes[1] = height;
	fit->data = data;
	fit->fdata = fdata;
	fit->binning_x = fit->binning_y = 1;
	if (nsamples == 1 || nsamples == 2) {
		fit->naxes[2] = 1;
		fit->naxis = 2;
		if (data) {
			fit->pdata[RLAYER] = fit->data;
			fit->pdata[GLAYER] = fit->data;
			fit->pdata[BLAYER] = fit->data;
		} else {
			fit->fpdata[RLAYER] = fit->fdata;
			fit->fpdata[GLAYER] = fit->fdata;
			fit->fpdata[BLAYER] = fit->fdata;
		}
	} else {
		fit->naxes[2] = 3;
		fit->naxis = 3;
		if (data) {
			fit->pdata[RLAYER] = fit->data;
			fit->pdata[GLAYER] = fit->data + npixels;
			fit->pdata[BLAYER] = fit->data + npixels * 2;
		} else {
			fit->fpdata[RLAYER] = fit->fdata;
			fit->fpdata[GLAYER] = fit->fdata + npixels;
			fit->fpdata[BLAYER] = fit->fdata + npixels * 2;
		}
	}
	switch (nbits) {
	case 8:
		fit->bitpix = BYTE_IMG;
		fit->type = DATA_USHORT;
		break;
	case 16:
		fit->bitpix = USHORT_IMG;
		fit->type = DATA_USHORT;
		mirrorx(fit, FALSE);
		break;
	case 32:
		fit->bitpix = FLOAT_IMG;
		fit->type = DATA_FLOAT;
		mirrorx(fit, FALSE);
	}
	fit->orig_bitpix = fit->bitpix;
	retval = nsamples;

	gchar *basename = g_path_get_basename(name);
	siril_log_message(_("Reading TIFF: %d-bit file %s, %ld layer(s), %ux%u pixels\n"),
						nbits, basename, fit->naxes[2], fit->rx, fit->ry);
	g_free(basename);

	return retval;
}

static void get_tif_data_from_ui(gchar **description, gchar **copyright) {
	if (!com.script && !com.headless) {
		/*******************************************************************
		 * If the user saves a tif from the graphical menu, he can set
		 * the Description and the Copyright of the Image
		 ******************************************************************/
		GtkTextView *description_txt_view = GTK_TEXT_VIEW(lookup_widget("Description_txt"));
		GtkTextView *copyright_txt_view = GTK_TEXT_VIEW(lookup_widget("Copyright_txt"));
		GtkTextBuffer *desbuf = gtk_text_view_get_buffer(description_txt_view);
		GtkTextBuffer *copybuf = gtk_text_view_get_buffer(copyright_txt_view);
		GtkTextIter itDebut;
		GtkTextIter itFin;

		gtk_text_buffer_get_start_iter(desbuf, &itDebut);
		gtk_text_buffer_get_end_iter(desbuf, &itFin);
		*description = gtk_text_buffer_get_text(desbuf, &itDebut, &itFin, TRUE);
		gtk_text_buffer_get_bounds(desbuf, &itDebut, &itFin);
		gtk_text_buffer_delete(desbuf, &itDebut, &itFin);
		gtk_text_buffer_get_start_iter(copybuf, &itDebut);
		gtk_text_buffer_get_end_iter(copybuf, &itFin);
		*copyright = gtk_text_buffer_get_text(copybuf, &itDebut, &itFin, TRUE);
		gtk_text_buffer_get_bounds(copybuf, &itDebut, &itFin);
		gtk_text_buffer_delete(copybuf, &itDebut, &itFin);
	}
}

/*** This function save the current image into a uncompressed 8- or 16-bit file *************/

int savetif(const char *name, fits *fit, uint16 bitspersample){
	int retval = 0;
	gchar *description, *copyright;

	gchar *filename = g_strdup(name);
	if (!ends_with(filename, ".tif") && (!ends_with(filename, ".tiff"))) {
		filename = str_append(&filename, ".tif");
	}

	TIFF* tif = Siril_TIFFOpen(filename, "w");

	if (tif == NULL) {
		siril_log_message(_("Siril cannot create TIFF file.\n"));
		free(filename);
		return 1;
	}
	const uint16 nsamples = (uint16) fit->naxes[2];
	const uint32 width = (uint32) fit->rx;
	const uint32 height = (uint32) fit->ry;
	
	get_tif_data_from_ui(&description, &copyright);

	/*******************************************************************/

	/* TIFF TAG FIELD */
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bitspersample);
	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
	TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tif, -1));
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, nsamples);
	TIFFSetField(tif, TIFFTAG_COMPRESSION, get_compression_mode());
	TIFFSetField(tif, TIFFTAG_IMAGEDESCRIPTION, description);
	TIFFSetField(tif, TIFFTAG_COPYRIGHT, copyright);
	TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(tif, TIFFTAG_MINSAMPLEVALUE, fit->mini);
	TIFFSetField(tif, TIFFTAG_MAXSAMPLEVALUE, fit->maxi);
	TIFFSetField(tif, TIFFTAG_SOFTWARE, PACKAGE " v" VERSION);

	if (nsamples == 1)
		TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	else if (nsamples == 3)
		TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
	else {
		TIFFClose(tif);
		siril_log_message(_("TIFF file has unexpected number of channels (not 1 or 3).\n"));
		free(filename);
		return 1;
	}

	WORD *gbuf[3] =	{ fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	float *gbuff[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER], fit->fpdata[BLAYER] };
	float norm = fit->orig_bitpix != BYTE_IMG ? UCHAR_MAX_SINGLE / USHRT_MAX_SINGLE : 1.f;

	mirrorx(fit, FALSE);

	switch (bitspersample) {
	case 8:
		siril_debug_print("Saving 8-bit TIFF file.\n");
		BYTE *buf8 = _TIFFmalloc(width * sizeof(unsigned char) * nsamples);

		for (uint32 row = 0; row < height; row++) {
			for (uint32 col = 0; col < width; col++) {
				for (uint16 n = 0; n < nsamples; n++) {
					buf8[col * nsamples + n] =
							(fit->type == DATA_USHORT) ?
									*gbuf[n]++ * norm :
									float_to_uchar_range(*gbuff[n]++);
				}
			}
			if (TIFFWriteScanline(tif, buf8, row, 0) < 0) {
				siril_debug_print("Error while writing in TIFF File.\n");
				retval = -1;
				break;
			}
		}
		_TIFFfree(buf8);
		break;
	case 16:
		siril_debug_print("Saving 16-bit TIFF file.\n");
		WORD *buf16 = _TIFFmalloc(width * sizeof(WORD) * nsamples);

		for (uint32 row = 0; row < height; row++) {
			for (uint32 col = 0; col < width; col++) {
				for (uint16 n = 0; n < nsamples; n++) {
					buf16[col * nsamples + n] =
							(fit->type == DATA_USHORT) ?
									*gbuf[n]++ * norm :
									float_to_ushort_range(*gbuff[n]++);
				}
			}
			if (TIFFWriteScanline(tif, buf16, row, 0) < 0) {
				siril_debug_print("Error while writing in TIFF File.\n");
				retval = -1;
				break;
			}
		}
		_TIFFfree(buf16);
		break;
	default:		// Should not happen
		retval = 1;
	}
	TIFFClose(tif);
	mirrorx(fit, FALSE);
	siril_log_message(_("Saving TIFF: %d-bit file %s, %ld layer(s), %ux%u pixels\n"),
			bitspersample, filename, nsamples, width, height);
	g_free(filename);
	return retval;
}
#endif	// HAVE_LIBTIFF


/********************* JPEG IMPORT AND EXPORT *********************/

#ifdef HAVE_LIBJPEG
int readjpg(const char* name, fits *fit){
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;

	FILE *f = g_fopen(name, "rb");
	if (f == NULL) {
		siril_log_message(_("Sorry but Siril cannot open the file: %s.\n"),	name);
		return -1;
	}
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, f);
	jpeg_read_header(&cinfo, TRUE);
	jpeg_start_decompress(&cinfo);

	const unsigned long npixels = cinfo.output_width * cinfo.output_height;
	WORD *data = malloc(npixels * sizeof(WORD) * 3);
	WORD *buf[3] = { data, data + npixels, data + npixels * 2 };
	int row_stride = cinfo.output_width * cinfo.output_components;
	JSAMPARRAY pJpegBuffer = (*cinfo.mem->alloc_sarray)((j_common_ptr) &cinfo, JPOOL_IMAGE,	row_stride, 1);

	while (cinfo.output_scanline < cinfo.output_height) {
		jpeg_read_scanlines(&cinfo, pJpegBuffer, 1);
		for (int i = 0; i < cinfo.output_width; i++) {
			*buf[RLAYER]++ = pJpegBuffer[0][cinfo.output_components * i + 0];
			*buf[GLAYER]++ = pJpegBuffer[0][cinfo.output_components * i + 1];
			*buf[BLAYER]++ = pJpegBuffer[0][cinfo.output_components * i + 2];
		}
	}

	fclose(f);
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	if (data != NULL) {
		clearfits(fit);
		fit->bitpix = fit->orig_bitpix = BYTE_IMG;
		if (cinfo.output_components == 1)
			fit->naxis = 2;
		else
			fit->naxis = 3;
		fit->rx = cinfo.output_width;
		fit->ry = cinfo.output_height;
		fit->naxes[0] = cinfo.output_width;
		fit->naxes[1] = cinfo.output_height;
		fit->naxes[2] = cinfo.output_components;
		fit->data = data;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data + npixels;
		fit->pdata[BLAYER] = fit->data + npixels * 2;
		fit->binning_x = fit->binning_y = 1;
		fit->type = DATA_USHORT;
	}
	mirrorx(fit, FALSE);
	gchar *basename = g_path_get_basename(name);
	siril_log_message(_("Reading JPG: file %s, %ld layer(s), %ux%u pixels\n"),
			basename, fit->naxes[2], fit->rx, fit->ry);
	g_free(basename);

	return cinfo.output_components;
}

int savejpg(const char *name, fits *fit, int quality){
	struct jpeg_compress_struct cinfo;    // Basic info for JPEG properties.
	struct jpeg_error_mgr jerr;           // In case of error.

	//## ALLOCATE AND INITIALIZE JPEG COMPRESSION OBJECT
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);

	char *filename = strdup(name);
	if (!ends_with(filename, ".jpg") && (!ends_with(filename, ".jpeg"))) {
		filename = str_append(&filename, ".jpg");
	}

	//## OPEN FILE FOR DATA DESTINATION:
	FILE *f = g_fopen(filename, "wb");
	if (f == NULL) {
		siril_log_message(_("Siril cannot create JPG file.\n"));
		free(filename);
		return 1;
	}
	jpeg_stdio_dest(&cinfo, f);

	//## SET PARAMETERS FOR COMPRESSION:
	cinfo.image_width  = fit->rx;   // |-- Image width and height in pixels.
	cinfo.image_height = fit->ry;   // |
	cinfo.input_components = fit->naxes[2];     // Number of color components per pixel.
	cinfo.in_color_space = (fit->naxes[2] == 3) ? JCS_RGB : JCS_GRAYSCALE; // Colorspace of input image as RGB.

	WORD *gbuf[3] =	{ fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	float *gbuff[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER], fit->fpdata[BLAYER] };

	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, quality, TRUE);

	//## CREATE IMAGE BUFFER TO WRITE FROM AND MODIFY THE IMAGE TO LOOK LIKE CHECKERBOARD:
	unsigned char *image_buffer = (unsigned char*) malloc(
			cinfo.image_width * cinfo.image_height * cinfo.num_components);

	float norm = (fit->orig_bitpix != BYTE_IMG ?
			UCHAR_MAX_SINGLE / USHRT_MAX_SINGLE : 1.f);

	for (int i = (cinfo.image_height - 1); i >= 0; i--) {
		for (int j = 0; j < cinfo.image_width; j++) {
			int pixelIdx = ((i * cinfo.image_width) + j) * cinfo.input_components;
			if (fit->type == DATA_USHORT) {
				WORD red = *gbuf[RLAYER]++;
				image_buffer[pixelIdx + 0] = round_to_BYTE(red * norm); // r |-- Set r,g,b components to
				if (cinfo.input_components == 3) {
					WORD green = *gbuf[GLAYER]++;
					WORD blue = *gbuf[BLAYER]++;
					image_buffer[pixelIdx + 1] = round_to_BYTE(green * norm); // g |   make this pixel
					image_buffer[pixelIdx + 2] = round_to_BYTE(blue * norm); // b |
				}
			} else {
				float red = *gbuff[RLAYER]++;
				image_buffer[pixelIdx + 0] = float_to_uchar_range(red); // r |-- Set r,g,b components to
				if (cinfo.input_components == 3) {
					float green = *gbuff[GLAYER]++;
					float blue = *gbuff[BLAYER]++;
					image_buffer[pixelIdx + 1] = float_to_uchar_range(green); // g |   make this pixel
					image_buffer[pixelIdx + 2] = float_to_uchar_range(blue); // b |
				}
			}
		}
	}
	//## START COMPRESSION:
	jpeg_start_compress(&cinfo, TRUE);
	int row_stride = cinfo.image_width * cinfo.input_components;        // JSAMPLEs per row in image_buffer

	JSAMPROW row_pointer[1];
	while (cinfo.next_scanline < cinfo.image_height) {
		row_pointer[0] = &image_buffer[cinfo.next_scanline * row_stride];
		(void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}
	// NOTE: jpeg_write_scanlines expects an array of pointers to scanlines.
	//       Here the array is only one element long, but you could pass
	//       more than one scanline at a time if that's more convenient.

	//## FINISH COMPRESSION AND CLOSE FILE:
	jpeg_finish_compress(&cinfo);
	fclose(f);
	jpeg_destroy_compress(&cinfo);
	free(image_buffer);
	siril_log_message(_("Saving JPG: file %s, quality=%d%%, %ld layer(s), %ux%u pixels\n"),
						filename, quality, fit->naxes[2], fit->rx, fit->ry);
	free(filename);
	return 0;
}
#endif	// HAVE_LIBJPEG


/********************* PNG IMPORT *********************/

#ifdef HAVE_LIBPNG
/* reads a PNG file and stores it in the fits argument.
 */
int readpng(const char *name, fits* fit) {
	png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL,	NULL);
	if (!png) {
		siril_log_message(_("Sorry but Siril cannot open the file: %s.\n"),	name);
		return -1;
	}

	png_infop info = png_create_info_struct(png);
	if (!info)
		return -1;

	if (setjmp(png_jmpbuf(png)))
		return -1;

	FILE *f = g_fopen(name, "rb");
	png_init_io(png, f);

	png_read_info(png, info);

	const int width = png_get_image_width(png, info);
	const int height = png_get_image_height(png, info);
	const unsigned long npixels = (unsigned long) width * (unsigned long) height;
	png_byte color_type = png_get_color_type(png, info);
	png_byte bit_depth = png_get_bit_depth(png, info);

	WORD *data = malloc(npixels * sizeof(WORD) * 3);
	WORD *buf[3] = { data, data + npixels, data + npixels * 2 };

	if (color_type == PNG_COLOR_TYPE_PALETTE)
		png_set_palette_to_rgb(png);

	// PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
	if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
		png_set_expand_gray_1_2_4_to_8(png);

	if (png_get_valid(png, info, PNG_INFO_tRNS))
		png_set_tRNS_to_alpha(png);

	// These color_type don't have an alpha channel then fill it with 0xff.
	if (color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_GRAY
			|| color_type == PNG_COLOR_TYPE_PALETTE)
		png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

	if (color_type == PNG_COLOR_TYPE_GRAY
			|| color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
		png_set_gray_to_rgb(png);

	png_read_update_info(png, info);

	png_bytep *row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
	for (int y = 0; y < height; y++) {
		row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png, info));
	}

	png_read_image(png, row_pointers);

	fclose(f);

	if (bit_depth == 16) {			//in 16-bit: it is stored as RRGGBB
		for (int y = height - 1; y > -1; y--) {
			png_byte* row = row_pointers[y];
			for (int x = 0; x < width; x++) {
				png_byte* ptr = &(row[x * 8]);
				*buf[RLAYER]++ = ptr[0] * (UCHAR_MAX + 1) + ptr[1];
				*buf[GLAYER]++ = ptr[2] * (UCHAR_MAX + 1) + ptr[3];
				*buf[BLAYER]++ = ptr[4] * (UCHAR_MAX + 1) + ptr[5];
			}
		}
	} else {
		for (int y = height - 1; y > -1; y--) {
			png_byte* row = row_pointers[y];
			for (int x = 0; x < width; x++) {
				png_byte* ptr = &(row[x * 4]);
				*buf[RLAYER]++ = ptr[0];
				*buf[GLAYER]++ = ptr[1];
				*buf[BLAYER]++ = ptr[2];
			}
		}
	}
	// We define the number of channel we have
	int nbplanes;
	switch (color_type) {
	case PNG_COLOR_TYPE_RGB_ALPHA:
	case PNG_COLOR_TYPE_RGB:
		nbplanes = 3;
		break;
	case PNG_COLOR_TYPE_GRAY_ALPHA:
	case PNG_COLOR_TYPE_GRAY:
		nbplanes = 1;
		break;
	default:
		nbplanes = 0;
	}

	// free allocated memory
	for (int y = 0; y < height; y++)
		free(row_pointers[y]);
	free(row_pointers);

	if (data != NULL) {
		clearfits(fit);
		fit->rx = width;
		fit->ry = height;
		fit->naxes[0] = width;
		fit->naxes[1] = height;
		fit->naxes[2] = nbplanes;
		if (nbplanes == 1)
			fit->naxis = 2;
		else
			fit->naxis = 3;
		fit->bitpix = (bit_depth == 16) ? USHORT_IMG : BYTE_IMG;
		fit->type = DATA_USHORT;
		fit->orig_bitpix = fit->bitpix;
		fit->data = data;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data + npixels;
		fit->pdata[BLAYER] = fit->data + npixels * 2;
		fit->binning_x = fit->binning_y = 1;

	}
	gchar *basename = g_path_get_basename(name);
	siril_log_message(_("Reading PNG: %d-bit file %s, %ld layer(s), %ux%u pixels\n"),
			bit_depth, basename, fit->naxes[2], fit->rx, fit->ry);
	g_free(basename);

	return nbplanes;
}

/* Code borrowed from SER Player */

// ------------------------------------------
// Save PNG image (colour)
// ------------------------------------------
static int32_t save_colour_file(const char *filename,
		const void *p_image_data, uint32_t width, uint32_t height,
		uint32_t bytes_per_sample) {
	int32_t ret = -1;
	png_image image; // The control structure used by libpng

	// Initialize the 'png_image' structure.
	memset(&image, 0, (sizeof image));
	image.version = PNG_IMAGE_VERSION;
	image.width = width;
	image.height = height;
	image.format = PNG_FORMAT_RGB;
	if (bytes_per_sample == 2) {
		image.format |= PNG_FORMAT_FLAG_LINEAR;
	}

	FILE *p_png_file = g_fopen(filename, "wb");
	int row_stride = image.width * -3;

	if (p_png_file != NULL) {
		png_image_write_to_stdio(&image, p_png_file, 0,
				(png_bytep) (p_image_data), row_stride,  // row_stride
				NULL);  // colormap

		ret = 0;
		fclose(p_png_file);
	}

	return ret;
}

// ------------------------------------------
// Save PNG image (mono B, G or R)
// ------------------------------------------
static int32_t save_mono_file(const char *filename, const void *p_image_data,
		uint32_t width, uint32_t height, uint32_t bytes_per_sample) {
	int32_t ret = -1;
	png_image image; // The control structure used by libpng

	// Initialize the 'png_image' structure.
	memset(&image, 0, (sizeof image));
	image.version = PNG_IMAGE_VERSION;
	image.width = width;
	image.height = height;

	if (bytes_per_sample == 1) {
		image.format = PNG_FORMAT_GRAY;  // 8-bit data
	} else {
		image.format = PNG_FORMAT_LINEAR_Y;  // 16-bit data
	}

	FILE *p_png_file = g_fopen(filename, "wb");
	int row_stride = image.width * -1;

	if (p_png_file != NULL) {
		png_image_write_to_stdio(&image, p_png_file, 0,  // convert_to_8bit
				(png_bytep) p_image_data, row_stride,  // row_stride
				NULL);  // colormap
		ret = 0;
		fclose(p_png_file);
	}

	return ret;
}

static WORD *convert_data(fits *image) {
	int ndata = image->rx * image->ry;
	int ch = image->naxes[2];

	WORD *buffer = malloc(ndata * ch * sizeof(WORD));
	for (int i = 0, j = 0; i < ndata * ch; i += ch, j++) {
		if (image->type == DATA_USHORT) {
			buffer[i + 0] = image->pdata[RLAYER][j];
			if (ch > 1) {
				buffer[i + 1] = image->pdata[GLAYER][j];
				buffer[i + 2] = image->pdata[BLAYER][j];
			}
		} else {
			buffer[i + 0] = float_to_ushort_range(image->fpdata[RLAYER][j]);
			if (ch > 1) {
				buffer[i + 1] = float_to_ushort_range(image->fpdata[GLAYER][j]);
				buffer[i + 2] = float_to_ushort_range(image->fpdata[BLAYER][j]);
			}
		}
	}
	return buffer;
}

static uint8_t *convert_data8(fits *image) {
	int ndata = image->rx * image->ry;
	const long ch = image->naxes[2];

	uint8_t *buffer = malloc(ndata * ch * sizeof(uint8_t));
	for (int i = 0, j = 0; i < ndata * ch; i += ch, j++) {
		if (image->type == DATA_USHORT) {
			buffer[i + 0] = (uint8_t) image->pdata[RLAYER][j];
			if (ch > 1) {
				buffer[i + 1] = (uint8_t) image->pdata[GLAYER][j];
				buffer[i + 2] = (uint8_t) image->pdata[BLAYER][j];
			}
		} else {
			buffer[i + 0] = float_to_uchar_range(image->fpdata[RLAYER][j]);
			if (ch > 1) {
				buffer[i + 1] = float_to_uchar_range(image->fpdata[GLAYER][j]);
				buffer[i + 2] = float_to_uchar_range(image->fpdata[BLAYER][j]);
			}
		}
	}
	return buffer;
}

int savepng(const char *name, fits *fit, uint32_t bytes_per_sample,
		gboolean is_colour) {
	int ret = -1;
	const uint32_t width = fit->rx;
	const uint32_t height = fit->ry;

	char *filename = strdup(name);
	if (!ends_with(filename, ".png")) {
		filename = str_append(&filename, ".png");
	}

	if (is_colour) {
		// Create colour PNG file
		if (bytes_per_sample == 2) {
			WORD *data = convert_data(fit);
			ret = save_colour_file(filename, data, width, height,
					bytes_per_sample);
			free(data);
		} else {
			uint8_t *data = convert_data8(fit);
			ret = save_colour_file(filename, data, width, height,
					bytes_per_sample);
			free(data);
		}
	} else {
		// Create monochrome PNG file
		if (bytes_per_sample == 2) {
			WORD *data = convert_data(fit);
			ret = save_mono_file(filename, data, width, height,
					bytes_per_sample);
		} else {
			uint8_t *data = convert_data8(fit);
			ret = save_mono_file(filename, data, width, height,
					bytes_per_sample);
			free(data);
		}
	}
	siril_log_message(_("Saving PNG: file %s, %ld layer(s), %ux%u pixels\n"),
			filename, fit->naxes[2], fit->rx, fit->ry);

	free(filename);
	return ret;
}

#endif	// HAVE_LIBPNG

/********************* RAW IMPORT *********************/
#ifdef HAVE_LIBRAW

static void get_FITS_date(time_t date, char *date_obs) {
	struct tm *t;
#ifdef HAVE_GMTIME_R
	struct tm t_;
#endif

#ifdef _WIN32
	t = gmtime (&date);
#else
#ifdef HAVE_GMTIME_R
	t = gmtime_r (&date, &t_);
#else
	t = gmtime(&date);
#endif /* HAVE_GMTIME_R */
#endif /* _WIN32 */

	/* If the gmtime() call has failed, "secs" is too big. */
	if (t) {
		g_snprintf(date_obs, FLEN_VALUE, "%04d-%02d-%02dT%02d:%02d:%02d",
				t->tm_year + 1900, t->tm_mon + 1, t->tm_mday, t->tm_hour, t->tm_min,
				t->tm_sec);
	}
}

#if (LIBRAW_MAJOR_VERSION == 0) && (LIBRAW_MINOR_VERSION < 18)
#define LIBRAW_FORMAT_1INCH 5
#endif

/* this is an estimation of the pixel size. Indeed, we cannot know
 * the real width resolution with libraw.
 * However, this approximation should be good enough.
 */
static float estimate_pixel_pitch(libraw_data_t *raw) {
	float s_width;

	switch (raw->lens.makernotes.CameraFormat) {
	case LIBRAW_FORMAT_APSC:
		if (!g_ascii_strncasecmp("Canon", raw->idata.make, 5))
			s_width = 22.3f;
		else
			s_width = 23.6f;
		break;
	case LIBRAW_FORMAT_FF:
		if (!g_ascii_strncasecmp("Sony", raw->idata.make, 4))
			s_width = 35.6f;
		else
			s_width = 36.0f;
		break;
	case LIBRAW_FORMAT_FT:
		s_width = 17.3f;
		break;
	case LIBRAW_FORMAT_APSH:
		s_width = 28.7f;
		break;
	case LIBRAW_FORMAT_1INCH:
		s_width = 13.2f;
		break;
	case LIBRAW_FORMAT_MF:
		s_width = 44.0f;
		break;
	default:
		s_width = 0.0f;
		break;
	}
//	printf("s_width=%f\n", s_width);
	float pitch = s_width / (float) raw->sizes.width * 1000.f;
	return roundf(pitch * 100.f) / 100.f;
}

static int siril_libraw_open_file(libraw_data_t* rawdata, const char *name) {
/* libraw_open_wfile is not defined for all windows compilers */
#if defined(_WIN32) && !defined(__MINGW32__) && defined(_MSC_VER) && (_MSC_VER > 1310)
	wchar_t *wname;

	wname = g_utf8_to_utf16(name, -1, NULL, NULL, NULL);
	if (wname == NULL) {
		return 1;
	}

	int ret = libraw_open_wfile(rawdata, wname);
	g_free(wname);
	return ret;
#elif defined(_WIN32)
	gchar *localefilename = g_win32_locale_filename_from_utf8(name);
	int ret = libraw_open_file(rawdata, localefilename);
	g_free(localefilename);
	return ret;
#else
	return(libraw_open_file(rawdata, name));
#endif
}

static int readraw(const char *name, fits *fit) {
	libraw_data_t *raw = libraw_init(0);

	int ret = siril_libraw_open_file(raw, name);
	if (ret) {
		siril_log_message(_("Error in libraw %s\n"), libraw_strerror(ret));
		libraw_recycle(raw);
		libraw_close(raw);
		return -1;
	}
	
	if (raw->other.shutter > 1.0)
		siril_log_message(_("Decoding %s %s file (ISO=%g, Exposure=%gs)\n"),
				raw->idata.make, raw->idata.model, raw->other.iso_speed, raw->other.shutter);
	else
		siril_log_message(_("Decoding %s %s file (ISO=%g, Exposure=1/%gs)\n"),
				raw->idata.make, raw->idata.model, raw->other.iso_speed, 1/raw->other.shutter);
		
	raw->params.output_bps = 16;						/* 16-bits files                           */
	raw->params.four_color_rgb = 0;						/* If == 1, interpolate RGB as four colors.*/
	raw->params.no_auto_bright = 1;						/* no auto_bright                          */
	raw->params.gamm[0] = 1.0 / com.raw_set.gamm[0];    /* Gamma curve set by the user             */
	raw->params.gamm[1] = com.raw_set.gamm[1];	                                         
	raw->params.bright = com.raw_set.bright;			/* Brightness                              */
	raw->params.user_flip = 0;							/* no flip                                 */
	raw->params.use_camera_wb = com.raw_set.use_camera_wb;
	raw->params.use_auto_wb = com.raw_set.use_auto_wb;
	if (com.raw_set.user_black == 1)
		raw->params.user_black = 0;						/* user black level equivalent to dcraw -k 0 */
	raw->params.output_color = 0;						/* output colorspace, 0=raw, 1=sRGB, 2=Adobe, 3=Wide, 4=ProPhoto, 5=XYZ*/
		
	if (!(com.raw_set.auto_mul)) { /* 4 multipliers (r,g,b,g) of the user's white balance.    */
		raw->params.user_mul[0] = (float) com.raw_set.mul[0];
		raw->params.user_mul[1] = raw->params.user_mul[3] = 1.0f;
		raw->params.user_mul[2] = (float) com.raw_set.mul[2];
		siril_log_message(_("Daylight multipliers: %f, %f, %f\n"),
				raw->params.user_mul[0], raw->params.user_mul[1],
				raw->params.user_mul[2]);
	} else {
		float mul[4]; /* 3 multipliers (r,g,b) from the camera white balance.  */
		mul[0] = raw->color.pre_mul[0] / raw->color.pre_mul[1];
		mul[1] = 1.0; /* raw->color.pre_mul[1]/raw->color.pre_mul[1]; */
		mul[2] = raw->color.pre_mul[2] / raw->color.pre_mul[1];
		mul[3] = raw->color.pre_mul[3] / raw->color.pre_mul[1];
		siril_log_message(_("Daylight multipliers: %f, %f, %f\n"), mul[0],
				mul[1], mul[2]);
	}

	if (raw->idata.filters == 9) {
		siril_log_color_message(_("XTRANS Sensor detected.\n"), "red");
	}

	switch (com.raw_set.user_qual) { /* Set interpolation                                        */
	case 0: /* bilinear interpolaton */
		raw->params.user_qual = 0;
		siril_log_message(_("Bilinear interpolation...\n"));
		break;
	case 2: /* VNG interpolaton */
		raw->params.user_qual = 1;
		siril_log_message(_("VNG interpolation...\n"));
		break;
	case 3: /* PPG interpolaton */
		raw->params.user_qual = 2;
		siril_log_message(_("PPG interpolation...\n"));
		break;
	default:
	case 1: /* AHD interpolaton */
		raw->params.user_qual = 3;
		siril_log_message(_("AHD interpolation...\n"));
		break;
	}
	
	const ushort width = raw->sizes.iwidth;
	const ushort height = raw->sizes.iheight;
	const float pitch = estimate_pixel_pitch(raw);
	const unsigned long npixels = width * height;

	WORD *data = malloc(npixels * sizeof(WORD) * 3);
	if (!data) {
		PRINT_ALLOC_ERR;
		return -1;
	}
	WORD *buf[3] = { data, data + npixels, data + npixels * 2 };
	ret = libraw_unpack(raw);
	if (ret) {
		printf("Error in libraw %s\n", libraw_strerror(ret));
		free(data);
		libraw_recycle(raw);
		libraw_close(raw);
		return -1;
	}

	ret = libraw_dcraw_process(raw);
	if (ret) {
		printf("Error in libraw %s\n", libraw_strerror(ret));
		free(data);
		libraw_recycle(raw);
		libraw_close(raw);
		return -1;
	}

	libraw_processed_image_t *image = libraw_dcraw_make_mem_image(raw, &ret);
	if (ret) {
		printf("Error in libraw %s\n", libraw_strerror(ret));
		free(data);
		libraw_dcraw_clear_mem(image);
		libraw_recycle(raw);
		libraw_close(raw);
		return -1;
	}

	int nbplanes = image->colors;
	if (nbplanes != 3) {
		free(data);
		libraw_dcraw_clear_mem(image);
		libraw_recycle(raw);
		libraw_close(raw);
		return -1;
	}
	// only for 16-bits because of endianness. Are there 8-bits RAW ???

	for (int i = 0; i < image->data_size; i += 6) {
		*buf[RLAYER]++ = (image->data[i + 0]) + (image->data[i + 1] << 8);
		*buf[GLAYER]++ = (image->data[i + 2]) + (image->data[i + 3] << 8);
		*buf[BLAYER]++ = (image->data[i + 4]) + (image->data[i + 5] << 8);
	}

	/*  Here we compute the correct size of the output image (imgdata.sizes.iwidth and imgdata.sizes.iheight) for the following cases:
    	- Files from Fuji cameras (with a 45-degree rotation)
    	- Files from cameras with non-square pixels
    	- Images shot by a rotated camera.
	 */
	libraw_adjust_sizes_info_only(raw);
	
	if (data != NULL) {
		clearfits(fit);
		fit->bitpix = fit->orig_bitpix = USHORT_IMG;
		fit->type = DATA_USHORT;
		fit->rx = (unsigned int) width;
		fit->ry = (unsigned int) height;
		fit->naxes[0] = (long) width;
		fit->naxes[1] = (long) height;
		fit->naxes[2] = nbplanes;
		if (nbplanes == 1)
			fit->naxis = 2;
		else
			fit->naxis = 3;
		fit->data = data;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data + npixels;
		fit->pdata[BLAYER] = fit->data + npixels * 2;
		fit->binning_x = fit->binning_y = 1;
		if (pitch > 0.f)
			fit->pixel_size_x = fit->pixel_size_y = pitch;
		if (raw->other.focal_len > 0.f)
			fit->focal_length = raw->other.focal_len;
		if (raw->other.iso_speed > 0.f)
			fit->iso_speed = raw->other.iso_speed;
		if (raw->other.shutter > 0.f)
			fit->exposure = raw->other.shutter;
		if (raw->other.aperture > 0.f)
			fit->aperture = raw->other.aperture;
		g_snprintf(fit->instrume, FLEN_VALUE, "%s %s", raw->idata.make,
				raw->idata.model);
		get_FITS_date(raw->other.timestamp, fit->date_obs);
		mirrorx(fit, FALSE);
	}

	libraw_dcraw_clear_mem(image);
	libraw_recycle(raw);
	libraw_close(raw);

	return nbplanes;
}

#define FC(filters, row, col) \
	(filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)

static const char filter[16][16] =
{ { 2,1,1,3,2,3,2,0,3,2,3,0,1,2,1,0 },
  { 0,3,0,2,0,1,3,1,0,1,1,2,0,3,3,2 },
  { 2,3,3,2,3,1,1,3,3,1,2,1,2,0,0,3 },
  { 0,1,0,1,0,2,0,2,2,0,3,0,1,3,2,1 },
  { 3,1,1,2,0,1,0,2,1,3,1,3,0,1,3,0 },
  { 2,0,0,3,3,2,3,1,2,0,2,0,3,2,2,1 },
  { 2,3,3,1,2,1,2,1,2,1,1,2,3,0,0,1 },
  { 1,0,0,2,3,0,0,3,0,3,0,3,2,1,2,3 },
  { 2,3,3,1,1,2,1,0,3,2,3,0,2,3,1,3 },
  { 1,0,2,0,3,0,3,2,0,1,1,2,0,1,0,2 },
  { 0,1,1,3,3,2,2,1,1,3,3,0,2,1,3,2 },
  { 2,3,2,0,0,1,3,0,2,0,1,2,3,0,1,0 },
  { 1,3,1,2,3,2,3,2,0,2,0,1,1,0,3,0 },
  { 0,2,0,3,1,0,0,1,1,3,3,2,3,2,2,1 },
  { 2,1,3,2,3,1,2,1,0,3,0,2,0,2,0,2 },
  { 0,3,1,0,0,2,0,3,2,1,3,1,1,3,1,3 } };

static int fcol(libraw_data_t *raw, int row, int col) {
	if (raw->idata.filters == 1)
		return filter[(row + raw->rawdata.sizes.top_margin) & 15][(col
				+ raw->rawdata.sizes.left_margin) & 15];
	if (raw->idata.filters == 9)
		return raw->idata.xtrans[(row + 6) % 6][(col + 6) % 6];
	return FC(raw->idata.filters, row, col);
}

static int readraw_in_cfa(const char *name, fits *fit) {
	libraw_data_t *raw = libraw_init(0);
	char pattern[FLEN_VALUE];

	int ret = siril_libraw_open_file(raw, name);
	if (ret) {
		printf("Error in libraw %s\n", libraw_strerror(ret));
		return -1;
	}

	ret = libraw_unpack(raw);
	if (ret) {
		printf("Error in libraw %s\n", libraw_strerror(ret));
		return -1;
	}

	/* This test checks if raw data exist. Sometimes it doesn't. This is
	 * the case for DNG built from lightroom for example */
	if ((void*) raw->rawdata.raw_image == 0x00
			&& (raw->rawdata.color3_image || raw->rawdata.color4_image)) {
		siril_log_message(_("Siril cannot open this file in CFA mode (no data available). "
				"Try to switch into RGB.\n"));
		return -1;
	}

	raw->params.user_flip = 0;				/* no flip                                 */
	raw->params.output_color = 0;			/* output colorspace, 0=raw, 1=sRGB, 2=Adobe, 3=Wide, 4=ProPhoto, 5=XYZ*/

	const ushort raw_width = raw->sizes.raw_width;
	const ushort raw_height = raw->sizes.raw_height;
	const ushort left_margin = raw->rawdata.sizes.left_margin;
	const ushort top_margin = raw->rawdata.sizes.top_margin;

	ushort width, height;

	if (raw->rawdata.ioparams.fuji_width) {
		const ushort right_margin = raw_width - raw->rawdata.ioparams.fuji_width
				- left_margin;
		width = raw_width - right_margin;
		height = raw_height;
	} else {
		width = raw->sizes.iwidth;
		height = raw->sizes.iheight;
	}

	float pitch = estimate_pixel_pitch(raw);
	const int npixels = width * height;
	
	if (raw->other.shutter > 0 && raw->other.shutter < 1)
		siril_log_message(_("Decoding %s %s file (ISO=%g, Exposure=1/%0.1f sec)\n"),
						raw->idata.make, raw->idata.model, raw->other.iso_speed, 1/raw->other.shutter);
	else
		siril_log_message(_("Decoding %s %s file (ISO=%g, Exposure=%0.1f sec)\n"),
						raw->idata.make, raw->idata.model, raw->other.iso_speed, raw->other.shutter);

	unsigned filters = raw->idata.filters;

	if (filters) {
		int fhigh = 2, fwide = 2;
		if ((filters ^ (filters >> 8)) & 0xff)
			fhigh = 4;
		if ((filters ^ (filters >> 16)) & 0xffff)
			fhigh = 8;
		if (filters == 1) /* Leaf Catchlight with 16x16 bayer matrix */
			fhigh = fwide = 16;
		if (filters == 9) /* Fuji X-Trans (6x6 matrix) */
			fhigh = fwide = 6;

		int j = 0;
		for (int i = 0; i < fhigh; i++) {
			for (int c = i /*&& (pattern[j++] = '/')*/ && 0; c < fwide; c++) {
				pattern[j++] = raw->idata.cdesc[fcol(raw, i, c)];
			}
		}
		pattern[j++] = '\0';
		siril_log_message(_("CFA pattern: %s\n"), pattern);
	}

	WORD *data = (WORD*) calloc(1, npixels * sizeof(WORD));
	if (!data) {
		PRINT_ALLOC_ERR;
		libraw_recycle(raw);
		libraw_close(raw);
		return -1;
	}

	WORD *buf = data;

	int offset = raw_width * top_margin + left_margin;

	int i = 0;
	for (int row = height - 1; row > -1; row--) {
		for (int col = 0; col < width; col++) {
			buf[i++] = raw->rawdata.raw_image[offset + col + (raw_width * row)];
		}
	}

	clearfits(fit);
	fit->bitpix = fit->orig_bitpix = USHORT_IMG;
	fit->type = DATA_USHORT;
	fit->rx = (unsigned int) (width);
	fit->ry = (unsigned int) (height);
	fit->naxes[0] = (long) (width);
	fit->naxes[1] = (long) (height);
	fit->naxes[2] = 1;
	fit->naxis = 2;
	fit->data = data;
	fit->pdata[RLAYER] = fit->data;
	fit->pdata[GLAYER] = fit->data;
	fit->pdata[BLAYER] = fit->data;
	fit->binning_x = fit->binning_y = 1;
	if (pitch > 0.f)
		fit->pixel_size_x = fit->pixel_size_y = pitch;
	if (raw->other.focal_len > 0.f)
		fit->focal_length = raw->other.focal_len;
	if (raw->other.iso_speed > 0.f)
		fit->iso_speed = raw->other.iso_speed;
	if (raw->other.shutter > 0.f)
		fit->exposure = raw->other.shutter;
	if (raw->other.aperture > 0.f)
		fit->aperture = raw->other.aperture;
	g_snprintf(fit->instrume, FLEN_VALUE, "%s %s", raw->idata.make,
			raw->idata.model);
	get_FITS_date(raw->other.timestamp, fit->date_obs);
	if (filters)
		g_snprintf(fit->bayer_pattern, FLEN_VALUE, "%s", pattern);


	libraw_recycle(raw);
	libraw_close(raw);
	return 1;
}

int open_raw_files(const char *name, fits *fit, int type) {
	int retvalue = 1;

	switch (type) {
	default:
	case 0:
		retvalue = readraw(name, fit);
		break;
	case 1:
		retvalue = readraw_in_cfa(name, fit);
		break;
	}
	if (retvalue >= 0) {
		gchar *basename = g_path_get_basename(name);
		siril_log_message(_("Reading RAW: file %s, %ld layer(s), %ux%u pixels\n"),
				basename, fit->naxes[2], fit->rx, fit->ry);
		g_free(basename);
	}
	return retvalue;
}
#endif

#ifdef HAVE_LIBHEIF
#define MAX_THUMBNAIL_SIZE thumbnail_size

struct HeifImage {
	uint32_t ID;
	char caption[100]; // image text (filled with resolution description)
	struct heif_image *thumbnail;
	int width, height;
};

static gboolean load_thumbnails(struct heif_context *heif, struct HeifImage *images) {
	int i;
	int numImages = heif_context_get_number_of_top_level_images(heif);

	// get list of all (top level) image IDs

	uint32_t *IDs = malloc(numImages * sizeof(uint32_t));
	heif_context_get_list_of_top_level_image_IDs(heif, IDs, numImages);

	// --- Load a thumbnail for each image.

	for (i = 0; i < numImages; i++) {

		images[i].ID = IDs[i];
		images[i].caption[0] = 0;
		images[i].thumbnail = NULL;

		// get image handle

		struct heif_image_handle *handle;
		struct heif_error err = heif_context_get_image_handle(heif, IDs[i],
				&handle);
		if (err.code) {
			g_printf("%s\n", err.message);
			continue;
		}

		// generate image caption

		int width = heif_image_handle_get_width(handle);
		int height = heif_image_handle_get_height(handle);

		if (heif_image_handle_is_primary_image(handle)) {
			sprintf(images[i].caption, "%dx%d (%s)", width, height,
					_("primary"));
		} else {
			sprintf(images[i].caption, "%dx%d", width, height);
		}

		// get handle to thumbnail image
		// if there is no thumbnail image, just the the image itself (will be scaled down later)

		struct heif_image_handle *thumbnail_handle;
		heif_item_id thumbnail_ID;

		int nThumbnails = heif_image_handle_get_list_of_thumbnail_IDs(handle,
				&thumbnail_ID, 1);

		if (nThumbnails > 0) {
			err = heif_image_handle_get_thumbnail(handle, thumbnail_ID,
					&thumbnail_handle);
			if (err.code) {
				g_printf("%s\n", err.message);
				continue;
			}
		} else {
			err = heif_context_get_image_handle(heif, IDs[i],
					&thumbnail_handle);
			if (err.code) {
				g_printf("%s\n", err.message);
				continue;
			}
		}

		// decode the thumbnail image

		struct heif_image *thumbnail_img;
		err = heif_decode_image(thumbnail_handle, &thumbnail_img,
				heif_colorspace_RGB, heif_chroma_interleaved_24bit,
				NULL);
		if (err.code) {
			g_printf("%s\n", err.message);
			continue;
		}

		// if thumbnail image size exceeds the maximum, scale it down

		int thumbnail_width = heif_image_handle_get_width(thumbnail_handle);
		int thumbnail_height = heif_image_handle_get_height(thumbnail_handle);

		if (thumbnail_width > MAX_THUMBNAIL_SIZE
				|| thumbnail_height > MAX_THUMBNAIL_SIZE) {

			// compute scaling factor to fit into a max sized box

			float factor_h = thumbnail_width / (float) MAX_THUMBNAIL_SIZE;
			float factor_v = thumbnail_height / (float) MAX_THUMBNAIL_SIZE;

			int new_width, new_height;

			if (factor_v > factor_h) {
				new_height = MAX_THUMBNAIL_SIZE;
				new_width = thumbnail_width / factor_v;
			} else {
				new_height = thumbnail_height / factor_h;
				new_width = MAX_THUMBNAIL_SIZE;
			}

			// scale the image

			struct heif_image *scaled_img = NULL;

			err = heif_image_scale_image(thumbnail_img,
					&scaled_img, new_width, new_height,
					NULL);
			if (err.code) {
				g_printf("%s\n", err.message);
				continue;
			}

			// release the old image and only keep the scaled down version

			heif_image_release(thumbnail_img);
			thumbnail_img = scaled_img;

			thumbnail_width = new_width;
			thumbnail_height = new_height;
		}

		heif_image_handle_release(thumbnail_handle);
		heif_image_handle_release(handle);

		// remember the HEIF thumbnail image (we need it for the GdkPixbuf)

		images[i].thumbnail = thumbnail_img;

		images[i].width = thumbnail_width;
		images[i].height = thumbnail_height;
	}

	return TRUE;
}

static gboolean heif_dialog(struct heif_context *heif, uint32_t *selected_image) {
	GtkWidget *dlg;
	GtkWidget *frame;
	gboolean run = FALSE;
	int i;

	int numImages = heif_context_get_number_of_top_level_images(heif);

	struct HeifImage *heif_images = malloc(numImages * sizeof(struct HeifImage));
	gboolean success = load_thumbnails(heif, heif_images);
	if (!success) {
		free(heif_images);
		return FALSE;
	}

	dlg = gtk_dialog_new_with_buttons(_("Load HEIF image content"),
			GTK_WINDOW(lookup_widget("control_window")), GTK_DIALOG_MODAL,
			_("_Cancel"), GTK_RESPONSE_CANCEL, _("_OK"), GTK_RESPONSE_OK, NULL);
	gtk_dialog_set_default_response(GTK_DIALOG(dlg), GTK_RESPONSE_OK);

	GtkContainer *content_area = GTK_CONTAINER(gtk_dialog_get_content_area(GTK_DIALOG(dlg)));
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 12);

	frame = gtk_frame_new(_("Select image"));
	gtk_container_add(content_area, GTK_WIDGET(frame));
	gtk_widget_show(frame);

// prepare list store with all thumbnails and caption

	GtkListStore *liststore;
	GtkTreeIter iter;

	liststore = gtk_list_store_new(2, G_TYPE_STRING, GDK_TYPE_PIXBUF);

	for (i = 0; i < numImages; i++) {
		gtk_list_store_append(liststore, &iter);
		gtk_list_store_set(liststore, &iter, 0, heif_images[i].caption, -1);

		int stride;
		const uint8_t *data = heif_image_get_plane_readonly(
				heif_images[i].thumbnail, heif_channel_interleaved, &stride);

		GdkPixbuf *pixbuf = gdk_pixbuf_new_from_data(data, GDK_COLORSPACE_RGB,
				FALSE, 8, heif_images[i].width, heif_images[i].height, stride,
				NULL, NULL);
		gtk_list_store_set(liststore, &iter, 1, pixbuf, -1);
	}

	GtkWidget *iconview = gtk_icon_view_new();
	gtk_icon_view_set_model((GtkIconView*) iconview, (GtkTreeModel*) liststore);
	gtk_icon_view_set_text_column((GtkIconView*) iconview, 0);
	gtk_icon_view_set_pixbuf_column((GtkIconView*) iconview, 1);
	gtk_icon_view_set_item_width((GtkIconView*) iconview, MAX_THUMBNAIL_SIZE);

	GtkWidget *scroll = gtk_scrolled_window_new(NULL, NULL);
	gtk_widget_set_size_request(scroll, -1, 400);
	g_object_set(scroll, "expand", TRUE, NULL);

	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll),	GTK_POLICY_NEVER, GTK_POLICY_ALWAYS);

	gtk_container_add(GTK_CONTAINER(frame), scroll);
	gtk_container_add(GTK_CONTAINER(scroll), iconview);

	gtk_widget_show(scroll);
	gtk_widget_show(iconview);

// pre-select the primary image

	int selected_idx = -1;
	for (i = 0; i < numImages; i++) {
		if (heif_images[i].ID == *selected_image) {
			selected_idx = i;
			break;
		}
	}

	if (selected_idx != -1) {
		GtkTreePath *path = gtk_tree_path_new_from_indices(selected_idx, -1);
		gtk_icon_view_select_path((GtkIconView*) iconview, path);
		gtk_tree_path_free(path);
	}

	gtk_widget_show(dlg);

	run = (gtk_dialog_run(GTK_DIALOG(dlg)) == GTK_RESPONSE_OK);

	if (run) {
		GList *selected_items = gtk_icon_view_get_selected_items(
				(GtkIconView*) iconview);

		if (selected_items) {
			GtkTreePath *path = (GtkTreePath*) (selected_items->data);
			gint *indices = gtk_tree_path_get_indices(path);

			*selected_image = heif_images[indices[0]].ID;

			g_list_free_full(selected_items,
					(GDestroyNotify) gtk_tree_path_free);
		}
	}

	gtk_widget_destroy(dlg);

// release thumbnail images

	for (i = 0; i < numImages; i++) {
		heif_image_release(heif_images[i].thumbnail);
	}

	free(heif_images);

	return run;
}

int readheif(const char* name, fits *fit, gboolean interactive){
	struct heif_error err;

	struct heif_context *ctx = heif_context_alloc();
	err = heif_context_read_from_file(ctx, name, NULL);
	if (err.code) {
		g_printf("%s\n", err.message);
		heif_context_free(ctx);
		return 1;
	}

	// analyze image content
	int num = heif_context_get_number_of_top_level_images(ctx);
	if (num == 0) {
		siril_log_message(_("Input file contains no readable images.\n"));
		heif_context_free(ctx);
		return 1;
	}

	  // get the primary image

	heif_item_id primary;

	err = heif_context_get_primary_image_ID(ctx, &primary);
	if (err.code) {
		g_printf("%s\n", err.message);
		heif_context_free(ctx);
		return 1;
	}

	// if primary image is no top level image or not present (invalid file), just take the first image

	if (!heif_context_is_top_level_image_ID(ctx, primary)) {
		int n = heif_context_get_list_of_top_level_image_IDs(ctx, &primary, 1);
		g_assert(n == 1);
	}

	heif_item_id selected_image = primary;

	if (num > 1) {
		if (!interactive) {
			siril_log_message(_("This is a sequence of %d images: "
					"loading the primary one.\n"), num);
		} else {
			if (!heif_dialog(ctx, &selected_image)) {
				heif_context_free(ctx);
				return 1;
			}
		}
	}

	// get a handle to the primary image
	struct heif_image_handle *handle;
	err = heif_context_get_image_handle(ctx, selected_image, &handle);
	if (err.code) {
		g_printf("%s\n", err.message);
		heif_context_free(ctx);
		return 1;
	}

	int has_alpha = heif_image_handle_has_alpha_channel(handle);

	struct heif_image *img = 0;
	err = heif_decode_image(handle, &img, heif_colorspace_RGB,
			has_alpha ? heif_chroma_interleaved_32bit :
			heif_chroma_interleaved_24bit, NULL);
	if (err.code) {
		g_printf("%s\n", err.message);
		heif_image_handle_release(handle);
		heif_context_free(ctx);
		return 1;
	}

	int stride;
	const uint8_t* udata = heif_image_get_plane_readonly(img, heif_channel_interleaved, &stride);
	const int width = heif_image_get_width(img, heif_channel_interleaved);
	const int height = heif_image_get_height(img, heif_channel_interleaved);

	long npixels = width * height;

	WORD *data = malloc(npixels * sizeof(WORD) * 3);
	WORD *buf[3] = { data, data + npixels, data + npixels * 2 };

	unsigned int nchannels = has_alpha ? 4 : 3;
	for (int row = 0; row < height; row += stride) {
		int nrow = (row + stride > height ? height - row : stride);
		for (int i = 0; i < width * nrow; i++) {
			*buf[RLAYER]++ = udata[i * nchannels + RLAYER];
			*buf[GLAYER]++ = udata[i * nchannels + GLAYER];
			*buf[BLAYER]++ = udata[i * nchannels + BLAYER];
		}
	}

	if (data != NULL) {
		clearfits(fit);
		fit->bitpix = fit->orig_bitpix = BYTE_IMG;
		fit->type = DATA_USHORT;
		fit->naxis = 3;
		fit->rx = width;
		fit->ry = height;
		fit->naxes[0] = fit->rx;
		fit->naxes[1] = fit->ry;
		fit->naxes[2] = 3;
		fit->data = data;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data + npixels;
		fit->pdata[BLAYER] = fit->data + npixels * 2;
		fit->binning_x = fit->binning_y = 1;
		mirrorx(fit, FALSE);
	}

	heif_image_handle_release(handle);
	heif_context_free(ctx);
	heif_image_release(img);
	gchar *basename = g_path_get_basename(name);
	siril_log_message(_("Reading HEIF: file %s, %ld layer(s), %ux%u pixels\n"),
			basename, fit->naxes[2], fit->rx, fit->ry);
	g_free(basename);

	return 0;
}
#endif
