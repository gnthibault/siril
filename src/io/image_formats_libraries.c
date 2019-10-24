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

#ifdef HAVE_CONFIG_H
#  include <config.h>
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

#include "core/siril.h"
#include "core/proto.h"
#include "algos/geometry.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "single_image.h"

/********************* TIFF IMPORT AND EXPORT *********************/

#ifdef HAVE_LIBTIFF

static int readtifstrip(TIFF* tif, uint32 width, uint32 height, uint16 nsamples, WORD **data) {
	unsigned int npixels;
	int i, j, scanline, retval = nsamples;
	WORD *buf;
	uint32 rowsperstrip;
	uint16 config;
	unsigned long nrow, row;

	TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
	TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);

	npixels = width * height;
	*data = malloc(npixels * sizeof(WORD) * nsamples);
	WORD *gbuf[3] = {*data, *data, *data};
	if (nsamples == 4) {
		siril_log_message(_("Alpha channel is ignored.\n"));
	}
	if ((nsamples == 3) || (nsamples == 4)) {
		gbuf[1] = *data + npixels;
		gbuf[2] = *data + npixels * 2;
	}

	scanline = TIFFScanlineSize(tif);
	buf = _TIFFmalloc(TIFFStripSize(tif));
	if (!buf) {
		PRINT_ALLOC_ERR;
		return -1;
	}
	for (row = 0; row < height; row += rowsperstrip){
		nrow = (row + rowsperstrip > height ? height - row : rowsperstrip);
		switch(config){
			case PLANARCONFIG_CONTIG:
				if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, row, 0), buf, nrow*scanline) < 0){
					siril_log_message(_("An unexpected error was encountered while trying to read the file.\n"));
					retval = -1;
					break;
				}
				for (i = 0; i < width*nrow; i++) {
					*gbuf[RLAYER]++ = buf[i*nsamples+0];
					if ((nsamples == 3) || (nsamples == 4)) {
						*gbuf[GLAYER]++ = buf[i*nsamples+1];
						*gbuf[BLAYER]++ = buf[i*nsamples+2];
					}
				}
				break;
			case PLANARCONFIG_SEPARATE:
				if (nsamples >= 3)		//don't need to read the alpha
					nsamples = 3;
				for (j=0; j<nsamples; j++){	//loop on the layer
					if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, row, j), buf, nrow*scanline) < 0){
						siril_log_message(_("An unexpected error was encountered while trying to read the file.\n"));
						retval = -1;
						break;
					}
					for (i = 0; i < width*nrow; i++)
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
	uint32 npixels;
	int retval = nsamples;

	npixels = width * height;
	*data = malloc(npixels * sizeof(WORD) * nsamples);
	WORD *gbuf[3] = {*data, *data, *data};
	if (nsamples == 4) {
		siril_log_message(_("Alpha channel is ignored.\n"));
	}
	if ((nsamples == 3) || (nsamples == 4)) {
		gbuf[1] = *data + npixels;
		gbuf[2] = *data + npixels * 2;
	}

	/* get the data */
	uint32* raster = (uint32*)_TIFFmalloc(npixels*sizeof(uint32));
	if (raster != NULL) {
		if (TIFFReadRGBAImage(tif, width, height, raster, 0)) {
			int i, j;
			for (j = 0; j < height; j++) {
				int istart = j * width;
				for (i = 0; i < width; i++) {
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
	uint32 height, width, npixels;
	uint16 nbits, nsamples, color;
	WORD *data;
	
	TIFF* tif = Siril_TIFFOpen(name, "r");

	if (!tif) {
		siril_log_message(_("Could not open the TIFF file %s\n"), name);
		return -1;
	}

	data = NULL;
	
	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
	TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);	
	TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &nbits);
	TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &color);
	TIFFGetField(tif, TIFFTAG_MINSAMPLEVALUE, &(fit->lo));
	TIFFGetField(tif, TIFFTAG_MAXSAMPLEVALUE, &(fit->hi));

	npixels = width * height;

	switch(nbits){
		case 8:
			/* High level functions in readtif8bits: should read every 8-bit TIFF file */
			retval = readtif8bits(tif, width, height, nsamples, &data);
			break;

		case 16:
			/* Lower level functions in readtifstrip: 
			 * maybe some particular TIFF files could not be read.
			 * No such files are known yet */
			retval = readtifstrip(tif, width, height, nsamples, &data);
			break;

		default :
			siril_log_message(_("Siril only works with 8/16-bit TIFF format.\n"));
			retval = -1;
	}
	TIFFClose(tif);
	if (retval <= 0) {
		if (data)
			free(data);
	} else if (data != NULL) {
		clearfits(fit);
		fit->rx = width;
		fit->ry = height;
		fit->naxes[0] = width;
		fit->naxes[1] = height;
		fit->data = data;
		fit->binning_x=fit->binning_y=1;
		if (nsamples == 1 || nsamples == 2) {
			fit->naxes[2] = 1;
			fit->naxis = 2;
			fit->pdata[RLAYER]=fit->data;
			fit->pdata[GLAYER]=fit->data;
			fit->pdata[BLAYER]=fit->data;
		} else {
			fit->naxes[2] = 3;
			fit->naxis = 3;
			fit->pdata[RLAYER]=fit->data;
			fit->pdata[GLAYER]=fit->data + npixels;
			fit->pdata[BLAYER]=fit->data + npixels * 2;
		}
		fit->bitpix = (nbits == 8) ? BYTE_IMG : USHORT_IMG;
		fit->orig_bitpix = fit->bitpix;
		retval = nsamples;
	}
	if (nbits==16) mirrorx(fit, FALSE);
	gchar *basename = g_path_get_basename(name);
	siril_log_message(_("Reading TIFF: %d-bit file %s, %ld layer(s), %ux%u pixels\n"),
						nbits, basename, fit->naxes[2], fit->rx, fit->ry);
	g_free(basename);

	return retval;
}

/*** This function save the current image into a uncompressed 8- or 16-bit file *************/

int savetif(const char *name, fits *fit, uint16 bitspersample){
	int retval = 0;
	char *filename;
	unsigned char *buf8;
	WORD *buf16;
	uint32 width, height, row, col, n;
	uint16 nsamples;
	double norm;

	mirrorx(fit, FALSE);

	filename = strdup(name);
	if (!ends_with(filename, ".tif") && (!ends_with(filename, ".tiff"))) {
		filename = str_append(&filename, ".tif");
	}

	TIFF* tif = Siril_TIFFOpen(filename, "w");

	if (tif == NULL) {
		siril_log_message(_("Siril cannot create TIFF file.\n"));
		free(filename);
		return 1;
	}
	nsamples = (uint16) fit->naxes[2];
	width = (uint32) fit->rx;
	height = (uint32) fit->ry;
	
	/*******************************************************************
	 * If the user saves a tif from the graphical menu, he can set
	 * the Description and the Copyright of the Image 
	 ******************************************************************/
	gchar *img_desc = NULL, *img_copy = NULL;
	GtkTextView *description = GTK_TEXT_VIEW(lookup_widget("Description_txt"));
	GtkTextView *copyright = GTK_TEXT_VIEW(lookup_widget("Copyright_txt"));
	GtkTextBuffer *desbuf = gtk_text_view_get_buffer(description);
	GtkTextBuffer *copybuf = gtk_text_view_get_buffer(copyright);
	GtkTextIter itDebut;
	GtkTextIter itFin;

	gtk_text_buffer_get_start_iter(desbuf, &itDebut);
	gtk_text_buffer_get_end_iter(desbuf, &itFin);
	if (desbuf)
		img_desc = gtk_text_buffer_get_text(desbuf, &itDebut, &itFin, TRUE);
	gtk_text_buffer_get_bounds(desbuf, &itDebut, &itFin);
	gtk_text_buffer_delete(desbuf, &itDebut, &itFin);
	gtk_text_buffer_get_start_iter(copybuf, &itDebut);
	gtk_text_buffer_get_end_iter(copybuf, &itFin);
	if (copybuf)
		img_copy = gtk_text_buffer_get_text(copybuf, &itDebut, &itFin, TRUE);
	gtk_text_buffer_get_bounds(copybuf, &itDebut, &itFin);
	gtk_text_buffer_delete(copybuf, &itDebut, &itFin);

	/*******************************************************************/

	/* TIFF TAG FIELD */
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bitspersample);
	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
	TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tif, -1));
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, nsamples);
	TIFFSetField(tif, TIFFTAG_COMPRESSION, get_compression_mode());
	TIFFSetField(tif, TIFFTAG_IMAGEDESCRIPTION, img_desc);
	TIFFSetField(tif, TIFFTAG_COPYRIGHT, img_copy);
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

	WORD *gbuf[3] =
			{ fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };

	switch (bitspersample) {
	case 8:
		buf8 = _TIFFmalloc(width * sizeof(unsigned char) * nsamples);
		norm = fit->orig_bitpix != BYTE_IMG ? UCHAR_MAX_DOUBLE / USHRT_MAX_DOUBLE : 1.0;
		for (row = 0; row < height; row++) {
			for (col = 0; col < width; col++) {
				for (n = 0; n < nsamples; n++) {
					/* UCHAR_MAX / USHRT_MAX is constant, it's 1/255
					 * This operation should be speed-up by doing a shift */
					buf8[col * nsamples + n] = (*gbuf[n]++) * norm;
				}
			}
			TIFFWriteScanline(tif, buf8, row, 0);
		}
		_TIFFfree(buf8);
		break;
	case 16:
		buf16 = _TIFFmalloc(width * sizeof(WORD) * nsamples);
		norm = (fit->orig_bitpix != BYTE_IMG) ? 1.0 : USHRT_MAX_DOUBLE / UCHAR_MAX_DOUBLE;

		for (row = 0; row < height; row++) {
			for (col = 0; col < width; col++) {
				for (n = 0; n < nsamples; n++)
					buf16[col * nsamples + n] = (*gbuf[n]++) * norm;
			}
			TIFFWriteScanline(tif, buf16, row, 0);
		}
		_TIFFfree(buf16);
		break;
	default:		// Should not happen
		retval = 1;
	}
	TIFFClose(tif);
	mirrorx(fit, FALSE);
	siril_log_message(
			_("Saving TIFF: %d-bit file %s, %ld layer(s), %ux%u pixels\n"),
			bitspersample, filename, fit->naxes[2], fit->rx, fit->ry);
	free(filename);
	return retval;
}
#endif	// HAVE_LIBTIFF


/********************* JPEG IMPORT AND EXPORT *********************/

#ifdef HAVE_LIBJPEG
int readjpg(const char* name, fits *fit){
	int i, npixels;
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	WORD *data;
	FILE *f;
	JSAMPARRAY pJpegBuffer;
	int row_stride;


	if ((f = g_fopen(name, "rb")) == NULL){
		siril_log_message(_("Sorry but Siril cannot open the file: %s.\n"), name);
		return -1;
	}
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, f);
	jpeg_read_header(&cinfo, TRUE);
	jpeg_start_decompress(&cinfo);
	npixels=cinfo.output_width*cinfo.output_height;
	data = malloc(npixels * sizeof(WORD) * 3);
	WORD *buf[3] = {data, data + npixels, data + npixels * 2};
	row_stride = cinfo.output_width * cinfo.output_components ;
	pJpegBuffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
	while (cinfo.output_scanline < cinfo.output_height) {
		jpeg_read_scanlines(&cinfo, pJpegBuffer, 1);
		for (i=0;i<cinfo.output_width;i++) {
			*buf[RLAYER]++ = pJpegBuffer[0][cinfo.output_components*i+0];
			*buf[GLAYER]++ = pJpegBuffer[0][cinfo.output_components*i+1];
			*buf[BLAYER]++ = pJpegBuffer[0][cinfo.output_components*i+2];
		}
	}
	fclose(f);
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	if (data != NULL){
		clearfits(fit);	
		fit->bitpix = fit->orig_bitpix = BYTE_IMG;
		if (cinfo.output_components==1)
			fit->naxis = 2;
		else
			fit->naxis = 3;
		fit->rx = cinfo.output_width;
		fit->ry = cinfo.output_height;
		fit->naxes[0] = cinfo.output_width;
		fit->naxes[1] = cinfo.output_height;
		fit->naxes[2] = cinfo.output_components;
		fit->data = data;
		fit->pdata[RLAYER]=fit->data;
		fit->pdata[GLAYER]=fit->data + npixels;
		fit->pdata[BLAYER]=fit->data + npixels * 2;
		fit->binning_x = fit->binning_y = 1;
	}
	mirrorx(fit, FALSE);
	gchar *basename = g_path_get_basename(name);
	siril_log_message(_("Reading JPG: file %s, %ld layer(s), %ux%u pixels\n"),
						basename, fit->naxes[2], fit->rx, fit->ry);
	g_free(basename);

	return cinfo.output_components;
}

int savejpg(const char *name, fits *fit, int quality){
	FILE *f;
	int i, j;
	WORD red, blue, green;
	struct jpeg_compress_struct cinfo;    // Basic info for JPEG properties.
	struct jpeg_error_mgr jerr;           // In case of error.
	JSAMPROW row_pointer[1];              // Pointer to JSAMPLE row[s].
	int row_stride;                       // Physical row width in image buffer.
	double norm;

	//## ALLOCATE AND INITIALIZE JPEG COMPRESSION OBJECT
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);

	char *filename = strdup(name);
	if (!ends_with(filename, ".jpg") && (!ends_with(filename, ".jpeg"))) {
		filename = str_append(&filename, ".jpg");
	}

	//## OPEN FILE FOR DATA DESTINATION:
	if ((f = g_fopen(filename, "wb")) == NULL) {
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

	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, quality, TRUE);

	//## CREATE IMAGE BUFFER TO WRITE FROM AND MODIFY THE IMAGE TO LOOK LIKE CHECKERBOARD:
	unsigned char *image_buffer = (unsigned char*) malloc(
			cinfo.image_width * cinfo.image_height * cinfo.num_components);

	norm = (fit->orig_bitpix != BYTE_IMG ?
			UCHAR_MAX_DOUBLE / USHRT_MAX_DOUBLE : 1.0);

	for (i = (cinfo.image_height - 1); i >= 0; i--) {
		for (j = 0; j < cinfo.image_width; j++) {
			int pixelIdx = ((i * cinfo.image_width) + j) * cinfo.input_components;
			red = *gbuf[RLAYER]++;
			image_buffer[pixelIdx + 0] = round_to_BYTE(red * norm);         // r |-- Set r,g,b components to
			if (cinfo.input_components == 3) {
				green = *gbuf[GLAYER]++;
				blue = *gbuf[BLAYER]++;
				image_buffer[pixelIdx + 1] = round_to_BYTE(green * norm); // g |   make this pixel
				image_buffer[pixelIdx + 2] = round_to_BYTE(blue * norm);  // b |
			}
		}
	}
	//## START COMPRESSION:
	jpeg_start_compress(&cinfo, TRUE);
	row_stride = cinfo.image_width * cinfo.input_components;        // JSAMPLEs per row in image_buffer

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
	int width, height, x, y, npixels, nbplanes;
	png_byte color_type;
	png_byte bit_depth;
	png_bytep *row_pointers;
	WORD *data;
	FILE *f;

	png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL,
			NULL);
	if (!png) {
		siril_log_message(_("Sorry but Siril cannot open the file: %s.\n"), name);
		return -1;
	}

	png_infop info = png_create_info_struct(png);
	if(!info) return -1;

	if(setjmp(png_jmpbuf(png))) return -1;

	f = g_fopen(name, "rb");
	png_init_io(png, f);

	png_read_info(png, info);

	width = png_get_image_width(png, info);
	height = png_get_image_height(png, info);
	npixels = width * height;
	color_type = png_get_color_type(png, info);
	bit_depth = png_get_bit_depth(png, info);

	data = malloc(npixels * sizeof(WORD) * 3);
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

	row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
	for (y = 0; y < height; y++) {
		row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png, info));
	}

	png_read_image(png, row_pointers);

	fclose(f);

	if (bit_depth == 16) {			//in 16-bit: it is stored as RRGGBB
		for (y = 0; y < height; y++) {
			png_byte* row = row_pointers[y];
			for (x = 0; x < width; x++) {
				png_byte* ptr = &(row[x * 8]);
				*buf[RLAYER]++ = ptr[0] * (UCHAR_MAX + 1) + ptr[1];
				*buf[GLAYER]++ = ptr[2] * (UCHAR_MAX + 1) + ptr[3];
				*buf[BLAYER]++ = ptr[4] * (UCHAR_MAX + 1) + ptr[5];
			}
		}
	} else {
		for (y = 0; y < height; y++) {
			png_byte* row = row_pointers[y];
			for (x = 0; x < width; x++) {
				png_byte* ptr = &(row[x * 4]);
				*buf[RLAYER]++ = ptr[0];
				*buf[GLAYER]++ = ptr[1];
				*buf[BLAYER]++ = ptr[2];
			}
		}
	}
	// We define the number of channel we have
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
	for (y = 0; y < height; y++)
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
		fit->orig_bitpix = fit->bitpix;
		fit->data = data;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data + npixels;
		fit->pdata[BLAYER] = fit->data + npixels * 2;
		fit->binning_x = fit->binning_y = 1;

	}
	mirrorx(fit, FALSE);
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
	int row_stride;
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
	row_stride = image.width * -3;

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
	int row_stride;
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
	row_stride = image.width * -1;

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
	int i, j;

	WORD *buffer = malloc(ndata * ch * sizeof(WORD));
	for (i = 0, j = 0; i < ndata * ch; i += ch, j++) {
		buffer[i + 0] = image->pdata[RLAYER][j];
		if (ch > 1) {
			buffer[i + 1] = image->pdata[GLAYER][j];
			buffer[i + 2] = image->pdata[BLAYER][j];
		}
	}
	return buffer;
}

static uint8_t *convert_data8(fits *image) {
	int ndata = image->rx * image->ry;
	int ch = image->naxes[2];
	int i, j;

	uint8_t *buffer = malloc(ndata * ch * sizeof(uint8_t));
	for (i = 0, j = 0; i < ndata * ch; i += ch, j++) {
		buffer[i + 0] = (uint8_t)image->pdata[RLAYER][j];
		if (ch > 1) {
			buffer[i + 1] = (uint8_t) image->pdata[GLAYER][j];
			buffer[i + 2] = (uint8_t) image->pdata[BLAYER][j];
		}
	}
	return buffer;
}

int savepng(const char *name, fits *fit, uint32_t bytes_per_sample,
		gboolean is_colour) {
	int ret = -1;
	uint32_t width = fit->rx;
	uint32_t height = fit->ry;

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
			ret = save_mono_file(filename, fit->data, width, height,
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
			s_width = 22.3;
		else
			s_width = 23.6;
		break;
	case LIBRAW_FORMAT_FF:
		if (!g_ascii_strncasecmp("Sony", raw->idata.make, 4))
			s_width = 35.6;
		else
			s_width = 36.0;
		break;
	case LIBRAW_FORMAT_FT:
		s_width = 17.3;
		break;
	case LIBRAW_FORMAT_APSH:
		s_width = 28.7;
		break;
	case LIBRAW_FORMAT_1INCH:
		s_width = 13.2;
		break;
	case LIBRAW_FORMAT_MF:
		s_width = 44.0;
		break;
	default:
		s_width = 0.0;
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
	ushort width, height;
	int npixels, i;
	double pitch;
	WORD *data = NULL;
	libraw_data_t *raw = libraw_init(0);
	libraw_processed_image_t *image = NULL;
	int ret;

	ret = siril_libraw_open_file(raw, name);
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
	
	width = raw->sizes.iwidth;
	height = raw->sizes.iheight;
	pitch = estimate_pixel_pitch(raw);

	npixels = width * height;

	data = malloc(npixels * sizeof(WORD) * 3);
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
	image = libraw_dcraw_make_mem_image(raw, &ret);
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

	for (i = 0; i < image->data_size; i += 6) {
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
	width = raw->sizes.iwidth;
	height = raw->sizes.iheight;
	
	if (data != NULL) {
		clearfits(fit);
		fit->bitpix = fit->orig_bitpix = USHORT_IMG;
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
	unsigned int i, j, c, col, row;
	float pitch;
	char pattern[FLEN_VALUE];
	ushort raw_width, raw_height, left_margin, top_margin;
	ushort width, height;
	int npixels;
	WORD *data = NULL;
	int ret;

	ret = siril_libraw_open_file(raw, name);

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

	raw_width = raw->sizes.raw_width;
	raw_height = raw->sizes.raw_height;

	left_margin = raw->rawdata.sizes.left_margin;
	top_margin = raw->rawdata.sizes.top_margin;

	if (raw->rawdata.ioparams.fuji_width) {
		ushort right_margin = raw_width - raw->rawdata.ioparams.fuji_width
				- left_margin;
		width = raw_width - right_margin;
		height = raw_height;
	} else {
		width = raw->sizes.iwidth;
		height = raw->sizes.iheight;
	}

	pitch = estimate_pixel_pitch(raw);
	npixels = width * height;
	
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

		j = 0;
		for (i = 0; i < fhigh; i++) {
			for (c = i /*&& (pattern[j++] = '/')*/ && 0; c < fwide; c++) {
				pattern[j++] = raw->idata.cdesc[fcol(raw, i, c)];
			}
		}
		pattern[j++] = '\0';
		siril_log_message(_("CFA pattern: %s\n"), pattern);
	}

	data = (WORD*) calloc(1, npixels * sizeof(WORD));
	if (!data) {
		PRINT_ALLOC_ERR;
		libraw_recycle(raw);
		libraw_close(raw);
		return -1;
	}

	WORD *buf = data;

	int offset = raw_width * top_margin + left_margin;

	i = 0;
	for (row = 0; row < height; row++) {
		for (col = 0; col < width; col++) {
			buf[i++] = raw->rawdata.raw_image[offset + col + (raw_width * row)];
		}
	}

	clearfits(fit);
	fit->bitpix = fit->orig_bitpix = USHORT_IMG;
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
	fit->maximum_pixel_value = (WORD) raw->color.maximum;
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
		mirrorx(fit, FALSE);
		gchar *basename = g_path_get_basename(name);
		siril_log_message(_("Reading RAW: file %s, %ld layer(s), %ux%u pixels\n"),
				basename, fit->naxes[2], fit->rx, fit->ry);
		g_free(basename);
	}
	return retvalue;
}
#endif
