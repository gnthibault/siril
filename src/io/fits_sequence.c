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
 *
 *
 * FITS sequences are not a sequence of FITS files but a FITS file containing a
 * sequence. It simply has as many elements in the third dimension as the
 * number of images in the sequence multiplied by the number of channels per
 * image. Given its use of the third dimension, it's sometimes called FITS cube.
 */

#include "core/siril.h"

#include "io/image_format_fits.h"
#include "gui/progress_and_log.h"
#include "core/siril_log.h"

#include "fits_sequence.h"

static int fitseq_write_image_for_writer(struct seqwriter_data *writer, fits *image, int index);

static int _find_hdus(fitsfile *fptr, int **hdus, int *nb_im) {
	int status = 0;
	int nb_hdu, ref_naxis = -1, ref_bitpix = 0, nb_images = 0;
	long ref_naxes[3] = { 0l };

	fits_get_num_hdus(fptr, &nb_hdu, &status);
	if (status || !nb_im)
		return 1;

	if (hdus) {
		*hdus = malloc(nb_hdu * sizeof(int));
		if (!*hdus) {
			PRINT_ALLOC_ERR;
			return 1;
		}
	}

	for (int i = 0; i < nb_hdu; i++) {
		status = 0;
		int type;
		if (fits_movabs_hdu(fptr, i + 1, &type, &status)) {
			report_fits_error(status);
			return -1;
		}

		if (type != IMAGE_HDU) continue;

		long naxes[3];
		int naxis;
		int bitpix;
		fits_get_img_param(fptr, 3, &bitpix, &naxis, naxes, &status);
		if (status) {
			report_fits_error(status);
			break;
		}

		if (naxis > 0) {
			if (ref_naxis == -1) {
				ref_naxis = naxis;
				ref_bitpix = bitpix;
				memcpy(ref_naxes, naxes, sizeof naxes);
				siril_debug_print("found reference HDU %ldx%ldx%d (%d)\n", naxes[0], naxes[1], naxis, bitpix);
			} else {
				if (naxis != ref_naxis || naxes[0] != ref_naxes[0] || naxes[1] != ref_naxes[1] || bitpix != ref_bitpix) {
					fprintf(stderr, "another image was found in the FITS file but does not has the same parameters as the first one\n");
					break;
				}
			}
			if (hdus)
				(*hdus)[nb_images] = i + 1;
			nb_images++;
		}
	}

	if (status) {
		if (hdus) {
			free(*hdus);
			*hdus = NULL;
		}
	}
	else {
		*nb_im = nb_images;
		siril_debug_print("found %d images with same params in the FITS sequence\n", nb_images);
		// we could realloc *hdus, but it's not much useful
	}
	return status;
}


// test if a file is a multi-extension FITS, a.k.a FITS cube or FITS sequence
int fitseq_is_fitseq(const char *filename, int *frames) {
	fitsfile *fptr;
	int status = 0;
	if (siril_fits_open_diskfile(&fptr, filename, READONLY, &status))
		return 0;

	int nb_images;
	status = _find_hdus(fptr, NULL, &nb_images);
	if (frames) *frames = nb_images;

	int status2 = 0;
	fits_close_file(fptr, &status2);
	return !status && nb_images > 1;
}

void fitseq_init_struct(fitseq *fitseq) {
	fitseq->filename = NULL;
	fitseq->bitpix = 0;
	fitseq->orig_bitpix = 0;
	fitseq->naxes[0] = 0;
	fitseq->frame_count = 0;
	fitseq->hdu_index = NULL;
	fitseq->fptr = NULL;
	fitseq->is_mt_capable = FALSE;
#ifdef _OPENMP
	fitseq->thread_fptr = NULL;
#endif
	fitseq->writer = NULL;
}

int fitseq_open(const char *filename, fitseq *fitseq) {
	if (fitseq->fptr) {
		fprintf(stderr, "FITS sequence: file already opened, or badly closed\n");
		return -1;
	}

	int status = 0;
	siril_fits_open_diskfile(&(fitseq->fptr), filename, READONLY, &status);
	if (status) {
		report_fits_error(status);
		siril_log_color_message(_("Cannot open FITS file %s\n"), "red", filename);
		return -1;
	}

	if (_find_hdus(fitseq->fptr, &fitseq->hdu_index, &fitseq->frame_count) || fitseq->frame_count <= 1) {
		siril_log_color_message(_("Cannot open FITS file %s: doesn't seem to be a FITS sequence\n"), "red", filename);
		return -1;
	}

	if (fits_movabs_hdu(fitseq->fptr, fitseq->hdu_index[0], NULL, &status)) {
		report_fits_error(status);
		return -1;
	}

	// we store the first image's dimensions in the struct
	int naxis;
	status = 0;
	fits_get_img_param(fitseq->fptr, 3, &(fitseq->bitpix), &naxis, fitseq->naxes, &status);
	if (status || naxis <= 1 || fitseq->naxes[0] == 0 || fitseq->naxes[1] == 0) {
		status = 0;
		fits_close_file(fitseq->fptr, &status);
		return -1;
	}
	if (naxis == 2)
		fitseq->naxes[2] = 1;

	manage_bitpix(fitseq->fptr, &(fitseq->bitpix), &(fitseq->orig_bitpix));

	if (fitseq->bitpix == LONGLONG_IMG) {
		siril_log_message(
				_("FITS images with 64 bits signed integer per pixel.channel are not supported.\n"));
		status = 0;
		fits_close_file(fitseq->fptr, &status);
		return -1;
	}

	fitseq->filename = strdup(filename);
	fitseq->is_mt_capable = FALSE;
	siril_debug_print("fitseq_open: sequence %s has %d frames, bitpix = %d, naxis = %d, naxes = { %ld, %ld, %ld }\n",
			filename, fitseq->frame_count, fitseq->bitpix, naxis,
			fitseq->naxes[0], fitseq->naxes[1], fitseq->naxes[2]);

#ifdef _OPENMP
	if (fits_is_reentrant()) {
		fitseq->is_mt_capable = TRUE;
		fprintf(stdout, "cfitsio was compiled with multi-thread support,"
				" parallel read of images will be possible\n");
	} else {
		fprintf(stdout, "cfitsio was compiled without multi-thread support,"
				" parallel read of images will be impossible\n");
		siril_log_message(_("Your version of cfitsio does not support multi-threading\n"));
	}
#endif
	return 0;
}

/* dest must be filled with zeros */
static int fitseq_read_frame_internal(fitseq *fitseq, int index, fits *dest, gboolean force_float, fitsfile *fptr) {
	if (!fptr) {
		return -1;
	}

	memcpy(dest->naxes, fitseq->naxes, sizeof fitseq->naxes);
	dest->naxis = fitseq->naxes[2] == 3 ? 3 : 2;
	dest->bitpix = fitseq->bitpix;
	dest->orig_bitpix = fitseq->orig_bitpix;
	dest->rx = dest->naxes[0];
	dest->ry = dest->naxes[1];
	dest->fptr = fptr;

	siril_debug_print("reading HDU %d (of %s)\n", fitseq->hdu_index[index], fitseq->filename);
	int status = 0;
	if (fits_movabs_hdu(fptr, fitseq->hdu_index[index], NULL, &status)) {
		report_fits_error(status);
		return -1;
	}

	read_fits_header(dest);	// stores useful header data in fit
	dest->header = copy_header(dest); // for display

	if (read_fits_with_convert(dest, fitseq->filename, force_float)) {
		return -1;
	}
	return 0;
}

int fitseq_read_frame(fitseq *fitseq, int index, fits *dest, gboolean force_float, int thread) {
	fitsfile *fptr = fitseq->fptr;
#ifdef _OPENMP
	if (thread >= 0 && fitseq->thread_fptr) {
		fptr = fitseq->thread_fptr[thread];
		siril_debug_print("fitseq: thread %d reading FITS image\n", thread);
	}
#endif
	return fitseq_read_frame_internal(fitseq, index, dest, force_float, fptr);
}

// we read a partial image and return it as fits
int fitseq_read_partial_fits(fitseq *fitseq, int layer, int index, fits *dest, const rectangle *area, gboolean do_photometry, int thread) {
	dest->type = get_data_type(fitseq->bitpix);
	if (dest->type == DATA_UNSUPPORTED) {
		siril_log_message(_("Unknown FITS data format in internal conversion\n"));
		return -1;
	}
	if (new_fit_image(&dest, area->w, area->h, 1, dest->type))
		return -1;
	fitsfile *fptr = fitseq->fptr;
#ifdef _OPENMP
	if (thread >= 0 && fitseq->thread_fptr)
		fptr = fitseq->thread_fptr[thread];
#endif
	dest->fptr = fptr;
	dest->bitpix = fitseq->bitpix;
	dest->orig_bitpix = fitseq->orig_bitpix;

	int status = 0;
	if (fits_movabs_hdu(fptr, fitseq->hdu_index[index], NULL, &status)) {
		report_fits_error(status);
		return -1;
	}

	if (do_photometry)
		fit_get_photometry_data(dest);

	status = internal_read_partial_fits(fptr, fitseq->naxes[1], fitseq->bitpix,
			dest->type == DATA_USHORT ? (void *)dest->data : (void *)dest->fdata,
			layer, area);
	return status;
}

// we read a partial image and return it as buffer
int fitseq_read_partial(fitseq *fitseq, int layer, int index, void *buffer, const rectangle *area, int thread) {
	if (area->x < 0 || area->y < 0 || area->x >= fitseq->naxes[0] || area->y >= fitseq->naxes[1]
			|| area->w <= 0 || area->h <= 0 || area->x + area->w > fitseq->naxes[0]
			|| area->y + area->h > fitseq->naxes[1]) {
		fprintf(stderr, "partial read from FITS file has been requested outside image bounds or with invalid size\n");
		return 1;
	}

	fitsfile *fptr = fitseq->fptr;
#ifdef _OPENMP
	if (thread >= 0 && fitseq->thread_fptr)
		fptr = fitseq->thread_fptr[thread];
#endif

	int status = 0;
	if (fits_movabs_hdu(fptr, fitseq->hdu_index[index], NULL, &status)) {
		report_fits_error(status);
		return -1;
	}

	if (internal_read_partial_fits(fptr, fitseq->naxes[1], fitseq->bitpix, buffer, layer, area))
		return 1;
	flip_buffer(fitseq->bitpix, buffer, area);
	return 0;
}

/* create a fits sequence with the given name into the given struct */
int fitseq_create_file(const char *filename, fitseq *fitseq, int frame_count) {
	g_unlink(filename); /* Delete old file if it already exists */
	fitseq_init_struct(fitseq);

	int status = 0;
	if (siril_fits_create_diskfile(&fitseq->fptr, filename, &status)) { /* create new FITS file */
		report_fits_error(status);
		return 1;
	}

	fitseq->filename = strdup(filename);
	fitseq->frame_count = frame_count;
	fitseq->writer = malloc(sizeof(struct seqwriter_data));
	fitseq->writer->write_image_hook = fitseq_write_image_for_writer;
	fitseq->writer->sequence = fitseq;
	siril_debug_print("Successfully created the FITS sequence file %s, for %d images, waiting for data\n",
			fitseq->filename, fitseq->frame_count);

	start_writer(fitseq->writer, frame_count);
	return 0;
}

static int fitseq_write_image_for_writer(struct seqwriter_data *writer, fits *image, int index) {
	fitseq *fitseq = (struct fits_sequence *)writer->sequence;
	int status = 0;
	if (fits_create_img(fitseq->fptr, image->bitpix,
				image->naxis, image->naxes, &status)) {
		report_fits_error(status);
		return 1;
	}

	image->fptr = fitseq->fptr;

	if (com.pref.comp.fits_enabled) {
		status = siril_fits_compress(image);
		if (status) {
			report_fits_error(status);
			return 1;
		}
	}

	return save_opened_fits(image); // warning: will change HDU
}

/* expected images (if a frame count is given on creation) MUST be notified in
 * all cases, even with a NULL image if there is in fact no image to write for
 * the index
 */
int fitseq_write_image(fitseq *fitseq, fits *image, int index) {
	if (!fitseq->fptr) {
		siril_log_color_message(_("Cannot save image in sequence not opened for writing\n"), "red");
		return 1;
	}
	siril_debug_print("FITS sequence %s pending image save %d\n", fitseq->filename, index);
	return seqwriter_append_write(fitseq->writer, image, index);
}

static int fitseq_destroy(fitseq *fitseq) {
	int retval = 0;
	if (fitseq->writer) {
		retval = stop_writer(fitseq->writer);
		free(fitseq->writer);
		fitseq->writer = NULL;
	}
	int status = 0;
	fits_close_file(fitseq->fptr, &status);
	if (fitseq->filename)
		free(fitseq->filename);
	return retval;
}

void fitseq_close_and_delete_file(fitseq *fitseq) {
	char *filename = fitseq->filename;
	fitseq->filename = NULL;
	fitseq_destroy(fitseq);
	siril_log_message(_("Removing failed FITS sequence file: %s\n"), filename);
	g_unlink(filename);
}

int fitseq_close_file(fitseq *fitseq) {
	return fitseq_destroy(fitseq);
}

// to call after open to read with several threads in the file
int fitseq_prepare_for_multiple_read(fitseq *fitseq) {
#ifdef _OPENMP
	if (fitseq->thread_fptr) return 0;
	if (!fitseq->is_mt_capable) return 0;
	guint num_proc = g_get_num_processors();
	fitseq->thread_fptr = malloc(num_proc * sizeof(fitsfile *));
	for (guint i = 0; i < num_proc; i++) {
		int status = 0;
		if (siril_fits_open_diskfile(&fitseq->thread_fptr[i], fitseq->filename, READONLY, &status)) {
			report_fits_error(status);
			return -1;
		}
	}
	siril_debug_print("initialized FITS sequence fd for %d threads reading\n", num_proc);
#endif
	return 0;
}

int fitseq_multiple_close(fitseq *fitseq) {
	int retval = 0;
#ifdef _OPENMP
	if (!fitseq->thread_fptr) return 0;
	guint num_proc = g_get_num_processors();
	for (guint i = 0; i < num_proc; i++) {
		int status = 0;
		fits_close_file(fitseq->thread_fptr[i], &status);
		if (status)
			retval = 1;
	}
	free(fitseq->thread_fptr);
	fitseq->thread_fptr = NULL;
	siril_debug_print("closing FITS sequence fd for %d threads\n", num_proc);
#endif
	return retval;
}


