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

#include "fits_sequence.h"
#include "image_format_fits.h"
#include "gui/progress_and_log.h"
#include "core/siril_log.h"

static void *write_worker(void *a);
static void notify_data_freed();

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
	fitseq->write_thread = NULL;
	fitseq->writes_queue = NULL;
	fitseq->frame_written= NULL;
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
	if (do_photometry)
		fit_get_photometry_data(dest);

	int status = 0;
	if (fits_movabs_hdu(fptr, fitseq->hdu_index[index], NULL, &status)) {
		report_fits_error(status);
		return -1;
	}

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
	if (frame_count > 0)
		fitseq->frame_written = calloc(frame_count, sizeof(gboolean));

	fitseq->writes_queue = g_async_queue_new();
	fitseq->write_thread = g_thread_new("processing", write_worker, fitseq);

	siril_debug_print("Successfully created the FITS sequence file %s, for %d images, waiting for data\n",
			fitseq->filename, fitseq->frame_count);

	return 0;
}

static void init_images(fitseq *fitseq, fits *example) {
	fitseq->bitpix = example->bitpix;
	memcpy(fitseq->naxes, example->naxes, sizeof fitseq->naxes);
}

struct _pending_write {
	fits *image;
	int index;
};

#define ABORT_TASK ((void *)0x66)

int fitseq_write_image(fitseq *fitseq, fits *image, int index) {
	if (!fitseq->fptr) {
		siril_log_color_message(_("Cannot save image in sequence not opened for writing\n"), "red");
		return 1;
	}
	siril_debug_print("FITS sequence %s pending image save %d\n", fitseq->filename, index);

	struct _pending_write *newtask = malloc(sizeof(struct _pending_write));
	newtask->image = image;
	newtask->index = index;
	g_async_queue_push(fitseq->writes_queue, newtask);
	return 0;
}

typedef enum {
	FITSEQ_OK = 0,
	FITSEQ_WRITE_ERROR,
	FITSEQ_INCOMPLETE
} fitseq_error;

static void *write_worker(void *a) {
	fitseq *fitseq = (struct fits_sequence *)a;
	fitseq_error retval = FITSEQ_OK;
	int nb_frames_written = 0, status;
	GList *next_images = NULL;

	do {
		struct _pending_write *task = NULL;
		GList *stored;
		for (stored = next_images; stored != NULL; stored = stored->next) {
			struct _pending_write *stored_task = (struct _pending_write *)stored->data;
			if (stored_task->index == nb_frames_written) {
				task = stored_task;
				next_images = g_list_delete_link(next_images, stored);
				siril_debug_print("fitseq write: image %d obtained from waiting list\n", task->index);
				break;
			}
		}

		if (!task) {	// if not in the waiting list, try to get it from processing threads
			do {
			siril_debug_print("fitseq write: waiting for message %d\n", nb_frames_written);
				task = g_async_queue_pop(fitseq->writes_queue);	// blocking
				if (fitseq->bitpix && (memcmp(task->image->naxes, fitseq->naxes, sizeof fitseq->naxes) ||
							task->image->bitpix != fitseq->bitpix)) {
					siril_log_color_message(_("Cannot add an image with different properties to an existing sequence.\n"), "red");
					retval = FITSEQ_WRITE_ERROR;
					break;
				}

				if (task == ABORT_TASK) {
					siril_debug_print("fitseq write: abort message\n");
					retval = FITSEQ_INCOMPLETE;
					break;
				}
				if (task->index >= 0 && task->index != nb_frames_written) {
					siril_debug_print("fitseq write: image %d put stored for later use\n", task->index);
					next_images = g_list_append(next_images, task);
					task = NULL;
				}
			} while (!task);
			siril_debug_print("fitseq write: image %d received\n", task->index);
		}
		if (retval == FITSEQ_INCOMPLETE)
			break;

		if (!fitseq->bitpix)
			init_images(fitseq, task->image);

		status = 0;
		if (fits_create_img(fitseq->fptr, task->image->bitpix,
					task->image->naxis, task->image->naxes, &status)) {
			report_fits_error(status);
			retval = FITSEQ_WRITE_ERROR;
			break;
		}

		siril_log_message(_("fitseq write: Saving FITS image %d, %ld layer(s), %ux%u pixels, %d bits\n"),
				task->index, task->image->naxes[2],
				task->image->rx, task->image->ry,
				task->image->type == DATA_FLOAT ? 32 : 16);
		task->image->fptr = fitseq->fptr;
		if (save_opened_fits(task->image)) // warning: will change HDU
			retval = FITSEQ_WRITE_ERROR;
		clearfits(task->image);

		if (retval != FITSEQ_WRITE_ERROR) {
			notify_data_freed();
			fitseq->frame_written[task->index] = TRUE;
			nb_frames_written++;
		}
		free(task->image);
		free(task);
	} while (retval == FITSEQ_OK &&
			(fitseq->frame_count <= 0 || nb_frames_written < fitseq->frame_count));

	if (retval == FITSEQ_OK) {
		int nb_hdu;
		status = 0;
		fits_get_num_hdus(fitseq->fptr, &nb_hdu, &status);
		// this should never happen
		if (nb_hdu != nb_frames_written) {
			siril_log_color_message(_("Frame number mismatch at the end of sequence saving (%d saved for %d written)\n"), "red", nb_hdu, nb_frames_written);
			retval = 1;
		}
	}
	else if (retval == FITSEQ_INCOMPLETE) {
		if (fitseq->frame_count < 0) {
			fitseq->frame_count = nb_frames_written;
			retval = FITSEQ_OK;
			siril_log_message(_("Successfully saved FITS sequence with %d images\n"), nb_frames_written);
		}
		// we don't know here if it's cancellation or end of
		// processing, so we can't call compact
	}

	siril_debug_print("fitseq writer exits with retval %d (0: ok, 1: error, 2: incomplete)\n", retval);
	return GINT_TO_POINTER(retval);
}

/*int fitseq_append_image_from_disk(fitseq *fitseq, fits *image) {
	int retval = 0;
	if (!image->fptr) {
		return -1;
	}

	fits_copy_hdu(image->fptr, fitseq->fptr, 0, &retval);
	return retval;
}*/

int fitseq_compact_file(fitseq *fitseq) {
	if (!fitseq->frame_written || fitseq->frame_count <= 0)
		return FITSEQ_OK;
	int nb_hdu;
	int status = 0;
	fits_get_num_hdus(fitseq->fptr, &nb_hdu, &status);
	siril_debug_print("compacting fitseq: %d HDU for %d frame count\n", nb_hdu, fitseq->frame_count);

	for (int i = 0, j = 0; i < fitseq->frame_count; i++) {
		if (!fitseq->frame_written[i]) {
			siril_debug_print("removing failed image %d from sequence\n", i);
			status = 0;
			if (fits_movabs_hdu(fitseq->fptr, j + 1, NULL, &status)) {
				siril_log_color_message(_("Failed to remove failed images from the sequence\n"), "red");
				return FITSEQ_WRITE_ERROR;
			}

			status = 0;
			if (fits_delete_hdu(fitseq->fptr, NULL, &status)) {
				siril_log_color_message(_("Failed to remove failed images from the sequence\n"), "red");
				return FITSEQ_WRITE_ERROR;
			}
		}
		else j++;
	}
	status = 0;
	fits_get_num_hdus(fitseq->fptr, &nb_hdu, &status);
	fitseq->frame_count = nb_hdu;
	siril_debug_print("after fitseq compaction: %d frame count\n", fitseq->frame_count);
	return FITSEQ_OK;
}

static int fitseq_destroy(fitseq *fitseq, gboolean compact) {
	int retval = 0;
	if (fitseq->write_thread) {
		g_async_queue_push(fitseq->writes_queue, ABORT_TASK);
		siril_debug_print("fitseq writing thread notified, waiting for exit...\n");
		gpointer ret = g_thread_join(fitseq->write_thread);
		fitseq->write_thread = NULL;
		g_async_queue_unref(fitseq->writes_queue);
		retval = GPOINTER_TO_INT(ret);
		siril_debug_print("fitseq writing thread joined (retval: %d)\n", retval);
		if (retval == FITSEQ_INCOMPLETE && compact)
			retval = fitseq_compact_file(fitseq);
		fitseq_set_max_active_blocks(0); // wake-up the callers
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
	fitseq_destroy(fitseq, FALSE);
	siril_log_message(_("Removing failed FITS sequence file: %s\n"), filename);
	g_unlink(filename);
}

void fitseq_close_file(fitseq *fitseq) {
	fitseq_destroy(fitseq, TRUE);
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

/* TO BE GENERALIZED TO ALL PROCESSING FUNCTIONS
 * FITS cannot be written by several threads at the same time. We still want to
 * read and process files in parallel and save the results into a FITS
 * sequence, so instead of writing in the file from each processing thread, we
 * queue the writes and a single thread, launched manually with the
 * write_worker function, writes to the file.
 * The problem with that is memory management. In most algorithms, we limit the
 * number of threads to match memory requirements, because each thread needs
 * memory to handle the image data. With the writes queued, memory is not freed
 * when the processing ends, but the thread is ready to process more, hence
 * allocate more. We have to pause the processing until the writing thread has
 * saved a result and freed the data, otherwise siril will go out of the memory
 * limits. In case the memory it larger than what the number of thread can
 * support, the threads won't be blocked until too many images are pending
 * write.
 * The code below counts the number of active memory blocks and provides a
 * waiting function.
 */

static int nb_blocks_active, configured_max_active_blocks;
static GCond pool_cond;
static GMutex pool_mutex;

void fitseq_set_max_active_blocks(int max) {
	siril_log_message(_("Number of images allowed in the FITS write queue: %d\n"), max);
	configured_max_active_blocks = max;
	nb_blocks_active = 0;
}

void fitseq_wait_for_memory() {
	if (configured_max_active_blocks <= 0)
		return;
	siril_debug_print("entering the wait function\n");
	g_mutex_lock(&pool_mutex);
	while (nb_blocks_active >= configured_max_active_blocks) {
		siril_debug_print("  waiting for free memory slot (%d active)\n", nb_blocks_active);
		g_cond_wait(&pool_cond, &pool_mutex);
	}
	nb_blocks_active++;
	siril_debug_print("got the slot!\n");
	g_mutex_unlock(&pool_mutex);
}

static void notify_data_freed() {
	g_mutex_lock(&pool_mutex);
	nb_blocks_active--;
	g_cond_signal(&pool_cond);
	g_mutex_unlock(&pool_mutex);
}

