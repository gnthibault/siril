/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2021 team free-astro (see more in AUTHORS file)
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

#include <assert.h>
#include <string.h>
#include <glib.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/sequence_filtering.h"
#include "core/OS_utils.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/seqwriter.h"
#include "io/fits_sequence.h"
#include "io/image_format_fits.h"
#include "algos/statistics.h"

// called in start_in_new_thread only
// works in parallel if the arg->parallel is TRUE for FITS or SER sequences
gpointer generic_sequence_worker(gpointer p) {
	struct generic_seq_args *args = (struct generic_seq_args *)p;
	struct timeval t_start, t_end;
	int frame;	// output frame index
	int input_idx;	// index of the frame being processed in the sequence
	int *index_mapping = NULL;
	int nb_frames, excluded_frames = 0, progress = 0;
	float nb_framesf;
	int abort = 0;	// variable for breaking out of loop

	assert(args);
	assert(args->seq);
	assert(args->image_hook);
	set_progress_bar_data(NULL, PROGRESS_RESET);
	gettimeofday(&t_start, NULL);

	if (args->nb_filtered_images > 0)
		nb_frames = args->nb_filtered_images;
	else {
		nb_frames = compute_nb_filtered_images(args->seq, args->filtering_criterion, args->filtering_parameter);
		args->nb_filtered_images = nb_frames;
		if (nb_frames <= 0) {
			siril_log_message(_("No image selected for processing, aborting\n"));
			args->retval = 1;
			goto the_end;
		}
	}
	nb_framesf = (float)nb_frames + 0.3f;	// leave margin for rounding errors and post processing
	args->retval = 0;

#ifdef _OPENMP
	if (args->max_thread < 1) {
		if (args->compute_mem_limits_hook)
			args->max_thread = args->compute_mem_limits_hook(args, FALSE);
		else args->max_thread = seq_compute_mem_limits(args, FALSE);
	}
	if (args->max_thread < 1) {
		args->retval = 1;
		goto the_end;
	}
	if (args->compute_mem_limits_hook)
		siril_log_message(_("%s: with the current memory and thread limits, up to %d thread(s) can be used\n"),
				args->description, args->max_thread);
#endif

	if (args->prepare_hook && args->prepare_hook(args)) {
		siril_log_message(_("Preparing sequence processing failed.\n"));
		args->retval = 1;
		goto the_end;
	}

	/* We have a sequence in which images can be filtered out. In order to
	 * distribute the workload fairly among all threads, the main iteration
	 * should not be on the list of images of the sequence, but on the
	 * actual list of selected images.
	 * Here we create this map of images to be processed, each cell of the
	 * array providing the image number in the input sequence. It cannot be
	 * done in parallel.
	 * This is mandatory for SER contiguous output. */
	if (args->filtering_criterion) {
		index_mapping = malloc(nb_frames * sizeof(int));
		if (!index_mapping) {
			PRINT_ALLOC_ERR;
			args->retval = 1;
			goto the_end;
		}
		for (input_idx = 0, frame = 0; input_idx < args->seq->number; input_idx++) {
			if (!args->filtering_criterion(args->seq, input_idx, args->filtering_parameter)) {
				continue;
			}
			index_mapping[frame++] = input_idx;
		}
		if (frame != nb_frames) {
			siril_log_message(_("Output index mapping failed (%d/%d).\n"), frame, nb_frames);
			args->retval = 1;
			goto the_end;
		}
	}

	if (args->has_output && !args->partial_image) {
		gint64 size;
		if (args->compute_size_hook)
			size = args->compute_size_hook(args, nb_frames);
		else size = seq_compute_size(args->seq, nb_frames, args->output_type);
		if (test_available_space(size)) {
			args->retval = 1;
			goto the_end;
		}
	}

	/* Output print of algorithm description */
	if (args->description) {
		gchar *desc = g_strdup_printf(_("%s: processing...\n"), args->description);
		siril_log_color_message(desc, "green");
		g_free(desc);
	}

	gboolean have_seqwriter = args->has_output &&
		((args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) ||
		 (args->force_ser_output || args->seq->type == SEQ_SER));
#ifdef _OPENMP
	omp_init_lock(&args->lock);
	if (have_seqwriter)
		omp_set_schedule(omp_sched_dynamic, 1);
	else omp_set_schedule(omp_sched_guided, 0);
#ifdef HAVE_FFMS2
	// we don't want to enable parallel processing for films, as ffms2 is not thread-safe
#pragma omp parallel for num_threads(args->max_thread) private(input_idx) schedule(runtime) \
	if(args->seq->type != SEQ_AVI && args->parallel && (args->seq->type == SEQ_SER || fits_is_reentrant()))
#else
#pragma omp parallel for num_threads(args->max_thread) private(input_idx) schedule(runtime) \
	if(args->parallel && (args->seq->type == SEQ_SER || fits_is_reentrant()))
#endif // HAVE_FFMS2
#endif // _OPENMP
	for (frame = 0; frame < nb_frames; frame++) {
		if (abort) continue;

		fits *fit = calloc(1, sizeof(fits));
		char filename[256];
		rectangle area = { .x = args->area.x, .y = args->area.y,
			.w = args->area.w, .h = args->area.h };

		if (!get_thread_run()) {
			abort = 1;
			continue;
		}
		if (index_mapping)
			input_idx = index_mapping[frame];
		else input_idx = frame;

		if (!seq_get_image_filename(args->seq, input_idx, filename)) {
			abort = 1;
			continue;
		}

		int thread_id = -1;
#ifdef _OPENMP
		thread_id = omp_get_thread_num();
		if (have_seqwriter) {
			seqwriter_wait_for_memory();
			if (abort) {
				seqwriter_release_memory();
				continue;
			}
		}
#endif

		if (args->partial_image) {
			// if we run in parallel, it will not be the same for all
			// and we don't want to overwrite the original anyway
			if (args->regdata_for_partial) {
				int shiftx = roundf_to_int(args->seq->regparam[args->layer_for_partial][input_idx].shiftx);
				int shifty = roundf_to_int(args->seq->regparam[args->layer_for_partial][input_idx].shifty);
				area.x -= shiftx;
				area.y += shifty;
			}

			// args->area may be modified in hooks
			enforce_area_in_image(&area, args->seq);
			if (seq_read_frame_part(args->seq,
						args->layer_for_partial,
						input_idx, fit, &area,
						args->get_photometry_data_for_partial, thread_id))
			{
				abort = 1;
				clearfits(fit);
				free(fit);
				continue;
			}
			/*char tmpfn[100];	// this is for debug purposes
			  sprintf(tmpfn, "/tmp/partial_%d.fit", input_idx);
			  savefits(tmpfn, fit);*/
		} else {
			// image is obtained bottom to top here, while it's in natural order for partial images!
			if (seq_read_frame(args->seq, input_idx, fit, args->force_float, thread_id)) {
				abort = 1;
				clearfits(fit);
				free(fit);
				continue;
			}
		}

		if (args->image_hook(args, frame, input_idx, fit, &area)) {
			if (args->stop_on_error)
				abort = 1;
			else {
				//args->seq->imgparam[frame].incl = FALSE;
#ifdef _OPENMP
#pragma omp critical
#endif
				{
					//if (args->nb_filtered_images > 0)
					//	args->nb_filtered_images--;
					//args->seq->selnum--;
					excluded_frames++;
				}
			}
			clearfits(fit);
			free(fit);
			// for seqwriter, we need to notify the failed frame
			if (have_seqwriter) {
				int retval;
				if (args->save_hook)
					retval = args->save_hook(args, frame, input_idx, NULL);
				else retval = generic_save(args, frame, input_idx, NULL);
				if (retval)
					abort = 1;
			}
			continue;
		}

		if (args->has_output) {
			int retval;
			if (args->save_hook)
				retval = args->save_hook(args, frame, input_idx, fit);
			else retval = generic_save(args, frame, input_idx, fit);
			if (retval) {
				abort = 1;
				clearfits(fit);
				free(fit);
				continue;
			}
		} else {
			/* save stats that may have been computed for the first
			 * time, but if fit has been modified for the new
			 * sequence, we shouldn't save it for the old one.
			 */
			save_stats_from_fit(fit, args->seq, input_idx);
		}

		if (!have_seqwriter) {
			clearfits(fit);
			free(fit);
		}

#ifdef _OPENMP
#pragma omp atomic
#endif
		progress++;
		gchar *msg = g_strdup_printf(_("%s. Processing image %d (%s)"), args->description, input_idx + 1, filename);
		set_progress_bar_data(msg, (float)progress / nb_framesf);
		g_free(msg);
	}

	if (abort) {
		set_progress_bar_data(_("Sequence processing failed. Check the log."), PROGRESS_RESET);
		siril_log_color_message(_("Sequence processing failed.\n"), "red");
		args->retval = abort;
	}
	else {
		if (excluded_frames) {
			set_progress_bar_data(_("Sequence processing partially succeeded. Check the log."), PROGRESS_RESET);
			siril_log_color_message(_("Sequence processing partially succeeded, with %d images that failed and that were temporarily excluded from the sequence.\n"), "salmon", excluded_frames);
		} else {
			set_progress_bar_data(_("Sequence processing succeeded."), PROGRESS_RESET);
			siril_log_color_message(_("Sequence processing succeeded.\n"), "green");
		}
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}

#ifdef _OPENMP
	omp_destroy_lock(&args->lock);
#endif
the_end:
	if (index_mapping) free(index_mapping);
	if (args->finalize_hook && args->finalize_hook(args)) {
		siril_log_message(_("Finalizing sequence processing failed.\n"));
		args->retval = 1;
	}

	if (args->already_in_a_thread) {
		if (args->idle_function)
			args->idle_function(args);
	} else {
		if (args->idle_function)
			siril_add_idle(args->idle_function, args);
		else siril_add_idle(end_generic_sequence, args);
	}
	return GINT_TO_POINTER(args->retval);
}

// defaut idle function (in GTK main thread) to run at the end of the generic sequence processing
gboolean end_generic_sequence(gpointer p) {
	struct generic_seq_args *args = (struct generic_seq_args *) p;

	if (args->has_output && args->load_new_sequence &&
			args->new_seq_prefix && !args->retval) {
		gchar *basename = g_path_get_basename(args->seq->seqname);
		gchar *seqname = g_strdup_printf("%s%s.seq", args->new_seq_prefix, basename);
		check_seq(0);
		update_sequences_list(seqname);
		free(seqname);
		g_free(basename);
	}
	
	free(p);
	return end_generic(NULL);
}

/* If for_writer is false, it computes how many images can be processed in
 * parallel, with regard to how many of them can fit in memory. It returns at
 * most com.max_thread.
 * If for_writer is true, it computes how many images can be stored in the
 * queue. It returns at most 3 times com.max_thread.
 */
int seq_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	unsigned int MB_per_image, MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, args->upscale_ratio, args->force_float, &MB_per_image, &MB_avail);
	if (limit == 0) {
		gchar *mem_per_image = g_format_size_full(MB_per_image * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);

		siril_log_color_message(_("%s: not enough memory to do this operation (%s required per image, %s considered available)\n"),
				"red", args->description, mem_per_image, mem_available);

		g_free(mem_per_image);
		g_free(mem_available);
	} else {
#ifdef _OPENMP
		if (for_writer) {
			int max_queue_size = com.max_thread * 3;
			if (limit > max_queue_size)
				limit = max_queue_size;
		}
		else if (limit > com.max_thread)
			limit = com.max_thread;
#else
		if (!for_writer)
			limit = 1;
		else if (limit > 3)
			limit = 3;
#endif
	}
	return limit;
}

int seq_prepare_hook(struct generic_seq_args *args) {
	int retval = 0;
	g_assert(args->has_output); // don't call this hook otherwise
	if (args->force_ser_output || args->seq->type == SEQ_SER) {
		gchar *dest;
		const char *ptr = strrchr(args->seq->seqname, G_DIR_SEPARATOR);
		if (ptr)
			dest = g_strdup_printf("%s%s.ser", args->new_seq_prefix, ptr + 1);
		else dest = g_strdup_printf("%s%s.ser", args->new_seq_prefix, args->seq->seqname);

		args->new_ser = malloc(sizeof(struct ser_struct));
		if (ser_create_file(dest, args->new_ser, TRUE, args->seq->ser_file)) {
			free(args->new_ser);
			args->new_ser = NULL;
			retval = 1;
		}
		g_free(dest);
	}
	else if (args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) {
		gchar *dest;
		const char *ptr = strrchr(args->seq->seqname, G_DIR_SEPARATOR);
		if (ptr)
			dest = g_strdup_printf("%s%s%s", args->new_seq_prefix, ptr + 1, com.pref.ext);
		else dest = g_strdup_printf("%s%s%s", args->new_seq_prefix, args->seq->seqname, com.pref.ext);

		args->new_fitseq = malloc(sizeof(fitseq));
		if (fitseq_create_file(dest, args->new_fitseq, args->nb_filtered_images)) {
			free(args->new_fitseq);
			args->new_fitseq = NULL;
			retval = 1;
		}
		g_free(dest);
	}
	else return 0;

	if (!retval)
		retval = seq_prepare_writer(args);
	return retval;
}

int seq_prepare_writer(struct generic_seq_args *args) {
	int limit = 0;
	if (args->compute_mem_limits_hook)
		limit = args->compute_mem_limits_hook(args, TRUE);
	else limit = seq_compute_mem_limits(args, TRUE); // the default

	if (limit == 0) {
		if (args->force_ser_output || args->seq->type == SEQ_SER) {
			ser_close_file(args->new_ser);
			free(args->new_ser);
			args->new_ser = NULL;
		}
		else if (args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) {
			fitseq_close_file(args->new_fitseq);
			free(args->new_fitseq);
			args->new_fitseq = NULL;
		}
		return 1;
	}
	seqwriter_set_max_active_blocks(limit);
	return 0;
}

int seq_finalize_hook(struct generic_seq_args *args) {
	int retval = 0;
	g_assert(args->has_output); // don't call this hook otherwise
	if ((args->force_ser_output || args->seq->type == SEQ_SER) && args->new_ser) {
		retval = ser_write_and_close(args->new_ser);
		free(args->new_ser);
	}
	else if ((args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) && args->new_fitseq) {
		retval = fitseq_close_file(args->new_fitseq);
		free(args->new_fitseq);
	}
	return retval;
}

/* In SER, all images must be in a contiguous sequence, so we use the out_index.
 * In FITS sequences, to keep track of image accross processings, we keep the
 * input file number all along (in_index is the index in the sequence, not the name).
 */
int generic_save(struct generic_seq_args *args, int out_index, int in_index, fits *fit) {
	if (args->force_ser_output || args->seq->type == SEQ_SER) {
		return ser_write_frame_from_fit(args->new_ser, fit, out_index);
	} else if (args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) {
		return fitseq_write_image(args->new_fitseq, fit, out_index);
	} else {
		char *dest = fit_sequence_get_image_filename_prefixed(args->seq,
				args->new_seq_prefix, in_index);
		fit->bitpix = fit->orig_bitpix;
		int retval = savefits(dest, fit);
		free(dest);
		return retval;
	}
}

/*****************************************************************************
 *      P R O C E S S I N G      T H R E A D      M A N A G E M E N T        *
 ****************************************************************************/

static gboolean thread_being_waited = FALSE;

// This function is reentrant. The pointer will be freed in the idle function,
// so it must be a proper pointer to an allocated memory chunk.
void start_in_new_thread(gpointer (*f)(gpointer), gpointer p) {
	g_mutex_lock(&com.mutex);
	if (com.run_thread || com.thread) {
		fprintf(stderr, "The processing thread is busy, stop it first.\n");
		g_mutex_unlock(&com.mutex);
		free(p);
		return;
	}

	com.run_thread = TRUE;
	g_mutex_unlock(&com.mutex);
	com.thread = g_thread_new("processing", f, p);
}

void start_in_reserved_thread(gpointer (*f)(gpointer), gpointer p) {
	g_mutex_lock(&com.mutex);
	if (com.thread) {
		fprintf(stderr, "The processing thread is busy, stop it first.\n");
		g_mutex_unlock(&com.mutex);
		free(p);
		return;
	}

	com.run_thread = TRUE;
	g_mutex_unlock(&com.mutex);
	com.thread = g_thread_new("processing", f, p);
}

gpointer waiting_for_thread() {
	gpointer retval = NULL;
	if (com.thread) {
		thread_being_waited = TRUE;
		retval = g_thread_join(com.thread);
	}
	com.thread = NULL;
	thread_being_waited = FALSE;
	set_thread_run(FALSE);	// do it anyway in case of wait without stop
	return retval;
}

void stop_processing_thread() {
	if (com.thread == NULL) {
		siril_debug_print("The processing thread is not running.\n");
		return;
	}

	set_thread_run(FALSE);
	if (!thread_being_waited)
		waiting_for_thread();
}

void set_thread_run(gboolean b) {
	g_mutex_lock(&com.mutex);
	com.run_thread = b;
	g_mutex_unlock(&com.mutex);
}

gboolean get_thread_run() {
	gboolean retval;
	g_mutex_lock(&com.mutex);
	retval = com.run_thread;
	g_mutex_unlock(&com.mutex);
	return retval;
}

// equivalent to atomic get and set if not running
gboolean reserve_thread() {
	gboolean retval;
	g_mutex_lock(&com.mutex);
	retval = !com.run_thread;
	if (retval)
		com.run_thread = TRUE;
	g_mutex_unlock(&com.mutex);
	return retval;
}

void unreserve_thread() {
	set_thread_run(FALSE);
}

/* should be called in a threaded function if nothing special has to be done at the end.
 * siril_add_idle(end_generic, NULL);
 */
gboolean end_generic(gpointer arg) {
	stop_processing_thread();
	
	set_cursor_waiting(FALSE);
	return FALSE;
}

/* wrapper for gdk_threads_add_idle that deactivates idle functions when
 * running in script/console mode */
guint siril_add_idle(GSourceFunc idle_function, gpointer data) {
	if (!com.script && !com.headless)
		return gdk_threads_add_idle(idle_function, data);
	return 0;
}

void wait_for_script_thread() {
	if (com.script_thread)
		g_thread_join(com.script_thread);
	com.script_thread = NULL;
}

void on_processes_button_cancel_clicked(GtkButton *button, gpointer user_data) {
	if (com.thread != NULL)
		siril_log_color_message(_("Process aborted by user\n"), "red");
	com.stop_script = TRUE;
	stop_processing_thread();
	wait_for_script_thread();
}

struct generic_seq_args *create_default_seqargs(sequence *seq) {
	struct generic_seq_args *args = calloc(1, sizeof(struct generic_seq_args));
	args->seq = seq;
	args->filtering_criterion = seq_filter_all;
	args->nb_filtered_images = seq->number;
	args->stop_on_error = TRUE;
	args->upscale_ratio = 1.0;
	args->parallel = TRUE;
	return args;
}
