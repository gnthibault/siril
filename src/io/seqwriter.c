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

#include "seqwriter.h"
#include "core/siril_log.h"
#include "io/image_format_fits.h"

typedef enum {
	SEQ_OK = 0,
	SEQ_WRITE_ERROR,
	SEQ_INCOMPLETE
} seq_error;

static void init_images(struct seqwriter_data *writer, fits *example) {
	writer->bitpix = example->bitpix;
	memcpy(writer->naxes, example->naxes, sizeof writer->naxes);
}

struct _pending_write {
	fits *image;
	int index;
};

#define ABORT_TASK ((void *)0x66)

static void notify_data_freed(struct seqwriter_data *writer, int index);

int seqwriter_append_write(struct seqwriter_data *writer, fits *image, int index) {
	struct _pending_write *newtask = malloc(sizeof(struct _pending_write));
	if (!newtask)
		return -1;
	newtask->image = image;
	newtask->index = index;
	g_async_queue_push(writer->writes_queue, newtask);
	return 0;
}

static void *write_worker(void *a) {
	struct seqwriter_data *writer = (struct seqwriter_data *)a;
	seq_error retval = SEQ_OK;
	int nb_frames_written = 0, current_index = 0;
	GList *next_images = NULL;

	do {
		struct _pending_write *task = NULL;
		GList *stored;
		for (stored = next_images; stored != NULL; stored = stored->next) {
			struct _pending_write *stored_task = (struct _pending_write *)stored->data;
			if (stored_task->index == current_index) {
				task = stored_task;
				next_images = g_list_delete_link(next_images, stored);
				siril_debug_print("writer: image %d obtained from waiting list\n", task->index);
				break;
			}
		}

		if (!task) {	// if not in the waiting list, try to get it from processing threads
			do {
				siril_debug_print("writer: waiting for message %d\n", current_index);
				task = g_async_queue_pop(writer->writes_queue);	// blocking
				if (task == ABORT_TASK) {
					siril_debug_print("writer: abort message\n");
					retval = SEQ_INCOMPLETE;
					break;
				}

				if (writer->bitpix && task->image &&
						(memcmp(task->image->naxes, writer->naxes, sizeof writer->naxes) ||
						 task->image->bitpix != writer->bitpix)) {
					siril_log_color_message(_("Cannot add an image with different properties to an existing sequence.\n"), "red");
					retval = SEQ_WRITE_ERROR;
					break;
				}

				if (task->index >= 0 && task->index != current_index) {
					siril_debug_print("writer: image %d put stored for later use\n", task->index);
					next_images = g_list_append(next_images, task);
					task = NULL;
				}
				else siril_debug_print("writer: image %d received\n", task->index);
			} while (!task);
		}
		if (!task)
			continue;
		if (retval == SEQ_INCOMPLETE)
			break;
		if (!task->image) {
			// failed image, hole in sequence, skip it
			siril_debug_print("writer: skipping image %d\n", task->index);
			notify_data_freed(writer, task->index);
			current_index++;
			writer->frame_count--;
			continue;
		}

		// from here on, we have a valid task and we will write an image
		if (!writer->bitpix)
			init_images(writer, task->image);

		siril_log_message(_("writer: Saving image %d, %ld layer(s), %ux%u pixels, %d bits\n"),
				task->index, task->image->naxes[2],
				task->image->rx, task->image->ry,
				task->image->type == DATA_FLOAT ? 32 : 16);

		retval = writer->write_image_hook(writer, task->image, nb_frames_written);
		clearfits(task->image);

		if (retval != SEQ_WRITE_ERROR) {
			notify_data_freed(writer, task->index);
			nb_frames_written++;
			current_index++;
		}
		free(task->image);
		free(task);
	} while (retval == SEQ_OK &&
			(writer->frame_count <= 0 || nb_frames_written < writer->frame_count));

	if (retval == SEQ_INCOMPLETE) {
		if (writer->frame_count <= 0) {
			writer->frame_count = nb_frames_written;
			retval = SEQ_OK;
			siril_log_message(_("Saved %d images in the sequence\n"), nb_frames_written);
		} else {
			siril_debug_print("writer: write aborted, expected %d images, got %d.\n",
					writer->frame_count, nb_frames_written);
		}
	}

	siril_debug_print("writer exits with retval %d (0: ok, 1: error, 2: incomplete)\n", retval);
	return GINT_TO_POINTER(retval);
}

/* frame_count can be unknown and nil or negative, otherwise, providing it will
 * give clearer output on completion of the sequence */
void start_writer(struct seqwriter_data *writer, int frame_count) {
	g_assert(writer->write_image_hook);
	g_assert(writer->sequence);
	writer->bitpix = 0;
	writer->naxes[0] = 0;
	writer->frame_count = frame_count;
	writer->writes_queue = g_async_queue_new();
	writer->write_thread = g_thread_new("writer", write_worker, writer);
}

int stop_writer(struct seqwriter_data *writer) {
	int retval = 0;
	if (writer->write_thread) {
		g_async_queue_push(writer->writes_queue, ABORT_TASK);
		siril_debug_print("writer thread notified, waiting for exit...\n");
		gpointer ret = g_thread_join(writer->write_thread);
		writer->write_thread = NULL;
		g_async_queue_unref(writer->writes_queue);
		retval = GPOINTER_TO_INT(ret);
		siril_debug_print("writer thread joined (retval: %d)\n", retval);
		seqwriter_set_max_active_blocks(0); // wake-up the callers
	}
	return retval;
}

/* FITS cannot be written by several threads at the same time. We still want to
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
 * limits. In case the memory is larger than what the number of thread can
 * support, the threads won't be blocked until too many images are pending
 * write.
 * The code below counts the number of active memory blocks and provides a
 * waiting function.
 */

static int nb_blocks_active, configured_max_active_blocks;
static int nb_outputs = 1;
static GCond pool_cond;
static GMutex pool_mutex;

/* here we keep the latest frame index written for each output sequence, to
 * synchronize in case of several output sequences for a single processing */
struct _outputs_struct {
	void *seq;
	int index;
};
static struct _outputs_struct *outputs;

void seqwriter_set_max_active_blocks(int max) {
	siril_log_message(_("Number of images allowed in the FITS write queue: %d (zero or less is unlimited)\n"), max);
	configured_max_active_blocks = max;
	nb_blocks_active = 0;
}

void seqwriter_wait_for_memory() {
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

static int get_output_for_seq(void *seq) {
	for (int i = 0; i < nb_outputs; i++) {
		if (!outputs[i].seq) {
			outputs[i].seq = seq;
			outputs[i].index = -1;
			return i;
		}
		if (outputs[i].seq == seq)
			return i;
	}
	siril_debug_print("### seqwriter get_output_for_seq: not found! should never happen ###\n");
	return -1;
}

static gboolean all_outputs_to_index(int index) {
	for (int i = 0; i < nb_outputs; i++) {
		if (!outputs[i].seq)
			return FALSE;
		if (outputs[i].index < index)
			return FALSE;
	}
	siril_debug_print("\tgot all outputs notified for index %d, signaling\n", index);
	return TRUE;
}

static void notify_data_freed(struct seqwriter_data *writer, int index) {
	g_mutex_lock(&pool_mutex);
	if (nb_outputs > 1) {
		int output_num = get_output_for_seq(writer->sequence);
		if (outputs[output_num].index + 1 != index) {
			fprintf(stderr, "inconsistent index in memory management (%d for expected %d)\n",
					outputs[output_num].index + 1, index);
		}
		outputs[output_num].index = index;
		if (!all_outputs_to_index(index)) {
			g_mutex_unlock(&pool_mutex);
			return;
		}
	}

	nb_blocks_active--;
	g_cond_signal(&pool_cond);
	g_mutex_unlock(&pool_mutex);
}

void seqwriter_set_number_of_outputs(int number_of_outputs) {
	siril_debug_print("seqwriter number of outputs: %d\n", number_of_outputs);
	nb_outputs = number_of_outputs;
	if (number_of_outputs > 1) {
		outputs = calloc(number_of_outputs, sizeof(struct _outputs_struct));
	} else {
		if (outputs)
			free(outputs);
		outputs = NULL;
	}
}
