#ifndef SEQ_WRITER_H
#define SEQ_WRITER_H

#include <glib.h>
#include "core/siril.h"

struct seqwriter_data {
	int bitpix;		// bitpix of the sequence
	long naxes[3];		// size of each dimension
	int frame_count;	// number of frames in the sequence, of HDU in the FITS

	GThread *write_thread;		// reads a script and executes its commands
	GAsyncQueue *writes_queue;	// the write tasks queue
	gint failed;

	int (*write_image_hook)(struct seqwriter_data *writer, fits *image, int index);
	void *sequence;
};

void start_writer(struct seqwriter_data *writer, int frame_count);
int stop_writer(struct seqwriter_data *writer, gboolean aborting);
int seqwriter_append_write(struct seqwriter_data *writer, fits *image, int index);

void seqwriter_set_max_active_blocks(int max);
void seqwriter_wait_for_memory();
void seqwriter_release_memory();
void seqwriter_set_number_of_outputs(int number_of_outputs);

#endif
