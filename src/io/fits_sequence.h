#ifndef _FITS_SEQUENCE_H
#define _FITS_SEQUENCE_H

#include <fitsio.h>
#include <glib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "core/siril.h"

struct fits_sequence {
	char *filename;
	int bitpix;		// bitpix of the sequence
	int orig_bitpix;	// bitpix of the data in the file
	long naxes[3];		// size of each dimension
	int frame_count;	// number of frames in the sequence, of HDU in the FITS
	int *hdu_index;		// HDU index for each frame (cfitsio indices, start at 1)

	fitsfile *fptr;		// cfitsio file descriptor.

	gboolean is_mt_capable;	// cfitsio has the option to use multi-threading
#ifdef _OPENMP
	fitsfile **thread_fptr;		// cfitsio file descriptor for each thread read only
#endif
	GThread *write_thread;		// reads a script and executes its commands
	GAsyncQueue *writes_queue;	// the write tasks queue
	gboolean *frame_written;	// keeping track of holes in the sequence
};

typedef struct fits_sequence fitseq;

void fitseq_init_struct(fitseq *fitseq);
int fitseq_is_fitseq(const char *filename, int *frames);

int fitseq_open(const char *filename, fitseq *fitseq);
int fitseq_read_frame(fitseq *fitseq, int index, fits *dest, gboolean force_float, int thread);
int fitseq_read_partial_fits(fitseq *fitseq, int layer, int index, fits *dest, const rectangle *area, gboolean do_photometry, int thread);
int fitseq_read_partial(fitseq *fitseq, int layer, int index, void *buffer, const rectangle *area, int thread);

int fitseq_create_file(const char *filename, fitseq *fitseq, int frame_count);
int fitseq_write_image(fitseq *fitseq, fits *image, int index);
void fitseq_close_and_delete_file(fitseq *fitseq);
void fitseq_close_file(fitseq *fitseq);

int fitseq_prepare_for_multiple_read(fitseq *fitseq);
int fitseq_multiple_close(fitseq *fitseq);

void fitseq_set_max_active_blocks(int max);
void fitseq_wait_for_memory();

#endif
