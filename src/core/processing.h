#ifndef _PROCESSING_H_
#define _PROCESSING_H_

#include "sequence_filtering.h"

/**
 *
 * \file processing.h
 * \brief
 *
 */

/** Main structure of the generic function */
struct generic_seq_args {
	/** sequence that will be processed */
	sequence *seq;

	/** process a partial image read from area instead of full-frame reading */
	gboolean partial_image;
	/** area of the partial image */
	rectangle area;
	/** in case of partial image reading, only one layer is read too */
	int layer_for_partial;
	/** in case of partial, we may use registration data to move the area */
	gboolean regdata_for_partial;
	/** flag to get photometry data */
	gboolean get_photometry_data_for_partial;

	/** filtering the images from the sequence, maybe we don't want them all */
	seq_image_filter filtering_criterion;
	/** filtering parameter */
	double filtering_parameter;
	/** if already known, the number of images after filtering, for smoother
	 *  progress report. < 1 is unknown */
	int nb_filtered_images;

	/** function called before iterating through the sequence */
	int (*prepare_hook)(struct generic_seq_args *);
	/** function called for each image with image index in sequence, number
	 *  of image currently processed and the image, area if partial */
	int (*image_hook)(struct generic_seq_args *, int, int, fits *, rectangle *);
	/** saving the processed image, the one passed to the image_hook, so
	 * in-place editing, if the image_hook succeeded. Only used if
	 * has_output, set to NULL to get default behaviour */
	int (*save_hook)(struct generic_seq_args *, int, int, fits *);
	/** function called after iterating through the sequence, or on
	 * clean-up, even in case of error */
	int (*finalize_hook)(struct generic_seq_args *);

	/** idle function to register at the end. If NULL, the default ending
	 *  that stops the thread is queued. Return false for single execution.
	 *  It should free its argument. */
	GSourceFunc idle_function;
	/** retval, useful for the idle_function, set by the worker */
	int retval;

	/** if false, ignore image_hook errors, unselect failing image from the
	    sequence and continue processing other images */
	gboolean stop_on_error;

	/** string description for progress and logs */
	const char *description;

	/** some processing may create a new image sequence */
	gboolean has_output;
	/** output files: prefix for the new sequence and automatic loading */
	const char *new_seq_prefix;
	/** flag to load or not a new sequence */
	gboolean load_new_sequence;
	/** flag to force output to be SER file */
	gboolean force_ser_output;
	/** new output SER if seq->type == SEQ_SER or force_ser_output (internal) */
	struct ser_struct *new_ser;

	/** user data: pointer to operation-specific data */
	void *user;

	/** if the generic sequence processing is run from an existing thread,
	 * the idle function is not executed in the GTK+ main thread but in this
	 * same thread. If this is false, the generic idle function is run. */
	gboolean already_in_a_thread;
	/** activate parallel execution */
	gboolean parallel;
#ifdef _OPENMP
	/** for in-hook synchronization (internal init, public use) */
	omp_lock_t lock;
#endif
};

gpointer generic_sequence_worker(gpointer p);
gboolean end_generic_sequence(gpointer p);

int ser_prepare_hook(struct generic_seq_args *args);
int ser_finalize_hook(struct generic_seq_args *args);
int generic_save(struct generic_seq_args *, int, int, fits *);

void start_in_new_thread(gpointer(*f)(gpointer p), gpointer p);
gpointer waiting_for_thread();
void stop_processing_thread();
void set_thread_run(gboolean b);
gboolean get_thread_run();

void start_in_reserved_thread(gpointer (*f)(gpointer), gpointer p);
gboolean reserve_thread();
void unreserve_thread();

gboolean end_generic(gpointer arg);
guint siril_add_idle(GSourceFunc idle_function, gpointer data);

#endif
