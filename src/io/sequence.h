#ifndef _SEQUENCE_H_
#define _SEQUENCE_H_

#include <stdint.h>
#include "../core/siril.h"

int	read_single_sequence(char *realname, int imagetype);
int	seqsetnum(int image_number);
int	check_seq(int recompute_stats);
int	check_only_one_film_seq(char* name);
int	seq_check_basic_data(sequence *seq, gboolean load_ref_into_gfit);
int	set_seq(const char *);
char *	seq_get_image_filename(sequence *seq, int index, char *name_buf);
int	seq_read_frame(sequence *seq, int index, fits *dest);
int	seq_read_frame_part(sequence *seq, int layer, int index, fits *dest, const rectangle *area, gboolean do_photometry);
int	seq_load_image(sequence *seq, int index, gboolean load_it);
int64_t seq_compute_size(sequence *seq, int nb_frames);
gboolean check_if_seq_exist(gchar *basename);
int	seq_open_image(sequence *seq, int index);
void	seq_close_image(sequence *seq, int index);
int	seq_opened_read_region(sequence *seq, int layer, int index, WORD *buffer, const rectangle *area);
void	set_fwhm_star_as_star_list(sequence *seq);
char *	fit_sequence_get_image_filename(sequence *seq, int index, char *name_buffer, gboolean add_fits_ext);
char *	fit_sequence_get_image_filename_prefixed(sequence *seq, const char *prefix, int index);
char *	get_possible_image_filename(sequence *seq, int image_number, char *name_buffer);
int	get_index_and_basename(const char *filename, char **basename, int *index, int *fixed);
void	remove_prefixed_sequence_files(sequence *seq, const char *prefix);
void	initialize_sequence(sequence *seq, gboolean is_zeroed);
void	free_sequence(sequence *seq, gboolean free_seq_too);
void	sequence_free_preprocessing_data(sequence *seq);
void	close_sequence(int loading_another);
gboolean sequence_is_loaded();

typedef enum {
	ORIGINAL_FRAME,
	FOLLOW_STAR_FRAME,
	REGISTERED_FRAME
} framing_mode;

/* Procedure signature for sequence processing.
 * Returns < 0 for an error that should stop the processing on the sequence.
 * Other return values are not used.
 * Processed data should be written in the sequence data directly. */
typedef int (*sequence_proc)(sequence *seq, int seq_layer, int frame_no, fits *fit, rectangle *source_area, void *arg);

#if 0
int sequence_processing(sequence *seq, sequence_proc process, int layer,
		gboolean run_in_thread, gboolean run_in_parallel,
		gboolean do_photometry, void *arg);
int	seqprocess_fwhm(sequence *seq, int seq_layer, int frame_no, fits *fit, rectangle *source_area, void *arg);
int	do_fwhm_sequence_processing(sequence *seq, int layer, gboolean print_psf, gboolean follow_star, gboolean run_in_thread, gboolean for_registration);
#endif
int	sequence_find_refimage(sequence *seq);
void	check_or_allocate_regparam(sequence *seq, int layer);
void	set_shifts(sequence *seq, int frame, int layer, float shiftx, float shifty, gboolean data_is_top_down);
sequence *create_internal_sequence(int size);
void	internal_sequence_set(sequence *seq, int index, fits *fit);
int	internal_sequence_find_index(sequence *seq, fits *fit);
fits	*internal_sequence_get(sequence *seq, int index);
gpointer crop_sequence(gpointer p);
gboolean sequence_is_rgb(sequence *seq);
void	enforce_area_in_image(rectangle *area, sequence *seq);

int seqpsf(sequence *seq, int layer, gboolean for_registration, gboolean regall,
		framing_mode framing, gboolean run_in_thread);

/* in export.c now */
void	update_export_crop_label();

#endif
