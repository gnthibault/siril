#ifndef STACKING_H_
#define STACKING_H_

#include "core/processing.h"

//#define STACK_DEBUG

/* the stacking method */
struct stacking_args;
typedef int (*stack_method)(struct stacking_args *args);

typedef struct normalization_coeff norm_coeff;

typedef enum {
	STACK_SUM,
	STACK_MEAN,
	STACK_MEDIAN,
	STACK_MAX,
	STACK_MIN,
} stackMethod;

/* TYPE OF SIGMA CLIPPING */
typedef enum {
	NO_REJEC,
	PERCENTILE,
	SIGMA,
	SIGMEDIAN,
	WINSORIZED,
	LINEARFIT,
} rejection;

/* TYPE OF NORMALIZATION */
typedef enum {
	NO_NORM,
	ADDITIVE,
	MULTIPLICATIVE,
	ADDITIVE_SCALING,
	MULTIPLICATIVE_SCALING,
} normalization;

/* identical to the combo box items */
typedef enum {
	ALL_IMAGES,
	SELECTED_IMAGES,
	BEST_PSF_IMAGES,
	BEST_ROUND_IMAGES,
	BEST_QUALITY_IMAGES
} stackType;

struct normalization_coeff {
	double *offset;
	double *mul;
	double *scale;
};

struct stacking_args {
	stack_method method;
	sequence *seq;
	int ref_image;	// takes precedences over seq->reference_image which may not be applicable
	seq_image_filter filtering_criterion;
	double filtering_parameter;
	int nb_images_to_stack; // calculated from the above, for display purposes
	int *image_indices;	// conversion between selected image indices and sequence image indices
	char *description;	// description of the filtering
	const char *output_filename;	// used in the idle function only
	gboolean output_overwrite;	// used in the idle function only
	struct timeval t_start;
	int retval;
	int max_number_of_rows;	/* number of rows that can be processed simultaneously,
				   function of max memory, image size and nb_images_to_stack */
	double sig[2];		/* low and high sigma rejection */
	rejection type_of_rejection;	/* type of rejection */
	normalization normalize;	/* type of normalization */
	norm_coeff coeff;		/* normalization data */
	gboolean force_norm;		/* TRUE = force normalization */
	gboolean norm_to_16;		/* normalize final image to 16bits */
	int reglayer;		/* layer used for registration data */
};

/* configuration from the command line */
struct stacking_configuration {
	struct timeval t_start;
	gchar *seqfile;
	gchar *result_file;
	stack_method method;
	double sig[2];
	gboolean force_no_norm;
	normalization norm;
	int number_of_loaded_sequences;
	float f_fwhm, f_fwhm_p, f_round, f_round_p, f_quality, f_quality_p; // on if >0
	gboolean filter_included;
};

void initialize_stacking_methods();

int stack_get_max_number_of_rows(sequence *seq, int nb_images_to_stack);

int stack_median(struct stacking_args *args);
int stack_mean_with_rejection(struct stacking_args *args);
int stack_addmax(struct stacking_args *args);
int stack_addmin(struct stacking_args *args);

void main_stack(struct stacking_args *args);
void clean_end_stacking(struct stacking_args *args);

void get_sequence_filtering_from_gui(seq_image_filter *filtering_criterion,
		double *filtering_parameter);
void update_stack_interface(gboolean dont_change_stack_type);


	/* normalization functions, normalize.c */

int do_normalization(struct stacking_args *args);


	/* median and mean functions */

struct _image_block {
	unsigned long channel, start_row, end_row, height;
};

/* pool of memory blocks for parallel processing */
struct _data_block {
	WORD **pix;	// buffer for a block on all images
	WORD *tmp;	// the actual single buffer for pix
	WORD *stack;	// the reordered stack for one pixel in all images
	int *rejected;  // 0 if pixel ok, 1 or -1 if rejected
	WORD *w_stack;	// stack for the winsorized rejection
	double *xf, *yf;// data for the linear fit rejection
};

int stack_open_all_files(struct stacking_args *args, int *bitpix, int *naxis, long *naxes, double *exposure, fits *fit);
int stack_create_result_fit(fits *fit, int bitpix, int naxis, long *naxes);
int stack_compute_parallel_blocks(struct _image_block **blocks, int max_number_of_rows,
		int nb_channels, long *naxes, long *largest_block_height,
		int *nb_parallel_stacks);
void stack_read_block_data(struct stacking_args *args, int use_regdata,
		struct _image_block *my_block, struct _data_block *data, long *naxes);
int find_refimage_in_indices(int *indices, int nb, int ref);

	/* up-scaling functions */

int upscale_sequence(struct stacking_args *args);
void remove_tmp_drizzle_files(struct stacking_args *args);

#endif
