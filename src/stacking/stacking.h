#ifndef STACKING_H_
#define STACKING_H_

#include "core/processing.h"

//#define STACK_DEBUG

/* the stacking method */
struct stacking_args;
typedef int (*stack_method)(struct stacking_args *args);

typedef struct normalization_coeff norm_coeff;

enum {
	ST_ALLOC_ERROR = -10,
	ST_SEQUENCE_ERROR = -2,
	ST_GENERIC_ERROR = -1,
	ST_OK = 0
};

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
	MAD,
	SIGMEDIAN,
	WINSORIZED,
	LINEARFIT,
	GESDT
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
	BEST_WPSF_IMAGES,
	BEST_ROUND_IMAGES,
	BEST_QUALITY_IMAGES
} stackType;

struct normalization_coeff {
	double *offset;
	double *mul;
	double *scale;
	double *poffset[3];
	double *pmul[3];
	double *pscale[3];
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
	float sig[2];		/* low and high sigma rejection or GESTD parameters */
	float *critical_value; /* index of critical_values for GESTD */
	rejection type_of_rejection;	/* type of rejection */
	normalization normalize;	/* type of normalization */
	norm_coeff coeff;		/* normalization data */
	gboolean force_norm;		/* TRUE = force normalization */
	gboolean output_norm;		/* normalize final image to the [0, 1] range */
	gboolean use_32bit_output;	/* output to 32 bit float */
	int reglayer;		/* layer used for registration data */

	gboolean apply_weight;			/* enable weights */
	double *weights; 				/* computed weights for each (layer,image)*/

	float (*sd_calculator)(const WORD *, const int); // internal, for ushort
	float (*mad_calculator)(const WORD *, const size_t, const double, gboolean) ; // internal, for ushort
};

/* configuration from the command line */
struct stacking_configuration {
	struct timeval t_start;
	gchar *seqfile;
	gchar *result_file;
	stack_method method;
	rejection type_of_rejection;	/* type of rejection */
	double sig[2];
	gboolean force_no_norm;
	gboolean output_norm;
	normalization norm;
	int number_of_loaded_sequences;
	float f_fwhm, f_fwhm_p, f_wfwhm, f_wfwhm_p, f_round, f_round_p, f_quality, f_quality_p; // on if >0
	gboolean filter_included;
	gboolean apply_weight;
};


/*
 * ESD test statistic data.
 */
struct outliers {
	float x;
	int i;
	int out;
};


void initialize_stacking_default();
void initialize_stacking_methods();
gboolean evaluate_stacking_should_output_32bits(stack_method method,
		sequence *seq, int nb_img_to_stack, gchar **err);

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

int check_G_values(float Gs, float Gc);
void confirm_outliers(struct outliers *out, int N, double median, int *rejected, guint64 rej[2]);

struct _image_block {
	long channel, start_row, end_row, height; // long matches naxes type
};

int stack_compute_parallel_blocks(struct _image_block **blocksptr, long max_number_of_rows,
		long naxes[3], int nb_threads, long *largest_block_height, int *nb_blocks);

/* pool of memory blocks for parallel processing */
struct _data_block {
	void **pix;	// buffer for a block on all images
	void *tmp;	// the actual single buffer for pix
	void *stack;	// the reordered stack for one pixel in all images
	int *rejected;  // 0 if pixel ok, 1 or -1 if rejected
	void *w_stack;	// stack for the winsorized rejection
	void *o_stack;  // original unordered stack
	float *xf, *yf, m_x, m_dx2;// data for the linear fit rejection
	int layer;	// to identify layer for normalization
};

int stack_open_all_files(struct stacking_args *args, int *bitpix, int *naxis, long *naxes, GList **date_time, fits *fit);
int find_refimage_in_indices(int *indices, int nb, int ref);

	/* up-scaling functions */

int upscale_sequence(struct stacking_args *args);
void remove_tmp_drizzle_files(struct stacking_args *args);


	/* rejection_float.c */

int apply_rejection_float(struct _data_block *data, int nb_frames, struct stacking_args *args, guint64 crej[2]);

#endif
