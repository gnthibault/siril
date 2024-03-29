#ifndef _REGISTRATION_H_
#define _REGISTRATION_H_

#include "core/siril.h"
#include "algos/PSF.h"
#include "core/processing.h"

#define NUMBER_OF_METHODS 7

struct registration_args;
typedef int (*registration_function)(struct registration_args *);

typedef enum {
	REQUIRES_NO_SELECTION,	// selection is not used
	REQUIRES_ANY_SELECTION,		// selection can be of any size and shape
	REQUIRES_SQUARED_SELECTION	// selection needs to be square-shaped
} selection_type;

typedef enum {
	REGTYPE_DEEPSKY, REGTYPE_PLANETARY
} registration_type;

typedef enum {
	PLANETARY_FULLDISK, PLANETARY_SURFACE
} planetary_type;

typedef enum {
	REG_PAGE_GLOBAL,
	REG_PAGE_COMET,
	REG_PAGE_3_STARS,
	REG_PAGE_MISC
} reg_notebook_page;

typedef enum {
	SHIFT_TRANSFORMATION,
	AFFINE_TRANSFORMATION,
	HOMOGRAPHY_TRANSFORMATION,
	FULLAFFINE_TRANSFORMATION,
} transformation_type;

/* arguments passed to registration functions */
struct registration_args {
	registration_function func;	// the registration function
	sequence *seq;			// the sequence to register
	int reference_image;		// reference image index
	gboolean process_all_frames;	// all frames of the sequence (opposite of selected frames)
	int layer;			// layer of images on which the registration is computed
	int retval;			// retval of func
	gboolean run_in_thread;		// true if the registration was run in a thread
	gboolean follow_star;		// follow star position between frames
	gboolean matchSelection;	// Match stars found in the seleciton of reference image
	rectangle selection;		// the selection rectangle
	gboolean x2upscale;		// apply an x2 upscale for pseudo drizzle
	gboolean cumul;			// cumul reg data with previous one
	int min_pairs;			// Minimum number of star pairs for success
	transformation_type type;   // Use affine transform  or homography

	/* data for generated sequence, for star alignment registration */
	gboolean translation_only;	// don't rotate images => no new sequence
	int new_total;                  // remaining images after registration
	imgdata *imgparam;		// imgparam for the new sequence
	regdata *regparam;		// regparam for the new sequence
	const gchar *prefix;		// prefix of the created sequence if any
	gboolean load_new_sequence;	// load the new sequence if success
	const gchar *new_seq_name;
	opencv_interpolation interpolation; // type of rotation interpolation
};

/* used to register a registration method */
struct registration_method {
	const char *name;
	registration_function method_ptr;
	selection_type sel;
	registration_type type;
};

struct registration_method *new_reg_method(const char *name, registration_function f,
		selection_type s, registration_type t); // for compositing
void initialize_registration_methods();
struct registration_method * get_selected_registration_method();
int register_shift_dft(struct registration_args *args);
int register_shift_fwhm(struct registration_args *args);
int register_star_alignment(struct registration_args *args);
int register_ecc(struct registration_args *args);
int register_comet(struct registration_args *regargs);
int register_3stars(struct registration_args *regargs);

pointf get_velocity();
void update_reg_interface(gboolean dont_change_reg_radio);
void compute_fitting_selection(rectangle *area, int hsteps, int vsteps, int preserve_square);
void get_the_registration_area(struct registration_args *reg_args,
		struct registration_method *method); // for compositing
void fill_comboboxregmethod();
gpointer register_thread_func(gpointer p);

/** getter */
int get_registration_layer(sequence *seq);


/**** star alignment (global and 3-star) registration ****/

struct star_align_data {
	struct registration_args *regargs;
	regdata *current_regdata;
	psf_star **refstars;
	int fitted_stars;
	BYTE *success;
	point ref;
};

regdata *star_align_get_current_regdata(struct registration_args *regargs);
int star_align_prepare_results(struct generic_seq_args *args);
int star_align_finalize_hook(struct generic_seq_args *args);

#endif
