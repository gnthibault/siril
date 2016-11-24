#ifndef _REGISTRATION_H_
#define _REGISTRATION_H_

#include "core/siril.h"

#define NUMBER_OF_METHOD 5
#define SUPER_SAMPLING 2

struct registration_args;
typedef int (*registration_function)(struct registration_args *);

/* arguments passed to registration functions */
struct registration_args {
	registration_function func;	// the registration function
	sequence *seq;			// the sequence to register
	gboolean process_all_frames;	// all frames of the sequence (opposite of selected frames)
	rectangle selection;		// the selection rectangle
	int layer;			// layer of images on which the registration is computed
	struct timeval t_start;		// start time of func
	int retval;			// retval of func
	gboolean run_in_thread;		// true if the registration was run in a thread
	const gchar *prefix;		// prefix of the created sequence if any
	gboolean follow_star;		// follow star position between frames
	gboolean load_new_sequence;	// load a new sequence if success
	gboolean matchSelection;	// Match stars found in the seleciton of reference image
	opencv_interpolation interpolation; // type of rotation interpolation
	gboolean translation_only;	// don't rotate images

	/* data for generated sequence, for star alignment registration */
	int new_total;
	imgdata *imgparam;
	regdata *regparam;
};

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
void update_reg_interface(gboolean dont_change_reg_radio);
void compute_fitting_selection(rectangle *area, int hsteps, int vsteps, int preserve_square);
void get_the_registration_area(struct registration_args *reg_args,
		struct registration_method *method); // for compositing
void fill_comboboxregmethod();

/** getter */
int get_registration_layer();

/* mouse behaviour */
typedef enum {
	MOUSE_ACTION_NONE,
	MOUSE_ACTION_SELECT_REG_AREA,
	MOUSE_ACTION_SELECT_PREVIEW1,
	MOUSE_ACTION_SELECT_PREVIEW2,
	MOUSE_ACTION_DRAW_SAMPLES,
} mouse_status_enum;

extern mouse_status_enum mouse_status;	// defined in registration_preview.c

#endif
