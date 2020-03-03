/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <float.h>
#include <sys/stat.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/OS_utils.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/PSF_list.h"
#include "gui/sequence_list.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/ser.h"
#include "registration/registration.h"
#include "algos/noise.h"
#include "algos/sorting.h"
#include "stacking/sum.h"
#include "opencv/opencv.h"

#include "stacking.h"

static struct stacking_args stackparam = {	// parameters passed to stacking
	NULL, NULL, -1, NULL, -1.0, 0, NULL, NULL, NULL, FALSE, { 0, 0 }, -1, 0,
	{ 0, 0 }, NO_REJEC, NO_NORM, { NULL, NULL, NULL}, FALSE, -1
};

#define MAX_FILTERS 5
static struct filtering_tuple stackfilters[MAX_FILTERS];

stack_method stacking_methods[] = {
	stack_summing_generic, stack_mean_with_rejection, stack_median, stack_addmax, stack_addmin
};

static gboolean end_stacking(gpointer p);
static int stack_addminmax(struct stacking_args *args, gboolean ismax);
static void stacking_args_deep_copy(struct stacking_args *from, struct stacking_args *to);
static void stacking_args_deep_free(struct stacking_args *args);


void initialize_stacking_methods() {
	GtkComboBoxText *stackcombo, *rejectioncombo;

	stackcombo = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(builder, "comboboxstack_methods"));
	rejectioncombo = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(builder, "comborejection"));
	gtk_combo_box_set_active(GTK_COMBO_BOX(stackcombo), com.stack.method);
	gtk_combo_box_set_active(GTK_COMBO_BOX(rejectioncombo), com.stack.rej_method);
}

gboolean evaluate_stacking_should_output_32bits(stack_method method, sequence *seq, int nb_img_to_stack) {
	if (method == stack_summing_generic) {
		if (seq->bitpix == BYTE_IMG)
			return nb_img_to_stack > 256;
		return TRUE;
	}
	if (method == stack_mean_with_rejection) {
		return TRUE;
	}
	if (method == stack_median) {
		return TRUE;
	}
	return seq->bitpix == FLOAT_IMG; // for min or max, only use it if input is already float
}

/******************************* ADDMIN AND ADDMAX STACKING ******************************
 * These methods are very close to summing stacking instead that the result
 * takes only the pixel if it is brighter (max) or dimmer (min) than the
 * previous one at the same coordinates.
 */
int stack_addmax(struct stacking_args *args) {
	return stack_addminmax(args, TRUE);
}

int stack_addmin(struct stacking_args *args) {
	return stack_addminmax(args, FALSE);
}

static int stack_addminmax(struct stacking_args *args, gboolean ismax) {
	int x, y, nx, ny, i, ii, j, shiftx, shifty, layer, reglayer;
	WORD *final_pixel[3];
	float *ffinal_pixel[3];
	double exposure=0.0;
	gboolean is_float;
	long nbdata = 0;
	char *tmpmsg, filename[256];
	int retval = 0, nb_frames, cur_nb = 0;
	fits fit = { 0 };

	/* should be pre-computed to display it in the stacking tab */
	nb_frames = args->nb_images_to_stack;
	reglayer = get_registration_layer(args->seq);

	if (nb_frames <= 1) {
		siril_log_message(_("No frame selected for stacking (select at least 2). Aborting.\n"));
		return -1;
	}

	final_pixel[0] = NULL;
	ffinal_pixel[0] = NULL;
	g_assert(args->seq->nb_layers == 1 || args->seq->nb_layers == 3);
	g_assert(nb_frames <= args->seq->number);

	for (j=0; j<args->seq->number; ++j){
		if (!get_thread_run()) {
			retval = -1;
			goto free_and_reset_progress_bar;
		}
		if (!args->filtering_criterion(args->seq, j, args->filtering_parameter)) {
			fprintf(stdout, "image %d is excluded from stacking\n", j);
			continue;
		}
		if (!seq_get_image_filename(args->seq, j, filename)) {
			retval = -1;
			goto free_and_reset_progress_bar;
		}
		tmpmsg = strdup(_("Processing image "));
		tmpmsg = str_append(&tmpmsg, filename);
		set_progress_bar_data(tmpmsg, (double)cur_nb/((double)nb_frames+1.));
		free(tmpmsg);

		cur_nb++;	// only used for progress bar

		if (seq_read_frame(args->seq, j, &fit)) {
			siril_log_message(_("Stacking: could not read frame, aborting\n"));
			retval = -3;
			goto free_and_reset_progress_bar;
		}

		g_assert(args->seq->nb_layers == 1 || args->seq->nb_layers == 3);
		g_assert(fit.naxes[2] == args->seq->nb_layers);

		/* first loaded image: init data structures for stacking */
		if (!nbdata) {
			is_float = fit.type == DATA_FLOAT;
			nbdata = fit.naxes[0] * fit.naxes[1];
			if (is_float) {
				if (ismax)
					ffinal_pixel[0] = calloc(nbdata * fit.naxes[2], sizeof(float));
				else {
					ffinal_pixel[0] = malloc(nbdata * fit.naxes[2] * sizeof(float));
					for (long k = 0; k < nbdata * fit.naxes[2]; k++)
						ffinal_pixel[0][k] = 1.0;
				}
				if (!ffinal_pixel[0]) {
					PRINT_ALLOC_ERR;
					retval = -2;
					goto free_and_reset_progress_bar;
				}
				if(args->seq->nb_layers == 3){
					ffinal_pixel[1] = ffinal_pixel[0] + nbdata;
					ffinal_pixel[2] = ffinal_pixel[1] + nbdata;
				}
			} else {
				if (ismax)
					final_pixel[0] = calloc(nbdata * fit.naxes[2], sizeof(WORD));
				else {
					final_pixel[0] = malloc(nbdata * fit.naxes[2] * sizeof(WORD));
					for (long k = 0; k < nbdata * fit.naxes[2]; k++)
						final_pixel[0][k] = USHRT_MAX;
				}
				if (!final_pixel[0]) {
					PRINT_ALLOC_ERR;
					retval = -2;
					goto free_and_reset_progress_bar;
				}
				if(args->seq->nb_layers == 3){
					final_pixel[1] = final_pixel[0] + nbdata;
					final_pixel[2] = final_pixel[1] + nbdata;
				}
			}
		} else if (fit.ry * fit.rx != nbdata) {
			siril_log_message(_("Stacking: image in sequence doesn't has the same dimensions\n"));
			retval = -3;
			goto free_and_reset_progress_bar;
		}

		/* load registration data for current image */
		if(reglayer != -1 && args->seq->regparam[reglayer]) {
			shiftx = round_to_int(args->seq->regparam[reglayer][j].shiftx * (float)args->seq->upscale_at_stacking);
			shifty = round_to_int(args->seq->regparam[reglayer][j].shifty * (float)args->seq->upscale_at_stacking);
		} else {
			shiftx = 0;
			shifty = 0;
		}
#ifdef STACK_DEBUG
		printf("stack image %d with shift x=%d y=%d\n", j, shiftx, shifty);
#endif

		/* Summing the exposure */
		exposure += fit.exposure;

		/* stack current image */
		i=0;	// index in final_pixel[0]
		for (y=0; y < fit.ry; ++y){
			for (x=0; x < fit.rx; ++x){
				nx = x - shiftx;
				ny = y - shifty;
				//printf("x=%d y=%d sx=%d sy=%d i=%d ii=%d\n",x,y,shiftx,shifty,i,ii);
				if (nx >= 0 && nx < fit.rx && ny >= 0 && ny < fit.ry) {
					ii = ny * fit.rx + nx;		// index in final_pixel[0] too
					//printf("shiftx=%d shifty=%d i=%d ii=%d\n",shiftx,shifty,i,ii);
					if (ii > 0 && ii < fit.rx * fit.ry){
						for(layer=0; layer<args->seq->nb_layers; ++layer){
							if (is_float) {
								float current_pixel = fit.fpdata[layer][ii];
								// we take the brightest pixel
								if ((ismax && current_pixel > ffinal_pixel[layer][i]) ||	
										// we take the darkest pixel
										(!ismax && current_pixel < ffinal_pixel[layer][i]))
									ffinal_pixel[layer][i] = current_pixel;
							} else {
								WORD current_pixel = fit.pdata[layer][ii];
								// we take the brightest pixel
								if ((ismax && current_pixel > final_pixel[layer][i]) ||	
										// we take the darkest pixel
										(!ismax && current_pixel < final_pixel[layer][i]))
									final_pixel[layer][i] = current_pixel;
							}
						}
					}
				}
				++i;
			}
		}
	}
	if (!get_thread_run()) {
		retval = -1;
		goto free_and_reset_progress_bar;
	}
	set_progress_bar_data(_("Finalizing stacking..."), (double)nb_frames/((double)nb_frames+1.));

	clearfits(&gfit);
	copyfits(&fit, &gfit, CP_FORMAT, 0);
	if (is_float) {
		gfit.fdata = ffinal_pixel[0];
		gfit.fpdata[RLAYER] = gfit.fdata;
		if (fit.naxes[2] == 3) {
			gfit.fpdata[GLAYER] = gfit.fdata + nbdata;
			gfit.fpdata[BLAYER] = gfit.fdata + 2 * nbdata;
		} else {
			gfit.fpdata[GLAYER] = gfit.fdata;
			gfit.fpdata[BLAYER] = gfit.fdata;
		}
	} else {
		gfit.data = final_pixel[0];
		gfit.pdata[RLAYER] = gfit.data;
		if (fit.naxes[2] == 3) {
			gfit.pdata[GLAYER] = gfit.data + nbdata;
			gfit.pdata[BLAYER] = gfit.data + 2 * nbdata;
		} else {
			gfit.pdata[GLAYER] = gfit.data;
			gfit.pdata[BLAYER] = gfit.data;
		}
	}

free_and_reset_progress_bar:
	if (retval) {
		set_progress_bar_data(_("Stacking failed. Check the log."), PROGRESS_RESET);
		siril_log_message(_("Stacking failed.\n"));
	} else {
		set_progress_bar_data(_("Stacking complete."), PROGRESS_DONE);
	}
	
	return retval;
}

/* the function that prepares the stacking and runs it */
void main_stack(struct stacking_args *args) {
	int nb_allowed_files;
	g_assert(args->ref_image >= 0 && args->ref_image < args->seq->number);

	/* first of all we need to check if we can process the files */
	if (args->seq->type == SEQ_REGULAR) {
		if (!allow_to_open_files(args->nb_images_to_stack, &nb_allowed_files)) {
			siril_log_message(_("Your system does not allow one to open more than %d files at the same time. "
						"You may consider either to enhance this limit (the method depends of "
						"your Operating System) or to convert your FITS sequence into a SER "
						"sequence before stacking, or to stack with the \"sum\" method.\n"),
					nb_allowed_files);
			args->retval = -1;
			return;
		}
	}

	siril_log_message(args->description);
	if (args->use_32bit_output)
		siril_log_message(_("Stacking result will be stored as a 32-bit image\n"));

	// 1. normalization
	if (do_normalization(args)) // does nothing if NO_NORM
		return;
	// 2. up-scale
	if (upscale_sequence(args)) // does nothing if args->seq->upscale_at_stacking <= 1.05
		return;
	// 3. stack
	args->max_number_of_rows = stack_get_max_number_of_rows(args->seq, args->nb_images_to_stack);
	args->retval = args->method(args);
}

/* the function that runs the thread. */
gpointer stack_function_handler(gpointer p) {
	struct stacking_args *args = (struct stacking_args *)p;

	main_stack(args);

	// 4. save result and clean-up
	siril_add_idle(end_stacking, args);
	return GINT_TO_POINTER(args->retval);
}


/* starts a summing operation using data stored in the stackparam structure
 * function is not reentrant but can be called again after it has returned and the thread is running */
static void start_stacking() {
	static GtkComboBox *method_combo = NULL, *rejec_combo = NULL, *norm_combo = NULL;
	static GtkEntry *output_file = NULL;
	static GtkToggleButton *overwrite = NULL, *force_norm = NULL;
	static GtkSpinButton *sigSpin[2] = {NULL, NULL};
	static GtkWidget *norm_to_max = NULL;

	if (method_combo == NULL) {
		method_combo = GTK_COMBO_BOX(gtk_builder_get_object(builder, "comboboxstack_methods"));
		output_file = GTK_ENTRY(gtk_builder_get_object(builder, "entryresultfile"));
		overwrite = GTK_TOGGLE_BUTTON(gtk_builder_get_object(builder, "checkbutoverwrite"));
		sigSpin[0] = GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button"));
		sigSpin[1] = GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button"));
		rejec_combo = GTK_COMBO_BOX(lookup_widget("comborejection"));
		norm_combo = GTK_COMBO_BOX(lookup_widget("combonormalize"));
		force_norm = GTK_TOGGLE_BUTTON(lookup_widget("checkforcenorm"));
		norm_to_max = lookup_widget("check_normalise_to_max");
	}

	if (get_thread_run()) {
		siril_log_message(_("Another task is already in progress, ignoring new request.\n"));
		return;
	}

	stackparam.sig[0] = (float) gtk_spin_button_get_value(sigSpin[0]);
	stackparam.sig[1] = (float) gtk_spin_button_get_value(sigSpin[1]);
	stackparam.type_of_rejection = gtk_combo_box_get_active(rejec_combo);
	stackparam.normalize = gtk_combo_box_get_active(norm_combo);
	stackparam.force_norm = gtk_toggle_button_get_active(force_norm);
	stackparam.norm_to_max= gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(norm_to_max)) && gtk_widget_is_visible(norm_to_max);
	stackparam.coeff.offset = NULL;
	stackparam.coeff.mul = NULL;
	stackparam.coeff.scale = NULL;
	stackparam.method =
			stacking_methods[gtk_combo_box_get_active(method_combo)];

	stackparam.use_32bit_output = evaluate_stacking_should_output_32bits(stackparam.method,
			&com.seq, stackparam.nb_images_to_stack);

	// ensure we have no normalization if not supported by the stacking method
	if (stackparam.method != stack_median && stackparam.method != stack_mean_with_rejection)
		stackparam.normalize = NO_NORM;
	stackparam.seq = &com.seq;
	stackparam.reglayer = get_registration_layer(&com.seq);
	siril_log_color_message(_("Stacking will use registration data of layer %d if some exist.\n"), "salmon", stackparam.reglayer);

	/* Do not display that cause it uses the generic function that already
	 * displays this text
	 */
	if (stackparam.method != &stack_summing_generic)
		siril_log_color_message(_("Stacking: processing...\n"), "red");
	gettimeofday(&stackparam.t_start, NULL);
	set_cursor_waiting(TRUE);

	stackparam.output_overwrite = gtk_toggle_button_get_active(overwrite);
	stackparam.output_filename = gtk_entry_get_text(output_file);

	/* Stacking. Result is in gfit if success */
	struct stacking_args *params = malloc(sizeof(struct stacking_args));
	stacking_args_deep_copy(&stackparam, params);
	start_in_new_thread(stack_function_handler, params);
}

static void _show_summary(struct stacking_args *args) {
	const char *norm_str, *rej_str;

	siril_log_message(_("Integration of %d images:\n"), args->nb_images_to_stack);

	/* Type of algorithm */
	if (args->method == &stack_mean_with_rejection) {
		siril_log_message(_("Pixel combination ......... average\n"));
	} else if (args->method == &stack_summing_generic) {
		siril_log_message(_("Pixel combination ......... normalized sum\n"));
	} else if (args->method == &stack_median) {
		siril_log_message(_("Pixel combination ......... median\n"));
	} else if (args->method == &stack_addmin) {
		siril_log_message(_("Pixel combination ......... minimum\n"));
	} else if (args->method == &stack_addmax) {
		siril_log_message(_("Pixel combination ......... maximum\n"));
	} else {
		siril_log_message(_("Pixel combination ......... none\n"));
	}

	/* Normalisation */
	if (args->method != &stack_mean_with_rejection &&
			args->method != &stack_median ) {
		norm_str = _("none");
	} else {
		switch (args->normalize) {
		default:
		case NO_NORM:
			norm_str = _("none");
			break;
		case ADDITIVE:
			norm_str = _("additive");
			break;
		case MULTIPLICATIVE:
			norm_str = _("multiplicative");
			break;
		case ADDITIVE_SCALING:
			norm_str = _("additive + scaling");
			break;
		case MULTIPLICATIVE_SCALING:
			norm_str = _("multiplicative + scaling");
			break;
		}
	}

	siril_log_message(_("Normalization ............. %s\n"), norm_str);

	/* Type of rejection */
	if (args->method != &stack_mean_with_rejection) {
		siril_log_message(_("Pixel rejection ........... none\n"));
		siril_log_message(_("Rejection parameters ...... none\n"));
	}
	else {

		switch (args->type_of_rejection) {
		default:
		case NO_REJEC:
			rej_str = _("none");
			break;
		case PERCENTILE:
			rej_str = _("percentile clipping");
			break;
		case SIGMA:
			rej_str = _("sigma clipping");
			break;
		case SIGMEDIAN:
			rej_str = _("median sigma clipping");
			break;
		case WINSORIZED:
			rej_str = _("Winsorized sigma clipping");
			break;
		case LINEARFIT:
			rej_str = _("linear fit clipping");
			break;
		}
		siril_log_message(_("Pixel rejection ........... %s\n"), rej_str);
		siril_log_message(_("Rejection parameters ...... low=%.3f high=%.3f\n"),
				args->sig[0], args->sig[1]);
	}
}

static void _show_bgnoise(gpointer p) {
	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return;
	}
	set_cursor_waiting(TRUE);

	struct noise_data *args = malloc(sizeof(struct noise_data));
	args->fit = com.uniq->fit;
	args->verbose = FALSE;
	args->use_idle = TRUE;
	memset(args->bgnoise, 0.0, sizeof(double[3]));

	start_in_new_thread(noise, args);
}

void clean_end_stacking(struct stacking_args *args) {
	if (!args->retval)
		_show_summary(args);
	remove_tmp_drizzle_files(args);
}

/* because this idle function is called after one of many stacking method
 * functions, it contains all generic wrap-up stuff instead of only graphical
 * operations. */
static gboolean end_stacking(gpointer p) {
	struct timeval t_end;
	struct stacking_args *args = (struct stacking_args *)p;
	fprintf(stdout, "Ending stacking idle function, retval=%d\n", args->retval);
	stop_processing_thread();	// can it be done here in case there is no thread?

	if (!args->retval) {
		clear_stars_list();
		/* check in com.seq, because args->seq may have been replaced */
		if (com.seq.upscale_at_stacking > 1.05)
			com.seq.current = SCALED_IMAGE;
		else com.seq.current = RESULT_IMAGE;
		/* Warning: the previous com.uniq is not freed, but calling
		 * close_single_image() will close everything before reopening it,
		 * which is quite slow */
		com.uniq = calloc(1, sizeof(single));
		com.uniq->comment = strdup(_("Stacking result image"));
		com.uniq->nb_layers = gfit.naxes[2];
		com.uniq->layers = calloc(com.uniq->nb_layers, sizeof(layer_info));
		com.uniq->fit = &gfit;
		/* Giving summary if average rejection stacking */
		_show_summary(args);
		/* Giving noise estimation (new thread) */
		_show_bgnoise(com.uniq->fit);

		/* save stacking result */
		if (args->output_filename != NULL && args->output_filename[0] != '\0') {
			GStatBuf st;
			if (!g_stat(args->output_filename, &st)) {
				int failed = !args->output_overwrite;
				if (!failed) {
					if (g_unlink(args->output_filename) == -1)
						failed = 1;
					if (!failed && savefits(args->output_filename, &gfit))
						failed = 1;
					if (!failed)
						com.uniq->filename = strdup(args->output_filename);
				}
				if (failed)
					com.uniq->filename = strdup(_("Unsaved stacking result"));
			}
			else {
				if (!savefits(args->output_filename, &gfit))
					com.uniq->filename = strdup(args->output_filename);
				else com.uniq->filename = strdup(_("Unsaved stacking result"));
			}
			display_filename();
			set_precision_switch(); // set precision on screen
		}
		/* remove tmp files if exist (Drizzle) */
		remove_tmp_drizzle_files(args);

		waiting_for_thread();		// bgnoise
		adjust_cutoff_from_updated_gfit();	// computes min and max
		set_sliders_value_to_gfit();
		initialize_display_mode();

		sliders_mode_set_state(com.sliders);
		set_cutoff_sliders_max_values();

		set_display_mode();

		/* update menus */
		update_MenuItem();

		if (com.seq.current == SCALED_IMAGE)
			adjust_vport_size_to_image();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
		sequence_list_change_current();
		update_stack_interface(TRUE);
	}

	set_cursor_waiting(FALSE);
	/* Do not display time for stack_summing_generic
	 * cause it uses the generic function that already
	 * displays the time
	 */
	if (args->method != &stack_summing_generic) {
		gettimeofday(&t_end, NULL);
		show_time(args->t_start, t_end);
	}
	stacking_args_deep_free(args);
	return FALSE;
}

void on_seqstack_button_clicked (GtkButton *button, gpointer user_data){
	control_window_switch_to_tab(OUTPUT_LOGS);
	start_stacking();
}

void on_comboboxstack_methods_changed (GtkComboBox *box, gpointer user_data) {
	static GtkNotebook* notebook = NULL;
	if (!notebook)
		notebook = GTK_NOTEBOOK(gtk_builder_get_object(builder, "notebook4"));
	com.stack.method = gtk_combo_box_get_active(box);

	gtk_notebook_set_current_page(notebook, com.stack.method);
	update_stack_interface(TRUE);
	writeinitfile();
}

void on_combonormalize_changed (GtkComboBox *box, gpointer user_data) {
	GtkWidget *widgetnormalize = lookup_widget("combonormalize");
	GtkWidget *force_norm = lookup_widget("checkforcenorm");
	gtk_widget_set_sensitive(force_norm,
			gtk_combo_box_get_active(GTK_COMBO_BOX(widgetnormalize)) != 0);
}


void on_comborejection_changed (GtkComboBox *box, gpointer user_data) {
	rejection type_of_rejection = gtk_combo_box_get_active(box);
	GtkLabel *label_rejection[2] = {NULL, NULL};

	if (!label_rejection[0]) {
		label_rejection[0] = GTK_LABEL(lookup_widget("label_low"));
		label_rejection[1] = GTK_LABEL(lookup_widget("label_high"));
	}
	/* set default values */
	switch (type_of_rejection) {
		case NO_REJEC:
			gtk_widget_set_visible(lookup_widget("stack_siglow_button"), FALSE);
			gtk_widget_set_visible(lookup_widget("stack_sighigh_button"), FALSE);
			gtk_widget_set_visible(lookup_widget("label_low"), FALSE);
			gtk_widget_set_visible(lookup_widget("label_high"), FALSE);
			break;
		case PERCENTILE :
			gtk_widget_set_visible(lookup_widget("stack_siglow_button"), TRUE);
			gtk_widget_set_visible(lookup_widget("stack_sighigh_button"), TRUE);
			gtk_widget_set_visible(lookup_widget("label_low"), TRUE);
			gtk_widget_set_visible(lookup_widget("label_high"), TRUE);
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button")), 0.0, 1.0);
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button")), 0.0, 1.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button")), 0.2);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button")), 0.1);
			gtk_label_set_text(label_rejection[0], _("Percentile low: "));
			gtk_label_set_text(label_rejection[1], _("Percentile high: "));
			break;
		case LINEARFIT:
			gtk_widget_set_visible(lookup_widget("stack_siglow_button"), TRUE);
			gtk_widget_set_visible(lookup_widget("stack_sighigh_button"), TRUE);
			gtk_widget_set_visible(lookup_widget("label_low"), TRUE);
			gtk_widget_set_visible(lookup_widget("label_high"), TRUE);
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button")), 0.0, 10.0);
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button")), 0.0, 10.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button")), 5.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button")), 5.0);
			gtk_label_set_text(label_rejection[0], _("Linear low: "));
			gtk_label_set_text(label_rejection[1], _("Linear high: "));
			break;
		default:
		case SIGMA:
		case WINSORIZED:
			gtk_widget_set_visible(lookup_widget("stack_siglow_button"), TRUE);
			gtk_widget_set_visible(lookup_widget("stack_sighigh_button"), TRUE);
			gtk_widget_set_visible(lookup_widget("label_low"), TRUE);
			gtk_widget_set_visible(lookup_widget("label_high"), TRUE);
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button")), 0.0, 10.0);
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button")), 0.0, 10.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button")), 4.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button")), 3.0);
			gtk_label_set_text(label_rejection[0], _("Sigma low: "));
			gtk_label_set_text(label_rejection[1], _("Sigma high: "));
	}
	com.stack.rej_method = gtk_combo_box_get_active(box);
	writeinitfile();
}

int stack_get_max_number_of_rows(sequence *seq, int nb_images_to_stack) {
	int max_memory = get_max_memory_in_MB();
	if (max_memory > 0) {
		siril_log_message(_("Using %d MB memory maximum for stacking\n"), max_memory);
		int elem_size = get_data_type(seq->bitpix) == DATA_FLOAT ? sizeof(float) : sizeof(WORD);
		uint64_t number_of_rows = (uint64_t)max_memory * BYTES_IN_A_MB /
			((uint64_t)seq->rx * nb_images_to_stack * elem_size * com.max_thread);
		// this is how many rows we can load in parallel from all images of the
		// sequence and be under the limit defined in config in megabytes.
		// We want to avoid having blocks larger than the half or they will decrease parallelism
		if (number_of_rows > seq->ry)
			return seq->ry;
		if (number_of_rows * 2 > seq->ry)
			return seq->ry / 2;
		return truncate_to_int32(number_of_rows);
	} else {
		siril_log_message(_("Not using limits on maximum memory for stacking\n"));
		return (seq->ry / 4) + 1;
	}
}

int find_refimage_in_indices(int *indices, int nb, int ref) {
	int i;
	for (i = 0; i < nb; i++) {
		if (indices[i] == ref)
			return i;
	}
	return -1;
}

/****************************************************************/

void on_stacksel_changed(GtkComboBox *widget, gpointer user_data) {
	update_stack_interface(TRUE);
}

void on_spinbut_percent_change(GtkSpinButton *spinbutton, gpointer user_data) {
	update_stack_interface(TRUE);
}

void on_filter_add1_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(lookup_widget("combofilter2"), TRUE);
	gtk_widget_set_visible(lookup_widget("stackspin2"), TRUE);
	gtk_widget_set_visible(lookup_widget("filter_add2"), TRUE);
	gtk_widget_set_visible(lookup_widget("filter_rem2"), TRUE);
	gtk_widget_set_visible(lookup_widget("labelfilter2"), TRUE);
	update_stack_interface(TRUE);
}

void on_filter_add2_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(lookup_widget("combofilter3"), TRUE);
	gtk_widget_set_visible(lookup_widget("stackspin3"), TRUE);
	gtk_widget_set_visible(lookup_widget("filter_rem3"), TRUE);
	gtk_widget_set_visible(lookup_widget("labelfilter3"), TRUE);
	update_stack_interface(TRUE);
}

void on_filter_rem2_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(lookup_widget("combofilter2"), FALSE);
	gtk_widget_set_visible(lookup_widget("stackspin2"), FALSE);
	gtk_widget_set_visible(lookup_widget("filter_add2"), FALSE);
	gtk_widget_set_visible(lookup_widget("filter_rem2"), FALSE);
	gtk_widget_set_visible(lookup_widget("labelfilter2"), FALSE);
	update_stack_interface(TRUE);
}

void on_filter_rem3_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(lookup_widget("combofilter3"), FALSE);
	gtk_widget_set_visible(lookup_widget("stackspin3"), FALSE);
	gtk_widget_set_visible(lookup_widget("filter_rem3"), FALSE);
	gtk_widget_set_visible(lookup_widget("labelfilter3"), FALSE);
	update_stack_interface(TRUE);
}

void get_sequence_filtering_from_gui(seq_image_filter *filtering_criterion,
		double *filtering_parameter) {
	int filter, guifilter, channel = 0, type;
	double percent = 0.0;
	static GtkComboBox *filter_combo[] = {NULL, NULL, NULL};
	static GtkAdjustment *stackadj[] = {NULL, NULL, NULL};
	static GtkWidget *spin[] = {NULL, NULL, NULL};
	if (!spin[0]) {
		spin[0] = lookup_widget("stackspin1");
		spin[1] = lookup_widget("stackspin2");
		spin[2] = lookup_widget("stackspin3");
		stackadj[0] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[0]));
		stackadj[1] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[1]));
		stackadj[2] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[2]));
		filter_combo[0] = GTK_COMBO_BOX(lookup_widget("combofilter1"));
		filter_combo[1] = GTK_COMBO_BOX(lookup_widget("combofilter2"));
		filter_combo[2] = GTK_COMBO_BOX(lookup_widget("combofilter3"));
	}
	for (filter = 0, guifilter = 0; guifilter < 3; guifilter++) {
		if (!gtk_widget_get_visible(GTK_WIDGET(filter_combo[guifilter]))) {
			continue;
		}

		type = gtk_combo_box_get_active(filter_combo[guifilter]);
		if (type != ALL_IMAGES && type != SELECTED_IMAGES) {
			channel = get_registration_layer(&com.seq);
			percent = gtk_adjustment_get_value(stackadj[guifilter]);
		}

		switch (type) {
			default:
			case ALL_IMAGES:
				stackfilters[filter].filter = seq_filter_all;
				stackfilters[filter].param = 0.0;
				gtk_widget_set_visible(spin[guifilter], FALSE);
				break;
			case SELECTED_IMAGES:
				stackfilters[filter].filter = seq_filter_included;
				stackfilters[filter].param = 0.0;
				gtk_widget_set_visible(spin[guifilter], FALSE);
				break;
			case BEST_PSF_IMAGES:
				stackfilters[filter].filter = seq_filter_fwhm;
				stackfilters[filter].param = compute_highest_accepted_fwhm(
						stackparam.seq, channel, percent);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				break;
			case BEST_WPSF_IMAGES:
				stackfilters[filter].filter = seq_filter_weighted_fwhm;
				stackfilters[filter].param = compute_highest_accepted_weighted_fwhm(
						stackparam.seq, channel, percent);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				break;
			case BEST_ROUND_IMAGES:
				stackfilters[filter].filter = seq_filter_roundness;
				stackfilters[filter].param = compute_lowest_accepted_roundness(
						stackparam.seq, channel, percent);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				break;
			case BEST_QUALITY_IMAGES:
				stackfilters[filter].filter = seq_filter_quality;
				stackfilters[filter].param = compute_lowest_accepted_quality(
						stackparam.seq, channel, percent);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				break;
		}
		filter++;
	}
	stackfilters[filter].filter = NULL;

	if (filter == 1) {
		*filtering_criterion = stackfilters[0].filter;
		*filtering_parameter = stackfilters[0].param;
	} else {
		*filtering_criterion = create_multiple_filter_from_list(stackfilters);
		*filtering_parameter = 0.0;
	}
}

static void update_filter_label() {
	static GtkComboBox *filter_combo[] = {NULL, NULL, NULL};
	static GtkLabel *filter_label[] = {NULL, NULL, NULL};
	gchar *filter_str;
	double param;
	int filter, type;

	if (!filter_combo[0]) {
		filter_combo[0] = GTK_COMBO_BOX(lookup_widget("combofilter1"));
		filter_combo[1] = GTK_COMBO_BOX(lookup_widget("combofilter2"));
		filter_combo[2] = GTK_COMBO_BOX(lookup_widget("combofilter3"));
		filter_label[0] = GTK_LABEL(lookup_widget("labelfilter1"));
		filter_label[1] = GTK_LABEL(lookup_widget("labelfilter2"));
		filter_label[2] = GTK_LABEL(lookup_widget("labelfilter3"));
	}

	for (filter = 0; filter < 3; filter++) {
		if (!gtk_widget_get_visible(GTK_WIDGET(filter_combo[filter]))) {
			break;
		}

		type = gtk_combo_box_get_active(filter_combo[filter]);
		param = stackfilters[filter].param;
		if (param == DBL_MIN || param == DBL_MAX || param == 0.0) {
			if (type == ALL_IMAGES || type == SELECTED_IMAGES)
				filter_str = g_strdup("");
			else filter_str = g_strdup("N/A");
		} else {
			switch (type) {
			default:
			case ALL_IMAGES:
			case SELECTED_IMAGES:
				filter_str = g_strdup("");
				break;
			case BEST_PSF_IMAGES:
			case BEST_WPSF_IMAGES:
				filter_str = g_strdup_printf("< %.2lf", param);
				break;
			case BEST_ROUND_IMAGES:
			case BEST_QUALITY_IMAGES:
				filter_str = g_strdup_printf("> %.3lf", param);
				break;
			}
		}
		gtk_label_set_text(filter_label[filter], filter_str);
		g_free(filter_str);
	}
}

/* Activates or not the stack button if there are 2 or more selected images,
 * all data related to stacking is set in stackparam, except the method itself,
 * determined at stacking start.
 */
void update_stack_interface(gboolean dont_change_stack_type) {
	static GtkWidget *go_stack = NULL,
			 *widgetnormalize = NULL, *force_norm = NULL, *norm_to_max = NULL;
	static GtkComboBox *method_combo = NULL, *filter_combo = NULL;
	static GtkLabel *result_label = NULL;
	gchar *labelbuffer;

	if(!go_stack) {
		go_stack = lookup_widget("gostack_button");
		filter_combo = GTK_COMBO_BOX(lookup_widget("combofilter1"));
		method_combo = GTK_COMBO_BOX(lookup_widget("comboboxstack_methods"));
		widgetnormalize = lookup_widget("combonormalize");
		force_norm = lookup_widget("checkforcenorm");
		norm_to_max = lookup_widget("check_normalise_to_max");
		result_label = GTK_LABEL(lookup_widget("stackfilter_label"));
	}
	if (!sequence_is_loaded()) {
		gtk_widget_set_sensitive(go_stack, FALSE);
		return;
	}
	stackparam.seq = &com.seq;

	if (!dont_change_stack_type && stackparam.seq->selnum < stackparam.seq->number) {
		g_signal_handlers_block_by_func(filter_combo, on_stacksel_changed, NULL);
		gtk_combo_box_set_active(filter_combo, SELECTED_IMAGES);
		g_signal_handlers_unblock_by_func(filter_combo, on_stacksel_changed, NULL);
	}

	switch (gtk_combo_box_get_active(method_combo)) {
	default:
	case STACK_SUM:
	case STACK_MAX:
	case STACK_MIN:
		gtk_widget_set_sensitive(widgetnormalize, FALSE);
		gtk_widget_set_sensitive(force_norm, FALSE);
		break;
	case STACK_MEAN:
	case STACK_MEDIAN:
		gtk_widget_set_sensitive(widgetnormalize, TRUE);
		gtk_widget_set_sensitive(force_norm,
				gtk_combo_box_get_active(GTK_COMBO_BOX(widgetnormalize)) != 0);
		gtk_widget_set_visible(norm_to_max, TRUE);
	}

	if (com.seq.reference_image == -1)
		com.seq.reference_image = sequence_find_refimage(&com.seq);
	stackparam.ref_image = com.seq.reference_image;

	get_sequence_filtering_from_gui(
			&stackparam.filtering_criterion, &stackparam.filtering_parameter);

	if (stackparam.description)
		free(stackparam.description);
	stackparam.description = describe_filter(stackparam.seq,
			stackparam.filtering_criterion, stackparam.filtering_parameter);

	update_filter_label();

	stackparam.nb_images_to_stack = compute_nb_filtered_images(&com.seq,
			stackparam.filtering_criterion, stackparam.filtering_parameter);
	labelbuffer = g_strdup_printf(_("Stacking %d images of the %d of the sequence"),
			stackparam.nb_images_to_stack, com.seq.number);
	gtk_label_set_text(result_label, labelbuffer);
	g_free(labelbuffer);

	if (stackparam.nb_images_to_stack >= 2) {
		stack_fill_list_of_unfiltered_images(&stackparam);
		gtk_widget_set_sensitive(go_stack, TRUE);
	} else {
		gtk_widget_set_sensitive(go_stack, FALSE);
	}
}

static void stacking_args_deep_copy(struct stacking_args *from, struct stacking_args *to) {
	memcpy(to, from, sizeof(struct stacking_args));
	// sequence is not duplicated
	to->image_indices = malloc(from->nb_images_to_stack * sizeof(int));
	memcpy(to->image_indices, from->image_indices, from->nb_images_to_stack * sizeof(int));
	to->description = strdup(from->description);
	// output_filename is not duplicated, can be changed until the last minute
}

static void stacking_args_deep_free(struct stacking_args *args) {
	free(args->image_indices);
	free(args->description);
	free(args);
}
