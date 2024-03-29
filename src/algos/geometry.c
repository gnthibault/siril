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

#include <string.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "algos/statistics.h"
#include "algos/siril_wcs.h"
#include "core/undo.h"
#include "core/processing.h"
#include "gui/utils.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/PSF_list.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/siril_preview.h"
#include "opencv/opencv.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"

#include "geometry.h"
#include "astrometry_solver.h"

/* this method rotates the image 180 degrees, useful after german mount flip.
 * fit->rx, fit->ry, fit->naxes[2] and fit->pdata[*] are required to be assigned correctly */
static void fits_rotate_pi_ushort(fits *fit) {
	int i, line, axis;
	WORD *line1, *line2, *src, *dst, swap;

	size_t line_size = fit->rx * sizeof(WORD);
	line1 = malloc(line_size);
	line2 = malloc(line_size);
	if (!line1 || !line2) {
		PRINT_ALLOC_ERR;
		return;
	}

	for (axis = 0; axis < fit->naxes[2]; axis++) {
		for (line = 0; line < fit->ry / 2; line++) {
			src = fit->pdata[axis] + line * fit->rx;
			dst = fit->pdata[axis] + (fit->ry - line - 1) * fit->rx;

			memcpy(line1, src, line_size);
			for (i = 0; i < fit->rx / 2; i++) {
				swap = line1[i];
				line1[i] = line1[fit->rx - i - 1];
				line1[fit->rx - i - 1] = swap;
			}
			memcpy(line2, dst, line_size);
			for (i = 0; i < fit->rx / 2; i++) {
				swap = line2[i];
				line2[i] = line2[fit->rx - i - 1];
				line2[fit->rx - i - 1] = swap;
			}
			memcpy(src, line2, line_size);
			memcpy(dst, line1, line_size);
		}
		if (fit->ry & 1) {
			/* swap the middle line */
			src = fit->pdata[axis] + line * fit->rx;
			for (i = 0; i < fit->rx / 2; i++) {
				swap = src[i];
				src[i] = src[fit->rx - i - 1];
				src[fit->rx - i - 1] = swap;
			}
		}
	}
	free(line1);
	free(line2);
}

static void fits_rotate_pi_float(fits *fit) {
	int i, line, axis;
	float *line1, *line2, *src, *dst, swap;

	size_t line_size = fit->rx * sizeof(float);
	line1 = malloc(line_size);
	line2 = malloc(line_size);
	if (!line1 || !line2) {
		PRINT_ALLOC_ERR;
		return;
	}

	for (axis = 0; axis < fit->naxes[2]; axis++) {
		for (line = 0; line < fit->ry / 2; line++) {
			src = fit->fpdata[axis] + line * fit->rx;
			dst = fit->fpdata[axis] + (fit->ry - line - 1) * fit->rx;

			memcpy(line1, src, line_size);
			for (i = 0; i < fit->rx / 2; i++) {
				swap = line1[i];
				line1[i] = line1[fit->rx - i - 1];
				line1[fit->rx - i - 1] = swap;
			}
			memcpy(line2, dst, line_size);
			for (i = 0; i < fit->rx / 2; i++) {
				swap = line2[i];
				line2[i] = line2[fit->rx - i - 1];
				line2[fit->rx - i - 1] = swap;
			}
			memcpy(src, line2, line_size);
			memcpy(dst, line1, line_size);
		}
		if (fit->ry & 1) {
			/* swap the middle line */
			src = fit->fpdata[axis] + line * fit->rx;
			for (i = 0; i < fit->rx / 2; i++) {
				swap = src[i];
				src[i] = src[fit->rx - i - 1];
				src[fit->rx - i - 1] = swap;
			}
		}
	}
	free(line1);
	free(line2);
}

static void fits_rotate_pi(fits *fit) {
	if (fit->type == DATA_USHORT) {
		fits_rotate_pi_ushort(fit);
	} else if (fit->type == DATA_FLOAT) {
		fits_rotate_pi_float(fit);
	}
}

void mirrorx_gui(fits *fit) {
		set_cursor_waiting(TRUE);
		undo_save_state(fit, _("Mirror X"));
		mirrorx(fit, TRUE);
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
		set_cursor_waiting(FALSE);
}

void mirrory_gui(fits *fit) {
		set_cursor_waiting(TRUE);
		undo_save_state(fit, _("Mirror Y"));
		mirrory(fit, TRUE);
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
		set_cursor_waiting(FALSE);
}

static void rotate_gui(fits *fit) {
	static GtkToggleButton *crop_rotation = NULL;
	double angle = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_rotation")));
	int interpolation = gtk_combo_box_get_active(
			GTK_COMBO_BOX(lookup_widget("combo_interpolation_rotation")));
	int cropped;

	if (crop_rotation == NULL) {
		crop_rotation = GTK_TOGGLE_BUTTON(
				lookup_widget("checkbutton_rotation_crop"));
	}
	cropped = gtk_toggle_button_get_active(crop_rotation);

	set_cursor_waiting(TRUE);
	undo_save_state(fit, _("Rotation (%.1lfdeg, cropped=%s)"), angle,
			cropped ? "TRUE" : "FALSE");
	verbose_rotate_image(fit, angle, interpolation, cropped);
	
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

/* These functions do not more than resize_gaussian and rotate_image
 * except for console outputs.
 * Indeed, siril_log_message seems not working in a cpp file */
int verbose_resize_gaussian(fits *image, int toX, int toY, int interpolation) {
	int retvalue;
	const char *str_inter;
	struct timeval t_start, t_end;

	switch (interpolation) {
	case OPENCV_NEAREST:
		str_inter = _("Nearest-Neighbor");
		break;
	default:
	case OPENCV_LINEAR:
		str_inter = _("Bilinear");
		break;
	case OPENCV_AREA:
		str_inter = _("Pixel Area Relation");
		break;
	case OPENCV_CUBIC:
		str_inter = _("Bicubic");
		break;
	case OPENCV_LANCZOS4:
		str_inter = _("Lanczos4");
		break;
	}

	siril_log_color_message(_("Resample (%s interpolation): processing...\n"),
			"green", str_inter);

	gettimeofday(&t_start, NULL);

	retvalue = cvResizeGaussian(image, toX, toY, interpolation);
	invalidate_WCS_keywords(image);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	return retvalue;
}

int verbose_rotate_image(fits *image, double angle, int interpolation,
		int cropped) {
	const char *str_inter;
	struct timeval t_start, t_end;

	switch (interpolation) {
	case -1:
		str_inter = _("No");
		break;
	case OPENCV_NEAREST:
		str_inter = _("Nearest-Neighbor");
		break;
	default:
	case OPENCV_LINEAR:
		str_inter = _("Bilinear");
		break;
	case OPENCV_AREA:
		str_inter = _("Pixel Area Relation");
		break;
	case OPENCV_CUBIC:
		str_inter = _("Bicubic");
		break;
	case OPENCV_LANCZOS4:
		str_inter = _("Lanczos4");
		break;
	}

	siril_log_color_message(
			_("Rotation (%s interpolation, angle=%g): processing...\n"), "green",
			str_inter, angle);
	gettimeofday(&t_start, NULL);

	point center = {image->rx / 2.0, image->ry / 2.0};

	cvRotateImage(image, center, angle, interpolation, cropped);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	if (image->wcslib) {
		rotate_astrometry_data(image, center, angle, cropped);
		load_WCS_from_memory(image);
	}

	return 0;
}

static void mirrorx_ushort(fits *fit, gboolean verbose) {
	int line, axis;
	WORD *swapline, *src, *dst;
	struct timeval t_start, t_end;

	if (verbose) {
		siril_log_color_message(_("Horizontal mirror: processing...\n"), "red");
		gettimeofday(&t_start, NULL);
	}

	size_t line_size = fit->rx * sizeof(WORD);
	swapline = malloc(line_size);
	if (!swapline) {
		PRINT_ALLOC_ERR;
		return;
	}

	for (axis = 0; axis < fit->naxes[2]; axis++) {
		for (line = 0; line < fit->ry / 2; line++) {
			src = fit->pdata[axis] + line * fit->rx;
			dst = fit->pdata[axis] + (fit->ry - line - 1) * fit->rx;

			memcpy(swapline, src, line_size);
			memcpy(src, dst, line_size);
			memcpy(dst, swapline, line_size);
		}
	}
	free(swapline);
	if (verbose) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}
}

static void mirrorx_float(fits *fit, gboolean verbose) {
	int line, axis;
	float *swapline, *src, *dst;
	struct timeval t_start, t_end;

	if (verbose) {
		siril_log_color_message(_("Horizontal mirror: processing...\n"), "green");
		gettimeofday(&t_start, NULL);
	}

	size_t line_size = fit->rx * sizeof(float);
	swapline = malloc(line_size);
	if (!swapline) {
		PRINT_ALLOC_ERR;
		return;
	}

	for (axis = 0; axis < fit->naxes[2]; axis++) {
		for (line = 0; line < fit->ry / 2; line++) {
			src = fit->fpdata[axis] + line * fit->rx;
			dst = fit->fpdata[axis] + (fit->ry - line - 1) * fit->rx;

			memcpy(swapline, src, line_size);
			memcpy(src, dst, line_size);
			memcpy(dst, swapline, line_size);
		}
	}
	free(swapline);
	if (verbose) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}
}

void mirrorx(fits *fit, gboolean verbose) {
	if (fit->type == DATA_USHORT) {
		mirrorx_ushort(fit, verbose);
	} else if (fit->type == DATA_FLOAT) {
		mirrorx_float(fit, verbose);
	}
	if (fit->wcslib) {
		flip_bottom_up_astrometry_data(fit);
		load_WCS_from_memory(fit);
	}
}

void mirrory(fits *fit, gboolean verbose) {
	struct timeval t_start, t_end;

	if (verbose) {
		siril_log_color_message(_("Vertical mirror: processing...\n"), "green");
		gettimeofday(&t_start, NULL);
	}

	fits_flip_top_to_bottom(fit);
	fits_rotate_pi(fit);

	if (verbose) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}

	if (fit->wcslib) {
		flip_left_right_astrometry_data(fit);
		load_WCS_from_memory(fit);
	}
}

/* inplace cropping of the image in fit
 * fit->data is not realloc, only fit->pdata points to a different area and
 * data is correctly written to this new area, which makes this function
 * quite dangerous to use when fit is used for something else afterwards.
 */
static void crop_ushort(fits *fit, rectangle *bounds) {
	int i, j, layer;
	int newnbdata;
	struct timeval t_start = { 0 }, t_end = { 0 };

	if (fit == &gfit) {
		siril_log_color_message(_("Crop: processing...\n"), "red");
		gettimeofday(&t_start, NULL);
	}

	newnbdata = bounds->w * bounds->h;
	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		WORD *from = fit->pdata[layer]
				+ (fit->ry - bounds->y - bounds->h) * fit->rx + bounds->x;
		fit->pdata[layer] = fit->data + layer * newnbdata;
		WORD *to = fit->pdata[layer];
		int stridefrom = fit->rx - bounds->w;

		for (i = 0; i < bounds->h; ++i) {
			for (j = 0; j < bounds->w; ++j) {
				*to++ = *from++;
			}
			from += stridefrom;
		}
	}
	fit->rx = fit->naxes[0] = bounds->w;
	fit->ry = fit->naxes[1] = bounds->h;

	if (fit == &gfit) {
		clear_stars_list();
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}
	invalidate_stats_from_fit(fit);
}

static void crop_float(fits *fit, rectangle *bounds) {
	int i, j, layer;
	int newnbdata;
	struct timeval t_start = { 0 }, t_end = { 0 };

	if (fit == &gfit) {
		siril_log_color_message(_("Crop: processing...\n"), "green");
		gettimeofday(&t_start, NULL);
	}

	newnbdata = bounds->w * bounds->h;
	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		float *from = fit->fpdata[layer]
				+ (fit->ry - bounds->y - bounds->h) * fit->rx + bounds->x;
		fit->fpdata[layer] = fit->fdata + layer * newnbdata;
		float *to = fit->fpdata[layer];
		int stridefrom = fit->rx - bounds->w;

		for (i = 0; i < bounds->h; ++i) {
			for (j = 0; j < bounds->w; ++j) {
				*to++ = *from++;
			}
			from += stridefrom;
		}
	}
	fit->rx = fit->naxes[0] = bounds->w;
	fit->ry = fit->naxes[1] = bounds->h;

	if (fit == &gfit) {
		clear_stars_list();
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}
	invalidate_stats_from_fit(fit);
}

int crop(fits *fit, rectangle *bounds) {
	point shift; //need to be computed before fit rx/ry are altered by crop
	shift.x = (double)(bounds->x);
	shift.y = fit->ry - (double)(bounds->h) - (double)(bounds->y) - 1; // for top-bottom flip
	if (fit->type == DATA_USHORT) {
		crop_ushort(fit, bounds);
	} else if (fit->type == DATA_FLOAT) {
		crop_float(fit, bounds);
	} else {
		return -1;
	}
	if (fit->wcslib) {
		crop_astrometry_data(fit, shift);
		load_WCS_from_memory(fit);
	}
	return 0;
}

/************************* CALLBACKS *************************************/

/**
 *  ROTATION
 */
void siril_rotate90() {
		set_cursor_waiting(TRUE);
		undo_save_state(&gfit, _("Rotation (90.0deg)"));
		verbose_rotate_image(&gfit, 90.0, -1, 0);	// fast rotation, no interpolation, no crop
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
		set_cursor_waiting(FALSE);
}

void siril_rotate270() {
		set_cursor_waiting(TRUE);
		undo_save_state(&gfit, _("Rotation (-90.0deg)"));
		verbose_rotate_image(&gfit, 270.0, -1, 0);// fast rotation, no interpolation, no crop
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
		set_cursor_waiting(FALSE);
}

void on_button_rotation_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("rotation_dialog");
}

void on_button_rotation_ok_clicked(GtkButton *button, gpointer user_data) {
	rotate_gui(&gfit);
}

/*************
 * RESAMPLE
 */

/* Resample */

void on_button_resample_ok_clicked(GtkButton *button, gpointer user_data) {
	if (confirm_delete_wcs_keywords(&gfit)) {
		double sample[2];
		sample[0] = gtk_spin_button_get_value(
				GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_X")));
		sample[1] = gtk_spin_button_get_value(
				GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_Y")));
		int interpolation = gtk_combo_box_get_active(
				GTK_COMBO_BOX(lookup_widget("combo_interpolation")));

		set_cursor_waiting(TRUE);
		int toX = round_to_int((sample[0] / 100.0) * gfit.rx);
		int toY = round_to_int((sample[1] / 100.0) * gfit.ry);
		undo_save_state(&gfit, _("Resample (%g - %g)"), sample[0] / 100.0,
				sample[1] / 100.0);
		verbose_resize_gaussian(&gfit, toX, toY, interpolation);
		
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
		set_cursor_waiting(FALSE);
	}
}

void on_button_resample_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("resample_dialog");
}

void on_spinbutton_resample_X_value_changed(GtkSpinButton *spinbutton,
		gpointer user_data) {
	GtkToggleButton *ratio = GTK_TOGGLE_BUTTON(
			lookup_widget("button_sample_ratio"));
	double xvalue = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_X")));

	if (gtk_toggle_button_get_active(ratio))
		gtk_spin_button_set_value(
				GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_Y")),
				xvalue);
}

void on_spinbutton_resample_Y_value_changed(GtkSpinButton *spinbutton,
		gpointer user_data) {
	GtkToggleButton *ratio = GTK_TOGGLE_BUTTON(
			lookup_widget("button_sample_ratio"));
	double yvalue = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_Y")));

	if (gtk_toggle_button_get_active(ratio))
		gtk_spin_button_set_value(
				GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_X")),
				yvalue);
}

void on_button_sample_ratio_toggled(GtkToggleButton *button, gpointer user_data) {
	double xvalue = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_X")));

	if (gtk_toggle_button_get_active(button))
		gtk_spin_button_set_value(
				GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_Y")),
				xvalue);
}

/**************
 * CROP
 */
void siril_crop() {
	undo_save_state(&gfit, _("Crop (x=%d, y=%d, w=%d, h=%d)"),
			com.selection.x, com.selection.y, com.selection.w,
			com.selection.h);
	if (is_preview_active()) {
		siril_message_dialog(GTK_MESSAGE_INFO, _("A live preview session is active"),
				_("It is impossible to crop the image when a filter with preview session is active. "
						"Please consider to close the filter dialog first."));
		return;
	}
	crop(&gfit, &com.selection);
	delete_selected_area();
	reset_display_offset();
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
}

gint64 crop_compute_size_hook(struct generic_seq_args *args, int nb_frames) {
	struct crop_sequence_data *c_args = (struct crop_sequence_data*) args->user;
	double ratio = (double)(c_args->area.h * c_args->area.w) / (double)(args->seq->rx * args->seq->ry);
	double fullseqsize = seq_compute_size(args->seq, nb_frames, args->output_type);
	return (gint64)(fullseqsize * ratio);
}

int crop_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_) {
	struct crop_sequence_data *c_args = (struct crop_sequence_data*) args->user;

	return crop(fit, &(c_args->area));
}

/* TODO: should we use the partial image? */
gpointer crop_sequence(struct crop_sequence_data *crop_sequence_data) {
	struct generic_seq_args *args = create_default_seqargs(crop_sequence_data->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = crop_sequence_data->seq->selnum;
	args->compute_size_hook = crop_compute_size_hook;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = seq_finalize_hook;
	args->image_hook = crop_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Crop Sequence");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->new_seq_prefix = crop_sequence_data->prefix;
	args->load_new_sequence = TRUE;
	args->user = crop_sequence_data;

	start_in_new_thread(generic_sequence_worker, args);

	return 0;
}

/*** GUI for crop sequence */
void on_crop_Apply_clicked(GtkButton *button, gpointer user_data) {
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

#ifdef HAVE_FFMS2
	if (com.seq.type == SEQ_AVI) {
		siril_log_message(_("Crop does not work with "
				"avi film. Please, convert your file to SER first.\n"));
		return;
	}
#endif
	if (com.seq.type == SEQ_INTERNAL) {
		siril_log_message(_("Not a valid sequence for cropping.\n"));
	}

	struct crop_sequence_data *args = malloc(sizeof(struct crop_sequence_data));

	GtkEntry *cropped_entry = GTK_ENTRY(lookup_widget("cropped_entry"));

	args->seq = &com.seq;
	memcpy(&args->area, &com.selection, sizeof(rectangle));
	args->prefix = gtk_entry_get_text(cropped_entry);

	set_cursor_waiting(TRUE);
	crop_sequence(args);
}

void on_crop_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("crop_dialog");
}
