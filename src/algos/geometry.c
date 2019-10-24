/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2019 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"
#include "algos/plateSolver.h"
#include "algos/statistics.h"
#include "core/undo.h"
#include "core/processing.h"
#include "gui/callbacks.h"
#include "gui/PSF_list.h"
#include "gui/dialogs.h"
#include "opencv/opencv.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "geometry.h"

#include <string.h>
#include <math.h>

/* this method rotates the image 180 degrees, useful after german mount flip.
 * fit->rx, fit->ry, fit->naxes[2] and fit->pdata[*] are required to be assigned correctly */
static void fits_rotate_pi(fits *fit) {
	int i, line, axis, line_size;
	WORD *line1, *line2, *src, *dst, swap;

	line_size = fit->rx * sizeof(WORD);
	line1 = malloc(line_size);
	line2 = malloc(line_size);

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

static void mirrorx_gui(fits *fit) {
		if (confirm_delete_wcs_keywords(fit)) {
			set_cursor_waiting(TRUE);
			undo_save_state(&gfit, "Processing: Mirror X");
			mirrorx(fit, TRUE);
			redraw(com.cvport, REMAP_ALL);
			redraw_previews();
			set_cursor_waiting(FALSE);
		}
}

static void mirrory_gui(fits *fit) {
		if (confirm_delete_wcs_keywords(fit)) {
			set_cursor_waiting(TRUE);
			undo_save_state(&gfit, "Processing: Mirror Y");
			mirrory(fit, TRUE);
			redraw(com.cvport, REMAP_ALL);
			redraw_previews();
			set_cursor_waiting(FALSE);
		}
}

static void rotate_gui(fits *fit) {
	if (confirm_delete_wcs_keywords(fit)) {
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
		undo_save_state(&gfit, "Processing: Rotation (%.1lfdeg, cropped=%s)", angle,
				cropped ? "TRUE" : "FALSE");
		verbose_rotate_image(fit, angle, interpolation, cropped);
		update_used_memory();
		adjust_vport_size_to_image();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
		set_cursor_waiting(FALSE);
	}
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
			"red", str_inter);

	gettimeofday(&t_start, NULL);

	retvalue = cvResizeGaussian(&gfit, toX, toY, interpolation);
	invalidate_WCS_keywords(&gfit);

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
			_("Rotation (%s interpolation, angle=%g): processing...\n"), "red",
			str_inter, angle);
	gettimeofday(&t_start, NULL);

	point center = {gfit.rx / 2.0, gfit.ry / 2.0};
	cvRotateImage(&gfit, center, angle, interpolation, cropped);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	invalidate_WCS_keywords(&gfit);

	return 0;
}

void mirrorx(fits *fit, gboolean verbose) {
	int line, axis, line_size;
	WORD *swapline, *src, *dst;
	struct timeval t_start, t_end;

	if (verbose) {
		siril_log_color_message(_("Horizontal mirror: processing...\n"), "red");
		gettimeofday(&t_start, NULL);
	}

	line_size = fit->rx * sizeof(WORD);
	swapline = malloc(line_size);

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
	invalidate_WCS_keywords(fit);
}

void mirrory(fits *fit, gboolean verbose) {
	struct timeval t_start, t_end;

	if (verbose) {
		siril_log_color_message(_("Vertical mirror: processing...\n"), "red");
		gettimeofday(&t_start, NULL);
	}

	fits_flip_top_to_bottom(fit);
	fits_rotate_pi(fit);

	if (verbose) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}
	invalidate_WCS_keywords(fit);
}

/* inplace cropping of the image in fit
 * fit->data is not realloc, only fit->pdata points to a different area and
 * data is correctly written to this new area, which makes this function
 * quite dangerous to use when fit is used for something else afterwards.
 */
int crop(fits *fit, rectangle *bounds) {
	int i, j, layer;
	int newnbdata;
	struct timeval t_start, t_end;

	memset(&t_start, 0, sizeof(struct timeval));
	memset(&t_end, 0, sizeof(struct timeval));

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
	invalidate_WCS_keywords(fit);

	return 0;
}

/************************* CALLBACKS *************************************/

/**
 *  ROTATION
 */
void on_menuitem_rotation90_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (confirm_delete_wcs_keywords(&gfit)) {
		static GtkToggleButton *crop_rotation = NULL;
		int cropped;

		if (crop_rotation == NULL) {
			crop_rotation = GTK_TOGGLE_BUTTON(
					lookup_widget("checkbutton_rotation_crop"));
		}
		cropped = gtk_toggle_button_get_active(crop_rotation);

		set_cursor_waiting(TRUE);
		undo_save_state(&gfit, "Processing: Rotation (90.0deg)");
		verbose_rotate_image(&gfit, 90.0, -1, cropped);	// fast rotation, no interpolation, no crop
		adjust_vport_size_to_image();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
		set_cursor_waiting(FALSE);
	}
}

void on_menuitem_rotation270_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (confirm_delete_wcs_keywords(&gfit)) {
		static GtkToggleButton *crop_rotation = NULL;
		int cropped;

		if (crop_rotation == NULL) {
			crop_rotation = GTK_TOGGLE_BUTTON(
					lookup_widget("checkbutton_rotation_crop"));
		}
		cropped = gtk_toggle_button_get_active(crop_rotation);

		set_cursor_waiting(TRUE);
		undo_save_state(&gfit, "Processing: Rotation (-90.0deg)");
		verbose_rotate_image(&gfit, 270.0, -1, cropped);// fast rotation, no interpolation, no crop
		adjust_vport_size_to_image();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
		set_cursor_waiting(FALSE);
	}
}

void on_menuitem_rotation_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded())
		siril_open_dialog("rotation_dialog");
}

void on_button_rotation_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("rotation_dialog");
}

void on_spinbutton_rotation_value_changed(GtkSpinButton *spin_button,
		gpointer user_data) {
	double angle = gtk_spin_button_get_value(spin_button);

	angle = fmod(angle, 360);
	while (angle < 0) {
		angle += 360.0;
	}
	gtk_spin_button_set_value(spin_button, angle);
}

void on_button_rotation_ok_clicked(GtkButton *button, gpointer user_data) {
	rotate_gui(&gfit);
}

/******
 * MIRROR
 */

void on_menuitem_mirrorx_activate(GtkMenuItem *menuitem, gpointer user_data) {
	mirrorx_gui(&gfit);
}

void on_mirrorx_button_clicked(GtkToolButton *button, gpointer user_data) {
	mirrorx_gui(&gfit);
}

void on_menuitem_mirrory_activate(GtkMenuItem *menuitem, gpointer user_data) {
	mirrory_gui(&gfit);
}

void on_mirrory_button_clicked(GtkToolButton *button, gpointer user_data) {
	mirrory_gui(&gfit);
}

/*************
 * RESAMPLE
 */

/* Resample */
void on_menuitem_resample_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded())
		siril_open_dialog("resample_dialog");
}

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
		undo_save_state(&gfit, "Processing: Resample (%g - %g)", sample[0] / 100.0,
				sample[1] / 100.0);
		verbose_resize_gaussian(&gfit, toX, toY, interpolation);
		update_used_memory();
		adjust_vport_size_to_image();
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
void on_menu_gray_crop_activate(GtkMenuItem *menuitem, gpointer user_data) {
	// if astrometry exists
	if (confirm_delete_wcs_keywords(&gfit)) {
		undo_save_state(&gfit, "Processing: Crop (x=%d, y=%d, w=%d, h=%d)",
				com.selection.x, com.selection.y, com.selection.w,
				com.selection.h);
		crop(&gfit, &com.selection);
		delete_selected_area();
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
		update_used_memory();
	}
}

void on_menu_gray_crop_seq_activate(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("crop_dialog");
}

/*** GUI for crop sequence */
void on_crop_Apply_clicked(GtkButton *button, gpointer user_data) {
	if (get_thread_run()) {
		siril_log_message(
				"Another task is already in progress, ignoring new request.\n");
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
	args->area = com.selection;
	args->prefix = gtk_entry_get_text(cropped_entry);

	set_cursor_waiting(TRUE);
	start_in_new_thread(crop_sequence, args);
}

void on_crop_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("crop_dialog");
}
