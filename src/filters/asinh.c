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

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/preview_timer.h"
#include "core/undo.h"

#include "asinh.h"

static gboolean asinh_rgb_space = FALSE;
static double asinh_stretch_value = 1.0, asinh_black_value = 0.0;
static fits asinh_gfit_backup;

static void asinh_startup() {
	copyfits(&gfit, &asinh_gfit_backup, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
}

static void asinh_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		copyfits(&asinh_gfit_backup, &gfit, CP_COPYA, -1);
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	} else {
		invalidate_stats_from_fit(&gfit);
		undo_save_state(&asinh_gfit_backup,
				"Processing: Asinh Transformation: (stretch=%6.1lf, bp=%7.5lf)",
				asinh_stretch_value, asinh_black_value);
	}
	clearfits(&asinh_gfit_backup);
	set_cursor_waiting(FALSE);
}

static int asinh_update_preview() {
	copyfits(&asinh_gfit_backup, &gfit, CP_COPYA, -1);
	asinhlut(&gfit, asinh_stretch_value, asinh_black_value, asinh_rgb_space);
	return 0;
}

static void asinh_recompute() {
	if (asinh_stretch_value == 0) {
		// sometimes this happens in gui
		return;
	}
	set_cursor_waiting(TRUE);
	asinh_update_preview();
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

int asinhlut(fits *fit, double beta, double offset, gboolean RGBspace) {
	int i;
	WORD *buf[3] = { fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	double norm, asinh_beta, factor_red, factor_green, factor_blue;

	norm = get_normalized_value(fit);
	asinh_beta = asinh(beta);
	factor_red = RGBspace ? 0.2126 : 0.3333;
	factor_green = RGBspace ? 0.7152 : 0.3333;
	factor_blue = RGBspace ? 0.0722 : 0.3333;

	if (fit->naxes[2] > 1) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->rx * 16)
#endif
		for (i = 0; i < fit->ry * fit->rx; i++) {
			double x, k;
			double r, g, b;

			r = (double) buf[RLAYER][i] / norm;
			g = (double) buf[GLAYER][i] / norm;
			b = (double) buf[BLAYER][i] / norm;

			x = factor_red * r + factor_green * g + factor_blue * b;

			k = (x == 0.0) ? 0.0 : asinh(beta * x) / (x * asinh_beta);

			buf[RLAYER][i] = round_to_WORD((r - offset) * k * norm);
			buf[GLAYER][i] = round_to_WORD((g - offset) * k * norm);
			buf[BLAYER][i] = round_to_WORD((b - offset) * k * norm);
		}
	} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, fit->rx * 16)
#endif
		for (i = 0; i < fit->ry * fit->rx; i++) {
			double x, k;
			x = buf[RLAYER][i] / norm;
			k = (x == 0.0) ? 0.0 : asinh(beta * x) / (x * asinh_beta);
			buf[RLAYER][i] = round_to_WORD((x - offset) * k * norm);
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

static void apply_asinh_changes() {
	gboolean status = (asinh_stretch_value != 1.0) || (asinh_black_value != 0.0) || asinh_rgb_space;
	asinh_close(!status);
}

void apply_asinh_cancel() {
	asinh_close(TRUE);
	siril_close_dialog("asinh_dialog");
}

/*** callbacks **/

void on_menuitem_asinh_activate(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("asinh_dialog");
}

void on_asinh_dialog_show(GtkWidget *widget, gpointer user_data) {
	GtkSpinButton *spin_stretch = GTK_SPIN_BUTTON(lookup_widget("spin_asinh"));
	GtkSpinButton *spin_black_p = GTK_SPIN_BUTTON(lookup_widget("black_point_spin_asinh"));
	GtkToggleButton *toggle_rgb = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_RGBspace"));

	asinh_startup();
	asinh_stretch_value = 1.0;
	asinh_black_value = 0.0;
	asinh_rgb_space = FALSE;
	gtk_toggle_button_set_active(toggle_rgb, asinh_rgb_space);
	gtk_spin_button_set_value(spin_stretch, asinh_stretch_value);
	gtk_spin_button_set_value(spin_black_p, asinh_black_value);
	gtk_spin_button_set_increments(spin_stretch, 0.001, 0.01);
}

void on_asinh_cancel_clicked(GtkButton *button, gpointer user_data) {
	apply_asinh_cancel();
}

void on_asinh_ok_clicked(GtkButton *button, gpointer user_data) {
	apply_asinh_changes();
	siril_close_dialog("asinh_dialog");
}

void on_asinh_dialog_close(GtkDialog *dialog, gpointer user_data) {
	apply_asinh_changes();
}

void on_asinh_RGBspace_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	asinh_rgb_space = gtk_toggle_button_get_active(togglebutton);
	asinh_recompute();
}

void on_asinh_undo_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *spin_stretch = GTK_SPIN_BUTTON(lookup_widget("spin_asinh"));
	GtkSpinButton *spin_black_p = GTK_SPIN_BUTTON(lookup_widget("black_point_spin_asinh"));
	GtkToggleButton *toggle_rgb = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_RGBspace"));
	asinh_stretch_value = 1.0;
	asinh_black_value = 0.0;
	asinh_rgb_space = FALSE;

	set_cursor_waiting(TRUE);
	g_signal_handlers_block_by_func(toggle_rgb, on_asinh_RGBspace_toggled, NULL);
	gtk_toggle_button_set_active(toggle_rgb, asinh_rgb_space);
	g_signal_handlers_unblock_by_func(toggle_rgb, on_asinh_RGBspace_toggled, NULL);

	gtk_spin_button_set_value(spin_stretch, asinh_stretch_value);
	gtk_spin_button_set_value(spin_black_p, asinh_black_value);

	copyfits(&asinh_gfit_backup, &gfit, CP_COPYA, -1);

	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

/*** adjusters **/
void on_spin_asinh_value_changed(GtkSpinButton *button, gpointer user_data) {
	asinh_stretch_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	notify_update((gpointer) param);
}

void on_black_point_spin_asinh_value_changed(GtkSpinButton *button, gpointer user_data) {
	asinh_black_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = asinh_update_preview;
	notify_update((gpointer) param);
}
