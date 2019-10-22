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

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
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
		undo_save_state(&asinh_gfit_backup, "Processing: Asinh Transformation: (stretch=%6.1lf, bp=%7.5lf)",
				asinh_stretch_value, asinh_black_value);
	}
	clearfits(&asinh_gfit_backup);
	set_cursor_waiting(FALSE);
}

static void asinh_recompute() {
	set_cursor_waiting(TRUE);
	copyfits(&asinh_gfit_backup, &gfit, CP_COPYA, -1);
	asinhlut(&gfit, asinh_stretch_value, asinh_black_value, asinh_rgb_space);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

int asinhlut(fits *fit, double beta, double offset, gboolean RGBspace) {
	int i, layer;
	WORD *buf[3] = { fit->pdata[RLAYER],
			fit->pdata[GLAYER], fit->pdata[BLAYER] };
	double norm;

	norm = get_normalized_value(fit);

	for (i = 0; i < fit->ry * fit->rx; i++) {
		double x, k;
		if (fit->naxes[2] > 1) {
			double r, g, b;

			r = (double) buf[RLAYER][i] / norm;
			g = (double) buf[GLAYER][i] / norm;
			b = (double) buf[BLAYER][i] / norm;
			/* RGB space */
			if (RGBspace)
				x = 0.2126 * r + 0.7152 * g + 0.0722 * b;
			else
				x = 0.3333 * r + 0.3333 * g + 0.3333 * b;
		} else {
			x = buf[RLAYER][i] / norm;
		}

		k = asinh(beta * x) / (x * asinh(beta));

		for (layer = 0; layer < fit->naxes[2]; ++layer) {
			double px = (double) buf[layer][i] / norm;
			px -= offset;
			px *= k;
			buf[layer][i] = round_to_WORD(px * norm);
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

void on_menuitem_asinh_activate(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("asinh_dialog");
}

void on_asinh_dialog_show(GtkWidget *widget, gpointer user_data) {
	asinh_startup();
	asinh_stretch_value = 1.0;
	asinh_black_value = 0.0;
	asinh_rgb_space = FALSE;
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_RGBspace")), asinh_rgb_space);
	gtk_range_set_value(GTK_RANGE(lookup_widget("scale_asinh")), asinh_stretch_value);
	gtk_range_set_value(GTK_RANGE(lookup_widget("black_point_asinh")), asinh_black_value);
}

void apply_asinh_cancel() {
	asinh_close(TRUE);
	siril_close_dialog("asinh_dialog");
}

void on_asinh_cancel_clicked(GtkButton *button, gpointer user_data) {
	apply_asinh_cancel();
}

static void apply_asinh_changes() {
	gboolean status = (asinh_stretch_value != 1.0) || (asinh_black_value != 0.0) || asinh_rgb_space;
	asinh_close(!status);
}
void on_asinh_ok_clicked(GtkButton *button, gpointer user_data) {
	apply_asinh_changes();
	siril_close_dialog("asinh_dialog");
}

void on_asinh_dialog_close(GtkDialog *dialog, gpointer user_data) {
	apply_asinh_changes();
}

gboolean on_scale_asinh_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	asinh_stretch_value = gtk_range_get_value(GTK_RANGE(widget));
	asinh_recompute();
	return FALSE;
}

gboolean on_scale_asinh_key_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	asinh_stretch_value = gtk_range_get_value(GTK_RANGE(widget));
	asinh_recompute();
	return FALSE;
}

gboolean on_black_point_asinh_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	asinh_black_value = gtk_range_get_value(GTK_RANGE(widget));
	asinh_recompute();
	return FALSE;
}

gboolean on_black_point_asinh_key_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	asinh_black_value = gtk_range_get_value(GTK_RANGE(widget));
	asinh_recompute();
	return FALSE;
}

void on_asinh_RGBspace_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	asinh_rgb_space = gtk_toggle_button_get_active(togglebutton);
	asinh_recompute();
}

void on_asinh_undo_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	asinh_stretch_value = 1.0;
	asinh_black_value = 0.0;
	asinh_rgb_space = FALSE;
	GtkToggleButton *check_button = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_RGBspace"));
	g_signal_handlers_block_by_func(check_button, on_asinh_RGBspace_toggled, NULL);
	gtk_toggle_button_set_active(check_button, asinh_rgb_space);
	g_signal_handlers_unblock_by_func(check_button, on_asinh_RGBspace_toggled, NULL);
	gtk_range_set_value(GTK_RANGE(lookup_widget("scale_asinh")), asinh_stretch_value);
	gtk_range_set_value(GTK_RANGE(lookup_widget("black_point_asinh")), asinh_black_value);
	copyfits(&asinh_gfit_backup, &gfit, CP_COPYA, -1);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}
