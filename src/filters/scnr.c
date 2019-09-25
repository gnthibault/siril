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
#include "core/undo.h"
#include "core/processing.h"
#include "algos/colors.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "gui/progress_and_log.h"
#include "gui/callbacks.h"
#include "gui/histogram.h"
#include "gui/dialogs.h"

#include "scnr.h"

// idle function executed at the end of the scnr processing
static gboolean end_scnr(gpointer p) {
	stop_processing_thread();
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	update_gfit_histogram_if_needed();
	set_cursor_waiting(FALSE);
	update_used_memory();

	return FALSE;
}

/* Subtractive Chromatic Noise Reduction.
 * No unprotected GTK+ calls can go there. */
gpointer scnr(gpointer p) {
	struct scnr_data *args = (struct scnr_data *) p;
	WORD *buf[3] = { args->fit->pdata[RLAYER], args->fit->pdata[GLAYER],
			args->fit->pdata[BLAYER] };
	double m;
	int nbdata = args->fit->rx * args->fit->ry;
	int i;
	struct timeval t_start, t_end;

	siril_log_color_message(_("SCNR: processing...\n"), "red");
	gettimeofday(&t_start, NULL);

	WORD norm = get_normalized_value(args->fit);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
#endif
	for (i = 0; i < nbdata; i++) {
		double red = (double) buf[RLAYER][i] / norm;
		double green = (double) buf[GLAYER][i] / norm;
		double blue = (double) buf[BLAYER][i] / norm;
		double x, y, z, L, a, b;

		if (args->preserve) {
			rgb_to_xyz(red, green, blue, &x, &y, &z);
			xyz_to_LAB(x, y, z, &L, &a, &b);
		}
		switch (args->type) {
		case 0:
			m = 0.5 * (red + blue);
			green = min(green, m);
			break;
		case 1:
			m = max(red, blue);
			green = min(green, m);
			break;
		case 2:
			m = max(red, blue);
			green = (green * (1.0 - args->amount) * (1.0 - m)) + (m * green);
			break;
		case 3:
			m = min(1.0, red + blue);
			green = (green * (1.0 - args->amount) * (1.0 - m)) + (m * green);
		}
		if (args->preserve) {
			double tmp;
			rgb_to_xyz(red, green, blue, &x, &y, &z);
			xyz_to_LAB(x, y, z, &tmp, &a, &b);
			LAB_to_xyz(L, a, b, &x, &y, &z);
			xyz_to_rgb(x, y, z, &red, &green, &blue);
		}
		buf[RLAYER][i] = round_to_WORD(red * (double) norm);
		buf[GLAYER][i] = round_to_WORD(green * (double) norm);
		buf[BLAYER][i] = round_to_WORD(blue * (double) norm);
	}

	invalidate_stats_from_fit(args->fit);
	free(args);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	siril_add_idle(end_scnr, NULL);
	return GINT_TO_POINTER(0);
}

void on_removegreen_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded() && isrgb(&gfit))
		siril_open_dialog("SCNR_dialog");
}

void on_SCNR_dialog_show(GtkWidget *widget, gpointer user_data) {
	GtkComboBox *comboscnr = GTK_COMBO_BOX(
			gtk_builder_get_object(builder, "combo_scnr"));
	int type = gtk_combo_box_get_active(comboscnr);

	if (type == -1)
		gtk_combo_box_set_active(comboscnr, 0);
}

void on_SCNR_Apply_clicked(GtkButton *button, gpointer user_data) {
	/* Type 0: Average Neutral protection
	 * Type 1: Maximum Neutral protection
	 */
	int type = gtk_combo_box_get_active(
			GTK_COMBO_BOX(gtk_builder_get_object(builder, "combo_scnr")));
	GtkToggleButton *light_button = GTK_TOGGLE_BUTTON(
			gtk_builder_get_object(builder, "preserve_light"));
	gboolean preserve = gtk_toggle_button_get_active(light_button);
	double amount = gtk_range_get_value(
			GTK_RANGE(gtk_builder_get_object(builder, "scale_scnr")));

	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return;
	}

	struct scnr_data *args = malloc(sizeof(struct scnr_data));
	undo_save_state(&gfit, "Processing: SCNR (type=%d, amount=%0.2lf, preserve=%s)",
			type, amount, preserve ? "true" : "false");

	args->fit = &gfit;
	args->type = type;
	args->amount = amount;
	args->preserve = preserve;
	set_cursor_waiting(TRUE);
	start_in_new_thread(scnr, args);
}

void on_SCNR_cancel_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("SCNR_dialog");
}

void on_combo_scnr_changed(GtkComboBoxText *box, gpointer user_data) {
	int type = gtk_combo_box_get_active(
			GTK_COMBO_BOX(gtk_builder_get_object(builder, "combo_scnr")));
	GtkScale *scale = GTK_SCALE(lookup_widget("scale_scnr"));
	GtkLabel *label = GTK_LABEL(lookup_widget("label56"));

	gtk_widget_set_sensitive(GTK_WIDGET(scale), type > 1);
	gtk_widget_set_sensitive(GTK_WIDGET(label), type > 1);
}

