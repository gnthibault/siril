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


#include "core/siril.h"
#include "core/proto.h"

#include "core/undo.h"
#include "core/processing.h"
#include "io/single_image.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"
#include "opencv/opencv.h"

#include "deconv.h"

gpointer RTdeconv(gpointer p) {
	struct deconv_data *args = (struct deconv_data *) p;
	struct timeval t_start, t_end;

	siril_log_color_message(_("Deconvolution: processing...\n"), "red");
	gettimeofday(&t_start, NULL);

	deconvolution(args);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	siril_add_idle(end_generic, args);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

/************************ GUI for deconvolution ***********************/
void on_menuitem_deconvolution_activate(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("deconvolution_dialog");
}

void on_deconvolution_cancel_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("deconvolution_dialog");
}

void on_deconvolution_reset_clicked(GtkButton *button, gpointer user_data) {
	GtkRange *threshold, *sigma, *corner_radius, *iterations;
	GtkToggleButton *auto_limit, *auto_threshold;

	threshold = GTK_RANGE(lookup_widget("scale_deconv_threshold"));
	sigma = GTK_RANGE(lookup_widget("scale_deconv_radius"));
	corner_radius = GTK_RANGE(lookup_widget("scale_deconv_corner"));
	iterations = GTK_RANGE(lookup_widget("scale_deconv_iterations"));
	auto_limit = GTK_TOGGLE_BUTTON(lookup_widget("toggle_deconv_auto"));
	auto_threshold = GTK_TOGGLE_BUTTON(lookup_widget("toggle_deconv_trheshold"));

	gtk_range_set_value(threshold, 20);
	gtk_range_set_value(sigma, 1.0);
	gtk_range_set_value(corner_radius, 0.0);
	gtk_range_set_value(iterations, 20);
	gtk_toggle_button_set_active(auto_limit, TRUE);
	gtk_toggle_button_set_active(auto_threshold, FALSE);
}

void on_toggle_deconv_trheshold_toggled(GtkToggleButton *button, gpointer user_data) {
	gtk_widget_set_sensitive(lookup_widget("scale_deconv_threshold"), !gtk_toggle_button_get_active(button));
	gtk_widget_set_sensitive(lookup_widget("spin_deconv_threshold"), !gtk_toggle_button_get_active(button));
}

void on_deconvolution_Apply_clicked(GtkButton *button, gpointer user_data) {
	GtkRange *threshold, *sigma, *corner_radius, *iterations;
	GtkToggleButton *auto_limit, *auto_threshold;

	threshold = GTK_RANGE(lookup_widget("scale_deconv_threshold"));
	sigma = GTK_RANGE(lookup_widget("scale_deconv_radius"));
	corner_radius = GTK_RANGE(lookup_widget("scale_deconv_corner"));
	iterations = GTK_RANGE(lookup_widget("scale_deconv_iterations"));
	auto_limit = GTK_TOGGLE_BUTTON(lookup_widget("toggle_deconv_auto"));
	auto_threshold = GTK_TOGGLE_BUTTON(lookup_widget("toggle_deconv_trheshold"));

	struct deconv_data *args = malloc(sizeof(struct deconv_data));

	set_cursor_waiting(TRUE);

	args->fit = &gfit;
	args->contrast_threshold = (size_t)gtk_range_get_value(threshold);
	args->sigma = gtk_range_get_value(sigma);
	args->corner_radius = gtk_range_get_value(corner_radius);
	args->iterations = (size_t)gtk_range_get_value(iterations);
	args->auto_limit = gtk_toggle_button_get_active(auto_limit);
	args->auto_contrast_threshold = gtk_toggle_button_get_active(auto_threshold);
	if (args->fit->type == DATA_USHORT) {
		args->clip = (args->fit->maxi <= 0) ? USHRT_MAX_DOUBLE : args->fit->maxi;
	} else {
		args->clip = (args->fit->maxi <= 0) ? USHRT_MAX_DOUBLE : args->fit->maxi * USHRT_MAX_DOUBLE;
	}
	undo_save_state(args->fit, "Processing: Deconv. (iter=%d, sig=%.3f)", args->iterations, args->sigma);

	start_in_new_thread(RTdeconv, args);
}
