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

#include "algos/statistics.h"
#include "core/undo.h"
#include "core/processing.h"
#include "io/single_image.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"
#include "gui/preview_timer.h"

#include "deconv.h"

static int deconv_threshold, deconv_iterations;
static double deconv_radius, deconv_boost;
static gboolean deconv_auto_thres, deconv_auto_iter;
static fits deconv_gfit_backup;

static void deconv_startup() {
	copyfits(&gfit, &deconv_gfit_backup, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	deconv_threshold = 20;
	deconv_iterations = 20;
	deconv_radius = 1.0;
	deconv_boost = 0.0;
	deconv_auto_thres = TRUE;
	deconv_auto_iter = TRUE;
}

static void deconv_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		copyfits(&deconv_gfit_backup, &gfit, CP_COPYA, -1);
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	} else {
		invalidate_stats_from_fit(&gfit);
		undo_save_state(&deconv_gfit_backup,
				"Processing: Deconv. (iter=%d, sig=%.3f)", deconv_iterations,
				deconv_radius);

	}
	clearfits(&deconv_gfit_backup);
	set_cursor_waiting(FALSE);
}

int deconv_update_preview() {
	copyfits(&deconv_gfit_backup, &gfit, CP_COPYA, -1);

	struct deconv_data *args = malloc(sizeof(struct deconv_data));

	set_cursor_waiting(TRUE);

	args->fit = &gfit;
	args->clip = args->fit->maxi;
	args->contrast_threshold = (size_t)deconv_threshold;
	args->sigma = deconv_radius;
	args->corner_radius = deconv_boost;
	args->iterations = (size_t)deconv_iterations;
	args->auto_limit = deconv_auto_iter;
	args->auto_contrast_threshold = deconv_auto_thres;

	start_in_new_thread(RTdeconv, args);
	return 0;
}


gpointer RTdeconv(gpointer p) {
	struct deconv_data *args = (struct deconv_data *) p;

	deconvolution(args);

	siril_add_idle(end_generic, args);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}


/** callbacks **/
void on_deconvolution_dialog_show(GtkWidget *widget, gpointer user_data) {
	deconv_startup();

	set_notify_block(TRUE);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_deconv_threshold")), deconv_threshold);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_deconv_radius")),	deconv_radius);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_deconv_boost")),	deconv_boost);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_deconv_iterations")),	deconv_iterations);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_deconv_trheshold")),	deconv_auto_iter);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_deconv_auto")),	deconv_auto_thres);
	set_notify_block(FALSE);

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = deconv_update_preview;
	notify_update((gpointer) param);
}


void on_menuitem_deconvolution_activate(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("deconvolution_dialog");
}

void on_deconvolution_cancel_clicked(GtkButton *button, gpointer user_data) {
	deconv_close(TRUE);
	siril_close_dialog("deconvolution_dialog");
}

void on_deconvolution_reset_clicked(GtkButton *button, gpointer user_data) {
	deconv_threshold = 20;
	deconv_iterations = 20;
	deconv_radius = 1.0;
	deconv_boost = 0.0;
	deconv_auto_thres = TRUE;
	deconv_auto_iter = TRUE;

	set_cursor_waiting(TRUE);
	set_notify_block(TRUE);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_deconv_threshold")), deconv_threshold);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_deconv_radius")),	deconv_radius);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_deconv_boost")),	deconv_boost);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_deconv_iterations")),	deconv_iterations);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_deconv_trheshold")),	deconv_auto_iter);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("toggle_deconv_auto")),	deconv_auto_thres);
	set_notify_block(FALSE);

	copyfits(&deconv_gfit_backup, &gfit, CP_COPYA, -1);

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = deconv_update_preview;
	notify_update((gpointer) param);
}

void on_deconvolution_Apply_clicked(GtkButton *button, gpointer user_data) {
	deconv_close(FALSE);
	siril_close_dialog("deconvolution_dialog");
}

void on_deconvolution_dialog_close(GtkDialog *dialog, gpointer user_data) {
	deconv_close(TRUE);
}

void on_spin_deconv_threshold_value_changed(GtkSpinButton *button, gpointer user_data) {
	deconv_threshold = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = deconv_update_preview;
	notify_update((gpointer) param);
}

void on_spin_deconv_radius_value_changed(GtkSpinButton *button, gpointer user_data) {
	deconv_radius = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = deconv_update_preview;
	notify_update((gpointer) param);
}

void on_spin_deconv_boost_value_changed(GtkSpinButton *button, gpointer user_data) {
	deconv_boost = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = deconv_update_preview;
	notify_update((gpointer) param);
}

void on_spin_deconv_iterations_value_changed(GtkSpinButton *button, gpointer user_data) {
	deconv_iterations = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = deconv_update_preview;
	notify_update((gpointer) param);
}

void on_toggle_deconv_trheshold_toggled(GtkToggleButton *button, gpointer user_data) {
	gtk_widget_set_sensitive(lookup_widget("scale_deconv_threshold"), !gtk_toggle_button_get_active(button));
	gtk_widget_set_sensitive(lookup_widget("spin_deconv_threshold"), !gtk_toggle_button_get_active(button));

	deconv_auto_thres = gtk_toggle_button_get_active(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = deconv_update_preview;
	notify_update((gpointer) param);
}

void on_toggle_deconv_auto_toggled(GtkToggleButton *button, gpointer user_data) {
	deconv_auto_iter = gtk_toggle_button_get_active(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = deconv_update_preview;
	notify_update((gpointer) param);
}
