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

#include "core/siril.h"
#include "core/proto.h"
#include "core/undo.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "algos/colors.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/callbacks.h"
#include "gui/histogram.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"

#include "saturation.h"

static gboolean satu_preserve_bkg;
static double satu_amount;
static int satu_hue_type;
static gboolean satu_show_preview;

static void satu_startup() {
	copy_gfit_to_backup();
	satu_amount = 0.0;
	satu_hue_type = 6;
}

static void satu_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		siril_preview_hide();
	} else {
		undo_save_state(get_preview_gfit_backup(),
				"Processing: Saturation enhancement (amount=%4.2lf)", satu_amount);
	}
	clear_backup();
	set_cursor_waiting(FALSE);
}


static void apply_satu_changes() {
	gboolean status = satu_amount != 0.0;
	satu_close(!status);
}

static int satu_update_preview() {
	if (get_thread_run()) {
		siril_debug_print(_("Another task is already in progress, ignoring new request.\n"));
		return 1;
	}

	set_cursor_waiting(TRUE);

	struct enhance_saturation_data *args = malloc(sizeof(struct enhance_saturation_data));

	switch (satu_hue_type) {
		case 0:		// Pink-Red to Red-Orange
			args->h_min = 346.0;
			args->h_max = 20.0;
			break;
		case 1:		// Orange-Brown to Yellow
			args->h_min = 21.0;
			args->h_max = 60.0;
			break;
		case 2:		// Yellow-Green to Green-Cyan
			args->h_min = 61.0;
			args->h_max = 200.0;
			break;
		case 3:		// Cyan
			args->h_min = 170.0;
			args->h_max = 200.0;
			break;
		case 4:		// Cyan-Blue to Blue-Magenta
			args->h_min = 201.0;
			args->h_max = 280.0;
			break;
		case 5:		// Magenta to Pink
			args->h_min = 281.0;
			args->h_max = 345.0;
			break;
		default:
		case 6:		// Global
			args->h_min = 0.0;
			args->h_max = 360.0;
	}

	args->input = satu_show_preview ? get_preview_gfit_backup() : &gfit;
	args->output = &gfit;
	args->coeff = satu_amount;
	args->preserve = satu_preserve_bkg;
	start_in_new_thread(enhance_saturation, args);

	return 0;
}

void on_satu_cancel_clicked(GtkButton *button, gpointer user_data) {
	satu_close(TRUE);
	siril_close_dialog("satu_dialog");
}

void on_satu_apply_clicked(GtkButton *button, gpointer user_data) {
	if (satu_show_preview == FALSE) {
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = satu_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}

	apply_satu_changes();
	siril_close_dialog("satu_dialog");
}

void on_satu_dialog_close(GtkDialog *dialog, gpointer user_data) {
	apply_satu_changes();
}

gpointer enhance_saturation_ushort(gpointer p) {
	struct enhance_saturation_data *args = (struct enhance_saturation_data *) p;
	double bg = 0;

	WORD *in[3] = { args->input->pdata[RLAYER], args->input->pdata[GLAYER],
		args->input->pdata[BLAYER] };
	WORD *out[3] = { args->output->pdata[RLAYER], args->output->pdata[GLAYER],
		args->output->pdata[BLAYER] };

	args->h_min /= 360.0;
	args->h_max /= 360.0;
	if (args->preserve) {
		imstats *stat = statistics(NULL, -1, args->input, GLAYER, NULL, STATS_BASIC, TRUE);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			siril_add_idle(end_generic, args);
			free(args);
			return GINT_TO_POINTER(1);
		}
		bg = stat->median + stat->sigma;
		bg /= stat->normValue;
		free_stats(stat);
	}

	size_t i, n = args->input->naxes[0] * args->input->naxes[1];
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
#endif
	for (i = 0; i < n; i++) {
		double h, s, l;
		double r = (double) in[RLAYER][i] / USHRT_MAX_DOUBLE;
		double g = (double) in[GLAYER][i] / USHRT_MAX_DOUBLE;
		double b = (double) in[BLAYER][i] / USHRT_MAX_DOUBLE;
		rgb_to_hsl(r, g, b, &h, &s, &l);
		if (l > bg) {
			if (args->h_min > args->h_max) {// Red case. TODO: find a nicer way to do it
				if (h >= args->h_min || h <= args->h_max) {
					s += (s * args->coeff);
				}
			} else {
				if (h >= args->h_min && h <= args->h_max) {
					s += (s * args->coeff);
				}
			}
			if (s < 0.0)
				s = 0.0;
			else if (s > 1.0)
				s = 1.0;
		}
		hsl_to_rgb(h, s, l, &r, &g, &b);
		out[RLAYER][i] = round_to_WORD(r * USHRT_MAX_DOUBLE);
		out[GLAYER][i] = round_to_WORD(g * USHRT_MAX_DOUBLE);
		out[BLAYER][i] = round_to_WORD(b * USHRT_MAX_DOUBLE);
	}
	invalidate_stats_from_fit(args->output);

	siril_add_idle(end_generic, args);
	free(args);

	return GINT_TO_POINTER(0);
}

static gpointer enhance_saturation_float(gpointer p) {
	struct enhance_saturation_data *args = (struct enhance_saturation_data *) p;
	float bg = 0;

	float *in[3] = { args->input->fpdata[RLAYER], args->input->fpdata[GLAYER],
		args->input->fpdata[BLAYER] };
	float *out[3] = { args->output->fpdata[RLAYER], args->output->fpdata[GLAYER],
		args->output->fpdata[BLAYER] };

	args->h_min /= 60.0;
	args->h_max /= 60.0;
	if (args->preserve) {
		imstats *stat = statistics(NULL, -1, args->input, GLAYER, NULL, STATS_BASIC, TRUE);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			siril_add_idle(end_generic, args);
			free(args);
			return GINT_TO_POINTER(1);
		}
		bg = stat->median + stat->sigma;
		bg /= stat->normValue;
		free_stats(stat);
	}

	float s_mult = 1.f + args->coeff;
	gboolean red_case = args->h_min > args->h_max;
	float h_min = args->h_min;
	float h_max = args->h_max;

	size_t i, n = args->input->naxes[0] * args->input->naxes[1];
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(dynamic, args->input->rx * 16)
#endif
	for (i = 0; i < n; i++) {
		float h, s, l;
		float r = in[RLAYER][i];
		float g = in[GLAYER][i];
		float b = in[BLAYER][i];

		rgb_to_hsl_float_sat(r, g, b, bg, &h, &s, &l);
		if (l > bg) {
			if (red_case) {// Red case. TODO: find a nicer way to do it
				if (h >= h_min || h <= h_max) {
					s *= s_mult;
				}
			} else {
				if (h >= h_min && h <= h_max) {
					s *= s_mult;
				}
			}
			s = s > 1.f ? 1.f : s;
			hsl_to_rgb_float_sat(h, s, l, &r, &g, &b);
		}
		out[RLAYER][i] = r;
		out[GLAYER][i] = g;
		out[BLAYER][i] = b;
	}
	invalidate_stats_from_fit(args->output);
	siril_add_idle(end_generic, args);
	free(args);

	return GINT_TO_POINTER(0);
}

gpointer enhance_saturation(gpointer p) {
	struct enhance_saturation_data *args = (struct enhance_saturation_data *) p;

	if (args->input->type == DATA_USHORT) {
		return enhance_saturation_ushort(args);
	} else if (args->input->type == DATA_FLOAT) {
		return enhance_saturation_float(args);
	}
	return GINT_TO_POINTER(-1);
}

/** callbacks **/

void on_menuitem_satu_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (!single_image_is_loaded() || !isrgb(&gfit))
		return;

	siril_open_dialog("satu_dialog");
}

void on_satu_dialog_show(GtkWidget *widget, gpointer user_data) {
	satu_startup();
	satu_amount = 0.0;
	satu_hue_type = 6;
	satu_preserve_bkg = TRUE;

	set_notify_block(TRUE);
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_saturation")), satu_hue_type);
	gtk_range_set_value(GTK_RANGE(lookup_widget("scale_satu")), satu_amount);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("preserve_bg")), satu_preserve_bkg);
	set_notify_block(FALSE);

	satu_show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("satu_preview")));
}

void on_preserve_bg_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	satu_preserve_bkg = gtk_toggle_button_get_active(togglebutton);

	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = 	satu_update_preview;
	param->show_preview = satu_show_preview;
	notify_update((gpointer) param);
}

void on_combo_saturation_changed(GtkComboBox* box, gpointer user_data) {
	satu_hue_type = gtk_combo_box_get_active(box);

	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = 	satu_update_preview;
	param->show_preview = satu_show_preview;
	notify_update((gpointer) param);
}

void on_satu_undo_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	satu_preserve_bkg = TRUE;
	satu_amount = 0.0;
	GtkToggleButton *check_button = GTK_TOGGLE_BUTTON(lookup_widget("preserve_bg"));

	set_notify_block(TRUE);
	gtk_toggle_button_set_active(check_button, satu_preserve_bkg);
	gtk_range_set_value(GTK_RANGE(lookup_widget("scale_satu")), satu_amount);
	set_notify_block(FALSE);

	copy_backup_to_gfit();
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void apply_satu_cancel() {
	satu_close(TRUE);
}

/*** adjusters **/
void on_spin_satu_value_changed(GtkSpinButton *button, gpointer user_data) {
	satu_amount = gtk_spin_button_get_value(button);

	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = 	satu_update_preview;
	param->show_preview = satu_show_preview;
	notify_update((gpointer) param);
}

void on_satu_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	if (satu_show_preview == TRUE) {
		/* if user click very fast */
		waiting_for_thread();
		siril_preview_hide();
	} else {
		copy_gfit_to_backup();

		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = satu_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
	satu_show_preview = !satu_show_preview;
}
