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

#include <gsl/gsl_fit.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/optimize_utils.h"
#include "core/undo.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "gui/message_dialog.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"

#include "linear_match.h"

static int find_linear_coeff_ushort(fits *target_fit, fits *reference_fit, double low,
		double high, double *a, double *b, gchar **error) {
	int count = 0;
	double c0, c1, cov00, cov01, cov11, sumsq;
	size_t ref_size = reference_fit->rx * reference_fit->ry;

	if (memcmp(target_fit->naxes, reference_fit->naxes, sizeof target_fit->naxes)) {
		gchar *err = siril_log_message(_("Images must have same dimensions.\n"));
		if (error) {
			*error = err;
		}
		return -1;
	}

	low *= USHRT_MAX_DOUBLE;
	high *= USHRT_MAX_DOUBLE;

	siril_log_color_message(_("Linear fit functions:\n"), "green");
	for (int channel = 0; channel < reference_fit->naxes[2]; channel++) {
		for (size_t i = 0; i < ref_size; i++) {
			if (inInterval(reference_fit->pdata[channel][i], low, high)) {
				count ++;
			}
		}

		double *x = malloc(count * sizeof(double));
		double *y = malloc(count * sizeof(double));
		for (size_t i = 0, j = 0; i < ref_size; i++) {
			if (inInterval(reference_fit->pdata[channel][i], low, high)) {
				y[j] = (double) reference_fit->pdata[channel][i] / USHRT_MAX_DOUBLE;
				if (target_fit->type == DATA_FLOAT) {
					x[j] = (double) target_fit->fpdata[channel][i];
				} else {
					x[j] = (double) target_fit->pdata[channel][i] / USHRT_MAX_DOUBLE;
				}
				j++;
			}
		}
		gsl_fit_linear(x, 1, y, 1, count, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
		siril_log_color_message("y_0 = %e + %e*x_0 (%d)\n", "blue", c0, c1, count);
		free(x);
		free(y);
		a[channel] = c1;
		b[channel] = c0;
	}
	return 0;
}

static int find_linear_coeff_float(fits *target_fit, fits *reference_fit, double low,
		double high, double *a, double *b, gchar **error) {
	int count = 0;
	double c0, c1, cov00, cov01, cov11, sumsq;
	size_t ref_size = reference_fit->rx * reference_fit->ry;

	if (memcmp(target_fit->naxes, reference_fit->naxes, sizeof target_fit->naxes)) {
		gchar *err = siril_log_message(_("Images must have same dimensions.\n"));
		if (error) {
			*error = err;
		}
		return -1;
	}

	siril_log_color_message(_("Linear fit functions:\n"), "green");
	for (int channel = 0; channel < reference_fit->naxes[2]; channel++) {
		for (size_t i = 0; i < ref_size; i++) {
			if (inInterval(reference_fit->fpdata[channel][i], low, high)) {
				count ++;
			}
		}

		double *x = malloc(count * sizeof(double));
		double *y = malloc(count * sizeof(double));
		for (size_t i = 0, j = 0; i < ref_size; i++) {
			if (inInterval(reference_fit->fpdata[channel][i], low, high)) {
				y[j] = (double) reference_fit->fpdata[channel][i];
				if (target_fit->type == DATA_FLOAT) {
					x[j] = (double) target_fit->fpdata[channel][i];
				} else {
					x[j] = (double) target_fit->pdata[channel][i] / USHRT_MAX_DOUBLE;
				}
				j++;
			}
		}
		gsl_fit_linear(x, 1, y, 1, count, &c0, &c1, &cov00, &cov01, &cov11,	&sumsq);
		siril_log_color_message("y_0 = %e + %e*x_0 (%d)\n", "blue", c0, c1, count);
		free(x);
		free(y);
		a[channel] = c1;
		b[channel] = c0;
	}
	return 0;
}

int find_linear_coeff(fits *target_fit, fits *reference_fit, double low,
		double high, double *a, double *b, gchar **error) {
	if (reference_fit->type == DATA_USHORT) {
		return find_linear_coeff_ushort(target_fit, reference_fit, low, high, a, b, error);
	} else if (reference_fit->type == DATA_FLOAT) {
		return find_linear_coeff_float(target_fit, reference_fit, low, high, a, b, error);
	}
	return 1;
}

static void apply_linear_to_fits_ushort(fits *fit, double *a, double *b) {
	size_t size = fit->rx * fit->ry;

	invalidate_stats_from_fit(&gfit);
	for (int channel = 0; channel < fit->naxes[2]; channel++) {
		for (size_t i = 0; i < size; i++) {
			fit->pdata[channel][i] = round_to_WORD(fit->pdata[channel][i] * a[channel]+ b[channel] * USHRT_MAX_DOUBLE);
		}
	}
}

static void apply_linear_to_fits_float(fits *fit, double *a, double *b) {
	size_t size = fit->rx * fit->ry * fit->naxes[2];

	invalidate_stats_from_fit(&gfit);
	for (int channel = 0; channel < fit->naxes[2]; channel++) {
		for (size_t i = 0; i < size; i++) {
			fit->fpdata[channel][i] = fit->fpdata[channel][i] * a[channel]
					+ b[channel];
		}
	}
}

void apply_linear_to_fits(fits *fit, double *a, double *b) {
	if (fit->type == DATA_USHORT) {
		apply_linear_to_fits_ushort(fit, a, b);
	} else if (fit->type == DATA_FLOAT) {
		apply_linear_to_fits_float(fit, a, b);
	}
}

static gchar *get_reference_filename() {
	GtkFileChooser *linearmatch_ref = GTK_FILE_CHOOSER(lookup_widget("reference_filechooser_linearmatch"));
	return gtk_file_chooser_get_filename(linearmatch_ref);
}

static gdouble get_high_rejection() {
	GtkSpinButton *button = GTK_SPIN_BUTTON(lookup_widget("spin_linearmatch_high"));

	return gtk_spin_button_get_value(button);
}

static gdouble get_low_rejection() {
	GtkSpinButton *button = GTK_SPIN_BUTTON(lookup_widget("spin_linearmatch_low"));

	return gtk_spin_button_get_value(button);
}

/*** callbacks **/

void on_menu_linearmatch_activate(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("linearmatch_dialog");
}

void on_linearmatch_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("linearmatch_dialog");
}

void on_linearmatch_apply_clicked(GtkButton *button, gpointer user_data) {
	if (!single_image_is_loaded()) return;
	gchar *filename = get_reference_filename();
	if (filename) {
		gchar *error;
		fits ref = { 0 };
		double a[3] = { 0.0 }, b[3] = { 0.0 };
		double low = get_low_rejection();
		double high = get_high_rejection();
		if (readfits(filename, &ref, NULL, gfit.type == DATA_FLOAT)) {
			g_free(filename);
			return;
		}
		g_free(filename);
		set_cursor_waiting(TRUE);
		undo_save_state(&gfit, "Linear Match");
		if (!find_linear_coeff(&gfit, &ref, low, high, a, b, &error)) {
			apply_linear_to_fits(&gfit, a, b);

			adjust_cutoff_from_updated_gfit();
			redraw(com.cvport, REMAP_ALL);
			redraw_previews();
		} else {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Cannot compute linear coefficients."),
					error);
		}
		clearfits(&ref);
		set_cursor_waiting(FALSE);
	}
}
