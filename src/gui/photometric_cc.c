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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_LIBCURL

#include <math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_version.h>

#include "photometric_cc.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/undo.h"
#include "algos/statistics.h"
#include "algos/photometry.h"
#include "algos/PSF.h"
#include "algos/plateSolver.h"
#include "algos/star_finder.h"
#include "gui/message_dialog.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"

enum {
	RED, GREEN, BLUE
};


static void initialize_photometric_cc_dialog() {
	GtkWidget *button_ips_ok, *button_cc_ok, *catalog_label, *catalog_box,
			*catalog_auto, *frame_cc_bkg, *frame_cc_norm;
	GtkWindow *parent;

	button_ips_ok = lookup_widget("buttonIPS_ok");
	button_cc_ok = lookup_widget("button_cc_ok");
	catalog_label = lookup_widget("GtkLabelCatalog");
	catalog_box = lookup_widget("ComboBoxIPSCatalog");
	catalog_auto = lookup_widget("GtkCheckButton_OnlineCat");
	frame_cc_bkg = lookup_widget("frame_cc_background");
	frame_cc_norm = lookup_widget("frame_cc_norm");

	parent = GTK_WINDOW(lookup_widget("ImagePlateSolver_Dial"));

	gtk_widget_set_visible(button_ips_ok, FALSE);
	gtk_widget_set_visible(button_cc_ok, TRUE);
	gtk_widget_set_visible(catalog_label, FALSE);
	gtk_widget_set_visible(catalog_box, FALSE);
	gtk_widget_set_visible(catalog_auto, FALSE);
	gtk_widget_set_visible(frame_cc_bkg, TRUE);
	gtk_widget_set_visible(frame_cc_norm, TRUE);

	gtk_window_set_title(parent, _("Photometric Color Calibration"));
}

static void start_photometric_cc() {
	struct plate_solver_data *args = malloc(sizeof(struct plate_solver_data));
	set_cursor_waiting(TRUE);

	args->for_photometry_cc = TRUE;
	fill_plate_solver_structure(args);

	start_in_new_thread(match_catalog, args);
}

static void read_photometry_cc_file(FILE *BV_file, fitted_PSF **stars, int *nb_stars) {
	char line[512];
	fitted_PSF *star;
	int i = 0;

	while (fgets(line, 512, BV_file) != NULL) {
		int tmp;
		star = malloc(sizeof(fitted_PSF));

		sscanf(line, "%d %lf %lf %lf\n", &tmp, &(star->xpos), &(star->ypos), &(star->BV));

		stars[i] = star;
		stars[i + 1] = NULL;
		i++;
	}
	*nb_stars = i;
}

static void bv2rgb(double *r, double *g, double *b, double bv) { // RGB <0,1> <- BV <-0.4,+2.0> [-]
	double t;
	*r = 0.0;
	*g = 0.0;
	*b = 0.0;
	if (bv < -0.4)
		bv = -0.4;
	if (bv > 2.0)
		bv = 2.0;
	if ((bv >= -0.40) && (bv < 0.00)) {
		t = (bv + 0.40) / (0.00 + 0.40);
		*r = 0.61 + (0.11 * t) + (0.1 * t * t);
	} else if ((bv >= 0.00) && (bv < 0.40)) {
		t = (bv - 0.00) / (0.40 - 0.00);
		*r = 0.83 + (0.17 * t);
	} else if ((bv >= 0.40) && (bv < 2.10)) {
		t = (bv - 0.40) / (2.10 - 0.40);
		*r = 1.00;
	}
	if ((bv >= -0.40) && (bv < 0.00)) {
		t = (bv + 0.40) / (0.00 + 0.40);
		*g = 0.70 + (0.07 * t) + (0.1 * t * t);
	} else if ((bv >= 0.00) && (bv < 0.40)) {
		t = (bv - 0.00) / (0.40 - 0.00);
		*g = 0.87 + (0.11 * t);
	} else if ((bv >= 0.40) && (bv < 1.60)) {
		t = (bv - 0.40) / (1.60 - 0.40);
		*g = 0.98 - (0.16 * t);
	} else if ((bv >= 1.60) && (bv < 2.00)) {
		t = (bv - 1.60) / (2.00 - 1.60);
		*g = 0.82 - (0.5 * t * t);
	}
	if ((bv >= -0.40) && (bv < 0.40)) {
		t = (bv + 0.40) / (0.40 + 0.40);
		*b = 1.00;
	} else if ((bv >= 0.40) && (bv < 1.50)) {
		t = (bv - 0.40) / (1.50 - 0.40);
		*b = 1.00 - (0.47 * t) + (0.1 * t * t);
	} else if ((bv >= 1.50) && (bv < 1.94)) {
		t = (bv - 1.50) / (1.94 - 1.50);
		*b = 0.63 - (0.6 * t * t);
	}
}

static int make_selection_around_a_star(fitted_PSF *stars, rectangle *area, fits *fit) {
	/* make a selection around the star */
	area->x = (int) (stars->xpos - com.phot_set.outer / 2);
	area->y = (int) (stars->ypos - com.phot_set.outer / 2);
	area->w = area->h = (int) com.phot_set.outer;

	/* Don't want stars to close of the edge */
	if (area->x + area->w >= fit->rx) {
		return 1;
	}
	if (area->x - area->w <= 0) {
		return 1;
	}
	if (area->y + area->h >= fit->ry) {
		return 1;
	}
	if (area->y + area->h <= 0) {
		return 1;
	}

	return 0;
}

static void get_white_balance_coeff(fitted_PSF **stars, int nb_stars, fits *fit, double kw[], int n_f) {
	int i = 0;
	int chan;

	double *data[3];

	data[RED] = calloc(sizeof(double), nb_stars);
	data[GREEN] = calloc(sizeof(double), nb_stars);
	data[BLUE] = calloc(sizeof(double), nb_stars);

	siril_log_message(_("Applying aperture photometry to %d stars.\n"), nb_stars);

	while (stars[i]) {
		rectangle area = { 0 };
		double flux[3] = { 0.0, 0.0, 0.0 };
		double kw_tmp[3] = { 0.0, 0.0, 0.0 };
		double r, g, b, bv;

		if (make_selection_around_a_star(stars[i], &area, fit)) {
			i++;
			continue;
		}

		for (chan = 0; chan < 3; chan ++) {
			fitted_PSF *photometry = psf_get_minimisation(fit, chan, &area, TRUE, FALSE);
			if (!photometry) {
				continue;
			}
			flux[chan] = pow(10, -0.4 * photometry->mag);
			free(photometry);
		}
		/* get r g b coefficient from bv color index */
		bv = stars[i]->BV;
		bv2rgb(&r, &g, &b, bv);

		/* get white balance coeff for current star */
		kw_tmp[RED] = flux[n_f] / flux[RED] * r;
		kw_tmp[GREEN] = flux[n_f] / flux[GREEN] * g;
		kw_tmp[BLUE] = flux[n_f] / flux[BLUE] * b;
		if (isnan(kw_tmp[RED]) || (isnan(kw_tmp[GREEN]) || (isnan(kw_tmp[BLUE])))
				|| isinf(kw_tmp[RED]) || (isinf(kw_tmp[GREEN]) || (isinf(kw_tmp[BLUE])))) {
			i++;
			continue;
		}

		data[RED][i] = kw_tmp[RED];
		data[GREEN][i] = kw_tmp[GREEN];
		data[BLUE][i] = kw_tmp[BLUE];

//		printf("%d: %.3lf %.3lf %.3lf (%.2lf/%.2lf)\n", i, flux[RED], flux[GREEN], flux[BLUE], stars[i]->xpos, stars[i]->ypos);
//		printf("%d: %.3lf %.3lf %.3lf\n", i, kw_tmp[RED], kw_tmp[GREEN], kw_tmp[BLUE]);

		i++;
	}

	/* sort in ascending order before using gsl_stats_trmean_from_sorted_data */
	gsl_sort(data[RED], 1, nb_stars);
	gsl_sort(data[GREEN], 1, nb_stars);
	gsl_sort(data[BLUE], 1, nb_stars);

#if (GSL_MAJOR_VERSION == 2) && (GSL_MINOR_VERSION > 4)
	double alpha = 0.3;

	kw[RED] = gsl_stats_trmean_from_sorted_data(alpha, data[RED], 1, nb_stars);
	kw[GREEN] = gsl_stats_trmean_from_sorted_data(alpha, data[GREEN], 1, nb_stars);
	kw[BLUE] = gsl_stats_trmean_from_sorted_data(alpha, data[BLUE], 1, nb_stars);
#else
	kw[RED] = gsl_stats_median_from_sorted_data(data[RED], 1, nb_stars);
	kw[GREEN] = gsl_stats_median_from_sorted_data(data[GREEN], 1, nb_stars);
	kw[BLUE] = gsl_stats_median_from_sorted_data(data[BLUE], 1, nb_stars);
#endif

	/* normalize factors */
	kw[RED] /= (kw[n_f]);
	kw[GREEN] /= (kw[n_f]);
	kw[BLUE] /= (kw[n_f]);
	siril_log_message(_("Color calibration factors:\n"));
	for (chan = 0; chan < 3; chan++) {
		siril_log_message("K%d: %5.3lf\n", chan, kw[chan]);
	}

	free(data[RED]);
	free(data[GREEN]);
	free(data[BLUE]);
}

static void get_background_coefficients(fits *fit, rectangle *area, coeff bg[]) {
	int chan;

	siril_log_message(_("Background reference:\n"));
	for (chan = 0; chan < 3; chan++) {
		imstats *stat = statistics(NULL, -1, fit, chan, area, STATS_BASIC);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return;
		}
		bg[chan].value = stat->median / stat->normValue;
		bg[chan].channel = chan;
		siril_log_message("B%d: %.5e\n", chan, bg[chan].value);
		free_stats(stat);
	}
}

static int calibrate_colors(fits *fit, double kw[], coeff bg[], double norm) {
	WORD *buf;
	int i, chan;

	for (chan = 0; chan < 3; chan++) {
		if (kw[chan] == 1.0)
			continue;
		WORD bgNorm = bg[chan].value * norm;
		buf = fit->pdata[chan];
		for (i = 0; i < fit->rx * fit->ry; ++i) {
			buf[i] = round_to_WORD((buf[i] - bgNorm) * kw[chan] + bgNorm);
		}
	}
	return 0;
}

/* This function equalize the background by giving equal value for all layers */
static void background_neutralize(fits* fit, coeff bg[], int n_f, double norm) {
	int chan, i;

	for (chan = 0; chan < 3; chan++) {
		int offset = (bg[chan].value - bg[n_f].value) * norm;
		siril_debug_print("offset: %d, %d\n", chan, offset);
		WORD *buf = fit->pdata[chan];
		for (i = 0; i < fit->rx * fit->ry; i++) {
			buf[i] = round_to_WORD(buf[i] - offset);
		}
	}
	invalidate_stats_from_fit(fit);
}

int struct_cmp(const void *a, const void *b) {
	coeff *a1 = (coeff *) a;
	coeff *a2 = (coeff*) b;
	if ((*a1).value > (*a2).value)
		return 1;
	else if ((*a1).value < (*a2).value)
		return -1;
	else
		return 0;
}

static int determine_chan_for_norm(coeff bg[], int n_f) {
	/* make copy bg coefficients */
	coeff tmp[3];
	memcpy(tmp, bg, 3 * sizeof(coeff));
	/* ascending order */
	qsort(tmp, 3, sizeof(tmp[0]), struct_cmp);

	if (n_f == 0) { /* on highest */
		return tmp[2].channel;
	} else if (n_f == 1) { /* on middle */
		return tmp[1].channel;
	} else { /* on lowest */
		return tmp[0].channel;
	}
}

static gboolean end_photometric_cc(gpointer p) {
	struct photometric_cc_data *args = (struct photometric_cc_data *) p;
	stop_processing_thread();

	free_fitted_stars(args->stars);
	fclose(args->BV_file);
	free(args);
	update_used_memory();

	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	return FALSE;
}

static gpointer photometric_cc(gpointer p) {
	struct photometric_cc_data *args = (struct photometric_cc_data *) p;
	double kw[3], norm;
	coeff bg[3];
	int nb_stars, chan;
	rectangle *bkg_sel;

	if (!args->bg_auto) {
		bkg_sel = &(args->bg_area);
	} else {
		bkg_sel = NULL;
	}

	read_photometry_cc_file(args->BV_file, args->stars, &nb_stars);

	get_background_coefficients(&gfit, bkg_sel, bg);
	chan = determine_chan_for_norm(bg, args->n_f);
	siril_log_message(_("Normalizing on %s channel\n"), (chan == 0) ? "red" : ((chan == 1) ? "green" : "blue"));
	get_white_balance_coeff(args->stars, nb_stars, &gfit, kw, chan);
	norm = (double) get_normalized_value(&gfit);

	calibrate_colors(&gfit, kw, bg, norm);
	background_neutralize(&gfit, bg, chan, norm);

	siril_add_idle(end_photometric_cc, args);
	return GINT_TO_POINTER(0);
}

static gboolean is_selection_ok() {
	static GtkSpinButton *selection_black_value[4] = { NULL, NULL, NULL, NULL };
	int width, height;

	if (!selection_black_value[0]) {
		selection_black_value[0] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_x"));
		selection_black_value[1] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_y"));
		selection_black_value[2] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_w"));
		selection_black_value[3] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_h"));
	}
	width = (int) gtk_spin_button_get_value(selection_black_value[2]);
	height = (int) gtk_spin_button_get_value(selection_black_value[3]);

	if ((!width) || (!height)) {
		return FALSE;
	}
	return TRUE;
}

static rectangle get_bkg_selection() {
	static GtkSpinButton *selection_black_value[4] = { NULL, NULL, NULL, NULL };
	rectangle black_selection;

	if (!selection_black_value[0]) {
		selection_black_value[0] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_x"));
		selection_black_value[1] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_y"));
		selection_black_value[2] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_w"));
		selection_black_value[3] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_h"));
	}

	black_selection.x = gtk_spin_button_get_value(selection_black_value[0]);
	black_selection.y = gtk_spin_button_get_value(selection_black_value[1]);
	black_selection.w = gtk_spin_button_get_value(selection_black_value[2]);
	black_selection.h = gtk_spin_button_get_value(selection_black_value[3]);

	return black_selection;
}

/****
 * PUBLIC FUNCTIONS
 */

FILE *open_bv_file(const gchar *mode) {
	gchar *filename = g_build_filename(g_get_tmp_dir(), "photometric_cc.dat", NULL);
	FILE *BV = g_fopen(filename, mode);
	g_free(filename);
	return BV;
}

int apply_photometric_cc() {
	fitted_PSF **stars;
	FILE *BV_file = NULL;
	GtkComboBox *norm_box;
	GtkToggleButton *auto_bkg;

	norm_box = GTK_COMBO_BOX(lookup_widget("combo_box_cc_norm"));
	auto_bkg = GTK_TOGGLE_BUTTON(lookup_widget("button_cc_bkg_auto"));

	undo_save_state(&gfit, "Processing: Photometric CC");
	invalidate_stats_from_fit(&gfit);

	set_cursor_waiting(TRUE);
	BV_file = open_bv_file("r+t");

	stars = malloc((MAX_STARS + 1) * sizeof(fitted_PSF *));
	if (stars == NULL) {
		printf("Memory allocation failed: apply_photometric_cc\n");
		return 1;
	}

	struct photometric_cc_data *args = malloc(sizeof(struct photometric_cc_data));
	set_cursor_waiting(TRUE);

	args->stars = stars;
	args->BV_file = BV_file;
	args->n_f = gtk_combo_box_get_active(norm_box);
	args->bg_area = get_bkg_selection();
	args->bg_auto = gtk_toggle_button_get_active(auto_bkg);

	start_in_new_thread(photometric_cc, args);

	return 0;
}

/*****
 * CALLBACKS FUNCTIONS
 */

void on_menuitemphotometriccalibration_activate() {
	initialize_photometric_cc_dialog();
	gtk_widget_show(lookup_widget("ImagePlateSolver_Dial"));
}

void on_button_cc_ok_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton *auto_bkg;

	auto_bkg = GTK_TOGGLE_BUTTON(lookup_widget("button_cc_bkg_auto"));

	if ((!gtk_toggle_button_get_active(auto_bkg)) && (!is_selection_ok())) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the background area"));
	} else
		start_photometric_cc();
}

void on_button_cc_bkg_auto_toggled(GtkToggleButton *button,
		gpointer user_data) {
	GtkWidget *box_cc_manual_bkg;

	box_cc_manual_bkg = lookup_widget("box_cc_manual_bkg");
	gtk_widget_set_sensitive(box_cc_manual_bkg, !gtk_toggle_button_get_active(button));
}

void on_button_cc_bkg_selection_clicked(GtkButton *button, gpointer user_data) {
	static GtkSpinButton *selection_bkg_value[4] = { NULL, NULL, NULL, NULL };

	if ((!com.selection.h) || (!com.selection.w)) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the background area"));
		return;
	}

	if (!selection_bkg_value[0]) {
		selection_bkg_value[0] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_x"));
		selection_bkg_value[1] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_y"));
		selection_bkg_value[2] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_w"));
		selection_bkg_value[3] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_h"));
	}

	gtk_spin_button_set_value(selection_bkg_value[0], com.selection.x);
	gtk_spin_button_set_value(selection_bkg_value[1], com.selection.y);
	gtk_spin_button_set_value(selection_bkg_value[2], com.selection.w);
	gtk_spin_button_set_value(selection_bkg_value[3], com.selection.h);
}

#endif