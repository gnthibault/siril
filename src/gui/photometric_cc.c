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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_LIBCURL

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_statistics.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/sleef.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/OS_utils.h"
#include "algos/sorting.h"
#include "algos/statistics.h"
#include "algos/photometry.h"
#include "algos/PSF.h"
#include "algos/plateSolver.h"
#include "algos/star_finder.h"
#include "gui/image_display.h"
#include "gui/message_dialog.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "gui/histogram.h"
#include "gui/dialogs.h"

#include "photometric_cc.h"

enum {
	RED, GREEN, BLUE
};

static void initialize_photometric_cc_dialog() {
	GtkWidget *button_ips_ok, *button_cc_ok, *catalog_label, *catalog_box_ips,
			*catalog_box_pcc, *catalog_auto, *frame_cc_bkg, *frame_cc_norm,
			*catalog_label_pcc;
	GtkWindow *parent;
	GtkAdjustment *selection_cc_black_adjustment[4];

	button_ips_ok = lookup_widget("buttonIPS_ok");
	button_cc_ok = lookup_widget("button_cc_ok");
	catalog_label = lookup_widget("GtkLabelCatalog");
	catalog_label_pcc = lookup_widget("GtkLabelCatalogPCC");
	catalog_box_ips = lookup_widget("ComboBoxIPSCatalog");
	catalog_box_pcc = lookup_widget("ComboBoxPCCCatalog");
	catalog_auto = lookup_widget("GtkCheckButton_OnlineCat");
	frame_cc_bkg = lookup_widget("frame_cc_background");
	frame_cc_norm = lookup_widget("frame_cc_norm");

	parent = GTK_WINDOW(lookup_widget("ImagePlateSolver_Dial"));

	selection_cc_black_adjustment[0] = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment_cc_bkg_x"));
	selection_cc_black_adjustment[1] = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment_cc_bkg_y"));
	selection_cc_black_adjustment[2] = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment_cc_bkg_w"));
	selection_cc_black_adjustment[3] = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment_cc_bkg_h"));

	gtk_widget_set_visible(button_ips_ok, FALSE);
	gtk_widget_set_visible(button_cc_ok, TRUE);
	gtk_widget_set_visible(catalog_label, FALSE);
	gtk_widget_set_visible(catalog_label_pcc, TRUE);
	gtk_widget_set_visible(catalog_box_ips, FALSE);
	gtk_widget_set_visible(catalog_box_pcc, TRUE);
	gtk_widget_set_visible(catalog_auto, FALSE);
	gtk_widget_set_visible(frame_cc_bkg, TRUE);
	gtk_widget_set_visible(frame_cc_norm, TRUE);

	gtk_window_set_title(parent, _("Photometric Color Calibration"));

	gtk_adjustment_set_upper(selection_cc_black_adjustment[0], gfit.rx);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[1], gfit.ry);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[2], gfit.rx);
	gtk_adjustment_set_upper(selection_cc_black_adjustment[3], gfit.ry);
	gtk_adjustment_set_value(selection_cc_black_adjustment[0], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[1], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[2], 0);
	gtk_adjustment_set_value(selection_cc_black_adjustment[3], 0);

}

static void start_photometric_cc() {
	struct plate_solver_data *args = malloc(sizeof(struct plate_solver_data));

	args->for_photometry_cc = TRUE;
	if (!fill_plate_solver_structure(args)) {
		set_cursor_waiting(TRUE);
		start_in_new_thread(match_catalog, args);
	}
}

static void read_photometry_cc_file(FILE *BV_file, fitted_PSF **stars, int *nb_stars) {
	char line[512];
	int i = 0;

	while (fgets(line, 512, BV_file) != NULL) {
		int tmp;
		fitted_PSF *star = malloc(sizeof(fitted_PSF));

		sscanf(line, "%d %lf %lf %lf\n", &tmp, &(star->xpos), &(star->ypos), &(star->BV));

		stars[i] = star;
		stars[i + 1] = NULL;
		i++;
	}
	*nb_stars = i;
}

static void bv2rgb(float *r, float *g, float *b, float bv) { // RGB <0,1> <- BV <-0.4,+2.0> [-]
	float t;
	*r = 0.f;
	*g = 0.f;
	*b = 0.f;
	if (bv < -0.4f)
		bv = -0.4f;
	if (bv > 2.f)
		bv = 2.f;
	if ((bv >= -0.4f) && (bv < 0.0f)) {
		t = (bv + 0.4f) / (0.f + 0.4f);
		*r = 0.61f + (0.11f * t) + (0.1f * t * t);
	} else if ((bv >= 0.f) && (bv < 0.4f)) {
		t = (bv - 0.0f) / (0.4f - 0.f);
		*r = 0.83f + (0.17f * t);
	} else if ((bv >= 0.4f) && (bv < 2.1f)) {
		*r = 1.f;
	}
	if ((bv >= -0.4f) && (bv < 0.f)) {
		t = (bv + 0.4f) / (0.f + 0.4f);
		*g = 0.7f + (0.07f * t) + (0.1f * t * t);
	} else if ((bv >= 0.f) && (bv < 0.4f)) {
		t = (bv - 0.f) / (0.4f - 0.f);
		*g = 0.87f + (0.11f * t);
	} else if ((bv >= 0.4f) && (bv < 1.6f)) {
		t = (bv - 0.4f) / (1.6f - 0.4f);
		*g = 0.98f - (0.16f * t);
	} else if ((bv >= 1.6f) && (bv < 2.f)) {
		t = (bv - 1.6f) / (2.f - 1.6f);
		*g = 0.82f - (0.5f * t * t);
	}
	if ((bv >= -0.4f) && (bv < 0.4f)) {
		*b = 1.f;
	} else if ((bv >= 0.4f) && (bv < 1.5f)) {
		t = (bv - 0.4f) / (1.5f - 0.4f);
		*b = 1.f - (0.47f * t) + (0.1f * t * t);
	} else if ((bv >= 1.5f) && (bv < 1.94f)) {
		t = (bv - 1.5f) / (1.94f - 1.5f);
		*b = 0.63f - (0.6f * t * t);
	}
}

static int make_selection_around_a_star(fitted_PSF *stars, rectangle *area, fits *fit) {
	/* make a selection around the star */
	area->x = round_to_int(stars->xpos - com.pref.phot_set.outer);
	area->y = round_to_int(stars->ypos - com.pref.phot_set.outer);
	area->w = area->h = round_to_int(com.pref.phot_set.outer * 2);

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
	if (area->y - area->h <= 0) {
		return 1;
	}

	return 0;
}

static float Qn0(const float sorted_data[], const size_t stride, const size_t n) {
	const size_t wsize = n * (n - 1) / 2;
	const size_t n_2 = n / 2;
	const size_t k = ((n_2 + 1) * n_2) / 2;
	size_t idx = 0;

	if (n < 2)
		return (0.0);

	float *work = malloc(wsize * sizeof(float));
	if (!work) {
		PRINT_ALLOC_ERR;
		return -1.0f;
	}

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j)
			work[idx++] = fabsf(sorted_data[i] - sorted_data[j]);
	}

	quickmedian_float(work, idx);
	float Qn = work[k - 1];

	free(work);
	return Qn;
}

static float siril_stats_trmean_from_sorted_data(const float trim,
		const float sorted_data[], const size_t stride, const size_t size) {
	if (trim >= 0.5f) {
		return (float) gsl_stats_float_median_from_sorted_data(sorted_data, stride, size);
	} else {
		size_t ilow = (size_t) floorf(trim * size);
		size_t ihigh = size - ilow - 1;
		float mean = 0.f;
		float k = 0.f;

		/* compute mean of middle samples in [ilow,ihigh] */
		for (size_t i = ilow; i <= ihigh; ++i) {
			float delta = sorted_data[i * stride] - mean;
			k += 1.f;
			mean += delta / k;
		}
		return mean;
	}
}

static float siril_stats_robust_mean(const float sorted_data[],
		const size_t stride, const size_t size) {
	float mx = (float) gsl_stats_float_median_from_sorted_data(sorted_data, stride, size);
	float qn0 = Qn0(sorted_data, 1, size);
	if (qn0 < 0)
		return -1.0f;
	float sx = 2.2219f * qn0;
	float *x, mean;
	int i, j;

	x = malloc(size * sizeof(float));
	if (!x) {
		PRINT_ALLOC_ERR;
		return -1.0f;
	}

	for (i = 0, j = 0; i < size; ++i) {
		if (fabsf(sorted_data[i] - (float) mx) < 3 * (float) sx) {
			x[j++] = sorted_data[i];
		}
	}
	/* not enough stars, try something anyway */
	if (j < 5) {
		mean = (float) siril_stats_trmean_from_sorted_data(0.3f, sorted_data, stride,
				size);
	} else {
		mean = (float) gsl_stats_float_mean(x, stride, j);
	}
	free(x);
	return mean;
}

static int get_white_balance_coeff(fitted_PSF **stars, int nb_stars, fits *fit, float kw[], int n_channel) {
	int i = 0, ngood = 0;
	gboolean no_phot = FALSE;
	int progress = 0;
	float *data[3];

	data[RED] = malloc(sizeof(float) * nb_stars);
	data[GREEN] = malloc(sizeof(float) * nb_stars);
	data[BLUE] = malloc(sizeof(float) * nb_stars);
	if (!data[RED] || !data[GREEN] || !data[BLUE]) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	/* initialize to DBL_MAX */
	for (int k = 0; k < nb_stars; k++) {
		data[RED][k] = FLT_MAX;
		data[GREEN][k] = FLT_MAX;
		data[BLUE][k] = FLT_MAX;
	}

	siril_log_message(_("Applying aperture photometry to %d stars.\n"), nb_stars);
	set_progress_bar_data(_("Photometry color calibration in progress..."), PROGRESS_RESET);

	while (stars[i]) {
		rectangle area = { 0 };
		float flux[3] = { 0.f, 0.f, 0.f };
		float r, g, b, bv;
		if (!(i % 16))	// every 16 iterations
			set_progress_bar_data(NULL, (double) progress / (double) nb_stars);
		progress++;

		if (make_selection_around_a_star(stars[i], &area, fit)) {
			i++;
			continue;
		}

		for (int chan = 0; chan < 3; chan ++) {
			fitted_PSF *photometry = psf_get_minimisation(fit, chan, &area, TRUE, FALSE, TRUE);
			if (!photometry || !photometry->phot_is_valid) {
				no_phot = TRUE;
				break;
			}
			flux[chan] = powf(10.f, -0.4f * (float) photometry->mag);
			free(photometry);
		}
		if (no_phot) {
			i++;
			no_phot = FALSE;
			continue;
		}
		/* get r g b coefficient from bv color index */
		bv = stars[i]->BV;
		bv2rgb(&r, &g, &b, bv);

		/* get Color calibration factors for current star */
		data[RED][i] = (flux[n_channel] / flux[RED]) * r;
		data[GREEN][i] = (flux[n_channel] / flux[GREEN]) * g;
		data[BLUE][i] = (flux[n_channel] / flux[BLUE]) * b;

		if (xisnanf(data[RED][i]) || xisnanf(data[GREEN][i]) || xisnanf(data[BLUE][i])) {
			data[RED][i] = FLT_MAX;
			data[GREEN][i] = FLT_MAX;
			data[BLUE][i] = FLT_MAX;
			i++;
			continue;
		}
		i++;
		ngood++;
	}
	int excl = nb_stars - ngood;
	gchar *str = ngettext("%d star excluded from the calculation\n", "%d stars excluded from the calculation\n", excl);
	str = g_strdup_printf(str, excl);
	siril_log_message(str);
	g_free(str);

	if (ngood == 0) {
		siril_log_message(_("No valid stars found.\n"));
		free(data[RED]);
		free(data[GREEN]);
		free(data[BLUE]);
		return 1;
	}
	/* sort in ascending order before using siril_stats_mean_from_linearFit
	 Hence, DBL_MAX are at the end of the tab */
	quickmedian_float(data[RED], nb_stars);
	quickmedian_float(data[GREEN], nb_stars);
	quickmedian_float(data[BLUE], nb_stars);

	/* we do not take into account DBL_MAX values */
	kw[RED] = siril_stats_robust_mean(data[RED], 1, ngood);
	kw[GREEN] = siril_stats_robust_mean(data[GREEN], 1, ngood);
	kw[BLUE] = siril_stats_robust_mean(data[BLUE], 1, ngood);
	if (kw[RED] < 0.f || kw[GREEN] < 0.f || kw[BLUE] < 0.f) {
		free(data[RED]);
		free(data[GREEN]);
		free(data[BLUE]);
		return 1;
	}

	/* normalize factors */
	kw[RED] /= (kw[n_channel]);
	kw[GREEN] /= (kw[n_channel]);
	kw[BLUE] /= (kw[n_channel]);
	siril_log_message(_("Color calibration factors:\n"));
	for (int chan = 0; chan < 3; chan++) {
		siril_log_message("K%d: %5.3lf\n", chan, kw[chan]);
	}

	free(data[RED]);
	free(data[GREEN]);
	free(data[BLUE]);

	return 0;
}

static void get_background_coefficients(fits *fit, rectangle *area, coeff bg[], gboolean verbose) {

	if (verbose) siril_log_message(_("Background reference:\n"));
	for (int chan = 0; chan < 3; chan++) {
		imstats *stat = statistics(NULL, -1, fit, chan, area, STATS_BASIC, TRUE);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return;
		}
		bg[chan].value = stat->median / stat->normValue;
		bg[chan].channel = chan;
		if (verbose) siril_log_message("B%d: %.5e\n", chan, bg[chan].value);
		free_stats(stat);
	}
}

static int apply_white_balance(fits *fit, float kw[]) {
	for (int chan = 0; chan < 3; chan++) {
		float scale = kw[chan];
		if (scale == 1.0) continue;

		size_t i, n = fit->naxes[0] * fit->naxes[1];
		if (fit->type == DATA_USHORT) {
			WORD *buf = fit->pdata[chan];
			for (i = 0; i < n; ++i) {
				buf[i] = roundf_to_WORD((float)buf[i] * scale);
			}
		}
		else if (fit->type == DATA_FLOAT) {
			float *buf = fit->fpdata[chan];
			for (i = 0; i < n; ++i) {
				buf[i] = buf[i] * scale;
			}
		}
		else return 1;
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

/* This function equalize the background by giving equal value for all layers */
static void background_neutralize(fits* fit, coeff bg[], int n_channel, double norm) {
	size_t i, n = fit->naxes[0] * fit->naxes[1];
	if (fit->type == DATA_USHORT) {
		for (int chan = 0; chan < 3; chan++) {
			float offset = (bg[chan].value - bg[n_channel].value) * norm;
			siril_debug_print("offset: %d, %f\n", chan, offset);
			WORD *buf = fit->pdata[chan];
			for (i = 0; i < n; ++i) {
				buf[i] = roundf_to_WORD((float)buf[i] - offset);
			}
		}
	}
	else if (fit->type == DATA_FLOAT) {
		for (int chan = 0; chan < 3; chan++) {
			float offset = bg[chan].value - bg[n_channel].value;
			siril_debug_print("offset: %d, %f\n", chan, offset);
			float *buf = fit->fpdata[chan];
			for (i = 0; i < n; ++i) {
				buf[i] = buf[i] - offset;
			}
		}
	}
	invalidate_stats_from_fit(fit);
}

static int cmp_coeff(const void *a, const void *b) {
	coeff *a1 = (coeff *) a;
	coeff *a2 = (coeff *) b;
	if ((*a1).value > (*a2).value)
		return 1;
	else if ((*a1).value < (*a2).value)
		return -1;
	else
		return 0;
}

static int determine_chan_for_norm(coeff bg[], int n_channel) {
	/* make a copy of bg coefficients because we don't
	 * want to sort original data */
	coeff tmp[3];
	memcpy(tmp, bg, 3 * sizeof(coeff));
	/* ascending order */
	qsort(tmp, 3, sizeof(tmp[0]), cmp_coeff);

	if (n_channel == 0) { /* on highest */
		return tmp[2].channel;
	} else if (n_channel == 1) { /* on middle */
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

	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	
	set_cursor_waiting(FALSE);
	return FALSE;
}

static gpointer photometric_cc(gpointer p) {
	struct photometric_cc_data *args = (struct photometric_cc_data *) p;
	float kw[3];
	coeff bg[3];
	int nb_stars, chan;
	rectangle *bkg_sel;

	if (!args->bg_auto) {
		bkg_sel = &(args->bg_area);
	} else {
		bkg_sel = NULL;
	}

	/* make sure parameters are initialized in a good way */
	if (com.pref.phot_set.outer == 0.0) {
		initialize_photometric_param();
	}

	read_photometry_cc_file(args->BV_file, args->stars, &nb_stars);

	get_background_coefficients(&gfit, bkg_sel, bg, FALSE);
	chan = determine_chan_for_norm(bg, args->n_channel);
	siril_log_message(_("Normalizing on %s channel.\n"), (chan == 0) ? _("red") : ((chan == 1) ? _("green") : _("blue")));
	int ret = get_white_balance_coeff(args->stars, nb_stars, &gfit, kw, chan);
	if (!ret) {
		double norm = get_normalized_value(&gfit);
		apply_white_balance(&gfit, kw);
		get_background_coefficients(&gfit, bkg_sel, bg, TRUE);
		background_neutralize(&gfit, bg, chan, norm);
		set_progress_bar_data(_("Photometric Color Calibration applied"), PROGRESS_DONE);
	} else {
		set_progress_bar_data(_("Photometric Color Calibration failed"), PROGRESS_DONE);
	}

	siril_add_idle(end_photometric_cc, args);
	return GINT_TO_POINTER(ret);
}

static gboolean is_selection_ok() {
	static GtkSpinButton *selection_black_value[4] = { NULL, NULL, NULL, NULL };

	if (!selection_black_value[0]) {
		selection_black_value[0] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_x"));
		selection_black_value[1] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_y"));
		selection_black_value[2] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_w"));
		selection_black_value[3] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_h"));
	}
	int width = (int) gtk_spin_button_get_value(selection_black_value[2]);
	int height = (int) gtk_spin_button_get_value(selection_black_value[3]);

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

	undo_save_state(&gfit, _("Photometric CC"));
	invalidate_stats_from_fit(&gfit);
	invalidate_gfit_histogram();

	set_cursor_waiting(TRUE);
	BV_file = open_bv_file("r+t");

	stars = malloc((MAX_STARS + 1) * sizeof(fitted_PSF *));
	if (stars == NULL) {
		PRINT_ALLOC_ERR;
		set_cursor_waiting(FALSE);
		fclose(BV_file);
		return 1;
	}

	struct photometric_cc_data *args = malloc(sizeof(struct photometric_cc_data));

	args->stars = stars;
	args->BV_file = BV_file;
	args->n_channel = gtk_combo_box_get_active(norm_box);
	args->bg_area = get_bkg_selection();
	args->bg_auto = gtk_toggle_button_get_active(auto_bkg);

	start_in_new_thread(photometric_cc, args);

	return 0;
}

int get_photometry_catalog() {
	GtkComboBox *box;
	int ret;

	box = GTK_COMBO_BOX(lookup_widget("ComboBoxPCCCatalog"));
	ret = gtk_combo_box_get_active(box);

	if (ret == 1) {
		return APASS;
	} else {
		return NOMAD;
	}
}

/*****
 * CALLBACKS FUNCTIONS
 */

void on_menuitemphotometriccalibration_activate() {
	initialize_photometric_cc_dialog();
	siril_open_dialog("ImagePlateSolver_Dial");
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
	static GtkSpinButton *selection_cc_bkg_value[4] = { NULL, NULL, NULL, NULL };

	if ((!com.selection.h) || (!com.selection.w)) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("There is no selection"),
				_("Make a selection of the background area"));
		return;
	}

	if (!selection_cc_bkg_value[0]) {
		selection_cc_bkg_value[0] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_x"));
		selection_cc_bkg_value[1] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_y"));
		selection_cc_bkg_value[2] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_w"));
		selection_cc_bkg_value[3] = GTK_SPIN_BUTTON(lookup_widget("spin_cc_bkg_h"));
	}

	gtk_spin_button_set_value(selection_cc_bkg_value[0], com.selection.x);
	gtk_spin_button_set_value(selection_cc_bkg_value[1], com.selection.y);
	gtk_spin_button_set_value(selection_cc_bkg_value[2], com.selection.w);
	gtk_spin_button_set_value(selection_cc_bkg_value[3], com.selection.h);
}

#endif
