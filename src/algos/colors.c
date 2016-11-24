/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2016 team free-astro (see more in AUTHORS file)
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
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include "core/siril.h"
#include "core/processing.h"
#include "gui/callbacks.h"
#include "io/single_image.h"
#include "gui/histogram.h"
#include "core/proto.h"
#include "algos/colors.h"
#include "core/undo.h"

/*
 * A Fast HSL-to-RGB Transform
 * by Ken Fishkin
 * from "Graphics Gems", Academic Press, 1990
 * */
/*
 *  * given h,s,l on [0..1],
 *   * return r,g,b on [0..1]
 *    */
void hsl_to_rgb(double h, double sl, double l, double * r, double * g,
		double * b) {
	double v;

	assert(h >= 0.0 && h <= 1.0);
	if (h >= 1.0) h -= 1.0;		// this code doesn't work for h = 1
	v = (l <= 0.5) ? (l * (1.0 + sl)) : (l + sl - l * sl);
	if (v <= 0) {
		*r = *g = *b = 0.0;
	} else {
		double m;
		double sv;
		int sextant;
		double fract, vsf, mid1, mid2;

		m = l + l - v;
		sv = (v - m) / v;
		h *= 6.0;
		sextant = h;
		fract = h - sextant;
		vsf = v * sv * fract;
		mid1 = m + vsf;
		mid2 = v - vsf;
		switch (sextant) {
		case 0:
			*r = v;
			*g = mid1;
			*b = m;
			break;
		case 1:
			*r = mid2;
			*g = v;
			*b = m;
			break;
		case 2:
			*r = m;
			*g = v;
			*b = mid1;
			break;
		case 3:
			*r = m;
			*g = mid2;
			*b = v;
			break;
		case 4:
			*r = mid1;
			*g = m;
			*b = v;
			break;
		case 5:
			*r = v;
			*g = m;
			*b = mid2;
			break;
		}
	}
}
/*
 *  * RGB-HSL transforms.
 *   * Ken Fishkin, Pixar Inc., January 1989.
 *    */

/*
 *  * given r,g,b on [0 ... 1],
 *   * return (h,s,l) on [0 ... 1]
 *    */
void rgb_to_hsl(double r, double g, double b, double *h, double *s, double *l) {
	double v;
	double m;
	double vm;
	double r2, g2, b2;

	v = MAX(r, g);
	v = MAX(v, b);
	m = MIN(r, g);
	m = MIN(m, b);
	*h = 0.0;
	*s = 0.0;	// init values

	if ((*l = (m + v) / 2.0) <= 0.0) {
		*l = 0.0;
		return;
	}
	if ((*s = vm = v - m) > 0.0) {
		*s /= (*l <= 0.5) ? (v + m) : (2.0 - v - m);
	} else
		return;

	r2 = (v - r) / vm;
	g2 = (v - g) / vm;
	b2 = (v - b) / vm;

	if (r == v)
		*h = (g == m ? 5.0 + b2 : 1.0 - g2);
	else if (g == v)
		*h = (b == m ? 1.0 + r2 : 3.0 - b2);
	else
		*h = (r == m ? 3.0 + g2 : 5.0 - r2);

	*h /= 6;
}

/* all variables are between 0 and 1. h takes 0 for grey */
void rgb_to_hsv(double r, double g, double b, double *h, double *s, double *v) {
	double cmax, cmin, delta;

	cmax = max(r, g);
	cmax = max(cmax, b);
	cmin = min(r, g);
	cmin = min(cmin, b);
	delta = cmax - cmin;
	*v = cmax;
	if (delta == 0.0) {
		*s = 0.0;
		*h = 0.0;
		return;
	}
	*s = delta / cmax;

	if (cmax == r)
		*h = (((g - b) / delta)) / 6.0;
	else if (cmax == g)
		*h = (((b - r) / delta) + 2.0) / 6.0;
	else
		*h = (((r - g) / delta) + 4.0) / 6.0;

	if (*h < 0.0)
		*h += 1.0;
}

void hsv_to_rgb(double h, double s, double v, double *r, double *g, double *b) {
	double p, q, t, f;
	int i;

	if (h >= 1.0)
		h -= 1.0;
	h *= 6.0;
	i = (int)h;
	f = h - (double)i;
	p = v * (1.0 - s);
	q = v * (1.0 - (s * f));
	t = v * (1.0 - (s * (1.0 - f)));

	switch (i) {
	case 0:
		*r = v;
		*g = t;
		*b = p;
		break;
	case 1:
		*r = q;
		*g = v;
		*b = p;
		break;
	case 2:
		*r = p;
		*g = v;
		*b = t;
		break;
	case 3:
		*r = p;
		*g = q;
		*b = v;
		break;
	case 4:
		*r = t;
		*g = p;
		*b = v;
		break;
	case 5:
	default:
		*r = v;
		*g = p;
		*b = q;
		break;
	}
}

void rgb_to_xyz(double r, double g, double b, double *x, double *y, double *z) {
	r = (r <= 0.04045) ? r / 12.92 : pow(((r + 0.055) / 1.055), 2.4);
	g = (g <= 0.04045) ? g / 12.92 : pow(((g + 0.055) / 1.055), 2.4);
	b = (b <= 0.04045) ? b / 12.92 : pow(((b + 0.055) / 1.055), 2.4);

	r *= 100;
	g *= 100;
	b *= 100;

	*x = 0.412453 * r + 0.357580 * g + 0.180423 * b;
	*y = 0.212671 * r + 0.715160 * g + 0.072169 * b;
	*z = 0.019334 * r + 0.119193 * g + 0.950227 * b;
}

void xyz_to_LAB(double x, double y, double z, double *L, double *a, double *b) {
	x /= 95.047;
	y /= 100.000;
	z /= 108.883;

	x = (x > 0.008856452) ? pow(x, 1 / 3.0) : (7.787037037 * x) + (16. / 116.);
	y = (y > 0.008856452) ? pow(y, 1 / 3.0) : (7.787037037 * y) + (16. / 116.);
	z = (z > 0.008856452) ? pow(z, 1 / 3.0) : (7.787037037 * z) + (16. / 116.);

	*L = (116.0 * y) - 16.0;
	*a = 500.0 * (x - y);
	*b = 200.0 * (y - z);
}

void LAB_to_xyz(double L, double a, double b, double *x, double *y, double *z) {
	*y = (L + 16.0) / 116.0;
	*x = a / 500.0 + (*y);
	*z = *y - b / 200.0;

	double x3, y3, z3;
	x3 = (*x) * (*x) * (*x);
	y3 = (*y) * (*y) * (*y);
	z3 = (*z) * (*z) * (*z);

	*x = (x3 > 0.008856452) ? x3 : (*x - 16. / 116.) / 7.787037037;
	*y = (y3 > 0.008856452) ? y3 : (*y - 16. / 116.) / 7.787037037;
	*z = (z3 > 0.008856452) ? z3 : (*z - 16. / 116.) / 7.787037037;

	*x *= 95.047;
	*y *= 100.000;
	*z *= 108.883;
}

void xyz_to_rgb(double x, double y, double z, double *r, double *g, double *b) {
	x /= 100.0;
	y /= 100.0;
	z /= 100.0;

	*r =  3.240479 * x - 1.537150 * y - 0.498535 * z;
	*g = -0.969256 * x + 1.875992 * y + 0.041556 * z;
	*b =  0.055648 * x - 0.204043 * y + 1.057311 * z;

	*r = (*r > 0.0031308) ? 1.055 * (pow(*r, (1 / 2.4))) - 0.055 : 12.92 * (*r);
	*g = (*g > 0.0031308) ? 1.055 * (pow(*g, (1 / 2.4))) - 0.055 : 12.92 * (*g);
	*b = (*b > 0.0031308) ? 1.055 * (pow(*b, (1 / 2.4))) - 0.055 : 12.92 * (*b);
}

// idle function executed at the end of the extract_channels processing
gboolean end_extract_channels(gpointer p) {
	struct extract_channels_data *args = (struct extract_channels_data *) p;
	if (args->process) {
		int i;

		stop_processing_thread();
		for (i = 0; i < 3; i++)
			save1fits16(args->channel[i], args->fit, i);
	}
	clearfits(args->fit);
	free(args);
	set_cursor_waiting(FALSE);
	update_used_memory();

	return FALSE;
}

gpointer extract_channels(gpointer p) {
	struct extract_channels_data *args = (struct extract_channels_data *) p;
	WORD *buf[3] = { args->fit->pdata[RLAYER], args->fit->pdata[GLAYER],
			args->fit->pdata[BLAYER] };
	int i;
	struct timeval t_start, t_end;
	args->process = TRUE;

	if (args->fit->naxes[2] != 3) {
		siril_log_message(
				_("Siril cannot extract layers. Make sure your image is in RGB mode.\n"));
		args->process = FALSE;
		gdk_threads_add_idle(end_extract_channels, args);
		return GINT_TO_POINTER(1);
	}

	siril_log_color_message(_("%s channel extraction: processing...\n"), "red",
			args->str_type);
	gettimeofday(&t_start, NULL);

	switch (args->type) {
	/* RGB space: nothing to do */
	case 0:
		break;
		/* HSL space */
	case 1:
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
#endif
		for (i = 0; i < args->fit->rx * args->fit->ry; i++) {
			double h, s, l;
			double r = (double) buf[RLAYER][i] / USHRT_MAX_DOUBLE;
			double g = (double) buf[GLAYER][i] / USHRT_MAX_DOUBLE;
			double b = (double) buf[BLAYER][i] / USHRT_MAX_DOUBLE;
			rgb_to_hsl(r, g, b, &h, &s, &l);
			buf[RLAYER][i] = round_to_WORD(h * 360.0);
			buf[GLAYER][i] = round_to_WORD(s * USHRT_MAX_DOUBLE);
			buf[BLAYER][i] = round_to_WORD(l * USHRT_MAX_DOUBLE);
		}
		break;
		/* HSV space */
	case 2:
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
#endif
		for (i = 0; i < args->fit->rx * args->fit->ry; i++) {
			double h, s, v;
			double r = (double) buf[RLAYER][i] / USHRT_MAX_DOUBLE;
			double g = (double) buf[GLAYER][i] / USHRT_MAX_DOUBLE;
			double b = (double) buf[BLAYER][i] / USHRT_MAX_DOUBLE;
			rgb_to_hsv(r, g, b, &h, &s, &v);
			buf[RLAYER][i] = round_to_WORD(h * 360.0);
			buf[GLAYER][i] = round_to_WORD(s * USHRT_MAX_DOUBLE);
			buf[BLAYER][i] = round_to_WORD(v * USHRT_MAX_DOUBLE);
		}
		break;
		/* CIE L*a*b */
	case 3:
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
#endif
		for (i = 0; i < args->fit->rx * args->fit->ry; i++) {
			double x, y, z, L, a, b;
			double red = (double) buf[RLAYER][i] / USHRT_MAX_DOUBLE;
			double green = (double) buf[GLAYER][i] / USHRT_MAX_DOUBLE;
			double blue = (double) buf[BLAYER][i] / USHRT_MAX_DOUBLE;
			rgb_to_xyz(red, green, blue, &x, &y, &z);
			xyz_to_LAB(x, y, z, &L, &a, &b);
			buf[RLAYER][i] = round_to_WORD(L / 100. * USHRT_MAX_DOUBLE);// 0 < L < 100
			buf[GLAYER][i] = round_to_WORD(
					((a + 128) / 255.) * USHRT_MAX_DOUBLE);	// -128 < a < 127
			buf[BLAYER][i] = round_to_WORD(
					((b + 128) / 255.) * USHRT_MAX_DOUBLE);	// -128 < b < 127
		}

	}
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	gdk_threads_add_idle(end_extract_channels, args);

	return GINT_TO_POINTER(0);
}

// idle function executed at the end of the enhance_saturation processing
gboolean end_enhance_saturation(gpointer p) {
	struct enhance_saturation_data *args = (struct enhance_saturation_data *) p;
	stop_processing_thread();
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	free(args);
	set_cursor_waiting(FALSE);
	update_used_memory();

	return FALSE;
}

gpointer enhance_saturation(gpointer p) {
	struct enhance_saturation_data *args = (struct enhance_saturation_data *) p;
	struct timeval t_start, t_end;
	double bg = 0;
	int i;

	if (!isrgb(args->fit)) {
		gdk_threads_add_idle(end_enhance_saturation, args);
		return GINT_TO_POINTER(1);
	}
	if (args->coeff == 0.0) {
		gdk_threads_add_idle(end_enhance_saturation, args);
		return GINT_TO_POINTER(1);
	}

	WORD *buf[3] = { args->fit->pdata[RLAYER], args->fit->pdata[GLAYER],
			args->fit->pdata[BLAYER] };

	siril_log_color_message(_("Saturation enhancement: processing...\n"), "red");
	gettimeofday(&t_start, NULL);

	args->h_min /= 360.0;
	args->h_max /= 360.0;
	if (args->preserve) {
		imstats *stat = statistics(args->fit, GLAYER, NULL, STATS_BASIC,
				STATS_ZERO_NULLCHECK);
		if (!stat) {
			siril_log_message(_("Error: no data computed.\n"));
			gdk_threads_add_idle(end_enhance_saturation, args);
			return GINT_TO_POINTER(1);
		}
		bg = stat->median + stat->sigma;
		bg /= stat->normValue;
		free(stat);

	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
#endif
	for (i = 0; i < args->fit->rx * args->fit->ry; i++) {
		double h, s, l;
		double r = (double) buf[RLAYER][i] / USHRT_MAX_DOUBLE;
		double g = (double) buf[GLAYER][i] / USHRT_MAX_DOUBLE;
		double b = (double) buf[BLAYER][i] / USHRT_MAX_DOUBLE;
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
		buf[RLAYER][i] = round_to_WORD(r * USHRT_MAX_DOUBLE);
		buf[GLAYER][i] = round_to_WORD(g * USHRT_MAX_DOUBLE);
		buf[BLAYER][i] = round_to_WORD(b * USHRT_MAX_DOUBLE);
	}
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	gdk_threads_add_idle(end_enhance_saturation, args);

	return GINT_TO_POINTER(0);
}

// idle function executed at the end of the scnr processing
gboolean end_scnr(gpointer p) {
	struct scnr_data *args = (struct scnr_data *) p;
	stop_processing_thread();
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	free(args);
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

	gettimeofday(&t_end, NULL);

	show_time(t_start, t_end);
	gdk_threads_add_idle(end_scnr, args);
	return GINT_TO_POINTER(0);
}

/****************** Color calibration ************************/
void on_button_bkg_selection_clicked(GtkButton *button, gpointer user_data) {
	static GtkSpinButton *selection_black_value[4] = { NULL, NULL, NULL, NULL };

	if ((!com.selection.h) || (!com.selection.w)) {
		show_dialog(_("Make a selection of the background area before"), "Warning",
				"gtk-dialog-warning");
		return;
	}

	if (!selection_black_value[0]) {
		selection_black_value[0] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_bkg_x"));
		selection_black_value[1] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_bkg_y"));
		selection_black_value[2] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_bkg_w"));
		selection_black_value[3] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_bkg_h"));
	}

	gtk_spin_button_set_value(selection_black_value[0], com.selection.x);
	gtk_spin_button_set_value(selection_black_value[1], com.selection.y);
	gtk_spin_button_set_value(selection_black_value[2], com.selection.w);
	gtk_spin_button_set_value(selection_black_value[3], com.selection.h);
}

void initialize_calibration_interface() {
	static GtkAdjustment *selection_black_adjustment[4] = { NULL, NULL, NULL,
			NULL };
	static GtkAdjustment *selection_white_adjustment[4] = { NULL, NULL, NULL,
			NULL };

	if (!selection_black_adjustment[0]) {
		selection_black_adjustment[0] = GTK_ADJUSTMENT(
				gtk_builder_get_object(builder, "adjustment_bkg_x"));
		selection_black_adjustment[1] = GTK_ADJUSTMENT(
				gtk_builder_get_object(builder, "adjustment_bkg_y"));
		selection_black_adjustment[2] = GTK_ADJUSTMENT(
				gtk_builder_get_object(builder, "adjustment_bkg_w"));
		selection_black_adjustment[3] = GTK_ADJUSTMENT(
				gtk_builder_get_object(builder, "adjustment_bkg_h"));
	}
	if (!selection_white_adjustment[0]) {
		selection_white_adjustment[0] = GTK_ADJUSTMENT(
				gtk_builder_get_object(builder, "adjustment_white_x"));
		selection_white_adjustment[1] = GTK_ADJUSTMENT(
				gtk_builder_get_object(builder, "adjustment_white_y"));
		selection_white_adjustment[2] = GTK_ADJUSTMENT(
				gtk_builder_get_object(builder, "adjustment_white_w"));
		selection_white_adjustment[3] = GTK_ADJUSTMENT(
				gtk_builder_get_object(builder, "adjustment_white_h"));
	}
	gtk_adjustment_set_upper(selection_black_adjustment[0], gfit.rx);
	gtk_adjustment_set_upper(selection_black_adjustment[1], gfit.ry);
	gtk_adjustment_set_upper(selection_black_adjustment[2], gfit.rx);
	gtk_adjustment_set_upper(selection_black_adjustment[3], gfit.ry);
	gtk_adjustment_set_value(selection_black_adjustment[0], 0);
	gtk_adjustment_set_value(selection_black_adjustment[1], 0);
	gtk_adjustment_set_value(selection_black_adjustment[2], 0);
	gtk_adjustment_set_value(selection_black_adjustment[3], 0);

	gtk_adjustment_set_upper(selection_white_adjustment[0], gfit.rx);
	gtk_adjustment_set_upper(selection_white_adjustment[1], gfit.ry);
	gtk_adjustment_set_upper(selection_white_adjustment[2], gfit.rx);
	gtk_adjustment_set_upper(selection_white_adjustment[3], gfit.ry);
	gtk_adjustment_set_value(selection_white_adjustment[0], 0);
	gtk_adjustment_set_value(selection_white_adjustment[1], 0);
	gtk_adjustment_set_value(selection_white_adjustment[2], 0);
	gtk_adjustment_set_value(selection_white_adjustment[3], 0);
}

/* This function equalize the background by giving equal value for all layers */
static void background_neutralize(fits* fit, rectangle black_selection) {
	int chan, i;
	imstats** stats;
	int ref = 0;

	assert(fit->naxes[2] == 3);

	stats = malloc(3 * sizeof(imstats *));
	for (chan = 0; chan < 3; chan++) {
		stats[chan] = statistics(fit, chan, &black_selection, STATS_BASIC,
				STATS_ZERO_NULLCHECK);
		if (!stats[chan]) {
			siril_log_message(_("Error: no data computed.\n"));
			return;
		}
		ref += stats[chan]->median;
	}
	ref /= 3;

	for (chan = 0; chan < 3; chan++) {
		int offset = stats[chan]->mean - ref;
		WORD *buf = fit->pdata[chan];
		for (i = 0; i < fit->rx * fit->ry; i++) {
			if (buf[i] < offset)
				buf[i] = 0;
			else
				buf[i] = (buf[i] - offset >= USHRT_MAX ? USHRT_MAX : buf[i] - offset);

		}
		free(stats[chan]);
	}
	free(stats);
}

void on_button_bkg_neutralization_clicked(GtkButton *button, gpointer user_data) {
	static GtkSpinButton *selection_black_value[4] = { NULL, NULL, NULL, NULL };
	rectangle black_selection;
	int width, height;

	if (!selection_black_value[0]) {
		selection_black_value[0] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_bkg_x"));
		selection_black_value[1] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_bkg_y"));
		selection_black_value[2] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_bkg_w"));
		selection_black_value[3] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_bkg_h"));
	}
	width = (int) gtk_spin_button_get_value(selection_black_value[2]);
	height = (int) gtk_spin_button_get_value(selection_black_value[3]);

	if ((!width) || (!height)) {
		show_dialog(_("Make a selection of the background area before"), "Warning",
				"gtk-dialog-warning");
		return;
	}
	black_selection.x = gtk_spin_button_get_value(selection_black_value[0]);
	black_selection.y = gtk_spin_button_get_value(selection_black_value[1]);
	black_selection.w = gtk_spin_button_get_value(selection_black_value[2]);
	black_selection.h = gtk_spin_button_get_value(selection_black_value[3]);

	undo_save_state("Processing: Background neutralization");

	set_cursor_waiting(TRUE);
	background_neutralize(&gfit, black_selection);
	delete_selected_area();

	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void on_button_white_selection_clicked(GtkButton *button, gpointer user_data) {
	static GtkSpinButton *selection_white_value[4] = { NULL, NULL, NULL, NULL };

	if (!selection_white_value[0]) {
		selection_white_value[0] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_white_x"));
		selection_white_value[1] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_white_y"));
		selection_white_value[2] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_white_w"));
		selection_white_value[3] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_white_h"));
	}

	if ((!com.selection.h) || (!com.selection.w)) {
		show_dialog(_("Make a selection of the background area before"), "Warning",
				"gtk-dialog-warning");
		return;
	}

	gtk_spin_button_set_value(selection_white_value[0], com.selection.x);
	gtk_spin_button_set_value(selection_white_value[1], com.selection.y);
	gtk_spin_button_set_value(selection_white_value[2], com.selection.w);
	gtk_spin_button_set_value(selection_white_value[3], com.selection.h);
}

static void get_coeff_for_wb(fits *fit, rectangle white, rectangle black,
		double kw[], double bg[], double *norm, double low, double high) {
	int chan, i, j, n;
	double tmp[3] = { 0.0, 0.0, 0.0 };
	WORD lo, hi;

	assert(fit->naxes[2] == 3);

	*norm = (double) get_normalized_value(fit);
	lo = round_to_WORD(low * (*norm));
	hi = round_to_WORD(high * (*norm));

	for (chan = 0; chan < 3; chan++) {
		n = 0;
		WORD *from = fit->pdata[chan] + (fit->ry - white.y - white.h) * fit->rx
				+ white.x;
		int stridefrom = fit->rx - white.w;

		for (i = 0; i < white.h; i++) {
			for (j = 0; j < white.w; j++) {
				if (*from > lo && *from < hi ) {
					kw[chan] += (double) *from / *norm;
					n++;
				}
				from++;
			}
			from += stridefrom;
		}
		if (n > 0)
			kw[chan] /= (double) n;
	}

	siril_log_message(_("Background reference:\n"));
	for (chan = 0; chan < 3; chan++) {
		imstats *stat = statistics(fit, chan, &black, STATS_BASIC, STATS_ZERO_NULLCHECK);
		if (!stat) {
			siril_log_message(_("Error: no data computed.\n"));
			return;
		}
		bg[chan] = stat->median / stat->normValue;
		siril_log_message("B%d : %.5e\n", chan, bg[chan]);
		free(stat);

	}

	siril_log_message(_("White reference:\n"));
	for (chan = 0; chan < 3; chan++) {
		siril_log_message("W%d : %.5e\n", chan, kw[chan]);
		kw[chan] = fabs(kw[chan] - bg[chan]);
	}

	int rc = (kw[0] > kw[1]) ? ((kw[0] > kw[2]) ? 0 : 2) :
					((kw[1] > kw[2]) ? 1 : 2);
	for (chan = 0; chan < 3; chan++) {
		if (chan == rc)
			tmp[chan] = 1;
		else
			tmp[chan] = kw[rc] / kw[chan];
	}

	siril_log_message(_("Color calibration factors:\n"));
	for (chan = 0; chan < 3; chan++) {
		kw[chan] = tmp[chan];
		siril_log_message("K%d : %5.3lf\n", chan, kw[chan]);
	}
}

static int calibrate(fits *fit, int layer, double kw, double bg, double norm) {
	WORD *buf;
	WORD bgNorm;
	int i;

	bgNorm = bg * norm;

	buf = fit->pdata[layer];
	for (i = 0; i < fit->rx * fit->ry; ++i) {
		buf[i] = round_to_WORD((buf[i] - bgNorm) * kw + bgNorm);
	}
	return 0;
}

static void white_balance(fits *fit, gboolean is_manual, rectangle white_selection,
		rectangle black_selection) {
	int chan;
	double norm, low, high;
	double kw[3] = { 0.0, 0.0, 0.0 };
	double bg[3] = { 0.0, 0.0, 0.0 };
	static GtkRange *scale_white_balance[3] = { NULL, NULL, NULL };
	static GtkRange *scaleLimit[2] = { NULL, NULL };

	if (scale_white_balance[RLAYER] == NULL) {
		scale_white_balance[RLAYER] = GTK_RANGE(lookup_widget("scale_r"));
		scale_white_balance[GLAYER] = GTK_RANGE(lookup_widget("scale_g"));
		scale_white_balance[BLAYER] = GTK_RANGE(lookup_widget("scale_b"));

		scaleLimit[0] = GTK_RANGE(lookup_widget("lowWhiteColorCalibScale"));
		scaleLimit[1] = GTK_RANGE(lookup_widget("upWhiteColorCalibScale"));
	}

	assert(fit->naxes[2] == 3);

	if (is_manual) {
		kw[RLAYER] = gtk_range_get_value(scale_white_balance[RLAYER]);
		kw[GLAYER] = gtk_range_get_value(scale_white_balance[GLAYER]);
		kw[BLAYER] = gtk_range_get_value(scale_white_balance[BLAYER]);
	} else {
		low = gtk_range_get_value(scaleLimit[0]);
		high = gtk_range_get_value(scaleLimit[1]);
		get_coeff_for_wb(fit, white_selection, black_selection, kw, bg, &norm, low, high);
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(chan) schedule(static)
#endif
	for (chan = 0; chan < 3; chan++) {
		if (kw[chan] == 1.0) continue;
		calibrate(fit, chan, kw[chan], bg[chan], norm);
	}
}

void on_calibration_apply_button_clicked(GtkButton *button, gpointer user_data) {
	rectangle black_selection, white_selection;
	static GtkSpinButton *selection_black_value[4] = { NULL, NULL, NULL, NULL };
	static GtkSpinButton *selection_white_value[4] = { NULL, NULL, NULL, NULL };
	struct timeval t_start, t_end;

	siril_log_color_message(_("Color Calibration: processing...\n"), "red");
	gettimeofday(&t_start, NULL);

	GtkToggleButton *manual = GTK_TOGGLE_BUTTON(
			lookup_widget("checkbutton_manual_calibration"));
	gboolean is_manual = gtk_toggle_button_get_active(manual);

	if (!selection_black_value[0]) {
		selection_black_value[0] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_bkg_x"));
		selection_black_value[1] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_bkg_y"));
		selection_black_value[2] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_bkg_w"));
		selection_black_value[3] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_bkg_h"));
	}

	if (!selection_white_value[0]) {
		selection_white_value[0] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_white_x"));
		selection_white_value[1] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_white_y"));
		selection_white_value[2] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_white_w"));
		selection_white_value[3] = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "spin_white_h"));
	}

	black_selection.x = gtk_spin_button_get_value(selection_black_value[0]);
	black_selection.y = gtk_spin_button_get_value(selection_black_value[1]);
	black_selection.w = gtk_spin_button_get_value(selection_black_value[2]);
	black_selection.h = gtk_spin_button_get_value(selection_black_value[3]);

	if (!black_selection.w || !black_selection.h) {
		show_dialog(_("Make a selection of the background area before"), "Warning",
				"gtk-dialog-warning");
		return;
	}

	white_selection.x = gtk_spin_button_get_value(selection_white_value[0]);
	white_selection.y = gtk_spin_button_get_value(selection_white_value[1]);
	white_selection.w = gtk_spin_button_get_value(selection_white_value[2]);
	white_selection.h = gtk_spin_button_get_value(selection_white_value[3]);

	if ((!white_selection.w || !white_selection.h) && !is_manual) {
		show_dialog(_("Make a selection of the background area before"), "Warning",
				"gtk-dialog-warning");
		return;
	}

	set_cursor_waiting(TRUE);
	undo_save_state("Processing: Color Calibration");
	white_balance(&gfit, is_manual, white_selection, black_selection);

	gettimeofday(&t_end, NULL);

	show_time(t_start, t_end);

	delete_selected_area();

	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void on_calibration_close_button_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("color_calibration"));
}

void on_checkbutton_manual_calibration_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	GtkWidget *manual_components = lookup_widget("grid25");
	gtk_widget_set_sensitive(manual_components,
			gtk_toggle_button_get_active(togglebutton));
}
