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
#include "core/arithm.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/OS_utils.h"
#include "algos/statistics.h"
#include "gui/image_display.h"
#include "gui/histogram.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "opencv/opencv.h"
#include "io/single_image.h"

#include "rgradient.h"

static void to_polar(int x, int y, point center, double *r, double *theta) {
	double dx = x - center.x;
	double dy = y - center.y;
	*r = sqrt(dx * dx + dy * dy);
	*theta = atan2(dy, dx);
}

static void to_cartesian(double r, double theta, point center, point *p) {
	p->x = center.x + r * cos(theta);
	p->y = center.y + r * sin(theta);
}

static gboolean end_rgradient_filter(gpointer p) {
	struct median_filter_data *args = (struct median_filter_data *) p;
	stop_processing_thread();
	adjust_vport_size_to_image();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);

	free(args);
	return FALSE;
}

gpointer rgradient_filter(gpointer p) {
	struct rgradient_filter_data *args = (struct rgradient_filter_data *) p;

	gboolean was_ushort;
	fits imA = { 0 }, imB = { 0 };
	point center = {args->xc, args->yc};
	int x, y, layer, retval = 0, cur_nb = 0;
	int w = args->fit->rx - 1;
	int h = args->fit->ry - 1;
	long n = args->fit->naxes[0] * args->fit->naxes[1] * args->fit->naxes[2];
	double dAlpha = M_PI / 180.0 * args->da;
	double total = (double)n;	// only used for progress bar
	set_progress_bar_data(_("Rotational gradient in progress..."), PROGRESS_RESET);
	was_ushort = args->fit->type == DATA_USHORT;

	/* convenient transformation to not inverse y sign */
	fits_flip_top_to_bottom(args->fit);

	retval = copyfits(args->fit, &imA, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	if (retval) { retval = 1; goto end_rgradient; }
	if (was_ushort) {
		float *newbuf = ushort_buffer_to_float(args->fit->data, n);
		if (!newbuf) { retval = 1; goto end_rgradient; }
		fit_replace_buffer(args->fit, newbuf, DATA_FLOAT);
	}

	retval = copyfits(&imA, &imB, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	if (retval) { retval = 1; goto end_rgradient; }

	retval = soper(args->fit, 2.0, OPER_MUL, TRUE);
	if (retval) { retval = 1; goto end_rgradient; }
	// args->fit will be float data after soper

	for (layer = 0; layer < args->fit->naxes[2]; layer++) {
		int i = 0;
		float *gbuf = args->fit->fpdata[layer];
		float *Abuf = imA.fpdata[layer];
		float *Bbuf = imB.fpdata[layer];
		for (y = 0; y < args->fit->ry; y++) {
			for (x = 0; x < args->fit->rx; x++) {
				double r, theta;
				point delta;

				if (!(i % 256))
					set_progress_bar_data(NULL, (double) cur_nb / total);

				to_polar(x, y, center, &r, &theta);

				// Positive differential
				to_cartesian(r - args->dR, theta + dAlpha, center, &delta);
				if (delta.x < 0)
					delta.x = fabs(delta.x);
				else if (delta.x > w)
					delta.x = 2 * w - delta.x;
				if (delta.y < 0)
					delta.y = fabs(delta.y);
				else if (delta.y > h)
					delta.y = 2 * h - delta.y;
				gbuf[i] -= Abuf[(int)delta.x + (int)delta.y * args->fit->rx];

				// Negative differential
				to_cartesian(r - args->dR, theta - dAlpha, center, &delta);
				if (delta.x < 0)
					delta.x = fabs(delta.x);
				else if (delta.x > w)
					delta.x = 2 * w - delta.x;
				if (delta.y < 0)
					delta.y = fabs(delta.y);
				else if (delta.y > h)
					delta.y = 2 * h - delta.y;
				gbuf[i] -= Bbuf[(int)delta.x + (int)delta.y * args->fit->rx];
				i++;
				cur_nb++;
			}
		}
	}

end_rgradient:
	fits_flip_top_to_bottom(args->fit);
	clearfits(&imA);
	clearfits(&imB);
	//if (was_ushort && !args->allow_32_bits) convert back to 16
	if (!retval)
		set_progress_bar_data(_("Rotational gradient complete."), PROGRESS_DONE);

	invalidate_stats_from_fit(args->fit);
	update_gfit_histogram_if_needed();
	siril_add_idle(end_rgradient_filter, args);

	return GINT_TO_POINTER(retval);
}

/// GUI

static double get_xc() {
	GtkEntry *entry;

	entry = GTK_ENTRY(lookup_widget("entry_rgradient_xc"));

	return atof(gtk_entry_get_text(entry));
}

static double get_yc() {
	GtkEntry *entry;

	entry = GTK_ENTRY(lookup_widget("entry_rgradient_yc"));

	return atof(gtk_entry_get_text(entry));
}

static double get_dR() {
	GtkRange *range;

	range = GTK_RANGE(lookup_widget("scale_radial_rgradient"));

	return gtk_range_get_value(range);
}

static double get_da() {
	GtkRange *range;

	range = GTK_RANGE(lookup_widget("scale_rot_rgradient"));

	return gtk_range_get_value(range);
}

////// CALLBACKS

void on_menuitem_rgradient_activate(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("rgradient_dialog");
}

void on_rgradient_cancel_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("rgradient_dialog");
}

void on_rgradient_Apply_clicked(GtkButton *button, gpointer user_data) {
	if (get_thread_run()) {
		siril_log_message(_("Another task is already in progress, ignoring new request.\n"));
		return;
	}

	if (!single_image_is_loaded()) return;

	struct rgradient_filter_data *args = malloc(sizeof(struct rgradient_filter_data));
	args->xc = get_xc();
	args->yc = get_yc();
	args->dR = get_dR();
	args->da = get_da();
	args->fit = &gfit;

	set_cursor_waiting(TRUE);

	undo_save_state(&gfit, "Processing: RGradient: (dR=%5.2lf, dA=%4.2lf, xc=%7.1lf, yc=%7.1lf)",
			args->dR, args->da, args->xc, args->yc);

	start_in_new_thread(rgradient_filter, args);
}
