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
#include "core/processing.h"
#include "core/undo.h"
#include "algos/statistics.h"
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
	stop_processing_thread();// can it be done here in case there is no thread?
	adjust_vport_size_to_image();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	update_used_memory();
	free(args);
	return FALSE;
}

gpointer rgradient_filter(gpointer p) {
	struct rgradient_filter_data *args = (struct rgradient_filter_data *) p;

	fits imA = { 0 };
	fits imB = { 0 };
	point center = {args->xc, args->yc};
	int x, y, layer, cur_nb = 0;
	int w = args->fit->rx - 1;
	int h = args->fit->ry - 1;
	double dAlpha = M_PI / 180.0 * args->da;

	double total = (double)(args->fit->rx * args->fit->ry * args->fit->naxes[2]);	// only used for progress bar
	set_progress_bar_data(_("Rotational gradient in progress..."), PROGRESS_RESET);

	/* convenient transformation to not inverse y sign */
	fits_flip_top_to_bottom(args->fit);

	copyfits(args->fit, &imA, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	copyfits(args->fit, &imB, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);

	soper(args->fit, 2.0, OPER_MUL);

	for (layer = 0; layer < args->fit->naxes[2]; layer++) {
		int i = 0;
		for (y = 0; y < args->fit->ry; y++) {
			for (x = 0; x < args->fit->rx; x++) {
				double r, theta;
				point delta;
				WORD *gbuf = args->fit->pdata[layer];
				WORD *Abuf = imA.pdata[layer];
				WORD *Bbuf = imB.pdata[layer];

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

	fits_flip_top_to_bottom(args->fit);
	set_progress_bar_data(_("Rotational gradient complete."), PROGRESS_DONE);

	clearfits(&imA);
	clearfits(&imB);
	invalidate_stats_from_fit(args->fit);
	update_gfit_histogram_if_needed();
	siril_add_idle(end_rgradient_filter, args);

	return GINT_TO_POINTER(0);
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
