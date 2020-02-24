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
#include "core/sleef.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/OS_utils.h"
#include "algos/statistics.h"
#include "algos/PSF.h"
#include "gui/image_display.h"
#include "gui/histogram.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "opencv/opencv.h"
#include "io/single_image.h"

#include "rgradient.h"

static void to_polar(int x, int y, point center, double *r, double *theta) {
	double dx = x - center.x;
	double dy = y - center.y;
	*r = sqrt(dx * dx + dy * dy);
	*theta = xatan2(dy, dx);
}

static void to_cartesian(double r, double theta, point center, point *p) {
    const double2 sincosval = xsincos(theta);
	p->x = center.x + r * sincosval.y;
	p->y = center.y + r * sincosval.x;
}

static gboolean end_rgradient_filter(gpointer p) {
	struct median_filter_data *args = (struct median_filter_data *) p;
	stop_processing_thread();// can it be done here in case there is no thread?
	adjust_vport_size_to_image();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	
	free(args);
	return FALSE;
}

gpointer rgradient_filter(gpointer p) {
	struct timeval t_start, t_end;

	siril_log_color_message(_("Rotational gradient: processing...\n"), "red");
	gettimeofday(&t_start, NULL);

	struct rgradient_filter_data *args = (struct rgradient_filter_data *) p;

	fits imA = { 0 };
	fits imB = { 0 };
	const point center = {args->xc, args->yc};
	const int w = args->fit->rx - 1;
	const int h = args->fit->ry - 1;
	const double dAlpha = M_PI / 180.0 * args->da;

	int cur_nb = 0;	// only used for progress bar
	const double total = args->fit->ry * args->fit->naxes[2];	// only used for progress bar
	set_progress_bar_data(_("Rotational gradient in progress..."), PROGRESS_RESET);

	/* convenient transformation to not inverse y sign */
	fits_flip_top_to_bottom(args->fit);

	copyfits(args->fit, &imA, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	copyfits(args->fit, &imB, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);

	soper(args->fit, 2.0, OPER_MUL);
    const double w2 = 2 * w;
    const double h2 = 2 * h;
    const double wd = w;
    const double hd = h;
    for (int layer = 0; layer < args->fit->naxes[2]; layer++) {
        WORD *gbuf = args->fit->pdata[layer];
        WORD *Abuf = imA.pdata[layer];
        WORD *Bbuf = imB.pdata[layer];
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int y = 0; y < args->fit->ry; y++) {
        	int i = y * args->fit->rx;
#ifdef _OPENMP
#pragma omp critical
#endif
			{
				set_progress_bar_data(NULL, cur_nb / total);
				cur_nb++;
			}
			for (int x = 0; x < args->fit->rx; x++) {
				double r, theta;
				point delta;

				to_polar(x, y, center, &r, &theta);

				// Positive differential
				to_cartesian(r - args->dR, theta + dAlpha, center, &delta);

				if (delta.x < 0)
					delta.x = fabs(delta.x);
				else if (delta.x > wd)
					delta.x = w2 - delta.x;

				if (delta.y < 0)
					delta.y = fabs(delta.y);
				else if (delta.y > hd)
					delta.y = h2 - delta.y;

				gbuf[i] -= Abuf[(int) delta.x + (int) delta.y * args->fit->rx];

				// Negative differential
				to_cartesian(r - args->dR, theta - dAlpha, center, &delta);
				if (delta.x < 0)
					delta.x = fabs(delta.x);
				else if (delta.x > wd)
					delta.x = w2 - delta.x;
				if (delta.y < 0)
					delta.y = fabs(delta.y);
				else if (delta.y > hd)
					delta.y = h2 - delta.y;
				gbuf[i] -= Bbuf[(int) delta.x + (int) delta.y * args->fit->rx];
				i++;
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

    gettimeofday(&t_end, NULL);
    show_time(t_start, t_end);

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

	if ((args->xc >= args->fit->rx) || (args->yc >= args->fit->ry)) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Wrong center coordinates"),
				_("The coordinates cannot be greater than the size of the image. Please change their values and retry."));
	} else {

	set_cursor_waiting(TRUE);

	undo_save_state(&gfit, "Processing: RGradient: (dR=%5.2lf, dA=%4.2lf, xc=%7.1lf, yc=%7.1lf)",
			args->dR, args->da, args->xc, args->yc);

	start_in_new_thread(rgradient_filter, args);
	}
}

void on_button_rgradient_selection_clicked(GtkButton *button, gpointer user_data) {
	if (com.selection.h && com.selection.w) {
		fitted_PSF *result = psf_get_minimisation(&gfit, 0, &com.selection, FALSE, TRUE, TRUE);
		gchar *x0 = g_strdup_printf("%.3lf", result->x0 + com.selection.x);
		gtk_entry_set_text(GTK_ENTRY(lookup_widget("entry_rgradient_xc")), x0);
		gchar *y0 = g_strdup_printf("%.3lf", com.selection.y + com.selection.h - result->y0);
		gtk_entry_set_text(GTK_ENTRY(lookup_widget("entry_rgradient_yc")), y0);

		g_free(x0);
		g_free(y0);
		free(result);
	}
}
