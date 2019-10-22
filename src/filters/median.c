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

#include "core/siril.h"
#include "core/proto.h"
#include "core/undo.h"
#include "core/processing.h"
#include "algos/statistics.h"
#include "algos/sorting.h"
#include "gui/progress_and_log.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "io/single_image.h"

#include "median.h"


void on_menuitem_medianfilter_activate(GtkMenuItem *menuitem,
		gpointer user_data) {
	if (single_image_is_loaded())
		siril_open_dialog("Median_dialog");
}


void on_Median_cancel_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("Median_dialog");
}

void on_Median_Apply_clicked(GtkButton *button, gpointer user_data) {
	int combo_size = gtk_combo_box_get_active(
			GTK_COMBO_BOX(
					gtk_builder_get_object(builder, "combo_ksize_median")));
	double amount = gtk_range_get_value(
			GTK_RANGE(gtk_builder_get_object(builder, "scale_median")));
	int iterations = round_to_int(gtk_spin_button_get_value(GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "median_button_iterations"))));

	if (get_thread_run()) {
		siril_log_message(
				_(	"Another task is already in progress, ignoring new request.\n"));
		return;
	}

	struct median_filter_data *args = malloc(sizeof(struct median_filter_data));

	switch (combo_size) {
	default:
	case 0:
		args->ksize = 3;
		break;
	case 1:
		args->ksize = 5;
		break;
	case 2:
		args->ksize = 7;
		break;
	case 3:
		args->ksize = 9;
		break;
	case 4:
		args->ksize = 11;
		break;
	case 5:
		args->ksize = 13;
		break;
	case 6:
		args->ksize = 15;
		break;
	}
	undo_save_state(&gfit, "Processing: Median Filter (filter=%dx%d px)",
			args->ksize, args->ksize);

	args->fit = &gfit;
	args->amount = amount;
	args->iterations = iterations;
	set_cursor_waiting(TRUE);
	start_in_new_thread(median_filter, args);

}

/*****************************************************************************
 *                      M E D I A N     F I L T E R                          *
 ****************************************************************************/

static gboolean end_median_filter(gpointer p) {
	struct median_filter_data *args = (struct median_filter_data *) p;
	stop_processing_thread();// can it be done here in case there is no thread?
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	update_used_memory();
	free(args);
	return FALSE;
}

/* The function smoothes an image using the median filter with the
 * ksize x ksize aperture. Each channel of a multi-channel image is
 * processed independently. In-place operation is supported. */
gpointer median_filter(gpointer p) {
	struct median_filter_data *args = (struct median_filter_data *) p;
	g_assert(args->ksize % 2 == 1 && args->ksize > 1);
	int i, x, y, xx, yy, layer, iter = 0;
	int nx = args->fit->rx;
	int ny = args->fit->ry;
	int radius = (args->ksize - 1) / 2;
	int ksize_squared = args->ksize * args->ksize;
	double norm = (double) get_normalized_value(args->fit);
	double cur = 0.0, total;
	g_assert(nx > 0 && ny > 0);

	struct timeval t_start, t_end;

	char *msg = siril_log_color_message(_("Median Filter: processing...\n"), "red");
	msg[strlen(msg) - 1] = '\0';
	set_progress_bar_data(msg, PROGRESS_RESET);
	gettimeofday(&t_start, NULL);

	WORD *data = calloc(ksize_squared, sizeof(WORD));
	if (data == NULL) {
		PRINT_ALLOC_ERR;
		siril_add_idle(end_median_filter, args);
		set_progress_bar_data(_("Median filter failed"), PROGRESS_DONE);
		return GINT_TO_POINTER(1);
	}

	do {
		for (layer = 0; layer < args->fit->naxes[2]; layer++) {
			/* FILL image upside-down */
			WORD **image = malloc(ny * sizeof(WORD *));
			if (image == NULL) {
				PRINT_ALLOC_ERR;
				siril_add_idle(end_median_filter, args);
				free(data);
				return GINT_TO_POINTER(1);
			}
			for (i = 0; i < ny; i++)
				image[ny - i - 1] = args->fit->pdata[layer] + i * nx;

			for (y = 0; y < ny; y++) {
				if (!get_thread_run())
					break;
				total = ny * args->fit->naxes[2] * args->iterations;
				if (!(y % 16))	// every 16 iterations
					set_progress_bar_data(NULL, cur / total);
				cur++;
				for (x = 0; x < nx; x++) {
					i = 0;
					for (yy = y - radius; yy <= y + radius; yy++) {
						for (xx = x - radius; xx <= x + radius; xx++) {
							WORD tmp;
							if (xx < 0 && yy >= 0) {
								if (yy >= ny)
									tmp = image[ny - 1][0];
								else
									tmp = image[yy][0];
							} else if (xx > 0 && yy <= 0) {
								if (xx >= nx)
									tmp = image[0][nx - 1];
								else
									tmp = image[0][xx];
							} else if (xx <= 0 && yy <= 0) {
								tmp = image[0][0];
							} else {
								if (xx >= nx && yy >= ny)
									tmp = image[ny - 1][nx - 1];
								else if (xx >= nx && yy < ny)
									tmp = image[yy][nx - 1];
								else if (xx < nx && yy >= ny)
									tmp = image[ny - 1][xx];
								else
									tmp = image[yy][xx];
							}
							data[i++] = tmp;
						}
					}
					WORD median = round_to_WORD(quickmedian(data,ksize_squared));
					double pixel = args->amount * (median / norm);
					pixel += (1.0 - args->amount)
							* ((double) image[y][x] / norm);
					image[y][x] = round_to_WORD(pixel * norm);
				}
			}
			free(image);
		}
		iter++;
	} while (iter < args->iterations && get_thread_run());
	invalidate_stats_from_fit(args->fit);
	free(data);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	set_progress_bar_data(_("Median filter applied"), PROGRESS_DONE);
	siril_add_idle(end_median_filter, args);

	return GINT_TO_POINTER(0);
}
