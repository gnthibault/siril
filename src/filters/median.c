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

#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/undo.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/optimize_utils.h"
#include "algos/statistics.h"
#include "algos/sorting.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "io/single_image.h"

#include "median.h"
#include "algos/median_fast.h"

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
 *                M E D I A N     I M A G E     F I L T E R S                *
 ****************************************************************************/

/* get the median of the neighbors of pixel (xx, yy), including itself if
 * include_self is TRUE. radius is 1 for a 3x3, 2 for a 5x5 and so on.
 * w and h are the size of the image passed in buf.
 */
double get_median_ushort(WORD *buf, const int xx, const int yy, const int w,
		const int h, int radius, gboolean is_cfa, gboolean include_self) {
	int n = 0, step = 1, x, y, ksize;
	WORD *values;
	double median;

	if (is_cfa) {
		step = 2;
		radius *= 2;
	}
	ksize = radius * 2 + 1;
	values = calloc(ksize * ksize, sizeof(WORD));

	for (y = yy - radius; y <= yy + radius; y += step) {
		for (x = xx - radius; x <= xx + radius; x += step) {
			if (y >= 0 && y < h && x >= 0 && x < w) {
				// ^ limit to image bounds ^
				// v exclude centre pixel v
				if (include_self || x != xx || y != yy) {
					values[n++] = buf[x + y * w];
				}
			}
		}
	}
	median = quickmedian(values, n);
	free(values);
	return median;
}

double get_median_float(float *buf, const int xx, const int yy, const int w,
		const int h, int radius, gboolean is_cfa, gboolean include_self) {
	int n = 0, step = 1, x, y, ksize;
	float *values;
	double median;

	if (is_cfa) {
		step = 2;
		radius *= 2;
	}
	ksize = radius * 2 + 1;
	values = calloc(ksize * ksize, sizeof(float));

	for (y = yy - radius; y <= yy + radius; y += step) {
		for (x = xx - radius; x <= xx + radius; x += step) {
			if (y >= 0 && y < h && x >= 0 && x < w) {
				// ^ limit to image bounds ^
				// v exclude centre pixel v
				if (include_self || x != xx || y != yy) {
					values[n++] = buf[x + y * w];
				}
			}
		}
	}
	median = quickmedian_float(values, n);
	free(values);
	return median;
}

float get_median_float_fast(float *buf, const int xx, const int yy, const int w,
		const int h, int radius) {

	int ksize = radius * 2 + 1;
	float values[ksize * ksize];

	int ystart = (yy - radius) < 0 ? 0 : yy - radius;
	int yend = (yy + radius) >= h ? h - 1 : yy + radius;
	int xstart = (xx - radius) < 0 ? 0 : xx - radius;
	int xend = (xx + radius) >= w ? w - 1 : xx + radius;
	int n = 0;
	for (int y = ystart; y <= yend; ++y) {
		for (int x = xstart; x <= xend; ++x) {
			// ^ limit to image bounds ^
			values[n++] = buf[x + y * w];
		}
	}
	return quickmedian_float(values, n);
}

double get_median_gsl(gsl_matrix *mat, const int xx, const int yy, const int w,
		const int h, int radius, gboolean is_cfa, gboolean include_self) {
	int n = 0, step = 1, x, y, ksize;
	double *values, median;

	if (is_cfa) {
		step = 2;
		radius *= 2;
	}
	ksize = radius * 2 + 1;
	values = calloc(ksize * ksize, sizeof(double));

	for (y = yy - radius; y <= yy + radius; y += step) {
		for (x = xx - radius; x <= xx + radius; x += step) {
			if (y >= 0 && y < h && x >= 0 && x < w) {
				// ^ limit to image bounds ^
				// v exclude centre pixel v
				if (include_self || x != xx || y != yy) {
					values[n++] = gsl_matrix_get(mat, y, x);
				}
			}
		}
	}
	median = quickmedian_double(values, n);
	free(values);
	return median;
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

	free(args);
	return FALSE;
}

static gpointer median_filter_ushort(gpointer p) {
	struct median_filter_data *args = (struct median_filter_data *)p;
	int progress = 0, x, y, layer, iter = 0;
	int nx = args->fit->rx;
	int ny = args->fit->ry;
	double total, norm = get_normalized_value(args->fit);
	struct timeval t_start, t_end;
	int radius = (args->ksize - 1) / 2;

	g_assert(args->ksize % 2 == 1 && args->ksize > 1);
	g_assert(nx > 0 && ny > 0);
	total = ny * args->fit->naxes[2] * args->iterations;

	char *msg = siril_log_color_message(_("Median Filter: processing...\n"), "red");
	msg[strlen(msg) - 1] = '\0';
	set_progress_bar_data(msg, PROGRESS_RESET);
	gettimeofday(&t_start, NULL);

	do {
		for (layer = 0; layer < args->fit->naxes[2]; layer++) {
			WORD *data = args->fit->pdata[layer];
			for (y = 0; y < ny; y++) {
				int pix_idx = y * nx;
				if (!get_thread_run()) break;
				if (!(y % 16))	// every 16 iterations
					set_progress_bar_data(NULL, (double)progress / total);
				progress++;
				for (x = 0; x < nx; x++) {
					double median = get_median_ushort(data, x, y, nx, ny, radius, FALSE, TRUE);
					if (args->amount != 1.0) {
						double pixel = args->amount * (median / norm);
						pixel += (1.0 - args->amount)
							* ((double)data[pix_idx] / norm);
						data[pix_idx] = round_to_WORD(pixel * norm);
					} else {
						data[pix_idx] = round_to_WORD(median);
					}
					pix_idx++;
				}
			}
		}
		iter++;
	} while (iter < args->iterations && get_thread_run());
	invalidate_stats_from_fit(args->fit);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	set_progress_bar_data(_("Median filter applied"), PROGRESS_DONE);
	siril_add_idle(end_median_filter, args);

	return GINT_TO_POINTER(0);
}

static gpointer median_filter_float(gpointer p) {
	struct median_filter_data *args = (struct median_filter_data *)p;
	int progress = 0;
	int nx = args->fit->rx;
	int ny = args->fit->ry;
	double total;
	struct timeval t_start, t_end;
	int radius = (args->ksize - 1) / 2;

	g_assert(args->ksize % 2 == 1 && args->ksize > 1);
	g_assert(nx > 0 && ny > 0);
	total = ny * args->fit->naxes[2] * args->iterations;

	char *msg = siril_log_color_message(_("Median Filter: processing...\n"), "red");
	msg[strlen(msg) - 1] = '\0';
	set_progress_bar_data(msg, PROGRESS_RESET);
	gettimeofday(&t_start, NULL);

	size_t alloc_size = args->fit->naxes[0] * args->fit->naxes[1] * sizeof(float);
	float *temp = calloc(1, alloc_size); // we need a temporary buffer
	if (!temp) {
		PRINT_ALLOC_ERR;
		return GINT_TO_POINTER(-1);
	}
	float amountf = args->amount;
	for (int layer = 0; layer < args->fit->naxes[2]; layer++) {
		for (int iter = 0; iter < args->iterations; ++iter) {
			float *dst = (iter % 2) ? args->fit->fpdata[layer] : temp;
			float *src = (iter % 2) ? temp : args->fit->fpdata[layer];
			// borders
			for (int y = 0; y < ny; y++) {
				if (y < radius || y >= ny - radius) {
					for (int x = 0; x < nx; x++) {
						if (x < radius || x >= nx - radius) {
							int pix_idx = y * nx + x;
							float median = get_median_float_fast(src, x, y, nx, ny, radius);
							dst[pix_idx] = intpf(amountf, median, src[pix_idx]);
							pix_idx++;
						}
					}
				}
			}
			if (args->ksize == 3) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,16) num_threads(com.max_thread)
#endif
				for (int y = 1; y < ny - 1; y++) {
					int pix_idx = y * nx + 1;
					int x = 1;
#ifdef __SSE2__
					for (; x <= nx - 4; x += 4) {
						__m128 medianv = median9sse(_mm_loadu_ps(&src[(y - 1) * nx + x - 1]),
								_mm_loadu_ps(&src[(y - 1) * nx + x]),
								_mm_loadu_ps(&src[(y - 1) * nx + x + 1]),
								_mm_loadu_ps(&src[y * nx + x - 1]),
								_mm_loadu_ps(&src[y * nx + x]),
								_mm_loadu_ps(&src[y * nx + x + 1]),
								_mm_loadu_ps(&src[(y + 1) * nx + x - 1]),
								_mm_loadu_ps(&src[(y + 1) * nx + x]),
								_mm_loadu_ps(&src[(y + 1) * nx + x + 1])
								);
						_mm_storeu_ps(&dst[pix_idx], intpsse(_mm_set1_ps(amountf), medianv, _mm_loadu_ps(&src[pix_idx])));
						pix_idx += 4;
					}
#endif
					for (; x < nx - 1; x++) {
						float median = median9f(src[(y - 1) * nx + x - 1],
								src[(y - 1) * nx + x],
								src[(y - 1) * nx + x + 1],
								src[y * nx + x - 1],
								src[y * nx + x],
								src[y * nx + x + 1],
								src[(y + 1) * nx + x - 1],
								src[(y + 1) * nx + x],
								src[(y + 1) * nx + x + 1]);
						dst[pix_idx] = intpf(amountf, median, src[pix_idx]);
						pix_idx++;
					}
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						++progress;
						if (!(progress % 32)) {
							set_progress_bar_data(NULL, (double)progress / total);
						}
					}
				}
			} else if (args->ksize == 5) {
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
				{
					float medbuf[25];
#ifdef __SSE2__
					__m128 medbufv[25];
#endif

#ifdef _OPENMP
#pragma omp for schedule(dynamic,16)
#endif
					for (int y = 2; y < ny - 2; y++) {
						int pix_idx = y * nx + 2;
						int x = 2;
#ifdef __SSE2__
						for (; x <= nx - 5; x += 4) {
							int nb = 0;
							for (int i = -2; i <= 2; ++i) {
								for (int j = -2; j <= 2; ++j) {
									_mm_storeu_ps((float*)&medbufv[nb++], _mm_loadu_ps(&src[(y + i) * nx + x + j]));
								}
							}
							__m128 medianv = median5x5sse(medbufv);
							__m128 amountv = _mm_set1_ps(amountf);
							_mm_storeu_ps(&dst[pix_idx], intpsse(amountv, medianv, _mm_loadu_ps(&src[pix_idx])));
							pix_idx += 4;
						}
#endif
						for (; x < nx - 2; x++) {
							int nb = 0;
							for (int i = -2; i <= 2 ; ++i) {
								for (int j = -2; j <= 2; ++j) {
									medbuf[nb++] = src[(y + i) * nx + x + j];
								}
							}
							float median = median5x5(medbuf);
							dst[pix_idx] = intpf(amountf, median, src[pix_idx]);
							pix_idx++;
						}
#ifdef _OPENMP
#pragma omp critical
#endif
						{
							++progress;
							if (!(progress % 32)) {
								set_progress_bar_data(NULL, (double)progress / total);
							}
						}
					}
				}
			} else if (args->ksize == 7) {
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
				{
					float medbuf[49];
#ifdef __SSE2__
					__m128 medbufv[49];
#endif

#ifdef _OPENMP
#pragma omp for schedule(dynamic,16)
#endif
					for (int y = 3; y < ny - 3; y++) {
						int pix_idx = y * nx + 3;
						int x = 3;
#ifdef __SSE2__
						for (; x < nx - 6; x += 4) {
							int nb = 0;
							for (int i = -3; i <= 3; ++i) {
								for (int j = -3; j <= 3; ++j) {
									_mm_storeu_ps((float*)&medbufv[nb++], _mm_loadu_ps(&src[(y + i) * nx + x + j]));
								}
							}
							__m128 medianv = median7x7sse(medbufv);
							__m128 amountv = _mm_set1_ps(amountf);
							_mm_storeu_ps(&dst[pix_idx], intpsse(amountv, medianv, _mm_loadu_ps(&src[pix_idx])));
							pix_idx += 4;
						}
#endif
						for (; x < nx - 3; x++) {
							int nb = 0;
							for (int i = -3; i <= 3 ; ++i) {
								for (int j = -3; j <= 3; ++j) {
									medbuf[nb++] = src[(y + i) * nx + x + j];
								}
							}
							float median = median7x7(medbuf);
							dst[pix_idx] = intpf(amountf, median, src[pix_idx]);
							pix_idx++;
						}
#ifdef _OPENMP
#pragma omp critical
#endif
						{
							++progress;
							if (!(progress % 32)) {
								set_progress_bar_data(NULL, (double)progress / total);
							}
						}
					}
				}
			} else if (args->ksize == 9) {
#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
				{
					float medbuf[81];
#ifdef __SSE2__
					__m128 medbufv[81];
#endif

#ifdef _OPENMP
#pragma omp for schedule(dynamic,16)
#endif
					for (int y = 4; y < ny - 4; y++) {
						int pix_idx = y * nx + 4;
						int x = 4;
#ifdef __SSE2__
						for (; x < nx - 7; x += 4) {
							int nb = 0;
							for (int i = -4; i <= 4; ++i) {
								for (int j = -4; j <= 4; ++j) {
									_mm_storeu_ps((float*)&medbufv[nb++], _mm_loadu_ps(&src[(y + i) * nx + x + j]));
								}
							}
							__m128 medianv = median9x9sse(medbufv);
							__m128 amountv = _mm_set1_ps(amountf);
							_mm_storeu_ps(&dst[pix_idx], intpsse(amountv, medianv, _mm_loadu_ps(&src[pix_idx])));
							pix_idx += 4;
						}
#endif
						for (; x < nx - 4; x++) {
							int nb = 0;
							for (int i = -4; i <= 4 ; ++i) {
								for (int j = -4; j <= 4; ++j) {
									medbuf[nb++] = src[(y + i) * nx + x + j];
								}
							}
							float median = median9x9(medbuf);
							dst[pix_idx] = intpf(amountf, median, src[pix_idx]);
							pix_idx++;
						}
#ifdef _OPENMP
#pragma omp critical
#endif
						{
							++progress;
							if (!(progress % 32)) {
								set_progress_bar_data(NULL, (double)progress / total);
							}
						}
					}
				}
			} else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic, 16)
#endif
				for (int y = 0; y < ny; y++) {
					int pix_idx = y * nx;
					for (int x = 0; x < nx; x++) {
						float median = get_median_float_fast(src, x, y, nx, ny, radius);
						dst[pix_idx] = intpf(amountf, median, src[pix_idx]);
						pix_idx++;
					}
#ifdef _OPENMP
#pragma omp critical
#endif
					{
						++progress;
						if (!(progress % 32)) {
							set_progress_bar_data(NULL, (double)progress / total);
						}
					}
				}
			}
		}
		if (args->iterations % 2) {
			// for odd number of iterations (1, 3, 5, ...) we have to copy the data back at the end
			float *dst = args->fit->fpdata[layer];
			float *src = temp;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread)
#endif
			for (int y = 0; y < ny; y++) {
				for (int x = 0; x < nx; x++) {
					dst[y * nx + x] = src[y * nx + x];
				}
			}
		}
	}
	free(temp);
	invalidate_stats_from_fit(args->fit);
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	set_progress_bar_data(_("Median filter applied"), PROGRESS_DONE);
	siril_add_idle(end_median_filter, args);

	return GINT_TO_POINTER(0);
}

/* The function smoothes an image using the median filter with the
 * ksize x ksize aperture. Each channel of a multi-channel image is
 * processed independently. In-place operation is supported. */
gpointer median_filter(gpointer p) {
	struct median_filter_data *args = (struct median_filter_data *)p;
	if (args->fit->type == DATA_USHORT)
		return median_filter_ushort(p);
	if (args->fit->type == DATA_FLOAT)
		return median_filter_float(p);
	siril_add_idle(end_median_filter, args);
	return GINT_TO_POINTER(1);
}
