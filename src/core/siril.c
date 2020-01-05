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

#include <gtk/gtk.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/arithm.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "gui/image_display.h"
#include "gui/histogram.h"
#include "gui/progress_and_log.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "algos/colors.h"
#include "algos/statistics.h"
#include "opencv/opencv.h"

int threshlo(fits *fit, int level) {
	int i, layer;

	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		if (fit->type == DATA_USHORT) {
			WORD *buf = fit->pdata[layer];
			for (i = 0; i < fit->rx * fit->ry; ++i) {
				*buf = max(level, *buf);
				buf++;
			}
		} else if (fit->type == DATA_FLOAT) {
			float *buf = fit->fpdata[layer];
			for (i = 0; i < fit->rx * fit->ry; ++i) {
				*buf = max(level, *buf);
				buf++;
			}
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int threshhi(fits *fit, int level) {
	int i, layer;

	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		if (fit->type == DATA_USHORT) {
			WORD *buf = fit->pdata[layer];
			for (i = 0; i < fit->rx * fit->ry; ++i) {
				*buf = min(level, *buf);
				buf++;
			}
		} else if (fit->type == DATA_FLOAT) {
			float *buf = fit->fpdata[layer];
			for (i = 0; i < fit->rx * fit->ry; ++i) {
				*buf = min(level, *buf);
				buf++;
			}
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int nozero(fits *fit, int level) {
	int i, layer;

	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		if (fit->type == DATA_USHORT) {
			WORD *buf = fit->pdata[layer];
			for (i = 0; i < fit->rx * fit->ry; ++i) {
				if (*buf == 0)
					*buf = level;
				buf++;
			}
		} else if (fit->type == DATA_FLOAT) {
			float *buf = fit->fpdata[layer];
			for (i = 0; i < fit->rx * fit->ry; ++i) {
				if (*buf == 0)
					*buf = level;
				buf++;
			}
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

/**********************************************************
 *
 */

int unsharp(fits *fit, double sigma, double amount, gboolean verbose) {
	struct timeval t_start, t_end;

	if (sigma <= 0.0)
		return 1;
	if (verbose) {
		siril_log_color_message(_("Unsharp: processing...\n"), "red");
		gettimeofday(&t_start, NULL);
	}
	cvUnsharpFilter(fit, sigma, amount);

	if (verbose) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}
	return 0;
}

/* This entropy function computes the entropy for the image in gfit for its
 * layer 'layer', in the area designated by area which can be NULL.
 * An optional imstats parameter can be used to provide the background and
 * sigma value, and when it is given, the entropy will only be computed for
 * pixels with values above background + 1 * sigma. It must be NULL otherwise.
 */
double entropy(fits *fit, int layer, rectangle *area, imstats *opt_stats) {
	double e = 0.0, threshold = 0.0;
	gsl_histogram *histo;
	size_t i, size, n;

	if (opt_stats && opt_stats->median >= 0.0 && opt_stats->sigma >= 0.0)
		threshold = opt_stats->median + 1 * opt_stats->sigma;

	if (area == NULL)
		histo = computeHisto(fit, layer);
	else
		histo = computeHisto_Selection(fit, layer, area);

	n = fit->rx * fit->ry;
	g_assert (n > 0);
	size = gsl_histogram_bins(histo);
	for (i = 0; i < size; i++) {
		double p = gsl_histogram_get(histo, i);
		if (p > threshold && p < size)
			e += (p / n) * log(n / p);
	}
	gsl_histogram_free(histo);

	return e;
}

static int loglut_ushort(fits *fit) {
	// This function maps fit with a log LUT
	int i, layer;
	WORD *buf[3] = { fit->pdata[RLAYER],
			fit->pdata[GLAYER], fit->pdata[BLAYER] };

	double norm = USHRT_MAX_DOUBLE / log(USHRT_MAX_DOUBLE);

	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		imstats *stat = statistics(NULL, -1, fit, layer, NULL, STATS_MINMAX);
		double min = stat->min;
		double wd = stat->max - stat->min;
		for (i = 0; i < fit->ry * fit->rx; i++) {
			float px = (float)buf[layer][i];
			buf[layer][i] = round_to_WORD(logf((px - min) / wd) * norm);
		}
		free_stats(stat);
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

static int loglut_float(fits *fit) {
	// This function maps fit with a log LUT
	int i, layer;
	float *buf[3] = { fit->fpdata[RLAYER],
			fit->fpdata[GLAYER], fit->fpdata[BLAYER] };

	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		imstats *stat = statistics(NULL, -1, fit, layer, NULL, STATS_MINMAX);
		double min = stat->min;
		double wd = stat->max - stat->min;
		for (i = 0; i < fit->ry * fit->rx; i++) {
			float px = buf[layer][i];
			buf[layer][i] = log1pf((px - min) / wd);
		}
		free_stats(stat);
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int loglut(fits *fit) {
	if (fit->type == DATA_USHORT) {
		return loglut_ushort(fit);
	} else if (fit->type == DATA_FLOAT) {
		return loglut_float(fit);
	}
	return -1;
}

int ddp(fits *a, int level, float coeff, float sigma) {
	fits fit = { 0 };
	// TODO: check retvals
	copyfits(a, &fit, CP_ALLOC | CP_COPYA | CP_FORMAT, 0);
	unsharp(&fit, sigma, 0, FALSE);
	soper(&fit, (double)level, OPER_ADD, TRUE);
	nozero(&fit, 1);
	siril_fdiv(a, &fit, level, TRUE);
	soper(a, (double)coeff, OPER_MUL, TRUE);
	clearfits(&fit);
	invalidate_stats_from_fit(a);
	return 0;
}

int visu(fits *fit, int low, int high) {
	if (low < 0 || low > USHRT_MAX || high < 1 || high > USHRT_MAX)
		return 1;
	if (single_image_is_loaded() && com.cvport < com.uniq->nb_layers) {
		com.uniq->layers[com.cvport].hi = high;
		com.uniq->layers[com.cvport].lo = low;
	} else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers) {
		com.seq.layers[com.cvport].hi = high;
		com.seq.layers[com.cvport].lo = low;
	} else
		return 1;
	set_cutoff_sliders_values();
	redraw(com.cvport, REMAP_ONLY);
	redraw_previews();
	return 0;
}

/* fill an image or selection with the value 'level' */
int fill(fits *fit, int level, rectangle *arearg) {
	int i, j, layer;
	rectangle area;

	if (arearg) {
		memcpy(&area, arearg, sizeof(rectangle));
	} else {
		if (com.selection.h && com.selection.w) {
			memcpy(&area, &com.selection, sizeof(rectangle));
		} else {
			area.w = fit->rx;
			area.h = fit->ry;
			area.x = 0;
			area.y = 0;
		}
	}
	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		if (fit->type == DATA_USHORT) {
			WORD *buf = fit->pdata[layer]
					+ (fit->ry - area.y - area.h) * fit->rx + area.x;
			int stridebuf = fit->rx - area.w;
			for (i = 0; i < area.h; ++i) {
				for (j = 0; j < area.w; ++j) {
					*buf++ = level;
				}
				buf += stridebuf;
			}
		} else if (fit->type == DATA_FLOAT) {
			float *buf = fit->fpdata[layer]
					+ (fit->ry - area.y - area.h) * fit->rx + area.x;
			int stridebuf = fit->rx - area.w;
			for (i = 0; i < area.h; ++i) {
				for (j = 0; j < area.w; ++j) {
					*buf++ = level;
				}
				buf += stridebuf;
			}
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

static int off_ushort(fits *fit, float level) {
	WORD *buf[3] = { fit->pdata[RLAYER], fit->pdata[GLAYER],
			fit->pdata[BLAYER] };
	int i, layer;
	g_assert(fit->naxes[2] <= 3);
	if (level == 0)
		return 0;
	if (level < -USHRT_MAX)
		level = -USHRT_MAX;
	else if (level > USHRT_MAX)
		level = USHRT_MAX;
	for (i = 0; i < fit->rx * fit->ry; ++i) {
		for (layer = 0; layer < fit->naxes[2]; ++layer) {
			float val = (float)buf[layer][i];
			buf[layer][i] = roundf_to_WORD(val + level);
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

static int off_float(fits *fit, float level) {
	float *buf[3] = { fit->fpdata[RLAYER], fit->fpdata[GLAYER],
			fit->fpdata[BLAYER] };
	int i, layer;
	g_assert(fit->naxes[2] <= 3);
	if (level == 0)
		return 0;
	if (level < -1.f)
		level = -1.f;
	else if (level > 1.f)
		level = 1.f;
	for (i = 0; i < fit->rx * fit->ry; ++i) {
		for (layer = 0; layer < fit->naxes[2]; ++layer) {
			float val = buf[layer][i];
			buf[layer][i] = val + level;
			buf[layer][i] = buf[layer][i] > 1.f ? 1.f : buf[layer][i];
			buf[layer][i] = buf[layer][i] < 0.f ? 0.f : buf[layer][i];
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int off(fits *fit, float level) {
	if (fit->type == DATA_USHORT) {
		return off_ushort(fit, level);
	} else if (fit->type == DATA_FLOAT) {
		return off_float(fit, level);
	}
	return -1;
}

/* computes the background value using the histogram and/or median value.
 * The argument layer can be -1 for automatic setting (= green for RGB) */
double background(fits* fit, int reqlayer, rectangle *selection) {
	int layer = RLAYER;
	double bg;

	if (reqlayer >= 0)
		layer = reqlayer;
	else if (isrgb(&gfit))
		layer = GLAYER;		//GLAYER is better to evaluate background

	imstats* stat = statistics(NULL, -1, fit, layer, selection, STATS_BASIC);
	if (!stat) {
		siril_log_message(_("Error: statistics computation failed.\n"));
		return 0.0;
	}
	bg = stat->median;
	free_stats(stat);
	return bg;
}

void show_FITS_header(fits *fit) {
	if (fit->header)
		show_data_dialog(fit->header, "FITS Header");
}

void compute_grey_flat(fits *fit) {
	double mean[4];
	double diag1, diag2, coeff1, coeff2;
	int config;

	/* compute means of 4 channels */
	compute_means_from_flat_cfa(fit, mean);

	/* compute coefficients */

	/* looking for green diagonal */
	diag1 = mean[0] / mean[3];
	diag2 = mean[1] / mean[2];

	/* BAYER_FILTER_RGGB
	 * BAYER_FILTER_BGGR */
	if (fabs(1 - diag1) < fabs(1 - diag2)) {
		coeff1 = mean[1] / mean[0];
		coeff2 = mean[2] / mean[3];
		config = 0;
	} /* BAYER_FILTER_GBRG
	 * BAYER_FILTER_GRBG */
	else {
		coeff1 = mean[0] / mean[1];
		coeff2 = mean[3] / mean[2];
		config = 1;
	}

	/* applies coefficients to cfa image */
	equalize_cfa_fit_with_coeffs(fit, coeff1, coeff2, config);

}
