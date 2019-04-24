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

#include <gtk/gtk.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_statistics.h>
#include <fitsio.h>
#include <complex.h>
#include <float.h>
#include <assert.h>
#include <libgen.h>
#include <fcntl.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "gui/histogram.h"
#include "gui/progress_and_log.h"
#include "gui/PSF_list.h"
#include "io/conversion.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/ser.h"
#include "algos/colors.h"
#include "algos/Def_Math.h"
#include "algos/Def_Wavelet.h"
#include "algos/cosmetic_correction.h"
#include "algos/statistics.h"
#include "algos/sorting.h"
#include "algos/plateSolver.h"
#include "opencv/opencv.h"

#define MAX_ITER 15
#define EPSILON 1E-4

/* this file contains all functions for image processing */

int threshlo(fits *fit, int level) {
	int i, layer;

	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		WORD *buf = fit->pdata[layer];
		for (i = 0; i < fit->rx * fit->ry; ++i) {
			*buf = max(level, *buf);
			buf++;
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int threshhi(fits *fit, int level) {
	int i, layer;

	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		WORD *buf = fit->pdata[layer];
		for (i = 0; i < fit->rx * fit->ry; ++i) {
			*buf = min(level, *buf);
			buf++;
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int nozero(fits *fit, int level) {
	int i, layer;

	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		WORD *buf = fit->pdata[layer];
		for (i = 0; i < fit->rx * fit->ry; ++i) {
			if (*buf == 0)
				*buf = level;
			buf++;
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

/*****************************************************************************
 *       S I R I L      A R I T H M E T I C      O P E R A T I O N S         *
 ****************************************************************************/

/* equivalent to (map simple_operation a), with simple_operation being
 * (lambda (pixel) (oper pixel scalar))
 * oper is a for addition, s for substraction (i for difference) and so on. */
int soper(fits *a, double scalar, char oper) {
	WORD *gbuf;
	int i, layer;
	int n = a->rx * a->ry;

	assert(n > 0);

	for (layer = 0; layer < a->naxes[2]; ++layer) {
		gbuf = a->pdata[layer];
		switch (oper) {
		case OPER_ADD:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD((double) gbuf[i] + scalar);
			}
			break;
		case OPER_SUB:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD((double) gbuf[i] - scalar);
			}
			break;
		case OPER_MUL:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD((double) gbuf[i] * scalar);
			}
			break;
		case OPER_DIV:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD((double) gbuf[i] / scalar);
			}
			break;
		}
	}
	invalidate_stats_from_fit(a);
	return 0;
}

/* applies operation of image a with image b, for all their layers:
 * a = a oper b
 * returns 0 on success */
int imoper(fits *a, fits *b, char oper) {
	int i, layer;

	if (a->rx != b->rx || a->ry != b->ry) {
		siril_log_message(
				_("imoper: images don't have the same size (w = %u|%u, h = %u|%u)\n"),
				a->rx, b->rx, a->ry, b->ry);
		return 1;
	}
	for (layer = 0; layer < a->naxes[2]; ++layer) {
		WORD *buf = b->pdata[layer];
		WORD *gbuf = a->pdata[layer];
		int n = a->rx * a->ry;
		switch (oper) {
		case OPER_ADD:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD(gbuf[i] + buf[i]);
			}
			break;
		case OPER_SUB:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD(gbuf[i] - buf[i]);
			}
			break;
		case OPER_MUL:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD(gbuf[i] * buf[i]);
			}
			break;
		case OPER_DIV:
			for (i = 0; i < n; ++i) {
				gbuf[i] = (buf[i] == 0) ? 0 : round_to_WORD(gbuf[i] / buf[i]);
			}
			break;
		}
	}
	invalidate_stats_from_fit(a);
	return 0;
}

int addmax(fits *a, fits *b) {
	WORD *gbuf[3] = { a->pdata[RLAYER], a->pdata[GLAYER], a->pdata[BLAYER] };
	WORD *buf[3] = { b->pdata[RLAYER], b->pdata[GLAYER], b->pdata[BLAYER] };
	gint i, layer;

	if (a->rx != b->rx || a->ry != b->ry || a->naxes[2] != b->naxes[2]) {
		siril_log_message(
				_("addmax: images don't have the same size (w = %d|%d, h = %d|%d, layers = %d|%d)\n"),
				a->rx, b->rx, a->ry, b->ry, a->naxes[2], b->naxes[2]);
		return 1;
	}
	assert(a->naxes[2] == 1 || a->naxes[2] == 3);

	for (layer = 0; layer < a->naxes[2]; ++layer) {
		for (i = 0; i < a->ry * a->rx; ++i) {
			if (buf[layer][i] > gbuf[layer][i])
				gbuf[layer][i] = buf[layer][i];
		}
	}
	invalidate_stats_from_fit(a);
	return 0;
}

/* If siril_fdiv is ok, function returns 0. If overflow, siril_fdiv returns 1*/
int siril_fdiv(fits *a, fits *b, float coef) {
	int i, layer;
	int retvalue = 0;
	double temp;

	if (a->rx != b->rx || a->ry != b->ry || a->naxes[2] != b->naxes[2]) {
		fprintf(stderr, "Wrong size or channel count: %u=%u? / %u=%u?\n", a->rx,
				b->rx, a->ry, b->ry);
		return -1;
	}
	for (layer = 0; layer < a->naxes[2]; ++layer) {
		WORD *buf = b->pdata[layer];
		WORD *gbuf = a->pdata[layer];
		for (i = 0; i < b->rx * b->ry; ++i) {
			if (buf[i] == 0)
				buf[i] = 1;		// avoid division by 0
			temp = ((double) coef * ((double) gbuf[i] / (double) buf[i]));
			if (temp > USHRT_MAX_DOUBLE)
				retvalue = 1;
			gbuf[i] = round_to_WORD(temp);
		}
	}
	invalidate_stats_from_fit(a);
	return retvalue;
}

/* normalized division a/b, stored in a, with max value equal to the original
 * max value of a, for each layer. */
int siril_ndiv(fits *a, fits *b) {
	double *div;
	int layer, i, nb_pixels;
	if (a->rx != b->rx || a->ry != b->ry || a->naxes[2] != b->naxes[2]) {
		fprintf(stderr,
				"Wrong size or channel count: %u=%u? / %u=%u?, %ld=%ld?\n",
				a->rx, b->rx, a->ry, b->ry, a->naxes[2], b->naxes[2]);
		return 1;
	}
	nb_pixels = a->rx * a->ry;
	div = malloc(nb_pixels * sizeof(double));

	for (layer = 0; layer < a->naxes[2]; ++layer) {
		double max = 0, norm;
		for (i = 0; i < nb_pixels; ++i) {
			if (!b->pdata[layer][i])
				div[i] = (double) a->pdata[layer][i];
			else
				div[i] = (double) a->pdata[layer][i]
						/ (double) b->pdata[layer][i];
			max = max(div[i], max);
		}
		norm = max / fit_get_max(a, layer);
		for (i = 0; i < nb_pixels; ++i) {
			a->pdata[layer][i] = round_to_WORD(div[i] / norm);
		}
	}

	invalidate_stats_from_fit(a);
	free(div);
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

/* takes the image in gfit, copies it in a temporary fit to shift it, and copy it back into gfit */
/* TODO: it can be done in the same, thus avoiding to allocate, it just needs to care
 * about the sign of sx and sy to avoid data overwriting in the same allocated space. */
int shift(int sx, int sy) {
	int x, y, nx, ny, i, ii, layer;
	fits tmpfit;
	copyfits(&(gfit), &tmpfit, CP_ALLOC | CP_FORMAT, 0);
	i = 0;
	/* the loop is the same than in composit() */
	for (y = 0; y < gfit.ry; ++y) {
		for (x = 0; x < gfit.rx; ++x) {
			nx = (x - sx);
			ny = (y - sy);
			//printf("x=%d y=%d sx=%d sy=%d i=%d ii=%d\n",x,y,shiftx,shifty,i,ii);
			if (nx >= 0 && nx < gfit.rx && ny >= 0 && ny < gfit.ry) {
				ii = ny * gfit.rx + nx;
				//printf("shiftx=%d shifty=%d i=%d ii=%d\n",shiftx,shifty,i,ii);
				if (ii > 0 && ii < gfit.rx * gfit.ry) {
					for (layer = 0; layer < gfit.naxes[2]; ++layer) {
						tmpfit.pdata[layer][i] = gfit.pdata[layer][ii];
					}
				}
			}
			++i;
		}
	}

	for (layer = 0; layer < gfit.naxes[2]; ++layer) {
		memcpy(gfit.pdata[layer], tmpfit.pdata[layer],
				gfit.rx * gfit.ry * sizeof(WORD));
	}
	free(tmpfit.data);
	invalidate_stats_from_fit(&gfit);

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
	assert (n > 0);
	size = gsl_histogram_bins(histo);
	for (i = 0; i < size; i++) {
		double p = gsl_histogram_get(histo, i);
		if (p > threshold && p < size)
			e += (p / n) * log(n / p);
	}
	gsl_histogram_free(histo);

	return e;
}

int loglut(fits *fit) {
	// This function maps fit with a log LUT
	int i, layer;
	WORD *buf[3] = { fit->pdata[RLAYER],
			fit->pdata[GLAYER], fit->pdata[BLAYER] };

	double norm = USHRT_MAX_DOUBLE / log(USHRT_MAX_DOUBLE);

	for (i = 0; i < fit->ry * fit->rx; i++) {
		for (layer = 0; layer < fit->naxes[2]; ++layer) {
			double px = (double)buf[layer][i];
			px = (px == 0.0) ? 1.0 : px;
			buf[layer][i] = round_to_WORD(norm * log(px));
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int asinhlut(fits *fit, double beta, double offset, gboolean RGBspace) {
	int i, layer;
	WORD *buf[3] = { fit->pdata[RLAYER],
			fit->pdata[GLAYER], fit->pdata[BLAYER] };
	double norm;

	norm = get_normalized_value(fit);

	for (i = 0; i < fit->ry * fit->rx; i++) {
		double x, k;
		if (fit->naxes[2] > 1) {
			double r, g, b;

			r = (double) buf[RLAYER][i] / norm;
			g = (double) buf[GLAYER][i] / norm;
			b = (double) buf[BLAYER][i] / norm;
			/* RGB space */
			if (RGBspace)
				x = 0.2126 * r + 0.7152 * g + 0.0722 * b;
			else
				x = 0.3333 * r + 0.3333 * g + 0.3333 * b;
		} else {
			x = buf[RLAYER][i] / norm;
		}

		k = asinh(beta * x) / (x * asinh(beta));

		for (layer = 0; layer < fit->naxes[2]; ++layer) {
			double px = (double) buf[layer][i] / norm;
			px -= offset;
			px *= k;
			buf[layer][i] = round_to_WORD(px * norm);
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int ddp(fits *a, int level, float coeff, float sigma) {
	fits fit;
	memset(&fit, 0, sizeof(fits));
	copyfits(a, &fit, CP_ALLOC | CP_COPYA | CP_FORMAT, 0);
	unsharp(&fit, sigma, 0, FALSE);
	soper(&fit, (double) level, OPER_ADD);
	nozero(&fit, 1);
	siril_fdiv(a, &fit, level);
	soper(a, (double) coeff, OPER_MUL);
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
	WORD *buf;
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
		buf = fit->pdata[layer] + (fit->ry - area.y - area.h) * fit->rx
				+ area.x;
		int stridebuf = fit->rx - area.w;
		for (i = 0; i < area.h; ++i) {
			for (j = 0; j < area.w; ++j) {
				*buf++ = level;
			}
			buf += stridebuf;
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}

int off(fits *fit, int level) {
	WORD *buf[3] =
			{ fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	int i, layer;
	assert(fit->naxes[2] <= 3);
	if (level == 0)
		return 0;
	if (level < -USHRT_MAX)
		level = -USHRT_MAX;
	else if (level > USHRT_MAX)
		level = USHRT_MAX;
	for (i = 0; i < fit->rx * fit->ry; ++i) {
		for (layer = 0; layer < fit->naxes[2]; ++layer) {
			WORD val = buf[layer][i];
			if ((level < 0 && val < -level))
				buf[layer][i] = 0;
			else if (level > 0 && val > USHRT_MAX - level)
				buf[layer][i] = USHRT_MAX;
			else
				buf[layer][i] = val + level;
		}
	}
	invalidate_stats_from_fit(fit);
	return 0;
}



/* This function fills the data in the lrgb image with LRGB information from l, r, g and b
 * images. Layers are not aligned, images need to be all of the same size.
 * It may be used in the command line, currently unused. */
int lrgb(fits *l, fits *r, fits *g, fits *b, fits *lrgb) {
	//
	// Combines l r g and b components into resulting lrgb
	// We transform each pixel from RGB to HSI,
	// then take I from the luminance l fits and
	// immediately step back to RGB to the working copy
	//
	guint x, y;
	gdouble rr, gg, bb, h, s, i/*, ps3, dps3, qps3, dpi*/;
	gint maxi;
	WORD *pr, *pg, *pb, *dr, *dg, *db, *pl;

	//
	// some stats used to normalize
	//
	if (image_find_minmax(r) || image_find_minmax(g) ||
			image_find_minmax(b) || image_find_minmax(l)) {
		siril_log_color_message(_("Could not compute normalization values for the images, aborting.\n"), "red");
		return -1;
	}
	maxi = max(r->maxi, max(g->maxi, b->maxi));
	//
	// initialize pointers
	//
	pr = r->data;
	pg = g->data;
	pb = b->data;
	pl = l->data;
	dr = lrgb->pdata[RLAYER];
	dg = lrgb->pdata[GLAYER];
	db = lrgb->pdata[BLAYER];
	//
	// some trigo constants
	// we stick to h in radians, not in degrees
	//
	//dpi=2*M_PI;
	//ps3=M_PI/3;
	//dps3=2*M_PI;
	//dps3=2*M_PI/3;
	//qps3=4*M_PI/3;
	//
	// Main loop
	//
	fprintf(stderr, "HSI->RGB %u %u\n", r->ry, r->rx);
	for (y = 0; y < r->ry; y++) {
		for (x = 0; x < r->rx; x++) {
			//
			// First normalize rgb to [0 1]
			//
			rr = (double) (*pr++) / maxi;
			gg = (double) (*pg++) / maxi;
			bb = (double) (*pb++) / maxi;

			rgb_to_hsl(rr, gg, bb, &h, &s, &i);
			//
			// replace luminance
			//
			i = *pl++ / (double) l->maxi;
			//
			// and back to RGB
			hsl_to_rgb(h, s, i, &rr, &gg, &bb);
			//
			// now denormalize and store
			//
			*dr++ = (WORD) (rr * maxi);
			*dg++ = (WORD) (gg * maxi);
			*db++ = (WORD) (bb * maxi);
		}
	}
	return 0;
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

/* This function computes wavelets with the number of Nbr_Plan and
 * extracts plan "Plan" in fit parameters */

int get_wavelet_layers(fits *fit, int Nbr_Plan, int Plan, int Type, int reqlayer) {
	int chan, start, end, retval = 0;
	wave_transf_des wavelet[3];

	assert(fit->naxes[2] <= 3);

	float *Imag = f_vector_alloc(fit->ry * fit->rx);
	if (Imag == NULL)
		return 1;

	if (reqlayer < 0 || reqlayer > 3) {
		start = 0;
		end = fit->naxes[2];
	}
	else {
		start = reqlayer;
		end = start + 1;
	}

	for (chan = start; chan < end; chan++) {
		int Nl, Nc;

		if (wavelet_transform(Imag, fit->ry, fit->rx, &wavelet[chan],
					Type, Nbr_Plan, fit->pdata[chan])) {
			retval = 1;
			break;
		}
		Nl = wavelet[chan].Nbr_Ligne;
		Nc = wavelet[chan].Nbr_Col;
		pave_2d_extract_plan(wavelet[chan].Pave.Data, Imag, Nl, Nc, Plan);
		reget_rawdata(Imag, Nl, Nc, fit->pdata[chan]);
		wave_io_free(&wavelet[chan]);
	}

	/* Free */
	free(Imag);
	return retval;
}

int extract_plans(fits *fit, int Nbr_Plan, int Type) {
	int i;

	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);

	for (i = 0; i < Nbr_Plan; i++) {
		char filename[256], msg[256];

		g_snprintf(filename, sizeof(filename), "layer%02d", i);
		snprintf(msg, 256, _("Extracting %s..."), filename);
		set_progress_bar_data(msg, (float)i / Nbr_Plan);
		get_wavelet_layers(fit, Nbr_Plan, i, Type, -1);
		if (savefits(filename, fit))
			break;
	}
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_DONE);
	return 0;
}

/*****************************************************************************
 *                      M E D I A N     F I L T E R                          *
 ****************************************************************************/

gboolean end_median_filter(gpointer p) {
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
	assert(args->ksize % 2 == 1 && args->ksize > 1);
	int i, x, y, xx, yy, layer, iter = 0;
	int nx = args->fit->rx;
	int ny = args->fit->ry;
	int radius = (args->ksize - 1) / 2;
	int ksize_squared = args->ksize * args->ksize;
	double norm = (double) get_normalized_value(args->fit);
	double cur = 0.0, total;
	assert(nx > 0 && ny > 0);

	struct timeval t_start, t_end;

	char *msg = siril_log_color_message(_("Median Filter: processing...\n"), "red");
	msg[strlen(msg) - 1] = '\0';
	set_progress_bar_data(msg, PROGRESS_RESET);
	gettimeofday(&t_start, NULL);

	WORD *data = calloc(ksize_squared, sizeof(WORD));
	if (data == NULL) {
		printf("median filter: error allocating data\n");
		siril_add_idle(end_median_filter, args);
		set_progress_bar_data(_("Median filter failed"), PROGRESS_DONE);
		return GINT_TO_POINTER(1);
	}

	do {
		for (layer = 0; layer < args->fit->naxes[2]; layer++) {
			/* FILL image upside-down */
			WORD **image = malloc(ny * sizeof(WORD *));
			if (image == NULL) {
				printf("median filter: error allocating data\n");
				siril_add_idle(end_median_filter, args);
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

static int fmul_layer(fits *a, int layer, float coeff) {
	WORD *buf;
	int i;

	if (coeff < 0.0)
		return 1;
	buf = a->pdata[layer];
	for (i = 0; i < a->rx * a->ry; ++i) {
		buf[i] = round_to_WORD(buf[i] * coeff);
	}
	invalidate_stats_from_fit(a);
	return 0;
}

/*****************************************************************************
 *      B A N D I N G      R E D U C T I O N      M A N A G E M E N T        *
 ****************************************************************************/

int banding_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_) {
	struct banding_data *banding_args = (struct banding_data *)args->user;
	return BandingEngine(fit, banding_args->sigma, banding_args->amount,
			banding_args->protect_highlights, banding_args->applyRotation);
}

void apply_banding_to_sequence(struct banding_data *banding_args) {
	struct generic_seq_args *args = malloc(sizeof(struct generic_seq_args));
	args->seq = &com.seq;
	args->partial_image = FALSE;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = com.seq.selnum;
	args->prepare_hook = ser_prepare_hook;
	args->finalize_hook = ser_finalize_hook;
	args->save_hook = NULL;
	args->image_hook = banding_image_hook;
	args->idle_function = NULL;
	args->stop_on_error = FALSE;
	args->description = _("Banding Reduction");
	args->has_output = TRUE;
	args->new_seq_prefix = banding_args->seqEntry;
	args->load_new_sequence = TRUE;
	args->force_ser_output = FALSE;
	args->user = banding_args;
	args->already_in_a_thread = FALSE;
	args->parallel = TRUE;

	banding_args->fit = NULL;	// not used here

	start_in_new_thread(generic_sequence_worker, args);
}

// idle function executed at the end of the BandingEngine processing
gboolean end_BandingEngine(gpointer p) {
	struct banding_data *args = (struct banding_data *) p;
	stop_processing_thread();// can it be done here in case there is no thread?
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	update_used_memory();
	free(args);
	return FALSE;
}

/*** Reduces Banding in Canon DSLR images.
 * This code come from CanonBandingReduction.js v0.9.1, a script of
 * PixInsight, originally written by Georg Viehoever and
 * distributed under the terms of the GNU General Public License ******/
gpointer BandingEngineThreaded(gpointer p) {
	struct banding_data *args = (struct banding_data *) p;
	struct timeval t_start, t_end;

	siril_log_color_message(_("Banding Reducing: processing...\n"), "red");
	gettimeofday(&t_start, NULL);

	int retval = BandingEngine(args->fit, args->sigma, args->amount, args->protect_highlights, args->applyRotation);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	siril_add_idle(end_BandingEngine, args);
	
	return GINT_TO_POINTER(retval);
}

int BandingEngine(fits *fit, double sigma, double amount, gboolean protect_highlights, gboolean applyRotation) {
	int chan, row, i, ret = 0;
	WORD *line, *fixline;
	double minimum = DBL_MAX, globalsigma = 0.0;
	fits *fiximage = NULL;
	double invsigma = 1.0 / sigma;

	if (applyRotation) {
		point center = {gfit.rx / 2.0, gfit.ry / 2.0};
		cvRotateImage(fit, center, 90.0, -1, OPENCV_LINEAR);
	}

	if (new_fit_image(&fiximage, fit->rx, fit->ry, fit->naxes[2]))
		return 1;

	for (chan = 0; chan < fit->naxes[2]; chan++) {
		imstats *stat = statistics(NULL, -1, fit, chan, NULL, STATS_BASIC | STATS_MAD);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		double background = stat->median;
		double *rowvalue = calloc(fit->ry, sizeof(double));
		if (rowvalue == NULL) {
			fprintf(stderr, "BandingEngine: error allocating data\n");
			free_stats(stat);
			return 1;
		}
		if (protect_highlights) {
			globalsigma = stat->mad * MAD_NORM;
		}
		free_stats(stat);
		for (row = 0; row < fit->ry; row++) {
			line = fit->pdata[chan] + row * fit->rx;
			WORD *cpyline = calloc(fit->rx, sizeof(WORD));
			if (cpyline == NULL) {
				fprintf(stderr, "BandingEngine: error allocating data\n");
				free(rowvalue);
				return 1;
			}
			memcpy(cpyline, line, fit->rx * sizeof(WORD));
			int n = fit->rx;
			double median;
			if (protect_highlights) {
				quicksort_s(cpyline, n);
				WORD reject = round_to_WORD(
						background + invsigma * globalsigma);
				for (i = fit->rx - 1; i >= 0; i--) {
					if (cpyline[i] < reject)
						break;
					n--;
				}
				median = gsl_stats_ushort_median_from_sorted_data(cpyline, 1, n);
			} else {
				median = round_to_WORD(quickmedian(cpyline, n));
			}

			rowvalue[row] = background - median;
			minimum = min(minimum, rowvalue[row]);
			free(cpyline);
		}
		for (row = 0; row < fit->ry; row++) {
			fixline = fiximage->pdata[chan] + row * fiximage->rx;
			for (i = 0; i < fit->rx; i++)
				fixline[i] = round_to_WORD(rowvalue[row] - minimum);
		}
		free(rowvalue);
	}
	for (chan = 0; chan < fit->naxes[2]; chan++)
		fmul_layer(fiximage, chan, amount);
	ret = imoper(fit, fiximage, OPER_ADD);

	invalidate_stats_from_fit(fit);
	clearfits(fiximage);
	if ((!ret) && applyRotation) {
		point center = {gfit.rx / 2.0, gfit.ry / 2.0};
		cvRotateImage(fit, center, -90.0, -1, OPENCV_LINEAR);
	}

	return ret;
}

/*****************************************************************************
 *       N O I S E     C O M P U T A T I O N      M A N A G E M E N T        *
 ****************************************************************************/

/* Based on Jean-Luc Starck and Fionn Murtagh (1998), Automatic Noise
 * Estimation from the Multiresolution Support, Publications of the
 * Royal Astronomical Society of the Pacific, vol. 110, pp. 193–199.
 * slow algorithm. For now it is replaced by faster one. BUT, we need to keep it
 * in case we need it -. */
int backgroundnoise(fits* fit, double sigma[]) {
	int layer, k;
	fits *waveimage = calloc(1, sizeof(fits));

	if (waveimage == NULL) {
		fprintf(stderr, "backgroundnoise: error allocating data\n");
		return 1;
	}

	copyfits(fit, waveimage, CP_ALLOC | CP_FORMAT | CP_COPYA, 0);
	cvComputeFinestScale(waveimage);

	for (layer = 0; layer < fit->naxes[2]; layer++) {
		double sigma0, mean, norm_val;
		double epsilon = 0.0;
		WORD lo, hi;
		WORD *buf = waveimage->pdata[layer];
		unsigned int i;
		unsigned int ndata = fit->rx * fit->ry;
		assert(ndata > 0);

		imstats *stat = statistics(NULL, -1, waveimage, layer, NULL, STATS_BASIC);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		sigma0 = stat->sigma;
		mean = stat->mean;
		norm_val = stat->normValue;
		free_stats(stat);

		WORD *array1 = calloc(ndata, sizeof(WORD));
		WORD *array2 = calloc(ndata, sizeof(WORD));
		if (array1 == NULL || array2 == NULL) {
			printf("backgroundnoise: Error allocating data\n");
			if (array1)
				free(array1);
			if (array2)
				free(array2);
			return 1;
		}
		WORD *set = array1, *subset = array2;
		memcpy(set, buf, ndata * sizeof(WORD));

		lo = round_to_WORD(LOW_BOUND * norm_val);
		hi = round_to_WORD(HIGH_BOUND * norm_val);

		sigma[layer] = sigma0;

		int n = 0;
		do {
			sigma0 = sigma[layer];
			for (i = 0, k = 0; i < ndata; i++) {
				if (set[i] >= lo && set[i] <= hi) {
					if (fabs(set[i] - mean) < 3.0 * sigma0) {
						subset[k++] = set[i];
					}
				}
			}
			ndata = k;
			sigma[layer] = gsl_stats_ushort_sd(subset, 1, ndata);
			set = subset;
			(set == array1) ? (subset = array2) : (subset = array1);
			if (ndata == 0) {
				free(array1);
				free(array2);
				siril_log_message(_("backgroundnoise: Error, no data computed\n"));
				sigma[layer] = 0.0;
				return 1;
			}
			n++;
			epsilon = fabs(sigma[layer] - sigma0) / sigma[layer];
		} while (epsilon > EPSILON && n < MAX_ITER);
		sigma[layer] *= SIGMA_PER_FWHM; // normalization
		sigma[layer] /= 0.974; // correct for 2% systematic bias
		if (n == MAX_ITER)
			siril_log_message(_("backgroundnoise: does not converge\n"));
		free(array1);
		free(array2);
	}
	clearfits(waveimage);
	invalidate_stats_from_fit(fit);

	return 0;
}

gboolean end_noise(gpointer p) {
	struct noise_data *args = (struct noise_data *) p;
	stop_processing_thread();
	set_cursor_waiting(FALSE);
	update_used_memory();
	if (args->verbose) {
		struct timeval t_end;
		gettimeofday(&t_end, NULL);
		show_time(args->t_start, t_end);
	}
	free(args);
	return FALSE;
}

gpointer noise(gpointer p) {
	struct noise_data *args = (struct noise_data *) p;
	int chan;
	args->retval = 0;

	if (args->verbose) {
		siril_log_color_message(_("Noise standard deviation: calculating...\n"),
				"red");
		gettimeofday(&args->t_start, NULL);
	}

	for (chan = 0; chan < args->fit->naxes[2]; chan++) {
		imstats *stat = statistics(NULL, -1, args->fit, chan, NULL, STATS_NOISE);
		if (!stat) {
			args->retval = 1;
			siril_log_message(_("Error: statistics computation failed.\n"));
			break;
		}
		args->bgnoise[chan] = stat->bgnoise;
		free_stats(stat);
	}

	if (!args->retval) {
		double norm = (double) get_normalized_value(args->fit);
		for (chan = 0; chan < args->fit->naxes[2]; chan++)
			siril_log_message(
					_("Background noise value (channel: #%d): %0.3lf (%.3e)\n"), chan,
					args->bgnoise[chan], args->bgnoise[chan] / norm);
	}

	int retval = args->retval;
	if (args->use_idle)
		siril_add_idle(end_noise, args);

	return GINT_TO_POINTER(retval);
}

gpointer LRdeconv(gpointer p) {
	struct RL_data *args = (struct RL_data *) p;
	struct timeval t_start, t_end;

	siril_log_color_message(_("Lucy-Richardson deconvolution: processing...\n"), "red");
	gettimeofday(&t_start, NULL);

	cvLucyRichardson(args->fit, args->sigma, args->iter);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	siril_add_idle(end_generic, args);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
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
