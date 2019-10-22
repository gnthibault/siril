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
#include "io/sequence.h"
#include "io/single_image.h"
#include "algos/colors.h"
#include "algos/Def_Math.h"
#include "algos/Def_Wavelet.h"
#include "algos/statistics.h"
#include "opencv/opencv.h"

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
		double dbuf, dgbuf;
		for (i = 0; i < a->rx * a->ry; ++i) {
			dbuf = (double) buf[i];
			dgbuf = (double) gbuf[i];
			switch (oper) {
			case OPER_ADD:
				gbuf[i] = round_to_WORD(dgbuf + dbuf);
				break;
			case OPER_SUB:
				gbuf[i] = round_to_WORD(dgbuf - dbuf);
				break;
			case OPER_MUL:
				gbuf[i] = round_to_WORD(dgbuf * dbuf);
				break;
			case OPER_DIV:
				gbuf[i] = (buf[i] == 0) ? 0 : round_to_WORD(dgbuf / dbuf);
				break;
			}
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
			if (temp > USHRT_MAX_DOUBLE) {
				siril_debug_print("OVERFLOW in FDIV: %lf\n", temp);
				retvalue = 1;
			}
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
	if (div == NULL) {
		PRINT_ALLOC_ERR;
		return 1;
	}

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
	if (Imag == NULL) {
		PRINT_ALLOC_ERR;
		return 1;
	}

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
