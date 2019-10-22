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

/* HOW STATISTICS WORK
 * Stats for an image are computed on request and are carried in an imstats
 * structure. Depending on the current mode, this structure is stored for
 * caching in the fits->stats or in the sequence->stats.
 * If it is stored in the fits, when it is disposed, it is copied in the
 * sequence if it belongs to one.
 * If it is stored in the sequence, when the sequence is disposed, it is saved
 * in the seqfile, in `M' fields.
 * When a sequence is loaded, previously computed stats are recovered that way.
 * All operations that need to access stats should do it with the statistics()
 * function at the bottom of this file which provides this abstraction.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <gsl/gsl_statistics.h>
#include "core/siril.h"
#include "core/proto.h"
#include "sorting.h"
#include "statistics.h"

static void stats_set_default_values(imstats *stat);

// copies the area of an image into the memory buffer data
static void select_area(fits *fit, WORD *data, int layer, rectangle *bounds) {
	int i, j, k = 0;

	WORD *from = fit->pdata[layer] + (fit->ry - bounds->y - bounds->h) * fit->rx
			+ bounds->x;
	int stridefrom = fit->rx - bounds->w;

	for (i = 0; i < bounds->h; ++i) {
		for (j = 0; j < bounds->w; ++j) {
			data[k] = *from++;
			k++;
		}
		from += stridefrom;
	}
}

/* For a univariate data set X1, X2, ..., Xn, the MAD is defined as the median
 * of the absolute deviations from the data's median:
 *  MAD = median (| Xi − median(X) |)
 */
static double siril_stats_ushort_mad(WORD* data, const size_t stride,
		const size_t n, const double m) {
	size_t i;
	double mad;
	int median = round_to_int(m);	// we use it on integer data anyway
	WORD *tmp = calloc(n, sizeof(WORD));

#pragma omp parallel for num_threads(com.max_thread) if(n > 10000) private(i) schedule(static)
	for (i = 0; i < n; i++) {
		int delta = data[i * stride] - median;
		tmp[i] = (WORD)abs(delta);
	}

	mad = histogram_median (tmp, n);
	free(tmp);
	return mad;
}

static double siril_stats_double_mad(const double* data, const size_t stride,
		const size_t n, const double median) {
	size_t i;
	double *tmp = calloc(n, sizeof(double));
	double mad;

#pragma omp parallel for num_threads(com.max_thread) if(n > 10000) private(i) schedule(static)
	for (i = 0; i < n; i++) {
		tmp[i] = fabs(data[i * stride] - median);
	}

	mad = histogram_median_double (tmp, n);
	free(tmp);
	return mad;
}

static double siril_stats_ushort_bwmv(const WORD* data, const size_t n,
		const double mad, const double median) {

	double bwmv = 0.0;
	double up = 0.0, down = 0.0;
	size_t i;

	if (mad > 0.0) {
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static) reduction(+:up,down)
		for (i = 0; i < n; i++) {
			double yi, ai, yi2;

			yi = ((double) data[i] - median) / (9 * mad);
			yi2 = yi * yi;
			ai = (fabs(yi) < 1.0) ? 1.0 : 0.0;

			up += ai * SQR((double ) data[i] - median) * SQR(SQR (1 - yi2));
			down += (ai * (1 - yi2) * (1 - 5 * yi2));
		}

		bwmv = n * (up / (down * down));
	}

	return bwmv;
}

static double siril_stats_double_bwmv(const double* data, const size_t n,
		const double mad, const double median) {

	double bwmv = 0.0;
	double up = 0.0, down = 0.0;
	size_t i;

	if (mad > 0.0) {
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static) reduction(+:up,down)
		for (i = 0; i < n; i++) {
			double yi, ai, yi2;

			yi = (data[i] - median) / (9 * mad);
			yi2 = yi * yi;
			ai = (fabs(yi) < 1.0) ? 1.0 : 0.0;

			up += ai * SQR(data[i] - median) * SQR(SQR (1 - yi2));
			down += (ai * (1 - yi2) * (1 - 5 * yi2));
		}
		bwmv = n * (up / (down * down));
	}

	return bwmv;
}

static int IKSS(double *data, int n, double *location, double *scale) {
	size_t i, j;
	double mad, s, s0, m, xlow, xhigh;

	quicksort_d(data, n);	// this sort is mandatory
	i = 0;
	j = n;
	s0 = 1;
	for (;;) {
		if (j - i < 1) {
			*location = *scale = 0;
			break;
		}
		m = gsl_stats_median_from_sorted_data(data + i, 1, j - i);
		mad = siril_stats_double_mad(data + i, 1, j - i, m);
		s = sqrt(siril_stats_double_bwmv(data + i, j - i, mad, m));
		if (s < 2E-23) {
			*location = m;
			*scale = 0;
			break;
		}
		if (((s0 - s) / s) < 10E-6) {
			*location = m;
			*scale = 0.991 * s;
			break;
		}
		s0 = s;
		xlow = m - 4 * s;
		xhigh = m + 4 * s;
		while (data[i] < xlow)
			i++;
		while (data[j - 1] > xhigh)
			j--;
	}
	return 0;
}

static WORD* reassign_to_non_null_data(WORD *data, long inputlen, long outputlen, int free_input) {
	int i, j = 0;
	WORD *ndata = malloc(outputlen * sizeof(WORD));
	if (!ndata) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	for (i = 0; i < inputlen; i++) {
		if (data[i] > 0) {
			if (j >= outputlen) {
				//assert(0);
				fprintf(stderr, "\n- stats MISMATCH in sizes (in: %ld, out: %ld), THIS IS A BUG: seqfile is wrong *********\n\n", inputlen, outputlen);
				break;
			}
			ndata[j] = data[i];
			j++;
		}
	}
	if (free_input)
		free(data);
	return ndata;
}

static void siril_stats_ushort_minmax(WORD *min_out, WORD *max_out,
		const WORD data[], const size_t stride, const size_t n) {
	/* finds the smallest and largest members of a dataset */

	if (n > 0 && data) {
		WORD min = data[0];
		WORD max = data[0];
		size_t i;

#pragma omp parallel for num_threads(com.max_thread) schedule(static) if(n > 10000) reduction(max:max) reduction(min:min)
		for (i = 0; i < n; i++) {
			WORD xi = data[i * stride];

			if (xi < min)
				min = xi;

			if (xi > max)
				max = xi;
		}

		*min_out = min;
		*max_out = max;
	}
}

/* this function tries to get the requested stats from the passed stats,
 * computes them and stores them in it if they have not already been */
static imstats* statistics_internal(fits *fit, int layer, rectangle *selection, int option, imstats *stats) {
	int nx, ny;
	WORD *data = NULL;
	int stat_is_local = 0, free_data = 0;
	imstats* stat = stats;
	// median is included in STATS_BASIC but required to compute other data
	int compute_median = (option & STATS_BASIC) || (option & STATS_AVGDEV) ||
		(option & STATS_MAD) || (option & STATS_BWMV);

	if (!stat) {
		allocate_stats(&stat);
		if (!stat) return NULL;
		stat_is_local = 1;
	}

	if (fit) {
		if (selection && selection->h > 0 && selection->w > 0) {
			nx = selection->w;
			ny = selection->h;
			data = calloc(nx * ny, sizeof(WORD));
			select_area(fit, data, layer, selection);
			free_data = 1;
		} else {
			nx = fit->rx;
			ny = fit->ry;
			data = fit->pdata[layer];
		}
		stat->total = nx * ny;
		if (stat->total == 0L) {
			if (stat_is_local) free(stat);
			return NULL;
		}
	}

	/* Calculation of min and max */
	if ((option & (STATS_MINMAX | STATS_BASIC)) && (stat->min < 0. || stat->max < 0.)) {
		// TODO 1.0: change to double
		WORD min = 0, max = 0, norm = 0;
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing minmax\n", stat, fit, layer);
		siril_stats_ushort_minmax(&min, &max, data, 1, stat->total);
		if (fit->bitpix == BYTE_IMG)
			norm = UCHAR_MAX;
		else norm = USHRT_MAX;
		stat->min = (double)min;
		stat->max = (double)max;
		stat->normValue = (double)norm;
	}

	/* Calculation of ngoodpix, mean, sigma and background noise */
	if ((option & (STATS_NOISE | STATS_BASIC)) && (stat->ngoodpix <= 0L || stat->mean < 0. ||
			stat->sigma < 0. || stat->bgnoise < 0.)) {
		int status = 0;
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing basic\n", stat, fit, layer);
		fits_img_stats_ushort(data, nx, ny, 1, 0, &stat->ngoodpix,
				NULL, NULL, &stat->mean, &stat->sigma, &stat->bgnoise,
				NULL, NULL, NULL, &status);
		if (status) {
			if (free_data) free(data);
			if (stat_is_local) free(stat);
			return NULL;
		}
	}
	if (stat->ngoodpix == 0L) {
		if (free_data) free(data);
		if (stat_is_local) free(stat);
		return NULL;
	}

	/* we exclude 0 if some computations remain to be done or copy data if
	 * median has to be computed */
	if (fit && (compute_median || ((option & STATS_IKSS) && stat->total != stat->ngoodpix))) {
		data = reassign_to_non_null_data(data, stat->total, stat->ngoodpix, free_data);
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;
		}
		free_data = 1;
	}

	/* Calculation of median */
	if (compute_median && stat->median < 0.) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing median\n", stat, fit, layer);
		stat->median = histogram_median(data, stat->ngoodpix);
	}

	/* Calculation of average absolute deviation from the median */
	if ((option & STATS_AVGDEV) && stat->avgDev < 0.) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing absdev\n", stat, fit, layer);
		stat->avgDev = gsl_stats_ushort_absdev_m(data, 1, stat->ngoodpix, stat->median);
	}

	/* Calculation of median absolute deviation */
	if (((option & STATS_MAD) || (option & STATS_BWMV)) && stat->mad < 0.) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing mad\n", stat, fit, layer);
		stat->mad = siril_stats_ushort_mad(data, 1, stat->ngoodpix, stat->median);
	}

	/* Calculation of Bidweight Midvariance */
	if ((option & STATS_BWMV) && stat->sqrtbwmv < 0.) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing bimid\n", stat, fit, layer);
		double bwmv = siril_stats_ushort_bwmv(data, stat->ngoodpix, stat->mad, stat->median);
		stat->sqrtbwmv = sqrt(bwmv);
	}

	/* Calculation of IKSS. Only used for stacking normalization */
	if ((option & STATS_IKSS) && (stat->location < 0. || stat->scale < 0.)) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing ikss\n", stat, fit, layer);
		long i;
		double *newdata = malloc(stat->ngoodpix * sizeof(double));
		if (!newdata) {
			if (stat_is_local) free(stat);
			if (free_data) free(data);
			PRINT_ALLOC_ERR;
			return NULL;
		}

		/* we convert in the [0, 1] range */
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
		for (i = 0; i < stat->ngoodpix; i++) {
			newdata[i] = (double)data[i] / stat->normValue;
		}
		IKSS(newdata, stat->ngoodpix, &stat->location, &stat->scale);
		/* go back to the original range */
		stat->location *= stat->normValue;
		stat->scale *= stat->normValue;
		free(newdata);
	}

	if (free_data) free(data);
	return stat;
}

/* Computes statistics on the given layer of the given opened image.
 * Mean, sigma and noise are computed with a cfitsio function rewritten here.
 * Min and max value, average deviation, MAD, Bidweight Midvariance and IKSS
 * are computed with gsl stats.
 *
 * If the selection is not null or empty, computed data is not stored and seq
 * is not used.
 * If seq is null (single image processing), image_index is ignored, data is
 * stored in the fit, which cannot be NULL.
 * If seq is non-null, fit can be null to check for cached data.
 * The return value, if non-null, may be freed only using the free_stats()
 * function because of the special rule of this object that has a reference
 * counter because it can be referenced in 3 different places.
 */
imstats* statistics(sequence *seq, int image_index, fits *fit, int layer, rectangle *selection, int option) {
	imstats *oldstat = NULL, *stat;
	if (selection && selection->h > 0 && selection->w > 0) {
		// we have a selection, don't store anything
		return statistics_internal(fit, layer, selection, option, NULL);
	} else if (!seq || image_index < 0) {
		// we have a single image, store in the fits
		if (fit->stats && fit->stats[layer]) {
			oldstat = fit->stats[layer];
			oldstat->_nb_refs++;
		}
		stat = statistics_internal(fit, layer, NULL, option, oldstat);
		if (!stat) {
			fprintf(stderr, "- stats failed for fit %p (%d)\n", fit, layer);
			if (oldstat) {
				stats_set_default_values(oldstat);
				oldstat->_nb_refs--;
			}
			return NULL;
		}
		if (!oldstat)
			add_stats_to_fit(fit, layer, stat);
		return stat;
	} else {
		// we have sequence data, store in the sequence
		if (seq->stats && seq->stats[layer]) {
			oldstat = seq->stats[layer][image_index];
			if (oldstat)	// can be NULL here
				oldstat->_nb_refs++;
		}
		stat = statistics_internal(fit, layer, NULL, option, oldstat);
		if (!stat) {
			if (fit)
				fprintf(stderr, "- stats failed for %d in seq (%d)\n",
						image_index, layer);
			if (oldstat) {
				stats_set_default_values(oldstat);
				oldstat->_nb_refs--;
			}
			return NULL;
		}
		if (!oldstat)
			add_stats_to_seq(seq, image_index, layer, stat);
		if (fit)
			add_stats_to_fit(fit, layer, stat);	// can be useful too
		return stat;
	}
}

int compute_means_from_flat_cfa(fits *fit, double mean[4]) {
	int row, col, c, i = 0;
	WORD *data;
	unsigned int width, height;
	unsigned int startx, starty;

	mean[0] = mean[1] = mean[2] = mean[3] = 0.0;
	data = fit->data;
	width = fit->rx;
	height = fit->ry;

	/* due to vignetting it is better to take an area in the
	 * center of the flat image
	 */
	startx = width / 3;
	starty = height / 3;

	/* make sure it does start at the beginning of CFA pattern */
	startx = startx % 2 == 0 ? startx : startx + 1;
	starty = starty % 2 == 0 ? starty : starty + 1;

	siril_debug_print("Computing stat in (%d, %d, %d, %d)\n", startx, starty,
			width - 1 - startx, height - 1 - starty);

	/* Compute mean of each channel */
	for (row = starty; row < height - 1 - starty; row += 2) {
		for (col = startx; col < width - 1 - startx; col += 2) {
			mean[0] += (double) data[col + row * width];
			mean[1] += (double) data[1 + col + row * width];
			mean[2] += (double) data[col + (1 + row) * width];
			mean[3] += (double) data[1 + col + (1 + row) * width];
			i++;
		}
	}

	for (c = 0; c < 4; c++) {
		mean[c] /= (double) i;
	}
	return 0;
}

/****************** statistics caching and data management *****************/

/* reference an imstats struct to a fits, creates the stats array if needed */
void add_stats_to_fit(fits *fit, int layer, imstats *stat) {
	if (!fit->stats)
		fit->stats = calloc(fit->naxes[2], sizeof(imstats *));
	if (fit->stats[layer]) {
		if (fit->stats[layer] != stat) {
			siril_debug_print("- stats %p in fit %p (%d) is being replaced\n", fit->stats[layer], fit, layer);
			free_stats(fit->stats[layer]);
		} else return;
	}
	fit->stats[layer] = stat;
	stat->_nb_refs++;
	siril_debug_print("- stats %p saved to fit %p (%d)\n", stat, fit, layer);
}

static void add_stats_to_stats(sequence *seq, int nb_layers, imstats ****stats, int image_index, int layer, imstats *stat) {
	if (!*stats)
		*stats = calloc(nb_layers, sizeof(imstats **));
	if (!(*stats)[layer])
		(*stats)[layer] = calloc(seq->number, sizeof(imstats *));
	if ((*stats)[layer][image_index]) {
		if ((*stats)[layer][image_index] != stat) {
			siril_debug_print("- stats %p, %d in seq (%d) is being replaced\n", (*stats)[layer][image_index], image_index, layer);
			free_stats((*stats)[layer][image_index]);
		} else return;
	}
	siril_debug_print("- stats %p, %d in seq (%d): saving data\n", stat, image_index, layer);
	(*stats)[layer][image_index] = stat;
	seq->needs_saving = TRUE;
	stat->_nb_refs++;
}

void add_stats_to_seq(sequence *seq, int image_index, int layer, imstats *stat) {
	add_stats_to_stats(seq, seq->nb_layers, &seq->stats, image_index, layer, stat);
}

void add_stats_to_seq_backup(sequence *seq, int image_index, int layer, imstats *stat) {
	add_stats_to_stats(seq, 3, &seq->stats_bkp, image_index, layer, stat);
}

/* saves cached stats from the fits to its sequence, and clears the cache of the fits */
void save_stats_from_fit(fits *fit, sequence *seq, int index) {
	int layer;
	if (!fit || !fit->stats || !seq || index < 0) return;
	for (layer = 0; layer < fit->naxes[2]; layer++) {
		if (fit->stats[layer])
			add_stats_to_seq(seq, index, layer, fit->stats[layer]);
		free_stats(fit->stats[layer]);
		fit->stats[layer] = NULL;
	}
}

/* fit must be already read from disk or have naxes set at least */
void copy_seq_stats_to_fit(sequence *seq, int index, fits *fit) {
	if (seq->stats) {
		int layer;
		for (layer = 0; layer < fit->naxes[2]; layer++) {
			if (seq->stats[layer] && seq->stats[layer][index]) {
				add_stats_to_fit(fit, layer, seq->stats[layer][index]);
				siril_debug_print("- stats %p, copied from seq %d (%d)\n", fit->stats[layer], index, layer);
			}
		}
	}
}

/* if image data has changed, use this to force recomputation of the stats */
void invalidate_stats_from_fit(fits *fit) {
	if (fit->stats) {
		int layer;
		for (layer = 0; layer < fit->naxes[2]; layer++) {
			siril_debug_print("- stats %p cleared from fit (%d)\n", fit->stats[layer], layer);
			free_stats(fit->stats[layer]);
			fit->stats[layer] = NULL;
		}
	}
}

/* if image data and image structure has changed, invalidate the complete stats data structure */
void full_stats_invalidation_from_fit(fits *fit) {
	if (fit->stats) {
		invalidate_stats_from_fit(fit);
		free(fit->stats);
		fit->stats = NULL;
	}
}

static void stats_set_default_values(imstats *stat) {
	stat->total = -1L;
	stat->ngoodpix = -1L;
	stat->mean = stat->avgDev = stat->median = stat->sigma = stat->bgnoise = stat->min = stat->max = stat->normValue = stat->mad = stat->sqrtbwmv = stat->location = stat->scale = -1.0;
}

/* allocates an imstat structure and initializes it with default values that
 * are used by the statistics() function.
 * Only use free_stats() to free the return value.
 * Increment the _nb_refs if a new reference to the struct's address is made. */
void allocate_stats(imstats **stat) {
	if (stat) {
		if (!*stat)
			*stat = malloc(sizeof(imstats));
		if (!*stat) { PRINT_ALLOC_ERR; return; } // OOM
		stats_set_default_values(*stat);
		(*stat)->_nb_refs = 1;
		siril_debug_print("- stats %p allocated\n", *stat);
	}
}

/* frees an imstats struct if there are no more references to it.
 * returns NULL if it was freed, the argument otherwise. */
imstats* free_stats(imstats *stat) {
	if (stat && stat->_nb_refs-- == 1) {
		siril_debug_print("- stats %p has no more refs, freed\n", stat);
		free(stat);
		return NULL;
	}
	if (stat)
		siril_debug_print("- stats %p has refs (%d)\n", stat, stat->_nb_refs);
	return stat;
}

/* calls free_stats on all stats of a sequence */
void clear_stats(sequence *seq, int layer) {
	if (seq->stats && seq->stats[layer]) {
		int i;
		for (i = 0; i < seq->number; i++) {
			if (seq->stats[layer][i]) {
				siril_debug_print("- stats %p freed from seq %d (%d)\n", seq->stats[layer][i], i, layer);
				free_stats(seq->stats[layer][i]);
				seq->stats[layer][i] = NULL;
			}
		}
	}
}
