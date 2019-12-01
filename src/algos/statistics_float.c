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
#include "gui/dialogs.h"
#include "sorting.h"
#include "statistics.h"

// copies the area of an image into the memory buffer data
static void select_area_float(fits *fit, float *data, int layer, rectangle *bounds) {
	int i, j, k = 0;

	float *from = fit->fpdata[layer] + (fit->ry - bounds->y - bounds->h) * fit->rx
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
static double siril_stats_float_mad(const float *data, const size_t stride,
		const size_t n, const double m) {
	size_t i;
	double mad;
	const float median = (float)m;
	float *tmp = calloc(n, sizeof(float));
	if (!tmp) {
		PRINT_ALLOC_ERR;
		return 0.0f;	// TODO: check return value
	}

#pragma omp parallel for num_threads(com.max_thread) if(n > 10000) private(i) schedule(static)
	for (i = 0; i < n; i++) {
		float delta = data[i * stride] - median;
		tmp[i] = fabsf(delta);
	}

	mad = histogram_median_float(tmp, n);
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

static double siril_stats_float_bwmv(const float* data, const size_t n,
		const double mad, const double median) {

	double bwmv = 0.0;
	double up = 0.0, down = 0.0;
	size_t i;

	if (mad > 0.0) {
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static) reduction(+:up,down)
		for (i = 0; i < n; i++) {
			double yi, ai, yi2;
			double i_med = (double)data[i] - median;

			yi = i_med / (9.0 * mad);
			yi2 = yi * yi;
			ai = (fabs(yi) < 1.0) ? 1.0 : 0.0;

			up += ai * SQR(i_med) * SQR(SQR (1 - yi2));
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

			yi = (data[i] - median) / (9.0 * mad);
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

static float* reassign_to_non_null_data_float(float *data, long inputlen, long outputlen, int free_input) {
	int i, j = 0;
	float *ndata = malloc(outputlen * sizeof(float));
	if (!ndata) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	for (i = 0; i < inputlen; i++) {
		if (data[i] > 0.0) {
			if (j >= outputlen) {
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

static void siril_stats_float_minmax(float *min_out, float *max_out,
		const float data[], const size_t stride, const size_t n) {
	/* finds the smallest and largest members of a dataset */

	if (n > 0 && data) {
		float min = data[0];
		float max = data[0];
		size_t i;

#pragma omp parallel for num_threads(com.max_thread) schedule(static) if(n > 10000) reduction(max:max) reduction(min:min)
		for (i = 0; i < n; i++) {
			float xi = data[i * stride];
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
imstats* statistics_internal_float(fits *fit, int layer, rectangle *selection, int option, imstats *stats) {
	int nx, ny;
	float *data = NULL;
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
			data = calloc(nx * ny, sizeof(float));
			select_area_float(fit, data, layer, selection);
			free_data = 1;
		} else {
			nx = fit->rx;
			ny = fit->ry;
			data = fit->fpdata[layer];
		}
		stat->total = nx * ny;
		if (stat->total == 0L) {
			if (stat_is_local) free(stat);
			return NULL;
		}
	}

	/* Calculation of min and max */
	if ((option & (STATS_MINMAX | STATS_BASIC)) && (stat->min < 0. || stat->max < 0.)) {
		float min = 0, max = 0;
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing minmax\n", stat, fit, layer);
		siril_stats_float_minmax(&min, &max, data, 1, stat->total);
		stat->min = (double)min;
		stat->max = (double)max;
		stat->normValue = 1.0;
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
		fits_img_stats_float(data, nx, ny, 1, 0.0f, &stat->ngoodpix,
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
		data = reassign_to_non_null_data_float(data, stat->total, stat->ngoodpix, free_data);
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
		stat->median = histogram_median_float(data, stat->ngoodpix);
	}

	/* Calculation of average absolute deviation from the median */
	if ((option & STATS_AVGDEV) && stat->avgDev < 0.) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing absdev\n", stat, fit, layer);
		stat->avgDev = gsl_stats_float_absdev_m(data, 1, stat->ngoodpix, stat->median);
	}

	/* Calculation of median absolute deviation */
	if (((option & STATS_MAD) || (option & STATS_BWMV)) && stat->mad < 0.) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing mad\n", stat, fit, layer);
		stat->mad = siril_stats_float_mad(data, 1, stat->ngoodpix, stat->median);
	}

	/* Calculation of Bidweight Midvariance */
	if ((option & STATS_BWMV) && stat->sqrtbwmv < 0.) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing bimid\n", stat, fit, layer);
		double bwmv = siril_stats_float_bwmv(data, stat->ngoodpix, stat->mad, stat->median);
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
			newdata[i] = (double)data[i];
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
