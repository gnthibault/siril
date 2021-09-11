/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2021 team free-astro (see more in AUTHORS file)
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
#include "gui/progress_and_log.h"

#define ACTIVATE_NULLCHECK_FLOAT 1

// copies the area of an image into the memory buffer data
static void select_area_float(fits *fit, float *data, int layer, rectangle *bounds) {
	int i, j, k = 0;

	float *from = fit->fpdata[layer] +
		(fit->ry - bounds->y - bounds->h) * fit->rx + bounds->x;
	int stridefrom = fit->rx - bounds->w;

	for (i = 0; i < bounds->h; ++i) {
		for (j = 0; j < bounds->w; ++j) {
			data[k++] = *from++;
		}
		from += stridefrom;
	}
}

// in this case, N is the number of frames, so int is fine
float siril_stats_float_sd(const float data[], const int N, float *m) {
	double accumulator = 0.0; // accumulating in double precision is important for accuracy
	for (int i = 0; i < N; ++i) {
		accumulator += data[i];
	}
	float mean = (float)(accumulator / N);
	accumulator = 0.0;
	for (int i = 0; i < N; ++i)
		accumulator += (data[i] - mean) * (data[i] - mean);

	if (m) *m = mean;

	return sqrtf((float)(accumulator / (N - 1)));
}

/* For a univariate data set X1, X2, ..., Xn, the MAD is defined as the median
 * of the absolute deviations from the data's median:
 *  MAD = median (| Xi − median(X) |)
 */
double siril_stats_float_mad(const float *data, const size_t n, const double m, gboolean multithread, float *buffer) {
	double mad;
	const float median = (float)m;
	float *tmp = buffer ? buffer : malloc(n * sizeof(float));
	if (!tmp) {
		PRINT_ALLOC_ERR;
		return 0.0;
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) if(multithread && n > 10000) schedule(static)
#endif
	for (size_t i = 0; i < n; i++) {
		tmp[i] = fabsf(data[i] - median);
	}

	mad = histogram_median_float(tmp, n, multithread);
	if (!buffer) {
	    free(tmp);
	}
	return mad;
}

static double siril_stats_float_bwmv(const float* data, const size_t n,
		const float mad, const float median, gboolean multithread) {

	double bwmv = 0.0;
	double up = 0.0, down = 0.0;

	if (mad > 0.f) {
        const float factor = 1.f / (9.f * mad);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) if(multithread) schedule(static) reduction(+:up,down)
#endif
		for (size_t i = 0; i < n; i++) {
			const float i_med = data[i] - median;

			const float yi = i_med * factor;
			const float yi2 = fabsf(yi) < 1.f ? yi * yi : 1.f;

			up += SQR(i_med * SQR (1 - yi2));
			down += (1 - yi2) * (1 - 5 * yi2);
		}

		bwmv = n * (up / (down * down));
	}

	return bwmv;
}

float siril_stats_trmean_from_sorted_data(const float trim,
		const float sorted_data[], const size_t stride, const size_t size) {
	if (trim >= 0.5f) {
		return (float) gsl_stats_float_median_from_sorted_data(sorted_data, stride, size);
	} else {
		size_t ilow = (size_t) floorf(trim * size);
		size_t ihigh = size - ilow - 1;
		float mean = 0.f;
		float k = 0.f;

		/* compute mean of middle samples in [ilow,ihigh] */
		for (size_t i = ilow; i <= ihigh; ++i) {
			float delta = sorted_data[i * stride] - mean;
			k += 1.f;
			mean += delta / k;
		}
		return mean;
	}
}

#if 0
int IKSS(float *data, size_t n, double *location, double *scale, gboolean multithread) {
	size_t i, j;
	double mad, s, s0, m;
	float xlow, xhigh;

	quicksort_f(data, n);	// this sort is mandatory
	i = 0;
	j = n;
	s0 = 1;
	float *buffer = malloc(n * sizeof(float));
	if (!buffer) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	for (;;) {
		if (j - i < 1) {
			*location = *scale = 0;
			break;
		}
		m = gsl_stats_float_median_from_sorted_data(data + i, 1, j - i);
		mad = siril_stats_float_mad(data + i, j - i, m, multithread, buffer);
		if (mad == 0.0f) {
			free(buffer);
			return 1;
		}
		s = sqrt(siril_stats_float_bwmv(data + i, j - i, mad, m, multithread));
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
	free(buffer);
	return 0;
}
#endif

int IKSSlite(float *data, size_t n, const float median, float mad, double *location, double *scale, gboolean multithread) {
	size_t i, kept;
	float xlow, xhigh;

	xlow = median - 6.0 * mad;
	xhigh = median + 6.0 * mad;

	kept = 0;
	// removing pixels outside of +/- 6mad from median
	for (i = 0; i < n; i++) {
		if ((data[i] >= xlow) && (data[i] <= xhigh)) {
			if (i != kept) {
				data[kept] = data[i];
			}
			kept++;
		}
	}
	if (kept == 0)
		return 1;

	*location = histogram_median_float(data, kept, multithread);
	mad = siril_stats_float_mad(data, kept, *location, multithread, NULL);
	if (mad == 0.0f)
		return 1;

	*scale = sqrt(siril_stats_float_bwmv(data, kept, mad, *location, multithread)) *.991;
	/* 0.991 factor is to keep consistency with IKSS scale */
	return 0;
}

static float* reassign_to_non_null_data_float(float *data, size_t inputlen, size_t outputlen, int free_input) {
	size_t i, j = 0;
	float *ndata = malloc(outputlen * sizeof(float));
	if (!ndata) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	for (i = 0; i < inputlen; i++) {
		if (data[i] != 0.f) {
			if (j >= outputlen) {
				fprintf(stderr, "\n- stats MISMATCH in sizes (in: %zu, out: %zu), THIS IS A BUG: seqfile is wrong *********\n\n", inputlen, outputlen);
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
		const float data[], const size_t n, gboolean multithread) {
	/* finds the smallest and largest members of a dataset */

	if (n > 0 && data) {
		float min = data[0];
		float max = data[0];
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) if(multithread && n > 10000) reduction(max:max) reduction(min:min)
#endif
		for (size_t i = 0; i < n; i++) {
			const float xi = data[i];
			min = min(min, xi);
			max = max(max, xi);
		}

		*min_out = min;
		*max_out = max;
	}
}


/* this function tries to get the requested stats from the passed stats,
 * computes them and stores them in it if they have not already been */
imstats* statistics_internal_float(fits *fit, int layer, rectangle *selection, int option, imstats *stats, int bitpix, gboolean multithread) {
	int nx, ny;
	float *data = NULL;
	int stat_is_local = 0, free_data = 0;
	imstats* stat = stats;
	// median is included in STATS_BASIC but required to compute other data
	int compute_median = (option & STATS_BASIC) || (option & STATS_AVGDEV) ||
		(option & STATS_MAD) || (option & STATS_BWMV) || (option & STATS_IKSS);

	if (!stat) {
		allocate_stats(&stat);
		if (!stat) return NULL;
		stat_is_local = 1;
	} else {
		atomic_int_incref(stat->_nb_refs);
	}

	if (fit) {
		if (selection && selection->h > 0 && selection->w > 0) {
			nx = selection->w;
			ny = selection->h;
			data = malloc(nx * ny * sizeof(float));
			if (!data) {
				PRINT_ALLOC_ERR;
				if (stat_is_local) free(stat);
				return NULL;
			}
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

	/* We need to convert stats to the bitdepth of the orginal images
	 * to make sure they are logged in a consistent manner and no mixing 16b/32b can occur in the seq file
	 * - if from fit, the original fit bitpix is passed as argument
	 * - if from cache, so from seq file, we take it from seq->bitpix and pass it as argument
	 * This is really important because we can compute float statistics from ushort images.
	 * Fixes https://gitlab.com/free-astro/siril/-/issues/700
	 */

	if (stat->normValue == NULL_STATS) {
		if (bitpix == FLOAT_IMG) stat->normValue = 1.f;
		else stat->normValue = (bitpix == BYTE_IMG) ? UCHAR_MAX_DOUBLE : USHRT_MAX_DOUBLE;
	}

	/* Calculation of min and max */
	if ((option & (STATS_MINMAX | STATS_BASIC)) && (stat->min == NULL_STATS || stat->max == NULL_STATS)) {
		float min = 0, max = 0;
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing minmax\n", stat, fit, layer);
		siril_stats_float_minmax(&min, &max, data, stat->total, multithread);
		stat->min = (double)min * stat->normValue;
		stat->max = (double)max * stat->normValue;
	}

	/* Calculation of ngoodpix, mean, sigma and background noise */
	if ((option & (STATS_SIGMEAN | STATS_BASIC)) && (stat->ngoodpix <= 0L || stat->mean == NULL_STATS ||
			stat->sigma == NULL_STATS || stat->bgnoise == NULL_STATS)) {
		int status = 0;
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing basic\n", stat, fit, layer);
		siril_fits_img_stats_float(data, nx, ny, ACTIVATE_NULLCHECK_FLOAT, 0.0f, &stat->ngoodpix,
				NULL, NULL, &stat->mean, &stat->sigma, &stat->bgnoise,
				NULL, NULL, NULL, multithread, &status);
		if (status) {
			if (free_data) free(data);
			if (stat_is_local) free(stat);
			siril_log_message("fits_img_stats_float failed\n");
			return NULL;
		}
		stat->mean *= stat->normValue;
		stat->sigma *= stat->normValue;
		stat->bgnoise *= stat->normValue;
	}
	if (stat->ngoodpix == 0L) {
		if (free_data) free(data);
		if (stat_is_local) free(stat);
		return NULL;
	}


	/* we exclude 0 if some computations remain to be done or copy data if
	 * median has to be computed */
	if (ACTIVATE_NULLCHECK_FLOAT && fit && compute_median && stat->total != stat->ngoodpix) {
		data = reassign_to_non_null_data_float(data, stat->total, stat->ngoodpix, free_data);
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;
		}
		free_data = 1;
	}


	/* Calculation of median */
	if (compute_median && stat->median == NULL_STATS) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing median\n", stat, fit, layer);
		stat->median = histogram_median_float(data, stat->ngoodpix, multithread) * stat->normValue;
	}

	/* Calculation of average absolute deviation from the median */
	if ((option & STATS_AVGDEV) && stat->avgDev == NULL_STATS) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing absdev\n", stat, fit, layer);
		stat->avgDev = gsl_stats_float_absdev_m(data, 1, stat->ngoodpix, stat->median / stat->normValue) * stat->normValue;
	}

	/* Calculation of median absolute deviation */
	if (((option & STATS_MAD) || (option & STATS_BWMV) || (option & STATS_IKSS)) && stat->mad == NULL_STATS) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing mad\n", stat, fit, layer);
		stat->mad = siril_stats_float_mad(data, stat->ngoodpix, stat->median / stat->normValue, multithread, NULL) * stat->normValue;
	}

	/* Calculation of Bidweight Midvariance */
	if ((option & STATS_BWMV) && stat->sqrtbwmv == NULL_STATS) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing bimid\n", stat, fit, layer);
		double bwmv = siril_stats_float_bwmv(data, stat->ngoodpix, stat->mad / stat->normValue, stat->median / stat->normValue, multithread);
		stat->sqrtbwmv = sqrt(bwmv) * stat->normValue;
	}

	/* Calculation of IKSS. Only used for stacking normalization */
	if ((option & STATS_IKSS) && (stat->location == NULL_STATS || stat->scale == NULL_STATS)) {
		if (!data) {
			if (stat_is_local) free(stat);
			return NULL;	// not in cache, don't compute
		}
		siril_debug_print("- stats %p fit %p (%d): computing ikss\n", stat, fit, layer);
		if (IKSSlite(data, stat->ngoodpix, stat->median / stat->normValue, stat->mad / stat->normValue, &stat->location, &stat->scale, multithread)) {
			if (stat_is_local) free(stat);
			if (free_data) free(data);
			return
			 NULL;
		}
		stat->location *= stat->normValue;
		stat->scale *= stat->normValue;
	}

	if (free_data) free(data);
	return stat;
}

int compute_means_from_flat_cfa_float(fits *fit, float mean[4]) {
	int row, col, c, i = 0;
	float *data;
	unsigned int width, height;
	unsigned int startx, starty;

	mean[0] = mean[1] = mean[2] = mean[3] = 0.0;
	data = fit->fdata;
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
			mean[0] += data[col + row * width];
			mean[1] += data[1 + col + row * width];
			mean[2] += data[col + (1 + row) * width];
			mean[3] += data[1 + col + (1 + row) * width];
			i++;
		}
	}

	for (c = 0; c < 4; c++) {
		mean[c] /= (float) i;
	}
	return 0;
}
