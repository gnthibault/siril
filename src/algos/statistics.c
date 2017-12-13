/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2017 team free-astro (see more in AUTHORS file)
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <gsl/gsl_statistics.h>
#include "core/siril.h"
#include "core/proto.h"

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

/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */

#define ELEM_SWAP_WORD(a,b) { register WORD t=(a);(a)=(b);(b)=t; }

static double siril_stats_ushort_median(WORD *arr, int n) {
	int low, high;
	int median;
	int middle, ll, hh;

	low = 0;
	high = n - 1;
	median = (low + high) / 2;
	for (;;) {
		if (high <= low) /* One element only */
			return (double) arr[median];

		if (high == low + 1) { /* Two elements only */
			if (arr[low] > arr[high])
				ELEM_SWAP_WORD(arr[low], arr[high]);
			return (double) arr[median];
		}

		/* Find median of low, middle and high items; swap into position low */
		middle = (low + high) / 2;
		if (arr[middle] > arr[high])
			ELEM_SWAP_WORD(arr[middle], arr[high]);
		if (arr[low] > arr[high])
			ELEM_SWAP_WORD(arr[low], arr[high]);
		if (arr[middle] > arr[low])
			ELEM_SWAP_WORD(arr[middle], arr[low]);

		/* Swap low item (now in position middle) into position (low+1) */
		ELEM_SWAP_WORD(arr[middle], arr[low + 1]);

		/* Nibble from each end towards middle, swapping items when stuck */
		ll = low + 1;
		hh = high;
		for (;;) {
			do
				ll++;
			while (arr[low] > arr[ll]);
			do
				hh--;
			while (arr[hh] > arr[low]);

			if (hh < ll)
				break;

			ELEM_SWAP_WORD(arr[ll], arr[hh]);
		}

		/* Swap middle item (in position low) back into correct position */
		ELEM_SWAP_WORD(arr[low], arr[hh]);

		/* Re-set active partition */
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
	}
	return -1;
}

#undef ELEM_SWAP_WORD

#define ELEM_SWAP_DOUBLE(a,b) { register double t=(a);(a)=(b);(b)=t; }

static double siril_stats_double_median(double *arr, int n) {
	int low, high;
	int median;
	int middle, ll, hh;

	low = 0;
	high = n - 1;
	median = (low + high) / 2;
	for (;;) {
		if (high <= low) /* One element only */
			return (double) arr[median];

		if (high == low + 1) { /* Two elements only */
			if (arr[low] > arr[high])
				ELEM_SWAP_DOUBLE(arr[low], arr[high]);
			return (double) arr[median];
		}

		/* Find median of low, middle and high items; swap into position low */
		middle = (low + high) / 2;
		if (arr[middle] > arr[high])
			ELEM_SWAP_DOUBLE(arr[middle], arr[high]);
		if (arr[low] > arr[high])
			ELEM_SWAP_DOUBLE(arr[low], arr[high]);
		if (arr[middle] > arr[low])
			ELEM_SWAP_DOUBLE(arr[middle], arr[low]);

		/* Swap low item (now in position middle) into position (low+1) */
		ELEM_SWAP_DOUBLE(arr[middle], arr[low + 1]);

		/* Nibble from each end towards middle, swapping items when stuck */
		ll = low + 1;
		hh = high;
		for (;;) {
			do
				ll++;
			while (arr[low] > arr[ll]);
			do
				hh--;
			while (arr[hh] > arr[low]);

			if (hh < ll)
				break;

			ELEM_SWAP_DOUBLE(arr[ll], arr[hh]);
		}

		/* Swap middle item (in position low) back into correct position */
		ELEM_SWAP_DOUBLE(arr[low], arr[hh]);

		/* Re-set active partition */
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
	}
	return -1;
}

#undef ELEM_SWAP_DOUBLE

/* For a univariate data set X1, X2, ..., Xn, the MAD is defined as the median
 * of the absolute deviations from the data's median:
 *  MAD = median (| Xi − median(X) |)
 */

static double siril_stats_ushort_mad(WORD* data, const size_t stride,
		const size_t n, const double m) {
	size_t i;
	double median;
	WORD *tmp;

	tmp = calloc(n, sizeof(WORD));

#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
	for (i = 0; i < n; i++) {
		const WORD delta = fabs(data[i * stride] - m);
		tmp[i] = delta;
	}
	median = siril_stats_ushort_median(tmp, n);
	free(tmp);

	return median;
}

static double siril_stats_double_mad(const double* data, const size_t stride,
		const size_t n, const double m) {
	size_t i;
	double median;
	double *tmp;

	tmp = calloc(n, sizeof(double));

#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
	for (i = 0; i < n; i++) {
		const double delta = fabs(data[i * stride] - m);
		tmp[i] = delta;
	}
	median = siril_stats_double_median(tmp, n);
	free(tmp);

	return median;
}

static double siril_stats_ushort_bwmv(const WORD* data, const size_t n,
		const double mad, const double median) {

	double bwmv = 0.0;
	double up = 0.0, down = 0.0;
	size_t i;

	if (mad > 0.0) {
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

	quicksort_d(data, n);
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

static WORD* reassign_data(WORD *data, int nb, long ngoodpix) {
	int i, j = 0;
	WORD *ndata = calloc(ngoodpix, sizeof(WORD));

	for (i = 0; i < nb; i++) {
		if (data[i] > 0) {
			ndata[j] = data[i];
			j++;
		}
	}
	free(data);
	return ndata;
}

/* computes statistics on the given layer of the given image. Mean, sigma and noise
 * are computed with a cfitsio function rewritten here.
 * min and max value, average deviation, MAD, Bidweight Midvariance and IKSS are computed with gsl stats.
 */
imstats* statistics(fits *fit, int layer, rectangle *selection, int option, int nullcheck) {
	double mean = 0.0;
	double median = 0.0;
	double sigma = 0.0;
	double noise = 0.0;
	long ngoodpix = 0L;
	double avgDev = 0.0;
	double mad = 0.0;
	double bwmv = 0.0;
	double location = 0.0, scale = 0.0;
	int status = 0;
	int nx, ny;
	WORD min = 0, max = 0, *data;
	size_t i, norm = USHRT_MAX;
	imstats* stat = NULL;

	if (selection && selection->h > 0 && selection->w > 0) {
		nx = selection->w;
		ny = selection->h;
		data = calloc(nx * ny, sizeof(WORD));
		select_area(fit, data, layer, selection);
	} else {
		nx = fit->rx;
		ny = fit->ry;
		data = calloc(nx * ny, sizeof(WORD));
		memcpy(data, fit->pdata[layer], nx * ny * sizeof(WORD));
	}

	/* Calculation of mean, sigma and noise */
	if (option & STATS_BASIC)
		fits_img_stats_ushort(data, nx, ny, nullcheck, 0, &ngoodpix, NULL, NULL, &mean,
			&sigma, &noise, NULL, NULL, NULL, &status);
	if (status) {
		free(data);
		return NULL;
	}
	if (ngoodpix == 0L) {
		free(data);
		return NULL;
	}

	/* we exclude 0 if nullcheck */
	if (nullcheck) {
		data = reassign_data(data, nx * ny, ngoodpix);
	}

	/* Calculation of median */
	if ((option & STATS_BASIC) || (option & STATS_AVGDEV)
			|| (option & STATS_MAD) || (option & STATS_BWMV))
		median = siril_stats_ushort_median(data, ngoodpix);

	if (option & STATS_BASIC) {
		gsl_stats_ushort_minmax(&min, &max, data, 1, ngoodpix);
		if (max <= UCHAR_MAX)
			norm = UCHAR_MAX;
		else
			norm = USHRT_MAX;
	}

	/* Calculation of average absolute deviation from the median */
	if (option & STATS_AVGDEV)
		avgDev = gsl_stats_ushort_absdev_m(data, 1, ngoodpix, median);

	/* Calculation of median absolute deviation */
	if ((option & STATS_MAD) || (option & STATS_BWMV))
		mad = siril_stats_ushort_mad(data, 1, ngoodpix, median);

	/* Calculation of Bidweight Midvariance */
	if (option & STATS_BWMV)
		bwmv = siril_stats_ushort_bwmv(data, ngoodpix, mad, median);

	/* Calculation of IKSS. Used for stacking */
	if (option & STATS_IKSS) {
		double *newdata = calloc(ngoodpix, sizeof(double));

		/* we convert in the [0, 1] range */
		for (i = 0; i < ngoodpix; i++) {
			newdata[i] = (double) data[i] / ((double) norm);
		}
		IKSS(newdata, ngoodpix, &location, &scale);
		/* go back to the original range */
		location *= ((double) norm);
		scale *= ((double) norm);
		free(newdata);
	}

	stat = malloc(sizeof(imstats));

	switch (layer) {
	case 0:
		if (fit->naxes[2] == 1)
			strcpy(stat->layername, "B&W");
		else
			strcpy(stat->layername, "Red");
		break;
	case 1:
		strcpy(stat->layername, "Green");
		break;
	case 2:
		strcpy(stat->layername, "Blue");
		break;
	}

	stat->total = (nx * ny);
	stat->ngoodpix = ngoodpix;
	stat->mean = mean;
	stat->avgDev = avgDev;
	stat->mad = mad;
	stat->median = median;
	stat->sigma = sigma;
	stat->bgnoise = noise;
	stat->min = (double) min;
	stat->max = (double) max;
	stat->sqrtbwmv = sqrt(bwmv);
	stat->location = location;
	stat->scale = scale;
	stat->normValue = (double) norm;

	free(data);
	return stat;
}
