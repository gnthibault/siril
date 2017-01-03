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
#include "gui/histogram.h"

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

static double siril_stats_ushort_median(gsl_histogram *histo, const size_t n, int nullcheck) {
	size_t i;
	size_t hist_size = gsl_histogram_bins(histo);
	double sum = 0.0;
	double median = 0.0;
	int zero = (nullcheck == 0) ? zero = 0 : 1;

	/* Get the median value without 0 value*/
	for (i = zero; i < hist_size; i++) {
		sum += gsl_histogram_get(histo, i);
		if (sum > ((double) n * 0.5)) {
			median = (double) i;
			break;	//we get out of the loop
		}
	}
	return median;
}

static double siril_stats_ushort_mad(const WORD* data, const size_t stride,
		const size_t n, const double m, int nullcheck) {
	size_t i;
	double median;
	gsl_histogram *histo;

	histo = gsl_histogram_alloc(USHRT_MAX + 1);
	gsl_histogram_set_ranges_uniform(histo, 0, USHRT_MAX);
	for (i = 0; i < n; i++) {
		const double delta = fabs(data[i * stride] - m);
		gsl_histogram_increment(histo, delta);
	}
	median = siril_stats_ushort_median(histo, n, nullcheck);
	gsl_histogram_free(histo);

	return median;
}

static double siril_stats_double_mad(const double* data, const size_t stride,
		const size_t n, const double m) {
	size_t i;
	double median;
	double *tmp;

	tmp = calloc(n, sizeof(double));

	for (i = 0; i < n; i++) {
		const double delta = fabs(data[i * stride] - m);
		tmp[i] = delta;
	}
	quicksort_d(tmp, n);
	median = gsl_stats_median_from_sorted_data(tmp, 1, n);
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

	for (i = 0; i < nb; i++)
		if (data[i] > 0) {
			ndata[j] = data[i];
			j++;
		}
	free(data);
	return ndata;
}

/* computes statistics on the given layer of the given image. It creates the
 * histogram to easily extract the median. mean, sigma and noise
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
	gsl_histogram* histo;
	size_t i, hist_size;
	imstats* stat = NULL;

	if (selection && selection->h > 0 && selection->w > 0) {
		nx = selection->w;
		ny = selection->h;
		data = calloc(nx * ny, sizeof(WORD));
		select_area(fit, data, layer, selection);
		histo = computeHisto_Selection(fit, layer, selection);
	} else {
		nx = fit->rx;
		ny = fit->ry;
		data = calloc(nx * ny, sizeof(WORD));
		histo = computeHisto(fit, layer);
		memcpy(data, fit->pdata[layer], nx * ny * sizeof(WORD));
	}
	hist_size = gsl_histogram_bins(histo);

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

	/* Calculation of median with histogram */
	if ((option & STATS_BASIC) || (option & STATS_AVGDEV)
			|| (option & STATS_MAD) || (option & STATS_BWMV))
		median = siril_stats_ushort_median(histo, ngoodpix, nullcheck);
	gsl_histogram_free(histo);

	if (nullcheck) {
		data = reassign_data(data, nx * ny, ngoodpix);
	}

	if (option & STATS_BASIC)
		gsl_stats_ushort_minmax(&min, &max, data, 1, ngoodpix);

	/* Calculation of average absolute deviation from the median */
	if (option & STATS_AVGDEV)
		avgDev = gsl_stats_ushort_absdev_m(data, 1, ngoodpix, median);

	/* Calculation of median absolute deviation */
	if ((option & STATS_MAD) || (option & STATS_BWMV))
		mad = siril_stats_ushort_mad(data, 1, ngoodpix, median, nullcheck);

	/* Calculation of Bidweight Midvariance */
	if (option & STATS_BWMV)
		bwmv = siril_stats_ushort_bwmv(data, ngoodpix, mad, median);

	/* Calculation of IKSS. Used for stacking */
	if (option & (STATS_IKSS)) {
		double *newdata = calloc(ngoodpix, sizeof(double));

		/* we convert in the [0, 1] range */
		for (i = 0; i < ngoodpix; i++) {
			newdata[i] = (double) data[i] / ((double) hist_size - 1);
		}
		IKSS(newdata, ngoodpix, &location, &scale);
		/* go back to the original range */
		location *= ((double) hist_size - 1);
		scale *= ((double) hist_size - 1);
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
	stat->normValue = (double) hist_size - 1;

	free(data);
	return stat;
}
