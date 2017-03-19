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

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "gui/callbacks.h"
#include "algos/PSF.h"
#include "algos/photometry.h"

#define hampel_a   1.7
#define hampel_b   3.4
#define hampel_c   8.5
#define sign(x,y)  ((y)>=0?fabs(x):-fabs(x))
#define epsilon(x) 0.00000001
#define maxit      50
#define min_sky    5
#define lo_data    0.0
#define hi_data    USHRT_MAX_DOUBLE

static void initializeParam() {
	com.phot_set.inner = 20;
	com.phot_set.outer = 30;
	com.phot_set.gain = 2.3;
}

static double hampel(double x) {
	if (x >= 0) {
		if (x < hampel_a)
			return x;
		if (x < hampel_b)
			return hampel_a;
		if (x < hampel_c)
			return hampel_a * (x - hampel_c) / (hampel_b - hampel_c);
	} else {
		if (x > -hampel_a)
			return x;
		if (x > -hampel_b)
			return -hampel_a;
		if (x > -hampel_c)
			return hampel_a * (x + hampel_c) / (hampel_b - hampel_c);
	}
	return 0.0;
}

static double dhampel(double x) {
	if (x >= 0) {
		if (x < hampel_a)
			return 1;
		if (x < hampel_b)
			return 0;
		if (x < hampel_c)
			return hampel_a / (hampel_b - hampel_c);
	} else {
		if (x > -hampel_a)
			return 1;
		if (x > -hampel_b)
			return 0;
		if (x > -hampel_c)
			return -hampel_a / (hampel_b - hampel_c);
	}
	return 0.0;
}

static double qmedD(int n, double *a)
/* Vypocet medianu algoritmem Quick Median (Wirth) */
{
	double w, x;
	int i, j;
	int k = ((n & 1) ? (n / 2) : ((n / 2) - 1));
	int l = 0;
	int r = n - 1;

	while (l < r) {
		x = a[k];
		i = l;
		j = r;
		do {
			while (a[i] < x)
				i++;
			while (x < a[j])
				j--;
			if (i <= j) {
				w = a[i];
				a[i] = a[j];
				a[j] = w;
				i++;
				j--;
			}
		} while (i <= j);
		if (j < k)
			l = i;
		if (k < i)
			r = j;
	}
	return a[k];
}

static int robustmean(int n, double *x, double *mean, double *stdev)
/* Newton's iterations */
{
	int i, it;
	double a, c, d, dt, r, s, sum1, sum2, sum3, psir;
	double *xx;

	if (n < 1) {
		if (mean)
			*mean = 0.0; /* a few data */
		if (stdev)
			*stdev = -1.0;
		return 1;
	}
	if (n == 1) { /* only one point, but correct case */
		if (mean)
			*mean = x[0];
		if (stdev)
			*stdev = 0.0;
		return 0;
	}

	/* initial values:
	 - median is the first approximation of location
	 - MAD/0.6745 is the first approximation of scale */
	xx = malloc(n * sizeof(double));
	memcpy(xx, x, n * sizeof(double));
	a = qmedD(n, xx);
	for (i = 0; i < n; i++)
		xx[i] = fabs(x[i] - a);
	s = qmedD(n, xx) / 0.6745;
	free(xx);

	/* almost identical points on input */
	if (fabs(s) < epsilon(s)) {
		if (mean)
			*mean = a;
		if (stdev) {
			double sum = 0.0;
			for (i = 0; i < n; i++)
				sum += (x[i] - a) * (x[i] - a);
			*stdev = sqrt(sum / n);
		}
		return 0;
	}

	/* corrector's estimation */
	dt = 0;
	c = s * s * n * n / (n - 1);
	for (it = 1; it <= maxit; it++) {
		sum1 = sum2 = sum3 = 0.0;
		for (i = 0; i < n; i++) {
			r = (x[i] - a) / s;
			psir = hampel(r);
			sum1 += psir;
			sum2 += dhampel(r);
			sum3 += psir * psir;
		}
		if (fabs(sum2) < epsilon(sum2))
			break;
		d = s * sum1 / sum2;
		a = a + d;
		dt = c * sum3 / (sum2 * sum2);
		if ((it > 2) && ((d * d < 1e-4 * dt) || (fabs(d) < 10.0 * epsilon(d))))
			break;
	}
	if (mean)
		*mean = a;
	if (stdev)
		*stdev = (dt > 0 ? sqrt(dt) : 0);
	return 0;
}

static photometry *phot_alloc(size_t size) {
	photometry *phot;

	phot = (photometry*) calloc((unsigned) size * sizeof(photometry), 1);
	if (phot == NULL) {
		printf("photometry: memory error\n");
	}
	return (phot);
}

static double getMagnitude(double intensity) {
	return -2.5 * log10(intensity);
}

static double getCameraGain() {
	return com.phot_set.gain;
}

static double getInnerRadius() {
	return com.phot_set.inner;
}

static double getOuterRadius() {
	return com.phot_set.outer;
}

static double getMagErr(double intensity, double area, int nsky, double skysig) {
	double skyvar, sigsq;
	double err1, err2, err3;
	double phpadu;

	skyvar = skysig * skysig; /* variance of the sky brightness */
	sigsq = skyvar / nsky; /* square of the standard error of the mean sky brightness */
	phpadu = getCameraGain();
	err1 = area * skyvar;
	err2 = intensity / phpadu;
	err3 = sigsq * area * area;

	return fmin(9.999, 1.0857 * sqrt(err1 + err2 + err3) / intensity);
}

/* Function that compute all photometric data. The result must be freed */
photometry *getPhotometryData(gsl_matrix* z, fitted_PSF *psf) {
	int width = z->size2;
	int height = z->size1;
	int n_sky = 0, ret;
	int x, y, x1, y1, x2, y2;
	double r1, r2, r, rmin_sq, appRadius;
	double xc, yc;
	double apmag = 0.0, mean = 0.0, stdev = 0.0, area = 0.0;
	double signalIntensity;
	double *data;
	photometry *phot;

	xc = psf->x0 - 1;
	yc = psf->y0 - 1;

	r1 = getInnerRadius();
	r2 = getOuterRadius();
	appRadius = sqrt(psf->sx / 2.0) * 2 * sqrt(log(2.0) * 2) + 0.5;
	if (appRadius >= r1) {
		/* Translator note: radii is plural for radius */
		siril_log_message(_("Inner and outer radii are too small. Please update values in setting box.\n"));
		return NULL;
	}

	x1 = xc - r2;
	if (x1 < 1)
		x1 = 1;
	x2 = xc + r2;
	if (x2 > width - 1)
		x2 = width - 1;
	y1 = yc - r2;
	if (y1 < 1)
		y1 = 1;
	y2 = yc + r2;
	if (y2 > height - 1)
		y2 = height - 1;

	r1 *= r1;
	r2 *= r2;
	rmin_sq = (appRadius - 0.5) * (appRadius - 0.5);

	data = calloc((y2 - y1) * (x2 - x1), sizeof(double));

	for (y = y1; y <= y2; ++y) {
		int yp = (y - yc) * (y - yc);
		for (x = x1; x <= x2; ++x) {
			r = yp + (x - xc) * (x - xc);
			double pixel = gsl_matrix_get(z, y, x);
			if (pixel > lo_data && pixel < hi_data) {
				double f = (r < rmin_sq ? 1 : appRadius - sqrt(r) + 0.5);
				if (f >= 0) {
					area += f;
					apmag += pixel * f;
				}
				/* annulus */
				if (r < r2 && r > r1) {
					data[n_sky] = pixel;
					n_sky++;
				}
			}
		}
	}
	if (area < 1) {
		free(data);
		return NULL;
	}

	if (n_sky < min_sky) {
		siril_log_message(_("Warning: There aren't enough pixels"
				" in the sky annulus. You need to make a larger selection.\n"));
		free(data);
		return NULL;
	}

	ret = robustmean(n_sky, data, &mean, &stdev);
	if (ret > 0) {
		free(data);
		return NULL;
	}

	phot = phot_alloc(1);
	if (phot) {
		signalIntensity = apmag - (area * mean);
		phot->mag = getMagnitude(signalIntensity);
		phot->s_mag = getMagErr(signalIntensity, area, n_sky, stdev);
	}

	free(data);
	return phot;
}

void on_button_reset_photometry_clicked(GtkButton *button, gpointer user_data) {
	double tmp;

	initializeParam();
	tmp = gfit.cvf;
	gfit.cvf = 0.0;
	set_GUI_photometry();
	gfit.cvf = tmp;
}

