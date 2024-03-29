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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_world_cs.h"
#include "algos/photometry.h"
#include "algos/sorting.h"
#include "algos/siril_wcs.h"
#include "filters/median.h"

#include "PSF.h"


#define MAX_ITER_NO_ANGLE  10		//Number of iteration in the minimization with no angle
#define MAX_ITER_ANGLE     10		//Number of iteration in the minimization with angle
#define EPSILON            0.001

const double radian_conversion = ((3600.0 * 180.0) / M_PI) / 1.0E3;

static gsl_matrix *removeHotPixels(gsl_matrix *in) {
	size_t width = in->size2;
	size_t height = in->size1;
	size_t x, y;
	gsl_matrix *out = gsl_matrix_alloc (in->size1, in->size2);
	if (!out) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	gsl_matrix_memcpy (out, in);
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			double a = get_median_gsl(in, x, y, width, height, 1, FALSE, FALSE);
			gsl_matrix_set(out, y, x, a);
		}
	}
	return out;
}

/* Compute initial values for the algorithm from data in the pixel value matrix */
static gsl_vector* psf_init_data(gsl_matrix* z, double bg) {
	gsl_vector * MaxV = gsl_vector_alloc(5);
	double max;
	size_t NbRows = z->size1;
	size_t NbCols = z->size2;
	size_t i, j;

	/* find maximum */
	/* first we remove hot pixels in the matrix */
	gsl_matrix *m_tmp = removeHotPixels(z);
	if (!m_tmp) return NULL;
	max = gsl_matrix_max(m_tmp);
	gsl_matrix_max_index(m_tmp, &i, &j);
	gsl_matrix_free(m_tmp);
	gsl_vector_set(MaxV, 0, j);
	gsl_vector_set(MaxV, 1, i);
	gsl_vector_set(MaxV, 2, max);

	size_t ii1 = (size_t) gsl_vector_get(MaxV, 1);
	size_t ii2 = (size_t) gsl_vector_get(MaxV, 1);
	size_t jj1 = (size_t) gsl_vector_get(MaxV, 0);
	size_t jj2 = (size_t) gsl_vector_get(MaxV, 0);
	size_t perm1 = (size_t) gsl_vector_get(MaxV, 1);
	size_t perm2 = (size_t) gsl_vector_get(MaxV, 0);

	while ((2.0 * (gsl_matrix_get(z, ii1, perm2) - bg)
			> (gsl_matrix_get(z, perm1, perm2) - bg)) && (ii1 < NbRows - 1.0)) {
		ii1++;
	}
	while ((2.0 * (gsl_matrix_get(z, ii2, perm2) - bg)
			> (gsl_matrix_get(z, perm1, perm2) - bg)) && (ii2 > 0)) {
		ii2--;
	}

	while ((2.0 * (gsl_matrix_get(z, perm1, jj1) - bg)
			> (gsl_matrix_get(z, perm1, perm2) - bg)) && (jj1 < NbCols - 1)) {
		jj1++;
	}
	while ((2.0 * (gsl_matrix_get(z, perm1, jj2) - bg)
			> (gsl_matrix_get(z, perm1, perm2) - bg)) && (jj2 > 0)) {
		jj2--;
	}
	gsl_vector_set(MaxV, 0, (jj1 + jj2 + 2) / 2.0);
	gsl_vector_set(MaxV, 1, (ii1 + ii2 + 2) / 2.0);
	gsl_vector_set(MaxV, 3, (size_t) (SQR(ii1 - ii2) / 4.0 / log(2.0)));
	gsl_vector_set(MaxV, 4, (size_t) (SQR(jj1 - jj2) / 4.0 / log(2.0)));

	return MaxV;
}

/* Basic magnitude computation. This is not really accurate, all pixels are
 * taken into account. But this fast function is used if the other one
 * failed and for star detection when magnitude is not needed.
 */
static double psf_get_mag(gsl_matrix* z, double B) {
	double intensity = 1.0;
	size_t NbRows = z->size1;
	size_t NbCols = z->size2;

	for (size_t i = 0; i < NbRows; i++) {
		for (size_t j = 0; j < NbCols; j++)
			intensity += gsl_matrix_get(z, i, j) - B;
	}
	return -2.5 * log10(intensity);
}

/* No angle */
static int psf_Gaussian_f(const gsl_vector * x, void *PSF_data, gsl_vector * f) {
	size_t NbRows = ((struct PSF_data *) PSF_data)->NbRows;
	size_t NbCols = ((struct PSF_data *) PSF_data)->NbCols;
	size_t i, j;
	double *y = ((struct PSF_data *) PSF_data)->y;
	double *sigma = ((struct PSF_data *) PSF_data)->sigma;
	double B = gsl_vector_get(x, 0);
	double A = gsl_vector_get(x, 1);
	double x0 = gsl_vector_get(x, 2);
	double y0 = gsl_vector_get(x, 3);
	double SX = gsl_vector_get(x, 4);
	double SY = gsl_vector_get(x, 5);
	double tmpx, tmpy, tmpc, sumres = 0.;

	for (i = 0; i < NbRows; i++) {
		for (j = 0; j < NbCols; j++) {
			tmpx = j + 1;
			tmpy = i + 1;
			tmpc = exp(-(SQR(tmpx-x0) / SX + SQR(tmpy-y0) / SY));
			gsl_vector_set(f, NbCols * i + j,
					(B + A * tmpc - y[NbCols * i + j]) / sigma[NbCols * i + j]);
			sumres += (B + A * tmpc - y[NbCols * i + j])
					* (B + A * tmpc - y[NbCols * i + j]);
		}
	}
	((struct PSF_data *) PSF_data)->rmse = sqrt(sumres / (NbRows * NbCols));
	return GSL_SUCCESS;
}

static int psf_Gaussian_df(const gsl_vector * x, void *PSF_data, gsl_matrix * J) {
	size_t NbRows = ((struct PSF_data *) PSF_data)->NbRows;
	size_t NbCols = ((struct PSF_data *) PSF_data)->NbCols;
	size_t i, j;
	double *sigma = ((struct PSF_data *) PSF_data)->sigma;
	double A = gsl_vector_get(x, 1);
	double x0 = gsl_vector_get(x, 2);
	double y0 = gsl_vector_get(x, 3);
	double SX = gsl_vector_get(x, 4);
	double SY = gsl_vector_get(x, 5);
	double tmpx, tmpy, tmpc, tmpd;

	for (i = 0; i < NbRows; i++) {
		for (j = 0; j < NbCols; j++) {
			tmpx = j + 1;
			tmpy = i + 1;
			double s = sigma[NbCols * i + j];
			tmpc = exp(-(SQR(tmpx-x0) / SX + SQR(tmpy-y0) / SY));
			gsl_matrix_set(J, NbCols * i + j, 0, 1. / s);
			gsl_matrix_set(J, NbCols * i + j, 1, tmpc / s);
			tmpd = A * tmpc * 2 * (tmpx - x0) / SX;
			gsl_matrix_set(J, NbCols * i + j, 2, tmpd / s);
			tmpd = A * tmpc * 2 * (tmpy - y0) / SY;
			gsl_matrix_set(J, NbCols * i + j, 3, tmpd / s);
			tmpd = A * tmpc * SQR(tmpx - x0) / SQR(SX);
			gsl_matrix_set(J, NbCols * i + j, 4, tmpd / s);
			tmpd = A * tmpc * SQR(tmpy - y0) / SQR(SY);
			gsl_matrix_set(J, NbCols * i + j, 5, tmpd / s);
		}
	}
	return GSL_SUCCESS;
}

int psf_Gaussian_fdf(const gsl_vector * x, void *PSF_data, gsl_vector * f,
		gsl_matrix * J) {
	psf_Gaussian_f(x, PSF_data, f);
	psf_Gaussian_df(x, PSF_data, J);
	return GSL_SUCCESS;
}

/* Angle */
static int psf_Gaussian_f_an(const gsl_vector * x, void *PSF_data,
		gsl_vector * f) {
	size_t NbRows = ((struct PSF_data *) PSF_data)->NbRows;
	size_t NbCols = ((struct PSF_data *) PSF_data)->NbCols;
	size_t n = ((struct PSF_data *) PSF_data)->n;
	size_t i, j;
	double *y = ((struct PSF_data *) PSF_data)->y;
	double *sigma = ((struct PSF_data *) PSF_data)->sigma;
	double B = gsl_vector_get(x, 0);
	double A = gsl_vector_get(x, 1);
	double x0 = gsl_vector_get(x, 2);
	double y0 = gsl_vector_get(x, 3);
	double SX = gsl_vector_get(x, 4);
	double SY = gsl_vector_get(x, 5);
	double alpha = gsl_vector_get(x, 6);
	double tmpx, tmpy, tmpc, sumres = 0.;

	for (i = 0; i < NbRows; i++) {
		for (j = 0; j < NbCols; j++) {
			tmpx = cos(alpha) * (j + 1 - x0) - sin(alpha) * (i + 1 - y0) + x0;
			tmpy = sin(alpha) * (j + 1 - x0) + cos(alpha) * (i + 1 - y0) + y0;
			tmpc = exp(-(SQR(tmpx-x0) / SX + SQR(tmpy-y0) / SY));
			gsl_vector_set(f, NbCols * i + j,
					(B + A * tmpc - y[NbCols * i + j]) / sigma[NbCols * i + j]);
			sumres += (B + A * tmpc - y[NbCols * i + j])
					* (B + A * tmpc - y[NbCols * i + j]);
		}
	}
	((struct PSF_data *) PSF_data)->rmse = sqrt(sumres / n);
	return GSL_SUCCESS;
}

static int psf_Gaussian_df_an(const gsl_vector * x, void *PSF_data,
		gsl_matrix * J) {
	size_t NbRows = ((struct PSF_data *) PSF_data)->NbRows;
	size_t NbCols = ((struct PSF_data *) PSF_data)->NbCols;
	size_t i, j;
	double *sigma = ((struct PSF_data *) PSF_data)->sigma;
	double A = gsl_vector_get(x, 1);
	double x0 = gsl_vector_get(x, 2);
	double y0 = gsl_vector_get(x, 3);
	double SX = gsl_vector_get(x, 4);
	double SY = gsl_vector_get(x, 5);
	double alpha = gsl_vector_get(x, 6);
	double tmpx, tmpy, tmpc, tmpd, tmpderxr, tmpderyr;
	double s;

	for (i = 0; i < NbRows; i++) {
		for (j = 0; j < NbCols; j++) {
			tmpx = cos(alpha) * (j + 1 - x0) - sin(alpha) * (i + 1 - y0) + x0;
			tmpy = sin(alpha) * (j + 1 - x0) + cos(alpha) * (i + 1 - y0) + y0;
			s = sigma[NbCols * i + j];
			tmpc = exp(-(SQR(tmpx-x0) / SX + SQR(tmpy-y0) / SY));
			gsl_matrix_set(J, NbCols * i + j, 0, 1. / s);
			gsl_matrix_set(J, NbCols * i + j, 1, tmpc / s);
			tmpd = A * tmpc * 2 * (tmpx - x0) / SX * cos(alpha);
			gsl_matrix_set(J, NbCols * i + j, 2, tmpd / s);
			tmpd = A * tmpc * 2 * (tmpy - y0) / SY * cos(alpha);
			gsl_matrix_set(J, NbCols * i + j, 3, tmpd / s);
			tmpd = A * tmpc * SQR(tmpx - x0) / SQR(SX);
			gsl_matrix_set(J, NbCols * i + j, 4, tmpd / s);
			tmpd = A * tmpc * SQR(tmpy - y0) / SQR(SY);
			gsl_matrix_set(J, NbCols * i + j, 5, tmpd / s);
			tmpderxr = -sin(alpha) * (j + 1 - x0) - cos(alpha) * (i + 1 - y0);
			tmpderyr = cos(alpha) * (j + 1 - x0) - sin(alpha) * (i + 1 - y0);
			tmpd = -A * tmpc
					* (2 * (tmpx - x0) / SX * tmpderxr
							+ 2 * (tmpy - y0) / SY * tmpderyr);
			gsl_matrix_set(J, NbCols * i + j, 6, tmpd / s);
		}
	}
	return GSL_SUCCESS;
}

static int psf_Gaussian_fdf_an(const gsl_vector * x, void *PSF_data,
		gsl_vector * f, gsl_matrix * J) {
	psf_Gaussian_f_an(x, PSF_data, f);
	psf_Gaussian_df_an(x, PSF_data, J);
	return GSL_SUCCESS;
}

/* The function returns the fitted parameters without angle. However it
 * returns NULL if the number of parameters is => to the pixel number.
 */
static psf_star *psf_minimiz_no_angle(gsl_matrix* z, double background) {
	size_t i, j;
	size_t NbRows = z->size1; //characteristics of the selection : height and width
	size_t NbCols = z->size2;
	const size_t p = 6;			// Number of parameters fitted
	const size_t n = NbRows * NbCols;
	gsl_vector *MaxV = psf_init_data(z, background);
	if (!MaxV)
		return NULL;
	int status;
	unsigned int iter = 0;
	gsl_matrix *covar = gsl_matrix_alloc(p, p);
	double *y = malloc(n * sizeof(double));
	double *sigma = malloc(n * sizeof(double));
	psf_star *psf = new_psf_star();
	if (!y || !sigma || !psf) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	struct PSF_data d = { n, y, sigma, NbRows, NbCols, 0 };
	gsl_multifit_function_fdf f;
	const gsl_rng_type * type;
	gsl_rng * r;
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	double x_init[] = { background, gsl_vector_get(MaxV, 2), gsl_vector_get(
			MaxV, 0), gsl_vector_get(MaxV, 1), gsl_vector_get(MaxV, 4),
			gsl_vector_get(MaxV, 3) };
	gsl_vector_view x = gsl_vector_view_array(x_init, p);

	gsl_rng_env_setup();

	type = gsl_rng_default;
	r = gsl_rng_alloc(type);

	f.f = &psf_Gaussian_f;
	f.df = &psf_Gaussian_df;
	f.fdf = &psf_Gaussian_fdf;
	f.n = n;
	if (n <= p) {
		free_psf(psf);
		free(y);
		free(sigma);
		return NULL;
	}
	f.p = p;
	f.params = &d;

	for (i = 0; i < NbRows; i++) {
		for (j = 0; j < NbCols; j++) {
			y[NbCols * i + j] = gsl_matrix_get(z, i, j);
			sigma[NbCols * i + j] = 1.0;
		}
	}

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc(T, n, p);
	gsl_multifit_fdfsolver_set(s, &f, &x.vector);

	do {
		iter++;
		status = gsl_multifit_fdfsolver_iterate(s);
		if (status)
			break;
		status = gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4);
	} while (status == GSL_CONTINUE && iter < MAX_ITER_NO_ANGLE);

#if HAVE_GSL_1
	gsl_multifit_covar(s->J, 0.0, covar);
#elif HAVE_GSL_2
	gsl_matrix * J = gsl_matrix_alloc(n, p);

	gsl_multifit_fdfsolver_jac(s, J);
	gsl_multifit_covar(J, 0.0, covar);

	gsl_matrix_free(J);
#endif

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))	//for now, errors are not displayed

	/* Output structure with parameters fitted */
	psf->B = FIT(0);
	psf->A = FIT(1);
	psf->x0 = FIT(2);
	psf->y0 = FIT(3);
	psf->sx = FIT(4);
	psf->sy = FIT(5);
	psf->fwhmx = sqrt(FIT(4) / 2.) * 2 * sqrt(log(2.) * 2);	//Set the real FWHMx with regards to the Sx parameter
	psf->fwhmy = sqrt(FIT(5) / 2.) * 2 * sqrt(log(2.) * 2);	//Set the real FWHMy with regards to the Sy parameter
	psf->fwhmx_arcsec = -1.0;
	psf->fwhmy_arcsec = -1.0;
	psf->angle = 0;	//The angle is not fitted here
	// Units
	psf->units = "px";
	// Magnitude
	psf->mag = psf_get_mag(z, psf->B);
	psf->phot = NULL;
	psf->phot_is_valid = FALSE;
	// RMSE
	psf->rmse = d.rmse;
	// absolute uncertainties
	psf->B_err = ERR(0) / FIT(0);
	psf->A_err = ERR(1) / FIT(1);
	psf->x_err = ERR(2) / FIT(2);
	psf->y_err = ERR(3) / FIT(3);
	psf->sx_err = ERR(4) / FIT(4);
	psf->sy_err = ERR(5) / FIT(5);
	psf->ang_err = 0;
	psf->xpos = 0;		// will be set by the peaker
	psf->ypos = 0;
	// we free the memory
	free(sigma);
	free(y);
	gsl_vector_free(MaxV);
	gsl_multifit_fdfsolver_free(s);
	gsl_matrix_free(covar);
	gsl_rng_free(r);
	return psf;
}

/* The function returns the fitted parameters with angle. However it returns
 * NULL if the number of parameters is => to the pixel number.
 * This should not happen because this case is already treated by the
 * minimiz_no_angle function */
static psf_star *psf_minimiz_angle(gsl_matrix* z, psf_star *psf, gboolean for_photometry, gboolean verbose) {
	size_t i, j;
	size_t NbRows = z->size1; //characteristics of the selection : height and width
	size_t NbCols = z->size2;
	const size_t p = 7;			// Number of parameters fitted
	const size_t n = NbRows * NbCols;
	g_assert (n > 0);
	int status;
	unsigned int iter = 0;
	psf_star *psf_angle = new_psf_star();
	gsl_matrix *covar = gsl_matrix_alloc(p, p);
	double *y = malloc(n * sizeof(double));
	double *sigma = malloc(n * sizeof(double));
	if (!psf_angle || !covar || !y || !sigma) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	struct PSF_data d = { n, y, sigma, NbRows, NbCols, 0 };
	gsl_multifit_function_fdf f_angle;
	double x_init[7] = { psf->B, psf->A, psf->x0, psf->y0, psf->sx, psf->sy, 0 };
	gsl_vector_view x = gsl_vector_view_array(x_init, p);
	const gsl_rng_type * type;
	gsl_rng * r;
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;

	gsl_rng_env_setup();

	type = gsl_rng_default;
	r = gsl_rng_alloc(type);

	f_angle.f = &psf_Gaussian_f_an;
	f_angle.df = &psf_Gaussian_df_an;
	f_angle.fdf = &psf_Gaussian_fdf_an;
	f_angle.n = n;
	f_angle.p = p;
	f_angle.params = &d;

	for (i = 0; i < NbRows; i++) {
		for (j = 0; j < NbCols; j++) {
			y[NbCols * i + j] = gsl_matrix_get(z, i, j);
			sigma[NbCols * i + j] = 1.0;
		}
	}

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc(T, n, p);
	gsl_multifit_fdfsolver_set(s, &f_angle, &x.vector);

	do {
		iter++;
		status = gsl_multifit_fdfsolver_iterate(s);
		if (status)
			break;
		status = gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4);
	} while (status == GSL_CONTINUE && iter < MAX_ITER_ANGLE);

#if HAVE_GSL_1
	gsl_multifit_covar(s->J, 0.0, covar);
#elif HAVE_GSL_2
	gsl_matrix * J;

	J = gsl_matrix_alloc(n, p);

	gsl_multifit_fdfsolver_jac(s, J);
	gsl_multifit_covar(J, 0.0, covar);

	gsl_matrix_free(J);
#endif

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))	//for now, errors are not displayed

	/*Output structure with parameters fitted */
	psf_angle->B = FIT(0);
	psf_angle->A = FIT(1);
	psf_angle->x0 = FIT(2);
	psf_angle->y0 = FIT(3);
	psf_angle->sx = FIT(4);
	psf_angle->sy = FIT(5);
	psf_angle->fwhmx = sqrt(FIT(4) / 2.0) * 2 * sqrt(log(2.0) * 2);	//Set the real FWHMx with regards to the Sx parameter
	psf_angle->fwhmy = sqrt(FIT(5) / 2.0) * 2 * sqrt(log(2.0) * 2);	//Set the real FWHMy with regards to the Sy parameter
	psf_angle->angle = -FIT(6) * 180.0 / M_PI;
	/* The angle must be => -90 and <= 90
	 * Otherwise, the solution may be degenerate
	 * and produce an angle > 90. So we're
	 * looking for the solution between
	 * the interval we want */
	while (fabs(psf_angle->angle) > 90.0) {
		if (psf_angle->angle > 0.0)
			psf_angle->angle -= 90.0;
		else
			psf_angle->angle += 90.0;
	}
	//Units
	psf_angle->units = "px";
	// Photometry
	if (for_photometry)
		psf_angle->phot = getPhotometryData(z, psf_angle, verbose);
	else {
		psf_angle->phot = NULL;
		psf_angle->phot_is_valid = FALSE;
	}
	// Magnitude
	if (psf_angle->phot != NULL) {
		psf_angle->mag = psf_angle->phot->mag;
		psf_angle->s_mag = psf_angle->phot->s_mag;
		psf_angle->SNR = psf_angle->phot->SNR;
		psf_angle->phot_is_valid = psf_angle->phot->valid;

	} else {
		psf_angle->mag = psf_get_mag(z, psf_angle->B);
		psf_angle->s_mag = 9.999;
		psf_angle->SNR = 0;
		psf_angle->phot_is_valid = FALSE;
	}
	//RMSE
	psf_angle->rmse = d.rmse;
	// absolute uncertainties
	psf_angle->B_err = ERR(0) / FIT(0);
	psf_angle->A_err = ERR(1) / FIT(1);
	psf_angle->x_err = ERR(2) / FIT(2);
	psf_angle->y_err = ERR(3) / FIT(3);
	psf_angle->sx_err = ERR(4) / FIT(4);
	psf_angle->sy_err = ERR(5) / FIT(5);
	psf_angle->ang_err = ERR(6) / FIT(6);

	//we free the memory
	free(sigma);
	free(y);
	gsl_multifit_fdfsolver_free(s);
	gsl_matrix_free(covar);
	gsl_rng_free(r);
	return psf_angle;
}

/******************************************************************************/

/* Returns the largest FWHM in pixels
 * The optional output parameter roundness is the ratio between the two axis FWHM */
double psf_get_fwhm(fits *fit, int layer, rectangle *selection, double *roundness) {
	psf_star *result = psf_get_minimisation(fit, layer, selection, FALSE, TRUE, TRUE);
	if (result == NULL) {
		*roundness = 0.0;
		return 0.0;
	}
	double retval;
	retval = result->fwhmx;
	if (roundness)
		*roundness = result->fwhmy / result->fwhmx;
	free_psf(result);
	return retval;
}

/* Computes the FWHM on data in the selection rectangle of image fit.
 * Selection rectangle is passed as third argument.
 * Return value is a structure, type psf_star, that has to be freed after use.
 */
psf_star *psf_get_minimisation(fits *fit, int layer, rectangle *area,
		gboolean for_photometry, gboolean verbose, gboolean multithread_stat) {
	int stridefrom, i, j;
	psf_star *result;
	double bg = background(fit, layer, area, multithread_stat);

	// fprintf(stdout, "background: %g\n", bg);
	gsl_matrix *z = gsl_matrix_alloc(area->h, area->w);
	stridefrom = fit->rx - area->w;

	// create the matrix with values from the selected rectangle
	if (fit->type == DATA_USHORT) {
		WORD *from = fit->pdata[layer] +
			(fit->ry - area->y - area->h) * fit->rx + area->x;

		for (i = 0; i < area->h; i++) {
			for (j = 0; j < area->w; j++) {
				gsl_matrix_set(z, i, j, (double)*from);
				from++;
			}
			from += stridefrom;
		}
	}
	else if (fit->type == DATA_FLOAT) {
		float *from = fit->fpdata[layer] +
			(fit->ry - area->y - area->h) * fit->rx + area->x;

		for (i = 0; i < area->h; i++) {
			for (j = 0; j < area->w; j++) {
				gsl_matrix_set(z, i, j, (double)*from);
				from++;
			}
			from += stridefrom;
		}
	}
	else {
		gsl_matrix_free(z);
		return NULL;
	}

	result = psf_global_minimisation(z, bg, TRUE, for_photometry, verbose);
	if (result) {
		fwhm_to_arcsec_if_needed(fit, result);
		result->layer = layer;
	}
	gsl_matrix_free(z);
	return result;
}

/* This function is the global minimisation. Every call to the minimisation
 * must come over here. It will check if the difference between Sx and Sy is
 * larger than or equal to 0.01 pixel.
 * In this case, Dynamic PSF fits additional angle parameter wich is the
 * rotation angle of the X axis with respect to the centroid coordinates: so,
 * by design we set Sx>Sy.
 * If the difference is smaller OR if fit_Angle is equal to FALSE (in the case
 * of the star_finder algorithm), no angle parameter is fitted.
 * The function returns NULL if values look bizarre.
 */
psf_star *psf_global_minimisation(gsl_matrix* z, double bg,
		gboolean fit_angle, gboolean for_photometry, gboolean verbose) {
	psf_star *psf;

	// To compute good starting values, we first compute with no angle
	if ((psf = psf_minimiz_no_angle(z, bg)) != NULL) {
		if (fit_angle) {
			/* This next check is to avoid possible angle divergence
			 * when sx and sy are too close (star is quite round).
			 * However, in this case we need to compute photometry if needed
			 */
			if ((fabs(psf->sx - psf->sy) < EPSILON)) {
				// Photometry
				if (for_photometry) {
					psf->phot = getPhotometryData(z, psf, verbose);
					if (psf->phot != NULL) {
						psf->mag = psf->phot->mag;
						psf->s_mag = psf->phot->s_mag;
						psf->SNR = psf->phot->SNR;
						psf->phot_is_valid = psf->phot->valid;
					}
				} else {
					psf->phot = NULL;
					psf->phot_is_valid = FALSE;
				}
			} else {
				psf_star *tmp_psf;
				if ((tmp_psf = psf_minimiz_angle(z, psf, for_photometry, verbose))
						== NULL) {
					free_psf(psf);
					return NULL;
				}
				free_psf(psf);
				psf = tmp_psf;
			}
		}
		// Solve symmetry problem in order to have Sx>Sy in any case !!!
		if (psf->sy > psf->sx) {
			SWAP(psf->sx, psf->sy);
			SWAP(psf->fwhmx, psf->fwhmy);
			if (fit_angle && psf->angle != 0.0) {
				if (psf->angle > 0.0)
					psf->angle -= 90.0;
				else psf->angle += 90.0;
			}
		}

		/* We quickly test the result. If it is bad we return NULL */
		if (!isfinite(psf->fwhmx) || !isfinite(psf->fwhmy) ||
				psf->fwhmx <= 0.0 || psf->fwhmy <= 0.0) {
			free_psf(psf);
			psf = NULL;
		}
	}
	/* When the first minimization gives NULL value, it's probably because the selected
	 * area was not big enough: we need more samples than parameters to fit the area */
	return psf;
}

void psf_display_result(psf_star *result, rectangle *area) {
	char *buffer, *coordinates;
	char *str;
	if (com.magOffset > 0.0)
		str = "true reduced";
	else
		str = "relative";

	double x = result->x0 + area->x;
	double y = area->y + area->h - result->y0;

	if (has_wcs(&gfit)) {
		double world_x, world_y;
		SirilWorldCS *world_cs;
		pix2wcs(&gfit, x, (double) gfit.ry - y, &world_x, &world_y);
		world_cs = siril_world_cs_new_from_a_d(world_x, world_y);
		if (world_cs) {
			gchar *ra = siril_world_cs_alpha_format(world_cs, "%02dh%02dm%02ds");
			gchar *dec = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%02d\"");
			coordinates = g_strdup_printf("x0=%0.2f px, y0=%0.2f px (%s , %s)", x, y, ra, dec);

			siril_world_cs_unref(world_cs);
			g_free(ra);
			g_free(dec);
		} else {
			coordinates = g_strdup_printf("x0=%0.2f px, y0=%0.2f px", x, y);
		}
	} else {
		coordinates = g_strdup_printf("x0=%0.2f px, y0=%0.2f px", x, y);
	}

	buffer = g_strdup_printf(_("PSF fit Result:\n"
			"%s\n"
			"FWHM X=%0.2f%s, FWHM Y=%0.2f%s\n"
			"Angle=%0.2f deg\n"
			"Background value=%0.6f\n"
			"Maximal intensity=%0.6f\n"
			"Magnitude (%s)=%0.2f\n"
			"SNR=%.1fdB\n"
			"RMSE=%.3e\n"),
			coordinates,
			result->fwhmx, result->units, result->fwhmy, result->units,
			result->angle,
			result->B,
			result->A,
			str,
			result->mag + com.magOffset,
			result->SNR,
			result->rmse);

	siril_log_message(buffer);
	g_free(buffer);
	g_free(coordinates);
}

#define _2_SQRT_2_LOG2 2.35482004503

/* If the pixel pitch and the focal length are known and filled in the 
 * setting box, we convert FWHM in pixel to arcsec by multiplying
 * the FWHM value with the sampling value */
void fwhm_to_arcsec_if_needed(fits* fit, psf_star *result) {

	if (!result) return;
	if (fit->focal_length <= 0.0 || fit->pixel_size_x <= 0.f
			|| fit->pixel_size_y <= 0.f || fit->binning_x <= 0
			|| fit->binning_y <= 0) {
		result->fwhmx_arcsec = -1.0;
		result->fwhmy_arcsec = -1.0;
		return;
	}

	double bin_X, bin_Y;
	double fwhmx, fwhmy;

	fwhmx = sqrt(result->sx * 0.5) * _2_SQRT_2_LOG2;
	fwhmy = sqrt(result->sy * 0.5) * _2_SQRT_2_LOG2;

	bin_X = fit->unbinned ? (double) fit->binning_x : 1.0;
	bin_Y = fit->unbinned ? (double) fit->binning_y : 1.0;

	result->fwhmx_arcsec = fwhmx * (radian_conversion * (double)fit->pixel_size_x / fit->focal_length) * bin_X;
	result->fwhmy_arcsec = fwhmy * (radian_conversion * (double)fit->pixel_size_y / fit->focal_length) * bin_Y;
	result->units = "\"";
}

void fwhm_to_pixels(psf_star *result) {
	result->fwhmx = sqrt(result->sx * 0.5) * _2_SQRT_2_LOG2;
	result->fwhmy = sqrt(result->sy * 0.5) * _2_SQRT_2_LOG2;
	result->units = "px";
}

// returns boolean if it was possible (true if arcsec)
gboolean get_fwhm_as_arcsec_if_possible(psf_star *star, double *fwhmx, double *fwhmy, char **unit) {
	if (!strcmp(star->units, "px")) {
		*fwhmx = star->fwhmx;
		*fwhmy = star->fwhmy;
		*unit = star->units;
		return FALSE;
	}
	if (star->fwhmx_arcsec <= 0.0) {
		fprintf(stderr, "FWHM wrongly stored as arcsec\n");
		star->units = "px";
		return get_fwhm_as_arcsec_if_possible(star, fwhmx, fwhmy, unit);
	}
	*fwhmx = star->fwhmx_arcsec;
	*fwhmy = star->fwhmy_arcsec;
	*unit = star->units;
	return TRUE;
}

double convert_single_fwhm_to_pixels(double fwhm, double s) {
	return sqrt(s * 0.5) * _2_SQRT_2_LOG2;
}

gboolean convert_single_fwhm_to_arcsec_if_possible(double fwhm, double bin, double px_size, double flenght, double *result) {
	double arcsec = fwhm * (radian_conversion * px_size / flenght) * bin;
	if (arcsec <= 0.0 || isnan(arcsec) || !isfinite(arcsec)) {
		*result = 0;
		return FALSE;
	}
	*result = arcsec;
	return TRUE;
}

psf_star *new_psf_star() {
	psf_star *star = malloc(sizeof(psf_star));
	star->phot = NULL;

	return star;
}

psf_star *duplicate_psf(psf_star *psf) {
	if (!psf)
		return NULL;
	psf_star *new_psf = new_psf_star();
	memcpy(new_psf, psf, sizeof(psf_star));
	if (psf->phot) {
		new_psf->phot = malloc(sizeof(photometry));
		memcpy(new_psf->phot, psf->phot, sizeof(photometry));
	} else {
		new_psf->phot = NULL;
	}
	return new_psf;
}

void free_psf(psf_star *psf) {
	if (psf->phot) free(psf->phot);
	free(psf);
}
