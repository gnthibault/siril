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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "core/siril.h"
#include "io/image_format_fits.h"

#ifdef HAVE_WCSLIB
#include <wcslib.h>
#include <wcsfix.h>
#endif

#include "siril_wcs.h"

#ifdef HAVE_WCSLIB
static GMutex wcs_mutex;
#endif

/* we force naxis to 2 */
#define NAXIS 2

gboolean has_wcs(fits *fit) {
#ifdef HAVE_WCSLIB
	return fit->wcslib != NULL;
#endif
	return FALSE;
}

void free_wcs(fits *fit) {
#ifdef HAVE_WCSLIB
	if (fit->wcslib) {
		if (!wcsfree(fit->wcslib))
			free(fit->wcslib);
		fit->wcslib = NULL;
	}
#endif
}

gboolean load_WCS_from_memory(fits *fit) {
#ifdef HAVE_WCSLIB
	int status;
	if (!fit->wcslib) {
		fit->wcslib = calloc(1, sizeof(struct wcsprm));
		fit->wcslib->flag = -1;
	}
	wcsinit(1, NAXIS, fit->wcslib, 0, 0, 0);

	const char CTYPE[2][9] = {"RA---TAN", "DEC--TAN"};

	double *cdij = fit->wcslib->cd;
	double *pcij = fit->wcslib->pc;
	for (int i = 0; i < NAXIS; i++) {
		for (int j = 0; j < NAXIS; j++) {
			*(cdij++) = fit->wcsdata.cd[i][j];
			*(pcij++) = fit->wcsdata.cd[i][j];
		}
	}

	for (int i = 0; i < NAXIS; i++) {
		fit->wcslib->crval[i] = fit->wcsdata.crval[i];
	}

	for (int i = 0; i < NAXIS; i++) {
		fit->wcslib->crota[i] = fit->wcsdata.crota[i];
	}

	for (int i = 0; i < NAXIS; i++) {
		fit->wcslib->crpix[i] = fit->wcsdata.crpix[i];
	}

	for (int i = 0; i < NAXIS; i++) {
		fit->wcslib->cdelt[i] = 1;
	}

	for (int i = 0; i < NAXIS; i++) {
		strcpy(fit->wcslib->ctype[i], &CTYPE[i][0]);
	}

	fit->wcslib->equinox = fit->wcsdata.equinox;
//	fit->wcslib->lonpole = 180;
	fit->wcslib->latpole = fit->wcsdata.crval[1];

	if ((status = wcsset(fit->wcslib)) != 0) {
		free_wcs(fit);
		siril_debug_print("wcsset error %d: %s.\n", status, wcs_errmsg[status]);
		return FALSE;
	}
	return TRUE;
#else
	return FALSE;
#endif
}


gboolean load_WCS_from_file(fits* fit) {
#ifdef HAVE_WCSLIB
	int status = 0, wcs_status = 0;
	char *header;
	struct wcsprm *data = NULL;
	int nkeyrec, nreject, nwcs;

	/* sanity check to avoid error in some strange files */
	if ((fit->wcsdata.crpix[0] == 0) && (fit->wcsdata.crpix[1] == 0)) {
		return FALSE;
	}

	if (fit->wcslib) {
		free_wcs(fit);
	}

	ffhdr2str(fit->fptr, 1, NULL, 0, &header, &nkeyrec, &status);
	if (status) {
		report_fits_error(status);
		return FALSE;
	}

	/** There is a bug with wcspih that it is not really thread-safe for wcslib version < 7.5.
	 * We need to lock it.
	 * TODO: check wcslib version ?*/
	g_mutex_lock(&wcs_mutex);
	wcs_status = wcspih(header, nkeyrec, 0, 0, &nreject, &nwcs, &data);

	if (wcs_status == 0) {
		for (int i = 0; i < nwcs; i++) {
			/* Find the master celestial WCS coordinates */
			struct wcsprm *prm = data + i;
//			/* ctype3 = 'RGB' fix */
//			if (prm->naxis == 3) {
//				cdfix(prm);
//			}
			wcsset(prm);
			if (prm->lng >= 0 && prm->lat >= 0
					&& (prm->alt[0] == '\0' || prm->alt[0] == ' ')) {
				int axes[2], nsub;
				nsub = 2;
				axes[0] = WCSSUB_LONGITUDE;
				axes[1] = WCSSUB_LATITUDE;
				fit->wcslib = (struct wcsprm*) calloc(1, sizeof(struct wcsprm));
				fit->wcslib->flag = -1;
				status = wcssub(1, prm, &nsub, axes, fit->wcslib);
				if (status == 0) {
					break;
				} else {
					siril_debug_print("wcssub error %d: %s.\n", status, wcs_errmsg[status]);
					wcsvfree(&nwcs, &fit->wcslib);
					fit->wcslib = NULL;
				}
			}
		}
	}
	wcsvfree(&nwcs, &data);
	g_mutex_unlock(&wcs_mutex);
	free(header);

	if (!fit->wcslib) {
		siril_debug_print("No world coordinate systems found.\n");
		wcsvfree(&nwcs, &fit->wcslib);
		fit->wcslib = NULL;
		return FALSE;
	}

	return TRUE;
#else
	return FALSE;
#endif
}

void pix2wcs(fits *fit, double x, double y, double *r, double *d) {
	*r = -1.0;
	*d = -1.0;
#ifdef HAVE_WCSLIB
	int status, stat[NWCSFIX];
	double imgcrd[NWCSFIX], phi, pixcrd[NWCSFIX], theta, world[NWCSFIX];

	pixcrd[0] = x;
	pixcrd[1] = y;

	status = wcsp2s(fit->wcslib, 1, 2, pixcrd, imgcrd, &phi, &theta, world, stat);
	if (status != 0)
		return;

	*r = world[0];
	*d = world[1];
#endif
}

void wcs2pix(fits *fit, double r, double d, double *x, double *y) {
	*x = -1.0;
	*y = -1.0;
#ifdef HAVE_WCSLIB
	int status, stat[NWCSFIX];
	double imgcrd[NWCSFIX], phi, pixcrd[NWCSFIX], theta, world[NWCSFIX];

	world[0] = r;
	world[1] = d;

	status = wcss2p(fit->wcslib, 1, 2, world, &phi, &theta, imgcrd, pixcrd, stat);
	if (status != 0)
		return;

	*x = pixcrd[0];
	*y = pixcrd[1];
#endif
}

/* get resolution in arcsec/pixel */
double get_wcs_image_resolution(fits *fit) {
	double resolution = -1.0;
#ifdef HAVE_WCSLIB
	if (fit->wcslib) {
		double res_x = sqrt(fit->wcslib->cd[0] * fit->wcslib->cd[0] + fit->wcslib->cd[2] * fit->wcslib->cd[2]);
		double res_y = sqrt(fit->wcslib->cd[1] * fit->wcslib->cd[1] + fit->wcslib->cd[3] * fit->wcslib->cd[3]);
		resolution = (res_x + res_y) * 0.5;
	}
#endif
	if (resolution <= 0.0) {
		if (gfit.focal_length >= 0.0 && gfit.pixel_size_x >= 0.0 && gfit.pixel_size_y == gfit.pixel_size_x)
			resolution = (RADCONV / gfit.focal_length * gfit.pixel_size_x) / 3600;
	}
	return resolution;
}

double *get_wcs_crval(fits *fit) {
#ifdef HAVE_WCSLIB
	static double ret[NWCSFIX] = { 0 };
	if (fit->wcslib) {
		for (int i = 0; i < NAXIS; i++) {
			ret[i] = fit->wcslib->crval[i];
		}
	}
	return ret;
#else
	return NULL;
#endif
}

