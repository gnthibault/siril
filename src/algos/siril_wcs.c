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

/* we force naxis to 2 */
#define NAXIS 2

static struct wcsprm *wcs = NULL;

gboolean has_wcs(){
	return wcs != NULL;
}

//static int parse_header_str(gchar *orig_header, gchar **header, int *nkeys, int *status) {
//	char keybuf[162], keyname[FLEN_KEYWORD], *headptr;
//	gchar **token = g_strsplit_set(orig_header, "\n", -1);
//	int totkeys = g_strv_length(token);
//
//	*nkeys = 0;
//
//	*header = (char*) calloc((totkeys + 1) * 80 + 1, 1);
//	if (!(*header)) {
//		PRINT_ALLOC_ERR;
//	}
//
//    headptr = *header;
//
//	/* read every keyword */
//	for (int i = 0; i < totkeys; i++) {
//		strcpy(keybuf, token[i]);
//		/* pad record with blanks so that it is at least 80 chars long */
//		strcat(keybuf,
//				"                                                                                ");
//
//		keyname[0] = '\0';
//		strncat(keyname, keybuf, 8); /* copy the keyword name */
//
//		if (!g_strcmp0("COMMENT ", keyname) || !g_strcmp0("HISTORY ", keyname)
//				|| !g_strcmp0("        ", keyname))
//			continue; /* skip this commentary keyword */
//
//		/* not in exclusion list, add this keyword to the string */
//		strcpy(headptr, keybuf);
//		headptr += 80;
//		(*nkeys)++;
//	}
//
//
//    /* add the END keyword */
//    strcpy(headptr,
//    "END                                                                             ");
//    headptr += 80;
//    (*nkeys)++;
//    g_strfreev(token);
//
//    *headptr = '\0';   /* terminate the header string */
//    /* minimize the allocated memory */
//    *header = (char *) realloc(*header, (*nkeys *80) + 1);
//    return *status;
//}

gboolean load_WCS_from_memory(fits *fit) {
#ifdef HAVE_WCSLIB
	int status;
	if (!has_wcs()) {
		wcs = calloc(1, sizeof(struct wcsprm));
		wcs->flag = -1;
	}
	wcsinit(1, NAXIS, wcs, 0, 0, 0);

	char CTYPE[2][9] = {"RA---TAN", "DEC--TAN"};

	double *cdij = wcs->cd;
	double *pcij = wcs->pc;
	for (int i = 0; i < NAXIS; i++) {
		for (int j = 0; j < NAXIS; j++) {
			*(cdij++) = fit->wcs.cd[i][j];
			*(pcij++) = fit->wcs.cd[i][j];
		}
	}

	for (int i = 0; i < NAXIS; i++) {
		wcs->crval[i] = fit->wcs.crval[i];
	}

	for (int i = 0; i < NAXIS; i++) {
		wcs->crota[i] = fit->wcs.crota[i];
	}

	for (int i = 0; i < NAXIS; i++) {
		wcs->crpix[i] = fit->wcs.crpix[i];
	}

	for (int i = 0; i < NAXIS; i++) {
		wcs->cdelt[i] = 1;
	}

	for (int i = 0; i < NAXIS; i++) {
		strcpy(wcs->ctype[i], &CTYPE[i][0]);
	}

	wcs->equinox = fit->wcs.equinox;
//	wcs->lonpole = 180;
	wcs->latpole = fit->wcs.crval[1];

	if ((status = wcsset(wcs)) != 0) {
		free_wcs();
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

	if (has_wcs()) {
		free_wcs();
	}

	ffhdr2str(fit->fptr, 1, NULL, 0, &header, &nkeyrec, &status);
	if (status) {
		report_fits_error(status);
		return FALSE;
	}

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
				wcs = (struct wcsprm*) calloc(1, sizeof(struct wcsprm));
				wcs->flag = -1;
				status = wcssub(1, prm, &nsub, axes, wcs);
				if (status == 0) {
					break;
				} else {
					siril_debug_print("wcssub error %d: %s.\n", status, wcs_errmsg[status]);
					wcsvfree(&nwcs, &wcs);
					wcs = NULL;
				}
			}
		}
	}
	wcsvfree(&nwcs, &data);
	free(header);

	if (!has_wcs()) {
		siril_debug_print("No world coordinate systems found.\n");
		wcsvfree(&nwcs, &wcs);
		wcs = NULL;
		return FALSE;
	}

	return TRUE;
#else
	return FALSE;
#endif
}

void pix2wcs(double x, double y, double *r, double *d) {
	*r = -1.0;
	*d = -1.0;
#ifdef HAVE_WCSLIB
	int status, stat[NWCSFIX];
	double imgcrd[NWCSFIX], phi, pixcrd[NWCSFIX], theta, world[NWCSFIX];

	pixcrd[0] = x;
	pixcrd[1] = y;

	status = wcsp2s(wcs, 1, 2, pixcrd, imgcrd, &phi, &theta, world, stat);
	if (status != 0)
		return;

	*r = world[0];
	*d = world[1];
#endif
}

void wcs2pix(double r, double d, double *x, double *y) {
	*x = -1.0;
	*y = -1.0;
#ifdef HAVE_WCSLIB
	int status, stat[NWCSFIX];
	double imgcrd[NWCSFIX], phi, pixcrd[NWCSFIX], theta, world[NWCSFIX];

	world[0] = r;
	world[1] = d;

	status = wcss2p(wcs, 1, 2, world, &phi, &theta, imgcrd, pixcrd, stat);
	if (status != 0)
		return;

	*x = pixcrd[0];
	*y = pixcrd[1];
#endif
}

/* get resolution in arcsec/pixel */
double get_wcs_image_resolution() {
	double resolution = -1.0;
#ifdef HAVE_WCSLIB
	if (has_wcs()) {
		double res_x = sqrt(wcs->cd[0] * wcs->cd[0] + wcs->cd[2] * wcs->cd[2]);
		double res_y = sqrt(wcs->cd[1] * wcs->cd[1] + wcs->cd[3] * wcs->cd[3]);
		resolution = (res_x + res_y) * 0.5;
	}
#endif
	if (resolution <= 0.0) {
		if (gfit.focal_length >= 0.0 && gfit.pixel_size_x >= 0.0 && gfit.pixel_size_y == gfit.pixel_size_x)
			resolution = (RADCONV / gfit.focal_length * gfit.pixel_size_x) / 3600;
	}
	return resolution;
}

double *get_wcs_crval() {
#ifdef HAVE_WCSLIB
	static double ret[NWCSFIX] = { 0 };
	if (has_wcs()) {
		for (int i = 0; i < NAXIS; i++) {
			ret[i] = wcs->crval[i];
		}
	}
	return ret;
#else
	return NULL;
#endif
}

void free_wcs() {
	// Clean up.
#ifdef HAVE_WCSLIB
	if (!wcsfree(wcs))
		free(wcs);
	wcs = NULL;
#endif
}
