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

#include <math.h>
#include "core/siril.h"
#include "core/proto.h"

#include "algos/PSF.h"
#include "core/command.h"
#include "core/siril_world_cs.h"
#include "io/sequence.h"
#include "algos/star_finder.h"
#include "algos/siril_wcs.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/image_display.h"
#include "gui/PSF_list.h"

static gchar *build_wcs_url(gchar *ra, gchar *dec) {
	if (!has_wcs(&gfit)) return NULL;

	double resolution = get_wcs_image_resolution(&gfit);

	gchar *tol = g_strdup_printf("%lf", resolution * 3600 * 15);

	GString *url = g_string_new("https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=");
	url = g_string_append(url, ra);
	url = g_string_append(url, dec);
	url = g_string_append(url, "&Radius=");
	url = g_string_append(url, tol);
	url = g_string_append(url, "&Radius.unit=arcsec");
	url = g_string_append(url, "#lab_basic");

	gchar *simbad_url = g_string_free(url, FALSE);
	gchar *cleaned_url = url_cleanup(simbad_url);

	g_free(tol);
	g_free(simbad_url);

	return cleaned_url;
}

void on_menu_gray_psf_activate(GtkMenuItem *menuitem, gpointer user_data) {
	gchar *msg, *coordinates, *url = NULL;
	fitted_PSF *result = NULL;
	int layer = match_drawing_area_widget(com.vport[com.cvport], FALSE);
	const char *str;

	if (layer == -1)
		return;
	if (!(com.selection.h && com.selection.w))
		return;
	if (com.selection.w > 300 || com.selection.h > 300) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Current selection is too large"),
				_("To determine the PSF, please make a selection around a star."));

		return;
	}
	result = psf_get_minimisation(&gfit, layer, &com.selection, TRUE, TRUE, TRUE);
	if (!result)
		return;

	if (com.magOffset > 0.0)
		str = "true reduced";
	else
		str = "relative";

	double x = result->x0 + com.selection.x;
	double y = com.selection.y + com.selection.h - result->y0;
	if (has_wcs(&gfit)) {
		double world_x, world_y;
		SirilWorldCS *world_cs;

		pix2wcs(&gfit, x, (double) gfit.ry - y, &world_x, &world_y);
		world_cs = siril_world_cs_new_from_a_d(world_x, world_y);
		if (world_cs) {
			gchar *ra = siril_world_cs_alpha_format(world_cs, "%02d %02d %.3lf");
			gchar *dec = siril_world_cs_delta_format(world_cs, "%c%02d %02d %.3lf");

			url = build_wcs_url(ra, dec);

			g_free(ra);
			g_free(dec);

			ra = siril_world_cs_alpha_format(world_cs, " %02dh%02dm%02ds");
			dec = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%02d\"");

			coordinates = g_strdup_printf("x0=%.2fpx\t%s J2000\n\t\ty0=%.2fpx\t%s J2000", x, ra, y, dec);

			g_free(ra);
			g_free(dec);
			siril_world_cs_unref(world_cs);
		} else {
			coordinates = g_strdup_printf("x0=%.2fpx\n\t\ty0=%.2fpx", x, y);
		}
	} else {
		coordinates = g_strdup_printf("x0=%.2fpx\n\t\ty0=%.2fpx", x, y);
	}

	double fwhmx, fwhmy;
	char *units;
	get_fwhm_as_arcsec_if_possible(result, &fwhmx, &fwhmy, &units);
	msg = g_strdup_printf(_("Centroid Coordinates:\n\t\t%s\n\n"
				"Full Width Half Maximum:\n\t\tFWHMx=%.2f%s\n\t\tFWHMy=%.2f%s\n\n"
				"Angle:\n\t\t%0.2fdeg\n\n"
				"Background Value:\n\t\tB=%.6f\n\n"
				"Maximal Intensity:\n\t\tA=%.6f\n\n"
				"Magnitude (%s):\n\t\tm=%.4f\u00B1%.4f\n\n"
				"RMSE:\n\t\tRMSE=%.3e"),
			coordinates, fwhmx, units, fwhmy, units, result->angle, result->B,
			result->A, str, result->mag + com.magOffset, result->s_mag, result->rmse);
	show_data_dialog(msg, "PSF Results", NULL, url);
	g_free(coordinates);
	g_free(msg);
	g_free(url);
	free(result);
}

void on_menu_gray_seqpsf_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (!sequence_is_loaded()) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("PSF for the sequence only applies on sequences"),
				_("Please load a sequence before trying to apply the PSF for the sequence."));
	} else {
		process_seq_psf(0);
	}
}
