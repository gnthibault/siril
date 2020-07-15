/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"
#include "algos/statistics.h"

#include "fix_xtrans_af.h"

supported_xtrans_list supported_xtrans[] =
		{
		// Camera Name      AF Pixels x,y,w,h        Sample x,y,w,h
		{ "Fujifilm X-T2", { 1510, 507, 3009, 3016 }, { 1992, 990,2048, 2048 } },
		{ "Fujifilm X-T20", { 1510, 507, 3009, 3016 }, { 1992, 990,2048, 2048 } }
};

static int get_nb_xtrans_supported() {
	return G_N_ELEMENTS(supported_xtrans);
}

static int get_model(const char *model) {
	int nb_xtrans = get_nb_xtrans_supported();
	for (int i = 0; i < nb_xtrans; i++) {
		if (!g_ascii_strcasecmp(model, supported_xtrans[i].model))
			return i;
	}
	return -1;
}

// This returns true if the pixel is a special auto focus pixel.
static int is_af_pixel(rectangle af, int x, int y) {
	// x=0, y=0 is the bottom left corner of the image.
	return (x >= af.x && x <= (af.x + af.w) && y >= af.y && y <= (af.y + af.h)
			&& (x + 2) % 3 == 0 && ((y + 5) % 12 == 0 || (y + 9) % 12 == 0));
}

static int subtract_fudge(fits *fit, rectangle af, float fudge) {
	int width = fit->rx;
	int height = fit->ry;

	if (fit->type == DATA_USHORT) {
		WORD *buf = fit->pdata[RLAYER];

		for (unsigned int y = 0; y < height; y++) {
			for (unsigned int x = 0; x < width; x++) {
				if (is_af_pixel(af, x, y)) {
					// This is an auto focus pixel.  Subtract the fudge.
					buf[x + y * width] -= roundf_to_WORD(fudge);
				}
			}
		}
	} else if (fit->type == DATA_FLOAT) {
		float *buf = fit->fpdata[RLAYER];

		for (unsigned int y = 0; y < height; y++) {
			for (unsigned int x = 0; x < width; x++) {
				if (is_af_pixel(af, x, y)) {
					// This is an auto focus pixel.  Subtract the fudge.
					buf[x + y * width] -= fudge;
				}
			}
		}
	}
	// the caller should call invalidate_stats_from_fit(fit);
	return 0;
}

int fix_xtrans_ac(fits *fit) {
	rectangle af, sam;

	int model = get_model(fit->instrume);
	if (model < 0) {
		siril_log_color_message(_("Fix X-Trans: Unknown camera %s, trying to read information from preferences.\n"), "red", fit->instrume);
		if (com.pref.xtrans_af.w != 0 && com.pref.xtrans_af.h != 0 && com.pref.xtrans_sample.w != 0 && com.pref.xtrans_sample.h != 0) {
			af = com.pref.xtrans_af;
			sam = com.pref.xtrans_sample;
		} else {
			siril_log_color_message(_("No information available in preferences.\n"), "red");
			return 1;
		}
	} else {
		af = supported_xtrans[model].af;
		sam = supported_xtrans[model].sample;
	}

	// non-focus pixels
	double nfsum = 0.0;
	float nfmean = 0.f;
	long nfcount = 0L;

	// auto focus pixels
	double afsum = 0.0;
	float afmean = 0.f;
	long afcount = 0L;

	// The fudge amount to apply to auto focus pixels. (computed)
	float fudge;

	WORD *buf = fit->pdata[RLAYER];
	float *fbuf = fit->fpdata[RLAYER];

	// Loop through sample rectangle and count/sum AF and non-AF pixels.
	for (unsigned int y = sam.y; y <= (sam.y + sam.h); y++) {
		for (unsigned int x = sam.x; x <= (sam.x + sam.w); x++) {
			float pixel =
					fit->type == DATA_FLOAT ?
							fbuf[x + y * fit->rx] :
							(float) buf[x + y * fit->rx];
			if (is_af_pixel(af, x, y)) {
				// This is an AF pixel.
				afcount++;
				afsum += (double) pixel;
			} else {
				// This is not an AF pixel.
				nfcount++;
				nfsum += (double) pixel;
			}
		}
	}

	// Make sure we have a valid sample.
	if (nfcount == 0 || afcount == 0) {
		siril_log_message(_("Failed to sample enough pixels.\n"));
		return -1.f;
	}

	// Compute averages and fudge amount.
	nfmean = nfsum / nfcount;
	afmean = afsum / afcount;
	fudge = afmean - nfmean;

	// Debug statements.  Remove later???
	siril_log_message(_("Reg Pixel Mean...   %.10f, Count: %ld\n"), nfmean, nfcount);
	siril_log_message(_("AF Pixel Mean....   %.10f, Count: %ld\n"), afmean, afcount);
	siril_log_message(_("AF Pixel Adjust..   %.10f\n"), fudge);

	// Stay FIT, Subtract the fudge!
	subtract_fudge(fit, af, fudge);

	invalidate_stats_from_fit(fit);

	return 0;
}
