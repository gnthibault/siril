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
#include "algos/statistics.h"
#include "io/conversion.h"
#include "io/image_format_fits.h"

#include "fix_xtrans_af.h"

supported_xtrans_list supported_xtrans[] =
		{
		// Camera Name      AF Pixels x,y,w,h        Sample x,y,w,h
		{ "Fujifilm X-T1",   { 1480, 997, 1972, 1313 }, { 1992, 990, 2048, 2048 } },
		{ "Fujifilm X-T2",   { 1510, 504, 3009, 3019 }, { 1992, 990, 2048, 2048 } },
		{ "Fujifilm X-T20",  { 1510, 504, 3009, 3019 }, { 1992, 990, 2048, 2048 } },
		{ "Fujifilm X-Pro2", { 1510, 504, 3009, 3019 }, { 1992, 990, 2048, 2048 } },
		{ "Fujifilm X-E3",   { 1510, 504, 3009, 3019 }, { 1992, 990, 2048, 2048 } },
		{ "Fujifilm X-H1",   { 1510, 504, 3009, 3019 }, { 1992, 990, 2048, 2048 } }
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

static void set_af_matrix(gchar *pattern, af_pixel_matrix af_matrix) {
	// If we don't find a match, we will not populate af_matrix.
	// af_pixel_matrix is [12][6].
	// Numbers are candidate green AF patterns.  G's are green.  Hyphens are red or blue.
	char matrix_str[72]="G0-G0-G3-G3---G--GG1-G1-G0-G0---G--GG2-G2-G1-G1---G--GG3-G3-G2-G2---G--G";

	// Swap with this to try the inverse.  We haven't found a sensor that works with this yet.
	//char matrix_str[72]="0G-0G-3G-3G---G--G1G-1G-0G-0G---G--G2G-2G-1G-1G---G--G3G-3G-2G-2G---G--G";

	// Start at new i in matrix_str until we find a match against pattern.
	for (int i = 0; i < 36; i += 6) {

		// Assume the pattern matches until we find it doesn't.
		int match = 1;

		// Attempt to match all 36 characters of pattern.
		for (int j = 0; j < 36; j++) {
			// If we don't match, we break.  Otherwise we keep looping.
			if ((pattern[j] == 'G' && matrix_str[j + i] == '-')
					|| (pattern[j] != 'G' && matrix_str[j + i] != '-')) {
				match = 0;
				break;
			}
		}

		// If we still match, then we found our string.
		// Create an af_matrix starting at offset = i.
		if (match) {
			for (int k = 0; k < 72; k++) {
				if (k < 72 - i) {
					af_matrix[k / 6][k % 6] = matrix_str[k + i];
				} else {
					af_matrix[k / 6][k % 6] = matrix_str[k + i - 72];
				}
			}

			// Print the matrix if debug is enabled.
			siril_debug_print("  %.6s\n", af_matrix[0]);
			siril_debug_print("  %.6s\n", af_matrix[1]);
			siril_debug_print("  %.6s\n", af_matrix[2]);
			siril_debug_print("  %.6s\n", af_matrix[3]);
			siril_debug_print("  %.6s\n", af_matrix[4]);
			siril_debug_print("  %.6s\n", af_matrix[5]);
			siril_debug_print("  %.6s\n", af_matrix[6]);
			siril_debug_print("  %.6s\n", af_matrix[7]);
			siril_debug_print("  %.6s\n", af_matrix[8]);
			siril_debug_print("  %.6s\n", af_matrix[9]);
			siril_debug_print("  %.6s\n", af_matrix[10]);
			siril_debug_print("  %.6s\n", af_matrix[11]);

			break;
		}
	}
}

// This returns the pixel type based on our AF matrix if we are within the AF rectangle.
// It returns an X if we are outside of the AF rectangle.
char get_pixel_type(rectangle af, int x, int y, af_pixel_matrix *af_matrix) {

	if (x >= af.x && x <= (af.x + af.w) && y >= af.y && y <= (af.y + af.h)) {
		// We are within the AF rectangle.
		// This is written assuming we don't know the size of the matrix.
		int matrix_cols = sizeof((*af_matrix)[0]);
		int matrix_rows = sizeof((*af_matrix)) / sizeof((*af_matrix)[0]);

		// This will return the corresponding pixel type.
		return (*af_matrix)[y % matrix_rows][x % matrix_cols];
	} else {
		// We are outside of the AF rectangle.
		return 'X';
	}
}

static int subtract_fudge(fits *fit, rectangle af, float fudge, af_pixel_matrix *af_matrix, char af_type ) {
	int width = fit->rx;
	int height = fit->ry;

	if (fit->type == DATA_USHORT) {
		WORD *buf = fit->pdata[RLAYER];
		WORD fudgew = (WORD)fudge; // truncated on purpose

		unsigned long total_fudgew = 0, total_pixels = 0;

		for (unsigned int y = 0; y < height; y++) {
			for (unsigned int x = 0; x < width; x++) {
				if (get_pixel_type(af, x, y, af_matrix) == af_type) {
					// This is an auto focus pixel.  Subtract the fudge.

					// Randomly add a 1 to some pixels to bring the average correction close to the computed value.
					WORD fudgew_rand = rand() / (float) RAND_MAX >= fudge - (float) fudgew ? fudgew : fudgew + 1;

					// Save for debugging.
					total_fudgew += fudgew_rand;
					total_pixels += 1;

					// Prevent driving the unsigned pixel negative.  Not a worry for normal use cases.
					if (fudgew_rand >= buf[x + y * width]) {
						buf[x + y * width] = 0;
					} else {
						buf[x + y * width] -= fudgew_rand;
					}
				}
			}
		}

		// Show the average integer adjustment.  Should be close to the calculated fudge by a few decimal places.
		siril_debug_print("XTRANS Integer Mean.... %.10lf\n", (double)total_fudgew/(double)total_pixels);


	} else if (fit->type == DATA_FLOAT) {
		float *buf = fit->fpdata[RLAYER];

		for (unsigned int y = 0; y < height; y++) {
			for (unsigned int x = 0; x < width; x++) {
				if (get_pixel_type(af, x, y, af_matrix) == af_type) {
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
	gboolean read_bottom_up = FALSE;

	int model = get_model(fit->instrume);
	if (model < 0) {
		siril_log_color_message(_("Fix X-Trans: Unknown camera %s, trying to read information from preferences.\n"), "red", fit->instrume);
		if (com.pref.xtrans_af.w != 0 && com.pref.xtrans_af.h != 0) {
			if (com.pref.xtrans_sample.w > fit->rx || com.pref.xtrans_sample.h > fit->ry) {
				siril_log_color_message(_("Sample box cannot be bigger than the image.\n"), "red");
				return 1;
			}
			if (com.pref.xtrans_af.w > fit->rx || com.pref.xtrans_af.h > fit->ry) {
				siril_log_color_message(_("AF box cannot be bigger than the image.\n"), "red");
				return 1;
			}
			af = com.pref.xtrans_af;
			if (com.pref.xtrans_sample.w != 0 && com.pref.xtrans_sample.h != 0) {
				sam = com.pref.xtrans_sample;
			} else {
				sam.x = 0;
				sam.y = 0;
				sam.w = fit->rx - 1;
				sam.h = fit->ry - 1;
			}
		} else {
			siril_log_color_message(_("No information available in preferences.\n"), "red");
			return 1;
		}
	} else {
		af = supported_xtrans[model].af;
		sam = supported_xtrans[model].sample;
	}


	// Flip the image so the xtrans pattern makes sense.
	// This matches logic in demosaicing.c.
	read_bottom_up = (com.pref.debayer.use_bayer_header
			&& !g_strcmp0(fit->row_order, "BOTTOM-UP"))
			|| (!com.pref.debayer.top_down);
	if (read_bottom_up) { fits_flip_top_to_bottom(fit); }


	// Struct for holding computations.
	struct af_type
	{
		// non-focus pixels
		double nfsum;
		float nfmean;
		long nfcount;

		// auto focus pixels
		double afsum;
		float afmean;
		long afcount;

		// The fudge amount to be computed.
		float fudge;

	} af_types[4] = { { 0.0, 0.f, 0L, 0.0, 0.f, 0L, 0.f }
	                , { 0.0, 0.f, 0L, 0.0, 0.f, 0L, 0.f }
	                , { 0.0, 0.f, 0L, 0.0, 0.f, 0L, 0.f }
	                , { 0.0, 0.f, 0L, 0.0, 0.f, 0L, 0.f } };

	// Variables for the winning pattern.
	float best_fudge = 0.f;
	char  best_af_type = '\0';


	// af_matrix is an RGB pattern where lowercase letters represent AF pixels and their color.
	af_pixel_matrix af_matrix = { 0 };
	set_af_matrix(fit->bayer_pattern, af_matrix);

	WORD *buf = fit->pdata[RLAYER];
	float *fbuf = fit->fpdata[RLAYER];

	if (af_matrix[0][0] == 0) {
		siril_log_color_message(_("This CFA pattern cannot be handled by fix_xtrans_ac.\n"), "red");
		return 1;
	}

	// Loop through sample rectangle and count/sum AF and non-AF pixels.
	for (unsigned int y = sam.y; y <= (sam.y + sam.h); y++) {
		for (unsigned int x = sam.x; x <= (sam.x + sam.w); x++) {
			double pixel = fit->type == DATA_FLOAT ?
							(double) fbuf[x + y * fit->rx] :
							(double) buf[x + y * fit->rx];

			switch (get_pixel_type(af, x, y, &af_matrix)) {
			case 'G': // This is a Green (non-AF) pixel.
				af_types[0].nfcount++;
				af_types[1].nfcount++;
				af_types[2].nfcount++;
				af_types[3].nfcount++;
				af_types[0].nfsum += pixel;
				af_types[1].nfsum += pixel;
				af_types[2].nfsum += pixel;
				af_types[3].nfsum += pixel;
				break;
			case '0':
				af_types[0].afcount++;
				af_types[1].nfcount++;
				af_types[2].nfcount++;
				af_types[3].nfcount++;
				af_types[0].afsum += pixel;
				af_types[1].nfsum += pixel;
				af_types[2].nfsum += pixel;
				af_types[3].nfsum += pixel;
				break;
			case '1':
				af_types[0].nfcount++;
				af_types[1].afcount++;
				af_types[2].nfcount++;
				af_types[3].nfcount++;
				af_types[0].nfsum += pixel;
				af_types[1].afsum += pixel;
				af_types[2].nfsum += pixel;
				af_types[3].nfsum += pixel;
				break;
			case '2':
				af_types[0].nfcount++;
				af_types[1].nfcount++;
				af_types[2].afcount++;
				af_types[3].nfcount++;
				af_types[0].nfsum += pixel;
				af_types[1].nfsum += pixel;
				af_types[2].afsum += pixel;
				af_types[3].nfsum += pixel;
				break;
			case '3':
				af_types[0].nfcount++;
				af_types[1].nfcount++;
				af_types[2].nfcount++;
				af_types[3].afcount++;
				af_types[0].nfsum += pixel;
				af_types[1].nfsum += pixel;
				af_types[2].nfsum += pixel;
				af_types[3].afsum += pixel;
				break;
			// case '-': // If we want to include R and B.
		 	// 	af_types[0].nfcount++;
		 	// 	af_types[1].nfcount++;
		 	// 	af_types[2].nfcount++;
		 	// 	af_types[3].nfcount++;
		 	// 	af_types[0].nfsum += pixel;
		 	// 	af_types[1].nfsum += pixel;
		 	// 	af_types[2].nfsum += pixel;
		 	// 	af_types[3].nfsum += pixel;
		 	// 	break;
			default:
				break;
			}
		}
	}

	for (int f = 0; f < 4; f++) {

		// Make sure we have a valid sample.
		if (af_types[f].nfcount == 0 || af_types[f].afcount == 0) {
			siril_log_message(_("Failed to sample enough pixels for AF type %d.\n"),f);
			return -1.f;
		}

		// Compute averages and fudge amount.
		af_types[f].nfmean = af_types[f].nfsum / af_types[f].nfcount;
		af_types[f].afmean = af_types[f].afsum / af_types[f].afcount;
		af_types[f].fudge = af_types[f].afmean - af_types[f].nfmean;

		// Debug statements.
		siril_debug_print("XTRANS %d non-AF Mean... %.10f (%ld pixels)\n", f, af_types[f].nfmean, af_types[f].nfcount);
		siril_debug_print("XTRANS %d AF Mean....... %.10f (%ld pixels)\n", f, af_types[f].afmean, af_types[f].afcount);
		siril_debug_print("XTRANS %d AF Adjust..... %.10f\n", f, af_types[f].fudge);

		// Pick the winner.
		if (fabsf(af_types[f].fudge) > fabsf(best_fudge)) {
			best_fudge = af_types[f].fudge;
			best_af_type = f+'0';  // single digit int to char
		}
	}

	// Debug for the winner.
	siril_debug_print("XTRANS Best Type %c .... %.10f\n", best_af_type, best_fudge);

	// Stay FIT, Subtract the fudge!
	subtract_fudge(fit, af, best_fudge, &af_matrix, best_af_type);

	// Flip the image back.
	if (read_bottom_up) { fits_flip_top_to_bottom(fit); }

	invalidate_stats_from_fit(fit);

	return 0;
}
