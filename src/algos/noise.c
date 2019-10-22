/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2019 team free-astro (see more in AUTHORS file)
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
#include <gsl/gsl_statistics.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "algos/statistics.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "opencv/opencv.h"

#include "noise.h"

#define MAX_ITER 15
#define EPSILON 1E-4

/*****************************************************************************
 *       N O I S E     C O M P U T A T I O N      M A N A G E M E N T        *
 ****************************************************************************/

/* Based on Jean-Luc Starck and Fionn Murtagh (1998), Automatic Noise
 * Estimation from the Multiresolution Support, Publications of the
 * Royal Astronomical Society of the Pacific, vol. 110, pp. 193–199.
 * slow algorithm. For now it is replaced by faster one. BUT, we need to keep it
 * in case we need it -. */
int backgroundnoise(fits* fit, double sigma[]) {
	int layer, k;
	fits *waveimage = calloc(1, sizeof(fits));

	if (waveimage == NULL) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	copyfits(fit, waveimage, CP_ALLOC | CP_FORMAT | CP_COPYA, 0);
	cvComputeFinestScale(waveimage);

	for (layer = 0; layer < fit->naxes[2]; layer++) {
		double sigma0, mean, norm_val;
		double epsilon = 0.0;
		WORD lo, hi;
		WORD *buf = waveimage->pdata[layer];
		unsigned int i;
		unsigned int ndata = fit->rx * fit->ry;
		g_assert(ndata > 0);

		imstats *stat = statistics(NULL, -1, waveimage, layer, NULL, STATS_BASIC);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		sigma0 = stat->sigma;
		mean = stat->mean;
		norm_val = stat->normValue;
		free_stats(stat);

		WORD *array1 = calloc(ndata, sizeof(WORD));
		WORD *array2 = calloc(ndata, sizeof(WORD));
		if (array1 == NULL || array2 == NULL) {
			PRINT_ALLOC_ERR;
			if (array1)
				free(array1);
			if (array2)
				free(array2);
			return 1;
		}
		WORD *set = array1, *subset = array2;
		memcpy(set, buf, ndata * sizeof(WORD));

		lo = round_to_WORD(LOW_BOUND * norm_val);
		hi = round_to_WORD(HIGH_BOUND * norm_val);

		sigma[layer] = sigma0;

		int n = 0;
		do {
			sigma0 = sigma[layer];
			for (i = 0, k = 0; i < ndata; i++) {
				if (set[i] >= lo && set[i] <= hi) {
					if (fabs(set[i] - mean) < 3.0 * sigma0) {
						subset[k++] = set[i];
					}
				}
			}
			ndata = k;
			sigma[layer] = gsl_stats_ushort_sd(subset, 1, ndata);
			set = subset;
			(set == array1) ? (subset = array2) : (subset = array1);
			if (ndata == 0) {
				free(array1);
				free(array2);
				siril_log_message(_("backgroundnoise: Error, no data computed\n"));
				sigma[layer] = 0.0;
				return 1;
			}
			n++;
			epsilon = fabs(sigma[layer] - sigma0) / sigma[layer];
		} while (epsilon > EPSILON && n < MAX_ITER);
		sigma[layer] *= SIGMA_PER_FWHM; // normalization
		sigma[layer] /= 0.974; // correct for 2% systematic bias
		if (n == MAX_ITER)
			siril_log_message(_("backgroundnoise: does not converge\n"));
		free(array1);
		free(array2);
	}
	clearfits(waveimage);
	invalidate_stats_from_fit(fit);

	return 0;
}

static gboolean end_noise(gpointer p) {
	struct noise_data *args = (struct noise_data *) p;
	stop_processing_thread();
	set_cursor_waiting(FALSE);
	update_used_memory();
	if (args->verbose) {
		struct timeval t_end;
		gettimeofday(&t_end, NULL);
		show_time(args->t_start, t_end);
	}
	free(args);
	return FALSE;
}

gpointer noise(gpointer p) {
	struct noise_data *args = (struct noise_data *) p;
	int chan;
	args->retval = 0;

	if (args->verbose) {
		siril_log_color_message(_("Noise standard deviation: calculating...\n"),
				"red");
		gettimeofday(&args->t_start, NULL);
	}

	for (chan = 0; chan < args->fit->naxes[2]; chan++) {
		imstats *stat = statistics(NULL, -1, args->fit, chan, NULL, STATS_NOISE);
		if (!stat) {
			args->retval = 1;
			siril_log_message(_("Error: statistics computation failed.\n"));
			break;
		}
		args->bgnoise[chan] = stat->bgnoise;
		free_stats(stat);
	}

	if (!args->retval) {
		double norm = (double) get_normalized_value(args->fit);
		for (chan = 0; chan < args->fit->naxes[2]; chan++)
			siril_log_message(
					_("Background noise value (channel: #%d): %0.3lf (%.3e)\n"), chan,
					args->bgnoise[chan], args->bgnoise[chan] / norm);
	}

	int retval = args->retval;
	if (args->use_idle)
		siril_add_idle(end_noise, args);

	return GINT_TO_POINTER(retval);
}


void on_menuitem_noise_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return;
	}

	set_cursor_waiting(TRUE);
	control_window_switch_to_tab(OUTPUT_LOGS);

	struct noise_data *args = malloc(sizeof(struct noise_data));
	args->fit = &gfit;
	args->verbose = TRUE;
	args->use_idle = TRUE;
	memset(args->bgnoise, 0.0, sizeof(double[3]));
	start_in_new_thread(noise, args);
}

