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
#include <string.h>
#include <gsl/gsl_statistics.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "algos/statistics.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"

#include "noise.h"

#define MAX_ITER 15
#define EPSILON 1E-4

static gboolean end_noise(gpointer p) {
	struct noise_data *args = (struct noise_data *) p;
	stop_processing_thread();
	set_cursor_waiting(FALSE);

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
	double norm = 1.0;
	args->retval = 0;

	if (args->verbose) {
		siril_log_color_message(_("Noise standard deviation: calculating...\n"),
				"green");
		gettimeofday(&args->t_start, NULL);
	}

	for (chan = 0; chan < args->fit->naxes[2]; chan++) {
		imstats *stat = statistics(NULL, -1, args->fit, chan, NULL, STATS_SIGMEAN, TRUE);
		if (!stat) {
			args->retval = 1;
			siril_log_message(_("Error: statistics computation failed.\n"));
			break;
		}
		args->bgnoise[chan] = stat->bgnoise;
		norm = stat->normValue;
		free_stats(stat);
	}

	if (!args->retval) {
		for (chan = 0; chan < args->fit->naxes[2]; chan++)
			if (args->fit->type == DATA_USHORT)
				siril_log_message(
						_("Background noise value (channel: #%d): %0.3f (%.3e)\n"),
						chan, args->bgnoise[chan], args->bgnoise[chan] / norm);
			else 
				siril_log_message(
						_("Background noise value (channel: #%d): %0.3f (%.3e)\n"),
						chan, args->bgnoise[chan] * USHRT_MAX_DOUBLE, args->bgnoise[chan]);
	}

	int retval = args->retval;
	if (args->use_idle)
		siril_add_idle(end_noise, args);

	return GINT_TO_POINTER(retval);
}


void evaluate_noise_in_image() {
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	set_cursor_waiting(TRUE);

	struct noise_data *args = malloc(sizeof(struct noise_data));
	args->fit = &gfit;
	args->verbose = TRUE;
	args->use_idle = TRUE;
	memset(args->bgnoise, 0.0, sizeof(double[3]));
	start_in_new_thread(noise, args);
}

