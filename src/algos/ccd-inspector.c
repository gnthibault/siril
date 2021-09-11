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

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/command.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/statistics.h"
#include "algos/sorting.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/progress_and_log.h"

#include "ccd-inspector.h"

static void draw_polygon(float rx, float ry, float m1, float m2, float m3, float m4, float mcentre) {
	float r1, r2, r3, r4;
	pointf c = { rx / 2.f, ry / 2.f };
	float m = (m1 + m2 + m3 + m4) / 4.f;
	float diag = sqrtf(rx * rx + ry * ry) / 4.f;

	/* now we compute the four radius. */
	r1 = diag * (((m1 - m) / m) + 1);
	r2 = diag * (((m2 - m) / m) + 1);
	r3 = diag * (((m3 - m) / m) + 1);
	r4 = diag * (((m4 - m) / m) + 1);

	com.tilt = malloc(sizeof(sensor_tilt));

	com.tilt->pt[0].x = c.x + (r1 * sin(7.0 * M_PI / 4.0));
	com.tilt->pt[0].y = ry - (c.y + (r1 * cos(7.0 * M_PI / 4.0)));
	com.tilt->fwhm[0] = m1;

	com.tilt->pt[1].x = c.x + (r2 * sin(M_PI / 4));
	com.tilt->pt[1].y = ry - (c.y + (r2 * cos(M_PI / 4.0)));
	com.tilt->fwhm[1] = m2;

	com.tilt->pt[2].x = c.x + (r3 * sin(5.0 * M_PI / 4));
	com.tilt->pt[2].y = ry - (c.y + (r3 * cos(5.0 * M_PI / 4.0)));
	com.tilt->fwhm[2] = m3;

	com.tilt->pt[3].x = c.x + (r4 * sin(3.0 * M_PI / 4));
	com.tilt->pt[3].y = ry - (c.y + (r4 * cos(3.0 * M_PI / 4.0)));
	com.tilt->fwhm[3] = m4;

	com.tilt->fwhm_centre = mcentre;

	redraw(com.cvport, REMAP_NONE);
}

void clear_sensor_tilt() {
	free(com.tilt);
	com.tilt = NULL;
}

int draw_sensor_tilt(fits *fit) {
	int i1 = 0, i2 = 0, i3 = 0, i4 = 0, ir1 = 0, ir2 = 0;
	int nbstars = 0;
	float best, worst;
	pointf center = {fit->rx / 2.f, fit->ry / 2.f };
	float r = sqrtf(center.x * center.x + center.y * center.y);
	float r1 = 0.25f * r;
	float r2 = 0.75f * r;
	int layer = com.cvport == RGB_VPORT ? GLAYER : com.cvport;

	delete_selected_area();

	psf_star **stars = peaker(fit, layer, &com.starfinder_conf, &nbstars, NULL, FALSE, FALSE);

	float *f = malloc(nbstars * sizeof(float));

	float *f1 = calloc(nbstars, sizeof(float));
	float *f2 = calloc(nbstars, sizeof(float));
	float *f3 = calloc(nbstars, sizeof(float));
	float *f4 = calloc(nbstars, sizeof(float));

	float *fr1 = calloc(nbstars, sizeof(float));
	float *fr2 = calloc(nbstars, sizeof(float));

	for (int i = 0; i < nbstars; i++) {
		float x = (float) stars[i]->xpos;
		float y = (float) stars[i]->ypos;

		/* global */
		f[i] = (float) (stars[i]->fwhmx + stars[i]->fwhmy) * 0.5f;
		/* check for 4 areas */
		if ((x < center.x) && (y < center.y)) {
			f1[i1++] = (float) (stars[i]->fwhmx + stars[i]->fwhmy) * 0.5f;
		} else if ((x > center.x) && (y < center.y)) {
			f2[i2++] = (float) (stars[i]->fwhmx + stars[i]->fwhmy) * 0.5f;
		} else if ((x < center.x) && (y > center.y)) {
			f3[i3++] = (float) (stars[i]->fwhmx + stars[i]->fwhmy) * 0.5f;
		} else if ((x > center.x) && (y > center.y)) {
			f4[i4++] = (float) (stars[i]->fwhmx + stars[i]->fwhmy) * 0.5f;
		}
		/* check for off-axis aberration */
		if (((x - center.x) * (x - center.x) + (y - center.y) * (y - center.y)) < (r1 * r1)) {
			fr1[ir1++] = (float) (stars[i]->fwhmx + stars[i]->fwhmy) * 0.5f;
		} else if (((x - center.x) * (x - center.x) + (y - center.y) * (y - center.y)) > (r2 * r2)) {
			fr2[ir2++] = (float) (stars[i]->fwhmx + stars[i]->fwhmy) * 0.5f;
		}
	}

	if ((i1 != 0) && (i2 != 0) && (i3 != 0) && (i4 != 0) && (ir1 != 0) && (ir2 != 0)) {
		quicksort_f(f, nbstars);
		quicksort_f(f1, i1);
		quicksort_f(f2, i2);
		quicksort_f(f3, i3);
		quicksort_f(f4, i4);
		quicksort_f(fr1, ir1);
		quicksort_f(fr2, ir2);

		float m = siril_stats_trmean_from_sorted_data(0.25, f, 1, nbstars);
		float m1 = siril_stats_trmean_from_sorted_data(0.25, f1, 1, i1);
		float m2 = siril_stats_trmean_from_sorted_data(0.25, f2, 1, i2);
		float m3 = siril_stats_trmean_from_sorted_data(0.25, f3, 1, i3);
		float m4 = siril_stats_trmean_from_sorted_data(0.25, f4, 1, i4);

		float mr1 = siril_stats_trmean_from_sorted_data(0.25, fr1, 1, ir1);
		float mr2 = siril_stats_trmean_from_sorted_data(0.25, fr2, 1, ir2);

		best = min(min(m1, m2), min(m3, m4));
		worst = max(max(m1, m2), max(m3, m4));

		draw_polygon((float) fit->rx, (float) fit->ry, m1, m2, m3, m4, mr1);
		siril_log_message(_("Stars: %d, Truncated mean[FWHM]: %.2f, Sensor tilt[FWHM]: %.2f, Off-axis aberration[FWHM]: %.2f\n"), nbstars, m, worst - best, mr2 - mr1);
	}

	free(f);

	free(f1);
	free(f2);
	free(f3);
	free(f4);

	free(fr1);
	free(fr2);

	free_fitted_stars(stars);

	return 0;
}

void on_tilt_button_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	confirm_peaker_GUI();
	draw_sensor_tilt(&gfit);
	set_cursor_waiting(FALSE);
}
