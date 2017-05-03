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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/siril.h"
#include "gui/callbacks.h"
#include "core/proto.h"
#include "registration/registration.h"
#include "io/sequence.h"
#include "io/single_image.h"

#define REGLAYER 0

static sequence *seq = NULL;		// the sequence of channels
static struct registration_method *reg_methods[3];


static int extract_fits() {
	copyfits(&gfit, wfit + 0, CP_ALLOC | CP_EXTRACT, RLAYER); /* r */
	copyfits(&gfit, wfit + 1, CP_ALLOC | CP_EXTRACT, GLAYER); /* g */
	copyfits(&gfit, wfit + 2, CP_ALLOC | CP_EXTRACT, BLAYER); /* b */

	return 0;
}

static void initialize_methods() {
	reg_methods[0] = new_reg_method(_("One star registration (deep-sky)"),
			&register_shift_fwhm, REQUIRES_ANY_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[1] = new_reg_method(_("Image pattern alignment (planetary/deep-sky)"),
			&register_shift_dft, REQUIRES_SQUARED_SELECTION, REGTYPE_PLANETARY);

	reg_methods[2] = NULL;
}

static int initialize_internal_rgb_sequence() {
	int i;
	if (seq) free_sequence(seq, TRUE);

	extract_fits();
	seq = create_internal_sequence(3);
	for (i = 0; i < 3; i++) {
		internal_sequence_set(seq, i, &wfit[i]);
	}
	seq->rx = gfit.rx;
	seq->ry = gfit.ry;

	return 0;
}

static void align_and_compose() {
	int x, y, i = 0;	// i is browsing the 1D buffer, i = y * rx + x
	for (y = 0; y < gfit.ry; ++y) {
		for (x = 0; x < gfit.rx; ++x) {
			int channel;
			for (channel = 0; channel < 3; channel++) {
				if (seq && seq->regparam) {
					WORD pixel;
					int realX = x - seq->regparam[REGLAYER][channel].shiftx;
					int realY = y - seq->regparam[REGLAYER][channel].shifty;
					if (realX < 0 || realX >= gfit.rx)
						pixel = 0;
					else if (realY < 0 || realY >= gfit.ry)
						pixel = 0;
					else
						pixel = wfit[channel].pdata[0][realX + gfit.rx * realY];
					gfit.pdata[channel][i] = pixel;
				}
			}
			i++;
		}
	}
}

int rgb_align(int m) {
	struct registration_args regargs;
	struct registration_method *method;

	initialize_methods();
	initialize_internal_rgb_sequence();

	/* align it */
	method = reg_methods[m];

	regargs.seq = seq;
	regargs.seq->nb_layers = 3;
	regargs.process_all_frames = TRUE;
	get_the_registration_area(&regargs, method);
	regargs.layer = REGLAYER;
	regargs.run_in_thread = FALSE;

	set_cursor_waiting(TRUE);
	set_progress_bar_data(NULL, PROGRESS_RESET);
	if (method->method_ptr(&regargs))
		set_progress_bar_data(_("Error in layers alignment."), PROGRESS_DONE);
	else
		set_progress_bar_data(_("Registration complete."), PROGRESS_DONE);
	set_cursor_waiting(FALSE);

	align_and_compose();

	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	return 0;
}
