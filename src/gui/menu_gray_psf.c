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
#include "io/sequence.h"
#include "algos/star_finder.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/image_display.h"
#include "gui/PSF_list.h"

void on_menu_gray_psf_activate(GtkMenuItem *menuitem, gpointer user_data) {
	fitted_PSF *result = NULL;
	int layer = match_drawing_area_widget(com.vport[com.cvport], FALSE);

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

	popup_psf_result(result, &com.selection);
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
