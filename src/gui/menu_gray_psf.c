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

#include "algos/PSF.h"
#include "core/command.h"
#include "io/sequence.h"
#include "algos/star_finder.h"

#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/image_display.h"
#include "gui/PSF_list.h"

void on_menu_gray_psf_activate(GtkMenuItem *menuitem, gpointer user_data) {
	gchar *msg;
	fitted_PSF *result = NULL;
	int layer = match_drawing_area_widget(com.vport[com.cvport], FALSE);
	char *str;

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
	msg = g_strdup_printf(_("Centroid Coordinates:\n\t\tx0=%.2fpx\n\t\ty0=%.2fpx\n\n"
				"Full Width Half Maximum:\n\t\tFWHMx=%.2f%s\n\t\tFWHMy=%.2f%s\n\n"
				"Angle:\n\t\t%0.2fdeg\n\n"
				"Background Value:\n\t\tB=%.6f\n\n"
				"Maximal Intensity:\n\t\tA=%.6f\n\n"
				"Magnitude (%s):\n\t\tm=%.4f\u00B1%.4f\n\n"
				"RMSE:\n\t\tRMSE=%.3e"), result->x0 + com.selection.x,
			com.selection.y + com.selection.h - result->y0, result->fwhmx,
			result->units, result->fwhmy, result->units, result->angle, result->B,
			result->A, str, result->mag + com.magOffset, result->s_mag, result->rmse);
	show_data_dialog(msg, "PSF Results");
	g_free(msg);
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

void on_menu_gray_pick_star_activate(GtkMenuItem *menuitem, gpointer user_data) {
	int layer = match_drawing_area_widget(com.vport[com.cvport], FALSE);
	int new_index;

	if (layer != -1) {
		if (!(com.selection.h && com.selection.w))
			return;
		if (com.selection.w > 300 || com.selection.h > 300) {
			siril_message_dialog(GTK_MESSAGE_WARNING, _("Current selection is too large"),
					_("To determine the PSF, please make a selection around a star."));
			return;
		}
		fitted_PSF *new_star = add_star(&gfit, layer, &new_index);
		if (new_star) {
			add_star_to_list(new_star);
			siril_open_dialog("stars_list_window");
		} else
			return;
	}
	redraw(com.cvport, REMAP_NONE);
}
