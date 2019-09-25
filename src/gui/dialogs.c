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

#include "core/siril.h"
#include "core/proto.h"
#include "gui/callbacks.h"
#include "gui/histogram.h"
#include "filters/asinh.h"
#include "filters/saturation.h"
#include "filters/wavelets.h"

#include "gui/dialogs.h"

static const SirilDialogEntry entries[] =
{
		{"asinh_dialog", IMAGE_PROCESSING_DIALOG, TRUE, apply_asinh_cancel},
		{"background_extraction_dialog", IMAGE_PROCESSING_DIALOG, FALSE, NULL},
		{"canon_fixbanding_dialog", IMAGE_PROCESSING_DIALOG, FALSE, NULL},
		{"CLAHE_dialog", IMAGE_PROCESSING_DIALOG, FALSE, NULL},
		{"composition_dialog", IMAGE_PROCESSING_DIALOG, FALSE, NULL},
		{"color_calibration", IMAGE_PROCESSING_DIALOG, FALSE, NULL},
		{"cosmetic_dialog", IMAGE_PROCESSING_DIALOG, FALSE, NULL},
		{"crop_dialog", IMAGE_PROCESSING_DIALOG, FALSE, NULL},
		{"deconvolution_dialog", IMAGE_PROCESSING_DIALOG, FALSE, NULL},
		{"dialog_FFT", IMAGE_PROCESSING_DIALOG, FALSE, NULL},
		{"extract_channel_dialog", OTHER_DIALOG, FALSE, NULL},
		{"extract_wavelets_layers_dialog", OTHER_DIALOG, FALSE, NULL},
		{"file_information", INFORMATION_DIALOG, FALSE, NULL},
		{"histogram_window", IMAGE_PROCESSING_DIALOG, TRUE, apply_histo_cancel},
		{"ImagePlateSolver_Dial", IMAGE_PROCESSING_DIALOG, FALSE, NULL},
		{"Median_dialog", IMAGE_PROCESSING_DIALOG, FALSE, NULL},
		{"resample_dialog", IMAGE_PROCESSING_DIALOG, FALSE, NULL},
		{"rgradient_dialog", IMAGE_PROCESSING_DIALOG, FALSE, NULL},
		{"rotation_dialog", IMAGE_PROCESSING_DIALOG, FALSE, NULL},
		{"satu_dialog", IMAGE_PROCESSING_DIALOG, TRUE, apply_satu_cancel},
		{"SCNR_dialog", IMAGE_PROCESSING_DIALOG,  FALSE, NULL},
		{"settings_window", INFORMATION_DIALOG, FALSE, NULL},
		{"split_cfa_dialog", OTHER_DIALOG, FALSE, NULL},
		{"stars_list_window", INFORMATION_DIALOG, FALSE, NULL},
		{"StatWindow", INFORMATION_DIALOG, FALSE, NULL},
		{"wavelets_dialog", IMAGE_PROCESSING_DIALOG, TRUE, apply_wavelets_cancel}
};

static SirilDialogEntry get_entry_by_id(gchar *id) {
	int i;
	SirilDialogEntry retvalue = {NULL, NO_DIALOG, FALSE, NULL};

	for (i = 0; i < G_N_ELEMENTS(entries); i++) {
		if (!g_ascii_strcasecmp(id, entries[i].identifier)) {
			return entries[i];
		}
	}
	return retvalue;
}

void siril_open_dialog(gchar *id) {
	int i;
	gint x, y;
	gboolean win_already_shown = FALSE;

	if (get_entry_by_id(id).type != INFORMATION_DIALOG) {
		for (i = 0; i < G_N_ELEMENTS(entries); i++) {
			GtkWidget *w = lookup_widget(entries[i].identifier);
			if (gtk_widget_get_visible(w) && entries[i].type != INFORMATION_DIALOG) {
				win_already_shown = TRUE;
				gtk_window_get_position(GTK_WINDOW(w), &x, &y);
				if (entries[i].has_preview)
					entries[i].apply_function();
				gtk_widget_hide(w);
				break;
			}
		}
	}
	GtkWindow *win = GTK_WINDOW(lookup_widget(id));
	if (win_already_shown && x >=0 && y >= 0) {
		gtk_window_move(win, x, y);
	} else {
		gtk_window_set_position (win, GTK_WIN_POS_CENTER_ON_PARENT);
	}
	gtk_window_set_transient_for(win, GTK_WINDOW(lookup_widget("main_window")));
	gtk_window_present_with_time(win, GDK_CURRENT_TIME);
}

void siril_close_dialog(gchar *id) {
	gtk_widget_hide(lookup_widget(id));
}

void siril_close_preview_dialogs() {
	int i;
	for (i = 0; i < G_N_ELEMENTS(entries); i++) {
		GtkWidget *w = lookup_widget(entries[i].identifier);
		if (gtk_widget_get_visible(w) && (entries[i].has_preview)) {
			entries[i].apply_function();
			gtk_widget_hide(w);
		}
	}

}
