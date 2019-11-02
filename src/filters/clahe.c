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

#include <opencv/opencv.h>
#include <opencv2/core/version.hpp>

#include "core/siril.h"
#include "core/proto.h"
#include "core/undo.h"
#include "core/processing.h"
#include "algos/colors.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"
#include "gui/callbacks.h"
#include "io/single_image.h"
#include "gui/message_dialog.h"

#include "clahe.h"

void on_menuitem_clahe_activate(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("CLAHE_dialog");
}

void on_clahe_cancel_clicked(GtkMenuItem *menuitem, gpointer user_data) {
	siril_close_dialog("CLAHE_dialog");
}

void on_clahe_Apply_clicked(GtkButton *button, gpointer user_data) {
	GtkRange *clip;
	GtkSpinButton *size;

	if (CV_MAJOR_VERSION < 3) {
		char *error = siril_log_message(_("Your version of opencv is "
				"too old for this feature. Please upgrade your system."));
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Upgrade your system"), error);
		return;
	}

	clip = GTK_RANGE(gtk_builder_get_object(builder, "scale_clahe"));
	size = GTK_SPIN_BUTTON(lookup_widget("clahe_tiles_size_spin"));

	struct CLAHE_data *args = malloc(sizeof(struct CLAHE_data));

	set_cursor_waiting(TRUE);

	args->fit = &gfit;
	args->clip = gtk_range_get_value(clip);
	args->tileSize = gtk_spin_button_get_value_as_int(size);

	undo_save_state(&gfit, "Processing: CLAHE (size=%d, clip=%.2f)",
			args->tileSize, args->clip);
	start_in_new_thread(clahe, args);
}

gpointer clahe(gpointer p) {
	struct CLAHE_data *args = (struct CLAHE_data *) p;
	struct timeval t_start, t_end;

	char *msg = siril_log_color_message(_("CLAHE: processing...\n"), "red");
	msg[strlen(msg) - 1] = '\0';
	set_progress_bar_data(msg, PROGRESS_PULSATE);
	gettimeofday(&t_start, NULL);

	cvClahe(args->fit, args->clip, args->tileSize);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	siril_add_idle(end_generic, args);
	set_progress_bar_data(_("CLAHE applied"), PROGRESS_DONE);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return GINT_TO_POINTER(0);
}
