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

#include <string.h>
#include <opencv2/core/version.hpp>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/statistics.h"
#include "core/processing.h"
#include "core/undo.h"
#include "algos/colors.h"
#include "opencv/opencv.h"
#include "gui/image_display.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"
#include "gui/callbacks.h"
#include "io/single_image.h"
#include "gui/message_dialog.h"
#include "gui/preview_timer.h"

#include "clahe.h"

static double clahe_limit_value;
static int clahe_tile_size;
static gboolean clahe_show_preview;

static void clahe_startup() {
	copy_gfit_to_backup();
	clahe_limit_value = 2.0;
	clahe_tile_size = 8;
}

static void clahe_close(gboolean revert) {
	set_cursor_waiting(TRUE);
	if (revert) {
		copy_backup_to_gfit();
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	} else {
		invalidate_stats_from_fit(&gfit);
		undo_save_state(get_preview_gfit_backup(),
				"Processing: CLAHE (size=%d, clip=%.2f)", clahe_tile_size, clahe_limit_value);
	}
	clear_backup();
	set_cursor_waiting(FALSE);
}

static int clahe_update_preview() {
	copy_backup_to_gfit();
	if (CV_MAJOR_VERSION < 3) {
		char *error = siril_log_message(_("Your version of opencv is "
				"too old for this feature. Please upgrade your system."));
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Upgrade your system"),
				error);
		return 1;
	}

	struct CLAHE_data *args = malloc(sizeof(struct CLAHE_data));

	set_cursor_waiting(TRUE);

	args->fit = &gfit;
	args->clip = clahe_limit_value;
	args->tileSize = clahe_tile_size;

	start_in_new_thread(clahe, args);
	return 0;
}

gpointer clahe(gpointer p) {
	struct CLAHE_data *args = (struct CLAHE_data*) p;

	cvClahe(args->fit, args->clip, args->tileSize);

	siril_add_idle(end_generic, args);
	return GINT_TO_POINTER(0);
}

void apply_clahe_cancel() {
	clahe_close(TRUE);
	siril_close_dialog("CLAHE_dialog");
}

/** callbacks **/

void on_menuitem_clahe_activate(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("CLAHE_dialog");
}

void on_clahe_cancel_clicked(GtkMenuItem *menuitem, gpointer user_data) {
	clahe_close(TRUE);
	siril_close_dialog("CLAHE_dialog");
}

void on_CLAHE_dialog_close(GtkDialog *dialog, gpointer user_data) {
	clahe_close(TRUE);
}

void on_clahe_undo_clicked(GtkButton *button, gpointer user_data) {
	GtkSpinButton *liit_value = GTK_SPIN_BUTTON(lookup_widget("spin_clahe"));
	GtkSpinButton *tiles_size = GTK_SPIN_BUTTON(lookup_widget("clahe_tiles_size_spin"));
	clahe_limit_value = 2.0;
	clahe_tile_size = 8;

	set_notify_block(TRUE);
	gtk_spin_button_set_value(liit_value, clahe_limit_value);
	gtk_spin_button_set_value(tiles_size, clahe_tile_size);
	set_notify_block(FALSE);

	copy_backup_to_gfit();

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = clahe_show_preview;
	notify_update((gpointer) param);
}

void on_clahe_Apply_clicked(GtkButton *button, gpointer user_data) {
	if (clahe_show_preview == FALSE) {
		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = clahe_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}

	clahe_close(FALSE);
	siril_close_dialog("CLAHE_dialog");
}

void on_CLAHE_dialog_show(GtkWidget *widget, gpointer user_data) {
	clahe_startup();

	set_notify_block(TRUE);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("clahe_tiles_size_spin")), clahe_tile_size);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spin_clahe")),	clahe_limit_value);
	set_notify_block(FALSE);

	clahe_show_preview = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("clahe_preview")));

	/* default parameters transform image, we need to update preview */
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = clahe_show_preview;
	notify_update((gpointer) param);
}

void on_clahe_preview_toggled(GtkToggleButton *button, gpointer user_data) {
	if (clahe_show_preview == TRUE) {
		copy_backup_to_gfit();
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	} else {
		copy_gfit_to_backup();

		update_image *param = malloc(sizeof(update_image));
		param->update_preview_fn = clahe_update_preview;
		param->show_preview = TRUE;
		notify_update((gpointer) param);
	}
	clahe_show_preview = !clahe_show_preview;
}

/** adjusters **/
void on_spin_clahe_value_changed(GtkSpinButton *button, gpointer user_data) {
	clahe_limit_value = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = clahe_show_preview;
	notify_update((gpointer) param);
}

void on_clahe_tiles_size_spin_value_changed(GtkSpinButton *button, gpointer user_data) {
	clahe_tile_size = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = clahe_update_preview;
	param->show_preview = clahe_show_preview;
	notify_update((gpointer) param);
}
