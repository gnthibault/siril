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

#include "siril_preview.h"

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/undo.h"
#include "gui/progress_and_log.h"
#include "gui/image_display.h"
#include "io/single_image.h"


#define PREVIEW_DELAY 200

static guint timer_id;
static gboolean notify_is_blocked;
static gboolean preview_is_active;
static fits preview_gfit_backup;

static gboolean update_preview(gpointer user_data) {
	update_image *im = (update_image*) user_data;

	if (notify_is_blocked == TRUE) return FALSE;

	if (im->show_preview) {
		siril_debug_print("update preview\n");
		set_cursor_waiting(TRUE);
		im->update_preview_fn();
	}

	waiting_for_thread(); // in case function is run in another thread
	set_progress_bar_data(NULL, PROGRESS_DONE);
	if (im->show_preview) {
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
		set_cursor_waiting(FALSE);
	}
	return FALSE;
}

static void free_struct(gpointer user_data) {
	update_image *im = (update_image*) user_data;

	timer_id = 0;
	free(im);
}

void copy_gfit_to_backup() {
	copyfits(&gfit, &preview_gfit_backup, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	preview_is_active = TRUE;
}

void copy_backup_to_gfit() {
	copyfits(&preview_gfit_backup, &gfit, CP_COPYA, -1);
}

fits *get_preview_gfit_backup() {
	return &preview_gfit_backup;
}

gboolean is_preview_active() {
	return preview_is_active;
}

void clear_backup() {
	clearfits(&preview_gfit_backup);
	preview_is_active = FALSE;
}

void set_notify_block(gboolean value) {
	notify_is_blocked = value;
}

void siril_preview_hide() {
	copy_backup_to_gfit();
	clear_backup();
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
}

void notify_update(gpointer user_data) {
	if (timer_id != 0) {
		g_source_remove(timer_id);
	}
	timer_id = g_timeout_add_full(G_PRIORITY_DEFAULT_IDLE,
			PREVIEW_DELAY, (GSourceFunc) update_preview, user_data,
			(GDestroyNotify) free_struct);
}
