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
#include "gui/progress_and_log.h"
#include "gui/image_display.h"
#include "io/single_image.h"

#include "preview_timer.h"

#define PREVIEW_DELAY 200

static guint timer_id;
static gboolean notify_is_blocked;

static gboolean update_preview(gpointer user_data) {
	update_image *im = (update_image*) user_data;

	if (notify_is_blocked == TRUE) return FALSE;

	set_cursor_waiting(TRUE);

	im->update_preview_fn();

	printf("update preview\n");
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	return FALSE;
}

static void free_struct(gpointer user_data) {
	update_image *im = (update_image*) user_data;

	timer_id = 0;
	free(im);
}

void set_notify_block(gboolean value) {
	notify_is_blocked = value;
}

void notify_update(gpointer user_data) {
	if (timer_id != 0) {
		g_source_remove(timer_id);
	}
	timer_id = g_timeout_add_full(G_PRIORITY_DEFAULT_IDLE,
			PREVIEW_DELAY, (GSourceFunc) update_preview, user_data,
			(GDestroyNotify) free_struct);
}
