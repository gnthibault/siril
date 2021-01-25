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
#include "core/siril_world_cs.h"
#include "gui/dialogs.h"
#include "gui/utils.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "algos/plateSolver.h"
#include "algos/annotate.h"
#include "algos/siril_wcs.h"

static gboolean parse_buffer(char *buffer) {
	char **token, **fields, *realname= NULL;
	point center;
	int nargs, i = 0;
	SirilWorldCS *world_cs = NULL;

	token = g_strsplit(buffer, "\n", -1);
	nargs = g_strv_length(token);

	while (i < nargs) {
		if (g_str_has_prefix(token[i], "%I.0 ")) {
				gchar *name = g_strstr_len(token[i], strlen(token[i]), "%I.0 ");
				realname = g_strdup(name + 5);

		} else if (g_str_has_prefix (token[i], "%J ")) {
			fields = g_strsplit(token[i], " ", -1);
			sscanf(fields[1], "%lf", &center.x);
			sscanf(fields[2], "%lf", &center.y);
			world_cs = siril_world_cs_new_from_a_d(center.x, center.y);

			g_strfreev(fields);
		}
		i++;
	}
	if (world_cs && realname) {
		add_object_in_catalogue(realname, world_cs);
		siril_world_cs_unref(world_cs);
	}

	g_strfreev(token);
	return (world_cs && realname);
}

void on_search_objects_entry_activate(GtkEntry *entry, gpointer user_data) {
	if (!has_wcs()) return;
	gchar *result = search_in_catalogs(gtk_entry_get_text(GTK_ENTRY(entry)));
	if (result) {
		if (parse_buffer(result)) {
			GtkToggleToolButton *button = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("annotate_button"));
			if (!gtk_toggle_tool_button_get_active(button)) {
				gtk_toggle_tool_button_set_active(button, TRUE);
			}
			redraw(com.cvport, REMAP_NONE);
			gtk_entry_set_text(GTK_ENTRY(entry), "");
			gtk_widget_hide(lookup_widget("search_objects"));
		}
	}
}
