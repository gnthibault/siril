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

#include <gtk/gtk.h>
#include "core/siril_actions.h"

static GActionEntry win_entries[] = {
	{ "snapshot", snapshot_action_activate },
	{ "close", close_action_activate },
	{ "undo", undo_action_activate },
	{ "redo", redo_action_activate },
	{ "scripts", scripts_action_activate },
	{ "updates", updates_action_activate },
	{ "full-screen", full_screen_activated},
	{ "hide-show-toolbar", toolbar_activate },
	{ "shortcuts", keyboard_shortcuts_activated},
	{ "cwd", cwd_action_activate },

	{ "conversion", tab_conversion_activate },
	{ "sequence", tab_sequence_activate },
	{ "registration", tab_registration_activate },
	{ "prepro", tab_prepro_activate },
	{ "plot", tab_plot_activate },
	{ "stacking", tab_stacking_activate },
	{ "logs", tab_logs_activate }
};

static GActionEntry image_entries[] = {
	{ "zoom-out", zoom_out_activate },
	{ "zoom-in", zoom_in_activate },
	{ "zoom-fit", zoom_fit_activate, NULL, "true", change_zoom_fit_state },
	{ "statistics", statistics_activate },
	{ "evaluate-noise", noise_activate },
	{ "astrometry", astrometry_activate },
	{ "image-information", astrometry_activate },
	{ "dyn-psf", dyn_psf_activate },
	{ "search-object", search_object_activate }
};

static void _siril_window_enable_action_group(GActionMap *map,
		const gchar **group, gboolean enable) {
	GAction *action;
	const gchar **it = group;

	for (it = group; *it != NULL; it++) {
		action = g_action_map_lookup_action(map, *it);
		if (G_LIKELY(action))
			g_simple_action_set_enabled(G_SIMPLE_ACTION(action), enable);
		else
			g_warning("Action not found in action group: %s", *it);
	}
}

void siril_window_enable_image_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *image_actions[] = {
		"zoom-out",
		"zoom-in",
		"zoom-fit",
		"statistics",
		"evaluate-noise",
		"astrometry",
		"image-information",
	    "dyn-psf",
        "search-object",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), image_actions, enable);
}

void siril_window_map_actions(GtkApplicationWindow *window) {
	g_action_map_add_action_entries(G_ACTION_MAP(window), win_entries, G_N_ELEMENTS(win_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), image_entries, G_N_ELEMENTS(image_entries), window);
}
