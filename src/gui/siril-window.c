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
		{ "bit-depth", NULL },
		{ "zoom-out", zoom_out_activate },
		{ "zoom-in", zoom_in_activate },
		{ "zoom-fit", zoom_fit_activate, NULL, "true", change_zoom_fit_state },
		{ "zoom-one", zoom_one_activate },
		{ "negative-view", negative_view_activate, NULL, "false", negative_view_state },
		{ "color-map", color_map_activate, NULL, "false", color_map_state },
		{ "snapshot", snapshot_action_activate },
		{ "fits-header", image_fits_header_activate },
		{ "statistics", statistics_activate },
		{ "evaluate-noise", noise_activate },
		{ "astrometry", astrometry_activate },
		{ "photometry", photometry_activate, NULL, "false", photometry_state },
		{ "image-information", image_information_activate },
		{ "dyn-psf", dyn_psf_activate },
		{ "annotate-object", annotate_object_activate, NULL, "false", annotate_object_state },
		{ "search-object", search_object_activate },
		{ "seq-list", seq_list_activate }
};

static GActionEntry selection_entries[] = {
		{ "pickstar", pick_star_activate },
		{ "psf", psf_activate },
		{ "crop", crop_activate }
};

static GActionEntry selection_sequence_entries[] = {
		{ "seq-psf", seq_psf_activate },
		{ "seq-crop", seq_crop_activate }
};

static GActionEntry rgb_processing_entries[] = {
		{ "remove-green-processing", remove_green_activate },
		{ "saturation-processing", saturation_activate },
		{ "color-calib-processing", color_calib_activate },
		{ "pcc-processing", pcc_activate },
		{ "split-channel-processing", split_channel_activate }
};

static GActionEntry any_processing_entries[] = {
		{ "negative-processing", negative_activate },
		{ "histo-processing", histo_activate },
		{ "fix-banding-processing", fix_banding_activate },
		{ "cosmetic-processing", cosmetic_activate },
		{ "background-extr-processing", background_extr_activate }
};

static GActionEntry any_mono_processing_entries[] = {
		{ "split-cfa-processing", split_cfa_activate }
};

static GActionEntry single_processing_entries[] = {
		{ "asinh-processing", asinh_activate },
		{ "deconvolution-processing", deconvolution_activate },
		{ "resample-processing", resample_activate },
		{ "rotation-processing", rotation_activate },
		{ "rotation90-processing", rotation90_activate },
		{ "rotation270-processing", rotation270_activate },
		{ "mirrorx-processing", mirrorx_activate },
		{ "mirrory-processing", mirrory_activate },
		{ "wavelets-processing", wavelets_activate },
		{ "split-wavelets-processing", split_wavelets_activate },
		{ "medianfilter-processing", medianfilter_activate },
		{ "rgradient-processing", rgradient_activate },
		{ "clahe-processing", clahe_activate },
		{ "linearmatch-processing", linearmatch_activate }
};

static GActionEntry none_processing_entries[] = {
		{ "fft-processing", fft_activate },
		{ "rgb-compositing-processing", rgb_compositing_activate }
};

static void _siril_window_enable_action_group(GActionMap *map,
		const gchar **group, gboolean enable) {
	GAction *action;

	for (const gchar **it = group; *it != NULL; it++) {
		action = g_action_map_lookup_action(map, *it);
		if (G_LIKELY(action))
			g_simple_action_set_enabled(G_SIMPLE_ACTION(action), enable);
		else
			g_warning("Action not found in action group: %s", *it);
	}
}

void siril_window_enable_image_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *image_actions[] = {
		"bit-depth",
		"zoom-out",
		"zoom-in",
		"zoom-fit",
		"zoom-one",
		"negative-view",
		"color-map",
		"snapshot",
		"statistics",
		"evaluate-noise",
		"astrometry",
		"photometry",
		"image-information",
	    "dyn-psf",
        "search-object",
		"seq-list",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), image_actions, enable);
}

void siril_window_enable_rgb_proc_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *rgb_processing_actions[] = {
		"remove-green-processing",
		"saturation-processing",
		"color-calib-processing",
		"pcc-processing",
		"split-channel-processing",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), rgb_processing_actions, enable);
}

void siril_window_enable_any_proc_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *any_processing_actions[] = {
		"negative-processing",
		"histo-processing",
		"fix-banding-processing",
		"cosmetic-processing",
		"background-extr-processing",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), any_processing_actions, enable);
}

void siril_window_enable_any_mono_proc_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *any_mono_processing_actions[] = {
		"split-cfa-processing",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), any_mono_processing_actions, enable);
}

void siril_window_enable_single_proc_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *single_processing_actions[] = {
		"asinh-processing",
		"deconvolution-processing",
		"resample-processing",
		"rotation-processing",
		"rotation90-processing",
		"rotation270-processing",
		"mirrorx-processing",
		"mirrory-processing",
		"wavelets-processing",
		"split-wavelets-processing",
		"medianfilter-processing",
		"rgradient-processing",
		"clahe-processing",
		"linearmatch-processing",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), single_processing_actions, enable);
}

void siril_window_enable_none_proc_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *none_processing_actions[] = {
		"fft-processing",
		"rgb-compositing-processing",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), none_processing_actions, enable);
}

void siril_window_enable_if_selection_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *selection_actions[] = {
		"pickstar",
		"psf",
		"crop",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), selection_actions, enable);
}

void siril_window_enable_if_selection_sequence_actions(GtkApplicationWindow *window, gboolean enable) {
	static const gchar *selection_sequebce_actions[] = {
		"seq-psf",
		"seq-crop",
		NULL,
	};
	_siril_window_enable_action_group(G_ACTION_MAP(window), selection_sequebce_actions, enable);
}

void siril_window_map_actions(GtkApplicationWindow *window) {
	g_action_map_add_action_entries(G_ACTION_MAP(window), win_entries, G_N_ELEMENTS(win_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), image_entries, G_N_ELEMENTS(image_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), rgb_processing_entries, G_N_ELEMENTS(rgb_processing_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), any_processing_entries, G_N_ELEMENTS(any_processing_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), any_mono_processing_entries, G_N_ELEMENTS(any_mono_processing_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), single_processing_entries, G_N_ELEMENTS(single_processing_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), none_processing_entries, G_N_ELEMENTS(none_processing_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), selection_entries, G_N_ELEMENTS(selection_entries), window);
	g_action_map_add_action_entries(G_ACTION_MAP(window), selection_sequence_entries, G_N_ELEMENTS(selection_sequence_entries), window);
}
