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
#ifndef SRC_CORE_SIRIL_ACTIONS_H_
#define SRC_CORE_SIRIL_ACTIONS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <glib.h>

void open_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void cwd_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void save_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void save_as_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void snapshot_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void undo_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void redo_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void quit_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void about_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void preferences_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void close_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void scripts_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) ;
void updates_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void full_screen_activated(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void keyboard_shortcuts_activated(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void tab_conversion_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void tab_sequence_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void tab_prepro_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void tab_registration_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void tab_plot_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void tab_stacking_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) ;
void tab_logs_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void zoom_in_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void zoom_out_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void zoom_one_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void zoom_fit_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void negative_view_state(GSimpleAction *action, GVariant *state, gpointer user_data);
void negative_view_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void color_map_state(GSimpleAction *action, GVariant *state, gpointer user_data);
void color_map_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void toolbar_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void change_zoom_fit_state(GSimpleAction *action, GVariant *state, gpointer user_data);
void astrometry_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void dyn_psf_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void search_object_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void annotate_object_state(GSimpleAction *action, GVariant *state, gpointer user_data);
void annotate_object_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void seq_list_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void statistics_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void noise_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void image_information_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void image_fits_header_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);

void remove_green_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void saturation_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void color_calib_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void pcc_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void split_channel_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void negative_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void histo_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void fix_banding_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void cosmetic_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void background_extr_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void asinh_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void deconvolution_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void resample_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void rotation_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void rotation90_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void rotation270_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void mirrorx_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void mirrory_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void wavelets_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void split_wavelets_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void medianfilter_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void rgradient_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void clahe_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void linearmatch_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void fft_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void rgb_compositing_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);
void split_cfa_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data);


#endif /* SRC_CORE_SIRIL_ACTIONS_H_ */
