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

#include "core/siril.h"
#include "core/proto.h"
#include "core/command.h"
#include "core/undo.h"
#include "core/siril_update.h"
#include "core/siril_cmd_help.h"
#include "algos/annotate.h"
#include "algos/colors.h"
#include "algos/noise.h"
#include "algos/geometry.h"
#include "algos/siril_wcs.h"
#include "algos/plateSolver.h"
#include "compositing/compositing.h"
#include "gui/about_dialog.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/histogram.h"
#include "gui/open_dialog.h"
#include "gui/message_dialog.h"
#include "gui/PSF_list.h"
#include "gui/save_dialog.h"
#include "gui/sequence_list.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/script_menu.h"
#include "gui/image_interactions.h"
#include "gui/image_display.h"
#include "gui/photometric_cc.h"

#include "siril_actions.h"

void open_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	header_open_button_clicked();
}

void cwd_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	cwd_btton_clicked();
}

void save_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_header_save_button_clicked();
}

void save_as_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_header_save_as_button_clicked();
}

void snapshot_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	on_header_snapshot_button_clicked();
}

void undo_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	set_cursor_waiting(TRUE);
	undo_display_data(UNDO);
	set_cursor_waiting(FALSE);
}

void redo_action_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	set_cursor_waiting(TRUE);
	undo_display_data(REDO);
	set_cursor_waiting(FALSE);
}

void quit_action_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	siril_quit();
}

void about_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_show_about_dialog();
}

void preferences_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("settings_window");
}

void close_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	process_close(0);
}

void scripts_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_get_on_script_pages();
}

void updates_action_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_check_updates(TRUE);
}

static gboolean is_extended = FALSE;

void full_screen_activated(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkWindow *window;
	GtkWidget *toolbarbox = lookup_widget("toolbarbox");
	GtkWidget *control_center_box = lookup_widget("control_center_box");
	GtkButton *button = GTK_BUTTON(lookup_widget("button_paned"));
	gboolean is_control_box_visible;
	gboolean is_fullscreen;

	window = GTK_WINDOW(GTK_APPLICATION_WINDOW(user_data));

	GdkWindow *gdk_window = gtk_widget_get_window(GTK_WIDGET(window));
	is_fullscreen = gdk_window_get_state(gdk_window) & GDK_WINDOW_STATE_FULLSCREEN;
	is_control_box_visible = gtk_widget_get_visible(control_center_box);

	if (is_fullscreen) {
		gtk_window_unfullscreen(window);
		if (is_extended)
			gtk_button_clicked(button);
	} else {
		gtk_window_fullscreen(window);
		if (is_control_box_visible) {
			gtk_button_clicked(button);
		}
		is_extended = is_control_box_visible;
	}
	gtk_widget_set_visible(toolbarbox, is_fullscreen);
}

void keyboard_shortcuts_activated(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkWindow *window;

	window = GTK_WINDOW(GTK_APPLICATION_WINDOW(user_data));

	siril_cmd_help_keyboard_shortcuts(window);
}

void tab_conversion_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(FILE_CONVERSION);
}

void tab_sequence_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(IMAGE_SEQ);
}

void tab_prepro_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(PRE_PROC);
}

void tab_registration_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(REGISTRATION);
}

void tab_plot_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(PLOT);
}

void tab_stacking_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(STACKING);
}

void tab_logs_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(OUTPUT_LOGS);
}

void toolbar_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkWidget *w = lookup_widget("toolbarbox");
	gtk_widget_set_visible(w, !gtk_widget_get_visible(w));
}

void change_zoom_fit_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	if (g_variant_get_boolean(state)) {
		com.zoom_value = ZOOM_FIT;
		reset_display_offset();
		redraw(com.cvport, REMAP_NONE);
	} else {
		com.zoom_value = get_zoom_val();
	}
	g_simple_action_set_state(action, state);
}

void zoom_fit_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void zoom_in_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	point center = get_center_of_vport();
	update_zoom(center.x, center.y, ZOOM_IN);
}

void zoom_out_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	point center = get_center_of_vport();
	update_zoom(center.x, center.y, ZOOM_OUT);
}

void zoom_one_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	update_zoom_fit_button();
	com.zoom_value = ZOOM_NONE;
	reset_display_offset();
	redraw(com.cvport, REMAP_NONE);
}

void negative_view_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	set_cursor_waiting(TRUE);
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	g_simple_action_set_state(action, state);
}

void negative_view_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void photometry_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	mouse_status = g_variant_get_boolean(state) ? MOUSE_ACTION_PHOTOMETRY : MOUSE_ACTION_SELECT_REG_AREA;
	g_simple_action_set_state(action, state);
	free(com.qphot);
	com.qphot = NULL;
	redraw(com.cvport, REMAP_NONE);
}

void photometry_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void color_map_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	set_cursor_waiting(TRUE);
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	g_simple_action_set_state(action, state);
}

void color_map_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void astrometry_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	open_astrometry_dialog();
}

void dyn_psf_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("stars_list_window");
}

void pick_star_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	pick_a_star();
}

void psf_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	psf_star *result = NULL;
	int layer = match_drawing_area_widget(com.vport[com.cvport], FALSE);

	if (layer == -1)
		return;
	if (!(com.selection.h && com.selection.w))
		return;
	result = psf_get_minimisation(&gfit, layer, &com.selection, TRUE, TRUE, TRUE);
	if (!result)
		return;

	popup_psf_result(result, &com.selection);
	free_psf(result);
}

void seq_psf_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	process_seq_psf(0);
}

void crop_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_crop();
}

void seq_crop_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("crop_dialog");
}

void search_object_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (has_wcs(&gfit))
		siril_open_dialog("search_objects");
}

void annotate_object_state(GSimpleAction *action, GVariant *state, gpointer user_data) {
	if (g_variant_get_boolean(state)) {
		if (has_wcs(&gfit)) {
			com.found_object = find_objects(&gfit);
		}
	} else {
		g_slist_free_full(com.found_object, (GDestroyNotify) free_object);
		com.found_object = NULL;
	}
	g_simple_action_set_state(action, state);
	redraw(com.cvport, REMAP_NONE);
}

void annotate_object_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GVariant *state;

	state = g_action_get_state(G_ACTION(action));
	g_action_change_state(G_ACTION(action), g_variant_new_boolean(!g_variant_get_boolean(state)));
	g_variant_unref(state);
}

void seq_list_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	if (gtk_widget_get_visible(lookup_widget("seqlist_dialog"))) {
		siril_close_dialog("seqlist_dialog");
	} else {
		gboolean confirm = TRUE;
		if (com.seq.current == RESULT_IMAGE) {
			confirm = siril_confirm_dialog(_("Save your changes before loading a frame of the sequence."),
					_("The image currently displayed is the result of the previous stack. "
							"If you load an image from the sequence, you might lose the entire process you performed on the image, "
							"but not the image itself. You need to save your data before doing this."),
					_("Load another image"));
		}
		if (confirm) {
			update_seqlist();
			siril_open_dialog("seqlist_dialog");
		}
	}
}

void statistics_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	set_cursor_waiting(TRUE);
	computeStat();
	siril_open_dialog("StatWindow");
	set_cursor_waiting(FALSE);
}

void noise_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	evaluate_noise_in_image();
}

void image_information_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("file_information");
}

void image_fits_header_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	show_FITS_header(&gfit);
}

/******* processing menu **************/

void remove_green_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("SCNR_dialog");
}

void saturation_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("satu_dialog");
}

void color_calib_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	initialize_calibration_interface();
	siril_open_dialog("color_calibration");
}

void pcc_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	initialize_photometric_cc_dialog();
	siril_open_dialog("ImagePlateSolver_Dial");
}

void split_channel_activate(GSimpleAction *action, GVariant *parameter,gpointer user_data) {
	siril_open_dialog("extract_channel_dialog");
}

void negative_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	negative_processing();
}

void histo_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	toggle_histogram_window_visibility();
}

void fix_banding_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("canon_fixbanding_dialog");
}

void cosmetic_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("cosmetic_dialog");
}

void background_extr_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("background_extraction_dialog");
}

void asinh_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("asinh_dialog");
}

void deconvolution_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("deconvolution_dialog");
}

void resample_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("resample_dialog");
}

void rotation_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("rotation_dialog");
}

void rotation90_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_rotate90();
}

void rotation270_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_rotate270();
}

void mirrorx_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	mirrorx_gui(&gfit);
}

void mirrory_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	mirrory_gui(&gfit);
}

void wavelets_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("wavelets_dialog");
}

void split_wavelets_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("extract_wavelets_layers_dialog");
}

void medianfilter_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("Median_dialog");
}

void rgradient_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("rgradient_dialog");
}

void clahe_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("CLAHE_dialog");
}

void linearmatch_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("linearmatch_dialog");
}

void fft_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	GtkFileChooserButton *magbutton, *phasebutton;

	magbutton = GTK_FILE_CHOOSER_BUTTON(lookup_widget("filechooser_mag"));
	phasebutton = GTK_FILE_CHOOSER_BUTTON(lookup_widget("filechooser_phase"));
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(magbutton), com.wd);
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(phasebutton), com.wd);
	siril_open_dialog("dialog_FFT");
}

void rgb_compositing_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	open_compositing_window();
}

void split_cfa_activate(GSimpleAction *action, GVariant *parameter, gpointer user_data) {
	siril_open_dialog("split_cfa_dialog");
}
