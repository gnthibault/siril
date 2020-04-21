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
#include "core/initfile.h"
#include "core/siril_language.h"
#include "algos/photometry.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/script_menu.h"
#include "gui/dialogs.h"
#include "gui/PSF_list.h"

#include "preferences.h"

#ifndef W_OK
#define W_OK 2
#endif

/*
 * Set swap dir to default
 */

static void reset_swapdir() {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_swap"));
	const gchar *dir;

	dir = g_get_tmp_dir();

	if (g_strcmp0(dir, com.pref.swap_dir)) {
		g_free(com.pref.swap_dir);
		com.pref.swap_dir = g_strdup(dir);
		gtk_file_chooser_set_filename(swap_dir, dir);
		writeinitfile();
	}
}

void update_libraw_and_debayer_interface() {
	/**********COLOR ADJUSTEMENT**************/
	com.pref.raw_set.bright = gtk_spin_button_get_value(	GTK_SPIN_BUTTON(lookup_widget("Brightness_spinbutton")));
	com.pref.raw_set.mul[0] = gtk_spin_button_get_value(	GTK_SPIN_BUTTON(lookup_widget("Red_spinbutton")));
	com.pref.raw_set.mul[2] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("Blue_spinbutton")));

	com.pref.raw_set.auto_mul = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_multipliers")));
	com.pref.raw_set.user_black = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_blackpoint")));

	/**************WHITE BALANCE**************/
	com.pref.raw_set.use_camera_wb = (int) gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_cam")));
	com.pref.raw_set.use_auto_wb = (int) gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_auto")));

	/********MATRIX INTERPOLATION**************/
	com.pref.raw_set.user_qual = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_dcraw_inter")));

	/********GAMMA CORRECTION**************/
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm0")))) {
		/* Linear Gamma Curve */
		com.pref.raw_set.gamm[0] = 1.0;
		com.pref.raw_set.gamm[1] = 1.0;
	} else if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm1")))) {
		/* BT.709 Gamma curve */
		com.pref.raw_set.gamm[0] = 2.222;
		com.pref.raw_set.gamm[1] = 4.5;
	} else {
		/* sRGB Gamma curve */
		com.pref.raw_set.gamm[0] = 2.40;
		com.pref.raw_set.gamm[1] = 12.92;
	}
	/* We write in config file */
	/*************SER**********************/
	com.pref.debayer.use_bayer_header = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_SER_use_header")));
	com.pref.debayer.compatibility = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_debayer_compatibility")));
	com.pref.debayer.xbayeroff = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(lookup_widget("xbayeroff_spin")));
	com.pref.debayer.ybayeroff = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(lookup_widget("ybayeroff_spin")));
	writeinitfile();
}

void on_button_reset_photometry_clicked(GtkButton *button, gpointer user_data) {

	initialize_photometric_param();
	double tmp = gfit.cvf;
	gfit.cvf = 0.0;
	set_GUI_photometry();
	gfit.cvf = tmp;
}

void update_photometry_interface() {
	com.pref.phot_set.gain = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinGain")));
	com.pref.phot_set.inner = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinInner")));
	com.pref.phot_set.outer = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinOuter")));
	com.pref.phot_set.minval = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinMinPhot")));
	com.pref.phot_set.maxval = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinMaxPhot")));
	writeinitfile();
}

void set_GUI_LIBRAW() {
	/**********COLOR ADJUSTEMENT**************/
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("Brightness_spinbutton")), com.pref.raw_set.bright);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("Red_spinbutton")), com.pref.raw_set.mul[0]);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("Blue_spinbutton")), com.pref.raw_set.mul[2]);

	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_multipliers")), com.pref.raw_set.auto_mul);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_blackpoint")), com.pref.raw_set.user_black);

	/**************WHITE BALANCE**************/
	if (com.pref.raw_set.use_camera_wb) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_cam")), com.pref.raw_set.use_camera_wb);
	}

	if (com.pref.raw_set.use_auto_wb) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_auto")), com.pref.raw_set.use_auto_wb);
	}

	/********MATRIX INTERPOLATION**************/
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_dcraw_inter")), com.pref.raw_set.user_qual);

	/********GAMMA CORRECTION**************/
	if (com.pref.raw_set.gamm[0] == 1.0 && com.pref.raw_set.gamm[1] == 1.0)
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm0")), TRUE);
	else if (com.pref.raw_set.gamm[0] == 2.222 && com.pref.raw_set.gamm[1] == 4.5)
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm1")), TRUE);
	else
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm2")), TRUE);

	/********** DEBAYER ******************/
	GtkComboBox *pattern = GTK_COMBO_BOX(lookup_widget("comboBayer_pattern"));
	GtkComboBox *inter = GTK_COMBO_BOX(lookup_widget("comboBayer_inter"));
	GtkToggleButton *compat = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_debayer_compatibility"));
	GtkToggleButton *use_header = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_SER_use_header"));
	GtkToggleButton *demosaicingButton = GTK_TOGGLE_BUTTON(lookup_widget("demosaicingButton"));
	GtkSpinButton *xbayer_spin = GTK_SPIN_BUTTON(lookup_widget("xbayeroff_spin"));
	GtkSpinButton *ybayer_spin = GTK_SPIN_BUTTON(lookup_widget("ybayeroff_spin"));
	gtk_combo_box_set_active(pattern, com.pref.debayer.bayer_pattern);
	gtk_combo_box_set_active(inter, com.pref.debayer.bayer_inter);
	gtk_toggle_button_set_active(compat, com.pref.debayer.compatibility);
	gtk_toggle_button_set_active(use_header, com.pref.debayer.use_bayer_header);
	gtk_toggle_button_set_active(demosaicingButton,	com.pref.debayer.open_debayer);
	gtk_spin_button_set_value(xbayer_spin, com.pref.debayer.xbayeroff);
	gtk_spin_button_set_value(ybayer_spin, com.pref.debayer.ybayeroff);
}

void on_checkbutton_cam_toggled(GtkButton *button, gpointer user_data) {
	GtkToggleButton *auto_button = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_auto"));
	GtkToggleButton *cam_button = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_cam"));

	if (gtk_toggle_button_get_active(auto_button)) {
		g_signal_handlers_block_by_func(auto_button, on_checkbutton_auto_toggled, NULL);
		gtk_toggle_button_set_active(auto_button, FALSE);
		g_signal_handlers_unblock_by_func(auto_button, on_checkbutton_auto_toggled, NULL);
		gtk_toggle_button_set_active(cam_button, TRUE);
	}
}

void on_checkbutton_auto_toggled(GtkButton *button, gpointer user_data) {
	GtkToggleButton *auto_button = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_auto"));
	GtkToggleButton *cam_button = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_cam"));

	if (gtk_toggle_button_get_active(cam_button)) {
		g_signal_handlers_block_by_func(cam_button, on_checkbutton_cam_toggled, NULL);
		gtk_toggle_button_set_active(cam_button, FALSE);
		g_signal_handlers_unblock_by_func(cam_button, on_checkbutton_cam_toggled, NULL);
		gtk_toggle_button_set_active(auto_button, TRUE);
	}
}

void on_checkbutton_auto_evaluate_toggled(GtkToggleButton *button,
		gpointer user_data) {
	GtkWidget *entry = (GtkWidget *)user_data;
	gtk_widget_set_sensitive(entry, !gtk_toggle_button_get_active(button));
}

void on_checkbutton_multipliers_toggled(GtkToggleButton *button,
		gpointer user_data) {
	gboolean active = gtk_toggle_button_get_active(button);

	gtk_widget_set_sensitive(lookup_widget("hbox8"), !active);
	gtk_widget_set_sensitive(lookup_widget("hbox11"), !active);
	if (active) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("Red_spinbutton")), 1.0);
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("Blue_spinbutton")), 1.0);
	}
}

void set_GUI_photometry() {
	if (gfit.cvf > 0.0) {
		com.pref.phot_set.gain = gfit.cvf;
	}
	if (com.pref.phot_set.gain > 0.0) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinGain")), com.pref.phot_set.gain);
	}
	if (com.pref.phot_set.inner > 0.0) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinInner")), com.pref.phot_set.inner);
	}
	if (com.pref.phot_set.outer > 0.0) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinOuter")), com.pref.phot_set.outer);
	}
	if (com.pref.phot_set.minval >= 0.0) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinMinPhot")), (gdouble) com.pref.phot_set.minval);
	}
	if (com.pref.phot_set.maxval >= 0.0) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinMaxPhot")), (gdouble) com.pref.phot_set.maxval);
	}
}

void initialize_path_directory() {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_swap"));
	if (com.pref.swap_dir && com.pref.swap_dir[0] != '\0') {
		gtk_file_chooser_set_filename (swap_dir, com.pref.swap_dir);
	} else {
		gtk_file_chooser_set_filename (swap_dir, g_get_tmp_dir());
	}
}

void set_libraw_settings_menu_available(gboolean activate) {
	if (!com.script) {
		gtk_widget_set_visible(lookup_widget("box_stack_p1"), activate);
	}
}

void on_comboBayer_pattern_changed(GtkComboBox* box, gpointer user_data) {
	com.pref.debayer.bayer_pattern = gtk_combo_box_get_active(box);
}

void on_comboBayer_inter_changed(GtkComboBox* box, gpointer user_data) {
	com.pref.debayer.bayer_inter = gtk_combo_box_get_active(box);
}

void on_button_reset_swap_clicked(GtkButton *button, gpointer user_data) {
	reset_swapdir();
}

void on_spinbutton_mem_ratio_value_changed(GtkSpinButton *button, gpointer user_data) {
	gdouble mem = gtk_spin_button_get_value(button);
	com.pref.stack.memory_ratio = mem;
	writeinitfile();
}

void on_spinbutton_mem_amount_value_changed(GtkSpinButton *button, gpointer user_data) {
	gdouble mem = gtk_spin_button_get_value(button);
	com.pref.stack.memory_amount = mem;
	writeinitfile();
}

void on_mem_radio_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkToggleButton *ratio = GTK_TOGGLE_BUTTON(lookup_widget("memfreeratio_radio")),
			*amount = GTK_TOGGLE_BUTTON(lookup_widget("memfixed_radio")),
			*unlimited = GTK_TOGGLE_BUTTON(lookup_widget("memunlimited_radio"));
	GtkWidget *ratio_spin = lookup_widget("spinbutton_mem_ratio"),
		  *amount_spin = lookup_widget("spinbutton_mem_amount");
	if (!gtk_toggle_button_get_active(togglebutton)) return;

	if (togglebutton == ratio) {
		com.pref.stack.mem_mode = 0;
		gtk_widget_set_sensitive(ratio_spin, TRUE);
		gtk_widget_set_sensitive(amount_spin, FALSE);
	} else if (togglebutton == amount) {
		com.pref.stack.mem_mode = 1;
		gtk_widget_set_sensitive(ratio_spin, FALSE);
		gtk_widget_set_sensitive(amount_spin, TRUE);
	} else if (togglebutton == unlimited) {
		com.pref.stack.mem_mode = 2;
		gtk_widget_set_sensitive(ratio_spin, FALSE);
		gtk_widget_set_sensitive(amount_spin, FALSE);
	}
}

void on_combobox_ext_changed(GtkComboBox *box, gpointer user_data) {
	if (com.pref.ext)
		free(com.pref.ext);

	com.pref.ext = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(box));
	writeinitfile();
	initialize_FITS_name_entries();
}

void on_combobox_type_changed(GtkComboBox *box, gpointer user_data) {
	com.pref.force_to_16bit = gtk_combo_box_get_active(box) == 0;
	writeinitfile();
}

void on_filechooser_swap_file_set(GtkFileChooserButton *fileChooser, gpointer user_data) {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(fileChooser);
	gchar *dir;

	dir = gtk_file_chooser_get_filename (swap_dir);

	if (g_access(dir, W_OK)) {
		gchar *msg = siril_log_color_message(_("You don't have permission to write in this directory: %s\n"), "red", dir);
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error"), msg);
		gtk_file_chooser_set_filename(swap_dir, com.pref.swap_dir);
		return;
	}

	g_free(com.pref.swap_dir);
	com.pref.swap_dir = dir;
	writeinitfile();
}

void on_rememberWindowsCheck_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	com.pref.remember_windows = gtk_toggle_button_get_active(togglebutton);
}

void on_show_preview_button_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkWidget *label = lookup_widget("thumbnails_label_size");
	GtkWidget *box = lookup_widget("thumbnails_box_size");

	com.pref.show_thumbnails = gtk_toggle_button_get_active(togglebutton);
	gtk_widget_set_sensitive(label, com.pref.show_thumbnails);
	gtk_widget_set_sensitive(box, com.pref.show_thumbnails);
}

void on_thumbnails_box_size_changed(GtkComboBoxText *box, gpointer user_data) {
	com.pref.thumbnail_size = (gtk_combo_box_get_active(GTK_COMBO_BOX(box)) == 0) ? 128 : 256;
}

void on_apply_settings_button_clicked(GtkButton *button, gpointer user_data) {
	update_libraw_and_debayer_interface();
	update_photometry_interface();
	update_language();
	fill_script_paths_list();
	refresh_stars_list(com.stars);
	save_main_window_state();
	siril_close_dialog("settings_window");
}

