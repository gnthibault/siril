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

#include <math.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/undo.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "gui/dialogs.h"
#include "io/single_image.h"
#include "algos/Def_Wavelet.h"

#include "wavelets.h"

static fits wavelets_gfit_backup;

static void reset_scale_w() {
	static GtkRange *range_w[6] = { NULL, NULL, NULL, NULL, NULL, NULL };
	int i;

	if (range_w[0] == NULL) {
		range_w[0] = GTK_RANGE(lookup_widget("scale_w0"));
		range_w[1] = GTK_RANGE(lookup_widget("scale_w1"));
		range_w[2] = GTK_RANGE(lookup_widget("scale_w2"));
		range_w[3] = GTK_RANGE(lookup_widget("scale_w3"));
		range_w[4] = GTK_RANGE(lookup_widget("scale_w4"));
		range_w[5] = GTK_RANGE(lookup_widget("scale_w5"));
	}

	for (i = 0; i < 6; i++) {
		gtk_range_set_value(range_w[i], 1.f);
	}
}

static void update_wavelets() {
	float scale[6];
	static GtkRange *range_w[6] = { NULL, NULL, NULL, NULL, NULL, NULL };
	int i;
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];
	const char *tmpdir;

	tmpdir = g_get_tmp_dir();

	if (range_w[0] == NULL) {
		range_w[0] = GTK_RANGE(lookup_widget("scale_w0"));
		range_w[1] = GTK_RANGE(lookup_widget("scale_w1"));
		range_w[2] = GTK_RANGE(lookup_widget("scale_w2"));
		range_w[3] = GTK_RANGE(lookup_widget("scale_w3"));
		range_w[4] = GTK_RANGE(lookup_widget("scale_w4"));
		range_w[5] = GTK_RANGE(lookup_widget("scale_w5"));
	}

	for (i = 0; i < 6; i++)
		scale[i] = (float) gtk_range_get_value(range_w[i]);

	set_cursor_waiting(TRUE);

	for (i = 0; i < gfit.naxes[2]; i++) {
		dir[i] = g_build_filename(tmpdir, File_Name_Transform[i], NULL);
		wavelet_reconstruct_file(dir[i], scale, gfit.pdata[i]);
		g_free(dir[i]);
	}

	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

static void wavelets_startup() {
	copyfits(&gfit, &wavelets_gfit_backup, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
}

void on_menuitem_wavelets_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded()) {
		siril_open_dialog("wavelets_dialog");
	}
}

void on_wavelets_dialog_show(GtkWidget *widget, gpointer user_data) {
	reset_scale_w();
	wavelets_startup();
}

gboolean on_scale_w0_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	update_wavelets();
	return FALSE;
}

gboolean on_scale_w1_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	update_wavelets();
	return FALSE;
}

gboolean on_scale_w2_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	update_wavelets();
	return FALSE;
}

gboolean on_scale_w3_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	update_wavelets();
	return FALSE;
}

gboolean on_scale_w4_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	update_wavelets();
	return FALSE;
}

gboolean on_scale_w5_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	update_wavelets();
	return FALSE;
}


gboolean on_scale_w0_key_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	update_wavelets();
	return FALSE;
}

gboolean on_scale_w1_key_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	update_wavelets();
	return FALSE;
}

gboolean on_scale_w2_key_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	update_wavelets();
	return FALSE;
}

gboolean on_scale_w3_key_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	update_wavelets();
	return FALSE;
}

gboolean on_scale_w4_key_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	update_wavelets();
	return FALSE;
}

gboolean on_scale_w5_key_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	update_wavelets();
	return FALSE;
}

void on_wavelets_dialog_hide(GtkWidget *widget, gpointer user_data) {
	gtk_widget_set_sensitive(lookup_widget("grid_w"), FALSE);
	gtk_widget_set_sensitive(lookup_widget("button_reset_w"), FALSE);
	clearfits(&wavelets_gfit_backup);
}


void on_button_reset_w_clicked(GtkButton *button, gpointer user_data) {
	reset_scale_w();
	update_wavelets();
}

void apply_wavelets_cancel() {
	if (gtk_widget_get_sensitive(lookup_widget("grid_w")) == TRUE) {
		reset_scale_w();
		update_wavelets();
	}
}

void on_button_ok_w_clicked(GtkButton *button, gpointer user_data) {
	if (gtk_widget_get_sensitive(lookup_widget("grid_w")) == TRUE) {
		update_wavelets();
		undo_save_state(&wavelets_gfit_backup, "Processing: Wavelets Transformation");
	}
	siril_close_dialog("wavelets_dialog");
}

void on_button_cancel_w_clicked(GtkButton *button, gpointer user_data) {
	apply_wavelets_cancel();
	siril_close_dialog("wavelets_dialog");
}

void on_button_compute_w_clicked(GtkButton *button, gpointer user_data) {
	int Type_Transform, Nbr_Plan, maxplan, mins, i;
	int nb_chan = gfit.naxes[2];
	float *Imag;
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];
	const char *tmpdir;

	tmpdir = g_get_tmp_dir();

	Nbr_Plan = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_plans_w")));
	Type_Transform = gtk_combo_box_get_active(
			GTK_COMBO_BOX(lookup_widget("combobox_type_w"))) + 1;

	mins = min(gfit.rx, gfit.ry);
	maxplan = log(mins) / log(2) - 2;

	if (Nbr_Plan > maxplan) {
		char *msg = siril_log_message(
				_("Wavelet: maximum number of plans for this image size is %d\n"),
				maxplan);
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Warning"), msg);
		Nbr_Plan = maxplan;
		gtk_spin_button_set_value(
				GTK_SPIN_BUTTON(lookup_widget("spinbutton_plans_w")), Nbr_Plan);
	}

	if (Type_Transform != TO_PAVE_LINEAR && Type_Transform != TO_PAVE_BSPLINE) {
		char *msg = siril_log_message(_("Wavelet: type must be %d or %d\n"),
		TO_PAVE_LINEAR, TO_PAVE_BSPLINE);
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Warning"), msg);
	}

	set_cursor_waiting(TRUE);

	Imag = (float *) malloc(gfit.rx * gfit.ry * sizeof(float));

	for (i = 0; i < nb_chan; i++) {
		dir[i] = malloc(strlen(tmpdir) + strlen(File_Name_Transform[i]) + 2);
		strcpy(dir[i], tmpdir);
		strcat(dir[i], G_DIR_SEPARATOR_S);
		strcat(dir[i], File_Name_Transform[i]);
		wavelet_transform_file(Imag, gfit.ry, gfit.rx, dir[i], Type_Transform, Nbr_Plan,
				gfit.pdata[i]);
		free(dir[i]);
	}

	free(Imag);
	Imag = NULL;
	gtk_widget_set_sensitive(lookup_widget("grid_w"), TRUE);
	gtk_widget_set_sensitive(lookup_widget("button_reset_w"), TRUE);
	set_cursor_waiting(FALSE);
	return;
}
