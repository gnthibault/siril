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

#include <math.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/undo.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "gui/dialogs.h"
#include "io/single_image.h"
#include "algos/Def_Wavelet.h"
#include "gui/preview_timer.h"

#include "wavelets.h"

static float wavelet_value[6];
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

static int update_wavelets() {
	int i;
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];
	const char *tmpdir;

	tmpdir = g_get_tmp_dir();

	set_cursor_waiting(TRUE);

	for (i = 0; i < gfit.naxes[2]; i++) {
		dir[i] = g_build_filename(tmpdir, File_Name_Transform[i], NULL);
		wavelet_reconstruct_file(dir[i], wavelet_value, gfit.pdata[i]);
		g_free(dir[i]);
	}
	return 0;
}

static void wavelets_startup() {
	int i;
	for (i = 0; i < 6; i++) {
		wavelet_value[i] = 1.f;
	}
	copyfits(&gfit, &wavelets_gfit_backup, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
}

/* This function computes wavelets with the number of Nbr_Plan and
 * extracts plan "Plan" in fit parameters */

int get_wavelet_layers(fits *fit, int Nbr_Plan, int Plan, int Type, int reqlayer) {
	int chan, start, end, retval = 0;
	wave_transf_des wavelet[3];

	g_assert(fit->naxes[2] <= 3);

	float *Imag;
	if (fit->type == DATA_USHORT) {
		Imag = f_vector_alloc(fit->ry * fit->rx);
		if (Imag == NULL) {
			PRINT_ALLOC_ERR;
			return 1;
		}
	}

	if (reqlayer < 0 || reqlayer > 3) {
		start = 0;
		end = fit->naxes[2];
	}
	else {
		start = reqlayer;
		end = start + 1;
	}

	for (chan = start; chan < end; chan++) {
		int Nl, Nc;

		if (fit->type == DATA_USHORT) {
			if (wavelet_transform(Imag, fit->ry, fit->rx, &wavelet[chan],
						Type, Nbr_Plan, fit->pdata[chan])) {
				retval = 1;
				break;
			}
		}
		else if (fit->type == DATA_FLOAT) {
			Imag = fit->fpdata[chan];
			if (wavelet_transform_float(Imag, fit->ry, fit->rx, &wavelet[chan],
						Type, Nbr_Plan)) {
				retval = 1;
				break;
			}
		}
		Nl = wavelet[chan].Nbr_Ligne;
		Nc = wavelet[chan].Nbr_Col;
		pave_2d_extract_plan(wavelet[chan].Pave.Data, Imag, Nl, Nc, Plan);
		if (fit->type == DATA_USHORT)
			reget_rawdata(Imag, Nl, Nc, fit->pdata[chan]);
		wave_io_free(&wavelet[chan]);
	}

	/* Free */
	if (fit->type == DATA_USHORT)
		free(Imag);
	return retval;
}

static gboolean end_wavelets_filter(gpointer p) {
	struct wavelets_filter_data *args = (struct wavelets_filter_data *) p;
	stop_processing_thread();// can it be done here in case there is no thread?
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_DONE);
	
	set_cursor_waiting(FALSE);
	free(args);
	return FALSE;
}

gpointer extract_plans(gpointer p) {
	int i;
	fits fit = { 0 };
	struct wavelets_filter_data *args = (struct wavelets_filter_data *) p;

	copyfits(args->fit, &fit, CP_ALLOC | CP_COPYA | CP_FORMAT, 0);

	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);

	for (i = 0; i < args->Nbr_Plan; i++) {
		gchar *filename, *msg;

		filename = g_strdup_printf("layer%02d", i);
		msg = g_strdup_printf(_("Extracting %s..."), filename);
		set_progress_bar_data(msg, (float)i / args->Nbr_Plan);
		get_wavelet_layers(&fit, args->Nbr_Plan, i, args->Type, -1);
		savefits(filename, &fit);
		g_free(filename);
		g_free(msg);
	}
	clearfits(&fit);
	siril_add_idle(end_wavelets_filter, args);
	return GINT_TO_POINTER(0);
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

void on_wavelets_dialog_hide(GtkWidget *widget, gpointer user_data) {
	gtk_widget_set_sensitive(lookup_widget("frame_wavelets"), FALSE);
	gtk_widget_set_sensitive(lookup_widget("button_reset_w"), FALSE);
	clearfits(&wavelets_gfit_backup);
}

void on_button_reset_w_clicked(GtkButton *button, gpointer user_data) {
	reset_scale_w();
	update_wavelets();
}

void apply_wavelets_cancel() {
	if (gtk_widget_get_sensitive(lookup_widget("frame_wavelets")) == TRUE) {
		reset_scale_w();
		update_wavelets();
	}
}

void on_button_ok_w_clicked(GtkButton *button, gpointer user_data) {
	if (gtk_widget_get_sensitive(lookup_widget("frame_wavelets")) == TRUE) {
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
	gtk_widget_set_sensitive(lookup_widget("frame_wavelets"), TRUE);
	gtk_widget_set_sensitive(lookup_widget("button_reset_w"), TRUE);
	set_cursor_waiting(FALSE);
	return;
}

/****************** GUI for Wavelet Layers Extraction *****************/

void on_menu_wavelet_separation_activate(GtkMenuItem *menuitem,
		gpointer user_data) {

	if (single_image_is_loaded()) {
		siril_open_dialog("extract_wavelets_layers_dialog");
	}
}

void on_button_extract_w_ok_clicked(GtkButton *button, gpointer user_data) {
	int Nbr_Plan, Type, maxplan, mins;
	static GtkSpinButton *Spin_Nbr_Plan = NULL;
	static GtkComboBox *Combo_Wavelets_Type = NULL;

	if (Spin_Nbr_Plan == NULL) {
		Spin_Nbr_Plan = GTK_SPIN_BUTTON(lookup_widget("spinbutton_extract_w"));
		Combo_Wavelets_Type = GTK_COMBO_BOX(
				lookup_widget("combo_interpolation_extract_w"));
	}

	Nbr_Plan = gtk_spin_button_get_value(Spin_Nbr_Plan);
	Type = gtk_combo_box_get_active(Combo_Wavelets_Type) + 1;// 1: linear, 2: bspline

	set_cursor_waiting(TRUE);
	mins = min(gfit.rx, gfit.ry);
	maxplan = log(mins) / log(2) - 2;

	if (Nbr_Plan > maxplan) {
		char *msg = siril_log_message(_("Wavelet: maximum number "
				"of plans for this image size is %d\n"), maxplan);
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Warning"), msg);
		set_cursor_waiting(FALSE);
		return;
	}

	struct wavelets_filter_data *args = malloc(sizeof(struct wavelets_filter_data));

	args->Type = Type;
	args->Nbr_Plan = Nbr_Plan;
	args->fit = &gfit;
	start_in_new_thread(extract_plans, args);
}

void on_button_extract_w_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("extract_wavelets_layers_dialog");
}

void on_spinbutton_plans_w_value_changed(GtkSpinButton *button, gpointer user_data) {
	int i;
	gint current_value = gtk_spin_button_get_value_as_int(button);
	for (i = 0; i < 6; i++) {
		gchar *tmp = g_strdup_printf("box_w%d", i);
		gtk_widget_set_visible(lookup_widget(tmp), current_value > i);
		g_free(tmp);
	}
}

void on_spin_w0_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[0] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	notify_update((gpointer) param);
}

void on_spin_w1_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[1] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	notify_update((gpointer) param);
}

void on_spin_w2_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[2] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	notify_update((gpointer) param);
}

void on_spin_w3_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[3] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	notify_update((gpointer) param);
}

void on_spin_w4_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[4] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	notify_update((gpointer) param);
}

void on_spin_w5_value_changed(GtkSpinButton *button, gpointer user_data) {
	wavelet_value[5] = gtk_spin_button_get_value(button);
	update_image *param = malloc(sizeof(update_image));
	param->update_preview_fn = update_wavelets;
	notify_update((gpointer) param);
}

