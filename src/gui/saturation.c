#include <stdlib.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/undo.h"
#include "core/processing.h"
#include "algos/colors.h"
#include "io/single_image.h"
#include "progress_and_log.h"
#include "callbacks.h"
#include "dialogs.h"
#include "saturation.h"

static gboolean satu_preserve_bkg = TRUE;
static double satu_amount = 0.0;
static int satu_hue_type = 6;
static fits satu_gfit_backup;

static void satu_startup() {
	copyfits(&gfit, &satu_gfit_backup, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
}

static void satu_close(gboolean revert) {
	if (revert) {
		copyfits(&satu_gfit_backup, &gfit, CP_COPYA, -1);
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	} else {
		undo_save_state(&satu_gfit_backup, "Processing: Saturation enhancement (amount=%4.2lf)", satu_amount);
	}
	clearfits(&satu_gfit_backup);
}

void on_satu_cancel_clicked(GtkButton *button, gpointer user_data) {
	satu_close(TRUE);
	siril_close_dialog("satu_dialog");
}

void on_satu_apply_clicked(GtkButton *button, gpointer user_data) {
	apply_satu_changes();
	siril_close_dialog("satu_dialog");
}

void on_satu_dialog_close(GtkDialog *dialog, gpointer user_data) {
	if (satu_amount != 0.0)
		apply_satu_changes();
}

void satu_recompute() {
	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return;
	}
	if (satu_amount == 0.0) return;
	set_cursor_waiting(TRUE);

	struct enhance_saturation_data *args = malloc(sizeof(struct enhance_saturation_data));

	switch (satu_hue_type) {
	case 0:		// Pink-Red to Red-Orange
		args->h_min = 346.0;
		args->h_max = 20.0;
		break;
	case 1:		// Orange-Brown to Yellow
		args->h_min = 21.0;
		args->h_max = 60.0;
		break;
	case 2:		// Yellow-Green to Green-Cyan
		args->h_min = 61.0;
		args->h_max = 200.0;
		break;
	case 3:		// Cyan
		args->h_min = 170.0;
		args->h_max = 200.0;
		break;
	case 4:		// Cyan-Blue to Blue-Magenta
		args->h_min = 201.0;
		args->h_max = 280.0;
		break;
	case 5:		// Magenta to Pink
		args->h_min = 281.0;
		args->h_max = 345.0;
		break;
	default:
	case 6:		// Global
		args->h_min = 0.0;
		args->h_max = 360.0;
	}

	args->input = &satu_gfit_backup;
	args->output = &gfit;
	args->coeff = satu_amount;
	args->preserve = satu_preserve_bkg;
	start_in_new_thread(enhance_saturation, args);
}

void on_menuitem_satu_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (!single_image_is_loaded() || !isrgb(&gfit))
		return;
	satu_startup();
	satu_amount = 0.0;
	satu_hue_type = 6;
	satu_preserve_bkg = TRUE;

	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_saturation")), satu_hue_type);
	gtk_range_set_value(GTK_RANGE(lookup_widget("scale_satu")), satu_amount);
	gtk_toggle_button_set_active(
			GTK_TOGGLE_BUTTON(lookup_widget("preserve_bg")), satu_preserve_bkg);
	siril_open_dialog("satu_dialog");
}

gboolean on_scale_satu_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	satu_amount = gtk_range_get_value(GTK_RANGE(widget));
	satu_recompute();
	return FALSE;
}

gboolean on_scale_satu_key_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	satu_amount = gtk_range_get_value(GTK_RANGE(widget));
	satu_recompute();
	return FALSE;
}

void on_preserve_bg_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	satu_preserve_bkg = gtk_toggle_button_get_active(togglebutton);
	satu_recompute();
}

void on_combo_saturation_changed(GtkComboBox* box, gpointer user_data) {
	satu_hue_type = gtk_combo_box_get_active(box);
	satu_recompute();
}

void on_satu_undo_clicked(GtkButton *button, gpointer user_data) {
	satu_preserve_bkg = TRUE;
	satu_amount = 0.0;
	GtkToggleButton *check_button = GTK_TOGGLE_BUTTON(lookup_widget("preserve_bg"));
	g_signal_handlers_block_by_func(check_button, on_preserve_bg_toggled, NULL);
	gtk_toggle_button_set_active(check_button, satu_preserve_bkg);
	g_signal_handlers_unblock_by_func(check_button, on_preserve_bg_toggled, NULL);
	gtk_range_set_value(GTK_RANGE(lookup_widget("scale_satu")), satu_amount);
	copyfits(&satu_gfit_backup, &gfit, CP_COPYA, -1);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
}

void apply_satu_changes() {
	gboolean status = satu_amount != 0.0;
	satu_close(!status);
}
