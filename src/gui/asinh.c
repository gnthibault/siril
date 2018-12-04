#include "core/siril.h"
#include "core/proto.h"
#include "io/single_image.h"
#include "callbacks.h"
#include "progress_and_log.h"
#include "core/undo.h"

static gboolean asinh_rgb_space = FALSE;
static double asinh_stretch_value = 1.0, asinh_black_value = 0.0;
static fits asinh_gfit_backup;

static void asinh_startup() {
	copyfits(&gfit, &asinh_gfit_backup, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
}

static void asinh_close(gboolean revert) {
	if (revert) {
		copyfits(&asinh_gfit_backup, &gfit, CP_COPYA, -1);
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	} else {
		undo_save_state(&asinh_gfit_backup, "Processing: Asinh Transformation: (stretch=%6.1lf, bp=%7.5lf)",
				asinh_stretch_value, asinh_black_value);
	}
	clearfits(&asinh_gfit_backup);
}

static void asinh_recompute() {
	set_cursor_waiting(TRUE);
	copyfits(&asinh_gfit_backup, &gfit, CP_COPYA, -1);
	asinhlut(&gfit, asinh_stretch_value, asinh_black_value, asinh_rgb_space);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void on_menuitem_asinh_activate(GtkMenuItem *menuitem, gpointer user_data) {
	asinh_startup();
	asinh_stretch_value = 1.0;
	asinh_black_value = 0.0;
	asinh_rgb_space = FALSE;
	gtk_toggle_button_set_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_RGBspace")), asinh_rgb_space);
	gtk_range_set_value(GTK_RANGE(lookup_widget("scale_asinh")), asinh_stretch_value);
	gtk_range_set_value(GTK_RANGE(lookup_widget("black_point_asinh")), asinh_black_value);
	gtk_widget_show(lookup_widget("asinh_dialog"));
}

void on_asinh_cancel_clicked(GtkButton *button, gpointer user_data) {
	asinh_close(TRUE);
	gtk_widget_hide(lookup_widget("asinh_dialog"));
}

void on_asinh_ok_clicked(GtkButton *button, gpointer user_data) {
	asinh_close(FALSE);
	gtk_widget_hide(lookup_widget("asinh_dialog"));
}

void on_asinh_dialog_close(GtkDialog *dialog, gpointer user_data) {
	asinh_close(FALSE);
}

gboolean on_scale_asinh_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	asinh_stretch_value = gtk_range_get_value(GTK_RANGE(widget));
	asinh_recompute();
	return FALSE;
}

gboolean on_black_point_asinh_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	asinh_black_value = gtk_range_get_value(GTK_RANGE(widget));
	asinh_recompute();
	return FALSE;
}

void on_asinh_RGBspace_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	asinh_rgb_space = gtk_toggle_button_get_active(togglebutton);
	asinh_recompute();
}
