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

#include "core/siril.h"
#include "core/proto.h"
#include "core/command.h"
#include "core/undo.h"
#include "core/siril_update.h"
#include "core/siril_cmd_help.h"
#include "gui/about_dialog.h"
#include "gui/callbacks.h"
#include "gui/open_dialog.h"
#include "gui/save_dialog.h"
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/script_menu.h"

#include "siril_actions.h"

void open_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	header_open_button_clicked();
}

void cwd_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	cwd_btton_clicked();
}

void save_as_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	on_header_save_button_clicked();
}

void undo_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	set_cursor_waiting(TRUE);
	undo_display_data(UNDO);
	set_cursor_waiting(FALSE);
	update_MenuItem();
}

void redo_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	set_cursor_waiting(TRUE);
	undo_display_data(REDO);
	set_cursor_waiting(FALSE);
	update_MenuItem();
}

void quit_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_quit();
}

void about_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_show_about_dialog();
}

void preferences_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_open_dialog("settings_window");
}

void close_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	process_close(0);
}

void scripts_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_get_on_script_pages();
}

#ifdef HAVE_LIBCURL
void updates_action_activate(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	siril_check_updates();
}
#endif

void full_screen_activated(GSimpleAction *action, GVariant *parameter,
		gpointer user_data) {
	GtkApplication *app;
	GtkWindow *window;
	GtkWidget *toolbarbox = lookup_widget("toolbarbox");
	gboolean is_fullscreen;

	app = GTK_APPLICATION(user_data);
	window = GTK_WINDOW(gtk_application_get_active_window(app));

	GdkWindow *gdk_window = gtk_widget_get_window(GTK_WIDGET(window));
	is_fullscreen = gdk_window_get_state(gdk_window) & GDK_WINDOW_STATE_FULLSCREEN;

	if (is_fullscreen) {
		gtk_window_unfullscreen(window);
	} else {
		gtk_window_fullscreen(window);
	}
	gtk_widget_set_visible(toolbarbox, is_fullscreen);
}

void keyboard_shortcuts_activated(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	GtkApplication *app;
	GtkWindow *window;

	app = GTK_APPLICATION(user_data);
	window = GTK_WINDOW(gtk_application_get_active_window(app));

	siril_cmd_help_keyboard_shortcuts(window);
}

void tab_conversion_activate(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(FILE_CONVERSION);
}

void tab_sequence_activate(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(IMAGE_SEQ);
}

void tab_prepro_activate(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(PRE_PROC);
}

void tab_registration_activate(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(REGISTRATION);
}

void tab_plot_activate(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(PLOT);
}

void tab_stacking_activate(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(STACKING);
}

void tab_logs_activate(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	control_window_switch_to_tab(OUTPUT_LOGS);
}

void toolbar_activate(GSimpleAction *action,
		GVariant *parameter, gpointer user_data) {
	GtkWidget *w = lookup_widget("toolbarbox");
	gtk_widget_set_visible(w, !gtk_widget_get_visible(w));
}
