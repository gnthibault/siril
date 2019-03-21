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
#include "gui/callbacks.h"

static gchar *build_timestamp() {
	GTimeVal time;

	g_get_current_time(&time);
	return g_time_val_to_iso8601(&time);
}

static void remove_not_valid_char(gchar *str, gchar c, gchar n) {
	gchar *s = str;
	while (*s) {
		if (*s == c) {
			*s = n;
		}
		s++;
	}
}

static void save_log_file(gchar *filename) {
	GtkTextBuffer *log;
	GtkTextView *tv;
	GtkTextIter start, end;
	gchar *str;
	FILE *f;

	tv = GTK_TEXT_VIEW(lookup_widget("output"));
	log = gtk_text_view_get_buffer(tv);
	gtk_text_buffer_get_bounds(log, &start, &end);
	str = gtk_text_buffer_get_text(log, &start, &end, FALSE);

	f = g_fopen(filename, "w");
	fprintf(f, "%s", str);
	fclose(f);
	g_free(str);
}

static SirilWidget *siril_file_chooser_save_log(GtkWindow *parent, GtkFileChooserAction action) {
#if (defined _WIN32) || (defined(__APPLE__) && defined(__MACH__))
	return gtk_file_chooser_native_new(_("Save File"), parent, action,
			_("_Save"), _("_Cancel"));
#else
	return gtk_file_chooser_dialog_new(_("Save File"), parent, action,
				_("_Cancel"), GTK_RESPONSE_CANCEL, _("_Save"), GTK_RESPONSE_ACCEPT,
				NULL);
#endif
}

static gint siril_dialog_run(SirilWidget *widgetdialog) {
#if (defined _WIN32) || (defined(__APPLE__) && defined(__MACH__))
	return gtk_native_dialog_run(GTK_NATIVE_DIALOG(widgetdialog));
#else
	return gtk_dialog_run(GTK_DIALOG(GTK_FILE_CHOOSER(widgetdialog)));
#endif
}

static void siril_widget_destroy(SirilWidget *widgetdialog) {
#if (defined _WIN32) || (defined(__APPLE__) && defined(__MACH__))
	g_object_unref(widgetdialog);
#else
	gtk_widget_destroy(widgetdialog);
#endif
}

static void set_filter(GtkFileChooser *dialog) {
	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, "log files (*.log)");
	gtk_file_filter_add_pattern(f, "*.log");
	gtk_file_chooser_add_filter(dialog, f);
	gtk_file_chooser_set_filter(dialog, f);
}

static void save_log_dialog() {
	SirilWidget *widgetdialog;
	GtkFileChooser *dialog = NULL;
	GtkWindow *control_window = GTK_WINDOW(lookup_widget("control_window"));
	gint res;
	gchar *filename;

	filename = build_timestamp();
	remove_not_valid_char(filename, ':', '.');
	filename = str_append(&filename, ".log");

	widgetdialog = siril_file_chooser_save_log(control_window, GTK_FILE_CHOOSER_ACTION_SAVE);
	dialog = GTK_FILE_CHOOSER(widgetdialog);
	gtk_file_chooser_set_current_folder(dialog, com.wd);
	gtk_file_chooser_set_select_multiple(dialog, FALSE);
	gtk_file_chooser_set_do_overwrite_confirmation(dialog, TRUE);
	gtk_file_chooser_set_current_name(dialog, filename);
	set_filter(dialog);

	res = siril_dialog_run(widgetdialog);
	if (res == GTK_RESPONSE_ACCEPT) {
		gchar *file = gtk_file_chooser_get_filename(dialog);
		save_log_file(file);
	}
	siril_widget_destroy(widgetdialog);
	g_free(filename);
}

static void export_siril_log() {
	save_log_dialog();
}


/************** Callbacks function ***********/

void on_clear_log_button_clicked(GtkButton *button, gpointer user_data) {
	process_clear(0);
}

void on_export_log_button_clicked(GtkButton *button, gpointer user_data) {
	export_siril_log();
}
