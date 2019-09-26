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

#if !GLIB_CHECK_VERSION(2,62,0)
/**
 * g_date_time_format_iso8601:
 * @datetime: A #GDateTime
 *
 * Format @datetime in [ISO 8601 format](https://en.wikipedia.org/wiki/ISO_8601),
 * including the date, time and time zone, and return that as a UTF-8 encoded
 * string.
 *
 * Returns: a newly allocated string formatted in ISO 8601 format
 *     or %NULL in the case that there was an error. The string
 *     should be freed with g_free().
 * Since: 2.62
 */
gchar* g_date_time_format_iso8601(GDateTime *datetime) {
	GString *outstr = NULL;
	gchar *main_date = NULL;
	time_t offset;

	/* Main date and time. */
	main_date = g_date_time_format(datetime, "%Y-%m-%dT%H:%M:%S");
	outstr = g_string_new(main_date);
	g_free(main_date);

	/* Timezone. Format it as `%:::z` unless the offset is zero, in which case
	 * we can simply use `Z`. */
	offset = g_date_time_get_utc_offset(datetime);

	if (offset == 0) {
		g_string_append_c(outstr, 'Z');
	} else {
		gchar *time_zone = g_date_time_format(datetime, "%:::z");
		g_string_append(outstr, time_zone);
		g_free(time_zone);
	}

	return g_string_free(outstr, FALSE);
}
#endif

static gchar* build_timestamp() {
#if !GLIB_CHECK_VERSION(2,26,0)
	GTimeVal time;

	g_get_current_time(&time);
	return g_time_val_to_iso8601(&time);
#else
	time_t t = time(NULL);
	GDateTime *dt = g_date_time_new_from_unix_utc(t);
	if (dt) {
		gchar *iso8601_string = g_date_time_format_iso8601(dt);
		g_date_time_unref(dt);
		return iso8601_string;
	} else {
		return "";
	}
#endif
}

static void replace_not_valid_char(gchar *str, gchar c, gchar n) {
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
	gchar *str, **token;
	guint nargs, i;
	FILE *f;

	tv = GTK_TEXT_VIEW(lookup_widget("output"));
	log = gtk_text_view_get_buffer(tv);
	gtk_text_buffer_get_bounds(log, &start, &end);
	str = gtk_text_buffer_get_text(log, &start, &end, FALSE);

	token = g_strsplit(str, "\n", -1);
	nargs = g_strv_length(token);

	f = g_fopen(filename, "w");

	for (i = 0; i < nargs; i++) {
		fprintf(f, "%s%s", token[i], SIRIL_EOL);
	}

	fclose(f);
	g_free(str);
	g_strfreev(token);
}

static void set_filter(GtkFileChooser *dialog) {
	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, _("Log files (*.log)"));
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
	replace_not_valid_char(filename, ':', '.');
	filename = str_append(&filename, ".log");

	widgetdialog = siril_file_chooser_save(control_window, GTK_FILE_CHOOSER_ACTION_SAVE);
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

		g_free(file);
	}
	siril_widget_destroy(widgetdialog);
	g_free(filename);
}

/************** Callbacks function ***********/

void on_clear_log_button_clicked(GtkButton *button, gpointer user_data) {
	process_clear(0);
}

void on_export_log_button_clicked(GtkButton *button, gpointer user_data) {
	save_log_dialog();
}
