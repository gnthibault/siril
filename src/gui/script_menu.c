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

#ifdef _WIN32
#include <windows.h>
/* Constant available since Shell32.dll 4.72 */
#ifndef CSIDL_APPDATA
#define CSIDL_APPDATA 0x001a
#endif
#endif
#include <string.h>
#include <locale.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/command.h" // for process_close()
#include "core/command_line_processor.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_app_dirs.h"
#include "gui/utils.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "algos/sorting.h"
#include "script_menu.h"

#define SCRIPT_EXT ".ssf"
#define CONFIRM_RUN_SCRIPTS _("You are about to use scripts. Running automatic scripts is something that is easy and generally it provides a nice image. However you have to keep in mind that scripts are not magic; automatic choices are made where human decision would probably be better. Also, every commands used in a script are available on the interface with a better parameter control.")

static GSList *initialize_script_paths(){
	GSList *list = NULL;
#ifdef _WIN32
	list = g_slist_prepend(list, g_build_filename(get_special_folder(CSIDL_APPDATA), "siril",
					"scripts", NULL));

	gchar *execpath = g_win32_get_package_installation_directory_of_module(NULL);

	list = g_slist_prepend(list, g_build_filename(execpath, "scripts", NULL));
	g_free(execpath);
#else
	list = g_slist_prepend(list, g_build_filename(siril_get_system_data_dir(), "scripts", NULL));
	list = g_slist_prepend(list, g_build_filename(g_get_home_dir(), ".siril", "scripts", NULL));
	list = g_slist_prepend(list, g_build_filename(g_get_home_dir(), "siril", "scripts", NULL));
#endif
	list = g_slist_reverse(list);
	return list;
}

static void add_path_to_gtkText(gchar *path) {
	static GtkTextBuffer *tbuf = NULL;
	static GtkTextView *text = NULL;
	GtkTextIter iter;

	if (!tbuf) {
		text = GTK_TEXT_VIEW(lookup_widget("GtkTxtScriptPath"));
		tbuf = gtk_text_view_get_buffer(text);
	}

	gtk_text_buffer_get_end_iter(tbuf, &iter);
	gtk_text_buffer_insert(tbuf, &iter, path, strlen(path));
	gtk_text_buffer_insert(tbuf, &iter, "\n", strlen("\n"));

	/* scroll to end */
	gtk_text_buffer_get_end_iter(tbuf, &iter);
	GtkTextMark *insert_mark = gtk_text_buffer_get_insert(tbuf);
	gtk_text_buffer_place_cursor(tbuf, &iter);
	gtk_text_view_scroll_to_mark(text, insert_mark, 0.0, TRUE, 0.0, 1.0);
}

static void clear_gtk_list() {
	GtkTextView *text = GTK_TEXT_VIEW(lookup_widget("GtkTxtScriptPath"));
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(text);
	GtkTextIter start_iter, end_iter;
	gtk_text_buffer_get_start_iter(tbuf, &start_iter);
	gtk_text_buffer_get_end_iter(tbuf, &end_iter);
	gtk_text_buffer_delete(tbuf, &start_iter, &end_iter);
}

static GSList *search_script(const char *path) {
	GSList *list = NULL;
	GDir *dir;
	GError *error = NULL;
	const gchar *file;

	dir = g_dir_open(path, 0, &error);
	if (!dir) {
		fprintf(stderr, "scripts: %s\n", error->message);
		g_clear_error(&error);
		return NULL;
	}
	while ((file = g_dir_read_name(dir)) != NULL) {
		if (g_str_has_suffix(file, SCRIPT_EXT)) {
			gchar *str = (gchar*) remove_ext_from_filename(file);

			list = g_slist_prepend(list, str);
		}
	}
	list = g_slist_sort(list, (GCompareFunc) strcompare);
	g_dir_close(dir);

	return list;
}

static void on_script_execution(GtkMenuItem *menuitem, gpointer user_data) {
	GError *error = NULL;

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	if (com.pref.save.warn_script) {
		gboolean dont_show_again;
		gboolean confirm = siril_confirm_dialog_and_remember(
				_("Please read me before using scripts"), CONFIRM_RUN_SCRIPTS, _("Run Script"), &dont_show_again);
		com.pref.save.warn_script = !dont_show_again;
		/* We do not use set_GUI_misc because some button state can be in an unsaved state if the
		 * preference dialog is opened
		 */
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("miscAskScript")), com.pref.save.warn_script);
		/* update config file */
		writeinitfile();
		if (!confirm) {
			return;
		}
	}

	if (com.script_thread)
		g_thread_join(com.script_thread);

	/* Switch to console tab */
	control_window_switch_to_tab(OUTPUT_LOGS);

	gchar *script_file = g_strdup_printf("%s%s", (gchar *) user_data, SCRIPT_EXT);
	GFile *file = g_file_new_for_path(script_file);

    GFileInfo *info;

	info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_SIZE,
			G_FILE_QUERY_INFO_NONE, NULL, &error);
	if (info) {
		GInputStream *input_stream = (GInputStream*) g_file_read(file, NULL, &error);

		if (input_stream == NULL) {
			if (error != NULL) {
				g_clear_error(&error);
				siril_log_message(_("File [%s] does not exist\n"), script_file);
			}

			g_free(script_file);
			g_object_unref(file);
			return;
		}
		/* ensure that everything is closed */
		process_close(0);
		/* Then, run script */
		siril_log_message(_("Starting script %s\n"), script_file);
		com.script_thread = g_thread_new("script", execute_script, input_stream);
	}

	g_free(script_file);
	g_object_unref(file);
}

int initialize_script_menu() {
	static GtkWidget *menuscript = NULL;
	GSList *list, *script, *s;
	GtkWidget *menu;
	gint nb_item = 0;

	if (!menuscript) {
		menuscript = lookup_widget("header_scripts_button");
	}
	
	script = set_list_to_preferences_dialog(com.pref.script_path);

	menu = gtk_menu_new();
	gtk_widget_hide(menuscript);

	for (s = script; s; s = s->next) {
		list = search_script(s->data);
		if (list) {
			GSList *l;
			if (!gtk_widget_get_visible(menuscript)) {
				gtk_widget_show(menuscript);
				gtk_menu_button_set_popup(GTK_MENU_BUTTON(menuscript), menu);
			}
			/* write separator but not for the first one */
			if (nb_item != 0) {
				GtkWidget *separator = gtk_separator_menu_item_new();
				gtk_menu_shell_append(GTK_MENU_SHELL(menu), separator);
				gtk_widget_show(separator);
			}
			siril_log_color_message(_("Searching scripts in: \"%s\"...\n"), "green", s->data);
			for (l = list; l; l = l->next) {
				nb_item ++;
				/* write an item per script file */
				GtkWidget *menu_item;

				menu_item = gtk_menu_item_new_with_label(l->data);
				gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
				gchar *full_path = g_build_filename(s->data, l->data,
						NULL);
				g_signal_connect(G_OBJECT(menu_item), "activate",
						G_CALLBACK(on_script_execution), (gchar * ) full_path);
				siril_log_message(_("Loading script: %s\n"), l->data);
				gtk_widget_show(menu_item);
			}
			g_slist_free_full(list, g_free);
		}
	}

	writeinitfile();

	return 0;
}

int refresh_scripts(gboolean update_list, gchar **error) {
	gchar *err = NULL;
	int retval;
	GSList *list = get_list_from_preferences_dialog();
	if (list == NULL) {
		err = siril_log_color_message(_("Cannot refresh the scripts if the list is empty.\n"), "red");
		retval = 1;
	} else {
		g_slist_free_full(com.pref.script_path, g_free);
		com.pref.script_path = list;
		retval = initialize_script_menu();
	}
	if (error) {
		*error = err;
	}
	return retval;
}

GSList *get_list_from_preferences_dialog() {
	GSList *list = NULL;
	static GtkTextBuffer *tbuf = NULL;
	static GtkTextView *text = NULL;
	GtkTextIter start, end;
	gchar *txt;
	gint i = 0;

	if (!tbuf) {
		text = GTK_TEXT_VIEW(lookup_widget("GtkTxtScriptPath"));
		tbuf = gtk_text_view_get_buffer(text);
	}
	gtk_text_buffer_get_bounds(tbuf, &start, &end);
	txt = gtk_text_buffer_get_text(tbuf, &start, &end, TRUE);
	if (txt) {
		gchar **token = g_strsplit(txt, "\n", -1);
		while (token[i]) {
			if (*token[i] != '\0')
				list = g_slist_prepend(list, g_strdup(token[i]));
			i++;
		}
		list = g_slist_reverse(list);
		g_strfreev(token);
	}

	return list;
}

GSList *set_list_to_preferences_dialog(GSList *list) {
	clear_gtk_list();
	if (list == NULL) {
		list = initialize_script_paths();
	}
	for (GSList *l = list; l; l = l->next) {
		add_path_to_gtkText((gchar *) l->data);
	}
	return list;
}

/* Get Scripts menu */

#define GET_SCRIPTS_URL "https://free-astro.org/index.php?title=Siril:scripts"

void siril_get_on_script_pages() {
	gboolean ret;
	const char *locale;
	const char *supported_languages[] = { "fr", NULL }; // en is NULL: default language
	gchar *lang = NULL;
	int i = 0;

	if (!g_strcmp0(com.pref.combo_lang, "")) {
		locale = setlocale(LC_MESSAGES, NULL);
	} else {
		locale = com.pref.combo_lang;
	}

	if (locale) {
		while (supported_languages[i]) {
			if (!strncmp(locale, supported_languages[i], 2)) {
				lang = g_strndup(locale, 2);
				break;
			}
			i++;
		}
	}
	gchar *url = g_build_path("/", GET_SCRIPTS_URL, lang, NULL);

#if GTK_CHECK_VERSION(3, 22, 0)
	GtkWidget* win = lookup_widget("control_window");
	ret = gtk_show_uri_on_window(GTK_WINDOW(GTK_APPLICATION_WINDOW(win)), url,
			gtk_get_current_event_time(), NULL);
#else
	ret = gtk_show_uri(gdk_screen_get_default(), url,
			gtk_get_current_event_time(), NULL);
#endif
	if (!ret) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Could not show link"),
				_("Please go to <a href=\""GET_SCRIPTS_URL"\">"GET_SCRIPTS_URL"</a> "
								"by copying the link."));
	}
	g_free(url);
	g_free(lang);
}

