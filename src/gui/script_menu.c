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

#ifdef _WIN32
#include <Windows.h>
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
#include "core/command.h"
#include "core/processing.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "algos/sorting.h"
#include "script_menu.h"

#define SCRIPT_EXT ".ssf"

static GSList *initialize_script_paths(){
	GSList *list = NULL;
#ifdef _WIN32
	list = g_slist_prepend(list, g_build_filename(get_special_folder(CSIDL_APPDATA), "siril",
					"scripts", NULL));

	gchar *execpath = g_win32_get_package_installation_directory_of_module(NULL);

	list = g_slist_prepend(list, g_build_filename(execpath, "scripts", NULL));
	g_free(execpath);
#else
	list = g_slist_prepend(list, g_build_filename(PACKAGE_DATA_DIR, "scripts", NULL));
	list = g_slist_prepend(list, g_build_filename(g_get_home_dir(), ".siril", "scripts", NULL));
	list = g_slist_prepend(list, g_build_filename(g_get_home_dir(), "siril", "scripts", NULL));
#endif
	list = g_slist_reverse(list);
	return list;
}

static GSList *get_list_from_textview() {
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

static void fill_gtkText(GSList *list) {
	while (list) {
		add_path_to_gtkText((gchar *) list->data);
		list = list->next;
	}
}

static GSList *search_script(const char *path) {
	GSList *list = NULL;
	GDir *dir;
	GError *error = NULL;
	const gchar *file;

	dir = g_dir_open(path, 0, &error);
	if (!dir) {
		fprintf(stderr, "scripts: %s\n", error->message);
		g_error_free(error);
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
	gchar *script_file;
	GString *str;

	if (get_thread_run()) {
		siril_log_message(_("Another task is already in "
				"progress, ignoring new request.\n"));
		return;
	}

	if (com.script_thread)
		g_thread_join(com.script_thread);

	/* append script file extension */
	str = g_string_new((gchar *) user_data);
	str = g_string_append(str, SCRIPT_EXT);
	script_file = g_string_free(str, FALSE);

	FILE* fp = g_fopen(script_file, "r");
	if (fp == NULL) {
		siril_log_message(_("File [%s] does not exist\n"), script_file);
		g_free(script_file);
		return;
	}
	/* Switch to console tab */
	control_window_switch_to_tab(OUTPUT_LOGS);
	/* ensure that everything is closed */
	process_close(0);
	/* Then, run script */
	siril_log_message(_("Starting script %s\n"), script_file);
	com.script_thread = g_thread_new("script", execute_script, fp);

	g_free(script_file);
}

int initialize_script_menu() {
	static GtkWidget *menuscript = NULL;
	GSList *list, *script, *s;
	GtkWidget *menu;
	gint nb_item = 0;

	if (!menuscript) {
		menuscript = lookup_widget("menuscript");
	}

	if (!com.script_path) {
		script = initialize_script_paths();
		com.script_path = script;
	} else {
		script = com.script_path;
	}
	fill_gtkText(script);

	menu = gtk_menu_new();

	for (s = script; s; s = s->next) {
		list = search_script(s->data);
		if (list) {
			GSList *l;
			gtk_widget_show(menuscript);
			gtk_menu_item_set_submenu(GTK_MENU_ITEM(menuscript), menu);
			/* write separator but not for the first one */
			if (nb_item != 0) {
				GtkWidget *separator = gtk_separator_menu_item_new();
				gtk_menu_shell_append(GTK_MENU_SHELL(menu), separator);
				gtk_widget_show(separator);
			}
			siril_log_message(_("Searching scripts in: \"%s\"...\n"), s->data);
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

void fill_script_paths_list() {
	GSList *list;

	g_slist_free_full(com.script_path, g_free);
	list = get_list_from_textview();
	com.script_path = list;
	writeinitfile();
}


/* Get Scripts menu */

#define GET_SCRIPTS_URL "https://free-astro.org/index.php?title=Siril:scripts"

void on_help_get_scripts_activate(GtkMenuItem *menuitem, gpointer user_data) {
	gboolean ret;
	const char *locale = setlocale(LC_MESSAGES, NULL);
	const char *supported_languages[] = { "el", "fr", "it", NULL }; // en is NULL: default language
	gchar *lang = NULL;
	int i = 0;

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
	ret = gtk_show_uri_on_window(GTK_WINDOW(win), url,
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

