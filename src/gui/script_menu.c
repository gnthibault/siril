/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2018 team free-astro (see more in AUTHORS file)
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
#include "core/siril.h"
#include "core/proto.h"
#include "core/command.h"
#include "core/processing.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "script_menu.h"

#define EXT ".ssf"

static GThread *script_thread = NULL;

static GSList *search_script(const char *path) {
	GSList *list = NULL;
	GDir *dir;
	GError *error = NULL;
	const gchar *file;

	dir = g_dir_open(path, 0, &error);
	if (!dir) {
		fprintf(stderr, "scripts: %s\n", error->message);
		return NULL;
	}
	while ((file = g_dir_read_name(dir)) != NULL) {
		if (g_str_has_suffix(file, EXT)) {
			char *str = remove_ext_from_filename(file);

			list = g_slist_append(list, str);
		}
	}
	list = g_slist_sort(list, (GCompareFunc) strcompare);
	g_dir_close(dir);

	return list;
}

void on_script_execution(GtkMenuItem *menuitem, gpointer user_data) {
	gchar *script_file;
	GString *str;

	if (get_thread_run()) {
		siril_log_message(_("Another task is already in "
				"progress, ignoring new request.\n"));
		return;
	}

	if (script_thread)
		g_thread_join(script_thread);

	/* append script file extension */
	str = g_string_new((gchar *) user_data);
	str = g_string_append(str, EXT);
	script_file = g_string_free(str, FALSE);

	FILE* fp = g_fopen(script_file, "r");
	if (fp == NULL) {
		siril_log_message(_("File [%s] does not exist\n"), script_file);
		g_free(script_file);
		return;
	}
	control_window_switch_to_tab(OUTPUT_LOGS);
	siril_log_message(_("Starting script %s\n"), script_file);
	script_thread = g_thread_new("script", execute_script, fp);

	g_free(script_file);
}

int initialize_script_menu() {
	GSList *list;
	GtkWidget *menuscript;
	gchar *home1_script;
	gchar *home2_script;
	gint i = 0, nb_item = 0;

#ifdef _WIN32
	wchar_t wFilename[MAX_PATH];

	home1_script = g_build_filename (get_special_folder (CSIDL_APPDATA),
			"siril", "scripts", NULL);
	home2_script = NULL;

	if (!GetModuleFileNameW(NULL, wFilename, MAX_PATH)) {
		fprintf(stderr, "initialize_script_menu error: %d\n", GetLastError());
	} else {
		gchar *fName = g_utf16_to_utf8(wFilename, -1, NULL, NULL, NULL);
		gchar *path = g_path_get_dirname(fName);
		path[strlen(path) - 4] = '\0';		/* remove "/bin" */
		home2_script = g_build_filename(path, "scripts", NULL);
		g_free(fName);
		g_free(path);
	}
#else
	home1_script = g_build_filename(g_get_home_dir(), ".siril", "scripts", NULL);
	home2_script = g_build_filename(g_get_home_dir(), "siril", "scripts", NULL);
#endif

	gchar *directories[] = { home1_script, home2_script, NULL };

	menuscript = lookup_widget("menuscript");
	GtkWidget *menu = gtk_menu_new();

	while (directories[i]) {
		list = search_script(directories[i]);
		if (list) {
			gtk_widget_show(menuscript);
			gtk_menu_item_set_submenu(GTK_MENU_ITEM(menuscript), menu);
			while (list) {
				nb_item ++;
				/* write separator but not for the first one */
				if (nb_item > 1) {
					GtkWidget *separator = gtk_separator_menu_item_new();
					gtk_menu_shell_append(GTK_MENU_SHELL(menu), separator);
					gtk_widget_show(separator);
				}
				/* write an item per script file */
				GtkWidget *menu_item;

				menu_item = gtk_menu_item_new_with_label(list->data);
				gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
				gchar *full_path = g_build_filename(directories[i], list->data,
				NULL);
				g_signal_connect(G_OBJECT(menu_item), "activate",
						G_CALLBACK(on_script_execution), (gchar * ) full_path);
				siril_log_message(_("Loading script: %s\n"), list->data);
				gtk_widget_show(menu_item);

				/* go to the next item */
				list = list->next;
			}
			g_slist_free_full(list, g_free);
		}
		i++;
	}
	g_free(home1_script);
	g_free(home2_script);
	return 0;
}
