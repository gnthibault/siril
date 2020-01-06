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
#include "gui/callbacks.h"
#include "core/siril_app_dirs.h"

#include "core/siril_cmd_help.h"

void siril_cmd_help_keyboard_shortcuts(GtkWindow *window) {
	static GtkWidget *shortcuts_window;
	char *shortcutfile = g_build_filename(siril_get_system_data_dir(), SHORTCUTS_UI, NULL);
	GError *err = NULL;

	if (shortcuts_window == NULL) {
		GtkBuilder *s_builder;

		s_builder = gtk_builder_new_from_file(shortcutfile);
		shortcuts_window = GTK_WIDGET(gtk_builder_get_object (s_builder, "shortcuts-siril"));

		g_signal_connect(shortcuts_window, "destroy",
				G_CALLBACK (gtk_widget_destroyed), &shortcuts_window);

		g_object_unref(s_builder);
	}

	if (GTK_WINDOW(window)
			!= gtk_window_get_transient_for(GTK_WINDOW(shortcuts_window))) {
		gtk_window_set_transient_for(GTK_WINDOW(shortcuts_window),
				GTK_WINDOW(window));
	}

	gtk_widget_show_all(shortcuts_window);
	gtk_window_present(GTK_WINDOW(shortcuts_window));
}
