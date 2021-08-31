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

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_app_dirs.h"

#include "siril_css.h"

#define CSS_FILENAME "siril.css"

static gchar *get_buffer_from_css_file(gchar *css) {
	GFile *file = g_file_new_for_path(css);
	GError *error = NULL;
	gchar *str;

	if (!g_file_load_contents(file, NULL, &str, NULL, NULL, &error)) {
		printf("Error loading %s: %s\n", g_file_peek_path(file), error->message);
		g_clear_error(&error);
	}
	return str;
}

/**
 * Loads the css sheet
 */
void load_css_style_sheet () {
	gchar *CSSFile;

	CSSFile = g_build_filename(siril_get_system_data_dir(), CSS_FILENAME, NULL);
	if (!g_file_test (CSSFile, G_FILE_TEST_EXISTS)) {
		g_error (_("Unable to load CSS style sheet file: %s. Please reinstall Siril\n"), CSSFile);
	}
	else {
		gchar *css_buffer = get_buffer_from_css_file(CSSFile);
		if (css_buffer) {
			/* make sure that scale is ok */
			if (com.pref.font_scale < 70.0) com.pref.font_scale = 100;

			GString *string = g_string_new(css_buffer);
			g_string_replace(string,
					"* { font-size: 1.0em; -gtk-icon-style: regular; }",
					"* { font-size: %lfem; -gtk-icon-style: %s; }", 1);

			gchar *new_buffer = g_string_free(string, FALSE);
			gchar *updated_css = g_strdup_printf(new_buffer, 1.0 + (com.pref.font_scale - 100.0) / 1000.0, com.pref.icon_symbolic ? "symbolic" : "regular");

			GtkCssProvider *css_provider = gtk_css_provider_new();
			GdkDisplay *display = gdk_display_get_default();
			GdkScreen *screen = gdk_display_get_default_screen(display);
			gtk_style_context_add_provider_for_screen(screen,
					GTK_STYLE_PROVIDER(css_provider),
					GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);

			gtk_css_provider_load_from_data(css_provider, updated_css, -1, NULL);

			g_fprintf(stdout, _("Successfully loaded '%s'\n"), CSSFile);
			g_free(new_buffer);
			g_free(css_buffer);
			g_free(updated_css);
			g_object_unref(css_provider);
		}
	}
	g_free(CSSFile);
}

