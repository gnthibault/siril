/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2017 team free-astro (see more in AUTHORS file)
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
#include "gui/utils.h"
#include "git-version.h"

#include "about_dialog.h"

static gchar **authors = (gchar *[] ) { "Vincent Hourdin <vh@free-astro.vinvin.tf>",
				"Cyril Richard <cyril@free-astro.org>", "François Meyer", NULL };

static gchar **documenters = (gchar *[] ) { "Laurent Rogé <l.roge@siril.org>", NULL };

static gchar **artists = (gchar *[] ) { "Maxime Oudoux <max.oudoux@gmail.com>", ""
				"Tobias Bernard <tbernard@gnome.org>",
				"Cyril Richard <cyril@free-astro.org>", NULL };

void siril_show_about_dialog() {
	GdkPixbuf *icon;
	GtkWindow *parent;
	gchar *copyright;
	gchar *version;
#ifdef SIRIL_UNSTABLE
	version = g_strdup_printf(_("%s\nThis is an unstable development release\n"
					"commit %s\n"), VERSION, SIRIL_GIT_VERSION_ABBREV);
#else
	version = g_strdup(VERSION);
#endif
	copyright = g_strdup_printf("Copyright © 2004-2011 François Meyer\n"
			"Copyright © 2012-%s team free-astro", SIRIL_GIT_LAST_COMMIT_YEAR);

	parent = GTK_WINDOW(GTK_APPLICATION_WINDOW(lookup_widget("control_window")));
	/* Create pixbuf from siril.svg file */
	gchar *image = g_build_filename(siril_get_system_data_dir(), "pixmaps", "siril.svg", NULL);
	icon = gdk_pixbuf_new_from_file_at_size(image, 256, 256, NULL);
	g_free(image);

	gtk_show_about_dialog(parent,
			"program-name", PACKAGE,
			"title", _("About Siril"),
			"logo", icon,
			"version", version,
			"copyright", copyright,
			"authors", authors,
			"documenters", documenters,
			"artists", artists,
			"comments", _("Astronomical image (pre-)processing program"),
			"translator-credits", _("translator-credits"),
			"website", PACKAGE_URL,
			"website-label", _("Visit the Siril website"),
			"license-type", GTK_LICENSE_GPL_3_0,
			NULL);

	g_free(copyright);
	g_free(version);
}
