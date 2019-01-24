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
#include "gui/callbacks.h"

static const gchar* copyright = ("Copyright © 2004-2011 François Meyer\n"
		"Copyright © 2012-2019 team free-astro");

static gchar **authors = (gchar *[] ) { "Vincent Hourdin <vh@free-astro.vinvin.tf>",
				"Cyril Richard <cyril@free-astro.org>", "François Meyer", NULL };

static gchar **documenters = (gchar *[] ) { "Laurent Rogé <l.roge@siril.org>", NULL };

static gchar **artists = (gchar *[] ) { "Maxime Oudoux <max.oudoux@gmail.com>",
				"Cyril Richard <cyril@free-astro.org>", NULL };

// translator names
static gchar *translator = N_("Cyril Richard <cyril@free-astro.org>\n"
		"Vincent Hourdin <vh@free-astro.vinvin.tf>");

void on_about_activate(GtkMenuItem *menuitem, gpointer user_data) {
	GdkPixbuf *icon;
	GtkWindow *parent;

	parent = GTK_WINDOW(lookup_widget("control_window"));
	icon = gtk_image_get_pixbuf(GTK_IMAGE(lookup_widget("pixmap1")));
	gtk_show_about_dialog(parent,
			"program-name", PACKAGE,
			"title", _("About Siril"),
			"logo", icon,
			"version", VERSION,
			"copyright", copyright,
			"authors", authors,
			"documenters", documenters,
			"artists", artists,
			"comments", _("Astronomical image (pre-)processing program"),
			"translator-credits", _(translator),
			"website", "https://siril.org",
			"website-label", _("Visit the Siril website"),
			"license-type", GTK_LICENSE_GPL_3_0,
			NULL);
}
