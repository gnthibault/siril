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

#include "utils.h"

struct _label_data {
	const char *label_name;
	char *text;
};

static gboolean set_label_text_idle(gpointer p) {
	struct _label_data *args = (struct _label_data *) p;
	GtkLabel *label = GTK_LABEL(lookup_widget(args->label_name));

	gtk_label_set_text(label, args->text);
	free(args->text);
	free(args);
	return FALSE;
}

void set_label_text_from_main_thread(const char *label_name, const char *text) {
	struct _label_data *data = malloc(sizeof(struct _label_data));
	data->label_name = label_name;
	data->text = strdup(text);
	gdk_threads_add_idle(set_label_text_idle, data);
}

GtkWidget* lookup_widget(const gchar *widget_name) {
	return GTK_WIDGET(gtk_builder_get_object(builder, widget_name));
}

void control_window_switch_to_tab(main_tabs tab) {
	GtkNotebook* notebook = GTK_NOTEBOOK(lookup_widget("notebook_center_box"));
	gtk_notebook_set_current_page(notebook, tab);
}

/**
 * Create a popover with icon and text
 * @param widget is the parent widget where the popover arises from
 * @param text will be shown in the popover
 * @return the GtkWidget of popover
 */
GtkWidget* popover_new(GtkWidget *widget, const gchar *text) {
	return popover_new_with_image(widget, text, NULL);
}

/**
 * Create a popover with icon and text
 * @param widget is the parent widget where the popover arises from
 * @param text will be shown in the popover
 * @param pixbuf will be shown in the popover
 * @return the GtkWidget of popover
 */
GtkWidget* popover_new_with_image(GtkWidget *widget, const gchar *text, GdkPixbuf *pixbuf) {
	GtkWidget *popover, *box, *image, *label;

	popover = gtk_popover_new(widget);
	label = gtk_label_new(NULL);
	box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);

	if (pixbuf) {
		float ratio = 1.0 * gdk_pixbuf_get_height(pixbuf) / gdk_pixbuf_get_width(pixbuf);
		int width = 128, height = 128 * ratio;
		GdkPixbuf *new_pixbuf = gdk_pixbuf_scale_simple(pixbuf, width, height,
				GDK_INTERP_BILINEAR);
		image = gtk_image_new_from_pixbuf(new_pixbuf);
		g_object_unref(new_pixbuf);
	} else {
		image = gtk_image_new_from_icon_name("dialog-information-symbolic",
					GTK_ICON_SIZE_DIALOG);
	}

	gtk_label_set_markup(GTK_LABEL(label), text);
	gtk_label_set_line_wrap(GTK_LABEL(label), TRUE);
	gtk_label_set_max_width_chars(GTK_LABEL(label), 64);

	gtk_box_pack_start(GTK_BOX(box), image, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(box), label, FALSE, FALSE, 0);
	gtk_container_add(GTK_CONTAINER(popover), box);

	/* make all sensitive in case where parent is not */
	gtk_widget_set_sensitive(label, TRUE);
	gtk_widget_set_sensitive(box, TRUE);
	gtk_widget_set_sensitive(popover, TRUE);

	gtk_widget_show_all(box);

	return popover;
}

/* size is in Bytes */
void set_GUI_MEM(guint64 used, const gchar *label) {
	if (com.headless)
		return;
	gchar *str;
	if (used > 0) {
		gchar *mem = g_format_size_full(used, G_FORMAT_SIZE_IEC_UNITS);
		str = g_strdup_printf(_("Mem: %s"), mem);
		g_free(mem);
	} else {
		str = g_strdup(_("Mem: N/A"));
	}
	set_label_text_from_main_thread(label, str);
	g_free(str);
}

void set_GUI_DiskSpace(gint64 space, const gchar *label) {
	if (com.headless)
		return;
	gchar *str;
	GtkStyleContext *context = gtk_widget_get_style_context(lookup_widget(label));
	gtk_style_context_remove_class(context, "label-info");

	if (space > 0) {
		if (space < 1073741824) { // we want to warn user of space is less than 1GiB
			gtk_style_context_add_class(context, "label-info");
		}
		gchar *mem = g_format_size_full(space, G_FORMAT_SIZE_IEC_UNITS);
		str = g_strdup_printf(_("Disk Space: %s"), mem);
		g_free(mem);
	} else {
		str = g_strdup(_("Disk Space: N/A"));
	}
	set_label_text_from_main_thread(label, str);
	g_free(str);
}

void set_suggested(GtkWidget *widget) {
	gtk_style_context_add_class(gtk_widget_get_style_context(widget),
			GTK_STYLE_CLASS_SUGGESTED_ACTION);
}

void unset_suggested(GtkWidget *widget) {
	gtk_style_context_remove_class(gtk_widget_get_style_context(widget),
			GTK_STYLE_CLASS_SUGGESTED_ACTION);
}
