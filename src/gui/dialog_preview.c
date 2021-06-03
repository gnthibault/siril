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
#include "core/processing.h"
#include "core/exif.h"
#include "gui/utils.h"
#include "gui/histogram.h"
#include "io/ser.h"
#include "io/image_format_fits.h"

#ifdef HAVE_LIBHEIF
#include <libheif/heif.h>
#endif

#include "dialog_preview.h"

static gboolean preview_allocated = FALSE; // flag needed when user load image before prevew was displayed.

struct _updta_preview_data {
	GtkFileChooser *file_chooser;
	gchar *filename;
	gchar *description;
	GdkPixbuf *pixbuf;
	GFileInfo *file_info;
	fileChooserPreview *preview;
};

struct _fileChooserPreview {
	GtkWidget *image;
	GtkWidget *name_label;
	GtkWidget *dim_label;
	GtkWidget *size_label;
};

static fileChooserPreview *new_preview_object() {
	fileChooserPreview *object = g_new(fileChooserPreview, 1);
	object->image = gtk_image_new();
	object->name_label = gtk_label_new(NULL);
	object->dim_label = gtk_label_new(NULL);
	object->size_label = gtk_label_new(NULL);
	return object;
}

static gboolean end_update_preview_cb(gpointer p) {
	struct _updta_preview_data *args = (struct _updta_preview_data *) p;

	stop_processing_thread();

	int bytes;
	const char *bytes_str;
	char *size_str = NULL;
	char *name_str = NULL;
	char *info_str = NULL;
	GFileType type;

	fileChooserPreview *preview = args->preview;

	if (!preview_allocated || !preview || !(GTK_IS_IMAGE((preview->image)))) {
		return FALSE;
	}

	name_str = g_path_get_basename(args->filename);

	if (!args->file_info) {
		return FALSE;
	}
	type = g_file_info_get_file_type(args->file_info);

	/* try to read file size */
	if (args->pixbuf && (bytes_str = gdk_pixbuf_get_option(args->pixbuf, "tEXt::Thumb::Size")) != NULL) {
		bytes = g_ascii_strtoll(bytes_str, NULL, 10);
		size_str = g_format_size(bytes);
	} else {
		if (type == G_FILE_TYPE_REGULAR) {
			size_str = g_format_size(g_file_info_get_size(args->file_info));
		}
	}

	/* load icon */
	if (type == G_FILE_TYPE_REGULAR && args->pixbuf) {
		gtk_image_set_from_pixbuf(GTK_IMAGE(preview->image), args->pixbuf);
		info_str = args->description;
	} else if (type == G_FILE_TYPE_DIRECTORY) {
		gtk_image_set_from_icon_name(GTK_IMAGE(preview->image), "folder", GTK_ICON_SIZE_DIALOG);
		gtk_image_set_pixel_size(GTK_IMAGE(preview->image), com.pref.thumbnail_size);
		info_str = g_strdup(_("Folder"));
	} else {
		image_type im_type = get_type_from_filename(args->filename);
		if (im_type == TYPEAVI || im_type == TYPESER ||
				(im_type == TYPEFITS && fitseq_is_fitseq(args->filename, NULL)))
			gtk_image_set_from_icon_name(GTK_IMAGE(preview->image), "video", GTK_ICON_SIZE_DIALOG);
		else gtk_image_set_from_icon_name(GTK_IMAGE(preview->image), "image", GTK_ICON_SIZE_DIALOG);
		gtk_image_set_pixel_size(GTK_IMAGE(preview->image), com.pref.thumbnail_size);
	}

	/* information strings */
	const char *format = "<span style=\"italic\">%s</span>";
	char *markup = g_markup_printf_escaped(format, name_str);
	gtk_label_set_markup(GTK_LABEL(preview->name_label), markup);
	gtk_label_set_ellipsize(GTK_LABEL(preview->name_label), PANGO_ELLIPSIZE_MIDDLE);
	gtk_label_set_width_chars(GTK_LABEL(preview->name_label), 25);
	gtk_label_set_max_width_chars(GTK_LABEL(preview->name_label), 25);

	gtk_label_set_text(GTK_LABEL(preview->dim_label), info_str);
	gtk_label_set_text(GTK_LABEL(preview->size_label), size_str);

	if (args->pixbuf)
		g_object_unref(args->pixbuf);
	g_free(markup);
	g_free(name_str);
	g_free(info_str);
	g_free(size_str);

	g_object_unref(args->file_info);
	g_free(args->filename);

	free(args);
	args = NULL;
	return FALSE;
}

static gpointer update_preview_cb_idle(gpointer p) {
	uint8_t *buffer = NULL;
	size_t size;
	char *mime_type = NULL;
	GdkPixbuf *pixbuf = NULL;
	image_type im_type;
	gboolean libheif_is_ok = FALSE;

	struct _updta_preview_data *args = (struct _updta_preview_data *) p;

	args->description = NULL;

	im_type = get_type_from_filename(args->filename);

	if (im_type == TYPEFITS) {
		/* try FITS file */
		pixbuf = get_thumbnail_from_fits(args->filename, &args->description);
	} else if (im_type == TYPESER) {
		pixbuf = get_thumbnail_from_ser(args->filename, &args->description);
	} else {
		if (im_type != TYPEUNDEF && !siril_get_thumbnail_exiv(args->filename, &buffer, &size,
				&mime_type)) {
			// Scale the image to the correct size
			GdkPixbuf *tmp;
			GdkPixbufLoader *loader = gdk_pixbuf_loader_new();
			if (!gdk_pixbuf_loader_write(loader, buffer, size, NULL))
				goto cleanup;
			// Calling gdk_pixbuf_loader_close forces the data to be parsed by the
			// loader. We must do this before calling gdk_pixbuf_loader_get_pixbuf.
			if (!gdk_pixbuf_loader_close(loader, NULL))
				goto cleanup;
			if (!(tmp = gdk_pixbuf_loader_get_pixbuf(loader)))
				goto cleanup;
			float ratio = 1.0 * gdk_pixbuf_get_height(tmp) / gdk_pixbuf_get_width(tmp);
			int width = com.pref.thumbnail_size, height = com.pref.thumbnail_size * ratio;
			pixbuf = gdk_pixbuf_scale_simple(tmp, width, height, GDK_INTERP_BILINEAR);
			args->description = siril_get_file_info(args->filename, pixbuf);

			cleanup: gdk_pixbuf_loader_close(loader, NULL);
			free(mime_type);
			free(buffer);
			g_object_unref(loader); // This should clean up tmp as well
		}

		/* if no pixbuf created try to directly read the file */
		/* libheif < 1.6.2 has a bug, therefore we can't open preview if libheif is too old
		 * bug fixed in https://github.com/strukturag/libheif/commit/fbd6d28e8604ecb53a2eb33b522a664b6bcabd0b*/
#ifdef HAVE_LIBHEIF
		libheif_is_ok = LIBHEIF_HAVE_VERSION(1, 6, 2);
#endif

		if (!pixbuf && (im_type != TYPEHEIF || libheif_is_ok)) {
			pixbuf = gdk_pixbuf_new_from_file_at_size(args->filename,
					com.pref.thumbnail_size, com.pref.thumbnail_size, NULL);
			args->description = siril_get_file_info(args->filename, pixbuf);
		}
	}

	args->pixbuf = pixbuf;
	siril_add_idle(end_update_preview_cb, args);
	return GINT_TO_POINTER(0);
}

static void update_preview_cb(GtkFileChooser *file_chooser, gpointer p) {
	gchar *uri;
	GFile *file;
	GFileInfo *file_info;
	fileChooserPreview *preview = (fileChooserPreview *)p;

	uri = gtk_file_chooser_get_preview_uri(file_chooser);
	if (uri == NULL) {
		gtk_file_chooser_set_preview_widget_active(file_chooser, FALSE);
		return;
	}

	file = g_file_new_for_uri(uri);
	file_info = g_file_query_info(file,
				       G_FILE_ATTRIBUTE_TIME_MODIFIED ","
				       G_FILE_ATTRIBUTE_STANDARD_TYPE ","
				       G_FILE_ATTRIBUTE_STANDARD_SIZE ","
				       G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
				       0, NULL, NULL);

	gtk_file_chooser_set_preview_widget_active(file_chooser, TRUE);

	struct _updta_preview_data *data = malloc(sizeof(struct _updta_preview_data));
	data->filename = g_file_get_path(file);
	data->file_info = file_info;
	data->file_chooser = file_chooser;
	data->preview = preview;

	g_free(uri);
	g_object_unref(file);

	start_in_new_thread(update_preview_cb_idle, data);
}

void siril_preview_free(fileChooserPreview *preview) {
	g_free(preview);
	preview_allocated = FALSE;
}

void siril_file_chooser_add_preview(GtkFileChooser *dialog, fileChooserPreview *preview) {
	if (com.pref.show_thumbnails) {
		GtkWidget *vbox;

		vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
		gtk_container_set_border_width(GTK_CONTAINER(vbox), 12);

		preview = new_preview_object();
		preview_allocated = TRUE;

		gtk_label_set_justify(GTK_LABEL(preview->name_label), GTK_JUSTIFY_CENTER);
		gtk_label_set_justify(GTK_LABEL(preview->dim_label), GTK_JUSTIFY_CENTER);
		gtk_label_set_justify(GTK_LABEL(preview->dim_label), GTK_JUSTIFY_CENTER);

		gtk_widget_set_size_request(preview->image, com.pref.thumbnail_size, com.pref.thumbnail_size);

		gtk_box_pack_start(GTK_BOX(vbox), preview->image, FALSE, TRUE, 0);
		gtk_box_pack_start(GTK_BOX(vbox), preview->name_label, FALSE, TRUE, 10);
		gtk_box_pack_start(GTK_BOX(vbox), preview->size_label, FALSE, TRUE, 0);
		gtk_box_pack_start(GTK_BOX(vbox), preview->dim_label, FALSE, TRUE, 0);

		gtk_widget_show_all(vbox);

		gtk_file_chooser_set_preview_widget(dialog, vbox);
		gtk_file_chooser_set_use_preview_label(dialog, FALSE);
		gtk_file_chooser_set_preview_widget_active(dialog, FALSE);

		g_signal_connect(dialog, "update-preview", G_CALLBACK(update_preview_cb), (gpointer)preview);
	}
}
