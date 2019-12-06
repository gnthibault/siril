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

/*
 * GTK File Chooser static functions
 */

#include <stdio.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/initfile.h"
#include "core/exif.h"
#include "algos/sorting.h"
#include "io/conversion.h"
#include "io/films.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "callbacks.h"
#include "progress_and_log.h"
#include "message_dialog.h"

#include "open_dialog.h"

#define preview_size 256

static fileChooserPreview preview;

static void gtk_filter_add(GtkFileChooser *file_chooser, const gchar *title,
		const gchar *pattern, gboolean set_default) {
	gchar **patterns;
	gint i;

	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, title);
	/* get the patterns */
	patterns = g_strsplit(pattern, ";", -1);
	for (i = 0; patterns[i] != NULL; i++)
		gtk_file_filter_add_pattern(f, patterns[i]);
	/* free the patterns */
	g_strfreev(patterns);
	gtk_file_chooser_add_filter(file_chooser, f);
	if (set_default)
		gtk_file_chooser_set_filter(file_chooser, f);
}

static void set_filters_dialog(GtkFileChooser *chooser, int whichdial) {
	GString *all_filter = NULL;
	gchar *fits_filter = "*.fit;*.FIT;*.fits;*.FITS;*.fts;*.FTS";
	gchar *netpbm_filter = "*.ppm;*.PPM;*.pnm;*.PNM;*.pgm;*.PGM";
	gchar *pic_filter = "*.pic;*.PIC";
	gchar *ser_filter = "*.ser;*.SER";
	if (whichdial != OD_CONVERT) {
		gtk_filter_add(chooser, _("FITS Files (*.fit, *.fits, *.fts)"),
				fits_filter, com.filter == TYPEFITS);
	} else {
		all_filter = g_string_new(fits_filter);
	}
	if (whichdial == OD_OPEN || whichdial == OD_CONVERT) {
#ifdef HAVE_LIBRAW
		/* RAW FILES */
		int nb_raw;
		char *raw;
		int i;

		nb_raw = get_nb_raw_supported();
		raw = calloc(sizeof(char), nb_raw * 12 + 1);// we assume the extension size of 3 char "*.xxx;*.XXX;" = 12
		for (i = 0; i < nb_raw; i++) {
			char *ext;
			gchar *upcase;

			upcase = g_ascii_strup(supported_raw[i].extension, strlen(supported_raw[i].extension));
			ext = g_strdup_printf("*.%s;*.%s;", supported_raw[i].extension,
					upcase);
			strcat(raw, ext);

			g_free(ext);
			g_free(upcase);
		}
		raw[strlen(raw) - 1] = '\0';
		if (whichdial != OD_CONVERT) {
			gtk_filter_add(chooser, _("RAW DSLR Camera Files"), raw,
				com.filter == TYPERAW);
		} else {
			all_filter = g_string_append(all_filter, ";");
			all_filter = g_string_append(all_filter, raw);
		}

		free(raw);

#endif
		/* GRAPHICS FILES */
		GString *s_supported_graph, *s_pattern;
		gchar *graphics_supported, *graphics_filter;

		s_supported_graph = g_string_new(_("Graphics Files (*.bmp"));
		s_pattern = g_string_new("*.bmp;*.BMP;");
#ifdef HAVE_LIBJPEG
		s_supported_graph = g_string_append(s_supported_graph, ", *.jpg, *.jpeg");
		s_pattern = g_string_append(s_pattern, "*.jpg;*.JPG;*.jpeg;*.JPEG;");
#endif

#ifdef HAVE_LIBPNG
		s_supported_graph = g_string_append(s_supported_graph, ", *.png");
		s_pattern = g_string_append(s_pattern, "*.png;*.PNG;");
#endif

#ifdef HAVE_LIBTIFF
		s_supported_graph = g_string_append(s_supported_graph, ", *.tif, *.tiff");
		s_pattern = g_string_append(s_pattern, "*.tif;*.TIF;*.tiff;*.TIFF");
#endif
		s_supported_graph = g_string_append(s_supported_graph, ")");

		graphics_supported = g_string_free(s_supported_graph, FALSE);
		graphics_filter = g_string_free(s_pattern, FALSE);
		if (whichdial != OD_CONVERT) {
			gtk_filter_add(chooser, graphics_supported, graphics_filter,
					com.filter == TYPEBMP || com.filter == TYPEJPG
							|| com.filter == TYPEPNG || com.filter == TYPETIFF);

			/* NETPBM FILES */
			gtk_filter_add(chooser, _("Netpbm Files (*.ppm, *.pnm, *.pgm)"),
					netpbm_filter, com.filter == TYPEPNM);
			/* IRIS FILES */
			gtk_filter_add(chooser, _("IRIS PIC Files (*.pic)"), pic_filter,
					com.filter == TYPEPIC);
			/* SER FILES */
			gtk_filter_add(chooser, _("SER files (*.ser)"), ser_filter,
					com.filter == TYPESER);
		} else {
			all_filter = g_string_append(all_filter, ";");
			all_filter = g_string_append(all_filter, graphics_filter);
			all_filter = g_string_append(all_filter, ";");
			all_filter = g_string_append(all_filter, netpbm_filter);
			all_filter = g_string_append(all_filter, ";");
			all_filter = g_string_append(all_filter, pic_filter);
			all_filter = g_string_append(all_filter, ";");
			all_filter = g_string_append(all_filter, ser_filter);
		}

#ifdef HAVE_FFMS2
		/* FILM FILES */
		int nb_film;
		char *film_filter;
		int j;

		nb_film = get_nb_film_ext_supported();
		film_filter = calloc(sizeof(char), nb_film * 14 + 1);// we assume the extension size of 4 char "*.xxxx;*.XXXX;" = 14
		for (j = 0; j < nb_film; j++) {
			char *ext;
			gchar *upcase;

			upcase = g_ascii_strup(supported_film[j].extension,
					strlen(supported_film[j].extension));
			ext = g_strdup_printf("*.%s;*.%s;", supported_film[j].extension,
					upcase);
			strcat(film_filter, ext);

			g_free(ext);
			g_free(upcase);
		}
		film_filter[strlen(film_filter) - 1] = '\0';

		if (whichdial != OD_CONVERT) {
		gtk_filter_add(chooser, _("Film Files (*.avi, *.mpg, ...)"), film_filter,
				com.filter == TYPEAVI);
		} else {
			all_filter = g_string_append(all_filter, ";");
			all_filter = g_string_append(all_filter, film_filter);
		}
		free(film_filter);
#endif
		g_free(graphics_supported);
		g_free(graphics_filter);

		if (whichdial == OD_CONVERT) {
			gchar *filter = g_string_free(all_filter, FALSE);

			gtk_filter_add(chooser, _("All supported files"), filter, TRUE);
			g_free(filter);
		}
	}
}

struct _updta_preview_data {
	GtkFileChooser *file_chooser;
	gchar *filename;
	GdkPixbuf *pixbuf;
	GFileInfo *file_info;
	gboolean have_preview;
};

static gboolean end_update_preview_cb(gpointer p) {
	struct _updta_preview_data *args = (struct _updta_preview_data *) p;
	stop_processing_thread();// can it be done here in case there is no thread?

	if (args->have_preview) {
		int bytes;
		int width;
		int height;
		const char *bytes_str;
		char *size_str = NULL;
		char *dim_str = NULL;

		/* try to read file size */
		bytes_str = gdk_pixbuf_get_option(args->pixbuf, "tEXt::Thumb::Size");
		if (bytes_str != NULL) {
			bytes = atoi(bytes_str);
			size_str = g_format_size(bytes);
		} else {
			size_str = g_format_size(g_file_info_get_size(args->file_info));
		}

		/* try to read image dimensions */
		GdkPixbufFormat *pixbuf_file_info = gdk_pixbuf_get_file_info(args->filename,
				&width, &height);

		if (pixbuf_file_info != NULL) {
			/* Pixel size of image: width x height in pixel */
			dim_str = g_strdup_printf("%d x %d %s", width, height,
					ngettext("pixel", "pixels", height));
		}

		gtk_label_set_text(GTK_LABEL(preview.dim_label), dim_str);
		gtk_label_set_text(GTK_LABEL(preview.size_label), size_str);
		gtk_image_set_from_pixbuf(GTK_IMAGE(preview.image), args->pixbuf);

		g_free(dim_str);
		g_free(size_str);
		g_object_unref(args->pixbuf);
	}

	gtk_file_chooser_set_preview_widget_active(args->file_chooser, args->have_preview);
	g_object_unref(args->file_info);
	g_free(args->filename);
	free(args);
	return FALSE;
}

static gpointer update_preview_cb_idle(gpointer p) {
	uint8_t *buffer = NULL;
	size_t size;
	char *mime_type = NULL;
	GdkPixbuf *pixbuf = NULL;
	gboolean have_preview = FALSE;

	struct _updta_preview_data *args = (struct _updta_preview_data *) p;

	if (!siril_exif_get_thumbnail(args->filename, &buffer, &size, &mime_type)) {
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
		float ratio = 1.0 * gdk_pixbuf_get_height(tmp)
				/ gdk_pixbuf_get_width(tmp);
		int width = preview_size, height = preview_size * ratio;
		pixbuf = gdk_pixbuf_scale_simple(tmp, width, height,
				GDK_INTERP_BILINEAR);

		have_preview = TRUE;

		cleanup: gdk_pixbuf_loader_close(loader, NULL);
		free(mime_type);
		free(buffer);
		g_object_unref(loader); // This should clean up tmp as well
	}

	if (!have_preview) {
		pixbuf = gdk_pixbuf_new_from_file_at_size(args->filename, preview_size, preview_size, NULL);
	}
	if(pixbuf != NULL) have_preview = TRUE;

	args->have_preview = have_preview;
	args->pixbuf = pixbuf;
	siril_add_idle(end_update_preview_cb, args);
	return GINT_TO_POINTER(0);
}

static void update_preview_cb(GtkFileChooser *file_chooser, gpointer p) {
	gchar *filename;
	char *uri;
	GFile *file;
	GFileInfo *file_info;

	uri = gtk_file_chooser_get_preview_uri (file_chooser);
	if (uri == NULL) {
		gtk_file_chooser_set_preview_widget_active (file_chooser, FALSE);
		return;
	}

	file = g_file_new_for_uri (uri);
	file_info = g_file_query_info (file,
				       G_FILE_ATTRIBUTE_TIME_MODIFIED ","
				       G_FILE_ATTRIBUTE_STANDARD_TYPE ","
				       G_FILE_ATTRIBUTE_STANDARD_SIZE ","
				       G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
				       0, NULL, NULL);

	filename = gtk_file_chooser_get_preview_filename(file_chooser);

	struct _updta_preview_data *data = malloc(sizeof(struct _updta_preview_data));
	data->filename = filename;
	data->file_info = file_info;
	data->file_chooser = file_chooser;

	start_in_new_thread(update_preview_cb_idle, data);
}

static void siril_file_chooser_add_preview(GtkWidget *widget) {
	if (com.show_preview) {
		GtkWidget *vbox;

		vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
		gtk_container_set_border_width(GTK_CONTAINER(vbox), 12);

		preview.image = gtk_image_new();
		preview.dim_label = gtk_label_new(NULL);
		preview.size_label = gtk_label_new(NULL);

		gtk_box_pack_start(GTK_BOX(vbox), preview.image, FALSE, TRUE, 0);
		gtk_box_pack_start(GTK_BOX(vbox), preview.dim_label, FALSE, TRUE, 0);
		gtk_box_pack_start(GTK_BOX(vbox), preview.size_label, FALSE, TRUE, 0);

		gtk_widget_show_all(vbox);

		gtk_file_chooser_set_preview_widget(GTK_FILE_CHOOSER(widget), vbox);
		gtk_file_chooser_set_use_preview_label(GTK_FILE_CHOOSER(widget), FALSE);
		gtk_file_chooser_set_preview_widget_active(GTK_FILE_CHOOSER(widget), FALSE);

		g_signal_connect(widget, "update-preview", G_CALLBACK(update_preview_cb), NULL);
	}
}

static void opendial(int whichdial) {
	SirilWidget *widgetdialog;
	GtkFileChooser *dialog = NULL;
	GtkWindow *control_window = GTK_WINDOW(lookup_widget("control_window"));
	gint res;

	if (!com.wd)
		return;

	switch (whichdial) {
	case OD_NULL:
		fprintf(stderr, "whichdial undefined, should not happen\n");
		return;
	case OD_FLAT:
	case OD_DARK:
	case OD_OFFSET:
		widgetdialog = siril_file_chooser_open(control_window, GTK_FILE_CHOOSER_ACTION_OPEN);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		gtk_file_chooser_set_current_folder(dialog, com.wd);
		gtk_file_chooser_set_select_multiple(dialog, FALSE);
		set_filters_dialog(dialog, whichdial);
		break;
	case OD_CWD:
		widgetdialog = siril_file_chooser_open(control_window, GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		gtk_file_chooser_set_current_folder(dialog, com.wd);
		gtk_file_chooser_set_select_multiple(dialog, FALSE);
		break;
	case OD_OPEN:
		widgetdialog = siril_file_chooser_open(control_window, GTK_FILE_CHOOSER_ACTION_OPEN);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		gtk_file_chooser_set_current_folder(dialog, com.wd);
		gtk_file_chooser_set_select_multiple(dialog, FALSE);
		set_filters_dialog(dialog, whichdial);
		siril_file_chooser_add_preview(widgetdialog);
		break;
	case OD_CONVERT:
		widgetdialog = siril_file_chooser_add(control_window, GTK_FILE_CHOOSER_ACTION_OPEN);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		gtk_file_chooser_set_current_folder(dialog, com.wd);
		gtk_file_chooser_set_select_multiple(dialog, TRUE);
		set_filters_dialog(dialog, whichdial);
	}

	if (!dialog)
		return;

	res = siril_dialog_run(widgetdialog);

	if (res == GTK_RESPONSE_ACCEPT) {
		GSList *list = NULL;
		gchar *filename, *err;
		gboolean anything_loaded;
		GtkFileChooser *chooser = GTK_FILE_CHOOSER(dialog);
		GtkEntry *flat_entry, *dark_entry, *bias_entry;
		GtkToggleButton *flat_button, *dark_button, *bias_button;
		GtkWidget *pbutton;

		filename = gtk_file_chooser_get_filename(chooser);
		if (!filename)
			return;

		pbutton  = lookup_widget("prepro_button");
		flat_entry = GTK_ENTRY(lookup_widget("flatname_entry"));
		dark_entry = GTK_ENTRY(lookup_widget("darkname_entry"));
		bias_entry = GTK_ENTRY(lookup_widget("offsetname_entry"));

		flat_button = GTK_TOGGLE_BUTTON(lookup_widget("useflat_button"));
		dark_button = GTK_TOGGLE_BUTTON(lookup_widget("usedark_button"));
		bias_button = GTK_TOGGLE_BUTTON(lookup_widget("useoffset_button"));

		anything_loaded = sequence_is_loaded() || single_image_is_loaded();

		switch (whichdial) {
		case OD_FLAT:
			gtk_entry_set_text(flat_entry, filename);
			gtk_toggle_button_set_active(flat_button, TRUE);
			gtk_widget_set_sensitive(pbutton, anything_loaded);
			break;

		case OD_DARK:
			gtk_entry_set_text(dark_entry, filename);
			gtk_toggle_button_set_active(dark_button, TRUE);
			gtk_widget_set_sensitive(pbutton, anything_loaded);
			break;

		case OD_OFFSET:
			gtk_entry_set_text(bias_entry, filename);
			gtk_toggle_button_set_active(bias_button, TRUE);
			gtk_widget_set_sensitive(pbutton, anything_loaded);
			break;

		case OD_CWD:
			if (!changedir(filename, &err)) {
				writeinitfile();
				set_GUI_CWD();
			} else {
				siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), err);
			}
			break;

		case OD_OPEN:
			set_cursor_waiting(TRUE);
			open_single_image(filename);
			set_cursor_waiting(FALSE);
			break;

		case OD_CONVERT:
			list = gtk_file_chooser_get_filenames(chooser);
			list = g_slist_sort(list, (GCompareFunc) strcompare);
			fill_convert_list(list);
			g_slist_free(list);
			break;
		}
		g_free(filename);
	}
	siril_widget_destroy(widgetdialog);
}

void on_darkfile_button_clicked(GtkButton *button, gpointer user_data) {
	opendial(OD_DARK);
}

void cwd_btton_clicked() {
	opendial(OD_CWD);
}

void on_offsetfile_button_clicked(GtkButton *button, gpointer user_data) {
	opendial(OD_OFFSET);
}

void on_flatfile_button_clicked(GtkButton *button, gpointer user_data) {
	opendial(OD_FLAT);
}

void header_open_button_clicked() {
	opendial(OD_OPEN);
}

void on_select_convert_button_clicked(GtkToolButton *button, gpointer user_data) {
	opendial(OD_CONVERT);
}

/** callback function for recent document opened
 *
 * @param chooser
 * @param user_data
 */
void on_open_recent_action_item_activated(GtkRecentChooser *chooser,
		gpointer user_data) {
	gchar *uri, *path;
	GError *error = NULL;

	uri = gtk_recent_chooser_get_current_uri(chooser);

	path = g_filename_from_uri(uri, NULL, &error);
	if (error) {
		g_warning("Could not convert uri \"%s\" to a local path: %s", uri,
				error->message);
		g_error_free(error);
		return;
	}

	open_single_image(path);

	g_free(uri);
	g_free(path);
}
