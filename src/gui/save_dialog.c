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

#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "gui/progress_and_log.h"
#include "gui/callbacks.h"
#include "gui/save_dialog.h"
#include "gui/message_dialog.h"
#include "io/sequence.h"
#include "io/single_image.h"


static image_type whichminisave = TYPEUNDEF;
static SirilWidget *saveDialog = NULL;

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

static void set_filters_save_dialog(GtkFileChooser *chooser) {
	gtk_filter_add(chooser, _("FITS Files (*.fit, *.fits, *.fts)"),
			"*.fit;*.FIT;*.fits;*.FITS;*.fts;*.FTS", com.filter == TYPEFITS);

	/* GRAPHICS FILES */
	/* BMP FILES */
	gtk_filter_add(chooser, _("BMP Files (*.bmp)"), "*.bmp;*.BMP",
			com.filter == TYPEBMP);
#ifdef HAVE_LIBJPEG
	gtk_filter_add(chooser, _("JPEG Files (*.jpg, *.jpeg)"),
			"*.jpg;*.JPG;*.jpeg;*.JPEG", com.filter == TYPEJPG);
#endif

#ifdef HAVE_LIBPNG
	gtk_filter_add(chooser, _("PNG Files (*.png)"), "*.png;*.PNG",
			com.filter == TYPEPNG);
#endif

#ifdef HAVE_LIBTIFF
	gtk_filter_add(chooser, _("TIFF Files (*.tif, *.tiff)"),
			"*.tif;*.TIF;*.tiff;*.TIFF", com.filter == TYPETIFF);
#endif

	/* NETPBM FILES */
	gtk_filter_add(chooser, _("Netpbm Files (*.ppm, *.pnm, *.pgm)"),
			"*.ppm;*.PPM;*.pnm;*.PNM;*.pgm;*.PGM", com.filter == TYPEPNM);
}

static int get_filetype(const gchar *filter) {
	gchar **string;
	int type = TYPEUNDEF;
	int i = 0;

	string = g_strsplit_set(filter, "*(),.", -1);

	while (string[i]) {
		if (!g_strcmp0(string[i], "fit")) {
			type = TYPEFITS;
			break;
		} else if (!g_strcmp0(string[i], "bmp")) {
			type = TYPEBMP;
			break;
		} else if (!g_strcmp0(string[i], "jpg")) {
			type = TYPEJPG;
			break;
		} else if (!g_strcmp0(string[i], "png")) {
			type = TYPEPNG;
			break;
		} else if (!g_strcmp0(string[i], "tif")) {
			type = TYPETIFF;
			break;
		} else if (!g_strcmp0(string[i], "ppm")) {
			type = TYPEPNM;
			break;
		}
		i++;
	}
	g_strfreev(string);

	return type;
}

static void set_programm_name_in_TIFF() {
	static GtkTextView *TIFF_txt = NULL;
	GtkTextBuffer *tbuf;
	GtkTextIter itStart, itEnd;
	gchar *copyright;

	if (TIFF_txt == NULL)
		TIFF_txt = GTK_TEXT_VIEW(lookup_widget("Copyright_txt"));

	tbuf = gtk_text_view_get_buffer(TIFF_txt);

	copyright = g_strdup_printf("%s v%s", PACKAGE, VERSION);
	copyright[0] = toupper(copyright[0]);			// convert siril to Siril

	gtk_text_buffer_get_bounds(tbuf, &itStart, &itEnd);
	gtk_text_buffer_delete(tbuf, &itStart, &itEnd);
	gtk_text_buffer_set_text(tbuf, copyright, strlen(copyright));

	g_free(copyright);
}

static void set_description_in_TIFF() {
	static GtkTextView *TIFF_txt = NULL;
	GtkTextBuffer *tbuf;
	GtkTextIter itStart, itEnd;
	int i;

	if (TIFF_txt == NULL)
		TIFF_txt = GTK_TEXT_VIEW(lookup_widget("Description_txt"));

	tbuf = gtk_text_view_get_buffer(TIFF_txt);

	gtk_text_buffer_get_bounds(tbuf, &itStart, &itEnd);
	gtk_text_buffer_delete(tbuf, &itStart, &itEnd);
	/* History already written in header */
	if (gfit.history) {
		GSList *list;
		for (list = gfit.history; list; list = list->next) {
			gtk_text_buffer_get_end_iter(tbuf, &itEnd);
			gtk_text_buffer_insert(tbuf, &itEnd, (gchar *)list->data, -1);
			gtk_text_buffer_get_end_iter(tbuf, &itEnd);
			gtk_text_buffer_insert(tbuf, &itEnd, "\n", 1);
		}
	}
	/* New history */
	if (com.history) {
		for (i = 0; i < com.hist_display; i++) {
			if (com.history[i].history[0] != '\0') {
				gtk_text_buffer_get_end_iter(tbuf, &itEnd);
				gtk_text_buffer_insert(tbuf, &itEnd, com.history[i].history, strlen(com.history[i].history));
				gtk_text_buffer_get_end_iter(tbuf, &itEnd);
				gtk_text_buffer_insert(tbuf, &itEnd, "\n", 1);
			}
		}
	}
}

static void prepare_savepopup(int type) {
	static GtkNotebook* notebookFormat = NULL;
	static GtkWidget *savepopup = NULL;
	static GtkWidget *savetxt = NULL;
	GtkWindow *parent;
	int tab;

	if (notebookFormat == NULL) {
		notebookFormat = GTK_NOTEBOOK(lookup_widget("notebookFormat"));
		savepopup = lookup_widget("savepopup");
		savetxt = lookup_widget("filenameframe");
	}
	parent = siril_get_active_window();
	if (!GTK_IS_WINDOW(parent)) {
		parent = GTK_WINDOW(lookup_widget("control_window"));
	}
	gtk_window_set_transient_for(GTK_WINDOW(savepopup),	parent);

	switch (type) {
	case TYPEBMP:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving BMP"));
		tab = PAGE_MISC;
		break;
	case TYPEPNG:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving PNG"));
		tab = PAGE_MISC;
		break;
	case TYPEPNM:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving Netpbm"));
		tab = PAGE_MISC;
		break;
	case TYPEJPG:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving JPG"));
		tab = PAGE_JPG;
		break;
	case TYPETIFF:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving TIFF"));
		set_programm_name_in_TIFF();
		set_description_in_TIFF();
		tab = PAGE_TIFF;
		break;
	default:
	case TYPEFITS:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving FITS"));
		tab = PAGE_FITS;
	}
	whichminisave = type;

	gtk_widget_set_visible(savetxt, FALSE);
	gtk_notebook_set_current_page(notebookFormat, tab);
}

static void init_dialog() {
	GtkWindow *parent;
	if (saveDialog == NULL) {
		parent = siril_get_active_window();
		saveDialog = siril_file_chooser_save(parent, GTK_FILE_CHOOSER_ACTION_SAVE);
	}
}

static void close_dialog() {
	if (saveDialog != NULL) {
		siril_widget_destroy(saveDialog);
		saveDialog = NULL;
	}
}

static gchar *get_filename() {
	gchar *basename;

	if (sequence_is_loaded() && com.seq.current != RESULT_IMAGE) {
		char fname[256];
		/* set the output file name default as the current image.jpg */
		seq_get_image_filename(&com.seq, com.seq.current, fname);
		basename = g_path_get_basename(fname);
	} else {
		basename = g_path_get_basename(com.uniq->filename);
	}

	return remove_ext_from_filename(basename);
}

static int save_dialog() {
	GtkFileChooser *chooser;
	GtkFileFilter *filter;
	GtkEntry *savetext;
	gint res;
	gchar *fname;

	init_dialog();

	chooser = GTK_FILE_CHOOSER(saveDialog);
	fname = get_filename();

	gtk_file_chooser_set_current_name(chooser, fname);
	gtk_file_chooser_set_do_overwrite_confirmation(chooser, TRUE);
	set_filters_save_dialog(chooser);
	g_free(fname);

	res = siril_dialog_run(saveDialog);
	if (res == GTK_RESPONSE_ACCEPT) {

		gchar *filename = gtk_file_chooser_get_filename(chooser);
		savetext = GTK_ENTRY(lookup_widget("savetxt"));
		gtk_entry_set_text(savetext, filename);
		g_free(filename);

		/* we get the type of filter */
		filter = gtk_file_chooser_get_filter(chooser);
		const gchar *str = gtk_file_filter_get_name(filter);
		prepare_savepopup(get_filetype(str));
		return res;
	}
	close_dialog();

	return res;
}

// idle function executed at the end of the Save Data processing
gboolean end_save(gpointer p) {
	struct savedial_data *args = (struct savedial_data *) p;
	if (args->retval)
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("File saving failed. Check the logs for more info."));

	gtk_entry_set_text(args->entry, "");
	gtk_widget_hide(lookup_widget("savepopup"));
	stop_processing_thread();
	set_cursor_waiting(FALSE);
	close_dialog();	// is this different from the hide above?
	update_MenuItem();
	update_used_memory();
	free(args);
	return FALSE;
}

static void initialize_data(gpointer p) {
	struct savedial_data *args = (struct savedial_data *) p;

	GtkToggleButton *fits_8 = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_save_fit8"));
	GtkToggleButton *fits_16s = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_save_fit16s"));
	GtkToggleButton *update_hilo = (GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_update_hilo")));
#ifdef HAVE_LIBJPEG
	GtkSpinButton *qlty_spin_button = GTK_SPIN_BUTTON(lookup_widget("quality_spinbutton"));
	args->quality = gtk_spin_button_get_value_as_int(qlty_spin_button);
#endif
#ifdef HAVE_LIBTIFF
	GtkToggleButton *BPS_Button = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton8bits"));
	args->bitspersamples = gtk_toggle_button_get_active(BPS_Button) ? 8 : 16;
#endif
	args->entry = GTK_ENTRY(lookup_widget("savetxt"));
	args->filename = gtk_entry_get_text(args->entry);

	if (gtk_toggle_button_get_active(fits_8))
		args->bitpix = BYTE_IMG;
	else if (gtk_toggle_button_get_active(fits_16s))
		args->bitpix = SHORT_IMG;
	else
		args->bitpix = USHORT_IMG;

	args->update_hilo = gtk_toggle_button_get_active(update_hilo);
}

static gpointer minisavedial(gpointer p) {
	struct savedial_data *args = (struct savedial_data *) p;
#ifdef HAVE_LIBPNG
	uint32_t bytes_per_sample;
#endif
	args->retval = 0;

	if (args->filename[0] != '\0') {
		switch (whichminisave) {
		case TYPEBMP:
			args->retval = savebmp(args->filename, &gfit);
			break;
#ifdef HAVE_LIBJPEG
		case TYPEJPG:
			args->retval = savejpg(args->filename, &gfit, args->quality);
			break;
#endif
#ifdef HAVE_LIBTIFF
		case TYPETIFF:
			args->retval = savetif(args->filename, &gfit, args->bitspersamples);
			break;
#endif
#ifdef HAVE_LIBPNG
		case TYPEPNG:
			bytes_per_sample = gfit.orig_bitpix != BYTE_IMG ? 2 : 1;
			args->retval = savepng(args->filename, &gfit, bytes_per_sample, gfit.naxes[2] == 3);
			break;
#endif
		case TYPEFITS:
			gfit.bitpix = args->bitpix;
			/* Check if MIPS-HI and MIPS-LO must be updated. If yes,
			 * Values are taken from the layer 0 */
			if (args->update_hilo) {
				if (sequence_is_loaded() && !single_image_is_loaded()) {
					gfit.hi = com.seq.layers[RLAYER].hi;
					gfit.lo = com.seq.layers[RLAYER].lo;
				} else {
					gfit.hi = com.uniq->layers[RLAYER].hi;
					gfit.lo = com.uniq->layers[RLAYER].lo;
				}
				if (gfit.orig_bitpix == BYTE_IMG
						&& (gfit.hi > UCHAR_MAX || gfit.lo > UCHAR_MAX)) {
					gfit.hi = UCHAR_MAX;
					gfit.lo = 0;
				} else if (gfit.orig_bitpix == SHORT_IMG
						&& (gfit.hi > SHRT_MAX || gfit.lo > SHRT_MAX)) {
					gfit.hi = UCHAR_MAX;
					gfit.lo = 0;
				}
				if (gfit.orig_bitpix == BYTE_IMG && gfit.bitpix != BYTE_IMG) {
					gfit.hi = USHRT_MAX;
					gfit.lo = 0;
				}
			}
			args->retval = savefits(args->filename, &gfit);
			break;
		case TYPEPNM:
			args->retval = saveNetPBM(args->filename, &gfit);
			break;
		default:
			siril_log_message(_("This type of file is not handled. Should not happen"));
			break;
		}
	}
	siril_add_idle(end_save, args);
	return NULL;
}

void on_menu_rgb_savefits_activate(GtkMenuItem *menuitem, gpointer user_data) {
	static GtkNotebook* notebookFormat = NULL;
	static GtkWidget *savepopup = NULL;
	GtkWidget *savetxt = lookup_widget("filenameframe");
	GtkToggleButton *b8bit = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_save_fit8"));
	GtkToggleButton *b16bitu = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_save_fit16"));
	GtkToggleButton *b16bits = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_save_fit16s"));

	if (notebookFormat == NULL) {
		notebookFormat = GTK_NOTEBOOK(
				gtk_builder_get_object(builder, "notebookFormat"));
		savepopup = lookup_widget("savepopup");
	}
	if (single_image_is_loaded() || sequence_is_loaded()) {
		switch(gfit.bitpix) {
		case BYTE_IMG:
			gtk_toggle_button_set_active(b8bit, TRUE);
			break;
		case SHORT_IMG:
			gtk_toggle_button_set_active(b16bits, TRUE);
			break;
		default:
			gtk_toggle_button_set_active(b16bitu, TRUE);
		}
		whichminisave = TYPEFITS;
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving FITS"));
		gtk_window_set_transient_for (GTK_WINDOW(savepopup), GTK_WINDOW(lookup_widget("rgb_window")));
		gtk_notebook_set_current_page(notebookFormat, PAGE_FITS);
		gtk_widget_set_visible(savetxt, TRUE);
		gtk_widget_show(savepopup);
	}
}

void on_menu_rgb_savetiff_activate(GtkMenuItem *menuitem, gpointer user_data) {
	static GtkNotebook* notebookFormat = NULL;
	static GtkWidget *savepopup = NULL;
	GtkWidget *savetxt = lookup_widget("filenameframe");


	if (notebookFormat == NULL) {
		notebookFormat = GTK_NOTEBOOK(
				gtk_builder_get_object(builder, "notebookFormat"));
		savepopup = lookup_widget("savepopup");
	}

	if (single_image_is_loaded() || sequence_is_loaded()) {
		whichminisave = TYPETIFF;
		set_programm_name_in_TIFF(); //Write "Siril Version X.Y in Copyright_Txt
		set_description_in_TIFF();
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving TIFF"));
		gtk_window_set_transient_for (GTK_WINDOW(savepopup), GTK_WINDOW(lookup_widget("rgb_window")));
		gtk_notebook_set_current_page(notebookFormat, PAGE_TIFF);
		gtk_widget_set_visible(savetxt, TRUE);
		gtk_widget_show(savepopup);
	}
}

void on_menu_rgb_savepng_activate(GtkMenuItem *menuitem, gpointer user_data) {
	static GtkNotebook* notebookFormat = NULL;
	static GtkWidget *savepopup = NULL;
	GtkWidget *savetxt = lookup_widget("filenameframe");


	if (notebookFormat == NULL) {
		notebookFormat = GTK_NOTEBOOK(
				gtk_builder_get_object(builder, "notebookFormat"));
		savepopup = lookup_widget("savepopup");
	}

	if (single_image_is_loaded() || sequence_is_loaded()) {
		whichminisave = TYPEPNG;
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving PNG"));
		gtk_window_set_transient_for (GTK_WINDOW(savepopup), GTK_WINDOW(lookup_widget("rgb_window")));
		gtk_notebook_set_current_page(notebookFormat, PAGE_MISC);
		gtk_widget_set_visible(savetxt, TRUE);
		gtk_widget_show(savepopup);
	}
}

void on_menu_rgb_save8ppm_activate(GtkMenuItem *menuitem, gpointer user_data) {
	static GtkNotebook* notebookFormat = NULL;
	static GtkWidget *savepopup = NULL;
	GtkWidget *savetxt = lookup_widget("filenameframe");

	if (notebookFormat == NULL) {
		notebookFormat = GTK_NOTEBOOK(
				gtk_builder_get_object(builder, "notebookFormat"));
		savepopup = lookup_widget("savepopup");
	}

	if (single_image_is_loaded() || sequence_is_loaded()) {
		whichminisave = TYPEPNM;
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving Netpbm"));
		gtk_window_set_transient_for (GTK_WINDOW(savepopup), GTK_WINDOW(lookup_widget("rgb_window")));
		gtk_notebook_set_current_page(notebookFormat, PAGE_MISC);
		gtk_widget_set_visible(savetxt, TRUE);
		gtk_widget_show(savepopup);
	}
}

void on_menu_rgb_savebmp_activate(GtkMenuItem *menuitem, gpointer user_data) {
	static GtkNotebook* notebookFormat = NULL;
	static GtkWidget *savepopup = NULL;
	GtkWidget *savetxt = lookup_widget("filenameframe");

	if (notebookFormat == NULL) {
		notebookFormat = GTK_NOTEBOOK(
				gtk_builder_get_object(builder, "notebookFormat"));
		savepopup = lookup_widget("savepopup");
	}

	if (single_image_is_loaded() || sequence_is_loaded()) {
		whichminisave = TYPEBMP;
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving BMP"));
		gtk_window_set_transient_for (GTK_WINDOW(savepopup), GTK_WINDOW(lookup_widget("rgb_window")));
		gtk_notebook_set_current_page(notebookFormat, PAGE_MISC);
		gtk_widget_set_visible(savetxt, TRUE);
		gtk_widget_show(savepopup);
	}
}

void on_menu_rgb_savejpg_activate(GtkMenuItem *menuitem, gpointer user_data) {
	static GtkNotebook* notebookFormat = NULL;
	static GtkWidget *savepopup = NULL;
	GtkWidget *savetxt = lookup_widget("filenameframe");

	if (notebookFormat == NULL) {
		notebookFormat = GTK_NOTEBOOK(
				gtk_builder_get_object(builder, "notebookFormat"));
		savepopup = lookup_widget("savepopup");
	}

	if (single_image_is_loaded() || sequence_is_loaded()) {
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving JPG"));
		gtk_window_set_transient_for (GTK_WINDOW(savepopup), GTK_WINDOW(lookup_widget("rgb_window")));
		if (sequence_is_loaded() && !single_image_is_loaded()) {
			char filename[256];
			/* set the output file name default as the current image.jpg */
			GtkEntry *entry = GTK_ENTRY(lookup_widget("savetxt"));
			seq_get_image_filename(&com.seq, com.seq.current, filename);
			gtk_entry_set_text(entry, filename);
		}
		whichminisave = TYPEJPG;
		gtk_notebook_set_current_page(notebookFormat, PAGE_JPG);
		gtk_widget_set_visible(savetxt, TRUE);
		gtk_widget_show(savepopup);
	}
}

void on_savetxt_changed(GtkEditable *editable, gpointer user_data) {
	GtkEntry *entry = GTK_ENTRY(editable);
	GtkWidget *button = lookup_widget("button_savepopup");

	const gchar *name = gtk_entry_get_text(entry);
	gtk_widget_set_sensitive(button, (*name != '\0'));
}

void on_button_savepopup_clicked(GtkButton *button, gpointer user_data) {
	struct savedial_data *args = malloc(sizeof(struct savedial_data));

	set_cursor_waiting(TRUE);
	initialize_data(args);
	start_in_new_thread(minisavedial, args);
}

void on_savetxt_activate(GtkEntry *entry, gpointer user_data) {
	struct savedial_data *args = malloc(sizeof(struct savedial_data));

	set_cursor_waiting(TRUE);
	initialize_data(args);
	start_in_new_thread(minisavedial, args);
}

void on_button_cancelpopup_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("savepopup"));
}

void on_save1_activate(GtkMenuItem *menuitem, gpointer user_data) {
	GtkWidget *savepopup = lookup_widget("savepopup");

	if (save_dialog() == GTK_RESPONSE_ACCEPT) {
		/* now it is not needed for some formats */
		if (whichminisave != TYPEBMP && whichminisave != TYPEPNG
				&& whichminisave != TYPEPNM) {
			close_dialog();
			GtkWindow *parent = siril_get_active_window();
			if (!GTK_IS_WINDOW(parent)) {
				parent = GTK_WINDOW(lookup_widget("control_window"));
			}
			gtk_window_set_transient_for(GTK_WINDOW(savepopup),	parent);
			gtk_widget_show(savepopup);
			gtk_window_present_with_time(GTK_WINDOW(savepopup), GDK_CURRENT_TIME);
		}
		else {
			struct savedial_data *args = malloc(sizeof(struct savedial_data));

			set_cursor_waiting(TRUE);
			initialize_data(args);
			start_in_new_thread(minisavedial, args);
		}
	}
}

void on_savepopup_show(GtkWidget *widget, gpointer user_data) {
	GtkScrolledWindow *scrolled_window = GTK_SCROLLED_WINDOW(lookup_widget("scrolledwindow3"));
	gint height, width;

	if (whichminisave == TYPETIFF) {
		width = 400;
		height = 100;
	} else {
		width = 100;
		height = 50;
	}

	gtk_scrolled_window_set_min_content_height(scrolled_window, height);
	gtk_scrolled_window_set_min_content_width(scrolled_window, width);
}
