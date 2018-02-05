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

#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "gui/progress_and_log.h"
#include "gui/callbacks.h"
#include "io/sequence.h"
#include "io/single_image.h"


static image_type whichminisave;

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
	int i=0;

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

static void prepare_savepopup(int type) {
	static GtkNotebook* notebookFormat = NULL;
	static GtkWidget *savepopup = NULL;
	static GtkWidget *savetxt = NULL;
	int tab;

	if (notebookFormat == NULL) {
		notebookFormat = GTK_NOTEBOOK(
				gtk_builder_get_object(builder, "notebookFormat"));
		savepopup = lookup_widget("savepopup");
		savetxt = lookup_widget("filenameframe");
	}

	switch (type) {
	case TYPEBMP:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving BMP"));
		tab = 3;
		break;
	case TYPEPNG:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving PNG"));
		tab = 3;
		break;
	case TYPEPNM:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving Netpbm"));
		tab = 3;
		break;
	case TYPEJPG:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving JPG"));
		tab = 1;
		break;
	case TYPETIFF:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving TIFF"));
		tab = 0;
		break;
	case TYPEFITS:
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving FITS"));
		tab = 2;
	}

	whichminisave = type;
	gtk_widget_set_visible(savetxt, FALSE);
	gtk_notebook_set_current_page(notebookFormat, tab);
}

static int save_dialog() {
	GtkWidget *dialog;
	GtkFileChooser *chooser;
	GtkFileFilter *filter;
	GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_SAVE;
	GtkWindow *parent_window = GTK_WINDOW(lookup_widget("control_window"));
	GtkEntry *savetext = GTK_ENTRY(lookup_widget("savetxt"));
	gint res;

	dialog = gtk_file_chooser_dialog_new(_("Save File"), parent_window, action,
			_("_Cancel"), GTK_RESPONSE_CANCEL, _("_Save"), GTK_RESPONSE_ACCEPT,
			NULL);

	chooser = GTK_FILE_CHOOSER(dialog);

	gtk_file_chooser_set_do_overwrite_confirmation(chooser, TRUE);
	set_filters_save_dialog(chooser);

	res = gtk_dialog_run(GTK_DIALOG(dialog));
	if (res == GTK_RESPONSE_ACCEPT) {

		gchar *filename = gtk_file_chooser_get_filename(chooser);
		gtk_entry_set_text(savetext, filename);
		g_free(filename);

		/* we get the type of filter */
		filter = gtk_file_chooser_get_filter(chooser);
		const gchar *str = gtk_file_filter_get_name(filter);
		prepare_savepopup(get_filetype(str));
	}

	gtk_widget_destroy(dialog);
	return res;
}

static void minisavedial(void) {
	const gchar *filename;
	GtkToggleButton *fits_8 = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_save_fit8"));
	GtkToggleButton *fits_16s = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_save_fit16s"));
#ifdef HAVE_LIBJPEG
	GtkSpinButton *qlty_spin_button = GTK_SPIN_BUTTON(lookup_widget("quality_spinbutton"));
	gint quality = gtk_spin_button_get_value_as_int(qlty_spin_button);
#endif
#ifdef HAVE_LIBTIFF
	gint bitspersamples = 16;
	GtkToggleButton *BPS_Button = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton8bits"));

	if (gtk_toggle_button_get_active(BPS_Button))
		bitspersamples = 8;
#endif
	GtkEntry *entry = GTK_ENTRY(lookup_widget("savetxt"));

	filename = gtk_entry_get_text(entry);
	if (filename[0] != '\0') {
		switch (whichminisave) {
		case TYPEBMP:
				savebmp(filename, &gfit);
			break;
#ifdef HAVE_LIBJPEG
		case TYPEJPG:
				savejpg(filename, &gfit, quality);
			break;
#endif
#ifdef HAVE_LIBTIFF
		case TYPETIFF:
				savetif(filename, &gfit, bitspersamples);
			break;
#endif
#ifdef HAVE_LIBPNG
		case TYPEPNG:
				savepng(filename, &gfit, 2, gfit.naxes[2] == 3);
			break;
#endif
		case TYPEFITS:
			if (gtk_toggle_button_get_active(fits_8))
				gfit.bitpix = BYTE_IMG;
			else if (gtk_toggle_button_get_active(fits_16s))
				gfit.bitpix = SHORT_IMG;
			else
				gfit.bitpix = USHORT_IMG;
			/* Check if MIPS-HI and MIPS-LO must be updated. If yes,
			 * Values are taken from the layer 0 */
			if (gtk_toggle_button_get_active(
					GTK_TOGGLE_BUTTON(
							lookup_widget(
									"checkbutton_update_hilo"))) == TRUE) {
				if (sequence_is_loaded() && !single_image_is_loaded()) {
					gfit.hi = com.seq.layers[RLAYER].hi;
					gfit.lo = com.seq.layers[RLAYER].lo;
				} else {
					gfit.hi = com.uniq->layers[RLAYER].hi;
					gfit.lo = com.uniq->layers[RLAYER].lo;
				}
				if (gfit.bitpix == BYTE_IMG
						&& (gfit.hi > UCHAR_MAX || gfit.lo > UCHAR_MAX)) {
					gfit.hi = UCHAR_MAX;
					gfit.lo = 0;
				} else if (gfit.bitpix == SHORT_IMG
						&& (gfit.hi > SHRT_MAX || gfit.lo > SHRT_MAX)) {
					gfit.hi = UCHAR_MAX;
					gfit.lo = 0;
				}

			}
				savefits(filename, &gfit);
			break;
		case TYPEPNM:
			saveNetPBM(filename, &gfit);
			break;
		default:
			siril_log_message(
					_("This type of file is not handled. Should not happen"));
			break;
		}
		gtk_widget_hide(lookup_widget("savepopup"));
		gtk_entry_set_text(entry, "");
	}
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
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving FITS"));
		gtk_widget_show_all(savepopup);
		gtk_notebook_set_current_page(notebookFormat, 2);
		gtk_widget_set_visible(savetxt, TRUE);
		whichminisave = TYPEFITS;
	}
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
}

static void Set_Programm_name_in_TIFF() {
	static GtkTextView *TIFF_txt = NULL;
	GtkTextBuffer *tbuf;
	GtkTextIter itDebut, itFin;
	char Copyright[64];

	if (TIFF_txt == NULL)
		TIFF_txt = GTK_TEXT_VIEW(lookup_widget("Copyright_txt"));

	tbuf = gtk_text_view_get_buffer(TIFF_txt);

	g_snprintf(Copyright, sizeof(Copyright), "%s v%s", PACKAGE, VERSION);
	Copyright[0] = toupper(Copyright[0]);			// convert siril to Siril

	gtk_text_buffer_get_bounds(tbuf, &itDebut, &itFin);
	gtk_text_buffer_delete(tbuf, &itDebut, &itFin);
	gtk_text_buffer_set_text(tbuf, Copyright, strlen(Copyright));
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
		Set_Programm_name_in_TIFF(); //Write "Siril Version X.Y in Copyright_Txt
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving TIFF"));
		gtk_widget_show_all(savepopup);
		gtk_notebook_set_current_page(notebookFormat, 0);
		gtk_widget_set_visible(savetxt, TRUE);
		whichminisave = TYPETIFF;
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
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving PNG"));
		gtk_widget_show_all(savepopup);
		gtk_notebook_set_current_page(notebookFormat, 3);
		gtk_widget_set_visible(savetxt, TRUE);
		whichminisave = TYPEPNG;
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
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving Netpbm"));
		gtk_widget_show_all(savepopup);
		gtk_notebook_set_current_page(notebookFormat, 3);
		gtk_widget_set_visible(savetxt, TRUE);
		whichminisave = TYPEPNM;
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
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving BMP"));
		gtk_widget_show_all(savepopup);
		gtk_notebook_set_current_page(notebookFormat, 3);
		gtk_widget_set_visible(savetxt, TRUE);
		whichminisave = TYPEBMP;
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
		if (sequence_is_loaded() && !single_image_is_loaded()) {
			char filename[256];
			/* set the output file name default as the current image.jpg */
			GtkEntry *entry = GTK_ENTRY(lookup_widget("savetxt"));
			seq_get_image_filename(&com.seq, com.seq.current, filename);
			gtk_entry_set_text(entry, filename);
		}
		gtk_widget_show_all(savepopup);
		gtk_notebook_set_current_page(notebookFormat, 1);
		gtk_widget_set_visible(savetxt, TRUE);
		whichminisave = TYPEJPG;
	}
}


gboolean on_savetxt_key_press_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {
	GtkWidget *button = lookup_widget("button_savepopup");
	gboolean handled = FALSE;

	switch (event->keyval) {
	case GDK_KEY_Return:
	case GDK_KEY_KP_Enter:
		handled = TRUE;
		gtk_widget_set_can_default(button, TRUE);
		gtk_widget_grab_focus(widget);
		set_cursor_waiting(TRUE);
		minisavedial();
		set_cursor_waiting(FALSE);
		break;
	}
	return handled;
}

void on_savetxt_changed(GtkEditable *editable, gpointer user_data) {
	GtkEntry *entry = GTK_ENTRY(editable);
	GtkWidget *button = lookup_widget("button_savepopup");

	const gchar *name = gtk_entry_get_text(entry);
	gtk_widget_set_sensitive(button, (*name != '\0'));
}

void on_button_savepopup_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	minisavedial();
	set_cursor_waiting(FALSE);
}

void on_button_cancelpopup_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("savepopup"));
}

void on_save1_activate(GtkMenuItem *menuitem, gpointer user_data) {
	static GtkWidget *savepopup = NULL;


	if (savepopup == NULL) {
		savepopup = lookup_widget("savepopup");
	}

	int res = save_dialog();
	if (res == GTK_RESPONSE_ACCEPT) {
		/* now it is not needed for some formats */
		if (whichminisave != TYPEBMP && whichminisave != TYPEPNG
				&& whichminisave != TYPEPNM)
			gtk_widget_show(savepopup);
		else {
			minisavedial();
		}
	}
}

