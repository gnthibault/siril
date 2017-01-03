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
#include "algos/PSF.h"
#include "registration/registration.h"	// for update_reg_interface
#include "stacking/stacking.h"	// for update_stack_interface

static gboolean fill_sequence_list_idle(gpointer p);

static const char *bg_colour[] = { "WhiteSmoke", "#1B1B1B" };
static const char *ref_bg_colour[] = { "Beige", "#4A4A39" };


static GtkListStore *list_store = NULL;

enum {
	COLUMN_IMNAME,		// string
	COLUMN_SHIFTX,		// int
	COLUMN_SHIFTY,		// int
	COLUMN_SELECTED,	// gboolean
	COLUMN_FWHM,		// converted to string, not pure double
	COLUMN_CURRENT,		// int weight, current file loaded, display IMNAME in bold
	COLUMN_REFERENCE,	// background color depending on the image being reference
	COLUMN_INDEX,		// int
	N_COLUMNS
};

void get_list_store() {
	if (list_store == NULL)
		list_store = GTK_LIST_STORE(gtk_builder_get_object(builder, "liststore1"));
}

/* Add an image to the list. If seq is NULL, the list is cleared. */
void add_image_to_sequence_list(sequence *seq, int index, int layer) {
	static GtkTreeSelection *selection = NULL;
	GtkTreeIter iter;
	char imname[256], fwhm_str[20];
	char *basename;
	int shiftx = -1, shifty = -1;

	get_list_store();
	if (!selection)
		selection = GTK_TREE_SELECTION(gtk_builder_get_object(builder, "treeview-selection1"));
	if (seq == NULL) {
		gtk_list_store_clear(list_store);
		return;		// just clear the list
	}
	if (seq->regparam && seq->regparam[layer]) {
		shiftx = seq->regparam[layer][index].shiftx;
		shifty = seq->regparam[layer][index].shifty;

		if (seq->regparam[layer][index].fwhm > 0.0f) {
			// Is it certain that FWHMX is > than FWHMY? The minimization seems to imply this.
			sprintf(fwhm_str, "%.3f", seq->regparam[layer][index].fwhm);
		} else if (seq->regparam[layer][index].quality >= 0.0) {
			sprintf(fwhm_str, "%.3f", seq->regparam[layer][index].quality);
		} else sprintf(fwhm_str, "N/A");
	} else sprintf(fwhm_str, "N/A");

	basename = g_path_get_basename(seq_get_image_filename(seq, index, imname));
	gtk_list_store_append (list_store, &iter);
	gtk_list_store_set (list_store, &iter,
			COLUMN_IMNAME, basename,
			COLUMN_SHIFTX, shiftx,
			COLUMN_SHIFTY, shifty,
			COLUMN_SELECTED, seq->imgparam[index].incl,
			COLUMN_FWHM, fwhm_str,
			COLUMN_CURRENT, index == seq->current ? 800 : 400,
			// weight value is 400 by default "normal":
			// http://developer.gnome.org/gtk3/stable/GtkCellRendererText.html#GtkCellRendererText--weight
			COLUMN_REFERENCE, index == seq->reference_image ?
			ref_bg_colour[com.have_dark_theme] : bg_colour[com.have_dark_theme],
			COLUMN_INDEX, index,
			-1);
	/* see example at http://developer.gnome.org/gtk3/3.5/GtkListStore.html */
	if (index == seq->current) {
		gtk_tree_selection_select_iter(selection, &iter);
	}
	g_free(basename);
}

struct _seq_list {
	sequence *seq;
	int layer;
};

/* called on sequence loading (set_seq), on layer tab change and on registration data update.
 * It is executed safely in the GTK thread. */
void fill_sequence_list(sequence *seq, int layer) {
	struct _seq_list *args;
	if (seq == NULL || layer >= seq->nb_layers) return;
	args = malloc(sizeof(struct _seq_list));
	args->seq = seq;
	args->layer = layer;
	gdk_threads_add_idle(fill_sequence_list_idle, args);
}

static gboolean fill_sequence_list_idle(gpointer p) {
	int i;
	struct _seq_list *args = (struct _seq_list *)p;
	add_image_to_sequence_list(NULL, 0, 0);	// clear  
	if (args->seq->number == 0) return FALSE;
	for (i=0; i<args->seq->number; i++) {
		add_image_to_sequence_list(args->seq, i, args->layer);
	}
	free(args);
	return FALSE;
}

void show_seqlist(GtkWidget *widget, gboolean show) {
	static gboolean was_extended = FALSE;
	if (!was_extended) {
		gint w, h;
		GtkWindow *window = GTK_WINDOW(lookup_widget("main_window"));
		gtk_window_get_size(window, &w, &h);
		gtk_window_resize(window, w+200, h);
		was_extended = TRUE;
	}
	gtk_paned_set_position(GTK_PANED(widget), show ? 200 : 0);
}

void on_toggle_show_seqlist_toggled(GtkToggleToolButton *togglebutton, gpointer user_data) {
	static GtkWidget *paned = NULL;
	if (!paned)
		paned = lookup_widget("paned1");
	show_seqlist(paned, gtk_toggle_tool_button_get_active(togglebutton));
}

int get_image_index_from_path(GtkTreePath *path) {
	GValue value = G_VALUE_INIT;
	gint index;
	GtkTreeIter iter;
	get_list_store();
	gtk_tree_model_get_iter(GTK_TREE_MODEL(list_store), &iter, path);
	gtk_tree_model_get_value(GTK_TREE_MODEL(list_store), &iter, COLUMN_INDEX, &value);
	index = g_value_get_int(&value);
	g_value_unset(&value);
	return index;
}

void on_seqlist_image_selection_toggled(GtkCellRendererToggle *cell_renderer,
		gchar *path, gpointer user_data) {
	gint index = get_image_index_from_path(gtk_tree_path_new_from_string(path));
	if (index < 0 || index >= com.seq.number) return;
	fprintf(stdout, "toggle selection index = %d\n", index);
	/*if (gtk_cell_renderer_toggle_get_active(cell_renderer) ==
			com.seq.imgparam[index].incl) {
		fprintf(stdout, "mismatch in selection toggle, what should I do?\n");
		return;
	}*/

	sequence_list_change_selection(path, !com.seq.imgparam[index].incl);
	siril_log_message(_("%s image %d in sequence %s\n"),
			com.seq.imgparam[index].incl ? _("excluding") : _("including"),
			index, com.seq.seqname);

	com.seq.imgparam[index].incl = !com.seq.imgparam[index].incl;
	if (com.seq.imgparam[index].incl)
		com.seq.selnum++;
	else 	com.seq.selnum--;
	adjust_exclude(index, TRUE);	// check or uncheck excluded checkbox in seq tab
	update_reg_interface(FALSE);
	update_stack_interface();
	writeseqfile(&com.seq);
	redraw(com.cvport, REMAP_NONE);
}

/* double click on an image -> open it */
void on_treeview1_row_activated(GtkTreeView *tree_view, GtkTreePath *path,
		GtkTreeViewColumn *column, gpointer user_data) {
	gint index = get_image_index_from_path(path);
	if (index < 0 || index >= com.seq.number) return;
	fprintf(stdout, "loading image %d\n", index);
	seq_load_image(&com.seq, index, &gfit, TRUE);
}

/****************** modification of the list store (tree model) ******************/

void sequence_list_change_selection(gchar *path, gboolean new_value) {
	GtkTreeIter iter;
	get_list_store();
	gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(list_store), &iter, path);
	gtk_list_store_set(list_store, &iter, COLUMN_SELECTED, new_value, -1);
}

void sequence_list_change_selection_index(int index) {
	sequence_list_change_selection(
			gtk_tree_path_to_string(gtk_tree_path_new_from_indices(index, -1)),
			com.seq.imgparam[index].incl);
}

void sequence_list_change_current() {
	GtkTreeIter iter;
	gboolean valid;
	gint row_count = 0;

	get_list_store();
	valid = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(list_store), &iter);
	while (valid) {
		gtk_list_store_set(list_store, &iter,
				COLUMN_CURRENT, (row_count == com.seq.current) ? 800 : 400,
				-1);
		row_count++;
		valid = gtk_tree_model_iter_next(GTK_TREE_MODEL(list_store), &iter);
	}
}

void sequence_list_change_reference() {
	GtkTreeIter iter;
	gboolean valid;
	gint row_count = 0;

	get_list_store();
	valid = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(list_store), &iter);
	while (valid) {
		gtk_list_store_set(list_store, &iter,
				COLUMN_REFERENCE,
				(row_count == com.seq.reference_image) ?
				ref_bg_colour[com.have_dark_theme] : bg_colour[com.have_dark_theme],
				-1);
		row_count++;
		valid = gtk_tree_model_iter_next(GTK_TREE_MODEL(list_store), &iter);
	}
}

void clear_sequence_list() {
	get_list_store();
	gtk_list_store_clear(list_store);
}
