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
#include "gui/dialogs.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "io/sequence.h"
#include "algos/PSF.h"
#include "registration/registration.h"	// for update_reg_interface
#include "stacking/stacking.h"	// for update_stack_interface

#include "sequence_list.h"

static gboolean fill_sequence_list_idle(gpointer p);

static const char *bg_colour[] = { "WhiteSmoke", "#1B1B1B" };
static const char *ref_bg_colour[] = { "Beige", "#4A4A39" };

struct _seq_list {
	sequence *seq;
	int layer;
};


static GtkListStore *list_store = NULL;

enum {
	COLUMN_IMNAME,		// string
	COLUMN_SHIFTX,		// int
	COLUMN_SHIFTY,		// int
	COLUMN_SELECTED,	// gboolean
	COLUMN_FWHM,		// gdouble
	COLUMN_CURRENT,		// int weight, current file loaded, display IMNAME in bold
	COLUMN_REFERENCE,	// background color depending on the image being reference
	COLUMN_INDEX,		// int
	N_COLUMNS
};

/******* Static functions **************/
static void fwhm_quality_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble quality;
	gchar buf[20];
	gtk_tree_model_get(model, iter, COLUMN_FWHM, &quality, -1);
	if (quality >= 0.0)
		g_snprintf(buf, sizeof(buf), "%.3f", quality);
	else
		g_strlcpy(buf, "N/A", sizeof(buf));
	g_object_set(renderer, "text", buf, NULL);
}

static void initialize_seqlist_dialog_combo() {
	if (!sequence_is_loaded()) return;

	GtkComboBoxText *seqcombo = GTK_COMBO_BOX_TEXT(lookup_widget("seqlist_dialog_combo"));
	g_signal_handlers_block_by_func(GTK_COMBO_BOX(seqcombo), on_seqlist_dialog_combo_changed, NULL);
	gtk_combo_box_text_remove_all(seqcombo);

	if (com.seq.nb_layers == 1) {
		gtk_combo_box_text_append_text(seqcombo, _("B&W channel"));
	} else {
		gtk_combo_box_text_append_text(seqcombo, _("Red channel"));
		gtk_combo_box_text_append_text(seqcombo, _("Green channel"));
		gtk_combo_box_text_append_text(seqcombo, _("Blue channel"));
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(seqcombo), 0);
	g_signal_handlers_unblock_by_func(GTK_COMBO_BOX(seqcombo), on_seqlist_dialog_combo_changed, NULL);
}

static void initialize_title() {
	GtkHeaderBar *bar = GTK_HEADER_BAR(lookup_widget("seqlistbar"));
	gchar *seq_basename;

	if (sequence_is_loaded())
		seq_basename = g_path_get_basename(com.seq.seqname);
	else
		seq_basename = g_strdup(_("No sequence loaded"));

	gtk_header_bar_set_subtitle(bar, seq_basename);

	g_free(seq_basename);
}

static void initialize_search_entry() {
	GtkTreeView *tree_view = GTK_TREE_VIEW(gtk_builder_get_object(builder, "treeview1"));
	GtkEntry *entry = GTK_ENTRY(lookup_widget("seqlistsearch"));

	gtk_entry_set_text (entry, "");
	gtk_tree_view_set_search_entry(tree_view, entry);
}

static void get_list_store() {
	if (list_store == NULL) {
		list_store = GTK_LIST_STORE(gtk_builder_get_object(builder, "liststore1"));

		GtkTreeViewColumn *col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(builder, "treeviewcolumn5"));
		GtkCellRenderer *cell = GTK_CELL_RENDERER(gtk_builder_get_object(builder, "cellrenderertext5"));
		gtk_tree_view_column_set_cell_data_func(col, cell, fwhm_quality_cell_data_function, NULL, NULL);
	}
}

/* Add an image to the list. If seq is NULL, the list is cleared. */
static void add_image_to_sequence_list(sequence *seq, int index, int layer) {
	static GtkTreeSelection *selection = NULL;
	GtkTreeIter iter;
	char imname[256];
	char *basename;
	int shiftx = -1, shifty = -1;
	double fwhm = -1.0;
	int color;

	get_list_store();
	if (!selection)
		selection = GTK_TREE_SELECTION(gtk_builder_get_object(builder, "treeview-selection1"));
	if (seq == NULL) {
		gtk_list_store_clear(list_store);
		return;		// just clear the list
	}
	if (seq->regparam && seq->regparam[layer]) {
		shiftx = roundf_to_int(seq->regparam[layer][index].shiftx);
		shifty = roundf_to_int(seq->regparam[layer][index].shifty);

		if (seq->regparam[layer][index].fwhm > 0.0f)
			fwhm = seq->regparam[layer][index].fwhm;
		else if (seq->regparam[layer][index].quality >= 0.0)
			fwhm = seq->regparam[layer][index].quality;
	}

	color = (com.combo_theme == 0) ? 1 : 0;

	basename = g_path_get_basename(seq_get_image_filename(seq, index, imname));
	gtk_list_store_append (list_store, &iter);
	gtk_list_store_set (list_store, &iter,
			COLUMN_IMNAME, basename,
			COLUMN_SHIFTX, shiftx,
			COLUMN_SHIFTY, shifty,
			COLUMN_SELECTED, seq->imgparam[index].incl,
			COLUMN_FWHM, fwhm,
			COLUMN_CURRENT, index == seq->current ? 800 : 400,
			// weight value is 400 by default "normal":
			// http://developer.gnome.org/gtk3/stable/GtkCellRendererText.html#GtkCellRendererText--weight
			COLUMN_REFERENCE, index == seq->reference_image ?
			ref_bg_colour[color] : bg_colour[color],
			COLUMN_INDEX, (index + 1),
			-1);
	/* see example at http://developer.gnome.org/gtk3/3.5/GtkListStore.html */
	if (index == seq->current) {
		gtk_tree_selection_select_iter(selection, &iter);
	}
	g_free(basename);
}

static void sequence_list_change_selection(gchar *path, gboolean new_value) {
	GtkTreeIter iter;
	get_list_store();
	gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(list_store), &iter, path);
	gtk_list_store_set(list_store, &iter, COLUMN_SELECTED, new_value, -1);
}

static GList *get_row_references_of_selected_rows(GtkTreeSelection *selection,
		GtkTreeModel *model) {
	GList *ref = NULL;
	GList *sel, *s;

	sel = gtk_tree_selection_get_selected_rows(selection, &model);

	for (s = sel; s; s = s->next) {
		GtkTreeRowReference *rowref = gtk_tree_row_reference_new(model,	(GtkTreePath *) s->data);
		ref = g_list_prepend(ref, rowref);
	}
	g_list_free_full(sel, (GDestroyNotify) gtk_tree_path_free);
	return ref;
}

static void unselect_select_frame_from_list(gboolean select) {
	GtkTreeSelection *selection;
	GtkTreeModel *model;
	GList *references, *list;
	GtkTreeView *tree_view;
	get_list_store();

	tree_view = GTK_TREE_VIEW(lookup_widget("treeview1"));
	model = gtk_tree_view_get_model(tree_view);
	selection = gtk_tree_view_get_selection(tree_view);
	references = get_row_references_of_selected_rows(selection, model);
	for (list = references; list; list = list->next) {
		GtkTreeIter iter;
		GtkTreePath *path = gtk_tree_row_reference_get_path((GtkTreeRowReference*)list->data);
		if (path) {
			if (gtk_tree_model_get_iter(model, &iter, path)) {
				gint idx;
				gtk_tree_model_get(model, &iter, COLUMN_INDEX, &idx, -1);
				idx --; // we count frames from 1 !!
				toggle_image_selection(idx);
			}
			gtk_tree_path_free(path);
		}
	}
	g_list_free(references);
}

/**** Callbacs *****/

void on_treeview1_row_activated(GtkTreeView *tree_view, GtkTreePath *path,
		GtkTreeViewColumn *column, gpointer user_data) {
	GtkTreeModel *tree_model;
	GtkTreeSelection *selection;
	GtkTreeIter iter;
	GList *list;

	tree_model = gtk_tree_view_get_model(tree_view);
	selection = gtk_tree_view_get_selection (tree_view);

	list = gtk_tree_selection_get_selected_rows(selection, &tree_model);
	if (g_list_length(list) == 1) {
		gint idx;
		GValue value = G_VALUE_INIT;
		gtk_tree_model_get_iter(tree_model, &iter, path);

		gtk_tree_model_get_value(tree_model, &iter, COLUMN_INDEX, &value);
		idx = g_value_get_int(&value) - 1;
		if (idx != com.seq.current) {
			fprintf(stdout, "loading image %d\n", idx);
			seq_load_image(&com.seq, idx, TRUE);
		}
		g_value_unset(&value);
	}
	g_list_free_full(list, (GDestroyNotify) gtk_tree_path_free);
}

void on_seqlist_dialog_combo_changed(GtkComboBoxText *widget, gpointer user_data) {
	int active = gtk_combo_box_get_active(GTK_COMBO_BOX(widget));
	if (active >= 0) {
		fill_sequence_list(&com.seq, active, FALSE);
	}
}

void initialize_seqlist() {
	initialize_title();
	initialize_seqlist_dialog_combo();
	initialize_search_entry();
}

/* called on sequence loading (set_seq), on layer tab change and on registration data update.
 * It is executed safely in the GTK thread if as_idle is true. */
void fill_sequence_list(sequence *seq, int layer, gboolean as_idle) {
	struct _seq_list *args;
	if (seq == NULL || layer >= seq->nb_layers) return;
	args = malloc(sizeof(struct _seq_list));
	args->seq = seq;
	args->layer = layer;
	if (as_idle)
		gdk_threads_add_idle(fill_sequence_list_idle, args);
	else fill_sequence_list_idle(args);
}

static gboolean fill_sequence_list_idle(gpointer p) {
	int i;
	struct _seq_list *args = (struct _seq_list *)p;
	add_image_to_sequence_list(NULL, 0, 0);	// clear
	if (args->seq->number > 0) {
		for (i = 0; i < args->seq->number; i++) {
			add_image_to_sequence_list(args->seq, i, args->layer);
		}
	}
	free(args);
	return FALSE;
}

void on_seqlist_button_clicked(GtkToolButton *button, gpointer user_data) {
	if (gtk_widget_get_visible(lookup_widget("seqlist_dialog"))) {
		siril_close_dialog("seqlist_dialog");
	} else {
		siril_open_dialog("seqlist_dialog");
	}
}

int get_image_index_from_path(GtkTreePath *path) {
	GValue value = G_VALUE_INIT;
	gint index;
	GtkTreeIter iter;
	get_list_store();
	gtk_tree_model_get_iter(GTK_TREE_MODEL(list_store), &iter, path);
	gtk_tree_model_get_value(GTK_TREE_MODEL(list_store), &iter, COLUMN_INDEX, &value);
	index = g_value_get_int(&value) - 1;
	g_value_unset(&value);
	return index;
}

void on_seqlist_image_selection_toggled(GtkCellRendererToggle *cell_renderer,
		gchar *path, gpointer user_data) {
	gint index = get_image_index_from_path(gtk_tree_path_new_from_string(path));
	if (index < 0 || index >= com.seq.number) return;
	fprintf(stdout, "toggle selection index = %d\n", index);

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
	update_stack_interface(FALSE);
	writeseqfile(&com.seq);
	redraw(com.cvport, REMAP_NONE);
}

void toggle_image_selection(int image_num) {
	gchar *msg;
	if (com.seq.imgparam[image_num].incl) {
		com.seq.imgparam[image_num].incl = FALSE;
		--com.seq.selnum;
		msg = g_strdup_printf(_("Image %d has been unselected from sequence\n"), image_num);
		if (image_num == com.seq.reference_image) {
			com.seq.reference_image = -1;
			sequence_list_change_reference();
			adjust_refimage(image_num);
		}
	} else {
		com.seq.imgparam[image_num].incl = TRUE;
		++com.seq.selnum;
		msg = g_strdup_printf(_("Image %d has been selected from sequence\n"), image_num);
	}
	siril_log_message(msg);
	g_free(msg);
	sequence_list_change_selection_index(image_num);
	update_reg_interface(FALSE);
	update_stack_interface(TRUE);
	adjust_exclude(image_num, TRUE);
	writeseqfile(&com.seq);
}

void on_selected_frames_unselect(GtkButton *button, gpointer user_data) {
	unselect_select_frame_from_list(FALSE);
}

void on_selected_frames_select(GtkButton *button, gpointer user_data) {
	unselect_select_frame_from_list(TRUE);
}

/****************** modification of the list store (tree model) ******************/

void sequence_list_change_selection_index(int index) {
	sequence_list_change_selection(
			gtk_tree_path_to_string(gtk_tree_path_new_from_indices(index, -1)),
			com.seq.imgparam[index].incl);
}

void sequence_list_change_current() {
	GtkTreeIter iter;
	gboolean valid;

	get_list_store();
	valid = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(list_store), &iter);
	while (valid) {
		gint index = 0;
		GValue value = G_VALUE_INIT;

		gtk_tree_model_get_value (GTK_TREE_MODEL(list_store), &iter, COLUMN_INDEX, &value);
		index = g_value_get_int(&value) - 1;
		g_value_unset(&value);
		gtk_list_store_set(list_store, &iter, COLUMN_CURRENT,
				(index == com.seq.current) ? 800 : 400, -1);
		valid = gtk_tree_model_iter_next(GTK_TREE_MODEL(list_store), &iter);
	}
}

void sequence_list_change_reference() {
	GtkTreeIter iter;
	gboolean valid;
	gint row_count = 0;
	int color;

	color = (com.combo_theme == 0) ? 1 : 0;

	get_list_store();
	valid = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(list_store), &iter);
	while (valid) {
		gtk_list_store_set(list_store, &iter,
				COLUMN_REFERENCE,
				(row_count == com.seq.reference_image) ?
				ref_bg_colour[color] : bg_colour[color], -1);
		row_count++;
		valid = gtk_tree_model_iter_next(GTK_TREE_MODEL(list_store), &iter);
	}
}

void clear_sequence_list() {
	get_list_store();
	gtk_list_store_clear(list_store);
}
