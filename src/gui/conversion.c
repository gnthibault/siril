/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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

#include <gdk/gdkkeysyms.h>
#include <gtk/gtk.h>
#include <string.h>

#include "gui/conversion.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "io/conversion.h"
#include "io/sequence.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "algos/sorting.h"

static gchar *destroot = NULL;
static GtkListStore *liststore_convert = NULL;
static GtkTreeView *tree_view = NULL;
static void check_for_conversion_form_completeness();

static void init_tree_view() {
	if (tree_view  == NULL)
		tree_view = GTK_TREE_VIEW(gtk_builder_get_object(builder, "treeview_convert"));
	g_assert(tree_view);
}

int count_converted_files() {
	init_tree_view();
	GtkTreeIter iter;
	GtkTreeModel *model = gtk_tree_view_get_model(tree_view);
	gboolean valid = gtk_tree_model_get_iter_first(model, &iter);
	
	int count = 0;
	while (valid) {
		gtk_tree_model_get(model, &iter, -1);
		valid = gtk_tree_model_iter_next (model, &iter);
		count++;
	}
	return count;
}

int count_selected_files() {
	init_tree_view();
	GtkTreeSelection *selection = gtk_tree_view_get_selection(tree_view);
	return gtk_tree_selection_count_selected_rows(selection);
}

static void initialize_convert() {
	GDir *dir;
	GError *error = NULL;
	gchar *file_data, *file_date;
	const gchar *index;
	static GtkEntry *startEntry = NULL;
	GtkTreeModel *model = NULL;
	GtkTreeIter iter;
	GList *list = NULL;
	
	init_tree_view();
	if (!startEntry) {
		startEntry = GTK_ENTRY(lookup_widget("startIndiceEntry"));
	}
	
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	if (file_name_has_invalid_chars(destroot)) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Invalid char"), _("Please remove invalid characters in the sequence name "
				"before trying to convert images into a new sequence again."));
		return;
	}

	if (g_file_test(destroot, G_FILE_TEST_EXISTS)) {
		char *title = siril_log_message(_("A file named %s already exists. "
				"Do you want to replace it?\n"), destroot);
		gboolean replace = siril_confirm_dialog(title, _("The file already exists. "
				"Replacing it will overwrite its contents."));
		if (!replace) return;
	}

	model = gtk_tree_view_get_model(tree_view);
	gboolean valid = gtk_tree_model_get_iter_first(model, &iter);
	if (!valid) return;	//The tree is empty

	//image_type first_type = TYPEUNDEF;
	gboolean /*several_type_of_files = FALSE, */no_sequence_to_convert = TRUE;
	gboolean there_is_an_image = FALSE, there_is_a_ser = FALSE;
	int count = 0;
	while (valid) {
		gtk_tree_model_get(model, &iter, COLUMN_FILENAME, &file_data,
				COLUMN_DATE, &file_date, -1);
		list = g_list_prepend(list, file_data);

		const char *src_ext = get_filename_ext(file_data);
		image_type type = get_type_for_extension(src_ext);
		if (type == TYPESER)
			there_is_a_ser = TRUE;
		if (type == TYPEAVI || type == TYPESER)
			no_sequence_to_convert = FALSE;
		else there_is_an_image = TRUE;
		//if (count == 0)
		//	first_type = type;
		//if (type != first_type)
		//	several_type_of_files = TRUE;

		valid = gtk_tree_model_iter_next(model, &iter);
		count++;
	}

	/* handle impossible cases */
	/* why is it forbidden? and shouldn't it be CONVDSTSER instead?
	 * apparently SER cannot be converted to SER, wouldn't it be nice to debayer them? */
	if ((convflags & CONVDEBAYER) && there_is_a_ser) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("A conflict has been detected."),
				_("The Debayer option is not allowed in SER conversion, please uncheck the option."));
		g_list_free_full(list, g_free);
		return;
	}
	if ((convflags & CONVMULTIPLE) && there_is_a_ser) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("A conflict has been detected."),
				_("The Multiple SER option is not allowed in SER conversion, please uncheck the option."));
		g_list_free_full(list, g_free);
		return;
	}
	if ((convflags & CONVMULTIPLE) && there_is_an_image) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("A conflict has been detected."),
				_("Creating multiple SER can only be done with only films as input."));
		g_list_free_full(list, g_free);
		return;
	}

	index = gtk_entry_get_text(startEntry);

	siril_log_color_message(_("Conversion: processing %d files...\n"), "green", count);
	
	set_cursor_waiting(TRUE);
	control_window_switch_to_tab(OUTPUT_LOGS);
	
	/* then, convert files to Siril's FITS format */
	if((dir = g_dir_open(com.wd, 0, &error)) == NULL){
		char *tmpmsg = siril_log_message(_("Conversion: error opening working directory %s.\n"), com.wd);
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), tmpmsg);
		fprintf(stderr, "Conversion: %s\n", error->message);
		g_error_free(error);
		g_list_free_full(list, g_free);
		set_cursor_waiting(FALSE);
		return ;
	}

	/* g_list_append() has to traverse the entire list to find the end,
	 * which is inefficient when adding multiple elements. A common idiom
	 * to avoid the inefficiency is to use g_list_prepend() and reverse the
	 * list with g_list_reverse() when all elements have been added. */
	list = g_list_reverse(list);
	/* convert the list to an array for parallel processing */
	char **files_to_convert = malloc(count * sizeof(char *));
	if (!files_to_convert) {
		PRINT_ALLOC_ERR;
		return;
	}
	GList *orig_list = list;
	for (int i = 0; i < count && list; list = list->next, i++)
		files_to_convert[i] = g_strdup(list->data);
	g_list_free_full(orig_list, g_free);

	struct _convert_data *args = malloc(sizeof(struct _convert_data));
	if (!args) {
		PRINT_ALLOC_ERR;
		return;
	}
	args->start = (atoi(index) <= 0 || atoi(index) >= 100000) ? 1 : atoi(index);
	args->dir = dir;
	args->list = files_to_convert;
	args->total = count;
	args->nb_converted = 0;
	args->compatibility = com.debayer.compatibility;
	args->command_line = FALSE;
	args->input_has_a_seq = !no_sequence_to_convert;
	args->destroot = g_strdup(destroot);
	gettimeofday(&(args->t_start), NULL);
	start_in_new_thread(convert_thread_worker, args);
	return;
}

void on_entry2_activate(GtkEntry *entry, gpointer user_data) {
	initialize_convert();
}

void on_convert_button_clicked(GtkButton *button, gpointer user_data) {
	initialize_convert();
}

static void add_convert_to_list(char *filename, GStatBuf st) {
	GtkTreeIter iter;
	char *date;

	date = ctime(&st.st_mtime);
	date[strlen(date) - 1] = 0;	// removing '\n' at the end of the string

	gtk_list_store_append(liststore_convert, &iter);
	gtk_list_store_set(liststore_convert, &iter, COLUMN_FILENAME, filename,	// copied in the store
			COLUMN_DATE, date, -1);
}
static void get_convert_list_store() {
	if (liststore_convert == NULL)
		liststore_convert = GTK_LIST_STORE(
				gtk_builder_get_object(builder, "liststore_convert"));
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

static void remove_selected_files_from_list() {
	GtkTreeSelection *selection;
	GtkTreeModel *model;
	GList *references, *list;

	init_tree_view();
	model = gtk_tree_view_get_model(tree_view);
	selection = gtk_tree_view_get_selection(tree_view);
	references = get_row_references_of_selected_rows(selection, model);
	for (list = references; list; list = list->next) {
		GtkTreeIter iter;
		GtkTreePath *path = gtk_tree_row_reference_get_path((GtkTreeRowReference*)list->data);
		if (path) {
			if (gtk_tree_model_get_iter(model, &iter, path)) {
				gtk_list_store_remove(liststore_convert, &iter);
			}
			gtk_tree_path_free(path);
		}
	}
	g_list_free(references);
	gtk_tree_selection_unselect_all(selection);
}

static gint sort_conv_tree(GtkTreeModel *model, GtkTreeIter *a, GtkTreeIter *b,
		gpointer user_data) {
	gchar *name_a, *name_b;
	gchar *collate_key1, *collate_key2;
	gint ret;

	gtk_tree_model_get(model, a, 0,	&name_a, -1);
	gtk_tree_model_get(model, b, 0,	&name_b, -1);

	collate_key1  = g_utf8_collate_key_for_filename(name_a, strlen(name_a));
	collate_key2  = g_utf8_collate_key_for_filename(name_b, strlen(name_b));

	ret = g_strcmp0(collate_key1, collate_key2);

	g_free(collate_key1);
	g_free(collate_key2);
	g_free(name_a);
	g_free(name_b);

	return ret;
}

void fill_convert_list(GSList *list) {
	GStatBuf st;
	GSList *l;

	get_convert_list_store();

	for (l = list; l; l = l->next) {
		char *filename;

		filename = (char *) l->data;
		if (g_stat(filename, &st) == 0) {
			add_convert_to_list(filename, st);
		}
		g_free(filename);
	}
	check_for_conversion_form_completeness();
}

void on_clear_convert_button_clicked(GtkToolButton *button, gpointer user_data) {
	get_convert_list_store();
	gtk_list_store_clear(liststore_convert);
	check_for_conversion_form_completeness();
}

void on_remove_convert_button_clicked(GtkToolButton *button, gpointer user_data) {
	remove_selected_files_from_list();
	check_for_conversion_form_completeness();
}

void on_treeview_convert_drag_data_received(GtkWidget *widget,
		GdkDragContext *context, gint x, gint y,
		GtkSelectionData *selection_data, guint info, guint time,
		gpointer user_data) {

	gchar **uris, **str;
	const guchar *data;
	GSList *list = NULL;
	gint bad_files = 0;

	if (info != 0)
		return;

	data = gtk_selection_data_get_data(selection_data);
	uris = g_uri_list_extract_uris((gchar *) data);

	for (str = uris; *str; str++) {
		GError *error = NULL;
		gchar *path = g_filename_from_uri(*str, NULL, &error);
		if (path) {
			const char *src_ext = get_filename_ext(path);
			if (src_ext) {
				if (get_type_for_extension(src_ext) == TYPEUNDEF) {
					bad_files++;
				} else {
					list = g_slist_prepend(list, path);
				}
			} else bad_files++;
		} else {
			fprintf(stderr, "Could not convert uri to local path: %s",
					error->message);
			bad_files++;
			g_error_free(error);
		}
	}
	list = g_slist_sort(list, (GCompareFunc) strcompare);
	fill_convert_list(list);
	if (bad_files) {
		char *msg = siril_log_message(_("%d %s while drag and drop\n"), bad_files,
				ngettext("file was ignored", "files were ignored", bad_files));
		siril_message_dialog(GTK_MESSAGE_INFO, msg,
				_("Files with unknown extension cannot be dropped in this area. "
						"Therefore they are ignored."));
	}
	g_strfreev(uris);
	g_slist_free(list);
}

gboolean on_treeview_convert_key_release_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {
	if (event->keyval == GDK_KEY_Delete || event->keyval == GDK_KEY_KP_Delete
			|| event->keyval == GDK_KEY_BackSpace) {
		remove_selected_files_from_list();
		check_for_conversion_form_completeness();
		return TRUE;
	}
	return FALSE;
}


static void check_for_conversion_form_completeness() {
	GtkTreeIter iter;
	GtkTreeModel *model = NULL;
	gboolean valid;
	GtkWidget *go_button = lookup_widget("convert_button");

	init_tree_view();
	model = gtk_tree_view_get_model(tree_view);
	valid = gtk_tree_model_get_iter_first(model, &iter);
	gtk_widget_set_sensitive (go_button, destroot && destroot[0] != '\0' && valid);

	/* we override the sort function in order to provide natural sort order */
	gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(model),
			COLUMN_FILENAME, (GtkTreeIterCompareFunc) sort_conv_tree, NULL,
			NULL);

	update_statusbar_convert();
}

/******************Callback functions*******************************************************************/


// TODO: put a red lining around the entry instead of removing bad chars
void insert_text_handler(GtkEntry *entry, const gchar *text, gint length,
		gint *position, gpointer data) {
	GtkEditable *editable = GTK_EDITABLE(entry);
	int i, count = 0;

	gchar *result = g_strndup(text, length);

	for (i = 0; i < length; i++) {
		if (is_forbiden_in_filename(text[i]))
			continue;
		result[count++] = text[i];
	}

	if (count > 0) {
		g_signal_handlers_block_by_func(G_OBJECT (editable),
				G_CALLBACK (insert_text_handler), data);
		gtk_editable_insert_text(editable, result, count, position);
		g_signal_handlers_unblock_by_func(G_OBJECT (editable),
				G_CALLBACK (insert_text_handler), data);
	}
	g_signal_stop_emission_by_name(G_OBJECT(editable), "insert_text");

	g_free(result);
}

void update_statusbar_convert() {
	GtkLabel *status_label = GTK_LABEL(lookup_widget("statuslabel_convert"));

	int nb_files = count_converted_files();
	if (nb_files == 0)
		gtk_label_set_text(status_label, " ");
	else {
		int selected = count_selected_files();
		gchar *str, *total;
		if (nb_files == 1) {
			str = g_strdup_printf(_("%d file loaded"), nb_files);
		} else {
			str = g_strdup_printf(_("%d files loaded"), nb_files);
		}
		if (selected == 0) {
			total = g_strdup(str);
		} else if (selected == 1) {
			total = g_strdup_printf(_("%d file selected, %s"), selected, str);
		} else {
			total = g_strdup_printf(_("%d files selected, %s"), selected, str);
		}
		gtk_label_set_text(status_label, total);
		g_free(str);
		g_free(total);
	}
}

void on_treeview_selection5_changed(GtkTreeSelection *treeselection,
		gpointer user_data) {
	update_statusbar_convert();
}

// truncates destroot if it's more than 120 characters, append a '_' if it
// doesn't end with one or a '-'. SER extensions are accepted and unmodified.
void on_convtoroot_changed(GtkEditable *editable, gpointer user_data){
	static GtkWidget *multiple_ser = NULL;
	const gchar *name = gtk_entry_get_text(GTK_ENTRY(editable));

	if (*name != 0) {
		if (!multiple_ser)
			multiple_ser = lookup_widget("multipleSER");
		if (destroot)
			g_free(destroot);

		destroot = g_str_to_ascii(name, NULL); // we want to avoid special char

		const char *ext = get_filename_ext(destroot);
		if (ext && !g_ascii_strcasecmp(ext, "ser")) {
			convflags |= CONVDSTSER;
			gtk_widget_set_visible(multiple_ser, TRUE);
			char *base = remove_ext_from_filename(destroot);
			if (check_if_seq_exist(base)) {
				set_icon_entry(GTK_ENTRY(editable), "gtk-dialog-warning");
			} else{
				set_icon_entry(GTK_ENTRY(editable), NULL);
			}
			free(base);
		} else {
			convflags &= ~CONVDSTSER;
			gtk_widget_set_visible(multiple_ser, FALSE);
			destroot = format_basename(destroot);
			if (check_if_seq_exist(destroot)) {
				set_icon_entry(GTK_ENTRY(editable), "gtk-dialog-warning");
			} else{
				set_icon_entry(GTK_ENTRY(editable), NULL);
			}
		}

	} else {
		set_icon_entry(GTK_ENTRY(editable), NULL);
		g_free(destroot);
		destroot = NULL;
	}
	check_for_conversion_form_completeness();
}

void check_debayer_button (GtkToggleButton *togglebutton) {
	static GtkToggleButton *but = NULL;
	if (!but) but = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton1"));
	if (gtk_toggle_button_get_active(togglebutton)) {
		set_debayer_in_convflags();
		gtk_toggle_button_set_active(but, TRUE);
		com.debayer.open_debayer = TRUE;
	}
	else {
		unset_debayer_in_convflags();	// used for conversion
		com.debayer.open_debayer = FALSE;	// used for image opening
	}
}

void on_demosaicing_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	check_debayer_button(togglebutton);
}

void on_multipleSER_toggled (GtkToggleButton *togglebutton, gpointer user_data) {
	if (gtk_toggle_button_get_active(togglebutton))
		convflags |= CONVMULTIPLE;
	else convflags &= ~CONVMULTIPLE;
}

void on_conv3planefit_toggled (GtkToggleButton *togglebutton, gpointer user_data){
	convflags |= CONV1X3;
	convflags &= ~(CONV3X1|CONV1X1);
}

void on_conv3_1plane_toggled (GtkToggleButton *togglebutton, gpointer user_data) {
	convflags |= CONV3X1;
	convflags &= ~(CONV1X1|CONV1X3);
}

void on_conv1_1plane_toggled (GtkToggleButton *togglebutton, gpointer user_data) {
	convflags |= CONV1X1;
	convflags &= ~(CONV3X1|CONV1X3);
}
