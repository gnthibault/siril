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

#include <gdk/gdkkeysyms.h>
#include <gtk/gtk.h>
#include <string.h>

#include "gui/conversion.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "io/conversion.h"
#include "io/sequence.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "algos/sorting.h"

static gchar *destroot = NULL;
static GtkListStore *liststore_convert = NULL;
static GtkTreeView *tree_view = NULL;
static GtkTreeModel *model = NULL;
static gboolean warning_is_displayed = FALSE;

static void check_for_conversion_form_completeness();
static void on_input_files_change();
static sequence_type get_activated_output_type();

static void init_widgets() {
	if (!tree_view) {
		tree_view = GTK_TREE_VIEW(gtk_builder_get_object(builder, "treeview_convert"));
		model = gtk_tree_view_get_model(tree_view);
		liststore_convert = GTK_LIST_STORE(
				gtk_builder_get_object(builder, "liststore_convert"));
	}
	g_assert(tree_view);
	g_assert(model);
	g_assert(liststore_convert);
}

int count_converted_files() {
	init_widgets();
	GtkTreeIter iter;
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
	init_widgets();
	GtkTreeSelection *selection = gtk_tree_view_get_selection(tree_view);
	return gtk_tree_selection_count_selected_rows(selection);
}

static void initialize_convert() {
	gchar *file_data, *file_date;
	GtkTreeIter iter;
	GList *list = NULL;
	
	init_widgets();
	
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
		else if (type == TYPEUNDEF) {
			char *title = siril_log_message(_("Filetype is not supported, cannot convert: %s\n"), src_ext);
			gchar *msg = g_strdup_printf(_("File extension '%s' is not supported.\n"
				"Verify that you typed the extension correctly.\n"
				"If so, you may need to install third-party software to enable "
				"this file type conversion, look at the README file.\n"
				"If the file type you are trying to load is listed in supported "
				"formats, you may notify the developers that the extension you are "
				"trying to use should be recognized for this type."), src_ext);
			siril_message_dialog(GTK_MESSAGE_ERROR, title, msg);
			return;

		}
		else there_is_an_image = TRUE;
		//if (count == 0)
		//	first_type = type;
		//if (type != first_type)
		//	several_type_of_files = TRUE;

		valid = gtk_tree_model_iter_next(model, &iter);
		count++;
	}

	sequence_type output_type = get_activated_output_type();
	int nb_allowed;
	if (!allow_to_open_files(count, &nb_allowed) && output_type == SEQ_REGULAR) {
		gboolean confirm = siril_confirm_dialog(_("Too many files are being converted."),
				_("You are about to convert a large amount of files into standard FITS files."
						"However, your OS limits the number of files that will be processed in the same time."
						"You may want to convert your input files into a FITS sequence."));
		if (!confirm) return;
	}

	gboolean multiple, debayer, symbolic_link;
	GtkToggleButton *toggle = GTK_TOGGLE_BUTTON(lookup_widget("multipleSER"));
	multiple = gtk_toggle_button_get_active(toggle);
	toggle = GTK_TOGGLE_BUTTON(lookup_widget("demosaicingButton"));
	debayer = gtk_toggle_button_get_active(toggle);
	toggle = GTK_TOGGLE_BUTTON(lookup_widget("convert_symlink"));
	symbolic_link = gtk_toggle_button_get_active(toggle);

	/* handle impossible cases */
	/* why is it forbidden?
	 * apparently SER cannot be converted to SER, wouldn't it be nice to debayer them? */
	if (debayer && there_is_a_ser) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("A conflict has been detected."),
				_("The Debayer option is not allowed in SER conversion, please uncheck the option."));
		g_list_free_full(list, g_free);
		return;
	}
	if (output_type == SEQ_REGULAR && debayer && symbolic_link) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("A conflict has been detected."),
				_("The symbolic link option is not allowed with the debayer one, please uncheck one option."));
		g_list_free_full(list, g_free);
		return;
	}
	if (multiple && there_is_a_ser) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("A conflict has been detected."),
				_("The Multiple SER option is not allowed in SER conversion, please uncheck the option."));
		g_list_free_full(list, g_free);
		return;
	}
	if (multiple && there_is_an_image) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("A conflict has been detected."),
				_("Creating multiple SER can only be done with only films as input."));
		g_list_free_full(list, g_free);
		return;
	}

	siril_log_color_message(_("Conversion: processing %d files...\n"), "green", count);
	
	set_cursor_waiting(TRUE);
	control_window_switch_to_tab(OUTPUT_LOGS);

	/* g_list_append() has to traverse the entire list to find the end,
	 * which is inefficient when adding multiple elements. A common idiom
	 * to avoid the inefficiency is to use g_list_prepend() and reverse the
	 * list with g_list_reverse() when all elements have been added. */
	list = g_list_reverse(list);
	/* convert the list to an array for parallel processing */
	char **files_to_convert = glist_to_array(list, &count);

	struct _convert_data *args = malloc(sizeof(struct _convert_data));
	if (!args) {
		PRINT_ALLOC_ERR;
		return;
	}
	if (output_type == SEQ_REGULAR) {
		GtkEntry *startEntry = GTK_ENTRY(lookup_widget("startIndiceEntry"));
		const gchar *index = gtk_entry_get_text(startEntry);
		args->start = (g_ascii_strtoll(index, NULL, 10) <= 0
						|| g_ascii_strtoll(index, NULL, 10) >= INDEX_MAX) ?	1 : g_ascii_strtoll(index, NULL, 10);
	}
	else args->start = 0;
	args->list = files_to_convert;
	args->total = count;
	args->nb_converted_files = 0;
	args->input_has_a_seq = !no_sequence_to_convert;
	args->destroot = g_strdup(destroot);
	args->debayer = debayer;
	args->make_link = symbolic_link;
	args->output_type = output_type;
	args->multiple_output = multiple;
	gettimeofday(&(args->t_start), NULL);
	start_in_new_thread(convert_thread_worker, args);
	return;
}

void on_convroot_entry_activate(GtkEntry *entry, gpointer user_data) {
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
	GList *references, *list;

	init_widgets();
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
	init_widgets();

	for (l = list; l; l = l->next) {
		char *filename;

		filename = (char *) l->data;
		if (g_stat(filename, &st) == 0) {
			add_convert_to_list(filename, st);
		}
		g_free(filename);
	}
	check_for_conversion_form_completeness();
	on_input_files_change();
}

void on_clear_convert_button_clicked(GtkToolButton *button, gpointer user_data) {
	init_widgets();
	gtk_list_store_clear(liststore_convert);
	check_for_conversion_form_completeness();
	on_input_files_change();
}

void on_remove_convert_button_clicked(GtkToolButton *button, gpointer user_data) {
	init_widgets();
	remove_selected_files_from_list();
	check_for_conversion_form_completeness();
	on_input_files_change();
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
			g_clear_error(&error);
		}
	}
	list = g_slist_sort(list, (GCompareFunc) strcompare);
	fill_convert_list(list);
	if (bad_files) {
		gchar *loc_str = ngettext("%d file was ignored while drag and drop\n",
				"%d files were ignored while drag and drop\n", bad_files);
		loc_str = g_strdup_printf(loc_str, bad_files);
		char *msg = siril_log_message(loc_str);
		siril_message_dialog(GTK_MESSAGE_INFO, msg,
				_("Files with unknown extension cannot be dropped in this area. "
						"Therefore they are ignored."));
		g_free(loc_str);
	}
	g_strfreev(uris);
	g_slist_free(list);
	on_input_files_change();
}

gboolean on_treeview_convert_key_release_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {
	if (event->keyval == GDK_KEY_Delete || event->keyval == GDK_KEY_KP_Delete
			|| event->keyval == GDK_KEY_BackSpace) {
		remove_selected_files_from_list();
		check_for_conversion_form_completeness();
		on_input_files_change();
		return TRUE;
	}
	return FALSE;
}

static void check_for_conversion_form_completeness() {
	GtkTreeIter iter;
	static GtkWidget *go_button = NULL;
	if (!go_button)
		go_button = lookup_widget("convert_button");

	init_widgets();
	gboolean valid = gtk_tree_model_get_iter_first(model, &iter);
	gtk_widget_set_sensitive(go_button, destroot && destroot[0] != '\0' && valid);
}

static void on_input_files_change() {
	/* we override the sort function in order to provide natural sort order */
	gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(model),
			COLUMN_FILENAME, (GtkTreeIterCompareFunc) sort_conv_tree, NULL,
			NULL);

	update_statusbar_convert();
}

/******************************* Callback functions ***********************************/

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
	static GtkLabel *status_label = NULL;
	if (!status_label)
		status_label = GTK_LABEL(lookup_widget("statuslabel_convert"));

	int nb_files = count_converted_files();
	if (nb_files == 0)
		gtk_label_set_text(status_label, " ");
	else {
		int selected = count_selected_files();
		gchar *str1, *total;

		str1 = ngettext("%d file loaded", "%d files loaded", nb_files);
		str1 = g_strdup_printf(str1, nb_files);
		if (selected == 0) {
			total = g_strdup(str1);
		} else {
			gchar *str2 = ngettext("%d file selected", "%d files selected", selected);
			str2 = g_strdup_printf(str2, selected);
			total = g_strdup_printf("%s, %s", str1, str2);
			g_free(str2);
		}
		gtk_label_set_text(status_label, total);
		g_free(str1);
		g_free(total);
	}
}

void on_treeview_selection5_changed(GtkTreeSelection *treeselection,
		gpointer user_data) {
	update_statusbar_convert();
}

void process_destroot(sequence_type output_type) {
	static GtkEntry *convroot_entry = NULL;
	if (!convroot_entry)
		convroot_entry = GTK_ENTRY(lookup_widget("convroot_entry"));

	const gchar *name = gtk_entry_get_text(convroot_entry);
	if (*name == '\0') {
		g_free(destroot);
		destroot = NULL;
		return;
	}
	destroot = g_str_to_ascii(name, NULL); // we want to avoid special char
	gboolean seq_exists = FALSE;
	if (output_type == SEQ_SER) {
		if (!g_str_has_suffix(destroot, ".ser"))
			str_append(&destroot, ".ser");
		seq_exists = check_if_seq_exist(destroot, FALSE);
	}
	else if (output_type == SEQ_FITSEQ) {
		if (!g_str_has_suffix(destroot, com.pref.ext))
			str_append(&destroot, com.pref.ext);
		seq_exists = check_if_seq_exist(destroot, FALSE);
	}
	else {
		destroot = format_basename(destroot, TRUE);
		seq_exists = check_if_seq_exist(destroot, TRUE);
	}

	if (seq_exists && !warning_is_displayed) {
		set_icon_entry(convroot_entry, "gtk-dialog-warning");
		warning_is_displayed = TRUE;
	}
	else if (!seq_exists && warning_is_displayed) {
		set_icon_entry(convroot_entry, NULL);
		warning_is_displayed = FALSE;
	}
}

/* 0: FITS images
 * 1: SER sequence
 * 2: FITS sequence
 */
static sequence_type get_activated_output_type() {
	static GtkComboBox *combo = NULL;
	if (!combo)
		combo = GTK_COMBO_BOX(gtk_builder_get_object(builder, "prepro_output_type_combo1"));
	return (sequence_type)gtk_combo_box_get_active(combo);
}

// truncates destroot if it's more than 120 characters, append a '_' if it
// doesn't end with one or a '-'
void on_convtoroot_changed(GtkEditable *editable, gpointer user_data){
	process_destroot(get_activated_output_type());
	check_for_conversion_form_completeness();
}

// used for global file opening
void on_demosaicing_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	com.pref.debayer.open_debayer = gtk_toggle_button_get_active(togglebutton);
}

void on_prepro_output_type_combo1_changed(GtkComboBox *combo, gpointer user_data) {
	static GtkWidget *multiple_ser = NULL, *convert_symlink = NULL;
	if (!multiple_ser) {
		multiple_ser = lookup_widget("multipleSER");
		convert_symlink = lookup_widget("convert_symlink");
	}

	sequence_type output = gtk_combo_box_get_active(combo);
	gtk_widget_set_visible(multiple_ser, output == SEQ_SER);
	gtk_widget_set_visible(convert_symlink, output == SEQ_REGULAR);
	process_destroot(output);
	check_for_conversion_form_completeness();
}
