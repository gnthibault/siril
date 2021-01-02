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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/OS_utils.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "core/processing.h"
#include "core/command_list.h"
#include "io/sequence.h"

#include "command_line_processor.h"

static void parse_line(char *myline, int len, int *nb) {
	int i = 0, wordnb = 0;
	char string_starter = '\0';	// quotes don't split words on spaces
	word[0] = NULL;

	do {
		while (i < len && isblank(myline[i]))
			i++;
		if (myline[i] == '"' || myline[i] == '\'')
			string_starter = myline[i++];
		if (myline[i] == '\0' || myline[i] == '\n' || myline[i] == '\r')
			break;
		word[wordnb++] = myline + i;	// the beginning of the word
		word[wordnb] = NULL;		// put next word to NULL
		do {
			i++;
			if (string_starter != '\0' && myline[i] == string_starter) {
				string_starter = '\0';
				break;
			}
		} while (i < len && (!isblank(myline[i]) || string_starter != '\0')
				&& myline[i] != '\r' && myline[i] != '\n');
		if (myline[i] == '\0')	// the end of the word and line (i == len)
			break;
		myline[i++] = '\0';		// the end of the word
	} while (wordnb < MAX_COMMAND_WORDS - 1);
	*nb = wordnb;
}

static void remove_trailing_cr(char *str) {
	if (str == NULL)
		return;
	int length = strlen(str);
	if (str[length - 1] == '\r')
		str[length - 1] = '\0';
}

static int execute_command(int wordnb) {
	// search for the command in the list
	if (word[0] == NULL) return 1;
	int i = G_N_ELEMENTS(commands);
	while (g_ascii_strcasecmp(commands[--i].name, word[0])) {
		if (i == 0) {
			siril_log_message(_("Unknown command: '%s' or not implemented yet\n"), word[0]);
			return 1 ;
		}
	}

	// verify argument count
	if (wordnb - 1 < commands[i].nbarg) {
		siril_log_message(_("Usage: %s\n"), commands[i].usage);
		return 1;
	}

	// verify if command is scriptable
	if (com.script) {
		if (!commands[i].scriptable) {
			siril_log_message(_("This command cannot be used in a script: %s\n"), commands[i].name);
			return 1;
		}
	}

	// process the command
	siril_log_color_message(_("Running command: %s\n"), "salmon", word[0]);
	fprintf(stderr, "%lu: running command %s\n", time(NULL), word[0]);
	return commands[i].process(wordnb);
}

static void update_log_icon(gboolean is_running) {
	GtkImage *image = GTK_IMAGE(lookup_widget("image_log"));
	if (is_running)
		gtk_image_set_from_icon_name(image, "gtk-yes", GTK_ICON_SIZE_LARGE_TOOLBAR);
	else
		gtk_image_set_from_icon_name(image, "gtk-no", GTK_ICON_SIZE_LARGE_TOOLBAR);
}

struct log_status_bar_idle_data {
	gchar *myline;
	int line;
};

static gboolean log_status_bar_idle_callback(gpointer p) {
	struct log_status_bar_idle_data *data = (struct log_status_bar_idle_data *) p;

	GtkStatusbar *statusbar_script = GTK_STATUSBAR(lookup_widget("statusbar_script"));
	gchar *status;
	gchar *newline;

	update_log_icon(TRUE);

	newline = g_strdup(data->myline);
	status = g_strdup_printf(_("Processing line %d: %s"), data->line, newline);

	gtk_statusbar_push(statusbar_script, 0, status);
	g_free(newline);
	g_free(status);
	g_free(data->myline);
	free(data);

	return FALSE;	// only run once
}

static void display_command_on_status_bar(int line, char *myline) {
	if (!com.headless) {
		struct log_status_bar_idle_data *data;

		data = malloc(sizeof(struct log_status_bar_idle_data));
		data->line = line;
		data->myline = myline ? g_strdup(myline) : NULL;
		gdk_threads_add_idle(log_status_bar_idle_callback, data);
	}
}

static void clear_status_bar() {
	if (!com.headless) {
		GtkStatusbar *bar = GTK_STATUSBAR(lookup_widget("statusbar_script"));
		gtk_statusbar_remove_all(bar, 0);
		update_log_icon(FALSE);
	}
}

static gboolean end_script(gpointer p) {
	clear_status_bar();
	set_GUI_CWD();
	
	set_cursor_waiting(FALSE);
	return FALSE;
}

gpointer execute_script(gpointer p) {
	GInputStream *input_stream = (GInputStream*) p;
	gboolean check_required = FALSE;
	gchar *buffer;
	int line = 0, retval = 0;
	int wordnb;
	int startmem, endmem;
	struct timeval t_start, t_end;

	com.script = TRUE;
	com.stop_script = FALSE;

	gettimeofday(&t_start, NULL);

	/* Now we want to save the cwd in order to come back after
	 * script execution
	 */
	gchar *saved_cwd = g_strdup(com.wd);
	startmem = get_available_memory() / BYTES_IN_A_MB;
	gsize length = 0;
	GDataInputStream *data_input = g_data_input_stream_new(input_stream);
	while ((buffer = g_data_input_stream_read_line_utf8(data_input, &length,
			NULL, NULL))) {
		++line;
		if (com.stop_script) {
			retval = 1;
			g_free (buffer);
			break;
		}
		/* Displays comments */
		if (buffer[0] == '#') {
			siril_log_color_message("%s\n", "blue", buffer);
			g_free (buffer);
			continue;
		}
		if (buffer[0] == '\0') {
			g_free (buffer);
			continue;
		}

		/* in Windows case, remove trailing CR */
		remove_trailing_cr(buffer);

		display_command_on_status_bar(line, buffer);
		parse_line(buffer, length, &wordnb);
		/* check for requires command */
		if (!g_ascii_strcasecmp(word[0], "requires")) {
			check_required = TRUE;
		} else {
			if (com.pref.script_check_requires && check_required == FALSE) {
				siril_log_color_message(_("The \"requires\" command is missing at the top of the script file."
						" This command is needed to check script compatibility.\n"), "red");
				retval = 1;
				g_free (buffer);
				break;
			}
		}
		if ((retval = execute_command(wordnb))) {
			siril_log_message(_("Error in line %d: '%s'.\n"), line, buffer);
			siril_log_message(_("Exiting batch processing.\n"));
			g_free (buffer);
			break;
		}
		if (waiting_for_thread()) {
			retval = 1;
			g_free (buffer);
			break;	// abort script on command failure
		}
		endmem = get_available_memory() / BYTES_IN_A_MB;
		siril_debug_print("End of command %s, memory difference: %d MB\n", word[0], startmem - endmem);
		startmem = endmem;
		memset(word, 0, sizeof word);
		g_free (buffer);
	}
	g_object_unref(data_input);
	g_object_unref(input_stream);
	com.script = FALSE;
	/* Now we want to restore the saved cwd */
	siril_change_dir(saved_cwd, NULL);
	writeinitfile();
	siril_add_idle(end_script, NULL);
	if (!retval) {
		siril_log_message(_("Script execution finished successfully.\n"));
		gettimeofday(&t_end, NULL);
		show_time_msg(t_start, t_end, _("Total execution time"));
	} else {
		char *msg = siril_log_message(_("Script execution failed.\n"));
		msg[strlen(msg) - 1] = '\0';
		set_progress_bar_data(msg, PROGRESS_DONE);
	}
	g_free(saved_cwd);
	fprintf(stderr, "Script thread exiting\n");
	return GINT_TO_POINTER(retval);
}

static void show_command_help_popup(GtkEntry *entry) {
	gchar *helper = NULL;

	const gchar *text = gtk_entry_get_text(entry);
	if (*text == '\0') {
		helper = g_strdup(_("Please enter an existing command before hitting this button"));
	} else {
		command *current = commands;
		gchar **command_line = g_strsplit_set(text, " ", -1);
		while (current->process) {
			if (!g_ascii_strcasecmp(current->name, command_line[0])) {
				gchar **token;

				token = g_strsplit_set(current->usage, " ", -1);
				GString *str = g_string_new(token[0]);
				str = g_string_prepend(str, "<span foreground=\"red\"><b>");
				str = g_string_append(str, "</b>");
				if (token[1] != NULL) {
					str = g_string_append(str,
							current->usage + strlen(token[0]));
				}
				str = g_string_append(str, "</span>\n\n\t");
				str = g_string_append(str, _(current->definition));
				str = g_string_append(str, "\n\n<b>");
				str = g_string_append(str, _("Can be used in a script: "));
				str = g_string_append(str, "<span foreground=\"red\">");
				if (current->scriptable) {
					str = g_string_append(str, _("YES"));
				} else {
					str = g_string_append(str, _("NO"));
				}
				str = g_string_append(str, "</span></b>");
				helper = g_string_free(str, FALSE);
				g_strfreev(token);
				break;
			}
			current++;
		}
		g_strfreev(command_line);
	}
	if (!helper) {
		helper = g_strdup(_("No help for this command"));
	}

	GtkWidget *popover = popover_new(lookup_widget("command"), helper);
#if GTK_CHECK_VERSION(3, 22, 0)
	gtk_popover_popup(GTK_POPOVER(popover));
#else
	gtk_widget_show(popover);
#endif
	g_free(helper);
}

int processcommand(const char *line) {
	int wordnb = 0;
	gchar *myline;
	GError *error = NULL;

	if (line[0] == '\0' || line[0] == '\n')
		return 0;
	if (line[0] == '@') { // case of files
		if (get_thread_run()) {
			PRINT_ANOTHER_THREAD_RUNNING;
			return 1;
		}
		if (com.script_thread)
			g_thread_join(com.script_thread);

		/* Switch to console tab */
		control_window_switch_to_tab(OUTPUT_LOGS);

		char filename[256];
		g_strlcpy(filename, line + 1, 250);
		expand_home_in_filename(filename, 256);

		GFile *file = g_file_new_for_path(filename);
		GInputStream *input_stream = (GInputStream *)g_file_read(file, NULL, &error);

		if (input_stream == NULL) {
			if (error != NULL) {
				g_clear_error(&error);
				siril_log_message(_("File [%s] does not exist\n"), filename);
			}

			g_object_unref(file);
			return 1;
		}
		/* ensure that everything is closed */
		process_close(0);
		/* Then, run script */
		siril_log_message(_("Starting script %s\n"), filename);
		com.script_thread = g_thread_new("script", execute_script, input_stream);
		g_object_unref(file);
	} else {
		/* Switch to console tab */
		control_window_switch_to_tab(OUTPUT_LOGS);

		myline = strdup(line);
		int len = strlen(line);
		parse_line(myline, len, &wordnb);
		if (execute_command(wordnb)) {
			siril_log_color_message(_("Command execution failed.\n"), "red");
			if (!com.script && !com.headless) {
				show_command_help_popup(GTK_ENTRY(lookup_widget("command")));
			}
			free(myline);
			return 1;
		}
		free(myline);
	}
	set_cursor_waiting(FALSE);
	return 0;
}

// loads the sequence from com.wd
sequence *load_sequence(const char *name, char **get_filename) {
	gchar *file = g_strdup(name);
	gchar *altfile = NULL;
	if (!g_str_has_suffix(name, ".seq")) {
		str_append(&file, ".seq");
		if (!g_str_has_suffix(name, "_"))
			altfile = g_strdup_printf("%s_.seq", name);
	}

	if (!is_readable_file(file) && (!altfile || !is_readable_file(altfile))) {
		if (check_seq(FALSE)) {
			siril_log_message(_("No sequence `%s' found.\n"), name);
			g_free(file);
			g_free(altfile);
			return NULL;
		}
	}

	sequence *seq;
	if ((seq = readseqfile(file))) {
		if (get_filename) {
			*get_filename = file;
			file = NULL; // do not free
		}
	}
	else if (altfile && (seq = readseqfile(altfile))) {
		if (get_filename) {
			*get_filename = altfile;
			altfile = NULL; // do not free
		}
	}
	if (!seq)
		siril_log_message(_("Loading sequence `%s' failed.\n"), name);
	else {
		if (seq_check_basic_data(seq, FALSE) == -1) {
			free(seq);
			seq = NULL;
		}
	}
	g_free(file);
	g_free(altfile);
	return seq;
}

/* callback functions */

#define COMPLETION_COLUMN 0

static gboolean on_match_selected(GtkEntryCompletion *widget, GtkTreeModel *model,
		GtkTreeIter *iter, gpointer user_data) {
	const gchar *cmd;
	GtkEditable *e = (GtkEditable *) gtk_entry_completion_get_entry(widget);
	gchar *s = gtk_editable_get_chars(e, 0, -1);
	gint cur_pos = gtk_editable_get_position(e);
	gint p = cur_pos;
	gchar *end;
	gint del_end_pos = -1;

	gtk_tree_model_get(model, iter, COMPLETION_COLUMN, &cmd, -1);

	end = s + cur_pos;

	if (end) {
		del_end_pos = end - s + 1;
	} else {
		del_end_pos = cur_pos;
	}

	gtk_editable_delete_text(e, 0, del_end_pos);
	gtk_editable_insert_text(e, cmd, -1, &p);
	gtk_editable_set_position(e, p);

	return TRUE;
}

static gboolean completion_match_func(GtkEntryCompletion *completion,
		const gchar *key, GtkTreeIter *iter, gpointer user_data) {
	gboolean res = FALSE;
	char *tag = NULL;

	if (*key == '\0') return FALSE;

	GtkTreeModel *model = gtk_entry_completion_get_model(completion);
	int column = gtk_entry_completion_get_text_column(completion);

	if (gtk_tree_model_get_column_type(model, column) != G_TYPE_STRING)
		return FALSE;

	gtk_tree_model_get(model, iter, column, &tag, -1);

	if (tag) {
		char *normalized = g_utf8_normalize(tag, -1, G_NORMALIZE_ALL);
		if (normalized) {
			char *casefold = g_utf8_casefold(normalized, -1);
			if (casefold) {
				res = g_strstr_len(casefold, -1, key) != NULL;
			}
			g_free(casefold);
		}
		g_free(normalized);
		g_free(tag);
	}

	return res;
}

void init_completion_command() {
	GtkEntryCompletion *completion = gtk_entry_completion_new();
	GtkListStore *model = gtk_list_store_new(1, G_TYPE_STRING);
	GtkTreeIter iter;
	GtkEntry *entry = GTK_ENTRY(lookup_widget("command"));

	gtk_entry_completion_set_text_column(completion, COMPLETION_COLUMN);
	gtk_entry_set_completion(entry, completion);
	gtk_entry_completion_set_inline_completion(completion, TRUE);
	gtk_entry_completion_set_popup_single_match(completion, FALSE);
	gtk_entry_completion_set_minimum_key_length(completion, 2);
	gtk_entry_completion_set_match_func(completion, completion_match_func, NULL, NULL);
	g_signal_connect(G_OBJECT(completion), "match-selected", G_CALLBACK(on_match_selected), NULL);

	/* Populate the completion database. */
	command *current = commands;

	while (current->process){
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter, COMPLETION_COLUMN, current->name, -1);
		current++;
	}
	gtk_entry_completion_set_model(completion, GTK_TREE_MODEL(model));
	g_object_unref(model);
}

void on_GtkCommandHelper_clicked(GtkButton *button, gpointer user_data) {
	show_command_help_popup((GtkEntry *)user_data);
}

/** Callbacks **/

/*
 * Command line history static function
 */

static void history_add_line(char *line) {
	if (!com.cmd_history) {
		com.cmd_hist_size = CMD_HISTORY_SIZE;
		com.cmd_history = calloc(com.cmd_hist_size, sizeof(const char*));
		com.cmd_hist_current = 0;
		com.cmd_hist_display = 0;
	}
	com.cmd_history[com.cmd_hist_current] = line;
	com.cmd_hist_current++;
	// circle at the end
	if (com.cmd_hist_current == com.cmd_hist_size)
		com.cmd_hist_current = 0;
	if (com.cmd_history[com.cmd_hist_current]) {
		free(com.cmd_history[com.cmd_hist_current]);
		com.cmd_history[com.cmd_hist_current] = NULL;
	}
	com.cmd_hist_display = com.cmd_hist_current;
}

/* handler for the single-line console */
gboolean on_command_key_press_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {
	const gchar *text;
	int handled = 0;
	static GtkEntry *entry = NULL;
	if (!entry)
		entry = GTK_ENTRY(widget);
	GtkEditable *editable = GTK_EDITABLE(entry);
	int entrylength = 0;

	switch (event->keyval) {
	case GDK_KEY_Return:
	case GDK_KEY_KP_Enter:
		handled = 1;
		text = gtk_entry_get_text(entry);
		history_add_line(strdup(text));
		if (!(processcommand(text))) {
			gtk_entry_set_text(entry, "");
			set_precision_switch();
		}
		break;
	case GDK_KEY_Up:
		handled = 1;
		if (!com.cmd_history)
			break;
		if (com.cmd_hist_display > 0) {
			if (com.cmd_history[com.cmd_hist_display - 1])
				--com.cmd_hist_display;
			// display previous entry
			gtk_entry_set_text(entry, com.cmd_history[com.cmd_hist_display]);
		} else if (com.cmd_history[com.cmd_hist_size - 1]) {
			// ring back, display previous
			com.cmd_hist_display = com.cmd_hist_size - 1;
			gtk_entry_set_text(entry, com.cmd_history[com.cmd_hist_display]);
		}
		entrylength = gtk_entry_get_text_length(entry);
		gtk_editable_set_position(editable, entrylength);
		break;
	case GDK_KEY_Down:
		handled = 1;
		if (!com.cmd_history)
			break;
		if (com.cmd_hist_display == com.cmd_hist_current)
			break;
		if (com.cmd_hist_display == com.cmd_hist_size - 1) {
			if (com.cmd_hist_current == 0) {
				// ring forward, end
				gtk_entry_set_text(entry, "");
				com.cmd_hist_display++;
			} else if (com.cmd_history[0]) {
				// ring forward, display next
				com.cmd_hist_display = 0;
				gtk_entry_set_text(entry, com.cmd_history[0]);
			}
		} else {
			if (com.cmd_hist_display == com.cmd_hist_current - 1) {
				// end
				gtk_entry_set_text(entry, "");
				com.cmd_hist_display++;
			} else if (com.cmd_history[com.cmd_hist_display + 1]) {
				// display next
				gtk_entry_set_text(entry,
						com.cmd_history[++com.cmd_hist_display]);
			}
		}
		entrylength = gtk_entry_get_text_length(entry);
		gtk_editable_set_position(editable, entrylength);
		break;
	case GDK_KEY_Page_Up:
	case GDK_KEY_Page_Down:
		handled = 1;
		// go to first and last in history
		break;
	}
	return (handled == 1);
}
