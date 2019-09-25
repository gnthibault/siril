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

#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "gui/callbacks.h"

#include "message_dialog.h"

static gboolean show_modal_dialog(gpointer p) {
	struct siril_dialog_data *args = (struct siril_dialog_data *) p;
	GtkTextIter itBeg, itEnd;
	GtkWidget *dialog;

	dialog = gtk_message_dialog_new(args->parent,
			GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
			GTK_MESSAGE_WARNING, GTK_BUTTONS_OK, "%s", args->primary_text);

	if (args->secondary_text)
		gtk_message_dialog_format_secondary_markup(GTK_MESSAGE_DIALOG(dialog),
				"%s", args->secondary_text);

	if (args->data) {
		GtkTextBuffer *buffer;
		GtkWidget *swindow;
		GtkWidget *tview;

		tview = gtk_text_view_new();
		buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(tview));
		gtk_text_buffer_get_bounds(buffer, &itBeg, &itEnd);
		gtk_text_buffer_delete(buffer, &itBeg, &itEnd);
		gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(tview), FALSE);
		gtk_text_view_set_editable(GTK_TEXT_VIEW(tview), FALSE);
		gtk_text_buffer_set_text(buffer, args->data, strlen(args->data));
		gtk_widget_set_halign(tview, GTK_ALIGN_FILL);
		gtk_widget_set_valign(tview, GTK_ALIGN_FILL);
		gtk_widget_set_margin_start(tview, 6);
		gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(tview), GTK_WRAP_WORD);

		swindow = gtk_scrolled_window_new(NULL, NULL);
		gtk_container_add(GTK_CONTAINER(swindow), tview);
		gtk_box_pack_end(GTK_BOX(gtk_dialog_get_content_area(GTK_DIALOG(dialog))),
				swindow, FALSE, FALSE, 0);
		gtk_widget_set_size_request(swindow, -1, 200);

		gtk_widget_show(swindow);
		gtk_widget_show(tview);
	}

	gtk_dialog_run(GTK_DIALOG(dialog));

	gtk_widget_destroy(dialog);
	g_free(args);
	return FALSE;
}

static gchar *strip_last_ret_char(gchar *str) {
	char *pch;
	int len;

	pch = strrchr(str, '\0');
	len = pch - str;

	if (str[len - 1] == '\n') {
		str[len - 1] = '\0';
	}
	return str;
}

/******* Public functions *****************/

void siril_message_dialog(GtkMessageType type, char *title, char *text) {
	/* headless has no GUI, so no dialog; script has a GUI but calls it from another thread
	 * so it's not safe to use dialogs in the calling thread, we just ignore it for now. */
	if (com.headless || com.script)
		return;	// show_dialog usually follows a siril_log_message() call
	struct siril_dialog_data *args = g_try_malloc(sizeof(struct siril_dialog_data));

	args->parent = siril_get_active_window();
	if (!GTK_IS_WINDOW(args->parent)) {
		args->parent = GTK_WINDOW(lookup_widget("control_window"));
	}
	/* first we want to remove the '\n' at the end of the title and text
	 * if message come from siril_log_message
	 */
	strip_last_ret_char(title);
	strip_last_ret_char(text);

	args->primary_text = title;
	args->secondary_text = text;
	args->data = NULL;
	args->type = type;
	show_modal_dialog(args);
}

struct message_data {
	GtkMessageType type;
	char *title;
	char *text;
};

static gboolean siril_message_dialog_idle(gpointer p) {
	struct message_data *data = (struct message_data *) p;
	siril_message_dialog(data->type, data->title, data->text);
	free(data->title);
	free(data->text);
	free(data);
	return FALSE;
}

void queue_message_dialog(GtkMessageType type, char *title, char *text) {
	if (com.headless || com.script)
		return;	// show_dialog usually follows a siril_log_message() call
	struct message_data *data = malloc(sizeof(struct message_data));
	data->type = type;
	data->title = strdup(title);
	data->text = strdup(text);
	gdk_threads_add_idle(siril_message_dialog_idle, data);
}

void siril_data_dialog(GtkMessageType type, char *title, char *text, gchar *data) {
	if (com.headless || com.script)
		return;	// show_dialog usually follows a siril_log_message() call
	struct siril_dialog_data *args = malloc(sizeof(struct siril_dialog_data));

	args->parent = siril_get_active_window();
	if (!GTK_IS_WINDOW(args->parent)) {
		args->parent = GTK_WINDOW(lookup_widget("control_window"));
	}
	/* first we want to remove the '\n' at the end of the title and text
	 * if message come from siril_log_message
	 */
	strip_last_ret_char(title);
	strip_last_ret_char(text);

	args->primary_text = title;
	args->secondary_text = text;
	args->data = data;
	args->type = type;
	show_modal_dialog(args);
}

gboolean siril_confirm_dialog(gchar *title, gchar *msg, gboolean show_checkbutton) {
	GtkWindow *parent;
	GtkWidget *dialog, *check = NULL;
	gint res;
	gboolean ok = FALSE;

	parent = siril_get_active_window();
	if (!GTK_IS_WINDOW(parent)) {
		/* could happend if the GtkWindow has been destroyed right after the call
		 * This is the case for chooser dialog */
		parent = GTK_WINDOW(lookup_widget("control_window"));
	}

	/* first we want to remove the '\n' at the end of the title and text
	 * if message come from siril_log_message
	 */
	strip_last_ret_char(title);
	strip_last_ret_char(msg);

	dialog = gtk_message_dialog_new(parent, GTK_DIALOG_MODAL,
			GTK_MESSAGE_QUESTION, GTK_BUTTONS_OK_CANCEL, "%s", title);
	gtk_message_dialog_format_secondary_text(GTK_MESSAGE_DIALOG(dialog), "%s", msg);
	if (show_checkbutton) {
		check = gtk_check_button_new_with_mnemonic(_("_Do not show this dialog again"));
		gtk_box_pack_end(GTK_BOX(gtk_dialog_get_content_area(GTK_DIALOG (dialog))),
				check, FALSE, FALSE, 0);
		gtk_widget_set_halign(check, GTK_ALIGN_START);
		gtk_widget_set_margin_start(check, 6);
		gtk_widget_show(check);
	}
	res = gtk_dialog_run(GTK_DIALOG(dialog));
	if (res == GTK_RESPONSE_OK) {
		ok = TRUE;
		if (check) {
			com.dontShowConfirm = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(check));
			set_GUI_misc();
			writeinitfile();
		}
	}
	gtk_widget_destroy(dialog);

	return ok;
}
