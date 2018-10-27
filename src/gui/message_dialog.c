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
#include "core/initfile.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"

static gboolean show_modal_dialog_idle(gpointer p) {
	struct siril_dialog_data *args = (struct siril_dialog_data *) p;
	GtkWidget *dialog;

	dialog = gtk_message_dialog_new(args->parent,
			GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
			GTK_MESSAGE_WARNING, GTK_BUTTONS_OK, "%s",
			args->primary_text);

	if (args->secondary_text)
		gtk_message_dialog_format_secondary_text(GTK_MESSAGE_DIALOG(dialog), "%s",
				args->secondary_text);

	gtk_dialog_run(GTK_DIALOG(dialog));

	gtk_widget_destroy(dialog);
	free(args);
	return FALSE;
}

static void strip_last_ret_char(gchar *str) {
	if (str[strlen(str) - 1] == '\n')
		str[strlen(str) - 1] = '\0';
}

/******* Public functions *****************/

void siril_message_dialog(GtkMessageType type, char *title, char *text) {
	if (com.headless)
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
	args->type = type;
	gdk_threads_add_idle(show_modal_dialog_idle, args);
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
