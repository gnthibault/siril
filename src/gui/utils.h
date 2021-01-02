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
#ifndef SRC_GUI_UTILS_H_
#define SRC_GUI_UTILS_H_

typedef enum {
	FILE_CONVERSION,
	IMAGE_SEQ,
	PRE_PROC,
	REGISTRATION,
	PLOT,
	STACKING,
	OUTPUT_LOGS
} main_tabs;


GtkWidget* lookup_widget (const gchar *widget_name);
void set_label_text_from_main_thread(const char *label_name, const char *text);
void control_window_switch_to_tab(main_tabs tab);
GtkWidget* popover_new(GtkWidget *widget, const gchar *text);
GtkWidget* popover_new_with_image(GtkWidget *widget, const gchar *text, GdkPixbuf *pixbuf);
void set_GUI_MEM(guint64 used, const gchar *label);
void set_GUI_DiskSpace(gint64 mem, const gchar *label);

#endif /* SRC_GUI_UTILS_H_ */
