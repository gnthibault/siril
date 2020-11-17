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
#ifndef SRC_CORE_OS_UTILS_H_
#define SRC_CORE_OS_UTILS_H_

#include <glib.h>

gboolean update_displayed_memory();
int test_available_space(gint64 req_size);
guint64 get_available_memory();
int get_max_memory_in_MB();
#ifdef _WIN32
gchar *get_special_folder(int csidl);
int ReconnectIO(int OpenNewConsole);
#endif
SirilWidget* siril_file_chooser_open(GtkWindow *parent,
		GtkFileChooserAction action);
SirilWidget* siril_file_chooser_add(GtkWindow *parent,
		GtkFileChooserAction action);
SirilWidget* siril_file_chooser_save(GtkWindow *parent,
		GtkFileChooserAction action);
gint siril_dialog_run(SirilWidget *widgetdialog);
void siril_widget_destroy(SirilWidget *widgetdialog);
gboolean allow_to_open_files(int nb_frames, int *nb_allowed_file);

#endif /* SRC_CORE_OS_UTILS_H_ */
