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
#ifndef SRC_GUI_SIRIL_PREVIEW_H_
#define SRC_GUI_SIRIL_PREVIEW_H_

#include "core/siril.h"

typedef struct update_preview_struct {
	gboolean show_preview;
	int (*update_preview_fn)(void);
} update_image;

void copy_gfit_to_backup();
void copy_backup_to_gfit();
fits *get_preview_gfit_backup();
gboolean is_preview_active();
void clear_backup();

void siril_preview_hide();

void set_notify_block(gboolean value);
void notify_update(gpointer user_data);

#endif /* SRC_GUI_SIRIL_PREVIEW_H_ */
