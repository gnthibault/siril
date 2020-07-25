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
#ifndef SRC_IO_FITS_SYMLINK_H_
#define SRC_IO_FITS_SYMLINK_H_

#ifdef _WIN32
#include <windows.h>
#endif

#include <glib.h>

gboolean test_if_symlink_is_ok();
gpointer symlink_thread_worker(gpointer p);
gboolean symlink_uniq_file(gchar *src_filename, gchar *dest_filename, gboolean allow_symlink);

#ifdef _WIN32
DWORD read_registre_value(LPTSTR lpKeyName, LPTSTR lpPolicyPath);
#endif

#endif /* SRC_IO_FITS_SYMLINK_H_ */
