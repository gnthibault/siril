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

#include "siril_app_dirs.h"

static const gchar *siril_share_dir = NULL;
static const gchar *siril_config_dir = NULL;
static const gchar *siril_startup_dir = NULL;
static const gchar *siril_locale_dir = NULL;

/* To set the data dir we are looking for the glade file */
static void search_for_data_dir() {
	gchar *path;

	/* First we are looking for in the package_data_dir */
#ifdef _WIN32 // in the case where the app is started with double click on seq file
	gchar *execname = g_win32_get_package_installation_directory_of_module(NULL);
	path = g_build_filename(execname, "share", PACKAGE, NULL);
	g_free(execname);
	if (g_file_test(path, G_FILE_TEST_IS_DIR)) {
		siril_share_dir = g_strdup(path);
	}
	g_free(path);
#else
	path = g_build_filename(PACKAGE_DATA_DIR, NULL);
	if (g_file_test(path, G_FILE_TEST_IS_DIR)) {
		siril_share_dir = g_strdup(path);
	}
	g_free(path);
#endif
	/* if not found we are looking for in the common dirs */
	if (siril_share_dir == NULL) {
		int i = 0;
		const gchar *const*system_data_dirs;

		system_data_dirs = g_get_system_data_dirs();

		do {
			path = g_build_filename(system_data_dirs[i], PACKAGE, NULL);
			gchar *gladefile = g_build_filename(path, GLADE_FILE, NULL);

			/* data dir is the dir when a glade file is found */
			if (g_file_test(gladefile, G_FILE_TEST_EXISTS)) {
				siril_share_dir = g_strdup(path);

				g_free(path);
				g_free(gladefile);
				break;
			}
			g_free(path);
			g_free(gladefile);

			i++;
		} while (system_data_dirs[i] != NULL);
	}
}

static void search_for_config_dir() {
	siril_config_dir = g_get_user_config_dir();
}

/** This function tries to set a startup directory. It first looks at the "Pictures" directory,
 *  then if it does not exist, the "Document" one, Finally, if it fails on some UNIX systems
 *  the dir is set to the home directory.
 *  @return a working directory path if success, NULL if error
 */

static void search_for_startup_dir() {
	const gchar *dir = NULL;
	gint i = 0;
	size_t size;

	size = sizeof(sdir) / sizeof(GUserDirectory);

	while (dir == NULL && i < size) {
		dir = g_get_user_special_dir(sdir[i]);
		i++;
	}
	/* Not every platform has a directory for these logical id */
	if (dir == NULL) {
		dir = g_get_home_dir();
	}
	if (dir)
		siril_startup_dir = g_strdup(dir);
}

/**
 * This function search for the locale dir
 * @return the locale dir
 */
static void search_for_locale_dir() {
#ifdef _WIN32
	gchar *win32_dir;

	win32_dir = g_win32_get_package_installation_directory_of_module (NULL);
	gchar *locale_dir = g_build_filename(win32_dir, "share", "locale", NULL);

	g_free(win32_dir);

	siril_locale_dir = locale_dir;
#else
	gchar *path = g_build_filename(LOCALEDIR, NULL);
		if (g_file_test(path, G_FILE_TEST_IS_DIR)) {
			siril_locale_dir = g_strdup(path);
		}
		g_free(path);
#endif
}

/** Public functions **/

const gchar* siril_get_locale_dir() {
	return siril_locale_dir;
}
const gchar* siril_get_startup_dir() {
	return siril_startup_dir;
}

const gchar* siril_get_system_data_dir() {
	return siril_share_dir;
}

const gchar* siril_get_config_dir() {
	return siril_config_dir;
}

void initialize_siril_directories() {
	search_for_locale_dir();
	search_for_startup_dir();
	search_for_config_dir();
	search_for_data_dir();
}

