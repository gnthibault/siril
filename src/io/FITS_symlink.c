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

#ifdef _WIN32
#undef _WIN32_WINNT
#define _WIN32_WINNT _WIN32_WINNT_WIN10
#include <windows.h>
#endif
#include <stdio.h>

#ifndef S_ISLNK
#define S_ISLNK(x) 0
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "io/sequence.h"
#include "io/conversion.h"
#include "io/image_format_fits.h"
#include "gui/progress_and_log.h"

#include "FITS_symlink.h"

#ifdef _WIN32
#define PATH_APPMODEUNLOCK      TEXT("SOFTWARE\\Microsoft\\Windows\\CurrentVersion\\AppModelUnlock" )
#define CLE_APPMODEUNLOCK_ADWDL TEXT("AllowDevelopmentWithoutDevLicense" )
#define CLE_APPMODEUNLOCK_AATA  TEXT("AllowAllTrustedApps" )

static DWORD ReadValue(LPTSTR lpKeyName, LPTSTR lpPolicyPath) {
	DWORD dwReturnKo = -1;
	HKEY hKey;
	LONG lResult;
	DWORD dwValue;
	DWORD dwSize = sizeof(dwValue);

	lResult = RegOpenKeyEx(HKEY_LOCAL_MACHINE, lpPolicyPath, 0, KEY_QUERY_VALUE,
			&hKey);
	if (lResult != ERROR_SUCCESS) {
		printf("RegOpenKeyEx KO\n");
		return dwReturnKo;
	}

	lResult = RegQueryValueEx(hKey, lpKeyName, 0, NULL, (LPBYTE) & dwValue,
			&dwSize);
	RegCloseKey(hKey);

	return (lResult == ERROR_SUCCESS) ? dwValue : dwReturnKo;
}
#endif

static gboolean end_symlink_idle(gpointer p) {
	struct _symlink_data *args = (struct _symlink_data *) p;
	struct timeval t_end;

	if (!args->retval && get_thread_run() && args->nb_linked_files > 1) {
		// load the sequence
		char *linked_seqname = NULL;
		linked_seqname = malloc(strlen(args->destroot) + 5);
		sprintf(linked_seqname, "%s.seq", args->destroot);
		check_seq(0);
		if (linked_seqname) {
			update_sequences_list(linked_seqname);
			free(linked_seqname);
		}
	}

	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_DONE);
	set_cursor_waiting(FALSE);
	gettimeofday(&t_end, NULL);
	show_time(args->t_start, t_end);
	stop_processing_thread();
	g_free(args->destroot);
	free(args);
	return FALSE;
}

gpointer symlink_thread_worker(gpointer p) {
	double progress = 0.0;
	struct _symlink_data *args = (struct _symlink_data *) p;
	unsigned int frame_index = 0;
#ifdef _WIN32
	gboolean allow_symlink = TRUE;
#endif
	gboolean symlink_is_ok = TRUE;

	args->nb_linked_files = 0;
	args->retval = 0;

#ifdef _WIN32
	// AllowDevelopmentWithoutDevLicense=1  and AllowAllTrustedApps = 1 if DevMode is enabled
	// AllowDevelopmentWithoutDevLicense=0  and AllowAllTrustedApps = 0 if DevMode is disabled
	DWORD cr = ReadValue(CLE_APPMODEUNLOCK_ADWDL, PATH_APPMODEUNLOCK);
	if (cr != 1 ) {
		siril_log_color_message(_("You should enable the Developer Mode in order to create symbolic links "
								"instead of simply copying files.\n"), "red");
		allow_symlink = FALSE;
		symlink_is_ok = FALSE;
	}
#endif

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(guided) \
	if(!args->input_has_a_seq && fits_is_reentrant())
	// we should run in parallel only when images are converted, not sequences
#endif
	for (int i = 0; i < args->total; i++) {
		if (args->retval || !get_thread_run()) {
			continue;
		}

		gchar *src_filename = args->list[i];
		const char *src_ext = get_filename_ext(src_filename);
		int index = args->input_has_a_seq ? frame_index : args->start + i;

		gchar *name = g_utf8_strrchr(src_filename, strlen(src_filename), G_DIR_SEPARATOR);
		gchar *msg_bar;
		if (name)
			msg_bar = g_strdup_printf(_("Making symbolic link %s..."), name + 1);
		else msg_bar = g_strdup_printf(_("Making symbolic link %s..."), src_filename);

		image_type imagetype = get_type_for_extension(src_ext);
		if (imagetype != TYPEFITS) {
			args->retval = 1;
			g_free(msg_bar);
			continue;
		} else {
			gchar *dest_filename = g_strdup_printf("%s%05d%s", args->destroot,
					index, com.pref.ext);
			/* remove symlink already existing to avoid error */
			GStatBuf dest_stat;
			if (g_lstat(dest_filename, &dest_stat) == 0) {
				g_unlink(dest_filename);
			}

#ifdef _WIN32
			if (allow_symlink) {
				wchar_t *wsrc, *wdst;

				wsrc = g_utf8_to_utf16(src_filename, -1, NULL, NULL, NULL);
				wdst = g_utf8_to_utf16(dest_filename, -1, NULL, NULL, NULL);

				if (CreateSymbolicLinkW(wdst, wsrc, SYMBOLIC_LINK_FLAG_ALLOW_UNPRIVILEGED_CREATE) == 0) {
					copy_fits_from_file(src_filename, dest_filename);
					symlink_is_ok = FALSE;
				}

				g_free(wsrc);
				g_free(wdst);
			} else {
				copy_fits_from_file(src_filename, dest_filename);
				symlink_is_ok = FALSE;
			}
#else
			if (symlink(src_filename, dest_filename) != 0) {
				copy_fits_from_file(src_filename, dest_filename);
				symlink_is_ok = FALSE;
			}
#endif

			g_free(dest_filename);
			frame_index++;
		}

#ifdef _OPENMP
#pragma omp atomic
#endif
		progress += 1.0;
		set_progress_bar_data(msg_bar, progress / ((double) args->total));
		g_free(msg_bar);
#ifdef _OPENMP
#pragma omp atomic
#endif
		args->nb_linked_files++;
	}

	g_dir_close(args->dir);
	for (int i = 0; i < args->total; i++)
		g_free(args->list[i]);
	if (args->retval)
		siril_log_message(_("%s ended with error, %d/%d input files done\n"),
				symlink_is_ok ? _("Symbolic link creation") : _("The copy of the files"), args->nb_linked_files, args->total);
	else {
		if (args->nb_linked_files == args->total)
			siril_log_message(_("%s succeeded, %d/%d input files done\n"),
					symlink_is_ok ? _("Symbolic link creation") : _("The copy of the files"), args->nb_linked_files, args->total);
		else siril_log_message(_("%s aborted, %d/%d input files done\n"),
				symlink_is_ok ? _("Symbolic link creation") : _("The copy of the files"), args->nb_linked_files, args->total);
	}
	siril_add_idle(end_symlink_idle, args);
	return NULL;
}
