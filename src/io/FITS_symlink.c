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

#ifdef _WIN32
#undef _WIN32_WINNT
#define _WIN32_WINNT _WIN32_WINNT_WIN10
#include <windows.h>
#else
#include <string.h>
#include <errno.h>
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

DWORD read_registre_value(LPTSTR lpKeyName, LPTSTR lpPolicyPath) {
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

gboolean test_if_symlink_is_ok() {
#ifdef _WIN32
	// AllowDevelopmentWithoutDevLicense=1  and AllowAllTrustedApps = 1 if DevMode is enabled
	// AllowDevelopmentWithoutDevLicense=0  and AllowAllTrustedApps = 0 if DevMode is disabled
	DWORD cr = read_registre_value(CLE_APPMODEUNLOCK_ADWDL, PATH_APPMODEUNLOCK);
	if (cr != 1 ) {
		siril_log_color_message(_("You should enable the Developer Mode in order to create symbolic links "
				"instead of simply copying files.\n"), "red");
		return FALSE;
	}
	return TRUE;
#else
	return TRUE;
#endif
}

int symlink_uniq_file(gchar *src_filename, gchar *dest_filename, gboolean allow_symlink) {
	int retval = 0;

	/* remove already existing file to avoid error */
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
			retval = copy_fits_from_file(src_filename, dest_filename);
		}

		g_free(wsrc);
		g_free(wdst);
	} else {
		retval = copy_fits_from_file(src_filename, dest_filename);
	}
#else
	static gboolean warned = FALSE;

	if (symlink(src_filename, dest_filename)) {
		char err[150];
		strerror_r(errno, err, 150);
		if (!warned) {
			siril_log_color_message(_("Symbolic link could not be made, copying the file. Error: %s\n"), "salmon", err);
			warned = TRUE;
		}
		retval = copy_fits_from_file(src_filename, dest_filename);
	}
#endif
	return retval;
}
