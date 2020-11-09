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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "core/siril_update.h"


#define DOMAIN_NAME "https://free-astro.org/"
#define GITLAB_URL "https://gitlab.com/free-astro/siril/"
#define TITLE_TAG_STRING "<a class=\"item-title"

static const gchar* gitlab_tags = GITLAB_URL"-/tags?sort=updated_desc";
static const gchar* gitlab_raw = GITLAB_URL"raw";
static const gchar* download_url = DOMAIN_NAME"index.php?title=Siril:";

/**
 * Check if the version is a patched version.
 * patched version are named like that x.y.z.patch where patch only contains digits.
 * if patch contains alpha char it is because that's an alpha or beta version. Not a patched one.
 * @param version version to be tested
 * @return 0 if the version is not patched. The version of the patch is returned otherwise.
 */
static gint check_for_patch(gchar *version) {
	gint i = 0;

	while (version[i]) {
		if (g_ascii_isalpha(version[i])) return 0;
		i++;
	}
	return (atoi(version));
}

static version_number get_current_version_number() {
	gchar **fullVersionNumber;
	version_number version;

	fullVersionNumber = g_strsplit_set(PACKAGE_VERSION, ".-", -1);
	version.major_version = atoi(fullVersionNumber[0]);
	version.minor_version = atoi(fullVersionNumber[1]);
	version.micro_version = atoi(fullVersionNumber[2]);
	version.patched_version = (fullVersionNumber[3] == NULL) ? 0 : check_for_patch(fullVersionNumber[3]);

	g_strfreev(fullVersionNumber);

	return version;
}

static version_number get_last_version_number(gchar *buffer) {
	gchar **token;
	gchar **fullVersionNumber;
	gchar *v = NULL;
	gint i, nargs;
	version_number version = { 0 };

	token = g_strsplit(buffer, "\n", -1);
	nargs = g_strv_length(token);

	i = 0;
	while (i < nargs) {
		if (g_str_has_prefix(g_strchug(token[i]), TITLE_TAG_STRING)) {
			/* get value between > and < */
			/* the first version number is the newer */
			strtok(token[i], ">");
			v = strtok(NULL, "<");
			break;
		}
		i++;
	}
	if (v) {
		g_fprintf(stdout, "last tagged version: %s\n", v);
		fullVersionNumber = g_strsplit_set(v, ".-", -1);

		version.major_version = atoi(fullVersionNumber[0]);
		version.minor_version = atoi(fullVersionNumber[1]);
		version.micro_version = atoi(fullVersionNumber[2]);
		version.patched_version = (fullVersionNumber[3] == NULL) ? 0 : atoi(fullVersionNumber[3]);

		g_strfreev(fullVersionNumber);
		g_strfreev(token);
	}
	return version;
}

/**
 * This function compare x1.y1.z1.patch1 vs x2.y2.z2.patch2
 * @param v1 First version number to be tested
 * @param v2 Second version number to be tested
 * @return -1 if v1 < v2, 1 if v1 > v2 and 0 if v1 is equal to v2
 */
static int compare_version(version_number v1, version_number v2) {
	if (v1.major_version < v2.major_version)
		return -1;
	else if (v1.major_version > v2.major_version)
		return 1;
	else {
		if (v1.minor_version < v2.minor_version)
			return -1;
		else if (v1.minor_version > v2.minor_version)
			return 1;
		else {
			if (v1.micro_version < v2.micro_version)
				return -1;
			else if (v1.micro_version > v2.micro_version)
				return 1;
			else {
				if (v1.patched_version < v2.patched_version)
					return -1;
				else if (v1.patched_version > v2.patched_version)
					return 1;
			}
		}
	}
	return 0;
}

static gchar *parse_changelog(gchar *changelog) {
	gchar **token;
	GString *strResult;
	gint nargs, i;

	token = g_strsplit(changelog, "\n", -1);
	nargs = g_strv_length(token);

	strResult = g_string_new(token[0]);
	strResult = g_string_append(strResult, "\n\n");
	/* we start at line 3 */
	i = 3;
	while (i < nargs && token[i][0] != '\0') {
		strResult = g_string_append(strResult, token[i]);
		strResult = g_string_append(strResult, "\n");
		i++;
	}
	g_strfreev(token);
	return g_string_free(strResult, FALSE);
}

static gchar *get_changelog(gint x, gint y, gint z, gint p) {
	GError *error = NULL;
	gchar *result = NULL;
	gchar *str;

	if (p != 0) {
		str = g_strdup_printf("/%d.%d.%d.%d/", x, y, z, p);
	} else {
		str = g_strdup_printf("/%d.%d.%d/", x, y, z);
	}
	GString *url = g_string_new(gitlab_raw);
	url = g_string_append(url, str);
	url = g_string_append(url, "ChangeLog");

	gchar *changelog_url = g_string_free(url, FALSE);
	GFile *file = g_file_new_for_uri(changelog_url);

	if (!g_file_load_contents(file, NULL, &result, NULL, NULL, &error)) {
		gchar *name = g_file_get_basename(file);
		printf("Error loading %s: %s\n", name, error->message);
		g_clear_error(&error);
		g_free(name);
	}

	g_free(changelog_url);
	g_free(str);
	g_object_unref(file);

	return result;
}

static gboolean end_update_idle(gpointer p) {
	char *msg = NULL;
	GtkMessageType message_type;
	gchar *changelog = NULL;
	gchar *data = NULL;
	version_number current_version, last_version_available;
	struct _update_data *args = (struct _update_data *) p;

	if (args->content) {
		last_version_available = get_last_version_number(args->content);
		current_version = get_current_version_number();
		gint x = last_version_available.major_version;
		gint y = last_version_available.minor_version;
		gint z = last_version_available.micro_version;
		gint patch = last_version_available.patched_version;
		if (compare_version(current_version, last_version_available) < 0) {
			msg = siril_log_message(_("New version is available. You can download it at "
							"<a href=\"%s%d.%d.%d\">%s%d.%d.%d</a>\n"),
					download_url, x, y, z, download_url, x, y, z);
			changelog = get_changelog(x, y, z, patch);
			data = parse_changelog(changelog);
			/* force the verbose variable */
			args->verbose = TRUE;
		} else if (compare_version(current_version, last_version_available)	> 0) {
			if (args->verbose)
				msg = siril_log_message(_("No update check: this is a development version\n"));
		} else {
			if (args->verbose)
				msg = siril_log_message(_("Siril is up to date\n"));
		}
		message_type = GTK_MESSAGE_INFO;
	} else {
		msg = args->msg;
		message_type = GTK_MESSAGE_ERROR;
	}
	set_cursor_waiting(FALSE);
	if (msg && args->verbose)
		siril_data_dialog(message_type, _("Software Update"), msg, data);

	/* free data */
	g_free(args->content);
	g_free(data);
	g_free(changelog);
	free(args);
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	stop_processing_thread();
	return FALSE;
}

static gpointer fetch_url(gpointer p) {
	struct _update_data *args = (struct _update_data *) p;
	GFile *file = g_file_new_for_uri(args->url);
	GError *error = NULL;

	if (!g_file_load_contents(file, NULL, &args->content, NULL, NULL, &error)) {
		gchar *name = g_file_get_basename(file);
		args->msg = siril_log_message("Error loading %s: %s\n", name, error->message);
		g_clear_error(&error);
		g_free(name);
	}
	gdk_threads_add_idle(end_update_idle, args);
	g_object_unref(file);

	return NULL;
}

void siril_check_updates(gboolean verbose) {
	struct _update_data *args;

	args = malloc(sizeof(struct _update_data));
	args->url = (gchar *)gitlab_tags;
	args->content = NULL;
	args->verbose = verbose;

	set_progress_bar_data(_("Looking for updates..."), PROGRESS_NONE);
	set_cursor_waiting(TRUE);
	start_in_new_thread(fetch_url, args);
}
