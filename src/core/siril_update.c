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

#ifdef HAVE_LIBCURL
#include <curl/curl.h>
#endif

#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "gui/utils.h"
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
static guint check_for_patch(gchar *version) {
	guint i = 0;

	while (version[i]) {
		if (g_ascii_isalpha(version[i])) return 0;
		i++;
	}
	return (g_ascii_strtoull(version, NULL, 10));
}

static version_number get_current_version_number() {
	gchar **fullVersionNumber;
	version_number version;

	fullVersionNumber = g_strsplit_set(PACKAGE_VERSION, ".-", -1);
	version.major_version = g_ascii_strtoull(fullVersionNumber[0], NULL, 10);
	version.minor_version = g_ascii_strtoull(fullVersionNumber[1], NULL, 10);
	version.micro_version = g_ascii_strtoull(fullVersionNumber[2], NULL, 10);
	version.patched_version = (fullVersionNumber[3] == NULL) ? 0 : check_for_patch(fullVersionNumber[3]);

	g_strfreev(fullVersionNumber);

	return version;
}

static version_number get_last_version_number(gchar *buffer) {
	gchar **token;
	gchar **fullVersionNumber;
	gchar *v = NULL;
	guint i, nargs;
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

		version.major_version = g_ascii_strtoull(fullVersionNumber[0], NULL, 10);
		version.minor_version = g_ascii_strtoull(fullVersionNumber[1], NULL, 10);
		version.micro_version = g_ascii_strtoull(fullVersionNumber[2], NULL, 10);
		version.patched_version = (fullVersionNumber[3] == NULL) ? 0 : g_ascii_strtoull(fullVersionNumber[3], NULL, 10);

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
	guint nargs, i;

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
		siril_log_message(_("Error loading url: %s: %s\n"), changelog_url, error->message);
		g_clear_error(&error);
	}

	g_free(changelog_url);
	g_free(str);
	g_object_unref(file);

	return result;
}

static gchar *check_version(struct _update_data *args, gchar **data) {
	gchar *changelog = NULL;
	gchar *msg = NULL;

	version_number last_version_available = get_last_version_number(args->content);
	version_number current_version = get_current_version_number();
	guint x = last_version_available.major_version;
	guint y = last_version_available.minor_version;
	guint z = last_version_available.micro_version;
	guint patch = last_version_available.patched_version;
	if (compare_version(current_version, last_version_available) < 0) {
		msg = siril_log_message(_("New version is available. You can download it at "
						"<a href=\"%s%d.%d.%d\">%s%d.%d.%d</a>\n"),
				download_url, x, y, z, download_url, x, y, z);
		changelog = get_changelog(x, y, z, patch);
		*data = parse_changelog(changelog);
		/* force the verbose variable */
		args->verbose = TRUE;
	} else if (compare_version(current_version, last_version_available)	> 0) {
		if (args->verbose)
			msg = siril_log_message(_("No update check: this is a development version\n"));
	} else {
		if (args->verbose)
			msg = siril_log_message(_("Siril is up to date\n"));
	}
	g_free(changelog);
	return msg;
}

// TODO: For now, to fix this bug https://gitlab.com/free-astro/siril/-/issues/604() we need to use GIO for Windows
#if defined HAVE_LIBCURL && !defined _WIN32
static const int DEFAULT_FETCH_RETRIES = 5;
static CURL *curl;

struct ucontent {
	gchar *data;
	size_t len;
};

static size_t cbk_curl(void *buffer, size_t size, size_t nmemb, void *userp) {
	size_t realsize = size * nmemb;
	struct ucontent *mem = (struct ucontent *) userp;

	mem->data = g_try_realloc(mem->data, mem->len + realsize + 1);

	memcpy(&(mem->data[mem->len]), buffer, realsize);
	mem->len += realsize;
	mem->data[mem->len] = 0;

	return realsize;
}

static void init() {
	if (!curl) {
		g_fprintf(stdout, "initializing CURL\n");
		curl_global_init(CURL_GLOBAL_ALL);
		curl = curl_easy_init();
	}

	if (!curl)
		exit(EXIT_FAILURE);
}

static void http_cleanup() {
	curl_easy_cleanup(curl);
	curl_global_cleanup();
	curl = NULL;
}

static gboolean end_update_idle(gpointer p) {
	char *msg = NULL;
	gchar *data = NULL;
	GtkMessageType message_type = GTK_MESSAGE_ERROR;
	struct _update_data *args = (struct _update_data *) p;

	if (args->content == NULL) {
		switch(args->code) {
		case 0:
			msg = siril_log_message(_("Unable to check updates! "
					"Please Check your network connection\n"));
			break;
		default:
			msg = siril_log_message(_("Unable to check updates! Error: %ld\n"),
					args->code);
		}
	} else {
		msg = check_version(args, &data);
		message_type = GTK_MESSAGE_INFO;
	}
	if (args->verbose) {
		set_cursor_waiting(FALSE);
		if (msg) {
			siril_data_dialog(message_type, _("Software Update"), msg, data);
		}
	}

	/* free data */
	g_free(args->content);
	g_free(data);
	free(args);
	http_cleanup();
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	stop_processing_thread();
	return FALSE;
}

static gpointer fetch_url(gpointer p) {
	struct ucontent *content;
	gchar *result;
	long code = -1L;
	int retries;
	unsigned int s;
	struct _update_data *args = (struct _update_data *) p;

	content = g_try_malloc(sizeof(struct ucontent));
	if (content == NULL) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	g_fprintf(stdout, "fetch_url(): %s\n", args->url);

	init();
	set_progress_bar_data(NULL, 0.1);

	result = NULL;

	retries = DEFAULT_FETCH_RETRIES;

	retrieve: content->data = g_malloc(1);
	content->data[0] = '\0';
	content->len = 0;

	curl_easy_setopt(curl, CURLOPT_URL, args->url);
	curl_easy_setopt(curl, CURLOPT_VERBOSE, 0L);
	curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, cbk_curl);
	curl_easy_setopt(curl, CURLOPT_WRITEDATA, content);
	curl_easy_setopt(curl, CURLOPT_USERAGENT, "siril/0.0");
	curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);

	CURLcode retval = curl_easy_perform(curl);
	if (retval == CURLE_OK) {
		if (retries == DEFAULT_FETCH_RETRIES) set_progress_bar_data(NULL, 0.4);
		curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &code);
		if (retries == DEFAULT_FETCH_RETRIES) set_progress_bar_data(NULL, 0.6);

		switch (code) {
		case 200:
			result = content->data;
			break;
		case 500:
		case 502:
		case 503:
		case 504:
			g_fprintf(stderr, "Fetch failed with code %ld for URL %s\n", code,
					args->url);

			if (retries && get_thread_run()) {
				double progress = (DEFAULT_FETCH_RETRIES - retries) / (double) DEFAULT_FETCH_RETRIES;
				progress *= 0.4;
				progress += 0.6;
				s = 2 * (DEFAULT_FETCH_RETRIES - retries) + 2;
				char *msg = siril_log_message(_("Error: %ld. Wait %us before retry\n"), code, s);
				msg[strlen(msg) - 1] = 0; /* we remove '\n' at the end */
				set_progress_bar_data(msg, progress);
				g_usleep(s * 1E6);

				g_free(content->data);
				retries--;
				goto retrieve;
			}

			break;
		default:
			g_fprintf(stderr, "Fetch failed with code %ld for URL %s\n", code,
					args->url);
		}
		g_fprintf(stderr, "Fetch succeeded with code %ld for URL %s\n", code,
				args->url);
		args->code = code;
	} else {
		siril_log_color_message(_("Cannot retrieve information from the update URL. Error: [%ld]\n"), "red", retval);
	}
	set_progress_bar_data(NULL, PROGRESS_DONE);

	if (!result)
		g_free(content->data);
	g_free(content);

	args->content = result;

	gdk_threads_add_idle(end_update_idle, args);
	return NULL;
}

void siril_check_updates(gboolean verbose) {
	struct _update_data *args;

	args = malloc(sizeof(struct _update_data));
	args->url = (gchar *)gitlab_tags;
	args->code = 0L;
	args->content = NULL;
	args->verbose = verbose;

	set_progress_bar_data(_("Looking for updates..."), PROGRESS_NONE);
	if (args->verbose)
		set_cursor_waiting(TRUE);
	start_in_new_thread(fetch_url, args);
}

#else

static void siril_check_updates_callback(GObject *source, GAsyncResult *result,
		gpointer user_data) {
	struct _update_data *args = (struct _update_data *) user_data;
	gchar *msg = NULL;
	GtkMessageType message_type = GTK_MESSAGE_ERROR;
	gchar *data = NULL;
	GError *error = NULL;

	if (g_file_load_contents_finish(G_FILE(source), result, &args->content,
			NULL, NULL, &error)) {
		msg = check_version(args, &data);
		message_type = GTK_MESSAGE_INFO;
	} else {
		gchar *uri = g_file_get_uri(G_FILE(source));
		g_printerr("%s: loading of %s failed: %s\n", G_STRFUNC,
				uri, error->message);
		msg = siril_log_message(_("Siril cannot access to %s\n"), uri);
		g_clear_error(&error);
		g_free(uri);
	}
	if (args->verbose) {
		set_cursor_waiting(FALSE);
		if (msg) {
			siril_data_dialog(message_type, _("Software Update"), msg, data);
		}
	}

	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);

	/* free data */
	g_free(args->content);
	g_free(data);
	free(args);
}

void siril_check_updates(gboolean verbose) {
	struct _update_data *args;

	args = malloc(sizeof(struct _update_data));
	args->url = (gchar *)gitlab_tags;
	args->content = NULL;
	args->verbose = verbose;

	GFile *file = g_file_new_for_uri(args->url);
	set_progress_bar_data(_("Looking for updates..."), PROGRESS_NONE);
	if (args->verbose)
		set_cursor_waiting(TRUE);

	g_file_load_contents_async(file, NULL, siril_check_updates_callback, args);

	g_object_unref(file);
}


#endif
