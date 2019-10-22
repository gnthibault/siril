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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_LIBCURL
#include <string.h>
#include <curl/curl.h>

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

static void http_cleanup() {
	curl_easy_cleanup(curl);
	curl_global_cleanup();
	curl = NULL;
}

static gchar *get_changelog(gint x, gint y, gint z, gint p) {
	struct ucontent *changelog;
	gchar *result = NULL;
	gchar str[20];
	gchar *changelog_url;
	long code;
	GString *url = g_string_new(gitlab_raw);

	changelog = g_try_malloc(sizeof(struct ucontent));
	if (changelog == NULL) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	changelog->data = g_malloc(1);
	changelog->data[0] = '\0';
	changelog->len = 0;

	if (p != 0) {
		g_snprintf(str, sizeof(str), "/%d.%d.%d.%d/", x, y, z, p);
	} else {
		g_snprintf(str, sizeof(str), "/%d.%d.%d/", x, y, z);
	}
	url = g_string_append(url, str);
	url = g_string_append(url, "ChangeLog");

	changelog_url = g_string_free(url, FALSE);
	curl_easy_setopt(curl, CURLOPT_URL, changelog_url);
	curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
	curl_easy_setopt(curl, CURLOPT_WRITEDATA, changelog);
	if (curl_easy_perform(curl) == CURLE_OK) {
		curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &code);
		if (code == 200) {
			result = g_strdup(changelog->data);
		}
	}
	if (!result)
		g_free(changelog->data);
	g_free(changelog);
	g_free(changelog_url);

	return result;
}

static gboolean end_update_idle(gpointer p) {
	gint ret;
	char *msg;
	gchar *changelog = NULL;
	gchar *data = NULL;
	version_number current_version, last_version_available;
	static GtkMessageType type[] = { GTK_MESSAGE_INFO, GTK_MESSAGE_ERROR };
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
		ret = 1;
	} else {
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
		} else if (compare_version(current_version, last_version_available) > 0) {
			msg = siril_log_message(_("No update check: this is a development version\n"));
		} else {
			msg = siril_log_message(_("Siril is up to date\n"));
		}
		ret = 0;
	}
	set_cursor_waiting(FALSE);
	siril_data_dialog(type[ret], _("Software Update"), msg, data);

	/* free data */
	g_free(args->content);
	g_free(data);
	g_free(changelog);
	free(args);
	http_cleanup();
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	stop_processing_thread();
	return FALSE;
}

static gpointer fetch_url(gpointer p) {
	struct ucontent *content;
	gchar *result;
	long code;
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

	if (curl_easy_perform(curl) == CURLE_OK) {
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
		args->code = code;
	}
	set_progress_bar_data(NULL, PROGRESS_DONE);

	if (!result)
		g_free(content->data);
	g_free(content);

	args->content = result;

	gdk_threads_add_idle(end_update_idle, args);
	return NULL;
}

void on_help_update_activate(GtkMenuItem *menuitem, gpointer user_data) {
	struct _update_data *args;

	args = malloc(sizeof(struct _update_data));
	args->url = (gchar *)gitlab_tags;
	args->code = 0L;
	args->content = NULL;

	set_progress_bar_data(_("Looking for updates..."), PROGRESS_NONE);
	set_cursor_waiting(TRUE);
	start_in_new_thread(fetch_url, args);
}

#endif
