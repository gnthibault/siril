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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <json-glib/json-glib.h>

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


#define DOMAIN "https://staging.siril.org/"
#define SIRIL_VERSIONS DOMAIN"siril_versions.json"
#define SIRIL_DOWNLOAD DOMAIN"download"
#define GITLAB_URL "https://gitlab.com/free-astro/siril/raw"


// taken from gimp
static gboolean siril_update_get_highest(JsonParser *parser,
		gchar **highest_version, gint64 *release_timestamp,
		gint *build_revision, gchar **build_comment) {
	JsonPath *path;
	JsonNode *result;
	JsonArray *versions;
	const gchar *platform;
	const gchar *path_str;
	const gchar *release_date = NULL;
	GError *error = NULL;
	gint i;

	g_return_val_if_fail(highest_version != NULL, FALSE);
	g_return_val_if_fail(release_timestamp != NULL, FALSE);
	g_return_val_if_fail(build_revision != NULL, FALSE);
	g_return_val_if_fail(build_comment != NULL, FALSE);

	*highest_version = NULL;
	*release_timestamp = 0;
	*build_revision = 0;
	*build_comment = NULL;

	path_str = "$['RELEASE'][*]";

	/* For Windows and macOS, let's look if installers are available.
	 * For other platforms, let's just look for source release.
	 */
	if (g_strcmp0(SIRIL_BUILD_PLATFORM_FAMILY, "windows") == 0
			|| g_strcmp0(SIRIL_BUILD_PLATFORM_FAMILY, "macos") == 0)
		platform = SIRIL_BUILD_PLATFORM_FAMILY;
	else
		platform = "source";

	path = json_path_new();
	/* Ideally we could just use Json path filters like this to
	 * retrieve only released binaries for a given platform:
	 * g_strdup_printf ("$['STABLE'][?(@.%s)]['version']", platform);
	 * json_array_get_string_element (result, 0);
	 * And that would be it! We'd have our last release for given
	 * platform.
	 * Unfortunately json-glib does not support filter syntax, so we
	 * end up looping through releases.
	 */
	if (!json_path_compile(path, path_str, &error)) {
		g_warning("%s: path compilation failed: %s\n", G_STRFUNC,
				error->message);
		g_clear_error(&error);
		g_object_unref(path);

		return FALSE;
	}
	result = json_path_match(path, json_parser_get_root(parser));
	if (!JSON_NODE_HOLDS_ARRAY(result)) {
		g_printerr("%s: match for \"%s\" is not a JSON array.\n",
		G_STRFUNC, path_str);
		g_object_unref(path);

		return FALSE;
	}

	versions = json_node_get_array(result);
	for (i = 0; i < (gint) json_array_get_length(versions); i++) {
		JsonObject *version;

		/* Note that we don't actually look for the highest version,
		 * but for the highest version for which a build for your
		 * platform (and optional build-id) is available.
		 *
		 * So we loop through the version list then the build array
		 * and break at first compatible release, since JSON arrays
		 * are ordered.
		 */
		version = json_array_get_object_element(versions, i);
		if (json_object_has_member(version, platform)) {
			JsonArray *builds;
			gint j;

			builds = json_object_get_array_member(version, platform);

			for (j = 0; j < (gint) json_array_get_length(builds); j++) {
				const gchar *build_id = NULL;
				JsonObject *build;

				build = json_array_get_object_element(builds, j);
				if (json_object_has_member(build, "build-id"))
					build_id = json_object_get_string_member (build, "build-id");
				if (g_strcmp0(build_id, "org.free_astro.siril") == 0
						|| g_strcmp0(platform, "source") == 0) {
					/* Release date is the build date if any set,
					 * otherwise the main version release date.
					 */
					if (json_object_has_member(build, "date"))
						release_date = json_object_get_string_member(build, "date");
					else
						release_date = json_object_get_string_member(version, "date");

					/* These are optional data. */
					if (json_object_has_member(build, "revision"))
						*build_revision = json_object_get_int_member(build, "revision");
					if (json_object_has_member(build, "comment"))
						*build_comment = g_strdup(json_object_get_string_member(build, "comment"));
					break;
				}
			}

			if (release_date) {
				*highest_version = g_strdup(json_object_get_string_member(version, "version"));
				break;
			}
		}
	}

	if (*highest_version && *release_date) {
		GDateTime *datetime;
		gchar *str;

		str = g_strdup_printf("%s 00:00:00Z", release_date);
		datetime = g_date_time_new_from_iso8601(str, NULL);
		g_free(str);

		if (datetime) {
			*release_timestamp = g_date_time_to_unix(datetime);
			g_date_time_unref(datetime);
		} else {
			/* JSON file data bug. */
			g_printerr(
					"%s: release date for version %s not properly formatted: %s\n",
					G_STRFUNC, *highest_version, release_date);

			g_clear_pointer(highest_version, g_free);
			g_clear_pointer(build_comment, g_free);
			*build_revision = 0;
		}
	}

	json_node_unref(result);
	g_object_unref(path);

	return (*highest_version != NULL);
}

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

static version_number get_last_version_number(gchar *version_str) {
	gchar **v;
	version_number version = { 0 };

	v = g_strsplit_set(version_str, ".-", -1);

	if (v[0])
		version.major_version = g_ascii_strtoull(v[0], NULL, 10);
	if (v[0] && v[1])
		version.minor_version = g_ascii_strtoull(v[1], NULL, 10);
	if (v[0] && v[1] && v[2])
		version.micro_version = g_ascii_strtoull(v[2], NULL, 10);
	if (v[0] && v[1] && v[2] && v[3] && v[4])
		version.patched_version = g_ascii_strtoull(v[3], NULL, 10);

	g_strfreev(v);
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
	GString *url = g_string_new(GITLAB_URL);
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

static gchar *check_version(gchar *version, gboolean *verbose, gchar **data) {
	gchar *changelog = NULL;
	gchar *msg = NULL;

	version_number last_version_available = get_last_version_number(version);
	version_number current_version = get_current_version_number();
	guint x = last_version_available.major_version;
	guint y = last_version_available.minor_version;
	guint z = last_version_available.micro_version;
	guint patch = last_version_available.patched_version;
	if (x == 0 && y == 0 && z == 0) {
		if (*verbose)
			msg = siril_log_message(_("No update check: cannot fetch version file\n"));
	} else {
		if (compare_version(current_version, last_version_available) < 0) {
			msg = siril_log_message(_("New version is available. You can download it at "
							"<a href=\"%s\">%s</a>\n"),
					SIRIL_DOWNLOAD, SIRIL_DOWNLOAD);
			changelog = get_changelog(x, y, z, patch);
			*data = parse_changelog(changelog);
			/* force the verbose variable */
			*verbose = TRUE;
		} else if (compare_version(current_version, last_version_available)	> 0) {
			if (*verbose)
				msg = siril_log_message(_("No update check: this is a development version\n"));
		} else {
			if (*verbose)
				msg = siril_log_message(_("Siril is up to date\n"));
		}
		g_free(changelog);
	}
	return msg;
}

// TODO: For now, to fix this bug https://gitlab.com/free-astro/siril/-/issues/604() we need to use GIO for Windows
#if defined HAVE_LIBCURL && !defined _WIN32

struct _update_data {
	gchar *url;
	long code;
	gchar *content;
	gboolean verbose;
};


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

static gchar *check_update_version(struct _update_data *args) {
	JsonParser *parser;
	gchar *last_version = NULL;
	gchar *build_comment = NULL;
	gint64 release_timestamp = 0;
	gint build_revision = 0;
	GError *error = NULL;
	gchar *msg = NULL;
	gchar *data = NULL;
	GtkMessageType message_type = GTK_MESSAGE_ERROR;

	parser = json_parser_new();
	if (!json_parser_load_from_data(parser, args->content, -1, &error)) {
//		g_printerr("%s: parsing of %s failed: %s\n", G_STRFUNC,
//				g_file_get_uri(G_FILE(source)), error->message);
		g_clear_object(&parser);
		g_clear_error(&error);

		return NULL;
	}

	siril_update_get_highest(parser, &last_version, &release_timestamp,	&build_revision, &build_comment);

	if (last_version) {
		g_fprintf(stdout, "Last available version: %s\n", last_version);

		msg = check_version(last_version, &(args->verbose), &data);
		message_type = GTK_MESSAGE_INFO;
	} else {
		msg = siril_log_message(_("Cannot fetch version file\n"));
	}

	g_clear_pointer(&last_version, g_free);
	g_clear_pointer(&build_comment, g_free);
	g_object_unref(parser);

	return msg;
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
		msg = check_update_version(args);
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
	args->url = SIRIL_VERSIONS;
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
	gboolean verbose = GPOINTER_TO_INT(user_data);
	char *file_contents = NULL;
	gsize file_length = 0;
	GError *error = NULL;
	gchar *msg = NULL;
	gchar *data = NULL;
	GtkMessageType message_type = GTK_MESSAGE_ERROR;

	if (g_file_load_contents_finish(G_FILE(source), result, &file_contents,
			&file_length,
			NULL, &error)) {
		JsonParser *parser;
		gchar *last_version = NULL;
		gchar *build_comment = NULL;
		gint64 release_timestamp = 0;
		gint build_revision = 0;

		parser = json_parser_new();
		if (!json_parser_load_from_data(parser, file_contents, file_length,
				&error)) {
			g_printerr("%s: parsing of %s failed: %s\n", G_STRFUNC,
					g_file_get_uri(G_FILE(source)), error->message);
			g_free(file_contents);
			g_clear_object(&parser);
			g_clear_error(&error);

			return;
		}

		siril_update_get_highest(parser, &last_version, &release_timestamp,	&build_revision, &build_comment);

		if (last_version) {
			g_fprintf(stdout, "Last available version: %s\n", last_version);

			msg = check_version(last_version, &verbose, &data);
			message_type = GTK_MESSAGE_INFO;
		} else {
			msg = siril_log_message(_("Cannot fetch version file\n"));
		}

		g_clear_pointer(&last_version, g_free);
		g_clear_pointer(&build_comment, g_free);
		g_object_unref(parser);
		g_free(file_contents);
	} else {
		g_printerr("%s: loading of %s failed: %s\n", G_STRFUNC,
				g_file_get_uri(G_FILE(source)), error->message);
		g_clear_error(&error);
		msg = siril_log_message(_("Cannot fetch version file\n"));
	}
	if (verbose) {
		set_cursor_waiting(FALSE);
		if (msg) {
			siril_data_dialog(message_type, _("Software Update"), msg, data);
		}
	}
	/* free data */
	g_free(data);
}

void siril_check_updates(gboolean verbose) {
	GFile *siril_versions;

	siril_versions = g_file_new_for_uri(SIRIL_VERSIONS);

	g_file_load_contents_async(siril_versions, NULL, siril_check_updates_callback, GINT_TO_POINTER(verbose));
	g_object_unref(siril_versions);
}
#endif
