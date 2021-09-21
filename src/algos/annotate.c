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

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_world_cs.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include "algos/siril_wcs.h"
#include "gui/image_display.h"

#include "annotate.h"

#define USER_CATALOGUE "user-catalogue.txt"

static GSList *siril_catalogue_list = NULL;
static gboolean show_catalog(int catalog);

static const gchar *cat[] = {
		"messier.txt",
		"ngc.txt",
		"ic.txt",
		"ldn.txt",
		"sh2.txt",
		"stars.txt"
};

struct _CatalogObjects {
	gchar *code;
	gdouble ra;
	gdouble dec;
	gdouble radius;
	gchar *name;
	gchar *alias;
	gint catalogue;
};

static CatalogObjects* new_catalog_object(const gchar *code, gdouble ra,
		gdouble dec, gdouble radius, const gchar *name, const gchar *alias,
		gint catalogue) {
	CatalogObjects *object = g_new(CatalogObjects, 1);

	object->code = g_strdup(code);
	object->ra = ra;
	object->dec = dec;
	object->radius = radius;
	object->name = g_strdup(name);
	object->alias = g_strdup(alias);
	object->catalogue = catalogue;
	return object;
}

static gboolean is_inside(fits *fit, double ra, double dec) {
	double x, y;

	wcs2pix(fit, ra, dec, &x, &y);
	return (x > 0 && x < fit->rx && y > 0 && y < fit->ry);
}

static gint object_compare(gconstpointer *a, gconstpointer *b) {
	const CatalogObjects *s1 = (const CatalogObjects *) a;
	const CatalogObjects *s2 = (const CatalogObjects *) b;

	if (!s1->alias) return 1;

	gchar **token = g_strsplit(s1->alias, "/", -1);
	guint nargs = g_strv_length(token);

	if (nargs == 1)
		return g_strcmp0(s1->alias, s2->code);
	else {
		for (int i = 0; i < nargs; i++) {
			if (!g_strcmp0(token[i], s2->code)) {
				return 0;
			}
		}
		return 1;
	}
}

static GSList *load_catalog(const gchar *filename, gint cat_index) {
	GFile *file;
	gchar *line;
	GSList *list = NULL;
	GError *error = NULL;

	file = g_file_new_for_path(filename);
	GInputStream *input_stream = (GInputStream *)g_file_read(file, NULL, &error);

	if (input_stream == NULL) {
		if (error != NULL) {
			g_clear_error(&error);
			siril_log_message(_("File [%s] does not exist\n"), g_file_peek_path(file));
		}
		g_object_unref(file);
		return NULL;
	}

	GDataInputStream *data_input = g_data_input_stream_new(input_stream);
	while ((line = g_data_input_stream_read_line_utf8(data_input, NULL,
				NULL, NULL))) {
		if (g_str_has_prefix (line, "Code")) {
			g_free(line);
			continue;
		}
		gchar **token = g_strsplit(line, ";", -1);
		guint nargs = g_strv_length(token);

		const gchar *code = NULL, *name = NULL, *alias = NULL;
		gdouble ra, dec, radius;

		/* mandatory tokens */
		code = token[0];
		ra = g_ascii_strtod(token[1], NULL) * 15.0;
		dec = g_strcmp0(token[2], "-") ? g_ascii_strtod(token[3], NULL) : g_ascii_strtod(token[3], NULL) * -1.0;
		radius = g_ascii_strtod(token[4], NULL) * 0.5;

		/* optional tokens */
		if (nargs > 5) {
			name = token[6];
			if (nargs > 6) {
				alias = token[7];
			}
		}

		CatalogObjects *object = new_catalog_object(code, ra, dec, radius, name, alias, cat_index);

		list = g_slist_prepend(list, (gpointer) object);

		g_strfreev(token);
		g_free(line);
	}
	list = g_slist_reverse(list);

	g_object_unref(data_input);
	g_object_unref(input_stream);
	g_object_unref(file);
	return list;
}

static void load_all_catalogues() {
	int cat_size = G_N_ELEMENTS(cat);

	for (int i = 0; i < cat_size; i++) {
		gchar *filename = g_build_filename(siril_get_system_data_dir(), "catalogue", cat[i], NULL);
		siril_catalogue_list = g_slist_concat(siril_catalogue_list, load_catalog(filename, i));

		g_free(filename);
	}
	/* load user catalogue */
	gchar *filename = g_build_filename(siril_get_config_dir(), PACKAGE, "catalogue", USER_CATALOGUE, NULL);
	siril_catalogue_list = g_slist_concat(siril_catalogue_list, load_catalog(filename, cat_size));

	g_free(filename);
}

static GSList *get_siril_catalogue_list() {
	return siril_catalogue_list;
}

static gboolean is_catalogue_loaded() {
	return siril_catalogue_list != NULL;
}

typedef struct {
	char *greek;			// Greek letter of stars
	char *latin;			// Greek letter written in Latin
} GreekLetters;

static GreekLetters convert_to_greek[] = {
        { "\u03b1", "alf" },
        { "\u03b2", "bet" },
        { "\u03b3", "gam" },
        { "\u03b4", "del" },
        { "\u03b5", "eps" },
        { "\u03b6", "zet" },
        { "\u03b7", "eta" },
        { "\u03b8", "tet" },
        { "\u03b9", "iot" },
        { "\u03ba", "kap" },
        { "\u03bb", "lam" },
        { "\u03bc", "mu." },
        { "\u03bd", "nu." },
        { "\u03be", "ksi" },
        { "\u03bf", "omi" },
        { "\u03c0", "pi." },
        { "\u03c1", "rho" },
        { "\u03c3", "sig" },
        { "\u03c4", "tau" },
        { "\u03c5", "ups" },
        { "\u03c6", "phi" },
        { "\u03c7", "chi" },
        { "\u03c8", "psi" },
        { "\u03c9", "ome" },
		{ NULL, NULL }
};

static gchar* replace_str(const gchar *s, const gchar *old, const gchar *new) {
	gchar *result;
	int i, cnt = 0;
	int newlen = strlen(new);
	int oldlen = strlen(old);

	// Counting the number of times old word
	// occur in the string
	for (i = 0; s[i] != '\0'; i++) {
		if (g_strstr_len(&s[i], -1, old) == &s[i]) {
			cnt++;

			// Jumping to index after the old word.
			i += oldlen - 1;
		}
	}

	// Making new string of enough length
	result = malloc(i + cnt * (newlen - oldlen) + 1);

	i = 0;
	while (*s) {
		// compare the substring with the result
		if (g_strstr_len(s, -1, old) == s) {
			strcpy(&result[i], new);
			i += newlen;
			s += oldlen;
		} else
			result[i++] = *s++;
	}

	result[i] = '\0';
	return result;
}

static void write_in_user_catalogue(CatalogObjects *object) {
	GError *error = NULL;
	GFile *file;
	GOutputStream *output_stream;
	gchar *root;

	/* First we test if root directory already exists */
	root = g_build_filename(siril_get_config_dir(), PACKAGE, "catalogue", NULL);

	if (!g_file_test(root, G_FILE_TEST_EXISTS)) {
		if (g_mkdir_with_parents(root, 0755) < 0) {
			siril_log_color_message(_("Cannot create output folder: %s\n"), "red", root);
			g_free(root);
			return;
		}
	}

	/* we write in the catalogue */
	file = g_file_new_build_filename(root, USER_CATALOGUE, NULL);
	g_free(root);

	output_stream = (GOutputStream *)g_file_append_to(file, G_FILE_CREATE_NONE, NULL, &error);
	if (output_stream == NULL) {
		if (error != NULL) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
		}
		g_object_unref(file);
		return;
	}
	gchar sign = object->dec < 0 ? '-' : '+';
	gchar *output_line = g_strdup_printf("%s;%lf;%c;%lf;;;;\n", object->code, object->ra / 15.0, sign, fabs(object->dec));

	g_output_stream_write_all(output_stream, output_line, strlen(output_line), NULL, NULL, NULL);

	g_free(output_line);
	g_object_unref(output_stream);
	g_object_unref(file);
}

GSList *find_objects(fits *fit) {
	if (!has_wcs(fit)) return NULL;
	GSList *targets = NULL;

	if (!is_catalogue_loaded())
		load_all_catalogues();
	GSList *list = get_siril_catalogue_list();

	for (GSList *l = list; l; l = l->next) {
		CatalogObjects *cur = (CatalogObjects *)l->data;
		if (!show_catalog(cur->catalogue)) continue;

		/* Search for objects in the image */
		if (is_inside(fit, cur->ra, cur->dec)) {
			if (!g_slist_find_custom(targets, cur, (GCompareFunc) object_compare)) {
				CatalogObjects *new_object = new_catalog_object(cur->code, cur->ra, cur->dec, cur->radius, cur->name, cur->alias, cur->catalogue);
				targets = g_slist_prepend(targets, new_object);
			}
		}
	}

	if (targets) {
		targets = g_slist_reverse(targets);
	}
	return targets;
}

void add_object_in_catalogue(gchar *code, SirilWorldCS *wcs) {
	int cat_size = G_N_ELEMENTS(cat);

	if (!is_catalogue_loaded())
		load_all_catalogues();

	CatalogObjects *new_object = new_catalog_object(code,
			siril_world_cs_get_alpha(wcs), siril_world_cs_get_delta(wcs), 0,
			NULL, NULL, cat_size);

	siril_catalogue_list = g_slist_append(siril_catalogue_list, new_object);
	write_in_user_catalogue(new_object);
}

gchar *get_catalogue_object_code(CatalogObjects *object) {
	gboolean found = FALSE;
	int i = 0;

	/* in case of stars we want to convert to Greek letter */
	while (convert_to_greek[i].latin) {
		gchar *latin_code = g_strstr_len(object->code, -1, convert_to_greek[i].latin);
		if (latin_code) {
			found = TRUE;
			break;
		}
		i++;
	}
	if (found) {
		gchar *code = g_strdup(replace_str(object->code, convert_to_greek[i].latin, convert_to_greek[i].greek));
		g_free(object->code);
		object->code = code;
	}
	return object->code;
}

gchar *get_catalogue_object_name(CatalogObjects *object) {
	return object->name;
}

gdouble get_catalogue_object_ra(CatalogObjects *object) {
	return object->ra;
}

gdouble get_catalogue_object_dec(CatalogObjects *object) {
	return object->dec;
}

gdouble get_catalogue_object_radius(CatalogObjects *object) {
	return object->radius;
}

void free_catalogue_object(CatalogObjects *object) {
	g_free(object->code);
	g_free(object->name);
	g_free(object);
}

void force_to_refresh_catalogue_list() {
	if (has_wcs(&gfit)) {
		if (com.found_object) {
			g_slist_free_full(com.found_object, (GDestroyNotify) free_catalogue_object);
		}
		com.found_object = find_objects(&gfit);
	}
}

static gboolean show_catalog(int catalog) {
	return com.pref.catalog[catalog];
}
