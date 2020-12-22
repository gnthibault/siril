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

#include <math.h>

#include "core/siril.h"
#include "core/siril_app_dirs.h"
#include "core/siril_log.h"
#include "gui/utils.h"
#include "algos/siril_wcs.h"
#include "gui/image_display.h"

#include "annotate.h"

static GSList *siril_catalogue_list = NULL;

/* set a tolerance for "same object" test, in degree */
#define TOLERANCE 20.0 / 3600.0;

static const gchar *cat[] = {
		"messier.txt",
		"ngc.txt",
		"ic.txt",
		"ldn.txt"
};

struct _CatalogObjects {
	gchar *code;
	gdouble ra;
	gdouble dec;
	gdouble radius;
	gchar *name;
};

static CatalogObjects *new_catalog_object(gchar *code, gdouble ra, gdouble dec, gdouble radius, gchar *name) {
	CatalogObjects *object = g_new(CatalogObjects, 1);
	object->code = g_strdup(code);
	object->ra = ra;
	object->dec = dec;
	object->radius = radius;
	object->name = g_strdup(name);
	return object;
}

gboolean is_inside(int circle_x, int circle_y, int rad, int x, int y) {
	// Compare radius of circle with distance
	// of its center from given point
	if ((x - circle_x) * (x - circle_x) + (y - circle_y) * (y - circle_y)
			<= rad * rad)
		return TRUE;
	else
		return FALSE;
}

static gboolean already_exist(GSList *list, double ra, double dec) {
	for (GSList *l = list; l; l = l->next) {

		gdouble cur_dec = ((CatalogObjects*) l->data)->dec;
		gdouble cur_ra = ((CatalogObjects*) l->data)->ra;

		double minDec = cur_dec - TOLERANCE;
		double maxDec = cur_dec + TOLERANCE;

		double minRa = cur_ra - TOLERANCE;
		double maxRa = cur_ra + TOLERANCE;

		/* compare */
		if (dec > minDec && dec < maxDec && ra > minRa && ra < maxRa) {
			return TRUE;
		}
	}
	return FALSE;
}

static GSList *load_catalog(const gchar *catalogue) {
	GFile *file;
	gchar *line;
	GSList *list = NULL;
	GError *error = NULL;

	file = g_file_new_build_filename(siril_get_system_data_dir(), "catalogue", catalogue, NULL);
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
		gchar **token = g_strsplit(line, "\t", -1);

		CatalogObjects *object = g_new(CatalogObjects, 1);
		object->code = g_strdup(token[0]);
		object->ra = g_ascii_strtod(token[1], NULL) * 15.0;
		object->dec = g_strcmp0(token[2], "-") ? g_ascii_strtod(token[3], NULL) : g_ascii_strtod(token[3], NULL) * -1.0;
		object->radius = g_ascii_strtod(token[4], NULL) * 0.5;
		object->name = g_strdup(token[6]);

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
	for (int i = 0; i < G_N_ELEMENTS(cat); i++) {
		siril_catalogue_list = g_slist_concat(siril_catalogue_list, load_catalog(cat[i]));
	}
}

static GSList *get_siril_catalogue_list() {
	return siril_catalogue_list;
}

static gboolean is_catalogue_loaded() {
	return siril_catalogue_list != NULL;
}

static GSList *find_objects(fits *fit) {
	if (!has_wcs()) return NULL;
	GSList *targets = NULL;
	gdouble x1, y1, x2, y2;
	double *crval;
	double resolution;

	crval = get_wcs_crval();
	resolution = get_wcs_image_resolution();

	if (crval == NULL) return NULL;
	if (crval[0] == 0.0 && crval[1] == 0.0) return NULL;
	if (resolution <= 0.0) return NULL;

	/* get center of the image */
	x1 = crval[0];
	y1 = crval[1];

	/* get radius of the fov */
	x2 = x1 + fit->rx * resolution;
	y2 = y1 + fit->ry * resolution;

	if (!is_catalogue_loaded())
		load_all_catalogues();
	GSList *list = get_siril_catalogue_list();

	for (GSList *l = list; l; l = l->next) {
		CatalogObjects *cur = (CatalogObjects *)l->data;

		/* Search for objects in the circle of radius defined by the image */
		if (is_inside(x1, y1, sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2)),
				cur->ra, cur->dec)) {
			if (!already_exist(targets, cur->ra, cur->dec)) {
				CatalogObjects *new_object = new_catalog_object(cur->code, cur->ra, cur->dec, cur->radius, cur->name);
				targets = g_slist_prepend(targets, new_object);
			}
		}
	}

	if (targets) {
		targets = g_slist_reverse(targets);
	}
	return targets;
}

gchar *get_catalogue_object_code(CatalogObjects *object) {
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

void free_object(CatalogObjects *object) {
	g_free(object->code);
	g_free(object->name);
	g_free(object);
}

/*** UI callbacks **///

void on_annotate_button_toggled(GtkToggleToolButton *togglebutton,
		gpointer user_data) {
	if (gtk_toggle_tool_button_get_active(togglebutton)) {
		if (has_wcs()) {
			com.found_object = find_objects(&gfit);
		}
	} else {
		g_slist_free_full(com.found_object, (GDestroyNotify) free_object);
		com.found_object = NULL;
	}
	redraw(com.cvport, REMAP_NONE);
}

void initialize_wcs_toggle_button() {
	gtk_toggle_tool_button_set_active(GTK_TOGGLE_TOOL_BUTTON(lookup_widget("annotate_button")), FALSE);
}
