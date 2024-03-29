/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2017 team free-astro (see more in AUTHORS file)
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#ifdef HAVE_LIBCURL
#include <curl/curl.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/sleef.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "core/siril_world_cs.h"
#include "core/undo.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "gui/photometric_cc.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/PSF_list.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/annotate.h"
#include "algos/siril_wcs.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "opencv/opencv.h"
#include "registration/registration.h"
#include "registration/matching/match.h"
#include "registration/matching/apply_match.h"
#include "registration/matching/misc.h"
#include "registration/matching/atpmatch.h"
#include "registration/matching/project_coords.h"

#include "astrometry_solver.h"

#define DOWNSAMPLE_FACTOR 0.25
#define CONV_TOLERANCE 1E-8

enum {
	COLUMN_RESOLVER,		// string
	COLUMN_NAME,		// string
	N_COLUMNS
};

#undef DEBUG           /* get some of diagnostic output */

typedef enum {
	RESOLVER_NED,
	RESOLVER_SIMBAD,
	RESOLVER_VIZIER,
	RESOLVER_NUMBER
} resolver;

struct object {
	gchar *name;
	double radius;
	int maxRecords;
	SirilWorldCS *world_cs;
	point imageCenter;
	gboolean south;
};

struct image_solved_struct {
	point size;
	SirilWorldCS *px_cat_center;
	SirilWorldCS *image_center;
	double crpix[2];
	double pixel_size, focal;
	Homography H;
};

void on_GtkTreeViewIPS_cursor_changed(GtkTreeView *tree_view,
		gpointer user_data);

static struct object platedObject[RESOLVER_NUMBER];
static GtkListStore *list_IPS = NULL;

static void initialize_ips_dialog() {
	GtkWidget *button_ips_ok, *button_cc_ok, *catalog_label, *catalog_box_ips,
			*catalog_box_pcc, *catalog_auto, *frame_cc_bkg, *frame_cc_norm,
			*catalog_label_pcc;
	GtkWindow *parent;

	button_ips_ok = lookup_widget("buttonIPS_ok");
	button_cc_ok = lookup_widget("button_cc_ok");
	catalog_label = lookup_widget("GtkLabelCatalog");
	catalog_label_pcc = lookup_widget("GtkLabelCatalogPCC");
	catalog_box_ips = lookup_widget("ComboBoxIPSCatalog");
	catalog_box_pcc = lookup_widget("ComboBoxPCCCatalog");
	catalog_auto = lookup_widget("GtkCheckButton_OnlineCat");
	frame_cc_bkg = lookup_widget("frame_cc_background");
	frame_cc_norm = lookup_widget("frame_cc_norm");

	parent = GTK_WINDOW(lookup_widget("ImagePlateSolver_Dial"));

	gtk_widget_set_visible(button_ips_ok, TRUE);
	gtk_widget_set_visible(button_cc_ok, FALSE);
	gtk_widget_set_visible(catalog_label, TRUE);
	gtk_widget_set_visible(catalog_label_pcc, FALSE);
	gtk_widget_set_visible(catalog_box_ips, TRUE);
	gtk_widget_set_visible(catalog_box_pcc, FALSE);
	gtk_widget_set_visible(catalog_auto, TRUE);
	gtk_widget_set_visible(frame_cc_bkg, FALSE);
	gtk_widget_set_visible(frame_cc_norm, FALSE);

	gtk_window_set_title(parent, _("Image Plate Solver"));
}

static void fov_in_DHMS(double var, gchar *fov) {
	int deg, decM;
	double decS;

	if (var < 0) {
		fprintf(stdout, "fov_in_DHMS: negative value, should not happened\n");
		return;
	}
	deg = (int) var;
	decM = abs((int) ((var - deg) * 60));
	decS = (fabs((var - deg) * 60) - decM) * 60;
	if (deg > 0)
		g_snprintf(fov, 256, "%02dd %02dm %.2lfs", deg, decM, decS);
	else if (decM > 0)
		g_snprintf(fov, 256, "%02d\' %.2lf\"", decM, decS);
	else if (decS > 0.0)
		g_snprintf(fov, 256, "%.2lf\"", decS);
}

static int parse_content_buffer(char *buffer, struct object *obj) {
	char **token, **fields;
	point center;
	int nargs, i = 0, resolver = -1;

	token = g_strsplit(buffer, "\n", -1);
	nargs = g_strv_length(token);

	while (i < nargs) {
		if (g_strrstr (token[i], "=NED")) {
			resolver = RESOLVER_NED;
		} else if (g_strrstr (token[i], "=Simbad")) {
			resolver = RESOLVER_SIMBAD;
		} else if (g_strrstr(token[i], "=VizieR")) {
			resolver = RESOLVER_VIZIER;
		} else if (g_str_has_prefix (token[i], "%J ")) {
			fields = g_strsplit(token[i], " ", -1);
			sscanf(fields[1], "%lf", &center.x);
			sscanf(fields[2], "%lf", &center.y);
			if (resolver != -1) {
				platedObject[resolver].world_cs = siril_world_cs_new_from_a_d(center.x, center.y);

				/* others */
				platedObject[resolver].imageCenter = center;
				platedObject[resolver].south = (center.y < 0.0);
			}
			g_strfreev(fields);
		} else if (g_str_has_prefix (token[i], "%I.0 ")) {
			if (resolver != -1) {
				gchar *name = g_strstr_len(token[i], strlen(token[i]), "%I.0 ");
				gchar *realname;
				realname = g_strdup(name + 5);
				platedObject[resolver].name = realname;
			}
		} else if (g_str_has_prefix (token[i], "%I NAME ")) {
			if (resolver != -1) {
				gchar *name = g_strstr_len(token[i], strlen(token[i]), "%I NAME ");
				gchar *realname;
				realname = g_strdup(name + 5 + 3);
				g_free(platedObject[resolver].name);
				platedObject[resolver].name = realname;
			}
		}
		i++;
	}

	g_strfreev(token);
	return 0;
}

static void free_Platedobject() {
	for (int i = 0; i < RESOLVER_NUMBER; i++) {
		if (platedObject[i].name) {
			siril_world_cs_unref(platedObject[i].world_cs);
			free(platedObject[i].name);
			platedObject[i].name = NULL;
		}
	}
}

static double get_focal() {
	GtkEntry *focal_entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_focal"));
	const gchar *value = gtk_entry_get_text(focal_entry);

	return g_ascii_strtod(value, NULL);
}

/* get pixel in µm */
static double get_pixel() {
	GtkEntry *pixel_entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_pixels"));
	const gchar *value = gtk_entry_get_text(pixel_entry);

	return g_ascii_strtod(value, NULL);
}

static double get_resolution(double focal, double pixel) {
	return RADCONV / focal * pixel;
}

/* get FOV in arcmin/px */
static double get_fov(double resolution, int image_size) {
	return (resolution * (double)image_size) / 60.0;
}

static double get_mag_limit(double fov) {
	GtkToggleButton *autobutton = GTK_TOGGLE_BUTTON(lookup_widget("GtkCheckButton_Mag_Limit"));
	if (gtk_toggle_button_get_active(autobutton)) {
		// Empiric formula for 1000 stars at 20 deg of galactic latitude
		double autoLimitMagnitudeFactor = 14.5;
		/* convert fov in degree */
		fov /= 60.0;
		double m = autoLimitMagnitudeFactor * pow(fov, -0.179);
		m = round(100 * min(20, max(7, m))) / 100;
		return m;
	} else {
		GtkSpinButton *magButton = GTK_SPIN_BUTTON(
				lookup_widget("GtkSpinIPS_Mag_Limit"));

		return gtk_spin_button_get_value(magButton);
	}
}

static SirilWorldCS *get_center_of_catalog() {
	GtkSpinButton *GtkSpinIPS_RA_h, *GtkSpinIPS_RA_m;
	GtkSpinButton *GtkSpinIPS_Dec_deg, *GtkSpinIPS_Dec_m;
	GtkEntry *GtkEntryIPS_RA_s, *GtkEntryIPS_Dec_s;
	GtkToggleButton *GtkCheckButtonIPS_S;

	/* get alpha center */
	GtkSpinIPS_RA_h = GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_RA_h"));
	GtkSpinIPS_RA_m = GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_RA_m"));
	GtkEntryIPS_RA_s = GTK_ENTRY(lookup_widget("GtkEntryIPS_RA_s"));

	gdouble hour = gtk_spin_button_get_value_as_int(GtkSpinIPS_RA_h);
	gdouble min = gtk_spin_button_get_value_as_int(GtkSpinIPS_RA_m);
	gdouble sec = g_ascii_strtod(gtk_entry_get_text(GtkEntryIPS_RA_s), NULL);

	/* get Dec center */
	GtkSpinIPS_Dec_deg = GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_Dec_deg"));
	GtkSpinIPS_Dec_m = GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_Dec_m"));
	GtkEntryIPS_Dec_s = GTK_ENTRY(lookup_widget("GtkEntryIPS_Dec_s"));
	GtkCheckButtonIPS_S = GTK_TOGGLE_BUTTON(lookup_widget("GtkCheckButtonIPS_S"));

	gdouble deg = gtk_spin_button_get_value_as_int(GtkSpinIPS_Dec_deg);
	gdouble m = gtk_spin_button_get_value_as_int(GtkSpinIPS_Dec_m);
	gdouble s = g_ascii_strtod(gtk_entry_get_text(GtkEntryIPS_Dec_s), NULL);
	if (gtk_toggle_button_get_active(GtkCheckButtonIPS_S)) {
		deg = -deg;
	}

	return siril_world_cs_new_from_ra_dec(hour, min, sec, deg, m, s);;
}

static gboolean is_detection_manual() {
	GtkToggleButton *button;

	button = GTK_TOGGLE_BUTTON(lookup_widget("checkButton_IPS_manual"));
	return gtk_toggle_button_get_active(button);
}

static gboolean flip_image_after_ps() {
	GtkToggleButton *button;

	button = GTK_TOGGLE_BUTTON(lookup_widget("checkButton_IPS_flip"));
	return gtk_toggle_button_get_active(button);
}

static gboolean is_downsample_activated() {
	GtkToggleButton *button;

	button = GTK_TOGGLE_BUTTON(lookup_widget("downsample_ips_button"));
	return gtk_toggle_button_get_active(button);
}

static gboolean is_cache_activated() {
	GtkToggleButton *button;

	button = GTK_TOGGLE_BUTTON(lookup_widget("use_cache_ips"));
	return gtk_toggle_button_get_active(button);
}

static gboolean is_autocrop_activated() {
	GtkToggleButton *button;

	button = GTK_TOGGLE_BUTTON(lookup_widget("autocrop_ips_button"));
	return gtk_toggle_button_get_active(button);
}

static gchar *get_catalog_url(SirilWorldCS *center, double mag_limit, double dfov, int type) {
	GString *url;
	gchar *coordinates;
	gchar *mag;
	gchar *fov;

	coordinates = g_strdup_printf("%f+%f", siril_world_cs_get_alpha(center), siril_world_cs_get_delta(center));
	mag = g_strdup_printf("%2.2lf", mag_limit);
	fov = g_strdup_printf("%2.1lf", dfov / 2.0);

	url = g_string_new("http://vizier.u-strasbg.fr/viz-bin/asu-tsv?-source=");
	switch (type) {
	case NOMAD:
		url = g_string_append(url, "NOMAD&-out.meta=-h-u-D&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-out=%20RAJ2000%20DEJ2000%20Vmag%20Bmag");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&Vmag=<");
		url = g_string_append(url, mag);
		break;
	default:
	case TYCHO2:
		url = g_string_append(url, "I/259/tyc2&-out.meta=-h-u-D&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-out=%20RAmdeg%20DEmdeg%20VTmag%20BTmag");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&VTmag=<");
		url = g_string_append(url, mag);
		break;
	case GAIA:
		url = g_string_append(url, "I/345/gaia2&-out.meta=-h-u-D&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-out=%20RAJ2000%20DEJ2000%20Gmag%20BPmag");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&Gmag=<");
		url = g_string_append(url, mag);
		break;
	case GAIAEDR3:
		url = g_string_append(url, "I/350/gaiaedr3&-out.meta=-h-u-D&-out.add=_r");
		url = g_string_append(url, "&-out=%20RAJ2000%20DEJ2000%20Gmag%20BPmag");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&Gmag=<");
		url = g_string_append(url, mag);
		break;
	case PPMXL:
		url = g_string_append(url, "I/317&-out.meta=-h-u-D&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-out=%20RAJ2000%20DEJ2000%20Jmag%20Hmag");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&Jmag=<");
		url = g_string_append(url, mag);
		break;
	case BRIGHT_STARS:
		url = g_string_append(url, "V/50/catalog&-out.meta=-h-u-D&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-out.add=_RAJ,_DEJ&-out=Vmag&-out=B-V");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&Vmag=<");
		url = g_string_append(url, mag);
		break;
	case APASS: // for photometry only
		url = g_string_append(url, "APASS&-out.meta=-h-u-D&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-out=%20RAJ2000%20DEJ2000%20Vmag%20Bmag");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&Vmag=<");
		url = g_string_append(url, mag);
		break;
	}

	g_free(coordinates);
	g_free(mag);
	g_free(fov);

	return g_string_free(url, FALSE);
}

#if defined HAVE_LIBCURL
/*****
 * HTTP functions
 ****/

static CURL *curl;
static const int DEFAULT_FETCH_RETRIES = 10;

struct ucontent {
	char *data;
	size_t len;
};

static void init() {
	if (!curl) {
		printf("initializing CURL\n");
		curl_global_init(CURL_GLOBAL_ALL);
		curl = curl_easy_init();
	}

	if (!curl)
		exit(EXIT_FAILURE);
}

static size_t cbk_curl(void *buffer, size_t size, size_t nmemb, void *userp) {
	size_t realsize = size * nmemb;
	struct ucontent *mem = (struct ucontent *) userp;

	mem->data = realloc(mem->data, mem->len + realsize + 1);

	memcpy(&(mem->data[mem->len]), buffer, realsize);
	mem->len += realsize;
	mem->data[mem->len] = 0;

	return realsize;
}

static char *fetch_url(const char *url) {
	struct ucontent *content = malloc(sizeof(struct ucontent));
	char *result, *error;
	long code;
	int retries;
	unsigned int s;

	siril_debug_print("fetch_url(): %s\n", url);

	init();

	result = NULL;

	retries = DEFAULT_FETCH_RETRIES;

	retrieve: content->data = malloc(1);
	content->data[0] = '\0';
	content->len = 0;

	curl_easy_setopt(curl, CURLOPT_URL, url);
	curl_easy_setopt(curl, CURLOPT_VERBOSE, 0);
	curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, cbk_curl);
	curl_easy_setopt(curl, CURLOPT_WRITEDATA, content);
	curl_easy_setopt(curl, CURLOPT_USERAGENT, PACKAGE_STRING);

	if (curl_easy_perform(curl) == CURLE_OK) {
		curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &code);

		switch (code) {
		case 200:
			result = content->data;
			break;
		case 500:
		case 502:
		case 503:
		case 504:
			printf("Fetch failed with code %ld for URL %s\n", code, url);

			if (retries) {
				s = 2 * (DEFAULT_FETCH_RETRIES - retries) + 2;
				printf("Wait %uds before retry\n", s);
				sleep(s);

				free(content->data);
				retries--;
				goto retrieve;
			}

			break;
		default:
			error = siril_log_message(_("Fetch failed with code %ld for URL %s\n"), code, url);
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), error);
		}
	}

	if (!result)
		free(content->data);

	free(content);

	return result;
}
#else
static gchar *fetch_url(const gchar *url) {
	GFile *file = g_file_new_for_uri(url);
	GError *error = NULL;
	gchar *content = NULL;

	siril_debug_print("fetch_url(): %s\n", url);

	if (!g_file_load_contents(file, NULL, &content, NULL, NULL, &error)) {
		siril_log_message(_("Error loading url: %s: %s\n"), url, error->message);
		g_clear_error(&error);
	}
	g_object_unref(file);
	return content;
}
#endif

static online_catalog get_online_catalog(double fov, double m) {
	GtkToggleButton *auto_button;
	int ret;

	auto_button = GTK_TOGGLE_BUTTON(lookup_widget("GtkCheckButton_OnlineCat"));
	if (gtk_toggle_button_get_active(auto_button)) {
		if (m <= 6.5) {
			ret = BRIGHT_STARS;
		} else if (fov > 180.0) {
			ret = NOMAD;
		} else if (fov < 30.0){
			ret = GAIAEDR3;
		} else {
			ret = PPMXL;
		}
		return ret;
	} else {
		GtkComboBox *box = GTK_COMBO_BOX(lookup_widget("ComboBoxIPSCatalog"));
		ret = gtk_combo_box_get_active(box);
		return (ret < 0 ? NOMAD : ret);
	}
}

static GFile *download_catalog(gboolean use_cache, online_catalog onlineCatalog, SirilWorldCS *catalog_center, double fov, double m) {
	gchar *buffer = NULL;
	GError *error = NULL;

	/* ------------------- get Vizier catalog in catalog.dat -------------------------- */

	/* check if catalogue already exist in cache */
	gchar *str = g_strdup_printf("cat-%d-%lf-%lf-%lf-%lf.cat",
			(int) onlineCatalog,
			siril_world_cs_get_alpha(catalog_center),
			siril_world_cs_get_delta(catalog_center),
			fov, m);
	siril_debug_print("Catalogue file: %s\n", str);
	GFile *file = g_file_new_build_filename(g_get_tmp_dir(), str, NULL);
	g_free(str);

	GOutputStream *output_stream = (GOutputStream*) g_file_create(file, G_FILE_CREATE_NONE, NULL, &error);

	if (output_stream == NULL) {
		if (error != NULL) {
			/* if file exist and user uses cache */
			if (error->code == G_IO_ERROR_EXISTS && use_cache) {
				siril_log_color_message(_("Using data in cache\n"), "salmon");
				g_clear_error(&error);
			} else {
				g_warning("%s\n", error->message);
				g_clear_error(&error);
				fprintf(stderr, "astrometry_solver: Cannot open catalogue\n");
				g_object_unref(file);
				return NULL;
			}
		}
	} else {
		gchar *url = get_catalog_url(catalog_center, m, fov, onlineCatalog);
		buffer = fetch_url(url);
		g_free(url);

		if (buffer) {
			if (!g_output_stream_write_all(output_stream, buffer, strlen(buffer), NULL, NULL, &error)) {
				g_warning("%s\n", error->message);
				g_clear_error(&error);
				g_free(buffer);
				g_object_unref(output_stream);
				g_object_unref(file);
				return NULL;
			}
			g_object_unref(output_stream);
			g_free(buffer);
		}
	}

	return file;
}


static gchar *project_catalog(GFile *catalogue_name, SirilWorldCS *catalog_center) {
	GError *error = NULL;
	gchar *foutput = NULL;
	/* -------------------------------------------------------------------------------- */

	/* --------- Project coords of Vizier catalog and save it into catalog.proj ------- */

	GFile *fproj = g_file_new_build_filename(g_get_tmp_dir(), "catalog.proj", NULL);

	/* We want to remove the file if already exisit */
	if (!g_file_delete(fproj, NULL, &error)
			&& !g_error_matches(error, G_IO_ERROR, G_IO_ERROR_NOT_FOUND)) {
		// deletion failed for some reason other than the file not existing:
		// so report the error
		g_warning("Failed to delete %s: %s", g_file_peek_path(fproj),
				error->message);
	}

	convert_catalog_coords(catalogue_name, catalog_center, fproj);
	foutput = g_file_get_path(fproj);
	g_object_unref(fproj);

	/* -------------------------------------------------------------------------------- */

//	solution.px_cat_center = siril_world_cs_ref(catalog_center);
	return foutput;
}

/*********
 *
 */

static void get_list_IPS() {
	if (list_IPS == NULL)
		list_IPS = GTK_LIST_STORE(gtk_builder_get_object(builder, "liststoreIPS"));
}

static void clear_all_objects() {
	gtk_list_store_clear(list_IPS);
}

static void add_object_to_list() {
	GtkTreeIter iter;

	get_list_IPS();
	clear_all_objects();

	if (platedObject[RESOLVER_NED].name) {
		gtk_list_store_append(list_IPS, &iter);
		gtk_list_store_set(list_IPS, &iter, COLUMN_RESOLVER, "NED", COLUMN_NAME,
				platedObject[RESOLVER_NED].name, -1);
	}

	if (platedObject[RESOLVER_SIMBAD].name) {
		gtk_list_store_append(list_IPS, &iter);
		gtk_list_store_set(list_IPS, &iter, COLUMN_RESOLVER, "Simbad",
				COLUMN_NAME, platedObject[RESOLVER_SIMBAD].name, -1);
	}

	if (platedObject[RESOLVER_VIZIER].name) {
		gtk_list_store_append(list_IPS, &iter);
		gtk_list_store_set(list_IPS, &iter, COLUMN_RESOLVER, "VizieR",
				COLUMN_NAME, platedObject[RESOLVER_VIZIER].name, -1);
	}
}

static void unselect_all_items() {
	GtkTreeSelection *selection;

	selection = GTK_TREE_SELECTION(gtk_builder_get_object(builder, "gtkselectionIPS"));
	gtk_tree_selection_unselect_all(selection);
}

static void update_coordinates(SirilWorldCS *world_cs) {
	gchar *RA_sec, *Dec_sec;
	gint ra_h, ra_m;
	gint dec_deg, dec_m;
	gdouble ra_s, dec_s;

	siril_world_cs_get_ra_hour_min_sec(world_cs, &ra_h, &ra_m, &ra_s);
	siril_world_cs_get_dec_deg_min_sec(world_cs, &dec_deg, &dec_m, &dec_s);

	RA_sec = g_strdup_printf("%6.4lf", ra_s);
	Dec_sec = g_strdup_printf("%6.4lf", dec_s);

	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("GtkCheckButtonIPS_S")), dec_deg < 0);

	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_RA_h")), ra_h);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_RA_m")), ra_m);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("GtkEntryIPS_RA_s")), RA_sec);

	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_Dec_deg")), abs(dec_deg));
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_Dec_m")), dec_m);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("GtkEntryIPS_Dec_s")), Dec_sec);

	g_free(RA_sec);
	g_free(Dec_sec);
}

static gboolean has_any_keywords() {
	return (gfit.focal_length > 0.0 ||
			gfit.pixel_size_x > 0.f ||
			gfit.pixel_size_y > 0.f ||
			(gfit.wcsdata.crval[0] > 0.0 && gfit.wcsdata.crval[1] != 0.0) ||
			(gfit.wcsdata.objctra[0] != '\0' && gfit.wcsdata.objctdec[0] != '\0') ||
			(gfit.wcsdata.ra != 0.0 && gfit.wcsdata.dec != 0.0));
}

static void update_coords() {
	SirilWorldCS *world_cs = NULL;
	fits *fit = &gfit;

	if (fit->wcsdata.ra != 0.0 && fit->wcsdata.dec != 0.0) {

		world_cs = siril_world_cs_new_from_a_d(fit->wcsdata.ra, fit->wcsdata.dec);

		update_coordinates(world_cs);
		unselect_all_items();
	} else if (fit->wcsdata.objctra[0] != '\0' && fit->wcsdata.objctdec[0] != '\0') {

		world_cs = siril_world_cs_new_from_objct_ra_dec(fit->wcsdata.objctra, fit->wcsdata.objctdec);

		update_coordinates(world_cs);
		unselect_all_items();
	} else if (fit->wcsdata.crval[0] != 0.0 && fit->wcsdata.crval[1] != 0.0) {

		// first transform coords to alpha and delta
		world_cs = siril_world_cs_new_from_a_d(fit->wcsdata.crval[0], fit->wcsdata.crval[1]);

		update_coordinates(world_cs);
		unselect_all_items();
	}
	if (world_cs)
		siril_world_cs_unref(world_cs);
}

static void update_pixel_size() {
	GtkEntry *entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_pixels"));
	float pixel;

	pixel = gfit.pixel_size_x > gfit.pixel_size_y ? gfit.pixel_size_x : gfit.pixel_size_y;

	if (pixel > 0.f) {
		gchar *cpixels = g_strdup_printf("%.2lf", (double) pixel);
		gtk_entry_set_text(entry, cpixels);
		g_free(cpixels);
	}
}

static void update_focal() {
	GtkEntry *entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_focal"));
	double focal;

	focal = gfit.focal_length;

	if (focal > 0.0) {
		gchar *cfocal = g_strdup_printf("%.1lf", focal);
		gtk_entry_set_text(entry, cfocal);
		g_free(cfocal);
	}
}

static void update_resolution_field() {
	GtkEntry *entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_resolution"));
	double res = get_resolution(get_focal(), get_pixel());
	gchar *cres;

	cres = g_strdup_printf("%1.3lf", res);
	gtk_entry_set_text(entry, cres);
	g_free(cres);
}

static void update_image_parameters_GUI() {
	/* update all fields. Resolution is updated as well
	 thanks to the Entry and combo changed signal
	  */
	update_focal();
	update_pixel_size();
	update_coords();
}

/* Extract CDELT from CD matrix.*/
static void extract_cdelt_from_cd(double cd1_1, double cd1_2, double cd2_1,
		double cd2_2, double *cdelt1, double *cdelt2) {
	int sign;
	if ((cd1_1 * cd2_2 - cd1_2 * cd2_1) >= 0)
		sign = +1;
	else
		sign = -1;

	*cdelt1 = sqrt((cd1_1 * cd1_1) + (cd2_1 * cd2_1)) * sign;
	*cdelt2 = sqrt((cd1_2 * cd1_2) + (cd2_2 * cd2_2));
}

static void print_platesolving_results(image_solved *image, gboolean downsample) {
	double rotation, det, scaleX, scaleY, resolution;
	double inliers;
	gchar *alpha, *delta;
	char field_x[256] = { 0 };
	char field_y[256] = { 0 };
	point fov;

	float factor = (downsample) ? DOWNSAMPLE_FACTOR : 1.0;
	Homography H = image->H;

	/* Matching information */
	gchar *str = ngettext("%d pair match.\n", "%d pair matches.\n", H.pair_matched);
	str = g_strdup_printf(str, H.pair_matched);
	siril_log_message(str);
	g_free(str);
	inliers = 1.0 - ((((double) H.pair_matched - (double) H.Inliers)) / (double) H.pair_matched);
	siril_log_message(_("Inliers:%*.3f\n"), 14, inliers);

	/* Plate Solving */
	scaleX = sqrt(H.h00 * H.h00 + H.h01 * H.h01);
	scaleY = sqrt(H.h10 * H.h10 + H.h11 * H.h11);
	resolution = (scaleX + scaleY) * 0.5 * factor; // we assume square pixels
	siril_log_message(_("Resolution:%*.3lf arcsec/px\n"), 11, resolution);

	/* rotation */
	rotation = atan2(H.h00 + H.h01, H.h10 + H.h11) * RADTODEG + 135.0;
	det = (H.h00 * H.h11 - H.h01 * H.h10); // determinant of rotation matrix (ad - bc)
	/* If the determinant of the top-left 2x2 rotation matrix is > 0
	 * the transformation is orientation-preserving. */

	if (det < 0)
		rotation = -90.0 - rotation;
	if (rotation < -180.0)
		rotation += 360.0;
	if (rotation > 180.0)
		rotation -= 360.0;
	siril_log_message(_("Rotation:%+*.2lf deg %s\n"), 12, rotation, det < 0.0 ? _("(flipped)") : "");

	fov.x = get_fov(resolution, image->size.x);
	fov.y = get_fov(resolution, image->size.y);
	siril_log_message(_("Focal:%*.2lf mm\n"), 15, RADCONV * image->pixel_size / resolution);
	siril_log_message(_("Pixel size:%*.2lf µm\n"), 10, image->pixel_size);
	fov_in_DHMS(fov.x / 60.0, field_x);
	fov_in_DHMS(fov.y / 60.0, field_y);
	siril_log_message(_("Field of view:    %s x %s\n"), field_x, field_y);

	alpha = siril_world_cs_alpha_format(image->image_center, " %02dh%02dm%02ds");
	delta = siril_world_cs_delta_format(image->image_center, "%c%02d°%02d\'%02d\"");
	siril_log_message(_("Image center: alpha: %s, delta: %s\n"), alpha, delta);

	g_free(alpha);
	g_free(delta);
}

static int read_NOMAD_catalog(GInputStream *stream, psf_star **cstars) {
	gchar *line;
	psf_star *star;

	int i = 0;

	GDataInputStream *data_input = g_data_input_stream_new(stream);
	while (i < MAX_STARS &&
			(line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {
		double r = 0.0, x = 0.0, y = 0.0, Vmag = 0.0, Bmag = 0.0;

		if (line[0] == COMMENT_CHAR) {
			g_free(line);
			continue;
		}
		if (is_blank(line)) {
			g_free(line);
			continue;
		}
		if (g_str_has_prefix(line, "---")) {
			g_free(line);
			continue;
		}
		int n = sscanf(line, "%lf %lf %lf %lf %lf", &r, &x, &y, &Vmag, &Bmag);

		star = new_psf_star();
		star->xpos = x;
		star->ypos = y;
		star->mag = Vmag;
		star->BV = n < 5 ? -99.9 : Bmag - Vmag;
		star->phot = NULL;
		cstars[i] = star;
		cstars[i + 1] = NULL;
		i++;
		g_free(line);
	}
	g_object_unref(data_input);
	sort_stars(cstars, i);
	siril_log_message(_("Catalog NOMAD size: %d objects\n"), i);
	return i;
}

static int read_TYCHO2_catalog(GInputStream *stream, psf_star **cstars) {
	gchar *line;
	psf_star *star;

	int i = 0;

	GDataInputStream *data_input = g_data_input_stream_new(stream);
	while (i < MAX_STARS &&
			(line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {
		double r = 0.0, x = 0.0, y = 0.0, Vmag = 0.0, Bmag = 0.0;

		if (line[0] == COMMENT_CHAR) {
			continue;
		}
		if (is_blank(line)) {
			continue;
		}
		if (g_str_has_prefix(line, "---")) {
			continue;
		}
		int n = sscanf(line, "%lf %lf %lf %lf %lf", &r, &x, &y, &Vmag, &Bmag);

		star = new_psf_star();
		star->xpos = x;
		star->ypos = y;
		star->mag = Vmag;
		star->BV = n < 5 ? -99.9 : Bmag - Vmag;
		star->phot = NULL;
		cstars[i] = star;
		cstars[i + 1] = NULL;
		i++;
	}
	g_object_unref(data_input);
	sort_stars(cstars, i);
	siril_log_message(_("Catalog TYCHO-2 size: %d objects\n"), i);
	return i;
}

static int read_GAIA_catalog(GInputStream *stream, psf_star **cstars) {
	gchar *line;
	psf_star *star;

	int i = 0;

	GDataInputStream *data_input = g_data_input_stream_new(stream);
	while (i < MAX_STARS &&
			(line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {
		double r = 0.0, x = 0.0, y = 0.0, Gmag = 0.0, BPmag = 0.0;

		if (line[0] == COMMENT_CHAR) {
			g_free(line);
			continue;
		}
		if (is_blank(line)) {
			g_free(line);
			continue;
		}
		if (g_str_has_prefix(line, "---")) {
			g_free(line);
			continue;
		}
		sscanf(line, "%lf %lf %lf %lf %lf", &r, &x, &y, &Gmag, &BPmag);

		star = new_psf_star();
		star->xpos = x;
		star->ypos = y;
		star->mag = Gmag;
		star->BV = -99.9;
		star->phot = NULL;
		cstars[i] = star;
		cstars[i + 1] = NULL;
		i++;
		g_free(line);
	}
	g_object_unref(data_input);
	sort_stars(cstars, i);
	siril_log_message(_("Catalog Gaia DR2 size: %d objects\n"), i);
	return i;
}

static int read_PPMXL_catalog(GInputStream *stream, psf_star **cstars) {
	gchar *line;
	psf_star *star;

	int i = 0;

	GDataInputStream *data_input = g_data_input_stream_new(stream);
	while (i < MAX_STARS &&
			(line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {
		double r = 0.0, x = 0.0, y = 0.0, Jmag = 0.0, Hmag = 0.0;

		if (line[0] == COMMENT_CHAR) {
			g_free(line);
			continue;
		}
		if (is_blank(line)) {
			g_free(line);
			continue;
		}
		if (g_str_has_prefix(line, "---")) {
			g_free(line);
			continue;
		}
		sscanf(line, "%lf %lf %lf %lf %lf", &r, &x, &y, &Jmag, &Hmag);

		star = new_psf_star();
		star->xpos = x;
		star->ypos = y;
		star->mag = Jmag;
		star->BV = -99.9;
		star->phot = NULL;
		cstars[i] = star;
		cstars[i + 1] = NULL;
		i++;
		g_free(line);
	}
	g_object_unref(data_input);
	sort_stars(cstars, i);
	siril_log_message(_("Catalog PPMXL size: %d objects\n"), i);
	return i;
}

static int read_BRIGHT_STARS_catalog(GInputStream *stream, psf_star **cstars) {
	gchar *line;
	psf_star *star;

	int i = 0;

	GDataInputStream *data_input = g_data_input_stream_new(stream);
	while (i < MAX_STARS &&
			(line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {
		double r = 0.0, x = 0.0, y = 0.0, Vmag = 0.0, BV = 0.0;

		if (line[0] == COMMENT_CHAR) {
			g_free(line);
			continue;
		}
		if (is_blank(line)) {
			g_free(line);
			continue;
		}
		if (g_str_has_prefix(line, "---")) {
			g_free(line);
			continue;
		}
		sscanf(line, "%lf %lf %lf %lf %lf", &r, &x, &y, &Vmag, &BV);

		star = new_psf_star();
		star->xpos = x;
		star->ypos = y;
		star->mag = Vmag;
		star->BV = BV;
		star->phot = NULL;
		cstars[i] = star;
		cstars[i + 1] = NULL;
		i++;
		g_free(line);
	}
	g_object_unref(data_input);
	sort_stars(cstars, i);
	siril_log_message(_("Catalog Bright stars size: %d objects\n"), i);
	return i;
}

static int read_APASS_catalog(GInputStream *stream, psf_star **cstars) {
	gchar *line;
	psf_star *star;

	int i = 0;

	GDataInputStream *data_input = g_data_input_stream_new(stream);
	while (i < MAX_STARS &&
			(line = g_data_input_stream_read_line_utf8(data_input, NULL, NULL, NULL))) {
		double r = 0.0, x = 0.0, y = 0.0, Vmag = 0.0, Bmag = 0.0;

		if (line[0] == COMMENT_CHAR) {
			g_free(line);
			continue;
		}
		if (is_blank(line)) {
			g_free(line);
			continue;
		}
		if (g_str_has_prefix(line, "---")) {
			g_free(line);
			continue;
		}
		int n = sscanf(line, "%lf %lf %lf %lf %lf", &r, &x, &y, &Vmag, &Bmag);

		star = new_psf_star();
		star->xpos = x;
		star->ypos = y;
		star->mag = Vmag;
		star->BV = n < 5 ? -99.9 : Bmag - Vmag;
		star->phot = NULL;
		cstars[i] = star;
		cstars[i + 1] = NULL;
		i++;
		g_free(line);
	}
	g_object_unref(data_input);
	sort_stars(cstars, i);
	siril_log_message(_("Catalog APASS size: %d objects\n"), i);
	return i;
}

static int read_catalog(GInputStream *stream, psf_star **cstars, int type) {
	switch (type) {
	default:
	case TYCHO2:
		return read_TYCHO2_catalog(stream, cstars);
	case NOMAD:
		return read_NOMAD_catalog(stream, cstars);
	case GAIA:
		return read_GAIA_catalog(stream, cstars);
	case GAIAEDR3:
		return read_GAIA_catalog(stream, cstars);
	case PPMXL:
		return read_PPMXL_catalog(stream, cstars);
	case BRIGHT_STARS:
		return read_BRIGHT_STARS_catalog(stream, cstars);
	case APASS:
		return read_APASS_catalog(stream, cstars);
	}
}

static TRANS H_to_linear_TRANS(Homography H) {
	TRANS trans;

	trans.order = AT_TRANS_LINEAR;

	trans.a = H.h02;
	trans.b = H.h00;
	trans.c = H.h01;
	trans.d = H.h12;
	trans.e = H.h10;
	trans.f = H.h11;

	return trans;
}

static gboolean check_affine_TRANS_sanity(TRANS trans) {
	double var1 = fabs(trans.b) - fabs(trans.f);
	double var2 = fabs(trans.c) - fabs(trans.e);
	siril_debug_print("abs(b+f)=%f et abs(c+e)=%f\n", var1, var2);

	return ((fabs(var1) < 0.3) && fabs(var2) < 0.3);
}

/* Determine if the image need to be flipped
 * 1 - is the image flipped? Need to compute the determinant of H
 * 2 - did the user asked for a flip?
 */
static gboolean flip_image(gboolean flip_image, Homography H) {
	double det = (H.h00 * H.h11 - H.h01 * H.h10); // determinant of rotation matrix (ad - bc)
	return flip_image && det < 0;
}

static gboolean end_plate_solver(gpointer p) {
	struct astrometry_data *args = (struct astrometry_data *) p;
	stop_processing_thread();

	if (args->downsample) {
		/* copying size of original fits */
		args->fit->naxes[0] = args->fit_backup->naxes[0];
		args->fit->naxes[1] = args->fit_backup->naxes[1];
		args->fit->rx = args->fit->naxes[0];
		args->fit->ry = args->fit->naxes[1];

		copyfits(args->fit_backup, args->fit, CP_ALLOC | CP_COPYA, -1);

		clearfits(args->fit_backup);
	}

	if (!args->manual)
		clear_stars_list();
	set_cursor_waiting(FALSE);

	if (args->ret) {
		char *title = siril_log_color_message(_("Plate Solving failed. "
				"The image could not be aligned with the reference stars.\n"), "red");
		if (!args->message) {
			args->message = g_strdup(_("This is usually because the initial parameters (pixel size, focal length, initial coordinates) "
					"are too far from the real metadata of the image.\n\n"
					"You could also try to look into another catalogue, or try to click on the \"Downsampling\" button, especially for image done with Drizzle.\n\n"
					"Finally, keep in mind that plate solving algorithm should only be applied on linear image."));
		}
		siril_message_dialog(GTK_MESSAGE_ERROR, title, args->message);
	} else {
		image_solved *solution = args->solution;
		print_platesolving_results(solution, args->downsample);

		/* update UI */
		update_image_parameters_GUI();
		set_GUI_CAMERA();
		update_coordinates(solution->image_center);
		siril_world_cs_unref(solution->px_cat_center);
		siril_world_cs_unref(solution->image_center);
		delete_selected_area();
		/* ****** */

		control_window_switch_to_tab(OUTPUT_LOGS);
		if (args->for_photometry_cc) {
			apply_photometric_cc();
		}
		if (flip_image(args->flip_image, solution->H)) {
			siril_log_color_message(_("Flipping image and updating astrometry data.\n"), "salmon");
			fits_flip_top_to_bottom(args->fit);
			flip_bottom_up_astrometry_data(args->fit);
			redraw(com.cvport, REMAP_ALL);
		} else {
			redraw(com.cvport, REMAP_NONE);
		}
		load_WCS_from_memory(args->fit);
		free(solution);
	}
	update_MenuItem();
	g_object_unref(args->catalog_name);
	g_free(args->catalogStars);
	g_free(args->message);
	free(args);
	
	return FALSE;
}

static void add_object_in_tree_view(const gchar *object) {
	struct object obj;
	GtkTreeView *GtkTreeViewIPS;

	GtkTreeViewIPS = GTK_TREE_VIEW(lookup_widget("GtkTreeViewIPS"));

	set_cursor_waiting(TRUE);

	gchar *result = search_in_catalogs(object);
	if (result) {
		free_Platedobject();
		parse_content_buffer(result, &obj);
		g_signal_handlers_block_by_func(GtkTreeViewIPS, on_GtkTreeViewIPS_cursor_changed, NULL);
		add_object_to_list();
		g_signal_handlers_unblock_by_func(GtkTreeViewIPS, on_GtkTreeViewIPS_cursor_changed, NULL);
		g_free(result);
	}
	set_cursor_waiting(FALSE);
}

static void start_image_plate_solve() {
	struct astrometry_data *args = malloc(sizeof(struct astrometry_data));

	args->for_photometry_cc = FALSE;
	if (!fill_plate_solver_structure(args)) {
		set_cursor_waiting(TRUE);
		start_in_new_thread(match_catalog, args);
	}
}

// From projected starlist and center (ra,dec), go back to original ra and dec
// All formulas from AIPS memo #27 III.A.ii
// https://library.nrao.edu/public/memos/aips/memos/AIPSM_027.pdf

static void deproject_starlist(int num_stars, s_star *star_list, double ra0, double dec0, int doASEC) {
	ra0 *= DEGTORAD;
	dec0 *= DEGTORAD;
	s_star *currstar;
	currstar = star_list;
	for (int i = 0; i < num_stars; i++) {
		double xi = currstar->x;
		double eta = currstar->y;
		if (doASEC > 0) {
			xi /= RADtoASEC;
			eta /= RADtoASEC;
		}
		double delta_ra = atan(xi / (cos(dec0) - eta * sin(dec0)));
		double ra = ra0 + delta_ra;
		double dec = atan(cos(delta_ra) * (eta * cos(dec0) + sin(dec0)) / (cos(dec0) - eta * sin(dec0)));
		currstar->x = ra / DEGTORAD;
		currstar->y = dec / DEGTORAD;
		currstar = currstar->next;
	}
}

// From starlist in (ra,dec) and center (ra,dec), project in "pixels" (in arcsec)
// All formulas from AIPS memo #27 III.A.i
// https://library.nrao.edu/public/memos/aips/memos/AIPSM_027.pdf

static void project_starlist(int num_stars, s_star *star_list, double ra0, double dec0, int doASEC) {
	double delta_ra;
	dec0 *= DEGTORAD;
	s_star *currstar;
	currstar = star_list;
	for (int i = 0; i < num_stars; i++) {
		double ra = currstar->x;
		double dec = currstar->y;
		if ((ra < 10) && (ra0 > 350)) {
			delta_ra = (ra + 360) - ra0;
		} else if ((ra > 350) && (ra0 < 10)) {
			delta_ra = (ra - 360) - ra0;
		} else {
			delta_ra = ra - ra0;
		}
		delta_ra *= DEGTORAD;
		dec *= DEGTORAD;

		/*
		 * let's transform from (delta_RA, delta_Dec) to (xi, eta),
		 */
		double xx = cos(dec) * sin(delta_ra);
		double yy = sin(dec0) * sin(dec) + cos(dec0) * cos(dec) * cos(delta_ra);
		double xi = (xx / yy);
		xx = cos(dec0) * sin(dec) - sin(dec0) * cos(dec) * cos(delta_ra);
		double eta = (xx / yy);

		if (doASEC > 0) {
			xi *= RADtoASEC;
			eta *= RADtoASEC;
		}
		currstar->x = xi;
		currstar->y = eta;
		currstar = currstar->next;
	}
}


/*****
 * CALLBACKS FUNCTIONS
 */

void on_GtkEntry_IPS_focal_changed(GtkEditable *editable, gpointer user_data) {
	update_resolution_field();
	com.pref.focal = g_ascii_strtod(gtk_editable_get_chars(editable, 0, -1), NULL);
}

void on_GtkEntry_IPS_pixels_changed(GtkEditable *editable, gpointer user_data) {
	update_resolution_field();
	com.pref.pitch = g_ascii_strtod(gtk_editable_get_chars(editable, 0, -1), NULL);
}

void on_GtkEntry_IPS_insert_text(GtkEntry *entry, const gchar *text, gint length,
		gint *position, gpointer data) {
	GtkEditable *editable = GTK_EDITABLE(entry);
	int i, count = 0;

	gchar *result = g_strndup(text, length);

	for (i = 0; i < length; i++) {
		if (!g_ascii_isdigit(text[i]) && text[i] != '.')
			continue;
		result[count++] = text[i];
	}

	if (count > 0) {
		g_signal_handlers_block_by_func(G_OBJECT (editable),
				G_CALLBACK (on_GtkEntry_IPS_insert_text), data);
		gtk_editable_insert_text(editable, result, count, position);
		g_signal_handlers_unblock_by_func(G_OBJECT (editable),
				G_CALLBACK (on_GtkEntry_IPS_insert_text), data);
	}
	g_signal_stop_emission_by_name(G_OBJECT(editable), "insert_text");

	g_free(result);
}

void on_info_menu_astrometry_clicked(GtkButton *button, gpointer user_data) {
	open_astrometry_dialog();
}

void on_buttonIPS_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("ImagePlateSolver_Dial");
}

void on_GtkTreeViewIPS_cursor_changed(GtkTreeView *tree_view,
		gpointer user_data) {

	GtkTreeModel *treeModel = gtk_tree_view_get_model(tree_view);
	GtkTreeSelection *selection = gtk_tree_view_get_selection (tree_view);
	GtkTreeIter iter;
	GValue value = G_VALUE_INIT;

	if (gtk_tree_model_get_iter_first(treeModel, &iter) == FALSE)
		return;	//The tree is empty
	if (gtk_tree_selection_get_selected(selection, &treeModel, &iter)) { //get selected item
		gtk_tree_model_get_value(treeModel, &iter, COLUMN_RESOLVER, &value);
		const gchar *res = g_value_get_string(&value);
		int selected_item;

		if (!g_strcmp0(res, "NED")) {
			selected_item = 0;
		} else if (!g_strcmp0(res, "Simbad")) {
			selected_item = 1;
		} else if (!g_strcmp0(res, "VizieR")) {
			selected_item = 2;
		} else {
			selected_item = -1;
		}

		if (selected_item >= 0) {
			update_coordinates(platedObject[selected_item].world_cs);
		}

		g_value_unset(&value);
	}
}

void on_GtkButton_IPS_metadata_clicked(GtkButton *button, gpointer user_data) {
	if (!has_any_keywords()) {
		char *msg = siril_log_message(_("There are no keywords stored in the FITS header.\n"));
		siril_message_dialog(GTK_MESSAGE_WARNING, _("No metadata"), msg);
	} else {
		update_image_parameters_GUI();
	}
}

void on_GtkButtonIPS_clicked(GtkButton *button, gpointer user_data) {
	GtkEntry *entry;

	entry = GTK_ENTRY(lookup_widget("GtkSearchIPS"));
	add_object_in_tree_view(gtk_entry_get_text(GTK_ENTRY(entry)));
}

void on_buttonIPS_ok_clicked(GtkButton *button, gpointer user_data) {
	start_image_plate_solve();
}

void on_GtkSearchIPS_activate(GtkEntry *entry, gpointer user_data) {
	add_object_in_tree_view(gtk_entry_get_text(GTK_ENTRY(entry)));
}

void on_GtkCheckButton_Mag_Limit_toggled(GtkToggleButton *button,
		gpointer user_data) {
	GtkWidget *spinmag;

	spinmag = lookup_widget("GtkSpinIPS_Mag_Limit");
	gtk_widget_set_sensitive(spinmag, !gtk_toggle_button_get_active(button));
}

void on_GtkCheckButton_OnlineCat_toggled(GtkToggleButton *button,
		gpointer user_data) {
	GtkWidget *combobox;

	combobox = lookup_widget("ComboBoxIPSCatalog");
	gtk_widget_set_sensitive(combobox, !gtk_toggle_button_get_active(button));
}

void print_updated_wcs_data(fits *fit) {
	/* debug output */
	siril_debug_print("****Updated WCS data*************\n");
	siril_debug_print("crpix1 = %*.12e\n", 20, fit->wcsdata.crpix[0]);
	siril_debug_print("crpix2 = %*.12e\n", 20, fit->wcsdata.crpix[1]);
	siril_debug_print("crval1 = %*.12e\n", 20, fit->wcsdata.crval[0]);
	siril_debug_print("crval2 = %*.12e\n", 20, fit->wcsdata.crval[1]);
	siril_debug_print("cdelt1 = %*.12e\n", 20, fit->wcsdata.cdelt[0]);
	siril_debug_print("cdelt2 = %*.12e\n", 20, fit->wcsdata.cdelt[1]);
	siril_debug_print("pc1_1  = %*.12e\n", 20, fit->wcsdata.pc[0][0]);
	siril_debug_print("pc1_2  = %*.12e\n", 20, fit->wcsdata.pc[0][1]);
	siril_debug_print("pc2_1  = %*.12e\n", 20, fit->wcsdata.pc[1][0]);
	siril_debug_print("pc2_2  = %*.12e\n", 20, fit->wcsdata.pc[1][1]);
	siril_debug_print("******************************************\n");
}

/******
 *
 * Public functions
 */

void flip_bottom_up_astrometry_data(fits *fit) {
	/* flip pc matrix */
	fit->wcsdata.pc[0][1] = -fit->wcsdata.pc[0][1];
	fit->wcsdata.pc[1][1] = -fit->wcsdata.pc[1][1];

	/* update crpix */
	fit->wcsdata.crpix[1] = fit->ry - fit->wcsdata.crpix[1];

	print_updated_wcs_data(fit);
}

void flip_left_right_astrometry_data(fits *fit) {
	/* flip pc matrix */
	fit->wcsdata.pc[0][0] = -fit->wcsdata.pc[0][0];
	fit->wcsdata.pc[1][0] = -fit->wcsdata.pc[1][0];

	/* update crpix */
	fit->wcsdata.crpix[0] = fit->rx - fit->wcsdata.crpix[0];

	print_updated_wcs_data(fit);
}

void rotate_astrometry_data(fits *fit, point center, double angle, gboolean cropped) {
	double pc1_1, pc1_2, pc2_1, pc2_2;
	point refpointout;

	const double2 sincosval = xsincos(angle * DEGTORAD);
	double sa, ca;
	sa = sincosval.x;
	ca = sincosval.y;

	pc1_1 =  ca * fit->wcsdata.pc[0][0] + sa * fit->wcsdata.pc[0][1];
	pc1_2 = -sa * fit->wcsdata.pc[0][0] + ca * fit->wcsdata.pc[0][1];
	pc2_1 =  ca * fit->wcsdata.pc[1][0] + sa * fit->wcsdata.pc[1][1];
	pc2_2 = -sa * fit->wcsdata.pc[1][0] + ca * fit->wcsdata.pc[1][1];

	point refpointin = {fit->wcsdata.crpix[0], fit->wcsdata.crpix[1]};
	cvRotateImageRefPoint(fit, center, angle, cropped, refpointin, &refpointout); 

	fit->wcsdata.pc[0][0] = pc1_1;
	fit->wcsdata.pc[0][1] = pc1_2;
	fit->wcsdata.pc[1][0] = pc2_1;
	fit->wcsdata.pc[1][1] = pc2_2;
	fit->wcsdata.crpix[0] = refpointout.x;
	fit->wcsdata.crpix[1] = refpointout.y;

	print_updated_wcs_data(fit);
}

void crop_astrometry_data(fits *fit, point shift) {
	fit->wcsdata.crpix[0] -= shift.x;
	fit->wcsdata.crpix[1] -= shift.y;

	print_updated_wcs_data(fit);
	load_WCS_from_memory(fit); //need to update before pix2wcs - will be called once again by the crop function to update ra/dec

	center2wcs(fit, &fit->wcsdata.ra, &fit->wcsdata.dec);
	if (fit->wcsdata.ra != -1.) update_coords(); //to have plate solve window well-centered in case of subsequent call
}

void wcs_cd_to_pc(double cd[2][2], double pc[2][2], double cdelt[2]) {
	extract_cdelt_from_cd(cd[0][0], cd[0][1], cd[1][0], cd[1][1], &cdelt[0], &cdelt[1]);

	pc[0][0] = cd[0][0] / cdelt[0];
	pc[0][1] = cd[0][1] / cdelt[0];
	pc[1][0] = cd[1][0] / cdelt[1];
	pc[1][1] = cd[1][1] / cdelt[1];
}

void wcs_pc_to_cd(double pc[2][2], double cdelt[2], double cd[2][2]) {
	cd[0][0] = pc[0][0] * cdelt[0];
	cd[0][1] = pc[0][1] * cdelt[0];
	cd[1][0] = pc[1][0] * cdelt[1];
	cd[1][1] = pc[1][1] * cdelt[1];
}

gpointer match_catalog(gpointer p) {
	struct astrometry_data *args = (struct astrometry_data *) p;
	GError *error = NULL;
	psf_star **cstars;
	int n_fit = 0, n_cat = 0;
	Homography H = { 0 };
	int nobj = AT_MATCH_CATALOG_NBRIGHT;
	int max_trials = 0;
	GFile *catalog;
	GInputStream *input_stream;
	s_star star_list_A, star_list_B;
	args->catalogStars = NULL;
	args->message = NULL;

	if (args->downsample) {
		args->fit_backup = calloc(1, sizeof(fits));
		copyfits(args->fit, args->fit_backup, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
		cvResizeGaussian(args->fit, DOWNSAMPLE_FACTOR * args->fit->rx, DOWNSAMPLE_FACTOR * args->fit->ry, OPENCV_AREA);
	}

	if (!args->manual) {
		com.stars = peaker(args->fit, 0, &com.starfinder_conf, &n_fit, &(args->solvearea), FALSE, FALSE); // TODO: use good layer
	} else {
		if (com.stars)
			while (com.stars[n_fit++]);
	}

	if (!com.stars || n_fit < AT_MATCH_STARTN_LINEAR) {
		args->message = g_strdup_printf(_("There are not enough stars picked in the image. "
				"At least %d stars are needed."), AT_MATCH_STARTN_LINEAR);
		siril_log_message("%s\n", args->message);
		args->ret = 1;
		siril_add_idle(end_plate_solver, args);
		return GINT_TO_POINTER(1);
	}
	if (fabs(args->xoffset) > 0.0 || fabs(args->yoffset) > 0.0 ) max_trials = 20; //retry to converge if solve is done at an offset from the center

	cstars = new_fitted_stars(MAX_STARS);
	if (cstars == NULL) {
		PRINT_ALLOC_ERR;
		args->ret = 1;
		siril_add_idle(end_plate_solver, args);
		return GINT_TO_POINTER(1);
	}

	/* project and open the file */
	args->catalogStars = project_catalog(args->catalog_name, args->cat_center);
	if (!args->catalogStars) {
		siril_world_cs_unref(args->cat_center);
		siril_message_dialog(GTK_MESSAGE_ERROR, _("No projection"), _("Cannot project the star catalog."));
		return GINT_TO_POINTER(1);
	}
	catalog = g_file_new_for_path(args->catalogStars);
	input_stream = (GInputStream*) g_file_read(catalog, NULL, &error);
	if (input_stream == NULL) {
		if (error != NULL) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
		}
		free_fitted_stars(cstars);
		args->ret = 1;
		siril_add_idle(end_plate_solver, args);
		g_object_unref(catalog);
		return GINT_TO_POINTER(1);
	}

	n_cat = read_catalog(input_stream, cstars, args->onlineCatalog);

	/* make sure that arrays are not too small
	* make  sure that the max of stars is BRIGHTEST_STARS */
	int n = min(min(n_fit, n_cat), BRIGHTEST_STARS);

	double scale_min = args->scale - 0.2;
	double scale_max = args->scale + 0.2;
	args->ret = 1;
	int attempt = 1;
	while (args->ret && attempt < NB_OF_MATCHING_TRY) {
		args->ret = new_star_match(com.stars, cstars, n, nobj, scale_min,
				scale_max, &H, args->for_photometry_cc,
				FULLAFFINE_TRANSFORMATION, &star_list_A, &star_list_B);
		if (attempt == 1) {
			scale_min = -1.0;
			scale_max = -1.0;
		} else {
			nobj += 50;
		}
		attempt++;
	}
	if (!args->ret) {
		double conv = DBL_MAX;
		image_solved *solution = malloc(sizeof(image_solved));

		memcpy(&solution->H, &H, sizeof(Homography));
		solution->px_cat_center = siril_world_cs_ref(args->cat_center);
		/* we only want to compare with linear function
		* Maybe one day we will apply match with homography matrix
		*/
		TRANS trans = H_to_linear_TRANS(H);
		if (check_affine_TRANS_sanity(trans)) {
			double ra0, dec0;
			point image_size = { args->fit->rx, args->fit->ry };

			solution->size = image_size;
			solution->crpix[0] = ((image_size.x - 1) / 2.0);
			solution->crpix[1] = ((image_size.y - 1) / 2.0);

			apply_match(solution->px_cat_center, solution->crpix, trans, &ra0, &dec0);
			int num_matched = H.pair_matched;
			int trial = 0;

			while ((conv > CONV_TOLERANCE) && (trial < max_trials) && (!args->ret)){
				double rainit = siril_world_cs_get_alpha(args->cat_center);
				double decinit = siril_world_cs_get_delta(args->cat_center);
				double orig_ra0 = ra0;
				double orig_dec0 = dec0;

				deproject_starlist(num_matched, &star_list_B, rainit, decinit, 1);
				siril_debug_print("Deprojecting from: alpha: %s, delta: %s\n",
						siril_world_cs_alpha_format(args->cat_center, "%02d %02d %.3lf"),
						siril_world_cs_delta_format(args->cat_center, "%c%02d %02d %.3lf"));
				args->cat_center = siril_world_cs_new_from_a_d(ra0, dec0);
				solution->px_cat_center = siril_world_cs_new_from_a_d(ra0, dec0);

				project_starlist(num_matched, &star_list_B, ra0, dec0, 1);
				siril_debug_print("Reprojecting to: alpha: %s, delta: %s\n",
						siril_world_cs_alpha_format(args->cat_center, "%02d %02d %.3lf"),
						siril_world_cs_delta_format(args->cat_center, "%c%02d %02d %.3lf"));
				solution->pixel_size = args->pixel_size;

				double scaleX = sqrt(solution->H.h00 * solution->H.h00 + solution->H.h01 * solution->H.h01);
				double scaleY = sqrt(solution->H.h10 * solution->H.h10 + solution->H.h11 * solution->H.h11);
				double resolution = (scaleX + scaleY) * 0.5; // we assume square pixels
				solution->focal = RADCONV * solution->pixel_size / resolution;
				siril_debug_print("Current focal: %0.2fmm\n", solution->focal);
				
				if (atPrepareHomography(num_matched, &star_list_A, num_matched, &star_list_B, &H, FALSE, FULLAFFINE_TRANSFORMATION)){
					siril_log_color_message(_("Updating homography failed.\n"), "red");
					args->ret = 1;
				}
				trans = H_to_linear_TRANS(H);
				memcpy(&solution->H, &H, sizeof(Homography));
				apply_match(solution->px_cat_center, solution->crpix, trans, &ra0, &dec0);

				conv = fabs((dec0 - orig_dec0) / orig_dec0) + fabs((ra0 - orig_ra0) / orig_ra0);

				trial += 1;
			} 
			solution->pixel_size = args->pixel_size;


			double scaleX = sqrt(solution->H.h00 * solution->H.h00 + solution->H.h01 * solution->H.h01);
			double scaleY = sqrt(solution->H.h10 * solution->H.h10 + solution->H.h11 * solution->H.h11);
			double resolution = (scaleX + scaleY) * 0.5; // we assume square pixels
			solution->focal = RADCONV * solution->pixel_size / resolution;

			solution->image_center = siril_world_cs_new_from_a_d(ra0, dec0);
			if (max_trials == 0) {
				siril_debug_print("Converged to: alpha: %0.8f, delta: %0.8f\n", ra0, dec0);
			} else if (trial == max_trials) {
				siril_debug_print("No convergence found: alpha: %0.8f, delta: %0.8f\n", ra0, dec0);
			} else {
				siril_debug_print("Converged to: alpha: %0.8f, delta: %0.8f at iteration #%d\n", ra0, dec0, trial);
			}

			double scalefactor = (args->downsample) ? 1.0 / DOWNSAMPLE_FACTOR : 1.0;
			if (args->downsample) {
				solution->size.x *= scalefactor;
				solution->size.y *= scalefactor;
				solution->focal *= scalefactor;
				solution->crpix[0] = ((image_size.x - 1) / 2.0);
				solution->crpix[1] = ((image_size.y - 1) / 2.0);
			}

			/* compute cd matrix */
			double ra7, dec7, delta_ra;

			/* first, convert center coordinates from deg to rad: */
			dec0 *= DEGTORAD;
			ra0 *= DEGTORAD;

			/* make 1 step in direction crpix1 */
			double crpix1[] = { solution->crpix[0] + 1.0 / scalefactor, solution->crpix[1] };
			apply_match(solution->px_cat_center, crpix1, trans, &ra7, &dec7);

			dec7 *= DEGTORAD;
			ra7 *= DEGTORAD;

			delta_ra = ra7 - ra0;
			if (delta_ra > +M_PI)
				delta_ra = 2.0 * M_PI - delta_ra;
			if (delta_ra < -M_PI)
				delta_ra = delta_ra - 2.0 * M_PI;
			double cd1_1 = (delta_ra) * cos(dec0) * RADTODEG;
			double cd2_1 = (dec7 - dec0) * RADTODEG;

			/* make 1 step in direction crpix2
			* WARNING: we use -1 because of the Y axis reversing */
			double crpix2[] = { solution->crpix[0], solution->crpix[1] - 1.0 / scalefactor };
			apply_match(solution->px_cat_center, crpix2, trans, &ra7, &dec7);

			dec7 *= DEGTORAD;
			ra7 *= DEGTORAD;

			delta_ra = ra7 - ra0;
			if (delta_ra > +M_PI)
				delta_ra = 2.0 * M_PI - delta_ra;
			if (delta_ra < -M_PI)
				delta_ra = delta_ra - 2.0 * M_PI;
			double cd1_2 = (delta_ra) * cos(dec0) * RADTODEG;
			double cd2_2 = (dec7 - dec0) * RADTODEG;

			// saving state for undo before modifying fit structure
			const char *undo_str = args->for_photometry_cc ? _("Photometric CC") : _("Plate Solve");
			fits *undo_fit = args->downsample ? args->fit_backup : args->fit;

			undo_save_state(undo_fit, undo_str);

			/**** Fill wcsdata fit structure ***/

			args->fit->wcsdata.equinox = 2000.0;
			args->fit->focal_length = solution->focal;
			args->fit->pixel_size_x = args->fit->pixel_size_y = solution->pixel_size;
			solution->crpix[0] *= scalefactor;
			solution->crpix[1] *= scalefactor;

			args->fit->wcsdata.crpix[0] = solution->crpix[0];
			args->fit->wcsdata.crpix[1] = solution->crpix[1];
			args->fit->wcsdata.crval[0] = ra0 * RADTODEG;
			args->fit->wcsdata.crval[1] = dec0 * RADTODEG;

			args->fit->wcsdata.ra = siril_world_cs_get_alpha(solution->image_center);
			args->fit->wcsdata.dec = siril_world_cs_get_delta(solution->image_center);

			args->fit->wcsdata.pltsolvd = TRUE;
			g_snprintf(args->fit->wcsdata.pltsolvd_comment, 21, "Siril internal solver");

			gchar *ra = siril_world_cs_alpha_format(solution->image_center, "%02d %02d %.3lf");
			gchar *dec = siril_world_cs_delta_format(solution->image_center, "%c%02d %02d %.3lf");

			g_sprintf(args->fit->wcsdata.objctra, "%s", ra);
			g_sprintf(args->fit->wcsdata.objctdec, "%s", dec);

			g_free(ra);
			g_free(dec);

			double cdelt1, cdelt2;

			extract_cdelt_from_cd(cd1_1, cd1_2, cd2_1, cd2_2, &cdelt1, &cdelt2);

			args->fit->wcsdata.cdelt[0] = cdelt1;
			args->fit->wcsdata.cdelt[1] = cdelt2;

			/* PC + CDELT seems to be the preferred approach
			 * according to Calabretta private discussion
			 *
			 *    |cd11 cd12|  = |cdelt1      0| * |pc11 pc12|
			 *    |cd21 cd22|    |0      cdelt2|   |pc21 pc22|
			 */

			args->fit->wcsdata.pc[0][0] = cd1_1 / cdelt1;
			args->fit->wcsdata.pc[0][1] = cd1_2 / cdelt1;
			args->fit->wcsdata.pc[1][0] = cd2_1 / cdelt2;
			args->fit->wcsdata.pc[1][1] = cd2_2 / cdelt2;

			siril_debug_print("****Solution found: WCS data*************\n");
			siril_debug_print("crpix1 = %*.12e\n", 20, solution->crpix[0]);
			siril_debug_print("crpix2 = %*.12e\n", 20, solution->crpix[1]);
			siril_debug_print("crval1 = %*.12e\n", 20, ra0 * RADTODEG);
			siril_debug_print("crval2 = %*.12e\n", 20, dec0 * RADTODEG);
			siril_debug_print("cdelt1 = %*.12e\n", 20, cdelt1);
			siril_debug_print("cdelt2 = %*.12e\n", 20, cdelt2);
			siril_debug_print("pc1_1  = %*.12e\n", 20, args->fit->wcsdata.pc[0][0]);
			siril_debug_print("pc1_2  = %*.12e\n", 20, args->fit->wcsdata.pc[0][1]);
			siril_debug_print("pc2_1  = %*.12e\n", 20, args->fit->wcsdata.pc[1][0]);
			siril_debug_print("pc2_2  = %*.12e\n", 20, args->fit->wcsdata.pc[1][1]);
			siril_debug_print("******************************************\n");
			
		} else {
			args->ret = 1;
		}
		args->solution = solution;
	}
	/* free data */
	if (n_cat > 0) free_fitted_stars(cstars);
	g_object_unref(input_stream);
	g_object_unref(catalog);
	siril_add_idle(end_plate_solver, args);
	return GINT_TO_POINTER(args->ret);
}


void open_astrometry_dialog() {
	if (single_image_is_loaded() || sequence_is_loaded()) {
		initialize_ips_dialog();
		siril_open_dialog("ImagePlateSolver_Dial");
	}
}

gchar *search_in_catalogs(const gchar *object) {
	GString *string_url;
	gchar *url, *result, *name;

	set_cursor_waiting(TRUE);

	name = g_utf8_strup(object, -1);

	string_url = g_string_new(CDSSESAME);
	string_url = g_string_append(string_url, "/-oI/A?");
	string_url = g_string_append(string_url, name);
	url = g_string_free(string_url, FALSE);

	gchar *cleaned_url = url_cleanup(url);

	result = fetch_url(cleaned_url);

	set_cursor_waiting(FALSE);
	g_free(cleaned_url);
	g_free(url);
	g_free(name);

	return result;
}

int fill_plate_solver_structure(struct astrometry_data *args) {
	double fov, px_size, scale, m, usedfov, scalefactor;
	SirilWorldCS *catalog_center;
	rectangle croparea = { 0 };

	px_size = get_pixel();
	scale = get_resolution(get_focal(), px_size);
	
	args->autocrop = is_autocrop_activated();
	args->manual = is_detection_manual();
	args->downsample = is_downsample_activated();
	args->use_cache = is_cache_activated();
	args->fit = &gfit;
	scalefactor = args->downsample ? DOWNSAMPLE_FACTOR : 1.0;
	
	if (!args->manual) {
		// first checking if there is a selection or if the full field is to be used
		if (com.selection.w != 0 && com.selection.h != 0) {
			memcpy(&croparea, &com.selection, sizeof(rectangle));
		} else {
			croparea.x = 0;
			croparea.y = 0;
			croparea.w = args->fit->rx;
			croparea.h = args->fit->ry;
		}
		fov = get_fov(scale, max(croparea.w, croparea.h));

		// then apply or not autocropping to 5deg (300 arcmin)
		usedfov = (args->autocrop) ?  min(fov, 300.) : fov;
		args->cropfactor = (usedfov < fov) ? usedfov / fov : 1.0; // to avoid cropping on rounding errors
		if (args->cropfactor != 1.0) {
			croparea.x += (int) ((croparea.w - croparea.w * args->cropfactor) / 2);
			croparea.y += (int) ((croparea.h - croparea.h * args->cropfactor) / 2);
			croparea.w = (int) (args->cropfactor * croparea.w);
			croparea.h = (int) (args->cropfactor * croparea.h);
           // TODO calc center offset if need
			siril_log_message(_("Auto-cropped factor: %.2f\n"), args->cropfactor);
		}
		siril_log_message(_("Solving on selected area: %d %d %d %d \n"), croparea.x, croparea.y, croparea.w, croparea.h);

		if (com.selection.w != 0 && com.selection.h != 0) {
			args->xoffset = (double) croparea.x + 0.5 * (double) croparea.w - 0.5 * (double) args->fit->rx;
			args->yoffset = (double) croparea.y + 0.5 * (double) croparea.h - 0.5 * (double) args->fit->ry;
		} else {
			args->xoffset = 0.0;
			args->yoffset = 0.0;
		}

		if (args->downsample) {
			croparea.w *= scalefactor;
			croparea.h *= scalefactor;
			croparea.x *= scalefactor;
			croparea.y *= scalefactor;
			args->xoffset *= scalefactor;
			args->yoffset *= scalefactor;
		}
	} else { //stars manual selection - use full field centered
		fov = get_fov(scale, max(args->fit->rx, args->fit->ry));
		usedfov = fov;
		args->cropfactor = 1.0 ;  //just in case
		args->xoffset = 0.0;
		args->yoffset = 0.0;
	}
	memcpy(&(args->solvearea), &croparea, sizeof(rectangle));
	m = get_mag_limit(usedfov * CROP_ALLOWANCE);
	catalog_center = get_center_of_catalog();

	if (siril_world_cs_get_alpha(catalog_center) == 0.0 && siril_world_cs_get_delta(catalog_center) == 0.0) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("No coordinates"), _("Please enter object coordinates."));
		return 1;
	}

	/* Filling structure */
	args->onlineCatalog = args->for_photometry_cc ? get_photometry_catalog() : get_online_catalog(usedfov * CROP_ALLOWANCE, m);
	GFile *catalog_name = download_catalog(args->use_cache, args->onlineCatalog, catalog_center, usedfov * CROP_ALLOWANCE, m);
	if (!catalog_name) {
		siril_world_cs_unref(catalog_center);
		siril_message_dialog(GTK_MESSAGE_ERROR, _("No catalog"), _("Cannot download the online star catalog."));
		return 1;
	}
	args->cat_center = catalog_center; //for projection later on
	args->catalog_name = catalog_name;
	args->scale = scale;
	args->pixel_size = px_size;
	args->flip_image = flip_image_after_ps();

	return 0;
}

gboolean confirm_delete_wcs_keywords(fits *fit) {
	gboolean erase = TRUE;

	if (fit->wcsdata.equinox > 0.0) {
		erase = siril_confirm_dialog(_("Astrometric solution detected"),
				_("The astrometric solution contained in "
				"the image will be erased by the geometric transformation and no undo "
				"will be possible."), _("Process"));
	}
	return erase;
}

void invalidate_WCS_keywords(fits *fit) {
	if (fit->wcsdata.equinox > 0.0) {
		memset(&fit->wcsdata, 0, sizeof(fit->wcsdata));
	}
	free_wcs(fit);
	if (!com.headless) {
		update_MenuItem();
	}
}

/** some getters and setters */

SirilWorldCS *get_image_solved_px_cat_center(image_solved *image) {
	return image->px_cat_center;
}

SirilWorldCS *get_image_solved_image_center(image_solved *image) {
	return image->image_center;
}

void set_focal_and_pixel_pitch() {
	GtkEntry *focal, *pitch;
	gchar *f_str, *p_str;

	focal = GTK_ENTRY(lookup_widget("GtkEntry_IPS_focal"));
	pitch = GTK_ENTRY(lookup_widget("GtkEntry_IPS_pixels"));

	f_str = g_strdup_printf("%.1lf", com.pref.focal);
	p_str = g_strdup_printf("%.2lf", com.pref.pitch);

	gtk_entry_set_text(focal, f_str);
	gtk_entry_set_text(pitch, p_str);

	g_free(f_str);
	g_free(p_str);
}
