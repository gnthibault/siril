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
#include "gui/message_dialog.h"
#include "core/processing.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"

#ifdef HAVE_LIBCURL

#include "gui/PSF_list.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/plateSolver.h"
#include "registration/matching/match.h"
#include "registration/matching/apply_match.h"
#include "registration/matching/misc.h"
#include "registration/matching/atpmatch.h"
#include "registration/matching/project_coords.h"

#define RADCONV ((3600.0 * 180.0) / M_PI) / 1.0E3

#undef DEBUG           /* get some of diagnostic output */

void on_GtkTreeViewIPS_cursor_changed(GtkTreeView *tree_view,
		gpointer user_data);

static struct object platedObject[RESOLVER_NUMBER];
static GtkListStore *list_IPS = NULL;
static int selectedItem = -1;
static image_solved is_result;

static RA convert_ra(double var) {
	RA ra;
	ra.hour = (int)(var / 15.0);
	ra.min = (int)(((var / 15.0) - ra.hour) * 60.0);
	ra.sec = ((((var / 15.0) - ra.hour) * 60.0) - ra.min) * 60.0;
	return ra;
}

static double ra_to_x(RA ra) {
	return ra.hour * 15.0 + ra.min * 15.0 / 60.0 + ra.sec * 15.0 / 3600.0;
}

static DEC convert_dec(double var) {
	DEC dec;
	dec.degree = (int) var;
	dec.min = abs((int) ((var - dec.degree) * 60.0));
	dec.sec = (fabs((var - dec.degree) * 60.0) - dec.min) * 60.0;
	return dec;
}

static double dec_to_y(DEC dec) {
	if (dec.degree > 0) {
		return ((dec.sec / 3600.0) + (dec.min / 60.0) + dec.degree);
	} else {
		return (-(dec.sec / 3600.0) - (dec.min / 60.0) + dec.degree);
	}
}

static void deg_to_HMS(double var, gchar *type, gchar *HMS) {
	if (!strncasecmp(type, "ra", 2)) {
		char rs = ' ';
		int raH, raM;
		double raS;

		if (var < 0) rs = '-';
		var = fabs(var);
		raH = (int)(var / 15.0);
		raM = (int)(((var / 15.0) - raH) * 60.0);
		raS = ((((var / 15.0) - raH) * 60.0) - raM) * 60.0;
		g_snprintf(HMS, 256, "%c%02d %02d %.3lf", rs, raH, raM, raS);
	} else if (!strncasecmp(type, "dec", 2)) {
		char ds = ' ';
		int deg, decM;
		double decS;

		if (var < 0) ds = '-';
		var = fabs(var);
		deg = (int) var;
		decM = abs((int) ((var - deg) * 60));
		decS = (fabs((var - deg) * 60) - decM) * 60;
		g_snprintf(HMS, 256, "%c%02d %02d %.3lf", ds, deg, decM, decS);
	}
}

static void fov_in_DHMS(double var, gchar *fov) {
	int deg, decM;
	double decS;

	if (var < 0) {
		fprintf(stdout, "fov_in_DHMS: negative value, should not happend\n");
		return;
	}
	deg = (int) var;
	decM = abs((int) ((var - deg) * 60));
	decS = (fabs((var - deg) * 60) - decM) * 60;
	if (deg > 0)
		g_snprintf(fov, 256, "%02dd %02d\' %.2lf\"", deg, decM, decS);
	else if (decM > 0)
		g_snprintf(fov, 256, "%02d\' %.2lf\"", decM, decS);
	else if (decS > 0.0)
		g_snprintf(fov, 256, "%.2lf\"", decS);
}

static int parse_curl_buffer(char *buffer, struct object *obj) {
	char **token, **fields;
	point center;
	int nargs, i = 0, resolver = -1;

	token = g_strsplit(buffer, "\n", -1);
	nargs = g_strv_length(token);

	while (i < nargs) {
		if (g_strrstr (token[i], "NED")) {
			resolver = RESOLVER_NED;
		}
		else if (g_strrstr (token[i], "Simbad")) {
			resolver = RESOLVER_SIMBAD;
		}
		else if (g_strrstr(token[i], "VizieR")) {
			resolver = RESOLVER_VIZIER;
		}
		else if (g_str_has_prefix (token[i], "%J ")) {
			fields = g_strsplit(token[i], " ", -1);
			sscanf(fields[1], "%lf", &center.x);
			sscanf(fields[2], "%lf", &center.y);
			if (resolver != -1) {
				/* RA coordinates */
				platedObject[resolver].RA = convert_ra(center.x);

				/* Dec coordinates */
				platedObject[resolver].Dec = convert_dec(center.y);

				/* others */
				platedObject[resolver].imageCenter = center;
				platedObject[resolver].south = (center.y < 0.0);
			}
			g_strfreev(fields);
		}
		else if (g_str_has_prefix (token[i], "%I.0 ")) {
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
	int i;
	for (i = 0; i < RESOLVER_NUMBER; i++) {
		if (platedObject[i].name) {
			free(platedObject[i].name);
			platedObject[i].name = NULL;
		}
	}
}

static double get_focal() {
	GtkEntry *focal_entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_focal"));
	const gchar *value = gtk_entry_get_text(focal_entry);

	return atof(value);
}

/* get pixel in µm */
static double get_pixel() {
	GtkEntry *pixel_entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_pixels"));
	const gchar *value = gtk_entry_get_text(pixel_entry);

	return atof(value);
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

static point get_center_of_catalog() {
	GtkSpinButton *GtkSpinIPS_RA_h, *GtkSpinIPS_RA_m;
	GtkSpinButton *GtkSpinIPS_Dec_deg, *GtkSpinIPS_Dec_m;
	GtkEntry *GtkEntryIPS_RA_s, *GtkEntryIPS_Dec_s;
	GtkToggleButton *GtkCheckButtonIPS_S;
	RA ra;
	DEC dec;
	point result;

	/* get RA center */
	GtkSpinIPS_RA_h = GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_RA_h"));
	GtkSpinIPS_RA_m = GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_RA_m"));
	GtkEntryIPS_RA_s = GTK_ENTRY(lookup_widget("GtkEntryIPS_RA_s"));

	ra.hour = gtk_spin_button_get_value_as_int(GtkSpinIPS_RA_h);
	ra.min = gtk_spin_button_get_value_as_int(GtkSpinIPS_RA_m);
	ra.sec = atof(gtk_entry_get_text(GtkEntryIPS_RA_s));

	/* get Dec center */
	GtkSpinIPS_Dec_deg = GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_Dec_deg"));
	GtkSpinIPS_Dec_m = GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_Dec_m"));
	GtkEntryIPS_Dec_s = GTK_ENTRY(lookup_widget("GtkEntryIPS_Dec_s"));
	GtkCheckButtonIPS_S = GTK_TOGGLE_BUTTON(lookup_widget("GtkCheckButtonIPS_S"));

	dec.degree = gtk_spin_button_get_value_as_int(GtkSpinIPS_Dec_deg);
	dec.min = gtk_spin_button_get_value_as_int(GtkSpinIPS_Dec_m);
	dec.sec = atof(gtk_entry_get_text(GtkEntryIPS_Dec_s));
	if (gtk_toggle_button_get_active(GtkCheckButtonIPS_S)) {
		dec.degree = -dec.degree;
	}

	/* convert */
	result.x = ra_to_x(ra);
	result.y = dec_to_y(dec);

	return result;
}

static gboolean is_detection_manual() {
	GtkToggleButton *button;

	button = GTK_TOGGLE_BUTTON(lookup_widget("checkButton_IPS_manual"));
	return gtk_toggle_button_get_active(button);
}

static gchar *get_catalog_url(point center, double mag_limit, double dfov, int type) {
	GString *url;
	gchar *coordinates;
	gchar *mag;
	gchar *fov;

	coordinates = g_strdup_printf("%lf+%lf", center.x, center.y);
	mag = g_strdup_printf("%2.2lf", mag_limit);
	fov = g_strdup_printf("%2.1lf", dfov / 2);

	url = g_string_new("http://vizier.u-strasbg.fr/viz-bin/asu-tsv?-source=");
	switch (type) {
	case NOMAD:
		url = g_string_append(url, "NOMAD&-out.meta=-h-u-D&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-out=%20RAJ2000%20DEJ2000%20Vmag%20Bmag");
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&Vmag=<");
		url = g_string_append(url, mag);
		break;
	default:
	case TYCHO2:
		url = g_string_append(url, "I/259/tyc2&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-c.u=deg&-out=%20RAmdeg%20DEmdeg%20VTmag%20BTmag");
		url = g_string_append(url, "&-out.max=200000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&VTmag=<");
		url = g_string_append(url, mag);
		break;
	}

	g_free(coordinates);
	g_free(mag);
	g_free(fov);

	return g_string_free(url, FALSE);
}

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
	char *result;
	long code;
	int retries;
	unsigned int s;

	printf("fetch_url(): %s\n", url);

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
			printf("Fetch failed with code %ld for URL %s\n", code, url);
		}
	}

	if (!result)
		free(content->data);

	free(content);

	return result;
}

static online_catalog get_online_catalog() {
	GtkComboBox *box;
	int ret;

	box = GTK_COMBO_BOX(lookup_widget("ComboBoxIPSCatalog"));
	ret = gtk_combo_box_get_active(box);
	return (ret < 0 ? TYCHO2 : ret);
}

static gchar *download_catalog(online_catalog onlineCatalog) {
	gchar *url;
	char *buffer = NULL;
	FILE *catalog = NULL;
	FILE *fproj = NULL;
	gchar *filename, *foutput;
	int index = selectedItem;
	point catalog_center;

	if (index < 0)
		return NULL;

	double fov = get_fov(get_resolution(get_focal(), get_pixel()),
			gfit.ry > gfit.rx ? gfit.ry : gfit.rx);

	catalog_center = get_center_of_catalog();

	/* ------------------- get Vizier catalog in catalog.dat -------------------------- */

	url = get_catalog_url(catalog_center, get_mag_limit(fov), fov,
			onlineCatalog);

	filename = g_build_filename(g_get_tmp_dir(), "catalog.dat", NULL);
	catalog = g_fopen(filename, "w+t");
	if (catalog == NULL) {
		fprintf(stderr, "plateSolver: Cannot open catalogue\n");
		return NULL;
	}
	buffer = fetch_url(url);
	fprintf(catalog, "%s", buffer);
	g_free(url);
	free(buffer);
	fclose(catalog);

	/* -------------------------------------------------------------------------------- */

	/* --------- Project coords of Vizier catalog and save it into catalog.proj ------- */

	foutput = g_build_filename(g_get_tmp_dir(), "catalog.proj", NULL);
	fproj = g_fopen(foutput, "w+t");
	if (fproj == NULL) {
		fprintf(stderr, "plateSolver: Cannot open fproj\n");
		g_free(foutput);
		return NULL;
	}

	convert_catalog_coords(filename, platedObject[index].imageCenter, fproj);
	fclose(fproj);

	/* -------------------------------------------------------------------------------- */

	is_result.px_cat_center = platedObject[index].imageCenter;

	g_free(filename);
	return foutput;
}

/*********
 *
 */

enum {
	COLUMN_RESOLVER,		// string
	COLUMN_NAME,		// string
	N_COLUMNS
};

static void get_list_IPS() {
	if (list_IPS == NULL)
		list_IPS = GTK_LIST_STORE(gtk_builder_get_object(builder, "liststoreIPS"));
}

static void clear_all_objects() {
	static GtkTreeSelection *selection;

	gtk_list_store_clear(list_IPS);
	selectedItem = -1;

	selection = GTK_TREE_SELECTION(gtk_builder_get_object(builder, "GtkTreeSelectionIPS"));

	gtk_tree_selection_unselect_all (selection);
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

static void update_coordinates() {
	gchar *RA_sec, *Dec_sec;

	int index = selectedItem;
	if (index < 0)
		return;

	RA_sec = g_strdup_printf("%6.4lf", platedObject[index].RA.sec);
	Dec_sec = g_strdup_printf("%6.4lf", platedObject[index].Dec.sec);

	gtk_toggle_button_set_active(
			GTK_TOGGLE_BUTTON(lookup_widget("GtkCheckButtonIPS_S")),
			platedObject[index].south);

	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_RA_h")),
			platedObject[index].RA.hour);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_RA_m")),
			platedObject[index].RA.min);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("GtkEntryIPS_RA_s")), RA_sec);

	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_Dec_deg")),
			abs(platedObject[index].Dec.degree));
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_Dec_m")),
			platedObject[index].Dec.min);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("GtkEntryIPS_Dec_s")), Dec_sec);

	g_free(RA_sec);
	g_free(Dec_sec);
}

static gboolean has_any_keywords() {
	return (gfit.focal_length > 0.0 ||
			gfit.pixel_size_x > 0.0 ||
			gfit.pixel_size_y > 0.0);
}

static void update_pixel_size() {
	GtkEntry *entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_pixels"));
	gchar *cpixels;
	double pixel;

	pixel = gfit.pixel_size_x > gfit.pixel_size_y ? gfit.pixel_size_x : gfit.pixel_size_y;

	if (pixel > 0.0) {
		cpixels = g_strdup_printf("%1.2lf", pixel);
		gtk_entry_set_text(entry, cpixels);
		g_free(cpixels);
	}
}

static void update_focal() {
	GtkEntry *entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_focal"));
	gchar *cfocal;
	double focal;

	focal = gfit.focal_length;

	if (focal > 0.0) {
		cfocal = g_strdup_printf("%.3lf", focal);
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

static void update_IPS_GUI() {
	/* update all fields. Resolution is updated as well
	 thanks to the Entry and combo changed signal
	  */
	update_focal();
	update_pixel_size();
}

static void update_gfit(image_solved image) {
	gfit.focal_length = image.focal;
	gfit.pixel_size_x = gfit.pixel_size_y = image.pixel_size;
	gfit.wcs.crpix1 = image.x;
	gfit.wcs.crpix2 = image.y;
	gfit.wcs.crval1 = image.ra;
	gfit.wcs.crval2 = image.dec;
	gfit.wcs.equinox = 2000;
	deg_to_HMS(image.ra, "ra", gfit.wcs.objctra);
	deg_to_HMS(image.dec, "dec", gfit.wcs.objctdec);
}

static void print_platesolving_results(Homography H, image_solved image) {
	double rotation, det, scaleX, scaleY;
	double inliers;
	char RA[256] = { 0 };
	char DEC[256] = { 0 };
	char field_x[256] = { 0 };
	char field_y[256] = { 0 };

	/* Matching information */
	siril_log_message(_("%d pair matches.\n"), H.pair_matched);
	inliers = 1.0 - ((((double) H.pair_matched - (double) H.Inliers)) / (double) H.pair_matched);
	siril_log_message(_("Inliers:%*.3f\n"), 14, inliers);

	/* Plate Solving */
	scaleX = sqrt(H.h00 * H.h00 + H.h01 * H.h01);
	scaleY = sqrt(H.h10 * H.h10 + H.h11 * H.h11);
	image.resolution = (scaleX + scaleY) * 0.5; // we assume square pixels
	siril_log_message(_("Resolution:%*.3lf arcsec/px\n"), 11, image.resolution);

	/* rotation */
	rotation = atan2(H.h00 + H.h01, H.h10 + H.h11) * 180 / M_PI + 135.0;
	det = (H.h00 * H.h11 - H.h01 * H.h10); // determinant of rotation matrix (ac - bd)
	/* If the determinant of the top-left 2x2 rotation matrix is > 0
	 * the transformation is orientation-preserving. */
	if (det < 0)
		rotation = -90 - rotation;
	rotation = - rotation;
	if (rotation < -180)
		rotation += 360;
	if (rotation > 180)
		rotation -= 360;
	siril_log_message(_("Rotation:%+*.2lf deg %s\n"), 12, rotation, det < 0 ? _("(flipped)") : "");
	image.focal = RADCONV * image.pixel_size / image.resolution;

	image.fov.x = get_fov(image.resolution, image.px_size.x);
	image.fov.y = get_fov(image.resolution, image.px_size.y);

	siril_log_message(_("Focal:%*.2lf mm\n"), 15, image.focal);
	siril_log_message(_("Pixel size:%*.2lf µm\n"), 10, image.pixel_size);
	fov_in_DHMS(image.fov.x / 60.0, field_x);
	fov_in_DHMS(image.fov.y / 60.0, field_y);
	siril_log_message(_("Field of view:    %s x %s\n"), field_x, field_y);
	deg_to_HMS(image.ra, "ra", RA);
	deg_to_HMS(image.dec, "dec", DEC);
	siril_log_message(_("Image center: RA: %s, DEC: %s\n"), RA, DEC);

   	update_gfit(image);
}

static int read_NOMAD_catalog(FILE *catalog, fitted_PSF **cstars, int shift_y) {
	char line[LINELEN];
	fitted_PSF *star;

	int i = 0;

	while (fgets(line, LINELEN, catalog) != NULL) {
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

		star = malloc(sizeof(fitted_PSF));
		star->xpos = x;
		star->ypos = y + shift_y; // shift_y is here because in get_stars of misc.c we have "image_size.y - s[i]->ypos"
		star->mag = Vmag;
		star->BV = n < 5 ? -99.9 : Bmag - Vmag;
		cstars[i] = star;
		cstars[i + 1] = NULL;
		i++;
	}
	sort_stars(cstars, i);

	return i;
}

static int read_TYCHO2_catalog(FILE *catalog, fitted_PSF **cstars, int shift_y) {
	char line[LINELEN] = {0};
	fitted_PSF *star;

	int i = 0;

	while (fgets(line, LINELEN, catalog) != NULL) {
		double r = 0.0, x = 0.0, y = 0.0, magA = 0.0, magB = 0.0;

		if (line[0] == COMMENT_CHAR) {
			continue;
		}
		if (is_blank(line)) {
			continue;
		}
		if (g_str_has_prefix(line, "---")) {
			continue;
		}
		if (g_str_has_prefix(line, "_r")) {
			continue;
		}
		if (g_str_has_prefix(line, "arcmin")) {
			continue;
		}
		int n = sscanf(line, "%lf %lf %lf %lf %lf", &r, &x, &y, &magA, &magB);

		star = malloc(sizeof(fitted_PSF));
		star->xpos = x;
		star->ypos = y + shift_y;
		star->mag = magA;
		cstars[i] = star;
		cstars[i + 1] = NULL;
		i++;

	}
	sort_stars(cstars, i);
	return i;
}

static int read_catalog(FILE *catalog, fitted_PSF **cstars, int shift_y,
		int type) {
	switch (type) {
	default:
	case TYCHO2:
		return read_TYCHO2_catalog(catalog, cstars, shift_y);
	case NOMAD:
		return read_NOMAD_catalog(catalog, cstars, shift_y);
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
	gboolean ok = FALSE;

	double var1 = fabs(trans.b/trans.f);
	double var2 = fabs(trans.c/trans.e);
	siril_debug_print("abs(b/f)=%lf et abs(c/e)=%lf\n", var1, var2);

	if (0.9 < var1 && var1 < 1.1) {
		if (0.9 < var2 && var2 < 1.1) {
			ok = TRUE;
		}
	}
	return ok;
}

static gboolean end_plate_solver(gpointer p) {
	struct plate_solver_data *args = (struct plate_solver_data *) p;
	stop_processing_thread();
	if (!args->manual)
		clear_stars_list();
	set_cursor_waiting(FALSE);

	if (args->ret) {
		char *title = siril_log_color_message(_("Plate Solving failed. "
				"The image could not be aligned with the reference stars.\n"), "red");
		if (!args->message) {
			args->message = g_strdup(_("This is usually because the initial parameters (pixel size, focal length, initial coordinates) "
					"are too far from the real metadata of the image."));
		}
		siril_message_dialog(GTK_MESSAGE_ERROR, title, args->message);
	} else {
		update_IPS_GUI();
		control_window_switch_to_tab(OUTPUT_LOGS);
	}
	g_free(args->catalogStars);
	g_free(args->message);
	free(args);
	update_used_memory();
	return FALSE;
}

static gpointer match_catalog(gpointer p) {
	struct plate_solver_data *args = (struct plate_solver_data *) p;
	FILE *catalog;
	fitted_PSF **cstars;
	int n_fit, n_cat, n, i = 0;
	int attempt = 1;
	point image_size = {args->fit->rx, args->fit->ry};
	Homography H = { 0 };
	int nobj = AT_MATCH_CATALOG_NBRIGHT;

	args->message = NULL;

	if (!args->manual) {
		com.stars = peaker(args->fit, 0, &com.starfinder_conf, &n_fit, NULL, FALSE); // TODO: use good layer
	} else {
		while (com.stars && com.stars[i]) {
			i++;
		}
		n_fit = i;
	}
	if (n_fit < AT_MATCH_STARTN_LINEAR) {
		args->message = g_strdup_printf(_("There are not enough stars picked in the image. "
				"At least %d stars are needed."), AT_MATCH_STARTN_LINEAR);
		siril_log_message("%s\n", args->message);
		args->ret = 1;
		siril_add_idle(end_plate_solver, args);
		return GINT_TO_POINTER(1);
	}

	cstars = malloc((MAX_STARS + 1) * sizeof(fitted_PSF *));
	if (cstars == NULL) {
		printf("Memory allocation failed: peaker\n");
		args->ret = 1;
		siril_add_idle(end_plate_solver, args);
		return GINT_TO_POINTER(1);
	}

	/* open the file */
	catalog = g_fopen(args->catalogStars, "r");
	if (catalog == NULL) {
		fprintf(stderr, "match_catalog: error opening file: %s\n", args->catalogStars);
		free(cstars);
		args->ret = 1;
		siril_add_idle(end_plate_solver, args);
		return GINT_TO_POINTER(1);
	}

	n_cat = read_catalog(catalog, cstars, image_size.y, args->onlineCatalog);

	/* make sure that arrays are not too small
	 * make  sure that the max of stars is BRIGHTEST_STARS */
	n = n_fit < n_cat ? n_fit : n_cat;
	n = n > BRIGHTEST_STARS ? BRIGHTEST_STARS : n;

	args->ret = 1;
	while (args->ret && attempt < NB_OF_MATCHING_TRY){
		args->ret = new_star_match(com.stars, cstars, n, nobj,
				args->scale - 0.2, args->scale + 0.2, &H, image_size);
		nobj += 50;
		attempt++;
	}
	if (!args->ret) {

		/* we only want to compare with linear function
		 * Maybe one day we will apply match with homography matrix
		 */
		TRANS trans = H_to_linear_TRANS(H);
		if (check_affine_TRANS_sanity(trans)) {

			is_result.px_size = image_size;
			is_result.x = image_size.x / 2.0;
			is_result.y = image_size.y / 2.0;
			is_result.pixel_size = args->pixel_size;

			apply_match(&is_result, trans);

			print_platesolving_results(H, is_result);
		} else {
			args->ret = 1;
		}

	}
	/* free data */
	if (n_cat > 0) free_fitted_stars(cstars);
	fclose(catalog);
	siril_add_idle(end_plate_solver, args);
	return GINT_TO_POINTER(args->ret);
}

static void url_encode(gchar **qry) {
	int new_string_length = 0;
	for (gchar *c = *qry; *c != '\0'; c++) {
		if (*c == ' ')
			new_string_length += 2;
		new_string_length++;
	}
	gchar *qstr = g_malloc((new_string_length + 1) * sizeof qstr[0]);
	gchar *c1, *c2;
	for (c1 = *qry, c2 = qstr; *c1 != '\0'; c1++) {
		if (*c1 == ' ') {
			c2[0] = '%';
			c2[1] = '2';
			c2[2] = '0';
			c2 += 3;
		} else {
			*c2 = *c1;
			c2++;
		}
	}
	*c2 = '\0';
	g_free(*qry);
	*qry = g_strdup(qstr);
	g_free(qstr);
}

static void search_object_in_catalogs(const gchar *object) {
	GString *url;
	gchar *gcurl, *result, *name;
	struct object obj;
	GtkTreeView *GtkTreeViewIPS;

	GtkTreeViewIPS = GTK_TREE_VIEW(lookup_widget("GtkTreeViewIPS"));

	set_cursor_waiting(TRUE);

	name = g_strdup(object);
	/* Removes leading and trailing whitespace */
	name = g_strstrip(name);
	/* replace whitespaces by %20 for html purposes */
	url_encode(&name);

	url = g_string_new(CDSSESAME);
	url = g_string_append(url, "/-oI/A?");
	url = g_string_append(url, name);
	gcurl = g_string_free(url, FALSE);

	result = fetch_url(gcurl);
	if (result) {
		parse_curl_buffer(result, &obj);
		g_signal_handlers_block_by_func(GtkTreeViewIPS, on_GtkTreeViewIPS_cursor_changed, NULL);
		add_object_to_list();
		g_signal_handlers_unblock_by_func(GtkTreeViewIPS, on_GtkTreeViewIPS_cursor_changed, NULL);
		free_Platedobject();
	}
	set_cursor_waiting(FALSE);
	g_free(gcurl);
	g_free(result);
	g_free(name);
}

static void start_image_plate_solve() {

	struct plate_solver_data *args = malloc(sizeof(struct plate_solver_data));
	set_cursor_waiting(TRUE);

	args->onlineCatalog = get_online_catalog();
	args->catalogStars = download_catalog(args->onlineCatalog);
	args->pixel_size = get_pixel();
	args->scale = get_resolution(get_focal(), get_pixel());
	args->manual = is_detection_manual();
	args->fit = &gfit;

	start_in_new_thread(match_catalog, args);
}

/*****
 * CALLBACKS FUNCTIONS
 */

void on_GtkEntry_IPS_changed(GtkEditable *editable, gpointer user_data) {
	update_resolution_field();
}

void on_menuitem_IPS_activate(GtkMenuItem *menuitem, gpointer user_data) {
	gtk_widget_show_all(lookup_widget("ImagePlateSolver_Dial"));
}

void on_buttonIPS_close_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("ImagePlateSolver_Dial"));
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
		if (!g_strcmp0(res, "NED")) {
			selectedItem = 0;
		} else if (!g_strcmp0(res, "Simbad")) {
			selectedItem = 1;
		} else if (!g_strcmp0(res, "VizieR")) {
			selectedItem = 2;
		} else {
			selectedItem = -1;
		}

		if (selectedItem >= 0)
			update_coordinates();

		g_value_unset(&value);
	}
}

void on_GtkButton_IPS_metadata_clicked(GtkButton *button, gpointer user_data) {
	if (!has_any_keywords()) {
		char *msg = siril_log_message(_("There are no keywords storred in the FITS header.\n"));
		siril_message_dialog(GTK_MESSAGE_WARNING, _("No metadata"), msg);
	} else {
		update_IPS_GUI();
	}
}

void on_GtkButtonIPS_clicked(GtkButton *button, gpointer user_data) {
	GtkEntry *entry;

	entry = GTK_ENTRY(lookup_widget("GtkSearchIPS"));
	search_object_in_catalogs(gtk_entry_get_text(GTK_ENTRY(entry)));
}

void on_buttonIPS_ok_clicked(GtkButton *button, gpointer user_data) {
	start_image_plate_solve();
}

void on_GtkSearchIPS_activate(GtkEntry *entry, gpointer user_data) {
	search_object_in_catalogs(gtk_entry_get_text(GTK_ENTRY(entry)));
}

void on_GtkCheckButton_Mag_Limit_toggled(GtkToggleButton *button,
		gpointer user_data) {
	GtkWidget *spinmag;

	spinmag = lookup_widget("GtkSpinIPS_Mag_Limit");
	gtk_widget_set_sensitive(spinmag, !gtk_toggle_button_get_active(button));
}
#endif

/******
 *
 * Public functions
 */

gboolean confirm_delete_wcs_keywords(fits *fit) {
	gboolean erase = TRUE;

	if (fit->wcs.equinox > 0) {
		erase = siril_confirm_dialog(_("Astrometric solution detected"), _("The astrometric solution contained in "
				"the image will be erased by the geometric transformation and no undo "
				"will be possible."), FALSE);
	}
	return erase;
}

void invalidate_WCS_keywords(fits *fit) {
	if (fit->wcs.equinox > 0) {
		fit->wcs.equinox = 0;
		fit->wcs.crpix1 = 0.0;
		fit->wcs.crpix2 = 0.0;
		fit->wcs.crval1 = 0.0;
		fit->wcs.crval2 = 0.0;
		memset(fit->wcs.objctra, 0, FLEN_VALUE);
		memset(fit->wcs.objctdec, 0, FLEN_VALUE);
	}
}
