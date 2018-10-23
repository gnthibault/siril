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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <curl/curl.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
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

enum {
	APPAS, NOMAD
};

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

static DEC convert_dec(double var) {
	DEC dec;
	dec.degree = (int) var;
	dec.min = fabs((int) ((var - dec.degree) * 60));
	dec.sec = (fabs((var - dec.degree) * 60) - dec.min) * 60;
	return dec;
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
		decM = fabs((int) ((var - deg) * 60));
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
	decM = fabs((int) ((var - deg) * 60));
	decS = (fabs((var - deg) * 60) - decM) * 60;
	g_snprintf(fov, 256, "%02dd %02d\' %.2lf\"", deg, decM, decS);
}

static int parse_curl_buffer(char *buffer, struct object *obj) {
	char **token, **fields;
	point center;
	gchar *object = NULL;
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
			resolver = RESOLVER_VISIER;
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


static gchar *get_catalog_url(point center, double mag_limit, double dfov, int type) {
	GString *url;
	gchar coordinates[256];
	gchar mag[6];
	gchar fov[10];

	g_snprintf(coordinates, sizeof(coordinates), "%lf+%lf", center.x, center.y);
	g_snprintf(mag, sizeof(mag), "%2.2lf", mag_limit);
	g_snprintf(fov, sizeof(fov), "%2.1lf", dfov / 2);

	url = g_string_new("http://vizier.u-strasbg.fr/viz-bin/asu-tsv?-source=");
	switch (type) {
	case NOMAD:
		url = g_string_append(url, "NOMAD&-out.meta=-h-u-D&-out.add=_r&-sort=_r");
		url = g_string_append(url, "&-out=%20RAJ2000%20DEJ2000%20Vmag%20Bmag");
		url = g_string_append(url, "&-c.rm=");
		url = g_string_append(url, fov);
		url = g_string_append(url, "&-out.max=5000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		url = g_string_append(url, "&Vmag=<");
		url = g_string_append(url, mag);
		break;
	default:
	case APPAS:
		url = g_string_append(url, "I/259/tyc2");
		url = g_string_append(url, "&-c.u=deg&-out=TYC1&-out=TYC2&-out=TYC3&-out=RAmdeg&-out=DEmdeg&-out=pmRA&-out=pmDE&-out=VTmag&-out=BTmag&-out=HIP");
		url = g_string_append(url, "&-out.max=5000");
		url = g_string_append(url, "&-c=");
		url = g_string_append(url, coordinates);
		// TODO: fov ?
		url = g_string_append(url, "&-c.r=51");
		url = g_string_append(url, "&VTmag=<");
		url = g_string_append(url, mag);
		break;
	}
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

static gchar *download_catalog() {
	gchar *url;
	char *buffer = NULL;
	FILE *catalog = NULL;
	FILE *fproj = NULL;
	gchar *filename, *foutput;
	int index = selectedItem;

	if (index < 0)
		return NULL;

	double fov = get_fov(get_resolution(get_focal(), get_pixel()),
			gfit.ry > gfit.rx ? gfit.ry : gfit.rx);

	/* ------------------- get Vizier catalog in catalog.dat -------------------------- */

	url = get_catalog_url(platedObject[index].imageCenter, get_mag_limit(fov), fov,
			NOMAD);

	filename = g_build_filename(g_get_tmp_dir(), "catalog.dat", NULL);
	catalog = g_fopen(filename, "w+t");
	if (catalog == NULL) {
		fprintf(stderr, "plateSolver: Cannot open catalog\n");
		g_free(foutput);
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
		gtk_list_store_set(list_IPS, &iter, COLUMN_RESOLVER, "SIMBAD",
				COLUMN_NAME, platedObject[RESOLVER_SIMBAD].name, -1);
	}

	if (platedObject[RESOLVER_VISIER].name) {
		gtk_list_store_append(list_IPS, &iter);
		gtk_list_store_set(list_IPS, &iter, COLUMN_RESOLVER, "VISIER",
				COLUMN_NAME, platedObject[RESOLVER_VISIER].name, -1);
	}
}


static void update_coordinates() {
	gchar RA_sec[7], Dec_sec[7];

	int index = selectedItem;
	if (index < 0)
		return;

	g_snprintf(RA_sec, sizeof(RA_sec), "%lf", platedObject[index].RA.sec);
	g_snprintf(Dec_sec, sizeof(Dec_sec), "%lf", platedObject[index].Dec.sec);

	gtk_toggle_button_set_active(
			GTK_TOGGLE_BUTTON(lookup_widget("GtkCheckButtonIPS_S")),
			platedObject[index].south);

	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_RA_h")),
			platedObject[index].RA.hour);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_RA_m")),
			platedObject[index].RA.min);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("GtkEntryIPS_RA_s")), RA_sec);

	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_Dec_deg")),
			fabs(platedObject[index].Dec.degree));
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_Dec_m")),
			platedObject[index].Dec.min);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("GtkEntryIPS_Dec_s")), Dec_sec);
}

static gboolean has_any_keywords() {
	return (gfit.focal_length > 0.0 ||
			gfit.pixel_size_x > 0.0 ||
			gfit.pixel_size_y > 0.0);
}

static void update_pixel_size() {
	GtkEntry *entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_pixels"));
	gchar cpixels[24];
	double pixel;

	pixel = gfit.pixel_size_x > gfit.pixel_size_y ? gfit.pixel_size_x : gfit.pixel_size_y;

	if (pixel > 0.0) {
		g_snprintf(cpixels, sizeof(cpixels), "%1.2lf", pixel);
		gtk_entry_set_text(entry, cpixels);
	}
}

static void update_focal() {
	GtkEntry *entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_focal"));
	gchar cfocal[24];
	double focal;

	focal = gfit.focal_length;

	if (focal > 0.0) {
		g_snprintf(cfocal, sizeof(cfocal), "%.3lf", focal);
		gtk_entry_set_text(entry, cfocal);
	}
}

static void update_resolution_field() {
	GtkEntry *entry = GTK_ENTRY(lookup_widget("GtkEntry_IPS_resolution"));
	double res = get_resolution(get_focal(), get_pixel());
	gchar cres[6];

	g_snprintf(cres, sizeof(cres), "%1.3lf", res);
	gtk_entry_set_text(entry, cres);
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
	gfit.crpix1 = image.px_center.x;
	gfit.crpix2 = image.px_center.y;
	gfit.crval1 = image.ra_center;
	gfit.crval2 = image.dec_center;
}

static void print_platesolving_results(Homography H, image_solved image) {
	double rotation, det, resolution, scaleX, scaleY, focal;
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
	resolution = (scaleX + scaleY) * 0.5; // we assume square pixels
	siril_log_message(_("Resolution:%*.3lf arcsec/px\n"), 11, resolution);
	rotation = atan2(H.h00 + H.h01, H.h10 + H.h11) * 180 / M_PI + 135.0;
	det = (H.h00 * H.h11 - H.h01 * H.h10); // determinant of rotation matrix (ac - bd)
	if (det < 0)
		rotation = -90 - rotation;
	rotation = - rotation;
	if (rotation < -180)
		rotation += 360;
	if (rotation > 180)
		rotation -= 360;
	siril_log_message(_("Rotation:%+*.2lf deg %s\n"), 12, rotation, det > 0 ? _("(flipped)") : "");
	focal = RADCONV * image.pixel_size / resolution;

	image.focal = focal;
	image.fov.x = get_fov(resolution, image.px_size.x);
	image.fov.y = get_fov(resolution, image.px_size.y);

	siril_log_message(_("Focal:%*.2lf mm\n"), 15, focal);
	siril_log_message(_("Pixel size:%*.2lf µm\n"), 10, image.pixel_size);
	fov_in_DHMS(image.fov.x / 60.0, field_x);
	fov_in_DHMS(image.fov.y / 60.0, field_y);
	siril_log_message(_("Field of view:    %s x %s\n"), field_x, field_y);
	deg_to_HMS(image.ra_center, "ra", RA);
	deg_to_HMS(image.dec_center, "dec", DEC);
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
		star->xpos = -x;
		star->ypos = y - shift_y;
		star->mag = Vmag;
		star->BV = n < 5 ? -99.9 : Bmag - Vmag;
		cstars[i] = star;
		cstars[i + 1] = NULL;
		i++;
	}
	sort_stars(cstars, i);

	return i;
}

static int read_APPAS_catalog(FILE *catalog, fitted_PSF **cstars, int shift_y) {
	char line[LINELEN];
	fitted_PSF *star;
	double x, y, magA, magB, pmRA, pmDec;
	int i = 0;
	int tmp;

	while (fgets(line, LINELEN, catalog) != NULL) {
		if (line[0] == COMMENT_CHAR) {
			continue;
		}
		if (is_blank(line)) {
			continue;
		}
		if (g_str_has_prefix(line, "---")) {
			continue;
		}
		if (g_str_has_prefix(line, "TYC")) {
			continue;
		}
		int n = sscanf(line, "%d %d %d %lf %lf %lf %lf %lf %lf %d", &tmp, &tmp, &tmp, &x, &y, &pmRA, &pmDec, &magA, &magB, &tmp);
		if (n < 10)
			continue;
		star = malloc(sizeof(fitted_PSF));
		star->xpos = x;
		star->ypos = -y + shift_y;
		star->mag = magA;
		cstars[i] = star;
		cstars[i + 1] = NULL;
		i++;
	}
	sort_stars(cstars, i);
	return i;
}

static int read_catalog(FILE *catalog, fitted_PSF **cstars, int shift_y, int type) {
	switch (type) {
	default:
	case APPAS:
		return read_APPAS_catalog(catalog, cstars, shift_y);
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

static gboolean end_plate_solver(gpointer p) {
	struct plate_solver_data *args = (struct plate_solver_data *) p;
	stop_processing_thread();
	clear_stars_list();
	set_cursor_waiting(FALSE);
	if (!args->ret)
		update_IPS_GUI();
	g_free(args->catalogStars);
	free(args);
	update_used_memory();
	return FALSE;
}

static gpointer match_catalog(gpointer p) {
	struct plate_solver_data *args = (struct plate_solver_data *) p;
	FILE *catalog;
	fitted_PSF **cstars;
	int n_fit, n_cat, n;
	point image_size = {args->fit->rx, args->fit->ry};
	Homography H = { 0 };
	int nobj = AT_MATCH_CATALOG_NBRIGHT;
	int attempt = 1;

	com.stars = peaker(args->fit, 0, &com.starfinder_conf, &n_fit, NULL, TRUE); // TODO: use good layer
	if (n_fit < AT_MATCH_MINPAIRS) {
		siril_log_message(_("Not enough stars.\n"));
		return GINT_TO_POINTER(1);
	}

	cstars = malloc((MAX_STARS + 1) * sizeof(fitted_PSF *));
	if (cstars == NULL) {
		printf("Memory allocation failed: peaker\n");
		return GINT_TO_POINTER(1);
	}

	/* open the file */
	catalog = g_fopen(args->catalogStars, "r");
	if (catalog == NULL) {
		fprintf(stderr, "match_catalog: error opening file: %s\n", args->catalogStars);
		free(cstars);
		return GINT_TO_POINTER(1);
	}

	n_cat = read_catalog(catalog, cstars, image_size.y, NOMAD);

	/* make sure that arrays are not too small
	 * make  sure that the max of stars is BRIGHTEST_STARS */
	n = n_fit < n_cat ? n_fit : n_cat;
	n = n > BRIGHTEST_STARS ? BRIGHTEST_STARS : n;

	args->ret = 1;
	while (args->ret && attempt < NB_OF_MATCHING_TRY){
		args->ret = new_star_match(com.stars, cstars, n, nobj, args->s - 0.2, args->s + 0.2, &H,
				image_size);
		nobj += 50;
		attempt++;
	}

	if (!args->ret) {
		/* we only want to compare with linear function
		 * Maybe one day we will apply match with homography matrix
		 */
		TRANS trans = H_to_linear_TRANS(H);

		is_result.px_size = image_size;
		is_result.px_center.x = image_size.x / 2.0;
		is_result.px_center.y = image_size.y / 2.0;
		is_result.pixel_size = args->pixel_size;

		apply_match(trans, &is_result);

		print_platesolving_results(H, is_result);

	} else {
		siril_log_color_message(_("Plate Solving failed. The image could not be aligned with the reference stars after %d attempts.\n"), "red", attempt);
		siril_log_color_message(_("This is usually because the initial parameters (pixel size, focal length, initial coordinates) "
				"are too far from the real metadata of the image.\n"), "red");
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

	args->catalogStars = download_catalog();
	args->pixel_size = get_pixel();
	args->s = get_resolution(get_focal(), get_pixel());
	args->fit = &gfit;

	start_in_new_thread(match_catalog, args);
}

void invalidate_WCS_keywords() {
	gfit.crpix1 = 0.0;
	gfit.crpix2 = 0.0;
	gfit.crval1 = 0.0;
	gfit.crval2 = 0.0;
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

	if (gtk_tree_model_get_iter_first(treeModel, &iter) == FALSE)
		return;	//The tree is empty
	if (gtk_tree_selection_get_selected(selection, &treeModel, &iter)) { //get selected item
		GtkTreePath *path = gtk_tree_model_get_path(treeModel, &iter);
		gint *index = gtk_tree_path_get_indices(path);
		if (index) {
			selectedItem = index[0];
			update_coordinates();
		}
		gtk_tree_path_free(path);
	}
}

void on_GtkButton_IPS_metadata_clicked(GtkButton *button, gpointer user_data) {
	if (!has_any_keywords()) {
		char *msg = siril_log_message(_("No keywords found in the header.\n"));
		show_dialog(msg, _("Warning"), "dialog-warning-symbolic");
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
