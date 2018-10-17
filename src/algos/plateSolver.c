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
#include <curl/curl.h>

#include "core/siril.h"
#include "core/proto.h"
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

#define BRIGHTEST_STARS 2500

#define RADCONV ((3600.0 * 180.0) / M_PI) / 1.0E3

#undef DEBUG           /* get some of diagnostic output */

enum {
	APPAS, NOMAD
};

static int get_mag_limit();
static double get_focal();
static double get_pixel();
static double get_resolution(double focal, double pixel);
static double get_fov();

static struct object platedObject[RESOLVER_NUMBER];
static GtkListStore *list_IPS = NULL;
static int selectedItem = -1;
static image_solved is_result;

static int parse_curl_buffer(char *buffer, struct object *obj) {
	char **token, **fields;
	point center;
	gchar *object = NULL, *name = NULL;
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
				double hour = center.x * 24.0 / 360.0;
				platedObject[resolver].RA.hour = (int) hour;
				double min = (hour - (int) hour) * 60.0;
				platedObject[resolver].RA.min = (int) min;
				double sec = (min - (int) min) * 60.0;
				platedObject[resolver].RA.sec = sec;

				if (center.y < 0.0)
					platedObject[resolver].south = TRUE;
				else
					platedObject[resolver].south = FALSE;

				double degree = fabs(center.y);
				platedObject[resolver].Dec.degree = (int) degree;
				double degmin = ((double) (degree - (int) degree) * 60.0);
				platedObject[resolver].Dec.min = (int) degmin;
				double degsec = (degmin - (int) degmin) * 60.0;
				platedObject[resolver].Dec.sec = degsec;

				platedObject[resolver].imageCenter = center;
			}
			g_strfreev(fields);
		}
		else if (g_str_has_prefix (token[i], "%I.0 ")) {
			if (resolver != -1) {
				name = g_strstr_len(token[i], strlen(token[i]), "%I.0 ");
				gchar *realname;
				realname = strdup(name + 5);
				platedObject[resolver].name = realname;
				name = NULL;
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

static gchar *get_catalog_url(point center, int mag_limit, double dfov, int type) {
	GString *url;
	gchar coordinates[256];
	gchar mag[3];
	gchar fov[10];

	g_snprintf(coordinates, sizeof(coordinates), "%lf+%lf", center.x, center.y);
	g_snprintf(mag, sizeof(mag), "%d", mag_limit);
	g_snprintf(fov, sizeof(fov), "%d", (int)dfov / 2);

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
//		url = g_string_append(url, "&Bmag=<");
//		url = g_string_append(url, mag);
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
	gchar *buffer;
	FILE *catalog = NULL;
	FILE *fproj = NULL;
	gchar *filename, *foutput;
	int index = selectedItem;

	if (index < 0)
		return NULL;

	double fov = get_fov(get_resolution(get_focal(), get_pixel()),
			gfit.ry > gfit.rx ? gfit.ry : gfit.rx);

	url = get_catalog_url(platedObject[index].imageCenter, get_mag_limit(), fov,
			NOMAD);
	buffer = fetch_url(url);

	filename = g_build_filename(g_get_tmp_dir(), "catalog.dat", NULL);
	catalog = g_fopen(filename, "w+t");
	if (catalog == NULL) {
		fprintf(stderr, "plateSolver: Cannot open catalog\n");
		g_free(foutput);
		return NULL;
	}
	fprintf(catalog, buffer);

	foutput = g_build_filename(g_get_tmp_dir(), "catalog.proj", NULL);
	fproj = g_fopen(foutput, "w+t");
	if (fproj == NULL) {
		fprintf(stderr, "plateSolver: Cannot open fproj\n");
		g_free(foutput);
		return NULL;
	}

	convert_catalog_coords(filename, platedObject[index].imageCenter, fproj);

	is_result.px_cat_center = platedObject[index].imageCenter;

	fclose(catalog);
	g_free(filename);

	/* free data */
	g_free(url);
	return foutput;
}

void http_cleanup() {
	curl_easy_cleanup(curl);
	curl_global_cleanup();
}

/*********
 *
 */

enum {
	COLUMN_RESOLVER,		// string
	COLUMN_NAME,		// string
	N_COLUMNS
};

// RA: 03 47 13.271  Dec: +24 19 17.77
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

static void get_list_IPS() {
	if (list_IPS == NULL)
		list_IPS = GTK_LIST_STORE(gtk_builder_get_object(builder, "liststoreIPS"));
}

static void add_object_to_list() {
	static GtkTreeSelection *selection = NULL;
	GtkTreeIter iter;

	get_list_IPS();
	if (!selection)
		selection = GTK_TREE_SELECTION(gtk_builder_get_object(builder, "GtkTreeSelectionIPS"));

	gtk_list_store_clear(list_IPS);

if (platedObject[RESOLVER_NED].name) {
		gtk_list_store_append(list_IPS, &iter);
		gtk_list_store_set(list_IPS, &iter, COLUMN_RESOLVER, "NEV", COLUMN_NAME,
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
			platedObject[index].Dec.degree);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_Dec_m")),
			platedObject[index].Dec.min);
	gtk_entry_set_text(GTK_ENTRY(lookup_widget("GtkEntryIPS_Dec_s")), Dec_sec);
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
	return (RADCONV / focal) * pixel;
}

/* get FOV in arcsec/px */
static double get_fov(double resolution, int image_size) {
	return (resolution * (double)image_size) / 60.0;
}


static int get_mag_limit() {
	GtkSpinButton *magButton = GTK_SPIN_BUTTON(lookup_widget("GtkSpinIPS_Mag_Limit"));

	return gtk_spin_button_get_value_as_int(magButton);
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

void update_IPS_GUI() {
	update_pixel_size();
	update_focal();
	update_resolution_field();
}

void on_GtkEntry_IPS_changed(GtkEditable *editable, gpointer user_data) {
	update_resolution_field();
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
	double rotation, resolution, scaleX, scaleY, focal, pixel;
	double inliers;
	char RA[256] = { 0 };
	char DEC[256] = { 0 };

	/* Matching information */
	siril_log_message(_("%d pair matches.\n"), H.pair_matched);
	inliers = 1.0 - ((((double) H.pair_matched - (double) H.Inliers)) / (double) H.pair_matched);
	siril_log_message(_("Inliers:%*.3f\n"), 14, inliers);

	/* Plate Solving */
	scaleX = sqrt(H.h00 * H.h00 + H.h01 * H.h01);
	scaleY = sqrt(H.h10 * H.h10 + H.h11 * H.h11);
	resolution = (scaleX + scaleY) * 0.5;
	siril_log_message(_("Resolution:%*.3f arcsec/px\n"), 11, resolution);
	rotation = atan2(H.h01, H.h00);
	siril_log_message(_("Rotation:%+*.2f deg\n"), 12, rotation * 180 / M_PI);
	pixel = get_pixel();
	focal = RADCONV * pixel / resolution;

	image.focal = focal;
	image.pixel_size = pixel;

	siril_log_message(_("Focal:%*.2f mm\n"), 15, focal);
	siril_log_message(_("Pixel size:%*.2f µm\n"), 10, pixel);
	deg_to_HMS(image.ra_center, "ra", RA);
	deg_to_HMS(image.dec_center, "dec", DEC);
	siril_log_message(_("Image center: RA: %s, DEC: %s\n"), RA, DEC);

	update_gfit(image);
	update_IPS_GUI();
}

static int read_NOMAD_catalog(FILE *catalog, fitted_PSF **cstars, int shift_y) {
	char line[LINELEN];
	fitted_PSF *psf;
	double r, x, y, Vmag, Bmag;
	int i = 0;

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
		int n = sscanf(line, "%lf %lf %lf %lf %lf", &r, &x, &y, &Vmag, &Bmag);
		if (n < 5)
			continue;
		psf = malloc(sizeof(fitted_PSF));
		psf->xpos = x;
		psf->ypos = -y + shift_y;
		psf->mag = Vmag;
		cstars[i] = psf;
		cstars[i + 1] = NULL;
		i++;
	}
	sort_stars(cstars, i);

	return i;
}

static int read_APPAS_catalog(FILE *catalog, fitted_PSF **cstars, int shift_y) {
	char line[LINELEN];
	fitted_PSF *psf;
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
		psf = malloc(sizeof(fitted_PSF));
		psf->xpos = x;
		psf->ypos = -y + shift_y;
		psf->mag = magA;
		cstars[i] = psf;
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

TRANS H_to_linear_TRANS(Homography H) {
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

static int match_catalog(gchar *catalogStars) {
	FILE *catalog;
	fitted_PSF **cstars;
	int n_fit, n_cat, n;
	point image_size = {gfit.rx, gfit.ry};
	Homography H = { 0 };

	com.stars = peaker(&gfit, 0, &com.starfinder_conf, &n_fit, NULL, TRUE); // TODO: use good layer
	if (n_fit < AT_MATCH_MINPAIRS) {
		siril_log_message(_("Not enough stars.\n"));
		return -1;
	}

	cstars = malloc((MAX_STARS + 1) * sizeof(fitted_PSF *));
	if (cstars == NULL) {
		printf("Memory allocation failed: peaker\n");
		return 1;
	}

	/* open the file */
	catalog = g_fopen(catalogStars, "r");
	if (catalog == NULL) {
		fprintf(stderr, "match_catalog: error opening file: %s\n", catalogStars);
		free(cstars);
		return 1;
	}

	n_cat = read_catalog(catalog, cstars, image_size.y, NOMAD);

	/* make sure that arrays are not too small
	 * make  sure that the max of stars is BRIGHTEST_STARS */
	n = n_fit < n_cat ? n_fit : n_cat;
	n = n > BRIGHTEST_STARS ? BRIGHTEST_STARS : n;

	int ret = new_star_match(com.stars, cstars, n, AT_MATCH_CATALOG_NBRIGHT, &H, image_size);

	if (!ret) {
		TRANS trans = H_to_linear_TRANS(H);

		is_result.px_center.x = image_size.x / 2.0;
		is_result.px_center.y = image_size.y / 2.0;
		apply_match(trans, &is_result);

		print_platesolving_results(H, is_result);
		update_IPS_GUI();
	} else {
		siril_log_color_message(_("Plate Solving failed. The image could not be aligned with the reference star field.\n"), "red");
		siril_log_color_message(_("This is usually because the initial parameters (pixel size, focal length, initial coordinates) "
				"are too far from the real metadata of the image.\n"), "red");
	}

	/* free data */
	if (n_cat > 0) free_fitted_stars(cstars);
	clear_stars_list();
	return 0;
}

static void search_object_in_catalogs(const gchar *object) {
	GString *url;
	gchar *gcurl, *result;
	struct object obj;

	set_cursor_waiting(TRUE);
	url = g_string_new(CDSSESAME);
	url = g_string_append(url, "/-oI/A?");
	url = g_string_append(url, object);
	gcurl = g_string_free(url, FALSE);
	result = fetch_url(gcurl);
	parse_curl_buffer(result, &obj);
	add_object_to_list();
	free_Platedobject();
	set_cursor_waiting(FALSE);
}

static void start_image_plate_solve() {
	gchar *catalog;

	set_cursor_waiting(TRUE);
	catalog = download_catalog();
	match_catalog(catalog);
	set_cursor_waiting(FALSE);

	g_free(catalog);
}

/*****
 * CALLBACKS FUNCTIONS
 */

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

	if (gtk_tree_model_get_iter_first(treeModel, &iter) == FALSE) return;	//The tree is empty
	if (gtk_tree_selection_get_selected(selection, &treeModel, &iter)){	//get selected item
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
	update_IPS_GUI();
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
