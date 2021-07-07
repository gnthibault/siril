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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/siril_world_cs.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "gui/PSF_list.h"
#include "gui/progress_and_log.h"
#include "algos/siril_wcs.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"

static GtkListStore *liststore_stars = NULL;

enum {
	COLUMN_CHANNEL,		// int
	COLUMN_B,			// gdouble
	COLUMN_A,			// gdouble
	COLUMN_X0,			// gdouble
	COLUMN_Y0,			// gdouble
	COLUMN_FWHMX,		// gdouble
	COLUMN_FWHMY,		// gdouble
	COLUMN_MAG,		    // gdouble
	COLUMN_ROUNDNESS,	// gdouble
	COLUMN_ANGLE,		// gdouble
	COLUMN_RMSE,		// gdouble
	N_COLUMNS
};

enum {					//different context_id of the GtkStatusBar
	COUNT_STATE
};

static gchar *units = "";

static void gdouble_fwhmx_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_FWHMX, &var, -1);
	buf = g_strdup_printf("%.2f%s", var, units);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void gdouble_fwhmy_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_FWHMY, &var, -1);
	buf = g_strdup_printf("%.2f%s", var, units);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void gdouble_x0_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_X0, &var, -1);
	buf = g_strdup_printf("%.2f", var);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void gdouble_y0_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_Y0, &var, -1);
	buf = g_strdup_printf("%.2f", var);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void gdouble_mag_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_MAG, &var, -1);
	buf = g_strdup_printf("%.2f", var);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void gdouble_r_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_ROUNDNESS, &var, -1);
	buf = g_strdup_printf("%.3f", var);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void gdouble_angle_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_ANGLE, &var, -1);
	if (var == 0.0)
		buf = g_strdup_printf("%s", "N/A");
	else
		buf = g_strdup_printf("%.2f", var);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void gdouble_rmse_cell_data_function(GtkTreeViewColumn *col,
		GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter,
		gpointer user_data) {
	gdouble var;
	gchar *buf;
	gtk_tree_model_get(model, iter, COLUMN_RMSE, &var, -1);
	buf = g_strdup_printf("%.2e", var);
	g_object_set(renderer, "text", buf, NULL);

	g_free(buf);
}

static void get_stars_list_store() {
	if (liststore_stars == NULL)
		liststore_stars = GTK_LIST_STORE(gtk_builder_get_object(builder, "liststore_stars"));

	GtkTreeViewColumn *col;
	GtkCellRenderer *cell;

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(builder, "treeviewcolumn7"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(builder, "cell_x0"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_x0_cell_data_function, NULL, NULL);

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(builder, "treeviewcolumn8"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(builder, "cell_y0"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_y0_cell_data_function, NULL, NULL);

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(builder, "treeviewcolumn9"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(builder, "cell_fwhmx"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_fwhmx_cell_data_function, NULL, NULL);

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(builder, "treeviewcolumn10"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(builder, "cell_fwhmy"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_fwhmy_cell_data_function, NULL, NULL);

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(builder, "treeviewcolumn_mag"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(builder, "cell_mag"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_mag_cell_data_function, NULL, NULL);

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(builder, "treeviewcolumn14"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(builder, "cell_r"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_r_cell_data_function, NULL, NULL);

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(builder, "treeviewcolumn6"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(builder, "cell_angle"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_angle_cell_data_function, NULL, NULL);

	col = GTK_TREE_VIEW_COLUMN(gtk_builder_get_object(builder, "treeviewcolumn15"));
	cell = GTK_CELL_RENDERER(gtk_builder_get_object(builder, "cell_rmse"));
	gtk_tree_view_column_set_cell_data_func(col, cell, gdouble_rmse_cell_data_function, NULL, NULL);
}

static void display_PSF(fitted_PSF **result) {
	if (result) {
		gchar *msg;
		int i = 0;
		double FWHMx = 0.0, FWHMy = 0.0, B = 0.0, A = 0.0, r = 0.0, angle = 0.0,
				rmse = 0.0;
		gboolean unit_is_arcsec;

		while (result[i]) {
			double fwhmx, fwhmy;
			char *unit;
			gboolean is_as = get_fwhm_as_arcsec_if_possible(result[i], &fwhmx, &fwhmy, &unit);
			if (i == 0)
				unit_is_arcsec = is_as;
			else if (is_as != unit_is_arcsec) {
				siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"),
						_("Stars FWHM must have the same units."));
				return;
			}

			B += result[i]->B;
			A += result[i]->A;
			FWHMx += fwhmx;
			FWHMy += fwhmy;
			angle += result[i]->angle;
			rmse += result[i]->rmse;
			i++;
		}
		if (i <= 0) return;
		/* compute average */
		B = B / (double)i;
		A = A / (double)i;
		FWHMx = FWHMx / (double)i;
		FWHMy = FWHMy / (double)i;
		r = FWHMy / FWHMx;
		angle = angle / (double)i;
		rmse = rmse / (double)i;

		msg = g_strdup_printf(_("Average Gaussian PSF\n\n"
				"N:\t%d stars\nB:\t%.6f\nA:\t%.6f\nFWHMx:\t%.2f%s\n"
				"FWHMy:\t%.2f%s\nr:\t%.3f\nAngle:\t%.2f deg\nrmse:\t%.3e\n"),
				i, B, A, FWHMx, result[0]->units, FWHMy,
				result[0]->units, r, angle, rmse);
		show_data_dialog(msg, _("Average Star Data"), "stars_list_window", NULL);
		g_free(msg);
	}
}

static gint get_index_of_selected_star(gdouble x, gdouble y) {
	int i = 0;

	while (com.stars && com.stars[i]) {
		if ((com.stars[i]->xpos == x) && (com.stars[i]->ypos == y)) {
			return i;
		}
		i++;
	}
	return -1;
}

static void display_status() {
	gchar *text;
	int i = 0;
	GtkStatusbar *statusbar;

	statusbar = GTK_STATUSBAR(lookup_widget("statusbar_PSF"));

	while (com.stars && com.stars[i])
		i++;
	if (com.selected_star == -1) {
		if (i > 0) {
			text = ngettext("%d star", "%d stars", i);
			text = g_strdup_printf(text, i);
		} else {
			text = g_strdup(" ");
		}
	} else {
		text = g_strdup_printf(_("Star %d of %d"), com.selected_star + 1, i);
	}
	gtk_statusbar_push(statusbar, COUNT_STATE, text);
	g_free(text);
}

static void remove_selected_star(int index) {
	GtkTreeSelection *selection = GTK_TREE_SELECTION(gtk_builder_get_object(builder, "treeview-selection"));
	GtkTreeModel *model = gtk_tree_view_get_model(GTK_TREE_VIEW(gtk_builder_get_object(builder, "Stars_stored")));
	GtkTreeIter iter;

	if (gtk_tree_selection_get_selected(selection, &model, &iter)) {
		gtk_list_store_remove(GTK_LIST_STORE(model), &iter);
		gtk_tree_selection_unselect_all(selection);

		remove_star(index);

		com.selected_star = -1;
		display_status();
	}
}

static void remove_all_stars(){
	clear_stars_list();
	com.selected_star = -1;
	display_status();
	redraw(com.cvport, REMAP_NONE);
}

static int save_list(gchar *filename) {
	int i = 0;
	if (!com.stars)
		return 1;
	GError *error = NULL;
	gboolean is_in_arcsec = FALSE;

	GFile *file = g_file_new_for_path(filename);
	GOutputStream *output_stream = (GOutputStream*) g_file_replace(file, NULL, FALSE,
			G_FILE_CREATE_NONE, NULL, &error);

	if (output_stream == NULL) {
		if (error != NULL) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
			fprintf(stderr, "save_list: Cannot save star list\n");
		}
		g_object_unref(file);
		return 1;
	}

	gchar *buffer = g_strdup_printf("star#\tlayer\tB\tA\tX\tY\tFWHMx [%s]\tFWHMy [%s]\tangle\tRMSE\tmag%s", com.stars[0]->units,com.stars[0]->units,SIRIL_EOL);
	if (!g_output_stream_write_all(output_stream, buffer, strlen(buffer), NULL, NULL, &error)) {
		g_warning("%s\n", error->message);
		g_free(buffer);
		g_clear_error(&error);
		g_object_unref(output_stream);
		g_object_unref(file);
		return 1;
	}
	g_free(buffer);
	if (com.stars[0]) {
		is_in_arcsec = com.stars[0]->fwhmx_arcsec > 0;
	}
	while (com.stars[i]) {
		if (is_in_arcsec) { 
			buffer = g_strdup_printf(
					"%d\t%d\t%10.6f\t%10.6f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%3.2f\t%10.3e\t%10.2f%s",
					i + 1, com.stars[i]->layer, com.stars[i]->B, com.stars[i]->A,
					com.stars[i]->xpos, com.stars[i]->ypos, com.stars[i]->fwhmx_arcsec,
					com.stars[i]->fwhmy_arcsec, com.stars[i]->angle, com.stars[i]->rmse, com.stars[i]->mag, SIRIL_EOL);
		} else {
			buffer = g_strdup_printf(
					"%d\t%d\t%10.6f\t%10.6f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%3.2f\t%10.3e\t%10.2f%s",
					i + 1, com.stars[i]->layer, com.stars[i]->B, com.stars[i]->A,
					com.stars[i]->xpos, com.stars[i]->ypos, com.stars[i]->fwhmx,
					com.stars[i]->fwhmy, com.stars[i]->angle, com.stars[i]->rmse, com.stars[i]->mag, SIRIL_EOL);
		}

		if (!g_output_stream_write_all(output_stream, buffer, strlen(buffer), NULL, NULL, &error)) {
			g_warning("%s\n", error->message);
			g_free(buffer);
			g_clear_error(&error);
			g_object_unref(output_stream);
			g_object_unref(file);
			return 1;
		}
		i++;
		g_free(buffer);
	}
	siril_log_message(_("The file %s has been created.\n"), filename);
	g_object_unref(file);

	return 0;
}

static void set_filter(GtkFileChooser *dialog) {
	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, _("Star list file (*.lst)"));
	gtk_file_filter_add_pattern(f, "*.lst");
	gtk_file_chooser_add_filter(dialog, f);
	gtk_file_chooser_set_filter(dialog, f);
}

static void save_stars_dialog() {
	SirilWidget *widgetdialog;
	GtkFileChooser *dialog = NULL;
	GtkWindow *parent = GTK_WINDOW(lookup_widget("stars_list_window"));
	gint res;

	widgetdialog = siril_file_chooser_save(parent, GTK_FILE_CHOOSER_ACTION_SAVE);
	dialog = GTK_FILE_CHOOSER(widgetdialog);
	gtk_file_chooser_set_current_folder(dialog, com.wd);
	gtk_file_chooser_set_select_multiple(dialog, FALSE);
	gtk_file_chooser_set_do_overwrite_confirmation(dialog, TRUE);
	gtk_file_chooser_set_current_name(dialog, "stars.lst");
	set_filter(dialog);

	res = siril_dialog_run(widgetdialog);
	if (res == GTK_RESPONSE_ACCEPT) {
		gchar *file = gtk_file_chooser_get_filename(dialog);
		save_list(file);

		g_free(file);
	}
	siril_widget_destroy(widgetdialog);
}

/********************* public ***********************/

void add_star_to_list(fitted_PSF *star) {
	static GtkTreeSelection *selection = NULL;
	GtkTreeIter iter;

	get_stars_list_store();
	if (!selection)
		selection = GTK_TREE_SELECTION(gtk_builder_get_object(builder, "treeview-selection"));
	if (star == NULL) {
		gtk_list_store_clear(liststore_stars);
		return;		// just clear the list
	}

	double fwhmx = star->fwhmx_arcsec < 0 ? star->fwhmx : star->fwhmx_arcsec;
	double fwhmy = star->fwhmy_arcsec < 0 ? star->fwhmy : star->fwhmy_arcsec;

	gtk_list_store_append (liststore_stars, &iter);
	gtk_list_store_set (liststore_stars, &iter,
			COLUMN_CHANNEL, star->layer,
			COLUMN_B, star->B,
			COLUMN_A, star->A,
			COLUMN_X0, star->xpos,
			COLUMN_Y0, star->ypos,
			COLUMN_FWHMX, fwhmx,
			COLUMN_FWHMY, fwhmy,
			COLUMN_MAG, star->mag,
			COLUMN_ROUNDNESS, fwhmy / fwhmx,
			COLUMN_ANGLE, star->angle,
			COLUMN_RMSE, star->rmse,
			-1);

	units = star->units;
	display_status();
}

void fill_stars_list(fits *fit, fitted_PSF **stars) {
	int i = 0;
	if (stars == NULL)
		return;
	add_star_to_list(NULL);	// clear

	while (stars[i]) {
		/* update units if needed */
		fwhm_to_arcsec_if_needed(fit, stars[i]);
		add_star_to_list(stars[i]);
		i++;
	}
	com.selected_star = -1;
}

void refresh_star_list(fitted_PSF **star){
	get_stars_list_store();
	gtk_list_store_clear(liststore_stars);
	fill_stars_list(&gfit, com.stars);
	redraw(com.cvport, REMAP_NONE);
}

void clear_stars_list() {
	if (com.stars) {
		if (!com.headless) {
			get_stars_list_store();
			gtk_list_store_clear(liststore_stars);
		}
		if (com.stars[0]) {
			/* freeing found stars. It must not be done when the only star in
			 * com.stars is the same as com.seq.imgparam[xxx].fwhm, as set in
			 * set_fwhm_star_as_star_list(), because it will be reused */
			if (com.stars[1] || !com.star_is_seqdata) {
				int i = 0;
				while (i < MAX_STARS && com.stars[i])
					free(com.stars[i++]);
			}
		}
		free(com.stars);
		com.stars = NULL;
	}
	com.star_is_seqdata = FALSE;
}

void pick_a_star() {
	int layer = match_drawing_area_widget(com.vport[com.cvport], FALSE);
	int new_index;

	if (layer != -1) {
		if (!(com.selection.h && com.selection.w))
			return;
		if (com.selection.w > 300 || com.selection.h > 300) {
			siril_message_dialog(GTK_MESSAGE_WARNING, _("Current selection is too large"),
					_("To determine the PSF, please make a selection around a star."));
			return;
		}
		fitted_PSF *new_star = add_star(&gfit, layer, &new_index);
		if (new_star) {
			add_star_to_list(new_star);
			siril_open_dialog("stars_list_window");
		} else
			return;
	}
	redraw(com.cvport, REMAP_NONE);
}

static gchar *build_wcs_url(gchar *ra, gchar *dec) {
	if (!has_wcs(&gfit)) return NULL;

	double resolution = get_wcs_image_resolution(&gfit);

	gchar *tol = g_strdup_printf("%lf", resolution * 3600 * 15);

	GString *url = g_string_new("https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=");
	url = g_string_append(url, ra);
	url = g_string_append(url, dec);
	url = g_string_append(url, "&Radius=");
	url = g_string_append(url, tol);
	url = g_string_append(url, "&Radius.unit=arcsec");
	url = g_string_append(url, "#lab_basic");

	gchar *simbad_url = g_string_free(url, FALSE);
	gchar *cleaned_url = url_cleanup(simbad_url);

	g_free(tol);
	g_free(simbad_url);

	return cleaned_url;
}

static const char *SNR_quality(double SNR) {
	if (SNR <= 0.0) return _("N/A");
	if (SNR <= 10.0) return _("Bad");
	if (SNR <= 15.0) return _("Poor");
	if (SNR <= 25.0) return _("Fair");
	if (SNR <= 40.0) return _("Good");
	else return _("Excellent");
}

void popup_psf_result(fitted_PSF *result) {
	gchar *msg, *coordinates, *url = NULL;
	const char *str;
	if (com.magOffset > 0.0)
		str = "true reduced";
	else
		str = "relative";

	double x = result->x0 + com.selection.x;
	double y = com.selection.y + com.selection.h - result->y0;
	if (has_wcs(&gfit)) {
		double world_x, world_y;
		SirilWorldCS *world_cs;

		pix2wcs(&gfit, x, (double) gfit.ry - y, &world_x, &world_y);
		world_cs = siril_world_cs_new_from_a_d(world_x, world_y);
		if (world_cs) {
			gchar *ra = siril_world_cs_alpha_format(world_cs, "%02d %02d %.3lf");
			gchar *dec = siril_world_cs_delta_format(world_cs, "%c%02d %02d %.3lf");

			url = build_wcs_url(ra, dec);

			g_free(ra);
			g_free(dec);

			ra = siril_world_cs_alpha_format(world_cs, " %02dh%02dm%02ds");
			dec = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%02d\"");

			coordinates = g_strdup_printf("x0=%.2fpx\t%s J2000\n\t\ty0=%.2fpx\t%s J2000", x, ra, y, dec);

			g_free(ra);
			g_free(dec);
			siril_world_cs_unref(world_cs);
		} else {
			coordinates = g_strdup_printf("x0=%.2fpx\n\t\ty0=%.2fpx", x, y);
		}
	} else {
		coordinates = g_strdup_printf("x0=%.2fpx\n\t\ty0=%.2fpx", x, y);
	}

	double fwhmx, fwhmy;
	char *units;
	get_fwhm_as_arcsec_if_possible(result, &fwhmx, &fwhmy, &units);
	msg = g_strdup_printf(_("Centroid Coordinates:\n\t\t%s\n\n"
				"Full Width Half Maximum:\n\t\tFWHMx=%.2f%s\n\t\tFWHMy=%.2f%s\n\n"
				"Angle:\n\t\t%0.2fdeg\n\n"
				"Background Value:\n\t\tB=%.6f\n\n"
				"Maximal Intensity:\n\t\tA=%.6f\n\n"
				"Magnitude (%s):\n\t\tm=%.4f\u00B1%.4f\n\n"
				"Signal-to-noise ratio:\n\t\tSNR=%.1fdB (%s)\n\n"
				"RMSE:\n\t\tRMSE=%.3e"),
			coordinates, fwhmx, units, fwhmy, units,
			result->angle, result->B, result->A, str,
			result->mag + com.magOffset, result->s_mag, result->SNR,
			SNR_quality(result->SNR), result->rmse);
	show_data_dialog(msg, "PSF Results", NULL, url);
	g_free(coordinates);
	g_free(msg);
	g_free(url);
}

/***************** callbacks ****************/

void on_treeview_cursor_changed(GtkTreeView *tree_view,
		gpointer user_data) {
	GtkTreeModel *treeModel = gtk_tree_view_get_model(tree_view);
	GtkTreeSelection *selection = gtk_tree_view_get_selection (tree_view);
	GtkTreeIter iter;
	GValue value_x = G_VALUE_INIT;
	GValue value_y = G_VALUE_INIT;

	if (gtk_tree_model_get_iter_first(treeModel, &iter) == FALSE)
		return;	//The tree is empty
	if (gtk_tree_selection_get_selected(selection, &treeModel, &iter)) { //get selected item
		gdouble x0, y0;

		gtk_tree_model_get_value(treeModel, &iter, COLUMN_X0, &value_x);
		x0 = g_value_get_double(&value_x);
		gtk_tree_model_get_value(treeModel, &iter, COLUMN_Y0, &value_y);
		y0 = g_value_get_double(&value_y);

		g_value_unset(&value_x);
		g_value_unset(&value_y);

		com.selected_star = get_index_of_selected_star(x0, y0);
		display_status();
		redraw(com.cvport, REMAP_NONE);
	}
}

void on_Stars_stored_key_release_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {
	if (event->keyval == GDK_KEY_Delete || event->keyval == GDK_KEY_KP_Delete
			|| event->keyval == GDK_KEY_BackSpace) {

		remove_selected_star(com.selected_star);
	}
}

void on_stars_list_window_hide(GtkWidget *object, gpointer user_data) {
	com.selected_star = -1;
}

void on_sum_button_clicked(GtkButton *button, gpointer user_data) {
	display_PSF(com.stars);
}

void on_remove_button_clicked(GtkButton *button, gpointer user_data) {
	remove_selected_star(com.selected_star);
}

void on_remove_all_button_clicked(GtkButton *button, gpointer user_data) {
	remove_all_stars();
}

void on_process_starfinder_button_clicked(GtkButton *button, gpointer user_data) {
	int nbstars = 0;
	int layer = com.cvport == RGB_VPORT ? GLAYER : com.cvport;
	if (!single_image_is_loaded() && !sequence_is_loaded()) {
		siril_log_color_message(_("Load an image first, aborted.\n"), "red");
		return;
	}
	set_cursor_waiting(TRUE);

	confirm_peaker_GUI(); //making sure the spin buttons values are read even without confirmation
	delete_selected_area();
	com.stars = peaker(&gfit, layer, &com.starfinder_conf, &nbstars, NULL, TRUE, FALSE);
	siril_log_message(_("Found %d stars in image, channel #%d\n"), nbstars, layer);
	if (com.stars)
		refresh_star_list(com.stars);
	set_cursor_waiting(FALSE);
}

void on_export_button_clicked(GtkButton *button, gpointer user_data) {
	if (com.stars) {
		save_stars_dialog();
	} else {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Nothing to export"),
				_("There are no stars in the list."));
	}
}

void on_stars_list_window_show(GtkWidget *widget, gpointer user_data) {
	update_peaker_GUI();
	fill_stars_list(&gfit, com.stars);
}

void on_button_stars_list_ok_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("stars_list_window");
}

