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

#include "core/siril.h"
#include "core/proto.h"
#include "gui/callbacks.h"
#include "algos/PSF.h"
#include "gui/PSF_list.h"
#include "algos/star_finder.h"

static GtkListStore *liststore_stars = NULL;

enum {
	COLUMN_CHANNEL,		// int
	COLUMN_B,			// gchar[]
	COLUMN_A,			// gchar[]
	COLUMN_X0,			// gchar[]
	COLUMN_Y0,			// gchar[]
	COLUMN_FWHMX,		// gchar[]
	COLUMN_FWHMY,		// gchar[]
	COLUMN_ROUNDNESS,	// gchar[]
	COLUMN_ANGLE,		// gchar[]
	COLUMN_RMSE,		// gchar[]
	N_COLUMNS
};

enum {					//different context_id of the GtkStatusBar
	COUNT_STATE
};

void get_stars_list_store() {
	if (liststore_stars == NULL)
		liststore_stars = GTK_LIST_STORE(gtk_builder_get_object(builder, "liststore_stars"));
}

void add_star_to_list(fitted_PSF *star) {
	static GtkTreeSelection *selection = NULL;
	GtkTreeIter iter;
	gchar B[16], A[16], xpos[16], ypos[16], angle[16],
		fwhmx[16], fwhmy[16], roundness[16], rmse[16];

	get_stars_list_store();
	if (!selection)
		selection = GTK_TREE_SELECTION(gtk_builder_get_object(builder, "treeview-selection"));
	if (star == NULL) {
		gtk_list_store_clear(liststore_stars);
		return;		// just clear the list
	}
	
	g_snprintf(B, sizeof(B), "%.6lf", star->B);
	g_snprintf(A, sizeof(A), "%.6lf", star->A);
	g_snprintf(xpos, sizeof(xpos), "%.2lf", star->xpos);
	g_snprintf(ypos, sizeof(ypos), "%.2lf", star->ypos);
	g_snprintf(fwhmx, sizeof(fwhmx), "%.2lf%s", star->fwhmx, star->units);
	g_snprintf(fwhmy, sizeof(fwhmy), "%.2lf%s", star->fwhmy, star->units);
	g_snprintf(roundness, sizeof(roundness), "%8.3lf", star->fwhmy/star->fwhmx);		//ALWAYS FWHMX > FWHMY
	g_snprintf(angle, sizeof(angle), "%.2lf", star->angle);
	g_snprintf(rmse, sizeof(rmse), "%.3e", star->rmse);
	gtk_list_store_append (liststore_stars, &iter);
	gtk_list_store_set (liststore_stars, &iter,
			COLUMN_CHANNEL, star->layer,
			COLUMN_B, B,
			COLUMN_A, A,
			COLUMN_X0, xpos,
			COLUMN_Y0, ypos,
			COLUMN_FWHMX, fwhmx,
			COLUMN_FWHMY, fwhmy,
			COLUMN_ROUNDNESS, roundness,
			COLUMN_ANGLE, angle,
			COLUMN_RMSE, rmse,
			-1);
}

void fill_stars_list(fitted_PSF **star) {
	int i=0;
	if (star == NULL) return;
	add_star_to_list(NULL);	// clear  
	while (star[i]) {
		add_star_to_list(star[i]);
		i++;
	}
	com.selected_star = -1;
	display_status();	//no stars selected
}

void refresh_stars_list(fitted_PSF **star){
	get_stars_list_store();
	gtk_list_store_clear(liststore_stars);
	fill_stars_list(com.stars);
	redraw(com.cvport, REMAP_NONE);
}

void clear_stars_list() {
	get_stars_list_store();
	gtk_list_store_clear(liststore_stars);
	if (com.stars) {
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

void display_PSF(fitted_PSF **result){
	
	if (result && result[0]) {
		char msg[512];
		int i = 0;
		double FWHMx=0., FWHMy=0., B=0., A=0., r=0., angle=0., rmse=0.;

		while (result[i]) {
			B+=result[i]->B;
			A+=result[i]->A;
			FWHMx+=result[i]->fwhmx;
			FWHMy+=result[i]->fwhmy;
			angle+=result[i]->angle;
			rmse+=result[i]->rmse;
			if (i>1){
				if (strcmp(result[i]->units,result[i-1]->units)){
					show_dialog(_("Stars must have the same units."), _("Error"), "gtk-dialog-error");
					return;
				}
			}
			i++;
		}
		B = B / (double)i;
		A = A / (double)i;
		FWHMx = FWHMx / (double)i;
		FWHMy = FWHMy / (double)i;
		r = FWHMy / FWHMx;
		angle = angle / (double)i;
		rmse = rmse / (double)i;
		g_snprintf(msg, sizeof(msg), _("Average Gaussian PSF\n\n"
				"N:\t%d stars\nB:\t%.6f\nA:\t%.6f\nFWHMx:\t%.2f%s\n"
				"FWHMy:\t%.2f%s\nr:\t%.3f\nAngle:\t%.2f deg\nrmse:\t%.3e\n"),
				i, B, A, FWHMx, result[0]->units, FWHMy, result[0]->units, r, angle, rmse);
		show_data_dialog(msg, "Average Star Data");
	}
}

/* This function returns the index of a selected row.
 * The index corresponds to the index of the com.stars.
 * If no rows are selected, the function returns the value -1.
 * If "remove_row" is true, then the selected row will
 * be removed from the list */
int get_index_of_selected_line(gboolean remove_row){
	GtkTreeSelection *selection = GTK_TREE_SELECTION(gtk_builder_get_object(builder, "treeview-selection"));
	GtkTreeModel *tree_stars = gtk_tree_view_get_model(GTK_TREE_VIEW(gtk_builder_get_object(builder, "Stars_stored")));
	GtkTreeIter iter;
	int retvalue = -1;
	
	if (com.stars){
		if (gtk_tree_model_get_iter_first(tree_stars, &iter) == FALSE) return retvalue;	//The tree is empty
		if (gtk_tree_selection_get_selected(selection, &tree_stars, &iter)){	//get selected item	
			GtkTreePath *path = gtk_tree_model_get_path(tree_stars, &iter);
			int *index = gtk_tree_path_get_indices(path);
			if (index) retvalue = index[0];
			gtk_tree_path_free(path);
		}
	}
	if (remove_row == TRUE && retvalue != -1) {
		gtk_list_store_remove(GTK_LIST_STORE(tree_stars), &iter);
		gtk_tree_selection_unselect_all(selection);
	}
	return retvalue;
}

void move_selected_line() {
	int index = get_index_of_selected_line(FALSE);
	
	com.selected_star = index;
	display_status();
	redraw(com.cvport, REMAP_NONE);
}

void remove_selected_line(){
	int index = get_index_of_selected_line(TRUE);
	
	remove_star(index);
	com.selected_star = -1;
	display_status();
}

void remove_all_lines(){
	clear_stars_list();
	com.selected_star = -1;
	display_status();
	redraw(com.cvport, REMAP_NONE);
}

void display_status(){
	char text[256];
	int nbstars=0;
	GtkStatusbar *statusbar = GTK_STATUSBAR(gtk_builder_get_object(builder, "statusbar_PSF"));
	
	while (com.stars && com.stars[nbstars]) nbstars++;
	if (com.selected_star == -1) g_snprintf(text, sizeof(text), " ");
	else g_snprintf(text, sizeof(text), _("Star %d of %d"), com.selected_star + 1, nbstars);
	gtk_statusbar_push(statusbar, COUNT_STATE, text);
}
