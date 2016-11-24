/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2016 team free-astro (see more in AUTHORS file)
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
#  include <config.h>
#endif
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "gui/callbacks.h"
#include "io/conversion.h"
#include "io/single_image.h"
#include "gui/PSF_list.h"
#include "gui/histogram.h"
#include "core/undo.h"

/* Closes and frees resources attached to the single image opened in gfit.
 * If a sequence is loaded and one of its images is displayed, nothing is done.
 */
void close_single_image() {
	if (sequence_is_loaded() && com.seq.current != RESULT_IMAGE)
		return;
	free_image_data();
	undo_flush();
}

/* frees resources when changing sequence or closing a single image
 * (image size may vary, so we reallocate) */
void free_image_data() {
	int i=0, vport;
	fprintf(stdout, "free_image_data() called, clearing loaded image\n");
	/* WARNING: single_image.fit references the actual fits image,
	 * shouldn't it be used here instead of gfit? */
	clearfits(&gfit);
	clear_stars_list();
	clear_histograms();
	delete_selected_area();
	clear_sampling_setting_box();	// clear focal and pixel pitch info

	for (vport=0; vport<MAXGRAYVPORT; vport++) {
		if (com.graybuf[vport]) {
			free(com.graybuf[vport]);
			com.graybuf[vport] = NULL;
		}
		if (com.surface[vport]) {
			cairo_surface_destroy(com.surface[vport]);
			com.surface[vport] = NULL;
		}
		com.surface_stride[vport] = 0;
		com.surface_height[vport] = 0;
	}
	activate_tab(RED_VPORT);
	if (com.rgbbuf) {
		free(com.rgbbuf);
		com.rgbbuf = NULL;
	}
	if (com.uniq) {
		if (com.uniq->filename)
			free(com.uniq->filename);
		if (com.uniq->comment)
			free(com.uniq->comment);
		free(com.uniq);
		com.uniq = NULL;
	}

	/* TODO: free alignment preview data */
	for (i=0; i<PREVIEW_NB; i++) {
		if (com.preview_surface[i]) {
			cairo_surface_destroy(com.preview_surface[i]);
			com.preview_surface[i] = NULL;
		}
	}
	if (com.refimage_surface) {
		cairo_surface_destroy(com.refimage_surface);
		com.refimage_surface = NULL;
	}
}

/* Reads an image from disk and stores it in the user allocated destination
 * fits. The realname_out argument can be NULL, and if not, it is set to the
 * real file name of the loaded file, since the given filename can be without
 * extension. This realname output has to be freed by user.
 */
int read_single_image(const char* filename, fits *dest, char **realname_out) {
	int retval;
	image_type imagetype;
	char *realname;

	realname = malloc(strlen(filename) + 10);
	retval = stat_file(filename, &imagetype, realname);
	if (retval) {
		char *msg = siril_log_message(_("Error opening image %s: file not found or not supported.\n"), filename);
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		set_cursor_waiting(FALSE);
		free(realname);
		return 1;
	}
	if (imagetype == TYPESER || imagetype == TYPEAVI) {
			/* Returns 3 if ok */
		retval = read_single_sequence(realname, imagetype);
	} else {
		retval = any_to_fits(imagetype, realname, dest);
		if (!retval)
			debayer_if_needed(imagetype, dest, com.debayer.compatibility);
	}
	if (retval != 0 && retval != 3)
		siril_log_message(_("Opening %s failed.\n"), realname);
	if (realname_out)
		*realname_out = realname;
	else free(realname);
	set_GUI_CAMERA();
	com.filter = (int)imagetype;
	return retval;
}

/* This function is used to load a single image, meaning outside a sequence,
 * whether a sequence is loaded or not, whether an image was already loaded or
 * not. The opened file is available in the usual global variable for current
 * image, gfit.
 */
int open_single_image(const char* filename) {
	int retval;
	char *realname;

	close_single_image();	// close the previous image and free resources

	retval = read_single_image(filename, &gfit, &realname);
	
	/* A single sequence has been loaded */
	if (retval == 3) {
		return 0;
	}
	if (retval == 2) {
		show_dialog(_("This file could not be opened because its extension is not supported.\n"), _("Error"), "gtk-dialog-error");
		return 1;
	}
	if (retval == 1) {
		show_dialog(_("There was an error when opening this image. See the log for more information."), _("Error"), "gtk-dialog-error");
		return 1;
	}

	if (sequence_is_loaded()) {
		char *basename = g_path_get_basename(com.seq.seqname);
		siril_log_message(_("Closing sequence %s\n"), basename);
		g_free(basename);
		clear_sequence_list();
		free_sequence(&(com.seq), FALSE);
		initialize_sequence(&com.seq, FALSE);
	}
	fprintf(stdout, "Loading image OK, now displaying\n");
	open_single_image_from_gfit(realname);
	return 0;
}

/* creates a single_image structure and displays a single image, found in gfit.
 * The argument, realname, should be the file path, and it must be allocated
 * on the heap, since destroying the single image will attempt to free it.
 */
void open_single_image_from_gfit(char *realname) {
	/* now initializing everything
	 * code based on seq_load_image or set_seq (sequence.c) */
	com.seq.current = UNRELATED_IMAGE;
	com.uniq = calloc(1, sizeof(single));
	com.uniq->filename = realname;
	com.uniq->nb_layers = gfit.naxes[2];
	com.uniq->layers = calloc(com.uniq->nb_layers, sizeof(layer_info));
	com.uniq->fit = &gfit;
	initialize_display_mode();

	image_find_minmax(&gfit, 0);
	init_layers_hi_and_lo_values(MIPSLOHI);		// If MIPS-LO/HI exist we load these values. If not it is min/max

	sliders_mode_set_state(com.sliders);
	set_cutoff_sliders_max_values();
	set_cutoff_sliders_values();

	set_display_mode();
	set_prepro_button_sensitiveness(); 			// enable or not the preprobutton
	adjust_sellabel();

	display_filename();	// display filename in gray window

	/* update menus */
	update_MenuItem();

	redraw(com.cvport, REMAP_ALL);
	update_used_memory();
	show_main_gray_window();
	adjust_vport_size_to_image();
	update_gfit_histogram_if_needed();
	if (gfit.naxes[2] == 3)
		show_rgb_window();
	else
		hide_rgb_window();
	close_tab();
}

/* searches the image for minimum and maximum pixel value, on each layer
 * the values are stored in fit->min[layer] and fit->max[layer] */
int image_find_minmax(fits *fit, int force_minmax){
	int i, layer;

	/* this should only be done once per image */
	if (fit->maxi != 0 && !force_minmax) return 1;

	/* search for min and max values in all layers */
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(layer, i) schedule(dynamic,1) if(fit->naxes[2] > 1)
#endif
	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		WORD *buf = fit->pdata[layer];
		fit->max[layer] = 0;
		fit->min[layer] = USHRT_MAX;

		for (i = 0; i < fit->rx * fit->ry; ++i) {
			fit->max[layer] = max(fit->max[layer], buf[i]);
			fit->min[layer] = min(fit->min[layer], buf[i]);
		}
		/* fprintf(stdout, "min and max (layer %d): %hu %hu\n",
				layer, fit->min[layer], fit->max[layer]); */
	}

	/* set the overall min and max values from layer values */
	fit->maxi = 0;
	fit->mini = USHRT_MAX;
	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		fit->maxi = max(fit->max[layer], fit->maxi);
		fit->mini = min(fit->min[layer], fit->mini);
	}
	return 0;
}

/* gfit has been loaded, now we copy the hi and lo values into the com.uniq or com.seq layers.
 * gfit.hi and gfit.lo may only be available in some FITS files; if they are not available, the
 * min and max value for the layer is used.
 * If gfit changed, its hi and lo values need to be updated, and they are taken from min and
 * max.
 */
void init_layers_hi_and_lo_values(sliders_mode force_minmax) {
	if (force_minmax == USER) return;
	int i, nb_layers;
	layer_info *layers=NULL;
	static GtkToggleButton *chainedbutton = NULL;
	gboolean is_chained;
	
	chainedbutton = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_chain"));
	is_chained = gtk_toggle_button_get_active(chainedbutton);

	if (com.uniq && com.uniq->layers && com.seq.current != RESULT_IMAGE) {
		nb_layers = com.uniq->nb_layers;
		layers = com.uniq->layers;
	} else if (sequence_is_loaded() && com.seq.layers) {
		nb_layers = com.seq.nb_layers;
		layers = com.seq.layers;
	} else {
		fprintf(stderr, "COULD NOT INIT HI AND LO VALUES\n");
		return;
	}
	for (i=0; i<nb_layers; i++) {
		if (gfit.hi == 0 || force_minmax == MINMAX) {
			com.sliders = MINMAX;
			if (!is_chained) {
				layers[i].hi = gfit.max[i];
				layers[i].lo = gfit.min[i];
			}
			else {
				layers[i].hi = gfit.maxi;
				layers[i].lo = gfit.mini;
			}
		} else {
			com.sliders = MIPSLOHI;
			layers[i].hi = gfit.hi;
			layers[i].lo = gfit.lo;
		}
	}
}

/* was level_adjust, to call when gfit changed and need min/max to be recomputed. */
void adjust_cutoff_from_updated_gfit() {
	image_find_minmax(&gfit, 1);
	update_gfit_histogram_if_needed();
	init_layers_hi_and_lo_values(com.sliders);
	set_cutoff_sliders_values();
}

void unique_free_preprocessing_data(single *uniq) {
	// free opened files
	if (uniq->ppprefix) {
		free(uniq->ppprefix);
		uniq->ppprefix = NULL;
	}
	if (uniq->offset) {
		clearfits(uniq->offset);
		free(uniq->offset);
		uniq->offset = NULL;
	}
	if (uniq->dark) {
		clearfits(uniq->dark);
		free(uniq->dark);
		uniq->dark = NULL;
	}
	if (uniq->flat) {
		clearfits(uniq->flat);
		free(uniq->flat);
		uniq->flat = NULL;
	}
}

int single_image_is_loaded() {
	return (com.uniq != NULL && com.uniq->nb_layers > 0);
}

