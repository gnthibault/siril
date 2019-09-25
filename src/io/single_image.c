/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2019 team free-astro (see more in AUTHORS file)
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
#include "algos/statistics.h"
#include "algos/background_extraction.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/message_dialog.h"
#include "io/conversion.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "gui/PSF_list.h"
#include "gui/histogram.h"
#include "gui/progress_and_log.h"
#include "core/undo.h"
#include "core/processing.h"

/* Closes and frees resources attached to the single image opened in gfit.
 * If a sequence is loaded and one of its images is displayed, nothing is done.
 */
void close_single_image() {
	if (sequence_is_loaded() && com.seq.current >= 0)
		return;
	fprintf(stdout, "MODE: closing single image\n");
	/* we need to close all dialogs in order to avoid bugs
	 * with previews
	 */
	if (!com.headless) {
		siril_close_preview_dialogs();
	}
	free_image_data();
	undo_flush();
}

/* frees resources when changing sequence or closing a single image
 * (image size may vary, so we reallocate) */
void free_image_data() {
	int vport;
	fprintf(stdout, "free_image_data() called, clearing loaded image\n");
	/* WARNING: single_image.fit references the actual fits image,
	 * shouldn't it be used here instead of gfit? */
	if (!single_image_is_loaded() && sequence_is_loaded())
		save_stats_from_fit(&gfit, &com.seq, com.seq.current);
	clearfits(&gfit);
	if (!com.headless) {
		clear_stars_list();
		delete_selected_area();
		clear_sampling_setting_box();	// clear focal and pixel pitch info
		free_background_sample_list(com.grad_samples);
		com.grad_samples = NULL;
	}
	clear_histograms();

	for (vport = 0; vport < MAXGRAYVPORT; vport++) {
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
	if (!com.headless)
		activate_tab(RED_VPORT);
	if (com.rgbbuf) {
		free(com.rgbbuf);
		com.rgbbuf = NULL;
	}
	if (com.uniq) {
		free(com.uniq->filename);
		free(com.uniq->comment);
		free(com.uniq->layers);
		free(com.uniq);
		com.uniq = NULL;
	}

	if (!com.headless) {
		/* free alignment preview data */
		int i;
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
}

static gboolean end_read_single_image(gpointer p) {
	set_GUI_CAMERA();
	set_GUI_photometry();
	return FALSE;
}

/* Reads an image from disk and stores it in the user allocated destination
 * fits. The realname_out argument can be NULL, and if not, it is set to the
 * real file name of the loaded file, since the given filename can be without
 * extension. This realname output has to be freed by user.
 */
int read_single_image(const char* filename, fits *dest, char **realname_out) {
	int retval;
	image_type imagetype;
	char *realname = NULL;

	retval = stat_file(filename, &imagetype, &realname);
	if (retval) {
		char *msg = siril_log_message(_("Error opening image %s: file not found or not supported.\n"), filename);
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error"), msg);
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
			debayer_if_needed(imagetype, dest, com.debayer.compatibility, FALSE, com.debayer.stretch);
	}
	if (retval != 0 && retval != 3)
		siril_log_message(_("Opening %s failed.\n"), realname);
	if (realname_out)
		*realname_out = realname;
	else
		free(realname);
	com.filter = (int) imagetype;
	siril_add_idle(end_read_single_image, NULL);
	return retval;
}

static gboolean end_open_single_image(gpointer arg) {
	open_single_image_from_gfit();
	return FALSE;
}

/* This function is used to load a single image, meaning outside a sequence,
 * whether a sequence is loaded or not, whether an image was already loaded or
 * not. The opened file is available in the usual global variable for current
 * image, gfit.
 */
int open_single_image(const char* filename) {
	int retval;
	char *realname;

	close_sequence(FALSE);	// closing a sequence if loaded
	close_single_image();	// close the previous image and free resources

	retval = read_single_image(filename, &gfit, &realname);
	
	/* A single sequence has been loaded */
	if (retval == 3) {
		return 0;
	}
	if (retval == 2) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error opening file"),
				_("This file could not be opened because "
						"its extension is not supported."));
		return 1;
	}
	if (retval == 1) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error opening file"),
				_("There was an error when opening this image. "
						"See the log for more information."));
		return 1;
	}

	fprintf(stdout, "Loading image OK, now displaying\n");

	/* Now initializing com struct */
	com.seq.current = UNRELATED_IMAGE;
	com.uniq = calloc(1, sizeof(single));
	com.uniq->filename = realname;
	com.uniq->nb_layers = gfit.naxes[2];
	com.uniq->layers = calloc(com.uniq->nb_layers, sizeof(layer_info));
	com.uniq->fit = &gfit;
	siril_add_idle(end_open_single_image, realname);
	return 0;
}

/* creates a single_image structure and displays a single image, found in gfit.
 */
void open_single_image_from_gfit() {
	/* now initializing everything
	 * code based on seq_load_image or set_seq (sequence.c) */

	initialize_display_mode();

	init_layers_hi_and_lo_values(MIPSLOHI);		// If MIPS-LO/HI exist we load these values. If not it is min/max

	sliders_mode_set_state(com.sliders);
	set_cutoff_sliders_max_values();
	set_cutoff_sliders_values();

	set_display_mode();
	update_prepro_interface(TRUE);
	adjust_sellabel();

	display_filename();	// display filename in gray window

	/* update menus */
	update_MenuItem();

	show_main_gray_window();
	close_tab();
	if (gfit.naxes[2] == 3) {
		show_rgb_window();
	} else {
		hide_rgb_window();
	}
	update_gfit_histogram_if_needed();
	adjust_vport_size_to_image();
	redraw(com.cvport, REMAP_ALL);
	update_used_memory();
}

/* searches the image for minimum and maximum pixel value, on each layer
 * the values are stored in fit->min[layer] and fit->max[layer] */
int image_find_minmax(fits *fit) {
	int layer;
	if (fit->maxi > 0.0)
		return 0;
	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		// calling statistics() saves stats in the fit already, we don't need
		// to use the returned handle
		free_stats(statistics(NULL, -1, fit, layer, NULL, STATS_MINMAX));
		if (!fit->stats[layer])
			return -1;
		fit->maxi = max(fit->maxi, fit->stats[layer]->max);
		fit->mini = min(fit->mini, fit->stats[layer]->min);
	}
	return 0;
}

static int fit_get_minmax(fits *fit, int layer) {
	// calling statistics() saves stats in the fit already, we don't need
	// to use the returned handle
	free_stats(statistics(NULL, -1, fit, layer, NULL, STATS_MINMAX));
	if (!fit->stats[layer])
		return -1;
	return 0;
}

double fit_get_max(fits *fit, int layer) {
	if (fit_get_minmax(fit, layer))
		return -1.0;
	return fit->stats[layer]->max;
}

double fit_get_min(fits *fit, int layer) {
	if (fit_get_minmax(fit, layer))
		return -1.0;
	return fit->stats[layer]->min;
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
				layers[i].hi = fit_get_max(&gfit, i);
				layers[i].lo = fit_get_min(&gfit, i);
			}
			else {
				image_find_minmax(&gfit);
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
	invalidate_stats_from_fit(&gfit);
	if (!com.script) {
		invalidate_gfit_histogram();
		compute_histo_for_gfit();
		init_layers_hi_and_lo_values(com.sliders);
		set_cutoff_sliders_values();
	}
}

int single_image_is_loaded() {
	return (com.uniq != NULL && com.uniq->nb_layers > 0);
}

