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

#include <stdlib.h>
#include <string.h>
#include "core/siril.h"
#include "core/proto.h"

#include "registration.h"
#include "algos/PSF.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "core/processing.h"
#include "opencv/opencv.h"
#include "gui/image_interactions.h"
#include "gui/image_display.h"
#include "gui/utils.h"

static int awaiting_star = 0;

static GtkWidget *three_buttons[3] = { 0 };

struct _3psf {
	fitted_PSF *stars[3];
};

static struct _3psf *results;
static int results_size;

// local functions
static int rotate_images(struct registration_args *regargs, regdata *current_regdata);

static void set_registration_ready(gboolean ready) {
	static GtkWidget *go_register = NULL;
	if (!go_register)
		go_register = lookup_widget("goregister_button");
	gtk_widget_set_sensitive(go_register, ready);
}

static void update_label(gchar* str) {
	static GtkLabel *labelreginfo = NULL;
	if (!labelreginfo)
		labelreginfo = GTK_LABEL(lookup_widget("labelregisterinfo"));
	gtk_label_set_text(labelreginfo, str);
}

static void update_icons(int idx, gboolean OK) {
	static GtkImage *image_3stars[3] = { NULL };

	if (!image_3stars[0]) {
		image_3stars[0] = GTK_IMAGE(lookup_widget("3stars-image1"));
		image_3stars[1] = GTK_IMAGE(lookup_widget("3stars-image2"));
		image_3stars[2] = GTK_IMAGE(lookup_widget("3stars-image3"));
	}
	gtk_image_set_from_icon_name(image_3stars[idx],
			OK ? "gtk-yes" : "gtk-no", GTK_ICON_SIZE_LARGE_TOOLBAR);

}

static void reset_icons() {
	for (int i = 0; i < 3; i++) {
		update_icons(i, FALSE);
	}
}

static gboolean _3stars_seqpsf_end(gpointer p) {
	/* the fun part, synchronizing the three threads */
	struct generic_seq_args *args = (struct generic_seq_args *)p;
	struct seqpsf_args *spsfargs = (struct seqpsf_args *)args->user;

	if (args->retval) {
		if (args->seq->current != 0)
			update_label(_("Make sure you load the first image"));
		else update_label(_("Star analysis failed"));
		goto psf_end;
	}

	GSList *iterator;
	for (iterator = spsfargs->list; iterator; iterator = iterator->next) {
		struct seqpsf_data *data = iterator->data;
		results[data->image_index].stars[awaiting_star - 1] = data->psf;
	}
	g_slist_free(spsfargs->list);
	int refimage = sequence_find_refimage(&com.seq);
	if (!results[refimage].stars[awaiting_star - 1]) {
		siril_log_color_message(_("The star was not found in the reference image. Change the selection or the reference image\n"), "red");
		for (int i = 0 ; i < com.seq.number; i++)
			results[i].stars[awaiting_star - 1] = NULL;
		goto psf_end;
	}

	unset_suggested(three_buttons[awaiting_star - 1]);
	int i;
	for (i = 0; i < 3 && results[args->seq->current].stars[i]; i++);
	if (i < 3) {
		set_suggested(three_buttons[i]);
		set_registration_ready(i == 2);
	} else {
		set_registration_ready(TRUE);
	}
	update_icons(awaiting_star - 1, TRUE);

	com.stars = realloc(com.stars, 4 * sizeof(fitted_PSF *)); // to be sure...
	com.stars[awaiting_star - 1] = duplicate_psf(results[args->seq->current].stars[awaiting_star - 1]);

psf_end:
	memset(&com.selection, 0, sizeof(rectangle));
	redraw(com.cvport, REMAP_NONE);

	awaiting_star = 0;
	free(args);
	free(spsfargs);
	return end_generic(NULL);
}

static void start_seqpsf() {
	struct seqpsf_args *spsfargs = malloc(sizeof(struct seqpsf_args));
	spsfargs->for_registration = TRUE; // if false, photometry is computed
	spsfargs->framing = FOLLOW_STAR_FRAME;
	spsfargs->list = NULL;	// GSList init is NULL
	struct generic_seq_args *args = calloc(1, sizeof(struct generic_seq_args));
	args->seq = &com.seq;
	args->partial_image = TRUE;
	args->layer_for_partial = get_registration_layer(&com.seq);
	args->regdata_for_partial = FALSE;
	args->get_photometry_data_for_partial = FALSE;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = com.seq.selnum;
	args->image_hook = seqpsf_image_hook;
	args->idle_function = _3stars_seqpsf_end;
	args->stop_on_error = FALSE;
	args->description = _("PSF on area for 3 stars");
	args->upscale_ratio = 1.0;
	args->user = spsfargs;
	args->already_in_a_thread = FALSE;
	args->parallel = FALSE;	// follow star implies not parallel
	memcpy(&args->area, &com.selection, sizeof(rectangle));
	if (!results) {
		results = calloc(com.seq.number, sizeof(struct _3psf));
		if (!results) {
			PRINT_ALLOC_ERR;
			return;
		}
		results_size = com.seq.number;
	}

	start_in_new_thread(generic_sequence_worker, args);
}

void on_select_star_button_clicked(GtkButton *button, gpointer user_data) {
	if (!three_buttons[0]) {
		three_buttons[0] = lookup_widget("pickstar1");
		three_buttons[1] = lookup_widget("pickstar2");
		three_buttons[2] = lookup_widget("pickstar3");
	}
	if (!com.selection.w || !com.selection.h) {
		update_label(_("draw a selection around the star"));
		return;
	}
	GtkWidget *widget = GTK_WIDGET(button);
	if (three_buttons[0] == widget)
		awaiting_star = 1;
	else if (three_buttons[1] == widget)
		awaiting_star = 2;
	else if (three_buttons[2] == widget)
		awaiting_star = 3;
	else {
	       fprintf(stderr, "unknown button clicked\n");
       	       return;
	}

	if (!com.stars)
		com.stars = calloc(4, sizeof(fitted_PSF *));

	start_seqpsf();
}

int register_3stars(struct registration_args *regargs) {
	int refimage = regargs->reference_image;
	if (!results[refimage].stars[0] || !results[refimage].stars[1]) {
		siril_log_color_message("Less than two stars were found in the reference image, try setting another as reference?\n", "red");
		return 1;
	}

	regdata *current_regdata = star_align_get_current_regdata(regargs);
	if (!current_regdata) return -2;

	/* set regparams for current sequence before closing it */
	for (int i = 0; i < regargs->seq->number; i++) {
		double sumx = 0.0, sumy = 0.0;
		int nb_stars = 0;
		if (results[i].stars[0]) {
			sumx += results[i].stars[0]->fwhmx;
			sumy += results[i].stars[0]->fwhmy;
			nb_stars++;
		}
		if (results[i].stars[1]) {
			sumx += results[i].stars[1]->fwhmx;
			sumy += results[i].stars[1]->fwhmy;
			nb_stars++;
		}
		if (results[i].stars[2]) {
			sumx += results[i].stars[2]->fwhmx;
			sumy += results[i].stars[2]->fwhmy;
			nb_stars++;
		}
		double fwhm = sumx / nb_stars;
		current_regdata[i].roundness = sumy / sumx;
		current_regdata[i].fwhm = fwhm;
		current_regdata[i].weighted_fwhm = fwhm; // TODO: compute it with nb_stars
	}

	return rotate_images(regargs, current_regdata);
}

/* image rotation sequence processing */
static int affine_transform_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *area) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	int refimage = regargs->reference_image;

	int nb_stars = 3;
	if (!results[refimage].stars[2])
		nb_stars = 2;
	if (nb_stars == 2 && (!results[in_index].stars[0] || !results[in_index].stars[1]))
		return 1;
	if (nb_stars == 2 || (nb_stars == 3 && results[in_index].stars[0] && results[in_index].stars[1] && results[in_index].stars[2])) {
		if (regargs->x2upscale || in_index != refimage) {
			pointf ref[3] = {
				{ results[refimage].stars[0]->xpos, results[refimage].stars[0]->ypos },
				{ results[refimage].stars[1]->xpos, results[refimage].stars[1]->ypos },
				{ 0 }
			};
			if (nb_stars == 3) {
				ref[2].x = results[refimage].stars[2]->xpos;
				ref[2].y = results[refimage].stars[2]->ypos;
			}
			pointf cur[3] = {
				{ results[in_index].stars[0]->xpos, results[in_index].stars[0]->ypos },
				{ results[in_index].stars[1]->xpos, results[in_index].stars[1]->ypos },
				{ 0 }
			};
			if (nb_stars == 3) {
				cur[2].x = results[in_index].stars[2]->xpos;
				cur[2].y = results[in_index].stars[2]->ypos;
			}

			if (cvAffineTransformation(fit, ref, cur, nb_stars, regargs->x2upscale, regargs->interpolation))
				return 1;
		}
	}
	else {
		int in_stars = (results[in_index].stars[0] != NULL) +
			(results[in_index].stars[1] != NULL) + (results[in_index].stars[2] != NULL);
		if (in_stars != 2)
			return 1;
		pointf ref[2] = { 0 };
		pointf cur[2] = { 0 };
		int star = 0;
		if (results[in_index].stars[0]) {
			ref[star].x = results[refimage].stars[0]->xpos;
			ref[star].y = results[refimage].stars[0]->ypos;
			cur[star].x = results[in_index].stars[0]->xpos;
			cur[star].y = results[in_index].stars[0]->ypos;
			star++;
		}
		if (results[in_index].stars[1]) {
			ref[star].x = results[refimage].stars[1]->xpos;
			ref[star].y = results[refimage].stars[1]->ypos;
			cur[star].x = results[in_index].stars[1]->xpos;
			cur[star].y = results[in_index].stars[1]->ypos;
			star++;
		}
		if (results[in_index].stars[2]) {
			ref[star].x = results[refimage].stars[2]->xpos;
			ref[star].y = results[refimage].stars[2]->ypos;
			cur[star].x = results[in_index].stars[2]->xpos;
			cur[star].y = results[in_index].stars[2]->ypos;
			star++;
		}

		if (cvAffineTransformation(fit, ref, cur, nb_stars, regargs->x2upscale, regargs->interpolation))
			return 1;
	}

	sadata->success[out_index] = 1;
	regargs->imgparam[out_index].filenum = args->seq->imgparam[in_index].filenum;
	regargs->imgparam[out_index].incl = SEQUENCE_DEFAULT_INCLUDE;
	regargs->regparam[out_index].fwhm = sadata->current_regdata[in_index].fwhm;
	regargs->regparam[out_index].weighted_fwhm = sadata->current_regdata[in_index].weighted_fwhm;
	regargs->regparam[out_index].roundness = sadata->current_regdata[in_index].roundness;

	if (regargs->x2upscale) {
		fit->pixel_size_x /= 2;
		fit->pixel_size_y /= 2;
		regargs->regparam[out_index].fwhm *= 2.0;
		regargs->regparam[out_index].weighted_fwhm *= 2.0;
	}
	return 0;
}

static int affine_transform_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) {
	unsigned int MB_per_image; int MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, args->upscale_ratio, FALSE, &MB_per_image, &MB_avail);
	if (limit > 0) {
		/* The registration memory consumption, n is image size:
		 * for monochrome images O(n), O(2n) for color
		 * Don't forget that this is in addition to the image being already read.
		 */
		unsigned int required;
		if (args->seq->nb_layers == 3)
			required = MB_per_image * 3;
		else required = MB_per_image;
		args->max_thread = MB_avail / required;
		siril_debug_print("Memory required per thread: %u MB, limiting to %d threads\n", required, args->max_thread);
	}

	if (limit == 0) {
		gchar *mem_per_image = g_format_size_full(MB_per_image * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);

		siril_log_color_message(_("%s: not enough memory to do this operation (%s required per image, %s considered available)\n"),
				"red", args->description, mem_per_image, mem_available);

		g_free(mem_per_image);
		g_free(mem_available);
	} else {
#ifdef _OPENMP
		if (for_writer) {
			/* we already have allowed limit [1, max_thread] for
			 * the non-writer case */
			limit -= com.max_thread;
			if (limit < 0) limit = 0;
			int max_queue_size = com.max_thread * 3;
			if (limit > max_queue_size)
				limit = max_queue_size;
		}
		else if (limit > com.max_thread)
			limit = com.max_thread;
#else
		if (for_writer) {
			limit--;
			if (limit < 0) limit = 0;
			if (limit > 3)
				limit = 3;
		} else {
			limit = 1;
		}
#endif
	}
	return limit;
}

static int rotate_images(struct registration_args *regargs, regdata *current_regdata) {
	struct generic_seq_args *args = create_default_seqargs(&com.seq);
	args->stop_on_error = FALSE;
	if (!regargs->process_all_frames) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = regargs->seq->selnum;
	}
	args->compute_mem_limits_hook = affine_transform_compute_mem_limits;
	args->prepare_hook = star_align_prepare_results;
	args->image_hook = affine_transform_hook;
	args->finalize_hook = star_align_finalize_hook;	// from global registration
	args->description = _("Creating the rotated image sequence");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->upscale_ratio = regargs->x2upscale ? 2.0 : 1.0;
	args->new_seq_prefix = regargs->prefix;
	args->load_new_sequence = TRUE;
	args->already_in_a_thread = TRUE;

	struct star_align_data *sadata = calloc(1, sizeof(struct star_align_data));
	if (!sadata) {
		free(args);
		return -1;
	}
	sadata->regargs = regargs;
	// we pass the regdata just to avoid recomputing it for the new sequence
	sadata->current_regdata = current_regdata; 
	args->user = sadata;

	generic_sequence_worker(args);
	
	for (int i = 0; i < results_size; i++) {
		for (int s = 0; s < 3; s++)
			if (results[i].stars[s])
				free(results[i].stars[s]);
	}
	free(results);
	results = NULL;
	reset_icons();
	for (int i = 0; i < 3; i++)
		unset_suggested(three_buttons[i]);
	set_suggested(three_buttons[0]);
	return args->retval;
}

