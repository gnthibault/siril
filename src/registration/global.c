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
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "algos/star_finder.h"
#include "algos/statistics.h"
#include "algos/PSF.h"
#include "gui/image_display.h"
#include "gui/PSF_list.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/image_format_fits.h"
#include "registration/registration.h"
#include "registration/matching/atpmatch.h"
#include "registration/matching/match.h"
#include "registration/matching/misc.h"
#include "opencv/opencv.h"


/* TODO:
 * check usage of openmp in functions called by these ones (to be disabled)
 * compact and clarify console output
 */

static void create_output_sequence_for_global_star(struct registration_args *args);
static void print_alignment_results(Homography H, int filenum, float FWHMx, float FWHMy, char *units);

regdata *star_align_get_current_regdata(struct registration_args *regargs) {
	regdata *current_regdata;
	if (regargs->seq->regparam[regargs->layer]) {
		siril_log_message(
				_("Recomputing already existing registration for this layer\n"));
		current_regdata = regargs->seq->regparam[regargs->layer];
		/* we reset all values as we may register different images */
		memset(current_regdata, 0, regargs->seq->number * sizeof(regdata));
	} else {
		current_regdata = calloc(regargs->seq->number, sizeof(regdata));
		if (current_regdata == NULL) {
			PRINT_ALLOC_ERR;
			return NULL;
		}
		regargs->seq->regparam[regargs->layer] = current_regdata;
	}
	return current_regdata;
}

int star_align_prepare_results(struct generic_seq_args *args) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;

	if (!regargs->translation_only) {
		// allocate destination sequence data
		regargs->imgparam = calloc(args->nb_filtered_images, sizeof(imgdata));
		regargs->regparam = calloc(args->nb_filtered_images, sizeof(regdata));
		if (!regargs->imgparam  || !regargs->regparam) {
			PRINT_ALLOC_ERR;
			return 1;
		}

		if (args->seq->type == SEQ_SER) {
			/* copied from seq_prepare_hook with one variation */
			args->new_ser = malloc(sizeof(struct ser_struct));
			if (!args->new_ser) {
				PRINT_ALLOC_ERR;
				return 1;
			}

			char dest[256];
			const char *ptr = strrchr(args->seq->seqname, G_DIR_SEPARATOR);
			if (ptr)
				snprintf(dest, 255, "%s%s.ser", regargs->prefix, ptr + 1);
			else
				snprintf(dest, 255, "%s%s.ser", regargs->prefix, args->seq->seqname);

			/* Here the last argument is NULL because we do not want copy SER file
			 * from the original. Indeed in the demosaicing case this would lead to
			 * a wrong file (B&W and not in RAW data). Moreover, header informations
			 * (like fps, local and UTC time, ...) have no sense now since some frames
			 * could be removed from the sequence.
			 */
			if (ser_create_file(dest, args->new_ser, TRUE, NULL)) {
				free(args->new_ser);
				args->new_ser = NULL;
				return 1;
			}

			if (seq_prepare_writer(args))
				return 1;
		}
		else if (args->seq->type == SEQ_FITSEQ) {
			if (seq_prepare_hook(args))
				return 1;
		}
	}

	sadata->success = calloc(args->nb_filtered_images, sizeof(BYTE));
	if (!sadata->success) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	return 0;
}


static int star_align_prepare_hook(struct generic_seq_args *args) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	float FWHMx, FWHMy;
	char *units;
	fits fit = { 0 };
	int i, nb_stars = 0;

	sadata->current_regdata = star_align_get_current_regdata(regargs);
	if (!sadata->current_regdata) return -2;
	
	/* first we're looking for stars in reference image */
	if (seq_read_frame(args->seq, regargs->reference_image, &fit, FALSE, -1)) {
		siril_log_message(_("Could not load reference image\n"));
		args->seq->regparam[regargs->layer] = NULL;
		free(sadata->current_regdata);
		return 1;
	}

	siril_log_color_message(_("Reference Image:\n"), "green");

	if (regargs->matchSelection && regargs->selection.w > 0 && regargs->selection.h > 0) {
		com.stars = peaker(&fit, regargs->layer, &com.starfinder_conf, &nb_stars, &regargs->selection, FALSE, TRUE);
	}
	else {
		com.stars = peaker(&fit, regargs->layer, &com.starfinder_conf, &nb_stars, NULL, FALSE, TRUE);
	}

	siril_log_message(_("Found %d stars in reference, channel #%d\n"), nb_stars, regargs->layer);


	if (!com.stars || nb_stars < AT_MATCH_MINPAIRS) {
		siril_log_message(
				_("There are not enough stars in reference image to perform alignment\n"));
		args->seq->regparam[regargs->layer] = NULL;
		free(sadata->current_regdata);
		clearfits(&fit);
		return 1;
	}
	if (!com.script && &com.seq == args->seq && com.seq.current == regargs->reference_image)
		queue_redraw(REMAP_NONE); // draw stars

	sadata->ref.x = fit.rx;
	sadata->ref.y = fit.ry;

	clearfits(&fit);

	if (regargs->x2upscale) {
		if (regargs->translation_only) {
			args->seq->upscale_at_stacking = 2.0;
		} else {
			sadata->ref.x *= 2.0;
			sadata->ref.y *= 2.0;
		}
	}
	else {
		if (regargs->translation_only) {
			args->seq->upscale_at_stacking = 1.0;
		}
	}

	/* we copy com.stars to refstars in case user take a look to another image of the sequence
	 * that would destroy com.stars
	 */
	i = 0;
	sadata->refstars = new_fitted_stars(MAX_STARS);
	if (!sadata->refstars) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	while (i < MAX_STARS && com.stars[i]) {
		fitted_PSF *tmp = malloc(sizeof(fitted_PSF));
		if (!tmp) {
			PRINT_ALLOC_ERR;
			sadata->refstars[i] = NULL;
			return 1;
		}
		memcpy(tmp, com.stars[i], sizeof(fitted_PSF));
		sadata->refstars[i] = tmp;
		sadata->refstars[i+1] = NULL;
		i++;
	}

	if (nb_stars >= MAX_STARS_FITTED) {
		sadata->fitted_stars = MAX_STARS_FITTED;
		siril_log_color_message(_("Reference Image: Limiting to %d brightest stars\n"), "green", MAX_STARS_FITTED);
	} else {
		sadata->fitted_stars = nb_stars;
	}
	FWHM_average(sadata->refstars, sadata->fitted_stars, &FWHMx, &FWHMy, &units);
	siril_log_message(_("FWHMx:%*.2f %s\n"), 12, FWHMx, units);
	siril_log_message(_("FWHMy:%*.2f %s\n"), 12, FWHMy, units);
	sadata->current_regdata[regargs->reference_image].roundness = FWHMy/FWHMx;
	sadata->current_regdata[regargs->reference_image].fwhm = FWHMx;
	sadata->current_regdata[regargs->reference_image].weighted_fwhm = FWHMx;

	return star_align_prepare_results(args);
}

/* reads the image, searches for stars in it, tries to match them with
 * reference stars, computes the homography matrix, applies it on the image,
 * possibly up-scales the image and stores registration data */
static int star_align_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	int nbpoints, nb_stars = 0;
	int retvalue;
	int nobj = 0;
	int attempt = 1;
	float FWHMx, FWHMy;
	char *units;
	Homography H = { 0 };
	int filenum = args->seq->imgparam[in_index].filenum;	// for display purposes

	if (regargs->translation_only) {
		/* if "translation only", we choose to initialize all frames
		 * to exclude status. If registration is ok, the status is
		 * set to include */
		args->seq->imgparam[out_index].incl = !SEQUENCE_DEFAULT_INCLUDE;
	}

	if (in_index != regargs->reference_image) {
		fitted_PSF **stars;
		if (args->seq->type == SEQ_SER || args->seq->type == SEQ_FITSEQ) {
			siril_log_color_message(_("Frame %d:\n"), "bold", filenum);
		}

		if (regargs->matchSelection && regargs->selection.w > 0 && regargs->selection.h > 0) {
			stars = peaker(fit, regargs->layer, &com.starfinder_conf, &nb_stars, &regargs->selection, FALSE, TRUE);
		}
		else {
			stars = peaker(fit, regargs->layer, &com.starfinder_conf, &nb_stars, NULL, FALSE, TRUE);
		}

		siril_log_message(_("Found %d stars in image %d, channel #%d\n"), nb_stars, filenum, regargs->layer);

		if (!stars || nb_stars < AT_MATCH_MINPAIRS) {
			siril_log_message(
					_("Not enough stars. Image %d skipped\n"), filenum);
			if (stars) free_fitted_stars(stars);
			return 1;
		}

		if (nb_stars >= sadata->fitted_stars) {
			if (nb_stars >= MAX_STARS_FITTED) {
				siril_log_color_message(_("Target Image: Limiting to %d brightest stars\n"), "green", MAX_STARS_FITTED);
			}
			nbpoints = sadata->fitted_stars;
		}
		else {
			nbpoints = nb_stars;
		}

		/* make a loop with different tries in order to align the two sets of data */
		double scale_min = 0.9;
		double scale_max = 1.1;
		retvalue = 1;
		while (retvalue && attempt < NB_OF_MATCHING_TRY){
			retvalue = new_star_match(stars, sadata->refstars, nbpoints, nobj, scale_min, scale_max, &H, FALSE);
			if (attempt == 1) {
				scale_min = -1.0;
				scale_max = -1.0;
			} else {
				nobj += 50;
			}
			attempt++;
		}
		if (retvalue) {
			siril_log_color_message(_("Cannot perform star matching: try #%d. Image %d skipped\n"),
					"red", attempt, filenum);
			free_fitted_stars(stars);
			return 1;
		}
		if (H.pair_matched < regargs->min_pairs) {
			siril_log_color_message(_("Not enough star pairs (%d): Image %d skipped\n"),
					"red", H.pair_matched, filenum);
			free_fitted_stars(stars);
			return 1;
		}

		FWHM_average(stars, nbpoints, &FWHMx, &FWHMy, &units);
#ifdef _OPENMP
#pragma omp critical
#endif
		print_alignment_results(H, filenum, FWHMx, FWHMy, units);

		sadata->current_regdata[in_index].roundness = FWHMy/FWHMx;
		sadata->current_regdata[in_index].fwhm =  FWHMx;
		sadata->current_regdata[in_index].weighted_fwhm = 2 * FWHMx
				* (((double) sadata->fitted_stars) - (double) H.pair_matched)
				/ (double) sadata->fitted_stars + FWHMx;

		if (!regargs->translation_only) {
			if (cvTransformImage(fit, sadata->ref.x, sadata->ref.y, H, regargs->x2upscale, regargs->interpolation)) {
				free_fitted_stars(stars);
				return 1;
			}
		}

		free_fitted_stars(stars);
	}
	else {
		if (regargs->x2upscale && !regargs->translation_only) {
			if (cvResizeGaussian(fit, fit->rx * 2, fit->ry * 2, OPENCV_NEAREST))
				return 1;
		}
	}

	if (!regargs->translation_only) {
		regargs->imgparam[out_index].filenum = args->seq->imgparam[in_index].filenum;
		regargs->imgparam[out_index].incl = SEQUENCE_DEFAULT_INCLUDE;
		regargs->regparam[out_index].fwhm = sadata->current_regdata[in_index].fwhm;	// not FWHMx because of the ref frame
		regargs->regparam[out_index].weighted_fwhm = sadata->current_regdata[in_index].weighted_fwhm;
		regargs->regparam[out_index].roundness = sadata->current_regdata[in_index].roundness;

		if (regargs->x2upscale) {
			fit->pixel_size_x /= 2;
			fit->pixel_size_y /= 2;
			regargs->regparam[out_index].fwhm *= 2.0;
			regargs->regparam[out_index].weighted_fwhm *= 2.0;
		}
	} else {
		set_shifts(args->seq, in_index, regargs->layer, (float) H.h02,
				(float) -H.h12, fit->top_down);
		args->seq->imgparam[out_index].incl = SEQUENCE_DEFAULT_INCLUDE;
	}
	sadata->success[out_index] = 1;
	return 0;
}

int star_align_finalize_hook(struct generic_seq_args *args) {
	struct star_align_data *sadata = args->user;
	struct registration_args *regargs = sadata->regargs;
	int i, failed = 0;

	// images may have been excluded but selnum wasn't updated
	fix_selnum(args->seq, FALSE);

	free_fitted_stars(sadata->refstars);

	if (!args->retval) {
		for (i = 0; i < args->nb_filtered_images; i++)
			if (!sadata->success[i])
				failed++;
		regargs->new_total = args->nb_filtered_images - failed;

		if (!regargs->translation_only) {
			if (failed) {
				// regargs->imgparam and regargs->regparam may have holes caused by images
				// that failed to be registered - compact them
				int j;
				for (i = 0, j = 0; i < regargs->new_total; i++, j++) {
					while (!sadata->success[j] && j < args->nb_filtered_images) j++;
					g_assert(sadata->success[j]);
					if (i != j) {
						regargs->imgparam[i] = regargs->imgparam[j];
						regargs->regparam[i] = regargs->regparam[j];
					}
				}
			}

			seq_finalize_hook(args);
		}
	} else {
		regargs->new_total = 0;
		free(args->seq->regparam[regargs->layer]);
		args->seq->regparam[regargs->layer] = NULL;

		// args->new_ser can be null if stars were not detected in the reference image
		// same as seq_finalize_hook but with file deletion
		if ((args->force_ser_output || args->seq->type == SEQ_SER) && args->new_ser) {
			ser_close_and_delete_file(args->new_ser);
			free(args->new_ser);
		}
		else if ((args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) && args->new_fitseq) {
			fitseq_close_and_delete_file(args->new_fitseq);
			free(args->new_fitseq);
		}
	}

	if (sadata->success) free(sadata->success);
	free(sadata);
	args->user = NULL;
	clear_stars_list();

	if (!args->retval) {
		siril_log_message(_("Registration finished.\n"));
		gchar *str = ngettext("%d image processed.\n", "%d images processed.\n", args->nb_filtered_images);
		str = g_strdup_printf(str, args->nb_filtered_images);
		siril_log_color_message(str, "green");
		siril_log_color_message(_("Total: %d failed, %d registered.\n"), "green", failed, regargs->new_total);

		g_free(str);
		if (!regargs->translation_only) {
			// explicit sequence creation to copy imgparam and regparam
			create_output_sequence_for_global_star(regargs);
			// will be loaded in the idle function if (load_new_sequence)
			regargs->load_new_sequence = TRUE; // only case where a new sequence must be loaded
		}
	}
	else {
		siril_log_message(_("Registration aborted.\n"));
	}
	return regargs->new_total == 0;
	// TODO: args is never freed because we don't call an end function for
	// this generic processing function. The register idle is called for
	// everything else, but does not know this pointer, and we cannot free
	// it here because it's still used in the generic processing function.
}

int star_align_compute_mem_limits(struct generic_seq_args *args, gboolean for_writer) { 
	unsigned int MB_per_orig_image, MB_per_scaled_image, MB_avail;
	int limit = compute_nb_images_fit_memory(args->seq, args->upscale_ratio, args->force_float,
			&MB_per_orig_image, &MB_per_scaled_image, &MB_avail);
	unsigned int required = MB_per_scaled_image;
	if (limit > 0) {
		/* The registration memory consumption, n is image size and m channel size.
		 * First, a threshold is computed for star pixel value, using statistics:
		 *	O(m), data is duplicated for median computation if
		 *	there are nil values, O(1) otherwise
		 * Then, still in peaker(), image is filtered using unsharp filter, duplicating
		 * the reference channel to act as input and output of the filter as float O(2m
		 * as float).
		 * Then, the image is rotated and upscaled by the generic function if enabled:
		 * cvTransformImage is O(n) in mem for unscaled, O(nscaled)=O(4m) for
		 * monochrome scaled and O(2nscaled)=O(21m) for color scaled
		 * All this is in addition to the image being already loaded, except for the
		 * color scaled image.
		 *
		 * Since these three operations are in sequence, we need room only for the
		 * largest.
		 * rotated color scaled float	mem needed
		 *       0     0      0     0	O(2m as float)
		 *       1     0      0     0	O(2m as float)
		 *       1     0      0     1	O(2m as float, same as 2n)
		 *       1     0      1     0	O(4m, same as 2m as float)
		 *       1     0      1     1	O(4m)
		 *       1     1      0     0	O(2m as float)
		 *       1     1      0     1	O(n)
		 *       1     1      1     0	O(8n or 2nscaled)
		 *       1     1      1     1	O(8n or 2nscaled)
		 */
		int is_color = args->seq->nb_layers == 3;
		int is_float = get_data_type(args->seq->bitpix) == DATA_FLOAT;
		int is_scaled = args->upscale_ratio == 2.0;
		unsigned int float_multiplier = is_float ? 1 : 2;
		unsigned int MB_per_float_image = MB_per_orig_image * float_multiplier;
		unsigned int MB_per_float_channel = is_color ? MB_per_float_image / 3 : MB_per_float_image;
		unsigned int MB_per_orig_channel = is_color ? MB_per_orig_image / 3 : MB_per_float_image;
		MB_per_float_channel = min(1, MB_per_float_channel);
		MB_per_orig_channel = min(1, MB_per_orig_channel);
		if (!args->has_output || (!is_scaled && (!is_color || !is_float))) {
			required = MB_per_orig_image + MB_per_float_channel * 2;
		}
		else if (args->has_output && !is_color && is_scaled) {
			required = MB_per_orig_image + 4 * MB_per_orig_channel;
		}
		else if (args->has_output && is_color && !is_scaled && is_float) {
			required = 2 * MB_per_orig_image;
		}
		else {
			required = 2 * MB_per_scaled_image;
		}
		int thread_limit = MB_avail / required;
		if (thread_limit > com.max_thread)
			thread_limit = com.max_thread;

		if (for_writer) {
			/* we allow the already allocated thread_limit images,
			 * plus how many images can be stored in what remains
			 * unused by the main processing */
			limit = thread_limit + (MB_avail - required * thread_limit) / MB_per_scaled_image;
		} else limit = thread_limit;
	}

	if (limit == 0) {
		gchar *mem_per_thread = g_format_size_full(required * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);
		gchar *mem_available = g_format_size_full(MB_avail * BYTES_IN_A_MB, G_FORMAT_SIZE_IEC_UNITS);

		siril_log_color_message(_("%s: not enough memory to do this operation (%s required per thread, %s considered available)\n"),
				"red", args->description, mem_per_thread, mem_available);

		g_free(mem_per_thread);
		g_free(mem_available);
	} else {
#ifdef _OPENMP
		if (for_writer) {
			int max_queue_size = com.max_thread * 3;
			if (limit > max_queue_size)
				limit = max_queue_size;
		}
		siril_debug_print("Memory required per thread: %u MB, per image: %u MB, limiting to %d %s\n",
				required, MB_per_scaled_image, limit, for_writer ? "images" : "threads");
#else
		if (!for_writer)
			limit = 1;
		else if (limit > 3)
			limit = 3;
#endif
	}
	return limit;
}

int register_star_alignment(struct registration_args *regargs) {
	struct generic_seq_args *args = create_default_seqargs(regargs->seq);
	if (!regargs->process_all_frames) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = regargs->seq->selnum;
	}
	args->compute_mem_limits_hook = star_align_compute_mem_limits;
	args->prepare_hook = star_align_prepare_hook;
	args->image_hook = star_align_image_hook;
	args->finalize_hook = star_align_finalize_hook;
	args->stop_on_error = FALSE;
	args->description = _("Global star registration");
	args->has_output = !regargs->translation_only;
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
	args->user = sadata;

	generic_sequence_worker(args);

	regargs->retval = args->retval;
	free(args);
	return regargs->retval;
}

static void create_output_sequence_for_global_star(struct registration_args *args) {
	sequence seq = { 0 };
	initialize_sequence(&seq, TRUE);

	/* we are not interested in the whole path */
	gchar *seqname = g_path_get_basename(args->seq->seqname);
	char *rseqname = malloc(
			strlen(args->prefix) + strlen(seqname) + 5);
	sprintf(rseqname, "%s%s.seq", args->prefix, seqname);
	g_free(seqname);
	g_unlink(rseqname);	// remove previous to overwrite
	args->new_seq_name = remove_ext_from_filename(rseqname);
	free(rseqname);
	seq.seqname = strdup(args->new_seq_name);
	seq.number = args->new_total;
	seq.selnum = args->new_total;
	seq.fixed = args->seq->fixed;
	seq.nb_layers = args->seq->nb_layers;
	seq.rx = args->seq->rx;
	seq.ry = args->seq->ry;
	seq.imgparam = args->imgparam;
	seq.regparam = calloc(seq.nb_layers, sizeof(regdata*));
	seq.regparam[args->layer] = args->regparam;
	seq.layers = calloc(seq.nb_layers, sizeof(layer_info));
	seq.beg = seq.imgparam[0].filenum;
	seq.end = seq.imgparam[seq.number-1].filenum;
	seq.type = args->seq->type;
	seq.current = -1;
	// don't copy from old sequence, it may not be the same image
	seq.reference_image = sequence_find_refimage(&seq);
	seq.needs_saving = TRUE;
	writeseqfile(&seq);
	free_sequence(&seq, FALSE);
}

static void print_alignment_results(Homography H, int filenum, float FWHMx, float FWHMy, char *units) {
	double rotation, scale, scaleX, scaleY;
	point shift;
	double inliers;

	/* Matching information */
	siril_log_color_message(_("Matching stars in image %d: done\n"), "green", filenum);
	gchar *str = ngettext("%d pair match.\n", "%d pair matches.\n", H.pair_matched);
	str = g_strdup_printf(str, H.pair_matched);
	siril_log_message(str);
	g_free(str);
	inliers = 1.0 - ((((double) H.pair_matched - (double) H.Inliers)) / (double) H.pair_matched);
	siril_log_message(_("Inliers:%*.3f\n"), 11, inliers);

	/* Scale */
	scaleX = sqrt(H.h00 * H.h00 + H.h01 * H.h01);
	scaleY = sqrt(H.h10 * H.h10 + H.h11 * H.h11);
	scale = (scaleX + scaleY) * 0.5;
	siril_log_message(_("scaleX:%*.3f\n"), 12, scaleX);
	siril_log_message(_("scaleY:%*.3f\n"), 12, scaleY);
	siril_log_message(_("scale:%*.3f\n"), 13, scale);

	/* Rotation */
	rotation = atan2(H.h01, H.h00) * 180 / M_PI;
	siril_log_message(_("rotation:%+*.3f deg\n"), 9, rotation);

	/* Translation */
	shift.x = -H.h02;
	shift.y = -H.h12;
	siril_log_message(_("dx:%+*.2f px\n"), 15, shift.x);
	siril_log_message(_("dy:%+*.2f px\n"), 15, shift.y);
	siril_log_message(_("FWHMx:%*.2f %s\n"), 12, FWHMx, units);
	siril_log_message(_("FWHMy:%*.2f %s\n"), 12, FWHMy, units);
}

