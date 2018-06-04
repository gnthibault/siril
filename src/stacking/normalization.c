#include <string.h>
#include "core/siril.h"
#include "algos/statistics.h"
#include "stacking.h"
#include "io/sequence.h"
#include "core/proto.h"
#include "gui/progress_and_log.h"
#include "gui/callbacks.h"

static int compute_normalization(struct stacking_args *args);

/* normalization: reading all images and making stats on their background level.
 * That's very long if not cached. */
int do_normalization(struct stacking_args *args) {
	if (args->normalize == NO_NORM) return 0;

	int nb_frames = args->nb_images_to_stack;

	args->coeff.offset = malloc(nb_frames * sizeof(double));
	args->coeff.mul = malloc(nb_frames * sizeof(double));
	args->coeff.scale = malloc(nb_frames * sizeof(double));
	if (!args->coeff.offset || !args->coeff.mul || !args->coeff.scale) {
		printf("allocation issue in stacking normalization\n");
		return -1;
	}

	if (compute_normalization(args))
		return -1;

	if (args->seq->needs_saving)	// if we had to compute new stats
		writeseqfile(args->seq);

	return 0;
}

/* scale0, mul0 and offset0 are output arguments when i = ref_image, input arguments otherwise */
static int _compute_normalization_for_image(struct stacking_args *args, int i, int ref_image,
		double *offset, double *mul, double *scale, normalization mode, double *scale0,
		double *mul0, double *offset0) {
	imstats *stat = NULL;
	int reglayer;

	reglayer = (args->reglayer == -1) ? 0 : args->reglayer;

	// try with no fit passed: fails if data is needed because data is not cached
	if (!(stat = statistics(args->seq, args->image_indices[i], NULL, reglayer, NULL, STATS_EXTRA))) {
		fits fit = { 0 };
		if (seq_read_frame(args->seq, args->image_indices[i], &fit)) {
			return 1;
		}
		// retry with the fit to compute it
		if (!(stat = statistics(args->seq, args->image_indices[i], &fit, reglayer, NULL, STATS_EXTRA)))
			return 1;
		if (args->seq->type != SEQ_INTERNAL)
			clearfits(&fit);
	}

	switch (mode) {
	default:
	case ADDITIVE_SCALING:
		scale[i] = stat->scale;
		if (i == ref_image)
			*scale0 = scale[ref_image];
		scale[i] = (scale[i] == 0) ? 1 : *scale0 / scale[i];
		/* no break */
	case ADDITIVE:
		offset[i] = stat->location;
		if (i == ref_image)
			*offset0 = offset[ref_image];
		offset[i] = scale[i] * offset[i] - *offset0;
		break;
	case MULTIPLICATIVE_SCALING:
		scale[i] = stat->scale;
		if (i == ref_image)
			*scale0 = scale[ref_image];
		scale[i] = (scale[i] == 0) ? 1 : *scale0 / scale[i];
		/* no break */
	case MULTIPLICATIVE:
		mul[i] = stat->location;
		if (i == ref_image)
			*mul0 = mul[ref_image];
		mul[i] = (mul[i] == 0) ? 1 : *mul0 / mul[i];
		break;
	}

	free_stats(stat);
	return 0;
}

static int compute_normalization(struct stacking_args *args) {
	int i, ref_image, ref_image_filtred_idx = -1, retval = 0,
			cur_nb = 1;
	double scale0, mul0, offset0;	// for reference frame
	char *tmpmsg;
	norm_coeff *coeff = &args->coeff;

	for (i = 0; i < args->nb_images_to_stack; i++) {
		coeff->offset[i] = 0.0;
		coeff->mul[i] = 1.0;
		coeff->scale[i] = 1.0;
	}
	scale0 = mul0 = offset0 = 0.0;
	if (args->normalize == NO_NORM)	// should never happen here
		return 0;

	tmpmsg = siril_log_message(_("Computing normalization...\n"));
	tmpmsg[strlen(tmpmsg) - 1] = '\0';
	set_progress_bar_data(tmpmsg, PROGRESS_RESET);

	// first, find the index of the ref image in the filtered image list
	ref_image = sequence_find_refimage(args->seq);
	for (i = 0; i < args->nb_images_to_stack; i++)
		if (args->image_indices[i] == ref_image) {
			ref_image_filtred_idx = i;
			break;
		}
	if (ref_image_filtred_idx == -1) {
		siril_log_color_message(_("The reference image is not in the selected set of images. "
				"Please choose another reference image.\n"), "red");
		siril_log_color_message(_("Normalisation skipped.\n"), "red");
		return 1;
	}

	/* We empty the cache if needed (force to recompute) */
	if (args->force_norm)
		clear_stats(args->seq, args->reglayer);

	// compute for the first image to have scale0 mul0 and offset0
	if (_compute_normalization_for_image(args,
				ref_image_filtred_idx, ref_image_filtred_idx,
				coeff->offset, coeff->mul, coeff->scale,
				args->normalize, &scale0, &mul0, &offset0)) {
		set_progress_bar_data(_("Normalization failed."), PROGRESS_NONE);
		return 1;
	}

	set_progress_bar_data(NULL, 1.0 / (double)args->nb_images_to_stack);

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static) if (args->seq->type == SEQ_SER || fits_is_reentrant())
#endif
	for (i = 0; i < args->nb_images_to_stack; ++i) {
		if (!retval && i != ref_image_filtred_idx) {
			if (!get_thread_run()) {
				retval = 1;
				continue;
			}
			if (_compute_normalization_for_image(args, i, ref_image_filtred_idx,
						coeff->offset, coeff->mul, coeff->scale,
						args->normalize, &scale0, &mul0, &offset0)) {
				retval = 1;
				continue;
			}
#ifdef _OPENMP
#pragma omp atomic
#endif
			cur_nb++;	// only used for progress bar
			set_progress_bar_data(NULL,
					(double)cur_nb / ((double)args->nb_images_to_stack));
		}
	}
	set_progress_bar_data(NULL, PROGRESS_DONE);
	return retval;
}

