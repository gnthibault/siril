#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "algos/statistics.h"
#include "stacking.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "gui/progress_and_log.h"
#include "gui/utils.h"

static int compute_normalization(struct stacking_args *args);

/* normalization: reading all images and making stats on their background level.
 * That's very long if not cached. */
int do_normalization(struct stacking_args *args) {
	if (args->normalize == NO_NORM) return ST_OK;

	int nb_frames = args->nb_images_to_stack;
	int nb_layers = args->seq->nb_layers;

	args->coeff.offset = malloc(nb_layers * nb_frames * sizeof(double));
	args->coeff.mul = malloc(nb_layers * nb_frames * sizeof(double));
	args->coeff.scale = malloc(nb_layers * nb_frames * sizeof(double));

	if (!args->coeff.offset || !args->coeff.mul || !args->coeff.scale) {
		PRINT_ALLOC_ERR;
		args->retval = ST_ALLOC_ERROR;
		return args->retval;
	}

	if (compute_normalization(args)) {
		args->retval = ST_GENERIC_ERROR;
		return args->retval;
	}

	if (args->seq->needs_saving)	// if we had to compute new stats
		writeseqfile(args->seq);

	return ST_OK;
}

/* scale0, mul0 and offset0 are output arguments when i = ref_image, input arguments otherwise */
static int _compute_normalization_for_image(struct stacking_args *args, int i, int ref_image,
		double **poffset, double **pmul, double **pscale, normalization mode, double *scale0,
		double *mul0, double *offset0, gboolean multithread, int thread_id) {
	imstats *stat = NULL;
	gboolean fit_is_open = FALSE;
	fits fit = { 0 };


	for (int layer = 0; layer < args->seq->nb_layers; ++layer) {
		// try with no fit passed: fails if data is needed because data is not cached
		if (!(stat = statistics(args->seq, args->image_indices[i], NULL, layer, NULL, STATS_NORM, multithread))) {
			if (!(fit_is_open)) {
				// read frames as float, it's faster to compute stats
				if (seq_read_frame(args->seq, args->image_indices[i], &fit, TRUE, thread_id)) {
					return ST_SEQUENCE_ERROR;
				}
				fit_is_open = TRUE; // to avoid opening fit more than once if RGB
			}
			// retry with the fit to compute it
			if (!(stat = statistics(args->seq, args->image_indices[i], &fit, layer, NULL, STATS_NORM, multithread)))
				return ST_GENERIC_ERROR;
		}

		/* deal with the conversion of offset to the correct range depending
           on what was returned by statistics . This is done on each layer of each frame
		   just in case stats were not consistent in the seq file. */
		double conversionfactor = 1.0;
		if (!(args->seq->bitpix == FLOAT_IMG) && (stat->normValue == 1.)) {
			conversionfactor = (args->seq->bitpix ==  BYTE_IMG) ? UCHAR_MAX_DOUBLE : USHRT_MAX_DOUBLE;
		}

		double sc = stat->scale;
		double loc = stat->location * conversionfactor;

		switch (mode) {
		default:
		case ADDITIVE_SCALING:
			pscale[layer][i] = sc;
			if (i == ref_image)
				scale0[layer] = pscale[layer][i];
			pscale[layer][i] = (pscale[layer][i] == 0) ? 1 : scale0[layer] / pscale[layer][i];
			/* no break */
		case ADDITIVE:
			poffset[layer][i] = loc;
			if (i == ref_image)
				offset0[layer] = poffset[layer][i];
			poffset[layer][i] = pscale[layer][i] * poffset[layer][i] - offset0[layer];
			break;
		case MULTIPLICATIVE_SCALING:
			pscale[layer][i] = sc;
			if (i == ref_image)
				scale0[layer] = pscale[layer][i];
			pscale[layer][i] = (pscale[layer][i] == 0) ? 1 : scale0[layer]  / pscale[layer][i];
			/* no break */
		case MULTIPLICATIVE:
			pmul[layer][i] = loc;
			if (i == ref_image)
				mul0[layer] = pmul[layer][i];
			pmul[layer][i] = (pmul[layer][i] == 0) ? 1 : mul0[layer] / pmul[layer][i];
			break;
		}
	}
	if (fit_is_open && args->seq->type != SEQ_INTERNAL)
		clearfits(&fit);
	free_stats(stat);
	return ST_OK;
}

static int normalization_get_max_number_of_threads(sequence *seq) {
	int max_memory_MB = get_max_memory_in_MB();
	/* The normalization memory consumption, n is image size and m channel size.
	 * It uses IKSS computation in stats, which can be done only on float data.
	 * IKSS requires computing the MAD, which requires its own copy of the data.
	 * The stats are computed successively on each channel.
	 * For DATA_USHORT, we have: the image O(n), rewrite of the channel without
	 * zeros O(m), a copy to float for IKSS O(2m), a copy of that for the MAD in
	 * IKSS O(2m).
	 * For DATA_FLOAT, we have: the image O(n), rewrite without zeros O(m),
	 * used directly for IKSS and a copy for MAD O(m).
	 */
	guint64 memory_per_image = seq->rx * seq->ry;
	if (get_data_type(seq->bitpix) == DATA_FLOAT)
		memory_per_image *= (seq->nb_layers + 2) * sizeof(float);
	else memory_per_image *= (seq->nb_layers + 1) * sizeof(WORD) + 2 * sizeof(float);
	unsigned int memory_per_image_MB = memory_per_image / BYTES_IN_A_MB;

	fprintf(stdout, "Memory per image: %u MB. Max memory: %d MB\n", memory_per_image_MB, max_memory_MB);

	if (memory_per_image_MB > max_memory_MB) {
		siril_log_color_message(_("Your system does not have enough memory to normalize images for stacking operation (%d MB free for %d MB required)\n"), "red", max_memory_MB, memory_per_image_MB);
		return 0;
	}

	int nb_threads = memory_per_image_MB ? max_memory_MB / memory_per_image_MB : 1;
	if (nb_threads > com.max_thread)
		nb_threads = com.max_thread;
	siril_log_message(_("With the current memory and thread (%d) limits, up to %d thread(s) can be used for sequence normalization\n"), com.max_thread, nb_threads);
	return nb_threads;
}

static int compute_normalization(struct stacking_args *args) {
	int i, ref_image_filtred_idx = -1, retval = 0, cur_nb = 1;
	double *scale0, *mul0, *offset0;	// for reference frame
	char *tmpmsg;
	norm_coeff *coeff = &args->coeff;
	int nb_layers = args->seq->nb_layers;
	int nb_frames = args->nb_images_to_stack;

	for (int layer = 0; layer < nb_layers; ++layer) {
		coeff->poffset[layer] = coeff->offset + layer * nb_frames;
		coeff->pmul[layer] = coeff->mul + layer * nb_frames;
		coeff->pscale[layer] = coeff->scale + layer * nb_frames;
	}

	for (i = 0; i < args->nb_images_to_stack; ++i) {
		for (int layer = 0; layer < nb_layers; ++layer) {
			coeff->poffset[layer][i] = 0.0;
			coeff->pmul[layer][i] = 1.0;
			coeff->pscale[layer][i] = 1.0;
		}
	}

	if (args->normalize == NO_NORM)	// should never happen here
		return 0;

	scale0 = malloc(nb_layers * sizeof(double));
	mul0 = malloc(nb_layers * sizeof(double));
	offset0 = malloc(nb_layers * sizeof(double));

	tmpmsg = siril_log_message(_("Computing normalization...\n"));
	tmpmsg[strlen(tmpmsg) - 1] = '\0';
	set_progress_bar_data(tmpmsg, PROGRESS_RESET);

	// first, find the index of the ref image in the filtered image list
	ref_image_filtred_idx = find_refimage_in_indices(args->image_indices,
			args->nb_images_to_stack, args->ref_image);
	if (ref_image_filtred_idx == -1) {
		siril_log_color_message(_("The reference image is not in the selected set of images. "
				"Please choose another reference image.\n"), "red");
		free(scale0);
		free(mul0);
		free(offset0);
		return ST_GENERIC_ERROR;
	}

	// check memory first
	const char *error_msg = (_("Normalization failed."));
	int nb_threads = normalization_get_max_number_of_threads(args->seq);
	if (nb_threads == 0) {
		set_progress_bar_data(error_msg, PROGRESS_NONE);
		free(scale0);
		free(mul0);
		free(offset0);
		return ST_GENERIC_ERROR;
	}

	/* We empty the cache if needed (force to recompute) */
	if (args->force_norm)
		clear_stats(args->seq, args->reglayer);

	// compute for the first image to have scale0 mul0 and offset0 for each layer

	if (_compute_normalization_for_image(args, ref_image_filtred_idx,
			ref_image_filtred_idx, coeff->poffset, coeff->pmul, coeff->pscale,
			args->normalize, scale0, mul0, offset0, TRUE, -1)) {
		siril_log_color_message(_("%s Check image %d first.\n"), "red", error_msg,
				ref_image_filtred_idx + 1);
		set_progress_bar_data(error_msg, PROGRESS_NONE);
		free(scale0);
		free(mul0);
		free(offset0);
		return ST_GENERIC_ERROR;
	}


	set_progress_bar_data(NULL, 1.0 / (double)args->nb_images_to_stack);

#ifdef _OPENMP
#pragma omp parallel for num_threads(nb_threads) private(i) schedule(guided) \
	if (args->seq->type == SEQ_SER || ((args->seq->type == SEQ_REGULAR || args->seq->type == SEQ_FITSEQ) && fits_is_reentrant()))
#endif

	for (int i = 0; i < args->nb_images_to_stack; ++i) {
		if (!retval && i != ref_image_filtred_idx) {
			if (!get_thread_run()) {
				retval = 1;
				continue;
			}
			int thread_id = -1;
#ifdef _OPENMP
			thread_id = omp_get_thread_num();
#endif
			if (_compute_normalization_for_image(args, i, ref_image_filtred_idx,
						coeff->poffset, coeff->pmul, coeff->pscale,
						args->normalize, scale0, mul0, offset0,
						FALSE, thread_id)) {
				siril_log_color_message(_("%s Check image %d first.\n"), "red",
						error_msg, args->image_indices[i] + 1);
				set_progress_bar_data(error_msg, PROGRESS_NONE);
				retval = 1;
				continue;
			}
#ifdef _OPENMP
#pragma omp atomic
#endif
			cur_nb++;	// only used for progress bar
			set_progress_bar_data(NULL,
					(double)cur_nb / (double)args->nb_images_to_stack );
		}
	}

	set_progress_bar_data(NULL, PROGRESS_DONE);
	free(scale0);
	free(mul0);
	free(offset0);
	return retval;
}

