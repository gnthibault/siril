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

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/stat.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics_ushort.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef MAC_INTEGRATION
#include "gtkmacintegration/gtkosxapplication.h"
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "gui/callbacks.h"
#include "io/single_image.h"
#include "registration/registration.h"
#include "stacking/stacking.h"
#include "algos/PSF.h"
#include "gui/PSF_list.h"
#include "gui/histogram.h"	// update_gfit_histogram_if_needed();
#include "io/ser.h"

#undef STACK_DEBUG

static struct stacking_args stackparam = {	// parameters passed to stacking
		NULL, NULL, NULL, -1.0, 0, NULL, { '\0' }, NULL, FALSE, { 0, 0 }, -1, 0, { 0, 0 }, NO_REJEC, NO_NORM, FALSE
};

static stack_method stacking_methods[] = {
	stack_summing, stack_mean_with_rejection, stack_median, stack_addmax, stack_addmin
};

static gboolean end_stacking(gpointer p);


/* pool of memory blocks for parallel processing */
struct _data_block {
	WORD **pix;	// buffer for a block on all images
	WORD *tmp;	// the actual single buffer for pix
	WORD *stack;	// the reordered stack for one pixel in all images
	int *rejected;  // 0 if pixel ok, 1 or -1 if rejected
};

void initialize_stacking_methods() {
	GtkComboBoxText *stackcombo, *rejectioncombo;

	stackcombo = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(builder, "comboboxstack_methods"));
	rejectioncombo = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(builder, "comborejection"));
	gtk_combo_box_set_active(GTK_COMBO_BOX(stackcombo), com.stack.method);
	gtk_combo_box_set_active(GTK_COMBO_BOX(rejectioncombo), com.stack.rej_method);
}

/* scale0, mul0 and offset0 are output arguments when i = ref_image, input arguments otherwise */
static int _compute_normalization_for_image(struct stacking_args *args, int i, int ref_image,
		double *offset, double *mul, double *scale, normalization mode, double *scale0,
		double *mul0, double *offset0) {
	imstats *stat = NULL;

	if (!(stat = seq_get_imstats(args->seq, args->image_indices[i], NULL, STATS_EXTRA))) {
		fits fit;
		memset(&fit, 0, sizeof(fits));
		if (seq_read_frame(args->seq, args->image_indices[i], &fit)) {
			return 1;
		}
		stat = seq_get_imstats(args->seq, args->image_indices[i], &fit, STATS_EXTRA);
		if (args->seq->type != SEQ_INTERNAL)
			clearfits(&fit);
	}

	switch (mode) {
	default:
	case ADDITIVE_SCALING:
		scale[i] = stat->scale;
		if (i == ref_image)
			*scale0 = scale[ref_image];
		scale[i] = *scale0 / scale[i];
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
		scale[i] = *scale0 / scale[i];
		/* no break */
	case MULTIPLICATIVE:
		mul[i] = stat->location;
		if (i == ref_image)
			*mul0 = mul[ref_image];
		mul[i] = *mul0 / mul[i];
		break;
	}
	return 0;
}

int compute_normalization(struct stacking_args *args, norm_coeff *coeff, normalization mode) {
	int i, ref_image, retval = 0, cur_nb = 1;
	double scale0, mul0, offset0;	// for reference frame
	char *tmpmsg;

	for (i = 0; i < args->nb_images_to_stack; i++) {
		coeff->offset[i] = 0.0;
		coeff->mul[i] = 1.0;
		coeff->scale[i] = 1.0;
	}
	scale0 = mul0 = offset0 = 0.0;
	if (mode == 0)
		return 0;

	tmpmsg = siril_log_message(_("Computing normalization...\n"));
	tmpmsg[strlen(tmpmsg) - 1] = '\0';
	set_progress_bar_data(tmpmsg, PROGRESS_RESET);

	if (args->seq->reference_image == -1)
		ref_image = 0;
	else ref_image = args->seq->reference_image;

	/* We empty the cache if needed (force to recompute) */
	if (args->force_norm) {
		for (i = 0; i < args->seq->number; i++) {
			if (args->seq->imgparam && args->seq->imgparam[i].stats) {
				free(args->seq->imgparam[i].stats);
				args->seq->imgparam[i].stats = NULL;
			}
		}
	}

	// compute for the first image to have scale0 mul0 and offset0
	if (_compute_normalization_for_image(args, ref_image, ref_image, coeff->offset, coeff->mul, coeff->scale, mode,
			&scale0, &mul0, &offset0)) {
		set_progress_bar_data(_("Normalization failed."), PROGRESS_NONE);
		return 1;
	}

	set_progress_bar_data(NULL, 1.0 / (double)args->nb_images_to_stack);

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static) if (args->seq->type == SEQ_SER || fits_is_reentrant())
#endif
	for (i = 0; i < args->nb_images_to_stack; ++i) {
		if (!retval && i != ref_image) {
			if (!get_thread_run()) {
				retval = 1;
				continue;
			}
			if (_compute_normalization_for_image(args, i, ref_image, coeff->offset, coeff->mul, coeff->scale,
					mode, &scale0, &mul0, &offset0)) {
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

/** STACK method **
 * This method takes several images and create a new being the sum of all
 * others (normalized to the maximum value of unsigned SHORT).
 */
int stack_summing(struct stacking_args *args) {
	int x, y, nx, ny, i, ii, j, shiftx, shifty, layer, reglayer;
	unsigned long *somme[3], *from, maxim = 0;
	WORD *to;
	double ratio;
	double exposure=0.0;
	unsigned int nbdata = 0;
	char filename[256];
	int retval = 0;
	int nb_frames, cur_nb = 0;
	fits *fit = &wfit[0];
	char *tmpmsg;
	memset(fit, 0, sizeof(fits));

	/* should be pre-computed to display it in the stacking tab */
	nb_frames = args->nb_images_to_stack;
	reglayer = get_registration_layer();

	if (nb_frames <= 1) {
		siril_log_message(_("No frame selected for stacking (select at least 2). Aborting.\n"));
		return -1;
	}

	somme[0] = NULL;
	assert(nb_frames <= args->seq->number);
	set_progress_bar_data(NULL, PROGRESS_RESET);

	for (j=0; j<args->seq->number; ++j){
		if (!get_thread_run()) {
			retval = -1;
			goto free_and_reset_progress_bar;
		}
		if (!args->filtering_criterion(args->seq, j, args->filtering_parameter)) {
			fprintf(stdout, "image %d is excluded from stacking\n", j);
			continue;
		}
		if (!seq_get_image_filename(args->seq, j, filename)) {
			retval = -1;
			goto free_and_reset_progress_bar;
		}
		tmpmsg = strdup(_("Processing image "));
		tmpmsg = str_append(&tmpmsg, filename);
		set_progress_bar_data(tmpmsg, (double)cur_nb/((double)nb_frames+1.));
		free(tmpmsg);

		cur_nb++;	// only used for progress bar

		if (seq_read_frame(args->seq, j, fit)) {
			siril_log_message(_("Stacking: could not read frame, aborting\n"));
			retval = -3;
			goto free_and_reset_progress_bar;
		}

		if (args->seq->nb_layers == -1) {
			/* sequence has not been opened before, this is set in set_seq.
			 * It happens with the stackall command that stacks a
			 * sequence right after readseqfile.
			 */
			args->seq->rx = fit->rx; args->seq->ry = fit->ry;
			args->seq->nb_layers = fit->naxes[2];
		}
		assert(args->seq->nb_layers == 1 || args->seq->nb_layers == 3);
		assert(fit->naxes[2] == args->seq->nb_layers);

		/* first loaded image: init data structures for stacking */
		if (!nbdata) {
			nbdata = fit->ry * fit->rx;
			somme[0] = calloc(nbdata, sizeof(unsigned long)*fit->naxes[2]);
			if (somme[0] == NULL){
				printf("Stacking: memory allocation failure\n");
				retval = -2;
				goto free_and_reset_progress_bar;
			}
			if(args->seq->nb_layers == 3){
				somme[1] = somme[0] + nbdata;	// index of green layer in somme[0]
				somme[2] = somme[0] + nbdata*2;	// index of blue layer in somme[0]
			}
			//~ "stacking operation\n");
		} else if (fit->ry * fit->rx != nbdata) {
			siril_log_message(_("Stacking: image in sequence doesn't has the same dimensions\n"));
			retval = -3;
			goto free_and_reset_progress_bar;
		}

		update_used_memory();

		/* load registration data for current image */
		if(reglayer != -1 && args->seq->regparam[reglayer]) {
			shiftx = args->seq->regparam[reglayer][j].shiftx;
			shifty = args->seq->regparam[reglayer][j].shifty;
		} else {
			shiftx = 0;
			shifty = 0;
		}
#ifdef STACK_DEBUG
		printf("Stack image %d with shift x=%d y=%d\n", j, shiftx, shifty);
#endif

		/* Summing the exposure */
		exposure += fit->exposure;

		/* stack current image */
		i=0;	// index in somme[0]
		for (y=0; y < fit->ry; ++y){
			for (x=0; x < fit->rx; ++x){
				nx = x - shiftx;
				ny = y - shifty;
				//printf("x=%d y=%d sx=%d sy=%d i=%d ii=%d\n",x,y,shiftx,shifty,i,ii);
				if (nx >= 0 && nx < fit->rx && ny >= 0 && ny < fit->ry) {
					ii = ny * fit->rx + nx;		// index in somme[0] too
					//printf("shiftx=%d shifty=%d i=%d ii=%d\n",shiftx,shifty,i,ii);
					if (ii > 0 && ii < fit->rx * fit->ry){
						for(layer=0; layer<args->seq->nb_layers; ++layer){
							WORD current_pixel = fit->pdata[layer][ii];
							somme[layer][i] += current_pixel;
							if (somme[layer][i] > maxim){
								maxim = somme[layer][i];
							}
						}
					}
				}
				++i;
			}
		}
	}
	set_progress_bar_data(_("Finalizing stacking..."), (double)nb_frames/((double)nb_frames + 1.));

	copyfits(fit, &gfit, CP_ALLOC|CP_FORMAT, 0);
	gfit.hi = round_to_WORD(maxim);
	gfit.bitpix = USHORT_IMG;
	gfit.exposure = exposure;

	if (maxim > USHRT_MAX)
		ratio = USHRT_MAX_DOUBLE / (double)maxim;
	else	ratio = 1.0;

	if (somme[0]) {
		assert(args->seq->nb_layers == 1 || args->seq->nb_layers == 3);
		for(layer=0; layer<args->seq->nb_layers; ++layer){
			from = somme[layer];
			to = gfit.pdata[layer];
			for (y=0; y < fit->ry * fit->rx; ++y) {
				if (ratio == 1.0)
					*to++ = round_to_WORD(*from++);
				else	*to++ = round_to_WORD((double)(*from++) * ratio);
			}
		}
	}

free_and_reset_progress_bar:
	if (somme[0]) free(somme[0]);
	if (retval) {
		set_progress_bar_data(_("Stacking failed. Check the log."), PROGRESS_RESET);
		siril_log_message(_("Stacking failed.\n"));
	} else {
		set_progress_bar_data(_("Stacking complete."), PROGRESS_DONE);
	}
	update_used_memory();
	return retval;
}

/******************************* MEDIAN STACKING ******************************
 * This is a bit special as median stacking requires all images to be in memory
 * So we dont use the generic readfits but directly the cfitsio routines, and
 * allocates as many pix tables as needed.
 * ****************************************************************************/
int stack_median(struct stacking_args *args) {
	int nb_frames;		/* number of frames actually used */
	int status;		/* CFITSIO status value MUST be initialized to zero for EACH call */
	int bitpix;
	int naxis, oldnaxis = -1, cur_nb = 0;
	long npixels_in_block, nbdata;
	long naxes[3], oldnaxes[3];
	int i;
	double exposure = 0.0;
	char filename[256], msg[256];
	int retval = 0;
	struct _data_block *data_pool = NULL;
	int pool_size = 1;
	fits *fit = &wfit[0];
	norm_coeff coeff;

	nb_frames = args->nb_images_to_stack;

	if (args->seq->type != SEQ_REGULAR && args->seq->type != SEQ_SER) {
		char *msg = siril_log_message(_("Median stacking is only supported for FITS images and SER sequences.\n"));
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		return -1;
	}
	if (nb_frames < 2) {
		siril_log_message(_("Select at least two frames for stacking. Aborting.\n"));
		return -1;
	}

	assert(nb_frames <= args->seq->number);
	set_progress_bar_data(NULL, PROGRESS_RESET);

	/* allocate data structures */
	oldnaxes[0] = oldnaxes[1] = oldnaxes[2] = 0;	// fix compiler warning
	naxes[0] = naxes[1] = 0; naxes[2] = 1;

	/* first loop: open all fits files and check they are of same size */
	if (args->seq->type == SEQ_REGULAR) {
		for (i=0; i<nb_frames; ++i) {
			int image_index = args->image_indices[i];	// image index in sequence
			if (!get_thread_run()) {
				retval = -1;
				goto free_and_close;
			}
			if (!fit_sequence_get_image_filename(args->seq, image_index, filename, TRUE))
				continue;

			snprintf(msg, 255, _("Median stack: opening image %s"), filename);
			msg[255] = '\0';
			set_progress_bar_data(msg, PROGRESS_NONE);

			/* open input images */
			if (seq_open_image(args->seq, image_index)) {
				retval = -1;
				goto free_and_close;
			}

			/* here we use the internal data of sequences, it's quite ugly, we should
			 * consider moving these tests in seq_open_image() or wrapping them in a
			 * sequence function */
			status = 0;
			fits_get_img_param(args->seq->fptr[image_index], 3, &bitpix, &naxis, naxes, &status);
			if (status) {
				fits_report_error(stderr, status); /* print error message */
				retval = status;
				goto free_and_close;
			}
			if (naxis > 3) {
				siril_log_message(_("Median stack error: images with > 3 dimensions "
						"are not supported\n"));
				retval = -1;
				goto free_and_close;
			}

			if(oldnaxis > 0) {
				if(naxis != oldnaxis ||
						oldnaxes[0] != naxes[0] ||
						oldnaxes[1] != naxes[1] ||
						oldnaxes[2] != naxes[2]) {
					siril_log_message(_("Median stack error: input images have "
							"different sizes\n"));
					retval = -2;
					goto free_and_close;
				}
			} else {
				oldnaxis = naxis;
				oldnaxes[0] = naxes[0];
				oldnaxes[1] = naxes[1];
				oldnaxes[2] = naxes[2];
			}

			/* exposure summing */
			double tmp;
			status = 0;
			/* and here we should provide a opened_fits_read_key for example */
			fits_read_key (args->seq->fptr[image_index], TDOUBLE, "EXPTIME", &tmp, NULL, &status);
			if (status || tmp <= 0.0) {
				status = 0;
				fits_read_key (args->seq->fptr[image_index], TDOUBLE, "EXPOSURE", &tmp, NULL, &status);
			}
			if (!status)
				exposure += tmp;
		}
		update_used_memory();
	}

	coeff.offset = malloc(nb_frames * sizeof(double));
	coeff.mul = malloc(nb_frames * sizeof(double));
	coeff.scale = malloc(nb_frames * sizeof(double));
	if (!coeff.offset || !coeff.mul || !coeff.scale) {
		printf("allocation issue in stacking normalization\n");
		retval = -1;
		goto free_and_close;
	}

	if (naxes[2] == 0)
		naxes[2] = 1;
	assert(naxes[2] <= 3);
	if (args->seq->type == SEQ_SER) {
		assert(args->seq->ser_file);
		naxes[0] = args->seq->ser_file->image_width;
		naxes[1] = args->seq->ser_file->image_height;
		ser_color type_ser = args->seq->ser_file->color_id;
		if (!com.debayer.open_debayer && type_ser != SER_RGB && type_ser != SER_BGR)
			type_ser = SER_MONO;
		naxes[2] = type_ser == SER_MONO ? 1 : 3;
		naxis = type_ser == SER_MONO ? 2 : 3;
		/* case of Super Pixel not handled yet */
		if (com.debayer.bayer_inter == BAYER_SUPER_PIXEL) {
			siril_log_message(_("Super-pixel is not handled yet for on the fly SER stacking\n")); /* TODO */
			retval = -1;
			goto free_and_close;
		}
	}
	if (naxes[0] == 0) {
		// no image has been loaded
		siril_log_message(_("Median stack error: uninitialized sequence\n"));
		retval = -2;
		goto free_and_close;
	}
	fprintf(stdout, "image size: %ldx%ld, %ld layers\n", naxes[0], naxes[1], naxes[2]);

	/* normalization: reading all images and making stats on their background level.
	 * That's very long if not cached. */
	if (compute_normalization(args, &coeff, args->normalize)) {
		retval = -1;
		goto free_and_close;
	}
	if (args->seq->needs_saving)	// if we had to compute new stats
		writeseqfile(args->seq);

	/* initialize result image */
	nbdata = naxes[0] * naxes[1];
	memset(fit, 0, sizeof(fits));
	fit->data = malloc(nbdata * naxes[2] * sizeof(WORD));
	if (!fit->data) {
		fprintf(stderr, "Memory allocation error for result\n");
		retval = -1;
		goto free_and_close;
	}
	fit->bitpix = USHORT_IMG;
	fit->naxes[0] = naxes[0];
	fit->naxes[1] = naxes[1];
	fit->naxes[2] = naxes[2];
	fit->rx = naxes[0];
	fit->ry = naxes[1];
	fit->naxis = naxis;
	fit->maxi = 0;
	if(fit->naxis == 3) {
		fit->pdata[RLAYER]=fit->data;
		fit->pdata[GLAYER]=fit->data + nbdata;
		fit->pdata[BLAYER]=fit->data + nbdata * 2;
	} else {
		fit->pdata[RLAYER]=fit->data;
		fit->pdata[GLAYER]=fit->data;
		fit->pdata[BLAYER]=fit->data;
	}
	update_used_memory();

	/* Define some useful constants */
	double total = (double)(naxes[2] * naxes[1] + 2);	// only used for progress bar

	int nb_threads;
#ifdef _OPENMP
	nb_threads = com.max_thread;
	if (args->seq->type == SEQ_REGULAR && fits_is_reentrant()) {
		fprintf(stdout, "cfitsio was compiled with multi-thread support,"
				" stacking will be executed by several cores\n");
	}
	if (args->seq->type == SEQ_REGULAR && !fits_is_reentrant()) {
		nb_threads = 1;
		fprintf(stdout, "cfitsio was compiled without multi-thread support,"
				" stacking will be executed on only one core\n");
		siril_log_message(_("Your version of cfitsio does not support multi-threading\n"));
	}
#else
	nb_threads = 1;
#endif

	int nb_channels = naxes[2];
	if (sequence_is_rgb(args->seq) && nb_channels != 3) {
		siril_log_message(_("Processing the sequence as RGB\n"));
		nb_channels = 3;
	}

	int size_of_stacks = args->max_number_of_rows / nb_threads; 
	if (size_of_stacks == 0)
		size_of_stacks = 1;
	/* Note: this size of stacks based on the max memory configured doesn't take into
	 * account memory for demosaicing if it applies.
	 * Now we compute the total number of "stacks" which are the independent areas where
	 * the stacking will occur. This will then be used to create the image areas. */
	long nb_parallel_stacks;
	int remainder;
	if (naxes[1] / size_of_stacks < 4) {
		/* We have enough RAM to process each channel with 4 threads.
		 * We should cut images at least in 4 on one channel to use enough threads,
		 * and if only one is available, it will use much less RAM for a small time overhead.
		 * Also, for slow data access like rotating drives or on-the-fly debayer,
		 * it feels more responsive this way.
		 */
		nb_parallel_stacks = 4 * nb_channels;
		size_of_stacks = naxes[1] / 4;
		remainder = naxes[1] % 4;
	} else {
		/* We don't have enough RAM to process a channel with all available threads */
		nb_parallel_stacks = naxes[1] * nb_channels / size_of_stacks;
		if (nb_parallel_stacks % nb_channels != 0
				|| (naxes[1] * nb_channels) % size_of_stacks != 0) {
			/* we need to take into account the fact that the stacks are computed for
			 * each channel, not for the total number of pixels. So it needs to be
			 * a factor of the number of channels.
			 */
			nb_parallel_stacks += nb_channels - (nb_parallel_stacks % nb_channels);
			size_of_stacks = naxes[1] * nb_channels / nb_parallel_stacks;
		}
		remainder = naxes[1] - (nb_parallel_stacks / nb_channels * size_of_stacks);
	}
	siril_log_message(_("We have %d parallel blocks of size %d (+%d) for stacking.\n"),
			nb_parallel_stacks, size_of_stacks, remainder);
	struct image_block {
		unsigned long channel, start_row, end_row, height;
	};
	long largest_block_height = 0;
	struct image_block *blocks = malloc(nb_parallel_stacks * sizeof(struct image_block));
	{
		long channel = 0, row = 0, end, j = 0;
		do {
			if (j >= nb_parallel_stacks) {
				siril_log_message(_("A bug has been found. "
						"Unable to split the image area into the correct processing blocks.\n"));
				retval = -1;
				goto free_and_close;
			}

			blocks[j].channel = channel;
			blocks[j].start_row = row;
			end = row + size_of_stacks - 1; 
			if (remainder > 0) {
				// just add one pixel from the remainder to the first blocks to
				// avoid having all of them in the last block
				end++;
				remainder--;
			}
			if (end >= naxes[1] - 1 ||	// end of the line
					(naxes[1] - end < size_of_stacks / 10)) { // not far from it
				end = naxes[1] - 1;
				row = 0;
				channel++;
				remainder = naxes[1] - (nb_parallel_stacks / nb_channels * size_of_stacks);
			} else {
				row = end + 1;
			}
			blocks[j].end_row = end;
			blocks[j].height = blocks[j].end_row - blocks[j].start_row + 1;
			if (largest_block_height < blocks[j].height) {
				largest_block_height = blocks[j].height;
			}
			fprintf(stdout, "Block %ld: channel %lu, from %lu to %lu (h = %lu)\n",
					j, blocks[j].channel, blocks[j].start_row,
					blocks[j].end_row, blocks[j].height);
			j++;
	
		} while (channel < nb_channels) ;
	}

	/* Allocate the buffers.
	 * We allocate as many as the number of threads, each thread will pick one of the buffers.
	 * Buffers are allocated to the largest block size calculated above.
	 */
#ifdef _OPENMP
	pool_size = nb_threads;
	assert(pool_size > 0);
#endif
	npixels_in_block = largest_block_height * naxes[0];
	fprintf(stdout, "allocating data for %d threads (each %'lu MB)\n", pool_size,
			(unsigned long) (nb_frames * npixels_in_block * sizeof(WORD)) / 1048576UL);
	data_pool = malloc(pool_size * sizeof(struct _data_block));
	for (i = 0; i < pool_size; i++) {
		int j;
		data_pool[i].pix = calloc(nb_frames, sizeof(WORD *));
		data_pool[i].tmp = calloc(nb_frames, npixels_in_block * sizeof(WORD));
		data_pool[i].stack = calloc(nb_frames, sizeof(WORD));
		if (!data_pool[i].pix || !data_pool[i].tmp || !data_pool[i].stack) {
			fprintf(stderr, "Memory allocation error on pix.\n");
			fprintf(stderr, "CHANGE MEMORY SETTINGS if stacking takes too much.\n");
			retval = -1;
			goto free_and_close;
		}
		for (j=0; j<nb_frames; ++j) {
			data_pool[i].pix[j] = data_pool[i].tmp + j * npixels_in_block;
		}
	}
	update_used_memory();

	siril_log_message(_("Starting stacking...\n"));
	set_progress_bar_data(_("Median stacking in progress..."), PROGRESS_RESET);

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static) if (args->seq->type == SEQ_SER || fits_is_reentrant())
#endif
	for (i = 0; i < nb_parallel_stacks; i++)
	{
		/**** Step 1: get allocated memory for the current thread ****/
		struct image_block *my_block = blocks+i;
		struct _data_block *data;
		int data_idx = 0, frame;
		long x, y;

		if (!get_thread_run()) retval = -1;
		if (retval) continue;
#ifdef _OPENMP
		data_idx = omp_get_thread_num();
#endif
		assert(data_idx < pool_size);
		//fprintf(stdout, "thread %d working on block %d gets data\n", data_idx, i);
		data = &data_pool[data_idx];

		/**** Step 2: load image data for the corresponding image block ****/
		/* area in C coordinates, starting with 0, not cfitsio coordinates. */
		rectangle area = {0, my_block->start_row, naxes[0], my_block->height};
		/* Read the block from all images, store them in pix[image] */
		for (frame = 0; frame < nb_frames; ++frame){
			if (!get_thread_run()) {
				retval = -1;
				break;
			}

			// reading pixels from current frame
			int success = seq_opened_read_region(args->seq, my_block->channel,
					args->image_indices[frame], data->pix[frame], &area);

			if (success < 0)
				retval = -1;

			if (retval) {
				siril_log_message(_("Error reading one of the image areas\n"));
				break;
			}
		}
		if (retval) continue;

		/**** Step 3: iterate over the y and x of the image block and stack ****/
		for (y = 0; y < my_block->height; y++)
		{
			/* index of the pixel in the result image
			 * we read line y, but we need to store it at
			 * ry - y - 1 to not have the image mirrored. */
			int pixel_idx = (naxes[1] - (my_block->start_row + y) - 1) * naxes[0]; 
			if (retval) break;

			// update progress bar
#ifdef _OPENMP
#pragma omp atomic
#endif
			cur_nb++;

			if (!get_thread_run()) {
				retval = -1;
				break;
			}
			set_progress_bar_data(NULL, (double)cur_nb/total);

			for (x = 0; x < naxes[0]; ++x){
				int ii;
				/* copy all images pixel values in the same row array `stack'
				 * to optimize caching and improve readability */
				for (ii=0; ii<nb_frames; ++ii) {
					double tmp = (double)data->pix[ii][y*naxes[0]+x] * coeff.scale[ii];
					switch (args->normalize) {
						default:
						case NO_NORM:		// no normalization (scale[ii] = 1, offset[ii] = 0, mul[ii] = 1)
						case ADDITIVE:		// additive (scale[ii] = 1, mul[ii] = 1)
						case ADDITIVE_SCALING:		// additive + scale (mul[ii] = 1)
							data->stack[ii] = round_to_WORD(tmp - coeff.offset[ii]);
							break;
						case MULTIPLICATIVE:		// multiplicative  (scale[ii] = 1, offset[ii] = 0)
						case MULTIPLICATIVE_SCALING:		// multiplicative + scale (offset[ii] = 0)
							data->stack[ii] = round_to_WORD(tmp * coeff.mul[ii]);
							break;
					}
				}
				quicksort_s(data->stack, nb_frames);
				fit->pdata[my_block->channel][pixel_idx] =
						gsl_stats_ushort_median_from_sorted_data(data->stack, 1, nb_frames);
				pixel_idx++;
			}
		}
	} /* end of loop over parallel stacks */

	if (retval)
		goto free_and_close;

	set_progress_bar_data(_("Finalizing stacking..."), PROGRESS_NONE);
	/* copy result to gfit if success */
	copyfits(fit, &gfit, CP_FORMAT, 0);
	if (gfit.data) free(gfit.data);
	gfit.data = fit->data;
	gfit.exposure = exposure;
	memcpy(gfit.pdata, fit->pdata, 3*sizeof(WORD *));

	fit->data = NULL;
	memset(fit->pdata, 0, 3*sizeof(WORD *));

free_and_close:
	fprintf(stdout, "free and close (%d)\n", retval);
	for (i=0; i<nb_frames; ++i) {
		seq_close_image(args->seq, args->image_indices[i]);
	}

	if (data_pool) {
		for (i=0; i<pool_size; i++) {
			if (data_pool[i].stack) free(data_pool[i].stack);
			if (data_pool[i].pix) free(data_pool[i].pix);
			if (data_pool[i].tmp) free(data_pool[i].tmp);
		}
		free(data_pool);
	}
	free(coeff.offset);
	free(coeff.mul);
	free(coeff.scale);
	if (retval) {
		/* if retval is set, gfit has not been modified */
		if (fit->data) free(fit->data);
		set_progress_bar_data(_("Median stacking failed. Check the log."), PROGRESS_RESET);
		siril_log_message(_("Stacking failed.\n"));
	} else {
		set_progress_bar_data(_("Median stacking complete."), PROGRESS_DONE);
		siril_log_message(_("Median stacking complete. %d have been stacked.\n"), nb_frames);
	}
	update_used_memory();
	return retval;
}


/** Addmax STACK method **
 * This method is very close to the summing one instead that the result
 * takes only the pixel if it is brighter than the previous one at the
 * same coordinates.
 */
int stack_addmax(struct stacking_args *args) {
	int x, y, nx, ny, i, ii, j, shiftx, shifty, layer, reglayer;
	WORD *final_pixel[3], *from, maxim = 0;
	WORD *to;
	double exposure=0.0;
	unsigned int nbdata = 0;
	char filename[256];
	int retval = 0;
	int nb_frames, cur_nb = 0;
	fits *fit = &wfit[0];
	char *tmpmsg;
	memset(fit, 0, sizeof(fits));

	/* should be pre-computed to display it in the stacking tab */
	nb_frames = args->nb_images_to_stack;
	reglayer = get_registration_layer();

	if (nb_frames <= 1) {
		siril_log_message(_("No frame selected for stacking (select at least 2). Aborting.\n"));
		return -1;
	}

	final_pixel[0] = NULL;
	assert(args->seq->nb_layers == 1 || args->seq->nb_layers == 3);
	assert(nb_frames <= args->seq->number);

	for (j=0; j<args->seq->number; ++j){
		if (!get_thread_run()) {
			retval = -1;
			goto free_and_reset_progress_bar;
		}
		if (!args->filtering_criterion(args->seq, j, args->filtering_parameter)) {
			fprintf(stdout, "image %d is excluded from stacking\n", j);
			continue;
		}
		if (!seq_get_image_filename(args->seq, j, filename)) {
			retval = -1;
			goto free_and_reset_progress_bar;
		}
		tmpmsg = strdup(_("Processing image "));
		tmpmsg = str_append(&tmpmsg, filename);
		set_progress_bar_data(tmpmsg, (double)cur_nb/((double)nb_frames+1.));
		free(tmpmsg);

		cur_nb++;	// only used for progress bar

		if (seq_read_frame(args->seq, j, fit)) {
			siril_log_message(_("Stacking: could not read frame, aborting\n"));
			retval = -3;
			goto free_and_reset_progress_bar;
		}

		assert(args->seq->nb_layers == 1 || args->seq->nb_layers == 3);
		assert(fit->naxes[2] == args->seq->nb_layers);

		/* first loaded image: init data structures for stacking */
		if (!nbdata) {
			nbdata = fit->ry * fit->rx;
			final_pixel[0] = calloc(nbdata, sizeof(WORD)*fit->naxes[2]);
			if (final_pixel[0] == NULL){
				printf("Stacking: memory allocation failure\n");
				retval = -2;
				goto free_and_reset_progress_bar;
			}
			if(args->seq->nb_layers == 3){
				final_pixel[1] = final_pixel[0] + nbdata;	// index of green layer in final_pixel[0]
				final_pixel[2] = final_pixel[0] + nbdata*2;	// index of blue layer in final_pixel[0]
			}
			//~ siril_log_message("Stacking: successfully allocated memory for "
					//~ "stacking operation\n");
		} else if (fit->ry * fit->rx != nbdata) {
			siril_log_message(_("Stacking: image in sequence doesn't has the same dimensions\n"));
			retval = -3;
			goto free_and_reset_progress_bar;
		}

		update_used_memory();

		/* load registration data for current image */
		if(reglayer != -1 && args->seq->regparam[reglayer]) {
			shiftx = args->seq->regparam[reglayer][j].shiftx;
			shifty = args->seq->regparam[reglayer][j].shifty;
		} else {
			shiftx = 0;
			shifty = 0;
		}
#ifdef STACK_DEBUG
		printf("stack image %d with shift x=%d y=%d\n", j, shiftx, shifty);
#endif

		/* Summing the exposure */
		exposure += fit->exposure;

		/* stack current image */
		i=0;	// index in final_pixel[0]
		for (y=0; y < fit->ry; ++y){
			for (x=0; x < fit->rx; ++x){
				nx = x - shiftx;
				ny = y - shifty;
				//printf("x=%d y=%d sx=%d sy=%d i=%d ii=%d\n",x,y,shiftx,shifty,i,ii);
				if (nx >= 0 && nx < fit->rx && ny >= 0 && ny < fit->ry) {
					ii = ny * fit->rx + nx;		// index in final_pixel[0] too
					//printf("shiftx=%d shifty=%d i=%d ii=%d\n",shiftx,shifty,i,ii);
					if (ii > 0 && ii < fit->rx * fit->ry){
						for(layer=0; layer<args->seq->nb_layers; ++layer){
							WORD current_pixel = fit->pdata[layer][ii];
							if (current_pixel > final_pixel[layer][i])	// we take the brighter pixel
								final_pixel[layer][i] = current_pixel;
							if (final_pixel[layer][i] > maxim){
								maxim = final_pixel[layer][i];
							}
						}
					}
				}
				++i;
			}
		}
	}
	if (!get_thread_run()) {
		retval = -1;
		goto free_and_reset_progress_bar;
	}
	set_progress_bar_data(_("Finalizing stacking..."), (double)nb_frames/((double)nb_frames+1.));

	copyfits(fit, &gfit, CP_ALLOC|CP_FORMAT, 0);
	gfit.hi = round_to_WORD(maxim);
	gfit.bitpix = USHORT_IMG;
	gfit.exposure = exposure;						// TODO : think if exposure has a sense here

	if (final_pixel[0]) {
		assert(args->seq->nb_layers == 1 || args->seq->nb_layers == 3);
		for (layer=0; layer<args->seq->nb_layers; ++layer){
			from = final_pixel[layer];
			to = gfit.pdata[layer];
			for (y=0; y < fit->ry * fit->rx; ++y) {
				*to++ = *from++;
			}
		}
	}

free_and_reset_progress_bar:
	if (final_pixel[0]) free(final_pixel[0]);
	if (retval) {
		set_progress_bar_data(_("Stacking failed. Check the log."), PROGRESS_RESET);
		siril_log_message(_("Stacking failed.\n"));
	} else {
		set_progress_bar_data(_("Stacking complete."), PROGRESS_DONE);
	}
	update_used_memory();
	return retval;
}

/** Addmax STACK method **
 * This method is very close to the summing one instead that the result
 * takes only the pixel if it is brighter than the previous one at the
 * same coordinates.
 */
int stack_addmin(struct stacking_args *args) {
	int x, y, nx, ny, i, ii, j, shiftx, shifty, layer, reglayer;
	WORD *final_pixel[3], *from;
	WORD *to, minim = USHRT_MAX;
	double exposure=0.0;
	unsigned int nbdata = 0;
	char filename[256];
	int retval = 0;
	int nb_frames, cur_nb = 0;
	fits *fit = &wfit[0];
	char *tmpmsg;
	memset(fit, 0, sizeof(fits));

	/* should be pre-computed to display it in the stacking tab */
	nb_frames = args->nb_images_to_stack;
	reglayer = get_registration_layer();

	if (nb_frames <= 1) {
		siril_log_message(_("No frame selected for stacking (select at least 2). Aborting.\n"));
		return -1;
	}

	final_pixel[0] = NULL;
	assert(args->seq->nb_layers == 1 || args->seq->nb_layers == 3);
	assert(nb_frames <= args->seq->number);

	for (j=0; j<args->seq->number; ++j){
		if (!get_thread_run()) {
			retval = -1;
			goto free_and_reset_progress_bar;
		}
		if (!args->filtering_criterion(args->seq, j, args->filtering_parameter)) {
			fprintf(stdout, "image %d is excluded from stacking\n", j);
			continue;
		}
		if (!seq_get_image_filename(args->seq, j, filename)) {
			retval = -1;
			goto free_and_reset_progress_bar;
		}
		tmpmsg = strdup(_("Processing image "));
		tmpmsg = str_append(&tmpmsg, filename);
		set_progress_bar_data(tmpmsg, (double)cur_nb/((double)nb_frames+1.));
		free(tmpmsg);

		cur_nb++;	// only used for progress bar

		if (seq_read_frame(args->seq, j, fit)) {
			siril_log_message(_("Stacking: could not read frame, aborting\n"));
			retval = -3;
			goto free_and_reset_progress_bar;
		}

		assert(args->seq->nb_layers == 1 || args->seq->nb_layers == 3);
		assert(fit->naxes[2] == args->seq->nb_layers);

		/* first loaded image: init data structures for stacking */
		if (!nbdata) {
			nbdata = fit->ry * fit->rx;
			final_pixel[0] = malloc(nbdata * fit->naxes[2] * sizeof(WORD));
			memset(final_pixel[0], USHRT_MAX, nbdata * fit->naxes[2] * sizeof(WORD));
			if (final_pixel[0] == NULL){
				printf("Stacking: memory allocation failure\n");
				retval = -2;
				goto free_and_reset_progress_bar;
			}
			if(args->seq->nb_layers == 3){
				final_pixel[1] = final_pixel[0] + nbdata;	// index of green layer in final_pixel[0]
				final_pixel[2] = final_pixel[0] + nbdata*2;	// index of blue layer in final_pixel[0]
			}
		} else if (fit->ry * fit->rx != nbdata) {
			siril_log_message(_("Stacking: image in sequence doesn't has the same dimensions\n"));
			retval = -3;
			goto free_and_reset_progress_bar;
		}

		update_used_memory();

		/* load registration data for current image */
		if(reglayer != -1 && args->seq->regparam[reglayer]) {
			shiftx = args->seq->regparam[reglayer][j].shiftx;
			shifty = args->seq->regparam[reglayer][j].shifty;
		} else {
			shiftx = 0;
			shifty = 0;
		}
#ifdef STACK_DEBUG
		printf("stack image %d with shift x=%d y=%d\n", j, shiftx, shifty);
#endif

		/* Summing the exposure */
		exposure += fit->exposure;

		/* stack current image */
		i=0;	// index in final_pixel[0]
		for (y=0; y < fit->ry; ++y){
			for (x=0; x < fit->rx; ++x){
				nx = x - shiftx;
				ny = y - shifty;
				//printf("x=%d y=%d sx=%d sy=%d i=%d ii=%d\n",x,y,shiftx,shifty,i,ii);
				if (nx >= 0 && nx < fit->rx && ny >= 0 && ny < fit->ry) {
					ii = ny * fit->rx + nx;		// index in final_pixel[0] too
					//printf("shiftx=%d shifty=%d i=%d ii=%d\n",shiftx,shifty,i,ii);
					if (ii > 0 && ii < fit->rx * fit->ry){
						for(layer=0; layer<args->seq->nb_layers; ++layer){
							WORD current_pixel = fit->pdata[layer][ii];
							if (current_pixel < final_pixel[layer][i])	// we take the darker pixel
								final_pixel[layer][i] = current_pixel;
							if (final_pixel[layer][i] < minim){
								minim = final_pixel[layer][i];
							}
						}
					}
				}
				++i;
			}
		}
	}
	if (!get_thread_run()) {
		retval = -1;
		goto free_and_reset_progress_bar;
	}
	set_progress_bar_data(_("Finalizing stacking..."), (double)nb_frames/((double)nb_frames+1.));

	copyfits(fit, &gfit, CP_ALLOC|CP_FORMAT, 0);
	gfit.hi = round_to_WORD(minim);
	gfit.bitpix = USHORT_IMG;
	gfit.exposure = exposure;						// TODO : think if exposure has a sense here

	if (final_pixel[0]) {
		assert(args->seq->nb_layers == 1 || args->seq->nb_layers == 3);
		for (layer=0; layer<args->seq->nb_layers; ++layer){
			from = final_pixel[layer];
			to = gfit.pdata[layer];
			for (y=0; y < fit->ry * fit->rx; ++y) {
				*to++ = *from++;
			}
		}
	}

free_and_reset_progress_bar:
	if (final_pixel[0]) free(final_pixel[0]);
	if (retval) {
		set_progress_bar_data(_("Stacking failed. Check the log."), PROGRESS_RESET);
		siril_log_message(_("Stacking failed.\n"));
	} else {
		set_progress_bar_data(_("Stacking complete."), PROGRESS_DONE);
	}
	update_used_memory();
	return retval;
}

static int percentile_clipping(WORD pixel, double sig[], double median, uint64_t rej[]) {
	double plow = sig[0];
	double phigh = sig[1];

	if ((median - (double)pixel) / median > plow) {
		rej[0]++;
		return -1;
	}
	else if (((double)pixel - median) / median > phigh) {
		rej[1]++;
		return 1;
	}
	else return 0;
}

/* Rejection of pixels, following sigma_(high/low) * sigma.
 * The function returns 0 if no rejections are required, 1 if it's a high
 * rejection and -1 for a low-rejection */
static int sigma_clipping(WORD pixel, double sig[], double sigma, double median, uint64_t rej[]) {
	double sigmalow = sig[0];
	double sigmahigh = sig[1];

	if (median - (double)pixel > sigmalow * sigma) {
		rej[0]++;
		return -1;
	}
	else if ((double)pixel - median > sigmahigh * sigma) {
		rej[1]++;
		return 1;
	}
	else return 0;
}

static int Winsorized(WORD *pixel, double m0, double m1) {
	if (*pixel < m0) *pixel = round_to_WORD(m0);
	else if (*pixel > m1) *pixel = round_to_WORD(m1);

	return 0;
}

static int line_clipping(WORD pixel, double sig[], double sigma, int i, double a, double b, uint64_t rej[]) {
	double sigmalow = sig[0];
	double sigmahigh = sig[1];

	if (((a * (double)i + b - (double)pixel) / sigma) > sigmalow) {
		rej[0]++;
		return -1;
	}
	else if ((((double)pixel - a * (double)i - b) / sigma) > sigmahigh) {
		rej[1]++;
		return 1;
	}
	else return 0;
}

static void remove_pixel(WORD *arr, int i, int N) {
	memmove(&arr[i], &arr[i + 1], (N - i - 1) * sizeof(*arr));
}

int stack_mean_with_rejection(struct stacking_args *args) {
	int nb_frames;		/* number of frames actually used */
	int status;		/* CFITSIO status value MUST be initialized to zero for EACH call */
	int reglayer;
	uint64_t irej[3][2] = {{0,0}, {0,0}, {0,0}};
	int bitpix;
	int naxis, oldnaxis = -1, cur_nb = 0;
	long npixels_in_block, nbdata;
	long naxes[3], oldnaxes[3];
	int i;
	double exposure = 0.0;
	char filename[256], msg[256];
	int retval = 0;
	struct _data_block *data_pool = NULL;
	int pool_size = 1;
	fits *fit = &wfit[0];
	norm_coeff coeff;

	nb_frames = args->nb_images_to_stack;
	reglayer = get_registration_layer();

	if (args->seq->type != SEQ_REGULAR && args->seq->type != SEQ_SER) {
		char *msg = siril_log_message(_("Rejection stacking is only supported for FITS images and SER sequences.\nUse \"Sum Stacking\" instead.\n"));
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		return -1;
	}
	if (nb_frames < 2) {
		siril_log_message(_("Select at least two frames for stacking. Aborting.\n"));
		return -1;
	}

	assert(nb_frames <= args->seq->number);
	set_progress_bar_data(NULL, PROGRESS_RESET);

	/* allocate data structures */
	oldnaxes[0] = oldnaxes[1] = oldnaxes[2] = 0;	// fix compiler warning
	naxes[0] = naxes[1] = 0; naxes[2] = 1;

	/* first loop: open all fits files and check they are of same size */
	if (args->seq->type == SEQ_REGULAR) {
		for (i=0; i<nb_frames; ++i) {
			int image_index = args->image_indices[i];	// image index in sequence
			if (!get_thread_run()) {
				retval = -1;
				goto free_and_close;
			}
			if (!fit_sequence_get_image_filename(args->seq, image_index, filename, TRUE))
				continue;

			snprintf(msg, 255, _("Rejection stack: opening image %s"), filename);
			msg[255] = '\0';
			set_progress_bar_data(msg, PROGRESS_NONE);

			/* open input images */
			if (seq_open_image(args->seq, image_index)) {
				retval = -1;
				goto free_and_close;
			}

			/* here we use the internal data of sequences, it's quite ugly, we should
			 * consider moving these tests in seq_open_image() or wrapping them in a
			 * sequence function */
			status = 0;
			fits_get_img_param(args->seq->fptr[image_index], 3, &bitpix, &naxis, naxes, &status);
			if (status) {
				fits_report_error(stderr, status); /* print error message */
				retval = status;
				goto free_and_close;
			}
			if (naxis > 3) {
				siril_log_message(_("Rejection stack error: images with > 3 dimensions "
						"are not supported\n"));
				retval = -1;
				goto free_and_close;
			}

			if(oldnaxis > 0) {
				if(naxis != oldnaxis ||
						oldnaxes[0] != naxes[0] ||
						oldnaxes[1] != naxes[1] ||
						oldnaxes[2] != naxes[2]) {
					siril_log_message(_("Rejection stack error: input images have "
							"different sizes\n"));
					retval = -2;
					goto free_and_close;
				}
			} else {
				oldnaxis = naxis;
				oldnaxes[0] = naxes[0];
				oldnaxes[1] = naxes[1];
				oldnaxes[2] = naxes[2];
			}

			/* exposure summing */
			double tmp;
			status = 0;
			/* and here we should provide a opened_fits_read_key for example */
			fits_read_key (args->seq->fptr[image_index], TDOUBLE, "EXPTIME", &tmp, NULL, &status);
			if (status || tmp <= 0.0) {
				status = 0;
				fits_read_key (args->seq->fptr[image_index], TDOUBLE, "EXPOSURE", &tmp, NULL, &status);
			}
			if (!status)
				exposure += tmp;
		}
		update_used_memory();
	}

	coeff.offset = malloc(nb_frames * sizeof(double));
	coeff.mul = malloc(nb_frames * sizeof(double));
	coeff.scale = malloc(nb_frames * sizeof(double));
	if (!coeff.offset || !coeff.mul || !coeff.scale) {
		printf("allocation issue in stacking normalization\n");
		retval = -1;
		goto free_and_close;
	}

	if (naxes[2] == 0)
		naxes[2] = 1;
	assert(naxes[2] <= 3);
	if (args->seq->type == SEQ_SER) {
		assert(args->seq->ser_file);
		naxes[0] = args->seq->ser_file->image_width;
		naxes[1] = args->seq->ser_file->image_height;
		ser_color type_ser = args->seq->ser_file->color_id;
		if (!com.debayer.open_debayer && type_ser != SER_RGB && type_ser != SER_BGR)
			type_ser = SER_MONO;
		naxes[2] = type_ser == SER_MONO ? 1 : 3;
		naxis = type_ser == SER_MONO ? 2 : 3;
		/* case of Super Pixel not handled yet */
		if (com.debayer.bayer_inter == BAYER_SUPER_PIXEL) {
			siril_log_message(_("Super-pixel is not handled yet for on the fly SER stacking\n"));
			retval = -1;
			goto free_and_close;
		}
	}
	if (naxes[0] == 0) {
		// no image has been loaded
		siril_log_message(_("Rejection stack error: uninitialized sequence\n"));
		retval = -2;
		goto free_and_close;
	}
	fprintf(stdout, "image size: %ldx%ld, %ld layers\n", naxes[0], naxes[1], naxes[2]);

	/* normalization: reading all images and making stats on their background level.
	 * That's very long if not cached. */
	if (compute_normalization(args, &coeff, args->normalize)) {
		retval = -1;
		goto free_and_close;
	}
	if (args->seq->needs_saving)	// if we had to compute new stats
		writeseqfile(args->seq);

	/* initialize result image */
	nbdata = naxes[0] * naxes[1];
	memset(fit, 0, sizeof(fits));
	fit->data = malloc(nbdata * naxes[2] * sizeof(WORD));
	if (!fit->data) {
		fprintf(stderr, "Memory allocation error for result\n");
		retval = -1;
		goto free_and_close;
	}
	fit->bitpix = USHORT_IMG;
	fit->naxes[0] = naxes[0];
	fit->naxes[1] = naxes[1];
	fit->naxes[2] = naxes[2];
	fit->rx = naxes[0];
	fit->ry = naxes[1];
	fit->naxis = naxis;
	fit->maxi = 0;
	if(fit->naxis == 3) {
		fit->pdata[RLAYER]=fit->data;
		fit->pdata[GLAYER]=fit->data + nbdata;
		fit->pdata[BLAYER]=fit->data + nbdata * 2;
	} else {
		fit->pdata[RLAYER]=fit->data;
		fit->pdata[GLAYER]=fit->data;
		fit->pdata[BLAYER]=fit->data;
	}
	update_used_memory();

	/* Define some useful constants */
	double total = (double)(naxes[2] * naxes[1] + 2);	// only used for progress bar

	int nb_threads;
#ifdef _OPENMP
	nb_threads = com.max_thread;
	if (args->seq->type == SEQ_REGULAR && fits_is_reentrant()) {
		fprintf(stdout, "cfitsio was compiled with multi-thread support,"
				" stacking will be executed by several cores\n");
	}
	if (args->seq->type == SEQ_REGULAR && !fits_is_reentrant()) {
		nb_threads = 1;
		fprintf(stdout, "cfitsio was compiled without multi-thread support,"
				" stacking will be executed on only one core\n");
		siril_log_message(_("Your version of cfitsio does not support multi-threading\n"));
	}
#else
	nb_threads = 1;
#endif

	int nb_channels = naxes[2];
	if (sequence_is_rgb(args->seq) && nb_channels != 3) {
		siril_log_message(_("Processing the sequence as RGB\n"));
		nb_channels = 3;
	}

	int size_of_stacks = args->max_number_of_rows / nb_threads; 
	if (size_of_stacks == 0)
		size_of_stacks = 1;
	/* Note: this size of stacks based on the max memory configured doesn't take into
	 * account memory for demosaicing if it applies.
	 * Now we compute the total number of "stacks" which are the independent areas where
	 * the stacking will occur. This will then be used to create the image areas. */
	long nb_parallel_stacks;
	int remainder;
	if (naxes[1] / size_of_stacks < 4) {
		/* We have enough RAM to process each channel with 4 threads.
		 * We should cut images at least in 4 on one channel to use enough threads,
		 * and if only one is available, it will use much less RAM for a small time overhead.
		 * Also, for slow data access like rotating drives or on-the-fly debayer,
		 * it feels more responsive this way.
		 */
		nb_parallel_stacks = 4 * nb_channels;
		size_of_stacks = naxes[1] / 4;
		remainder = naxes[1] % 4;
	} else {
		/* We don't have enough RAM to process a channel with all available threads */
		nb_parallel_stacks = naxes[1] * nb_channels / size_of_stacks;
		if (nb_parallel_stacks % nb_channels != 0
				|| (naxes[1] * nb_channels) % size_of_stacks != 0) {
			/* we need to take into account the fact that the stacks are computed for
			 * each channel, not for the total number of pixels. So it needs to be
			 * a factor of the number of channels.
			 */
			nb_parallel_stacks += nb_channels - (nb_parallel_stacks % nb_channels);
			size_of_stacks = naxes[1] * nb_channels / nb_parallel_stacks;
		}
		remainder = naxes[1] - (nb_parallel_stacks / nb_channels * size_of_stacks);
	}
	siril_log_message(_("We have %d parallel blocks of size %d (+%d) for stacking.\n"),
			nb_parallel_stacks, size_of_stacks, remainder);
	struct image_block {
		unsigned long channel, start_row, end_row, height;
	};
	long largest_block_height = 0;
	struct image_block *blocks = malloc(nb_parallel_stacks * sizeof(struct image_block));
	{
		long channel = 0, row = 0, end, j = 0;
		do {
			if (j >= nb_parallel_stacks) {
				siril_log_message(_("A bug has been found. "
						"Unable to split the image area into the correct processing blocks.\n"));
				retval = -1;
				goto free_and_close;
			}

			blocks[j].channel = channel;
			blocks[j].start_row = row;
			end = row + size_of_stacks - 1; 
			if (remainder > 0) {
				// just add one pixel from the remainder to the first blocks to
				// avoid having all of them in the last block
				end++;
				remainder--;
			}
			if (end >= naxes[1] - 1 ||	// end of the line
					(naxes[1] - end < size_of_stacks / 10)) { // not far from it
				end = naxes[1] - 1;
				row = 0;
				channel++;
				remainder = naxes[1] - (nb_parallel_stacks / nb_channels * size_of_stacks);
			} else {
				row = end + 1;
			}
			blocks[j].end_row = end;
			blocks[j].height = blocks[j].end_row - blocks[j].start_row + 1;
			if (largest_block_height < blocks[j].height) {
				largest_block_height = blocks[j].height;
			}
			fprintf(stdout, "Block %ld: channel %lu, from %lu to %lu (h = %lu)\n",
					j, blocks[j].channel, blocks[j].start_row,
					blocks[j].end_row, blocks[j].height);
			j++;
	
		} while (channel < nb_channels) ;
	}

	/* Allocate the buffers.
	 * We allocate as many as the number of threads, each thread will pick one of the buffers.
	 * Buffers are allocated to the largest block size calculated above.
	 */
#ifdef _OPENMP
	pool_size = nb_threads;
	assert(pool_size > 0);
#endif
	npixels_in_block = largest_block_height * naxes[0];
	fprintf(stdout, "allocating data for %d threads (each %'lu MB)\n", pool_size,
			(unsigned long) (nb_frames * npixels_in_block * sizeof(WORD)) / 1048576UL);
	data_pool = malloc(pool_size * sizeof(struct _data_block));
	for (i = 0; i < pool_size; i++) {
		int j;
		data_pool[i].pix = malloc(nb_frames * sizeof(WORD *));
		data_pool[i].tmp = malloc(nb_frames * npixels_in_block * sizeof(WORD));
		data_pool[i].stack = malloc(nb_frames * sizeof(WORD));
		data_pool[i].rejected = calloc(nb_frames, sizeof(int));
		if (!data_pool[i].pix || !data_pool[i].tmp || !data_pool[i].stack || !data_pool[i].rejected) {
			fprintf(stderr, "Memory allocation error on pix.\n");
			fprintf(stderr, "CHANGE MEMORY SETTINGS if stacking takes too much.\n");
			retval = -1;
			goto free_and_close;
		}
		for (j=0; j<nb_frames; ++j) {
			data_pool[i].pix[j] = data_pool[i].tmp + j * npixels_in_block;
		}
	}
	update_used_memory();

	siril_log_message(_("Starting stacking...\n"));
	set_progress_bar_data(_("Rejection stacking in progress..."), PROGRESS_RESET);

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static) if (args->seq->type == SEQ_SER || fits_is_reentrant())
#endif
	for (i = 0; i < nb_parallel_stacks; i++)
	{
		/**** Step 1: get allocated memory for the current thread ****/
		struct image_block *my_block = blocks+i;
		struct _data_block *data;
		int data_idx = 0, frame;
		long x, y;

		if (!get_thread_run()) retval = -1;
		if (retval) continue;
#ifdef _OPENMP
		data_idx = omp_get_thread_num();
#endif
		assert(data_idx < pool_size);
		//fprintf(stdout, "thread %d working on block %d gets data\n", data_idx, i);
		data = &data_pool[data_idx];

		/**** Step 2: load image data for the corresponding image block ****/
		/* Read the block from all images, store them in pix[image] */
		for (frame = 0; frame < nb_frames; ++frame){
			int shifty = 0;
			gboolean clear = FALSE, readdata = TRUE;
			long offset = 0;
			if (!get_thread_run()) {
				retval = -1;
				break;
			}
			/* area in C coordinates, starting with 0, not cfitsio coordinates. */
			rectangle area = {0, my_block->start_row, naxes[0], my_block->height};

			/* Load registration data for current image and modify area.
			 * Here, only the y shift is managed. If possible, the remaining part
			 * of the original area is read, the rest is filled with zeros. The x
			 * shift is managed in the main loop after the read. */
			if (reglayer != -1 && args->seq->regparam[reglayer]) {
				shifty = args->seq->regparam[reglayer][args->image_indices[frame]].shifty;
				if (area.y + area.h - 1 + shifty < 0 || area.y + shifty >= naxes[1]) {
					// entirely outside image below or above: all black pixels
					clear = TRUE; readdata = FALSE;
				} else if (area.y + shifty < 0) {
					/* we read only the bottom part of the area here, which
					 * requires an offset in pix */
					clear = TRUE;
					area.h += area.y + shifty;	// cropping the height
					offset = naxes[0] * (area.y - shifty);	// positive
					area.y = 0;
				} else if (area.y + area.h - 1 + shifty >= naxes[1]) {
					/* we read only the upper part of the area here */
					clear = TRUE;
					area.y += shifty;
					area.h += naxes[1] - (area.y + area.h);
				} else {
					area.y += shifty;
				}
			}

			if (clear) {
				/* we are reading outside an image, fill with
				 * zeros and attempt to read lines that fit */
				memset(data->pix[frame], 0, npixels_in_block * sizeof(WORD));
			}

			if (readdata) {
				// reading pixels from current frame
				int success = seq_opened_read_region(args->seq, my_block->channel,
						args->image_indices[frame], data->pix[frame]+offset, &area);

				if (success < 0)
					retval = -1;

				if (retval) {
					siril_log_message(_("Error reading one of the image areas\n"));
					break;
				}
			}
		}
		if (retval) continue;

		/**** Step 3: iterate over the y and x of the image block and stack ****/
		for (y = 0; y < my_block->height; y++)
		{
			/* index of the pixel in the result image
			 * we read line y, but we need to store it at
			 * ry - y - 1 to not have the image mirrored. */
			int pdata_idx = (naxes[1] - (my_block->start_row + y) - 1) * naxes[0]; 
			/* index of the line in the read data, data->pix[frame] */
			int pix_idx = y * naxes[0];
			if (retval) break;

			// update progress bar
#ifdef _OPENMP
#pragma omp atomic
#endif
			cur_nb++;

			if (!get_thread_run()) {
				retval = -1;
				break;
			}
			set_progress_bar_data(NULL, (double)cur_nb/total);

			double sigma = -1.0;
			uint64_t crej[2] = {0, 0};
	
			for (x = 0; x < naxes[0]; ++x){
				/* copy all images pixel values in the same row array `stack'
				 * to optimize caching and improve readability */
				for (frame = 0; frame < nb_frames; ++frame) {
					int shiftx = 0;
					if (reglayer != -1 && args->seq->regparam[reglayer]) {
						shiftx = args->seq->regparam[reglayer][args->image_indices[frame]].shiftx;
					}
					if (shiftx && (x - shiftx >= naxes[0] || x - shiftx < 0)) {
						/* outside bounds, images are black. We could
						 * also set the background value instead, if available */
						data->stack[frame] = 0;
					}
					else {
						double tmp;
						switch (args->normalize) {
						default:
						case NO_NORM:		// no normalization (scale[frame] = 1, offset[frame] = 0, mul[frame] = 1)
							data->stack[frame] = data->pix[frame][pix_idx+x-shiftx];
							/* it's faster if we don't convert it to double
							 * to make identity operations */
							break;
						case ADDITIVE:		// additive (scale[frame] = 1, mul[frame] = 1)
						case ADDITIVE_SCALING:		// additive + scale (mul[frame] = 1)
							tmp = (double)data->pix[frame][pix_idx+x-shiftx] * coeff.scale[frame];
							data->stack[frame] = round_to_WORD(tmp - coeff.offset[frame]);
							break;
						case MULTIPLICATIVE:		// multiplicative  (scale[frame] = 1, offset[frame] = 0)
						case MULTIPLICATIVE_SCALING:		// multiplicative + scale (offset[frame] = 0)
							tmp = (double)data->pix[frame][pix_idx+x-shiftx] * coeff.scale[frame];
							data->stack[frame] = round_to_WORD(tmp * coeff.mul[frame]);
							break;
						}
					}
				}

				int N = nb_frames;// N is the number of pixels kept from the current stack
				double median;
				int n, j, r = 0;
				switch (args->type_of_rejection) {
				case PERCENTILE:
					quicksort_s(data->stack, N);
					median = gsl_stats_ushort_median_from_sorted_data(data->stack, 1, N);
					for (frame = 0; frame < N; frame++) {
						data->rejected[frame] =	percentile_clipping(data->stack[frame], args->sig, median, crej);
					}
					for (frame = 0, j = 0; frame < N; frame++, j++) {
						if (data->rejected[j] != 0 && N > 1) {
							remove_pixel(data->stack, frame, N);
							frame--;
							N--;
						}
					}
					break;
				case SIGMA:
					do {
						sigma = gsl_stats_ushort_sd(data->stack, 1, N);
						quicksort_s(data->stack, N);
						median = gsl_stats_ushort_median_from_sorted_data(data->stack, 1, N);
						n = 0;
						for (frame = 0; frame < N; frame++) {
							data->rejected[frame] =	sigma_clipping(data->stack[frame], args->sig, sigma, median, crej);
							if (data->rejected[frame])
								r++;
							if (N - r <= 4) break;
						}
						for (frame = 0, j = 0; frame < N - n; frame++, j++) {
							if (data->rejected[j] != 0) {
								remove_pixel(data->stack, frame, N - n);
								n++;
								frame--;
							}
						}
						N = N - n;
					} while (n > 0 && N > 3);
					break;
				case SIGMEDIAN:
					do {
						sigma = gsl_stats_ushort_sd(data->stack, 1, N);
						quicksort_s(data->stack, N);
						median = gsl_stats_ushort_median_from_sorted_data(data->stack, 1, N);
						n = 0;
						for (frame = 0; frame < N; frame++) {
							if (sigma_clipping(data->stack[frame], args->sig, sigma, median, crej)) {
								data->stack[frame] = round_to_WORD(median);
								n++;
							}
						}
					} while (n > 0 && N > 3);
					break;
				case WINSORIZED:
					do {
						double sigma0;
						sigma = gsl_stats_ushort_sd(data->stack, 1, N);
						quicksort_s(data->stack, N);
						median = gsl_stats_ushort_median_from_sorted_data(data->stack, 1, N);
						WORD *w_stack = malloc(N * sizeof(WORD));
						memcpy(w_stack, data->stack, N * sizeof(WORD));
						do {
							int jj;
							double m0 = median - 1.5 * sigma;
							double m1 = median + 1.5 * sigma;
							for (jj = 0; jj < N; jj++)
								Winsorized(&w_stack[jj], m0, m1);
							quicksort_s(w_stack, N);
							median = gsl_stats_ushort_median_from_sorted_data(w_stack, 1, N);
							sigma0 = sigma;
							sigma = 1.134 * gsl_stats_ushort_sd(w_stack, 1, N);
						} while ((fabs(sigma - sigma0) / sigma0) > 0.0005);
						free(w_stack);
						n = 0;
						for (frame = 0; frame < N; frame++) {
							data->rejected[frame] = sigma_clipping(
									data->stack[frame], args->sig, sigma,
									median, crej);
							if (data->rejected[frame] != 0)
								r++;
							if (N - r <= 4) break;

						}
						for (frame = 0, j = 0; frame < N - n; frame++, j++) {
							if (data->rejected[j] != 0) {
								remove_pixel(data->stack, frame, N - n);
								frame--;
								n++;
							}
						}
						N = N - n;
					} while (n > 0 && N > 3);
					break;
				case LINEARFIT:
					do {
						double *xf = malloc(N * sizeof(double));
						double *yf = malloc(N * sizeof(double));
						double a, b, cov00, cov01, cov11, sumsq;
						quicksort_s(data->stack, N);
						for (frame = 0; frame < N; frame++) {
							xf[frame] = (double) frame;
							yf[frame] = (double) data->stack[frame];
						}
						gsl_fit_linear(xf, 1, yf, 1, N, &b, &a, &cov00, &cov01, &cov11, &sumsq);
						sigma = 0.0;
						for (frame=0; frame < N; frame++)
							sigma += (fabs((double)data->stack[frame] - (a*(double)frame + b)));
						sigma /= (double)N;
						n = 0;
						for (frame = 0; frame < N; frame++) {
							data->rejected[frame] =
									line_clipping(data->stack[frame], args->sig, sigma, frame, a, b, crej);
							if (data->rejected[frame] != 0)
								r++;
							if (N - r <= 4) break;
						}
						for (frame = 0, j = 0; frame < N - n; frame++, j++) {
							if (data->rejected[j] != 0) {
								remove_pixel(data->stack, frame, N - n);
								frame--;
								n++;
							}
						}
						N = N - n;
						free(xf);
						free(yf);
					} while (n > 0 && N > 3);
					break;
				default:
				case NO_REJEC:
					;		// Nothing to do, no rejection
				}

				double sum = 0.0;
				for (frame = 0; frame < N; ++frame) {
					sum += data->stack[frame];
				}
				fit->pdata[my_block->channel][pdata_idx++] = round_to_WORD(sum/(double)N);
			} // end of for x
#ifdef _OPENMP
#pragma omp critical
#endif
			{
				irej[my_block->channel][0] += crej[0];
				irej[my_block->channel][1] += crej[1];
			}

		} // end of for y
	} /* end of loop over parallel stacks */

	if (retval)
		goto free_and_close;

	set_progress_bar_data(_("Finalizing stacking..."), PROGRESS_NONE);
	double nb_tot = (double) naxes[0] * naxes[1] * nb_frames;
	long channel;
	for (channel = 0; channel < naxes[2]; channel++) {
		siril_log_message(_("Pixel rejection in channel #%d: %.3lf%% - %.3lf%%\n"),
				channel, irej[channel][0] / (nb_tot) * 100.0,
				irej[channel][1] / (nb_tot) * 100.0);
	}

	/* copy result to gfit if success */
	copyfits(fit, &gfit, CP_FORMAT, 0);
	if (gfit.data) free(gfit.data);
	gfit.data = fit->data;
	gfit.exposure = exposure;
	memcpy(gfit.pdata, fit->pdata, 3*sizeof(WORD *));

	fit->data = NULL;
	memset(fit->pdata, 0, 3*sizeof(WORD *));

free_and_close:
	fprintf(stdout, "free and close (%d)\n", retval);
	for (i=0; i<nb_frames; ++i) {
		seq_close_image(args->seq, args->image_indices[i]);
	}

	if (data_pool) {
		for (i=0; i<pool_size; i++) {
			if (data_pool[i].stack) free(data_pool[i].stack);
			if (data_pool[i].pix) free(data_pool[i].pix);
			if (data_pool[i].tmp) free(data_pool[i].tmp);
			if (data_pool[i].rejected) free(data_pool[i].rejected);
		}
		free(data_pool);
	}
	free(coeff.offset);
	free(coeff.mul);
	free(coeff.scale);
	if (retval) {
		/* if retval is set, gfit has not been modified */
		if (fit->data) free(fit->data);
		set_progress_bar_data(_("Rejection stacking failed. Check the log."), PROGRESS_RESET);
		siril_log_message(_("Stacking failed.\n"));
	} else {
		set_progress_bar_data(_("Rejection stacking complete."), PROGRESS_DONE);
	}
	update_used_memory();
	return retval;
}

/* the function that runs the thread. Easier to do the simple indirection than
 * changing all return values and adding the idle everywhere. */
gpointer stack_function_handler(gpointer p) {
	struct stacking_args *args = (struct stacking_args *)p;
	args->retval = args->method(p);
	gdk_threads_add_idle(end_stacking, args);
	return GINT_TO_POINTER(args->retval);	// not used anyway
}

/* starts a summing operation using data stored in the stackparam structure
 * function is not reentrant but can be called again after it has returned and the thread is running */
void start_stacking() {
	static GtkComboBox *method_combo = NULL, *rejec_combo = NULL, *norm_combo = NULL;
	static GtkEntry *output_file = NULL;
	static GtkToggleButton *overwrite = NULL, *force_norm = NULL;
	static GtkSpinButton *sigSpin[2] = {NULL, NULL};
	int max_memory;		// maximum memory to use in MB

	if (method_combo == NULL) {
		method_combo = GTK_COMBO_BOX(gtk_builder_get_object(builder, "comboboxstack_methods"));
		output_file = GTK_ENTRY(gtk_builder_get_object(builder, "entryresultfile"));
		overwrite = GTK_TOGGLE_BUTTON(gtk_builder_get_object(builder, "checkbutoverwrite"));
		sigSpin[0] = GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button"));
		sigSpin[1] = GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button"));
		rejec_combo = GTK_COMBO_BOX(lookup_widget("comborejection"));
		norm_combo = GTK_COMBO_BOX(lookup_widget("combonormalize"));
		force_norm = GTK_TOGGLE_BUTTON(lookup_widget("checkforcenorm"));
	}

	if (get_thread_run()) {
		siril_log_message(_("Another task is already in progress, ignoring new request.\n"));
		return;
	}

	stackparam.sig[0] = gtk_spin_button_get_value(sigSpin[0]);
	stackparam.sig[1] = gtk_spin_button_get_value(sigSpin[1]);
	stackparam.type_of_rejection = gtk_combo_box_get_active(rejec_combo);
	stackparam.normalize = gtk_combo_box_get_active(norm_combo);
	stackparam.force_norm = gtk_toggle_button_get_active(force_norm);

	stackparam.method =
			stacking_methods[gtk_combo_box_get_active(method_combo)];
	stackparam.seq = &com.seq;
	max_memory = (int) (com.stack.memory_percent
			* (double) get_available_memory_in_MB());
	siril_log_message(_("Using %d MB memory maximum for stacking\n"), max_memory);
	uint64_t number_of_rows = (uint64_t)max_memory * 1048576L /
		((uint64_t)com.seq.rx * stackparam.nb_images_to_stack * sizeof(WORD) * com.max_thread);
	// this is how many rows we can load in parallel from all images of the
	// sequence and be under the limit defined in config in megabytes.
	// We want to avoid having blocks larger than the half or they will decrease parallelism
	if (number_of_rows > com.seq.ry)
		stackparam.max_number_of_rows = com.seq.ry;
	else if (number_of_rows * 2 > com.seq.ry)
		stackparam.max_number_of_rows = com.seq.ry / 2;
	else stackparam.max_number_of_rows = number_of_rows;

	siril_log_color_message(_("Stacking: processing...\n"), "red");
	gettimeofday(&stackparam.t_start, NULL);
	set_cursor_waiting(TRUE);
	siril_log_message(stackparam.description);

	stackparam.output_overwrite = gtk_toggle_button_get_active(overwrite);
	stackparam.output_filename = gtk_entry_get_text(output_file);

	/* Stacking. Result is in gfit if success */
	start_in_new_thread(stack_function_handler, &stackparam);
}

static void _show_summary(struct stacking_args *args) {
	const char *norm_str, *rej_str;

	siril_log_message(_("Integration of %d images:\n"), args->nb_images_to_stack);

	/* Type of algorithm */
	if (args->method == &stack_mean_with_rejection) {
		siril_log_message(_("Pixel combination ......... average\n"));
	}
	else if (args->method == &stack_summing) {
		siril_log_message(_("Pixel combination ......... normalized sum\n"));
	}
	else if (args->method == &stack_median) {
		siril_log_message(_("Pixel combination ......... median\n"));
	}
	else if (args->method == &stack_addmin) {
		siril_log_message(_("Pixel combination ......... minimum\n"));
	}
	else if (args->method == &stack_addmax) {
		siril_log_message(_("Pixel combination ......... maximum\n"));
	}
	else {
		siril_log_message(_("Pixel combination ......... none\n"));
	}

	/* Normalisation */
	if (args->method != &stack_mean_with_rejection) {
		norm_str = "none";
	} else {
		switch (args->normalize) {
		default:
		case NO_NORM:
			norm_str = _("none");
			break;
		case ADDITIVE:
			norm_str = _("additive");
			break;
		case MULTIPLICATIVE:
			norm_str = _("multiplicative");
			break;
		case ADDITIVE_SCALING:
			norm_str = _("additive + scaling");
			break;
		case MULTIPLICATIVE_SCALING:
			norm_str = _("multiplicative + scaling");
			break;
		}
	}

	siril_log_message(_("Normalization ............. %s\n"), norm_str);

	/* Type of rejection */
	if (args->method != &stack_mean_with_rejection) {
		siril_log_message(_("Pixel rejection ........... none\n"));
		siril_log_message(_("Rejection parameters ...... none\n"));
	}
	else {

		switch (args->type_of_rejection) {
		default:
		case NO_REJEC:
			rej_str = _("none");
			break;
		case PERCENTILE:
			rej_str = _("percentile clipping");
			break;
		case SIGMA:
			rej_str = _("sigma clipping");
			break;
		case SIGMEDIAN:
			rej_str = _("median sigma clipping");
			break;
		case WINSORIZED:
			rej_str = _("Winsorized sigma clipping");
			break;
		case LINEARFIT:
			rej_str = _("linear fit clipping");
			break;
		}
		siril_log_message(_("Pixel rejection ........... %s\n"), rej_str);
		siril_log_message(_("Rejection parameters ...... low=%.3f high=%.3f\n"),
				args->sig[0], args->sig[1]);
	}
}

static void _show_bgnoise(gpointer p) {
	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return;
	}

	struct noise_data *args = malloc(sizeof(struct noise_data));

	set_cursor_waiting(TRUE);

	args->fit = com.uniq->fit;
	args->verbose = FALSE;
	memset(args->bgnoise, 0.0, sizeof(double[3]));
	start_in_new_thread(noise, args);
}

static gboolean end_stacking(gpointer p) {
	struct timeval t_end;
	struct stacking_args *args = (struct stacking_args *)p;
	fprintf(stdout, "Ending stacking idle function, retval=%d\n", args->retval);
	stop_processing_thread();	// can it be done here in case there is no thread?
	if (!args->retval) {
		clear_stars_list();
		com.seq.current = RESULT_IMAGE;
		/* Warning: the previous com.uniq is not freed, but calling
		 * close_single_image() will close everything before reopening it,
		 * which is quite slow */
		com.uniq = calloc(1, sizeof(single));
		com.uniq->comment = strdup("Stacking result image");
		com.uniq->nb_layers = gfit.naxes[2];
		com.uniq->layers = calloc(com.uniq->nb_layers, sizeof(layer_info));
		com.uniq->fit = &gfit;
		com.uniq->fit->maxi = 0;	// force to recompute min/max
		/* Giving summary if average rejection stacking */
		_show_summary(args);
		/* Giving noise estimation */
		_show_bgnoise(com.uniq->fit);
		stop_processing_thread();

		/* save result */
		if (args->output_filename != NULL && args->output_filename[0] != '\0') {
			struct stat st;
			if (!stat(args->output_filename, &st)) {
				int failed = !args->output_overwrite;
				if (!failed) {
					if (unlink(args->output_filename) == -1)
						failed = 1;
					if (!failed && savefits(args->output_filename, &gfit))
						failed = 1;
					if (!failed)
						com.uniq->filename = strdup(args->output_filename);
				}
				if (failed)
					com.uniq->filename = strdup(_("Unsaved stacking result"));
			}
			else {
				if (!savefits(args->output_filename, &gfit))
					com.uniq->filename = strdup(args->output_filename);
				else com.uniq->filename = strdup(_("Unsaved stacking result"));
			}
			display_filename();
		}
		initialize_display_mode();

		adjust_cutoff_from_updated_gfit();
		set_sliders_value_to_gfit();
		initialize_display_mode();

		sliders_mode_set_state(com.sliders);
		set_cutoff_sliders_max_values();

		set_display_mode();

		/* update menus */
		update_MenuItem();

		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
		sequence_list_change_current();
	}

	set_cursor_waiting(FALSE);
#ifdef MAC_INTEGRATION
	GtkosxApplication *osx_app = gtkosx_application_get();
	gtkosx_application_attention_request(osx_app, INFO_REQUEST);
	g_object_unref (osx_app);
#endif
	gettimeofday (&t_end, NULL);
	show_time(args->t_start, t_end);
	return FALSE;
}

void on_seqstack_button_clicked (GtkButton *button, gpointer user_data){
	control_window_switch_to_tab(OUTPUT_LOGS);
	start_stacking();
}

void on_comboboxstack_methods_changed (GtkComboBox *box, gpointer user_data) {
	static GtkNotebook* notebook = NULL;
	if (!notebook)
		notebook = GTK_NOTEBOOK(gtk_builder_get_object(builder, "notebook4"));
	com.stack.method = gtk_combo_box_get_active(box);

	gtk_notebook_set_current_page(notebook, com.stack.method);
	update_stack_interface();
	writeinitfile();
}

void on_comborejection_changed (GtkComboBox *box, gpointer user_data) {
	rejection type_of_rejection = gtk_combo_box_get_active(box);
	GtkLabel *label_rejection[2] = {NULL, NULL};

	if (!label_rejection[0]) {
		label_rejection[0] = GTK_LABEL(lookup_widget("label120"));
		label_rejection[1] = GTK_LABEL(lookup_widget("label122"));
	}
	/* set default values */
	switch (type_of_rejection) {
		case NO_REJEC:
			gtk_widget_set_sensitive(lookup_widget("stack_siglow_button"), FALSE);
			gtk_widget_set_sensitive(lookup_widget("stack_sighigh_button"), FALSE);
			break;
		case PERCENTILE :
			gtk_widget_set_sensitive(lookup_widget("stack_siglow_button"), TRUE);
			gtk_widget_set_sensitive(lookup_widget("stack_sighigh_button"), TRUE);
			gtk_spin_button_set_range (GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button")), 0.0, 1.0);
			gtk_spin_button_set_range (GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button")), 0.0, 1.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button")), 0.2);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button")), 0.1);
			gtk_label_set_text (label_rejection[0], "Percentile low: ");
			gtk_label_set_text (label_rejection[1], "Percentile high: ");
			break;
		case LINEARFIT:
			gtk_widget_set_sensitive(lookup_widget("stack_siglow_button"), TRUE);
			gtk_widget_set_sensitive(lookup_widget("stack_sighigh_button"), TRUE);
			gtk_spin_button_set_range (GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button")), 0.0, 10.0);
			gtk_spin_button_set_range (GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button")), 0.0, 10.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button")), 5.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button")), 2.5);
			gtk_label_set_text (label_rejection[0], "Linear low: ");
			gtk_label_set_text (label_rejection[1], "Linear high: ");
			break;
		default:
		case SIGMA:
		case WINSORIZED:
			gtk_widget_set_sensitive(lookup_widget("stack_siglow_button"), TRUE);
			gtk_widget_set_sensitive(lookup_widget("stack_sighigh_button"), TRUE);
			gtk_spin_button_set_range (GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button")), 0.0, 10.0);
			gtk_spin_button_set_range (GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button")), 0.0, 10.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button")), 4.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button")), 3.0);
			gtk_label_set_text (label_rejection[0], "Sigma low: ");
			gtk_label_set_text (label_rejection[1], "Sigma high: ");
	}
	com.stack.rej_method = gtk_combo_box_get_active(box);
	writeinitfile();
}


/******************* IMAGE FILTERING CRITERIA *******************/

/* a criterion exists for each image filtering method, and is called in a
 * processing to verify if an image should be included or not.
 * These functions have the same signture, defined in stacking.h as
 * stack_filter, and return 1 if the image is included and 0 if not.
 * The functions are also called in a loop to determine the number of images to
 * be processed.
 */
int stack_filter_all(sequence *seq, int nb_img, double any) {
	return 1;
}

int stack_filter_included(sequence *seq, int nb_img, double any) {
	return (seq->imgparam[nb_img].incl);
}

/* filter for deep-sky */
int stack_filter_fwhm(sequence *seq, int nb_img, double max_fwhm) {
	int layer;
	if (!seq->regparam) return 0;
	layer = get_registration_layer();
	if (layer == -1) return 0;
	if (!seq->regparam[layer]) return 0;
	if (seq->imgparam[nb_img].incl && seq->regparam[layer][nb_img].fwhm > 0.0f)
		return (seq->regparam[layer][nb_img].fwhm <= max_fwhm);
	else return 0;
}

/* filter for planetary */
int stack_filter_quality(sequence *seq, int nb_img, double max_quality) {
	int layer;
	if (!seq->regparam) return 0;
	layer = get_registration_layer();
	if (layer == -1) return 0;
	if (!seq->regparam[layer]) return 0;
	if (seq->imgparam[nb_img].incl && seq->regparam[layer][nb_img].quality > 0.0)
		return (seq->regparam[layer][nb_img].quality >= max_quality);
	else return 0;
}

/* browse the images to konw how many fit the criterion, from global data */
int compute_nb_filtered_images() {
	int i, count = 0;
	if (!sequence_is_loaded()) return 0;
	for (i=0; i<com.seq.number; i++) {
		if (stackparam.filtering_criterion(
					&com.seq, i,
					stackparam.filtering_parameter))
			count++;
	}
	return count;
}

/* fill the image_indices mapping for the args->image_indices array, which has
 * to be already allocated to the correct size at least */
void fill_list_of_unfiltered_images(struct stacking_args *args) {
	int i, j;
	for (i=0, j=0; i<args->seq->number; i++) {
		if (args->filtering_criterion(
					&com.seq, i,
					args->filtering_parameter)) {
			args->image_indices[j] = i;
			j++;
		}
	}
	assert(j <= args->nb_images_to_stack);
}

/****************************************************************/

/* For a sequence of images with PSF registration data and a percentage of
 * images to include in a processing, computes the highest FWHM value accepted.
 */
double compute_highest_accepted_fwhm(double percent) {
	int i, layer;
	double *val = malloc(com.seq.number * sizeof(double));
	double highest_accepted;
	layer = get_registration_layer();
	if (layer == -1 || !com.seq.regparam || !com.seq.regparam[layer]) {
		free(val);
		return 0.0;
	}
	// copy values
	for (i=0; i<com.seq.number; i++) {
		if (com.seq.regparam[layer][i].fwhm <= 0.0f) {
			siril_log_message(_("Error in highest FWHM accepted for sequence processing: some images don't have this kind of information available\n"));
			free(val);
			return 0.0;
		}
		val[i] = com.seq.regparam[layer][i].fwhm;
	}

	//sort values
	quicksort_d(val, com.seq.number);
	/*fprintf(stdout, "sorted values:\n");
	  for (i=0; i<com.seq.number; i++)
	  fprintf(stdout, "%g ", val[i]);
	  fputc('\n', stdout);*/

	// get highest accepted
	highest_accepted = val[(int) (percent * (double) com.seq.number / 100.0)];
	free(val);
	return highest_accepted;
}

/* For a sequence of images with quality registration data and a percentage of
 * images to include in a processing, computes the highest quality value accepted.
 */
double compute_highest_accepted_quality(double percent) {
	int i, layer;
	double *val = malloc(com.seq.number * sizeof(double));
	double highest_accepted;
	layer = get_registration_layer();
	if (layer == -1 || !com.seq.regparam || !com.seq.regparam[layer]) {
		free(val);
		return 0.0;
	}
	// copy values
	for (i=0; i<com.seq.number; i++) {
		if (com.seq.imgparam[i].incl && com.seq.regparam[layer][i].quality < 0.0) {
			siril_log_message(_("Error in highest quality accepted for sequence processing: some images don't have this kind of information available for channel #%d.\n"), layer);
			free(val);
			return 0.0;
		}
		else val[i] = com.seq.regparam[layer][i].quality;
	}

	//sort values
	quicksort_d(val, com.seq.number);

	// get highest accepted
	highest_accepted = val[(int) ((100.0 - percent) * (double) com.seq.number
			/ 100.0)];
	free(val);
	return highest_accepted;
}


/* Activates or not the stack button if there are 2 or more selected images,
 * all data related to stacking is set in stackparam, except the method itself,
 * determined at stacking start.
 */
void update_stack_interface() {	// was adjuststackspin
	static GtkAdjustment *stackadj = NULL;
	static GtkWidget *go_stack = NULL, *stack[] = {NULL, NULL}, *widgetnormalize=NULL;
	static GtkComboBox *stack_type = NULL, *method_combo = NULL;
	double percent;
	char labelbuffer[256];

	if(!stackadj) {
		go_stack = lookup_widget("gostack_button");
		stack[0] = lookup_widget("stackspin");
		stackadj = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(stack[0]));
		stack[1] = lookup_widget("label27");
		stack_type = GTK_COMBO_BOX(lookup_widget("comboboxstacksel"));
		method_combo = GTK_COMBO_BOX(lookup_widget("comboboxstack_methods"));
		widgetnormalize = lookup_widget("combonormalize");

	}
	if (!sequence_is_loaded()) return;
	stackparam.seq = &com.seq;

	switch (gtk_combo_box_get_active(method_combo)) {
		default:
		case 0:
		case 3:
		case 4:
			gtk_widget_set_sensitive(widgetnormalize, FALSE);
			break;
		case 1:
		case 2:
			gtk_widget_set_sensitive(widgetnormalize, TRUE);
	}

	switch (gtk_combo_box_get_active(stack_type)) {
		case 0:
			stackparam.filtering_criterion = stack_filter_all;
			stackparam.nb_images_to_stack = com.seq.number;
			sprintf(stackparam.description, _("Stacking all images in the sequence (%d)\n"), com.seq.number);
			gtk_widget_set_sensitive(stack[0], FALSE);
			gtk_widget_set_sensitive(stack[1], FALSE);
			break;
		case 1:
			stackparam.filtering_criterion = stack_filter_included;
			stackparam.nb_images_to_stack = com.seq.selnum;
			sprintf(stackparam.description, _("Stacking only selected images in the sequence (%d)\n"), com.seq.selnum);
			gtk_widget_set_sensitive(stack[0], FALSE);
			gtk_widget_set_sensitive(stack[1], FALSE);
			break;
		case 2:
			/* we should check if the sequence has this kind of data
			 * available before allowing the option to be selected. */
			percent = gtk_adjustment_get_value(stackadj);
			stackparam.filtering_criterion = stack_filter_fwhm;
			stackparam.filtering_parameter = compute_highest_accepted_fwhm(percent);
			stackparam.nb_images_to_stack = compute_nb_filtered_images();
			sprintf(stackparam.description, _("Stacking images of the sequence with a FWHM lower or equal than %g (%d)\n"),
					stackparam.filtering_parameter,
					stackparam.nb_images_to_stack);
			gtk_widget_set_sensitive(stack[0], TRUE);
			gtk_widget_set_sensitive(stack[1], TRUE);
			if (stackparam.filtering_parameter > 0.0)
				sprintf(labelbuffer, _("Based on FWHM < %.2f (%d images)"), stackparam.filtering_parameter, stackparam.nb_images_to_stack);
			else
				sprintf(labelbuffer, _("Based on FWHM"));
			gtk_label_set_text(GTK_LABEL(stack[1]), labelbuffer);
			break;

		case 3:
			percent = gtk_adjustment_get_value(stackadj);
			stackparam.filtering_criterion = stack_filter_quality;
			stackparam.filtering_parameter = compute_highest_accepted_quality(percent);
			stackparam.nb_images_to_stack = compute_nb_filtered_images();
			sprintf(stackparam.description, _("Stacking images of the sequence with a quality higher or equal than %g (%d)\n"),
					stackparam.filtering_parameter,
					stackparam.nb_images_to_stack);
			gtk_widget_set_sensitive(stack[0], TRUE);
			gtk_widget_set_sensitive(stack[1], TRUE);
			if (stackparam.filtering_parameter > 0.0)
				sprintf(labelbuffer, _("Based on quality > %.2f (%d images)"), stackparam.filtering_parameter, stackparam.nb_images_to_stack);
			else
				sprintf(labelbuffer, _("Based on quality"));
			gtk_label_set_text(GTK_LABEL(stack[1]), labelbuffer);
			break;

		default:	// could it be -1?
			fprintf(stderr, "unexpected value from the stack type combo box\n");
			stackparam.nb_images_to_stack = 0;
	}

	if (stackparam.nb_images_to_stack >= 2) {
		if (stackparam.image_indices) free(stackparam.image_indices);
		stackparam.image_indices = malloc(stackparam.nb_images_to_stack * sizeof(int));
		fill_list_of_unfiltered_images(&stackparam);
		gtk_widget_set_sensitive(go_stack, TRUE);
	} else {
		gtk_widget_set_sensitive(go_stack, FALSE);
	}
}

void on_stacksel_changed(GtkComboBox *widget, gpointer user_data) {
	update_stack_interface();
}

void on_spinbut_percent_change(GtkSpinButton *spinbutton, gpointer user_data) {
	update_stack_interface();
}
