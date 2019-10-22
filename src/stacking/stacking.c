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

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics_ushort.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef MAC_INTEGRATION
#include <gtkosxapplication.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "gui/PSF_list.h"
#include "gui/histogram.h"	// update_gfit_histogram_if_needed();
#include "io/sequence.h"
#include "io/single_image.h"
#include "io/ser.h"
#include "registration/registration.h"
#include "algos/PSF.h"
#include "algos/noise.h"
#include "algos/sorting.h"
#include "stacking.h"
#include "sum.h"
#include "opencv/opencv.h"

static struct stacking_args stackparam = {	// parameters passed to stacking
	NULL, NULL, -1, NULL, -1.0, 0, NULL, NULL, NULL, FALSE, { 0, 0 }, -1, 0,
	{ 0, 0 }, NO_REJEC, NO_NORM, { NULL, NULL, NULL}, FALSE, -1
};

#define MAX_FILTERS 5
static struct filtering_tuple stackfilters[MAX_FILTERS];

stack_method stacking_methods[] = {
	stack_summing_generic, stack_mean_with_rejection, stack_median, stack_addmax, stack_addmin
};

static gboolean end_stacking(gpointer p);
static int stack_addminmax(struct stacking_args *args, gboolean ismax);
static void stacking_args_deep_copy(struct stacking_args *from, struct stacking_args *to);
static void stacking_args_deep_free(struct stacking_args *args);


void initialize_stacking_methods() {
	GtkComboBoxText *stackcombo, *rejectioncombo;

	stackcombo = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(builder, "comboboxstack_methods"));
	rejectioncombo = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(builder, "comborejection"));
	gtk_combo_box_set_active(GTK_COMBO_BOX(stackcombo), com.stack.method);
	gtk_combo_box_set_active(GTK_COMBO_BOX(rejectioncombo), com.stack.rej_method);
}

static void normalize_to16bit(int bitpix, double *mean) {
	switch(bitpix) {
	case BYTE_IMG:
		*mean *= (USHRT_MAX_DOUBLE / UCHAR_MAX_DOUBLE);
		break;
	default:
	case SHORT_IMG:
	case USHORT_IMG:
		; // do nothing
	}
}

/******************************* MEDIAN STACKING ******************************
 * Median stacking requires all images to be in memory, so we dont use the
 * generic readfits but directly the cfitsio routines, and allocates as many
 * pix tables as needed.
 * Median stacking does not use registration data, as it's generally used for
 * preprocessing master file creation.
 * ****************************************************************************/
int stack_median(struct stacking_args *args) {
	int nb_frames;		/* number of frames actually used */
	int bitpix, i, naxis, cur_nb = 0, retval = 0, pool_size = 1;
	long npixels_in_block, naxes[3];
	double exposure;
	struct _data_block *data_pool = NULL;
	struct _image_block *blocks = NULL;
	fits fit = { 0 };

	nb_frames = args->nb_images_to_stack;
	naxes[0] = naxes[1] = 0; naxes[2] = 1;

	if (nb_frames < 2) {
		siril_log_message(_("Select at least two frames for stacking. Aborting.\n"));
		return -1;
	}
	g_assert(nb_frames <= args->seq->number);
	set_progress_bar_data(NULL, PROGRESS_RESET);

	/* first loop: open all fits files and check they are of same size */
	if ((retval = stack_open_all_files(args, &bitpix, &naxis, naxes, &exposure, &fit))) {
		goto free_and_close;
	}

	if (naxes[0] == 0) {
		// no image has been loaded
		siril_log_message(_("Median stack error: uninitialized sequence\n"));
		retval = -2;
		goto free_and_close;
	}
	fprintf(stdout, "image size: %ldx%ld, %ld layers\n", naxes[0], naxes[1], naxes[2]);

	/* initialize result image */
	if ((retval = stack_create_result_fit(&fit, bitpix, naxis, naxes))) {
		goto free_and_close;
	}
	if (args->norm_to_16 || fit.orig_bitpix != BYTE_IMG) {
		fit.bitpix = USHORT_IMG;
		if (args->norm_to_16)
			fit.orig_bitpix = USHORT_IMG;
	}

	/* Define some useful constants */
	double total = (double)(naxes[2] * naxes[1] + 2);	// only used for progress bar

	int nb_threads;
#ifdef _OPENMP
	nb_threads = com.max_thread;
	if (nb_threads > 1 && args->seq->type == SEQ_REGULAR) {
		if (fits_is_reentrant()) {
			fprintf(stdout, "cfitsio was compiled with multi-thread support,"
					" stacking will be executed by several cores\n");
		} else {
			nb_threads = 1;
			fprintf(stdout, "cfitsio was compiled without multi-thread support,"
					" stacking will be executed on only one core\n");
			siril_log_message(_("Your version of cfitsio does not support multi-threading\n"));
		}
	}
#else
	nb_threads = 1;
#endif

	int nb_channels = naxes[2];
	if (sequence_is_rgb(args->seq) && nb_channels != 3) {
		siril_log_message(_("Processing the sequence as RGB\n"));
		nb_channels = 3;
	}


	long largest_block_height;
	int nb_blocks;
	/* Compute parallel processing data: the data blocks, later distributed to threads */
	if ((retval = stack_compute_parallel_blocks(&blocks, args->max_number_of_rows, nb_channels,
					naxes, &largest_block_height, &nb_blocks))) {
		goto free_and_close;
	}

	/* Allocate the buffers.
	 * We allocate as many as the number of threads, each thread will pick one of the buffers.
	 * Buffers are allocated to the largest block size calculated above.
	 */
#ifdef _OPENMP
	pool_size = nb_threads;
	g_assert(pool_size > 0);
#endif
	npixels_in_block = largest_block_height * naxes[0];
	g_assert(npixels_in_block > 0);
	fprintf(stdout, "allocating data for %d threads (each %'lu MB)\n", pool_size,
			(unsigned long) (nb_frames * npixels_in_block * sizeof(WORD)) / BYTES_IN_A_MB);
	data_pool = calloc(pool_size, sizeof(struct _data_block));
	for (i = 0; i < pool_size; i++) {
		int j;
		data_pool[i].pix = calloc(nb_frames, sizeof(WORD *));
		data_pool[i].tmp = calloc(nb_frames, npixels_in_block * sizeof(WORD));
		data_pool[i].stack = calloc(nb_frames, sizeof(WORD));
		if (!data_pool[i].pix || !data_pool[i].tmp || !data_pool[i].stack) {
			PRINT_ALLOC_ERR;
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
#pragma omp parallel for num_threads(nb_threads) private(i) schedule(dynamic) if (nb_threads > 1 && (args->seq->type == SEQ_SER || fits_is_reentrant()))
#endif
	for (i = 0; i < nb_blocks; i++)
	{
		/**** Step 1: get allocated memory for the current thread ****/
		struct _image_block *my_block = blocks+i;
		struct _data_block *data;
		int data_idx = 0, frame;
		long x, y;

		if (!get_thread_run()) retval = -1;
		if (retval) continue;
#ifdef _OPENMP
		data_idx = omp_get_thread_num();
#ifdef STACK_DEBUG
		fprintf(stdout, "Thread %d takes block %d.\n", data_idx, i);
#endif
#endif
		data = &data_pool[data_idx];

		/**** Step 2: load image data for the corresponding image block ****/
		stack_read_block_data(args, 0, my_block, data, naxes);

		/**** Step 3: iterate over the y and x of the image block and stack ****/
		for (y = 0; y < my_block->height; y++)
		{
			/* index of the pixel in the result image
			 * we read line y, but we need to store it at
			 * ry - y - 1 to not have the image mirrored. */
			int pixel_idx = (naxes[1] - (my_block->start_row + y) - 1) * naxes[0]; 
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
			if (!(cur_nb % 16))	// every 16 iterations
				set_progress_bar_data(NULL, (double)cur_nb/total);

			for (x = 0; x < naxes[0]; ++x){
				/* copy all images pixel values in the same row array `stack'
				 * to optimize caching and improve readability */
				for (frame = 0; frame < nb_frames; ++frame) {
					double tmp;
					switch (args->normalize) {
						default:
						case NO_NORM:
							// no normalization (scale[frame] = 1, offset[frame] = 0, mul[frame] = 1)
							data->stack[frame] = data->pix[frame][pix_idx+x];
							/* it's faster if we don't convert it to double
							 * to make identity operations */
							break;
						case ADDITIVE:
							// additive (scale[frame] = 1, mul[frame] = 1)
						case ADDITIVE_SCALING:
							// additive + scale (mul[frame] = 1)
							tmp = (double)data->pix[frame][pix_idx+x] * args->coeff.scale[frame];
							data->stack[frame] = round_to_WORD(tmp - args->coeff.offset[frame]);
							break;
						case MULTIPLICATIVE:
							// multiplicative  (scale[frame] = 1, offset[frame] = 0)
						case MULTIPLICATIVE_SCALING:
							// multiplicative + scale (offset[frame] = 0)
							tmp = (double)data->pix[frame][pix_idx+x] * args->coeff.scale[frame];
							data->stack[frame] = round_to_WORD(tmp * args->coeff.mul[frame]);
							break;
					}
				}
				double median = quickmedian(data->stack, nb_frames);
				if (args->norm_to_16) {
					normalize_to16bit(bitpix, &median);
				}
				fit.pdata[my_block->channel][pixel_idx] = round_to_WORD(median);
				pixel_idx++;
			}
		}
	} /* end of loop over parallel stacks */

	if (retval)
		goto free_and_close;

	set_progress_bar_data(_("Finalizing stacking..."), (double)cur_nb/total);
	/* copy result to gfit if success */
	clearfits(&gfit);
	copyfits(&fit, &gfit, CP_FORMAT, 0);
	gfit.data = fit.data;
	for (i = 0; i < fit.naxes[2]; i++)
		gfit.pdata[i] = fit.pdata[i];

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
	if (blocks) free(blocks);
	if (args->coeff.offset) free(args->coeff.offset);
	if (args->coeff.mul) free(args->coeff.mul);
	if (args->coeff.scale) free(args->coeff.scale);
	if (retval) {
		/* if retval is set, gfit has not been modified */
		if (fit.data) free(fit.data);
		set_progress_bar_data(_("Median stacking failed. Check the log."), PROGRESS_RESET);
		siril_log_message(_("Stacking failed.\n"));
	} else {
		set_progress_bar_data(_("Median stacking complete."), PROGRESS_DONE);
		siril_log_message(_("Median stacking complete. %d have been stacked.\n"), nb_frames);
	}
	update_used_memory();
	return retval;
}


/******************************* ADDMIN AND ADDMAX STACKING ******************************
 * These methods are very close to summing stacking instead that the result
 * takes only the pixel if it is brighter (max) or dimmer (min) than the
 * previous one at the same coordinates.
 */
int stack_addmax(struct stacking_args *args) {
	return stack_addminmax(args, TRUE);
}

int stack_addmin(struct stacking_args *args) {
	return stack_addminmax(args, FALSE);
}

static int stack_addminmax(struct stacking_args *args, gboolean ismax) {
	int x, y, nx, ny, i, ii, j, shiftx, shifty, layer, reglayer;
	WORD *final_pixel[3], *from, *to, minmaxim = ismax ? 0 : USHRT_MAX;;
	double exposure=0.0;
	unsigned int nbdata = 0;
	char filename[256];
	int retval = 0;
	int nb_frames, cur_nb = 0;
	fits fit;
	char *tmpmsg;
	memset(&fit, 0, sizeof(fits));

	/* should be pre-computed to display it in the stacking tab */
	nb_frames = args->nb_images_to_stack;
	reglayer = get_registration_layer(args->seq);

	if (nb_frames <= 1) {
		siril_log_message(_("No frame selected for stacking (select at least 2). Aborting.\n"));
		return -1;
	}

	final_pixel[0] = NULL;
	g_assert(args->seq->nb_layers == 1 || args->seq->nb_layers == 3);
	g_assert(nb_frames <= args->seq->number);

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

		if (seq_read_frame(args->seq, j, &fit)) {
			siril_log_message(_("Stacking: could not read frame, aborting\n"));
			retval = -3;
			goto free_and_reset_progress_bar;
		}

		g_assert(args->seq->nb_layers == 1 || args->seq->nb_layers == 3);
		g_assert(fit.naxes[2] == args->seq->nb_layers);

		/* first loaded image: init data structures for stacking */
		if (!nbdata) {
			nbdata = fit.ry * fit.rx;
			final_pixel[0] = malloc(nbdata * fit.naxes[2] * sizeof(WORD));
			memset(final_pixel[0], ismax ? 0 : USHRT_MAX, nbdata * fit.naxes[2] * sizeof(WORD));
			if (final_pixel[0] == NULL){
				printf("Stacking: memory allocation failure\n");
				retval = -2;
				goto free_and_reset_progress_bar;
			}
			if(args->seq->nb_layers == 3){
				final_pixel[1] = final_pixel[0] + nbdata;	// index of green layer in final_pixel[0]
				final_pixel[2] = final_pixel[0] + nbdata*2;	// index of blue layer in final_pixel[0]
			}
		} else if (fit.ry * fit.rx != nbdata) {
			siril_log_message(_("Stacking: image in sequence doesn't has the same dimensions\n"));
			retval = -3;
			goto free_and_reset_progress_bar;
		}

		update_used_memory();

		/* load registration data for current image */
		if(reglayer != -1 && args->seq->regparam[reglayer]) {
			shiftx = round_to_int(args->seq->regparam[reglayer][j].shiftx * args->seq->upscale_at_stacking);
			shifty = round_to_int(args->seq->regparam[reglayer][j].shifty * args->seq->upscale_at_stacking);
		} else {
			shiftx = 0;
			shifty = 0;
		}
#ifdef STACK_DEBUG
		printf("stack image %d with shift x=%d y=%d\n", j, shiftx, shifty);
#endif

		/* Summing the exposure */
		exposure += fit.exposure;

		/* stack current image */
		i=0;	// index in final_pixel[0]
		for (y=0; y < fit.ry; ++y){
			for (x=0; x < fit.rx; ++x){
				nx = x - shiftx;
				ny = y - shifty;
				//printf("x=%d y=%d sx=%d sy=%d i=%d ii=%d\n",x,y,shiftx,shifty,i,ii);
				if (nx >= 0 && nx < fit.rx && ny >= 0 && ny < fit.ry) {
					ii = ny * fit.rx + nx;		// index in final_pixel[0] too
					//printf("shiftx=%d shifty=%d i=%d ii=%d\n",shiftx,shifty,i,ii);
					if (ii > 0 && ii < fit.rx * fit.ry){
						for(layer=0; layer<args->seq->nb_layers; ++layer){
							WORD current_pixel = fit.pdata[layer][ii];
							if ((ismax && current_pixel > final_pixel[layer][i]) ||	// we take the brighter pixel
									(!ismax && current_pixel < final_pixel[layer][i]))	// we take the darker pixel
								final_pixel[layer][i] = current_pixel;
							if ((ismax && final_pixel[layer][i] > minmaxim) ||
									(!ismax && final_pixel[layer][i] < minmaxim)){
								minmaxim = final_pixel[layer][i];
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

	copyfits(&fit, &gfit, CP_ALLOC|CP_FORMAT, 0);
	gfit.hi = round_to_WORD(minmaxim);
	gfit.bitpix = USHORT_IMG;

	if (final_pixel[0]) {
		g_assert(args->seq->nb_layers == 1 || args->seq->nb_layers == 3);
		for (layer=0; layer<args->seq->nb_layers; ++layer){
			from = final_pixel[layer];
			to = gfit.pdata[layer];
			for (y=0; y < fit.ry * fit.rx; ++y) {
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


/******************************* REJECTION STACKING ******************************
 * The functions below are those managing the rejection, the stacking code is
 * after and similar to median but takes into account the registration data and
 * does a different operation to keep the final pixel values.
 *********************************************************************************/
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

static void Winsorize(WORD *pixel, double m0, double m1) {
	if (*pixel < m0) *pixel = round_to_WORD(m0);
	else if (*pixel > m1) *pixel = round_to_WORD(m1);
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

int stack_mean_with_rejection(struct stacking_args *args) {
	int nb_frames;		/* number of frames actually used */
	uint64_t irej[3][2] = {{0,0}, {0,0}, {0,0}};
	int bitpix;
	int naxis, cur_nb = 0;
	long npixels_in_block;
	long naxes[3];
	int i;
	double exposure = 0.0;
	int retval = 0;
	struct _data_block *data_pool = NULL;
	int pool_size = 1;
	fits fit = { 0 };
	struct _image_block *blocks = NULL;
	regdata *layerparam = NULL;

	nb_frames = args->nb_images_to_stack;
	naxes[0] = naxes[1] = 0; naxes[2] = 1;

	if (nb_frames < 2) {
		siril_log_message(_("Select at least two frames for stacking. Aborting.\n"));
		return -1;
	}
	g_assert(nb_frames <= args->seq->number);

	if (args->reglayer < 0)
		fprintf(stderr, "No registration layer passed, ignoring regdata!\n");
	else layerparam = args->seq->regparam[args->reglayer];

	set_progress_bar_data(NULL, PROGRESS_RESET);

	/* first loop: open all fits files and check they are of same size */
	if ((retval = stack_open_all_files(args, &bitpix, &naxis, naxes, &exposure, &fit))) {
		goto free_and_close;
	}

	if (naxes[0] == 0) {
		// no image has been loaded
		siril_log_message(_("Rejection stack error: uninitialized sequence\n"));
		retval = -2;
		goto free_and_close;
	}
	fprintf(stdout, "image size: %ldx%ld, %ld layers\n", naxes[0], naxes[1], naxes[2]);

	/* initialize result image */
	if ((retval = stack_create_result_fit(&fit, bitpix, naxis, naxes))) {
		goto free_and_close;
	}
	if (args->norm_to_16 || fit.orig_bitpix != BYTE_IMG) {
		fit.bitpix = USHORT_IMG;
		if (args->norm_to_16)
			fit.orig_bitpix = USHORT_IMG;
	}

	/* Define some useful constants */
	double total = (double)(naxes[2] * naxes[1] + 2);	// only used for progress bar

	int nb_threads;
#ifdef _OPENMP
	nb_threads = com.max_thread;
	if (nb_threads > 1 && args->seq->type == SEQ_REGULAR) {
		if (fits_is_reentrant()) {
			fprintf(stdout, "cfitsio was compiled with multi-thread support,"
					" stacking will be executed by several cores\n");
		} else {
			nb_threads = 1;
			fprintf(stdout, "cfitsio was compiled without multi-thread support,"
					" stacking will be executed on only one core\n");
			siril_log_message(_("Your version of cfitsio does not support multi-threading\n"));
		}
	}
#else
	nb_threads = 1;
#endif

	int nb_channels = naxes[2];
	if (sequence_is_rgb(args->seq) && nb_channels != 3) {
		siril_log_message(_("Processing the sequence as RGB\n"));
		nb_channels = 3;
	}

	long largest_block_height;
	int nb_blocks;
	/* Compute parallel processing data: the data blocks, later distributed to threads */
	if ((retval = stack_compute_parallel_blocks(&blocks, args->max_number_of_rows, nb_channels,
					naxes, &largest_block_height, &nb_blocks))) {
		goto free_and_close;
	}

	/* Allocate the buffers.
	 * We allocate as many as the number of threads, each thread will pick one of the buffers.
	 * Buffers are allocated to the largest block size calculated above.
	 */
#ifdef _OPENMP
	pool_size = nb_threads;
	g_assert(pool_size > 0);
#endif
	npixels_in_block = largest_block_height * naxes[0];
	g_assert(npixels_in_block > 0);

	fprintf(stdout, "allocating data for %d threads (each %'lu MB)\n", pool_size,
			(unsigned long) (nb_frames * npixels_in_block * sizeof(WORD)) / BYTES_IN_A_MB);
	data_pool = calloc(pool_size, sizeof(struct _data_block));
	for (i = 0; i < pool_size; i++) {
		int j;
		data_pool[i].pix = malloc(nb_frames * sizeof(WORD *));
		data_pool[i].tmp = malloc(nb_frames * npixels_in_block * sizeof(WORD));
		data_pool[i].stack = malloc(nb_frames * sizeof(WORD));
		data_pool[i].rejected = calloc(nb_frames, sizeof(int));
		if (!data_pool[i].pix || !data_pool[i].tmp || !data_pool[i].stack || !data_pool[i].rejected) {
			PRINT_ALLOC_ERR;
			fprintf(stderr, "CHANGE MEMORY SETTINGS if stacking takes too much.\n");
			retval = -1;
			goto free_and_close;
		}
		if (args->type_of_rejection == WINSORIZED) {
			data_pool[i].w_stack = malloc(nb_frames * sizeof(WORD));
			if (!data_pool[i].w_stack) {
				PRINT_ALLOC_ERR;
				fprintf(stderr, "CHANGE MEMORY SETTINGS if stacking takes too much.\n");
				retval = -1;
				goto free_and_close;
			}
		}

		if (args->type_of_rejection == LINEARFIT) {
			data_pool[i].xf = malloc(nb_frames * sizeof(double));
			data_pool[i].yf = malloc(nb_frames * sizeof(double));
			if (!data_pool[i].xf || !data_pool[i].yf) {
				PRINT_ALLOC_ERR;
				fprintf(stderr, "CHANGE MEMORY SETTINGS if stacking takes too much.\n");
				retval = -1;
				goto free_and_close;
			}
		}

		for (j=0; j<nb_frames; ++j) {
			data_pool[i].pix[j] = data_pool[i].tmp + j * npixels_in_block;
		}
	}
	update_used_memory();

	siril_log_message(_("Starting stacking...\n"));
	set_progress_bar_data(_("Rejection stacking in progress..."), PROGRESS_RESET);

#ifdef _OPENMP
#pragma omp parallel for num_threads(nb_threads) private(i) schedule(dynamic) if (nb_threads > 1 && (args->seq->type == SEQ_SER || fits_is_reentrant()))
#endif
	for (i = 0; i < nb_blocks; i++)
	{
		/**** Step 1: get allocated memory for the current thread ****/
		struct _image_block *my_block = blocks+i;
		struct _data_block *data;
		int data_idx = 0;
		long x, y;

		if (!get_thread_run()) retval = -1;
		if (retval) continue;
#ifdef _OPENMP
		data_idx = omp_get_thread_num();
#ifdef STACK_DEBUG
		struct timeval thread_start;
		gettimeofday(&thread_start, NULL);
		fprintf(stdout, "Thread %d takes block %d.\n", data_idx, i);
#endif
#endif
		//fprintf(stdout, "thread %d working on block %d gets data\n", data_idx, i);
		data = &data_pool[data_idx];

		/**** Step 2: load image data for the corresponding image block ****/
		stack_read_block_data(args, 1, my_block, data, naxes);

#if defined _OPENMP && defined STACK_DEBUG
		{
			struct timeval thread_mid;
			int min, sec;
			gettimeofday(&thread_mid, NULL);
			get_min_sec_from_timevals(thread_start, thread_mid, &min, &sec);
			fprintf(stdout, "Thread %d loaded block %d after %d min %02d s.\n\n",
					data_idx, i, min, sec);
		}
#endif

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
			if (!(cur_nb % 16))	// every 16 iterations
				set_progress_bar_data(NULL, (double)cur_nb/total);

			double sigma = -1.0;
			uint64_t crej[2] = {0, 0};
	
			for (x = 0; x < naxes[0]; ++x){
				int frame;
				/* copy all images pixel values in the same row array `stack'
				 * to optimize caching and improve readability */
				for (frame = 0; frame < nb_frames; ++frame) {
					int shiftx = 0;
					if (layerparam) {
						shiftx = round_to_int(
								layerparam[args->image_indices[frame]].shiftx *
								args->seq->upscale_at_stacking);
					}

					if (shiftx && (x - shiftx >= naxes[0] || x - shiftx < 0)) {
						/* outside bounds, images are black. We could
						 * also set the background value instead, if available */
						data->stack[frame] = 0;
					}
					else {
						WORD pixel = data->pix[frame][pix_idx+x-shiftx];
						double tmp;
						switch (args->normalize) {
						default:
						case NO_NORM:
							// no normalization (scale[frame] = 1, offset[frame] = 0, mul[frame] = 1)
							data->stack[frame] = pixel;
							/* it's faster if we don't convert it to double
							 * to make identity operations */
							break;
						case ADDITIVE:
							// additive (scale[frame] = 1, mul[frame] = 1)
						case ADDITIVE_SCALING:
							// additive + scale (mul[frame] = 1)
							tmp = (double)pixel * args->coeff.scale[frame];
							data->stack[frame] = round_to_WORD(tmp - args->coeff.offset[frame]);
							break;
						case MULTIPLICATIVE:
							// multiplicative  (scale[frame] = 1, offset[frame] = 0)
						case MULTIPLICATIVE_SCALING:
							// multiplicative + scale (offset[frame] = 0)
							tmp = (double)pixel * args->coeff.scale[frame];
							data->stack[frame] = round_to_WORD(tmp * args->coeff.mul[frame]);
							break;
						}
					}
				}

				int N = nb_frames;// N is the number of pixels kept from the current stack
				double median;
				int pixel, output, changed, n, r = 0;
				switch (args->type_of_rejection) {
				case PERCENTILE:
					median = quickmedian (data->stack, N);
					for (frame = 0; frame < N; frame++) {
						data->rejected[frame] =	percentile_clipping(data->stack[frame], args->sig, median, crej);
					}

					for (pixel = 0, output = 0; pixel < N; pixel++) {
						if (!data->rejected[pixel]) {
							// copy only if there was a rejection
							if (pixel != output)
								data->stack[output] = data->stack[pixel];
							output++;
						}
					}
					N = output;
					break;
				case SIGMA:
					do {
						sigma = gsl_stats_ushort_sd(data->stack, 1, N);
						median = quickmedian (data->stack, N);
						for (frame = 0; frame < N; frame++) {
							data->rejected[frame] =	sigma_clipping(data->stack[frame], args->sig, sigma, median, crej);
							if (data->rejected[frame])
								r++;
							if (N - r <= 4) break;
						}
						for (pixel = 0, output = 0; pixel < N; pixel++) {
							if (!data->rejected[pixel]) {
								// copy only if there was a rejection
								if (pixel != output)
									data->stack[output] = data->stack[pixel];
								output++;
							}
						}
						changed = N != output;
						N = output;
					} while (changed && N > 3);
					break;
				case SIGMEDIAN:
					do {
						sigma = gsl_stats_ushort_sd(data->stack, 1, N);
						median = quickmedian (data->stack, N);
						n = 0;
						for (frame = 0; frame < N; frame++) {
							if (sigma_clipping(data->stack[frame], args->sig, sigma, median, crej)) {
								data->stack[frame] = median;
								n++;
							}
						}
					} while (n > 0 && N > 3);
					break;
				case WINSORIZED:
					do {
						double sigma0;
						sigma = gsl_stats_ushort_sd(data->stack, 1, N);
						median = quickmedian (data->stack, N);
						memcpy(data->w_stack, data->stack, N * sizeof(WORD));
						do {
							int jj;
							double m0 = median - 1.5 * sigma;
							double m1 = median + 1.5 * sigma;
							for (jj = 0; jj < N; jj++)
								Winsorize(data->w_stack+jj, m0, m1);
							median = quickmedian (data->w_stack, N);
							sigma0 = sigma;
							sigma = 1.134 * gsl_stats_ushort_sd(data->w_stack, 1, N);
						} while ((fabs(sigma - sigma0) / sigma0) > 0.0005);
						for (frame = 0; frame < N; frame++) {
							data->rejected[frame] = sigma_clipping(
									data->stack[frame], args->sig, sigma,
									median, crej);
							if (data->rejected[frame] != 0)
								r++;
							if (N - r <= 4) break;

						}
						for (pixel = 0, output = 0; pixel < N; pixel++) {
							if (!data->rejected[pixel]) {
								// copy only if there was a rejection
								if (pixel != output)
									data->stack[output] = data->stack[pixel];
								output++;
							}
						}
						changed = N != output;
						N = output;
					} while (changed && N > 3);
					break;
				case LINEARFIT:
					do {
						double a, b, cov00, cov01, cov11, sumsq;
						quicksort_s(data->stack, N);
						for (frame = 0; frame < N; frame++) {
							data->xf[frame] = (double)frame;
							data->yf[frame] = (double)data->stack[frame];
						}
						gsl_fit_linear(data->xf, 1, data->yf, 1, N, &b, &a, &cov00, &cov01, &cov11, &sumsq);
						sigma = 0.0;
						for (frame = 0; frame < N; frame++)
							sigma += (fabs((double)data->stack[frame] - (a*(double)frame + b)));
						sigma /= (double)N;
						for (frame = 0; frame < N; frame++) {
							data->rejected[frame] =
									line_clipping(data->stack[frame], args->sig, sigma, frame, a, b, crej);
							if (data->rejected[frame] != 0)
								r++;
							if (N - r <= 4) break;
						}
						for (pixel = 0, output = 0; pixel < N; pixel++) {
							if (!data->rejected[pixel]) {
								// copy only if there was a rejection
								if (pixel != output)
									data->stack[output] = data->stack[pixel];
								output++;
							}
						}
						changed = N != output;
						N = output;
					} while (changed && N > 3);
					break;
				default:
				case NO_REJEC:
					;		// Nothing to do, no rejection
				}

				int64_t sum = 0L;
				double mean;
				for (frame = 0; frame < N; ++frame) {
					sum += data->stack[frame];
				}
				mean = sum / (double)N;
				if (args->norm_to_16) {
					normalize_to16bit(bitpix, &mean);
				}
				fit.pdata[my_block->channel][pdata_idx++] = round_to_WORD(mean);
			} // end of for x
#ifdef _OPENMP
#pragma omp critical
#endif
			{
				irej[my_block->channel][0] += crej[0];
				irej[my_block->channel][1] += crej[1];
			}

		} // end of for y
#if defined _OPENMP && defined STACK_DEBUG
		{
			struct timeval thread_end;
			int min, sec;
			gettimeofday(&thread_end, NULL);
			get_min_sec_from_timevals(thread_start, thread_end, &min, &sec);
			fprintf(stdout, "Thread %d finishes block %d after %d min %02d s.\n",
					data_idx, i, min, sec);
		}
#endif
	} /* end of loop over parallel stacks */

	if (retval)
		goto free_and_close;

	set_progress_bar_data(_("Finalizing stacking..."), (double)cur_nb/total);
	double nb_tot = (double) naxes[0] * naxes[1] * nb_frames;
	long channel;
	for (channel = 0; channel < naxes[2]; channel++) {
		siril_log_message(_("Pixel rejection in channel #%d: %.3lf%% - %.3lf%%\n"),
				channel, irej[channel][0] / (nb_tot) * 100.0,
				irej[channel][1] / (nb_tot) * 100.0);
	}

	/* copy result to gfit if success */
	clearfits(&gfit);
	copyfits(&fit, &gfit, CP_FORMAT, 0);
	gfit.exposure = exposure;
	gfit.data = fit.data;
	for (i = 0; i < fit.naxes[2]; i++)
		gfit.pdata[i] = fit.pdata[i];

free_and_close:
	fprintf(stdout, "free and close (%d)\n", retval);
	for (i = 0; i < nb_frames; ++i) {
		seq_close_image(args->seq, args->image_indices[i]);
	}

	if (data_pool) {
		for (i=0; i<pool_size; i++) {
			if (data_pool[i].stack) free(data_pool[i].stack);
			if (data_pool[i].pix) free(data_pool[i].pix);
			if (data_pool[i].tmp) free(data_pool[i].tmp);
			if (data_pool[i].rejected) free(data_pool[i].rejected);
			if (data_pool[i].w_stack) free(data_pool[i].w_stack);
			if (data_pool[i].xf) free(data_pool[i].xf);
			if (data_pool[i].yf) free(data_pool[i].yf);
		}
		free(data_pool);
	}
	if (blocks) free(blocks);
	if (args->coeff.offset) free(args->coeff.offset);
	if (args->coeff.mul) free(args->coeff.mul);
	if (args->coeff.scale) free(args->coeff.scale);
	if (retval) {
		/* if retval is set, gfit has not been modified */
		if (fit.data) free(fit.data);
		set_progress_bar_data(_("Rejection stacking failed. Check the log."), PROGRESS_RESET);
		siril_log_message(_("Stacking failed.\n"));
	} else {
		set_progress_bar_data(_("Rejection stacking complete."), PROGRESS_DONE);
	}
	update_used_memory();
	return retval;
}

/* the function that prepares the stacking and runs it */
void main_stack(struct stacking_args *args) {
	int nb_allowed_files;
	assert(args->ref_image >= 0 && args->ref_image < args->seq->number);

	/* first of all we need to check if we can process the files */
	if (args->seq->type == SEQ_REGULAR) {
		if (!allow_to_open_files(args->nb_images_to_stack, &nb_allowed_files)) {
			siril_log_message(_("Your system does not allow one to open more than %d files at the same time. "
						"You may consider either to enhance this limit (the method depends of "
						"your Operating System) or to convert your FITS sequence into a SER "
						"sequence before stacking, or to stack with the \"sum\" method.\n"),
					nb_allowed_files);
			args->retval = -1;
			return;
		}
	}

	siril_log_message(args->description);

	// 1. normalization
	if (do_normalization(args)) // does nothing if NO_NORM
		return;
	// 2. up-scale
	if (upscale_sequence(args)) // does nothing if args->seq->upscale_at_stacking <= 1.05
		return;
	// 3. stack
	args->max_number_of_rows = stack_get_max_number_of_rows(args->seq, args->nb_images_to_stack);
	args->retval = args->method(args);
}

/* the function that runs the thread. */
gpointer stack_function_handler(gpointer p) {
	struct stacking_args *args = (struct stacking_args *)p;

	main_stack(args);

	// 4. save result and clean-up
	siril_add_idle(end_stacking, args);
	return GINT_TO_POINTER(args->retval);
}


/* starts a summing operation using data stored in the stackparam structure
 * function is not reentrant but can be called again after it has returned and the thread is running */
static void start_stacking() {
	static GtkComboBox *method_combo = NULL, *rejec_combo = NULL, *norm_combo = NULL;
	static GtkEntry *output_file = NULL;
	static GtkToggleButton *overwrite = NULL, *force_norm = NULL;
	static GtkSpinButton *sigSpin[2] = {NULL, NULL};
	static GtkWidget *norm_to_16 = NULL;

	if (method_combo == NULL) {
		method_combo = GTK_COMBO_BOX(gtk_builder_get_object(builder, "comboboxstack_methods"));
		output_file = GTK_ENTRY(gtk_builder_get_object(builder, "entryresultfile"));
		overwrite = GTK_TOGGLE_BUTTON(gtk_builder_get_object(builder, "checkbutoverwrite"));
		sigSpin[0] = GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button"));
		sigSpin[1] = GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button"));
		rejec_combo = GTK_COMBO_BOX(lookup_widget("comborejection"));
		norm_combo = GTK_COMBO_BOX(lookup_widget("combonormalize"));
		force_norm = GTK_TOGGLE_BUTTON(lookup_widget("checkforcenorm"));
		norm_to_16 = lookup_widget("check_normalise_to_16b");
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
	stackparam.norm_to_16 = gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(norm_to_16)) && gtk_widget_is_visible(norm_to_16);
	stackparam.coeff.offset = NULL;
	stackparam.coeff.mul = NULL;
	stackparam.coeff.scale = NULL;
	stackparam.method =
			stacking_methods[gtk_combo_box_get_active(method_combo)];
	// ensure we have no normalization if not supported by the stacking method
	if (stackparam.method != stack_median && stackparam.method != stack_mean_with_rejection)
		stackparam.normalize = NO_NORM;
	stackparam.seq = &com.seq;
	stackparam.reglayer = get_registration_layer(&com.seq);
	siril_log_color_message(_("Stacking will use registration data of layer %d if some exist.\n"), "salmon", stackparam.reglayer);

	/* Do not display that cause it uses the generic function that already
	 * displays this text
	 */
	if (stackparam.method != &stack_summing_generic)
		siril_log_color_message(_("Stacking: processing...\n"), "red");
	gettimeofday(&stackparam.t_start, NULL);
	set_cursor_waiting(TRUE);

	stackparam.output_overwrite = gtk_toggle_button_get_active(overwrite);
	stackparam.output_filename = gtk_entry_get_text(output_file);

	/* Stacking. Result is in gfit if success */
	struct stacking_args *params = malloc(sizeof(struct stacking_args));
	stacking_args_deep_copy(&stackparam, params);
	start_in_new_thread(stack_function_handler, params);
}

static void _show_summary(struct stacking_args *args) {
	const char *norm_str, *rej_str;

	siril_log_message(_("Integration of %d images:\n"), args->nb_images_to_stack);

	/* Type of algorithm */
	if (args->method == &stack_mean_with_rejection) {
		siril_log_message(_("Pixel combination ......... average\n"));
	} else if (args->method == &stack_summing_generic) {
		siril_log_message(_("Pixel combination ......... normalized sum\n"));
	} else if (args->method == &stack_median) {
		siril_log_message(_("Pixel combination ......... median\n"));
	} else if (args->method == &stack_addmin) {
		siril_log_message(_("Pixel combination ......... minimum\n"));
	} else if (args->method == &stack_addmax) {
		siril_log_message(_("Pixel combination ......... maximum\n"));
	} else {
		siril_log_message(_("Pixel combination ......... none\n"));
	}

	/* Normalisation */
	if (args->method != &stack_mean_with_rejection &&
			args->method != &stack_median ) {
		norm_str = _("none");
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
	set_cursor_waiting(TRUE);

	struct noise_data *args = malloc(sizeof(struct noise_data));
	args->fit = com.uniq->fit;
	args->verbose = FALSE;
	args->use_idle = TRUE;
	memset(args->bgnoise, 0.0, sizeof(double[3]));

	start_in_new_thread(noise, args);
}

void clean_end_stacking(struct stacking_args *args) {
	if (!args->retval)
		_show_summary(args);
	remove_tmp_drizzle_files(args);
}

/* because this idle function is called after one of many stacking method
 * functions, it contains all generic wrap-up stuff instead of only graphical
 * operations. */
static gboolean end_stacking(gpointer p) {
	struct timeval t_end;
	struct stacking_args *args = (struct stacking_args *)p;
	fprintf(stdout, "Ending stacking idle function, retval=%d\n", args->retval);
	stop_processing_thread();	// can it be done here in case there is no thread?

	if (!args->retval) {
		clear_stars_list();
		/* check in com.seq, because args->seq may have been replaced */
		if (com.seq.upscale_at_stacking > 1.05)
			com.seq.current = SCALED_IMAGE;
		else com.seq.current = RESULT_IMAGE;
		/* Warning: the previous com.uniq is not freed, but calling
		 * close_single_image() will close everything before reopening it,
		 * which is quite slow */
		com.uniq = calloc(1, sizeof(single));
		com.uniq->comment = strdup(_("Stacking result image"));
		com.uniq->nb_layers = gfit.naxes[2];
		com.uniq->layers = calloc(com.uniq->nb_layers, sizeof(layer_info));
		com.uniq->fit = &gfit;
		/* Giving summary if average rejection stacking */
		_show_summary(args);
		/* Giving noise estimation (new thread) */
		_show_bgnoise(com.uniq->fit);

		/* save stacking result */
		if (args->output_filename != NULL && args->output_filename[0] != '\0') {
			GStatBuf st;
			if (!g_stat(args->output_filename, &st)) {
				int failed = !args->output_overwrite;
				if (!failed) {
					if (g_unlink(args->output_filename) == -1)
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
		/* remove tmp files if exist (Drizzle) */
		remove_tmp_drizzle_files(args);

		waiting_for_thread();		// bgnoise
		adjust_cutoff_from_updated_gfit();	// computes min and max
		set_sliders_value_to_gfit();
		initialize_display_mode();

		sliders_mode_set_state(com.sliders);
		set_cutoff_sliders_max_values();

		set_display_mode();

		/* update menus */
		update_MenuItem();

		if (com.seq.current == SCALED_IMAGE)
			adjust_vport_size_to_image();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
		sequence_list_change_current();
		update_stack_interface(TRUE);
	}

	set_cursor_waiting(FALSE);
#ifdef MAC_INTEGRATION
	GtkosxApplication *osx_app = gtkosx_application_get();
	gtkosx_application_attention_request(osx_app, INFO_REQUEST);
	g_object_unref (osx_app);
#endif
	/* Do not display time for stack_summing_generic
	 * cause it uses the generic function that already
	 * displays the time
	 */
	if (args->method != &stack_summing_generic) {
		gettimeofday(&t_end, NULL);
		show_time(args->t_start, t_end);
	}
	stacking_args_deep_free(args);
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
	update_stack_interface(TRUE);
	writeinitfile();
}

void on_combonormalize_changed (GtkComboBox *box, gpointer user_data) {
	GtkWidget *widgetnormalize = lookup_widget("combonormalize");
	GtkWidget *force_norm = lookup_widget("checkforcenorm");
	gtk_widget_set_sensitive(force_norm,
			gtk_combo_box_get_active(GTK_COMBO_BOX(widgetnormalize)) != 0);
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
			gtk_label_set_text (label_rejection[0], _("Percentile low: "));
			gtk_label_set_text (label_rejection[1], _("Percentile high: "));
			break;
		case LINEARFIT:
			gtk_widget_set_sensitive(lookup_widget("stack_siglow_button"), TRUE);
			gtk_widget_set_sensitive(lookup_widget("stack_sighigh_button"), TRUE);
			gtk_spin_button_set_range (GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button")), 0.0, 10.0);
			gtk_spin_button_set_range (GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button")), 0.0, 10.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("stack_siglow_button")), 5.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("stack_sighigh_button")), 5.0);
			gtk_label_set_text (label_rejection[0], _("Linear low: "));
			gtk_label_set_text (label_rejection[1], _("Linear high: "));
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
			gtk_label_set_text (label_rejection[0], _("Sigma low: "));
			gtk_label_set_text (label_rejection[1], _("Sigma high: "));
	}
	com.stack.rej_method = gtk_combo_box_get_active(box);
	writeinitfile();
}

int stack_get_max_number_of_rows(sequence *seq, int nb_images_to_stack) {
	int max_memory = get_max_memory_in_MB();
	if (max_memory > 0) {
		siril_log_message(_("Using %d MB memory maximum for stacking\n"), max_memory);
		uint64_t number_of_rows = (uint64_t)max_memory * BYTES_IN_A_MB /
			((uint64_t)seq->rx * nb_images_to_stack * sizeof(WORD) * com.max_thread);
		// this is how many rows we can load in parallel from all images of the
		// sequence and be under the limit defined in config in megabytes.
		// We want to avoid having blocks larger than the half or they will decrease parallelism
		if (number_of_rows > seq->ry)
			return seq->ry;
		if (number_of_rows * 2 > seq->ry)
			return seq->ry / 2;
		return truncate_to_int32(number_of_rows);
	} else {
		siril_log_message(_("Not using limits on maximum memory for stacking\n"));
		return (seq->ry / 4) + 1;
	}
}

int find_refimage_in_indices(int *indices, int nb, int ref) {
	int i;
	for (i = 0; i < nb; i++) {
		if (indices[i] == ref)
			return i;
	}
	return -1;
}

/****************************************************************/

void on_stacksel_changed(GtkComboBox *widget, gpointer user_data) {
	update_stack_interface(TRUE);
}

void on_spinbut_percent_change(GtkSpinButton *spinbutton, gpointer user_data) {
	update_stack_interface(TRUE);
}

void on_filter_add1_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(lookup_widget("combofilter2"), TRUE);
	gtk_widget_set_visible(lookup_widget("stackspin2"), TRUE);
	gtk_widget_set_visible(lookup_widget("filter_add2"), TRUE);
	gtk_widget_set_visible(lookup_widget("filter_rem2"), TRUE);
	gtk_widget_set_visible(lookup_widget("labelfilter2"), TRUE);
	update_stack_interface(TRUE);
}

void on_filter_add2_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(lookup_widget("combofilter3"), TRUE);
	gtk_widget_set_visible(lookup_widget("stackspin3"), TRUE);
	gtk_widget_set_visible(lookup_widget("filter_rem3"), TRUE);
	gtk_widget_set_visible(lookup_widget("labelfilter3"), TRUE);
	update_stack_interface(TRUE);
}

void on_filter_rem2_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(lookup_widget("combofilter2"), FALSE);
	gtk_widget_set_visible(lookup_widget("stackspin2"), FALSE);
	gtk_widget_set_visible(lookup_widget("filter_add2"), FALSE);
	gtk_widget_set_visible(lookup_widget("filter_rem2"), FALSE);
	gtk_widget_set_visible(lookup_widget("labelfilter2"), FALSE);
	update_stack_interface(TRUE);
}

void on_filter_rem3_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_set_visible(lookup_widget("combofilter3"), FALSE);
	gtk_widget_set_visible(lookup_widget("stackspin3"), FALSE);
	gtk_widget_set_visible(lookup_widget("filter_rem3"), FALSE);
	gtk_widget_set_visible(lookup_widget("labelfilter3"), FALSE);
	update_stack_interface(TRUE);
}

void get_sequence_filtering_from_gui(seq_image_filter *filtering_criterion,
		double *filtering_parameter) {
	int filter, guifilter, channel = 0, type;
	double percent = 0.0;
	static GtkComboBox *filter_combo[] = {NULL, NULL, NULL};
	static GtkAdjustment *stackadj[] = {NULL, NULL, NULL};
	static GtkWidget *spin[] = {NULL, NULL, NULL};
	if (!spin[0]) {
		spin[0] = lookup_widget("stackspin1");
		spin[1] = lookup_widget("stackspin2");
		spin[2] = lookup_widget("stackspin3");
		stackadj[0] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[0]));
		stackadj[1] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[1]));
		stackadj[2] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[2]));
		filter_combo[0] = GTK_COMBO_BOX(lookup_widget("combofilter1"));
		filter_combo[1] = GTK_COMBO_BOX(lookup_widget("combofilter2"));
		filter_combo[2] = GTK_COMBO_BOX(lookup_widget("combofilter3"));
	}
	for (filter = 0, guifilter = 0; guifilter < 3; guifilter++) {
		if (!gtk_widget_get_visible(GTK_WIDGET(filter_combo[guifilter]))) {
			continue;
		}

		type = gtk_combo_box_get_active(filter_combo[guifilter]);
		if (type != ALL_IMAGES && type != SELECTED_IMAGES) {
			channel = get_registration_layer(&com.seq);
			percent = gtk_adjustment_get_value(stackadj[guifilter]);
		}

		switch (type) {
			default:
			case ALL_IMAGES:
				stackfilters[filter].filter = seq_filter_all;
				stackfilters[filter].param = 0.0;
				gtk_widget_set_visible(spin[guifilter], FALSE);
				break;
			case SELECTED_IMAGES:
				stackfilters[filter].filter = seq_filter_included;
				stackfilters[filter].param = 0.0;
				gtk_widget_set_visible(spin[guifilter], FALSE);
				break;
			case BEST_PSF_IMAGES:
				stackfilters[filter].filter = seq_filter_fwhm;
				stackfilters[filter].param = compute_highest_accepted_fwhm(
						stackparam.seq, channel, percent);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				break;
			case BEST_ROUND_IMAGES:
				stackfilters[filter].filter = seq_filter_roundness;
				stackfilters[filter].param = compute_lowest_accepted_roundness(
						stackparam.seq, channel, percent);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				break;
			case BEST_QUALITY_IMAGES:
				stackfilters[filter].filter = seq_filter_quality;
				stackfilters[filter].param = compute_lowest_accepted_quality(
						stackparam.seq, channel, percent);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				break;
		}
		filter++;
	}
	stackfilters[filter].filter = NULL;

	if (filter == 1) {
		*filtering_criterion = stackfilters[0].filter;
		*filtering_parameter = stackfilters[0].param;
	} else {
		*filtering_criterion = create_multiple_filter_from_list(stackfilters);
		*filtering_parameter = 0.0;
	}
}

static void update_filter_label() {
	static GtkComboBox *filter_combo[] = {NULL, NULL, NULL};
	static GtkLabel *filter_label[] = {NULL, NULL, NULL};
	gchar *filter_str;
	double param;
	int filter, type;

	if (!filter_combo[0]) {
		filter_combo[0] = GTK_COMBO_BOX(lookup_widget("combofilter1"));
		filter_combo[1] = GTK_COMBO_BOX(lookup_widget("combofilter2"));
		filter_combo[2] = GTK_COMBO_BOX(lookup_widget("combofilter3"));
		filter_label[0] = GTK_LABEL(lookup_widget("labelfilter1"));
		filter_label[1] = GTK_LABEL(lookup_widget("labelfilter2"));
		filter_label[2] = GTK_LABEL(lookup_widget("labelfilter3"));
	}

	for (filter = 0; filter < 3; filter++) {
		if (!gtk_widget_get_visible(GTK_WIDGET(filter_combo[filter]))) {
			break;
		}

		type = gtk_combo_box_get_active(filter_combo[filter]);
		param = stackfilters[filter].param;
		if (param == DBL_MIN || param == DBL_MAX || param == 0.0) {
			if (type == ALL_IMAGES || type == SELECTED_IMAGES)
				filter_str = g_strdup("");
			else filter_str = g_strdup("N/A");
		} else {
			switch (type) {
			default:
			case ALL_IMAGES:
			case SELECTED_IMAGES:
				filter_str = g_strdup("");
				break;
			case BEST_PSF_IMAGES:
				filter_str = g_strdup_printf("< %.2lf", param);
				break;
			case BEST_ROUND_IMAGES:
			case BEST_QUALITY_IMAGES:
				filter_str = g_strdup_printf("> %.3lf", param);
				break;
			}
		}
		gtk_label_set_text(filter_label[filter], filter_str);
		g_free(filter_str);
	}
}

/* Activates or not the stack button if there are 2 or more selected images,
 * all data related to stacking is set in stackparam, except the method itself,
 * determined at stacking start.
 */
void update_stack_interface(gboolean dont_change_stack_type) {
	static GtkWidget *go_stack = NULL,
			 *widgetnormalize = NULL, *force_norm = NULL, *norm_to_16 = NULL;
	static GtkComboBox *method_combo = NULL, *filter_combo = NULL;
	static GtkLabel *result_label = NULL;
	gchar *labelbuffer;

	if(!go_stack) {
		go_stack = lookup_widget("gostack_button");
		filter_combo = GTK_COMBO_BOX(lookup_widget("combofilter1"));
		method_combo = GTK_COMBO_BOX(lookup_widget("comboboxstack_methods"));
		widgetnormalize = lookup_widget("combonormalize");
		force_norm = lookup_widget("checkforcenorm");
		norm_to_16 = lookup_widget("check_normalise_to_16b");
		result_label = GTK_LABEL(lookup_widget("stackfilter_label"));
	}
	if (!sequence_is_loaded()) {
		gtk_widget_set_sensitive(go_stack, FALSE);
		return;
	}
	stackparam.seq = &com.seq;

	if (!dont_change_stack_type && stackparam.seq->selnum < stackparam.seq->number) {
		g_signal_handlers_block_by_func(filter_combo, on_stacksel_changed, NULL);
		gtk_combo_box_set_active(filter_combo, SELECTED_IMAGES);
		g_signal_handlers_unblock_by_func(filter_combo, on_stacksel_changed, NULL);
	}

	switch (gtk_combo_box_get_active(method_combo)) {
	default:
	case STACK_SUM:
	case STACK_MAX:
	case STACK_MIN:
		gtk_widget_set_sensitive(widgetnormalize, FALSE);
		gtk_widget_set_sensitive(force_norm, FALSE);
		break;
	case STACK_MEAN:
	case STACK_MEDIAN:
		gtk_widget_set_sensitive(widgetnormalize, TRUE);
		gtk_widget_set_sensitive(force_norm,
				gtk_combo_box_get_active(GTK_COMBO_BOX(widgetnormalize)) != 0);
		gtk_widget_set_visible(norm_to_16, stackparam.seq->bitpix == BYTE_IMG);
	}

	if (com.seq.reference_image == -1)
		com.seq.reference_image = sequence_find_refimage(&com.seq);
	stackparam.ref_image = com.seq.reference_image;

	get_sequence_filtering_from_gui(
			&stackparam.filtering_criterion, &stackparam.filtering_parameter);

	if (stackparam.description)
		free(stackparam.description);
	stackparam.description = describe_filter(stackparam.seq,
			stackparam.filtering_criterion, stackparam.filtering_parameter);

	update_filter_label();

	stackparam.nb_images_to_stack = compute_nb_filtered_images(&com.seq,
			stackparam.filtering_criterion, stackparam.filtering_parameter);
	labelbuffer = g_strdup_printf(_("Stacking %d images of the %d of the sequence"),
			stackparam.nb_images_to_stack, com.seq.number);
	gtk_label_set_text(result_label, labelbuffer);
	g_free(labelbuffer);

	if (stackparam.nb_images_to_stack >= 2) {
		stack_fill_list_of_unfiltered_images(&stackparam);
		gtk_widget_set_sensitive(go_stack, TRUE);
	} else {
		gtk_widget_set_sensitive(go_stack, FALSE);
	}
}

static void stacking_args_deep_copy(struct stacking_args *from, struct stacking_args *to) {
	memcpy(to, from, sizeof(struct stacking_args));
	// sequence is not duplicated
	to->image_indices = malloc(from->nb_images_to_stack * sizeof(int));
	memcpy(to->image_indices, from->image_indices, from->nb_images_to_stack * sizeof(int));
	to->description = strdup(from->description);
	// output_filename is not duplicated, can be changed until the last minute
}

static void stacking_args_deep_free(struct stacking_args *args) {
	free(args->image_indices);
	free(args->description);
	free(args);
}
