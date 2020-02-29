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

#include <string.h>
#include <math.h>
#include <gsl/gsl_statistics_ushort.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "stacking.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "gui/progress_and_log.h"
#include "algos/sorting.h"
#include "stacking/siril_fit_linear.h"

static int stack_mean_or_median(struct stacking_args *args, gboolean is_mean);

/*************************** MEDIAN AND MEAN STACKING **************************
 * Median and mean stacking requires all images to be in memory, so we don't
 * use the generic readfits() but directly the cfitsio routines to open them
 * and seq_opened_read_region() to read data randomly from them.
 * Since all data of all images cannot fit in memory, a divide and conqueer
 * strategy is used, where each thread reads and processes only a part of the
 * image, which size is computed depending on the available memory, image size
 * and thread number.
 *
 * Median stacking does not use registration data, as it's generally used for
 * preprocessing master file creation. Mean stacking however does.
 *
 * The difference between median and mean stacking is that once we have pixel
 * data for all images, in the first case the result is the median of all
 * values, in the other some values can be rejected and the average of the
 * remaining ones is used.
 * ****************************************************************************/

int stack_open_all_files(struct stacking_args *args, int *bitpix, int *naxis, long *naxes, double *exposure, fits *fit) {
	char msg[256], filename[256];
	int i, status, oldbitpix = 0, oldnaxis = -1, nb_frames = args->nb_images_to_stack;
	long oldnaxes[3] = { 0 };
	*exposure = 0.0;

	if (args->seq->type == SEQ_REGULAR) {
		for (i = 0; i < nb_frames; ++i) {
			int image_index = args->image_indices[i];	// image index in sequence
			if (!get_thread_run()) {
				return 1;
			}
			if (!fit_sequence_get_image_filename(args->seq, image_index, filename, TRUE))
				continue;

			snprintf(msg, 255, _("Opening image %s for stacking"), filename);
			msg[255] = '\0';
			set_progress_bar_data(msg, PROGRESS_NONE);

			/* open input images */
			if (seq_open_image(args->seq, image_index)) {
				return 1;
			}

			/* here we use the internal data of sequences, it's quite ugly, we should
			 * consider moving these tests in seq_open_image() or wrapping them in a
			 * sequence function */
			status = 0;
			fits_get_img_param(args->seq->fptr[image_index], 3, bitpix, naxis, naxes, &status);
			if (status) {
				fits_report_error(stderr, status); /* print error message */
				return 1;
			}
			if (*naxis > 3) {
				siril_log_message(_("Stacking error: images with > 3 dimensions "
						"are not supported\n"));
				return 1;
			}

			if (oldnaxis > 0) {
				if (*naxis != oldnaxis ||
						oldnaxes[0] != naxes[0] ||
						oldnaxes[1] != naxes[1] ||
						oldnaxes[2] != naxes[2]) {
					siril_log_message(_("Stacking error: input images have "
							"different sizes\n"));
					return 2;
				}
			} else {
				oldnaxis = *naxis;
				oldnaxes[0] = naxes[0];
				oldnaxes[1] = naxes[1];
				oldnaxes[2] = naxes[2];
			}

			if (oldbitpix > 0) {
				if (*bitpix != oldbitpix) {
					siril_log_message(_("Stacking error: input images have "
							"different precision\n"));
					return 2;
				}
			} else {
				oldbitpix = *bitpix;
			}

			/* exposure summing */
			*exposure += get_exposure_from_fitsfile(args->seq->fptr[image_index]);

			/* We copy metadata from reference to the final fit */
			if (args->seq->type == SEQ_REGULAR && image_index == args->ref_image)
				import_metadata_from_fitsfile(args->seq->fptr[image_index], fit);
		}

		if (naxes[2] == 0)
			naxes[2] = 1;
		g_assert(naxes[2] <= 3);
	}

	else if (args->seq->type == SEQ_SER) {
		g_assert(args->seq->ser_file);
		naxes[0] = args->seq->ser_file->image_width;
		naxes[1] = args->seq->ser_file->image_height;
		ser_color type_ser = args->seq->ser_file->color_id;
		*bitpix = (args->seq->ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8) ? BYTE_IMG : USHORT_IMG;
		if (!com.debayer.open_debayer && type_ser != SER_RGB && type_ser != SER_BGR)
			type_ser = SER_MONO;
		naxes[2] = type_ser == SER_MONO ? 1 : 3;
		*naxis = type_ser == SER_MONO ? 2 : 3;
		/* case of Super Pixel not handled yet */
		if (com.debayer.open_debayer && com.debayer.bayer_inter == BAYER_SUPER_PIXEL) {
			siril_log_message(_("Super-pixel is not handled yet for on the fly SER stacking\n"));
			return 1;
		}
	}
	else {
		siril_log_message(_("Rejection stacking is only supported for FITS images and SER sequences.\nUse \"Sum Stacking\" instead.\n"));
		return 2;
	}

	return 0;
}

int stack_compute_parallel_blocks(struct _image_block **blocksptr, int max_number_of_rows,
		int nb_channels, long *naxes, long *largest_block_height,
		int *nb_blocks, int nb_threads) {
	int size_of_stacks = max_number_of_rows;
	if (size_of_stacks == 0)
		size_of_stacks = 1;
	/* Note: this size of stacks based on the max memory configured doesn't take into
	 * account memory for demosaicing if it applies.
	 * Now we compute the total number of "stacks" which are the independent areas where
	 * the stacking will occur. This will then be used to create the image areas. */
	int remainder;
	if (naxes[1] / size_of_stacks < nb_threads) {
		/* We have enough RAM to process each channel with nb_threads threads.
		 * We should cut images at least in nb_threads on one channel to use enough threads,
		 * and if only one is available, it will use much less RAM for a small time overhead.
		 * Also, for slow data access like rotating drives or on-the-fly debayer,
		 * it feels more responsive this way.
		 */
        int mult = 1;
        // calculate mult so nb_blocks will be a multiple of nb_channels * nb_threads
        while ((mult * nb_channels) % nb_threads) {
            ++mult;
        }

		*nb_blocks = mult * nb_channels;
		size_of_stacks = naxes[1] / mult;
		remainder = naxes[1] % mult;
	} else {
		/* We don't have enough RAM to process a channel with all available threads */
		*nb_blocks = naxes[1] * nb_channels / size_of_stacks;
		if (*nb_blocks % nb_channels != 0
				|| (naxes[1] * nb_channels) % size_of_stacks != 0) {
			/* we need to take into account the fact that the stacks are computed for
			 * each channel, not for the total number of pixels. So it needs to be
			 * a factor of the number of channels.
			 */
			*nb_blocks += nb_channels - (*nb_blocks % nb_channels);
			size_of_stacks = naxes[1] * nb_channels / *nb_blocks;
		}
		remainder = naxes[1] - (*nb_blocks / nb_channels * size_of_stacks);
	}
	siril_log_message(_("We have %d parallel blocks of size %d (+%d) for stacking.\n"),
			*nb_blocks, size_of_stacks, remainder);

	*largest_block_height = 0;
	long channel = 0, row = 0, end, j = 0;
	*blocksptr = malloc(*nb_blocks * sizeof(struct _image_block));
	if (!*blocksptr) return 1;
	struct _image_block *blocks = *blocksptr;
	do {
		if (j >= *nb_blocks) {
			siril_log_message(_("A bug has been found. Unable to split the image "
						"area into the correct processing blocks.\n"));
			return 1;
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
			remainder = naxes[1] - (*nb_blocks / nb_channels * size_of_stacks);
		} else {
			row = end + 1;
		}
		blocks[j].end_row = end;
		blocks[j].height = blocks[j].end_row - blocks[j].start_row + 1;
		if (*largest_block_height < blocks[j].height) {
			*largest_block_height = blocks[j].height;
		}
		fprintf(stdout, "Block %ld: channel %lu, from %lu to %lu (h = %lu)\n",
				j, blocks[j].channel, blocks[j].start_row,
				blocks[j].end_row, blocks[j].height);
		j++;

	} while (channel < nb_channels) ;

	return 0;
}

static void stack_read_block_data(struct stacking_args *args, int use_regdata,
		struct _image_block *my_block, struct _data_block *data,
		long *naxes, data_type itype) {

	int frame, ielem_size = itype == DATA_FLOAT ? sizeof(float) : sizeof(WORD);
	/* Read the block from all images, store them in pix[image] */
	for (frame = 0; frame < args->nb_images_to_stack; ++frame){
		gboolean clear = FALSE, readdata = TRUE;
		long offset = 0;
		/* area in C coordinates, starting with 0, not cfitsio coordinates. */
		rectangle area = {0, my_block->start_row, naxes[0], my_block->height};

		if (!get_thread_run()) {
			return;
		}
		if (use_regdata && args->reglayer >= 0) {
			/* Load registration data for current image and modify area.
			 * Here, only the y shift is managed. If possible, the remaining part
			 * of the original area is read, the rest is filled with zeros. The x
			 * shift is managed in the main loop after the read. */
			regdata *layerparam = args->seq->regparam[args->reglayer];
			if (layerparam) {
				int shifty = round_to_int(
						layerparam[args->image_indices[frame]].shifty *
						args->seq->upscale_at_stacking);
#ifdef STACK_DEBUG
				fprintf(stdout, "shifty for image %d: %d\n", args->image_indices[frame], shifty);
#endif
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
#ifdef STACK_DEBUG
			else fprintf(stderr, "NO REGPARAM\n");
#endif

			if (clear) {
				/* we are reading outside an image, fill with
				 * zeros and attempt to read lines that fit */
				memset(data->pix[frame], 0, my_block->height * naxes[0] * ielem_size);
			}
		}

		if (!use_regdata || readdata) {
			// reading pixels from current frame
			void *buffer;
			if (itype == DATA_FLOAT)
				buffer = ((float*)data->pix[frame])+offset;
			else 	buffer = ((WORD *)data->pix[frame])+offset;
			int retval = seq_opened_read_region(args->seq, my_block->channel,
					args->image_indices[frame], buffer, &area);
			if (retval) {
#ifdef _OPENMP
				int tid = omp_get_thread_num();
				if (tid == 0)
#endif
					siril_log_message(_("Error reading one of the image areas\n"));
				break;
			}
		}
	}
}

static void normalize_to16bit(int bitpix, double *mean) {
	switch(bitpix) {
		case BYTE_IMG:
			*mean *= (USHRT_MAX_DOUBLE / UCHAR_MAX_DOUBLE);
			break;
		default:
			; // do nothing
	}
}

/******************************* REJECTION STACKING ******************************
 * The functions below are those managing the rejection, the stacking code is
 * after and similar to median but takes into account the registration data and
 * does a different operation to keep the final pixel values.
 *********************************************************************************/
static int percentile_clipping(WORD pixel, float sig[], float median, uint64_t rej[]) {
	double plow = sig[0];
	double phigh = sig[1];

	if ((median - pixel) / median > plow) {
		rej[0]++;
		return -1;
	}
	else if ((pixel - median) / median > phigh) {
		rej[1]++;
		return 1;
	}
	return 0;
}

/* Rejection of pixels, following sigma_(high/low) * sigma.
 * The function returns 0 if no rejections are required, 1 if it's a high
 * rejection and -1 for a low-rejection */
static int sigma_clipping(WORD pixel, float sig[], float sigma, float median, uint64_t rej[]) {
	float sigmalow = sig[0];
	float sigmahigh = sig[1];

	if (median - pixel > sigmalow * sigma) {
		rej[0]++;
		return -1;
	}
	else if (pixel - median > sigmahigh * sigma) {
		rej[1]++;
		return 1;
	}
	return 0;
}

static void Winsorize(WORD *pixel, float m0, float m1) {
	if (*pixel < m0) *pixel = roundf_to_WORD(m0);
	else if (*pixel > m1) *pixel = roundf_to_WORD(m1);
}

static int line_clipping(WORD pixel, float sig[], float sigma, int i, float a, float b, uint64_t rej[]) {
	float sigmalow = sig[0];
	float sigmahigh = sig[1];

	if (((a * (float)i + b - pixel) / sigma) > sigmalow) {
		rej[0]++;
		return -1;
	}
	else if (((pixel - a * (float)i - b) / sigma) > sigmahigh) {
		rej[1]++;
		return 1;
	}
	return 0;
}

static int apply_rejection_ushort(struct _data_block *data, int nb_frames, struct stacking_args *args, uint64_t crej[2]) {
	int N = nb_frames;	// N is the number of pixels kept from the current stack
	double median;
	int frame, pixel, output, changed, n, r = 0;
	int firstloop = 1;

	WORD *stack = (WORD *)data->stack;
	WORD *w_stack = (WORD *)data->w_stack;
	int *rejected = (int *)data->rejected;

	/* prepare median and check that the stack is not mostly zero */
	switch (args->type_of_rejection) {
		case PERCENTILE:
		case SIGMA:
		case SIGMEDIAN:
		case WINSORIZED:
			median = quickmedian (stack, N);
			if (median == 0.0)
				return 0;
			break;
		default:
			break;
	}

	switch (args->type_of_rejection) {
		case PERCENTILE:
			for (frame = 0; frame < N; frame++) {
				rejected[frame] = percentile_clipping(stack[frame], args->sig, (float) median, crej);
			}

			for (pixel = 0, output = 0; pixel < N; pixel++) {
				if (!rejected[pixel]) {
					// copy only if there was a rejection
					if (pixel != output)
						stack[output] = stack[pixel];
					output++;
				}
			}
			N = output;
			break;
		case SIGMA:
			do {
				float sigma = (float) gsl_stats_ushort_sd(stack, 1, N);
				if (!firstloop)
					median = quickmedian (stack, N);
				else firstloop = 0;
				for (frame = 0; frame < N; frame++) {
					if (N - r <= 4) {
						// no more rejections
						rejected[frame] = 0;
					} else {
						rejected[frame] = sigma_clipping(stack[frame], args->sig, sigma, (float) median, crej);
						if (rejected[frame])
							r++;
					}
				}
				for (pixel = 0, output = 0; pixel < N; pixel++) {
					if (!rejected[pixel]) {
						// copy only if there was a rejection
						if (pixel != output)
							stack[output] = stack[pixel];
						output++;
					}
				}
				changed = N != output;
				N = output;
			} while (changed && N > 3);
			break;
		case SIGMEDIAN:
			do {
				float sigma = (float) gsl_stats_ushort_sd(stack, 1, N);
				if (!firstloop)
					median = quickmedian (stack, N);
				else firstloop = 0;
				n = 0;
				for (frame = 0; frame < N; frame++) {
					if (sigma_clipping(stack[frame], args->sig, sigma, (float) median, crej)) {
						stack[frame] = median;
						n++;
					}
				}
			} while (n > 0);
			break;
		case WINSORIZED:
			do {
				double sigma0;
				float sigma = (float) gsl_stats_ushort_sd(stack, 1, N);
				if (!firstloop)
					median = quickmedian (stack, N);
				else firstloop = 0;
				memcpy(w_stack, stack, N * sizeof(WORD));
				do {
					int jj;
					double m0 = median - 1.5 * sigma;
					double m1 = median + 1.5 * sigma;
					for (jj = 0; jj < N; jj++)
						Winsorize(w_stack+jj, m0, m1);
					median = quickmedian (w_stack, N);
					sigma0 = sigma;
					sigma = 1.134 * gsl_stats_ushort_sd(w_stack, 1, N);
				} while ((fabs(sigma - sigma0) / sigma0) > 0.0005);
				for (frame = 0; frame < N; frame++) {
					if (N - r <= 4) {
						// no more rejections
						rejected[frame] = 0;
					} else {
						rejected[frame] = sigma_clipping(stack[frame], args->sig, sigma, (float) median, crej);
						if (rejected[frame] != 0)
							r++;
					}

				}
				for (pixel = 0, output = 0; pixel < N; pixel++) {
					if (!rejected[pixel]) {
						// copy only if there was a rejection
						if (pixel != output)
							stack[output] = stack[pixel];
						output++;
					}
				}
				changed = N != output;
				N = output;
			} while (changed && N > 3);
			break;
		case LINEARFIT:
			do {
				float a, b;
				quicksort_s(stack, N);
				for (frame = 0; frame < N; frame++) {
					data->xf[frame] = (float)frame;
					data->yf[frame] = (float)stack[frame];
				}
				siril_fit_linear(data->xf, data->yf, N, &b, &a);
				float sigma = 0.f;
				for (frame = 0; frame < N; frame++)
					sigma += (fabsf((float)stack[frame] - (a*(float)frame + b)));
				sigma /= (float)N;
				for (frame = 0; frame < N; frame++) {
					if (N - r <= 4) {
						// no more rejections
						rejected[frame] = 0;
					} else {
						rejected[frame] =
							line_clipping(stack[frame], args->sig, sigma, frame, a, b, crej);
						if (rejected[frame] != 0)
							r++;
					}
				}
				for (pixel = 0, output = 0; pixel < N; pixel++) {
					if (!rejected[pixel]) {
						// copy only if there was a rejection
						if (pixel != output)
							stack[output] = stack[pixel];
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
	return N;
}

int stack_mean_with_rejection(struct stacking_args *args) {
	return stack_mean_or_median(args, TRUE);
}

int stack_median(struct stacking_args *args) {
	return stack_mean_or_median(args, FALSE);
}

static int stack_mean_or_median(struct stacking_args *args, gboolean is_mean) {
	int nb_frames;		/* number of frames actually used */
	int bitpix, i, naxis, ielem_size, cur_nb = 0, retval = 0, pool_size = 1;
	long npixels_in_block, naxes[3];
	double exposure = 0.0;
	struct _data_block *data_pool = NULL;
	struct _image_block *blocks = NULL;
	data_type itype;	// input data type
	fits fit = { 0 }, *fptr; // output result
	fits ref = { 0 }; // reference image, used to get metadata back
	// data for mean/rej only
	uint64_t irej[3][2] = {{0,0}, {0,0}, {0,0}};
	regdata *layerparam = NULL;
	gboolean use_regdata = is_mean;	// TODO see the other TODO at the top

	nb_frames = args->nb_images_to_stack;
	naxes[0] = naxes[1] = 0; naxes[2] = 1;

	if (nb_frames < 2) {
		siril_log_message(_("Select at least two frames for stacking. Aborting.\n"));
		return -1;
	}
	g_assert(nb_frames <= args->seq->number);

	if (use_regdata) {
		if (args->reglayer < 0)
			fprintf(stderr, "No registration layer passed, ignoring regdata!\n");
		else layerparam = args->seq->regparam[args->reglayer];
	}

	set_progress_bar_data(NULL, PROGRESS_RESET);

	/* first loop: open all fits files and check they are of same size */
	if ((retval = stack_open_all_files(args, &bitpix, &naxis, naxes, &exposure, &ref))) {
		goto free_and_close;
	}

	itype = get_data_type(bitpix);

	if (naxes[0] == 0) {
		// no image has been loaded
		siril_log_message(_("Rejection stack error: uninitialized sequence\n"));
		retval = -2;
		goto free_and_close;
	}
	fprintf(stdout, "image size: %ldx%ld, %ld layers\n", naxes[0], naxes[1], naxes[2]);

	/* initialize result image */
	// TODO: ensure norm_to_16 is not activated when  use_32bit_output is, because it's useless
	fptr = &fit;
	if ((retval = new_fit_image(&fptr, naxes[0], naxes[1], naxes[2],
					args->use_32bit_output ? DATA_FLOAT : DATA_USHORT))) {
		goto free_and_close;
	}
	copy_fits_metadata(&ref, fptr);
	clearfits(&ref);
	if (!args->use_32bit_output && (args->norm_to_16 || fit.orig_bitpix != BYTE_IMG)) {
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
					naxes, &largest_block_height, &nb_blocks, nb_threads))) {
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
	ielem_size = itype == DATA_FLOAT ? sizeof(float) : sizeof(WORD);

	fprintf(stdout, "allocating data for %d threads (each %'lu MB)\n", pool_size,
			(unsigned long)(nb_frames * npixels_in_block * ielem_size) / BYTES_IN_A_MB);
	data_pool = calloc(pool_size, sizeof(struct _data_block));
	for (i = 0; i < pool_size; i++) {
		int j;
		data_pool[i].pix = malloc(nb_frames * sizeof(void *));
		data_pool[i].tmp = malloc(nb_frames * npixels_in_block * ielem_size);
		data_pool[i].stack = malloc(nb_frames * ielem_size);
		if (is_mean)
			data_pool[i].rejected = calloc(nb_frames, sizeof(int));
		if (!data_pool[i].pix || !data_pool[i].tmp || !data_pool[i].stack || (is_mean && !data_pool[i].rejected)) {
			PRINT_ALLOC_ERR;
			fprintf(stderr, "CHANGE MEMORY SETTINGS if stacking takes too much.\n");
			retval = -1;
			goto free_and_close;
		}

		if (is_mean) {
			if (args->type_of_rejection == WINSORIZED) {
				data_pool[i].w_stack = malloc(nb_frames * ielem_size);
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
		}

		for (j=0; j<nb_frames; ++j) {
			if (itype == DATA_FLOAT)
				data_pool[i].pix[j] = ((float*)data_pool[i].tmp) + j * npixels_in_block;
			else 	data_pool[i].pix[j] = ((WORD *)data_pool[i].tmp) + j * npixels_in_block;
		}
	}

	siril_log_message(_("Starting stacking...\n"));
	if (is_mean)
		set_progress_bar_data(_("Rejection stacking in progress..."), PROGRESS_RESET);
	else	set_progress_bar_data(_("Median stacking in progress..."), PROGRESS_RESET);

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
		data = &data_pool[data_idx];

		/**** Step 2: load image data for the corresponding image block ****/
	    stack_read_block_data(args, use_regdata, my_block, data, naxes, itype);

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
			int line_idx = y * naxes[0];
			uint64_t crej[2] = {0, 0};
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
				int frame;
				/* copy all images pixel values in the same row array `stack'
				 * to optimize caching and improve readability */
				for (frame = 0; frame < nb_frames; ++frame) {
					int pix_idx = line_idx + x;
					if (use_regdata) {
						int shiftx = 0;
						if (layerparam) {
							shiftx = round_to_int(
									layerparam[args->image_indices[frame]].shiftx *
									args->seq->upscale_at_stacking);
						}

						if (shiftx && (x - shiftx >= naxes[0] || x - shiftx < 0)) {
							/* outside bounds, images are black. We could
							 * also set the background value instead, if available */
							if (itype == DATA_FLOAT)
								((float*)data->stack)[frame] = 0.0f;
							else	((WORD *)data->stack)[frame] = 0;
							continue;
						}

						pix_idx -= shiftx;
					}

					WORD pixel; float fpixel;
					if (itype == DATA_FLOAT)
						fpixel = ((float*)data->pix[frame])[pix_idx];
					else	pixel  = ((WORD *)data->pix[frame])[pix_idx];
					double tmp;
					switch (args->normalize) {
						default:
						case NO_NORM:
							// no normalization (scale[frame] = 1, offset[frame] = 0, mul[frame] = 1)
							if (itype == DATA_FLOAT)
								((float*)data->stack)[frame] = fpixel;
							else	((WORD *)data->stack)[frame] = pixel;
							/* it's faster if we don't convert it to double
							 * to make identity operations */
							break;
						case ADDITIVE:
							// additive (scale[frame] = 1, mul[frame] = 1)
						case ADDITIVE_SCALING:
							// additive + scale (mul[frame] = 1)
							if (itype == DATA_FLOAT) {
								tmp = fpixel * args->coeff.scale[frame];
								((float*)data->stack)[frame] = (float)(tmp - args->coeff.offset[frame]);
							} else {
								tmp = (double)pixel * args->coeff.scale[frame];
								((WORD *)data->stack)[frame] = round_to_WORD(tmp - args->coeff.offset[frame]);
							}
							break;
						case MULTIPLICATIVE:
							// multiplicative  (scale[frame] = 1, offset[frame] = 0)
						case MULTIPLICATIVE_SCALING:
							// multiplicative + scale (offset[frame] = 0)
							if (itype == DATA_FLOAT) {
								tmp = fpixel * args->coeff.scale[frame];
								((float*)data->stack)[frame] = (float)(tmp * args->coeff.mul[frame]);
							} else {
								tmp = (double)pixel * args->coeff.scale[frame];
								((WORD *)data->stack)[frame] = round_to_WORD(tmp * args->coeff.mul[frame]);
							}
							break;
					}
				}

				double result; // resulting pixel value, either mean or median
				if (is_mean) {
					double mean;
					if (itype == DATA_USHORT) {
						int kept_pixels = apply_rejection_ushort(data, nb_frames, args, crej);
						if (kept_pixels == 0)
							mean = 0.0;
						else {
							int64_t sum = 0L;
							for (frame = 0; frame < kept_pixels; ++frame) {
								sum += ((WORD *)data->stack)[frame];
							}
							mean = sum / (double)kept_pixels;
						}
					} else {
						int kept_pixels = apply_rejection_float(data, nb_frames, args, crej);
						if (kept_pixels == 0)
							mean = 0.0;
						else {
							double sum = 0.0;
							for (frame = 0; frame < kept_pixels; ++frame) {
								sum += ((float*)data->stack)[frame];
							}
							mean = sum / (double)kept_pixels;
						}
					}
					result = mean;
				} else {
					if (itype == DATA_USHORT)
						result = quickmedian(data->stack, nb_frames);
					else 	result = quickmedian_float(data->stack, nb_frames);
				}
				if (args->norm_to_16) {
					normalize_to16bit(bitpix, &result);
				}
				if (args->use_32bit_output) {
					if (itype == DATA_USHORT)
						fit.fpdata[my_block->channel][pdata_idx] = min(double_ushort_to_float_range(result), 1.f);
					else	fit.fpdata[my_block->channel][pdata_idx] = min((float)result, 1.f);
				} else {
					fit.pdata[my_block->channel][pdata_idx] = round_to_WORD(result);
				}
				pdata_idx++;
			} // end of for x

			if (is_mean && args->type_of_rejection != NO_REJEC) {
#ifdef _OPENMP
#pragma omp critical
#endif
				{
					irej[my_block->channel][0] += crej[0];
					irej[my_block->channel][1] += crej[1];
				}
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
	if (is_mean) {
		double nb_tot = (double)(naxes[0] * naxes[1] * nb_frames);
		long channel;
		for (channel = 0; channel < naxes[2]; channel++) {
			siril_log_message(_("Pixel rejection in channel #%d: %.3lf%% - %.3lf%%\n"),
					channel, irej[channel][0] / nb_tot * 100.0,
					irej[channel][1] / nb_tot * 100.0);
		}
	}

	/* copy result to gfit if success */
	clearfits(&gfit);
	copyfits(&fit, &gfit, CP_FORMAT, 0);
	if (args->use_32bit_output) {
		gfit.fdata = fit.fdata;
		for (i = 0; i < fit.naxes[2]; i++)
			gfit.fpdata[i] = fit.fpdata[i];
	} else {
		gfit.data = fit.data;
		for (i = 0; i < fit.naxes[2]; i++)
			gfit.pdata[i] = fit.pdata[i];
	}
	if (is_mean)
		gfit.exposure = exposure;

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
		if (fit.fdata) free(fit.fdata);
		if (is_mean)
			set_progress_bar_data(_("Rejection stacking failed. Check the log."), PROGRESS_RESET);
		else	set_progress_bar_data(_("Median stacking failed. Check the log."), PROGRESS_RESET);
		siril_log_message(_("Stacking failed.\n"));
	} else {
		if (is_mean) {
			set_progress_bar_data(_("Rejection stacking complete."), PROGRESS_DONE);
			siril_log_message(_("Rejection stacking complete. %d images have been stacked.\n"), nb_frames);
		} else {
			set_progress_bar_data(_("Median stacking complete."), PROGRESS_DONE);
			siril_log_message(_("Median stacking complete. %d imageshave been stacked.\n"), nb_frames);
		}
	}
	return retval;
}

