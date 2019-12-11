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

#include <string.h>
#include <math.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics_ushort.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "stacking.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "gui/progress_and_log.h"
#include "algos/sorting.h"

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
		int *nb_blocks) {
	int size_of_stacks = max_number_of_rows;
	if (size_of_stacks == 0)
		size_of_stacks = 1;
	/* Note: this size of stacks based on the max memory configured doesn't take into
	 * account memory for demosaicing if it applies.
	 * Now we compute the total number of "stacks" which are the independent areas where
	 * the stacking will occur. This will then be used to create the image areas. */
	int remainder;
	if (naxes[1] / size_of_stacks < 4) {
		/* We have enough RAM to process each channel with 4 threads.
		 * We should cut images at least in 4 on one channel to use enough threads,
		 * and if only one is available, it will use much less RAM for a small time overhead.
		 * Also, for slow data access like rotating drives or on-the-fly debayer,
		 * it feels more responsive this way.
		 */
		*nb_blocks = 4 * nb_channels;
		size_of_stacks = naxes[1] / 4;
		remainder = naxes[1] % 4;
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

void stack_read_block_data(struct stacking_args *args, int use_regdata,
		struct _image_block *my_block, struct _data_block *data, long *naxes) {

	int frame;
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
				memset(data->pix[frame], 0, my_block->height * naxes[0] * sizeof(WORD));
			}
		}

		if (!use_regdata || readdata) {
			// reading pixels from current frame
			int retval = seq_opened_read_region(args->seq, my_block->channel,
					args->image_indices[frame], data->pix[frame]+offset, &area);
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

int apply_rejection_ushort(struct _data_block *data, int nb_frames, struct stacking_args *args, uint64_t crej[2]) {
	int N = nb_frames;// N is the number of pixels kept from the current stack
	double median, sigma = -1.0;
	int frame, pixel, output, changed, n, r = 0;

	WORD *stack = (WORD *)data->stack;
	WORD *w_stack = (WORD *)data->w_stack;
	WORD *rejected = (WORD *)data->rejected;

	switch (args->type_of_rejection) {
		case PERCENTILE:
			median = quickmedian (stack, N);
			for (frame = 0; frame < N; frame++) {
				rejected[frame] = percentile_clipping(stack[frame], args->sig, median, crej);
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
				sigma = gsl_stats_ushort_sd(stack, 1, N);
				median = quickmedian (stack, N);
				for (frame = 0; frame < N; frame++) {
					rejected[frame] = sigma_clipping(stack[frame], args->sig, sigma, median, crej);
					if (rejected[frame])
						r++;
					if (N - r <= 4) break;
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
				sigma = gsl_stats_ushort_sd(stack, 1, N);
				median = quickmedian (stack, N);
				n = 0;
				for (frame = 0; frame < N; frame++) {
					if (sigma_clipping(stack[frame], args->sig, sigma, median, crej)) {
						stack[frame] = median;
						n++;
					}
				}
			} while (n > 0 && N > 3);
			break;
		case WINSORIZED:
			do {
				double sigma0;
				sigma = gsl_stats_ushort_sd(stack, 1, N);
				median = quickmedian (stack, N);
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
					rejected[frame] = sigma_clipping(
							stack[frame], args->sig, sigma,
							median, crej);
					if (rejected[frame] != 0)
						r++;
					if (N - r <= 4) break;

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
				double a, b, cov00, cov01, cov11, sumsq;
				quicksort_s(stack, N);
				for (frame = 0; frame < N; frame++) {
					data->xf[frame] = (double)frame;
					data->yf[frame] = (double)stack[frame];
				}
				gsl_fit_linear(data->xf, 1, data->yf, 1, N, &b, &a, &cov00, &cov01, &cov11, &sumsq);
				sigma = 0.0;
				for (frame = 0; frame < N; frame++)
					sigma += (fabs((double)stack[frame] - (a*(double)frame + b)));
				sigma /= (double)N;
				for (frame = 0; frame < N; frame++) {
					rejected[frame] =
						line_clipping(stack[frame], args->sig, sigma, frame, a, b, crej);
					if (rejected[frame] != 0)
						r++;
					if (N - r <= 4) break;
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
