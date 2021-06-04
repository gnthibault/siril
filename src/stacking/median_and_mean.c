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

#include <string.h>
#include <math.h>
#include <gsl/gsl_statistics_ushort.h>
#include <gsl/gsl_cdf.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/siril_date.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/image_format_fits.h"
#include "gui/progress_and_log.h"
#include "algos/sorting.h"
#include "algos/statistics.h"
#include "stacking/stacking.h"
#include "stacking/siril_fit_linear.h"

typedef struct {
	GDateTime *date_obs;
	gdouble exposure;
} DateEvent;

static void free_list_date(gpointer data) {
	DateEvent *item = data;

	g_date_time_unref(item->date_obs);
	g_slice_free(DateEvent, item);
}

static DateEvent* new_item(GDateTime *dt, gdouble exposure) {
	DateEvent *item;

	item = g_slice_new(DateEvent);
	item->exposure = exposure;
	item->date_obs = dt;

	return item;
}

static gint list_date_compare(gconstpointer *a, gconstpointer *b) {
	const DateEvent *dt1 = (const DateEvent *) a;
	const DateEvent *dt2 = (const DateEvent *) b;

	return g_date_time_compare(dt1->date_obs, dt2->date_obs);
}

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

int stack_open_all_files(struct stacking_args *args, int *bitpix, int *naxis, long *naxes, GList **list_date, fits *fit) {
	char msg[256], filename[256];
	int status, oldbitpix = 0, oldnaxis = -1, nb_frames = args->nb_images_to_stack;
	long oldnaxes[3] = { 0 };

	if (args->seq->type == SEQ_REGULAR) {
		for (int i = 0; i < nb_frames; ++i) {
			int image_index = args->image_indices[i];	// image index in sequence
			if (!get_thread_run()) {
				return ST_GENERIC_ERROR;
			}
			if (!fit_sequence_get_image_filename(args->seq, image_index, filename, TRUE))
				continue;

			snprintf(msg, 255, _("Opening image %s for stacking"), filename);
			msg[255] = '\0';
			set_progress_bar_data(msg, PROGRESS_NONE);

			/* open input images */
			if (seq_open_image(args->seq, image_index)) {
				return ST_SEQUENCE_ERROR;
			}

			/* here we use the internal data of sequences, it's quite ugly, we should
			 * consider moving these tests in seq_open_image() or wrapping them in a
			 * sequence function */
			status = 0;
			fits_get_img_param(args->seq->fptr[image_index], 3, bitpix, naxis, naxes, &status);
			if (status) {
				fits_report_error(stderr, status); /* print error message */
				return ST_SEQUENCE_ERROR;
			}
			if (*naxis > 3) {
				siril_log_message(_("Stacking error: images with > 3 dimensions "
						"are not supported\n"));
				return ST_SEQUENCE_ERROR;
			}

			if (oldnaxis > 0) {
				if (*naxis != oldnaxis ||
						oldnaxes[0] != naxes[0] ||
						oldnaxes[1] != naxes[1] ||
						oldnaxes[2] != naxes[2]) {
					siril_log_message(_("Stacking error: input images have "
							"different sizes\n"));
					return ST_SEQUENCE_ERROR;
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
					return ST_SEQUENCE_ERROR;
				}
			} else {
				oldbitpix = *bitpix;
			}

			gdouble current_exp;
			GDateTime *dt = NULL;

			get_date_data_from_fitsfile(args->seq->fptr[image_index], &dt, &current_exp);
			if (dt)
				*list_date = g_list_prepend(*list_date, new_item(dt, current_exp));

			/* We copy metadata from reference to the final fit */
			if (image_index == args->ref_image)
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
		if (!com.pref.debayer.open_debayer && type_ser != SER_RGB && type_ser != SER_BGR)
			type_ser = SER_MONO;
		naxes[2] = type_ser == SER_MONO ? 1 : 3;
		*naxis = type_ser == SER_MONO ? 2 : 3;
		/* case of Super Pixel not handled yet */
		if (com.pref.debayer.open_debayer && com.pref.debayer.bayer_inter == BAYER_SUPER_PIXEL) {
			siril_log_message(_("Super-pixel is not handled yet for on the fly SER stacking\n"));
			return ST_GENERIC_ERROR;
		}

		import_metadata_from_serfile(args->seq->ser_file, fit);
		for (int frame = 0; frame < args->seq->number; frame++) {
			GDateTime *dt = ser_read_frame_date(args->seq->ser_file, frame);
			if (dt)
				*list_date = g_list_prepend(*list_date,	new_item(dt, 0.0));
		}
	}
	else if (args->seq->type == SEQ_FITSEQ) {
		g_assert(args->seq->fitseq_file);
		memcpy(naxes, args->seq->fitseq_file->naxes, sizeof args->seq->fitseq_file->naxes);
		*naxis = naxes[2] == 3 ? 3 : 2;
		*bitpix = args->seq->fitseq_file->bitpix;

		for (int frame = 0; frame < args->seq->number; frame++) {
			if (fitseq_set_current_frame(args->seq->fitseq_file, frame)) {
				siril_log_color_message(_("There was an error opening frame %d for stacking\n"), "red", frame);
				return ST_SEQUENCE_ERROR;
			}
			gdouble current_exp;
			GDateTime *dt = NULL;

			get_date_data_from_fitsfile(args->seq->fitseq_file->fptr, &dt, &current_exp);
			if (dt)
				*list_date = g_list_prepend(*list_date, new_item(dt, current_exp));

			/* We copy metadata from reference to the final fit */
			if (frame == args->ref_image)
				import_metadata_from_fitsfile(args->seq->fitseq_file->fptr, fit);
		}
	} else {
		siril_log_message(_("Rejection stacking is only supported for FITS images/sequences and SER sequences.\nUse \"Sum Stacking\" instead.\n"));
		return ST_SEQUENCE_ERROR;
	}

	return ST_OK;
}

/* The number of blocks must be divisible by the number of channels or they won't be
 * nearly the same size. It must also be divisible by the number of threads, or possibly
 * be close to being when there are many. This favors memory over threads, but since they
 * all start reading at the same time their execution will likely shift, so it may be
 * better to use all available memory.
 */
static int refine_blocks_candidate(int nb_threads, int nb_channels, int minimum_blocks) {
	// we assume that minimum_blocks, the candidate, is at least equal to the number
	// of threads
	int factor_of = nb_channels;
	if (nb_threads < 4) {
		// only allow factors of nb_threads
		if (factor_of != 1 && nb_threads % factor_of == 0)
			factor_of = nb_threads;
		else factor_of *= nb_threads;
		return round_to_ceiling_multiple(minimum_blocks, factor_of);
	}
	// allow 1 minus the factor for 4 - 7 threads
	// allow 3 minus the factor for 8 and more threads
	int minus_factor_allowed = nb_threads < 8 ? 1 : 3;
	int candidate = round_to_ceiling_multiple(minimum_blocks, factor_of);
	do {
		int rem = candidate % nb_threads;
		if (rem == 0 || rem >= (nb_threads - minus_factor_allowed))
			return candidate;
		candidate += factor_of;
	} while (1);
	return candidate;
}

/* median or mean stacking require that the value of each pixel from all images
 * is available. We cannot load all images in memory and it's too slow to read
 * one pixel at a time in all images, so we prepare blocks.
 * Blocks are a part of a channel that will be read from all images, and the
 * stacking will be done on each pixel of the block before going to the next.
 * Parallelization of work is done by assigning some blocks to a thread, which
 * it will process sequentially.
 * To improve load distribution, blocks should be small enough to allow all
 * threads to work but as big as possible for the available memory.
 */
int stack_compute_parallel_blocks(struct _image_block **blocksptr, long max_number_of_rows,
		long naxes[3], int nb_threads, long *largest_block_height, int *nb_blocks) {
	int candidate = nb_threads;	// candidate number of blocks
	while ((max_number_of_rows * candidate) / nb_threads < naxes[1] * naxes[2])
		candidate++;
	candidate = refine_blocks_candidate(nb_threads, (naxes[2] == 3L) ? 3 : 1, candidate);

	*nb_blocks = candidate;
	long height_of_blocks = naxes[1] * naxes[2] / candidate;
	int remainder = naxes[1] % (candidate / naxes[2]);
	siril_log_message(_("We have %d parallel blocks of size %d (+%d) for stacking.\n"),
			*nb_blocks, height_of_blocks, remainder);

	*largest_block_height = 0;
	long channel = 0, row = 0, end, j = 0;
	*blocksptr = malloc(*nb_blocks * sizeof(struct _image_block));
	if (!*blocksptr) {
		PRINT_ALLOC_ERR;
		return ST_ALLOC_ERROR;
	}
	struct _image_block *blocks = *blocksptr;
	do {
		if (j >= *nb_blocks) {
			siril_log_message(_("A bug has been found. Unable to split the image "
						"area into the correct processing blocks.\n"));
			return ST_GENERIC_ERROR;
		}

		blocks[j].channel = channel;
		blocks[j].start_row = row;
		end = row + height_of_blocks - 1;
		if (remainder > 0) {
			// just add one pixel from the remainder to the first blocks to
			// avoid having all of them in the last block
			end++;
			remainder--;
		}
		if (end >= naxes[1] - 1 ||	// end of the line
				(naxes[1] - end < height_of_blocks / 10)) { // not far from it
			end = naxes[1] - 1;
			row = 0;
			channel++;
			remainder = naxes[1] - (*nb_blocks / naxes[2] * height_of_blocks);
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

	} while (channel < naxes[2]) ;

	return ST_OK;
}

static void stack_read_block_data(struct stacking_args *args, int use_regdata,
		struct _image_block *my_block, struct _data_block *data,
		long *naxes, data_type itype, int thread_id) {

	int ielem_size = itype == DATA_FLOAT ? sizeof(float) : sizeof(WORD);
	/* store the layer info to retrieve normalization coeffs*/
	data->layer = (int)my_block->channel;
	/* Read the block from all images, store them in pix[image] */
	for (int frame = 0; frame < args->nb_images_to_stack; ++frame) {
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
					args->image_indices[frame], buffer, &area, thread_id);
			if (retval) {
#ifdef _OPENMP
				int tid = omp_get_thread_num();
				if (tid == 0)
#endif
					siril_log_color_message(_("Error reading one of the image areas\n"), "red");
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

static void norm_to_0_1_range(fits *fit) {
	float mini = fit->fdata[0];
	float maxi = fit->fdata[0];
	long n = fit->naxes[0] * fit->naxes[1] * fit->naxes[2];

	/* search for min / max */
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static) reduction(max:maxi) reduction(min:mini)
#endif
	for (int i = 1; i < n; i++) {
		float tmp = fit->fdata[i];
		if (tmp < mini)
			mini = tmp;
		if (tmp > maxi)
			maxi = tmp;
	}
	/* normalize to [0, 1] range */
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (int i = 0; i < n; i++) {
		fit->fdata[i] = (fit->fdata[i] - mini) / (maxi - mini);
	}
}

/******************************* REJECTION STACKING ******************************
 * The functions below are those managing the rejection, the stacking code is
 * after and similar to median but takes into account the registration data and
 * does a different operation to keep the final pixel values.
 *********************************************************************************/

static int percentile_clipping(WORD pixel, float sig[], float median, guint64 rej[]) {
	float plow = sig[0];
	float phigh = sig[1];

	if ((median - (float) pixel) / median > plow) {
		rej[0]++;
		return -1;
	}
	else if (((float)pixel - median) / median > phigh) {
		rej[1]++;
		return 1;
	}
	return 0;
}

/* Rejection of pixels, following sigma_(high/low) * sigma.
 * The function returns 0 if no rejections are required, 1 if it's a high
 * rejection and -1 for a low-rejection */
static int sigma_clipping(WORD pixel, float sig[], float sigma, float median, guint64 rej[]) {
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

static void Winsorize(WORD *pixel, WORD m0, WORD m1, int N) {
	for (int j = 0; j < N; ++j) {
		pixel[j] = pixel[j] < m0 ? m0 : pixel[j];
		pixel[j] = pixel[j] > m1 ? m1 : pixel[j];
	}
}

static int line_clipping(WORD pixel, float sig[], float sigma, int i, float a, float b, guint64 rej[]) {
	float sigmalow = sig[0];
	float sigmahigh = sig[1];

	if (a * i + b - pixel > sigma * sigmalow) {
		rej[0]++;
		return -1;
	} else if (pixel - a * i - b > sigma * sigmahigh) {
		rej[1]++;
		return 1;
	}
	return 0;
}

static void remove_element(WORD *array, int index, int array_length) {
	for (int i = index; i < array_length - 1; i++)
		array[i] = array[i + 1];
}

static float siril_stats_ushort_sd(const WORD data[], const int N, float *m) {
	double accumulator = 0.0; // accumulating in double precision is important for accuracy
	for (int i = 0; i < N; ++i) {
		accumulator += data[i];
	}
	float mean = (float)(accumulator / N);
	accumulator = 0.0;
	for (int i = 0; i < N; ++i)
		accumulator += (data[i] - mean) * (data[i] - mean);

	if (m) *m = mean;

	return sqrtf((float)(accumulator / (N - 1)));
}

static void grubbs_stat(WORD *stack, int N, float *GCal, int *max_ind) {
	float avg_y;
	float sd = siril_stats_ushort_sd(stack, N, &avg_y);

	/* data are sorted */
	float max_of_deviations = avg_y - stack[0];
	float md2 = stack[N - 1] - avg_y;

	if (md2 > max_of_deviations) {
		max_of_deviations = md2;
		*max_ind = N - 1;
	} else {
		*max_ind = 0;
	}
	*GCal = max_of_deviations / sd;
}

int check_G_values(float Gs, float Gc) {
	return (Gs > Gc);
}

void confirm_outliers(struct outliers *out, int N, double median, int *rejected,
		guint64 rej[2]) {
	int i = N - 1;

	while (!out[i].out && i >= 0) {
		i--;
	}
	for (int j = i; j >= 0; j--) {
		out[j].out = 1;
		if (out[j].x >= median) {
			rejected[out[j].i] = 1;
			rej[1]++;
		} else {
			rejected[out[j].i] = -1;
			rej[0]++;
		}
	}
}

static int apply_rejection_ushort(struct _data_block *data, int nb_frames, struct stacking_args *args, guint64 crej[2]) {
	int N = nb_frames;	// N is the number of pixels kept from the current stack
	float median = 0.f;
	int pixel, output, changed, n, r = 0;
	int firstloop = 1;
	int kept = 0, removed = 0;
	WORD *stack = (WORD *)data->stack;
	WORD *w_stack = (WORD *)data->w_stack;
	int *rejected = (int *)data->rejected;
	WORD *o_stack = (WORD *)data->o_stack;

	memcpy(o_stack, stack, N * sizeof(WORD)); /* making a copy of unsorted stack to apply weights (before the median sorts in place)*/

	/* remove null pixels */
	for (int frame = 0; frame < N; frame++) {
		if (stack[frame] > 0) {
			if (frame != kept) {
				stack[kept] = stack[frame];
			}
			kept++;
		}
	}
	/* Preventing problems
	   0: should not happen but just in case.
	   1 or 2: no need to reject */
	if (kept <= 2) {
		return kept;
	}
	removed = N - kept;
	N = kept;

	/* prepare median and check that the stack is not mostly zero */
	switch (args->type_of_rejection) {
		case PERCENTILE:
		case SIGMA:
		case SIGMEDIAN:
		case WINSORIZED:
			median = quickmedian(stack, N);
			if (median == 0.f)
				return 0;
			break;
		default:
			break;
	}

	switch (args->type_of_rejection) {
		case PERCENTILE:
			for (int frame = 0; frame < N; frame++) {
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
				const float sigma = args->sd_calculator(stack, N);
				if (!firstloop) {
					median = quickmedian(stack, N);
				} else {
					firstloop = 0;
				}
				for (int frame = 0; frame < N; frame++) {
					if (N - r <= 4) {
					// no more rejections
						rejected[frame] = 0;
					} else {
						rejected[frame] = sigma_clipping(stack[frame], args->sig, sigma, median, crej);
						if (rejected[frame]) {
							r++;
						}
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
				const float sigma = args->sd_calculator(stack, N);
				if (!firstloop) {
					median = quickmedian (stack, N);
				} else {
					firstloop = 0;
				}
				n = 0;
				for (int frame = 0; frame < N; frame++) {
					if (sigma_clipping(stack[frame], args->sig, sigma, median, crej)) {
						stack[frame] = median;
						n++;
					}
				}
			} while (n > 0);
			break;
		case WINSORIZED:
			do {
				float sigma0;
				float sigma = args->sd_calculator(stack, N);
				if (!firstloop)
					median = quickmedian (stack, N);
				else firstloop = 0;
				memcpy(w_stack, stack, N * sizeof(WORD));
				do {
					Winsorize(w_stack, roundf_to_WORD(median - 1.5f * sigma), roundf_to_WORD(median + 1.5f * sigma), N);
					sigma0 = sigma;
					sigma = 1.134f * args->sd_calculator(w_stack, N);
				} while (fabs(sigma - sigma0) > sigma0 * 0.0005f);
				for (int frame = 0; frame < N; frame++) {
					if (N - r <= 4) {
						// no more rejections
						rejected[frame] = 0;
					} else {
						rejected[frame] = sigma_clipping(stack[frame], args->sig, sigma, median, crej);
						if (rejected[frame] != 0)
							r++;
					}

				}
				for (pixel = 0, output = 0; pixel < N; pixel++) {
					if (!rejected[pixel]) {
						// copy only if there was a rejection
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
				quicksort_s(stack, N);
				for (int frame = 0; frame < N; frame++) {
					data->yf[frame] = (float)stack[frame];
				}
				float a, b;
				siril_fit_linear(data->xf, data->yf, data->m_x, data->m_dx2, N, &b, &a);
				float sigma = 0.f;
				for (int frame = 0; frame < N; frame++)
					sigma += fabsf(stack[frame] - (a * frame + b));
				sigma /= (float)N;
				for (int frame = 0; frame < N; frame++) {
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
		case GESDT:
			/* Normaly The algorithm does not need to play with sorted data.
			 * But our implementation (after the rejection) needs to be sorted.
			 * So we do it, and by the way we get the median value. Indeed, by design
			 * this algorithm does not have low and high representation of rejection.
			 * We define:
			 * - cold pixel: rejected < median
			 * - hot pixel: rejected > median
			 */

			quicksort_s(stack, N);
			median = gsl_stats_ushort_median_from_sorted_data(stack, 1, N);

			int max_outliers = (int) nb_frames * args->sig[0];
			if (removed >= max_outliers) { /* more than max allowable have been already been removed, should not reject anymore*/
				return kept;
			}
			max_outliers -= removed;
			struct outliers *out = malloc(max_outliers * sizeof(struct outliers));

			memcpy(w_stack, stack, N * sizeof(WORD));
			memset(rejected, 0, N * sizeof(int));

			for (int iter = 0, size = N; iter < max_outliers; iter++, size--) {
				float Gstat;
				int max_index = 0;

				grubbs_stat(w_stack, size, &Gstat, &max_index);
				out[iter].out = check_G_values(Gstat, args->critical_value[iter + removed]);
				out[iter].x = w_stack[max_index];
				out[iter].i = max_index;
				remove_element(w_stack, max_index, size);
			}
			confirm_outliers(out, max_outliers, median, rejected, crej);
			free(out);

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
		default:
		case NO_REJEC:
			;		// Nothing to do, no rejection
	}
	return N;
}

static double mean_and_reject(struct stacking_args *args, struct _data_block *data,
		int stack_size, data_type itype, guint64 crej[2]) {
	double mean;

	int layer = data->layer;
	if (itype == DATA_USHORT) {
		int kept_pixels = apply_rejection_ushort(data, stack_size, args, crej);
		if (kept_pixels == 0)
			mean = quickmedian(data->stack, stack_size);
		else {
			if (args->apply_weight) {
				double *pweights = args->weights + layer * stack_size;
				WORD pmin = 65535, pmax = 0; /* min and max computed here instead of rejection step to avoid dealing with too many particular cases */
				for (int frame = 0; frame < kept_pixels; ++frame) {
					if (pmin > ((WORD*)data->stack)[frame]) pmin = ((WORD*)data->stack)[frame];
					if (pmax < ((WORD*)data->stack)[frame]) pmax = ((WORD*)data->stack)[frame];
				}
				double sum = 0.0;
				double norm = 0.0;
				WORD val;

				for (int frame = 0; frame < stack_size; ++frame) {
					val = ((WORD*)data->o_stack)[frame];
					if (val >= pmin && val <= pmax && val > 0) {
						sum += (float)val * pweights[frame];
						norm += pweights[frame];
					}
				}
				mean = sum / norm;
			} else {
				gint64 sum = 0L;
				for (int frame = 0; frame < kept_pixels; ++frame) {
					sum += ((WORD *)data->stack)[frame];
				}
				mean = sum / (double)kept_pixels;
			}
		}
	} else {
		int kept_pixels = apply_rejection_float(data, stack_size, args, crej);
		if (kept_pixels == 0)
			mean = quickmedian_float(data->stack, stack_size);
		else {
			if (args->apply_weight) {
				double *pweights = args->weights + layer * stack_size;
				float pmin = 10000.0, pmax = -10000.0; /* min and max computed here instead of rejection step to avoid dealing with too many particular cases */
				for (int frame = 0; frame < kept_pixels; ++frame) {
					if (pmin > ((float*)data->stack)[frame]) pmin = ((float*)data->stack)[frame];
					if (pmax < ((float*)data->stack)[frame]) pmax = ((float*)data->stack)[frame];
				}
				double sum = 0.0;
				double norm = 0.0;
				float val;
				for (int frame = 0; frame < stack_size; ++frame) {
					val = ((float*)data->o_stack)[frame];
					if (val >= pmin && val <= pmax && val != 0.f) {
						sum += val * pweights[frame];
						norm += pweights[frame];
					}
				}
			mean = sum / norm;
			} else {
				double sum = 0.0;
				for (int frame = 0; frame < kept_pixels; ++frame) {
					sum += ((float*)data->stack)[frame];
				}
				mean = sum / (double)kept_pixels;
			}
		}
	}
	return mean;
}

int stack_mean_with_rejection(struct stacking_args *args) {
	return stack_mean_or_median(args, TRUE);
}

int stack_median(struct stacking_args *args) {
	return stack_mean_or_median(args, FALSE);
}

static int compute_weights(struct stacking_args *args) {
	int nb_frames = args->nb_images_to_stack;
	int nb_layers = args->seq->nb_layers;
	int idx = 0;

	args->weights = malloc(nb_layers * nb_frames * sizeof(double));
	double *pweights[3];

	for (int layer = 0; layer < nb_layers; ++layer) {
		double norm = 0.0;
		pweights[layer] = args->weights + layer * nb_frames;
		for (int i = 0; i < args->nb_images_to_stack; ++i) {
			idx = args->image_indices[i];
			pweights[layer][i] = 1.f /
				(args->coeff.pscale[layer][i] * args->coeff.pscale[layer][i] *
				 args->seq->stats[layer][idx]->bgnoise * args->seq->stats[layer][idx]->bgnoise);
			norm += pweights[layer][i];
		}
		norm /= (double) nb_frames;

		for (int i = 0; i < args->nb_images_to_stack; i++) {
			pweights[layer][i] /= norm;
		}
	}
	return ST_OK;
}

static void compute_date_time_keywords(GList *list_date, fits *fit) {
	if (list_date) {
		gdouble exposure = 0.0;
		GDateTime *date_obs;
		gdouble start, end;
		GList *list;
		/* First we want to sort the list */
		list_date = g_list_sort(list_date, (GCompareFunc) list_date_compare);
		/* Then we compute the sum of exposure */
		for (list = list_date; list; list = list->next) {
			exposure += ((DateEvent *)list->data)->exposure;
		}

		/* go to the first stacked image and get needed values */
		list_date = g_list_first(list_date);
		date_obs = g_date_time_ref(((DateEvent *)list_date->data)->date_obs);
		start = date_time_to_Julian(((DateEvent *)list_date->data)->date_obs);

		/* go to the last stacked image and get needed values
		 * This time we need to add the exposure to the date_obs
		 * to exactly retrieve the end of the exposure
		 */
		list_date = g_list_last(list_date);
		gdouble last_exp = ((DateEvent *)list_date->data)->exposure;
		GDateTime *last_date = ((DateEvent *)list_date->data)->date_obs;
		GDateTime *corrected_last_date = g_date_time_add_seconds(last_date, (gdouble) last_exp);

		end = date_time_to_Julian(corrected_last_date);

		g_date_time_unref(corrected_last_date);

		/* we address the computed values to the keywords */
		fit->exposure = exposure;
		fit->date_obs = date_obs;
		fit->expstart = start;
		fit->expend = end;
	}
}

/* How many rows fit in memory, based on image size, number and available memory.
 * It returns at most the total number of rows of the image (naxes[1] * naxes[2]) */
static long stack_get_max_number_of_rows(long naxes[3], data_type type, int nb_images_to_stack) {
	int max_memory = get_max_memory_in_MB();
	long total_nb_rows = naxes[1] * naxes[2];

	siril_log_message(_("Using %d MB memory maximum for stacking\n"), max_memory);
	int elem_size = type == DATA_FLOAT ? sizeof(float) : sizeof(WORD);
	guint64 number_of_rows = (guint64)max_memory * BYTES_IN_A_MB /
		((guint64)naxes[0] * nb_images_to_stack * elem_size);
	// this is how many rows we can load in parallel from all images of the
	// sequence and be under the limit defined in config in megabytes.
	if (total_nb_rows < number_of_rows)
		return total_nb_rows;
	return (long)number_of_rows;
}

static int stack_mean_or_median(struct stacking_args *args, gboolean is_mean) {
	int bitpix, i, naxis, cur_nb = 0, retval = ST_OK, pool_size = 1;
	long naxes[3];
	struct _data_block *data_pool = NULL;
	struct _image_block *blocks = NULL;
	fits fit = { 0 }; // output result
	fits ref = { 0 }; // reference image, used to get metadata back
	// data for mean/rej only
	guint64 irej[3][2] = {{0,0}, {0,0}, {0,0}};
	regdata *layerparam = NULL;
	gboolean use_regdata = is_mean;


	int nb_frames = args->nb_images_to_stack; // number of frames actually used
	naxes[0] = naxes[1] = 0; naxes[2] = 1;

	if (nb_frames < 2) {
		siril_log_message(_("Select at least two frames for stacking. Aborting.\n"));
		return ST_GENERIC_ERROR;
	} else if (nb_frames < 3 && is_mean && args->type_of_rejection == GESDT) {
		siril_log_message(_("The Generalized Extreme Studentized Deviate Test needs at least three frames for stacking. Aborting.\n"));
		return ST_GENERIC_ERROR;
	}
	g_assert(nb_frames <= args->seq->number);

	if (use_regdata) {
		if (args->reglayer < 0) {
			siril_log_message(_("No registration layer passed, ignoring registration data!\n"));
			use_regdata = FALSE;
		}
		else layerparam = args->seq->regparam[args->reglayer];
	}

	set_progress_bar_data(NULL, PROGRESS_RESET);

	/* first loop: open all fits files and check they are of same size */
	GList *list_date = NULL;
	if ((retval = stack_open_all_files(args, &bitpix, &naxis, naxes, &list_date, &ref))) {
		goto free_and_close;
	}

	if (naxes[0] == 0) {
		// no image has been loaded
		siril_log_color_message(_("Rejection stack error: uninitialized sequence\n"), "red");
		retval = ST_SEQUENCE_ERROR;
		goto free_and_close;
	}
	if (naxes[0] != args->seq->rx || naxes[1] != args->seq->ry) {
		siril_log_color_message(_("Rejection stack error: sequence has wrong image size (%dx%d for sequence, %ldx%ld for images)\n"), "red", args->seq->rx, args->seq->ry, naxes[0], naxes[1]);
		retval = ST_SEQUENCE_ERROR;
		goto free_and_close;
	}
	if (sequence_is_rgb(args->seq) && naxes[2] != 3) {
		siril_log_message(_("Processing the sequence as RGB\n"));
		naxes[2] = 3;
	}
	fprintf(stdout, "image size: %ldx%ld, %ld layers\n", naxes[0], naxes[1], naxes[2]);

	/* initialize result image */
	fits *fptr = &fit;
	if ((retval = new_fit_image(&fptr, naxes[0], naxes[1], naxes[2],
					args->use_32bit_output ? DATA_FLOAT : DATA_USHORT))) {
		goto free_and_close;
	}
	copy_fits_metadata(&ref, fptr);
	clearfits(&ref);
	if (!args->use_32bit_output && (args->output_norm || fit.orig_bitpix != BYTE_IMG)) {
		fit.bitpix = USHORT_IMG;
		if (args->output_norm)
			fit.orig_bitpix = USHORT_IMG;
	}

	/* manage threads */
	int nb_threads;
#ifdef _OPENMP
	nb_threads = com.max_thread;
	if (nb_threads > 1 && (args->seq->type == SEQ_REGULAR || args->seq->type == SEQ_FITSEQ)) {
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
#ifdef HAVE_FFMS2
	if (args->seq->type == SEQ_AVI) {
		siril_log_color_message(_("Stacking a film will work only on one core and will be slower than if you convert it to SER\n"), "salmon");
		nb_threads = 1;
	}
#endif // HAVE_FFMS2
#else
	nb_threads = 1;
#endif

	/* manage memory */
	long largest_block_height;
	int nb_blocks;
	data_type itype = get_data_type(bitpix);
	long max_number_of_rows = stack_get_max_number_of_rows(naxes, itype, args->nb_images_to_stack);
	/* Compute parallel processing data: the data blocks, later distributed to threads */
	if ((retval = stack_compute_parallel_blocks(&blocks, max_number_of_rows, naxes, nb_threads,
					&largest_block_height, &nb_blocks))) {
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
	size_t npixels_in_block = largest_block_height * naxes[0];
	g_assert(npixels_in_block > 0);
	int ielem_size = itype == DATA_FLOAT ? sizeof(float) : sizeof(WORD);

	fprintf(stdout, "allocating data for %d threads (each %'lu MB)\n", pool_size,
			(unsigned long)(nb_frames * npixels_in_block * ielem_size) / BYTES_IN_A_MB);
	data_pool = calloc(pool_size, sizeof(struct _data_block));
	size_t bufferSize = ielem_size * nb_frames * (npixels_in_block + 1ul) + 4ul; // buffer for tmp and stack, added 4 byte for alignment
	if (is_mean) {
		bufferSize += nb_frames * sizeof(int); // for rejected
		bufferSize += ielem_size * nb_frames; // for o_stack
		if (args->type_of_rejection == WINSORIZED) {
			bufferSize += ielem_size * nb_frames; // for w_frame
		} else if (args->type_of_rejection == GESDT) {
			bufferSize += ielem_size * nb_frames; // for w_frame
			bufferSize += sizeof(float) * (int) floor(nb_frames * args->sig[0]); //and GCritical
		} else if (args->type_of_rejection == LINEARFIT) {
			bufferSize += 2 * sizeof(float) * nb_frames; // for xc and yc
		}
	}
	for (i = 0; i < pool_size; i++) {
		data_pool[i].pix = malloc(nb_frames * sizeof(void *));
		data_pool[i].tmp = malloc(bufferSize);
		if (!data_pool[i].pix || !data_pool[i].tmp) {
			PRINT_ALLOC_ERR;
			gchar *available = g_format_size_full(get_available_memory(), G_FORMAT_SIZE_IEC_UNITS);
			fprintf(stderr, "Cannot allocate %zu (free memory: %s)\n", bufferSize / BYTES_IN_A_MB, available);
			fprintf(stderr, "CHANGE MEMORY SETTINGS if stacking takes too much.\n");

			g_free(available);
			retval = ST_ALLOC_ERROR;
			goto free_and_close;
		}
		data_pool[i].stack = (void*) ((char*) data_pool[i].tmp
				+ nb_frames * npixels_in_block * ielem_size);
		if (is_mean) {
			size_t offset = (size_t)ielem_size * nb_frames * (npixels_in_block + 1);
			int temp = offset % sizeof(int);
			if (temp > 0) { // align buffer
				offset += sizeof(int) - temp;
			}
			data_pool[i].rejected = (int*)((char*)data_pool[i].tmp + offset);
			data_pool[i].o_stack = (void*)((char*)data_pool[i].rejected + sizeof(int) * nb_frames);

			if (args->type_of_rejection == WINSORIZED) {
				data_pool[i].w_stack = (void*)((char*)data_pool[i].o_stack + ielem_size * nb_frames);
			} else if (args->type_of_rejection == GESDT) {
				data_pool[i].w_stack = (void*)((char*)data_pool[i].o_stack + ielem_size * nb_frames);
				int max_outliers = (int) floor(nb_frames * args->sig[0]);
				args->critical_value = malloc(max_outliers * sizeof(float));
				for (int j = 0, size = nb_frames; j < max_outliers; j++, size--) {
					float t_dist = gsl_cdf_tdist_Pinv(1 - args->sig[1] / (2 * size), size - 2);
					float numerator = (size - 1) * t_dist;
					float denominator = sqrtf(size) * sqrtf(size - 2 + (t_dist * t_dist));
					args->critical_value[j] = numerator / denominator;
				}
			} else if (args->type_of_rejection == LINEARFIT) {
				data_pool[i].xf = (float (*)) ((char*)data_pool[i].o_stack + ielem_size * nb_frames);
				data_pool[i].yf = data_pool[i].xf + nb_frames;
				// precalculate some stuff
				data_pool[i].m_x = (nb_frames - 1) * 0.5f;
				data_pool[i].m_dx2 = 0.f;
				for (int j = 0; j < nb_frames; ++j) {
					const float dx = j - data_pool[i].m_x;
					data_pool[i].xf[j] = 1.f / (j + 1);
					data_pool[i].m_dx2 += (dx * dx - data_pool[i].m_dx2)
							* data_pool[i].xf[j];
				}
				data_pool[i].m_dx2 = 1.f / data_pool[i].m_dx2;
			}
		}

		for (int j = 0; j < nb_frames; ++j) {
			if (itype == DATA_FLOAT)
				data_pool[i].pix[j] = ((float*)data_pool[i].tmp) + j * npixels_in_block;
			else 	data_pool[i].pix[j] = ((WORD *)data_pool[i].tmp) + j * npixels_in_block;
		}
	}

	if (itype == DATA_USHORT)
		args->sd_calculator = nb_frames < 65536 ? siril_stats_ushort_sd_32 : siril_stats_ushort_sd_64;

	if (args->apply_weight) {
		siril_log_message(_("Computing weights...\n"));
		retval = compute_weights(args);
		if(retval) {
			retval = ST_GENERIC_ERROR;
			goto free_and_close;
		}
	}

	siril_log_message(_("Starting stacking...\n"));
	if (is_mean)
		set_progress_bar_data(_("Rejection stacking in progress..."), PROGRESS_RESET);
	else	set_progress_bar_data(_("Median stacking in progress..."), PROGRESS_RESET);
	double total = (double)(naxes[2] * naxes[1] + 2); // for progress bar

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

		if (!get_thread_run()) retval = ST_GENERIC_ERROR;
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
		stack_read_block_data(args, use_regdata, my_block, data, naxes, itype, data_idx);

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
		int layer = my_block->channel;
		for (y = 0; y < my_block->height; y++)
		{
			/* index of the pixel in the result image
			 * we read line y, but we need to store it at
			 * ry - y - 1 to not have the image mirrored. */
			size_t pdata_idx = (naxes[1] - (my_block->start_row + y) - 1) * naxes[0];
			/* index of the line in the read data, data->pix[frame] */
			size_t line_idx = y * naxes[0];
			guint64 crej[2] = {0, 0};
			if (retval) break;

			// update progress bar
#ifdef _OPENMP
#pragma omp atomic
#endif
			cur_nb++;

			if (!get_thread_run()) {
				retval = ST_GENERIC_ERROR;
				break;
			}
			if (!(cur_nb % 16))	// every 16 iterations
				set_progress_bar_data(NULL, (double)cur_nb/total);

			for (x = 0; x < naxes[0]; ++x) {
				/* copy all images pixel values in the same row array `stack'
				 * to optimize caching and improve readability */
				for (int frame = 0; frame < nb_frames; ++frame) {
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
							else ((WORD *)data->stack)[frame] = 0;
							continue;
						}

						pix_idx -= shiftx;
					}

					WORD pixel = 0; float fpixel = 0.f;
					if (itype == DATA_FLOAT)
						fpixel = ((float*) data->pix[frame])[pix_idx];
					else
						pixel = ((WORD*) data->pix[frame])[pix_idx];
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
								if (fpixel != 0.f) { // do not normalize null pixels to detect them later
									tmp = fpixel * args->coeff.pscale[layer][frame];
									((float*)data->stack)[frame] = (float)(tmp - args->coeff.poffset[layer][frame]);
								} else {
									((float*)data->stack)[frame] = 0.f;
								}
							} else {
								if (pixel > 0) { // do not normalize null pixels to detect them later
									tmp = (double)pixel * args->coeff.pscale[layer][frame];
									((WORD *)data->stack)[frame] = round_to_WORD(tmp - args->coeff.poffset[layer][frame]);
								} else {
									((WORD *)data->stack)[frame] = 0;
								}
							}
							break;
						case MULTIPLICATIVE:
							// multiplicative  (scale[frame] = 1, offset[frame] = 0)
						case MULTIPLICATIVE_SCALING:
							// multiplicative + scale (offset[frame] = 0)
							if (itype == DATA_FLOAT) {
								tmp = fpixel * args->coeff.pscale[layer][frame];
								((float*)data->stack)[frame] = (float)(tmp * args->coeff.pmul[layer][frame]);
							} else {
								tmp = (double)pixel * args->coeff.pscale[layer][frame];
								((WORD *)data->stack)[frame] = round_to_WORD(tmp * args->coeff.pmul[layer][frame]);
							}
							break;
					}
				}

				double result; // resulting pixel value, either mean or median
				if (is_mean) {
					result = mean_and_reject(args, data, nb_frames, itype, crej);
				} else {
					if (itype == DATA_USHORT)
						result = quickmedian(data->stack, nb_frames);
					else 	result = quickmedian_float(data->stack, nb_frames);
				}

				if (args->use_32bit_output) {
					if (itype == DATA_USHORT)
						fit.fpdata[my_block->channel][pdata_idx] = min(double_ushort_to_float_range(result), 1.f);
					else	fit.fpdata[my_block->channel][pdata_idx] = min((float)result, 1.f);
				} else {
					/* in case of 8bit data we may want to normalize to 16bits */
					if (args->output_norm) {
						normalize_to16bit(bitpix, &result);
					}
					fit.pdata[my_block->channel][pdata_idx] = round_to_WORD(result);
				}
				pdata_idx++;
			} // end of for x

			if (is_mean && args->type_of_rejection != NO_REJEC) {
#ifdef _OPENMP
#pragma omp atomic
#endif
				irej[my_block->channel][0] += crej[0];
#ifdef _OPENMP
#pragma omp atomic
#endif
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
	if (is_mean) {
		double nb_tot = (double) naxes[0] * (double) naxes[1] * (double) nb_frames;
		for (long channel = 0; channel < naxes[2]; channel++) {
			siril_log_message(_("Pixel rejection in channel #%d: %.3lf%% - %.3lf%%\n"),
					channel, (double) irej[channel][0] / nb_tot * 100.0,
					(double) irej[channel][1] / nb_tot * 100.0);
		}
	}

	/* copy result to gfit if success */
	clearfits(&gfit);
	copyfits(&fit, &gfit, CP_FORMAT, -1);

	if (args->use_32bit_output) {
		gfit.fdata = fit.fdata;
		for (i = 0; i < fit.naxes[2]; i++)
			gfit.fpdata[i] = fit.fpdata[i];

		if (args->output_norm) {
			norm_to_0_1_range(&gfit);
		}
	} else {
		gfit.data = fit.data;
		for (i = 0; i < fit.naxes[2]; i++)
			gfit.pdata[i] = fit.pdata[i];
	}

	compute_date_time_keywords(list_date, &gfit);

free_and_close:
	fprintf(stdout, "free and close (%d)\n", retval);
	for (i = 0; i < nb_frames; ++i) {
		seq_close_image(args->seq, args->image_indices[i]);
	}

	if (data_pool) {
		for (i=0; i<pool_size; i++) {
			if (data_pool[i].pix) free(data_pool[i].pix);
			if (data_pool[i].tmp) free(data_pool[i].tmp);
		}
		free(data_pool);
	}
	g_list_free_full(list_date, (GDestroyNotify) free_list_date);
	if (blocks) free(blocks);
	if (args->normalize) {
		free(args->coeff.offset);
		free(args->coeff.scale);
		free(args->coeff.mul);
	}

	if (args->weights) free(args->weights);
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
			siril_log_message(_("Median stacking complete. %d images have been stacked.\n"), nb_frames);
		}
	}
	return retval;
}

