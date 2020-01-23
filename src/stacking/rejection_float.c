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
#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics_float.h>
#include <stdint.h>

#include "core/siril.h"
#include "stacking.h"
#include "algos/sorting.h"

static int percentile_clipping(float pixel, double sig[], double median, uint64_t rej[]) {
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
	return 0;
}

/* Rejection of pixels, following sigma_(high/low) * sigma.
 * The function returns 0 if no rejections are required, 1 if it's a high
 * rejection and -1 for a low-rejection */
static int sigma_clipping(float pixel, double sig[], double sigma, double median, uint64_t rej[]) {
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
	return 0;
}

static void Winsorize(float *pixel, double m0, double m1) {
	if (*pixel < m0) *pixel = (float)m0;
	else if (*pixel > m1) *pixel = (float)m1;
}

static int line_clipping(float pixel, double sig[], double sigma, int i, double a, double b, uint64_t rej[]) {
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
	return 0;
}

int apply_rejection_float(struct _data_block *data, int nb_frames, struct stacking_args *args, uint64_t crej[2]) {
	int N = nb_frames;	// N is the number of pixels kept from the current stack
	double median, sigma = -1.0;
	int frame, pixel, output, changed, n, r = 0;
	int firstloop = 1;

	float *stack = (float *)data->stack;
	float *w_stack = (float *)data->w_stack;
	float *rejected = (float *)data->rejected;

	/* prepare median and check that the stack is not mostly zero */
	switch (args->type_of_rejection) {
		case PERCENTILE:
		case SIGMA:
		case SIGMEDIAN:
		case WINSORIZED:
			median = quickmedian_float (stack, N);
			if (median == 0.0)
				return 0;
			break;
		default:
			break;
	}

	switch (args->type_of_rejection) {
		case PERCENTILE:
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
				sigma = gsl_stats_float_sd(stack, 1, N);
				if (!firstloop)
					median = quickmedian_float (stack, N);
				else firstloop = 0;
				for (frame = 0; frame < N; frame++) {
					if (N - r <= 4) {
						// no more rejections
						rejected[frame] = 0;
					} else {
						rejected[frame] = sigma_clipping(stack[frame], args->sig, sigma, median, crej);
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
				sigma = gsl_stats_float_sd(stack, 1, N);
				if (!firstloop)
					median = quickmedian_float (stack, N);
				else firstloop = 0;
				n = 0;
				for (frame = 0; frame < N; frame++) {
					if (sigma_clipping(stack[frame], args->sig, sigma, median, crej)) {
						stack[frame] = median;
						n++;
					}
				}
			} while (n > 0);
			break;
		case WINSORIZED:
			do {
				double sigma0;
				sigma = gsl_stats_float_sd(stack, 1, N);
				if (!firstloop)
					median = quickmedian_float (stack, N);
				else firstloop = 0;
				memcpy(w_stack, stack, N * sizeof(float));
				do {
					int jj;
					double m0 = median - 1.5 * sigma;
					double m1 = median + 1.5 * sigma;
					for (jj = 0; jj < N; jj++)
						Winsorize(w_stack+jj, m0, m1);
					median = quickmedian_float (w_stack, N);
					sigma0 = sigma;
					sigma = 1.134 * gsl_stats_float_sd(w_stack, 1, N);
				} while ((fabs(sigma - sigma0) / sigma0) > 0.0005);
				for (frame = 0; frame < N; frame++) {
					if (N - r <= 4) {
						// no more rejections
						rejected[frame] = 0;
					} else {
						rejected[frame] = sigma_clipping(
								stack[frame], args->sig, sigma,
								median, crej);
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
				double a, b, cov00, cov01, cov11, sumsq;
				quicksort_f(stack, N);
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


