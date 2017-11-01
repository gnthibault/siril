/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2017 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/processing.h"
#include "core/proto.h"		// FITS functions
#include "stacking.h"

struct sum_stacking_data {
	unsigned long *sum[3];	// the new image's channels
	double exposure;	// sum of the exposures
	int reglayer;		// layer used for registration data
};

int sum_stacking_prepare_hook(struct generic_seq_args *args) {
	struct sum_stacking_data *ssdata = args->user;
	unsigned int nbdata = args->seq->ry * args->seq->rx;
	ssdata->sum[0] = calloc(nbdata, sizeof(unsigned long)*args->seq->nb_layers);
	if (ssdata->sum[0] == NULL){
		printf("Stacking: memory allocation failure\n");
		return -1;
	}
	if(args->seq->nb_layers == 3){
		ssdata->sum[1] = ssdata->sum[0] + nbdata;	// index of green layer in sum[0]
		ssdata->sum[2] = ssdata->sum[0] + nbdata*2;	// index of blue layer in sum[0]
	} else {
		ssdata->sum[1] = NULL;
		ssdata->sum[2] = NULL;
	}

	ssdata->exposure = 0.0;
	return 0;
}

int sum_stacking_image_hook(struct generic_seq_args *args, int i, fits *fit, rectangle *_) {
	struct sum_stacking_data *ssdata = args->user;
	int shiftx, shifty, nx, ny, x, y, ii, layer;
	int pixel = 0;	// index in sum[0]

#ifdef _OPENMP
#pragma omp atomic
#endif
	ssdata->exposure += fit->exposure;
	
	if (ssdata->reglayer != -1 && args->seq->regparam[ssdata->reglayer]) {
		shiftx = args->seq->regparam[ssdata->reglayer][i].shiftx;
		shifty = args->seq->regparam[ssdata->reglayer][i].shifty;
	} else {
		shiftx = 0;
		shifty = 0;
	}

	for (y=0; y < fit->ry; ++y){
		for (x=0; x < fit->rx; ++x){
			nx = x - shiftx;
			ny = y - shifty;
			if (nx >= 0 && nx < fit->rx && ny >= 0 && ny < fit->ry) {
				// we have data for this pixel
				ii = ny * fit->rx + nx;		// index in source image
				if (ii > 0 && ii < fit->rx * fit->ry){
					for(layer=0; layer<args->seq->nb_layers; ++layer) {
#ifdef _OPENMP
#pragma omp atomic
#endif
						ssdata->sum[layer][pixel] += fit->pdata[layer][ii];
					}
				}
			}
			++pixel;
		}
	}
	return 0;
}

// convert the result and store it into gfit
int sum_stacking_finalize_hook(struct generic_seq_args *args) {
	struct sum_stacking_data *ssdata = args->user;
	unsigned long max = 0L;	// max value of the image's channels
	unsigned int i, nbdata;
	int layer;

	nbdata = args->seq->ry * args->seq->rx * args->seq->nb_layers;
	// find the max first
	#pragma omp parallel for reduction(max:max)
	for (i=0; i < nbdata; ++i)
		if (ssdata->sum[0][i] > max)
			max = ssdata->sum[0][i];

	clearfits(&gfit);
	new_fit_image(&gfit, args->seq->rx, args->seq->ry, args->seq->nb_layers);
	gfit.hi = round_to_WORD(max);
	gfit.bitpix = USHORT_IMG;
	gfit.exposure = ssdata->exposure;

	double ratio = 1.0;
	if (max > USHRT_MAX)
		ratio = USHRT_MAX_DOUBLE / (double)max;

	nbdata = args->seq->ry * args->seq->rx;
	for (layer=0; layer<args->seq->nb_layers; ++layer){
		unsigned long* from = ssdata->sum[layer];
		WORD *to = gfit.pdata[layer];
		for (i=0; i < nbdata; ++i) {
			if (ratio == 1.0)
				*to++ = round_to_WORD(*from++);
			else	*to++ = round_to_WORD((double)(*from++) * ratio);
		}
	}

	free(ssdata->sum[0]);
	return 0;
}


int stack_summing_generic(struct stacking_args *stackargs) {
	struct generic_seq_args *args = malloc(sizeof(struct generic_seq_args));
	args->seq = stackargs->seq;
	args->partial_image = FALSE;
	args->filtering_criterion = stackargs->filtering_criterion;
	args->filtering_parameter = stackargs->filtering_parameter;
	args->nb_filtered_images = stackargs->nb_images_to_stack;
	args->prepare_hook = sum_stacking_prepare_hook;
	args->image_hook = sum_stacking_image_hook;
	args->save_hook = NULL;
	args->finalize_hook = sum_stacking_finalize_hook;
	args->idle_function = NULL;
	args->description = "sum stacking";
	args->has_output = FALSE;
	args->already_in_a_thread = TRUE;
	args->parallel = TRUE;

	struct sum_stacking_data *ssdata = malloc(sizeof(struct sum_stacking_data));
	ssdata->reglayer = stackargs->reglayer;
	args->user = ssdata;

	//start_in_new_thread(generic_sequence_worker, args);
	generic_sequence_worker(args);
	return args->retval;
}

