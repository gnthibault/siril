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

#include "core/siril.h"
#include "core/proto.h"
#include "gui/progress_and_log.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "registration/registration.h"

#include "stacking/stacking.h"

static int stack_addminmax(struct stacking_args *args, gboolean ismax);


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
	WORD *final_pixel[3];
	float *ffinal_pixel[3];
	double exposure = 0.0;
	gboolean is_float = TRUE; // init only for warning
	size_t nbdata = 0;
	char filename[256];
	int retval = ST_OK, nb_frames, cur_nb = 0;
	fits fit = { 0 };

	/* should be pre-computed to display it in the stacking tab */
	nb_frames = args->nb_images_to_stack;
	int reglayer = get_registration_layer(args->seq);

	if (nb_frames <= 1) {
		siril_log_message(_("No frame selected for stacking (select at least 2). Aborting.\n"));
		return ST_GENERIC_ERROR;
	}

	final_pixel[0] = NULL;
	ffinal_pixel[0] = NULL;
	g_assert(args->seq->nb_layers == 1 || args->seq->nb_layers == 3);
	g_assert(nb_frames <= args->seq->number);

	for (int j = 0; j < args->seq->number; ++j) {
		if (!get_thread_run()) {
			retval = ST_GENERIC_ERROR;
			goto free_and_reset_progress_bar;
		}
		if (!args->filtering_criterion(args->seq, j, args->filtering_parameter)) {
			fprintf(stdout, "image %d is excluded from stacking\n", j);
			continue;
		}
		if (!seq_get_image_filename(args->seq, j, filename)) {
			retval = ST_GENERIC_ERROR;
			goto free_and_reset_progress_bar;
		}
		gchar *tmpmsg = g_strdup_printf(_("Processing image %s"), filename);
		set_progress_bar_data(tmpmsg, (double) cur_nb / ((double) nb_frames + 1.));
		g_free(tmpmsg);

		cur_nb++;	// only used for progress bar

		if (seq_read_frame(args->seq, j, &fit, FALSE, -1)) {
			siril_log_message(_("Stacking: could not read frame, aborting\n"));
			retval = ST_SEQUENCE_ERROR;
			goto free_and_reset_progress_bar;
		}

		g_assert(args->seq->nb_layers == 1 || args->seq->nb_layers == 3);
		g_assert(fit.naxes[2] == args->seq->nb_layers);

		/* first loaded image: init data structures for stacking */
		if (!nbdata) {
			is_float = fit.type == DATA_FLOAT;
			nbdata = fit.naxes[0] * fit.naxes[1];
			if (is_float) {
				if (ismax)
					ffinal_pixel[0] = calloc(nbdata * fit.naxes[2], sizeof(float));
				else {
					ffinal_pixel[0] = malloc(nbdata * fit.naxes[2] * sizeof(float));
					for (long k = 0; k < nbdata * fit.naxes[2]; k++)
						ffinal_pixel[0][k] = 1.0;
				}
				if (!ffinal_pixel[0]) {
					PRINT_ALLOC_ERR;
					retval = ST_ALLOC_ERROR;
					goto free_and_reset_progress_bar;
				}
				if (args->seq->nb_layers == 3) {
					ffinal_pixel[1] = ffinal_pixel[0] + nbdata;
					ffinal_pixel[2] = ffinal_pixel[1] + nbdata;
				}
			} else {
				if (ismax)
					final_pixel[0] = calloc(nbdata * fit.naxes[2], sizeof(WORD));
				else {
					final_pixel[0] = malloc(nbdata * fit.naxes[2] * sizeof(WORD));
					for (long k = 0; k < nbdata * fit.naxes[2]; k++)
						final_pixel[0][k] = USHRT_MAX;
				}
				if (!final_pixel[0]) {
					PRINT_ALLOC_ERR;
					retval = ST_ALLOC_ERROR;
					goto free_and_reset_progress_bar;
				}
				if (args->seq->nb_layers == 3) {
					final_pixel[1] = final_pixel[0] + nbdata;
					final_pixel[2] = final_pixel[1] + nbdata;
				}
			}
		} else if (fit.ry * fit.rx != nbdata) {
			siril_log_message(_("Stacking: image in sequence doesn't has the same dimensions\n"));
			retval = ST_SEQUENCE_ERROR;
			goto free_and_reset_progress_bar;
		}

		/* load registration data for current image */
		int shiftx, shifty;
		if (reglayer != -1 && args->seq->regparam[reglayer]) {
			shiftx = round_to_int(args->seq->regparam[reglayer][j].shiftx * (float) args->seq->upscale_at_stacking);
			shifty = round_to_int(args->seq->regparam[reglayer][j].shifty * (float) args->seq->upscale_at_stacking);
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
		size_t i = 0;	// index in final_pixel[0]
		for (int y = 0; y < fit.ry; ++y) {
			for (int x = 0; x < fit.rx; ++x) {
				int nx = x - shiftx;
				int ny = y - shifty;
				//printf("x=%d y=%d sx=%d sy=%d i=%d ii=%d\n",x,y,shiftx,shifty,i,ii);
				if (nx >= 0 && nx < fit.rx && ny >= 0 && ny < fit.ry) {
					size_t ii = ny * fit.rx + nx;		// index in final_pixel[0] too
					//printf("shiftx=%d shifty=%d i=%d ii=%d\n",shiftx,shifty,i,ii);
					if (ii > 0 && ii < fit.rx * fit.ry) {
						for (int layer = 0; layer < args->seq->nb_layers; ++layer) {
							if (is_float) {
								float current_pixel = fit.fpdata[layer][ii];
								// we take the brightest pixel
								if ((ismax && current_pixel > ffinal_pixel[layer][i]) ||
										// we take the darkest pixel
										(!ismax && current_pixel < ffinal_pixel[layer][i]))
									ffinal_pixel[layer][i] = current_pixel;
							} else {
								WORD current_pixel = fit.pdata[layer][ii];
								// we take the brightest pixel
								if ((ismax && current_pixel > final_pixel[layer][i]) ||
										// we take the darkest pixel
										(!ismax && current_pixel < final_pixel[layer][i]))
									final_pixel[layer][i] = current_pixel;
							}
						}
					}
				}
				++i;
			}
		}
		clearfits(&fit);
	}
	if (!get_thread_run()) {
		retval = ST_GENERIC_ERROR;
		goto free_and_reset_progress_bar;
	}
	set_progress_bar_data(_("Finalizing stacking..."), (double) nb_frames / ((double) nb_frames + 1.));

	clearfits(&gfit);
	fits *result = &gfit;
	if (is_float) {
		if (new_fit_image_with_data(&result, args->seq->rx, args->seq->ry, args->seq->nb_layers, DATA_FLOAT, ffinal_pixel[0]))
			return ST_GENERIC_ERROR;
	} else {
		if (new_fit_image_with_data(&result, args->seq->rx, args->seq->ry, args->seq->nb_layers, DATA_USHORT, final_pixel[0]))
			return ST_GENERIC_ERROR;
	}

	/* We copy metadata from reference to the final fit */
	if (args->seq->type == SEQ_REGULAR) {
		int ref = args->ref_image;
		if (!seq_open_image(args->seq, ref)) {
			import_metadata_from_fitsfile(args->seq->fptr[ref], &gfit);
			seq_close_image(args->seq, ref);
		}
	}

free_and_reset_progress_bar:
	if (retval) {
		set_progress_bar_data(_("Stacking failed. Check the log."), PROGRESS_RESET);
		siril_log_message(_("Stacking failed.\n"));
	} else {
		set_progress_bar_data(_("Stacking complete."), PROGRESS_DONE);
	}

	return retval;
}
