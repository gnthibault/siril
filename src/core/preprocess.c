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
#include <math.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "algos/statistics.h"
#include "filters/cosmetic_correction.h"
#include "gui/progress_and_log.h"
#include "io/sequence.h"
#include "io/conversion.h"
#include "io/ser.h"

#include "preprocess.h"

static double evaluateNoiseOfCalibratedImage(fits *fit, fits *dark, double k) {
	double noise = 0.0;
	fits dark_tmp = { 0 }, fit_tmp = { 0 };
	int chan, ret = 0;
	rectangle area = { 0 };

	/* square of 512x512 in the center of the image */
	int size = 512;
	area.x = (fit->rx - size) / 2;
	area.y = (fit->ry - size) / 2;
	area.w = size;
	area.h = size;

	copyfits(dark, &dark_tmp, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	copyfits(fit, &fit_tmp, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);

	soper(&dark_tmp, k, OPER_MUL);
	ret = imoper(&fit_tmp, &dark_tmp, OPER_SUB);
	if (ret) {
		clearfits(&dark_tmp);
		clearfits(&fit_tmp);
		return -1.0;
	}

	for (chan = 0; chan < fit->naxes[2]; chan++) {
		imstats *stat = statistics(NULL, -1, &fit_tmp, chan, &area, STATS_BASIC);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 0.0;	// not -1?
		}
		noise += stat->sigma;
		free_stats(stat);
	}
	clearfits(&dark_tmp);
	clearfits(&fit_tmp);

	return noise;
}

#define GR ((sqrt(5.0) - 1.0) / 2.0)

static double goldenSectionSearch(fits *raw, fits *dark, double a, double b,
		double tol) {
	double c, d;
	double fc, fd;
	int iter = 0;

	c = b - GR * (b - a);
	d = a + GR * (b - a);
	// TODO: check return values of evaluateNoiseOfCalibratedImage
	fc = evaluateNoiseOfCalibratedImage(raw, dark, c);
	fd = evaluateNoiseOfCalibratedImage(raw, dark, d);
	do {
		siril_debug_print("Iter: %d (%1.2lf, %1.2lf)\n", ++iter, c, d);
		if (fc < 0.0 || fd < 0.0)
			return -1.0;
		if (fc < fd) {
			b = d;
			d = c;
			fd = fc;
			c = b - GR * (b - a);
			fc = evaluateNoiseOfCalibratedImage(raw, dark, c);
		} else {
			a = c;
			c = d;
			fc = fd;
			d = a + GR * (b - a);
			fd = evaluateNoiseOfCalibratedImage(raw, dark, d);

		}
	} while (fabs(c - d) > tol);
	return ((b + a) / 2.0);
}


static int preprocess(fits *raw, struct preprocessing_data *args) {
	int ret = 0;
	if (args->use_bias) {
		ret = imoper(raw, args->bias, OPER_SUB);
	}

	/* if dark optimization, the master-dark has already been subtracted */
	if (!ret && args->use_dark && !args->use_dark_optim) {
		ret = imoper(raw, args->dark, OPER_SUB);
	}

	if (!ret && args->use_flat) {
		// return value is an error if there is an overflow, but it is usual
		// for now, so don't treat as error
		/*ret =*/ siril_fdiv(raw, args->flat, args->normalisation);
	}

	return ret;
}

static int darkOptimization(fits *raw, fits *dark, fits *offset) {
	double k0;
	double lo = 0.0, up = 2.0;
	int ret = 0;
	fits dark_tmp = { 0 };

	if (raw->rx != dark->rx || raw->ry != dark->ry) {
		fprintf(stderr, "image size mismatch\n");
		return -1;
	}

	if (copyfits(dark, &dark_tmp, CP_ALLOC | CP_COPYA | CP_FORMAT, 0))
		return 1;

	/* Minimization of background noise to find better k */
	invalidate_stats_from_fit(raw);
	k0 = goldenSectionSearch(raw, &dark_tmp, lo, up, 1E-3);
	if (k0 < 0.0) {
		ret = -1;
	} else {
		siril_log_message(_("Dark optimization: k0=%.3lf\n"), k0);
		/* Multiply coefficient to master-dark */
		soper(&dark_tmp, k0, OPER_MUL);
		ret = imoper(raw, &dark_tmp, OPER_SUB);
	}
	clearfits(&dark_tmp);
	return ret;
}

static int prepro_prepare_hook(struct generic_seq_args *args) {
	struct preprocessing_data *prepro = args->user;

	if (prepro->seq) {
		// checking disk space: removing old sequence and computing free space
		remove_prefixed_sequence_files(args->seq, prepro->ppprefix);

		int64_t size = seq_compute_size(args->seq, args->seq->number);
		if (prepro->debayer)
			size *= 3;
		if (test_available_space(size))
			return 1;

		// handling SER
		if (ser_prepare_hook(args))
			return 1;
	}

	// precompute flat levels
	if (prepro->use_flat) {
		if (prepro->equalize_cfa) {
			compute_grey_flat(prepro->flat);
		}
		if (prepro->autolevel) {
			imstats *stat = statistics(NULL, -1, prepro->flat, RLAYER, NULL, STATS_BASIC);
			if (!stat) {
				siril_log_message(_("Error: statistics computation failed.\n"));
				return 1;
			}
			prepro->normalisation = stat->mean;
			siril_log_message(_("Normalisation value auto evaluated: %.2lf\n"),
					prepro->normalisation);
			free_stats(stat);
		}
	}

	// proceed to cosmetic correction
	if (prepro->use_cosmetic_correction && prepro->use_dark) {
		if (prepro->dark->naxes[2] == 1) {
			prepro->dev = find_deviant_pixels(prepro->dark, prepro->sigma,
					&(prepro->icold), &(prepro->ihot));
			siril_log_message(_("%ld pixels corrected (%ld + %ld)\n"),
					prepro->icold + prepro->ihot, prepro->icold, prepro->ihot);
		} else
			siril_log_message(_("Darkmap cosmetic correction "
						"is only supported with single channel images\n"));
	}
	return 0;
}

static int prepro_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_) {
	struct preprocessing_data *prepro = args->user;
	if (prepro->use_dark_optim && prepro->use_dark) {
		if (darkOptimization(fit, prepro->dark, prepro->bias))
			return 1;
	}

	if (preprocess(fit, prepro))
		return 1;

	if (prepro->use_cosmetic_correction && prepro->use_dark && prepro->dark->naxes[2] == 1)
		cosmeticCorrection(fit, prepro->dev, prepro->icold + prepro->ihot, prepro->is_cfa);

	if (prepro->debayer) {
		if (!prepro->seq || prepro->seq->type == SEQ_REGULAR) {
			// not for SER because it is done on-the-fly
			debayer_if_needed(TYPEFITS, fit,
					prepro->compatibility, TRUE, prepro->stretch_cfa);
		}
	}

	return 0;
}

static void clear_preprocessing_data(struct preprocessing_data *prepro) {
	if (prepro->use_bias && prepro->bias)
		clearfits(prepro->bias);
	if (prepro->use_dark && prepro->dark)
		clearfits(prepro->dark);
	if (prepro->use_flat && prepro->flat)
		clearfits(prepro->flat);
}

static int prepro_finalize_hook(struct generic_seq_args *args) {
	int retval = ser_finalize_hook(args);
	struct preprocessing_data *prepro = args->user;
	clear_preprocessing_data(prepro);
	free(args->user);
	return retval;
}

gpointer prepro_worker(gpointer p) {
	gpointer retval = generic_sequence_worker(p);

	struct generic_seq_args *args = (struct generic_seq_args *)p;
	free_sequence(args->seq, TRUE);
	free(args);
	return retval;
}

void start_sequence_preprocessing(struct preprocessing_data *prepro, gboolean from_script) {
	struct generic_seq_args *args = malloc(sizeof(struct generic_seq_args));
	args->seq = prepro->seq;
	args->partial_image = FALSE;
	args->filtering_criterion = seq_filter_all;
	args->nb_filtered_images = prepro->seq->number;
	args->prepare_hook = prepro_prepare_hook;
	args->image_hook = prepro_image_hook;
	args->save_hook = NULL;
	args->finalize_hook = prepro_finalize_hook;
	args->idle_function = NULL;
	args->stop_on_error = TRUE;
	args->description = _("Preprocessing");
	args->has_output = TRUE;
	args->new_seq_prefix = prepro->ppprefix;
	args->load_new_sequence = TRUE;
	args->force_ser_output = FALSE;
	args->parallel = TRUE;
	args->user = prepro;

	if (from_script) {
		args->already_in_a_thread = TRUE;
		start_in_new_thread(prepro_worker, args);
	} else {
		args->already_in_a_thread = FALSE;
		start_in_new_thread(generic_sequence_worker, args);
	}
}

/********** SINGLE IMAGE ************/
int preprocess_single_image(struct preprocessing_data *args) {
	gchar *dest_filename, *msg;
	fits fit = { 0 };
	int ret = 0;

	msg = g_strdup_printf(_("Pre-processing image %s"), com.uniq->filename);
	set_progress_bar_data(msg, 0.5);
	g_free(msg);
	struct generic_seq_args generic = { .user = args };

	copyfits(com.uniq->fit, &fit, CP_ALLOC | CP_FORMAT | CP_COPYA, 0);
	copy_fits_metadata(com.uniq->fit, &fit);

	ret = prepro_prepare_hook(&generic);
	if (!ret)
		ret = prepro_image_hook(&generic, 0, 0, &fit, NULL);
	clear_preprocessing_data(args);

	if (!ret) {
		gchar *filename = g_path_get_basename(com.uniq->filename);
		char *filename_noext = remove_ext_from_filename(filename);
		g_free(filename);
		dest_filename = g_strdup_printf("%s%s%s", args->ppprefix, filename_noext, com.ext);
		msg = g_strdup_printf(_("Saving image %s"), filename_noext);
		set_progress_bar_data(msg, PROGRESS_NONE);
		ret = savefits(dest_filename, &fit);

		if (!ret) {
			// open the new image?
			copyfits(&fit, com.uniq->fit, CP_ALLOC | CP_FORMAT | CP_COPYA, 0);
			if (com.uniq->nb_layers != fit.naxes[2]) {
				com.uniq->nb_layers = fit.naxes[2];
				com.uniq->layers = realloc(com.uniq->layers, com.uniq->nb_layers * sizeof(layer_info));
			}
			if (com.uniq->filename)
				free(com.uniq->filename);
			com.uniq->filename = strdup(dest_filename);
		}

		clearfits(&fit);
		free(filename_noext);
		g_free(dest_filename);
		g_free(msg);
	}

	return ret;
}
