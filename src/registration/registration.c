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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <gtk/gtk.h>

#include "core/siril.h"
#include "core/proto.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/message_dialog.h"
#include "gui/plot.h"
#include "gui/progress_and_log.h"
#include "gui/sequence_list.h"
#include "core/initfile.h"
#include "core/OS_utils.h"
#include "registration/registration.h"
#include "registration/matching/misc.h"
#include "registration/matching/match.h"
#include "registration/matching/atpmatch.h"
#include "stacking/stacking.h"
#include "algos/PSF.h"
#include "gui/PSF_list.h"
#include "algos/quality.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/single_image.h"
#include "io/image_format_fits.h"
#include "opencv/opencv.h"
#include "opencv/ecc/ecc.h"

#undef DEBUG

static char *tooltip_text[] = { N_("<b>One Star Registration</b>: This is the simplest method to register deep-sky images. "
		"Because only one star is concerned for register, images are aligned using shifting "
		"(at a fraction of pixel). No rotation or scaling are performed. "
		"Shifts at pixel precision are saved in seq file."),
		N_("<b>Two or Three Stars Registration</b>: This method looks like the one star registration except one need to select "
		"two or three stars. This is very useful for field with a few stars."),
		N_("<b>Global Star Alignment</b>: This is a more powerful and accurate algorithm (but also slower) "
		"to perform deep-sky images. The global matching is based on triangle similarity method for automatically "
		"identify common stars in each image. "
		"A new sequence is created with the prefix of your choice (r_ by default)."),
		N_("<b>Image Pattern Alignment</b>: This is a simple registration by translation method "
		"using cross correlation in the spatial domain. This method is fast and is used to register "
		"planetary movies. It can also be used for some deep-sky images registration. "
		"Shifts at pixel precision are saved in seq file."),
		N_("<b>Enhanced Correlation Coefficient Maximization</b>: It is based on the enhanced correlation "
		"coefficient maximization algorithm. This method is more complex and slower than Image Pattern Alignment "
		"but no selection is required. It is good for moon surface images registration. Only translation is taken "
		"into account yet."),
		N_("<b>Comet/Asteroid Registration</b>: This algorithm is dedicated to the comet and asteroid registration. It is necessary to have timestamps "
		"stored in FITS header and to load a sequence of star aligned images. This methods makes a translation of a certain number of pixels depending on "
		"the timestamp of each images and the global shift of the object between the first and the last image.")
};
/* callback for the selected area event */
void _reg_selected_area_callback() {
	if (!com.headless)
		update_reg_interface(TRUE);
}

static struct registration_method *reg_methods[NUMBER_OF_METHODS];
static gboolean end_register_idle(gpointer p);

struct registration_method *new_reg_method(const char *name, registration_function f,
		selection_type s, registration_type t) {
	struct registration_method *reg = malloc(sizeof(struct registration_method));
	reg->name = strdup(name);
	reg->method_ptr = f;
	reg->sel = s;
	reg->type = t;
	return reg;
}

void initialize_registration_methods() {
	GtkComboBoxText *regcombo;
	int i = 0, j = 0;
	GString *tip;
	gchar *ctip;

	reg_methods[i++] = new_reg_method(_("One Star Registration (deep-sky)"),
			&register_shift_fwhm, REQUIRES_ANY_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[i++] = new_reg_method(_("Two or Three Stars Registration (deep-sky)"),
			&register_3stars, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[i++] = new_reg_method(_("Global Star Alignment (deep-sky)"),
			&register_star_alignment, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[i++] = new_reg_method(_("Image Pattern Alignment (planetary - full disk)"),
			&register_shift_dft, REQUIRES_SQUARED_SELECTION, REGTYPE_PLANETARY);
	reg_methods[i++] = new_reg_method(_("Enhanced Correlation Coefficient (planetary - surfaces)"),
			&register_ecc, REQUIRES_NO_SELECTION, REGTYPE_PLANETARY);
	reg_methods[i++] = new_reg_method(_("Comet/Asteroid Registration"),
			&register_comet, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
	reg_methods[i] = NULL;

	tip = g_string_new ("");
	for (j = 0; j < i; j ++) {
		tip = g_string_append(tip, _(tooltip_text[j]));
		if (j < i - 1)
			tip = g_string_append(tip, "\n\n");
	}
	ctip = g_string_free (tip, FALSE);
	gtk_widget_set_tooltip_markup(lookup_widget("comboboxregmethod"), ctip);
	g_free(ctip);

	/* fill comboboxregmethod */
	regcombo = GTK_COMBO_BOX_TEXT(
			gtk_builder_get_object(builder, "comboboxregmethod"));
	gtk_combo_box_text_remove_all(regcombo);
	i = 0;
	while (reg_methods[i] != NULL) {
		gtk_combo_box_text_append_text(regcombo, reg_methods[i]->name);
		siril_log_message(_("Loading registration method: %s\n"),
				reg_methods[i]->name);
		i++;
	}
	if (i > 0) {
		gtk_combo_box_set_active(GTK_COMBO_BOX(regcombo), com.reg_settings);
	}

	/* register to the new area selected event */
	register_selection_update_callback(_reg_selected_area_callback);
}

struct registration_method *get_selected_registration_method() {
	GtkComboBoxText *regcombo = GTK_COMBO_BOX_TEXT(
			gtk_builder_get_object(builder, "comboboxregmethod"));
	int index = 0;

	gchar *text = gtk_combo_box_text_get_active_text (regcombo);
	while (reg_methods[index] && text != NULL) {
		if (!strcmp(reg_methods[index]->name, text))
			break;
		index++;
	}
	g_free(text);
	return reg_methods[index];
}

static void normalizeQualityData(struct registration_args *args, double q_min, double q_max) {
	int frame;
	double diff = q_max - q_min;

	/* this case occurs when all images but one are excluded */
	if (diff == 0) {
		q_min = 0;
		diff = q_max;
	}

	for (frame = 0; frame < args->seq->number; ++frame) {
		args->seq->regparam[args->layer][frame].quality -= q_min;
		args->seq->regparam[args->layer][frame].quality /= diff;
		/* if thread has been manually stopped, some values will be < 0 */
		if ((args->seq->regparam[args->layer][frame].quality < 0)
				|| isnan(args->seq->regparam[args->layer][frame].quality))
			args->seq->regparam[args->layer][frame].quality = -1.0;
	}
}

/* Calculate shift in images to be aligned with the reference image, using
 * discrete Fourrier transform on a square selected area and matching the
 * phases.
 */
int register_shift_dft(struct registration_args *args) {
	fits fit_ref = { 0 };
	unsigned int frame, size, sqsize, j;
	fftwf_complex *ref, *in, *out, *convol;
	fftwf_plan p, q;
	int ret;
	int abort = 0;
	float nb_frames, cur_nb;
	int ref_image;
	regdata *current_regdata;
	double q_max = 0, q_min = DBL_MAX;
	int q_index = -1;

	/* the selection needs to be squared for the DFT */
	assert(args->selection.w == args->selection.h);
	size = args->selection.w;
	sqsize = size * size;

	if (args->process_all_frames)
		nb_frames = (float) args->seq->number;
	else
		nb_frames = (float) args->seq->selnum;

	if (args->seq->regparam[args->layer]) {
		siril_log_message(
				_("Recomputing already existing registration for this layer\n"));
		current_regdata = args->seq->regparam[args->layer];
		/* we reset all values as we may register different images */
		memset(current_regdata, 0, args->seq->number * sizeof(regdata));
	} else {
		current_regdata = calloc(args->seq->number, sizeof(regdata));
		if (current_regdata == NULL) {
			PRINT_ALLOC_ERR;
			return -2;
		}
		args->seq->regparam[args->layer] = current_regdata;
	}

	/* loading reference frame */
	ref_image = sequence_find_refimage(args->seq);

	set_progress_bar_data(
			_("Register DFT: loading and processing reference frame"),
			PROGRESS_NONE);
	ret = seq_read_frame_part(args->seq, args->layer, ref_image, &fit_ref,
			&args->selection, FALSE, -1);

	if (ret) {
		siril_log_message(
				_("Register: could not load first image to register, aborting.\n"));
		args->seq->regparam[args->layer] = NULL;
		free(current_regdata);
		clearfits(&fit_ref);
		return ret;
	}

	ref = fftwf_malloc(sizeof(fftwf_complex) * sqsize);
	in = fftwf_malloc(sizeof(fftwf_complex) * sqsize);
	out = fftwf_malloc(sizeof(fftwf_complex) * sqsize);
	convol = fftwf_malloc(sizeof(fftwf_complex) * sqsize);

	gchar* wisdomFile = g_build_filename(g_get_user_cache_dir(), "siril_fftw.wisdom", NULL);

	// test for available wisdom
	p = fftwf_plan_dft_2d(size, size, ref, out, FFTW_FORWARD, FFTW_WISDOM_ONLY);
	if (!p) {
		// no wisdom available, load wisdom from file
		fftwf_import_wisdom_from_filename(wisdomFile);
		// test again for wisdom
		p = fftwf_plan_dft_2d(size, size, ref, out, FFTW_FORWARD, FFTW_WISDOM_ONLY);
		if (!p) {
			// build plan with FFTW_MEASURE
			p = fftwf_plan_dft_2d(size, size, ref, out, FFTW_FORWARD, FFTW_MEASURE);
			// save the wisdom
			fftwf_export_wisdom_to_filename(wisdomFile);
		}
	}

	q = fftwf_plan_dft_2d(size, size, convol, out, FFTW_BACKWARD, FFTW_WISDOM_ONLY | FFTW_DESTROY_INPUT);
	if (!q) {
		// no wisdom available, load wisdom from file
		fftwf_import_wisdom_from_filename(wisdomFile);
		// test again for wisdom
		q = fftwf_plan_dft_2d(size, size, convol, out, FFTW_BACKWARD, FFTW_WISDOM_ONLY | FFTW_DESTROY_INPUT);
		if (!q) {
			// build plan with FFTW_MEASURE
			q = fftwf_plan_dft_2d(size, size, convol, out, FFTW_BACKWARD, FFTW_MEASURE | FFTW_DESTROY_INPUT);
			// save the wisdom
			fftwf_export_wisdom_to_filename(wisdomFile);
		}
	}
	g_free(wisdomFile);
	fftwf_free(out); // was needed to build the plan, can be freed now
	fftwf_free(convol); // was needed to build the plan, can be freed now

	// copying image selection into the fftw data
	if (fit_ref.type == DATA_USHORT) {
		for (j = 0; j < sqsize; j++)
			ref[j] = (float)fit_ref.data[j];
	}
	else if (fit_ref.type == DATA_FLOAT) {
		for (j = 0; j < sqsize; j++)
			ref[j] = fit_ref.fdata[j];
	}

	// We don't need fit_ref anymore, we can destroy it.
	current_regdata[ref_image].quality = QualityEstimate(&fit_ref, args->layer);
	clearfits(&fit_ref);
	fftwf_execute_dft(p, ref, in); /* repeat as needed */

	fftwf_free(ref); // not needed anymore

	set_shifts(args->seq, ref_image, args->layer, 0.0, 0.0, FALSE);

	q_min = q_max = current_regdata[ref_image].quality;
	q_index = ref_image;

	cur_nb = 0.f;

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(guided) \
	if (args->seq->type == SEQ_SER || ((args->seq->type == SEQ_REGULAR || args->seq->type == SEQ_FITSEQ) && fits_is_reentrant()))
#endif
	for (frame = 0; frame < args->seq->number; ++frame) {
		if (abort) continue;
		if (args->run_in_thread && !get_thread_run()) {
			abort = 1;
			continue;
		}
		if (frame == ref_image) continue;
		if (!args->process_all_frames && !args->seq->imgparam[frame].incl)
			continue;

		fits fit = { 0 };
		char tmpmsg[1024], tmpfilename[256];

		seq_get_image_filename(args->seq, frame, tmpfilename);
		g_snprintf(tmpmsg, 1024, _("Register: processing image %s"), tmpfilename);
		set_progress_bar_data(tmpmsg, PROGRESS_NONE);
		int thread_id = -1;
#ifdef _OPENMP
		thread_id = omp_get_thread_num();
#endif
		if (!seq_read_frame_part(args->seq, args->layer, frame, &fit,
						&args->selection, FALSE, thread_id)) {
			int x;
			fftwf_complex *img = fftwf_malloc(sizeof(fftwf_complex) * sqsize);
			fftwf_complex *out2 = fftwf_malloc(sizeof(fftwf_complex) * sqsize);

			// copying image selection into the fftw data
			if (fit.type == DATA_USHORT) {
				for (x = 0; x < sqsize; x++)
					img[x] = (float)fit.data[x];
			}
			if (fit.type == DATA_FLOAT) {
				for (x = 0; x < sqsize; x++)
					img[x] = fit.fdata[x];
			}

			current_regdata[frame].quality = QualityEstimate(&fit, args->layer);
			// after this call, fit data is dead

#ifdef _OPENMP
#pragma omp critical
#endif
			{
				double qual = current_regdata[frame].quality;
				if (qual > q_max) {
					q_max = qual;
					q_index = frame;
				}
				q_min = min(q_min, qual);
			}

			fftwf_execute_dft(p, img, out2); /* repeat as needed */

			fftwf_complex *convol2 = img; // reuse buffer

			for (x = 0; x < sqsize; x++) {
				convol2[x] = in[x] * conjf(out2[x]);
			}

			fftwf_execute_dft(q, convol2, out2); /* repeat as needed */

			int shift = 0;
			for (x = 1; x < sqsize; ++x) {
				if (crealf(out2[x]) > crealf(out2[shift])) {
					shift = x;
					// break or get last value?
				}
			}
			int shifty = shift / size;
			int shiftx = shift % size;
			if (shifty > size / 2) {
				shifty -= size;
			}
			if (shiftx > size / 2) {
				shiftx -= size;
			}

			set_shifts(args->seq, frame, args->layer, (float)shiftx, (float)shifty,
					fit.top_down);

			// We don't need fit anymore, we can destroy it.
			clearfits(&fit);

			/* shiftx and shifty are the x and y values for translation that
			 * would make this image aligned with the reference image.
			 * WARNING: the y value is counted backwards, since the FITS is
			 * stored down from up.
			 */
#ifdef DEBUG
			fprintf(stderr,
					"reg: frame %d, shiftx=%f shifty=%f quality=%g\n",
					args->seq->imgparam[frame].filenum,
					current_regdata[frame].shiftx, current_regdata[frame].shifty,
					current_regdata[frame].quality);
#endif
#ifdef _OPENMP
#pragma omp atomic
#endif
			cur_nb += 1.f;
			set_progress_bar_data(NULL, cur_nb / nb_frames);
			fftwf_free(img);
			fftwf_free(out2);
		} else {
			//report_fits_error(ret, error_buffer);
			args->seq->regparam[args->layer] = NULL;
			free(current_regdata);
			abort = ret = 1;
			continue;
		}
	}

	fftwf_destroy_plan(p);
	fftwf_destroy_plan(q);
	fftwf_free(in);
	if (!ret) {
		if (args->x2upscale)
			args->seq->upscale_at_stacking = 2.0;
		else
			args->seq->upscale_at_stacking = 1.0;
		normalizeQualityData(args, q_min, q_max);

		siril_log_message(_("Registration finished.\n"));
		siril_log_color_message(_("Best frame: #%d.\n"), "bold", q_index + 1);
	} else {
		free(args->seq->regparam[args->layer]);
		args->seq->regparam[args->layer] = NULL;
	}
	return ret;
}

/* register images: calculate shift in images to be aligned with the reference image;
 * images are not modified, only shift parameters are saved in regparam in the sequence.
 * layer is the layer on which the registration will be done, green by default (set in siril_init())
 */
int register_shift_fwhm(struct registration_args *args) {
	int frame, ref_image;
	float nb_frames, cur_nb = 0.f;
	double reference_xpos, reference_ypos;
	double fwhm_min = DBL_MAX;
	int fwhm_index = -1;
	regdata *current_regdata;

	framing_mode framing = ORIGINAL_FRAME;
	if (args->follow_star)
		framing = FOLLOW_STAR_FRAME;

	/* First and longest step: get the minimization data on one star for all
	 * images to register, which provides FWHM but also star coordinates */
	// TODO: detect that it was already computed, and don't do it again
	// -> should be done at a higher level and passed in the args
	if (seqpsf(args->seq, args->layer, TRUE, args->process_all_frames, framing, FALSE))
		return 1;

	// regparam is managed in seqpsf idle function already
	current_regdata = args->seq->regparam[args->layer];

	if (args->process_all_frames)
		nb_frames = (float)args->seq->number;
	else nb_frames = (float)args->seq->selnum;

	/* loading reference frame */
	ref_image = sequence_find_refimage(args->seq);
	if (!current_regdata[ref_image].fwhm_data) {
		siril_log_message(
				_("Registration PSF: failed to compute PSF for reference frame at least\n"));
		return -1;
	}
	reference_xpos = current_regdata[ref_image].fwhm_data->xpos;
	reference_ypos = current_regdata[ref_image].fwhm_data->ypos;

	fwhm_min = current_regdata[ref_image].fwhm_data->fwhmx;

	fwhm_index = ref_image;

	/* Second step: align image by aligning star coordinates together */
	for (frame = 0; frame < args->seq->number; frame++) {
		double tmp;
		if (args->run_in_thread && !get_thread_run())
			break;
		if (!args->process_all_frames && !args->seq->imgparam[frame].incl)
			continue;
		if (frame == ref_image || !current_regdata[frame].fwhm_data) {
			set_shifts(args->seq, frame, args->layer, 0.0, 0.0, FALSE);
			continue;
		}
		if (current_regdata[frame].fwhm < fwhm_min
				&& current_regdata[frame].fwhm > 0.0) {
			fwhm_min = current_regdata[frame].fwhm;
			fwhm_index = frame;
		}
		tmp = reference_xpos - current_regdata[frame].fwhm_data->xpos;
		current_regdata[frame].shiftx = tmp;
		tmp = current_regdata[frame].fwhm_data->ypos - reference_ypos;
		current_regdata[frame].shifty = tmp;

		fprintf(stderr, "reg: file %d, shiftx=%f shifty=%f\n",
				args->seq->imgparam[frame].filenum,
				current_regdata[frame].shiftx, current_regdata[frame].shifty);
		cur_nb += 1.f;
		set_progress_bar_data(NULL, cur_nb / nb_frames);
	}

	if (args->x2upscale)
		args->seq->upscale_at_stacking = 2.0;
	else
		args->seq->upscale_at_stacking = 1.0;

	siril_log_message(_("Registration finished.\n"));
	siril_log_color_message(_("Best frame: #%d with fwhm=%.3g.\n"), "bold",
			fwhm_index + 1, fwhm_min);
	return 0;
}

int register_ecc(struct registration_args *args) {
	int frame, ref_image, failed = 0;
	float nb_frames, cur_nb;
	regdata *current_regdata;
	fits ref = { 0 };
	double q_max = 0, q_min = DBL_MAX;
	int q_index = -1;
	int abort = 0;

	if (args->seq->regparam[args->layer]) {
		current_regdata = args->seq->regparam[args->layer];
		/* we reset all values as we may register different images */
		memset(current_regdata, 0, args->seq->number * sizeof(regdata));
	} else {
		current_regdata = calloc(args->seq->number, sizeof(regdata));
		if (current_regdata == NULL) {
			PRINT_ALLOC_ERR;
			return -2;
		}
		args->seq->regparam[args->layer] = current_regdata;
	}

	if (args->process_all_frames)
		nb_frames = (float)args->seq->number;
	else nb_frames = (float)args->seq->selnum;

	/* loading reference frame */
	ref_image = sequence_find_refimage(args->seq);

	if (seq_read_frame(args->seq, ref_image, &ref, FALSE, -1)) {
		siril_log_message(_("Could not load reference image\n"));
		args->seq->regparam[args->layer] = NULL;
		free(current_regdata);
		return 1;
	}
	current_regdata[ref_image].quality = QualityEstimate(&ref, args->layer);
	/* we make sure to free data in the destroyed fit */
	clearfits(&ref);
	/* Ugly code: as QualityEstimate destroys fit we need to reload it */
	seq_read_frame(args->seq, ref_image, &ref, FALSE, -1);
	image_find_minmax(&ref);
	q_min = q_max = current_regdata[ref_image].quality;
	q_index = ref_image;

	/* then we compare to other frames */
	if (args->process_all_frames)
		args->new_total = args->seq->number;
	else args->new_total = args->seq->selnum;

	cur_nb = 0.f;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(guided) \
	if (args->seq->type == SEQ_SER || ((args->seq->type == SEQ_REGULAR || args->seq->type == SEQ_FITSEQ) && fits_is_reentrant()))
#endif
	for (frame = 0; frame < args->seq->number; frame++) {
		if (abort) continue;
		if (args->run_in_thread && !get_thread_run()) {
			abort = 1;
			continue;
		}
		if (frame == ref_image) continue;
		if (!args->process_all_frames && !args->seq->imgparam[frame].incl)
			continue;

		set_shifts(args->seq, frame, args->layer, 0.0, 0.0, FALSE);
		fits im = { 0 };
		char tmpmsg[1024], tmpfilename[256];

		seq_get_image_filename(args->seq, frame, tmpfilename);
		g_snprintf(tmpmsg, 1024, _("Register: processing image %s"),
				tmpfilename);
		set_progress_bar_data(tmpmsg, PROGRESS_NONE);
		int thread_id = -1;
#ifdef _OPENMP
		thread_id = omp_get_thread_num();
#endif
		if (!seq_read_frame(args->seq, frame, &im, FALSE, thread_id)) {
			reg_ecc reg_param = { 0 };
			image_find_minmax(&im);

			if (findTransform(&ref, &im, args->layer, &reg_param)) {
				siril_log_message(
						_("Cannot perform ECC alignment for frame %d\n"),
						frame + 1);
				/* We exclude this frame */
				args->seq->imgparam[frame].incl = FALSE;
				current_regdata[frame].quality = 0.0;
				args->seq->selnum--;
#ifdef _OPENMP
#pragma omp atomic
#endif
				++failed;
				clearfits(&im);
				continue;
			}

			current_regdata[frame].quality = QualityEstimate(&im,
					args->layer);
#ifdef _OPENMP
#pragma omp critical
#endif
			{
				double qual = current_regdata[frame].quality;
				if (qual > q_max) {
					q_max = qual;
					q_index = frame;
				}
				q_min = min(q_min, qual);
			}

			set_shifts(args->seq, frame, args->layer, -reg_param.dx,
					-reg_param.dy, im.top_down);
#ifdef _OPENMP
#pragma omp atomic
#endif
			cur_nb += 1.f;
			set_progress_bar_data(NULL, cur_nb / nb_frames);
			clearfits(&im);
		}
	}

	if (args->x2upscale)
		args->seq->upscale_at_stacking = 2.0;
	else
		args->seq->upscale_at_stacking = 1.0;

	normalizeQualityData(args, q_min, q_max);
	clearfits(&ref);

	siril_log_message(_("Registration finished.\n"));
	if (failed) {
		gchar *str = ngettext("%d file was ignored and excluded\n", "%d files were ignored and excluded\n", failed);
		str = g_strdup_printf(str, failed);
		siril_log_color_message(str, "red");
		g_free(str);
	}
	siril_log_color_message(_("Best frame: #%d.\n"), "bold", q_index + 1);

	return 0;
}

void on_comboboxregmethod_changed(GtkComboBox *box, gpointer user_data) {
	int index = 0;
	gchar *text = gtk_combo_box_text_get_active_text (GTK_COMBO_BOX_TEXT(box));

	while (reg_methods[index] && text != NULL) {
		if (!strcmp(reg_methods[index]->name, text))
			break;
		index++;
	}
	g_free(text);

	com.reg_settings = index;
	update_reg_interface(TRUE);
	writeinitfile();
}

void on_comboreg_transfo_changed(GtkComboBox *box, gpointer user_data) {
	GtkAdjustment *register_minpairs = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "register_minpairs"));
	double val = gtk_adjustment_get_value(register_minpairs);

	switch(gtk_combo_box_get_active(box)) {
	case SHIFT_TRANSFORMATION:
	case AFFINE_TRANSFORMATION:
		gtk_adjustment_set_lower (register_minpairs, 3);
		break;
	case HOMOGRAPHY_TRANSFORMATION:
		gtk_adjustment_set_lower (register_minpairs, 4);
		if (val < 4)
			gtk_adjustment_set_value(register_minpairs, 4);
		break;
	default:
		printf("on_comboreg_transfo_changed: Value not handled.\n");
	}
}

/* for now, the sequence argument is used only when executing a script */
int get_registration_layer(sequence *seq) {
	if (!com.script) {
		GtkComboBox *registbox = GTK_COMBO_BOX(lookup_widget("comboboxreglayer"));
		int reglayer = gtk_combo_box_get_active(registbox);
		if (!seq || !seq->regparam || seq->nb_layers < 0 || seq->nb_layers <= reglayer)
			return -1;
		return reglayer;
	} else {
		// find first available regdata
		if (!seq || !seq->regparam || seq->nb_layers < 0)
			return -1;
		int i;
		for (i = 0; i < seq->nb_layers; i++)
			if (seq->regparam[i])
				return i;
		return -1;
	}
}

/* Selects the "register all" or "register selected" according to the number of
 * selected images, if argument is false.
 * Verifies that enough images are selected and an area is selected.
 */
void update_reg_interface(gboolean dont_change_reg_radio) {
	static GtkWidget *go_register = NULL, *follow = NULL, *cumul_data = NULL;
	static GtkLabel *labelreginfo = NULL;
	static GtkToggleButton *reg_all = NULL, *reg_sel = NULL;
	static GtkNotebook *notebook_reg = NULL;
	int nb_images_reg; /* the number of images to register */
	struct registration_method *method;
	gboolean selection_is_done;

	if (!go_register) {
		go_register = lookup_widget("goregister_button");
		follow = lookup_widget("followStarCheckButton");
		reg_all = GTK_TOGGLE_BUTTON(lookup_widget("regallbutton"));
		reg_sel = GTK_TOGGLE_BUTTON(lookup_widget("regselbutton"));
		labelreginfo = GTK_LABEL(lookup_widget("labelregisterinfo"));
		notebook_reg = GTK_NOTEBOOK(lookup_widget("notebook_registration"));
		cumul_data = lookup_widget("check_button_comet");
	}

	if (!dont_change_reg_radio) {
		if (com.seq.selnum < com.seq.number)
			gtk_toggle_button_set_active(reg_sel, TRUE);
		else
			gtk_toggle_button_set_active(reg_all, TRUE);
	}

	selection_is_done = (com.selection.w > 0 && com.selection.h > 0);

	/* initialize default */
	gtk_notebook_set_current_page(notebook_reg, REG_PAGE_MISC);
	gtk_widget_set_visible(cumul_data, FALSE);

	/* getting the selected registration method */
	method = get_selected_registration_method();

	/* number of registered image */
	nb_images_reg = gtk_toggle_button_get_active(reg_all) ? com.seq.number : com.seq.selnum;

	if (method && nb_images_reg > 1 && (selection_is_done || method->sel == REQUIRES_NO_SELECTION)) {
		if (method->method_ptr == &register_star_alignment) {
			gtk_notebook_set_current_page(notebook_reg, REG_PAGE_GLOBAL);
		} else if (method->method_ptr == &register_comet) {
			gtk_notebook_set_current_page(notebook_reg, REG_PAGE_COMET);
		} else if (method->method_ptr == &register_3stars) {
			gtk_notebook_set_current_page(notebook_reg, REG_PAGE_3_STARS);
		}
		gtk_widget_set_visible(follow, method->method_ptr == &register_shift_fwhm);
		gtk_widget_set_visible(cumul_data, method->method_ptr == &register_comet);
		if (method->method_ptr == &register_3stars && com.seq.current != 0)
			gtk_label_set_text(labelreginfo, _("Make sure you load the first image"));
		else if (gfit.naxes[2] == 1 && gfit.bayer_pattern[0] != '\0')
			gtk_label_set_text(labelreginfo, _("Debayer the sequence for registration"));
		else gtk_label_set_text(labelreginfo, "");
		// the 3 stars method has special GUI requirements
		if (method->method_ptr != &register_3stars) {
			gtk_widget_set_sensitive(go_register, TRUE);
		}
	} else {
		gtk_widget_set_sensitive(go_register, FALSE);
		if (nb_images_reg <= 1 && !selection_is_done) {
			if (sequence_is_loaded()) {
				if (method && method->sel == REQUIRES_NO_SELECTION) {
					gtk_label_set_text(labelreginfo, _("Select images in the sequence"));
				} else {
					gtk_label_set_text(labelreginfo, _("Select an area in image first, and select images in the sequence"));
				}
			}
			else {
				gtk_label_set_text(labelreginfo, _("Load a sequence first"));
			}
		} else if (nb_images_reg <= 1) {
			gtk_label_set_text(labelreginfo, _("Select images in the sequence"));
		} else {
			gtk_label_set_text(labelreginfo, _("Select an area in image first"));
		}
	}
}

/* try to maximize the area within the image size (based on gfit)
 * hsteps and vsteps are used to resize the selection zone when it is larger than the image
 * they must be at least 2 */
void compute_fitting_selection(rectangle *area, int hsteps, int vsteps, int preserve_square) {
	//fprintf(stdout, "function entry: %d,%d,\t%dx%d\n", area->x, area->y, area->w, area->h);
	if (area->x >= 0 && area->x + area->w <= gfit.rx && area->y >= 0
			&& area->y + area->h <= gfit.ry)
		return;

	if (area->x < 0) {
		area->x++;
		if (area->x + area->w > gfit.rx) {
			/* reduce area */
			area->w -= hsteps;
			if (preserve_square) {
				area->h -= vsteps;
				area->y++;
			}
		}
	} else if (area->x + area->w > gfit.rx) {
		area->x--;
		if (area->x < 0) {
			/* reduce area */
			area->x++;
			area->w -= hsteps;
			if (preserve_square) {
				area->h -= vsteps;
				area->y++;
			}
		}
	}

	if (area->y < 0) {
		area->y++;
		if (area->y + area->h > gfit.ry) {
			/* reduce area */
			area->h -= hsteps;
			if (preserve_square) {
				area->w -= vsteps;
				area->x++;
			}
		}
	} else if (area->y + area->h > gfit.ry) {
		area->y--;
		if (area->y < 0) {
			/* reduce area */
			area->y++;
			area->h -= vsteps;
			if (preserve_square) {
				area->w -= hsteps;
				area->x++;
			}
		}
	}

	return compute_fitting_selection(area, hsteps, vsteps, preserve_square);
}

void get_the_registration_area(struct registration_args *reg_args,
		struct registration_method *method) {
	int max;
	switch (method->sel) {
		/* even in the case of REQUIRES_NO_SELECTION selection is needed for MatchSelection of starAlignment */
		case REQUIRES_NO_SELECTION:
		case REQUIRES_ANY_SELECTION:
			memcpy(&reg_args->selection, &com.selection, sizeof(rectangle));
			break;
		case REQUIRES_SQUARED_SELECTION:
			/* Passed arguments are X,Y of the center of the square and the size of
			 * the square. */
			if (com.selection.w > com.selection.h)
				max = com.selection.w;
			else
				max = com.selection.h;

			reg_args->selection.x = com.selection.x + com.selection.w / 2 - max / 2;
			reg_args->selection.w = max;
			reg_args->selection.y = com.selection.y + com.selection.h / 2 - max / 2;
			reg_args->selection.h = max;
			compute_fitting_selection(&reg_args->selection, 2, 2, 1);

			/* save it back to com.selection do display it properly */
			memcpy(&com.selection, &reg_args->selection, sizeof(rectangle));
			fprintf(stdout, "final area: %d,%d,\t%dx%d\n", reg_args->selection.x,
					reg_args->selection.y, reg_args->selection.w,
					reg_args->selection.h);
			redraw(com.cvport, REMAP_NONE);
			break;
	}
}

/* callback for 'Go register' button, GTK thread */
void on_seqregister_button_clicked(GtkButton *button, gpointer user_data) {
	struct registration_args *reg_args;
	struct registration_method *method;
	char *msg;
	GtkToggleButton *regall, *follow, *matchSel, *no_translate, *x2upscale,
			*cumul;
	GtkComboBox *cbbt_layers;
	GtkComboBoxText *ComboBoxRegInter, *ComboBoxTransfo;
	GtkSpinButton *minpairs;

	if (!reserve_thread()) {	// reentrant from here
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	if (!com.seq.regparam) {
		fprintf(stderr, "regparam should have been created before\n");
		// means that a call to seq_check_basic_data() or
		// check_or_allocate_regparam() is missing somewhere else
		unreserve_thread();
		return;
	}

	method = get_selected_registration_method();

	if (com.selection.w <= 0 && com.selection.h <= 0
			&& method->sel != REQUIRES_NO_SELECTION) {
		msg = siril_log_message(
				_("All prerequisites are not filled for registration. Select a rectangle first.\n"));
		siril_message_dialog( GTK_MESSAGE_WARNING, _("Warning"), msg);
		unreserve_thread();
		return;
	}

	reg_args = calloc(1, sizeof(struct registration_args));

	control_window_switch_to_tab(OUTPUT_LOGS);

	/* filling the arguments for registration */
	regall = GTK_TOGGLE_BUTTON(lookup_widget("regallbutton"));
	follow = GTK_TOGGLE_BUTTON(lookup_widget("followStarCheckButton"));
	matchSel = GTK_TOGGLE_BUTTON(lookup_widget("checkStarSelect"));
	no_translate = GTK_TOGGLE_BUTTON(lookup_widget("regTranslationOnly"));
	x2upscale = GTK_TOGGLE_BUTTON(lookup_widget("upscaleCheckButton"));
	cbbt_layers = GTK_COMBO_BOX(lookup_widget("comboboxreglayer"));
	ComboBoxRegInter = GTK_COMBO_BOX_TEXT(lookup_widget("ComboBoxRegInter"));
	cumul = GTK_TOGGLE_BUTTON(lookup_widget("check_button_comet"));
	minpairs = GTK_SPIN_BUTTON(lookup_widget("spinbut_minpairs"));
	ComboBoxTransfo = GTK_COMBO_BOX_TEXT(lookup_widget("comboreg_transfo"));

	reg_args->func = method->method_ptr;
	reg_args->seq = &com.seq;
	reg_args->reference_image = sequence_find_refimage(&com.seq);
	reg_args->process_all_frames = gtk_toggle_button_get_active(regall);
	reg_args->follow_star = gtk_toggle_button_get_active(follow);
	reg_args->matchSelection = gtk_toggle_button_get_active(matchSel);
	reg_args->translation_only = gtk_toggle_button_get_active(no_translate);
	reg_args->x2upscale = gtk_toggle_button_get_active(x2upscale);
	reg_args->cumul = gtk_toggle_button_get_active(cumul);
	reg_args->prefix = gtk_entry_get_text(GTK_ENTRY(lookup_widget("regseqname_entry")));
	reg_args->min_pairs = gtk_spin_button_get_value_as_int(minpairs);
	reg_args->type = gtk_combo_box_get_active(GTK_COMBO_BOX(ComboBoxTransfo));


	/* We check that available disk space is enough when:
	 * - activating the subpixel alignment, which requires generating a new
	 *   sequence with bigger images
	 * - using global star registration with rotation enabled, also generating a
	 *   new sequence */
	if (reg_args->x2upscale ||
			(method->method_ptr == register_star_alignment &&
			 !reg_args->translation_only)) {
		// first, remove the files that we are about to create
		remove_prefixed_sequence_files(reg_args->seq, reg_args->prefix);

		int nb_frames = reg_args->process_all_frames ? reg_args->seq->number : reg_args->seq->selnum;
		gint64 size = seq_compute_size(reg_args->seq, nb_frames, get_data_type(reg_args->seq->bitpix));
		if (reg_args->x2upscale)
			size *= 4;
		if (test_available_space(size) > 0) {
			free(reg_args);
			unreserve_thread();
			return;
		}
	} else if (method->method_ptr == register_comet) {
		pointf velocity = get_velocity();
		if ((velocity.x == 0.0 && velocity.y == 0.0)
				|| isinf(velocity.x) || isinf(velocity.y)) {
			msg = siril_log_color_message(_("The object is not moving, please check your registration data.\n"), "red");
			siril_message_dialog( GTK_MESSAGE_WARNING, _("Warning"), msg);
			free(reg_args);
			unreserve_thread();
			return;
		}
	}
	/* getting the selected registration layer from the combo box. The value is the index
	 * of the selected line, and they are in the same order than layers so there should be
	 * an exact matching between the two */
	reg_args->layer = gtk_combo_box_get_active(cbbt_layers);
	reg_args->interpolation = gtk_combo_box_get_active(GTK_COMBO_BOX(ComboBoxRegInter));
	get_the_registration_area(reg_args, method);	// sets selection
	reg_args->run_in_thread = TRUE;
	reg_args->load_new_sequence = FALSE; // only TRUE for global registration. Will be updated in this case

	clear_stars_list(); //to avoid problems with com.stars later on in the process

	msg = siril_log_color_message(_("Registration: processing using method: %s\n"),
			"green", method->name);
	msg[strlen(msg) - 1] = '\0';
	set_cursor_waiting(TRUE);
	set_progress_bar_data(msg, PROGRESS_RESET);

	start_in_reserved_thread(register_thread_func, reg_args);
}

// worker thread function for the registration
gpointer register_thread_func(gpointer p) {
	struct registration_args *args = (struct registration_args *) p;
	int retval;

	args->retval = args->func(args);

	if (args->seq->reference_image == -1) {
		// set new reference image: should we do it all the time?
		// also done in generated sequence in global.c
		args->seq->reference_image = sequence_find_refimage(args->seq);
	}
	writeseqfile(args->seq);
	retval = args->retval;
	if (!siril_add_idle(end_register_idle, args)) {
		free_sequence(args->seq, TRUE);
		free(args);
	}
	return GINT_TO_POINTER(retval);
}

// end of registration, GTK thread. Executed when started from the GUI and in
// the graphical command line but not from a script (headless mode)
static gboolean end_register_idle(gpointer p) {
	struct registration_args *args = (struct registration_args *) p;
	stop_processing_thread();

	if (!args->retval) {
		if (!args->load_new_sequence) {
			int chan = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("seqlist_dialog_combo")));
			update_seqlist();
			fill_sequence_list(args->seq, chan, FALSE);
			set_layers_for_registration();	// update display of available reg data
		}
		else {
			check_seq(0);
			update_sequences_list(args->new_seq_name);
		}
	}
	set_progress_bar_data(_("Registration complete."), PROGRESS_DONE);

	drawPlot();
	update_stack_interface(TRUE);
	adjust_sellabel();

	set_cursor_waiting(FALSE);

	free(args);
	return FALSE;
}
