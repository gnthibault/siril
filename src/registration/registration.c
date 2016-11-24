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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>
#include <assert.h>
#include <float.h>
#include <gtk/gtk.h>
#ifdef MAC_INTEGRATION
#include "gtkmacintegration/gtkosxapplication.h"
#endif

#include "core/siril.h"
#include "gui/callbacks.h"
#include "gui/plot.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "registration/registration.h"
#include "registration/matching/misc.h"
#include "registration/matching/match.h"
#include "registration/matching/atpmatch.h"
#include "algos/star_finder.h"
#include "stacking/stacking.h"
#include "algos/PSF.h"
#include "gui/PSF_list.h"
#include "algos/quality.h"
#include "io/ser.h"
#ifdef HAVE_OPENCV
#include "opencv/opencv.h"
#include "opencv/ecc/ecc.h"
#endif

#define MAX_STARS_FITTED 2000
#undef DEBUG

static char *tooltip_text[] = { N_("One Star Registration: This is the simplest method to register deep-sky images. "
		"Because only one star is concerned for register, images are aligned using shifting "
		"(at a fraction of pixel). No rotation or scaling are performed. "
		"Shifts at pixel precision are saved in seq file."),
#ifdef HAVE_OPENCV
		N_("Global Star Alignment: This is a more powerful and accurate algorithm (but also slower) "
		"to perform deep-sky images. The global matching is based on triangle similarity method for automatically "
		"identify common stars in each image. "
		"A new sequence is created with the prefix of your choice (r_ by default)."),
#endif
		N_("Image Pattern Alignment: This is a simple registration by translation method "
		"using cross correlation in the spatial domain. This method is fast and is used to register "
		"planetary movies. It can also be used for some deep-sky images registration. "
		"Shifts at pixel precision are saved in seq file.")
#ifdef HAVE_OPENCV
,
		N_("Enhanced Correlation Coefficient Maximization: It is based on the enhanced correlation "
		"coefficient maximization algorithm. This method is more complex and slower than Image Pattern Alignment "
		"but no selection is required. It is good for moon surface images registration. Only translation is taken "
		"into account yet.")
#endif
};
/* callback for the selected area event */
void _reg_selected_area_callback() {
	update_reg_interface(TRUE);
}

static struct registration_method *reg_methods[NUMBER_OF_METHOD];
static gpointer register_thread_func(gpointer p);
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
#ifdef HAVE_OPENCV
	reg_methods[i++] = new_reg_method(_("Global Star Alignment (deep-sky)"),
			&register_star_alignment, REQUIRES_NO_SELECTION, REGTYPE_DEEPSKY);
#endif
	reg_methods[i++] = new_reg_method(_("Image Pattern Alignment (planetary - full disk)"),
			&register_shift_dft, REQUIRES_SQUARED_SELECTION, REGTYPE_PLANETARY);
#ifdef HAVE_OPENCV
	reg_methods[i++] = new_reg_method(_("Enhanced Correlation Coefficient (planetary - surfaces)"),
			&register_ecc, REQUIRES_NO_SELECTION, REGTYPE_PLANETARY);
#endif
	//if (theli_is_available())
	//	reg_methods[i++] = new_reg_method("theli", register_theli, REQUIRES_NO_SELECTION);
	reg_methods[i] = NULL;

	tip = g_string_new ("");
	for (j = 0; j < i; j ++) {
		g_string_append(tip, _(tooltip_text[j]));
		if (j < i - 1)
			g_string_append(tip, "\n\n");
	}
	ctip = g_string_free (tip, FALSE);
	gtk_widget_set_tooltip_text(lookup_widget("comboboxregmethod"), ctip);
	g_free(ctip);

	/* fill comboboxregmethod */
	regcombo = GTK_COMBO_BOX_TEXT(
			gtk_builder_get_object(builder, "comboboxregmethod"));
	gtk_combo_box_text_remove_all(regcombo);
	i = 0;
	while (reg_methods[i] != NULL) {
		gtk_combo_box_text_append_text(regcombo, reg_methods[i]->name);
		siril_log_message(_("Added a registration method: %s\n"),
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

	for (frame = 0; frame < args->seq->number; ++frame) {
		if (args->run_in_thread && !get_thread_run()) {
			break;
		}
		if (!args->process_all_frames && !args->seq->imgparam[frame].incl)
			continue;

		args->seq->regparam[args->layer][frame].quality -= q_min;
		args->seq->regparam[args->layer][frame].quality /= (q_max - q_min);
	}
}

/* register images: calculate shift in images to be aligned with the reference image;
 * images are not modified, only shift parameters are saved in regparam in the sequence.
 * layer is the layer on which the registration will be done, green by default (set in siril_init())
 */
int register_shift_dft(struct registration_args *args) {
	fits fit_ref, fit;
	int frame, size, sqsize;
	fftw_complex *ref, *in, *out, *convol;
	fftw_plan p, q;
	int ret, j;
	int plan;
	int abort = 0;
	float nb_frames, cur_nb;
	int ref_image;
	regdata *current_regdata;
	rectangle full_area;	// the area to use after getting image_part
	double q_max = 0, q_min = DBL_MAX;
	int q_index = -1;

	/* the selection needs to be squared for the DFT */
	assert(args->selection.w == args->selection.h);
	size = args->selection.w;
	full_area.x = full_area.y = 0;
	full_area.h = full_area.w = size;
	sqsize = size * size;

	if (args->process_all_frames)
		nb_frames = (float) args->seq->number;
	else
		nb_frames = (float) args->seq->selnum;

	if (!args->seq->regparam) {
		fprintf(stderr, "regparam should have been created before\n");
		return -1;
	}
	if (args->seq->regparam[args->layer]) {
		siril_log_message(
				_("Recomputing already existing registration for this layer\n"));
		current_regdata = args->seq->regparam[args->layer];
	} else {
		current_regdata = calloc(args->seq->number, sizeof(regdata));
		if (current_regdata == NULL) {
			printf("Error allocating registration data\n");
			return -2;
		}
	}

	/* loading reference frame */
	if (args->seq->reference_image == -1)
		ref_image = 0;
	else
		ref_image = args->seq->reference_image;

	set_progress_bar_data(
			_("Register DFT: loading and processing reference frame"),
			PROGRESS_NONE);
	memset(&fit_ref, 0, sizeof(fits));
	ret = seq_read_frame_part(args->seq, args->layer, ref_image, &fit_ref,
			&args->selection);

	if (ret) {
		siril_log_message(
				_("Register: could not load first image to register, aborting.\n"));
		free(current_regdata);
		clearfits(&fit_ref);
		return ret;
	}

	ref = fftw_malloc(sizeof(fftw_complex) * sqsize);
	in = fftw_malloc(sizeof(fftw_complex) * sqsize);
	out = fftw_malloc(sizeof(fftw_complex) * sqsize);
	convol = fftw_malloc(sizeof(fftw_complex) * sqsize);

	if (nb_frames > 200.f)
		plan = FFTW_MEASURE;
	else
		plan = FFTW_ESTIMATE;

	p = fftw_plan_dft_2d(size, size, ref, out, FFTW_FORWARD, plan);
	q = fftw_plan_dft_2d(size, size, convol, out, FFTW_BACKWARD, plan);

	// copying image selection into the fftw data
	for (j = 0; j < sqsize; j++)
		ref[j] = (double) fit_ref.data[j];

	// We don't need fit anymore, we can destroy it.
	current_regdata[ref_image].quality = QualityEstimate(&fit_ref, args->layer, QUALTYPE_NORMAL);
	clearfits(&fit_ref);
	fftw_execute_dft(p, ref, in); /* repeat as needed */
	current_regdata[ref_image].shiftx = 0;
	current_regdata[ref_image].shifty = 0;

	q_min = q_max = current_regdata[ref_image].quality;
	q_index = ref_image;

	cur_nb = 0.f;

	memset(&fit, 0, sizeof(fits));
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) firstprivate(fit) schedule(static) \
	if((args->seq->type == SEQ_REGULAR && fits_is_reentrant()) || args->seq->type == SEQ_SER)
#endif
	for (frame = 0; frame < args->seq->number; ++frame) {
		if (!abort) {
			if (args->run_in_thread && !get_thread_run()) {
				abort = 1;
				continue;
			}
			if (frame == ref_image)
				continue;
			if (!args->process_all_frames && !args->seq->imgparam[frame].incl)
				continue;

			char tmpmsg[1024], tmpfilename[256];

			seq_get_image_filename(args->seq, frame, tmpfilename);
			g_snprintf(tmpmsg, 1024, _("Register: processing image %s\n"),
					tmpfilename);
			set_progress_bar_data(tmpmsg, PROGRESS_NONE);
			if (!(seq_read_frame_part(args->seq, args->layer, frame, &fit,
					&args->selection))) {

				int x;
				fftw_complex *img = fftw_malloc(sizeof(fftw_complex) * sqsize);
				fftw_complex *out2 = fftw_malloc(sizeof(fftw_complex) * sqsize);

				// copying image selection into the fftw data
				for (x = 0; x < sqsize; x++)
					img[x] = (double) fit.data[x];

				// We don't need fit anymore, we can destroy it.
				current_regdata[frame].quality = QualityEstimate(&fit, args->layer,
						QUALTYPE_NORMAL);

				clearfits(&fit);

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

				fftw_execute_dft(p, img, out2); /* repeat as needed */

				fftw_complex *convol2 = fftw_malloc(sizeof(fftw_complex) * sqsize);

				for (x = 0; x < sqsize; x++) {
					convol2[x] = in[x] * conj(out2[x]);
				}

				fftw_execute_dft(q, convol2, out2); /* repeat as needed */
				fftw_free(convol2);

				int shift = 0;
				for (x = 1; x < sqsize; ++x) {
					if (creal(out2[x]) > creal(out2[shift])) {
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

				current_regdata[frame].shiftx = shiftx;
				current_regdata[frame].shifty = shifty;

				/* shiftx and shifty are the x and y values for translation that
				 * would make this image aligned with the reference image.
				 * WARNING: the y value is counted backwards, since the FITS is
				 * stored down from up.
				 */
#ifdef DEBUG
				fprintf(stderr,
						"reg: frame %d, shiftx=%d shifty=%d quality=%g\n",
						args->seq->imgparam[frame].filenum,
						current_regdata[frame].shiftx, current_regdata[frame].shifty,
						current_regdata[frame].quality);
#endif
#ifdef _OPENMP
#pragma omp atomic
#endif
				cur_nb += 1.f;
				set_progress_bar_data(NULL, cur_nb / nb_frames);
				fftw_free(img);
				fftw_free(out2);
			} else {
				//report_fits_error(ret, error_buffer);
				if (current_regdata == args->seq->regparam[args->layer])
					args->seq->regparam[args->layer] = NULL;
				free(current_regdata);
				abort = ret = 1;
				continue;
			}
		}
	}

	fftw_destroy_plan(p);
	fftw_destroy_plan(q);
	fftw_free(in);
	fftw_free(out);
	fftw_free(ref);
	fftw_free(convol);
	if (!ret) {
		args->seq->regparam[args->layer] = current_regdata;
		normalizeQualityData(args, q_min, q_max);
		update_used_memory();
		siril_log_message(_("Registration finished.\n"));
		siril_log_color_message(_("Best frame: #%d.\n"), "bold", q_index);
	}
	return ret;
}

/* register images: calculate shift in images to be aligned with the reference image;
 * images are not modified, only shift parameters are saved in regparam in the sequence.
 * layer is the layer on which the registration will be done, green by default (set in siril_init())
 */
int register_shift_fwhm(struct registration_args *args) {
	int frame, ref_image;
	float nb_frames, cur_nb;
	double reference_xpos, reference_ypos;
	double fwhm_min = DBL_MAX;
	int fwhm_index = -1;
	regdata *current_regdata;
	/* First and longest step: get the minimization data on one star for all
	 * images to register, which provides FWHM but also star coordinates */
	// TODO: detect that it was already computed, and don't do it again
	// -> should be done at a higher level and passed in the args
	if (do_fwhm_sequence_processing(args->seq, args->layer, TRUE, args->follow_star,
			args->run_in_thread, TRUE))	// stores in regparam
		return 1;

	current_regdata = args->seq->regparam[args->layer];

	if (args->process_all_frames)
		nb_frames = (float) args->seq->number;
	else
		nb_frames = (float) args->seq->selnum;

	/* loading reference frame */
	if (args->seq->reference_image == -1)
		ref_image = 0;
	else
		ref_image = args->seq->reference_image;
	if (!current_regdata[ref_image].fwhm_data) {
		siril_log_message(
				_("Registration PSF: failed to compute PSF for reference frame at least\n"));
		if (current_regdata != args->seq->regparam[args->layer])
			free(current_regdata);
		return -1;
	}
	reference_xpos = current_regdata[ref_image].fwhm_data->xpos;
	reference_ypos = current_regdata[ref_image].fwhm_data->ypos;

	fwhm_min = current_regdata[ref_image].fwhm_data->fwhmx;

	fwhm_index = ref_image;

	/* Second step: align image by aligning star coordinates together */
	for (frame = 0, cur_nb = 0.f; frame < args->seq->number; frame++) {
		double tmp;
		if (args->run_in_thread && !get_thread_run())
			break;
		if (!args->process_all_frames && !args->seq->imgparam[frame].incl)
			continue;
		if (frame == ref_image || !current_regdata[frame].fwhm_data) {
			current_regdata[frame].shiftx = 0;
			current_regdata[frame].shifty = 0;
			continue;
		}
		if (current_regdata[frame].fwhm < fwhm_min
				&& current_regdata[frame].fwhm > 0.0) {
			fwhm_min = current_regdata[frame].fwhm;
			fwhm_index = frame;
		}
		tmp = reference_xpos - current_regdata[frame].fwhm_data->xpos;
		current_regdata[frame].shiftx = round_to_int(tmp);
		tmp = current_regdata[frame].fwhm_data->ypos - reference_ypos;
		current_regdata[frame].shifty = round_to_int(tmp);

		/* shiftx and shifty are the x and y values for translation that
		 * would make this image aligned with the reference image.
		 * WARNING: the y value is counted backwards, since the FITS is
		 * stored down from up.
		 */
		fprintf(stderr, "reg: file %d, shiftx=%d shifty=%d\n",
				args->seq->imgparam[frame].filenum,
				current_regdata[frame].shiftx, current_regdata[frame].shifty);
		cur_nb += 1.f;
		set_progress_bar_data(NULL, cur_nb / nb_frames);
	}

	args->seq->regparam[args->layer] = current_regdata;
	update_used_memory();
	siril_log_message(_("Registration finished.\n"));
	siril_log_color_message(_("Best frame: #%d with fwhm=%.3g.\n"), "bold",
			fwhm_index, fwhm_min);
	return 0;
}

#ifdef HAVE_OPENCV

static void _print_result(TRANS *trans, float FWHMx, float FWHMy) {
	double rotation, scale;
	point shift;

	switch (trans->order) {
	case AT_TRANS_LINEAR:
		rotation = -atan2(trans->c, trans->b);
		shift.x = trans->a;
		shift.y = -trans->d;
		scale = sqrt(trans->b * trans->b + trans->c * trans->c);
		siril_log_color_message(_("Matching stars: done\n"), "green");
		siril_log_message(_("%d pair matches.\n"), trans->nr);
		siril_log_message(_("scale:%*.3f\n"), 13, scale);
		siril_log_message(_("rotation:%+*.2f deg\n"), 9, rotation * 180 / M_PI);
		siril_log_message(_("dx:%+*.2f px\n"), 15, shift.x);
		siril_log_message(_("dy:%+*.2f px\n"), 15, shift.y);
		siril_log_message(_("FWHMx:%*.2f px\n"), 12, FWHMx);
		siril_log_message(_("FWHMy:%*.2f px\n"), 12, FWHMy);
		break;
	default:
		siril_log_color_message(_("Not handled yet\n"), "red");
	}
}

int register_star_alignment(struct registration_args *args) {
	int frame, ref_image, ret, i;
	int abort = 0;
	int fitted_stars, failed = 0, skipped;
	float nb_frames, cur_nb;
	float FWHMx, FWHMy;
	fitted_PSF **stars;
	TRANS trans;
	regdata *current_regdata;
	starFinder sf;
	fits fit;
	struct ser_struct *new_ser = NULL;

	memset(&fit, 0, sizeof(fits));
	memset(&sf, 0, sizeof(starFinder));

	if (!args->seq->regparam) {
		fprintf(stderr, "regparam should have been created before\n");
		return -1;
	}
	if (args->seq->regparam[args->layer]) {
		siril_log_message(
				_("Recomputing already existing registration for this layer\n"));
		current_regdata = args->seq->regparam[args->layer];
		/* we reset all values as we built another sequence */
		memset(current_regdata, 0, args->seq->number * sizeof(regdata));
	} else {
		current_regdata = calloc(args->seq->number, sizeof(regdata));
		if (current_regdata == NULL) {
			printf("Error allocating registration data\n");
			return -2;
		}
	}

	if (args->process_all_frames)
		nb_frames = (float) args->seq->number;
	else
		nb_frames = (float) args->seq->selnum;

	/* loading reference frame */
	if (args->seq->reference_image == -1)
		ref_image = 0;
	else ref_image = args->seq->reference_image;

	/* first we're looking for stars in reference image */
	ret = seq_read_frame(args->seq, ref_image, &fit);
	if (ret) {
		siril_log_message(_("Could not load reference image\n"));
		if (current_regdata == args->seq->regparam[args->layer])
			args->seq->regparam[args->layer] = NULL;
		free(current_regdata);
		return 1;
	}
	siril_log_color_message(_("Reference Image:\n"), "green");

	if ((com.selection.w != 0) && (com.selection.h != 0) && args->matchSelection) {
		com.stars = peaker(&fit, args->layer, &sf, &com.selection);
	}
	else {
		com.stars = peaker(&fit, args->layer, &sf, NULL);
	}
	if (sf.nb_stars < AT_MATCH_MINPAIRS) {
		siril_log_message(
				_("There are not enough stars in reference image to perform alignment\n"));
		if (current_regdata == args->seq->regparam[args->layer])
			args->seq->regparam[args->layer] = NULL;
		free(current_regdata);
		return 1;
	}
	redraw(com.cvport, REMAP_NONE); // draw stars
#ifdef DEBUG
		FILE *pfile;

		pfile = fopen("ref.txt", "w+");
		fprintf(pfile, "REFERENCE IMAGE\n");
		for (i = 0; i < MAX_STARS_FITTED; i++) {
			fprintf(pfile, "%.3lf\t%.3lf\t%.3lf\n",
					com.stars[i]->xpos, com.stars[i]->ypos, com.stars[i]->mag);
		}
		fclose(pfile);
#endif
	fitted_stars = (sf.nb_stars > MAX_STARS_FITTED) ? MAX_STARS_FITTED : sf.nb_stars;
	FWHM_average(com.stars, &FWHMx, &FWHMy, fitted_stars);
	siril_log_message(_("FWHMx:%*.2f px\n"), 12, FWHMx);
	siril_log_message(_("FWHMy:%*.2f px\n"), 12, FWHMy);
	current_regdata[ref_image].fwhm = FWHMx;

	/* then we compare to other frames */
	if (args->process_all_frames)
		args->new_total = args->seq->number;
	else args->new_total = args->seq->selnum;
	args->imgparam = calloc(args->new_total, sizeof(imgdata));
	args->regparam = calloc(args->new_total, sizeof(regdata));

	if (args->seq->type == SEQ_SER) {
		char dest[256];

		new_ser = malloc(sizeof(struct ser_struct));

		const char *ptr = strrchr(args->seq->seqname, '/');
		if (ptr)
			snprintf(dest, 255, "%s%s.ser", args->prefix, ptr + 1);
		else
			snprintf(dest, 255, "%s%s.ser", args->prefix, args->seq->seqname);

		/* Here the last argument is NULL because we do not want copy SER file
		 * from the original. Indeed in the demosaicing case this would lead to
		 * a wrong file (B&W and not in RAW data). Moreover, header informations
		 * (like fps, local and UTC time, ...) have no sense now since some frames
		 * could be removed from the sequence.
		 */
		ser_create_file(dest, new_ser, TRUE, NULL);
	}

	skipped = 0;
	for (frame = 0, cur_nb = 0.f; frame < args->seq->number; frame++) {
		if (!abort) {
			if (args->run_in_thread && !get_thread_run()) {
				abort = 1;
				continue;
			}
			if (!args->process_all_frames && !args->seq->imgparam[frame].incl) {
				skipped++;
				continue;
			}

			ret = seq_read_frame(args->seq, frame, &fit);
			if (!ret) {
				char dest[256], filename[256];
				int nbpoints;

				if (frame != ref_image) {
					if (args->seq->type == SEQ_SER) {
						siril_log_color_message(_("Frame %d:\n"), "bold", frame);
					}
					stars = peaker(&fit, args->layer, &sf, NULL);
					if (sf.nb_stars < AT_MATCH_MINPAIRS) {
						siril_log_message(
								_("Not enough stars. Image %d skipped\n"), frame);
						args->new_total--;
						failed++;
						continue;
					}

#ifdef DEBUG
					FILE *pfile2;

					pfile2 = fopen("im.txt", "w+");
					fprintf(pfile2, "IMAGE %d\n", frame);
					for (i = 0; i < MAX_STARS_FITTED; i++) {
						fprintf(pfile2, "%.3lf\t%.3lf\t%.3lf\n", stars[i]->xpos,
								stars[i]->ypos, stars[i]->mag);
					}
					fprintf(pfile2, "\n\n", frame);
					fclose(pfile2);

#endif

					nbpoints = (sf.nb_stars < fitted_stars) ?
									sf.nb_stars : fitted_stars;

					if (star_match(stars, com.stars, nbpoints, &trans)) {
						siril_log_color_message(_("Cannot perform star matching. Image %d skipped\n"),
								"red", frame);
						args->new_total--;
						failed++;
						i = 0;
						while (i < MAX_STARS && stars[i])
							free(stars[i++]);
						free(stars);
						continue;
					}

					FWHM_average(stars, &FWHMx, &FWHMy, nbpoints);
					_print_result(&trans, FWHMx, FWHMy);
					current_regdata[frame].fwhm = FWHMx;

					if (!args->translation_only) {
						fits_flip_top_to_bottom(&fit);	// this is because in cvTransformImage, rotation center point is at (0, 0)

						/* An alternative method for generating an improved rotation: Basically we apply the operation to an image
						 * that is at least twice (or more) the size of the final image size wanted. After rotating the image,
						 * the image is resized down to its final size so as to produce a very sharp lines,
						 * edges, and much cleaner looking fonts. */
						cvResizeGaussian(&fit, fit.rx * SUPER_SAMPLING, fit.ry * SUPER_SAMPLING, OPENCV_CUBIC);
						trans.a *= SUPER_SAMPLING;
						trans.d *= SUPER_SAMPLING;
						cvTransformImage(&fit, trans, args->interpolation);
						cvResizeGaussian(&fit, fit.rx / SUPER_SAMPLING, fit.ry / SUPER_SAMPLING, OPENCV_CUBIC);

						fits_flip_top_to_bottom(&fit);
					}

					i = 0;
					while (i < MAX_STARS && stars[i])
						free(stars[i++]);
					free(stars);
				}

				if (!args->translation_only) {
					if (args->seq->type == SEQ_SER) {
						ser_write_frame_from_fit(new_ser, &fit,
								frame - failed - skipped);
						args->imgparam[frame - failed - skipped].filenum =
							frame - failed - skipped;
					} else {
						fit_sequence_get_image_filename(args->seq, frame, filename, TRUE);
						snprintf(dest, 256, "%s%s", args->prefix, filename);
						savefits(dest, &fit);
						args->imgparam[frame - failed - skipped].filenum = args->seq->imgparam[frame].filenum;
					}
					args->imgparam[frame - failed - skipped].incl = SEQUENCE_DEFAULT_INCLUDE;
					args->regparam[frame - failed - skipped].fwhm = current_regdata[frame].fwhm;	// not FWHMx because of the ref frame
				} else {
					current_regdata[frame].shiftx = trans.a;
					current_regdata[frame].shifty = -trans.d;
				}

				cur_nb += 1.f;
				set_progress_bar_data(NULL, cur_nb / nb_frames);
			}
		}
	}
	if (args->seq->type == SEQ_SER) {
		ser_write_and_close(new_ser);
		free(new_ser);
	}
	args->seq->regparam[args->layer] = current_regdata;
	update_used_memory();
	if (!abort) {
		siril_log_message(_("Registration finished.\n"));
		siril_log_color_message(_("%d images processed.\n"), "green",
				args->new_total + failed);
		siril_log_color_message(_("Total: %d failed, %d registered.\n"), "green",
				failed, args->new_total);

		args->load_new_sequence = !args->translation_only;

	}
	else {
		siril_log_message(_("Registration aborted.\n"));
		args->load_new_sequence = FALSE;
	}

	return args->new_total == 0;
}

int register_ecc(struct registration_args *args) {
	int frame, ref_image, ret, failed = 0;
	float nb_frames, cur_nb;
	regdata *current_regdata;
	fits ref, im;
	double q_max = 0, q_min = DBL_MAX;
	int q_index = -1;
	int abort = 0;

	if (!args->seq->regparam) {
		fprintf(stderr, "regparam should have been created before\n");
		return -1;
	}
	if (args->seq->regparam[args->layer]) {
		current_regdata = args->seq->regparam[args->layer];
		/* we reset all values as we built another sequence */
		memset(current_regdata, 0, args->seq->number * sizeof(regdata));
	} else {
		current_regdata = calloc(args->seq->number, sizeof(regdata));
		if (current_regdata == NULL) {
			printf("Error allocating registration data\n");
			return -2;
		}
	}

	if (args->process_all_frames)
		nb_frames = (float) args->seq->number;
	else
		nb_frames = (float) args->seq->selnum;

	/* loading reference frame */
	if (args->seq->reference_image == -1)
		ref_image = 0;
	else
		ref_image = args->seq->reference_image;

	memset(&ref, 0, sizeof(fits));

	/* first we're looking for stars in reference image */
	ret = seq_read_frame(args->seq, ref_image, &ref);
	if (ret) {
		siril_log_message(_("Could not load reference image\n"));
		if (current_regdata == args->seq->regparam[args->layer])
			args->seq->regparam[args->layer] = NULL;
		free(current_regdata);
		return 1;
	}
	current_regdata[ref_image].quality = QualityEstimate(&ref, args->layer, QUALTYPE_NORMAL);
	/* we make sure to free data in the destroyed fit */
	clearfits(&ref);
	/* Ugly code: as QualityEstimate destroys fit we need to reload it */
	seq_read_frame(args->seq, ref_image, &ref);
	q_min = q_max = current_regdata[ref_image].quality;
	q_index = ref_image;

	/* then we compare to other frames */
	if (args->process_all_frames)
		args->new_total = args->seq->number;
	else args->new_total = args->seq->selnum;

	cur_nb = 0.f;

	memset(&im, 0, sizeof(fits));
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) firstprivate(im) schedule(static) \
	if((args->seq->type == SEQ_REGULAR && fits_is_reentrant()) || args->seq->type == SEQ_SER)
#endif
	for (frame = 0; frame < args->seq->number; frame++) {
		if (!abort) {
			if (args->run_in_thread && !get_thread_run()) {
				abort = 1;
				continue;
			}
			if (!args->process_all_frames && !args->seq->imgparam[frame].incl)
				continue;
			current_regdata[frame].shiftx = 0;
			current_regdata[frame].shifty = 0;

			char tmpmsg[1024], tmpfilename[256];

			seq_get_image_filename(args->seq, frame, tmpfilename);
			g_snprintf(tmpmsg, 1024, _("Register: processing image %s\n"),
					tmpfilename);
			set_progress_bar_data(tmpmsg, PROGRESS_NONE);

			if (frame != ref_image) {

				ret = seq_read_frame(args->seq, frame, &im);
				if (!ret) {
					reg_ecc reg_param;
					memset(&reg_param, 0, sizeof(reg_ecc));

					if (findTransform(&ref, &im, args->layer, &reg_param)) {
						siril_log_message(
								_("Cannot perform ECC alignment for frame %d\n"),
								frame);
						/* We exclude this frame */
						com.seq.imgparam[frame].incl = FALSE;
#ifdef _OPENMP
#pragma omp atomic
#endif
						++failed;
						clearfits(&im);
						continue;
					}

					current_regdata[frame].quality = QualityEstimate(&im,
							args->layer, QUALTYPE_NORMAL);

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

					current_regdata[frame].shiftx = -round_to_int(reg_param.dx);
					current_regdata[frame].shifty = -round_to_int(reg_param.dy);

#ifdef _OPENMP
#pragma omp atomic
#endif
					cur_nb += 1.f;
					set_progress_bar_data(NULL, cur_nb / nb_frames);
					clearfits(&im);
				}
			}
		}
	}
	args->seq->regparam[args->layer] = current_regdata;
	normalizeQualityData(args, q_min, q_max);
	clearfits(&ref);
	update_used_memory();
	siril_log_message(_("Registration finished.\n"));
	if (failed)
		siril_log_color_message(_("%d frames were excluded.\n"), "red", failed);
	siril_log_color_message(_("Best frame: #%d.\n"), "bold", q_index);

	return 0;
}

#endif

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

int get_registration_layer() {
	int reglayer;
	GtkComboBox *registbox = GTK_COMBO_BOX(lookup_widget("comboboxreglayer"));

	if (!sequence_is_loaded())
		return -1;
	reglayer = gtk_combo_box_get_active(registbox);

	return reglayer;
}

/* Selects the "register all" or "register selected" according to the number of
 * selected images, if argument is false.
 * Verifies that enough images are selected and an area is selected.
 */
void update_reg_interface(gboolean dont_change_reg_radio) {
	static GtkWidget *go_register = NULL, *newSequence = NULL, *follow = NULL;
	static GtkLabel *labelreginfo = NULL;
	static GtkToggleButton *reg_all = NULL, *reg_sel = NULL;
	int nb_images_reg; /* the number of images to register */
	struct registration_method *method;

	if (!go_register) {
		go_register = lookup_widget("goregister_button");
		newSequence = lookup_widget("box29");
		follow = lookup_widget("followStarCheckButton");
		reg_all = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "regallbutton"));
		reg_sel = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "regselbutton"));
		labelreginfo = GTK_LABEL(
				gtk_builder_get_object(builder, "labelregisterinfo"));
	}

	if (!dont_change_reg_radio) {
		if (com.seq.selnum < com.seq.number)
			gtk_toggle_button_set_active(reg_sel, TRUE);
		else
			gtk_toggle_button_set_active(reg_all, TRUE);
	}

	/* getting the selected registration method */
	method = get_selected_registration_method();

	if (gtk_toggle_button_get_active(reg_all))
		nb_images_reg = com.seq.number;
	else
		nb_images_reg = com.seq.selnum;
	if (method && ((nb_images_reg > 1 && com.selection.w > 0 && com.selection.h > 0)
			|| (nb_images_reg > 1 && method->sel == REQUIRES_NO_SELECTION))) {
		gtk_widget_set_sensitive(go_register, TRUE);
		gtk_label_set_text(labelreginfo, "");
#ifdef HAVE_OPENCV
		gtk_widget_set_visible(newSequence, method->method_ptr == &register_star_alignment);
#endif
		gtk_widget_set_visible(follow, method->method_ptr == &register_shift_fwhm);
	} else {
		gtk_widget_set_sensitive(go_register, FALSE);
		gtk_widget_set_visible(newSequence, FALSE);
		if (nb_images_reg <= 1 && com.selection.w <= 0 && com.selection.h <= 0)
			if (!sequence_is_loaded())
				gtk_label_set_text(labelreginfo, _("Load a sequence first."));
			else
				gtk_label_set_text(labelreginfo,
					_("Select an area in image first, and select images in the sequence."));
		else if (nb_images_reg <= 1)
			gtk_label_set_text(labelreginfo, _("Select images in the sequence."));
		else
			gtk_label_set_text(labelreginfo, _("Select an area in image first."));
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
	case REQUIRES_NO_SELECTION:
		break;
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
	GtkToggleButton *regall, *follow, *matchSel, *no_translate;
	GtkComboBox *cbbt_layers;
	GtkComboBoxText *ComboBoxRegInter;

	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return;
	}

	/* getting the selected registration method */
	method = get_selected_registration_method();

	if (com.selection.w <= 0 && com.selection.h <= 0
			&& method->sel != REQUIRES_NO_SELECTION) {
		msg = siril_log_message(
				_("All prerequisites are not filled for registration. Select a rectangle first.\n"));
		show_dialog(msg, _("Warning"), "gtk-dialog-warning");
		return;
	}
	// TODO: check for reentrance

	reg_args = malloc(sizeof(struct registration_args));

	control_window_switch_to_tab(OUTPUT_LOGS);

	/* filling the arguments for registration */
	reg_args->seq = &com.seq;
	regall = GTK_TOGGLE_BUTTON(lookup_widget("regallbutton"));
	follow = GTK_TOGGLE_BUTTON(lookup_widget("followStarCheckButton"));
	matchSel = GTK_TOGGLE_BUTTON(lookup_widget("checkStarSelect"));
	no_translate = GTK_TOGGLE_BUTTON(lookup_widget("regTranslationOnly"));
	reg_args->process_all_frames = gtk_toggle_button_get_active(regall);
	reg_args->follow_star = gtk_toggle_button_get_active(follow);
	reg_args->matchSelection = gtk_toggle_button_get_active(matchSel);
	reg_args->translation_only = gtk_toggle_button_get_active(no_translate);
	/* getting the selected registration layer from the combo box. The value is the index
	 * of the selected line, and they are in the same order than layers so there should be
	 * an exact matching between the two */
	cbbt_layers = GTK_COMBO_BOX(
			gtk_builder_get_object(builder, "comboboxreglayer"));
	reg_args->layer = gtk_combo_box_get_active(cbbt_layers);
	ComboBoxRegInter = GTK_COMBO_BOX_TEXT(
			gtk_builder_get_object(builder, "ComboBoxRegInter"));
	reg_args->interpolation = gtk_combo_box_get_active(GTK_COMBO_BOX(ComboBoxRegInter));
	get_the_registration_area(reg_args, method);
	reg_args->func = method->method_ptr;
	reg_args->run_in_thread = TRUE;
	reg_args->prefix = gtk_entry_get_text(
			GTK_ENTRY(gtk_builder_get_object(builder, "regseqname_entry")));

	msg = siril_log_color_message(_("Registration: processing using method: %s\n"),
			"red", method->name);
	msg[strlen(msg) - 1] = '\0';
	gettimeofday(&(reg_args->t_start), NULL);
	set_cursor_waiting(TRUE);
	set_progress_bar_data(msg, PROGRESS_RESET);

	start_in_new_thread(register_thread_func, reg_args);
}

// worker thread
static gpointer register_thread_func(gpointer p) {
	struct registration_args *args = (struct registration_args *) p;
	args->retval = args->func(args);
	gdk_threads_add_idle(end_register_idle, args);
	return GINT_TO_POINTER(args->retval);	// not used anyway
}

// end of registration, GTK thread
static gboolean end_register_idle(gpointer p) {
	struct timeval t_end;
	struct registration_args *args = (struct registration_args *) p;
	stop_processing_thread();
	if (!args->retval) {
		writeseqfile(args->seq);
		fill_sequence_list(args->seq, com.cvport);
		set_layers_for_registration();	// update display of available reg data

		/* Load new sequence. Only star alignment method uses new sequence. */
#ifdef HAVE_OPENCV
		if (args->func == &register_star_alignment) {
			if (args->load_new_sequence) {
				sequence *seq;
				if (!(seq = malloc(sizeof(sequence)))) {
					fprintf(stderr, "could not allocate new sequence\n");
					goto failed_end;
				}
				initialize_sequence(seq, FALSE);

				/* we are not interested in the whole path */
				gchar *seqname = g_path_get_basename (com.seq.seqname);
				char *rseqname = malloc(
						strlen(args->prefix) + strlen(seqname) + 5);

				sprintf(rseqname, "%s%s.seq", args->prefix, seqname);
				g_free(seqname);
				unlink(rseqname);	// remove previous to overwrite
				//check_seq(0);		// search for the new sequence
				char *newname = remove_ext_from_filename(rseqname);
				seq->seqname = newname;
				seq->number = args->new_total;
				seq->selnum = args->new_total;
				seq->fixed = args->seq->fixed;
				seq->nb_layers = args->seq->nb_layers;
				seq->rx = args->seq->rx;
				seq->ry = args->seq->ry;
				seq->imgparam = args->imgparam;
				seq->regparam = calloc(seq->nb_layers, sizeof(regdata*));
				seq->regparam[args->layer] = args->regparam;
				seq->layers = calloc(seq->nb_layers, sizeof(layer_info));
				seq->beg = seq->imgparam[0].filenum;
				seq->end = seq->imgparam[seq->number-1].filenum;
				seq->type = args->seq->type;
				seq->ser_file = args->seq->ser_file;
				seq->current = -1;
				seq->needs_saving = TRUE;
				writeseqfile(seq);

				free_sequence(args->seq, FALSE);	// probably com.seq

				/* copy the new over com.seq to leave it in a usable state,
				 * even if it will be freed in update_sequences_list to be
				 * reloaded from the file */
				memcpy(&com.seq, seq, sizeof(sequence));

				update_sequences_list(rseqname);
				free(rseqname);
			}
			clear_stars_list();
		}
#endif
	}
	set_progress_bar_data(_("Registration complete."), PROGRESS_DONE);
	drawPlot();
#ifdef HAVE_OPENCV
failed_end:
#endif
	update_stack_interface();
	set_cursor_waiting(FALSE);
#ifdef MAC_INTEGRATION
	GtkosxApplication *osx_app = gtkosx_application_get();
	gtkosx_application_attention_request(osx_app, INFO_REQUEST);
	g_object_unref (osx_app);
#endif
	gettimeofday(&t_end, NULL);
	show_time(args->t_start, t_end);
	free(args);
	return FALSE;
}
