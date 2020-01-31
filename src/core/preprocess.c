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

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/arithm.h"
#include "core/processing.h"
#include "core/OS_utils.h"
#include "algos/statistics.h"
#include "filters/cosmetic_correction.h"
#include "gui/callbacks.h"
#include "gui/histogram.h"
#include "gui/progress_and_log.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "io/conversion.h"
#include "io/ser.h"

#include "preprocess.h"

static double evaluateNoiseOfCalibratedImage(fits *fit, fits *dark,
		double k, gboolean allow_32bit_output) {
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

	soper(&dark_tmp, k, OPER_MUL, allow_32bit_output);
	ret = imoper(&fit_tmp, &dark_tmp, OPER_SUB, allow_32bit_output);
	if (ret) {
		clearfits(&dark_tmp);
		clearfits(&fit_tmp);
		return -1.0;
	}

	for (chan = 0; chan < fit->naxes[2]; chan++) {
		imstats *stat = statistics(NULL, -1, &fit_tmp, chan, &area, STATS_BASIC, FALSE);
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
		double tol, gboolean allow_32bits) {
	double c, d;
	double fc, fd;
	int iter = 0;

	c = b - GR * (b - a);
	d = a + GR * (b - a);
	// TODO: check return values of evaluateNoiseOfCalibratedImage
	fc = evaluateNoiseOfCalibratedImage(raw, dark, c, allow_32bits);
	fd = evaluateNoiseOfCalibratedImage(raw, dark, d, allow_32bits);
	do {
		siril_debug_print("Iter: %d (%1.2lf, %1.2lf)\n", ++iter, c, d);
		if (fc < 0.0 || fd < 0.0)
			return -1.0;
		if (fc < fd) {
			b = d;
			d = c;
			fd = fc;
			c = b - GR * (b - a);
			fc = evaluateNoiseOfCalibratedImage(raw, dark, c, allow_32bits);
		} else {
			a = c;
			c = d;
			fc = fd;
			d = a + GR * (b - a);
			fd = evaluateNoiseOfCalibratedImage(raw, dark, d, allow_32bits);

		}
	} while (fabs(c - d) > tol);
	return ((b + a) / 2.0);
}


static int preprocess(fits *raw, struct preprocessing_data *args) {
	int ret = 0;
	if (args->use_bias) {
		ret = imoper(raw, args->bias, OPER_SUB, args->allow_32bit_output);
	}

	/* if dark optimization, the master-dark has already been subtracted */
	if (!ret && args->use_dark && !args->use_dark_optim) {
		ret = imoper(raw, args->dark, OPER_SUB, args->allow_32bit_output);
	}

	if (!ret && args->use_flat) {
		// return value is an error if there is an overflow, but it is usual
		// for now, so don't treat as error
		/*ret =*/ siril_fdiv(raw, args->flat, args->normalisation, args->allow_32bit_output);
	}

	return ret;
}

static int darkOptimization(fits *raw, struct preprocessing_data *args) {
	double k0;
	double lo = 0.0, up = 2.0;
	int ret = 0;
	fits *dark = args->dark;
	fits dark_tmp = { 0 };

	if (memcmp(raw->naxes, dark->naxes, sizeof raw->naxes)) {
		siril_log_message(_("imoper: images must have same dimensions\n"));
		return 1;
	}

	if (copyfits(dark, &dark_tmp, CP_ALLOC | CP_COPYA | CP_FORMAT, 0))
		return 1;

	/* Minimization of background noise to find better k */
	invalidate_stats_from_fit(raw);
	k0 = goldenSectionSearch(raw, &dark_tmp, lo, up, 1E-3, args->allow_32bit_output);
	if (k0 < 0.0) {
		ret = -1;
	} else {
		siril_log_message(_("Dark optimization: k0=%.3lf\n"), k0);
		/* Multiply coefficient to master-dark */
		ret = soper(&dark_tmp, k0, OPER_MUL, args->allow_32bit_output);
		if (!ret)
			ret = imoper(raw, &dark_tmp, OPER_SUB, args->allow_32bit_output);
	}
	clearfits(&dark_tmp);
	return ret;
}

static int prepro_prepare_hook(struct generic_seq_args *args) {
	struct preprocessing_data *prepro = args->user;

	if (prepro->seq) {
		// checking disk space: removing old sequence and computing free space
		remove_prefixed_sequence_files(args->seq, prepro->ppprefix);

		int64_t size = seq_compute_size(args->seq, args->seq->number,
				prepro->use_flat ? DATA_FLOAT : DATA_USHORT);
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
			imstats *stat = statistics(NULL, -1, prepro->flat, RLAYER, NULL, STATS_BASIC, FALSE);
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
		if (darkOptimization(fit, prepro))
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

static void test_for_master_files(struct preprocessing_data *args) {
	GtkToggleButton *tbutton;
	GtkEntry *entry;

	tbutton = GTK_TOGGLE_BUTTON(lookup_widget("useoffset_button"));
	if (gtk_toggle_button_get_active(tbutton)) {
		const char *filename;
		entry = GTK_ENTRY(lookup_widget("offsetname_entry"));
		filename = gtk_entry_get_text(entry);
		if (filename[0] == '\0') {
			gtk_toggle_button_set_active(tbutton, FALSE);
		} else {
			const char *error = NULL;
			set_progress_bar_data(_("Opening offset image..."), PROGRESS_NONE);
			args->bias = calloc(1, sizeof(fits));
			if (!readfits(filename, args->bias, NULL)) {
				if (args->bias->naxes[2] != gfit.naxes[2]) {
					error = _("NOT USING OFFSET: number of channels is different");
				} else if (args->bias->naxes[0] != gfit.naxes[0] ||
						args->bias->naxes[1] != gfit.naxes[1]) {
					error = _("NOT USING OFFSET: image dimensions are different");
				} else {
					args->use_bias = TRUE;
				}

			} else error = _("NOT USING OFFSET: cannot open the file\n");
			if (error) {
				siril_log_message("%s\n", error);
				set_progress_bar_data(error, PROGRESS_DONE);
				free(args->bias);
				gtk_entry_set_text(entry, "");
				args->use_bias = FALSE;
			}
		}
	}

	tbutton = GTK_TOGGLE_BUTTON(lookup_widget("usedark_button"));
	if (gtk_toggle_button_get_active(tbutton)) {
		const char *filename;
		entry = GTK_ENTRY(lookup_widget("darkname_entry"));
		filename = gtk_entry_get_text(entry);
		if (filename[0] == '\0') {
			gtk_toggle_button_set_active(tbutton, FALSE);
		} else {
			const char *error = NULL;
			set_progress_bar_data(_("Opening dark image..."), PROGRESS_NONE);
			args->dark = calloc(1, sizeof(fits));
			if (!readfits(filename, args->dark, NULL)) {
				if (args->dark->naxes[2] != gfit.naxes[2]) {
					error = _("NOT USING DARK: number of channels is different");
				} else if (args->dark->naxes[0] != gfit.naxes[0] ||
						args->dark->naxes[1] != gfit.naxes[1]) {
					error = _("NOT USING DARK: image dimensions are different");
				} else {
					args->use_dark = TRUE;
				}

			} else error = _("NOT USING DARK: cannot open the file\n");
			if (error) {
				siril_log_message("%s\n", error);
				set_progress_bar_data(error, PROGRESS_DONE);
				free(args->dark);
				gtk_entry_set_text(entry, "");
				args->use_dark = FALSE;
			}
		}
		// dark optimization
		tbutton = GTK_TOGGLE_BUTTON(lookup_widget("checkDarkOptimize"));
		args->use_dark_optim = gtk_toggle_button_get_active(tbutton);

		// cosmetic correction
		tbutton = GTK_TOGGLE_BUTTON(lookup_widget("cosmEnabledCheck"));
		args->use_cosmetic_correction = gtk_toggle_button_get_active(tbutton);
	}

	tbutton = GTK_TOGGLE_BUTTON(lookup_widget("useflat_button"));
	if (gtk_toggle_button_get_active(tbutton)) {
		const char *filename;
		entry = GTK_ENTRY(lookup_widget("flatname_entry"));
		filename = gtk_entry_get_text(entry);
		if (filename[0] == '\0') {
			gtk_toggle_button_set_active(tbutton, FALSE);
		} else {
			const char *error = NULL;
			set_progress_bar_data(_("Opening flat image..."), PROGRESS_NONE);
			args->flat = calloc(1, sizeof(fits));
			if (!readfits(filename, args->flat, NULL)) {
				if (args->flat->naxes[2] != gfit.naxes[2]) {
					error = _("NOT USING FLAT: number of channels is different");
				} else if (args->flat->naxes[0] != gfit.naxes[0] ||
						args->flat->naxes[1] != gfit.naxes[1]) {
					error = _("NOT USING FLAT: image dimensions are different");
				} else {
					args->use_flat = TRUE;
				}

			} else error = _("NOT USING FLAT: cannot open the file\n");
			if (error) {
				siril_log_message("%s\n", error);
				set_progress_bar_data(error, PROGRESS_DONE);
				free(args->flat);
				gtk_entry_set_text(entry, "");
				args->use_flat = FALSE;
			}
		}
	}
}

void on_prepro_button_clicked(GtkButton *button, gpointer user_data) {
	struct preprocessing_data *args;
	GtkEntry *entry;
	GtkWidget *autobutton;
	GtkToggleButton *CFA, *debayer, *equalize_cfa, *compatibility, *stretch_cfa;
	GtkSpinButton *sigHot, *sigCold;

	if (!single_image_is_loaded() && !sequence_is_loaded())
		return;

	if (!single_image_is_loaded() && get_thread_run()) {
		siril_log_message(_("Another task is already in "
					"progress, ignoring new request.\n"));
		return;
	}

	args = calloc(1, sizeof(struct preprocessing_data));
	test_for_master_files(args);
	siril_log_color_message(_("Preprocessing...\n"), "red");
	gettimeofday(&args->t_start, NULL);
	args->autolevel = TRUE;
	args->normalisation = 1.0f;	// will be updated anyway

	// set output filename (preprocessed file name prefix)
	entry = GTK_ENTRY(lookup_widget("preproseqname_entry"));
	args->ppprefix = gtk_entry_get_text(entry);

	autobutton = lookup_widget("checkbutton_auto_evaluate");
	args->autolevel = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(autobutton));
	if (!args->autolevel) {
		GtkEntry *norm_entry = GTK_ENTRY(lookup_widget("entry_flat_norm"));
		args->normalisation = atof(gtk_entry_get_text(norm_entry));
	}

	CFA = GTK_TOGGLE_BUTTON(lookup_widget("cosmCFACheck"));
	debayer = GTK_TOGGLE_BUTTON(lookup_widget("checkButton_pp_dem"));
	compatibility = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_debayer_compatibility"));
	stretch_cfa = GTK_TOGGLE_BUTTON(lookup_widget("stretch_CFA_to16_button"));
	equalize_cfa = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_equalize_cfa"));
	sigHot = GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeHot"));
	sigCold = GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeCold"));

	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkSigCold"))))
		args->sigma[0] = gtk_spin_button_get_value(sigCold);
	else args->sigma[0] = -1.0;

	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkSigHot"))))
		args->sigma[1] = gtk_spin_button_get_value(sigHot);
	else args->sigma[1] = -1.0;

	args->is_cfa = gtk_toggle_button_get_active(CFA);
	args->compatibility = gtk_toggle_button_get_active(compatibility);
	args->debayer = gtk_toggle_button_get_active(debayer);
	args->stretch_cfa =  gtk_toggle_button_get_active(stretch_cfa);
	args->equalize_cfa = gtk_toggle_button_get_active(equalize_cfa);

	/****/

	if (sequence_is_loaded()) {
		args->is_sequence = TRUE;
		args->seq = &com.seq;
		set_cursor_waiting(TRUE);
		control_window_switch_to_tab(OUTPUT_LOGS);
		start_sequence_preprocessing(args, FALSE);
	} else {
		int retval;
		args->is_sequence = FALSE;
		set_cursor_waiting(TRUE);
		control_window_switch_to_tab(OUTPUT_LOGS);

		retval = preprocess_single_image(args);

		free(args);

		if (retval)
			set_progress_bar_data(_("Error in preprocessing."), PROGRESS_NONE);
		else {
			set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
			invalidate_gfit_histogram();
			open_single_image_from_gfit();
		}
		set_cursor_waiting(FALSE);
	}
}

void on_GtkButtonEvaluateCC_clicked(GtkButton *button, gpointer user_data) {
	GtkEntry *entry;
	GtkLabel *label[2];
	GtkWidget *widget[2];
	const char *filename;
	gchar *str[2];
	double sig[2];
	long icold = 0L, ihot = 0L;
	double rate, total;
	fits fit = { 0 };

	set_cursor_waiting(TRUE);
	sig[0] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeColdBox")));
	sig[1] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeHotBox")));
	widget[0] = lookup_widget("GtkLabelColdCC");
	widget[1] = lookup_widget("GtkLabelHotCC");
	label[0] = GTK_LABEL(lookup_widget("GtkLabelColdCC"));
	label[1] = GTK_LABEL(lookup_widget("GtkLabelHotCC"));
	entry = GTK_ENTRY(lookup_widget("darkname_entry"));
	filename = gtk_entry_get_text(entry);
	if (readfits(filename, &fit, NULL)) {
		str[0] = g_markup_printf_escaped(_("<span foreground=\"red\">ERROR</span>"));
		str[1] = g_markup_printf_escaped(_("<span foreground=\"red\">ERROR</span>"));
		gtk_label_set_markup(label[0], str[0]);
		gtk_label_set_markup(label[1], str[1]);
		set_cursor_waiting(FALSE);
		return;
	}
	count_deviant_pixels(&fit, sig, &icold, &ihot);
	total = fit.rx * fit.ry;
	clearfits(&fit);
	rate = (double)icold / total;
	/* 1% of cold pixels seems to be a reasonable limit */
	if (rate > 0.01) {
		str[0] = g_markup_printf_escaped(_("<span foreground=\"red\">Cold: %ld px</span>"), icold);
		gtk_widget_set_tooltip_text(widget[0], _("This value may be to high. Please, consider to change sigma value or uncheck the box."));
	}
	else {
		str[0] = g_markup_printf_escaped(_("Cold: %ld px"), icold);
		gtk_widget_set_tooltip_text(widget[0], "");
	}
	gtk_label_set_markup(label[0], str[0]);

	rate = (double)ihot / total;
	/* 1% of hot pixels seems to be a reasonable limit */
	if (rate > 0.01) {
		str[1] = g_markup_printf_escaped(_("<span foreground=\"red\">Hot: %ld px</span>"), ihot);
		gtk_widget_set_tooltip_text(widget[1], _("This value may be to high. Please, consider to change sigma value or uncheck the box."));
	}
	else {
		str[1] = g_markup_printf_escaped(_("Hot: %ld px"), ihot);
		gtk_widget_set_tooltip_text(widget[1], "");
	}
	gtk_label_set_markup(label[1], str[1]);
	g_free(str[0]);
	g_free(str[1]);
	set_cursor_waiting(FALSE);
}
