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
#include <float.h>
#include <gsl/gsl_statistics.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/undo.h"
#include "core/processing.h"
#include "algos/statistics.h"
#include "algos/sorting.h"
#include "gui/progress_and_log.h"
#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "opencv/opencv.h"

#include "banding.h"

/*****************************************************************************
 *      B A N D I N G      R E D U C T I O N      M A N A G E M E N T        *
 ****************************************************************************/

int banding_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_) {
	struct banding_data *banding_args = (struct banding_data *)args->user;
	return BandingEngine(fit, banding_args->sigma, banding_args->amount,
			banding_args->protect_highlights, banding_args->applyRotation);
}

void apply_banding_to_sequence(struct banding_data *banding_args) {
	struct generic_seq_args *args = malloc(sizeof(struct generic_seq_args));
	args->seq = &com.seq;
	args->partial_image = FALSE;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = com.seq.selnum;
	args->prepare_hook = ser_prepare_hook;
	args->finalize_hook = ser_finalize_hook;
	args->save_hook = NULL;
	args->image_hook = banding_image_hook;
	args->idle_function = NULL;
	args->stop_on_error = FALSE;
	args->description = _("Banding Reduction");
	args->has_output = TRUE;
	args->new_seq_prefix = banding_args->seqEntry;
	args->load_new_sequence = TRUE;
	args->force_ser_output = FALSE;
	args->user = banding_args;
	args->already_in_a_thread = FALSE;
	args->parallel = TRUE;

	banding_args->fit = NULL;	// not used here

	start_in_new_thread(generic_sequence_worker, args);
}

// idle function executed at the end of the BandingEngine processing
gboolean end_BandingEngine(gpointer p) {
	struct banding_data *args = (struct banding_data *) p;
	stop_processing_thread();// can it be done here in case there is no thread?
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	update_used_memory();
	free(args);
	return FALSE;
}

static int fmul_layer(fits *a, int layer, float coeff) {
	WORD *buf;
	int i;

	if (coeff < 0.0)
		return 1;
	buf = a->pdata[layer];
	for (i = 0; i < a->rx * a->ry; ++i) {
		buf[i] = round_to_WORD(buf[i] * coeff);
	}
	invalidate_stats_from_fit(a);
	return 0;
}

/*** Reduces Banding in Canon DSLR images.
 * This code come from CanonBandingReduction.js v0.9.1, a script of
 * PixInsight, originally written by Georg Viehoever and
 * distributed under the terms of the GNU General Public License ******/
gpointer BandingEngineThreaded(gpointer p) {
	struct banding_data *args = (struct banding_data *) p;
	struct timeval t_start, t_end;

	siril_log_color_message(_("Banding Reducing: processing...\n"), "red");
	gettimeofday(&t_start, NULL);

	int retval = BandingEngine(args->fit, args->sigma, args->amount, args->protect_highlights, args->applyRotation);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	siril_add_idle(end_BandingEngine, args);

	return GINT_TO_POINTER(retval);
}

int BandingEngine(fits *fit, double sigma, double amount, gboolean protect_highlights, gboolean applyRotation) {
	int chan, row, i, ret = 0;
	WORD *line, *fixline;
	double minimum = DBL_MAX, globalsigma = 0.0;
	fits *fiximage = NULL;
	double invsigma = 1.0 / sigma;

	if (applyRotation) {
		point center = {gfit.rx / 2.0, gfit.ry / 2.0};
		cvRotateImage(fit, center, 90.0, -1, OPENCV_LINEAR);
	}

	if (new_fit_image(&fiximage, fit->rx, fit->ry, fit->naxes[2]))
		return 1;

	for (chan = 0; chan < fit->naxes[2]; chan++) {
		imstats *stat = statistics(NULL, -1, fit, chan, NULL, STATS_BASIC | STATS_MAD);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}
		double background = stat->median;
		double *rowvalue = calloc(fit->ry, sizeof(double));
		if (rowvalue == NULL) {
			PRINT_ALLOC_ERR;
			free_stats(stat);
			return 1;
		}
		if (protect_highlights) {
			globalsigma = stat->mad * MAD_NORM;
		}
		free_stats(stat);
		for (row = 0; row < fit->ry; row++) {
			line = fit->pdata[chan] + row * fit->rx;
			WORD *cpyline = calloc(fit->rx, sizeof(WORD));
			if (cpyline == NULL) {
				PRINT_ALLOC_ERR;
				free(rowvalue);
				return 1;
			}
			memcpy(cpyline, line, fit->rx * sizeof(WORD));
			int n = fit->rx;
			double median;
			if (protect_highlights) {
				quicksort_s(cpyline, n);
				WORD reject = round_to_WORD(
						background + invsigma * globalsigma);
				for (i = fit->rx - 1; i >= 0; i--) {
					if (cpyline[i] < reject)
						break;
					n--;
				}
				median = gsl_stats_ushort_median_from_sorted_data(cpyline, 1, n);
			} else {
				median = round_to_WORD(quickmedian(cpyline, n));
			}

			rowvalue[row] = background - median;
			minimum = min(minimum, rowvalue[row]);
			free(cpyline);
		}
		for (row = 0; row < fit->ry; row++) {
			fixline = fiximage->pdata[chan] + row * fiximage->rx;
			for (i = 0; i < fit->rx; i++)
				fixline[i] = round_to_WORD(rowvalue[row] - minimum);
		}
		free(rowvalue);
	}
	for (chan = 0; chan < fit->naxes[2]; chan++)
		fmul_layer(fiximage, chan, amount);
	ret = imoper(fit, fiximage, OPER_ADD);

	invalidate_stats_from_fit(fit);
	clearfits(fiximage);
	if ((!ret) && applyRotation) {
		point center = {gfit.rx / 2.0, gfit.ry / 2.0};
		cvRotateImage(fit, center, -90.0, -1, OPENCV_LINEAR);
	}

	return ret;
}

/***************** GUI for Canon Banding Reduction ********************/

void on_menuitem_fixbanding_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (sequence_is_loaded()) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkBandingSeq")), TRUE);
	}
	else if (single_image_is_loaded()) {
		// not a processing result
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkBandingSeq")), FALSE);
	}
	else
		return;
	siril_open_dialog("canon_fixbanding_dialog");
}

void on_button_ok_fixbanding_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("canon_fixbanding_dialog");
}

void on_button_apply_fixbanding_clicked(GtkButton *button, gpointer user_data) {
	static GtkRange *range_amount = NULL;
	static GtkRange *range_invsigma = NULL;
	static GtkToggleButton *toggle_protect_highlights_banding = NULL,
		*vertical = NULL, *seq = NULL;
	static GtkEntry *bandingSeqEntry = NULL;
	double amount, invsigma;
	gboolean protect_highlights;

	if (get_thread_run()) {
		siril_log_message(
				_(	"Another task is already in progress, ignoring new request.\n"));
		return;
	}

	struct banding_data *args = malloc(sizeof(struct banding_data));

	if (range_amount == NULL) {
		range_amount = GTK_RANGE(lookup_widget("scale_fixbanding_amount"));
		range_invsigma = GTK_RANGE(lookup_widget("scale_fixbanding_invsigma"));
		toggle_protect_highlights_banding = GTK_TOGGLE_BUTTON(
				lookup_widget("checkbutton_fixbanding"));
		vertical = GTK_TOGGLE_BUTTON(lookup_widget("checkBandingVertical"));
		seq = GTK_TOGGLE_BUTTON(lookup_widget("checkBandingSeq"));
		bandingSeqEntry = GTK_ENTRY(lookup_widget("entryBandingSeq"));
	}
	amount = gtk_range_get_value(range_amount);
	invsigma = gtk_range_get_value(range_invsigma);
	protect_highlights = gtk_toggle_button_get_active(
			toggle_protect_highlights_banding);

	if (!protect_highlights)
		undo_save_state(&gfit, "Processing: Canon Banding Reduction (amount=%.2lf)", amount);
	else
		undo_save_state(&gfit, "Processing: Canon Banding Reduction (amount=%.2lf, Protect=TRUE, invsigma=%.2lf)",
				amount, invsigma);

	args->fit = &gfit;
	args->protect_highlights = protect_highlights;
	args->amount = amount;
	args->sigma = invsigma;
	args->applyRotation = gtk_toggle_button_get_active(vertical);
	args->seqEntry = gtk_entry_get_text(bandingSeqEntry);
	set_cursor_waiting(TRUE);

	if (gtk_toggle_button_get_active(seq) && sequence_is_loaded()) {
		if (args->seqEntry && args->seqEntry[0] == '\0')
			args->seqEntry = "unband_";
		apply_banding_to_sequence(args);
	} else {
		start_in_new_thread(BandingEngineThreaded, args);
	}
}

void on_checkbutton_fixbanding_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	static GtkWidget *bandingHighlightBox = NULL;
	gboolean is_active;

	if (bandingHighlightBox == NULL) {
		bandingHighlightBox = lookup_widget("bandingHighlightBox");
	}

	is_active = gtk_toggle_button_get_active(togglebutton);
	gtk_widget_set_sensitive(bandingHighlightBox, is_active);
}
