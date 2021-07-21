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
#include "core/siril_date.h"
#include "core/processing.h"
#include "algos/PSF.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "registration.h"

static pointf velocity = { 0.f, 0.f };
static GDateTime *t_of_image_1 = NULL;
static GDateTime *t_of_image_2 = NULL;
static pointf pos_of_image1 = { 0 };
static pointf pos_of_image2 = { 0 };

static pointf compute_velocity(GDateTime *t1, GDateTime *t2, pointf d1, pointf d2) {
	float delta_t;
	pointf delta_d, px_per_hour = { 0.f, 0.f };

	if (t1 && t2) {
		delta_t = g_date_time_difference(t2, t1);
		delta_d.x = d2.x - d1.x;
		delta_d.y = d2.y - d1.y;

		px_per_hour.x = delta_d.x / delta_t * 3600000000.f;
		px_per_hour.y = delta_d.y / delta_t * 3600000000.f;
	}

	return px_per_hour;
}

static int get_comet_shift(GDateTime *ref, GDateTime *img, pointf px_per_hour, pointf *reg) {
	float delta_t;

	if (img && ref) {
		delta_t = (float) g_date_time_difference(img, ref);
		delta_t /= 3600000000.f;
		reg->x = delta_t * px_per_hour.x;
		reg->y = delta_t * px_per_hour.y;
	}
	return 0;
}

static void update_velocity() {
	GtkLabel *label = GTK_LABEL(lookup_widget("label1_comet"));

	velocity = compute_velocity(t_of_image_1, t_of_image_2, pos_of_image1, pos_of_image2);

	gchar *v_txt = g_strdup_printf("Δx: %.2lf, Δy: %.2lf", velocity.x, -velocity.y);
	gtk_label_set_text(label, v_txt);

	g_free(v_txt);
}

static void update_entry1(pointf pt) {
	GtkEntry *entry_x = GTK_ENTRY(lookup_widget("entry1_x_comet"));
	GtkEntry *entry_y = GTK_ENTRY(lookup_widget("entry1_y_comet"));
	gchar *txt_x, *txt_y;

	txt_x = g_strdup_printf("%7.2f", pt.x);
	txt_y = g_strdup_printf("%7.2f", pt.y);

	gtk_entry_set_text(entry_x, txt_x);
	gtk_entry_set_text(entry_y, txt_y);

	g_free(txt_x);
	g_free(txt_y);
}

static void update_entry2(pointf pt) {
	GtkEntry *entry_x = GTK_ENTRY(lookup_widget("entry2_x_comet"));
	GtkEntry *entry_y = GTK_ENTRY(lookup_widget("entry2_y_comet"));
	gchar *txt_x, *txt_y;

	txt_x = g_strdup_printf("%7.2f", pt.x);
	txt_y = g_strdup_printf("%7.2f", pt.y);

	gtk_entry_set_text(entry_x, txt_x);
	gtk_entry_set_text(entry_y, txt_y);

	g_free(txt_x);
	g_free(txt_y);
}

static int get_reglayer() {
	GtkComboBox *cbbt_layers = GTK_COMBO_BOX(lookup_widget("comboboxreglayer"));

	return gtk_combo_box_get_active(cbbt_layers);
}

void on_button1_comet_clicked(GtkButton *button, gpointer p) {
	psf_star *result = NULL;
	int layer = get_reglayer();

	if (com.selection.h && com.selection.w) {
		set_cursor_waiting(TRUE);
		result = psf_get_minimisation(&gfit, layer, &com.selection, FALSE, FALSE, TRUE);
		if (result) {
			pos_of_image1.x = result->x0 + com.selection.x;
			pos_of_image1.y = com.selection.y + com.selection.h - result->y0;
			free_psf(result);
			if (!gfit.date_obs) {
				siril_message_dialog(GTK_MESSAGE_ERROR,
						_("There is no timestamp stored in the file"),
						_("Siril cannot perform the registration without date information in the file."));
			} else {
				if (t_of_image_1) {
					g_date_time_unref(t_of_image_1);
				}
				t_of_image_1 = g_date_time_ref(gfit.date_obs);
				if (!t_of_image_1) {
					siril_message_dialog(GTK_MESSAGE_ERROR,
							_("Unable to convert DATE-OBS to a valid date"),
							_("Siril cannot convert the DATE-OBS keyword into a valid date needed in the alignment."));
				}
				update_entry1(pos_of_image1);
			}
		}
		set_cursor_waiting(FALSE);
	}
}

void on_button2_comet_clicked(GtkButton *button, gpointer p) {
	psf_star *result = NULL;
	int layer = get_reglayer();

	if (com.selection.h && com.selection.w) {
		set_cursor_waiting(TRUE);
		result = psf_get_minimisation(&gfit, layer, &com.selection, FALSE, FALSE, TRUE);
		if (result) {
			pos_of_image2.x = result->x0 + com.selection.x;
			pos_of_image2.y = com.selection.y + com.selection.h - result->y0;
			free_psf(result);
			if (!gfit.date_obs) {
				siril_message_dialog(GTK_MESSAGE_ERROR,
						_("There is no timestamp stored in the file"),
						_("Siril cannot perform the registration without date information in the file."));
			} else {
				if (t_of_image_2) {
					g_date_time_unref(t_of_image_2);
				}
				t_of_image_2 = g_date_time_ref(gfit.date_obs);
				if (!t_of_image_2) {
					siril_message_dialog(GTK_MESSAGE_ERROR,
							_("Unable to convert DATE-OBS to a valid date"),
							_("Siril cannot convert the DATE-OBS keyword into a valid date needed in the alignment."));
				}
				update_entry2(pos_of_image2);
			}
		}
		set_cursor_waiting(FALSE);
	}
}

void on_entry_comet_changed(GtkEditable *editable, gpointer user_data) {
	GtkEntry *entry1_x = GTK_ENTRY(lookup_widget("entry1_x_comet"));
	GtkEntry *entry1_y = GTK_ENTRY(lookup_widget("entry1_y_comet"));
	GtkEntry *entry2_x = GTK_ENTRY(lookup_widget("entry2_x_comet"));
	GtkEntry *entry2_y = GTK_ENTRY(lookup_widget("entry2_y_comet"));

	pos_of_image1.x = g_ascii_strtod(gtk_entry_get_text(entry1_x), NULL);
	pos_of_image1.y = g_ascii_strtod(gtk_entry_get_text(entry1_y), NULL);
	pos_of_image2.x = g_ascii_strtod(gtk_entry_get_text(entry2_x), NULL);
	pos_of_image2.y = g_ascii_strtod(gtk_entry_get_text(entry2_y), NULL);

	update_velocity();
}

pointf get_velocity() {
	return velocity;
}

/***** generic moving object registration *****/

struct comet_align_data {
	struct registration_args *regargs;
	regdata *current_regdata;
	GDateTime *reference_date;
};

static int comet_align_prepare_hook(struct generic_seq_args *args) {
	struct comet_align_data *cadata = args->user;
	struct registration_args *regargs = cadata->regargs;
	int ref_image;
	fits ref = { 0 };

	if (args->seq->regparam[regargs->layer]) {
		cadata->current_regdata = args->seq->regparam[regargs->layer];
	} else {
		cadata->current_regdata = calloc(args->seq->number, sizeof(regdata));
		if (cadata->current_regdata == NULL) {
			PRINT_ALLOC_ERR;
			return -2;
		}
		args->seq->regparam[regargs->layer] = cadata->current_regdata;
	}

	/* loading reference frame */
	ref_image = sequence_find_refimage(args->seq);

	// TODO: reading only the date can be made in a less resource-consuming way
	if (seq_read_frame(args->seq, ref_image, &ref, FALSE, -1)) {
		siril_log_message(_("Could not load reference image\n"));
		args->seq->regparam[regargs->layer] = NULL;
		free(cadata->current_regdata);
		return 1;
	}
	cadata->reference_date = g_date_time_ref(ref.date_obs);
	clearfits(&ref);

	if (regargs->x2upscale)
		args->seq->upscale_at_stacking = 2.0;
	else args->seq->upscale_at_stacking = 1.0;
	return 0;
}

static int comet_align_image_hook(struct generic_seq_args *args, int out_index, int in_index, fits *fit, rectangle *_) {
	struct comet_align_data *cadata = args->user;
	struct registration_args *regargs = cadata->regargs;
	pointf reg = { 0.f, 0.f };

	if (!regargs->cumul) {
		/* data initialization */
		set_shifts(args->seq, in_index, regargs->layer, reg.x, reg.y, FALSE);
	}

	get_comet_shift(cadata->reference_date, fit->date_obs, velocity, &reg);

	/* get_comet_shift does not car about orientation of image */
	set_shifts(args->seq, in_index, regargs->layer, -reg.x, reg.y, FALSE);
	return 0;
}

static int comet_align_finalize_hook(struct generic_seq_args *args) {
	struct comet_align_data *cadata = args->user;
	struct registration_args *regargs = cadata->regargs;

	if (args->retval) {
		free(args->seq->regparam[regargs->layer]);
		args->seq->regparam[regargs->layer] = NULL;
	}

	if (cadata->reference_date)
		g_date_time_unref(cadata->reference_date);

	free(cadata);
	cadata = NULL;
	return 0;
}

int register_comet(struct registration_args *regargs) {
	struct generic_seq_args *args = create_default_seqargs(regargs->seq);
	/* we don't need to read image data, for simplicity we just read one
	 * pixel from it, making sure the header is read */
	args->partial_image = TRUE;
	args->area.x = 0; args->area.y = 0;
	args->area.w = 1; args->area.h = 1;
	args->layer_for_partial = 0;
	args->get_photometry_data_for_partial = TRUE;

	if (!regargs->process_all_frames) {
		args->filtering_criterion = seq_filter_included;
		args->nb_filtered_images = regargs->seq->selnum;
	}
	args->prepare_hook = comet_align_prepare_hook;
	args->image_hook = comet_align_image_hook;
	args->finalize_hook = comet_align_finalize_hook;
	args->description = _("Moving object registration");
	args->already_in_a_thread = TRUE;

	struct comet_align_data *cadata = calloc(1, sizeof(struct comet_align_data));
	if (!cadata) {
		free(args);
		return -1;
	}
	cadata->regargs = regargs;
	args->user = cadata;

	generic_sequence_worker(args);
	return args->retval;
}

