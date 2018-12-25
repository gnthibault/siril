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

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "algos/PSF.h"
#include "io/sequence.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "registration.h"


static point velocity = { 0.0, 0.0 };
static time_t t_of_image_1 = { 0 };
static time_t t_of_image_2 = { 0 };
static point pos_of_image1 = { 0 };
static point pos_of_image2 = { 0 };

static time_t FITS_date_key_to_sec(char *date) {
	struct tm timeinfo = { };
	time_t t;
	int year = 0, month = 0, day = 0, hour = 0, min = 0;
	float sec = 0.0;

	if (date[0] == '\0')
		return -1;

	sscanf(date, "%04d-%02d-%02dT%02d:%02d:%f", &year, &month, &day, &hour,
			&min, &sec);

	timeinfo.tm_year = year - 1900;
	timeinfo.tm_mon = month - 1;
	timeinfo.tm_mday = day;
	timeinfo.tm_hour = hour;
	timeinfo.tm_min = min;
	timeinfo.tm_sec = (int) sec;

	// Hopefully these are not needed
	timeinfo.tm_wday = 0;
	timeinfo.tm_yday = 0;
	timeinfo.tm_isdst = -1;

	/* get local time from timeinfo* */
	t = mktime(&timeinfo);

	return t;
}

static point compute_velocity(time_t t1, time_t t2, point d1, point d2) {
	double delta_t;
	point delta_d, px_per_hour;

	delta_t = (double) t2 - (double) t1;
	delta_d.x = d2.x - d1.x;
	delta_d.y = d2.y - d1.y;

	px_per_hour.x = delta_d.x / delta_t * 3600.0;
	px_per_hour.y = delta_d.y / delta_t * 3600.0;

	return px_per_hour;
}

static int get_comet_shift(fits *ref, fits *img, point px_per_hour, float *reg_x, float *reg_y) {
	time_t t_ref, t_img;
	double delta_t;

	t_ref = FITS_date_key_to_sec(ref->date_obs);
	t_img = FITS_date_key_to_sec(img->date_obs);

	delta_t = (double) t_img - (double) t_ref;
	delta_t /= 3600.0;
	*reg_x = delta_t * px_per_hour.x;
	*reg_y = delta_t * px_per_hour.y;

	return 0;
}

static void update_velocity() {
	GtkLabel *label = GTK_LABEL(lookup_widget("label1_comet"));

	velocity = compute_velocity(t_of_image_1, t_of_image_2, pos_of_image1, pos_of_image2);

	gchar *v_txt = g_strdup_printf("Δx: %.2lf, Δy: %.2lf", velocity.x, -velocity.y);
	gtk_label_set_text(label, v_txt);

	g_free(v_txt);
}

static void update_entry1(double x, double y) {
	GtkEntry *entry_x = GTK_ENTRY(lookup_widget("entry1_x_comet"));
	GtkEntry *entry_y = GTK_ENTRY(lookup_widget("entry1_y_comet"));
	gchar *txt_x, *txt_y;

	txt_x = g_strdup_printf("%7.2lf", x);
	txt_y = g_strdup_printf("%7.2lf", y);

	gtk_entry_set_text(entry_x, txt_x);
	gtk_entry_set_text(entry_y, txt_y);

	g_free(txt_x);
	g_free(txt_y);
}

static void update_entry2(double x, double y) {
	GtkEntry *entry_x = GTK_ENTRY(lookup_widget("entry2_x_comet"));
	GtkEntry *entry_y = GTK_ENTRY(lookup_widget("entry2_y_comet"));
	gchar *txt_x, *txt_y;

	txt_x = g_strdup_printf("%7.2lf", x);
	txt_y = g_strdup_printf("%7.2lf", y);

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
	fitted_PSF *result = NULL;
	int layer = get_reglayer();

	if (com.selection.h && com.selection.w) {
		set_cursor_waiting(TRUE);
		result = psf_get_minimisation(&gfit, layer, &com.selection, FALSE);
		if (result) {
			pos_of_image1.x = result->x0 + com.selection.x;
			pos_of_image1.y = com.selection.y + com.selection.h - result->y0;
			free(result);
			if (gfit.date_obs[0] == '\0') {
				siril_message_dialog(GTK_MESSAGE_ERROR,
						_("There is no timestamp stored in the file"),
						_("Siril cannot perform the registration without date information in the file."));
			} else {
				t_of_image_1 = FITS_date_key_to_sec(gfit.date_obs);
				update_entry1(pos_of_image1.x, pos_of_image1.y);
			}
		}
		set_cursor_waiting(FALSE);
	}
}

void on_button2_comet_clicked(GtkButton *button, gpointer p) {
	fitted_PSF *result = NULL;
	int layer = get_reglayer();

	if (com.selection.h && com.selection.w) {
		set_cursor_waiting(TRUE);
		result = psf_get_minimisation(&gfit, layer, &com.selection, FALSE);
		if (result) {
			pos_of_image2.x = result->x0 + com.selection.x;
			pos_of_image2.y = com.selection.y + com.selection.h - result->y0;
			free(result);
			if (gfit.date_obs[0] == '\0') {
				siril_message_dialog(GTK_MESSAGE_ERROR,
						_("There is no timestamp stored in the file"),
						_("Siril cannot perform the registration without date information in the file."));
			} else {
				t_of_image_2 = FITS_date_key_to_sec(gfit.date_obs);
				update_entry2(pos_of_image2.x, pos_of_image2.y);
			}
		}
		set_cursor_waiting(FALSE);
	}
}

int register_comet(struct registration_args *args) {
	int frame, ref_image, ret;
	int abort = 0;
	float cur_nb;
	regdata *current_regdata;
	fits ref = { 0 }, im = { 0 };

	if (args->seq->regparam[args->layer]) {
		current_regdata = args->seq->regparam[args->layer];
	} else {
		current_regdata = calloc(args->seq->number, sizeof(regdata));
		if (current_regdata == NULL) {
			fprintf(stderr, "Error allocating registration data\n");
			return -2;
		}
	}

	/* loading reference frame */
	ref_image = sequence_find_refimage(args->seq);

	ret = seq_read_frame(args->seq, ref_image, &ref);
	if (ret) {
		siril_log_message(_("Could not load reference image\n"));
		args->seq->regparam[args->layer] = NULL;
		free(current_regdata);
		return 1;
	}

	/* then we compare to other frames */
	if (args->process_all_frames)
		args->new_total = args->seq->number;
	else args->new_total = args->seq->selnum;

	cur_nb = 0.f;

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
			char tmpmsg[1024], tmpfilename[256];
			float reg_x, reg_y;

			if (!args->process_all_frames && !args->seq->imgparam[frame].incl)
				continue;

			if (!args->cumul) {
				/* data initialization */
				current_regdata[frame].shiftx = 0.0;
				current_regdata[frame].shifty = 0.0;
			}

			seq_get_image_filename(args->seq, frame, tmpfilename);
			g_snprintf(tmpmsg, 1024, _("Register: processing image %s"),
					tmpfilename);
			set_progress_bar_data(tmpmsg, PROGRESS_NONE);
			ret = seq_read_frame(args->seq, frame, &im);
			if (!ret) {
				get_comet_shift(&ref, &im, velocity, &reg_x, &reg_y);

				current_regdata[frame].shiftx += -reg_x;
				current_regdata[frame].shifty += reg_y;

			}
			cur_nb += 1.f;
			set_progress_bar_data(NULL, cur_nb / args->new_total);
			clearfits(&im);
		}
	}
	args->seq->regparam[args->layer] = current_regdata;
	if (args->x2upscale)
		args->seq->upscale_at_stacking = 2.0;
	else
		args->seq->upscale_at_stacking = 1.0;

	clearfits(&ref);
	update_used_memory();
	siril_log_message(_("Registration finished.\n"));

	return 0;
}

void on_entry_comet_changed(GtkEditable *editable, gpointer user_data) {
	GtkEntry *entry1_x = GTK_ENTRY(lookup_widget("entry1_x_comet"));
	GtkEntry *entry1_y = GTK_ENTRY(lookup_widget("entry1_y_comet"));
	GtkEntry *entry2_x = GTK_ENTRY(lookup_widget("entry2_x_comet"));
	GtkEntry *entry2_y = GTK_ENTRY(lookup_widget("entry2_y_comet"));

	pos_of_image1.x = atof(gtk_entry_get_text(entry1_x));
	pos_of_image1.y = atof(gtk_entry_get_text(entry1_y));
	pos_of_image2.x = atof(gtk_entry_get_text(entry2_x));
	pos_of_image2.y = atof(gtk_entry_get_text(entry2_y));

	update_velocity();
}
