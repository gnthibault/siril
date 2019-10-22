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

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <assert.h>
#include <fftw3.h>

#include "gui/callbacks.h"
#include "gui/dialogs.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "algos/statistics.h"

#include "fft.h"

enum {
	TYPE_CENTERED,
	TYPE_REGULAR
};

static void fft_to_spectra(fftw_complex *frequency_repr, double *as, double *ps,
		double *maxi, int nbdata) {
	unsigned int i;
	*maxi = 0.0;

	for (i = 0; i < nbdata; i++) {
		double r = creal(frequency_repr[i]);
		double im = cimag(frequency_repr[i]);
		as[i] = hypot(r, im);
		ps[i] = atan2(im, r);
		*maxi = max(as[i], *maxi);
	}
}

static void fft_to_freq(fftw_complex *frequency_repr, double *as, double *ps,
		int nbdata) {
	unsigned int i;

	for (i = 0; i < nbdata; i++) {
		frequency_repr[i] = as[i] * (cos(ps[i]) + I * sin(ps[i]));
	}
}

static void change_symmetry_forward(unsigned int width, unsigned int height,
		unsigned int i, unsigned int j, unsigned int *x, unsigned int *y) {

	*x = i + width / 2;
	if (*x >= width) {
		*x = *x - width;
	}
	*y = j + height / 2;
	if (*y >= height) {
		*y = *y - height;
	}
}

static void change_symmetry_backward(unsigned int width, unsigned int height,
		unsigned int i, unsigned int j, unsigned int *x, unsigned int *y) {

	*x = i + width - width / 2;
	if (*x >= width) {
		*x = *x - width;
	}
	*y = j + height - height / 2;
	if (*y >= height) {
		*y = *y - height;
	}
}

static void centered(WORD *buf, unsigned int width, unsigned int height,
		int type) {
	unsigned int i, j;

	WORD *temp = malloc(width * height * sizeof(WORD));
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			unsigned int x = i;
			unsigned int y = j;
			if (type == FFTW_FORWARD) {
				change_symmetry_forward(width, height, i, j, &x, &y);
			} else {
				change_symmetry_backward(width, height, i, j, &x, &y);
			}

			temp[j * width + i] = buf[y * width + x];
		}
	}

	memcpy(buf, temp, sizeof(WORD) * width * height);
	free(temp);
}

static void normalisation_spectra(unsigned int w, unsigned int h, double *modul, double *phase,
		WORD *abuf, WORD *pbuf, double maxi) {
	unsigned int i;

	for (i = 0; i < h * w; i++) {
		pbuf[i] = round_to_WORD(((phase[i] + M_PI) * USHRT_MAX_DOUBLE / (2 * M_PI)));
		abuf[i] = round_to_WORD((modul[i] * USHRT_MAX_DOUBLE / maxi));
	}
}

static void save_dft_information_in_gfit(fits *fit) {
	int i;

	strcpy(gfit.dft.ord, fit->dft.type);
	strcpy(gfit.dft.ord, fit->dft.ord);
	for (i = 0; i < fit->naxes[2]; i++)
		gfit.dft.norm[i] = fit->dft.norm[i];

}

static void FFTD(fits *fit, fits *x, fits *y, int type_order, int layer) {
	WORD *xbuf = x->pdata[layer];
	WORD *ybuf = y->pdata[layer];
	WORD *gbuf = fit->pdata[layer];
	unsigned int i;
	unsigned int width = fit->rx, height = fit->ry;
	double maxi;
	int nbdata = width * height;

	fftw_complex *spatial_repr = fftw_malloc(sizeof(fftw_complex) * nbdata);
	if (!spatial_repr) {
		PRINT_ALLOC_ERR;
		return;
	}
	fftw_complex *frequency_repr = fftw_malloc(sizeof(fftw_complex) * nbdata);
	if (!frequency_repr) {
		PRINT_ALLOC_ERR;
		fftw_free(spatial_repr);
		return;
	}

	/* copying image selection into the fftw data */
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static) if(nbdata > 15000)
#endif
	for (i = 0; i < nbdata; i++) {
		spatial_repr[i] = (double) gbuf[i];
	}

	/* we run the Fourier Transform */
	fftw_plan p = fftw_plan_dft_2d(height, width, spatial_repr, frequency_repr,
			FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);

	/* we compute modulus and phase */
	double *modul = malloc(nbdata * sizeof(double));
	double *phase = malloc(nbdata * sizeof(double));

	fft_to_spectra(frequency_repr, modul, phase, &maxi, nbdata);

	//We normalize the modulus and the phase
	normalisation_spectra(width, height, modul, phase, xbuf, ybuf, maxi);
	if (type_order == TYPE_CENTERED) {
		strcpy(x->dft.ord, "CENTERED");
		centered(xbuf, width, height, FFTW_FORWARD);
		centered(ybuf, width, height, FFTW_FORWARD);
	} else {
		strcpy(x->dft.ord, "REGULAR");
	}
	strcpy(y->dft.ord, x->dft.ord);
	x->dft.norm[layer] = maxi / USHRT_MAX_DOUBLE;

	free(modul);
	free(phase);
	fftw_destroy_plan(p);
	fftw_free(spatial_repr);
	fftw_free(frequency_repr);
}

static void FFTI(fits *fit, fits *xfit, fits *yfit, int type_order, int layer) {
	WORD *xbuf = xfit->pdata[layer];
	WORD *ybuf = yfit->pdata[layer];
	WORD *gbuf = fit->pdata[layer];
	unsigned int i;
	unsigned int width = xfit->rx;
	unsigned int height = xfit->ry;
	int nbdata = width * height;

	double *modul = calloc(1, nbdata * sizeof(double));
	double *phase = calloc(1, nbdata * sizeof(double));

	if (type_order == TYPE_CENTERED) {
		centered(xbuf, width, height, FFTW_BACKWARD);
		centered(ybuf, width, height, FFTW_BACKWARD);
	}

	for (i = 0; i < height * width; i++) {
		modul[i] = (double) xbuf[i] * (xfit->dft.norm[layer]);
		phase[i] = (double) ybuf[i] * (2 * M_PI / USHRT_MAX_DOUBLE);
		phase[i] -= M_PI;
	}

	fftw_complex* spatial_repr = fftw_malloc(sizeof(fftw_complex) * nbdata);
	if (!spatial_repr) {
		PRINT_ALLOC_ERR;
		free(modul);
		free(phase);
		return;
	}

	fftw_complex* frequency_repr = fftw_malloc(sizeof(fftw_complex) * nbdata);
	if (!frequency_repr) {
		PRINT_ALLOC_ERR;
		fftw_free(spatial_repr);
		free(modul);
		free(phase);
		return;
	}

	fft_to_freq(frequency_repr, modul, phase, nbdata);

	fftw_plan p = fftw_plan_dft_2d(height, width, frequency_repr, spatial_repr,
			FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);

	for (i = 0; i < nbdata; i++) {
		double pxl = creal(spatial_repr[i]) / nbdata;
		gbuf[i] = round_to_WORD(pxl);
	}
	delete_selected_area();
	invalidate_stats_from_fit(fit);

	free(modul);
	free(phase);
	fftw_destroy_plan(p);
	fftw_free(spatial_repr);
	fftw_free(frequency_repr);
}

// idle function executed at the end of the fourier_transform processing
static gboolean end_fourier_transform(gpointer p) {
	struct fft_data *args = (struct fft_data *)p;
	stop_processing_thread();
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	free(args->type);
	free(args);
	set_cursor_waiting(FALSE);
	update_used_memory();

	return FALSE;
}

gpointer fourier_transform(gpointer p) {
	struct fft_data *args = (struct fft_data *) p;
	unsigned int width = args->fit->rx;
	unsigned int height = args->fit->ry;
	int chan;
	struct timeval t_start, t_end;
	fits *tmp = NULL, *tmp1 = NULL, *tmp2 = NULL;

	siril_log_color_message(_("Fourier Transform: processing...\n"), "red");
	gettimeofday(&t_start, NULL);
	args->retval = 0;

	//type must be either "ffti" or "fftd"
	switch (args->type[3]) {
	default:
	case 'd':
	case 'D':
		if (new_fit_image(&tmp1, width, height, args->fit->naxes[2]) ||
				new_fit_image(&tmp2, width, height, args->fit->naxes[2])) {
			args->retval = 1;
			goto end;
		}

		for (chan = 0; chan < args->fit->naxes[2]; chan++)
			FFTD(args->fit, tmp1, tmp2, args->type_order, chan);
		strcpy(tmp1->dft.type, "SPECTRUM");
		if (savefits(args->modulus, tmp1)) {
			args->retval = 1;
			goto end;
		}
		strcpy(tmp2->dft.type, "PHASE");
		if (savefits(args->phase, tmp2)) {
			args->retval = 1;
			goto end;
		}

		/* We display the modulus on screen */
		if (copyfits(tmp1, &gfit, CP_ALLOC | CP_FORMAT | CP_COPYA, 0)) {
			args->retval = 1;
			goto end;
		}

		/* we copy the header informations */
		save_dft_information_in_gfit(tmp1);
		break;
	case 'i':
	case 'I':
		tmp = calloc(1, sizeof(fits));
		if (!tmp || readfits(args->modulus, tmp, NULL)) {
			PRINT_ALLOC_ERR;
			args->retval = 1;
			goto end;
		}
		tmp1 = calloc(1, sizeof(fits));
		if (!tmp1 || readfits(args->phase, tmp1, NULL)) {
			PRINT_ALLOC_ERR;
			args->retval = 1;
			goto end;
		}
		if (tmp->dft.ord[0] == 'C')		// CENTERED
			args->type_order = TYPE_CENTERED;
		else if (tmp->dft.ord[0] == 'R')	// REGULAR
			args->type_order = TYPE_REGULAR;
		else {
			args->retval = 1;
			siril_log_message(_("There is something wrong in your files\n"));
			goto end;
		}
		new_fit_image(&tmp2, width, height, tmp->naxes[2]);
		for (chan = 0; chan < args->fit->naxes[2]; chan++)
			FFTI(tmp2, tmp, tmp1, args->type_order, chan);
		/* We display the result on screen */
		if (copyfits(tmp2, &gfit, CP_ALLOC | CP_FORMAT | CP_COPYA, 0)) {
			args->retval = 1;
			goto end;
		}
	}

end:
	invalidate_stats_from_fit(args->fit);
	if (tmp)  { clearfits(tmp);  free(tmp);  }
	if (tmp1) { clearfits(tmp1); free(tmp1); }
	if (tmp2) { clearfits(tmp2); free(tmp2); }

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	siril_add_idle(end_fourier_transform, args);

	return GINT_TO_POINTER(args->retval);
}

/************************* GUI for FFT ********************************/

void on_button_fft_apply_clicked(GtkButton *button, gpointer user_data) {
	const char *mag, *phase;
	char *type = NULL, page;
	int type_order = -1;
	static GtkToggleButton *order = NULL;
	static GtkNotebook* notebookFFT = NULL;

	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return;
	}

	if (notebookFFT == NULL) {
		notebookFFT = GTK_NOTEBOOK(
				gtk_builder_get_object(builder, "notebook_fft"));
		order = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "fft_centered"));
	}

	page = gtk_notebook_get_current_page(notebookFFT);

	if (page == 0) {
		if (sequence_is_loaded()) {
			char *msg = siril_log_message(_("FFT does not work with sequences !\n"));
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), msg);
			set_cursor_waiting(FALSE);
			return;
		}
		if (!single_image_is_loaded()) {
			char *msg = siril_log_message(_("Open an image first !\n"));
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), msg);
			set_cursor_waiting(FALSE);
			return;
		}

		GtkEntry *entry_mag = GTK_ENTRY(lookup_widget("fftd_mag_entry"));
		GtkEntry *entry_phase = GTK_ENTRY(lookup_widget("fftd_phase_entry"));

		type_order = !gtk_toggle_button_get_active(order);
		type = strdup("fftd");
		mag = gtk_entry_get_text(entry_mag);
		phase = gtk_entry_get_text(entry_phase);
	} else {
		type = strdup("ffti");
		mag = gtk_file_chooser_get_filename(
				GTK_FILE_CHOOSER(lookup_widget("filechooser_mag")));
		phase = gtk_file_chooser_get_filename(
				GTK_FILE_CHOOSER(lookup_widget("filechooser_phase")));

		if (mag == NULL || phase == NULL) {
			char *msg = siril_log_message(_("Select magnitude and phase before !\n"));
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), msg);
			set_cursor_waiting(FALSE);
			free(type);
			return;
		}
		open_single_image(mag);
	}

	if ((mag != NULL) && (phase != NULL)) {
		set_cursor_waiting(TRUE);
		struct fft_data *args = malloc(sizeof(struct fft_data));
		args->fit = &gfit;
		args->type = type;
		args->modulus = mag;
		args->phase = phase;
		args->type_order = type_order;
		set_cursor_waiting(TRUE);
		start_in_new_thread(fourier_transform, args);
	} else {
		free(type);
	}
}

void on_button_fft_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("dialog_FFT");
}

void on_menuitem_fft_activate(GtkMenuItem *menuitem, gpointer user_data) {
	GtkFileChooserButton *magbutton, *phasebutton;

	magbutton = GTK_FILE_CHOOSER_BUTTON(lookup_widget("filechooser_mag"));
	phasebutton = GTK_FILE_CHOOSER_BUTTON(lookup_widget("filechooser_phase"));
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(magbutton), com.wd);
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(phasebutton), com.wd);
	siril_open_dialog("dialog_FFT");
}

