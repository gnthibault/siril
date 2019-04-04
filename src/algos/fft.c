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
#include "gui/progress_and_log.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "io/single_image.h"
#include "algos/statistics.h"
#include "algos/fft.h"

enum {
	TYPE_CENTERED,
	TYPE_REGULAR
};

static void fft_to_spectra(fits* fit, fftw_complex *frequency_repr, double *as,
		double *ps, int nbdata) {
	unsigned int i;

	for (i = 0; i < nbdata; i++) {
		double r = creal(frequency_repr[i]);
		double im = cimag(frequency_repr[i]);
		as[i] = hypot(r, im);
		ps[i] = atan2(im, r);
	}
}

static void fft_to_freq(fits* fit, fftw_complex *frequency_repr, double *as, double *ps, int nbdata) {
	unsigned int i;

	for (i = 0; i < nbdata; i++) {
		frequency_repr[i] = as[i] * (cos(ps[i]) + I * sin(ps[i]));
	}
}

static void change_symmetry(fits* fit, unsigned int i, unsigned int j, unsigned int *x,
		unsigned int *y) {
	int width = fit->rx;
	int height = fit->ry;
	if (i < width / 2 && j < height / 2) {
		*x = i + width / 2;
		*y = j + height / 2;
	}
	if (i >= width / 2 && j < height / 2) {
		*x = i - width / 2;
		*y = j + height / 2;
	}
	if (i < width / 2 && j >= height / 2) {
		*x = i + width / 2;
		*y = j - height / 2;
	}
	if (i >= width / 2 && j >= height / 2) {
		*x = i - width / 2;
		*y = j - height / 2;
	}
}

static void normalisation_spectra(fits* fit, double *modulus, double* phase,
		WORD *abuf, WORD *pbuf, int type_order) {
	unsigned int i, j;
	int width = fit->rx;
	int height = fit->ry;

	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			unsigned int x = i;
			unsigned int y = j;

			if (type_order == TYPE_CENTERED)
				change_symmetry(fit, i, j, &x, &y);
			pbuf[j * width + i] = round_to_WORD(((phase[y * width + x] + M_PI) * USHRT_MAX_DOUBLE
							/ (2 * M_PI)));
			abuf[j * width + i] = round_to_WORD((modulus[y * width + x] / width / height));
		}
	}
}

static void save_dft_information_in_gfit(fits *fit) {
	strcpy(gfit.dft.ord, fit->dft.type);
	strcpy(gfit.dft.ord, fit->dft.ord);
}

static void FFTD(fits *fit, fits *x, fits *y, int type_order, int layer) {
	WORD *xbuf = x->pdata[layer];
	WORD *ybuf = y->pdata[layer];
	WORD *gbuf = fit->pdata[layer];
	unsigned int i;
	int width = fit->rx, height = fit->ry;
	int nbdata = width * height;

	fftw_complex *spatial_repr = fftw_malloc(sizeof(fftw_complex) * nbdata);
	if (!spatial_repr) {
		return;
	}
	fftw_complex *frequency_repr = fftw_malloc(sizeof(fftw_complex) * nbdata);
	if (!frequency_repr) {
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
	double *modulus = malloc(nbdata * sizeof(double));
	double *phase = malloc(nbdata * sizeof(double));

	fft_to_spectra(fit, frequency_repr, modulus, phase, nbdata);

	//We normalize the modulus and the phase
	normalisation_spectra(fit, modulus, phase, xbuf, ybuf, type_order);
	if (type_order == TYPE_CENTERED)
		strcpy(x->dft.ord, "CENTERED");
	else
		strcpy(x->dft.ord, "REGULAR");
	strcpy(y->dft.ord, x->dft.ord);

	free(modulus);
	free(phase);
	fftw_destroy_plan(p);
	fftw_free(spatial_repr);
	fftw_free(frequency_repr);
}

static void FFTI(fits *fit, fits *xfit, fits *yfit, int type_order, int layer) {
	WORD *xbuf = xfit->pdata[layer];
	WORD *ybuf = yfit->pdata[layer];
	WORD *gbuf = fit->pdata[layer];
	unsigned int i, j;
	int width = xfit->rx;
	int height = xfit->ry;
	int nbdata = width * height;

	double *modulus = calloc(1, nbdata * sizeof(double));
	double *phase = calloc(1, nbdata * sizeof(double));

	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			unsigned int x = i;
			unsigned int y = j;

			if (type_order == TYPE_CENTERED)
				change_symmetry(fit, i, j, &x, &y);
			modulus[j * width + i] = (double) xbuf[y * width + x] * (width * height);
			phase[j * width + i] = (double) ybuf[y * width + x]	* (2 * M_PI / USHRT_MAX_DOUBLE);
			phase[j * width + i] -= M_PI;
		}
	}

	fftw_complex* spatial_repr = fftw_malloc(sizeof(fftw_complex) * nbdata);
	if (!spatial_repr) {
		return;
	}

	fftw_complex* frequency_repr = fftw_malloc(sizeof(fftw_complex) * nbdata);
	if (!frequency_repr) {
		fftw_free(spatial_repr);
		return;
	}

	fft_to_freq(fit, frequency_repr, modulus, phase, nbdata);

	fftw_plan p = fftw_plan_dft_2d(height, width, frequency_repr, spatial_repr,
			FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);

	for (i = 0; i < nbdata; i++) {
		double pxl = creal(spatial_repr[i]) / nbdata;
		gbuf[i] = round_to_WORD(pxl);
	}
	delete_selected_area();
	invalidate_stats_from_fit(fit);

	free(modulus);
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
	int width = args->fit->rx;
	int height = args->fit->ry;
	int ndata = width * height, chan;
	struct timeval t_start, t_end;
	WORD *from[3], *to[3];
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
			args->retval = 1;
			goto end;
		}
		tmp1 = calloc(1, sizeof(fits));
		if (!tmp1 || readfits(args->phase, tmp1, NULL)) {
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
		new_fit_image(&args->fit, width, height, tmp->naxes[2]);
		for (chan = 0; chan < args->fit->naxes[2]; chan++) {
			from[chan] = tmp2->pdata[chan];
			to[chan] = args->fit->pdata[chan];
			memcpy(to[chan], from[chan], ndata * sizeof(WORD));
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
