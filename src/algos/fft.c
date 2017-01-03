/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2017 team free-astro (see more in AUTHORS file)
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
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "io/single_image.h"
#include "algos/fft.h"

void fft_to_spectra(fits* fit, fftw_complex *frequency_repr, double *as,
		double *ps) {
	unsigned int i;
	int nbdata = fit->rx * fit->ry;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static) if(nbdata > 3000)
#endif
	for (i = 0; i < nbdata; i++) {
		as[i] = hypot(creal(frequency_repr[i]), cimag(frequency_repr[i]));
		ps[i] = atan2(cimag(frequency_repr[i]), creal(frequency_repr[i]));
	}
}

void fft_to_fr(fits* fit, fftw_complex *frequency_repr, double *as, double *ps) {
	unsigned int i;
	int nbdata = fit->rx * fit->ry;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static) if(nbdata > 3000)
#endif
	for (i = 0; i < nbdata; i++) {
		frequency_repr[i] = as[i] * (cos(ps[i]) + I * sin(ps[i]));
	}
}

void change_symmetry(fits* fit, unsigned int i, unsigned int j, unsigned int *x,
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

double normalisation_spectra(fits* fit, double *modulus, double* phase,
		WORD *abuf, WORD *pbuf, int type_order) {
	unsigned int i, j;
	int width = fit->rx;
	int height = fit->ry;
	double max_m = 0.0;

	for (i = 0; i < width * height; i++)
		max_m = max(max_m, modulus[i]);

	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			unsigned int x = i;
			unsigned int y = j;

			if (type_order == 0)
				change_symmetry(fit, i, j, &x, &y);
			pbuf[j * width + i] = round_to_WORD(
					((phase[y * width + x] + M_PI) * USHRT_MAX_DOUBLE
							/ (2 * M_PI)));
			abuf[j * width + i] = round_to_WORD(
					(modulus[y * width + x] * USHRT_MAX_DOUBLE / max_m));
		}
	}
	return max_m / USHRT_MAX_DOUBLE;
}

void save_dft_information_in_gfit(fits *fit) {
	int i;

	strcpy(gfit.dft_ord, fit->dft_type);
	strcpy(gfit.dft_ord, fit->dft_ord);
	for (i = 0; i < fit->naxes[2]; i++)
		gfit.dft_norm[i] = fit->dft_norm[i];
	gfit.dft_rx = fit->dft_rx;
	gfit.dft_ry = fit->dft_ry;
}

void FFTD(fits *fit, fits *x, fits *y, int type_order, int layer) {
	WORD *xbuf = x->pdata[layer];
	WORD *ybuf = y->pdata[layer];
	WORD *gbuf = fit->pdata[layer];
	unsigned int i;
	int width = fit->rx, height = fit->ry;
	int nbdata = width * height;

	assert(nbdata);
	fftw_complex *spatial_repr = (fftw_complex*) fftw_malloc(
			sizeof(fftw_complex) * nbdata);
	fftw_complex *frequency_repr = (fftw_complex*) fftw_malloc(
			sizeof(fftw_complex) * nbdata);

	/* copying image selection into the fftw data */
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static) if(nbdata > 15000)
#endif
	for (i = 0; i < nbdata; i++)
		spatial_repr[i] = (double) gbuf[i];

	/* we run the Fourier Transform */
	fftw_plan p = fftw_plan_dft_2d(width, height, spatial_repr, frequency_repr,
			FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);

	/* we compute modulus and phase */
	double *modulus = calloc(1, nbdata * sizeof(double));
	double *phase = calloc(1, nbdata * sizeof(double));

	fft_to_spectra(fit, frequency_repr, modulus, phase);

	//We normalize the modulus and the phase
	x->dft_norm[layer] = normalisation_spectra(fit, modulus, phase, xbuf, ybuf,
			type_order);
	if (type_order == 0)
		strcpy(x->dft_ord, "CENTERED");
	else
		strcpy(x->dft_ord, "REGULAR");
	strcpy(y->dft_ord, x->dft_ord);

	free(modulus);
	free(phase);
	fftw_destroy_plan(p);
	fftw_free(spatial_repr);
	fftw_free(frequency_repr);
}

void FFTI(fits *fit, fits *xfit, fits *yfit, int type_order, int layer) {
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
			if (type_order == 0)
				change_symmetry(fit, i, j, &x, &y);
			modulus[j * width + i] = (double) xbuf[y * width + x]
					* (xfit->dft_norm[layer]);
			phase[j * width + i] = (double) ybuf[y * width + x]
					* (2 * M_PI / USHRT_MAX_DOUBLE);
			phase[j * width + i] -= M_PI;
		}
	}

	fftw_complex* spatial_repr = (fftw_complex*) fftw_malloc(
			sizeof(fftw_complex) * nbdata);
	fftw_complex* frequency_repr = (fftw_complex*) fftw_malloc(
			sizeof(fftw_complex) * nbdata);

	fft_to_fr(fit, frequency_repr, modulus, phase);

	fftw_plan p = fftw_plan_dft_2d(width, height, frequency_repr, spatial_repr,
			FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);

	for (i = 0; i < nbdata; i++) {
		double pxl = creal(spatial_repr[i]) / nbdata;
		gbuf[i] = round_to_WORD(pxl);
	}
	delete_selected_area();

	free(modulus);
	free(phase);
	fftw_destroy_plan(p);
	fftw_free(spatial_repr);
	fftw_free(frequency_repr);
}

// idle function executed at the end of the fourier_transform processing
gboolean end_fourier_transform(gpointer p) {
	struct fft_data *args = (struct fft_data *) p;
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

	siril_log_color_message(_("Fourier Transform: processing...\n"), "red");
	gettimeofday(&t_start, NULL);

	//type must be either "ffti" or "fftd"
	switch (args->type[3]) {
	default:
	case 'd':
	case 'D':
		/* We transform the image in a squared picture */
		if (args->fit->rx != args->fit->ry) {
			int size = max(width, height);
			new_fit_image(&wfit[0], size, size, args->fit->naxes[2]);
			for (chan = 0; chan < args->fit->naxes[2]; chan++) {
				from[chan] = args->fit->pdata[chan];
				to[chan] = wfit[0].pdata[chan];
				memcpy(to[chan], from[chan], ndata * sizeof(WORD));
			}
			copyfits(&wfit[0], args->fit, CP_ALLOC | CP_FORMAT | CP_COPYA, 0);
		}
		/* ******************************************* */
		new_fit_image(&wfit[1], args->fit->rx, args->fit->ry,
				args->fit->naxes[2]);
		new_fit_image(&wfit[2], args->fit->rx, args->fit->ry,
				args->fit->naxes[2]);
		for (chan = 0; chan < args->fit->naxes[2]; chan++)
			FFTD(args->fit, &wfit[1], &wfit[2], args->type_order, chan);
		/* we save the original size in the FITS HEADER */
		wfit[1].dft_rx = wfit[2].dft_rx = width;
		wfit[1].dft_ry = wfit[2].dft_ry = height;
		strcpy(wfit[1].dft_type, "SPECTRUM");
		savefits(args->modulus, &wfit[1]);
		strcpy(wfit[2].dft_type, "PHASE");
		savefits(args->phase, &wfit[2]);

		/* We display the modulus on screen */
		copyfits(&wfit[1], &gfit, CP_ALLOC | CP_FORMAT | CP_COPYA, 0);

		/* we copy the header informations */
		save_dft_information_in_gfit(&wfit[1]);
		break;
	case 'i':
	case 'I':
		if (readfits(args->modulus, &wfit[0], NULL)) {
			gdk_threads_add_idle(end_fourier_transform, args);
			return GINT_TO_POINTER(1);
		}
		if (readfits(args->phase, &wfit[1], NULL)) {
			gdk_threads_add_idle(end_fourier_transform, args);
			return GINT_TO_POINTER(1);
		}
		if (wfit[0].dft_ord[0] == 'C')		// CENTERED
			args->type_order = 0;
		else if (wfit[0].dft_ord[0] == 'R')	// REGULAR
			args->type_order = 1;
		else {
			siril_log_message(_("There is something wrong in your files\n"));
			gdk_threads_add_idle(end_fourier_transform, args);
			return GINT_TO_POINTER(1);
		}
		new_fit_image(&wfit[2], wfit[0].rx, wfit[0].ry, wfit[0].naxes[2]);
		for (chan = 0; chan < args->fit->naxes[2]; chan++)
			FFTI(&wfit[2], &wfit[0], &wfit[1], args->type_order, chan);
		new_fit_image(args->fit, wfit[0].dft_rx, wfit[0].dft_ry,
				wfit[0].naxes[2]);
		for (chan = 0; chan < args->fit->naxes[2]; chan++) {
			from[chan] = wfit[2].pdata[chan];
			to[chan] = args->fit->pdata[chan];
			memcpy(to[chan], from[chan],
					wfit[0].dft_rx * wfit[0].dft_ry * sizeof(WORD));
		}
	}
	clearfits(&wfit[0]);
	clearfits(&wfit[1]);
	clearfits(&wfit[2]);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	gdk_threads_add_idle(end_fourier_transform, args);

	return GINT_TO_POINTER(0);
}
