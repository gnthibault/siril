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

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_version.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/OS_utils.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "algos/statistics.h"
#include "algos/geometry.h"
#include "algos/sorting.h"
#include "gui/callbacks.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "gui/dialogs.h"
#include "registration/registration.h"	// for mouse_status
#include "background_extraction.h"


#define NPARAM_POLY4 15		// Number of parameters used with 4rd order
#define NPARAM_POLY3 10		// Number of parameters used with 3rd order
#define NPARAM_POLY2 6		// Number of parameters used with 2nd order
#define NPARAM_POLY1 3		// Number of parameters used with 1nd order

#define SAMPLE_SIZE 25

//C contains background function
#define C(i) (gsl_vector_get(c,(i)))

static double poly_4(gsl_vector *c, double x, double y) {
	double value = C(0) + C(1) * x + C(2) * y + (x * x) * C(3) + (x * y) * C(4)
			+ (y * y) * C(5) + (x * x) * x * C(6) + (x * x) * y * C(7)
			+ x * (y * y) * C(8)+ (y * y) * y * C(9) + (x * x) * (x * x) * C(10)
			+ (x * x) * (x * y) * C(11) + (x * x) * (y * y) * C(12)
			+ (x * y) * (y * y) * C(13) + (y * y) * (y * y) * C(14);

	return (value);
}

static double poly_3(gsl_vector *c, double x, double y) {
	double value = C(0) + C(1) * x + C(2) * y + (x * x) * C(3) + (y * x) * C(4)
			+ (y * y) * C(5) + (x * x) * x * C(6) + (x * x) * y * C(7) + x * (y * y) * C(8)
			+ (y * y) * y * C(9);

	return (value);
}

static double poly_2(gsl_vector *c, double x, double y) {
	double value = C(0) + C(1) * x + C(2) * y + C(3) * (x * x) + C(4) * (y * x)
			+ C(5) * (y * y);

	return (value);
}

static double poly_1(gsl_vector *c, double x, double y) {
	double value = C(0) + C(1) * x + C(2) * y;

	return (value);
}

static double *computeBackground(GSList *list, int channel, unsigned int width, unsigned int height, poly_order order, gchar **err) {
	size_t k = 0;
	double chisq, pixel;
	double row, col;
	gsl_matrix *J, *cov;
	gsl_vector *y, *w, *c;
	GSList *l;

	guint n = g_slist_length(list);

	int nbParam;
	switch (order) {
	case POLY_1:
		nbParam = NPARAM_POLY1;
		break;
	case POLY_2:
		nbParam = NPARAM_POLY2;
		break;
	case POLY_3:
		nbParam = NPARAM_POLY3;
		break;
	case POLY_4:
	default:
		nbParam = NPARAM_POLY4;
	}

	if (n < nbParam) {
		*err = siril_log_message(_("There are not enough background samples. "
				"The background to be extracted cannot be computed.\n"));
		return NULL;
	}

	// J is the Jacobian
	// y contains data (pixel intensity)
	J = gsl_matrix_calloc(n, nbParam);
	y = gsl_vector_calloc(n);
	w = gsl_vector_calloc(n);
	c = gsl_vector_calloc(nbParam);
	cov = gsl_matrix_calloc(nbParam, nbParam);

	for (l = list; l; l = l->next) {
		background_sample *sample = (background_sample *) l->data;

		col = sample->position.x;
		row = sample->position.y;
		pixel = sample->median[channel];
		// here, it is a bit sketchy in the sense that if there is not value to report in a box (because the threshold is too
		// low for example, then I just skip the initialization of J and y. gsl automatically discard the non assigned values
		// during the minimization. I tested it with Matlab and it works fine. The results agree.
		if (pixel < 0)
			continue;

		gsl_matrix_set(J, k, 0, 1.0);
		gsl_matrix_set(J, k, 1, col);
		gsl_matrix_set(J, k, 2, row);

		if (order != POLY_1) {
			gsl_matrix_set(J, k, 3, col * col);
			gsl_matrix_set(J, k, 4, col * row);
			gsl_matrix_set(J, k, 5, row * row);
		}

		if (order == POLY_3 || order == POLY_4) {
			gsl_matrix_set(J, k, 6, col * col * col);
			gsl_matrix_set(J, k, 7, col * col * row);
			gsl_matrix_set(J, k, 8, col * row * row);
			gsl_matrix_set(J, k, 9, row * row * row);
		}

		if (order == POLY_4) {
			gsl_matrix_set(J, k, 10, col * col * col * col);
			gsl_matrix_set(J, k, 11, col * col * col * row);
			gsl_matrix_set(J, k, 12, col * col * row * row);
			gsl_matrix_set(J, k, 13, col * row * row * row);
			gsl_matrix_set(J, k, 14, row * row * row * row);
		}

		gsl_vector_set(y, k, pixel);
		gsl_vector_set(w, k, 1.0);

		k++;
	}

	// Must turn off error handler or it aborts on error
	gsl_set_error_handler_off();

	gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, nbParam);
	int status = gsl_multifit_wlinear(J, w, y, c, cov, &chisq, work);
	if (status != GSL_SUCCESS) {
		*err = siril_log_message("GSL multifit error: %s\n", gsl_strerror(status));
		gsl_matrix_free(J);
		gsl_vector_free(y);
		gsl_vector_free(w);
		gsl_vector_free(c);
		gsl_matrix_free(cov);
		return NULL;
	}

	// Calculation of the background with the same dimension that the input matrix.
	double *background = malloc(height * width * sizeof(double));
	if (!background) {
		PRINT_ALLOC_ERR;
		*err = _("Out of memory - aborting");
		return NULL;
	}

	for (unsigned int i = 0; i < height; i++) {
		for (unsigned int j = 0; j < width; j++) {
			switch (order) {
			case POLY_1:
				pixel = poly_1(c, (double) j, (double) i);
				break;
			case POLY_2:
				pixel = poly_2(c, (double) j, (double) i);
				break;
			case POLY_3:
				pixel = poly_3(c, (double) j, (double) i);
				break;
			default:
			case POLY_4:
				pixel = poly_4(c, (double) j, (double) i);
			}
			background[j + i * width] = pixel;
		}
	}

	/* free memory */
	gsl_multifit_linear_free(work);
	gsl_matrix_free(J);
	gsl_vector_free(y);
	gsl_vector_free(w);
	gsl_vector_free(c);
	gsl_matrix_free(cov);

	return background;
}

static background_sample *get_sample(float *buf, const int xx,
		const int yy, const int w, const int h) {
	int radius, x, y;
	double *data;
	size_t size = SAMPLE_SIZE * SAMPLE_SIZE;
	background_sample *sample = (background_sample *) g_malloc(sizeof(background_sample));
	if (!sample) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	radius = (int) (SAMPLE_SIZE / 2);

	int n = 0;
	data = calloc(size, sizeof(double));
	if (!data) {
		free(sample);
		PRINT_ALLOC_ERR;
		return NULL;
	}
	for (y = yy - radius; y <= yy + radius; y ++) {
		for (x = xx - radius; x <= xx + radius; x ++) {
			if (y >= 0 && y < h) {
				if (x >= 0 && x < w) {
					data[n++] = (double)buf[x + y * w];
				}
			}
		}
	}
	gsl_stats_minmax(&sample->min, &sample->max, data, 1, size);
	sample->mean = gsl_stats_mean(data, 1, size);
	sample->median[RLAYER] = quickmedian_double(data, size);
	sample->median[GLAYER] = sample->median[BLAYER] = sample->median[RLAYER];
	sample->position.x = xx;
	sample->position.y = yy;
	sample->size = SAMPLE_SIZE;
	sample->valid = TRUE;

	free(data);
	return sample;
}

static double get_sample_median(double *buf, const int xx,
		const int yy, const int w, const int h) {
	int radius, x, y, n;
	double *data, median;
	size_t size = SAMPLE_SIZE * SAMPLE_SIZE;

	radius = get_sample_radius();

	n = 0;
	data = calloc(size, sizeof(double));
	if (!data) {
		PRINT_ALLOC_ERR;
		return -1.0;
	}
	for (y = yy - radius; y <= yy + radius; y ++) {
		for (x = xx - radius; x <= xx + radius; x ++) {
			if (y >= 0 && y < h) {
				if (x >= 0 && x < w) {
					data[n++] = buf[x + y * w];
				}
			}
		}
	}
	median = quickmedian_double(data, size);

	free(data);
	return median;
}

static long dither(long max) {
	unsigned long
	// max <= RAND_MAX < ULONG_MAX, so this is okay.
	num_bins = (unsigned long) max + 1, num_rand = (unsigned long) RAND_MAX + 1,
			bin_size = num_rand / num_bins, defect = num_rand % num_bins;

	long x;

	do {
		x = rand();
	}
	// This is carefully written not to overflow
	while (num_rand - defect <= (unsigned long) x);

	// Truncated division is intentional
	return x / bin_size;
}

static double *convert_fits_to_img(fits *fit, int channel, gboolean add_dither) {
	size_t i, n = fit->naxes[0] * fit->naxes[1];
	double *image = malloc(n * sizeof(double));
	if (!image) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	if (add_dither) {
		/* initialize random seed: */
		srand(time(NULL));
	}

	mirrorx(fit, FALSE);
	/*  copy data to new array and normalize pixel data */
	for (i = 0; i < n; i++) {
		if (fit->type == DATA_USHORT) {
			image[i] = (double) fit->pdata[channel][i] / USHRT_MAX_DOUBLE;
		} else if (fit->type == DATA_FLOAT) {
			image[i] = (double) fit->fpdata[channel][i];
		}
		if (add_dither) {
			/* add dithering in order to avoid colour banding */
			double dithering = (dither(999) * 1E-7);
			image[i] += dithering;
		}
	}
	mirrorx(fit, FALSE);
	return image;
}

static float *convert_fits_to_luminance(fits *fit) {
	size_t i, n = fit->naxes[0] * fit->naxes[1];
	/* allocating memory to image */
	float *image = malloc(n * sizeof(float));
	if (!image) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	mirrorx(fit, FALSE);

	for (i = 0; i < n; i++) {
		if (fit->naxes[2] > 1) {
			float r, g, b;
			if (fit->type == DATA_USHORT) {
				r = (float) fit->pdata[RLAYER][i] / USHRT_MAX_SINGLE;
				g = (float) fit->pdata[GLAYER][i] / USHRT_MAX_SINGLE;
				b = (float) fit->pdata[BLAYER][i] / USHRT_MAX_SINGLE;
			} else if (fit->type == DATA_FLOAT) {
				r = fit->fpdata[RLAYER][i];
				g = fit->fpdata[GLAYER][i];
				b = fit->fpdata[BLAYER][i];
			} else return NULL;
			image[i] = 0.2126f * r + 0.7152f * g + 0.0722f * b;
		} else {
			if (fit->type == DATA_USHORT) {
				image[i] = (float) fit->pdata[RLAYER][i] / USHRT_MAX_SINGLE;
			} else if (fit->type == DATA_FLOAT) {
				image[i] = fit->fpdata[RLAYER][i];
			}
		}
	}

	mirrorx(fit, FALSE);

	return image;
}

static void convert_img_to_fits(double *image, fits *fit, int channel) {
	size_t i, n = fit->naxes[0] * fit->naxes[1];
	mirrorx(fit, FALSE);

	if (fit->type == DATA_USHORT) {
		WORD *buf = fit->pdata[channel];
		for (i = 0; i < n; i++) {
			buf[i] = round_to_WORD(image[i] * USHRT_MAX_DOUBLE);
		}
	} else if (fit->type == DATA_FLOAT) {
		float *buf = fit->fpdata[channel];
		for (i = 0; i < n; i++) {
			buf[i] = (float)image[i];
		}
	}

	mirrorx(fit, FALSE);
}

static double siril_stats_mad(const float data[], const size_t stride,
		const size_t n, float work[]) {
	float median, mad;
	size_t i;

	/* copy input data to work */
	for (i = 0; i < n; ++i)
		work[i] = (float) data[i * stride];

	/* compute median of input data using double version */
	median = histogram_median_float(work, n, TRUE);

	/* compute absolute deviations from median */
	for (i = 0; i < n; ++i)
		work[i] = fabsf((float) data[i * stride] - median);

	mad = histogram_median_float(work, n, TRUE);

	return mad;
}

static GSList *generate_samples(fits *fit, int nb_per_line, double tolerance, size_t size) {
	int nx = fit->rx;
	int ny = fit->ry;
	int dist, starty, startx;
	unsigned int x, y;
	float median, mad0, *work;
	size_t radius;
	size_t n = fit->naxes[0] * fit->naxes[1];
	GSList *list = NULL;

	float *image = convert_fits_to_luminance(fit);

	work = malloc(n * sizeof(float));
	if (!work) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	dist = (int) (nx / nb_per_line);
	radius = size / 2;
	startx = ((nx - size) % dist) / 2;
	starty = ((ny - size) % dist) / 2;
	mad0 = siril_stats_mad(image, 1, nx * ny, work);
	median = histogram_median_float(image, nx * ny, TRUE);

	for (y = starty; y <= ny - radius; y = y + dist) {
		for (x = startx; x <= nx - radius; x = x + dist) {
			background_sample *sample = get_sample(image, x, y, nx, ny);
			if (sample->median[RLAYER] > 0.0
					&& sample->median[RLAYER] <= (mad0 * exp(tolerance)) + median) {
				list = g_slist_prepend(list, sample);
			} else {
				g_free(sample);
			}
		}
	}
	list = g_slist_reverse(list);
	free(image);
	free(work);

	return list;
}

static GSList *update_median_for_rgb_samples(GSList *orig, fits *fit) {
	GSList *list;
	int nx = fit->rx;
	int ny = fit->ry;
	int channel;
	double *rgb[3];

	rgb[RLAYER] = convert_fits_to_img(fit, RLAYER, FALSE);
	rgb[GLAYER] = convert_fits_to_img(fit, GLAYER, FALSE);
	rgb[BLAYER] = convert_fits_to_img(fit, BLAYER, FALSE);

	for (list = orig; list; list = list->next) {
		background_sample *sample = (background_sample *)list->data;
		for (channel = 0; channel < fit->naxes[2]; channel++) {
			sample->median[channel] = get_sample_median(rgb[channel],
					sample->position.x, sample->position.y, nx, ny);
		}
	}

	free(rgb[RLAYER]);
	free(rgb[GLAYER]);
	free(rgb[BLAYER]);

	return orig;
}

static poly_order get_poly_order() {
	GtkComboBox *combo_box_poly_order;

	combo_box_poly_order = GTK_COMBO_BOX(lookup_widget("box_background_order"));
	return gtk_combo_box_get_active(combo_box_poly_order);
}

static int get_correction_type() {
	GtkComboBox *combo_box_correction;

	combo_box_correction = GTK_COMBO_BOX(lookup_widget("box_background_correction"));
	return gtk_combo_box_get_active(combo_box_correction);
}

static int get_nb_samples_per_line() {
	GtkSpinButton *nb_samples = GTK_SPIN_BUTTON(lookup_widget("spin_background_nb_samples"));

	return gtk_spin_button_get_value_as_int(nb_samples);
}

static double get_tolerance_value() {
	GtkRange *tol = GTK_RANGE(lookup_widget("scale_background_nb_samples"));

	return gtk_range_get_value(tol);
}

static void remove_gradient(double *img, double *background, size_t ndata, int type) {
	size_t i;
	double mean;

	mean = gsl_stats_mean(img, 1, ndata);

	switch (type) {
	default:
	case 0: // Subtraction
		for (i = 0; i < ndata; i++) {
			img[i] -= background[i];
			img[i] += mean;
		}
		break;
	case 1: // Division
		for (i = 0; i < ndata; i++) {
			img[i] /= background[i];
			img[i] *= mean;
		}
	}
}

/************* PUBLIC FUNCTIONS *************/

int get_sample_radius() {
	return (int) (SAMPLE_SIZE / 2);
}

void free_background_sample_list(GSList *list) {
	if (list == NULL) return;
	g_slist_free_full(list, g_free);
}

GSList *add_background_sample(GSList *orig, fits *fit, point pt) {
	GSList *list;
	int nx = fit->rx;
	int ny = fit->ry;
	float *image;

	image = convert_fits_to_luminance(fit);

	list = orig;

	background_sample *sample = get_sample(image, pt.x, pt.y, nx, ny);
	list = g_slist_append(list, sample);

	free(image);

	return list;
}

GSList *remove_background_sample(GSList *orig, fits *fit, point pt) {
	GSList *list;
	float *image;

	image = convert_fits_to_luminance(fit);

	for (list = orig; list; list = list->next) {
		background_sample *sample;
		double radius;
		double dx;
		double dy;

		sample = (background_sample *)list->data;
		dx = pt.x - sample->position.x;
		dy = pt.y - sample->position.y;
		radius = sqrt(dx * dx + dy * dy);

		if (radius <= sample->size * 2) {
			orig = g_slist_remove(orig, sample);
			g_free((background_sample *) sample);
			break;
		}
	}
	free(image);

	return orig;
}

void generate_background_samples(int nb_of_samples, double tolerance) {
	set_cursor_waiting(TRUE);

	free_background_sample_list(com.grad_samples);
	com.grad_samples = generate_samples(&gfit, nb_of_samples, tolerance, SAMPLE_SIZE);
	if (gfit.naxes[2] > 1) {
		com.grad_samples = update_median_for_rgb_samples(com.grad_samples,
				&gfit);
	}

	redraw(com.cvport, REMAP_ALL);
	set_cursor_waiting(FALSE);
}

void remove_gradient_from_image(int correction, poly_order degree) {
	double *background, *image[3] = {0};
	gchar *error;

	for (int channel = 0; channel < gfit.naxes[2]; channel++) {
		/* compute background */
		image[channel] = convert_fits_to_img(&gfit, channel, TRUE);
		background = computeBackground(com.grad_samples, channel, gfit.rx, gfit.ry, degree, &error);
		if (background == NULL) {
			if (error && !com.script) {
				siril_message_dialog(GTK_MESSAGE_ERROR, _("Not enough samples."), error);
			}

			free(image[channel]);
			set_cursor_waiting(FALSE);
			return;
		}
		/* remove background */
		const char *c_name = vport_number_to_name(channel);
		siril_log_message(_("Background extraction from channel %s.\n"), c_name);
		remove_gradient(image[channel], background,	gfit.naxes[0] * gfit.naxes[1], correction);
		convert_img_to_fits(image[channel], &gfit, channel);

		/* free memory */
		free(image[channel]);
		free(background);
	}
}

/** Apply for sequence **/

static int background_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_) {
	struct background_data *b_args = (struct background_data*) args->user;

	double *background, *image[3] = {0};
	gchar *error;

	GSList *samples = generate_samples(fit, b_args->nb_of_samples, b_args->tolerance, SAMPLE_SIZE);
	if (samples == NULL) {
		return 1;
	}

	if (fit->naxes[2] > 1) {
		samples = update_median_for_rgb_samples(samples, fit);
	}

	for (int channel = 0; channel < fit->naxes[2]; channel++) {
		/* compute background */
		image[channel] = convert_fits_to_img(fit, channel, TRUE);
		background = computeBackground(samples, channel, fit->rx, fit->ry, b_args->degree, &error);
		if (background == NULL) {
			if (error) {
				siril_log_message(error);
			}

			set_cursor_waiting(FALSE);
			free(image[channel]);
			free_background_sample_list(samples);
			return 1;
		}
		/* remove background */
		remove_gradient(image[channel], background, fit->naxes[0] * fit->naxes[1], b_args->correction);
		convert_img_to_fits(image[channel], fit, channel);

		/* free memory */
		free(image[channel]);
		free(background);
	}
	free_background_sample_list(samples);

	return 0;
}
void apply_background_extraction_to_sequence(struct background_data *background_args) {
	struct generic_seq_args *args = malloc(sizeof(struct generic_seq_args));
	args->seq = background_args->seq;
	args->force_float = FALSE;
	args->partial_image = FALSE;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = background_args->seq->selnum;
	args->prepare_hook = ser_prepare_hook;
	args->finalize_hook = ser_finalize_hook;
	args->save_hook = NULL;
	args->image_hook = background_image_hook;
	args->idle_function = NULL;
	args->stop_on_error = FALSE;
	args->description = _("Background Extraction");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->new_seq_prefix = background_args->seqEntry;
	args->load_new_sequence = TRUE;
	args->force_ser_output = FALSE;
	args->user = background_args;
	args->already_in_a_thread = FALSE;
	args->parallel = TRUE;

	background_args->fit = NULL;	// not used here

	start_in_new_thread(generic_sequence_worker, args);
}

/************* CALLBACKS *************/

void on_menuitem_background_extraction_activate(GtkMenuItem *menuitem,
		gpointer user_data) {
	siril_open_dialog("background_extraction_dialog");
}

void on_background_generate_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	int nb_of_samples;
	double tolerance;

	nb_of_samples = get_nb_samples_per_line();
	tolerance = get_tolerance_value();

	generate_background_samples(nb_of_samples, tolerance);

	redraw(com.cvport, REMAP_ALL);
	set_cursor_waiting(FALSE);
}

void on_background_clear_all_clicked(GtkButton *button, gpointer user_data) {
	free_background_sample_list(com.grad_samples);
	com.grad_samples = NULL;

	redraw(com.cvport, REMAP_ALL);
	set_cursor_waiting(FALSE);
}

void on_background_ok_button_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);

	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkBkgSeq")))
			&& sequence_is_loaded()) {
		struct background_data *args = malloc(sizeof(struct background_data));

		args->nb_of_samples = get_nb_samples_per_line();
		args->tolerance = get_tolerance_value();
		args->correction = get_correction_type();
		args->degree = get_poly_order();
		if (args->degree > POLY_1) {
			int confirm = siril_confirm_dialog(_("Polynomial order seems too high."),
					_("You are about to process a sequence of preprocessed files with "
							"a polynomial degree greater than 1. This is unlikely because such "
							"gradients are often linear and a correction with a polynomial "
							"function of degree 1 is probably enough. Click OK to confirm or "
							"cancel to change your mind."));
			if (!confirm) {
				free(args);
				set_cursor_waiting(FALSE);
				return;
			}
		}

		args->seqEntry = gtk_entry_get_text(GTK_ENTRY(lookup_widget("entryBkgSeq")));
		if (args->seqEntry && args->seqEntry[0] == '\0')
			args->seqEntry = "bkg_";
		args->seq = &com.seq;
		apply_background_extraction_to_sequence(args);
	} else {
		if (com.grad_samples == NULL) {
			set_cursor_waiting(FALSE);
			return;
		}

		int correction = get_correction_type();
		poly_order degree = get_poly_order();
		undo_save_state(&gfit, "Processing: Background extraction (Correction: %s)",
				correction ? "Division" : "Subtraction");
		remove_gradient_from_image(correction, degree);

		invalidate_stats_from_fit(&gfit);
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);

	}
	set_cursor_waiting(FALSE);
}

void on_background_close_button_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("background_extraction_dialog");
}

void on_background_extraction_dialog_hide(GtkWidget *widget, gpointer user_data) {
	free_background_sample_list(com.grad_samples);
	com.grad_samples = NULL;
	mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	redraw(com.cvport, REMAP_ALL);
}

void on_background_extraction_dialog_show(GtkWidget *widget, gpointer user_data) {
	mouse_status = MOUSE_ACTION_DRAW_SAMPLES;
}
