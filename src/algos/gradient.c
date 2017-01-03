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
#include <assert.h>
#include <string.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#include "core/siril.h"
#include "core/proto.h"
#include "gui/callbacks.h"
#include "algos/gradient.h"
#include "registration/registration.h"	// for mouse_status

#define NPARAM_POLY4 15		// Number of parameters used with 4rd order
#define NPARAM_POLY3 10		// Number of parameters used with 3rd order
#define NPARAM_POLY2 6		// Number of parameters used with 2nd order
#define NPARAM_POLY1 3		// Number of parameters used with 1nd order

//C contains background function
#define C(i) (gsl_vector_get(c,(i)))

static double poly_4(gsl_vector * c, int x, int y) {
	double value = C(0) * 1 + C(1) * x + C(2) * y + C(3) * x * x + C(4) * y * x
			+ C(5) * y * y +
			C(6) * x * x * x + C(7) * x * x * y + C(8) * x * y * y
			+ C(9) * y * y * y +
			C(10) * x * x * x * x + C(11) * x * x * x * y
			+ C(12) * x * x * y * y + C(13) * x * y * y * y
			+ C(14) * y * y * y * y;

	return (value < 0 ? 0.0 : value);
}

static double poly_3(gsl_vector * c, int x, int y) {
	double value = C(0) * 1 + C(1) * x + C(2) * y + C(3) * x * x + C(4) * y * x
			+ C(5) * y * y +
			C(6) * x * x * x + C(7) * x * x * y + C(8) * x * y * y
			+ C(9) * y * y * y;

	return (value < 0 ? 0.0 : value);
}

static double poly_2(gsl_vector * c, int x, int y) {
	double value = C(0) * 1 + C(1) * x + C(2) * y + C(3) * x * x + C(4) * y * x
			+ C(5) * y * y;

	return (value < 0 ? 0.0 : value);
}

static double poly_1(gsl_vector * c, int x, int y) {
	double value = C(0) * 1 + C(1) * x + C(2) * y;

	return (value < 0 ? 0.0 : value);
}


static int buildBoxesAutomatically(gsl_vector *MatR, newBackground *bkg, int layer) {
	size_t i, k;
	size_t row, col, inc_row, inc_col;
	size_t inc = 0;
	double tmpRow, tmpCol;
	double tolerance = bkg->tolerance;
	double deviation = bkg->deviation;
	double unbalance = bkg->unbalance;
	size_t box = bkg->box, midbox = box * 0.5;
	size_t boxPerRow = bkg->boxPerRow;
	size_t boxPerCol = bkg->boxPerCol;
	size_t height = bkg->row;
	size_t width = bkg->col;

	assert(box > 0);

	bkg->NbBoxes = boxPerRow * boxPerCol;
	clearSamples();
	com.grad = malloc(bkg->NbBoxes * sizeof(gradient));

	if (((bkg->order == POLY_1) && (bkg->NbBoxes < NPARAM_POLY1))
			|| ((bkg->order == POLY_2) && (bkg->NbBoxes < NPARAM_POLY2))
			|| ((bkg->order == POLY_3) && (bkg->NbBoxes < NPARAM_POLY3))
			|| ((bkg->order == POLY_4) && (bkg->NbBoxes < NPARAM_POLY4))) {
		return -1;
	}
	gsl_vector *vecRow = gsl_vector_alloc(boxPerCol);
	gsl_vector *vecCol = gsl_vector_alloc(boxPerRow);

	tmpRow = (double) midbox - 1.0;
	tmpCol = (double) midbox - 1.0;

	for (i = 0; i < boxPerCol; i++) {
		gsl_vector_set(vecRow, i, tmpRow);
		tmpRow += (double) ((height - 2 * midbox) / (boxPerCol - 1));
	}
	for (i = 0; i < boxPerRow; i++) {
		gsl_vector_set(vecCol, i, tmpCol);
		tmpCol += (double) ((width - 2 * midbox) / (boxPerRow - 1));
	}

	bkg->meshRow = gsl_vector_alloc(bkg->NbBoxes);
	bkg->meshCol = gsl_vector_alloc(bkg->NbBoxes);
	bkg->meshVal = gsl_vector_alloc(bkg->NbBoxes);

	for (row = 0; row < boxPerCol; row++) {
		size_t start_row = round(gsl_vector_get(vecRow, row) - midbox + 1);
		for (col = 0; col < boxPerRow; col++) {
			double *data_box = calloc(box * box, sizeof(double));
			size_t start_col = round(gsl_vector_get(vecCol, col) - midbox + 1);
			for (inc_row = 0, k = 0; inc_row < box; inc_row++) {
				for (inc_col = 0; inc_col < box; inc_col++) {
					data_box[k++] = gsl_vector_get(MatR,
							((size_t) start_row + inc_row) * width + start_col
									+ inc_col);
				}
			}
			double sigma = gsl_stats_sd(data_box, 1, k);
			gsl_sort(data_box, 1, k);
			double median = gsl_stats_median_from_sorted_data(data_box, 1, k);

			for (inc_row = 0, k = 0; inc_row < box; inc_row++) {
				for (inc_col = 0; inc_col < box; inc_col++) {
					double current_pixel = gsl_vector_get(MatR,
							(size_t) (start_row + inc_row) * width + start_col
									+ inc_col);
					if (current_pixel > (tolerance * sigma + median)) {
						gsl_vector_set(MatR,
								(size_t) (start_row + inc_row) * width
										+ start_col + inc_col, median);
					}
					data_box[k++] = gsl_vector_get(MatR,
							(size_t) (start_row + inc_row) * width + start_col
									+ inc_col);
				}
			}
			gsl_sort(data_box, 1, k);
			double value = gsl_stats_median_from_sorted_data(data_box, 1, k);

			gsl_vector_set(bkg->meshVal, inc, value);

			/* Drawing boxes */
			com.grad[inc].centre.x = gsl_vector_get(vecCol, col) + midbox;
			com.grad[inc].centre.y = height - gsl_vector_get(vecRow, row)
					+ midbox;
			gsl_vector_set(bkg->meshRow, inc, gsl_vector_get(vecRow, row));
			gsl_vector_set(bkg->meshCol, inc, gsl_vector_get(vecCol, col));
			inc++;
			free(data_box);
		}
	}

	double *data = calloc(bkg->NbBoxes, sizeof(double));
	for (i = 0; i < bkg->NbBoxes; i++)
		data[i] = gsl_vector_get(bkg->meshVal, i);

	gsl_sort(data, 1, bkg->NbBoxes);
	double median = gsl_stats_median_from_sorted_data(data, 1, bkg->NbBoxes);
	double sigma = gsl_stats_sd(data, 1, bkg->NbBoxes);

	for (i = 0; i < bkg->NbBoxes; i++) {
		double pixel = gsl_vector_get(bkg->meshVal, i);
		if (((pixel - median) / sigma > deviation)
				|| ((median - pixel) / sigma > (deviation * unbalance)))
			gsl_vector_set(bkg->meshVal, i, -1);
		com.grad[i].boxvalue[layer] = gsl_vector_get(bkg->meshVal, i);
	}
	free(data);
	return 0;
}

static gsl_matrix *computeBackground(newBackground *bkg) {
	size_t n, i, j;
	size_t inc = 0;
	double chisq, pixel_value, tmpRow, tmpCol;
	size_t height = bkg->row;
	size_t width = bkg->col;
	gsl_matrix *J, *cov;
	gsl_vector *y, *w, *c;

	n = bkg->NbBoxes;

	int nbParam;
	switch (bkg->order) {
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

	// J is the Jacobian
	// y contains data (pixel intensity)
	J = gsl_matrix_calloc(n, nbParam);
	y = gsl_vector_calloc(n);
	w = gsl_vector_calloc(n);
	c = gsl_vector_calloc(nbParam);
	cov = gsl_matrix_calloc(nbParam, nbParam);

	for (inc = 0; inc < n; inc++) {
		tmpCol = gsl_vector_get(bkg->meshCol, inc);
		tmpRow = gsl_vector_get(bkg->meshRow, inc);
		pixel_value = gsl_vector_get(bkg->meshVal, inc);
		// here, it is a bit sketchy in the sense that if there is not value to report in a box (because the threshold is too
		// low for example, then I just skip the initialization of J and y. gsl automatically discard the non assigned values
		// during the minimization. I tested it with Matlab and it works fine. The results agree.
		if (pixel_value < 0)
			continue;

		gsl_matrix_set(J, inc, 0, 1.0);
		gsl_matrix_set(J, inc, 1, tmpCol);
		gsl_matrix_set(J, inc, 2, tmpRow);

		if (bkg->order != POLY_1) {
			gsl_matrix_set(J, inc, 3, tmpCol * tmpCol);
			gsl_matrix_set(J, inc, 4, tmpCol * tmpRow);
			gsl_matrix_set(J, inc, 5, tmpRow * tmpRow);
		}

		if (bkg->order == POLY_3 || bkg->order == POLY_4) {
			gsl_matrix_set(J, inc, 6, tmpCol * tmpCol * tmpCol);
			gsl_matrix_set(J, inc, 7, tmpCol * tmpCol * tmpRow);
			gsl_matrix_set(J, inc, 8, tmpCol * tmpRow * tmpRow);
			gsl_matrix_set(J, inc, 9, tmpRow * tmpRow * tmpRow);
		}

		if (bkg->order == POLY_4) {
			gsl_matrix_set(J, inc, 10, tmpCol * tmpCol * tmpCol * tmpCol);
			gsl_matrix_set(J, inc, 11, tmpCol * tmpCol * tmpCol * tmpRow);
			gsl_matrix_set(J, inc, 12, tmpCol * tmpCol * tmpRow * tmpRow);
			gsl_matrix_set(J, inc, 13, tmpCol * tmpRow * tmpRow * tmpRow);
			gsl_matrix_set(J, inc, 14, tmpRow * tmpRow * tmpRow * tmpRow);
		}

		gsl_vector_set(y, inc, pixel_value);
		gsl_vector_set(w, inc, 1.0);
	}

	gsl_multifit_linear_workspace *work;

	work = gsl_multifit_linear_alloc(n,	nbParam);
	gsl_multifit_wlinear(J, w, y, c, cov, &chisq, work);
	gsl_multifit_linear_free(work);
	gsl_matrix_free(J);
	gsl_vector_free(y);
	gsl_vector_free(w);
	gsl_vector_free(bkg->meshVal);
	gsl_vector_free(bkg->meshRow);
	gsl_vector_free(bkg->meshCol);

// Calculation of the background with the same dimension that the input matrix.
	gsl_matrix *bkgMatrix = gsl_matrix_alloc(height, width);

	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			switch (bkg->order) {
			case POLY_1:
				pixel_value = poly_1(c, j, i);
				break;
			case POLY_2:
				pixel_value = poly_2(c, j, i);
				break;
			case POLY_3:
				pixel_value = poly_3(c, j, i);
				break;
			case POLY_4:
			default:		//default should not happen
				pixel_value = poly_4(c, j, i);
			}
			gsl_matrix_set(bkgMatrix, i, j, pixel_value);
		}
	}
	return bkgMatrix;
}

static int extractBackgroundAuto(fits *imgfit, fits *bkgfit, newBackground *bkg) {
	size_t ndata, i, j;
	WORD *buf = imgfit->pdata[bkg->layer];

	ndata = (size_t) (imgfit->rx * imgfit->ry);
	gsl_vector *MatR = gsl_vector_alloc(ndata);

	for (i = 0; i < ndata; i++) {
		gsl_vector_set(MatR, i, (double) buf[i]);
	}

	com.grad_nb_boxes = bkg->boxPerCol * bkg->boxPerRow;
	com.grad_size_boxes = bkg->box;

	if (buildBoxesAutomatically(MatR, bkg, bkg->layer))
		return 1;		// not enough samples
	gsl_vector_free(MatR);

	gsl_matrix *bkgMatrix = computeBackground(bkg);

	if (imgfit->naxes[2] > 1)
		copyfits(imgfit, bkgfit, CP_ALLOC | CP_FORMAT | CP_EXPAND, bkg->layer);
	else
		copyfits(imgfit, bkgfit, CP_ALLOC | CP_FORMAT | CP_COPYA, 0);

	WORD *tbuf = bkgfit->pdata[bkg->layer];
	for (i = 0; i < bkg->row; i++) {
		for (j = 0; j < bkg->col; j++)
			*tbuf++ = (WORD) gsl_matrix_get(bkgMatrix, i, j);
	}

	siril_log_message(_("Channel #%d: background extraction done.\n"), bkg->layer);
	gsl_matrix_free(bkgMatrix);
	return 0;
}

static int extractBackgroundManual(fits *imgfit, fits *bkgfit, newBackground *bkg) {
	gsl_matrix *bkgMatrix;
	size_t i, j;

	int n = bkg->NbBoxes = com.grad_nb_boxes;

	if (((bkg->order == POLY_1) && (n < NPARAM_POLY1))
			|| ((bkg->order == POLY_2) && (n < NPARAM_POLY2))
			|| ((bkg->order == POLY_3) && (n < NPARAM_POLY3))
			|| ((bkg->order == POLY_4) && (n < NPARAM_POLY4))) {
		return -1;
	}

	bkg->meshRow = gsl_vector_alloc(n);
	bkg->meshCol = gsl_vector_alloc(n);
	bkg->meshVal = gsl_vector_alloc(n);

	for (i = 0; i < n; i++) {
		gsl_vector_set(bkg->meshRow, i,
				imgfit->ry - com.grad[i].centre.y + (bkg->box * 0.5));
		gsl_vector_set(bkg->meshCol, i,
				com.grad[i].centre.x - (bkg->box * 0.5));
		gsl_vector_set(bkg->meshVal, i, com.grad[i].boxvalue[bkg->layer]);
	}

	bkgMatrix = computeBackground(bkg);

	if (imgfit->naxes[2] > 1)
		copyfits(imgfit, bkgfit, CP_ALLOC | CP_FORMAT | CP_EXPAND, bkg->layer);
	else
		copyfits(imgfit, bkgfit, CP_ALLOC | CP_FORMAT | CP_COPYA, 0);

	WORD *tbuf = bkgfit->pdata[bkg->layer];
	for (i = 0; i < bkg->row; i++) {
		for (j = 0; j < bkg->col; j++)
			*tbuf++ = (WORD) gsl_matrix_get(bkgMatrix, i, j);
	}

	siril_log_message(_("Channel #%d: background extraction done.\n"), bkg->layer);
	gsl_matrix_free(bkgMatrix);
	return 0;
}

void clearSamples() {
	if (com.grad) {
		free(com.grad);
		com.grad = NULL;
	}
}

void bkgExtractBackground(fits *fit, gboolean automatic) {
	int layer;
	struct timeval t_start, t_end;
	GtkComboBox *comboBkgPolyOrder;
	GtkSpinButton *spinBkgSizeBox, *spinBkgInterval, *spinBkgTolerance;
	GtkSpinButton *spinBkgDeviation, *spinBkgUnbalance;
	newBackground bkg;

	comboBkgPolyOrder = GTK_COMBO_BOX(gtk_builder_get_object(builder, "combo_polyorder"));
	spinBkgSizeBox = GTK_SPIN_BUTTON(lookup_widget("spinbutton_bkg_sizebox"));
	/* Automatic mode */
	spinBkgInterval = GTK_SPIN_BUTTON(lookup_widget("spinbutton_bkg_Box_sep"));
	spinBkgTolerance = GTK_SPIN_BUTTON(lookup_widget("spinbutton_bkg_tolerance"));
	spinBkgDeviation = GTK_SPIN_BUTTON(lookup_widget("spinbutton_bkg_deviation"));
	spinBkgUnbalance = GTK_SPIN_BUTTON(lookup_widget("spinbutton_bkg_unbalance"));

	siril_log_color_message(_("Background extraction: processing...\n"), "red");
	gettimeofday(&t_start, NULL);

	bkg.order = gtk_combo_box_get_active(comboBkgPolyOrder);
	bkg.box = (size_t) gtk_spin_button_get_value(spinBkgSizeBox) * 2;
	bkg.row = (size_t) gfit.ry;
	bkg.col = (size_t) gfit.rx;

	for (layer = 0; layer < com.uniq->nb_layers; layer++) {
		bkg.layer = layer;
		if (automatic) {
			int interval = (size_t) gtk_spin_button_get_value(spinBkgInterval);

			bkg.tolerance = gtk_spin_button_get_value(spinBkgTolerance);
			bkg.deviation = gtk_spin_button_get_value(spinBkgDeviation);
			bkg.unbalance = gtk_spin_button_get_value(spinBkgUnbalance);

			bkg.boxPerRow = (size_t) ((double) bkg.col / ((double) bkg.box + interval - 1));
			bkg.boxPerCol = (size_t) ((double) bkg.row / ((double) bkg.box + interval - 1));

			if (extractBackgroundAuto(&gfit, fit, &bkg)) {
				siril_log_message(_("Insufficient background samples.\n"));
				return;
			}
		}
		else {
			if (extractBackgroundManual(&gfit, fit, &bkg)) {
				siril_log_message(_("Insufficient background samples.\n"));
				return;
			}
		}
	}
	gtk_widget_set_sensitive(lookup_widget("frame_bkg_tools"), TRUE);
	gtk_widget_set_sensitive(lookup_widget("button_bkg_correct"), TRUE);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
}

static void remove_pixel(double *arr, int i, int N) {
	memmove(&arr[i], &arr[i + 1], (N - i - 1) * sizeof(*arr));
}

double get_value_from_box(fits *fit, point box, size_t size, int layer) {
	double value = -1.0;
	double sigma, median;
	int i, j, stridefrom, data;
	WORD *from;
	rectangle area;
	memset(&area, 0, sizeof(rectangle));
	double *databox = calloc(size * size, sizeof(double));

	area.w = area.h = size;
	area.x = box.x - size / 2.0;
	area.y = box.y - size / 2.0;

	data = size * size;

	// create the matrix with values from the selected rectangle
	from = fit->pdata[layer] + (fit->ry - area.y - area.h) * fit->rx
			+ area.x;
	stridefrom = fit->rx - area.w;

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			databox[j + i * size] = (double) *from;
			from++;
		}
		from += stridefrom;
	}
	sigma = gsl_stats_sd(databox, 1, data);
	gsl_sort(databox, 1, data);
	median = gsl_stats_median_from_sorted_data(databox, 1, data);

	for (i = 0; i < data; i++) {
		double tmp = databox[i];

		if (tmp > (sigma + median) || tmp < (median - sigma)) {
			remove_pixel(databox, i, data);
			data--;
			i--;
		}
	}
	value = gsl_stats_median_from_sorted_data(databox, 1, data);

	free(databox);
	return value;
}

void update_bkg_interface() {
	static GtkToggleButton *bgkManButton = NULL;
	GtkWidget *frame24 = GTK_WIDGET(lookup_widget("frame24"));
	GtkWidget *frame23 = GTK_WIDGET(lookup_widget("frame23"));
	GtkWidget *label44 = GTK_WIDGET(lookup_widget("label44"));
	GtkWidget *spinSeparation = GTK_WIDGET(lookup_widget("spinbutton_bkg_Box_sep"));
	GtkWidget *drawBoxes = GTK_WIDGET(lookup_widget("checkbutton_bkg_boxes"));
	gboolean manButton;

	clearSamples();

	if (bgkManButton == NULL) {
		bgkManButton = GTK_TOGGLE_BUTTON(lookup_widget("bkgButtonManual"));
	}
	manButton = gtk_toggle_button_get_active(bgkManButton);
	if (!manButton) {
		gtk_widget_set_sensitive(frame23, TRUE);
		gtk_widget_set_sensitive(frame24, TRUE);
		gtk_widget_set_sensitive(label44, TRUE);
		gtk_widget_set_sensitive(spinSeparation, TRUE);
		gtk_widget_set_sensitive(drawBoxes, TRUE);
		mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	}
	else {
		gtk_widget_set_sensitive(frame23, FALSE);
		gtk_widget_set_sensitive(frame24, FALSE);
		gtk_widget_set_sensitive(label44, FALSE);
		gtk_widget_set_sensitive(spinSeparation, FALSE);
		gtk_widget_set_sensitive(drawBoxes, FALSE);
		mouse_status = MOUSE_ACTION_DRAW_SAMPLES;
	}
}
