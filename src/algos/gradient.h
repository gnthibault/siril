#ifndef GRAD_H
#define GRAD_H

#include <gsl/gsl_vector.h>

#define NB_MAX_OF_SAMPLES 2000

typedef struct {
	size_t row, col;
	size_t box;
	size_t NbBoxes;
	size_t boxPerRow;
	size_t boxPerCol;
	double tolerance;
	double deviation;
	double unbalance;
	gsl_vector *meshVal; /* mean of the pixel intensity in the box. Pixels above threshold are rejected of the mean calcul */
	gsl_vector *meshCol; /* column coordinates of the box center */
	gsl_vector *meshRow; /* row coordinates of the box center */
	poly_order order;
	int layer;
} newBackground;

void clearSamples();
void bkgExtractBackground(fits *fit, gboolean automatic);
double get_value_from_box(fits *fit, point box, size_t size, int layer);
void update_bkg_interface();

#endif
