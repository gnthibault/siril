/*
 * photometry.h
 *
 *  Created on: 9 mars 2017
 *      Author: cyril
 */

#ifndef SRC_ALGOS_PHOTOMETRY_H_
#define SRC_ALGOS_PHOTOMETRY_H_

#include <gsl/gsl_matrix.h>

struct photometry_struct {
	double mag; // magnitude
	double s_pxl; // variance of mean sky level
	double s_mag; // magnitude uncertainty
	double Int; // Intensity
	double B_mean; // background value
	size_t n_sky; // number of pixels in sky
	size_t area; // number of pixels in aperture
};
typedef struct photometry_struct photometry;

photometry *getPhotometricData(gsl_matrix* z, fitted_PSF *psf);

#endif /* SRC_ALGOS_PHOTOMETRY_H_ */
