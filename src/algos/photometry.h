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
	double s_mag; // magnitude uncertainty
};
typedef struct photometry_struct photometry;

photometry *getPhotometryData(gsl_matrix* z, fitted_PSF *psf);

#endif /* SRC_ALGOS_PHOTOMETRY_H_ */
