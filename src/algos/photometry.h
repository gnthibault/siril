#ifndef SRC_ALGOS_PHOTOMETRY_H_
#define SRC_ALGOS_PHOTOMETRY_H_

#include <gsl/gsl_matrix.h>

struct photometry_struct {
	double mag; // magnitude
	double s_mag; // magnitude uncertainty
	gboolean valid; // TRUE if no pixel outside of the range
};
typedef struct photometry_struct photometry;

photometry *getPhotometryData(gsl_matrix* z, fitted_PSF *psf, gboolean verbose);
void initialize_photometric_param();

#endif /* SRC_ALGOS_PHOTOMETRY_H_ */
