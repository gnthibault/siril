#ifndef PSF_H_
#define PSF_H_

#include <gsl/gsl_matrix.h>
#include <algos/photometry.h>

//in siril.h: typedef struct fwhm_struct psf_star;

struct fwhm_struct {
	double B; /* average sky background value */
	double A; /* amplitude */
	double x0, y0; /* coordinates of the peak */
	double sx, sy; /* Size of the fitted function on the x and y axis in PSF coordinates */
	double fwhmx, fwhmy; /* FWHM in x and y axis */
	double fwhmx_arcsec, fwhmy_arcsec; /* FWHM in x and y axis in arc second */
	double angle; /* angle of the axis x,y with respect to the image's */
	double mag; /* magnitude of the star : this parameter is not fitted but calculated with the vector G and the parameter B */
	double s_mag; /* error on the magnitude */
	double SNR; /* SNR of the star */
	photometry *phot; /* photometry data */
	gboolean phot_is_valid; /* valid if computed by photometry and no saturated pixel detected */
	double xpos, ypos; /* position of the star in the image, not set by Minimization */
	double rmse; /* RMSE of the minimization */
	double BV; /* only use in BV calibration */

	/* uncertainties */
	double B_err;
	double A_err;
	double x_err, y_err;
	double sx_err, sy_err;
	double ang_err;
	int layer;
	char* units;
};

struct PSF_data {
	size_t n;
	double *y;
	double *sigma;
	size_t NbRows;
	size_t NbCols;
	double rmse;
};

double psf_get_fwhm(fits *fit, int layer, rectangle *selection, double *roundness);
psf_star *psf_get_minimisation(fits *, int, rectangle *, gboolean, gboolean, gboolean);
psf_star *psf_global_minimisation(gsl_matrix *, double, gboolean, gboolean, gboolean);
void psf_display_result(psf_star *, rectangle *);
void fwhm_to_arcsec_if_needed(fits*, psf_star*);
void fwhm_to_pixels(psf_star *result);
gboolean get_fwhm_as_arcsec_if_possible(psf_star *star, double *fwhmx, double *fwhmy, char **unit);
double convert_single_fwhm_to_pixels(double fwhm, double s);
gboolean convert_single_fwhm_to_arcsec_if_possible(double fwhm, double bin, double px_size, double flenght, double *result);

psf_star *new_psf_star();
psf_star *duplicate_psf(psf_star *);
void free_psf(psf_star *psf);

#endif
