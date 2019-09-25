#ifndef _COLORS_H_
#define _COLORS_H_

/* extract_channels data from GUI */
struct extract_channels_data {
	gboolean process;
	fits *fit;
	const char *channel[3];
	int type;
	const char* str_type;
};

void rgb_to_hsl(double, double, double, double *, double *, double *);
void hsl_to_rgb(double, double, double, double *, double *, double *);
void rgb_to_hsv(double, double, double, double *, double *, double *);
void hsv_to_rgb(double, double, double, double *, double *, double *);
void rgb_to_xyz(double, double, double, double *, double *, double *);
void xyz_to_LAB(double, double, double, double *, double *, double *);
void LAB_to_xyz(double, double, double, double *, double *, double *);
void xyz_to_rgb(double, double, double, double *, double *, double *);
double BV_to_T(double BV);

int equalize_cfa_fit_with_coeffs(fits *fit, double coeff1, double coeff2, int config);

gpointer extract_channels(gpointer p);

void initialize_calibration_interface();

#endif
