#ifndef _COLORS_H_
#define _COLORS_H_

/* scnr data from GUI */
struct scnr_data {
	fits *fit;
	int type;
	double amount;
	gboolean preserve;
};

/* extract_channels data from GUI */
struct extract_channels_data {
	gboolean process;
	fits *fit;
	const char *channel[3];
	int type;
	const char* str_type;
};

/* color saturation data from GUI */
struct enhance_saturation_data {
	fits *fit;
	double coeff, h_min, h_max;
	gboolean preserve;
};

void rgb_to_hsl(double, double, double, double *, double *, double *);
void hsl_to_rgb(double, double, double, double *, double *, double *);
void rgb_to_hsv(double, double, double, double *, double *, double *);
void hsv_to_rgb(double, double, double, double *, double *, double *);
void rgb_to_xyz(double, double, double, double *, double *, double *);
void xyz_to_LAB(double, double, double, double *, double *, double *);
void LAB_to_xyz(double, double, double, double *, double *, double *);
void xyz_to_rgb(double, double, double, double *, double *, double *);

gpointer extract_channels(gpointer p);
gpointer enhance_saturation(gpointer p);
gpointer scnr(gpointer p);

void initialize_calibration_interface();

#endif
