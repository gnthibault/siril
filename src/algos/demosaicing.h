#ifndef _DEMOSAICING_H
#define _DEMOSAICING_H

struct split_cfa_data {
	fits *fit;
	sequence *seq;
	const gchar *seqEntry;
};

WORD *debayer_buffer(WORD *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern);
int debayer(fits*, interpolation_method, sensor_pattern pattern);

#ifdef __cplusplus
extern "C" {
#endif
int retrieve_Bayer_pattern(fits *fit, sensor_pattern *pattern);
WORD *debayer_buffer_superpixel_ushort(WORD *buf, int *width, int *height, sensor_pattern pattern);
float *debayer_buffer_superpixel_float(float *buf, int *width, int *height, sensor_pattern pattern);
#ifdef __cplusplus
}
#endif

void get_debayer_area(const rectangle *area, rectangle *debayer_area,
		const rectangle *image_area, int *debayer_offset_x,
		int *debayer_offset_y);
int extractHa_ushort(fits *in, fits *Ha, sensor_pattern pattern);
int extractHa_float(fits *in, fits *Ha, sensor_pattern pattern);
void apply_extractHa_to_sequence(struct split_cfa_data *split_cfa_args);
int extractHaOIII_ushort(fits *in, fits *Ha, fits *OIII, sensor_pattern pattern);
int extractHaOIII_float(fits *in, fits *Ha, fits *OIII, sensor_pattern pattern);
void apply_extractHaOIII_to_sequence(struct split_cfa_data *split_cfa_args);
int split_cfa_ushort(fits *in, fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3);
int split_cfa_float(fits *in, fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3);
void apply_split_cfa_to_sequence(struct split_cfa_data *split_cfa_args);

#ifdef __cplusplus
extern "C" {
#endif
/* from demosaicing_rtp.cpp */
WORD *debayer_buffer_new_ushort(WORD *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, unsigned int xtrans[6][6]);

float *debayer_buffer_new_float(float *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, unsigned int xtrans[6][6]);
#ifdef __cplusplus
}
#endif

#endif
