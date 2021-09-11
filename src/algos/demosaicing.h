#ifndef _DEMOSAICING_H
#define _DEMOSAICING_H

extern const char *filter_pattern[];
int retrieveBayerPatternFromChar(const char *bayer);

WORD *debayer_buffer(WORD *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, int bit_depth);
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

sensor_pattern get_bayer_pattern(fits *fit);
void clear_Bayer_information(fits *fit);

void get_debayer_area(const rectangle *area, rectangle *debayer_area,
		const rectangle *image_area, int *debayer_offset_x,
		int *debayer_offset_y);

#ifdef __cplusplus
extern "C" {
#endif
/* from demosaicing_rtp.cpp */
WORD *debayer_buffer_new_ushort(WORD *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, unsigned int xtrans[6][6], int bit_depth);

float *debayer_buffer_new_float(float *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, unsigned int xtrans[6][6]);
#ifdef __cplusplus
}
#endif

#endif
