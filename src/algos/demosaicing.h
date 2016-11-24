#ifndef _DEMOSAICING_H
#define _DEMOSAICING_H

int bayer_Bilinear(const WORD *, WORD *, int, int, sensor_pattern);
int bayer_NearestNeighbor(const WORD *, WORD *, int, int, sensor_pattern);
int bayer_VNG(const WORD *, WORD *, int, int, sensor_pattern);
int bayer_AHD(const WORD *, WORD *, int, int, sensor_pattern);
int super_pixel(const WORD*, WORD*, int, int, sensor_pattern);
WORD *debayer_buffer(WORD *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern);
int debayer(fits*, interpolation_method);
void get_debayer_area(const rectangle *area, rectangle *debayer_area,
		const rectangle *image_area, int *debayer_offset_x,
		int *debayer_offset_y);

#endif
