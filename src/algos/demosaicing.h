#ifndef _DEMOSAICING_H
#define _DEMOSAICING_H

WORD *debayer_buffer(WORD *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, int xtrans[6][6]);
int debayer(fits*, interpolation_method, gboolean);
void get_debayer_area(const rectangle *area, rectangle *debayer_area,
		const rectangle *image_area, int *debayer_offset_x,
		int *debayer_offset_y);

#endif
