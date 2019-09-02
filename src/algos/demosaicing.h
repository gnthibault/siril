#ifndef _DEMOSAICING_H
#define _DEMOSAICING_H

struct split_cfa_data {
	fits *fit;
	sequence *seq;
	const gchar *seqEntry;
};

WORD *debayer_buffer(WORD *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, int xtrans[6][6]);
int debayer(fits*, interpolation_method, gboolean);
void get_debayer_area(const rectangle *area, rectangle *debayer_area,
		const rectangle *image_area, int *debayer_offset_x,
		int *debayer_offset_y);
int split_cfa(fits *in, fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3);
void apply_split_cfa_to_sequence(struct split_cfa_data *split_cfa_args);

#endif
