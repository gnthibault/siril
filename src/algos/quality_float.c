/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2019 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

/** This code comes from PIPP https://sites.google.com/site/astropipp/ */

#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "core/siril.h"
#include "core/proto.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "algos/quality.h"

static float SubSample(float *ptr, int img_wid, int x_size, int y_size);
static float *_smooth_image_float(float *buf, int width, int height);
static double Gradient(float *buf, int width, int height);

// -------------------------------------------------------
// Method to estimate quality.
// Runs on the complete layer and destroys it.
// -------------------------------------------------------
double QualityEstimate_float(fits *fit, int layer) {
	int width = fit->rx;
	int height = fit->ry;
	int x1, y1;
	int subsample, region_w, region_h;
	int i, j, n, x, y, max, x_inc;
	int x_samples, y_samples, y_last;
	float *buffer, *buf, maxp[MAXP];
	float *new_image = NULL;
	double q, dval = 0.0;

	buffer = fit->fpdata[layer];

	// dimensions of the region we want to analyse
	x1 = 0; y1 = 0;
	region_w = width - 1;
	region_h = height - 1;

	// Allocate the intermediate buffer. Will be 16bpp greyscale
	buf = calloc(region_w * region_h, sizeof(float));

	subsample = QSUBSAMPLE_MIN;
	while (subsample <= QSUBSAMPLE_MAX) {
		float* ptr;

		// Number of h & v pixels in subimage
		x_samples = region_w / subsample;
		y_samples = region_h / subsample;

		if (x_samples < 2 || y_samples < 2) {
			break;
		}

		y_last = y1 + (y_samples - 1) * subsample; // second last row of subsampled output
		x_inc = subsample;

		for (i = 0; i < MAXP; ++i) {
			maxp[i] = 0;
		}

		// First row - ignore histo-stretch
		y = y1;
		n = 0;
		ptr = buffer + (y * width + x1);
		for (x = 0; x < x_samples; ++x, ptr += x_inc) {
			buf[n++] = SubSample(ptr, width, subsample, subsample);
		}

		// Rows 1 .. y_last-1 additional histo-stretch code
		for (y += subsample; y < y_last; y += subsample) {
			ptr = buffer + (y * width + x1);
			for (x = 0; x < x_samples; ++x, ptr += x_inc) {
				float v = SubSample(ptr, width, subsample, subsample);

				if (v > maxp[2] && v < 65530) {
					int slot;
					if (v > maxp[0]) {
						slot = 0;
					} else if (v > maxp[1]) {
						slot = 1;
					} else {
						slot = 2;
					}

					for (j = MAXP - 1; j > slot; --j) {
						maxp[j] = maxp[j - 1];
						maxp[j] = v;
					}
				}

				buf[n++] = v;
			}
		}

		// Last row, ignore histo-stretch
		ptr = buffer + (y * width + x1);
		for (x = 0; x < x_samples; ++x, ptr += x_inc) {
			buf[n++] = SubSample(ptr, width, subsample, subsample);
		}

		// Average the bottom half brightest pixels to get the real max
		// to reduce noise effects
		j = MAXP / 2;
		for (i = j, max = 0; i < MAXP; ++i) {
			max += maxp[i];
		}

		// Test idea - reduce quality if histogram peak is lower
		max /= (MAXP - j);

		// 3x3 smoothing
		new_image = _smooth_image_float(buf, x_samples, y_samples);

		q = Gradient(new_image, x_samples, y_samples);

		dval += (q * ((QSUBSAMPLE_MIN * QSUBSAMPLE_MIN)
					/ (subsample * subsample)));
		//printf("dval val : %f\n", dval);

		do {
			subsample += QSUBSAMPLE_INC;
		} while (width / subsample == x_samples &&
				height / subsample == y_samples);
	}

	dval = sqrt(dval);

	/*
	 double histo_val;
	 histo_val = histo_quality(colour);
	 histo_val = histo_val / (50 * 255 * 255);
	 dval = dval * histo_val;
	 */
	free(buf);
	return dval;
}

/*
 * Subsample a region starting at *ptr of size X size pixels.
 */
static float SubSample(float *ptr, int img_wid, int x_size, int y_size) {
	int x, y;
	float val = 0;
	for (y = 0; y < y_size; ++y) {
		for (x = 0; x < x_size; x++) {
			val += ptr[x];
		}
		ptr += img_wid;
	}
	return val / (float)(x_size * y_size);
}

static double Gradient(float *buf, int width, int height) {
	int pixels;
	int x, y;
	int yborder = (int) ((double) height * QMARGIN) + 1;
	int xborder = (int) ((double) width * QMARGIN) + 1;
	double d1, d2;
	double val, avg = 0;
	float threshold = THRESHOLD_FLOAT;
	unsigned char *map = calloc(width * height, sizeof(unsigned char));

	// pass 1 locate all pixels > threshold and flag the 3x3 region
	// around them for inclusion in the algorithm
	pixels = 0;
	avg = 0;
	for (y = yborder; y < height - yborder; ++y) {
		int o = y * width + xborder;
		for (x = xborder; x < width - xborder; ++x, ++o) {
			if (buf[o] >= threshold) {
				map[o - width - 1] = map[o - width] = map[o - width + 1] = 1;
				map[o - 1] = map[o] = map[o + 1] = 1;
				map[o + width - 1] = map[o + width] = map[o + width + 1] = 1;
				++pixels;
				avg += buf[o];
			}
		}
	}

	// Average of the significant pixels
	if (!pixels) {
		val = -1.0;
		goto end;
	}

	avg /= (double)pixels;

	val = 0.0;
	pixels = 0;

	for (y = yborder; y < height - yborder; ++y) {
		int o = y * width + xborder;      // start of row y
		for (x = xborder; x < width - xborder; ++x, ++o)
			if (map[o]) {
				// Pixel differences
				d1 = buf[o] - buf[o+1];
				d2 = buf[o] - buf[o+width];
				val += (d1 * d1 + d2 * d2);
				pixels++;
			}
	}

	val = val / (double)pixels; // normalise value to per-pixel
	val = val / 10.0;
end:
	free(map);

	return val;
}

/* 3*3 averaging convolution filter, does nothing on the edges */
static float *_smooth_image_float(float *buf, int width, int height) {
	float *new_buff = calloc(width * height * 2, sizeof(float));
	int x, y;

	for (y = 1; y < height - 1; ++y) {
		int o = y * width + 1;
		for (x = 1; x < width - 1; ++x, ++o) {
			float v = 0.0f;
			v += buf[o - width - 1];
			v += buf[o - width];
			v += buf[o - width + 1];

			v += buf[o - 1];
			v += buf[o];
			v += buf[o + 1];

			v += buf[o + width - 1];
			v += buf[o + width];
			v += buf[o + width + 1];

			new_buff[o] = v / 9.0f;
		}
	}

	return new_buff;
}
