/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "algos/quality.h"

static double QualityEstimate_ushort(fits *fit, int layer);
static int32_t SubSample(WORD *ptr, int img_wid, int x_size, int y_size);
static void _smooth_image_16(unsigned short *buf, int width, int height);
static double Gradient(WORD *buf, int width, int height);

double QualityEstimate(fits *fit, int layer) {
	if (fit->type == DATA_USHORT)
		return QualityEstimate_ushort(fit, layer);
	if (fit->type == DATA_FLOAT)
		return QualityEstimate_float(fit, layer);
	return -1.0;
}

// -------------------------------------------------------
// Method to estimate quality.
// Runs on the complete layer and destroys it.
// -------------------------------------------------------
static double QualityEstimate_ushort(fits *fit, int layer) {
	int width = fit->rx;
	int height = fit->ry;
	int x1, y1;
	int subsample, region_w, region_h;
	int i, j, n, x, y, max, x_inc;
	int x_samples, y_samples, y_last;
	WORD *buffer, *buf, maxp[MAXP];
	double q, dval = 0.0;

	buffer = fit->pdata[layer];

	// dimensions of the region we want to analyse
	x1 = 0; y1 = 0;
	region_w = width - 1;
	region_h = height - 1;

	// Allocate the intermediate buffer. Will be 16bpp greyscale
	buf = calloc((region_w / QSUBSAMPLE_MIN + 1) * (region_h / QSUBSAMPLE_MIN + 1), sizeof(WORD));
	if (!buf) {
		PRINT_ALLOC_ERR;
		return -1.0;
	}
	subsample = QSUBSAMPLE_MIN;
	while (subsample <= QSUBSAMPLE_MAX) {
		WORD* ptr;

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
				WORD v = SubSample(ptr, width, subsample, subsample);

				if (v > maxp[2] && v < 65530) {
					int slot;
					if (v > maxp[0]) {
						slot = 0;
					} else if (v > maxp[1]) {
						slot = 1;
					} else {
						slot = 2;
					}

					for (j = MAXP - 1; j > slot; --j)
						maxp[j] = maxp[j - 1];
					maxp[j] = v;
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

		// Stretch histogram
		if (max > 0) {
			double mult = 60000.0 / (double)max;
			for (i = 0; i < n; ++i) {
				unsigned int v = buf[i];
				v = (unsigned int)((double)v * mult);
				if (v > 65535)
					v = 65535;
				buf[i] = v;
			}
		}

		// 3x3 smoothing
		_smooth_image_16(buf, x_samples, y_samples);

#ifdef DEBUG
		/*******************************/
		char filename[1024];
		FILE *out;
		sprintf(filename, "sample_%d.ppm", subsample);
		out = g_fopen(filename, "wb");
		if (out == NULL)
			printf("Cannot write subsampled image %d\n", subsample);
		else {
			fprintf(out, "P5\n%d %d\n255\n", x_samples, y_samples);
			for (i = 0; i < n; ++i)
				putc(buf[i] >> 8, out);
			fclose(out);
		}
		/*********************************/
#endif
		q = Gradient(buf, x_samples, y_samples);

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
static int32_t SubSample(WORD *ptr, int img_wid, int x_size, int y_size) {
	int x, y, val = 0;
	for (y = 0; y < y_size; ++y) {
		for (x = 0; x < x_size; x++) {
			val += ptr[x];
		}
		ptr += img_wid;
	}
	return round_to_WORD((double)val / (double)(x_size * y_size));
}

static double Gradient(WORD *buf, int width, int height) {
	int pixels;
	int x, y;
	int yborder = (int) ((double) height * QMARGIN) + 1;
	int xborder = (int) ((double) width * QMARGIN) + 1;
	double d1, d2;
	double val;
	int threshold = THRESHOLD_USHRT;
	unsigned char *map = calloc(width * height, sizeof(unsigned char));
	if (!map) {
		PRINT_ALLOC_ERR;
		return -1.0;
	}

	// pass 1 locate all pixels > threshold and flag the 3x3 region
	// around them for inclusion in the algorithm
	pixels = 0;
	for (y = yborder; y < height - yborder; ++y) {
		int o = y * width + xborder;
		for (x = xborder; x < width - xborder; ++x, ++o) {
			if (buf[o] >= threshold) {
				map[o - width - 1] = map[o - width] = map[o - width + 1] = 1;
				map[o - 1] = map[o] = map[o + 1] = 1;
				map[o + width - 1] = map[o + width] = map[o + width + 1] = 1;
				++pixels;
			}
		}
	}

	// Average of the significant pixels
	if (!pixels) {
		val = -1.0;
		goto end;
	}

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
static void _smooth_image_16(unsigned short *buf, int width,
		int height) {
	unsigned short lineBuffer[2][width];
	// copy first line to lineBuffer
	for (int x = 0; x < width; ++x) {
		lineBuffer[0][x] = buf[x];
	}
	int prevLine = 0;
	int currLine = 1;
	for (int y = 1; y < height - 1; ++y) {
		int o = y * width;
		// copy current line to lineBuffer
		for (int x = 0; x < width; ++x) {
			lineBuffer[currLine][x] = buf[o + x];
		}
		++o; // increment because we start at x = 1
		for (int x = 1; x < width - 1; ++x, ++o) {
			const unsigned int v = (lineBuffer[prevLine][x - 1] + lineBuffer[prevLine][x])
				+ (lineBuffer[prevLine][x + 1] + lineBuffer[currLine][x - 1]) + (lineBuffer[currLine][x] + lineBuffer[currLine][x + 1])
				+ (buf[o + width - 1] + buf[o + width])
				+ buf[o + width + 1];
			buf[o] = v / 9;
		}
		// swap lineBuffers
		prevLine ^= 1;
		currLine ^= 1;
	}
}

// Scan the region given by (x1,y1) - (x2,y2) and return the barycentre (centre of brightness)

// For a pixel to be counted it's orthogonal neighbors must all be above the threshhold. This
// stops hot pixels and isolated pixels from counting.
// Also, a horizontal gap of 3 or more pixels will cause the last counted pixel to be un-counted.
//

int BlankImageCount = 0;
int MinPixels = 50;

static int _FindCentre_Barycentre_ushort(fits *fit, int x1, int y1, int x2, int y2,
		float *x_avg, float *y_avg) {
	int img_width = fit->rx;
	int img_height = fit->ry;
	int x, y;
	int count = 0;	// count of significant pixels
	float x_total = 0.f, y_total = 0.F;
	float RealThreshHold;

	// must prevent scanning near the edge due to the extended tests below that look
	// +/- 1 pixel above and below
	if (x1 < 1)
		x1 = 1;
	if (y1 < 1)
		y1 = 1;
	if (x2 >= img_width - 1)
		x2 = img_width - 2;
	if (y2 >= img_height - 1)
		y2 = img_height - 2;

	if (get_normalized_value(fit) == UCHAR_MAX_DOUBLE)
		RealThreshHold = THRESHOLD_UCHAR;
	else	RealThreshHold = THRESHOLD_USHRT;

	for (y = y1; y <= y2; ++y) {
		unsigned short *iptr = fit->data + y * img_width + x1;
		for (x = x1; x <= x2; ++x, ++iptr) {
			if (*iptr >= RealThreshHold && *(iptr - 1) >= RealThreshHold
					&& *(iptr + 1) >= RealThreshHold
					&& *(iptr - img_width) >= RealThreshHold
					&& *(iptr + img_width) >= RealThreshHold) {
				x_total += x;
				y_total += y;
				count++;
			}
		}
	}

	if (count == 0) {
		printf("[no image] ");
		if (BlankImageCount >= 0)
			++BlankImageCount;
		return 1;
	}

	if (count < MinPixels) {
		printf("[Not enough pixels. Found %d, require %d] ", count, MinPixels);
		if (BlankImageCount >= 0)
			++BlankImageCount;
		return 1;
	}

	if (count > 0) {
		*x_avg = (x_total / (float) count + 0.5f);
		*y_avg = (y_total / (float) count + 0.5f);
		BlankImageCount = 0;
	}

	*y_avg = img_height - *y_avg;

	return 0;
}

static int _FindCentre_Barycentre_float(fits *fit, int x1, int y1, int x2, int y2,
		float *x_avg, float *y_avg) {
	int img_width = fit->rx;
	int img_height = fit->ry;
	int x, y;
	int count = 0;	// count of significant pixels
	float x_total = 0.f, y_total = 0.f;
	float RealThreshHold;

	// must prevent scanning near the edge due to the extended tests below that look
	// +/- 1 pixel above and below
	if (x1 < 1)
		x1 = 1;
	if (y1 < 1)
		y1 = 1;
	if (x2 >= img_width - 1)
		x2 = img_width - 2;
	if (y2 >= img_height - 1)
		y2 = img_height - 2;

	RealThreshHold = THRESHOLD_FLOAT;

	for (y = y1; y <= y2; ++y) {
		float *iptr = fit->fdata + y * img_width + x1;
		for (x = x1; x <= x2; ++x, ++iptr) {
			if (*iptr >= RealThreshHold && *(iptr - 1) >= RealThreshHold
					&& *(iptr + 1) >= RealThreshHold
					&& *(iptr - img_width) >= RealThreshHold
					&& *(iptr + img_width) >= RealThreshHold) {
				x_total += x;
				y_total += y;
				count++;
			}
		}
	}

	if (count == 0) {
		printf("[no image] ");
		if (BlankImageCount >= 0)
			++BlankImageCount;
		return 1;
	}

	if (count < MinPixels) {
		printf("[Not enough pixels. Found %d, require %d] ", count, MinPixels);
		if (BlankImageCount >= 0)
			++BlankImageCount;
		return 1;
	}

	if (count > 0) {
		*x_avg = (x_total / (float) count + 0.5f);
		*y_avg = (y_total / (float) count + 0.5f);
		BlankImageCount = 0;
	}

	*y_avg = img_height - *y_avg;

	return 0;
}

// find the centre of brightness of the whole image
int FindCentre(fits *fit, float *x_avg, float *y_avg) {
	int x1 = 2;
	int x2 = fit->rx - 3;
	int y1 = 0;
	int y2 = fit->ry - 1;

	if (fit->type == DATA_USHORT) {
		return _FindCentre_Barycentre_ushort(fit, x1, y1, x2, y2, x_avg, y_avg);
	} else {
		return _FindCentre_Barycentre_float(fit, x1, y1, x2, y2, x_avg, y_avg);
	}
}
