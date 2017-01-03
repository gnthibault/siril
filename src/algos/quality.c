/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2017 team free-astro (see more in AUTHORS file)
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
#include "algos/quality.h"

static int32_t SubSample(WORD *ptr, int img_wid, int x_size, int y_size);

static unsigned short *_smooth_image_16(unsigned short *buf, int width,
		int height);

static double Gradient(WORD *buf, int width, int height, int qtype);

// -------------------------------------------------------
// Method to estimate quality.
// Runs on the complete layer and destroys it.
// What is qtype?
// -------------------------------------------------------
double QualityEstimate(fits *fit, int layer, int qtype) {
	WORD *buffer;
	int width = fit->rx;
	int height = fit->ry;
	int x1, y1;
	int subsample, region_w, region_h;
	int i, j, n, x, y, max, maxp[MAXP], x_inc;
	int x_samples, y_samples, y_last;
	WORD *buf;
	WORD *new_image = NULL;
	double mult, q, dval = 0.0;

	buffer = fit->pdata[layer];

#if 0
	if (qtype == QUALTYPE_HISTO) {
		double size_d = (double) ((m_y_planet_max - m_y_planet_min)
				* (m_x_planet_max - m_x_planet_min));
		return histo_quality(colour) / size_d;
	}

	if (qtype == QUALTYPE_TOTAL) {
		return average_pixel_quality(colour);
	}
#endif

	/* dimensions of the region we want to analyse
	 * We take all the tmpfit
	 */
	region_w = width - 1;
	region_h = height - 1;

	x1 = 0;
	y1 = 0;

	// Allocate the intermediate buffer. Will be 16bpp greyscale
	buf = calloc(region_w * region_h, sizeof(WORD));

	subsample = QSUBSAMPLE_MIN;
	while (subsample <= QSUBSAMPLE_MAX) {
		WORD* ptr;

		/*
		 * Number of h & v pixels in subimage
		 */
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
				register int v = SubSample(ptr, width, subsample, subsample);

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

		// Stretch histogram
		if (max > 0) {
			mult = (double) 60000 / (double) max;
			for (i = 0; i < n; ++i) {
				unsigned int v = buf[i];
				v = (unsigned int) ((double) v * mult);
				if (v > 65535)
					v = 65535;
				buf[i] = v;
			}
		}

		// 3x3 smoothing
		new_image = _smooth_image_16(buf, x_samples, y_samples);

#ifdef DEBUG
		/*******************************/
		char filename[1024];
		FILE *out;

		sprintf(filename, "sample_%d.ppm", subsample);
		out = fopen(filename, "wb");
		if (out == NULL)
		printf("Cannot write subsampled image %d\n", subsample);
		else {
			fprintf(out, "P5\n%d %d\n255\n", x_samples, y_samples);
			for (i = 0; i < n; ++i)
			putc(new_image[i] >> 8, out);
			fclose(out);
		}
		/*********************************/
#endif
		q = Gradient(new_image, x_samples, y_samples, qtype);

		if (qtype == QUALTYPE_NINOX) {
			dval += q;
		} else {
			dval += (q
					* ((QSUBSAMPLE_MIN * QSUBSAMPLE_MIN)
							/ (subsample * subsample))); // New scheme
		}

		//printf("dval val : %lf\n", dval);

		do {
			subsample += QSUBSAMPLE_INC;
		} while (width / subsample == x_samples
				&& height / subsample == y_samples);
	}

	if (qtype == QUALTYPE_NORMAL || qtype == QUALTYPE_NINOX) {
		dval = sqrt(dval);
	}

	/*	if (qtype == QUALTYPE_NORMAL) {
	 double histo_val;
	 histo_val = histo_quality(colour);
	 histo_val = histo_val / (50 * 255 * 255);
	 dval = dval * histo_val;
	 }
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

	return val / (x_size * y_size);
}

static double Gradient(WORD *buf, int width, int height, int qtype) {
	int pixels;
	int x, y;
	int yborder = (int) ((double) height * QMARGIN) + 1;
	int xborder = (int) ((double) width * QMARGIN) + 1;
	double d1, d2;
	double val, avg = 0;
	int threshhold = (THRESHOLD) << 8;
	unsigned char *map = calloc(width * height, sizeof(unsigned char));

	// pass 1 locate all pixels > threshhold and flag the 3x3 region
	// around them for inclusion in the algorithm
	pixels = 0;
	avg = 0;
	for (y = yborder; y < height - yborder; ++y) {
		int o = y * width + xborder;
		for (x = xborder; x < width - xborder; ++x, ++o) {
			if (buf[o] >= threshhold) {
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

	avg /= (double) pixels;

	val = 0;
	pixels = 0;

	if (qtype == QUALTYPE_NINOX) {
		for (y = yborder; y < height - yborder; ++y) {
			int o = y * width + xborder;      // start of row y
			for (x = xborder; x < width - xborder; ++x, ++o)
				if (map[o]) {
					// Pixel differences
					d1 = buf[o];
					d2 = buf[o];

					d1 = fabs(d1 - (int) (buf[o + 1]));
					d2 = fabs(d2 - (int) (buf[o + width]));

					val += d1 + d2;
					pixels++;
				}
		}
	} else {
		for (y = yborder; y < height - yborder; ++y) {
			int o = y * width + xborder;      // start of row y
			for (x = xborder; x < width - xborder; ++x, ++o)
				if (map[o]) {
					// Pixel differences
					d1 = buf[o];
					d2 = buf[o];

					d1 = d1 - (int) (buf[o + 1]);
					d2 = d2 - (int) (buf[o + width]);

					val += (d1 * d1 + d2 * d2);
					pixels++;
				}
		}
	}

	if (qtype == QUALTYPE_NINOX) {
		val = val / (double) pixels; // normalise value to per-pixel
		val = (val * 50) / avg;    // arbitrary increase +
								   // normalise differences to compensate for dim/bright images
	} else {
		val = val / (double) pixels; // normalise value to per-pixel
		val = val / 10;
	}

	end:

	return val;
}

/* 3*3 averaging convolution filter */
static unsigned short *_smooth_image_16(unsigned short *buf, int width,
		int height) {
	unsigned short *new_buff = calloc(width * height * 2, sizeof(WORD));
	int x, y;

	for (y = 1; y < height - 1; ++y) {
		int o = y * width + 1;
		for (x = 1; x < width - 1; ++x, ++o) {
			unsigned int v = (unsigned int) buf[o];
			v += (unsigned int) buf[o - width - 1];
			v += (unsigned int) buf[o - width];
			v += (unsigned int) buf[o - width + 1];

			v += (unsigned int) buf[o - 1];
			v += (unsigned int) buf[o + 1];

			v += (unsigned int) buf[o + width - 1];
			v += (unsigned int) buf[o + width];
			v += (unsigned int) buf[o + width + 1];

			new_buff[o] = v / 9;
		}
	}

	return new_buff;
}

// Scan the region given by (x1,y1) - (x2,y2) and return the barycentre (centre of brightness)

// For a pixel to be counted it's orthogonal neighbors must all be above the threshhold. This
// stops hot pixels and isolated pixels from counting.
// Also, a horizontal gap of 3 or more pixels will cause the last counted pixel to be un-counted.
//

int BlankImageCount = 0;
int MinPixels = 50;

static int _FindCentre_Barycentre(fits *fit, int x1, int y1, int x2, int y2,
		double *x_avg, double *y_avg) {
	int img_width = fit->rx;
	int img_height = fit->ry;
	int x, y;
	int count = 0;	// count of significant pixels
	double x_total = 0, y_total = 0;
	double RealThreshHold;

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

	if (get_normalized_value(fit) == UCHAR_MAX)
		RealThreshHold = THRESHOLD;
	else
		RealThreshHold = THRESHOLD * 256;

	for (y = y1; y <= y2; ++y) {
		int rowcount = 0;

		unsigned short *iptr = fit->data + y * img_width + x1;
		for (x = x1; x <= x2; ++x, ++iptr) {
			if (*iptr >= RealThreshHold && *(iptr - 1) >= RealThreshHold
					&& *(iptr + 1) >= RealThreshHold
					&& *(iptr - img_width) >= RealThreshHold
					&& *(iptr + img_width) >= RealThreshHold) {
				x_total += x;
				y_total += y;
				count++;
				rowcount++;
			}
		}
	}

	if (count == 0) {
		printf("[no image] ");
		if (BlankImageCount >= 0)
			++BlankImageCount;
		return (0);
	}

	if (count < MinPixels) {
		printf("[Not enough pixels. Found %d, require %d] ", count, MinPixels);
		if (BlankImageCount >= 0)
			++BlankImageCount;
		return (0);
	}

	if (count > 0) {
		*x_avg = ((double) x_total / (double) count + 0.5);
		*y_avg = ((double) y_total / (double) count + 0.5);
		BlankImageCount = 0;
	}

	*y_avg = img_height - *y_avg;

	return (1);
}

static int _FindCentre(fits *fit, int x1, int y1, int x2, int y2, double *x_avg,
		double *y_avg) {

	return _FindCentre_Barycentre(fit, x1, y1, x2, y2, x_avg, y_avg);

	return 0;
}

// find the centre of brightness of the whole image
int FindCentre(fits *fit, double *x_avg, double *y_avg) {
	int x1 = 2;
	int x2 = fit->rx - 3;
	int y1 = 0;
	int y2 = fit->ry - 1;

	return _FindCentre(fit, x1, y1, x2, y2, x_avg, y_avg);
}

/* Aperture algorithm for determining quality of small objects (round, featureless).
 Uses global value QF_APERTURE_RADIUS as the radius of a circle centered on the target.
 Create a second circle around that so that the annulus between them has the same area (ie
 radius of outer circle = QF_APERTURE_RADIUS * 1.414
 Then calculate difference in adu count between the inner circle and annulus. Higher difference
 values indicate better quality
 Code comes from NINOX
 */

double Aperture(fits *fit, int layer) {
	int width = fit->rx;
	int height = fit->ry;
	WORD *ubuf = (WORD *) fit->pdata[layer];
	int x, y, x1, y1, x2, y2, r1, r2, r;
	double xc, yc;
	double total1, total2, rval;

	xc = yc = 0;

	FindCentre(fit, &xc, &yc);

	r1 = QF_APERTURE_RADIUS;
	r2 = r1 * 1.414 + 0.5;

	x1 = xc - r2;
	if (x1 < 1)
		x1 = 1;
	x2 = xc + r2;
	if (x2 > width - 1)
		x2 = width - 1;
	y1 = yc - r2;
	if (y1 < 1)
		y1 = 1;
	y2 = yc + r2;
	if (y2 > height - 1)
		y2 = height - 1;

	r1 *= r1;
	r2 *= r2;

	total1 = total2 = 0;

	for (y = y1; y <= y2; ++y) {
		int o = y * width + x1;
		int yp = (y - yc) * (y - yc);
		for (x = x1; x <= x2; ++x, ++o) {
			r = yp + (x - xc) * (x - xc);
			if (r <= r2)
				total2 += ubuf[o];
			if (r <= r1) {
				total1 += ubuf[o];
			}
		}
	}

	if (total2 == 0)
		total2 = 1;
	rval = total1 / total2;
	return rval;
}

