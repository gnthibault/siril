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
 *
 * This file is a copy of demosaicing.c with WORD replaced by float, and
 * related adjustments.
 */

#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/command.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "gui/dialogs.h"
#include "io/sequence.h"
#include "algos/demosaicing.h"
#include "algos/statistics.h"


/** Calculate the bayer pattern color from the row and column **/
static inline int FC(const size_t row, const size_t col, const uint32_t filters) {
	return filters >> (((row << 1 & 14) + (col & 1)) << 1) & 3;
}

/* width and height are sizes of the original image */
static int super_pixel(const float *buf, float *newbuf, int width, int height,
		sensor_pattern pattern) {
	int i, col, row;
	float tmp;

	i = 0;
	for (row = 0; row < height - 1; row += 2) {
		for (col = 0; col < width - 1; col += 2) {
			switch (pattern) {
			default:
			case BAYER_FILTER_RGGB:
				newbuf[i + 0] = buf[col + row * width];
				tmp = buf[1 + col + row * width];
				tmp += buf[col + (1 + row) * width];
				tmp *= 0.5f;
				newbuf[i + 1] = tmp;
				newbuf[i + 2] = buf[1 + col + (1 + row) * width];
				break;
			case BAYER_FILTER_BGGR:
				newbuf[i + 2] = buf[col + row * width];
				tmp = buf[1 + col + row * width];
				tmp += buf[(col + row * width) + width];
				tmp *= 0.5f;
				newbuf[i + 1] = tmp;
				newbuf[i + 0] = buf[(1 + col + row * width) + width];
				break;
			case BAYER_FILTER_GBRG:
				newbuf[i + 2] = buf[1 + col + row * width];
				newbuf[i + 0] = buf[(col + row * width) + width];
				tmp = buf[col + row * width];
				tmp += buf[(1 + col + row * width) + width];
				tmp *= 0.5f;
				newbuf[i + 1] = tmp;
				break;
			case BAYER_FILTER_GRBG:
				newbuf[i + 0] = buf[1 + col + row * width];
				newbuf[i + 2] = buf[(col + row * width) + width];
				tmp = buf[col + row * width];
				tmp += buf[(1 + col + row * width) + width];
				tmp *= 0.5f;
				newbuf[i + 1] = tmp;
				break;
			}
			i += 3;
		}
	}
	return 0;
}

/***************************************************
 * 
 * Written by Damien Douxchamps and Frederic Devernay
 * The original VNG and AHD Bayer decoding are from Dave Coffin's DCRAW.
 * https://code.google.com/p/gst-plugins-elphel/
 * 
 * *************************************************/

static void ClearBorders(float *rgb, int sx, int sy, int w) {
	int i, j;

	/* black edges: */
	i = 3 * sx * w - 1;
	j = 3 * sx * sy - 1;
	while (i >= 0) {
		rgb[i--] = 0.0f;
		rgb[j--] = 0.0f;
	}

	int low = sx * (w - 1) * 3 - 1 + w * 3;
	i = low + sx * (sy - w * 2 + 1) * 3;
	while (i > low) {
		j = 6 * w;
		while (j > 0) {
			rgb[i--] = 0.0F;
			j--;
		}
		i -= (sx - 2 * w) * 3;
	}
}

/* OpenCV's Bayer decoding */
static int bayer_Bilinear(const float *bayer, float *rgb, int sx, int sy,
		sensor_pattern tile) {
	const int bayerStep = sx;
	const int rgbStep = 3 * sx;
	int width = sx;
	int height = sy;
	int blue = tile == BAYER_FILTER_BGGR || tile == BAYER_FILTER_GBRG ? -1 : 1;
	int start_with_green = tile == BAYER_FILTER_GBRG
			|| tile == BAYER_FILTER_GRBG;

	if (tile > BAYER_FILTER_MAX || tile < BAYER_FILTER_MIN)
		return -1;

	ClearBorders(rgb, sx, sy, 1);
	rgb += rgbStep + 3 + 1;
	height -= 2;
	width -= 2;

	for (; height--; bayer += bayerStep, rgb += rgbStep) {
		float t0, t1;
		const float *bayerEnd = bayer + width;

		if (start_with_green) {
			t0 = (bayer[1] + bayer[bayerStep * 2 + 1] + 1) * 0.5f;
			t1 = (bayer[bayerStep] + bayer[bayerStep + 2] + 1) * 0.5f;
			rgb[-blue] = t0;
			rgb[0] = bayer[bayerStep + 1];
			rgb[blue] = t1;
			bayer++;
			rgb += 3;
		}

		if (blue > 0) {
			for (; bayer <= bayerEnd - 2; bayer += 2, rgb += 6) {
				t0 = (bayer[0] + bayer[2] + bayer[bayerStep * 2]
						+ bayer[bayerStep * 2 + 2] + 2) * 0.25f;
				t1 = (bayer[1] + bayer[bayerStep] + bayer[bayerStep + 2]
						+ bayer[bayerStep * 2 + 1] + 2) * 0.25f;
				rgb[-1] = t0;
				rgb[0] = t1;
				rgb[1] = bayer[bayerStep + 1];

				t0 = (bayer[2] + bayer[bayerStep * 2 + 2] + 1) * 0.5f;
				t1 = (bayer[bayerStep + 1] + bayer[bayerStep + 3] + 1) * 0.5f;
				rgb[2] = t0;
				rgb[3] = bayer[bayerStep + 2];
				rgb[4] = t1;
			}
		} else {
			for (; bayer <= bayerEnd - 2; bayer += 2, rgb += 6) {
				t0 = (bayer[0] + bayer[2] + bayer[bayerStep * 2]
						+ bayer[bayerStep * 2 + 2] + 2) * 0.25f;
				t1 = (bayer[1] + bayer[bayerStep] + bayer[bayerStep + 2]
						+ bayer[bayerStep * 2 + 1] + 2) * 0.25f;
				rgb[1] = t0;
				rgb[0] = t1;
				rgb[-1] = bayer[bayerStep + 1];

				t0 = (bayer[2] + bayer[bayerStep * 2 + 2] + 1) * 0.5f;
				t1 = (bayer[bayerStep + 1] + bayer[bayerStep + 3] + 1) * 0.5f;
				rgb[4] = t0;
				rgb[3] = bayer[bayerStep + 2];
				rgb[2] = t1;
			}
		}

		if (bayer < bayerEnd) {
			t0 = (bayer[0] + bayer[2] + bayer[bayerStep * 2]
					+ bayer[bayerStep * 2 + 2] + 2) * 0.25f;
			t1 = (bayer[1] + bayer[bayerStep] + bayer[bayerStep + 2]
					+ bayer[bayerStep * 2 + 1] + 2) * 0.25f;
			rgb[-blue] = t0;
			rgb[0] = t1;
			rgb[blue] = bayer[bayerStep + 1];
			bayer++;
			rgb += 3;
		}

		bayer -= width;
		rgb -= width * 3;

		blue = -blue;
		start_with_green = !start_with_green;
	}

	return 0;
}

/* insprired by OpenCV's Bayer decoding */
static int bayer_NearestNeighbor(const float *bayer, float *rgb, int sx, int sy,
		sensor_pattern tile) {
	const int bayerStep = sx;
	const int rgbStep = 3 * sx;
	int width = sx;
	int height = sy;
	int blue = tile == BAYER_FILTER_BGGR || tile == BAYER_FILTER_GBRG ? -1 : 1;
	int start_with_green = tile == BAYER_FILTER_GBRG
			|| tile == BAYER_FILTER_GRBG;
	int i, iinc, imax;
	if ((tile > BAYER_FILTER_MAX) || (tile < BAYER_FILTER_MIN))
		return 1;
	/* add black border */
	imax = sx * sy * 3;
	for (i = sx * (sy - 1) * 3; i < imax; i++) {
		rgb[i] = 0;
	}
	iinc = (sx - 1) * 3;
	for (i = (sx - 1) * 3; i < imax; i += iinc) {
		rgb[i++] = 0;
		rgb[i++] = 0;
		rgb[i++] = 0;
	}
	rgb += 1;
	height -= 1;
	width -= 1;
	for (; height--; bayer += bayerStep, rgb += rgbStep) {
		const float *bayerEnd = bayer + width;
		if (start_with_green) {
			rgb[-blue] = bayer[1];
			rgb[0] = bayer[bayerStep + 1];
			rgb[blue] = bayer[bayerStep];
			bayer++;
			rgb += 3;
		}
		if (blue > 0) {
			for (; bayer <= bayerEnd - 2; bayer += 2, rgb += 6) {
				rgb[-1] = bayer[0];
				rgb[0] = bayer[1];
				rgb[1] = bayer[bayerStep + 1];
				rgb[2] = bayer[2];
				rgb[3] = bayer[bayerStep + 2];
				rgb[4] = bayer[bayerStep + 1];
			}
		} else {
			for (; bayer <= bayerEnd - 2; bayer += 2, rgb += 6) {
				rgb[1] = bayer[0];
				rgb[0] = bayer[1];
				rgb[-1] = bayer[bayerStep + 1];
				rgb[4] = bayer[2];
				rgb[3] = bayer[bayerStep + 2];
				rgb[2] = bayer[bayerStep + 1];
			}
		}
		if (bayer < bayerEnd) {
			rgb[-blue] = bayer[0];
			rgb[0] = bayer[1];
			rgb[blue] = bayer[bayerStep + 1];
			bayer++;
			rgb += 3;
		}
		bayer -= width;
		rgb -= width * 3;
		blue = -blue;
		start_with_green = !start_with_green;
	}
	return 0;
}

#define ABSOLU(x) (((int)(x) ^ ((int)(x) >> 31)) - ((int)(x) >> 31))

static int bayer_VNG(const float *bayer, float *dst, int sx, int sy,
		sensor_pattern pattern) {
	const signed char bayervng_terms[] = { -2, -2, +0, -1, 0, 0x01, -2, -2,
			+0, +0, 1, 0x01, -2, -1, -1, +0, 0, 0x01, -2, -1, +0, -1, 0, 0x02, -2,
			-1, +0, +0, 0, 0x03, -2, -1, +0, +1, 1, 0x01, -2, +0, +0, -1, 0, 0x06,
			-2, +0, +0, +0, 1, 0x02, -2, +0, +0, +1, 0, 0x03, -2, +1, -1, +0, 0,
			0x04, -2, +1, +0, -1, 1, 0x04, -2, +1, +0, +0, 0, 0x06, -2, +1, +0, +1,
			0, 0x02, -2, +2, +0, +0, 1, 0x04, -2, +2, +0, +1, 0, 0x04, -1, -2, -1,
			+0, 0, 0x80, -1, -2, +0, -1, 0, 0x01, -1, -2, +1, -1, 0, 0x01, -1, -2,
			+1, +0, 1, 0x01, -1, -1, -1, +1, 0, 0x88, -1, -1, +1, -2, 0, 0x40, -1,
			-1, +1, -1, 0, 0x22, -1, -1, +1, +0, 0, 0x33, -1, -1, +1, +1, 1, 0x11,
			-1, +0, -1, +2, 0, 0x08, -1, +0, +0, -1, 0, 0x44, -1, +0, +0, +1, 0,
			0x11, -1, +0, +1, -2, 1, 0x40, -1, +0, +1, -1, 0, 0x66, -1, +0, +1, +0,
			1, 0x22, -1, +0, +1, +1, 0, 0x33, -1, +0, +1, +2, 1, 0x10, -1, +1, +1,
			-1, 1, 0x44, -1, +1, +1, +0, 0, 0x66, -1, +1, +1, +1, 0, 0x22, -1, +1,
			+1, +2, 0, 0x10, -1, +2, +0, +1, 0, 0x04, -1, +2, +1, +0, 1, 0x04, -1,
			+2, +1, +1, 0, 0x04, +0, -2, +0, +0, 1, 0x80, +0, -1, +0, +1, 1, 0x88,
			+0, -1, +1, -2, 0, 0x40, +0, -1, +1, +0, 0, 0x11, +0, -1, +2, -2, 0,
			0x40, +0, -1, +2, -1, 0, 0x20, +0, -1, +2, +0, 0, 0x30, +0, -1, +2, +1,
			1, 0x10, +0, +0, +0, +2, 1, 0x08, +0, +0, +2, -2, 1, 0x40, +0, +0, +2,
			-1, 0, 0x60, +0, +0, +2, +0, 1, 0x20, +0, +0, +2, +1, 0, 0x30, +0, +0,
			+2, +2, 1, 0x10, +0, +1, +1, +0, 0, 0x44, +0, +1, +1, +2, 0, 0x10, +0,
			+1, +2, -1, 1, 0x40, +0, +1, +2, +0, 0, 0x60, +0, +1, +2, +1, 0, 0x20,
			+0, +1, +2, +2, 0, 0x10, +1, -2, +1, +0, 0, 0x80, +1, -1, +1, +1, 0,
			0x88, +1, +0, +1, +2, 0, 0x08, +1, +0, +2, -1, 0, 0x40, +1, +0, +2, +1,
			0, 0x10 }, bayervng_chood[] = { -1, -1, -1, 0, -1, +1, 0, +1, +1, +1,
			+1, 0, +1, -1, 0, -1 };
	const int height = sy, width = sx;
	const signed char *cp;
	float (*brow[5])[3], *pix; /* [FD] */
	int code[8][2][320], *ip, gval[8], gmin, gmax, sum[4];
	int row, col, x, y, x1, x2, y1, y2, t, weight, grads, color, diag;
	int g, diff, thold, num, c;
	unsigned long int filters; /* [FD] */

	/* first, use bilinear bayer decoding */

	bayer_Bilinear(bayer, dst, sx, sy, pattern);

	switch (pattern) {
	case BAYER_FILTER_BGGR:
		filters = 0x16161616;
		break;
	case BAYER_FILTER_GRBG:
		filters = 0x61616161;
		break;
	case BAYER_FILTER_RGGB:
		filters = 0x94949494;
		break;
	case BAYER_FILTER_GBRG:
		filters = 0x49494949;
		break;
	default:
		return -1;
	}

	for (row = 0; row < 8; row++) { /* Precalculate for VNG */
		for (col = 0; col < 2; col++) {
			ip = code[row][col];
			for (cp = bayervng_terms, t = 0; t < 64; t++) {
				y1 = *cp++;
				x1 = *cp++;
				y2 = *cp++;
				x2 = *cp++;
				weight = *cp++;
				grads = *cp++;
				color = FC(row + y1, col + x1, filters);
				if (FC(row+y2,col+x2, filters) != (unsigned int) color)
					continue;
				diag = (FC(row,col+1, filters) == (unsigned int) color
						&& FC(row+1,col, filters) == (unsigned int) color) ? 2 : 1;
				if (abs(y1 - y2) == diag && abs(x1 - x2) == diag)
					continue;
				*ip++ = (y1 * width + x1) * 3 + color; /* [FD] */
				*ip++ = (y2 * width + x2) * 3 + color; /* [FD] */
				*ip++ = weight;
				for (g = 0; g < 8; g++)
					if (grads & 1 << g)
						*ip++ = g;
				*ip++ = -1;
			}
			*ip++ = INT_MAX;
			for (cp = bayervng_chood, g = 0; g < 8; g++) {
				y = *cp++;
				x = *cp++;
				*ip++ = (y * width + x) * 3; /* [FD] */
				color = FC(row, col, filters);
				if (FC(row+y,col+x, filters) != (unsigned int) color
						&& FC(row+y*2,col+x*2, filters) == (unsigned int) color)
					*ip++ = (y * width + x) * 6 + color; /* [FD] */
				else
					*ip++ = 0;
			}
		}
	}
	brow[4] = calloc(width * 3, sizeof **brow);
	//merror (brow[4], "vng_interpolate()");
	for (row = 0; row < 3; row++)
		brow[row] = brow[4] + row * width;
	for (row = 2; row < height - 2; row++) { /* Do VNG interpolation */
		for (col = 2; col < width - 2; col++) {
			pix = dst + (row * width + col) * 3; /* [FD] */
			ip = code[row & 7][col & 1];
			memset(gval, 0, sizeof gval);
			while ((g = ip[0]) != INT_MAX) { /* Calculate gradients */
				diff = ABSOLU(pix[g] - pix[ip[1]]) << ip[2];
				gval[ip[3]] += diff;
				ip += 5;
				if ((g = ip[-1]) == -1)
					continue;
				gval[g] += diff;
				while ((g = *ip++) != -1)
					gval[g] += diff;
			}
			ip++;
			gmin = gmax = gval[0]; /* Choose a threshold */
			for (g = 1; g < 8; g++) {
				if (gmin > gval[g])
					gmin = gval[g];
				if (gmax < gval[g])
					gmax = gval[g];
			}
			if (gmax == 0) {
				memcpy(brow[2][col], pix, 3 * sizeof *dst); /* [FD] */
				continue;
			}
			thold = gmin + (gmax * 0.5f);
			memset(sum, 0, sizeof sum);
			color = FC(row, col, filters);
			for (num = g = 0; g < 8; g++, ip += 2) { /* Average the neighbors */
				if (gval[g] <= thold) {
					for (c = 0; c < 3; c++) /* [FD] */
						if (c == color && ip[1])
							sum[c] += (pix[c] + pix[ip[1]]) * 0.5f;
						else
							sum[c] += pix[ip[0] + c];
					num++;
				}
			}
			for (c = 0; c < 3; c++) { /* [FD] Save to buffer */
				float tmp = pix[color];
				if (c != color)
					tmp += (sum[c] - sum[color]) / num;
				//~ CLIP16(t,brow[2][col][c], 16); /* [FD] */
				brow[2][col][c] = tmp; /* [FD] */
			}
		}
		if (row > 3) /* Write buffer to image */
			memcpy(dst + 3 * ((row - 2) * width + 2), brow[0] + 2,
					(width - 4) * 3 * sizeof *dst); /* [FD] */
		for (g = 0; g < 4; g++)
			brow[(g - 1) & 3] = brow[g];
	}
	memcpy(dst + 3 * ((row - 2) * width + 2), brow[0] + 2,
			(width - 4) * 3 * sizeof *dst);
	memcpy(dst + 3 * ((row - 1) * width + 2), brow[1] + 2,
			(width - 4) * 3 * sizeof *dst);
	free(brow[4]);

	return 0;
}

/* AHD interpolation ported from dcraw to libdc1394 by Samuel Audet */
static gboolean ahd_inited = FALSE; /* WARNING: not multi-processor safe */

#define LIM(x,min,max) MAX(min,MIN(x,max))
#define ULIM(x,y,z) ((y) < (z) ? LIM(x,y,z) : LIM(x,z,y))

static const float xyz_rgb[3][3] = { /* XYZ from RGB */
	{ 0.412453f, 0.357580f, 0.180423f },
	{ 0.212671f, 0.715160f, 0.072169f },
	{ 0.019334f, 0.119193f, 0.950227f } };
static const float d65_white[3] = { 0.950456, 1, 1.088754 };

static void cam_to_cielab(float cam[3], float lab[3]) /* [SA] */
{
	float xyz[3];
	static float xyz_cam[3][4];

	if (cam == NULL) {	// this is the init of this function...
		int i, j;
		for (i = 0; i < 3; i++)
			for (j = 0; j < 3; j++) /* [SA] */
				xyz_cam[i][j] = xyz_rgb[i][j] / d65_white[i]; /* [SA] */
	} else {
		int c;

		xyz[0] = xyz[1] = xyz[2] = 0.5;
		for (c = 0; c < 3; c++)
		{ /* [SA] */
			xyz[0] += xyz_cam[0][c] * cam[c];
			xyz[1] += xyz_cam[1][c] * cam[c];
			xyz[2] += xyz_cam[2][c] * cam[c];
		}
		xyz[0] = cbrtf(xyz[0]); /* [SA] */
		xyz[1] = cbrtf(xyz[1]); /* [SA] */
		xyz[2] = cbrtf(xyz[2]); /* [SA] */
		/* TODO: LAB has int ranges, and the algorithms using it use this range too */
		lab[0] = 116.0f * xyz[1] - 16.0f;
		lab[1] = 500.0f * (xyz[0] - xyz[1]);
		lab[2] = 200.0f * (xyz[1] - xyz[2]);
	}
}

/*
 Adaptive Homogeneity-Directed interpolation is based on
 the work of Keigo Hirakawa, Thomas Parks, and Paul Lee.
 */
#define TS 256 /* Tile Size */

static int bayer_AHD(const float *bayer, float *dst, int sx, int sy,
		sensor_pattern pattern) {
	int i, j, top, left, row, col, tr, tc, fc, c, d, hm[2];
	/* the following has the same type as the image */
	float (*pix)[3], (*rix)[3]; /* [SA] */
	static const int dir[4] = { -1, 1, -TS, TS };
	unsigned ldiff[2][4], abdiff[2][4], leps, abeps;
	float flab[3];
	float (*rgb)[TS][TS][3]; /* [SA] */
	short (*lab)[TS][TS][3];
	char (*homo)[TS][TS], *buffer;
	/* start - new code for libdc1394 */
	uint32_t filters;
	const int height = sy, width = sx;
	int x, y;
	float val;

	if (ahd_inited == FALSE) {
		/* WARNING: this might not be multi-processor safe */
		cam_to_cielab(NULL, NULL);
		ahd_inited = TRUE;
	}

	switch (pattern) {
	case BAYER_FILTER_BGGR:
		filters = 0x16161616;
		break;
	case BAYER_FILTER_GRBG:
		filters = 0x61616161;
		break;
	case BAYER_FILTER_RGGB:
		filters = 0x94949494;
		break;
	case BAYER_FILTER_GBRG:
		filters = 0x49494949;
		break;
	default:
		return -1;
	}

	/* fill-in destination with known exact values */
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			int channel = FC(y, x, filters);
			dst[(y * width + x) * 3 + channel] = bayer[y * width + x];
		}
	}
	/* end - new code for libdc1394 */

	/* start - code from border_interpolate(int border) */
	{
		int border = 3;
		unsigned row, col, y, x, f, c;
		float sum[8];

		for (row = 0; row < (unsigned int) height; row++)
			for (col = 0; col < (unsigned int) width; col++) {
				if (col == (unsigned int) border && row >= (unsigned int) border
						&& row < (unsigned int) height - (unsigned int) border)
					col = width - border;
				memset(sum, 0, sizeof sum);
				for (y = row - 1; y != row + 2; y++)
					for (x = col - 1; x != col + 2; x++)
						if (y < (unsigned int) height
								&& x < (unsigned int) width) {
							f = FC(y, x, filters);
							sum[f] += dst[(y * width + x) * 3 + f]; /* [SA] */
							sum[f + 4] += 1.0f;	// why?
						}
				f = FC(row, col, filters);
				for (c = 0; c < 3; c++)
					if (c != f && sum[c + 4]) /* [SA] */
						dst[(row * width + col) * 3 + c] = sum[c] / sum[c + 4]; /* [SA] */
			}
	}
	/* end - code from border_interpolate(int border) */

	buffer = (char *) malloc(26 * TS * TS); /* 1664 kB */
	/* merror (buffer, "ahd_interpolate()"); */
	rgb = (float (*)[TS][TS][3]) buffer; /* [SA] */
	lab = (short (*)[TS][TS][3]) (buffer + 12 * TS * TS);
	homo = (char (*)[TS][TS]) (buffer + 24 * TS * TS);

	for (top = 0; top < height; top += TS - 6)
		for (left = 0; left < width; left += TS - 6) {
			memset(rgb, 0, 12 * TS * TS);

			/* Interpolate green horizontally and vertically: */
			for (row = ((top < 2) ? 2 : top);
					row < top + TS && row < height - 2; row++) {
				col = left + (FC(row,left, filters) == 1);
				if (col < 2)
					col += 2;
				for (fc = FC(row, col, filters); col < left + TS && col < width - 2;
						col += 2) {
					pix = (float(*)[3]) dst + (row * width + col); /* [SA] */
					val = ((pix[-1][1] + pix[0][fc] + pix[1][1]) * 2.0f
							- pix[-2][fc] - pix[2][fc]) * 0.25f;
					rgb[0][row - top][col - left][1] = ULIM(val, pix[-1][1],
							pix[1][1]);
					val = ((pix[-width][1] + pix[0][fc] + pix[width][1]) * 2.0f
							- pix[-2 * width][fc] - pix[2 * width][fc]) * 0.25f;
					rgb[1][row - top][col - left][1] = ULIM(val, pix[-width][1],
							pix[width][1]);
				}
			}
			/* Interpolate red and blue, and convert to CIELab: */
			for (d = 0; d < 2; d++)
				for (row = top + 1; row < top + TS - 1 && row < height - 1;
						row++)
					for (col = left + 1; col < left + TS - 1 && col < width - 1;
							col++) {
						pix = (float (*)[3]) dst + (row * width + col); /* [SA] */
						rix = &rgb[d][row - top][col - left];
						if ((c = 2 - FC(row, col, filters)) == 1) {
							c = FC(row + 1, col, filters);
							val = pix[0][1]
									+ ((pix[-1][2 - c] + pix[1][2 - c]
											- rix[-1][1] - rix[1][1]) * 0.5f);
							rix[0][2 - c] = val; /* [SA] */
							val = pix[0][1]
									+ ((pix[-width][c] + pix[width][c]
											- rix[-TS][1] - rix[TS][1]) * 0.5f);
						} else
							val = rix[0][1]
									+ ((pix[-width - 1][c] + pix[-width + 1][c]
											+ pix[+width - 1][c]
											+ pix[+width + 1][c]
											- rix[-TS - 1][1] - rix[-TS + 1][1]
											- rix[+TS - 1][1] - rix[+TS + 1][1]
											+ 1) * 0.25f);
						rix[0][c] = val; /* [SA] */
						c = FC(row, col, filters);
						rix[0][c] = pix[0][c];
						cam_to_cielab(rix[0], flab);
						for (c = 0; c < 3; c++)
							lab[d][row - top][col - left][c] = 64 * flab[c];
					}
			/* Build homogeneity maps from the CIELab images: */
			memset(homo, 0, 2 * TS * TS);
			for (row = top + 2; row < top + TS - 2 && row < height; row++) {
				tr = row - top;
				for (col = left + 2; col < left + TS - 2 && col < width;
						col++) {
					tc = col - left;
					for (d = 0; d < 2; d++)
						for (i = 0; i < 4; i++)
							ldiff[d][i] = ABSOLU(
									lab[d][tr][tc][0]
											- lab[d][tr][tc + dir[i]][0]);
					leps = MIN(MAX(ldiff[0][0],ldiff[0][1]),
							MAX(ldiff[1][2],ldiff[1][3]));
					for (d = 0; d < 2; d++)
						for (i = 0; i < 4; i++)
							if (i * 0.5f == d || ldiff[d][i] <= leps)
								abdiff[d][i] =
										SQR(
												lab[d][tr][tc][1]
														- lab[d][tr][tc + dir[i]][1]) + SQR(lab[d][tr][tc][2] -lab[d][tr][tc+dir[i]][2]);
					abeps = MIN(MAX(abdiff[0][0],abdiff[0][1]),
							MAX(abdiff[1][2],abdiff[1][3]));
					for (d = 0; d < 2; d++)
						for (i = 0; i < 4; i++)
							if (ldiff[d][i] <= leps && abdiff[d][i] <= abeps)
								homo[d][tr][tc]++;
				}
			}
			/* Combine the most homogenous pixels for the final result: */
			for (row = top + 3; row < top + TS - 3 && row < height - 3; row++) {
				tr = row - top;
				for (col = left + 3; col < left + TS - 3 && col < width - 3;
						col++) {
					tc = col - left;
					for (d = 0; d < 2; d++)
						for (hm[d] = 0, i = tr - 1; i <= tr + 1; i++)
							for (j = tc - 1; j <= tc + 1; j++)
								hm[d] += homo[d][i][j];
					if (hm[0] != hm[1])
						for (c = 0; c < 3; c++)
							dst[(row * width + col) * 3 + c] = 
									rgb[hm[1] > hm[0]][tr][tc][c]; /* [SA] */
					else
						for (c = 0; c < 3; c++)
							dst[(row * width + col) * 3 + c] =
									(rgb[0][tr][tc][c] + rgb[1][tr][tc][c])
											* 0.5f; /* [SA] */
				}
			}
		}
	free(buffer);

	return 0;
}

#define fcol(row, col) xtrans[(row) % 6][(col) % 6]

/* Code from RAWTherapee:
 * It is a simple algorithm. Certainly not the best (probably the worst) but it works yet. */
static int fast_xtrans_interpolate(const float *bayer, float *dst, int sx, int sy, int xtrans[6][6]) {
	uint32_t filters = 9;
	const int height = sy, width = sx;
	int row;

	/* start - code from border_interpolate(int border) */
	{
		int border = 1;
		unsigned row, col, y, x, f, c;
	       	float sum[8];

		for (row = 0; row < (unsigned int) height; row++)
			for (col = 0; col < (unsigned int) width; col++) {
				if (col == (unsigned int) border && row >= (unsigned int) border
						&& row < (unsigned int) height - (unsigned int) border)
					col = width - border;
				memset(sum, 0, sizeof sum);
				for (y = row - 1; y != row + 2; y++)
					for (x = col - 1; x != col + 2; x++)
						if (y < (unsigned int) height
								&& x < (unsigned int) width) {
							f = FC(y, x, filters);
							sum[f] += dst[(y * width + x) * 3 + f]; /* [SA] */
							sum[f + 4] += 1.0f;
						}
				f = FC(row, col, filters);
				for (c = 0; c < 3; c++)
					if (c != f && sum[c + 4]) /* [SA] */
						dst[(row * width + col) * 3 + c] = sum[c] / sum[c + 4]; /* [SA] */
			}
	}
	/* end - code from border_interpolate(int border) */
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
	for (row = 1; row < (height - 1); row ++) {
		int col;
		for (col = 1; col < (width - 1); col ++) {
			float sum[3] = { 0.f };

			int v, h;
			for (v = -1; v <= 1; v++) {
				for (h = -1; h <= 1; h++) {
					sum[fcol(row + v, col + h)] += bayer[(col + h) + (row + v) * width];
				}
			}

			switch (fcol(row, col)) {
			case 0: /* Red */
				dst[(row * width + col) * 3 + 0] = bayer[col + row * width];
				dst[(row * width + col) * 3 + 1] = sum[1] * 0.2f;
				dst[(row * width + col) * 3 + 2] = sum[2] * 0.33333333f;
				break;

			case 1: /* Green */
				dst[(row * width + col) * 3 + 0] = sum[0] * 0.5f;
				dst[(row * width + col) * 3 + 1] = bayer[col + row * width];
				dst[(row * width + col) * 3 + 2] = sum[2] * 0.5f;
				break;

			case 2: /* Blue */
				dst[(row * width + col) * 3 + 0] = sum[0] * 0.33333333f;
				dst[(row * width + col) * 3 + 1] = sum[1] * 0.2f;
				dst[(row * width + col) * 3 + 2] = bayer[col + row * width];
				break;
			}

		}
	}
	return 0;
}

#undef fcol

/**
 * debayer a buffer of a given size into a newly allocated and returned buffer,
 * using the given bayer pattern and interpolation
 *
 * @param buf original RAW data
 * @param width width of image
 * @param height height of image
 * @param interpolation type of interpolation used for demosaicing algorithm
 * @param pattern type of pattern used for demosaicing algorithm
 * @param xtrans this array is only used in case of FUJI XTRANS RAWs. Can be NULL.
 * @return a new buffer of demosaiced data
 */
float *debayer_buffer_float(float *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, int xtrans[6][6]) {
	float *newbuf;
	long npixels;
	int retval;
	switch (interpolation) {
		case BAYER_BILINEAR:
		case BAYER_NEARESTNEIGHBOR:
		case BAYER_VNG:
		case BAYER_AHD:
		case XTRANS:
			npixels = (*width) * (*height);
			break;
		case BAYER_SUPER_PIXEL:
			npixels = (*width / 2 + *width % 2) * (*height / 2 + *height % 2);
			break;
	}
	newbuf = calloc(3, npixels * sizeof(float));
	if (newbuf == NULL) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	switch (interpolation) {
	case BAYER_BILINEAR:
		retval = bayer_Bilinear(buf, newbuf, *width, *height, pattern);
		break;
	case BAYER_NEARESTNEIGHBOR:
		retval = bayer_NearestNeighbor(buf, newbuf, *width, *height, pattern);
		break;
	default:
	case BAYER_VNG:
		retval = bayer_VNG(buf, newbuf, *width, *height, pattern);
		break;
	case BAYER_AHD:
		retval = bayer_AHD(buf, newbuf, *width, *height, pattern);
		break;
	case BAYER_SUPER_PIXEL:
		retval = super_pixel(buf, newbuf, *width, *height, pattern);
		*width = *width / 2 + *width % 2;
		*height = *height / 2 + *height % 2;
		break;
	case XTRANS:
		if (xtrans == NULL)
			return NULL;
		retval = fast_xtrans_interpolate(buf, newbuf, *width, *height, xtrans);
		break;
	}
	if (retval) {
		free(newbuf);
		return NULL;
	}
	return newbuf;
}

/* This function retrieve the xtrans matrix from the FITS header
 */
static int retrieveXTRANSPattern(char *bayer, int xtrans[6][6]) {
	int x, y, i = 0;
	int len;

	len = strlen(bayer);

	if (len == 36) {
		for (x = 0; x < 6; x++) {
			for (y = 0; y < 6; y++) {
				switch (bayer[i]) {
				default:
				case 'R':
					xtrans[x][y] = 0;
					break;
				case 'G':
					xtrans[x][y] = 1;
					break;
				case 'B':
					xtrans[x][y] = 2;
					break;
				}
				i++;
			}
		}
	}
	return 0;
}

int debayer_float(fits* fit, interpolation_method interpolation, gboolean stretch_cfa) {
	int i, j;
	int width = fit->rx;
	int height = fit->ry;
	int npixels;
	float *buf = fit->fdata;
	float *newbuf;
	int xtrans[6][6] = { 0 };
	int xbayeroff = 0;
	int ybayeroff = 0;

	if (interpolation == XTRANS)
		retrieveXTRANSPattern(fit->bayer_pattern, xtrans);
	full_stats_invalidation_from_fit(fit);

	if (!com.debayer.use_bayer_header) {
		xbayeroff = com.debayer.xbayeroff;
		ybayeroff = com.debayer.ybayeroff;
	} else {
		xbayeroff = fit->bayer_xoffset;
		ybayeroff = fit->bayer_yoffset;
	}

	if (xbayeroff == 1) {
		buf += width;
		height--;
	}

	if (ybayeroff == 1) {
		buf++;
	}

	newbuf = debayer_buffer_float(buf, &width, &height, interpolation,
			com.debayer.bayer_pattern, xtrans);
	if (newbuf == NULL) {
		return 1;
	}
	npixels = width * height;

	// usual color RGBRGB format to fits RRGGBB format
	fit->fdata = realloc(fit->fdata, 3 * npixels * sizeof(float));
	fit->naxes[0] = width;
	fit->naxes[1] = height;
	fit->naxes[2] = 3;
	fit->naxis = 3;
	fit->rx = width;
	fit->ry = height;
	fit->fpdata[RLAYER] = fit->fdata;
	fit->fpdata[GLAYER] = fit->fdata + npixels;
	fit->fpdata[BLAYER] = fit->fdata + npixels * 2;
	fit->bitpix = fit->orig_bitpix;
	for (i = 0, j = 0; j < npixels; i += 3, j++) {
		float r = newbuf[i + RLAYER];
		float g = newbuf[i + GLAYER];
		float b = newbuf[i + BLAYER];
		/* TODO:
		if (stretch_cfa && fit->maximum_pixel_value) {
			r = (r / (double) fit->maximum_pixel_value);
			g = (g / (double) fit->maximum_pixel_value);
			b = (b / (double) fit->maximum_pixel_value);
		}*/
		fit->fpdata[RLAYER][j] = r;
		fit->fpdata[GLAYER][j] = g;
		fit->fpdata[BLAYER][j] = b;
	}
	free(newbuf);
	return 0;
}

