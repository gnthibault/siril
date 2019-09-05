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

#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "algos/demosaicing.h"
#include "algos/statistics.h"


/** Calculate the bayer pattern color from the row and column **/
static inline int FC(const size_t row, const size_t col, const uint32_t filters) {
	return filters >> (((row << 1 & 14) + (col & 1)) << 1) & 3;
}

/* width and height are sizes of the original image */
static int super_pixel(const WORD *buf, WORD *newbuf, int width, int height,
		sensor_pattern pattern) {
	int i, col, row;
	double tmp;

	i = 0;
	for (row = 0; row < height - 1; row += 2) {
		for (col = 0; col < width - 1; col += 2) {
			switch (pattern) {
			default:
			case BAYER_FILTER_RGGB:
				newbuf[i + 0] = buf[col + row * width];
				tmp = (double) buf[1 + col + row * width];
				tmp += (double) buf[col + (1 + row) * width];
				tmp *= 0.5;
				newbuf[i + 1] = round_to_WORD(tmp);
				newbuf[i + 2] = buf[1 + col + (1 + row) * width];
				break;
			case BAYER_FILTER_BGGR:
				newbuf[i + 2] = buf[col + row * width];
				tmp = (double) buf[1 + col + row * width];
				tmp += (double) buf[(col + row * width) + width];
				tmp *= 0.5;
				newbuf[i + 1] = round_to_WORD(tmp);
				newbuf[i + 0] = buf[(1 + col + row * width) + width];
				break;
			case BAYER_FILTER_GBRG:
				newbuf[i + 2] = buf[1 + col + row * width];
				newbuf[i + 0] = buf[(col + row * width) + width];
				tmp = (double) buf[col + row * width];
				tmp += (double) buf[(1 + col + row * width) + width];
				tmp *= 0.5;
				newbuf[i + 1] = round_to_WORD(tmp);
				break;
			case BAYER_FILTER_GRBG:
				newbuf[i + 0] = buf[1 + col + row * width];
				newbuf[i + 2] = buf[(col + row * width) + width];
				tmp = (double) buf[col + row * width];
				tmp += (double) buf[(1 + col + row * width) + width];
				tmp *= 0.5;
				newbuf[i + 1] = round_to_WORD(tmp);
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


static void ClearBorders(WORD *rgb, int sx, int sy, int w) {
	int i, j;

	/* black edges: */
	i = 3 * sx * w - 1;
	j = 3 * sx * sy - 1;
	while (i >= 0) {
		rgb[i--] = 0;
		rgb[j--] = 0;
	}

	int low = sx * (w - 1) * 3 - 1 + w * 3;
	i = low + sx * (sy - w * 2 + 1) * 3;
	while (i > low) {
		j = 6 * w;
		while (j > 0) {
			rgb[i--] = 0;
			j--;
		}
		i -= (sx - 2 * w) * 3;
	}
}

/* OpenCV's Bayer decoding */
static int bayer_Bilinear(const WORD *bayer, WORD *rgb, int sx, int sy,
		sensor_pattern tile) {
	const int bayerStep = sx;
	const int rgbStep = 3 * sx;
	int width = sx;
	int height = sy;
	int blue = tile == BAYER_FILTER_BGGR || tile == BAYER_FILTER_GBRG ? -1 : 1;
	int start_with_green = tile == BAYER_FILTER_GBRG
			|| tile == BAYER_FILTER_GRBG;

	if ((tile > BAYER_FILTER_MAX) || (tile < BAYER_FILTER_MIN))
		return -1;

	ClearBorders(rgb, sx, sy, 1);
	rgb += rgbStep + 3 + 1;
	height -= 2;
	width -= 2;

	for (; height--; bayer += bayerStep, rgb += rgbStep) {
		int t0, t1;
		const WORD *bayerEnd = bayer + width;

		if (start_with_green) {
			t0 = (bayer[1] + bayer[bayerStep * 2 + 1] + 1) >> 1;
			t1 = (bayer[bayerStep] + bayer[bayerStep + 2] + 1) >> 1;
			rgb[-blue] = round_to_WORD(t0);
			rgb[0] = bayer[bayerStep + 1];
			rgb[blue] = round_to_WORD(t1);
			bayer++;
			rgb += 3;
		}

		if (blue > 0) {
			for (; bayer <= bayerEnd - 2; bayer += 2, rgb += 6) {
				t0 = (bayer[0] + bayer[2] + bayer[bayerStep * 2]
						+ bayer[bayerStep * 2 + 2] + 2) >> 2;
				t1 = (bayer[1] + bayer[bayerStep] + bayer[bayerStep + 2]
						+ bayer[bayerStep * 2 + 1] + 2) >> 2;
				rgb[-1] = round_to_WORD(t0);
				rgb[0] = round_to_WORD(t1);
				rgb[1] = bayer[bayerStep + 1];

				t0 = (bayer[2] + bayer[bayerStep * 2 + 2] + 1) >> 1;
				t1 = (bayer[bayerStep + 1] + bayer[bayerStep + 3] + 1) >> 1;
				rgb[2] = round_to_WORD(t0);
				rgb[3] = bayer[bayerStep + 2];
				rgb[4] = round_to_WORD(t1);
			}
		} else {
			for (; bayer <= bayerEnd - 2; bayer += 2, rgb += 6) {
				t0 = (bayer[0] + bayer[2] + bayer[bayerStep * 2]
						+ bayer[bayerStep * 2 + 2] + 2) >> 2;
				t1 = (bayer[1] + bayer[bayerStep] + bayer[bayerStep + 2]
						+ bayer[bayerStep * 2 + 1] + 2) >> 2;
				rgb[1] = round_to_WORD(t0);
				rgb[0] = round_to_WORD(t1);
				rgb[-1] = bayer[bayerStep + 1];

				t0 = (bayer[2] + bayer[bayerStep * 2 + 2] + 1) >> 1;
				t1 = (bayer[bayerStep + 1] + bayer[bayerStep + 3] + 1) >> 1;
				rgb[4] = round_to_WORD(t0);
				rgb[3] = bayer[bayerStep + 2];
				rgb[2] = round_to_WORD(t1);
			}
		}

		if (bayer < bayerEnd) {
			t0 = (bayer[0] + bayer[2] + bayer[bayerStep * 2]
					+ bayer[bayerStep * 2 + 2] + 2) >> 2;
			t1 = (bayer[1] + bayer[bayerStep] + bayer[bayerStep + 2]
					+ bayer[bayerStep * 2 + 1] + 2) >> 2;
			rgb[-blue] = round_to_WORD(t0);
			rgb[0] = round_to_WORD(t1);
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
static int bayer_NearestNeighbor(const WORD *bayer, WORD *rgb, int sx, int sy,
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
		const WORD *bayerEnd = bayer + width;
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

static int bayer_VNG(const WORD *bayer, WORD *dst, int sx, int sy,
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
	/* the following has the same type as the image */
	uint16_t (*brow[5])[3], *pix; /* [FD] */
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
			thold = gmin + (gmax >> 1);
			memset(sum, 0, sizeof sum);
			color = FC(row, col, filters);
			for (num = g = 0; g < 8; g++, ip += 2) { /* Average the neighbors */
				if (gval[g] <= thold) {
					for (c = 0; c < 3; c++) /* [FD] */
						if (c == color && ip[1])
							sum[c] += (pix[c] + pix[ip[1]]) >> 1;
						else
							sum[c] += pix[ip[0] + c];
					num++;
				}
			}
			for (c = 0; c < 3; c++) { /* [FD] Save to buffer */
				t = pix[color];
				if (c != color)
					t += (sum[c] - sum[color]) / num;
				//~ CLIP16(t,brow[2][col][c], 16); /* [FD] */
				brow[2][col][c] = round_to_WORD(t); /* [FD] */
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
#define CLIPOUT(x) LIM(x,0,255)
#define CLIPOUT16(x,bits) LIM(x,0,((1<<bits)-1))

static const double xyz_rgb[3][3] = { /* XYZ from RGB */
{ 0.412453, 0.357580, 0.180423 }, { 0.212671, 0.715160, 0.072169 }, { 0.019334,
		0.119193, 0.950227 } };
static const float d65_white[3] = { 0.950456, 1, 1.088754 };

static void cam_to_cielab(uint16_t cam[3], float lab[3]) /* [SA] */
{
	float xyz[3];
	static float cbrt[0x10000], xyz_cam[3][4];

	if (cam == NULL) {
		int i, j;

		for (i = 0; i < 0x10000; i++) {
			float r = i / 65535.0;
			cbrt[i] = r > 0.008856 ? pow(r, 1 / 3.0) : 7.787 * r + 16 / 116.0;
		}
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
		xyz[0] = cbrt[round_to_WORD(xyz[0])]; /* [SA] */
		xyz[1] = cbrt[round_to_WORD(xyz[1])]; /* [SA] */
		xyz[2] = cbrt[round_to_WORD(xyz[2])]; /* [SA] */
		lab[0] = 116 * xyz[1] - 16;
		lab[1] = 500 * (xyz[0] - xyz[1]);
		lab[2] = 200 * (xyz[1] - xyz[2]);
	}
}

/*
 Adaptive Homogeneity-Directed interpolation is based on
 the work of Keigo Hirakawa, Thomas Parks, and Paul Lee.
 */
#define TS 256 /* Tile Size */

static int bayer_AHD(const WORD *bayer, WORD *dst, int sx, int sy,
		sensor_pattern pattern) {
	int i, j, top, left, row, col, tr, tc, fc, c, d, val, hm[2];
	/* the following has the same type as the image */
	uint16_t (*pix)[3], (*rix)[3]; /* [SA] */
	static const int dir[4] = { -1, 1, -TS, TS };
	unsigned ldiff[2][4], abdiff[2][4], leps, abeps;
	float flab[3];
	uint16_t (*rgb)[TS][TS][3]; /* [SA] */
	short (*lab)[TS][TS][3];
	char (*homo)[TS][TS], *buffer;
	/* start - new code for libdc1394 */
	uint32_t filters;
	const int height = sy, width = sx;
	int x, y;

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
		unsigned row, col, y, x, f, c, sum[8];

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
							sum[f + 4]++;
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
	rgb = (uint16_t (*)[TS][TS][3]) buffer; /* [SA] */
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
					pix = (uint16_t (*)[3]) dst + (row * width + col); /* [SA] */
					val = ((pix[-1][1] + pix[0][fc] + pix[1][1]) * 2
							- pix[-2][fc] - pix[2][fc]) >> 2;
					rgb[0][row - top][col - left][1] = ULIM(val, pix[-1][1],
							pix[1][1]);
					val = ((pix[-width][1] + pix[0][fc] + pix[width][1]) * 2
							- pix[-2 * width][fc] - pix[2 * width][fc]) >> 2;
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
						pix = (uint16_t (*)[3]) dst + (row * width + col); /* [SA] */
						rix = &rgb[d][row - top][col - left];
						if ((c = 2 - FC(row, col, filters)) == 1) {
							c = FC(row + 1, col, filters);
							val = pix[0][1]
									+ ((pix[-1][2 - c] + pix[1][2 - c]
											- rix[-1][1] - rix[1][1]) >> 1);
							rix[0][2 - c] = round_to_WORD(val); /* [SA] */
							val = pix[0][1]
									+ ((pix[-width][c] + pix[width][c]
											- rix[-TS][1] - rix[TS][1]) >> 1);
						} else
							val = rix[0][1]
									+ ((pix[-width - 1][c] + pix[-width + 1][c]
											+ pix[+width - 1][c]
											+ pix[+width + 1][c]
											- rix[-TS - 1][1] - rix[-TS + 1][1]
											- rix[+TS - 1][1] - rix[+TS + 1][1]
											+ 1) >> 2);
						rix[0][c] = round_to_WORD((double) val); /* [SA] */
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
							if (i >> 1 == d || ldiff[d][i] <= leps)
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
							dst[(row * width + col) * 3 + c] = round_to_WORD(
									rgb[hm[1] > hm[0]][tr][tc][c]); /* [SA] */
					else
						for (c = 0; c < 3; c++)
							dst[(row * width + col) * 3 + c] = round_to_WORD(
									(rgb[0][tr][tc][c] + rgb[1][tr][tc][c])
											>> 1); /* [SA] */
				}
			}
		}
	free(buffer);

	return 0;
}

#define fcol(row, col) xtrans[(row) % 6][(col) % 6]

/* Code from RAWTherapee:
 * It is a simple algorithm. Certainly not the best (probably the worst) but it works yet. */
static int fast_xtrans_interpolate(const WORD *bayer, WORD *dst, int sx, int sy,
		sensor_pattern pattern, int xtrans[6][6]) {
	uint32_t filters = 9;
	const int height = sy, width = sx;
	int row;

	/* start - code from border_interpolate(int border) */
	{
		int border = 1;
		unsigned row, col, y, x, f, c, sum[8];

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
							sum[f + 4]++;
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
WORD *debayer_buffer(WORD *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, int xtrans[6][6]) {
	WORD *newbuf;
	int npixels;

	switch (interpolation) {
	case BAYER_BILINEAR:
		npixels = (*width) * (*height);
		newbuf = calloc(1, 3 * npixels * sizeof(WORD));
		if (newbuf == NULL) {
			PRINT_ALLOC_ERR;
			return NULL;
		}
		//siril_log_message("Bilinear interpolation...\n");
		bayer_Bilinear(buf, newbuf, *width, *height, pattern);
		break;
	case BAYER_NEARESNEIGHBOR:
		npixels = (*width) * (*height);
		newbuf = calloc(1, 3 * npixels * sizeof(WORD));
		if (newbuf == NULL) {
			PRINT_ALLOC_ERR;
			return NULL;
		}
		//siril_log_message("Nearest Neighbor interpolation...\n");
		bayer_NearestNeighbor(buf, newbuf, *width, *height, pattern);
		break;
	default:
	case BAYER_VNG:
		npixels = (*width) * (*height);
		newbuf = calloc(1, 3 * npixels * sizeof(WORD));
		if (newbuf == NULL) {
			PRINT_ALLOC_ERR;
			return NULL;
		}
		//siril_log_message("VNG interpolation...\n");
		bayer_VNG(buf, newbuf, *width, *height, pattern);
		break;
	case BAYER_AHD:
		npixels = (*width) * (*height);
		newbuf = calloc(1, 3 * npixels * sizeof(WORD));
		if (newbuf == NULL) {
			PRINT_ALLOC_ERR;
			return NULL;
		}
		//siril_log_message("AHD interpolation...\n");
		bayer_AHD(buf, newbuf, *width, *height, pattern);
		break;
	case BAYER_SUPER_PIXEL:
		npixels = (*width / 2 + *width % 2) * (*height / 2 + *height % 2);
		newbuf = calloc(1, 3 * npixels * sizeof(WORD));
		if (newbuf == NULL) {
			PRINT_ALLOC_ERR;
			return NULL;
		}
		//siril_log_message("Super Pixel interpolation...\n");
		super_pixel(buf, newbuf, *width, *height, pattern);
		*width = *width / 2 + *width % 2;
		*height = *height / 2 + *height % 2;
		break;
	case XTRANS:
		if (xtrans == NULL)
			return NULL;
		npixels = (*width) * (*height);
		newbuf = calloc(1, 3 * npixels * sizeof(WORD));
		if (newbuf == NULL) {
			PRINT_ALLOC_ERR;
			return NULL;
		}
		fast_xtrans_interpolate(buf, newbuf, *width, *height, pattern, xtrans);
		break;
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

int debayer(fits* fit, interpolation_method interpolation, gboolean stretch_cfa) {
	int i, j;
	int width = fit->rx;
	int height = fit->ry;
	int npixels;
	WORD *buf = fit->data;
	WORD *newbuf;
	int xtrans[6][6] = { 0 };

	retrieveXTRANSPattern(fit->bayer_pattern, xtrans);
	full_stats_invalidation_from_fit(fit);

	if (fit->bayer_yoffset == 1) {
		buf += width;
		height--;
	}

	if (fit->bayer_xoffset == 1) {
		buf++;
	}

	newbuf = debayer_buffer(buf, &width, &height, interpolation,
			com.debayer.bayer_pattern, xtrans);
	if (newbuf == NULL) {
		return 1;
	}
	npixels = width * height;

	// usual color RGBRGB format to fits RRGGBB format
	fit->data = realloc(fit->data, 3 * npixels * sizeof(WORD));
	fit->naxes[0] = width;
	fit->naxes[1] = height;
	fit->naxes[2] = 3;
	fit->naxis = 3;
	fit->rx = width;
	fit->ry = height;
	fit->pdata[RLAYER] = fit->data;
	fit->pdata[GLAYER] = fit->data + npixels;
	fit->pdata[BLAYER] = fit->data + npixels * 2;
	fit->bitpix = fit->orig_bitpix;
	for (i = 0, j = 0; j < npixels; i += 3, j++) {
		double r = (double) newbuf[i + RLAYER];
		double g = (double) newbuf[i + GLAYER];
		double b = (double) newbuf[i + BLAYER];
		if (stretch_cfa && fit->maximum_pixel_value) {
			double norm = fit->bitpix == 8 ? UCHAR_MAX_DOUBLE : USHRT_MAX_DOUBLE;
			r = (r / (double) fit->maximum_pixel_value) * norm;
			g = (g / (double) fit->maximum_pixel_value) * norm;
			b = (b / (double) fit->maximum_pixel_value) * norm;
		}
		fit->pdata[RLAYER][j] =
				(fit->bitpix == 8) ? round_to_BYTE(r) : round_to_WORD(r);
		fit->pdata[GLAYER][j] =
				(fit->bitpix == 8) ? round_to_BYTE(g) : round_to_WORD(g);
		fit->pdata[BLAYER][j] =
				(fit->bitpix == 8) ? round_to_BYTE(b) : round_to_WORD(b);
	}
	free(newbuf);
	return 0;
}

/* From an area, get the area corresponding to the debayer data for all colors,
 * the dashed area below.
 * 0 1 2 3 4 5
 * - - - - - -
 * - - - - - -
 * - - G R - -
 * - - B G - -
 * - - - - - -
 * - - - - - -
 *
 * area is the requested area of an image (simplified as GRBG above)
 * debayer_area is the result of this function, the area with enough pixels to
 *	have a valid debayer
 * image_area is the size of the image, to avoid going out of bounds
 * debayer_offset_x and y are the offset that need to be applied to the debayer
 *	data to find the original area (between 0 and 3).
 */
void get_debayer_area(const rectangle *area, rectangle *debayer_area,
		const rectangle *image_area, int *debayer_offset_x,
		int *debayer_offset_y) {
	int right, bottom;	// temp debayer negative offsets

	/* left side */
	if (area->x & 1)
		*debayer_offset_x = 3;
	else
		*debayer_offset_x = 2;
	if (area->x - *debayer_offset_x < 0) {
		debayer_area->x = 0;
		*debayer_offset_x = area->x;
	} else {
		debayer_area->x = area->x - *debayer_offset_x;
	}

	/* right side */
	int xend = area->x + area->w - 1;
	if (xend & 1)
		right = 2;
	else
		right = 3;
	if (xend + right >= image_area->w) {
		right = image_area->w - xend - 1;
	}
	debayer_area->w = area->w + (area->x - debayer_area->x) + right;

	/* top */
	if (area->y & 1)
		*debayer_offset_y = 3;
	else
		*debayer_offset_y = 2;
	if (area->y - *debayer_offset_y < 0) {
		debayer_area->y = 0;
		*debayer_offset_y = area->y;
	} else {
		debayer_area->y = area->y - *debayer_offset_y;
	}

	/* bottom */
	int yend = area->y + area->h - 1;
	if (yend & 1)
		bottom = 2;
	else
		bottom = 3;
	if (yend + bottom >= image_area->h) {
		bottom = image_area->h - yend - 1;
	}
	debayer_area->h = area->h + (area->y - debayer_area->y) + bottom;

	assert(debayer_area->x < image_area->w);
	assert(debayer_area->y < image_area->h);
	assert(debayer_area->h > 2);
	assert(debayer_area->w > 2);
}

int split_cfa(fits *in, fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3) {
	int width = in->rx;
	int height = in->ry;
	int j, row, col;

	if (strlen(in->bayer_pattern) > 4) {
		siril_log_message(_("Split CFA does not work on non-Bayer filter camera images!\n"));
		return 1;
	}

	width = width / 2 + width % 2;
	height = height / 2 + height % 2;

	if (new_fit_image(&cfa0, width, height, 1)) {
		return 1;
	}

	if (new_fit_image(&cfa1, width, height, 1)) {
		return 1;
	}

	if (new_fit_image(&cfa2, width, height, 1)) {
		return 1;
	}

	if (new_fit_image(&cfa3, width, height, 1)) {
		return 1;
	}

	double c0, c1, c2, c3;
	j = 0;

	for (row = 0; row < in->ry - 1; row += 2) {
		for (col = 0; col < in->rx - 1; col += 2) {
			/* not c0, c1, c2 and c3 because of the read orientation */
			c1 = in->data[col + row * in->rx];
			c3 = in->data[1 + col + row * in->rx];
			c0 = in->data[col + (1 + row) * in->rx];
			c2 = in->data[1 + col + (1 + row) * in->rx];

			cfa0->data[j] =
					(in->bitpix == 8) ? round_to_BYTE(c0) : round_to_WORD(c0);
			cfa1->data[j] =
					(in->bitpix == 8) ? round_to_BYTE(c1) : round_to_WORD(c1);
			cfa2->data[j] =
					(in->bitpix == 8) ? round_to_BYTE(c2) : round_to_WORD(c2);
			cfa3->data[j] =
					(in->bitpix == 8) ? round_to_BYTE(c3) : round_to_WORD(c3);
			j++;
		}
	}

	return 0;
}

int split_cfa_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_) {
	fits f_cfa0 = { 0 };
	fits f_cfa1 = { 0 };
	fits f_cfa2 = { 0 };
	fits f_cfa3 = { 0 };

	struct split_cfa_data *cfa_args = (struct split_cfa_data *) args->user;

	gchar *cfa0 = g_strdup_printf("%s0_%s_%05d%s", cfa_args->seqEntry, cfa_args->seq->seqname, o, com.ext);
	gchar *cfa1 = g_strdup_printf("%s1_%s_%05d%s", cfa_args->seqEntry, cfa_args->seq->seqname, o, com.ext);
	gchar *cfa2 = g_strdup_printf("%s2_%s_%05d%s", cfa_args->seqEntry, cfa_args->seq->seqname, o, com.ext);
	gchar *cfa3 = g_strdup_printf("%s3_%s_%05d%s", cfa_args->seqEntry, cfa_args->seq->seqname, o, com.ext);

	int ret = split_cfa(fit, &f_cfa0, &f_cfa1, &f_cfa2, &f_cfa3);
	if (ret) {
		g_free(cfa0);
		g_free(cfa1);
		g_free(cfa2);
		g_free(cfa3);
		clearfits(&f_cfa0);
		clearfits(&f_cfa1);
		clearfits(&f_cfa2);
		clearfits(&f_cfa3);
		return ret;
	}

	save1fits16(cfa0, &f_cfa0, 0);
	save1fits16(cfa1, &f_cfa1, 0);
	save1fits16(cfa2, &f_cfa2, 0);
	save1fits16(cfa3, &f_cfa3, 0);

	g_free(cfa0);
	g_free(cfa1);
	g_free(cfa2);
	g_free(cfa3);

	clearfits(&f_cfa0);
	clearfits(&f_cfa1);
	clearfits(&f_cfa2);
	clearfits(&f_cfa3);

	return 0;
}

void apply_split_cfa_to_sequence(struct split_cfa_data *split_cfa_args) {
	struct generic_seq_args *args = malloc(sizeof(struct generic_seq_args));
	args->seq = split_cfa_args->seq;
	args->partial_image = FALSE;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = split_cfa_args->seq->selnum;
	args->prepare_hook = ser_prepare_hook;
	args->finalize_hook = ser_finalize_hook;
	args->save_hook = NULL;
	args->image_hook = split_cfa_image_hook;
	args->idle_function = NULL;
	args->stop_on_error = TRUE;
	args->description = _("Split CFA");
	args->has_output = FALSE;
	args->new_seq_prefix = split_cfa_args->seqEntry;
	args->load_new_sequence = FALSE;
	args->force_ser_output = FALSE;
	args->user = split_cfa_args;
	args->already_in_a_thread = FALSE;
	args->parallel = TRUE;

	split_cfa_args->fit = NULL;	// not used here

	start_in_new_thread(generic_sequence_worker, args);
}
