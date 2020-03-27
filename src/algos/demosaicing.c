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

#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/command.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"
#include "gui/dialogs.h"
#include "io/sequence.h"
#include "algos/demosaicing.h"
#include "algos/statistics.h"

#define USE_SIRIL_DEBAYER FALSE

/** Calculate the bayer pattern color from the row and column **/
static inline int FC(const size_t row, const size_t col, const uint32_t filters) {
	return filters >> (((row << 1 & 14) + (col & 1)) << 1) & 3;
}

/* width and height are sizes of the original image */
static void super_pixel_ushort(const WORD *buf, WORD *newbuf, int width, int height,
		sensor_pattern pattern) {
	long i = 0;
	for (int row = 0; row < height - 1; row += 2) {
		for (int col = 0; col < width - 1; col += 2) {
			float tmp;
			switch (pattern) {
			default:
			case BAYER_FILTER_RGGB:
				newbuf[i + 0] = buf[col + row * width];
				tmp = buf[1 + col + row * width];
				tmp += buf[col + (1 + row) * width];
				newbuf[i + 1] = round_to_WORD(tmp * 0.5f);
				newbuf[i + 2] = buf[1 + col + (1 + row) * width];
				break;
			case BAYER_FILTER_BGGR:
				newbuf[i + 2] = buf[col + row * width];
				tmp = buf[1 + col + row * width];
				tmp += buf[(col + row * width) + width];
				newbuf[i + 1] = round_to_WORD(tmp * 0.5f);
				newbuf[i + 0] = buf[(1 + col + row * width) + width];
				break;
			case BAYER_FILTER_GBRG:
				newbuf[i + 2] = buf[1 + col + row * width];
				newbuf[i + 0] = buf[(col + row * width) + width];
				tmp = buf[col + row * width];
				tmp += buf[(1 + col + row * width) + width];
				newbuf[i + 1] = round_to_WORD(tmp * 0.5f);
				break;
			case BAYER_FILTER_GRBG:
				newbuf[i + 0] = buf[1 + col + row * width];
				newbuf[i + 2] = buf[(col + row * width) + width];
				tmp = buf[col + row * width];
				tmp += buf[(1 + col + row * width) + width];
				newbuf[i + 1] = round_to_WORD(tmp * 0.5f);
				break;
			}
			i += 3;
		}
	}
}

/* width and height are sizes of the original image */
static void super_pixel_float(const float *buf, float *newbuf, int width, int height,
		sensor_pattern pattern) {
	long i = 0;
	for (int row = 0; row < height - 1; row += 2) {
		for (int col = 0; col < width - 1; col += 2) {
			float tmp;
			switch (pattern) {
			default:
			case BAYER_FILTER_RGGB:
				newbuf[i + 0] = buf[col + row * width];
				tmp = buf[1 + col + row * width];
				tmp += buf[col + (1 + row) * width];
				newbuf[i + 1] = tmp * 0.5f;
				newbuf[i + 2] = buf[1 + col + (1 + row) * width];
				break;
			case BAYER_FILTER_BGGR:
				newbuf[i + 2] = buf[col + row * width];
				tmp = buf[1 + col + row * width];
				tmp += buf[(col + row * width) + width];
				newbuf[i + 1] = tmp * 0.5f;
				newbuf[i + 0] = buf[(1 + col + row * width) + width];
				break;
			case BAYER_FILTER_GBRG:
				newbuf[i + 2] = buf[1 + col + row * width];
				newbuf[i + 0] = buf[(col + row * width) + width];
				tmp = buf[col + row * width];
				tmp += buf[(1 + col + row * width) + width];
				newbuf[i + 1] = tmp * 0.5f;
				break;
			case BAYER_FILTER_GRBG:
				newbuf[i + 0] = buf[1 + col + row * width];
				newbuf[i + 2] = buf[(col + row * width) + width];
				tmp = buf[col + row * width];
				tmp += buf[(1 + col + row * width) + width];
				newbuf[i + 1] = tmp * 0.5f;
				break;
			}
			i += 3;
		}
	}
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

	if (tile > BAYER_FILTER_MAX || tile < BAYER_FILTER_MIN)
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
	WORD (*brow[5])[3], *pix; /* [FD] */
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

static const float xyz_rgb[3][3] = { /* XYZ from RGB */
	{ 0.412453f, 0.357580f, 0.180423f },
	{ 0.212671f, 0.715160f, 0.072169f },
	{ 0.019334f, 0.119193f, 0.950227f } };
/* TODO: is it wise to use a D65 here? */
static const float d65_white[3] = { 0.950456f, 1.0f, 1.088754f };
/* TODO: store the precomputation of xyz_rgb * d65 instead of running a silly init */

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
		unsigned int xtrans[6][6]) {
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

static WORD *debayer_buffer_siril(WORD *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern, unsigned int xtrans[6][6]) {
	WORD *newbuf;
	size_t npixels;
	int retval;
	switch (interpolation) {
		default:
			npixels = (*width) * (*height);
			break;
		case BAYER_SUPER_PIXEL:
			npixels = (*width / 2 + *width % 2) * (*height / 2 + *height % 2);
			break;
	}
	newbuf = calloc(3, npixels * sizeof(WORD));
	if (!newbuf) {
		PRINT_ALLOC_ERR;
		return NULL;
	}

	switch (interpolation) {
	case BAYER_BILINEAR:
		retval = bayer_Bilinear(buf, newbuf, *width, *height, pattern);
		break;
//	case BAYER_NEARESTNEIGHBOR:
//		retval = bayer_NearestNeighbor(buf, newbuf, *width, *height, pattern);
//		break;
	default:
	case BAYER_VNG:
		retval = bayer_VNG(buf, newbuf, *width, *height, pattern);
		break;
	case BAYER_AHD:
		retval = bayer_AHD(buf, newbuf, *width, *height, pattern);
		break;
	case BAYER_SUPER_PIXEL:
		super_pixel_ushort(buf, newbuf, *width, *height, pattern);
		*width = *width / 2 + *width % 2;
		*height = *height / 2 + *height % 2;
		retval = 0;
		break;
	case XTRANS:
		if (!xtrans)
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

WORD *debayer_buffer_superpixel_ushort(WORD *buf, int *width, int *height, sensor_pattern pattern) {
	int new_rx = *width / 2 + *width % 2;
	int new_ry = *height / 2 + *height % 2;
	size_t npixels = new_rx * new_ry;
	WORD *newbuf = malloc(3 * npixels * sizeof(WORD));
	if (!newbuf) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	super_pixel_ushort(buf, newbuf, *width, *height, pattern);
	*width = new_rx;
	*height = new_ry;
	return newbuf;
}

/**
 * debayer a buffer of a given size into a newly allocated and returned buffer,
 * using the given bayer pattern and interpolation (only used for SER demosaicing)
 *
 * @param buf original RAW data
 * @param width width of image
 * @param height height of image
 * @param interpolation type of interpolation used for demosaicing algorithm
 * @param pattern type of pattern used for demosaicing algorithm
 * @return a new buffer of demosaiced data
 */
WORD *debayer_buffer(WORD *buf, int *width, int *height,
		interpolation_method interpolation, sensor_pattern pattern) {
	if (USE_SIRIL_DEBAYER)
		return debayer_buffer_siril(buf, width, height, interpolation, pattern, NULL);
	return debayer_buffer_new_ushort(buf, width, height, interpolation, pattern, NULL);
}

float *debayer_buffer_superpixel_float(float *buf, int *width, int *height, sensor_pattern pattern) {
	int new_rx = *width / 2 + *width % 2;
	int new_ry = *height / 2 + *height % 2;
	size_t npixels = new_rx * new_ry;
	float *newbuf = malloc(3 * npixels * sizeof(float));
	if (!newbuf) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	super_pixel_float(buf, newbuf, *width, *height, pattern);
	*width = new_rx;
	*height = new_ry;
	return newbuf;
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

/* This function retrieve the xtrans matrix from the FITS header */
int retrieveXTRANSPattern(char *bayer, unsigned int xtrans[6][6]) {
	int x, y, i = 0;

	if (strlen(bayer) != 36) {
		siril_log_color_message(_("FITS header does not contain a proper XTRANS pattern, demosaicing cannot be done"), "red");
		return 1;
	}

	for (x = 0; x < 6; x++) {
		for (y = 0; y < 6; y++) {
			switch (bayer[i]) {
				default:	// shouldn't default be an error?
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
	return 0;
}

static int debayer_ushort(fits *fit, interpolation_method interpolation, sensor_pattern pattern) {
	size_t i, j, npixels = fit->naxes[0] * fit->naxes[1];
	int width = fit->rx;
	int height = fit->ry;
	WORD *buf = fit->data;
	int xbayeroff = 0, ybayeroff = 0;

	unsigned int xtrans[6][6];
	if (interpolation == XTRANS) {
		retrieveXTRANSPattern(fit->bayer_pattern, xtrans);
	}

	if (!com.debayer.use_bayer_header) {
		xbayeroff = com.debayer.xbayeroff;
		ybayeroff = com.debayer.ybayeroff;
	} else {
		xbayeroff = fit->bayer_xoffset;
		ybayeroff = fit->bayer_yoffset;
	}

	if (xbayeroff == 1) {
		switch (pattern) {
		case BAYER_FILTER_RGGB:
			pattern = BAYER_FILTER_GRBG;
			break;
		case BAYER_FILTER_BGGR:
			pattern = BAYER_FILTER_GBRG;
			break;
		case BAYER_FILTER_GBRG:
			pattern = BAYER_FILTER_BGGR;
			break;
		case BAYER_FILTER_GRBG:
			pattern = BAYER_FILTER_RGGB;
			break;
		default:
		case BAYER_FILTER_NONE:
			return 1;
		}
	}

	if (ybayeroff == 1) {
		switch (pattern) {
		case BAYER_FILTER_RGGB:
			pattern = BAYER_FILTER_GBRG;
			break;
		case BAYER_FILTER_BGGR:
			pattern = BAYER_FILTER_GRBG;
			break;
		case BAYER_FILTER_GBRG:
			pattern = BAYER_FILTER_RGGB;
			break;
		case BAYER_FILTER_GRBG:
			pattern = BAYER_FILTER_BGGR;
			break;
		default:
		case BAYER_FILTER_NONE:
			return 1;
		}
	}

	if (USE_SIRIL_DEBAYER) {
		WORD *newbuf = debayer_buffer_siril(buf, &width, &height, interpolation, pattern, xtrans);
		if (!newbuf)
			return 1;

		fit_debayer_buffer(fit, newbuf);
		// size might have changed, in case of superpixel
		fit->naxes[0] = width;
		fit->naxes[1] = height;
		fit->rx = width;
		fit->ry = height;
		fit->bitpix = fit->orig_bitpix;
		// color RGBRGB format to fits RRGGBB format
		for (i = 0, j = 0; j < npixels; i += 3, j++) {
			double r = (double) newbuf[i + RLAYER];
			double g = (double) newbuf[i + GLAYER];
			double b = (double) newbuf[i + BLAYER];
			fit->pdata[RLAYER][j] =
					(fit->bitpix == 8) ? round_to_BYTE(r) : round_to_WORD(r);
			fit->pdata[GLAYER][j] =
					(fit->bitpix == 8) ? round_to_BYTE(g) : round_to_WORD(g);
			fit->pdata[BLAYER][j] =
					(fit->bitpix == 8) ? round_to_BYTE(b) : round_to_WORD(b);
		}
	} else {
		// use librtprocess debayer
		WORD *newbuf = debayer_buffer_new_ushort(buf, &width, &height, interpolation, pattern, xtrans);
		if (!newbuf)
			return 1;

		fit_debayer_buffer(fit, newbuf);
	}

	return 0;
}

static int debayer_float(fits* fit, interpolation_method interpolation, sensor_pattern pattern) {
	int width = fit->rx;
	int height = fit->ry;
	float *buf = fit->fdata;
	int xbayeroff = 0, ybayeroff = 0;

	unsigned int xtrans[6][6];
	if (interpolation == XTRANS) {
		retrieveXTRANSPattern(fit->bayer_pattern, xtrans);
	}

	if (!com.debayer.use_bayer_header) {
		xbayeroff = com.debayer.xbayeroff;
		ybayeroff = com.debayer.ybayeroff;
	} else {
		xbayeroff = fit->bayer_xoffset;
		ybayeroff = fit->bayer_yoffset;
	}

	if (xbayeroff == 1) {
		switch (pattern) {
		case BAYER_FILTER_RGGB:
			pattern = BAYER_FILTER_GRBG;
			break;
		case BAYER_FILTER_BGGR:
			pattern = BAYER_FILTER_GBRG;
			break;
		case BAYER_FILTER_GBRG:
			pattern = BAYER_FILTER_BGGR;
			break;
		case BAYER_FILTER_GRBG:
			pattern = BAYER_FILTER_RGGB;
			break;
		default:
		case BAYER_FILTER_NONE:
			return 1;
		}
	}

	if (ybayeroff == 1) {
		switch (pattern) {
		case BAYER_FILTER_RGGB:
			pattern = BAYER_FILTER_GBRG;
			break;
		case BAYER_FILTER_BGGR:
			pattern = BAYER_FILTER_GRBG;
			break;
		case BAYER_FILTER_GBRG:
			pattern = BAYER_FILTER_RGGB;
			break;
		case BAYER_FILTER_GRBG:
			pattern = BAYER_FILTER_BGGR;
			break;
		default:
		case BAYER_FILTER_NONE:
			return 1;
		}
	}

	float *newbuf = debayer_buffer_new_float(buf, &width, &height, interpolation, pattern, xtrans);
	if (!newbuf)
		return 1;

	fit_debayer_buffer(fit, newbuf);
	return 0;
}

int debayer(fits* fit, interpolation_method interpolation, sensor_pattern pattern) {
	if (fit->type == DATA_USHORT)
		return debayer_ushort(fit, interpolation, pattern);
	else if (fit->type == DATA_FLOAT)
		return debayer_float(fit, interpolation, pattern);
	else return -1;
}

int split_cfa_ushort(fits *in, fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3) {
	int width = in->rx;
	int height = in->ry;
	int j, row, col;

	if (strlen(in->bayer_pattern) > 4) {
		siril_log_message(_("Split CFA does not work on non-Bayer filter camera images!\n"));
		return 1;
	}

	width = width / 2 + width % 2;
	height = height / 2 + height % 2;

	if (new_fit_image(&cfa0, width, height, 1, DATA_USHORT) ||
			new_fit_image(&cfa1, width, height, 1, DATA_USHORT) ||
			new_fit_image(&cfa2, width, height, 1, DATA_USHORT) ||
			new_fit_image(&cfa3, width, height, 1, DATA_USHORT)) {
		return 1;
	}

	WORD c0, c1, c2, c3;
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

int split_cfa_float(fits *in, fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3) {
	int width = in->rx;
	int height = in->ry;
	int j, row, col;

	if (strlen(in->bayer_pattern) > 4) {
		siril_log_message(_("Split CFA does not work on non-Bayer filter camera images!\n"));
		return 1;
	}

	width = width / 2 + width % 2;
	height = height / 2 + height % 2;

	if (new_fit_image(&cfa0, width, height, 1, DATA_FLOAT) ||
			new_fit_image(&cfa1, width, height, 1, DATA_FLOAT) ||
			new_fit_image(&cfa2, width, height, 1, DATA_FLOAT) ||
			new_fit_image(&cfa3, width, height, 1, DATA_FLOAT)) {
		return 1;
	}

	float c0, c1, c2, c3;
	j = 0;

	for (row = 0; row < in->ry - 1; row += 2) {
		for (col = 0; col < in->rx - 1; col += 2) {
			/* not c0, c1, c2 and c3 because of the read orientation */
			c1 = in->fdata[col + row * in->rx];
			c3 = in->fdata[1 + col + row * in->rx];
			c0 = in->fdata[col + (1 + row) * in->rx];
			c2 = in->fdata[1 + col + (1 + row) * in->rx];

			cfa0->fdata[j] = c0;
			cfa1->fdata[j] = c1;
			cfa2->fdata[j] = c2;
			cfa3->fdata[j] = c3;
			j++;
		}
	}

	return 0;
}

int split_cfa_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_) {
	int ret = 1;
	struct split_cfa_data *cfa_args = (struct split_cfa_data *) args->user;

	fits f_cfa0 = { 0 }, f_cfa1 = { 0 }, f_cfa2 = { 0 }, f_cfa3 = { 0 };

	gchar *cfa0 = g_strdup_printf("%s0_%s_%05d%s", cfa_args->seqEntry, cfa_args->seq->seqname, o, com.ext);
	gchar *cfa1 = g_strdup_printf("%s1_%s_%05d%s", cfa_args->seqEntry, cfa_args->seq->seqname, o, com.ext);
	gchar *cfa2 = g_strdup_printf("%s2_%s_%05d%s", cfa_args->seqEntry, cfa_args->seq->seqname, o, com.ext);
	gchar *cfa3 = g_strdup_printf("%s3_%s_%05d%s", cfa_args->seqEntry, cfa_args->seq->seqname, o, com.ext);

	if (fit->type == DATA_USHORT) {
		if (!(ret = split_cfa_ushort(fit, &f_cfa0, &f_cfa1, &f_cfa2, &f_cfa3))) {
			ret = save1fits16(cfa0, &f_cfa0, 0) ||
				save1fits16(cfa1, &f_cfa1, 0) ||
				save1fits16(cfa2, &f_cfa2, 0) ||
				save1fits16(cfa3, &f_cfa3, 0);
		}
	}
	else if (fit->type == DATA_FLOAT) {
		if (!(ret = split_cfa_float(fit, &f_cfa0, &f_cfa1, &f_cfa2, &f_cfa3))) {
			ret = save1fits32(cfa0, &f_cfa0, 0) ||
				save1fits32(cfa1, &f_cfa1, 0) ||
				save1fits32(cfa2, &f_cfa2, 0) ||
				save1fits32(cfa3, &f_cfa3, 0);
		}
	}

	g_free(cfa0); g_free(cfa1);
	g_free(cfa2); g_free(cfa3);
	clearfits(&f_cfa0); clearfits(&f_cfa1);
	clearfits(&f_cfa2); clearfits(&f_cfa3);
	return ret;
}

void apply_split_cfa_to_sequence(struct split_cfa_data *split_cfa_args) {
	struct generic_seq_args *args = malloc(sizeof(struct generic_seq_args));
	args->seq = split_cfa_args->seq;
	args->force_float = FALSE;
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

/******* SPLIT CFA ******************************/

void on_menu_slpitcfa_activate(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("split_cfa_dialog");
}

void on_split_cfa_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("split_cfa_dialog");
}

void on_split_cfa_apply_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton *seq = GTK_TOGGLE_BUTTON(lookup_widget("checkSplitCFASeq"));
	GtkEntry *entrySplitCFA;

	entrySplitCFA = GTK_ENTRY(lookup_widget("entrySplitCFA"));

	if (gtk_toggle_button_get_active(seq) && sequence_is_loaded()) {
		struct split_cfa_data *args = malloc(sizeof(struct split_cfa_data));

		set_cursor_waiting(TRUE);
		args->seq = &com.seq;
		args->seqEntry = gtk_entry_get_text(entrySplitCFA);
		if (args->seqEntry && args->seqEntry[0] == '\0')
			args->seqEntry = "CFA_";
		apply_split_cfa_to_sequence(args);
	} else {
		process_split_cfa(0);
	}
}
