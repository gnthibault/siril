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

#include "core/siril.h"
#include "core/proto.h"
#include "algos/statistics.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "io/single_image.h"


/*****************************************************************************
 *       S I R I L      A R I T H M E T I C      O P E R A T I O N S         *
 ****************************************************************************/

/* equivalent to (map simple_operation a), with simple_operation being
 * (lambda (pixel) (oper pixel scalar))
 * oper is a for addition, s for substraction (i for difference) and so on. */
int soper(fits *a, double scalar, char oper) {
	WORD *gbuf;
	int i, layer;
	int n = a->rx * a->ry;

	g_assert(n > 0);

	for (layer = 0; layer < a->naxes[2]; ++layer) {
		gbuf = a->pdata[layer];
		switch (oper) {
		case OPER_ADD:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD((double) gbuf[i] + scalar);
			}
			break;
		case OPER_SUB:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD((double) gbuf[i] - scalar);
			}
			break;
		case OPER_MUL:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD((double) gbuf[i] * scalar);
			}
			break;
		case OPER_DIV:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD((double) gbuf[i] / scalar);
			}
			break;
		}
	}
	invalidate_stats_from_fit(a);
	return 0;
}

/* applies operation of image a with image b, for all their layers:
 * a = a oper b
 * returns 0 on success */
int imoper(fits *a, fits *b, char oper) {
	int i, layer;

	if (a->rx != b->rx || a->ry != b->ry) {
		siril_log_message(
				_("imoper: images don't have the same size (w = %u|%u, h = %u|%u)\n"),
				a->rx, b->rx, a->ry, b->ry);
		return 1;
	}
	for (layer = 0; layer < a->naxes[2]; ++layer) {
		WORD *buf = b->pdata[layer];
		WORD *gbuf = a->pdata[layer];
		double dbuf, dgbuf;
		for (i = 0; i < a->rx * a->ry; ++i) {
			dbuf = (double) buf[i];
			dgbuf = (double) gbuf[i];
			switch (oper) {
			case OPER_ADD:
				gbuf[i] = round_to_WORD(dgbuf + dbuf);
				break;
			case OPER_SUB:
				gbuf[i] = round_to_WORD(dgbuf - dbuf);
				break;
			case OPER_MUL:
				gbuf[i] = round_to_WORD(dgbuf * dbuf);
				break;
			case OPER_DIV:
				gbuf[i] = (buf[i] == 0) ? 0 : round_to_WORD(dgbuf / dbuf);
				break;
			}
		}
	}
	invalidate_stats_from_fit(a);
	return 0;
}

int addmax(fits *a, fits *b) {
	WORD *gbuf[3] = { a->pdata[RLAYER], a->pdata[GLAYER], a->pdata[BLAYER] };
	WORD *buf[3] = { b->pdata[RLAYER], b->pdata[GLAYER], b->pdata[BLAYER] };
	gint i, layer;

	if (a->rx != b->rx || a->ry != b->ry || a->naxes[2] != b->naxes[2]) {
		siril_log_message(
				_("addmax: images don't have the same size (w = %d|%d, h = %d|%d, layers = %d|%d)\n"),
				a->rx, b->rx, a->ry, b->ry, a->naxes[2], b->naxes[2]);
		return 1;
	}
	g_assert(a->naxes[2] == 1 || a->naxes[2] == 3);

	for (layer = 0; layer < a->naxes[2]; ++layer) {
		for (i = 0; i < a->ry * a->rx; ++i) {
			if (buf[layer][i] > gbuf[layer][i])
				gbuf[layer][i] = buf[layer][i];
		}
	}
	invalidate_stats_from_fit(a);
	return 0;
}

/* a = coef * a / b
 * a is expected as USHORT, returned as FLOAT. b is expected as USHORT.
 * If overflow, siril_fdiv returns 1*/
int siril_fdiv(fits *a, fits *b, float coef) {
	int i, retvalue = 0;
	unsigned long nbdata;
	float *newdata;
	double coefd = (double)coef;

	if (a->type != DATA_USHORT) {
		siril_log_message(_("siril_fdiv: not yet working with 32-bit input data."));
		return -1;
	}
	if (a->rx != b->rx || a->ry != b->ry || a->naxes[2] != b->naxes[2]) {
		fprintf(stderr, "Wrong size or channel count: %u=%u? / %u=%u?\n", a->rx,
				b->rx, a->ry, b->ry);
		return -1;
	}

	nbdata = a->rx * a->ry * a->naxes[2];
	newdata = malloc(nbdata * sizeof(float));
	WORD *abuf = a->data;
	WORD *bbuf = b->data;
	float *fbbuf = b->fdata;

	for (i = 0; i < nbdata; ++i) {
		double val, denominator = 1.0;
		if (b->type == DATA_USHORT && bbuf[i] != 0)
			denominator = (double)bbuf[i];
		else if (b->type == DATA_FLOAT && fbbuf[i] != 0.0f)
			denominator = (double)fbbuf[i];
		if (coef != 1.0f)
			val = coefd * (double)abuf[i] / denominator;
		else val = (double)abuf[i] / denominator;
		if (val > USHRT_MAX_DOUBLE) {
			siril_debug_print("OVERFLOW in FDIV: %lf\n", val);
			retvalue = 1;
		}

		// normalize to float data range [0, 1]
		newdata[i] = (float)(val / USHRT_MAX_DOUBLE);
		if (newdata[i] > 1.0f) newdata[i] = 1.0f;
		//if (newdata[i] < 0.0f) newdata[i] = 0.0f;
	}

	fit_replace_buffer(a, newdata, DATA_FLOAT);
	return retvalue;
}

/* normalized division a/b, stored in a, with max value equal to the original
 * max value of a, for each layer. */
int siril_ndiv(fits *a, fits *b) {
	double *div;
	int layer, i, nb_pixels;
	if (a->rx != b->rx || a->ry != b->ry || a->naxes[2] != b->naxes[2]) {
		fprintf(stderr,
				"Wrong size or channel count: %u=%u? / %u=%u?, %ld=%ld?\n",
				a->rx, b->rx, a->ry, b->ry, a->naxes[2], b->naxes[2]);
		return 1;
	}
	nb_pixels = a->rx * a->ry;
	div = malloc(nb_pixels * sizeof(double));
	if (div == NULL) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	for (layer = 0; layer < a->naxes[2]; ++layer) {
		double max = 0, norm;
		for (i = 0; i < nb_pixels; ++i) {
			if (!b->pdata[layer][i])
				div[i] = (double) a->pdata[layer][i];
			else
				div[i] = (double) a->pdata[layer][i]
						/ (double) b->pdata[layer][i];
			max = max(div[i], max);
		}
		norm = max / fit_get_max(a, layer);
		for (i = 0; i < nb_pixels; ++i) {
			a->pdata[layer][i] = round_to_WORD(div[i] / norm);
		}
	}

	invalidate_stats_from_fit(a);
	free(div);
	return 0;
}
