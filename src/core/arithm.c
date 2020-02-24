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

#include "core/siril.h"
#include "core/proto.h"
#include "algos/statistics.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "io/single_image.h"


/*****************************************************************************
 *       S I R I L      A R I T H M E T I C      O P E R A T I O N S         *
 ****************************************************************************/

static int soper_ushort_to_ushort(fits *a, double scalar, image_operator oper) {
	WORD *data;
	long i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	g_assert(n > 0);
	data = a->data;

	switch (oper) {
		case OPER_ADD:
			for (i = 0; i < n; ++i) {
				data[i] = round_to_WORD((double)data[i] + scalar);
			}
			break;
		case OPER_SUB:
			for (i = 0; i < n; ++i) {
				data[i] = round_to_WORD((double)data[i] - scalar);
			}
			break;
		case OPER_MUL:
			for (i = 0; i < n; ++i) {
				data[i] = round_to_WORD((double)data[i] * scalar);
			}
			break;
		case OPER_DIV:
			for (i = 0; i < n; ++i) {
				data[i] = round_to_WORD((double)data[i] / scalar);
			}
			break;
	}
	invalidate_stats_from_fit(a);
	return 0;
}

static int soper_ushort_to_float(fits *a, double scalar, image_operator oper) {
	WORD *data;
	float *result;
	long i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	g_assert(n > 0);
	data = a->data;
	result = malloc(n * sizeof(float));
	if (!result) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	switch (oper) {
		case OPER_ADD:
			for (i = 0; i < n; ++i) {
				result[i] = double_ushort_to_float_range((double)data[i] + scalar);
			}
			break;
		case OPER_SUB:
			for (i = 0; i < n; ++i) {
				result[i] = double_ushort_to_float_range((double)data[i] - scalar);
			}
			break;
		case OPER_MUL:
			for (i = 0; i < n; ++i) {
				result[i] = double_ushort_to_float_range((double)data[i] * scalar);
			}
			break;
		case OPER_DIV:
			for (i = 0; i < n; ++i) {
				result[i] = double_ushort_to_float_range((double)data[i] / scalar);
			}
			break;
	}
	fit_replace_buffer(a, result, DATA_FLOAT);
	return 0;
}

static int soper_float(fits *a, double scalar, image_operator oper) {
	float *data;
	long i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	g_assert(n > 0);
	data = a->fdata;

	switch (oper) {
		case OPER_ADD:
			for (i = 0; i < n; ++i) {
				data[i] = (float)((double)data[i] + scalar);
			}
			break;
		case OPER_SUB:
			for (i = 0; i < n; ++i) {
				data[i] = (float)((double)data[i] - scalar);
			}
			break;
		case OPER_MUL:
			for (i = 0; i < n; ++i) {
				data[i] = (float)((double)data[i] * scalar);
			}
			break;
		case OPER_DIV:
			for (i = 0; i < n; ++i) {
				data[i] = (float)((double)data[i] / scalar);
			}
			break;
	}
	invalidate_stats_from_fit(a);
	return 0;
}

/* equivalent to (map simple_operation a), with simple_operation being
 * (lambda (pixel) (oper pixel scalar))
 * If oper is ADD or SUB, it will be applied on the original type of the image,
 * in case conversion to 32 bits is required.
 */
int soper(fits *a, double scalar, image_operator oper, gboolean conv_to_float) {
	if (a->type == DATA_USHORT) {
		if (conv_to_float)
			return soper_ushort_to_float(a, scalar, oper);
		return soper_ushort_to_ushort(a, scalar, oper);
	}
	if (a->type == DATA_FLOAT)
		return soper_float(a, scalar, oper);
	return 1;
}

static int imoper_ushort_to_ushort(fits *a, fits *b, image_operator oper, float factor) {
	long i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];

	if (memcmp(a->naxes, b->naxes, sizeof a->naxes)) {
		siril_log_message(_("imoper: images must have same dimensions\n"));
		return 1;
	}

	WORD *abuf = a->data, *bbuf = b->data;
	if (oper == OPER_DIV) {
		for (i = 0; i < n; ++i) {
			if (bbuf[i] == 0)
				abuf[i] = 0;
			else {
				float aval = (float)abuf[i];
				float bval = (float)bbuf[i];
				if (factor != 1.0f)
					abuf[i] = round_to_WORD(factor * (aval / bval));
				else abuf[i] = round_to_WORD(aval / bval);
			}
		}
	} else {
		for (i = 0; i < n; ++i) {
			int aval = (int)abuf[i];
			int bval = (int)bbuf[i];
			switch (oper) {
				case OPER_ADD:
					abuf[i] = truncate_to_WORD(aval + bval);
					break;
				case OPER_SUB:
					abuf[i] = truncate_to_WORD(aval - bval);
					break;
				case OPER_MUL:
					abuf[i] = truncate_to_WORD(aval * bval);
					break;
				case OPER_DIV:	// handled above
					break;
			}
			if (factor != 1.0f)
				abuf[i] = round_to_WORD(factor * (float)abuf[i]);
		}
	}
	invalidate_stats_from_fit(a);
	return 0;
}

int imoper_to_float(fits *a, fits *b, image_operator oper, float factor) {
	long i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];
	float *result;

	if (memcmp(a->naxes, b->naxes, sizeof a->naxes)) {
		siril_log_message(_("imoper: images must have same dimensions\n"));
		return 1;
	}

	if (a->type == DATA_FLOAT) {
		result = a->fdata;
		if (factor > 1.0) {
			factor /= USHRT_MAX_SINGLE;
		}
	} else {
		result = malloc(n * sizeof(float));
		if (!result) {
			PRINT_ALLOC_ERR;
			return 1;
		}
	}

	for (i = 0; i < n; ++i) {
		float aval = a->type == DATA_USHORT ? ushort_to_float_range(a->data[i]) : a->fdata[i];
		float bval = b->type == DATA_USHORT ? ushort_to_float_range(b->data[i]) : b->fdata[i];
		switch (oper) {
			case OPER_ADD:
				result[i] = aval + bval;
				break;
			case OPER_SUB:
				result[i] = aval - bval;
				break;
			case OPER_MUL:
				result[i] = aval * bval;
				break;
			case OPER_DIV:
				printf("aval=%g\tbval=%g\tfactor=%g\n", aval, bval, factor);
				if (bval == 0.0f)
					result[i] = 0.0f;
				else result[i] = aval / bval;
		}
		if (factor != 1.0f)
			result[i] *= factor;
		if (result[i] > 1.0f)	// should we truncate by default?
			result[i] = 1.0f;
	}
	if (a->type == DATA_USHORT) {
		fit_replace_buffer(a, result, DATA_FLOAT);
	} else invalidate_stats_from_fit(a);
	return 0;
}

/* applies operation of image a with image b, for all their layers:
 * a = factor * a oper b
 * returns 0 on success */
static int imoper_with_factor(fits *a, fits *b, image_operator oper, float factor, gboolean allow_32bits) {
	// ushort result can only be forced when both input images are ushort
	if (allow_32bits)
		return imoper_to_float(a, b, oper, factor);
	else {
		if (a->type == DATA_USHORT && b->type == DATA_USHORT)
			return imoper_ushort_to_ushort(a, b, oper, 1.0f);
		siril_log_color_message(_("image operations can only be kept 16 bits if the two input images are 16 bits too. Aborting.\n"), "red");
	}
	return 1;
}

int imoper(fits *a, fits *b, image_operator oper, gboolean allow_32bits) {
	return imoper_with_factor(a, b, oper, 1.0f, allow_32bits);
}

/* a = coef * a / b
 * a is expected as USHORT, returned as FLOAT. b is expected as USHORT.
 * If overflow, siril_fdiv returns 1*/
int siril_fdiv(fits *a, fits *b, float coef, gboolean allow_32bits) {
	return imoper_with_factor(a, b, OPER_DIV, coef, allow_32bits);
}

// a = max(a, b)
int addmax(fits *a, fits *b) {
	long i, n = a->naxes[0] * a->naxes[1] * a->naxes[2];

	if (memcmp(a->naxes, b->naxes, sizeof a->naxes)) {
		siril_log_message(_("addmax: images must have same dimensions\n"));
		return 1;
	}
	if (a->type != b->type) {
		siril_log_message(_("addmax: images must have same data type\n"));
		return 1;
	}
	g_assert(a->naxes[2] == 1 || a->naxes[2] == 3);

	if (a->type == DATA_USHORT) {
		WORD *abuf = a->data, *bbuf = b->data;
		for (i = 0; i < n; ++i) {
			if (bbuf[i] > abuf[i])
				abuf[i] = bbuf[i];
		}
	} else {
		float *abuf = a->fdata, *bbuf = b->fdata;
		for (i = 0; i < n; ++i) {
			if (bbuf[i] > abuf[i])
				abuf[i] = bbuf[i];
		}
	}
	invalidate_stats_from_fit(a);
	return 0;
}

