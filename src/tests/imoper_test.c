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

#include <criterion/criterion.h>
#include <math.h>
#include "core/siril.h"
/** I include c file to test static function!! */
#include "core/arithm.c"
/**********************************************/
#include "io/image_format_fits.h"

cominfo com;	// the main data struct
GtkBuilder *builder = NULL;	// get widget references anywhere
fits gfit;	// currently loaded image

/* required links:
 * new_fit_image (io/image_format_fits.o),
 * (io/image_format_fits.o),
 * imoper (core/arithm.o),
 * truncate_to_WORD (core/utils.o)
 * roundf_to_WORD (core/utiles.o)
 * invalidate_stats_from_fit (not important here, dummy.c ?)
 * optional links:
 * siril_log_color_message (dummy.c ?)
 */

static void set_ushort_data(fits *fit, WORD *data, int length) {
	memcpy(fit->data, data, length * sizeof(WORD));
}

static WORD *alloc_data(const WORD *from, int length) {
	WORD *data = malloc(5 * sizeof(WORD));
	memcpy(data, from, length * sizeof(WORD));
	return data;
}

static float *alloc_fdata(const float *from, int length) {
	float *data = malloc(5 * sizeof(float));
	memcpy(data, from, length * sizeof(float));
	return data;
}

void test_ushort() {
	fits *a = NULL, *b = NULL;
	int size = 5;

	WORD origa[] = { 0, 1, 2, 1000, 65535 };
	WORD origb[] = { 2, 2, 2, 2, 2 };

	WORD *dataa = alloc_data(origa, size);
	WORD *datab = alloc_data(origb, size);

	new_fit_image_with_data(&a, size, 1, 1, DATA_USHORT, dataa);
	new_fit_image_with_data(&b, size, 1, 1, DATA_USHORT, datab);

	/* with factor = 1 first */
	int retval = imoper(a, b, OPER_ADD, FALSE);
	cr_assert(!retval, "imoper ADD failed");
	cr_expect_eq(a->data[0], 2);
	cr_expect_eq(a->data[1], 3);
	cr_expect_eq(a->data[2], 4);
	cr_expect_eq(a->data[3], 1002);
	cr_expect_eq(a->data[4], 65535);

	set_ushort_data(a, origa, size);
	retval = imoper(a, b, OPER_SUB, FALSE);
	cr_assert(!retval, "imoper SUB failed");
	cr_expect_eq(a->data[0], 0);
	cr_expect_eq(a->data[1], 0);
	cr_expect_eq(a->data[2], 0);
	cr_expect_eq(a->data[3], 998);
	cr_expect_eq(a->data[4], 65533);

	set_ushort_data(a, origa, size);
	retval = imoper(a, b, OPER_DIV, FALSE);
	cr_assert(!retval, "imoper DIV failed");
	cr_expect_eq(a->data[0], 0);
	cr_expect_eq(a->data[1], 1);
	cr_expect_eq(a->data[2], 1);
	cr_expect_eq(a->data[3], 500);
	cr_expect_eq(a->data[4], 32768);

	set_ushort_data(a, origa, size);
	retval = imoper(a, b, OPER_MUL, FALSE);
	cr_assert(!retval, "imoper MUL failed");
	cr_expect_eq(a->data[0], 0);
	cr_expect_eq(a->data[1], 2);
	cr_expect_eq(a->data[2], 4);
	cr_expect_eq(a->data[3], 2000);
	cr_expect_eq(a->data[4], 65535);

	/* now with factor != 1 */
	set_ushort_data(a, origa, size);
	retval = imoper_with_factor(a, b, OPER_ADD, 2.0f, FALSE);
	cr_assert(!retval, "imoper ADD with factor failed");
	cr_expect_eq(a->data[0], 4);
	cr_expect_eq(a->data[1], 6);
	cr_expect_eq(a->data[2], 8);
	cr_expect_eq(a->data[3], 2004);
	cr_expect_eq(a->data[4], 65535);

	set_ushort_data(a, origa, size);
	retval = imoper_with_factor(a, b, OPER_SUB, 2.0f, FALSE);
	cr_assert(!retval, "imoper SUB with factor failed");
	cr_expect_eq(a->data[0], 0);
	cr_expect_eq(a->data[1], 0);
	cr_expect_eq(a->data[2], 0);
	cr_expect_eq(a->data[3], 1996);
	cr_expect_eq(a->data[4], 65535);

	set_ushort_data(a, origa, size);
	retval = imoper_with_factor(a, b, OPER_DIV, 3.0f, FALSE);
	cr_assert(!retval, "imoper DIV with factor failed");
	cr_expect_eq(a->data[0], 0);
	cr_expect_eq(a->data[1], 2);
	cr_expect_eq(a->data[2], 3);
	cr_expect_eq(a->data[3], 1500);
	cr_expect_eq(a->data[4], 65535);

	set_ushort_data(a, origa, size);
	retval = imoper_with_factor(a, b, OPER_MUL, 0.5f, FALSE);
	cr_assert(!retval, "imoper MUL with factor failed");
	cr_expect_eq(a->data[0], origa[0]);
	cr_expect_eq(a->data[1], origa[1]);
	cr_expect_eq(a->data[2], origa[2]);
	cr_expect_eq(a->data[3], origa[3]);
	cr_expect_eq(a->data[4], origa[4]);

	/* with float output and no factor: this frees a->data each time */
	//WORD origa[] = { 0, 1, 2, 1000, 65535 };
	//WORD origb[] = { 2, 2, 2, 2, 2 };
	set_ushort_data(a, origa, size);
	retval = imoper(a, b, OPER_ADD, TRUE);
	cr_assert(!retval, "imoper ADD to 32 bits failed");
	cr_expect_float_eq(a->fdata[0], 2.0f*INV_USHRT_MAX_SINGLE, 1e-6);
	cr_expect_float_eq(a->fdata[1], 3.0f*INV_USHRT_MAX_SINGLE, 1e-6);
	cr_expect_float_eq(a->fdata[2], 4.0f*INV_USHRT_MAX_SINGLE, 1e-6);
	cr_expect_float_eq(a->fdata[3], 1002.0f*INV_USHRT_MAX_SINGLE, 1e-6);
	cr_expect_float_eq(a->fdata[4], 1.0f, 1e-7);

	clearfits(a);
	dataa = alloc_data(origa, size);
	new_fit_image_with_data(&a, size, 1, 1, DATA_USHORT, dataa);
	retval = imoper(a, b, OPER_SUB, TRUE);
	cr_assert(!retval, "imoper SUB to 32 bits failed");
	cr_expect_float_eq(a->fdata[0], -2.0f*INV_USHRT_MAX_SINGLE, 1e-6, "expected -2 as float, got %f", a->fdata[0]);
	cr_expect_float_eq(a->fdata[1], -1.0f*INV_USHRT_MAX_SINGLE, 1e-6, "expected -1 as float, got %f", a->fdata[1]);
	cr_expect_float_eq(a->fdata[2], 0.0f, 1e-7, "expected 0, got %f", a->fdata[2]);
	cr_expect_float_eq(a->fdata[3], 998.f*INV_USHRT_MAX_SINGLE, 1e-6, "expected 998 as float, got %f", a->fdata[3]);
	cr_expect_float_eq(a->fdata[4], 65533.f*INV_USHRT_MAX_SINGLE, 1e-6, "expected 65535 as float, got %f", a->fdata[4]);

	clearfits(a);
	dataa = alloc_data(origa, size);
	new_fit_image_with_data(&a, size, 1, 1, DATA_USHORT, dataa);
	retval = imoper(a, b, OPER_DIV, TRUE);
	cr_assert(!retval, "imoper DIV to 32 bits failed");
	cr_expect_float_eq(a->fdata[0], 0.0f, 1e-7, "expected 0, got %f", a->fdata[0]);
	cr_expect_float_eq(a->fdata[1], 0.5f, 1e-6, "expected 0.5, got %f", a->fdata[1]);
	cr_expect_float_eq(a->fdata[2], 1.0f, 1e-6, "expected 1, got %f", a->fdata[2]);
	cr_expect_float_eq(a->fdata[3], 1.0f, 1e-6, "expected 1, got %f", a->fdata[3]);
	cr_expect_float_eq(a->fdata[4], 1.0f, 1e-6, "expected 1, got %f", a->fdata[4]);

	/* the expected results for division and multiplication are unclear in
	 * that case, for example, 1 * 65535 converted to float would be
	 * .000015259022 because 65535 becomes 1 and 1 becomes .000015259022 */
	/*clearfits(a);
	dataa = alloc_data(origa, size);
	new_fit_image_with_data(&a, size, 1, 1, DATA_USHORT, dataa);
	retval = imoper(a, b, OPER_MUL, TRUE);
	cr_assert(!retval, "imoper MUL to 32 bits failed");
	cr_expect_float_eq(a->fdata[0], 0.0f, 1e-7, "expected 0, got %f", a->fdata[0]);
	cr_expect_float_eq(a->fdata[1], 2.0f*INV_USHRT_MAX_SINGLE, 1e-6, "expected 2 as float, got %f", a->fdata[1]);
	cr_expect_float_eq(a->fdata[2], 4.0f*INV_USHRT_MAX_SINGLE, 1e-6, "expected 4 as float, got %f", a->fdata[2]);
	cr_expect_float_eq(a->fdata[3], 2000.0f*INV_USHRT_MAX_SINGLE, 1e-6, "expected 2000 as float, got %f", a->fdata[3]);
	cr_expect_float_eq(a->fdata[4], 1.0f, 1e-6, "expected 1, got %f", a->fdata[4]);*/

	free(datab);
}

void test_ushort_float() {
	fits *a = NULL, *b = NULL;
	int size = 5;

	WORD origa[] = { 0, 1, 2, 1000, 65535 };
	float origb[] = { 0.1f, 0.1f, 0.1f, 0.1f, 0.1f };

	WORD *dataa = alloc_data(origa, size);
	float *datab = alloc_fdata(origb, size);

	new_fit_image_with_data(&a, size, 1, 1, DATA_USHORT, dataa);
	new_fit_image_with_data(&b, size, 1, 1, DATA_FLOAT, datab);

	/* with factor = 1 first */
	int retval = imoper(a, b, OPER_ADD, FALSE);
	cr_assert(!retval, "imoper failed");
	cr_expect_eq(a->data[0], 6553);
	cr_expect_eq(a->data[1], 6554);
	cr_expect_eq(a->data[2], 6555);
	cr_expect_eq(a->data[3], 7553);
	cr_expect_eq(a->data[4], 65535);

	set_ushort_data(a, origa, size);
	retval = imoper(a, b, OPER_SUB, FALSE);
	cr_assert(!retval, "imoper failed");
	cr_expect_eq(a->data[0], 0);
	cr_expect_eq(a->data[1], 0);
	cr_expect_eq(a->data[2], 0);
	cr_expect_eq(a->data[3], 0);
	cr_expect_eq(a->data[4], 58982);

	set_ushort_data(a, origa, size);
	retval = imoper(a, b, OPER_DIV, FALSE);
	cr_assert(!retval, "imoper failed");
	cr_expect_eq(a->data[0], 0);
	cr_expect_eq(a->data[1], 0);
	cr_expect_eq(a->data[2], 0);
	cr_expect_eq(a->data[3], 0);
	cr_expect_eq(a->data[4], 10);

	set_ushort_data(a, origa, size);
	retval = imoper(a, b, OPER_MUL, FALSE);
	cr_assert(!retval, "imoper failed");
	cr_expect_eq(a->data[0], 0);
	cr_expect_eq(a->data[1], 6554); // TODO: not consistent with addition: truncate vs. roundf_to_WORD
	cr_expect_eq(a->data[2], 13107);// TODO: not consistent with previous line: truncate vs. roundf_to_WORD
	cr_expect_eq(a->data[3], 65535);
	cr_expect_eq(a->data[4], 65535);


	free(dataa);
	free(datab);
}

Test(arithmetics, ushort) { test_ushort(); }
Test(arithmetics, float) { test_ushort_float(); }
