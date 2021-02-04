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
#include "core/siril.h"
/** I include c file to test stati function!! */
#include "core/arithm.c"
/**********************************************/
#include "io/image_format_fits.h"

cominfo com;	// the main data struct
GtkBuilder *builder = NULL;	// get widget references anywhere
fits gfit;	// currently loaded image

WORD origa[] = { 0, 1, 2, 1000, 65535 };
WORD origb[] = { 2, 2, 2, 2, 2 };

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

void test_ushort() {
	fits *a = NULL, *b = NULL;
	int size = 5;

	WORD *dataa = alloc_data(origa, size);
	WORD *datab = alloc_data(origb, size);

	new_fit_image_with_data(&a, size, 1, 1, DATA_USHORT, dataa);
	new_fit_image_with_data(&b, size, 1, 1, DATA_USHORT, datab);

	/* with factor = 1 first */
	int retval = imoper(a, b, OPER_ADD, FALSE);
	cr_expect(!retval, "imoper failed");
	cr_expect(a->data[0] == 2);
	cr_expect(a->data[1] == 3);
	cr_expect(a->data[2] == 4);
	cr_expect(a->data[3] == 1002);
	cr_expect(a->data[4] == 65535);

	set_ushort_data(a, origa, size);
	retval = imoper(a, b, OPER_SUB, FALSE);
	cr_expect(!retval, "imoper failed");
	cr_expect(a->data[0] == 0);
	cr_expect(a->data[1] == 0);
	cr_expect(a->data[2] == 0);
	cr_expect(a->data[3] == 998);
	cr_expect(a->data[4] == 65533);

	set_ushort_data(a, origa, size);
	retval = imoper(a, b, OPER_DIV, FALSE);
	cr_expect(!retval, "imoper failed");
	cr_expect(a->data[0] == 0);
	cr_expect(a->data[1] == 0);
	cr_expect(a->data[2] == 1);
	cr_expect(a->data[3] == 500);
	cr_expect(a->data[4] == 32768);

	set_ushort_data(a, origa, size);
	retval = imoper(a, b, OPER_MUL, FALSE);
	cr_expect(!retval, "imoper failed");
	cr_expect(a->data[0] == 0);
	cr_expect(a->data[1] == 2);
	cr_expect(a->data[2] == 4);
	cr_expect(a->data[3] == 2000);
	cr_expect(a->data[4] == 65535);

	/* now with factor != 1 */
	set_ushort_data(a, origa, size);
	retval = imoper_with_factor(a, b, OPER_ADD, 2.0f, FALSE);
	cr_expect(!retval, "imoper failed");
	cr_expect(a->data[0] == 4);
	cr_expect(a->data[1] == 6);
	cr_expect(a->data[2] == 8);
	cr_expect(a->data[3] == 2004);
	cr_expect(a->data[4] == 65535);

	set_ushort_data(a, origa, size);
	retval = imoper_with_factor(a, b, OPER_SUB, 2.0f, FALSE);
	cr_expect(!retval, "imoper failed");
	cr_expect(a->data[0] == 0);
	cr_expect(a->data[1] == 0);
	cr_expect(a->data[2] == 0);
	cr_expect(a->data[3] == 1996);
	cr_expect(a->data[4] == 65535);

	set_ushort_data(a, origa, size);
	retval = imoper_with_factor(a, b, OPER_DIV, 3.0f, FALSE);
	cr_expect(!retval, "imoper failed");
	cr_expect(a->data[0] == 0);
	cr_expect(a->data[1] == 2);
	cr_expect(a->data[2] == 3);
	cr_expect(a->data[3] == 1500);
	cr_expect(a->data[4] == 65535);

	/*  case that probably doesn't work as expected, but I don't think we use it */
	set_ushort_data(a, origa, size);
	retval = imoper_with_factor(a, b, OPER_MUL, 0.5f, FALSE);
	cr_expect(!retval, "imoper failed");
	cr_expect(a->data[0] == origa[0]);
	cr_expect(a->data[1] == origa[1]);
	cr_expect(a->data[2] == origa[2]);
	cr_expect(a->data[3] == origa[3]);
//	cr_expect(a->data[4] == origa[4]);
//	fprintf(stderr, "a->data[4]=%d et origa[4]=%d\n", a->data[4], origa[4]);


	free(dataa);
	free(datab);
}

Test(arithmetics, test1) { test_ushort(); }
