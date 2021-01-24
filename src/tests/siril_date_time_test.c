//#define WITH_MAIN
#ifndef WITH_MAIN
#include <criterion/criterion.h>
#endif

#include "../core/siril.h"
#include "core/siril_date.h"
#include <stdio.h>

#define INPUT_TIME    G_GUINT64_CONSTANT(637232717926133380)
#define SER_TIME_1970 G_GUINT64_CONSTANT(621355968000000000) // 621.355.968.000.000.000 ticks between 1st Jan 0001 and 1st Jan 1970.

cominfo com;	// the main data struct

/**
 *  Test consistency of siril_date_time functions
 *   */
int test_consistency() {
	GDateTime *dt = ser_timestamp_to_date_time(INPUT_TIME);
	guint64 output = date_time_to_ser_timestamp(dt);

	/**
	 *  ser timestamp precision is down 0.1 microsecond while our
	 *  structure is accurate down to 1 microsecond
	 */
	gint64 retval = INPUT_TIME - output;
	cr_expect(retval == 0, "Failed with retval=%lu", retval);

	gchar *str = date_time_to_FITS_date(dt);
	GDateTime *new_dt = FITS_date_to_date_time(str);

	cr_expect(g_date_time_equal(dt,new_dt), "date_time are not equal");

	g_date_time_unref(dt);
	g_date_time_unref(new_dt);
	return 0;
}

Test(check_date, test1) { cr_assert(!test_consistency()); }
