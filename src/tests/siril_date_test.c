//#define WITH_MAIN
#ifndef WITH_MAIN
#include <criterion/criterion.h>
#endif

#include "core/siril_date.h"

#define UNDER_US      G_GUINT64_CONSTANT(7)
#define INPUT_TIME    G_GUINT64_CONSTANT(637232717926133380) + UNDER_US
#define SER_TIME_1970 G_GUINT64_CONSTANT(621355968000000000) // 621.355.968.000.000.000 ticks between 1st Jan 0001 and 1st Jan 1970.

/**
 *  Test consistency of siril_date_time functions
 *   */
int test_consistency() {
	GDateTime *dt1, *dt2, *dt3, *dt4;
	guint64 diff, ts;
	gchar *date_str;

	dt1 = ser_timestamp_to_date_time(INPUT_TIME);
	ts = date_time_to_ser_timestamp(dt1);

	/**
	 *  ser timestamp precision is down 0.1 microsecond while our
	 *  structure is accurate down to 1 microsecond
	 */
	diff = INPUT_TIME - ts;
	cr_expect(diff == UNDER_US, "Failed with retval=%lu", diff);

	dt2 = g_date_time_new_from_iso8601 ("2016-11-31T22:10:42Z", NULL);
	ts = date_time_to_ser_timestamp(dt2);
	dt3 = ser_timestamp_to_date_time(ts);
	cr_expect(g_date_time_equal(dt2, dt3), "date_time from ser are not equal");

	/**
	 *  Test FITS date time consistency
	 */
	date_str = date_time_to_FITS_date(dt2);
	dt4 = FITS_date_to_date_time(date_str);
	cr_expect(g_date_time_equal(dt2, dt4), "date_time from FITS are not equal");

	g_date_time_unref(dt1);
	g_date_time_unref(dt2);
	g_date_time_unref(dt3);
	g_date_time_unref(dt4);
	return 0;
}

Test(check_date, test1) { cr_assert(!test_consistency()); }
