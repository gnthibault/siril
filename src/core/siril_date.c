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
#include <glib.h>
#include <glib/gprintf.h>

#include "siril_date.h"

#define SER_TIME_1970 G_GUINT64_CONSTANT(621355968000000000) // 621.355.968.000.000.000 ticks between 1st Jan 0001 and 1st Jan 1970.

#if !GLIB_CHECK_VERSION(2,62,0)
/**
 * g_date_time_format_iso8601:
 * @datetime: A #GDateTime
 *
 * Format @datetime in [ISO 8601 format](https://en.wikipedia.org/wiki/ISO_8601),
 * including the date, time and time zone, and return that as a UTF-8 encoded
 * string.
 *
 * Since GLib 2.66, this will output to sub-second precision if needed.
 *
 * Returns: a newly allocated string formatted in ISO 8601 format
 *     or %NULL in the case that there was an error. The string
 *     should be freed with g_free().
 * Since: 2.62
 */
static gchar* g_date_time_format_iso8601(GDateTime *datetime) {
	GString *outstr = NULL;
	gchar *main_date = NULL;
	time_t offset;
	gchar *format = "%Y-%m-%dT%H:%M:%S";

	/* if datetime has sub-second non-zero values below the second precision we
	 * should print them as well */
	if (g_date_time_get_microsecond(datetime) % G_TIME_SPAN_SECOND != 0)
		format = "%Y-%m-%dT%H:%M:%S.%f";

	/* Main date and time. */
	main_date = g_date_time_format(datetime, format);
	outstr = g_string_new(main_date);
	g_free(main_date);

	/* Timezone. Format it as `%:::z` unless the offset is zero, in which case
	 * we can simply use `Z`. */
	offset = g_date_time_get_utc_offset(datetime);

	if (offset == 0) {
		g_string_append_c(outstr, 'Z');
	} else {
		gchar *time_zone = g_date_time_format(datetime, "%:::z");
		g_string_append(outstr, time_zone);
		g_free(time_zone);
	}

	return g_string_free(outstr, FALSE);
}
#endif

/**
 * From a datetime it computes the Julian date needed in photometry
 * (code borrowed from muniwin)
 * @param dt timestamp in GDateTime format
 * @return the Julian date
 */
gdouble date_time_to_Julian(GDateTime *dt) {
	gdouble jd1;
	gboolean before;
	gint d1, d2;
	gint year, month, day;
	gint hour, min, sec, ms;

	if (!dt) return 0;

	g_date_time_get_ymd(dt, &year, &month, &day);

	/* Check date and time */
	if (day <= 0 || year <= 0 || month <= 0)
		return 0;

	/* Compute Julian date from input citizen year, month and day. */
	/* Tested for YEAR>0 except 1582-10-07/15 */
	if (year > 1582) {
		before = FALSE;
	} else if (year < 1582) {
		before = TRUE;
	} else if (month > 10) {
		before = FALSE;
	} else if (month < 10) {
		before = TRUE;
	} else if (day >= 15) {
		before = FALSE;
	} else {
		before = TRUE;
	}
	if (month <= 2) {
		d1 = (gint) (365.25 * (year - 1));
		d2 = (gint) (30.6001 * (month + 13));
	} else {
		d1 = (gint) (365.25 * (year));
		d2 = (gint) (30.6001 * (month + 1));
	}

	hour = g_date_time_get_hour(dt);
	min = g_date_time_get_minute(dt);
	sec = g_date_time_get_second(dt);
	ms = g_date_time_get_microsecond(dt) * 0.001;

	jd1 = 1720994.5 + d1 + d2 + day;
	jd1 += 1.0 * hour / 24;
	jd1 += 1.0 * min / 1440.0;
	jd1 += 1.0 * sec / 86400.0;
	jd1 += 1.0 * ms / 86400000.0;

	if (before) {
		return jd1;
	} else {
		return jd1 + 2 - (year / 100) + (year / 400);
	}
}

/**
 * Build filename in the iso8601 format
 * @return a newly allocated string formatted in ISO 8601 format
 *     or %NULL in the case that there was an error. The string
 *     should be freed with g_free().
 */
gchar* build_timestamp_filename() {
	GDateTime *dt = g_date_time_new_now_local();
	gchar *iso8601_string = NULL;
	if (dt) {
		iso8601_string = g_date_time_format_iso8601(dt);
		g_date_time_unref(dt);
	}

	return iso8601_string;
}

/**
 * SER timestamp are converted to GDateTime
 * @param SER timestamp in us
 * @return a newly allocated GDateTime that should be freed
 * with g_date_time_unref().
 */
GDateTime *ser_timestamp_to_date_time(guint64 timestamp) {
	GDateTime *dt, *new_dt = NULL;
	guint64 t1970_ms = (timestamp - SER_TIME_1970) / 10000;
	gint64 secs = t1970_ms / 1000;
	gint ms = t1970_ms % 1000;

	dt = g_date_time_new_from_unix_utc(secs);
	if (dt) {
		/* add milliseconds */
		new_dt = g_date_time_add_seconds(dt, (gdouble) ms * 0.001);

		g_date_time_unref(dt);
	}

	return new_dt;
}

/**
 * From a char * in FITS format to a GDateTime
 * @param date
 * @return a GDateTime or NULL if date is not in the right format
 */
GDateTime *FITS_date_to_date_time(char *date) {
	gint year = 0, month = 0, day = 0, hour = 0, min = 0;
	gdouble sec = 0.0;

	if (date[0] == '\0')
		return NULL;

	if (sscanf(date, "%04d-%02d-%02dT%02d:%02d:%lf", &year, &month, &day, &hour, &min, &sec) != 6) {
		return NULL;
	}
	GTimeZone *tz = g_time_zone_new_utc();
	GDateTime *new_date = g_date_time_new(tz, year, month, day, hour, min, sec);
	g_time_zone_unref(tz);

	return new_date;
}

/**
 *
 * @param date, a GDateTIme
 * @return a newly allocated string formatted as expected in
 * FITS header: "%Y-%m-%dT%H:%M:%S.%f" or NULL in the case that
 * there was an error  The string should be freed with g_free().
 */
gchar *date_time_to_FITS_date(GDateTime *date) {
	return g_date_time_format(date, "%Y-%m-%dT%H:%M:%S.%f");
}
