#ifndef _CONVERSION_GUI_H
#define _CONVERSION_GUI_H

#include <glib.h>

enum {
	COLUMN_FILENAME,	// string
	COLUMN_SIZE,		// string
	COLUMN_SIZE_INT64,	// gint64
	COLUMN_DATE,		// string
	COLUMN_DATE_UNIX,	// guint64
	N_COLUMNS_CONVERT
};

void fill_convert_list(GSList *list);

#endif
