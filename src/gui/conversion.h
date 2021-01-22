#ifndef _CONVERSION_GUI_H
#define _CONVERSION_GUI_H

#include <glib.h>

enum {
	COLUMN_FILENAME,	// string
	COLUMN_SIZE,		// string
	COLUMN_DATE,		// string
	N_COLUMNS_CONVERT
};

void fill_convert_list(GSList *list);

#endif
