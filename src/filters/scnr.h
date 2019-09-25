#ifndef SRC_FILTERS_SCNR_H_
#define SRC_FILTERS_SCNR_H_

#include "core/siril.h"

/* scnr data from GUI */
struct scnr_data {
	fits *fit;
	int type;
	double amount;
	gboolean preserve;
};


gpointer scnr(gpointer p);

#endif /* SRC_FILTERS_SCNR_H_ */
