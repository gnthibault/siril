#ifndef SRC_ALGOS_RGRADIENT_H_
#define SRC_ALGOS_RGRADIENT_H_

#include "core/siril.h"

/* rgradient filter data from GUI */
struct rgradient_filter_data {
	fits *fit;
	double xc, yc, dR, da;
};


gpointer rgradient_filter(gpointer p);

#endif /* SRC_ALGOS_RGRADIENT_H_ */
