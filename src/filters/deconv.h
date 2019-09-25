#ifndef SRC_FILTERS_DECONV_H_
#define SRC_FILTERS_DECONV_H_

#include "core/siril.h"

/* Lucy-Richardson data from GUI */
struct RL_data {
	fits *fit;
	double sigma;
	int iter;
};

gpointer LRdeconv(gpointer p);


#endif /* SRC_FILTERS_DECONV_H_ */
