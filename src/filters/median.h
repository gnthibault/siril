#ifndef SRC_FILTERS_MEDIAN_H_
#define SRC_FILTERS_MEDIAN_H_

#include <glib.h>

/* median filter data from GUI */
struct median_filter_data {
	fits *fit;
	int ksize;
	double amount;
	int iterations;
};

gpointer median_filter(gpointer p);

#endif /* SRC_FILTERS_MEDIAN_H_ */
