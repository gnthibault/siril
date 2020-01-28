#ifndef SRC_ALGOS_DECONVOLUTION_H_
#define SRC_ALGOS_DECONVOLUTION_H_

#include "core/siril.h"

/* deconvolution filter data from GUI */
struct deconvolution_filter_data {
	fits *fit;
	size_t contrastThreshold;
	double sigma, deconvSigmaOffset;
	size_t iterations;
	gboolean iterCheck, showMask;
};

gpointer deconvolution(gpointer p);
//gpointer rgradient_filter(gpointer p);

#endif /* SRC_ALGOS_DECONVOLUTION_H_ */
