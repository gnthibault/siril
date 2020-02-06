#ifndef SRC_FILTERS_DECONV_H_
#define SRC_FILTERS_DECONV_H_

#include "core/siril.h"

/* Lucy-Richardson data from GUI */
struct deconv_data {
	fits *fit;
	double clip;
	double sigma;
	double corner_radius;
	size_t contrast_threshold;
	size_t iterations;
	gboolean auto_limit;
	gboolean auto_contrast_threshold;
};

void apply_deconv_cancel();
gpointer RTdeconv(gpointer p);
#ifdef __cplusplus
extern "C" {
#endif
gpointer deconvolution(gpointer p);
#ifdef __cplusplus
}
#endif


#endif /* SRC_FILTERS_DECONV_H_ */
