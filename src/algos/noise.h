#ifndef SRC_ALGOS_NOISE_H_
#define SRC_ALGOS_NOISE_H_

#include "core/siril.h"

/* Noise data from GUI */
struct noise_data {
	gboolean verbose;
	gboolean use_idle;
	fits *fit;
	double bgnoise[3];
	struct timeval t_start;
	int retval;
};

gpointer noise(gpointer p);
void evaluate_noise_in_image();

#endif /* SRC_ALGOS_NOISE_H_ */
