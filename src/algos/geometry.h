#ifndef SRC_ALGOS_GEOMETRY_H_
#define SRC_ALGOS_GEOMETRY_H_

#include "core/siril.h"

/* crop sequence data from GUI */
struct crop_sequence_data {
	sequence *seq;
	rectangle area;
	const char *prefix;
	int retvalue;
};

int verbose_resize_gaussian(fits *, int, int, int);
int verbose_rotate_image(fits *, double, int, int);

void mirrorx(fits *fit, gboolean verbose);
void mirrory(fits *fit, gboolean verbose);

int crop(fits *fit, rectangle *bounds);
gpointer crop_sequence(struct crop_sequence_data *crop_sequence_data);

#endif /* SRC_ALGOS_GEOMETRY_H_ */
