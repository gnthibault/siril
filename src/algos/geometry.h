#ifndef SRC_ALGOS_GEOMETRY_H_
#define SRC_ALGOS_GEOMETRY_H_

#include "core/siril.h"

int	verbose_resize_gaussian(fits *, int, int, int);
int	verbose_rotate_image(fits *, double, int, int);

void 	mirrorx(fits *fit, gboolean verbose);
void 	mirrory(fits *fit, gboolean verbose);

int	crop(fits *fit, rectangle *bounds);

#endif /* SRC_ALGOS_GEOMETRY_H_ */
