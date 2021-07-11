#ifndef SIRIL_OPENCV_H_
#define SIRIL_OPENCV_H_

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include "registration/matching/misc.h"
#include "registration/matching/atpmatch.h"
#include "gui/progress_and_log.h"

int cvResizeGaussian(fits *, int, int, int);

int cvRotateImage(fits *, point, double, int, int);

int cvAffineTransformation(fits *image, pointf *refpoints, pointf *curpoints, int nb_points,
		gboolean upscale2x, int interpolation);

unsigned char *cvCalculH(s_star *star_array_img,
		struct s_star *star_array_ref, int n, Homography *H);


int cvTransformImage(fits *image, unsigned int width, unsigned int height, Homography Hom, gboolean upscale2x, int interpolation);

int cvUnsharpFilter(fits*, double, double);

int cvClahe(fits *image, double clip_limit, int size);

#ifdef __cplusplus
}
#endif

#endif
