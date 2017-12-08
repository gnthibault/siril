#ifndef SIRIL_OPENCV_H_
#define SIRIL_OPENCV_H_

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#ifdef HAVE_OPENCV

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include "registration/matching/misc.h"
#include "registration/matching/atpmatch.h"

int cvResizeGaussian(fits *, int, int, int);
int cvResizeGaussian_data8(uint8_t *dataIn, int rx, int ry, uint8_t *dataOut,
		int toX, int toY, int chan, int interpolation);
int cvRotateImage(fits *, double, int, int);
int cvCalculH(s_star *star_array_img,
		struct s_star *star_array_ref, int n, Homography *H);
int cvTransformImage(fits *, point, Homography, int);
int cvUnsharpFilter(fits*, double, double);
int cvComputeFinestScale(fits *image);
int cvLucyRichardson(fits *image, double sigma, int iterations);
#ifdef __cplusplus
}
#endif

#endif	/* HAVE_OPENCV */

#endif
