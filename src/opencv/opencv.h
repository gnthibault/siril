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
int cvResizeGaussian_data8(uint8_t *dataIn, int rx, int ry, uint8_t *dataOut,
		int toX, int toY, int chan, int interpolation);
int cvTranslateImage(fits *image, point shift, int interpolation);
int cvRotateImage(fits *, point, double, int, int);
unsigned char *cvCalculH(s_star *star_array_img,
		struct s_star *star_array_ref, int n, Homography *H);
int cvApplyScaleToH(Homography *H1, double scale);
int cvTransformImage(fits *image, long width, long height, Homography Hom, int interpolation);
int cvUnsharpFilter(fits*, double, double);
int cvComputeFinestScale(fits *image);
int cvLucyRichardson(fits *image, double sigma, int iterations);
int cvClahe(fits *image, double clip_limit, int size);
#ifdef __cplusplus
}
#endif

#endif
