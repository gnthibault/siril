#ifndef _ECC_H_
#define _ECC_H_

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#ifdef HAVE_OPENCV

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Enhanced Correlation Coefficient (ECC) Image Alignment Algorithm
 * http://xanthippi.ceid.upatras.gr/people/evangelidis/ecc
 *
 * Authors: Georgios Evangelidis
 * e-mail: evagelid@ceid.upatras.gr
 * e-mail: georgios.evangelidis@inria.fr
 *
 * Copyright (2010): G.D.Evangelidis
 *
 * For details, see the paper:
 * G.D. Evangelidis, E.Z. Psarakis, "Parametric Image Alignment using
 * Enhanced Correlation Coefficient Maximization, IEEE Transaction on
 * Pattern Analysis & Machine Intelligence, vol. 30, no. 10, 2008
 */


typedef enum {
	WARP_MODE_TRANSLATION,
	WARP_MODE_EUCLIDEAN,
	WARP_MODE_AFFINE,
	WARP_MODE_HOMOGRAPHY
} WARP_MODE;

typedef  struct reg_ecc_struct reg_ecc;

struct reg_ecc_struct {
	float dx;
	float dy;
};

int findTransform(fits *reference, fits *image, int layer, reg_ecc *reg_param);

#ifdef __cplusplus
}
#endif

#endif	/* HAVE_OPENCV */


#endif //_ECC_H_
