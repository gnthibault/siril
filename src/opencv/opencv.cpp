/*
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2015 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 *
 * Useful links about OpenCV:
 * http://docs.opencv.org/modules/core/doc/intro.html
 * http://docs.opencv.org/modules/imgproc/doc/geometric_transformations.html#resize
 */
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#ifdef HAVE_OPENCV
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "core/siril.h"
#include "core/proto.h"
#include "registration/matching/misc.h"
#include "opencv.h"
#include "opencv/ecc/ecc.h"

using namespace cv;

static WORD *fits_to_bgrbgr(fits *image) {
	int ndata = image->rx * image->ry;
	WORD *bgrbgr = new WORD[ndata * 3];
	for (int i = 0, j = 0; i < ndata * 3; i += 3, j++) {
		bgrbgr[i + 0] = image->pdata[BLAYER][j];
		bgrbgr[i + 1] = image->pdata[GLAYER][j];
		bgrbgr[i + 2] = image->pdata[RLAYER][j];
	}
	return bgrbgr;
}

int cvResizeGaussian_data8(uint8_t *dataIn, int rx, int ry, uint8_t *dataOut,
		int toX, int toY, int chan, int interpolation) {
	int mode = (chan == 1 ? CV_8UC1 : CV_8UC3);

	Mat in(ry, rx, mode, dataIn);
	Mat out(toY, toX, mode);

	resize(in, out, out.size(), 0, 0, interpolation);

	for (int i = 0; i < toX * toY * chan; i++)
		dataOut[i] = (uint8_t) out.data[i];

	in = Mat();
	out = Mat();
	return 0;
}

/* resizes image to the sizes toX * toY, and stores it back in image */
int cvResizeGaussian(fits *image, int toX, int toY, int interpolation) {
	assert(image->data);
	assert(image->rx);

	WORD *bgrbgr = fits_to_bgrbgr(image);

	Mat in(image->ry, image->rx, CV_16UC3, bgrbgr);
	Mat out(toY, toX, CV_16UC3);

	resize(in, out, out.size(), 0, 0, interpolation);

	image->rx = toX;
	image->naxes[0] = toX;
	image->ry = toY;
	image->naxes[1] = toY;
	WORD *newdata = (WORD*) realloc(image->data,
			toX * toY * sizeof(WORD) * image->naxes[2]);
	if (!newdata) {
		free(newdata);
		return 1;
	}
	image->data = newdata;
	Mat channel[3];
	split(out, channel);

	memcpy(image->data, channel[2].data, toX * toY * sizeof(WORD));
	if (image->naxes[2] == 3) {
		memcpy(image->data + toX * toY, channel[1].data,
				toX * toY * sizeof(WORD));
		memcpy(image->data + toX * toY * 2, channel[0].data,
				toX * toY * sizeof(WORD));
	}

	if (image->naxes[2] == 1) {
		image->pdata[0] = image->data;
		image->pdata[1] = image->data;
		image->pdata[2] = image->data;
	} else {
		image->pdata[0] = image->data;
		image->pdata[1] = image->data + (toX * toY);
		image->pdata[2] = image->data + (toX * toY) * 2;
	}
	delete[] bgrbgr;
	in = Mat();
	out = Mat();
	return 0;
}

/* Rotate an image with the angle "angle" */
int cvRotateImage(fits *image, double angle, int interpolation, int cropped) {
	assert(image->data);
	assert(image->rx);
	assert(image->ry);

	int ndata = image->rx * image->ry;

	WORD *bgrbgr = fits_to_bgrbgr(image);

	Mat in(image->ry, image->rx, CV_16UC3, bgrbgr);
	Mat out(image->ry, image->rx, CV_16UC3);

	if ((fmod(angle, 90.0) == 0) && interpolation == -1) {	// fast rotation
		transpose(in, out);
		if (angle == 90.0)
			flip(out, out, 0);
		else // 270, -90
			flip(out, out, 1);
	} else {
		Point2f pt(in.cols / 2.0, in.rows / 2.0);// We take the center of the image. Should we pass this in function parameters ?
		Mat r = getRotationMatrix2D(pt, angle, 1.0);
		if (cropped == 1) {
			warpAffine(in, out, r, in.size(), interpolation);
		} else {

			// determine bounding rectangle
			Rect frame = RotatedRect(pt, in.size(), angle).boundingRect();
			// adjust transformation matrix
			r.at<double>(0, 2) += frame.width / 2.0 - pt.x;
			r.at<double>(1, 2) += frame.height / 2.0 - pt.y;

			warpAffine(in, out, r, frame.size(), interpolation);
			ndata = out.cols * out.rows;
			WORD *newdata = (WORD*) realloc(image->data,
					ndata * image->naxes[2] * sizeof(WORD));
			if (!newdata) {
				free(newdata);
				return 1;
			}
			image->data = newdata;
		}
	}
	Mat channel[3];
	split(out, channel);

	memcpy(image->data, channel[2].data, ndata * sizeof(WORD));
	if (image->naxes[2] == 3) {
		memcpy(image->data + ndata, channel[1].data, ndata * sizeof(WORD));
		memcpy(image->data + ndata * 2, channel[0].data, ndata * sizeof(WORD));
	}

	if (image->naxes[2] == 1) {
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data;
		image->pdata[BLAYER] = image->data;
	} else {
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data + ndata;
		image->pdata[BLAYER] = image->data + ndata * 2;
	}
	image->rx = out.cols;
	image->ry = out.rows;
	image->naxes[0] = image->rx;
	image->naxes[1] = image->ry;

	/* free data */
	delete[] bgrbgr;
	in = Mat();
	out = Mat();
	return 0;
}

int cvTransformImage(fits *image, TRANS trans, int interpolation) {
	assert(image->data);
	assert(image->rx);
	assert(image->ry);

	int ndata = image->rx * image->ry;

	WORD *bgrbgr = fits_to_bgrbgr(image);

	Mat in(image->ry, image->rx, CV_16UC3, bgrbgr);
	Mat out(image->ry, image->rx, CV_16UC3);

	double angle = -atan2(trans.c, trans.b);
	double s = sqrt(trans.b * trans.b + trans.c * trans.c);

    // http://en.wikipedia.org/wiki/Transformation_matrix#Affine_transformations
	Mat transform = Mat::eye(2, 3, CV_64FC1);
	transform.at<double>(0, 0) = s * cos(angle);
	transform.at<double>(0, 1) = -s * sin(angle);
	transform.at<double>(1, 0) = s * sin(angle);
	transform.at<double>(1, 1) = s * cos(angle);
	transform.at<double>(0, 2) = trans.a;	// shift dx
	transform.at<double>(1, 2) = trans.d;	// shift dy

	warpAffine(in, out, transform, in.size(), interpolation);

	Mat channel[3];
	split(out, channel);

	memcpy(image->data, channel[2].data, ndata * sizeof(WORD));
	if (image->naxes[2] == 3) {
		memcpy(image->data + ndata, channel[1].data, ndata * sizeof(WORD));
		memcpy(image->data + ndata * 2, channel[0].data, ndata * sizeof(WORD));
	}

	if (image->naxes[2] == 1) {
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data;
		image->pdata[BLAYER] = image->data;
	} else {
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data + ndata;
		image->pdata[BLAYER] = image->data + ndata * 2;
	}
	image->rx = out.cols;
	image->ry = out.rows;
	image->naxes[0] = image->rx;
	image->naxes[1] = image->ry;

	/* free data */
	delete[] bgrbgr;
	in = Mat();
	out = Mat();
	transform = Mat();
	return 0;
}

int cvUnsharpFilter(fits* image, double sigma, double amount) {
	assert(image->data);
	assert(image->rx);
	int type = CV_16U;
	if (image->naxes[2] != 1)
		type = CV_16UC3;

	Mat in(image->ry, image->rx, type, image->data);
	Mat out, contrast;
	GaussianBlur(in, out, Size(), sigma);
	if (fabs(amount) > 0.0) {
		Mat sharpened = in * (1 + amount) + out * (-amount);
		out = sharpened.clone();
		sharpened.release();
	}

	memcpy(image->data, out.data,
			image->rx * image->ry * sizeof(WORD) * image->naxes[2]);
	if (image->naxes[2] == 1) {
		image->pdata[0] = image->data;
		image->pdata[1] = image->data;
		image->pdata[2] = image->data;
	} else {
		image->pdata[0] = image->data;
		image->pdata[1] = image->data + (image->rx * image->ry);
		image->pdata[2] = image->data + (image->rx * image->ry) * 2;
	}
	in.release();
	return 0;
}

int cvComputeFinestScale(fits *image) {
	assert(image->data);
	assert(image->rx);
	assert(image->ry);

	int ndata = image->rx * image->ry;

	WORD *bgrbgr = fits_to_bgrbgr(image);

	Mat in(image->ry, image->rx, CV_16UC3, bgrbgr);
	Mat out(image->ry, image->rx, CV_16UC3);
	blur(in, out, Size(3, 3));
	out = in - out;

	Mat channel[3];
	split(out, channel);

	memcpy(image->data, channel[2].data, ndata * sizeof(WORD));
	if (image->naxes[2] == 3) {
		memcpy(image->data + ndata, channel[1].data, ndata * sizeof(WORD));
		memcpy(image->data + ndata * 2, channel[0].data, ndata * sizeof(WORD));
	}

	if (image->naxes[2] == 1) {
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data;
		image->pdata[BLAYER] = image->data;
	} else {
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data + ndata;
		image->pdata[BLAYER] = image->data + ndata * 2;
	}
	image->rx = out.cols;
	image->ry = out.rows;
	image->naxes[0] = image->rx;
	image->naxes[1] = image->ry;

	/* free data */
	delete[] bgrbgr;
	in = Mat();
	out = Mat();

	return 0;
}

#endif
