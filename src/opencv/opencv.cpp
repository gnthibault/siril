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
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "opencv2/core/version.hpp"
#if CV_MAJOR_VERSION == 2
#include "opencv/findHomography/calib3d.hpp"
#else
#if CV_MAJOR_VERSION == 4
#define CV_RANSAC FM_RANSAC
#endif
#include <opencv2/calib3d.hpp>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "registration/matching/misc.h"
#include "registration/matching/atpmatch.h"
#include "opencv.h"
#include "opencv/ecc/ecc.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "algos/statistics.h"
#ifdef __cplusplus
}
#endif


#define defaultRANSACReprojThreshold 3

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

	in.release();
	out.release();
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

	WORD *newdata = (WORD*) realloc(image->data,
			toX * toY * sizeof(WORD) * image->naxes[2]);
	if (!newdata) {
		free(newdata);
		return 1;
	}
	image->rx = toX;
	image->naxes[0] = toX;
	image->ry = toY;
	image->naxes[1] = toY;

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
	in.release();
	out.release();
	channel[0].release();
	channel[1].release();
	channel[2].release();
	invalidate_stats_from_fit(image);
	return 0;
}

int cvTranslateImage(fits *image, point shift, int interpolation) {
	assert(image->data);
	assert(image->rx);
	assert(image->ry);

	int ndata = image->rx * image->ry;

	WORD *bgrbgr = fits_to_bgrbgr(image);

	Mat in(image->ry, image->rx, CV_16UC3, bgrbgr);
	Mat out(image->ry, image->rx, CV_16UC3);

	Mat M = (Mat_<double>(2,3) << 1, 0, shift.x, 0, 1, shift.y);

	warpAffine(in, out, M, in.size(), interpolation);

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
	invalidate_stats_from_fit(image);

	/* free data */
	delete[] bgrbgr;
	in.release();
	out.release();
	channel[0].release();
	channel[1].release();
	channel[2].release();
	return 0;
}

/* Rotate an image with the angle "angle" */
int cvRotateImage(fits *image, point center, double angle, int interpolation, int cropped) {
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
		Point2f pt(center.x, center.y);
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
	invalidate_stats_from_fit(image);

	/* free data */
	delete[] bgrbgr;
	in.release();
	out.release();
	channel[0].release();
	channel[1].release();
	channel[2].release();
	return 0;
}

static void convert_H_to_MatH(Homography *from, Mat &to) {
	to.at<double>(0, 0) = from->h00;
	to.at<double>(0, 1) = from->h01;
	to.at<double>(0, 2) = from->h02;
	to.at<double>(1, 0) = from->h10;
	to.at<double>(1, 1) = from->h11;
	to.at<double>(1, 2) = from->h12;
	to.at<double>(2, 0) = from->h20;
	to.at<double>(2, 1) = from->h21;
	to.at<double>(2, 2) = from->h22;
}

static void convert_MatH_to_H(Mat from, Homography *to) {
	to->h00 = from.at<double>(0, 0);
	to->h01 = from.at<double>(0, 1);
	to->h02 = from.at<double>(0, 2);
	to->h10 = from.at<double>(1, 0);
	to->h11 = from.at<double>(1, 1);
	to->h12 = from.at<double>(1, 2);
	to->h20 = from.at<double>(2, 0);
	to->h21 = from.at<double>(2, 1);
    to->h22 = from.at<double>(2, 2);
}

unsigned char *cvCalculH(s_star *star_array_img,
		struct s_star *star_array_ref, int n, Homography *Hom) {

	std::vector<Point2f> ref;
	std::vector<Point2f> img;
	Mat mask;
	unsigned char *ret;
	int i;

	/* build vectors with lists of stars. */
	for (i = 0; i < n; i++) {
		ref.push_back(Point2f(star_array_ref[i].x, star_array_ref[i].y));
		img.push_back(Point2f(star_array_img[i].x, star_array_img[i].y));
	}

	Mat H = findHomography(img, ref, CV_RANSAC, defaultRANSACReprojThreshold, mask);
	if (countNonZero(H) < 1) {
		return NULL;
	}
	Hom->Inliers = countNonZero(mask);
	ret = (unsigned char *) malloc(n * sizeof(unsigned char));
	for (i = 0; i < n; i++) {
		ret[i] = mask.at<uchar>(i);
	}

	convert_MatH_to_H(H, Hom);

	mask.release();
	H.release();
	return ret;
}

/**
 * Apply the upscale to H. Same value of x and y.
 * @param H1
 * @param scale
 * @return 0
 */
int cvApplyScaleToH(Homography *H1, double s) {
	Mat H = Mat::eye(3, 3, CV_64FC1);
	Mat S = Mat::eye(3, 3, CV_64FC1);

	convert_H_to_MatH(H1, H);

	/* we define Scale Matrix S
	 *
	 *     s   0   0
	 * S = 0   s   0
	 *     0   0   1
	 *
	 */
	S.at<double>(0,0) = s;
	S.at<double>(1,1) = s;

	/* We apply transform */
	Mat result = S * H * S.inv();

	convert_MatH_to_H(result, H1);

	H.release();
	result.release();
	return 0;
}

int cvTransformImage(fits *image, point ref, Homography Hom, int interpolation) {
	assert(image->data);
	assert(image->rx);
	assert(image->ry);

	int ndata = ref.x * ref.y;
	WORD *newdata;
	WORD *bgrbgr = fits_to_bgrbgr(image);

	Mat in(image->ry, image->rx, CV_16UC3, bgrbgr);
	Mat out(ref.y, ref.x, CV_16UC3);
	Mat H = Mat::eye(3, 3, CV_64FC1);

	convert_H_to_MatH(&Hom, H);

	warpPerspective(in, out, H, Size(ref.x, ref.y), interpolation);

	Mat channel[3];
	split(out, channel);

	if (image->ry != ref.y || image->rx != ref.x) {
		newdata = (WORD*) realloc(image->data,
				ref.x * ref.y * sizeof(WORD) * image->naxes[2]);
		if (!newdata) {
			free(newdata);
			return 1;
		}
		image->data = newdata;
	}

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
	invalidate_stats_from_fit(image);

	/* free data */
	delete[] bgrbgr;
	in.release();
	out.release();
	channel[0].release();
	channel[1].release();
	channel[2].release();
	H.release();
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
	invalidate_stats_from_fit(image);
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
	invalidate_stats_from_fit(image);

	/* free data */
	delete[] bgrbgr;
	channel[0].release();
	channel[1].release();
	channel[2].release();
	in.release();
	out.release();

	return 0;
}

static Mat RLTikh_deconvolution(Mat observed, Mat psf, double mu, int iterations) {

	Mat deconv = observed.clone();

	// Iterate
	for (int i = 0; i < iterations; i++) {
		set_progress_bar_data(NULL, (double) i / iterations);

		// Temporary matrix
		Mat ratio;
		sepFilter2D(deconv, ratio, deconv.depth(), psf, psf, Point(-1, -1), 0,
				BORDER_REFLECT);

		divide(observed, ratio, ratio);

		sepFilter2D(ratio, ratio, ratio.depth(), psf, psf, Point(-1, -1), 0,
				BORDER_REFLECT);

		// TV Regularization
		Mat denom;
		Laplacian(deconv, denom, deconv.depth(), 1, 1, 0, BORDER_REFLECT);
		denom = 1.0 - 2.0 * mu * denom;
		divide(ratio, denom, ratio);

		// Apply iteration on the estimate
		multiply(deconv, ratio, deconv);
	}

	return deconv;
}

#define KERNEL_SIZE_FACTOR 6

int cvLucyRichardson(fits *image, double sigma, int iterations) {
	assert(image->data);
	assert(image->rx);
	assert(image->ry);

	int ndata = image->rx * image->ry;
	int ksize = (int)((KERNEL_SIZE_FACTOR * sigma) + 0.5);
	ksize = ksize % 2 != 0 ? ksize : ksize + 1;

	WORD *bgrbgr = fits_to_bgrbgr(image);

	Mat in(image->ry, image->rx, CV_16UC3, bgrbgr);
	Mat out(image->ry, image->rx, CV_16UC3);
	set_progress_bar_data(_("Deconvolution..."), PROGRESS_NONE);

	// From here on, use 64-bit floats
	// Convert original_image to float
	Mat float_image;
	in.convertTo(float_image, CV_64FC3);
	float_image *= 1.0 / 65535.0;

	Mat psf = getGaussianKernel(ksize, sigma);

	double mu = 0.01;

	Mat estimation = RLTikh_deconvolution(float_image, psf, mu, iterations);
	estimation *= 65535.0;

	estimation.convertTo(out, CV_16UC3);

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

	set_progress_bar_data(_("Deconvoltuion applied"), PROGRESS_DONE);
	/* free data */
	delete[] bgrbgr;
	in.release();
	out.release();
	channel[0].release();
	channel[1].release();
	channel[2].release();
	invalidate_stats_from_fit(image);

	return 0;

}
