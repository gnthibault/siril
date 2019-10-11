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

static BYTE *fits8_to_bgrbgr(fits *image) {
	int ndata = image->rx * image->ry;
	BYTE *bgrbgr = new BYTE[ndata * 3];
	for (int i = 0, j = 0; i < ndata * 3; i += 3, j++) {
		bgrbgr[i + 0] = (BYTE)image->pdata[BLAYER][j];
		bgrbgr[i + 1] = (BYTE)image->pdata[GLAYER][j];
		bgrbgr[i + 2] = (BYTE)image->pdata[RLAYER][j];
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

	// preparing data
	Mat in, out;
	WORD *bgrbgr = NULL;

	if (image->naxes[2] == 1) {
		in = Mat(image->ry, image->rx, CV_16UC1, image->data);
		out = Mat(toY, toX, CV_16UC1);
	}
	else if (image->naxes[2] == 3) {
		bgrbgr = fits_to_bgrbgr(image);
		in = Mat(image->ry, image->rx, CV_16UC3, bgrbgr);
		out = Mat(toY, toX, CV_16UC3);
	}
	else {
		siril_log_message("Image resizing is not supported for images with %d channels\n", image->naxes[2]);
		return -1;
	}

	// OpenCV function
	resize(in, out, out.size(), 0, 0, interpolation);

	// saving result
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

	unsigned int dataSize = toX * toY;
	if (image->naxes[2] == 1) {
		memcpy(image->data, out.data, dataSize * sizeof(WORD));
		image->pdata[0] = image->data;
		image->pdata[1] = image->data;
		image->pdata[2] = image->data;
	}
	else {
		std::vector<Mat> channel(3);
		split(out, channel);
		memcpy(image->data, channel[2].data, dataSize * sizeof(WORD));
		memcpy(image->data + dataSize, channel[1].data, dataSize * sizeof(WORD));
		memcpy(image->data + dataSize * 2, channel[0].data, dataSize * sizeof(WORD));

		image->pdata[0] = image->data;
		image->pdata[1] = image->data + dataSize;
		image->pdata[2] = image->data + dataSize * 2;

		channel[0].release();
		channel[1].release();
		channel[2].release();
		delete[] bgrbgr;
	}
	in.release();
	out.release();
	invalidate_stats_from_fit(image);
	return 0;
}

int cvTranslateImage(fits *image, point shift, int interpolation) {
	assert(image->data);
	assert(image->rx);
	assert(image->ry);

	int ndata = image->rx * image->ry;

	WORD *bgrbgr = NULL;
	Mat in, out;

	if (image->naxes[2] == 1) {
		in = Mat(image->ry, image->rx, CV_16UC1, image->data);
		out = Mat(image->ry, image->rx, CV_16UC1);
	}
	else if (image->naxes[2] == 3) {
		bgrbgr = fits_to_bgrbgr(image);
		in = Mat(image->ry, image->rx, CV_16UC3, bgrbgr);
		out = Mat(image->ry, image->rx, CV_16UC3);
	}
	else {
		siril_log_message(_("Translate is not supported for images with %d channels\n"), image->naxes[2]);
		return -1;
	}

	Mat M = (Mat_<double>(2,3) << 1, 0, shift.x, 0, 1, shift.y);

	warpAffine(in, out, M, in.size(), interpolation);

	std::vector<Mat> channel(3);
	split(out, channel);

	memcpy(image->data, channel[2].data, ndata * sizeof(WORD));
	if (image->naxes[2] == 3) {
		memcpy(image->data + ndata, channel[1].data, ndata * sizeof(WORD));
		memcpy(image->data + ndata * 2, channel[0].data, ndata * sizeof(WORD));
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data + ndata;
		image->pdata[BLAYER] = image->data + ndata * 2;
	} else {
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data;
		image->pdata[BLAYER] = image->data;
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

	WORD *bgrbgr = NULL;
	Mat in, out;

	if (image->naxes[2] == 1) {
		in = Mat(image->ry, image->rx, CV_16UC1, image->data);
		out = Mat(image->ry, image->rx, CV_16UC1);
	}
	else if (image->naxes[2] == 3) {
		bgrbgr = fits_to_bgrbgr(image);
		in = Mat(image->ry, image->rx, CV_16UC3, bgrbgr);
		out = Mat(image->ry, image->rx, CV_16UC3);
	}
	else {
		siril_log_message(_("Rotate is not supported for images with %d channels\n"), image->naxes[2]);
		return -1;
	}

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
	std::vector<Mat> channel(3);
	split(out, channel);

	memcpy(image->data, channel[2].data, ndata * sizeof(WORD));
	if (image->naxes[2] == 3) {
		memcpy(image->data + ndata, channel[1].data, ndata * sizeof(WORD));
		memcpy(image->data + ndata * 2, channel[0].data, ndata * sizeof(WORD));
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data + ndata;
		image->pdata[BLAYER] = image->data + ndata * 2;
	} else {
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data;
		image->pdata[BLAYER] = image->data;
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
	unsigned char *ret = NULL;
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
	if (n > 0) {
		ret = (unsigned char *) malloc(n * sizeof(unsigned char));
		for (i = 0; i < n; i++) {
			ret[i] = mask.at<uchar>(i);
		}
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

// transform an image using the homography.
int cvTransformImage(fits *image, long width, long height, Homography Hom, int interpolation) {
	assert(image->data);
	assert(image->rx);
	assert(image->ry);

	// preparing data
	Mat in, out;
	WORD *bgrbgr = NULL;
	WORD *newdata;

	if (image->naxes[2] == 1) {
		in = Mat(image->ry, image->rx, CV_16UC1, image->data);
		out = Mat(height, width, CV_16UC1);
	}
	else if (image->naxes[2] == 3) {
		bgrbgr = fits_to_bgrbgr(image);
		in = Mat(image->ry, image->rx, CV_16UC3, bgrbgr);
		out = Mat(height, width, CV_16UC3);
	}
	else {
		siril_log_message(_("Transformation is not supported for images with %d channels\n"), image->naxes[2]);
		return -1;
	}

	Mat H = Mat::eye(3, 3, CV_64FC1);
	convert_H_to_MatH(&Hom, H);

	// OpenCV function
	warpPerspective(in, out, H, Size(width, height), interpolation);

	// saving result
	long ndata = height * width;

	if (image->ry != height || image->rx != width) {
		newdata = (WORD*) realloc(image->data, ndata * sizeof(WORD) * image->naxes[2]);
		if (!newdata) {
			PRINT_ALLOC_ERR;
			return 1;
		}
		image->data = newdata;
	}

	if (image->naxes[2] == 1) {
		memcpy(image->data, out.data, ndata * sizeof(WORD));
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data;
		image->pdata[BLAYER] = image->data;
	} else {
		std::vector<Mat> channel(3);
		split(out, channel);
		memcpy(image->data, channel[2].data, ndata * sizeof(WORD));
		memcpy(image->data + ndata, channel[1].data, ndata * sizeof(WORD));
		memcpy(image->data + ndata * 2, channel[0].data, ndata * sizeof(WORD));

		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data + ndata;
		image->pdata[BLAYER] = image->data + ndata * 2;


		channel[0].release();
		channel[1].release();
		channel[2].release();
		delete[] bgrbgr;
	}
	image->rx = image->naxes[0] = out.cols;
	image->ry = image->naxes[1] = out.rows;
	H.release();
	in.release();
	out.release();
	invalidate_stats_from_fit(image);
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

	WORD *bgrbgr = NULL;
	Mat in, out;

	if (image->naxes[2] == 1) {
		in = Mat(image->ry, image->rx, CV_16UC1, image->data);
		out = Mat(image->ry, image->rx, CV_16UC1);
	}
	else if (image->naxes[2] == 3) {
		bgrbgr = fits_to_bgrbgr(image);
		in = Mat(image->ry, image->rx, CV_16UC3, bgrbgr);
		out = Mat(image->ry, image->rx, CV_16UC3);
	}
	else {
		siril_log_message(_("Deconvolution is not supported for images with %d channels\n"), image->naxes[2]);
		return -1;
	}

	blur(in, out, Size(3, 3));
	out = in - out;

	std::vector<Mat> channel(3);
	split(out, channel);

	if (image->naxes[2] == 3) {
		memcpy(image->data, channel[2].data, ndata * sizeof(WORD));
		memcpy(image->data + ndata, channel[1].data, ndata * sizeof(WORD));
		memcpy(image->data + ndata * 2, channel[0].data, ndata * sizeof(WORD));
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data + ndata;
		image->pdata[BLAYER] = image->data + ndata * 2;
	} else {
		memcpy(image->data, channel[0].data, ndata * sizeof(WORD));
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data;
		image->pdata[BLAYER] = image->data;
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
	Mat in, out;

	WORD *bgrbgr = NULL;

	if (image->naxes[2] == 1) {
		in = Mat(image->ry, image->rx, CV_16UC1, image->data);
		out = Mat(image->ry, image->rx, CV_16UC1);
	}
	else if (image->naxes[2] == 3) {
		bgrbgr = fits_to_bgrbgr(image);
		in = Mat(image->ry, image->rx, CV_16UC3, bgrbgr);
		out = Mat(image->ry, image->rx, CV_16UC3);
	}
	else {
		siril_log_message(_("Deconvolution is not supported for images with %d channels\n"), image->naxes[2]);
		return -1;
	}

	set_progress_bar_data(_("Deconvolution..."), PROGRESS_NONE);

	// From here on, use 64-bit floats
	// Convert original_image to float
	Mat float_image;
	in.convertTo(float_image, CV_64FC3);
	float_image *= 1.0 / USHRT_MAX_DOUBLE;

	Mat psf = getGaussianKernel(ksize, sigma);

	double mu = 0.01;

	Mat estimation = RLTikh_deconvolution(float_image, psf, mu, iterations);
	estimation *= USHRT_MAX_DOUBLE;

	estimation.convertTo(out, CV_16UC3);

	std::vector<Mat> channel(3);
	split(out, channel);

	if (image->naxes[2] == 3) {
		memcpy(image->data, channel[2].data, ndata * sizeof(WORD));
		memcpy(image->data + ndata, channel[1].data, ndata * sizeof(WORD));
		memcpy(image->data + ndata * 2, channel[0].data, ndata * sizeof(WORD));
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data + ndata;
		image->pdata[BLAYER] = image->data + ndata * 2;
	} else {
		memcpy(image->data, channel[0].data, ndata * sizeof(WORD));
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data;
		image->pdata[BLAYER] = image->data;
	}

	image->rx = out.cols;
	image->ry = out.rows;
	image->naxes[0] = image->rx;
	image->naxes[1] = image->ry;

	set_progress_bar_data(_("Deconvolution applied"), PROGRESS_DONE);
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

/* Work on grey images. If image is in RGB it must be first converted
 * in CieLAB. Then, only the first channel is applied
 */
int cvClahe(fits *image, double clip_limit, int size) {
	assert(image->data);
	assert(image->rx);
	assert(image->ry);

	int ndata = image->rx * image->ry;

	// preparing data
	Mat in, out;

	Ptr<CLAHE> clahe = createCLAHE();
	clahe->setClipLimit(clip_limit);
	clahe->setTilesGridSize(Size(size, size));

	if (image->naxes[2] == 3) {
		Mat lab_image;
		std::vector<Mat> lab_planes(3);
		BYTE *bgrbgr8;
		WORD *bgrbgr;

		switch (image->bitpix) {
		case BYTE_IMG:
			bgrbgr8 = fits8_to_bgrbgr(image);
			in = Mat(image->ry, image->rx, CV_8UC3, bgrbgr8);
			out = Mat();
			// convert the RGB color image to Lab
			cvtColor(in, lab_image, COLOR_BGR2Lab);

			// Extract the L channel
			split(lab_image, lab_planes); // now we have the L image in lab_planes[0]

			// apply the CLAHE algorithm to the L channel
			clahe->apply(lab_planes[0], lab_planes[0]);

			// Merge the color planes back into an Lab image
			merge(lab_planes, lab_image);

			// convert back to RGB
			cvtColor(lab_image, out, COLOR_Lab2BGR);
			out.convertTo(out, CV_16UC3, 1.0);

			delete[] bgrbgr8;

			break;
		default:
		case USHORT_IMG:
			bgrbgr = fits_to_bgrbgr(image);
			in = Mat(image->ry, image->rx, CV_16UC3, bgrbgr);
			in.convertTo(in, CV_32F, 1.0 / USHRT_MAX_DOUBLE);
			out = Mat();

			// convert the RGB color image to Lab
			cvtColor(in, lab_image, COLOR_BGR2Lab);

			// Extract the L channel
			split(lab_image, lab_planes); // now we have the L image in lab_planes[0]

			// apply the CLAHE algorithm to the L channel (does not work with 32F images)
			lab_planes[0].convertTo(lab_planes[0], CV_16U,	USHRT_MAX_DOUBLE / 100.0);
			clahe->apply(lab_planes[0], lab_planes[0]);
			lab_planes[0].convertTo(lab_planes[0], CV_32F, 100.0 / USHRT_MAX_DOUBLE);

			// Merge the color planes back into an Lab image
			merge(lab_planes, lab_image);

			// convert back to RGB
			cvtColor(lab_image, out, COLOR_Lab2BGR);
			out.convertTo(out, CV_16UC3, USHRT_MAX_DOUBLE);

			delete[] bgrbgr;
		}

	} else {
		in = Mat(image->ry, image->rx, CV_16UC1, image->data);
		out = Mat();
		switch (image->bitpix) {
		case BYTE_IMG:
			in.convertTo(in, CV_8U, 1.0);
			clahe->apply(in, out);
			out.convertTo(out, CV_16UC3, 1.0);
			break;
		default:
		case USHORT_IMG:
			clahe->apply(in, out);
		}
	}

	std::vector<Mat> channel(3);
	split(out, channel);

	if (image->naxes[2] == 3) {
		memcpy(image->data, channel[2].data, ndata * sizeof(WORD));
		memcpy(image->data + ndata, channel[1].data, ndata * sizeof(WORD));
		memcpy(image->data + ndata * 2, channel[0].data, ndata * sizeof(WORD));
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data + ndata;
		image->pdata[BLAYER] = image->data + ndata * 2;
	} else {
		memcpy(image->data, channel[0].data, ndata * sizeof(WORD));
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data;
		image->pdata[BLAYER] = image->data;
	}

	image->rx = out.cols;
	image->ry = out.rows;
	image->naxes[0] = image->rx;
	image->naxes[1] = image->ry;

	/* free data */
	in.release();
	out.release();
	channel[0].release();
	channel[1].release();
	channel[2].release();
	invalidate_stats_from_fit(image);

	return 0;
}
