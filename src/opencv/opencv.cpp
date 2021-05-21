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

static WORD *fits_to_bgrbgr_ushort(fits *image) {
	size_t ndata = image->rx * image->ry * 3;
	WORD *bgrbgr = new WORD[ndata];
	for (size_t i = 0, j = 0; i < ndata; i += 3, j++) {
		bgrbgr[i + 0] = image->pdata[BLAYER][j];
		bgrbgr[i + 1] = image->pdata[GLAYER][j];
		bgrbgr[i + 2] = image->pdata[RLAYER][j];
	}
	return bgrbgr;
}

static float *fits_to_bgrbgr_float(fits *image) {
	size_t ndata = image->rx * image->ry * 3;
	float *bgrbgr = new float[ndata];
	for (size_t i = 0, j = 0; i < ndata; i += 3, j++) {
		bgrbgr[i + 0] = image->fpdata[BLAYER][j];
		bgrbgr[i + 1] = image->fpdata[GLAYER][j];
		bgrbgr[i + 2] = image->fpdata[RLAYER][j];
	}
	return bgrbgr;
}

static BYTE *fits8_to_bgrbgr(fits *image) {
	size_t ndata = image->rx * image->ry * 3;
	BYTE *bgrbgr = new BYTE[ndata];
	for (size_t i = 0, j = 0; i < ndata; i += 3, j++) {
		bgrbgr[i + 0] = (BYTE)image->pdata[BLAYER][j];
		bgrbgr[i + 1] = (BYTE)image->pdata[GLAYER][j];
		bgrbgr[i + 2] = (BYTE)image->pdata[RLAYER][j];
	}
	return bgrbgr;
}

int cvResizeGaussian_uchar(uint8_t *dataIn, int rx, int ry, uint8_t *dataOut,
		int toX, int toY, int chan, int interpolation) {
	int mode = (chan == 1 ? CV_8UC1 : CV_8UC3);

	Mat in(ry, rx, mode, dataIn);
	Mat out(toY, toX, mode);

	resize(in, out, out.size(), 0, 0, interpolation);

	size_t ndata = toX * toY * chan;
	for (size_t i = 0; i < ndata; i++)
		dataOut[i] = (uint8_t) out.data[i];

	in.release();
	out.release();
	return 0;
}

/* resizes image to the sizes toX * toY, and stores it back in image */
static int cvResizeGaussian_ushort(fits *image, int toX, int toY, int interpolation) {
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
		bgrbgr = fits_to_bgrbgr_ushort(image);
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
	size_t nbpixels = toX * toY;
	size_t dataSize = nbpixels * sizeof(WORD);
	WORD *newdata = (WORD*) realloc(image->data, dataSize * image->naxes[2]);
	if (!newdata) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	image->rx = toX;
	image->naxes[0] = toX;
	image->ry = toY;
	image->naxes[1] = toY;
	image->data = newdata;

	if (image->naxes[2] == 1) {
		memcpy(image->data, out.data, dataSize);
		image->pdata[0] = image->data;
		image->pdata[1] = image->data;
		image->pdata[2] = image->data;
	}
	else {
		std::vector<Mat> channel(3);
		split(out, channel);
		memcpy(image->data, channel[2].data, dataSize);
		memcpy(image->data + nbpixels, channel[1].data, dataSize);
		memcpy(image->data + nbpixels * 2, channel[0].data, dataSize);

		image->pdata[0] = image->data;
		image->pdata[1] = image->data + nbpixels;
		image->pdata[2] = image->data + nbpixels * 2;

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

static int cvResizeGaussian_float(fits *image, int toX, int toY, int interpolation) {
	assert(image->fdata);
	assert(image->rx);

	// preparing data
	Mat in, out;
	float *bgrbgr = NULL;

	if (image->naxes[2] == 1) {
		in = Mat(image->ry, image->rx, CV_32FC1, image->fdata);
		out = Mat(toY, toX, CV_32FC1);
	}
	else if (image->naxes[2] == 3) {
		bgrbgr = fits_to_bgrbgr_float(image);
		in = Mat(image->ry, image->rx, CV_32FC3, bgrbgr);
		out = Mat(toY, toX, CV_32FC3);
	}
	else {
		siril_log_message("Image resizing is not supported for images with %d channels\n", image->naxes[2]);
		return -1;
	}

	// OpenCV function
	resize(in, out, out.size(), 0, 0, interpolation);

	// saving result
	size_t nbpixels = toX * toY;
	size_t dataSize = nbpixels * sizeof(float);
	float *newdata = (float*) realloc(image->fdata, dataSize * image->naxes[2]);
	if (!newdata) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	image->rx = toX;
	image->naxes[0] = toX;
	image->ry = toY;
	image->naxes[1] = toY;
	image->fdata = newdata;

	if (image->naxes[2] == 1) {
		memcpy(image->fdata, out.data, dataSize);
		image->fpdata[0] = image->fdata;
		image->fpdata[1] = image->fdata;
		image->fpdata[2] = image->fdata;
	}
	else {
		std::vector<Mat> channel(3);
		split(out, channel);
		memcpy(image->fdata, channel[2].data, dataSize);
		memcpy(image->fdata + nbpixels, channel[1].data, dataSize);
		memcpy(image->fdata + nbpixels * 2, channel[0].data, dataSize);

		image->fpdata[0] = image->fdata;
		image->fpdata[1] = image->fdata + nbpixels;
		image->fpdata[2] = image->fdata + nbpixels * 2;

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

/* resizes image to the sizes toX * toY, and stores it back in image */
int cvResizeGaussian(fits *image, int toX, int toY, int interpolation) {
	if (image->type == DATA_USHORT)
		return cvResizeGaussian_ushort(image, toX, toY, interpolation);
	if (image->type == DATA_FLOAT)
		return cvResizeGaussian_float(image, toX, toY, interpolation);
	return -1;
}

/* Rotate an image with the angle "angle" */
int cvRotateImage_ushort(fits *image, point center, double angle, int interpolation, int cropped) {
	assert(image->data);
	assert(image->rx);
	assert(image->ry);

	size_t ndata;
	WORD *bgrbgr = NULL;
	Mat in, out;

	if (image->naxes[2] == 1) {
		in = Mat(image->ry, image->rx, CV_16UC1, image->data);
		out = Mat(image->ry, image->rx, CV_16UC1);
	}
	else if (image->naxes[2] == 3) {
		bgrbgr = fits_to_bgrbgr_ushort(image);
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
		ndata = image->naxes[0] * image->naxes[1];
	} else {
		Point2f pt(center.x, center.y);
		Mat r = getRotationMatrix2D(pt, angle, 1.0);
		if (cropped == 1) {
			warpAffine(in, out, r, in.size(), interpolation);
			ndata = image->naxes[0] * image->naxes[1];
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
				PRINT_ALLOC_ERR;
				return 1;
			}
			image->data = newdata;
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

int cvRotateImage_float(fits *image, point center, double angle, int interpolation, int cropped) {
	assert(image->fdata);
	assert(image->rx);
	assert(image->ry);

	size_t ndata;
	float *bgrbgr = NULL;
	Mat in, out;

	if (image->naxes[2] == 1) {
		in = Mat(image->ry, image->rx, CV_32FC1, image->fdata);
		out = Mat(image->ry, image->rx, CV_32FC1);
	}
	else if (image->naxes[2] == 3) {
		bgrbgr = fits_to_bgrbgr_float(image);
		in = Mat(image->ry, image->rx, CV_32FC3, bgrbgr);
		out = Mat(image->ry, image->rx, CV_32FC3);
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
		ndata = image->naxes[0] * image->naxes[1];
	} else {
		Point2f pt(center.x, center.y);
		Mat r = getRotationMatrix2D(pt, angle, 1.0);
		if (cropped == 1) {
			warpAffine(in, out, r, in.size(), interpolation);
			ndata = image->naxes[0] * image->naxes[1];
		} else {

			// determine bounding rectangle
			Rect frame = RotatedRect(pt, in.size(), angle).boundingRect();
			// adjust transformation matrix
			r.at<double>(0, 2) += frame.width / 2.0 - pt.x;
			r.at<double>(1, 2) += frame.height / 2.0 - pt.y;

			warpAffine(in, out, r, frame.size(), interpolation);
			ndata = out.cols * out.rows;
			float *newdata = (float*) realloc(image->fdata, ndata * image->naxes[2] * sizeof(float));
			if (!newdata) {
				PRINT_ALLOC_ERR;
				return 1;
			}
			image->fdata = newdata;
		}
	}
	std::vector<Mat> channel(3);
	split(out, channel);

	if (image->naxes[2] == 3) {
		memcpy(image->fdata, channel[2].data, ndata * sizeof(float));
		memcpy(image->fdata + ndata, channel[1].data, ndata * sizeof(float));
		memcpy(image->fdata + ndata * 2, channel[0].data, ndata * sizeof(float));
		image->fpdata[RLAYER] = image->fdata;
		image->fpdata[GLAYER] = image->fdata + ndata;
		image->fpdata[BLAYER] = image->fdata + ndata * 2;
	} else {
		memcpy(image->fdata, channel[0].data, ndata * sizeof(float));
		image->fpdata[RLAYER] = image->fdata;
		image->fpdata[GLAYER] = image->fdata;
		image->fpdata[BLAYER] = image->fdata;
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

int cvRotateImage(fits *image, point center, double angle, int interpolation, int cropped) {
	if (image->type == DATA_USHORT)
		return cvRotateImage_ushort(image, center, angle, interpolation, cropped);
	if (image->type == DATA_FLOAT)
		return cvRotateImage_float(image, center, angle, interpolation, cropped);
	return -1;
}

int cvAffineTransformation(fits *image, pointf *refpoints, pointf *curpoints, int nb_points, gboolean upscale2x, int interpolation) {
	// see https://docs.opencv.org/3.4/d4/d61/tutorial_warp_affine.html
	std::vector<Point2f> ref;
	std::vector<Point2f> cur;

	/* build vectors with lists of 3 stars. */
	for (int i = 0; i < nb_points; i++) {
		ref.push_back(Point2f(refpoints[i].x, image->ry - refpoints[i].y - 1));
		if (upscale2x)
			cur.push_back(Point2f(curpoints[i].x * 0.5f, (image->ry - curpoints[i].y - 1) * 0.5f));
		else cur.push_back(Point2f(curpoints[i].x, image->ry - curpoints[i].y - 1));
	}

	Mat m = estimateAffinePartial2D(cur, ref);
	//std::cout << m << std::endl;

	/* test that m is not a zero matrix */
	if (countNonZero(m) < 1) {
		siril_log_color_message(_("Singular Matrix. Cannot compute Affine Transformation.\n"), "red");
		return -1;
	}

	/* convert image */
	Mat in, out;
	WORD *bgrbgr_ushort = NULL;
	float *bgrbgr_float = NULL;
	int target_rx = image->rx, target_ry = image->ry;
	if (upscale2x) {
		target_rx *= 2;
		target_ry *= 2;
	}

	if (image->naxes[2] != 1 && image->naxes[2] != 3) {
		siril_log_message(_("Transformation is not supported for images with %ld channels\n"), image->naxes[2]);
		return -1;
	}
	if (image->type == DATA_USHORT) {
		if (image->naxes[2] == 1) {
			in = Mat(image->ry, image->rx, CV_16UC1, image->data);
			out = Mat(target_ry, target_rx, CV_16UC1);
		}
		else if (image->naxes[2] == 3) {
			bgrbgr_ushort = fits_to_bgrbgr_ushort(image);
			in = Mat(image->ry, image->rx, CV_16UC3, bgrbgr_ushort);
			out = Mat(target_ry, target_rx, CV_16UC3);
		}
	}
	else if (image->type == DATA_FLOAT) {
		if (image->naxes[2] == 1) {
			in = Mat(image->ry, image->rx, CV_32FC1, image->fdata);
			out = Mat(target_ry, target_rx, CV_32FC1);
		}
		else if (image->naxes[2] == 3) {
			bgrbgr_float = fits_to_bgrbgr_float(image);
			in = Mat(image->ry, image->rx, CV_32FC3, bgrbgr_float);
			out = Mat(target_ry, target_rx, CV_32FC3);
		}
	}
	else return -1;

	warpAffine(in, out, m, out.size(), interpolation, BORDER_TRANSPARENT);

	/* store result */
	std::vector<Mat> channel(3);
	split(out, channel);

	size_t ndata = target_rx * target_ry;
	if (image->type == DATA_USHORT) {
		size_t data_size = ndata * sizeof(WORD);
		WORD *newdata = (WORD*) realloc(image->data, data_size * image->naxes[2]);
		if (!newdata) {
			PRINT_ALLOC_ERR;
			return 1;
		}
		image->data = newdata;

		if (image->naxes[2] == 3) {
			memcpy(image->data, channel[2].data, data_size);
			memcpy(image->data + ndata, channel[1].data, data_size);
			memcpy(image->data + ndata * 2, channel[0].data, data_size);
			image->pdata[RLAYER] = image->data;
			image->pdata[GLAYER] = image->data + ndata;
			image->pdata[BLAYER] = image->data + ndata * 2;
		} else {
			memcpy(image->data, channel[0].data, data_size);
			image->pdata[RLAYER] = image->data;
			image->pdata[GLAYER] = image->data;
			image->pdata[BLAYER] = image->data;
		}
	} else {
		size_t data_size = ndata * sizeof(float);
		float *newdata = (float *) realloc(image->fdata, data_size * image->naxes[2]);
		if (!newdata) {
			PRINT_ALLOC_ERR;
			return 1;
		}
		image->fdata = newdata;


		if (image->naxes[2] == 3) {
			memcpy(image->fdata, channel[2].data, data_size);
			memcpy(image->fdata + ndata, channel[1].data, data_size);
			memcpy(image->fdata + ndata * 2, channel[0].data, data_size);
			image->fpdata[RLAYER] = image->fdata;
			image->fpdata[GLAYER] = image->fdata + ndata;
			image->fpdata[BLAYER] = image->fdata + ndata * 2;
		} else {
			memcpy(image->fdata, channel[0].data, data_size);
			image->fpdata[RLAYER] = image->fdata;
			image->fpdata[GLAYER] = image->fdata;
			image->fpdata[BLAYER] = image->fdata;
		}
	}
	image->rx = out.cols;
	image->ry = out.rows;
	image->naxes[0] = image->rx;
	image->naxes[1] = image->ry;
	invalidate_stats_from_fit(image);

	/* free data */
	if (bgrbgr_ushort) delete[] bgrbgr_ushort;
	if (bgrbgr_float) delete[] bgrbgr_float;
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

static int cvTransformImage_ushort(fits *image, long width, long height, Homography Hom, int interpolation) {
	// preparing data
	Mat in, out;
	WORD *bgrbgr = NULL;
	WORD *newdata;
	assert(image->data);

	if (image->naxes[2] == 1) {
		in = Mat(image->ry, image->rx, CV_16UC1, image->data);
		out = Mat(height, width, CV_16UC1);
	}
	else if (image->naxes[2] == 3) {
		bgrbgr = fits_to_bgrbgr_ushort(image);
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
	warpPerspective(in, out, H, Size(width, height), interpolation, BORDER_TRANSPARENT);

	// saving result
	size_t ndata = height * width;
	if (image->naxes[1] != height || image->naxes[0] != width) {
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

static int cvTransformImage_float(fits *image, long width, long height, Homography Hom, int interpolation) {
	// preparing data
	Mat in, out;
	float *bgrbgr = NULL;
	float *newdata;
	assert(image->fdata);

	if (image->naxes[2] == 1) {
		in = Mat(image->ry, image->rx, CV_32FC1, image->fdata);
		out = Mat(height, width, CV_32FC1);
	}
	else if (image->naxes[2] == 3) {
		bgrbgr = fits_to_bgrbgr_float(image);
		in = Mat(image->ry, image->rx, CV_32FC3, bgrbgr);
		out = Mat(height, width, CV_32FC3);
	}
	else {
		siril_log_message(_("Transformation is not supported for images with %d channels\n"), image->naxes[2]);
		return -1;
	}

	Mat H = Mat::eye(3, 3, CV_64FC1);
	convert_H_to_MatH(&Hom, H);

	// OpenCV function
	warpPerspective(in, out, H, Size(width, height), interpolation, BORDER_TRANSPARENT);

	// saving result
	size_t ndata = height * width;
	if (image->naxes[1] != height || image->naxes[0] != width) {
		newdata = (float*) realloc(image->fdata, ndata * sizeof(float) * image->naxes[2]);
		if (!newdata) {
			PRINT_ALLOC_ERR;
			return 1;
		}
		image->fdata = newdata;
	}

	if (image->naxes[2] == 1) {
		memcpy(image->fdata, out.data, ndata * sizeof(float));
		image->fpdata[RLAYER] = image->fdata;
		image->fpdata[GLAYER] = image->fdata;
		image->fpdata[BLAYER] = image->fdata;
	} else {
		std::vector<Mat> channel(3);
		split(out, channel);
		memcpy(image->fdata, channel[2].data, ndata * sizeof(float));
		memcpy(image->fdata + ndata, channel[1].data, ndata * sizeof(float));
		memcpy(image->fdata + ndata * 2, channel[0].data, ndata * sizeof(float));

		image->fpdata[RLAYER] = image->fdata;
		image->fpdata[GLAYER] = image->fdata + ndata;
		image->fpdata[BLAYER] = image->fdata + ndata * 2;

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

// transform an image using the homography.
int cvTransformImage(fits *image, long width, long height, Homography Hom, int interpolation) {
	assert(image->rx);
	assert(image->ry);
	if (image->type == DATA_USHORT)
		return cvTransformImage_ushort(image, width, height, Hom, interpolation);
	if (image->type == DATA_FLOAT)
		return cvTransformImage_float(image, width, height, Hom, interpolation);
	return -1;
}

static int cvUnsharpFilter_ushort(fits* image, double sigma, double amount) {
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

	size_t nbpixels = image->naxes[0] * image->naxes[1];
	memcpy(image->data, out.data,
			nbpixels * sizeof(WORD) * image->naxes[2]);
	if (image->naxes[2] == 1) {
		image->pdata[0] = image->data;
		image->pdata[1] = image->data;
		image->pdata[2] = image->data;
	} else {
		image->pdata[0] = image->data;
		image->pdata[1] = image->data + nbpixels;
		image->pdata[2] = image->data + nbpixels * 2;
	}
	in.release();
	invalidate_stats_from_fit(image);
	return 0;
}

static int cvUnsharpFilter_float(fits* image, double sigma, double amount) {
	assert(image->fdata);
	assert(image->rx);
	int type = CV_32F;
	if (image->naxes[2] != 1)
		type = CV_32FC3;

	Mat in(image->ry, image->rx, type, image->fdata);
	Mat out, contrast;
	GaussianBlur(in, out, Size(), sigma);
	if (fabs(amount) > 0.0) {
		Mat sharpened = in * (1 + amount) + out * (-amount);
		out = sharpened.clone();
		sharpened.release();
	}

	size_t nbpixels = image->naxes[0] * image->naxes[1];
	memcpy(image->fdata, out.data,
			nbpixels * sizeof(float) * image->naxes[2]);
	if (image->naxes[2] == 1) {
		image->fpdata[0] = image->fdata;
		image->fpdata[1] = image->fdata;
		image->fpdata[2] = image->fdata;
	} else {
		image->fpdata[0] = image->fdata;
		image->fpdata[1] = image->fdata + nbpixels;
		image->fpdata[2] = image->fdata + nbpixels * 2;
	}
	in.release();
	invalidate_stats_from_fit(image);
	return 0;
}

int cvUnsharpFilter(fits* image, double sigma, double amount) {
	if (image->type == DATA_USHORT)
		return cvUnsharpFilter_ushort(image, sigma, amount);
	if (image->type == DATA_FLOAT)
		return cvUnsharpFilter_float(image, sigma, amount);
	return -1;
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

static int cvLucyRichardson_ushort(fits *image, double sigma, int iterations) {
	assert(image->data);
	assert(image->rx);
	assert(image->ry);

	int ksize = (int)((KERNEL_SIZE_FACTOR * sigma) + 0.5);
	ksize = ksize % 2 != 0 ? ksize : ksize + 1;
	Mat in, out;

	WORD *bgrbgr = NULL;

	if (image->naxes[2] == 1) {
		in = Mat(image->ry, image->rx, CV_16UC1, image->data);
		out = Mat(image->ry, image->rx, CV_16UC1);
	}
	else if (image->naxes[2] == 3) {
		bgrbgr = fits_to_bgrbgr_ushort(image);
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
	in.convertTo(float_image, CV_32FC3);
	float_image *= 1.0 / USHRT_MAX_DOUBLE;

	Mat psf = getGaussianKernel(ksize, sigma);

	double mu = 0.01;

	Mat estimation = RLTikh_deconvolution(float_image, psf, mu, iterations);
	estimation *= USHRT_MAX_DOUBLE;

	estimation.convertTo(out, CV_16UC3);

	std::vector<Mat> channel(3);
	split(out, channel);

	size_t nbpixels = image->naxes[0] * image->naxes[1];
	size_t ndata = nbpixels * sizeof(WORD);
	if (image->naxes[2] == 3) {
		memcpy(image->data, channel[2].data, ndata);
		memcpy(image->data + nbpixels, channel[1].data, ndata);
		memcpy(image->data + nbpixels * 2, channel[0].data, ndata);
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data + nbpixels;
		image->pdata[BLAYER] = image->data + nbpixels* 2;
	} else {
		memcpy(image->data, channel[0].data, ndata);
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

static int cvLucyRichardson_float(fits *image, double sigma, int iterations) {
	assert(image->fdata);
	assert(image->rx);
	assert(image->ry);

	int ksize = (int)((KERNEL_SIZE_FACTOR * sigma) + 0.5);
	ksize = ksize % 2 != 0 ? ksize : ksize + 1;
	Mat in, out;

	float *bgrbgr = NULL;

	if (image->naxes[2] == 1) {
		in = Mat(image->ry, image->rx, CV_32FC1, image->fdata);
	}
	else if (image->naxes[2] == 3) {
		bgrbgr = fits_to_bgrbgr_float(image);
		in = Mat(image->ry, image->rx, CV_32FC3, bgrbgr);
	}
	else {
		siril_log_message(_("Deconvolution is not supported for images with %d channels\n"), image->naxes[2]);
		return -1;
	}

	set_progress_bar_data(_("Deconvolution..."), PROGRESS_NONE);


	in.convertTo(in, CV_32FC3);

	Mat psf = getGaussianKernel(ksize, sigma);

	double mu = 0.01;

	out = RLTikh_deconvolution(in, psf, mu, iterations);

	std::vector<Mat> channel(3);
	split(out, channel);

	size_t nbpixels = image->naxes[0] * image->naxes[1];
	size_t ndata = nbpixels * sizeof(float);
	if (image->naxes[2] == 3) {
		memcpy(image->fdata, channel[2].data, ndata);
		memcpy(image->fdata + nbpixels, channel[1].data, ndata);
		memcpy(image->fdata + nbpixels * 2, channel[0].data, ndata);
		image->fpdata[RLAYER] = image->fdata;
		image->fpdata[GLAYER] = image->fdata + nbpixels;
		image->fpdata[BLAYER] = image->fdata + nbpixels* 2;
	} else {
		memcpy(image->fdata, channel[0].data, ndata);
		image->fpdata[RLAYER] = image->fdata;
		image->fpdata[GLAYER] = image->fdata;
		image->fpdata[BLAYER] = image->fdata;
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

int cvLucyRichardson(fits *image, double sigma, int iterations) {
	assert(image->rx);
	assert(image->ry);
	if (image->type == DATA_USHORT)
		return cvLucyRichardson_ushort(image, sigma, iterations);
	if (image->type == DATA_FLOAT)
		return cvLucyRichardson_float(image, sigma, iterations);
	return -1;
}


/* Work on grey images. If image is in RGB it must be first converted
 * in CieLAB. Then, only the first channel is applied
 */
static int cvClahe_ushort(fits *image, double clip_limit, int size) {
	assert(image->data);
	assert(image->rx);
	assert(image->ry);

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
				bgrbgr = fits_to_bgrbgr_ushort(image);
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

	size_t nbpixels = image->naxes[0] * image->naxes[1];
	size_t ndata = nbpixels * sizeof(WORD);
	if (image->naxes[2] == 3) {
		memcpy(image->data, channel[2].data, ndata);
		memcpy(image->data + nbpixels, channel[1].data, ndata);
		memcpy(image->data + nbpixels * 2, channel[0].data, ndata);
		image->pdata[RLAYER] = image->data;
		image->pdata[GLAYER] = image->data + nbpixels;
		image->pdata[BLAYER] = image->data + nbpixels* 2;
	} else {
		memcpy(image->data, channel[0].data, ndata);
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

static int cvClahe_float(fits *image, double clip_limit, int size) {
	assert(image->fdata);
	assert(image->rx);
	assert(image->ry);

	// preparing data
	Mat in, out;

	Ptr<CLAHE> clahe = createCLAHE();
	clahe->setClipLimit(clip_limit);
	clahe->setTilesGridSize(Size(size, size));

	if (image->naxes[2] == 3) {
		Mat lab_image;
		std::vector<Mat> lab_planes(3);
		float *bgrbgr;

		bgrbgr = fits_to_bgrbgr_float(image);
		in = Mat(image->ry, image->rx, CV_32FC3, bgrbgr);
		out = Mat();

		// convert the RGB color image to Lab
		cvtColor(in, lab_image, COLOR_BGR2Lab);

		// Extract the L channel
		split(lab_image, lab_planes); // now we have the L image in lab_planes[0]

		// apply the CLAHE algorithm to the L channel (does not work with 32F images)
		lab_planes[0].convertTo(lab_planes[0], CV_16U, USHRT_MAX_DOUBLE / 100.0);
		clahe->apply(lab_planes[0], lab_planes[0]);
		lab_planes[0].convertTo(lab_planes[0], CV_32F, 100.0 / USHRT_MAX_DOUBLE);

		// Merge the color planes back into an Lab image
		merge(lab_planes, lab_image);

		// convert back to RGB
		cvtColor(lab_image, out, COLOR_Lab2BGR);
		out.convertTo(out, CV_32FC3, 1.0);

		delete[] bgrbgr;

	} else {
		in = Mat(image->ry, image->rx, CV_32FC1, image->fdata);
		out = Mat();

		in.convertTo(in, CV_16U, USHRT_MAX_DOUBLE);
		clahe->apply(in, out);
		out.convertTo(out, CV_32F, 1.0 / USHRT_MAX_DOUBLE);
	}

	std::vector<Mat> channel(3);
	split(out, channel);

	size_t nbpixels = image->naxes[0] * image->naxes[1];
	size_t ndata = nbpixels * sizeof(float);
	if (image->naxes[2] == 3) {
		memcpy(image->fdata, channel[2].data, ndata);
		memcpy(image->fdata + nbpixels, channel[1].data, ndata);
		memcpy(image->fdata + nbpixels * 2, channel[0].data, ndata);
		image->fpdata[RLAYER] = image->fdata;
		image->fpdata[GLAYER] = image->fdata + nbpixels;
		image->fpdata[BLAYER] = image->fdata + nbpixels* 2;
	} else {
		memcpy(image->fdata, channel[0].data, ndata);
		image->fpdata[RLAYER] = image->fdata;
		image->fpdata[GLAYER] = image->fdata;
		image->fpdata[BLAYER] = image->fdata;
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

int cvClahe(fits *image, double clip_limit, int size) {
	assert(image->rx);
	assert(image->ry);
	if (image->type == DATA_USHORT)
		return cvClahe_ushort(image, clip_limit, size);
	if (image->type == DATA_FLOAT)
		return cvClahe_float(image, clip_limit, size);
	return -1;
}


