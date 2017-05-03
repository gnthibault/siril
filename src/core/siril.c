/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2017 team free-astro (see more in AUTHORS file)
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
 */

#include <gtk/gtk.h>
#ifdef MAC_INTEGRATION
#include "gtkmacintegration/gtkosxapplication.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_statistics.h>
#include <fitsio.h>
#include <complex.h>
#include <float.h>
#include <assert.h>
#include <libgen.h>
#include <fcntl.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "gui/callbacks.h"
#include "algos/colors.h"
#include "gui/histogram.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "algos/gradient.h"
#include "gui/PSF_list.h"
#include "opencv/opencv.h"
#include "algos/Def_Math.h"
#include "algos/Def_Wavelet.h"
#include "algos/cosmetic_correction.h"
#include "io/ser.h"

#define MAX_ITER 15
#define EPSILON 1E-4

/* this file contains all functions for image processing */

int threshlo(fits *fit, int level) {
	int i, layer;

	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		WORD *buf = fit->pdata[layer];
		for (i = 0; i < fit->rx * fit->ry; ++i) {
			*buf = max(level, *buf);
			buf++;
		}
	}
	return 0;
}

int threshhi(fits *fit, int level) {
	int i, layer;

	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		WORD *buf = fit->pdata[layer];
		for (i = 0; i < fit->rx * fit->ry; ++i) {
			*buf = min(level, *buf);
			buf++;
		}
	}
	return 0;
}

int nozero(fits *fit, int level) {
	int i, layer;

	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		WORD *buf = fit->pdata[layer];
		for (i = 0; i < fit->rx * fit->ry; ++i) {
			if (*buf == 0)
				*buf = level;
			buf++;
		}
	}
	return 0;
}

/*****************************************************************************
 *       S I R I L      A R I T H M E T I C      O P E R A T I O N S         *
 ****************************************************************************/

/* equivalent to (map simple_operation a), with simple_operation being
 * (lambda (pixel) (oper pixel scalar))
 * oper is a for addition, s for substraction (i for difference) and so on. */
int soper(fits *a, double scalar, char oper) {
	WORD *gbuf;
	int i, layer;
	int n = a->rx * a->ry;

	assert(n > 0);

	for (layer = 0; layer < a->naxes[2]; ++layer) {
		gbuf = a->pdata[layer];
		switch (oper) {
		case OPER_ADD:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD((double) gbuf[i] + scalar);
			}
			break;
		case OPER_SUB:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD((double) gbuf[i] - scalar);
			}
			break;
		case OPER_MUL:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD((double) gbuf[i] * scalar);
			}
			break;
		case OPER_DIV:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD((double) gbuf[i] / scalar);
			}
			break;
		}
	}
	return 0;
}

/* applies operation of image a with image b, for all their layers:
 * a = a oper b
 * returns 0 on success */
int imoper(fits *a, fits *b, char oper) {
	int i, layer;

	if (a->rx != b->rx || a->ry != b->ry) {
		siril_log_message(
				_("imoper: images don't have the same size (w = %u|%u, h = %u|%u)\n"),
				a->rx, b->rx, a->ry, b->ry);
		return 1;
	}
	for (layer = 0; layer < a->naxes[2]; ++layer) {
		WORD *buf = b->pdata[layer];
		WORD *gbuf = a->pdata[layer];
		int n = a->rx * a->ry;
		switch (oper) {
		case OPER_ADD:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD(gbuf[i] + buf[i]);
			}
			break;
		case OPER_SUB:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD(gbuf[i] - buf[i]);
			}
			break;
		case OPER_MUL:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD(gbuf[i] * buf[i]);
			}
			break;
		case OPER_DIV:
			for (i = 0; i < n; ++i) {
				gbuf[i] = round_to_WORD(gbuf[i] / buf[i]);
			}
			break;
		}
	}
	return 0;
}

/* This function applies a subtraction but contrary to Sub in imoper
 * it will use double type format.
 */
int sub_background(fits* image, fits* background, int layer) {
	double *pxl_image, *pxl_bkg, min = DBL_MAX;
	WORD *image_buf = image->pdata[layer];
	WORD *bkg_buf = background->pdata[layer];
	size_t i, ndata;

	if ((image->rx) != (background->rx) || ((image->ry) != (background->ry))) {
		char *msg = siril_log_message(
				_("Images don't have the same size (w = %d|%d, h = %d|%d)\n"),
				image->rx, background->rx, image->ry, background->ry);
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		return 1;
	}
	ndata = image->rx * image->ry;

	// First step we convert data, apply the subtraction and search for minimum
	pxl_image = malloc(sizeof(double) * ndata);
	pxl_bkg = malloc(sizeof(double) * ndata);
	for (i = 0; i < ndata; i++) {
		pxl_image[i] = (double) image_buf[i] / USHRT_MAX_DOUBLE;
		pxl_bkg[i] = (double) bkg_buf[i] / USHRT_MAX_DOUBLE;
		pxl_image[i] -= pxl_bkg[i];
		min = min(pxl_image[i], min);
	}
	image_buf = image->pdata[layer];
	// Second we apply an offset to the result and re-convert the data
	for (i = 0; i < ndata; i++) {
		pxl_image[i] += fabs(min);
		image_buf[i] = round_to_WORD(pxl_image[i] * USHRT_MAX_DOUBLE);
	}

	// We free memory
	free(pxl_image);
	free(pxl_bkg);
	return 0;
}

int addmax(fits *a, fits *b) {
	WORD *gbuf[3] = { a->pdata[RLAYER], a->pdata[GLAYER], a->pdata[BLAYER] };
	WORD *buf[3] = { b->pdata[RLAYER], b->pdata[GLAYER], b->pdata[BLAYER] };
	gint i, layer;

	if (a->rx != b->rx || a->ry != b->ry || a->naxes[2] != b->naxes[2]) {
		siril_log_message(
				_("addmax: images don't have the same size (w = %d|%d, h = %d|%d, layers = %d|%d)\n"),
				a->rx, b->rx, a->ry, b->ry, a->naxes[2], b->naxes[2]);
		return 1;
	}
	assert(a->naxes[2] == 1 || a->naxes[2] == 3);

	for (layer = 0; layer < a->naxes[2]; ++layer) {
		for (i = 0; i < a->ry * a->rx; ++i) {
			if (buf[layer][i] > gbuf[layer][i])
				gbuf[layer][i] = buf[layer][i];
		}
	}
	return 0;
}

/* If fdiv is ok, function returns 0. If overflow, fdiv returns 1*/
int fdiv(fits *a, fits *b, float coef) {
	int i, layer;
	int retvalue = 0;
	double temp;

	if (a->rx != b->rx || a->ry != b->ry || a->naxes[2] != b->naxes[2]) {
		fprintf(stderr, "Wrong size or channel count: %u=%u? / %u=%u?\n", a->rx,
				b->rx, a->ry, b->ry);
		return -1;
	}
	for (layer = 0; layer < a->naxes[2]; ++layer) {
		WORD *buf = b->pdata[layer];
		WORD *gbuf = a->pdata[layer];
		for (i = 0; i < b->rx * b->ry; ++i) {
			if (buf[i] == 0)
				buf[i] = 1;		// avoid division by 0
			temp = ((double) coef * ((double) gbuf[i] / (double) buf[i]));
			if (temp > USHRT_MAX_DOUBLE)
				retvalue = 1;
			gbuf[i] = round_to_WORD(temp);
		}
	}
	return retvalue;
}

/* normalized division a/b, stored in a, with max value equal to the original
 * max value of a, for each layer. */
int ndiv(fits *a, fits *b) {
	double *div;
	int layer, i, nb_pixels;
	if (a->rx != b->rx || a->ry != b->ry || a->naxes[2] != b->naxes[2]) {
		fprintf(stderr,
				"Wrong size or channel count: %u=%u? / %u=%u?, %ld=%ld?\n",
				a->rx, b->rx, a->ry, b->ry, a->naxes[2], b->naxes[2]);
		return 1;
	}
	nb_pixels = a->rx * a->ry;
	div = malloc(nb_pixels * sizeof(double));

	for (layer = 0; layer < a->naxes[2]; ++layer) {
		double max = 0, norm;
		for (i = 0; i < nb_pixels; ++i) {
			if (!b->pdata[layer][i])
				div[i] = (double) a->pdata[layer][i];
			else
				div[i] = (double) a->pdata[layer][i]
						/ (double) b->pdata[layer][i];
			max = max(div[i], max);
		}
		norm = max / (double) a->max[layer];
		for (i = 0; i < nb_pixels; ++i) {
			a->pdata[layer][i] = round_to_WORD(div[i] / norm);
		}
	}

	free(div);
	return 0;
}

/**********************************************************
 *
 */

#ifdef HAVE_OPENCV
int unsharp(fits *fit, double sigma, double amount, gboolean verbose) {
	struct timeval t_start, t_end;

	if (sigma <= 0.0)
		return 1;
	if (verbose) {
		siril_log_color_message(_("Unsharp: processing...\n"), "red");
		gettimeofday(&t_start, NULL);
	}
	cvUnsharpFilter(fit, sigma, amount);

	if (verbose) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}
	return 0;
}

#else
double gaussienne(double sigma, int size, double *gauss) {
	double s2, s;
	int i,j,n;

	n=size/2;
	s2=sigma*sigma;
	s=(double)0;
	for (i=0;i<size;++i) {
		for (j=0;j<size;++j) {
			s+=gauss[i*size+j]=exp(-((double)(i-n)*(i-n) + (double)(j-n)*(j-n))/2/s2);
			//~ fprintf(stderr,"%d:%d %f %f \n", i, j, s, gauss[i*size+j]);
		}
	}
	return s;
}

int unsharp(fits *fit, double sigma, double mult, gboolean verbose) {
	// if fabs(mult) > 0.01, unsharp computes the unsharp mask
	// else unsharp computes a gaussian filter

	double normalize,g, *gauss/*, *gaussbuf*/;
	WORD *buf, *gbuf;
	int size, ss2, stride, i, j, k, l, layer;
	struct timeval t_start, t_end;

	if (verbose) {
		siril_log_color_message("Unsharp: processing...\n", "red");
		gettimeofday (&t_start, NULL);
	}

	//	fprintf(stderr,"gfitrx: %d mult :%f\n", gfit.rx, mult);
	//~ size=(int)(4*sigma);
	size = (int)(2*(((sigma - 0.8) / 0.3) + 1));// heuristic. Just to be homogeneous with opencv
	if (!(size%2))
	++size;

	ss2=size/2;
	gauss=malloc(size*size*sizeof(double));
	if (!gauss) {
		perror("unsharp.c: Alloc error for gauss");
		return -1;
	}
	//	fprintf(stderr,"size: %d sigma:%f\n", size, sigma);

	normalize=gaussienne(sigma, size, gauss);
	//	fprintf(stderr,"gfitrx: %d ss2:%d\n", gfit.rx, ss2);
	wfit[4].data=(WORD *)calloc(1, fit->rx*fit->ry*sizeof(WORD));

	wfit[4].rx=fit->rx;
	wfit[4].ry=fit->ry;
	wfit[4].lo=fit->lo;
	wfit[4].hi=fit->hi;
	stride=wfit[4].rx-size;
	for (layer=0; layer<fit->naxes[2]; ++layer) {
		memcpy(wfit[4].data, fit->pdata[layer], fit->rx*fit->ry*sizeof(WORD));

		buf=fit->pdata[layer]+ss2+ss2*fit->rx;
		//	fprintf(stderr,"gfitrx: %d ss2:%d\n", gfit.rx, ss2);
		for (i=ss2;i<fit->ry-ss2;++i) {
			for (j=ss2;j<fit->rx-ss2;++j) {
				g=(double)0;
				gbuf=wfit[4].data+(i-ss2)*fit->rx+j-ss2;
				//~ gaussbuf=gauss;
				for (k=0;k<size;++k) {
					for(l=0;l<size;++l) {
						g+=(*gbuf++)*(gauss[k*size+l]);
					}
					gbuf+=stride;
				}
				*(buf++)=g/normalize;
			}
			buf+=ss2+ss2;
		}

		buf=fit->pdata[layer];
		gbuf=wfit[4].data;
		if (fabs(mult)>0.0) {
			for (i=0; i<fit->rx * fit->ry; i++) {
				//~ double tmp = gbuf[i] + mult * (gbuf[i] - buf[i]);
				double tmp = gbuf[i] * (1.0 + mult) + buf[i] * (- mult);
				if (tmp < 0.0) buf[i] = 0;
				else if (tmp > USHRT_MAX_DOUBLE) buf[i] = USHRT_MAX;
				else buf[i] = (WORD) tmp;
			}
		}
	}
	free(gauss);
	clearfits(&wfit[4]);
	if (verbose) {
		gettimeofday (&t_end, NULL);
		show_time(t_start, t_end);
	}
	return 0;
}
#endif

/* inplace cropping of the image in fit
 * fit->data is not realloc, only fit->pdata points to a different area and
 * data is correctly written to this new area, which makes this function
 * quite dangerous to use when fit is used for something else afterwards.
 */
int crop(fits *fit, rectangle *bounds) {
	int i, j, layer;
	int newnbdata;
	struct timeval t_start, t_end;

	memset(&t_start, 0, sizeof(struct timeval));
	memset(&t_end, 0, sizeof(struct timeval));

	if (fit == &gfit) {
		siril_log_color_message(_("Crop: processing...\n"), "red");
		gettimeofday(&t_start, NULL);
	}

	newnbdata = bounds->w * bounds->h;
	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		WORD *from = fit->pdata[layer]
				+ (fit->ry - bounds->y - bounds->h) * fit->rx + bounds->x;
		fit->pdata[layer] = fit->data + layer * newnbdata;
		WORD *to = fit->pdata[layer];
		int stridefrom = fit->rx - bounds->w;

		for (i = 0; i < bounds->h; ++i) {
			for (j = 0; j < bounds->w; ++j) {
				*to++ = *from++;
			}
			from += stridefrom;
		}
	}
	fit->rx = fit->naxes[0] = bounds->w;
	fit->ry = fit->naxes[1] = bounds->h;

	if (fit == &gfit) {
		clear_stars_list();
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}

	return 0;
}

/* takes the image in gfit, copies it in a temporary fit to shift it, and copy it back into gfit */
/* TODO: it can be done in the same, thus avoiding to allocate, it just needs to care
 * about the sign of sx and sy to avoid data overwriting in the same allocated space. */
int shift(int sx, int sy) {
	int x, y, nx, ny, i, ii, layer;
	fits tmpfit;
	copyfits(&(gfit), &tmpfit, CP_ALLOC | CP_FORMAT, 0);
	i = 0;
	/* the loop is the same than in composit() */
	for (y = 0; y < gfit.ry; ++y) {
		for (x = 0; x < gfit.rx; ++x) {
			nx = (x - sx);
			ny = (y - sy);
			//printf("x=%d y=%d sx=%d sy=%d i=%d ii=%d\n",x,y,shiftx,shifty,i,ii);
			if (nx >= 0 && nx < gfit.rx && ny >= 0 && ny < gfit.ry) {
				ii = ny * gfit.rx + nx;
				//printf("shiftx=%d shifty=%d i=%d ii=%d\n",shiftx,shifty,i,ii);
				if (ii > 0 && ii < gfit.rx * gfit.ry) {
					for (layer = 0; layer < gfit.naxes[2]; ++layer) {
						tmpfit.pdata[layer][i] = gfit.pdata[layer][ii];
					}
				}
			}
			++i;
		}
	}

	for (layer = 0; layer < gfit.naxes[2]; ++layer) {
		memcpy(gfit.pdata[layer], tmpfit.pdata[layer],
				gfit.rx * gfit.ry * sizeof(WORD));
	}
	free(tmpfit.data);

	return 0;
}

#if 0
int rshift2(char *genname, char *outname, int number, char *shiftfile) {
	int j,shiftx,shifty,count,n;
	char line[256];
	FILE *sf=NULL;

	if (shiftfile!=NULL) {
		sf=fopen(shiftfile,"r");
		if (sf==NULL) {
			siril_log_message("rshift2: could not open shift file %s\n", shiftfile);
			return 0;
		}
		fgets(line,255,sf);
		while(line[0]=='#') {
			fgets(line,255,sf);
		}
	}
	shiftx=shifty=0;

	for (j=1;j<=number;++j) {
		if(sf!=NULL) {
			count=sscanf(line,"%d %d %d",&n,&shiftx,&shifty);
			if(count!=3) {
				siril_log_message("rshift2: format error in shift file %s\n", shiftfile);
				return 0;
			}
			fgets(line,255,sf);
		}
		else {
			shiftx=0;
			shifty=0;
			n=j;
		}
		fprintf(stderr,"rshift2 %d %d %d\n", n, shiftx,shifty);
		siril_log_message("Processing image %s%d%s", genname, n, com.ext);

		buildfilename(genname,n);
		if(readfits(com.formname, &(gfit), com.formname)) {
			fprintf(stderr,"missed %s com.formanme",com.formname);
			return 1;
		}
		shift(shiftx, shifty);
		buildfilename(outname,n);
		savefits(com.formname,&gfit);
	}
	if(sf!=NULL) {
		fclose(sf);
	}
	return 0;
}
#endif

/* This entropy function computes the entropy for the image in gfit for its
 * layer 'layer', in the area designated by area which can be NULL.
 * An optional imstats parameter can be used to provide the background and
 * sigma value, and when it is given, the entropy will only be computed for
 * pixels with values above background + 1 * sigma. It must be NULL otherwise.
 */
double entropy(fits *fit, int layer, rectangle *area, imstats *opt_stats) {
	double e = 0.0, threshold = 0.0;
	gsl_histogram *histo;
	size_t i, size, n;

	if (opt_stats && opt_stats->median >= 0.0 && opt_stats->sigma >= 0.0)
		threshold = opt_stats->median + 1 * opt_stats->sigma;

	if (area == NULL)
		histo = computeHisto(fit, layer);
	else
		histo = computeHisto_Selection(fit, layer, area);

	n = fit->rx * fit->ry;
	assert (n > 0);
	size = gsl_histogram_bins(histo);
	for (i = 0; i < size; i++) {
		double p = gsl_histogram_get(histo, i);
		if (p > threshold && p < size)
			e += (p / n) * log(n / p);
	}
	gsl_histogram_free(histo);

	return e;
}

int loglut(fits *fit, int dir) {
	// This function maps fit with a log LUT
	int i, layer;
	gdouble normalisation, temp;
	WORD *buf[3] =
			{ fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	assert(fit->naxes[2] <= 3);

	normalisation = USHRT_MAX_DOUBLE / log(USHRT_MAX_DOUBLE);

	for (i = 0; i < fit->ry * fit->rx; i++) {
		for (layer = 0; layer < fit->naxes[2]; ++layer) {
			temp = buf[layer][i] + 1;
			if (dir == LOG)
				buf[layer][i] = normalisation * log(temp);
			else
				buf[layer][i] = exp(temp / normalisation);
		}
	}
	return 0;
}

double contrast(fits* fit, int layer) {
	int i;
	WORD *buf = fit->pdata[layer];
	double contrast = 0.0;
	imstats *stat = statistics(fit, layer, &com.selection, STATS_BASIC, STATS_ZERO_NULLCHECK);
	if (!stat) {
		siril_log_message(_("Error: no data computed.\n"));
		return -1.0;
	}
	double mean = stat->mean;
	free(stat);

	for (i = 0; i < fit->rx * fit->ry; i++)
		contrast += SQR((double )buf[i] - mean);
	contrast /= (fit->rx * fit->ry);
	return contrast;
}

int ddp(fits *a, int level, float coeff, float sigma) {
	copyfits(a, &wfit[0], CP_ALLOC | CP_COPYA | CP_FORMAT, 0);
	unsharp(&wfit[0], sigma, 0, FALSE);
	soper(&wfit[0], (double) level, OPER_ADD);
	nozero(&wfit[0], 1);
	fdiv(a, &wfit[0], level);
	soper(a, (double) coeff, OPER_MUL);
	clearfits(&wfit[0]);
	return 0;
}

int visu(fits *fit, int low, int high) {
	if (low < 0 || low > USHRT_MAX || high < 1 || high > USHRT_MAX)
		return 1;
	if (single_image_is_loaded() && com.cvport < com.uniq->nb_layers) {
		com.uniq->layers[com.cvport].hi = high;
		com.uniq->layers[com.cvport].lo = low;
	} else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers) {
		com.seq.layers[com.cvport].hi = high;
		com.seq.layers[com.cvport].lo = low;
	} else
		return 1;
	set_cutoff_sliders_values();
	redraw(com.cvport, REMAP_ONLY);
	redraw_previews();
	return 0;
}

/* fill an image or selection with the value 'level' */
int fill(fits *fit, int level, rectangle *arearg) {
	WORD *buf;
	int i, j, layer;
	rectangle area;

	if (arearg) {
		memcpy(&area, arearg, sizeof(rectangle));
	} else {
		if (com.selection.h && com.selection.w) {
			memcpy(&area, &com.selection, sizeof(rectangle));
		} else {
			area.w = fit->rx;
			area.h = fit->ry;
			area.x = 0;
			area.y = 0;
		}
	}
	for (layer = 0; layer < fit->naxes[2]; ++layer) {
		buf = fit->pdata[layer] + (fit->ry - area.y - area.h) * fit->rx
				+ area.x;
		int stridebuf = fit->rx - area.w;
		for (i = 0; i < area.h; ++i) {
			for (j = 0; j < area.w; ++j) {
				*buf++ = level;
			}
			buf += stridebuf;
		}
	}
	return 0;
}

int off(fits *fit, int level) {
	WORD *buf[3] =
			{ fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	int i, layer;
	assert(fit->naxes[2] <= 3);
	if (level == 0)
		return 0;
	if (level < -USHRT_MAX)
		level = -USHRT_MAX;
	else if (level > USHRT_MAX)
		level = USHRT_MAX;
	for (i = 0; i < fit->rx * fit->ry; ++i) {
		for (layer = 0; layer < fit->naxes[2]; ++layer) {
			WORD val = buf[layer][i];
			if ((level < 0 && val < -level))
				buf[layer][i] = 0;
			else if (level > 0 && val > USHRT_MAX - level)
				buf[layer][i] = USHRT_MAX;
			else
				buf[layer][i] = val + level;
		}
	}
	return 0;
}

void mirrorx(fits *fit, gboolean verbose) {
	int line, axis, line_size;
	WORD *swapline, *src, *dst;
	struct timeval t_start, t_end;

	if (verbose) {
		siril_log_color_message(_("Horizontal mirror: processing...\n"), "red");
		gettimeofday(&t_start, NULL);
	}

	line_size = fit->rx * sizeof(WORD);
	swapline = malloc(line_size);

	for (axis = 0; axis < fit->naxes[2]; axis++) {
		for (line = 0; line < fit->ry / 2; line++) {
			src = fit->pdata[axis] + line * fit->rx;
			dst = fit->pdata[axis] + (fit->ry - line - 1) * fit->rx;

			memcpy(swapline, src, line_size);
			memcpy(src, dst, line_size);
			memcpy(dst, swapline, line_size);
		}
	}
	free(swapline);
	if (verbose) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}
}

void mirrory(fits *fit, gboolean verbose) {
	struct timeval t_start, t_end;

	if (verbose) {
		siril_log_color_message(_("Vertical mirror: processing...\n"), "red");
		gettimeofday(&t_start, NULL);
	}

	fits_flip_top_to_bottom(fit);
	fits_rotate_pi(fit);

	if (verbose) {
		gettimeofday(&t_end, NULL);
		show_time(t_start, t_end);
	}
}

/* this method rotates the image 180 degrees, useful after german mount flip.
 * fit->rx, fit->ry, fit->naxes[2] and fit->pdata[*] are required to be assigned correctly */
void fits_rotate_pi(fits *fit) {
	int i, line, axis, line_size;
	WORD *line1, *line2, *src, *dst, swap;

	line_size = fit->rx * sizeof(WORD);
	line1 = malloc(line_size);
	line2 = malloc(line_size);

	for (axis = 0; axis < fit->naxes[2]; axis++) {
		for (line = 0; line < fit->ry / 2; line++) {
			src = fit->pdata[axis] + line * fit->rx;
			dst = fit->pdata[axis] + (fit->ry - line - 1) * fit->rx;

			memcpy(line1, src, line_size);
			for (i = 0; i < fit->rx / 2; i++) {
				swap = line1[i];
				line1[i] = line1[fit->rx - i - 1];
				line1[fit->rx - i - 1] = swap;
			}
			memcpy(line2, dst, line_size);
			for (i = 0; i < fit->rx / 2; i++) {
				swap = line2[i];
				line2[i] = line2[fit->rx - i - 1];
				line2[fit->rx - i - 1] = swap;
			}
			memcpy(src, line2, line_size);
			memcpy(dst, line1, line_size);
		}
		if (fit->ry & 1) {
			/* swap the middle line */
			src = fit->pdata[axis] + line * fit->rx;
			for (i = 0; i < fit->rx / 2; i++) {
				swap = src[i];
				src[i] = src[fit->rx - i - 1];
				src[fit->rx - i - 1] = swap;
			}
		}
	}
	free(line1);
	free(line2);
}

/* This function fills the data in the lrgb image with LRGB information from l, r, g and b
 * images. Layers are not aligned, images need to be all of the same size.
 * It may be used in the command line, currently unused. */
int lrgb(fits *l, fits *r, fits *g, fits *b, fits *lrgb) {
	//
	// Combines l r g and b components into resulting lrgb
	// We transform each pixel from RGB to HSI,
	// then take I from the luminance l fits and
	// immediately step back to RGB to the working copy
	//
	guint x, y;
	gdouble rr, gg, bb, h, s, i/*, ps3, dps3, qps3, dpi*/;
	gint maxi;
	WORD *pr, *pg, *pb, *dr, *dg, *db, *pl;

	//
	// some stats used to normalize
	//
	image_find_minmax(r, 0);
	image_find_minmax(g, 0);
	image_find_minmax(b, 0);
	maxi = max(r->maxi, max(g->maxi, b->maxi));
	image_find_minmax(l, 0);
	//
	// initialize pointers
	//
	pr = r->data;
	pg = g->data;
	pb = b->data;
	pl = l->data;
	dr = lrgb->pdata[RLAYER];
	dg = lrgb->pdata[GLAYER];
	db = lrgb->pdata[BLAYER];
	//
	// some trigo constants
	// we stick to h in radians, not in degrees
	//
	//dpi=2*M_PI;
	//ps3=M_PI/3;
	//dps3=2*M_PI;
	//dps3=2*M_PI/3;
	//qps3=4*M_PI/3;
	//
	// Main loop
	//
	fprintf(stderr, "HSI->RGB %u %u\n", r->ry, r->rx);
	for (y = 0; y < r->ry; y++) {
		for (x = 0; x < r->rx; x++) {
			//
			// First normalize rgb to [0 1]
			//
			rr = (double) (*pr++) / maxi;
			gg = (double) (*pg++) / maxi;
			bb = (double) (*pb++) / maxi;

			rgb_to_hsl(rr, gg, bb, &h, &s, &i);
			//
			// replace luminance
			//
			i = *pl++ / (double) l->maxi;
			//
			// and back to RGB
			hsl_to_rgb(h, s, i, &rr, &gg, &bb);
			//
			// now denormalize and store
			//
			*dr++ = (WORD) (rr * maxi);
			*dg++ = (WORD) (gg * maxi);
			*db++ = (WORD) (bb * maxi);
		}
	}
	return 0;
}

static double evaluateNoiseOfCalibratedImage(fits *fit, fits *dark, double k) {
	double noise = 0;
	fits *dark_tmp;
	fits *fit_tmp;
	int chan;

	dark_tmp = calloc(1, sizeof(fits));
	fit_tmp = calloc(1, sizeof(fits));

	new_fit_image(dark_tmp, dark->rx, dark->ry, 1);
	new_fit_image(fit_tmp, fit->rx, fit->ry, 1);

	copyfits(dark, dark_tmp, CP_ALLOC | CP_EXTRACT, 0);
	copyfits(fit, fit_tmp, CP_ALLOC | CP_EXTRACT, 0);

	soper(dark_tmp, k, OPER_MUL);
	imoper(fit_tmp, dark_tmp, OPER_SUB);

	for (chan = 0; chan < fit->naxes[2]; chan++) {
		imstats *stat = NULL;
		stat = statistics(fit_tmp, chan, NULL, STATS_BASIC, STATS_ZERO_NULLCHECK);
		if (!stat) {
			siril_log_message(_("Error: no data computed.\n"));
			return 0.0;
		}
		noise += stat->bgnoise;
		//printf("noise=%lf, k=%lf\n", noise, k);
		free(stat);
	}
	clearfits(dark_tmp);
	clearfits(fit_tmp);

	return noise;
}

#define GR ((sqrt(5) - 1) / 2)

static double goldenSectionSearch(fits *brut, fits *dark, double a, double b,
		double tol) {
	double c, d;
	double fc, fd;

	c = b - GR * (b - a);
	d = a + GR * (b - a);
	do {
		fc = evaluateNoiseOfCalibratedImage(brut, dark, c);
		fd = evaluateNoiseOfCalibratedImage(brut, dark, d);
		if (fc < fd) {
			b = d;
			d = c;
			c = b - GR * (b - a);
		} else {
			a = c;
			c = d;
			d = a + GR * (b - a);
		}
	} while (fabs(c - d) > tol);
	return ((b + a) / 2);
}

static int preprocess(fits *brut, fits *offset, fits *dark, fits *flat, float level) {

	if (com.preprostatus & USE_OFFSET) {
		imoper(brut, offset, OPER_SUB);
	}

	/* if dark optimization, the master-dark has already been subtracted */
	if ((com.preprostatus & USE_DARK) && !(com.preprostatus & USE_OPTD)) {
		imoper(brut, dark, OPER_SUB);
	}

	if (com.preprostatus & USE_FLAT) {
		fdiv(brut, flat, level);
	}

	return 0;
}

static int darkOptimization(fits *brut, fits *dark, fits *offset) {
	double k;
	double lo = 0.0;
	double up = 2.0;

	fits *dark_tmp = calloc(1, sizeof(fits));
	new_fit_image(dark_tmp, dark->rx, dark->ry, 1);
	copyfits(dark, dark_tmp, CP_ALLOC | CP_EXTRACT, 0);

	/* Minimization of background noise to find better k */
	k = goldenSectionSearch(brut, dark_tmp, lo, up, 1E-3);

	siril_log_message(_("Dark optimization: %.3lf\n"), k);
	/* Multiply coefficient to master-dark */
	if (com.preprostatus & USE_OFFSET)
		imoper(dark_tmp, offset, OPER_SUB);
	soper(dark_tmp, k, OPER_MUL);
	imoper(brut, dark_tmp, OPER_SUB);

	clearfits(dark_tmp);

	return 0;
}

// idle function executed at the end of the sequence preprocessing
static gboolean end_sequence_prepro(gpointer p) {
	struct preprocessing_data *args = (struct preprocessing_data *) p;
	struct timeval t_end;
	fprintf(stdout, "Ending sequence prepro idle function, retval=%d\n",
			args->retval);
	stop_processing_thread();// can it be done here in case there is no thread?
	set_cursor_waiting(FALSE);
	gettimeofday(&t_end, NULL);
	show_time(args->t_start, t_end);
	update_used_memory();
	if (!args->retval && !single_image_is_loaded()) {
		// load the new sequence
		char *ppseqname = malloc(
				strlen(com.seq.ppprefix) + strlen(com.seq.seqname) + 5);
		sprintf(ppseqname, "%s%s.seq", com.seq.ppprefix, com.seq.seqname);
		check_seq(0);
		update_sequences_list(ppseqname);
		free(ppseqname);
	}
	sequence_free_preprocessing_data(&com.seq);
#ifdef MAC_INTEGRATION
	GtkosxApplication *osx_app = gtkosx_application_get();
	gtkosx_application_attention_request(osx_app, INFO_REQUEST);
	g_object_unref (osx_app);
#endif
	free(args);
	return FALSE;
}

/* doing the preprocessing. No unprotected GTK+ calls can go there.
 * returns 1 on error */
gpointer seqpreprocess(gpointer p) {
	char dest_filename[256], msg[256];
	fits *dark, *offset, *flat;
	struct preprocessing_data *args = (struct preprocessing_data *) p;

	if (single_image_is_loaded()) {
		dark = com.uniq->dark;
		offset = com.uniq->offset;
		flat = com.uniq->flat;
	} else if (sequence_is_loaded()) {
		dark = com.seq.dark;
		offset = com.seq.offset;
		flat = com.seq.flat;
	} else
		return GINT_TO_POINTER(1);

	if (com.preprostatus & USE_FLAT) {
		if (args->autolevel) {
			/* TODO: evaluate the layer to apply but generally RLAYER is a good choice.
			 * Indeed, if it is image from APN, CFA picture are in black & white */
			imstats *stat = statistics(flat, RLAYER, NULL, STATS_BASIC, STATS_ZERO_NULLCHECK);
			if (!stat) {
				siril_log_message(_("Error: no data computed.\n"));
				return GINT_TO_POINTER(1);
			}
			args->normalisation = stat->mean;
			siril_log_message(_("Normalisation value auto evaluated: %.2lf\n"),
					args->normalisation);
			free(stat);
		}
	}

	if (single_image_is_loaded()) {
		snprintf(msg, 255, _("Pre-processing image %s"), com.uniq->filename);
		msg[255] = '\0';
		set_progress_bar_data(msg, 0.5);

		if ((com.preprostatus & USE_OPTD) && (com.preprostatus & USE_DARK))
			darkOptimization(com.uniq->fit, dark, offset);

		preprocess(com.uniq->fit, offset, dark, flat, args->normalisation);

		if ((com.preprostatus & USE_COSME) && (com.preprostatus & USE_DARK)) {
			if (dark->naxes[2] == 1) {
			/* Cosmetic correction */
			long icold, ihot;
			deviant_pixel *dev = find_deviant_pixels(dark, args->sigma, &icold, &ihot);
			siril_log_message(_("%ld pixels corrected (%ld + %ld)\n"),
					icold + ihot, icold, ihot);
			cosmeticCorrection(com.uniq->fit, dev, icold + ihot, args->is_cfa);
			if (dev)
				free(dev);
			}
			else
				siril_log_message(_("Darkmap cosmetic correction "
						"is only supported with single channel images\n"));
		}

		gchar *filename = g_path_get_basename(com.uniq->filename);
		char *filename_noext = remove_ext_from_filename(filename);
		snprintf(dest_filename, 255, "%s%s", com.uniq->ppprefix, filename_noext);
		dest_filename[255] = '\0';
		snprintf(msg, 255, _("Saving image %s"), filename_noext);
		msg[255] = '\0';
		set_progress_bar_data(msg, PROGRESS_NONE);
		savefits(dest_filename, com.uniq->fit);
		g_free(filename);
		free(filename_noext);
	} else {	// sequence
		struct ser_struct *new_ser_file = NULL;
		char source_filename[256];
		int i;
		long icold = 0L, ihot = 0L;
		deviant_pixel *dev = NULL;

		// creating a SER file if the input data is SER
		if (com.seq.type == SEQ_SER) {
			char new_ser_filename[256];
			new_ser_file = calloc(1, sizeof(struct ser_struct));
			snprintf(new_ser_filename, 255, "%s%s", com.seq.ppprefix, com.seq.ser_file->filename);
			ser_create_file(new_ser_filename, new_ser_file, TRUE, com.seq.ser_file);
		}

		if ((com.preprostatus & USE_COSME) && (com.preprostatus & USE_DARK)) {
			if (dark->naxes[2] == 1) {
				dev = find_deviant_pixels(dark, args->sigma, &icold, &ihot);
				siril_log_message(_("%ld pixels corrected (%ld + %ld)\n"),
						icold + ihot, icold, ihot);
			} else
				siril_log_message(_("Darkmap cosmetic correction "
						"is only supported with single channel images\n"));
		}

		fits *fit = calloc(1, sizeof(fits));

		for (i = 0; i < com.seq.number; i++) {
			if (!get_thread_run())
				break;
			seq_get_image_filename(&com.seq, i, source_filename);
			snprintf(msg, 255, _("Loading and pre-processing image %d/%d (%s)"),
					i + 1, com.seq.number, source_filename);
			msg[255] = '\0';
			set_progress_bar_data(msg,
					(double) (i + 1) / (double) com.seq.number);
			if (seq_read_frame(&com.seq, i, fit)) {
				snprintf(msg, 255, _("Could not read one of the raw files: %s."
						" Aborting preprocessing."), source_filename);
				msg[255] = '\0';
				set_progress_bar_data(msg, PROGRESS_RESET);
				args->retval = 1;
				if (com.seq.type == SEQ_SER) {
					ser_close_file(new_ser_file);
					free(new_ser_file);
				}
				gdk_threads_add_idle(end_sequence_prepro, args);
				return GINT_TO_POINTER(1);
			}
			if ((com.preprostatus & USE_OPTD) && (com.preprostatus & USE_DARK))
				darkOptimization(fit, dark, offset);

			preprocess(fit, offset, dark, flat, args->normalisation);

			if ((com.preprostatus & USE_COSME) && (com.preprostatus & USE_DARK) && (dark->naxes[2] == 1))
				cosmeticCorrection(fit, dev, icold + ihot, args->is_cfa);

			snprintf(dest_filename, 255, "%s%s", com.seq.ppprefix,
					source_filename);
			dest_filename[255] = '\0';
			snprintf(msg, 255, "Saving image %d/%d (%s)", i + 1, com.seq.number,
					dest_filename);
			if (com.seq.type == SEQ_SER) {
				ser_write_frame_from_fit(new_ser_file, fit, i);
			} else {
				savefits(dest_filename, fit);
			}
			clearfits(fit);
		}
		free(fit);
		// closing SER file if it applies
		if (com.seq.type == SEQ_SER && (new_ser_file != NULL)) {
			close(new_ser_file->fd);
			free(new_ser_file);
			new_ser_file = NULL;
		}
		set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
		if (dev) free(dev);
	}
	args->retval = 0;
	gdk_threads_add_idle(end_sequence_prepro, args);
	return GINT_TO_POINTER(0);
}

/* computes the background value using the histogram and/or median value.
 * The argument layer can be -1 for automatic setting (= green for RGB) */
double background(fits* fit, int reqlayer, rectangle *selection) {
	int layer = RLAYER;
	double bg;

	if (reqlayer >= 0)
		layer = reqlayer;
	else if (isrgb(&gfit))
		layer = GLAYER;		//GLAYER is better to evaluate background

	imstats* stat = statistics(fit, layer, selection, STATS_BASIC, STATS_ZERO_NULLCHECK);
	if (!stat) {
		siril_log_message(_("Error: no data computed.\n"));
		return 0.0;
	}
	bg = stat->median;

	free(stat);
	stat = NULL;
	return bg;
}

void show_FITS_header(fits *fit) {
	if (fit->header)
		show_data_dialog(fit->header, "FITS Header");
}

#ifdef HAVE_OPENCV
/* These functions do not more than resize_gaussian and rotate_image
 * except for console outputs. 
 * Indeed, siril_log_message seems not working in a cpp file */
int verbose_resize_gaussian(fits *image, int toX, int toY, int interpolation) {
	int retvalue;
	const char *str_inter;
	struct timeval t_start, t_end;

	switch (interpolation) {
	case OPENCV_NEAREST:
		str_inter = _("Nearest-Neighbor");
		break;
	default:
	case OPENCV_LINEAR:
		str_inter = _("Bilinear");
		break;
	case OPENCV_AREA:
		str_inter = _("Pixel Area Relation");
		break;
	case OPENCV_CUBIC:
		str_inter = _("Bicubic");
		break;
	case OPENCV_LANCZOS4:
		str_inter = _("Lanczos4");
		break;
	}

	siril_log_color_message(_("Resample (%s interpolation): processing...\n"),
			"red", str_inter);

	gettimeofday(&t_start, NULL);

	retvalue = cvResizeGaussian(&gfit, toX, toY, interpolation);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	return retvalue;
}

int verbose_rotate_image(fits *image, double angle, int interpolation,
		int cropped) {
	const char *str_inter;
	struct timeval t_start, t_end;

	switch (interpolation) {
	case -1:
		str_inter = _("No");
		break;
	case OPENCV_NEAREST:
		str_inter = _("Nearest-Neighbor");
		break;
	default:
	case OPENCV_LINEAR:
		str_inter = _("Bilinear");
		break;
	case OPENCV_AREA:
		str_inter = _("Pixel Area Relation");
		break;
	case OPENCV_CUBIC:
		str_inter = _("Bicubic");
		break;
	case OPENCV_LANCZOS4:
		str_inter = _("Lanczos4");
		break;
	}

	siril_log_color_message(
			_("Rotation (%s interpolation, angle=%g): processing...\n"), "red",
			str_inter, angle);
	gettimeofday(&t_start, NULL);

	cvRotateImage(&gfit, angle, interpolation, cropped);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);

	return 0;
}

#endif

/* This function computes wavelets with the number of Nbr_Plan and
 * extracts plan "Plan" in fit parameters */

int get_wavelet_layers(fits *fit, int Nbr_Plan, int Plan, int Type, int reqlayer) {
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];
	const char *tmpdir;
	int chan, start, end;
	wave_transf_des Wavelet[3];

	assert(fit->naxes[2] <= 3);
	tmpdir = g_get_tmp_dir();

	float *Imag = f_vector_alloc(fit->ry * fit->rx);
	if (Imag == NULL)
		return 1;

	if (reqlayer < 0 || reqlayer > 3) {
		start = 0;
		end = fit->naxes[2];
	}
	else {
		start = reqlayer;
		end = start + 1;
	}
	for (chan = start; chan < end; chan++) {
		int Nl, Nc;

		dir[chan] = malloc(
				strlen(tmpdir) + strlen(File_Name_Transform[chan]) + 2);
		strcpy(dir[chan], tmpdir);
		strcat(dir[chan], "/");
		strcat(dir[chan], File_Name_Transform[chan]);
		if (wavelet_transform_file(Imag, fit->ry, fit->rx, dir[chan], Type, Nbr_Plan,
				fit->pdata[chan])) {
			free((char *) Imag);
			free(dir[chan]);
			return 1;
		}
		if (wave_io_read(dir[chan], &Wavelet[chan])) {
			free((char *) Imag);
			free(dir[chan]);
			return 1;
		}
		Nl = Wavelet[chan].Nbr_Ligne;
		Nc = Wavelet[chan].Nbr_Col;
		pave_2d_extract_plan(Wavelet[chan].Pave.Data, Imag, Nl, Nc, Plan);
		reget_rawdata(Imag, Nl, Nc, fit->pdata[chan]);
		free(dir[chan]);
	}

	/* Free */
	if (Imag)
		free((char *) Imag);
	for (chan = start; chan < end; chan++) {
		wave_io_free(&Wavelet[chan]);
	}
	return 0;
}

/*****************************************************************************
 *                      M E D I A N     F I L T E R                          *
 ****************************************************************************/

gboolean end_median_filter(gpointer p) {
	struct median_filter_data *args = (struct median_filter_data *) p;
	stop_processing_thread();// can it be done here in case there is no thread?
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	update_used_memory();
	free(args);
	return FALSE;
}

/* The function smoothes an image using the median filter with the
 * ksize x ksize aperture. Each channel of a multi-channel image is 
 * processed independently. In-place operation is supported. */
gpointer median_filter(gpointer p) {
	struct median_filter_data *args = (struct median_filter_data *) p;
	assert(args->ksize % 2 == 1 && args->ksize > 1);
	int i, x, y, xx, yy, layer, iter = 0;
	int nx = args->fit->rx;
	int ny = args->fit->ry;
	int radius = (args->ksize - 1) / 2;
	double norm = (double) get_normalized_value(args->fit);

	assert(nx > 0 && ny > 0);

	struct timeval t_start, t_end;

	siril_log_color_message(_("Median Filter: processing...\n"), "red");
	gettimeofday(&t_start, NULL);

	do {
		if (args->iterations != 1)
			siril_log_message(_("Iteration #%d...\n"), iter + 1);
		for (layer = 0; layer < com.uniq->nb_layers; layer++) {
			/* FILL image upside-down */
			WORD **image = malloc(ny * sizeof(WORD *));
			if (image == NULL) {
				printf("median filter: error allocating data\n");
				gdk_threads_add_idle(end_median_filter, args);
				return GINT_TO_POINTER(1);
			}
			for (i = 0; i < ny; i++)
				image[ny - i - 1] = args->fit->pdata[layer] + i * nx;

			for (y = 0; y < ny; y++) {
				if (!get_thread_run())
					break;
				for (x = 0; x < nx; x++) {
					WORD *data = calloc(args->ksize * args->ksize,
							sizeof(WORD));
					if (data == NULL) {
						printf("median filter: error allocating data\n");
						free(image);
						gdk_threads_add_idle(end_median_filter, args);
						return GINT_TO_POINTER(1);
					}
					i = 0;
					for (yy = y - radius; yy <= y + radius; yy++) {
						for (xx = x - radius; xx <= x + radius; xx++) {
							WORD tmp;
							if (xx < 0 && yy >= 0) {
								if (yy >= ny)
									tmp = image[ny - 1][00];
								else
									tmp = image[yy][0];
							} else if (xx > 0 && yy <= 0) {
								if (xx >= nx)
									tmp = image[00][nx - 1];
								else
									tmp = image[0][xx];
							} else if (xx <= 0 && yy <= 0) {
								tmp = image[0][0];
							} else {
								if (xx >= nx && yy >= ny)
									tmp = image[ny - 1][nx - 1];
								else if (xx >= nx && yy < ny)
									tmp = image[yy][nx - 1];
								else if (xx < nx && yy >= ny)
									tmp = image[ny - 1][xx];
								else
									tmp = image[yy][xx];
							}
							data[i++] = tmp;
						}
					}
					quicksort_s(data, args->ksize * args->ksize);
					WORD median = round_to_WORD(gsl_stats_ushort_median_from_sorted_data(data, 1, args->ksize * args->ksize));
					double pixel = args->amount * (median / norm);
					pixel += (1.0 - args->amount)
							* ((double) image[y][x] / norm);
					image[y][x] = round_to_WORD(pixel * norm);
					free(data);
				}
			}
			free(image);
		}
		iter++;
	} while (iter < args->iterations && get_thread_run());
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	gdk_threads_add_idle(end_median_filter, args);

	return GINT_TO_POINTER(0);
}

static int fmul_layer(fits *a, int layer, float coeff) {
	WORD *buf;
	int i;

	if (coeff < 0.0)
		return 1;
	buf = a->pdata[layer];
	for (i = 0; i < a->rx * a->ry; ++i) {
		buf[i] = round_to_WORD(buf[i] * coeff);
	}
	return 0;
}

/*****************************************************************************
 *      B A N D I N G      R E D U C T I O N      M A N A G E M E N T        *
 ****************************************************************************/

int banding_image_hook(struct generic_seq_args *args, int i, fits *fit, rectangle *_) {
	struct banding_data *banding_args = (struct banding_data *)args->user;
	return BandingEngine(fit, banding_args->sigma, banding_args->amount,
			banding_args->protect_highlights, banding_args->applyRotation);
}

void apply_banding_to_sequence(struct banding_data *banding_args) {
	struct generic_seq_args *args = malloc(sizeof(struct generic_seq_args));
	args->seq = &com.seq;
	args->partial_image = FALSE;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = com.seq.selnum;
	args->prepare_hook = ser_prepare_hook;
	args->finalize_hook = ser_finalize_hook;
	args->save_hook = NULL;
	args->image_hook = banding_image_hook;
	args->idle_function = NULL;
	args->description = "Banding Reduction";
	args->has_output = TRUE;
	args->new_seq_prefix = banding_args->seqEntry;
	args->load_new_sequence = TRUE;
	args->force_ser_output = FALSE;
	args->user = banding_args;
	args->already_in_a_thread = FALSE;
	args->parallel = TRUE;

	banding_args->fit = NULL;	// not used here

	start_in_new_thread(generic_sequence_worker, args);
}

// idle function executed at the end of the BandingEngine processing
gboolean end_BandingEngine(gpointer p) {
	struct banding_data *args = (struct banding_data *) p;
	stop_processing_thread();// can it be done here in case there is no thread?
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	update_used_memory();
	free(args);
	return FALSE;
}

/*** Reduces Banding in Canon DSLR images.
 * This code come from CanonBandingReduction.js v0.9.1, a script of
 * PixInsight, originally written by Georg Viehoever and
 * distributed under the terms of the GNU General Public License ******/
gpointer BandingEngineThreaded(gpointer p) {
	struct banding_data *args = (struct banding_data *) p;
	struct timeval t_start, t_end;

	siril_log_color_message(_("Banding Reducing: processing...\n"), "red");
	gettimeofday(&t_start, NULL);

	int retval = BandingEngine(args->fit, args->sigma, args->amount, args->protect_highlights, args->applyRotation);

	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	gdk_threads_add_idle(end_BandingEngine, args);
	
	return GINT_TO_POINTER(retval);
}

int BandingEngine(fits *fit, double sigma, double amount, gboolean protect_highlights, gboolean applyRotation) {
	int chan, row, i;
	WORD *line, *fixline;
	double minimum = DBL_MAX, globalsigma = 0.0;
	fits *fiximage;
	double invsigma = 1.0 / sigma;

	fiximage = calloc(1, sizeof(fits));
	if (fiximage == NULL) {
		fprintf(stderr, "BandingEngine: error allocating data\n");
		return 1;
	}

	if (applyRotation) {
#ifdef HAVE_OPENCV
		cvRotateImage(fit, 90.0, -1, OPENCV_LINEAR);
#else
		siril_log_message(_("Rotation is only possible when Siril has been compiled with OpenCV support.\n"));
		free(fiximage);
		return 1;
#endif
	}

	new_fit_image(fiximage, fit->rx, fit->ry, fit->naxes[2]);

	for (chan = 0; chan < fit->naxes[2]; chan++) {
		imstats *stat = statistics(fit, chan, NULL, STATS_BASIC | STATS_MAD, STATS_ZERO_NULLCHECK);
		if (!stat) {
			siril_log_message(_("Error: no data computed.\n"));
			return 1;
		}
		double background = stat->median;
		double *rowvalue = calloc(fit->ry, sizeof(double));
		if (rowvalue == NULL) {
			fprintf(stderr, "BandingEngine: error allocating data\n");
			free(stat);
			return 1;
		}
		if (protect_highlights) {
			globalsigma = stat->mad * MAD_NORM;
		}
		for (row = 0; row < fit->ry; row++) {
			line = fit->pdata[chan] + row * fit->rx;
			WORD *cpyline = calloc(fit->rx, sizeof(WORD));
			if (cpyline == NULL) {
				fprintf(stderr, "BandingEngine: error allocating data\n");
				free(stat);
				free(rowvalue);
				return 1;
			}
			memcpy(cpyline, line, fit->rx * sizeof(WORD));
			int n = fit->rx;
			quicksort_s(cpyline, n);
			if (protect_highlights) {
				WORD reject = round_to_WORD(
						background + invsigma * globalsigma);
				for (i = fit->rx - 1; i >= 0; i--) {
					if (cpyline[i] < reject)
						break;
					n--;
				}
			}

			double median = gsl_stats_ushort_median_from_sorted_data(cpyline, 1, n);
			rowvalue[row] = background - median;
			minimum = min(minimum, rowvalue[row]);
			free(cpyline);
		}
		for (row = 0; row < fit->ry; row++) {
			fixline = fiximage->pdata[chan] + row * fiximage->rx;
			for (i = 0; i < fit->rx; i++)
				fixline[i] = round_to_WORD(rowvalue[row] - minimum);
		}
		free(rowvalue);
		free(stat);
	}
	for (chan = 0; chan < fit->naxes[2]; chan++)
		fmul_layer(fiximage, chan, amount);
	imoper(fit, fiximage, OPER_ADD);

	clearfits(fiximage);
#ifdef HAVE_OPENCV
	if (applyRotation)
		cvRotateImage(fit, -90.0, -1, OPENCV_LINEAR);
#endif
	return 0;
}

/*****************************************************************************
 *       N O I S E     C O M P U T A T I O N      M A N A G E M E N T        *
 ****************************************************************************/

/* Based on Jean-Luc Starck and Fionn Murtagh (1998), Automatic Noise
 * Estimation from the Multiresolution Support, Publications of the
 * Royal Astronomical Society of the Pacific, vol. 110, pp. 193–199.
 * slow algorithm. For now it is replaced by faster one. BUT, we need to keep it
 * in case we need it -. */
int backgroundnoise(fits* fit, double sigma[]) {
	int layer, k;
	fits *waveimage = calloc(1, sizeof(fits));

	if (waveimage == NULL) {
		fprintf(stderr, "backgroundnoise: error allocating data\n");
		return 1;
	}

	copyfits(fit, waveimage, CP_ALLOC | CP_FORMAT | CP_COPYA, 0);
#ifdef HAVE_OPENCV	// a bit faster
	cvComputeFinestScale(waveimage);
#else
	if (get_wavelet_layers(waveimage, 4, 0, TO_PAVE_BSPLINE, -1)) {
		siril_log_message(_("Siril cannot evaluate the noise in the image\n"));
		clearfits(waveimage);
		return 1;
	}
#endif

	for (layer = 0; layer < fit->naxes[2]; layer++) {
		imstats *stat = statistics(waveimage, layer, NULL, STATS_BASIC, STATS_ZERO_NULLCHECK);
		if (!stat) {
			siril_log_message(_("Error: no data computed.\n"));
			return 1;
		}
		double sigma0 = stat->sigma;
		double mean = stat->mean;
		double epsilon = 0.0;
		WORD lo, hi;
		WORD *buf = waveimage->pdata[layer];
		unsigned int i;
		unsigned int ndata = fit->rx * fit->ry;
		assert(ndata > 0);
		WORD *array1 = calloc(ndata, sizeof(WORD));
		WORD *array2 = calloc(ndata, sizeof(WORD));
		if (array1 == NULL || array2 == NULL) {
			printf("backgroundnoise: Error allocating data\n");
			if (array1)
				free(array1);
			if (array2)
				free(array2);
			free(stat);
			return 1;
		}
		WORD *set = array1, *subset = array2;
		memcpy(set, buf, ndata * sizeof(WORD));

		lo = round_to_WORD(LOW_BOUND * stat->normValue);
		hi = round_to_WORD(HIGH_BOUND * stat->normValue);

		sigma[layer] = sigma0;

		int n = 0;
		do {
			sigma0 = sigma[layer];
			for (i = 0, k = 0; i < ndata; i++) {
				if (set[i] >= lo && set[i] <= hi) {
					if (fabs(set[i] - mean) < 3.0 * sigma0) {
						subset[k++] = set[i];
					}
				}
			}
			ndata = k;
			sigma[layer] = gsl_stats_ushort_sd(subset, 1, ndata);
			set = subset;
			(set == array1) ? (subset = array2) : (subset = array1);
			if (ndata == 0) {
				free(array1);
				free(array2);
				free(stat);
				siril_log_message(_("backgroundnoise: Error, no data computed\n"));
				sigma[layer] = 0.0;
				return 1;
			}
			n++;
			epsilon = fabs(sigma[layer] - sigma0) / sigma[layer];
		} while (epsilon > EPSILON && n < MAX_ITER);
		sigma[layer] *= SIGMA_PER_FWHM; // normalization
		sigma[layer] /= 0.974; // correct for 2% systematic bias
		if (n == MAX_ITER)
			siril_log_message(_("backgroundnoise: does not converge\n"));
		free(array1);
		free(array2);
		free(stat);
	}
	clearfits(waveimage);

	return 0;
}

gboolean end_noise(gpointer p) {
	struct noise_data *args = (struct noise_data *) p;
	stop_processing_thread();// can it be done here in case there is no thread?
	int chan, nb_chan;
	struct timeval t_end;

	nb_chan = args->fit->naxes[2];

	double norm = (double) get_normalized_value(args->fit);
	for (chan = 0; chan < nb_chan; chan++)
		siril_log_message(
				_("Background noise value (channel: #%d): %0.3lf (%.3e)\n"), chan,
				args->bgnoise[chan], args->bgnoise[chan] / norm);
	set_cursor_waiting(FALSE);
	update_used_memory();
	if (args->verbose) {
		gettimeofday(&t_end, NULL);
		show_time(args->t_start, t_end);
	}
	free(args);
	return FALSE;
}

gpointer noise(gpointer p) {
	struct noise_data *args = (struct noise_data *) p;
	int  chan;

	if (args->verbose) {
		siril_log_color_message(_("Noise standard deviation: calculating...\n"),
				"red");
		gettimeofday(&args->t_start, NULL);
	}

/*	if (backgroundnoise(args->fit, args->bgnoise)) {
		gdk_threads_add_idle(end_noise, args);
		return GINT_TO_POINTER(1);
	}
	*/
	for (chan = 0; chan < args->fit->naxes[2]; chan++) {
		imstats *stat = statistics(args->fit, chan, NULL, STATS_BASIC, STATS_ZERO_NULLCHECK);
		if (!stat) {
			siril_log_message(_("Error: no data computed.\n"));
			gdk_threads_add_idle(end_noise, args);
			return GINT_TO_POINTER(1);
		}
		args->bgnoise[chan] = stat->bgnoise;
		free(stat);
	}

	gdk_threads_add_idle(end_noise, args);

	return GINT_TO_POINTER(0);
}


