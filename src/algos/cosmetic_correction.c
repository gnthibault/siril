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

#include <stdio.h>
#include <string.h>
#include <gsl/gsl_statistics_ushort.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "gui/callbacks.h"
#include "io/single_image.h"
#include "io/ser.h"
#include "algos/cosmetic_correction.h"


static WORD getMedian5x5(WORD *buf, const int xx, const int yy, const int w,
		const int h, gboolean is_cfa) {
	int step, radius, x, y;
	WORD *value, median;

	if (is_cfa) {
		step = 2;
		radius = 4;
	}
	else {
		step = 1;
		radius = 2;
	}

	int n = 0;
	int start;
	value = calloc(24, sizeof(WORD));
	for (y = yy - radius; y <= yy + radius; y += step) {
		for (x = xx - radius; x <= xx + radius; x += step) {
			if (y >= 0 && y < h) {
				if (x >= 0 && x < w) {
					if ((x != xx) || (y != yy)) {
						value[n++] = buf[x + y * w];
					}
				}
			}
		}
	}
	start = 24 - n - 1;
	quicksort_s(value, 24);
	median = round_to_WORD(gsl_stats_ushort_median_from_sorted_data(value + start, 1, n));
	free(value);
	return median;
}


static WORD *getAverage3x3Line(WORD *buf, const int yy, const int w, const int h,
		gboolean is_cfa) {
	int step, radius, x, xx, y;
	WORD *cpyline;

	if (is_cfa)
		step = radius = 2;
	else
		step = radius = 1;

	cpyline = calloc(w, sizeof(WORD));
	for (xx = 0; xx < w; ++xx) {
		int n = 0;
		double value = 0;
		for (y = yy - radius; y <= yy + radius; y += step) {
			if (y != yy) {	// we skip the line
				for (x = xx - radius; x <= xx + radius; x += step) {
					if (y >= 0 && y < h) {
						if (x >= 0 && x < w) {
							value += (double) buf[x + y * w];
							n++;
						}
					}
				}
			}
		}
		cpyline[xx] = round_to_WORD(value / n);
	}
	return cpyline;
}


static WORD getAverage3x3(WORD *buf, const int xx, const int yy, const int w,
		const int h, gboolean is_cfa) {
	int step, radius, x, y;
	double value = 0;

	if (is_cfa)
		step = radius = 2;
	else
		step = radius = 1;

	int n = 0;
	for (y = yy - radius; y <= yy + radius; y += step) {
		for (x = xx - radius; x <= xx + radius; x += step) {
			if (y >= 0 && y < h) {
				if (x >= 0 && x < w) {
					if ((x != xx) || (y != yy)) {
						value += (double) buf[x + y * w];
						n++;
					}
				}
			}
		}
	}
	return round_to_WORD(value / n);
}

long count_deviant_pixels(fits *fit, double sig[2], long *icold, long *ihot) {
	int i;
	WORD *buf = fit->pdata[RLAYER];
	imstats *stat;
	double sigma, median, thresHot, thresCold;

	/** statistics **/
	stat = statistics(fit, RLAYER, NULL, STATS_BASIC, STATS_ZERO_NULLCHECK);
	if (!stat) {
		siril_log_message(_("Error: no data computed.\n"));
		return 0L;
	}
	sigma = stat->sigma;
	median = stat->median;

	if (sig[0] == -1.0) {	// flag for no cold detection
		thresCold = -1.0;
	}
	else {
		double val = median - (sig[0] * sigma);
		thresCold = (val > 0) ? val : 0.0;
	}
	if (sig[1] == -1.0) {	// flag for no hot detection
		thresHot = USHRT_MAX_DOUBLE + 1;
	}
	else {
		double val = median + (sig[1] * sigma);
		thresHot = (val > USHRT_MAX_DOUBLE) ? USHRT_MAX_DOUBLE : val;
	}

	free(stat);

	/** We count deviant pixels **/
	*icold = 0;
	*ihot = 0;
	for (i = 0; i < fit->rx * fit->ry; i++) {
		if (buf[i] >= thresHot) (*ihot)++;
		else if (buf[i] <= thresCold) (*icold)++;
	}

	return (*icold + *ihot);
}


/* Gives a list of point p containing deviant pixel coordinates
 * p MUST be freed after the call
 * if cold == -1 or hot == -1, this is a flag to not compute cold or hot
 */
deviant_pixel *find_deviant_pixels(fits *fit, double sig[2], long *icold, long *ihot) {
	int x, y, i;
	WORD *buf = fit->pdata[RLAYER];
	imstats *stat;
	double sigma, median, thresHot, thresCold;
	deviant_pixel *dev;

	/** statistics **/
	stat = statistics(fit, RLAYER, NULL, STATS_BASIC, STATS_ZERO_NULLCHECK);
	if (!stat) {
		siril_log_message(_("Error: no data computed.\n"));
		return NULL;
	}
	sigma = stat->sigma;
	median = stat->median;

	if (sig[0] == -1.0) {	// flag for no cold detection
		thresCold = -1.0;
	}
	else {
		double val = median - (sig[0] * sigma);
		thresCold = (val > 0) ? val : 0.0;
	}
	if (sig[1] == -1.0) {	// flag for no hot detection
		thresHot = USHRT_MAX_DOUBLE + 1;
	}
	else {
		double val = median + (sig[1] * sigma);
		thresHot = (val > USHRT_MAX_DOUBLE) ? USHRT_MAX_DOUBLE : val;
	}

	free(stat);

	/** First we count deviant pixels **/
	*icold = 0;
	*ihot = 0;
	for (i = 0; i < fit->rx * fit->ry; i++) {
		if (buf[i] >= thresHot) (*ihot)++;
		else if (buf[i] <= thresCold) (*icold)++;
	}

	/** Second we store deviant pixels in p*/
	int n = (*icold) + (*ihot);
	if (n <= 0) return NULL;
	dev = calloc(n, sizeof(deviant_pixel));
	i = 0;
	for (y = 0; y < fit->ry; y++) {
		for (x = 0; x < fit->rx; x++) {
			double pixel = (double) buf[x + y * fit->rx];
			if (pixel >= thresHot) {
				dev[i].p.x = x;
				dev[i].p.y = y;
				dev[i].type = HOT_PIXEL;
				i++;
			}
			else if (pixel <= thresCold) {
				dev[i].p.x = x;
				dev[i].p.y = y;
				dev[i].type = COLD_PIXEL;
				i++;
			}
		}
	}
	return dev;
}

int cosmeticCorrOnePoint(fits *fit, deviant_pixel dev, gboolean is_cfa) {
	WORD *buf = fit->pdata[RLAYER];		// Cosmetic correction, as developed here, is only used on 1-channel images
	WORD newpixel;
	int width = fit->rx;
	int height = fit->ry;
	int x = (int) dev.p.x;
	int y = (int) dev.p.y;

	if (dev.type == COLD_PIXEL)
		newpixel = getMedian5x5(buf, x, y, width, height, is_cfa);
	else
		newpixel = getAverage3x3(buf, x, y, width, height, is_cfa);

	buf[x + y * fit->rx] = newpixel;
	return 0;
}

int cosmeticCorrOneLine(fits *fit, deviant_pixel dev, gboolean is_cfa) {
	WORD *buf = fit->pdata[RLAYER];
	WORD *line, *newline;
	int width = fit->rx;
	int height = fit->ry;
	int row = (int) dev.p.y;

	line = buf + row * width;
	newline = getAverage3x3Line(buf, row, width, height, is_cfa);
	memcpy(line, newline, width * sizeof(WORD));

	free(newline);

	return 0;
}

int cosmeticCorrection(fits *fit, deviant_pixel *dev, int size, gboolean is_cfa) {
	int i;
	WORD *buf = fit->pdata[RLAYER];		// Cosmetic correction, as developed here, is only used on 1-channel images
	int width = fit->rx;
	int height = fit->ry;

	for (i = 0; i < size; i++) {
		WORD newPixel;
		int xx = (int) dev[i].p.x;
		int yy = (int) dev[i].p.y;

		if (dev[i].type == COLD_PIXEL)
			newPixel = getMedian5x5(buf, xx, yy, width, height, is_cfa);
		else
			newPixel = getAverage3x3(buf, xx, yy, width, height, is_cfa);

		buf[xx + yy * width] = newPixel;
	}
	return 0;
}

/**** Autodetect *****/
int cosmetic_image_hook(struct generic_seq_args *args, int i, fits *fit, rectangle *_) {
	struct cosmetic_data *c_args = (struct cosmetic_data *) args->user;
	int retval, chan;
	/* Count variables, icold and ihot, need to be local in order to be parallelized */
	long icold, ihot;

	icold = ihot = 0L;
	for (chan = 0; chan < fit->naxes[2]; chan++) {
		retval = autoDetect(fit, chan, c_args->sigma, &icold, &ihot,
				c_args->amount, c_args->is_cfa);
		if (retval)
			return retval;
	}
	siril_log_color_message(_("Image %d: %ld pixels corrected (%ld + %ld)\n"),
			"bold", i, icold + ihot, icold, ihot);
	return 0;
}

void apply_cosmetic_to_sequence(struct cosmetic_data *cosme_args) {
	struct generic_seq_args *args = malloc(sizeof(struct generic_seq_args));
	args->seq = &com.seq;
	args->partial_image = FALSE;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = com.seq.selnum;
	args->prepare_hook = ser_prepare_hook;
	args->finalize_hook = ser_finalize_hook;
	args->save_hook = NULL;
	args->image_hook = cosmetic_image_hook;
	args->idle_function = NULL;
	args->description = "Cosmetic Correction";
	args->has_output = TRUE;
	args->new_seq_prefix = cosme_args->seqEntry;
	args->load_new_sequence = TRUE;
	args->force_ser_output = FALSE;
	args->user = cosme_args;
	args->already_in_a_thread = FALSE;
	args->parallel = TRUE;

	cosme_args->fit = NULL;	// not used here

	start_in_new_thread(generic_sequence_worker, args);
}

// idle function executed at the end of the BandingEngine processing
gboolean end_autoDetect(gpointer p) {
	struct cosmetic_data *args = (struct cosmetic_data *) p;

	stop_processing_thread();// can it be done here in case there is no thread?
	siril_log_message(_("%ld pixels corrected (%ld + %ld)\n"),
			args->icold + args->ihot, args->icold, args->ihot);
	set_cursor_waiting(FALSE);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	update_used_memory();
	free(args);
	return FALSE;
}

gpointer autoDetectThreaded(gpointer p) {
	struct cosmetic_data *args = (struct cosmetic_data *) p;
	struct timeval t_start, t_end;
	int retval = 0, chan;
	long icold, ihot;

	siril_log_color_message(_("Cosmetic Correction: processing...\n"), "red");
	gettimeofday(&t_start, NULL);

	icold = ihot = 0L;
	for (chan = 0; chan < args->fit->naxes[2]; chan++) {
		retval = autoDetect(args->fit, chan, args->sigma, &icold, &ihot,
				args->amount, args->is_cfa);
		if (retval)
			break;
	}
	args->icold = icold;
	args->ihot = ihot;
	gettimeofday(&t_end, NULL);
	show_time(t_start, t_end);
	gdk_threads_add_idle(end_autoDetect, args);

	return GINT_TO_POINTER(retval);
}

/* this is an autodetect algorithm. Cold and hot pixels
 *  are corrected in the same time */
int autoDetect(fits *fit, int layer, double sig[2], long *icold, long *ihot, double amount,
		gboolean is_cfa) {
	int x, y;
	int width = fit->rx;
	int height = fit->ry;
	double bkg, avgDev;
	double f0 = amount;
	double f1 = 1 - f0;
	imstats *stat;

	/* XXX: if cfa, stats are irrelevant. We should compute them taking
	 * into account the Bayer pattern */
	stat = statistics(fit, layer, NULL, STATS_BASIC | STATS_AVGDEV, STATS_ZERO_NULLCHECK);
	if (!stat) {
		siril_log_message(_("Error: no data computed.\n"));
		return 1;
	}
	bkg = stat->median;
	avgDev = stat->avgDev;
	free(stat);

	WORD *buf = fit->pdata[layer];

	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			WORD pixel = buf[x + y * width];
			WORD a = getAverage3x3(buf, x, y, width, height, is_cfa);
			WORD m = getMedian5x5(buf, x, y, width, height, is_cfa);

			/* Hot autodetect */
			if (sig[1] != -1.0) {
				double k1 = avgDev;
				double k2 = k1 / 2;
				double k3 = sig[1] * k1;
				if ((a < bkg + k2) && (pixel > bkg + k1) && (pixel > m + k3)) {
					(*ihot)++;
					buf[x + y * width] = a * f0 + pixel * f1;
				}
			}

			/* Cold autodetect */
			if (sig[0] != -1.0) {
				double k = avgDev * sig[0];
				if (((pixel + k) < bkg) && ((pixel + k) < m)) {
					(*icold)++;
					buf[x + y * width] = m * f0 + pixel * f1;
				}
			}
		}
	}
	return 0;
}
