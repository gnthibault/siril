/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2019 team free-astro (see more in AUTHORS file)
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

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "algos/statistics.h"
#include "algos/cosmetic_correction.h"
#include "gui/progress_and_log.h"
#include "io/sequence.h"
#include "io/conversion.h"
#include "io/ser.h"

static double evaluateNoiseOfCalibratedImage(fits *fit, fits *dark, double k) {
	double noise = 0.0;
	fits dark_tmp = { 0 }, fit_tmp = { 0 };
	int chan, ret = 0;
	rectangle area = { 0 };

	/* square of 512x512 in the center of the image */
	int size = 512;
	area.x = (fit->rx - size) / 2;
	area.y = (fit->ry - size) / 2;
	area.w = size;
	area.h = size;

	copyfits(dark, &dark_tmp, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);
	copyfits(fit, &fit_tmp, CP_ALLOC | CP_COPYA | CP_FORMAT, -1);

	soper(&dark_tmp, k, OPER_MUL);
	ret = imoper(&fit_tmp, &dark_tmp, OPER_SUB);
	if (ret) {
		clearfits(&dark_tmp);
		clearfits(&fit_tmp);
		return -1.0;
	}

	for (chan = 0; chan < fit->naxes[2]; chan++) {
		imstats *stat = statistics(NULL, -1, &fit_tmp, chan, &area, STATS_BASIC);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 0.0;
		}
		noise += stat->sigma;
		free_stats(stat);
	}
	clearfits(&dark_tmp);
	clearfits(&fit_tmp);

	return noise;
}

// TODO: is this computed as double?
#define GR ((sqrt(5) - 1) / 2)

static double goldenSectionSearch(fits *brut, fits *dark, double a, double b,
		double tol) {
	double c, d;
	double fc, fd;
	int iter = 0;

	c = b - GR * (b - a);
	d = a + GR * (b - a);
	fc = evaluateNoiseOfCalibratedImage(brut, dark, c);
	fd = evaluateNoiseOfCalibratedImage(brut, dark, d);
	do {
		siril_debug_print("Iter: %d (%1.2lf, %1.2lf)\n", ++iter, c, d);
		if (fc < 0.0 || fd < 0.0)
			return -1.0;
		if (fc < fd) {
			b = d;
			d = c;
			fd = fc;
			c = b - GR * (b - a);
			fc = evaluateNoiseOfCalibratedImage(brut, dark, c);
		} else {
			a = c;
			c = d;
			fc = fd;
			d = a + GR * (b - a);
			fd = evaluateNoiseOfCalibratedImage(brut, dark, d);

		}
	} while (fabs(c - d) > tol);
	return ((b + a) / 2.0);
}


static int preprocess(fits *brut, fits *offset, fits *dark, fits *flat, float level) {
	int ret = 0;

	if (com.preprostatus & USE_OFFSET) {
		ret = imoper(brut, offset, OPER_SUB);
		if (ret)
			return ret;
	}

	/* if dark optimization, the master-dark has already been subtracted */
	if ((com.preprostatus & USE_DARK) && !(com.preprostatus & USE_OPTD)) {
		ret = imoper(brut, dark, OPER_SUB);
		if (ret)
			return ret;
	}

	if (com.preprostatus & USE_FLAT) {
		siril_fdiv(brut, flat, level);
	}

	return 0;
}

static int darkOptimization(fits *brut, fits *dark, fits *offset) {
	double k0;
	double lo = 0.0;
	double up = 2.0;
	int ret = 0;
	fits dark_tmp = { 0 };

	if (brut->rx != dark->rx || brut->ry != dark->ry) {
		return -1;
	}

	copyfits(dark, &dark_tmp, CP_ALLOC | CP_COPYA | CP_FORMAT, 0);

	/* Minimization of background noise to find better k */
	invalidate_stats_from_fit(brut);
	k0 = goldenSectionSearch(brut, &dark_tmp, lo, up, 1E-3);
	if (k0 < 0.0) {
		ret = -1;
	} else {
		siril_log_message(_("Dark optimization: k0=%.3lf\n"), k0);
		/* Multiply coefficient to master-dark */
		soper(&dark_tmp, k0, OPER_MUL);
		ret = imoper(brut, &dark_tmp, OPER_SUB);
	}
	clearfits(&dark_tmp);
	return ret;
}

static void prepro_cleanup(struct preprocessing_data *args) {
	free_sequence(args->seq, TRUE);
	free(args);
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
	if (args->is_sequence) {
		if (!args->retval) {
			// load the new sequence
			char *ppseqname = malloc(
					strlen(args->seq->ppprefix) + strlen(args->seq->seqname) + 5);
			sprintf(ppseqname, "%s%s.seq", args->seq->ppprefix,
					args->seq->seqname);
			check_seq(0);
			update_sequences_list(ppseqname);
			free(ppseqname);
		}
		sequence_free_preprocessing_data(args->seq);
		free(args->seq->ppprefix);
	}
#ifdef MAC_INTEGRATION
	GtkosxApplication *osx_app = gtkosx_application_get();
	gtkosx_application_attention_request(osx_app, INFO_REQUEST);
	g_object_unref (osx_app);
#endif
	free(args);
	return FALSE;
}

int preprocess_single_image(struct preprocessing_data *args) {
	char dest_filename[256], msg[256];
	fits *dark, *offset, *flat;
	fits fit = { 0 };
	int ret = 0;
	dark = args->dark;
	offset = args->offset;
	flat = args->flat;
	snprintf(msg, 255, _("Pre-processing image %s"), com.uniq->filename);
	msg[255] = '\0';
	set_progress_bar_data(msg, 0.5);

	copyfits(com.uniq->fit, &fit, CP_ALLOC | CP_FORMAT | CP_COPYA, 0);
	copy_fits_metadata(com.uniq->fit, &fit);

	if ((com.preprostatus & USE_OPTD) && (com.preprostatus & USE_DARK)) {
		ret = darkOptimization(&fit, dark, offset);
		if (ret) {
			set_progress_bar_data(msg, PROGRESS_NONE);
			clearfits(&fit);
			return 1;
		}
	}

	ret = preprocess(&fit, offset, dark, flat, args->normalisation);
	if (ret) {
		set_progress_bar_data(msg, PROGRESS_NONE);
		clearfits(&fit);
		return 1;
	}

	if ((com.preprostatus & USE_COSME) && (com.preprostatus & USE_DARK)) {
		if (dark->naxes[2] == 1) {
			/* Cosmetic correction */
			long icold, ihot;
			deviant_pixel *dev = find_deviant_pixels(dark, args->sigma, &icold, &ihot);
			siril_log_message(_("%ld pixels corrected (%ld + %ld)\n"),
					icold + ihot, icold, ihot);
			cosmeticCorrection(&fit, dev, icold + ihot, args->is_cfa);
			if (dev)
				free(dev);
		}
		else
			siril_log_message(_("Darkmap cosmetic correction "
						"is only supported with single channel images\n"));
	}

	if (args->debayer) {
		debayer_if_needed(TYPEFITS, &fit, args->compatibility, TRUE, args->stretch_cfa);
	}

	gchar *filename = g_path_get_basename(com.uniq->filename);
	char *filename_noext = remove_ext_from_filename(filename);
	snprintf(dest_filename, 255, "%s%s", com.uniq->ppprefix, filename_noext);
	dest_filename[255] = '\0';
	snprintf(msg, 255, _("Saving image %s"), filename_noext);
	msg[255] = '\0';
	set_progress_bar_data(msg, PROGRESS_NONE);
	args->retval = savefits(dest_filename, &fit);
	clearfits(&fit);
	g_free(filename);
	free(filename_noext);
	return 0;
}

/* doing the preprocessing. No unprotected GTK+ calls can go there.
 * returns 1 on error */
gpointer seqpreprocess(gpointer p) {
	char dest_filename[256], msg[256];
	fits *dark, *offset, *flat;
	int ret = 0;
	struct preprocessing_data *args = (struct preprocessing_data *) p;

	dark = args->dark;
	offset = args->offset;
	flat = args->flat;
	args->retval = 0;

	// remove old sequence
	if (args->is_sequence) {
		char *ppseqname = malloc(
				strlen(args->seq->ppprefix) + strlen(args->seq->seqname) + 5);
		sprintf(ppseqname, "%s%s.seq", args->seq->ppprefix, args->seq->seqname);
		unlink(ppseqname);
		free(ppseqname);
	}

	if (com.preprostatus & USE_FLAT) {
		if (args->equalize_cfa) {
			compute_grey_flat(flat);
		}
		if (args->autolevel) {
			imstats *stat = statistics(NULL, -1, flat, RLAYER, NULL, STATS_BASIC);
			if (!stat) {
				siril_log_message(_("Error: statistics computation failed.\n"));
				return GINT_TO_POINTER(1);
			}
			args->normalisation = stat->mean;
			siril_log_message(_("Normalisation value auto evaluated: %.2lf\n"),
					args->normalisation);
			free_stats(stat);
		}
	}

	if (!args->is_sequence) {
		preprocess_single_image(args);
	} else {	// sequence
		struct ser_struct *new_ser_file = NULL;
		char source_filename[256];
		int i;
		long icold = 0L, ihot = 0L;
		fits fit = { 0 };
		deviant_pixel *dev = NULL;

		int64_t size = seq_compute_size(args->seq, args->seq->number);
		if (args->debayer) size *= 3;
		if (test_available_space(size))
			return GINT_TO_POINTER(1);

		// creating a SER file if the input data is SER
		if (args->seq->type == SEQ_SER) {
			char new_ser_filename[256];
			new_ser_file = calloc(1, sizeof(struct ser_struct));
			snprintf(new_ser_filename, 255, "%s%s", args->seq->ppprefix, args->seq->ser_file->filename);
			if (ser_create_file(new_ser_filename, new_ser_file, TRUE, args->seq->ser_file)) {
				free(new_ser_file);
				new_ser_file = NULL;
				return GINT_TO_POINTER(1);
			}
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

		/* allocating memory to new fits */
		for (i = 0; i < args->seq->number; i++) {
			if (!get_thread_run())
				break;
			seq_get_image_filename(args->seq, i, source_filename);
			snprintf(msg, 255, _("Loading and pre-processing image %d/%d (%s)"),
					i + 1, args->seq->number, source_filename);
			msg[255] = '\0';
			set_progress_bar_data(msg,
					(double) (i + 1) / (double) args->seq->number);
			if (seq_read_frame(args->seq, i, &fit)) {
				snprintf(msg, 255, _("Could not read one of the raw files: %s."
						" Aborting preprocessing."), source_filename);
				msg[255] = '\0';
				set_progress_bar_data(msg, PROGRESS_RESET);
				args->retval = 1;
				break;
			}
			if ((com.preprostatus & USE_OPTD) && (com.preprostatus & USE_DARK)) {
				ret = darkOptimization(&fit, dark, offset);
				if (ret) {
					args->retval = 1;
					set_progress_bar_data(msg, PROGRESS_NONE);
					clearfits(&fit);
					break;
				}
			}

			ret = preprocess(&fit, offset, dark, flat, args->normalisation);
			if (ret) {
				args->retval = 1;
				set_progress_bar_data(msg, PROGRESS_NONE);
				clearfits(&fit);
				break;
			}

			if ((com.preprostatus & USE_COSME) && (com.preprostatus & USE_DARK) && (dark->naxes[2] == 1))
				cosmeticCorrection(&fit, dev, icold + ihot, args->is_cfa);

			if (args->debayer && args->seq->type == SEQ_REGULAR) {
				debayer_if_needed(TYPEFITS, &fit, args->compatibility, TRUE, args->stretch_cfa);
			}

			snprintf(dest_filename, 255, "%s%s", args->seq->ppprefix,
					source_filename);
			dest_filename[255] = '\0';
			snprintf(msg, 255, "Saving image %d/%d (%s)", i + 1, args->seq->number,
					dest_filename);
			if (args->seq->type == SEQ_SER) {
				args->retval = ser_write_frame_from_fit(new_ser_file, &fit, i);
			} else {
				args->retval = savefits(dest_filename, &fit);
			}
			clearfits(&fit);
			if (args->retval) {
				set_progress_bar_data(msg, PROGRESS_RESET);
				break;
			}
		}
		// closing SER file if it applies
		if (args->seq->type == SEQ_SER && (new_ser_file != NULL)) {
			ser_write_and_close(new_ser_file);
			free(new_ser_file);
			new_ser_file = NULL;
		}
		set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
		if (dev) free(dev);
	}
	if (!siril_add_idle(end_sequence_prepro, args))
		prepro_cleanup(args);
	return GINT_TO_POINTER(args->retval);
}

