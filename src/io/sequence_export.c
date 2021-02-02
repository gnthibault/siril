/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2021 team free-astro (see more in AUTHORS file)
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

#include <string.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_date.h"
#include "sequence.h"
#include "ser.h"
#include "stacking/stacking.h"
#include "registration/registration.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "io/image_format_fits.h"
#ifdef HAVE_FFMS2
#include "io/films.h"
#endif
#include "avi_pipp/avi_writer.h"
#ifdef HAVE_FFMPEG
#include "io/mp4_output.h"
#endif
#include "algos/geometry.h"

/* same order as in the combo box 'comboExport' */
typedef enum {
	EXPORT_FITS,
	EXPORT_FITSEQ,
	EXPORT_TIFF,
	EXPORT_SER,
	EXPORT_AVI,
	EXPORT_MP4,
	EXPORT_WEBM
} export_format;

struct exportseq_args {
	sequence *seq;
	seq_image_filter filtering_criterion;
	double filtering_parameter;

	char *basename;
	export_format output;
	gboolean normalize;

	int film_fps;		// has to be int for avi or ffmpeg
	int film_quality;	// [1, 5], for mp4 and webm

	gboolean resample;
	int32_t dest_width, dest_height;

	gboolean crop;
	rectangle crop_area;
};


/* Used for avi exporter, creates buffer as BGRBGR, using GUI cursors lo and hi */
static uint8_t *fits_to_uint8(fits *fit) {
	uint8_t *data;
	int w, h, i, j, channel, step;
	WORD lo, hi;

	w = fit->rx;
	h = fit->ry;
	channel = fit->naxes[2];
	step = (channel == 3 ? 2 : 0);
	float slope = compute_slope(&lo, &hi);

	data = malloc(w * h * channel * sizeof(uint8_t));
	for (i = 0, j = 0; i < w * h * channel; i += channel, j++) {
		float pixel = fit->pdata[RLAYER][j] - lo < 0 ? 0.0f : (float)(fit->pdata[RLAYER][j] - lo);
		data[i + step] = (uint8_t) roundf_to_BYTE(pixel * slope);
		if (channel > 1) {
			pixel = fit->pdata[GLAYER][j] - lo < 0 ? 0.0f : (float)(fit->pdata[GLAYER][j] - lo);
			data[i + 1] = (uint8_t) roundf_to_BYTE(pixel * slope);
			pixel = fit->pdata[BLAYER][j] - lo < 0 ? 0.0f : (float)(fit->pdata[BLAYER][j] - lo);
			data[i + 2 - step] = (uint8_t) roundf_to_BYTE(pixel * slope);
		}
	}
	return data;
}

static gpointer export_sequence(gpointer ptr) {
	int retval = 0, cur_nb = 0;
	unsigned int out_width, out_height, in_width, in_height;
	uint8_t *data;
	fits *destfit = NULL;	// must be declared before any goto!
	char filename[256], dest[256];
	struct ser_struct *ser_file = NULL;
	fitseq *fitseq_file = NULL;
	GSList *timestamp = NULL;
	char *filter_descr;
	GDateTime *strTime;
#ifdef HAVE_FFMPEG
	struct mp4_struct *mp4_file = NULL;
#endif
	struct exportseq_args *args = (struct exportseq_args *)ptr;
	norm_coeff coeff = { 0 };

	int reglayer = get_registration_layer(args->seq);
	siril_log_message(_("Using registration information from layer %d to export sequence\n"), reglayer);
	if (args->crop) {
		in_width  = args->crop_area.w;
		in_height = args->crop_area.h;
	} else {
		in_width  = args->seq->rx;
		in_height = args->seq->ry;
	}

	if (args->resample) {
		out_width = args->dest_width;
		out_height = args->dest_height;
		if (out_width == in_width && out_height == in_height)
			args->resample = FALSE;
	} else {
		out_width = in_width;
		out_height = in_height;
	}

	gboolean have_seqwriter = args->output == EXPORT_FITSEQ || args->output == EXPORT_SER;
	if (have_seqwriter)
		seqwriter_set_max_active_blocks(3);

	int output_bitpix = USHORT_IMG;

	/* possible output formats: FITS images, FITS cube, TIFF, SER, AVI, MP4, WEBM */
	// create the sequence file for single-file sequence formats
	switch (args->output) {
		case EXPORT_FITS:
			output_bitpix = args->seq->bitpix;
			break;
		case EXPORT_FITSEQ:
			fitseq_file = malloc(sizeof(fitseq));
			snprintf(dest, 256, "%s%s", args->basename, com.pref.ext);
			if (fitseq_create_file(dest, fitseq_file, -1)) {
				free(fitseq_file);
				fitseq_file = NULL;
				retval = -1;
				goto free_and_reset_progress_bar;
			}
			output_bitpix = args->seq->bitpix;
			break;
		case EXPORT_TIFF:
			output_bitpix = args->seq->bitpix;
			// limit to 16 bits for TIFF
			if (output_bitpix == FLOAT_IMG)
				output_bitpix = USHORT_IMG;
			break;

		case EXPORT_SER:
			ser_file = malloc(sizeof(struct ser_struct));
			snprintf(dest, 256, "%s.ser", args->basename);
			if (ser_create_file(dest, ser_file, TRUE, args->seq->ser_file)) {
				free(ser_file);
				ser_file = NULL;
				retval = -1;
				goto free_and_reset_progress_bar;
			}
			output_bitpix = args->seq->bitpix;
			if (output_bitpix == FLOAT_IMG)
				output_bitpix = USHORT_IMG;
			break;

		case EXPORT_AVI:
			snprintf(dest, 256, "%s.avi", args->basename);
			int32_t avi_format;

			if (args->seq->nb_layers == 1)
				avi_format = AVI_WRITER_INPUT_FORMAT_MONOCHROME;
			else avi_format = AVI_WRITER_INPUT_FORMAT_COLOUR;

			if (avi_file_create(dest, out_width, out_height, avi_format,
					AVI_WRITER_CODEC_DIB, args->film_fps)) {
				siril_log_color_message(_("AVI file `%s' could not be created\n"), "red", dest);
				retval = -1;
				goto free_and_reset_progress_bar;
			}

			output_bitpix = BYTE_IMG;
			break;

		case EXPORT_MP4:
		case EXPORT_WEBM:
#ifndef HAVE_FFMPEG
			siril_log_message(_("MP4 output is not supported because siril was not compiled with ffmpeg support.\n"));
			retval = -1;
			goto free_and_reset_progress_bar;
#else
			/* resampling is managed by libswscale */
			snprintf(dest, 256, "%s.%s", args->basename,
					args->output == EXPORT_MP4 ? "mp4" : "webm");

			if (in_width % 32 || out_height % 2 || out_width % 2) {
				siril_log_message(_("Film output needs to have a width that is a multiple of 32 and an even height, resizing selection.\n"));
				if (in_width % 32) in_width = (in_width / 32) * 32 + 32;
				if (in_height % 2) in_height++;
				if (args->crop) {
					args->crop_area.w = in_width;
					args->crop_area.h = in_height;
				} else {
					args->crop = TRUE;
					args->crop_area.x = 0;
					args->crop_area.y = 0;
					args->crop_area.w = in_width;
					args->crop_area.h = in_height;
				}
				compute_fitting_selection(&args->crop_area, 32, 2, 0);
				memcpy(&com.selection, &args->crop_area, sizeof(rectangle));
				fprintf(stdout, "final input area: %d,%d,\t%dx%d\n",
						args->crop_area.x, args->crop_area.y,
						args->crop_area.w, args->crop_area.h);
				in_width = args->crop_area.w;
				in_height = args->crop_area.h;
				if (!args->resample) {
					out_width = in_width;
					out_height = in_height;
				} else {
					if (out_width % 2) out_width++;
					if (out_height % 2) out_height++;
				}
			}

			mp4_file = mp4_create(dest, out_width, out_height, args->film_fps, args->seq->nb_layers, args->film_quality, in_width, in_height);
			if (!mp4_file) {
				retval = -1;
				goto free_and_reset_progress_bar;
			}
			output_bitpix = BYTE_IMG;
			break;
#endif
	}

	if (output_bitpix == FLOAT_IMG)
		output_bitpix = com.pref.force_to_16bit ? USHORT_IMG : FLOAT_IMG;

	int nb_frames = compute_nb_filtered_images(args->seq,
			args->filtering_criterion, args->filtering_parameter);
	filter_descr = describe_filter(args->seq, args->filtering_criterion,
			args->filtering_parameter);
	siril_log_message(filter_descr);
	g_free(filter_descr);

	if (args->normalize) {
		struct stacking_args stackargs = { 0 };
		stackargs.force_norm = FALSE;
		stackargs.seq = args->seq;
		stackargs.filtering_criterion = args->filtering_criterion;
		stackargs.filtering_parameter = args->filtering_parameter;
		stackargs.nb_images_to_stack = nb_frames; 
		stackargs.normalize = ADDITIVE_SCALING;
		stackargs.reglayer = reglayer;

		// build image indices used by normalization
		if (stack_fill_list_of_unfiltered_images(&stackargs))
			goto free_and_reset_progress_bar;

		do_normalization(&stackargs);
		coeff.offset = stackargs.coeff.offset;
		coeff.scale = stackargs.coeff.scale;
		free(stackargs.coeff.mul);

		// and dispose them, because we don't need them anymore
		free(stackargs.image_indices);
	}

	long naxes[3];
	size_t nbpix = 0;
	set_progress_bar_data(NULL, PROGRESS_RESET);

	for (int i = 0, skipped = 0; i < args->seq->number; ++i) {
		if (!get_thread_run()) {
			retval = -1;
			goto free_and_reset_progress_bar;
		}
		if (!args->filtering_criterion(args->seq, i, args->filtering_parameter)) {
			siril_log_message(_("image %d is excluded from export\n"), i);
			skipped++;
			continue;
		}

		if (have_seqwriter)
			seqwriter_wait_for_memory();

		if (!seq_get_image_filename(args->seq, i, filename)) {
			seqwriter_release_memory();
			retval = -1;
			goto free_and_reset_progress_bar;
		}
		gchar *tmpmsg = g_strdup_printf(_("Processing image %s"), filename);
		set_progress_bar_data(tmpmsg, (double)cur_nb / (double)nb_frames);
		g_free(tmpmsg);

		/* we read the full frame */
		fits fit = { 0 };
		if (seq_read_frame(args->seq, i, &fit, FALSE, -1)) {
			seqwriter_release_memory();
			siril_log_message(_("Export: could not read frame, aborting\n"));
			retval = -3;
			goto free_and_reset_progress_bar;
		}

		/* destfit is allocated to the full size. Data will be copied from fit,
		 * image buffers are duplicated. It will be cropped after the copy if
		 * needed */
		if (!nbpix || have_seqwriter) {
			if (!nbpix) {
				memcpy(naxes, fit.naxes, sizeof naxes);
				nbpix = fit.naxes[0] * fit.naxes[1];
			}
			else {
				if (memcmp(naxes, fit.naxes, sizeof naxes)) {
					fprintf(stderr, "An image of the sequence doesn't have the same dimensions\n");
					retval = -3;
					clearfits(&fit);
					seqwriter_release_memory();
					goto free_and_reset_progress_bar;
				}
			}
			data_type dest_type = output_bitpix == FLOAT_IMG ? DATA_FLOAT : DATA_USHORT;
			destfit = NULL;	// otherwise it's overwritten by new_fit_image()
			if (new_fit_image(&destfit, fit.rx, fit.ry, fit.naxes[2], dest_type)) {
				retval = -1;
				clearfits(&fit);
				seqwriter_release_memory();
				goto free_and_reset_progress_bar;
			}
			destfit->bitpix = output_bitpix;
			destfit->orig_bitpix = output_bitpix;
		}
		else if (memcmp(naxes, fit.naxes, sizeof naxes)) {
			fprintf(stderr, "An image of the sequence doesn't have the same dimensions\n");
			retval = -3;
			clearfits(&fit);
			seqwriter_release_memory();
			goto free_and_reset_progress_bar;
		}
		else {
			/* we don't have a seq writer or it's not the first frame, we can reuse destfit */
			if (destfit->type == DATA_FLOAT) {
				memset(destfit->fdata, 0, nbpix * fit.naxes[2] * sizeof(float));
				if (args->crop) {
					/* reset destfit damaged by the crop function */
					if (fit.naxes[2] == 3) {
						destfit->fpdata[1] = destfit->fdata + nbpix;
						destfit->fpdata[2] = destfit->fdata + nbpix * 2;
					}
					destfit->rx = destfit->naxes[0] = fit.rx;
					destfit->ry = destfit->naxes[1] = fit.ry;
				}
			} else {
				memset(destfit->data, 0, nbpix * fit.naxes[2] * sizeof(WORD));
				if (args->crop) {
					/* reset destfit damaged by the crop function */
					if (fit.naxes[2] == 3) {
						destfit->pdata[1] = destfit->data + nbpix;
						destfit->pdata[2] = destfit->data + nbpix * 2;
					}
					destfit->rx = destfit->naxes[0] = fit.rx;
					destfit->ry = destfit->naxes[1] = fit.ry;
				}
			}
		}
		/* we copy the header */
		copy_fits_metadata(&fit, destfit);

		int shiftx, shifty;
		/* load registration data for current image */
		if (reglayer != -1 && args->seq->regparam[reglayer]) {
			shiftx = roundf_to_int(args->seq->regparam[reglayer][i].shiftx);
			shifty = roundf_to_int(args->seq->regparam[reglayer][i].shifty);
		} else {
			shiftx = 0;
			shifty = 0;
		}

		/* fill the image with shifted data and normalization */
		for (int layer = 0; layer < fit.naxes[2]; ++layer) {
			for (int y = 0; y < fit.ry; ++y) {
				for (int x = 0; x < fit.rx; ++x) {
					int nx = x + shiftx;
					int ny = y + shifty;
					if (nx >= 0 && nx < fit.rx && ny >= 0 && ny < fit.ry) {
						if (fit.type == DATA_USHORT) {
							WORD pixel = fit.pdata[layer][x + y * fit.rx];
							if (args->normalize) {
								double tmp = (double) pixel;
								tmp *= coeff.scale[i];
								tmp -= coeff.offset[i];
								pixel = round_to_WORD(tmp);
							}
							destfit->pdata[layer][nx + ny * fit.rx] = pixel;
						}
						else if (fit.type == DATA_FLOAT) {
							float pixel = fit.fpdata[layer][x + y * fit.rx];
							if (args->normalize) {
								pixel *= (float) coeff.scale[i];
								pixel -= (float) coeff.offset[i];
							}
							if (destfit->type == DATA_FLOAT) {
								destfit->fpdata[layer][nx + ny * fit.rx] = pixel;
							} else {
								destfit->pdata[layer][nx + ny * fit.rx] = roundf_to_WORD(pixel * USHRT_MAX_SINGLE);
							}
						} else {
							retval = -1;
							clearfits(&fit);
							seqwriter_release_memory();
							goto free_and_reset_progress_bar;
						}
						// for 8 bit output, destfit will be transformed later with a linear scale
					}
				}
			}
		}
		clearfits(&fit);

		if (args->crop) {
			crop(destfit, &args->crop_area);
		}

		switch (args->output) {
			case EXPORT_FITS:
				snprintf(dest, 255, "%s%05d%s", args->basename, i + 1, com.pref.ext);
				retval = savefits(dest, destfit);
				break;
			case EXPORT_FITSEQ:
				retval = fitseq_write_image(fitseq_file, destfit, i - skipped);
				break;
#ifdef HAVE_LIBTIFF
			case EXPORT_TIFF:
				snprintf(dest, 255, "%s%05d", args->basename, i + 1);
				retval = savetif(dest, destfit, 16);
				break;
#endif
			case EXPORT_SER:
				if (destfit->date_obs) {
					strTime = g_date_time_ref(destfit->date_obs);
					timestamp = g_slist_append(timestamp, strTime);
				}
				retval = ser_write_frame_from_fit(ser_file, destfit, i - skipped);
				break;
			case EXPORT_AVI:
				data = fits_to_uint8(destfit);
				retval = avi_file_write_frame(0, data);
				break;
#ifdef HAVE_FFMPEG
			case EXPORT_MP4:
			case EXPORT_WEBM:
				// an equivalent to fits_to_uint8 is called in there (fill_rgb_image)...
				retval = mp4_add_frame(mp4_file, destfit);
				break;
#endif
		}
		if (retval) {
			seqwriter_release_memory();
			goto free_and_reset_progress_bar;
		}
		cur_nb++;
	}

free_and_reset_progress_bar:
	if (destfit && !have_seqwriter)
		clearfits(destfit);
	if (args->normalize) {
		free(coeff.offset);
		free(coeff.scale);
	}
	// close the sequence file for single-file sequence formats
	switch (args->output) {
		case EXPORT_FITSEQ:
			if (fitseq_file) {
				fitseq_close_file(fitseq_file);
				free(fitseq_file);
			}
			break;
		case EXPORT_SER:
			if (ser_file) {
				if (timestamp)
					ser_convertTimeStamp(ser_file, timestamp);
				ser_write_and_close(ser_file);
				free(ser_file);
			}
			g_slist_free_full(timestamp, (GDestroyNotify) g_date_time_unref);
			break;
		case EXPORT_AVI:
			avi_file_close(0);
			break;
		case EXPORT_MP4:
		case EXPORT_WEBM:
#ifdef HAVE_FFMPEG
			if (mp4_file) {
				mp4_close(mp4_file);
				free(mp4_file);
			}
#endif
			break;
		default:
			break;
	}

	if (retval) {
		set_progress_bar_data(_("Sequence export failed. Check the log."), PROGRESS_RESET);
		siril_log_message(_("Sequence export failed\n"));
	}
	else {
		set_progress_bar_data(_("Sequence export succeeded."), PROGRESS_RESET);
		siril_log_message(_("Sequence export succeeded.\n"));
	}

	free(args->basename);
	free(args);
	args = NULL;
	siril_add_idle(end_generic, args);
	return NULL;
}

void on_buttonExportSeq_clicked(GtkButton *button, gpointer user_data) {
	int selected = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("comboExport")));
	const char *bname = gtk_entry_get_text(GTK_ENTRY(lookup_widget("entryExportSeq")));
	gboolean normalize = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("exportNormalize")));
	struct exportseq_args *args;

	if (bname[0] == '\0') return;
	if (selected == -1) return;

	args = malloc(sizeof(struct exportseq_args));
	args->seq = &com.seq;
	get_sequence_filtering_from_gui(&args->filtering_criterion, &args->filtering_parameter);
	args->basename = g_str_to_ascii(bname, NULL);
	args->output = (export_format)selected;
	args->normalize = normalize;
	args->resample = FALSE;
	args->crop = com.selection.w && com.selection.h;
	if (args->crop)
		memcpy(&args->crop_area, &com.selection, sizeof(rectangle));

	if (args->output == EXPORT_AVI || args->output == EXPORT_MP4 || args->output == EXPORT_WEBM) {
		GtkEntry *fpsEntry = GTK_ENTRY(lookup_widget("entryAviFps"));
		args->film_fps = round_to_int(g_ascii_strtod(gtk_entry_get_text(fpsEntry), NULL));
		if (args->film_fps <= 0) args->film_fps = 1;
	}
	if (args->output == EXPORT_MP4 || args->output == EXPORT_WEBM) {
		GtkAdjustment *adjQual = GTK_ADJUSTMENT(gtk_builder_get_object(builder,"adjustment3"));
		GtkToggleButton *checkResize = GTK_TOGGLE_BUTTON(lookup_widget("checkAviResize"));
		args->film_quality = (int)gtk_adjustment_get_value(adjQual);
		args->resample = gtk_toggle_button_get_active(checkResize);
		if (args->resample) {
			GtkEntry *widthEntry = GTK_ENTRY(lookup_widget("entryAviWidth"));
			GtkEntry *heightEntry = GTK_ENTRY(lookup_widget("entryAviHeight"));
			args->dest_width = g_ascii_strtoll(gtk_entry_get_text(widthEntry), NULL, 10);
			args->dest_height = g_ascii_strtoll(gtk_entry_get_text(heightEntry), NULL, 10);
			if (args->dest_height == 0 || args->dest_width == 0) {
				siril_log_message(_("Width or height cannot be null. Not resizing.\n"));
				gtk_toggle_button_set_active(checkResize, FALSE);
				args->resample = FALSE;
			} else if (args->dest_height == args->seq->ry && args->dest_width == args->seq->rx) {
				gtk_toggle_button_set_active(checkResize, FALSE);
				args->resample = FALSE;
			}
		}
	}
	else if (args->output == EXPORT_FITS || args->output == EXPORT_TIFF) {
		// add a trailing '_' for multiple-files sequences
		args->basename = format_basename(args->basename, TRUE);
	}
	// Display a useful warning because I always forget to remove selection
	if (args->crop) {
		gboolean confirm = siril_confirm_dialog(_("Export cropped sequence?"),
				_("An active selection was detected. The exported sequence will only contain data within the drawn selection. You can confirm the crop or cancel it. "
						"If you choose to click on cancel, the exported sequence will contain all data."), _("Confirm Crop"));
		args->crop = confirm;
	}

	set_cursor_waiting(TRUE);
	start_in_new_thread(export_sequence, args);
}

void on_comboExport_changed(GtkComboBox *box, gpointer user_data) {
	GtkWidget *avi_options = lookup_widget("boxAviOptions");
	GtkWidget *checkAviResize = lookup_widget("checkAviResize");
	GtkWidget *quality = lookup_widget("exportQualScale");
	int output_type = gtk_combo_box_get_active(box);
	gtk_widget_set_visible(avi_options, output_type >= EXPORT_AVI);
	gtk_widget_set_visible(quality, output_type >= EXPORT_MP4);
	gtk_widget_set_sensitive(checkAviResize, output_type >= EXPORT_MP4);
}

void on_checkAviResize_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkWidget *heightEntry = lookup_widget("entryAviHeight");
	GtkWidget *widthEntry = lookup_widget("entryAviWidth");
	gtk_widget_set_sensitive(heightEntry, gtk_toggle_button_get_active(togglebutton));
	gtk_widget_set_sensitive(widthEntry, gtk_toggle_button_get_active(togglebutton));
}

void update_export_crop_label() {
	static GtkLabel *label = NULL;
	if (!label) 
		label = GTK_LABEL(lookup_widget("exportLabel"));
	if (com.selection.w && com.selection.h)
		gtk_label_set_text(label, _("Cropping to selection"));
	else gtk_label_set_text(label, _("Select area to crop"));
}

void on_entryExportSeq_changed(GtkEditable *editable, gpointer user_data){
	gchar *name = (gchar *)gtk_entry_get_text(GTK_ENTRY(editable));
	if (*name != 0) {
		if (check_if_seq_exist(name, !g_str_has_suffix(name, ".ser"))) {
			set_icon_entry(GTK_ENTRY(editable), "gtk-dialog-warning");
		} else {
			set_icon_entry(GTK_ENTRY(editable), NULL);
		}
	} else {
		set_icon_entry(GTK_ENTRY(editable), NULL);
	}
}

void on_entryAviHeight_changed(GtkEditable *editable, gpointer user_data);

void on_entryAviWidth_changed(GtkEditable *editable, gpointer user_data) {
	double ratio, width, height;
	gchar *c_height;
	GtkEntry *heightEntry = GTK_ENTRY(lookup_widget("entryAviHeight"));

	if (com.selection.w && com.selection.h) return;
	ratio = (double) com.seq.ry / (double) com.seq.rx;
	width = g_ascii_strtod(gtk_entry_get_text(GTK_ENTRY(editable)),NULL);
	height = ratio * width;
	c_height = g_strdup_printf("%d", (int)(height));

	g_signal_handlers_block_by_func(heightEntry, on_entryAviHeight_changed, NULL);
	gtk_entry_set_text(heightEntry, c_height);
	g_signal_handlers_unblock_by_func(heightEntry, on_entryAviHeight_changed, NULL);
	g_free(c_height);
}

void on_entryAviHeight_changed(GtkEditable *editable, gpointer user_data) {
	double ratio, width, height;
	gchar *c_width;
	GtkEntry *widthEntry = GTK_ENTRY(lookup_widget("entryAviWidth"));

	if (com.selection.w && com.selection.h) return;
	ratio = (double) com.seq.rx / (double) com.seq.ry;
	height = g_ascii_strtod(gtk_entry_get_text(GTK_ENTRY(editable)), NULL);
	width = ratio * height;
	c_width = g_strdup_printf("%d", (int)(width));

	g_signal_handlers_block_by_func(widthEntry, on_entryAviWidth_changed, NULL);
	gtk_entry_set_text(widthEntry, c_width);
	g_signal_handlers_unblock_by_func(widthEntry, on_entryAviWidth_changed, NULL);
	g_free(c_width);
}
