#include <string.h>

#include "algos/geometry.h"
#include "core/siril.h"
#include "core/proto.h"
#include "sequence.h"
#include "ser.h"
#include "stacking/stacking.h"
#include "registration/registration.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#ifdef HAVE_FFMS2
#include "io/films.h"
#endif
#include "avi_pipp/avi_writer.h"
#include "opencv/opencv.h"
#ifdef HAVE_FFMPEG
#include "io/mp4_output.h"
#endif

struct exportseq_args {
	sequence *seq;
	char *basename;
	int convflags;
	gboolean normalize;
//	int gif_delay, gif_loops;
	double avi_fps;
	int quality;	// [1, 5], for mp4 and webm
	gboolean resize;
	int32_t dest_width, dest_height;
	gboolean crop;
	rectangle crop_area;
	seq_image_filter filtering_criterion;
	double filtering_parameter;
};

enum {
	EXPORT_FITS,
	EXPORT_TIFF,
	EXPORT_SER,
	EXPORT_AVI,
	EXPORT_MP4,
	EXPORT_WEBM
};

/* Used for avi exporter */
static uint8_t *fits_to_uint8(fits *fit) {
	uint8_t *data;
	int w, h, i, j, channel, step;
	float pente;
	WORD lo, hi;

	w = fit->rx;
	h = fit->ry;
	channel = fit->naxes[2];
	step = (channel == 3 ? 2 : 0);
	pente = computePente(&lo, &hi);

	data = malloc(w * h * channel * sizeof(uint8_t));
	for (i = 0, j = 0; i < w * h * channel; i += channel, j++) {
		data[i + step] = (uint8_t) round_to_BYTE(((double) fit->pdata[RLAYER][j] * pente));
		if (channel > 1) {
			data[i + 1] = (uint8_t) round_to_BYTE(((double) fit->pdata[GLAYER][j] * pente));
			data[i + 2 - step] = (uint8_t) round_to_BYTE(((double) fit->pdata[BLAYER][j] * pente));
		}
	}
	return data;
}

static gpointer export_sequence(gpointer ptr) {
	int i, x, y, nx, ny, shiftx, shifty, layer, retval = 0, reglayer,
	    nb_layers, skipped, nb_frames, cur_nb = 0;
	unsigned int out_width, out_height, in_width, in_height, nbdata = 0;
	uint8_t *data;
	fits fit = { 0 };
	fits destfit = { 0 };
	char filename[256], dest[256];
	struct ser_struct *ser_file = NULL;
	GSList *timestamp = NULL;
	char *filter_descr;
	gchar *strTime;
#ifdef HAVE_FFMPEG
	struct mp4_struct *mp4_file = NULL;
#endif
	struct exportseq_args *args = (struct exportseq_args *)ptr;
	norm_coeff coeff = { 0 };

	reglayer = get_registration_layer(args->seq);
	siril_log_message(_("Using registration information from layer %d to export sequence\n"), reglayer);
	if (args->crop) {
		in_width  = args->crop_area.w;
		in_height = args->crop_area.h;
	} else {
		in_width  = args->seq->rx;
		in_height = args->seq->ry;
	}

	if (args->resize) {
		out_width = args->dest_width;
		out_height = args->dest_height;
		if (out_width == in_width && out_height == in_height)
			args->resize = FALSE;
	} else {
		out_width = in_width;
		out_height = in_height;
	}

	switch (args->convflags) {
		case TYPESER:
			/* image size is not known here, no problem for crop or resize */
			ser_file = malloc(sizeof(struct ser_struct));
			snprintf(dest, 256, "%s.ser", args->basename);
			if (ser_create_file(dest, ser_file, TRUE, args->seq->ser_file))
				siril_log_message(_("Creating the SER file failed, aborting.\n"));
			break;

		case TYPEAVI:
			/* image size is set here, resize is managed by opencv when
			 * writing frames, we don't need crop size here */
			snprintf(dest, 256, "%s.avi", args->basename);
			int32_t avi_format;

			if (args->seq->nb_layers == 1)
				avi_format = AVI_WRITER_INPUT_FORMAT_MONOCHROME;
			else avi_format = AVI_WRITER_INPUT_FORMAT_COLOUR;

			avi_file_create(dest, out_width, out_height, avi_format,
					AVI_WRITER_CODEC_DIB, args->avi_fps);
			break;

		case TYPEMP4:
		case TYPEWEBM:
#ifndef HAVE_FFMPEG
			siril_log_message(_("MP4 output is not supported because siril was not compiled with ffmpeg support.\n"));
			retval = -1;
			goto free_and_reset_progress_bar;
#else
			/* image size is set here, resize is managed by ffmpeg so it also
			 * needs to know the input image size after crop */
			snprintf(dest, 256, "%s.%s", args->basename,
					args->convflags == TYPEMP4 ? "mp4" : "webm");
			if (args->avi_fps <= 0) args->avi_fps = 25;

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
				if (!args->resize) {
					out_width = in_width;
					out_height = in_height;
				} else {
					if (out_width % 2) out_width++;
					if (out_height % 2) out_height++;
				}
			}

			mp4_file = mp4_create(dest, out_width, out_height, args->avi_fps, args->seq->nb_layers, args->quality, in_width, in_height);
			if (!mp4_file) {
				retval = -1;
				goto free_and_reset_progress_bar;
			}
#endif
			break;
	}

	nb_frames = compute_nb_filtered_images(args->seq,
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

	set_progress_bar_data(NULL, PROGRESS_RESET);
	for (i = 0, skipped = 0; i < args->seq->number; ++i) {
		if (!get_thread_run()) {
			retval = -1;
			goto free_and_reset_progress_bar;
		}
		if (!args->filtering_criterion(args->seq, i, args->filtering_parameter)) {
			siril_log_message(_("image %d is excluded from export\n"), i);
			skipped++;
			continue;
		}

		if (!seq_get_image_filename(args->seq, i, filename)) {
			retval = -1;
			goto free_and_reset_progress_bar;
		}
		char *tmpmsg = strdup(_("Processing image "));
		tmpmsg = str_append(&tmpmsg, filename);
		set_progress_bar_data(tmpmsg, (double)cur_nb / (double)nb_frames);
		free(tmpmsg);

		if (seq_read_frame(args->seq, i, &fit)) {
			siril_log_message(_("Export: could not read frame, aborting\n"));
			retval = -3;
			goto free_and_reset_progress_bar;
		}

		if (!nbdata) {
			/* destfit is allocated to the real size because of the possible
			 * shifts and of the inplace cropping. FITS data is copied from
			 * fit, image buffers are duplicated. */
			memcpy(&destfit, &fit, sizeof(fits));
			destfit.header = NULL;
			destfit.fptr = NULL;
			nbdata = fit.rx * fit.ry;
			destfit.data = calloc(nbdata * fit.naxes[2], sizeof(WORD));
			destfit.stats = NULL;
			if (!destfit.data) {
				PRINT_ALLOC_ERR;
				retval = -1;
				goto free_and_reset_progress_bar;
			}
			destfit.pdata[0] = destfit.data;
			if (fit.naxes[2] == 1) {
				destfit.pdata[1] = destfit.data;
				destfit.pdata[2] = destfit.data;
			} else {
				destfit.pdata[1] = destfit.data + nbdata;
				destfit.pdata[2] = destfit.data + nbdata * 2;
			}
			nb_layers = fit.naxes[2];
		}
		else if (fit.ry * fit.rx != nbdata || nb_layers != fit.naxes[2]) {
			fprintf(stderr, "An image of the sequence doesn't have the same dimensions\n");
			retval = -3;
			goto free_and_reset_progress_bar;
		}
		else {
			/* we want copy the header */
			// TODO: why not use copyfits here?
			copy_fits_metadata(&fit, &destfit);
			memset(destfit.data, 0, nbdata * fit.naxes[2] * sizeof(WORD));
			if (args->crop) {
				/* reset destfit damaged by the crop function */
				if (fit.naxes[2] == 3) {
					destfit.pdata[1] = destfit.data + nbdata;
					destfit.pdata[2] = destfit.data + nbdata * 2;
				}
				destfit.rx = destfit.naxes[0] = fit.rx;
				destfit.ry = destfit.naxes[1] = fit.ry;
			}
		}

		/* load registration data for current image */
		if (reglayer != -1 && args->seq->regparam[reglayer]) {
			shiftx = roundf_to_int(args->seq->regparam[reglayer][i].shiftx);
			shifty = roundf_to_int(args->seq->regparam[reglayer][i].shifty);
		} else {
			shiftx = 0;
			shifty = 0;
		}

		/* fill the image with shift data and normalization */
		for (layer=0; layer<fit.naxes[2]; ++layer) {
			for (y=0; y < fit.ry; ++y){
				for (x=0; x < fit.rx; ++x){
					nx = x + shiftx;
					ny = y + shifty;
					if (nx >= 0 && nx < fit.rx && ny >= 0 && ny < fit.ry) {
						if (args->normalize) {
							double tmp = fit.pdata[layer][x + y * fit.rx];
							tmp *= coeff.scale[i];
							tmp -= coeff.offset[i];
							destfit.pdata[layer][nx + ny * fit.rx] = round_to_WORD(tmp);
						} else {
							destfit.pdata[layer][nx + ny * fit.rx] = fit.pdata[layer][x + y * fit.rx];
						}
					}
				}
			}
		}

		if (args->crop) {
			crop(&destfit, &args->crop_area);
		}

		switch (args->convflags) {
			case TYPEFITS:
				snprintf(dest, 255, "%s%05d", args->basename, i);
				if (savefits(dest, &destfit)) {
					retval = -1;
					goto free_and_reset_progress_bar;
				}
				break;
			case TYPETIFF:
				snprintf(dest, 255, "%s%05d%s", args->basename, i, com.ext);
				if (savetif(dest, &destfit, 16)) {
					retval = -1;
					goto free_and_reset_progress_bar;
				}

				break;
			case TYPESER:
				strTime = strdup(destfit.date_obs);
				timestamp = g_slist_append (timestamp, strTime);
				if (ser_write_frame_from_fit(ser_file, &destfit, i - skipped))
					siril_log_message(
							_("Error while converting to SER (no space left?)\n"));
				break;
			case TYPEAVI:
				data = fits_to_uint8(&destfit);

				if (args->resize) {
					uint8_t *newdata = malloc(out_width * out_height * destfit.naxes[2]);
					cvResizeGaussian_data8(data, destfit.rx, destfit.ry, newdata,
							out_width, out_height, destfit.naxes[2], OPENCV_CUBIC);
					avi_file_write_frame(0, newdata);
					free(newdata);
				}
				else
					avi_file_write_frame(0, data);
				free(data);
				break;
#ifdef HAVE_FFMPEG
			case TYPEMP4:
			case TYPEWEBM:
				mp4_add_frame(mp4_file, &destfit);
				break;
#endif
		}
		cur_nb++;
		clearfits(&fit);
	}

free_and_reset_progress_bar:
	clearfits(&fit);	// in case of goto
	free(destfit.data);
	if (args->normalize) {
		free(coeff.offset);
		free(coeff.scale);
	}
	if (args->convflags == TYPESER) {
		ser_convertTimeStamp(ser_file, timestamp);
		ser_write_and_close(ser_file);
		free(ser_file);
		g_slist_free_full(timestamp, g_free);
	}
	else if (args->convflags == TYPEAVI) {
		avi_file_close(0);
	}
#ifdef HAVE_FFMPEG
	else if (mp4_file && (args->convflags == TYPEMP4 || args->convflags == TYPEWEBM)) {
		mp4_close(mp4_file);
		free(mp4_file);
	}
#endif

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
	struct exportseq_args *args;
	GtkToggleButton *exportNormalize, *checkResize;
	GtkEntry *fpsEntry, *widthEntry, *heightEntry;
	GtkAdjustment *adjQual;

	if (bname[0] == '\0') return;
	if (selected == -1) return;

	args = malloc(sizeof(struct exportseq_args));
	args->basename = g_str_to_ascii(bname, NULL);
	args->seq = &com.seq;
	exportNormalize = GTK_TOGGLE_BUTTON(lookup_widget("exportNormalize"));
	args->normalize = gtk_toggle_button_get_active(exportNormalize);
	args->crop = com.selection.w && com.selection.h;
	args->resize = FALSE;
	if (args->crop)
		memcpy(&args->crop_area, &com.selection, sizeof(rectangle));

	// filtering
	get_sequence_filtering_from_gui(&args->filtering_criterion, &args->filtering_parameter);

	// format
	switch (selected) {
	case EXPORT_FITS:
		args->convflags = TYPEFITS;
		args->basename = format_basename(args->basename);
		break;
	case EXPORT_TIFF:
		args->convflags = TYPETIFF;
		args->basename = format_basename(args->basename);
		break;
	case EXPORT_SER:
		args->convflags = TYPESER;
		break;
	case EXPORT_AVI:
	case EXPORT_MP4:
	case EXPORT_WEBM:
		fpsEntry = GTK_ENTRY(lookup_widget("entryAviFps"));
		args->avi_fps = atoi(gtk_entry_get_text(fpsEntry));
		widthEntry = GTK_ENTRY(lookup_widget("entryAviWidth"));
		args->dest_width = atof(gtk_entry_get_text(widthEntry));
		heightEntry = GTK_ENTRY(lookup_widget("entryAviHeight"));
		args->dest_height = atof(gtk_entry_get_text(heightEntry));
		checkResize = GTK_TOGGLE_BUTTON(lookup_widget("checkAviResize"));
		adjQual = GTK_ADJUSTMENT(gtk_builder_get_object(builder,"adjustment3"));
		args->quality = (int)gtk_adjustment_get_value(adjQual);

		if (args->dest_height == 0 || args->dest_width == 0) {
			siril_log_message(_("Width or height cannot be null. Not resizing.\n"));
			args->resize = FALSE;
			gtk_toggle_button_set_active(checkResize, FALSE);
		} else if (args->dest_height == args->seq->ry && args->dest_width == args->seq->rx) {
			args->resize = FALSE;
			gtk_toggle_button_set_active(checkResize, FALSE);
		} else {
			args->resize = gtk_toggle_button_get_active(checkResize);
		}
		args->convflags = TYPEAVI;
		if (selected == EXPORT_MP4)
			args->convflags = TYPEMP4;
		else if (selected == EXPORT_WEBM)
			args->convflags = TYPEWEBM;
		break;
	default:
		free(args);
		return;
	}
	set_cursor_waiting(TRUE);
	start_in_new_thread(export_sequence, args);
}

void on_comboExport_changed(GtkComboBox *box, gpointer user_data) {
	GtkWidget *avi_options = lookup_widget("boxAviOptions");
	GtkWidget *checkAviResize = lookup_widget("checkAviResize");
	GtkWidget *quality = lookup_widget("exportQualScale");
	gtk_widget_set_visible(avi_options, gtk_combo_box_get_active(box) >= 3);
	gtk_widget_set_visible(quality, gtk_combo_box_get_active(box) >= 4);
	gtk_widget_set_sensitive(checkAviResize, TRUE);

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

		const char *ext = get_filename_ext(name);
		if (ext && !g_ascii_strcasecmp(ext, "ser")) {
			if (check_if_seq_exist(name)) {
				set_icon_entry(GTK_ENTRY(editable), "gtk-dialog-warning");
			} else{
				set_icon_entry(GTK_ENTRY(editable), NULL);
			}
		} else {
			if (check_if_seq_exist(name)) {
				set_icon_entry(GTK_ENTRY(editable), "gtk-dialog-warning");
			} else{
				set_icon_entry(GTK_ENTRY(editable), NULL);
			}
		}

	} else {
		set_icon_entry(GTK_ENTRY(editable), NULL);
	}
}
