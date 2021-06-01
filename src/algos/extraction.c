#include <string.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/command.h"
#include "algos/demosaicing.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "gui/dialogs.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "extraction.h"

static void update_sampling_information(fits *fit) {
	clear_Bayer_information(fit);

	fit->pixel_size_x *= 2;
	fit->pixel_size_y *= 2;
}

int extractHa_ushort(fits *in, fits *Ha, sensor_pattern pattern) {
	int width = in->rx / 2, height = in->ry / 2;

	if (strlen(in->bayer_pattern) > 4) {
		siril_log_message(_("Extract_Ha does not work on non-Bayer filter camera images!\n"));
		return 1;
	}
	if (new_fit_image(&Ha, width, height, 1, DATA_USHORT)) {
		return 1;
	}

	int j = 0;

	for (int row = 0; row < in->ry - 1; row += 2) {
		for (int col = 0; col < in->rx - 1; col += 2) {
			WORD c0 = in->data[col + row * in->rx];
			WORD c1 = in->data[1 + col + row * in->rx];
			WORD c2 = in->data[col + (1 + row) * in->rx];
			WORD c3 = in->data[1 + col + (1 + row) * in->rx];

			switch(pattern) {
			case BAYER_FILTER_RGGB:
				Ha->data[j] = (in->bitpix == 8) ? round_to_BYTE(c0) : c0;
				break;
			case BAYER_FILTER_BGGR:
				Ha->data[j] = (in->bitpix == 8) ? round_to_BYTE(c3) : c3;
				break;
			case BAYER_FILTER_GRBG:
				Ha->data[j] = (in->bitpix == 8) ? round_to_BYTE(c1) : c1;
				break;
			case BAYER_FILTER_GBRG:
				Ha->data[j] = (in->bitpix == 8) ? round_to_BYTE(c2) : c2;
				break;
			default:
				printf("Should not happen.\n");
				return 1;
			}
			j++;
		}
	}

	/* We update FITS keywords */
	copy_fits_metadata(in, Ha);
	update_sampling_information(Ha);

	return 0;
}

int extractHa_float(fits *in, fits *Ha, sensor_pattern pattern) {
	int width = in->rx / 2, height = in->ry / 2;

	if (strlen(in->bayer_pattern) > 4) {
		siril_log_message(_("Extract_Ha does not work on non-Bayer filter camera images!\n"));
		return 1;
	}
	if (new_fit_image(&Ha, width, height, 1, DATA_FLOAT)) {
		return 1;
	}

	int j = 0;

	for (int row = 0; row < in->ry - 1; row += 2) {
		for (int col = 0; col < in->rx - 1; col += 2) {
			float c0 = in->fdata[col + row * in->rx];
			float c1 = in->fdata[1 + col + row * in->rx];
			float c2 = in->fdata[col + (1 + row) * in->rx];
			float c3 = in->fdata[1 + col + (1 + row) * in->rx];

			switch(pattern) {
			case BAYER_FILTER_RGGB:
				Ha->fdata[j] = c0;
				break;
			case BAYER_FILTER_BGGR:
				Ha->fdata[j] = c3;
				break;
			case BAYER_FILTER_GRBG:
				Ha->fdata[j] = c1;
				break;
			case BAYER_FILTER_GBRG:
				Ha->fdata[j] = c2;
				break;
			default:
				printf("Should not happen.\n");
				return 1;
			}
			j++;
		}
	}

	/* We update FITS keywords */
	copy_fits_metadata(in, Ha);
	update_sampling_information(Ha);

	return 0;
}

sensor_pattern get_bayer_pattern(fits *fit) {
	/* Get Bayer informations from header if available */
	sensor_pattern tmp_pattern = com.pref.debayer.bayer_pattern;
	if (com.pref.debayer.use_bayer_header) {
		sensor_pattern bayer;
		bayer = retrieveBayerPatternFromChar(fit->bayer_pattern);

		if (bayer <= BAYER_FILTER_MAX) {
			if (bayer != tmp_pattern) {
				if (bayer == BAYER_FILTER_NONE) {
					siril_log_color_message(_("No Bayer pattern found in the header file.\n"), "salmon");
				}
				else {
					siril_log_color_message(_("Bayer pattern found in header (%s) is different"
								" from Bayer pattern in settings (%s). Overriding settings.\n"),
							"salmon", filter_pattern[bayer], filter_pattern[com.pref.debayer.bayer_pattern]);
					tmp_pattern = bayer;
				}
			}
		} else {
			siril_log_message(_("XTRANS pattern not supported for this feature.\n"));
			return 1;
		}
	}
	if (tmp_pattern >= BAYER_FILTER_MIN && tmp_pattern <= BAYER_FILTER_MAX) {
		siril_log_message(_("Filter Pattern: %s\n"),
				filter_pattern[tmp_pattern]);
	}

	retrieve_Bayer_pattern(fit, &tmp_pattern);
	return tmp_pattern;
}

int extractHa_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_) {
	int ret = 1;
	fits f_Ha = { 0 };
	sensor_pattern pattern = get_bayer_pattern(fit);

	if (fit->type == DATA_USHORT)
		ret = extractHa_ushort(fit, &f_Ha, pattern);
	else if (fit->type == DATA_FLOAT)
		ret = extractHa_float(fit, &f_Ha, pattern);
	else return 1;
	if (!ret) {
		clearfits(fit);
		memcpy(fit, &f_Ha, sizeof(fits));
	}
	return ret;
}

void apply_extractHa_to_sequence(struct split_cfa_data *split_cfa_args) {
	struct generic_seq_args *args = create_default_seqargs(split_cfa_args->seq);
	args->seq = split_cfa_args->seq;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = split_cfa_args->seq->selnum;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = seq_finalize_hook;
	args->image_hook = extractHa_image_hook;
	args->description = _("Extract Ha");
	args->has_output = TRUE;
	args->new_seq_prefix = split_cfa_args->seqEntry;
	args->load_new_sequence = TRUE;
	args->force_ser_output = FALSE;
	args->user = split_cfa_args;

	split_cfa_args->fit = NULL;	// not used here

	start_in_new_thread(generic_sequence_worker, args);
}

int extractGreen_ushort(fits *in, fits *green, sensor_pattern pattern) {
	int width = in->rx / 2, height = in->ry / 2;

	if (strlen(in->bayer_pattern) > 4) {
		siril_log_message(_("Extract_Green does not work on non-Bayer filter camera images!\n"));
		return 1;
	}
	if (new_fit_image(&green, width, height, 1, DATA_USHORT))
		return 1;

	int j = 0;

	for (int row = 0; row < in->ry - 1; row += 2) {
		for (int col = 0; col < in->rx - 1; col += 2) {
			WORD c0, c1, c2, c3;
			switch(pattern) {
			case BAYER_FILTER_RGGB:
			case BAYER_FILTER_BGGR:
				c1 = in->data[1 + col + row * in->rx];
				c2 = in->data[col + (1 + row) * in->rx];
				green->data[j] = (c1 + c2) / 2;
				break;
			case BAYER_FILTER_GRBG:
			case BAYER_FILTER_GBRG:
				c0 = in->data[col + row * in->rx];
				c3 = in->data[1 + col + (1 + row) * in->rx];
				green->data[j] = (c0 + c3) / 2;
				break;
			default:
				printf("Should not happen.\n");
				return 1;
			}
			j++;
		}
	}

	/* We update FITS keywords */
	copy_fits_metadata(in, green);
	update_sampling_information(green);

	return 0;
}

int extractGreen_float(fits *in, fits *green, sensor_pattern pattern) {
	int width = in->rx / 2, height = in->ry / 2;

	if (strlen(in->bayer_pattern) > 4) {
		siril_log_message(_("Extract_Green does not work on non-Bayer filter camera images!\n"));
		return 1;
	}
	if (new_fit_image(&green, width, height, 1, DATA_FLOAT))
		return 1;

	int j = 0;

	for (int row = 0; row < in->ry - 1; row += 2) {
		for (int col = 0; col < in->rx - 1; col += 2) {
			float c0, c1, c2, c3;
			switch(pattern) {
			case BAYER_FILTER_RGGB:
			case BAYER_FILTER_BGGR:
				c1 = in->fdata[1 + col + row * in->rx];
				c2 = in->fdata[col + (1 + row) * in->rx];
				green->fdata[j] = (c1 + c2) * 0.5f;
				break;
			case BAYER_FILTER_GRBG:
			case BAYER_FILTER_GBRG:
				c0 = in->fdata[col + row * in->rx];
				c3 = in->fdata[1 + col + (1 + row) * in->rx];
				green->fdata[j] = (c0 + c3) * 0.5f;
				break;
			default:
				printf("Should not happen.\n");
				return 1;
			}
			j++;
		}
	}

	/* We update FITS keywords */
	copy_fits_metadata(in, green);
	update_sampling_information(green);

	return 0;
}

int extractGreen_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_) {
	int ret = 1;
	fits f_Ha = { 0 };
	sensor_pattern pattern = get_bayer_pattern(fit);

	if (fit->type == DATA_USHORT)
		ret = extractGreen_ushort(fit, &f_Ha, pattern);
	else if (fit->type == DATA_FLOAT)
		ret = extractGreen_float(fit, &f_Ha, pattern);
	else return 1;
	if (!ret) {
		clearfits(fit);
		memcpy(fit, &f_Ha, sizeof(fits));
	}
	return ret;
}

void apply_extractGreen_to_sequence(struct split_cfa_data *split_cfa_args) {
	struct generic_seq_args *args = create_default_seqargs(split_cfa_args->seq);
	args->seq = split_cfa_args->seq;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = split_cfa_args->seq->selnum;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = seq_finalize_hook;
	args->image_hook = extractGreen_image_hook;
	args->description = _("Extract Green");
	args->has_output = TRUE;
	args->new_seq_prefix = split_cfa_args->seqEntry;
	args->load_new_sequence = TRUE;
	args->force_ser_output = FALSE;
	args->user = split_cfa_args;

	split_cfa_args->fit = NULL;	// not used here

	start_in_new_thread(generic_sequence_worker, args);
}

int extractHaOIII_ushort(fits *in, fits *Ha, fits *OIII, sensor_pattern pattern) {
	int width = in->rx / 2, height = in->ry / 2;

	if (strlen(in->bayer_pattern) > 4) {
		siril_log_message(_("Extract_Ha does not work on non-Bayer filter camera images!\n"));
		return 1;
	}
	if (new_fit_image(&Ha, width, height, 1, DATA_USHORT) ||
			new_fit_image(&OIII, width, height, 1, DATA_USHORT)) {
		return 1;
	}

	int j = 0;

	for (int row = 0; row < in->ry - 1; row += 2) {
		for (int col = 0; col < in->rx - 1; col += 2) {
			WORD c0 = in->data[col + row * in->rx];
			WORD c1 = in->data[1 + col + row * in->rx];
			WORD c2 = in->data[col + (1 + row) * in->rx];
			WORD c3 = in->data[1 + col + (1 + row) * in->rx];

			switch(pattern) {
			case BAYER_FILTER_RGGB:
				Ha->data[j] = (in->bitpix == 8) ? round_to_BYTE(c0) : c0;
				OIII->data[j] = (c1 + c2 + c3) / 3;
				break;
			case BAYER_FILTER_BGGR:
				Ha->data[j] = (in->bitpix == 8) ? round_to_BYTE(c3) : c3;
				OIII->data[j] = (c1 + c2 + c0) / 3;
				break;
			case BAYER_FILTER_GRBG:
				Ha->data[j] = (in->bitpix == 8) ? round_to_BYTE(c1) : c1;
				OIII->data[j] = (c0 + c2 + c3) / 3;
				break;
			case BAYER_FILTER_GBRG:
				Ha->data[j] = (in->bitpix == 8) ? round_to_BYTE(c2) : c2;
				OIII->data[j] = (c1 + c0 + c3) / 3;
				break;
			default:
				printf("Should not happen.\n");
				return 1;
			}
			j++;
		}
	}

	/* We update FITS keywords */
	copy_fits_metadata(in, Ha);
	update_sampling_information(Ha);

	copy_fits_metadata(in, OIII);
	update_sampling_information(OIII);

	return 0;
}

int extractHaOIII_float(fits *in, fits *Ha, fits *OIII, sensor_pattern pattern) {
	int width = in->rx / 2, height = in->ry / 2;

	if (strlen(in->bayer_pattern) > 4) {
		siril_log_message(_("Extract_HaOIII does not work on non-Bayer filter camera images!\n"));
		return 1;
	}
	if (new_fit_image(&Ha, width, height, 1, DATA_FLOAT) ||
			new_fit_image(&OIII, width, height, 1, DATA_FLOAT)) {
		return 1;
	}

	int j = 0;

	for (int row = 0; row < in->ry - 1; row += 2) {
		for (int col = 0; col < in->rx - 1; col += 2) {
			float c0 = in->fdata[col + row * in->rx];
			float c1 = in->fdata[1 + col + row * in->rx];
			float c2 = in->fdata[col + (1 + row) * in->rx];
			float c3 = in->fdata[1 + col + (1 + row) * in->rx];

			switch(pattern) {
			case BAYER_FILTER_RGGB:
				Ha->fdata[j] = c0;
				OIII->fdata[j] = (c1 + c2 + c3) / 3;
				break;
			case BAYER_FILTER_BGGR:
				Ha->fdata[j] = c3;
				OIII->fdata[j] = (c1 + c2 + c0) / 3;
				break;
			case BAYER_FILTER_GRBG:
				Ha->fdata[j] = c1;
				OIII->fdata[j] = (c0 + c2 + c3) / 3;
				break;
			case BAYER_FILTER_GBRG:
				Ha->fdata[j] = c2;
				OIII->fdata[j] = (c1 + c0 + c3) / 3;
				break;
			default:
				printf("Should not happen.\n");
				return 1;
			}
			j++;
		}
	}

	/* We update FITS keywords */
	copy_fits_metadata(in, Ha);
	update_sampling_information(Ha);

	copy_fits_metadata(in, OIII);
	update_sampling_information(OIII);

	return 0;
}

struct _double_split {
	int index;
	fits *ha;
	fits *oiii;
};

int extractHaOIII_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_) {
	int ret = 1;
	struct split_cfa_data *cfa_args = (struct split_cfa_data *) args->user;

	sensor_pattern pattern = get_bayer_pattern(fit);

	/* Demosaic and store images for write */
	struct _double_split *double_data = malloc(sizeof(struct _double_split));
	double_data->ha = calloc(1, sizeof(fits));
	double_data->oiii = calloc(1, sizeof(fits));
	double_data->index = o;

	if (fit->type == DATA_USHORT) {
		ret = extractHaOIII_ushort(fit, double_data->ha, double_data->oiii, pattern);
	}
	else if (fit->type == DATA_FLOAT) {
		ret = extractHaOIII_float(fit, double_data->ha, double_data->oiii, pattern);
	}

	if (ret) {
		clearfits(double_data->ha);
		clearfits(double_data->oiii);
		free(double_data);
	} else {
#ifdef _OPENMP
		omp_set_lock(&args->lock);
#endif
		cfa_args->processed_images = g_list_append(cfa_args->processed_images, double_data);
#ifdef _OPENMP
		omp_unset_lock(&args->lock);
#endif
		siril_debug_print("Ha-OIII: processed images added to the save list (%d)\n", o);
	}
	return ret;
}

static int dual_prepare(struct generic_seq_args *args) {
	struct split_cfa_data *cfa_args = (struct split_cfa_data *) args->user;
	// we call the generic prepare twice with different prefixes
	args->new_seq_prefix = "Ha_";
	if (seq_prepare_hook(args))
		return 1;
	// but we copy the result between each call
	cfa_args->new_ser_ha = args->new_ser;
	cfa_args->new_fitseq_ha = args->new_fitseq;

	args->new_seq_prefix = "OIII_";
	if (seq_prepare_hook(args))
		return 1;
	cfa_args->new_ser_oiii = args->new_ser;
	cfa_args->new_fitseq_oiii = args->new_fitseq;

	args->new_seq_prefix = NULL;
	args->new_ser = NULL;
	args->new_fitseq = NULL;

	seqwriter_set_number_of_outputs(2);
	return 0;
}

static int dual_finalize(struct generic_seq_args *args) {
	struct split_cfa_data *cfa_args = (struct split_cfa_data *) args->user;
	args->new_ser = cfa_args->new_ser_ha;
	args->new_fitseq = cfa_args->new_fitseq_ha;
	int retval = seq_finalize_hook(args);
	cfa_args->new_ser_ha = NULL;
	cfa_args->new_fitseq_ha = NULL;

	args->new_ser = cfa_args->new_ser_oiii;
	args->new_fitseq = cfa_args->new_fitseq_oiii;
	retval = seq_finalize_hook(args) || retval;
	cfa_args->new_ser_oiii = NULL;
	cfa_args->new_fitseq_oiii = NULL;
	seqwriter_set_number_of_outputs(1);
	return retval;
}

static int dual_save(struct generic_seq_args *args, int out_index, int in_index, fits *fit) {
	struct split_cfa_data *cfa_args = (struct split_cfa_data *) args->user;
	struct _double_split *double_data = NULL;
	// images are passed from the image_hook to the save in a list, because
	// there are two, which is unsupported by the generic arguments
#ifdef _OPENMP
	omp_set_lock(&args->lock);
#endif
	GList *list = cfa_args->processed_images;
	while (list) {
		if (((struct _double_split *)list->data)->index == out_index) {
			double_data = list->data;
			break;
		}
		list = g_list_next(cfa_args->processed_images);
	}
	if (double_data)
		cfa_args->processed_images = g_list_remove(cfa_args->processed_images, double_data);
#ifdef _OPENMP
	omp_unset_lock(&args->lock);
#endif
	if (!double_data) {
		siril_log_color_message(_("Image %d not found for writing\n"), "red", in_index);
		return 1;
	}

	siril_debug_print("Ha-OIII: images to be saved (%d)\n", out_index);
	if (double_data->ha->naxes[0] == 0 || double_data->oiii->naxes[0] == 0) {
		siril_debug_print("empty data\n");
		return 1;
	}

	int retval1, retval2;
	if (args->force_ser_output || args->seq->type == SEQ_SER) {
		retval1 = ser_write_frame_from_fit(cfa_args->new_ser_ha, double_data->ha, out_index);
		retval2 = ser_write_frame_from_fit(cfa_args->new_ser_oiii, double_data->oiii, out_index);
		clearfits(double_data->ha);
		clearfits(double_data->oiii);
	} else if (args->force_fitseq_output || args->seq->type == SEQ_FITSEQ) {
		retval1 = fitseq_write_image(cfa_args->new_fitseq_ha, double_data->ha, out_index);
		retval2 = fitseq_write_image(cfa_args->new_fitseq_oiii, double_data->oiii, out_index);
		// the two fits are freed by the writing thread
		if (!retval1 && !retval2) {
			/* special case because it's not done in the generic */
			clearfits(fit);
			free(fit);
		}
	} else {
		char *dest = fit_sequence_get_image_filename_prefixed(args->seq, "Ha_", in_index);
		if (fit->type == DATA_USHORT) {
			retval1 = save1fits16(dest, double_data->ha, RLAYER);
		} else {
			retval1 = save1fits32(dest, double_data->ha, RLAYER);
		}
		free(dest);
		dest = fit_sequence_get_image_filename_prefixed(args->seq, "OIII_", in_index);
		if (fit->type == DATA_USHORT) {
			retval2 = save1fits16(dest, double_data->oiii, RLAYER);
		} else {
			retval2 = save1fits32(dest, double_data->oiii, RLAYER);
		}
		free(dest);
		clearfits(double_data->ha);
		clearfits(double_data->oiii);
	}
	free(double_data);
	return retval1 || retval2;
}

void apply_extractHaOIII_to_sequence(struct split_cfa_data *split_cfa_args) {
	struct generic_seq_args *args = create_default_seqargs(split_cfa_args->seq);
	args->seq = split_cfa_args->seq;
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = split_cfa_args->seq->selnum;
	args->prepare_hook = dual_prepare;
	args->finalize_hook = dual_finalize;
	args->save_hook = dual_save;
	args->image_hook = extractHaOIII_image_hook;
	args->description = _("Extract Ha and OIII");
	args->has_output = TRUE;
	args->output_type = get_data_type(args->seq->bitpix);
	args->upscale_ratio = 1.23;	// sqrt(1.5), for memory management
	args->new_seq_prefix = NULL;
	args->user = split_cfa_args;

	split_cfa_args->fit = NULL;	// not used here

	start_in_new_thread(generic_sequence_worker, args);
}

int split_cfa_ushort(fits *in, fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3) {
	int width = in->rx / 2, height = in->ry / 2;

	if (strlen(in->bayer_pattern) > 4) {
		siril_log_message(_("Split CFA does not work on non-Bayer filter camera images!\n"));
		return 1;
	}
	if (new_fit_image(&cfa0, width, height, 1, DATA_USHORT) ||
			new_fit_image(&cfa1, width, height, 1, DATA_USHORT) ||
			new_fit_image(&cfa2, width, height, 1, DATA_USHORT) ||
			new_fit_image(&cfa3, width, height, 1, DATA_USHORT)) {
		return 1;
	}

	int j = 0;

	for (int row = 0; row < in->ry - 1; row += 2) {
		for (int col = 0; col < in->rx - 1; col += 2) {
			/* not c0, c1, c2 and c3 because of the read orientation */
			WORD c1 = in->data[col + row * in->rx];
			WORD c3 = in->data[1 + col + row * in->rx];
			WORD c0 = in->data[col + (1 + row) * in->rx];
			WORD c2 = in->data[1 + col + (1 + row) * in->rx];

			cfa0->data[j] = (in->bitpix == 8) ? round_to_BYTE(c0) : c0;
			cfa1->data[j] = (in->bitpix == 8) ? round_to_BYTE(c1) : c1;
			cfa2->data[j] = (in->bitpix == 8) ? round_to_BYTE(c2) : c2;
			cfa3->data[j] = (in->bitpix == 8) ? round_to_BYTE(c3) : c3;
			j++;
		}
	}

	copy_fits_metadata(in, cfa0);
	copy_fits_metadata(in, cfa1);
	copy_fits_metadata(in, cfa2);
	copy_fits_metadata(in, cfa3);

	/* we remove Bayer header because not needed now */
	clear_Bayer_information(cfa0);
	clear_Bayer_information(cfa1);
	clear_Bayer_information(cfa2);
	clear_Bayer_information(cfa3);

	return 0;
}

int split_cfa_float(fits *in, fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3) {
	int width = in->rx / 2, height = in->ry / 2;

	if (strlen(in->bayer_pattern) > 4) {
		siril_log_message(_("Split CFA does not work on non-Bayer filter camera images!\n"));
		return 1;
	}
	if (new_fit_image(&cfa0, width, height, 1, DATA_FLOAT) ||
			new_fit_image(&cfa1, width, height, 1, DATA_FLOAT) ||
			new_fit_image(&cfa2, width, height, 1, DATA_FLOAT) ||
			new_fit_image(&cfa3, width, height, 1, DATA_FLOAT)) {
		return 1;
	}

	int j = 0;

	for (int row = 0; row < in->ry - 1; row += 2) {
		for (int col = 0; col < in->rx - 1; col += 2) {
			/* not c0, c1, c2 and c3 because of the read orientation */
			float c1 = in->fdata[col + row * in->rx];
			float c3 = in->fdata[1 + col + row * in->rx];
			float c0 = in->fdata[col + (1 + row) * in->rx];
			float c2 = in->fdata[1 + col + (1 + row) * in->rx];

			cfa0->fdata[j] = c0;
			cfa1->fdata[j] = c1;
			cfa2->fdata[j] = c2;
			cfa3->fdata[j] = c3;
			j++;
		}
	}

	copy_fits_metadata(in, cfa0);
	copy_fits_metadata(in, cfa1);
	copy_fits_metadata(in, cfa2);
	copy_fits_metadata(in, cfa3);

	/* we remove Bayer header because not needed now */
	clear_Bayer_information(cfa0);
	clear_Bayer_information(cfa1);
	clear_Bayer_information(cfa2);
	clear_Bayer_information(cfa3);

	return 0;
}

int split_cfa_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_) {
	int ret = 1;
	struct split_cfa_data *cfa_args = (struct split_cfa_data *) args->user;

	fits f_cfa0 = { 0 }, f_cfa1 = { 0 }, f_cfa2 = { 0 }, f_cfa3 = { 0 };

	gchar *cfa0 = g_strdup_printf("%s0_%s%05d%s", cfa_args->seqEntry, cfa_args->seq->seqname, o, com.pref.ext);
	gchar *cfa1 = g_strdup_printf("%s1_%s%05d%s", cfa_args->seqEntry, cfa_args->seq->seqname, o, com.pref.ext);
	gchar *cfa2 = g_strdup_printf("%s2_%s%05d%s", cfa_args->seqEntry, cfa_args->seq->seqname, o, com.pref.ext);
	gchar *cfa3 = g_strdup_printf("%s3_%s%05d%s", cfa_args->seqEntry, cfa_args->seq->seqname, o, com.pref.ext);

	if (fit->type == DATA_USHORT) {
		if (!(ret = split_cfa_ushort(fit, &f_cfa0, &f_cfa1, &f_cfa2, &f_cfa3))) {
			ret = save1fits16(cfa0, &f_cfa0, RLAYER) ||
				save1fits16(cfa1, &f_cfa1, RLAYER) ||
				save1fits16(cfa2, &f_cfa2, RLAYER) ||
				save1fits16(cfa3, &f_cfa3, RLAYER);
		}
	}
	else if (fit->type == DATA_FLOAT) {
		if (!(ret = split_cfa_float(fit, &f_cfa0, &f_cfa1, &f_cfa2, &f_cfa3))) {
			ret = save1fits32(cfa0, &f_cfa0, RLAYER) ||
				save1fits32(cfa1, &f_cfa1, RLAYER) ||
				save1fits32(cfa2, &f_cfa2, RLAYER) ||
				save1fits32(cfa3, &f_cfa3, RLAYER);
		}
	}

	g_free(cfa0); g_free(cfa1);
	g_free(cfa2); g_free(cfa3);
	clearfits(&f_cfa0); clearfits(&f_cfa1);
	clearfits(&f_cfa2); clearfits(&f_cfa3);
	return ret;
}

void apply_split_cfa_to_sequence(struct split_cfa_data *split_cfa_args) {
	struct generic_seq_args *args = create_default_seqargs(split_cfa_args->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = split_cfa_args->seq->selnum;
	args->image_hook = split_cfa_image_hook;
	args->description = _("Split CFA");
	args->new_seq_prefix = split_cfa_args->seqEntry;
	args->user = split_cfa_args;

	split_cfa_args->fit = NULL;	// not used here

	start_in_new_thread(generic_sequence_worker, args);
}

/******* SPLIT CFA ******************************/

void on_split_cfa_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("split_cfa_dialog");
}

void on_split_cfa_apply_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton *seq = GTK_TOGGLE_BUTTON(lookup_widget("checkSplitCFASeq"));
	GtkEntry *entrySplitCFA = GTK_ENTRY(lookup_widget("entrySplitCFA"));
	gint method = gtk_combo_box_get_active(GTK_COMBO_BOX(lookup_widget("combo_split_cfa_method")));

	if (gtk_toggle_button_get_active(seq) && sequence_is_loaded()) {
		struct split_cfa_data *args = calloc(1, sizeof(struct split_cfa_data));

		set_cursor_waiting(TRUE);
		args->seq = &com.seq;
		args->seqEntry = gtk_entry_get_text(entrySplitCFA);
		switch (method) {
			case 0:
				if (args->seqEntry && args->seqEntry[0] == '\0')
					args->seqEntry = "CFA_";
				apply_split_cfa_to_sequence(args);
				break;
			case 1:
				if (args->seqEntry && args->seqEntry[0] == '\0')
					args->seqEntry = "Ha_";
				apply_extractHa_to_sequence(args);
				break;
			case 2:
				apply_extractHaOIII_to_sequence(args);
				break;
			case 3:
				if (args->seqEntry && args->seqEntry[0] == '\0')
					args->seqEntry = "Green_";
				apply_extractGreen_to_sequence(args);
				break;
			default:
				fprintf(stderr, "unhandled case!\n");
		}
	} else {
		switch (method) {
			case 0:
				process_split_cfa(0);
				break;
			case 1:
				process_extractHa(0);
				break;
			case 2:
				process_extractHaOIII(0);
				break;
			case 3:
				process_extractGreen(0);
				break;
			default:
				fprintf(stderr, "unhandled case!\n");
		}
	}
}

void on_combo_split_cfa_method_changed(GtkComboBox *box, gpointer user_data) {
	GtkWidget *w = lookup_widget("label10");
	GtkWidget *txt = lookup_widget("entrySplitCFA");
	gint method = gtk_combo_box_get_active(box);

	gtk_widget_set_sensitive(w, method != 2);
	gtk_widget_set_sensitive(txt, method != 2);
	switch (method) {
		case 0:
			gtk_entry_set_text(GTK_ENTRY(txt), "CFA_");
			break;
		case 1:
			gtk_entry_set_text(GTK_ENTRY(txt), "Ha_");
			break;
		case 3:
			gtk_entry_set_text(GTK_ENTRY(txt), "Green_");
			break;
	}
}

