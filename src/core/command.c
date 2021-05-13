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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_histogram.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <opencv2/core/version.hpp>
#include <glib.h>
#include <libgen.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/arithm.h"
#include "core/undo.h"
#include "core/initfile.h"
#include "core/preprocess.h"
#include "core/processing.h"
#include "core/sequence_filtering.h"
#include "core/OS_utils.h"
#include "io/conversion.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "gui/utils.h"
#include "gui/callbacks.h"
#include "gui/PSF_list.h"
#include "gui/histogram.h"
#include "gui/plot.h"
#include "gui/progress_and_log.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/linear_match.h"
#include "gui/sequence_list.h"
#include "gui/siril_preview.h"
#include "gui/script_menu.h"
#include "filters/asinh.h"
#include "filters/banding.h"
#include "filters/clahe.h"
#include "filters/cosmetic_correction.h"
#include "filters/deconv.h"
#include "filters/median.h"
#include "filters/fft.h"
#include "filters/rgradient.h"
#include "filters/saturation.h"
#include "filters/scnr.h"
#include "filters/wavelets.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/Def_Math.h"
#include "algos/Def_Wavelet.h"
#include "algos/background_extraction.h"
#include "algos/demosaicing.h"
#include "algos/colors.h"
#include "algos/quality.h"
#include "algos/noise.h"
#include "algos/statistics.h"
#include "algos/sorting.h"
#include "algos/siril_wcs.h"
#include "algos/geometry.h"
#include "opencv/opencv.h"
#include "stacking/stacking.h"
#include "stacking/sum.h"
#include "registration/registration.h"
#include "registration/matching/match.h"
#include "algos/fix_xtrans_af.h"
#include "algos/annotate.h"

#include "command.h"
#include "command_def.h"
#include "command_list.h"
#include "command_line_processor.h"

#define PRINT_NOT_FOR_SEQUENCE siril_log_message(_("This command cannot be applied on a sequence.\n"))
#define PRINT_NOT_FOR_SINGLE siril_log_message(_("This command can only be used when a sequence is loaded.\n"))

char *word[MAX_COMMAND_WORDS];	// NULL terminated

int process_load(int nb){
	char filename[256];
	
	strncpy(filename, word[1], 250);
	filename[250] = '\0';
	
	for (int i = 1; i < nb - 1; ++i) {
		strcat(filename, " ");
		strcat(filename, word[i + 1]);
	}
	expand_home_in_filename(filename, 256);
	int retval = open_single_image(filename);
	return (retval < 0);
}

int process_satu(int nb){
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	struct enhance_saturation_data *args = malloc(sizeof(struct enhance_saturation_data));
	
	args->coeff = g_ascii_strtod(word[1], NULL);

	args->input = &gfit;
	args->output = &gfit;
	args->h_min = 0.0;
	args->h_max = 360.0;
	args->background_factor = 1.0;

	enhance_saturation(args);

	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();

	return 0;
}

int process_save(int nb){
	if (sequence_is_loaded() && !single_image_is_loaded()) {
		gfit.hi = com.seq.layers[RLAYER].hi;
		gfit.lo = com.seq.layers[RLAYER].lo;
	}
	else if (single_image_is_loaded()) {
		gfit.hi = com.uniq->layers[RLAYER].hi;
		gfit.lo = com.uniq->layers[RLAYER].lo;
	} else {
		return 1;
	}

	gchar *filename = g_strdup(word[1]);
	set_cursor_waiting(TRUE);
	int retval = savefits(filename, &(gfit));
	set_precision_switch();
	set_cursor_waiting(FALSE);
	g_free(filename);
	return retval;
}

int process_savebmp(int nb){
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	gchar *filename = g_strdup_printf("%s.bmp", word[1]);

	set_cursor_waiting(TRUE);
	savebmp(filename, &(gfit));
	set_cursor_waiting(FALSE);
	g_free(filename);
	return 0;
}

#ifdef HAVE_LIBJPEG
int process_savejpg(int nb){
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	int quality = 100;
	
	if ((nb == 3) && g_ascii_strtoull(word[2], NULL, 10) <= 100 && g_ascii_strtoull(word[2], NULL, 10) > 0)
		quality = g_ascii_strtoull(word[2], NULL, 10);

	gchar *filename = g_strdup_printf("%s.jpg", word[1]);

	set_cursor_waiting(TRUE);
	savejpg(filename, &gfit, quality);
	set_cursor_waiting(FALSE);
	g_free(filename);
	return 0;
}
#endif

#ifdef HAVE_LIBPNG
int process_savepng(int nb){

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	gchar *filename = g_strdup_printf("%s.png", word[1]);

	set_cursor_waiting(TRUE);
	uint32_t bytes_per_sample = gfit.orig_bitpix != BYTE_IMG ? 2 : 1;
	savepng(filename, &gfit, bytes_per_sample, gfit.naxes[2] == 3);
	set_cursor_waiting(FALSE);
	g_free(filename);
	return 0;
}
#endif

#ifdef HAVE_LIBTIFF
int process_savetif(int nb){
	uint16 bitspersample = 16;

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	if (strcasecmp(word[0], "savetif8") == 0)
		bitspersample = 8;
	else if (strcasecmp(word[0], "savetif32") == 0)
		bitspersample = 32;
	gchar *filename = g_strdup_printf("%s.tif", word[1]);
	set_cursor_waiting(TRUE);
	savetif(filename, &gfit, bitspersample);
	set_cursor_waiting(FALSE);
	g_free(filename);
	return 0;
}
#endif

int process_savepnm(int nb){
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	saveNetPBM(word[1], &gfit);
	return 0;
}

int process_imoper(int nb){
	fits fit = { 0 };
	if (!single_image_is_loaded()) return 1;
	if (readfits(word[1], &fit, NULL, !com.pref.force_to_16bit)) return -1;

	image_operator oper;
	switch (word[0][1]) {
		case 'a':
		case 'A':
			oper = OPER_ADD;
			break;
		case 's':
		case 'S':
			oper = OPER_SUB;
			break;
		case 'm':
		case 'M':
			oper = OPER_MUL;
			break;
		case 'd':
		case 'D':
			oper = OPER_DIV;
			break;
		default:
			siril_log_color_message(_("Could not understand the requested operator\n"), "red");
			clearfits(&fit);
			return 1;
	}
	int retval = imoper(&gfit, &fit, oper, !com.pref.force_to_16bit);

	clearfits(&fit);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return retval;
}

int process_addmax(int nb){
	fits fit = { 0 };

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	if (readfits(word[1], &fit, NULL, gfit.type == DATA_FLOAT))
		return -1;
	if (addmax(&gfit, &fit) == 0) {
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	}
	clearfits(&fit);
	return 0;
}

int process_fdiv(int nb){
	// combines an image division and a scalar multiplication.
	fits fit = { 0 };
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	float norm = g_ascii_strtod(word[2], NULL);
	if (readfits(word[1], &fit, NULL, !com.pref.force_to_16bit)) return -1;
	siril_fdiv(&gfit, &fit, norm, TRUE);

	clearfits(&fit);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_fmul(int nb){
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	float coeff = g_ascii_strtod(word[1], NULL);
	if (coeff <= 0.f) {
		siril_log_message(_("Multiplying by a coefficient less than or equal to 0 is not possible.\n"));
		return 1;
	}
	soper(&gfit, coeff, OPER_MUL, TRUE);

	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_entropy(int nb){
	rectangle area;
	float e;

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	if (com.selection.w > 0 && com.selection.h > 0) {
		memcpy(&area, &com.selection, sizeof(rectangle));
		e = entropy(&gfit, com.cvport, &area, NULL);
	}
	else {
		e = entropy(&gfit, com.cvport, NULL, NULL);
	}
	siril_log_message(_("Entropy: %.3f\n"), e);
	return 0;
}

int process_gauss(int nb){
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	unsharp(&(gfit), g_ascii_strtod(word[1], NULL), 0.0, TRUE);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_grey_flat(int nb) {
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	compute_grey_flat(&gfit);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();

	return 0;
}

int process_rl(int nb) {
	double sigma, corner;
	int iter;

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	if (!com.script)
		control_window_switch_to_tab(OUTPUT_LOGS);
	sigma = g_ascii_strtod(word[1], NULL);
	corner = g_ascii_strtod(word[2], NULL);
	iter = g_ascii_strtoull(word[3], NULL, 10);

	if (sigma < 0.4 || sigma > 2.0) {
		siril_log_message(_("Sigma must be between [0.4, 2.0]\n"));
		return 1;
	}

	if (corner < -0.5 || corner > 0.5) {
		siril_log_message(_("Corner radius boost must be between [0.5, 0.5]\n"));
		return 1;
	}

	if (iter <= 0) {
		siril_log_message(_("Number of iterations must be > 0.\n"));
		return 1;
	}

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	struct deconv_data *args = malloc(sizeof(struct deconv_data));

	args->fit = &gfit;
	if (args->fit->type == DATA_USHORT) {
		args->clip = (args->fit->maxi <= 0) ? USHRT_MAX_DOUBLE : args->fit->maxi;
	} else {
		args->clip = (args->fit->maxi <= 0) ? USHRT_MAX_DOUBLE : args->fit->maxi * USHRT_MAX_DOUBLE;
	}
	args->auto_contrast_threshold = TRUE;
	args->sigma = sigma;
	args->corner_radius = corner;
	args->iterations = (size_t)iter;
	args->auto_limit = TRUE;

	set_cursor_waiting(TRUE);

	start_in_new_thread(RTdeconv, args);

	return 0;
}

int process_unsharp(int nb) {
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	unsharp(&(gfit), g_ascii_strtod(word[1], NULL), g_ascii_strtod(word[2], NULL), TRUE);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_crop(int nb) {
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}
	if (is_preview_active()) {
		siril_log_message(_("It is impossible to crop the image when a filter with preview session is active. "
						"Please consider to close the filter dialog first.\n"));
		return 1;
	}


	rectangle area;
	if ((!com.selection.h) || (!com.selection.w)) {
		if (nb == 5) {
			if (g_ascii_strtoull(word[1], NULL, 10) < 0 || g_ascii_strtoull(word[2], NULL, 10) < 0) {
				siril_log_message(_("Crop: x and y must be positive values.\n"));
				return 1;
			}
			if (g_ascii_strtoull(word[3], NULL, 10) <= 0 || g_ascii_strtoull(word[4], NULL, 10) <= 0) {
				siril_log_message(_("Crop: width and height must be greater than 0.\n"));
				return 1;
			}
			if (g_ascii_strtoull(word[1], NULL, 10) + g_ascii_strtoull(word[3], NULL, 10) > gfit.rx || g_ascii_strtoull(word[2], NULL, 10) + g_ascii_strtoull(word[4], NULL, 10) > gfit.ry) {
				siril_log_message(_("Crop: width and height, respectively, must be less than %d and %d.\n"), gfit.rx,gfit.ry);
				return 1;
			}
			area.x = g_ascii_strtoull(word[1], NULL, 10);
			area.y = g_ascii_strtoull(word[2], NULL, 10);
			area.w = g_ascii_strtoull(word[3], NULL, 10);
			area.h = g_ascii_strtoull(word[4], NULL, 10);
		}
		else {
			siril_log_message(_("Crop: select a region or provide x, y, width, height\n"));
			return 1;
		}
	} else {
		memcpy(&area, &com.selection, sizeof(rectangle));
	}

	crop(&gfit, &area);
	delete_selected_area();
	reset_display_offset();
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	
	return 0;
}

int process_cd(int nb) {
	char filename[256];
	int retval;

	g_strlcpy(filename, word[1], 250);

	expand_home_in_filename(filename, 256);
	retval = siril_change_dir(filename, NULL);
	if (!retval) {
		writeinitfile();
		if (!com.script) {
			set_GUI_CWD();
		}
	}
	return retval;
}

int process_wrecons(int nb) {
	int i;
	float coef[7];
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];
	const char *tmpdir;
	int nb_chan;

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	nb_chan = gfit.naxes[2];

	g_assert(nb_chan == 1 || nb_chan == 3);

	tmpdir = g_get_tmp_dir();

	for (i = 0; i < nb - 1; ++i) {
		coef[i] = g_ascii_strtod(word[i + 1], NULL);
	}

	for (i = 0; i < nb_chan; i++) {
		dir[i] = g_build_filename(tmpdir, File_Name_Transform[i], NULL);
		if (gfit.type == DATA_USHORT) {
			wavelet_reconstruct_file(dir[i], coef, gfit.pdata[i]);
		} else if (gfit.type == DATA_FLOAT) {
			wavelet_reconstruct_file_float(dir[i], coef, gfit.fpdata[i]);
		}
		else return 1;
		g_free(dir[i]);
	}

	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_wavelet(int nb) {
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];
	const char* tmpdir;
	int Type_Transform, Nbr_Plan, maxplan, mins, chan, nb_chan;

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	tmpdir = g_get_tmp_dir();

	Nbr_Plan = g_ascii_strtoull(word[1], NULL, 10);
	Type_Transform = g_ascii_strtoull(word[2], NULL, 10);
	
	nb_chan = gfit.naxes[2];
	g_assert(nb_chan <= 3);

	mins = min (gfit.rx, gfit.ry);
	maxplan = log(mins) / log(2) - 2;

	if ( Nbr_Plan > maxplan ){
		siril_log_message(_("Wavelet: maximum number of plans for this image size is %d\n"),
				maxplan);
		return 1;
	}

	if(Type_Transform != TO_PAVE_LINEAR && Type_Transform != TO_PAVE_BSPLINE){
		siril_log_message(_("Wavelet: type must be %d or %d\n"), TO_PAVE_LINEAR, TO_PAVE_BSPLINE);
		return 1;
	}

	if (gfit.type == DATA_USHORT) {
		float *Imag = f_vector_alloc(gfit.rx * gfit.ry);
		if (!Imag) {
			PRINT_ALLOC_ERR;
			return 1;
		}

		for (chan = 0; chan < nb_chan; chan++) {
			dir[chan] = g_build_filename(tmpdir, File_Name_Transform[chan],	NULL);
			wavelet_transform_file(Imag, gfit.ry, gfit.rx, dir[chan],
					Type_Transform, Nbr_Plan, gfit.pdata[chan]);
			g_free(dir[chan]);
		}

		free(Imag);
	} else if (gfit.type == DATA_FLOAT) {
		for (chan = 0; chan < nb_chan; chan++) {
			dir[chan] = g_build_filename(tmpdir, File_Name_Transform[chan],	NULL);
			wavelet_transform_file_float(gfit.fpdata[chan], gfit.ry, gfit.rx, dir[chan],
					Type_Transform, Nbr_Plan);
			g_free(dir[chan]);
		}
	}
	else return 1;
	return 0;
}

int process_log(int nb){
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	loglut(&gfit);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_linear_match(int nb) {
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}
	fits ref = { 0 };
	double a[3] = { 0.0 }, b[3] = { 0.0 };
	double low = g_ascii_strtod(word[2], NULL);
	double high = g_ascii_strtod(word[3], NULL);
	if (readfits(word[1], &ref, NULL, gfit.type == DATA_FLOAT))
		return 1;
	set_cursor_waiting(TRUE);
	if (!find_linear_coeff(&gfit, &ref, low, high, a, b, NULL)) {
		apply_linear_to_fits(&gfit, a, b);

		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	}
	clearfits(&ref);
	set_cursor_waiting(FALSE);
	return 0;
}

int process_asinh(int nb) {
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	double beta = g_ascii_strtod(word[1], NULL);

	asinhlut(&gfit, beta, 0, FALSE);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_clahe(int nb) {
	double clip_limit;
	int size;

	if (CV_MAJOR_VERSION < 3) {
		siril_log_message(_("Your version of opencv is "
				"too old for this feature. Please upgrade your system."));
		return 1;
	}

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	if (!com.script)
		control_window_switch_to_tab(OUTPUT_LOGS);
	clip_limit = g_ascii_strtod(word[1], NULL);

	if (clip_limit <= 0.0) {
		siril_log_message(_("Clip limit must be > 0.\n"));
		return 1;
	}

	size = g_ascii_strtoull(word[2], NULL, 10);

	if (size <= 0.0) {
		siril_log_message(_("Tile size must be > 0.\n"));
		return 1;
	}

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	struct CLAHE_data *args = malloc(sizeof(struct CLAHE_data));

	args->fit = &gfit;
	args->clip = clip_limit;
	args->tileSize = size;

	set_cursor_waiting(TRUE);

	start_in_new_thread(clahe, args);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();

	return 0;
}

#ifndef _WIN32
int process_ls(int nb){
	struct dirent **list;
	gchar *path = NULL;
	
	/* If a path is given in argument */
	if (nb > 1) {
		if (word[1][0] != '\0') {
			/* Absolute path */
			if (word[1][0] == G_DIR_SEPARATOR || word[1][0] == '~') {
				char filename[256];
				
				g_strlcpy(filename, word[1], 250);
				filename[250] = '\0';
				expand_home_in_filename(filename, 256);
				path = g_build_filename(filename, NULL);
			}
			/* Relative path */
			else {
				path = g_build_filename(com.wd, word[1], NULL);
			}
		}
		/* Should not happen */
		else {
			printf("Cannot list files in %s\n", word[1]);
			return 1;
		}
	}
	/* No paths are given in argument */
	else {
		if (!com.wd) {
			siril_log_message(_("Cannot list files, set working directory first.\n"));
			return 1;
		}
		path = g_strdup(com.wd);
	}
	if (path == NULL) {
		siril_log_message(_("Siril cannot open the directory.\n"));
		return 1;
	}

	int i, n;

	n = scandir(path, &list, 0, alphasort);
	if (n < 0) {
		perror("scandir");
		siril_log_message(_("Siril cannot open the directory.\n"));
		g_free(path);
		return 1;
	}

	/* List the entries */
	for (i = 0; i < n; ++i) {
		GStatBuf entrystat;
		gchar *filename;
		const char *ext;
		if (list[i]->d_name[0] == '.')
			continue; /* no hidden files */

		filename = g_build_filename(path, list[i]->d_name, NULL);

		if (g_lstat(filename, &entrystat)) {
			perror("stat");
			g_free(filename);
			break;
		}
		g_free(filename);
		if (S_ISLNK(entrystat.st_mode)) {
			siril_log_color_message(_("Link: %s\n"), "bold", list[i]->d_name);
			continue;
		}
		if (S_ISDIR(entrystat.st_mode)) {
			siril_log_color_message(_("Directory: %s\n"), "green",
					list[i]->d_name);
			continue;
		}
		ext = get_filename_ext(list[i]->d_name);
		if (!ext)
			continue;
		image_type type = get_type_for_extension(ext);
		if (type != TYPEUNDEF) {
			if (type == TYPEAVI || type == TYPESER)
				siril_log_color_message(_("Sequence: %s\n"), "salmon",
						list[i]->d_name);
			else if (type == TYPEFITS)
				siril_log_color_message(_("Image: %s\n"), "plum", list[i]->d_name);
			else
				siril_log_color_message(_("Image: %s\n"), "red", list[i]->d_name);
		} else if (!strncmp(ext, "seq", 3))
			siril_log_color_message(_("Sequence: %s\n"), "blue", list[i]->d_name);
	}
	for (i = 0; i < n; i++)
		free(list[i]);
	free(list);
	siril_log_message(_("********* END OF THE LIST *********\n"));
	g_free(path);

	return 0;
}
#endif

int process_merge(int nb) {
	int retval = 0, nb_seq = nb-2;
	if (!com.wd) {
		siril_log_message(_("Merge: no working directory set.\n"));
		set_cursor_waiting(FALSE);
		return 1;
	}
	char *dest_dir = strdup(com.wd);
	sequence **seqs = calloc(nb_seq, sizeof(sequence *));
	GList *list = NULL;
	for (int i = 0; i < nb_seq; i++) {
		char *seqpath1 = strdup(word[i + 1]), *seqpath2 = strdup(word[i + 1]);
		char *dir = g_path_get_dirname(seqpath1);
		char *seqname = g_path_get_basename(seqpath2);
#ifdef _WIN32
		gchar **token = g_strsplit(dir, "/", -1);
		g_free(dir);
		dir = g_strjoinv(G_DIR_SEPARATOR_S, token);
		g_strfreev(token);
#endif
		if (dir[0] != '\0' && !(dir[0] == '.' && dir[1] == '\0'))
			siril_change_dir(dir, NULL);
		if (!(seqs[i] = load_sequence(seqname, NULL))) {
			siril_log_message(_("Could not open sequence `%s' for merging\n"), word[i + 1]);
			retval = 1;
			free(seqpath1); free(seqpath2);	g_free(seqname); g_free(dir);
			goto merge_clean_up;
		}
		g_free(seqname);
		if (seq_check_basic_data(seqs[i], FALSE) < 0) {
			siril_log_message(_("Sequence `%s' is invalid, could not merge\n"), word[i + 1]);
			retval = 1;
			free(seqpath1); free(seqpath2); g_free(dir);
			goto merge_clean_up;
		}

		if (i != 0 && (seqs[i]->rx != seqs[0]->rx ||
					seqs[i]->ry != seqs[0]->ry ||
					seqs[i]->nb_layers != seqs[0]->nb_layers ||
					seqs[i]->bitpix != seqs[0]->bitpix ||
					seqs[i]->type != seqs[0]->type)) {
			siril_log_message(_("All sequences must be the same format for merging. Sequence `%s' is different\n"), word[i + 1]);
			retval = 1;
			free(seqpath1); free(seqpath2); g_free(dir);
			goto merge_clean_up;
		}

		if (seqs[i]->type == SEQ_REGULAR) {
			// we need to build the list of files
			char filename[256];
			for (int image = 0; image < seqs[i]->number; image++) {
				fit_sequence_get_image_filename(seqs[i], image, filename, TRUE);
				list = g_list_append(list, g_build_filename(dir, filename, NULL));
			}
		}
		free(seqpath1); free(seqpath2); g_free(dir);
		siril_change_dir(dest_dir, NULL);	// they're all relative to this one
	}

	char *outseq_name;
	struct ser_struct out_ser;
	struct _convert_data *args;
	fitseq out_fitseq;
	switch (seqs[0]->type) {
		case SEQ_REGULAR:
			// use the conversion, it makes symbolic links or copies as a fallback
			args = malloc(sizeof(struct _convert_data));
			args->start = 0;
			args->total = 0; // init to get it from glist_to_array()
			args->list = glist_to_array(list, &args->total);
			args->destroot = format_basename(word[nb - 1], FALSE);
			args->input_has_a_seq = FALSE;
			args->input_has_a_film = FALSE;
			args->debayer = FALSE;
			args->multiple_output = FALSE;
			args->output_type = SEQ_REGULAR;
			args->make_link = TRUE;
			gettimeofday(&(args->t_start), NULL);
			start_in_new_thread(convert_thread_worker, args);
			break;

		case SEQ_SER:
			if (g_str_has_suffix(word[nb - 1], ".ser"))
				outseq_name = g_strdup(word[nb - 1]);
			else outseq_name = g_strdup_printf("%s.ser", word[nb - 1]);
			if (ser_create_file(outseq_name, &out_ser, TRUE, seqs[0]->ser_file)) {
				siril_log_message(_("Failed to create the output SER file `%s'\n"), word[nb - 1]);
				retval = 1;
				goto merge_clean_up;
			}
			free(outseq_name);
			seqwriter_set_max_active_blocks(2);
			int written_frames = 0;
			for (int i = 0; i < nb_seq; i++) {
				for (unsigned int frame = 0; frame < seqs[i]->number; frame++) {
					seqwriter_wait_for_memory();
					fits *fit = calloc(1, sizeof(fits));
					if (ser_read_frame(seqs[i]->ser_file, frame, fit, FALSE)) {
						siril_log_message(_("Failed to read frame %d from input sequence `%s'\n"), frame, word[i + 1]);
						retval = 1;
						seqwriter_release_memory();
						ser_close_and_delete_file(&out_ser);
						goto merge_clean_up;
					}

					if (ser_write_frame_from_fit(&out_ser, fit, written_frames)) {
						siril_log_message(_("Failed to write frame %d in merged sequence\n"), written_frames);
						retval = 1;
						seqwriter_release_memory();
						ser_close_and_delete_file(&out_ser);
						goto merge_clean_up;
					}
					written_frames++;
				}
			}
			if (ser_write_and_close(&out_ser)) {
				siril_log_message(_("Error while finalizing the merged sequence\n"));
				retval = 1;
			}
			break;

		case SEQ_FITSEQ:
			if (g_str_has_suffix(word[nb - 1], com.pref.ext))
				outseq_name = g_strdup(word[nb - 1]);
			else outseq_name = g_strdup_printf("%s%s", word[nb - 1], com.pref.ext);
			if (fitseq_create_file(outseq_name, &out_fitseq, -1)) {
				siril_log_message(_("Failed to create the output SER file `%s'\n"), word[nb - 1]);
				retval = 1;
				goto merge_clean_up;
			}
			free(outseq_name);
			seqwriter_set_max_active_blocks(2);
			written_frames = 0;
			for (int i = 0; i < nb_seq; i++) {
				for (unsigned int frame = 0; frame < seqs[i]->number; frame++) {
					seqwriter_wait_for_memory();
					fits *fit = calloc(1, sizeof(fits));
					if (fitseq_read_frame(seqs[i]->fitseq_file, frame, fit, FALSE, -1)) {
						siril_log_message(_("Failed to read frame %d from input sequence `%s'\n"), frame, word[i + 1]);
						retval = 1;
						seqwriter_release_memory();
						fitseq_close_and_delete_file(&out_fitseq);
						goto merge_clean_up;
					}

					if (fitseq_write_image(&out_fitseq, fit, written_frames)) {
						siril_log_message(_("Failed to write frame %d in merged sequence\n"), written_frames);
						retval = 1;
						seqwriter_release_memory();
						fitseq_close_and_delete_file(&out_fitseq);
						goto merge_clean_up;
					}
					written_frames++;
				}
			}
			if (fitseq_close_file(&out_fitseq)) {
				siril_log_message(_("Error while finalizing the merged sequence\n"));
				retval = 1;
			}
			break;
		default:
			siril_log_message(_("This type of sequence cannot be created by Siril, aborting the merge\n"));
			retval = 1;
	}

merge_clean_up:
	for (int i = 0; i < nb_seq; i++) {
		if (seqs[i])
			free_sequence(seqs[i], TRUE);
	}
	free(seqs);
	siril_change_dir(dest_dir, NULL);
	free(dest_dir);
	return retval;
}

int	process_mirrorx(int nb){
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	mirrorx(&gfit, TRUE);
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int	process_mirrory(int nb){
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	mirrory(&gfit, TRUE);
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_mtf(int nb) {
	float lo = g_ascii_strtod(word[1], NULL);
	float mid = g_ascii_strtod(word[2], NULL);
	float hi = g_ascii_strtod(word[3], NULL);

	mtf_with_parameters(&gfit, lo, mid, hi);

	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_resample(int nb) {
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	double factor = g_ascii_strtod(word[1], NULL);
	if (factor > 5.0) {
		siril_log_message(_("The scaling factor must be less than 5.0\n"));
		return 1;
	}
	int toX = round_to_int(factor * gfit.rx);
	int toY = round_to_int(factor * gfit.ry);
	
	set_cursor_waiting(TRUE);
	verbose_resize_gaussian(&gfit, toX, toY, OPENCV_LINEAR);
	
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	return 0;
}

int process_rgradient(int nb) {
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	struct rgradient_filter_data *args = malloc(sizeof(struct rgradient_filter_data));
	args->xc = g_ascii_strtod(word[1], NULL);
	args->yc = g_ascii_strtod(word[2], NULL);
	args->dR = g_ascii_strtod(word[3], NULL);
	args->da = g_ascii_strtod(word[4], NULL);
	args->fit = &gfit;

	if ((args->xc >= args->fit->rx) || (args->yc >= args->fit->ry)) {
		siril_log_message(_("The coordinates cannot be greater than the size of the image. "
				"Please change their values and retry.\n"));
	} else {

		set_cursor_waiting(TRUE);

		start_in_new_thread(rgradient_filter, args);
	}
	return 0;
}

int process_rotate(int nb) {
	double degree;
	
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	set_cursor_waiting(TRUE);
	degree = g_ascii_strtod(word[1], NULL);
	verbose_rotate_image(&gfit, degree, OPENCV_LINEAR, 1);	//INTER_LINEAR
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	return 0;
}

int process_rotatepi(int nb){
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	verbose_rotate_image(&gfit, 180.0, OPENCV_LINEAR, 1);

	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_set_mag(int nb) {
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	int layer = match_drawing_area_widget(com.vport[com.cvport], FALSE);
	double mag = g_ascii_strtod(word[1], NULL);

	if (layer != -1) {

		if (com.selection.w > 300 || com.selection.h > 300){
			siril_log_message(_("Current selection is too large. To determine the PSF, please make a selection around a single star.\n"));
			return 1;
		}
		if (com.selection.w <= 0 || com.selection.h <= 0){
			siril_log_message(_("Select an area first\n"));
			return 1;
		}
		fitted_PSF *result = psf_get_minimisation(&gfit, layer, &com.selection, TRUE, TRUE, TRUE);
		if (result) {
			com.magOffset = mag - result->mag;
			siril_log_message(_("Relative magnitude: %.3lf, "
					"True reduced magnitude: %.3lf, "
					"Offset: %.3lf\n"), result->mag, mag, com.magOffset);
			free(result);
		}
	}
	return 0;
}

int process_set_ref(int nb) {
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq)
		return 1;

	int n = g_ascii_strtoull(word[2], NULL, 10) - 1;
	if (n < 0 || n > seq->number) {
		siril_log_message(_("The reference image must be set between 1 and %d\n"), seq->number);
		return 1;
	}

	seq->reference_image = n;
	// a reference image should not be excluded to avoid confusion
	if (!seq->imgparam[seq->current].incl) {
		seq->imgparam[seq->current].incl = TRUE;
	}

	writeseqfile(seq);

	return 0;
}

int process_unset_mag(int nb) {
	com.magOffset = 0.0;
	return 0;
}

int process_set_mag_seq(int nb) {
	if (!sequence_is_loaded()) {
		PRINT_NOT_FOR_SINGLE;
		return 1;
	}
	double mag = g_ascii_strtod(word[1], NULL);
	int i;
	for (i = 0; i < MAX_SEQPSF && com.seq.photometry[i]; i++);
	com.seq.reference_star = i - 1;
	if (i == 0) {
		siril_log_message(_("Run a PSF for the sequence first (see seqpsf)\n"));
		return 1;
	}
	com.seq.reference_mag = mag;
	siril_log_message(_("Reference magnitude has been set for star %d to %f and will be computed for each image\n"), i - 1, mag);
	drawPlot();
	return 0;
}

int process_set_ext(int nb) {
	if (word[1]) {
		GString *str = NULL;

		if ((g_ascii_strncasecmp(word[1], "fit", 3))
				&& (g_ascii_strncasecmp(word[1], "fts", 3))
				&& (g_ascii_strncasecmp(word[1], "fits", 4))) {
			siril_log_message(_("FITS extension unknown: %s\n"), word[1]);
		}

		free(com.pref.ext);
		str = g_string_new(".");
		str = g_string_append(str, word[1]);
		str = g_string_ascii_down(str);
		com.pref.ext = g_string_free(str, FALSE);
		writeinitfile();
	}

	return 0;
}

int process_set_findstar(int nb) {
	double sigma = g_ascii_strtod(word[1], NULL);
	double roundness = g_ascii_strtod(word[2], NULL);
	int retval = 0;

	if (sigma >= 0.05 && roundness >= 0 && roundness <= 0.9) {
		com.starfinder_conf.sigma = sigma;
		com.starfinder_conf.roundness = roundness;
	} else {
		siril_log_message(_("Wrong parameter values. Sigma must be >= 0.05 and roundness between 0 and 0.9.\n"));
		retval = 1;
	}
	return retval;
}

int process_unset_mag_seq(int nb) {
	if (!sequence_is_loaded()) {
		PRINT_NOT_FOR_SINGLE;
		return 1;
	}
	com.seq.reference_star = -1;
	com.seq.reference_mag = -1001.0;
	siril_log_message(_("Reference magnitude unset for sequence\n"));
	drawPlot();
	return 0;
}

int process_psf(int nb){
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	int layer = match_drawing_area_widget(com.vport[com.cvport], FALSE);
	if (layer != -1) {

		if (com.selection.w > 300 || com.selection.h > 300){
			siril_log_message(_("Current selection is too large. To determine the PSF, please make a selection around a single star.\n"));
			return 1;
		}
		if (com.selection.w <= 0 || com.selection.h <= 0){
			siril_log_message(_("Select an area first\n"));
			return 1;
		}
		fitted_PSF *result = psf_get_minimisation(&gfit, layer, &com.selection, TRUE, TRUE, TRUE);
		if (result) {
			psf_display_result(result, &com.selection);
			free(result);
		}
	}
	return 0;
}

int process_seq_psf(int nb) {
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}
	if (com.selection.w > 300 || com.selection.h > 300){
		siril_log_message(_("Current selection is too large. To determine the PSF, please make a selection around a single star.\n"));
		return 1;
	}
	if (com.selection.w <= 0 || com.selection.h <= 0){
		siril_log_message(_("Select an area first\n"));
		return 1;
	}

	int layer = match_drawing_area_widget(com.vport[com.cvport], FALSE);
	if (sequence_is_loaded() && layer != -1) {
		framing_mode framing = REGISTERED_FRAME;
		if (framing == REGISTERED_FRAME && !com.seq.regparam[layer])
			framing = ORIGINAL_FRAME;
		if (framing == ORIGINAL_FRAME) {
			GtkToggleButton *follow = GTK_TOGGLE_BUTTON(lookup_widget("followStarCheckButton"));
			if (gtk_toggle_button_get_active(follow))
				framing = FOLLOW_STAR_FRAME;
		}
		siril_log_message(_("Running the PSF on the loaded sequence, layer %d\n"), layer);
		seqpsf(&com.seq, layer, FALSE, FALSE, framing, TRUE);
		return 0;
	}
	else {
		PRINT_NOT_FOR_SINGLE;
		return 1;
	}
}

int process_seq_crop(int nb) {
	if (!sequence_is_loaded()) {
		PRINT_NOT_FOR_SINGLE;
		return 1;
	}
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	rectangle area;

	int startoptargs = 5;

	if ((!com.selection.h) || (!com.selection.w)) {
		if (nb >= startoptargs) {
			if (g_ascii_strtoull(word[1], NULL, 10) < 0 || g_ascii_strtoull(word[2], NULL, 10) < 0) {
				siril_log_message(_("Crop: x and y must be positive values.\n"));
				return 1;
			}
			if (g_ascii_strtoull(word[3], NULL, 10) <= 0 || g_ascii_strtoull(word[4], NULL, 10) <= 0) {
				siril_log_message(_("Crop: width and height must be greater than 0.\n"));
				return 1;
			}
			area.x = g_ascii_strtoull(word[1], NULL, 10);
			area.y = g_ascii_strtoull(word[2], NULL, 10);
			area.w = g_ascii_strtoull(word[3], NULL, 10);
			area.h = g_ascii_strtoull(word[4], NULL, 10);
		} else {
			siril_log_message(_("Crop: select a region or provide x, y, width, height\n"));
			return 1;
		}
	} else {
		memcpy(&area, &com.selection, sizeof(rectangle));
	}


	if (g_ascii_strtoull(word[3], NULL, 10) > com.seq.rx || g_ascii_strtoull(word[4], NULL, 10) > com.seq.ry) {
		siril_log_message(_("Crop: width and height, respectively, must be less than %d and %d.\n"),
				com.seq.rx, com.seq.ry);
		return 1;
	}

	struct crop_sequence_data *args = malloc(sizeof(struct crop_sequence_data));

	args->seq = &com.seq;
	args->area = area;
	args->prefix = "cropped_";
	
	if (nb > startoptargs) {
		for (int i = startoptargs; i < nb; i++) {
			if (word[i]) {
				if (g_str_has_prefix(word[i], "-prefix=")) {
					char *current = word[i], *value;
					value = current + 8;
					if (value[0] == '\0') {
						siril_log_message(_("Missing argument to %s, aborting.\n"),	current);
						return 1;
					}
					args->prefix = strdup(value);
				}
			}
		}
	}

	set_cursor_waiting(TRUE);

	crop_sequence(args);
	return 0;
}

int process_bg(int nb){
	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;
	WORD us_bg;

	for (int layer = 0; layer < gfit.naxes[2]; layer++) {
		double bg = background(&gfit, layer, &com.selection, TRUE);
		if (gfit.type == DATA_USHORT) {
			us_bg = round_to_WORD(bg);
			bg = bg / get_normalized_value(&gfit);
		} else if (gfit.type == DATA_FLOAT) {
			us_bg = float_to_ushort_range(bg);
		} else return 1;
		siril_log_message(_("Background value (channel: #%d): %d (%lf)\n"), layer, us_bg, bg);
	}
	return 0;
}

int process_bgnoise(int nb){
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

	struct noise_data *args = malloc(sizeof(struct noise_data));

	if (!com.script) {
		control_window_switch_to_tab(OUTPUT_LOGS);
	}

	args->fit = &gfit;
	args->verbose = TRUE;
	args->use_idle = TRUE;
	memset(args->bgnoise, 0.0, sizeof(double[3]));
	set_cursor_waiting(TRUE);

	start_in_new_thread(noise, args);
	return 0;
}

int process_histo(int nb){
	GError *error = NULL;
	int nlayer = g_ascii_strtoull(word[1], NULL, 10);
	const gchar* clayer;
	
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	if (nlayer>3 || nlayer <0)
		return 1;
	gsl_histogram *histo = computeHisto(&gfit, nlayer);
	if (!isrgb(&gfit))
		clayer = "bw";		//if B&W
	else
		clayer = vport_number_to_name(nlayer);
	gchar *filename = g_strdup_printf("histo_%s.dat", clayer);

	GFile *file = g_file_new_for_path(filename);
	g_free(filename);

	GOutputStream *output_stream = (GOutputStream*) g_file_replace(file, NULL, FALSE,
			G_FILE_CREATE_NONE, NULL, &error);

	if (output_stream == NULL) {
		if (error != NULL) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
			fprintf(stderr, "Cannot save histo\n");
		}
		g_object_unref(file);
		return 1;
	}
	for (size_t i = 0; i < USHRT_MAX + 1; i++) {
		gchar *buffer = g_strdup_printf("%zu %d\n", i, (int) gsl_histogram_get (histo, i));

		if (!g_output_stream_write_all(output_stream, buffer, strlen(buffer), NULL, NULL, &error)) {
			g_warning("%s\n", error->message);
			g_free(buffer);
			g_clear_error(&error);
			g_object_unref(output_stream);
			g_object_unref(file);
			return 1;
		}
		g_free(buffer);
	}

	siril_log_message(_("The file %s has been created for the %s layer.\n"), g_file_peek_path(file), clayer);

	g_object_unref(output_stream);
	g_object_unref(file);
	gsl_histogram_free(histo);
	return 0;
}

int process_thresh(int nb){
	int lo, hi;

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	lo = g_ascii_strtoull(word[1], NULL, 10);
	if (lo < 0 || lo > USHRT_MAX) {
		siril_log_message(_("replacement value is out of range (0 - %d)\n"), USHRT_MAX);
		return 1;
	}
	hi = g_ascii_strtoull(word[2], NULL, 10);
	if (hi < 0 || hi > USHRT_MAX) {
		siril_log_message(_("replacement value is out of range (0 - %d)\n"), USHRT_MAX);
		return 1;
	}
	threshlo(&gfit, lo);
	threshhi(&gfit, hi);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_threshlo(int nb){
	int lo;

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	lo = g_ascii_strtoull(word[1], NULL, 10);
	if (lo < 0 || lo > USHRT_MAX) {
		siril_log_message(_("replacement value is out of range (0 - %d)\n"), USHRT_MAX);
		return 1;
	}
	threshlo(&gfit, lo);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_threshhi(int nb){
	int hi;

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	hi = g_ascii_strtoull(word[1], NULL, 10);
	if (hi < 0 || hi > USHRT_MAX) {
		siril_log_message(_("replacement value is out of range (0 - %d)\n"), USHRT_MAX);
		return 1;
	}
	threshhi(&gfit, hi);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_neg(int nb) {
	set_cursor_waiting(TRUE);
	pos_to_neg(&gfit);
	update_gfit_histogram_if_needed();
	invalidate_stats_from_fit(&gfit);
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	return 0;
}

int process_nozero(int nb){
	int level;

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	level = g_ascii_strtoull(word[1], NULL, 10);
	if (level < 0 || level > USHRT_MAX) {
		siril_log_message(_("replacement value is out of range (0 - %d)\n"), USHRT_MAX);
		return 1;
	}
	nozero(&gfit, (WORD)level);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_ddp(int nb){
	// combines an image division and a scalar multiplication.
	float coeff, sigma;
	unsigned level;

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	level = g_ascii_strtoull(word[1], NULL, 10);
	coeff = g_ascii_strtod(word[2], NULL);
	sigma = g_ascii_strtod(word[3], NULL);
	ddp(&gfit, level, coeff, sigma);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_new(int nb){
	int width, height, layers;
	
	width = g_ascii_strtod(word[1], NULL);
	height = g_ascii_strtod(word[2], NULL);
	layers = g_ascii_strtoull(word[3], NULL, 10);
	if (layers != 1 && layers != 3) {
		siril_log_message(_("Number of layers MUST be 1 or 3\n"));
		return 1;
	}
	if (!height || !width) return 1;

	close_single_image();
	close_sequence(FALSE);

	fits *fit = &gfit;
	if (new_fit_image(&fit, width, height, layers, DATA_FLOAT))
		return 1;
	memset(gfit.fdata, 0, width * height * layers * sizeof(float));

	com.seq.current = UNRELATED_IMAGE;
	com.uniq = calloc(1, sizeof(single));
	com.uniq->filename = strdup(_("new empty image"));
	com.uniq->fileexist = FALSE;
	com.uniq->nb_layers = gfit.naxes[2];
	com.uniq->layers = calloc(com.uniq->nb_layers, sizeof(layer_info));
	com.uniq->fit = &gfit;

	open_single_image_from_gfit();
	return 0;
}

int process_visu(int nb){
	int low, high;
	
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	low = g_ascii_strtoull(word[1], NULL, 10);
	high = g_ascii_strtoull(word[2], NULL, 10);
	if ((high > USHRT_MAX) || (low < 0)) {
		siril_log_message(_("Values must be positive and less than %d.\n"), USHRT_MAX);
		return 1;
	}
	visu(&gfit, low, high);
	return 0;
}

int process_fill2(int nb){
	int level = g_ascii_strtoull(word[1], NULL, 10);
	rectangle area;

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	if ((!com.selection.h) || (!com.selection.w)) {
		if (nb == 6) {
			area.x = g_ascii_strtoull(word[2], NULL, 10);
			area.y = g_ascii_strtoull(word[3], NULL, 10);
			area.w = g_ascii_strtoull(word[4], NULL, 10);
			area.h = g_ascii_strtoull(word[5], NULL, 10);
			if ((area.w + area.x > gfit.rx) || (area.h + area.y > gfit.ry)) {
				siril_log_message(_("Wrong parameters.\n"));
				return 1;
			}
		}
		else {
			siril_log_message(_("Fill2: select a region or provide x, y, width, height\n"));
			return 1;
		}
	} else {
		memcpy(&area, &com.selection, sizeof(rectangle));
	}
	fill(&gfit, level, &area);
	area.x = gfit.rx - area.x - area.w;
	area.y = gfit.ry - area.y - area.h;
	fill(&gfit, level, &area);
	redraw(com.cvport, REMAP_ALL);
	return 0;
}

int process_findstar(int nb){
	int nbstars = 0;

	int layer = com.cvport == RGB_VPORT ? GLAYER : com.cvport;

	delete_selected_area();
	com.stars = peaker(&gfit, layer, &com.starfinder_conf, &nbstars, NULL, TRUE);
	siril_log_message(_("Found %d stars in image, channel #%d\n"), nbstars, layer);
	if (com.stars)
		refresh_star_list(com.stars);
	return 0;
}

int process_findhot(int nb){
	GError *error = NULL;
	long icold, ihot;
	gchar type;

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	if (gfit.naxes[2] != 1) {
		siril_log_message(_("find_hot must be applied on an one-channel master-dark frame"));
		return 1;
	}
	double sig[2];
	sig[0] = g_ascii_strtod(word[2], NULL);
	sig[1] = g_ascii_strtod(word[3], NULL);

	deviant_pixel *dev = find_deviant_pixels(&gfit, sig, &icold, &ihot, FALSE);
	siril_log_message(_("%ld cold and %ld hot pixels\n"), icold, ihot);

	gchar *filename = g_strdup_printf("%s.lst", word[1]);
	GFile *file = g_file_new_for_path(filename);
	g_free(filename);

	GOutputStream *output_stream = (GOutputStream*) g_file_replace(file, NULL, FALSE,
			G_FILE_CREATE_NONE, NULL, &error);

	if (output_stream == NULL) {
		if (error != NULL) {
			g_warning("%s\n", error->message);
			g_clear_error(&error);
			fprintf(stderr, "Cannot open file: %s\n", filename);
		}
		g_object_unref(file);
		return 1;
	}

	for (int i = 0; i < icold + ihot; i++) {
		int y = gfit.ry - (int) dev[i].p.y - 1;  /* FITS is stored bottom to top */
		if (dev[i].type == HOT_PIXEL)
			type = 'H';
		else
			type = 'C';
		gchar *buffer = g_strdup_printf("P %d %d %c\n", (int) dev[i].p.x, y, type);
		if (!g_output_stream_write_all(output_stream, buffer, strlen(buffer), NULL, NULL, &error)) {
			g_warning("%s\n", error->message);
			g_free(buffer);
			g_clear_error(&error);
			g_object_unref(output_stream);
			g_object_unref(file);
			return 1;
		}
		g_free(buffer);
	}

	free(dev);
	g_object_unref(output_stream);
	g_object_unref(file);

	return 0;
}

int process_fix_xtrans(int nb) {
	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

	fix_xtrans_ac(&gfit);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	return 0;
}

int process_cosme(int nb) {
	GError *error = NULL;
	deviant_pixel dev;
	gchar *filename;
	double dirty;
	int is_cfa, i = 0, retval = 0;
	int nb_tokens;
	gchar *line;
	char type;

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	if (!g_str_has_suffix(word[1], ".lst")) {
		filename = g_strdup_printf("%s.lst", word[1]);
	} else {
		filename = g_strdup(word[1]);
	}
	GFile *file = g_file_new_for_path(filename);
	g_free(filename);

	GInputStream *input_stream = (GInputStream *)g_file_read(file, NULL, &error);

	if (input_stream == NULL) {
		if (error != NULL) {
			g_clear_error(&error);
			siril_log_message(_("File [%s] does not exist\n"), filename);
		}

		g_object_unref(file);
		return 1;
	}

	is_cfa = (word[0][5] == '_') ? 1 : 0;

	GDataInputStream *data_input = g_data_input_stream_new(input_stream);
	while ((line = g_data_input_stream_read_line_utf8(data_input, NULL,
				NULL, NULL))) {
		++i;
		switch (line[0]) {
		case '#': // comments.
			g_free(line);
			continue;
			break;
		case 'P':
			nb_tokens = sscanf(line + 2, "%lf %lf %c", &dev.p.x, &dev.p.y, &type);
			if (nb_tokens != 2 && nb_tokens != 3) {
				fprintf(stderr, "cosmetic correction: "
						"cosme file format error at line %d: %s", i, line);
				retval = 1;
				g_free(line);
				continue;
			}
			if (nb_tokens == 2)	{
				type = 'H';
			}
			if (type == 'H')
				dev.type = HOT_PIXEL;
			else
				dev.type = COLD_PIXEL;
			dev.p.y = gfit.ry - dev.p.y - 1;  /* FITS are stored bottom to top */
			cosmeticCorrOnePoint(&gfit, dev, is_cfa);
			break;
		case 'L':
			nb_tokens = sscanf(line + 2, "%lf %lf %c", &dev.p.y, &dirty, &type);
			if (nb_tokens != 2 && nb_tokens != 3) {
				fprintf(stderr, "cosmetic correction: "
						"cosme file format error at line %d: %s\n", i, line);
				retval = 1;
				g_free(line);
				continue;
			}
			dev.type = HOT_PIXEL; // we force it
			dev.p.y = gfit.ry - dev.p.y - 1; /* FITS are stored bottom to top */
			cosmeticCorrOneLine(&gfit, dev, is_cfa);
			break;
		case 'C':
			nb_tokens = sscanf(line + 2, "%lf %lf %c", &dev.p.y, &dirty, &type);
			if (nb_tokens != 2 && nb_tokens != 3) {
				fprintf(stderr, "cosmetic correction: "
						"cosme file format error at line %d: %s\n", i, line);
				retval = 1;
				g_free(line);
				continue;
			}
			point center = {gfit.rx / 2.0, gfit.ry / 2.0};
			dev.type = HOT_PIXEL; // we force it
			dev.p.y = gfit.rx - dev.p.y - 1; /* FITS are stored bottom to top */
			cvRotateImage(&gfit, center, 90.0, -1, OPENCV_LINEAR);
			cosmeticCorrOneLine(&gfit, dev, is_cfa);
			cvRotateImage(&gfit, center, -90.0, -1, OPENCV_LINEAR);

			break;
		default:
			fprintf(stderr, _("cosmetic correction: "
					"cosme file format error at line %d: %s\n"), i, line);
			retval = 1;
		}
		g_free(line);
	}

	g_object_unref(input_stream);
	g_object_unref(file);
	if (retval)
		siril_log_color_message(_("There were some errors, please check your input file.\n"), "salmon");

	invalidate_stats_from_fit(&gfit);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_fmedian(int nb){
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}
	
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	struct median_filter_data *args = malloc(sizeof(struct median_filter_data));
	args->ksize = g_ascii_strtoull(word[1], NULL, 10);
	args->amount = g_ascii_strtod(word[2], NULL);
	args->iterations = 1;
	
	if (!(args->ksize & 1) || args->ksize < 2 || args->ksize > 15) {
		siril_log_message(_("The size of the kernel MUST be odd and in the range [3, 15].\n"));
		free(args);
		return 1;
	}
	if (args->amount < 0.0 || args->amount > 1.0) {
		siril_log_message(_("Modulation value MUST be between 0 and 1\n"));
		free(args);
		return 1;
	}
	args->fit = &gfit;

	set_cursor_waiting(TRUE);

	start_in_new_thread(median_filter, args);
	
	return 0;
}

/* The name of this command should be COG in english but this choice
 * was done to be consistent with IRIS
 */
int process_cdg(int nb) {
	float x_avg, y_avg;

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	if (!FindCentre(&gfit, &x_avg, &y_avg)) {
		siril_log_message(_("Center of gravity coordinates are (%.3lf, %.3lf)\n"), x_avg, y_avg);
		return 0;
	}
	return 1;
}

int process_clear(int nb) {
	if (com.script) return 0;
	GtkTextView *text = GTK_TEXT_VIEW(lookup_widget("output"));
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(text);
	GtkTextIter start_iter, end_iter;
	gtk_text_buffer_get_start_iter(tbuf, &start_iter);
	gtk_text_buffer_get_end_iter(tbuf, &end_iter);
	gtk_text_buffer_delete(tbuf, &start_iter, &end_iter);
	return 0;
}

int process_clearstar(int nb){
	clear_stars_list();
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_NONE);
	redraw_previews();
	return 0;
}

int process_close(int nb) {
	close_sequence(FALSE);
	close_single_image();
	if (!com.script) {
		update_MenuItem();
		reset_plot(); // reset all plots
		close_tab();	//close Green and Blue Tab if a 1-layer sequence is loaded
		
	}
	return 0;
}

int process_fill(int nb){
	int level;
	rectangle area;
	
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	if ((!com.selection.h) || (!com.selection.w)) {
		if (nb == 6) {
			area.x = g_ascii_strtoull(word[2], NULL, 10);
			area.y = g_ascii_strtoull(word[3], NULL, 10);
			area.w = g_ascii_strtoull(word[4], NULL, 10);
			area.h = g_ascii_strtoull(word[5], NULL, 10);
			if ((area.w + area.x > gfit.rx) || (area.h + area.y > gfit.ry)) {
				siril_log_message(_("Wrong parameters.\n"));
				return 1;
			}
		}
		else {
			area.w = gfit.rx; area.h = gfit.ry;
			area.x = 0; area.y = 0;
		}
	} else {
		memcpy(&area, &com.selection, sizeof(rectangle));
	}
	level = g_ascii_strtoull(word[1], NULL, 10);
	fill(&gfit, level, &area);
	redraw(com.cvport, REMAP_ALL);
	return 0;
}

int process_offset(int nb){
	int level;
	
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	level = g_ascii_strtod(word[1], NULL);
	off(&gfit, level);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

/* The version in command line is a minimal version
 * Only neutral type are available (no amount needed), 
 * then we always preserve the lightness */
int process_scnr(int nb){
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}
	if (gfit.naxes[2] == 1) return 1;

	struct scnr_data *args = malloc(sizeof(struct scnr_data));
	
	args->type = g_ascii_strtoull(word[1], NULL, 10);
	args->fit = &gfit;
	args->amount = 0.0;
	args->preserve = TRUE;

	set_cursor_waiting(TRUE);

	start_in_new_thread(scnr, args);

	return 0;
}

int process_fft(int nb){
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	if (sequence_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	struct fft_data *args = malloc(sizeof(struct fft_data));
	
	args->fit = &gfit;
	args->type = strdup(word[0]);
	args->modulus = strdup(word[1]);
	args->phase = strdup(word[2]);
	args->type_order = 0;
	
	set_cursor_waiting(TRUE);

	start_in_new_thread(fourier_transform, args);
	
	return 0;
}

int process_fixbanding(int nb) {
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	struct banding_data *args = malloc(sizeof(struct banding_data));

	args->amount = g_ascii_strtod(word[1], NULL);
	args->sigma = g_ascii_strtod(word[2], NULL);
	args->protect_highlights = TRUE;
	args->fit = &gfit;

	set_cursor_waiting(TRUE);

	start_in_new_thread(BandingEngineThreaded, args);
	
	return 0;
}


int process_subsky(int nb) {
	gboolean is_sequence;
	sequence *seq = NULL;
	int degree = 0;

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	is_sequence = (word[0][2] == 'q');

	if (is_sequence) {
		seq = load_sequence(word[1], NULL);
		if (!seq)
			return 1;
		degree = g_ascii_strtoull(word[2], NULL, 10);
	} else {
		if (!single_image_is_loaded()) return 1;
		degree = g_ascii_strtoull(word[1], NULL, 10);
	}

	if (degree < 1 || degree > 4) {
		siril_log_message(_("Polynomial degree order must be within the [1, 4] range.\n"));
		return 1;
	}

	set_cursor_waiting(TRUE);

	if (is_sequence) {
		struct background_data *args = malloc(sizeof(struct background_data));

		args->seq = seq;
		args->nb_of_samples = 20;
		args->tolerance = 1.0;
		args->correction = 0; //subtraction
		args->seqEntry = "bkg_";
		args->degree = (poly_order) (degree - 1);
		args->dither = TRUE;

		int startoptargs = 3;
		if (nb > startoptargs) {
			for (int i = startoptargs; i < nb; i++) {
				if (word[i]) {
					if (g_str_has_prefix(word[i], "-prefix=")) {
						char *current = word[i], *value;
						value = current + 8;
						if (value[0] == '\0') {
							siril_log_message(_("Missing argument to %s, aborting.\n"),	current);
							return 1;
						}
						args->seqEntry = strdup(value);
					}
				}
			}
		}


		apply_background_extraction_to_sequence(args);
	} else {
		generate_background_samples(20, 1.0);
		remove_gradient_from_image(0, (poly_order) (degree - 1));
		free_background_sample_list(com.grad_samples);
		com.grad_samples = NULL;

		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
	}
	set_cursor_waiting(FALSE);

	return 0;
}


int process_findcosme(int nb) {
	gboolean is_sequence;
	sequence *seq = NULL;
	int i = 0;

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	is_sequence = (word[0][0] == 's');

	if (is_sequence) {
		seq = load_sequence(word[1], NULL);
		if (!seq)
			return 1;
		i++;
	} else {
		if (!single_image_is_loaded()) return 1;
	}

	struct cosmetic_data *args = malloc(sizeof(struct cosmetic_data));

	args->seq = seq;
	args->sigma[0] = g_ascii_strtod(word[1 + i], NULL);
	args->sigma[1] = g_ascii_strtod(word[2 + i], NULL);
	args->is_cfa = (word[0][10] == '_' || word[0][13] == '_');	// find_cosme_cfa or seqfind_cosme_cfa
	args->amount = 1.0;
	args->fit = &gfit;

	set_cursor_waiting(TRUE);

	if (is_sequence) {
		args->seqEntry = "cc_";
		args->multithread = FALSE;

		int startoptargs = i + 3;
		int nb_command_max = i + 4;
		if (nb > startoptargs) {
			for (int j = startoptargs; j < nb_command_max; j++) {
				if (word[j]) {
					if (g_str_has_prefix(word[j], "-prefix=")) {
						char *current = word[j], *value;
						value = current + 8;
						if (value[0] == '\0') {
							siril_log_message(_("Missing argument to %s, aborting.\n"), current);
							return 1;
						}
						args->seqEntry = strdup(value);
					}
				}
			}
		}
		apply_cosmetic_to_sequence(args);
	} else {
		args->multithread = TRUE;
		start_in_new_thread(autoDetectThreaded, args);
	}

	return 0;
}

int select_unselect(gboolean select) {
	if (!sequence_is_loaded()) {
		siril_log_message(_("Use this command to select images in a sequence, load a sequence first.\n"));
		return 1;
	}
	int from = g_ascii_strtoull(word[1], NULL, 10);
	int to = g_ascii_strtoull(word[2], NULL, 10);
	if (from < 0 || from >= com.seq.number) {
		siril_log_message(_("The first argument must be between 0 and the number of images minus one.\n"));
		return 1;
	}
	gboolean current_updated = FALSE;
	for (int i = from; i <= to; i++) {
		if (i >= com.seq.number) break;
		if (com.seq.imgparam[i].incl != select) {
			com.seq.imgparam[i].incl = select;
			if (!com.headless)
				sequence_list_change_selection_index(i, i);
			if (select)
				com.seq.selnum++;
			else	com.seq.selnum--;
			if (i == com.seq.current)
				current_updated = TRUE;
		}
		if (!select && com.seq.reference_image == i) {
			com.seq.reference_image = -1;
			if (!com.headless) {
				sequence_list_change_reference();
				adjust_refimage(com.seq.current);
			}
		}
	}

	if (!com.headless) {
		if (current_updated) {
			redraw(com.cvport, REMAP_NONE);
			drawPlot();
			adjust_sellabel();
		}
		update_reg_interface(FALSE);
		adjust_sellabel();
	}
	writeseqfile(&com.seq);
	siril_log_message(_("Selection update finished, %d images are selected in the sequence\n"), com.seq.selnum);

	return 0;
}

int process_select(int nb){
	return select_unselect(TRUE);
}

int process_unselect(int nb){
	return select_unselect(FALSE);
}

int process_split(int nb){
	if (!single_image_is_loaded()) {
		PRINT_NOT_FOR_SEQUENCE;
		return 1;
	}

	if (!isrgb(&gfit)) {
		siril_log_message(_("Siril cannot split layers. Make sure your image is in RGB mode.\n"));
		return 1;
	}

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	struct extract_channels_data *args = malloc(sizeof(struct extract_channels_data));
	if (!args) {
		PRINT_ALLOC_ERR;
		return 1;
	}

	args->type = 0;
	args->str_type = _("RGB");

	args->channel[0] = g_strdup_printf("%s%s", word[1], com.pref.ext);
	args->channel[1] = g_strdup_printf("%s%s", word[2], com.pref.ext);
	args->channel[2] = g_strdup_printf("%s%s", word[3], com.pref.ext);

	args->fit = calloc(1, sizeof(fits));
	set_cursor_waiting(TRUE);
	if (copyfits(&gfit, args->fit, CP_ALLOC | CP_COPYA | CP_FORMAT, -1)) {
		siril_log_message(_("Could not copy the input image, aborting.\n"));
		free(args->fit);
		free(args->channel[0]);
		free(args->channel[1]);
		free(args->channel[2]);
		free(args);
		return 1;
	}
	copy_fits_metadata(&gfit, args->fit);
	start_in_new_thread(extract_channels, args);
	return 0;
}

int process_split_cfa(int nb) {
	if (isrgb(&gfit)) {
		siril_log_message(_("Siril cannot split CFA channel. Make sure your image is in CFA mode.\n"));
		return 1;
	}
	char *filename = NULL;
	int ret = 1;

	fits f_cfa0 = { 0 }, f_cfa1 = { 0 }, f_cfa2 = { 0 }, f_cfa3 = { 0 };

	if (sequence_is_loaded() && !single_image_is_loaded()) {
		filename = g_path_get_basename(com.seq.seqname);
	}
	else {
		if (com.uniq->filename != NULL) {
			char *tmp = remove_ext_from_filename(com.uniq->filename);
			filename = g_path_get_basename(tmp);
			free(tmp);
		}
	}

	gchar *cfa0 = g_strdup_printf("CFA0_%s%s", filename, com.pref.ext);
	gchar *cfa1 = g_strdup_printf("CFA1_%s%s", filename, com.pref.ext);
	gchar *cfa2 = g_strdup_printf("CFA2_%s%s", filename, com.pref.ext);
	gchar *cfa3 = g_strdup_printf("CFA3_%s%s", filename, com.pref.ext);

	if (gfit.type == DATA_USHORT) {
		if (!(ret = split_cfa_ushort(&gfit, &f_cfa0, &f_cfa1, &f_cfa2, &f_cfa3))) {
			ret = save1fits16(cfa0, &f_cfa0, 0) ||
				save1fits16(cfa1, &f_cfa1, 0) ||
				save1fits16(cfa2, &f_cfa2, 0) ||
				save1fits16(cfa3, &f_cfa3, 0);
		}
	}
	else if (gfit.type == DATA_FLOAT) {
		if (!(ret = split_cfa_float(&gfit, &f_cfa0, &f_cfa1, &f_cfa2, &f_cfa3))) {
			ret = save1fits32(cfa0, &f_cfa0, 0) ||
				save1fits32(cfa1, &f_cfa1, 0) ||
				save1fits32(cfa2, &f_cfa2, 0) ||
				save1fits32(cfa3, &f_cfa3, 0);
		}
	}

	g_free(cfa0); g_free(cfa1);
	g_free(cfa2); g_free(cfa3);
	clearfits(&f_cfa0); clearfits(&f_cfa1);
	clearfits(&f_cfa2); clearfits(&f_cfa3);
	free(filename);
	return ret;
}

int process_extractHa(int nb) {
	if (isrgb(&gfit)) {
		siril_log_message(_("Siril cannot split CFA channel. Make sure your image is in CFA mode.\n"));
		return 1;
	}
	char *filename = NULL;
	int ret = 1;

	fits f_Ha = { 0 };

	if (sequence_is_loaded() && !single_image_is_loaded()) {
		filename = g_path_get_basename(com.seq.seqname);
	}
	else {
		if (com.uniq->filename != NULL) {
			char *tmp = remove_ext_from_filename(com.uniq->filename);
			filename = g_path_get_basename(tmp);
			free(tmp);
		}
	}

	/* Get Bayer informations from header if available */
	sensor_pattern tmp_pattern = com.pref.debayer.bayer_pattern;
	if (com.pref.debayer.use_bayer_header) {
		sensor_pattern bayer;
		bayer = retrieveBayerPatternFromChar(gfit.bayer_pattern);

		if (bayer <= BAYER_FILTER_MAX) {
			if (bayer != tmp_pattern) {
				if (bayer == BAYER_FILTER_NONE) {
					siril_log_color_message(_("No Bayer pattern found in the header file.\n"), "red");
				}
				else {
					siril_log_color_message(_("Bayer pattern found in header (%s) is different"
								" from Bayer pattern in settings (%s). Overriding settings.\n"),
							"salmon", filter_pattern[bayer], filter_pattern[com.pref.debayer.bayer_pattern]);
					tmp_pattern = bayer;
				}
			}
		} else {
			siril_log_message(_("XTRANS pattern not handled for this feature.\n"));
			return 1;
		}
	}
	if (tmp_pattern >= BAYER_FILTER_MIN && tmp_pattern <= BAYER_FILTER_MAX) {
		siril_log_message(_("Filter Pattern: %s\n"),
				filter_pattern[tmp_pattern]);
	}

	retrieve_Bayer_pattern(&gfit, &tmp_pattern);

	gchar *Ha = g_strdup_printf("Ha_%s%s", filename, com.pref.ext);
	if (gfit.type == DATA_USHORT) {
		if (!(ret = extractHa_ushort(&gfit, &f_Ha, tmp_pattern))) {
			ret = save1fits16(Ha, &f_Ha, 0);
		}
	}
	else if (gfit.type == DATA_FLOAT) {
		if (!(ret = extractHa_float(&gfit, &f_Ha, tmp_pattern))) {
			ret = save1fits32(Ha, &f_Ha, 0);
		}
	} else return 1;

	g_free(Ha);
	clearfits(&f_Ha);
	free(filename);
	return ret;
}

int process_extractHaOIII(int nb) {
	if (isrgb(&gfit)) {
		siril_log_message(_("Siril cannot split CFA channel. Make sure your image is in CFA mode.\n"));
		return 1;
	}
	char *filename = NULL;
	int ret = 1;

	fits f_Ha = { 0 }, f_OIII = { 0 };

	if (sequence_is_loaded() && !single_image_is_loaded()) {
		filename = g_path_get_basename(com.seq.seqname);
	}
	else {
		if (com.uniq->filename != NULL) {
			char *tmp = remove_ext_from_filename(com.uniq->filename);
			filename = g_path_get_basename(tmp);
			free(tmp);
		}
	}

	/* Get Bayer informations from header if available */
	sensor_pattern tmp_pattern = com.pref.debayer.bayer_pattern;
	if (com.pref.debayer.use_bayer_header) {
		sensor_pattern bayer;
		bayer = retrieveBayerPatternFromChar(gfit.bayer_pattern);

		if (bayer <= BAYER_FILTER_MAX) {
			if (bayer != tmp_pattern) {
				if (bayer == BAYER_FILTER_NONE) {
					siril_log_color_message(_("No Bayer pattern found in the header file.\n"), "red");
				}
				else {
					siril_log_color_message(_("Bayer pattern found in header (%s) is different"
								" from Bayer pattern in settings (%s). Overriding settings.\n"),
							"salmon", filter_pattern[bayer], filter_pattern[com.pref.debayer.bayer_pattern]);
					tmp_pattern = bayer;
				}
			}
		} else {
			siril_log_message(_("XTRANS pattern not handled for this feature.\n"));
			return 1;
		}
	}
	if (tmp_pattern >= BAYER_FILTER_MIN && tmp_pattern <= BAYER_FILTER_MAX) {
		siril_log_message(_("Filter Pattern: %s\n"),
				filter_pattern[tmp_pattern]);
	}

	retrieve_Bayer_pattern(&gfit, &tmp_pattern);

	gchar *Ha = g_strdup_printf("Ha_%s%s", filename, com.pref.ext);
	gchar *OIII = g_strdup_printf("OIII_%s%s", filename, com.pref.ext);
	if (gfit.type == DATA_USHORT) {
		if (!(ret = extractHaOIII_ushort(&gfit, &f_Ha, &f_OIII, tmp_pattern))) {
			ret = save1fits16(Ha, &f_Ha, 0) ||
					save1fits16(OIII, &f_OIII, 0);
		}
	}
	else if (gfit.type == DATA_FLOAT) {
		if (!(ret = extractHaOIII_float(&gfit, &f_Ha, &f_OIII, tmp_pattern))) {
			ret = save1fits32(Ha, &f_Ha, 0) ||
					save1fits16(OIII, &f_OIII, 0);
		}
	} else return 1;

	g_free(Ha);
	g_free(OIII);
	clearfits(&f_Ha);
	clearfits(&f_OIII);
	free(filename);
	return ret;
}

int process_seq_mtf(int nb) {
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq)
		return 1;

	struct mtf_data *args = malloc(sizeof(struct mtf_data));

	args->seq = seq;
	args->fit = &gfit;
	args->seqEntry = "mtf_";
	args->lo = g_ascii_strtod(word[2], NULL);
	args->mid = g_ascii_strtod(word[3], NULL);
	args->hi = g_ascii_strtod(word[4], NULL);

	int startoptargs = 5;
	if (nb > startoptargs) {
		for (int i = startoptargs; i < nb; i++) {
			if (word[i]) {
				if (g_str_has_prefix(word[i], "-prefix=")) {
					char *current = word[i], *value;
					value = current + 8;
					if (value[0] == '\0') {
						siril_log_message(_("Missing argument to %s, aborting.\n"),	current);
						return 1;
					}
					args->seqEntry = strdup(value);
				}
			}
		}
	}

	set_cursor_waiting(TRUE);

	apply_mtf_to_sequence(args);

	return 0;
}

int process_seq_split_cfa(int nb) {
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	sequence *seq = load_sequence(word[1], NULL);
	if (!seq)
		return 1;

	if (seq->nb_layers > 1) {
		siril_log_message(_("Siril cannot split CFA channel. Make sure your image is in CFA mode.\n"));
		return 1;
	}

	struct split_cfa_data *args = calloc(1, sizeof(struct split_cfa_data));

	args->seq = seq;
	args->fit = &gfit;
	args->seqEntry = "CFA_"; // propose to default to "CFA" for consistency of output names with single image split_cfa

	int startoptargs = 2;
	if (nb > startoptargs) {
		for (int i = startoptargs; i < nb; i++) {
			if (word[i]) {
				if (g_str_has_prefix(word[i], "-prefix=")) {
					char *current = word[i], *value;
					value = current + 8;
					if (value[0] == '\0') {
						siril_log_message(_("Missing argument to %s, aborting.\n"),	current);
						return 1;
					}
					args->seqEntry = strdup(value);
				}
			}
		}
	}

	set_cursor_waiting(TRUE);

	apply_split_cfa_to_sequence(args);

	return 0;
}

int process_seq_extractHa(int nb) {
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	sequence *seq = load_sequence(word[1], NULL);
	if (!seq)
		return 1;

	if (seq->nb_layers > 1) {
		siril_log_message(_("Siril cannot split CFA channel. Make sure your image is in CFA mode.\n"));
		return 1;
	}

	struct split_cfa_data *args = calloc(1, sizeof(struct split_cfa_data));

	args->seq = seq;
	args->seqEntry = "Ha_";

	int startoptargs = 2;
	if (nb > startoptargs) {
		for (int i = startoptargs; i < nb; i++) {
			if (word[i]) {
				if (g_str_has_prefix(word[i], "-prefix=")) {
					char *current = word[i], *value;
					value = current + 8;
					if (value[0] == '\0') {
						siril_log_message(_("Missing argument to %s, aborting.\n"),	current);
						return 1;
					}
					args->seqEntry = strdup(value);
				}
			}
		}
	}

	set_cursor_waiting(TRUE);

	apply_extractHa_to_sequence(args);

	return 0;
}

int process_seq_extractHaOIII(int nb) {
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	sequence *seq = load_sequence(word[1], NULL);
	if (!seq)
		return 1;

	if (seq->nb_layers > 1) {
		siril_log_message(_("Siril cannot split CFA channel. Make sure your image is in CFA mode.\n"));
		return 1;
	}

	struct split_cfa_data *args = calloc(1, sizeof(struct split_cfa_data));

	args->seq = seq;
	args->seqEntry = ""; // not used

	set_cursor_waiting(TRUE);

	apply_extractHaOIII_to_sequence(args);

	return 0;
}

int process_stat(int nb){
	int nplane;
	int layer;
	char layername[6];

	nplane = gfit.naxes[2];

	for (layer = 0; layer < nplane; layer++) {
		imstats* stat = statistics(NULL, -1, &gfit, layer, &com.selection, STATS_MAIN, TRUE);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}

		switch (layer) {
			case 0:
				if (nplane == 1)
					strcpy(layername, "B&W");
				else
					strcpy(layername, "Red");
				break;
			case 1:
				strcpy(layername, "Green");
				break;
			case 2:
				strcpy(layername, "Blue");
				break;
		}

		if (gfit.type == DATA_USHORT) {
			siril_log_message(
					_("%s layer: Mean: %0.1lf, Median: %0.1lf, Sigma: %0.1lf, "
							"AvgDev: %0.1lf, Min: %0.1lf, Max: %0.1lf\n"),
					layername, stat->mean, stat->median, stat->sigma,
					stat->avgDev, stat->min, stat->max);
		} else {
			siril_log_message(
					_("%s layer: Mean: %0.1lf, Median: %0.1lf, Sigma: %0.1lf, "
							"AvgDev: %0.1lf, Min: %0.1lf, Max: %0.1lf\n"),
					layername, stat->mean * USHRT_MAX_DOUBLE,
					stat->median * USHRT_MAX_DOUBLE,
					stat->sigma * USHRT_MAX_DOUBLE,
					stat->avgDev * USHRT_MAX_DOUBLE,
					stat->min * USHRT_MAX_DOUBLE, stat->max * USHRT_MAX_DOUBLE);
		}
		free_stats(stat);
	}
	return 0;
}

int process_seq_stat(int nb) {
	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	sequence *seq = load_sequence(word[1], NULL);
	if (!seq)
		return 1;

	struct stat_data *args = calloc(1, sizeof(struct stat_data));

	args->seq = seq;
	args->seqEntry = ""; // not used
	args->csv_name = g_strdup(word[2]);

	if (word[3] && !g_strcmp0(word[3], "main")) {
		args->option = STATS_MAIN;
	} else {
		args->option = STATS_BASIC;
	}

	set_cursor_waiting(TRUE);

	apply_stats_to_sequence(args);

	return 0;
}

int process_convertraw(int nb) {
	GDir *dir;
	GError *error = NULL;
	const gchar *file;
	GList *list = NULL;
	int idx = 1;
	gchar *destroot = g_strdup(word[1]);
	sequence_type output = SEQ_REGULAR;
	gboolean debayer = FALSE;

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	if (!com.wd) {
		siril_log_message(_("Conversion: no working directory set.\n"));
		return 1;
	}

	for (int i = 2; i < nb; i++) {
		char *current = word[i], *value;
		if (!strcmp(current, "-debayer")) {
			debayer = TRUE;
		} else if (!strcmp(current, "-fitseq")) {
			output = SEQ_FITSEQ;
			if (!g_str_has_suffix(destroot, com.pref.ext))
				str_append(&destroot, com.pref.ext);
		} else if (!strcmp(current, "-ser")) {
			output = SEQ_SER;
			if (!g_str_has_suffix(destroot, ".ser"))
				str_append(&destroot, ".ser");
		} else if (g_str_has_prefix(current, "-start=")) {
			value = current + 7;
			idx = (g_ascii_strtoull(value, NULL, 10) <= 0 || g_ascii_strtoull(value, NULL, 10) >= INDEX_MAX) ? 1 : g_ascii_strtoull(value, NULL, 10);
		} else if (g_str_has_prefix(current, "-out=")) {
			value = current + 5;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				return 1;
			}
			if (!g_file_test(value, G_FILE_TEST_EXISTS)) {
				if (g_mkdir_with_parents(value, 0755) < 0) {
					siril_log_color_message(_("Cannot create output folder: %s\n"), "red", value);
					return 1;
				}
			}
			gchar *filename = g_build_filename(value, destroot, NULL);
			g_free(destroot);
			destroot = filename;
		}
	}

	if ((dir = g_dir_open(com.wd, 0, &error)) == NULL){
		siril_log_message(_("Conversion: error opening working directory %s.\n"), com.wd);
		fprintf (stderr, "Conversion: %s\n", error->message);
		g_clear_error(&error);
		return 1;
	}

	int count = 0;
	while ((file = g_dir_read_name(dir)) != NULL) {
		const char *ext = get_filename_ext(file);
		if (!ext)
			continue;
		image_type type = get_type_for_extension(ext);
		if (type == TYPERAW) {
			if (output == SEQ_SER && !g_ascii_strcasecmp(ext, "raf") && !debayer) {
				siril_log_message(_("FujiFilm XTRANS sensors are not supported by SER v2 (CFA-style) standard. You may use FITS sequences instead."));
				g_list_free_full(list, g_free);
				return 1;
			}
			list = g_list_append(list, g_build_filename(com.wd, file, NULL));
			count++;
		}
	}
	g_dir_close(dir);
	if (!count) {
		siril_log_message(_("No RAW files were found for conversion\n"));
		g_list_free_full(list, g_free);
		return 1;
	}
	/* sort list */
	list = g_list_sort(list, (GCompareFunc) strcompare);
	/* convert the list to an array for parallel processing */
	char **files_to_convert = glist_to_array(list, &count);

	siril_log_color_message(_("Conversion: processing %d RAW files...\n"), "green", count);

	set_cursor_waiting(TRUE);
	if (!com.script)
		control_window_switch_to_tab(OUTPUT_LOGS);

	struct _convert_data *args = malloc(sizeof(struct _convert_data));
	args->start = idx;
	args->list = files_to_convert;
	args->total = count;
	if (output == SEQ_REGULAR)
		args->destroot = format_basename(destroot, TRUE);
	else
		args->destroot = destroot;
	args->input_has_a_seq = FALSE;
	args->input_has_a_film = FALSE;
	args->debayer = debayer;
	args->output_type = output;
	args->multiple_output = FALSE;
	args->make_link = FALSE;
	gettimeofday(&(args->t_start), NULL);
	start_in_new_thread(convert_thread_worker, args);
	return 0;
}

int process_link(int nb) {
	GDir *dir;
	GError *error = NULL;
	const gchar *file;
	GList *list = NULL;
	int idx = 1;
	gchar *destroot = g_strdup(word[1]);

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	for (int i = 2; i < nb; i++) {
		char *current = word[i], *value;
		if (g_str_has_prefix(current, "-start=")) {
			value = current + 7;
			idx = (g_ascii_strtoull(value, NULL, 10) <= 0 ||
					g_ascii_strtoull(value, NULL, 10) >= INDEX_MAX) ?
				1 : g_ascii_strtoull(value, NULL, 10);
		} else if (g_str_has_prefix(current, "-out=")) {
			value = current + 5;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				return 1;
			}
			if (!g_file_test(value, G_FILE_TEST_EXISTS)) {
				if (g_mkdir_with_parents(value, 0755) < 0) {
					siril_log_color_message(_("Cannot create output folder: %s\n"), "red", value);
					return 1;
				}
			}
			gchar *filename = g_build_filename(value, destroot, NULL);
			g_free(destroot);
			destroot = filename;
		}
	}

	if ((dir = g_dir_open(com.wd, 0, &error)) == NULL){
		siril_log_message(_("Link: error opening working directory %s.\n"), com.wd);
		fprintf (stderr, "Link: %s\n", error->message);
		g_clear_error(&error);
		set_cursor_waiting(FALSE);
		return 1;
	}

	int count = 0;
	while ((file = g_dir_read_name(dir)) != NULL) {
		const char *ext = get_filename_ext(file);
		if (!ext)
			continue;
		image_type type = get_type_for_extension(ext);
		if (type == TYPEFITS) {
			list = g_list_append(list, g_build_filename(com.wd, file, NULL));
			count++;
		}
	}
	g_dir_close(dir);
	if (!count) {
		siril_log_message(_("No FITS files were found for link\n"));
		return 1;
	}
	/* sort list */
	list = g_list_sort(list, (GCompareFunc) strcompare);
	/* convert the list to an array for parallel processing */
	char **files_to_link = glist_to_array(list, &count);

	gchar *str = ngettext("Link: processing %d FITS file...\n", "Link: processing %d FITS files...\n", count);
	str = g_strdup_printf(str, count);
	siril_log_color_message(str, "green");
	g_free(str);

	set_cursor_waiting(TRUE);
	if (!com.script)
		control_window_switch_to_tab(OUTPUT_LOGS);

	if (!com.wd) {
		siril_log_message(_("Link: no working directory set.\n"));
		set_cursor_waiting(FALSE);
		return 1;
	}

	struct _convert_data *args = malloc(sizeof(struct _convert_data));
	args->start = idx;
	args->list = files_to_link;
	args->total = count;
	args->destroot = format_basename(destroot, TRUE);
	args->input_has_a_seq = FALSE;
	args->input_has_a_film = FALSE;
	args->debayer = FALSE;
	args->multiple_output = FALSE;
	args->output_type = SEQ_REGULAR; // fallback if symlink does not work
	args->make_link = TRUE;
	gettimeofday(&(args->t_start), NULL);
	start_in_new_thread(convert_thread_worker, args);

	return 0;
}

int process_convert(int nb) {
	GDir *dir;
	GError *error = NULL;
	const gchar *file;
	GList *list = NULL;
	int idx = 1;
	gboolean debayer = FALSE;
	gboolean make_link = TRUE;
	sequence_type output = SEQ_REGULAR;
	gchar *destroot = g_strdup(word[1]);

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}

	for (int i = 2; i < nb; i++) {
		char *current = word[i], *value;
		if (!strcmp(current, "-debayer")) {
			debayer = TRUE;
			make_link = FALSE;
		} else if (!strcmp(current, "-fitseq")) {
			output = SEQ_FITSEQ;
			if (!g_str_has_suffix(destroot, com.pref.ext))
				str_append(&destroot, com.pref.ext);
		} else if (!strcmp(current, "-ser")) {
			output = SEQ_SER;
			if (!g_str_has_suffix(destroot, ".ser"))
				str_append(&destroot, ".ser");
		} else if (g_str_has_prefix(current, "-start=")) {
			value = current + 7;
			idx = (g_ascii_strtoull(value, NULL, 10) <= 0 || g_ascii_strtoull(value, NULL, 10) >= INDEX_MAX) ?
				1 : g_ascii_strtoull(value, NULL, 10);
		} else if (g_str_has_prefix(current, "-out=")) {
			value = current + 5;
			if (value[0] == '\0') {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				return 1;
			}
			if (!g_file_test(value, G_FILE_TEST_EXISTS)) {
				if (g_mkdir_with_parents(value, 0755) < 0) {
					siril_log_color_message(_("Cannot create output folder: %s\n"), "red", value);
					return 1;
				}
			}
			gchar *filename = g_build_filename(value, destroot, NULL);
			g_free(destroot);
			destroot = filename;
		}
	}

	if ((dir = g_dir_open(com.wd, 0, &error)) == NULL){
		siril_log_message(_("Convert: error opening working directory %s.\n"), com.wd);
		fprintf (stderr, "Convert: %s\n", error->message);
		g_clear_error(&error);
		set_cursor_waiting(FALSE);
		return 1;
	}

	int count = 0;
	while ((file = g_dir_read_name(dir)) != NULL) {
		const char *ext = get_filename_ext(file);
		if (!ext)
			continue;
		image_type type = get_type_for_extension(ext);
		if (type != TYPEUNDEF && type != TYPEAVI && type != TYPESER) {
			list = g_list_append(list, g_build_filename(com.wd, file, NULL));
			count++;
		}
	}
	g_dir_close(dir);
	if (!count) {
		siril_log_message(_("No files were found for convert\n"));
		return 1;
	}
	/* sort list */
	list = g_list_sort(list, (GCompareFunc) strcompare);
	/* convert the list to an array for parallel processing */
	char **files_to_link = glist_to_array(list, &count);

	gchar *str = ngettext("Convert: processing %d FITS file...\n", "Convert: processing %d FITS files...\n", count);
	str = g_strdup_printf(str, count);
	siril_log_color_message(str, "green");
	g_free(str);

	set_cursor_waiting(TRUE);
	if (!com.script)
		control_window_switch_to_tab(OUTPUT_LOGS);

	if (!com.wd) {
		siril_log_message(_("Convert: no working directory set.\n"));
		set_cursor_waiting(FALSE);
		return 1;
	}

	struct _convert_data *args = malloc(sizeof(struct _convert_data));
	args->start = idx;
	args->list = files_to_link;
	args->total = count;
	if (output == SEQ_REGULAR)
		args->destroot = format_basename(destroot, TRUE);
	else
		args->destroot = destroot;
	args->input_has_a_seq = FALSE;
	args->input_has_a_film = FALSE;
	args->debayer = debayer;
	args->multiple_output = FALSE;
	args->output_type = output;
	args->make_link = make_link;
	gettimeofday(&(args->t_start), NULL);
	start_in_new_thread(convert_thread_worker, args);

	return 0;
}

int process_register(int nb) {
	struct registration_args *reg_args;
	struct registration_method *method;
	char *msg;
	int i;

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return 1;
	}
	sequence *seq = load_sequence(word[1], NULL);
	if (!seq)
		return 1;

	/* getting the selected registration method */
	method = malloc(sizeof(struct registration_method));
	method->name = strdup(_("Global Star Alignment (deep-sky)"));
	method->method_ptr = &register_star_alignment;
	method->sel = REQUIRES_NO_SELECTION;
	method->type = REGTYPE_DEEPSKY;

	reg_args = calloc(1, sizeof(struct registration_args));

	if (!com.script)
		control_window_switch_to_tab(OUTPUT_LOGS);

	/* filling the arguments for registration */
	reg_args->func = method->method_ptr;
	reg_args->seq = seq;
	reg_args->reference_image = sequence_find_refimage(seq);
	reg_args->process_all_frames = TRUE;
	reg_args->follow_star = FALSE;
	reg_args->matchSelection = FALSE;
	reg_args->translation_only = FALSE;
	reg_args->x2upscale = FALSE;
	reg_args->prefix = "r_";
	reg_args->min_pairs = AT_MATCH_MINPAIRS;

	/* check for options */
	for (i = 2; i < nb; i++) {
		if (word[i]) {
			if (!strcmp(word[i], "-drizzle")) {
				reg_args->x2upscale = TRUE;
			} else if (!strcmp(word[i], "-norot")) {
				reg_args->translation_only = TRUE;
			} else if (g_str_has_prefix(word[i], "-prefix=")) {
				char *current = word[i], *value;
				value = current + 8;
				if (value[0] == '\0') {
					siril_log_message(_("Missing argument to %s, aborting.\n"), current);
					return 1;
				}
				reg_args->prefix = strdup(value);
			} else if (g_str_has_prefix(word[i], "-minpairs=")) {
				char *current = word[i], *value;
				value = current + 10;
				if (value[0] == '\0') {
					siril_log_message(_("Missing argument to %s, aborting.\n"), current);
					return 1;
				}
				if (g_ascii_strtoull(value, NULL, 10) < AT_MATCH_MINPAIRS) {
					gchar *str = ngettext("%d smaller than minimum allowable star pairs: %d, aborting.\n", "%d smaller than minimum allowable star pairs: %d, aborting.\n",
							g_ascii_strtoull(value, NULL, 10));
					str = g_strdup_printf(str, g_ascii_strtoull(value, NULL, 10), AT_MATCH_MINPAIRS);
					siril_log_message(str);
					g_free(str);

					return 1;
				}
				reg_args->min_pairs = g_ascii_strtoull(value, NULL, 10);
			}
		}
	}

	// testing free space
	if (reg_args->x2upscale ||
			(method->method_ptr == register_star_alignment &&
			 !reg_args->translation_only)) {
		// first, remove the files that we are about to create
		remove_prefixed_sequence_files(reg_args->seq, reg_args->prefix);

		int nb_frames = reg_args->process_all_frames ? reg_args->seq->number : reg_args->seq->selnum;
		int64_t size = seq_compute_size(reg_args->seq, nb_frames, get_data_type(seq->bitpix));
		if (reg_args->x2upscale)
			size *= 4;
		if (test_available_space(size) > 0) {
			free(reg_args);
			free(method);
			return 1;
		}
	}

	/* getting the selected registration layer from the combo box. The value is the index
	 * of the selected line, and they are in the same order than layers so there should be
	 * an exact matching between the two */
	reg_args->layer = (reg_args->seq->nb_layers == 3) ? 1 : 0;
	reg_args->interpolation = OPENCV_CUBIC;
	get_the_registration_area(reg_args, method);	// sets selection
	reg_args->run_in_thread = TRUE;
	reg_args->load_new_sequence = FALSE;	// don't load it for command line execution

	msg = siril_log_color_message(
			_("Registration: processing using method: %s\n"), "green",
			method->name);
	free(method);
	msg[strlen(msg) - 1] = '\0';
	set_progress_bar_data(msg, PROGRESS_RESET);

	set_cursor_waiting(TRUE);

	start_in_new_thread(register_thread_func, reg_args);
	return 0;
}

// parse normalization and filters from the stack command line, starting at word `first'
static int parse_stack_command_line(struct stacking_configuration *arg, int first, gboolean norm_allowed, gboolean out_allowed) {
	while (word[first]) {
		char *current = word[first], *value;
		if (!strcmp(current, "-nonorm") || !strcmp(current, "-no_norm"))
			arg->force_no_norm = TRUE;
		else if (!strcmp(current, "-output_norm")) {
			arg->output_norm = TRUE;
		} else if (!strcmp(current, "-weighted")) {
			if (arg->method != stack_mean_with_rejection) {
				siril_log_message(_("Weighting is allowed only with average stacking, ignoring.\n"));
			} else if (arg->norm == NO_NORM) {
				siril_log_message(_("Weighting is allowed only if normalization has been activated, ignoring.\n"));
			} else{
				arg->apply_weight = TRUE;
			}
		} else if (g_str_has_prefix(current, "-norm=")) {
			if (!norm_allowed) {
				siril_log_message(_("Normalization options are not allowed in this context, ignoring.\n"));
			} else {
				value = current + 6;
				if (!strcmp(value, "add"))
					arg->norm = ADDITIVE;
				else if (!strcmp(value, "addscale"))
					arg->norm = ADDITIVE_SCALING;
				else if (!strcmp(value, "mul"))
					arg->norm = MULTIPLICATIVE;
				else if (!strcmp(value, "mulscale"))
					arg->norm = MULTIPLICATIVE_SCALING;
			}
		} else if (g_str_has_prefix(current, "-filter-fwhm=")) {
			value = strchr(current, '=') + 1;
			if (value[0] != '\0') {
				char *end;
				float val = strtof(value, &end);
				if (end == value) {
					siril_log_message(_("Could not parse argument `%s' to the filter `%s', aborting.\n"), value, current);
					return 1;
				}
				if (*end == '%')
					arg->f_fwhm_p = val;
				else arg->f_fwhm = val;
			} else {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				return 1;
			}
		} else if (g_str_has_prefix(current, "-filter-wfwhm=")) {
			value = strchr(current, '=') + 1;
			if (value[0] != '\0') {
				char *end;
				float val = strtof(value, &end);
				if (end == value) {
					siril_log_message(_("Could not parse argument `%s' to the filter `%s', aborting.\n"), value, current);
					return 1;
				}
				if (*end == '%')
					arg->f_wfwhm_p = val;
				else arg->f_wfwhm = val;
			} else {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				return 1;
			}
		} else if (g_str_has_prefix(current, "-filter-round=") ||
				g_str_has_prefix(current, "-filter-roundness=")) {
			value = strchr(current, '=') + 1;
			if (value[0] != '\0') {
				char *end;
				float val = strtof(value, &end);
				if (end == value) {
					siril_log_message(_("Could not parse argument `%s' to the filter `%s', aborting.\n"), value, current);
					return 1;
				}
				if (*end == '%')
					arg->f_round_p = val;
				else arg->f_round = val;
			} else {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				return 1;
			}
		} else if (g_str_has_prefix(current, "-filter-qual=") ||
				g_str_has_prefix(current, "-filter-quality=")) {
			value = strchr(current, '=') + 1;
			if (value[0] != '\0') {
				char *end;
				float val = strtof(value, &end);
				if (end == value) {
					siril_log_message(_("Could not parse argument `%s' to the filter `%s', aborting.\n"), value, current);
					return 1;
				}
				if (*end == '%')
					arg->f_quality_p = val;
				else arg->f_quality = val;
			} else {
				siril_log_message(_("Missing argument to %s, aborting.\n"), current);
				return 1;
			}
		} else if (g_str_has_prefix(current, "-filter-incl") ||
				g_str_has_prefix(current, "-filter-included")) {
			arg->filter_included = TRUE;
		} else if (g_str_has_prefix(current, "-out=")) {
			if (out_allowed) {
				value = current + 5;
				if (value[0] == '\0') {
					siril_log_message(_("Missing argument to %s, aborting.\n"), current);
					return 1;
				}
				arg->result_file = strdup(value);
			}
			else {
				siril_log_message(_("Output filename option is not allowed in this context, ignoring.\n"));
			}
		}
		else {
			siril_log_message(_("Unexpected argument to stacking `%s', aborting.\n"), current);
			return 1;
		}
		first++;
	}
	return 0;
}

static int stack_one_seq(struct stacking_configuration *arg) {
	int retval = -1;
	sequence *seq = readseqfile(arg->seqfile);
	if (seq != NULL) {
		struct stacking_args args = { 0 };
		gchar *error = NULL;
		if (seq_check_basic_data(seq, FALSE) == -1) {
			free(seq);
			return 1;
		}
		siril_log_message(_("Stacking sequence %s\n"), seq->seqname);
		args.seq = seq;
		args.ref_image = sequence_find_refimage(seq);
		// the three below: used only if method is average w/ rejection
		if (arg->method == stack_mean_with_rejection && (arg->sig[0] != 0.0 || arg->sig[1] != 0.0)) {
			args.sig[0] = arg->sig[0];
			args.sig[1] = arg->sig[1];
			args.type_of_rejection = arg->type_of_rejection;
		} else {
			args.type_of_rejection = NO_REJEC;
			siril_log_message(_("Not using rejection for stacking\n"));
		}
		args.coeff.offset = NULL;
		args.coeff.mul = NULL;
		args.coeff.scale = NULL;
		if (!arg->force_no_norm &&
				(arg->method == stack_median || arg->method == stack_mean_with_rejection))
			args.normalize = arg->norm;
		else args.normalize = NO_NORM;
		args.method = arg->method;
		args.force_norm = FALSE;
		args.output_norm = arg->output_norm;
		args.reglayer = args.seq->nb_layers == 1 ? 0 : 1;
		args.apply_weight = arg->apply_weight;

		// manage filters
		if (convert_stack_data_to_filter(arg, &args) ||
				setup_filtered_data(&args)) {
			free_sequence(seq, TRUE);
			return 1;
		}
		args.description = describe_filter(seq, args.filtering_criterion, args.filtering_parameter);
		args.use_32bit_output = evaluate_stacking_should_output_32bits(args.method,
			args.seq, args.nb_images_to_stack, &error);
		if (error) {
			siril_log_color_message(error, "red");
			free_sequence(seq, TRUE);
			return 1;
		}

		if (!arg->result_file) {
			char filename[256];
			char *suffix = g_str_has_suffix(seq->seqname, "_") ||
				g_str_has_suffix(seq->seqname, "-") ? "" : "_";
			snprintf(filename, 256, "%s%sstacked%s",
					seq->seqname, suffix, com.pref.ext);
			arg->result_file = strdup(filename);
		}

		main_stack(&args);

		retval = args.retval;
		clean_end_stacking(&args);
		free_sequence(seq, TRUE);
		free(args.image_indices);
		free(args.description);

		if (!retval) {
			struct noise_data noise_args = { .fit = &gfit, .verbose = FALSE, .use_idle = FALSE };
			noise(&noise_args);
			if (savefits(arg->result_file, &gfit))
				siril_log_color_message(_("Could not save the stacking result %s\n"),
						"red", arg->result_file);
			clearfits(&gfit);
			++arg->number_of_loaded_sequences;
		}
		else if (!get_thread_run()) return -1;

	} else {
		siril_log_message(_("No sequence `%s' found.\n"), arg->seqfile);
	}
	return retval;
}

static gpointer stackall_worker(gpointer garg) {
	GDir *dir;
	GError *error = NULL;
	const gchar *file;
	struct timeval t_end;
	struct stacking_configuration *arg = (struct stacking_configuration *)garg;
	gboolean was_in_script = com.script;
	com.script = TRUE;

	siril_log_message(_("Looking for sequences in current working directory...\n"));
	if (check_seq(FALSE) || (dir = g_dir_open(com.wd, 0, &error)) == NULL) {
		siril_log_message(_("Error while searching sequences or opening the directory.\n"));
		if (error) {
			fprintf(stderr, "stackall: %s\n", error->message);
			g_clear_error(&error);
		}
		siril_add_idle(end_generic, NULL);
		return NULL;
	}

	siril_log_message(_("Starting stacking of found sequences...\n"));
	while ((file = g_dir_read_name(dir)) != NULL) {
		if (g_str_has_suffix(file, ".seq")) {
			arg->seqfile = strdup(file);
			stack_one_seq(arg);

			g_free(arg->result_file);
			arg->result_file = NULL;
			g_free(arg->seqfile);
		}
	}

	siril_log_message(_("Stacked %d sequences successfully.\n"), arg->number_of_loaded_sequences);
	gettimeofday(&t_end, NULL);
	show_time(arg->t_start, t_end);
	g_dir_close(dir);
	free(arg);
	com.script = was_in_script;
	siril_add_idle(end_generic, NULL);
	return NULL;
}

int process_stackall(int nb) {
	struct stacking_configuration *arg;

	arg = calloc(1, sizeof(struct stacking_configuration));
	arg->f_fwhm = -1.f; arg->f_fwhm_p = -1.f; arg->f_round = -1.f;
	arg->f_round_p = -1.f; arg->f_quality = -1.f; arg->f_quality_p = -1.f;
	arg->filter_included = FALSE; arg->norm = NO_NORM; arg->force_no_norm = FALSE;
	arg->apply_weight = FALSE;

	// stackall { sum | min | max } [-filter-fwhm=value[%]] [-filter-wfwhm=value[%]] [-filter-round=value[%]] [-filter-quality=value[%]] [-filter-incl[uded]]
	// stackall { med | median } [-nonorm, norm=] [-filter-incl[uded]]
	// stackall { rej | mean } sigma_low sigma_high [-nonorm, norm=] [-filter-fwhm=value[%]] [-filter-round=value[%]] [-filter-quality=value[%]] [-filter-incl[uded]] [-weighted]
	if (!word[1]) {
		arg->method = stack_summing_generic;
	} else {
		int start_arg_opt = 2;
		gboolean allow_norm = FALSE;
		if (!strcmp(word[1], "sum")) {
			arg->method = stack_summing_generic;
		} else if (!strcmp(word[1], "max")) {
			arg->method = stack_addmax;
		} else if (!strcmp(word[1], "min")) {
			arg->method = stack_addmin;
		} else if (!strcmp(word[1], "med") || !strcmp(word[1], "median")) {
			arg->method = stack_median;
			allow_norm = TRUE;
		} else if (!strcmp(word[1], "rej") || !strcmp(word[1], "mean")) {
			int shift = 1;
			if (!strcmp(word[3], "p") || !strcmp(word[3], "percentile")) {
				arg->type_of_rejection = PERCENTILE;
			} else if (!strcmp(word[3], "s") || !strcmp(word[3], "sigma")) {
				arg->type_of_rejection = SIGMA;
			} else if (!strcmp(word[3], "m") || !strcmp(word[3], "median")) {
				arg->type_of_rejection = SIGMEDIAN;
			} else if (!strcmp(word[3], "l") || !strcmp(word[3], "linear")) {
				arg->type_of_rejection = LINEARFIT;
			} else if (!strcmp(word[3], "w") || !strcmp(word[3], "winsorized")) {
				arg->type_of_rejection = WINSORIZED;
			} else if (!strcmp(word[3], "g") || !strcmp(word[3], "generalized")) {
				arg->type_of_rejection = GESDT;
			} else {
				arg->type_of_rejection = WINSORIZED;
				shift = 0;
			}
			if (!word[2 + shift] || !word[3 + shift] || (arg->sig[0] = g_ascii_strtod(word[2 + shift], NULL)) < 0.0
					|| (arg->sig[1] = g_ascii_strtod(word[3 + shift], NULL)) < 0.0) {
				siril_log_color_message(_("The average stacking with rejection requires two extra arguments: sigma low and high.\n"), "red");
				goto failure;
			}
			arg->method = stack_mean_with_rejection;
			start_arg_opt = 4 + shift;
			allow_norm = TRUE;
		}
		else {
			siril_log_color_message(_("Stacking method type '%s' is invalid\n"), "red", word[2]);
			goto failure;
		}
		if (parse_stack_command_line(arg, start_arg_opt, allow_norm, FALSE))
			goto failure;
	}
	set_cursor_waiting(TRUE);
	if (!com.headless)
		control_window_switch_to_tab(OUTPUT_LOGS);
	gettimeofday(&arg->t_start, NULL);

	start_in_new_thread(stackall_worker, arg);
	return 0;

failure:
	g_free(arg->result_file);
	g_free(arg->seqfile);
	free(arg);
	return 1;
}

static gpointer stackone_worker(gpointer garg) {
	int retval = 0;
	struct timeval t_end;
	struct stacking_configuration *arg = (struct stacking_configuration *)garg;
	gboolean was_in_script = com.script;
	com.script = TRUE;

	retval = stack_one_seq(arg);

	if (retval) {
		if (retval == ST_ALLOC_ERROR) {
			siril_log_message(_("It looks like there is a memory allocation error, change memory settings and try to fix it.\n"));
		}
	} else {
		siril_log_message(_("Stacked sequence successfully.\n"));
	}

	gettimeofday(&t_end, NULL);
	show_time(arg->t_start, t_end);

	g_free(arg->result_file);
	g_free(arg->seqfile);
	free(arg);
	com.script = was_in_script;
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}

int process_stackone(int nb) {
	struct stacking_configuration *arg;

	arg = calloc(1, sizeof(struct stacking_configuration));
	arg->f_fwhm = -1.f; arg->f_fwhm_p = -1.f; arg->f_round = -1.f;
	arg->f_round_p = -1.f; arg->f_quality = -1.f; arg->f_quality_p = -1.f;
	arg->filter_included = FALSE; arg->norm = NO_NORM; arg->force_no_norm = FALSE;
	arg->apply_weight = FALSE;

	sequence *seq = load_sequence(word[1], &arg->seqfile);
	if (!seq)
		goto failure;

	// stack seqfilename { sum | min | max } [-filter-fwhm=value[%]] [-filter-wfwhm=value[%]] [-filter-round=value[%]] [-filter-quality=value[%]] [-filter-incl[uded]] -out=result_filename
	// stack seqfilename { med | median } [-nonorm, norm=] [-filter-incl[uded]] -out=result_filename
	// stack seqfilename { rej | mean } sigma_low sigma_high [-nonorm, norm=] [-filter-fwhm=value[%]] [-filter-round=value[%]] [-filter-quality=value[%]] [-filter-incl[uded]] [-weighted] -out=result_filename
	if (!word[2]) {
		arg->method = stack_summing_generic;
	} else {
		int start_arg_opt = 3;
		gboolean allow_norm = FALSE;
		if (!strcmp(word[2], "sum")) {
			arg->method = stack_summing_generic;
		} else if (!strcmp(word[2], "max")) {
			arg->method = stack_addmax;
		} else if (!strcmp(word[2], "min")) {
			arg->method = stack_addmin;
		} else if (!strcmp(word[2], "med") || !strcmp(word[2], "median")) {
			arg->method = stack_median;
			allow_norm = TRUE;
		} else if (!strcmp(word[2], "rej") || !strcmp(word[2], "mean")) {
			int shift = 1;
			if (!strcmp(word[3], "p") || !strcmp(word[3], "percentile")) {
				arg->type_of_rejection = PERCENTILE;
			} else if (!strcmp(word[3], "s") || !strcmp(word[3], "sigma")) {
				arg->type_of_rejection = SIGMA;
			} else if (!strcmp(word[3], "m") || !strcmp(word[3], "median")) {
				arg->type_of_rejection = SIGMEDIAN;
			} else if (!strcmp(word[3], "l") || !strcmp(word[3], "linear")) {
				arg->type_of_rejection = LINEARFIT;
			} else if (!strcmp(word[3], "w") || !strcmp(word[3], "winsorized")) {
				arg->type_of_rejection = WINSORIZED;
			} else if (!strcmp(word[3], "g") || !strcmp(word[3], "generalized")) {
				arg->type_of_rejection = GESDT;
			} else {
				arg->type_of_rejection = WINSORIZED;
				shift = 0;
			}
			if (!word[3 + shift] || !word[4 + shift] || (arg->sig[0] = g_ascii_strtod(word[3 + shift], NULL)) < 0.0
					|| (arg->sig[1] = g_ascii_strtod(word[4 + shift], NULL)) < 0.0) {
				siril_log_color_message(_("The average stacking with rejection requires two extra arguments: sigma low and high.\n"), "red");
				goto failure;
			}
			arg->method = stack_mean_with_rejection;
			start_arg_opt = 5 + shift;
			allow_norm = TRUE;
		}
		else {
			siril_log_color_message(_("Stacking method type '%s' is invalid\n"), "red", word[2]);
			goto failure;
		}
		if (parse_stack_command_line(arg, start_arg_opt, allow_norm, TRUE))
			goto failure;
	}
	set_cursor_waiting(TRUE);
	gettimeofday(&arg->t_start, NULL);
	if (!com.headless)
		control_window_switch_to_tab(OUTPUT_LOGS);

	start_in_new_thread(stackone_worker, arg);
	return 0;

failure:
	g_free(arg->result_file);
	g_free(arg->seqfile);
	free(arg);
	return 1;
}

int process_preprocess(int nb) {
	struct preprocessing_data *args;
	int i, retvalue = 0;

	if (word[1][0] == '\0') {
		return -1;
	}

	sequence *seq = load_sequence(word[1], NULL);
	if (!seq)
		return 1;

	args = calloc(1, sizeof(struct preprocessing_data));
	args->ppprefix = "pp_";
	
	/* checking for options */
	for (i = 2; i < nb; i++) {
		if (word[i]) {
			if (g_str_has_prefix(word[i], "-bias=")) {
				args->bias = calloc(1, sizeof(fits));
				if (!readfits(word[i] + 6, args->bias, NULL, !com.pref.force_to_16bit)) {
					args->use_bias = TRUE;
				} else {
					retvalue = 1;
					free(args->bias);
					break;
				}
			} else if (g_str_has_prefix(word[i], "-dark=")) {
				args->dark = calloc(1, sizeof(fits));
				if (!readfits(word[i] + 6, args->dark, NULL, !com.pref.force_to_16bit)) {
					args->use_dark = TRUE;
					args->use_cosmetic_correction = TRUE;
				} else {
					retvalue = 1;
					free(args->dark);
					break;
				}
			} else if (g_str_has_prefix(word[i], "-flat=")) {
				args->flat = calloc(1, sizeof(fits));
				if (!readfits(word[i] + 6, args->flat, NULL, !com.pref.force_to_16bit)) {
					args->use_flat = TRUE;
				} else {
					retvalue = 1;
					free(args->flat);
					break;
				}
			} else if (g_str_has_prefix(word[i], "-prefix=")) {
				char *current = word[i], *value;
				value = current + 8;
				if (value[0] == '\0') {
					siril_log_message(_("Missing argument to %s, aborting.\n"), current);
					retvalue = 1;
					break;
				}
				args->ppprefix = strdup(value);
			} else if (!strcmp(word[i], "-opt")) {
				args->use_dark_optim = TRUE;
			} else if (!strcmp(word[i], "-fix_xtrans")) {
				args->fix_xtrans = TRUE;
			} else if (!strcmp(word[i], "-cfa")) {
				args->is_cfa = TRUE;
			} else if (!strcmp(word[i], "-debayer")) {
				args->debayer = TRUE;
			} else if (!strcmp(word[i], "-stretch")) {
				siril_log_message(_("-stretch option is now deprecated.\n")); // TODO. Should we keep it only for compatibility?
			} else if (!strcmp(word[i], "-flip")) {
				siril_log_message(_("-flip option is now deprecated.\n")); // TODO. Should we keep it only for compatibility?
			} else if (!strcmp(word[i], "-equalize_cfa")) {
				args->equalize_cfa = TRUE;
			} else if (!strcmp(word[i], "-fitseq")) {
				args->output_seqtype = SEQ_FITSEQ;
			}
		}
	}

	if (retvalue) {
		free(args);
		return -1;
	}

	siril_log_color_message(_("Preprocessing...\n"), "green");
	gettimeofday(&args->t_start, NULL);
	args->seq = seq;
	args->is_sequence = TRUE;
	args->autolevel = TRUE;
	args->normalisation = 1.0f;	// will be updated anyway
	args->sigma[0] = -1.00; /* cold pixels: it is better to deactivate it */
	args->sigma[1] =  3.00; /* hot pixels */
	args->allow_32bit_output = (args->seq->type == SEQ_REGULAR
			|| args->seq->type == SEQ_FITSEQ) && !com.pref.force_to_16bit;

	// start preprocessing
	set_cursor_waiting(TRUE);

	start_sequence_preprocessing(args);
	return 0;
}

int process_set_32bits(int nb) {
	com.pref.force_to_16bit = word[0][3] == '1';
	if (com.pref.force_to_16bit)
		siril_log_message(_("16-bit per channel in processed images mode is active\n"));
	else siril_log_message(_("32-bit per channel in processed images mode is active\n"));
	writeinitfile();
	if (!com.headless)
		set_GUI_misc();
	return 0;
}

int process_set_compress(int nb) {
	gboolean compress = g_ascii_strtoull(word[1], NULL, 10) == 1;
	int method = 0;
	double q = 16.0, hscale= 4.0;

	if (compress) {
		if (!word[2] || !word[3] || (!g_str_has_prefix(word[2], "-type="))) {
			siril_log_message(_("Please specify the type of compression and quantization value.\n"));
			return 1;
		}
		gchar *comp = NULL;
		if (!g_ascii_strncasecmp(word[2] + 6, "rice", 4)) {
			method = RICE_COMP;
			comp = g_strdup("rice");
		} else if (!g_ascii_strncasecmp(word[2] + 6, "gzip1", 5)) {
			method = GZIP1_COMP;
			comp = g_strdup("GZIP1");
		} else if (!g_ascii_strncasecmp(word[2] + 6, "gzip2", 5))  {
			method = GZIP2_COMP;
			comp = g_strdup("GZIP2");
		} else if (!g_ascii_strncasecmp(word[2] + 6, "hcompress", 9)) {
			method = HCOMPRESS_COMP;
			if (!word[4]) {
				siril_log_message(_("Please specify the value of hcompress scale factor.\n"));
				g_free(comp);
				return 1;
			}
			hscale = g_ascii_strtod(word[4], NULL);
			comp = g_strdup_printf("hcompress (scale factor = %.2lf) ", hscale);
		} else {
			siril_log_message(_("Wrong type of compression. Choices are rice, gzip1, gzip2 or hcompress\n"));
			return 1;
		}
		if (!word[3]) {
			siril_log_message(_("Please specify the value of quantization.\n"));
			g_free(comp);
			return 1;
		}
		q = g_ascii_strtod(word[3], NULL);
		if (q == 0.0 && (method == RICE_COMP || (method == HCOMPRESS_COMP))) {
			siril_log_message(_("Quantization can only be equal to 0 for GZIP1 and GZIP2 algorithms.\n"));
			return 1;
		}
		siril_log_message(_("Compression enabled with the %s algorithm and a quantization value of %.2lf\n"), comp, q);
		g_free(comp);
	} else {
		siril_log_message(_("No compression enabled.\n"));
	}
	com.pref.comp.fits_enabled = compress;
	com.pref.comp.fits_method = method;
	com.pref.comp.fits_quantization = q;
	com.pref.comp.fits_hcompress_scale = hscale;
	if (!com.headless)
		set_GUI_compression();
	writeinitfile();
	return 0;
}

#ifdef _OPENMP
int process_set_cpu(int nb){
	int proc_in, proc_out, proc_max;

	proc_in = g_ascii_strtoull(word[1], NULL, 10);
	proc_max = omp_get_num_procs();
	if (proc_in > proc_max || proc_in < 1) {
		siril_log_message(_("Number of logical processors MUST be greater "
				"than 0 and lower or equal to %d.\n"), proc_max);
		return 1;
	}
	omp_set_num_threads(proc_in);

#pragma omp parallel
	{
		proc_out = omp_get_num_threads();
	}

	gchar *str = ngettext("Using now %d logical processor\n", "Using now %d logical processors\n", proc_out);
	str = g_strdup_printf(str, proc_out);
	siril_log_message(str);
	g_free(str);

	com.max_thread = proc_out;
	if (!com.headless)
		update_spinCPU(0);

	return 0;
}
#endif

int process_set_mem(int nb){
	double ratio = g_ascii_strtod(word[1], NULL);
	if (ratio < 0.05 || ratio > 4.0) {
		siril_log_message(_("The accepted range for the ratio of memory used for stacking is [0.05, 4], with values below the available memory recommended\n"));
		return 1;
	}
	if (ratio > 1.0) {
		siril_log_message(_("Setting the ratio of memory used for stacking above 1 will require the use of on-disk memory, which can be very slow and is unrecommended (%g requested)\n"), ratio);
	}
	com.pref.stack.memory_ratio = ratio;
	writeinitfile();
	siril_log_message(_("Usable memory for stacking changed to %g\n"), ratio);
	if (!com.headless)
		set_GUI_misc();
	return 0;
}

int process_help(int nb){
	control_window_switch_to_tab(OUTPUT_LOGS);

	command *current = commands;
	siril_log_message(_("********* LIST OF AVAILABLE COMMANDS *********\n"));
	while (current->process) {
		siril_log_message("%s\n", current->usage);
		current++;
	}
	siril_log_message(_("********* END OF THE LIST *********\n"));
	return 0;
}

int process_exit(int nb){
	gtk_main_quit();
	return 0;
}

int process_extract(int nb) {
	int Nbr_Plan, maxplan, mins;
	
	if (!single_image_is_loaded()) return 1;

	Nbr_Plan = g_ascii_strtoull(word[1], NULL, 10);

	mins = min (gfit.rx, gfit.ry);
	maxplan = log(mins) / log(2) - 2;

	if ( Nbr_Plan > maxplan ){
		siril_log_message(_("Wavelet: maximum number of plans for this image size is %d\n"),
				maxplan);
		return 1;
	}
	
	struct wavelets_filter_data *args = malloc(sizeof(struct wavelets_filter_data));

	args->Type = TO_PAVE_BSPLINE;
	args->Nbr_Plan = Nbr_Plan;
	args->fit = &gfit;
	start_in_new_thread(extract_plans, args);

	return 0;
}

int process_reloadscripts(int nb){
	return refresh_scripts(FALSE, NULL);
}


int process_requires(int nb) {
	gchar **version, **required;
	gint major, minor, micro;
	gint req_major, req_minor, req_micro;

	version = g_strsplit(PACKAGE_VERSION, ".", 3);
	required = g_strsplit(word[1], ".", 3);

	if (g_strv_length(required) != 3) {
		siril_log_color_message(_("Required version is not correct.\n"), "red");

		g_strfreev(version);
		g_strfreev(required);
		return 1;
	}

	major = g_ascii_strtoull(version[0], NULL, 10);
	minor = g_ascii_strtoull(version[1], NULL, 10);
	micro = g_ascii_strtoull(version[2], NULL, 10);

	req_major = g_ascii_strtoull(required[0], NULL, 10);
	req_minor = g_ascii_strtoull(required[1], NULL, 10);
	req_micro = g_ascii_strtoull(required[2], NULL, 10);

	g_strfreev(version);
	g_strfreev(required);

	if ((major > req_major || (major == req_major && minor > req_minor)
			|| (major == req_major && minor == req_minor && micro >= req_micro))) {
		// no need to output something in script conditions
		if (!com.script) {
			siril_log_message(_("The required version of Siril is ok.\n"));
		}
		return 0;
	} else {
		if (!com.script) {
			siril_log_color_message(_("A newer version of Siril is required, please update your version.\n"), "red");
		} else {
			siril_log_color_message(_("The script you are executing requires a newer version of Siril to run (%s), aborting.\n"), "red", word[1]);
		}
		return 1;
	}
}
