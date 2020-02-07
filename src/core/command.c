/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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
#include "io/sequence.h"
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "gui/PSF_list.h"
#include "gui/histogram.h"
#include "gui/plot.h"
#include "gui/progress_and_log.h"
#include "gui/image_display.h"
#include "gui/image_interactions.h"
#include "gui/sequence_list.h"
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
#include "algos/demosaicing.h"
#include "algos/quality.h"
#include "algos/noise.h"
#include "algos/statistics.h"
#include "algos/sorting.h"
#include "algos/geometry.h"
#include "opencv/opencv.h"
#include "stacking/stacking.h"
#include "stacking/sum.h"
#include "registration/registration.h"
#include "registration/matching/match.h"

#include "command.h"
#include "command_def.h"
#include "command_list.h"
#include "command_line_processor.h"

char *word[MAX_COMMAND_WORDS];	// NULL terminated


int process_load(int nb){
	char filename[256];
	int retval, i;
	
	strncpy(filename, word[1], 250);
	filename[250] = '\0';
	
	for (i = 1; i < nb - 1; ++i) {
		strcat(filename, " ");
		strcat(filename, word[i + 1]);
	}
	expand_home_in_filename(filename, 256);
	retval = open_single_image(filename);
	return retval;
}

int process_satu(int nb){
	if (get_thread_run()) {
		siril_log_message(_("Another task is already in progress, ignoring new request.\n"));
		return 1;
	}
	if (!single_image_is_loaded()) return 1;

	struct enhance_saturation_data *args = malloc(sizeof(struct enhance_saturation_data));
	
	args->coeff = atof(word[1]);
	if (args->coeff == 0.0) args->coeff = 1.0;

	args->input = &gfit;
	args->output = &gfit;
	args->h_min = 0.0;
	args->h_max = 360.0;
	args->preserve = TRUE;

	set_cursor_waiting(TRUE);
	start_in_new_thread(enhance_saturation, args);

	return 0;
}

int process_save(int nb){
	gchar *filename;
	int retval;
	
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

	filename = g_strdup(word[1]);
	set_cursor_waiting(TRUE);
	retval = savefits(filename, &(gfit));
	set_precision_switch();
	set_cursor_waiting(FALSE);
	g_free(filename);
	return retval;
}

int process_savebmp(int nb){
	gchar *filename;
	
	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

	filename = g_strdup_printf("%s.bmp", word[1]);

	set_cursor_waiting(TRUE);
	savebmp(filename, &(gfit));
	set_cursor_waiting(FALSE);
	g_free(filename);
	return 0;
}

#ifdef HAVE_LIBJPEG
int process_savejpg(int nb){
	gchar *filename;

	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

	int quality = 100;
	
	if ((nb == 3) && atoi(word[2]) <= 100 && atoi(word[2]) > 0)
		quality = atoi(word[2]);

	filename = g_strdup_printf("%s.jpg", word[1]);

	set_cursor_waiting(TRUE);
	savejpg(filename, &gfit, quality);
	set_cursor_waiting(FALSE);
	g_free(filename);
	return 0;
}
#endif

#ifdef HAVE_LIBPNG
int process_savepng(int nb){
	gchar *filename;
	uint32_t bytes_per_sample;

	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

	filename = g_strdup_printf("%s.png", word[1]);

	set_cursor_waiting(TRUE);
	bytes_per_sample = gfit.orig_bitpix != BYTE_IMG ? 2 : 1;
	savepng(filename, &gfit, bytes_per_sample, gfit.naxes[2] == 3);
	set_cursor_waiting(FALSE);
	g_free(filename);
	return 0;
}
#endif

#ifdef HAVE_LIBTIFF
int process_savetif(int nb){
	gchar *filename;
	uint16 bitspersample = 16;

	if (!(single_image_is_loaded() || sequence_is_loaded()))
		return 1;

	if (strcasecmp(word[0], "savetif8") == 0)
		bitspersample = 8;
	filename = g_strdup_printf("%s.tif", word[1]);
	set_cursor_waiting(TRUE);
	savetif(filename, &gfit, bitspersample);
	set_cursor_waiting(FALSE);
	g_free(filename);
	return 0;
}
#endif

int process_savepnm(int nb){
	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

	saveNetPBM(word[1], &gfit);
	return 0;
}

int process_imoper(int nb){
	fits fit = { 0 };
	if (!single_image_is_loaded()) return 1;
	if (readfits(word[1], &fit, NULL)) return -1;
	imoper(&gfit, &fit, word[0][1], TRUE);

	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_addmax(int nb){
	fits fit = { 0 };

	if (!single_image_is_loaded()) return 1;

	if (readfits(word[1], &fit, NULL))
		return -1;
	if (addmax(&gfit, &fit) == 0) {
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	}
	return 0;
}

int process_fdiv(int nb){
	// combines an image division and a scalar multiplication.
	float norm;
	fits fit = { 0 };
	if (!single_image_is_loaded()) return 1;

	norm = atof(word[2]);
	if (readfits(word[1], &fit, NULL)) return -1;
	siril_fdiv(&gfit, &fit, norm, TRUE);

	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_fmul(int nb){
	float coeff;
	if (!single_image_is_loaded()) return 1;

	coeff = atof(word[1]);
	if (coeff <= 0.0) {
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
	double e;

	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

	if (com.selection.w > 0 && com.selection.h > 0) {
		memcpy(&area, &com.selection, sizeof(rectangle));
		e = entropy(&gfit, com.cvport, &area, NULL);
	}
	else {
		e = entropy(&gfit, com.cvport, NULL, NULL);
	}
	siril_log_message(_("Entropy: %.3lf\n"), e);
	return 0;
}

int process_gauss(int nb){
	if (!single_image_is_loaded()) return 1;

	unsharp(&(gfit), atof(word[1]), 0.0, TRUE);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_grey_flat(int nb) {
	if (!single_image_is_loaded()) return 1;

	compute_grey_flat(&gfit);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();

	return 0;
}

int process_rl(int nb) {
	double sigma, corner;
	int threshold, iter;

	if (!single_image_is_loaded()) return 1;

	if (!com.script)
		control_window_switch_to_tab(OUTPUT_LOGS);
	threshold = atoi(word[1]);
	sigma = atof(word[2]);
	corner = atof(word[3]);
	iter = atoi(word[4]);

	if (threshold < 0 || threshold > 200) {
		siril_log_message(_("threshold must be between [0, 200].\n"));
		return 1;
	}

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
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return 1;
	}

	struct deconv_data *args = malloc(sizeof(struct deconv_data));

	args->fit = &gfit;
	args->contrast_threshold = (size_t)threshold;
	args->sigma = sigma;
	args->corner_radius = corner;
	args->iterations = (size_t)iter;
	args->auto_limit = TRUE;

	set_cursor_waiting(TRUE);

	start_in_new_thread(RTdeconv, args);

	return 0;
}

int process_unsharp(int nb) {
	if (!single_image_is_loaded()) return 1;

	unsharp(&(gfit), atof(word[1]), atof(word[2]), TRUE);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_crop(int nb) {
	if (!single_image_is_loaded()) return 1;

	rectangle area;
	if ((!com.selection.h) || (!com.selection.w)) {
		if (nb == 5) {
			if (atoi(word[1]) < 0 || atoi(word[2]) < 0) {
				siril_log_message(_("Crop: x and y must be positive values.\n"));
				return 1;
			}
			if (atoi(word[3]) <= 0 || atoi(word[4]) <= 0) {
				siril_log_message(_("Crop: width and height must be greater than 0.\n"));
				return 1;
			}
			if (atoi(word[1]) + atoi(word[3]) > gfit.rx || atoi(word[2]) + atoi(word[4]) > gfit.ry) {
				siril_log_message(_("Crop: width and height, respectively, must be less than %d and %d.\n"), gfit.rx,gfit.ry);
				return 1;
			}
			area.x = atoi(word[1]);
			area.y = atoi(word[2]);
			area.w = atoi(word[3]);
			area.h = atoi(word[4]);
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
	retval = changedir(filename, NULL);
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

	if (!single_image_is_loaded()) return 1;

	nb_chan = gfit.naxes[2];

	g_assert(nb_chan == 1 || nb_chan == 3);

	tmpdir = g_get_tmp_dir();

	for (i = 0; i < nb - 1; ++i) {
		coef[i] = atof(word[i + 1]);
	}

	for (i = 0; i < nb_chan; i++) {
		dir[i] = g_build_filename(tmpdir, File_Name_Transform[i], NULL);
		wavelet_reconstruct_file(dir[i], coef, gfit.pdata[i]);
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
	float *Imag;

	if (!single_image_is_loaded()) return 1;

	tmpdir = g_get_tmp_dir();

	Nbr_Plan = atoi(word[1]);
	Type_Transform = atoi(word[2]);
	
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

	Imag = f_vector_alloc (gfit.rx * gfit.ry);
	
	for (chan = 0; chan < nb_chan; chan++) {
		dir[chan] = g_build_filename(tmpdir, File_Name_Transform[chan], NULL);
		wavelet_transform_file (Imag, gfit.ry, gfit.rx, dir[chan], Type_Transform, Nbr_Plan, gfit.pdata[chan]);
		g_free(dir[chan]);
	}
	
	free (Imag);
	return 0;
}

int process_log(int nb){
	if (!single_image_is_loaded()) return 1;

	loglut(&gfit);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_asinh(int nb) {
	if (!single_image_is_loaded()) return 1;

	double beta = atof(word[1]);

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

	if (!single_image_is_loaded()) return 1;

	if (!com.script)
		control_window_switch_to_tab(OUTPUT_LOGS);
	clip_limit = atof(word[1]);

	if (clip_limit <= 0.0) {
		siril_log_message(_("Clip limit must be > 0.\n"));
		return 1;
	}

	size = atoi(word[2]);

	if (size <= 0.0) {
		siril_log_message(_("Tile size must be > 0.\n"));
		return 1;
	}

	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return 1;
	}

	struct CLAHE_data *args = malloc(sizeof(struct CLAHE_data));

	args->fit = &gfit;
	args->clip = clip_limit;

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

int	process_mirrorx(int nb){
	if (!single_image_is_loaded()) return 1;

	mirrorx(&gfit, TRUE);
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int	process_mirrory(int nb){
	if (!single_image_is_loaded()) return 1;

	mirrory(&gfit, TRUE);
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_resample(int nb) {
	if (!single_image_is_loaded()) return 1;

	double factor = atof(word[1]);
	if (factor > 5.0) {
		siril_log_message(_("The scaling factor must be less than 5.0\n"));
		return 1;
	}
	int toX = round_to_int(factor * gfit.rx);
	int toY = round_to_int(factor * gfit.ry);
	
	set_cursor_waiting(TRUE);
	verbose_resize_gaussian(&gfit, toX, toY, OPENCV_LINEAR);
	
	adjust_vport_size_to_image();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	return 0;
}

int process_rgradient(int nb) {
	if (get_thread_run()) {
		siril_log_message(_("Another task is already in progress, ignoring new request.\n"));
		return 1;
	}

	if (!single_image_is_loaded()) return 1;

	struct rgradient_filter_data *args = malloc(sizeof(struct rgradient_filter_data));
	args->xc = atof(word[1]);
	args->yc = atof(word[2]);
	args->dR = atof(word[3]);
	args->da = atof(word[4]);
	args->fit = &gfit;

	set_cursor_waiting(TRUE);

	start_in_new_thread(rgradient_filter, args);
	return 0;
}

int process_rotate(int nb) {
	double degree;
	
	if (!single_image_is_loaded()) return 1;

	set_cursor_waiting(TRUE);
	degree = atof(word[1]);
	verbose_rotate_image(&gfit, degree, OPENCV_LINEAR, 1);	//INTER_LINEAR
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	return 0;
}

int process_rotatepi(int nb){
	if (!single_image_is_loaded()) return 1;

	verbose_rotate_image(&gfit, 180.0, OPENCV_LINEAR, 1);

	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_set_mag(int nb) {
	if (!single_image_is_loaded()) return 1;

	int layer = match_drawing_area_widget(com.vport[com.cvport], FALSE);
	double mag = atof(word[1]);

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
	if (single_image_is_loaded()) return 1;
	int n = atoi(word[1]) - 1;
	if (n < 0 || n > com.seq.number) {
		siril_log_message(_("The reference image must be set between 1 and %d\n"), com.seq.number);
		return 1;
	}
	if (sequence_is_loaded()) {
		free_reference_image();

		com.seq.reference_image = n;
		test_and_allocate_reference_image(-1);
		// a reference image should not be excluded to avoid confusion
		if (!com.seq.imgparam[com.seq.current].incl) {
			toggle_image_selection(com.seq.current);
		}

		if (!com.script) {
			sequence_list_change_reference();
			update_stack_interface(FALSE);// get stacking info and enable the Go button
			adjust_sellabel();	// reference image is named in the label
			drawPlot();		// update plots
		}
		writeseqfile(&com.seq);
	}

	return 0;
}

int process_unset_mag(int nb) {
	com.magOffset = 0.0;
	return 0;
}

int process_set_mag_seq(int nb) {
	if (!sequence_is_loaded()) {
		siril_log_message(_("This command can be used only when a sequence is loaded\n"));
		return 1;
	}
	double mag = atof(word[1]);
	int i;
	for (i = 0; i < MAX_SEQPSF && com.seq.photometry[i]; i++);
	com.seq.reference_star = i-1;
	if (i == 0) {
		siril_log_message(_("Run a PSF for the sequence first (see seqpsf)\n"));
		return 1;
	}
	com.seq.reference_mag = mag;
	siril_log_message(_("Reference magnitude has been set for star %d to %f and will be computed for each image\n"), i-1, mag);
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

		free(com.ext);
		str = g_string_new(".");
		str = g_string_append(str, word[1]);
		str = g_string_ascii_down(str);
		com.ext = g_string_free(str, FALSE);
		writeinitfile();
	}

	return 0;
}

int process_set_findstar(int nb) {
	double sigma = atof(word[1]);
	double roundness = atof(word[2]);
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
		siril_log_message(_("This command can be used only when a sequence is loaded\n"));
		return 1;
	}
	com.seq.reference_star = -1;
	com.seq.reference_mag = -1001.0;
	siril_log_message(_("Reference magnitude unset for sequence\n"));
	drawPlot();
	return 0;
}

int process_psf(int nb){
	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

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
		siril_log_message(_("Another task is already in progress, ignoring new request.\n"));
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
		siril_log_message(_("This command can be used only when a sequence is loaded\n"));
		return 1;
	}
}

int process_seq_crop(int nb) {
	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return 1;
	}

	rectangle area;
	sequence *seq;
	gchar *file;

	if ((!com.selection.h) || (!com.selection.w)) {
		if (nb == 6) {
			if (atoi(word[2]) < 0 || atoi(word[3]) < 0) {
				siril_log_message(_("Crop: x and y must be positive values.\n"));
				return 1;
			}
			if (atoi(word[4]) <= 0 || atoi(word[5]) <= 0) {
				siril_log_message(_("Crop: width and height must be greater than 0.\n"));
				return 1;
			}
			area.x = atoi(word[2]);
			area.y = atoi(word[3]);
			area.w = atoi(word[4]);
			area.h = atoi(word[5]);
		}
		else {
			siril_log_message(_("Crop: select a region or provide x, y, width, height\n"));
			return 1;
		}
	} else {
		memcpy(&area, &com.selection, sizeof(rectangle));
	}

	file = g_strdup(word[1]);
	if (!ends_with(file, ".seq")) {
		str_append(&file, ".seq");
	}

	if (!existseq(file)) {
		if (check_seq(FALSE)) {
			siril_log_message(_("No sequence `%s' found.\n"), file);
			return 1;
		}
	}
	seq = readseqfile(file);
	if (seq == NULL) {
		siril_log_message(_("No sequence `%s' found.\n"), file);
		return 1;
	}
	if (seq_check_basic_data(seq, FALSE) == -1) {
		free(seq);
		return 1;
	}
	if (atoi(word[4]) > seq->rx || atoi(word[5]) > seq->ry) {
		siril_log_message(_("Crop: width and height, respectively, must be less than %d and %d.\n"),
				seq->rx, seq->ry);
		return 1;
	}

	struct crop_sequence_data *args = malloc(sizeof(struct crop_sequence_data));

	args->seq = seq;
	args->area = area;
	args->prefix = "cropped_";

	set_cursor_waiting(TRUE);

	start_in_new_thread(crop_sequence, args);
	return 0;
}

int process_bg(int nb){
	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

	WORD bg = round_to_WORD(background(&gfit, -1, &com.selection, TRUE));
	siril_log_message(_("Background value: %d\n"), bg);
	return 0;
}

int process_bgnoise(int nb){
	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return 1;
	}

	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

	struct noise_data *args = malloc(sizeof(struct noise_data));

	if (!com.script) {
		control_window_switch_to_tab(OUTPUT_LOGS);
		set_cursor_waiting(TRUE);
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
	size_t i;
	int nlayer = atoi(word[1]);
	char* clayer;
	char name [20];
	
	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

	if (nlayer>3 || nlayer <0)
		return 1;
	gsl_histogram* histo = computeHisto(&gfit, nlayer);
	if (!isrgb(&gfit)) clayer = strdup("bw");		//if B&W
	else clayer = vport_number_to_name(nlayer);
	snprintf(name, 20, "histo_%s.dat",clayer);

	FILE *f = g_fopen(name, "w");

	if (f == NULL) {
		free(clayer);
		return 1;
	}
	for (i = 0; i < USHRT_MAX + 1; i++)
		fprintf(f, "%zu %d\n", i, (int) gsl_histogram_get (histo, i));
	fclose(f);
	gsl_histogram_free(histo);
	siril_log_message(_("The file %s has been created for the %s layer.\n"), name, clayer);
	free(clayer);
	return 0;
}

int process_thresh(int nb){
	int lo, hi;

	if (!single_image_is_loaded()) return 1;

	lo = atoi(word[1]);
	hi = atoi(word[2]);
	threshlo(&gfit, lo);
	threshhi(&gfit, hi);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_threshlo(int nb){
	int lo;

	if (!single_image_is_loaded()) return 1;

	lo = atoi(word[1]);
	threshlo(&gfit, lo);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_threshhi(int nb){
	int hi;

	if (!single_image_is_loaded()) return 1;

	hi = atoi(word[1]);
	threshhi(&gfit, hi);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_nozero(int nb){
	int level;

	if (!single_image_is_loaded()) return 1;

	level = atoi(word[1]);
	nozero(&gfit, level);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_ddp(int nb){
	// combines an image division and a scalar multiplication.
	float coeff, sigma;
	unsigned level;

	if (!single_image_is_loaded()) return 1;

	level = atoi(word[1]);
	coeff = atof(word[2]);
	sigma = atof(word[3]);
	ddp(&gfit, level, coeff, sigma);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_new(int nb){
	int width, height, layers;
	
	width = atof(word[1]);
	height = atof(word[2]);
	layers = atoi(word[3]);
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
	com.uniq->nb_layers = gfit.naxes[2];
	com.uniq->layers = calloc(com.uniq->nb_layers, sizeof(layer_info));
	com.uniq->fit = &gfit;

	open_single_image_from_gfit();
	return 0;
}

int process_visu(int nb){
	int low, high;
	
	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

	low = atoi(word[1]);
	high = atoi(word[2]);
	if ((high > USHRT_MAX) || (low < 0)) {
		siril_log_message(_("Values must be positive and less than %d.\n"), USHRT_MAX);
		return 1;
	}
	visu(&gfit, low, high);
	return 0;
}

int process_fill2(int nb){
	int level = atoi(word[1]);
	rectangle area;

	if (!single_image_is_loaded()) return 1;

	if ((!com.selection.h) || (!com.selection.w)) {
		if (nb == 6) {
			area.x = atoi(word[2]);
			area.y = atoi(word[3]);
			area.w = atoi(word[4]);
			area.h = atoi(word[5]);
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
	int layer = RLAYER;

	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

	if (isrgb(&gfit)) layer = GLAYER;
	delete_selected_area();
	com.stars = peaker(&gfit, layer, &com.starfinder_conf, &nbstars, NULL, TRUE);
	siril_log_message(_("Found %d stars in image, channel #%d\n"), nbstars, layer);
	if (com.stars)
		refresh_stars_list(com.stars);
	return 0;
}

int process_findhot(int nb){
	long icold, ihot;
	char filename[256];
	int i;
	char type;

	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

	if (gfit.naxes[2] != 1) {
		siril_log_message(_("find_hot must be applied on an one-channel master-dark frame"));
		return 1;
	}
	double sig[2];
	sig[0] = atof(word[2]);
	sig[1] = atof(word[3]);

	deviant_pixel *dev = find_deviant_pixels(&gfit, sig, &icold, &ihot);
	siril_log_message(_("%ld cold and %ld hot pixels\n"), icold, ihot);

	sprintf(filename, "%s.lst", word[1]);
	FILE *cosme_file = g_fopen(filename, "w");
	if (cosme_file == NULL) {
		siril_log_message(_("Cannot open file: %s\n"), filename);
		free(dev);
		return 1;
	}

	for (i = 0; i < icold + ihot; i++) {
		int y = gfit.ry - (int) dev[i].p.y - 1;  /* FITS is stored bottom to top */
		if (dev[i].type == HOT_PIXEL)
			type = 'H';
		else
			type = 'C';
		fprintf(cosme_file, "P %d %d %c\n", (int) dev[i].p.x, y, type);
	}

	free(dev);
	fclose(cosme_file);

	return 0;
}

int process_cosme(int nb) {
	FILE* cosme_file = NULL;
	deviant_pixel dev;
	char *filename;
	double dirty;
	int is_cfa, i = 0, retval = 0;
	int nb_tokens;
	char line[64];
	char type;

	if (!single_image_is_loaded()) return 1;

	if (!ends_with(word[1], ".lst")) {
		filename = g_strdup_printf("%s.lst", word[1]);
	} else {
		filename = g_strdup(word[1]);
	}
	cosme_file = g_fopen(filename, "r");
	if (cosme_file == NULL) {
		siril_log_message(_("Cannot open file: %s\n"), filename);
		g_free(filename);
		return 1;
	}
	g_free(filename);
	if (word[0][5] == '_')
		is_cfa = 1;
	else
		is_cfa = 0;

	while (fgets(line, 63, cosme_file)) {
		++i;
		switch (line[0]) {
		case '#': // comments.
			continue;
			break;
		case 'P':
			nb_tokens = sscanf(line + 2, "%lf %lf %c", &dev.p.x, &dev.p.y, &type);
			if (nb_tokens != 2 && nb_tokens != 3) {
				fprintf(stderr, "cosmetic correction: "
						"cosme file format error at line %d: %s", i, line);
				retval = 1;
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
	}

	fclose(cosme_file);
	if (retval)
		siril_log_message(_("There were some errors, please check your input file.\n"));

	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_fmedian(int nb){
	if (get_thread_run()) {
		siril_log_message(_("Another task is already in progress, ignoring new request.\n"));
		return 1;
	}
	
	if (!single_image_is_loaded()) return 1;

	struct median_filter_data *args = malloc(sizeof(struct median_filter_data));
	args->ksize = atoi(word[1]);
	args->amount = atof(word[2]);
	args->iterations = 1;
	
	if (!(args->ksize & 1) || args->ksize < 2) {
		siril_log_message(_("The size of the kernel MUST be odd and greater than 1.\n"));
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
	double x_avg, y_avg;

	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

	if (!FindCentre(&gfit, &x_avg, &y_avg)) {
		y_avg = gfit.ry - y_avg;	// FITS are stored bottom to top
		siril_log_message(_("Center of gravity coordinates are (%.3lf, %.3lf)\n"), x_avg, y_avg);
		return 0;
	}
	return 1;
}

int process_clear(int nb) {
	if (com.script) return 0;
	GtkTextView *text = GTK_TEXT_VIEW(gtk_builder_get_object(builder, "output"));
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(text);
	GtkTextIter start_iter, end_iter;
	gtk_text_buffer_get_start_iter(tbuf, &start_iter);
	gtk_text_buffer_get_end_iter(tbuf, &end_iter);
	gtk_text_buffer_delete(tbuf, &start_iter, &end_iter);
	return 0;
}

int process_clearstar(int nb){
	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

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
	
	if (!single_image_is_loaded()) return 1;

	if ((!com.selection.h) || (!com.selection.w)) {
		if (nb == 6) {
			area.x = atoi(word[2]);
			area.y = atoi(word[3]);
			area.w = atoi(word[4]);
			area.h = atoi(word[5]);
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
	level = atoi(word[1]);
	fill(&gfit, level, &area);
	redraw(com.cvport, REMAP_ALL);
	return 0;
}

int process_offset(int nb){
	int level;
	
	if (!single_image_is_loaded()) return 1;

	level = atof(word[1]);
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
		siril_log_message(_("Another task is already in progress, ignoring new request.\n"));
		return 1;
	}

	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;
	if (gfit.naxes[2] == 1) return 1;

	struct scnr_data *args = malloc(sizeof(struct scnr_data));
	
	args->type = atoi(word[1]);
	args->fit = &gfit;
	args->amount = 0.0;
	args->preserve = TRUE;

	set_cursor_waiting(TRUE);

	start_in_new_thread(scnr, args);

	return 0;
}

int process_fft(int nb){
	if (get_thread_run()) {
		siril_log_message(_("Another task is already in progress, ignoring new request.\n"));
		return 1;
	}

	if (sequence_is_loaded()) {
		siril_log_message(_("FFT does not work with sequences\n"));
		return 1;
	}

	if (!single_image_is_loaded()) return 1;

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
		siril_log_message(_("Another task is already in progress, ignoring new request.\n"));
		return 1;
	}

	if (!single_image_is_loaded()) return 1;

	struct banding_data *args = malloc(sizeof(struct banding_data));

	args->amount = atof(word[1]);
	args->sigma = atof(word[2]);
	args->protect_highlights = TRUE;
	args->fit = &gfit;

	set_cursor_waiting(TRUE);

	start_in_new_thread(BandingEngineThreaded, args);
	
	return 0;
}


int process_findcosme(int nb) {
	gboolean is_sequence;
	sequence *seq = NULL;
	int i = 0;

	if (get_thread_run()) {
		siril_log_message(_("Another task is "
				"already in progress, ignoring new request.\n"));
		return 1;
	}

	is_sequence = (word[0][0] == 's');

	if (is_sequence) {
		gchar *file = g_strdup(word[1]);
		if (!ends_with(file, ".seq")) {
			str_append(&file, ".seq");
		}

		if (!existseq(file)) {
			if (check_seq(FALSE)) {
				siril_log_message(_("No sequence `%s' found.\n"), file);
				return 1;
			}
		}
		seq = readseqfile(file);
		if (seq == NULL) {
			siril_log_message(_("No sequence `%s' found.\n"), file);
			return 1;
		}
		if (seq_check_basic_data(seq, FALSE) == -1) {
			free(seq);
			return 1;
		}
		i++;
	} else {
		if (!single_image_is_loaded()) return 1;
	}

	struct cosmetic_data *args = malloc(sizeof(struct cosmetic_data));

	args->seq = seq;
	args->sigma[0] = atof(word[1 + i]);
	args->sigma[1] = atof(word[2 + i]);
	args->is_cfa = (word[0][10] == '_' || word[0][13] == '_');	// find_cosme_cfa or seqfind_cosme_cfa
	args->fit = &gfit;

	set_cursor_waiting(TRUE);

	if (is_sequence) {
		args->seqEntry = "cc_";
		apply_cosmetic_to_sequence(args);
	} else {
		start_in_new_thread(autoDetectThreaded, args);
	}

	return 0;
}

int select_unselect(gboolean select) {
	if (!sequence_is_loaded()) {
		siril_log_message(_("Use this command to select images in a sequence, load a sequence first.\n"));
		return 1;
	}
	int from = atoi(word[1]);
	int to = atoi(word[2]);
	if (from < 0 || from >= com.seq.number) {
		siril_log_message(_("The first argument must be between 0 and the number of images minus one.\n"));
		return 1;
	}
	int i;
	gboolean current_updated = FALSE;
	for (i=from; i<=to; i++) {
		if (i >= com.seq.number) break;
		if (com.seq.imgparam[i].incl != select) {
			com.seq.imgparam[i].incl = select;
			if (!com.headless)
				sequence_list_change_selection_index(i);
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
	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

	if (!isrgb(&gfit)) {
		siril_log_message(_("Siril cannot split layers. Make sure your image is in RGB mode.\n"));
		return 1;
	}
	gchar *R = g_strdup_printf("%s%s", word[1], com.ext);
	gchar *G = g_strdup_printf("%s%s", word[2], com.ext);
	gchar *B = g_strdup_printf("%s%s", word[3], com.ext);

	save1fits16(R, &gfit, RLAYER);
	save1fits16(G, &gfit, GLAYER);
	save1fits16(B, &gfit, BLAYER);

	g_free(R);
	g_free(G);
	g_free(B);
	return 0;
}

int process_split_cfa(int nb) {
	if (isrgb(&gfit)) {
		siril_log_message(_("Siril cannot split CFA channel. Make sure your image is in CFA mode.\n"));
		return 1;
	}
	char *filename = NULL;

	fits f_cfa0 = { 0 };
	fits f_cfa1 = { 0 };
	fits f_cfa2 = { 0 };
	fits f_cfa3 = { 0 };

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

	gchar *cfa0 = g_strdup_printf("CFA0_%s%s", filename, com.ext);
	gchar *cfa1 = g_strdup_printf("CFA1_%s%s", filename, com.ext);
	gchar *cfa2 = g_strdup_printf("CFA2_%s%s", filename, com.ext);
	gchar *cfa3 = g_strdup_printf("CFA3_%s%s", filename, com.ext);

	int ret = split_cfa(&gfit, &f_cfa0, &f_cfa1, &f_cfa2, &f_cfa3);
	if (ret) {
		g_free(cfa0);
		g_free(cfa1);
		g_free(cfa2);
		g_free(cfa3);
		clearfits(&f_cfa0);
		clearfits(&f_cfa1);
		clearfits(&f_cfa2);
		clearfits(&f_cfa3);
		return ret;
	}

	save1fits16(cfa0, &f_cfa0, 0);
	save1fits16(cfa1, &f_cfa1, 0);
	save1fits16(cfa2, &f_cfa2, 0);
	save1fits16(cfa3, &f_cfa3, 0);

	g_free(cfa0);
	g_free(cfa1);
	g_free(cfa2);
	g_free(cfa3);

	clearfits(&f_cfa0);
	clearfits(&f_cfa1);
	clearfits(&f_cfa2);
	clearfits(&f_cfa3);

	free(filename);
	return 0;
}

int process_seq_split_cfa(int nb) {
	sequence *seq = NULL;

	if (get_thread_run()) {
		siril_log_message(_("Another task is "
				"already in progress, ignoring new request.\n"));
		return 1;
	}

	gchar *file = g_strdup(word[1]);
	if (!ends_with(file, ".seq")) {
		str_append(&file, ".seq");
	}

	if (!existseq(file)) {
		if (check_seq(FALSE)) {
			siril_log_message(_("No sequence `%s' found.\n"), file);
			return 1;
		}
	}
	seq = readseqfile(file);
	if (seq == NULL) {
		siril_log_message(_("No sequence `%s' found.\n"), file);
		return 1;
	}
	if (seq_check_basic_data(seq, FALSE) == -1) {
		free(seq);
		return 1;
	}

	struct split_cfa_data *args = malloc(sizeof(struct split_cfa_data));

	args->seq = seq;
	args->fit = &gfit;
	args->seqEntry = "CFA_";

	set_cursor_waiting(TRUE);

	apply_split_cfa_to_sequence(args);

	return 0;
}


int process_stat(int nb){
	int nplane;
	int layer;
	char layername[6];

	if (!(single_image_is_loaded() || sequence_is_loaded())) return 1;

	nplane = gfit.naxes[2];

	for (layer = 0; layer < nplane; layer++) {
		imstats* stat = statistics(NULL, -1, &gfit, layer, &com.selection, STATS_MAIN, TRUE);
		if (!stat) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 1;
		}

		switch (layer) {
			case 0:
				if (gfit.naxes[2] == 1)
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

		siril_log_message(
				_("%s layer: Mean: %0.1lf, Median: %0.1lf, Sigma: %0.1lf, "
						"AvgDev: %0.1lf, Min: %0.1lf, Max: %0.1lf\n"),
				layername, stat->mean, stat->median, stat->sigma,
				stat->avgDev, stat->min, stat->max);
		free_stats(stat);
	}
	return 0;
}

int process_convertraw(int nb) {
	GDir *dir;
	GError *error = NULL;
	const gchar *file;
	GList *list = NULL;

	struct timeval t_start;

	if (get_thread_run()) {
		siril_log_message(_("Another task is "
				"already in progress, ignoring new request.\n"));
		return 1;
	}

	if (word[2]) {
		if (!strcmp(word[2], "-debayer")) {
			set_debayer_in_convflags();
		}
	}

	if((dir = g_dir_open(com.wd, 0, &error)) == NULL){
		siril_log_message(_("Conversion: error opening working directory %s.\n"), com.wd);
		fprintf (stderr, "Conversion: %s\n", error->message);
		g_error_free(error);
		set_cursor_waiting(FALSE);
		return 1;
	}

	while ((file = g_dir_read_name(dir)) != NULL) {
		const char *ext = get_filename_ext(file);
		if (!ext)
			continue;
		image_type type = get_type_for_extension(ext);
		if (type == TYPERAW) {
			list = g_list_append (list, g_strdup(file));
		}
	}
	/* sort list */
	list = g_list_sort(list, (GCompareFunc) strcompare);

	siril_log_color_message(_("Conversion: processing...\n"), "red");
	gettimeofday(&t_start, NULL);

	set_cursor_waiting(TRUE);
	if (!com.script)
		control_window_switch_to_tab(OUTPUT_LOGS);

	/* then, convert files to Siril's FITS format */
	struct _convert_data *args;
	set_cursor_waiting(TRUE);
	if (!com.wd) {
		siril_log_message(_("Conversion: no working directory set.\n"));
		set_cursor_waiting(FALSE);
		return 1;
	}

	args = malloc(sizeof(struct _convert_data));
	args->start = 1;
	args->dir = dir;
	args->list = list;
	args->total = g_list_length(list);
	args->nb_converted = 0;
	args->t_start.tv_sec = t_start.tv_sec;
	args->t_start.tv_usec = t_start.tv_usec;
	args->compatibility = FALSE;	// not used here
	args->command_line = TRUE;
	args->destroot = g_strdup(word[1]);
	start_in_new_thread(convert_thread_worker, args);
	return 0;
}

int process_register(int nb) {
	struct registration_args *reg_args;
	struct registration_method *method;
	char *msg;
	int i;

	if (get_thread_run()) {
		siril_log_message(_("Another task is "
				"already in progress, ignoring new request.\n"));
		return 1;
	}

	gchar *file = g_strdup(word[1]);
	if (!ends_with(file, ".seq")) {
		str_append(&file, ".seq");
	}

	if (!existseq(file)) {
		if (check_seq(FALSE)) {
			siril_log_message(_("No sequence `%s' found.\n"), file);
			g_free(file);
			return 1;
		}
	}
	sequence *seq = readseqfile(file);
	if (seq == NULL) {
		siril_log_message(_("No sequence `%s' found.\n"), file);
		g_free(file);
		return 1;
	}
	if (seq_check_basic_data(seq, FALSE) == -1) {
		g_free(file);
		free(seq);
		return 1;
	}

	g_free(file);
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

	/* check for options */
	for (i = 2; i < 4; i++) {
		if (word[i]) {
			if (!strcmp(word[i], "-drizzle")) {
				reg_args->x2upscale = TRUE;
			} else if (!strcmp(word[i], "-norot")) {
				reg_args->translation_only = TRUE;
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
		int64_t size = seq_compute_size(reg_args->seq, nb_frames, DATA_USHORT);
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
			_("Registration: processing using method: %s\n"), "red",
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
		else if (g_str_has_prefix(current, "-norm=")) {
			if (!norm_allowed) {
				siril_log_message(_("Normalization options are not allowed in this context, ignoring.\n"));
			} else {
				value = current+6;
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
			args.type_of_rejection = WINSORIZED;
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
		args.norm_to_16 = TRUE;
		args.reglayer = args.seq->nb_layers == 1 ? 0 : 1;

		// manage filters
		if (convert_stack_data_to_filter(arg, &args) ||
				setup_filtered_data(&args)) {
			free_sequence(seq, TRUE);
			return 1;
		}
		args.description = describe_filter(seq, args.filtering_criterion, args.filtering_parameter);
		args.use_32bit_output = evaluate_stacking_should_output_32bits(args.method,
			args.seq, args.nb_images_to_stack);

		if (!arg->result_file) {
			char filename[256];
			char *suffix = ends_with(seq->seqname, "_") ? "" :
				(ends_with(com.seq.seqname, "-") ? "" : "_");
			snprintf(filename, 256, "%s%sstacked%s",
					seq->seqname, suffix, com.ext);
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
			g_error_free(error);
		}
		siril_add_idle(end_generic, NULL);
		return NULL;
	}

	siril_log_message(_("Starting stacking of found sequences...\n"));
	while ((file = g_dir_read_name(dir)) != NULL) {
		if (ends_with(file, ".seq")) {
			arg->seqfile = strdup(file);
			stack_one_seq(arg);
		}
	}

	siril_log_message(_("Stacked %d sequences successfully.\n"), arg->number_of_loaded_sequences);
	gettimeofday(&t_end, NULL);
	show_time(arg->t_start, t_end);
	g_dir_close(dir);
	g_free(arg->result_file);
	g_free(arg->seqfile);
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

	// stackall { sum | min | max } [-filter-fwhm=value[%]] [-filter-wfwhm=value[%]] [-filter-round=value[%]] [-filter-quality=value[%]] [-filter-incl[uded]]
	// stackall { med | median } [-nonorm, norm=] [-filter-incl[uded]]
	// stackall { rej | mean } sigma_low sigma_high [-nonorm, norm=] [-filter-fwhm=value[%]] [-filter-round=value[%]] [-filter-quality=value[%]] [-filter-incl[uded]]
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
			if (!word[2] || !word[3] || (arg->sig[0] = atof(word[2])) < 0.0
					|| (arg->sig[1] = atof(word[3])) < 0.0) {
				siril_log_message(_("The average stacking with rejection uses the Winsorized "
										"rejection here and requires two extra arguments: sigma low and high.\n"));
				goto failure;
			}
			arg->method = stack_mean_with_rejection;
			start_arg_opt = 4;
			allow_norm = TRUE;
		}
		else {
			siril_log_message(_("Stacking method type '%s' is invalid\n"), word[1]);
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

	if (!retval)
		siril_log_message(_("Stacked sequence successfully.\n"));

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
	gchar *file;

	arg = calloc(1, sizeof(struct stacking_configuration));
	arg->f_fwhm = -1.f; arg->f_fwhm_p = -1.f; arg->f_round = -1.f;
	arg->f_round_p = -1.f; arg->f_quality = -1.f; arg->f_quality_p = -1.f;
	arg->filter_included = FALSE; arg->norm = NO_NORM; arg->force_no_norm = FALSE;

	file = g_strdup(word[1]);
	if (!ends_with(file, ".seq")) {
		str_append(&file, ".seq"); // reallocs file
	}

	if (!existseq(file)) {
		if (check_seq(FALSE)) {
			siril_log_message(_("No sequence `%s' found.\n"), file);
			goto failure;
		}
	}
	arg->seqfile = file;

	// stack seqfilename { sum | min | max } [-filter-fwhm=value[%]] [-filter-wfwhm=value[%]] [-filter-round=value[%]] [-filter-quality=value[%]] [-filter-incl[uded]] -out=result_filename
	// stack seqfilename { med | median } [-nonorm, norm=] [-filter-incl[uded]] -out=result_filename
	// stack seqfilename { rej | mean } sigma_low sigma_high [-nonorm, norm=] [-filter-fwhm=value[%]] [-filter-round=value[%]] [-filter-quality=value[%]] [-filter-incl[uded]] -out=result_filename
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
			if (!word[3] || !word[4] || (arg->sig[0] = atof(word[3])) < 0.0
					|| (arg->sig[1] = atof(word[4])) < 0.0) {
				siril_log_message(_("The average stacking with rejection uses the Winsorized "
										"rejection here and requires two extra arguments: sigma low and high.\n"));
				goto failure;
			}
			arg->method = stack_mean_with_rejection;
			start_arg_opt = 5;
			allow_norm = TRUE;
		}
		else {
			siril_log_message(_("Stacking method type '%s' is invalid\n"), word[2]);
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
	int nb_command_max = 10;
	gchar *file;
	int i, retvalue = 0;

	if (word[1][0] == '\0') {
		return -1;
	}

	file = g_strdup(word[1]);
	if (!ends_with(file, ".seq")) {
		str_append(&file, ".seq");
	}

	if (!existseq(file)) {
		if (check_seq(FALSE)) {
			siril_log_message(_("No sequence `%s' found.\n"), file);
			g_free(file);
			return 1;
		}
	}
	sequence *seq = readseqfile(file);
	if (seq == NULL) {
		siril_log_message(_("No sequence `%s' found.\n"), file);
		g_free(file);
		return 1;
	}
	g_free(file);
	if (seq_check_basic_data(seq, FALSE) == -1) {
		free(seq);
		return 1;
	}

	args = calloc(1, sizeof(struct preprocessing_data));

	/* checking for options */
	for (i = 2; i < nb_command_max; i++) {
		if (word[i]) {
			if (g_str_has_prefix(word[i], "-bias=")) {
				args->bias = calloc(1, sizeof(fits));
				if (!readfits(word[i] + 6, args->bias, NULL)) {
					args->use_bias = TRUE;
				} else {
					retvalue = 1;
					free(args->bias);
					break;
				}
			} else if (g_str_has_prefix(word[i], "-dark=")) {
				args->dark = calloc(1, sizeof(fits));
				if (!readfits(word[i] + 6, args->dark, NULL)) {
					args->use_dark = TRUE;
					args->use_cosmetic_correction = TRUE;
				} else {
					retvalue = 1;
					free(args->dark);
					break;
				}
			} else if (g_str_has_prefix(word[i], "-flat=")) {
				args->flat = calloc(1, sizeof(fits));
				if (!readfits(word[i] + 6, args->flat, NULL)) {
					args->use_flat = TRUE;
				} else {
					retvalue = 1;
					free(args->flat);
					break;
				}
			} else if (!strcmp(word[i], "-opt")) {
				args->use_dark_optim = TRUE;
			} else if (!strcmp(word[i], "-cfa")) {
				args->is_cfa = TRUE;
			} else if (!strcmp(word[i], "-debayer")) {
				args->debayer = TRUE;
			} else if (!strcmp(word[i], "-stretch")) {
				args->stretch_cfa = TRUE;
			} else if (!strcmp(word[i], "-flip")) {
				args->compatibility = TRUE;
			} else if (!strcmp(word[i], "-equalize_cfa")) {
				args->equalize_cfa = TRUE;
			}
		}
	}

	if (retvalue) {
		free(args);
		return -1;
	}

	siril_log_color_message(_("Preprocessing...\n"), "red");
	gettimeofday(&args->t_start, NULL);
	args->seq = seq;
	args->is_sequence = TRUE;
	args->autolevel = TRUE;
	args->normalisation = 1.0f;	// will be updated anyway
	args->sigma[0] = -1.00; /* cold pixels: it is better to deactive it */
	args->sigma[1] =  3.00; /* hot pixels */
	args->ppprefix = "pp_";
	args->allow_32bit_output = args->seq->type == SEQ_REGULAR;

	// start preprocessing
	set_cursor_waiting(TRUE);

	start_sequence_preprocessing(args, TRUE);
	return 0;
}


#ifdef _OPENMP
int process_set_cpu(int nb){
	int proc_in, proc_out, proc_max;

	proc_in = atoi(word[1]);
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
	siril_log_message(_("Using now %d logical processors\n"), proc_out);
	com.max_thread = proc_out;
	if (!com.headless)
		update_spinCPU(0);

	return 0;
}
#endif

int process_set_mem(int nb){
	double ratio = atof(word[1]);
	if (ratio < 0.05 || ratio > 4.0) {
		siril_log_message(_("The accepted range for the ratio of memory used for stacking is [0.05, 4], with values below the available memory recommended\n"));
		return 1;
	}
	if (ratio > 1.0) {
		siril_log_message(_("Setting the ratio of memory used for stacking above 1 will require the use of on-disk memory, which can be very slow and is unrecommended (%g requested)\n"), ratio);
	}
	com.stack.memory_ratio = ratio;
	if (!writeinitfile())
		siril_log_message(_("Usable memory for stacking changed to %g\n"), ratio);
	return 0;
}

int process_help(int nb){
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

	Nbr_Plan = atoi(word[1]);

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
