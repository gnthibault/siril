/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2016 team free-astro (see more in AUTHORS file)
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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <dirent.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_histogram.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "core/siril.h"
#include "core/command.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/processing.h"
#include "io/conversion.h"
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "gui/PSF_list.h"
#include "gui/histogram.h"
#include "gui/plot.h"
#include "algos/colors.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/Def_Math.h"
#include "algos/Def_Wavelet.h"
#include "algos/gradient.h"
#include "algos/fft.h"
#include "algos/quality.h"
#include "algos/cosmetic_correction.h"
#include "stacking/stacking.h"
#include "registration/registration.h"

#ifdef HAVE_OPENCV
#include "opencv/opencv.h"
#endif

static char *word[MAX_COMMAND_WORDS];

command commande[] = {
	/* name,	nbarg,	usage,			function pointer */
	{"addmax",	1,	"addmax filename",	process_addmax},
	
	{"bg", 0, "bg", process_bg},
	{"bgnoise", 0, "bgnoise", process_bgnoise},
	
	{"cd", 1, "cd directory (define the working directory)", process_cd},
	{"cdg", 0, "cdg", process_cdg},
	{"clearstar", 0, "clearstar", process_clearstar},
	{"contrast", 0, "contrast", process_contrast},
	{"cosme", 1, "cosme [filename].lst", process_cosme},
	{"cosme_cfa", 1, "cosme [filename].lst", process_cosme},
	{"crop", 0, "crop [x y width height]", process_crop}, 

	{"ddp", 3, "ddp level coef sigma", process_ddp}, 
	
	{"entropy", 0, "entropy", process_entropy},
	{"exit", 0, "exit", process_exit},
	{"extract", 1, "extract NbPlans", process_extract},
	
	{"fdiv", 2, "fdiv filename scalar", process_fdiv},
	{"fftd", 2, "fftd magnitude phase", process_fft},
	{"ffti", 2, "ffti magnitude phase", process_fft},
	{"fill", 1, "fill value", process_fill},
	{"fill2", 1, "fill2 value [x y width height]", process_fill2},
	{"find_hot", 3, "find_hot filename cold_sigma hot_sigma", process_findhot},
	{"find_cosme", 2, "find_cosme cold_sigma hot_sigma", process_findcosme},
	{"find_cosme_cfa", 2, "find_cosme_cfa cold_sigma hot_sigma", process_findcosme},
	{"findstar", 0, "findstar", process_findstar},
	{"fmedian", 2, "fmedian ksize modulation", process_fmedian},
	{"fmul", 1, "fmul scalar", process_fmul},
	{"fixbanding", 2, "fixbanding amount sigma", process_fixbanding},

	
	{"gauss", 1, "gauss sigma ", process_gauss},	
	//~ {"gauss2", 1, "gauss sigma", process_gauss2},

	{"help", 0, "help", process_help},	
	{"histo", 1, "histo layer (layer=0, 1, 2 with 0: red, 1: green, 2: blue)", process_histo},
	
	/* commands oper filename and curent image */
	{"iadd", 1, "add filename", process_imoper}, 
	{"idiv", 1, "idiv filename", process_imoper},
	{"imul", 1, "imul filename", process_imoper}, 
	{"isub", 1, "isub filename", process_imoper},
	
	{"load", 1, "load filename.[ext]", process_load}, 
	// specific loads are not required, but could be used to force the
	// extension to a higher priority in case two files with same basename
	// exist (stat_file() manages that priority order for now).
	{"log", 0, "log ", process_log}, /* logarifies current image */
	{"ls", 0, "ls ", process_ls},
	
	{"mirrorx", 0, "mirrorx", process_mirrorx},
	{"mirrory", 0, "mirrory", process_mirrory},
	
	{"new", 3, "new width height nb_layers", process_new},
	{"nozero", 1, "nozero level (replaces null values by level)", process_nozero}, /* replaces null values by level */
	
	{"offset", 1, "offset value", process_offset},
	
	{"psf", 0, "psf", process_psf},
	
#ifdef HAVE_OPENCV
	{"resample", 1, "resample factor", process_resample},
#endif	
	{"rmgreen", 1, "rmgreen type", process_scnr},
#ifdef HAVE_OPENCV
	{"rotate", 1, "rotate angle", process_rotate},
#endif
	{"rotatePi", 0, "rotatePi", process_rotatepi},
	
	{"satu", 1, "satu coeff ", process_satu}, 
	{"save", 1, "save filename (save current image in fit)", process_save}, 
	{"savebmp", 1, "savebmp filename (save display image in bmp)", process_savebmp}, 
#ifdef HAVE_LIBJPEG
	{"savejpg", 1, "savejpg filename [quality] (save current display in jpg)", process_savejpg},
#endif
	{"savepnm", 1, "savepnm filename (save current image in Netpbm)", process_savepnm},
#ifdef HAVE_LIBTIFF
	{"savetif", 1, "savetif filename (save current image in tif 16bits)", process_savetif},
	{"savetif8", 1, "savetif8 filename (save current image in tif 8bits)", process_savetif},
#endif
	{"select", 2, "select from to", process_select},
	{"seqcrop", 0, "seqcrop", process_seq_crop},
	{"seqfind_cosme", 2, "seqfind_cosme cold_sigma hot_sigma", process_findcosme},
	{"seqfind_cosme_cfa", 2, "seqfind_cosme_cfa cold_sigma hot_sigma", process_findcosme},
	{"seqpsf", 0, "seqpsf", process_seq_psf},
#ifdef _OPENMP
	{"setcpu", 1, "setcpu number", process_set_cpu},
#endif
	{"setmag", 1, "setmag magnitude", process_set_mag},
	{"setmagseq", 1, "setmagseq magnitude", process_set_mag_seq},
	{"split", 3, "split R G B", process_split},
	{"stat", 0, "stat", process_stat},
	{"stackall", 0, "stackall", process_stackall},
	
	{"threshlo", 1, "threshlo level", process_threshlo},
	{"threshhi", 1, "threshi level", process_threshhi}, 
	{"thresh", 2, "thresh hi lo (threshes hi and lo)", process_thresh}, /* threshes hi and lo */
	
	/* unsharp masking of current image or genname sequence */
	{"unselect", 2, "unselect from to", process_unselect},
	{"unsharp", 2, "unsharp sigma multi", process_unsharp},
	{"unsetmag", 0, "unsetmag", process_unset_mag},
	{"unsetmagseq", 0, "unsetmagseq", process_unset_mag_seq},
//	{"unsharp2", 5, "unsharp2 sigma multi src dest number", process_unsharp2},

	{"visu", 2, "visu low high", process_visu},
	
	/* wavelet transform in nbr_plan plans */ 
	{"wavelet", 1, "wavelet nbr_plan type (1=linear 2=spline)", process_wavelet},
	/* reconstruct from wavelet transform and weighs plans with c1, c2, c3... */ 
	{"wrecons", 2, "wrecons c1 c2 c3 ...", process_wrecons},
	
	{"",0,"",0}
};

int process_load(int nb){
	char filename[256];
	int retval, i;
	
	strncpy(filename, word[1], 250);
	filename[250] = '\0';	
	
	for (i=1; i<nb-1; ++i){
		strcat(filename, " ");
		strcat(filename, word[i+1]);
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
	struct enhance_saturation_data *args = malloc(sizeof(struct enhance_saturation_data));
	
	args->coeff = atof(word[1]);
	if (args->coeff == 0.0) args->coeff = 1.0;

	args->fit = &gfit;
	args->h_min = 0.0;
	args->h_max = 360.0;
	args->preserve = TRUE;
	set_cursor_waiting(TRUE);
	start_in_new_thread(enhance_saturation, args);
	
	return 0;
}

int process_save(int nb){
	char filename[256];
	
	if (sequence_is_loaded() && !single_image_is_loaded()) {
		gfit.hi = com.seq.layers[RLAYER].hi;
		gfit.lo = com.seq.layers[RLAYER].lo;
	}
	else {
		gfit.hi = com.uniq->layers[RLAYER].hi;
		gfit.lo = com.uniq->layers[RLAYER].lo;
	}


	sprintf(filename, "%s", word[1]);
	set_cursor_waiting(TRUE);
	savefits(filename, &(gfit));
	set_cursor_waiting(FALSE);
	return 0;
}

int process_savebmp(int nb){
	char filename[256];
	
	sprintf(filename, "%s", strcat(word[1], ".bmp"));
	set_cursor_waiting(TRUE);
	savebmp(filename, &(gfit));
	set_cursor_waiting(FALSE);
	return 0;
}

#ifdef HAVE_LIBJPEG
int process_savejpg(int nb){
	char filename[256];
	int quality = 100;
	
	if ((nb == 3) && atoi(word[2]) <= 100 && atoi(word[2]) > 0)
		quality=atoi(word[2]);
	strcpy(filename, word[1]);
	strcat(filename, ".jpg");
	set_cursor_waiting(TRUE);
	savejpg(filename, &gfit, quality);
	set_cursor_waiting(FALSE);
	return 0;
}
#endif

#ifdef HAVE_LIBTIFF
int process_savetif(int nb){
	char filename[256];
	uint16 bitspersample = 16;
	
	if (strcasecmp(word[0],"savetif8")==0) bitspersample=8;
	sprintf(filename,"%s", strcat(word[1],".tif"));
	set_cursor_waiting(TRUE);
	savetif(filename, &gfit, bitspersample);
	set_cursor_waiting(FALSE);
	return 0;
}
#endif

int process_savepnm(int nb){
	char filename[256];
	
	switch(gfit.naxes[2]){
		case 1:
			sprintf(filename,"%s", strcat(word[1],".pgm"));
				set_cursor_waiting(TRUE);
				savepgm(filename, &(gfit));
				set_cursor_waiting(FALSE);
		break;
		case 3:
			sprintf(filename,"%s", strcat(word[1],".ppm"));
			set_cursor_waiting(TRUE);
			saveppm(filename, &(gfit));
			set_cursor_waiting(FALSE);
		break;
		
		default:	/*Should not happen */
			return 1;
		}
	return 0;	
}

int process_imoper(int nb){
	clearfits(&(wfit[4]));
	readfits(word[1], &(wfit[4]), NULL);
	imoper(&gfit, &wfit[4], word[0][1]);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_addmax(int nb){
	clearfits(&(wfit[4]));
	readfits(word[1], &(wfit[4]), NULL);
	if (addmax(&gfit, &wfit[4])==0) {
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	}
	return 0;
}

int process_fdiv(int nb){
	// combines an image division and a scalar multiplication.
	float norm;

	norm = atof(word[2]);
	clearfits(&(wfit[4]));
	readfits(word[1], &(wfit[4]), NULL);
	fdiv(&gfit,&wfit[4], norm);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_fmul(int nb){
	float coeff;

	coeff = atof(word[1]);
	if (coeff <= 0.0) {
		siril_log_message(_("Multiplying by a coefficient less than or equal to 0 is not possible.\n"));
		return 1;
	}
	soper(&gfit, coeff, OPER_MUL);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_entropy(int nb){
	rectangle area;
	double e;

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
	unsharp(&(gfit), atof(word[1]), (double)0, TRUE);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_unsharp(int nb){
	unsharp(&(gfit), atof(word[1]), atof(word[2]), TRUE);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_crop(int nb){
	rectangle area;
	if ((!com.selection.h) || (!com.selection.w)) {
		if (nb==5){
			if (atoi(word[1])<0 || atoi(word[2])<0){
				siril_log_message(_("Crop: x and y must be positive values.\n"));
				return 1;
			}			
			if (atoi(word[3])<=0 || atoi(word[4])<=0){
				siril_log_message(_("Crop: width and height must be greater than 0.\n"));
				return 1;
			}
			if (atoi(word[3])>gfit.rx || atoi(word[4])>gfit.ry){
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
	update_used_memory();
	return 0;
}

int process_cd(int nb){
	char filename[256];
	int retval, i;
	
	strncpy(filename, word[1], 250);
	filename[250] = '\0';
	
	for (i=1; i<nb-1; ++i){
		strcat(filename, " ");
		strcat(filename, word[i+1]);
	}
	
	expand_home_in_filename(filename, 256);
	if (!(retval=changedir(filename)))
		writeinitfile();
	return retval;
}

int process_wrecons(int nb){
	int i;
	float coef[7];
	char *File_Name_Transform[3] = {"r_rawdata.wave", "g_rawdata.wave", "b_rawdata.wave"}, *dir[3];
	const char *tmpdir;
	int nb_chan = gfit.naxes[2];
	
	assert(nb_chan == 1 || nb_chan == 3);
		
 	tmpdir = g_get_tmp_dir();
 	
	for (i=0;i<nb-1;++i){
		coef[i] = atof(word[i+1]);
	}

	for (i=0; i < nb_chan; i++) {
		dir[i] = malloc(strlen(tmpdir) + strlen(File_Name_Transform[i]) + 2);
		strcpy(dir[i], tmpdir);
		strcat(dir[i], "/");
		strcat(dir[i], File_Name_Transform[i]);
		wavelet_reconstruct_file (dir[i], coef, gfit.pdata[i]);
		free(dir[i]);
	}

	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}	

int process_wavelet(int nb){
	char *File_Name_Transform[3] = {"r_rawdata.wave", "g_rawdata.wave", "b_rawdata.wave"}, *dir[3];
	const char* tmpdir;
	int Type_Transform, Nbr_Plan, maxplan, mins, chan, nb_chan;
	float *Imag;
	
 	tmpdir = g_get_tmp_dir();
	
	Nbr_Plan = atoi(word[1]);
	Type_Transform = atoi(word[2]);
	
	nb_chan = gfit.naxes[2];
	assert(nb_chan <= 3);

	mins = min (gfit.rx, gfit.ry);
	maxplan = log(mins) / log(2) - 2;

	if ( Nbr_Plan > maxplan ){
		siril_log_message(_("Wavelet: maximum number of plans for this image size is %d\n"),
				maxplan);
		return 1;
	}

	if(Type_Transform != TO_PAVE_LINEAR && Type_Transform !=TO_PAVE_BSPLINE){
		siril_log_message(_("Wavelet: type must be %d or %d\n"), TO_PAVE_LINEAR, TO_PAVE_BSPLINE);
		return 1;
	}

	Imag = f_vector_alloc (gfit.rx * gfit.ry);
	
	for (chan = 0; chan < nb_chan; chan++) {
		dir[chan] = malloc(strlen(tmpdir) + strlen(File_Name_Transform[chan]) + 2);
		strcpy(dir[chan], tmpdir);
		strcat(dir[chan], "/");
		strcat(dir[chan], File_Name_Transform[chan]);
		wavelet_transform_file (Imag, gfit.ry, gfit.rx, dir[chan], Type_Transform, Nbr_Plan, gfit.pdata[chan]);
		free(dir[chan]);
	}
	
	free (Imag);
	return 0;
}

int process_log(int nb){
	loglut(&gfit, LOG);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_ls(int nb){
	struct dirent **list;
	char filename[256], *path;
	int i, n;
	
	filename[0]='\0';
	/* If a path is given in argument */
	if (nb>1){
		if (word[1][0]!='\0'){	
			/* Absolute path */
			if(word[1][0]=='/' || word[1][0]=='~'){
				strncpy(filename, word[1], 250);
				filename[250] = '\0';
				expand_home_in_filename(filename, 256);
			}
			/* Relative path */
			else {
				strcpy(filename, com.wd);
				strcat(filename, "/");
				strcat(filename, word[1]);
			}	
			path = strdup(filename);
		}
		/* Should not happen */
		else {
			printf("Cannot list files in %s\n", filename);
			return 1;
		}
	}
	/* No paths are given in argument */
	else {
		if (!com.wd) {
			siril_log_message(_("Cannot list files, set working directory first.\n"));
			return 1;
		}
		path = strdup(com.wd);
	}
	if (!path) {
		siril_log_message(_("Siril cannot open the directory.\n"));
		return 1;
	}

	n = scandir(path, &list, 0, alphasort);
	if (n < 0)
		perror("scandir");

	/* List the entries */
	for (i = 0; i < n; ++i) {
		struct stat entrystat;
		char file_path[256];
		const char *ext;
		if (list[i]->d_name[0] == '.')
			continue; /* no hidden files */
		if (filename[0] != '\0')
			sprintf(file_path, "%s/%s", filename, list[i]->d_name);
		else {
			sprintf(file_path, "%s", list[i]->d_name);
		}
		if (lstat(file_path, &entrystat)) {
			perror("stat");
			break;
		}
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
		} else if (!strncmp(ext, "seq", 4))
			siril_log_color_message(_("Sequence: %s\n"), "blue", list[i]->d_name);
	}
	siril_log_message(_("********* END OF THE LIST *********\n"));
	for (i = 0; i < n; i++)
		free(list[i]);
	free(list);
	free(path);

	return 0;
}

int	process_mirrorx(int nb){
	mirrorx(&gfit, TRUE);
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int	process_mirrory(int nb){
	mirrory(&gfit, TRUE);
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

#ifdef HAVE_OPENCV
int process_resample(int nb) {
	double factor = atof(word[1]);
	if (factor > 5.0) {
		siril_log_message(_("The scaling factor must be less than 5.0\n"));
		return 1;
	}
	int toX = round_to_int(factor * gfit.rx);
	int toY = round_to_int(factor * gfit.ry);
	
	set_cursor_waiting(TRUE);
	verbose_resize_gaussian(&gfit, toX, toY, OPENCV_LINEAR);
	update_used_memory();
	adjust_vport_size_to_image();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();	
	set_cursor_waiting(FALSE);
	return 0;
}


int process_rotate(int nb) {
	double degree;
	
	set_cursor_waiting(TRUE);
	degree = atof(word[1]);
	verbose_rotate_image(&gfit, degree, OPENCV_LINEAR, 1);	//INTER_LINEAR
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
	return 0;
}
#endif

int process_rotatepi(int nb){
#ifdef HAVE_OPENCV
	verbose_rotate_image(&gfit, 180.0, OPENCV_LINEAR, 1);
#else
	fits_rotate_pi(&gfit);
#endif
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;	
}

int process_set_mag(int nb) {
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
		fitted_PSF *result = psf_get_minimisation(&gfit, layer, &com.selection);
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
		fitted_PSF *result = psf_get_minimisation(&gfit, layer, &com.selection);
		if (result) {
			psf_display_result(result, &com.selection);
			free(result);
		}
	}
	return 0;
}

gboolean end_seqpsf(gpointer arg) {
	set_layers_for_registration();	// update display of available reg data
	if (com.seq.needs_saving)
		writeseqfile(&com.seq);
	drawPlot();
	notify_new_photometry();
	return end_generic(arg);
}

void *_psf_thread(void *arg) {
	int layer = (intptr_t) arg;
	do_fwhm_sequence_processing(&com.seq, layer, TRUE, TRUE, TRUE, FALSE);
	gdk_threads_add_idle(end_seqpsf, NULL);
	return NULL;
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
		siril_log_message(_("Running the PSF on the loaded sequence, layer %d\n"), layer);
		siril_log_message(_("Results will be displayed at the end of the processing, on the console output, in the following form:\n"));
		start_in_new_thread(_psf_thread, (void *)(intptr_t)layer);
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

	if (com.selection.w != 0 || com.selection.h != 0)
		return 1;

	struct crop_sequence_data *args = malloc(sizeof(struct crop_sequence_data));

	args->seq = &com.seq;
	args->area = &com.selection;
	args->prefix = "cropped_";

	set_cursor_waiting(TRUE);
	start_in_new_thread(crop_sequence, args);
	return 0;
}

int process_bg(int nb){
	WORD bg = round_to_WORD(background(&gfit, -1, &com.selection));
	siril_log_message(_("Background value: %d\n"), bg);
	return 0;
}

int process_bgnoise(int nb){
	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return 1;
	}

	struct noise_data *args = malloc(sizeof(struct noise_data));

	set_cursor_waiting(TRUE);
	control_window_switch_to_tab(OUTPUT_LOGS);

	args->fit = &gfit;
	args->verbose = TRUE;
	memset(args->bgnoise, 0.0, sizeof(double[3]));
	start_in_new_thread(noise, args);
	return 0;
}

int process_histo(int nb){
	size_t i;
	int nlayer = atoi(word[1]);
	char* clayer;
	char name [20];
	
	if (nlayer>3 || nlayer <0)
		return 1;
	gsl_histogram* histo = computeHisto(&gfit, nlayer);
	if (!isrgb(&gfit)) clayer = strdup("bw");		//if B&W
	else clayer = vport_number_to_name(nlayer);
	snprintf(name, 20, "histo_%s.dat",clayer);

	FILE *f = fopen(name, "w");

	if (f == NULL) {
		free(clayer);
		return 1;
	}
	for (i=0; i < USHRT_MAX + 1; i++)
		fprintf(f, "%zu %d\n", i, (int) gsl_histogram_get (histo, i));
	fclose(f);
	gsl_histogram_free(histo);
	siril_log_message(_("The file %s has been created for the %s layer.\n"), name, clayer);
	return 0;
}

int process_thresh(int nb){
	int lo, hi;

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

	lo = atoi(word[1]);
	threshlo(&gfit, lo);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_threshhi(int nb){
	int hi;

	hi = atoi(word[1]);
	threshhi(&gfit, hi);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_nozero(int nb){
	int level;

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

	new_fit_image(&gfit, width, height, layers);

	open_single_image_from_gfit(strdup(_("new empty image")));
	return 0;
}

int process_visu(int nb){
	int low, high;
	
	low = atoi(word[1]);
	high = atoi(word[2]);
	if ((high>USHRT_MAX) || (low<0)){
		siril_log_message(_("Values must be positive and less than %d.\n"), USHRT_MAX);
		return 1;		
	}
	visu(&gfit, low, high);
	return 0;
}

int process_fill2(int nb){
	int level=atoi(word[1]);
	rectangle area;
	if (!com.drawn || com.drawing){		// TODO: what's that test?
		if (nb==6){
			area.x = atoi(word[2]);
			area.y = atoi(word[3]);
			area.w = atoi(word[4]);
			area.h = atoi(word[5]);
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
	int layer = RLAYER;
	starFinder sf;

	memset(&sf, 0, sizeof(starFinder));
	
	if (!single_image_is_loaded()) return 0;
	if (isrgb(&gfit)) layer = GLAYER;
	delete_selected_area();
	com.stars = peaker(&gfit, layer, &sf, NULL);
	refresh_stars_list(com.stars);
	return 0;
}

int process_findhot(int nb){
	long icold, ihot;
	char filename[256];
	int i;
	char type;
	if (gfit.naxes[2] != 1) {
		siril_log_message(_("find_hot must be applied on an one-channel master-dark frame"));
		return 1;
	}
	double sig[2];
	sig[0] = atof(word[2]);
	sig[1] = atof(word[3]);

	deviant_pixel *dev = find_deviant_pixels(&gfit, sig, &icold, &ihot);
	siril_log_message(_("%ld cold and %ld hot pixels\n"), icold, ihot);

	FILE* cosme_file = NULL;
	sprintf(filename, "%s.lst", word[1]);
	cosme_file = fopen(filename, "w");
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
	double dirty;
	int is_cfa, i = 0, retval = 0;
	int nb_tokens;
	char line[64];
	char type;

	if (!ends_with(word[1], ".lst"))
		strcat(word[1], ".lst");
	cosme_file = fopen(word[1], "r");
	if (cosme_file == NULL) {
		siril_log_message(_("Cannot open file: %s\n"), word[1]);
		return 1;
	}
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
#ifdef HAVE_OPENCV
			nb_tokens = sscanf(line + 2, "%lf %lf %c", &dev.p.y, &dirty, &type);
			if (nb_tokens != 2 && nb_tokens != 3) {
				fprintf(stderr, "cosmetic correction: "
						"cosme file format error at line %d: %s\n", i, line);
				retval = 1;
				continue;
			}
			dev.type = HOT_PIXEL; // we force it
			dev.p.y = gfit.rx - dev.p.y - 1; /* FITS are stored bottom to top */
			cvRotateImage(&gfit, 90.0, -1, OPENCV_LINEAR);
			cosmeticCorrOneLine(&gfit, dev, is_cfa);
			cvRotateImage(&gfit, -90.0, -1, OPENCV_LINEAR);
#else
			siril_log_message(_("Opencv need to be compiled to remove bad column.\n"));
			retval = 1;
#endif
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

	FindCentre(&gfit, &x_avg, &y_avg);
	y_avg = gfit.ry - y_avg;	// FITS are stored bottom to top
	siril_log_message(_("Center of gravity coordinates are (%.3lf, %.3lf)\n"), x_avg, y_avg);

	return 0;
}

int process_clearstar(int nb){
	clear_stars_list();
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_NONE);
	redraw_previews();
	return 0;
}

int process_contrast(int nb){
	int layer;
	double result[gfit.naxes[2]], value=0;

	for (layer = 0; layer < gfit.naxes[2]; layer++)
		result[layer] = contrast(&gfit, layer);
	for (layer = 0; layer < gfit.naxes[2]; layer++)
		value += result[layer];
	value /= gfit.naxes[2];
	
	siril_log_message(_("Contrast: %lf\n"), value);
	return 0;
}

int process_fill(int nb){	
	int level;
	rectangle area;
	
	if (!com.drawn || com.drawing){		// TODO: what's that test?
		if (nb==6){
			area.x = atoi(word[2]);
			area.y = atoi(word[3]);
			area.w = atoi(word[4]);
			area.h = atoi(word[5]);
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
	
	level = atoi(word[1]);
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
	if (sequence_is_loaded()) {
		siril_log_message(_("FFT does not work with sequences\n"));
		return 1;
	}
	if (get_thread_run()) {
		siril_log_message(_("Another task is already in progress, ignoring new request.\n"));
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
		siril_log_message(_("Another task is already in progress, ignoring new request.\n"));
		return 1;
	}
	
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

	if (!sequence_is_loaded() && !single_image_is_loaded())
		return 1;
	struct cosmetic_data *args = malloc(sizeof(struct cosmetic_data));

	args->sigma[0] = atof(word[1]);
	args->sigma[1] = atof(word[2]);
	if (word[0][10] == '_' || word[0][13] == '_') {	// find_cosme_cfa or seqfind_cosme_cfa
		args->is_cfa = TRUE;
	}
	else {
		args->is_cfa = FALSE;
	}
	args->fit = &gfit;
	set_cursor_waiting(TRUE);
	if (word[0][0] == 's' && sequence_is_loaded()) {
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
			sequence_list_change_selection_index(i);
			if (select)
				com.seq.selnum++;
			else	com.seq.selnum--;
			if (i == com.seq.current)
				current_updated = TRUE;
		}
	}

	if (current_updated) {
		adjust_exclude(com.seq.current, TRUE);
	}

	update_reg_interface(FALSE);
	adjust_sellabel();
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
	char R[256], G[256], B[256];
	
	if (!isrgb(&gfit)) {
		siril_log_message(_("Siril cannot split layers. Make sure your image is in RGB mode.\n"));
		return 1;
	}
	sprintf(R, "%s%s", word[1], com.ext);
	sprintf(G, "%s%s", word[2], com.ext);
	sprintf(B, "%s%s", word[3], com.ext);
	save1fits16(R, &gfit, RLAYER);
	save1fits16(G, &gfit, GLAYER);
	save1fits16(B, &gfit, BLAYER);
	return 0;
}

int process_stat(int nb){
	int nplane = gfit.naxes[2];
	int layer;

	for (layer = 0; layer < nplane; layer++) {
		imstats* stat = statistics(&gfit, layer, &com.selection, STATS_MAIN,
		STATS_ZERO_NULLCHECK);
		if (!stat) {
			siril_log_message(_("Error: no data computed.\n"));
			return 1;
		}
		siril_log_message(
				_("%s layer: Mean: %0.1lf, Median: %0.1lf, Sigma: %0.1lf, "
						"AvgDev: %0.1lf, Min: %0.1lf, Max: %0.1lf\n"),
				stat->layername, stat->mean, stat->median, stat->sigma,
				stat->avgDev, stat->min, stat->max);
		free(stat);
		stat = NULL;
	}
	return 0;
}

gpointer stackall_worker(gpointer args) {
	DIR *dir;
	struct dirent *file;
	int number_of_loaded_sequences = 0, retval = 0;

	siril_log_message(_("Looking for sequences in current working directory...\n"));
	if (check_seq(0) || (dir = opendir(com.wd)) == NULL) {
		siril_log_message(_("Error while searching sequences or opening the directory.\n"));
		com.wd[0] = '\0';
		gdk_threads_add_idle(end_generic, NULL);
		return NULL;
	}
	siril_log_message(_("Starting stacking of found sequences...\n"));
	while ((file = readdir(dir)) != NULL) {
		char *suf;

		if ((suf = strstr(file->d_name, ".seq")) && strlen(suf) == 4) {
			sequence *seq = readseqfile(file->d_name);
			if (seq != NULL) {
				char filename[256];
				struct stacking_args args;
				//args.method = stack_summing;
				args.seq = seq;
				args.filtering_criterion = stack_filter_all;
				args.nb_images_to_stack = seq->number;
				snprintf(filename, 256, "%s%sstacked%s", seq->seqname,
						ends_with(seq->seqname, "_") ?
								"" : (ends_with(com.seq.seqname, "-") ? "" : "_"),
						com.ext);
				gettimeofday(&args.t_start, NULL);

				retval = stack_summing(&args);
				if (savefits(filename, &gfit))
					siril_log_message(_("Could not save the stacking result %s\n"),
							filename);

				free_sequence(seq, TRUE);
				++number_of_loaded_sequences;
				if (retval) break;
			}
		}
	}
	closedir(dir);
	siril_log_message(_("Stacked %d sequences %s.\n"), number_of_loaded_sequences,
			retval ? _("with errors") : _("successfully"));
	gdk_threads_add_idle(end_generic, NULL);
	return NULL;
}

int process_stackall(int nb) {
	start_in_new_thread(stackall_worker, NULL);
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

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		proc_out = omp_get_num_threads();
	}
	siril_log_message(_("Using now %d logical processors\n"), proc_out);
	com.max_thread = proc_out;
	update_spinCPU(0);

	return 0;
}
#endif

int process_help(int nb){
	command *current = commande;
	siril_log_message(_("********* LIST OF AVAILABLE COMMANDS *********\n"));
	while(current->process){
		siril_log_message("%s\n", current->usage);
		current++;
	}
	siril_log_message(_("********* END OF THE LIST *********\n"));
	return 0;
}

int process_exit(int nb){
	exit(0);
}

int process_extract(int nb) {
	int Nbr_Plan, maxplan, mins, i;
	
	Nbr_Plan = atoi(word[1]);

	mins = min (gfit.rx, gfit.ry);
	maxplan = log(mins) / log(2) - 2;

	if ( Nbr_Plan > maxplan ){
		siril_log_message(_("Wavelet: maximum number of plans for this image size is %d\n"),
				maxplan);
		return 1;
	}
	fits *fit = calloc(1, sizeof(fits));
	copyfits(&gfit, fit, CP_ALLOC | CP_COPYA | CP_FORMAT, 0);
	
	for (i=0; i < Nbr_Plan; i++) {
		char filename[256];
		
		sprintf(filename, "layer%02d", i);
		get_wavelet_layers(fit, Nbr_Plan, i, TO_PAVE_BSPLINE, -1);
		savefits(filename, fit);
	}
	clearfits(fit);
	update_used_memory();
	return 0;
}

static void parseLine(char *myline, int len, int *nb) {
	int i = 0, wordnb = 0;
	char string_starter = '\0';	// quotes don't split words on spaces
	word[0] = NULL;

	do {
		while (i < len && isblank(myline[i]))
			i++;
		if (myline[i] == '"' || myline[i] == '\'')
			string_starter = myline[i++];
		if (myline[i] == '\0' || myline[i] == '\n')
			break;
		word[wordnb++] = myline + i;	// the beginning of the word
		word[wordnb] = NULL;		// put next word to NULL
		do {
			i++;
			if (string_starter != '\0' && myline[i] == string_starter) {
				string_starter = '\0';
				break;
			}
		} while (i < len && (!isblank(myline[i]) || string_starter != '\0')
				&& myline[i] != '\n');
		if (myline[i] == '\0')	// the end of the word and line (i == len)
			break;
		myline[i++] = '\0';		// the end of the word
	} while (wordnb < MAX_COMMAND_WORDS - 1);
	*nb = wordnb;
}

static int executeCommand(int wordnb) {
	int i;
	// search for the command in the list
	if (word[0] == NULL) return 1;
	i = sizeof(commande)/sizeof(command);
	while (strcasecmp (commande[--i].name, word[0])) {
		if (i == 0) {
			siril_log_message(_("Unknown command: '%s' or not implemented yet\n"), word[0]);
			return 1 ;
		}
	}

	// verify argument count
	if(wordnb - 1 < commande[i].nbarg) {
		siril_log_message(_("Usage: %s\n"), commande[i].usage);
		return 1;
	}

	// process the command
	commande[i].process(wordnb);
	return 0;
}

int processcommand(const char *line) {
	int wordnb = 0, len, i = 0;
	char *myline;

	if (line[0] == '\0' || line[0] == '\n')
		return 0;
	if (line[0] == '@') { // case of files
		FILE * fp;
		char * linef = NULL;
		size_t lenf = 0;
		ssize_t read;

		fp = fopen(line + 1, "r");
		if (fp == NULL) {
			siril_log_message(_("File [%s] does not exist\n"), line + 1);
			return 1;
		}
		while ((read = getline(&linef, &lenf, fp)) != -1) {
			++i;
			if (linef[0] == '#') continue;	// comments
			if (linef[0] == '\0' || linef[0] == '\n')
				continue;
			myline = strdup(linef);
			parseLine(myline, read, &wordnb);
			if (executeCommand(wordnb)) {
				siril_log_message(_("Error in line: %d. Exiting batch processing\n"), i);
				free(myline);
				return 1;
			}
			free(myline);
		}

		fclose(fp);
		free(linef);
	} else {
		myline = strdup(line);
		len = strlen(line);
		parseLine(myline, len, &wordnb);
		if (executeCommand(wordnb)) {
			return 1;
		}
		free(myline);
	}
	return 0;
}
