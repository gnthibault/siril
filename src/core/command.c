/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2018 team free-astro (see more in AUTHORS file)
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

#ifdef MAC_INTEGRATION
#include "gtkmacintegration/gtkosxapplication.h"
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
#include "core/undo.h"
#include "core/initfile.h"
#include "core/processing.h"
#include "io/conversion.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "gui/PSF_list.h"
#include "gui/histogram.h"
#include "gui/plot.h"
#include "gui/progress_and_log.h"
#include "algos/colors.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/Def_Math.h"
#include "algos/Def_Wavelet.h"
#include "algos/gradient.h"
#include "algos/fft.h"
#include "algos/quality.h"
#include "algos/cosmetic_correction.h"
#include "algos/statistics.h"
#include "stacking/stacking.h"
#include "stacking/sum.h"
#include "registration/registration.h"
#include "opencv/opencv.h"

static char *word[MAX_COMMAND_WORDS];	// NULL terminated

static GThread *script_thread = NULL;

command commande[] = {
	/* name,	nbarg,	usage,			function pointer */
	{"addmax",	1,	"addmax filename",	process_addmax, N_("Computes a new image IMG with IMG_1 and IMG_2. The pixel of IMG_1 is replaced by the pixel at the same coordinates of IMG_2 if the intensity of 2 is greater than 1")},
	
	{"bg", 0, "bg", process_bg, N_("Returns the background level of the image loaded in memory")},
	{"bgnoise", 0, "bgnoise", process_bgnoise, N_("Returns the background noise level")},
	
	{"cd", 1, "cd directory", process_cd, N_("Set the new current working directory. The argument \"directory\" can contain the ~ token, expanded as the home directory, directories with spaces in the name can be protected using single or double quotes")},
	{"cdg", 0, "cdg", process_cdg, N_("Returns the coordinates of the center of gravity of the image")},
	{"clearstar", 0, "clearstar", process_clearstar, N_("Clear all the stars saved in memory and displayed on the screen")},
	{"close", 0, "close", process_close, N_("Properly closes the opened image and the opened sequence, if any")},
	{"cosme", 1, "cosme [filename].lst", process_cosme, N_("Apply the local mean to a set of pixels on the in-memory image (cosmetic correction). The coordinate of this pixels are in an ASCII file [.lst file]. COSME is adapted to correct residual hot and cold pixels after preprocessing")},
	{"cosme_cfa", 1, "cosme_cfa [filename].lst", process_cosme, N_("Same function that COSME but applying to RAW CFA images")},
	{"crop", 0, "crop [x y width height]", process_crop, N_("Crops the current image within the rectangle previously selected")},

	{"ddp", 3, "ddp level coef sigma", process_ddp, N_("Performs a DDP (digital development processing) as described first by Kunihiko Okano. This implementation is the one described in IRIS. It combines a linear distribution on low levels (below \"level\") and a non-linear on high levels. It uses a Gaussian filter of sigma \"sigma\" and multiplies the resulting image by \"coef\". The typical values for \"sigma\" are included between 0.7 and 2")},
	
	{"entropy", 0, "entropy", process_entropy, N_("Computes the entropy of the opened image on the displayed layer, only in the selected area if one has been selected or in the whole image else. The entropy is one way of measuring the noise or the details in an image")},
	{"exit", 0, "exit", process_exit, N_("Quits the application")},
	{"extract", 1, "extract NbPlans", process_extract, N_("Extracts \"NbPlans\" planes of wavelet domain")},
	
	{"fdiv", 2, "fdiv filename scalar", process_fdiv, N_("Divides the image in memory by the image given in argument. The resulting image is multiplied by the value of the \"scalar\" argument. See also IDIV")},
	{"fftd", 2, "fftd modulus phase", process_fft, N_("Applies a Fast Fourier Transform to the image loaded in memory. \"Modulus\" and \"phase\" given in argument are saved in FITS files")},
	{"ffti", 2, "ffti modulus phase", process_fft, N_("Retrieves corrected image applying an inverse transformation. The \"modulus\" and \"phase\" used are the files given in argument")},
	{"fill", 1, "fill value", process_fill, N_("Fills the whole current image (or selection) with pixels having the \"value\" intensity")},
	{"fill2", 1, "fill2 value [x y width height]", process_fill2, N_("Same command than FILL but this is a symmetric fill of a region defined by the mouse. Used to process an image in the Fourier (FFT) domain")},
	{"find_cosme", 2, "find_cosme cold_sigma hot_sigma", process_findcosme, N_("Applies an automatic detection of cold and hot pixels following the thresholds written in arguments")},
	{"find_cosme_cfa", 2, "find_cosme_cfa cold_sigma hot_sigma", process_findcosme, N_("Same command than FIND_COSME but for monochromatic CFA images")},
	{"find_hot", 3, "find_hot filename cold_sigma hot_sigma", process_findhot, N_("Provides a list file \"filename\" (format text) in the working directory which contains the coordinates of the pixels which have an intensity \"hot_sigma\" times higher and \"cold_sigma\" lower than standard deviation. We generally use this command on a master-dark file")},
	{"findstar", 0, "findstar", process_findstar, N_("Detects stars having a level greater than a threshold computed by Siril. The algorithm is based on the publication of Mighell, K. J. 1999, in ASP Conf. Ser., Vol. 172, Astronomical Data Analysis Software and Systems VIII, eds. D. M. Mehringer, R. L. Plante, & D. A. Roberts (San Francisco: ASP), 317. After that, a PSF is applied and Siril rejects all detected structures that don't fulfill a set of prescribed detection criteria. Finaly, a circle is drawn around detected stars. See also the command CLEARSTAR")},
	{"fmedian", 2, "fmedian ksize modulation", process_fmedian, N_("Performs a median filter of size \"ksize\" x \"ksize\" (\"ksize\" MUST be odd) to the original image with a modulation parameter \"modulation\". The output pixel is computed as : out=mod x m + (1 − mod) x in, where m is the median-filtered pixel value. A modulation's value of 1 will apply no modulation")},
	{"fmul", 1, "fmul scalar", process_fmul, N_("Multiplies the loaded image by the \"scalar\" given in argument")},
	{"fixbanding", 2, "fixbanding amount sigma", process_fixbanding, N_("Tries to remove the canon banding. Argument \"amount\" define the amount of correction. \"Sigma\" defines a protection level of the algorithm, higher sigma gives higher protection")},

	{"gauss", 1, "gauss sigma", process_gauss, N_("Performs a Gaussian filter with the given \"sigma\"")},

	{"help", 0, "help", process_help, N_("Gives the available commands")},
	{"histo", 1, "histo channel (channel=0, 1, 2 with 0: red, 1: green, 2: blue)", process_histo, N_("Calculates the histogram of the image channel in memory and produces file histo_[channel name].dat in the working directory")},
	
	/* commands oper filename and curent image */
	{"iadd", 1, "iadd filename", process_imoper, N_("Adds the image in memory to the image given in argument")},
	{"idiv", 1, "idiv filename", process_imoper, N_("Divides the image in memory by the image given in argument. See also FDIV")},
	{"imul", 1, "imul filename", process_imoper, N_("Multiplies the image in memory by the image given in argument")},
	{"isub", 1, "isub filename", process_imoper, N_("Subtracts the image in memory by the image given in argument")},
	
	{"load", 1, "load filename.[ext]", process_load, N_("Loads the image \"filename\"; it first attempts to load \"filename\", then \"filename\".fit and finally \"filename\".fits and after, all supported format, aborting if none of these are found. These scheme is applicable to every Siril command implying reading files. Fits headers MIPS-HI and MIPS-LO are read and their values given to the current viewing levels. Writing a known extension at the end of \"filename\" will load the image \"filename\".ext: this is used when numerous files have the same name but not the same extension")},
	// specific loads are not required, but could be used to force the
	// extension to a higher priority in case two files with same basename
	// exist (stat_file() manages that priority order for now).
	{"log", 0, "log", process_log, N_("Computes and applies a logarithmic scale to the current image")}, /* logarifies current image */
#ifndef _WIN32
	{"ls", 0, "ls", process_ls, N_("Lists files and directories in the working directory")},
#endif
	
	{"mirrorx", 0, "mirrorx", process_mirrorx, N_("Rotates the image around a vertical axis")},
	{"mirrory", 0, "mirrory", process_mirrory, N_("Rotates the image around an horizontal axis")},
	
	{"new", 3, "new width height nb_channel", process_new, N_("Creates a new image filled with zeros with a size of \"width\" x \"height\". The image is in 16-bit format, and it contains \"nb_channel\" channels, \"nb_channel\" being 1 or 3. It is not saved, but displayed and can be saved afterwards")},
	{"nozero", 1, "nozero level (replaces null values by level)", process_nozero, N_("Replaces null values by \"level\" values. Useful before an idiv or fdiv operation")}, /* replaces null values by level */
	
	{"offset", 1, "offset value", process_offset, N_("Adds the constant \"value\" to the current image. This constant can take a negative value. As Siril uses unsigned FITS files, if the intensity of the pixel become negative its value is replaced by 0 and by 65535 (for a 16-bit file) if the pixel intensity overflows")},
	
	{"preprocess", 1, "preprocess sequencename [-bias=, -dark=, -flat=] [-cfa] [-debayer]", process_preprocess, N_("Preprocesses the sequence \"sequencename\" using bias, dark and flat given in argument. It is possible to specify if images are CFA for cosmetic correction purposes with the option \"-cfa\" and also to demosaic images at the end of the process with \"-debayer\"")},
	{"psf", 0, "psf", process_psf, N_("Performs a PSF (Point Spread Function) on the selected star")},
	
	{"register", 1, "register sequence [-drizzle]", process_register, N_("Performs geometric transforms on images of the sequence given in arguement so that they may be superimposed on the reference image. The output sequence name starts with the prefix \"r_\". Using stars for registration, this algorithm only works with deepsky images. The option \"-drizzle\" performs a x2 drizzle on the images")},
	{"resample", 1, "resample factor", process_resample, N_("Resamples image with a factor \"factor\"")},
	{"rl", 2, "rl iterations sigma", process_rl, N_("Restores an image using the Richardson-Lucy method. \"Iterations\" is the number of iterations to be performed (typically between 10 and 50). \"Sigma\" is the size of the kernel to be applied")},
	{"rmgreen", 1, "rmgreen type", process_scnr, N_("Chromatic noise reduction filter. It removes green noise in the current image. This filter is based on PixInsight's SCNR Average Neutral algorithm and it is the same filter used by HLVG plugin in Photoshop. \"Type\"=1 stands for Average Neutral Protection, while \"type\"=2 stands for Maximum Neutral Protection")},
	{"rotate", 1, "rotate degree", process_rotate, N_("Rotates the image of an angle of \"degree\" value")},
	{"rotatePi", 0, "rotatePi", process_rotatepi, N_("Rotates the image of an angle of 180° around its center. This is equivalent to the command \"ROTATE 180\" or \"ROTATE -180\"")},
	
	{"satu", 1, "satu coeff", process_satu, N_("Enhances the global saturation of the image. Try iteratively to obtain best results")},
	{"save", 1, "save filename", process_save, N_("Saves current image to \"filename\".fit. Fits headers MIPS-HI and MIPS-LO are added with values corresponding to the current viewing levels")},
	{"savebmp", 1, "savebmp filename", process_savebmp, N_("Saves current image under the form of a bitmap file with 8-bit per channel: \"filename\".bmp (BMP 24-bit). This function is like a screenshot of what you see with the levels updated. This is very usefull to share an image in the bitmap format")},
#ifdef HAVE_LIBJPEG
	{"savejpg", 1, "savejpg filename [quality]", process_savejpg, N_("Saves current image into a JPG file: \"filename\".jpg. You have the possibility to adjust the quality of the compression. A value 100 for \"quality\" parameter offers best fidelity while a low value increases the compression ratio. If no value is specified, it holds a value of 100")},
#endif
#ifdef HAVE_LIBPNG
	{"savepng", 1, "savepng filename", process_savepng, N_("Saves current image into a PNG file: \"filename\".png")},
#endif
	{"savepnm", 1, "savepnm filename", process_savepnm, N_("Saves current image under the form of a Netpbm file format with 16-bit per channel. The extension of the output will be \"filename\".ppm for RGB image and \"filename\".pgm for gray-level image")},
#ifdef HAVE_LIBTIFF
	{"savetif", 1, "savetif filename", process_savetif, N_("Saves current image under the form of a uncompressed TIFF file with 16-bit per channel: \"filename\".tif")},
	{"savetif8", 1, "savetif8 filename", process_savetif, N_("Same command than SAVE_TIF but the output file is saved in 8-bit per channel: \"filename\".tif")},
#endif
	{"select", 2, "select from to", process_select, N_("This command allows easy mass selection of images in the loaded sequence (\"from\" - \"to\", to included)")},
	{"seqcrop", 0, "seqcrop", process_seq_crop, N_("Crops the loaded sequence")},
	{"seqfind_cosme", 3, "seqfind_cosme sequencename cold_sigma hot_sigma", process_findcosme, N_("Same command than FIND_COSME but for the sequence \"sequencename\"")},
	{"seqfind_cosme_cfa", 3, "seqfind_cosme_cfa sequencename cold_sigma hot_sigma", process_findcosme, N_("Same command than FIND_COSME_CFA but for the sequence \"sequencename\"")},
	{"seqpsf", 0, "seqpsf", process_seq_psf, N_("Same command than PSF but works for sequences. Results are dumped in the console in a form that can be used to produce brightness variation curves")},
#ifdef _OPENMP
	{"setcpu", 1, "setcpu number", process_set_cpu, N_("Defines the number of processing threads used for calculation. Can be as high as the number of virtual threads existing on the system, which is the number of CPU cores or twice this number if hyperthreading (Intel HT) is available")},
#endif
	{"setmag", 1, "setmag magnitude", process_set_mag, N_("Calibrates the magnitude by selecting a star and giving the known apparent magnitude. All PSF computations will return the calibrated apparent magnitude afterwards, instead of an apparent magnitude relative to ADU values. To reset the magnitude constant see UNSETMAG")},
	{"setmagseq", 1, "setmagseq magnitude", process_set_mag_seq, N_("This command is only valid after having run SEQPSF or its graphical counterpart (select the area around a star and launch the PSF analysis for the sequence, it will appear in the graphs). This command has the same goal as SETMAG but recomputes the reference magnitude for each image of the sequence where the reference star has been found. When running the command, the last star that has been analysed will be considered as the reference star. Displaying the magnitude plot before typing the command makes it easy to understand. To reset the reference star and magnitude offset, see UNSETMAGSEQ")},
	{"split", 3, "split R G B", process_split, N_("Splits the color image into three distincts files (one for each color) and save them in \"r\" \"g\" and \"b\" file")},
	{"stack", 1, "stack sequencename [type] [sigma low] [sigma high] [-nonorm, norm=]", process_stackone, N_("Stacks the \"sequencename\" sequence, using options. The allowed types are: sum, max, min, med or median, and rej or mean that requires the use of additional arguments \"sigma low\" and \"high\" used for the Winsorized sigma clipping rejection algorithm (cannot be changed from here). If no argument other than the sequence name is provided, sum stacking is assumed")},
	{"stackall", 0, "stackall", process_stackall, N_("Opens all sequences in the CWD and stacks them with the optionally specified stacking type or with sum stacking. See STACK command for options description")},
	{"stat", 0, "stat", process_stat, N_("Returns global statistics of the current image. If a selection is made, the command returns statistics within the selection")},
	
	{"threshlo", 1, "threshlo level", process_threshlo, N_("Replaces values below \"level\" with \"level\"")},
	{"threshhi", 1, "threshi level", process_threshhi, N_("Replaces values above \"level\" with \"level\"")},
	{"thresh", 2, "thresh lo hi", process_thresh, N_("Replaces values below \"lo\" with \"lo\" and values above \"hi\" with \"hi\"")}, /* threshes hi and lo */
	
	/* unsharp masking of current image or genname sequence */
	{"unselect", 2, "unselect from to", process_unselect, N_("Allows easy mass unselection of images in the loaded sequence (\"from\" - \"to\"). See SELECT")},
	{"unsetmag", 0, "unsetmag", process_unset_mag, N_("Reset the magnitude calibration to 0. See SETMAG")},
	{"unsetmagseq", 0, "unsetmagseq", process_unset_mag_seq, N_("Resets the magnitude calibration and reference star for the sequence. See SETMAGSEQ")},
	{"unsharp", 2, "unsharp sigma multi", process_unsharp, N_("Applies to the working image an unsharp mask with sigma \"sigma\" and coefficient \"multi\"")},
	{"visu", 2, "visu low high", process_visu, N_("Displays an image with \"low\" and \"high\" as the low and high threshold")},
	
	/* wavelet transform in nbr_plan plans */ 
	{"wavelet", 1, "wavelet nbr_plan type", process_wavelet, N_("Computes the wavelet transform on \"nbr_plan\" plans using linear (type=1) or bspline (type=2) version of the 'a trous' algorithm. The result is stored in a file as a structure containing the planes, ready for weighted reconstruction with WRECONS")},
	/* reconstruct from wavelet transform and weighs plans with c1, c2, c3... */ 
	{"wrecons", 2, "wrecons c1 c2 c3 ...", process_wrecons, N_("Reconstructs to current image from the planes previously computed with wavelets and weighted with coefficients \"c1\", \"c2\", ..., \"cn\" according to the number of planes used for wavelet transform")},
	
	{"",0,"",0, ""}
};

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
	else if (single_image_is_loaded()) {
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

#ifdef HAVE_LIBPNG
int process_savepng(int nb){
	char filename[256];

	strcpy(filename, word[1]);
	strcat(filename, ".png");
	set_cursor_waiting(TRUE);
	savepng(filename, &gfit, 2, gfit.naxes[2] == 3);
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
	saveNetPBM(word[1], &gfit);
	return 0;
}

int process_imoper(int nb){
	fits fit = { 0 };

	if (readfits(word[1], &fit, NULL))
		return -1;
	imoper(&gfit, &fit, word[0][1]);
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	return 0;
}

int process_addmax(int nb){
	fits fit = { 0 };

	if (readfits(word[1], &fit, NULL))
		return -1;
	if (addmax(&gfit, &fit)==0) {
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

	norm = atof(word[2]);
	if (readfits(word[1], &fit, NULL))
		return -1;
	fdiv(&gfit, &fit, norm);
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

int process_rl(int nb) {

	double sigma;
	int iter;

	if (!com.headless)
		control_window_switch_to_tab(OUTPUT_LOGS);
	iter = atoi(word[1]);
	sigma = atof(word[2]);
	if (iter <= 0) {
		siril_log_message(_("Number of iterations must be > 0.\n"));
		return 1;
	}

	if (sigma <= 0) {
		siril_log_message(_("Sigma must be > 0.\n"));
		return 1;
	}

	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return 1;
	}

	struct RL_data *args = malloc(sizeof(struct RL_data));

	set_cursor_waiting(TRUE);

	args->fit = &gfit;
	args->sigma = sigma;
	args->iter = iter;

	start_in_new_thread(LRdeconv, args);

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

int process_cd(int nb) {
	char filename[256];
	int retval;

	g_strlcpy(filename, word[1], 250);

	expand_home_in_filename(filename, 256);
	retval = changedir(filename, NULL);
	if (!retval) {
		writeinitfile();
	}
	return retval;
}

int process_wrecons(int nb) {
	int i;
	float coef[7];
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];
	const char *tmpdir;
	int nb_chan = gfit.naxes[2];

	assert(nb_chan == 1 || nb_chan == 3);

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
	loglut(&gfit, LOG);
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
		struct stat entrystat;
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
		} else if (!strncmp(ext, "seq", 4))
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

int process_rotatepi(int nb){
	verbose_rotate_image(&gfit, 180.0, OPENCV_LINEAR, 1);

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
		fitted_PSF *result = psf_get_minimisation(&gfit, layer, &com.selection, TRUE);
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
		fitted_PSF *result = psf_get_minimisation(&gfit, layer, &com.selection, TRUE);
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
	if (!com.headless)
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

	FILE *f = g_fopen(name, "w");

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

	fits *fit = &gfit;
	if (new_fit_image(&fit, width, height, layers))
		return 1;
	memset(gfit.data, 0, width*height*layers*sizeof(WORD));

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
	double dirty;
	int is_cfa, i = 0, retval = 0;
	int nb_tokens;
	char line[64];
	char type;

	if (!ends_with(word[1], ".lst"))
		strcat(word[1], ".lst");
	cosme_file = g_fopen(word[1], "r");
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

int process_close(int nb) {
	free_image_data();
	close_sequence(TRUE);
	undo_flush();
	hide_rgb_window();
	hide_gray_window();
	reset_plot(); // reset all plots
	close_tab();	//close Green and Blue Tab if a 1-layer sequence is loaded
	update_used_memory();
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

		if (!existseq(file))
			check_seq(FALSE);
		seq = readseqfile(file);
		if (seq == NULL) {
			siril_log_message(_("No sequence %s found.\n"), file);
			return 1;
		}
		seq_check_basic_data(seq, FALSE);
		i++;
	} else {
		if (!single_image_is_loaded())
			return 1;
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
	char layername[6];

	for (layer = 0; layer < nplane; layer++) {
		imstats* stat = statistics(NULL, -1, &gfit, layer, &com.selection, STATS_MAIN);
		if (!stat) {
			siril_log_message(_("Error: no data computed.\n"));
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

int process_register(int nb) {
	struct registration_args *reg_args;
	struct registration_method *method;
	char *msg;

	if (get_thread_run()) {
		siril_log_message(_("Another task is "
				"already in progress, ignoring new request.\n"));
		return 1;
	}

	gchar *file = g_strdup(word[1]);
	if (!ends_with(file, ".seq")) {
		str_append(&file, ".seq");
	}

	if (!existseq(file))
		check_seq(FALSE);
	sequence *seq = readseqfile(file);
	if (seq == NULL) {
		siril_log_message(_("No sequence %s found.\n"), file);
		return 1;
	}
	seq_check_basic_data(seq, FALSE);

	/* getting the selected registration method */
	struct registration_method *reg = malloc(sizeof(struct registration_method));
	reg->name = strdup(_("Global Star Alignment (deep-sky)"));
	reg->method_ptr = &register_star_alignment;
	reg->sel = REQUIRES_NO_SELECTION;
	reg->type = REGTYPE_DEEPSKY;

	method = reg;

	reg_args = calloc(1, sizeof(struct registration_args));

	if (!com.headless)
		control_window_switch_to_tab(OUTPUT_LOGS);

	/* filling the arguments for registration */
	reg_args->func = method->method_ptr;
	reg_args->seq = seq;
	reg_args->process_all_frames = TRUE;
	reg_args->follow_star = FALSE;
	reg_args->matchSelection = FALSE;
	reg_args->translation_only = FALSE;
	if (word[2] && (!strcmp(word[2], "-drizzle")))
		reg_args->x2upscale = TRUE;
	else
		reg_args->x2upscale = FALSE;
	/* Here we should test available free disk space for Drizzle operation */
	if (reg_args->x2upscale) {
		double size = seq_compute_size(reg_args->seq);
		double diff = test_available_space(size * 4.0); //FIXME: 4 is only ok for x2 Drizzle
		if (diff < 0.0) {
			msg = siril_log_message(_("Not enough disk space to "
					"perform Drizzle operation!\n"));
			free(reg_args);
			return 1;
		}
	}
	/* getting the selected registration layer from the combo box. The value is the index
	 * of the selected line, and they are in the same order than layers so there should be
	 * an exact matching between the two */
	reg_args->layer = (reg_args->seq->nb_layers == 3) ? 1 : 0;
	reg_args->interpolation = OPENCV_LINEAR;
	get_the_registration_area(reg_args, method);	// sets selection
	reg_args->run_in_thread = TRUE;
	reg_args->prefix = "r_";
	reg_args->load_new_sequence = FALSE;	// don't load it for command line execution

	msg = siril_log_color_message(
			_("Registration: processing using method: %s\n"), "red",
			method->name);
	msg[strlen(msg) - 1] = '\0';
	gettimeofday(&(reg_args->t_start), NULL);
	set_progress_bar_data(msg, PROGRESS_RESET);

	start_in_new_thread(register_thread_func, reg_args);
	return 0;
}

struct _stackall_data {
	const gchar *file;
	stack_method method;
	double sig[2];
	gboolean force_no_norm;
	normalization norm;
	int number_of_loaded_sequences;
};

static int stack_one_seq(struct _stackall_data *arg) {
	int retval = -1;
	sequence *seq = readseqfile(arg->file);
	if (seq != NULL) {
		char filename[256];
		struct stacking_args args;
		if (seq_check_basic_data(seq, FALSE) == -1) {
			free(seq);
			return 1;
		}
		siril_log_message(_("Stacking sequence %s\n"), seq->seqname);
		args.seq = seq;
		args.filtering_criterion = stack_filter_all;
		args.filtering_parameter = 0.0;
		args.nb_images_to_stack = seq->number;
		args.image_indices = malloc(seq->number * sizeof(int));
		gettimeofday(&args.t_start, NULL);
		args.max_number_of_rows = stack_get_max_number_of_rows(seq, seq->number);
		// the three below: used only if method is average w/ rejection
		args.sig[0] = arg->sig[0];
		args.sig[1] = arg->sig[1];
		args.type_of_rejection = WINSORIZED;
		args.coeff.offset = NULL;
		args.coeff.mul = NULL;
		args.coeff.scale = NULL;
		if (!arg->force_no_norm &&
				(arg->method == stack_median || arg->method == stack_mean_with_rejection))
			args.normalize = arg->norm;
		else args.normalize = NO_NORM;
		args.force_norm = FALSE;
		args.norm_to_16 = TRUE;
		args.reglayer = args.seq->nb_layers == 1 ? 0 : 1;
		stack_fill_list_of_unfiltered_images(&args);

		char *suffix = ends_with(seq->seqname, "_") ? "" :
		       	(ends_with(com.seq.seqname, "-") ? "" : "_");
		snprintf(filename, 256, "%s%sstacked%s",
				seq->seqname, suffix, com.ext);

		// 1. normalization
		do_normalization(&args);	// does nothing if NO_NORM
		// 2. up-scale
		upscale_sequence(&args); // does nothing if args->seq->upscale_at_stacking <= 1.05
		// 3. stack
		retval = arg->method(&args);

		// TODO, from end_stacking:
		//_show_summary(args);
		//noise(&gfit);
		//remove_tmp_drizzle_files(args, TRUE);

		free_sequence(seq, TRUE);
		free(args.image_indices);
		if (!retval) {
			if (savefits(filename, &gfit))
				siril_log_color_message(_("Could not save the stacking result %s\n"),
						"red", filename);
			++arg->number_of_loaded_sequences;
		}
		else if (!get_thread_run()) return -1;

	} else {
		siril_log_message(_("No sequence %s found.\n"), arg->file);
	}
	return retval;
}

static gpointer stackall_worker(gpointer garg) {
	GDir *dir;
	GError *error = NULL;
	const gchar *file;
	struct _stackall_data *arg = (struct _stackall_data *)garg;

	if (!com.headless)
		control_window_switch_to_tab(OUTPUT_LOGS);
	siril_log_message(_("Looking for sequences in current working directory...\n"));
	if (check_seq(0) || (dir = g_dir_open(com.wd, 0, &error)) == NULL) {
		siril_log_message(_("Error while searching sequences or opening the directory.\n"));
		fprintf (stderr, "stackall: %s\n", error->message);
		com.wd[0] = '\0';
		siril_add_idle(end_generic, NULL);
		return NULL;
	}
	siril_log_message(_("Starting stacking of found sequences...\n"));
	arg->number_of_loaded_sequences = 0;
	while ((file = g_dir_read_name(dir)) != NULL) {
		char *suf;

		if ((suf = strstr(file, ".seq")) && strlen(suf) == 4) {
			arg->file = file;
			stack_one_seq(arg);
		}
	}
	g_dir_close(dir);
	free(arg);
	siril_log_message(_("Stacked %d sequences successfully.\n"), arg->number_of_loaded_sequences);
	siril_add_idle(end_generic, NULL);
	return NULL;
}

int process_stackall(int nb) {
	struct _stackall_data *arg = malloc(sizeof (struct _stackall_data));
	arg->force_no_norm = FALSE;
	arg->norm = ADDITIVE_SCALING;

	if (!word[1] || !strcmp(word[1], "sum"))
		arg->method = stack_summing_generic;
	else if (!strcmp(word[1], "max"))
		arg->method = stack_addmax;
	else if (!strcmp(word[1], "min"))
		arg->method = stack_addmin;
	else if (!strcmp(word[1], "med") || !strcmp(word[1], "median")) {
		arg->method = stack_median;
		if (word[2] && (!strcmp(word[2], "-nonorm") || !strcmp(word[2], "-no_norm")))
			arg->force_no_norm = TRUE;
		else if (word[2] && g_str_has_prefix(word[2], "-norm=")) {
			if (!strcmp(word[2], "-norm=add")) {
				arg->norm = ADDITIVE;
			} else if (!strcmp(word[2], "-norm=addscale")) {
				arg->norm = ADDITIVE_SCALING;
			} else if (!strcmp(word[2], "-norm=mul")) {
				arg->norm = MULTIPLICATIVE;
			} else if (!strcmp(word[2], "-norm=mulscale")) {
				arg->norm = MULTIPLICATIVE_SCALING;
			}
		}
	} else if (!strcmp(word[1], "rej") || !strcmp(word[1], "mean")) {
		if (!word[2] || !word[3] ||
				(arg->sig[0] = atof(word[2])) < 0.001 ||
				(arg->sig[1] = atof(word[3])) < 0.001) {
			siril_log_message(_("The average stacking with rejection uses the Winsorized "
					"rejection here and requires two extra arguments: sigma low and high.\n"));
			free(arg);
			return 1;
		}
		arg->method = stack_mean_with_rejection;
		if (word[4] && (!strcmp(word[4], "-nonorm") || !strcmp(word[4], "-no_norm")))
			arg->force_no_norm = TRUE;
		else if (word[4] && g_str_has_prefix(word[4], "-norm=")) {
			if (!strcmp(word[4], "-norm=add")) {
				arg->norm = ADDITIVE;
			} else if (!strcmp(word[4], "-norm=addscale")) {
				arg->norm = ADDITIVE_SCALING;
			} else if (!strcmp(word[4], "-norm=mul")) {
				arg->norm = MULTIPLICATIVE;
			} else if (!strcmp(word[4], "-norm=mulscale")) {
				arg->norm = MULTIPLICATIVE_SCALING;
			}
		}
	}
	else {
		siril_log_message(_("The provided type of stacking is unknown (%s).\n"), word[1]);
		free(arg);
		return 1;
	}

	start_in_new_thread(stackall_worker, arg);
	return 0;
}

static gpointer stackone_worker(gpointer garg) {
	char *suf;
	int retval = 0;
	struct _stackall_data *arg = (struct _stackall_data *)garg;

	siril_log_message(_("Looking for sequences in current working directory...\n"));
	if (check_seq(0)) {
		siril_log_message(_("Error while searching sequences.\n"));
		com.wd[0] = '\0';
		siril_add_idle(end_generic, NULL);
		return GINT_TO_POINTER(1);
	}

	if ((suf = strstr(arg->file, ".seq")) && strlen(suf) == 4) {
		retval = stack_one_seq(arg);
	}
	free(arg);
	if (!retval)
		siril_log_message(_("Stacked sequence successfully.\n"));
	siril_add_idle(end_generic, NULL);
	return NULL;
}

int process_stackone(int nb) {
	struct _stackall_data *arg = malloc(sizeof (struct _stackall_data));
	gchar *file;

	if (word[1][0] == '\0') {
		free(arg);
		return -1;
	}

	arg->force_no_norm = FALSE;

	file = g_strdup(word[1]);
	if (!ends_with(file, ".seq")) {
		str_append(&file, ".seq");
	}

	if (!existseq(file))
		check_seq(FALSE);
	sequence *seq = readseqfile(file);
	if (seq == NULL) {
		siril_log_message(_("No sequence %s found.\n"), file);
		free(arg);
		return 1;
	}
	seq_check_basic_data(seq, FALSE);

	arg->file = file;
	if (!word[2] || !strcmp(word[2], "sum"))
		arg->method = stack_summing_generic;
	else if (!strcmp(word[2], "max"))
		arg->method = stack_addmax;
	else if (!strcmp(word[2], "min"))
		arg->method = stack_addmin;
	else if (!strcmp(word[2], "med") || !strcmp(word[2], "median")) {
		arg->method = stack_median;
		if (word[3] && (!strcmp(word[3], "-nonorm") || !strcmp(word[3], "-no_norm")))
			arg->force_no_norm = TRUE;
		else if (word[3] && g_str_has_prefix(word[3], "-norm=")) {
			if (!strcmp(word[3], "-norm=add")) {
				arg->norm = ADDITIVE;
			} else if (!strcmp(word[3], "-norm=addscale")) {
				arg->norm = ADDITIVE_SCALING;
			} else if (!strcmp(word[3], "-norm=mul")) {
				arg->norm = MULTIPLICATIVE;
			} else if (!strcmp(word[3], "-norm=mulscale")) {
				arg->norm = MULTIPLICATIVE_SCALING;
			}
		}
	} else if (!strcmp(word[2], "rej") || !strcmp(word[2], "mean")) {
		if (!word[3] || !word[4] ||
				(arg->sig[0] = atof(word[3])) < 0.001 ||
				(arg->sig[1] = atof(word[4])) < 0.001) {
			siril_log_message(_("The average stacking with rejection uses the Winsorized "
					"rejection here and requires two extra arguments: sigma low and high.\n"));
			free(arg);
			return 1;
		}
		arg->method = stack_mean_with_rejection;
		if (word[5] && (!strcmp(word[5], "-nonorm") || !strcmp(word[5], "-no_norm")))
			arg->force_no_norm = TRUE;
		else if (word[5] && g_str_has_prefix(word[5], "-norm=")) {
			if (!strcmp(word[5], "-norm=add")) {
				arg->norm = ADDITIVE;
			} else if (!strcmp(word[5], "-norm=addscale")) {
				arg->norm = ADDITIVE_SCALING;
			} else if (!strcmp(word[5], "-norm=mul")) {
				arg->norm = MULTIPLICATIVE;
			} else if (!strcmp(word[5], "-norm=mulscale")) {
				arg->norm = MULTIPLICATIVE_SCALING;
			}
		}
	}
	else {
		siril_log_message(_("The provided type of stacking is unknown (%s).\n"), word[2]);
		free(arg);
		return 1;
	}

	start_in_new_thread(stackone_worker, arg);
	return 0;
}

// preprocess sequencename -bias= -dark= -flat= -cfa -debayer
int process_preprocess(int nb) {
	struct preprocessing_data *args = malloc(sizeof(struct preprocessing_data));

	com.preprostatus = 0;
	gboolean is_cfa = FALSE;
	gboolean do_debayer = FALSE;
	gchar *file;
	fits *master_bias = NULL;
	fits *master_dark = NULL;
	fits *master_flat = NULL;
	int i, retvalue = 0;

	if (word[1][0] == '\0') {
		free(args);
		return -1;
	}

	file = g_strdup(word[1]);
	if (!ends_with(file, ".seq")) {
		str_append(&file, ".seq");
	}

	if (!existseq(file))
		check_seq(FALSE);
	sequence *seq = readseqfile(file);
	if (seq == NULL) {
		siril_log_message(_("No sequence %s found.\n"), file);
		free(args);
		return 1;
	}
	seq_check_basic_data(seq, FALSE);

	for (i = 2; i < 6; i++) {
		if (word[i]) {
			if (g_str_has_prefix(word[i], "-bias=")) {
				master_bias = calloc(1, sizeof(fits));
				if (!readfits(word[i] + 6, master_bias, NULL)) {
					com.preprostatus |= USE_OFFSET;
					seq->offset = master_bias;
				} else {
					retvalue = 1;
					break;
				}
			} else if (g_str_has_prefix(word[i], "-dark=")) {
				master_dark = calloc(1, sizeof(fits));
				if (!readfits(word[i] + 6, master_dark, NULL)) {
					com.preprostatus |= USE_DARK;
					com.preprostatus |= USE_COSME;
					seq->dark = master_dark;
				} else {
					retvalue = 1;
					break;
				}
			} else if (g_str_has_prefix(word[i], "-flat=")) {
				master_flat = calloc(1, sizeof(fits));
				if (!readfits(word[i] + 6, master_flat, NULL)) {
					com.preprostatus |= USE_FLAT;
					seq->flat = master_flat;
				} else {
					retvalue = 1;
					break;
				}
			} else if (!strcmp(word[i], "-cfa")) {
				is_cfa = TRUE;
			}  else if (!strcmp(word[i], "-debayer")) {
				do_debayer = TRUE;
			}
		}
	}

	if (retvalue || !com.preprostatus) {
		if (master_bias) free(master_bias);
		if (master_dark) free(master_dark);
		if (master_flat) free(master_flat);
		free(args);
		return -1;
	}

	siril_log_color_message(_("Preprocessing...\n"), "red");
	gettimeofday(&args->t_start, NULL);

	/* Get parameters */
	args->seq = seq;
	args->autolevel = TRUE;
	args->normalisation = 1.0f;	// will be updated anyway

	args->sigma[0] = -1.00; /* cold pixels */
	args->sigma[1] =  3.00; /* hot poxels */

	args->compatibility = FALSE;

	args->debayer = do_debayer;
	args->is_cfa = is_cfa;

	args->offset = args->seq->offset;
	args->dark = args->seq->dark;
	args->flat = args->seq->flat;
	args->is_sequence = TRUE;

	/****/

	// sequence, executed in a background thread
	args->seq->ppprefix = strdup("pp_");

	// start preprocessing
	start_in_new_thread(seqpreprocess, args);

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
	undo_flush();
	exit(EXIT_SUCCESS);
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
	siril_log_color_message(_("Running command: %s\n"), "salmon", word[0]);
	return commande[i].process(wordnb);
}

gpointer execute_script(gpointer p) {
	FILE *fp = (FILE *)p;
	ssize_t read;
	char *linef;
	int line = 0, retval = 0;
	struct timeval t_start, t_end;

	com.headless = TRUE;
	com.stop_script = FALSE;
	gettimeofday(&t_start, NULL);
#if (_POSIX_C_SOURCE < 200809L)
	linef = calloc(256, sizeof(char));
	while (fgets(linef, 256, fp)) {
		read = strlen(linef) + 1;
#else
	size_t lenf = 0;
	linef = NULL;
	while ((read = getline(&linef, &lenf, fp)) != -1) {
#endif
		++line;
		if (com.stop_script) {
			retval = 1;
			break;
		}
		if (linef[0] == '#') {
			siril_log_color_message(linef, "blue");
			continue;
		}
		if (linef[0] == '\0' || linef[0] == '\n')
			continue;
		int wordnb;
		char *myline = strdup(linef);
		parseLine(myline, read, &wordnb);
		if (executeCommand(wordnb)) {
			siril_log_message(_("Error in line %d. Exiting batch processing\n"), line);
			free(myline);
			retval = 1;
			break;
		}
		free(myline);
		if (waiting_for_thread())
			break;	// abort script on command failure
	}
	free(linef);
	fclose(fp);
	com.headless = FALSE;
	com.stop_script = FALSE;
	if (!retval) {
		siril_log_message(_("Script execution finished successfully.\n"));
		gettimeofday(&t_end, NULL);
		show_time_msg(t_start, t_end, _("Total execution time"));
	}
	return GINT_TO_POINTER(retval);
}

int processcommand(const char *line) {
	int wordnb = 0, len;
	char *myline;

	if (line[0] == '\0' || line[0] == '\n')
		return 0;
	if (line[0] == '@') { // case of files
		if (get_thread_run()) {
			siril_log_message(_("Another task is already in progress, ignoring new request.\n"));
			return 1;
		}
		if (script_thread)
			g_thread_join(script_thread);
		FILE* fp = g_fopen(line + 1, "r");
		if (fp == NULL) {
			siril_log_message(_("File [%s] does not exist\n"), line + 1);
			return 1;
		}
		control_window_switch_to_tab(OUTPUT_LOGS);
		siril_log_message(_("Starting script %s\n"), line + 1);
		script_thread = g_thread_new("script", execute_script, fp);
	} else {
		myline = strdup(line);
		len = strlen(line);
		parseLine(myline, len, &wordnb);
		if (executeCommand(wordnb)) {
			free(myline);
			return 1;
		}
		free(myline);
	}
	return 0;
}

/* callback functions */

#define COMPLETION_COLUMN 0

static gboolean on_match_selected(GtkEntryCompletion *widget, GtkTreeModel *model,
		GtkTreeIter *iter, gpointer user_data) {
	const gchar *cmd;
	GtkEditable *e = (GtkEditable *) gtk_entry_completion_get_entry(widget);
	gchar *s = gtk_editable_get_chars(e, 0, -1);
	gint cur_pos = gtk_editable_get_position(e);
	gint p = cur_pos;
	gchar *end;
	gint del_end_pos = -1;

	gtk_tree_model_get(model, iter, COMPLETION_COLUMN, &cmd, -1);

	end = s + cur_pos;

	if (end) {
		del_end_pos = end - s + 1;
	} else {
		del_end_pos = cur_pos;
	}

	gtk_editable_delete_text(e, 0, del_end_pos);
	gtk_editable_insert_text(e, cmd, -1, &p);
	gtk_editable_set_position(e, p);

	return TRUE;
}

static gboolean completion_match_func(GtkEntryCompletion *completion,
		const gchar *key, GtkTreeIter *iter, gpointer user_data) {
	gboolean res = FALSE;
	char *tag = NULL;
	GtkTreeModel *model = gtk_entry_completion_get_model(completion);
	int column = gtk_entry_completion_get_text_column(completion);

	if (gtk_tree_model_get_column_type(model, column) != G_TYPE_STRING)
		return FALSE;

	gtk_tree_model_get(model, iter, column, &tag, -1);

	if (tag) {
		char *normalized = g_utf8_normalize(tag, -1, G_NORMALIZE_ALL);
		if (normalized) {
			char *casefold = g_utf8_casefold(normalized, -1);
			if (casefold) {
				res = g_strstr_len(casefold, -1, key) != NULL;
			}
			g_free(casefold);
		}
		g_free(normalized);
		g_free(tag);
	}

	return res;
}

void init_completion_command() {
	GtkEntryCompletion *completion = gtk_entry_completion_new();
	GtkListStore *model = gtk_list_store_new(1, G_TYPE_STRING);
	GtkTreeIter iter;
	gint i;
	GtkEntry *entry = GTK_ENTRY(lookup_widget("command"));

	gtk_entry_completion_set_text_column(completion, COMPLETION_COLUMN);
	gtk_entry_set_completion(entry, completion);
	gtk_entry_completion_set_inline_completion(completion, TRUE);
	gtk_entry_completion_set_popup_single_match(completion, FALSE);
	gtk_entry_completion_set_minimum_key_length(completion, 2);
	gtk_entry_completion_set_match_func(completion, completion_match_func, NULL, NULL);
	g_signal_connect(G_OBJECT(completion), "match-selected", G_CALLBACK(on_match_selected), NULL);

	/* Populate the completion database. */
	for (i = 0; i < sizeof(commande) / sizeof(command); i++) {
		gtk_list_store_append(model, &iter);
		gtk_list_store_set(model, &iter, COMPLETION_COLUMN, commande[i].name, -1);
	}
	gtk_entry_completion_set_model(completion, GTK_TREE_MODEL(model));
	g_object_unref(model);
}

void on_GtkCommandHelper_clicked(GtkButton *button, gpointer user_data) {
	GtkEntry *entry;
	GString *str;
	GtkWidget *popover;
	gchar **command_line;
	const gchar *text;
	gchar *helper = NULL;
	gint i;

	entry = GTK_ENTRY(lookup_widget("command"));
	text = gtk_entry_get_text(entry);
	if (*text != 0) {
		command_line = g_strsplit_set(text, " ", -1);

		for (i = 0; i < sizeof(commande) / sizeof(command); i++) {
			if (!g_ascii_strcasecmp(commande[i].name, command_line[0])) {
				gchar **token;

				token = g_strsplit_set(commande[i].usage, " ", -1);
				str = g_string_new(token[0]);
				str = g_string_prepend(str, "<span foreground=\"red\"><b>");
				str = g_string_append(str, "</b>");
				if (token[1] != NULL) {
					str = g_string_append(str, commande[i].usage + strlen(token[0]));
				}
				str = g_string_append(str, "</span>\n\n\t");
				str = g_string_append(str, _(commande[i].definition));
				helper = g_string_free(str, FALSE);
				g_strfreev(token);
			}
		}
		if (!helper) {
			helper = g_strdup(_("No help for this command"));
		}

		g_strfreev(command_line);

		popover = popover_new(lookup_widget("command"), helper);
#if GTK_MAJOR_VERSION >= 3 && GTK_MINOR_VERSION < 22
		gtk_widget_show(popover);
#else
		gtk_popover_popup(GTK_POPOVER(popover));
#endif
		g_free(helper);
	}
}
