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
#ifdef _WIN32
#include <windows.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>

#include "conversion.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/processing.h"
#include "io/films.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "algos/demosaicing.h"

static unsigned int supported_filetypes = 0;	// initialized by initialize_converters()

// only used by the converter and its GUI
unsigned int convflags = CONV1X3;	// default

// NULL-terminated array, initialized by initialize_converters(), used only by stat_file
char *supported_extensions[MAX_EXTENSIONS];

supported_raw_list supported_raw[] = {
	{"dng",	"Adobe"},
	{"mos",	"Aptus"},
	{"cr2",	"Canon"},
	{"crw",	"Canon"},
	{"bay",	"Casio"},		// Not tested
	{"erf",	"Epson"},
	{"raf",	"Fuji"},
	{"3fr",	"Hasselblad"},
	{"kdc",	"Kodak"},
	{"dcr",	"Kodak"},
	{"mef",	"Mamiya"},
	{"mrw",	"Minolta"},
	{"nef",	"Nikon"},
	{"nrw",	"Nikon"},
	{"orf",	"Olympus"},
	{"raw",	"Leica"},
	{"rw2",	"Panasonic"},
	{"pef",	"Pentax"},
	{"ptx",	"Pentax"},		// Not tested
	{"x3f",	"Sigma"},		// Not supported yet
	{"srw",	"Samsung"},
	{"arw",	"Sony"}
};

char *filter_pattern[] = {
	"RGGB",
	"BGGR",
	"GBRG",
	"GRBG",
	"RBGBRGGGRGGBGGBGGRBRGRBGGGBGGRGGRGGB", /* XTRANS */
	"GBGGRGRGRBGBGBGGRGGRGGBGBGBRGRGRGGBG",
	"GGRGGBGGBGGRBRGRBGGGBGGRGGRGGBRBGBRG"
};

static int film_conversion(const char *src_filename, int index, unsigned int *added_frames, struct ser_struct *ser_file, struct _convert_data *args);
static int ser_conversion(const char *src_filename, int index, unsigned int *added_frames, struct ser_struct *ser_file, struct _convert_data *args);

int get_nb_raw_supported() {
	return G_N_ELEMENTS(supported_raw);
}

/* This function is used with command line only */
void list_format_available() {
	puts("=======================================================");
	puts("[            Supported image file formats             ]");
	puts("=======================================================");
	puts("FITS\t(*.fit, *.fits, *.fts)");
	puts("BMP\t(*.bmp)");
	puts("NetPBM\t(*.ppm, *.pgm, *.pnm)");
	puts("PIC\t(*.pic)");
#ifdef HAVE_LIBRAW
	printf("RAW\t(");
	int i, nb_raw;

	nb_raw = get_nb_raw_supported();
	for (i = 0; i < nb_raw; i++) {
		printf("*.%s", supported_raw[i].extension);
		if (i != nb_raw - 1)
			printf(", ");
	}
	printf(")\n");
#endif

#ifdef HAVE_LIBTIFF
	puts("TIFF\t(*.tif, *.tiff)");
#endif
#ifdef HAVE_LIBJPEG
	puts("JPEG\t(*.jpg, *.jpeg)");
#endif
#ifdef HAVE_LIBPNG
	puts("PNG\t(*.png)");
#endif
#ifdef HAVE_LIBHEIF
	puts("HEIF\t(*.heic, *.heif)");
#endif
}

// input is destroot
static char *create_sequence_filename(const char *destroot, int counter, char *output, int outsize) {
	const char *ext = get_filename_ext(destroot);
	if (ext) {
		/* we need to insert a number before the extension */
		gchar *the_ext = g_strdup(ext);
		gchar *the_root = g_strdup(destroot);
		if (the_ext == the_root) {
			g_snprintf(output, outsize, "%s", destroot);
			g_free(the_ext);
			g_free(the_root);
			return output;
		}
		the_root[ext-destroot-1] = '\0';
		gchar last_char = the_root[strlen(the_root)-1];
		if (last_char == '-' || last_char == '_')
			g_snprintf(output, outsize, "%s%05d.%s", the_root, counter, the_ext);
		else 	g_snprintf(output, outsize, "%s_%05d.%s", the_root, counter, the_ext);
		g_free(the_ext);
		g_free(the_root);
	} else {
		/* create the file name with destroot_number */
		g_snprintf(output, outsize, "%s%05d", destroot, counter);
	}
	return output;
}

/* This function sets all default values of libraw settings in the com.raw_set
 * struct, as defined in the glade file.
 * When the ini file is read, the values of com.raw_set are overwritten, but if the
 * file is missing, like the first time Siril is launched, we don't want to have the
 * GUI states reset to zero by set_GUI_LIBRAW() because the data in com.raw_set had
 * not been initialized with the default GUI values (= initialized to 0).
 */
static void initialize_libraw_settings() {
	com.raw_set.bright = 1.0;		// brightness
	com.raw_set.mul[0] = 1.0;		// multipliers: red
	com.raw_set.mul[1] = 1.0;		// multipliers: green, not used because always equal to 1
	com.raw_set.mul[2] = 1.0;		// multipliers: blue
	com.raw_set.auto_mul = 1;		// multipliers are Either read from file, or calculated on the basis of file data, or taken from hardcoded constants
	com.raw_set.user_black = 0;		// black point correction
	com.raw_set.use_camera_wb = 0;	// if possible, use the white balance from the camera.
	com.raw_set.use_auto_wb = 0;		// use automatic white balance obtained after averaging over the entire image
	com.raw_set.user_qual = 1;		// type of interpolation. AHD by default
	com.raw_set.gamm[0] = 1.0;		// gamma curve: linear by default
	com.raw_set.gamm[1] = 1.0;
}

static void initialize_ser_debayer_settings() {
	com.debayer.open_debayer = FALSE;
	com.debayer.use_bayer_header = TRUE;
	com.debayer.compatibility = FALSE;
	com.debayer.bayer_pattern = BAYER_FILTER_RGGB;
	com.debayer.bayer_inter = BAYER_RCD;
	com.debayer.xbayeroff= 0;
	com.debayer.ybayeroff= 0;
}

static gboolean end_convert_idle(gpointer p) {
	struct _convert_data *args = (struct _convert_data *) p;
	struct timeval t_end;

	if (!args->retval && get_thread_run() && args->nb_converted > 1) {
		// load the sequence
		char *ppseqname = malloc(strlen(args->destroot) + 5);
		sprintf(ppseqname, "%s.seq", args->destroot);
		check_seq(0);
		update_sequences_list(ppseqname);
		free(ppseqname);
	}

	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_DONE);
	set_cursor_waiting(FALSE);
	gettimeofday(&t_end, NULL);
	show_time(args->t_start, t_end);
	stop_processing_thread();
	g_free(args->destroot);
	free(args);
	return FALSE;
}

/* from a fits object, save to file or files, based on the channel policy from convflags */
static int save_to_target_fits(fits *fit, const char *dest_filename) {
	if (convflags & CONV3X1) {	// an RGB image to 3 fits, one for each channel
		char filename[130];

		if (fit->naxis != 3) {
			siril_log_message(_("Saving to 3 FITS files cannot be done because the source image does not have three channels\n"));
			return 1;
		}
		sprintf(filename, "r_%s", dest_filename);
		if (save1fits16(filename, fit, RLAYER)) {
			printf("tofits: save1fit8 error, CONV3X1\n");
			return 1;
		}
		sprintf(filename, "g_%s", dest_filename);
		if (save1fits16(filename, fit, GLAYER)) {
			printf("tofits: save1fit8 error, CONV3X1\n");
			return 1;
		}
		sprintf(filename, "b_%s", dest_filename);
		if (save1fits16(filename, fit, BLAYER)) {
			printf("tofits: save1fit8 error, CONV3X1\n");
			return 1;
		}
	} else if (convflags & CONV1X1) { // a single FITS to convert from an RGB grey image
		if (save1fits16(dest_filename, fit, RLAYER)) {
			printf("tofits: save1fit8 error, CONV1X1\n");
			return 1;
		}
	} else {			// normal FITS save, any format
		if (savefits(dest_filename, fit)) {
			printf("tofits: savefit error, CONV1X3\n");
			return 1;
		}
	}
	return 0;
}

/* open the file with path source from any image type and load it into a new FITS object */
static fits *any_to_new_fits(image_type imagetype, const char *source, gboolean compatibility) {
	int retval = 0;
	fits *tmpfit = calloc(1, sizeof(fits));

	retval = any_to_fits(imagetype, source, tmpfit, FALSE, FALSE);

	if (!retval)
		retval = debayer_if_needed(imagetype, tmpfit, compatibility, FALSE);

	if (retval) {
		clearfits(tmpfit);
		free(tmpfit);
		return NULL;
	}

	return tmpfit;
}

int retrieveBayerPattern(char *bayer) {
	int i;

	for (i = 0; i < (sizeof(filter_pattern) / sizeof(char *)); i++) {
		if (g_ascii_strcasecmp(bayer, filter_pattern[i]) == 0) {
			return i;
		}
	}
	return BAYER_FILTER_NONE;
}

/**************************Public functions***********************************************************/

/* initialize converters (utilities used for different image types importing) *
 * updates the label listing the supported input file formats, and modifies the
 * list of file types used in convflags */
gchar *initialize_converters() {
	GString *string;
	gchar *text;
	int count_ext = 0;

	/* internal converters */
	supported_filetypes |= TYPEBMP;
	string = g_string_new(_("BMP images, "));
	supported_filetypes |= TYPEPIC;
	string = g_string_append(string, _("PIC images (IRIS), "));
	supported_filetypes |= TYPEPNM;
	string = g_string_append(string, _("PGM and PPM binary images"));

	supported_extensions[count_ext++] = ".fit";
	supported_extensions[count_ext++] = ".fits";
	supported_extensions[count_ext++] = ".fts";
	supported_extensions[count_ext++] = ".bmp";
	supported_extensions[count_ext++] = ".ppm";
	supported_extensions[count_ext++] = ".pgm";
	supported_extensions[count_ext++] = ".pnm";
	supported_extensions[count_ext++] = ".pic";

	initialize_ser_debayer_settings();	// below in the file

#ifdef HAVE_LIBRAW
	int i, nb_raw;

	supported_filetypes |= TYPERAW;
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("RAW images"));
	initialize_libraw_settings();	// below in the file

	nb_raw = get_nb_raw_supported();
	for (i = 0; i < nb_raw; i++) {
		// add the '.' to the extension
		char *ext = malloc(strlen(supported_raw[i].extension) + 2 * sizeof(char));
		sprintf(ext, ".%s", supported_raw[i].extension);
		supported_extensions[count_ext+i] = ext;
	}
	count_ext += nb_raw;
#endif
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("FITS-CFA images"));

#ifdef HAVE_FFMS2
	supported_filetypes |= TYPEAVI;
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("Films"));
#endif

	supported_filetypes |= TYPESER;
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("SER sequences"));

	/* library converters (detected by configure) */
#ifdef HAVE_LIBTIFF
	supported_filetypes |= TYPETIFF;
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("TIFF images"));
	supported_extensions[count_ext++] = ".tif";
	supported_extensions[count_ext++] = ".tiff";
#endif

#ifdef HAVE_LIBJPEG
	supported_filetypes |= TYPEJPG;
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("JPG images"));
	supported_extensions[count_ext++] = ".jpg";
	supported_extensions[count_ext++] = ".jpeg";
#endif

#ifdef HAVE_LIBPNG
	supported_filetypes |= TYPEPNG;
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("PNG images"));
	supported_extensions[count_ext++] = ".png";
#endif

#ifdef HAVE_LIBHEIF
	supported_filetypes |= TYPEHEIF;
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("HEIF images"));
	supported_extensions[count_ext++] = ".heic";
	supported_extensions[count_ext++] = ".heif";
#endif
	supported_extensions[count_ext++] = NULL;		// NULL-terminated array

	string = g_string_append(string, ".");
	text = g_string_free(string, FALSE);

	siril_log_message(_("Supported file types: %s\n"), text);
	return text;
}

int check_for_raw_extensions(const char *extension) {
	int i, nb_raw;
	nb_raw = get_nb_raw_supported();
	for (i = 0; i < nb_raw; i++) {
		if (!g_ascii_strcasecmp(extension, supported_raw[i].extension))
			return 0;
	}
	return 1;
}

/* returns the image_type for the extension without the dot, only if it is supported by
 * the current instance of Siril. */
image_type get_type_for_extension(const char *extension) {
	if ((supported_filetypes & TYPEBMP) && !g_ascii_strcasecmp(extension, "bmp")) {
		return TYPEBMP;
	} else if ((supported_filetypes & TYPEJPG) &&
			(!g_ascii_strcasecmp(extension, "jpg") || !g_ascii_strcasecmp(extension, "jpeg"))) {
		return TYPEJPG;
	} else if ((supported_filetypes & TYPEHEIF) &&
			(!g_ascii_strcasecmp(extension, "heic") || !g_ascii_strcasecmp(extension, "heif"))) {
		return TYPEHEIF;
	} else if ((supported_filetypes & TYPETIFF) &&
			(!g_ascii_strcasecmp(extension, "tif") || !g_ascii_strcasecmp(extension, "tiff"))) {
		return TYPETIFF;
	} else if ((supported_filetypes & TYPEPNG) && !g_ascii_strcasecmp(extension, "png")) {
		return TYPEPNG;
	} else if ((supported_filetypes & TYPEPNM) &&
			(!g_ascii_strcasecmp(extension, "pnm") || !g_ascii_strcasecmp(extension, "ppm") ||
			 !g_ascii_strcasecmp(extension, "pgm"))) {
		return TYPEPNM;
	} else if ((supported_filetypes & TYPEPIC) && !g_ascii_strcasecmp(extension, "pic")){
		return TYPEPIC;
	} else if ((supported_filetypes & TYPERAW) && !check_for_raw_extensions(extension)) {
		return TYPERAW;
#ifdef HAVE_FFMS2
		// check_for_film_extensions is undefined without FFMS2
	} else if ((supported_filetypes & TYPEAVI) && !check_for_film_extensions(extension)) {
		return TYPEAVI;
#endif
	} else if ((supported_filetypes & TYPESER) && !g_ascii_strcasecmp(extension, "ser")) {
		return TYPESER;
	} else if (!g_ascii_strcasecmp(extension, "fit") || !g_ascii_strcasecmp(extension, "fits") ||
			!g_ascii_strcasecmp(extension, "fts")) {
		return TYPEFITS;
	}
	return TYPEUNDEF; // not recognized or not supported
}

gpointer convert_thread_worker(gpointer p) {
	char dest_filename[128], msg_bar[256];
	double progress = 0.0;
	struct ser_struct ser_file;
	struct _convert_data *args = (struct _convert_data *) p;
	unsigned int frame_index = 0;
	int i;

	if (convflags & CONVDSTSER) {
		if (convflags & CONV3X1) {
			siril_log_color_message(_("SER output will take precedence over the one-channel per image creation option.\n"), "salmon");
			convflags &= ~CONV3X1;
		}

		if (!(convflags & CONVMULTIPLE)) {
			if (ser_create_file(args->destroot, &ser_file, TRUE, NULL)) {
				siril_log_message(_("Creating the SER file failed, aborting.\n"));
				goto clean_exit;
			}
		}
	}

	args->retval = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(dynamic,3) \
	if(!args->input_has_a_seq && ((convflags & CONVDSTSER) || fits_is_reentrant()))
	// we should run in parallel only when images are converted, not sequences
#endif
	for (i = 0; i < args->total; i++) {
		if (args->retval || !get_thread_run()) {
			continue;
		}

		gchar *src_filename = args->list[i];
		const char *src_ext = get_filename_ext(src_filename);
		image_type imagetype;
		int index = args->input_has_a_seq ? frame_index : args->start + i;

		gchar *name = g_utf8_strrchr(src_filename, strlen(src_filename), G_DIR_SEPARATOR);
		if (name)
			g_snprintf(msg_bar, 256, _("Converting %s..."), name + 1);
		else g_snprintf(msg_bar, 256, _("Converting %s..."), src_filename);

		imagetype = get_type_for_extension(src_ext);
		com.filter = (int) imagetype;
		if (imagetype == TYPEUNDEF) {
			char msg[512];
			char *title = siril_log_message(_("Filetype is not supported, cannot convert: %s\n"), src_ext);
			g_snprintf(msg, 512, _("File extension '%s' is not supported.\n"
				"Verify that you typed the extension correctly.\n"
				"If so, you may need to install third-party software to enable "
				"this file type conversion, look at the README file.\n"
				"If the file type you are trying to load is listed in supported "
				"formats, you may notify the developers that the extension you are "
				"trying to use should be recognized for this type."), src_ext);
			siril_message_dialog(GTK_MESSAGE_ERROR, title, msg);
			args->retval = 1;
			continue;
		}

		if (imagetype == TYPEAVI) {
			unsigned int added_frames;
			if (film_conversion(src_filename, index, &added_frames, &ser_file, args))
				args->retval = 1;
			frame_index += added_frames;
		}
		else if (imagetype == TYPESER) {
			unsigned int added_frames;
			if (ser_conversion(src_filename, index, &added_frames, &ser_file, args))
				args->retval = 1;
			frame_index += added_frames;
		}
		else {	// single image
			fits *fit = any_to_new_fits(imagetype, src_filename, args->compatibility);
			if (fit) {
				if (convflags & CONVDSTSER) {
					if (convflags & CONV1X1)
						keep_first_channel_from_fits(fit);
					if (ser_write_frame_from_fit(&ser_file, fit, i)) {
						siril_log_message(_("Error while converting to SER (no space left?)\n"));
						args->retval = 1;
					}
				} else {
					g_snprintf(dest_filename, 128, "%s%05d", args->destroot, index);
					if (save_to_target_fits(fit, dest_filename)) {
						siril_log_message(_("Error while converting to FITS (no space left?)\n"));
						args->retval = 1;
					}
				}
				clearfits(fit);
				free(fit);
			}
			frame_index++;
		}

#ifdef _OPENMP
#pragma omp atomic
#endif
		progress += 1.0;
		set_progress_bar_data(msg_bar, progress/((double)args->total));
#ifdef _OPENMP
#pragma omp atomic
#endif
		args->nb_converted++;
	}

clean_exit:
	if (convflags & CONVDSTSER) {
		if (!(convflags & CONVMULTIPLE)) {
			if (args->retval)
				ser_close_and_delete_file(&ser_file);
			else ser_write_and_close(&ser_file);
		}
	}
	if (args->command_line) {
		unset_debayer_in_convflags();
	}
	g_dir_close(args->dir);
	for (i = 0; i < args->total; i++)
		g_free(args->list[i]);
	if (args->retval)
		siril_log_message(_("Conversion ended with error, %d/%d input files converted\n"), args->nb_converted, args->total);
	else {
		if (args->nb_converted == args->total)
			siril_log_message(_("Conversion succeeded, %d/%d input files converted\n"), args->nb_converted, args->total);
		else siril_log_message(_("Conversion aborted, %d/%d input files converted\n"), args->nb_converted, args->total);
	}
	siril_add_idle(end_convert_idle, args);
	return NULL;
}

int debayer_if_needed(image_type imagetype, fits *fit, gboolean compatibility, gboolean force_debayer) {
	int retval = 0;
	sensor_pattern tmp;
	/* Siril's FITS are stored bottom-up, debayering will give wrong results.
	 * So before demosacaing we need to flip the image with fits_flip_top_to_bottom().
	 * But sometimes FITS are created by acquisition software top-down, in that case
	 * the user can indicate it ('compatibility') and we don't flip for debayer.
	 */
	if (imagetype == TYPEFITS && (((convflags & CONVDEBAYER) && !force_debayer) || force_debayer)) {
		tmp = com.debayer.bayer_pattern;
		if (fit->naxes[2] != 1) {
			siril_log_message(_("Cannot perform debayering on image with more than one channel\n"));
			return retval;
		}
		if (!compatibility)
			fits_flip_top_to_bottom(fit);
		/* Get Bayer informations from header if available */
		if (com.debayer.use_bayer_header) {
			sensor_pattern bayer;
			bayer = retrieveBayerPattern(fit->bayer_pattern);

			if (bayer <= BAYER_FILTER_MAX) {
				if (bayer != com.debayer.bayer_pattern) {
					if (bayer == BAYER_FILTER_NONE) {
						siril_log_color_message(_("No Bayer pattern found in the header file.\n"), "red");
					}
					else {
						siril_log_color_message(_("Bayer pattern found in header (%s) is different"
								" from Bayer pattern in settings (%s). Overriding settings.\n"),
								"red", filter_pattern[bayer], filter_pattern[com.debayer.bayer_pattern]);
						com.debayer.bayer_pattern = bayer;
					}
				}
			} else {
				com.debayer.bayer_pattern = XTRANS_FILTER;
				com.debayer.bayer_inter = XTRANS;
				siril_log_color_message(_("XTRANS Sensor detected. Using special algorithm.\n"), "red");
			}
		}
		if (com.debayer.bayer_pattern >= BAYER_FILTER_MIN
				&& com.debayer.bayer_pattern <= BAYER_FILTER_MAX) {
			siril_log_message(_("Filter Pattern: %s\n"), filter_pattern[com.debayer.bayer_pattern]);
		}

		if (debayer(fit, com.debayer.bayer_inter, com.debayer.bayer_pattern)) {
			siril_log_message(_("Cannot perform debayering\n"));
			retval = -1;
		} else {
			if (!compatibility)
				fits_flip_top_to_bottom(fit);
		}
		com.debayer.bayer_pattern = tmp;
	}
	return retval;
}

static int film_conversion(const char *src_filename, int index, unsigned int *added_frames, struct ser_struct *ser_file, struct _convert_data *args) {
	// we need to do a semi-recursive thing here,
	// thankfully it's only one level deep
	*added_frames = 0;
#ifdef HAVE_FFMS2
	int frame, ser_frames = index;
	int retval = 0;
	struct film_struct film_file;
	if (film_open_file(src_filename, &film_file) != FILM_SUCCESS) {
		siril_log_message(_("Error while opening film %s, aborting.\n"), src_filename);
		return 1;
	}
	if (convflags & CONVMULTIPLE) {
		char dest_filename[128];
		if (ser_create_file(create_sequence_filename(args->destroot, index, dest_filename, 128),
					ser_file, TRUE, NULL)) {
			siril_log_message(_("Creating the SER file failed, aborting.\n"));
			return 1;
		}
		ser_frames = 0;
	}
	for (frame = 0; frame < film_file.frame_count; frame++) {
		if (!get_thread_run() || retval) {
			break;
		}
		// read frame from the film
		fits fit = { 0 };
		if (film_read_frame(&film_file, frame, &fit) != FILM_SUCCESS) {
			siril_log_message(_("Error while reading frame %d from %s, aborting.\n"),
					frame, src_filename);
			retval = 1;
			break;	// no need to call clearfits
		}

		// save to the destination file
		if (convflags & CONVDSTSER) {
			if (convflags & CONV1X1)
				keep_first_channel_from_fits(&fit);
			if (ser_write_frame_from_fit(ser_file, &fit, frame + ser_frames)) {
				siril_log_message(_("Error while converting to SER (no space left?)\n"));
				retval = 1;
			}
		} else {
			char dest_filename[128];
			g_snprintf(dest_filename, 128, "%s%05d", args->destroot, index++);
			if (save_to_target_fits(&fit, dest_filename)) {
				siril_log_message(_("Error while converting to FITS (no space left?)\n"));
				retval = 1;
			}
		}
		clearfits(&fit);
	}
	if (convflags & CONVMULTIPLE) {
		if (retval)
			ser_close_and_delete_file(ser_file);
		else ser_write_and_close(ser_file);
	}
	*added_frames = frame;
	return retval;
#else
	return 1;
#endif
}

static int ser_conversion(const char *src_filename, int index, unsigned int *added_frames, struct ser_struct *ser_file, struct _convert_data *args) {
	int frame, ser_frames = index;
	int retval = 0;
	*added_frames = 0;

	struct ser_struct tmp_ser;
	ser_init_struct(&tmp_ser);
	siril_log_message(_("Reading %s\n"), src_filename);
	if (ser_open_file(src_filename, &tmp_ser)) {
		siril_log_message(_("Error while opening ser file %s, aborting.\n"), src_filename);
		return 1;
	}
	if (args->nb_converted > 0 && (convflags & CONVDSTSER)) {
		if (tmp_ser.image_height != ser_file->image_height
				|| tmp_ser.image_width != ser_file->image_width) {
			siril_log_color_message(_("Input SER files must have the same size to be joined.\n"), "red");
			return 1;
		}
	}
	/* We want to copy header informations from the first SER file */
	if (ser_file && ser_frames == 0) {
		/* TODO: this code should be in io/ser.c, common with ser_create_file() */
		ser_file->lu_id = tmp_ser.lu_id;
		ser_file->date = tmp_ser.date;
		ser_file->date_utc = tmp_ser.date_utc;
		ser_file->file_id = strdup(tmp_ser.file_id);
		memcpy(ser_file->observer, tmp_ser.observer, 40);
		memcpy(ser_file->instrument, tmp_ser.instrument, 40);
		memcpy(ser_file->telescope, tmp_ser.telescope, 40);
		ser_file->byte_pixel_depth = tmp_ser.byte_pixel_depth;
		ser_file->color_id = tmp_ser.color_id;
	}
	set_progress_bar_data(NULL, PROGRESS_PULSATE);
	for (frame = 0; frame < tmp_ser.frame_count; frame++) {
		if (!get_thread_run() || retval) {
			break;
		}
		// read frame from the sequence
		fits fit = { 0 };
		if (ser_read_frame(&tmp_ser, frame, &fit)) {
			siril_log_message(_("Error while reading frame %d from %s, aborting.\n"),
					frame, src_filename);
			retval = 1;
			break;	// no need to call clearfits
		}

		// save to the destination file
		if (convflags & CONVDSTSER) {
			if (convflags & CONV1X1)
				keep_first_channel_from_fits(&fit);
			if (ser_write_frame_from_fit(ser_file, &fit, frame + ser_frames)) {
				siril_log_message(_("Error while converting to SER (no space left?)\n"));
				retval = 1;
			}
		} else {
			char dest_filename[128];
			g_snprintf(dest_filename, 128, "%s%05d", args->destroot, index++);
			if (save_to_target_fits(&fit, dest_filename)) {
				siril_log_message(_("Error while converting to FITS (no space left?)\n"));
				retval = 1;
			}
		}
		clearfits(&fit);
	}
	*added_frames = frame;
	ser_close_file(&tmp_ser);
	return retval;
}


#ifdef _WIN32
char* g_real_path(const char *source) {
	HANDLE hFile;
	DWORD maxchar = 2048;
	TCHAR *FilePath;
	gchar *gFilePath;

	if (!(GetFileAttributesA(source) & FILE_ATTRIBUTE_REPARSE_POINT)) { /* Ce n'est pas un lien symbolique , je sors */
		return NULL;
	}

	FilePath = malloc(maxchar + 1);
	if (!FilePath) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	FilePath[0] = 0;

	hFile = CreateFile(source, GENERIC_READ, FILE_SHARE_READ, NULL,
			OPEN_EXISTING, 0, NULL);
	if (hFile == INVALID_HANDLE_VALUE) {
		free(FilePath);
		return NULL;
	}
	GetFinalPathNameByHandleA(hFile, FilePath, maxchar, 0);
	gFilePath = g_locale_to_utf8(FilePath + 4, -1, NULL, NULL, NULL); // +4 = enleve les 4 caracteres du prefixe "//?/"
	CloseHandle(hFile);
	return gFilePath;
}
#endif

/* open the file with path source from any image type and load it into the given FITS object */
int any_to_fits(image_type imagetype, const char *source, fits *dest, gboolean interactive, gboolean force_float) {
	int retval = 0;

	switch (imagetype) {
		case TYPEFITS:
			retval = (readfits(source, dest, NULL, force_float) != 0);
			break;
		case TYPEBMP:
			retval = (readbmp(source, dest) < 0);
			break;
		case TYPEPIC:
			retval = (readpic(source, dest) < 0);
			break;
#ifdef HAVE_LIBTIFF
		case TYPETIFF:
			retval = (readtif(source, dest) < 0);
			break;
#endif
		case TYPEPNM:
			retval = (import_pnm_to_fits(source, dest) < 0);
			break;
#ifdef HAVE_LIBJPEG
		case TYPEJPG:
			retval = (readjpg(source, dest) < 0);
			break;		
#endif
#ifdef HAVE_LIBHEIF
		case TYPEHEIF:
			retval = (readheif(source, dest, interactive) < 0);
			break;
#endif
#ifdef HAVE_LIBPNG
		case TYPEPNG:
			retval = (readpng(source, dest) < 0);
			break;
#endif
#ifdef HAVE_LIBRAW
		case TYPERAW:
			{
				const char *src = source ;
#ifdef _WIN32
				char *rsrc = g_real_path( source ) ;
				if ( rsrc != NULL )
				{
					src  = rsrc;					
				}
#endif
				retval = (open_raw_files(src , dest, !(convflags & CONVDEBAYER)) < 0);
#ifdef _WIN32
				if ( rsrc != NULL )
				{
					g_free( rsrc ) ;					
				}
#endif
			}
			break;
#endif
		case TYPESER:
		case TYPEAVI:
			siril_log_message(_("Requested converting a sequence file to single FITS image, should not happen\n"));
			retval = 1;
			break;
		case TYPEUNDEF:
		default:	// when the ifdefs are not compiled, default happens!
			siril_log_message(_("Error opening %s: file type not supported.\n"), source);
			retval = 1;
	}

	return retval;
}

void set_debayer_in_convflags() {
	convflags |= CONVDEBAYER;
}

void unset_debayer_in_convflags() {
	convflags &= ~CONVDEBAYER;
}

