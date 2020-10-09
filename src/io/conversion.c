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

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

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

#include "core/siril.h"
#include "core/proto.h"
#include "core/OS_utils.h"
#include "core/processing.h"
#include "io/films.h"
#include "io/fits_sequence.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/seqwriter.h"
#include "io/ser.h"
#include "io/FITS_symlink.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "algos/demosaicing.h"
#include "conversion.h"

#ifdef HAVE_LIBRAW
#include <libraw/libraw_version.h>
#endif

static unsigned int supported_filetypes = 0;	// initialized by initialize_converters()

// NULL-terminated array, initialized by initialize_converters(), used only by stat_file
char *supported_extensions[MAX_EXTENSIONS];

supported_raw_list supported_raw[] = {
	{"dng",	"Adobe"},
	{"mos",	"Aptus"},
	{"cr2",	"Canon"},
	{"crw",	"Canon"},
#ifdef HAVE_LIBRAW
#if LIBRAW_VERSION > LIBRAW_MAKE_VERSION(0, 19, 5)
	{"cr3",	"Canon"},
#endif
#endif
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

/* XTRANS */
	"GGRGGB" // ----> XTRANS_1
	"GGBGGR"
	"BRGRBG"
	"GGBGGR"
	"GGRGGB"
	"RBGBRG",

	"RBGBRG" // ----> XTRANS_2
	"GGRGGB"
	"GGBGGR"
	"BRGRBG"
	"GGBGGR"
	"GGRGGB",

	"GRGGBG"
	"BGBRGR"
	"GRGGBG"
	"GBGGRG"
	"RGTBGB"
	"GBGGRG",

	"GBGGRG"
	"RGRBGB"
	"GBGGRG"
	"GRGGBG"
	"BGBRGR"
	"GRGGBG"
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
		the_root[ext - destroot - 1] = '\0';
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
	com.pref.raw_set.bright = 1.0;		// brightness
	com.pref.raw_set.mul[0] = 1.0;		// multipliers: red
	com.pref.raw_set.mul[1] = 1.0;		// multipliers: green, not used because always equal to 1
	com.pref.raw_set.mul[2] = 1.0;		// multipliers: blue
	com.pref.raw_set.auto_mul = 1;		// multipliers are Either read from file, or calculated on the basis of file data, or taken from hardcoded constants
	com.pref.raw_set.user_black = 0;		// black point correction
	com.pref.raw_set.use_camera_wb = 0;	// if possible, use the white balance from the camera.
	com.pref.raw_set.use_auto_wb = 0;		// use automatic white balance obtained after averaging over the entire image
	com.pref.raw_set.user_qual = 2;		// type of interpolation. VNG by default
	com.pref.raw_set.gamm[0] = 1.0;		// gamma curve: linear by default
	com.pref.raw_set.gamm[1] = 1.0;
}

static void initialize_ser_debayer_settings() {
	com.pref.debayer.open_debayer = FALSE;
	com.pref.debayer.use_bayer_header = TRUE;
	com.pref.debayer.top_down = TRUE;
	com.pref.debayer.bayer_pattern = BAYER_FILTER_RGGB;
	com.pref.debayer.bayer_inter = BAYER_RCD;
	com.pref.debayer.xbayeroff = 0;
	com.pref.debayer.ybayeroff = 0;
}

static gboolean end_convert_idle(gpointer p) {
	struct _convert_data *args = (struct _convert_data *) p;
	struct timeval t_end;

	if (!args->retval && get_thread_run() && args->nb_converted_files > 1) {
		// load the sequence
		char *converted_seqname = NULL;
		if (args->output_type != SEQ_REGULAR) {
			int extidx = get_extension_index(args->destroot);
			if (extidx) {
				converted_seqname = malloc(extidx + 5);
				strncpy(converted_seqname, args->destroot, extidx);
				strcpy(converted_seqname+extidx, ".seq");
			}
		} else {
			converted_seqname = malloc(strlen(args->destroot) + 5);
			sprintf(converted_seqname, "%s.seq", args->destroot);
		}
		check_seq(0);
		if (converted_seqname) {
			update_sequences_list(converted_seqname);
			free(converted_seqname);
		}
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

/* open the file with path source from any image type and load it into a new FITS object */
static fits *any_to_new_fits(image_type imagetype, const char *source, gboolean debayer) {
	int retval = 0;
	fits *tmpfit = calloc(1, sizeof(fits));

	retval = any_to_fits(imagetype, source, tmpfit, FALSE, FALSE, debayer);

	if (!retval)
		retval = debayer_if_needed(imagetype, tmpfit, debayer);

	if (retval) {
		clearfits(tmpfit);
		free(tmpfit);
		return NULL;
	}

	return tmpfit;
}

/**************************Public functions***********************************************************/

int retrieveBayerPatternFromChar(char *bayer) {
	for (int i = 0; i < G_N_ELEMENTS(filter_pattern); i++) {
		if (g_ascii_strcasecmp(bayer, filter_pattern[i]) == 0) {
			return i;
		}
	}
	return BAYER_FILTER_NONE;
}

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

/* the list of files must be checked for unsupported file types before calling this
 * destroot must contain the file name with extension in case of SER or FITS sequence output. */
gpointer convert_thread_worker(gpointer p) {
	double progress = 0.0;
	struct ser_struct ser_file;
	fitseq fitseq_file;
	struct _convert_data *args = (struct _convert_data *) p;
	unsigned int frame_index = 0;

	args->nb_converted_files = 0;
	args->retval = 0;

	gboolean allow_symlink = args->output_type == SEQ_REGULAR && test_if_symlink_is_ok();

	if (args->output_type == SEQ_SER) {
		if (!args->multiple_output) {
			ser_init_struct(&ser_file);
			if (ser_create_file(args->destroot, &ser_file, TRUE, NULL)) {
				siril_log_message(_("Creating the SER file failed, aborting.\n"));
				goto clean_exit;
			}
		}
	}
	else if (args->output_type == SEQ_FITSEQ) {
		if (fitseq_create_file(args->destroot, &fitseq_file,
					args->input_has_a_seq ? -1 : args->total)) {
			siril_log_message(_("Creating the FITS sequence file failed, aborting.\n"));
			args->retval = 1;
			goto clean_exit;
		}

		int limit = 0;
#ifdef _OPENMP
		/* we don't know the image size here, max * 2 + 1 may be ok for most users */
		limit = com.max_thread * 2 + 1;
#endif
		seqwriter_set_max_active_blocks(limit);
	}
	else g_assert(args->output_type == SEQ_REGULAR);

	if (args->multiple_output && args->output_type != SEQ_SER)
		args->multiple_output = FALSE;

	gboolean have_seqwriter = args->output_type == SEQ_FITSEQ || args->output_type == SEQ_SER;
#ifdef _OPENMP
	if (have_seqwriter)
		omp_set_schedule(omp_sched_dynamic, 1);
	else omp_set_schedule(omp_sched_static, 0);
#pragma omp parallel for num_threads(com.max_thread) schedule(runtime) \
	if(!args->input_has_a_seq && (have_seqwriter || fits_is_reentrant()))
	// we should run in parallel only when images are converted, not sequences
#endif
	for (int i = 0; i < args->total; i++) {
		if (args->retval || !get_thread_run()) {
			continue;
		}

		gchar *src_filename = args->list[i];
		const char *src_ext = get_filename_ext(src_filename);
		int index = args->input_has_a_seq ? frame_index : args->start + i;

		gchar *name = g_utf8_strrchr(src_filename, strlen(src_filename), G_DIR_SEPARATOR);
		gchar *msg_bar;
		if (name)
			msg_bar = g_strdup_printf(_("Converting %s..."), name + 1);
		else msg_bar = g_strdup_printf(_("Converting %s..."), src_filename);

		image_type imagetype = get_type_for_extension(src_ext);
		com.filter = (int) imagetype;
		if (imagetype == TYPEUNDEF) {
			args->retval = 1;
			g_free(msg_bar);
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
			if (imagetype == TYPEFITS && args->make_link && allow_symlink) {
				gchar *dest_filename = g_strdup_printf("%s%05d%s", args->destroot, index,
						com.pref.ext);
				symlink_uniq_file(src_filename, dest_filename, allow_symlink);
				g_free(dest_filename);
			} else {
				if (have_seqwriter) {
					seqwriter_wait_for_memory();
				}

				fits *fit = any_to_new_fits(imagetype, src_filename, args->debayer);
				if (args->output_type == SEQ_SER) {
					if (ser_write_frame_from_fit(&ser_file, fit, i)) {
						siril_log_message(_("Error while converting to SER (no space left?)\n"));
						args->retval = 1;
					}
				} else if (args->output_type == SEQ_FITSEQ) {
					if (fitseq_write_image(&fitseq_file, fit, i)) {
						siril_log_message(_("Error while converting to SER (no space left?)\n"));
						args->retval = 1;
					}
				} else {
					if (fit) {
						gchar *dest_filename = g_strdup_printf("%s%05d", args->destroot, index);
						if (savefits(dest_filename, fit)) {
							siril_log_message(_("Error while converting to FITS (no space left?)\n"));
							args->retval = 1;
						}
						clearfits(fit);
						free(fit);
						g_free(dest_filename);
					}
				}
			}
			frame_index++;
		}

#ifdef _OPENMP
#pragma omp atomic
#endif
		progress += 1.0;
		set_progress_bar_data(msg_bar, progress / ((double) args->total));
		g_free(msg_bar);
#ifdef _OPENMP
#pragma omp atomic
#endif
		args->nb_converted_files++;
	}

clean_exit:
	if (args->output_type == SEQ_SER) {
		if (!args->multiple_output) {
			if (args->retval)
				ser_close_and_delete_file(&ser_file);
			else ser_write_and_close(&ser_file);
		}
	}
	else if (args->output_type == SEQ_FITSEQ) {
		if (args->retval)
			fitseq_close_and_delete_file(&fitseq_file);
		else fitseq_close_file(&fitseq_file);
	}
	for (int i = 0; i < args->total; i++)
		g_free(args->list[i]);
	free(args->list);
	if (args->retval)
		siril_log_message(_("Conversion ended with error, %d/%d input files converted\n"), args->nb_converted_files, args->total);
	else {
		if (args->nb_converted_files == args->total)
			siril_log_message(_("Conversion succeeded, %d/%d input files converted\n"), args->nb_converted_files, args->total);
		else siril_log_message(_("Conversion aborted, %d/%d input files converted\n"), args->nb_converted_files, args->total);
	}
	siril_add_idle(end_convert_idle, args);
	return NULL;
}

// debayers the image if it's a FITS image and if debayer is activated globally
// or if the force argument is passed
int debayer_if_needed(image_type imagetype, fits *fit, gboolean force_debayer) {
	if (imagetype != TYPEFITS || (!com.pref.debayer.open_debayer && !force_debayer))
		return 0;

	/* Siril's FITS are stored bottom-up, debayering will give wrong results.
	 * So before demosacaing we need to flip the image with fits_flip_top_to_bottom().
	 * But sometimes FITS are created by acquisition software top-down, in that case
	 * the user can indicate it ('compatibility') and we don't flip for debayer.
	 */
	if (fit->naxes[2] != 1) {
		siril_log_message(_("Cannot perform debayering on image with more than one channel\n"));
		return 0;
	}

	/* Get Bayer informations from header if available */
	sensor_pattern tmp_pattern = com.pref.debayer.bayer_pattern;
	interpolation_method tmp_algo = com.pref.debayer.bayer_inter;
	if (com.pref.debayer.use_bayer_header) {
		sensor_pattern bayer;
		bayer = retrieveBayerPatternFromChar(fit->bayer_pattern);

		if (bayer <= BAYER_FILTER_MAX) {
			if (bayer != tmp_pattern) {
				if (bayer == BAYER_FILTER_NONE) {
					siril_log_color_message(_("No Bayer pattern found in the header file.\n"), "red");
				}
				else {
					siril_log_color_message(_("Bayer pattern found in header (%s) is different"
								" from Bayer pattern in settings (%s). Overriding settings.\n"),
							"red", filter_pattern[bayer], filter_pattern[com.pref.debayer.bayer_pattern]);
					tmp_pattern = bayer;
				}
			}
		} else {
			tmp_pattern = bayer;
			tmp_algo = XTRANS;
			siril_log_color_message(_("XTRANS Sensor detected. Using special algorithm.\n"), "green");
		}
	}
	if (tmp_pattern >= BAYER_FILTER_MIN && tmp_pattern <= BAYER_FILTER_MAX) {
		siril_log_message(_("Filter Pattern: %s\n"),
				filter_pattern[tmp_pattern]);
	}

	int retval = 0;
	if (debayer(fit, tmp_algo, tmp_pattern)) {
		siril_log_message(_("Cannot perform debayering\n"));
		retval = -1;
	}
	return retval;
}

static int film_conversion(const char *src_filename, int index,
		unsigned int *added_frames, struct ser_struct *ser_file,
		struct _convert_data *args) {
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
	if (args->multiple_output) {
		char dest_filename[128];
		ser_init_struct(ser_file);
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
		if (args->output_type == SEQ_SER) {
			if (ser_write_frame_from_fit(ser_file, &fit, frame + ser_frames)) {
				siril_log_message(_("Error while converting to SER (no space left?)\n"));
				retval = 1;
			}
		} else {
			// TODO: handle FITS sequence
			char dest_filename[128];
			g_snprintf(dest_filename, 128, "%s%05d", args->destroot, index++);
			if (savefits(dest_filename, &fit)) {
				siril_log_message(_("Error while converting to FITS (no space left?)\n"));
				retval = 1;
			}
		}
		clearfits(&fit);
	}
	if (args->multiple_output) {
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

static int ser_conversion(const char *src_filename, int index,
		unsigned int *added_frames, struct ser_struct *ser_file,
		struct _convert_data *args) {
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
	if (args->nb_converted_files > 0 && args->output_type == SEQ_SER) {
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
		if (args->output_type == SEQ_SER) {
			if (ser_write_frame_from_fit(ser_file, &fit, frame + ser_frames)) {
				siril_log_message(_("Error while converting to SER (no space left?)\n"));
				retval = 1;
			}
		} else {
			// TODO: handle FITS sequence
			char dest_filename[128];
			g_snprintf(dest_filename, 128, "%s%05d", args->destroot, index++);
			if (savefits(dest_filename, &fit)) {
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

	wchar_t *wsource = g_utf8_to_utf16(source, -1, NULL, NULL, NULL);
	if ( wsource == NULL ) {
		return NULL ;
	}

	if (!(GetFileAttributesW(wsource) & FILE_ATTRIBUTE_REPARSE_POINT)) { /* Ce n'est pas un lien symbolique , je sors */
		g_free(wsource);
		return NULL;
	}

	wchar_t *wFilePath = g_new(wchar_t, maxchar + 1);
	if (!wFilePath) {
		PRINT_ALLOC_ERR;
		g_free(wsource);
		return NULL;
	}
	wFilePath[0] = 0;

	hFile = CreateFileW(wsource, GENERIC_READ, FILE_SHARE_READ, NULL,
			OPEN_EXISTING, 0, NULL);
	if (hFile == INVALID_HANDLE_VALUE) {
		g_free(wFilePath);
		g_free(wsource);
		return NULL;
	}

	GetFinalPathNameByHandleW(hFile, wFilePath, maxchar, 0);

	gchar *gFilePath = g_utf16_to_utf8(wFilePath + 4, -1, NULL, NULL, NULL); // +4 = enleve les 4 caracteres du prefixe "//?/"
	g_free(wsource);
	g_free(wFilePath);
	CloseHandle(hFile);
	return gFilePath;
}
#endif

/* open the file with path source from any image type and load it into the given FITS object */
int any_to_fits(image_type imagetype, const char *source, fits *dest,
		gboolean interactive, gboolean force_float, gboolean debayer) {
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
			retval = (readtif(source, dest, force_float) < 0);
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
			/* need to retrieve all return values */
			retval = readheif(source, dest, interactive);
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
				char *rsrc = g_real_path(source) ;
				if (rsrc != NULL) {
					src  = rsrc;					
				}
#endif
				retval = (open_raw_files(src , dest, debayer) < 0);
#ifdef _WIN32
				if (rsrc != NULL) {
					g_free(rsrc);
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

