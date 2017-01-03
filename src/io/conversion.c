/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2017 team free-astro (see more in AUTHORS file)
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

#include <gdk/gdkkeysyms.h>
#include <gtk/gtk.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <dirent.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "io/conversion.h"
#include "io/films.h"
#include "io/ser.h"
#include "gui/callbacks.h"
#include "algos/demosaicing.h"

#define MAX_OF_EXTENSIONS 50	// actual size of supported_extensions

static gchar *destroot = NULL;
static unsigned int convflags = CONV1X3;	// default
static unsigned int supported_filetypes = 0;	// initialized by initialize_converters()

// NULL-terminated array, initialized by initialize_converters(), used only by stat_file
char **supported_extensions;

supported_raw_list supported_raw[] = {
	{"dng",	"Adobe", BAYER_FILTER_RGGB},
	{"mos",	"Aptus", BAYER_FILTER_RGGB},
	{"cr2",	"Canon", BAYER_FILTER_RGGB},
	{"crw",	"Canon", BAYER_FILTER_RGGB},
	{"bay",	"Casio", BAYER_FILTER_NONE},		// Not tested
	{"erf",	"Epson", BAYER_FILTER_RGGB},
	{"raf",	"Fuji", BAYER_FILTER_GBRG},		// Not really supported, specially XTRANS
	{"3fr",	"Hasselblad", BAYER_FILTER_GRBG},	// GRBG, RGGB		
	{"kdc",	"Kodak", BAYER_FILTER_GRBG},
	{"dcr",	"Kodak", BAYER_FILTER_GRBG},
	{"mef",	"Mamiya", BAYER_FILTER_RGGB},
	{"mrw",	"Minolta", BAYER_FILTER_RGGB},
	{"nef",	"Nikon", BAYER_FILTER_RGGB},
	{"nrw",	"Nikon", BAYER_FILTER_RGGB},
	{"orf",	"Olympus", BAYER_FILTER_GRBG},
	{"raw",	"Leica", BAYER_FILTER_RGGB},
	{"rw2",	"Panasonic", BAYER_FILTER_BGGR},
	{"pef",	"Pentax", BAYER_FILTER_BGGR},
	{"ptx",	"Pentax", BAYER_FILTER_NONE},		// Not tested
	{"x3f",	"Sigma", BAYER_FILTER_NONE},		// Not supported yet
	{"srw",	"Samsung", BAYER_FILTER_BGGR},
	{"arw",	"Sony", BAYER_FILTER_RGGB}
};

char *filter_pattern[] = {
	"RGGB",
	"BGGR",
	"GBRG",
	"GRBG"
};

static gpointer convert_thread_worker(gpointer p);
static gboolean end_convert_idle(gpointer p);


int get_nb_raw_supported() {
	return sizeof(supported_raw) / sizeof(supported_raw_list);
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
		printf("*.%s",supported_raw[i].extension);
		if (i != nb_raw - 1) printf(", ");
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
}

void check_for_conversion_form_completeness() {
	static GtkTreeView *tree_convert = NULL;
	GtkTreeIter iter;
	GtkTreeModel *model = NULL;
	gboolean valid;
	GtkWidget *go_button = lookup_widget("convert_button");
	
	if (tree_convert == NULL)
		tree_convert = GTK_TREE_VIEW(gtk_builder_get_object(builder, "treeview_convert"));
	
	model = gtk_tree_view_get_model(tree_convert);
	valid = gtk_tree_model_get_iter_first(model, &iter);
	gtk_widget_set_sensitive (go_button, destroot && destroot[0] != '\0' && valid);
	update_statusbar_convert();
}

// input is destroot
char *create_sequence_filename(int counter, char *output, int outsize) {
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

// truncates destroot if it's more than 120 characters, append a '_' if it
// doesn't end with one or a '-'. SER extensions are accepted and unmodified.
void on_convtoroot_changed (GtkEditable *editable, gpointer user_data){
	static GtkWidget *multiple_ser = NULL;
	const gchar *name = gtk_entry_get_text(GTK_ENTRY(editable));
	if (!multiple_ser)
		multiple_ser = lookup_widget("multipleSER");
	if (destroot) g_free(destroot);
	destroot = g_strdup(name);

	const char *ext = get_filename_ext(destroot);
	if (ext && !g_ascii_strcasecmp(ext, "ser")) {
		convflags |= CONVDSTSER;
		gtk_widget_set_visible(multiple_ser, TRUE);
		return;
	}
	gtk_widget_set_visible(multiple_ser, FALSE);

	destroot = format_basename(destroot);

	check_for_conversion_form_completeness();
}

void on_demosaicing_toggled (GtkToggleButton *togglebutton, gpointer user_data) {
	static GtkToggleButton *but = NULL;
	if (!but) but = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton1"));
	if (gtk_toggle_button_get_active(togglebutton)) {
		convflags |= CONVDEBAYER;
		gtk_toggle_button_set_active(but, TRUE);
		com.debayer.open_debayer = TRUE;
	}
	else {
		convflags &= ~CONVDEBAYER;	// used for conversion
		com.debayer.open_debayer = FALSE;	// used for image opening
	}
}

void on_multipleSER_toggled (GtkToggleButton *togglebutton, gpointer user_data) {
	if (gtk_toggle_button_get_active(togglebutton))
		convflags |= CONVMULTIPLE;
	else convflags &= ~CONVMULTIPLE;
}

void on_conv3planefit_toggled (GtkToggleButton *togglebutton, gpointer user_data){
	convflags |= CONV1X3;
	convflags &= ~(CONV3X1|CONV1X1);
}

void on_conv3_1plane_toggled (GtkToggleButton *togglebutton, gpointer user_data) {
	convflags |= CONV3X1;
	convflags &= ~(CONV1X1|CONV1X3);
}

void on_conv1_1plane_toggled (GtkToggleButton *togglebutton, gpointer user_data) {
	convflags |= CONV1X1;
	convflags &= ~(CONV3X1|CONV1X3);
}

/*************************************************************************************/

/* This function sets all default values of libraw settings in the com.raw_set
 * struct, as defined in the glade file.
 * When the ini file is read, the values of com.raw_set are overwritten, but if the
 * file is missing, like the first time Siril is launched, we don't want to have the
 * GUI states reset to zero by set_GUI_LIBRAW() because the data in com.raw_set had
 * not been initialized with the default GUI values (= initialized to 0).
 */
void initialize_libraw_settings() {
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

void initialize_ser_debayer_settings() {
	com.debayer.open_debayer = FALSE;
	com.debayer.use_bayer_header = TRUE;
	com.debayer.compatibility = FALSE;
	com.debayer.bayer_pattern = BAYER_FILTER_RGGB;
	com.debayer.bayer_inter = BAYER_VNG;
}
 
/* initialize converters (utilities used for different image types importing) *
 * updates the label listing the supported input file formats, and modifies the
 * list of file types used in convflags */
void initialize_converters() {
	GtkLabel *label_supported;
	GString *string;
	gchar *text;
	int count_ext = 0;

	string = g_string_new("\t");
	/* internal converters */
	supported_filetypes |= TYPEBMP;
	g_string_append(string, _("BMP images, "));
	supported_filetypes |= TYPEPIC;
	g_string_append(string, _("PIC images (IRIS), "));
	supported_filetypes |= TYPEPNM;
	g_string_append(string, _("PGM and PPM binary images"));
		
	supported_extensions = malloc(MAX_OF_EXTENSIONS * sizeof(char*));
	/* internal extensions */
	if (supported_extensions == NULL) {
		fprintf(stderr, "initialize_converters: error allocating data\n");
		return;
	}
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
	g_string_append(string, ", ");
	g_string_append(string, _("RAW images"));
	set_libraw_settings_menu_available(TRUE);	// enable libraw settings
	initialize_libraw_settings();	// below in the file
	
	nb_raw = get_nb_raw_supported();
	for (i = 0; i < nb_raw; i++) {
		supported_extensions[count_ext+i] = malloc(strlen(supported_raw[i].extension) + 2 * sizeof (char));
		strcpy(supported_extensions[count_ext+i], ".");
		strcat(supported_extensions[count_ext+i], supported_raw[i].extension);
	}
	count_ext += nb_raw;
#else
	set_libraw_settings_menu_available(FALSE);	// disable libraw settings
#endif
	g_string_append(string, ", ");
	g_string_append(string, _("FITS-CFA images"));

#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
	supported_filetypes |= TYPEAVI;
	g_string_append(string, ", ");
	g_string_append(string, _("Films"));
#endif

	supported_filetypes |= TYPESER;
	g_string_append(string, ", ");
	g_string_append(string, _("SER sequences"));

	/* library converters (detected by configure) */
#ifdef HAVE_LIBTIFF
	supported_filetypes |= TYPETIFF;
	g_string_append(string, ", ");
	g_string_append(string, _("TIFF images"));
	supported_extensions[count_ext++] = ".tif";
	supported_extensions[count_ext++] = ".tiff";
#endif

#ifdef HAVE_LIBJPEG
	supported_filetypes |= TYPEJPG;
	g_string_append(string, ", ");
	g_string_append(string, _("JPG images"));
	supported_extensions[count_ext++] = ".jpg";
	supported_extensions[count_ext++] = ".jpeg";
#endif

#ifdef HAVE_LIBPNG
	supported_filetypes |= TYPEPNG;
	g_string_append(string, ", ");
	g_string_append(string, _("PNG images"));
	supported_extensions[count_ext++] = ".png";
#endif
	supported_extensions[count_ext++] = NULL;		// NULL-terminated array

	g_string_append(string, ".");
	label_supported = GTK_LABEL(gtk_builder_get_object(builder, "label_supported_types"));
	text = g_string_free(string, FALSE);
	gtk_label_set_text(label_supported, text);
	siril_log_message(_("Supported file types: %s\n"), text + 1);
	g_free(text);
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
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
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

int count_selected_files() {
	static GtkTreeView *tree_convert = NULL;
	GtkTreeModel *model = NULL;
	GtkTreeIter iter;
	gboolean valid;
	int count = 0;
	
	if (tree_convert == NULL)
		tree_convert = GTK_TREE_VIEW(gtk_builder_get_object(builder, "treeview_convert"));

	model = gtk_tree_view_get_model(tree_convert);
	valid = gtk_tree_model_get_iter_first(model, &iter);
	
	while (valid) {
		gchar *file_name, *file_date;
		gtk_tree_model_get (model, &iter, COLUMN_FILENAME, &file_name,
											COLUMN_DATE, &file_date,
											-1);
		valid = gtk_tree_model_iter_next (model, &iter);
		count ++;
	}
	return count;
}

struct _convert_data {
	struct timeval t_start;
	DIR *dir;
	GList *list;
	int start;
	int total;
	int nb_converted;
	gboolean compatibility;
};

void on_convert_button_clicked(GtkButton *button, gpointer user_data) {
	DIR *dir;
	gchar *file_data, *file_date;
	const gchar *indice;
	static GtkTreeView *tree_convert = NULL;
	static GtkEntry *startEntry = NULL;
	GtkTreeModel *model = NULL;
	GtkTreeIter iter;
	gboolean valid;
	GList *list = NULL;
	int count = 0;
	
	if (tree_convert == NULL) {
		tree_convert = GTK_TREE_VIEW(gtk_builder_get_object(builder, "treeview_convert"));
		startEntry = GTK_ENTRY(gtk_builder_get_object(builder, "startIndiceEntry"));
	}

	struct timeval t_start;
	
	if (get_thread_run()) {
		siril_log_message(_("Another task is already in progress, ignoring new request.\n"));
		return;
	}

	model = gtk_tree_view_get_model(tree_convert);
	valid = gtk_tree_model_get_iter_first(model, &iter);
	if (valid == FALSE) return;	//The tree is empty
	
	while (valid) {
		gtk_tree_model_get (model, &iter, COLUMN_FILENAME, &file_data,
				COLUMN_DATE, &file_date,
				-1);
		list = g_list_append (list, file_data);
		valid = gtk_tree_model_iter_next (model, &iter);
		count ++;
	}
	
	indice = gtk_entry_get_text(startEntry);

	siril_log_color_message(_("Conversion: processing...\n"), "red");
	gettimeofday(&t_start, NULL);
	
	set_cursor_waiting(TRUE);
	control_window_switch_to_tab(OUTPUT_LOGS);
	
	/* then, convert files to Siril's FITS format */
	struct _convert_data *args;
	set_cursor_waiting(TRUE);
	char *tmpmsg;
	if (!com.wd) {
		tmpmsg = siril_log_message(_("Conversion: no working directory set.\n"));
		show_dialog(tmpmsg, _("Warning"), "gtk-dialog-warning");
		set_cursor_waiting(FALSE);
		return;
	}
	if((dir = opendir(com.wd)) == NULL){
		tmpmsg = siril_log_message(_("Conversion: error opening working directory %s.\n"), com.wd);
		show_dialog(tmpmsg, _("Error"), "gtk-dialog-error");
		set_cursor_waiting(FALSE);
		return ;
	}

	args = malloc(sizeof(struct _convert_data));
	args->start = (atof(indice) == 0 || atof(indice) > USHRT_MAX) ? 1 : atof(indice);
	args->dir = dir;
	args->list = list;
	args->total = count;
	args->nb_converted = 0;
	args->t_start.tv_sec = t_start.tv_sec;
	args->t_start.tv_usec = t_start.tv_usec;
	args->compatibility = com.debayer.compatibility;
	start_in_new_thread(convert_thread_worker, args);
	return;
}

static gpointer convert_thread_worker(gpointer p) {
	char dest_filename[128], msg_bar[256];
	int indice;
	double progress = 0.0;
	struct ser_struct *ser_file = NULL;
	struct _convert_data *args = (struct _convert_data *) p;
	GList *list;
	
	list = g_list_first(args->list);
	indice = args->start;

	if (convflags & CONVDSTSER) {
		if (convflags & CONV3X1) {
			siril_log_color_message(_("SER output will take precedence over the one-channel per image creation option.\n"), "salmon");
			convflags &= ~CONV3X1;
		} else {
			ser_file = malloc(sizeof(struct ser_struct));
			if (!(convflags & CONVMULTIPLE)) {
				if (ser_create_file(destroot, ser_file, TRUE, NULL)) {
					siril_log_message(_("Creating the SER file failed, aborting.\n"));
					goto clean_exit;
				}
			}
		}
	}

	while (list) {
		gchar *src_filename = (gchar *)list->data;
		const char *src_ext = get_filename_ext(src_filename);
		image_type imagetype;

		if (!get_thread_run()) {
			break;
		}
		gchar *name = g_utf8_strrchr(src_filename, strlen(src_filename), '/');
		if (name)
			g_snprintf(msg_bar, 256, _("Converting %s..."), name + 1);
		else g_snprintf(msg_bar, 256, _("Converting %s..."), src_filename);

		imagetype = get_type_for_extension(src_ext);
		if (imagetype == TYPEUNDEF) {
			char msg[512];
			siril_log_message(_("FILETYPE IS NOT SUPPORTED, CANNOT CONVERT: %s\n"), src_ext);
			g_snprintf(msg, 512, _("File extension '%s' is not supported.\n"
				"Verify that you typed the extension correctly.\n"
				"If so, you may need to install third-party software to enable "
				"this file type conversion, look at the README file.\n"
				"If the file type you are trying to load is listed in supported "
				"formats, you may notify the developpers that the extension you are "
				"trying to use should be recognized for this type."), src_ext);
			show_dialog(msg, _("Error"), "gtk-dialog-error");
			break;	// avoid 100 error popups
		}

		if (imagetype == TYPEAVI) {
			// we need to do a semi-recursive thing here,
			// thankfully it's only one level deep
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
			int frame;
			fits *fit = calloc(1, sizeof(fits));
			struct film_struct film_file;
			if (film_open_file(src_filename, &film_file) != FILM_SUCCESS) {
				siril_log_message(_("Error while opening film %s, aborting.\n"), src_filename);
				break;
			}
			if (convflags & CONVMULTIPLE) {
				if (ser_create_file(create_sequence_filename(indice++, dest_filename, 128),
							ser_file, TRUE, NULL)) {
					siril_log_message(_("Creating the SER file failed, aborting.\n"));
					goto clean_exit;
				}
			}
			for (frame = 0; frame < film_file.frame_count; frame++) {
				// read frame from the film
				if (film_read_frame(&film_file, frame, fit) != FILM_SUCCESS) {
					siril_log_message(_("Error while reading frame %d from %s, aborting.\n"),
							frame, src_filename);
					goto clean_exit;
				}

				// save to the destination file
				if (convflags & CONVDSTSER) {
					if (convflags & CONV1X1)
						keep_first_channel_from_fits(fit);
					if (ser_write_frame_from_fit(ser_file, fit, frame)) {
						siril_log_message(_("Error while converting to SER (no space left?)\n"));
						goto clean_exit;
					}
				} else {
					g_snprintf(dest_filename, 128, "%s%05d", destroot, indice++);
					if (save_to_target_fits(fit, dest_filename)) {
						siril_log_message(_("Error while converting to FITS (no space left?)\n"));
						goto clean_exit;
					}
				}
				clearfits(fit);
			}
			if (convflags & CONVMULTIPLE) {
				ser_write_and_close(ser_file);
			}
			free(fit);
#endif
		}
		else if (imagetype == TYPESER) {
			siril_log_message(_("Converting from SER is not yet supported\n"));
			break;
		}
		else {	// single image
			fits *fit = any_to_new_fits(imagetype, src_filename, args->compatibility);
			if (convflags & CONVDSTSER) {
				if (convflags & CONV1X1)
					keep_first_channel_from_fits(fit);
				if (ser_write_frame_from_fit(ser_file, fit, args->nb_converted)) {
					siril_log_message(_("Error while converting to SER (no space left?)\n"));
					break;
				}
			} else {
				g_snprintf(dest_filename, 128, "%s%05d", destroot, indice++);
				if (save_to_target_fits(fit, dest_filename)) {
					siril_log_message(_("Error while converting to FITS (no space left?)\n"));
					break;
				}
			}
			clearfits(fit);
			free(fit);
		}

		set_progress_bar_data(msg_bar, progress/((double)args->total));
		progress += 1.0;
		args->nb_converted++;

		list = g_list_next(list);
	}

clean_exit:
	if (convflags & CONVDSTSER) {
		if (!(convflags & CONVMULTIPLE))
			ser_write_and_close(ser_file);
		free(ser_file);
	}

	gdk_threads_add_idle(end_convert_idle, args);
	return NULL;
}

static gboolean end_convert_idle(gpointer p) {
	struct _convert_data *args = (struct _convert_data *) p;
	struct timeval t_end;
	
	if (get_thread_run() && args->nb_converted > 1) {
		// load the sequence
		char *ppseqname = malloc(strlen(destroot) + 5);
		sprintf(ppseqname, "%s.seq", destroot);
		check_seq(0);
		update_sequences_list(ppseqname);
		free(ppseqname);
	}
	update_used_memory();
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_DONE);
	set_cursor_waiting(FALSE);
	gettimeofday(&t_end, NULL);
	show_time(args->t_start, t_end);
	stop_processing_thread();
	g_list_free_full(args->list, g_free);
	free(args);
	return FALSE;
}

/* from a fits object, save to file or files, based on the channel policy from convflags */
int save_to_target_fits(fits *fit, const char *dest_filename) {
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

static int retrieveBayerPattern(char *bayer) {
	int i;

	for (i = BAYER_FILTER_MIN; i < BAYER_FILTER_MAX; i++) {
		if (g_ascii_strcasecmp(bayer, filter_pattern[i]) == 0) {
			return i;
		}
	}
	return BAYER_FILTER_NONE;
}

int debayer_if_needed(image_type imagetype, fits *fit, gboolean compatibility) {
	int retval = 0;
	sensor_pattern tmp;
	/* What the hell?
	 * Siril's FITS are stored bottom to top, debayering will throw 
	 * wrong results. So before demosacaing we need to transforme the image
	 * with fits_flip_top_to_bottom() function */
	if (imagetype == TYPEFITS && (convflags & CONVDEBAYER)) {
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
		}
		if (com.debayer.bayer_pattern >= 0)
			siril_log_message(_("Filter Pattern: %s\n"), filter_pattern[com.debayer.bayer_pattern]);

		if (debayer(fit, com.debayer.bayer_inter)) {
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

/* open the file with path source from any image type and load it into a new FITS object */
fits *any_to_new_fits(image_type imagetype, const char *source, gboolean compatibility) {
	int retval = 0;
	fits *tmpfit = calloc(1, sizeof(fits));

	retval = any_to_fits(imagetype, source, tmpfit);

	if (!retval)
		retval = debayer_if_needed(imagetype, tmpfit, compatibility);

	if (retval) {
		clearfits(tmpfit);
		free(tmpfit);
		return NULL;
	}

	return tmpfit;
}

/* open the file with path source from any image type and load it into the given FITS object */
int any_to_fits(image_type imagetype, const char *source, fits *dest) {
	int retval = 0;

	switch (imagetype) {
		case TYPEFITS:
			retval = (readfits(source, dest, NULL) != 0);
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
#ifdef HAVE_LIBPNG
		case TYPEPNG:
			retval = (readpng(source, dest) < 0);
			break;
#endif
#ifdef HAVE_LIBRAW
		case TYPERAW:
			retval = (open_raw_files(source, dest, !(convflags & CONVDEBAYER)) < 0);
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

