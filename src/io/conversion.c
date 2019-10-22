/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2019 team free-astro (see more in AUTHORS file)
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
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "io/conversion.h"
#include "io/films.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "gui/callbacks.h"
#include "gui/message_dialog.h"
#include "gui/progress_and_log.h"
#include "algos/demosaicing.h"
#include "algos/sorting.h"

#define MAX_OF_EXTENSIONS 50	// actual size of supported_extensions

static gchar *destroot = NULL;
static unsigned int convflags = CONV1X3;	// default
static unsigned int supported_filetypes = 0;	// initialized by initialize_converters()

// NULL-terminated array, initialized by initialize_converters(), used only by stat_file
char **supported_extensions;

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

enum {
	COLUMN_CONV_FILENAME,		// gchar[]
	COLUMN_CONV_DATE,			// gchar[]
	N_COLUMNS
};

static gboolean end_convert_idle(gpointer p);


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

static gint sort_conv_tree(GtkTreeModel *model, GtkTreeIter *a, GtkTreeIter *b,
		gpointer user_data) {
	gchar *name_a, *name_b;
	gchar *collate_key1, *collate_key2;
	gint ret;

	gtk_tree_model_get(model, a, 0,	&name_a, -1);
	gtk_tree_model_get(model, b, 0,	&name_b, -1);

	collate_key1  = g_utf8_collate_key_for_filename(name_a, strlen(name_a));
	collate_key2  = g_utf8_collate_key_for_filename(name_b, strlen(name_b));

	ret = g_strcmp0(collate_key1, collate_key2);

	g_free(collate_key1);
	g_free(collate_key2);
	g_free(name_a);
	g_free(name_b);

	return ret;
}

// input is destroot
static char *create_sequence_filename(int counter, char *output, int outsize) {
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
	com.debayer.stretch = TRUE;
	com.debayer.compatibility = FALSE;
	com.debayer.bayer_pattern = BAYER_FILTER_RGGB;
	com.debayer.bayer_inter = BAYER_VNG;
}

static gboolean end_convert_idle(gpointer p) {
	struct _convert_data *args = (struct _convert_data *) p;
	struct timeval t_end;

	if (get_thread_run() && args->nb_converted > 1) {
		// load the sequence
		char *ppseqname = malloc(strlen(args->destroot) + 5);
		sprintf(ppseqname, "%s.seq", args->destroot);
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
static fits *any_to_new_fits(image_type imagetype, const char *source, gboolean compatibility, gboolean stretch_cfa) {
	int retval = 0;
	fits *tmpfit = calloc(1, sizeof(fits));

	retval = any_to_fits(imagetype, source, tmpfit);

	if (!retval)
		retval = debayer_if_needed(imagetype, tmpfit, compatibility, FALSE, stretch_cfa);

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

static void check_for_conversion_form_completeness() {
	static GtkTreeView *tree_convert = NULL;
	GtkTreeIter iter;
	GtkTreeModel *model = NULL;
	gboolean valid;
	GtkWidget *go_button = lookup_widget("convert_button");

	if (tree_convert == NULL)
		tree_convert = GTK_TREE_VIEW(lookup_widget("treeview_convert"));

	model = gtk_tree_view_get_model(tree_convert);
	valid = gtk_tree_model_get_iter_first(model, &iter);
	gtk_widget_set_sensitive (go_button, destroot && destroot[0] != '\0' && valid);

	/* we override the sort function in order to provide natural sort order */
	gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(model),
			COLUMN_CONV_FILENAME, (GtkTreeIterCompareFunc) sort_conv_tree, NULL,
			NULL);

	update_statusbar_convert();
}

static void unset_debayer_in_convflags() {
	convflags &= ~CONVDEBAYER;
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
	string = g_string_new("BMP images, ");
	supported_filetypes |= TYPEPIC;
	string = g_string_append(string, _("PIC images (IRIS), "));
	supported_filetypes |= TYPEPNM;
	string = g_string_append(string, _("PGM and PPM binary images"));
		
	supported_extensions = malloc(MAX_OF_EXTENSIONS * sizeof(char*));
	/* internal extensions */
	if (supported_extensions == NULL) {
		PRINT_ALLOC_ERR;
		return NULL;
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
	string = g_string_append(string, ", ");
	string = g_string_append(string, _("RAW images"));
	initialize_libraw_settings();	// below in the file
	
	nb_raw = get_nb_raw_supported();
	for (i = 0; i < nb_raw; i++) {
		supported_extensions[count_ext+i] = malloc(strlen(supported_raw[i].extension) + 2 * sizeof (char));
		strcpy(supported_extensions[count_ext+i], ".");
		strcat(supported_extensions[count_ext+i], supported_raw[i].extension);
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

int count_converted_files() {
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
		gtk_tree_model_get(model, &iter, COLUMN_FILENAME, &file_name,
				COLUMN_DATE, &file_date, -1);
		valid = gtk_tree_model_iter_next (model, &iter);
		count ++;
		g_free(file_name);
		g_free(file_date);
	}
	return count;
}

int count_selected_files() {
	GtkTreeView *tree_view = GTK_TREE_VIEW(lookup_widget("treeview_convert"));
	GtkTreeSelection *selection = gtk_tree_view_get_selection(tree_view);

	return gtk_tree_selection_count_selected_rows(selection);
}

static void initialize_convert() {
	GDir *dir;
	GError *error = NULL;
	gchar *file_data, *file_date;
	const gchar *indice;
	static GtkTreeView *tree_convert = NULL;
	static GtkEntry *startEntry = NULL;
	GtkTreeModel *model = NULL;
	GtkTreeIter iter;
	gboolean valid, several_type_of_files = FALSE;
	image_type imagetype = TYPEUNDEF;
	GList *list = NULL;
	int count = 0;
	
	if (tree_convert == NULL) {
		tree_convert = GTK_TREE_VIEW(lookup_widget("treeview_convert"));
		startEntry = GTK_ENTRY(lookup_widget("startIndiceEntry"));
	}

	struct timeval t_start;
	
	if (get_thread_run()) {
		siril_log_message(_("Another task is already in progress, ignoring new request.\n"));
		return;
	}

	/* test if forbidden chars exist */
	char *forbid_char = strchr(destroot, '/');
	if (forbid_char == NULL) {
		forbid_char = strchr(destroot, '\\');
	}
	if (forbid_char != NULL) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Invalid char"), _("Please remove invalid char in the sequence name "
				"before trying to convert images into a new sequence again."));
		return;
	}

	if (g_file_test(destroot, G_FILE_TEST_EXISTS)) {
		char *title = siril_log_message(_("A file named %s already exists. "
				"Do you want to replace it?\n"), destroot);
		gboolean replace = siril_confirm_dialog(title, _("The file already exists. "
				"Replacing it will overwrite its contents."), FALSE);
		if (!replace) return;
	}

	model = gtk_tree_view_get_model(tree_convert);
	valid = gtk_tree_model_get_iter_first(model, &iter);
	if (valid == FALSE) return;	//The tree is empty
	
	while (valid) {
		gtk_tree_model_get(model, &iter, COLUMN_FILENAME, &file_data,
				COLUMN_DATE, &file_date, -1);
		list = g_list_prepend(list, file_data);

		const char *src_ext = get_filename_ext(file_data);
		if (count != 0) {
			if (imagetype != get_type_for_extension(src_ext)) {
				several_type_of_files = TRUE;
			}
		}
		imagetype = get_type_for_extension(src_ext);
		valid = gtk_tree_model_iter_next(model, &iter);
		count++;
	}

	if ((convflags & CONVDEBAYER) && imagetype == TYPESER && !several_type_of_files) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("A conflict has been detected."),
				_("The Debayer option is not allowed in SER conversion, please uncheck the option."));
		g_list_free_full(list, g_free);
		set_cursor_waiting(FALSE);
		return;
	}
	if ((convflags & CONVMULTIPLE) && imagetype == TYPESER && !several_type_of_files) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("A conflict has been detected."),
				_("The Multiple SER option is not allowed in SER conversion, please uncheck the option."));
		g_list_free_full(list, g_free);
		set_cursor_waiting(FALSE);
		return;
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
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Warning"), tmpmsg);
		g_list_free_full(list, g_free);
		set_cursor_waiting(FALSE);
		return;
	}
	if((dir = g_dir_open(com.wd, 0, &error)) == NULL){
		tmpmsg = siril_log_message(_("Conversion: error opening working directory %s.\n"), com.wd);
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), tmpmsg);
		fprintf (stderr, "Conversion: %s\n", error->message);
		g_error_free(error);
		g_list_free_full(list, g_free);
		set_cursor_waiting(FALSE);
		return ;
	}

	/* g_list_append() has to traverse the entire list to find the end,
	 * which is inefficient when adding multiple elements. A common idiom
	 * to avoid the inefficiency is to use g_list_prepend() and reverse the
	 * list with g_list_reverse() when all elements have been added. */
	list = g_list_reverse(list);

	args = malloc(sizeof(struct _convert_data));
	args->start = (atof(indice) == 0 || atof(indice) > USHRT_MAX) ? 1 : atof(indice);
	args->dir = dir;
	args->list = list;
	args->total = count;
	args->nb_converted = 0;
	args->t_start.tv_sec = t_start.tv_sec;
	args->t_start.tv_usec = t_start.tv_usec;
	args->compatibility = com.debayer.compatibility;
	args->stretch_cfa = com.debayer.stretch;
	args->command_line = FALSE;
	args->several_type_of_files = several_type_of_files;
	args->destroot = g_strdup(destroot);
	start_in_new_thread(convert_thread_worker, args);
	return;
}

void on_entry2_activate(GtkEntry *entry, gpointer user_data) {
	initialize_convert();
}


void on_convert_button_clicked(GtkButton *button, gpointer user_data) {
	initialize_convert();
}

gpointer convert_thread_worker(gpointer p) {
	char dest_filename[128], msg_bar[256];
	int indice;
	int ser_frames = 0;
	double progress = 0.0;
	struct ser_struct *ser_file = NULL;
	struct _convert_data *args = (struct _convert_data *) p;
	GList *list;
	
	indice = args->start;

	if (convflags & CONVDSTSER) {
		if (convflags & CONV3X1) {
			siril_log_color_message(_("SER output will take precedence over the one-channel per image creation option.\n"), "salmon");
			convflags &= ~CONV3X1;
		}

		ser_file = malloc(sizeof(struct ser_struct));
		if (!(convflags & CONVMULTIPLE)) {
			if (ser_create_file(args->destroot, ser_file, TRUE, NULL)) {
				siril_log_message(_("Creating the SER file failed, aborting.\n"));
				goto clean_exit;
			}
		}
	}

	for (list = args->list; list; list = list->next) {
		gchar *src_filename = (gchar *)list->data;
		const char *src_ext = get_filename_ext(src_filename);
		image_type imagetype;

		if (!get_thread_run()) {
			break;
		}
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
			break;	// avoid 100 error popups
		}

		if (imagetype == TYPEAVI) {
			// we need to do a semi-recursive thing here,
			// thankfully it's only one level deep
#ifdef HAVE_FFMS2
			int frame;
			fits *fit = calloc(1, sizeof(fits));
			struct film_struct film_file;
			if (film_open_file(src_filename, &film_file) != FILM_SUCCESS) {
				siril_log_message(_("Error while opening film %s, aborting.\n"), src_filename);
				clearfits(fit);
				free(fit);
				break;
			}
			if (convflags & CONVMULTIPLE) {
				if (ser_create_file(create_sequence_filename(indice++, dest_filename, 128),
							ser_file, TRUE, NULL)) {
					siril_log_message(_("Creating the SER file failed, aborting.\n"));
					clearfits(fit);
					free(fit);
					goto clean_exit;
				}
			}
			for (frame = 0; frame < film_file.frame_count; frame++) {
				if (!get_thread_run()) {
					break;
				}
				// read frame from the film
				if (film_read_frame(&film_file, frame, fit) != FILM_SUCCESS) {
					siril_log_message(_("Error while reading frame %d from %s, aborting.\n"),
							frame, src_filename);
					clearfits(fit);
					free(fit);
					goto clean_exit;
				}

				// save to the destination file
				if (convflags & CONVDSTSER) {
					if (convflags & CONV1X1)
						keep_first_channel_from_fits(fit);
					if (ser_write_frame_from_fit(ser_file, fit, frame)) {
						siril_log_message(_("Error while converting to SER (no space left?)\n"));
						clearfits(fit);
						free(fit);
						goto clean_exit;
					}
				} else {
					g_snprintf(dest_filename, 128, "%s%05d", args->destroot, indice++);
					if (save_to_target_fits(fit, dest_filename)) {
						siril_log_message(_("Error while converting to FITS (no space left?)\n"));
						clearfits(fit);
						free(fit);
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
			if (args->several_type_of_files) {
				siril_log_message(_("Joining SER files is only possible with a list "
						"only containing SER files. Please, remove non SER files.\n"));
				break;
			}
			int frame;
			fits *fit = calloc(1, sizeof(fits));
			struct ser_struct tmp_ser;
			ser_init_struct(&tmp_ser);
			if (ser_open_file(src_filename, &tmp_ser)) {
				siril_log_message(_("Error while opening ser file %s, aborting.\n"), src_filename);
				clearfits(fit);
				free(fit);
				break;
			}
			if (args->nb_converted > 0 && (convflags & CONVDSTSER)) {
				if (tmp_ser.image_height != ser_file->image_height
						|| tmp_ser.image_width != ser_file->image_width) {
					siril_log_color_message(_("Input SER files must have the same size to be joined.\n"), "red");
					clearfits(fit);
					free(fit);
					break;
				}
			}
			set_progress_bar_data(msg_bar, PROGRESS_PULSATE);
			for (frame = 0; frame < tmp_ser.frame_count; frame++) {
				if (!get_thread_run()) {
					break;
				}
				// read frame from the film
				if (ser_read_frame(&tmp_ser, frame, fit)) {
					siril_log_message(_("Error while reading frame %d from %s, aborting.\n"),
							frame, src_filename);
					clearfits(fit);
					free(fit);
					goto clean_exit;
				}

				// save to the destination file
				if (convflags & CONVDSTSER) {
					if (convflags & CONV1X1)
						keep_first_channel_from_fits(fit);
					if (ser_write_frame_from_fit(ser_file, fit, frame + ser_frames)) {
						siril_log_message(_("Error while converting to SER (no space left?)\n"));
						clearfits(fit);
						free(fit);
						goto clean_exit;
					}
				} else {
					g_snprintf(dest_filename, 128, "%s%05d", args->destroot, indice++);
					if (save_to_target_fits(fit, dest_filename)) {
						siril_log_message(_("Error while converting to FITS (no space left?)\n"));
						clearfits(fit);
						free(fit);
						goto clean_exit;
					}
				}
				clearfits(fit);
			}
			ser_frames += frame;
			ser_close_file(&tmp_ser);
			free(fit);
		}
		else {	// single image
			fits *fit = any_to_new_fits(imagetype, src_filename, args->compatibility, args->stretch_cfa);
			if (fit) {
				if (convflags & CONVDSTSER) {
					if (convflags & CONV1X1)
						keep_first_channel_from_fits(fit);
					if (ser_write_frame_from_fit(ser_file, fit, args->nb_converted)) {
						siril_log_message(_("Error while converting to SER (no space left?)\n"));
						break;
					}
				} else {
					g_snprintf(dest_filename, 128, "%s%05d", args->destroot, indice++);
					if (save_to_target_fits(fit, dest_filename)) {
						siril_log_message(_("Error while converting to FITS (no space left?)\n"));
						break;
					}
				}
				clearfits(fit);
				free(fit);
			}
		}

		set_progress_bar_data(msg_bar, progress/((double)args->total));
		progress += 1.0;
		args->nb_converted++;
	}

clean_exit:
	if (convflags & CONVDSTSER) {
		if (!(convflags & CONVMULTIPLE))
			ser_write_and_close(ser_file);
		free(ser_file);
	}
	if (args->command_line) {
		unset_debayer_in_convflags();
	}
	g_list_free_full(args->list, g_free);
	g_dir_close(args->dir);
	siril_add_idle(end_convert_idle, args);
	return NULL;
}

int debayer_if_needed(image_type imagetype, fits *fit, gboolean compatibility, gboolean force_debayer, gboolean stretch_cfa) {
	int retval = 0;
	sensor_pattern tmp;
	/* What the hell?
	 * Siril's FITS are stored bottom to top, debayering will throw 
	 * wrong results. So before demosacaing we need to transforme the image
	 * with fits_flip_top_to_bottom() function */
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
			} else { /* FIXME: XTRANS CASE. TESTED FOR ONE FILE */
				com.debayer.bayer_pattern = XTRANS_FILTER;
				com.debayer.bayer_inter = XTRANS;
				siril_log_color_message(_("XTRANS Sensor detected. Using special algorithm.\n"), "red");
			}
		}
		if (com.debayer.bayer_pattern >= BAYER_FILTER_MIN
				&& com.debayer.bayer_pattern <= BAYER_FILTER_MAX) {
			siril_log_message(_("Filter Pattern: %s\n"), filter_pattern[com.debayer.bayer_pattern]);
		}

		if (stretch_cfa && fit->maximum_pixel_value) {
			siril_log_message(_("The FITS file is being normalized to 16-bit\n"));
		}

		if (debayer(fit, com.debayer.bayer_inter, stretch_cfa)) {
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

/**************** Conversion tree managment ***********************************/
/*
 * Main conversion list static functions
 */

static GtkListStore *liststore_convert = NULL;

static void add_convert_to_list(char *filename, GStatBuf st) {
	GtkTreeIter iter;
	char *date;

	date = ctime(&st.st_mtime);
	date[strlen(date) - 1] = 0;	// removing '\n' at the end of the string

	gtk_list_store_append(liststore_convert, &iter);
	gtk_list_store_set(liststore_convert, &iter, COLUMN_FILENAME, filename,	// copied in the store
			COLUMN_DATE, date, -1);
}
static void get_convert_list_store() {
	if (liststore_convert == NULL)
		liststore_convert = GTK_LIST_STORE(
				gtk_builder_get_object(builder, "liststore_convert"));
}

static GList *get_row_references_of_selected_rows(GtkTreeSelection *selection,
		GtkTreeModel *model) {
	GList *ref = NULL;
	GList *sel, *s;

	sel = gtk_tree_selection_get_selected_rows(selection, &model);

	for (s = sel; s; s = s->next) {
		GtkTreeRowReference *rowref = gtk_tree_row_reference_new(model,	(GtkTreePath *) s->data);
		ref = g_list_prepend(ref, rowref);
	}
	g_list_free_full(sel, (GDestroyNotify) gtk_tree_path_free);
	return ref;
}

static void remove_selected_files_from_list() {
	GtkTreeSelection *selection;
	GtkTreeModel *model;
	GList *references, *list;
	GtkTreeView *tree_view;

	tree_view = GTK_TREE_VIEW(lookup_widget("treeview_convert"));
	model = gtk_tree_view_get_model(tree_view);
	selection = gtk_tree_view_get_selection(tree_view);
	references = get_row_references_of_selected_rows(selection, model);
	for (list = references; list; list = list->next) {
		GtkTreeIter iter;
		GtkTreePath *path = gtk_tree_row_reference_get_path((GtkTreeRowReference*)list->data);
		if (path) {
			if (gtk_tree_model_get_iter(model, &iter, path)) {
				gtk_list_store_remove(liststore_convert, &iter);
			}
			gtk_tree_path_free(path);
		}
	}
	g_list_free(references);
	gtk_tree_selection_unselect_all(selection);
}

void fill_convert_list(GSList *list) {
	GStatBuf st;
	GSList *l;

	get_convert_list_store();

	for (l = list; l; l = l->next) {
		char *filename;

		filename = (char *) l->data;
		if (g_stat(filename, &st) == 0) {
			add_convert_to_list(filename, st);
		}
		g_free(filename);
	}
	check_for_conversion_form_completeness();
}

void on_clear_convert_button_clicked(GtkToolButton *button, gpointer user_data) {
	get_convert_list_store();
	gtk_list_store_clear(liststore_convert);
	check_for_conversion_form_completeness();
}

void on_remove_convert_button_clicked(GtkToolButton *button, gpointer user_data) {
	remove_selected_files_from_list();
	check_for_conversion_form_completeness();
}

void on_treeview_convert_drag_data_received(GtkWidget *widget,
		GdkDragContext *context, gint x, gint y,
		GtkSelectionData *selection_data, guint info, guint time,
		gpointer user_data) {

	gchar **uris, **str;
	const guchar *data;
	GSList *list = NULL;
	gint bad_files = 0;

	if (info != 0)
		return;

	data = gtk_selection_data_get_data(selection_data);
	uris = g_uri_list_extract_uris((gchar *) data);

	for (str = uris; *str; str++) {
		GError *error = NULL;
		gchar *path = g_filename_from_uri(*str, NULL, &error);
		if (path) {
			const char *src_ext = get_filename_ext(path);
			if (src_ext) {
				if (get_type_for_extension(src_ext) == TYPEUNDEF) {
					bad_files++;
				} else {
					list = g_slist_prepend(list, path);
				}
			} else bad_files++;
		} else {
			fprintf(stderr, "Could not convert uri to local path: %s",
					error->message);
			bad_files++;
			g_error_free(error);
		}
	}
	list = g_slist_sort(list, (GCompareFunc) strcompare);
	fill_convert_list(list);
	if (bad_files) {
		char *msg = siril_log_message(_("%d file(s) were ignored while drag and drop\n"), bad_files);
		siril_message_dialog(GTK_MESSAGE_INFO, msg,
				_("Files with unknown extension cannot be dropped in this area. "
						"Therefore they are ignored."));
	}
	g_strfreev(uris);
	g_slist_free(list);
}

gboolean on_treeview_convert_key_release_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {
	if (event->keyval == GDK_KEY_Delete || event->keyval == GDK_KEY_KP_Delete
			|| event->keyval == GDK_KEY_BackSpace) {
		remove_selected_files_from_list();
		check_for_conversion_form_completeness();
		return TRUE;
	}
	return FALSE;
}


/******************Callback functions*******************************************************************/

static gchar forbidden_char[] = { '/', '\\' };

static gboolean is_forbiden(gchar c) {
	int i;

	for (i = 0; i < G_N_ELEMENTS(forbidden_char); i++) {
		if (c == forbidden_char[i]) {
			return TRUE;
		}
	}
	return FALSE;
}

void insert_text_handler(GtkEntry *entry, const gchar *text, gint length,
		gint *position, gpointer data) {
	GtkEditable *editable = GTK_EDITABLE(entry);
	int i, count = 0;

	gchar *result = g_strndup(text, length);

	for (i = 0; i < length; i++) {
		if (is_forbiden(text[i]))
			continue;
		result[count++] = text[i];
	}

	if (count > 0) {
		g_signal_handlers_block_by_func(G_OBJECT (editable),
				G_CALLBACK (insert_text_handler), data);
		gtk_editable_insert_text(editable, result, count, position);
		g_signal_handlers_unblock_by_func(G_OBJECT (editable),
				G_CALLBACK (insert_text_handler), data);
	}
	g_signal_stop_emission_by_name(G_OBJECT(editable), "insert_text");

	g_free(result);
}

void on_treeview_selection5_changed(GtkTreeSelection *treeselection,
		gpointer user_data) {
	update_statusbar_convert();
}

// truncates destroot if it's more than 120 characters, append a '_' if it
// doesn't end with one or a '-'. SER extensions are accepted and unmodified.
void on_convtoroot_changed(GtkEditable *editable, gpointer user_data){
	static GtkWidget *multiple_ser = NULL;
	const gchar *name = gtk_entry_get_text(GTK_ENTRY(editable));

	if (*name != 0) {
		if (!multiple_ser)
			multiple_ser = lookup_widget("multipleSER");
		if (destroot)
			g_free(destroot);

		destroot = g_str_to_ascii(name, NULL); // we want to avoid special char

		const char *ext = get_filename_ext(destroot);
		if (ext && !g_ascii_strcasecmp(ext, "ser")) {
			convflags |= CONVDSTSER;
			gtk_widget_set_visible(multiple_ser, TRUE);
			char *base = remove_ext_from_filename(destroot);
			if (check_if_seq_exist(base)) {
				set_icon_entry(GTK_ENTRY(editable), "gtk-dialog-warning");
			} else{
				set_icon_entry(GTK_ENTRY(editable), NULL);
			}
			free(base);
		} else {
			convflags &= ~CONVDSTSER;
			gtk_widget_set_visible(multiple_ser, FALSE);
			destroot = format_basename(destroot);
			if (check_if_seq_exist(destroot)) {
				set_icon_entry(GTK_ENTRY(editable), "gtk-dialog-warning");
			} else{
				set_icon_entry(GTK_ENTRY(editable), NULL);
			}
		}

	} else {
		set_icon_entry(GTK_ENTRY(editable), NULL);
		g_free(destroot);
		destroot = NULL;
	}
	check_for_conversion_form_completeness();
}

void on_demosaicing_toggled (GtkToggleButton *togglebutton, gpointer user_data) {
	static GtkToggleButton *but = NULL;
	if (!but) but = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton1"));
	if (gtk_toggle_button_get_active(togglebutton)) {
		set_debayer_in_convflags();
		gtk_toggle_button_set_active(but, TRUE);
		com.debayer.open_debayer = TRUE;
	}
	else {
		unset_debayer_in_convflags();	// used for conversion
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
