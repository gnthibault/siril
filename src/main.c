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

#define MAIN
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <gtk/gtk.h>
#ifdef MAC_INTEGRATION
#include <gtkosxapplication.h>
#endif
#ifdef _WIN32
#include <windows.h>
#include <tchar.h>
#include <io.h>
#include <fcntl.h>
#endif
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <locale.h>
#if (defined(__APPLE__) && defined(__MACH__))
#include <stdlib.h>
#include <libproc.h>
#endif
#include <getopt.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/command.h"
#include "core/pipe.h"
#include "core/undo.h"
#include "core/signals.h"
#include "algos/star_finder.h"
#include "algos/photometry.h"
#include "io/sequence.h"
#include "io/conversion.h"
#include "io/single_image.h"
#include "gui/callbacks.h"
#include "gui/script_menu.h"
#include "gui/progress_and_log.h"
#include "registration/registration.h"
#include "stacking/stacking.h"

#define GLADE_FILE "siril3.glade"

/* the global variables of the whole project */
cominfo com;	// the main data struct
fits gfit;	// currently loaded image
GtkBuilder *builder = NULL;	// get widget references anywhere

#ifdef MAC_INTEGRATION
static gboolean osx_open_file(GtkosxApplication *osx_app, gchar *path,
		gpointer data) {
	if (path != NULL) {
		const char *ext = get_filename_ext(path);
		if (ext && !strncmp(ext, "seq", 3)) {
			gchar *sequence_dir = g_path_get_dirname(path);
			if (!changedir(sequence_dir, NULL)) {
				if (check_seq(FALSE)) {
					siril_log_message(_("No sequence `%s' found.\n"), path);
				} else {
					set_seq(path);
				}
				g_free(sequence_dir);
			}
		} else {
			open_single_image(path);
		}
		return FALSE;
	}
	return TRUE;
}

static void gui_add_osx_to_app_menu(GtkosxApplication *osx_app, const gchar *item_name, gint index) {
	GtkWidget *item;

	item = lookup_widget(item_name);
	if (GTK_IS_MENU_ITEM(item))
		gtkosx_application_insert_app_menu_item(osx_app, GTK_WIDGET(item), index);
}

static void set_osx_integration(GtkosxApplication *osx_app) {
	GtkWidget *menubar = lookup_widget("menubar1");
	GtkWidget *file_quit_menu_item = lookup_widget("exit");
	GtkWidget *help_menu = lookup_widget("help1");
	GtkWidget *window_menu = lookup_widget("menuitemWindows");
	GtkWidget *sep;
	GdkPixbuf *icon;
	gchar *icon_path;
	
	g_signal_connect(osx_app, "NSApplicationOpenFile", G_CALLBACK(osx_open_file), NULL);

	gtk_widget_hide(menubar);

	gtkosx_application_set_menu_bar(osx_app, GTK_MENU_SHELL(menubar));
	gtkosx_application_set_window_menu(osx_app, GTK_MENU_ITEM(window_menu));

	gui_add_osx_to_app_menu(osx_app, "help_item1", 0);
	gui_add_osx_to_app_menu(osx_app, "help_get_scripts", 1);
	gui_add_osx_to_app_menu(osx_app, "help_update", 2);
	sep = gtk_separator_menu_item_new();
	gtkosx_application_insert_app_menu_item(osx_app, sep, 2);
	gui_add_osx_to_app_menu(osx_app, "settings", 3);
	sep = gtk_separator_menu_item_new();
	gtkosx_application_insert_app_menu_item(osx_app, sep, 4);

	gtk_widget_hide(file_quit_menu_item);
	gtk_widget_hide(help_menu);
	
	icon_path = g_build_filename(com.app_path, "pixmaps", "siril.svg", NULL);
	icon = gdk_pixbuf_new_from_file(icon_path, NULL);
	gtkosx_application_set_dock_icon_pixbuf(osx_app, icon);
		
	gtkosx_application_ready(osx_app);
	g_free(icon_path);
}

#endif

#ifdef _WIN32
/* origine du source: https://stackoverflow.com/questions/24171017/win32-console-application-that-can-open-windows */
static int ReconnectIO(int OpenNewConsole) {
	int hConHandle;
	HANDLE lStdHandle;
	FILE *fp;
	int MadeConsole;

	MadeConsole = 0;
	if (!AttachConsole(ATTACH_PARENT_PROCESS)) {
		if (!OpenNewConsole)
			return 0;

		MadeConsole = 1;
		if (!AllocConsole())
			return 0;
	}

	// STDOUT to the console
	lStdHandle = GetStdHandle(STD_OUTPUT_HANDLE);
	hConHandle = _open_osfhandle((intptr_t) lStdHandle, _O_TEXT);
	fp = _fdopen(hConHandle, "w");
	*stdout = *fp;
	setvbuf( stdout, NULL, _IONBF, 0);

	// STDIN to the console
	lStdHandle = GetStdHandle(STD_INPUT_HANDLE);
	hConHandle = _open_osfhandle((intptr_t) lStdHandle, _O_TEXT);
	fp = _fdopen(hConHandle, "r");
	*stdin = *fp;
	setvbuf( stdin, NULL, _IONBF, 0);

	// STDERR to the console
	lStdHandle = GetStdHandle(STD_ERROR_HANDLE);
	hConHandle = _open_osfhandle((intptr_t) lStdHandle, _O_TEXT);
	fp = _fdopen(hConHandle, "w");
	*stderr = *fp;
	setvbuf(stderr, NULL, _IONBF, 0);

	return MadeConsole;
}
#endif

static char *siril_sources[] = {
#ifdef _WIN32
	"../share/siril/",
#elif (defined(__APPLE__) && defined(__MACH__))
	"/tmp/siril/Contents/Resources/share/siril/",
#endif
	PACKAGE_DATA_DIR"/",
	"/usr/share/siril/",
	"/usr/local/share/siril/",
	""
};

static gchar *main_option_directory = NULL;
static gchar *main_option_script = NULL;
static gchar *main_option_initfile = NULL;
static gboolean main_option_pipe = FALSE;

static gboolean _print_version_and_exit(const gchar *option_name,
		const gchar *value, gpointer data, GError **error) {
	g_print("%s %s\n", PACKAGE, VERSION);
	exit(EXIT_SUCCESS);
	return TRUE;
}

static gboolean _print_list_of_formats_and_exit(const gchar *option_name,
		const gchar *value, gpointer data, GError **error) {
	list_format_available();
	exit(EXIT_SUCCESS);
	return TRUE;
}

static GOptionEntry main_option[] = {
		{ "directory", 'd', 0, G_OPTION_ARG_FILENAME, &main_option_directory, N_("changing the current working directory as the argument"), NULL },
		{ "script", 's', 0, G_OPTION_ARG_FILENAME, &main_option_script, N_("run the siril commands script in console mode"), NULL },
		{ "initfile", 'i', 0, G_OPTION_ARG_FILENAME, &main_option_initfile, N_("load configuration from file name instead of the default configuration file"), NULL },
		{ "pipe", 'p', 0, G_OPTION_ARG_NONE, &main_option_pipe, N_("run in console mode with command and log stream through named pipes"), NULL },
		{ "format", 'f', G_OPTION_FLAG_NO_ARG, G_OPTION_ARG_CALLBACK, _print_list_of_formats_and_exit, N_("print all supported image file formats (depending on installed libraries)" ), NULL },
		{ "version", 'v', G_OPTION_FLAG_NO_ARG, G_OPTION_ARG_CALLBACK, _print_version_and_exit, N_("print the application’s version"), NULL},
		{ NULL },
};

static gchar *load_glade_file() {
	int i = 0;

	/* try to load the glade file, from the sources defined above */
	builder = gtk_builder_new();

	do {
		GError *err = NULL;
		gchar *gladefile;

		gladefile = g_build_filename(siril_sources[i], GLADE_FILE, NULL);
		if (gtk_builder_add_from_file(builder, gladefile, &err)) {
			fprintf(stdout, _("Successfully loaded '%s'\n"), gladefile);
			g_free(gladefile);
			break;
		}
		fprintf (stderr, _("%s. Looking into another directory...\n"), err->message);
		g_error_free(err);
		g_free(gladefile);
		i++;
	} while (i < G_N_ELEMENTS(siril_sources));
	if (i == G_N_ELEMENTS(siril_sources)) {
#ifdef _WIN32 // in the case where the app is started with double click on seq file
		gchar *execname = g_win32_get_package_installation_directory_of_module(NULL);
		gchar *prefix = g_build_filename(execname, "share", PACKAGE, NULL);
		g_free(execname);
		gchar *gladefile = g_build_filename(prefix, GLADE_FILE, NULL);
		if (gtk_builder_add_from_file(builder, gladefile, NULL)) {
			fprintf(stdout, _("Successfully loaded '%s'\n"), gladefile);
			g_free(gladefile);
			return prefix;
		}
		g_free(gladefile);
#endif
		fprintf(stderr, _("%s was not found or contains errors, cannot render GUI. Exiting.\n"), GLADE_FILE);
		exit(EXIT_FAILURE);
	}
	/* get back to the saved working directory */
	return g_strdup(siril_sources[i]);
}

static void global_initialization() {
	com.cvport = RED_VPORT;
	com.show_excluded = TRUE;
	com.selected_star = -1;
	com.star_is_seqdata = FALSE;
	com.stars = NULL;
	com.uniq = NULL;
	com.color = NORMAL_COLOR;
	int i;

	for (i = 0; i < MAXVPORT; i++)
		com.buf_is_dirty[i] = TRUE;
	memset(&com.selection, 0, sizeof(rectangle));
	memset(com.layers_hist, 0, sizeof(com.layers_hist));
	/* initialize the com struct and zoom level */
	com.sliders = MINMAX;
	com.zoom_value = ZOOM_DEFAULT;
	com.stack.mem_mode = 0;
	com.stack.memory_ratio = 0.9;
	com.stack.memory_amount = 4.0;
	com.app_path = NULL;
}

static void init_num_procs() {
	/* Get CPU number and set the number of threads */
#ifdef _OPENMP
	int num_proc = (int) g_get_num_processors();
	int omp_num_proc = omp_get_num_procs();
	if (num_proc != omp_num_proc) {
		siril_log_message(_("Questionable parallel processing efficiency - openmp reports %d processors. "
				"Possibly broken opencv/openblas installation.\n"), omp_num_proc);
	}
	siril_log_message(_("Parallel processing %s: Using %d logical processor(s).\n"), _("enabled"), com.max_thread = num_proc);
#else
	siril_log_message(_("Parallel processing %s: Using %d logical processor.\n"), _("disabled"), com.max_thread = 1);
#endif

	if (!com.headless) {
		update_spinCPU(com.max_thread);
	}
}

int main(int argc, char *argv[]) {
	GError *error = NULL;
	GOptionContext *ctx;
	gchar **args;
	gchar *startup_cwd = NULL;
	gboolean forcecwd = FALSE;
	gchar *cwd_forced = NULL;
	gchar *dir;

	g_setenv("LC_NUMERIC", "C", TRUE); // avoid possible bugs using french separator ","

	/* for translation */

	dir = get_siril_locale_dir();
	setlocale (LC_ALL, "");
	bindtextdomain(PACKAGE, dir);
	bind_textdomain_codeset(PACKAGE, "UTF-8");
	textdomain(PACKAGE);
	g_free(dir);

	memset(&com, 0, sizeof(struct cominf));	// needed?
	com.initfile = NULL;

	signals_init();

#ifdef _WIN32
	args = g_win32_get_command_line();
#else
	args = g_strdupv(argv);
#endif

	ctx = g_option_context_new(_("[FILE…]"));
	g_option_context_add_main_entries(ctx, main_option, PACKAGE);
	g_option_context_add_group(ctx, gtk_get_option_group(FALSE));

	if (!g_option_context_parse_strv(ctx, &args, &error)) {
		gchar *help_msg;

		help_msg = g_strdup_printf(_("Run “%s --help” to see a full "
				"list of available command line "
				"options."), args[0]);
		g_printerr("%s\n%s\n", error->message, help_msg);
		g_error_free(error);
		g_free(help_msg);
		g_option_context_free(ctx);
		g_strfreev(args);

		return 1;
	}
	g_option_context_free(ctx);

	if (main_option_script || main_option_pipe) {
		com.script = TRUE;
		com.headless = TRUE;
		/* need to force cwd to the current dir if no option -d */
		if (!forcecwd) {
			cwd_forced = g_strdup(g_get_current_dir());
			forcecwd = TRUE;
		}
	}

	if (main_option_initfile) {
		com.initfile = g_strdup(main_option_initfile);
	}

	if (main_option_directory) {
		if (!g_path_is_absolute(main_option_directory)) {
			cwd_forced = g_build_filename(g_get_current_dir(), main_option_directory, NULL);
		} else {
			cwd_forced = g_strdup(main_option_directory);
		}
		forcecwd = TRUE;
	}

	global_initialization();

	if (!com.headless) {
		gtk_init(&argc, &args);
	}

	siril_log_color_message(_("Welcome to %s v%s\n"), "bold", PACKAGE, VERSION);

	/***************
	 *  initialization of some parameters that need to be done before
	 * checkinitfile
	 ***************/
	/* initialize converters (utilities used for different image types importing) */
	gchar *supported_files = initialize_converters();
	/* initialize photometric variables */
	initialize_photometric_param();
	/* initialize peaker variables */
	init_peaker_default();
	/* initialize sequence-related stuff */
	initialize_sequence(&com.seq, TRUE);

	/* set default CWD, and load init file
	 * checkinitfile will load all saved parameters
	 * */
	com.wd = siril_get_startup_dir();
	startup_cwd = g_get_current_dir();
	if (checkinitfile()) {
		fprintf(stderr,	_("Could not load or create settings file, exiting.\n"));
		g_strfreev(args);
		exit(EXIT_FAILURE);
	}

	if (!com.headless) {
		/* load preferred theme */
		load_prefered_theme(com.combo_theme);
		/* Load glade file */
		gchar *path = load_glade_file();
		if (g_path_is_absolute(path)) {
			com.app_path = g_strdup(path);
		} else {
			com.app_path = g_build_filename(g_get_current_dir(), path, NULL);
		}
		g_free(path);
		/* load the css sheet for general style */
		load_css_style_sheet();
	}

	if (changedir(com.wd, NULL))
		com.wd = siril_get_startup_dir();

	if (!com.headless) {
		gtk_builder_connect_signals (builder, NULL);
		initialize_all_GUI(supported_files);
	}
	g_free(supported_files);

	init_num_procs();

	/* handling OS-X integration */
#ifdef MAC_INTEGRATION
	GtkosxApplication *osx_app = gtkosx_application_get();
	if (!com.headless) {
		set_osx_integration(osx_app);
	}
#endif //MAC_INTEGRATION

	/* open image, or sequence in argument, changing dir to be in its directory too */
	if (g_strv_length(args) > 1) {
		const char *ext = get_filename_ext(args[1]);
		if (ext && !strncmp(ext, "seq", 4)) {
			gchar *sequence_dir = g_path_get_dirname(args[1]);
			if (!changedir(sequence_dir, NULL)) {
				if (check_seq(FALSE)) {
					siril_log_message(_("No sequence `%s' found.\n"), args[1]);
				} else {
					set_seq(args[1]);
				}
				g_free(sequence_dir);
			}
		} else {
			if (startup_cwd) {
				changedir(startup_cwd, NULL);
			}
			open_single_image(args[1]);
			if (!forcecwd) {
				gchar *image_dir = g_path_get_dirname(args[1]);
				changedir(image_dir, NULL);
				g_free(image_dir);
			}
		}
	}
	g_free(startup_cwd);
	g_strfreev(args);

	if (forcecwd && cwd_forced) {
		changedir(cwd_forced, NULL);
		g_free(cwd_forced);
	}

	if (!com.script) {
		set_GUI_CWD();
	}

	if (com.headless) {
		if (main_option_script) {
			FILE* fp = g_fopen(main_option_script, "r");
			if (fp == NULL) {
				siril_log_message(_("File [%s] does not exist\n"), main_option_script);
				exit(EXIT_FAILURE);
			}
#ifdef _WIN32			
			ReconnectIO(1);
#endif
			if (execute_script(fp)) {
				exit(EXIT_FAILURE);
			}
		}
		else {
			pipe_start();
			read_pipe(NULL);
		}
	}
	else gtk_main();

	/* quit Siril */
	close_sequence(FALSE);	// closing a sequence if loaded
	close_single_image();	// close the previous image and free resources
	pipe_stop();		// close the pipes and their threads
	g_free(com.app_path);
#ifdef MAC_INTEGRATION
	g_object_unref(osx_app);
#endif //MAC_INTEGRATION
	return 0;
}
