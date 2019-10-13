/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2015 team free-astro (see more in AUTHORS file)
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

#include <errno.h>
#include <libconfig.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <gio/gio.h>
#include <glib/gstdio.h>
#ifdef _WIN32
#define DATADIR datadir
/* Constant available since Shell32.dll 4.72 */
#ifndef CSIDL_APPDATA
#define CSIDL_APPDATA 0x001a
#endif
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"

#define CFG_FILE "siril.cfg"

static const char *keywords[] = { "working-directory", "libraw-settings",
		"debayer-settings", "prepro-settings", "registration-settings",
		"stacking-settings", "photometry-settings", "misc-settings" };

static int readinitfile() {
	config_t config;
	const char *dir = NULL, *swap_dir = NULL, *extension = NULL;
	GSList *list = NULL;

	if (!com.initfile)
		return 1;

	config_init(&config);

	if (config_read_file(&config, com.initfile) == CONFIG_FALSE)
		return 1;
	siril_log_message(_("Loading init file: '%s'\n"), com.initfile);

	/* Keeping the up-scaled files poses a few problems with sequence
	 * filtering changing and user comprehension, so for now it can only be
	 * enabled by uncommenting the following line. */
	//com.cache_upscaled = TRUE;

	/* Working directory */
	if (config_lookup_string(&config, keywords[WD], &dir)) {
		free(com.wd);
		com.wd = g_strdup(dir);
	}

	/* Libraw setting */
	config_setting_t *raw_setting = config_lookup(&config, keywords[RAW]);
	if (raw_setting) {
		config_setting_lookup_float(raw_setting, "mul_0", &com.raw_set.mul[0]);
		config_setting_lookup_float(raw_setting, "mul_2", &com.raw_set.mul[2]);
		config_setting_lookup_float(raw_setting, "bright", &com.raw_set.bright);
		config_setting_lookup_int(raw_setting, "auto", &com.raw_set.auto_mul);
		config_setting_lookup_int(raw_setting, "cam_wb",
				&com.raw_set.use_camera_wb);
		config_setting_lookup_int(raw_setting, "auto_wb",
				&com.raw_set.use_auto_wb);
		config_setting_lookup_int(raw_setting, "user_qual",
				&com.raw_set.user_qual);
		config_setting_lookup_float(raw_setting, "gamm_0",
				&com.raw_set.gamm[0]);
		config_setting_lookup_float(raw_setting, "gamm_1",
				&com.raw_set.gamm[1]);
		config_setting_lookup_int(raw_setting, "user_black",
				&com.raw_set.user_black);
	}

	/* Debayer setting */
	config_setting_t *debayer_setting = config_lookup(&config, keywords[BAY]);
	if (debayer_setting) {
		config_setting_lookup_bool(debayer_setting, "ser_use_bayer_header",
				&com.debayer.use_bayer_header);
		config_setting_lookup_int(debayer_setting, "pattern",
				&com.debayer.bayer_pattern);
		config_setting_lookup_bool(debayer_setting, "compatibility",
						&com.debayer.compatibility);
		int inter;
		config_setting_lookup_int(debayer_setting, "inter", &inter);
		com.debayer.bayer_inter = inter;
		config_setting_lookup_bool(debayer_setting, "stretch",
						&com.debayer.stretch);
	}

	/* Preprocessing settings */
	config_setting_t *prepro_setting = config_lookup(&config, keywords[PRE]);
	if (prepro_setting) {
		config_setting_lookup_bool(prepro_setting, "cfa", &com.prepro_cfa);
		config_setting_lookup_bool(prepro_setting, "equalize_cfa", &com.prepro_equalize_cfa);
	}

	/* Registration setting */
	config_setting_t *reg_setting = config_lookup(&config, keywords[REG]);
	if (reg_setting) {
		config_setting_lookup_int(reg_setting, "method", &com.reg_settings);
	}

	/* Stacking setting */
	config_setting_t *stack_setting = config_lookup(&config, keywords[STK]);
	if (stack_setting) {
		config_setting_lookup_int(stack_setting, "method", &com.stack.method);
		config_setting_lookup_int(stack_setting, "rejection",
				&com.stack.rej_method);
		config_setting_lookup_int(stack_setting, "normalisation",
				&com.stack.normalisation_method);

		int mode = 0;
		config_setting_lookup_int(stack_setting, "mem_mode",
				&mode);
		com.stack.mem_mode = mode;
		config_setting_lookup_float(stack_setting, "maxmem",
				&com.stack.memory_ratio);
		config_setting_lookup_float(stack_setting, "maxmem_gb",
				&com.stack.memory_amount);
	}
	if (com.stack.mem_mode < 0 || com.stack.mem_mode > 2)
		com.stack.mem_mode = RATIO;
	if (com.stack.memory_ratio <= 0.05)
		com.stack.memory_ratio = 0.9;

	/* Photometry setting */
	config_setting_t *photometry_setting = config_lookup(&config, keywords[PTM]);
	if (photometry_setting) {
		config_setting_lookup_float(photometry_setting, "gain", &com.phot_set.gain);
		config_setting_lookup_float(photometry_setting, "inner-radius", &com.phot_set.inner);
		config_setting_lookup_float(photometry_setting, "outer-radius", &com.phot_set.outer);
		config_setting_lookup_int(photometry_setting, "minval", &com.phot_set.minval);
		config_setting_lookup_int(photometry_setting, "maxval", &com.phot_set.maxval);
	}

	/* Misc setting */
	config_setting_t *misc_setting = config_lookup(&config, keywords[MISC]);
	if (misc_setting) {
		config_setting_lookup_bool(misc_setting, "confirm", &com.dontShowConfirm);
		config_setting_lookup_int(misc_setting, "theme", &com.combo_theme);
		config_setting_lookup_bool(misc_setting, "remember_winpos", &com.remember_windows);
		config_setting_lookup_string(misc_setting, "swap_directory", &swap_dir);
		config_setting_lookup_string(misc_setting, "extension", &extension);

		misc_setting = config_lookup(&config, "misc-settings.scripts_paths");
		if (misc_setting != NULL) {
			unsigned int count = config_setting_length(misc_setting);
			unsigned int i;
			const char *tmp = NULL;

			for (i = 0; i < count; ++i) {
				tmp = config_setting_get_string_elem(misc_setting, i);
				list = g_slist_append(list, g_strdup(tmp));
			}
		}
		misc_setting = config_lookup(&config, "misc-settings.main_w_pos");
		if (misc_setting != NULL) {
			com.main_w_pos.x = config_setting_get_int_elem(misc_setting, 0);
			com.main_w_pos.y = config_setting_get_int_elem(misc_setting, 1);
			com.main_w_pos.w = config_setting_get_int_elem(misc_setting, 2);
			com.main_w_pos.h = config_setting_get_int_elem(misc_setting, 3);
		}

		misc_setting = config_lookup(&config, "misc-settings.rgb_w_pos");
		if (misc_setting != NULL) {
			com.rgb_w_pos.x = config_setting_get_int_elem(misc_setting, 0);
			com.rgb_w_pos.y = config_setting_get_int_elem(misc_setting, 1);
			com.rgb_w_pos.w = config_setting_get_int_elem(misc_setting, 2);
			com.rgb_w_pos.h = config_setting_get_int_elem(misc_setting, 3);
		}

	}
	if (swap_dir && swap_dir[0] != '\0') {
		if (com.swap_dir)
			g_free(com.swap_dir);
		com.swap_dir = g_strdup(swap_dir);
	} else {
		const char* sw_dir = g_get_tmp_dir();
		com.swap_dir = g_strdup(sw_dir);
	}
	if (extension && extension[0] != '\0') {
		if (com.ext)
			free(com.ext);
		com.ext = strdup(extension);
	} else {
		com.ext = strdup(".fit");
	}
	com.script_path = list;
	config_destroy(&config);
	return 0;
}

static void _save_wd(config_t *config, config_setting_t *root) {
	config_setting_t *directory;

	directory = config_setting_add(root, keywords[WD], CONFIG_TYPE_STRING);
	config_setting_set_string(directory, com.wd);
}

static void _save_libraw(config_t *config, config_setting_t *root) {
	config_setting_t *libraw_group, *raw_setting;

	libraw_group = config_setting_add(root, keywords[RAW], CONFIG_TYPE_GROUP);

	raw_setting = config_setting_add(libraw_group, "mul_0", CONFIG_TYPE_FLOAT);
	config_setting_set_float(raw_setting, com.raw_set.mul[0]);

	raw_setting = config_setting_add(libraw_group, "mul_2", CONFIG_TYPE_FLOAT);
	config_setting_set_float(raw_setting, com.raw_set.mul[2]);

	raw_setting = config_setting_add(libraw_group, "bright", CONFIG_TYPE_FLOAT);
	config_setting_set_float(raw_setting, com.raw_set.bright);

	raw_setting = config_setting_add(libraw_group, "auto", CONFIG_TYPE_INT);
	config_setting_set_int(raw_setting, com.raw_set.auto_mul);

	raw_setting = config_setting_add(libraw_group, "cam_wb", CONFIG_TYPE_INT);
	config_setting_set_int(raw_setting, com.raw_set.use_camera_wb);

	raw_setting = config_setting_add(libraw_group, "auto_wb", CONFIG_TYPE_INT);
	config_setting_set_int(raw_setting, com.raw_set.use_auto_wb);

	raw_setting = config_setting_add(libraw_group, "user_qual",
	CONFIG_TYPE_INT);
	config_setting_set_int(raw_setting, com.raw_set.user_qual);

	raw_setting = config_setting_add(libraw_group, "gamm_0", CONFIG_TYPE_FLOAT);
	config_setting_set_float(raw_setting, com.raw_set.gamm[0]);

	raw_setting = config_setting_add(libraw_group, "gamm_1", CONFIG_TYPE_FLOAT);
	config_setting_set_float(raw_setting, com.raw_set.gamm[1]);

	raw_setting = config_setting_add(libraw_group, "user_black",
	CONFIG_TYPE_INT);
	config_setting_set_int(raw_setting, com.raw_set.user_black);
}

static void _save_debayer(config_t *config, config_setting_t *root) {
	config_setting_t *debayer_group, *debayer_setting;

	debayer_group = config_setting_add(root, keywords[BAY], CONFIG_TYPE_GROUP);

	debayer_setting = config_setting_add(debayer_group, "ser_use_bayer_header",
			CONFIG_TYPE_BOOL);
	config_setting_set_bool(debayer_setting, com.debayer.use_bayer_header);

	debayer_setting = config_setting_add(debayer_group, "pattern",
			CONFIG_TYPE_INT);
	config_setting_set_int(debayer_setting, com.debayer.bayer_pattern);

	debayer_setting = config_setting_add(debayer_group, "compatibility",
			CONFIG_TYPE_BOOL);
	config_setting_set_bool(debayer_setting, com.debayer.compatibility);

	debayer_setting = config_setting_add(debayer_group, "inter",
			CONFIG_TYPE_INT);
	config_setting_set_int(debayer_setting, com.debayer.bayer_inter);

	debayer_setting = config_setting_add(debayer_group, "stretch",
			CONFIG_TYPE_BOOL);
	config_setting_set_bool(debayer_setting, com.debayer.stretch);
}

static void _save_preprocessing(config_t *config, config_setting_t *root) {
	config_setting_t *prepro_group, *prepro_setting;

	prepro_group = config_setting_add(root, keywords[PRE], CONFIG_TYPE_GROUP);

	prepro_setting = config_setting_add(prepro_group, "cfa", CONFIG_TYPE_BOOL);
	config_setting_set_bool(prepro_setting, com.prepro_cfa);

	prepro_setting = config_setting_add(prepro_group, "equalize_cfa", CONFIG_TYPE_BOOL);
	config_setting_set_bool(prepro_setting, com.prepro_equalize_cfa);
}

static void _save_registration(config_t *config, config_setting_t *root) {
	config_setting_t *reg_group, *reg_setting;

	reg_group = config_setting_add(root, keywords[REG], CONFIG_TYPE_GROUP);

	reg_setting = config_setting_add(reg_group, "method", CONFIG_TYPE_INT);
	config_setting_set_int(reg_setting, com.reg_settings);
}

static void _save_stacking(config_t *config, config_setting_t *root) {
	config_setting_t *stk_group, *stk_setting;

	stk_group = config_setting_add(root, keywords[STK], CONFIG_TYPE_GROUP);

	stk_setting = config_setting_add(stk_group, "method", CONFIG_TYPE_INT);
	config_setting_set_int(stk_setting, com.stack.method);

	stk_setting = config_setting_add(stk_group, "rejection", CONFIG_TYPE_INT);
	config_setting_set_int(stk_setting, com.stack.rej_method);

	stk_setting = config_setting_add(stk_group, "mem_mode", CONFIG_TYPE_INT);
	config_setting_set_int(stk_setting, com.stack.mem_mode);

	stk_setting = config_setting_add(stk_group, "maxmem", CONFIG_TYPE_FLOAT);
	config_setting_set_float(stk_setting, com.stack.memory_ratio);

	stk_setting = config_setting_add(stk_group, "maxmem_gb", CONFIG_TYPE_FLOAT);
	config_setting_set_float(stk_setting, com.stack.memory_amount);
}

static void _save_photometry(config_t *config, config_setting_t *root) {
	config_setting_t *photometry_group, *photometry_setting;

	photometry_group = config_setting_add(root, keywords[PTM], CONFIG_TYPE_GROUP);

	photometry_setting = config_setting_add(photometry_group, "gain", CONFIG_TYPE_FLOAT);
	config_setting_set_float(photometry_setting, com.phot_set.gain);
	photometry_setting = config_setting_add(photometry_group, "inner-radius", CONFIG_TYPE_FLOAT);
	config_setting_set_float(photometry_setting, com.phot_set.inner);
	photometry_setting = config_setting_add(photometry_group, "outer-radius", CONFIG_TYPE_FLOAT);
	config_setting_set_float(photometry_setting, com.phot_set.outer);
	photometry_setting = config_setting_add(photometry_group, "minval", CONFIG_TYPE_INT);
	config_setting_set_int(photometry_setting, com.phot_set.minval);
	photometry_setting = config_setting_add(photometry_group, "maxval", CONFIG_TYPE_INT);
	config_setting_set_int(photometry_setting, com.phot_set.maxval);
}

static void _save_misc(config_t *config, config_setting_t *root) {
	config_setting_t *misc_group, *misc_setting;
	GSList *list = com.script_path;

	misc_group = config_setting_add(root, keywords[MISC], CONFIG_TYPE_GROUP);

	misc_setting = config_setting_add(misc_group, "swap_directory",
			CONFIG_TYPE_STRING);
	config_setting_set_string(misc_setting, com.swap_dir);

	misc_setting = config_setting_add(misc_group, "extension",
			CONFIG_TYPE_STRING);
	config_setting_set_string(misc_setting, com.ext);

	misc_setting = config_setting_add(misc_group, "confirm", CONFIG_TYPE_BOOL);
	config_setting_set_bool(misc_setting, com.dontShowConfirm);

	misc_setting = config_setting_add(misc_group, "theme", CONFIG_TYPE_INT);
	config_setting_set_int(misc_setting, com.combo_theme);

	misc_setting = config_setting_add(misc_group, "remember_winpos",
			CONFIG_TYPE_BOOL);
	config_setting_set_bool(misc_setting, com.remember_windows);

	misc_setting = config_setting_add(misc_group, "scripts_paths",
			CONFIG_TYPE_LIST);
	while (list) {
		config_setting_set_string_elem(misc_setting, -1, (char *)list->data);
		list = list->next;
	}
	misc_setting = config_setting_add(misc_group, "main_w_pos",
			CONFIG_TYPE_LIST);
	config_setting_set_int_elem(misc_setting, -1, com.main_w_pos.x);
	config_setting_set_int_elem(misc_setting, -1, com.main_w_pos.y);
	config_setting_set_int_elem(misc_setting, -1, com.main_w_pos.w);
	config_setting_set_int_elem(misc_setting, -1, com.main_w_pos.h);
	misc_setting = config_setting_add(misc_group, "rgb_w_pos",
			CONFIG_TYPE_LIST);
	config_setting_set_int_elem(misc_setting, -1, com.rgb_w_pos.x);
	config_setting_set_int_elem(misc_setting, -1, com.rgb_w_pos.y);
	config_setting_set_int_elem(misc_setting, -1, com.rgb_w_pos.w);
	config_setting_set_int_elem(misc_setting, -1, com.rgb_w_pos.h);
}

static int siril_config_write_file(config_t *config, const char *filename) {
	gchar *fname = get_locale_filename(filename);
	int ret = config_write_file(config, fname);
	g_free(fname);
	return ret;
}

/**
 * Public functions
 */

int writeinitfile() {
	config_t config;
	config_setting_t *root;

	config_init(&config);
	root = config_root_setting(&config);

	_save_wd(&config, root);
	_save_libraw(&config, root);
	_save_debayer(&config, root);
	_save_preprocessing(&config, root);
	_save_registration(&config, root);
	_save_stacking(&config, root);
	_save_photometry(&config, root);
	_save_misc(&config, root);

	if (!siril_config_write_file(&config, com.initfile)) {
		fprintf(stderr, "Error while writing file.\n");
		config_destroy(&config);
		return 1;
	}
	config_destroy(&config);
	return 0;
}

int checkinitfile() {
	gchar *home = NULL;
	GStatBuf sts;

	// try to read the file given on command line
	if (com.initfile && !readinitfile()) {
		return 0;
	}

	// no file given on command line, set initfile to default location
#ifdef _WIN32
	home = g_build_filename (get_special_folder (CSIDL_APPDATA),
			"siril", NULL);
	if( g_mkdir_with_parents( home, 1 ) == 0 ) {
		fprintf( stderr, "Created homefolder %s!\n", home );
	} else {
		fprintf( stderr, "Failed to create homefolder %s!\n", com.initfile );
	}
#else
	const gchar *tmp = g_get_home_dir();
	if (tmp == NULL) {
		fprintf(stderr,
				"Could not get the environment variable $HOME, no config file.\n");
		return 1;
	}
	home = g_strdup(tmp);
#endif

#if (defined(__APPLE__) && defined(__MACH__))
	fprintf(stderr, "Creating initfile in Application Support.\n");
	gchar *homefolder;
	homefolder = g_build_filename(g_get_home_dir(),
			"Library", "Application Support", "siril", NULL);
	if (g_mkdir_with_parents(homefolder, 0755) == 0) {
		com.initfile = g_build_filename(homefolder, CFG_FILE, NULL);
		fprintf(stderr, "The initfile name is %s.\n", com.initfile);
	} else {
		fprintf(stderr, "Failed to create homefolder %s.\n", homefolder);
	}

#elif defined (_WIN32)
	com.initfile = g_build_filename(home, CFG_FILE, NULL);
#else
	com.initfile = g_new(gchar, strlen(home) + 20);
	sprintf(com.initfile, "%s/.siril/%s", home, CFG_FILE);
#endif
	if (readinitfile()) {	// couldn't read it
		char filename[255];

		// if that fails, check and create the default ini file
#if (defined(__APPLE__) && defined(__MACH__))
		snprintf(filename, 255, "%s", homefolder);
		g_free(homefolder);
#elif defined (_WIN32)
		snprintf(filename, 255, "%s", home);
#else
		snprintf(filename, 255, "%s/.siril", home);
#endif
		g_free(home);
		if (g_stat(filename, &sts) != 0) {
			if (errno == ENOENT) {
				if (g_mkdir(filename, 0755)) {
					fprintf(stderr, "Could not create dir %s, please check\n",
							filename);
					return 1;
				}
				com.swap_dir = g_strdup(g_get_tmp_dir());
				com.ext = strdup(".fit");
				return (writeinitfile());
			}
		}

		if (!(S_ISDIR(sts.st_mode))) {
			fprintf(stderr,
					"There is a file named %s, that is not a directory.\n"
							"Remove or rename it first\n", filename);
			return 1;
		}

		com.swap_dir = g_strdup(g_get_tmp_dir());
		com.ext = strdup(".fit");
		return (writeinitfile());
	}
	g_free(home);

	return 0;
}
