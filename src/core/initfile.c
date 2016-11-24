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

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "gui/callbacks.h"

#define CFG_FILE "siril.cfg"

static const char *keywords[] = { "working-directory", "libraw-settings",
		"debayer-settings", "prepro-settings", "registration-settings",
		"stacking-settings", "misc-settings" };

static int readinitfile() {
	config_t config;
	const char *dir = NULL, *swap_dir = NULL, *extension = NULL;

	if (!com.initfile)
		return 1;

	config_init(&config);

	if (config_read_file(&config, com.initfile) == CONFIG_FALSE)
		return 1;
	siril_log_message(_("Loading init file: '%s'\n"), com.initfile);

	/* Working directory */
	if (config_lookup_string(&config, keywords[WD], &dir)) {
		if (changedir(dir)) {
			siril_log_message(
					_("Reverting current working directory to startup directory, the saved directory is not available anymore\n"));
			set_GUI_CWD();
			writeinitfile();
		}
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
	}

	/* Preprocessing settings */
	config_setting_t *prepro_setting = config_lookup(&config, keywords[PRE]);
	if (prepro_setting) {
		config_setting_lookup_bool(prepro_setting, "cfa", &com.prepro_cfa);
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
		config_setting_lookup_float(stack_setting, "maxmem",
				&com.stack.memory_percent);
	}
	if (com.stack.memory_percent <= 0.0001)
		com.stack.memory_percent = 0.9;

	/* Misc setting */
	config_setting_t *misc_setting = config_lookup(&config, keywords[MISC]);
	if (misc_setting) {
		config_setting_lookup_bool(misc_setting, "confirm",
				&com.dontShowConfirm);
		config_setting_lookup_bool(misc_setting, "darktheme",
				&com.have_dark_theme);
		config_setting_lookup_string(misc_setting, "swap_directory", &swap_dir);
		config_setting_lookup_string(misc_setting, "extension", &extension);
		set_GUI_misc();
	}
	if (swap_dir && swap_dir[0] != '\0') {
		if (com.swap_dir)
			free(com.swap_dir);
		com.swap_dir = strdup(swap_dir);
	} else {
		const char* sw_dir = g_get_tmp_dir();
		com.swap_dir = strdup(sw_dir);
	}
	if (extension && extension[0] != '\0') {
		if (com.ext)
			free(com.ext);
		com.ext = strdup(extension);
	} else {
		com.ext = strdup(".fit");
	}
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
}

static void _save_preprocessing(config_t *config, config_setting_t *root) {
	config_setting_t *prepro_group, *prepro_setting;

	prepro_group = config_setting_add(root, keywords[PRE], CONFIG_TYPE_GROUP);

	prepro_setting = config_setting_add(prepro_group, "cfa", CONFIG_TYPE_BOOL);
	config_setting_set_bool(prepro_setting, com.prepro_cfa);
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

	stk_setting = config_setting_add(stk_group, "maxmem", CONFIG_TYPE_FLOAT);
	config_setting_set_float(stk_setting, com.stack.memory_percent);
}

static void _save_misc(config_t *config, config_setting_t *root) {
	config_setting_t *misc_group, *misc_setting;

	misc_group = config_setting_add(root, keywords[MISC], CONFIG_TYPE_GROUP);

	misc_setting = config_setting_add(misc_group, "swap_directory",
			CONFIG_TYPE_STRING);
	config_setting_set_string(misc_setting, com.swap_dir);

	misc_setting = config_setting_add(misc_group, "extension",
			CONFIG_TYPE_STRING);
	config_setting_set_string(misc_setting, com.ext);

	misc_setting = config_setting_add(misc_group, "confirm", CONFIG_TYPE_BOOL);
	config_setting_set_bool(misc_setting, com.dontShowConfirm);

	misc_setting = config_setting_add(misc_group, "darktheme",
			CONFIG_TYPE_BOOL);
	config_setting_set_bool(misc_setting, com.have_dark_theme);
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
	_save_misc(&config, root);

	if (!config_write_file(&config, com.initfile)) {
		fprintf(stderr, "Error while writing file.\n");
		config_destroy(&config);
		return 1;
	}
	config_destroy(&config);

	return 0;
}

int checkinitfile() {
	char *home;
	struct stat sts;

	// try to read the file given on command line
	if (com.initfile && !readinitfile()) {
		return 0;
	}

	// no file given on command line, set initfile to default location
	if ((home = getenv("HOME")) == NULL) {
		fprintf(stderr,
				"Could not get the environment variable $HOME, no config file.\n");
		return 1;
	}
	com.initfile = malloc(strlen(home) + 20);
	sprintf(com.initfile, "%s/.siril/%s", home, CFG_FILE);
	if (readinitfile()) {	// couldn't read it
		char filename[255];

		set_GUI_CWD();
		// if that fails, check and create the default ini file
		snprintf(filename, 255, "%s/.siril", home);
		if (stat(filename, &sts) != 0) {
			if (errno == ENOENT) {
				if (mkdir(filename, 0755)) {
					fprintf(stderr, "Could not create dir %s, please check\n",
							filename);
					return 1;
				}
				com.swap_dir = strdup(g_get_tmp_dir());
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

		com.swap_dir = strdup(g_get_tmp_dir());
		com.ext = strdup(".fit");
		return (writeinitfile());
	}
	return 0;
}
