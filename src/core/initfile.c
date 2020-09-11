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

#include <libconfig.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_app_dirs.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"

#include "initfile.h"

#define CONFIG_FILE "siril.config"

static const char *keywords[] = { "working-directory", "libraw-settings",
		"debayer-settings", "prepro-settings", "registration-settings",
		"stacking-settings", "photometry-settings", "misc-settings", "compression-settings" };

static int readinitfile() {
	config_t config;
	const char *dir = NULL;
	GSList *list = NULL;

	if (!com.initfile)
		return 1;

	config_init(&config);

	if (config_read_file(&config, com.initfile) == CONFIG_FALSE) {
		return 1;
	}
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
		config_setting_lookup_float(raw_setting, "mul_0", &com.pref.raw_set.mul[0]);
		config_setting_lookup_float(raw_setting, "mul_2", &com.pref.raw_set.mul[2]);
		config_setting_lookup_float(raw_setting, "bright", &com.pref.raw_set.bright);
		config_setting_lookup_int(raw_setting, "auto", &com.pref.raw_set.auto_mul);
		config_setting_lookup_int(raw_setting, "cam_wb", &com.pref.raw_set.use_camera_wb);
		config_setting_lookup_int(raw_setting, "auto_wb", &com.pref.raw_set.use_auto_wb);
		config_setting_lookup_int(raw_setting, "user_qual", &com.pref.raw_set.user_qual);
		config_setting_lookup_float(raw_setting, "gamm_0", &com.pref.raw_set.gamm[0]);
		config_setting_lookup_float(raw_setting, "gamm_1", &com.pref.raw_set.gamm[1]);
		config_setting_lookup_int(raw_setting, "user_black", &com.pref.raw_set.user_black);
	}

	/* Debayer setting */
	config_setting_t *debayer_setting = config_lookup(&config, keywords[BAY]);
	if (debayer_setting) {
		config_setting_lookup_bool(debayer_setting, "ser_use_bayer_header", &com.pref.debayer.use_bayer_header);
		config_setting_lookup_int(debayer_setting, "pattern", &com.pref.debayer.bayer_pattern);
		config_setting_lookup_bool(debayer_setting, "debayer_top_down", &com.pref.debayer.top_down);
		config_setting_lookup_int(debayer_setting, "debayer_algo", (int*)&com.pref.debayer.bayer_inter);
		config_setting_lookup_int(debayer_setting, "xbayeroff", &com.pref.debayer.xbayeroff);
		config_setting_lookup_int(debayer_setting, "ybayeroff", &com.pref.debayer.ybayeroff);

	}

	/* Preprocessing settings */
	config_setting_t *prepro_setting = config_lookup(&config, keywords[PRE]);
	if (prepro_setting) {
		const char *bias = NULL, *dark = NULL, *flat = NULL;

		config_setting_lookup_bool(prepro_setting, "cfa", &com.pref.prepro_cfa);
		config_setting_lookup_bool(prepro_setting, "equalize_cfa", &com.pref.prepro_equalize_cfa);
		config_setting_lookup_bool(prepro_setting, "fix_xtrans", &com.pref.fix_xtrans);
		config_setting_lookup_string(prepro_setting, "bias_lib", &bias);
		com.pref.prepro_bias_lib = g_strdup(bias);
		config_setting_lookup_bool(prepro_setting, "use_bias_lib", &com.pref.use_bias_lib);
		config_setting_lookup_string(prepro_setting, "dark_lib", &dark);
		com.pref.prepro_dark_lib = g_strdup(dark);
		config_setting_lookup_bool(prepro_setting, "use_dark_lib", &com.pref.use_dark_lib);
		config_setting_lookup_string(prepro_setting, "flat_lib", &flat);
		com.pref.prepro_flat_lib = g_strdup(flat);
		config_setting_lookup_bool(prepro_setting, "use_flat_lib", &com.pref.use_flat_lib);
		prepro_setting = config_lookup(&config, "prepro-settings.xtrans_af");
		if (prepro_setting != NULL) {
			com.pref.xtrans_af.x = config_setting_get_int_elem(prepro_setting, 0);
			com.pref.xtrans_af.y = config_setting_get_int_elem(prepro_setting, 1);
			com.pref.xtrans_af.w = config_setting_get_int_elem(prepro_setting, 2);
			com.pref.xtrans_af.h = config_setting_get_int_elem(prepro_setting, 3);
		}
		prepro_setting = config_lookup(&config, "prepro-settings.xtrans_sample");
		if (prepro_setting != NULL) {
			com.pref.xtrans_sample.x = config_setting_get_int_elem(prepro_setting, 0);
			com.pref.xtrans_sample.y = config_setting_get_int_elem(prepro_setting, 1);
			com.pref.xtrans_sample.w = config_setting_get_int_elem(prepro_setting, 2);
			com.pref.xtrans_sample.h = config_setting_get_int_elem(prepro_setting, 3);
		}
	}

	/* Registration setting */
	config_setting_t *reg_setting = config_lookup(&config, keywords[REG]);
	if (reg_setting) {
		config_setting_lookup_int(reg_setting, "method", &com.reg_settings);
	}

	/* Stacking setting */
	config_setting_t *stack_setting = config_lookup(&config, keywords[STK]);
	if (stack_setting) {
		config_setting_lookup_int(stack_setting, "method", &com.pref.stack.method);
		config_setting_lookup_int(stack_setting, "rejection", &com.pref.stack.rej_method);
		config_setting_lookup_int(stack_setting, "normalisation", &com.pref.stack.normalisation_method);
		config_setting_lookup_float(stack_setting, "sigma_low", &com.pref.stack.sigma_low);
		config_setting_lookup_float(stack_setting, "sigma_high", &com.pref.stack.sigma_high);
		config_setting_lookup_float(stack_setting, "linear_low", &com.pref.stack.linear_low);
		config_setting_lookup_float(stack_setting, "linear_high", &com.pref.stack.linear_high);
		config_setting_lookup_float(stack_setting, "percentile_low", &com.pref.stack.percentile_low);
		config_setting_lookup_float(stack_setting, "percentile_high", &com.pref.stack.percentile_high);


		config_setting_lookup_int(stack_setting, "mem_mode", (int*)&com.pref.stack.mem_mode);
		config_setting_lookup_float(stack_setting, "maxmem", &com.pref.stack.memory_ratio);
		config_setting_lookup_float(stack_setting, "maxmem_gb",	&com.pref.stack.memory_amount);
	}
	if (com.pref.stack.mem_mode < 0 || com.pref.stack.mem_mode > 2)
		com.pref.stack.mem_mode = RATIO;
	if (com.pref.stack.memory_ratio <= 0.05)
		com.pref.stack.memory_ratio = 0.9;

	/* FITS compression setting */
	config_setting_t *comp_setting = config_lookup(&config, keywords[CMP]);
	if (comp_setting) {
		config_setting_lookup_bool(comp_setting, "fits_enabled", &com.pref.comp.fits_enabled);
		config_setting_lookup_int(comp_setting, "fits_method", &com.pref.comp.fits_method);
		config_setting_lookup_float(comp_setting, "fits_quantization", &com.pref.comp.fits_quantization);
		config_setting_lookup_float(comp_setting, "fits_hcompress_scale", &com.pref.comp.fits_hcompress_scale);
	}

	/* Photometry setting */
	config_setting_t *photometry_setting = config_lookup(&config, keywords[PTM]);
	if (photometry_setting) {
		config_setting_lookup_float(photometry_setting, "gain", &com.pref.phot_set.gain);
		config_setting_lookup_float(photometry_setting, "inner-radius", &com.pref.phot_set.inner);
		config_setting_lookup_float(photometry_setting, "outer-radius", &com.pref.phot_set.outer);
		config_setting_lookup_int(photometry_setting, "minval", &com.pref.phot_set.minval);
		config_setting_lookup_int(photometry_setting, "maxval", &com.pref.phot_set.maxval);
	}

	/* Misc setting */
	config_setting_t *misc_setting = config_lookup(&config, keywords[MISC]);
	if (misc_setting) {
		int type;
		const char *swap_dir = NULL, *extension = NULL, *lang = NULL, *copyright = NULL;

		config_setting_lookup_bool(misc_setting, "first_start_1_0_0", &com.pref.first_start);
		config_setting_lookup_bool(misc_setting, "confirm_quit", &com.pref.save.quit);
		config_setting_lookup_bool(misc_setting, "confirm_script", &com.pref.save.script);
		config_setting_lookup_bool(misc_setting, "show_thumbnails", &com.pref.show_thumbnails);
		config_setting_lookup_int(misc_setting, "thumbnail_size", &com.pref.thumbnail_size);
		config_setting_lookup_int(misc_setting, "theme", &com.pref.combo_theme);
		config_setting_lookup_string(misc_setting, "lang", &lang);
		com.pref.combo_lang = g_strdup(lang);
		config_setting_lookup_bool(misc_setting, "remember_winpos", &com.pref.remember_windows);
		config_setting_lookup_bool(misc_setting, "is_maximized", &com.pref.is_maximized);
		config_setting_lookup_string(misc_setting, "swap_directory", &swap_dir);
		com.pref.swap_dir = g_strdup(swap_dir);
		config_setting_lookup_string(misc_setting, "extension", &extension);
		com.pref.ext = g_strdup(extension);
		config_setting_lookup_int(misc_setting, "FITS_type", &type);
		com.pref.force_to_16bit = (type == 0);
		config_setting_lookup_int(misc_setting, "selection_guides", &com.pref.selection_guides);
		config_setting_lookup_string(misc_setting, "copyright", &copyright);
		com.pref.copyright = g_strdup(copyright);
		config_setting_lookup_bool(misc_setting, "check_update", &com.pref.check_update);

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
			com.pref.main_w_pos.x = config_setting_get_int_elem(misc_setting, 0);
			com.pref.main_w_pos.y = config_setting_get_int_elem(misc_setting, 1);
			com.pref.main_w_pos.w = config_setting_get_int_elem(misc_setting, 2);
			com.pref.main_w_pos.h = config_setting_get_int_elem(misc_setting, 3);
		}
	}
	com.pref.script_path = list;
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
	config_setting_set_float(raw_setting, com.pref.raw_set.mul[0]);

	raw_setting = config_setting_add(libraw_group, "mul_2", CONFIG_TYPE_FLOAT);
	config_setting_set_float(raw_setting, com.pref.raw_set.mul[2]);

	raw_setting = config_setting_add(libraw_group, "bright", CONFIG_TYPE_FLOAT);
	config_setting_set_float(raw_setting, com.pref.raw_set.bright);

	raw_setting = config_setting_add(libraw_group, "auto", CONFIG_TYPE_INT);
	config_setting_set_int(raw_setting, com.pref.raw_set.auto_mul);

	raw_setting = config_setting_add(libraw_group, "cam_wb", CONFIG_TYPE_INT);
	config_setting_set_int(raw_setting, com.pref.raw_set.use_camera_wb);

	raw_setting = config_setting_add(libraw_group, "auto_wb", CONFIG_TYPE_INT);
	config_setting_set_int(raw_setting, com.pref.raw_set.use_auto_wb);

	raw_setting = config_setting_add(libraw_group, "user_qual",	CONFIG_TYPE_INT);
	config_setting_set_int(raw_setting, com.pref.raw_set.user_qual);

	raw_setting = config_setting_add(libraw_group, "gamm_0", CONFIG_TYPE_FLOAT);
	config_setting_set_float(raw_setting, com.pref.raw_set.gamm[0]);

	raw_setting = config_setting_add(libraw_group, "gamm_1", CONFIG_TYPE_FLOAT);
	config_setting_set_float(raw_setting, com.pref.raw_set.gamm[1]);

	raw_setting = config_setting_add(libraw_group, "user_black", CONFIG_TYPE_INT);
	config_setting_set_int(raw_setting, com.pref.raw_set.user_black);
}

static void _save_debayer(config_t *config, config_setting_t *root) {
	config_setting_t *debayer_group, *debayer_setting;

	debayer_group = config_setting_add(root, keywords[BAY], CONFIG_TYPE_GROUP);

	debayer_setting = config_setting_add(debayer_group, "ser_use_bayer_header",	CONFIG_TYPE_BOOL);
	config_setting_set_bool(debayer_setting, com.pref.debayer.use_bayer_header);

	debayer_setting = config_setting_add(debayer_group, "pattern", CONFIG_TYPE_INT);
	config_setting_set_int(debayer_setting, com.pref.debayer.bayer_pattern);

	debayer_setting = config_setting_add(debayer_group, "debayer_top_down", CONFIG_TYPE_BOOL);
	config_setting_set_bool(debayer_setting, com.pref.debayer.top_down);

	debayer_setting = config_setting_add(debayer_group, "debayer_algo", CONFIG_TYPE_INT);
	config_setting_set_int(debayer_setting, com.pref.debayer.bayer_inter);

	debayer_setting = config_setting_add(debayer_group, "xbayeroff", CONFIG_TYPE_INT);
	config_setting_set_int(debayer_setting, com.pref.debayer.xbayeroff);
	debayer_setting = config_setting_add(debayer_group, "ybayeroff", CONFIG_TYPE_INT);
	config_setting_set_int(debayer_setting, com.pref.debayer.ybayeroff);
}

static void _save_preprocessing(config_t *config, config_setting_t *root) {
	config_setting_t *prepro_group, *prepro_setting;

	prepro_group = config_setting_add(root, keywords[PRE], CONFIG_TYPE_GROUP);

	prepro_setting = config_setting_add(prepro_group, "cfa", CONFIG_TYPE_BOOL);
	config_setting_set_bool(prepro_setting, com.pref.prepro_cfa);

	prepro_setting = config_setting_add(prepro_group, "equalize_cfa", CONFIG_TYPE_BOOL);
	config_setting_set_bool(prepro_setting, com.pref.prepro_equalize_cfa);

	prepro_setting = config_setting_add(prepro_group, "fix_xtrans", CONFIG_TYPE_BOOL);
	config_setting_set_bool(prepro_setting, com.pref.fix_xtrans);

	prepro_setting = config_setting_add(prepro_group, "xtrans_af",	CONFIG_TYPE_LIST);
	config_setting_set_int_elem(prepro_setting, -1, com.pref.xtrans_af.x);
	config_setting_set_int_elem(prepro_setting, -1, com.pref.xtrans_af.y);
	config_setting_set_int_elem(prepro_setting, -1, com.pref.xtrans_af.w);
	config_setting_set_int_elem(prepro_setting, -1, com.pref.xtrans_af.h);

	prepro_setting = config_setting_add(prepro_group, "xtrans_sample",	CONFIG_TYPE_LIST);
	config_setting_set_int_elem(prepro_setting, -1, com.pref.xtrans_sample.x);
	config_setting_set_int_elem(prepro_setting, -1, com.pref.xtrans_sample.y);
	config_setting_set_int_elem(prepro_setting, -1, com.pref.xtrans_sample.w);
	config_setting_set_int_elem(prepro_setting, -1, com.pref.xtrans_sample.h);

	prepro_setting = config_setting_add(prepro_group, "bias_lib", CONFIG_TYPE_STRING);
	config_setting_set_string(prepro_setting, com.pref.prepro_bias_lib);
	prepro_setting = config_setting_add(prepro_group, "use_bias_lib", CONFIG_TYPE_BOOL);
	config_setting_set_bool(prepro_setting, com.pref.use_bias_lib);

	prepro_setting = config_setting_add(prepro_group, "dark_lib", CONFIG_TYPE_STRING);
	config_setting_set_string(prepro_setting, com.pref.prepro_dark_lib);
	prepro_setting = config_setting_add(prepro_group, "use_dark_lib", CONFIG_TYPE_BOOL);
	config_setting_set_bool(prepro_setting, com.pref.use_dark_lib);

	prepro_setting = config_setting_add(prepro_group, "flat_lib", CONFIG_TYPE_STRING);
	config_setting_set_string(prepro_setting, com.pref.prepro_flat_lib);
	prepro_setting = config_setting_add(prepro_group, "use_flat_lib", CONFIG_TYPE_BOOL);
	config_setting_set_bool(prepro_setting, com.pref.use_flat_lib);
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
	config_setting_set_int(stk_setting, com.pref.stack.method);

	stk_setting = config_setting_add(stk_group, "rejection", CONFIG_TYPE_INT);
	config_setting_set_int(stk_setting, com.pref.stack.rej_method);

	stk_setting = config_setting_add(stk_group, "sigma_low", CONFIG_TYPE_FLOAT);
	config_setting_set_float(stk_setting, com.pref.stack.sigma_low);

	stk_setting = config_setting_add(stk_group, "sigma_high", CONFIG_TYPE_FLOAT);
	config_setting_set_float(stk_setting, com.pref.stack.sigma_high);

	stk_setting = config_setting_add(stk_group, "linear_low", CONFIG_TYPE_FLOAT);
	config_setting_set_float(stk_setting, com.pref.stack.linear_low);

	stk_setting = config_setting_add(stk_group, "linear_high", CONFIG_TYPE_FLOAT);
	config_setting_set_float(stk_setting, com.pref.stack.linear_high);

	stk_setting = config_setting_add(stk_group, "percentile_low", CONFIG_TYPE_FLOAT);
	config_setting_set_float(stk_setting, com.pref.stack.percentile_low);

	stk_setting = config_setting_add(stk_group, "percentile_high", CONFIG_TYPE_FLOAT);
	config_setting_set_float(stk_setting, com.pref.stack.percentile_high);

	stk_setting = config_setting_add(stk_group, "mem_mode", CONFIG_TYPE_INT);
	config_setting_set_int(stk_setting, com.pref.stack.mem_mode);

	stk_setting = config_setting_add(stk_group, "maxmem", CONFIG_TYPE_FLOAT);
	config_setting_set_float(stk_setting, com.pref.stack.memory_ratio);

	stk_setting = config_setting_add(stk_group, "maxmem_gb", CONFIG_TYPE_FLOAT);
	config_setting_set_float(stk_setting, com.pref.stack.memory_amount);
}

static void _save_comp(config_t *config, config_setting_t *root) {
	config_setting_t *cmp_group, *cmp_setting;

	cmp_group = config_setting_add(root, keywords[CMP], CONFIG_TYPE_GROUP);

	cmp_setting = config_setting_add(cmp_group, "fits_enabled", CONFIG_TYPE_BOOL);
	config_setting_set_bool(cmp_setting, com.pref.comp.fits_enabled);

	cmp_setting = config_setting_add(cmp_group, "fits_method", CONFIG_TYPE_INT);
	config_setting_set_int(cmp_setting, com.pref.comp.fits_method);

	cmp_setting = config_setting_add(cmp_group, "fits_quantization", CONFIG_TYPE_FLOAT);
	config_setting_set_float(cmp_setting, com.pref.comp.fits_quantization);

	cmp_setting = config_setting_add(cmp_group, "fits_hcompress_scale", CONFIG_TYPE_FLOAT);
	config_setting_set_float(cmp_setting, com.pref.comp.fits_hcompress_scale);

}

static void _save_photometry(config_t *config, config_setting_t *root) {
	config_setting_t *photometry_group, *photometry_setting;

	photometry_group = config_setting_add(root, keywords[PTM], CONFIG_TYPE_GROUP);

	photometry_setting = config_setting_add(photometry_group, "gain", CONFIG_TYPE_FLOAT);
	config_setting_set_float(photometry_setting, com.pref.phot_set.gain);
	photometry_setting = config_setting_add(photometry_group, "inner-radius", CONFIG_TYPE_FLOAT);
	config_setting_set_float(photometry_setting, com.pref.phot_set.inner);
	photometry_setting = config_setting_add(photometry_group, "outer-radius", CONFIG_TYPE_FLOAT);
	config_setting_set_float(photometry_setting, com.pref.phot_set.outer);
	photometry_setting = config_setting_add(photometry_group, "minval", CONFIG_TYPE_INT);
	config_setting_set_int(photometry_setting, com.pref.phot_set.minval);
	photometry_setting = config_setting_add(photometry_group, "maxval", CONFIG_TYPE_INT);
	config_setting_set_int(photometry_setting, com.pref.phot_set.maxval);
}

static void _save_misc(config_t *config, config_setting_t *root) {
	config_setting_t *misc_group, *misc_setting;
	GSList *list = com.pref.script_path;

	misc_group = config_setting_add(root, keywords[MISC], CONFIG_TYPE_GROUP);

	misc_setting = config_setting_add(misc_group, "swap_directory", CONFIG_TYPE_STRING);
	config_setting_set_string(misc_setting, com.pref.swap_dir);

	misc_setting = config_setting_add(misc_group, "first_start_1_0_0", CONFIG_TYPE_BOOL);
	config_setting_set_bool(misc_setting, com.pref.first_start);

	misc_setting = config_setting_add(misc_group, "extension", CONFIG_TYPE_STRING);
	config_setting_set_string(misc_setting, com.pref.ext);

	misc_setting = config_setting_add(misc_group, "FITS_type", CONFIG_TYPE_INT);
	config_setting_set_int(misc_setting, com.pref.force_to_16bit ? 0 : 1);

	misc_setting = config_setting_add(misc_group, "selection_guides", CONFIG_TYPE_INT);
	config_setting_set_int(misc_setting, com.pref.selection_guides);

	misc_setting = config_setting_add(misc_group, "copyright", CONFIG_TYPE_STRING);
	config_setting_set_string(misc_setting, com.pref.copyright);

	misc_setting = config_setting_add(misc_group, "confirm_quit", CONFIG_TYPE_BOOL);
	config_setting_set_bool(misc_setting, com.pref.save.quit);

	misc_setting = config_setting_add(misc_group, "confirm_script", CONFIG_TYPE_BOOL);
	config_setting_set_bool(misc_setting, com.pref.save.script);

	misc_setting = config_setting_add(misc_group, "show_thumbnails", CONFIG_TYPE_BOOL);
	config_setting_set_bool(misc_setting, com.pref.show_thumbnails);

	misc_setting = config_setting_add(misc_group, "thumbnail_size", CONFIG_TYPE_INT);
	config_setting_set_int(misc_setting, com.pref.thumbnail_size);

	misc_setting = config_setting_add(misc_group, "theme", CONFIG_TYPE_INT);
	config_setting_set_int(misc_setting, com.pref.combo_theme);

	misc_setting = config_setting_add(misc_group, "lang", CONFIG_TYPE_STRING);
	config_setting_set_string(misc_setting, com.pref.combo_lang);

	misc_setting = config_setting_add(misc_group, "remember_winpos", CONFIG_TYPE_BOOL);
	config_setting_set_bool(misc_setting, com.pref.remember_windows);

	misc_setting = config_setting_add(misc_group, "scripts_paths", CONFIG_TYPE_LIST);
	while (list) {
		config_setting_set_string_elem(misc_setting, -1, (char *)list->data);
		list = list->next;
	}
	misc_setting = config_setting_add(misc_group, "main_w_pos",	CONFIG_TYPE_LIST);
	config_setting_set_int_elem(misc_setting, -1, com.pref.main_w_pos.x);
	config_setting_set_int_elem(misc_setting, -1, com.pref.main_w_pos.y);
	config_setting_set_int_elem(misc_setting, -1, com.pref.main_w_pos.w);
	config_setting_set_int_elem(misc_setting, -1, com.pref.main_w_pos.h);

	misc_setting = config_setting_add(misc_group, "is_maximized", CONFIG_TYPE_BOOL);
	config_setting_set_bool(misc_setting, com.pref.is_maximized);

	misc_setting = config_setting_add(misc_group, "check_update", CONFIG_TYPE_BOOL);
	config_setting_set_bool(misc_setting, com.pref.check_update);
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
	_save_comp(&config, root);
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
	/* First we try to read the file given on command line */
	if (com.initfile && !readinitfile()) {
		return 0;
	}
	/* no file given on command line, set initfile to default location */
	gchar *pathname = g_build_filename(siril_get_config_dir(), PACKAGE, NULL);
	gchar *config_file = g_build_filename(pathname, CONFIG_FILE, NULL);
	if (!g_file_test(config_file, G_FILE_TEST_EXISTS)) {
		if (g_mkdir_with_parents(pathname, 0755) == 0) {
			g_fprintf(stderr, "Created config dir %s\n", pathname);
		} else {
			g_fprintf(stderr, "Failed to create config dir %s!\n", pathname);
			g_free(pathname);
			return 1;
		}
	}
	g_free(pathname);

	com.initfile = g_strdup(config_file);

	if (readinitfile()) {
		/* init file does not exist, so we create it */
		return writeinitfile();
	}
	return 0;
}
