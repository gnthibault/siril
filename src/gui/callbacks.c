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

#include <gtk/gtk.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdarg.h>
#include <math.h>	// for M_PI
#include <libgen.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/undo.h"
#include "core/preprocess.h"
#include "core/command.h"
#include "core/command_line_processor.h"
#include "core/siril_app_dirs.h"
#include "algos/demosaicing.h"
#include "algos/colors.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/Def_Wavelet.h"
#include "algos/sorting.h"
#include "algos/background_extraction.h"
#include "io/conversion.h"
#include "io/films.h"
#include "io/sequence.h"
#include "io/ser.h"
#include "io/single_image.h"
#include "registration/registration.h"
#include "stacking/stacking.h"
#include "compositing/compositing.h"
#include "compositing/align_rgb.h"
#include "opencv/opencv.h"
#include "image_display.h"
#include "image_interactions.h"

#include "callbacks.h"
#include "message_dialog.h"
#include "PSF_list.h"
#include "sequence_list.h"
#include "histogram.h"
#include "script_menu.h"
#include "progress_and_log.h"
#include "dialogs.h"
#include "plot.h"

layer_info predefined_layers_colors[] = {
	/* name, lambda, lo, hi, c/over, c/under, mode */
	{ N_("Luminance"), 0.0, 0, 0, FALSE, FALSE, NORMAL_DISPLAY }, // no color, undefined value is <0
	{ N_("Red"), 650.0, 0, 0, FALSE, FALSE, NORMAL_DISPLAY }, // approx. of the middle of the color
	{ N_("Green"), 530.0, 0, 0, FALSE, FALSE, NORMAL_DISPLAY },	// approx. of the middle of the color
	{ N_("Blue"), 450.0, 0, 0, FALSE, FALSE, NORMAL_DISPLAY }// approx. of the middle of the color
};

/*****************************************************************************
 *                    S T A T I C      F U N C T I O N S                     *
 ****************************************************************************/
/*
 * Memory label static functions
 */

struct _label_data {
	const char *label_name;
	char *text;
};

static gboolean set_label_text_idle(gpointer p) {
	struct _label_data *args = (struct _label_data *) p;
	GtkLabel *label = GTK_LABEL(
			gtk_builder_get_object(builder, args->label_name));
	gtk_label_set_text(label, args->text);
	free(args->text);
	free(args);
	return FALSE;
}

static void set_label_text_from_main_thread(const char *label_name, const char *text) {
	struct _label_data *data = malloc(sizeof(struct _label_data));
	data->label_name = label_name;
	data->text = strdup(text);
	gdk_threads_add_idle(set_label_text_idle, data);
}

/* enables or disables the "display reference" checkbox in registration preview */
void enable_view_reference_checkbox(gboolean status) {
	static GtkToggleButton *check_display_ref = NULL;
	static GtkWidget *widget = NULL, *labelRegRef = NULL;
	if (check_display_ref == NULL) {
		check_display_ref = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_displayref"));
		widget = GTK_WIDGET(check_display_ref);
		labelRegRef = lookup_widget("labelRegRef");
	}
	if (status && gtk_widget_get_sensitive(widget))
		return;	// may be already enabled but deactivated by user, don't force it again
	gtk_widget_set_sensitive(widget, status);
	gtk_widget_set_visible(labelRegRef, !status);
	gtk_toggle_button_set_active(check_display_ref, status);
}

void set_viewer_mode_widgets_sensitive(gboolean sensitive) {
	static GtkWidget *scalemax = NULL;
	static GtkWidget *scalemin = NULL;
	static GtkWidget *entrymin = NULL;
	static GtkWidget *entrymax = NULL;
	static GtkWidget *minmax = NULL;
	static GtkWidget *hilo = NULL;
	static GtkWidget *user = NULL;

	if (!scalemax) {
		scalemax = lookup_widget("scalemax");
		scalemin = lookup_widget("scalemin");
		entrymin = lookup_widget("min_entry");
		entrymax = lookup_widget("max_entry");
		minmax = lookup_widget("radiobutton_minmax");
		hilo = lookup_widget("radiobutton_hilo");
		user = lookup_widget("radiobutton_user");
	}
	gtk_widget_set_sensitive(scalemax, sensitive);
	gtk_widget_set_sensitive(scalemin, sensitive);
	gtk_widget_set_sensitive(entrymin, sensitive);
	gtk_widget_set_sensitive(entrymax, sensitive);
	gtk_widget_set_sensitive(minmax, sensitive);
	gtk_widget_set_sensitive(hilo, sensitive);
	gtk_widget_set_sensitive(user, sensitive);
}

/* vport can be -1 if the correct viewport should be tested */
void test_and_allocate_reference_image(int vport) {
	static GtkComboBox *cbbt_layers = NULL;
	if (cbbt_layers == NULL) {
		cbbt_layers = GTK_COMBO_BOX(lookup_widget("comboboxreglayer"));
	}
	if (vport == -1)
		vport = gtk_combo_box_get_active(cbbt_layers);

	if (sequence_is_loaded() && com.seq.current == com.seq.reference_image
			&& gtk_combo_box_get_active(cbbt_layers) == vport) {
		/* this is the registration layer and the reference frame,
		 * save the buffer for alignment preview */
		if (!com.refimage_regbuffer || !com.refimage_surface) {
			guchar *oldbuf = com.refimage_regbuffer;
			com.refimage_regbuffer = realloc(com.refimage_regbuffer,
					com.surface_stride[vport] * gfit.ry * sizeof(guchar));
			if (com.refimage_regbuffer == NULL) {
				PRINT_ALLOC_ERR;
				if (oldbuf)
					free(oldbuf);
				return;
			}

			if (com.refimage_surface)
				cairo_surface_destroy(com.refimage_surface);
			com.refimage_surface = cairo_image_surface_create_for_data(
					com.refimage_regbuffer, CAIRO_FORMAT_RGB24, gfit.rx,
					gfit.ry, com.surface_stride[vport]);
			if (cairo_surface_status(com.refimage_surface)
					!= CAIRO_STATUS_SUCCESS) {
				fprintf(stderr,
						"Error creating the Cairo image surface for the reference image.\n");
				cairo_surface_destroy(com.refimage_surface);
				com.refimage_surface = NULL;
			} else {
				fprintf(stdout,
						"Saved the reference frame buffer for alignment preview.\n");
				enable_view_reference_checkbox(TRUE);
			}
		}
		memcpy(com.refimage_regbuffer, com.graybuf[vport],
				com.surface_stride[vport] * gfit.ry * sizeof(guchar));
		cairo_surface_flush(com.refimage_surface);
		cairo_surface_mark_dirty(com.refimage_surface);
	}
}

/*
 * Update FWHM UNITS static function
 */

static void update_fwhm_units_ok() {
	GtkWidget *label_ok = GTK_WIDGET(lookup_widget("label_ok"));
	gboolean update;

	update = gfit.focal_length > 0.0 && gfit.pixel_size_x > 0.0f && gfit.pixel_size_y > 0.0f;

	gtk_widget_set_visible(label_ok, update);
}

/*
 * Set swap dir to default
 */

static void reset_swapdir() {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_swap"));
	const gchar *dir;

	dir = g_get_tmp_dir();

	if (g_strcmp0(dir, com.swap_dir)) {
		g_free(com.swap_dir);
		com.swap_dir = g_strdup(dir);
		gtk_file_chooser_set_filename(swap_dir, dir);
		writeinitfile();
	}
}

/*****************************************************************************
 *                    P U B L I C      F U N C T I O N S                     *
 ****************************************************************************/

GtkWidget* lookup_widget(const gchar *widget_name) {
	return GTK_WIDGET(gtk_builder_get_object(builder, widget_name));
}

static void update_theme_button(const gchar *button_name, const gchar *path) {
	gchar *image;
	GtkWidget *w_image;

	image = g_build_filename(siril_get_system_data_dir(), "pixmaps", path, NULL);
	w_image = gtk_image_new_from_file(image);
	gtk_widget_show(w_image);
	gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(lookup_widget(button_name)), w_image);
	gtk_widget_show(lookup_widget(button_name));

	g_free(image);
}

static void update_icons_to_theme(gboolean is_dark) {
	siril_debug_print("Loading %s theme...\n", is_dark ? "dark" : "light");
	if (is_dark) {
		update_theme_button("rotate90_anticlock_button", "rotate-acw_dark.png");
		update_theme_button("rotate90_clock_button", "rotate-cw_dark.png");
		update_theme_button("mirrorx_button", "mirrorx_dark.png");
		update_theme_button("mirrory_button", "mirrory_dark.png");

		update_theme_button("process_starfinder_button", "starfinder_dark.png");
		update_theme_button("sum_button", "sum_dark.png");
		update_theme_button("export_button", "export_dark.png");

		update_theme_button("histoToolAutoStretch", "mtf_dark.png");
	} else {
		update_theme_button("rotate90_anticlock_button", "rotate-acw.png");
		update_theme_button("rotate90_clock_button", "rotate-cw.png");
		update_theme_button("mirrorx_button", "mirrorx.png");
		update_theme_button("mirrory_button", "mirrory.png");

		update_theme_button("process_starfinder_button", "starfinder.png");
		update_theme_button("sum_button", "sum.png");
		update_theme_button("export_button", "export.png");

		update_theme_button("histoToolAutoStretch", "mtf.png");
	}
}

void on_combo_theme_changed(GtkComboBox *box, gpointer user_data) {
	GtkSettings *settings;

	com.combo_theme = gtk_combo_box_get_active(box);

	settings = gtk_settings_get_default();
	g_object_set(settings, "gtk-application-prefer-dark-theme", com.combo_theme == 0, NULL);
	update_icons_to_theme(com.combo_theme == 0);
}

static void initialize_theme_GUI() {
	GtkComboBox *box;

	box = GTK_COMBO_BOX(lookup_widget("combo_theme"));

	g_signal_handlers_block_by_func(box, on_combo_theme_changed, NULL);
	gtk_combo_box_set_active(box, com.combo_theme);
	g_signal_handlers_unblock_by_func(box, on_combo_theme_changed, NULL);
	update_icons_to_theme(com.combo_theme == 0);
}

void load_prefered_theme(gint theme) {
	GtkSettings *settings;

	settings = gtk_settings_get_default();

	g_object_set(settings, "gtk-application-prefer-dark-theme", com.combo_theme == 0, NULL);
}

void set_sliders_value_to_gfit() {
	static GtkAdjustment *adj1 = NULL, *adj2 = NULL;

	if (adj1 == NULL) {
		adj1 = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment1"));// scalemax
		adj2 = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment2"));// scalemin
	}

	gfit.hi = gtk_adjustment_get_value(adj1);
	gfit.lo = gtk_adjustment_get_value(adj2);
}

/* Sets maximum value for contrast scales. Minimum is always 0.
 * Should be done on first image load of a sequence and when single images are loaded.
 * Max value is taken from gfit.maxi, recomputed if not present.
 */
void set_cutoff_sliders_max_values() {
	static GtkAdjustment *adj1 = NULL, *adj2 = NULL;
	gdouble max_val;
	if (adj1 == NULL) {
		adj1 = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment1"));// scalemax
		adj2 = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment2"));// scalemin
	}
	fprintf(stdout, _("Setting MAX value for cutoff sliders adjustments\n"));
	/* set max value for range according to number of bits of original image
	 * We should use gfit.bitpix for this, but it's currently always USHORT_IMG.
	 * Since 0.9.8 we have orig_bitpix, but it's not filled for SER and other images.
	 */
	max_val = (double)get_normalized_value(&gfit);
	gtk_adjustment_set_upper(adj1, max_val);
	gtk_adjustment_set_upper(adj2, max_val);
}

/* Sets the value scalemin and scalemax sliders so that they match the hi and
 * lo values stored for the current viewport.
 * It also sets the 'cut' checkbox status.
 * The function should be called on every tab change and file opening */
void set_cutoff_sliders_values() {
	gchar buffer[10];
	static GtkAdjustment *adjmin = NULL, *adjmax = NULL;
	static GtkEntry *maxentry = NULL, *minentry = NULL;
	static GtkToggleButton *cutmax = NULL;
	WORD hi, lo;
	gboolean cut_over;
	if (adjmin == NULL) {
		adjmax = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment1")); // scalemax
		adjmin = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment2")); // scalemin
		maxentry = GTK_ENTRY(gtk_builder_get_object(builder, "max_entry"));
		minentry = GTK_ENTRY(gtk_builder_get_object(builder, "min_entry"));
		cutmax = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "checkcut_max"));
	}
	if (single_image_is_loaded() && com.cvport < com.uniq->nb_layers &&
			com.uniq->layers && com.seq.current != RESULT_IMAGE) {
		hi = com.uniq->layers[com.cvport].hi;
		lo = com.uniq->layers[com.cvport].lo;
		cut_over = com.uniq->layers[com.cvport].cut_over;
	}
	/* When com.seq.current == RESULT_IMAGE we take sequence values */
	else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers
			&& com.seq.layers) {
		hi = com.seq.layers[com.cvport].hi;
		lo = com.seq.layers[com.cvport].lo;
		cut_over = com.seq.layers[com.cvport].cut_over;
	} else
		return;	// there should be no other normal cases
	fprintf(stdout, _("setting ranges scalemin=%d, scalemax=%d\n"), lo, hi);
	WORD maxvalue = get_normalized_value(&gfit);
	gtk_adjustment_set_lower(adjmin, 0.0);
	gtk_adjustment_set_lower(adjmax, 0.0);
	gtk_adjustment_set_upper(adjmin, (gdouble) maxvalue);
	gtk_adjustment_set_upper(adjmax, (gdouble) maxvalue);
	gtk_adjustment_set_value(adjmin, (gdouble) lo);
	gtk_adjustment_set_value(adjmax, (gdouble) hi);
	g_snprintf(buffer, 6, "%u", hi);
	g_signal_handlers_block_by_func(maxentry, on_max_entry_changed, NULL);
	gtk_entry_set_text(maxentry, buffer);
	g_signal_handlers_unblock_by_func(maxentry, on_max_entry_changed, NULL);
	g_snprintf(buffer, 6, "%u", lo);
	g_signal_handlers_block_by_func(minentry, on_min_entry_changed, NULL);
	gtk_entry_set_text(minentry, buffer);
	g_signal_handlers_unblock_by_func(minentry, on_min_entry_changed, NULL);
	gtk_toggle_button_set_active(cutmax, cut_over);
}

/* Sets the display mode combo box to the value stored in the relevant struct.
 * The operation is purely graphical. */
void set_display_mode() {
	static GtkComboBox *modecombo = NULL;
	display_mode mode;

	if (!modecombo)
		modecombo = GTK_COMBO_BOX(lookup_widget("combodisplay"));

	if (single_image_is_loaded() && com.cvport < com.uniq->nb_layers && com.uniq->layers
			&& com.seq.current != RESULT_IMAGE)
		mode = com.uniq->layers[com.cvport].rendering_mode;
	else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers
			&& com.seq.layers)
		mode = com.seq.layers[com.cvport].rendering_mode;
	else
		return;

	g_signal_handlers_block_by_func(modecombo, on_combodisplay_changed, NULL);
	gtk_combo_box_set_active(modecombo, mode);
	g_signal_handlers_unblock_by_func(modecombo, on_combodisplay_changed, NULL);
}

/* fill the label indicating how many images are selected in the gray and
 * which one is the reference image, at the bottom of the main window */
int adjust_sellabel() {
	static GtkLabel *global_label = NULL;
	char bufferglobal[256];
	gchar *seq_basename = NULL;

	if (global_label == NULL) {
		global_label = GTK_LABEL(lookup_widget("labelseq"));
	}

	if (sequence_is_loaded()) {
		seq_basename = g_path_get_basename(com.seq.seqname);

		if (com.seq.reference_image != -1) {
			char format[150];
			if (com.seq.fixed <= 1) {
				g_snprintf(format, sizeof(format),
						_("<%%s.seq>: %%d images selected out of %%d, reference image is %%d"));
			} else {
				g_snprintf(format, sizeof(format),
						_("<%%s.seq>: %%d images selected out of %%d, reference image is %%.%dd"),
						com.seq.fixed);
			}
		}
		g_snprintf(bufferglobal, sizeof(bufferglobal), _("%s, %d images selected"),
				seq_basename, com.seq.selnum);
		//gtk_widget_set_sensitive(lookup_widget("goregister_button"), com.seq.selnum>0?TRUE:FALSE);
	} else {
		g_snprintf(bufferglobal, sizeof(bufferglobal), _("- none -"));
		gtk_widget_set_sensitive(lookup_widget("goregister_button"), FALSE);
	}

	gtk_label_set_text(global_label, bufferglobal);
	g_free(seq_basename);
	return 0;
}

void set_icon_entry(GtkEntry *entry, gchar *string) {
	const gchar *text = NULL;

	gtk_entry_set_icon_from_icon_name(entry, GTK_ENTRY_ICON_SECONDARY, string);
	if (string) {
		text = _("This sequence name already exists!! Please change the name before converting.");
	}
	gtk_entry_set_icon_tooltip_text (entry, GTK_ENTRY_ICON_SECONDARY, text);
}

void update_MenuItem() {
	gboolean is_a_single_image_loaded;		/* An image is loaded. Not a sequence or only the result of stacking process */
	gboolean is_a_singleRGB_image_loaded;	/* A RGB image is loaded. Not a sequence or only the result of stacking process */
	gboolean any_image_is_loaded;			/* Something is loaded. Single image or Sequence */

	is_a_singleRGB_image_loaded = isrgb(&gfit) && (!sequence_is_loaded()
			|| (sequence_is_loaded() && (com.seq.current == RESULT_IMAGE
					|| com.seq.current == SCALED_IMAGE)));

	is_a_single_image_loaded = single_image_is_loaded()	&& (!sequence_is_loaded()
			|| (sequence_is_loaded() && (com.seq.current == RESULT_IMAGE
					|| com.seq.current == SCALED_IMAGE)));

	any_image_is_loaded = single_image_is_loaded() || sequence_is_loaded();

	/* toolbar button */
	gtk_widget_set_sensitive(lookup_widget("GtkToolMainBar"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("header_undo_button"), is_undo_available());
	gtk_widget_set_sensitive(lookup_widget("header_redo_button"), is_redo_available());
	/* File Menu */
	gtk_widget_set_sensitive(lookup_widget("header_save_button"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("info_menu_headers"), any_image_is_loaded && gfit.header != NULL);
	gtk_widget_set_sensitive(lookup_widget("info_menu_informations"), is_a_single_image_loaded);

	/* Image processing Menu */
	gtk_widget_set_sensitive(lookup_widget("removegreen"), is_a_singleRGB_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_satu"), is_a_singleRGB_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menu_negative"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitemcalibration"), is_a_singleRGB_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitemphotometriccalibration"), is_a_singleRGB_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menu_channel_separation"), is_a_singleRGB_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menu_slpitcfa"), any_image_is_loaded && !isrgb(&gfit));
	gtk_widget_set_sensitive(lookup_widget("menuitem_histo"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_asinh"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_fixbanding"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_cosmetic"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_deconvolution"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_resample"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_rotation"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_rotation90"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_rotation270"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_mirrorx"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_mirrory"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_background_extraction"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_wavelets"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menu_wavelet_separation"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_medianfilter"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_rgradient"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_clahe"), is_a_single_image_loaded);

	/* Image information menu */
	gtk_widget_set_sensitive(lookup_widget("info_menu_noise_estimation"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("info_menu_statistics"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("info_menu_astrometry"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("info_menu_dynamic_psf"), any_image_is_loaded);
#ifdef HAVE_LIBCURL
	/* Updates check */
	gtk_widget_set_visible(lookup_widget("main_menu_updates"), TRUE);
	/* Astrometry tool */
	gtk_widget_set_visible(lookup_widget("info_menu_astrometry"), TRUE);
#endif
}

void sliders_mode_set_state(sliders_mode sliders) {
	GtkToggleButton *radiobutton;	// Must not be static
	char *str[] =
	{ "radiobutton_hilo", "radiobutton_minmax", "radiobutton_user" };
	void *func[] = { on_radiobutton_hilo_toggled, on_radiobutton_minmax_toggled,
		on_radiobutton_user_toggled };

	radiobutton = GTK_TOGGLE_BUTTON(lookup_widget(str[sliders]));

	g_signal_handlers_block_by_func(radiobutton, func[sliders], NULL);
	gtk_toggle_button_set_active(radiobutton, TRUE);
	g_signal_handlers_unblock_by_func(radiobutton, func[sliders], NULL);
}

/* When rendering settings are chained, they need to be copied to other layers
 * when modified on the current layer. This procedure does that. It can be
 * called whenever a value has changed, or when the chaning has been enabled,
 * to synchronize all data.
 * Current view port is com.cvport, and data is stored in layer_info structs
 * Synchronized data: hi and lo cursors, cut over box, rendering mode.
 * DOES NOT REMAP/REDRAW.
 *
 * from_GUI: TRUE if get values from the GUI, FALSE if get the values from structs.
 * Returns 1 if chained, 0 if not.
 */
int copy_rendering_settings_when_chained(gboolean from_GUI) {
	static GtkToggleButton *chainedbutton = NULL;
	static GtkRange *range_lo = NULL, *range_hi = NULL;
	static GtkComboBox *modecombo = NULL;
	static GtkToggleButton *cutmax = NULL;

	gboolean is_chained;
	display_mode mode;
	WORD lo, hi;
	gboolean cut_over;
	int i, nb_layers;
	layer_info *layers = NULL;

	if (!chainedbutton) {	// init widgets
		chainedbutton = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_chain"));
		modecombo = GTK_COMBO_BOX(lookup_widget("combodisplay"));
		range_lo = GTK_RANGE(gtk_builder_get_object(builder, "scalemin"));
		range_hi = GTK_RANGE(gtk_builder_get_object(builder, "scalemax"));
		cutmax = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "checkcut_max"));
	}

	is_chained = gtk_toggle_button_get_active(chainedbutton);
	if (com.cvport == RGB_VPORT && !is_chained) return 0;
	int cvport = com.cvport == RGB_VPORT ? 0 : com.cvport;

	if (single_image_is_loaded() &&
			cvport < com.uniq->nb_layers && com.uniq->layers &&
			com.seq.current != RESULT_IMAGE) {
		layers = com.uniq->layers;
		nb_layers = com.uniq->nb_layers;
	} else if (sequence_is_loaded() && cvport < com.seq.nb_layers
			&& com.seq.layers) {
		layers = com.seq.layers;
		nb_layers = com.seq.nb_layers;
	} else
		return 0;

	if (from_GUI) {
		int raw_mode = gtk_combo_box_get_active(modecombo);
		/* update values in the layer_info for cvport */
		layers[cvport].rendering_mode =
			raw_mode >= 0 ? raw_mode : NORMAL_DISPLAY;
		layers[cvport].lo = round_to_WORD(gtk_range_get_value(range_lo));
		layers[cvport].hi = round_to_WORD(gtk_range_get_value(range_hi));
		layers[cvport].cut_over = gtk_toggle_button_get_active(cutmax);
	}
	if (!is_chained)
		return 0;
	mode = layers[cvport].rendering_mode;
	lo = layers[cvport].lo;
	hi = layers[cvport].hi;
	cut_over = layers[cvport].cut_over;

	for (i = 0; i < nb_layers; i++) {
		if (i == cvport)
			continue;
		layers[i].rendering_mode = mode;
		layers[i].lo = lo;
		layers[i].hi = hi;
		layers[i].cut_over = cut_over;
	}

	return 1;
}

void update_prepro_interface(gboolean allow_debayer) {
	static GtkToggleButton *udark = NULL, *uoffset = NULL, *uflat = NULL,
			       *checkAutoEvaluate = NULL, *pp_debayer = NULL;
	if (udark == NULL) {
		udark = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "usedark_button"));
		uoffset = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "useoffset_button"));
		uflat = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "useflat_button"));
		checkAutoEvaluate = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "checkbutton_auto_evaluate"));
		pp_debayer = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "checkButton_pp_dem"));
	}

	gtk_widget_set_sensitive(lookup_widget("prepro_button"),
			(sequence_is_loaded() || single_image_is_loaded())
			&& (gtk_toggle_button_get_active(udark)
				|| gtk_toggle_button_get_active(uoffset)
				|| gtk_toggle_button_get_active(uflat)));
	gtk_widget_set_sensitive(lookup_widget("grid24"),
			gtk_toggle_button_get_active(udark));
	gtk_widget_set_sensitive(lookup_widget("checkDarkOptimize"),
			gtk_toggle_button_get_active(udark));
	gtk_widget_set_sensitive(lookup_widget("GtkBoxFlat"),
			gtk_toggle_button_get_active(uflat));
	gtk_widget_set_sensitive(lookup_widget("entry_flat_norm"),
			gtk_toggle_button_get_active(uflat)
			&& !gtk_toggle_button_get_active(checkAutoEvaluate));

	gtk_widget_set_sensitive(lookup_widget("checkButton_pp_dem"),
			allow_debayer && (gtk_toggle_button_get_active(udark)
				|| gtk_toggle_button_get_active(uoffset)
				|| gtk_toggle_button_get_active(uflat)));
	if (!allow_debayer) {
		gtk_toggle_button_set_active(pp_debayer, FALSE);
	}
}

void clear_sampling_setting_box() {
	GtkComboBox *binning = GTK_COMBO_BOX(
			gtk_builder_get_object(builder, "combobinning"));
	GtkEntry* focal_entry = GTK_ENTRY(lookup_widget("focal_entry"));
	GtkEntry* pitchX_entry = GTK_ENTRY(lookup_widget("pitchX_entry"));
	GtkEntry* pitchY_entry = GTK_ENTRY(lookup_widget("pitchY_entry"));

	gtk_entry_set_text(focal_entry, "");
	gtk_entry_set_text(pitchX_entry, "");
	gtk_entry_set_text(pitchY_entry, "");
	gtk_combo_box_set_active(binning, 0);
}

void update_libraw_and_debayer_interface() {
	/**********COLOR ADJUSTEMENT**************/
	com.raw_set.bright = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("Brightness_spinbutton")));
	com.raw_set.mul[0] = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("Red_spinbutton")));
	com.raw_set.mul[2] = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("Blue_spinbutton")));

	com.raw_set.auto_mul = gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_multipliers")));
	com.raw_set.user_black = gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_blackpoint")));

	/**************WHITE BALANCE**************/
	com.raw_set.use_camera_wb = (int) gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_cam")));
	com.raw_set.use_auto_wb = (int) gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_auto")));

	/********MATRIX INTERPOLATION**************/
	com.raw_set.user_qual = gtk_combo_box_get_active(
			GTK_COMBO_BOX(lookup_widget("combo_dcraw_inter")));

	/********GAMMA CORRECTION**************/
	if (gtk_toggle_button_get_active(
				GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm0")))) {
		/* Linear Gamma Curve */
		com.raw_set.gamm[0] = 1.0;
		com.raw_set.gamm[1] = 1.0;
	} else if (gtk_toggle_button_get_active(
				GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm1")))) {
		/* BT.709 Gamma curve */
		com.raw_set.gamm[0] = 2.222;
		com.raw_set.gamm[1] = 4.5;
	} else {
		/* sRGB Gamma curve */
		com.raw_set.gamm[0] = 2.40;
		com.raw_set.gamm[1] = 12.92;
	}
	/* We write in config file */
	/*************SER**********************/
	com.debayer.use_bayer_header = gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_SER_use_header")));
	com.debayer.compatibility = gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_debayer_compatibility")));
	com.debayer.stretch = gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("stretch_CFA_to16_button")));
	com.debayer.xbayeroff = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(lookup_widget("xbayeroff_spin")));
	com.debayer.ybayeroff = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(lookup_widget("ybayeroff_spin")));
	writeinitfile();
}

void update_photometry_interface() {
	com.phot_set.gain = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinGain")));
	com.phot_set.inner = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinInner")));
	com.phot_set.outer = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinOuter")));
	com.phot_set.minval = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinMinPhot")));
	com.phot_set.maxval = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinMaxPhot")));
	writeinitfile();
}

char *vport_number_to_name(int vport) {
	switch (vport) {
		case RED_VPORT:
			return strdup("red");
		case GREEN_VPORT:
			return strdup("green");
		case BLUE_VPORT:
			return strdup("blue");
		case RGB_VPORT:
			return strdup("rgb");
	}
	return NULL;
}

int match_drawing_area_widget(GtkWidget *drawing_area, gboolean allow_rgb) {
	/* could be done with a for i=0 loop, to get rid of these defines */
	if (drawing_area == com.vport[RED_VPORT])
		return RED_VPORT;
	if (drawing_area == com.vport[GREEN_VPORT])
		return GREEN_VPORT;
	else if (drawing_area == com.vport[BLUE_VPORT])
		return BLUE_VPORT;
	else if (allow_rgb && drawing_area == com.vport[RGB_VPORT])
		return RGB_VPORT;
	return -1;
}

void calculate_fwhm(GtkWidget *widget) {
	/* calculate and display FWHM */
	int layer = match_drawing_area_widget(widget, FALSE);
	if (layer != -1) {
		char buf[64], label_name[16];
		char *layer_name = vport_number_to_name(layer);
		GtkLabel *label;
		if (com.selection.w && com.selection.h) {// Now we don't care about the size of the sample. Minimization checks that
			if (com.selection.w < 300 && com.selection.h < 300) {
				double roundness;
				double fwhm_val;

				fwhm_val = psf_get_fwhm(&gfit, layer, &roundness);
				g_snprintf(buf, sizeof(buf), "fwhm = %.2f, r = %.2f", fwhm_val,
						roundness);
			} else
				g_snprintf(buf, sizeof(buf), _("fwhm: selection is too large"));
		} else {
			g_snprintf(buf, sizeof(buf), _("fwhm: no selection"));
		}
		g_snprintf(label_name, sizeof(label_name), "labelfwhm%s", layer_name);
		free(layer_name);
		label = GTK_LABEL(gtk_builder_get_object(builder, label_name));
		gtk_label_set_text(label, buf);
	}
}

/* displays the opened image file name in the layers window.
 * if a unique file is loaded, its details are used instead of any sequence data
 */
void display_filename() {
	GtkLabel *fn_label;
	int nb_layers;
	char *str, *filename;
	if (com.uniq) {	// unique image
		filename = com.uniq->filename;
		nb_layers = com.uniq->nb_layers;
	} else {	// sequence
		filename = malloc(256);
		seq_get_image_filename(&com.seq, com.seq.current, filename);
		nb_layers = com.seq.nb_layers;
	}
	fn_label = GTK_LABEL(gtk_builder_get_object(builder, "labelfilename_red"));
	str = g_strdup_printf(_("%s (channel 0)"), filename);
	gtk_label_set_text(fn_label, str);
	g_free(str);

	if (nb_layers == 3) {	//take in charge both sequence and single image
		fn_label = GTK_LABEL(
				gtk_builder_get_object(builder, "labelfilename_green"));
		str = g_strdup_printf(_("%s (channel 1)"), filename);
		gtk_label_set_text(fn_label, str);
		g_free(str);

		fn_label = GTK_LABEL(
				gtk_builder_get_object(builder, "labelfilename_blue"));
		str = g_strdup_printf(_("%s (channel 2)"), filename);
		gtk_label_set_text(fn_label, str);
		g_free(str);

	}
	if (!com.uniq) {
		free(filename);
	}
}

/* set available layers in the layer list of registration */
void set_layers_for_assign() {
	int i;
	if (!com.seq.layers)
		return;
	for (i = 0; i < com.seq.nb_layers; i++) {
		char layer[100];
		if (!com.seq.layers[i].name) {
			if (com.seq.nb_layers == 1) {
				com.seq.layers[i].name = strdup(
						_(predefined_layers_colors[i].name));
				com.seq.layers[i].wavelength =
					predefined_layers_colors[i].wavelength;
			} else if (com.seq.nb_layers == 3) {
				com.seq.layers[i].name = strdup(
						_(predefined_layers_colors[i + 1].name));
				com.seq.layers[i].wavelength =
					predefined_layers_colors[i + 1].wavelength;
			} else {
				com.seq.layers[i].name = strdup(_("Unassigned"));
				com.seq.layers[i].wavelength = -1.0;
			}
		}
		g_snprintf(layer, sizeof(layer), "%d: %s", i, com.seq.layers[i].name);
	}
}

/* updates the combo box of registration layers to reflect data availability */
void set_layers_for_registration() {
	static GtkComboBoxText *cbbt_layers = NULL;
	int i;
	int reminder;

	if (cbbt_layers == NULL)
		cbbt_layers = GTK_COMBO_BOX_TEXT(
				gtk_builder_get_object(builder, "comboboxreglayer"));
	reminder = gtk_combo_box_get_active(GTK_COMBO_BOX(cbbt_layers));

	gtk_combo_box_text_remove_all(cbbt_layers);
	for (i = 0; i < com.seq.nb_layers; i++) {
		char layer[100];
		if (com.seq.layers[i].name)
			g_snprintf(layer, sizeof(layer), "%d: %s", i,
					com.seq.layers[i].name);
		else
			g_snprintf(layer, sizeof(layer), _("%d: not affected yet"), i);
		if (com.seq.regparam[i]) {
			strcat(layer, " (*)");
			if (reminder == -1)	// set as default selection
				reminder = i;
		}
		gtk_combo_box_text_append_text(cbbt_layers, layer);
	}
	/* First initialization */
	if (reminder == -1) {
		if (com.seq.nb_layers == 3)
			gtk_combo_box_set_active(GTK_COMBO_BOX(cbbt_layers), 1);
		else
			gtk_combo_box_set_active(GTK_COMBO_BOX(cbbt_layers), 0);
	}
	/* Already initialized or default selection to channel with data */
	else
		gtk_combo_box_set_active(GTK_COMBO_BOX(cbbt_layers), reminder);
}

void show_data_dialog(char *text, char *title) {
	GtkTextView *tv = GTK_TEXT_VIEW(lookup_widget("data_txt"));
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(tv);
	GtkTextIter itDebut;
	GtkTextIter itFin;

	gtk_text_buffer_get_bounds(tbuf, &itDebut, &itFin);
	gtk_text_buffer_delete(tbuf, &itDebut, &itFin);
	gtk_text_buffer_set_text(tbuf, text, strlen(text));
	gtk_window_set_title(GTK_WINDOW(lookup_widget("data_dialog")), title);

	gtk_widget_show(lookup_widget("data_dialog"));
}

/**
 * Get the active window on toplevels
 * @return the GtkWindow activated
 */
GtkWindow *siril_get_active_window() {
	GtkWindow *win = NULL;
	GList *list, *l;

	list = gtk_window_list_toplevels();

	for (l = list; l; l = l->next) {
		if (gtk_window_is_active((GtkWindow *) l->data)) {
			win = (GtkWindow *) l->data;
			break;
		}
	}
	g_list_free(list);
	return win;
}

static void initialize_FITS_name_entries() {
	GtkEntry *moffset, *mdark, *mflat, *final_stack;
	gchar *str[4];
	gint i;

	moffset = GTK_ENTRY(lookup_widget("offsetname_entry"));
	mdark = GTK_ENTRY(lookup_widget("darkname_entry"));
	mflat = GTK_ENTRY(lookup_widget("flatname_entry"));
	final_stack = GTK_ENTRY(lookup_widget("entryresultfile"));

	str[0] = g_strdup_printf("master-offset%s", com.ext);
	str[1] = g_strdup_printf("master-dark%s", com.ext);
	str[2] = g_strdup_printf("master-flat%s", com.ext);
	str[3] = g_strdup_printf("stack_result%s", com.ext);

	gtk_entry_set_text(moffset, str[0]);
	gtk_entry_set_text(mdark, str[1]);
	gtk_entry_set_text(mflat, str[2]);
	gtk_entry_set_text(final_stack, str[3]);

	for (i = 0; i < 4; i++)
		g_free(str[i]);
}

/* when a sequence is loaded, the processing (stacking) output file name is
 * modified to include the name of the sequence */
void set_output_filename_to_sequence_name() {
	static GtkEntry *output_file = NULL;
	gchar *msg;
	if (!output_file)
		output_file = GTK_ENTRY(lookup_widget("entryresultfile"));
	if (!com.seq.seqname || *com.seq.seqname == '\0')
		return;
	msg = g_strdup_printf("%s%sstacked%s", com.seq.seqname,
			ends_with(com.seq.seqname, "_") ?
			"" : (ends_with(com.seq.seqname, "-") ? "" : "_"), com.ext);
	gtk_entry_set_text(output_file, msg);

	g_free(msg);
}

void adjust_refimage(int n) {
	static GtkWidget *ref_butt2 = NULL;
	if (ref_butt2 == NULL) {
		ref_butt2 = lookup_widget("refframe2");
	}

	g_signal_handlers_block_by_func(ref_butt2, on_ref_frame_toggled, NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ref_butt2), com.seq.reference_image == n);
	g_signal_handlers_unblock_by_func(ref_butt2, on_ref_frame_toggled, NULL);
}

void close_tab() {
	GtkNotebook* Color_Layers = GTK_NOTEBOOK(lookup_widget("notebook1"));
	GtkWidget* page;

	if (com.seq.nb_layers == 1 || gfit.naxes[2] == 1) {
		page = gtk_notebook_get_nth_page(Color_Layers, RGB_VPORT);
		gtk_widget_hide(page);
		page = gtk_notebook_get_nth_page(Color_Layers, GREEN_VPORT);
		gtk_widget_hide(page);
		page = gtk_notebook_get_nth_page(Color_Layers, BLUE_VPORT);
		gtk_widget_hide(page);
		page = gtk_notebook_get_nth_page(Color_Layers, RED_VPORT);
		gtk_notebook_set_tab_label_text(Color_Layers, page, _("B&W channel"));
	} else {
		page = gtk_notebook_get_nth_page(Color_Layers, RED_VPORT);
		gtk_notebook_set_tab_label_text(Color_Layers, page, _("Red channel"));
		page = gtk_notebook_get_nth_page(Color_Layers, GREEN_VPORT);
		gtk_widget_show(page);
		page = gtk_notebook_get_nth_page(Color_Layers, BLUE_VPORT);
		gtk_widget_show(page);
		page = gtk_notebook_get_nth_page(Color_Layers, RGB_VPORT);
		gtk_widget_show(page);
	}
}

void activate_tab(int vport) {
	GtkNotebook* notebook = GTK_NOTEBOOK(lookup_widget("notebook1"));
	gtk_notebook_set_current_page(notebook, vport);
	// com.cvport is set in the event handler for changed page
}

void control_window_switch_to_tab(main_tabs tab) {
	GtkNotebook* notebook = GTK_NOTEBOOK(lookup_widget("notebook2"));
	gtk_notebook_set_current_page(notebook, tab);
}

void update_statusbar_convert() {
	GtkLabel *status_label = GTK_LABEL(lookup_widget("statuslabel_convert"));

	int nb_files = count_converted_files();
	if (nb_files == 0)
		gtk_label_set_text(status_label, " ");
	else {
		int selected = count_selected_files();
		gchar *str, *total;
		if (nb_files == 1) {
			str = g_strdup_printf(_("%d file loaded"), nb_files);
		} else {
			str = g_strdup_printf(_("%d files loaded"), nb_files);
		}
		if (selected == 0) {
			total = g_strdup(str);
		} else if (selected == 1) {
			total = g_strdup_printf(_("%d file selected, %s"), selected, str);
		} else {
			total = g_strdup_printf(_("%d files selected, %s"), selected, str);
		}
		gtk_label_set_text(status_label, total);
		g_free(str);
		g_free(total);
	}
}

void update_spinCPU(int max) {
	static GtkSpinButton *spin_cpu = NULL;

	if (spin_cpu == NULL) {
		spin_cpu = GTK_SPIN_BUTTON(lookup_widget("spinCPU"));
	}
	if (max > 0) {
		gtk_spin_button_set_range (spin_cpu, 1, (gdouble) max);
	}
	gtk_spin_button_set_value (spin_cpu, (gdouble) com.max_thread);
}

/*****************************************************************************
 *             I N I T I A L I S A T I O N      F U N C T I O N S            *
 ****************************************************************************/

static void add_accelerator(GtkApplication *app, const gchar *action_name,
		const gchar *accel) {
	const gchar *vaccels[] = { accel, NULL };

	gtk_application_set_accels_for_action(app, action_name, vaccels);
}

static void load_accels() {
	GApplication *application = g_application_get_default();

	add_accelerator(GTK_APPLICATION(application), "app.quit", "<Primary>Q");
	add_accelerator(GTK_APPLICATION(application), "app.preferences", "<Primary>P");
	add_accelerator(GTK_APPLICATION(application), "app.open", "<Primary>O");
	add_accelerator(GTK_APPLICATION(application), "app.undo", "<Primary>Z");
#ifdef _WIN32
	add_accelerator(GTK_APPLICATION(application), "app.redo", "<Primary>Y");
#else
	add_accelerator(GTK_APPLICATION(application), "app.redo", "<Primary><Shift>Z");
#endif
	add_accelerator(GTK_APPLICATION(application), "app.save_as", "<Primary><Shift>S");
	add_accelerator(GTK_APPLICATION(application), "app.cwd", "<Primary>D");
}

/* Initialize the combobox when loading new single_image */
void initialize_display_mode() {
	static GtkComboBox *modecombo = NULL;
	static GtkToggleButton *chainedbutton = NULL;
	display_mode mode;
	int i, raw_mode;

	if (!modecombo)
		modecombo = GTK_COMBO_BOX(lookup_widget("combodisplay"));
	raw_mode = gtk_combo_box_get_active(modecombo);
	/* Check if never initialized. In this case the mode is set to linear */
	if (raw_mode == -1)
		mode = NORMAL_DISPLAY;
	else
		mode = raw_mode;
	/* The mode is applyed for each layer */
	if (single_image_is_loaded() && com.cvport < com.uniq->nb_layers
			&& com.seq.current != RESULT_IMAGE) {
		for (i = 0; i < com.uniq->nb_layers; i++)
			com.uniq->layers[i].rendering_mode = mode;
	} else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers) {
		for (i = 0; i < com.seq.nb_layers; i++)
			com.seq.layers[i].rendering_mode = mode;
	}
	/* In the case where the layer were unchained, we chaine it */
	if (!chainedbutton)
		chainedbutton = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_chain"));
	if (!gtk_toggle_button_get_active(chainedbutton)) {
		g_signal_handlers_block_by_func(chainedbutton, on_checkchain_toggled,
				NULL);
		gtk_toggle_button_set_active(chainedbutton, TRUE);
		g_signal_handlers_unblock_by_func(chainedbutton, on_checkchain_toggled,
				NULL);
	}
}

void set_GUI_CWD() {
	if (!com.wd)
		return;
	gchar *str;
	GtkLabel *label = GTK_LABEL(lookup_widget("labelcwd"));
	GtkHeaderBar *bar = GTK_HEADER_BAR(lookup_widget("headebar"));

	gtk_label_set_text(label, com.wd);

	str = g_strdup_printf("Siril-%s", VERSION);
	gtk_header_bar_set_title(bar , str);
	gtk_header_bar_set_subtitle(bar, com.wd);

	g_free(str);
}

void set_GUI_misc() {
	GtkToggleButton *ToggleButton;
	GtkSpinButton *memory_percent, *memory_amount;

	ToggleButton = GTK_TOGGLE_BUTTON(lookup_widget("miscAskQuit"));
	gtk_toggle_button_set_active(ToggleButton, com.dontShowConfirm);
	ToggleButton = GTK_TOGGLE_BUTTON(lookup_widget("rememberWindowsCheck"));
	gtk_toggle_button_set_active(ToggleButton, com.remember_windows);

	memory_percent = GTK_SPIN_BUTTON(lookup_widget("spinbutton_mem_ratio"));
	gtk_spin_button_set_value(memory_percent, com.stack.memory_ratio);
	memory_amount = GTK_SPIN_BUTTON(lookup_widget("spinbutton_mem_amount"));
	gtk_spin_button_set_value(memory_amount, com.stack.memory_amount);

	GtkToggleButton *modes[3] = { GTK_TOGGLE_BUTTON(lookup_widget("memfreeratio_radio")),
		GTK_TOGGLE_BUTTON(lookup_widget("memfixed_radio")),
		GTK_TOGGLE_BUTTON(lookup_widget("memunlimited_radio")) };
	gtk_toggle_button_set_active(modes[com.stack.mem_mode], TRUE);
}

/* size is in kiB */
void set_GUI_MEM(unsigned long size) {
	char *str;
	if (size != 0)
		str = g_strdup_printf(_("Mem: %ldMB"), size / 1024);
	else
		str = g_strdup(_("Mem: N/A"));
	set_label_text_from_main_thread("labelmem", str);
	g_free(str);
}

void set_GUI_DiskSpace(int64_t space) {
	gchar *str;
	if (space > 0) {
		gchar *mem = pretty_print_memory(space);
		str = g_strdup_printf(_("Disk Space: %s"), mem);
		g_free(mem);
	} else
		str = g_strdup(_("Disk Space: N/A"));
	set_label_text_from_main_thread("labelFreeSpace", str);
	g_free(str);
}

static void initialize_preprocessing() {
	GtkToggleButton *cfaButton, *eqButton;

	cfaButton = GTK_TOGGLE_BUTTON(lookup_widget("cosmCFACheck"));
	gtk_toggle_button_set_active(cfaButton, com.prepro_cfa);
	eqButton = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_equalize_cfa"));
	gtk_toggle_button_set_active(eqButton, com.prepro_equalize_cfa);

	update_prepro_interface(FALSE);
}

static void set_libraw_settings_menu_available(gboolean activate) {
	if (!com.script) {
		GtkNotebook *notebook = GTK_NOTEBOOK(lookup_widget("notebook3"));
		GtkWidget *widget = gtk_notebook_get_nth_page(notebook, 0);

		gtk_widget_set_visible(widget, activate);
	}
}

void set_GUI_CAMERA() {
	GtkComboBox *binning = GTK_COMBO_BOX(lookup_widget("combobinning"));

	if (gfit.focal_length) {
		gchar *focal = g_strdup_printf("%.3lf", gfit.focal_length);
		gtk_entry_set_text(GTK_ENTRY(lookup_widget("focal_entry")), focal);
		g_free(focal);
	}
	if (gfit.pixel_size_x) {
		gchar *pitchX = g_strdup_printf("%.2lf", gfit.pixel_size_x);
		gtk_entry_set_text(GTK_ENTRY(lookup_widget("pitchX_entry")), pitchX);
		g_free(pitchX);
	}
	if (gfit.pixel_size_y) {
		gchar *pitchY = g_strdup_printf("%.2lf", gfit.pixel_size_y);
		gtk_entry_set_text(GTK_ENTRY(lookup_widget("pitchY_entry")), pitchY);
		g_free(pitchY);
	}

	if (!gfit.binning_x || !gfit.binning_y) {
		gtk_combo_box_set_active(binning, 0);
	}
	/* squared binning */
	else if (gfit.binning_x == gfit.binning_y)
		gtk_combo_box_set_active(binning, (gint) gfit.binning_x - 1);
	else {
		short coeff =
			gfit.binning_x > gfit.binning_y ?
			gfit.binning_x / gfit.binning_y :
			gfit.binning_y / gfit.binning_x;
		switch (coeff) {
			case 2:
				gtk_combo_box_set_active(binning, 4);
				break;
			case 3:
				gtk_combo_box_set_active(binning, 5);
				break;
			default:
				siril_log_message(_("This binning is not handled yet\n"));
		}
	}
}

static void set_GUI_LIBRAW() {

	/**********COLOR ADJUSTEMENT**************/
	gtk_spin_button_set_value(
			GTK_SPIN_BUTTON(lookup_widget("Brightness_spinbutton")),
			com.raw_set.bright);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("Red_spinbutton")),
			com.raw_set.mul[0]);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("Blue_spinbutton")),
			com.raw_set.mul[2]);

	gtk_toggle_button_set_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_multipliers")),
			com.raw_set.auto_mul);
	gtk_toggle_button_set_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_blackpoint")),
			com.raw_set.user_black);

	/**************WHITE BALANCE**************/
	if (com.raw_set.use_camera_wb) {
		gtk_toggle_button_set_active(
				GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_cam")),
				com.raw_set.use_camera_wb);
	}

	if (com.raw_set.use_auto_wb) {
		gtk_toggle_button_set_active(
				GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_auto")),
				com.raw_set.use_auto_wb);
	}

	/********MATRIX INTERPOLATION**************/
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_dcraw_inter")),
			com.raw_set.user_qual);

	/********GAMMA CORRECTION**************/
	if (com.raw_set.gamm[0] == 1.0 && com.raw_set.gamm[1] == 1.0)
		gtk_toggle_button_set_active(
				GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm0")), TRUE);
	else if (com.raw_set.gamm[0] == 2.222 && com.raw_set.gamm[1] == 4.5)
		gtk_toggle_button_set_active(
				GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm1")), TRUE);
	else
		gtk_toggle_button_set_active(
				GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm2")), TRUE);

	/********** DEBAYER ******************/
	GtkComboBox *pattern = GTK_COMBO_BOX(lookup_widget("comboBayer_pattern"));
	GtkComboBox *inter = GTK_COMBO_BOX(lookup_widget("comboBayer_inter"));
	GtkToggleButton *compat = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_debayer_compatibility"));
	GtkToggleButton *use_header = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_SER_use_header"));
	GtkToggleButton *stretch_cfa = GTK_TOGGLE_BUTTON(lookup_widget("stretch_CFA_to16_button"));
	GtkToggleButton *demosaicingButton = GTK_TOGGLE_BUTTON(lookup_widget("demosaicingButton"));
	GtkSpinButton *xbayer_spin = GTK_SPIN_BUTTON(lookup_widget("xbayeroff_spin"));
	GtkSpinButton *ybayer_spin = GTK_SPIN_BUTTON(lookup_widget("ybayeroff_spin"));
	gtk_combo_box_set_active(pattern, com.debayer.bayer_pattern);
	gtk_combo_box_set_active(inter, com.debayer.bayer_inter);
	gtk_toggle_button_set_active(compat, com.debayer.compatibility);
	gtk_toggle_button_set_active(use_header, com.debayer.use_bayer_header);
	gtk_toggle_button_set_active(demosaicingButton,	com.debayer.open_debayer);
	gtk_toggle_button_set_active(stretch_cfa, com.debayer.stretch);
	gtk_spin_button_set_value(xbayer_spin, com.debayer.xbayeroff);
	gtk_spin_button_set_value(ybayer_spin, com.debayer.ybayeroff);
}

void set_GUI_photometry() {
	if (gfit.cvf > 0.0)
		com.phot_set.gain = gfit.cvf;
	if (com.phot_set.gain > 0.0) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinGain")),
				com.phot_set.gain);
	}
	if (com.phot_set.inner > 0.0) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinInner")),
				com.phot_set.inner);
	}
	if (com.phot_set.outer > 0.0) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinOuter")),
				com.phot_set.outer);
	}
	if (com.phot_set.minval > 0.0) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinMinPhot")),
				com.phot_set.minval);
	}
	if (com.phot_set.maxval > 0.0) {
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("spinMaxPhot")),
				com.phot_set.maxval);
	}
}

static void initialize_path_directory() {
	GtkFileChooser *swap_dir;

	swap_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_swap"));
	if (com.swap_dir && com.swap_dir[0] != '\0')
		gtk_file_chooser_set_filename (swap_dir, com.swap_dir);
	else
		gtk_file_chooser_set_filename (swap_dir, g_get_tmp_dir());
}

static void initialize_scrollbars() {
	int i;
	char *vport_names[] = { "r", "g", "b", "rgb" };
	char window_name[32];

	for (i = 0; i < sizeof(vport_names) / sizeof(char *); i++) {
		sprintf(window_name, "scrolledwindow%s", vport_names[i]);
		GtkScrolledWindow *win = GTK_SCROLLED_WINDOW(gtk_builder_get_object(builder, window_name));
		com.hadj[i] = gtk_scrolled_window_get_hadjustment(win);
		g_signal_connect(com.hadj[i], "value-changed",
				G_CALLBACK(scrollbars_hadjustment_changed_handler), NULL);
		com.vadj[i] = gtk_scrolled_window_get_vadjustment(win);
		g_signal_connect(com.vadj[i], "value-changed",
				G_CALLBACK(scrollbars_vadjustment_changed_handler), NULL);
	}
}

static GtkTargetEntry drop_types[] = {
	{ "text/uri-list", 0, 0 }
};

void initialize_all_GUI(gchar *supported_files) {
	/* initializing internal structures with widgets (drawing areas) */
	com.vport[RED_VPORT] = lookup_widget("drawingarear");
	com.vport[GREEN_VPORT] = lookup_widget("drawingareag");
	com.vport[BLUE_VPORT] = lookup_widget("drawingareab");
	com.vport[RGB_VPORT] = lookup_widget("drawingareargb");
	com.preview_area[0] = lookup_widget("drawingarea_preview1");
	com.preview_area[1] = lookup_widget("drawingarea_preview2");
	initialize_image_display();
	initialize_scrollbars();
	init_mouse();

	/* Keybord Shortcuts */
	load_accels();

	/* Select combo boxes that trigger some text display or other things */
	gtk_combo_box_set_active(GTK_COMBO_BOX(gtk_builder_get_object(builder, "comboboxstack_methods")), 0);

	GtkLabel *label_supported = GTK_LABEL(lookup_widget("label_supported_types"));
	gtk_label_set_text(label_supported, supported_files);

	adjust_sellabel();

	/* initialize theme */
	initialize_theme_GUI();

	/* initialize menu gui */
	update_MenuItem();
	initialize_script_menu();

	/* initialize command completion */
	init_completion_command();

	/* initialize preprocessing */
	initialize_preprocessing();

	/* initialize registration methods */
	initialize_registration_methods();

	/* initialize stacking methods */
	initialize_stacking_methods();

	/* register some callbacks */
	register_selection_update_callback(update_export_crop_label);

	/* initialization of the binning parameters */
	GtkComboBox *binning = GTK_COMBO_BOX(lookup_widget("combobinning"));
	gtk_combo_box_set_active(binning, 0);

	/* initialization of some paths */
	initialize_path_directory();

	/* initialization of default FITS extension */
	GtkComboBox *box = GTK_COMBO_BOX(lookup_widget("combobox_ext"));
	gtk_combo_box_set_active_id(box, com.ext);
	initialize_FITS_name_entries();

	initialize_log_tags();

	/* support for converting files by dragging onto the GtkTreeView */
	gtk_drag_dest_set(lookup_widget("treeview_convert"),
			GTK_DEST_DEFAULT_MOTION, drop_types, G_N_ELEMENTS(drop_types),
			GDK_ACTION_COPY);

	set_GUI_CWD();
	set_GUI_misc();
	set_GUI_photometry();
	init_peaker_GUI();
#ifdef HAVE_LIBRAW
	set_libraw_settings_menu_available(TRUE);	// enable libraw settings
	set_GUI_LIBRAW();
#else
	set_libraw_settings_menu_available(FALSE);	// disable libraw settings
#endif
	update_spinCPU(com.max_thread);
	update_used_memory();
}

/*****************************************************************************
 *      P U B L I C      C A L L B A C K      F U N C T I O N S              *
 ****************************************************************************/

void on_register_all_toggle(GtkToggleButton *togglebutton, gpointer user_data) {
	update_reg_interface(TRUE);
}

/* when the cursor moves, update the value displayed in the textbox and save it
 * in the related layer_info. Does not change display until cursor is released. */
void on_minscale_changed(GtkRange *range, gpointer user_data) {
	static GtkEntry *minentry = NULL;
	gchar *buffer;
	if (minentry == NULL)
		minentry = GTK_ENTRY(gtk_builder_get_object(builder, "min_entry"));

	if (single_image_is_loaded() && com.seq.current < RESULT_IMAGE) {
		int value = (int) gtk_range_get_value(range);
		if (com.cvport < com.uniq->nb_layers)
			com.uniq->layers[com.cvport].lo = value;
		buffer = g_strdup_printf("%u", value);
	} else if (sequence_is_loaded()) {
		int value = (int) gtk_range_get_value(range);
		if (com.cvport < com.seq.nb_layers)
			com.seq.layers[com.cvport].lo = value;
		buffer = g_strdup_printf("%u", value);

	} else return;
	g_signal_handlers_block_by_func(minentry, on_min_entry_changed, NULL);
	gtk_entry_set_text(minentry, buffer);
	g_signal_handlers_unblock_by_func(minentry, on_min_entry_changed, NULL);
	g_free(buffer);
}

/* when the cursor moves, update the value displayed in the textbox and save it
 * in the related layer_info. Does not change display until cursor is released. */
void on_maxscale_changed(GtkRange *range, gpointer user_data) {
	static GtkEntry *maxentry = NULL;
	gchar *buffer;
	if (maxentry == NULL)
		maxentry = GTK_ENTRY(gtk_builder_get_object(builder, "max_entry"));

	if (single_image_is_loaded() && com.seq.current < RESULT_IMAGE) {
		int value = (int) gtk_range_get_value(range);
		if (com.cvport < com.uniq->nb_layers)
			com.uniq->layers[com.cvport].hi = value;
		buffer = g_strdup_printf("%u", value);
	} else if (sequence_is_loaded()) {
		int value = (int) gtk_range_get_value(range);
		if (com.cvport < com.seq.nb_layers)
			com.seq.layers[com.cvport].hi = value;
		buffer = g_strdup_printf("%u", value);
	} else return;
	g_signal_handlers_block_by_func(maxentry, on_max_entry_changed, NULL);
	gtk_entry_set_text(maxentry, buffer);
	g_signal_handlers_unblock_by_func(maxentry, on_max_entry_changed, NULL);
	g_free(buffer);
}

gboolean on_minscale_release(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	if (com.sliders != USER) {
		com.sliders = USER;
		sliders_mode_set_state(com.sliders);
	}
	if (copy_rendering_settings_when_chained(TRUE))
		redraw(com.cvport, REMAP_ALL);
	else
		redraw(com.cvport, REMAP_ONLY);
	redraw_previews();
	return FALSE;
}

gboolean on_maxscale_release(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	if (com.sliders != USER) {
		com.sliders = USER;
		sliders_mode_set_state(com.sliders);
	}
	if (copy_rendering_settings_when_chained(TRUE))
		redraw(com.cvport, REMAP_ALL);
	else
		redraw(com.cvport, REMAP_ONLY);
	redraw_previews();
	return FALSE;
}

/* a checkcut checkbox was toggled. Update the layer_info and others if chained. */
void on_checkcut_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	if (copy_rendering_settings_when_chained(TRUE))
		redraw(com.cvport, REMAP_ALL);
	else
		redraw(com.cvport, REMAP_ONLY);
	redraw_previews();
}

void on_cosmEnabledCheck_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkWidget *CFA, *SigHot, *SigCold, *checkHot, *checkCold, *evaluateButton;
	gboolean is_active;

	CFA = lookup_widget("cosmCFACheck");
	SigHot = lookup_widget("spinSigCosmeHot");
	SigCold = lookup_widget("spinSigCosmeCold");
	checkHot = lookup_widget("checkSigHot");
	checkCold = lookup_widget("checkSigCold");
	evaluateButton = lookup_widget("GtkButtonEvaluateCC");

	is_active = gtk_toggle_button_get_active(button);

	gtk_widget_set_sensitive(CFA, is_active);
	gtk_widget_set_sensitive(SigHot, is_active);
	gtk_widget_set_sensitive(SigCold, is_active);
	gtk_widget_set_sensitive(checkHot, is_active);
	gtk_widget_set_sensitive(checkCold, is_active);
	gtk_widget_set_sensitive(evaluateButton, is_active);
}

void on_cosmCFACheck_toggled(GtkToggleButton *button, gpointer user_data) {
	com.prepro_cfa = gtk_toggle_button_get_active(button);
	writeinitfile();
}

void on_checkbutton_equalize_cfa_toggled(GtkToggleButton *button, gpointer user_data) {
	com.prepro_equalize_cfa = gtk_toggle_button_get_active(button);
	writeinitfile();
}

void on_GtkButtonEvaluateCC_clicked(GtkButton *button, gpointer user_data) {
	GtkEntry *entry;
	GtkLabel *label[2];
	GtkWidget *widget[2];
	const char *filename;
	gchar *str[2];
	double sig[2];
	long icold = 0L, ihot = 0L;
	double rate, total;
	fits fit = { 0 };

	set_cursor_waiting(TRUE);
	sig[0] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeColdBox")));
	sig[1] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeHotBox")));
	widget[0] = lookup_widget("GtkLabelColdCC");
	widget[1] = lookup_widget("GtkLabelHotCC");
	label[0] = GTK_LABEL(lookup_widget("GtkLabelColdCC"));
	label[1] = GTK_LABEL(lookup_widget("GtkLabelHotCC"));
	entry = GTK_ENTRY(lookup_widget("darkname_entry"));
	filename = gtk_entry_get_text(entry);
	if (readfits(filename, &fit, NULL)) {
		str[0] = g_markup_printf_escaped(_("<span foreground=\"red\">ERROR</span>"));
		str[1] = g_markup_printf_escaped(_("<span foreground=\"red\">ERROR</span>"));
		gtk_label_set_markup(label[0], str[0]);
		gtk_label_set_markup(label[1], str[1]);
		set_cursor_waiting(FALSE);
		return;
	}
	count_deviant_pixels(&fit, sig, &icold, &ihot);
	total = fit.rx * fit.ry;
	clearfits(&fit);
	rate = (double)icold / total;
	/* 1% of cold pixels seems to be a reasonable limit */
	if (rate > 0.01) {
		str[0] = g_markup_printf_escaped(_("<span foreground=\"red\">Cold: %ld px</span>"), icold);
		gtk_widget_set_tooltip_text(widget[0], _("This value may be to high. Please, consider to change sigma value or uncheck the box."));
	}
	else {
		str[0] = g_markup_printf_escaped(_("Cold: %ld px"), icold);
		gtk_widget_set_tooltip_text(widget[0], "");
	}
	gtk_label_set_markup(label[0], str[0]);

	rate = (double)ihot / total;
	/* 1% of hot pixels seems to be a reasonable limit */
	if (rate > 0.01) {
		str[1] = g_markup_printf_escaped(_("<span foreground=\"red\">Hot: %ld px</span>"), ihot);
		gtk_widget_set_tooltip_text(widget[1], _("This value may be to high. Please, consider to change sigma value or uncheck the box."));
	}
	else {
		str[1] = g_markup_printf_escaped(_("Hot: %ld px"), ihot);
		gtk_widget_set_tooltip_text(widget[1], "");
	}
	gtk_label_set_markup(label[1], str[1]);
	g_free(str[0]);
	g_free(str[1]);
	set_cursor_waiting(FALSE);
}

void on_info_menu_headers_clicked(GtkButton *button, gpointer user_data) {
	show_FITS_header(&gfit);
}

void on_apply_settings_button_clicked(GtkButton *button, gpointer user_data) {
	update_libraw_and_debayer_interface();
	update_photometry_interface();
	fill_script_paths_list();
	refresh_stars_list(com.stars);
	save_main_window_state();
	siril_close_dialog("settings_window");
}

void on_focal_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar* focal_entry = gtk_entry_get_text(GTK_ENTRY(editable));
	gfit.focal_length = atof(focal_entry);
	update_fwhm_units_ok();
}

void on_pitchX_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar* pitchX_entry = gtk_entry_get_text(GTK_ENTRY(editable));
	gfit.pixel_size_x = (float) atof(pitchX_entry);
	update_fwhm_units_ok();
}

void on_pitchY_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar* pitchY_entry = gtk_entry_get_text(GTK_ENTRY(editable));
	gfit.pixel_size_y = (float) atof(pitchY_entry);
	update_fwhm_units_ok();
}

void on_combobinning_changed(GtkComboBox *box, gpointer user_data) {
	gint index = gtk_combo_box_get_active(box);

	switch (index) {
		case 0:
		case 1:
		case 2:
		case 3:
			gfit.binning_x = gfit.binning_y = (short) index + 1;
			break;
		case 4:
			gfit.binning_x = 1;
			gfit.binning_x = 2;
			break;
		case 5:
			gfit.binning_x = 1;
			gfit.binning_y = 3;
			break;
		default:
			fprintf(stderr, "Should not happen\n");
	}
}

void on_info_menu_informations_clicked(GtkButton *button, gpointer user_data) {
	siril_open_dialog("file_information");
}

void on_file_information_close_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("file_information"));
}


void on_toggleButtonUnbinned_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkWidget *box;

	box = lookup_widget("box_binning_info");

	gfit.unbinned = gtk_toggle_button_get_active(button);
	gtk_widget_set_sensitive(box, gfit.unbinned);
}

void on_button_clear_sample_clicked(GtkButton *button, gpointer user_data) {
	clear_sampling_setting_box();
}

void on_comboBayer_pattern_changed(GtkComboBox* box, gpointer user_data) {
	com.debayer.bayer_pattern = gtk_combo_box_get_active(box);
}

void on_comboBayer_inter_changed(GtkComboBox* box, gpointer user_data) {
	com.debayer.bayer_inter = gtk_combo_box_get_active(box);
}

void on_checkbutton_cam_toggled(GtkButton *button, gpointer user_data) {
	GtkToggleButton *auto_button = GTK_TOGGLE_BUTTON(
			lookup_widget("checkbutton_auto"));
	GtkToggleButton *cam_button = GTK_TOGGLE_BUTTON(
			lookup_widget("checkbutton_cam"));

	if (gtk_toggle_button_get_active(auto_button)) {
		g_signal_handlers_block_by_func(auto_button,
				on_checkbutton_auto_toggled, NULL);
		gtk_toggle_button_set_active(auto_button, FALSE);
		g_signal_handlers_unblock_by_func(auto_button,
				on_checkbutton_auto_toggled, NULL);
		gtk_toggle_button_set_active(cam_button, TRUE);
	}
}

void on_checkbutton_auto_toggled(GtkButton *button, gpointer user_data) {
	GtkToggleButton *auto_button = GTK_TOGGLE_BUTTON(
			lookup_widget("checkbutton_auto"));
	GtkToggleButton *cam_button = GTK_TOGGLE_BUTTON(
			lookup_widget("checkbutton_cam"));

	if (gtk_toggle_button_get_active(cam_button)) {
		g_signal_handlers_block_by_func(cam_button, on_checkbutton_cam_toggled,
				NULL);
		gtk_toggle_button_set_active(cam_button, FALSE);
		g_signal_handlers_unblock_by_func(cam_button,
				on_checkbutton_cam_toggled, NULL);
		gtk_toggle_button_set_active(auto_button, TRUE);
	}
}

void on_checkbutton_auto_evaluate_toggled(GtkToggleButton *button,
		gpointer user_data) {
	GtkWidget *entry = lookup_widget("entry_flat_norm");

	gtk_widget_set_sensitive(entry, !gtk_toggle_button_get_active(button));
}

void on_checkbutton_multipliers_toggled(GtkToggleButton *button,
		gpointer user_data) {
	gboolean active;

	active = gtk_toggle_button_get_active(button);

	gtk_widget_set_sensitive(lookup_widget("hbox8"), !active);
	gtk_widget_set_sensitive(lookup_widget("hbox11"), !active);
	if (active) {
		gtk_spin_button_set_value(
				GTK_SPIN_BUTTON(lookup_widget("Red_spinbutton")), 1.0);
		gtk_spin_button_set_value(
				GTK_SPIN_BUTTON(lookup_widget("Blue_spinbutton")), 1.0);
	}
}

void on_filechooser_swap_file_set(GtkFileChooserButton *fileChooser, gpointer user_data) {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(fileChooser);
	gchar *dir;

	dir = gtk_file_chooser_get_filename (swap_dir);

	if (g_access (dir, W_OK)) {
		gchar *msg = siril_log_color_message(_("You don't have permission to write in this directory: %s\n"), "red", dir);
		siril_message_dialog( GTK_MESSAGE_ERROR, _("Error"), msg);
		gtk_file_chooser_set_filename(swap_dir, com.swap_dir);
		return;
	}

	g_free(com.swap_dir);
	com.swap_dir = dir;
	writeinitfile();
}

void on_button_reset_swap_clicked(GtkButton *button, gpointer user_data) {
	reset_swapdir();
}

void on_spinbutton_mem_ratio_value_changed(GtkSpinButton *button, gpointer user_data) {
	gdouble mem = gtk_spin_button_get_value(button);
	com.stack.memory_ratio = mem;
	writeinitfile();
}

void on_spinbutton_mem_amount_value_changed(GtkSpinButton *button, gpointer user_data) {
	gdouble mem = gtk_spin_button_get_value(button);
	com.stack.memory_amount = mem;
	writeinitfile();
}

void on_mem_radio_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkToggleButton *ratio = GTK_TOGGLE_BUTTON(lookup_widget("memfreeratio_radio")),
			*amount = GTK_TOGGLE_BUTTON(lookup_widget("memfixed_radio")),
			*unlimited = GTK_TOGGLE_BUTTON(lookup_widget("memunlimited_radio"));
	GtkWidget *ratio_spin = lookup_widget("spinbutton_mem_ratio"),
		  *amount_spin = lookup_widget("spinbutton_mem_amount");
	if (!gtk_toggle_button_get_active(togglebutton)) return;

	if (togglebutton == ratio) {
		com.stack.mem_mode = 0;
		gtk_widget_set_sensitive(ratio_spin, TRUE);
		gtk_widget_set_sensitive(amount_spin, FALSE);
	} else if (togglebutton == amount) {
		com.stack.mem_mode = 1;
		gtk_widget_set_sensitive(ratio_spin, FALSE);
		gtk_widget_set_sensitive(amount_spin, TRUE);
	} else if (togglebutton == unlimited) {
		com.stack.mem_mode = 2;
		gtk_widget_set_sensitive(ratio_spin, FALSE);
		gtk_widget_set_sensitive(amount_spin, FALSE);
	}
}

void on_combobox_ext_changed(GtkComboBox *box, gpointer user_data) {
	if (com.ext)
		free(com.ext);

	com.ext = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(box));
	writeinitfile();
	initialize_FITS_name_entries();
}

static rectangle get_window_position(GtkWindow *window) {
	gint x, y, w, h;
	rectangle rec = { 0 };

	gtk_window_get_position(window, &x, &y);
	gtk_window_get_size(window, &w, &h);
	rec.x = x;
	rec.y = y;
	rec.w = w;
	rec.h = h;
	return rec;
}

void save_main_window_state() {
	if (!com.script && com.remember_windows) {
		GtkWidget *main_w = lookup_widget("control_window");
		com.main_w_pos = get_window_position(GTK_WINDOW(main_w));
		com.is_maximized = gtk_window_is_maximized(GTK_WINDOW(main_w));
	}
}

void load_main_window_state() {
	GtkWidget *win;

	win = lookup_widget("control_window");
	int x = com.main_w_pos.x;
	int y = com.main_w_pos.y;
	int w = com.main_w_pos.w;
	int h = com.main_w_pos.h;
	if (com.remember_windows && w >= 0 && h >= 0) {
		if (com.is_maximized) {
			gtk_window_maximize(GTK_WINDOW(lookup_widget("control_window")));
		} else {
			gtk_window_move(GTK_WINDOW(win), x, y);
			gtk_window_resize(GTK_WINDOW(win), w, h);
		}
	}
}

gboolean on_control_window_configure_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {

	save_main_window_state();
	if (com.zoom_value == -1)
		adjust_vport_size_to_image();
	return FALSE;
}

gboolean on_control_window_window_state_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {

	save_main_window_state();
	return FALSE;
}

void gtk_main_quit() {
	writeinitfile();		// save settings (like window positions)
	close_sequence(FALSE);	// save unfinished business
	close_single_image();	// close the previous image and free resources
	g_slist_free_full(com.script_path, g_free);
	exit(EXIT_SUCCESS);
}

void siril_quit() {
	if (com.dontShowConfirm) {
		gtk_main_quit();
	}
	gboolean quit = siril_confirm_dialog(_("Closing application"), _("Are you sure you want to quit?"), TRUE);
	if (quit) {
		gtk_main_quit();
	} else {
		fprintf(stdout, "Staying on the application.\n");
	}
}

/* We give one signal event by toggle button to fix a bug. Without this solution
 * the toggled function was called 2 times
 * one call to select new button, one call to unselect previous one */
void on_radiobutton_minmax_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	if (gtk_toggle_button_get_active(togglebutton)) {
		com.sliders = MINMAX;
		init_layers_hi_and_lo_values(com.sliders);
		set_cutoff_sliders_values();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	}
}

void on_radiobutton_hilo_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	if (gtk_toggle_button_get_active(togglebutton)) {
		com.sliders = MIPSLOHI;
		init_layers_hi_and_lo_values(com.sliders);
		set_cutoff_sliders_values();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	}
}

void on_radiobutton_user_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	if (gtk_toggle_button_get_active(togglebutton)) {
		com.sliders = USER;
		init_layers_hi_and_lo_values(com.sliders);
		set_cutoff_sliders_values();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	}
}

void on_neg_button_clicked(GtkToolButton *button, gpointer user_data) {
	int tmp;
	static GtkToggleButton *chainedbutton = NULL;
	gboolean is_chained;

	set_cursor_waiting(TRUE);

	chainedbutton = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_chain"));
	is_chained = gtk_toggle_button_get_active(chainedbutton);

	/* swaps values of hi and lo and redraw */
	if (!is_chained) {
		if (single_image_is_loaded() && com.cvport < com.uniq->nb_layers
				&& com.seq.current != RESULT_IMAGE) {
			tmp = com.uniq->layers[com.cvport].hi;
			com.uniq->layers[com.cvport].hi = com.uniq->layers[com.cvport].lo;
			com.uniq->layers[com.cvport].lo = tmp;
		} else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers) {
			tmp = com.seq.layers[com.cvport].hi;
			com.seq.layers[com.cvport].hi = com.seq.layers[com.cvport].lo;
			com.seq.layers[com.cvport].lo = tmp;
		} else
			return;
		set_cutoff_sliders_values();
		redraw(com.cvport, REMAP_ONLY);	// sliders are only set for cvport
	} else {
		int i;
		if (single_image_is_loaded() && com.seq.current != RESULT_IMAGE) {
			for (i = 0; i < com.uniq->nb_layers; i++) {
				tmp = com.uniq->layers[i].hi;
				com.uniq->layers[i].hi = com.uniq->layers[i].lo;
				com.uniq->layers[i].lo = tmp;
			}
		} else if (sequence_is_loaded()) {
			for (i = 0; i < com.seq.nb_layers; i++) {
				tmp = com.seq.layers[i].hi;
				com.seq.layers[i].hi = com.seq.layers[i].lo;
				com.seq.layers[i].lo = tmp;
			}
		} else
			return;
		set_cutoff_sliders_values();
		redraw(com.cvport, REMAP_ALL);
	}
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void on_colormap_button_toggled(GtkToggleToolButton *togglebutton,
		gpointer user_data) {
	set_cursor_waiting(TRUE);
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

/* Callback for the display mode change */
void on_combodisplay_changed(GtkComboBox *widget, gpointer user_data) {
	if (copy_rendering_settings_when_chained(TRUE))
		redraw(com.cvport, REMAP_ALL);
	else
		redraw(com.cvport, REMAP_ONLY);
	redraw_previews();
}

void on_checkchain_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	if (copy_rendering_settings_when_chained(FALSE))
		redraw(com.cvport, REMAP_ALL);
}

void on_max_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar *txt;
	int value;

	txt = gtk_entry_get_text(GTK_ENTRY(editable));
	if (g_ascii_isalnum(txt[0])) {

		value = atoi(txt);

		if (com.sliders != USER) {
			com.sliders = USER;
			sliders_mode_set_state(com.sliders);
		}
		if (single_image_is_loaded() && com.cvport < com.uniq->nb_layers
				&& com.seq.current != RESULT_IMAGE)
			com.uniq->layers[com.cvport].hi = value;
		else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers)
			com.seq.layers[com.cvport].hi = value;
		else
			return;

		set_cutoff_sliders_values();

		if (copy_rendering_settings_when_chained(FALSE))
			redraw(com.cvport, REMAP_ALL);
		else
			redraw(com.cvport, REMAP_ONLY);
		redraw_previews();
	}
}

gboolean on_max_entry_focus_out_event(GtkWidget *widget, gpointer user_data) {
	const gchar *txt;
	int len, i;
	gboolean isalnum = TRUE;

	txt = gtk_entry_get_text(GTK_ENTRY(widget));
	len = gtk_entry_get_text_length(GTK_ENTRY(widget));
	for (i = 0; i < len; i++)
		if (!g_ascii_isalnum(txt[i])) {
			isalnum = FALSE;
			break;
		}
	if (isalnum == FALSE || len == 0)
		gtk_entry_set_text(GTK_ENTRY(widget), "65535");
	return FALSE;
}

gboolean on_min_entry_focus_out_event(GtkWidget *widget, gpointer user_data) {
	const gchar *txt;
	int len, i;
	gboolean isalnum = TRUE;

	txt = gtk_entry_get_text(GTK_ENTRY(widget));
	len = gtk_entry_get_text_length(GTK_ENTRY(widget));
	for (i = 0; i < len; i++)
		if (!g_ascii_isalnum(txt[i])) {
			isalnum = FALSE;
			break;
		}
	if (isalnum == FALSE || len == 0)
		gtk_entry_set_text(GTK_ENTRY(widget), "0");
	return FALSE;
}

void on_min_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar *txt;
	int value;

	txt = gtk_entry_get_text(GTK_ENTRY(editable));
	if (g_ascii_isalnum(txt[0])) {

		value = atoi(txt);

		if (com.sliders != USER) {
			com.sliders = USER;
			sliders_mode_set_state(com.sliders);
		}
		if (single_image_is_loaded() && com.cvport < com.uniq->nb_layers
				&& com.seq.current != RESULT_IMAGE)
			com.uniq->layers[com.cvport].lo = value;
		else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers)
			com.seq.layers[com.cvport].lo = value;
		else
			return;
		set_cutoff_sliders_values();
		if (copy_rendering_settings_when_chained(FALSE))
			redraw(com.cvport, REMAP_ALL);
		else
			redraw(com.cvport, REMAP_ONLY);
		redraw_previews();
	}
}

static void test_for_master_files(struct preprocessing_data *args) {
	GtkToggleButton *tbutton;
	GtkEntry *entry;

	tbutton = GTK_TOGGLE_BUTTON(lookup_widget("useoffset_button"));
	if (gtk_toggle_button_get_active(tbutton)) {
		const char *filename;
		entry = GTK_ENTRY(lookup_widget("offsetname_entry"));
		filename = gtk_entry_get_text(entry);
		if (filename[0] == '\0') {
			gtk_toggle_button_set_active(tbutton, FALSE);
		} else {
			const char *error = NULL;
			set_progress_bar_data(_("Opening offset image..."), PROGRESS_NONE);
			args->bias = calloc(1, sizeof(fits));
			if (!readfits(filename, args->bias, NULL)) {
				if (args->bias->naxes[2] != gfit.naxes[2]) {
					error = _("NOT USING OFFSET: number of channels is different");
				} else if (args->bias->naxes[0] != gfit.naxes[0] ||
						args->bias->naxes[1] != gfit.naxes[1]) {
					error = _("NOT USING OFFSET: image dimensions are different");
				} else {
					args->use_bias = TRUE;
				}

			} else error = _("NOT USING OFFSET: cannot open the file\n");
			if (error) {
				siril_log_message("%s\n", error);
				set_progress_bar_data(error, PROGRESS_DONE);
				free(args->bias);
				gtk_entry_set_text(entry, "");
				args->use_bias = FALSE;
			}
		}
	}

	tbutton = GTK_TOGGLE_BUTTON(lookup_widget("usedark_button"));
	if (gtk_toggle_button_get_active(tbutton)) {
		const char *filename;
		entry = GTK_ENTRY(lookup_widget("darkname_entry"));
		filename = gtk_entry_get_text(entry);
		if (filename[0] == '\0') {
			gtk_toggle_button_set_active(tbutton, FALSE);
		} else {
			const char *error = NULL;
			set_progress_bar_data(_("Opening dark image..."), PROGRESS_NONE);
			args->dark = calloc(1, sizeof(fits));
			if (!readfits(filename, args->dark, NULL)) {
				if (args->dark->naxes[2] != gfit.naxes[2]) {
					error = _("NOT USING DARK: number of channels is different");
				} else if (args->dark->naxes[0] != gfit.naxes[0] ||
						args->dark->naxes[1] != gfit.naxes[1]) {
					error = _("NOT USING DARK: image dimensions are different");
				} else {
					args->use_dark = TRUE;
				}

			} else error = _("NOT USING DARK: cannot open the file\n");
			if (error) {
				siril_log_message("%s\n", error);
				set_progress_bar_data(error, PROGRESS_DONE);
				free(args->dark);
				gtk_entry_set_text(entry, "");
				args->use_dark = FALSE;
			}
		}
		// dark optimization
		tbutton = GTK_TOGGLE_BUTTON(lookup_widget("checkDarkOptimize"));
		args->use_dark_optim = gtk_toggle_button_get_active(tbutton);

		// cosmetic correction
		tbutton = GTK_TOGGLE_BUTTON(lookup_widget("cosmEnabledCheck"));
		args->use_cosmetic_correction = gtk_toggle_button_get_active(tbutton);
	}

	tbutton = GTK_TOGGLE_BUTTON(lookup_widget("useflat_button"));
	if (gtk_toggle_button_get_active(tbutton)) {
		const char *filename;
		entry = GTK_ENTRY(lookup_widget("flatname_entry"));
		filename = gtk_entry_get_text(entry);
		if (filename[0] == '\0') {
			gtk_toggle_button_set_active(tbutton, FALSE);
		} else {
			const char *error = NULL;
			set_progress_bar_data(_("Opening flat image..."), PROGRESS_NONE);
			args->flat = calloc(1, sizeof(fits));
			if (!readfits(filename, args->flat, NULL)) {
				if (args->flat->naxes[2] != gfit.naxes[2]) {
					error = _("NOT USING FLAT: number of channels is different");
				} else if (args->flat->naxes[0] != gfit.naxes[0] ||
						args->flat->naxes[1] != gfit.naxes[1]) {
					error = _("NOT USING FLAT: image dimensions are different");
				} else {
					args->use_flat = TRUE;
				}

			} else error = _("NOT USING FLAT: cannot open the file\n");
			if (error) {
				siril_log_message("%s\n", error);
				set_progress_bar_data(error, PROGRESS_DONE);
				free(args->flat);
				gtk_entry_set_text(entry, "");
				args->use_flat = FALSE;
			}
		}
	}
}

void on_prepro_button_clicked(GtkButton *button, gpointer user_data) {
	struct preprocessing_data *args;
	GtkEntry *entry;
	GtkWidget *autobutton;
	GtkToggleButton *CFA, *debayer, *equalize_cfa, *compatibility, *stretch_cfa;
	GtkSpinButton *sigHot, *sigCold;

	if (!single_image_is_loaded() && !sequence_is_loaded())
		return;

	if (!single_image_is_loaded() && get_thread_run()) {
		siril_log_message(_("Another task is already in "
					"progress, ignoring new request.\n"));
		return;
	}

	args = calloc(1, sizeof(struct preprocessing_data));
	test_for_master_files(args);
	siril_log_color_message(_("Preprocessing...\n"), "red");
	gettimeofday(&args->t_start, NULL);
	args->autolevel = TRUE;
	args->normalisation = 1.0f;	// will be updated anyway

	// set output filename (preprocessed file name prefix)
	entry = GTK_ENTRY(lookup_widget("preproseqname_entry"));
	args->ppprefix = gtk_entry_get_text(entry);

	autobutton = lookup_widget("checkbutton_auto_evaluate");
	args->autolevel = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(autobutton));
	if (!args->autolevel) {
		GtkEntry *norm_entry = GTK_ENTRY(lookup_widget("entry_flat_norm"));
		args->normalisation = atof(gtk_entry_get_text(norm_entry));
	}

	CFA = GTK_TOGGLE_BUTTON(lookup_widget("cosmCFACheck"));
	debayer = GTK_TOGGLE_BUTTON(lookup_widget("checkButton_pp_dem"));
	compatibility = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_debayer_compatibility"));
	stretch_cfa = GTK_TOGGLE_BUTTON(lookup_widget("stretch_CFA_to16_button"));
	equalize_cfa = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_equalize_cfa"));
	sigHot = GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeHot"));
	sigCold = GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeCold"));

	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkSigCold"))))
		args->sigma[0] = gtk_spin_button_get_value(sigCold);
	else args->sigma[0] = -1.0;

	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkSigHot"))))
		args->sigma[1] = gtk_spin_button_get_value(sigHot);
	else args->sigma[1] = -1.0;

	args->is_cfa = gtk_toggle_button_get_active(CFA);
	args->compatibility = gtk_toggle_button_get_active(compatibility);
	args->debayer = gtk_toggle_button_get_active(debayer);
	args->stretch_cfa =  gtk_toggle_button_get_active(stretch_cfa);
	args->equalize_cfa = gtk_toggle_button_get_active(equalize_cfa);

	/****/

	if (sequence_is_loaded()) {
		args->is_sequence = TRUE;
		args->seq = &com.seq;
		set_cursor_waiting(TRUE);
		control_window_switch_to_tab(OUTPUT_LOGS);
		start_sequence_preprocessing(args, FALSE);
	} else {
		int retval;
		args->is_sequence = FALSE;
		set_cursor_waiting(TRUE);
		control_window_switch_to_tab(OUTPUT_LOGS);

		retval = preprocess_single_image(args);

		free(args);

		if (retval)
			set_progress_bar_data(_("Error in preprocessing."), PROGRESS_NONE);
		else {
			set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
			invalidate_gfit_histogram();
			open_single_image_from_gfit();
		}
		set_cursor_waiting(FALSE);
	}
}

void on_regTranslationOnly_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkWidget *Algo, *Prefix;

	Algo = lookup_widget("ComboBoxRegInter");
	Prefix = lookup_widget("regseqname_entry");

	gtk_widget_set_sensitive(Algo, !gtk_toggle_button_get_active(togglebutton));
	gtk_widget_set_sensitive(Prefix, !gtk_toggle_button_get_active(togglebutton));
}

void on_seqproc_entry_changed(GtkComboBox *widget, gpointer user_data) {
	gchar *name = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(widget));
	if (name && name[0] != '\0') {
		gchar *msg;
		gchar *type;

		set_cursor_waiting(TRUE);
		const char *ext = get_filename_ext(name);
		if (!strcmp(ext, "ser")) {
			name[strlen(name) - 1] = 'q';
			type = " SER";
#ifdef HAVE_FFMS2
		} else if (!check_for_film_extensions(ext)) {
			int len = strlen(ext);
			strncpy(name + strlen(name) - len - 1, "seq", len + 1);
			type = " AVI";
#endif
		} else
			type = "";
		msg = g_strdup_printf(_("Selected %s sequence %s..."), type, name);
		set_progress_bar_data(msg, PROGRESS_DONE);
		set_seq(name);
		set_cursor_waiting(FALSE);
		set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
		g_free(msg);
	}
	g_free(name);
}

/* signal handler for the gray window layer change */
void on_notebook1_switch_page(GtkNotebook *notebook, GtkWidget *page,
		guint page_num, gpointer user_data) {
	com.cvport = page_num;
	set_cutoff_sliders_values();// load the previous known values for sliders
	set_display_mode();		// change the mode in the combo box if needed
	redraw(com.cvport, REMAP_ONLY);
	calculate_fwhm(com.vport[com.cvport]);
}

struct checkSeq_filter_data {
	int force;
	int retvalue;
	GtkToggleButton *forceButton;
};

static gboolean end_checkSeq(gpointer p) {
	struct checkSeq_filter_data *args = (struct checkSeq_filter_data *) p;
	stop_processing_thread();

	/* it's better to uncheck the force button each time it is used */
	if (args->force)
		gtk_toggle_button_set_active(args->forceButton, FALSE);
	if (args->retvalue)
		update_sequences_list(NULL);
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
	set_cursor_waiting(FALSE);
	free(args);
	update_used_memory();

	return FALSE;
}

static gpointer checkSeq(gpointer p) {
	struct checkSeq_filter_data *args = (struct checkSeq_filter_data *) p;

	if (!check_seq(args->force))
		args->retvalue = 1;
	siril_add_idle(end_checkSeq, args);
	return GINT_TO_POINTER(0);
}

void on_checkseqbutton_clicked(GtkButton *button, gpointer user_data) {
	static GtkToggleButton *forceButton = NULL;
	int force;

	if (forceButton == NULL) {
		forceButton = GTK_TOGGLE_BUTTON(lookup_widget("checkforceseq"));
	}
	force = gtk_toggle_button_get_active(forceButton);

	if (get_thread_run()) {
		siril_log_message(_("Another task is already "
					"in progress, ignoring new request.\n"));
		return;
	}

	set_cursor_waiting(TRUE);
	set_progress_bar_data(_("Searching for sequences in "
				"the current working directory..."), PROGRESS_NONE);

	struct checkSeq_filter_data *args = malloc(sizeof(struct checkSeq_filter_data));

	args->force = force;
	args->forceButton = forceButton;
	args->retvalue = 0;
	set_cursor_waiting(TRUE);
	start_in_new_thread(checkSeq, args);
}

void on_confirmDontShowButton_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {

	com.dontShowConfirm = gtk_toggle_button_get_active(togglebutton);
	set_GUI_misc();
	writeinitfile();
}

void on_button_data_ok_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("data_dialog"));
}

//void on_combozoom_changed(GtkComboBox *widget, gpointer user_data) {
//	gint active = gtk_combo_box_get_active(widget);
//	switch (active) {
//		case 0: /* 16:1 */
//			com.zoom_value = 16.;
//			break;
//		case 1: /* 8:1 */
//			com.zoom_value = 8.;
//			break;
//		case 2: /* 4:1 */
//			com.zoom_value = 4.;
//			break;
//		case 3: /* 2:1 */
//			com.zoom_value = 2.;
//			break;
//		case -1:
//		case 4: /* 1:1 */
//			com.zoom_value = 1.;
//			break;
//		case 5: /* 1:2 */
//			com.zoom_value = .5;
//			break;
//		case 6: /* 1:4 */
//			com.zoom_value = .25;
//			break;
//		case 7: /* 1:8 */
//			com.zoom_value = .125;
//			break;
//		case 8: /* fit to window */
//			com.zoom_value = -1.;
//			break;
//	}
//	fprintf(stdout, "zoom is now %f\n", com.zoom_value);
//	adjust_vport_size_to_image();
//	redraw(com.cvport, REMAP_NONE);
//}

void on_comboboxreglayer_changed(GtkComboBox *widget, gpointer user_data) {
	if (gtk_combo_box_get_active(widget) == -1)
		return;
	free_reference_image();
	update_stack_interface(TRUE);
}

void scrollbars_hadjustment_changed_handler(GtkAdjustment *adjustment,
		gpointer user_data) {
	int i;
	double value = gtk_adjustment_get_value(adjustment);

	for (i = 0; i < MAXVPORT; i++) {
		if (com.hadj[i] != adjustment) {
			gtk_adjustment_set_value(com.hadj[i], value);
		}
	}
}

void scrollbars_vadjustment_changed_handler(GtkAdjustment *adjustment,
		gpointer user_data) {
	int i;
	double value = gtk_adjustment_get_value(adjustment);
	for (i = 0; i < MAXVPORT; i++) {
		if (com.vadj[i] != adjustment) {
			gtk_adjustment_set_value(com.vadj[i], value);
		}
	}
}

void on_spinCPU_value_changed (GtkSpinButton *spinbutton, gpointer user_data) {
	com.max_thread = gtk_spin_button_get_value_as_int(spinbutton);
}

void on_rememberWindowsCheck_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	com.remember_windows = gtk_toggle_button_get_active(togglebutton);
}


void on_menu_rgb_align_select(GtkMenuItem *menuitem, gpointer user_data) {
	gboolean sel_is_drawn = ((com.selection.w > 0.0) && (com.selection.h > 0.0));

	gtk_widget_set_sensitive(lookup_widget("rgb_align_dft"), sel_is_drawn);
	gtk_widget_set_sensitive(lookup_widget("rgb_align_psf"), sel_is_drawn);
}

void on_rgb_align_dft_activate(GtkMenuItem *menuitem, gpointer user_data) {
	undo_save_state(&gfit, "Processing: RGB alignment (DFT)");
	rgb_align(1);
}

void on_rgb_align_psf_activate(GtkMenuItem *menuitem, gpointer user_data) {
	undo_save_state(&gfit, "Processing: RGB alignment (PSF)");
	rgb_align(0);
}

void on_gotoStacking_button_clicked(GtkButton *button, gpointer user_data) {
	control_window_switch_to_tab(STACKING);
}

gboolean on_right_panel_image_button_press_event(GtkWidget *button,
		GdkEventButton *event, gpointer user_data) {
	if (event->button == 1) {
		static gboolean panel_is_extended = TRUE;
		GtkImage *image = GTK_IMAGE(lookup_widget("right_panel_image"));
		GtkWidget *widget = gtk_paned_get_child2(GTK_PANED(lookup_widget("main_panel")));

		gtk_widget_set_visible(widget, !panel_is_extended);

		if (!panel_is_extended) {
			gtk_image_set_from_icon_name(image, "pan-end-symbolic",
					GTK_ICON_SIZE_BUTTON);
		} else {
			gtk_image_set_from_icon_name(image, "pan-start-symbolic",
					GTK_ICON_SIZE_BUTTON);
		}
		panel_is_extended = !panel_is_extended;
	}
	return TRUE;
}
