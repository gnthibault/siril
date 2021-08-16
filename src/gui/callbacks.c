/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2021 team free-astro (see more in AUTHORS file)
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

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/undo.h"
#include "core/command.h"
#include "core/command_line_processor.h"
#include "core/siril_app_dirs.h"
#include "core/siril_language.h"
#include "core/OS_utils.h"
#include "algos/siril_wcs.h"
#include "algos/star_finder.h"
#include "io/conversion.h"
#include "io/films.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "registration/registration.h"
#include "stacking/stacking.h"
#include "compositing/align_rgb.h"
#include "preferences.h"
#include "image_display.h"
#include "image_interactions.h"
#include "single_image.h"
#include "callbacks.h"

#include "algos/astrometry_solver.h"
#include "utils.h"
#include "plot.h"
#include "message_dialog.h"
#include "PSF_list.h"
#include "histogram.h"
#include "script_menu.h"
#include "progress_and_log.h"
#include "dialogs.h"
#include "fix_xtrans_af.h"
#include "siril_intro.h"
#include "siril_preview.h"
#include "siril-window.h"

layer_info predefined_layers_colors[] = {
	/* name, lambda, lo, hi, c/over, c/under, mode */
	{ N_("Luminance"), 0.0, 0, 0, FALSE, FALSE, LINEAR_DISPLAY }, // no color, undefined value is <0
	{ N_("Red"), 650.0, 0, 0, FALSE, FALSE, LINEAR_DISPLAY }, // approx. of the middle of the color
	{ N_("Green"), 530.0, 0, 0, FALSE, FALSE, LINEAR_DISPLAY },	// approx. of the middle of the color
	{ N_("Blue"), 450.0, 0, 0, FALSE, FALSE, LINEAR_DISPLAY }// approx. of the middle of the color
};

/*****************************************************************************
 *                    S T A T I C      F U N C T I O N S                     *
 ****************************************************************************/
/*
 * Memory label static functions
 */

void set_viewer_mode_widgets_sensitive(gboolean sensitive) {
	GtkWidget *scalemax = lookup_widget("scalemax");
	GtkWidget *scalemin = lookup_widget("scalemin");
	GtkWidget *entrymin = lookup_widget("min_entry");
	GtkWidget *entrymax = lookup_widget("max_entry");
	GtkWidget *minmax = lookup_widget("radiobutton_minmax");
	GtkWidget *hilo = lookup_widget("radiobutton_hilo");
	GtkWidget *user = lookup_widget("radiobutton_user");

	gtk_widget_set_sensitive(scalemax, sensitive);
	gtk_widget_set_sensitive(scalemin, sensitive);
	gtk_widget_set_sensitive(entrymin, sensitive);
	gtk_widget_set_sensitive(entrymax, sensitive);
	gtk_widget_set_sensitive(minmax, sensitive);
	gtk_widget_set_sensitive(hilo, sensitive);
	gtk_widget_set_sensitive(user, sensitive);
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
		update_theme_button("annotate_button", "astrometry_dark.svg");
		update_theme_button("photometry_button", "photometry_dark.svg");

		update_theme_button("rotate90_anticlock_button", "rotate-acw_dark.svg");
		update_theme_button("rotate90_clock_button", "rotate-cw_dark.svg");
		update_theme_button("mirrorx_button", "mirrorx_dark.svg");
		update_theme_button("mirrory_button", "mirrory_dark.svg");

		update_theme_button("process_starfinder_button", "starfinder_dark.svg");
		update_theme_button("tilt_button", "tilt_dark.svg");
		update_theme_button("sum_button", "sum_dark.svg");
		update_theme_button("export_button", "export_dark.svg");

		update_theme_button("histoToolAutoStretch", "mtf_dark.svg");
	} else {
		update_theme_button("annotate_button", "astrometry.svg");
		update_theme_button("photometry_button", "photometry.svg");

		update_theme_button("rotate90_anticlock_button", "rotate-acw.svg");
		update_theme_button("rotate90_clock_button", "rotate-cw.svg");
		update_theme_button("mirrorx_button", "mirrorx.svg");
		update_theme_button("mirrory_button", "mirrory.svg");

		update_theme_button("process_starfinder_button", "starfinder.svg");
		update_theme_button("tilt_button", "tilt.svg");
		update_theme_button("sum_button", "sum.svg");
		update_theme_button("export_button", "export.svg");

		update_theme_button("histoToolAutoStretch", "mtf.svg");
	}
}

void siril_set_theme(int active) {
	GtkSettings *settings = gtk_settings_get_default();
	g_object_set(settings, "gtk-application-prefer-dark-theme", active == 0, NULL);
	update_icons_to_theme(active == 0);
}

void on_combo_theme_changed(GtkComboBox *box, gpointer user_data) {
	int active = gtk_combo_box_get_active(box);

	siril_set_theme(active);
}

static void initialize_theme_GUI() {
	GtkComboBox *box;

	box = GTK_COMBO_BOX(lookup_widget("combo_theme"));

	g_signal_handlers_block_by_func(box, on_combo_theme_changed, NULL);
	gtk_combo_box_set_active(box, com.pref.combo_theme);
	g_signal_handlers_unblock_by_func(box, on_combo_theme_changed, NULL);
	update_icons_to_theme(com.pref.combo_theme == 0);
}

void load_prefered_theme(gint theme) {
	GtkSettings *settings;

	settings = gtk_settings_get_default();

	g_object_set(settings, "gtk-application-prefer-dark-theme", com.pref.combo_theme == 0, NULL);
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
	/* set max value for range according to number of bits of original image
	 * We should use gfit.bitpix for this, but it's currently always USHORT_IMG.
	 * Since 0.9.8 we have orig_bitpix, but it's not filled for SER and other images.
	 */

	max_val = (gfit.type == DATA_FLOAT ? USHRT_MAX_DOUBLE : (double)get_normalized_value(&gfit));
	siril_debug_print(_("Setting MAX value for cutoff sliders adjustments (%f)\n"), max_val);
	gtk_adjustment_set_lower(adj1, 0.0);
	gtk_adjustment_set_lower(adj2, 0.0);
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
	int vport;
	gboolean cut_over;
	if (adjmin == NULL) {
		adjmax = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment1")); // scalemax
		adjmin = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment2")); // scalemin
		maxentry = GTK_ENTRY(gtk_builder_get_object(builder, "max_entry"));
		minentry = GTK_ENTRY(gtk_builder_get_object(builder, "min_entry"));
		cutmax = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "checkcut_max"));
	}
	vport = com.cvport;
	if (com.cvport ==  RGB_VPORT) vport = GREEN_VPORT;
	if (single_image_is_loaded() && vport < com.uniq->nb_layers &&
			com.uniq->layers && com.seq.current != RESULT_IMAGE) {
		hi = com.uniq->layers[vport].hi;
		lo = com.uniq->layers[vport].lo;
		cut_over = com.uniq->layers[vport].cut_over;
	}
	/* When com.seq.current == RESULT_IMAGE we take sequence values */
	else if (sequence_is_loaded() && vport < com.seq.nb_layers
			&& com.seq.layers) {
		hi = com.seq.layers[vport].hi;
		lo = com.seq.layers[vport].lo;
		cut_over = com.seq.layers[vport].cut_over;
	} else
		return;	// there should be no other normal cases
	siril_debug_print(_("Setting ranges scalemin=%d, scalemax=%d\n"), lo, hi);
	gtk_adjustment_set_value(adjmin, (gdouble)lo);
	gtk_adjustment_set_value(adjmax, (gdouble)hi);
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

void on_menu_display_selection_done(GtkMenuShell *menushell, gpointer user_data) {
	GtkWidget *w = gtk_menu_get_active(GTK_MENU(menushell));
	const GtkWidget *const lbl = gtk_bin_get_child(GTK_BIN(w));
	const char *const text = gtk_label_get_text(GTK_LABEL(lbl));

	gtk_label_set_text((GtkLabel *)user_data, text);

	if (copy_rendering_settings_when_chained(TRUE))
		redraw(com.cvport, REMAP_ALL);
	else
		redraw(com.cvport, REMAP_ONLY);
	redraw_previews();
}

/* Sets the display mode combo box to the value stored in the relevant struct.
 * The operation is purely graphical. */
void set_display_mode() {
	static GtkMenu *display_menu = NULL;
	display_mode mode;
	int vport;

	if (!display_menu) {
		display_menu = GTK_MENU(lookup_widget("menu_display"));
	}

	vport = com.cvport;
	if (com.cvport == RGB_VPORT) vport = GREEN_VPORT;
	if (single_image_is_loaded() && vport < com.uniq->nb_layers && com.uniq->layers
			&& com.seq.current != RESULT_IMAGE)
		mode = com.uniq->layers[vport].rendering_mode;
	else if (sequence_is_loaded() && vport < com.seq.nb_layers
			&& com.seq.layers)
		mode = com.seq.layers[vport].rendering_mode;
	else
		return;

	g_signal_handlers_block_by_func(display_menu, on_menu_display_selection_done, NULL);
	gtk_menu_set_active(display_menu, mode);
	g_signal_handlers_unblock_by_func(display_menu, on_menu_display_selection_done, NULL);
}

/* fill the label indicating how many images are selected in the gray and
 * which one is the reference image, at the bottom of the main window */
void adjust_sellabel() {
	static GtkLabel *global_label = NULL, *label_name_of_seq = NULL;
	gchar *buffer_global, *buffer_title;

	if (global_label == NULL) {
		global_label = GTK_LABEL(lookup_widget("labelseq"));
		label_name_of_seq = GTK_LABEL(lookup_widget("label_name_of_seq"));
	}

	if (sequence_is_loaded()) {
		gchar *seq_basename = g_path_get_basename(com.seq.seqname);

		buffer_global = g_strdup_printf(_("%s, %d/%d images selected"),
				seq_basename, com.seq.selnum, com.seq.number);
		buffer_title = g_strdup(_("Sequence:"));

		g_free(seq_basename);
	} else if (single_image_is_loaded()) {
		gchar *filename = g_path_get_basename(com.uniq->filename);

		buffer_global = g_strdup_printf("%s", filename);
		buffer_title = g_strdup(_("Image:"));

		g_free(filename);
	} else {
		buffer_global = g_strdup(_("- none -"));
		buffer_title = NULL;

		gtk_widget_set_sensitive(lookup_widget("goregister_button"), FALSE);
	}

	gtk_label_set_text(label_name_of_seq, buffer_title);
	gtk_label_set_text(global_label, buffer_global);

	g_free(buffer_global);
	g_free(buffer_title);
}

void set_icon_entry(GtkEntry *entry, gchar *string) {

	gtk_entry_set_icon_from_icon_name(entry, GTK_ENTRY_ICON_SECONDARY, string);
	if (string) {
		gchar *text = g_strdup(_("This sequence name already exists!! "
				"Please change the name before converting."));
		gtk_entry_set_icon_tooltip_text (entry, GTK_ENTRY_ICON_SECONDARY, text);

		g_free(text);
	}
}

void update_MenuItem() {
	GtkApplicationWindow *app_win = GTK_APPLICATION_WINDOW(lookup_widget("control_window"));
	gboolean is_a_single_image_loaded;		/* An image is loaded. Not a sequence or only the result of stacking process */
	gboolean is_a_singleRGB_image_loaded;	/* A RGB image is loaded. Not a sequence or only the result of stacking process */
	gboolean any_image_is_loaded;			/* Something is loaded. Single image or Sequence */

	is_a_single_image_loaded = single_image_is_loaded()	&& (!sequence_is_loaded() || (sequence_is_loaded() && (com.seq.current == RESULT_IMAGE || com.seq.current == SCALED_IMAGE)));
	is_a_singleRGB_image_loaded = isrgb(&gfit) && single_image_is_loaded();
	any_image_is_loaded = single_image_is_loaded() || sequence_is_loaded();

	/* some toolbar buttons */
	gtk_widget_set_sensitive(lookup_widget("toolbarbox"), any_image_is_loaded);

	gboolean enable_button = any_image_is_loaded && has_wcs(&gfit);
	GAction *action_annotate = g_action_map_lookup_action(G_ACTION_MAP(app_win), "annotate-object");
	g_simple_action_set_enabled(G_SIMPLE_ACTION(action_annotate), enable_button);
	/* untoggle if disabled */
	if (!enable_button) {
		GVariant *state = g_action_get_state(action_annotate);
		if (g_variant_get_boolean(g_action_get_state(action_annotate))) {
			g_action_change_state(action_annotate, g_variant_new_boolean(FALSE));
		}
		g_variant_unref(state);
	}

	/* undo and redo */
	GAction *action_undo = g_action_map_lookup_action(G_ACTION_MAP(app_win), "undo");
	GAction *action_redo = g_action_map_lookup_action(G_ACTION_MAP(app_win), "redo");
	g_simple_action_set_enabled(G_SIMPLE_ACTION (action_undo), is_undo_available());
	g_simple_action_set_enabled(G_SIMPLE_ACTION (action_redo), is_redo_available());

	set_undo_redo_tooltip();

	/* save and save as */
	GAction *action_save = g_action_map_lookup_action (G_ACTION_MAP(gtk_window_get_application(GTK_WINDOW(app_win))), "save");
	GAction *action_save_as = g_action_map_lookup_action (G_ACTION_MAP(gtk_window_get_application(GTK_WINDOW(app_win))), "save-as");
	g_simple_action_set_enabled(G_SIMPLE_ACTION(action_save), is_a_single_image_loaded && com.uniq->fileexist);
	g_simple_action_set_enabled(G_SIMPLE_ACTION(action_save_as), any_image_is_loaded);

	/* popup menu */

	/* fits header */
	GAction *action_fits_header = g_action_map_lookup_action (G_ACTION_MAP(app_win), "fits-header");
	g_simple_action_set_enabled(G_SIMPLE_ACTION(action_fits_header), any_image_is_loaded && gfit.header != NULL);
	/* search object */
	GAction *action_search_objectr = g_action_map_lookup_action (G_ACTION_MAP(app_win), "search-object");
	g_simple_action_set_enabled(G_SIMPLE_ACTION(action_search_objectr), any_image_is_loaded && has_wcs(&gfit));
	/* selection is needed */
	siril_window_enable_if_selection_actions(app_win, com.selection.w && com.selection.h);
	/* selection and sequence is needed */
	siril_window_enable_if_selection_sequence_actions(app_win, com.selection.w && com.selection.h && sequence_is_loaded());

	/* Image processing Menu */

	/* single RGB image is needed */
	siril_window_enable_rgb_proc_actions(app_win, is_a_singleRGB_image_loaded);
	/* any image is needed */
	siril_window_enable_any_proc_actions(app_win, any_image_is_loaded);
	/* any mono image is needed */
	siril_window_enable_any_mono_proc_actions(app_win, any_image_is_loaded && !isrgb(&gfit));
	/* only single image needed */
	siril_window_enable_single_proc_actions(app_win, is_a_single_image_loaded);
	/* no images needed */
	siril_window_enable_none_proc_actions(app_win, TRUE);

	/* Image information menu */
	siril_window_enable_image_actions(app_win, any_image_is_loaded);
}

void sliders_mode_set_state(sliders_mode sliders) {
	GtkToggleButton *radiobutton;	// Must not be static
	gchar *str[] = { "radiobutton_hilo", "radiobutton_minmax", "radiobutton_user" };
	void *func[] = { on_radiobutton_hilo_toggled, on_radiobutton_minmax_toggled, on_radiobutton_user_toggled };

	radiobutton = GTK_TOGGLE_BUTTON(lookup_widget(str[sliders]));

	g_signal_handlers_block_by_func(radiobutton, func[sliders], NULL);
	gtk_toggle_button_set_active(radiobutton, TRUE);
	g_signal_handlers_unblock_by_func(radiobutton, func[sliders], NULL);
}

display_mode get_display_mode_from_menu() {
	GtkWidget *w = gtk_menu_get_active(GTK_MENU(lookup_widget("menu_display")));
	if (w == lookup_widget("log_item"))
		return LOG_DISPLAY;
	else if (w == lookup_widget("square_root_item"))
		return SQRT_DISPLAY;
	else if (w == lookup_widget("squared_item"))
		return SQUARED_DISPLAY;
	else if (w == lookup_widget("asinh_item"))
		return ASINH_DISPLAY;
	else if (w == lookup_widget("auto_item"))
		return STF_DISPLAY;
	else if (w == lookup_widget("histo_item"))
		return HISTEQ_DISPLAY;
	else
		return LINEAR_DISPLAY;
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
	static GtkToggleButton *cutmax = NULL;

	gboolean is_chained;
	display_mode mode;
	WORD lo, hi;
	gboolean cut_over;
	int i, nb_layers;
	layer_info *layers = NULL;

	if (!chainedbutton) {	// init widgets
		chainedbutton = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_chain"));
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
		int raw_mode = get_display_mode_from_menu();
		/* update values in the layer_info for cvport */
		layers[cvport].rendering_mode =
			raw_mode >= 0 ? raw_mode : LINEAR_DISPLAY;
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
			       *checkAutoEvaluate = NULL;
	static GtkWidget *prepro_button = NULL, *cosme_grid = NULL, *dark_optim = NULL;
       	static GtkWidget *equalize = NULL, *auto_eval = NULL, *flat_norm = NULL;
       	static GtkWidget *debayer = NULL, *fix_xtrans = NULL;
	static GtkComboBox *output_type = NULL;
	if (udark == NULL) {
		udark = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "usedark_button"));
		uoffset = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "useoffset_button"));
		uflat = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "useflat_button"));
		checkAutoEvaluate = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "checkbutton_auto_evaluate"));
		output_type = GTK_COMBO_BOX(
				gtk_builder_get_object(builder, "prepro_output_type_combo"));
		prepro_button = lookup_widget("prepro_button");
		cosme_grid = lookup_widget("grid24");
		dark_optim = lookup_widget("checkDarkOptimize");
		equalize = lookup_widget("checkbutton_equalize_cfa");
		auto_eval = lookup_widget("checkbutton_auto_evaluate");
		flat_norm = lookup_widget("entry_flat_norm");
		debayer = lookup_widget("checkButton_pp_dem");
		fix_xtrans = lookup_widget("fix_xtrans_af");
	}

	gtk_widget_set_sensitive(prepro_button,
			(sequence_is_loaded() || single_image_is_loaded())
			&& (gtk_toggle_button_get_active(udark)
				|| gtk_toggle_button_get_active(uoffset)
				|| gtk_toggle_button_get_active(uflat)));
	gtk_widget_set_sensitive(cosme_grid, gtk_toggle_button_get_active(udark));
	gtk_widget_set_sensitive(dark_optim, gtk_toggle_button_get_active(udark));
	gtk_widget_set_sensitive(equalize, gtk_toggle_button_get_active(uflat));
	gtk_widget_set_sensitive(auto_eval, gtk_toggle_button_get_active(uflat));
	gtk_widget_set_sensitive(flat_norm,
			gtk_toggle_button_get_active(uflat) &&
			!gtk_toggle_button_get_active(checkAutoEvaluate));

	gtk_widget_set_sensitive(debayer, allow_debayer && gtk_widget_get_sensitive(prepro_button));
	gtk_widget_set_sensitive(fix_xtrans, gtk_toggle_button_get_active(udark) || gtk_toggle_button_get_active(uoffset));

	gtk_widget_set_sensitive(GTK_WIDGET(output_type), sequence_is_loaded());
	int type = com.seq.type;
	if (com.seq.type < 0 || com.seq.type > SEQ_FITSEQ)
		type = SEQ_REGULAR;
	gtk_combo_box_set_active(output_type, type);
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

const char *vport_number_to_name(int vport) {
	switch (vport) {
		case RED_VPORT:
			return (_("red"));
		case GREEN_VPORT:
			return (_("green"));
		case BLUE_VPORT:
			return (_("blue"));
		case RGB_VPORT:
			return (_("rgb"));
	}
	return NULL;
}

const char *untranslated_vport_number_to_name(int vport) {
	switch (vport) {
		case RED_VPORT:
			return ("red");
		case GREEN_VPORT:
			return ("green");
		case BLUE_VPORT:
			return ("blue");
		case RGB_VPORT:
			return ("rgb");
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

void update_display_selection() {
	if (com.cvport == RGB_VPORT) return;
	static const gchar *label_selection[] = { "labelselection_red", "labelselection_green", "labelselection_blue", "labelselection_rgb"};
	static gchar selection_buffer[256] = { 0 };

	if (com.selection.w && com.selection.h) {
		g_sprintf(selection_buffer, _("W: %dpx H: %dpx ratio: %.4f"), com.selection.w, com.selection.h,
			(double)com.selection.w / (double)com.selection.h);
		gtk_label_set_text(GTK_LABEL(lookup_widget(label_selection[com.cvport])), selection_buffer);
	} else {
		gtk_label_set_text(GTK_LABEL(lookup_widget(label_selection[com.cvport])), "");
	}
}

void update_display_fwhm() {
	if (com.cvport == RGB_VPORT) return;
	static const gchar *label_fwhm[] = { "labelfwhm_red", "labelfwhm_green", "labelfwhm_blue", "labelfwhm_rgb"};
	static gchar fwhm_buffer[256] = { 0 };

	if (com.selection.w && com.selection.h) {// Now we don't care about the size of the sample. Minimization checks that
		if (com.selection.w < 300 && com.selection.h < 300) {
			double roundness;
			double fwhm_val = psf_get_fwhm(&gfit, com.cvport, &com.selection, &roundness);
			g_sprintf(fwhm_buffer, _("fwhm: %.2f px, r: %.2f"), fwhm_val, roundness);
		} else
			g_sprintf(fwhm_buffer, _("fwhm: N/A"));
	} else {
		g_sprintf(fwhm_buffer, _("fwhm: N/A"));
	}
	gtk_label_set_text(GTK_LABEL(lookup_widget(label_fwhm[com.cvport])), fwhm_buffer);
}

/* displays the opened image file name in the layers window.
 * if a unique file is loaded, its details are used instead of any sequence data
 */
void display_filename() {
	GtkLabel *fn_label;
	int nb_layers, vport;
	char *filename;
	if (com.uniq) {	// unique image
		filename = com.uniq->filename;
		nb_layers = com.uniq->nb_layers;
	} else {	// sequence
		filename = malloc(256);
		seq_get_image_filename(&com.seq, com.seq.current, filename);
		nb_layers = com.seq.nb_layers;
	}
	vport = nb_layers > 1 ? nb_layers + 1 : nb_layers;

	gchar *	base_name = g_path_get_basename(filename);
	for (int channel = 0; channel < vport; channel++) {
		gchar *c = g_strdup_printf("labelfilename_%s", untranslated_vport_number_to_name(channel));
		fn_label = GTK_LABEL(lookup_widget(c));
		gtk_label_set_text(fn_label, base_name);
		g_free(c);
	}

	if (!com.uniq) {
		free(filename);
	}
	g_free(base_name);
}

void on_precision_item_toggled(GtkCheckMenuItem *checkmenuitem, gpointer user_data) {
	if (!single_image_is_loaded()) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Cannot convert a sequence file"),
				_("A sequence file cannot be converted to 32 bits. This operation can only be done on a single file."));
	} else {
		size_t ndata = gfit.naxes[0] * gfit.naxes[1] * gfit.naxes[2];

		if (gfit.type == DATA_FLOAT) {
			gboolean convert = siril_confirm_dialog(_("Precision loss"),
					_("Converting the image from 32 bits to 16 bits may lead to a loss of numerical accuracy. "
							"Getting back to 32 bits will not recover this loss.\n"
							"Are you sure you want to convert your data?"), _("Convert to 16 bits"));
			if (convert) {
				if (is_preview_active())
					fit_replace_buffer(get_preview_gfit_backup(), float_buffer_to_ushort(gfit.fdata, ndata), DATA_USHORT);
				fit_replace_buffer(&gfit, float_buffer_to_ushort(gfit.fdata, ndata), DATA_USHORT);
				invalidate_gfit_histogram();
				update_gfit_histogram_if_needed();
				redraw(com.cvport, REMAP_ALL);
			}
		} else if (gfit.type == DATA_USHORT) {
			if (is_preview_active())
				fit_replace_buffer(get_preview_gfit_backup(), ushort_buffer_to_float(gfit.data, ndata), DATA_FLOAT);
			fit_replace_buffer(&gfit, ushort_buffer_to_float(gfit.data, ndata), DATA_FLOAT);
			invalidate_gfit_histogram();
			update_gfit_histogram_if_needed();
			redraw(com.cvport, REMAP_ALL);
		}
	}
	set_precision_switch();
}

void set_precision_switch() {
	if (!com.script) {
		GtkLabel *label = GTK_LABEL(lookup_widget("precision_button_name"));
		GtkCheckMenuItem *float_button = GTK_CHECK_MENU_ITEM(lookup_widget("32bits_item"));
		GtkCheckMenuItem *ushort_button = GTK_CHECK_MENU_ITEM(lookup_widget("16bits_item"));

		gtk_label_set_text(label, gfit.type == DATA_USHORT ? _("16 bits") : _("32 bits"));
		g_signal_handlers_block_by_func(float_button, on_precision_item_toggled, NULL);
		gtk_check_menu_item_set_active(float_button, gfit.type == DATA_FLOAT);
		gtk_check_menu_item_set_active(ushort_button, gfit.type == DATA_USHORT);
		g_signal_handlers_unblock_by_func(float_button,	on_precision_item_toggled, NULL);
	}
}

/* set available layers in the layer list of registration */
void set_layers_for_assign() {
	int i;
	if (!com.seq.layers)
		return;
	for (i = 0; i < com.seq.nb_layers; i++) {
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
		gchar *layer;
		if (com.seq.layers[i].name)
			layer = g_strdup_printf("%d: %s", i, com.seq.layers[i].name);
		else
			layer = g_strdup_printf(_("%d: not affected yet"), i);
		if (com.seq.regparam[i]) {
			str_append(&layer,  " (*)");
			if (reminder == -1)	// set as default selection
				reminder = i;
		}
		gtk_combo_box_text_append_text(cbbt_layers, layer);

		g_free(layer);
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

void show_data_dialog(char *text, char *title, gchar *parent, gchar *url) {
	GtkTextView *tv = GTK_TEXT_VIEW(lookup_widget("data_txt"));
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(tv);
	GtkTextIter itDebut;
	GtkTextIter itFin;
	GtkWindow *win = GTK_WINDOW(lookup_widget("data_dialog"));

	gtk_text_buffer_get_bounds(tbuf, &itDebut, &itFin);
	gtk_text_buffer_delete(tbuf, &itDebut, &itFin);
	gtk_text_buffer_set_text(tbuf, text, strlen(text));
	gtk_window_set_title(win, title);

	gtk_widget_show(lookup_widget("data_dialog"));
	gtk_widget_set_visible(lookup_widget("data_extra_button"), url != NULL);
	if (url) {
		gtk_link_button_set_uri((GtkLinkButton *)lookup_widget("data_extra_button"), url);
	}
	if (parent) {
		gtk_window_set_transient_for(win, GTK_WINDOW(lookup_widget(parent)));
	} else {
		gtk_window_set_transient_for(win, GTK_WINDOW(lookup_widget("control_window")));
	}
}

/**
 * Get the active window on toplevels
 * @return the GtkWindow activated
 */
GtkWindow *siril_get_active_window() {
	GtkWindow *win = NULL;
	GList *list = gtk_window_list_toplevels();

	for (GList *l = list; l; l = l->next) {
		if (gtk_window_is_active((GtkWindow *) l->data)) {
			win = (GtkWindow *) l->data;
			break;
		}
	}
	g_list_free(list);
	return win;
}

void initialize_FITS_name_entries() {
	GtkEntry *mbias, *mdark, *mflat, *final_stack;
	gchar *str[4] = { NULL };

	mbias = GTK_ENTRY(lookup_widget("offsetname_entry"));
	mdark = GTK_ENTRY(lookup_widget("darkname_entry"));
	mflat = GTK_ENTRY(lookup_widget("flatname_entry"));
	final_stack = GTK_ENTRY(lookup_widget("entryresultfile"));

	if (com.pref.prepro_bias_lib && (g_file_test(com.pref.prepro_bias_lib, G_FILE_TEST_EXISTS))) {
		if (com.pref.use_bias_lib) {
			str[0] = g_strdup_printf("%s", com.pref.prepro_bias_lib);
		}
	}

	if (com.pref.prepro_bias_synth) {
		if (com.pref.use_bias_synth) {
			str[0] = g_strdup_printf("%s", com.pref.prepro_bias_synth);
		}
	}

	if (com.pref.prepro_dark_lib && (g_file_test(com.pref.prepro_dark_lib, G_FILE_TEST_EXISTS))) {
		if (com.pref.use_dark_lib) {
			str[1] = g_strdup_printf("%s", com.pref.prepro_dark_lib);
		}
	}

	if (com.pref.prepro_flat_lib && (g_file_test(com.pref.prepro_flat_lib, G_FILE_TEST_EXISTS))) {
		if (com.pref.use_flat_lib) {
			str[2] = g_strdup_printf("%s", com.pref.prepro_flat_lib);
		}

	}
	str[3] = g_strdup_printf("stack_result%s", com.pref.ext);

	const char *t_bias = gtk_entry_get_text(mbias);
	if (!str[0] && t_bias[0] == '\0') {
		str[0] = g_strdup_printf("master-bias%s", com.pref.ext);
	}
	if (str[0]) {
		gtk_entry_set_text(mbias, str[0]);
	}

	const char *t_dark = gtk_entry_get_text(mdark);
	if (!str[1] && t_dark[0] == '\0') {
		str[1] = g_strdup_printf("master-dark%s", com.pref.ext);
	}
	if (str[1]) {
		gtk_entry_set_text(mdark, str[1]);
	}

	const char *t_flat = gtk_entry_get_text(mflat);
	if (!str[2] && t_flat[0] == '\0') {
		str[2] = g_strdup_printf("master-flat%s", com.pref.ext);
	}
	if (str[2]) {
		gtk_entry_set_text(mflat, str[2]);
	}

	const char *t_stack = gtk_entry_get_text(final_stack);
	if (t_stack[0] == '\0') {
		gtk_entry_set_text(final_stack, str[3]);
	}

	for (int i = 0; i < 4; i++)
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
			g_str_has_suffix(com.seq.seqname, "_") ?
			"" : (g_str_has_suffix(com.seq.seqname, "-") ? "" : "_"), com.pref.ext);
	gtk_entry_set_text(output_file, msg);

	g_free(msg);
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
		gtk_notebook_set_tab_label_text(Color_Layers, page, _("B&W"));
	} else {
		page = gtk_notebook_get_nth_page(Color_Layers, RED_VPORT);
		gtk_notebook_set_tab_label_text(Color_Layers, page, _("Red"));
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
	if (gtk_notebook_get_current_page(notebook) != vport)
		gtk_notebook_set_current_page(notebook, vport);
	// com.cvport is set in the event handler for changed page
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

static void load_accels() {
	static const gchar * const accelmap[] = {
		"app.quit",                   "<Primary>q", NULL,
		"app.preferences",            "<Primary>p", NULL,
		"app.open",                   "<Primary>o", NULL,
		"app.save",                   "<Primary>s", NULL,
		"app.save-as",                "<Primary><Shift>s", NULL,

		"win.undo",                   "<Primary>z", NULL,
		"win.redo",                   "<Primary><Shift>z", NULL,
		"win.close",                  "<Primary>w", NULL,
		"win.cwd",                    "<Primary>d", NULL,
		"win.full-screen",            "<Primary>f", NULL,

		"win.conversion",             "F1", NULL,
		"win.sequence",               "F2", NULL,
		"win.prepro",                 "F3", NULL,
		"win.registration",           "F4", NULL,
		"win.plot",                   "F5", NULL,
		"win.stacking",               "F6", NULL,
		"win.logs",                   "F7", NULL,

		"win.hide-show-toolbar",      "<Primary>T", NULL,

		"win.zoom-out",               "<Primary>minus", "<Primary>KP_Subtract", NULL,
		"win.zoom-in",                "<Primary>plus", "<Primary>KP_Add", NULL,
		"win.zoom-fit",               "<Primary>0", NULL,

		"win.search-object",          "<Primary>slash", NULL,
		"win.astrometry",             "<Primary><Shift>a", NULL,
		"win.pickstar",               "<Primary>space", NULL,
		"win.dyn-psf",                "<Primary>F6", NULL,

		"win.negative-processing",    "<Primary>i", NULL,
		"win.rotation90-processing",  "<Primary>Right", NULL,
		"win.rotation270-processing", "<Primary>Left", NULL,
		"win.mirrorx-processing",     "<Primary>Up", NULL,
		"win.mirrory-processing",     "<Primary>Down", NULL,

		NULL /* Terminating NULL */
	};

	set_accel_map(accelmap);
}

void set_accel_map(const gchar * const *accelmap) {
	GApplication *application = g_application_get_default();

	for (const gchar *const *it = accelmap; it[0]; it += g_strv_length((gchar**) it) + 1) {
		gtk_application_set_accels_for_action(GTK_APPLICATION(application), it[0], &it[1]);
	}
}

/* Initialize the combobox when loading new single_image */
void initialize_display_mode() {
	static GtkToggleButton *chainedbutton = NULL;
	display_mode mode;
	int i;

	int raw_mode = get_display_mode_from_menu();

	/* Check if never initialized. In this case the mode is set to linear */
	if (raw_mode == -1)
		mode = LINEAR_DISPLAY;
	else
		mode = raw_mode;
	/* The mode is applyed for each layer */
	if (single_image_is_loaded() && com.seq.current != RESULT_IMAGE) {
		for (i = 0; i < com.uniq->nb_layers; i++)
			com.uniq->layers[i].rendering_mode = mode;
	} else if (sequence_is_loaded()) {
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
	GtkHeaderBar *bar = GTK_HEADER_BAR(lookup_widget("headerbar"));

	gchar *str = g_strdup_printf("Siril-%s", VERSION);
	gtk_header_bar_set_title(bar , str);

	gchar *truncated_wd = siril_truncate_str(com.wd, 50);
	gtk_header_bar_set_subtitle(bar, truncated_wd);

	g_free(str);
	g_free(truncated_wd);
}

void set_GUI_misc() {
	GtkToggleButton *ToggleButton;
	GtkSpinButton *memory_percent, *memory_amount, *font_size;

	ToggleButton = GTK_TOGGLE_BUTTON(lookup_widget("miscAskQuit"));
	gtk_toggle_button_set_active(ToggleButton, com.pref.save.quit);
	ToggleButton = GTK_TOGGLE_BUTTON(lookup_widget("miscAskUpdateStartup"));
	gtk_toggle_button_set_active(ToggleButton, com.pref.check_update);
	ToggleButton = GTK_TOGGLE_BUTTON(lookup_widget("miscAskScript"));
	gtk_toggle_button_set_active(ToggleButton, com.pref.save.warn_script);
	ToggleButton = GTK_TOGGLE_BUTTON(lookup_widget("script_check_version"));
	gtk_toggle_button_set_active(ToggleButton, com.pref.script_check_requires);
	ToggleButton = GTK_TOGGLE_BUTTON(lookup_widget("show_preview_button"));
	gtk_toggle_button_set_active(ToggleButton, com.pref.show_thumbnails);
	GtkComboBox *thumb_box = GTK_COMBO_BOX(lookup_widget("thumbnails_box_size"));
	gtk_combo_box_set_active(thumb_box, com.pref.thumbnail_size == 256 ? 1: 0);
	ToggleButton = GTK_TOGGLE_BUTTON(lookup_widget("rememberWindowsCheck"));
	gtk_toggle_button_set_active(ToggleButton, com.pref.remember_windows);
	font_size = GTK_SPIN_BUTTON(lookup_widget("pref_fontsize"));
	gtk_spin_button_set_value(font_size, com.pref.font_scale);
	ToggleButton = GTK_TOGGLE_BUTTON(lookup_widget("pref_iconstyle"));
	gtk_toggle_button_set_active(ToggleButton, com.pref.icon_symbolic);

	memory_percent = GTK_SPIN_BUTTON(lookup_widget("spinbutton_mem_ratio"));
	gtk_spin_button_set_value(memory_percent, com.pref.stack.memory_ratio);
	memory_amount = GTK_SPIN_BUTTON(lookup_widget("spinbutton_mem_amount"));
	gtk_spin_button_set_value(memory_amount, com.pref.stack.memory_amount);

	GtkToggleButton *modes[2] = { GTK_TOGGLE_BUTTON(lookup_widget("memfreeratio_radio")),
		GTK_TOGGLE_BUTTON(lookup_widget("memfixed_radio"))};
	gtk_toggle_button_set_active(modes[com.pref.stack.mem_mode], TRUE);

	/* initialization of default FITS extension and type */
	GtkComboBox *combobox_type = GTK_COMBO_BOX(lookup_widget("combobox_type"));
	gtk_combo_box_set_active(combobox_type, com.pref.force_to_16bit ? 0 : 1);
	GtkComboBox *fit_ext = GTK_COMBO_BOX(lookup_widget("combobox_ext"));
	gtk_combo_box_set_active_id(fit_ext, com.pref.ext);
}

static void initialize_preprocessing() {
	GtkToggleButton *cfaButton, *eqButton, *xtransButton;

	cfaButton = GTK_TOGGLE_BUTTON(lookup_widget("cosmCFACheck"));
	gtk_toggle_button_set_active(cfaButton, com.pref.prepro_cfa);
	eqButton = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_equalize_cfa"));
	gtk_toggle_button_set_active(eqButton, com.pref.prepro_equalize_cfa);
	xtransButton = GTK_TOGGLE_BUTTON(lookup_widget("fix_xtrans_af"));
	gtk_toggle_button_set_active(xtransButton, com.pref.fix_xtrans);

	update_prepro_interface(FALSE);
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

static void initialize_scrollbars() {
	int i;
	char *vport_names[] = { "r", "g", "b", "rgb" };
	char *window_name;

	for (i = 0; i < sizeof(vport_names) / sizeof(char *); i++) {
		window_name = g_strdup_printf("scrolledwindow%s", vport_names[i]);
		GtkScrolledWindow *win = GTK_SCROLLED_WINDOW(gtk_builder_get_object(builder, window_name));
		com.hadj[i] = gtk_scrolled_window_get_hadjustment(win);
		g_signal_connect(com.hadj[i], "value-changed",
				G_CALLBACK(scrollbars_hadjustment_changed_handler), NULL);
		com.vadj[i] = gtk_scrolled_window_get_vadjustment(win);
		g_signal_connect(com.vadj[i], "value-changed",
				G_CALLBACK(scrollbars_vadjustment_changed_handler), NULL);

		g_free(window_name);
	}
}

static GtkTargetEntry drop_types[] = {
	{ "text/uri-list", 0, 0 }
};

static gboolean on_control_window_configure_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	save_main_window_state();
	return FALSE;
}

static gboolean on_control_window_window_state_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	save_main_window_state();
	return FALSE;
}

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

	/* populate language combo */
	siril_language_fill_combo(com.pref.combo_lang);

	/* Keybord Shortcuts */
	load_accels();

	/* Select combo boxes that trigger some text display or other things */
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("comboboxstack_methods")), 0);

	GtkLabel *label_supported = GTK_LABEL(lookup_widget("label_supported_types"));
	gtk_label_set_text(label_supported, supported_files);

	adjust_sellabel();

	/* initialize theme */
	initialize_theme_GUI();

	/* map all actions for main window */
	siril_window_map_actions(GTK_APPLICATION_WINDOW(lookup_widget("control_window")));

	/* initialize menu gui */
	update_MenuItem();
	initialize_script_menu();

	/* initialize command processor */
	init_command();

	/* initialize preprocessing */
	initialize_preprocessing();

	/* initialize registration methods */
	initialize_registration_methods();

	/* initialize stacking methods */
	initialize_stacking_methods();

	/* set all preferences button in the dialog */
	set_preferences_dialog_from_global();

	/* set focal and pixel pitch */
	set_focal_and_pixel_pitch();

	initialize_FITS_name_entries();

	initialize_log_tags();

	/* register some callbacks */
	register_selection_update_callback(update_export_crop_label);
	register_selection_update_callback(update_display_selection);
	register_selection_update_callback(update_display_fwhm);

	/* support for converting files by dragging onto the GtkTreeView */
	gtk_drag_dest_set(lookup_widget("treeview_convert"),
			GTK_DEST_DEFAULT_MOTION, drop_types, G_N_ELEMENTS(drop_types),
			GDK_ACTION_COPY);

	siril_drag_single_image_set_dest();

	set_GUI_CWD();
	siril_log_message(_("Default FITS extension is set to %s\n"), com.pref.ext);

	init_peaker_GUI();
#ifdef HAVE_LIBRAW
	set_libraw_settings_menu_available(TRUE);	// enable libraw settings
#else
	set_libraw_settings_menu_available(FALSE);	// disable libraw settings
#endif
	update_spinCPU(com.max_thread);

	if (com.pref.first_start) {
		com.pref.first_start = FALSE;
		writeinitfile();

		gchar *ver = g_strdup_printf(_("Welcome to %s"), PACKAGE_STRING);

		int ret = siril_confirm_dialog(ver,
				_("Hello, this is the first time you use this new version of Siril. Please, have a seat and take the time "
						"to watch the short introduction we have prepared for you. "
						"Be aware you can replay this introduction at any times in the Miscellaneous tab of the preferences dialog box."),
						_("See Introduction"));
		if (ret)
			start_intro_script();

		g_free(ver);
	}

	/* every 0.5sec update memory display */
	g_timeout_add(500, update_displayed_memory, NULL);

	/* now that everything is loaded we can connect these signals
	 * Doing it in the glade file is a bad idea because they are called too many times during loading */
	g_signal_connect(lookup_widget("control_window"), "configure-event", G_CALLBACK(on_control_window_configure_event), NULL);
	g_signal_connect(lookup_widget("control_window"), "window-state-event", G_CALLBACK(on_control_window_window_state_event), NULL);
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
	GtkEntry *minentry;
	gchar *buffer;

	minentry = (GtkEntry *)user_data;

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
	GtkEntry *maxentry;
	gchar *buffer;

	maxentry = (GtkEntry *)user_data;

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

void on_focal_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar* focal_entry = gtk_entry_get_text(GTK_ENTRY(editable));
	gfit.focal_length = g_ascii_strtod(focal_entry, NULL);
	drawPlot();
}

void on_pitchX_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar* pitchX_entry = gtk_entry_get_text(GTK_ENTRY(editable));
	gfit.pixel_size_x = (float) g_ascii_strtod(pitchX_entry, NULL);
	drawPlot();
}

void on_pitchY_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar* pitchY_entry = gtk_entry_get_text(GTK_ENTRY(editable));
	gfit.pixel_size_y = (float) g_ascii_strtod(pitchY_entry, NULL);
	drawPlot();
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
	drawPlot();
}

void on_file_information_close_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("file_information"));
}


void on_toggleButtonUnbinned_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkWidget *box;

	box = (GtkWidget *)user_data;

	gfit.unbinned = gtk_toggle_button_get_active(button);
	gtk_widget_set_sensitive(box, gfit.unbinned);
}

void on_button_clear_sample_clicked(GtkButton *button, gpointer user_data) {
	clear_sampling_setting_box();
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
	static GtkWidget *main_w = NULL;
	if (!com.script && com.pref.remember_windows) {
		if (!main_w)
			main_w = lookup_widget("control_window");
		com.pref.main_w_pos = get_window_position(GTK_WINDOW(GTK_APPLICATION_WINDOW(main_w)));
		com.pref.is_maximized = gtk_window_is_maximized(GTK_WINDOW(GTK_APPLICATION_WINDOW(main_w)));
	}
}

void load_main_window_state() {
	GtkWidget *win = lookup_widget("control_window");
	GdkRectangle workarea = { 0 };

	gdk_monitor_get_workarea(gdk_display_get_primary_monitor(gdk_display_get_default()), &workarea);

	int w = com.pref.main_w_pos.w;
	int h = com.pref.main_w_pos.h;

	int x = CLAMP(com.pref.main_w_pos.x, 0, workarea.width - w);
	int y = CLAMP(com.pref.main_w_pos.y, 0, workarea.height - h);

	if (com.pref.remember_windows && w > 0 && h > 0) {
		if (com.pref.is_maximized) {
			gtk_window_maximize(GTK_WINDOW(GTK_APPLICATION_WINDOW(win)));
		} else {
			gtk_window_move(GTK_WINDOW(GTK_APPLICATION_WINDOW(win)), x, y);
			gtk_window_resize(GTK_WINDOW(GTK_APPLICATION_WINDOW(win)), w, h);
		}
	}
}

void gtk_main_quit() {
	writeinitfile();		// save settings (like window positions)
	close_sequence(FALSE);	// save unfinished business
	close_single_image();	// close the previous image and free resources
	g_slist_free_full(com.pref.script_path, g_free);
	exit(EXIT_SUCCESS);
}

void siril_quit() {
	if (com.pref.save.quit) {
		gtk_main_quit();
	}
	gboolean quit = siril_confirm_dialog_and_remember(_("Closing application"),
			_("Are you sure you want to quit?"), _("Exit"), &com.pref.save.quit);
	if (quit) {
		set_GUI_misc();
		writeinitfile();
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

void on_checkchain_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	if (copy_rendering_settings_when_chained(FALSE))
		redraw(com.cvport, REMAP_ALL);
}

void on_max_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar *txt = gtk_entry_get_text(GTK_ENTRY(editable));
	if (g_ascii_isalnum(txt[0])) {

		guint value = g_ascii_strtoull(txt, NULL, 10);

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
	gboolean isalnum = TRUE;

	const gchar *txt = gtk_entry_get_text(GTK_ENTRY(widget));
	int len = gtk_entry_get_text_length(GTK_ENTRY(widget));
	for (int i = 0; i < len; i++)
		if (!g_ascii_isalnum(txt[i])) {
			isalnum = FALSE;
			break;
		}
	if (isalnum == FALSE || len == 0)
		gtk_entry_set_text(GTK_ENTRY(widget), "65535");
	return FALSE;
}

gboolean on_min_entry_focus_out_event(GtkWidget *widget, gpointer user_data) {
	gboolean isalnum = TRUE;

	const gchar *txt = gtk_entry_get_text(GTK_ENTRY(widget));
	int len = gtk_entry_get_text_length(GTK_ENTRY(widget));
	for (int i = 0; i < len; i++)
		if (!g_ascii_isalnum(txt[i])) {
			isalnum = FALSE;
			break;
		}
	if (isalnum == FALSE || len == 0)
		gtk_entry_set_text(GTK_ENTRY(widget), "0");
	return FALSE;
}

void on_min_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar *txt = gtk_entry_get_text(GTK_ENTRY(editable));
	if (g_ascii_isalnum(txt[0])) {

		guint value = g_ascii_strtoull(txt, NULL, 10);

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

void on_regTranslationOnly_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkWidget *Algo = lookup_widget("ComboBoxRegInter");
	GtkWidget *Prefix = lookup_widget("regseqname_entry");

	gtk_widget_set_sensitive(Algo, !gtk_toggle_button_get_active(togglebutton));
	gtk_widget_set_sensitive(Prefix, !gtk_toggle_button_get_active(togglebutton));
}

void on_seqproc_entry_changed(GtkComboBox *widget, gpointer user_data) {
	gchar *name = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(widget));
	if (name && name[0] != '\0') {
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
		gchar *msg = g_strdup_printf(_("Selected %s sequence %s..."), type, name);
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
	update_display_selection();	// update the dimensions of the selection when switching page
	update_display_fwhm();
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
	GtkToggleButton *forceButton = (GtkToggleButton *)user_data;
	int force = gtk_toggle_button_get_active(forceButton);

	if (get_thread_run()) {
		PRINT_ANOTHER_THREAD_RUNNING;
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

void on_button_data_ok_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("data_dialog"));
}

void on_comboboxreglayer_changed(GtkComboBox *widget, gpointer user_data) {
	if (gtk_combo_box_get_active(widget) == -1)
		return;
	free_reference_image();
	update_stack_interface(TRUE);
}

void scrollbars_hadjustment_changed_handler(GtkAdjustment *adjustment,
		gpointer user_data) {
	double value = gtk_adjustment_get_value(adjustment);

	for (int i = 0; i < MAXVPORT; i++) {
		if (com.hadj[i] != adjustment) {
			gtk_adjustment_set_value(com.hadj[i], value);
		}
	}
}

void scrollbars_vadjustment_changed_handler(GtkAdjustment *adjustment,
		gpointer user_data) {
	double value = gtk_adjustment_get_value(adjustment);
	for (int i = 0; i < MAXVPORT; i++) {
		if (com.vadj[i] != adjustment) {
			gtk_adjustment_set_value(com.vadj[i], value);
		}
	}
}

void on_spinCPU_value_changed (GtkSpinButton *spinbutton, gpointer user_data) {
	com.max_thread = gtk_spin_button_get_value_as_int(spinbutton);
}

void on_menu_rgb_align_select(GtkMenuItem *menuitem, gpointer user_data) {
	gboolean sel_is_drawn = ((com.selection.w > 0.0) && (com.selection.h > 0.0));

	gtk_widget_set_sensitive(lookup_widget("rgb_align_dft"), sel_is_drawn);
	gtk_widget_set_sensitive(lookup_widget("rgb_align_psf"), sel_is_drawn);
}

void on_rgb_align_dft_activate(GtkMenuItem *menuitem, gpointer user_data) {
	undo_save_state(&gfit, _("RGB alignment (DFT)"));
	rgb_align(1);
}

void on_rgb_align_psf_activate(GtkMenuItem *menuitem, gpointer user_data) {
	undo_save_state(&gfit, _("RGB alignment (PSF)"));
	rgb_align(0);
}

void on_gotoStacking_button_clicked(GtkButton *button, gpointer user_data) {
	control_window_switch_to_tab(STACKING);
}

void on_button_paned_clicked(GtkButton *button, gpointer user_data) {
	static gboolean is_extended = TRUE;
	GtkPaned *paned = (GtkPaned*) user_data;
	GtkImage *image = GTK_IMAGE(gtk_bin_get_child(GTK_BIN(button)));
	GtkWidget *widget = gtk_paned_get_child2(paned);

	gtk_widget_set_visible(widget, !is_extended);

	if (!is_extended) {
		gtk_image_set_from_icon_name(image, "pan-end-symbolic",
				GTK_ICON_SIZE_BUTTON);
	} else {
		gtk_image_set_from_icon_name(image, "pan-start-symbolic",
				GTK_ICON_SIZE_BUTTON);
	}
	is_extended = !is_extended;
}
