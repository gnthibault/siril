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
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <string.h>
#include <assert.h>
#include <stdlib.h>

#include "compositing/compositing.h"
#include "core/siril.h"
#include "core/proto.h"
#include "algos/colors.h"
#include "io/single_image.h"
#include "gui/PSF_list.h"
#include "gui/callbacks.h"
#include "registration/registration.h"
#include "compositing/filters.h"
#include "stacking/stacking.h"
#ifdef HAVE_OPENCV
#include "opencv/opencv.h"
#endif

#undef DEBUG

static int compositing_loaded = 0;

typedef enum {
	HSL,
	HSV,
	CIELAB
} coloring_type_enum;

static coloring_type_enum coloring_type = HSL;

/* The result is stored in gfit.
 * gfit.rx and gfit.ry are the reference 1x1 binning output image size. */

#define MAX_LAYERS 8

/* the structure storing information for each layer to be composed
 * (one layer = one source image) and one associated colour */
typedef struct {
	/* widgets data */
	GtkButton *remove_button;
	GtkDrawingArea *color_w;	// the simulated color chooser
	GtkFileChooserButton *chooser;	// the file choosers
	GtkLabel *label;		// the labels
	GtkSpinButton *spinbutton_x;	// the X spin button
	GtkSpinButton *spinbutton_y;	// the Y spin button

	/* useful data */
	GdkRGBA color;			// real color of the layer
	GdkRGBA saturated_color;	// saturated color of the layer
	fits the_fit;			// the fits for layers
} layer;

/* the list of layers. It is dynamic in content but fixed in size.
 * The first element is reserved for luminance and cannot be removed. */
static layer *layers[MAX_LAYERS+1];	// NULL terminated
static int layers_count = 1;		// the glade has only luminance

static int luminance_mode = 0;		// 0 if luminance is not used

static struct registration_method *reg_methods[3];

static sequence *seq = NULL;		// the sequence of layers, for alignments and normalization
static norm_coeff *coeff = NULL;	// the normalization coefficients

/* special case of the color associated to luminance */
void set_luminance(GdkRGBA *rgba) { rgba->red = -42.0; }
int is_luminance(GdkRGBA *rgba) { return (rgba->red == -42.0); }

static GdkRGBA list_of_12_palette_colors[12];
static const char *list_of_12_color_names[12] = {
	"#ff0000", "#7f0000",
	"#00ff00", "#007f00",
	"#0000ff", "#00007f",
	"#ffff00", "#7f7f00",
	"#ff00ff", "#7f007f",
	"#00ffff", "#007f7f"
};

static GtkColorChooserDialog *color_dialog = NULL;
static int current_layer_color_choosing = 0;
static int color_quick_edit = 0;
static GdkRGBA qe_ref_color;
static GtkEntry *wl_entry = NULL;
static GtkComboBoxText *box = NULL;

static GtkGrid *grid_layers = NULL;

static GtkButton *add_button = NULL;

/******* internal functions *******/
void add_the_layer_add_button();
void grid_add_row(int layer, int index, int first_time);
void grid_remove_row(int layer, int free_the_row);
int has_fit(int layer);
void update_compositing_interface();
WORD get_composition_pixel_value(int fits_index, int reg_layer, int x, int y);
void increment_pixel_components_from_layer_value(int fits_index, GdkRGBA *rgbpixel, WORD layer_pixel_value);
void increment_pixel_components_from_layer_saturated_value(int fits_index, GdkRGBA *rgbpixel, WORD layer_pixel_value);
void colors_align_and_compose();		// the rgb procedure
void luminance_and_colors_align_and_compose();	// the lrgb procedure
void color_has_been_updated(int layer);
void update_color_from_saturation(int layer, double newl);
void rgb_pixel_limiter(GdkRGBA *pixel);
void clear_pixel(GdkRGBA *pixel);
void update_result(int and_refresh);
void populate_filter_lists();
void coeff_clear();

/* callbacks for programatic GTK */
void on_layer_remove(GtkButton *button, gpointer user_data);
gboolean on_color_button_press_event(GtkDrawingArea *widget, GdkEventButton *event, gpointer user_data);
gboolean on_color_button_release_event(GtkDrawingArea *widget, GdkEventButton *event, gpointer user_data);
gboolean on_color_button_motion_event(GtkWidget *widget, GdkEventMotion *event, gpointer user_data);
gboolean draw_layer_color(GtkDrawingArea *widget, cairo_t *cr, gpointer data);
void on_filechooser_file_set(GtkFileChooserButton *widget, gpointer user_data);

/********************************************************/

/* the compositing menu callback */
void on_menu_compositing_activate(GtkMenuItem *menuitem, gpointer user_data) {
	open_compositing_window();
}

/* creates a new row with all widgets and bindings at the row index in the
 * layers grid. Indices start at 0, but row 0 holds only one label, and row 1 is
 * reserved to the luminance layer. */
layer *create_layer(int index) {
	layer *ret = malloc(sizeof(layer));

	/* create the widgets and set properties and signals */
	//ret->remove_button = GTK_BUTTON(gtk_button_new_from_icon_name("list-remove", GTK_ICON_SIZE_BUTTON));
	ret->remove_button = GTK_BUTTON(gtk_button_new());
	gtk_button_set_image(ret->remove_button,
			gtk_image_new_from_icon_name("list-remove", GTK_ICON_SIZE_BUTTON));
	g_object_ref(G_OBJECT(ret->remove_button));	// don't destroy it on removal from grid
	g_signal_connect(ret->remove_button, "clicked", G_CALLBACK(on_layer_remove), NULL);

	ret->color_w = GTK_DRAWING_AREA(gtk_drawing_area_new());
	gtk_widget_set_events(GTK_WIDGET(ret->color_w),
			GDK_BUTTON_MOTION_MASK | GDK_BUTTON_PRESS_MASK |
			GDK_BUTTON_RELEASE_MASK);
	g_signal_connect(GTK_WIDGET(ret->color_w), "button-release-event", G_CALLBACK(on_color_button_release_event), NULL);
	g_signal_connect(GTK_WIDGET(ret->color_w), "button-press-event", G_CALLBACK(on_color_button_press_event), NULL);
	g_signal_connect(GTK_WIDGET(ret->color_w), "motion-notify-event", G_CALLBACK(on_color_button_motion_event), NULL);
	g_signal_connect(ret->color_w, "draw", G_CALLBACK(draw_layer_color), NULL);
	g_object_ref(G_OBJECT(ret->color_w));	// don't destroy it on removal from grid

	ret->chooser = GTK_FILE_CHOOSER_BUTTON(
			gtk_file_chooser_button_new(
				"Select source image", GTK_FILE_CHOOSER_ACTION_OPEN));
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(ret->chooser), com.wd);
	gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(ret->chooser),
			GTK_FILE_FILTER(gtk_builder_get_object(builder, "filefilter1")));
			gtk_file_chooser_button_set_width_chars(ret->chooser, 16);
	g_signal_connect(ret->chooser, "file-set", G_CALLBACK(on_filechooser_file_set), NULL);
	g_object_ref(G_OBJECT(ret->chooser));	// don't destroy it on removal from grid

	ret->label = GTK_LABEL(gtk_label_new(_("not loaded")));
	g_object_ref(G_OBJECT(ret->label));	// don't destroy it on removal from grid

	ret->spinbutton_x = GTK_SPIN_BUTTON(gtk_spin_button_new_with_range(-1000.0, 1000.0, 1.0));
	gtk_spin_button_set_value(ret->spinbutton_x, 0.0);
	gtk_widget_set_sensitive(GTK_WIDGET(ret->spinbutton_x), FALSE);
	g_object_ref(G_OBJECT(ret->spinbutton_x));	// don't destroy it on removal from grid

	ret->spinbutton_y = GTK_SPIN_BUTTON(gtk_spin_button_new_with_range(-1000.0, 1000.0, 1.0));
	gtk_spin_button_set_value(ret->spinbutton_y, 0.0);
	gtk_widget_set_sensitive(GTK_WIDGET(ret->spinbutton_y), FALSE);
	g_object_ref(G_OBJECT(ret->spinbutton_y));	// don't destroy it on removal from grid

	/* set other layer data */
	memset(&ret->the_fit, 0, sizeof(fits));
	assert(index >= 2);	// 1 is luminance
	if (index <= 7)		// copy default RGB colours
		memcpy(&ret->color,
				&list_of_12_palette_colors[(index-2)*2],
				sizeof(GdkRGBA));
	else clear_pixel(&ret->color);
	clear_pixel(&ret->saturated_color);
	return ret;
}

/* callback of the '+' button that is clicked to add a layer in the list */
void on_layer_add(GtkButton *button, gpointer user_data) {
	++layers_count;

	/* move down the plus button */
	gtk_container_remove(GTK_CONTAINER(grid_layers), GTK_WIDGET(add_button));
	if (layers_count < MAX_LAYERS)
		add_the_layer_add_button();

	/* add the new layer */
	layers[layers_count-1] = create_layer(layers_count);
	layers[layers_count] = NULL;

	grid_add_row(layers_count-1, layers_count, 1);

	color_has_been_updated(layers_count-1);

	coeff_clear();
}

/* adds the '+' button at the bottom of the list. Creates the trailing grid row too */
void add_the_layer_add_button() {
	int first_time = 0;
	if (!add_button) {
		first_time = 1;
		add_button = GTK_BUTTON(gtk_button_new());
		gtk_button_set_image(add_button,
				gtk_image_new_from_icon_name("list-add", GTK_ICON_SIZE_BUTTON));
		g_object_ref(G_OBJECT(add_button));	// don't destroy it on removal from grid
		g_signal_connect(add_button, "clicked", G_CALLBACK(on_layer_add), NULL);
	}

	gtk_grid_attach(grid_layers, GTK_WIDGET(add_button), 0, layers_count+1, 1, 1);
	if (first_time)
		gtk_widget_show(GTK_WIDGET(add_button));
}

/* callback of the '-' button that is clicked to remove a layer in the list */
void on_layer_remove(GtkButton *button, gpointer user_data) {
	int layer, refresh = 0;
	for (layer = 1; layers[layer]; layer++)
		if (layers[layer]->remove_button == button)
			break;
	if (!layers[layer]) return;

	if (layer != MAX_LAYERS-1) {
		// the add button is not present if we're at the maximum number of layers
		gtk_container_remove(GTK_CONTAINER(grid_layers), GTK_WIDGET(add_button));
	}

	if (has_fit(layer)) {
		clearfits(&layers[layer]->the_fit);
		refresh = 1;
	}
	grid_remove_row(layer, 1);	// remove these widgets for good
	free(layers[layer]);

	do {
		layers[layer] = layers[layer+1];
		grid_remove_row(layer, 0);		// switch rows
		grid_add_row(layer, layer+1, 0);
		layer++;
	} while (layers[layer]) ;

	--layers_count;

	add_the_layer_add_button();

	coeff_clear();

	if (refresh)
		update_result(1);
}

void grid_remove_row(int layer, int free_the_row) {
	GtkContainer *cont = GTK_CONTAINER(grid_layers);
	if (!layers[layer]) return;
	gtk_container_remove(cont, GTK_WIDGET(layers[layer]->remove_button));
	gtk_container_remove(cont, GTK_WIDGET(layers[layer]->color_w));
	gtk_container_remove(cont, GTK_WIDGET(layers[layer]->chooser));
	gtk_container_remove(cont, GTK_WIDGET(layers[layer]->label));
	gtk_container_remove(cont, GTK_WIDGET(layers[layer]->spinbutton_x));
	gtk_container_remove(cont, GTK_WIDGET(layers[layer]->spinbutton_y));
	if (free_the_row) {
		g_object_unref(G_OBJECT(layers[layer]->remove_button));	// don't destroy it on removal from grid
		g_object_unref(G_OBJECT(layers[layer]->color_w));	// don't destroy it on removal from grid
		g_object_unref(G_OBJECT(layers[layer]->chooser));	// don't destroy it on removal from grid
		g_object_unref(G_OBJECT(layers[layer]->label));	// don't destroy it on removal from grid
		g_object_unref(G_OBJECT(layers[layer]->spinbutton_x));	// don't destroy it on removal from grid
		g_object_unref(G_OBJECT(layers[layer]->spinbutton_y));	// don't destroy it on removal from grid
	}
}

void grid_add_row(int layer, int index, int first_time) {
	if (!layers[layer]) return;
	gtk_grid_attach(grid_layers, GTK_WIDGET(layers[layer]->remove_button),	0, index, 1, 1);
	gtk_grid_attach(grid_layers, GTK_WIDGET(layers[layer]->color_w),	1, index, 1, 1);
	gtk_grid_attach(grid_layers, GTK_WIDGET(layers[layer]->chooser),	2, index, 1, 1);
	gtk_grid_attach(grid_layers, GTK_WIDGET(layers[layer]->label),		3, index, 1, 1);
	gtk_grid_attach(grid_layers, GTK_WIDGET(layers[layer]->spinbutton_x),	4, index, 1, 1);
	gtk_grid_attach(grid_layers, GTK_WIDGET(layers[layer]->spinbutton_y),	5, index, 1, 1);

	if (first_time) {
		gtk_widget_show(GTK_WIDGET(layers[layer]->remove_button));
		gtk_widget_show(GTK_WIDGET(layers[layer]->color_w));
		gtk_widget_show(GTK_WIDGET(layers[layer]->chooser));
		gtk_widget_show(GTK_WIDGET(layers[layer]->label));
		gtk_widget_show(GTK_WIDGET(layers[layer]->spinbutton_x));
		gtk_widget_show(GTK_WIDGET(layers[layer]->spinbutton_y));
	}
}

/* load all glade data, connect signals, configure the dynamic objects of the
 * composition window and make it visible */
void open_compositing_window() {
	int i;
	if (!compositing_loaded) {
		register_selection_update_callback(update_compositing_interface);

		gtk_builder_connect_signals (builder, NULL);

		/* parse the default palette */
		for (i=0; i<sizeof(list_of_12_color_names)/sizeof(const char*); i++)
			gdk_rgba_parse(&list_of_12_palette_colors[i], list_of_12_color_names[i]);
		color_dialog = GTK_COLOR_CHOOSER_DIALOG(gtk_builder_get_object(builder, "colorchooserdialog"));
		gtk_color_chooser_add_palette(GTK_COLOR_CHOOSER(color_dialog),
				GTK_ORIENTATION_VERTICAL, 2, 12, list_of_12_palette_colors);
		populate_filter_lists();
		wl_entry = GTK_ENTRY(gtk_builder_get_object(builder, "entry_wavelength"));
		box = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(builder, "comboboxtext_filters"));


		/* allocate default layers and populate widget data */
		grid_layers = GTK_GRID(gtk_builder_get_object(builder, "grid_layers"));
		add_the_layer_add_button();

		layers[0] = calloc(1, sizeof(layer));
		layers[0]->chooser = GTK_FILE_CHOOSER_BUTTON(gtk_builder_get_object(builder, "filechooser_lum"));
		gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(layers[0]->chooser), com.wd);
		layers[0]->label = GTK_LABEL(gtk_builder_get_object(builder, "label_lum"));
		layers[0]->spinbutton_x = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spinbutton_lum_x"));
		layers[0]->spinbutton_y = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spinbutton_lum_y"));
		set_luminance(&layers[0]->color);	// special case of luminance

		for (i=1; i<4; i++) {
			/* Create the three default layers */
			on_layer_add(NULL, NULL);
		}
		layers[i] = NULL;

		/* the list below depends on the content of the glade file. It
		 * should be done in the same way as in registration.c, but it
		 * woud be easier if the two glades are merged. */
		reg_methods[0] = new_reg_method(_("One star registration (deep-sky)"), &register_shift_fwhm, REQUIRES_ANY_SELECTION, REGTYPE_DEEPSKY);
		reg_methods[1] = new_reg_method(_("Image pattern alignment (planetary/deep-sky)"), &register_shift_dft, REQUIRES_SQUARED_SELECTION, REGTYPE_PLANETARY);

		reg_methods[2] = NULL;
		update_compositing_interface();
		/* fill compositing_align_method_combo */
		GtkComboBoxText *aligncombo = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(builder, "compositing_align_method_combo"));
		gtk_combo_box_text_remove_all(aligncombo);
		i = 0;
		while (reg_methods[i] != NULL) {
			gtk_combo_box_text_append_text(aligncombo, reg_methods[i]->name);
			i++;
		}
		if (i > 0) {
			gtk_combo_box_set_active(GTK_COMBO_BOX(aligncombo), com.reg_settings);
		}

		compositing_loaded = 1;
	} else {
		/* not the first load, update the CWD just in case it changed in the meantime */
		i = 0;
		do {
			gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(layers[i]->chooser), com.wd);
			i++;
		} while (layers[i]) ;
	}

	if (compositing_loaded == 1)
		gtk_widget_show(lookup_widget("composition_dialog"));
}

/* returns true if the layer number layer has a loaded FITS image */
int has_fit(int layer) {
	return (layers[layer] && layers[layer]->the_fit.rx != 0);
}

/* returns true if only one image is loaded */
int number_of_images_loaded() {
	int i, count = 0;
	for (i=0; layers[i]; i++)
		if (has_fit(i))
			count++;
	return count;
}

/* returns true if all layers have an image loaded */
int all_images_loaded() {
	int i;
	for (i=0; layers[i]; i++) {
		if (i == 0 && !luminance_mode) continue;
		if (has_fit(i))
			return 0;
	}
	return 1;
}

/* returns true if none of the colour layers have an image loaded */
int no_color_available() {	// don't test luminance
	int i;
	for (i=1; layers[i]; i++) {
		if (has_fit(i))
			return 0;
	}
	return 1;
}

/* the 'enable luminance' checkbox callback */
void on_composition_use_lum_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	luminance_mode = gtk_toggle_button_get_active(togglebutton);
	if (has_fit(0) && number_of_images_loaded() >= 1)
		update_result(1);
}

/* callback for the file chooser's file selection: try to load the pointed file, allocate the
 * destination image if this is the first image, and update the result. */
void on_filechooser_file_set(GtkFileChooserButton *widget, gpointer user_data) {
	int layer, retval;
	char buf[48], *filename;

	for (layer = 0; layers[layer]; layer++)
		if (layers[layer]->chooser == widget)
			break;
	if (!layers[layer]) return;	// not found
	/* layer is the number of the opened layer */

	filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(widget));
	if (!filename) return;
	if ((retval = read_single_image(filename, &layers[layer]->the_fit, NULL))) {
		gtk_label_set_text(layers[layer]->label, _("ERROR"));
	} else {
		/* Force first tab to be Red and not B&W if an image was already loaded */
		GtkNotebook* Color_Layers = GTK_NOTEBOOK(gtk_builder_get_object(builder, "notebook1"));
		GtkWidget *page = gtk_notebook_get_nth_page(Color_Layers, RED_VPORT);
		gtk_notebook_set_tab_label_text(GTK_NOTEBOOK(Color_Layers), page, _("Red channel"));

		if (number_of_images_loaded() > 1 &&
				(gfit.rx != layers[layer]->the_fit.rx ||
				 gfit.ry != layers[layer]->the_fit.ry)) {
			/* TODO: handle the binning cases. Values should be stored
			 * in, or even taken from fit->binning_x and binning_y */
			//if (layers[layer]->the_fit.binning_x > 1 ||
			//layers[layer]->the_fit.binning_y > 1)
#ifdef HAVE_OPENCV
			if (gfit.rx < layers[layer]->the_fit.rx ||
					gfit.ry < layers[layer]->the_fit.ry) {
				siril_log_message(_("The first loaded image should have the greatest sizes for now\n"));
				sprintf(buf, _("NOT OK %ux%u"), layers[layer]->the_fit.rx, layers[layer]->the_fit.ry);
				gtk_label_set_text(layers[layer]->label, buf);
				retval = 1;
			} else {
				siril_log_message(_("Resizing the loaded image from %dx%d to %dx%d\n"),
						layers[layer]->the_fit.rx,
						layers[layer]->the_fit.ry, gfit.rx, gfit.ry);
				sprintf(buf, _("OK upscaled from %ux%u"),
						layers[layer]->the_fit.rx, layers[layer]->the_fit.ry);
				cvResizeGaussian(&layers[layer]->the_fit, gfit.rx, gfit.ry, OPENCV_LINEAR); // BILINEAR
				image_find_minmax(&layers[layer]->the_fit, 0);
				gtk_label_set_text(layers[layer]->label, buf);
			}
#else
			siril_log_message(_("You need to install opencv to compose images with different sizes\n"));
			sprintf(buf, _("NOT OK %ux%u"), layers[layer]->the_fit.rx, layers[layer]->the_fit.ry);
			gtk_label_set_text(layers[layer]->label, buf);
			retval = 1;
#endif
		}
		else if (!retval) {
			image_find_minmax(&layers[layer]->the_fit, 0);
			sprintf(buf, _("OK %ux%u"), layers[layer]->the_fit.rx, layers[layer]->the_fit.ry);
			gtk_label_set_text(layers[layer]->label, buf);
		}
	}
	g_free(filename);

	/* special case of luminance selected */
	if (layer == 0) {
		GtkToggleButton *lum_button = GTK_TOGGLE_BUTTON(gtk_builder_get_object(builder, "composition_use_lum"));
		g_signal_handlers_block_by_func(lum_button, on_composition_use_lum_toggled, NULL);
		gtk_toggle_button_set_active(lum_button, !retval);
		g_signal_handlers_unblock_by_func(lum_button, on_composition_use_lum_toggled, NULL);
		luminance_mode = !retval;
	}
	if (retval) {
		clearfits(&layers[layer]->the_fit);
		return;
	}

	/* create the new result image if it's the first opened image */
	if (number_of_images_loaded() == 1) {
		close_single_image();
		copyfits(&layers[layer]->the_fit, &gfit, CP_ALLOC | CP_INIT | CP_FORMAT | CP_EXPAND, -1);
	}

	update_result(0);

	/* open the single image.
	 * code taken from stacking.c:start_stacking() and read_single_image() */
	if (number_of_images_loaded() == 1) {
		clear_stars_list();
		com.seq.current = UNRELATED_IMAGE;
		com.uniq = calloc(1, sizeof(single));
		com.uniq->comment = strdup("Compositing result image");
		com.uniq->filename = strdup(_("Unsaved compositing result"));
		com.uniq->nb_layers = gfit.naxes[2];
		com.uniq->layers = calloc(com.uniq->nb_layers, sizeof(layer_info));
		com.uniq->fit = &gfit;
		display_filename();
		sliders_mode_set_state(com.sliders);

		image_find_minmax(&gfit, 0);
		init_layers_hi_and_lo_values(MIPSLOHI);
		set_cutoff_sliders_max_values();
		set_cutoff_sliders_values();
		set_display_mode();
		redraw(com.cvport, REMAP_ALL);
		update_used_memory();
		show_main_gray_window();
		show_rgb_window();
		sequence_list_change_current();
	}
	else {
		update_MenuItem();
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
	}
}

void create_the_internal_sequence() {
	int i, j, nb_layers;
	if (seq) free_sequence(seq, TRUE);

	nb_layers = number_of_images_loaded();
	if (nb_layers == 0 || nb_layers == 1) {
		char *msg = siril_log_message(_("You must at least load two layers before!\n"));
		show_dialog(msg, _("Warning"), "gtk-dialog-warning");
		seq = NULL;
		return;
	}
	if (luminance_mode) {
		seq = create_internal_sequence(nb_layers);
		for (j = 0, i=0; i<layers_count; i++) {
			if (has_fit(i)) {
				internal_sequence_set(seq, j, &layers[i]->the_fit);
				j++;
			}
		}
	} else {
		if (has_fit(0))
			nb_layers--;
		seq = create_internal_sequence(nb_layers);
		for (j=0, i=1; i<layers_count; i++) {
			if (has_fit(i)) {
				internal_sequence_set(seq, j, &layers[i]->the_fit);
				j++;
			}
		}
	}

	seq->rx = gfit.rx;
	seq->ry = gfit.ry;
}

/* start alignming the layers: create an 'internal' sequence and run the selected method on it */
void on_button_align_clicked(GtkButton *button, gpointer user_data) {
	int i, j;
	struct registration_args regargs;
	struct registration_method *method;
	char *msg;
	GtkComboBox *regcombo;

	create_the_internal_sequence();

	/* align it */
	regcombo = GTK_COMBO_BOX(gtk_builder_get_object(builder, "compositing_align_method_combo"));
	method = reg_methods[gtk_combo_box_get_active(regcombo)];

	regargs.seq = seq;
	regargs.process_all_frames = TRUE;
	get_the_registration_area(&regargs, method);
	regargs.layer = 0;	// TODO: fix with dynamic layers list
	regargs.run_in_thread = FALSE;

	msg = siril_log_message(_("Starting registration using method: %s\n"), method->name);
	msg[strlen(msg)-1] = '\0';
	set_cursor_waiting(TRUE);
	set_progress_bar_data(msg, PROGRESS_RESET);
	if (method->method_ptr(&regargs))
		set_progress_bar_data(_("Error in layers alignment."), PROGRESS_DONE);
	else set_progress_bar_data(_("Registration complete."), PROGRESS_DONE);
	set_cursor_waiting(FALSE);

	/* display the values */
	if (luminance_mode) {
		for (j=0, i=0; i<layers_count; i++) {
			if (has_fit(i)) {
				gtk_spin_button_set_value(layers[i]->spinbutton_x, seq->regparam[0][j].shiftx);
				gtk_spin_button_set_value(layers[i]->spinbutton_y, seq->regparam[0][j].shifty);
				j++;
			}
		}
	} else {
		for (j=0, i=1; i<layers_count; i++) {
			if (has_fit(i)) {
				gtk_spin_button_set_value(layers[i]->spinbutton_x, seq->regparam[0][j].shiftx);
				gtk_spin_button_set_value(layers[i]->spinbutton_y, seq->regparam[0][j].shifty);
				j++;
			}
		}
	}

	/* align the image and display it.
	 * Layers are aligned against the reference layer, with zeros where there is not data */
	update_result(1);
}

WORD get_normalized_pixel_value(int fits_index, WORD layer_pixel_value) {
	double tmp = (double)layer_pixel_value;
	if (!has_fit(0))
		fits_index--;
	tmp *= coeff->scale[fits_index];
	tmp -= coeff->offset[fits_index];
	return round_to_WORD(tmp);
}

/* get the pixel value at coordinates x,y for image in layers[fits_index]->the_fit.
 * x and y are given in buffer coordinates, not image coordinates.
 * Handles (not yet - binning and) registration offset */
WORD get_composition_pixel_value(int fits_index, int reg_layer, int x, int y) {
	int realX = x, realY = y;
	if (seq && seq->regparam && reg_layer < seq->number && reg_layer >= 0) {
		// all images have one layer, hence the 0 below
		realX = x - seq->regparam[0][reg_layer].shiftx;
		if (realX < 0 || realX >= gfit.rx) return (WORD)0;
		realY = y - seq->regparam[0][reg_layer].shifty;
		if (realY < 0 || realY >= gfit.ry) return (WORD)0;
	}
	WORD pixel_value = layers[fits_index]->the_fit.pdata[0][realX + realY * gfit.rx];
	if (coeff) {
		// normalization
		pixel_value = get_normalized_pixel_value(fits_index, pixel_value);
	}
	return pixel_value;
}

/* increments the color values in rgbpixel from the pixel value for a particular
 * layer. GdkRGBA values are stored in the [0, 1] interval. */
void increment_pixel_components_from_layer_value(int fits_index, GdkRGBA *rgbpixel, WORD layer_pixel_value) {
	GdkRGBA *layer_color = &layers[fits_index]->color;
	rgbpixel->red += layer_color->red * layer_pixel_value / USHRT_MAX_DOUBLE;
	rgbpixel->green += layer_color->green * layer_pixel_value / USHRT_MAX_DOUBLE;
	rgbpixel->blue += layer_color->blue * layer_pixel_value / USHRT_MAX_DOUBLE;
}

/* increments the color values in rgbpixel from the saturated pixel value for a
 * particular layer. GdkRGBA values are stored in the [0, 1] interval. */
void increment_pixel_components_from_layer_saturated_value(int fits_index, GdkRGBA *rgbpixel, WORD layer_pixel_value) {
	GdkRGBA *layer_color = &layers[fits_index]->saturated_color;
	rgbpixel->red += layer_color->red * layer_pixel_value / USHRT_MAX_DOUBLE;
	rgbpixel->green += layer_color->green * layer_pixel_value / USHRT_MAX_DOUBLE;
	rgbpixel->blue += layer_color->blue * layer_pixel_value / USHRT_MAX_DOUBLE;
}

/* called when selection changed */
void update_compositing_interface() {
	if (!builder) return;
	GtkLabel *label = GTK_LABEL(gtk_builder_get_object(builder, "label_msg"));
	GtkComboBox *combo = GTK_COMBO_BOX(gtk_builder_get_object(builder, "compositing_align_method_combo"));
	//int ref_layer = gtk_combo_box_get_active(GTK_COMBO_BOX(gtk_builder_get_object(builder, "compositing_align_layer_combo")));
	int sel_method = gtk_combo_box_get_active(combo);
	/* select default method as function of selection size */
	if (sel_method == -1 && com.selection.w > 0 && com.selection.h > 0) {
		if (com.selection.w > 180 || com.selection.h > 180)
			gtk_combo_box_set_active(combo, 0);	// activate DFT
		else gtk_combo_box_set_active(combo, 1);	// activate FWHM
	}

	if (com.selection.w <= 0 && com.selection.h <= 0) {
		gtk_label_set_text(label, _("An image area must be selected for align"));
		gtk_widget_set_sensitive(lookup_widget("button_align"), FALSE);
	/*} else if (ref_layer == -1 || (!luminance_mode && ref_layer == 0)) {
		gtk_label_set_text(label, "A reference layer must be selected for align");
		gtk_widget_set_sensitive(lookup_widget("button_align"), FALSE);*/
	} else if (number_of_images_loaded() < 2) {
		gtk_label_set_text(label, _("At least 2 channels must be loaded for align"));
		gtk_widget_set_sensitive(lookup_widget("button_align"), FALSE);
	} else {
		gtk_label_set_text(label, "");
		gtk_widget_set_sensitive(lookup_widget("button_align"), TRUE);
	}

}

/* callback for changes of the selected reference layer */
void on_compositing_align_layer_combo_changed(GtkComboBox *widget, gpointer user_data) {
	update_compositing_interface();
}


void on_composition_combo_coloringtype_changed(GtkComboBox *widget, gpointer user_data) {
	coloring_type = gtk_combo_box_get_active(widget);
	update_result(1);
}

/* Image composition without luminance. Used for RGB composition for example.
 * Result is in gfit. */
void colors_align_and_compose() {
	int x, y, i = 0;	// i is browsing the 1D buffer, i = y * rx + x
	if (no_color_available()) return;
	fprintf(stdout, "colour layers only composition\n");
	for (y = 0; y < gfit.ry; ++y) {
		for (x = 0; x < gfit.rx; ++x) {
			int layer;
			GdkRGBA pixel;
			clear_pixel(&pixel);
			for (layer = 1; layers[layer]; layer++) {
				if (has_fit(layer)) {
					int reg_layer = seq ? internal_sequence_find_index(seq, &layers[layer]->the_fit) : -1;
					WORD layer_value = get_composition_pixel_value(layer, reg_layer, x, y);
					if (layer_value != (WORD)0)
						increment_pixel_components_from_layer_value(layer, &pixel, layer_value);
				}
			}

			rgb_pixel_limiter(&pixel);
			gfit.pdata[RLAYER][i] = round_to_WORD(pixel.red * USHRT_MAX_DOUBLE);
			gfit.pdata[GLAYER][i] = round_to_WORD(pixel.green * USHRT_MAX_DOUBLE);
			gfit.pdata[BLAYER][i] = round_to_WORD(pixel.blue * USHRT_MAX_DOUBLE);
			i++;
		}
	}
}

/* This function fills the data in the gfit image with LRGB information from
 * layers[*]->the_fit images. Layers are aligned with registration data, no
 * binning yet. */
void luminance_and_colors_align_and_compose() {
	/* Each pixel is transformed from RGB to HSI, I is replaced by the
	 * luminance layer's value and transformed back to RGB. */
	guint x, y;
	assert(has_fit(0));

	if (no_color_available()) {
		/* luminance only: we copy its data to all result layers */
		int i;
		unsigned int nbdata = gfit.rx * gfit.ry;
		fprintf(stdout, "luminance-only, no composition\n");
		for (i=0; i<3; i++)
			memcpy(gfit.pdata[i], layers[0]->the_fit.data, nbdata*sizeof(WORD));
		return;
	}
	fprintf(stdout, "luminance-enabled composition\n");

	double norm = (double)(layers[0]->the_fit.maxi);

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(y,x) schedule(static)
#endif
	for (y = 0; y < gfit.ry; y++) {
		for (x = 0; x < gfit.rx; x++) {
			int layer;
			gdouble h, s, i;
			gdouble X, Y, Z;
			gdouble a, b;
			/* get color information */
			GdkRGBA pixel;
			clear_pixel(&pixel);
			for (layer = 1; layers[layer]; layer++) {
				if (has_fit(layer)) {
					WORD layer_value = get_composition_pixel_value(layer, layer, x, y);
					if (layer_value != (WORD)0)
						increment_pixel_components_from_layer_value(layer, &pixel, layer_value);
				}
			}
			rgb_pixel_limiter(&pixel);

			switch (coloring_type) {
			case HSL:
				rgb_to_hsl(pixel.red, pixel.green, pixel.blue, &h, &s, &i);
				/* add luminance by replacing it in the HSI */
				i = (double) get_composition_pixel_value(0, 0, x, y) / norm;
				/* converting back to RGB */
				hsl_to_rgb(h, s, i, &pixel.red, &pixel.green, &pixel.blue);
				break;
			case HSV:
				rgb_to_hsv(pixel.red, pixel.green, pixel.blue, &h, &s, &i);
				/* add luminance by replacing it in the HSI */
				i = (double) get_composition_pixel_value(0, 0, x, y) / norm;
				/* converting back to RGB */
				hsv_to_rgb(h, s, i, &pixel.red, &pixel.green, &pixel.blue);
				break;
			case CIELAB:
				rgb_to_xyz(pixel.red, pixel.green, pixel.blue, &X, &Y, &Z);
				xyz_to_LAB(X, Y, Z, &i, &a, &b);
				i = (double) get_composition_pixel_value(0, 0, x, y) / norm;
				i *= 100.0;		// 0 < L < 100
				LAB_to_xyz(i, a, b, &X, &Y, &Z);
				xyz_to_rgb(X, Y, Z, &pixel.red, &pixel.green, &pixel.blue);
				break;
			}

			rgb_pixel_limiter(&pixel);

			/* and store in gfit */
			int dst_index = y * gfit.rx + x;
			gfit.pdata[RLAYER][dst_index] = round_to_WORD(pixel.red * USHRT_MAX_DOUBLE);
			gfit.pdata[GLAYER][dst_index] = round_to_WORD(pixel.green * USHRT_MAX_DOUBLE);
			gfit.pdata[BLAYER][dst_index] = round_to_WORD(pixel.blue * USHRT_MAX_DOUBLE);
		}
	}
}

void on_compositing_cancel_clicked(GtkButton *button, gpointer user_data){
	gtk_widget_hide(lookup_widget("composition_dialog"));
}

/* When summing all layers to get the RGB values for one pixel, it may overflow.
 * This procedure defines what happens in that case. */
void rgb_pixel_limiter(GdkRGBA *pixel) {
#ifdef DEBUG
	if (pixel->red > 1.2 || pixel->green > 1.2 || pixel->blue > 1.2)
		fprintf(stdout, "large overflow %g,%g,%g\n", pixel->red,
				pixel->green, pixel->blue);
#endif
	if (pixel->red >= 1.0)
		pixel->red = 1.0;
	if (pixel->green >= 1.0)
		pixel->green = 1.0;
	if (pixel->blue >= 1.0)
		pixel->blue = 1.0;
}

/* initializes a GdkRGBA to black */
void clear_pixel(GdkRGBA *pixel) {
	pixel->red = 0.0;
	pixel->green = 0.0;
	pixel->blue = 0.0;
	pixel->alpha = 1.0;
}

/* recompute the layer composition and optionnally refresh the displayed result image */
void update_result(int and_refresh) {
	if (luminance_mode && has_fit(0)) {
		luminance_and_colors_align_and_compose();
	} else {
		colors_align_and_compose();
	}
	if (and_refresh) {
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
	}
}

/****************** colour management ******************/

// update the saturated color from the new real colour
void color_has_been_updated(int layer) {
	double h, s, v;
	GdkRGBA *real = &layers[layer]->color;
	GdkRGBA *satu = &layers[layer]->saturated_color;
	rgb_to_hsv(real->red, real->green, real->blue, &h,&s,&v);
	printf("%d: saturation: %g, light: %g\n", layer, s, v);
	/* in HSL, the actual saturated pure colour happens at l=0.5 and s=1 */
	s = 1.0; v = 1.0;
	hsv_to_rgb(h,s,v, &satu->red, &satu->green, &satu->blue);
	printf("%d: r: %g, g: %g, b: %g\n", layer, satu->red, satu->green, satu->blue);
}

// update a real colour from the saturated colour with a new lightness
void update_color_from_saturation(int layer, double newl) {
	double h, s, v;
	GdkRGBA *real = &layers[layer]->color;
	GdkRGBA *satu = &layers[layer]->saturated_color;
	rgb_to_hsv(satu->red, satu->green, satu->blue, &h,&s,&v);
	hsv_to_rgb(h,s,newl, &real->red, &real->green, &real->blue);
}

void on_colordialog_response(GtkColorChooserDialog *chooser, gint response_id, gpointer user_data) {
	/* this callback is called on any action of the dialog, and must be
	 * filtered according to the response_id. List is here:
	 * https://developer.gnome.org/gtk3/stable/GtkDialog.html#GTK-RESPONSE-NONE:CAPS
	 */
	if (response_id == GTK_RESPONSE_DELETE_EVENT ||
			response_id == GTK_RESPONSE_CANCEL ||
			response_id == GTK_RESPONSE_CLOSE) {
		current_layer_color_choosing = 0;
		gtk_widget_hide(GTK_WIDGET(chooser));
		gtk_editable_delete_text(GTK_EDITABLE(wl_entry), 0, -1);
		gtk_combo_box_set_active(GTK_COMBO_BOX(box), -1);
		return;
	}

	if (current_layer_color_choosing > 0 && layers[current_layer_color_choosing]) {
		gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(chooser), &layers[current_layer_color_choosing]->color);
		color_has_been_updated(current_layer_color_choosing);
		gtk_widget_queue_draw(GTK_WIDGET(layers[current_layer_color_choosing]->color_w));
		gtk_widget_hide(GTK_WIDGET(chooser));
		gtk_editable_delete_text(GTK_EDITABLE(wl_entry), 0, -1);
		gtk_combo_box_set_active(GTK_COMBO_BOX(box), -1);
		if (has_fit(current_layer_color_choosing))
			update_result(1);
	}
}

/* managing the drawing area that displays the color and opens the color dialog when clicked */
gboolean draw_layer_color(GtkDrawingArea *widget, cairo_t *cr, gpointer data) {
	int layer, w, h;
	for (layer = 0; layers[layer]; layer++)
		if (layers[layer]->color_w == widget)
			break;
	if (!layers[layer]) return FALSE;

	w = gtk_widget_get_allocated_width(GTK_WIDGET(widget));
	h = gtk_widget_get_allocated_height(GTK_WIDGET(widget));
	cairo_set_source_rgb(cr, layers[layer]->color.red,
			layers[layer]->color.green, layers[layer]->color.blue);
	//cairo_rectangle(cr, (double)w*0.33, 1, (double)w*0.33, (double)h-2.0);
	cairo_rectangle(cr, 1, 1, w-2.0, h-2.0);
	cairo_fill(cr);
	return FALSE;
}

/* click on the colored area: button press, only configure for quick color edit with
 * right button */
gboolean on_color_button_press_event(GtkDrawingArea *widget, GdkEventButton *event, gpointer user_data) {
	int layer;
	for (layer = 0; layers[layer]; layer++)
		if (layers[layer]->color_w == widget)
			break;
	if (!layers[layer]) return FALSE;
	if (event->button == 3) {	// right click
		current_layer_color_choosing = layer;
		color_quick_edit = 1;
		memcpy(&qe_ref_color, &layers[layer]->color, sizeof(GdkRGBA));
		return TRUE;
	}
	return FALSE;
}

/* click on the colored area: open the color chooser dialog */
gboolean on_color_button_release_event(GtkDrawingArea *widget, GdkEventButton *event, gpointer user_data) {
	int layer;
	for (layer = 0; layers[layer]; layer++)
		if (layers[layer]->color_w == widget)
			break;
	if (!layers[layer]) return FALSE;

	if (event->button == 1) {	// left click
		current_layer_color_choosing = layer;
		gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(color_dialog), &layers[layer]->color);
		gtk_editable_delete_text(GTK_EDITABLE(wl_entry), 0, -1);
		gtk_combo_box_set_active(GTK_COMBO_BOX(box), -1);
		gtk_widget_show(GTK_WIDGET(color_dialog));
	} else if (event->button == 3) {	// right click
		if (has_fit(current_layer_color_choosing))
			update_result(1);
		current_layer_color_choosing = 0;
		gtk_editable_delete_text(GTK_EDITABLE(wl_entry), 0, -1);
		gtk_combo_box_set_active(GTK_COMBO_BOX(box), -1);
	}
	color_quick_edit = 0;
	return TRUE;
}

gboolean on_color_button_motion_event(GtkWidget *widget, GdkEventMotion *event, gpointer user_data) {
	if (color_quick_edit) {	// right click
		double h, s, v;
		//fprintf(stdout, "%g\n", event->x);
		rgb_to_hsv(qe_ref_color.red, qe_ref_color.green,
				qe_ref_color.blue, &h,&s,&v);
		h += event->x / 600.0;
		v -= event->y / 600.0;
		while (h < 0.0) h += 1.0;
		while (h > 1.0) h -= 1.0;
		if (v < 0.0) v = 0.0;
		if (v > 1.0) v = 1.0;
		hsv_to_rgb(h,s,v, &layers[current_layer_color_choosing]->color.red,
				&layers[current_layer_color_choosing]->color.green,
				&layers[current_layer_color_choosing]->color.blue);
		color_has_been_updated(current_layer_color_choosing);
		gtk_widget_queue_draw(GTK_WIDGET(layers[current_layer_color_choosing]->color_w));
	}
	return FALSE;
}
/*******************************************************/

/* fill the combo box containing filter names */
void populate_filter_lists() {
	GtkComboBoxText *box;
	int i, nb_filters;
	box = GTK_COMBO_BOX_TEXT(gtk_builder_get_object(builder, "comboboxtext_filters"));
	nb_filters = get_nb_narrow_filters();
	gtk_combo_box_text_remove_all(box);
	for (i=0; i<nb_filters; i++)
		gtk_combo_box_text_append_text(box, narrow_band_filters[i].name);
	//gtk_combo_box_set_active(GTK_COMBO_BOX(box), -1);
}

/* the combo box containing filter names has one item selected: update colours */
void on_filter_changed(GtkComboBox *widget, gpointer user_data) {
	gint active = gtk_combo_box_get_active(widget);
	char wl_text[20];
	if (active == -1) return;
	sprintf(wl_text, "%g", narrow_band_filters[active].wavelength);
	gtk_entry_set_text(wl_entry, wl_text);
}

void on_wavelength_changed(GtkEditable *editable, gpointer user_data){
	GdkRGBA color;
	double wavelength = atof(gtk_entry_get_text(GTK_ENTRY(editable)));
	if (wavelength < 380.0 || wavelength > 780.0) return;
	wavelength_to_RGB(wavelength, &color);
	gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(color_dialog), &color);
}

void on_compositing_reset_clicked(GtkButton *button, gpointer user_data){
	int i;
	static GtkWidget *main_window = NULL;
	static GtkWidget *rgb_window = NULL;

	if (main_window == NULL) {
		main_window = lookup_widget("main_window");
		rgb_window = lookup_widget("rgb_window");
	}

	if (com.uniq) {
		close_single_image();
	}
	if (has_fit(0))
		clearfits(&layers[0]->the_fit);
	for (i=1; layers[i]; i++) {
		if (has_fit(i))
			clearfits(&layers[i]->the_fit);
		grid_remove_row(i, 1);
		free(layers[i]);
		layers[i] = NULL;
		layers_count--;
	}
	// already done in on_layer_add
	//gtk_container_remove(GTK_CONTAINER(grid_layers), GTK_WIDGET(add_button));

	/* Reset GtkFileChooserButton Luminance */
	GtkFileChooserButton *lum = GTK_FILE_CHOOSER_BUTTON(gtk_builder_get_object(builder, "filechooser_lum"));
	gtk_file_chooser_unselect_all (GTK_FILE_CHOOSER (lum));

	if (seq) {
		free_sequence(seq, TRUE);
		seq = NULL;
	}

	for (i=1; i<4; i++) {
		/* Create the three default layers */
		on_layer_add(NULL, NULL);
	}
	layers[i] = NULL;

	gtk_widget_hide(GTK_WIDGET(color_dialog));
	current_layer_color_choosing = 0;

	luminance_mode = 0;
	GtkToggleButton *lum_button = GTK_TOGGLE_BUTTON(gtk_builder_get_object(builder, "composition_use_lum"));
	gtk_toggle_button_set_active(lum_button, 0);
	gtk_label_set_text(layers[0]->label, _("not loaded"));


	update_compositing_interface();
	open_compositing_window();	// update the CWD just in case

	/* Close all oppened windows */
	if (gtk_widget_get_visible(main_window))
		gtk_widget_hide(main_window);
	if (gtk_widget_get_visible(rgb_window))
		gtk_widget_hide(rgb_window);

	update_used_memory();
}

/* Reduce brightness of colours associated to layers so that they never overflow on composition.
 * Algorithm: take the max of the composition and normalize brightness with
 * this max. It has to be done three times because it's the same layers that
 * can act on the resulting RGB channels.
 * This algorithm doesn't give the optimal answer, which could be found
 * iteratively, but it should never give an overflow.
 */
void autoadjust(int force_redraw) {
	int layer, nb_images_red = 0, nb_images_green = 0, nb_images_blue = 0;
	GdkRGBA max_pixel;

	set_cursor_waiting(TRUE);
	clear_pixel(&max_pixel);
	/* sum the max per channel */
	/* should we assume that fits mini and maxi are correct? */
	for (layer = 1; layers[layer]; layer++) {
		if (has_fit(layer)) {
			WORD max_value = layers[layer]->the_fit.maxi;
			if (coeff)
				max_value = get_normalized_pixel_value(layer, max_value);
			increment_pixel_components_from_layer_saturated_value(
					layer, &max_pixel, max_value);

			if (layers[layer]->color.red > 0.0) nb_images_red++;
			if (layers[layer]->color.green > 0.0) nb_images_green++;
			if (layers[layer]->color.blue > 0.0) nb_images_blue++;
		}
	}

	if (max_pixel.red <= 1.0 && max_pixel.green <= 1.0 && max_pixel.blue <= 1.0) {
		if (force_redraw) {
			siril_log_message(_("No overflow with the current colours, redrawing only\n"));
			update_result(1);
		} else {
			siril_log_message(_("Nothing to adjust, no overflow\n"));
			set_cursor_waiting(FALSE);
		}
		return;
	}

	/* update the real colours of layers from their saturated colour, based
	 * on how much each colour of the composition overflows */
	// amounts of normalization to be done on each layer's image for each channel
	double to_redistribute_red = (max_pixel.red - 1.0) / (double)nb_images_red;
	double to_redistribute_green = (max_pixel.green - 1.0) / (double)nb_images_green;
	double to_redistribute_blue = (max_pixel.blue - 1.0) / (double)nb_images_blue;
	for (layer = 1; layers[layer]; layer++) {
		if (has_fit(layer)) {
			double to_redistribute = 0.0;	// for this layer

			if (layers[layer]->color.red > 0.0 && to_redistribute_red > 0.0) {
				to_redistribute = to_redistribute_red;
			}
			/* for each layer, we check if a channel requires the
			 * current layer to be readjusted and we take the more
			 * severe value of all requirements */

			if (layers[layer]->color.green > 0.0 && to_redistribute_green > 0.0) {
				if (to_redistribute_green > to_redistribute)
					to_redistribute = to_redistribute_green;
			}

			if (layers[layer]->color.blue > 0.0 && to_redistribute_blue > 0.0) {
				if (to_redistribute_blue > to_redistribute)
					to_redistribute = to_redistribute_blue;
			}

			siril_log_message(_("Readjusting layer %d to %g times bright\n"),
					layer, 1.0-to_redistribute);
			/* to_redistribute here is the maximum reduction we
			 * need to give to the layer */
			update_color_from_saturation(layer, 1.0 - to_redistribute);
		}
	}

	/* redraw colours and composition */
	for (layer = 1; layers[layer]; layer++) {
		gtk_widget_queue_draw(GTK_WIDGET(layers[layer]->color_w));
	}
	update_result(1);
	set_cursor_waiting(FALSE);
}

void on_compositing_autoadjust_clicked(GtkButton *button, gpointer user_data){
	autoadjust(0);
}

/* Normalization functions, not used anymore */

void coeff_clear() {
	if (coeff) {
		free(coeff->offset);
		free(coeff->scale);
		free(coeff->mul);
		free(coeff);
		coeff = NULL;
	}
}

void on_composition_rgbcolor_clicked(GtkButton *button, gpointer user_data){
	initialize_calibration_interface();
	gtk_widget_show(lookup_widget("color_calibration"));
}


