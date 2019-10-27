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
#include "core/command.h"	// for processcommand()
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

#include "callbacks.h"
#include "message_dialog.h"
#include "PSF_list.h"
#include "histogram.h"
#include "script_menu.h"
#include "progress_and_log.h"
#include "dialogs.h"
#include "plot.h"

static gboolean is_shift_on = FALSE;

layer_info predefined_layers_colors[] = {
		/* name, lambda, lo, hi, c/over, c/under, mode */
		{ N_("Luminance"), 0.0, 0, 0, FALSE, FALSE, NORMAL_DISPLAY }, // no color, undefined value is <0
		{ N_("Red"), 650.0, 0, 0, FALSE, FALSE, NORMAL_DISPLAY }, // approx. of the middle of the color
		{ N_("Green"), 530.0, 0, 0, FALSE, FALSE, NORMAL_DISPLAY },	// approx. of the middle of the color
		{ N_("Blue"), 450.0, 0, 0, FALSE, FALSE, NORMAL_DISPLAY }// approx. of the middle of the color
};

/* remap index data, an index for each layer */
static BYTE *remap_index[MAXGRAYVPORT];
static float last_pente[MAXGRAYVPORT];
static display_mode last_mode[MAXGRAYVPORT];

/* STF (auto-stretch) data */
static gboolean stfComputed;	// Flag to know if STF parameters are available
static double stfShadows, stfHighlights, stfM;

/*****************************************************************************
 *                    S T A T I C      F U N C T I O N S                     *
 ****************************************************************************/

/* this function calculates the "fit to window" zoom values, given the window
 * size in argument and the image size in gfit.
 * Should not be called before displaying the main gray window when using zoom to fit */
static double get_zoom_val() {
	int window_width, window_height;
	double wtmp, htmp;
	static GtkWidget *scrolledwin = NULL;
	if (scrolledwin == NULL)
		scrolledwin = lookup_widget("scrolledwindowr");
	if (com.zoom_value > 0.)
		return com.zoom_value;
	/* else if zoom is < 0, it means fit to window */
	window_width = gtk_widget_get_allocated_width(scrolledwin);
	window_height = gtk_widget_get_allocated_height(scrolledwin);
	if (gfit.rx == 0 || gfit.ry == 0 || window_height <= 1 || window_width <= 1)
		return 1.0;
	wtmp = (double) window_width / (double) gfit.rx;
	htmp = (double) window_height / (double) gfit.ry;
	//fprintf(stdout, "computed fit to window zooms: %f, %f\n", wtmp, htmp);
	return min(wtmp, htmp);
}

/*
 * Returns the modifier mask. For Linux it is Control key, but for Apple - OS X
 * it should be Apple Key -> Mod2
 */
static GdkModifierType get_default_modifier() {
	GdkDisplay *display = gdk_display_get_default();
	GdkKeymap *keymap = gdk_keymap_get_for_display(display);
	GdkModifierType primary, real;

	g_return_val_if_fail(GDK_IS_KEYMAP (keymap), 0);

	/* Retrieve the real modifier mask */
	real = gdk_keymap_get_modifier_mask(keymap,
			GDK_MODIFIER_INTENT_PRIMARY_ACCELERATOR);

	primary = real;

	/* We need to translate the real modifiers into a virtual modifier
	 (like Super, Meta, etc.).
	 The following call adds the virtual modifiers for each real modifier
	 defined in primary.
	 */
	gdk_keymap_add_virtual_modifiers(keymap, &primary);

	if (primary != real) {
		/* In case the virtual and real modifiers are different, we need to
		 remove the real modifier from the result, and keep only the
		 virtual one.
		 */
		primary &= ~real;
	}
	return primary;
}

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

/*
 * Image display static functions
 */

static gboolean inimage(GdkEvent *event) {
	double zoom = get_zoom_val();
	return ((GdkEventButton*) event)->x > 0
			&& ((GdkEventButton*) event)->x < gfit.rx * zoom
			&& ((GdkEventButton*) event)->y > 0
			&& ((GdkEventButton*) event)->y < gfit.ry * zoom;
}

static void draw_empty_image(cairo_t *cr, guint width, guint height) {
	// black image with red square
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_rectangle(cr, 0, 0, width, height);
	cairo_fill(cr);
	cairo_set_source_rgb(cr, 0.3, 0, 0);
	cairo_rectangle(cr, 100, 70, 50, 50);
	cairo_fill(cr);
}


static void remaprgb(void) {
	guchar *dst;
	guchar *bufr, *bufg, *bufb;
	gint i, j;
	int nbdata;

	fprintf(stderr, "remaprgb\n");
	if (!isrgb(&gfit))
		return;

	// allocate if not already done or the same size
	if (cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, gfit.rx)
			!= com.surface_stride[RGB_VPORT]
			|| gfit.ry != com.surface_height[RGB_VPORT]
			|| !com.surface[RGB_VPORT] || !com.rgbbuf) {
		guchar *oldbuf = com.rgbbuf;
		fprintf(stderr, "RGB display buffers and surface (re-)allocation\n");
		com.surface_stride[RGB_VPORT] = cairo_format_stride_for_width(
				CAIRO_FORMAT_RGB24, gfit.rx);
		com.surface_height[RGB_VPORT] = gfit.ry;
		com.rgbbuf = realloc(com.rgbbuf,
				com.surface_stride[RGB_VPORT] * gfit.ry * sizeof(guchar));
		if (com.rgbbuf == NULL) {
			PRINT_ALLOC_ERR;
			if (oldbuf)
				free(oldbuf);
			return;
		}
		if (com.surface[RGB_VPORT])
			cairo_surface_destroy(com.surface[RGB_VPORT]);
		com.surface[RGB_VPORT] = cairo_image_surface_create_for_data(com.rgbbuf,
				CAIRO_FORMAT_RGB24, gfit.rx, gfit.ry,
				com.surface_stride[RGB_VPORT]);
		if (cairo_surface_status(com.surface[RGB_VPORT])
				!= CAIRO_STATUS_SUCCESS) {
			fprintf(stderr,
					"Error creating the Cairo image surface for the RGB image\n");
			cairo_surface_destroy(com.surface[RGB_VPORT]);
			com.surface[RGB_VPORT] = NULL;
			return;
		}
	}
	// WARNING : this assumes that R, G and B buffers are already allocated and mapped
	// it seems ok, but one can probably imagine situations where it segfaults
	bufr = com.graybuf[RED_VPORT];
	bufg = com.graybuf[GREEN_VPORT];
	bufb = com.graybuf[BLUE_VPORT];
	if (bufr == NULL || bufg == NULL || bufb == NULL) {
		fprintf(stderr, "remaprgb: gray buffers not allocated for display\n");
		return;
	}
	dst = com.rgbbuf;	// index is j
	nbdata = gfit.rx * gfit.ry * 4;	// source images are 32-bit RGBA

	for (i = 0, j = 0; i < nbdata; i += 4) {
		dst[j++] = bufb[i];
		dst[j++] = bufg[i];
		dst[j++] = bufr[i];
		j++;		// alpha padding
	}

	// flush to ensure all writing to the image was done and redraw the surface
	cairo_surface_flush(com.surface[RGB_VPORT]);
	cairo_surface_mark_dirty(com.surface[RGB_VPORT]);
}

static void set_viewer_mode_widgets_sensitive(gboolean sensitive) {
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

static int make_index_for_current_display(display_mode mode, WORD lo, WORD hi,
		int vport);
static int make_index_for_rainbow(BYTE index[][3]);

/* enables or disables the "display reference" checkbox in registration preview */
static void enable_view_reference_checkbox(gboolean status) {
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

/* vport can be -1 if the correct viewport should be tested */
static void test_and_allocate_reference_image(int vport) {
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

static void remap(int vport) {
	// This function maps fit data with a linear LUT between lo and hi levels
	// to the buffer to be displayed; display only is modified
	guint y;
	BYTE *dst, *index, rainbow_index[UCHAR_MAX + 1][3];
	WORD *src, hi, lo;
	display_mode mode;
	color_map color;
	gboolean do_cut_over, inverted;

	fprintf(stderr, "remap %d\n", vport);
	if (vport == RGB_VPORT) {
		remaprgb();
		return;
	}

	int no_data = 0;
	if (single_image_is_loaded()) {
	       if (vport >= com.uniq->nb_layers)
		       no_data = 1;
	}
	else if (sequence_is_loaded()) {
		if (vport >= com.seq.nb_layers)
			no_data = 1;
	}
	else no_data = 1;
	if (no_data) {
		fprintf(stderr, "vport is out of bounds or data is not loaded yet\n");
		return;
	}

	// allocate if not already done or the same size
	if (cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, gfit.rx) !=
			com.surface_stride[vport] ||
			gfit.ry != com.surface_height[vport] ||
			!com.surface[vport] ||
			!com.graybuf[vport]) {
		guchar *oldbuf = com.graybuf[vport];
		fprintf(stderr, "Gray display buffers and surface (re-)allocation\n");
		if (gfit.rx == 0 || gfit.ry == 0) {
			fprintf(stderr, "gfit has a zero size, must not happen!\n");
			return;
		}
		com.surface_stride[vport] = cairo_format_stride_for_width(
				CAIRO_FORMAT_RGB24, gfit.rx);
		com.surface_height[vport] = gfit.ry;
		com.graybuf[vport] = realloc(com.graybuf[vport],
				com.surface_stride[vport] * gfit.ry * sizeof(guchar));
		if (com.graybuf[vport] == NULL) {
			PRINT_ALLOC_ERR;
			if (oldbuf)
				free(oldbuf);
			return;
		}
		if (com.surface[vport])
			cairo_surface_destroy(com.surface[vport]);
		com.surface[vport] = cairo_image_surface_create_for_data(
				com.graybuf[vport], CAIRO_FORMAT_RGB24, gfit.rx, gfit.ry,
				com.surface_stride[vport]);
		if (cairo_surface_status(com.surface[vport]) != CAIRO_STATUS_SUCCESS) {
			fprintf(stderr,
					"Error creating the Cairo image surface for vport %d\n",
					vport);
			cairo_surface_destroy(com.surface[vport]);
			com.surface[vport] = NULL;
			return;
		}
	}
	if (single_image_is_loaded() && com.seq.current != RESULT_IMAGE) {
		mode = com.uniq->layers[vport].rendering_mode;
		hi = com.uniq->layers[vport].hi;
		lo = com.uniq->layers[vport].lo;
		do_cut_over = com.uniq->layers[vport].cut_over;
	} else if (sequence_is_loaded() && vport < com.seq.nb_layers) {
		// the check above is needed because there may be a different
		// number of channels between the unique image and the sequence
		mode = com.seq.layers[vport].rendering_mode;
		hi = com.seq.layers[vport].hi;
		lo = com.seq.layers[vport].lo;
		do_cut_over = com.seq.layers[vport].cut_over;
	} else {
		fprintf(stderr, "BUG in unique image remap\n");
		return;
	}

	if (lo > hi) {
		// negative display
		WORD tmp = hi;
		hi = lo;
		lo = tmp;
		inverted = TRUE;
	} else
		inverted = FALSE;

	if (mode == HISTEQ_DISPLAY) {
		double hist_sum;
		double nb_pixels;
		size_t hist_nb_bins;
		size_t i;
		gsl_histogram *histo;

		compute_histo_for_gfit();
		histo = com.layers_hist[vport];
		hist_nb_bins = gsl_histogram_bins(histo);
		/*if (hist_nb_bins <= USHRT_MAX) {
		 fprintf(stderr, "Error remapping: histogram is not the correct size\n");
		 return;
		 }*/
		nb_pixels = (double) (gfit.rx * gfit.ry);
		// build the remap_index
		if (!remap_index[vport])
			remap_index[vport] = malloc(USHRT_MAX + 1);

		remap_index[vport][0] = 0;
		hist_sum = gsl_histogram_get(histo, 0);
		for (i = 1; i < hist_nb_bins; i++) {
			hist_sum += gsl_histogram_get(histo, i);
			remap_index[vport][i] = round_to_BYTE(
					(hist_sum / nb_pixels) * UCHAR_MAX_DOUBLE);
		}

		last_mode[vport] = mode;
		set_viewer_mode_widgets_sensitive(FALSE);
	} else {
		// for all other modes, the index can be reused
		make_index_for_current_display(mode, lo, hi, vport);
		if (mode == STF_DISPLAY)
			set_viewer_mode_widgets_sensitive(FALSE);
		else
			set_viewer_mode_widgets_sensitive(TRUE);
	}

	src = gfit.pdata[vport];
	/* Siril's FITS are stored bottom to top, so mapping needs to revert data order */
	dst = com.graybuf[vport];

	color = gtk_toggle_tool_button_get_active(
			GTK_TOGGLE_TOOL_BUTTON(lookup_widget("colormap_button")));

	if (color == RAINBOW_COLOR)
		make_index_for_rainbow(rainbow_index);
	index = remap_index[vport];

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(y) schedule(static)
#endif
	for (y = 0; y < gfit.ry; y++) {
		guint x;
		for (x = 0; x < gfit.rx; x++) {
			guint src_index = y * gfit.rx + x;
			BYTE dst_pixel_value;
			WORD tmp_pixel_value;
			if (mode == HISTEQ_DISPLAY || mode == STF_DISPLAY)	// special case, no lo & hi
				dst_pixel_value = index[src[src_index]];
			else if (do_cut_over && src[src_index] > hi)	// cut
				dst_pixel_value = 0;
			else {
				if (src[src_index] - lo < 0)
					tmp_pixel_value = 0;
				else
					tmp_pixel_value = src[src_index] - lo;
				dst_pixel_value = index[tmp_pixel_value];
			}
			if (inverted)
				dst_pixel_value = UCHAR_MAX - dst_pixel_value;

			guint dst_index = ((gfit.ry - 1 - y) * gfit.rx + x) * 4;
			switch (color) {
				default:
				case NORMAL_COLOR:
					dst[dst_index++] = dst_pixel_value;
					dst[dst_index++] = dst_pixel_value;
					dst[dst_index++] = dst_pixel_value;
					break;
				case RAINBOW_COLOR:
					dst[dst_index++] = rainbow_index[dst_pixel_value][0];
					dst[dst_index++] = rainbow_index[dst_pixel_value][1];
					dst[dst_index++] = rainbow_index[dst_pixel_value][2];
			}
		}
	}

	// flush to ensure all writing to the image was done and redraw the surface
	cairo_surface_flush(com.surface[vport]);
	cairo_surface_mark_dirty(com.surface[vport]);

	test_and_allocate_reference_image(vport);
}

static int make_index_for_current_display(display_mode mode, WORD lo, WORD hi,
		int vport) {
	float pente;
	int i;
	BYTE *index;
	double pxl;
	if (mode == STF_DISPLAY) {
		if (!stfComputed) {
			stfM = findMidtonesBalance(&gfit, &stfShadows, &stfHighlights);
			stfComputed = TRUE;
		}
	}

	/* initialization of data required to build the remap_index */
	switch (mode) {
	case NORMAL_DISPLAY:
		pente = UCHAR_MAX_SINGLE / (float) (hi - lo);
		break;
	case LOG_DISPLAY:
		pente = fabsf(UCHAR_MAX_SINGLE / logf(((float) (hi - lo)) * 0.1f));
		break;
	case SQRT_DISPLAY:
		pente = UCHAR_MAX_SINGLE / sqrtf((float) (hi - lo));
		break;
	case SQUARED_DISPLAY:
		pente = UCHAR_MAX_SINGLE / SQR((float )(hi - lo));
		break;
	case ASINH_DISPLAY:
		pente = UCHAR_MAX_SINGLE / asinhf(((float) (hi - lo)) * 0.001f);
		break;
	case STF_DISPLAY:
		pente = UCHAR_MAX_SINGLE;
		break;
	default:
		return 1;
	}
	if ((mode != HISTEQ_DISPLAY && mode != STF_DISPLAY) && pente == last_pente[vport]
			&& mode == last_mode[vport]) {
		fprintf(stdout, "Re-using previous remap_index\n");
		return 0;
	}
	fprintf(stdout, "Rebuilding remap_index\n");

	/************* Building the remap_index **************/
	if (!remap_index[vport]) {
		remap_index[vport] = malloc(USHRT_MAX + 1);
		if (!remap_index[vport]) {
			fprintf(stderr,
					"allocation error in remap_index, aborting remap\n");
			return 1;
		}
	}
	index = remap_index[vport];

	for (i = 0; i <= USHRT_MAX; i++) {
		switch (mode) {
		case LOG_DISPLAY:
			// ln(5.56*10^110) = 255
			if (i < 10)
				index[i] = 0; /* avoid null and negative values */
			else
				index[i] = round_to_BYTE(logf((float) i / 10.f) * pente); //10.f is arbitrary: good matching with ds9
			break;
		case SQRT_DISPLAY:
			// sqrt(2^16) = 2^8
			index[i] = round_to_BYTE(sqrtf((float) i) * pente);
			break;
		case SQUARED_DISPLAY:
			// pow(2^4,2) = 2^8
			index[i] = round_to_BYTE(SQR((float)i) * pente);
			break;
		case ASINH_DISPLAY:
			// asinh(2.78*10^110) = 255
			index[i] = round_to_BYTE(asinhf((float) i / 1000.f) * pente); //1000.f is arbitrary: good matching with ds9, could be asinhf(a*Q*i)/Q
			break;
		case NORMAL_DISPLAY:
			index[i] = round_to_BYTE((float) i * pente);
			break;
		case STF_DISPLAY:
			pxl = (gfit.orig_bitpix == BYTE_IMG ?
					(double) i / UCHAR_MAX_DOUBLE :
					(double) i / USHRT_MAX_DOUBLE);
			index[i] = round_to_BYTE((float) (MTF(pxl, stfM, stfShadows, stfHighlights)) * pente);
			break;
		default:
			return 1;
		}
		// check for maximum overflow, given that df/di > 0. Should not happen with round_to_BYTE
		if (index[i] == UCHAR_MAX)
			break;
	}
	if (i != USHRT_MAX + 1) {
		/* no more computation needed, just fill with max value */
		for (++i; i <= USHRT_MAX; i++)
			index[i] = UCHAR_MAX;
	}

	last_pente[vport] = pente;
	last_mode[vport] = mode;
	return 0;
}

static int make_index_for_rainbow(BYTE index[][3]) {
	int i;
	double h, s, v, r, g, b;

	for (i = 0; i < UCHAR_MAX + 1; i++) {
		r = g = b = (double) i / UCHAR_MAX_DOUBLE;
		rgb_to_hsv(r, g, b, &h, &s, &v);
		double off = 300.0 / 360.0;  /* Arbitrary: we want h from 300 to 0 deg */
		h = (off - (double) i * (off / UCHAR_MAX_DOUBLE));
		s = 1.;
		v = 1.; /* Saturation and Value are set to 100%  */
		hsv_to_rgb(h, s, v, &r, &g, &b);
		index[i][0] = round_to_BYTE(r * UCHAR_MAX_DOUBLE);
		index[i][1] = round_to_BYTE(g * UCHAR_MAX_DOUBLE);
		index[i][2] = round_to_BYTE(b * UCHAR_MAX_DOUBLE);
	}
	return 0;
}

/*
 * Free reference image
 */

static void free_reference_image() {
	fprintf(stdout, "Purging previously saved reference frame data.\n");
	if (com.refimage_regbuffer) {
		free(com.refimage_regbuffer);
		com.refimage_regbuffer = NULL;
	}
	if (com.refimage_surface) {
		cairo_surface_destroy(com.refimage_surface);
		com.refimage_surface = NULL;
	}
	if (com.seq.reference_image == -1)
		enable_view_reference_checkbox(FALSE);
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

/*
 * Command line history static function
 */

static void history_add_line(char *line) {
	if (!com.cmd_history) {
		com.cmd_hist_size = CMD_HISTORY_SIZE;
		com.cmd_history = calloc(com.cmd_hist_size, sizeof(const char*));
		com.cmd_hist_current = 0;
		com.cmd_hist_display = 0;
	}
	com.cmd_history[com.cmd_hist_current] = line;
	com.cmd_hist_current++;
	// circle at the end
	if (com.cmd_hist_current == com.cmd_hist_size)
		com.cmd_hist_current = 0;
	if (com.cmd_history[com.cmd_hist_current]) {
		free(com.cmd_history[com.cmd_hist_current]);
		com.cmd_history[com.cmd_hist_current] = NULL;
	}
	com.cmd_hist_display = com.cmd_hist_current;
}

/*
 * Selection static functions
 */

/* selection zone event management */
#define MAX_CALLBACKS_PER_EVENT 10
static selection_update_callback _registered_callbacks[MAX_CALLBACKS_PER_EVENT];
static int _nb_registered_callbacks = 0;

void register_selection_update_callback(selection_update_callback f) {
	if (_nb_registered_callbacks < MAX_CALLBACKS_PER_EVENT) {
		_registered_callbacks[_nb_registered_callbacks] = f;
		_nb_registered_callbacks++;
	}
}

void unregister_selection_update_callback(selection_update_callback f) {
	int i;
	for (i = 0; i < _nb_registered_callbacks; ++i) {
		if (_registered_callbacks[i] == f) {
			_registered_callbacks[i] =
					_registered_callbacks[_nb_registered_callbacks];
			_registered_callbacks[_nb_registered_callbacks] = NULL;
			_nb_registered_callbacks--;
			return;
		}
	}
}

// send the events
static void new_selection_zone() {
	int i;
	fprintf(stdout, "selection: %d,%d,\t%dx%d\n", com.selection.x, com.selection.y,
			com.selection.w, com.selection.h);
	for (i = 0; i < _nb_registered_callbacks; ++i) {
		_registered_callbacks[i]();
	}
	redraw(com.cvport, REMAP_NONE);
}

void delete_selected_area() {
	memset(&com.selection, 0, sizeof(rectangle));
	new_selection_zone();
}

/*
 * MISC static functions
 */

static void toggle_image_selection(int image_num) {
	gchar *msg;
	if (com.seq.imgparam[image_num].incl) {
		com.seq.imgparam[image_num].incl = FALSE;
		--com.seq.selnum;
		msg = g_strdup_printf(_("Image %d has been unselected from sequence\n"), image_num);
		if (image_num == com.seq.reference_image) {
			com.seq.reference_image = -1;
			sequence_list_change_reference();
			adjust_refimage(image_num);
		}
	} else {
		com.seq.imgparam[image_num].incl = TRUE;
		++com.seq.selnum;
		msg = g_strdup_printf(_("Image %d has been selected from sequence\n"), image_num);
	}
	siril_log_message(msg);
	g_free(msg);
	sequence_list_change_selection_index(image_num);
	update_reg_interface(FALSE);
	update_stack_interface(TRUE);
	adjust_exclude(image_num, TRUE);
	writeseqfile(&com.seq);
}

/* method handling all include or all exclude from a sequence */
static void sequence_setselect_all(gboolean include_all) {
	int i;

	if (!com.seq.imgparam)
		return;
	for (i = 0; i < com.seq.number; ++i) {
		if (com.seq.imgparam[i].incl != include_all) {
			com.seq.imgparam[i].incl = include_all;
			sequence_list_change_selection_index(i);
		}
	}
	if (include_all) {
		com.seq.selnum = com.seq.number;
		siril_log_message(_("Selected all images from sequence\n"));
	} else {
		com.seq.selnum = 0;
		com.seq.reference_image = -1;
		siril_log_message(_("Unselected all images from sequence\n"));
		sequence_list_change_reference();
		adjust_refimage(com.seq.current);
	}
	adjust_exclude(com.seq.current, TRUE);
	update_reg_interface(FALSE);
	update_stack_interface(TRUE);
	writeseqfile(&com.seq);
}

/* RGB popup menu */

static void do_popup_rgbmenu(GtkWidget *my_widget, GdkEventButton *event) {
	static GtkMenu *menu = NULL;

	if (!menu) {
		menu = GTK_MENU(gtk_builder_get_object(builder, "menurgb"));
		gtk_menu_attach_to_widget(GTK_MENU(menu), my_widget, NULL);
	}

#if GTK_CHECK_VERSION(3, 22, 0)
	gtk_menu_popup_at_pointer(GTK_MENU(menu), NULL);
#else
	int button, event_time;

	if (event) {
		button = event->button;
		event_time = event->time;
	} else {
		button = 0;
		event_time = gtk_get_current_event_time();
	}

	gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL, button,
			event_time);
#endif
}

/* Gray popup menu */

static void do_popup_graymenu(GtkWidget *my_widget, GdkEventButton *event) {
	static GtkMenu *menu = NULL;
	gboolean selected;

	gboolean is_a_single_image_loaded = single_image_is_loaded()	&& (!sequence_is_loaded()
			|| (sequence_is_loaded() && (com.seq.current == RESULT_IMAGE
					|| com.seq.current == SCALED_IMAGE)));

	if (!menu) {
		menu = GTK_MENU(gtk_builder_get_object(builder, "menugray"));
		gtk_menu_attach_to_widget(GTK_MENU(menu), my_widget, NULL);
	}

	selected = com.selection.w && com.selection.h;
	gtk_widget_set_sensitive(lookup_widget("undo_item1"), is_undo_available());
	gtk_widget_set_sensitive(lookup_widget("redo_item1"), is_redo_available());
	gtk_widget_set_sensitive(lookup_widget("menu_gray_psf"), selected);
	gtk_widget_set_sensitive(lookup_widget("menu_gray_seqpsf"), selected);
	gtk_widget_set_sensitive(lookup_widget("menu_gray_pick_star"), selected);
	gtk_widget_set_sensitive(lookup_widget("menu_gray_crop"), selected && is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menu_gray_crop_seq"), selected && sequence_is_loaded());

#if GTK_CHECK_VERSION(3, 22, 0)
	gtk_menu_popup_at_pointer(GTK_MENU(menu), NULL);
#else
	int button, event_time;

	if (event) {
		button = event->button;
		event_time = event->time;
	} else {
		button = 0;
		event_time = gtk_get_current_event_time();
	}

	gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL, button, event_time);
#endif
}

/*****************************************************************************
 *                    P U B L I C      F U N C T I O N S                     *
 ****************************************************************************/

GtkWidget* lookup_widget(const gchar *widget_name) {
	return GTK_WIDGET(gtk_builder_get_object(builder, widget_name));
}

GtkWidget *popover_new(GtkWidget *widget, const gchar *text) {
	GtkWidget *popover, *box, *image, *label;

	popover = gtk_popover_new(widget);
	label = gtk_label_new(NULL);
	box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	image = gtk_image_new_from_icon_name("dialog-information-symbolic",
			GTK_ICON_SIZE_DIALOG);

	gtk_label_set_markup(GTK_LABEL(label), text);
	gtk_label_set_line_wrap(GTK_LABEL(label), TRUE);
	gtk_label_set_max_width_chars(GTK_LABEL(label), 64);

	gtk_box_pack_start(GTK_BOX(box), image, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(box), label, FALSE, FALSE, 0);
	gtk_container_add(GTK_CONTAINER(popover), box);

	gtk_widget_show_all(box);

	return popover;
}

static void update_theme_button(const gchar *button_name, const gchar *path) {
	gchar *images;

	images = g_build_filename(com.app_path, "pixmaps", path, NULL);
	gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(lookup_widget(button_name)),
			gtk_image_new_from_file(images));
	gtk_widget_show_all(lookup_widget(button_name));

	g_free(images);
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
	if ((com.combo_theme == 1) || (com.have_dark_theme && com.combo_theme == 0)) {
		com.want_dark = TRUE;
	} else {
		com.want_dark = FALSE;
	}

	settings = gtk_settings_get_default();
	g_object_set(settings, "gtk-application-prefer-dark-theme", com.want_dark, NULL);
	update_icons_to_theme(com.want_dark);
}

static void initialize_theme_GUI() {
	GtkComboBox *box;

	box = GTK_COMBO_BOX(lookup_widget("combo_theme"));

	g_signal_handlers_block_by_func(box, on_combo_theme_changed, NULL);
	gtk_combo_box_set_active(box, com.combo_theme);
	g_signal_handlers_unblock_by_func(box, on_combo_theme_changed, NULL);
	update_icons_to_theme(com.want_dark);
}

void load_prefered_theme(gint theme) {
	GtkSettings *settings;

	settings = gtk_settings_get_default();
	g_object_get(settings, "gtk-application-prefer-dark-theme", &com.have_dark_theme,
				NULL);

	if ((theme == 1) || (com.have_dark_theme && theme == 0)) {
		com.want_dark = TRUE;
	} else {
		com.want_dark = FALSE;
	}

	g_object_set(settings, "gtk-application-prefer-dark-theme", com.want_dark, NULL);
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
	gtk_adjustment_set_upper(adjmin, (double) maxvalue);
	gtk_adjustment_set_upper(adjmax, (double) maxvalue);
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

/* sets the maximum value for the spin button and display the initial file name */
int seqsetnum(int image_number) {
	GtkSpinButton *spin;
	GtkAdjustment *adj;
	if (com.seq.number <= 0 || image_number >= com.seq.number)
		return 1;
	spin = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "imagenumber_spin"));
	adj = gtk_spin_button_get_adjustment(spin);

	gtk_adjustment_set_upper(adj, (gdouble) com.seq.number - 1);
	gtk_adjustment_set_value(adj, (gdouble) image_number);	// 0 is default
	display_image_number(image_number);
	//on_imagenumberspin_output(GTK_SPIN_BUTTON(spin), NULL);	// redraw the real file number instead of 0
	return 0;
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

/* called when an image is included or excluded from the sequence.
 * it sets the toggle button "exclude_button" so it matches the inclusion state.
 * n is the image number in the sequence
 * changed indicates if a value was changed and if the display needs to be refreshed.
 */
void adjust_exclude(int n, gboolean changed) {
	static GtkWidget *excl_butt = NULL;
	if (!com.seq.imgparam || n < 0 || n >= com.seq.number)
		return;
	if (excl_butt == NULL)
		excl_butt = lookup_widget("exclude_button");

	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(excl_butt))
			== com.seq.imgparam[n].incl) {
		g_signal_handlers_block_by_func(excl_butt, on_excludebutton_toggled,
				NULL);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(excl_butt),
				!com.seq.imgparam[n].incl);
		g_signal_handlers_unblock_by_func(excl_butt, on_excludebutton_toggled,
				NULL);
	}

	if (changed) {
		redraw(com.cvport, REMAP_NONE);
		drawPlot();
		adjust_sellabel();
	}
}

/* fill the label indicating how many images are selected in the gray and
 * which one is the reference image, at the bottom of the main window */
int adjust_sellabel() {
	static GtkLabel *local_label = NULL, *global_label = NULL;
	char bufferlocal[256], bufferglobal[256];
	gchar *seq_basename = NULL;

	if (local_label == NULL) {
		local_label = GTK_LABEL(lookup_widget("imagesel_label"));
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
			g_snprintf(bufferlocal, sizeof(bufferlocal), format,
					seq_basename, com.seq.selnum, com.seq.number,
					com.seq.imgparam[com.seq.reference_image].filenum);

		} else {
			g_snprintf(bufferlocal, sizeof(bufferlocal),
					_("<%s.seq>: %d images selected out of %d, no reference image set"),
					seq_basename, com.seq.selnum, com.seq.number);
		}
		g_snprintf(bufferglobal, sizeof(bufferglobal), _("%s, %d images selected"),
				seq_basename, com.seq.selnum);
		//gtk_widget_set_sensitive(lookup_widget("goregister_button"), com.seq.selnum>0?TRUE:FALSE);
	} else {
		g_snprintf(bufferlocal, sizeof(bufferlocal), _("No sequence"));
		g_snprintf(bufferglobal, sizeof(bufferglobal), _("- none -"));
		gtk_widget_set_sensitive(lookup_widget("goregister_button"), FALSE);
	}

	gtk_label_set_text(local_label, bufferlocal);
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
	gboolean any_RGB_image_is_loaded;		/* Some RGB data are loaded. Single image or Sequence */

	is_a_singleRGB_image_loaded = isrgb(&gfit) && (!sequence_is_loaded()
			|| (sequence_is_loaded() && (com.seq.current == RESULT_IMAGE
					|| com.seq.current == SCALED_IMAGE)));

	is_a_single_image_loaded = single_image_is_loaded()	&& (!sequence_is_loaded()
			|| (sequence_is_loaded() && (com.seq.current == RESULT_IMAGE
					|| com.seq.current == SCALED_IMAGE)));

	any_image_is_loaded = single_image_is_loaded() || sequence_is_loaded();

	any_RGB_image_is_loaded = isrgb(&gfit) && (single_image_is_loaded() || sequence_is_loaded());

	/* File Menu */
	gtk_widget_set_sensitive(lookup_widget("save1"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menu_FITS_header"), any_image_is_loaded && gfit.header != NULL);
	gtk_widget_set_sensitive(lookup_widget("menu_file_information"), is_a_single_image_loaded);

	/* Edit Menu */
	gtk_widget_set_sensitive(lookup_widget("undo_item"), is_undo_available());
	gtk_widget_set_sensitive(lookup_widget("redo_item"), is_redo_available());

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

	/* Analysis Menu */
	gtk_widget_set_sensitive(lookup_widget("menuitem_noise"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_stat"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_IPS"), any_image_is_loaded);

	/* Windows Menu */
	gtk_widget_set_sensitive(lookup_widget("menuitemgray"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitemcolor"), any_RGB_image_is_loaded);

#ifdef HAVE_LIBCURL
	/* Updates check */
	gtk_widget_set_visible(lookup_widget("help_update"), TRUE);
	/* Astrometry tool */
	gtk_widget_set_visible(lookup_widget("menuitem_IPS"), TRUE);
#endif
}

gboolean redraw(int vport, int doremap) {
	if (com.script) return FALSE;
	GtkWidget *widget;

	if (vport >= MAXVPORT) {
		fprintf(stderr, _("redraw: maximum number of layers supported is %d"
				" (current image has %d).\n"), MAXVPORT, vport);
		return FALSE;
	}
	widget = com.vport[vport];

	switch (vport) {
	case RED_VPORT:
	case BLUE_VPORT:
	case GREEN_VPORT:
		if (doremap == REMAP_ONLY) {
			remap(vport);
		} else if (doremap == REMAP_ALL) {
			stfComputed = FALSE;
			int i;
//#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)		//probably causes crashes in HESTEQ_MODE
			for (i = 0; i < gfit.naxes[2]; i++) {
				remap(i);
			}
		}
		gtk_widget_queue_draw(widget);
		if (gfit.naxes[2] == 1)
			break;
		/* no break */
	case RGB_VPORT:
		if (gfit.naxis == 3) {
			if (doremap != REMAP_NONE) {
				remaprgb();
			}
			widget = com.vport[RGB_VPORT];
			gtk_widget_queue_draw(widget);
		}
		break;
	default:
		fprintf(stderr, "redraw: unknown viewport number %d\n", vport);
		break;
	}
	//fprintf(stdout, "end of redraw\n");
	return FALSE;
}

static gboolean redraw_idle(gpointer p) {
	redraw(com.cvport, GPOINTER_TO_INT(p)); // draw stars
	return FALSE;
}

void queue_redraw(int doremap) {	// request a redraw from another thread
	siril_add_idle(redraw_idle, GINT_TO_POINTER(doremap));
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
	if (single_image_is_loaded() &&
			com.cvport < com.uniq->nb_layers && com.uniq->layers &&
			com.seq.current != RESULT_IMAGE) {
		layers = com.uniq->layers;
		nb_layers = com.uniq->nb_layers;
	} else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers
			&& com.seq.layers) {
		layers = com.seq.layers;
		nb_layers = com.seq.nb_layers;
	} else
		return 0;

	if (from_GUI) {
		int raw_mode = gtk_combo_box_get_active(modecombo);
		/* update values in the layer_info for cvport */
		layers[com.cvport].rendering_mode =
				raw_mode >= 0 ? raw_mode : NORMAL_DISPLAY;
		layers[com.cvport].lo = round_to_WORD(gtk_range_get_value(range_lo));
		layers[com.cvport].hi = round_to_WORD(gtk_range_get_value(range_hi));
		layers[com.cvport].cut_over = gtk_toggle_button_get_active(cutmax);
	}
	if (!is_chained)
		return 0;
	mode = layers[com.cvport].rendering_mode;
	lo = layers[com.cvport].lo;
	hi = layers[com.cvport].hi;
	cut_over = layers[com.cvport].cut_over;

	for (i = 0; i < nb_layers; i++) {
		if (i == com.cvport)
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

void update_libraw_interface() {
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
	char str[64], *filename;
	if (com.uniq) {	// unique image
		filename = com.uniq->filename;
		nb_layers = com.uniq->nb_layers;
	} else {	// sequence
		filename = malloc(256);
		seq_get_image_filename(&com.seq, com.seq.current, filename);
		nb_layers = com.seq.nb_layers;
	}
	fn_label = GTK_LABEL(gtk_builder_get_object(builder, "labelfilename_red"));
	g_snprintf(str, sizeof(str), _("%s (channel 0)"), filename);
	gtk_label_set_text(fn_label, str);
	if (nb_layers == 3) {	//take in charge both sequence and single image
		fn_label = GTK_LABEL(
				gtk_builder_get_object(builder, "labelfilename_green"));
		g_snprintf(str, sizeof(str), _("%s (channel 1)"), filename);
		gtk_label_set_text(fn_label, str);
		fn_label = GTK_LABEL(
				gtk_builder_get_object(builder, "labelfilename_blue"));
		g_snprintf(str, sizeof(str), _("%s (channel 2)"), filename);
		gtk_label_set_text(fn_label, str);
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

void display_image_number(int index) {
	static GtkSpinButton *spin = NULL;
	if (!spin)
		spin = GTK_SPIN_BUTTON(lookup_widget("imagenumber_spin"));
	char text[16];
	char format[10];
	if (com.seq.fixed <= 1)
		g_snprintf(format, sizeof(format), "%%d");
	else
		g_snprintf(format, sizeof(format), "%%.%dd", com.seq.fixed);
	g_snprintf(text, sizeof(text), format, com.seq.imgparam[index].filenum);
	gtk_entry_set_text(GTK_ENTRY(spin), text);
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

	gtk_widget_show_all(lookup_widget("data_dialog"));
}

void show_main_gray_window() {
	GtkCheckMenuItem *graycheck = GTK_CHECK_MENU_ITEM(lookup_widget("menuitemgray"));
	GtkWidget *win;

	win = lookup_widget("main_window");
	int x = com.main_w_pos.x;
	int y = com.main_w_pos.y;
	int w = com.main_w_pos.w;
	int h = com.main_w_pos.h;

	gtk_check_menu_item_set_active(graycheck, TRUE);
	if (com.remember_windows && w >= 0 && h >= 0) {
		gtk_window_move(GTK_WINDOW(win), x, y);
		gtk_window_resize(GTK_WINDOW(win), w, h);
		gtk_window_present_with_time(GTK_WINDOW(lookup_widget("control_window")), GDK_CURRENT_TIME);
	}
	gtk_widget_show(win);
}

void show_rgb_window() {
	GtkCheckMenuItem *rgbcheck = GTK_CHECK_MENU_ITEM(lookup_widget("menuitemcolor"));
	GtkWidget *win;

	win = lookup_widget("rgb_window");
	int x = com.rgb_w_pos.x;
	int y = com.rgb_w_pos.y;
	int w = com.rgb_w_pos.w;
	int h = com.rgb_w_pos.h;

	gtk_check_menu_item_set_active(rgbcheck, TRUE);
	if (com.remember_windows && w >= 0 && h >= 0) {
		gtk_window_move(GTK_WINDOW(win), x, y);
		gtk_window_resize(GTK_WINDOW(win), w, h);
		gtk_window_present_with_time(GTK_WINDOW(lookup_widget("control_window")), GDK_CURRENT_TIME);
	}
	gtk_widget_show(win);
}

void hide_rgb_window() {
	gtk_widget_hide(lookup_widget("rgb_window"));
}

void hide_gray_window() {
	gtk_widget_hide(lookup_widget("main_window"));
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

static void zoomcombo_update_display_for_zoom() {
	static GtkComboBox *zoomcombo = NULL;
	static double indexes[] = { 16., 8., 4., 2., 1., .5, .25, .125, /*.0625, */
	-1. };
	int i;
	char *msg;

	if (zoomcombo == NULL)
		zoomcombo = GTK_COMBO_BOX(lookup_widget("combozoom"));
	for (i = 0; i < sizeof(indexes) / sizeof(double); i++) {
		if (indexes[i] == com.zoom_value) {
			g_signal_handlers_block_by_func(zoomcombo, on_combozoom_changed,
					NULL);
			gtk_combo_box_set_active(zoomcombo, i);
			g_signal_handlers_unblock_by_func(zoomcombo, on_combozoom_changed,
					NULL);
			return;
		}
	}
	msg = siril_log_message(
			_("Unknown zoom_value value, what is the current zoom?\n"));
	siril_message_dialog( GTK_MESSAGE_ERROR, _("Error"), msg);
}

static void initialize_FITS_name_entries() {
	GtkEntry *moffset, *mdark, *mflat, *final_stack;
	GString *str[4];
	gchar *txt[4];
	gint i;

	moffset = GTK_ENTRY(lookup_widget("offsetname_entry"));
	mdark = GTK_ENTRY(lookup_widget("darkname_entry"));
	mflat = GTK_ENTRY(lookup_widget("flatname_entry"));
	final_stack = GTK_ENTRY(lookup_widget("entryresultfile"));

	str[0] = g_string_new("master-offset");
	str[1] = g_string_new("master-dark");
	str[2] = g_string_new("master-flat");
	str[3] = g_string_new("stack_result");

	for (i = 0; i < 4; i++) {
		str[i] = g_string_append(str[i], com.ext);
		txt[i] = g_string_free(str[i], FALSE);
	}
	gtk_entry_set_text(moffset, txt[0]);
	gtk_entry_set_text(mdark, txt[1]);
	gtk_entry_set_text(mflat, txt[2]);
	gtk_entry_set_text(final_stack, txt[3]);

	for (i = 0; i < 4; i++)
		g_free(txt[i]);
}

void adjust_vport_size_to_image() {
	if (com.script) return;
	int vport;
	// make GtkDrawingArea the same size than the image
	// http://developer.gnome.org/gtk3/3.4/GtkWidget.html#gtk-widget-set-size-request
	double zoom = get_zoom_val();
	int w, h;
	if (zoom <= 0)
		return;
	w = (int) (((double) gfit.rx) * zoom);
	h = (int) (((double) gfit.ry) * zoom);
	for (vport = 0; vport < MAXVPORT; vport++)
		gtk_widget_set_size_request(com.vport[vport], w, h);
	siril_debug_print("set new vport size (%d, %d)\n", w, h);
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
	static GtkWidget *ref_butt = NULL;
	if (ref_butt == NULL)
		ref_butt = lookup_widget("refframe");

	g_signal_handlers_block_by_func(ref_butt, on_ref_frame_toggled, NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ref_butt), com.seq.reference_image == n);
	g_signal_handlers_unblock_by_func(ref_butt, on_ref_frame_toggled, NULL);
}

void close_tab() {
	GtkNotebook* Color_Layers = GTK_NOTEBOOK(lookup_widget("notebook1"));
	GtkWidget* page;

	if (com.seq.nb_layers == 1 || gfit.naxes[2] == 1) {
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
static void add_accelerator_to_tooltip(GtkWidget *widget, guint key, GdkModifierType mod) {
	gchar *text = gtk_widget_get_tooltip_text(widget);
	gchar *accel_str = gtk_accelerator_get_label(key, mod);
	gchar *tip = g_strdup_printf("%s (%s)", text, accel_str);

	gtk_widget_set_tooltip_text(widget, tip);
	g_free(accel_str);
	g_free(tip);
	g_free(text);
}

static void initialize_shortcuts() {
	/* activate accelerators (keyboard shortcut in GTK language) */
	static GtkAccelGroup *accel = NULL;
	GdkModifierType mod = get_default_modifier();

	if (accel == NULL) {
		accel = GTK_ACCEL_GROUP(gtk_builder_get_object(builder, "accelgroup1"));
	}
	/* EXIT */
	gtk_widget_add_accelerator(lookup_widget("exit"), "activate", accel,
	GDK_KEY_q, mod, GTK_ACCEL_VISIBLE);
	/* UNDO */
	gtk_widget_add_accelerator(lookup_widget("undo_item"), "activate", accel,
	GDK_KEY_z, mod, GTK_ACCEL_VISIBLE);
	gtk_widget_add_accelerator(lookup_widget("undo_item1"), "activate", accel,
	GDK_KEY_z, mod, GTK_ACCEL_VISIBLE);
	/* REDO */
#ifdef _WIN32
	gtk_widget_add_accelerator(lookup_widget("redo_item"), "activate", accel,
	GDK_KEY_y, mod, GTK_ACCEL_VISIBLE);
	gtk_widget_add_accelerator(lookup_widget("redo_item1"), "activate", accel,
	GDK_KEY_y, mod, GTK_ACCEL_VISIBLE);
#else
	gtk_widget_add_accelerator(lookup_widget("redo_item"), "activate", accel,
	GDK_KEY_z, mod | GDK_SHIFT_MASK, GTK_ACCEL_VISIBLE);
	gtk_widget_add_accelerator(lookup_widget("redo_item1"), "activate", accel,
	GDK_KEY_z, mod | GDK_SHIFT_MASK, GTK_ACCEL_VISIBLE);
#endif
	/* PREFERENCES */
	gtk_widget_add_accelerator(lookup_widget("settings"), "activate", accel,
	GDK_KEY_p, mod, GTK_ACCEL_VISIBLE);
	/* OPEN */
	gtk_widget_add_accelerator(lookup_widget("open1"), "activate", accel,
	GDK_KEY_o, mod, GTK_ACCEL_VISIBLE);
	/* SAVE */
	gtk_widget_add_accelerator(lookup_widget("save1"), "activate", accel,
	GDK_KEY_s, mod, GTK_ACCEL_VISIBLE);
	/* NEGATIVE */
	gtk_widget_add_accelerator(lookup_widget("menu_negative"), "activate", accel,
	GDK_KEY_i, mod, GTK_ACCEL_VISIBLE);
	/* OPEN WD */
	gtk_widget_add_accelerator(lookup_widget("cwd_button"), "clicked", accel,
	GDK_KEY_d, mod, GTK_ACCEL_VISIBLE);
	add_accelerator_to_tooltip(lookup_widget("cwd_button"), GDK_KEY_d, mod);
}

void initialize_remap() {
	int i;
	for (i = 0; i < MAXGRAYVPORT; i++) {
		remap_index[i] = NULL;
		last_pente[i] = 0.f;
		last_mode[i] = HISTEQ_DISPLAY;
		// only HISTEQ mode always computes the index, it's a good initializer here
	}
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

	gtk_label_set_text(label, com.wd);

	str = g_strdup_printf("%s v%s - %s", PACKAGE, VERSION, com.wd);
	gtk_window_set_title(GTK_WINDOW(lookup_widget("main_window")), str);
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

void set_libraw_settings_menu_available(gboolean activate) {
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
	gtk_combo_box_set_active(pattern, com.debayer.bayer_pattern);
	gtk_combo_box_set_active(inter, com.debayer.bayer_inter);
	gtk_toggle_button_set_active(compat, com.debayer.compatibility);
	gtk_toggle_button_set_active(use_header, com.debayer.use_bayer_header);
	gtk_toggle_button_set_active(demosaicingButton,	com.debayer.open_debayer);
	gtk_toggle_button_set_active(stretch_cfa, com.debayer.stretch);

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
  {"text/uri-list", 0, 0}
};

void initialize_all_GUI(gchar *supported_files) {
	/* initializing internal structures with widgets (drawing areas) */
	com.vport[RED_VPORT] = lookup_widget("drawingarear");
	com.vport[GREEN_VPORT] = lookup_widget("drawingareag");
	com.vport[BLUE_VPORT] = lookup_widget("drawingareab");
	com.vport[RGB_VPORT] = lookup_widget("drawingareargb");
	com.preview_area[0] = lookup_widget("drawingarea_preview1");
	com.preview_area[1] = lookup_widget("drawingarea_preview2");
	initialize_remap();
	initialize_scrollbars();
	init_mouse();

	/* Keybord Shortcuts */
	initialize_shortcuts();

	/* Select combo boxes that trigger some text display or other things */
	gtk_combo_box_set_active(GTK_COMBO_BOX(gtk_builder_get_object(builder, "comboboxstack_methods")), 0);
	zoomcombo_update_display_for_zoom();

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
	g_object_ref(G_OBJECT(lookup_widget("main_window"))); // don't destroy it on removal
	g_object_ref(G_OBJECT(lookup_widget("rgb_window")));  // don't destroy it on removal

	update_used_memory();
}

/*****************************************************************************
 *      P U B L I C      C A L L B A C K      F U N C T I O N S              *
 ****************************************************************************/

void on_register_all_toggle(GtkToggleButton *togglebutton, gpointer user_data) {
	update_reg_interface(TRUE);
}

/* callback for GtkDrawingArea, draw event
 * see http://developer.gnome.org/gtk3/3.2/GtkDrawingArea.html
 * http://developer.gnome.org/gdk-pixbuf/stable/gdk-pixbuf-Image-Data-in-Memory.html
 * http://www.cairographics.org/manual/
 * http://www.cairographics.org/manual/cairo-Image-Surfaces.html#cairo-image-surface-create-for-data
 */
gboolean redraw_drawingarea(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int image_width, image_height, window_width, window_height;
	int vport, i = 0;
	double zoom;

	// we need to identify which vport is being redrawn
	vport = match_drawing_area_widget(widget, TRUE);
	if (vport == -1) {
		fprintf(stderr, "Could not find the vport for the draw callback\n");
		return TRUE;
	}

	window_width = gtk_widget_get_allocated_width(widget);
	window_height = gtk_widget_get_allocated_height(widget);
	zoom = get_zoom_val();
	image_width = (int) (((double) window_width) / zoom);
	image_height = (int) (((double) window_height) / zoom);

	/* draw the RGB and gray images */
	if (vport == RGB_VPORT) {
		if (com.rgbbuf) {
			cairo_scale(cr, zoom, zoom);
			cairo_set_source_surface(cr, com.surface[RGB_VPORT], 0, 0);
			cairo_paint(cr);
		} else {
			fprintf(stdout, "RGB buffer is empty, drawing black image\n");
			draw_empty_image(cr, window_width, window_height);
		}
	} else {
		if (com.graybuf[vport]) {
			cairo_scale(cr, zoom, zoom);
			cairo_set_source_surface(cr, com.surface[vport], 0, 0);
			cairo_paint(cr);
		} else {
			fprintf(stdout, "Buffer %d is empty, drawing black image\n", vport);
			draw_empty_image(cr, window_width, window_height);
		}
	}

	/* draw the selection rectangle */
	if (com.selection.w > 0 && com.selection.h > 0) {
		static double dash_format[] = { 4.0, 2.0 };
		/* fprintf(stdout, "drawing the selection rectangle (%d,%d) (%d,%d)\n",
		 com.selection.x, com.selection.y, com.selection.w, com.selection.h); */
		cairo_set_line_width(cr, 0.8 / zoom);
		cairo_set_dash(cr, dash_format, 2, 0);
		cairo_set_source_rgb(cr, 0.8, 1.0, 0.8);
		cairo_rectangle(cr, (double) com.selection.x, (double) com.selection.y,
				(double) com.selection.w, (double) com.selection.h);
		cairo_stroke(cr);
	}

	/* draw detected stars and highlight the selected star */
	if (com.stars) {
		/* com.stars is a NULL-terminated array */
		cairo_set_dash(cr, NULL, 0, 0);
		cairo_set_source_rgba(cr, 1.0, 0.4, 0.0, 0.9);
		cairo_set_line_width(cr, 1.5 / zoom);

		while (com.stars[i]) {
			// by design Sx>Sy, we redefine FWHM to be sure to have the value in px
			double size = sqrt(com.stars[i]->sx / 2.) * 2 * sqrt(log(2.) * 3);

			if (i == com.selected_star) {
				// We draw horizontal and vertical lines to show the star
				cairo_set_line_width(cr, 2.0 / zoom);
				cairo_set_source_rgba(cr, 0.0, 0.4, 1.0, 0.6);

				cairo_move_to(cr, com.stars[i]->xpos, 0);
				cairo_line_to(cr, com.stars[i]->xpos, image_height);
				cairo_stroke(cr);
				cairo_move_to(cr, 0, com.stars[i]->ypos);
				cairo_line_to(cr, image_width, com.stars[i]->ypos);
				cairo_stroke(cr);

				cairo_set_source_rgba(cr, 1.0, 0.4, 0.0, 0.9);
				cairo_set_line_width(cr, 1.5 / zoom);
			}
			cairo_arc(cr, com.stars[i]->xpos, com.stars[i]->ypos, size, 0., 2. * M_PI);
			cairo_stroke(cr);
			i++;
		}
	}

	if (sequence_is_loaded() && com.seq.current >= 0) {
		/* draw seqpsf stars */
		for (i = 0; i < MAX_SEQPSF && com.seq.photometry[i]; i++) {
			cairo_set_dash(cr, NULL, 0, 0);
			cairo_set_source_rgba(cr, com.seq.photometry_colors[i][0],
					com.seq.photometry_colors[i][1],
					com.seq.photometry_colors[i][2], 1.0);
			cairo_set_line_width(cr, 2.0 / zoom);
			fitted_PSF *the_psf = com.seq.photometry[i][com.seq.current];
			if (the_psf) {
				double size = the_psf->sx;
				cairo_arc(cr, the_psf->xpos, the_psf->ypos, size, 0., 2. * M_PI);
				cairo_stroke(cr);
				cairo_arc(cr, the_psf->xpos, the_psf->ypos, com.phot_set.inner, 0.,
						2. * M_PI);
				cairo_stroke(cr);
				cairo_arc(cr, the_psf->xpos, the_psf->ypos, com.phot_set.outer, 0.,
						2. * M_PI);
				cairo_stroke(cr);
				cairo_select_font_face(cr, "Purisa", CAIRO_FONT_SLANT_NORMAL,
						CAIRO_FONT_WEIGHT_BOLD);
				cairo_set_font_size(cr, 40);
				cairo_move_to(cr, the_psf->xpos + com.phot_set.outer + 5, the_psf->ypos);
				if (i == 0) {
					cairo_show_text(cr, "V");
				}
				else {
					char tmp[2];
					g_snprintf(tmp, 2, "%d", i);
					cairo_show_text(cr, tmp);
				}
				cairo_stroke(cr);
			}
		}

		/* draw a cross on excluded images */
		if (com.seq.imgparam && com.seq.current >= 0 &&
				!com.seq.imgparam[com.seq.current].incl) {
			if (image_width > gfit.rx)
				image_width = gfit.rx;
			if (image_height > gfit.ry)
				image_height = gfit.ry;
			cairo_set_dash(cr, NULL, 0, 0);
			cairo_set_source_rgb(cr, 1.0, 0.8, 0.7);
			cairo_set_line_width(cr, 2.0 / zoom);
			cairo_move_to(cr, 0, 0);
			cairo_line_to(cr, image_width, image_height);
			cairo_move_to(cr, 0, image_height);
			cairo_line_to(cr, image_width, 0);
			cairo_stroke(cr);
		}

		/* draw preview rectangles for the manual registration */
		for (i = 0; i < PREVIEW_NB; i++) {
			if (com.seq.previewX[i] >= 0) {
				int textX, textY;
				char text[3];
				cairo_set_line_width(cr, 0.5 / zoom);
				cairo_set_source_rgb(cr, 0.1, 0.6, 0.0);
				cairo_rectangle(cr,
						com.seq.previewX[i] - com.seq.previewW[i] / 2,
						com.seq.previewY[i] - com.seq.previewH[i] / 2,
						com.seq.previewW[i], com.seq.previewH[i]);
				cairo_stroke(cr);

				textX = com.seq.previewX[i] - com.seq.previewW[i] / 2;
				if (textX < 0)
					textX += com.seq.previewW[i] - 20;
				else
					textX += 15;
				textY = com.seq.previewY[i] - com.seq.previewH[i] / 2;
				if (textY < 0)
					textY += com.seq.previewH[i] - 15;
				else
					textY += 20;
				g_snprintf(text, sizeof(text), "%d", i + 1);

				cairo_set_font_size(cr, 12.0 / zoom);
				cairo_move_to(cr, textX, textY);
				cairo_show_text(cr, text);
			}
		}
	}

	/* draw background removal gradient selection boxes */
	GSList *list;
	for (list = com.grad_samples; list; list = list->next) {
		background_sample *sample = (background_sample *)list->data;
		if (sample->valid) {
			int radius = (int) (sample->size / 2);
			cairo_set_line_width(cr, 1.5 / zoom);
			cairo_set_source_rgba(cr, 0.2, 1.0, 0.3, 1.0);
			cairo_rectangle(cr, sample->position.x - radius, sample->position.y - radius,
					sample->size, sample->size);
			cairo_stroke(cr);
		}
	}
	return FALSE;
}

/* when the cursor moves, update the value displayed in the textbox and save it
 * in the related layer_info. Does not change display until cursor is released. */
void on_minscale_changed(GtkRange *range, gpointer user_data) {
	static GtkEntry *minentry = NULL;
	gchar *buffer;
	if (minentry == NULL)
		minentry = GTK_ENTRY(gtk_builder_get_object(builder, "min_entry"));

	if (single_image_is_loaded() && com.seq.current < RESULT_IMAGE &&
			com.cvport < com.uniq->nb_layers) {
		com.uniq->layers[com.cvport].lo = (int) gtk_range_get_value(range);
		buffer = g_strdup_printf("%u", com.uniq->layers[com.cvport].lo);
	} else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers) {
		com.seq.layers[com.cvport].lo = (int) gtk_range_get_value(range);
		buffer = g_strdup_printf("%u", com.seq.layers[com.cvport].lo);

	} else return;
	g_signal_handlers_block_by_func(minentry, on_min_entry_changed, NULL);
	gtk_entry_set_text(minentry, buffer);
	g_signal_handlers_unblock_by_func(minentry, on_min_entry_changed, NULL);
	g_free(buffer);
}

gboolean on_minscale_release(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	if (com.sliders != USER) {
		com.sliders = USER;
		sliders_mode_set_state(com.sliders);
	}
	if (copy_rendering_settings_when_chained(FALSE))
		redraw(com.cvport, REMAP_ALL);
	else
		redraw(com.cvport, REMAP_ONLY);
	redraw_previews();
	return FALSE;
}

/* when the cursor moves, update the value displayed in the textbox and save it
 * in the related layer_info. Does not change display until cursor is released. */
void on_maxscale_changed(GtkRange *range, gpointer user_data) {
	static GtkEntry *maxentry = NULL;
	gchar *buffer;
	if (maxentry == NULL)
		maxentry = GTK_ENTRY(gtk_builder_get_object(builder, "max_entry"));

	if (single_image_is_loaded() && com.seq.current < RESULT_IMAGE &&
			com.cvport < com.uniq->nb_layers) {
		com.uniq->layers[com.cvport].hi = (int) gtk_range_get_value(range);
		buffer = g_strdup_printf("%u", com.uniq->layers[com.cvport].hi);
	} else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers) {
		com.seq.layers[com.cvport].hi = (int) gtk_range_get_value(range);
		buffer = g_strdup_printf("%u", com.seq.layers[com.cvport].hi);
	} else return;
	g_signal_handlers_block_by_func(maxentry, on_max_entry_changed, NULL);
	gtk_entry_set_text(maxentry, buffer);
	g_signal_handlers_unblock_by_func(maxentry, on_max_entry_changed, NULL);
	g_free(buffer);
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

void on_settings_activate(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("settings_window");
}

void on_menu_FITS_header_activate(GtkMenuItem *menuitem, gpointer user_data) {
	show_FITS_header(&gfit);
}

void on_apply_settings_button_clicked(GtkButton *button, gpointer user_data) {
	update_libraw_interface();
	update_photometry_interface();
	fill_script_paths_list();
	refresh_stars_list(com.stars);
	save_all_windows_position();
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

void on_menu_file_information_activate(GtkMenuItem *menuitem, gpointer user_data) {
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

void save_all_windows_position() {
	if (!com.script && com.remember_windows) {
		GtkWidget *main_w = lookup_widget("main_window");
		GtkWidget *rgb_w = lookup_widget("rgb_window");

		if (gtk_widget_get_visible(main_w))
			com.main_w_pos = get_window_position(GTK_WINDOW(main_w));
		if (gtk_widget_get_visible(rgb_w))
			com.rgb_w_pos = get_window_position(GTK_WINDOW(rgb_w));
	}
}

gboolean on_rgb_window_configure_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {

	save_all_windows_position();
	if (com.zoom_value == -1)
		adjust_vport_size_to_image();
	return FALSE;
}

gboolean on_main_window_configure_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {

	save_all_windows_position();
	if (com.zoom_value == -1)
		adjust_vport_size_to_image();
	return FALSE;
}

void gtk_main_quit() {
	writeinitfile();		// save settings (like window positions)
	close_sequence(FALSE);	// save unfinished business
	close_single_image();	// close the previous image and free resources
	g_slist_free_full(com.script_path, g_free);
	if (!com.headless) {
		g_object_unref(G_OBJECT(lookup_widget("main_window")));
		g_object_unref(G_OBJECT(lookup_widget("rgb_window")));
	}
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

void on_exit_activate(GtkMenuItem *menuitem, gpointer user_data) {
	siril_quit();
}

/* handler for the single-line console */
gboolean on_command_key_press_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {
	const gchar *text;
	int handled = 0;
	static GtkEntry *entry = NULL;
	if (!entry)
		entry = GTK_ENTRY(widget);
	GtkEditable *editable = GTK_EDITABLE(entry);
	int entrylength = 0;

	switch (event->keyval) {
	case GDK_KEY_Return:
	case GDK_KEY_KP_Enter:
		handled = 1;
		text = gtk_entry_get_text(entry);
		history_add_line(strdup(text));
		if (!(processcommand(text)))
			gtk_entry_set_text(entry, "");
		break;
	case GDK_KEY_Up:
		handled = 1;
		if (!com.cmd_history)
			break;
		if (com.cmd_hist_display > 0) {
			if (com.cmd_history[com.cmd_hist_display - 1])
				--com.cmd_hist_display;
			// display previous entry
			gtk_entry_set_text(entry, com.cmd_history[com.cmd_hist_display]);
		} else if (com.cmd_history[com.cmd_hist_size - 1]) {
			// ring back, display previous
			com.cmd_hist_display = com.cmd_hist_size - 1;
			gtk_entry_set_text(entry, com.cmd_history[com.cmd_hist_display]);
		}
		entrylength = gtk_entry_get_text_length(entry);
		gtk_editable_set_position(editable, entrylength);
		break;
	case GDK_KEY_Down:
		handled = 1;
		if (!com.cmd_history)
			break;
		if (com.cmd_hist_display == com.cmd_hist_current)
			break;
		if (com.cmd_hist_display == com.cmd_hist_size - 1) {
			if (com.cmd_hist_current == 0) {
				// ring forward, end
				gtk_entry_set_text(entry, "");
				com.cmd_hist_display++;
			} else if (com.cmd_history[0]) {
				// ring forward, display next
				com.cmd_hist_display = 0;
				gtk_entry_set_text(entry, com.cmd_history[0]);
			}
		} else {
			if (com.cmd_hist_display == com.cmd_hist_current - 1) {
				// end
				gtk_entry_set_text(entry, "");
				com.cmd_hist_display++;
			} else if (com.cmd_history[com.cmd_hist_display + 1]) {
				// display next
				gtk_entry_set_text(entry,
						com.cmd_history[++com.cmd_hist_display]);
			}
		}
		entrylength = gtk_entry_get_text_length(entry);
		gtk_editable_set_position(editable, entrylength);
		break;
	case GDK_KEY_Page_Up:
	case GDK_KEY_Page_Down:
		handled = 1;
		// go to first and last in history
		break;
	}
	return (handled == 1);
}

/* mouse callbacks */
double marge_size = 10;

static gboolean is_over_the_left_side_of_sel(double zoomedX, double zoomedY, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = marge_size / zoom;
	if ((zoomedX > com.selection.x - s && zoomedX < com.selection.x + s)) {
		if (zoomedY > com.selection.y - s
				&& zoomedY < com.selection.y + com.selection.h + s)
			return TRUE;
	}

	return FALSE;
}

static gboolean is_over_the_right_side_of_sel(double zoomedX, double zoomedY, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = marge_size / zoom;
	if ((zoomedX > com.selection.x + com.selection.w - s
			&& zoomedX < com.selection.x + com.selection.w + s)) {
		if (zoomedY > com.selection.y - s
				&& zoomedY < com.selection.y + com.selection.h + s)
			return TRUE;
	}
	return FALSE;
}

static gboolean is_over_the_bottom_of_sel(double zoomedX, double zoomedY, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = marge_size / zoom;
	if ((zoomedY > com.selection.y + com.selection.h - s
			&& zoomedY < com.selection.y + com.selection.h + s)) {
		if (zoomedX > com.selection.x - s
				&& zoomedX < com.selection.x + com.selection.w + s)
			return TRUE;
	}
	return FALSE;
}

static gboolean is_over_the_top_of_sel(double zoomedX, double zoomedY, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = marge_size / zoom;
	if ((zoomedY > com.selection.y - s && zoomedY < com.selection.y + s)) {
		if (zoomedX > com.selection.x - s
				&& zoomedX < com.selection.x + com.selection.w + s)
			return TRUE;
	}
	return FALSE;
}

gboolean on_drawingarea_button_press_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	if (inimage((GdkEvent *) event)) {
		/* click on RGB image */
		if (widget == com.vport[RGB_VPORT]) {
			if (event->button == 3) {
				do_popup_rgbmenu(widget, event);
				return TRUE;
			}
			return FALSE;
		}

		/* else, click on gray image */
		if (event->button == 1) {	// left click
			if (mouse_status == MOUSE_ACTION_SELECT_REG_AREA) {
				if (com.drawing) {
					com.drawing = FALSE;
				} else {
					double zoom = get_zoom_val();
					if (is_over_the_left_side_of_sel(event->x / zoom,
							event->y / zoom, zoom)) {
						com.drawing = TRUE;
						com.startX = com.selection.x + com.selection.w;
						com.freezeY = TRUE;
						com.freezeX = FALSE;
					} else if (is_over_the_right_side_of_sel(event->x / zoom,
							event->y / zoom, zoom)) {
						com.drawing = TRUE;
						com.startX = com.selection.x;
						com.freezeY = TRUE;
						com.freezeX = FALSE;
					} else if (is_over_the_bottom_of_sel(event->x / zoom,
							event->y / zoom, zoom)) {
						com.drawing = TRUE;
						com.startY = com.selection.y;
						com.freezeY = FALSE;
						com.freezeX = TRUE;
					} else if (is_over_the_top_of_sel(event->x / zoom,
							event->y / zoom, zoom)) {
						com.drawing = TRUE;
						com.startY = com.selection.y + com.selection.h;
						com.freezeY = FALSE;
						com.freezeX = TRUE;
					} else {
						com.drawing = TRUE;
						com.startX = event->x / zoom;
						com.startY = event->y / zoom;
						com.selection.h = 0;
						com.selection.w = 0;
						com.freezeX = com.freezeY = FALSE;
					}
				}
				gtk_widget_queue_draw(widget);
			} else if (mouse_status == MOUSE_ACTION_DRAW_SAMPLES) {
				double zoom = get_zoom_val();
				point pt;
				int radius = get_sample_radius();

				pt.x = (event->x / zoom) - radius;
				pt.y = (event->y / zoom) - radius;

				if (pt.x + radius <= gfit.rx && pt.y + radius <= gfit.ry
						&& pt.x - radius >= 0 && pt.y - radius >= 0) {
					com.grad_samples = add_background_sample(com.grad_samples, &gfit, pt);

					redraw(com.cvport, REMAP_NONE);
					redraw_previews();
				}
			}
		} else if (event->button == 3) {	// right click
			if (mouse_status == MOUSE_ACTION_DRAW_SAMPLES) {
				double zoom = get_zoom_val();
				point pt;
				int radius = (int) (25 / 2);

				pt.x = (event->x / zoom) - radius;
				pt.y = (event->y / zoom) - radius;

				if (pt.x + radius <= gfit.rx && pt.y + radius <= gfit.ry
						&& pt.x - radius >= 0 && pt.y - radius >= 0) {
					com.grad_samples = remove_background_sample(com.grad_samples, &gfit, pt);

					redraw(com.cvport, REMAP_NONE);
					redraw_previews();
				}
			}
		}
	}
	return FALSE;
}

gboolean on_drawingarea_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	double zoom = get_zoom_val();
	gdouble zoomedX, zoomedY;

	if (inimage((GdkEvent *) event)) {
		zoomedX = event->x / zoom;
		zoomedY = event->y / zoom;
	} else {
		if (event->x < 0)
			zoomedX = 0.0;
		else if (event->x > gfit.rx * zoom)
			zoomedX = gfit.rx;
		else
			zoomedX = event->x / zoom;
		if (event->y < 0)
			zoomedY = 0.0;
		else if (event->y > gfit.ry * zoom)
			zoomedY = gfit.ry;
		else
			zoomedY = event->y / zoom;
	}
	if (event->button == 1) {
		if (com.drawing && mouse_status == MOUSE_ACTION_SELECT_REG_AREA) {
			com.drawing = FALSE;
			/* finalize selection rectangle coordinates */
			if (!com.freezeX) {
				if (zoomedX > com.startX) {
					com.selection.x = com.startX;
					com.selection.w = zoomedX - com.selection.x;
				} else {
					com.selection.x = zoomedX;
					com.selection.w = com.startX - zoomedX;
				}
			}
			if (!com.freezeY) {
				if (zoomedY > com.startY) {
					com.selection.y = com.startY;
					com.selection.h = zoomedY - com.selection.y;
				} else {
					com.selection.y = zoomedY;
					com.selection.h = com.startY - zoomedY;
				}
			}
			/* we have a new rectangular selection zone,
			 * or an unselection (empty zone) */
			new_selection_zone();

			/* calculate and display FWHM - not in event
			 * callbacks because it's in the same file and
			 * requires a special argument */
			calculate_fwhm(widget);
		} else if (mouse_status == MOUSE_ACTION_SELECT_PREVIEW1) {
			set_preview_area(0, zoomedX, zoomedY);
			mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
			// redraw to get the position of the new preview area
			gtk_widget_queue_draw(widget);
		} else if (mouse_status == MOUSE_ACTION_SELECT_PREVIEW2) {
			set_preview_area(1, zoomedX, zoomedY);
			mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
			gtk_widget_queue_draw(widget);
		}
		is_shift_on = FALSE;
	} else if (event->button == 2) {
		if (inimage((GdkEvent *) event)) {
			double dX, dY, w, h;

			dX = 1.5 * com.phot_set.outer;
			dY = dX;
			w = 3 * com.phot_set.outer;
			h = w;

			if ((dX <= zoomedX) && (dY <= zoomedY)
					&& (zoomedX - dX + w < gfit.rx)
					&& (zoomedY - dY + h < gfit.ry)) {

				com.selection.x = zoomedX - dX;
				com.selection.y = zoomedY - dY;
				com.selection.w = w;
				com.selection.h = h;

				new_selection_zone();
				calculate_fwhm(widget);
			}
		}

	} else if (event->button == 3) {
		if (mouse_status != MOUSE_ACTION_DRAW_SAMPLES) {
			do_popup_graymenu(widget, NULL);
		}
	}
	return FALSE;
}

gboolean on_drawingarea_motion_notify_event(GtkWidget *widget,
		GdkEventMotion *event, gpointer user_data) {
	//static int delay = 5;
	char label[32] = "labeldensity";
	fits *fit = &(gfit);
	double zoom = get_zoom_val();
	gint zoomedX = 0, zoomedY = 0;

	if (inimage((GdkEvent *) event)) {
		char buffer[45];
		char format[25], *format_base = "x: %%.%dd y: %%.%dd = %%.%dd";
		int coords_width = 3, val_width = 3;
		zoomedX = (gint) (event->x / zoom);
		zoomedY = (gint) (event->y / zoom);
		if (fit->rx >= 1000 || fit->ry >= 1000)
			coords_width = 4;
		if (fit->hi >= 1000)
			val_width = 4;
		if (fit->hi >= 10000)
			val_width = 5;
		g_snprintf(format, sizeof(format), format_base, coords_width,
				coords_width, val_width);
		g_snprintf(buffer, sizeof(buffer), format, zoomedX, zoomedY,
				fit->pdata[com.cvport][fit->rx * (fit->ry - zoomedY - 1)
						+ zoomedX]);
		/* TODO: fix to use the new function vport_number_to_name() */
		if (widget == com.vport[RED_VPORT])
			strcat(label, "r");
		else if (widget == com.vport[GREEN_VPORT])
			strcat(label, "g");
		else if (widget == com.vport[BLUE_VPORT])
			strcat(label, "b");
		else
			return FALSE;
		gtk_label_set_text(GTK_LABEL(gtk_builder_get_object(builder, label)),
				buffer);
	}

	if (com.drawing) {	// with button 1 down
		if (!inimage((GdkEvent *) event)) {
			set_cursor("crosshair");
			if (event->x < 0)
				zoomedX = 0;
			else if (event->x > gfit.rx * zoom)
				zoomedX = gfit.rx;
			else
				zoomedX = round_to_int(event->x / zoom);
			if (event->y < 0)
				zoomedY = 0;
			else if (event->y > gfit.ry * zoom)
				zoomedY = gfit.ry;
			else
				zoomedY = round_to_int(event->y / zoom);
		}
		if (!com.freezeX) {
			if (zoomedX > com.startX) {
				com.selection.x = com.startX;
				com.selection.w = zoomedX - com.selection.x;
			} else {
				com.selection.x = zoomedX;
				com.selection.w = com.startX - zoomedX;
			}
		}
		if (!com.freezeY) {
			if (zoomedY > com.startY) {
				com.selection.y = com.startY;
				if (is_shift_on)
					com.selection.h = com.selection.w;
				else
					com.selection.h = zoomedY - com.selection.y;
			} else {
				com.selection.y = zoomedY;
				if (is_shift_on)
					com.selection.h = com.selection.w;
				else
					com.selection.h = com.startY - zoomedY;
			}
		}
		gtk_widget_queue_draw(widget);
	}
	if (inimage((GdkEvent *) event)) {
		if (mouse_status == MOUSE_ACTION_DRAW_SAMPLES) {
			set_cursor("cell");
		} else {
			if (!com.drawing) {
				if (is_over_the_left_side_of_sel(zoomedX, zoomedY, zoom)) {
					set_cursor("w-resize");
				} else if (is_over_the_right_side_of_sel(zoomedX, zoomedY, zoom)) {
					set_cursor("e-resize");
				} else if (is_over_the_bottom_of_sel(zoomedX, zoomedY, zoom)) {
					set_cursor("s-resize");
				} else if (is_over_the_top_of_sel(zoomedX, zoomedY, zoom)) {
					set_cursor("n-resize");
				} else {
					set_cursor("crosshair");
				}
			}
		}
	}

	return FALSE;
}

void on_drawingarea_leave_notify_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	/* trick to get default cursor */
	set_cursor_waiting(FALSE);
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

gboolean on_main_window_key_press_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {

	GtkWidget *WidgetFocused = gtk_window_get_focus(
			GTK_WINDOW(lookup_widget("main_window")));

	if ((WidgetFocused != lookup_widget("max_entry")
			&& WidgetFocused != lookup_widget("min_entry"))) {
		/*** ZOOM shortcuts ***/
		double oldzoom;
		//fprintf(stdout, "drawing area key event\n");
		oldzoom = com.zoom_value;
		is_shift_on = FALSE;

		switch (event->keyval) {
		case GDK_KEY_plus:
		case GDK_KEY_KP_Add:
			if (oldzoom < 0)
				com.zoom_value = 1.0;
			else
				com.zoom_value = min(ZOOM_MAX, oldzoom * 2.0);
			break;
		case GDK_KEY_minus:
		case GDK_KEY_KP_Subtract:
			if (oldzoom < 0)
				com.zoom_value = 1.0;
			else
				com.zoom_value = max(ZOOM_MIN, oldzoom / 2.0);
			break;
		case GDK_KEY_equal:
		case GDK_KEY_KP_Multiply:
			com.zoom_value = 1.0;
			break;
		case GDK_KEY_KP_0:
		case GDK_KEY_0:
			com.zoom_value = -1.0;
			break;
		case GDK_KEY_Shift_L:
		case GDK_KEY_Shift_R:
			is_shift_on = TRUE;
			break;
		default:
			//~ fprintf(stdout, "No bind found for key '%x'.\n", event->keyval);
			break;
		}
		if (com.zoom_value != oldzoom) {
			fprintf(stdout, _("new zoom value: %f\n"), com.zoom_value);
			zoomcombo_update_display_for_zoom();
			adjust_vport_size_to_image();
			redraw(com.cvport, REMAP_NONE);
		}
	}
	return FALSE;
}

void on_excludebutton_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	if (!sequence_is_loaded() || com.seq.current < 0) {
		return;
	}
	toggle_image_selection(com.seq.current);
}

struct load_img_data {
	int index, do_display;
};

static gboolean seq_load_image_as_idle(gpointer p) {
	struct load_img_data *data = (struct load_img_data *)p;
	if (com.seq.current != data->index) // avoid double click
		seq_load_image(&com.seq, data->index, data->do_display);
	free(data);
	return FALSE;
}

/* Returns TRUE if the value has been displayed */
gboolean on_imagenumber_spin_output(GtkSpinButton *spin, gpointer user_data) {
	static GtkAdjustment *adjustment = NULL;
	int index, do_display;
	if (adjustment == NULL)
		adjustment = gtk_spin_button_get_adjustment(spin);
	index = (int) gtk_adjustment_get_value(adjustment);

	//fprintf(stdout, "spinchanged output: index=%d\n", index);
	if (!sequence_is_loaded()) {
		return FALSE;
	}
	if (index > com.seq.number || com.seq.current == index) {// already done for this one
		//fprintf(stderr, "index already processed or data not init. Current: %d, idx: %d\n",
		//		com.seq.current, index);
		return TRUE;
	}
	//fprintf(stdout, "SPINCHANGED: index=%d\n", index);

	do_display = (com.seq.imgparam[index].incl || com.show_excluded);
	struct load_img_data *data = malloc(sizeof(struct load_img_data));
	data->index = index;
	data->do_display = do_display;
	//return !seq_load_image(&com.seq, index, do_display);
	gdk_threads_add_idle(seq_load_image_as_idle, data);
	return TRUE;
}

/* for the spin button to be able to display number which are not the actual value of
 * the adjustment, the output callback is used to modify the way they are displayed,
 * but the input callback is also required to do the opposite operation, i.e. convert
 * real image number into the number in the sequence, which is the value of the adj. */
gboolean on_imagenumberspin_input(GtkSpinButton *spin, gdouble *new_val,
		gpointer p) {
	const char *imgname = gtk_entry_get_text(GTK_ENTRY(spin));
	int imgname_int = atoi(imgname);
	int i;
	if (!sequence_is_loaded()) {
		//fprintf(stderr, "No sequence loaded\n");
		return FALSE;
	}
	i = com.seq.current;
	//fprintf(stdout, "current index in input is %d (looking for %d)\n", i, imgname_int);
	if (com.seq.imgparam[i].filenum == imgname_int) {
		*new_val = (gdouble) i;
		//fprintf(stdout, "found at 0 %d\n", i);
		return TRUE;
	} else if (i > 0 && com.seq.imgparam[i - 1].filenum == imgname_int) {
		*new_val = (gdouble) (i - 1);
		//fprintf(stdout, "found at -1 %d\n", i-1);
		return TRUE;
	} else if (i < com.seq.number - 1
			&& com.seq.imgparam[i + 1].filenum == imgname_int) {
		*new_val = (gdouble) (i + 1);
		//fprintf(stdout, "found at +1 %d\n", i+1);
		return TRUE;
	} else {
		/* no luck with neighbours, sweep it all */
		for (i = 0; i < com.seq.number; i++) {
			if (com.seq.imgparam[i].filenum == imgname_int) {
				*new_val = (gdouble) i;
				//fprintf(stdout, "sweep found at %d\n", i);
				return TRUE;
			}
		}
	}
	return GTK_INPUT_ERROR;
}

/* Returns:	TRUE to stop other handlers from being invoked for the event.
 *		FALSE to propagate the event further. */
gboolean on_imagenumberspin_key_release_event(GtkWidget *widget,
		GdkEventKey *event, gpointer user_data) {
	static GtkAdjustment *adj = NULL;
	int n;

	if (!sequence_is_loaded()) {
		//fprintf(stderr, "No sequence loaded\n");
		return TRUE;
	}
	if (!adj)
		adj = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(widget));
	n = (int) gtk_adjustment_get_value(adj);
	if (n > com.seq.number) {
		return TRUE;
	}
	if (event->keyval == GDK_KEY_space) {
		toggle_image_selection(n);
	}
	return FALSE;
}

void on_seqexcludeall_button_clicked(GtkButton *button, gpointer user_data) {
	gboolean exclude_all;

	exclude_all = siril_confirm_dialog(_("Exclude all images ?"),
			_("This erases previous image selection and there's no possible undo."), FALSE);
	if (exclude_all) {
		sequence_setselect_all(FALSE);
	}
}

void on_seqselectall_button_clicked(GtkButton *button, gpointer user_data) {
	gboolean select_all;

	select_all = siril_confirm_dialog(_("Include all images ?"),
			_("This erases previous image selection and there's no possible undo."), FALSE);
	if (select_all) {
		sequence_setselect_all(TRUE);
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

void on_showexcluded_button_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	com.show_excluded = gtk_toggle_button_get_active(togglebutton);
}

void on_ref_frame_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	if (sequence_is_loaded()) {
		free_reference_image();
		if ((gtk_toggle_button_get_active(togglebutton) == FALSE)) {
			if (com.seq.reference_image == com.seq.current)
				com.seq.reference_image = -1;
		} else {
			com.seq.reference_image = com.seq.current;
			test_and_allocate_reference_image(-1);
			// a reference image should not be excluded to avoid confusion
			if (!com.seq.imgparam[com.seq.current].incl) {
				toggle_image_selection(com.seq.current);
			}
		}
		sequence_list_change_reference();
		update_stack_interface(FALSE);// get stacking info and enable the Go button
		adjust_sellabel();	// reference image is named in the label
		writeseqfile(&com.seq);
		drawPlot();		// update plots
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
	fill_sequence_list(&com.seq, com.cvport, FALSE);
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

void on_menuitemgray_toggled(GtkCheckMenuItem *checkmenuitem,
		gpointer user_data) {
	if (gtk_check_menu_item_get_active(checkmenuitem))
		gtk_widget_show_all(lookup_widget("main_window"));
	else
		gtk_widget_hide(lookup_widget("main_window"));
}

void on_menuitemcolor_toggled(GtkCheckMenuItem *checkmenuitem,
		gpointer user_data) {
	if (gtk_check_menu_item_get_active(checkmenuitem))
		gtk_widget_show_all(lookup_widget("rgb_window"));
	else
		gtk_widget_hide(lookup_widget("rgb_window"));
}

gboolean rgb_area_popup_menu_handler(GtkWidget *widget) {
	do_popup_rgbmenu(widget, NULL);
	return TRUE;
}

void on_rgb_window_hide(GtkWidget *object, gpointer user_data) {
	GtkCheckMenuItem *rgbcheck = GTK_CHECK_MENU_ITEM(
			gtk_builder_get_object(builder, "menuitemcolor"));
	gtk_check_menu_item_set_active(rgbcheck, FALSE);
}

void on_gray_window_hide(GtkWidget *object, gpointer user_data) {
	GtkCheckMenuItem *graycheck = GTK_CHECK_MENU_ITEM(
			gtk_builder_get_object(builder, "menuitemgray"));
	gtk_check_menu_item_set_active(graycheck, FALSE);
}

void on_combozoom_changed(GtkComboBox *widget, gpointer user_data) {
	gint active = gtk_combo_box_get_active(widget);
	switch (active) {
	case 0: /* 16:1 */
		com.zoom_value = 16.;
		break;
	case 1: /* 8:1 */
		com.zoom_value = 8.;
		break;
	case 2: /* 4:1 */
		com.zoom_value = 4.;
		break;
	case 3: /* 2:1 */
		com.zoom_value = 2.;
		break;
	case -1:
	case 4: /* 1:1 */
		com.zoom_value = 1.;
		break;
	case 5: /* 1:2 */
		com.zoom_value = .5;
		break;
	case 6: /* 1:4 */
		com.zoom_value = .25;
		break;
	case 7: /* 1:8 */
		com.zoom_value = .125;
		break;
	case 8: /* fit to window */
		com.zoom_value = -1.;
		break;
	}
	fprintf(stdout, "zoom is now %f\n", com.zoom_value);
	adjust_vport_size_to_image();
	redraw(com.cvport, REMAP_NONE);
}

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

void on_menuitem_stat_activate(GtkMenuItem *menuitem, gpointer user_data) {
	set_cursor_waiting(TRUE);
	computeStat();
	siril_open_dialog("StatWindow");
	set_cursor_waiting(FALSE);
}

/**********************************************************************/

void on_menu_channel_separation_activate(GtkMenuItem *menuitem,
		gpointer user_data) {
	if (single_image_is_loaded() && isrgb(&gfit))
		siril_open_dialog("extract_channel_dialog");
}

void on_menuitemcalibration_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded() && isrgb(&gfit)) {
		initialize_calibration_interface();
		siril_open_dialog("color_calibration");
	}
}

void on_extract_channel_button_close_clicked(GtkButton *button,
		gpointer user_data) {
	siril_close_dialog("extract_channel_dialog");
}

void on_combo_extract_colors_changed(GtkComboBox *box, gpointer user_data) {
	switch(gtk_combo_box_get_active(box)) {
	default:
	case 0: // RGB
		gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c1")), _("Red: "));
		gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c2")), _("Green: "));
		gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c3")), _("Blue: "));
		break;
	case 1: // HSL
		gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c1")), _("Hue: "));
		gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c2")), _("Saturation: "));
		gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c3")), _("Lightness: "));
		break;
	case 2: // HSV
		gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c1")), _("Hue: "));
		gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c2")), _("Saturation: "));
		gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c3")), _("Value: "));
		break;
	case 3: // CIE L*a*b*
		gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c1")), "L*: ");
		gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c2")), "a*: ");
		gtk_label_set_text(GTK_LABEL(lookup_widget("label_extract_c3")), "b*: ");
	}
}

void on_extract_channel_button_ok_clicked(GtkButton *button, gpointer user_data) {
	static GtkEntry *channel_extract_entry[3] = { NULL, NULL, NULL };
	static GtkComboBox *combo_extract_channel = NULL;

	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return;
	}

	struct extract_channels_data *args = malloc(sizeof(struct extract_channels_data));
	if (args == NULL) {
		PRINT_ALLOC_ERR;
		return;
	}

	if (combo_extract_channel == NULL) {
		combo_extract_channel = GTK_COMBO_BOX(
				lookup_widget("combo_extract_colors"));
		channel_extract_entry[0] = GTK_ENTRY(
				lookup_widget("Ch1_extract_channel_entry"));
		channel_extract_entry[1] = GTK_ENTRY(
				lookup_widget("Ch2_extract_channel_entry"));
		channel_extract_entry[2] = GTK_ENTRY(
				lookup_widget("Ch3_extract_channel_entry"));
	}

	args->type = gtk_combo_box_get_active(combo_extract_channel);
	args->str_type = gtk_combo_box_get_active_id(combo_extract_channel);

	args->channel[0] = gtk_entry_get_text(channel_extract_entry[0]);
	args->channel[1] = gtk_entry_get_text(channel_extract_entry[1]);
	args->channel[2] = gtk_entry_get_text(channel_extract_entry[2]);

	if ((args->channel[0][0] != '\0') && (args->channel[1][0] != '\0')
			&& (args->channel[2][0] != '\0')) {
		args->fit = calloc(1, sizeof(fits));
		set_cursor_waiting(TRUE);
		copyfits(&gfit, args->fit, CP_ALLOC | CP_COPYA | CP_FORMAT, 0);
		start_in_new_thread(extract_channels, args);
	}
	else {
		free(args);
	}
}

/******************* POPUP GRAY MENU *******************************/

void on_menu_gray_psf_activate(GtkMenuItem *menuitem, gpointer user_data) {
	gchar *msg;
	fitted_PSF *result = NULL;
	int layer = match_drawing_area_widget(com.vport[com.cvport], FALSE);
	char *str;

	if (layer == -1)
		return;
	if (!(com.selection.h && com.selection.w))
		return;
	if (com.selection.w > 300 || com.selection.h > 300) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Current selection is too large"),
				_("To determine the PSF, please make a selection around a star."));

		return;
	}
	result = psf_get_minimisation(&gfit, layer, &com.selection, TRUE, TRUE);
	if (!result)
		return;

	if (com.magOffset > 0.0)
		str = "true reduced";
	else
		str = "relative";
	msg = g_strdup_printf(_("Centroid Coordinates:\n\t\tx0=%.2fpx\n\t\ty0=%.2fpx\n\n"
					"Full Width Half Maximum:\n\t\tFWHMx=%.2f%s\n\t\tFWHMy=%.2f%s\n\n"
					"Angle:\n\t\t%0.2fdeg\n\n"
					"Background Value:\n\t\tB=%.6f\n\n"
					"Maximal Intensity:\n\t\tA=%.6f\n\n"
					"Magnitude (%s):\n\t\tm=%.4f\u00B1%.4f\n\n"
					"RMSE:\n\t\tRMSE=%.3e"), result->x0 + com.selection.x,
			com.selection.y + com.selection.h - result->y0, result->fwhmx,
			result->units, result->fwhmy, result->units, result->angle, result->B,
			result->A, str, result->mag + com.magOffset, result->s_mag, result->rmse);
	show_data_dialog(msg, "PSF Results");
	g_free(msg);
	free(result);
}

void on_menu_gray_seqpsf_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (!sequence_is_loaded()) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("PSF for the sequence only applies on sequences"),
					_("Please load a sequence before trying to apply the PSF for the sequence."));
	} else {
		process_seq_psf(0);
	}
}

void on_menu_gray_pick_star_activate(GtkMenuItem *menuitem, gpointer user_data) {
	int layer = match_drawing_area_widget(com.vport[com.cvport], FALSE);
	int new_index;
	GtkCheckMenuItem *PSFcheck = GTK_CHECK_MENU_ITEM(
			gtk_builder_get_object(builder, "menuitemPSF"));
	GtkWidget *window = lookup_widget("stars_list_window");

	if (layer != -1) {
		if (!(com.selection.h && com.selection.w))
			return;
		if (com.selection.w > 300 || com.selection.h > 300) {
			siril_message_dialog(GTK_MESSAGE_WARNING, _("Current selection is too large"),
							_("To determine the PSF, please make a selection around a star."));
			return;
		}
		fitted_PSF *new_star = add_star(&gfit, layer, &new_index);
		if (new_star) {
			add_star_to_list(new_star);
			if (!(gtk_widget_get_visible(window)))//We open the stars_list_window
				gtk_widget_show_all(window);
			gtk_check_menu_item_set_active(PSFcheck, TRUE);
		} else
			return;
	}
	redraw(com.cvport, REMAP_NONE);
}

void on_menu_gray_stat_activate(GtkMenuItem *menuitem, gpointer user_data) {
	computeStat();
	siril_open_dialog("StatWindow");
}


/****************** GUI for Wavelet Layers Extraction *****************/

void on_menu_wavelet_separation_activate(GtkMenuItem *menuitem,
		gpointer user_data) {

	if (single_image_is_loaded()) {
		siril_open_dialog("extract_wavelets_layers_dialog");
	}
}

void on_button_extract_w_ok_clicked(GtkButton *button, gpointer user_data) {
	fits fit = { 0 };
	int Nbr_Plan, Type, maxplan, mins;
	static GtkSpinButton *Spin_Nbr_Plan = NULL;
	static GtkComboBox *Combo_Wavelets_Type = NULL;

	if (Spin_Nbr_Plan == NULL) {
		Spin_Nbr_Plan = GTK_SPIN_BUTTON(lookup_widget("spinbutton_extract_w"));
		Combo_Wavelets_Type = GTK_COMBO_BOX(
				lookup_widget("combo_interpolation_extract_w"));
	}

	Nbr_Plan = gtk_spin_button_get_value(Spin_Nbr_Plan);
	Type = gtk_combo_box_get_active(Combo_Wavelets_Type) + 1;// 1: linear, 2: bspline

	set_cursor_waiting(TRUE);
	mins = min(gfit.rx, gfit.ry);
	maxplan = log(mins) / log(2) - 2;

	if (Nbr_Plan > maxplan) {
		char *msg = siril_log_message(_("Wavelet: maximum number "
				"of plans for this image size is %d\n"), maxplan);
		siril_message_dialog(GTK_MESSAGE_WARNING, _("Warning"), msg);
		set_cursor_waiting(FALSE);
		return;
	}

	copyfits(&gfit, &fit, CP_ALLOC | CP_COPYA | CP_FORMAT, 0);

	extract_plans(&fit, Nbr_Plan, Type);

	clearfits(&fit);
	update_used_memory();
	set_cursor_waiting(FALSE);
}

void on_button_extract_w_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("extract_wavelets_layers_dialog");
}

/******* SPLIT CFA ******************************/

void on_menu_slpitcfa_activate(GtkMenuItem *menuitem, gpointer user_data) {
	siril_open_dialog("split_cfa_dialog");
}

void on_split_cfa_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("split_cfa_dialog");
}

void on_split_cfa_apply_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton *seq = GTK_TOGGLE_BUTTON(lookup_widget("checkSplitCFASeq"));
	GtkEntry *entrySplitCFA;

	entrySplitCFA = GTK_ENTRY(lookup_widget("entrySplitCFA"));

	if (gtk_toggle_button_get_active(seq) && sequence_is_loaded()) {
		struct split_cfa_data *args = malloc(sizeof(struct split_cfa_data));

		set_cursor_waiting(TRUE);
		args->seq = &com.seq;
		args->seqEntry = gtk_entry_get_text(entrySplitCFA);
		if (args->seqEntry && args->seqEntry[0] == '\0')
					args->seqEntry = "CFA_";
		apply_split_cfa_to_sequence(args);
	} else {
		process_split_cfa(0);
	}
}

void on_spinCPU_value_changed (GtkSpinButton *spinbutton, gpointer user_data) {
	com.max_thread = (int)gtk_spin_button_get_value(spinbutton);
}

void on_undo_item_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded() && is_undo_available()) {
		set_cursor_waiting(TRUE);
		undo_display_data(UNDO);
		set_cursor_waiting(FALSE);
	}

	/* update menus */
	update_MenuItem();
}

void on_redo_item_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded() && is_redo_available()) {
		set_cursor_waiting(TRUE);
		undo_display_data(REDO);
		set_cursor_waiting(FALSE);
	}

	/* update menus */
	update_MenuItem();
}

void on_undo_item1_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded() && is_undo_available()) {
		set_cursor_waiting(TRUE);
		undo_display_data(UNDO);
		set_cursor_waiting(FALSE);
	}
}

void on_redo_item1_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded() && is_redo_available()) {
		set_cursor_waiting(TRUE);
		undo_display_data(REDO);
		set_cursor_waiting(FALSE);
	}
}

void on_rememberWindowsCheck_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	com.remember_windows = gtk_toggle_button_get_active(togglebutton);
}

void on_entryAviWidth_changed(GtkEditable *editable, gpointer user_data) {
	double ratio, width, height;
	gchar *c_height;
	GtkEntry *heightEntry = GTK_ENTRY(lookup_widget("entryAviHeight"));

	if (com.selection.w && com.selection.h) return;
	ratio = (double) com.seq.ry / (double) com.seq.rx;
	width = atof(gtk_entry_get_text(GTK_ENTRY(editable)));
	height = ratio * width;
	c_height = g_strdup_printf("%d", (int)(height));

	g_signal_handlers_block_by_func(heightEntry, on_entryAviHeight_changed, NULL);
	gtk_entry_set_text(heightEntry, c_height);
	g_signal_handlers_unblock_by_func(heightEntry, on_entryAviHeight_changed, NULL);
	g_free(c_height);
}

void on_entryAviHeight_changed(GtkEditable *editable, gpointer user_data) {
	double ratio, width, height;
	gchar *c_width;
	GtkEntry *widthEntry = GTK_ENTRY(lookup_widget("entryAviWidth"));

	if (com.selection.w && com.selection.h) return;
	ratio = (double) com.seq.rx / (double) com.seq.ry;
	height = atof(gtk_entry_get_text(GTK_ENTRY(editable)));
	width = ratio * height;
	c_width = g_strdup_printf("%d", (int)(width));

	g_signal_handlers_block_by_func(widthEntry, on_entryAviWidth_changed, NULL);
	gtk_entry_set_text(widthEntry, c_width);
	g_signal_handlers_unblock_by_func(widthEntry, on_entryAviWidth_changed, NULL);
	g_free(c_width);
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
