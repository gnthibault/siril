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

#include <math.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_app_dirs.h"
#include "core/processing.h"
#include "algos/colors.h"
#include "algos/background_extraction.h"
#include "algos/PSF.h"
#include "algos/siril_wcs.h"
#include "algos/annotate.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
#include "histogram.h"
#include "git-version.h"

#include "image_display.h"

typedef struct draw_data {
	cairo_t *cr;
	int vport;
	double zoom;
	cairo_filter_t filter;
	guint image_width, image_height;
	guint window_width, window_height;
} draw_data_t;

/* remap index data, an index for each layer */
static BYTE *remap_index[MAXGRAYVPORT];
static float last_pente[MAXGRAYVPORT];
static display_mode last_mode[MAXGRAYVPORT];

/* STF (auto-stretch) data */
static gboolean stfComputed;	// Flag to know if STF parameters are available
static float stfShadows, stfHighlights, stfM;

static void remaprgb(void) {
	guint32 *dst;
	const guint32 *bufr, *bufg, *bufb;
	gint i;
	int nbdata;

	siril_debug_print("remaprgb\n");
	if (!isrgb(&gfit))
		return;

	// allocate if not already done or the same size
	if (cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, gfit.rx)
			!= com.surface_stride[RGB_VPORT]
			|| gfit.ry != com.surface_height[RGB_VPORT]
			|| !com.surface[RGB_VPORT] || !com.rgbbuf) {
		guchar *oldbuf = com.rgbbuf;
		siril_debug_print("RGB display buffers and surface (re-)allocation\n");
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
			siril_debug_print("Error creating the Cairo image surface for the RGB image\n");
			cairo_surface_destroy(com.surface[RGB_VPORT]);
			com.surface[RGB_VPORT] = NULL;
			return;
		}
	}
	// WARNING : this assumes that R, G and B buffers are already allocated and mapped
	// it seems ok, but one can probably imagine situations where it segfaults
	bufr = (const guint32*) com.graybuf[RED_VPORT];
	bufg = (const guint32*) com.graybuf[GREEN_VPORT];
	bufb = (const guint32*) com.graybuf[BLUE_VPORT];
	if (bufr == NULL || bufg == NULL || bufb == NULL) {
		siril_debug_print("remaprgb: gray buffers not allocated for display\n");
		return;
	}
	dst = (guint32*) com.rgbbuf;	// index is j
	nbdata = gfit.rx * gfit.ry;	// source images are 32-bit RGBA

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread)
#endif
	for (i = 0; i < nbdata; ++i) {
		dst[i] = (bufr[i] & 0xFF0000) | (bufg[i] & 0xFF00) | (bufb[i] & 0xFF);
	}

	// flush to ensure all writing to the image was done and redraw the surface
	cairo_surface_flush(com.surface[RGB_VPORT]);
	cairo_surface_mark_dirty(com.surface[RGB_VPORT]);
}

static int make_index_for_current_display(display_mode mode, WORD lo, WORD hi,
		int vport);

static int make_index_for_rainbow(BYTE index[][3]);

static void remap(int vport) {
	// This function maps fit data with a linear LUT between lo and hi levels
	// to the buffer to be displayed; display only is modified
	guint y;
	BYTE *dst, *index, rainbow_index[UCHAR_MAX + 1][3];
	WORD *src, hi, lo;
	float *fsrc;
	display_mode mode;
	color_map color;
	gboolean do_cut_over, inverted;

	siril_debug_print("remap %d\n", vport);
	if (vport == RGB_VPORT) {
		remaprgb();
		return;
	}

	int no_data = 0;
	if (single_image_is_loaded()) {
		if (vport >= com.uniq->nb_layers)
			no_data = 1;
	} else if (sequence_is_loaded()) {
		if (vport >= com.seq.nb_layers)
			no_data = 1;
	}
	else no_data = 1;
	if (gfit.type == DATA_UNSUPPORTED)
		no_data = 1;
	if (no_data) {
		siril_debug_print("vport is out of bounds or data is not loaded yet\n");
		return;
	}

	// allocate if not already done or the same size
	if (cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, gfit.rx)
			!= com.surface_stride[vport] || gfit.ry != com.surface_height[vport]
			|| !com.surface[vport] || !com.graybuf[vport]) {
		guchar *oldbuf = com.graybuf[vport];
		siril_debug_print("Gray display buffers and surface (re-)allocation\n");
		if (gfit.rx == 0 || gfit.ry == 0) {
			siril_debug_print("gfit has a zero size, must not happen!\n");
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
			siril_debug_print("Error creating the Cairo image surface for vport %d\n",
					vport);
			cairo_surface_destroy(com.surface[vport]);
			com.surface[vport] = NULL;
			return;
		}
	}
	inverted = gtk_toggle_tool_button_get_active(
			GTK_TOGGLE_TOOL_BUTTON(lookup_widget("neg_button")));

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
		siril_debug_print("BUG in unique image remap\n");
		return;
	}

	if (mode == HISTEQ_DISPLAY) {
		if (gfit.type == DATA_USHORT) {
			double hist_sum, nb_pixels;
			size_t i, hist_nb_bins;
			gsl_histogram *histo;

			compute_histo_for_gfit();
			histo = com.layers_hist[vport];
			hist_nb_bins = gsl_histogram_bins(histo);
			nb_pixels = (double)(gfit.rx * gfit.ry);

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
		}
		else if (gfit.type == DATA_FLOAT) {
			double hist_sum, nb_pixels;
			size_t i, hist_nb_bins;
			gsl_histogram *histo;

			compute_histo_for_gfit();
			histo = com.layers_hist[vport];
			hist_nb_bins = gsl_histogram_bins(histo);
			nb_pixels = (double)(gfit.rx * gfit.ry);

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
		}
		set_viewer_mode_widgets_sensitive(FALSE);
	} else {
		// for all other modes and ushort data, the index can be reused
		if (mode == STF_DISPLAY && !stfComputed) {
			stfM = findMidtonesBalance(&gfit, &stfShadows, &stfHighlights);
			stfComputed = TRUE;
		}
		make_index_for_current_display(mode, lo, hi, vport);
		if (mode == STF_DISPLAY)
			set_viewer_mode_widgets_sensitive(FALSE);
		else
			set_viewer_mode_widgets_sensitive(TRUE);
	}

	src = gfit.pdata[vport];
	fsrc = gfit.fpdata[vport];
	dst = com.graybuf[vport];

	color = gtk_toggle_tool_button_get_active(
			GTK_TOGGLE_TOOL_BUTTON(lookup_widget("colormap_button")));

	if (color == RAINBOW_COLOR)
		make_index_for_rainbow(rainbow_index);
	index = remap_index[vport];

	gboolean special_mode = (mode == HISTEQ_DISPLAY || mode == STF_DISPLAY);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(y) schedule(static)
#endif
	for (y = 0; y < gfit.ry; y++) {
		guint x;
		guint src_i = y * gfit.rx;
		guint dst_i = ((gfit.ry - 1 - y) * gfit.rx) * 4;
		for (x = 0; x < gfit.rx; ++x, ++src_i, dst_i += 2) {
			guint src_index = y * gfit.rx + x;
			BYTE dst_pixel_value = 0;
			if (gfit.type == DATA_USHORT) {
				if (special_mode) // special case, no lo & hi
					dst_pixel_value = index[src[src_index]];
				else if (do_cut_over && src[src_index] > hi)	// cut
					dst_pixel_value = 0;
				else {
					dst_pixel_value = index[src[src_index] - lo < 0 ? 0 : src[src_index] - lo];
				}
			} else if (gfit.type == DATA_FLOAT) {
				if (special_mode) // special case, no lo & hi
					dst_pixel_value = index[roundf_to_WORD(fsrc[src_index] * USHRT_MAX_SINGLE)];
				else if (do_cut_over && roundf_to_WORD(fsrc[src_index] * USHRT_MAX_SINGLE) > hi)	// cut
					dst_pixel_value = 0;
				else {
					dst_pixel_value = index[
						roundf_to_WORD(fsrc[src_index] * USHRT_MAX_SINGLE) - lo < 0 ? 0 :
							roundf_to_WORD(fsrc[src_index] * USHRT_MAX_SINGLE) - lo];
				}
			}

			dst_pixel_value = inverted ? UCHAR_MAX - dst_pixel_value : dst_pixel_value;

			// Siril's FITS are stored bottom to top, so mapping needs to revert data order
			guint dst_index = ((gfit.ry - 1 - y) * gfit.rx + x) * 4;
			switch (color) {
				default:
				case NORMAL_COLOR:
					*(guint32*)(dst + dst_index) = dst_pixel_value << 16 | dst_pixel_value << 8 | dst_pixel_value;
					break;
				case RAINBOW_COLOR:
					*(guint32*)(dst + dst_index) = rainbow_index[dst_pixel_value][0] << 16 | rainbow_index[dst_pixel_value][1] << 8 | rainbow_index[dst_pixel_value][2];
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
	float slope;
	int i;
	BYTE *index;
	float pxl;

	/* initialization of data required to build the remap_index */
	switch (mode) {
		case NORMAL_DISPLAY:
			slope = UCHAR_MAX_SINGLE / (float) (hi - lo);
			break;
		case LOG_DISPLAY:
			slope = fabsf(UCHAR_MAX_SINGLE / logf((float)(hi - lo) * 0.1f));
			break;
		case SQRT_DISPLAY:
			slope = UCHAR_MAX_SINGLE / sqrtf((float)(hi - lo));
			break;
		case SQUARED_DISPLAY:
			slope = UCHAR_MAX_SINGLE / SQR((float)(hi - lo));
			break;
		case ASINH_DISPLAY:
			slope = UCHAR_MAX_SINGLE / asinhf((float)(hi - lo) * 0.001f);
			break;
		case STF_DISPLAY:
			slope = UCHAR_MAX_SINGLE;
			break;
		default:
			return 1;
	}
	if ((mode != HISTEQ_DISPLAY && mode != STF_DISPLAY) && slope == last_pente[vport]
			&& mode == last_mode[vport]) {
		siril_debug_print("Re-using previous remap_index\n");
		return 0;
	}
	siril_debug_print("Rebuilding remap_index\n");

	/************* Building the remap_index **************/
	if (!remap_index[vport]) {
		remap_index[vport] = malloc(USHRT_MAX + 1);
		if (!remap_index[vport]) {
			PRINT_ALLOC_ERR;
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
					index[i] = round_to_BYTE(logf((float) i / 10.f) * slope); //10.f is arbitrary: good matching with ds9
				break;
			case SQRT_DISPLAY:
				// sqrt(2^16) = 2^8
				index[i] = round_to_BYTE(sqrtf((float) i) * slope);
				break;
			case SQUARED_DISPLAY:
				// pow(2^4,2) = 2^8
				index[i] = round_to_BYTE(SQR((float)i) * slope);
				break;
			case ASINH_DISPLAY:
				// asinh(2.78*10^110) = 255
				index[i] = round_to_BYTE(asinhf((float) i / 1000.f) * slope); //1000.f is arbitrary: good matching with ds9, could be asinhf(a*Q*i)/Q
				break;
			case NORMAL_DISPLAY:
				index[i] = round_to_BYTE((float) i * slope);
				break;
			case STF_DISPLAY:
				pxl = (gfit.orig_bitpix == BYTE_IMG ?
						(float) i / UCHAR_MAX_SINGLE :
						(float) i / USHRT_MAX_SINGLE);
				index[i] = round_to_BYTE((MTF(pxl, stfM, stfShadows, stfHighlights)) * slope);
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

	last_pente[vport] = slope;
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

static void draw_empty_image(const draw_data_t* dd) {
	cairo_t *cr = dd->cr;
	guint width = dd->window_width;
	guint height = dd->window_height;

	cairo_rectangle(cr, 0, 0, width, height);
	cairo_fill(cr);

	int pix_w = gdk_pixbuf_get_width(com.siril_pix);
	int pix_h = gdk_pixbuf_get_height(com.siril_pix);

	gdk_cairo_set_source_pixbuf(cr, com.siril_pix, (width - pix_w) / 2, (height - pix_h) / 2);
	cairo_paint(cr);
	cairo_fill(cr);

#ifdef SIRIL_UNSTABLE
	{
		GtkWidget *widget = lookup_widget("drawingareargb");
		GtkStyleContext *context = gtk_widget_get_style_context(widget);
		GtkStateFlags state = gtk_widget_get_state_flags(widget);
		PangoLayout *layout;
		gchar *msg;
		GtkAllocation allocation;
		gdouble scale;
		GdkRGBA color;
		gint w, h;

		layout = gtk_widget_create_pango_layout(widget, NULL);

		msg = g_strdup_printf(_("<big>Unstable Development Version</big>\n\n"
				"<small>commit <tt>%s</tt></small>\n\n"
				"<small>Please test bugs against "
				"latest git master branch\n"
				"before reporting them.</small>"),
				SIRIL_GIT_VERSION_ABBREV);
		pango_layout_set_markup(layout, msg, -1);
		g_free(msg);
		pango_layout_set_alignment(layout, PANGO_ALIGN_CENTER);

		pango_layout_get_pixel_size(layout, &w, &h);
		gtk_widget_get_allocation(widget, &allocation);

		scale = MIN(((gdouble ) allocation.width / 2.0) / (gdouble ) w,
				((gdouble ) allocation.height / 2.0) / (gdouble ) h / 2);

		gtk_style_context_get_color(context, state, &color);
		gdk_cairo_set_source_rgba(cr, &color);

		cairo_move_to(cr, (allocation.width - (w * scale)) / 2,
				(allocation.height - (h * scale)) / 2 - pix_h);

		cairo_scale(cr, scale, scale);

		pango_cairo_show_layout(cr, layout);

		g_object_unref(layout);
	}
#endif /* SIRIL_UNSTABLE */
}

static void draw_vport(const draw_data_t* dd) {
	cairo_set_source_surface(dd->cr, com.surface[dd->vport], 0, 0);
	cairo_pattern_set_filter(cairo_get_source(dd->cr), dd->filter);
	cairo_paint(dd->cr);
}

/* This method setup the viewport coordinates and draw the main image
 */
static void draw_main_image(const draw_data_t* dd) {
	if ((dd->vport == RGB_VPORT && com.rgbbuf) || com.graybuf[dd->vport]) {
		cairo_transform(dd->cr, &com.display_matrix);
		draw_vport(dd);
	} else {
		// For empty image, coordinates are untouched
		draw_empty_image(dd);
	}
}

static void draw_selection(const draw_data_t* dd) {
	if (com.selection.w > 0 && com.selection.h > 0) {
		cairo_t *cr = dd->cr;
		static double dash_format[] = { 4.0, 2.0 };
		cairo_set_line_width(cr, 1.5 / dd->zoom);
		cairo_set_dash(cr, dash_format, 2, 0);
		cairo_set_source_rgb(cr, 0.8, 1.0, 0.8);
		cairo_rectangle(cr, (double) com.selection.x, (double) com.selection.y,
						(double) com.selection.w, (double) com.selection.h);
		cairo_stroke(cr);

		// display a grid when the selection is being made / modified, when it is big enough
		if (com.pref.selection_guides > 1 && com.drawing && com.selection.w > 40 / dd->zoom && com.selection.h > 40 / dd->zoom) {
			cairo_set_line_width(cr, 0.4 / dd->zoom);
			cairo_set_dash(cr, NULL, 0, 0);
			for (int i = 1; i < com.pref.selection_guides; i++) {
				int x = com.selection.x + com.selection.w * i / com.pref.selection_guides;
				int y = com.selection.y + com.selection.h * i / com.pref.selection_guides;
				cairo_move_to(cr, x, com.selection.y);
				cairo_line_to(cr, x, com.selection.y + com.selection.h);
				cairo_move_to(cr, com.selection.x, y);
				cairo_line_to(cr, com.selection.x + com.selection.w, y);
			}
			cairo_stroke(cr);
		}

		// display a mini cross when the selection is being dragged
		if (com.freezeX && com.freezeY) {
			cairo_set_line_width(cr, 1.0 / dd->zoom);
			point selection_center = { com.selection.x + com.selection.w / 2,
				com.selection.y + com.selection.h / 2 };
			cairo_move_to(cr, selection_center.x, selection_center.y - 2 / dd->zoom);
			cairo_line_to(cr, selection_center.x, selection_center.y + 2 / dd->zoom);
			cairo_move_to(cr, selection_center.x - 2 / dd->zoom, selection_center.y);
			cairo_line_to(cr, selection_center.x + 2 / dd->zoom, selection_center.y);
			cairo_stroke(cr);
		}
	}
}

static void draw_stars(const draw_data_t* dd) {
	cairo_t *cr = dd->cr;
	int i = 0;

	if (com.stars && !com.script) {
		/* com.stars is a NULL-terminated array */
		cairo_set_dash(cr, NULL, 0, 0);
		cairo_set_source_rgba(cr, 1.0, 0.4, 0.0, 0.9);
		cairo_set_line_width(cr, 1.5 / dd->zoom);

		while (com.stars[i]) {
			double size = com.stars[i]->sx;
			if (com.stars[i]->units[0] == '\"') {
				size = convert_single_fwhm_to_pixels(com.stars[i]->fwhmx, com.stars[i]->sx);
			}

			if (i == com.selected_star) {
				// We draw horizontal and vertical lines to show the star
				cairo_set_line_width(cr, 2.0 / dd->zoom);
				cairo_set_source_rgba(cr, 0.0, 0.4, 1.0, 0.6);

				cairo_move_to(cr, com.stars[i]->xpos, 0);
				cairo_line_to(cr, com.stars[i]->xpos, dd->image_height);
				cairo_stroke(cr);
				cairo_move_to(cr, 0, com.stars[i]->ypos);
				cairo_line_to(cr, dd->image_width, com.stars[i]->ypos);
				cairo_stroke(cr);

				cairo_set_source_rgba(cr, 1.0, 0.4, 0.0, 0.9);
				cairo_set_line_width(cr, 1.5 / dd->zoom);
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
			cairo_set_line_width(cr, 2.0 / dd->zoom);
			fitted_PSF *the_psf = com.seq.photometry[i][com.seq.current];
			if (the_psf) {
				double size = the_psf->sx;
				cairo_arc(cr, the_psf->xpos, the_psf->ypos, size, 0., 2. * M_PI);
				cairo_stroke(cr);
				cairo_arc(cr, the_psf->xpos, the_psf->ypos, com.pref.phot_set.inner, 0.,
						  2. * M_PI);
				cairo_stroke(cr);
				cairo_arc(cr, the_psf->xpos, the_psf->ypos, com.pref.phot_set.outer, 0.,
						  2. * M_PI);
				cairo_stroke(cr);
				cairo_select_font_face(cr, "Purisa", CAIRO_FONT_SLANT_NORMAL,
									   CAIRO_FONT_WEIGHT_BOLD);
				cairo_set_font_size(cr, 40);
				cairo_move_to(cr, the_psf->xpos + com.pref.phot_set.outer + 5, the_psf->ypos);
				if (i == 0) {
					cairo_show_text(cr, "V");
				}
				else {
					gchar *tmp;
					tmp = g_strdup_printf("%d", i);
					cairo_show_text(cr, tmp);

					g_free(tmp);
				}
				cairo_stroke(cr);
			}
		}

		/* draw a cross on excluded images */
		if (com.seq.imgparam && com.seq.current >= 0 &&
				!com.seq.imgparam[com.seq.current].incl) {
			int w = dd->image_width > gfit.rx ? gfit.rx : dd->image_width;
			int h = dd->image_height > gfit.ry ? gfit.ry : dd->image_height;
			cairo_set_dash(cr, NULL, 0, 0);
			cairo_set_source_rgb(cr, 1.0, 0.8, 0.7);
			cairo_set_line_width(cr, 2.0 / dd->zoom);
			cairo_move_to(cr, 0, 0);
			cairo_line_to(cr, w, h);
			cairo_move_to(cr, 0, h);
			cairo_line_to(cr, w, 0.0);
			cairo_stroke(cr);
		}

		/* draw preview rectangles for the manual registration */
		for (i = 0; i < PREVIEW_NB; i++) {
			if (com.seq.previewX[i] >= 0) {
				int textX, textY;
				gchar *text;
				cairo_set_line_width(cr, 1.0 / dd->zoom);
				cairo_set_source_rgb(cr, 0.1, 0.6, 0.0);
				cairo_rectangle(cr,
						com.seq.previewX[i] - com.seq.previewW[i] / 2,
						com.seq.previewY[i] - com.seq.previewH[i] / 2,
						com.seq.previewW[i], com.seq.previewH[i]);
				cairo_stroke(cr);

				textX = com.seq.previewX[i] - com.seq.previewW[i] / 2;
				textX += 0.1 * com.seq.previewW[i];

				textY = com.seq.previewY[i] - com.seq.previewH[i] / 2;
				textY += 0.1 * com.seq.previewH[i];

				text = g_strdup_printf("%d", i + 1);

				cairo_set_font_size(cr, 12.0 / dd->zoom);
				cairo_move_to(cr, textX, textY);
				cairo_show_text(cr, text);
				g_free(text);
			}
		}
	}
}

static void draw_brg_boxes(const draw_data_t* dd) {
	GSList *list;
	for (list = com.grad_samples; list; list = list->next) {
		background_sample *sample = (background_sample *)list->data;
		if (background_sample_is_valid(sample)) {
			int radius = (int) (background_sample_get_size(sample) / 2);
			point position = background_sample_get_position(sample);
			cairo_set_line_width(dd->cr, 1.5 / dd->zoom);
			cairo_set_source_rgba(dd->cr, 0.2, 1.0, 0.3, 1.0);
			cairo_rectangle(dd->cr, position.x - radius - 1, position.y - radius,
					radius * 2, radius * 2);
			cairo_stroke(dd->cr);
		}
	}
}

static gdouble x_circle(gdouble x, gdouble radius) {
	return x + radius * cos(315 * M_PI / 180);
}

static gdouble y_circle(gdouble y, gdouble radius) {
	return y + radius * sin(315 * M_PI / 180);
}

static void draw_annotates(const draw_data_t* dd) {
	if (!com.found_object) return;
	cairo_t *cr = dd->cr;
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_set_source_rgba(cr, 0.5, 1.0, 0.3, 0.9);
	cairo_set_line_width(cr, 1.5 / dd->zoom);
	GSList *list;
	for (list = com.found_object; list; list = list->next) {
		CatalogObjects *object = (CatalogObjects *)list->data;
		gdouble radius = get_catalogue_object_radius(object);
		gdouble world_x = get_catalogue_object_ra(object);
		gdouble world_y = get_catalogue_object_dec(object);
		gchar *code = get_catalogue_object_code(object);
		gdouble resolution = get_wcs_image_resolution();
		gdouble x, y;
		gdouble size = 18 * (com.pref.font_scale / 100.0);

		if (resolution <= 0) return;

		radius = radius / resolution / 60.0;

		wcs2pix(world_x, world_y, &x, &y);
		y = gfit.ry - y;

		if (x > 0 && x < gfit.rx && y > 0 && y < gfit.ry) {
			point offset = {10, -10};
			if (radius < 0) {
				// objects we don't have an accurate location (LdN, Sh2)
			} else if (radius > 5) {
				cairo_arc(cr, x, y, radius, 0., 2. * M_PI);
				cairo_stroke(cr);
				cairo_move_to(cr, x_circle(x, radius), y_circle(y, radius));
				offset.x = x_circle(x, radius * 1.3) - x;
				offset.y = y_circle(y, radius * 1.3) - y;
				cairo_line_to(cr, offset.x + x, offset.y + y);
			} else {
				/* it is punctual */
				cairo_move_to(cr, x, y - 20);
				cairo_line_to(cr, x, y - 10);
				cairo_stroke(cr);
				cairo_move_to(cr, x, y + 20);
				cairo_line_to(cr, x, y + 10);
				cairo_stroke(cr);
				cairo_move_to(cr, x - 20, y);
				cairo_line_to(cr, x - 10, y);
				cairo_stroke(cr);
				cairo_move_to(cr, x + 20, y);
				cairo_line_to(cr, x + 10, y);
				cairo_stroke(cr);
			}
			if (code) {
				cairo_set_font_size(cr, size / dd->zoom);
				cairo_move_to(cr, x + offset.x, y + offset.y);
				cairo_show_text(cr, code);
				cairo_stroke(cr);
			}
		}
	}
}

static gboolean redraw_idle(gpointer p) {
	redraw(com.cvport, GPOINTER_TO_INT(p)); // draw stars
	return FALSE;
}

void initialize_image_display() {
	int i;
	for (i = 0; i < MAXGRAYVPORT; i++) {
		remap_index[i] = NULL;
		last_pente[i] = 0.f;
		last_mode[i] = HISTEQ_DISPLAY;
		// only HISTEQ mode always computes the index, it's a good initializer here
	}

	/* Create pixbuf from siril.svg file */
	gchar *image = g_build_filename(siril_get_system_data_dir(), "pixmaps", "siril.svg", NULL);
	com.siril_pix = gdk_pixbuf_new_from_file(image, NULL);
	g_free(image);

	cairo_matrix_init_identity(&com.display_matrix);
}

/* this function calculates the "fit to window" zoom values, given the window
 * size in argument and the image size in gfit.
 * Should not be called before displaying the main gray window when using zoom to fit */
double get_zoom_val() {
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
	return min(wtmp, htmp);
}

void adjust_vport_size_to_image() {
	if (com.script) return;
	double zoom = get_zoom_val();
	if (zoom <= 0.0) return;
	/* Init display matrix from current display state */
	cairo_matrix_init(&com.display_matrix,
					  zoom, 0, 0, zoom,
					  com.display_offset.x,
					  com.display_offset.y);
	/* Compute the inverse display matrix used for coordinate transformation */
	com.image_matrix = com.display_matrix;
	cairo_matrix_invert(&com.image_matrix);
}

static const gchar *label_zoom[] = { "labelzoom_red", "labelzoom_green", "labelzoom_blue", "labelzoom_rgb"};

static gboolean set_label_zoom_text_idle(gpointer p) {
	const gchar *txt = (const gchar *) p;
	GtkLabel *label = GTK_LABEL(lookup_widget(label_zoom[com.cvport]));

	gtk_label_set_text(label, txt);
	return FALSE;
}

static void update_zoom_label(gdouble zoom) {
	static gchar zoom_buffer[256] = { 0 };
	if ((single_image_is_loaded() || sequence_is_loaded()) && com.cvport < RGB_VPORT) {
		if (zoom < 0) {
			zoom = get_zoom_val();
		}
		g_sprintf(zoom_buffer, "%d%%", (int) (zoom * 100.0));
		gdk_threads_add_idle(set_label_zoom_text_idle, zoom_buffer);
	}
}

void redraw(int vport, int doremap) {
	if (com.script) return;

	update_zoom_label(com.zoom_value);

	if (vport >= MAXVPORT) {
		siril_debug_print(_("redraw: maximum number of layers supported is %d"
					" (current image has %d).\n"), MAXVPORT, vport);
		return;
	}
	GtkWidget *widget = com.vport[vport];

	if (doremap == REMAP_ALL) {
		stfComputed = FALSE;
		for (int i = 0; i < gfit.naxes[2]; i++) {
			remap(i);
		}
		if (gfit.naxis == 3)
			remaprgb();
		gtk_widget_queue_draw(widget);
		return;
	}

	switch (vport) {
		case RED_VPORT:
		case BLUE_VPORT:
		case GREEN_VPORT:
			if (doremap == REMAP_ONLY) {
				remap(vport);
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
			siril_debug_print("redraw: unknown viewport number %d\n", vport);
			break;
	}
}

void queue_redraw(int doremap) {
	// request a redraw from another thread
	siril_add_idle(redraw_idle, GINT_TO_POINTER(doremap));
}


/* callback for GtkDrawingArea, draw event
 * see http://developer.gnome.org/gtk3/3.2/GtkDrawingArea.html
 * http://developer.gnome.org/gdk-pixbuf/stable/gdk-pixbuf-Image-Data-in-Memory.html
 * http://www.cairographics.org/manual/
 * http://www.cairographics.org/manual/cairo-Image-Surfaces.html#cairo-image-surface-create-for-data
 */
gboolean redraw_drawingarea(GtkWidget *widget, cairo_t *cr, gpointer data) {
	draw_data_t dd;

	// we need to identify which vport is being redrawn
	dd.vport = match_drawing_area_widget(widget, TRUE);
	if (dd.vport == -1) {
		fprintf(stderr, "Could not find the vport for the draw callback\n");
		return TRUE;
	}

	/* catch and compute rendering data */
	dd.cr = cr;
	dd.window_width = gtk_widget_get_allocated_width(widget);
	dd.window_height = gtk_widget_get_allocated_height(widget);
	dd.zoom = get_zoom_val();
	dd.image_width = gfit.rx;
	dd.image_height = gfit.ry;
	dd.filter = (dd.zoom < 1.0) ? CAIRO_FILTER_GOOD : CAIRO_FILTER_FAST;

	adjust_vport_size_to_image();

	cairo_save(cr);

	/* RGB or gray images */
	draw_main_image(&dd);

	/* selection rectangle */
	draw_selection(&dd);

	/* detected stars and highlight the selected star */
	draw_stars(&dd);

	/* detected objects */
	draw_annotates(&dd);

	/* background removal gradient selection boxes */
	draw_brg_boxes(&dd);

	cairo_restore(cr);

	return FALSE;
}

point get_center_of_vport() {
	GtkWidget *widget = lookup_widget("drawingarear");

	guint window_width = gtk_widget_get_allocated_width(widget);
	guint window_height = gtk_widget_get_allocated_height(widget);

	point center = { window_width / 2, window_height / 2 };

	return center;
}

void add_image_and_label_to_cairo(cairo_t *cr, int vport) {
	draw_data_t dd;

	GtkWidget *widget = lookup_widget("drawingarear");

	dd.vport = vport;
	dd.cr = cr;
	dd.window_width = gtk_widget_get_allocated_width(widget);
	dd.window_height = gtk_widget_get_allocated_height(widget);
	dd.zoom = get_zoom_val();
	dd.image_width = gfit.rx;
	dd.image_height = gfit.ry;
	dd.filter = (dd.zoom < 1.0) ? CAIRO_FILTER_GOOD : CAIRO_FILTER_FAST;

	/* RGB or gray images */
	draw_main_image(&dd);
	/* detected objects */
	draw_annotates(&dd);
}
