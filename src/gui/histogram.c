/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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

#include <gsl/gsl_histogram.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/siril_app_dirs.h"
#include "algos/statistics.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "gui/image_display.h"
#include "gui/callbacks.h"	// for lookup_widget()
#include "gui/progress_and_log.h"
#include "gui/dialogs.h"
#include "gui/siril_preview.h"
#include "core/undo.h"

#include "histogram.h"

#define shadowsClipping -2.80f /* Shadows clipping point measured in sigma units from the main histogram peak. */
#define targetBackground 0.25f /* final "luminance" of the image for autostretch in the [0,1] range */
#define GRADIENT_HEIGHT 12

#undef HISTO_DEBUG

/* The gsl_histogram, documented here:
 * http://linux.math.tifr.res.in/manuals/html/gsl-ref-html/gsl-ref_21.html
 * is able to group values into bins, it does not need to handle all values
 * separately. That is useful for display purposes, but not currently used.
 */

// margin between axis and drawing area border
// colors of layers histograms		R	G	B	RGB
static double histo_color_r[] = { 1.0, 0.0, 0.0, 0.0 };
static double histo_color_g[] = { 0.0, 1.0, 0.0, 0.0 };
static double histo_color_b[] = { 0.0, 0.0, 1.0, 0.0 };
static float graph_height = 0.f;	// the max value of all bins
static guint64 clipped[] = { 0, 0 };

static GtkToggleToolButton *toggles[MAXVPORT] = { NULL };
static GtkToggleToolButton *toggleGrid = NULL, *toggleCurve = NULL;

/* the original histogram, used as starting point of each computation */
static gsl_histogram *hist_backup[MAXVPORT] = { NULL, NULL, NULL, NULL };

static float _midtones, _shadows, _highlights;

static gboolean _click_on_histo = FALSE;
static ScaleType _type_of_scale;

static void apply_mtf_to_fits(fits *from, fits *to);
static void set_histogram(gsl_histogram *histo, int layer);

static int get_width_of_histo() {
	return gtk_widget_get_allocated_width(lookup_widget("drawingarea_histograms"));
}

static int get_height_of_histo() {
	return gtk_widget_get_allocated_height(lookup_widget("drawingarea_histograms"));
}

static void clear_hist_backup() {
	if (hist_backup[0]) {
		for (int i = 0; i < gfit.naxes[2]; i++) {
			gsl_histogram_free(hist_backup[i]);
			hist_backup[i] = NULL;
		}
	}
}

static void histo_startup() {
	copy_gfit_to_backup();
	// also get the backup histogram
	compute_histo_for_gfit();
	for (int i = 0; i < gfit.naxes[2]; i++)
		hist_backup[i] = gsl_histogram_clone(com.layers_hist[i]);
}

static void histo_close(gboolean revert) {
	int i;
	if (revert) {
		set_cursor_waiting(TRUE);

		for (i = 0; i < gfit.naxes[2]; i++) {
			set_histogram(hist_backup[i], i);
			hist_backup[i] = NULL;
		}
		copy_backup_to_gfit();
		adjust_cutoff_from_updated_gfit();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
		set_cursor_waiting(FALSE);
	}

	// free data
	clear_backup();
	clear_hist_backup();
}

static void histo_recompute() {
	set_cursor("progress");
	copy_backup_to_gfit();

	apply_mtf_to_fits(get_preview_gfit_backup(), &gfit);
	// com.layers_hist should be good, update_histo_mtf() is always called before

	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
}

static void _init_clipped_pixels() {
	clipped[0] = 0;
	clipped[1] = 0;
}

static void _initialize_clip_text() {
	static GtkEntry *clip_high = NULL, *clip_low = NULL;

	if (clip_high == NULL) {
		clip_high = GTK_ENTRY(lookup_widget("clip_highlights"));
		clip_low = GTK_ENTRY(lookup_widget("clip_shadows"));
	}
	gtk_entry_set_text(clip_low, "0.000%");
	gtk_entry_set_text(clip_high, "0.000%");
}

static void _update_entry_text() {
	GtkEntry *histoMidEntry = GTK_ENTRY(lookup_widget("histoMidEntry"));
	GtkEntry *histoShadEntry = GTK_ENTRY(lookup_widget("histoShadEntry"));
	GtkEntry *histoHighEntry = GTK_ENTRY(lookup_widget("histoHighEntry"));
	gchar *buffer;

	buffer = g_strdup_printf("%.7f", _shadows);
	gtk_entry_set_text(histoShadEntry, buffer);
	g_free(buffer);
	buffer = g_strdup_printf("%.7f", _highlights);
	gtk_entry_set_text(histoHighEntry, buffer);
	g_free(buffer);
	buffer = g_strdup_printf("%.7f", _midtones);
	gtk_entry_set_text(histoMidEntry, buffer);
	g_free(buffer);
}

static void _update_clipped_pixels(size_t data) {
	static GtkEntry *clip_high = NULL, *clip_low = NULL;
	double tmp;
	char buffer[16];

	if (clip_high == NULL) {
		clip_high = GTK_ENTRY(lookup_widget("clip_highlights"));
		clip_low = GTK_ENTRY(lookup_widget("clip_shadows"));
	}
	tmp = (double)clipped[1] * 100.0 / (double)data;
	g_snprintf(buffer, sizeof(buffer), "%.3f%%", tmp);
	gtk_entry_set_text(clip_high, buffer);
	tmp = (double)clipped[0] * 100.0 / (double)data;
	g_snprintf(buffer, sizeof(buffer), "%.3f%%", tmp);
	gtk_entry_set_text(clip_low, buffer);
}

static int is_histogram_visible() {
	GtkWidget *window = lookup_widget("histogram_dialog");
	return gtk_widget_get_visible(window);
}

static void init_toggles() {
	if (!toggles[0]) {
		toggles[0] = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("histoToolRed"));
		toggles[1] = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("histoToolGreen"));
		toggles[2] = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("histoToolBlue"));
		toggles[3] = NULL;
		toggleGrid = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("histoToolGrid"));
		toggleCurve = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("histoToolCurve"));
	}
}

// sets the channel names of the toggle buttons in the histogram window, based on
// the number of layers of gfit
static void set_histo_toggles_names() {
	gchar *image;

	init_toggles();

	if (gfit.naxis == 2) {
		gtk_widget_set_tooltip_text(GTK_WIDGET(toggles[0]), _("Gray channel"));
		GtkWidget *w;
		if (com.pref.combo_theme == 0) {
			image = g_build_filename(siril_get_system_data_dir(), "pixmaps", "monochrome_dark.png", NULL);
			w = gtk_image_new_from_file(image);
		} else {
			image = g_build_filename(siril_get_system_data_dir(), "pixmaps", "monochrome.png", NULL);
			w = gtk_image_new_from_file(image);
		}
		gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(toggles[0]), w);
		gtk_widget_show(w);
		gtk_toggle_tool_button_set_active(toggles[0], TRUE);
		gtk_widget_set_visible(GTK_WIDGET(toggles[1]), FALSE);
		gtk_widget_set_visible(GTK_WIDGET(toggles[2]), FALSE);
		/* visible has no effect in GTK+ 3.12, trying sensitive too
		 * Yes it does. The solution is to call the window (widget)
		 * with gtk_widget_show and not gtk_widget_show_all */
		gtk_widget_set_sensitive(GTK_WIDGET(toggles[1]), FALSE);
		gtk_widget_set_sensitive(GTK_WIDGET(toggles[2]), FALSE);
		if (toggles[3])
			gtk_widget_set_visible(GTK_WIDGET(toggles[3]), FALSE);

	} else {
		image = g_build_filename(siril_get_system_data_dir(), "pixmaps", "r.png", NULL);
		gtk_widget_set_tooltip_text(GTK_WIDGET(toggles[0]), _("Red channel"));
		GtkWidget *w = gtk_image_new_from_file(image);
		gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(toggles[0]), w);
		gtk_widget_show(w);
		gtk_toggle_tool_button_set_active(toggles[0], TRUE);
		gtk_toggle_tool_button_set_active(toggles[1], TRUE);
		gtk_toggle_tool_button_set_active(toggles[2], TRUE);
		gtk_widget_set_sensitive(GTK_WIDGET(toggles[1]), TRUE);
		gtk_widget_set_sensitive(GTK_WIDGET(toggles[2]), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(toggles[1]), TRUE);
		gtk_widget_set_visible(GTK_WIDGET(toggles[2]), TRUE);
		if (toggles[3]) {
			gtk_widget_set_visible(GTK_WIDGET(toggles[3]), TRUE);
			gtk_toggle_tool_button_set_active(toggles[3], TRUE);
		}
	}
	g_free(image);
}

static double get_histoZoomValueH() {
	static GtkAdjustment *histoAdjZoomH = NULL;
	if (!histoAdjZoomH)
		histoAdjZoomH = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "histoAdjZoomH"));

	return gtk_adjustment_get_value(histoAdjZoomH);
}

static double get_histoZoomValueV() {
	static GtkAdjustment *histoAdjZoomV = NULL;
	if (!histoAdjZoomV)
		histoAdjZoomV = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "histoAdjZoomV"));

	return gtk_adjustment_get_value(histoAdjZoomV);
}

static void adjust_histogram_vport_size() {
	GtkWidget *drawarea, *vport;
	int targetW, targetH, cur_width, cur_height;
	double zoomH = get_histoZoomValueH();
	double zoomV = get_histoZoomValueV();

	drawarea = lookup_widget("drawingarea_histograms");
	vport = lookup_widget("viewport1");

	cur_width = gtk_widget_get_allocated_width(vport);
	cur_height = gtk_widget_get_allocated_height(vport);
	targetW = (int) (((double)cur_width) * zoomH);
	targetH = (int) (((double)cur_height) * zoomV);
	gtk_widget_set_size_request(drawarea, targetW, targetH);
#ifdef HISTO_DEBUG
	fprintf(stdout, "Histo vport size (%d, %d)\n", targetW, targetH);
#endif
}

size_t get_histo_size(fits *fit) {
	if (fit->type == DATA_USHORT)
		return (size_t)get_normalized_value(fit);
	return (size_t)USHRT_MAX;
}

// create a new hitogram object for the passed fit and layer
gsl_histogram* computeHisto(fits *fit, int layer) {
	g_assert(layer < 3);
	size_t i, ndata, size;

	size = get_histo_size(fit);
	gsl_histogram *histo = gsl_histogram_alloc(size + 1);
	gsl_histogram_set_ranges_uniform(histo, 0, fit->type == DATA_FLOAT ? 1.0 + 1.0 / size : size + 1);
	ndata = fit->naxes[0] * fit->naxes[1];

#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread)
#endif
	{
		gsl_histogram *histo_thr = gsl_histogram_alloc(size + 1);
		gsl_histogram_set_ranges_uniform(histo_thr, 0, fit->type == DATA_FLOAT ? 1.0 + 1.0 / size : size + 1);

		if (fit->type == DATA_USHORT) {
			WORD *buf = fit->pdata[layer];
#ifdef _OPENMP
#pragma omp for private(i) schedule(static)
#endif
			for (i = 0; i < ndata; i++) {
				gsl_histogram_increment(histo_thr, (double) buf[i]);
			}
		} else if (fit->type == DATA_FLOAT) {
			float *buf = fit->fpdata[layer];
#ifdef _OPENMP
#pragma omp for private(i) schedule(static)
#endif
			for (i = 0; i < ndata; i++) {
				gsl_histogram_increment(histo_thr, (double) buf[i]);
			}
		}
#ifdef _OPENMP
#pragma omp critical
#endif
		{
			gsl_histogram_add(histo, histo_thr);
		}
		gsl_histogram_free(histo_thr);
	}

	return histo;
}

static void draw_curve(cairo_t *cr, int width, int height) {
	// draw curve
	int k;
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_set_line_width(cr, 1.0);
	cairo_set_source_rgb(cr, .9, .9, .9);

	for (k = 0; k < width + 1; k++) {
		float x = k / (float) width;
		float y = MTF(x, _midtones, _shadows, _highlights);
		cairo_line_to(cr, k, height * (1 - y));
	}
	cairo_stroke(cr);
}

static void draw_grid(cairo_t *cr, int width, int height) {
	double dash_format[] = { 1.0, 1.0 };

	cairo_set_line_width(cr, 1.0);
	cairo_set_source_rgb(cr, 0.4, 0.4, 0.4);
	// quarters in solid, eights in dashed line
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_move_to(cr, width * 0.25, 0);
	cairo_line_to(cr, width * 0.25, height);
	cairo_move_to(cr, width * 0.5, 0);
	cairo_line_to(cr, width * 0.5, height);
	cairo_move_to(cr, width * 0.75, 0);
	cairo_line_to(cr, width * 0.75, height);

	cairo_set_line_width(cr, 1.0);
	cairo_move_to(cr, 0, height * 0.25);
	cairo_line_to(cr, width, height * 0.25);
	cairo_move_to(cr, 0, height * 0.5);
	cairo_line_to(cr, width, height * 0.5);
	cairo_move_to(cr, 0, height * 0.75);
	cairo_line_to(cr, width, height * 0.75);

	cairo_stroke(cr);

	cairo_set_line_width(cr, 1.0);
	cairo_set_dash(cr, dash_format, 2, 0);
	cairo_move_to(cr, width * 0.125, 0);
	cairo_line_to(cr, width * 0.125, height);
	cairo_move_to(cr, width * 0.375, 0);
	cairo_line_to(cr, width * 0.375, height);
	cairo_move_to(cr, width * 0.625, 0);
	cairo_line_to(cr, width * 0.625, height);
	cairo_move_to(cr, width * 0.875, 0);
	cairo_line_to(cr, width * 0.875, height);

	cairo_set_line_width(cr, 1.0);
	cairo_move_to(cr, 0, height * 0.125);
	cairo_line_to(cr, width, height * 0.125);
	cairo_move_to(cr, 0, height * 0.375);
	cairo_line_to(cr, width, height * 0.375);
	cairo_move_to(cr, 0, height * 0.625);
	cairo_line_to(cr, width, height * 0.625);
	cairo_move_to(cr, 0, height * 0.875);
	cairo_line_to(cr, width, height * 0.875);
	cairo_stroke(cr);

}

// erase image and redraw the background color and grid
static void erase_histo_display(cairo_t *cr, int width, int height) {
	gboolean drawGrid, drawCurve;
	// clear all with background color
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_rectangle(cr, 0, 0, width, height);
	cairo_fill(cr);

	// draw grid
	init_toggles();
	drawGrid = gtk_toggle_tool_button_get_active(toggleGrid);
	drawCurve = gtk_toggle_tool_button_get_active(toggleCurve);
	if (drawGrid)
		draw_grid(cr, width, height);
	if (drawCurve)
		draw_curve(cr, width, height);
}

static gboolean is_log_scale() {
	static GtkToggleButton *HistoCheckLogButton = NULL;

	if (HistoCheckLogButton == NULL)
		HistoCheckLogButton = GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckLogButton"));
	return (gtk_toggle_button_get_active(HistoCheckLogButton));
}


static void draw_gradient(cairo_t *cr, int width, int height) {
	cairo_pattern_t *pat;

	pat = cairo_pattern_create_linear(0.0, 0.0, width, 0.0);
	cairo_pattern_add_color_stop_rgb(pat, 0, 0, 0, 0);
	cairo_pattern_add_color_stop_rgb(pat, width, 1, 1, 1);
	cairo_rectangle(cr, 0, height - GRADIENT_HEIGHT, width, GRADIENT_HEIGHT);
	cairo_set_source(cr, pat);
	cairo_fill(cr);

	cairo_pattern_destroy(pat);
}

static void draw_slider(cairo_t *cr, int width, int height, int xpos) {
	if (xpos > width / 2) {
		cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
	} else {
		cairo_set_source_rgb(cr, 0.9, 0.9, 0.9);
	}
	cairo_move_to(cr, -10 + xpos, height);
	cairo_line_to(cr, 10 + xpos, height);
	cairo_line_to(cr, xpos, height - GRADIENT_HEIGHT);
	cairo_line_to(cr, -10 + xpos, height);
	cairo_stroke(cr);
}

static void display_scale(cairo_t *cr, int width, int height) {
	draw_gradient(cr, width, height);
	float delta = ((_highlights - _shadows) * _midtones) + _shadows;
	draw_slider(cr, width, height, _shadows * width);
	draw_slider(cr, width, height, (delta) * width);
	draw_slider(cr, width, height, _highlights * width);
}

static void display_histo(gsl_histogram *histo, cairo_t *cr, int layer, int width,
		int height, double zoomH, double zoomV) {
	if (width <= 0) return;
	int current_bin;
	size_t norm = gsl_histogram_bins(histo) - 1;

	float vals_per_px = (float)norm / (float)width;	// size of a bin
	size_t i, nb_orig_bins = gsl_histogram_bins(histo);

	// We need to store the binned histogram in order to find the binned maximum
	static gfloat *displayed_values = NULL;
	static int nb_bins_allocated = 0;
	/* we create a bin for each pixel in the displayed width.
	 * nb_bins_allocated is thus equal to the width of the image */
	if (nb_bins_allocated != width) {
		gfloat *tmp;
		nb_bins_allocated = width;
		tmp = realloc(displayed_values, nb_bins_allocated * sizeof(gfloat));
		if (!tmp) {
			if (displayed_values != NULL) {
				g_free(displayed_values);
				displayed_values = NULL;
			}
			PRINT_ALLOC_ERR;
			histo_close(TRUE);
			return;
		}
		displayed_values = tmp;
		memset(displayed_values, 0, nb_bins_allocated);
	}

	if (gfit.naxis == 2)
		cairo_set_source_rgb(cr, 255.0, 255.0, 255.0);
	else
		cairo_set_source_rgb(cr, histo_color_r[layer], histo_color_g[layer],
				histo_color_b[layer]);
	cairo_set_dash(cr, NULL, 0, 0);
	cairo_set_line_width(cr, 1.5);

	// first loop builds the bins and finds the maximum
	i = 0;
	current_bin = 0;
	do {
		float bin_val = 0.f;
		while (i < nb_orig_bins
				&& (float)i / vals_per_px <= (float)current_bin + 0.5f) {
			bin_val += (float)gsl_histogram_get(histo, i);
			i++;
		}
		if (is_log_scale() && bin_val != 0.f) {
			bin_val = logf(bin_val);
		}
		displayed_values[current_bin] = bin_val;
		if (bin_val > graph_height)	// check for maximum
			graph_height = bin_val;
		current_bin++;
	} while (i < nb_orig_bins && current_bin < nb_bins_allocated);
	for (i = 0; i < nb_bins_allocated; i++) {
		float bin_height = height - height * displayed_values[i] / graph_height;
		cairo_line_to(cr, i, bin_height);
	}
	cairo_stroke(cr);
}

static void apply_mtf_to_fits(fits *from, fits *to) {

	g_assert(from->naxes[2] == 1 || from->naxes[2] == 3);
	const size_t ndata = from->naxes[0] * from->naxes[1] * from->naxes[2];
	g_assert(from->type == to->type);

	if (from->type == DATA_USHORT) {
		float norm = (float)get_normalized_value(from);
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0; i < ndata; i++) {
			float pxl = (float)from->data[i] / norm;
			float mtf = MTF(pxl, _midtones, _shadows, _highlights);
			to->data[i] = round_to_WORD(mtf * norm);
		}
	}
	else if (from->type == DATA_FLOAT) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) schedule(static)
#endif
		for (size_t i = 0; i < ndata; i++) {
			to->fdata[i] = MTF(from->fdata[i], _midtones, _shadows, _highlights);
		}
	}
	else return;

	invalidate_stats_from_fit(to);
}

static void apply_mtf_to_histo(gsl_histogram *histo, float norm,
		float m, float lo, float hi) {

	size_t int_norm = (size_t)norm;
	gsl_histogram *mtf_histo = gsl_histogram_alloc(int_norm + 1);

	gsl_histogram_set_ranges_uniform(mtf_histo, 0, norm);

// #pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static) // disabled because of ISSUE #136 (https://free-astro.org/bugs/view.php?id=136)
	for (size_t i = 0; i < int_norm + 1; i++) {
		WORD mtf;
		float binval = gsl_histogram_get(histo, i);
		float pxl = ((float)i / norm);
		guint64 clip[2] = { 0, 0 };

		if (i < round_to_WORD(lo * norm)) {
			pxl = lo;
			clip[0] += binval;
		} else if (i > round_to_WORD(hi * norm)) {
			pxl = hi;
			clip[1] += binval;
		}
		mtf = round_to_WORD(MTF(pxl, m, lo, hi) * norm);
		gsl_histogram_accumulate(mtf_histo, (double)mtf, (double)binval);
//#ifdef _OPENMP
//#pragma omp critical
//#endif
		{
			clipped[0] += clip[0];
			clipped[1] += clip[1];
		}
	}
	gsl_histogram_memcpy(histo, mtf_histo);
	gsl_histogram_free(mtf_histo);
}

static void reset_cursors_and_values() {
	_shadows = 0.f;
	_midtones = 0.5f;
	_highlights = 1.0f;
	graph_height = 0.f;

	_init_clipped_pixels();
	_initialize_clip_text();
	_update_entry_text();
	update_gfit_histogram_if_needed();
}

static void queue_window_redraw() {
	static GtkWidget *drawarea = NULL;
	if (!drawarea)
		drawarea = lookup_widget("drawingarea_histograms");
	gtk_widget_queue_draw(drawarea);
}

static void update_histo_mtf() {
	float norm = (float)gsl_histogram_bins(com.layers_hist[0]) - 1;

	_init_clipped_pixels();
	for (long i = 0; i < gfit.naxes[2]; i++) {
		gsl_histogram_memcpy(com.layers_hist[i], hist_backup[i]);
		apply_mtf_to_histo(com.layers_hist[i], norm, _midtones, _shadows, _highlights);
	}
	size_t data = gfit.naxes[0] * gfit.naxes[1] * gfit.naxes[2];
	_update_clipped_pixels(data);
	queue_window_redraw();
}

static void set_histogram(gsl_histogram *histo, int layer) {
	g_assert(layer >= 0 && layer < MAXVPORT);
	if (com.layers_hist[layer])
		gsl_histogram_free(com.layers_hist[layer]);
	com.layers_hist[layer] = histo;
}

static gboolean on_gradient(GdkEvent *event, int width, int height) {
	return ((GdkEventButton*) event)->x > 0
			&& ((GdkEventButton*) event)->x < width
			&& ((GdkEventButton*) event)->y > height - GRADIENT_HEIGHT
			&& ((GdkEventButton*) event)->y < height;
}

/*
 * Public functions
 */

void mtf_with_parameters(fits *fit, float lo, float mid, float hi) {
	_shadows = lo;
	_midtones = mid;
	_highlights = hi;
	apply_mtf_to_fits(fit, fit);
}

gsl_histogram* computeHisto_Selection(fits* fit, int layer,
		rectangle *selection) {
	g_assert(layer < 3);

	size_t size = get_histo_size(fit);
	gsl_histogram* histo = gsl_histogram_alloc(size + 1);
	gsl_histogram_set_ranges_uniform(histo, 0, fit->type == DATA_FLOAT ? 1.0 : size);
	size_t stridefrom = fit->rx - selection->w;

	if (fit->type == DATA_USHORT) {
		WORD *from = fit->pdata[layer] + (fit->ry - selection->y - selection->h) * fit->rx
			+ selection->x;
		for (size_t i = 0; i < selection->h; i++) {
			for (size_t j = 0; j < selection->w; j++) {
				gsl_histogram_increment(histo, (double)*from);
				from++;
			}
			from += stridefrom;
		}
	}
	else if (fit->type == DATA_FLOAT) {
		float *from = fit->fpdata[layer] + (fit->ry - selection->y - selection->h) * fit->rx
			+ selection->x;
		for (size_t i = 0; i < selection->h; i++) {
			for (size_t j = 0; j < selection->w; j++) {
				gsl_histogram_increment(histo, (double)*from);
				from++;
			}
			from += stridefrom;
		}
	}
	return histo;
}

void compute_histo_for_gfit() {
	int nb_layers = 3;
	if (gfit.naxis == 2)
		nb_layers = 1;
	for (int i = 0; i < nb_layers; i++) {
		if (!com.layers_hist[i])
			set_histogram(computeHisto(&gfit, i), i);
	}
	set_histo_toggles_names();
}

void invalidate_gfit_histogram() {
	for (int layer = 0; layer < MAXVPORT; layer++) {
		set_histogram(NULL, layer);
	}
}

void update_gfit_histogram_if_needed() {
	if (is_histogram_visible()) {
		compute_histo_for_gfit();
		queue_window_redraw();
	}
}

void clear_histograms() {
	for (int i = 0; i < MAXVPORT; i++) {
		set_histogram(NULL, i);
	}
	// TODO: call histo_close? should it be done by the caller?
}

float MTF(float x, float m, float lo, float hi) {
	if (x <= lo)
		return 0.f;
	if (x >= hi)
		return 1.f;

	float xp = (x - lo) / (hi - lo);

	return ((m - 1.f) * xp) / (((2.f * m - 1.f) * xp) - m);
}

float findMidtonesBalance(fits *fit, float *shadows, float *highlights) {
	float c0 = 0.0, c1 = 0.0;
	float m = 0.0;
	int i, n, invertedChannels = 0;
	imstats *stat[3];

	n = fit->naxes[2];

	for (i = 0; i < n; ++i) {
		stat[i] = statistics(NULL, -1, fit, i, NULL, STATS_BASIC | STATS_MAD, TRUE);
		if (!stat[i]) {
			siril_log_message(_("Error: statistics computation failed.\n"));
			return 0.f;
		}

		if (stat[i]->median / stat[i]->normValue > 0.5)
			++invertedChannels;
	}

	if (invertedChannels < n) {
		for (i = 0; i < n; ++i) {
			float median, mad, normValue;

			normValue = (float)stat[i]->normValue;
			median = (float) stat[i]->median / normValue;
			mad = (float) stat[i]->mad / normValue * (float)MAD_NORM;
			/* this is a guard to avoid breakdown point */
			if (mad == 0.f) mad = 0.001f;

			c0 += median + shadowsClipping * mad;
			m += median;
		}
		c0 /= (float) n;
		if (c0 < 0.f) c0 = 0.f;
		float m2 = m / (float) n - c0;
		m = MTF(m2, targetBackground, 0.f, 1.f);
		*shadows = c0;
		*highlights = 1.0;
	} else {
		for (i = 0; i < n; ++i) {
			float median, mad, normValue;

			normValue = (float) stat[i]->normValue;
			median = (float) stat[i]->median / normValue;
			mad = (float) stat[i]->mad / normValue * (float)MAD_NORM;
			/* this is a guard to avoid breakdown point */
			if (mad == 0.f) mad = 0.001f;

			m += median;
			c1 += median - shadowsClipping * mad;
		}
		c1 /= (float) n;
		if (c1 > 1.f) c1 = 1.f;
		float m2 = c1 - m / (float) n;
		m = 1.f - MTF(m2, targetBackground, 0.f, 1.f);
		*shadows = 0.f;
		*highlights = c1;

	}
	for (i = 0; i < n; ++i)
		free_stats(stat[i]);
	return m;
}

static int mtf_image_hook(struct generic_seq_args *args, int o, int i, fits *fit,
		rectangle *_) {
	struct mtf_data *m_args = (struct mtf_data*) args->user;

	mtf_with_parameters(fit, m_args->lo, m_args->mid, m_args->hi);
	return 0;
}

/* Callback functions */

gboolean redraw_histo(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int i, width, height;
	double zoomH, zoomV;

	init_toggles();
	width = get_width_of_histo();
	height = get_height_of_histo();
#ifdef HISTO_DEBUG
	fprintf(stdout, "histogram redraw\n");
	fprintf(stdout, "w = %d and h = %d\n", width, height);
#endif
	zoomH = get_histoZoomValueH();
	zoomV = get_histoZoomValueV();

	if (height == 1)
		return FALSE;
	erase_histo_display(cr, width, height - GRADIENT_HEIGHT);
	graph_height = 0.0;
	for (i = 0; i < MAXVPORT; i++) {
		if (com.layers_hist[i]
				&& (!toggles[i] || gtk_toggle_tool_button_get_active(toggles[i])))
			display_histo(com.layers_hist[i], cr, i, width, height - GRADIENT_HEIGHT, zoomH, zoomV);
	}
	display_scale(cr, width, height);
	return FALSE;
}

void on_histo_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	queue_window_redraw();
}

void on_histogram_window_show(GtkWidget *object, gpointer user_data) {
	histo_startup();
	_initialize_clip_text();
	reset_cursors_and_values();
	compute_histo_for_gfit();
}

void on_button_histo_close_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	reset_cursors_and_values();
	histo_close(TRUE);
	set_cursor_waiting(FALSE);
	siril_close_dialog("histogram_dialog");
}

void on_button_histo_reset_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	reset_cursors_and_values();
	histo_close(TRUE);
	histo_startup();
	set_cursor_waiting(FALSE);
}

gboolean on_scale_key_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	set_cursor_waiting(FALSE);
	return FALSE;
}

void on_button_histo_apply_clicked(GtkButton *button, gpointer user_data) {
	if ((_midtones == 0.5f) && (_shadows == 0.f) && (_highlights == 1.f)) {
		return;
	}

	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkMTFSeq")))
			&& sequence_is_loaded()) {
		/* Apply to the whole sequence */
		struct mtf_data *args = malloc(sizeof(struct mtf_data));

		args->lo = _shadows;
		args->mid = _midtones;
		args->hi = _highlights;
		args->seqEntry = gtk_entry_get_text(GTK_ENTRY(lookup_widget("entryMTFSeq")));
		if (args->seqEntry && args->seqEntry[0] == '\0')
			args->seqEntry = "mtf_";
		args->seq = &com.seq;
		/* here it is a bit tricky.
		 * It is better to first close the window as it is a liveview tool
		 * TODO: could we improve this behavior?
		 */
		reset_cursors_and_values();
		histo_close(TRUE);
		siril_close_dialog("histogram_dialog");

		/* apply the process */
		apply_mtf_to_sequence(args);

	} else {
		// the apply button resets everything after recomputing with the current values
		histo_recompute();
		// partial cleanup
		siril_debug_print("Applying histogram (mid=%.3f, lo=%.3f, hi=%.3f)\n",
				_midtones, _shadows, _highlights);
		undo_save_state(get_preview_gfit_backup(),
				_("Histogram Transf. (mid=%.3f, lo=%.3f, hi=%.3f)"), _midtones, _shadows,
				_highlights);

		clear_backup();
		clear_hist_backup();
		// reinit
		histo_startup();
		reset_cursors_and_values();

		set_cursor("default");
	}
}

void apply_histo_cancel() {
	set_cursor_waiting(TRUE);
	reset_cursors_and_values();
	histo_close(TRUE);
	set_cursor_waiting(FALSE);
}

void on_histoZoom100_clicked(GtkButton *button, gpointer user_data) {
	GtkAdjustment *histoAdjZoomH, *histoAdjZoomV;

	histoAdjZoomH = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "histoAdjZoomH"));
	histoAdjZoomV = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "histoAdjZoomV"));

	gtk_adjustment_set_value(histoAdjZoomH, 1.0);
	gtk_adjustment_set_value(histoAdjZoomV, 1.0);
}

void on_histoSpinZoom_value_changed(GtkRange *range, gpointer user_data) {
	adjust_histogram_vport_size();
	queue_window_redraw();
}

void on_histoToolAutoStretch_clicked(GtkToolButton *button, gpointer user_data) {
	float m, shadows = 0.f, highlights = 0.f;

	set_cursor_waiting(TRUE);
	/* we always apply this function on original data */
	m = findMidtonesBalance(get_preview_gfit_backup(), &shadows, &highlights);
	_shadows = shadows;
	_midtones = m;
	_highlights = 1.f;

	_update_entry_text();
	update_histo_mtf();
	histo_recompute();
	set_cursor_waiting(FALSE);
}

void on_menuitem_histo_activate(GtkMenuItem *menuitem, gpointer user_data) {
	set_cursor_waiting(TRUE);

	siril_open_dialog("histogram_dialog");
	set_cursor_waiting(FALSE);
}

void toggle_histogram_window_visibility(GtkToolButton *button, gpointer user_data) {
	if (gtk_widget_get_visible((GtkWidget *)user_data)) {
		set_cursor_waiting(TRUE);
		reset_cursors_and_values();
		histo_close(TRUE);
		set_cursor_waiting(FALSE);
		siril_close_dialog("histogram_dialog");
	} else {
		on_menuitem_histo_activate(NULL, NULL);
	}
}

gboolean on_drawingarea_histograms_motion_notify_event(GtkWidget *widget, GdkEventMotion *event,
		gpointer user_data) {

	int width = get_width_of_histo();
	int height = get_height_of_histo();

	if (on_gradient((GdkEvent *) event, width, height)) {
		set_cursor("grab");
	} else {
		set_cursor("default");
	}

	if (_click_on_histo) {
		gdouble xpos = ((GdkEventButton*) event)->x / (gdouble) width;
		gchar *buffer = NULL;
		GtkEntry *histoMidEntry = GTK_ENTRY(lookup_widget("histoMidEntry"));
		GtkEntry *histoShadEntry = GTK_ENTRY(lookup_widget("histoShadEntry"));
		GtkEntry *histoHighEntry = GTK_ENTRY(lookup_widget("histoHighEntry"));

		if (xpos < 0.0)
			xpos = 0.0;
		if (xpos > 1.0)
			xpos = 1.0;

		switch (_type_of_scale) {
		case SCALE_LOW:
			if ((float)xpos > _highlights) {
				_shadows = _highlights;
			} else {
				_shadows = (float)xpos;
			}
			buffer = g_strdup_printf("%.7f", _shadows);
			gtk_entry_set_text(histoShadEntry, buffer);
			break;

		case SCALE_MID:
			if (_highlights == _shadows) {
				_midtones = _highlights;
			} else {
				_midtones = ((float)xpos - _shadows) / (_highlights - _shadows);
			}
			if (_midtones > 1.f) _midtones = 1.f;
			if (_midtones < 0.f) _midtones = 0.f;
			buffer = g_strdup_printf("%.7f", _midtones);
			gtk_entry_set_text(histoMidEntry, buffer);
			break;

		case SCALE_HI:
			if ((float)xpos < _shadows) {
				_shadows =_highlights;
			} else {
				_highlights = (float)xpos;
			}
			buffer = g_strdup_printf("%.7f", _highlights);
			gtk_entry_set_text(histoHighEntry, buffer);
			break;
		}
		set_cursor("grabbing");
		update_histo_mtf();
		g_free(buffer);
	}
	return FALSE;
}

void on_drawingarea_histograms_leave_notify_event(GtkWidget *widget,
		GdkEvent *event, gpointer user_data) {
	set_cursor("default");
}

gboolean on_drawingarea_histograms_button_press_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	int width = get_width_of_histo();
	int height = get_height_of_histo();

	if (on_gradient((GdkEvent *) event, width, height)) {
		float delta = ((_highlights - _shadows) * _midtones) + _shadows;

		_click_on_histo = TRUE;
		gdouble xpos = ((GdkEventButton*) event)->x / (gdouble) width;
		if (fabsf((float)xpos - _highlights) < fabsf((float)xpos - _shadows) && fabsf((float)xpos - _highlights) < fabsf((float)xpos - delta)) {
			_type_of_scale = SCALE_HI;
		} else if (fabsf((float)xpos - _shadows) < fabsf((float)xpos - delta) && fabsf((float)xpos - _shadows) < fabsf((float)xpos - _highlights)) {
			_type_of_scale = SCALE_LOW;
		} else if ((_shadows == _highlights) && (_shadows > 0.f)) {
			_type_of_scale = SCALE_LOW;
		} else if ((_shadows == _highlights) && (_shadows == 0.f)) {
			_type_of_scale = SCALE_HI;
		} else {
			_type_of_scale = SCALE_MID;
		}
		set_cursor("grabbing");
	}

	return FALSE;
}

gboolean on_drawingarea_histograms_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	set_cursor("default");
	if (_click_on_histo) {
		_click_on_histo = FALSE;
		set_cursor_waiting(TRUE);
		update_histo_mtf();
		histo_recompute();
		set_cursor_waiting(FALSE);
	}
	return FALSE;
}

void on_histoMidEntry_activate(GtkEntry *entry, gpointer user_data) {
	float mid = atof(gtk_entry_get_text(entry));
	if (mid <= _shadows) mid = _shadows;
	if (mid >= _highlights) mid = _highlights;
	_midtones = mid;
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	histo_recompute();
	gchar *str = g_strdup_printf("%8.7f", mid);
	gtk_entry_set_text(entry, str);
	g_free(str);
	set_cursor_waiting(FALSE);
}

void on_histoShadEntry_activate(GtkEntry *entry, gpointer user_data) {
	float lo = atof(gtk_entry_get_text(entry));
	if (lo <= 0.f) lo = 0.f;
	if (lo >= _highlights) lo = _highlights;
	_shadows = lo;
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	histo_recompute();
	gchar *str = g_strdup_printf("%8.7f", lo);
	gtk_entry_set_text(entry, str);
	g_free(str);
	set_cursor_waiting(FALSE);
}

void on_histoHighEntry_activate(GtkEntry *entry, gpointer user_data) {
	float hi = atof(gtk_entry_get_text(entry));
	if (hi <= _shadows) hi = _shadows;
	if (hi >= 1.f) hi = 1.f;
	_highlights = hi;
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	histo_recompute();
	gchar *str = g_strdup_printf("%8.7f", hi);
	gtk_entry_set_text(entry, str);
	g_free(str);
	set_cursor_waiting(FALSE);
}

void apply_mtf_to_sequence(struct mtf_data *mtf_args) {
	struct generic_seq_args *args = create_default_seqargs(mtf_args->seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = mtf_args->seq->selnum;
	args->prepare_hook = seq_prepare_hook;
	args->finalize_hook = seq_finalize_hook;
	args->image_hook = mtf_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("Midtone Transfer Function");
	args->has_output = TRUE;
	args->new_seq_prefix = mtf_args->seqEntry;
	args->load_new_sequence = TRUE;
	args->user = mtf_args;

	mtf_args->fit = NULL;	// not used here

	start_in_new_thread(generic_sequence_worker, args);
}
