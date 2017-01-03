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

#include <gsl/gsl_histogram.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "core/siril.h"
#include "core/proto.h"
#include "io/single_image.h"
#include "gui/histogram.h"
#include "gui/callbacks.h"	// for lookup_widget()
#include "core/undo.h"

#define shadowsClipping -2.80 /* Shadows clipping point measured in sigma units from the main histogram peak. */
#define targetBackground 0.25 /* final "luminance" of the image for autostretch in the [0,1] range */
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
static double graph_height = 0.;
static gsl_histogram *histCpy[MAXVPORT] = { NULL, NULL, NULL, NULL };
static uint64_t clipped[] = { 0, 0 };

static GtkToggleToolButton *toggles[MAXVPORT] = { NULL };
static GtkToggleToolButton *toggleGrid = NULL;

static void get_sliders_values(double *m, double *lo, double *hi) {
	static GtkRange *scale_transfert_function[3] = { NULL, NULL, NULL };

	if (scale_transfert_function[1] == NULL) {
		scale_transfert_function[0] = GTK_RANGE(lookup_widget("scale_shadows"));
		scale_transfert_function[1] = GTK_RANGE(
				lookup_widget("scale_midtones"));
		scale_transfert_function[2] = GTK_RANGE(
				lookup_widget("scale_highlights"));
	}

	*m = gtk_range_get_value(scale_transfert_function[1]);
	*lo = gtk_range_get_value(scale_transfert_function[0]);
	*hi = gtk_range_get_value(scale_transfert_function[2]);
}

static void _init_clipped_pixels() {
	clipped[0] = 0;
	clipped[1] = 0;
}

static void _initialize_clip_text() {
	static GtkEntry *clip_high = NULL, *clip_low = NULL;

	if (clip_high == NULL) {
		clip_high = GTK_ENTRY(
				gtk_builder_get_object(builder, "clip_highlights"));
		clip_low = GTK_ENTRY(gtk_builder_get_object(builder, "clip_shadows"));
	}
	gtk_entry_set_text(clip_low, "0.000%");
	gtk_entry_set_text(clip_high, "0.000%");
}

static void _update_clipped_pixels(int data) {
	static GtkEntry *clip_high = NULL, *clip_low = NULL;
	double tmp;
	char buffer[16];

	if (clip_high == NULL) {
		clip_high = GTK_ENTRY(lookup_widget("clip_highlights"));
		clip_low = GTK_ENTRY(lookup_widget("clip_shadows"));
	}
	tmp = (double) clipped[1] * 100.0 / data;
	g_snprintf(buffer, sizeof(buffer), "%.3lf%%", tmp);
	gtk_entry_set_text(clip_high, buffer);
	tmp = (double) clipped[0] * 100.0 / data;
	g_snprintf(buffer, sizeof(buffer), "%.3lf%%", tmp);
	gtk_entry_set_text(clip_low, buffer);

}

static int is_histogram_visible() {
	GtkWidget *window = lookup_widget("histogram_window");
	return gtk_widget_get_visible(window);
}

gsl_histogram* computeHisto(fits* fit, int layer) {
	assert(layer < 3);
	size_t i, ndata, size;
	WORD *buf;

	size = (size_t) get_normalized_value(fit);
	gsl_histogram* histo = gsl_histogram_alloc(size + 1);
	gsl_histogram_set_ranges_uniform(histo, 0, size);

	buf = fit->pdata[layer];
	ndata = fit->rx * fit->ry;
//#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static) // cause errors !!!
	for (i = 0; i < ndata; i++) {
		gsl_histogram_increment(histo, (double) buf[i]);
	}
	return histo;
}

gsl_histogram* computeHisto_Selection(fits* fit, int layer,
		rectangle *selection) {
	assert(layer < 3);
	WORD *from;
	size_t stridefrom, i, j, size;

	size = (size_t) get_normalized_value(fit);
	gsl_histogram* histo = gsl_histogram_alloc(size + 1);
	gsl_histogram_set_ranges_uniform(histo, 0, size);

	from = fit->pdata[layer] + (fit->ry - selection->y - selection->h) * fit->rx
			+ selection->x;
	stridefrom = fit->rx - selection->w;
	for (i = 0; i < selection->h; i++) {
		for (j = 0; j < selection->w; j++) {
			gsl_histogram_increment(histo, (double) *from);
			from++;
		}
		from += stridefrom;
	}
	return histo;
}

//compute histogram for all pixel values below backgroud value
gsl_histogram* histo_bg(fits* fit, int layer, double bg) {
	unsigned int i, ndata;
	WORD *buf;
	if (layer >= 3)
		return NULL;
	size_t size = (size_t) bg;
	gsl_histogram* histo = gsl_histogram_alloc(size + 1);
	gsl_histogram_set_ranges_uniform(histo, 0, size);

	buf = fit->pdata[layer];
	ndata = fit->rx * fit->ry;
	for (i = 0; i < ndata; i++) {
		if (buf[i] <= (WORD) bg)
			gsl_histogram_increment(histo, (double) buf[i]);
	}
	return histo;
}

void compute_histo_for_gfit(int force) {
	int nb_layers = 3;
	int i;
	if (gfit.naxis == 2)
		nb_layers = 1;
	for (i = 0; i < nb_layers; i++) {
		if (force || !com.layers_hist[i])
			set_histogram(computeHisto(&gfit, i), i);
		if (histCpy[i])
			gsl_histogram_free(histCpy[i]);
		histCpy[i] = gsl_histogram_clone(com.layers_hist[i]);
	}
	set_histo_toggles_names();
}

void update_gfit_histogram_if_needed() {
	static GtkWidget *drawarea = NULL;
	if (!drawarea) {
		drawarea = lookup_widget("drawingarea_histograms");
	}
	if (is_histogram_visible()) {
		compute_histo_for_gfit(1);
		gtk_widget_queue_draw(drawarea);
	}
}

void set_histogram(gsl_histogram *histo, int layer) {
	assert(layer >= 0 && layer < MAXVPORT);
	if (com.layers_hist[layer])
		gsl_histogram_free(com.layers_hist[layer]);
	com.layers_hist[layer] = histo;
}

void clear_histograms() {
	int i;
	for (i = 0; i < MAXVPORT; i++) {
		set_histogram(NULL, i);
		if (histCpy[i])
			gsl_histogram_free(histCpy[i]);
		histCpy[i] = NULL;
	}
}

void init_toggles() {
	if (!toggles[0]) {
		toggles[0] = GTK_TOGGLE_TOOL_BUTTON(
				gtk_builder_get_object(builder, "histoToolRed"));
		toggles[1] = GTK_TOGGLE_TOOL_BUTTON(
				gtk_builder_get_object(builder, "histoToolGreen"));
		toggles[2] = GTK_TOGGLE_TOOL_BUTTON(
				gtk_builder_get_object(builder, "histoToolBlue"));
		toggles[3] = NULL;
	}
	toggleGrid = GTK_TOGGLE_TOOL_BUTTON(
			gtk_builder_get_object(builder, "histoToolGrid"));
}

// sets the channel names of the toggle buttons in the histogram window, based on
// the number of layers of gfit
void set_histo_toggles_names() {
	init_toggles();

	if (gfit.naxis == 2) {
		const char* test = gtk_tool_button_get_label(GTK_TOOL_BUTTON(toggles[0]));
		if (strcmp(test, "K")) {
			gtk_tool_button_set_label(GTK_TOOL_BUTTON(toggles[0]), "K");
			gtk_widget_set_tooltip_text(GTK_WIDGET(toggles[0]), "Gray channel");
		}
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
		const char* test = gtk_tool_button_get_label(GTK_TOOL_BUTTON(toggles[0]));
		if (strcmp(test, "R")) {
			gtk_tool_button_set_label(GTK_TOOL_BUTTON(toggles[0]), "R");
			gtk_widget_set_tooltip_text(GTK_WIDGET(toggles[0]), "Red channel");
		}
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
}

static double get_histoZoomValueH() {
	static GtkAdjustment *histoAdjZoomH = NULL;
	if (!histoAdjZoomH)
		histoAdjZoomH = GTK_ADJUSTMENT(
				gtk_builder_get_object(builder, "histoAdjZoomH"));

	return gtk_adjustment_get_value(histoAdjZoomH);
}

static double get_histoZoomValueV() {
	static GtkAdjustment *histoAdjZoomV = NULL;
	if (!histoAdjZoomV)
		histoAdjZoomV = GTK_ADJUSTMENT(
				gtk_builder_get_object(builder, "histoAdjZoomV"));

	return gtk_adjustment_get_value(histoAdjZoomV);
}

static void adjust_histogram_vport_size() {
	GtkWidget *drawarea, *vport;
	int targetW, targetH;
	double zoomH = get_histoZoomValueH();
	double zoomV = get_histoZoomValueV();

	drawarea = lookup_widget("drawingarea_histograms");
	vport = lookup_widget("viewport1");
	int cur_width = gtk_widget_get_allocated_width(vport);
	int cur_height = gtk_widget_get_allocated_height(vport);
	targetW = (int) (((double) cur_width) * zoomH);
	targetH = (int) (((double) cur_height) * zoomV);
	gtk_widget_set_size_request(drawarea, targetW, targetH);
#ifdef HISTO_DEBUG
	fprintf(stdout, "Histo vport size (%d, %d)\n", targetW, targetH);
#endif
}

gboolean redraw_histo(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int i, width, height;
	double zoomH, zoomV;

	init_toggles();
	width = gtk_widget_get_allocated_width(widget);
	height = gtk_widget_get_allocated_height(widget);
#ifdef HISTO_DEBUG
	fprintf(stdout, "histogram redraw\n");
	fprintf(stdout, "w = %d and h = %d\n", width, height);
#endif
	zoomH = get_histoZoomValueH();
	zoomV = get_histoZoomValueV();

	if (histCpy[0] == NULL) {
		for (i = 0; i < gfit.naxes[2]; i++)
			histCpy[i] = gsl_histogram_clone(com.layers_hist[i]);
	}

	if (height == 1)
		return FALSE;
	erase_histo_display(cr, width, height);
	graph_height = 0.0;
	for (i = 0; i < MAXVPORT; i++) {
		if (com.layers_hist[i]
				&& (!toggles[i] || gtk_toggle_tool_button_get_active(toggles[i])))
			display_histo(com.layers_hist[i], cr, i, width, height, zoomH, zoomV);
	}
	return FALSE;
}

void on_histo_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	static GtkWidget *drawarea = NULL;
	if (!drawarea) {
		drawarea = lookup_widget("drawingarea_histograms");
	}
	gtk_widget_queue_draw(drawarea);
}

#if 0
static void draw_MTF_curve(cairo_t *cr, int width, int height) {
	double m, lo, hi, y;

	get_sliders_values(&m, &lo, &hi);

	if (m != 0.0 && m != 1.0) {

		cairo_set_line_width(cr, 1.0);
		cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);

		y = MTF(0.5, m);
//		cairo_curve_to(cr, lo * width, height, (((hi - lo) / 2.0 + lo) / 0.5) * m * width,
//				((1 - y) * height), hi * width, 0.0);
		cairo_curve_to(cr, 0, height, 0.01 * width, 0.01 * height, hi * width, 0.0);
//		printf("test : %.9lf et %.9lf\n", 1 - y,
//				(((hi - lo) / 2.0 + lo) / 0.5) * m);

		cairo_stroke(cr);
	}
}
#endif

static void draw_grid(cairo_t *cr, int width, int height) {
	double dash_format[] = { 1.0, 1.0};

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
void erase_histo_display(cairo_t *cr, int width, int height) {
	gboolean drawGrid;
	// clear all with background color
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_rectangle(cr, 0, 0, width, height);
	cairo_fill(cr);

	// draw curve // disabled because bugged
#if 0
		draw_MTF_curve(cr, width, height);
#endif
	// draw grid
	drawGrid = gtk_toggle_tool_button_get_active(toggleGrid);
	if (drawGrid)
		draw_grid(cr, width, height);
}

static gboolean is_log_scale() {
	static GtkToggleButton *HistoCheckLogButton = NULL;

	if (HistoCheckLogButton == NULL)
		HistoCheckLogButton = GTK_TOGGLE_BUTTON(lookup_widget("HistoCheckLogButton"));
	return (gtk_toggle_button_get_active(HistoCheckLogButton));
}

void display_histo(gsl_histogram *histo, cairo_t *cr, int layer, int width,
		int height, double zoomH, double zoomV) {
	int current_bin;
	size_t norm = gsl_histogram_bins(histo) - 1;

	float vals_per_px = (float) (norm) / (float) width;	// size of a bin
	size_t i, nb_orig_bins = gsl_histogram_bins(histo);

	// We need to store the binned histogram in order to find the binned maximum
	static double *displayed_values = NULL;
	static int nb_bins_allocated = 0;
	/* we create a bin for each pixel in the displayed width.
	 * nb_bins_allocated is thus equal to the width of the image */
	if (nb_bins_allocated != width) {
		double *tmp;
		nb_bins_allocated = width;
		tmp = realloc(displayed_values, nb_bins_allocated * sizeof(double));
		if (!tmp) {
			if (displayed_values)
				free(displayed_values);
			fprintf(stderr, "Failed to reallocate histogram bins\n");
			return;
		}
		displayed_values = tmp;
		memset(displayed_values, 0, nb_bins_allocated);
	}
	assert(displayed_values);

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
		double bin_val = 0.0;
		while (i < nb_orig_bins
				&& (float) i / vals_per_px <= (float) current_bin + 0.5f) {
			bin_val += gsl_histogram_get(histo, i);
			i++;
		}
		if (is_log_scale()) {
			bin_val = (bin_val == 0) ? bin_val : log(bin_val);
		}
		displayed_values[current_bin] = bin_val;
		if (bin_val > graph_height)	// check for maximum
			graph_height = bin_val;
		current_bin++;
	} while (i < nb_orig_bins && current_bin < nb_bins_allocated);
	for (i = 0; i < nb_bins_allocated; i++) {
		double bin_height = height - height * displayed_values[i] / graph_height;
		cairo_line_to(cr, i, bin_height);
	}
	cairo_stroke(cr);
}

void on_histogram_window_show(GtkWidget *object, gpointer user_data) {
//	register_selection_update_callback(_histo_on_selection_changed);
	_initialize_clip_text();
}

void on_histogram_window_hide(GtkWidget *object, gpointer user_data) {
//	unregister_selection_update_callback(_histo_on_selection_changed);
}

void on_button_histo_close_clicked() {
	graph_height = 0.;
	gtk_widget_hide(lookup_widget("histogram_window"));
	reset_curors_and_values();
}

void reset_curors_and_values() {
	gtk_range_set_value(
			GTK_RANGE(gtk_builder_get_object(builder, "scale_midtones")), 0.5);
	gtk_range_set_value(
			GTK_RANGE(gtk_builder_get_object(builder, "scale_shadows")), 0.0);
	gtk_range_set_value(
			GTK_RANGE(gtk_builder_get_object(builder, "scale_highlights")),
			1.0);
	_init_clipped_pixels();
	_initialize_clip_text();
	update_gfit_histogram_if_needed();
}

void on_button_histo_reset_clicked() {
	set_cursor_waiting(TRUE);
	reset_curors_and_values();
	set_cursor_waiting(FALSE);
}

void apply_mtf_to_fits(fits *fit, double m, double lo, double hi) {
	double pente;
	int i, chan, nb_chan, ndata;
	WORD *buf[3] =
			{ fit->pdata[RLAYER], fit->pdata[GLAYER], fit->pdata[BLAYER] };
	WORD norm = get_normalized_value(fit);

	assert(fit->naxes[2] == 1 || fit->naxes[2] == 3);
	nb_chan = fit->naxes[2];
	ndata = fit->rx * fit->ry;

	pente = 1.0 / (hi - lo);

	undo_save_state("Processing: Histogram Transformation "
			"(mid=%.3lf, low=%.3lf, high=%.3lf)", m, lo, hi);

	for (chan = 0; chan < nb_chan; chan++) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)
#endif
		for (i = 0; i < ndata; i++) {
			double pxl = ((double) buf[chan][i] / (double) norm);
			pxl = (pxl - lo < 0.0) ? 0.0 : pxl - lo;
			pxl *= pente;
			buf[chan][i] = round_to_WORD(MTF(pxl, m) * (double) norm);
		}
	}
}

gboolean on_scale_key_release_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	set_cursor_waiting(TRUE);
	update_histo_mtf();
	set_cursor_waiting(FALSE);
	return FALSE;
}

void on_button_histo_apply_clicked(GtkButton *button, gpointer user_data) {
	double m, lo, hi;
	int i;

	get_sliders_values(&m, &lo, &hi);

	set_cursor_waiting(TRUE);
	apply_mtf_to_fits(&gfit, m, lo, hi);
	_init_clipped_pixels();
	update_gfit_histogram_if_needed();
	for (i = 0; i < gfit.naxes[2]; i++)
		gsl_histogram_memcpy(histCpy[i], com.layers_hist[i]);
	reset_curors_and_values();

	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();

	set_cursor_waiting(FALSE);
}

double MTF(double x, double m) {
	double out;

	if (m == 0.0)
		out = 0.0;
	else if (m == 0.5)
		out = x;
	else if (m == 1.0)
		out = 1.0;
	else {
		out = ((m - 1.0) * x) / (((2.0 * m - 1.0) * x) - m);
	}
	return out;
}

void apply_mtf_to_histo(gsl_histogram *histo, double norm, double m, double lo,
		double hi) {
	gsl_histogram *mtf_histo;
	unsigned short i;
	double pente = 1.0 / (hi - lo);

	mtf_histo = gsl_histogram_alloc((size_t) norm + 1);
	gsl_histogram_set_ranges_uniform(mtf_histo, 0, norm);

// #pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static) // disabled because of ISSUE #136 (https://free-astro.org/bugs/view.php?id=136)
	for (i = 0; i < round_to_WORD(norm); i++) {
		WORD mtf;
		double binval = gsl_histogram_get(histo, i);
		double pxl = ((double) i / norm);
		uint64_t clip[2] = { 0, 0 };

		if (i < round_to_WORD(lo * norm)) {
			pxl = lo;
			clip[0] += binval;
		} else if (i > round_to_WORD(hi * norm)) {
			pxl = hi;
			clip[1] += binval;
		}
		pxl -= lo;
		pxl *= pente;
		mtf = round_to_WORD(MTF(pxl, m) * norm);
		gsl_histogram_accumulate(mtf_histo, mtf, binval);
#ifdef _OPENMP
#pragma omp critical
#endif
		{
			clipped[0] += clip[0];
			clipped[1] += clip[1];
		}
	}
	gsl_histogram_memcpy(histo, mtf_histo);
	gsl_histogram_free(mtf_histo);
}

void update_histo_mtf() {
	double lo, hi, m;
	unsigned int i, data = 0;
	static GtkRange *scale_transfert_function[3] = { NULL, NULL, NULL };
	static GtkWidget *drawarea = NULL;
	double norm = (double) gsl_histogram_bins(com.layers_hist[0]) - 1;

	if (scale_transfert_function[0] == NULL) {
		scale_transfert_function[0] = GTK_RANGE(lookup_widget("scale_shadows"));
		scale_transfert_function[1] = GTK_RANGE(
				lookup_widget("scale_midtones"));
		scale_transfert_function[2] = GTK_RANGE(
				lookup_widget("scale_highlights"));
	}

	if (drawarea == NULL) {
		drawarea = lookup_widget("drawingarea_histograms");
	}
	get_sliders_values(&m, &lo, &hi);

	gtk_range_set_range(scale_transfert_function[0], 0.0, hi);
	gtk_range_set_range(scale_transfert_function[2], lo, 1.0);

	_init_clipped_pixels();
	for (i = 0; i < gfit.naxes[2]; i++) {
		gsl_histogram_memcpy(com.layers_hist[i], histCpy[i]);
		apply_mtf_to_histo(com.layers_hist[i], norm, m, lo, hi);
	}
	data = gfit.rx * gfit.ry * gfit.naxes[2];
	_update_clipped_pixels(data);

	/* redraw the histogram */
	gtk_widget_queue_draw(drawarea);
}

double findMidtonesBalance(fits *fit, double *shadows, double *highlights) {
	double c0 = 0.0, c1 = 0.0;
	double m = 0.0;
	int i, n, invertedChannels = 0;
	imstats *stat[3];

	n = fit->naxes[2];

	for (i = 0; i < n; ++i) {
		stat[i] = statistics(fit, i, NULL, STATS_BASIC | STATS_MAD, STATS_ZERO_NULLCHECK);
		if (!stat[i]) {
			siril_log_message(_("Error: no data computed.\n"));
			return 0.0;
		}

		if (stat[i]->median / stat[i]->normValue > 0.5)
			++invertedChannels;
	}

	if (invertedChannels < n) {
		for (i = 0; i < n; ++i) {
			double median, mad, normValue;

			normValue = stat[i]->normValue;
			median = stat[i]->median / normValue;
			mad = stat[i]->mad / normValue * MAD_NORM;

			c0 += median + shadowsClipping * mad;
			m += median;
		}
		c0 /= n;
		double m2 = m / n - c0;
		m = MTF(m2, targetBackground);
		*shadows = c0;
		*highlights = 1.0;
	} else {
		for (i = 0; i < n; ++i) {
			double median, mad, normValue;

			normValue = stat[i]->normValue;
			median = stat[i]->median / normValue;
			mad = stat[i]->mad / normValue * MAD_NORM;

			m += median;
			c1 += median - shadowsClipping * mad;
		}
		c1 /= n;
		double m2 = c1 - m / n;
		m = 1.0 - MTF(m2, targetBackground);
		*shadows = 0.0;
		*highlights = c1;

	}
	for (i = 0; i < n; ++i)
		free(stat[i]);
	return m;
}

void on_histoMidEntry_changed(GtkEditable *editable, gpointer user_data);
void on_histoShadEntry_changed(GtkEditable *editable, gpointer user_data);
void on_histoHighEntry_changed(GtkEditable *editable, gpointer user_data);

void on_scale_midtones_value_changed(GtkRange *range, gpointer user_data) {
	static GtkEntry *histoMidEntry = NULL;
	char buffer[10];
	double value;

	value = gtk_range_get_value(range);

	if (histoMidEntry == NULL)
		histoMidEntry = GTK_ENTRY(
				gtk_builder_get_object(builder, "histoMidEntry"));

	g_snprintf(buffer, 10, "%.7lf", value);
	g_signal_handlers_block_by_func(histoMidEntry, on_histoMidEntry_changed,
			NULL);
	gtk_entry_set_text(histoMidEntry, buffer);
	g_signal_handlers_unblock_by_func(histoMidEntry, on_histoMidEntry_changed,
			NULL);

	update_histo_mtf();
}

void on_scale_shadows_value_changed(GtkRange *range, gpointer user_data) {
	static GtkEntry *histoShadEntry = NULL;
	char buffer[10];
	double value;

	value = gtk_range_get_value(range);

	if (histoShadEntry == NULL)
		histoShadEntry = GTK_ENTRY(
				gtk_builder_get_object(builder, "histoShadEntry"));

	g_snprintf(buffer, 10, "%.7lf", value);

	g_signal_handlers_block_by_func(histoShadEntry, on_histoShadEntry_changed, NULL);
	gtk_entry_set_text(histoShadEntry, buffer);
	g_signal_handlers_unblock_by_func(histoShadEntry, on_histoShadEntry_changed, NULL);

	update_histo_mtf();
}

void on_scale_highlights_value_changed(GtkRange *range, gpointer user_data) {
	static GtkEntry *histoHighEntry = NULL;
	char buffer[10];
	double value;

	value = gtk_range_get_value(range);

	if (histoHighEntry == NULL)
		histoHighEntry = GTK_ENTRY(
				gtk_builder_get_object(builder, "histoHighEntry"));

	g_snprintf(buffer, 10, "%.7lf", value);

	g_signal_handlers_block_by_func(histoHighEntry, on_histoHighEntry_changed, NULL);
	gtk_entry_set_text(histoHighEntry, buffer);
	g_signal_handlers_unblock_by_func(histoHighEntry, on_histoHighEntry_changed, NULL);

	update_histo_mtf();
}

void on_histoMidEntry_changed(GtkEditable *editable, gpointer user_data) {
	double value;
	GtkRange *MidRange;

	value = atof(gtk_entry_get_text(GTK_ENTRY(editable)));
	if (value < 0.0) {
		gtk_entry_set_text(GTK_ENTRY(editable), "0.0000000");
		return;
	}

	if (value > 1.0) {
		gtk_entry_set_text(GTK_ENTRY(editable), "1.0000000");
		value = 0;
	}

	MidRange = GTK_RANGE(GTK_SCALE(gtk_builder_get_object(builder, "scale_midtones")));

	g_signal_handlers_block_by_func(MidRange, on_scale_midtones_value_changed, NULL);
	gtk_range_set_value(MidRange, value);
	g_signal_handlers_unblock_by_func(MidRange, on_scale_midtones_value_changed, NULL);
	update_histo_mtf();
}

void on_histoShadEntry_changed(GtkEditable *editable, gpointer user_data) {
	double value;
	GtkAdjustment *Shadows;
	GtkRange *ShadRange;

	value = atof(gtk_entry_get_text(GTK_ENTRY(editable)));
	if (value < 0.0) {
		gtk_entry_set_text(GTK_ENTRY(editable), "0.0000000");
		return;
	}

	if (value > 1.0) {
		gtk_entry_set_text(GTK_ENTRY(editable), "1.0000000");
		value = 0;
	}

	set_cursor_waiting(TRUE);
	Shadows = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adj_shadows"));
	ShadRange = GTK_RANGE(GTK_SCALE(gtk_builder_get_object(builder, "scale_shadows")));

	g_signal_handlers_block_by_func(ShadRange, on_scale_shadows_value_changed, NULL);
	gtk_adjustment_set_value(Shadows, value);
	g_signal_handlers_unblock_by_func(ShadRange, on_scale_shadows_value_changed, NULL);
	update_histo_mtf();
	set_cursor_waiting(FALSE);
}

void on_histoHighEntry_changed(GtkEditable *editable, gpointer user_data) {
	double value;
	GtkAdjustment *Highligts;
	GtkRange *HighRange;

	value = atof(gtk_entry_get_text(GTK_ENTRY(editable)));
	if (value < 0.0) {
		gtk_entry_set_text(GTK_ENTRY(editable), "0.0000000");
		return;
	}

	if (value > 1.0) {
		gtk_entry_set_text(GTK_ENTRY(editable), "1.0000000");
		value = 0;
	}

	set_cursor_waiting(TRUE);

	Highligts = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adj_highlights"));
	HighRange = GTK_RANGE(GTK_SCALE(gtk_builder_get_object(builder, "scale_highlights")));

	g_signal_handlers_block_by_func(HighRange, on_scale_highlights_value_changed, NULL);
	gtk_adjustment_set_value(Highligts, value);
	g_signal_handlers_unblock_by_func(HighRange, on_scale_highlights_value_changed, NULL);
	update_histo_mtf();
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
	static GtkWidget *drawarea = NULL;

	if (!drawarea) {
		drawarea = lookup_widget("drawingarea_histograms");
	}
	adjust_histogram_vport_size();
	gtk_widget_queue_draw(drawarea);
}

void on_histoToolAutoStretch_clicked(GtkToolButton *button, gpointer user_data) {
	double m, shadows = 0.0, highlights = 0.0;
	int i;
	unsigned int data;
	static GtkWidget *drawarea = NULL;
	double norm = (double) gsl_histogram_bins(com.layers_hist[0]) - 1;
	static GtkRange *scale_transfert_function[3] = { NULL, NULL, NULL };

	if (scale_transfert_function[1] == NULL) {
		scale_transfert_function[0] = GTK_RANGE(lookup_widget("scale_shadows"));
		scale_transfert_function[1] = GTK_RANGE(
				lookup_widget("scale_midtones"));
		scale_transfert_function[2] = GTK_RANGE(
				lookup_widget("scale_highlights"));
	}

	if (drawarea == NULL) {
		drawarea = lookup_widget("drawingarea_histograms");
	}

	set_cursor_waiting(TRUE);

	m = findMidtonesBalance(&gfit, &shadows, &highlights);
	_init_clipped_pixels();
	for (i = 0; i < gfit.naxes[2]; i++) {
		gsl_histogram_memcpy(com.layers_hist[i], histCpy[i]);
		apply_mtf_to_histo(com.layers_hist[i], norm, m, shadows, 1.0);
	}
	data = gfit.rx * gfit.ry * gfit.naxes[2];
	gtk_range_set_value(scale_transfert_function[1], m);
	gtk_range_set_value(scale_transfert_function[0], shadows);
	gtk_range_set_value(scale_transfert_function[2], highlights);
	_update_clipped_pixels(data);

	/* redraw the histogram */
	gtk_widget_queue_draw(drawarea);

	set_cursor_waiting(FALSE);
}
