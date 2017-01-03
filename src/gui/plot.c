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

#include "plot.h"

#include <cairo.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/siril.h"
#include "core/proto.h"
#include "gui/callbacks.h"
#include "kplot.h"
#include "algos/PSF.h"

static GtkWidget *drawingPlot = NULL, *sourceCombo = NULL, *combo = NULL,
		 *buttonExport = NULL, *buttonClearAll = NULL, *buttonClearLatest = NULL;
static pldata *plot_data;
static struct kpair ref;
static gboolean is_fwhm = FALSE, use_photometry = FALSE, requires_color_update = FALSE;
static char *ylabel = NULL;
static enum photmetry_source selected_source = ROUNDNESS;

static void update_ylabel();
static void set_colors(struct kplotcfg *cfg);
void on_GtkEntryCSV_changed(GtkEditable *editable, gpointer user_data);

static pldata *alloc_plot_data(int size) {
	pldata *plot = malloc(sizeof(pldata));
	if (!plot) return NULL;
	plot->data = malloc(size * sizeof(struct kpair));
	if (!plot->data) { free(plot); return NULL; }
	plot->nb = size;
	plot->next = NULL;
	return plot;
}

static void build_registration_dataset(sequence *seq, int layer, int ref_image, pldata *plot) {
	int i, j;

	for (i = 0, j = 0; i < plot->nb; i++) {
		if (!seq->imgparam[i].incl) continue;
		plot->data[j].x = (double) i;
		plot->data[j].y = is_fwhm ?
						seq->regparam[layer][i].fwhm :
						seq->regparam[layer][i].quality;
		j++;
	}
	plot->nb = j;

	ref.x = (double) ref_image;
	ref.y = is_fwhm ?
			seq->regparam[layer][ref_image].fwhm :
			seq->regparam[layer][ref_image].quality;

}

static void build_photometry_dataset(sequence *seq, int dataset, int size, int ref_image, pldata *plot) {
	int i, j;
	double offset = -1001.0;
	fitted_PSF **psfs = seq->photometry[dataset], *ref_psf;
	if (seq->reference_star >= 0 && !seq->photometry[seq->reference_star])
		seq->reference_star = -1;

	for (i = 0, j = 0; i < size; i++) {
		if (!seq->imgparam[i].incl) continue;
		if (psfs[i]) {
			plot->data[j].x = (double)i;
			switch (selected_source) {
				case ROUNDNESS:
					plot->data[j].y = psfs[i]->fwhmy / psfs[i]->fwhmx;
					break;
				case FWHM:
					plot->data[j].y = psfs[i]->fwhmx;
					break;
				case AMPLITUDE:
					plot->data[j].y = psfs[i]->A;
					break;
				case MAGNITUDE:
					plot->data[j].y = psfs[i]->mag;
					if (seq->reference_star >= 0) {
						/* we have a reference star for the sequence,
						 * with photometry data */
						ref_psf = seq->photometry[seq->reference_star][i];
						if (ref_psf)
							offset = seq->reference_mag - ref_psf->mag;
					}
					else if (com.magOffset > 0.0)
						offset = com.magOffset;

					/* apply the absolute apparent magnitude offset */
					if (offset > -1000.0)
						plot->data[j].y += offset;
					break;
				case BACKGROUND:
					plot->data[j].y = psfs[i]->B;
					break;
				case X_POSITION:
					plot->data[j].y = psfs[i]->xpos;
					break;
				case Y_POSITION:
					plot->data[j].y = psfs[i]->ypos;
					break;
			}
		}

		/* we'll just take the reference image point from the last data set rendered */
		if (i == ref_image) {
			ref.x = (double) ref_image;
			ref.y = plot->data[j].y;
		}

		j++;
	}
	plot->nb = j;
}

static int exportCSV(pldata *plot, sequence *seq) {
	int i, j, ret = 0;
	const gchar *file;
	gchar *filename, *msg;
	GtkEntry *EntryCSV;

	EntryCSV = GTK_ENTRY(lookup_widget("GtkEntryCSV"));
	file = gtk_entry_get_text(EntryCSV);
	if (file && file[0] != '\0') {
		filename = g_strndup(file, strlen(file) + 5);
		g_strlcat(filename, ".csv", strlen(file) + 5);
		FILE *csv = fopen(filename, "w");
		if (csv == NULL) {
			ret = 1;
		} else {
			if (use_photometry) {
				pldata *tmp_plot = plot;
					for (i = 0, j = 0; i < plot->nb; i++) {
						if (!seq->imgparam[i].incl)
							continue;
						int x = 0;
						fprintf(csv, "%g", tmp_plot->data[j].x);
						while (x < MAX_SEQPSF && seq->photometry[x]){
							fprintf(csv, ", %g", tmp_plot->data[j].y);
							tmp_plot = tmp_plot->next;
							++x;
						}
						fprintf(csv, "\n");
						tmp_plot = plot;
						j++;
					}
			} else {
				for (i = 0, j = 0; i < plot->nb; i++) {
					if (!seq->imgparam[i].incl)
						continue;
					fprintf(csv, "%g, %g\n", plot->data[j].x, plot->data[j].y);
					j++;
				}
			}
			g_free(filename);
			fclose(csv);
		}
	}
	if (!ret) {
		msg = siril_log_message(_("%s.csv has been saved.\n"), file);
		show_dialog(msg, _("Information"), "gtk-dialog-info");
	}
	else {
		show_dialog(_("Something went wrong while saving plot"), _("Error"), "gtk-dialog-error");
	}
	return 0;
}

static void free_plot_data() {
	pldata *plot = plot_data;
	while (plot) {
		pldata *next = plot->next;
		if (plot->data)
			free(plot->data);
		free(plot);
		plot = next;
	}
	plot_data = NULL;
}

void on_plotSourceCombo_changed(GtkComboBox *box, gpointer user_data) {
	use_photometry = gtk_combo_box_get_active(GTK_COMBO_BOX(box));
	gtk_widget_set_visible(combo, use_photometry);
	drawPlot();
}

void on_GtkEntryCSV_changed(GtkEditable *editable, gpointer user_data) {
	const gchar *txt;
	if (!buttonExport) return;
	txt = gtk_entry_get_text(GTK_ENTRY(editable));
	gtk_widget_set_sensitive(buttonExport, txt[0] != '\0' && plot_data);
}

void reset_plot() {
	free_plot_data();
	if (sourceCombo) {
		gtk_combo_box_set_active(GTK_COMBO_BOX(sourceCombo), 0);
		gtk_widget_set_visible(sourceCombo, FALSE);
		gtk_widget_set_visible(combo, FALSE);
		gtk_widget_set_sensitive(buttonExport, FALSE);
		gtk_widget_set_sensitive(buttonClearLatest, FALSE);
		gtk_widget_set_sensitive(buttonClearAll, FALSE);
	}
}

void drawPlot() {
	int i, ref_image, layer = 0;
	sequence *seq;

	if (drawingPlot == NULL) {
		drawingPlot = lookup_widget("DrawingPlot");
		combo = lookup_widget("plotCombo");
		sourceCombo = lookup_widget("plotSourceCombo");
		buttonExport = lookup_widget("ButtonSaveCSV");
		buttonClearAll  = lookup_widget("clearAllPhotometry");
		buttonClearLatest = lookup_widget("clearLastPhotometry");
	}

	seq = &com.seq;
	if (plot_data)
		free_plot_data();

	if (seq->reference_image == -1)
		ref_image = 0;
	else ref_image = seq->reference_image;

	if (use_photometry) {
		// photometry data display
		pldata *plot;
		update_ylabel();

		plot = alloc_plot_data(seq->number);
		plot_data = plot;
		for (i = 0; i < MAX_SEQPSF && seq->photometry[i]; i++) {
			if (i > 0) {
				plot->next = alloc_plot_data(seq->number);
				plot = plot->next;
			}
			
			build_photometry_dataset(seq, i, seq->number, ref_image, plot);
		}
	} else {
		// registration data display
		if (!(seq->regparam))
			return;

		for (i = 0; i < seq->nb_layers; i++) {
			if (com.seq.regparam[i]) {
				layer = i;
				break;
			}
		}
		if ((!seq->regparam[layer]))
			return;

		if (seq->regparam[layer][ref_image].fwhm > 0.0f) {
			is_fwhm = TRUE;
			ylabel = _("FWHM");
		} else if (seq->regparam[layer][ref_image].quality >= 0.0) {
			is_fwhm = FALSE;
			ylabel = _("Quality");
		} else return;

		/* building data array */
		plot_data = alloc_plot_data(seq->number);

		build_registration_dataset(seq, layer, ref_image, plot_data);
	}
	on_GtkEntryCSV_changed(GTK_EDITABLE(lookup_widget("GtkEntryCSV")), NULL);
	gtk_widget_queue_draw(drawingPlot);
}

void on_ButtonSaveCSV_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	exportCSV(plot_data, &com.seq);
	set_cursor_waiting(FALSE);
}

void free_photometry_set(sequence *seq, int set) {
	int j;
	for (j = 0; j < seq->number; j++) {
		if (seq->photometry[set][j])
			free(seq->photometry[set][j]);
	}
	free(seq->photometry[set]);
	seq->photometry[set] = NULL;
}

void on_clearLatestPhotometry_clicked(GtkButton *button, gpointer user_data) {
	int i;
	for (i = 0; i < MAX_SEQPSF && com.seq.photometry[i]; i++);
	if (i != 0) {
		i--;
		free_photometry_set(&com.seq, i);
	}
	if (i == 0)
		reset_plot();
	drawPlot();
}

void on_clearAllPhotometry_clicked(GtkButton *button, gpointer user_data) {
	int i;
	for (i = 0; i < MAX_SEQPSF && com.seq.photometry[i]; i++) {
		free_photometry_set(&com.seq, i);
	}
	reset_plot();
	drawPlot();
}

gboolean on_DrawingPlot_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
	guint width, height, i, j;
	double mean/*, sigma*/;
	int min, max, nb_graphs = 0;
	struct kpair *avg;
	struct kplotcfg	 cfgplot;
	struct kdatacfg	 cfgdata;
	struct kdata *d1, *ref_d, *mean_d;
	struct kplot *p;

	if (plot_data) {
		pldata *plot = plot_data;
		d1 = ref_d = mean_d = NULL;

		kplotcfg_defaults(&cfgplot);
		kdatacfg_defaults(&cfgdata);
		set_colors(&cfgplot);
		cfgplot.xaxislabel = _("Frames");
		cfgplot.yaxislabel = ylabel;
		cfgplot.yaxislabelrot = M_PI_2 * 3.0;
		//cfgplot.y2axislabel = _("Sigma");
		cfgplot.xticlabelpad = cfgplot.yticlabelpad = 10.0;
		cfgdata.point.radius = 10;

		p = kplot_alloc(&cfgplot);

		// data plots
		while (plot) {
			d1 = kdata_array_alloc(plot->data, plot->nb);
			kplot_attach_data(p, d1, ((plot_data->nb <= 100) ? KPLOT_LINESPOINTS : KPLOT_LINES), NULL);
			plot = plot->next;
			nb_graphs++;
		}

		/* mean and sigma */
		mean = kdata_ymean(d1);
		//sigma = kdata_ystddev(d1);
		min = plot_data->data[0].x;
		/* if reference is ploted, we take it as maximum if it is */
		max = (plot_data->data[plot_data->nb-1].x > ref.x) ? plot_data->data[plot_data->nb-1].x + 1: ref.x + 1;

		if (nb_graphs == 1) {
			avg = calloc(max - min, sizeof(struct kpair));
			j = min;
			for (i = 0; i < max - min; i++) {
				avg[i].x = (double) j;
				avg[i].y = mean;
				++j;
			}

			mean_d = kdata_array_alloc(avg, max - min);
			kplot_attach_data(p, mean_d, KPLOT_LINES, NULL);	// mean plot
			free(avg);
		}

		ref_d = kdata_array_alloc(&ref, 1);

		kplot_attach_data(p, ref_d, KPLOT_POINTS, &cfgdata);	// ref image dot

		width = gtk_widget_get_allocated_width(widget);
		height = gtk_widget_get_allocated_height(widget);

		cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
		cairo_rectangle(cr, 0.0, 0.0, width, height);
		cairo_fill(cr);
		kplot_draw(p, width, height, cr);

		/* copy graph colours for star highlight */
		if (requires_color_update) {
			for (i = 0; i < cfgplot.clrsz; i++) {
				com.seq.photometry_colors[i][0] = cfgplot.clrs[i].rgba[0];
				com.seq.photometry_colors[i][1] = cfgplot.clrs[i].rgba[1];
				com.seq.photometry_colors[i][2] = cfgplot.clrs[i].rgba[2];
			}
			redraw(com.cvport, REMAP_ONLY);
			requires_color_update = FALSE;
		}

		kplot_free(p);
		kdata_destroy(d1);
		kdata_destroy(ref_d);
		if (mean_d)
			kdata_destroy(mean_d);
	}
	return FALSE;
}

void on_plotCombo_changed(GtkComboBox *box, gpointer user_data) {
	drawPlot();
}

static void update_ylabel() {
	selected_source = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
	switch (selected_source) {
		case ROUNDNESS:
			ylabel = _("Star roundness (1 is round)");
			break;
		case FWHM:
			ylabel = _("FWHM");
			break;
		case AMPLITUDE:
			ylabel = _("Amplitude");
			break;
		case MAGNITUDE:
			if (com.magOffset > 0.0 || com.seq.reference_star >= 0)
				ylabel = _("Star magnitude (absolute)");
			else ylabel = _("Star magnitude (relative, use setmag)");
			break;
		case BACKGROUND:
			ylabel = _("Background value");
			break;
		case X_POSITION:
			ylabel = _("Star position on X axis");
			break;
		case Y_POSITION:
			ylabel = _("Star position on Y axis");
			break;
	}
}

void notify_new_photometry() {
	requires_color_update = TRUE;
	gtk_widget_set_visible(sourceCombo, TRUE);
	gtk_combo_box_set_active(GTK_COMBO_BOX(sourceCombo), 1);
	gtk_widget_set_sensitive(buttonClearLatest, TRUE);
	gtk_widget_set_sensitive(buttonClearAll, TRUE);
}

static void set_colors(struct kplotcfg *cfg) {
	int i;
	cfg->clrsz = MAX_SEQPSF;
	cfg->clrs = calloc(cfg->clrsz, sizeof(struct kplotccfg));
	for (i = 0; i < cfg->clrsz; i++) {
		cfg->clrs[i].type = KPLOTCTYPE_RGBA;
		cfg->clrs[i].rgba[3] = 1.0;
	}
	cfg->clrs[0].rgba[0] = 0x94 / 255.0;
	cfg->clrs[0].rgba[1] = 0x04 / 255.0;
	cfg->clrs[0].rgba[2] = 0xd3 / 255.0;
	cfg->clrs[1].rgba[0] = 0x00 / 255.0;
	cfg->clrs[1].rgba[1] = 0x9e / 255.0;
	cfg->clrs[1].rgba[2] = 0x73 / 255.0;
	cfg->clrs[2].rgba[0] = 0x56 / 255.0;
	cfg->clrs[2].rgba[1] = 0xb4 / 255.0;
	cfg->clrs[2].rgba[2] = 0xe9 / 255.0;
	cfg->clrs[3].rgba[0] = 0xe6 / 255.0;
	cfg->clrs[3].rgba[1] = 0x9f / 255.0;
	cfg->clrs[3].rgba[2] = 0x00 / 255.0;
	cfg->clrs[4].rgba[0] = 0xf0 / 255.0;
	cfg->clrs[4].rgba[1] = 0xe4 / 255.0;
	cfg->clrs[4].rgba[2] = 0x42 / 255.0;
	cfg->clrs[5].rgba[0] = 0x00 / 255.0;
	cfg->clrs[5].rgba[1] = 0x72 / 255.0;
	cfg->clrs[5].rgba[2] = 0xb2 / 255.0;
	cfg->clrs[6].rgba[0] = 0xe5 / 255.0;
	cfg->clrs[6].rgba[1] = 0x1e / 255.0;
	cfg->clrs[6].rgba[2] = 0x10 / 255.0;
}

