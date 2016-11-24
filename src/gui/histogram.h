#ifndef _HIST_H_
#define _HIST_H_

gsl_histogram* computeHisto(fits*, int);
gsl_histogram* computeHisto_Selection(fits*, int, rectangle *);
gsl_histogram* histo_bg(fits*, int, double);
void compute_histo_for_gfit(int force);
void update_gfit_histogram_if_needed();
void set_histogram(gsl_histogram *histo, int layer);
void clear_histograms();
void set_histo_toggles_names();
void erase_histo_display(cairo_t *cr, int width, int height);
void display_histo(gsl_histogram *histo, cairo_t *cr, int layer, int width,
		int height, double ZoomH, double zoomV);
void reset_curors_and_values();
void update_histo_mtf();
void apply_mtf_to_fits(fits *fit, double m, double lo, double hi);
void apply_mtf_to_histo(gsl_histogram *histo, double norm, double m, double lo,
		double hi);
double MTF(double x, double m);
double findMidtonesBalance(fits *fit, double *shadows, double *highlights) ;
void on_scale_midtones_value_changed(GtkRange *range, gpointer user_data);
void on_scale_shadows_value_changed(GtkRange *range, gpointer user_data);
void on_scale_highlights_value_changed(GtkRange *range, gpointer user_data);

#endif
