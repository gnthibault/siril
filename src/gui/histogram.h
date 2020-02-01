#ifndef _HIST_H_
#define _HIST_H_

typedef enum {
	SCALE_LOW,
	SCALE_MID,
	SCALE_HI
} ScaleType;

gsl_histogram* computeHisto(fits*, int);
gsl_histogram* computeHisto_Selection(fits*, int, rectangle *);
void compute_histo_for_gfit();
void invalidate_gfit_histogram();
void update_gfit_histogram_if_needed();
void clear_histograms();
float MTF(float x, float m, float lo, float hi);
float findMidtonesBalance(fits *fit, float *shadows, float *highlights);
void apply_histo_cancel();

void on_histoMidEntry_changed(GtkEditable *editable, gpointer user_data);
void on_histoShadEntry_changed(GtkEditable *editable, gpointer user_data);
void on_histoHighEntry_changed(GtkEditable *editable, gpointer user_data);

#endif
