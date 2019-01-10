/* functions and variables defined only for linking */

#include "../core/siril.h"
#include "../core/pipe.h"

/* the global variables of the whole project (replacing main.c) */
cominfo com;	// the main data struct
fits gfit;	// currently loaded image
GtkBuilder *builder;	// get widget references anywhere
char **supported_extensions;

gboolean sequence_is_loaded() {
        fprintf(stderr, "ERROR: calling undefined function sequence_is_loaded\n");
	return FALSE;
}

void free_sequence(sequence *seq, gboolean free_seq_too) {
        fprintf(stderr, "ERROR: calling undefined function free_sequence\n");
}

void close_sequence(int loading_another) {
        fprintf(stderr, "ERROR: calling undefined function close_sequence\n");
}

sequence * readseqfile(const char *name){
        fprintf(stderr, "ERROR: calling undefined function readseqfile\n");
	return NULL;
}

#ifndef HAVE_STATS
imstats* free_stats(imstats *stat) {
        fprintf(stderr, "ERROR: calling undefined function free_stats\n");
	return NULL;
}

void add_stats_to_fit(fits *fit, int layer, imstats *stat) {
        fprintf(stderr, "ERROR: calling undefined function add_stats_to_fit\n");
}
#else
int fits_img_stats_ushort(WORD *array, long nx, long ny, int nullcheck, WORD
		nullvalue, long *ngoodpix, WORD *minvalue, WORD *maxvalue,
		double *mean, double *sigma, double *noise1, double *noise2,
		double *noise3, double *noise5, int *status)
{
        fprintf(stderr, "ERROR: calling undefined function pipe_send_message\n");
	return -1;
}
#endif

int pipe_send_message(pipe_message msgtype, pipe_verb verb, const char *arg) {
        fprintf(stderr, "ERROR: calling undefined function pipe_send_message\n");
	return 0;
}

double fit_get_min(fits *fit, int layer) {
        fprintf(stderr, "ERROR: calling undefined function fit_get_min\n");
	return 0.0;
}

double fit_get_max(fits *fit, int layer) {
        fprintf(stderr, "ERROR: calling undefined function fit_get_max\n");
	return 1.0;
}

int image_find_minmax(fits *fit) {
        fprintf(stderr, "ERROR: calling undefined function image_find_minmax\n");
	return 0;
}

image_type get_type_for_extension(const char *extension) {
        fprintf(stderr, "ERROR: calling undefined function get_type_for_extension\n");
	return TYPEFITS;
}

int single_image_is_loaded() {
        fprintf(stderr, "ERROR: calling undefined function single_image_is_loaded\n");
	return 0;
}

void set_GUI_MEM(unsigned long size) {
        fprintf(stderr, "ERROR: calling undefined function set_GUI_MEM\n");
}

void set_GUI_DiskSpace(double mem) {
        fprintf(stderr, "ERROR: calling undefined function set_GUI_DiskSpace\n");
}

void set_GUI_CWD() {
        fprintf(stderr, "ERROR: calling undefined function set_GUI_CWD\n");
}
