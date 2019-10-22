#ifndef PROTO_H_
#define PROTO_H_
#include "core/siril.h"
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#include <fftw3.h>
#include <gsl/gsl_histogram.h>
#include <stdint.h>

#ifdef HAVE_LIBTIFF
#define uint64 uint64_hack_
#define int64 int64_hack_
#include <tiffio.h>
#undef uint64
#undef int64
#endif

/****************** image_format_fits.h ******************/
int	readfits(const char *filename, fits *fit, char *realname);
double get_exposure_from_fitsfile(fitsfile *fptr);
int import_metadata_from_fitsfile(fitsfile *fptr, fits *to);
void	clearfits(fits *);
int	readfits_partial(const char *filename, int layer, fits *fit, const rectangle *area, gboolean read_date);
int	read_opened_fits_partial(sequence *seq, int layer, int index, WORD *buffer, const rectangle *area);
int 	savefits(const char *, fits *);
int	copyfits(fits *from, fits *to, unsigned char oper, int layer);
int	copy_fits_metadata(fits *from, fits *to);
int	save1fits16(const char *filename, fits *fit, int layer);
int siril_fits_open_diskfile(fitsfile **fptr, const char *filename, int iomode, int *status);

void	rgb24bit_to_fits48bit(unsigned char *rgbbuf, fits *fit, gboolean inverted);
void	rgb8bit_to_fits16bit(unsigned char *graybuf, fits *fit);
void	rgb48bit_to_fits48bit(WORD *rgbbuf, fits *fit, gboolean inverted, gboolean change_endian);

void	fits_flip_top_to_bottom(fits *fit);
void	extract_region_from_fits(fits *from, int layer, fits *to, const rectangle *area);
int 	new_fit_image(fits **fit, int width, int height, int nblayer);
void	keep_first_channel_from_fits(fits *fit);

/****************** image_formats_internal.h ******************/
/* BMP */
int 	readbmp(const char *, fits *);
int 	savebmp(const char *, fits *);

/* PNM */
int 	import_pnm_to_fits(const char *filename, fits *fit);
int	saveNetPBM(const char *name, fits *fit);

/* PIC */
struct pic_struct {
	unsigned long magic;
	unsigned short width;
	unsigned short height;
	unsigned short bin[6];
	unsigned short nbplane;
	unsigned short hi;
	unsigned short lo;
	char *date;
	char *time;

	// internal stuff
	FILE *file;
};
int	readpic(const char *name, fits *fit);

/****************** image_formats_libraries.h ******************/
#ifdef HAVE_LIBTIFF
int readtif(const char *name, fits *fit);
int savetif(const char *name, fits *fit, uint16 bitspersample);
#endif

#ifdef HAVE_LIBJPEG
int readjpg(const char*, fits *);
int savejpg(const char *, fits *, int);
#endif

#ifdef HAVE_LIBPNG
int readpng(const char*, fits *);
int savepng(const char *filename, fits *fit, uint32_t bytes_per_sample,
		gboolean is_colour);
#endif

#ifdef HAVE_LIBRAW
int open_raw_files(const char *, fits *, int);
#endif

/****************** utils.h ******************/
int	round_to_int(double x);
int	roundf_to_int(float x);
WORD	round_to_WORD(double x);
BYTE	round_to_BYTE(double x);
BYTE	conv_to_BYTE(double x);
int	truncate_to_int32(uint64_t x);
gboolean isrgb(fits *fit);
char *f2utf8(const char *filename);
gboolean ends_with(const char *str, const char *ending);
int	get_extension_index(const char *filename);
int	is_readable_file(const char *filename);
int	stat_file(const char *filename2, image_type *type, char **realname);
const char *get_filename_ext(const char *filename);

gchar *siril_get_startup_dir();
int	changedir(const char *dir, gchar **err);
gchar *get_locale_filename(const gchar *path);
int	update_sequences_list(const char *sequence_name_to_select);
void	expand_home_in_filename(char *filename, int size);
WORD	get_normalized_value(fits*);
void	read_and_show_textfile(char*, char*);
void	swap_param(double *, double *);
char*	remove_ext_from_filename(const char *basename);
char*	str_append(char** data, const char* newdata);
char*	format_basename(char *root);
float	computePente(WORD *lo, WORD *hi);
void	load_css_style_sheet();
double	encodeJD(dateTime dt);

/**************** OS_utils.h *****************/
gchar *get_siril_locale_dir();
void	update_used_memory();
gchar *pretty_print_memory(int64_t bytes);
int test_available_space(int64_t req_size);
int	get_available_memory_in_MB();
int	get_max_memory_in_MB();
#ifdef _WIN32
gchar *get_special_folder(int csidl);
#endif
SirilWidget *siril_file_chooser_open(GtkWindow *parent, GtkFileChooserAction action);
SirilWidget *siril_file_chooser_add(GtkWindow *parent, GtkFileChooserAction action);
SirilWidget *siril_file_chooser_save(GtkWindow *parent, GtkFileChooserAction action);
gint siril_dialog_run(SirilWidget *widgetdialog);
void siril_widget_destroy(SirilWidget *widgetdialog);
gboolean allow_to_open_files(int nb_frames, int *nb_allowed_file);

/****************** quantize.h ***************/
int fits_img_stats_ushort(WORD *array, long nx, long ny, int nullcheck,
		WORD nullvalue, long *ngoodpix, WORD *minvalue, WORD *maxvalue,
		double *mean, double *sigma, double *noise1, double *noise2, double *noise3,
		double *noise5, int *status);

/****************** siril.h ******************/
/* crop sequence data from GUI */
struct crop_sequence_data {
	sequence *seq;
	rectangle area;
	const char *prefix;
	int retvalue;
};

/* Noise data from GUI */
struct noise_data {
	gboolean verbose;
	gboolean use_idle;
	fits *fit;
	double bgnoise[3];
	struct timeval t_start;
	int retval;
};

int 	threshlo(fits *fit, int level);
int 	threshhi(fits *fit, int level);
int 	nozero(fits *fit, int level);
int	soper(fits *a, double scalar, char oper);
int	imoper(fits *a, fits *b, char oper);
int 	addmax(fits *a, fits *b);
int	siril_fdiv(fits *a, fits *b, float scalar);
int siril_ndiv(fits *a, fits *b);
double 	gaussienne(double sigma, int size, double *gauss);
int 	unsharp(fits *,double sigma, double mult, gboolean verbose);
int 	shift(int sx, int sy);
double entropy(fits *fit, int layer, rectangle *area, imstats *opt_stats);
int 	loglut(fits *fit);
int	ddp(fits *a, int lev, float coef, float sig);
int	visu(fits *fit, int low, int high);
int	fill(fits *fit, int level, rectangle *arearg);
int 	off(fits *a, int level);
int	lrgb(fits *l, fits *r, fits *g, fits *b, fits *lrgb);
double	background(fits* fit, int reqlayer, rectangle *selection);
void	show_FITS_header(fits *);
double gauss_cvf(double p);
int get_wavelet_layers(fits *fit, int Nbr_Plan, int Plan, int Type, int reqlayer);
int extract_plans(fits *fit, int Nbr_Plan, int Type);
void compute_grey_flat(fits *fit);

/****************** seqfile.h ******************/
sequence * readseqfile(const char *name);
int	writeseqfile(sequence *seq);
gboolean existseq(const char *name);
int	buildseqfile(sequence *seq, int force_recompute);

/****************** registration_preview.h ******************/
void	redraw_previews();
void	set_preview_area(int preview_area, int centerX, int centerY);
void	init_mouse();
void	adjust_reginfo();
void	on_spinbut_shift_value_change(GtkSpinButton *spinbutton, gpointer user_data);

/****************** sequence_list.h ******************/
void	sequence_list_change_selection(gchar *path, gboolean new_value);
void	sequence_list_change_selection_index(int index);
void	sequence_list_change_current();
void	sequence_list_change_reference();
void	fill_sequence_list(sequence *seq, int layer, gboolean as_idle);
void	clear_sequence_list();

/****************** statistics_list.h ******************/
void computeStat();

#endif
