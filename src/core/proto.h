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
char*	list_header(fits *fit);
void	clearfits(fits *);
void	report_fits_error(int status);
int	readfits_partial(const char *filename, int layer, fits *fit, const rectangle *area);
int	read_opened_fits_partial(sequence *seq, int layer, int index, WORD *buffer, const rectangle *area);
int 	savefits(const char *, fits *);
void 	save_fits_header(fits *);
int	copyfits(fits *from, fits *to, unsigned char oper, int layer);
int copy_header(fits *from, fits *to);
int	save1fits16(const char *filename, fits *fit, int layer);

void	rgb24bit_to_fits48bit(unsigned char *rgbbuf, fits *fit, gboolean inverted);
void	rgb8bit_to_fits16bit(unsigned char *graybuf, fits *fit);
void	rgb48bit_to_fits48bit(WORD *rgbbuf, fits *fit, gboolean inverted, gboolean change_endian);

void	fits_flip_top_to_bottom(fits *fit);
void	extract_region_from_fits(fits *from, int layer, fits *to, const rectangle *area);
int 	new_fit_image(fits *fit, int width, int height, int nblayer);
void	keep_first_channel_from_fits(fits *fit);

/****************** image_formats_internal.h ******************/
/* BMP */
int 	readbmp(const char *, fits *);
int 	savebmp(const char *, fits *);
int	bmp32tofits48(unsigned char *rvb, int rx, int ry, fits *fitr, gboolean inverted);
int	bmp24tofits48(unsigned char *rvb, int rx, int ry, fits *fitr);
int	bmp8tofits(unsigned char *rvb, int rx, int ry, fits *fitr);

/* PNM */
int 	import_pnm_to_fits(const char *filename, fits *fit);
int	saveppm(const char *name, fits *fit);
int	savepgm(const char *name, fits *fit);

/* PIC */
struct pic_struct {
	unsigned short magic[2];
	unsigned short width;
	unsigned short height;
	unsigned short bin[6];
	unsigned short nbplane;
	unsigned short hi;
	unsigned short lo;
	char *date;			
	char *time;		

	// internal stuff
	int fd;
};
int	pictofit(WORD *buf, fits *fit);
int	pictofitrgb(WORD *buf, fits *fit);
int	readpic(const char *name, fits *fit);

/****************** image_formats_libraries.h ******************/
#ifdef HAVE_LIBTIFF
int	readtif8bits(TIFF* tif, uint32 width, uint32 height, uint16 nsamples, WORD **data);
int	readtifstrip(TIFF* tif, uint32 width, uint32 height, uint16 nsamples, WORD **data);
int 	readtif(const char *name, fits *fit);
int	savetif(const char *name, fits *fit, uint16 bitspersample);
#endif

#ifdef HAVE_LIBJPEG
int	readjpg(const char* , fits *);
int	savejpg(char *, fits *, int);
#endif

#ifdef HAVE_LIBPNG
int	readpng(const char* , fits *);
#endif

#ifdef HAVE_LIBRAW
int readraw(const char *, fits *);
int readraw_in_cfa(const char *, fits *);
int open_raw_files(const char *, fits *, int);
#endif

/****************** utils.h ******************/
int	round_to_int(double x);
WORD	round_to_WORD(double x);
BYTE	round_to_BYTE(double x);
BYTE	conv_to_BYTE(double x);
gboolean isrgb(fits *fit);
gboolean ends_with(const char *str, const char *ending);
int	get_extension_index(const char *filename);
int	is_readable_file(const char *filename);
int	stat_file(const char *filename2, image_type *type, char *realname);
const char *get_filename_ext(const char *filename);

int	changedir(const char *dir);
int	update_sequences_list(const char *sequence_name_to_select);
void	update_used_memory();
int	get_available_memory_in_MB();
void	expand_home_in_filename(char *filename, int size);
WORD	get_normalized_value(fits*);
void	read_and_show_textfile(char*, char*);
void	swap_param(double *, double *);
void	quicksort_d (double *a, int n);
void	quicksort_s (WORD *a, int n);
char*	remove_ext_from_filename(const char *basename);
char*	str_append(char** data, const char* newdata);
char*	format_basename(char *root);
float	computePente(WORD *lo, WORD *hi);
void load_css_style_sheet (char *path);

/****************** quantize.h ***************/
int fits_img_stats_ushort(WORD *array, long nx, long ny, int nullcheck,
		WORD nullvalue, long *ngoodpix, WORD *minvalue, WORD *maxvalue,
		double *mean, double *sigma, double *noise1, double *noise2, double *noise3,
		double *noise5, int *status);

/****************** siril.h ******************/
/* crop sequence data from GUI */
struct crop_sequence_data {
	sequence *seq;
	rectangle *area;
	const char *prefix;
	int retvalue;
};

/* median filter data from GUI */
struct median_filter_data {
	fits *fit;
	int ksize;
	double amount;
	int iterations;
};

/* Banding data from GUI */
struct banding_data {
	fits *fit;
	double sigma;
	double amount;
	gboolean protect_highlights;
	gboolean applyRotation;
	const gchar *seqEntry;
};

/* Noise data from GUI */
struct noise_data {
	gboolean verbose;
	fits *fit;
	double bgnoise[3];
	struct timeval t_start;
};

int 	threshlo(fits *fit, int level);
int 	threshhi(fits *fit, int level);
int 	nozero(fits *fit, int level);
int	soper(fits *a, double scalar, char oper);
int	imoper(fits *a, fits *b, char oper);
int sub_background(fits* image, fits* background, int layer);
int 	addmax(fits *a, fits *b);
int	fdiv(fits *a, fits *b, float scalar);
int ndiv(fits *a, fits *b);
int fmul(fits *a, float coeff);
double 	gaussienne(double sigma, int size, double *gauss);
int 	unsharp(fits *,double sigma, double mult, gboolean verbose);
int	crop(fits *fit, rectangle *bounds);
int 	shift(int sx, int sy);
double entropy(fits *fit, int layer, rectangle *area, imstats *opt_stats);
double contrast(fits* fit, int layer) ;
int 	loglut(fits *fit, int dir);
int	ddp(fits *a, int lev, float coef, float sig);
int	visu(fits *fit, int low, int high);
int	fill(fits *fit, int level, rectangle *arearg);
int 	off(fits *a, int level);
void 	mirrorx(fits *fit, gboolean verbose);
void 	mirrory(fits *fit, gboolean verbose);
void 	fits_rotate_pi(fits *fit);
int	lrgb(fits *l, fits *r, fits *g, fits *b, fits *lrgb);
gpointer seqpreprocess(gpointer empty);
void	initialize_preprocessing();
double	background(fits* fit, int reqlayer, rectangle *selection);
int backgroundnoise(fits* fit, double sigma[]);
imstats* statistics(fits *, int, rectangle *, int, int);
void	show_FITS_header(fits *);
#ifdef HAVE_OPENCV
int	verbose_resize_gaussian(fits *, int, int, int);
int	verbose_rotate_image(fits *, double, int, int);
#endif
double gauss_cvf(double p);
int get_wavelet_layers(fits *fit, int Nbr_Plan, int Plan, int Type, int reqlayer);
gpointer median_filter(gpointer p);
void apply_banding_to_sequence(struct banding_data *banding_args);
gpointer BandingEngineThreaded(gpointer p);
int BandingEngine(fits *fit, double sigma, double amount, gboolean protect_highlights, gboolean applyRotation);
gpointer noise(gpointer p);

/****************** sequence.h ******************/
int	read_single_sequence(char *realname, int imagetype);
int	seqsetnum(int image_number);
int	check_seq(int force);
int	check_only_one_film_seq(char* name);
int	set_seq(const char *);
char *	seq_get_image_filename(sequence *seq, int index, char *name_buf);
int	seq_read_frame(sequence *seq, int index, fits *dest);
int	seq_read_frame_part(sequence *seq, int layer, int index, fits *dest, const rectangle *area);
int	seq_load_image(sequence *seq, int index, fits *dest, gboolean load_it);
int	seq_open_image(sequence *seq, int index);
void	seq_close_image(sequence *seq, int index);
int	seq_opened_read_region(sequence *seq, int layer, int index, WORD *buffer, const rectangle *area);
void	set_fwhm_star_as_star_list(sequence *seq);
void	set_fwhm_star_as_star_list_with_layer(sequence *seq, int layer);
char *	fit_sequence_get_image_filename(sequence *seq, int index, char *name_buffer, gboolean add_fits_ext);
char *	get_possible_image_filename(sequence *seq, int image_number, char *name_buffer);
int	get_index_and_basename(const char *filename, char **basename, int *index, int *fixed);
void	initialize_sequence(sequence *seq, gboolean is_zeroed);
void	free_sequence(sequence *seq, gboolean free_seq_too);
void	sequence_free_preprocessing_data(sequence *seq);
gboolean sequence_is_loaded();
int	sequence_processing(sequence *seq, sequence_proc process, int layer, gboolean run_in_thread, gboolean run_in_parallel, void *arg);
int	seqprocess_fwhm(sequence *seq, int seq_layer, int frame_no, fits *fit, rectangle *source_area, void *arg);
int	do_fwhm_sequence_processing(sequence *seq, int layer, gboolean print_psf, gboolean follow_star, gboolean run_in_thread, gboolean for_registration);
void	check_or_allocate_regparam(sequence *seq, int layer);
sequence *create_internal_sequence(int size);
void	internal_sequence_set(sequence *seq, int index, fits *fit);
int	internal_sequence_find_index(sequence *seq, fits *fit);
gpointer crop_sequence(gpointer p);
gboolean sequence_is_rgb(sequence *seq);
imstats* seq_get_imstats(sequence *seq, int index, fits *the_image, int option);
void	enforce_area_in_image(rectangle *area, sequence *seq);
void	update_export_crop_label();

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
void	fill_sequence_list(sequence *seq, int layer);
void	clear_sequence_list();

/****************** statistics_list.h ******************/
void computeStat();

#endif
