#ifndef PROTO_H_
#define PROTO_H_
#include <stdint.h>
#include "core/siril.h"
#ifdef HAVE_LIBTIFF
#define uint64 uint64_hack_
#define int64 int64_hack_
#include <tiffio.h>
#undef uint64
#undef int64
#endif

#ifdef __cplusplus
extern "C" {
#endif

/****************** image_formats_internal.h ******************/
/* BMP */
int readbmp(const char*, fits*);
int savebmp(const char*, fits*);

/* PNM */
int import_pnm_to_fits(const char *filename, fits *fit);
int saveNetPBM(const char *name, fits *fit);

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
int readpic(const char *name, fits *fit);

/****************** image_formats_libraries.h ******************/
#ifdef HAVE_LIBTIFF
int readtif(const char *name, fits *fit, gboolean force_float);
int savetif(const char *name, fits *fit, uint16 bitspersample);
#endif

#ifdef HAVE_LIBJPEG
int readjpg(const char*, fits*);
int savejpg(const char*, fits*, int);
#endif

#ifdef HAVE_LIBPNG
int readpng(const char*, fits*);
int savepng(const char *filename, fits *fit, uint32_t bytes_per_sample,
		gboolean is_colour);
#endif

#ifdef HAVE_LIBRAW
int open_raw_files(const char*, fits*, gboolean);
#endif

#ifdef HAVE_LIBHEIF
int readheif(const char* name, fits *fit, gboolean interactive);
#endif

/****************** utils.h ******************/
int round_to_int(double x);
int roundf_to_int(float x);
WORD round_to_WORD(double x);
BYTE round_to_BYTE(double x);
BYTE roundf_to_BYTE(float f);
WORD roundf_to_WORD(float f);
BYTE conv_to_BYTE(double x);
int truncate_to_int32(uint64_t x);
WORD truncate_to_WORD(int x);
float set_float_in_interval(float val, float low, float high);
double set_double_in_interval(double val, double low, double high);
float ushort_to_float_range(WORD w);
float uchar_to_float_range(BYTE w);
float double_ushort_to_float_range(double d);
WORD float_to_ushort_range(float f);
BYTE float_to_uchar_range(float f);
float ushort_to_float_bitpix(fits *fit, WORD value);
WORD *float_buffer_to_ushort(float *buffer, size_t ndata);
float *uchar_buffer_to_float(BYTE *buffer, size_t ndata);
float *ushort_buffer_to_float(WORD *buffer, size_t ndata);
float *ushort8_buffer_to_float(WORD *buffer, size_t ndata);
uint16_t change_endianness16(uint16_t x);
uint16_t cpu_to_le16(uint16_t x);
uint16_t cpu_to_be16(uint16_t x);
uint16_t le16_to_cpu(uint16_t x);
uint16_t be16_to_cpu(uint16_t x);
uint32_t change_endianness32(uint32_t x);
uint32_t cpu_to_le32(uint32_t x);
uint32_t cpu_to_be32(uint32_t x);
uint32_t le32_to_cpu(uint32_t x);
uint32_t be32_to_cpu(uint32_t x);
uint64_t change_endianness64(uint64_t x);
uint64_t cpu_to_le64(uint64_t x);
uint64_t cpu_to_be64(uint64_t x);
uint64_t le64_to_cpu(uint64_t x);
uint64_t be64_to_cpu(uint64_t x);
gboolean isrgb(fits *fit);
gboolean ends_with(const char *str, const char *ending);
int get_extension_index(const char *filename);
image_type get_type_from_filename(const gchar *filename);
int is_readable_file(const char *filename);
gboolean is_forbiden_in_filename(gchar c);
gboolean file_name_has_invalid_chars(const char *name);
int stat_file(const char *filename2, image_type *type, char **realname);
const char* get_filename_ext(const char *filename);

int changedir(const char *dir, gchar **err);
gchar* get_locale_filename(const gchar *path);
int update_sequences_list(const char *sequence_name_to_select);
void expand_home_in_filename(char *filename, int size);
double get_normalized_value(fits*);
void swap_param(double*, double*);
char* remove_ext_from_filename(const char *basename);
gchar* str_append(char **data, const char *newdata);
char* format_basename(char *root);
float compute_slope(WORD *lo, WORD *hi);
void load_css_style_sheet();
double encodeJD(dateTime dt);
gchar *siril_get_file_info(const gchar *filename, GdkPixbuf *pixbuf);
gchar *siril_truncate_str(gchar *str, gint size);
GtkWidget* popover_new(GtkWidget *widget, const gchar *text);

/****************** quantize.h ***************/
int siril_fits_img_stats_ushort(WORD *array, long nx, long ny, int nullcheck,
		WORD nullvalue, long *ngoodpix, WORD *minvalue, WORD *maxvalue,
		double *mean, double *sigma, double *noise1, double *noise2,
		double *noise3, double *noise5, gboolean multithread, int *status);

int siril_fits_img_stats_float(float *array, long nx, long ny, int nullcheck,
		float nullvalue, long *ngoodpix, float *minvalue, float *maxvalue,
		double *mean, double *sigma, double *noise1, double *noise2,
		double *noise3, double *noise5, gboolean multithread, int *status);

/****************** siril.h ******************/

int threshlo(fits *fit, WORD level);
int threshhi(fits *fit, WORD level);
int nozero(fits *fit, WORD level);
int unsharp(fits*, double sigma, double mult, gboolean verbose);
float entropy(fits *fit, int layer, rectangle *area, imstats *opt_stats);
int loglut(fits *fit);
int ddp(fits *a, int lev, float coef, float sig);
int visu(fits *fit, int low, int high);
int fill(fits *fit, int level, rectangle *arearg);
int off(fits *a, float level);
double background(fits *fit, int reqlayer, rectangle *selection, gboolean multithread);
void show_FITS_header(fits*);
void compute_grey_flat(fits *fit);

/****************** seqfile.h ******************/
sequence* readseqfile(const char *name);
int writeseqfile(sequence *seq);
gboolean existseq(const char *name);
int buildseqfile(sequence *seq, int force_recompute);

/****************** registration_preview.h ******************/
void redraw_previews();
void set_preview_area(int preview_area, int centerX, int centerY);
void init_mouse();
void adjust_reginfo();
void on_spinbut_shift_value_change(GtkSpinButton *spinbutton,
		gpointer user_data);
void test_and_allocate_reference_image(int vport);
void enable_view_reference_checkbox(gboolean status);

/****************** statistics_list.h ******************/
void computeStat();

/****************** siril_log.h ******************/
char* siril_log_message(const char* format, ...);
char* siril_log_color_message(const char* format, const char* color, ...);

#ifdef __cplusplus
}
#endif

#endif
