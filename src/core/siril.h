#ifndef SIRIL_H
#define SIRIL_H
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <glib.h>
#include <glib/gstdio.h>
#include <glib/gprintf.h>
#include <gtk/gtk.h>
#include <fitsio.h>	// fitsfile
#include <gsl/gsl_histogram.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <libintl.h>

#define _(String) gettext (String)
#define gettext_noop(String) String
#define N_(String) gettext_noop (String)


#ifdef SIRIL_OUTPUT_DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

/* https://stackoverflow.com/questions/1644868/define-macro-for-debug-printing-in-c */
#define siril_debug_print(fmt, ...) \
   do { if (DEBUG_TEST) fprintf(stdout, fmt, ##__VA_ARGS__); } while (0)

#define PRINT_ALLOC_ERR fprintf(stderr, "Out of memory in %s (%s:%d) - aborting\n", __func__, __FILE__, __LINE__)

#undef max
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#undef min
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define SWAP(a,b)  { double temp = (a); (a) = (b); (b) = temp; }

#define SQR(x) ((x)*(x))

#define USHRT_MAX_DOUBLE ((double)USHRT_MAX)
#define USHRT_MAX_SINGLE ((float)USHRT_MAX)
#define UCHAR_MAX_DOUBLE ((double)UCHAR_MAX)
#define UCHAR_MAX_SINGLE ((float)UCHAR_MAX)

#define BYTES_IN_A_MB 1048576	// 1024

#define SEQUENCE_DEFAULT_INCLUDE TRUE	// select images by default

typedef unsigned char BYTE;		// default type for image display data
typedef unsigned short WORD;		// default type for internal image data

#define MAX_COMMAND_WORDS 16		// max number of words to split in command line input

#define MAX_SEQPSF 7			// max number of stars for which seqpsf can be run

#define CMD_HISTORY_SIZE 50		// size of the command line history

#define ZOOM_MAX	16.0
#define ZOOM_MIN	0.0625
#define ZOOM_NONE	1.0
#define ZOOM_FIT	-1.0	// or any value < 0
#define ZOOM_DEFAULT	ZOOM_FIT

#define LOW_BOUND  0.00002
#define HIGH_BOUND 0.99998

/* Some statistic constants */
#define SIGMA_PER_FWHM 2.35482
#define AVGDEV_NORM 1.2533
#define MAD_NORM 1.4826
#define BWMV_NORM 0.9901
#define PBMV_NORM 0.9709
#define SN_NORM 1.1926
#define QN_NORM 2.2191

typedef struct _SirilDialogEntry SirilDialogEntry;

/* used for open and savedialog */
#if (defined _WIN32) || (defined(__APPLE__) && defined(__MACH__))
#define SirilWidget GtkFileChooserNative
#define SIRIL_EOL "\r\n"
#else
#define SirilWidget GtkWidget
#define SIRIL_EOL "\n"
#endif

/* when requesting an image redraw, it can be asked to remap its data before redrawing it.
 * REMAP_NONE	doesn't remaps the data,
 * REMAP_ONLY	remaps only the current viewport (color channel) and the mixed (RGB) image
 * REMAP_ALL	remaps all view ports, useful when all the colors come from the same file.
 */
#define REMAP_NONE	0
#define REMAP_ONLY	1
#define REMAP_ALL	2

typedef enum {
	TYPEUNDEF = (1 << 1),
	TYPEFITS = (1 << 2),
	TYPETIFF = (1 << 3),
	TYPEBMP = (1 << 4),
	TYPEPNG = (1 << 5),
	TYPEJPG = (1 << 6),
	TYPEPNM = (1 << 7),
	TYPEPIC = (1 << 8),
	TYPERAW = (1 << 9),
	TYPEAVI = (1 << 10),
	TYPESER = (1 << 11),
	TYPEMP4 = (1 << 12),
	TYPEWEBM = (1 << 13),
} image_type;

/* indices of the image data layers */
#define BW_LAYER 	0
#define RLAYER		0
#define GLAYER		1
#define BLAYER		2
#define RGB_LAYER 	3

/* indices of the viewports (graphical elements) */
#define BW_VPORT 	0
#define RED_VPORT 	0
#define GREEN_VPORT 	1
#define BLUE_VPORT 	2
#define RGB_VPORT 	3
#define MAXGRAYVPORT 	3	// 3 gray vports supported only (R, G, B)
#define MAXCOLORVPORT	1	// 1 color vport supported only (RGB)
#define MAXVPORT 	MAXGRAYVPORT + MAXCOLORVPORT

/* defines for copyfits actions */
#define CP_INIT		0x01	// initialize data array with 0s
#define CP_ALLOC	0x02	// reallocs data array
#define CP_COPYA	0x04	// copy data array content
#define CP_FORMAT	0x08	// copy size and bitpix info
#define CP_EXTRACT	0x10	// extract a 16bit plane from a 48 bit fit
#define CP_EXPAND	0x20	// expands a 16bit fits to a 48bit one.

/* processing */
#define CONVDEBAYER (1 << 1)
/* SER flags */
#define CONVDSTFITS (1 << 2)	// assumed as default
#define CONVDSTSER (1 << 3)
#define CONVMULTIPLE (1 << 4)
/* channel conversion type */
#define CONV1X3	(1 << 6)	// assumed as default
#define CONV3X1	(1 << 7)
#define CONV1X1	(1 << 8)

/* operations on image data */
#define OPER_ADD 'a'
#define OPER_SUB 's'
#define OPER_MUL 'm'
#define OPER_DIV 'd'

#define PREVIEW_NB 2

/* special values for com.seq.current, the currently loaded image of the
 * sequence. Negative values can be used to indicate that a specific image is
 * loaded while there is a sequence or not in com.seq */
#define RESULT_IMAGE -1		// used as current image index in a sequence
				// when the result of a processing is displayed
#define UNRELATED_IMAGE -2	// image loaded while a sequence was loaded too
#define SCALED_IMAGE -3		// the single image has a different size than
				// the loaded sequence

#define MAX_STARS 200000		// maximum length of com.stars

typedef struct imdata imgdata;
typedef struct registration_data regdata;
typedef struct layer_info_struct layer_info;
typedef struct sequ sequence;
typedef struct single_image single;
typedef struct wcs_struct wcs_info;
typedef struct dft_struct dft_info;
typedef struct ffit fits;
typedef struct libraw_config libraw;
typedef struct phot_config phot;
typedef struct stack_config stackconf;
typedef struct cominf cominfo;
typedef struct image_stats imstats;
typedef struct rectangle_struct rectangle;
typedef struct point_struct point;
typedef struct historic_struct historic;
typedef struct dateTime_struct dateTime;
typedef struct fwhm_struct fitted_PSF;
typedef struct star_finder_struct star_finder_params;

/* global structures */

/* ORDER OF POLYNOMES */
typedef enum {
	POLY_1,
	POLY_2,
	POLY_3,
	POLY_4,
} poly_order;

typedef enum {
	NORMAL_DISPLAY,	
	LOG_DISPLAY,
	SQRT_DISPLAY,
	SQUARED_DISPLAY,
	ASINH_DISPLAY,
	STF_DISPLAY,
	HISTEQ_DISPLAY
} display_mode;			// used in the layer_info_struct below
#define DISPLAY_MODE_MAX HISTEQ_DISPLAY

typedef enum {
	NORMAL_COLOR,	
	RAINBOW_COLOR
} color_map;

typedef enum {
	MIPSLOHI,
	MINMAX,
	USER
} sliders_mode;

typedef enum {
	FILE_CONVERSION,
	IMAGE_SEQ,
	PRE_PROC,
	REGISTRATION,
	PLOT,
	STACKING,
	OUTPUT_LOGS
} main_tabs;

typedef enum {
	BAYER_BILINEAR,
	BAYER_NEARESNEIGHBOR,
	BAYER_VNG,
	BAYER_AHD,
	BAYER_SUPER_PIXEL,
	XTRANS
} interpolation_method;

typedef enum {
	OPENCV_NEAREST = 0,
	OPENCV_LINEAR = 1,
	OPENCV_AREA = 2,
	OPENCV_CUBIC = 3,
	OPENCV_LANCZOS4 = 4,
	OPENCV_INTER_MAX = 7
} opencv_interpolation;

typedef enum {
    BAYER_FILTER_RGGB,
    BAYER_FILTER_BGGR,
    BAYER_FILTER_GBRG,
    BAYER_FILTER_GRBG,
	XTRANS_FILTER,
    BAYER_FILTER_NONE = -1		//case where pattern is undefined or untested
} sensor_pattern ;
#define BAYER_FILTER_MIN BAYER_FILTER_RGGB
#define BAYER_FILTER_MAX BAYER_FILTER_GRBG

struct layer_info_struct {
	char *name;			// name of the layer (a filter name)
	double wavelength;		// the wavelength of the filter, in nanometres
	WORD lo, hi;			// the values of the cutoff sliders
	//WORD min, max;		// the min and max values of the sliders
	gboolean cut_over, cut_under;	// display values over hi or under lo as negative
	display_mode rendering_mode;	// defaults to NORMAL_DISPLAY
};

typedef enum { SEQ_REGULAR, SEQ_SER,
#ifdef HAVE_FFMS2
	SEQ_AVI,
#endif
	SEQ_INTERNAL
} sequence_type;

/* image data, exists once for each image */
struct imdata {
	int filenum;		/* real file index in the sequence, i.e. for mars9.fit = 9 */
	gboolean incl;		/* selected in the sequence, included for future processings? */
	char *date_obs;		/* date of the observation, processed and copied from the header */
};

/* registration data, exists once for each image and each layer */
struct registration_data {
	float shiftx, shifty;	// we could have a subpixel precision, but is it needed? saved
	fitted_PSF *fwhm_data;	// used in PSF/FWHM registration, not saved
	float fwhm;		// copy of fwhm->fwhmx, used as quality indicator, saved data
	float roundness;	// fwhm->fwhmy / fwhm->fwhmx, 0 when uninit, ]0, 1] when set 
	double quality;
};

/* see explanation about sequence and single image management in io/sequence.c */

struct sequ {
	char *seqname;		// name of the sequence, as in name.seq
	int number;		// number of images in the sequence
	int selnum;		// number of selected images in the sequence
	int fixed;		// fixed length of image index in filename (like %3d)
	int nb_layers;		// number of layers embedded in each image file, -1 if unknown
	unsigned int rx;	// first image width
	unsigned int ry;	// first image height
	int bitpix;		// image pixel format, from fits
	double data_max; // data_max used for cdata conversion
	layer_info *layers;	// info about layers, may be null if nb_layers is unknown
	int reference_image;	// reference image for registration
	imgdata *imgparam;	// a structure for each image of the sequence
	regdata **regparam;	// *regparam[nb_layers], may be null if nb_layers is unknown
	imstats ***stats;	// statistics of the images for each layer, may be null too
	/* in the case of a CFA sequence, depending on the opening mode, we cannot store
	 * and use everything that was in the seqfile, so we back them up here */
	regdata **regparam_bkp;	// *regparam[3], null if nothing to back up
	imstats ***stats_bkp;	// statistics of the images for 3 layers, may be null too

	/* beg and end are used prior to imgparam allocation, hence their usefulness */
	int beg;		// imgparam[0]->filenum
	int end;		// imgparam[number-1]->filenum
	double exposure;	// exposure of frames (we assume they are all identical)

	sequence_type type;
	struct ser_struct *ser_file;
	gboolean cfa_opened_monochrome;	// in case the CFA SER was opened in monochrome mode
#ifdef HAVE_FFMS2
	struct film_struct *film_file;
	char *ext;		// extension of video, NULL if not video
#endif
	fits **internal_fits;	// for INTERNAL sequences: images references. Length: number
	fitsfile **fptr;	// file descriptors for open-mode operations
#ifdef _OPENMP
	omp_lock_t *fd_lock;	// locks for open-mode threaded operations
#endif

	int current;		// file number currently loaded in gfit (or displayed)

	/* registration previsualisation and manual alignment data */
	int previewX[PREVIEW_NB], previewY[PREVIEW_NB];	// center, -1 is uninitialized value
	int previewW[PREVIEW_NB], previewH[PREVIEW_NB];	// 0 is uninitialized value

	double upscale_at_stacking;// up-scale factor during stacking (see #215)
	
	gboolean needs_saving;	// a dirty flag for the sequence, avoid saving it too often

	fitted_PSF **photometry[MAX_SEQPSF];// psf for multiple stars for all images
	int reference_star;	// reference star for apparent magnitude (index of photometry)
	double reference_mag;	// reference magnitude for the reference star
	double photometry_colors[MAX_SEQPSF][3]; // colors for each photometry curve
};

/* this struct is used to manage data associated with a single image loaded, outside a sequence */
struct single_image {
	char *filename;		// the name of the file
	char *comment;		// comment on how the file got there (user load, result...)
	int nb_layers;		// number of layers embedded in each image file
	layer_info *layers;	// info about layers
	fits *fit;		// the fits is still gfit, but a reference doesn't hurt
};

struct wcs_struct {
	unsigned int equinox;
	double crpix1, crpix2;
	double crval1, crval2;
	double cdelt1, cdelt2;
	double crota1, crota2;
	char objctra[FLEN_VALUE];
	char objctdec[FLEN_VALUE];
};

struct dft_struct {
	double norm[3];			// Normalization value
	char type[FLEN_VALUE];		// spectrum, phase
	char ord[FLEN_VALUE];		// regular, centered
};

struct ffit {
	unsigned int rx;	// image width	(naxes[0])
	unsigned int ry;	// image height	(naxes[1])
	int bitpix;		// current bitpix of loaded data
	int orig_bitpix;	// original bitpix of the file
	/* bitpix can take the following values:
	 * BYTE_IMG	(8-bit byte pixels, 0 - 255)
	 * SHORT_IMG	(16 bit signed integer pixels)	
	 * USHORT_IMG	(16 bit unsigned integer pixels)	(used by Siril, quite off-standard)
	 * LONG_IMG	(32-bit integer pixels)
	 * FLOAT_IMG	(32-bit floating point pixels)
	 * DOUBLE_IMG	(64-bit floating point pixels)
	 * http://heasarc.nasa.gov/docs/software/fitsio/quick/node9.html
	 */
	int naxis;		// number of dimensions of the image
	long naxes[3];		// size of each dimension
	/* naxes[0] is rx, naxes[1] is ry
	 * Then, for gray images, naxes[2] is unused in FITS but set to 1, and naxis is 2.
	 * For RGB images, naxes[2] is 3 and naxis is 3.
	 * */

	/* data obtained from the FITS file */
	char *header;	// entire header of the FITS file. NULL for non-FITS file.
	WORD lo;	// MIPS-LO key in FITS file, which is "Lower visualization cutoff"
	WORD hi;	// MIPS-HI key in FITS file, which is "Upper visualization cutoff"
	double data_max; // used to check if 32b float is between 0 and 1
	WORD maximum_pixel_value; // value obtained from libraw, Maximum pixel value. Calculated from the data for most cameras, hardcoded for others.
	float pixel_size_x, pixel_size_y;	// XPIXSZ and YPIXSZ keys
	unsigned int binning_x, binning_y;		// XBINNING and YBINNING keys
	gboolean unbinned;
	char date_obs[FLEN_VALUE];		// YYYY-MM-DDThh:mm:ss observation start, UT
	char date[FLEN_VALUE];		// YYYY-MM-DDThh:mm:ss creation of file, UT
	char instrume[FLEN_VALUE];		// INSTRUME key
	char telescop[FLEN_VALUE];		// TELESCOP key
	char observer[FLEN_VALUE];		// OBSERVER key
	char bayer_pattern[FLEN_VALUE];	// BAYERPAT key Bayer Pattern if available
	int bayer_xoffset, bayer_yoffset;
	/* data obtained from FITS or RAW files */
	double focal_length, iso_speed, exposure, aperture, ccd_temp;
	double cvf; // Conversion factor (e-/adu)

	/* Plate Solving data */
	wcs_info wcs;

	/* data used in the Fourier space */
	dft_info dft;
	
	/* data computed or set by Siril */
	imstats **stats;	// stats of fit for each layer, null if naxes[2] is unknown
	double mini, maxi;	// min and max of the stats->max[3]

	fitsfile *fptr;		// file descriptor. Only used for file read and write.
	WORD *data;		// 16-bit image data (depending on image type)
	WORD *pdata[3];		// pointers on data, per layer data access (RGB)

	gboolean top_down;	// image data is stored top-down, normally false for FITS, true for SER

	GSList *history;	// Former HISTORY comments of FITS file
};

/* This structure is used for all the elements in the box libraw_settings.
 * Don't forget to update conversion.c:initialize_libraw_settings() data when
 * modifying the glade settings */
struct libraw_config {
	double mul[3], bright;					// Color  & brightness adjustement mul[0] = red, mul[1] = green = 1, mul[2] = blue
	int auto_mul, use_camera_wb, use_auto_wb;		// White Balance parameters
	int user_qual;						// Index of the Matrix interpolation set in dcraw, 0: bilinear, 1: VNG, 2: PPG, 3: AHD
	int user_black;						// black point correction
	double gamm[2];						// Gamma correction
};

/* This structure is used for storing all parameters used in photometry module */
struct phot_config {
	double gain;	// A/D converter gain in electrons per ADU
	double inner;	// Inner radius of the annulus used to measure local background.
	double outer;	// Outer radius of the annulus used to measure local background.
	int minval, maxval;
};

struct debayer_config {
	gboolean open_debayer;			// debayer images being opened
	gboolean use_bayer_header;		// use the pattern given in the file header
	sensor_pattern bayer_pattern;		// user-defined Bayer pattern
	interpolation_method bayer_inter;	// interpolation method for non-libraw debayer
	gboolean compatibility;				// ensure KSTARS compatibility if TRUE
	gboolean stretch;                  // stretch DSLR CFA data to 16-bit if wanted
};

struct stack_config {
	int method;				// 0=sum, 1=median, 2=average, 3=pixel max, 4=pixel min - Use to save preferences in the init file
	int normalisation_method;
	int rej_method;
	enum { RATIO, AMOUNT, UNLIMITED } mem_mode; // mode of memory management
	double memory_ratio;			// ratio of available memory to use for stacking (and others)
	double memory_amount;			// amount of memory in GB to use for stacking (and others)
};

struct rectangle_struct {
	int x, y, w, h;
};

struct point_struct {
	double x, y;
};

struct gradient_struct {
	point centre;
	double boxvalue[3];
};

struct historic_struct {
	char *filename;
	char history[FLEN_VALUE];
	int rx, ry;
};

struct dateTime_struct {
	int year;
	int month;
	int day;
	int hour;
	int min;
	int sec;
	int ms;
};

struct star_finder_struct {
	int radius;
	double sigma;
	double roundness;
};

/* The global data structure of siril, the only with gfit and the gtk builder,
 * declared in main.c */
struct cominf {
	/* current version of GTK, through GdkPixmap, doesn't handle gray images, so
	 * graybufs are the same size than the rgbbuf with 3 times the same value */
	guchar *graybuf[MAXGRAYVPORT];	// one B/W display buffer per viewport (R,G,B)
	guchar *rgbbuf;			// one rgb display buffer
	/* cairo image surface related data */
	int surface_stride[MAXVPORT];	// allocated stride
	int surface_height[MAXVPORT];	// allocated height
	cairo_surface_t *surface[MAXVPORT];
	gboolean buf_is_dirty[MAXVPORT];// dirtyness of each buffer (= need to redraw)
	
	/* Color map */
	color_map color;

	GtkWidget *vport[MAXVPORT];	// one drawingarea per layer, one rgb drawingarea
	int cvport;			// current viewport, index in the list vport above
	GtkAdjustment *hadj[MAXVPORT];	// adjustments of vport scrollbars
	GtkAdjustment *vadj[MAXVPORT];	// adjustments of vport scrollbars
	sliders_mode sliders;		// 0: min/max, 1: MIPS-LO/HI, 2: user
	gboolean prepro_cfa;	// Use to save type of sensor for cosmetic correction in preprocessing
	gboolean prepro_equalize_cfa;  // Use to save if flat will be equalized in preprocessing
	gboolean show_excluded;		// show excluded images in sequences
	double zoom_value;		// 1.0 is normal zoom, use get_zoom_val() to access it

	/* positions of all windows */
	gboolean remember_windows;
	rectangle main_w_pos;
	rectangle rgb_w_pos;

	/* selection rectangle for registration, FWHM, PSF */
	gboolean drawing;		// true if the rectangle is being set (clicked motion)
	gint startX, startY;		// where the mouse was originally clicked to
	gboolean freezeX, freezeY;
	rectangle selection;		// coordinates of the selection rectangle

	/* alignment preview data */
	//guchar *preview_buf[PREVIEW_NB];
	cairo_surface_t *preview_surface[PREVIEW_NB];
	GtkWidget *preview_area[PREVIEW_NB];
	guchar *refimage_regbuffer;	// the graybuf[registration_layer] of the reference image
	cairo_surface_t *refimage_surface;

	gchar *wd;			// working directory, where images and sequences are
	gchar *initfile;	// the path of the init file
	gchar *app_path;	// the path of the application
	
	char *ext;		// FITS extension used in SIRIL

	int reg_settings;		// Use to save registration method in the init file
	
	gboolean dontShowConfirm;

	gboolean have_dark_theme;	// global theme is dark
	gint combo_theme;           // value of the combobox theme
	gboolean want_dark;			// User want dark theme for siril

	stackconf stack;

	gboolean cache_upscaled;	// keep up-scaled files for 'drizzle' (only used by developers)
	
	int filter;			// file extension filter for open/save dialogs

	/* history of the command line. This is a circular buffer (cmd_history)
	 * of size cmd_hist_size, position to be written is cmd_hist_current and
	 * position being browser for display of the history is cmd_hist_display.
	 */
	char **cmd_history;		// the history of the command line
	int cmd_hist_size;		// allocated size
	int cmd_hist_current;		// current command index
	int cmd_hist_display;		// displayed command index

	/* history of operations */
	historic *history;			// the history of all operations
	int hist_size;			// allocated size
	int hist_current;		// current index
	int hist_display;		// displayed index
	gchar *swap_dir;		// swap directory
	GSList *script_path;	// script path directories

	libraw raw_set;			// the libraw settings
	struct debayer_config debayer;	// debayer settings
	phot phot_set;          // photometry settings

	sequence seq;			// currently loaded sequence	TODO: *seq
	single *uniq;			// currently loaded image, if outside sequence

	gsl_histogram *layers_hist[MAXVPORT]; // current image's histograms

	star_finder_params starfinder_conf;	// star finder settings, from GUI or init file
	fitted_PSF **stars;		// list of stars detected in the current image
	gboolean star_is_seqdata;	// the only star in stars belongs to seq, don't free it
	int selected_star;		// current selected star in the GtkListStore
	double magOffset;		// offset to reduce the real magnitude, single image
	
	GSList *grad_samples;

	int max_thread;			// maximum of thread used for parallel execution

	GThread *thread;		// the thread for processing
	GMutex mutex;			// a mutex we use for this thread
	gboolean run_thread;		// the main thread loop condition

	gboolean headless;		// pure console, no GUI
	gboolean script;		// scripts execution
	gboolean stop_script;		// abort script execution
	GThread *script_thread;		// reads a script and executes its commands
};

/* this structure is used to characterize the statistics of the image */
struct image_stats {
	long total,	// number of pixels
	     ngoodpix;	// number of non-zero pixels
	double mean, median, sigma, avgDev, mad, sqrtbwmv,
	       location, scale, min, max, normValue, bgnoise;

	int _nb_refs;	// reference counting for data management
};

typedef struct Homo {
	double h00, h01, h02;
	double h10, h11, h12;
	double h20, h21, h22;
	int pair_matched;
	int Inliers;
} Homography;

#if 0
/* TODO: this structure aims to allow the composition of several 1-channel images and make
 * more easy the management of RGB compositing */
typedef struct image_layer_struct image_layer;
struct image_layer_struct {
	char		*layer_name;		/* the name of the layer (the color or band name) */
	fits		*fit;			/* fits data of the layer */
	int		naxis;			/* this image is naxis in fits */
	guchar		*graybuf;		/* mapped image for display purposes */
	int		stride;			/* Cairo data width */
	cairo_surface_t	*surface;		/* Cairo image surface */
	GtkWidget	*vport;			/* the viewport */
	double		wavelength;		/* the wavelength associated with the channel */
	guchar		rmap, gmap, bmap;	/* mapping to rgb colors, initialized function of the wavelength */
	unsigned int	hi, lo;			/* same as fits_data->{hi,lo} but for display/compositing purposes */
	char		*filename;		/* the filename of the fits_data file */
	int		naxis;			/* axis number of the fits file filename, may not be 0 if it's an RGB fits for example */
};
#endif

#ifndef MAIN
extern GtkBuilder *builder;	// get widget references anywhere
extern cominfo com;		// the main data struct
extern fits gfit;		// currently loaded image
extern char **supported_extensions;
extern char *filter_pattern[];
#endif

#endif /*SIRIL */
