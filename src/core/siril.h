#ifndef SIRIL_H
#define SIRIL_H
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <gtk/gtk.h>
#include <fitsio.h>	// fitsfile
#include <gsl/gsl_histogram.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <libintl.h>

// PATH_MAX is not available on Hurd at least
#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#define _(String) gettext (String)
#define gettext_noop(String) String
#define N_(String) gettext_noop (String)

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
 #define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#define SQR(x) ((x)*(x))

#define USHRT_MAX_DOUBLE ((double)USHRT_MAX)
#define USHRT_MAX_SINGLE ((float)USHRT_MAX)
#define UCHAR_MAX_DOUBLE ((double)UCHAR_MAX)
#define UCHAR_MAX_SINGLE ((float)UCHAR_MAX)

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

#define STATS_BASIC		(1 << 1)	// median, mean, sigma, noise
#define STATS_AVGDEV	(1 << 2)	// average absolute deviation
#define STATS_MAD		(1 << 3)	// median absolute deviation
#define STATS_BWMV		(1 << 5)	// bidweight midvariance
#define STATS_MAIN		STATS_BASIC | STATS_AVGDEV | STATS_MAD | STATS_BWMV
#define STATS_IKSS		(1 << 6)	// Iterative K-sigma Estimator of Location and Scale. Take time, needed only for stacking
#define STATS_EXTRA		STATS_MAIN | STATS_IKSS

#define	STATS_ZERO_NONE 0
#define	STATS_ZERO_NULLCHECK (!STATS_ZERO_NONE)


/* when requesting an image redraw, it can be asked to remap its data before redrawing it.
 * REMAP_NONE	doesn't remaps the data,
 * REMAP_ONLY	remaps only the current viewport (color channel) and the mixed (RGB) image
 * REMAP_ALL	remaps all view ports, useful when all the colors come from the same file.
 */
#define REMAP_NONE	0
#define REMAP_ONLY	1
#define REMAP_ALL	2

enum {
	COLUMN_FILENAME,		// string
	COLUMN_DATE,		// string
	N_COLUMNS_CONVERT
};


typedef enum {
	TYPEUNDEF=(1 << 1),
	TYPEFITS= (1 << 2),
	TYPETIFF= (1 << 3),
	TYPEBMP = (1 << 4),
	TYPEPNG = (1 << 5),
	TYPEJPG = (1 << 6),
	TYPEPNM = (1 << 7),
	TYPEPIC = (1 << 8),
	TYPERAW = (1 << 9),
	TYPEAVI = (1 << 10),
	TYPESER = (1 << 11),
	TYPEMP4 = (1 << 12),
	TYPEWEBM= (1 << 13),
} image_type;

#define USE_DARK	0x01
#define USE_FLAT	0x02
#define USE_OFFSET	0x04
#define USE_COSME   0x08	/* cosmetic correction */
#define USE_OPTD    0x10	/* dark optimization */

/* cookies for the file chooser */
#define OD_NULL 	0
#define OD_FLAT 	1
#define OD_DARK 	2
#define OD_OFFSET 	3
#define OD_CWD 		4
#define OD_OPEN 	5
#define OD_CONVERT 	6

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
#define MAXGRAYVPORT 	3			// 3 gray vports supported only (R, G, B)
#define MAXCOLORVPORT	1			// 1 color vport supported only (RGB)
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

#define RESULT_IMAGE -1		// used as current image index in a sequence
				// when the result of a processing is displayed
#define UNRELATED_IMAGE -2	// image loaded while a sequence was loaded too

#define MAX_STARS 50000

/* constants for loglut function */
#define LOG 1
#define EXP -1

#define PROGRESS_NONE -2.0		// don't update the progress bar value
#define PROGRESS_PULSATE -1.0		// pulsate the progress bar
#define PROGRESS_RESET 0.0		// reset the progress bar
#define PROGRESS_DONE 1.0		// fill the progress bar
#define PROGRESS_TEXT_RESET ""		// reset the progress bar's text

typedef struct imdata imgdata;
typedef struct registration_data regdata;
typedef struct layer_info_struct layer_info;
typedef struct sequ sequence;
typedef struct single_image single;
typedef struct ffit fits;
typedef struct libraw_config libraw;
typedef struct stack_config stackconf;
typedef struct cominf cominfo;
typedef struct image_stats imstats;
typedef struct rectangle_struct rectangle;
typedef struct point_struct point;
typedef struct gradient_struct gradient;
typedef struct historic_struct historic;
typedef struct fwhm_struct fitted_PSF;

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
	BAYER_SUPER_PIXEL
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
    BAYER_FILTER_NONE = -1		//case where bayer pattern is undefined or untested
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
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
	SEQ_AVI,
#endif
	SEQ_INTERNAL
} sequence_type;

/* image data, exists once for each image */
struct imdata {
	int filenum;		/* real file index in the sequence, i.e. for mars9.fit = 9 */
	gboolean incl;		/* selected in the sequence, included for future processings? */
	//double eval;		/* an evaluation of the quality of the image, not yet used */
	/* for deep-sky images, the quality is the FWHM of a star or an average of FWHM of stars,
	 * for planetary and deep-sky images, it's computed from entropy, contrast or other eval functions. */
	imstats *stats;		/* statistics of the image, used as a cache for full image first
				   channel statistics, particularly used in stacking normalization.
				   NULL if not available, this is stored in seqfile if non NULL */
};

/* Procedure signature for sequence processing.
 * Returns < 0 for an error that should stop the processing on the sequence.
 * Other return values are not used.
 * Processed data should be written in the sequence data directly. */
typedef int (*sequence_proc)(sequence *seq, int seq_layer, int frame_no, fits *fit, rectangle *source_area, void *arg);

/* preprocessing data from GUI */
struct preprocessing_data {
	struct timeval t_start;
	gboolean autolevel;
	double sigma[2];
	gboolean is_cfa;
	float normalisation;
	int retval;
};

/* registration data, exists once for each image and each layer */
struct registration_data {
	int shiftx, shifty;	// we could have a subpixel precision, but is it needed? saved
	float rot_centre_x, rot_centre_y;	// coordinates for the rotation centre, saved
	float angle;		// angle for the rotation, saved
	fitted_PSF *fwhm_data;	// used in PSF/FWHM registration, not saved
	float fwhm;		// copy of fwhm->fwhmx, used as quality indicator, saved data
	//double entropy;		// used in DFT registration, saved data
	// entropy could be replaced by a more general quality indicator at
	// some point. It's the only double value stored in the .seq files, others are single.
	double quality;
};

struct sequ {
	char *seqname;		// name of the sequence, as in name.seq
	int number;		// number of images in the sequence
	int selnum;		// number of selected images in the sequence
	int fixed;		// fixed length of image index in filename (like %3d)
	int nb_layers;		// number of layers embedded in each image file
	unsigned int rx;	// first image width
	unsigned int ry;	// first image height
	layer_info *layers;	// info about layers
	int reference_image;	// reference image for registration
	imgdata *imgparam;	// a structure for each image of the sequence
	regdata **regparam;	// *regparam[nb_layers]
	/* beg and end are used prior to imgparam allocation, hence their usefulness */
	int beg;		// imgparam[0]->filenum
	int end;		// imgparam[number-1]->filenum

	/* registration previsualisation and manual alignment data */
	int previewX[PREVIEW_NB], previewY[PREVIEW_NB];	// center, -1 is uninitialized value
	int previewW[PREVIEW_NB], previewH[PREVIEW_NB];	// 0 is uninitialized value

	sequence_type type;
	struct ser_struct *ser_file;
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
	struct film_struct *film_file;
	char *ext;		// extension of video, NULL if not video
#endif
	fits **internal_fits;	// for INTERNAL sequences: images references. Length: number
	fitsfile **fptr;	// file descriptors for open-mode operations
#ifdef _OPENMP
	omp_lock_t *fd_lock;	// locks for open-mode threaded operations
#endif

	fits *offset;		// the image containing offset data
	fits *dark;		// the image containing dark data
	fits *flat;		// the image containing flat data
	char *ppprefix;		// prefix for filename output of preprocessing
	int current;		// file number currently loaded in wfit (or displayed)
	//struct registration_method reg_method;	// is it the right place for that?
	
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

	/* enabling pre-processing on a single image */
	fits *offset;		// the image containing offset data
	fits *dark;		// the image containing dark data
	fits *flat;		// the image containing flat data
	char *ppprefix;		// prefix for filename output of preprocessing
};

struct ffit {
	unsigned int rx;	// image width	(naxes[0])
	unsigned int ry;	// image height	(naxes[1])
	int bitpix;
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
	WORD lo;	// MIPS-LO key in FITS file, which is "Lower visualization cutoff"
	WORD hi;	// MIPS-HI key in FITS file, which is "Upper visualization cutoff"
	float pixel_size_x, pixel_size_y;	// XPIXSZ and YPIXSZ keys
	unsigned int binning_x, binning_y;		// XBINNING and YBINNING keys
	char date_obs[FLEN_VALUE];		// YYYY-MM-DDThh:mm:ss observation start, UT
	char date[FLEN_VALUE];		// YYYY-MM-DDThh:mm:ss creation of file, UT
	char instrume[FLEN_VALUE];		// INSTRUME key
	char telescop[FLEN_VALUE];		// TELESCOP key
	char observer[FLEN_VALUE];		// OBSERVER key
	char bayer_pattern[FLEN_VALUE];	// BAYERPAT key Bayer Pattern if available
	/* data obtained from FITS or RAW files */
	double focal_length, iso_speed, exposure, aperture, ccd_temp;

	/* data used in the Fourier space */
	double dft_norm[3];			// Normalization value
	char dft_type[FLEN_VALUE];		// spectrum, phase
	char dft_ord[FLEN_VALUE];		// regular, centered
	unsigned int dft_rx, dft_ry;		// padding: original value of picture size
	
	/* data computed or set by Siril */
	unsigned short min[3];	// min for all layers
	unsigned short max[3];	// max for all layers
	unsigned short maxi;	// max of the max[3]
	unsigned short mini;	// min of the min[3]

	fitsfile *fptr;		// file descriptor. Only used for file read and write.
	WORD *data;		// 16-bit image data (depending on image type)
	WORD *pdata[3];		// pointers on data, per layer data access (RGB)
	char *header;		// entire header of the FITS file. NULL for non-FITS file.
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

struct debayer_config {
	gboolean open_debayer;			// debayer images being opened
	gboolean use_bayer_header;		// use the pattern given in the file header
	sensor_pattern bayer_pattern;		// user-defined Bayer pattern
	interpolation_method bayer_inter;	// interpolation method for non-libraw debayer
	gboolean compatibility;				// ensure KSTARS compatibility if TRUE
};

struct stack_config {
	int method;				// 0=sum, 1=median, 2=average, 3=pixel max - Use to save preferences in the init file
	int normalisation_method;
	int rej_method;
	double memory_percent;			// percent of available memory to use for stacking
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
	gboolean leveldrag;		// middle click being dragged if true
	int preprostatus;
	gboolean prepro_cfa;	// Use to save type of sensor for cosmetic correction in preprocessing
	gboolean show_excluded;		// show excluded images in sequences
	double zoom_value;		// 1.0 is normal zoom, use get_zoom_val() to access it

	/* selection rectangle for registration, FWHM, PSF */
	gboolean drawn;			// true if the selection rectangle has been drawn TO REMOVE!
	gboolean drawing;		// true if the rectangle is being set (clicked motion)
	gint startX, startY;		// where the mouse was originally clicked to
	rectangle selection;		// coordinates of the selection rectangle

	/* alignment preview data */
	//guchar *preview_buf[PREVIEW_NB];
	cairo_surface_t *preview_surface[PREVIEW_NB];
	GtkWidget *preview_area[PREVIEW_NB];
	guchar *refimage_regbuffer;	// the graybuf[registration_layer] of the reference image
	cairo_surface_t *refimage_surface;

	char *wd;			// working directory, where images and sequences are
	char *initfile;			// the path of the init file
	
	char *ext;		// FITS extension used in SIRIL
	int len_ext;

	int reg_settings;		// Use to save registration method in the init file
	
	gboolean dontShowConfirm;

	gboolean have_dark_theme;	// we have a dark theme, use bright colours

	stackconf stack;
	
	int filter;

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
	char *swap_dir;

	libraw raw_set;			// the libraw settings
	struct debayer_config debayer;	// debayer settings

	sequence seq;			// currently loaded sequence	TODO: *seq
	single *uniq;			// currently loaded image, if outside sequence

	gsl_histogram *layers_hist[MAXVPORT]; // current image's histograms

	fitted_PSF **stars;		// list of stars detected in the current image
	gboolean star_is_seqdata;	// the only star in stars belongs to seq, don't free it
	int selected_star;		// current selected star in the GtkListStore
	double magOffset;		// offset to reduce the real magnitude, single image
	
	gradient *grad;
	int grad_nb_boxes, grad_size_boxes;
	gboolean grad_boxes_drawn;

	GThread *thread;		// the thread for processing
	GMutex mutex;			// a mutex we use for this thread
	gboolean run_thread;		// the main thread loop condition
	int max_thread;			// maximum of thread used
};

/* this structure is used to characterize the statistics of the image */
struct image_stats {
	long total, ngoodpix;
	double mean, avgDev, median, sigma, bgnoise, min, max, normValue, mad, sqrtbwmv,
			location, scale;
	char layername[6];
};

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
extern fits wfit[5];		// used for temp files, can probably be replaced by local variables
extern char **supported_extensions;
extern char *filter_pattern[];
#endif

#endif /*SIRIL */
