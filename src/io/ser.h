#ifndef _SER_H_
#define _SER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "core/siril.h"
#include "io/seqwriter.h"

/* This file is part of Siril, https://free-astro.org/
 *
 * WARNING: the code in this file and its .c counterpart will not work properly
 * on big endian systems.
 */

#if defined(__unix__) || (defined(OS_OSX))
#include <sys/param.h>		// define or not BSD macro
#endif

#ifdef __linux__
#define fseek64 fseeko  // Linux
#define ftell64 ftello  // Linux
#elif defined (OS_OSX)
#define fseek64 fseeko  // OS X
#define ftell64 ftello  // OS X
#elif defined(BSD)
#define fseek64 fseeko  // DragonFly BSD, FreeBSD, OpenBSD, NetBSD
#define ftell64 ftello  // DragonFly BSD, FreeBSD, OpenBSD, NetBSD
#elif defined (__FreeBSD_kernel__) && defined (__GLIBC__)
#define fseek64 fseeko  // kFreeBSD
#define ftell64 ftello  // kFreeBSD
#elif defined (__gnu_hurd__)
#define fseek64 fseeko  // GNU/Hurd
#define ftell64 ftello  // GNU/Hurd
#elif defined(__CYGWIN__)
#define fseek64 fseeko  // CYGWIN
#define ftell64 ftello  // CYGWIN
#else
#define fseek64 _fseeki64  // Windows
#define ftell64 _ftelli64  // Windows
#endif


#define SER_HEADER_LEN 178

typedef enum {
	SER_MONO = 0,
	SER_BAYER_RGGB = 8,
	SER_BAYER_GRBG = 9,
	SER_BAYER_GBRG = 10,
	SER_BAYER_BGGR = 11,
	SER_BAYER_CYYM = 16,
	SER_BAYER_YCMY = 17,
	SER_BAYER_YMCY = 18,
	SER_BAYER_MYYC = 19,
	SER_RGB = 100,	// SER v3
	SER_BGR = 101	// SER v3
} ser_color;

/* Endianness of the frame data for 16-bit images */
typedef enum {
//	SER_BIG_ENDIAN = 0 /* = FALSE */, SER_LITTLE_ENDIAN = 1 /* = TRUE */
	/* For an unknown reason, several of the first programs to support SER
	 * disrespect the specification regarding the endianness flag. The specification
	 * states that a boolean value is used for the LittleEndian header, and they
	 * use it as a BigEndian header, with 0 for little-endian and 1 for big-endian.
	 * Consequently, to not break compatibility with these first implementations,
	 * later programs, like Siril and GoQat, have also decided to implement this
	 * header in opposite meaning to the specification. */
	SER_LITTLE_ENDIAN = 0, SER_BIG_ENDIAN = 1
} ser_endian;

typedef enum {
	SER_PIXEL_DEPTH_8 = 1, SER_PIXEL_DEPTH_16 = 2
} ser_pixdepth;


/* this struct does not reflect the exact header:
 * - strings are nul-terminated, which makes them one character larger than in
 *   the header
 * - the integer types of the header are little endian, which may not be the
 *   case when compiling the struct
 */
struct ser_struct {
	char *file_id;			// 14 bytes (0)
	int lu_id;			// 4	(14)
	ser_color color_id;		// 4	(18)
	ser_endian little_endian;	// 4	(22)
	int image_width;		// 4	(26)
	int image_height;		// 4	(30)
	int bit_pixel_depth;		// 4	(34)
	unsigned int frame_count;	// 4	(38)
	char observer[40];		// 40	(42)
	char instrument[40];		// 40	(82)
	char telescope[40];		// 40	(122)
	guint64 date;		// 8	(162)
	guint64 date_utc;	// 8 (170)

	/* timestamps (not in the header, timestamps are in trailer) */
	guint64 *ts;			// total timestamps
	int ts_alloc;			// allocated number of elements in ts
	guint64 ts_min, ts_max;// min and max timestamp
	double fps;				// frame rate

	gint64 filesize;			// size of the file

	// internal representations of header data
	ser_pixdepth byte_pixel_depth;	// more useful representation of the bit_pixel_depth
	unsigned int number_of_planes;	// derived from the color_id
	FILE *file;
	char *filename;
#ifdef _OPENMP
	omp_lock_t fd_lock, ts_lock;
#endif

	struct seqwriter_data *writer;
};

gboolean ser_is_cfa(struct ser_struct *ser_file);
void ser_convertTimeStamp(struct ser_struct *ser_file, GSList *timestamp);
void ser_init_struct(struct ser_struct *ser_file);
void ser_display_info(struct ser_struct *ser_file);
int ser_open_file(const char *filename, struct ser_struct *ser_file);
int ser_close_and_delete_file(struct ser_struct *ser_file);
int ser_write_and_close(struct ser_struct *ser_file);
int ser_create_file(const char *filename, struct ser_struct *ser_file, gboolean overwrite, struct ser_struct *copy_from);
int ser_close_file(struct ser_struct *ser_file);
int ser_metadata_as_fits(struct ser_struct *ser_file, fits *fit);
int ser_read_frame(struct ser_struct *ser_file, int frame_no, fits *fit, gboolean force_float);
int ser_read_opened_partial_fits(struct ser_struct *ser_file, int layer,
		int frame_no, fits *fit, const rectangle *area);
int ser_read_opened_partial(struct ser_struct *ser_file, int layer,
		int frame_no, WORD *buffer, const rectangle *area);
int ser_write_frame_from_fit(struct ser_struct *ser_file, fits *fit, int frame);
gint64 ser_compute_file_size(struct ser_struct *ser_file, int nb_frames);
int import_metadata_from_serfile(struct ser_struct *ser_file, fits *to);
GdkPixbuf* get_thumbnail_from_ser(char *filename, gchar **descr);
GDateTime *ser_read_frame_date(struct ser_struct *ser_file, int frame_no);

#endif

