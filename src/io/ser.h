#ifndef _SER_H_
#define _SER_H_

#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* This file is part of Siril, https://free-astro.org/
 *
 * WARNING: the code in this file and its .c counterpart will not work properly
 * on big endian systems.
 */

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
	char observer[40];			// 40	(42)
	char instrument[40];		// 40	(82)
	char telescope[40];		// 40	(122)
	uint64_t date;	// 8	(162)
	uint64_t date_utc;	// 8 (170)

	/* timestamps (not in the header, timestamps are in trailer) */
	uint64_t *ts;			// total timestamps
	uint64_t ts_min, ts_max;// min and max timestamp
	double fps;				// frame rate

	off_t filesize;			// size of the file

	// internal representations of header data
	ser_pixdepth byte_pixel_depth;	// more useful representation of the bit_pixel_depth
	unsigned int number_of_planes;	// derived from the color_id
	int fd;
	char *filename;
#ifdef _OPENMP
	omp_lock_t fd_lock;
#endif
};

void ser_convertTimeStamp(struct ser_struct *ser_file, GSList *timestamp);
void ser_init_struct(struct ser_struct *ser_file);
void ser_display_info(struct ser_struct *ser_file);
int ser_open_file(char *filename, struct ser_struct *ser_file);
int ser_write_and_close(struct ser_struct *ser_file);
int ser_create_file(const char *filename, struct ser_struct *ser_file, gboolean overwrite, struct ser_struct *copy_from);
int ser_close_file(struct ser_struct *ser_file);
int ser_read_frame(struct ser_struct *ser_file, int frame_no, fits *fit);
int ser_read_opened_partial(struct ser_struct *ser_file, int layer,
		int frame_no, WORD *buffer, const rectangle *area);
int ser_write_frame_from_fit(struct ser_struct *ser_file, fits *fit, int frame);

#endif

