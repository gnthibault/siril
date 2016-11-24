/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2016 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * WARNING: the code in this file and its header will not work properly
 * on big endian systems.
 */

#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "core/siril.h"
#include "core/proto.h"
#include "gui/callbacks.h"
#include "algos/demosaicing.h"
#include "io/ser.h"

/* 62135596800 sec from year 0001 to 01 janv. 1970 00:00:00 GMT */
static const uint64_t epochTicks = 621355968000000000UL;
static const uint64_t ticksPerSecond = 10000000;

static int ser_write_header(struct ser_struct *ser_file);

/* Given an SER timestamp, return a char string representation
 * MUST be freed
 */
static char *ser_timestamp(uint64_t timestamp) {
	char *str = malloc(64);
	uint64_t t1970_ms = (timestamp - epochTicks) / 10000;
	time_t secs = t1970_ms / 1000;
	int ms = t1970_ms % 1000;
	struct tm *t;

	t = gmtime(&secs);

	sprintf(str, "%04d-%02d-%02d %02d:%02d:%02d.%03d", t->tm_year + 1900,
			t->tm_mon + 1, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec, ms);

	return str;
}

/* Output SER timestamp */
static int display_date(uint64_t timestamp, char *txt) {
	if (timestamp == 0)
		return -1;

	char *str = ser_timestamp(timestamp);
	printf("%s%s\n", txt, str);
	free(str);
	return 0;
}

/* Comes from http://linux.die.net/man/3/timegm
 * This is not thread-safe
 */
static time_t __timegm(struct tm *tm) {
	time_t ret;
	char *tz;

	tz = getenv("TZ");
	setenv("TZ", "", 1);
	tzset();
	ret = mktime(tm);
	if (tz)
		setenv("TZ", tz, 1);
	else
		unsetenv("TZ");
	tzset();
	return ret;
}

/* Convert FITS keyword DATE in a UNIX time format
 * DATE match this pattern: 1900-01-01T00:00:00
 */
static int FITS_date_key_to_Unix_time(char *date, uint64_t *utc,
		uint64_t *local) {
	struct tm timeinfo = { };
	time_t ut, t;
	int year = 0, month = 0, day = 0, hour = 0, min = 0, sec = 0;

	if (date[0] == '\0')
		return -1;

	sscanf(date, "%04d-%02d-%02dT%02d:%02d:%02d", &year, &month, &day, &hour,
			&min, &sec);

	timeinfo.tm_year = year - 1900;
	timeinfo.tm_mon = month - 1;
	timeinfo.tm_mday = day;
	timeinfo.tm_hour = hour;
	timeinfo.tm_min = min;
	timeinfo.tm_sec = sec;

	// Hopefully these are not needed
	timeinfo.tm_wday = 0;
	timeinfo.tm_yday = 0;
	timeinfo.tm_isdst = -1;

	/* get UTC time from timeinfo* */
	ut = __timegm(&timeinfo);
	ut *= ticksPerSecond;
	ut += epochTicks;
	*utc = (uint64_t) ut;

	/* get local time from timeinfo* */
	t = mktime(&timeinfo);
	t *= ticksPerSecond;
	t += epochTicks;
	*local = (uint64_t) t;

	return 0;
}

static char *convert_color_id_to_char(ser_color color_id) {
	switch (color_id) {
	case SER_MONO:
		return "MONO";
	case SER_BAYER_RGGB:
		return "RGGB";
	case SER_BAYER_BGGR:
		return "BGGR";
	case SER_BAYER_GBRG:
		return "GBRG";
	case SER_BAYER_GRBG:
		return "GRBG";
	case SER_BAYER_CYYM:
		return "CYYM";
	case SER_BAYER_YCMY:
		return "YCMY";
	case SER_BAYER_YMCY:
		return "YMCY";
	case SER_BAYER_MYYC:
		return "MYYC";
	case SER_RGB:
		return "RGB";
	case SER_BGR:
		return "BGR";
	default:
		return "";
	}
}

static int ser_read_timestamp(struct ser_struct *ser_file) {
	int frame_size, i;
	gboolean timestamps_in_order = TRUE;
	uint64_t previous_ts = 0L;
	off_t filesize;

	filesize = ser_file->filesize;
	frame_size = ser_file->image_width * ser_file->image_height * ser_file->number_of_planes;
	off_t offset = SER_HEADER_LEN + (off_t) frame_size * (off_t) ser_file->byte_pixel_depth
		* (off_t) ser_file->frame_count;

	/* Check if file is large enough to have timestamps */
	if (filesize >= offset + (8 * ser_file->frame_count)) {
		ser_file->ts = calloc(8, ser_file->frame_count);

		// Seek to start of timestamps
		for (i = 0; i < ser_file->frame_count; i++) {
			if ((off_t) -1 == lseek(ser_file->fd, offset + (i * 8), SEEK_SET)) {
				return -1;
			}
			char timestamp[8];

			memset(timestamp, 0, sizeof(timestamp));
			if (8 == read(ser_file->fd, timestamp, 8)) {
				memcpy(&ser_file->ts[i], timestamp, 8);
			} else {
				// No valid frames per second value can be calculated
				ser_file->fps = -1.0;
				return 0;
			}
		}

		/* Check order of Timestamps */
		uint64_t *ts_ptr = ser_file->ts;
		uint64_t min_ts = *ts_ptr;
		uint64_t max_ts = *ts_ptr;

		for (i = 0; i < ser_file->frame_count; i++) {
			if (*ts_ptr < previous_ts) {
				// Timestamps are not in order
				timestamps_in_order = FALSE;
			}
			previous_ts = *ts_ptr;
			// Keep track of maximum timestamp value
			if (*ts_ptr > max_ts) {
				max_ts = *ts_ptr;
			}
			// Keep track of minimum timestamp value
			if (*ts_ptr < min_ts) {
				min_ts = *ts_ptr;
			}
			ts_ptr++;
		}

		if (timestamps_in_order) {
			if (min_ts == max_ts) {
				printf("Timestamps are all identical\n");
			} else {
				printf("Timestamps are all in order\n");
			}
		} else {
			printf("Timestamps are not in order\n");
		}

		ser_file->ts_min = min_ts;
		ser_file->ts_max = max_ts;
		uint64_t diff_ts = (ser_file->ts_max - ser_file->ts_min) / 1000.0; // Now in units of 100 us
		if (diff_ts > 0) {
			// There is a positive time difference between first and last timestamps
			// We can calculate a frames per second value
			ser_file->fps = ((double) (ser_file->frame_count - 1) * 10000)
				/ (double) diff_ts;
		} else {
			// No valid frames per second value can be calculated
			ser_file->fps = -1.0;
		}
	} else {
		printf("No timestamps stored in the file.\n");
		ser_file->fps = -1.0;
	}
	return 0;
}

static int ser_fix_broken_file(struct ser_struct *ser_file) {
	int frame_count_calculated;
	off_t filesize = ser_file->filesize;

	siril_log_message(_("Trying to fix broken SER file...\n"));
	int frame_size = ser_file->image_width * ser_file->image_height;

	if (ser_file->color_id == SER_RGB || ser_file->color_id == SER_BGR) {
		frame_size *= 3;  // Color images have twice as many samples
	}

	if (ser_file->bit_pixel_depth > 8) {
		frame_size *= 2;  // Greater than 8-bit data has 2 bytes per pixel rather than one
	}

	filesize -= SER_HEADER_LEN;  // Remove header size from file size
	frame_count_calculated = filesize / frame_size;

	return frame_count_calculated;
}

static int ser_read_header(struct ser_struct *ser_file) {
	char header[SER_HEADER_LEN];

	if (!ser_file || ser_file->fd <= 0)
		return -1;

	/* Get file size */
	ser_file->filesize = lseek(ser_file->fd, 0, SEEK_END);
	if (ser_file->filesize == -1) {
		perror("seek");
		return -1;
	}
	lseek(ser_file->fd, 0, SEEK_SET);

	/* Read header (size of 178) */
	if (SER_HEADER_LEN != read(ser_file->fd, header, sizeof(header))) {
		perror("read");
		return -1;
	}

	// modify this to support big endian
	memcpy(&ser_file->lu_id, header + 14, 28);	// read all integers

	memcpy(&ser_file->date, header + 162, 8);
	memcpy(&ser_file->date_utc, header + 170, 8);

	// strings
	ser_file->file_id = strndup(header, 14);

	memcpy(ser_file->observer, header + 42, 40);
	memcpy(ser_file->instrument, header + 82, 40);
	memcpy(ser_file->telescope, header + 122, 40);

	/* internal representations of header data */
	if (ser_file->bit_pixel_depth <= 8)
		ser_file->byte_pixel_depth = SER_PIXEL_DEPTH_8;
	else ser_file->byte_pixel_depth = SER_PIXEL_DEPTH_16;

	if (ser_file->color_id == SER_RGB || ser_file->color_id == SER_BGR)
		ser_file->number_of_planes = 3;
	else
		ser_file->number_of_planes = 1;

/* In some cases, oacapture, firecapture, ... crash before writing frame_count
 * data. Here we try to get the calculated frame count which has not been written
 * in the header. Then we fix the SER file
 */
	if (ser_file->frame_count == 0) {
		ser_file->frame_count = ser_fix_broken_file(ser_file);

		if (ser_file->frame_count > 0) {
			if (ser_write_header(ser_file) == 0)
				siril_log_message(_("SER file has been fixed...\n"));
		}
	}

	ser_read_timestamp(ser_file);

	return 0;
}

static int ser_write_timestamp(struct ser_struct *ser_file) {
	int frame_size, i;

	if (ser_file->frame_count > 0) {
		if (ser_file->ts) {
			// Seek to start of timestamps
			frame_size = ser_file->image_width * ser_file->image_height
					* ser_file->number_of_planes;
			off_t offset = SER_HEADER_LEN
					+ (off_t) frame_size * (off_t) ser_file->byte_pixel_depth
							* (off_t) ser_file->frame_count;
			for (i = 0; i < ser_file->frame_count; i++) {
				if ((off_t) -1
						== lseek(ser_file->fd, offset + (i * 8), SEEK_SET)) {
					return -1;
				}
				char timestamp[8];

				memset(timestamp, 0, sizeof(timestamp));
				memcpy(timestamp, &ser_file->ts[i], 8);
				if (8 != write(ser_file->fd, timestamp, 8)) {
					perror("WriteTimetamps:");
					return -1;
				}
			}
		}
	}
	return 0;
}

static int ser_write_header(struct ser_struct *ser_file) {
	char header[SER_HEADER_LEN];

	if (!ser_file || ser_file->fd <= 0)
		return -1;
	if ((off_t) -1 == lseek(ser_file->fd, 0, SEEK_SET)) {
		perror("seek");
		return -1;
	}

	// modify this to support big endian
	memset(header, 0, sizeof(header));
	memcpy(header, ser_file->file_id, 14);
	memcpy(header + 14, &ser_file->lu_id, 28);
	memcpy(header + 42, ser_file->observer, 40);
	memcpy(header + 82, ser_file->instrument, 40);
	memcpy(header + 122, ser_file->telescope, 40);
	memcpy(header + 162, &ser_file->date, 8);
	memcpy(header + 170, &ser_file->date_utc, 8);

	if (sizeof(header) != write(ser_file->fd, header, sizeof(header))) {
		perror("write");
		return 1;
	}

	ser_write_timestamp(ser_file);

	return 0;
}

/* populate fields that are not already set in ser_create_file */
static void ser_write_header_from_fit(struct ser_struct *ser_file, fits *fit) {
	ser_file->image_width = fit->rx;
	ser_file->image_height = fit->ry;
	fprintf(stdout, "setting SER image size as %dx%d\n", fit->rx, fit->ry);
	if (fit->naxes[2] == 1) {
		// maybe it's read as CFA...
		ser_file->color_id = SER_MONO;
	} else if (fit->naxes[2] == 3) {
		ser_file->color_id = SER_RGB;
	}
	if (ser_file->color_id == SER_RGB || ser_file->color_id == SER_BGR)
		ser_file->number_of_planes = 3;
	else ser_file->number_of_planes = 1;

	if (fit->bitpix == BYTE_IMG) {
		ser_file->byte_pixel_depth = SER_PIXEL_DEPTH_8;
		ser_file->bit_pixel_depth = 8;
	} else if (fit->bitpix == USHORT_IMG || fit->bitpix == SHORT_IMG) {
		ser_file->byte_pixel_depth = SER_PIXEL_DEPTH_16;
		ser_file->bit_pixel_depth = 16;
	} else {
		siril_log_message(_("Writing to SER files from larger than 16-bit FITS images is not yet implemented\n"));
	}
	if (fit->instrume[0] != 0) {
		memset(ser_file->instrument, 0, 40);
		memcpy(ser_file->instrument, fit->instrume, 40);
	}
	if (fit->observer[0] != 0) {
		memset(ser_file->observer, 0, 40);
		memcpy(ser_file->observer, fit->observer, 40);
	}
	if (fit->instrume[0] != 0) {
		memset(ser_file->telescope, 0, 40);
		memcpy(ser_file->telescope, fit->telescop, 40);
	}
	int ret = FITS_date_key_to_Unix_time(fit->date_obs, &ser_file->date_utc, &ser_file->date);
	if (ret == -1)
		FITS_date_key_to_Unix_time(fit->date, &ser_file->date_utc, &ser_file->date);
}

static int retrieveSERBayerPattern(ser_color pattern) {

	switch (pattern) {
	case SER_BAYER_RGGB:
		return BAYER_FILTER_RGGB;
		break;
	case SER_BAYER_BGGR:
		return BAYER_FILTER_BGGR;
		break;
	case SER_BAYER_GBRG:
		return BAYER_FILTER_GBRG;
		break;
	case SER_BAYER_GRBG:
		return BAYER_FILTER_GRBG;
		break;
	default:
		return BAYER_FILTER_NONE;
	}
}

/* once a buffer (data) has been acquired from the file, with frame_size pixels
 * read in it, depending on ser_file's endianess and pixel depth, data is
 * reorganized to match Siril's data format . */
void ser_manage_endianess_and_depth(struct ser_struct *ser_file, WORD *data, int frame_size) {
	WORD pixel;
	int i;
	if (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8) {
		// inline conversion to 16 bit
		for (i = frame_size - 1; i >= 0; i--)
			data[i] = (WORD) (((BYTE*)data)[i]);
	} else if (ser_file->little_endian == SER_BIG_ENDIAN) {
		// inline conversion
		for (i = frame_size - 1; i >= 0; i--) {
			pixel = data[i];
			pixel = (pixel >> 8) | (pixel << 8);
			data[i] = pixel;
		}
	}
}

/*
 * Public functions
 */

void ser_display_info(struct ser_struct *ser_file) {
	char *color = convert_color_id_to_char(ser_file->color_id);

	fprintf(stdout, "=========== SER file info ==============\n");
	fprintf(stdout, "file id: %s\n", ser_file->file_id);
	fprintf(stdout, "lu id: %d\n", ser_file->lu_id);
	fprintf(stdout, "little endian: %d\n", ser_file->little_endian);
	fprintf(stdout, "sensor type: %s\n", color);
	fprintf(stdout, "image size: %d x %d (%d bits)\n", ser_file->image_width,
			ser_file->image_height, ser_file->bit_pixel_depth);
	fprintf(stdout, "frame count: %u\n", ser_file->frame_count);
	fprintf(stdout, "observer: %.40s\n", ser_file->observer);
	fprintf(stdout, "instrument: %.40s\n", ser_file->instrument);
	fprintf(stdout, "telescope: %.40s\n", ser_file->telescope);
	display_date(ser_file->date, "local time: ");
	display_date(ser_file->date_utc, "UTC time: ");
	fprintf(stdout, "fps: %.3lf\n", ser_file->fps);
	fprintf(stdout, "========================================\n");
}

int ser_write_and_close(struct ser_struct *ser_file) {
	ser_write_header(ser_file);	// writes the header
	return ser_close_file(ser_file);// closes, frees and zeroes
}

/* ser_file must be allocated */
int ser_create_file(const char *filename, struct ser_struct *ser_file,
		gboolean overwrite, struct ser_struct *copy_from) {
	if (overwrite)
		unlink(filename);
	if ((ser_file->fd = open(filename, O_CREAT | O_RDWR,
			S_IWRITE | S_IREAD)) == -1) {
		perror("open SER file for creation");
		return 1;
	}

	if (copy_from) {
		memcpy(&ser_file->lu_id, &copy_from->lu_id, 28);
		memcpy(&ser_file->date, &copy_from->date, 8);
		memcpy(&ser_file->date_utc, &copy_from->date_utc, 8);
		ser_file->file_id = strdup(copy_from->file_id);
		memcpy(ser_file->observer, copy_from->observer, 40);
		memcpy(ser_file->instrument, copy_from->instrument, 40);
		memcpy(ser_file->telescope, copy_from->telescope, 40);
		ser_file->byte_pixel_depth = copy_from->byte_pixel_depth;
		ser_file->number_of_planes = copy_from->number_of_planes;

		int i;
		if (copy_from->ts) {
			ser_file->ts = calloc(8, copy_from->frame_count);
			for (i = 0; i < copy_from->frame_count; i++)
				ser_file->ts[i] = copy_from->ts[i];
		} else {
			ser_file->ts = NULL;
			ser_file->fps = -1.0;
		}
		/* we write the header now, but it should be written again
		 * before closing in case the number of the image in the new
		 * SER changes from the copied SER */
		ser_write_header(ser_file);
	} else {	// new SER
		ser_file->file_id = strdup("Made by Siril");
		ser_file->lu_id = 0;
		ser_file->little_endian = SER_LITTLE_ENDIAN; // what will it do on big endian machine?
		memset(ser_file->observer, 0, 40);
		memset(ser_file->instrument, 0, 40);
		memset(ser_file->telescope, 0, 40);
		ser_file->ts = NULL;
		memset(&ser_file->date, 0, 8);
		memset(&ser_file->date_utc, 0, 8);
		ser_file->number_of_planes = 0;	// used as an indicator of new SER

		/* next operation should be ser_write_frame_from_fit, which writes with no
		 * seek and expects to be after the header */
		if ((off_t) -1 == lseek(ser_file->fd, SER_HEADER_LEN, SEEK_SET)) {
			perror("seek");
			return -1;
		}
	}
	ser_file->filename = strdup(filename);
	ser_file->frame_count = 0;	// incremented on image add
#ifdef _OPENMP
	omp_init_lock(&ser_file->fd_lock);
#endif
	siril_log_message(_("Created SER file %s\n"), filename);
	return 0;
}

int ser_open_file(char *filename, struct ser_struct *ser_file) {
	if (ser_file->fd > 0) {
		fprintf(stderr, "SER: file already opened, or badly closed\n");
		return -1;
	}
	ser_file->fd = open(filename, O_RDWR); // now we can fix broken file, so not O_RDONLY anymore
	if (ser_file->fd == -1) {
		perror("SER file open");
		return -1;
	}
	if (ser_read_header(ser_file)) {
		fprintf(stderr, "SER: reading header failed, closing file %s\n",
				filename);
		ser_close_file(ser_file);
		return -1;
	}
	ser_file->filename = strdup(filename);

#ifdef _OPENMP
	omp_init_lock(&ser_file->fd_lock);
#endif
	return 0;
}

int ser_close_file(struct ser_struct *ser_file) {
	int retval = 0;
	if (!ser_file)
		return retval;
	if (ser_file->fd > 0) {
		retval = close(ser_file->fd);
		ser_file->fd = -1;
	}
	if (ser_file->file_id)
		free(ser_file->file_id);
	if (ser_file->ts)
		free(ser_file->ts);
	if (ser_file->filename)
		free(ser_file->filename);
#ifdef _OPENMP
	omp_destroy_lock(&ser_file->fd_lock);
#endif
	ser_init_struct(ser_file);
	return retval;
}

void ser_init_struct(struct ser_struct *ser_file) {
	memset(ser_file, 0, sizeof(struct ser_struct));
}

/* frame number starts at 0 */
int ser_read_frame(struct ser_struct *ser_file, int frame_no, fits *fit) {
	int retval, frame_size, i, j, swap = 0;
	off_t offset;
	WORD *olddata, *tmp;
	if (!ser_file || ser_file->fd <= 0 || !ser_file->number_of_planes ||
			!fit || frame_no < 0 || frame_no >= ser_file->frame_count)
		return -1;
	frame_size = ser_file->image_width * ser_file->image_height
			* ser_file->number_of_planes;
	olddata = fit->data;
	if ((fit->data = realloc(fit->data, frame_size * sizeof(WORD))) == NULL) {
		fprintf(stderr, "ser_read: error realloc %s %d\n", ser_file->filename,
				frame_size);
		if (olddata)
			free(olddata);
		return -1;
	}

	offset = SER_HEADER_LEN + (off_t)frame_size *
		(off_t)ser_file->byte_pixel_depth * (off_t)frame_no;
	/*fprintf(stdout, "offset is %lu (frame %d, %d pixels, %d-byte)\n", offset,
	 frame_no, frame_size, ser_file->pixel_bytedepth);*/
#ifdef _OPENMP
	omp_set_lock(&ser_file->fd_lock);
#endif
	if ((off_t) -1 == lseek(ser_file->fd, offset, SEEK_SET)) {
#ifdef _OPENMP
		omp_unset_lock(&ser_file->fd_lock);
#endif
		return -1;
	}
	retval = read(ser_file->fd, fit->data, frame_size * ser_file->byte_pixel_depth);
#ifdef _OPENMP
	omp_unset_lock(&ser_file->fd_lock);
#endif
	if (retval != frame_size * ser_file->byte_pixel_depth)
		return -1;

	ser_manage_endianess_and_depth(ser_file, fit->data, frame_size);

	fit->bitpix = (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8) ? BYTE_IMG : USHORT_IMG;

	/* If the user checks the SER CFA box, the video is opened in B&W
	 * RGB and BGR are not coming from raw data. In consequence CFA does
	 * not exist for these kind of cam */
	ser_color type_ser = ser_file->color_id;
	if (!com.debayer.open_debayer && type_ser != SER_RGB && type_ser != SER_BGR)
		type_ser = SER_MONO;

	switch (type_ser) {
	case SER_MONO:
		fit->naxis = 2;
		fit->naxes[0] = fit->rx = ser_file->image_width;
		fit->naxes[1] = fit->ry = ser_file->image_height;
		fit->naxes[2] = 1;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data;
		fit->pdata[BLAYER] = fit->data;
		break;
	case SER_BAYER_RGGB:
	case SER_BAYER_BGGR:
	case SER_BAYER_GBRG:
	case SER_BAYER_GRBG:
		fit->naxes[0] = fit->rx = ser_file->image_width;
		fit->naxes[1] = fit->ry = ser_file->image_height;
		fit->naxes[2] = 3;
		/* Get Bayer informations from header if available */
		sensor_pattern sensortmp;
		sensortmp = com.debayer.bayer_pattern;
		if (com.debayer.use_bayer_header) {
			sensor_pattern bayer;
			bayer = retrieveSERBayerPattern(type_ser);
			if (bayer != com.debayer.bayer_pattern) {
				if (bayer == BAYER_FILTER_NONE) {
					siril_log_color_message(_("No Bayer pattern found in the header file.\n"), "red");
				}
				else {
					siril_log_color_message(_("Bayer pattern found in header (%s) is different"
							" from Bayer pattern in settings (%s). Overriding settings.\n"),
							"red", filter_pattern[bayer], filter_pattern[com.debayer.bayer_pattern]);
					com.debayer.bayer_pattern = bayer;
				}
			}
		}
		debayer(fit, com.debayer.bayer_inter);
		com.debayer.bayer_pattern = sensortmp;
		break;
	case SER_BGR:
		swap = 2;
		/* no break */
	case SER_RGB:
		tmp = malloc(frame_size * sizeof(WORD));
		memcpy(tmp, fit->data, sizeof(WORD) * frame_size);
		fit->naxes[0] = fit->rx = ser_file->image_width;
		fit->naxes[1] = fit->ry = ser_file->image_height;
		fit->naxes[2] = 3;
		fit->naxis = 3;
		fit->pdata[RLAYER] = fit->data;
		fit->pdata[GLAYER] = fit->data + fit->rx * fit->ry;
		fit->pdata[BLAYER] = fit->data + fit->rx * fit->ry * 2;
		for (i = 0, j = 0; j < fit->rx * fit->ry; i += 3, j++) {
			fit->pdata[0 + swap][j] = tmp[i + RLAYER];
			fit->pdata[1][j] = tmp[i + GLAYER];
			fit->pdata[2 - swap][j] = tmp[i + BLAYER];
		}
		free(tmp);
		break;
	case SER_BAYER_CYYM:
	case SER_BAYER_YCMY:
	case SER_BAYER_YMCY:
	case SER_BAYER_MYYC:
	default:
		siril_log_message(_("This type of Bayer pattern is not handled yet.\n"));
		return -1;
	}
	fits_flip_top_to_bottom(fit);
	return 0;
}

/* read an area of an image in an opened SER sequence */
int ser_read_opened_partial(struct ser_struct *ser_file, int layer,
		int frame_no, WORD *buffer, const rectangle *area) {
	off_t offset;
	int frame_size, read_size, retval, xoffset, yoffset, x, y, color_offset;
	ser_color type_ser;
	WORD *rawbuf, *demosaiced_buf, *rgb_buf;
	rectangle debayer_area, image_area;
	sensor_pattern sensortmp;

	if (!ser_file || ser_file->fd <= 0 || frame_no < 0
			|| frame_no >= ser_file->frame_count)
		return -1;
	frame_size = ser_file->image_width * ser_file->image_height *
		ser_file->number_of_planes * ser_file->byte_pixel_depth;

	/* If the user checks the SER CFA box, the video is opened in B&W
	 * RGB and BGR are not coming from raw data. In consequence CFA does
	 * not exist for these kind of cam */
	type_ser = ser_file->color_id;
	if (!com.debayer.open_debayer && type_ser != SER_RGB && type_ser !=
			SER_BGR)
		type_ser = SER_MONO;

	switch (type_ser) {
	case SER_MONO:
		offset = SER_HEADER_LEN + (off_t)frame_size * (off_t)frame_no +	// requested frame
			(off_t)(area->y * ser_file->image_width + area->x)
						* ser_file->byte_pixel_depth;	// requested area
#ifdef _OPENMP
		omp_set_lock(&ser_file->fd_lock);
#endif
		if ((off_t) -1 == lseek(ser_file->fd, offset, SEEK_SET)) {
#ifdef _OPENMP
			omp_unset_lock(&ser_file->fd_lock);
#endif
			return -1;
		}
		read_size = area->w * area->h * ser_file->byte_pixel_depth;
		retval = read(ser_file->fd, buffer, read_size);
#ifdef _OPENMP
		omp_unset_lock(&ser_file->fd_lock);
#endif
		if (retval != read_size)
			return -1;

		ser_manage_endianess_and_depth(ser_file, buffer, area->w * area->h);
		break;

	case SER_BAYER_RGGB:
	case SER_BAYER_BGGR:
	case SER_BAYER_GBRG:
	case SER_BAYER_GRBG:
		/* SER v2: RGB images obtained from demosaicing.
		 * Original is monochrome, we demosaic it in an area slightly larger than the
		 * requested area, giving 3 channels in form of RGBRGBRGB buffers, and finally
		 * we extract one of the three channels and crop it to the requested area. */

		/* Get Bayer informations from header if available */
		sensortmp = com.debayer.bayer_pattern;
		if (com.debayer.use_bayer_header) {
			sensor_pattern bayer;
			bayer = retrieveSERBayerPattern(type_ser);
			if (bayer != com.debayer.bayer_pattern) {
				if (bayer == BAYER_FILTER_NONE) {
					siril_log_color_message(_("No Bayer pattern found in the header file.\n"), "red");
				}
				else {
					siril_log_color_message(_("Bayer pattern found in header (%s) is different"
							" from Bayer pattern in settings (%s). Overriding settings.\n"),
							"red", filter_pattern[bayer], filter_pattern[com.debayer.bayer_pattern]);
					com.debayer.bayer_pattern = bayer;
				}
			}
		}
		if (layer < 0 || layer >= 3) {
			siril_log_message(_("For a demosaiced image, layer has to be R, G or B (0 to 2).\n"));
			return -1;
		}

		image_area = (rectangle) { .x = 0, .y = 0,
			.w = ser_file->image_width, .h = ser_file->image_height };
		get_debayer_area(area, &debayer_area, &image_area, &xoffset, &yoffset);

		offset = SER_HEADER_LEN + frame_size * frame_no +	// requested frame
				(debayer_area.y * ser_file->image_width + debayer_area.x)
						* ser_file->byte_pixel_depth;	// requested area
#ifdef _OPENMP
		omp_set_lock(&ser_file->fd_lock);
#endif
		if ((off_t) -1 == lseek(ser_file->fd, offset, SEEK_SET)) {
#ifdef _OPENMP
			omp_unset_lock(&ser_file->fd_lock);
#endif
			return -1;
		}

		// allocating a buffer for WORD because it's going to be converted in-place
		rawbuf = malloc(debayer_area.w * debayer_area.h * sizeof(WORD));
		if (!rawbuf) {
#ifdef _OPENMP
			omp_unset_lock(&ser_file->fd_lock);
#endif
			siril_log_message(_("Out of memory - aborting\n"));
			return -1;
		}
		read_size = debayer_area.w * debayer_area.h * ser_file->byte_pixel_depth;
		retval = read(ser_file->fd, rawbuf, read_size);
#ifdef _OPENMP
		omp_unset_lock(&ser_file->fd_lock);
#endif
		if (retval != read_size) {
			free(rawbuf);
			perror("read");
			return -1;
		}

		ser_manage_endianess_and_depth(ser_file, rawbuf, debayer_area.w * debayer_area.h);

		demosaiced_buf = debayer_buffer(rawbuf, &debayer_area.w,
				&debayer_area.h, com.debayer.bayer_inter,
				com.debayer.bayer_pattern);
		free(rawbuf);
		if (demosaiced_buf == NULL) {
			return -1;
		}

		/* area is the destination area.
		 * debayer_area is the demosaiced buf area.
		 * xoffset and yoffset are the x,y offsets of area in the debayer area.
		 */
		for (y = 0; y < area->h; y++) {
			for (x = 0; x < area->w; x++) {
				buffer[y*area->w + x] = demosaiced_buf[(yoffset+y)*debayer_area.w*3 + xoffset+x*3 + layer]; 
			}
		}

		free(demosaiced_buf);
		com.debayer.bayer_pattern = sensortmp;
		break;
	case SER_BGR:
	case SER_RGB:
		assert(ser_file->number_of_planes == 3);

		offset = SER_HEADER_LEN + frame_size * frame_no +	// requested frame
			(area->y * ser_file->image_width + area->x) *
			ser_file->byte_pixel_depth * 3;	// requested area
#ifdef _OPENMP
		omp_set_lock(&ser_file->fd_lock);
#endif
		if ((off_t) -1 == lseek(ser_file->fd, offset, SEEK_SET)) {
#ifdef _OPENMP
			omp_unset_lock(&ser_file->fd_lock);
#endif
			return -1;
		}

		read_size = area->w * area->h * ser_file->byte_pixel_depth * 3;
		// allocating a buffer for WORD because it's going to be converted in-place
		rgb_buf = malloc(area->w * area->h * 3 * sizeof(WORD));
		if (!rgb_buf) {
#ifdef _OPENMP
			omp_unset_lock(&ser_file->fd_lock);
#endif
			siril_log_message(_("Out of memory - aborting\n"));
			return -1;
		}
		retval = read(ser_file->fd, rgb_buf, read_size);
#ifdef _OPENMP
		omp_unset_lock(&ser_file->fd_lock);
#endif
		if (retval != read_size) {
			free(rgb_buf);
			perror("read");
			return -1;
		}

		ser_manage_endianess_and_depth(ser_file, rgb_buf, area->w * area->h * 3);

		color_offset = layer;
		if (type_ser == SER_BGR) {
			color_offset = 2 - layer;
		}

		for (y = 0; y < area->h; y++) {
			for (x = 0; x < area->w; x++) {
				buffer[y*area->w + x] = rgb_buf[y*area->w*3 + x*3 + color_offset]; 
			}
		}
		free(rgb_buf);
		break;
	default:
		siril_log_message(_("This type of Bayer pattern is not handled yet.\n"));
		return -1;
	}

	return 0;
}

int ser_write_frame_from_fit(struct ser_struct *ser_file, fits *fit, int frame_no) {
	int frame_size, pixel, plane, dest;
	int ret, retval = 0;
	off_t offset;
	BYTE *data8 = NULL;			// for 8-bit files
	WORD *data16 = NULL;		// for 16-bit files

	if (!ser_file || ser_file->fd <= 0 || !fit)
		return -1;
	if (ser_file->number_of_planes == 0) {
		// adding first frame of a new sequence, use it to populate the header
		ser_write_header_from_fit(ser_file, fit);
	}
	if (fit->rx != ser_file->image_width || fit->ry != ser_file->image_height) {
		siril_log_message(_("Trying to add an image of different size in a SER\n"));
		return 1;
	}

	fits_flip_top_to_bottom(fit);
	frame_size = ser_file->image_width * ser_file->image_height *
		ser_file->number_of_planes;

	offset = SER_HEADER_LEN	+ (off_t) frame_size *
			(off_t) ser_file->byte_pixel_depth * (off_t) frame_no;

	if (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8)
		data8 = malloc(frame_size * ser_file->byte_pixel_depth);
	else
		data16 = malloc(frame_size * ser_file->byte_pixel_depth);

	for (plane = 0; plane < ser_file->number_of_planes; plane++) {
		dest = plane;
		for (pixel = 0; pixel < ser_file->image_width * ser_file->image_height;
				pixel++) {
			if (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8)
				data8[dest] = (BYTE)(fit->pdata[plane][pixel]);
			else {
				if (ser_file->little_endian == SER_BIG_ENDIAN)
					data16[dest] = (fit->pdata[plane][pixel] >> 8 | fit->pdata[plane][pixel] << 8);
				else
					data16[dest] = fit->pdata[plane][pixel];
			}
			dest += ser_file->number_of_planes;
		}
	}

#ifdef _OPENMP
	omp_set_lock(&ser_file->fd_lock);
#endif
	if ((off_t)-1 == lseek(ser_file->fd, offset, SEEK_SET)) {
#ifdef _OPENMP
		omp_unset_lock(&ser_file->fd_lock);
#endif
		perror("seek");
		retval = -1;
		goto free_and_quit;
	}

	if (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8) {
		ret = write(ser_file->fd, data8, frame_size * ser_file->byte_pixel_depth);
#ifdef _OPENMP
		omp_unset_lock(&ser_file->fd_lock);
#endif
		if (ret != frame_size * ser_file->byte_pixel_depth) {
			perror("write image in SER");
			retval = 1;
			goto free_and_quit;
		}
	} else {
		ret = write(ser_file->fd, data16, frame_size * ser_file->byte_pixel_depth);
#ifdef _OPENMP
		omp_unset_lock(&ser_file->fd_lock);
#endif
		if (ret != frame_size * ser_file->byte_pixel_depth) {
			perror("write image in SER");
			retval = 1;
			goto free_and_quit;
		}
	}

#ifdef _OPENMP
#pragma omp atomic
#endif
	ser_file->frame_count++;
	retval = 0;

free_and_quit:
	if (data8) free(data8);
	if (data16) free(data16);
	return retval;
}
