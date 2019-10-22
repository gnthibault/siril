/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2019 team free-astro (see more in AUTHORS file)
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
#ifdef _WIN32
#include <io.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "algos/demosaicing.h"
#include "io/ser.h"

static gboolean user_warned = FALSE;

/* 62135596800 sec from year 0001 to 01 janv. 1970 00:00:00 GMT */
static const uint64_t epochTicks = 621355968000000000UL;
static const uint64_t ticksPerSecond = 10000000;

static int ser_write_header(struct ser_struct *ser_file);

/* Given a SER timestamp, return a char string representation
 * MUST be freed
 */
static char *ser_timestamp(uint64_t timestamp) {
	char *str = malloc(64);
	uint64_t t1970_ms = (timestamp - epochTicks) / 10000;
	time_t secs = t1970_ms / 1000;
	int ms = t1970_ms % 1000;
	struct tm *t;
#ifdef HAVE_GMTIME_R
	struct tm t_;
#endif

#ifdef _WIN32
	t = gmtime (&secs);
#else
#ifdef HAVE_GMTIME_R
	t = gmtime_r (&secs, &t_);
#else
	t = gmtime(&secs);
#endif /* HAVE_GMTIME_R */
#endif /* _WIN32 */

	/* If the gmtime() call has failed, "secs" is too big. */
	if (t == NULL) {
		free(str);
		return NULL;
	}

	sprintf(str, "%04d-%02d-%02dT%02d:%02d:%02d.%03d", t->tm_year + 1900,
			t->tm_mon + 1, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec, ms);

	return str;
}

/* Output SER timestamp */
static int display_date(uint64_t timestamp, char *txt) {
	if (timestamp == 0)
		return -1;

	char *str = ser_timestamp(timestamp);
	if (str) {
		fprintf(stdout, "%s%s\n", txt, str);
		free(str);
	}
	return 0;
}

static time_t mktime_utc(struct tm *tm) {
	time_t retval;

#ifndef HAVE_TIMEGM
	static const gint days_before[] = { 0, 31, 59, 90, 120, 151, 181, 212, 243,
			273, 304, 334 };
#endif

#ifndef HAVE_TIMEGM
	if (tm->tm_mon < 0 || tm->tm_mon > 11)
		return (time_t) -1;

	retval = (tm->tm_year - 70) * 365;
	retval += (tm->tm_year - 68) / 4;
	retval += days_before[tm->tm_mon] + tm->tm_mday - 1;

	if (tm->tm_year % 4 == 0 && tm->tm_mon < 2)
		retval -= 1;

	retval = ((((retval * 24) + tm->tm_hour) * 60) + tm->tm_min) * 60
			+ tm->tm_sec;
#else
	retval = timegm (tm);
#endif /* !HAVE_TIMEGM */

	return retval;
}

/* Convert FITS keyword DATE in a UNIX time format
 * DATE match this pattern: 1900-01-01T00:00:00
 */
static int FITS_date_key_to_Unix_time(char *date, uint64_t *utc,
		uint64_t *local) {
	struct tm timeinfo = { };
	time_t ut, t;
	int year = 0, month = 0, day = 0, hour = 0, min = 0, ms = 0;
	float sec = 0.0;

	if (date[0] == '\0')
		return -1;

	sscanf(date, "%04d-%02d-%02dT%02d:%02d:%f", &year, &month, &day, &hour,
			&min, &sec);

	timeinfo.tm_year = year - 1900;
	timeinfo.tm_mon = month - 1;
	timeinfo.tm_mday = day;
	timeinfo.tm_hour = hour;
	timeinfo.tm_min = min;
	timeinfo.tm_sec = (int) sec;
	ms = ((int) (sec * 1000) % 1000);

	// Hopefully these are not needed
	timeinfo.tm_wday = 0;
	timeinfo.tm_yday = 0;
	timeinfo.tm_isdst = -1;

	/* get UTC time from timeinfo* */
	ut = mktime_utc(&timeinfo);
	ut *= ticksPerSecond;
	ut += epochTicks;
	ut += ms * 10000;
	*utc = (uint64_t) ut;

	/* get local time from timeinfo* */
	t = mktime(&timeinfo);
	t *= ticksPerSecond;
	t += epochTicks;
	t += ms * 10000;
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

/* reads timestamps from the trailer of the file and stores them in ser_file->ts */
static int ser_read_timestamp(struct ser_struct *ser_file) {
	int i;
	gboolean timestamps_in_order = TRUE;
	uint64_t previous_ts = 0L;
	int64_t frame_size;

	ser_file->fps = -1.0;	// will be calculated from the timestamps

	if (!ser_file->frame_count || ser_file->image_width <= 0 ||
			ser_file->image_height <= 0 || ser_file->byte_pixel_depth <= 0 ||
			!ser_file->number_of_planes)
		return 0;

	frame_size = ser_file->image_width *
		ser_file->image_height * ser_file->number_of_planes;
	int64_t offset = SER_HEADER_LEN + frame_size *
		(int64_t)ser_file->byte_pixel_depth * (int64_t)ser_file->frame_count;
	/* Check if file is large enough to have timestamps */
	if (ser_file->filesize >= offset + (8 * ser_file->frame_count)) {
		ser_file->ts = calloc(8, ser_file->frame_count);
		ser_file->ts_alloc = ser_file->frame_count;

		// Seek to start of timestamps
		for (i = 0; i < ser_file->frame_count; i++) {
			if ((int64_t)-1 == fseek64(ser_file->file, offset+(i*8), SEEK_SET))
				return -1;

			if (8 != fread(&ser_file->ts[i], 1, 8, ser_file->file))
				return 0;
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
			if (min_ts == max_ts)
				fprintf(stdout, _("Warning: timestamps in the SER sequence are all identical.\n"));
			else fprintf(stdout, _("Timestamps in the SER sequence are correctly ordered.\n"));
		} else {
			fprintf(stdout, _("Warning: timestamps in the SER sequence are not in the correct order.\n"));
		}

		ser_file->ts_min = min_ts;
		ser_file->ts_max = max_ts;
		double diff_ts = (ser_file->ts_max - ser_file->ts_min) / 1000.0;
		// diff_ts now in units of 100 us or ten thousandths of a second
		if (diff_ts > 0.0) {
			// There is a positive time difference between first and last
			// timestamps, we can calculate a frames per second value
			ser_file->fps = (ser_file->frame_count - 1) * 10000.0 / diff_ts;
		}
	} else {
		fprintf(stdout, _("Warning: no timestamps stored in the SER sequence.\n"));
	}
	return 0;
}

static int ser_recompute_frame_count(struct ser_struct *ser_file) {
	int frame_count_calculated;
	int64_t filesize = ser_file->filesize;

	siril_log_message(_("Trying to fix broken SER file...\n"));
	int64_t frame_size = ser_file->image_width * ser_file->image_height;
	if (frame_size == 0)
		return 0;

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

	if (!ser_file || ser_file->file == NULL)
		return -1;

	/* Get file size */
	fseek64(ser_file->file, 0, SEEK_END);
	ser_file->filesize = ftell64(ser_file->file);
	fseek64(ser_file->file, 0, SEEK_SET);
	if (ser_file->filesize == -1) {
		perror("seek");
		return -1;
	}

	/* Read header (size of 178) */
	if (SER_HEADER_LEN != fread(header, 1, sizeof header, ser_file->file)) {
		perror("fread");
		return -1;
	}

	// modify this to support big endian
	memcpy(&ser_file->lu_id, header + 14, 28);	// read all integers

	memcpy(&ser_file->date, header + 162, 8);
	memcpy(&ser_file->date_utc, header + 170, 8);

	// strings
	ser_file->file_id = g_strndup(header, 14);

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

	/* In some cases, oacapture, firecapture, ... crash before writing
	 * frame_count data. Here we try to get the calculated frame count
	 * which has not been written in the header. Then we fix the SER file
	 */
	if (ser_file->frame_count == 0) {
		ser_file->frame_count = ser_recompute_frame_count(ser_file);

		if (ser_file->frame_count > 0) {
			if (ser_write_header(ser_file) == 0)
				siril_log_message(_("SER file has been fixed...\n"));
		}
	}

	ser_read_timestamp(ser_file);

	return 0;
}

static int ser_write_timestamps(struct ser_struct *ser_file) {
	int i;
	int64_t frame_size;

	if (!ser_file->frame_count || ser_file->image_width <= 0 ||
			ser_file->image_height <= 0 || ser_file->byte_pixel_depth <= 0 ||
			!ser_file->number_of_planes)
		return -1;

	if (ser_file->ts) {
		// Seek to start of timestamps
		frame_size = ser_file->image_width * ser_file->image_height
			* ser_file->number_of_planes;
		int64_t offset = SER_HEADER_LEN + frame_size * 
			(int64_t)ser_file->byte_pixel_depth * (int64_t)ser_file->frame_count;

		for (i = 0; i < ser_file->frame_count; i++) {
			if (i >= ser_file->ts_alloc)
				break;
			if ((int64_t)-1 == fseek64(ser_file->file, offset+(i*8), SEEK_SET)) {
				return -1;
			}
			if (8 != fwrite(&ser_file->ts[i], 1, 8, ser_file->file)) {
				perror("write timestamps:");
				return -1;
			}
		}
	}
	return 0;
}

/* (over)write the header of the opened file on the disk */
static int ser_write_header(struct ser_struct *ser_file) {
	char header[SER_HEADER_LEN];

	if (!ser_file || ser_file->file == NULL)
		return -1;
	if ((int64_t) -1 == fseek64(ser_file->file, 0, SEEK_SET)) {
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

	if (sizeof(header) != fwrite(header, 1, sizeof(header), ser_file->file)) {
		perror("write");
		return 1;
	}
	return 0;
}

/* populate fields that are not already set in ser_create_file */
static void ser_write_header_from_fit(struct ser_struct *ser_file, fits *fit) {
	ser_file->image_width = fit->rx;
	ser_file->image_height = fit->ry;
	fprintf(stdout, "setting SER image size as %dx%d\n", fit->rx, fit->ry);
	// already managed during creation for monochrome formats
	if (fit->naxes[2] == 3) {
		ser_file->color_id = SER_RGB;
	}
	if (ser_file->color_id == SER_RGB)
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

	if (FITS_date_key_to_Unix_time(fit->date_obs, &ser_file->date_utc, &ser_file->date) == -1)
		FITS_date_key_to_Unix_time(fit->date, &ser_file->date_utc, &ser_file->date);
}

static int get_SER_Bayer_Pattern(ser_color pattern) {
	switch (pattern) {
	case SER_BAYER_RGGB:
		return BAYER_FILTER_RGGB;
	case SER_BAYER_BGGR:
		return BAYER_FILTER_BGGR;
	case SER_BAYER_GBRG:
		return BAYER_FILTER_GBRG;
	case SER_BAYER_GRBG:
		return BAYER_FILTER_GRBG;
	default:
		return BAYER_FILTER_NONE;
	}
}

/* once a buffer (data) has been acquired from the file, with frame_size pixels
 * read in it, depending on ser_file's endianess and pixel depth, data is
 * reorganized to match Siril's data format . */
static void ser_manage_endianess_and_depth(struct ser_struct *ser_file,
		WORD *data, int64_t frame_size) {
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

static int ser_alloc_ts(struct ser_struct *ser_file, int frame_no) {
	int retval = 0;
#ifdef _OPENMP
	omp_set_lock(&ser_file->ts_lock);
#endif
	if (ser_file->ts_alloc <= frame_no) {
		uint64_t *new = realloc(ser_file->ts, (frame_no + 1) * 2 * sizeof(uint64_t));
		if (!new) {
			PRINT_ALLOC_ERR;
			retval = 1;
		} else {
			ser_file->ts = new;
			ser_file->ts_alloc = (frame_no + 1) * 2;
		}
	}
#ifdef _OPENMP
	omp_unset_lock(&ser_file->ts_lock);
#endif
	return retval;
}

/*
 * Public functions
 */

gboolean ser_is_cfa(struct ser_struct *ser_file) {
	return ser_file && (ser_file->color_id == SER_BAYER_RGGB || 
			ser_file->color_id == SER_BAYER_GRBG || 
			ser_file->color_id == SER_BAYER_GBRG || 
			ser_file->color_id == SER_BAYER_BGGR); 
	// SER_BAYER_CYYM SER_BAYER_YCMY SER_BAYER_YMCY SER_BAYER_MYYC are not
	// supported yet so returning false for them here is good
}

/* set the timestamps of the ser_file using a list of timestamps in string form */
void ser_convertTimeStamp(struct ser_struct *ser_file, GSList *timestamp) {
	int i = 0;
	if (ser_file->ts)
		free(ser_file->ts);
	ser_file->ts = calloc(8, ser_file->frame_count);
	ser_file->ts_alloc = ser_file->frame_count;

	GSList *t = timestamp;
	while (t && i < ser_file->frame_count) {
		uint64_t utc, local;
		FITS_date_key_to_Unix_time(t->data, &utc, &local);
		t = t->next;
		memcpy(&ser_file->ts[i], &utc, 8);
		i++;
	}
}

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
	if (ser_file == NULL) return -1;
	if (!ser_file->frame_count) {
		siril_log_color_message(_("The SER sequence is being created with no image in it.\n"), "red");
		char *filename = ser_file->filename;
		ser_file->filename = NULL;
		ser_close_file(ser_file);// closes, frees and zeroes
		g_unlink(filename);
		return -1;
	}
	ser_write_header(ser_file);	// writes the header
	ser_write_timestamps(ser_file);	// writes the trailer
	return ser_close_file(ser_file);// closes, frees and zeroes
}

/* calling ser_write_frame_from_fit() with image indices that do not cover a
 * contiguous range will pose some problems since there will be holes in images
 * data and the indices of frames in the file will be incorrect. To solve this
 * and still allow a parallel processing that can safely fail to be done with
 * SER output, we compact the frames that have been identified as failed in the
 * processing. The provided array contains booleans that inform about the
 * success of the processing and the presence of the frame with a given index
 * in the file. nb_frames is the size of this array and the last index + 1 of
 * the frames added to the file.
 * This function is not thread-safe. */
int ser_compact_file(struct ser_struct *ser_file, unsigned char *successful_frames, int nb_frames) {
	int64_t offseti, offsetj, frame_size;
	int i, j;
	unsigned char *buffer = NULL;
	frame_size = ser_file->image_width * ser_file->image_height *
		ser_file->number_of_planes * ser_file->byte_pixel_depth;

	// frame_count should be fine because it's incremented only when adding
	// one, but the real number of images for the file size if nb_frames
	for (i = 0, j = 0; i < ser_file->frame_count; i++, j++) {
		while (!successful_frames[j] && j < nb_frames) j++;
		if (i != j) {
			// move image j to i
			if (!buffer) {
				buffer = malloc(frame_size);
				if (!buffer) {
					PRINT_ALLOC_ERR;
					return 1;
				}
				siril_log_message(_("Compacting SER file after parallel output to it...\n"));
			}
			offseti = SER_HEADER_LEN + frame_size * (int64_t)i;
			offsetj = SER_HEADER_LEN + frame_size * (int64_t)j;

			if ((int64_t)-1 == fseek64(ser_file->file, offsetj, SEEK_SET)) {
				perror("seek");
				free(buffer);
				return 1;
			}
			if (fread(buffer, 1, frame_size, ser_file->file) != frame_size) {
				perror("fread");
				free(buffer);
				return 1;
			}
			if ((int64_t)-1 == fseek64(ser_file->file, offseti, SEEK_SET)) {
				perror("seek");
				free(buffer);
				return 1;
			}
			if (fwrite(buffer, 1, frame_size, ser_file->file) != frame_size) {
				perror("fwrite");
				free(buffer);
				return 1;
			}

			ser_file->ts[i] = ser_file->ts[j];
		}
	}

	free(buffer);
	return 0;
}

/* ser_file must be allocated
 * the file is created with no image size, the first image added will set it. */
int ser_create_file(const char *filename, struct ser_struct *ser_file,
		gboolean overwrite, struct ser_struct *copy_from) {
	if (overwrite)
		g_unlink(filename);
	if ((ser_file->file = g_fopen(filename, "w+b")) == NULL) {
		perror("open SER file for creation");
		return 1;
	}

	ser_file->filename = strdup(filename);
	ser_file->ts = NULL;
	ser_file->ts_alloc = 0;
	ser_file->fps = -1.0;
	ser_file->frame_count = 0;	// incremented on image add

	if (copy_from) {
		memcpy(&ser_file->lu_id, &copy_from->lu_id, 12);
		memset(&ser_file->image_width, 0, 16);
		memcpy(&ser_file->date, &copy_from->date, 8);
		memcpy(&ser_file->date_utc, &copy_from->date_utc, 8);
		ser_file->file_id = strdup(copy_from->file_id);
		memcpy(ser_file->observer, copy_from->observer, 40);
		memcpy(ser_file->instrument, copy_from->instrument, 40);
		memcpy(ser_file->telescope, copy_from->telescope, 40);
		ser_file->byte_pixel_depth = copy_from->byte_pixel_depth;
		ser_file->number_of_planes = 0;	// used as an indicator of new SER

		if (copy_from->ts && copy_from->frame_count > 0) {
			ser_file->ts = calloc(8, copy_from->frame_count);
			ser_file->ts_alloc = copy_from->frame_count;
		}
		/* we write the header now, but it should be written again
		 * before closing in case the number of the image in the new
		 * SER changes from the copied SER */
		ser_write_header(ser_file);
	} else {	// new SER
		ser_file->file_id = strdup("LUCAM-RECORDER");
		ser_file->lu_id = 0;
		ser_file->color_id = SER_MONO;	// this is 0
		ser_file->little_endian = SER_LITTLE_ENDIAN; // what will it do on big endian machine?
		memset(ser_file->observer, 0, 40);
		memset(ser_file->instrument, 0, 40);
		memset(ser_file->telescope, 0, 40);
		memset(&ser_file->date, 0, 8);
		memset(&ser_file->date_utc, 0, 8);
		ser_file->number_of_planes = 0;	// used as an indicator of new SER
	}
#ifdef _OPENMP
	omp_init_lock(&ser_file->fd_lock);
	omp_init_lock(&ser_file->ts_lock);
#endif
	siril_log_message(_("Created SER file %s\n"), filename);
	return 0;
}

int ser_open_file(const char *filename, struct ser_struct *ser_file) {
	if (ser_file->file) {
		fprintf(stderr, "SER: file already opened, or badly closed\n");
		return -1;
	}
	ser_file->file = g_fopen(filename, "r+b"); // now we can fix broken file, so not O_RDONLY anymore
	if (ser_file->file == NULL) {
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
	omp_init_lock(&ser_file->ts_lock);
#endif
	return 0;
}

int ser_close_file(struct ser_struct *ser_file) {
	int retval = 0;
	if (!ser_file)
		return -1;
	if (ser_file->file) {
		retval = fclose(ser_file->file);
		ser_file->file = NULL;
	}
	if (ser_file->file_id)
		free(ser_file->file_id);
	if (ser_file->ts)
		free(ser_file->ts);
	if (ser_file->filename)
		free(ser_file->filename);
#ifdef _OPENMP
	omp_destroy_lock(&ser_file->fd_lock);
	omp_destroy_lock(&ser_file->ts_lock);
#endif
	ser_init_struct(ser_file);
	return retval;
}

void ser_init_struct(struct ser_struct *ser_file) {
	memset(ser_file, 0, sizeof(struct ser_struct));
}

/* reads a frame on an already opened SER sequence.
 * frame number starts at 0 */
int ser_read_frame(struct ser_struct *ser_file, int frame_no, fits *fit) {
	int retval = 0, i, j, swap = 0;
	int64_t offset, frame_size;
	size_t read_size;
	WORD *olddata, *tmp;
	if (!ser_file || ser_file->file == NULL || !ser_file->number_of_planes ||
			!fit || frame_no < 0 || frame_no >= ser_file->frame_count)
		return -1;

	frame_size = ser_file->image_width * ser_file->image_height *
			ser_file->number_of_planes;
	read_size = frame_size * ser_file->byte_pixel_depth;

	olddata = fit->data;
	if ((fit->data = realloc(fit->data, frame_size * sizeof(WORD))) == NULL) {
		PRINT_ALLOC_ERR;
		if (olddata)
			free(olddata);
		return -1;
	}

	offset = SER_HEADER_LEN	+ frame_size *
		(int64_t)ser_file->byte_pixel_depth * (int64_t)frame_no;
	/*fprintf(stdout, "offset is %lu (frame %d, %d pixels, %d-byte)\n", offset,
	 frame_no, frame_size, ser_file->pixel_bytedepth);*/
#ifdef _OPENMP
	omp_set_lock(&ser_file->fd_lock);
#endif
	if ((int64_t)-1 == fseek64(ser_file->file, offset, SEEK_SET)) {
		perror("fseek in SER");
		retval = -1;
	} else {
		if (fread(fit->data, 1, read_size, ser_file->file) != read_size)
			retval = -1;
	}
#ifdef _OPENMP
	omp_unset_lock(&ser_file->fd_lock);
#endif
	if (retval)
		return -1;

	ser_manage_endianess_and_depth(ser_file, fit->data, frame_size);

	fit->bitpix = (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8) ? BYTE_IMG : USHORT_IMG;
	fit->orig_bitpix = fit->bitpix;

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
			bayer = get_SER_Bayer_Pattern(type_ser);
			if (bayer != com.debayer.bayer_pattern) {
				if (bayer == BAYER_FILTER_NONE  && user_warned == FALSE) {
					siril_log_color_message(_("No Bayer pattern found in the header file.\n"), "red");
				}
				else {
					if (user_warned == FALSE) {
						siril_log_color_message(_("Bayer pattern found in header (%s) is different"
								" from Bayer pattern in settings (%s). Overriding settings.\n"),
								"red", filter_pattern[bayer], filter_pattern[com.debayer.bayer_pattern]);
					}
					com.debayer.bayer_pattern = bayer;
				}
				user_warned = TRUE;
			}
		}
		/* for performance consideration (and many others) we force the interpolation algorithm
		 * to be BAYER_BILINEAR
		 */
		debayer(fit, BAYER_BILINEAR, FALSE);
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

	/* copy the SER timestamp to the fits */
	if (ser_file->ts) {
		char *timestamp = ser_timestamp(ser_file->ts[frame_no]);
		if (timestamp) {
			g_snprintf(fit->date_obs, FLEN_VALUE, "%s", timestamp);
			free(timestamp);
		}
	}

	fits_flip_top_to_bottom(fit);
	fit->top_down = FALSE;
	return 0;
}

/* multi-type cropping, works in constant space if needed */
#define crop_area_from_lines(BUFFER_TYPE) { \
	int x, y, src, dst = 0; \
	BUFFER_TYPE *inbuf = (BUFFER_TYPE *)read_buffer; \
	BUFFER_TYPE *out = (BUFFER_TYPE *)outbuf; \
	for (y = 0; y < area->h; y++) { \
		src = y * ser_file->image_width + area->x; \
		for (x = 0; x < area->w; x++) \
			out[dst++] = inbuf[src++]; \
	} \
}

/* multi-type RGB reordering, works in constant space if needed */
#define crop_area_from_color_lines(BUFFER_TYPE) { \
	int x, y, src, dst = 0; \
	BUFFER_TYPE *inbuf = (BUFFER_TYPE *)read_buffer; \
	BUFFER_TYPE *out = (BUFFER_TYPE *)outbuf; \
	int color_offset; \
	if (ser_file->color_id == SER_BGR) { \
		color_offset = 2 - layer; \
	} else { \
		color_offset = layer; \
	} \
	for (y = 0; y < area->h; y++) { \
		src = (y * ser_file->image_width + area->x) * 3 + color_offset; \
		for (x = 0; x < area->w; x++) { \
			out[dst++] = inbuf[src]; \
			src += 3; \
		} \
	} \
}

/* reading an area from a SER frame, for one layer only, either layer == -1 for
 * monochrome and debayer, or 0-2 for color.
 * the area is read in one read(2) call to limit the number of syscalls, for a
 * full-width area of same height as requested, then cropped horizontally to
 * get the requested area.
 * This function is the first one of siril to handle two different data types
 * (BYTE and WORD) for the same algorithm! This uses VIPS-style macros.
 * */
static int read_area_from_image(struct ser_struct *ser_file, const int frame_no,
		WORD *outbuf, const rectangle *area, const int layer) {
	int64_t offset, frame_size;
	int retval = 0;
	WORD *read_buffer;
	size_t read_size = ser_file->image_width * area->h * ser_file->byte_pixel_depth;
	if (layer != -1) read_size *= 3;
	if (layer != -1 || area->w != ser_file->image_width) {
		// allocated space is probably not enough to
		// store whole lines or RGB data
		read_buffer = malloc(read_size);
	}
	else read_buffer = outbuf;

	frame_size = ser_file->image_width * ser_file->image_height *
		ser_file->number_of_planes * ser_file->byte_pixel_depth;

#ifdef _OPENMP
	omp_set_lock(&ser_file->fd_lock);
#endif
	// we read the full-stride rectangle that contains the requested area
	offset = SER_HEADER_LEN + frame_size * frame_no +	// requested frame
		area->y * ser_file->image_width *
		ser_file->byte_pixel_depth * (layer != -1 ? 3 : 1);	// requested area

	if ((int64_t)-1 == fseek64(ser_file->file, offset, SEEK_SET)) {
		perror("fseek in SER");
		retval = -1;
	} else {
		if (fread(read_buffer, 1, read_size, ser_file->file) != read_size) {
			retval = -1;
		}
	}
#ifdef _OPENMP
	omp_unset_lock(&ser_file->fd_lock);
#endif
	if (!retval) {
		if (area->w != ser_file->image_width) {
			// here we crop x-wise our area
			if (layer != -1) {
				/* reorder the RGBRGB to RRGGBB and crop */
				if (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8) {
					crop_area_from_color_lines(BYTE);
				} else if (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_16) {
					crop_area_from_color_lines(WORD);
				}
			} else {
				/* just crop */
				if (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8) {
					crop_area_from_lines(BYTE);
				} else if (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_16) {
					crop_area_from_lines(WORD);
				}
			}
		} else if (layer != -1) {
			/* just reorder RGB data, the crop function works too */
			if (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8) {
				crop_area_from_color_lines(BYTE);
			} else if (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_16) {
				crop_area_from_color_lines(WORD);
			}
		}
	}
	if (layer != -1 || area->w != ser_file->image_width)
		free(read_buffer);
	return retval;
}

/* read an area of an image in an opened SER sequence */
int ser_read_opened_partial(struct ser_struct *ser_file, int layer,
		int frame_no, WORD *buffer, const rectangle *area) {
	int xoffset, yoffset, x, y;
	ser_color type_ser;
	WORD *rawbuf, *demosaiced_buf;
	rectangle debayer_area, image_area;
	sensor_pattern sensortmp;

	if (!ser_file || ser_file->file == NULL || frame_no < 0
			|| frame_no >= ser_file->frame_count)
		return -1;

	type_ser = ser_file->color_id;
	if (!com.debayer.open_debayer &&
			type_ser != SER_RGB && type_ser != SER_BGR)
		type_ser = SER_MONO;

	switch (type_ser) {
	case SER_MONO:
		if (read_area_from_image(ser_file, frame_no, buffer, area, -1))
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
			bayer = get_SER_Bayer_Pattern(type_ser);
			if (bayer != com.debayer.bayer_pattern) {
				if (bayer == BAYER_FILTER_NONE && user_warned == FALSE) {
					siril_log_color_message(_("No Bayer pattern found in the header file.\n"), "red");
				}
				else {
					if (user_warned == FALSE) {
						siril_log_color_message(_("Bayer pattern found in header (%s) is different"
								" from Bayer pattern in settings (%s). Overriding settings.\n"),
								"red", filter_pattern[bayer], filter_pattern[com.debayer.bayer_pattern]);
					}
					com.debayer.bayer_pattern = bayer;
				}
				user_warned = TRUE;
			}
		}
		if (layer < 0 || layer >= 3) {
			siril_log_message(_("For a demosaiced image, layer has to be R, G or B (0 to 2).\n"));
			return -1;
		}

		image_area = (rectangle) { .x = 0, .y = 0,
			.w = ser_file->image_width, .h = ser_file->image_height };
		get_debayer_area(area, &debayer_area, &image_area, &xoffset, &yoffset);

		// allocating a buffer for WORD because it's going to be converted in-place
		rawbuf = malloc(debayer_area.w * debayer_area.h * sizeof(WORD));
		if (!rawbuf) {
			PRINT_ALLOC_ERR;
			return -1;
		}
		if (read_area_from_image(ser_file, frame_no, rawbuf, &debayer_area, -1)) {
			free(rawbuf);
			return -1;
		}
		ser_manage_endianess_and_depth(ser_file, rawbuf, debayer_area.w * debayer_area.h);

		/* for performance consideration (and many others) we force the interpolation algorithm
		 * to be BAYER_BILINEAR
		 */
		demosaiced_buf = debayer_buffer(rawbuf, &debayer_area.w,
				&debayer_area.h, BAYER_BILINEAR, com.debayer.bayer_pattern,
				NULL);
		free(rawbuf);
		if (!demosaiced_buf)
			return -1;

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
		if (read_area_from_image(ser_file, frame_no, buffer, area, layer))
			return -1;
		ser_manage_endianess_and_depth(ser_file, buffer, area->w * area->h);
		break;
	default:
		siril_log_message(_("This type of Bayer pattern is not handled yet.\n"));
		return -1;
	}

	return 0;
}

int ser_read_opened_partial_fits(struct ser_struct *ser_file, int layer,
		int frame_no, fits *fit, const rectangle *area) {
	if (new_fit_image(&fit, area->w, area->h, 1))
		return -1;
	fit->top_down = TRUE;
	if (ser_file->ts) {
		char *timestamp = ser_timestamp(ser_file->ts[frame_no]);
		if (timestamp) {
			g_snprintf(fit->date_obs, FLEN_VALUE, "%s", timestamp);
			free(timestamp);
		}
	}
	return ser_read_opened_partial(ser_file, layer, frame_no, fit->pdata[0], area);
}

int ser_write_frame_from_fit(struct ser_struct *ser_file, fits *fit, int frame_no) {
	int pixel, plane, dest;
	int ret, retval = 0;
	int64_t offset, frame_size;
	BYTE *data8 = NULL;			// for 8-bit files
	WORD *data16 = NULL;		// for 16-bit files

	if (!ser_file || ser_file->file == NULL || !fit)
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

	offset = SER_HEADER_LEN	+ frame_size *
			(int64_t)ser_file->byte_pixel_depth * (int64_t)frame_no;

	if (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8) {
		data8 = malloc(frame_size * ser_file->byte_pixel_depth);
		if (!data8) return -1;
	} else {
		data16 = malloc(frame_size * ser_file->byte_pixel_depth);
		if (!data16) return -1;
	}

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
	if ((int64_t)-1 == fseek64(ser_file->file, offset, SEEK_SET)) {
#ifdef _OPENMP
		omp_unset_lock(&ser_file->fd_lock);
#endif
		perror("seek");
		retval = -1;
		goto free_and_quit;
	}

	if (ser_file->byte_pixel_depth == SER_PIXEL_DEPTH_8) {
		ret = fwrite(data8, 1, frame_size * ser_file->byte_pixel_depth, ser_file->file);
#ifdef _OPENMP
		omp_unset_lock(&ser_file->fd_lock);
#endif
		if (ret != frame_size * ser_file->byte_pixel_depth) {
			perror("write image in SER");
			retval = 1;
			goto free_and_quit;
		}
	} else {
		ret = fwrite(data16, 1, frame_size * ser_file->byte_pixel_depth, ser_file->file);
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

	if (!ser_alloc_ts(ser_file, frame_no)) {
		uint64_t utc, local;
		FITS_date_key_to_Unix_time(fit->date_obs, &utc, &local);
		ser_file->ts[frame_no] = utc;
	}

free_and_quit:
	if (data8) free(data8);
	if (data16) free(data16);
	return retval;
}

int64_t ser_compute_file_size(struct ser_struct *ser_file, int nb_frames) {
	int64_t frame_size, size = ser_file->filesize;

	if (nb_frames != ser_file->frame_count) {
		frame_size = (size - SER_HEADER_LEN) / ser_file->frame_count;
		size = SER_HEADER_LEN + frame_size * nb_frames;
	}
	/* SER can be demosaiced on the fly on creation.
	 * TODO: Is this the good test? */
	if (ser_is_cfa(ser_file) && com.debayer.open_debayer)
		size *= 3;
	return size;
}
