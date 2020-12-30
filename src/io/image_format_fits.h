#ifndef _IMAGE_FORMAT_FITS_H
#define _IMAGE_FORMAT_FITS_H

#include <fitsio.h>
#include "core/siril.h"

/****************** image_format_fits.h ******************/
void read_fits_header(fits *fit);
char *copy_header(fits *fit);
data_type get_data_type(int bitpix);
void fit_get_photometry_data(fits *fit);
int readfits(const char *filename, fits *fit, char *realname, gboolean force_float);
void get_date_data_from_fitsfile(fitsfile *fptr, GDateTime **dt, double *exposure);
int import_metadata_from_fitsfile(fitsfile *fptr, fits *to);
void clearfits(fits*);
int readfits_partial(const char *filename, int layer, fits *fit,
		const rectangle *area, gboolean read_date);
void flip_buffer(int bitpix, void *buffer, const rectangle *area);
int read_opened_fits_partial(sequence *seq, int layer, int index, void *buffer,
		const rectangle *area);
int siril_fits_compress(fits *f);
int save_opened_fits(fits *f);
int savefits(const char*, fits*);
int copyfits(fits *from, fits *to, unsigned char oper, int layer);
int copy_fits_metadata(fits *from, fits *to);
int copy_fits_from_file(char *source, char *destination);
int save1fits16(const char *filename, fits *fit, int layer);
int save1fits32(const char *filename, fits *fit, int layer);
int siril_fits_open_diskfile(fitsfile **fptr, const char *filename, int iomode,
		int *status);

void rgb24bit_to_fits48bit(unsigned char *rgbbuf, fits *fit, gboolean inverted);
void rgb8bit_to_fits16bit(unsigned char *graybuf, fits *fit);
void rgb48bit_to_fits48bit(WORD *rgbbuf, fits *fit, gboolean inverted,
		gboolean change_endian);

void fits_flip_top_to_bottom(fits *fit);
void extract_region_from_fits(fits *from, int layer, fits *to,
		const rectangle *area);
int new_fit_image(fits **fit, int width, int height, int nblayer, data_type type);
int new_fit_image_with_data(fits **fit, int width, int height, int nblayer, data_type type, void *data);
void fit_replace_buffer(fits *fit, void *newbuf, data_type newtype);
void fit_debayer_buffer(fits *fit, void *newbuf);

void keep_first_channel_from_fits(fits *fit);
GdkPixbuf* get_thumbnail_from_fits(char *filename, gchar **descr);

// internal read of FITS file, for FITS images and FITS sequences
void manage_bitpix(fitsfile *fptr, int *bitpix, int *orig_bitpix);
int read_fits_with_convert(fits* fit, const char* filename, gboolean force_float);
int internal_read_partial_fits(fitsfile *fptr, unsigned int ry,
		int bitpix, void *dest, int layer, const rectangle *area);
int siril_fits_create_diskfile(fitsfile **fptr, const char *filename, int *status);
void save_fits_header(fits *fit);
void report_fits_error(int status);

#endif
