/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2021 team free-astro (see more in AUTHORS file)
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
 */

/**
 *
 * \file utils.c
 * \brief Misc. function utilities.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#ifdef _WIN32
#include <windows.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_app_dirs.h"
#include "core/exif.h"
#include "io/conversion.h"
#include "io/ser.h"
#include "io/sequence.h"
#include "gui/utils.h"
#include "gui/progress_and_log.h"
#include "io/single_image.h"

/**
 * Round double value to an integer
 * @param x value to round
 * @return an integer
 */
int round_to_int(double x) {
	if (x <= INT_MIN + 0.5) return INT_MIN;
	if (x >= INT_MAX - 0.5) return INT_MAX;
	if (x >= 0.0)
		return (int)(x + 0.5);
	return (int)(x - 0.5);
}

/**
 * Round float value to an integer
 * @param x value to round
 * @return an integer
 */
int roundf_to_int(float x) {
	if (x <= INT_MIN + 0.5f) return INT_MIN;
	if (x >= INT_MAX - 0.5f) return INT_MAX;
	if (x >= 0.0f)
		return (int)(x + 0.5f);
	return (int)(x - 0.5f);
}

/**
 * Round double value to a WORD
 * @param x value to round
 * @return a WORD
 */
WORD round_to_WORD(double x) {
	if (x <= 0.0)
		return (WORD)0;
	if (x > USHRT_MAX_DOUBLE)
		return USHRT_MAX;
	return (WORD)(x + 0.5);
}

/**
 * Round double value to a BYTE
 * @param x value to round
 * @return a BYTE
 */
BYTE round_to_BYTE(double x) {
	if (x <= 0.0)
		return (BYTE)0;
	if (x > UCHAR_MAX_DOUBLE)
		return UCHAR_MAX;
	return (BYTE)(x + 0.5);
}

/**
 * Round float value to a BYTE
 * @param f value to round
 * @return a truncated and rounded BYTE
 */
BYTE roundf_to_BYTE(float f) {
	if (f <= 0.5f) return 0;
	if (f >= UCHAR_MAX - 0.5f) return UCHAR_MAX;
	return (BYTE)(f + 0.5f);
}

/**
 * Round float value to a WORD
 * @param f value to round
 * @return a truncated and rounded WORD
 */
WORD roundf_to_WORD(float f) {
	if (f <= 0.5f) return 0;
	if (f >= USHRT_MAX - 0.5f) return USHRT_MAX;
	return (WORD)(f + 0.5f);
}

/**
 * Compute a ceiling factor
 * @param x the number to test
 * @param factor the factor
 * @return x if it is a factor of factor or the next factor
 */
int round_to_ceiling_multiple(int x, int factor) {
	if (x % factor == 0)
		return x;
	return (x / factor + 1) * factor;
}

/**
 * convert double value to a BYTE
 * @param x value to convert
 * @return a BYTE
 */
BYTE conv_to_BYTE(double x) {
	if (x == 0.0)
		return (BYTE)0;
	if (x == USHRT_MAX_DOUBLE)
		return UCHAR_MAX;
	x = ((x / USHRT_MAX_DOUBLE) * UCHAR_MAX_DOUBLE);
	return((BYTE)(x));
}

/**
 * truncate a 64 bit unsigned int to a 32 bit signed int
 * @param x value to truncate
 * @return an int
 */
int truncate_to_int32(uint64_t x) {
	if (x > (uint64_t)INT_MAX)
		return INT_MAX;
	return (int)x;
}

WORD truncate_to_WORD(int x) {
	if (x < 0)
		return 0;
	if (x > USHRT_MAX)
		return USHRT_MAX;
	return (WORD)x;
}

/**
 * Clamp an integer value in the interval given by [low, high]
 * @param val value to be checked
 * @param low low value of the interval
 * @param high high value of the interval
 * @return a new value set in the [low, high] interval
 */
int set_int_in_interval(int val, int low, int high) {
	return max(low, min(val, high));
}

/**
 * Clamp a float value in the interval given by [low, high]
 * @param val value to be checked
 * @param low low value of the interval
 * @param high high value of the interval
 * @return a new value set in the [low, high] interval
 */
float set_float_in_interval(float val, float low, float high) {
	return max(low, min(val, high));
}

/**
 * Clamp a double value in the interval given by [low, high]
 * @param val value to be checked
 * @param low low value of the interval
 * @param high high value of the interval
 * @return a new value set in the [low, high] interval
 */
double set_double_in_interval(double val, double low, double high) {
	return max(low, min(val, high));
}

/**
 * convert an unsigned short value to siril's representation of float values [0, 1]
 * @param w value to convert
 * @return the float equivalent
 */
float ushort_to_float_range(WORD w) {
	return (float)w * INV_USHRT_MAX_SINGLE;
}

/**
 * convert an unsigned char value to siril's representation of float values [0, 1]
 * @param w value to convert
 * @return the float equivalent
 */
float uchar_to_float_range(BYTE w) {
	return (float)w * INV_UCHAR_MAX_SINGLE;
}

/**
 * convert an double value from the unsigned short range to siril's representation
 * of float values [0, 1]
 * @param d value to convert
 * @return the float equivalent
 */
float double_ushort_to_float_range(double d) {
	return (float)d * INV_USHRT_MAX_SINGLE;
}

/**
 * convert a siril float [0, 1] to an unsigned short
 * @param f value to convert
 * @return the unsigned short equivalent
 */
WORD float_to_ushort_range(float f) {
	return roundf_to_WORD(f * USHRT_MAX_SINGLE);
}

/**
 * convert a siril float [0, 1] to an unsigned char
 * @param f value to convert
 * @return the unsigned char equivalent
 */
BYTE float_to_uchar_range(float f) {
	return roundf_to_BYTE(f * UCHAR_MAX_SINGLE);
}

/**
 * convert the pixel value of an image to a float [0, 1] normalized using bitpix
 * value depending on btpix
 * @param fit the image the data is from
 * @return a float [0, 1] value for the given integer value
 */
float ushort_to_float_bitpix(fits *fit, WORD value) {
	float fval = (float)value;
	return fit->orig_bitpix == BYTE_IMG ?
		fval * INV_UCHAR_MAX_SINGLE :
		fval * INV_USHRT_MAX_SINGLE;
}

/**
 * convert a float type buffer into a WORD buffer
 * @param buffer in float
 * @param ndata
 * @return
 */
WORD *float_buffer_to_ushort(float *buffer, size_t ndata) {
	WORD *buf = malloc(ndata * sizeof(WORD));
	if (!buf) {
		PRINT_ALLOC_ERR;
	} else {
		for (size_t i = 0; i < ndata; i++) {
			buf[i] = float_to_ushort_range(buffer[i]);
		}
	}
	return buf;
}

/**
 * convert a BYTE type buffer into a float buffer
 * @param buffer in BYTE
 * @param ndata
 * @return
 */
float *uchar_buffer_to_float(BYTE *buffer, size_t ndata) {
	float *buf = malloc(ndata * sizeof(float));
	if (!buf) {
		PRINT_ALLOC_ERR;
	} else {
		for (size_t i = 0; i < ndata; i++) {
			buf[i] = uchar_to_float_range(buffer[i]);
		}
	}
	return buf;
}

/**
 * convert a WORD type buffer into a float buffer
 * @param buffer in WORD
 * @param ndata
 * @return
 */
float *ushort_buffer_to_float(WORD *buffer, size_t ndata) {
	float *buf = malloc(ndata * sizeof(float));
	if (!buf) {
		PRINT_ALLOC_ERR;
	} else {
		for (size_t i = 0; i < ndata; i++) {
			buf[i] = ushort_to_float_range(buffer[i]);
		}
	}
	return buf;
}

/**
 * convert a WORD type buffer representing 8bit data into a float buffer
 * @param buffer in WORD
 * @param ndata
 * @return
 */
float *ushort8_buffer_to_float(WORD *buffer, size_t ndata) {
	float *buf = malloc(ndata * sizeof(float));
	if (!buf) {
		PRINT_ALLOC_ERR;
	} else {
		for (size_t i = 0; i < ndata; i++) {
			buf[i] = uchar_to_float_range((BYTE) buffer[i]);
		}
	}
	return buf;
}

/**
 * change endianness of a 16 bit unsigned int
 * @param x value to convert
 * @return byte-swapped value
 */
uint16_t change_endianness16(uint16_t x) {
    return (x >> 8) | (x << 8);
}

/**
 * convert a 16 bit unsigned int in CPU byte order to little endian
 * @param x value to convert
 * @return little endian value
 */
uint16_t cpu_to_le16(uint16_t x) {
#ifdef __BIG_ENDIAN__
    return change_endianness16(x);
#else
    return x;
#endif
}

/**
 * convert a 16 bit unsigned int in CPU byte order to big endian
 * @param x value to convert
 * @return big endian value
 */
uint16_t cpu_to_be16(uint16_t x) {
#ifdef __BIG_ENDIAN__
    return x;
#else
    return change_endianness16(x);
#endif
}

/**
 * convert a 16 bit unsigned int from little endian to CPU byte order
 * @param x little endian value to convert
 * @return value
 */
uint16_t le16_to_cpu(uint16_t x) {
    return cpu_to_le16(x);
}

/**
 * convert a 16 bit unsigned int from big endian to CPU byte order
 * @param x big endian value to convert
 * @return value
 */
uint16_t be16_to_cpu(uint16_t x) {
    return cpu_to_be16(x);
}

/**
 * change endianness of a 32 bit unsigned int
 * @param x value to convert
 * @return byte-swapped value
 */
uint32_t change_endianness32(uint32_t x) {
    return (x >> 24) | ((x & 0xFF0000) >> 8) | ((x & 0xFF00) << 8) | (x << 24);
}

/**
 * convert a 32 bit unsigned int in CPU byte order to little endian
 * @param x value to convert
 * @return little endian value
 */
uint32_t cpu_to_le32(uint32_t x) {
#ifdef __BIG_ENDIAN__
    return change_endianness32(x);
#else
    return x;
#endif
}

/**
 * convert a 32 bit unsigned int in CPU byte order to big endian
 * @param x value to convert
 * @return big endian value
 */
uint32_t cpu_to_be32(uint32_t x) {
#ifdef __BIG_ENDIAN__
    return x;
#else
    return change_endianness32(x);
#endif
}

/**
 * convert a 32 bit unsigned int from little endian to CPU byte order
 * @param x little endian value to convert
 * @return value
 */
uint32_t le32_to_cpu(uint32_t x) {
    return cpu_to_le32(x);
}

/**
 * convert a 32 bit unsigned int from big endian to CPU byte order
 * @param x big endian value to convert
 * @return value
 */
uint32_t be32_to_cpu(uint32_t x) {
    return cpu_to_be32(x);
}

/**
 * change endianness of a 64 bit unsigned int
 * @param x value to convert
 * @return byte-swapped value
 */
uint64_t change_endianness64(uint64_t x) {
    return
        (x >> 56)
        | ((x & 0xFF000000000000) >> 40)
        | ((x & 0xFF0000000000) >> 24)
        | ((x & 0xFF00000000) >> 8)
        | ((x & 0xFF000000) << 8)
        | ((x & 0xFF0000) << 24)
        | ((x & 0xFF00) << 40)
        | (x << 56);
}

/**
 * convert a 64 bit unsigned int in CPU byte order to little endian
 * @param x value to convert
 * @return little endian value
 */
uint64_t cpu_to_le64(uint64_t x) {
#ifdef __BIG_ENDIAN__
    return change_endianness64(x);
#else
    return x;
#endif
}

/**
 * convert a 64 bit unsigned int in CPU byte order to big endian
 * @param x value to convert
 * @return big endian value
 */
uint64_t cpu_to_be64(uint64_t x) {
#ifdef __BIG_ENDIAN__
    return x;
#else
    return change_endianness64(x);
#endif
}

/**
 * convert a 64 bit unsigned int from little endian to CPU byte order
 * @param x little endian value to convert
 * @return value
 */
uint64_t le64_to_cpu(uint64_t x) {
    return cpu_to_le64(x);
}

/**
 * convert a 64 bit unsigned int from big endian to CPU byte order
 * @param x big endian value to convert
 * @return value
 */
uint64_t be64_to_cpu(uint64_t x) {
    return cpu_to_be64(x);
}

/**
 * Test if fit has 3 channels
 * @param fit input FITS image
 * @return TRUE if fit image has 3 channels
 */
gboolean isrgb(fits *fit) {
	return (fit->naxis == 3);
}

/**
 *  Searches for an extension '.something' in filename from the end
 *  @param filename input filename or path
 *  @return the index of the first '.' found
 */
int get_extension_index(const char *filename) {
	if (filename == NULL || filename[0] == '\0')
		return -1;
	for (int i = strlen(filename) - 1; i > 0; i--) {
		if (filename[i] == '\\' || filename[i] == '/')
			break;
		if (filename[i] == '.')
			return i;
	}
	return -1;
}

/**
 * Get the extension of a file, without the dot.
 * @param filename input filename or path
 * @return extension pointed from the filename itself or NULL
 */
const char *get_filename_ext(const char *filename) {
	gchar *basename;
	int len;
	const char *dot, *p;

	basename = g_path_get_basename(filename);
	len = strlen(filename) - strlen(basename);
	g_free(basename);

	p = filename + len;
	dot = strrchr(p, '.');
	if (!dot || dot == p) {
		return NULL;
	}
	return dot + 1;
}

/**
 *
 * @param filename
 * @return the type of the file from its filename
 */
image_type get_type_from_filename(const gchar *filename) {
	const char *ext = get_filename_ext(filename);
	if (!ext)
		return TYPEUNDEF;
	return get_type_for_extension(ext);
}

/**
 * Removes extension of the filename or path
 * @param filename file path with extension
 * @return newly allocated filename without extension
 */
char *remove_ext_from_filename(const char *filename) {
	char *file = NULL;
	int ext_index = -1;

	for (int i = strlen(filename) - 1; i > 0; i--) {
		if (filename[i] == '\\' || filename[i] == '/')
			break;
		if (filename[i] == '.') {
			ext_index = i;
			break;
		}
	}
	if (ext_index == -1)
		return strdup(filename);

	file = malloc(ext_index + 1);
	strncpy(file, filename, ext_index);
	file[ext_index] = '\0';
	return file;
}

/**
 * Replaces the extension of a file name or path
 * @param path the original path
 * @param new_ext the next extension to put
 * @return a new string with the new extension
 */
gchar *replace_ext(const char *path, const char *new_ext) {
	int idx = get_extension_index(path);
	gchar *retval = g_strdup(path);
	if (idx != -1)
		retval[idx] = '\0';
	return str_append(&retval, new_ext);
}

/**
 * Check is a string contains directory separators and thus represent a path
 * @param file the string to test
 * @return true if it contains a separator
 */
gboolean string_is_a_path(const char *file) {
	int len = strlen(file);
	for (int i = 0; i < len; i++) {
		if (file[i] == '\\' || file[i] == '/')
			return TRUE;
	}
	return FALSE;
}

/**
 * Tests whether the given file is either regular or a symlink
 * @param filename input
 * @return 1 if file is readable (not actually opened to verify)
 */
int is_readable_file(const char *filename) {
	GStatBuf sts;
	if (g_lstat(filename, &sts))
		return 0;
	if (S_ISREG (sts.st_mode)
#ifndef _WIN32
			|| S_ISLNK(sts.st_mode)
#else
		|| (GetFileAttributesA(filename) & FILE_ATTRIBUTE_REPARSE_POINT )
#endif
	)
		return 1;
	return 0;
}

static gchar forbidden_char[] = { '/', '\\', '"', '\'' };

gboolean is_forbiden_in_filename(gchar c) {
	for (int i = 0; i < G_N_ELEMENTS(forbidden_char); i++) {
		if (c == forbidden_char[i])
			return TRUE;
	}
	return FALSE;
}

gboolean file_name_has_invalid_chars(const char *name) {
	if (!name)
		return TRUE;	// NULL is kind of invalid
	for (int i = 0; i < strlen(name); i++)
		if (is_forbiden_in_filename(name[i]))
			return TRUE;
	return FALSE;
}

/** Tests if filename is the canonical name of a known file type
 *  If filename contains an extension, only this file name is tested, else all
 *  extensions are tested for the file name until one is found.
 * @param[in] filename the filename to test for.
 * @param[in] type is set according to the result of the test.
 * @param[out] realname (optionnal) is set according to the found file name: it
 *  must be freed with when no longer needed.
 * @return 0 if sucess, 1 if error
 */
int stat_file(const char *filename, image_type *type, char **realname) {
	int k;
	const char *ext;
	*type = TYPEUNDEF;	// default value

	/* check for an extension in filename and isolate it, including the . */
	if (filename[0] == '\0')
		return 1;

	ext = get_filename_ext(filename);
	/* if filename has an extension, we only test for it */
	if (ext) {
		if (is_readable_file(filename)) {
			if (realname)
				*realname = strdup(filename);
			*type = get_type_for_extension(ext);
			return 0;
		}
		return 1;
	}

	/* else, we can test various file extensions */
	/* first we test lowercase, then uppercase */
	for (k = 0; k < 2; k++) {
		int i = 0;
		while (supported_extensions[i]) {
			GString *testName = g_string_new(filename);
			if (k == 0) {
				testName = g_string_append(testName, supported_extensions[i]);
			} else {
				gchar *tmp = g_ascii_strup(supported_extensions[i],
						strlen(supported_extensions[i]));
				testName = g_string_append(testName, tmp);
				g_free(tmp);
			}
			gchar *name = g_string_free(testName, FALSE);

			if (is_readable_file(name)) {
				*type = get_type_for_extension(supported_extensions[i] + 1);
				assert(*type != TYPEUNDEF);
				if (realname)
					*realname = strdup(name);
				g_free(name);
				return 0;
			}
			i++;
			g_free(name);
		}
	}
	return 1;
}

static gchar* siril_canonicalize_filename(const gchar *filename,
		const gchar *relative_to)
{
#if GLIB_CHECK_VERSION(2,58,0)
	return g_canonicalize_filename(filename, relative_to);
}
#else
/**
 * g_canonicalize_filename:
 * @filename: (type filename): the name of the file
 * @relative_to: (type filename) (nullable): the relative directory, or %NULL
 * to use the current working directory
 *
 * Gets the canonical file name from @filename. All triple slashes are turned into
 * single slashes, and all `..` and `.`s resolved against @relative_to.
 *
 * Symlinks are not followed, and the returned path is guaranteed to be absolute.
 *
 * If @filename is an absolute path, @relative_to is ignored. Otherwise,
 * @relative_to will be prepended to @filename to make it absolute. @relative_to
 * must be an absolute path, or %NULL. If @relative_to is %NULL, it'll fallback
 * to g_get_current_dir().
 *
 * This function never fails, and will canonicalize file paths even if they don't
 * exist.
 *
 * No file system I/O is done.
 *
 * Returns: (type filename) (transfer full): a newly allocated string with the
 * canonical file path
 * Since: 2.58
 */
  gchar *canon, *start, *p, *q;
  guint i;

  g_return_val_if_fail (relative_to == NULL || g_path_is_absolute (relative_to), NULL);

  if (!g_path_is_absolute (filename))
    {
      gchar *cwd_allocated = NULL;
      const gchar  *cwd;

      if (relative_to != NULL)
        cwd = relative_to;
      else
        cwd = cwd_allocated = g_get_current_dir ();

      canon = g_build_filename (cwd, filename, NULL);
      g_free (cwd_allocated);
    }
  else
    {
      canon = g_strdup (filename);
    }

  start = (char *)g_path_skip_root (canon);

  if (start == NULL)
    {
      /* This shouldn't really happen, as g_get_current_dir() should
         return an absolute pathname, but bug 573843 shows this is
         not always happening */
      g_free (canon);
      return g_build_filename (G_DIR_SEPARATOR_S, filename, NULL);
    }

  /* POSIX allows double slashes at the start to
   * mean something special (as does windows too).
   * So, "//" != "/", but more than two slashes
   * is treated as "/".
   */
  i = 0;
  for (p = start - 1;
       (p >= canon) &&
         G_IS_DIR_SEPARATOR (*p);
       p--)
    i++;
  if (i > 2)
    {
      i -= 1;
      start -= i;
      memmove (start, start+i, strlen (start+i) + 1);
    }

  /* Make sure we're using the canonical dir separator */
  p++;
  while (p < start && G_IS_DIR_SEPARATOR (*p))
    *p++ = G_DIR_SEPARATOR;

  p = start;
  while (*p != 0)
    {
      if (p[0] == '.' && (p[1] == 0 || G_IS_DIR_SEPARATOR (p[1])))
        {
          memmove (p, p+1, strlen (p+1)+1);
        }
      else if (p[0] == '.' && p[1] == '.' && (p[2] == 0 || G_IS_DIR_SEPARATOR (p[2])))
        {
          q = p + 2;
          /* Skip previous separator */
          p = p - 2;
          if (p < start)
            p = start;
          while (p > start && !G_IS_DIR_SEPARATOR (*p))
            p--;
          if (G_IS_DIR_SEPARATOR (*p))
            *p++ = G_DIR_SEPARATOR;
          memmove (p, q, strlen (q)+1);
        }
      else
        {
          /* Skip until next separator */
          while (*p != 0 && !G_IS_DIR_SEPARATOR (*p))
            p++;

          if (*p != 0)
            {
              /* Canonicalize one separator */
              *p++ = G_DIR_SEPARATOR;
            }
        }

      /* Remove additional separators */
      q = p;
      while (*q && G_IS_DIR_SEPARATOR (*q))
        q++;

      if (p != q)
        memmove (p, q, strlen (q) + 1);
    }

  /* Remove trailing slashes */
  if (p > start && G_IS_DIR_SEPARATOR (*(p-1)))
    *(p-1) = 0;

  return canon;
}
#endif

/** Try to change the CWD to the argument, absolute or relative.
 *  If success, the new CWD is written to com.wd
 *  @param[in] dir absolute or relative path we want to set as cwd
 *  @param[out] err error message when return value is different of 1. Can be NULL if message is not needed.
 *  @return 0 if success, any other values for error
 */
int siril_change_dir(const char *dir, gchar **err) {
	gchar *error = NULL;
	int retval = 0;
	char *new_dir = NULL;

	if (dir == NULL || dir[0] == '\0') {
		error = siril_log_message(_("Unknown error\n"));
		retval = -1;
	} else if (!g_file_test(dir, G_FILE_TEST_EXISTS)) {
		error = siril_log_message(_("'%s' No such file or directory\n"), dir);
		retval = 2;
	} else if (!g_file_test(dir, G_FILE_TEST_IS_DIR)) {
		error = siril_log_message(_("'%s' is not a directory\n"), dir);
		retval = 3;
	} else if (g_access(dir, W_OK)) {
		error = siril_log_color_message(_("You don't have permission "
				"to write in this directory: '%s'\n"), "red", dir);
		retval = 4;
	} else {
		/* sequences are invalidate when cwd is changed */
		close_sequence(FALSE);
		if (!g_chdir(dir)) {
			/* do we need to search for sequences in the directory now? We still need to
			 * press the check seq button to display the list, and this is also done there. */
			/* check_seq();
			   update_sequence_list();*/
			// Don't follow symbolic links
			if (g_path_is_absolute(dir)) {
				new_dir = g_memdup(dir, strlen(dir) + 1);
				g_free(com.wd);
				com.wd = new_dir;
			} else {
				new_dir = siril_canonicalize_filename(dir, com.wd);
				g_free(com.wd);
				com.wd = new_dir;
			}

		  siril_log_message(_("Setting CWD (Current Working Directory) to '%s'\n"), com.wd);
		  retval = 0;
		} else {
			error = siril_log_message(_("Could not change directory to '%s'.\n"), dir);
			retval = 1;
		}
	}
	if (err) {
		*err = error;
	}
	return retval;
}
/**
 * If Windows OS, converts a filename from UTF-8 to the system codepage. Do nothing on other system
 * @param path to convert
 * @return converted filename
 */
gchar *get_locale_filename(const gchar *path) {
	gchar *str;
#ifdef _WIN32
	str = g_win32_locale_filename_from_utf8(path);
	if (!str) {
		siril_log_color_message("Conversion of the filename to system codepage failed. Please consider removing all wide chars.\n", "red");
		return g_strdup(path);
	}
#else // _WIN32
	str = g_strdup(path);
#endif // _WIN32
	return str;
}

#ifdef _WIN32
static int ListSequences(const gchar *sDir, const char *sequence_name_to_select,
		GtkComboBoxText *seqcombo, int *index_of_seq_to_load) {
	WIN32_FIND_DATAW fdFile;
	HANDLE hFind = NULL;
	char sPath[2048];
	char filename[256];
	int number_of_loaded_sequences = 0;
	wchar_t *wpath;

	//Specify a file mask. *.seq = We want only seq file!
	sprintf(sPath, "%s\\*.seq", sDir);

	wpath = g_utf8_to_utf16(sPath, -1, NULL, NULL, NULL);
	if (wpath == NULL)
		return 1;

	if ((hFind = FindFirstFileW(wpath, &fdFile)) == INVALID_HANDLE_VALUE) {
		fprintf(stderr, "Path not found: [%s]\n", sDir);
		g_free(wpath);
		return 1;
	}
	g_free(wpath);

	do {
		//Is the entity a File or Folder?
		if (!(fdFile.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
			gchar *cFileName = g_utf16_to_utf8(fdFile.cFileName, -1, NULL, NULL, NULL);
			if (cFileName == NULL) {
				return 1;
			}
			sequence *seq = readseqfile(cFileName);
			if (seq != NULL) {
				strncpy(filename, cFileName, 255);
				free_sequence(seq, TRUE);
				gtk_combo_box_text_append_text(seqcombo, filename);
				if (sequence_name_to_select
						&& !strncmp(filename, sequence_name_to_select,
								strlen(filename))) {
					*index_of_seq_to_load = number_of_loaded_sequences;
				}
				++number_of_loaded_sequences;
			}
			g_free(cFileName);
		}
	} while (FindNextFileW(hFind, &fdFile)); //Find the next file.

	FindClose(hFind);

	return number_of_loaded_sequences;
}
#endif

/** This method populates the sequence combo box with the sequences found in the CWD.
 *  If only one sequence is found, or if a sequence whose name matches the
 *  possibly NULL argument is found, it is automatically selected, which triggers
 *  its loading
 *  @param sequence_name_to_select the name of the input sequence
 *  @return 0 if success
 */
int update_sequences_list(const char *sequence_name_to_select) {
	GtkComboBoxText *seqcombo;
	int number_of_loaded_sequences = 0;
	int index_of_seq_to_load = -1;
	char *seqname = NULL;

	// clear the previous list
	seqcombo = GTK_COMBO_BOX_TEXT(lookup_widget("sequence_list_combobox"));
	gtk_combo_box_text_remove_all(seqcombo);

	if (sequence_name_to_select) {
	       if (g_str_has_suffix(sequence_name_to_select, ".seq"))
		       seqname = strdup(sequence_name_to_select);
	       else {
		       seqname = malloc(strlen(sequence_name_to_select) + 5);
		       sprintf(seqname, "%s.seq", sequence_name_to_select);
	       }
	}

#ifdef _WIN32
	number_of_loaded_sequences = ListSequences(com.wd, seqname, seqcombo, &index_of_seq_to_load);
#else
	struct dirent **list;
	int i, n;

	n = scandir(com.wd, &list, 0, alphasort);
	if (n < 0)
		perror("scandir");

	for (i = 0; i < n; ++i) {
		char *suf;

		if ((suf = strstr(list[i]->d_name, ".seq")) && strlen(suf) == 4) {
			sequence *seq = readseqfile(list[i]->d_name);
			if (seq != NULL) {
				free_sequence(seq, TRUE);
				char *filename = list[i]->d_name;
				gtk_combo_box_text_append_text(seqcombo, filename);
				if (seqname && !strcmp(filename, seqname))
					index_of_seq_to_load = number_of_loaded_sequences;
				++number_of_loaded_sequences;
			}
		}
	}
	for (i = 0; i < n; i++)
		free(list[i]);
	free(list);
#endif

	if (seqname) free(seqname);

	if (!number_of_loaded_sequences) {
		fprintf(stderr, "No valid sequence found in CWD.\n");
		return -1;
	} else {
		fprintf(stdout, "Loaded %d %s\n", number_of_loaded_sequences,
				ngettext("sequence", "sequences", number_of_loaded_sequences));
	}

	if (number_of_loaded_sequences > 1 && index_of_seq_to_load < 0) {
		gtk_combo_box_popup(GTK_COMBO_BOX(seqcombo));
	} else if (index_of_seq_to_load >= 0)
		gtk_combo_box_set_active(GTK_COMBO_BOX(seqcombo), index_of_seq_to_load);
	else
		gtk_combo_box_set_active(GTK_COMBO_BOX(seqcombo), 0);
	return 0;
}

/**
 * Expands the ~ in filenames
 * @param[in] filename input filename
 * @param[in] size maximum size of the filename
 */
void expand_home_in_filename(char *filename, int size) {
	if (filename[0] == '~' && filename[1] == '\0')
		strcat(filename, G_DIR_SEPARATOR_S);
	int len = strlen(filename);
	if (len < 2)
		return;		// not very necessary now with the first line
	if (filename[0] == '~' && filename[1] == G_DIR_SEPARATOR) {
		const gchar *homepath = g_get_home_dir();
		int j, homelen = strlen(homepath);
		if (len + homelen > size - 1) {
			siril_log_message(_("Filename is too long, not expanding it\n"));
			return;
		}
		for (j = len; j > 0; j--)		// edit in place
			filename[j + homelen - 1] = filename[j];
		// the -1 above is tricky: it's the removal of the ~ character from
		// the original string
		memcpy(filename, homepath, homelen);
	}
}

/**
 * Tries to get normalized value of a fit image. Make assumption that
 * an image with no values greater than 2^8 comes from 8-bit images
 * @param fit input FITS image
 * @return 255 or 65535 if 8- or 16-bit image
 */
double get_normalized_value(fits *fit) {
	if (fit->type == DATA_USHORT) {
		image_find_minmax(fit);
		if (fit->maxi <= UCHAR_MAX_DOUBLE)
			return UCHAR_MAX_DOUBLE;
		return USHRT_MAX_DOUBLE;
	}
	if (fit->type == DATA_FLOAT) {
		return 1.0;
	}
	return -1.0;
}

/**
 * append a string to the end of an existing string
 * @param data original string
 * @param newdata suffix to add
 * @return a new string that should be freed when no longer needed
 */
gchar* str_append(gchar** data, const gchar* newdata) {
	gchar* p;
	int len = (*data ? strlen(*data) : 0);
	if ((p = g_try_realloc(*data, len + strlen(newdata) + 1)) == NULL) {
		g_free(p);
		PRINT_ALLOC_ERR;
		return NULL;
	}
	*data = p;
	g_strlcpy(*data + len, newdata, len + strlen(newdata));
	return *data;
}

/**
 * Cut a base name to 120 characters and add a trailing underscore if needed.
 * WARNING: may return a newly allocated string and free the argument
 * @param root the original base name
 * @param can_free allow root to be freed in case a new string is allocated
 * @return a string ending with trailing underscore
 */
char *format_basename(char *root, gboolean can_free) {
	int len = strlen(root);
	if (len > 120) {
		root[120] = '\0';
		len = 120;
	}
	if (root[len - 1] == '-' || root[len - 1] == '_') {
		return root;
	}

	char *appended = malloc(len + 2);
	sprintf(appended, "%s_", root);
	if (can_free)
		free(root);
	return appended;
}

/**
 * Computes slope using low and high values
 * @param lo low value
 * @param hi high value
 * @return the computed slope
 */
float compute_slope(WORD *lo, WORD *hi) {
	if (sequence_is_loaded() && !single_image_is_loaded()) {
		*hi = com.seq.layers[RLAYER].hi;
		*lo = com.seq.layers[RLAYER].lo;
	}
	else {
		*hi = com.uniq->layers[RLAYER].hi;
		*lo = com.uniq->layers[RLAYER].lo;
	}
	return UCHAR_MAX_SINGLE / (float) (*hi - *lo);
}

/**
 * Try to get file info, i.e width and height
 * @param filename name of the file
 * @param pixbuf
 * @return a newly allocated and formated string containing dimension information or NULL
 */
gchar* siril_get_file_info(const gchar *filename, GdkPixbuf *pixbuf) {
	int width, height;
	int n_channel = 0;

	GdkPixbufFormat *pixbuf_file_info = gdk_pixbuf_get_file_info(filename,
			&width, &height);

	if (pixbuf) {
		n_channel = gdk_pixbuf_get_n_channels(pixbuf);
	}

	if (pixbuf_file_info != NULL) {
		/* Pixel size of image: width x height in pixel */
		return g_strdup_printf("%d x %d %s\n%d %s", width, height,
				ngettext("pixel", "pixels", height), n_channel,
				ngettext("channel", "channels", n_channel));
	}
	return NULL;
}

/**
 * Truncate a string str to not exceed an length of size
 * @param str the string to be truncated
 * @param size maximum size of the string
 * @return the truncated size starting by "..." and followed by "/"
 * if possible
 */
gchar *siril_truncate_str(gchar *str, gint size) {
	GString *trunc_str = g_string_new(str);
	gint len = strlen(str);

	if (len > size) {
		gint pos = len - size;
		/* locate first "/" */
		char *ptr = strchr(str + pos, G_DIR_SEPARATOR);
		if (ptr != NULL) {
			pos = ptr - str;
		}
		trunc_str = g_string_erase(trunc_str, 0, pos);
		trunc_str = g_string_prepend(trunc_str, "...");
	}
	return g_string_free(trunc_str, FALSE);
}

/**
 *
 * @param list
 * @param arg_count
 * @return
 */
char **glist_to_array(GList *list, int *arg_count) {
	int count;
	if (arg_count && *arg_count > 0)
		count = *arg_count;
	else {
		count = g_list_length(list);
		if (arg_count)
			*arg_count = count;
	}
	char **array = malloc(count * sizeof(char *));
	if (!array) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	GList *orig_list = list;
	for (int i = 0; i < count && list; list = list->next, i++)
		array[i] = g_strdup(list->data);
	g_list_free_full(orig_list, g_free);
	return array;
}

/**
 *
 * @param uri_string
 * @return a new allocated string
 */
gchar* url_cleanup(const gchar *uri_string) {
	GString *copy;
	const gchar *end;

	/* Skip leading whitespace */
	while (g_ascii_isspace(*uri_string))
		uri_string++;

	/* Ignore trailing whitespace */
	end = uri_string + strlen(uri_string);
	while (end > uri_string && g_ascii_isspace(*(end - 1)))
		end--;

	/* Copy the rest, encoding unencoded spaces and stripping other whitespace */
	copy = g_string_sized_new(end - uri_string);
	while (uri_string < end) {
		if (*uri_string == ' ')
			g_string_append(copy, "%20");
		else if (g_ascii_isspace(*uri_string))
			; // @suppress("Suspicious semicolon")
		else
			g_string_append_c(copy, *uri_string);
		uri_string++;
	}

	return g_string_free(copy, FALSE);
}
