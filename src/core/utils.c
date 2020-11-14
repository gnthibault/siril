/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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
#include "gui/callbacks.h"
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
 *  Looks whether the string str ends with ending. This is case insensitive
 *  @param str the string to check
 *  @param ending the suffix to look for
 *  @return TRUE if str ends with ending
 */
gboolean ends_with(const char *str, const char *ending) {
	if (!str || str[0] == '\0')
		return FALSE;
	if (!ending || ending[0] == '\0')
		return TRUE;
	int ending_len = strlen(ending);
	int str_len = strlen(str);
	if (ending_len > str_len)
		return FALSE;
	return !strncasecmp(str + str_len - ending_len, ending, ending_len);
}

/**
 *  Searches for an extension '.something' in filename from the end
 *  @param filename input filename
 *  @return the index of the first '.' found
 */
int get_extension_index(const char *filename) {
	int i;
	if (filename == NULL || filename[0] == '\0')
		return -1;
	i = strlen(filename) - 1;
	do {
		if (filename[i] == '.')
			return i;
		i--;
	} while (i > 0);
	return -1;
}

/**
 * Get the extension of a file, without the dot.
 * @param filename input filename
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
#if GLIB_CHECK_VERSION(2,58,0)
				new_dir = g_canonicalize_filename(dir, com.wd);
#else
    			new_dir = g_build_filename(dir, com.wd, NULL);
#endif
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
	       if (ends_with(sequence_name_to_select, ".seq"))
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
 * Removes extension of the filename
 * @param filename file path with extension
 * @return newly allocated filename without extension
 */
char *remove_ext_from_filename(const char *filename) {
	size_t filelen;
	const char *p;
	char *file = NULL;

	p = strrchr(filename, '.');

	if (p == NULL) {
		file = malloc(1);
		file[0] = '\0';
		return file;
	}

	filelen = p - filename;
	file = malloc(filelen + 1);
	strncpy(file, filename, filelen);
	file[filelen] = '\0';

	return file;
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
 * Create a popover with icon and text
 * @param widget is the parent widget where the popover arises from
 * @param text will be shown in the popover
 * @return the GtkWidget of popover
 */
GtkWidget* popover_new(GtkWidget *widget, const gchar *text) {
	GtkWidget *popover, *box, *image, *label;

	popover = gtk_popover_new(widget);
	label = gtk_label_new(NULL);
	box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	image = gtk_image_new_from_icon_name("dialog-information-symbolic",
			GTK_ICON_SIZE_DIALOG);

	gtk_label_set_markup(GTK_LABEL(label), text);
	gtk_label_set_line_wrap(GTK_LABEL(label), TRUE);
	gtk_label_set_max_width_chars(GTK_LABEL(label), 64);

	gtk_box_pack_start(GTK_BOX(box), image, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(box), label, FALSE, FALSE, 0);
	gtk_container_add(GTK_CONTAINER(popover), box);

	/* make all sensitive in case where parent is not */
	gtk_widget_set_sensitive(label, TRUE);
	gtk_widget_set_sensitive(box, TRUE);
	gtk_widget_set_sensitive(popover, TRUE);

	gtk_widget_show_all(box);

	return popover;
}

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
