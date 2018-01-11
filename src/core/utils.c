/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2018 team free-astro (see more in AUTHORS file)
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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <dirent.h>
#ifdef WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#endif
#if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#include <sys/param.h>		// define or not BSD macro
#endif
#if (defined(__APPLE__) && defined(__MACH__))
#include <mach/task.h>
#include <mach/mach_init.h>
#include <mach/mach_types.h>
#include <mach/mach_host.h>
#include <sys/sysctl.h>
#include <mach/vm_statistics.h>
#endif
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "io/conversion.h"
#include "io/sequence.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "io/single_image.h"

#ifdef HAVE_SYS_STATVFS_H
#include <sys/statvfs.h>
#endif
#if HAVE_SYS_VFS_H
#include <sys/vfs.h>
#elif HAVE_SYS_MOUNT_H
#if HAVE_SYS_PARAM_H
#include <sys/param.h>
#endif
#include <sys/mount.h>
#endif

/**
 * Round double value to an integer
 * @param x value to round
 * @return an integer
 */
int round_to_int(double x) {
	if (x <= INT_MIN + 0.5) return INT_MIN;
	if (x >= INT_MAX - 0.5) return INT_MAX;
	if (x >= 0.0)
		return (int) (x + 0.5);
	return (int) (x - 0.5);
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
		return (int) (x + 0.5f);
	return (int) (x - 0.5f);
}

/**
 * Round double value to a WORD
 * @param x value to round
 * @return a WORD
 */
WORD round_to_WORD(double x) {
	if (x <= 0.0)
		return (WORD) 0;
	if (x > USHRT_MAX_DOUBLE)
		return USHRT_MAX;
	return (WORD) (x + 0.5);
}

/**
 * Round double value to a BYTE
 * @param x value to round
 * @return a BYTE
 */
BYTE round_to_BYTE(double x) {
	if (x <= 0.0)
		return (BYTE) 0;
	if (x > UCHAR_MAX_DOUBLE)
		return UCHAR_MAX;
	return (BYTE) (x + 0.5);
}

/**
 * convert double value to a BYTE
 * @param x value to convert
 * @return a BYTE
 */
BYTE conv_to_BYTE(double x) {
	if (x == 0.0)
		return (BYTE) 0;
	if (x == USHRT_MAX_DOUBLE)
		return UCHAR_MAX;
	x = ((x / USHRT_MAX_DOUBLE) * UCHAR_MAX_DOUBLE);
	return((BYTE) (x));
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
 * Converts a string which is in the encoding used by GLib for filenames into a UTF-8 string.
 * Note that on Windows GLib uses UTF-8 for filenames; on other platforms,
 * this function indirectly depends on the current locale.
 * g_free the result when not needed anymore.
 * @param filename input filename
 * @return The converted string, or "<charset conversion error>" on an error.
 */
char *f2utf8(const char *filename) {
	char *utf8;

	if (!(utf8 = g_filename_to_utf8(filename, -1, NULL, NULL, NULL)))
		utf8 = g_strdup("<charset conversion error>");

	return (utf8);
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
	const char *dot = strrchr(filename, '.');
	if (!dot || dot == filename)
		return NULL;
	return dot + 1;
}

/**
 * Tests whether the given file is either regular or a symlink
 * @param filename input
 * @return 1 if file is readable (not actually opened to verify)
 */
int is_readable_file(const char *filename) {
	struct stat sts;
	if (g_stat(filename, &sts))
		return 0;
	if (S_ISREG (sts.st_mode)
#ifndef WIN32
			|| S_ISLNK(sts.st_mode)
#endif
	)
		return 1;
	return 0;
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

static GUserDirectory sdir[] = { G_USER_DIRECTORY_PICTURES,
		G_USER_DIRECTORY_DOCUMENTS };
/** This function tries to set a startup directory. It first looks at the "Pictures" directory,
 *  then if it does not exist, the "Document" one, Finally, if it fails on some UNIX systems
 *  the dir is set to the home directory.
 *  @return a working directory path if success, NULL if error
 */
gchar *siril_get_startup_dir() {
	const gchar *dir = NULL;
	gchar *startup_dir = NULL;
	gint i = 0;
	size_t size;

	size = sizeof(sdir) / sizeof(GUserDirectory);

	while (dir == NULL && i < size) {
		dir = g_get_user_special_dir(sdir[i]);
		i++;
	}
	/* Not every platform has a directory for these logical id */
	if (dir == NULL) {
		dir = g_get_home_dir();
	}
	if (dir)
		startup_dir = g_strdup(dir);
	return startup_dir;
}

/** Try to change the CWD to the argument, absolute or relative.
 *  If success, the new CWD is written to com.wd
 *  @param[in] dir absolute or relative path we want to set as cwd
 *  @param[out] err error message when return value is different of 1. Can be NULL if message is not needed.
 *  @return 0 if success, any other values for error
 */
int changedir(const char *dir, gchar **err) {
	gchar *error = NULL;
	int retval = 0;

	if (dir == NULL || dir[0] == '\0') {
		error = siril_log_message(_("Unknown error\n"));
		retval = -1;
	} else if (!g_file_test(dir, G_FILE_TEST_EXISTS)) {
		error = siril_log_message(_("No such file or directory\n"));
		retval = 2;
	} else if (!g_file_test(dir, G_FILE_TEST_IS_DIR)) {
		error = siril_log_message(_("\"%s\" is not a directory\n"), dir);
		retval = 3;
	} else if (g_access(dir, W_OK)) {
		error = siril_log_color_message(_("You don't have permission "
				"to write in this directory: %s\n"), "red", dir);
		retval = 4;
	} else {
		if (!g_chdir(dir)) {
			/* do we need to search for sequences in the directory now? We still need to
			 * press the check seq button to display the list, and this is also done there. */
			/* check_seq();
			 update_sequence_list();*/
			g_free(com.wd);
			com.wd = g_get_current_dir();
			siril_log_message(_("Setting CWD (Current "
					"Working Directory) to '%s'\n"), com.wd);
			set_GUI_CWD();
			update_used_memory();
			retval = 0;
		} else {
			error = siril_log_message(_("Could not change "
					"directory to '%s'.\n"), dir);
			retval = 1;
		}
	}
	if (err) {
		*err = error;
	}
	return retval;
}

#ifdef WIN32
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
	struct dirent **list;
	int number_of_loaded_sequences = 0;
	int index_of_seq_to_load = -1;

	// clear the previous list
	seqcombo = GTK_COMBO_BOX_TEXT(
			gtk_builder_get_object(builder, "sequence_list_combobox"));
	gtk_combo_box_text_remove_all(seqcombo);

#ifdef WIN32
	number_of_loaded_sequences = ListSequences(com.wd, sequence_name_to_select, seqcombo, &index_of_seq_to_load);
#else
	int i, n;
	char filename[256];

	n = scandir(com.wd, &list, 0, alphasort);
	if (n < 0)
		perror("scandir");

	for (i = 0; i < n; ++i) {
		char *suf;

		if ((suf = strstr(list[i]->d_name, ".seq")) && strlen(suf) == 4) {
			sequence *seq = readseqfile(list[i]->d_name);
			if (seq != NULL) {
				strncpy(filename, list[i]->d_name, 255);
				free_sequence(seq, TRUE);
				gtk_combo_box_text_append_text(seqcombo, filename);
				if (sequence_name_to_select
						&& !strncmp(filename, sequence_name_to_select,
								strlen(filename))) {
					index_of_seq_to_load = number_of_loaded_sequences;
				}
				++number_of_loaded_sequences;
			}
		}
	}
	for (i = 0; i < n; i++)
		free(list[i]);
	free(list);
#endif

	if (!number_of_loaded_sequences) {
		fprintf(stderr, "No valid sequence found in CWD.\n");
		return -1;
	} else {
		fprintf(stdout, "Loaded %d sequence(s)\n", number_of_loaded_sequences);
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
 * Find the space remaining in a directory, in bytes. A double for >32bit
 * problem avoidance. <0 for error.
 * @param name the path of the directory to be tested
 * @return the disk space remaining in bytes, or a value less than 0 if error
 */
#ifdef HAVE_SYS_STATVFS_H
static double find_space(const gchar *name) {
	struct statvfs st;
	double sz;

	if (statvfs (name, &st))
		/* Set to error value.
		 */
		sz = -1;
	else
		sz = (double) st.f_frsize * st.f_bavail;

	return (sz);
}
#elif (HAVE_SYS_VFS_H || HAVE_SYS_MOUNT_H)
static double find_space(const gchar *name) {
	struct statfs st;
	double sz;

	if (statfs (name, &st))
		sz = -1;
	else
		sz = (double) st.f_bsize * st.f_bavail;

	return (sz);
}
#elif defined WIN32
static double find_space(const gchar *name) {
	ULARGE_INTEGER avail;
	double sz;

	gchar *localdir = g_path_get_dirname(name);
	wchar_t *wdirname = g_utf8_to_utf16(localdir, -1, NULL, NULL, NULL);

	if (!GetDiskFreeSpaceExW(wdirname, &avail, NULL, NULL))
		sz = -1;
	else
		sz = (double) avail.QuadPart;

	g_free(localdir);
	g_free(wdirname);
	return (sz);
}
#else
static double find_space(const gchar *name) {
	return (-1);
}
#endif /*HAVE_SYS_STATVFS_H*/

#if defined(__linux__)
static unsigned long update_used_RAM_memory() {
	unsigned long size, resident, share, text, lib, data, dt;
	static int page_size_in_k = 0;
	const char* statm_path = "/proc/self/statm";
	FILE *f = fopen(statm_path, "r");

	if (page_size_in_k == 0) {
		page_size_in_k = getpagesize() / 1024;
	}
	if (!f) {
		perror(statm_path);
		return 0UL;
	}
	if (7 != fscanf(f, "%lu %lu %lu %lu %lu %lu %lu",
			&size, &resident, &share, &text, &lib, &data, &dt)) {
		perror(statm_path);
		fclose(f);
		return 0UL;
	}
	fclose(f);
	return (resident * page_size_in_k);
}
#elif (defined(__APPLE__) && defined(__MACH__))
static unsigned long update_used_RAM_memory() {
	struct task_basic_info t_info;

	mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
	task_info(current_task(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
	return ((unsigned long) t_info.resident_size / 1024UL);
}
#elif defined(BSD) /* BSD (DragonFly BSD, FreeBSD, OpenBSD, NetBSD). In fact, it could work with linux */
static unsigned long update_used_RAM_memory() {
	struct rusage usage;

	getrusage(RUSAGE_SELF, &usage);
	return ((unsigned long) usage.ru_maxrss);
}
#elif defined(WIN32) /* Windows */
static unsigned long update_used_RAM_memory() {
    PROCESS_MEMORY_COUNTERS memCounter;
    
	if (GetProcessMemoryInfo(GetCurrentProcess(), &memCounter, sizeof(memCounter)))
        return (memCounter.WorkingSetSize / 1024UL);
	return 0UL;
}
#else
static unsigned long update_used_RAM_memory() {
	return 0UL;
}
#endif

/**
 * Updates RAM memory used by siril, available free disk space
 * and displays information on the control window.
 */
void update_used_memory() {
	unsigned long ram;
	double freeDisk;

	ram = update_used_RAM_memory();
	freeDisk = find_space(com.wd);
	/* update GUI */
	set_GUI_MEM(ram);
	set_GUI_DiskSpace(freeDisk);
}

/**
 * Test if there is enough free disk space by returning the difference
 * between available free disk space and the size given in parameters
 * @param seq_size size to be tested
 * @return a value greater than 0 if there is enough disk space, a value
 * less than 0 otherwise. The function returns -1 if an error occurs.
 */
double test_available_space(double seq_size) {
	double freeDisk;

	freeDisk = find_space(com.wd);
	if ((freeDisk < 0) || (seq_size < 0)) {
		return -1;
	}
	return (freeDisk - seq_size);
}

/**
 * Gets available memory for stacking process
 * @return available memory in MB, 2048 if it fails.
 */
#if defined(__linux__)
int get_available_memory_in_MB() {
	int mem = 2048; /* this is the default value if we can't retrieve any values */
	FILE* fp = fopen("/proc/meminfo", "r");
	if (fp != NULL) {
		size_t bufsize = 1024 * sizeof(char);
		char* buf = (char*) malloc(bufsize);
		long value = -1L;
		while (getline(&buf, &bufsize, fp) >= 0) {
			if (strncmp(buf, "MemAvailable", 12) != 0)
				continue;
			sscanf(buf, "%*s%ld", &value);
			break;
		}
		fclose(fp);
		free((void*) buf);
		if (value != -1L)
			mem = (int) (value / 1024L);
	}
	return mem;
}
#elif (defined(__APPLE__) && defined(__MACH__))
int get_available_memory_in_MB() {
	int mem = 2048; /* this is the default value if we can't retrieve any values */
	vm_size_t page_size;
	mach_port_t mach_port;
	mach_msg_type_number_t count;
	vm_statistics64_data_t vm_stats;

	mach_port = mach_host_self();
	count = sizeof(vm_stats) / sizeof(natural_t);
	if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
			KERN_SUCCESS == host_statistics64(mach_port, HOST_VM_INFO,
					(host_info64_t)&vm_stats, &count))	{

		int64_t unused_memory = ((int64_t)vm_stats.free_count +
				(int64_t)vm_stats.inactive_count +
				(int64_t)vm_stats.wire_count) * (int64_t)page_size;

		mem = (int) ((unused_memory) / (1024 * 1024));
	}
	return mem;
}
#elif defined(BSD) /* BSD (DragonFly BSD, FreeBSD, OpenBSD, NetBSD). ----------- */
int get_available_memory_in_MB() {
	int mem = 2048; /* this is the default value if we can't retrieve any values */
	FILE* fp = fopen("/var/run/dmesg.boot", "r");
	if (fp != NULL) {
		size_t bufsize = 1024 * sizeof(char);
		char* buf = (char*) malloc(bufsize);
		long value = -1L;
		while (getline(&buf, &bufsize, fp) >= 0) {
			if (strncmp(buf, "avail memory", 12) != 0)
				continue;
			sscanf(buf, "%*s%*s%*s%ld", &value);
			break;
		}
		fclose(fp);
		free((void*) buf);
		if (value != -1L)
			mem = (int) (value / 1024L);
	}
	return mem;
}
#elif defined(WIN32) /* Windows */
int get_available_memory_in_MB() {
	int mem = 2048; /* this is the default value if we can't retrieve any values */
	MEMORYSTATUSEX memStatusEx = {0};
	memStatusEx.dwLength = sizeof(MEMORYSTATUSEX);
	const DWORD dwMBFactor = 1024 * 1024;
	DWORDLONG dwTotalPhys = memStatusEx.ullTotalPhys / dwMBFactor;
	if (dwTotalPhys > 0)
		mem = (int) dwTotalPhys;
	return mem;
}
#else
int get_available_memory_in_MB() {
	printf("Siril failed to get available free RAM memory\n");
	return 2048;
}
#endif

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
		strncpy(filename, homepath, homelen);
	}
}

/**
 * Tries to get normalized value of a fit image. Make assumption that
 * an image with no values greater than 2^8 comes from 8-bit images
 * @param fit input FITS image
 * @return 255 or 65535 if 8- or 16-bit image
 */
WORD get_normalized_value(fits *fit) {
	image_find_minmax(fit, 0);
	if (fit->maxi <= UCHAR_MAX)
		return UCHAR_MAX;
	return USHRT_MAX;
}

/**
 * This function reads a text file and displays it in the
 * show_data_dialog
 * @param path filename to display
 * @param title text shown as dialog title
 */
void read_and_show_textfile(char *path, char *title) {
	char line[64] = "";
	char txt[1024] = "";

	FILE *f = g_fopen(path, "r");
	if (!f) {
		show_dialog(_("File not found"), _("Error"), "gtk-dialog-error");
		return;
	}
	while (fgets(line, sizeof(line), f) != NULL)
		strcat(txt, line);
	show_data_dialog(txt, title);
	fclose(f);
}

/**
 * Switch the two parameters of the function:
 * Useful in Dynamic PSF (PSF.c)
 * @param a first parameter to switch
 * @param b second parameter to switch
 */
void swap_param(double *a, double *b) {
	double tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
}

/**
 * In-place quick sort of array of double a of size n
 * @param a array to sort
 * @param n size of the array
 */
void quicksort_d(double *a, int n) {
	if (n < 2)
		return;
	double p = a[n / 2];
	double *l = a;
	double *r = a + n - 1;
	while (l <= r) {
		if (*l < p) {
			l++;
			continue;
		}
		if (*r > p) {
			r--;
			continue; // we need to check the condition (l <= r) every time we change the value of l or r
		}
		double t = *l;
		*l++ = *r;
		*r-- = t;
	}
	quicksort_d(a, r - a + 1);
	quicksort_d(l, a + n - l);
}

/**
 * In-place quick sort of array of WORD a of size n
 * @param a array to sort
 * @param n size of the array
 */
void quicksort_s(WORD *a, int n) {
	if (n < 2)
		return;
	WORD p = a[n / 2];
	WORD *l = a;
	WORD *r = a + n - 1;
	while (l <= r) {
		if (*l < p) {
			l++;
			continue;
		}
		if (*r > p) {
			r--;
			continue; // we need to check the condition (l <= r) every time we change the value of l or r
		}
		WORD t = *l;
		*l++ = *r;
		*r-- = t;
	}
	quicksort_s(a, r - a + 1);
	quicksort_s(l, a + n - l);
}

/**
 * Removes extension of the filename
 * @param filename file path with extension
 * @return filename without extension
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
char* str_append(char** data, const char* newdata) {
	char* p;
	int len = (*data ? strlen(*data) : 0);
	if ((p = realloc(*data, len + strlen(newdata) + 1)) == NULL) {
		free(p);
		printf("str_append: error allocating data\n");
		return NULL;
	}
	*data = p;
	strcpy(*data + len, newdata);
	return *data;
}

/**
 * Cut a base name to 120 characters and add a trailing underscore if needed.
 * WARNING: may return a newly allocated string and free the argument
 * @param root the original base name
 * @return a string ending with trailing underscore
 */
char *format_basename(char *root) {
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
	free(root);
	return appended;
}

/**
 * Computes slop using low and high values
 * @param lo low value
 * @param hi high value
 * @return the computed slope
 */
float computePente(WORD *lo, WORD *hi) {
	float pente;

	if (sequence_is_loaded() && !single_image_is_loaded()) {
		*hi = com.seq.layers[RLAYER].hi;
		*lo = com.seq.layers[RLAYER].lo;
	}
	else {
		*hi = com.uniq->layers[RLAYER].hi;
		*lo = com.uniq->layers[RLAYER].lo;
	}

	pente = UCHAR_MAX_SINGLE / (float) (*hi - *lo);

	return pente;
}

static const gchar *checking_css_filename() {
	printf(_("Checking GTK version ... GTK-%d.%d\n"), GTK_MAJOR_VERSION, GTK_MINOR_VERSION);
	if ((GTK_MAJOR_VERSION >= 3) && (GTK_MINOR_VERSION >= 20))
		return "gtk.css";
	else if ((GTK_MAJOR_VERSION >= 3) && (GTK_MINOR_VERSION < 20))
		return "gtk_old.css";
	else {
		return NULL;
	}
}

/**
 * Loads the css sheet
 * @param path path of the file being loaded
 */
void load_css_style_sheet (char *path) {
	GtkCssProvider *css_provider;
	GdkDisplay *display;
	GdkScreen *screen;
	gchar *CSSFile;
	const gchar *css_filename;

	css_filename = checking_css_filename();
	if (css_filename == NULL) {
		printf(_("The version of GTK does not match requirements: (GTK-%d.%d)\n"), GTK_MAJOR_VERSION, GTK_MINOR_VERSION);
		exit(1);
	}

	CSSFile = g_build_filename (path, css_filename, NULL);
	if (!g_file_test (CSSFile, G_FILE_TEST_EXISTS)) {
		g_error (_("Unable to load CSS style sheet file: %s. Please reinstall %s\n"), CSSFile, PACKAGE);
	}
	else {
		css_provider = gtk_css_provider_new();
		display = gdk_display_get_default();
		screen = gdk_display_get_default_screen(display);
		gtk_style_context_add_provider_for_screen(screen,
				GTK_STYLE_PROVIDER(css_provider),
				GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);
		gtk_css_provider_load_from_path(css_provider, CSSFile, NULL);
		fprintf(stdout, _("Successfully loaded '%s'\n"), CSSFile);
		g_object_unref (css_provider);
	}
	g_free(CSSFile);
}

/**
 * From a datetime it computes the Julian date needed in photometry
 * (code borrowed from muniwin)
 * @param dt timestamp in datetime format
 * @return the Julian date
 */
double encodeJD(dateTime dt) {
	double jd1;
	int before, d1, d2;

	/* Check date and time */
	if (dt.day <= 0 || dt.year <= 0 || dt.month <= 0)
		return 0;

	/* Compute Julian date from input citizen year, month and day. */
	/* Tested for YEAR>0 except 1582-10-07/15 */
	if (dt.year > 1582) {
		before = 0;
	} else if (dt.year < 1582) {
		before = 1;
	} else if (dt.month > 10) {
		before = 0;
	} else if (dt.month < 10) {
		before = 1;
	} else if (dt.day >= 15) {
		before = 0;
	} else {
		before = 1;
	}
	if (dt.month <= 2) {
		d1 = (int) (365.25 * (dt.year - 1));
		d2 = (int) (30.6001 * (dt.month + 13));
	} else {
		d1 = (int) (365.25 * (dt.year));
		d2 = (int) (30.6001 * (dt.month + 1));
	}
	jd1 = 1720994.5 + d1 + d2 + dt.day;
	jd1 += 1.0 * dt.hour / 24;
	jd1 += 1.0 * dt.min / 1440.0;
	jd1 += 1.0 * dt.sec / 86400.0;
	jd1 += 1.0 * dt.ms / 86400000.0;

	if (before) {
		if (dt.year < 0)
			return jd1 - 1;
		else
			return jd1;
	} else {
		return jd1 + 2 - (dt.year / 100) + (dt.year / 400);
	}
}
