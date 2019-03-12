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
 */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

// TODO: remove unused
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <dirent.h>
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#include <direct.h>
#include <shlobj.h>
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
#include <string.h>
#include <assert.h>
#include <fitsio.h>
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

#include "core/siril.h"
#include "gui/callbacks.h"

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
#elif defined _WIN32
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

#if defined(__linux__) || defined(__CYGWIN__)
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
#elif defined(_WIN32) /* Windows */
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
#if defined(__linux__) || defined(__CYGWIN__)
int get_available_memory_in_MB() {
	int mem = 2048; /* this is the default value if we can't retrieve any values */
	FILE* fp = fopen("/proc/meminfo", "r");
	if (fp != NULL) {
		size_t bufsize = 1024 * sizeof(char);
		gchar *buf = g_new(gchar, bufsize);
		long value = -1L;
		while (getline(&buf, &bufsize, fp) >= 0) {
			if (strncmp(buf, "MemAvailable", 12) != 0)
				continue;
			sscanf(buf, "%*s%ld", &value);
			break;
		}
		fclose(fp);
		g_free(buf);
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
		gchar *buf = g_new(gchar, bufsize);
		long value = -1L;
		while (getline(&buf, &bufsize, fp) >= 0) {
			if (strncmp(buf, "avail memory", 12) != 0)
				continue;
			sscanf(buf, "%*s%*s%*s%ld", &value);
			break;
		}
		fclose(fp);
		g_free(buf);
		if (value != -1L)
			mem = (int) (value / 1024L);
	}
	return mem;
}
#elif defined(_WIN32) /* Windows */
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
	fprintf(stderr, "Siril failed to get available free RAM memory\n");
	return 2048;
}
#endif

/**
 *
 * @param filename
 * @param size
 */
#ifdef _WIN32
/* stolen from gimp which in turn stole it from glib 2.35 */
gchar *get_special_folder(int csidl) {
	wchar_t path[MAX_PATH + 1];
	HRESULT hr;
	LPITEMIDLIST pidl = NULL;
	BOOL b;
	gchar *retval = NULL;

	hr = SHGetSpecialFolderLocation(NULL, csidl, &pidl);
	if (hr == S_OK) {
		b = SHGetPathFromIDListW(pidl, path);
		if (b)
			retval = g_utf16_to_utf8(path, -1, NULL, NULL, NULL);
		CoTaskMemFree(pidl);
	}

	return retval;
}
#endif

/**
 * Check how many files a process can have open and try to extend the limit if possible.
 * The max files depends of the Operating System and of cfitsio (NMAXFILES)
 * @param nb_frames number of file processed
 * @param nb_allowed_file the maximum of file that can be opened
 * @return TRUE if the system can open all the files, FALSE otherwise
 */
gboolean allow_to_open_files(int nb_frames, int *nb_allowed_file) {
	int open_max, maxfile, MAX_NO_FILE_CFITSIO, MAX_NO_FILE;
	float version;

	/* get the limit of cfitsio */
	fits_get_version(&version);
	MAX_NO_FILE_CFITSIO = (version < 3.45) ? 1000 : 10000;

	/* get the OS limit and extend it if possible */
#ifdef _WIN32
	MAX_NO_FILE = min(MAX_NO_FILE_CFITSIO, 2048);
	open_max = _getmaxstdio();
	if (open_max < MAX_NO_FILE) {
		/* extend the limit to 2048 if possible
		 * 2048 is the maximum on WINDOWS */
		_setmaxstdio(MAX_NO_FILE);
		open_max = _getmaxstdio();
	}
#else
	struct rlimit rlp;

/* we first set the limit to the CFITSIO limit */
	MAX_NO_FILE = MAX_NO_FILE_CFITSIO;
	if (getrlimit(RLIMIT_NOFILE, &rlp) == 0) {
		MAX_NO_FILE = (rlp.rlim_max == RLIM_INFINITY) ?
						MAX_NO_FILE_CFITSIO : rlp.rlim_max;

		if (rlp.rlim_cur != RLIM_INFINITY) {
			open_max = rlp.rlim_cur;
			MAX_NO_FILE = min(MAX_NO_FILE_CFITSIO, MAX_NO_FILE);
			if (open_max < MAX_NO_FILE) {
				rlp.rlim_cur = MAX_NO_FILE;
				/* extend the limit to NMAXFILES if possible */
				int retval = setrlimit(RLIMIT_NOFILE, &rlp);
				if (!retval) {
					getrlimit(RLIMIT_NOFILE, &rlp);
					open_max = rlp.rlim_cur;
				}
			}
		} else { // no soft limits
			open_max = MAX_NO_FILE;
		}
	} else {
		open_max = sysconf(_SC_OPEN_MAX); // if no success with getrlimit, try with sysconf
	}
#endif // _WIN32

	maxfile = min(open_max, MAX_NO_FILE);
	siril_debug_print("Maximum of files that will be opened=%d\n", maxfile);
	*nb_allowed_file = maxfile;

	return nb_frames < maxfile;
}

