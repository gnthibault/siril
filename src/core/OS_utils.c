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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <fcntl.h>

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#include <direct.h>
#include <shlobj.h>
#include <tchar.h>
#include <io.h>
#include <fcntl.h>
#else
#include <sys/resource.h>
#endif
#if defined(__unix__) || defined(OS_OSX)
#include <sys/param.h>		// define or not BSD macro
#endif
#ifdef OS_OSX
#include <mach/task.h>
#include <mach/mach_init.h>
#include <mach/mach_types.h>
#include <mach/mach_host.h>
#include <sys/sysctl.h>
#include <mach/vm_statistics.h>
#endif
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
#include "core/proto.h"
#include "gui/callbacks.h"
#include "gui/progress_and_log.h"
#include "gui/message_dialog.h"

#include "OS_utils.h"

/**
 * Find the space remaining in a directory, in bytes.
 * @param name the path of the directory to be tested
 * @return the disk space remaining in bytes, or a negative value if error
 */
#ifdef OS_X
static int64_t find_space(const gchar *name) {
	NSError *error;

	NSDictionary* fileAttributes = [[NSFileManager defaultManager] attributesOfFileSystemForPath:@name
	                                                                                   error:&error];
	int64_t freeSpace = [[fileAttributes objectForKey:NSFileSystemFreeSize] longLongValue];
	return freeSpace;
}
#elif HAVE_SYS_STATVFS_H
static int64_t find_space(const gchar *name) {
	struct statvfs st;
	int64_t available;
	if (statvfs (name, &st))
		return -1LL;
	available = st.f_bavail;        // force 64 bits
	return available * st.f_frsize;
}
#elif (HAVE_SYS_VFS_H || HAVE_SYS_MOUNT_H)
static int64_t find_space(const gchar *name) {
	struct statfs st;
	int64_t available;
	if (statfs (name, &st))
		return -1LL;
	available = st.f_bavail;        // force 64 bits
        return available * st.f_bsize;
}
#elif defined _WIN32
static int64_t find_space(const gchar *name) {
	ULARGE_INTEGER avail;
	int64_t sz;

	gchar *localdir = g_path_get_dirname(name);
	wchar_t *wdirname = g_utf8_to_utf16(localdir, -1, NULL, NULL, NULL);

	if (!GetDiskFreeSpaceExW(wdirname, &avail, NULL, NULL))
		sz = -1;
	else
		sz = avail.QuadPart;

	g_free(localdir);
	g_free(wdirname);
	return sz;
}
#else
static int64_t find_space(const gchar *name) {
	return -1LL;
}
#endif /*HAVE_SYS_STATVFS_H*/

#if defined(__linux__) || defined(__CYGWIN__)
static unsigned long long update_used_RAM_memory() {
	static gboolean initialized = FALSE;
	static long page_size;
	static gint fd = -1;
	gchar buffer[128];
	gint size;
	unsigned long long resident;
	unsigned long long shared;

	if (!initialized) {
		page_size = getpagesize() / 1024L;

		if (page_size > 0)
			fd = g_open("/proc/self/statm", O_RDONLY);

		initialized = TRUE;
	}

	if (fd < 0)
		return 0ULL;

	if (lseek(fd, 0, SEEK_SET))
		return 0ULL;

	size = read(fd, buffer, sizeof(buffer) - 1);

	if (size <= 0)
		return 0ULL;

	buffer[size] = '\0';

	if (sscanf(buffer, "%*u %llu %llu", &resident, &shared) != 2)
		return 0ULL;

	return (unsigned long long) (resident /*- shared*/) * page_size;
}
#elif defined(OS_OSX)
static unsigned long long update_used_RAM_memory() {
	struct task_basic_info t_info;

	mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
	task_info(current_task(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
	return ((unsigned long long) t_info.resident_size / 1024ULL);
}
#elif defined(BSD) /* BSD (DragonFly BSD, FreeBSD, OpenBSD, NetBSD). In fact, it could work with linux */
static unsigned long long update_used_RAM_memory() {
	struct rusage usage;

	getrusage(RUSAGE_SELF, &usage);
	return ((unsigned long long) usage.ru_maxrss);
}
#elif defined(_WIN32) /* Windows */
static unsigned long long update_used_RAM_memory() {
	PROCESS_MEMORY_COUNTERS memCounter;

	if (GetProcessMemoryInfo(GetCurrentProcess(), &memCounter, sizeof(memCounter)))
		return (memCounter.WorkingSetSize / 1024ULL);
	return 0ULL;
}
#else
static unsigned long long update_used_RAM_memory() {
	return 0ULL;
}
#endif


/**
 * Updates RAM memory used by siril, available free disk space
 * and displays information on the control window.
 * @return always return TRUE
 */
gboolean update_displayed_memory() {
	set_GUI_MEM(update_used_RAM_memory());
	set_GUI_DiskSpace((double)find_space(com.wd));
	return TRUE;
}

/**
 * From a number of bytes in input, returns a string (to be freed) comprehensible by a
 * human for this size, for example 1.5G instead of 1500000000
 * @param bytes
 * @return a formated string showing memory
 */
gchar *pretty_print_memory(int64_t bytes) {
	const char *units[] = { "", "k", "M", "G", "T", "P", "E", "Z", "Y" };
	int i = 0;
	double mem = (double)bytes;
	while (mem >= 1000.0 && i < sizeof units) {
		mem = mem / 1024.0;
		i++;
	}
	return g_strdup_printf("%.1f%s", mem, units[i]);
}

/**
 * Test if there is enough free disk space by returning the difference
 * in bytes between available free disk space and the size given as parameter
 * @param req_size available space to be tested
 * @return 0 if there is enough disk space, 1 otherwise, -1 on error.
 */
int test_available_space(int64_t req_size) {
	int64_t free_space = find_space(com.wd);
	if (free_space < 0 || req_size <= 0)
		return -1;

	if (req_size > free_space) {
		gchar *avail = pretty_print_memory(free_space);
		gchar *required = pretty_print_memory(req_size);
		gchar *missing = pretty_print_memory(req_size - free_space);
		char *msg = siril_log_message(_("Not enough free disk space to perform this operation: "
					"%sB available for %sB needed (missing %sB)\n"),
				avail, required, missing);
		queue_message_dialog(GTK_MESSAGE_ERROR, _("Not enough disk space"), msg);
		g_free(avail);
		g_free(required);
		g_free(missing);
		return 1;
	}
	siril_debug_print("Tested free space ok: %ld for %ld MB free\n",
			(long)(req_size / BYTES_IN_A_MB), (long)(free_space / BYTES_IN_A_MB));
	return 0;
}

/**
 * Gets available memory for stacking process
 * @return available memory in MB, 0 if it fails.
 */
#if defined(__linux__) || defined(__CYGWIN__)
int get_available_memory_in_MB() {
	static gboolean initialized = FALSE;
	static int64_t last_check_time = 0;
	static gint fd;
	static uint64_t available;
	static gboolean has_available = FALSE;
	int64_t time;

	if (!initialized) {
		fd = open("/proc/meminfo", O_RDONLY);

		initialized = TRUE;
	}

	if (fd < 0)
		return 0;

	/* we don't have a config option for limiting the swap size, so we simply
	 * return the free space available on the filesystem containing the swap
	 */

	time = g_get_monotonic_time();

	if (time - last_check_time >= G_TIME_SPAN_SECOND) {
		gchar buffer[512];
		gint size;
		gchar *str;

		last_check_time = time;

		has_available = FALSE;

		if (lseek(fd, 0, SEEK_SET))
			return 0;

		size = read(fd, buffer, sizeof(buffer) - 1);

		if (size <= 0)
			return 0;

		buffer[size] = '\0';

		str = strstr(buffer, "MemAvailable:");

		if (!str)
			return 0;

		available = strtoull(str + 13, &str, 0);

		if (!str)
			return 0;

		for (; *str; str++) {
			if (*str == 'k') {
				available <<= 10;
				break;
			} else if (*str == 'M') {
				available <<= 20;
				break;
			}
		}

		if (!*str)
			return 0;

		has_available = TRUE;
	}

	if (!has_available)
		return 0;

	return (int) (available / (uint64_t)BYTES_IN_A_MB);
}
#elif defined(OS_OSX)
int get_available_memory_in_MB() {
	int mem = 0; /* this is the default value if we can't retrieve any values */
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

		mem = (int) ((unused_memory) / BYTES_IN_A_MB);
	}
	return mem;
}
#elif defined(BSD) /* BSD (DragonFly BSD, FreeBSD, OpenBSD, NetBSD). ----------- */
int get_available_memory_in_MB() {
	int mem = 0; /* this is the default value if we can't retrieve any values */
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
	int mem = 0; /* this is the default value if we can't retrieve any values */
	MEMORYSTATUSEX memStatusEx = { 0 };
	memStatusEx.dwLength = sizeof(MEMORYSTATUSEX);
	if (GlobalMemoryStatusEx(&memStatusEx)) {
		mem = (int) (memStatusEx.ullAvailPhys / BYTES_IN_A_MB);
	}
	return mem;
}
#else
int get_available_memory_in_MB() {
	fprintf(stderr, "Siril failed to get available free RAM memory\n");
	return 0;
}
#endif

/**
 * Get max memory depending on memory management mode
 * @return return the max memory, and -1 for unlimited
 */
int get_max_memory_in_MB() {
	switch (com.stack.mem_mode) {
		default:
		case RATIO:
			return round_to_int(com.stack.memory_ratio *
					(double)get_available_memory_in_MB());

		case AMOUNT:
			return round_to_int(com.stack.memory_amount * 1024.0);
		case UNLIMITED:
			return -1;
	}
}

/**
 *
 * @param filename
 * @param size
 */
#ifdef _WIN32
/* stolen from gimp which in turn stole it from glib 2.35 */
gchar* get_special_folder(int csidl) {
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
	MAX_NO_FILE_CFITSIO = (version < 3.45f) ? 1000 : 10000;

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

SirilWidget *siril_file_chooser_open(GtkWindow *parent, GtkFileChooserAction action) {
	gchar *title;
	SirilWidget *w;
	if (action == GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER) {
		title = g_strdup(_("Select Folder"));
	} else {
		title = g_strdup(_("Open File"));
	}
	w = gtk_file_chooser_dialog_new(title, parent, action, _("_Cancel"),
			GTK_RESPONSE_CANCEL, _("_Open"), GTK_RESPONSE_ACCEPT,
			NULL);
	g_free(title);
	return w;
}

SirilWidget *siril_file_chooser_add(GtkWindow *parent, GtkFileChooserAction action) {
	return gtk_file_chooser_dialog_new(_("Add Files"), parent, action,
			_("_Cancel"), GTK_RESPONSE_CANCEL, _("_Add"), GTK_RESPONSE_ACCEPT,
			NULL);
}

SirilWidget *siril_file_chooser_save(GtkWindow *parent, GtkFileChooserAction action) {
	return gtk_file_chooser_dialog_new(_("Save File"), parent, action,
			_("_Cancel"), GTK_RESPONSE_CANCEL, _("_Save"), GTK_RESPONSE_ACCEPT,
			NULL);
}

gint siril_dialog_run(SirilWidget *widgetdialog) {
	return gtk_dialog_run(GTK_DIALOG(GTK_FILE_CHOOSER(widgetdialog)));
}

void siril_widget_destroy(SirilWidget *widgetdialog) {
	gtk_widget_destroy(widgetdialog);
}

#ifdef _WIN32
/* origin of sources: https://stackoverflow.com/questions/24171017/win32-console-application-that-can-open-windows */
int ReconnectIO(int OpenNewConsole) {
	int hConHandle;
	HANDLE lStdHandle;
	FILE *fp;
	int MadeConsole;

	MadeConsole = 0;
	if (!AttachConsole(ATTACH_PARENT_PROCESS)) {
		if (!OpenNewConsole)
			return 0;

		MadeConsole = 1;
		if (!AllocConsole())
			return 0;
	}

	// STDOUT to the console
	lStdHandle = GetStdHandle(STD_OUTPUT_HANDLE);
	hConHandle = _open_osfhandle((intptr_t) lStdHandle, _O_TEXT);
	fp = _fdopen(hConHandle, "w");
	*stdout = *fp;
	setvbuf( stdout, NULL, _IONBF, 0);

	// STDIN to the console
	lStdHandle = GetStdHandle(STD_INPUT_HANDLE);
	hConHandle = _open_osfhandle((intptr_t) lStdHandle, _O_TEXT);
	fp = _fdopen(hConHandle, "r");
	*stdin = *fp;
	setvbuf( stdin, NULL, _IONBF, 0);

	// STDERR to the console
	lStdHandle = GetStdHandle(STD_ERROR_HANDLE);
	hConHandle = _open_osfhandle((intptr_t) lStdHandle, _O_TEXT);
	fp = _fdopen(hConHandle, "w");
	*stderr = *fp;
	setvbuf(stderr, NULL, _IONBF, 0);

	return MadeConsole;
}
#endif
