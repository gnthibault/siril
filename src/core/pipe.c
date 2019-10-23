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

/* This file manages the external command stream to siril, a named pipe */

#define PIPE_NAME_R "siril_command.in"
#define PIPE_NAME_W "siril_command.out"
#define PIPE_PATH_R "/tmp/" PIPE_NAME_R  // TODO: use g_get_tmp_dir()
#define PIPE_PATH_W "/tmp/" PIPE_NAME_W  // TODO: use g_get_tmp_dir()
#define PIPE_MSG_SZ 512	// max input command length

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <signal.h>
#ifdef _WIN32
// doc at: https://docs.microsoft.com/en-us/windows/desktop/ipc/named-pipes
// samples from: https://docs.microsoft.com/en-us/windows/desktop/ipc/using-pipes
// With windows, pipes can be bidirectional. To keep common code between windows and the
//	rest of the world, we still use unidirectional modes for windows named pipes.
// It's also different in the philosophy of pipe creation, where a pipe is created only
//	for one client and has to be created again to serve another.
#include <windows.h>
#include <conio.h>
#include <tchar.h>
#else
#include <sys/select.h>
#endif

#include "core/siril.h"
#include "pipe.h"
#include "command.h"
//#include "processing.h"
	void stop_processing_thread();	// avoid including everything
	gpointer waiting_for_thread();
	gboolean get_thread_run();
#include "gui/progress_and_log.h"

#ifdef _WIN32
LPTSTR lpszPipename_r = TEXT("\\\\.\\pipe\\" PIPE_NAME_R);
LPTSTR lpszPipename_w = TEXT("\\\\.\\pipe\\" PIPE_NAME_W);
HANDLE hPipe_r = INVALID_HANDLE_VALUE;
HANDLE hPipe_w = INVALID_HANDLE_VALUE;
#else
static int pipe_fd_r = -1;
static int pipe_fd_w = -1;
#endif

static GThread *pipe_thread_w, *worker_thread;
static int pipe_active;
static GCond write_cond, read_cond;
static GMutex write_mutex, read_mutex;
static GList *command_list, *pending_writes;
// ^ could use GQueue instead since it's used as a queue, avoids the cells memory leak

static void sigpipe_handler(int signum) { }	// do nothing

int pipe_create() {
#ifdef _WIN32
	if (hPipe_w != INVALID_HANDLE_VALUE || hPipe_r != INVALID_HANDLE_VALUE)
		return 0;

	hPipe_w = CreateNamedPipe(
			lpszPipename_w,           // pipe name
			PIPE_ACCESS_OUTBOUND,     // write access
			PIPE_TYPE_MESSAGE |       // message type pipe
			PIPE_READMODE_MESSAGE |   // message-read mode
			PIPE_WAIT,                // blocking mode
			PIPE_UNLIMITED_INSTANCES, // max. instances
			3*PIPE_MSG_SZ,            // output buffer size
			0,                        // input buffer size
			0,                        // client time-out
			NULL);                    // default security attribute
	if (hPipe_w == INVALID_HANDLE_VALUE)
	{
		siril_log_message(_("Output pipe creation failed with error %d\n"), GetLastError());
		return -1;
	}

	hPipe_r = CreateNamedPipe(
			lpszPipename_r,           // pipe name
			PIPE_ACCESS_INBOUND,      // read access
			PIPE_TYPE_MESSAGE |       // message type pipe
			PIPE_READMODE_MESSAGE |   // message-read mode
			PIPE_WAIT,                // blocking mode
			PIPE_UNLIMITED_INSTANCES, // max. instances
			3*PIPE_MSG_SZ,            // output buffer size
			0,                        // input buffer size
			0,                        // client time-out
			NULL);                    // default security attribute
	if (hPipe_r == INVALID_HANDLE_VALUE)
	{
		siril_log_message(_("Input pipe creation failed with error %d\n"), GetLastError());
		return -1;
	}
#else
	if (pipe_fd_r >= 0 || pipe_fd_w > 0) return 0;

	struct sigaction sa;
	sa.sa_handler = sigpipe_handler;
	sigemptyset(&sa.sa_mask);
	sa.sa_flags = SA_RESTART; /* Restart functions if
				     interrupted by handler */
	if (sigaction(SIGPIPE, &sa, NULL) == -1) {
		perror("sigaction");
		return -1;
	}

	struct stat st;
	if (stat(PIPE_PATH_R, &st)) {
		if (mkfifo(PIPE_PATH_R, 0666)) {
			siril_log_message(_("Could not create the named pipe "PIPE_PATH_R"\n"));
			perror("mkfifo");
			return -1;
		}
	}
	else if (!S_ISFIFO(st.st_mode)) {
		siril_log_message(_("The named pipe file " PIPE_PATH_R " already exists but is not a fifo, cannot create or open\n"));
		return -1;
	}

	if (stat(PIPE_PATH_W, &st)) {
		if (mkfifo(PIPE_PATH_W, 0666)) {
			siril_log_message(_("Could not create the named pipe "PIPE_PATH_W"\n"));
			perror("mkfifo");
			return -1;
		}
	}
	else if (!S_ISFIFO(st.st_mode)) {
		siril_log_message(_("The named pipe file " PIPE_PATH_W " already exists but is not a fifo, cannot create or open\n"));
		return -1;
	}
#endif
	return 0;
}

static int pipe_write(const char *string) {
#ifdef _WIN32
	int length;
	DWORD  retval ;
	if (hPipe_w == INVALID_HANDLE_VALUE)
		return -1;
	length = strlen(string);
	BOOL result = WriteFile(hPipe_w, string, length, &retval, NULL);
	if(result && retval == length)
		return 0;
	int err = GetLastError();
	if (err == ERROR_BROKEN_PIPE) {
		fprintf(stderr, "Output stream disconnected.\n");
		return 1;
	}
	else if (err == ERROR_NO_DATA) {
		fprintf(stderr, "Output stream closed on receiving side.\n");
		return 1;
	}
	else {
		fprintf(stderr, "Error writing to output stream; error code was 0x%08x.\n", err);
		return 1;
	}
#else
	int length, retval;
	if (pipe_fd_w <= 0)
		return -1;
	length = strlen(string);
	retval = write(pipe_fd_w, string, length);
	// buffer full, short writes and disconnections are treated as errors
	return retval != length;
#endif
}

int pipe_send_message(pipe_message msgtype, pipe_verb verb, const char *arg) {
#ifdef _WIN32
	if (hPipe_w == INVALID_HANDLE_VALUE) return -1;
#else
	if (pipe_fd_w <= 0) return -1;
#endif
	char *msg = NULL;

	switch (msgtype) {
		case PIPE_LOG:
			msg = malloc(strlen(arg) + 6);
			sprintf(msg, "log: %s", arg);
			break;
		case PIPE_STATUS:
			msg = malloc((arg ? strlen(arg) : 0) + 20);
			switch (verb) {
				case PIPE_STARTING:
					sprintf(msg, "status: starting %s", arg);
					break;
				case PIPE_SUCCESS:
					sprintf(msg, "status: success %s", arg);
					break;
				case PIPE_ERROR:
					sprintf(msg, "status: error %s", arg);
					break;
				case PIPE_EXIT:
					sprintf(msg, "status: exit\n");
					break;
				case PIPE_NA:
					free(msg);
					return -1;
			}
			break;
		case PIPE_PROGRESS:
			msg = strdup(arg);
			break;
	}

	if (msg) {
		g_mutex_lock(&write_mutex);
		pending_writes = g_list_append(pending_writes, msg);

		g_cond_signal(&write_cond);
		g_mutex_unlock(&write_mutex);
	}
	return 0;
}

int enqueue_command(char *command) {
	if (!strncmp(command, "cancel", 6))
		return 1;
	if ((command[0] >= 'a' && command[0] <= 'z') ||
			(command[0] >= 'A' && command[0] <= 'Z')) {
		g_mutex_lock(&read_mutex);
		command_list = g_list_append(command_list, command);
		g_cond_signal(&read_cond);
		g_mutex_unlock(&read_mutex);
	}
	return 0;
}

void empty_command_queue() {
	g_mutex_lock(&read_mutex);
	while (command_list) {
		free(command_list->data);
		command_list = g_list_next(command_list);
	}
	g_mutex_unlock(&read_mutex);
}

void *read_pipe(void *p) {
#ifdef _WIN32
	do {
		/* try to open the pipe */
		// will block until the other end is opened
		if (!ConnectNamedPipe(hPipe_r, NULL) && GetLastError() != ERROR_PIPE_CONNECTED) {
			siril_log_message(_("Could not open the named pipe\n"));
			break;
		}

		fprintf(stdout, "opened read pipe\n");
		/* now, try to read from it */
		int bufstart = 0;
		DWORD len;
		char buf[PIPE_MSG_SZ];
		do
		{
			// Read from the pipe.
			BOOL fSuccess = ReadFile(
					hPipe_r,         // pipe handle
					buf+bufstart,    // buffer to receive reply
					PIPE_MSG_SZ-1-bufstart, // size of buffer
					&len,            // number of bytes read
					NULL);           // not overlapped

			if ((fSuccess || GetLastError() == ERROR_MORE_DATA) && len > 0) {
				int i = 0, nbnl = 0;
				buf[len] = '\0';
				while (i < len && buf[i] != '\0') {
					if (buf[i] == '\n')
						nbnl++;
					i++;
				}
				if (nbnl == 0) {
					pipe_send_message(PIPE_STATUS, PIPE_ERROR, _("command too long or malformed\n"));
					fSuccess = FALSE;
				}

				if (fSuccess) {
					/* we have several commands in the buffer, we need to
					 * cut them, enqueue them and prepare next buffer for
					 * incomplete commands */
					char backup_char;
					char *command = NULL ;
					for (i = 0; i < len && buf[i] != '\0'; i++) {
						if (buf[i] == '\n') {
							backup_char = buf[i + 1];
							buf[i + 1] = '\0';
							command = strdup(buf+bufstart);
							buf[i + 1] = backup_char;
							bufstart = i + 1;

							if (enqueue_command(command)) {
								fSuccess = FALSE;
								break;
							}
						}
					}
					if (bufstart == i)
						bufstart = 0;
					else memcpy(buf, buf+bufstart, len-bufstart);
				}

			}
			if (!fSuccess && GetLastError() != ERROR_MORE_DATA) {
				fprintf(stdout, "closed read pipe\n");
				CloseHandle(hPipe_r);
				hPipe_r = INVALID_HANDLE_VALUE;
				empty_command_queue();
				if (get_thread_run()) {
					stop_processing_thread();
					pipe_send_message(PIPE_STATUS, PIPE_ERROR, _("command interrupted\n"));
				}
				break;
			}
		} while (1);
	} while (pipe_active);
#else
	do {
		// open will block until the other end is opened
		fprintf(stdout, "read pipe waiting to be opened...\n");
		if ((pipe_fd_r = open(PIPE_PATH_R, O_RDONLY)) == -1) {
			siril_log_message(_("Could not open the named pipe\n"));
			perror("open");
			break;
		}
		fprintf(stdout, "opened read pipe\n");

		int bufstart = 0;
		char buf[PIPE_MSG_SZ];
		do {
			int retval;
			fd_set rfds;
			FD_ZERO(&rfds);
			FD_SET(pipe_fd_r, &rfds);

			retval = select(pipe_fd_r+1, &rfds, NULL, NULL, NULL);
			if (retval == 1) {
				char *command;
				int len = read(pipe_fd_r, buf+bufstart, PIPE_MSG_SZ-1-bufstart);
				if (len == -1 || len == 0)
					retval = -1;
				else {
					int i = 0, nbnl = 0;
					buf[len] = '\0';
					while (i < len && buf[i] != '\0') {
						if (buf[i] == '\n')
							nbnl++;
						i++;
					}
					if (nbnl == 0) {
						pipe_send_message(PIPE_STATUS, PIPE_ERROR, _("command too long or malformed\n"));
						retval = -1;
					}

					if (retval == 1) {
						/* we have several commands in the buffer, we need to
						 * cut them, enqueue them and prepare next buffer for
						 * incomplete commands */
						char backup_char;
						for (i = 0; i < len && buf[i] != '\0'; i++) {
							if (buf[i] == '\n') {
								backup_char = buf[i + 1];
								buf[i + 1] = '\0';
								command = strdup(buf+bufstart);
								buf[i + 1] = backup_char;
								bufstart = i + 1;

								if (enqueue_command(command)) {
									retval = -1;
									free(command);
									break;
								}
							}
						}
						if (bufstart == i)
							bufstart = 0;
						else memcpy(buf, buf+bufstart, len-bufstart);
					}
				}
			}
			if (retval <= 0) {
				fprintf(stdout, "closed read pipe\n");
				close(pipe_fd_r);
				pipe_fd_r = -1;
				empty_command_queue();
				if (get_thread_run()) {
					stop_processing_thread();
					pipe_send_message(PIPE_STATUS, PIPE_ERROR, _("command interrupted\n"));
				}
				break;
			}
		} while (1);
	} while (pipe_active);
#endif

	return GINT_TO_POINTER(pipe_active ? -1 : 0);
}

void *process_commands(void *p) {
	while (pipe_active) {
		char *command;
		g_mutex_lock(&read_mutex);
		while (!command_list && pipe_active) {
			fprintf(stdout, "waiting for commands to be read from the pipe\n");
			g_cond_wait(&read_cond, &read_mutex);
		}
		if (!pipe_active) {
			g_mutex_unlock(&read_mutex);
			break;
		}

		command = (char*)command_list->data;
		command_list = g_list_next(command_list);
		g_mutex_unlock(&read_mutex);

		pipe_send_message(PIPE_STATUS, PIPE_STARTING, command);
		if (processcommand(command))
			pipe_send_message(PIPE_STATUS, PIPE_ERROR, command);
		else pipe_send_message(PIPE_STATUS, PIPE_SUCCESS, command);

		free(command);
		if (waiting_for_thread())
			empty_command_queue();
	}
	return NULL;
}

static void *write_pipe(void *p) {
	do {
		fprintf(stdout, "write pipe waiting to be opened...\n");
#ifdef _WIN32
		// will block until the other end is opened
		if (!ConnectNamedPipe(hPipe_w, NULL) && GetLastError() != ERROR_PIPE_CONNECTED) {
			siril_log_message(_("Could not open the named pipe\n"));
			break;
		}
#else
		// open will block until the other end is opened
		if ((pipe_fd_w = open(PIPE_PATH_W, O_WRONLY)) == -1) {
			siril_log_message(_("Could not open the named pipe\n"));
			perror("open");
			break;
		}
#endif
		fprintf(stdout, "opened write pipe\n");

		do {
			char *msg;
			// wait for messages to write
			g_mutex_lock(&write_mutex);
			while (!pending_writes && pipe_active)
				g_cond_wait(&write_cond, &write_mutex);
			if (!pipe_active) {
				g_mutex_unlock(&write_mutex);
				break;
			}

			msg = (char *)pending_writes->data;
			pending_writes = g_list_next(pending_writes);
			g_mutex_unlock(&write_mutex);

			if (pipe_write(msg)) {
#ifdef _WIN32
				CloseHandle(hPipe_w);
				hPipe_w = INVALID_HANDLE_VALUE;
#else
				fprintf(stdout, "closed write pipe\n");
				close(pipe_fd_w);
				pipe_fd_w = -1;
#endif
				free(msg);
				break;
			}
			free(msg);
		} while (1);
	} while (pipe_active);
	return GINT_TO_POINTER(-1);
}

/* not reentrant */
int pipe_start() {
	if (pipe_active)
		return 0;
	if (pipe_create())
		return -1;

	pipe_active = 1;
	worker_thread = g_thread_new("worker", process_commands, NULL);
	pipe_thread_w = g_thread_new("pipe writer", write_pipe, NULL);
	return 0;
}

/* not working, not used: blocked open calls are not signaled,
 * and blocked write_cond throws a deadlock error */
void pipe_stop() {
	fprintf(stdout, "closing pipes\n");
	g_mutex_lock(&read_mutex);
	g_mutex_lock(&write_mutex);
	pipe_active = 0;
#ifdef _WIN32
	if (hPipe_r != INVALID_HANDLE_VALUE)
		CloseHandle(hPipe_r);
	hPipe_r = INVALID_HANDLE_VALUE;
	if (hPipe_w != INVALID_HANDLE_VALUE)
		CloseHandle(hPipe_w);
	hPipe_w = INVALID_HANDLE_VALUE;
#else
	if (pipe_fd_r >= 0)
		close(pipe_fd_r);
	pipe_fd_r = -1;
	if (pipe_fd_w > 0)
		close(pipe_fd_w);
	pipe_fd_w = -1;
#endif
	g_cond_signal(&write_cond);
	g_cond_signal(&read_cond);
	g_mutex_unlock(&read_mutex);
	g_mutex_unlock(&write_mutex);
	if (pipe_thread_w)
		g_thread_join(pipe_thread_w);
	if (worker_thread)
		g_thread_join(worker_thread);
}
