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

/* This file manages the external command stream to siril, a named pipe */

#define PIPE_NAME_R "siril_command.in"
#define PIPE_NAME_W "siril_command.out"
#define PIPE_PATH_R "/tmp/" PIPE_NAME_R  // TODO: use g_get_tmp_dir()
#define PIPE_PATH_W "/tmp/" PIPE_NAME_W  // TODO: use g_get_tmp_dir()
#define PIPE_MSG_SZ 512

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <signal.h>
#include <sys/select.h>

#include "core/siril.h"
#include "pipe.h"
#include "command.h"
//#include "processing.h"
	void stop_processing_thread();	// avoid including everything
	gpointer waiting_for_thread();
#include "gui/progress_and_log.h"

/* TODO:
 * use PIPE_STATUS
 * check the mutexes around the conditional variables and fix the deadlock on stop
 * windows version
 */

#ifdef _WIN32
// https://docs.microsoft.com/en-us/windows/desktop/ipc/named-pipes
#else
static int pipe_fd_r = -1;
static int pipe_fd_w = -1;
#endif

static char pipe_buf_o[PIPE_MSG_SZ];
static GThread *pipe_thread_w, *worker_thread;
static int pipe_active;
static GCond write_cond, read_cond;
static GMutex write_mutex, read_mutex;
static GList *command_list;

static void sigpipe_handler(int signum) { }	// do nothing

int pipe_create() {
	/* using commands from a thread that is not the GUI thread is only
	 * supported when this is toggled, otherwise changing the cursor and
	 * possibly other graphical operations are badly done. */
	//com.script = TRUE;
#ifdef _WIN32
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

/* not reentrant, arg needs the newline for logs but not for others */
int pipe_send_message(pipe_message msgtype, pipe_verb verb, const char *arg) {
	if (pipe_fd_w <= 0) return -1;

	g_mutex_lock(&write_mutex);
	switch (msgtype) {
		case PIPE_LOG:
			snprintf(pipe_buf_o, PIPE_MSG_SZ, "log: %s", arg);
			break;
		case PIPE_STATUS:
			switch (verb) {
				case PIPE_STARTING:
					snprintf(pipe_buf_o, PIPE_MSG_SZ, "status: starting %s\n", arg);
					break;
				case PIPE_SUCCESS:
					snprintf(pipe_buf_o, PIPE_MSG_SZ, "status: success %s\n", arg);
					break;
				case PIPE_ERROR:
					snprintf(pipe_buf_o, PIPE_MSG_SZ, "status: error %s\n", arg);
					break;
				case PIPE_EXIT:
					snprintf(pipe_buf_o, PIPE_MSG_SZ, "status: exit\n");
					break;
				case PIPE_NA:
					return -1;
			}
			break;
		case PIPE_PROGRESS:
			snprintf(pipe_buf_o, PIPE_MSG_SZ, "%s", arg);
			break;

	}

	g_cond_signal(&write_cond);
	g_mutex_unlock(&write_mutex);
	//return pipe_write(pipe_buf_o);
	return 0;
}

#ifdef _WIN32
#else
void *read_pipe_old(void *p) {
	char pipe_buf_i[PIPE_MSG_SZ];

	do {
		// open will block until the other end is opened
		fprintf(stdout, "read pipe waiting to be opened...\n");
		if ((pipe_fd_r = open(PIPE_PATH_R, O_RDONLY)) == -1) {
			siril_log_message(_("Could not open the named pipe\n"));
			perror("open");
			break;
		}
		fprintf(stdout, "opened read pipe\n");

		FILE *fd = fdopen(pipe_fd_r, "r");
		if (!fd) { // should never happen
			perror("fdopen");
			break;
		}

		/*do {
			if (!fgets(pipe_buf_i, PIPE_MSG_SZ, fd)) {
				fprintf(stdout, "closed read pipe\n");
				fclose(fd);
				pipe_fd_r = -1;
				break;
			}

			processcommand(pipe_buf_i);

		} while (1);*/
		// even easier and safer:
		execute_script(fd);
	} while (pipe_active);

	return GINT_TO_POINTER(pipe_active ? -1 : 0);
}

void *read_pipe(void *p) {
	do {
		// open will block until the other end is opened
		fprintf(stdout, "read pipe waiting to be opened...\n");
		if ((pipe_fd_r = open(PIPE_PATH_R, O_RDONLY)) == -1) {
			siril_log_message(_("Could not open the named pipe\n"));
			perror("open");
			break;
		}
		fprintf(stdout, "opened read pipe\n");

		FILE *fd = fdopen(pipe_fd_r, "r");
		if (!fd) { // should never happen
			perror("fdopen");
			break;
		}

		do {
			int retval;
			fd_set rfds;
			FD_ZERO(&rfds);
			FD_SET(pipe_fd_r, &rfds);

			retval = select(pipe_fd_r+1, &rfds, NULL, NULL, NULL);
			if (retval == 1) {
				char *command = malloc(PIPE_MSG_SZ);
				if (!fgets(command, PIPE_MSG_SZ, fd)) {
					retval = -1;
				}

				if (retval == 1) {
					if (!strncmp(command, "cancel", 6))
						retval = -1;
					else if (command[0] > 'a' && command[0] < 'z') {
						g_mutex_lock(&read_mutex);
						command_list = g_list_append(command_list,
								command);
						g_cond_signal(&read_cond);
						g_mutex_unlock(&read_mutex);
					}
				}
			}
			if (retval <= 0) {
				fprintf(stdout, "closed read pipe\n");
				fclose(fd);
				pipe_fd_r = -1;
				stop_processing_thread();
				break;
			}
		} while (1);
	} while (pipe_active);

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

		processcommand(command);

		free(command);
		if (waiting_for_thread()) {
			g_mutex_lock(&read_mutex);
			command_list = NULL;
			g_mutex_unlock(&read_mutex);
		}
	}
	return NULL;
}

static void *write_pipe(void *p) {
	do {
		// open will block until the other end is opened
		fprintf(stdout, "write pipe waiting to be opened...\n");
		if ((pipe_fd_w = open(PIPE_PATH_W, O_WRONLY)) == -1) {
			siril_log_message(_("Could not open the named pipe\n"));
			perror("open");
			break;
		}
		fprintf(stdout, "opened write pipe\n");

		do {
			// wait for messages to write
			g_mutex_lock(&write_mutex);
			while (pipe_buf_o[0] == '\0')
				g_cond_wait(&write_cond, &write_mutex);
			g_mutex_unlock(&write_mutex);

			if (pipe_write(pipe_buf_o)) {
				fprintf(stdout, "closed write pipe\n");
				close(pipe_fd_w);
				pipe_fd_w = -1;
				break;
			}
			pipe_buf_o[0] = '\0';
		} while (1);
	} while (pipe_active);
	return GINT_TO_POINTER(-1);
}

#endif

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

void pipe_stop() {
	fprintf(stdout, "closing pipes\n");
	pipe_active = 0;
	if (pipe_fd_r >= 0)
		close(pipe_fd_r);
	if (pipe_fd_w > 0)
		close(pipe_fd_w);
	pipe_buf_o[0] = '\1';
	g_cond_signal(&write_cond);
	g_cond_signal(&read_cond);
	if (pipe_thread_w)
		g_thread_join(pipe_thread_w);
	if (worker_thread)
		g_thread_join(worker_thread);
}
