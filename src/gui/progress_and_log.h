#ifndef _PROGRESSLOG_H
#define _PROGRESSLOG_H

#include <sys/time.h>
#include <gtk/gtk.h>

#define PROGRESS_NONE -2.0		// don't update the progress bar value
#define PROGRESS_PULSATE -1.0		// pulsate the progress bar
#define PROGRESS_RESET 0.0		// reset the progress bar
#define PROGRESS_DONE 1.0		// fill the progress bar
#define PROGRESS_TEXT_RESET ""		// reset the progress bar's text

void initialize_log_tags();

char* siril_log_message(const char* format, ...);
char* siril_log_color_message(const char* format, const char* color, ...);

void set_progress_bar_data(const char *text, double percent);

void show_time(struct timeval, struct timeval);
void show_time_msg(struct timeval t_start, struct timeval t_end, const char *msg);
void get_min_sec_from_timevals(struct timeval t_start, struct timeval t_end,
		int *min, int *sec);

void set_cursor_waiting(gboolean waiting);
void set_cursor(const gchar* cursor_name);

#endif
