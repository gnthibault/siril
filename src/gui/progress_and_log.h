#ifndef _PROGRESSLOG_H
#define _PROGRESSLOG_H

#include <sys/time.h>

#define PROGRESS_NONE -2.0		// don't update the progress bar value
#define PROGRESS_PULSATE -1.0		// pulsate the progress bar
#define PROGRESS_RESET 0.0		// reset the progress bar
#define PROGRESS_DONE 1.0		// fill the progress bar
#define PROGRESS_TEXT_RESET ""		// reset the progress bar's text

char* siril_log_message(const char* format, ...);
char* siril_log_color_message(const char* format, const char* color, ...);


void set_progress_bar_data(const char *text, double percent);
// deprecated progress bar functions below
void progress_bar_set_text(const char *text);
void progress_bar_reset_ready();
void progress_bar_set_percent(double percent);


void show_time(struct timeval, struct timeval);

#endif
