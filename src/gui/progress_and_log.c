#include "progress_and_log.h"
#include <gtk/gtk.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include "callbacks.h"

/*
 * Progress bar static functions
 */

// should be static but is still used in preprocessing
void progress_bar_set_text(const char *text) {
	static GtkProgressBar *pbar = NULL;
	if (pbar == NULL)
		pbar = GTK_PROGRESS_BAR(
				gtk_builder_get_object(builder, "progressbar1"));
	/* It will not happen that text is NULL here, because it's
	 * catched by set_progress_bar_data() */
	if (!text || text[0] == '\0')
		text = "Ready.";
	gtk_progress_bar_set_text(pbar, text);
}

// should be static but is still used in preprocessing
void progress_bar_reset_ready() {
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
}

struct progress_bar_idle_data {
	char *progress_bar_text;
	double progress_bar_percent;
};

// should be static but is still used in preprocessing
/* http://developer.gnome.org/gtk3/3.4/GtkProgressBar.html */
void progress_bar_set_percent(double percent) {
	static GtkProgressBar *pbar = NULL;
	if (pbar == NULL)
		pbar = GTK_PROGRESS_BAR(
				gtk_builder_get_object(builder, "progressbar1"));
	if (percent == PROGRESS_PULSATE) {
#ifdef _OPENMP
#pragma omp critical
#endif
		gtk_progress_bar_pulse(pbar);
	}
	else {
		assert(percent >= 0.0 && percent <= 1.0);
		gtk_progress_bar_set_fraction(pbar, percent);
	}
}

static gboolean progress_bar_idle_callback(gpointer p) {
	struct progress_bar_idle_data *data = (struct progress_bar_idle_data *) p;

	if (data->progress_bar_text) {
		progress_bar_set_text(data->progress_bar_text);
		free(data->progress_bar_text);
	}
	if (data->progress_bar_percent != PROGRESS_NONE)
		progress_bar_set_percent(data->progress_bar_percent);
	free(data);
	return FALSE;	// only run once
}

// Thread-safe progress bar update.
// text can be NULL, percent can be -1 for pulsating, -2 for nothing, or between 0 and 1 for percent
void set_progress_bar_data(const char *text, double percent) {
	struct progress_bar_idle_data *data;
	g_mutex_lock(&com.mutex);
	//fprintf(stdout, "progress: %s, %g\n", text ? text : "NULL", percent);
	data = malloc(sizeof(struct progress_bar_idle_data));
	data->progress_bar_text = text ? strdup(text) : NULL;
	data->progress_bar_percent = percent;
	assert(percent == PROGRESS_PULSATE || percent == PROGRESS_NONE ||
			(percent >= 0.0 && percent <= 1.0));
	gdk_threads_add_idle(progress_bar_idle_callback, data);
	g_mutex_unlock(&com.mutex);
}

/*
 * Siril log message static functions
 */

struct log_message {
	char *timestamp;
	char *message;
	const char* color;
};

// The main thread internal function that does the printing.
static gboolean idle_messaging(gpointer p) {
	static GtkTextBuffer *tbuf = NULL;
	static GtkTextView *text = NULL;
	GtkTextIter iter;
	struct log_message *log = (struct log_message *) p;

	if (!tbuf) {
		text = GTK_TEXT_VIEW(gtk_builder_get_object(builder, "output"));
		tbuf = gtk_text_view_get_buffer(text);
	}

	if (log->message[0] == '\n' && log->message[1] == '\0') {
		gtk_text_buffer_get_start_iter(tbuf, &iter);
		gtk_text_buffer_insert(tbuf, &iter, log->message, strlen(log->message));
		free(log);
		return FALSE;
	}

	gtk_text_buffer_get_end_iter(tbuf, &iter);
	gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter, log->timestamp,
			strlen(log->timestamp), "bold", NULL);

	if (log->color == NULL)
		gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter, log->message,
				strlen(log->message), "normal", NULL);
	else
		gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter, log->message,
				strlen(log->message), log->color, NULL);

	/* scroll to end */
	gtk_text_buffer_get_end_iter(tbuf, &iter);
	GtkTextMark *insert_mark = gtk_text_buffer_get_insert(tbuf);
	gtk_text_buffer_place_cursor(tbuf, &iter);
	gtk_text_view_scroll_to_mark(text, insert_mark, 0.0, TRUE, 0.0, 1.0);
	gtk_widget_queue_draw(GTK_WIDGET(text));

	free(log->timestamp);
	free(log->message);
	free(log);
	return FALSE;
}

/* This function writes a message on Siril's console/log. It is not thread safe.
 * There is a limit in number of characters that it is able to write in one call: 1023.
 * Return value is the string printed from arguments, or NULL if argument was empty or
 * only newline. It is an allocated string and must not be freed. It can be
 * reused until next call to this function.
 */
static char* siril_log_internal(const char* format, const char* color, va_list arglist) {
	static char *msg = NULL;
	struct tm *now;
	time_t now_sec;
	char timestamp[30];
	struct log_message *new_msg;

	if (msg == NULL) {
		msg = malloc(1024);
		msg[1023] = '\0';
	}

	vsnprintf(msg, 1023, format, arglist);

	if (msg == NULL || msg[0] == '\0')
		return NULL;

	if (msg[0] == '\n' && msg[1] == '\0') {
		fputc('\n', stdout);
		new_msg = malloc(sizeof(struct log_message));
		new_msg->timestamp = NULL;
		new_msg->message = "\n";
		new_msg->color = NULL;
		gdk_threads_add_idle(idle_messaging, new_msg);
		return NULL;
	}

	fprintf(stdout, "log: %s", msg);
	now_sec = time(NULL);
	now = localtime(&now_sec);
	g_snprintf(timestamp, sizeof(timestamp), "%.2d:%.2d:%.2d: ", now->tm_hour,
			now->tm_min, now->tm_sec);

	new_msg = malloc(sizeof(struct log_message));
	new_msg->timestamp = strdup(timestamp);
	new_msg->message = strdup(msg);
	new_msg->color = color;
	gdk_threads_add_idle(idle_messaging, new_msg);

	return msg;
}

char* siril_log_message(const char* format, ...) {
	va_list args;
	va_start(args, format);
	g_mutex_lock(&com.mutex);
	char *msg = siril_log_internal(format, NULL, args);
	g_mutex_unlock(&com.mutex);
	va_end(args);
	return msg;
}

char* siril_log_color_message(const char* format, const char* color, ...) {
	va_list args;
	va_start(args, color);
	g_mutex_lock(&com.mutex);
	char *msg = siril_log_internal(format, color, args);
	g_mutex_unlock(&com.mutex);
	va_end(args);
	return msg;
}

void show_time(struct timeval t_start, struct timeval t_end) {
	float time = (float) (((t_end.tv_sec - t_start.tv_sec) * 1000000L
			+ t_end.tv_usec) - t_start.tv_usec) / 1000000.0;

	if (time > 60.f) {
		int min = (int) time / 60;
		int sec = (int) time % 60 + 1;
		siril_log_color_message(_("Execution time: %d min %.2d s.\n"), "green",
				min, sec);
	} else if (time < 1.f) {
		float ms = time * 1000.f;
		siril_log_color_message(_("Execution time: %.2f ms.\n"), "green", ms);
	} else {
		siril_log_color_message(_("Execution time: %.2f s.\n"), "green", time);
	}
}

