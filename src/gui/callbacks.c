/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2016 team free-astro (see more in AUTHORS file)
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

#include <gtk/gtk.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/time.h>
#include <dirent.h>
#include <stdarg.h>
#include <math.h>	// for M_PI
#include <libgen.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/initfile.h"
#include "core/undo.h"
#include "gui/callbacks.h"
#include "gui/PSF_list.h"
#include "gui/histogram.h"
#include "algos/colors.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/fft.h"
#include "algos/Def_Wavelet.h"
#include "algos/cosmetic_correction.h"
#include "algos/gradient.h"
#include "io/conversion.h"
#include "io/films.h"
#include "io/ser.h"
#include "io/single_image.h"
#include "core/command.h"	// for processcommand()
#include "registration/registration.h"
#include "stacking/stacking.h"
#include "compositing/compositing.h"
#include "compositing/align_rgb.h"
#include "plot.h"
#ifdef HAVE_OPENCV
#include "opencv/opencv.h"
#endif

static enum {
	CD_NULL, CD_INCALL, CD_EXCALL, CD_QUIT
} confirm;
static gboolean is_shift_on = FALSE;

layer_info predefined_layers_colors[] = {
		/* name, lambda, lo, hi, c/over, c/under, mode */
		{ "Luminance", 0.0, 0, 0, FALSE, FALSE, NORMAL_DISPLAY }, // no color, undefined value is <0
		{ "Red", 650.0, 0, 0, FALSE, FALSE, NORMAL_DISPLAY }, // approx. of the middle of the color
		{ "Green", 530.0, 0, 0, FALSE, FALSE, NORMAL_DISPLAY },	// approx. of the middle of the color
		{ "Blue", 450.0, 0, 0, FALSE, FALSE, NORMAL_DISPLAY }// approx. of the middle of the color
};

/* remap index data, an index for each layer */
static BYTE *remap_index[MAXGRAYVPORT];
static float last_pente[MAXGRAYVPORT];
static display_mode last_mode[MAXGRAYVPORT];

/*****************************************************************************
 *                    S T A T I C      F U N C T I O N S                     *
 ****************************************************************************/

/* this function calculates the "fit to window" zoom values, given the window
 * size in argument and the image size in gfit.
 * Should not be called before displaying the main gray window when using zoom to fit */
static double get_zoom_val() {
	int window_width, window_height;
	double wtmp, htmp;
	static GtkWidget *scrolledwin = NULL;
	if (scrolledwin == NULL)
		scrolledwin = lookup_widget("scrolledwindowr");
	if (com.zoom_value > 0.)
		return com.zoom_value;
	/* else if zoom is < 0, it means fit to window */
	window_width = gtk_widget_get_allocated_width(scrolledwin);
	window_height = gtk_widget_get_allocated_height(scrolledwin);
	if (gfit.rx == 0 || gfit.ry == 0 || window_height <= 1 || window_width <= 1)
		return 1.0;
	wtmp = (double) window_width / (double) gfit.rx;
	htmp = (double) window_height / (double) gfit.ry;
	//fprintf(stdout, "computed fit to window zooms: %f, %f\n", wtmp, htmp);
	return min(wtmp, htmp);
}

/*
 * Progress bar static functions
 */

static void progress_bar_set_text(const char *text) {
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

static void progress_bar_reset_ready() {
	set_progress_bar_data(PROGRESS_TEXT_RESET, PROGRESS_RESET);
}

struct progress_bar_idle_data {
	char *progress_bar_text;
	double progress_bar_percent;
};

/* http://developer.gnome.org/gtk3/3.4/GtkProgressBar.html */
static void progress_bar_set_percent(double percent) {
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

/*
 * Returns the modifier mask. For Linux it is Control key, but for Apple - OS X
 * it should be Apple Key -> Mod2
 */

static GdkModifierType get_default_modifier() {
	GdkDisplay *display = gdk_display_get_default();
	GdkKeymap *keymap = gdk_keymap_get_for_display(display);
	GdkModifierType primary, real;

	g_return_val_if_fail(GDK_IS_KEYMAP (keymap), 0);

	/* Retrieve the real modifier mask */
	real = gdk_keymap_get_modifier_mask(keymap,
			GDK_MODIFIER_INTENT_PRIMARY_ACCELERATOR);

	primary = real;

	/* We need to translate the real modifiers into a virtual modifier
	 (like Super, Meta, etc.).
	 The following call adds the virtual modifiers for each real modifier
	 defined in primary.
	 */
	gdk_keymap_add_virtual_modifiers(keymap, &primary);

	if (primary != real) {
		/* In case the virtual and real modifiers are different, we need to
		 remove the real modifier from the result, and keep only the
		 virtual one.
		 */
		primary &= ~real;
	}
	return primary;
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

/*
 * Dialog window text static functions
 */

struct _dialog_data {
	const char *text;
	const char *title;
	const char *icon;
};

static gboolean show_dialog_idle(gpointer p) {
	static GtkLabel *label = NULL;
	struct _dialog_data *args = (struct _dialog_data *) p;
	GtkImage *image = GTK_IMAGE(lookup_widget("image1"));
	if (label == NULL) {
		label = GTK_LABEL(gtk_builder_get_object(builder, "labeldialog1"));
	}
	gtk_window_set_title(GTK_WINDOW(lookup_widget("dialog1")), args->title);
	gtk_image_set_from_icon_name(image, args->icon, GTK_ICON_SIZE_DIALOG);
	gtk_label_set_text(label, args->text);
	gtk_widget_show(lookup_widget("dialog1"));
	gtk_window_present (GTK_WINDOW(lookup_widget("dialog1")));
	free(args);
	return FALSE;
}

/*
 * Wavelet static functions
 */

static void reset_scale_w() {
	static GtkSpinButton *spin_w[6] = { NULL, NULL, NULL, NULL, NULL, NULL };
	int i;

	if (spin_w[0] == NULL) {
		spin_w[0] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w0"));
		spin_w[1] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w1"));
		spin_w[2] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w2"));
		spin_w[3] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w3"));
		spin_w[4] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w4"));
		spin_w[5] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w5"));
	}

	for (i = 0; i < 6; i++) {
		g_signal_handlers_block_by_func(spin_w[i], on_spin_w_changed, NULL);
		gtk_spin_button_set_value(spin_w[i], 1.f);
		g_signal_handlers_unblock_by_func(spin_w[i], on_spin_w_changed, NULL);
	}

	gtk_widget_set_sensitive(lookup_widget("button_apply_w"), FALSE);
}

static void update_wavelets() {
	float scale[6];
	static GtkSpinButton *spin_w[6] = { NULL, NULL, NULL, NULL, NULL, NULL };
	int i;
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];
	const char *tmpdir;

	tmpdir = g_get_tmp_dir();

	if (spin_w[0] == NULL) {
		spin_w[0] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w0"));
		spin_w[1] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w1"));
		spin_w[2] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w2"));
		spin_w[3] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w3"));
		spin_w[4] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w4"));
		spin_w[5] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w5"));
	}

	for (i = 0; i < 6; i++)
		scale[i] = (float) gtk_spin_button_get_value(spin_w[i]);

	set_cursor_waiting(TRUE);

	for (i = 0; i < gfit.naxes[2]; i++) {
		dir[i] = malloc(strlen(tmpdir) + strlen(File_Name_Transform[i]) + 2);
		strcpy(dir[i], tmpdir);
		strcat(dir[i], "/");
		strcat(dir[i], File_Name_Transform[i]);
		wavelet_reconstruct_file(dir[i], scale, gfit.pdata[i]);
		free(dir[i]);
	}

	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

/*
 * Memory label static functions
 */

struct _label_data {
	const char *label_name;
	char *text;
};

static gboolean set_label_text_idle(gpointer p) {
	struct _label_data *args = (struct _label_data *) p;
	GtkLabel *label = GTK_LABEL(
			gtk_builder_get_object(builder, args->label_name));
	gtk_label_set_text(label, args->text);
	free(args->text);
	free(args);
	return FALSE;
}

static void set_label_text_from_main_thread(const char *label_name, const char *text) {
	struct _label_data *data = malloc(sizeof(struct _label_data));
	data->label_name = label_name;
	data->text = strdup(text);
	gdk_threads_add_idle(set_label_text_idle, data);
}

/*
 * Image display static functions
 */

static gboolean inimage(GdkEvent *event) {
	double zoom = get_zoom_val();
	return ((GdkEventButton*) event)->x > 0
			&& ((GdkEventButton*) event)->x < gfit.rx * zoom
			&& ((GdkEventButton*) event)->y > 0
			&& ((GdkEventButton*) event)->y < gfit.ry * zoom;
}

static void draw_empty_image(cairo_t *cr, guint width, guint height) {
	// black image with red square
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_rectangle(cr, 0, 0, width, height);
	cairo_fill(cr);
	cairo_set_source_rgb(cr, 0.3, 0, 0);
	cairo_rectangle(cr, 100, 70, 50, 50);
	cairo_fill(cr);
}


static void remaprgb(void) {
	guchar *dst;
	guchar *bufr, *bufg, *bufb;
	gint i, j;
	int nbdata;

	fprintf(stderr, "remaprgb\n");
	if (!isrgb(&gfit))
		return;

	// allocate if not already done or the same size
	if (cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, gfit.rx)
			!= com.surface_stride[RGB_VPORT]
			|| gfit.ry != com.surface_height[RGB_VPORT]
			|| !com.surface[RGB_VPORT] || !com.rgbbuf) {
		guchar *oldbuf = com.rgbbuf;
		fprintf(stderr, "RGB display buffers and surface (re-)allocation\n");
		com.surface_stride[RGB_VPORT] = cairo_format_stride_for_width(
				CAIRO_FORMAT_RGB24, gfit.rx);
		com.surface_height[RGB_VPORT] = gfit.ry;
		com.rgbbuf = realloc(com.rgbbuf,
				com.surface_stride[RGB_VPORT] * gfit.ry * sizeof(guchar));
		if (com.rgbbuf == NULL) {
			fprintf(stderr,
					"Could not allocate memory for RGB buffer (out of memory?)\n");
			if (oldbuf)
				free(oldbuf);
			return;
		}
		if (com.surface[RGB_VPORT])
			cairo_surface_destroy(com.surface[RGB_VPORT]);
		com.surface[RGB_VPORT] = cairo_image_surface_create_for_data(com.rgbbuf,
				CAIRO_FORMAT_RGB24, gfit.rx, gfit.ry,
				com.surface_stride[RGB_VPORT]);
		if (cairo_surface_status(com.surface[RGB_VPORT])
				!= CAIRO_STATUS_SUCCESS) {
			fprintf(stderr,
					"Error creating the Cairo image surface for the RGB image\n");
			cairo_surface_destroy(com.surface[RGB_VPORT]);
			com.surface[RGB_VPORT] = NULL;
			return;
		}
	}
	// WARNING : this assumes that R, G and B buffers are already allocated and mapped
	// it seems ok, but one can probably imagine situations where it segfaults
	bufr = com.graybuf[RED_VPORT];
	bufg = com.graybuf[GREEN_VPORT];
	bufb = com.graybuf[BLUE_VPORT];
	if (bufr == NULL || bufg == NULL || bufb == NULL) {
		fprintf(stderr, "remaprgb: gray buffers not allocated for display\n");
		return;
	}
	dst = com.rgbbuf;	// index is j
	nbdata = gfit.rx * gfit.ry * 4;	// source images are 32-bit RGBA

	for (i = 0, j = 0; i < nbdata; i += 4) {
		dst[j++] = bufb[i];
		dst[j++] = bufg[i];
		dst[j++] = bufr[i];
		j++;		// alpha padding
	}

	// flush to ensure all writing to the image was done and redraw the surface
	cairo_surface_flush(com.surface[RGB_VPORT]);
	cairo_surface_mark_dirty(com.surface[RGB_VPORT]);
}

static void set_viewer_mode_widgets_sensitive(gboolean sensitive) {
	static GtkWidget *scalemax = NULL;
	static GtkWidget *scalemin = NULL;
	static GtkWidget *entrymin = NULL;
	static GtkWidget *entrymax = NULL;
	static GtkWidget *minmax = NULL;
	static GtkWidget *hilo = NULL;
	static GtkWidget *user = NULL;

	if (!scalemax) {
		scalemax = lookup_widget("scalemax");
		scalemin = lookup_widget("scalemin");
		entrymin = lookup_widget("min_entry");
		entrymax = lookup_widget("max_entry");
		minmax = lookup_widget("radiobutton_minmax");
		hilo = lookup_widget("radiobutton_hilo");
		user = lookup_widget("radiobutton_user");
	}
	gtk_widget_set_sensitive(scalemax, sensitive);
	gtk_widget_set_sensitive(scalemin, sensitive);
	gtk_widget_set_sensitive(entrymin, sensitive);
	gtk_widget_set_sensitive(entrymax, sensitive);
	gtk_widget_set_sensitive(minmax, sensitive);
	gtk_widget_set_sensitive(hilo, sensitive);
	gtk_widget_set_sensitive(user, sensitive);
}

static int make_index_for_current_display(display_mode mode, WORD lo, WORD hi,
		int vport);
static int make_index_for_rainbow(BYTE index[][3]);

/* enables or disables the "display reference" checkbox in registration preview */
static void enable_view_reference_checkbox(gboolean status) {
	static GtkToggleButton *check_display_ref = NULL;
	static GtkWidget *widget = NULL, *labelRegRef = NULL;
	if (check_display_ref == NULL) {
		check_display_ref = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "checkbutton_displayref"));
		widget = GTK_WIDGET(check_display_ref);
		labelRegRef = lookup_widget("labelRegRef");
	}
	if (status && gtk_widget_get_sensitive(widget))
		return;	// may be already enabled but deactivated by user, don't force it again
	gtk_widget_set_sensitive(widget, status);
	gtk_widget_set_visible(labelRegRef, !status);
	gtk_toggle_button_set_active(check_display_ref, status);
}

/* vport can be -1 if the correct viewport should be tested */
static void test_and_allocate_reference_image(int vport) {
	static GtkComboBox *cbbt_layers = NULL;
	if (cbbt_layers == NULL) {
		cbbt_layers = GTK_COMBO_BOX(
				gtk_builder_get_object(builder, "comboboxreglayer"));
	}
	if (vport == -1)
		vport = gtk_combo_box_get_active(cbbt_layers);

	if (sequence_is_loaded() && com.seq.current == com.seq.reference_image
			&& gtk_combo_box_get_active(cbbt_layers) == vport) {
		/* this is the registration layer and the reference frame,
		 * save the buffer for alignment preview */
		if (!com.refimage_regbuffer || !com.refimage_surface) {
			guchar *oldbuf = com.refimage_regbuffer;
			com.refimage_regbuffer = realloc(com.refimage_regbuffer,
					com.surface_stride[vport] * gfit.ry * sizeof(guchar));
			if (com.refimage_regbuffer == NULL) {
				fprintf(stderr,
						"Could not allocate memory for the reference image buffer\n");
				if (oldbuf)
					free(oldbuf);
				return;
			}

			if (com.refimage_surface)
				cairo_surface_destroy(com.refimage_surface);
			com.refimage_surface = cairo_image_surface_create_for_data(
					com.refimage_regbuffer, CAIRO_FORMAT_RGB24, gfit.rx,
					gfit.ry, com.surface_stride[vport]);
			if (cairo_surface_status(com.refimage_surface)
					!= CAIRO_STATUS_SUCCESS) {
				fprintf(stderr,
						"Error creating the Cairo image surface for the reference image.\n");
				cairo_surface_destroy(com.refimage_surface);
				com.refimage_surface = NULL;
			} else {
				fprintf(stdout,
						"Saved the reference frame buffer for alignment preview.\n");
				enable_view_reference_checkbox(TRUE);
			}
		}
		memcpy(com.refimage_regbuffer, com.graybuf[vport],
				com.surface_stride[vport] * gfit.ry * sizeof(guchar));
		cairo_surface_flush(com.refimage_surface);
		cairo_surface_mark_dirty(com.refimage_surface);
	}
}

static void remap(int vport) {
	// This function maps fit data with a linear LUT between lo and hi levels
	// to the buffer to be displayed; display only is modified
	guint x, y;
	BYTE *dst, *index, rainbow_index[UCHAR_MAX + 1][3];
	WORD *src, hi, lo;
	display_mode mode;
	color_map color;
	gboolean do_cut_over, inverted;

	fprintf(stderr, "remap %d\n", vport);
	if (vport == RGB_VPORT) {
		remaprgb();
		return;
	}

	int no_data = 0;
	if (single_image_is_loaded()) {
	       if (vport >= com.uniq->nb_layers)
		       no_data = 1;
	}
	else if (sequence_is_loaded()) {
		if (vport >= com.seq.nb_layers)
			no_data = 1;
	}
	else no_data = 1;
	if (no_data) {
		fprintf(stderr, "vport is out of bounds or data is not loaded yet\n");
		return;
	}

	// allocate if not already done or the same size
	if (cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, gfit.rx) !=
			com.surface_stride[vport] ||
			gfit.ry != com.surface_height[vport] ||
			!com.surface[vport] ||
			!com.graybuf[vport]) {
		guchar *oldbuf = com.graybuf[vport];
		fprintf(stderr, "Gray display buffers and surface (re-)allocation\n");
		if (gfit.rx == 0 || gfit.ry == 0) {
			fprintf(stderr, "gfit has a zero size, must not happen!\n");
			return;
		}
		com.surface_stride[vport] = cairo_format_stride_for_width(
				CAIRO_FORMAT_RGB24, gfit.rx);
		com.surface_height[vport] = gfit.ry;
		com.graybuf[vport] = realloc(com.graybuf[vport],
				com.surface_stride[vport] * gfit.ry * sizeof(guchar));
		if (com.graybuf[vport] == NULL) {
			fprintf(stderr,
					"Could not allocate memory for gray buffer %d (out of memory?)\n",
					vport);
			if (oldbuf)
				free(oldbuf);
			return;
		}
		if (com.surface[vport])
			cairo_surface_destroy(com.surface[vport]);
		com.surface[vport] = cairo_image_surface_create_for_data(
				com.graybuf[vport], CAIRO_FORMAT_RGB24, gfit.rx, gfit.ry,
				com.surface_stride[vport]);
		if (cairo_surface_status(com.surface[vport]) != CAIRO_STATUS_SUCCESS) {
			fprintf(stderr,
					"Error creating the Cairo image surface for vport %d\n",
					vport);
			cairo_surface_destroy(com.surface[vport]);
			com.surface[vport] = NULL;
			return;
		}
	}
	if (single_image_is_loaded() && com.seq.current != RESULT_IMAGE) {
		mode = com.uniq->layers[vport].rendering_mode;
		hi = com.uniq->layers[vport].hi;
		lo = com.uniq->layers[vport].lo;
		do_cut_over = com.uniq->layers[vport].cut_over;
	} else if (sequence_is_loaded() && vport < com.seq.nb_layers) {
		// the check above is needed because there may be a different
		// number of channels between the unique image and the sequence
		mode = com.seq.layers[vport].rendering_mode;
		hi = com.seq.layers[vport].hi;
		lo = com.seq.layers[vport].lo;
		do_cut_over = com.seq.layers[vport].cut_over;
	} else {
		fprintf(stderr, "BUG in unique image remap\n");
		return;
	}

	if (lo > hi) {
		// negative display
		WORD tmp = hi;
		hi = lo;
		lo = tmp;
		inverted = TRUE;
	} else
		inverted = FALSE;

	if (mode == HISTEQ_DISPLAY) {
		double hist_sum;
		double nb_pixels;
		size_t hist_nb_bins;
		size_t i;
		gsl_histogram *histo = NULL;

		compute_histo_for_gfit(1);
		histo = com.layers_hist[vport];
		hist_nb_bins = gsl_histogram_bins(histo);
		/*if (hist_nb_bins <= USHRT_MAX) {
		 fprintf(stderr, "Error remapping: histogram is not the correct size\n");
		 return;
		 }*/
		nb_pixels = (double) (gfit.rx * gfit.ry);
		// build the remap_index
		if (!remap_index[vport])
			remap_index[vport] = malloc(USHRT_MAX + 1);

		remap_index[vport][0] = 0;
		hist_sum = gsl_histogram_get(histo, 0);
		for (i = 1; i < hist_nb_bins; i++) {
			hist_sum += gsl_histogram_get(histo, i);
			remap_index[vport][i] = round_to_BYTE(
					(hist_sum / nb_pixels) * UCHAR_MAX_DOUBLE);
		}

		last_mode[vport] = mode;
		set_viewer_mode_widgets_sensitive(FALSE);
	} else {
		// for all other modes, the index can be reused
		make_index_for_current_display(mode, lo, hi, vport);
		if (mode == STF_DISPLAY)
			set_viewer_mode_widgets_sensitive(FALSE);
		else
			set_viewer_mode_widgets_sensitive(TRUE);
	}

	src = gfit.pdata[vport];
	/* Siril's FITS are stored bottom to top, so mapping needs to revert data order */
	dst = com.graybuf[vport];

	color = gtk_toggle_tool_button_get_active(
			GTK_TOGGLE_TOOL_BUTTON(lookup_widget("colormap_button")));

	if (color == RAINBOW_COLOR)
		make_index_for_rainbow(rainbow_index);
	index = remap_index[vport];

#ifdef _OPENMP
#pragma omp parallel for num_threads(com.max_thread) private(y,x) schedule(static)
#endif
	for (y = 0; y < gfit.ry; y++) {
		for (x = 0; x < gfit.rx; x++) {
			guint src_index = y * gfit.rx + x;
			BYTE dst_pixel_value;
			WORD tmp_pixel_value;
			if (mode == HISTEQ_DISPLAY || mode == STF_DISPLAY)	// special case, no lo & hi
				dst_pixel_value = index[src[src_index]];
			else if (do_cut_over && src[src_index] > hi)	// cut
				dst_pixel_value = 0;
			else {
				if (src[src_index] - lo < 0)
					tmp_pixel_value = 0;
				else
					tmp_pixel_value = src[src_index] - lo;
				dst_pixel_value = index[tmp_pixel_value];
			}
			if (inverted)
				dst_pixel_value = UCHAR_MAX - dst_pixel_value;

			guint dst_index = ((gfit.ry - 1 - y) * gfit.rx + x) * 4;
			switch (color) {
				default:
				case NORMAL_COLOR:
					dst[dst_index++] = dst_pixel_value;
					dst[dst_index++] = dst_pixel_value;
					dst[dst_index++] = dst_pixel_value;
					break;
				case RAINBOW_COLOR:
					dst[dst_index++] = rainbow_index[dst_pixel_value][0];
					dst[dst_index++] = rainbow_index[dst_pixel_value][1];
					dst[dst_index++] = rainbow_index[dst_pixel_value][2];
			}
		}
	}

	// flush to ensure all writing to the image was done and redraw the surface
	cairo_surface_flush(com.surface[vport]);
	cairo_surface_mark_dirty(com.surface[vport]);

	test_and_allocate_reference_image(vport);
}

static int make_index_for_current_display(display_mode mode, WORD lo, WORD hi,
		int vport) {
	float pente;
	int i;
	BYTE *index;
	double m = 0.0;
	double pxl, shadows = 0.0, highlights = 0.0;
	if (mode == STF_DISPLAY)
		m = findMidtonesBalance(&gfit, &shadows, &highlights);

	/* initialization of data required to build the remap_index */
	switch (mode) {
	case NORMAL_DISPLAY:
		pente = UCHAR_MAX_SINGLE / (float) (hi - lo);
		break;
	case LOG_DISPLAY:
		pente = fabsf(UCHAR_MAX_SINGLE / logf(((float) (hi - lo)) * 0.1f));
		break;
	case SQRT_DISPLAY:
		pente = UCHAR_MAX_SINGLE / sqrtf((float) (hi - lo));
		break;
	case SQUARED_DISPLAY:
		pente = UCHAR_MAX_SINGLE / SQR((float )(hi - lo));
		break;
	case ASINH_DISPLAY:
		pente = UCHAR_MAX_SINGLE / asinhf(((float) (hi - lo)) * 0.001f);
		break;
	case STF_DISPLAY:
		pente = UCHAR_MAX_SINGLE;
		break;
	default:
		return 1;
	}
	if ((mode != HISTEQ_DISPLAY && mode != STF_DISPLAY) && pente == last_pente[vport]
			&& mode == last_mode[vport]) {
		fprintf(stdout, "Re-using previous remap_index\n");
		return 0;
	}
	fprintf(stdout, "Rebuilding remap_index\n");

	/************* Building the remap_index **************/
	if (!remap_index[vport]) {
		remap_index[vport] = malloc(USHRT_MAX + 1);
		if (!remap_index[vport]) {
			fprintf(stderr,
					"allocation error in remap_index, aborting remap\n");
			return 1;
		}
	}
	index = remap_index[vport];

	for (i = 0; i <= USHRT_MAX; i++) {
		switch (mode) {
		case LOG_DISPLAY:
			// ln(5.56*10^110) = 255
			if (i < 10)
				index[i] = 0; /* avoid null and negative values */
			else
				index[i] = round_to_BYTE(logf((float) i / 10.f) * pente); //10.f is arbitrary: good matching with ds9
			break;
		case SQRT_DISPLAY:
			// sqrt(2^16) = 2^8
			index[i] = round_to_BYTE(sqrtf((float) i) * pente);
			break;
		case SQUARED_DISPLAY:
			// pow(2^4,2) = 2^8
			index[i] = round_to_BYTE(SQR((float)i) * pente);
			break;
		case ASINH_DISPLAY:
			// asinh(2.78*10^110) = 255
			index[i] = round_to_BYTE(asinhf((float) i / 1000.f) * pente); //1000.f is arbitrary: good matching with ds9, could be asinhf(a*Q*i)/Q
			break;
		case NORMAL_DISPLAY:
			index[i] = round_to_BYTE((float) i * pente);
			break;
		case STF_DISPLAY:
			pxl = (gfit.bitpix == BYTE_IMG ?
					(double) i / UCHAR_MAX_DOUBLE :
					(double) i / USHRT_MAX_DOUBLE);
			pxl = (pxl - shadows < 0.0) ? 0.0 : pxl - shadows;
			pxl /= (highlights - shadows);
			index[i] = round_to_BYTE((float) (MTF(pxl, m)) * pente);
			break;
		default:
			return 1;
		}
		// check for maximum overflow, given that df/di > 0. Should not happen with round_to_BYTE
		if (index[i] == UCHAR_MAX)
			break;
	}
	if (i != USHRT_MAX + 1) {
		/* no more computation needed, just fill with max value */
		for (++i; i <= USHRT_MAX; i++)
			index[i] = UCHAR_MAX;
	}

	last_pente[vport] = pente;
	last_mode[vport] = mode;
	return 0;
}

static int make_index_for_rainbow(BYTE index[][3]) {
	int i;
	double h, s, v, r, g, b;

	for (i = 0; i < UCHAR_MAX + 1; i++) {
		r = g = b = (double) i / UCHAR_MAX_DOUBLE;
		rgb_to_hsv(r, g, b, &h, &s, &v);
		double off = 300.0 / 360.0;  /* Arbitrary: we want h from 300 to 0 deg */
		h = (off - (double) i * (off / UCHAR_MAX_DOUBLE));
		s = 1.;
		v = 1.; /* Saturation and Value are set to 100%  */
		hsv_to_rgb(h, s, v, &r, &g, &b);
		index[i][0] = round_to_BYTE(r * UCHAR_MAX_DOUBLE);
		index[i][1] = round_to_BYTE(g * UCHAR_MAX_DOUBLE);
		index[i][2] = round_to_BYTE(b * UCHAR_MAX_DOUBLE);
	}
	return 0;
}

/*
 * Free reference image
 */

static void free_reference_image() {
	fprintf(stdout, "Purging previously saved reference frame data.\n");
	if (com.refimage_regbuffer) {
		free(com.refimage_regbuffer);
		com.refimage_regbuffer = NULL;
	}
	if (com.refimage_surface) {
		cairo_surface_destroy(com.refimage_surface);
		com.refimage_surface = NULL;
	}
	enable_view_reference_checkbox(FALSE);
}

/*
 * Main conversion list static functions
 */

static GtkListStore *liststore_convert = NULL;

static void add_convert_to_list(char *filename, struct stat st) {
	GtkTreeIter iter;
	char *date;

	date = ctime(&st.st_mtime);
	date[strlen(date) - 1] = 0;	// removing '\n' at the end of the string

	gtk_list_store_append(liststore_convert, &iter);
	gtk_list_store_set(liststore_convert, &iter, COLUMN_FILENAME, filename,	// copied in the store
			COLUMN_DATE, date, -1);
}
static void get_convert_list_store() {
	if (liststore_convert == NULL)
		liststore_convert = GTK_LIST_STORE(
				gtk_builder_get_object(builder, "liststore_convert"));
}

static void fill_convert_list(GSList *list) {
	struct stat st;

	get_convert_list_store();

	while (list) {
		char *filename;

		filename = (char *) list->data;
		if (stat(filename, &st) == 0) {
			add_convert_to_list(filename, st);
			list = list->next;
		} else
			break;	// no infinite loop
		g_free(filename);
	}
	check_for_conversion_form_completeness();
}

/*
 * GTK File Chooser static functions
 */

static int whichdial;

static void gtk_filter_add(GtkFileChooser *file_chooser, const gchar *title,
		const gchar *pattern, gboolean set_default) {
	gchar **patterns;
	gint i;

	GtkFileFilter *f = gtk_file_filter_new();
	gtk_file_filter_set_name(f, title);
	/* get the patterns */
	patterns = g_strsplit(pattern, ";", -1);
	for (i = 0; patterns[i] != NULL; i++)
		gtk_file_filter_add_pattern(f, patterns[i]);
	/* free the patterns */
	g_strfreev(patterns);
	gtk_file_chooser_add_filter(file_chooser, f);
	if (set_default)
		gtk_file_chooser_set_filter(file_chooser, f);
}

static void set_filters_dialog(GtkFileChooser *chooser) {
	gtk_filter_add(chooser, _("FITS Files (*.fit, *.fits, *.fts)"),
			"*.fit;*.FIT;*.fits;*.FITS;*.fts;*.FTS",
			com.filter == TYPEFITS);
	if (whichdial == OD_OPEN || whichdial == OD_CONVERT) {
#ifdef HAVE_LIBRAW
		/* RAW FILES */
		int nb_raw;
		char *raw;
		int i;

		nb_raw = get_nb_raw_supported();
		raw = calloc(sizeof(char), nb_raw * 12 + 1);// we assume the extension size of 3 char "*.xxx;*.XXX;" = 12
		for (i = 0; i < nb_raw; i++) {
			char ext[20];
			gchar *upcase;

			upcase = g_ascii_strup(supported_raw[i].extension, strlen(supported_raw[i].extension));
			g_snprintf(ext, sizeof(ext), "*.%s;*.%s;",
					supported_raw[i].extension, upcase);
			strcat(raw, ext);

			g_free(upcase);
		}
		gtk_filter_add(chooser, _("RAW DSLR Camera Files"), raw,
				com.filter == TYPERAW);
		free(raw);
#endif
		/*GRAPHICS FILES*/
		char graphics_supported[256], pattern[256];
		g_snprintf(graphics_supported, sizeof(graphics_supported),
				_("Graphics Files (*.bmp"));
		g_snprintf(pattern, sizeof(pattern), "*.bmp;*.BMP;");
#ifdef HAVE_LIBJPEG
		strcat(graphics_supported, ", *.jpg, *.jpeg");
		strcat(pattern, "*.jpg;*.JPG;*.jpeg;*.JPEG;");
#endif

#ifdef HAVE_LIBPNG
		strcat(graphics_supported, ", *.png");
		strcat(pattern, "*.png;*.PNG;");
#endif

#ifdef HAVE_LIBTIFF
		strcat(graphics_supported, ", *.tif, *.tiff");
		strcat(pattern, "*.tif;*.TIF;*.tiff;*.TIFF");
#endif
		strcat(graphics_supported, ")");
		gtk_filter_add(chooser, graphics_supported, pattern,
				com.filter == TYPEBMP || com.filter == TYPEJPG
						|| com.filter == TYPEPNG || com.filter == TYPETIFF);

		/*NETPBM FILES*/
		gtk_filter_add(chooser, _("Netpbm Files (*.ppm, *.pnm, *.pgm)"),
				"*.ppm;*.PPM;*.pnm:*.PNM;*.pgm;*.PGM", com.filter == TYPEPNM);
		/*IRIS FILES*/
		gtk_filter_add(chooser, _("IRIS PIC Files (*.pic)"), "*.pic;*.PIC",
				com.filter == TYPEPIC);
		/* SER FILES */
		gtk_filter_add(chooser, _("SER files (*.ser)"), "*.ser;*.SER",
				com.filter == TYPESER);

#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
		/* FILM FILES */
		int nb_film;
		char *film;
		int j;

		nb_film = get_nb_film_ext_supported();
		film = calloc(sizeof(char), nb_film * 14 + 1);// we assume the extension size of 4 char "*.xxxx;*.XXXX;" = 14
		for (j = 0; j < nb_film; j++) {
			char ext[20];
			gchar *upcase;

			upcase = g_ascii_strup(supported_film[j].extension,
					strlen(supported_film[j].extension));
			g_snprintf(ext, sizeof(ext), "*.%s;*.%s;",
					supported_film[j].extension, upcase);
			strcat(film, ext);

			g_free(upcase);
		}
		gtk_filter_add(chooser, _("Film Files (*.avi, *.mpg, ...)"), film,
				com.filter == TYPEAVI);
		free(film);
#endif
	}
}

static void opendial(void) {
	GtkWidget *widgetdialog = NULL;
	GtkFileChooser *dialog = NULL;
	gint res;
	GtkWindow *main_window = GTK_WINDOW(
			gtk_builder_get_object(builder, "main_window"));
	GtkWindow *control_window = GTK_WINDOW(
			gtk_builder_get_object(builder, "control_window"));

	if (!com.wd)
		return;

	switch (whichdial) {
	case OD_NULL:
		fprintf(stderr, "whichdial undefined, should not happen\n");
		return;
	case OD_FLAT:
	case OD_DARK:
	case OD_OFFSET:
		widgetdialog = gtk_file_chooser_dialog_new(_("Open File"), control_window,
				GTK_FILE_CHOOSER_ACTION_OPEN, _("_Cancel"), GTK_RESPONSE_CANCEL,
				_("_Open"), GTK_RESPONSE_ACCEPT,
				NULL);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		gtk_file_chooser_set_current_folder(dialog, com.wd);
		gtk_file_chooser_set_select_multiple(dialog, FALSE);
		set_filters_dialog(dialog);
		break;
	case OD_CWD:
		widgetdialog = gtk_file_chooser_dialog_new(_("Open File"), control_window,
				GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER, _("_Cancel"),
				GTK_RESPONSE_CANCEL, _("_Open"), GTK_RESPONSE_ACCEPT,
				NULL);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		gtk_file_chooser_set_current_folder(dialog, com.wd);
		gtk_file_chooser_set_select_multiple(dialog, FALSE);
		break;
	case OD_OPEN:
		widgetdialog = gtk_file_chooser_dialog_new(_("Open File"), main_window,
				GTK_FILE_CHOOSER_ACTION_OPEN, _("_Cancel"), GTK_RESPONSE_CANCEL,
				_("_Open"), GTK_RESPONSE_ACCEPT,
				NULL);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		gtk_file_chooser_set_current_folder(dialog, com.wd);
		gtk_file_chooser_set_select_multiple(dialog, FALSE);
		set_filters_dialog(dialog);
		break;
	case OD_CONVERT:
		widgetdialog = gtk_file_chooser_dialog_new(_("Open File"), control_window,
				GTK_FILE_CHOOSER_ACTION_OPEN, _("_Cancel"), GTK_RESPONSE_CANCEL,
				_("_Open"), GTK_RESPONSE_ACCEPT,
				NULL);
		dialog = GTK_FILE_CHOOSER(widgetdialog);
		gtk_file_chooser_set_current_folder(dialog, com.wd);
		gtk_file_chooser_set_select_multiple(dialog, TRUE);
		set_filters_dialog(dialog);
	}

	if (!dialog)
		return;
	res = gtk_dialog_run(GTK_DIALOG(dialog));

	if (res == GTK_RESPONSE_ACCEPT) {
		GSList *list = NULL;
		char *filename;
		GtkFileChooser *chooser = GTK_FILE_CHOOSER(dialog);
		filename = gtk_file_chooser_get_filename(chooser);

		if (!(filename))
			return;

		switch (whichdial) {
		case OD_FLAT:
			gtk_entry_set_text(
					GTK_ENTRY(
							gtk_builder_get_object(builder, "flatname_entry")),
					filename);
			gtk_toggle_button_set_active(
					GTK_TOGGLE_BUTTON(
							gtk_builder_get_object(builder, "useflat_button")),
					TRUE);
			if (sequence_is_loaded() || single_image_is_loaded())
				gtk_widget_set_sensitive(lookup_widget("prepro_button"), TRUE);
			break;

		case OD_DARK:
			gtk_entry_set_text(
					GTK_ENTRY(
							gtk_builder_get_object(builder, "darkname_entry")),
					filename);
			gtk_toggle_button_set_active(
					GTK_TOGGLE_BUTTON(
							gtk_builder_get_object(builder, "usedark_button")),
					TRUE);
			if (sequence_is_loaded() || single_image_is_loaded())
				gtk_widget_set_sensitive(lookup_widget("prepro_button"), TRUE);
			break;

		case OD_OFFSET:
			gtk_entry_set_text(
					GTK_ENTRY(
							gtk_builder_get_object(builder,
									"offsetname_entry")), filename);
			gtk_toggle_button_set_active(
					GTK_TOGGLE_BUTTON(
							gtk_builder_get_object(builder,
									"useoffset_button")),
					TRUE);
			if (sequence_is_loaded() || single_image_is_loaded())
				gtk_widget_set_sensitive(lookup_widget("prepro_button"), TRUE);
			break;

		case OD_CWD:
			if (!changedir(filename))
				writeinitfile();
			break;

		case OD_OPEN:
			set_cursor_waiting(TRUE);
			open_single_image(filename);
			set_cursor_waiting(FALSE);
			break;
		case OD_CONVERT:
			list = gtk_file_chooser_get_filenames(chooser);
			fill_convert_list(list);
			g_slist_free(list);
			break;
		}
		whichdial = OD_NULL;
		g_free(filename);
	}
	gtk_widget_destroy(widgetdialog);
}

static void Set_Programm_name_in_TIFF() {
	static GtkTextView *TIFF_txt = NULL;
	GtkTextBuffer *tbuf;
	GtkTextIter itDebut, itFin;
	char Copyright[64];

	if (TIFF_txt == NULL)
		TIFF_txt = GTK_TEXT_VIEW(lookup_widget("Copyright_txt"));

	tbuf = gtk_text_view_get_buffer(TIFF_txt);

	g_snprintf(Copyright, sizeof(Copyright), "%s v%s", PACKAGE, VERSION);
	Copyright[0] = toupper(Copyright[0]);			// convert siril to Siril

	gtk_text_buffer_get_bounds(tbuf, &itDebut, &itFin);
	gtk_text_buffer_delete(tbuf, &itDebut, &itFin);
	gtk_text_buffer_set_text(tbuf, Copyright, strlen(Copyright));
}

static int savedial(char *filename, const gchar *title, const gchar *pattern) {
	GtkWidget *dialog;
	GtkFileChooser *chooser;
	GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_SAVE;
	GtkWindow *parent_window = GTK_WINDOW(
			gtk_builder_get_object(builder, "savepopup"));
	gint res, retval = 0;

	dialog = gtk_file_chooser_dialog_new(_("Save File"), parent_window, action,
			_("_Cancel"), GTK_RESPONSE_CANCEL, _("_Save"), GTK_RESPONSE_ACCEPT,
			NULL);

	chooser = GTK_FILE_CHOOSER(dialog);

	gtk_file_chooser_set_do_overwrite_confirmation(chooser, TRUE);
	gtk_filter_add(chooser, title, pattern, FALSE);
	gtk_file_chooser_set_filename(chooser, filename);
	gtk_file_chooser_set_current_name(chooser, (filename));

	res = gtk_dialog_run(GTK_DIALOG(dialog));
	if (res == GTK_RESPONSE_ACCEPT) {
		char *new_filename;

		new_filename = gtk_file_chooser_get_filename(chooser);
		strcpy(filename, new_filename);
		g_free(new_filename);
		retval = 1;
	}

	gtk_widget_destroy(dialog);
	return retval;
}

static image_type whichminisave;

static void minisavedial(void) {
	const gchar *name;
	gchar filename[256];
	GtkToggleButton *fits_8 = GTK_TOGGLE_BUTTON(
			lookup_widget("radiobutton_save_fit8"));
	GtkToggleButton *fits_16s = GTK_TOGGLE_BUTTON(
			lookup_widget("radiobutton_save_fit16s"));
#ifdef HAVE_LIBJPEG
	GtkSpinButton *qlty_spin_button = GTK_SPIN_BUTTON(
			lookup_widget("quality_spinbutton"));
	gint quality = gtk_spin_button_get_value_as_int(qlty_spin_button);
#endif
#ifdef HAVE_LIBTIFF
	gint bitspersamples = 16;
	GtkToggleButton *BPS_Button = GTK_TOGGLE_BUTTON(
			lookup_widget("radiobutton8bits"));

	if (gtk_toggle_button_get_active(BPS_Button))
		bitspersamples = 8;
#endif
	GtkEntry *entry = GTK_ENTRY(lookup_widget("savetxt"));

	name = gtk_entry_get_text(entry);
	if (name[0] != '\0') {
		int nplanes;

		strcpy(filename, name);
		switch (whichminisave) {
		case TYPEBMP:
			strcat(filename, ".bmp");
			if (savedial(filename, _("BMP Files"), "*.bmp;*.BMP"))
				savebmp(filename, &gfit);
			break;
#ifdef HAVE_LIBJPEG
		case TYPEJPG:
			strcat(filename, ".jpg");
			if (savedial(filename, _("JPEG Files"), "*.jpg;*.JPG;*.jpeg;*.JPEG"))
				savejpg(filename, &gfit, quality);
			break;
#endif
#ifdef HAVE_LIBTIFF
		case TYPETIFF:
			strcat(filename, ".tif");
			if (savedial(filename, _("TIFF Files"), "*.tif;*.TIF;*.tiff;*.TIFF"))
				savetif(filename, &gfit, bitspersamples);
			break;
#endif
		case TYPEFITS:
			if (gtk_toggle_button_get_active(fits_8))
				gfit.bitpix = BYTE_IMG;
			else if (gtk_toggle_button_get_active(fits_16s))
				gfit.bitpix = SHORT_IMG;
			else
				gfit.bitpix = USHORT_IMG;
			/* Check if MIPS-HI and MIPS-LO must be updated. If yes,
			 * Values are taken from the layer 0 */
			if (gtk_toggle_button_get_active(
					GTK_TOGGLE_BUTTON(
							lookup_widget(
									"checkbutton_update_hilo"))) == TRUE) {
				if (sequence_is_loaded() && !single_image_is_loaded()) {
					gfit.hi = com.seq.layers[RLAYER].hi;
					gfit.lo = com.seq.layers[RLAYER].lo;
				} else {
					gfit.hi = com.uniq->layers[RLAYER].hi;
					gfit.lo = com.uniq->layers[RLAYER].lo;
				}
				if (gfit.bitpix == BYTE_IMG
						&& (gfit.hi > UCHAR_MAX || gfit.lo > UCHAR_MAX)) {
					gfit.hi = UCHAR_MAX;
					gfit.lo = 0;
				} else if (gfit.bitpix == SHORT_IMG
						&& (gfit.hi > SHRT_MAX || gfit.lo > SHRT_MAX)) {
					gfit.hi = UCHAR_MAX;
					gfit.lo = 0;
				}

			}
			strcat(filename, ".fit");
			if (savedial(filename, _("FITS Files"),
					"*.fit;*.FIT;*.fts;*.FTS;*.fits;*.FITS"))
				savefits(filename, &gfit);
			break;
		case TYPEPNM:
			nplanes = gfit.naxes[2];
			if (nplanes == 1) {
				strcat(filename, ".pgm");
				if (savedial(filename, _("NetPBM Files"), "*.pgm;*.PGM"))
					savepgm(filename, &gfit);
			} else if (nplanes == 3) {
				strcat(filename, ".ppm");
				if (savedial(filename, _("NetPBM Files"), "*.ppm;*.PPM"))
					saveppm(filename, &gfit);
			} else
				return;			//should not happen
			break;
		default:
			siril_log_message(
					_("This type of file is not handled. Should not happen"));
			break;
		}
		gtk_widget_hide(lookup_widget("savepopup"));
		gtk_entry_set_text(entry, "");
	}
}

/*
 * Update FWHMP static function
 */

static void update_fwhm_units_ok() {
	GtkWidget *label_ok = GTK_WIDGET(lookup_widget("label_ok"));

	gtk_widget_set_visible(label_ok,
			gfit.focal_length > 0.0 && gfit.pixel_size_x > 0.0f
					&& gfit.pixel_size_y > 0.0f);
}

/*
 * Set swap dir to default
 */

static void reset_swapdir() {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(lookup_widget("filechooser_swap"));
	GtkLabel *label = GTK_LABEL(lookup_widget("label_swap_dir"));
	const char *dir;

	dir = g_get_tmp_dir();

	if (strcmp(dir, com.swap_dir)) {
		if (com.swap_dir)
			free(com.swap_dir);
		com.swap_dir = strdup(dir);
		gtk_file_chooser_set_filename(swap_dir, dir);
		gtk_label_set_text(label, dir);
		writeinitfile();
	}
}

/*
 * Command line history static function
 */

static void history_add_line(char *line) {
	if (!com.cmd_history) {
		com.cmd_hist_size = CMD_HISTORY_SIZE;
		com.cmd_history = calloc(com.cmd_hist_size, sizeof(const char*));
		com.cmd_hist_current = 0;
		com.cmd_hist_display = 0;
	}
	com.cmd_history[com.cmd_hist_current] = line;
	com.cmd_hist_current++;
	// circle at the end
	if (com.cmd_hist_current == com.cmd_hist_size)
		com.cmd_hist_current = 0;
	if (com.cmd_history[com.cmd_hist_current]) {
		free(com.cmd_history[com.cmd_hist_current]);
		com.cmd_history[com.cmd_hist_current] = NULL;
	}
	com.cmd_hist_display = com.cmd_hist_current;
}

/*
 * Selection static functions
 */

/* selection zone event management */
#define MAX_CALLBACKS_PER_EVENT 10
static selection_update_callback _registered_callbacks[MAX_CALLBACKS_PER_EVENT];
static int _nb_registered_callbacks = 0;

void register_selection_update_callback(selection_update_callback f) {
	if (_nb_registered_callbacks < MAX_CALLBACKS_PER_EVENT) {
		_registered_callbacks[_nb_registered_callbacks] = f;
		_nb_registered_callbacks++;
	}
}

void unregister_selection_update_callback(selection_update_callback f) {
	int i;
	for (i = 0; i < _nb_registered_callbacks; ++i) {
		if (_registered_callbacks[i] == f) {
			_registered_callbacks[i] =
					_registered_callbacks[_nb_registered_callbacks];
			_registered_callbacks[_nb_registered_callbacks] = NULL;
			_nb_registered_callbacks--;
			return;
		}
	}
}

// send the events
static void new_selection_zone() {
	int i;
	fprintf(stdout, "selection: %d,%d,\t%dx%d\n", com.selection.x, com.selection.y,
			com.selection.w, com.selection.h);
	for (i = 0; i < _nb_registered_callbacks; ++i) {
		_registered_callbacks[i]();
	}
}

void delete_selected_area() {
	memset(&com.selection, 0, sizeof(rectangle));
	new_selection_zone();
}

/*
 * MISC static functions
 */

static void toggle_image_selection(int image_num) {
	char msg[60];
	if (com.seq.imgparam[image_num].incl) {
		com.seq.imgparam[image_num].incl = FALSE;
		--com.seq.selnum;
		g_snprintf(msg, sizeof(msg),
				_("Image %d has been unselected from sequence\n"), image_num);
	} else {
		com.seq.imgparam[image_num].incl = TRUE;
		++com.seq.selnum;
		g_snprintf(msg, sizeof(msg),
				_("Image %d has been selected from sequence\n"), image_num);
	}
	siril_log_message(msg);
	sequence_list_change_selection_index(image_num);
	update_reg_interface(FALSE);
	adjust_exclude(image_num, TRUE);
	writeseqfile(&com.seq);
}

static int get_index_in_predefined_colors_for_wavelength(double wl) {
	int i;
	for (i = 0; i < sizeof(predefined_layers_colors) / sizeof(layer_info);
			i++) {
		if (predefined_layers_colors[i].wavelength == wl)
			return i;
	}
	return -1;
}

/* method handling all include or all exclude from a sequence */
static void sequence_setselect_all(gboolean include_all) {
	int i;

	if (!com.seq.imgparam)
		return;
	for (i = 0; i <= com.seq.number; ++i) {
		if (com.seq.imgparam[i].incl != include_all) {
			com.seq.imgparam[i].incl = include_all;
			sequence_list_change_selection_index(i);
		}
	}
	if (include_all) {
		com.seq.selnum = com.seq.number;
		siril_log_message(_("Selected all images from sequence\n"));
	} else {
		com.seq.selnum = 0;
		siril_log_message(_("Unselected all images from sequence\n"));
	}
	adjust_exclude(com.seq.current, TRUE);
	update_reg_interface(FALSE);
	writeseqfile(&com.seq);
}

/* RGB popup menu */

static void do_popup_rgbmenu(GtkWidget *my_widget, GdkEventButton *event) {
	static GtkMenu *menu = NULL;

	if (!menu) {
		menu = GTK_MENU(gtk_builder_get_object(builder, "menurgb"));
		gtk_menu_attach_to_widget(GTK_MENU(menu), my_widget, NULL);
	}

#if GTK_MAJOR_VERSION >= 3 && GTK_MINOR_VERSION < 22
	int button, event_time;

	if (event) {
		button = event->button;
		event_time = event->time;
	} else {
		button = 0;
		event_time = gtk_get_current_event_time();
	}

	gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL, button,
			event_time);
#else
	gtk_menu_popup_at_pointer(GTK_MENU(menu), NULL);
#endif
}

/* Gray popup menu */

static void do_popup_graymenu(GtkWidget *my_widget, GdkEventButton *event) {
	static GtkMenu *menu = NULL;
	gboolean selected;
	gboolean is_a_single_image_loaded = single_image_is_loaded() && (!sequence_is_loaded()
				|| (sequence_is_loaded() && com.seq.current == RESULT_IMAGE));

	if (!menu) {
		menu = GTK_MENU(gtk_builder_get_object(builder, "menugray"));
		gtk_menu_attach_to_widget(GTK_MENU(menu), my_widget, NULL);
	}

	selected = com.selection.w && com.selection.h;
	gtk_widget_set_sensitive(lookup_widget("undo_item1"), is_undo_available());
	gtk_widget_set_sensitive(lookup_widget("redo_item1"), is_redo_available());
	gtk_widget_set_sensitive(lookup_widget("menu_gray_psf"), selected);
	gtk_widget_set_sensitive(lookup_widget("menu_gray_seqpsf"), selected);
	gtk_widget_set_sensitive(lookup_widget("menu_gray_pick_star"), selected);
	gtk_widget_set_sensitive(lookup_widget("menu_gray_crop"), selected && is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menu_gray_crop_seq"), selected && sequence_is_loaded());

#if GTK_MAJOR_VERSION >= 3 && GTK_MINOR_VERSION < 22
	int button, event_time;

	if (event) {
		button = event->button;
		event_time = event->time;
	} else {
		button = 0;
		event_time = gtk_get_current_event_time();
	}

	gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL, button, event_time);
#else
	gtk_menu_popup_at_pointer(GTK_MENU(menu), NULL);
#endif
}

/*****************************************************************************
 *                    P U B L I C      F U N C T I O N S                     *
 ****************************************************************************/

GtkWidget* lookup_widget(const gchar *widget_name) {
	return GTK_WIDGET(gtk_builder_get_object(builder, widget_name));
}

void set_sliders_value_to_gfit() {
	static GtkAdjustment *adj1 = NULL, *adj2 = NULL;

	if (adj1 == NULL) {
		adj1 = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment1"));// scalemax
		adj2 = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment2"));// scalemin
	}

	gfit.hi = gtk_adjustment_get_value(adj1);
	gfit.lo = gtk_adjustment_get_value(adj2);
}

/* Sets maximum value for contrast scales. Minimum is always 0.
 * Should be done on first image load of a sequence and when single images are loaded.
 * Max value is taken from gfit.maxi, recomputed if not present.
 *
 */
void set_cutoff_sliders_max_values() {
	static GtkAdjustment *adj1 = NULL, *adj2 = NULL;
	gdouble max_val;
	if (adj1 == NULL) {
		adj1 = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment1"));// scalemax
		adj2 = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment2"));// scalemin
	}
	fprintf(stdout, _("Setting MAX value for cutoff sliders adjustments\n"));
	/* set max value for range according to number of bits of original image
	 * We should use gfit.bitpix for this, but it's currently always USHORT_IMG */
	/*if (gfit.bitpix == BYTE_IMG)
	 max_val = UCHAR_MAX_DOUBLE;
	 else max_val = USHRT_MAX_DOUBLE;*/
	if (gfit.maxi == 0)
		image_find_minmax(&gfit, 0);
	if (gfit.maxi <= UCHAR_MAX)
		max_val = UCHAR_MAX_DOUBLE;
	else
		max_val = USHRT_MAX_DOUBLE;

	gtk_adjustment_set_upper(adj1, max_val);
	gtk_adjustment_set_upper(adj2, max_val);
}

/* Sets the value scalemin and scalemax sliders so that they match the hi and
 * lo values stored for the current viewport.
 * It also sets the 'cut' checkbox status.
 * The function should be called on every tab change and file opening */
void set_cutoff_sliders_values() {
	gchar buffer[10];
	static GtkAdjustment *adjmin = NULL, *adjmax = NULL;
	static GtkEntry *maxentry = NULL, *minentry = NULL;
	static GtkToggleButton *cutmax = NULL;
	WORD hi, lo;
	gboolean cut_over;
	if (adjmin == NULL) {
		adjmax = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment1")); // scalemax
		adjmin = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjustment2")); // scalemin
		maxentry = GTK_ENTRY(gtk_builder_get_object(builder, "max_entry"));
		minentry = GTK_ENTRY(gtk_builder_get_object(builder, "min_entry"));
		cutmax = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "checkcut_max"));
	}
	if (single_image_is_loaded() && com.cvport < com.uniq->nb_layers &&
			com.uniq->layers && com.seq.current != RESULT_IMAGE) {
		hi = com.uniq->layers[com.cvport].hi;
		lo = com.uniq->layers[com.cvport].lo;
		cut_over = com.uniq->layers[com.cvport].cut_over;
	}
	/* When com.seq.current == RESULT_IMAGE we take sequence values */
	else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers
			&& com.seq.layers) {
		hi = com.seq.layers[com.cvport].hi;
		lo = com.seq.layers[com.cvport].lo;
		cut_over = com.seq.layers[com.cvport].cut_over;
	} else
		return;	// there should be no other normal cases
	fprintf(stdout, _("setting ranges scalemin=%d, scalemax=%d\n"), lo, hi);
	WORD maxvalue = get_normalized_value(&gfit);
	gtk_adjustment_set_lower(adjmin, 0.0);
	gtk_adjustment_set_lower(adjmax, 0.0);
	gtk_adjustment_set_upper(adjmin, (double) maxvalue);
	gtk_adjustment_set_upper(adjmax, (double) maxvalue);
	gtk_adjustment_set_value(adjmin, (gdouble) lo);
	gtk_adjustment_set_value(adjmax, (gdouble) hi);
	g_snprintf(buffer, 6, "%u", hi);
	g_signal_handlers_block_by_func(maxentry, on_max_entry_changed, NULL);
	gtk_entry_set_text(maxentry, buffer);
	g_signal_handlers_unblock_by_func(maxentry, on_max_entry_changed, NULL);
	g_snprintf(buffer, 6, "%u", lo);
	g_signal_handlers_block_by_func(minentry, on_min_entry_changed, NULL);
	gtk_entry_set_text(minentry, buffer);
	g_signal_handlers_unblock_by_func(minentry, on_min_entry_changed, NULL);
	gtk_toggle_button_set_active(cutmax, cut_over);
}

/* sets the maximum value for the spin button and display the initial file name */
int seqsetnum(int image_number) {
	GtkSpinButton *spin;
	GtkAdjustment *adj;
	if (com.seq.number <= 0 || image_number >= com.seq.number)
		return 1;
	spin = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "imagenumber_spin"));
	adj = gtk_spin_button_get_adjustment(spin);

	gtk_adjustment_set_upper(adj, (gdouble) com.seq.number - 1);
	gtk_adjustment_set_value(adj, (gdouble) image_number);	// 0 is default
	display_image_number(image_number);
	//on_imagenumberspin_output(GTK_SPIN_BUTTON(spin), NULL);	// redraw the real file number instead of 0
	return 0;
}

/* Sets the display mode combo box to the value stored in the relevant struct.
 * The operation is purely graphical. */
void set_display_mode() {
	static GtkComboBox *modecombo = NULL;
	display_mode mode;

	if (!modecombo)
		modecombo = GTK_COMBO_BOX(lookup_widget("combodisplay"));

	if (single_image_is_loaded() && com.cvport < com.uniq->nb_layers && com.uniq->layers
	&& com.seq.current != RESULT_IMAGE)
		mode = com.uniq->layers[com.cvport].rendering_mode;
	else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers
			&& com.seq.layers)
		mode = com.seq.layers[com.cvport].rendering_mode;
	else
		return;

	g_signal_handlers_block_by_func(modecombo, on_combodisplay_changed, NULL);
	gtk_combo_box_set_active(modecombo, mode);
	g_signal_handlers_unblock_by_func(modecombo, on_combodisplay_changed, NULL);
}

/* called when an image is included or excluded from the sequence.
 * it sets the toggle button "exclude_button" so it matches the inclusion state.
 * n is the image number in the sequence
 * changed indicates if a value was changed and if the display needs to be refreshed.
 */
void adjust_exclude(int n, gboolean changed) {
	static GtkWidget *excl_butt = NULL;
	if (!com.seq.imgparam || n < 0 || n >= com.seq.number)
		return;
	if (excl_butt == NULL)
		excl_butt = lookup_widget("exclude_button");

	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(excl_butt))
			== com.seq.imgparam[n].incl) {
		g_signal_handlers_block_by_func(excl_butt, on_excludebutton_toggled,
				NULL);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(excl_butt),
				!com.seq.imgparam[n].incl);
		g_signal_handlers_unblock_by_func(excl_butt, on_excludebutton_toggled,
				NULL);
	}

	if (changed) {
		redraw(com.cvport, REMAP_NONE);
		drawPlot();
		adjust_sellabel();
	}
}

/* fill the label indicating how many images are selected in the gray and
 * which one is the reference image, at the bottom of the main window */
int adjust_sellabel() {
	static GtkLabel *local_label = NULL, *global_label = NULL;
	char bufferlocal[256], bufferglobal[256];
	if (local_label == NULL) {
		local_label = GTK_LABEL(lookup_widget("imagesel_label"));
		global_label = GTK_LABEL(lookup_widget("labelseq"));
	}
	if (sequence_is_loaded()) {
		if (com.seq.reference_image != -1) {
			char format[150];
			if (com.seq.fixed <= 1) {
				g_snprintf(format, sizeof(format),
						_("<%%s.seq>: %%d images selected out of %%d, reference image is %%d"));
			} else {
				g_snprintf(format, sizeof(format),
						_("<%%s.seq>: %%d images selected out of %%d, reference image is %%.%dd"),
						com.seq.fixed);
			}
			g_snprintf(bufferlocal, sizeof(bufferlocal), format,
					com.seq.seqname, com.seq.selnum, com.seq.number,
					com.seq.imgparam[com.seq.reference_image].filenum);

		} else {
			g_snprintf(bufferlocal, sizeof(bufferlocal),
					_("<%s.seq>: %d images selected out of %d, no reference image set"),
					com.seq.seqname, com.seq.selnum, com.seq.number);
		}
		g_snprintf(bufferglobal, sizeof(bufferglobal), _("%s, %d images selected"),
				com.seq.seqname, com.seq.selnum);
		//gtk_widget_set_sensitive(lookup_widget("goregister_button"), com.seq.selnum>0?TRUE:FALSE);
	} else {
		g_snprintf(bufferlocal, sizeof(bufferlocal), _("No sequence"));
		g_snprintf(bufferglobal, sizeof(bufferglobal), _("- none -"));
		gtk_widget_set_sensitive(lookup_widget("goregister_button"), FALSE);
	}

	gtk_label_set_text(local_label, bufferlocal);
	gtk_label_set_text(global_label, bufferglobal);
	return 0;
}

void update_MenuItem() {
	gboolean is_a_single_image_loaded;		/* An image is loaded. Not a sequence or only the result of stacking process */
	gboolean is_a_singleRGB_image_loaded;	/* A RGB image is laoded. Not a sequence or only the result of stacking process */
	gboolean any_image_is_loaded;			/* Something is loaded. Single image or Sequence */
	gboolean any_RGB_image_is_loaded;		/* Some RGB data are loaded. Single image or Sequence */

	is_a_singleRGB_image_loaded = isrgb(&gfit) && (!sequence_is_loaded()
			|| (sequence_is_loaded() && com.seq.current == RESULT_IMAGE));

	is_a_single_image_loaded = single_image_is_loaded()	&& (!sequence_is_loaded()
			|| (sequence_is_loaded() && com.seq.current == RESULT_IMAGE));

	any_image_is_loaded = single_image_is_loaded() || sequence_is_loaded();

	any_RGB_image_is_loaded = isrgb(&gfit) && (single_image_is_loaded() || sequence_is_loaded());

	/* File Menu */
	gtk_widget_set_sensitive(lookup_widget("menu_save_fits"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menu_save_tiff"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menu_save_bmp"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menu_save_jpg"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menu_save_pbm"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menu_FITS_header"), any_image_is_loaded && gfit.header != NULL);

	/* Edit Menu */
	gtk_widget_set_sensitive(lookup_widget("undo_item"), is_undo_available());
	gtk_widget_set_sensitive(lookup_widget("redo_item"), is_redo_available());

	/* Image processing Menu */
	gtk_widget_set_sensitive(lookup_widget("removegreen"), is_a_singleRGB_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_satu"), is_a_singleRGB_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitemcalibration"), is_a_singleRGB_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menu_channel_separation"), is_a_singleRGB_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_histo"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_fixbanding"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_cosmetic"), any_image_is_loaded);
#ifdef HAVE_OPENCV
	gtk_widget_set_sensitive(lookup_widget("menuitem_resample"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_rotation"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_rotation90"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_rotation270"), is_a_single_image_loaded);
#else
	gtk_widget_set_sensitive(lookup_widget("menuitem_resample"), FALSE);
	gtk_widget_set_sensitive(lookup_widget("menuitem_rotation"), FALSE);
	gtk_widget_set_sensitive(lookup_widget("menuitem_rotation90"), FALSE);
	gtk_widget_set_sensitive(lookup_widget("menuitem_rotation270"), FALSE);
#endif
	gtk_widget_set_sensitive(lookup_widget("menuitem_mirrorx"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_mirrory"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_bkg_extraction"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_wavelets"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menu_wavelet_separation"), is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_medianfilter"), is_a_single_image_loaded);

	/* Analysis Menu */
	gtk_widget_set_sensitive(lookup_widget("menuitem_noise"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitem_stat"), any_image_is_loaded);

	/* Windows Menu */
	gtk_widget_set_sensitive(lookup_widget("menuitemgray"), any_image_is_loaded);
	gtk_widget_set_sensitive(lookup_widget("menuitemcolor"), any_RGB_image_is_loaded);
}

gboolean redraw(int vport, int doremap) {
	GtkWidget *widget;

	if (vport >= MAXVPORT) {
		fprintf(stderr, _("redraw: maximum number of layers supported is %d"
				" (current image has %d).\n"), MAXVPORT, vport);
		return FALSE;
	}
	widget = com.vport[vport];

	switch (vport) {
	case RED_VPORT:
	case BLUE_VPORT:
	case GREEN_VPORT:
		if (doremap == REMAP_ONLY) {
			remap(vport);
		} else if (doremap == REMAP_ALL) {
			int i;
//#pragma omp parallel for num_threads(com.max_thread) private(i) schedule(static)		//probably causes crashes in HESTEQ_MODE
			for (i = 0; i < gfit.naxes[2]; i++) {
				remap(i);
			}
		}
		gtk_widget_queue_draw(widget);
		if (gfit.naxes[2] == 1)
			break;
		/* no break */
	case RGB_VPORT:
		if (gfit.naxis == 3) {
			if (doremap != REMAP_NONE) {
				remaprgb();
			}
			widget = com.vport[RGB_VPORT];
			gtk_widget_queue_draw(widget);
		}
		break;
	default:
		fprintf(stderr, "redraw: unknown viewport number %d\n", vport);
		break;
	}
	//fprintf(stdout, "end of redraw\n");
	com.drawn = FALSE;
	return FALSE;
}

void sliders_mode_set_state(sliders_mode sliders) {
	GtkToggleButton *radiobutton;	// Must not be static
	char *str[] =
			{ "radiobutton_hilo", "radiobutton_minmax", "radiobutton_user" };
	void *func[] = { on_radiobutton_hilo_toggled, on_radiobutton_minmax_toggled,
			on_radiobutton_user_toggled };

	radiobutton = GTK_TOGGLE_BUTTON(
			gtk_builder_get_object(builder, str[sliders]));

	g_signal_handlers_block_by_func(radiobutton, func[sliders], NULL);
	gtk_toggle_button_set_active(radiobutton, TRUE);
	g_signal_handlers_unblock_by_func(radiobutton, func[sliders], NULL);
}

/* When rendering settings are chained, they need to be copied to other layers
 * when modified on the current layer. This procedure does that. It can be
 * called whenever a value has changed, or when the chaning has been enabled,
 * to synchronize all data.
 * Current view port is com.cvport, and data is stored in layer_info structs
 * Synchronized data: hi and lo cursors, cut over box, rendering mode.
 * DOES NOT REMAP/REDRAW.
 *
 * from_GUI: TRUE if get values from the GUI, FALSE if get the values from structs.
 * Returns 1 if chained, 0 if not.
 */
int copy_rendering_settings_when_chained(gboolean from_GUI) {
	static GtkToggleButton *chainedbutton = NULL;
	static GtkRange *range_lo = NULL, *range_hi = NULL;
	static GtkComboBox *modecombo = NULL;
	static GtkToggleButton *cutmax = NULL;

	gboolean is_chained;
	display_mode mode;
	WORD lo, hi;
	gboolean cut_over;
	int i, nb_layers;
	layer_info *layers = NULL;

	if (!chainedbutton) {	// init widgets
		chainedbutton = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_chain"));
		modecombo = GTK_COMBO_BOX(lookup_widget("combodisplay"));
		range_lo = GTK_RANGE(gtk_builder_get_object(builder, "scalemin"));
		range_hi = GTK_RANGE(gtk_builder_get_object(builder, "scalemax"));
		cutmax = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "checkcut_max"));
	}

	is_chained = gtk_toggle_button_get_active(chainedbutton);
	if (single_image_is_loaded() &&
			com.cvport < com.uniq->nb_layers && com.uniq->layers &&
			com.seq.current != RESULT_IMAGE) {
		layers = com.uniq->layers;
		nb_layers = com.uniq->nb_layers;
	} else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers
			&& com.seq.layers) {
		layers = com.seq.layers;
		nb_layers = com.seq.nb_layers;
	} else
		return 0;

	if (from_GUI) {
		int raw_mode = gtk_combo_box_get_active(modecombo);
		/* update values in the layer_info for cvport */
		layers[com.cvport].rendering_mode =
				raw_mode >= 0 ? raw_mode : NORMAL_DISPLAY;
		layers[com.cvport].lo = round_to_WORD(gtk_range_get_value(range_lo));
		layers[com.cvport].hi = round_to_WORD(gtk_range_get_value(range_hi));
		layers[com.cvport].cut_over = gtk_toggle_button_get_active(cutmax);
	}
	if (!is_chained)
		return 0;
	mode = layers[com.cvport].rendering_mode;
	lo = layers[com.cvport].lo;
	hi = layers[com.cvport].hi;
	cut_over = layers[com.cvport].cut_over;

	for (i = 0; i < nb_layers; i++) {
		if (i == com.cvport)
			continue;
		layers[i].rendering_mode = mode;
		layers[i].lo = lo;
		layers[i].hi = hi;
		layers[i].cut_over = cut_over;
	}

	return 1;
}

void set_prepro_button_sensitiveness() {
	static GtkToggleButton *udark = NULL, *uoffset = NULL, *uflat = NULL,
			*checkAutoEvaluate = NULL;
	if (udark == NULL) {
		udark = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "usedark_button"));
		uoffset = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "useoffset_button"));
		uflat = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "useflat_button"));
		checkAutoEvaluate = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "checkbutton_auto_evaluate"));
	}

	gtk_widget_set_sensitive(lookup_widget("prepro_button"),
			(sequence_is_loaded() || single_image_is_loaded())
					&& (gtk_toggle_button_get_active(udark)
							|| gtk_toggle_button_get_active(uoffset)
							|| gtk_toggle_button_get_active(uflat)));
	gtk_widget_set_sensitive(lookup_widget("grid24"),
			gtk_toggle_button_get_active(udark));
	gtk_widget_set_sensitive(lookup_widget("checkDarkOptimize"),
			gtk_toggle_button_get_active(udark));
	gtk_widget_set_sensitive(lookup_widget("checkbutton_auto_evaluate"),
			gtk_toggle_button_get_active(uflat));
	gtk_widget_set_sensitive(lookup_widget("entry_flat_norm"),
			gtk_toggle_button_get_active(uflat)
					&& !gtk_toggle_button_get_active(checkAutoEvaluate));
}

void clear_sampling_setting_box() {
	GtkComboBox *binning = GTK_COMBO_BOX(
			gtk_builder_get_object(builder, "combobinning"));
	GtkEntry* focal_entry = GTK_ENTRY(lookup_widget("focal_entry"));
	GtkEntry* pitchX_entry = GTK_ENTRY(lookup_widget("pitchX_entry"));
	GtkEntry* pitchY_entry = GTK_ENTRY(lookup_widget("pitchY_entry"));

	gtk_entry_set_text(focal_entry, "");
	gtk_entry_set_text(pitchX_entry, "");
	gtk_entry_set_text(pitchY_entry, "");
	gtk_combo_box_set_active(binning, 0);
}

void update_libraw_interface() {
	/**********COLOR ADJUSTEMENT**************/
	com.raw_set.bright = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("Brightness_spinbutton")));
	com.raw_set.mul[0] = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("Red_spinbutton")));
	com.raw_set.mul[2] = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("Blue_spinbutton")));

	com.raw_set.auto_mul = gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_multipliers")));
	com.raw_set.user_black = gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_blackpoint")));

	/**************WHITE BALANCE**************/
	com.raw_set.use_camera_wb = (int) gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_cam")));
	com.raw_set.use_auto_wb = (int) gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_auto")));

	/********MATRIX INTERPOLATION**************/
	com.raw_set.user_qual = gtk_combo_box_get_active(
			GTK_COMBO_BOX(lookup_widget("combo_dcraw_inter")));

	/********GAMMA CORRECTION**************/
	if (gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm0"))) == TRUE) {
		/* Linear Gamma Curve */
		com.raw_set.gamm[0] = 1.0;
		com.raw_set.gamm[1] = 1.0;
	} else if (gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm1"))) == TRUE) {
		/* BT.709 Gamma curve */
		com.raw_set.gamm[0] = 2.222;
		com.raw_set.gamm[1] = 4.5;
	} else {
		/* sRGB Gamma curve */
		com.raw_set.gamm[0] = 2.40;
		com.raw_set.gamm[1] = 12.92;
	}
	/* We write in config file */
	/*************SER**********************/
	com.debayer.use_bayer_header = gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_SER_use_header")));
	com.debayer.compatibility = gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_debayer_compatibility")));
	writeinitfile();
}

char *vport_number_to_name(int vport) {
	switch (vport) {
	case RED_VPORT:
		return strdup("red");
	case GREEN_VPORT:
		return strdup("green");
	case BLUE_VPORT:
		return strdup("blue");
	case RGB_VPORT:
		return strdup("rgb");
	}
	return NULL;
}

int match_drawing_area_widget(GtkWidget *drawing_area, gboolean allow_rgb) {
	/* could be done with a for i=0 loop, to get rid of these defines */
	if (drawing_area == com.vport[RED_VPORT])
		return RED_VPORT;
	if (drawing_area == com.vport[GREEN_VPORT])
		return GREEN_VPORT;
	else if (drawing_area == com.vport[BLUE_VPORT])
		return BLUE_VPORT;
	else if (allow_rgb && drawing_area == com.vport[RGB_VPORT])
		return RGB_VPORT;
	return -1;
}

void calculate_fwhm(GtkWidget *widget) {
	/* calculate and display FWHM */
	int layer = match_drawing_area_widget(widget, FALSE);
	if (layer != -1) {
		char buf[64], label_name[16];
		char *layer_name = vport_number_to_name(layer);
		GtkLabel *label;
		if (com.selection.w && com.selection.h) {// Now we don't care about the size of the sample. Minimization checks that
			if (com.selection.w < 300 && com.selection.h < 300) {
				double roundness;
				double fwhm_val;

				fwhm_val = psf_get_fwhm(&gfit, layer, &roundness);
				g_snprintf(buf, sizeof(buf), "fwhm = %.2f, r = %.2f", fwhm_val,
						roundness);
			} else
				g_snprintf(buf, sizeof(buf), _("fwhm: selection is too large"));
		} else {
			g_snprintf(buf, sizeof(buf), _("fwhm: no selection"));
		}
		g_snprintf(label_name, sizeof(label_name), "labelfwhm%s", layer_name);
		free(layer_name);
		label = GTK_LABEL(gtk_builder_get_object(builder, label_name));
		gtk_label_set_text(label, buf);
	}
}

/* displays the opened image file name in the layers window.
 * if a unique file is loaded, its details are used instead of any sequence data
 */
void display_filename() {
	GtkLabel *fn_label;
	int nb_layers;
	char str[64], *filename;
	if (com.uniq) {	// unique image
		filename = com.uniq->filename;
		nb_layers = com.uniq->nb_layers;
	} else {	// sequence
		filename = malloc(256);
		seq_get_image_filename(&com.seq, com.seq.current, filename);
		nb_layers = com.seq.nb_layers;
	}
	fn_label = GTK_LABEL(gtk_builder_get_object(builder, "labelfilename_red"));
	gchar *name = g_path_get_basename(filename);
	g_snprintf(str, sizeof(str), _("%s (channel 0)"), name);
	gtk_label_set_text(fn_label, str);
	if (nb_layers == 3) {	//take in charge both sequence and single image
		fn_label = GTK_LABEL(
				gtk_builder_get_object(builder, "labelfilename_green"));
		g_snprintf(str, sizeof(str), _("%s (channel 1)"), name);
		gtk_label_set_text(fn_label, str);
		fn_label = GTK_LABEL(
				gtk_builder_get_object(builder, "labelfilename_blue"));
		g_snprintf(str, sizeof(str), _("%s (channel 2)"), name);
		gtk_label_set_text(fn_label, str);
	}
	g_free(name);
}

/* set available layers in the layer list of registration */
void set_layers_for_assign() {
	int i;
	if (!com.seq.layers)
		return;
	for (i = 0; i < com.seq.nb_layers; i++) {
		char layer[100];
		if (!com.seq.layers[i].name) {
			if (com.seq.nb_layers == 1) {
				com.seq.layers[i].name = strdup(
						predefined_layers_colors[i].name);
				com.seq.layers[i].wavelength =
						predefined_layers_colors[i].wavelength;
			} else if (com.seq.nb_layers == 3) {
				com.seq.layers[i].name = strdup(
						predefined_layers_colors[i + 1].name);
				com.seq.layers[i].wavelength =
						predefined_layers_colors[i + 1].wavelength;
			} else {
				com.seq.layers[i].name = strdup("Unassigned");
				com.seq.layers[i].wavelength = -1.0;
			}
		}
		g_snprintf(layer, sizeof(layer), "%d: %s", i, com.seq.layers[i].name);
	}
}

void set_layers_for_registration() {
	static GtkComboBoxText *cbbt_layers = NULL;
	int i;
	int reminder;

	if (cbbt_layers == NULL)
		cbbt_layers = GTK_COMBO_BOX_TEXT(
				gtk_builder_get_object(builder, "comboboxreglayer"));
	reminder = gtk_combo_box_get_active(GTK_COMBO_BOX(cbbt_layers));
	gtk_combo_box_text_remove_all(cbbt_layers);
	for (i = 0; i < com.seq.nb_layers; i++) {
		char layer[100];
		if (com.seq.layers[i].name)
			g_snprintf(layer, sizeof(layer), "%d: %s", i,
					com.seq.layers[i].name);
		else
			g_snprintf(layer, sizeof(layer), _("%d: not affected yet"), i);
		if (com.seq.regparam[i]) {
			// calculate average quality for this layer
			strcat(layer, " (*)");
		}
		gtk_combo_box_text_append_text(cbbt_layers, layer);
	}
	/* First initialization */
	if (reminder == -1) {
		if (com.seq.nb_layers == 3)
			gtk_combo_box_set_active(GTK_COMBO_BOX(cbbt_layers), 1);
		else
			gtk_combo_box_set_active(GTK_COMBO_BOX(cbbt_layers), 0);
	}
	/* Already initialized */
	else
		gtk_combo_box_set_active(GTK_COMBO_BOX(cbbt_layers), reminder);
}

void display_image_number(int index) {
	static GtkSpinButton *spin = NULL;
	if (!spin)
		spin = GTK_SPIN_BUTTON(
				gtk_builder_get_object(builder, "imagenumber_spin"));
	char text[16];
	char format[10];
	if (com.seq.fixed <= 1)
		g_snprintf(format, sizeof(format), "%%d");
	else
		g_snprintf(format, sizeof(format), "%%.%dd", com.seq.fixed);
	g_snprintf(text, sizeof(text), format, com.seq.imgparam[index].filenum);
	gtk_entry_set_text(GTK_ENTRY(spin), text);
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

/* sets text in the label and displays the dialog window 1 */
void show_dialog(const char *text, const char *title, const char *icon) {
	struct _dialog_data *args = malloc(sizeof(struct _dialog_data));
	args->text = text;
	args->title = title;
	args->icon = icon;
	gdk_threads_add_idle(show_dialog_idle, args);
}

void show_data_dialog(char *text, char *title) {
	GtkTextView *tv = GTK_TEXT_VIEW(lookup_widget("data_txt"));
	GtkTextBuffer *tbuf = gtk_text_view_get_buffer(tv);
	GtkTextIter itDebut;
	GtkTextIter itFin;

	gtk_text_buffer_get_bounds(tbuf, &itDebut, &itFin);
	gtk_text_buffer_delete(tbuf, &itDebut, &itFin);
	gtk_text_buffer_set_text(tbuf, text, strlen(text));
	gtk_window_set_title(GTK_WINDOW(lookup_widget("data_dialog")), title);

	gtk_widget_show_all(lookup_widget("data_dialog"));
}

void show_main_gray_window() {
	GtkCheckMenuItem *graycheck = GTK_CHECK_MENU_ITEM(
			gtk_builder_get_object(builder, "menuitemgray"));
	gtk_check_menu_item_set_active(graycheck, TRUE);
	gtk_widget_show_all(lookup_widget("main_window"));
	gtk_window_present(GTK_WINDOW(lookup_widget("main_window")));
}

void show_rgb_window() {
	GtkCheckMenuItem *rgbcheck = GTK_CHECK_MENU_ITEM(
			gtk_builder_get_object(builder, "menuitemcolor"));
	gtk_check_menu_item_set_active(rgbcheck, TRUE);
	gtk_widget_show_all(lookup_widget("rgb_window"));
}

void hide_rgb_window() {
	/* unchecking the menu item is done in the window destruction callback */
	/*GtkCheckMenuItem *rgbcheck =
	 GTK_CHECK_MENU_ITEM(gtk_builder_get_object(builder, "menuitemcolor"));
	 gtk_check_menu_item_set_active(rgbcheck, FALSE);*/
	gtk_widget_hide(lookup_widget("rgb_window"));
}

void set_cursor_waiting(gboolean waiting) {
	GdkCursor *cursor;
	static GdkCursor *clock = NULL;
	GdkDisplay *display;
	GdkScreen *screen;
	GList *list;

	display = gdk_display_get_default ();

	if (clock == NULL)
		clock = gdk_cursor_new_for_display(display, GDK_WATCH);

	screen = gdk_screen_get_default();
	list = gdk_screen_get_toplevel_windows(screen);

	if (waiting) {
		cursor = clock;
	} else {
		cursor = NULL;
	}
	while (list) {
		GdkWindow *window = GDK_WINDOW(list->data);
		gdk_window_set_cursor(window, cursor);
		gdk_display_sync(gdk_window_get_display(window));
		gdk_flush();
		list = g_list_next(list);
	}
	g_free(list);
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

void zoomcombo_update_display_for_zoom() {
	static GtkComboBox *zoomcombo = NULL;
	static double indexes[] = { 16., 8., 4., 2., 1., .5, .25, .125, /*.0625, */
	-1. };
	int i;
	char *msg;

	if (zoomcombo == NULL)
		zoomcombo = GTK_COMBO_BOX(gtk_builder_get_object(builder, "combozoom"));
	for (i = 0; i < sizeof(indexes) / sizeof(double); i++) {
		if (indexes[i] == com.zoom_value) {
			g_signal_handlers_block_by_func(zoomcombo, on_combozoom_changed,
					NULL);
			gtk_combo_box_set_active(zoomcombo, i);
			g_signal_handlers_unblock_by_func(zoomcombo, on_combozoom_changed,
					NULL);
			return;
		}
	}
	msg = siril_log_message(
			_("Unknown zoom_value value, what is the current zoom?\n"));
	show_dialog(msg, _("Error"), "gtk-dialog-error");
}

void initialize_FITS_name_entries() {
	GtkEntry *moffset, *mdark, *mflat, *final_stack;
	GString *str[4];
	gchar *txt[4];
	gint i;

	moffset = GTK_ENTRY(lookup_widget("offsetname_entry"));
	mdark = GTK_ENTRY(lookup_widget("darkname_entry"));
	mflat = GTK_ENTRY(lookup_widget("flatname_entry"));
	final_stack = GTK_ENTRY(lookup_widget("entryresultfile"));

	str[0] = g_string_new("master-offset");
	str[1] = g_string_new("master-dark");
	str[2] = g_string_new("master-flat");
	str[3] = g_string_new("stack_result");

	for (i = 0; i < 4; i++) {
		g_string_append(str[i], com.ext);
		txt[i] = g_string_free(str[i], FALSE);
	}
	gtk_entry_set_text(moffset, txt[0]);
	gtk_entry_set_text(mdark, txt[1]);
	gtk_entry_set_text(mflat, txt[2]);
	gtk_entry_set_text(final_stack, txt[3]);

	for (i = 0; i < 4; i++)
		g_free(txt[i]);
}

void adjust_vport_size_to_image() {
	int vport;
	// make GtkDrawingArea the same size than the image
	// http://developer.gnome.org/gtk3/3.4/GtkWidget.html#gtk-widget-set-size-request
	double zoom = get_zoom_val();
	int w, h;
	if (zoom <= 0)
		return;
	w = (int) (((double) gfit.rx) * zoom);
	h = (int) (((double) gfit.ry) * zoom);
	for (vport = 0; vport < MAXVPORT; vport++)
		gtk_widget_set_size_request(com.vport[vport], w, h);
	fprintf(stdout, "set new vport size (%d, %d)\n", w, h);
}

/* when a sequence is loaded, the processing (stacking) output file name is
 * modified to include the name of the sequence */
void set_output_filename_to_sequence_name() {
	static GtkEntry *output_file = NULL;
	gchar msg[256];
	if (!output_file)
		output_file = GTK_ENTRY(
				gtk_builder_get_object(builder, "entryresultfile"));
	if (!com.seq.seqname || *com.seq.seqname == '\0')
		return;
	g_snprintf(msg, sizeof(msg), "%s%sstacked%s", com.seq.seqname,
			ends_with(com.seq.seqname, "_") ?
					"" : (ends_with(com.seq.seqname, "-") ? "" : "_"), com.ext);
	gtk_entry_set_text(output_file, msg);
}

void adjust_refimage(int n) {
	static GtkWidget *ref_butt = NULL;
	if (ref_butt == NULL)
		ref_butt = lookup_widget("refframe");

	//fprintf(stdout, "adjust refimage: %d (ref is %d)\n", n, com.seq.reference_image);
	g_signal_handlers_block_by_func(ref_butt, on_ref_frame_toggled, NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ref_butt), com.seq.reference_image == n);
	g_signal_handlers_unblock_by_func(ref_butt, on_ref_frame_toggled, NULL);
}

void close_tab() {
	GtkNotebook* Color_Layers = GTK_NOTEBOOK(
			gtk_builder_get_object(builder, "notebook1"));
	GtkWidget* page;

	if (com.seq.nb_layers == 1 || gfit.naxes[2] == 1) {
		page = gtk_notebook_get_nth_page(Color_Layers, GREEN_VPORT);
		gtk_widget_hide(page);
		page = gtk_notebook_get_nth_page(Color_Layers, BLUE_VPORT);
		gtk_widget_hide(page);
		page = gtk_notebook_get_nth_page(Color_Layers, RED_VPORT);
		gtk_notebook_set_tab_label_text(Color_Layers, page, _("B&W channel"));
	} else {
		page = gtk_notebook_get_nth_page(Color_Layers, RED_VPORT);
		gtk_notebook_set_tab_label_text(Color_Layers, page, _("Red channel"));
		page = gtk_notebook_get_nth_page(Color_Layers, GREEN_VPORT);
		gtk_widget_show(page);
		page = gtk_notebook_get_nth_page(Color_Layers, BLUE_VPORT);
		gtk_widget_show(page);
	}
}

void activate_tab(int vport) {
	GtkNotebook* notebook = GTK_NOTEBOOK(
			gtk_builder_get_object(builder, "notebook1"));
	gtk_notebook_set_current_page(notebook, vport);
	// com.cvport is set in the event handler for changed page
}

void control_window_switch_to_tab(main_tabs tab) {
	GtkNotebook* notebook = GTK_NOTEBOOK(
			gtk_builder_get_object(builder, "notebook2"));
	gtk_notebook_set_current_page(notebook, tab);
}

void update_statusbar_convert() {
	GtkLabel *status_label = GTK_LABEL(
			gtk_builder_get_object(builder, "statuslabel_convert"));

	int nb_files = count_selected_files();
	if (nb_files == 0)
		gtk_label_set_text(status_label, " ");
	else {
		char str[64];
		g_snprintf(str, sizeof(str), _("%d files loaded"), nb_files);
		gtk_label_set_text(status_label, str);
	}
}

void update_spinCPU(int max) {
	static GtkSpinButton * spin_cpu = NULL;

	if (spin_cpu == NULL) {
		spin_cpu = GTK_SPIN_BUTTON(lookup_widget("spinCPU"));
	}
	if (max > 0) {
		gtk_spin_button_set_range (spin_cpu, 1, (gdouble) max);
	}
	gtk_spin_button_set_value (spin_cpu, (gdouble) com.max_thread);
}

/*****************************************************************************
 *             I N I T I A L I S A T I O N      F U N C T I O N S            *
 ****************************************************************************/

void initialize_shortcuts() {
	/* activate accelerators (keyboard shortcut in GTK language) */
	static GtkAccelGroup *accel = NULL;

	if (accel == NULL) {
		accel = GTK_ACCEL_GROUP(gtk_builder_get_object(builder, "accelgroup1"));
	}
	/* EXIT */
	gtk_widget_add_accelerator(lookup_widget("exit"), "activate", accel,
	GDK_KEY_q, get_default_modifier(), GTK_ACCEL_VISIBLE);
	/* UNDO */
	gtk_widget_add_accelerator(lookup_widget("undo_item"), "activate", accel,
	GDK_KEY_z, get_default_modifier(), GTK_ACCEL_VISIBLE);
	gtk_widget_add_accelerator(lookup_widget("undo_item1"), "activate", accel,
	GDK_KEY_z, get_default_modifier(), GTK_ACCEL_VISIBLE);
	/* REDO */
	gtk_widget_add_accelerator(lookup_widget("redo_item"), "activate", accel,
	GDK_KEY_z, get_default_modifier() | GDK_SHIFT_MASK, GTK_ACCEL_VISIBLE);
	gtk_widget_add_accelerator(lookup_widget("redo_item1"), "activate", accel,
	GDK_KEY_z, get_default_modifier() | GDK_SHIFT_MASK, GTK_ACCEL_VISIBLE);
	/* OPEN */
	gtk_widget_add_accelerator(lookup_widget("open1"), "activate", accel,
	GDK_KEY_o, get_default_modifier(), GTK_ACCEL_VISIBLE);
	/* SAVE */
	gtk_widget_add_accelerator(lookup_widget("menu_save_fits"), "activate", accel,
	GDK_KEY_s, get_default_modifier(), GTK_ACCEL_VISIBLE);
	gtk_widget_add_accelerator(lookup_widget("menu_rgb_savefits"), "activate", accel,
	GDK_KEY_s, get_default_modifier(), GTK_ACCEL_VISIBLE);
	gtk_widget_add_accelerator(lookup_widget("menu_save_tiff"), "activate", accel,
	GDK_KEY_t, get_default_modifier(), GTK_ACCEL_VISIBLE);
	gtk_widget_add_accelerator(lookup_widget("menu_rgb_savetiff"), "activate", accel,
	GDK_KEY_t, get_default_modifier(), GTK_ACCEL_VISIBLE);
	gtk_widget_add_accelerator(lookup_widget("menu_save_bmp"), "activate", accel,
	GDK_KEY_b, get_default_modifier(), GTK_ACCEL_VISIBLE);
	gtk_widget_add_accelerator(lookup_widget("menu_rgb_savebmp"), "activate", accel,
	GDK_KEY_b, get_default_modifier(), GTK_ACCEL_VISIBLE);
	gtk_widget_add_accelerator(lookup_widget("menu_save_jpg"), "activate", accel,
	GDK_KEY_j, get_default_modifier(), GTK_ACCEL_VISIBLE);
	gtk_widget_add_accelerator(lookup_widget("menu_rgb_savejpg"), "activate", accel,
	GDK_KEY_j, get_default_modifier(), GTK_ACCEL_VISIBLE);
	gtk_widget_add_accelerator(lookup_widget("menu_save_pbm"), "activate", accel,
	GDK_KEY_p, get_default_modifier(), GTK_ACCEL_VISIBLE);
	gtk_widget_add_accelerator(lookup_widget("menu_rgb_save8ppm"), "activate", accel,
	GDK_KEY_p, get_default_modifier(), GTK_ACCEL_VISIBLE);
}

void initialize_remap() {
	int i;
	for (i = 0; i < MAXGRAYVPORT; i++) {
		remap_index[i] = NULL;
		last_pente[i] = 0.f;
		last_mode[i] = HISTEQ_DISPLAY;
		// only HISTEQ mode always computes the index, it's a good initializer here
	}
}

/* Initialize the combobox when loading new single_image */
void initialize_display_mode() {
	static GtkComboBox *modecombo = NULL;
	static GtkToggleButton *chainedbutton = NULL;
	display_mode mode;
	int i, raw_mode;

	if (!modecombo)
		modecombo = GTK_COMBO_BOX(lookup_widget("combodisplay"));
	raw_mode = gtk_combo_box_get_active(modecombo);
	/* Check if never initialized. In this case the mode is set to linear */
	if (raw_mode == -1)
		mode = NORMAL_DISPLAY;
	else
		mode = raw_mode;
	/* The mode is applyed for each layer */
	if (single_image_is_loaded() && com.cvport < com.uniq->nb_layers
	&& com.seq.current != RESULT_IMAGE) {
		for (i = 0; i < com.uniq->nb_layers; i++)
			com.uniq->layers[i].rendering_mode = mode;
	} else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers) {
		for (i = 0; i < com.seq.nb_layers; i++)
			com.seq.layers[i].rendering_mode = mode;
	}
	/* In the case where the layer were unchained, we chaine it */
	if (!chainedbutton)
		chainedbutton = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_chain"));
	if (!gtk_toggle_button_get_active(chainedbutton)) {
		g_signal_handlers_block_by_func(chainedbutton, on_checkchain_toggled,
				NULL);
		gtk_toggle_button_set_active(chainedbutton, TRUE);
		g_signal_handlers_unblock_by_func(chainedbutton, on_checkchain_toggled,
				NULL);
	}
}

void set_GUI_CWD() {
	if (!com.wd)
		return;
	GtkLabel *label = GTK_LABEL(lookup_widget("labelcwd"));
	gtk_label_set_text(label, com.wd);
}

void set_GUI_misc() {
	GtkToggleButton *ToggleButton;

	ToggleButton = GTK_TOGGLE_BUTTON(lookup_widget("miscAskQuit"));
	gtk_toggle_button_set_active(ToggleButton, com.dontShowConfirm);
	ToggleButton = GTK_TOGGLE_BUTTON(lookup_widget("darkThemeCheck"));
	gtk_toggle_button_set_active(ToggleButton, com.have_dark_theme);
}

/* size is in kiB */
void set_GUI_MEM(unsigned long size) {
	char str[20];
	if (size != 0)
		g_snprintf(str, sizeof(str), _("Mem: %ldMB"), size / 1024);
	else
		g_snprintf(str, sizeof(str), _("Mem: N/A"));
	set_label_text_from_main_thread("labelmem", str);
}

void initialize_preprocessing() {
	GtkToggleButton *ToggleButton;

	ToggleButton = GTK_TOGGLE_BUTTON(lookup_widget("cosmCFACheck"));
	gtk_toggle_button_set_active(ToggleButton, com.prepro_cfa);
}

void set_libraw_settings_menu_available(gboolean activate) {
	GtkNotebook *notebook = GTK_NOTEBOOK(lookup_widget("notebook3"));
	GtkWidget *widget = gtk_notebook_get_nth_page (notebook, 0);

	gtk_widget_set_visible(widget, activate);
}

void set_GUI_CAMERA() {
	GtkComboBox *binning = GTK_COMBO_BOX(
			gtk_builder_get_object(builder, "combobinning"));

	if (gfit.focal_length) {
		char focal[8];
		g_snprintf(focal, sizeof(focal), "%g", gfit.focal_length);
		gtk_entry_set_text(GTK_ENTRY(lookup_widget("focal_entry")), focal);
	}
	if (gfit.pixel_size_x) {
		char pitchX[8];
		g_snprintf(pitchX, sizeof(pitchX), "%g", gfit.pixel_size_x);
		gtk_entry_set_text(GTK_ENTRY(lookup_widget("pitchX_entry")), pitchX);
	}
	if (gfit.pixel_size_y) {
		char pitchY[8];
		g_snprintf(pitchY, sizeof(pitchY), "%g", gfit.pixel_size_y);
		gtk_entry_set_text(GTK_ENTRY(lookup_widget("pitchY_entry")), pitchY);
	}

	if (!gfit.binning_x || !gfit.binning_y) {
		gtk_combo_box_set_active(binning, 0);
	}
	/* squared binning */
	else if (gfit.binning_x == gfit.binning_y)
		gtk_combo_box_set_active(binning, (gint) gfit.binning_x - 1);
	else {
		short coeff =
				gfit.binning_x > gfit.binning_y ?
						gfit.binning_x / gfit.binning_y :
						gfit.binning_y / gfit.binning_x;
		switch (coeff) {
		case 2:
			gtk_combo_box_set_active(binning, 4);
			break;
		case 3:
			gtk_combo_box_set_active(binning, 5);
			break;
		default:
			siril_log_message(_("This binning is not handled yet\n"));
		}
	}
}

void set_GUI_LIBRAW() {

	/**********COLOR ADJUSTEMENT**************/
	gtk_spin_button_set_value(
			GTK_SPIN_BUTTON(lookup_widget("Brightness_spinbutton")),
			com.raw_set.bright);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("Red_spinbutton")),
			com.raw_set.mul[0]);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lookup_widget("Blue_spinbutton")),
			com.raw_set.mul[2]);

	gtk_toggle_button_set_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_multipliers")),
			com.raw_set.auto_mul);
	gtk_toggle_button_set_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_blackpoint")),
			com.raw_set.user_black);

	/**************WHITE BALANCE**************/
	if (com.raw_set.use_camera_wb) {
		gtk_toggle_button_set_active(
				GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_cam")),
				com.raw_set.use_camera_wb);
	}

	if (com.raw_set.use_auto_wb) {
		gtk_toggle_button_set_active(
				GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_auto")),
				com.raw_set.use_auto_wb);
	}

	/********MATRIX INTERPOLATION**************/
	gtk_combo_box_set_active(GTK_COMBO_BOX(lookup_widget("combo_dcraw_inter")),
			com.raw_set.user_qual);

	/********GAMMA CORRECTION**************/
	if (com.raw_set.gamm[0] == 1.0 && com.raw_set.gamm[1] == 1.0)
		gtk_toggle_button_set_active(
				GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm0")), TRUE);
	else if (com.raw_set.gamm[0] == 2.222 && com.raw_set.gamm[1] == 4.5)
		gtk_toggle_button_set_active(
				GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm1")), TRUE);
	else
		gtk_toggle_button_set_active(
				GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_gamm2")), TRUE);

	/********** DEBAYER ******************/
	GtkComboBox *pattern = GTK_COMBO_BOX(lookup_widget("comboBayer_pattern"));
	GtkComboBox *inter = GTK_COMBO_BOX(lookup_widget("comboBayer_inter"));
	GtkToggleButton *compat = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_debayer_compatibility"));
	GtkToggleButton *use_header = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_SER_use_header"));
	GtkToggleButton *demosaicingButton = GTK_TOGGLE_BUTTON(lookup_widget("demosaicingButton"));
	gtk_combo_box_set_active(pattern, com.debayer.bayer_pattern);
	gtk_combo_box_set_active(inter, com.debayer.bayer_inter);
	gtk_toggle_button_set_active(compat, com.debayer.compatibility);
	gtk_toggle_button_set_active(use_header, com.debayer.use_bayer_header);
	gtk_toggle_button_set_active(demosaicingButton,	com.debayer.open_debayer);
}

/*****************************************************************************
 *      P U B L I C      C A L L B A C K      F U N C T I O N S              *
 ****************************************************************************/

void on_register_all_toggle(GtkToggleButton *togglebutton, gpointer user_data) {
	update_reg_interface(TRUE);
}

/* callback for GtkDrawingArea, draw event
 * see http://developer.gnome.org/gtk3/3.2/GtkDrawingArea.html
 * http://developer.gnome.org/gdk-pixbuf/stable/gdk-pixbuf-Image-Data-in-Memory.html
 * http://www.cairographics.org/manual/
 * http://www.cairographics.org/manual/cairo-Image-Surfaces.html#cairo-image-surface-create-for-data
 */
gboolean redraw_drawingarea(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int image_width, image_height, window_width, window_height;
	int vport, i = 0;
	double zoom;

	// we need to identify which vport is being redrawn
	vport = match_drawing_area_widget(widget, TRUE);
	if (vport == -1) {
		fprintf(stderr, "Could not find the vport for the draw callback\n");
		return TRUE;
	}

	window_width = gtk_widget_get_allocated_width(widget);
	window_height = gtk_widget_get_allocated_height(widget);
	zoom = get_zoom_val();
	image_width = (int) (((double) window_width) / zoom);
	image_height = (int) (((double) window_height) / zoom);

	/* draw the RGB and gray images */
	if (vport == RGB_VPORT) {
		if (com.rgbbuf) {
			cairo_scale(cr, zoom, zoom);
			cairo_set_source_surface(cr, com.surface[RGB_VPORT], 0, 0);
			cairo_paint(cr);
		} else {
			fprintf(stdout, "RGB buffer is empty, drawing black image\n");
			draw_empty_image(cr, window_width, window_height);
		}
	} else {
		if (com.graybuf[vport]) {
			cairo_scale(cr, zoom, zoom);
			cairo_set_source_surface(cr, com.surface[vport], 0, 0);
			cairo_paint(cr);
		} else {
			fprintf(stdout, "Buffer %d is empty, drawing black image\n", vport);
			draw_empty_image(cr, window_width, window_height);
		}
	}

	/* draw the selection rectangle */
	if (com.selection.w > 0 && com.selection.h > 0) {
		static double dash_format[] = { 4.0, 2.0 };
		/* fprintf(stdout, "drawing the selection rectangle (%d,%d) (%d,%d)\n",
		 com.selection.x, com.selection.y, com.selection.w, com.selection.h); */
		cairo_set_line_width(cr, 0.8 / zoom);
		cairo_set_dash(cr, dash_format, 2, 0);
		cairo_set_source_rgb(cr, 0.8, 1.0, 0.8);
		cairo_rectangle(cr, (double) com.selection.x, (double) com.selection.y,
				(double) com.selection.w, (double) com.selection.h);
		cairo_stroke(cr);
	}

	/* draw detected stars and highlight the selected star */
	if (com.stars) {
		/* com.stars is a NULL-terminated array */
		cairo_set_dash(cr, NULL, 0, 0);
		cairo_set_source_rgba(cr, 1.0, 0.4, 0.0, 0.9);
		cairo_set_line_width(cr, 1.5);

		while (com.stars[i]) {
			// by design Sx>Sy, we redefine FWHM to be sure to have the value in px
			double size = sqrt(com.stars[i]->fwhmx / 2.) * 2 * sqrt(log(2.) * 3);

			if (i == com.selected_star) {
				cairo_set_line_width(cr, 2);
				cairo_set_source_rgba(cr, 0.0, 0.4, 1.0, 0.6);
				cairo_rectangle(cr, com.stars[i]->xpos - 1.5 * size,
						com.stars[i]->ypos - 1.5 * size, 3 * size, 3 * size);
				cairo_stroke(cr);

				cairo_set_line_width(cr, 1.5/zoom);
				cairo_set_source_rgba(cr, 1.0, 0.4, 0.0, 0.9);
			}
			cairo_arc(cr, com.stars[i]->xpos, com.stars[i]->ypos, size, 0., 2. * M_PI);
			cairo_stroke(cr);
			i++;
		}
	}

	if (sequence_is_loaded()) {
		/* draw seqpsf stars */
		for (i = 0; i < MAX_SEQPSF && com.seq.photometry[i]; i++) {
			cairo_set_dash(cr, NULL, 0, 0);
			cairo_set_source_rgba(cr, com.seq.photometry_colors[i][0],
					com.seq.photometry_colors[i][1],
					com.seq.photometry_colors[i][2], 1.0);
			cairo_set_line_width(cr, 2.0/zoom);
			fitted_PSF *the_psf = com.seq.photometry[i][com.seq.current];
			if (the_psf) {
				double size = sqrt(the_psf->fwhmx / 2.) * 2 * sqrt(log(2.) * 4);
				cairo_arc(cr, the_psf->xpos, the_psf->ypos, size, 0., 2. * M_PI);
				cairo_stroke(cr);
			}
		}

		/* draw a cross on excluded images */
		if (com.seq.imgparam && com.seq.current >= 0 &&
				!com.seq.imgparam[com.seq.current].incl) {
			if (image_width > gfit.rx)
				image_width = gfit.rx;
			if (image_height > gfit.ry)
				image_height = gfit.ry;
			cairo_set_dash(cr, NULL, 0, 0);
			cairo_set_source_rgb(cr, 1.0, 0.8, 0.7);
			cairo_set_line_width(cr, 2.0 / zoom);
			cairo_move_to(cr, 0, 0);
			cairo_line_to(cr, image_width, image_height);
			cairo_move_to(cr, 0, image_height);
			cairo_line_to(cr, image_width, 0);
			cairo_stroke(cr);
		}

		/* draw preview rectangles for the manual registration */
		for (i = 0; i < PREVIEW_NB; i++) {
			if (com.seq.previewX[i] >= 0) {
				int textX, textY;
				char text[3];
				cairo_set_line_width(cr, 0.5 / zoom);
				cairo_set_source_rgb(cr, 0.1, 0.6, 0.0);
				cairo_rectangle(cr,
						com.seq.previewX[i] - com.seq.previewW[i] / 2,
						com.seq.previewY[i] - com.seq.previewH[i] / 2,
						com.seq.previewW[i], com.seq.previewH[i]);
				cairo_stroke(cr);

				textX = com.seq.previewX[i] - com.seq.previewW[i] / 2;
				if (textX < 0)
					textX += com.seq.previewW[i] - 20;
				else
					textX += 15;
				textY = com.seq.previewY[i] - com.seq.previewH[i] / 2;
				if (textY < 0)
					textY += com.seq.previewH[i] - 15;
				else
					textY += 20;
				g_snprintf(text, sizeof(text), "%d", i + 1);

				cairo_set_font_size(cr, 12.0 / zoom);
				cairo_move_to(cr, textX, textY);
				cairo_show_text(cr, text);
			}
		}
	}

	/* draw background removal gradient selection boxes */
	if (com.grad && com.grad_boxes_drawn) {
		int i = 0;
		while (i < com.grad_nb_boxes) {
			if (com.grad[i].boxvalue[0] != -1.0) {
				cairo_set_line_width(cr, 1.5);
				cairo_set_source_rgba(cr, 0.2, 1.0, 0.3, 1.0);
				cairo_rectangle(cr, com.grad[i].centre.x - com.grad_size_boxes,
						com.grad[i].centre.y - com.grad_size_boxes,
						com.grad_size_boxes, com.grad_size_boxes);
				cairo_stroke(cr);
			}
			i++;
		}
	}
	return FALSE;
}

/* when the cursor moves, update the value displayed in the textbox and save it
 * in the related layer_info. Does not change display until cursor is released. */
void on_minscale_changed(GtkRange *range, gpointer user_data) {
	static GtkEntry *minentry = NULL;
	char buffer[10];
	if (minentry == NULL)
		minentry = GTK_ENTRY(gtk_builder_get_object(builder, "min_entry"));

	if (single_image_is_loaded() && com.seq.current < RESULT_IMAGE &&
			com.cvport < com.uniq->nb_layers) {
		com.uniq->layers[com.cvport].lo = (int) gtk_range_get_value(range);
		g_snprintf(buffer, 6, "%u", com.uniq->layers[com.cvport].lo);
	} else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers) {
		com.seq.layers[com.cvport].lo = (int) gtk_range_get_value(range);
		g_snprintf(buffer, 6, "%u", com.seq.layers[com.cvport].lo);
	} else return;
	g_signal_handlers_block_by_func(minentry, on_min_entry_changed, NULL);
	gtk_entry_set_text(minentry, buffer);
	g_signal_handlers_unblock_by_func(minentry, on_min_entry_changed, NULL);
}

gboolean on_minscale_release(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	if (com.sliders != USER) {
		com.sliders = USER;
		sliders_mode_set_state(com.sliders);
	}
	if (copy_rendering_settings_when_chained(FALSE))
		redraw(com.cvport, REMAP_ALL);
	else
		redraw(com.cvport, REMAP_ONLY);
	redraw_previews();
	return FALSE;
}

/* when the cursor moves, update the value displayed in the textbox and save it
 * in the related layer_info. Does not change display until cursor is released. */
void on_maxscale_changed(GtkRange *range, gpointer user_data) {
	static GtkEntry *maxentry = NULL;
	char buffer[10];
	if (maxentry == NULL)
		maxentry = GTK_ENTRY(gtk_builder_get_object(builder, "max_entry"));

	if (single_image_is_loaded() && com.seq.current < RESULT_IMAGE &&
			com.cvport < com.uniq->nb_layers) {
		com.uniq->layers[com.cvport].hi = (int) gtk_range_get_value(range);
		g_snprintf(buffer, 6, "%u", com.uniq->layers[com.cvport].hi);
	} else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers) {
		com.seq.layers[com.cvport].hi = (int) gtk_range_get_value(range);
		g_snprintf(buffer, 6, "%u", com.seq.layers[com.cvport].hi);
	} else return;
	g_signal_handlers_block_by_func(maxentry, on_max_entry_changed, NULL);
	gtk_entry_set_text(maxentry, buffer);
	g_signal_handlers_unblock_by_func(maxentry, on_max_entry_changed, NULL);
}

gboolean on_maxscale_release(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {
	if (com.sliders != USER) {
		com.sliders = USER;
		sliders_mode_set_state(com.sliders);
	}
	if (copy_rendering_settings_when_chained(TRUE))
		redraw(com.cvport, REMAP_ALL);
	else
		redraw(com.cvport, REMAP_ONLY);
	redraw_previews();
	return FALSE;
}

/* a checkcut checkbox was toggled. Update the layer_info and others if chained. */
void on_checkcut_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	if (copy_rendering_settings_when_chained(TRUE))
		redraw(com.cvport, REMAP_ALL);
	else
		redraw(com.cvport, REMAP_ONLY);
	redraw_previews();
}

void on_darkfile_button_clicked(GtkButton *button, gpointer user_data) {
	whichdial = OD_DARK;
	opendial();
}

void on_cwd_btton_clicked(GtkButton *button, gpointer user_data) {
	whichdial = OD_CWD;
	opendial();
}

void on_offsetfile_button_clicked(GtkButton *button, gpointer user_data) {
	whichdial = OD_OFFSET;
	opendial();
}

void on_flatfile_button_clicked(GtkButton *button, gpointer user_data) {
	whichdial = OD_FLAT;
	opendial();
}

void on_open1_activate(GtkMenuItem *menuitem, gpointer user_data) {
	whichdial = OD_OPEN;
	opendial();
}

void on_cosmEnabledCheck_toggled(GtkToggleButton *button, gpointer user_data) {
	GtkWidget *CFA, *SigHot, *SigCold, *checkHot, *checkCold, *evaluateButton;
	gboolean is_active;

	CFA = lookup_widget("cosmCFACheck");
	SigHot = lookup_widget("spinSigCosmeHot");
	SigCold = lookup_widget("spinSigCosmeCold");
	checkHot = lookup_widget("checkSigHot");
	checkCold = lookup_widget("checkSigCold");
	evaluateButton = lookup_widget("GtkButtonEvaluateCC");

	is_active = gtk_toggle_button_get_active(button);

	gtk_widget_set_sensitive(CFA, is_active);
	gtk_widget_set_sensitive(SigHot, is_active);
	gtk_widget_set_sensitive(SigCold, is_active);
	gtk_widget_set_sensitive(checkHot, is_active);
	gtk_widget_set_sensitive(checkCold, is_active);
	gtk_widget_set_sensitive(evaluateButton, is_active);
}

void on_cosmCFACheck_toggled(GtkToggleButton *button, gpointer user_data) {
	com.prepro_cfa = gtk_toggle_button_get_active(button);
	writeinitfile();
}

void on_GtkButtonEvaluateCC_clicked(GtkButton *button, gpointer user_data) {
	GtkEntry *entry;
	GtkLabel *label[2];
	GtkWidget *widget[2];
	const char *filename;
	char *str[2];
	double sig[2];
	long icold = 0L, ihot = 0L;

	set_cursor_waiting(TRUE);
	sig[0] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeColdBox")));
	sig[1] = gtk_spin_button_get_value(GTK_SPIN_BUTTON(lookup_widget("spinSigCosmeHotBox")));
	widget[0] = lookup_widget("GtkLabelColdCC");
	widget[1] = lookup_widget("GtkLabelHotCC");
	label[0] = GTK_LABEL(lookup_widget("GtkLabelColdCC"));
	label[1] = GTK_LABEL(lookup_widget("GtkLabelHotCC"));
	entry = GTK_ENTRY(lookup_widget("darkname_entry"));
	filename = gtk_entry_get_text(entry);
	if (filename) {
		int ret = readfits(filename, &(wfit[4]), NULL);
		if (!ret) {
			count_deviant_pixels(&(wfit[4]), sig, &icold, &ihot);
		}
	}
	if (icold > 10000) {
		str[0] = g_markup_printf_escaped(_("<span foreground=\"red\">Cold: %ld px</span>"), icold);
		gtk_widget_set_tooltip_text(widget[0], _("This value may be to high. Please, consider to change sigma value or uncheck the box."));
	}
	else {
		str[0] = g_markup_printf_escaped(_("Cold: %ld px"), icold);
		gtk_widget_set_tooltip_text(widget[0], "");
	}
	gtk_label_set_markup(label[0], str[0]);

	if (ihot > 10000) {
		str[1] = g_markup_printf_escaped(_("<span foreground=\"red\">Hot: %ld px</span>"), ihot);
		gtk_widget_set_tooltip_text(widget[1], _("This value may be to high. Please, consider to change sigma value or uncheck the box."));
	}
	else {
		str[1] = g_markup_printf_escaped(_("Hot: %ld px"), ihot);
		gtk_widget_set_tooltip_text(widget[1], "");
	}
	gtk_label_set_markup(label[1], str[1]);
	g_free(str[0]);
	g_free(str[1]);
	set_cursor_waiting(FALSE);
}

void on_settings_activate(GtkMenuItem *menuitem, gpointer user_data) {
	gtk_widget_show(lookup_widget("settings_window"));
}

void on_menu_FITS_header_activate(GtkMenuItem *menuitem, gpointer user_data) {
	show_FITS_header(&gfit);
}

void on_close_settings_button_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("settings_window"));
}

void on_focal_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar* focal_entry = gtk_entry_get_text(GTK_ENTRY(editable));
	gfit.focal_length = atof(focal_entry);
	update_fwhm_units_ok();
}

void on_pitchX_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar* pitchX_entry = gtk_entry_get_text(GTK_ENTRY(editable));
	gfit.pixel_size_x = (float)atof(pitchX_entry);
	update_fwhm_units_ok();
}

void on_pitchY_entry_changed(GtkEditable *editable, gpointer user_data) {
	const gchar* pitchY_entry = gtk_entry_get_text(GTK_ENTRY(editable));
	gfit.pixel_size_y = (float)atof(pitchY_entry);
	update_fwhm_units_ok();
}

void on_button_clear_sample_clicked(GtkButton *button, gpointer user_data) {
	clear_sampling_setting_box();
}

void on_comboBayer_pattern_changed(GtkComboBox* box, gpointer user_data) {
	com.debayer.bayer_pattern = gtk_combo_box_get_active(box);
}

void on_comboBayer_inter_changed(GtkComboBox* box, gpointer user_data) {
	com.debayer.bayer_inter = gtk_combo_box_get_active(box);
}

void on_checkbutton_cam_toggled(GtkButton *button, gpointer user_data) {
	GtkToggleButton *auto_button = GTK_TOGGLE_BUTTON(
			lookup_widget("checkbutton_auto"));
	GtkToggleButton *cam_button = GTK_TOGGLE_BUTTON(
			lookup_widget("checkbutton_cam"));

	if (gtk_toggle_button_get_active(auto_button) == TRUE) {
		g_signal_handlers_block_by_func(auto_button,
				on_checkbutton_auto_toggled, NULL);
		gtk_toggle_button_set_active(auto_button, FALSE);
		g_signal_handlers_unblock_by_func(auto_button,
				on_checkbutton_auto_toggled, NULL);
		gtk_toggle_button_set_active(cam_button, TRUE);
	}
}

void on_checkbutton_auto_toggled(GtkButton *button, gpointer user_data) {
	GtkToggleButton *auto_button = GTK_TOGGLE_BUTTON(
			lookup_widget("checkbutton_auto"));
	GtkToggleButton *cam_button = GTK_TOGGLE_BUTTON(
			lookup_widget("checkbutton_cam"));

	if (gtk_toggle_button_get_active(cam_button) == TRUE) {
		g_signal_handlers_block_by_func(cam_button, on_checkbutton_cam_toggled,
				NULL);
		gtk_toggle_button_set_active(cam_button, FALSE);
		g_signal_handlers_unblock_by_func(cam_button,
				on_checkbutton_cam_toggled, NULL);
		gtk_toggle_button_set_active(auto_button, TRUE);
	}
}

void on_checkbutton_auto_evaluate_toggled(GtkToggleButton *button,
		gpointer user_data) {
	GtkWidget *entry = lookup_widget("entry_flat_norm");

	gtk_widget_set_sensitive(entry, !gtk_toggle_button_get_active(button));
}

void on_settings_window_hide(GtkWidget *widget, gpointer user_data) {
	update_libraw_interface();
}

void on_combobinning_changed(GtkComboBox *box, gpointer user_data) {
	gint index = gtk_combo_box_get_active(box);

	switch (index) {
	case 0:
	case 1:
	case 2:
	case 3:
		gfit.binning_x = gfit.binning_y = (short) index + 1;
		break;
	case 4:
		gfit.binning_x = 1;
		gfit.binning_x = 2;
		break;
	case 5:
		gfit.binning_x = 1;
		gfit.binning_y = 3;
		break;
	default:
		fprintf(stderr, "Should not happen\n");
	}
}

void on_checkbutton_multipliers_toggled(GtkToggleButton *button,
		gpointer user_data) {
	gboolean active;

	active = gtk_toggle_button_get_active(button);

	gtk_widget_set_sensitive(lookup_widget("hbox8"), !active);
	gtk_widget_set_sensitive(lookup_widget("hbox11"), !active);
	if (active) {
		gtk_spin_button_set_value(
				GTK_SPIN_BUTTON(lookup_widget("Red_spinbutton")), 1.0);
		gtk_spin_button_set_value(
				GTK_SPIN_BUTTON(lookup_widget("Blue_spinbutton")), 1.0);
	}
}

void on_filechooser_swap_file_set(GtkFileChooserButton *fileChooser, gpointer user_data) {
	GtkFileChooser *swap_dir = GTK_FILE_CHOOSER(fileChooser);
	GtkLabel *label = GTK_LABEL(lookup_widget("label_swap_dir"));
	char *dir;

	dir = gtk_file_chooser_get_filename (swap_dir);

	if (com.swap_dir)
		free(com.swap_dir);
	com.swap_dir = (char*) dir;
	gtk_label_set_text (label, com.swap_dir);
	writeinitfile();
}

void on_button_reset_swap_clicked(GtkButton *button, gpointer user_data) {
	reset_swapdir();
}

void on_spinbutton_mem_value_changed(GtkSpinButton *button, gpointer user_data) {
	gdouble mem;

	mem = gtk_spin_button_get_value(button);
	com.stack.memory_percent = mem;
	writeinitfile();
}

void on_combobox_ext_changed(GtkComboBox *box, gpointer user_data) {

	if (com.ext)
		free(com.ext);

	com.ext = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(box));
	com.len_ext = strlen(com.ext);
	writeinitfile();
	initialize_FITS_name_entries();
}

void gtk_main_quit() {
	GtkWidget *widget = lookup_widget("confirmlabel");
	GtkWidget *dontShow = lookup_widget("confirmDontShowButton");

	confirm = CD_QUIT;
	if (!com.dontShowConfirm) {
		gtk_widget_set_visible(dontShow, TRUE);
		gtk_label_set_text(GTK_LABEL(widget),
				_("Are you sure you want to quit ?"));
		gtk_widget_show(lookup_widget("confirm_dialog"));
	} else {
		undo_flush();
		exit(EXIT_SUCCESS);
	}
}

void on_exit_activate(GtkMenuItem *menuitem, gpointer user_data) {
	gtk_main_quit();
}

/* handler for the single-line console */
gboolean on_command_key_press_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {
	const gchar *text;
	int handled = 0;
	static GtkEntry *entry = NULL;
	if (!entry)
		entry = GTK_ENTRY(widget);
	GtkEditable *editable = GTK_EDITABLE(entry);
	int entrylength = 0;

	switch (event->keyval) {
	case GDK_KEY_Return:
	case GDK_KEY_KP_Enter:
		handled = 1;
		text = gtk_entry_get_text(entry);
		history_add_line(strdup(text));
		if (!(processcommand(text)))
			gtk_entry_set_text(entry, "");
		break;
	case GDK_KEY_Up:
		handled = 1;
		if (!com.cmd_history)
			break;
		if (com.cmd_hist_display > 0) {
			if (com.cmd_history[com.cmd_hist_display - 1])
				--com.cmd_hist_display;
			// display previous entry
			gtk_entry_set_text(entry, com.cmd_history[com.cmd_hist_display]);
		} else if (com.cmd_history[com.cmd_hist_size - 1]) {
			// ring back, display previous
			com.cmd_hist_display = com.cmd_hist_size - 1;
			gtk_entry_set_text(entry, com.cmd_history[com.cmd_hist_display]);
		}
		entrylength = gtk_entry_get_text_length(entry);
		gtk_editable_set_position(editable, entrylength);
		break;
	case GDK_KEY_Down:
		handled = 1;
		if (!com.cmd_history)
			break;
		if (com.cmd_hist_display == com.cmd_hist_current)
			break;
		if (com.cmd_hist_display == com.cmd_hist_size - 1) {
			if (com.cmd_hist_current == 0) {
				// ring forward, end
				gtk_entry_set_text(entry, "");
				com.cmd_hist_display++;
			} else if (com.cmd_history[0]) {
				// ring forward, display next
				com.cmd_hist_display = 0;
				gtk_entry_set_text(entry, com.cmd_history[0]);
			}
		} else {
			if (com.cmd_hist_display == com.cmd_hist_current - 1) {
				// end
				gtk_entry_set_text(entry, "");
				com.cmd_hist_display++;
			} else if (com.cmd_history[com.cmd_hist_display + 1]) {
				// display next
				gtk_entry_set_text(entry,
						com.cmd_history[++com.cmd_hist_display]);
			}
		}
		entrylength = gtk_entry_get_text_length(entry);
		gtk_editable_set_position(editable, entrylength);
		break;
	case GDK_KEY_Page_Up:
	case GDK_KEY_Page_Down:
		handled = 1;
		// go to first and last in history
		break;
	}
	return (handled == 1);
}

/* mouse callbacks */
gboolean on_drawingarea_button_press_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	if (inimage((GdkEvent *) event)) {
		/* click on RGB image */
		if (widget == com.vport[RGB_VPORT]) {
			if (event->button == 3) {
				do_popup_rgbmenu(widget, event);
				return TRUE;
			}
			return FALSE;
		}

		/* else, click on gray image */
		if (event->button == 1) {	// left click
			if (mouse_status == MOUSE_ACTION_SELECT_REG_AREA) {
				if (com.drawing) {
					com.drawing = FALSE;
				} else {
					double zoom = get_zoom_val();
					com.drawing = TRUE;
					com.startX = event->x / zoom;
					com.startY = event->y / zoom;
					com.selection.h = 0;
					com.selection.w = 0;
				}
				gtk_widget_queue_draw(widget);
			}
			else if (mouse_status == MOUSE_ACTION_DRAW_SAMPLES) {
				double zoom = get_zoom_val();
				if (!com.grad) {
					com.grad = malloc(NB_MAX_OF_SAMPLES * sizeof(gradient));
					com.grad_boxes_drawn = TRUE;
					com.grad_nb_boxes = 0;
				}
				int i = com.grad_nb_boxes;
				if (i < NB_MAX_OF_SAMPLES) {
					int layer;
					GtkSpinButton *size = GTK_SPIN_BUTTON(
							lookup_widget("spinbutton_bkg_sizebox"));
					point pt;
					int midbox;

					midbox = (size_t) gtk_spin_button_get_value(size);
					com.grad_size_boxes = midbox * 2;
					pt.x = (event->x / zoom);
					pt.y = (event->y / zoom);

					if (pt.x + midbox <= gfit.rx && pt.y + midbox <= gfit.ry
							&& pt.x - midbox >= 0 && pt.y - midbox >= 0) {
						com.grad[i].centre.x = pt.x + midbox;
						com.grad[i].centre.y = pt.y + midbox;
						for (layer = 0; layer < gfit.naxes[2]; layer++)
							com.grad[i].boxvalue[layer] = get_value_from_box(
									&gfit, pt, com.grad_size_boxes, layer);
						com.grad_nb_boxes++;
						redraw(com.cvport, REMAP_NONE);
						redraw_previews();
					}
				}
			}
		} else if (event->button == 2) {	// middle click

		}
	}
	return FALSE;
}

gboolean on_drawingarea_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	if (inimage((GdkEvent *) event)) {
		double zoom = get_zoom_val();
		gdouble zoomedX, zoomedY;
		zoomedX = event->x / zoom;
		zoomedY = event->y / zoom;
		if (event->button == 1) {
			if (com.drawing && mouse_status == MOUSE_ACTION_SELECT_REG_AREA) {
				com.drawing = FALSE;
				/* finalize selection rectangle coordinates */
				if (zoomedX > com.startX) {
					com.selection.x = com.startX;
					com.selection.w = zoomedX - com.selection.x;
				} else {
					com.selection.x = zoomedX;
					com.selection.w = com.startX - zoomedX;
				}
				if (zoomedY > com.startY) {
					com.selection.y = com.startY;
					com.selection.h = zoomedY - com.selection.y;
				} else {
					com.selection.y = zoomedY;
					com.selection.h = com.startY - zoomedY;
				}
				/* we have a new rectangular selection zone,
				 * or an unselection (empty zone) */
				new_selection_zone();

				/* calculate and display FWHM - not in event
				 * callbacks because it's in the same file and
				 * requires a special argument */
				calculate_fwhm(widget);
				com.drawn = TRUE;
			} else if (mouse_status == MOUSE_ACTION_SELECT_PREVIEW1) {
				set_preview_area(0, zoomedX, zoomedY);
				mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
				// redraw to get the position of the new preview area
				gtk_widget_queue_draw(widget);
			} else if (mouse_status == MOUSE_ACTION_SELECT_PREVIEW2) {
				set_preview_area(1, zoomedX, zoomedY);
				mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
				gtk_widget_queue_draw(widget);
			}
			is_shift_on = FALSE;
		} else if (event->button == 2) {
			com.leveldrag = FALSE;
		} else if (event->button == 3) {
			do_popup_graymenu(widget, NULL);
		}
	}
	return FALSE;
}

gboolean on_drawingarea_motion_notify_event(GtkWidget *widget,
		GdkEventMotion *event, gpointer user_data) {
	//static int delay = 5;
	gchar label[32] = "labeldensity";
	fits *fit = &(gfit);

	if (inimage((GdkEvent *) event)) {
		if (com.leveldrag) {	// with button 2 down
			/* FIXME: CODE DISABLED
			 fit->lo = (event->y * (double)fit->max[com.cvport] / 2 / (double)fit->ry);
			 fit->hi = (event->x * (double)fit->max[com.cvport] * 2 / (double)fit->rx);
			 --delay;
			 if(delay<=0) {
			 level_adjust(fit, FALSE);
			 delay = 5;
			 }*/
		} else {
			char buffer[45];
			char format[25], *format_base = "x: %%.%dd y: %%.%dd = %%.%dd";
			int coords_width = 3, val_width = 3;
			double zoom = get_zoom_val();
			gint zoomedX, zoomedY;
			zoomedX = (gint) (event->x / zoom);
			zoomedY = (gint) (event->y / zoom);
			if (fit->rx >= 1000 || fit->ry >= 1000)
				coords_width = 4;
			if (fit->hi >= 1000)
				val_width = 4;
			if (fit->hi >= 10000)
				val_width = 5;
			g_snprintf(format, sizeof(format), format_base, coords_width,
					coords_width, val_width);
			g_snprintf(buffer, sizeof(buffer), format, zoomedX, zoomedY,
					fit->pdata[com.cvport][fit->rx * (fit->ry - zoomedY - 1)
							+ zoomedX]);
			/* TODO: fix to use the new function vport_number_to_name() */
			if (widget == com.vport[RED_VPORT])
				strcat(label, "r");
			else if (widget == com.vport[GREEN_VPORT])
				strcat(label, "g");
			else if (widget == com.vport[BLUE_VPORT])
				strcat(label, "b");
			else
				return FALSE;
			gtk_label_set_text(
					GTK_LABEL(gtk_builder_get_object(builder, label)), buffer);

			if (com.drawing) {	// with button 1 down
				if (zoomedX > com.startX) {
					com.selection.x = com.startX;
					com.selection.w = zoomedX - com.selection.x;
				} else {
					com.selection.x = zoomedX;
					com.selection.w = com.startX - zoomedX;
				}

				if (zoomedY > com.startY) {
					com.selection.y = com.startY;
					if (is_shift_on)
						com.selection.h = com.selection.w;
					else
						com.selection.h = zoomedY - com.selection.y;
				} else {
					com.selection.y = zoomedY;
					if (is_shift_on)
						com.selection.h = com.selection.w;
					else
						com.selection.h = com.startY - zoomedY;
				}
				gtk_widget_queue_draw(widget);
			}
		}
	}
	return FALSE;
}

void on_drawingarea_entry_notify_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	static GdkWindow *window = NULL;
	static GdkDisplay *display = NULL;
	static GdkCursor *cross = NULL;
	if (!window) {
		window = gtk_widget_get_window(lookup_widget("main_window"));
		display = gdk_window_get_display(window);
		cross = gdk_cursor_new_for_display(display, GDK_CROSSHAIR);
	}
	gdk_window_set_cursor(window, cross);
}

void on_drawingarea_leave_notify_event(GtkWidget *widget, GdkEvent *event, gpointer user_data) {
	static GdkWindow *window = NULL;
	if (!window) {
		window = gtk_widget_get_window(lookup_widget("main_window"));
	}
	gdk_window_set_cursor(window, NULL);
}

/* We give one signal event by toggle button to fix a bug. Without this solution
 * the toggled function was called 2 times
 * one call to select new button, one call to unselect previous one */
void on_radiobutton_minmax_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	if (gtk_toggle_button_get_active(togglebutton) == TRUE) {
		com.sliders = MINMAX;
		init_layers_hi_and_lo_values(com.sliders);
		set_cutoff_sliders_values();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	}
}

void on_radiobutton_hilo_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	if (gtk_toggle_button_get_active(togglebutton) == TRUE) {
		com.sliders = MIPSLOHI;
		init_layers_hi_and_lo_values(com.sliders);
		set_cutoff_sliders_values();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	}
}

void on_radiobutton_user_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	if (gtk_toggle_button_get_active(togglebutton) == TRUE) {
		com.sliders = USER;
		init_layers_hi_and_lo_values(com.sliders);
		set_cutoff_sliders_values();
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
	}
}

void on_neg_button_clicked(GtkToolButton *button, gpointer user_data) {
	int tmp;
	static GtkToggleButton *chainedbutton = NULL;
	gboolean is_chained;

	set_cursor_waiting(TRUE);

	chainedbutton = GTK_TOGGLE_BUTTON(lookup_widget("checkbutton_chain"));
	is_chained = gtk_toggle_button_get_active(chainedbutton);

	/* swaps values of hi and lo and redraw */
	if (!is_chained) {
		if (single_image_is_loaded() && com.cvport < com.uniq->nb_layers
		&& com.seq.current != RESULT_IMAGE) {
			tmp = com.uniq->layers[com.cvport].hi;
			com.uniq->layers[com.cvport].hi = com.uniq->layers[com.cvport].lo;
			com.uniq->layers[com.cvport].lo = tmp;
		} else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers) {
			tmp = com.seq.layers[com.cvport].hi;
			com.seq.layers[com.cvport].hi = com.seq.layers[com.cvport].lo;
			com.seq.layers[com.cvport].lo = tmp;
		} else
			return;
		set_cutoff_sliders_values();
		redraw(com.cvport, REMAP_ONLY);	// sliders are only set for cvport
	} else {
		int i;
		if (single_image_is_loaded() && com.seq.current != RESULT_IMAGE) {
			for (i = 0; i < com.uniq->nb_layers; i++) {
				tmp = com.uniq->layers[i].hi;
				com.uniq->layers[i].hi = com.uniq->layers[i].lo;
				com.uniq->layers[i].lo = tmp;
			}
		} else if (sequence_is_loaded()) {
			for (i = 0; i < com.seq.nb_layers; i++) {
				tmp = com.seq.layers[i].hi;
				com.seq.layers[i].hi = com.seq.layers[i].lo;
				com.seq.layers[i].lo = tmp;
			}
		} else
			return;
		set_cutoff_sliders_values();
		redraw(com.cvport, REMAP_ALL);
	}
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void on_colormap_button_toggled(GtkToggleToolButton *togglebutton,
		gpointer user_data) {
	set_cursor_waiting(TRUE);
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

/* Callback for the display mode change */
void on_combodisplay_changed(GtkComboBox *widget, gpointer user_data) {
	if (copy_rendering_settings_when_chained(TRUE))
		redraw(com.cvport, REMAP_ALL);
	else
		redraw(com.cvport, REMAP_ONLY);
	redraw_previews();
}

void on_checkchain_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	if (copy_rendering_settings_when_chained(FALSE))
		redraw(com.cvport, REMAP_ALL);
}

void on_mirrorx_button_clicked(GtkToolButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	mirrorx(&gfit, TRUE);
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void on_mirrory_button_clicked(GtkToolButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	mirrory(&gfit, TRUE);
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void on_max_entry_changed(GtkEditable *editable, gpointer user_data) {
	int value = atoi(gtk_entry_get_text(GTK_ENTRY(editable)));

	if (com.sliders != USER) {
		com.sliders = USER;
		sliders_mode_set_state(com.sliders);
	}
	if (single_image_is_loaded() && com.cvport < com.uniq->nb_layers
	&& com.seq.current != RESULT_IMAGE)
		com.uniq->layers[com.cvport].hi = value;
	else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers)
		com.seq.layers[com.cvport].hi = value;
	else
		return;
	set_cutoff_sliders_values();
	if (copy_rendering_settings_when_chained(FALSE))
		redraw(com.cvport, REMAP_ALL);
	else
		redraw(com.cvport, REMAP_ONLY);
	redraw_previews();
}

void on_min_entry_changed(GtkEditable *editable, gpointer user_data) {
	int value = atoi(gtk_entry_get_text(GTK_ENTRY(editable)));

	if (com.sliders != USER) {
		com.sliders = USER;
		sliders_mode_set_state(com.sliders);
	}
	if (single_image_is_loaded() && com.cvport < com.uniq->nb_layers
	&& com.seq.current != RESULT_IMAGE)
		com.uniq->layers[com.cvport].lo = value;
	else if (sequence_is_loaded() && com.cvport < com.seq.nb_layers)
		com.seq.layers[com.cvport].lo = value;
	else
		return;
	set_cutoff_sliders_values();
	if (copy_rendering_settings_when_chained(FALSE))
		redraw(com.cvport, REMAP_ALL);
	else
		redraw(com.cvport, REMAP_ONLY);
	redraw_previews();
}

gboolean on_main_window_key_press_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {
	//fprintf(stdout, "main window key event\n");

	/* FIXME: there is a propagation problem somewhere, we have to do it by hand */
	return on_drawingarea_key_press_event(widget, event, user_data);
}

static const gchar* copyright = N_("Copyright © 2004-2011 François Meyer\n"
		"Copyright © 2012-2016 team free-astro");

static gchar **authors = (gchar *[] ) { "Vincent Hourdin <vh@free-astro.vinvin.tf>",
				"Cyril Richard <cyril@free-astro.org>", "François Meyer", NULL };

static gchar **documenters = (gchar *[] ) { "Laurent Roge <siril.doc@orange.fr>", NULL };

static gchar **artists = (gchar *[] ) { "Coralie Monnier",
				"Cyril Richard <cyril@free-astro.org>", NULL };

// translator names
static gchar *translator = N_("Cyril Richard <cyril@free-astro.org>\n"
		"Vincent Hourdin <vh@free-astro.vinvin.tf>");

void on_about_activate(GtkMenuItem *menuitem, gpointer user_data) {
	GdkPixbuf *icon;
	GtkWindow *parent;

	parent = GTK_WINDOW(lookup_widget("control_window"));
	icon = gtk_image_get_pixbuf(GTK_IMAGE(lookup_widget("pixmap1")));
	gtk_show_about_dialog(parent,
			"program-name", PACKAGE,
			"title", _("About siril"),
			"logo", icon,
			"version", VERSION,
			"copyright", _(copyright),
			"authors", authors,
			"documenters", documenters,
			"artists", artists,
			"comments", _("Astronomical image (pre-)processing program"),
			"translator-credits", _(translator),
			"website", "https://free-astro.org/index.php/Siril",
			"website-label", _("Visit the Siril website"),
			"license-type", GTK_LICENSE_GPL_3_0,
			NULL);
}

void on_excludebutton_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	if (!com.seq.imgparam) {
		return;
	}
	toggle_image_selection(com.seq.current);
}

void on_layer_assign_selected(GtkComboBox *widget, gpointer user_data) {
	static GtkEntry *entry_name = NULL, *entry_wl = NULL;
	static GtkComboBox *cbbt_colors = NULL, *cbbt_sel = NULL;
	int layer, predef_idx;
	char wl[30];
	if (entry_name == NULL) {
		entry_name = GTK_ENTRY(gtk_builder_get_object(builder, "entrycolor"));
		entry_wl = GTK_ENTRY(gtk_builder_get_object(builder, "entrywavelen"));
		cbbt_colors = GTK_COMBO_BOX(
				gtk_builder_get_object(builder, "cbbt_colors"));
		cbbt_sel = widget;
	}
	layer = gtk_combo_box_get_active(cbbt_sel);
	if (layer == -1)
		return;

	predef_idx = get_index_in_predefined_colors_for_wavelength(
			com.seq.layers[layer].wavelength);
	if (predef_idx >= 0)
		gtk_combo_box_set_active(cbbt_colors, predef_idx);

	gtk_entry_set_text(entry_name, com.seq.layers[layer].name);
	if (com.seq.layers[layer].wavelength > 0.0)
		g_snprintf(wl, sizeof(wl), _("%g nm"), com.seq.layers[layer].wavelength);
	else
		g_snprintf(wl, sizeof(wl), _("undefined"));
	gtk_entry_set_text(entry_wl, wl);
}

/* Returns :	TRUE if the value has been displayed */
gboolean on_imagenumberspin_output(GtkSpinButton *spin, gpointer user_data) {
	static GtkAdjustment *adjustment = NULL;
	int index, do_display;
	if (adjustment == NULL)
		adjustment = gtk_spin_button_get_adjustment(spin);
	index = (int) gtk_adjustment_get_value(adjustment);

	//fprintf(stdout, "spinchanged output: index=%d\n", index);
	if (!sequence_is_loaded()) {
		// no sequence has been loaded yet
		//fprintf(stderr, "No sequence loaded\n");
		return FALSE;
	}
	if (index > com.seq.number || com.seq.current == index) {// already done for this one
	//fprintf(stderr, "index already processed or data not init. Current: %d, idx: %d\n",
	//		com.seq.current, index);
		return TRUE;
	}
	//fprintf(stdout, "SPINCHANGED: index=%d\n", index);

	do_display = (com.seq.imgparam[index].incl || com.show_excluded);
	return !seq_load_image(&com.seq, index, &gfit, do_display);
}

/* for the spin button to be able to display number which are not the actual value of
 * the adjustment, the output callback is used to modify the way they are displayed,
 * but the input callback is also required to do the opposite operation, i.e. convert
 * real image number into the number in the sequence, which is the value of the adj. */
gboolean on_imagenumberspin_input(GtkSpinButton *spin, gdouble *new_val,
		gpointer p) {
	const char *imgname = gtk_entry_get_text(GTK_ENTRY(spin));
	int imgname_int = atoi(imgname);
	int i;
	if (!sequence_is_loaded()) {
		//fprintf(stderr, "No sequence loaded\n");
		return FALSE;
	}
	i = com.seq.current;
	//fprintf(stdout, "current index in input is %d (looking for %d)\n", i, imgname_int);
	if (com.seq.imgparam[i].filenum == imgname_int) {
		*new_val = (gdouble) i;
		//fprintf(stdout, "found at 0 %d\n", i);
		return TRUE;
	} else if (i > 0 && com.seq.imgparam[i - 1].filenum == imgname_int) {
		*new_val = (gdouble) (i - 1);
		//fprintf(stdout, "found at -1 %d\n", i-1);
		return TRUE;
	} else if (i < com.seq.number - 1
			&& com.seq.imgparam[i + 1].filenum == imgname_int) {
		*new_val = (gdouble) (i + 1);
		//fprintf(stdout, "found at +1 %d\n", i+1);
		return TRUE;
	} else {
		/* no luck with neighbours, sweep it all */
		for (i = 0; i < com.seq.number; i++) {
			if (com.seq.imgparam[i].filenum == imgname_int) {
				*new_val = (gdouble) i;
				//fprintf(stdout, "sweep found at %d\n", i);
				return TRUE;
			}
		}
	}
	return GTK_INPUT_ERROR;
}

/* Returns:	TRUE to stop other handlers from being invoked for the event.
 *		FALSE to propagate the event further. */
gboolean on_imagenumberspin_key_release_event(GtkWidget *widget,
		GdkEventKey *event, gpointer user_data) {
	static GtkAdjustment *adj = NULL;
	int n;

	if (!sequence_is_loaded()) {
		//fprintf(stderr, "No sequence loaded\n");
		return TRUE;
	}
	if (!adj)
		adj = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(widget));
	n = (int) gtk_adjustment_get_value(adj);
	if (n > com.seq.number) {
		return TRUE;
	}
	if (event->keyval == GDK_KEY_space) {
		toggle_image_selection(n);
	}
	return FALSE;
}

void on_seqexcludeall_button_clicked(GtkButton *button, gpointer user_data) {
	GtkWidget *widget = lookup_widget("confirmlabel");
	GtkWidget *dontShow = lookup_widget("confirmDontShowButton");

	confirm = CD_EXCALL;
	gtk_label_set_text(GTK_LABEL(widget),
			_("Exclude all images ?\n (this erases previous image selection\n ... and there's no undo)"));
	gtk_widget_set_visible(dontShow, FALSE);
	gtk_widget_show(lookup_widget("confirm_dialog"));
}

void on_seqselectall_button_clicked(GtkButton *button, gpointer user_data) {
	GtkWidget *widget = lookup_widget("confirmlabel");
	GtkWidget *dontShow = lookup_widget("confirmDontShowButton");
	confirm = CD_INCALL;
	gtk_label_set_text(GTK_LABEL(widget),
			_("Include all images ?\n (this erases previous image selection\n ... and there's no undo)"));
	gtk_widget_set_visible(dontShow, FALSE);
	gtk_widget_show(lookup_widget("confirm_dialog"));
}

void on_prepro_button_clicked(GtkButton *button, gpointer user_data) {
	GtkToggleButton *tbutton;
	GtkEntry *entry;
	com.preprostatus = 0;

	if (!single_image_is_loaded() && !sequence_is_loaded())
		return;

	if (!single_image_is_loaded() && get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return;
	}

	// if dark selected
	tbutton = GTK_TOGGLE_BUTTON(
			gtk_builder_get_object(builder, "usedark_button"));
	if (gtk_toggle_button_get_active(tbutton) == TRUE) {
		const char *filename;
		entry = GTK_ENTRY(gtk_builder_get_object(builder, "darkname_entry"));
		filename = gtk_entry_get_text(entry);
		if (filename[0] == '\0') {
			gtk_toggle_button_set_active(tbutton, FALSE);
		} else {
			fits *dark_fit;
			progress_bar_set_text(_("Opening dark image..."));
			dark_fit = calloc(1, sizeof(fits));
			if (readfits(filename, dark_fit, NULL)) {
				siril_log_message(_("NOT USING DARK: cannot open the file\n"));
				free(dark_fit);
				gtk_entry_set_text(entry, "");
			} else {
				if (dark_fit->naxes[2] != gfit.naxes[2]) {
					const char *msg = _("NOT USING DARK: number of channels is different");
					siril_log_message("%s\n", msg);
					progress_bar_set_text(msg);
					free(dark_fit);
					gtk_entry_set_text(entry, "");
				}
				else {
					com.preprostatus |= USE_DARK;
					if (single_image_is_loaded())
						com.uniq->dark = dark_fit;
					else
						com.seq.dark = dark_fit;
				}
			}
		}
		// if dark optimization selected
		tbutton = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "checkDarkOptimize"));
		if (gtk_toggle_button_get_active(tbutton) == TRUE) {
			com.preprostatus |= USE_OPTD;
		}

		// if cosmetic correction selected
		tbutton = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "cosmEnabledCheck"));
		if (gtk_toggle_button_get_active(tbutton) == TRUE) {
			com.preprostatus |= USE_COSME;
		}
	}

	// if flat selected
	tbutton = GTK_TOGGLE_BUTTON(
			gtk_builder_get_object(builder, "useflat_button"));
	if (gtk_toggle_button_get_active(tbutton) == TRUE) {
		const char *filename;
		entry = GTK_ENTRY(gtk_builder_get_object(builder, "flatname_entry"));
		filename = gtk_entry_get_text(entry);
		if (filename[0] == '\0') {
			gtk_toggle_button_set_active(tbutton, FALSE);
		} else {
			fits *flat_fit;
			progress_bar_set_text(_("Opening flat image..."));
			flat_fit = calloc(1, sizeof(fits));
			if (readfits(filename, flat_fit, NULL)) {
				siril_log_message(_("NOT USING FLAT: cannot open the file\n"));
				free(flat_fit);
				gtk_entry_set_text(entry, "");
			} else {
				if (flat_fit->naxes[2] != gfit.naxes[2]) {
					const char *msg = _("NOT USING FLAT: number of channels is different");
					siril_log_message("%s\n", msg);
					progress_bar_set_text(msg);
					free(flat_fit);
					gtk_entry_set_text(entry, "");
				}
				else {
					com.preprostatus |= USE_FLAT;
					if (single_image_is_loaded())
						com.uniq->flat = flat_fit;
					else
						com.seq.flat = flat_fit;
				}
			}
		}
	}

	// if offset selected
	tbutton = GTK_TOGGLE_BUTTON(
			gtk_builder_get_object(builder, "useoffset_button"));
	if (gtk_toggle_button_get_active(tbutton) == TRUE) {
		const char *filename;
		entry = GTK_ENTRY(gtk_builder_get_object(builder, "offsetname_entry"));
		filename = gtk_entry_get_text(entry);
		if (filename[0] == '\0') {
			gtk_toggle_button_set_active(tbutton, FALSE);
		} else {
			fits *bias_fit;
			progress_bar_set_text(_("Opening offset image..."));
			bias_fit = calloc(1, sizeof(fits));
			if (readfits(filename, bias_fit, NULL)) {
				siril_log_message(_("NOT USING OFFSET: cannot open the file\n"));
				free(bias_fit);
				gtk_entry_set_text(entry, "");
			} else {
				if (bias_fit->naxes[2] != gfit.naxes[2]) {
					const char *msg = _("NOT USING OFFSET: number of channels is different");
					siril_log_message("%s\n", msg);
					progress_bar_set_text(msg);
					free(bias_fit);
					gtk_entry_set_text(entry, "");
				}
				else {
					com.preprostatus |= USE_OFFSET;
					if (single_image_is_loaded())
						com.uniq->offset = bias_fit;
					else
						com.seq.offset = bias_fit;
				}
			}
		}
	}

	if (com.preprostatus == 0)
		return;

	// set output filename (preprocessed file name prefix)
	entry = GTK_ENTRY(gtk_builder_get_object(builder, "preproseqname_entry"));
	struct preprocessing_data *args = malloc(sizeof(struct preprocessing_data));
	siril_log_color_message(_("Preprocessing...\n"), "red");
	gettimeofday(&args->t_start, NULL);

	/* Get parameters */
	GtkWidget *autobutton = lookup_widget("checkbutton_auto_evaluate");
	args->autolevel = gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(autobutton));
	if (args->autolevel) {
		args->normalisation = 1.0f;	// will be updated anyway
	} else {
		GtkEntry *norm_entry = GTK_ENTRY(
				gtk_builder_get_object(builder, "entry_flat_norm"));
		args->normalisation = atof(gtk_entry_get_text(norm_entry));
	}
	GtkToggleButton *CFA;
	GtkSpinButton *sigHot, *sigCold;

	CFA = GTK_TOGGLE_BUTTON(gtk_builder_get_object(builder,"cosmCFACheck"));
	sigHot = GTK_SPIN_BUTTON(gtk_builder_get_object(builder,"spinSigCosmeHot"));
	sigCold = GTK_SPIN_BUTTON(gtk_builder_get_object(builder,"spinSigCosmeCold"));

	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkSigCold"))))
		args->sigma[0] = gtk_spin_button_get_value(sigCold);
	else
		args->sigma[0] = -1.0;

	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkSigHot"))))
		args->sigma[1] = gtk_spin_button_get_value(sigHot);
	else
		args->sigma[1] = -1.0;

	args->is_cfa = gtk_toggle_button_get_active(CFA);

	/****/

	if (single_image_is_loaded()) {
		int success = 0;

		com.uniq->ppprefix = strdup(gtk_entry_get_text(entry));
		// start preprocessing
		set_cursor_waiting(TRUE);
		control_window_switch_to_tab(OUTPUT_LOGS);
		success = seqpreprocess(args) == GPOINTER_TO_INT(0);
		if (success)
			progress_bar_reset_ready();
		else
			progress_bar_set_percent(0.0);
		// end_sequence_prepro is also executed there by seqpreprocess

		free(com.uniq->ppprefix);
		com.uniq->ppprefix = NULL;
		unique_free_preprocessing_data(com.uniq);
	} else {	// sequence, executed in a background thread
		if (com.seq.ppprefix)
			free(com.seq.ppprefix);
		com.seq.ppprefix = strdup(gtk_entry_get_text(entry));

		// start preprocessing
		set_cursor_waiting(TRUE);
		control_window_switch_to_tab(OUTPUT_LOGS);
		start_in_new_thread(seqpreprocess, args);
	}
}

void on_showexcluded_button_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	com.show_excluded = gtk_toggle_button_get_active(togglebutton);
}

void on_ref_frame_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	free_reference_image();
	if ((gtk_toggle_button_get_active(togglebutton) == FALSE)) {
		if (com.seq.reference_image == com.seq.current)
			com.seq.reference_image = -1;
	} else {
		com.seq.reference_image = com.seq.current;
		test_and_allocate_reference_image(-1);
	}
	sequence_list_change_reference();
	adjust_sellabel();	// reference image is named in the label
	writeseqfile(&com.seq);
	drawPlot();		// update plots
}


void on_regTranslationOnly_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	GtkWidget *Algo, *Prefix;

	Algo = lookup_widget("ComboBoxRegInter");
	Prefix = lookup_widget("regseqname_entry");

	gtk_widget_set_sensitive(Algo, !gtk_toggle_button_get_active(togglebutton));
	gtk_widget_set_sensitive(Prefix, !gtk_toggle_button_get_active(togglebutton));
}

void on_seqproc_entry_changed(GtkComboBox *widget, gpointer user_data) {
	gchar *name = gtk_combo_box_text_get_active_text(
			GTK_COMBO_BOX_TEXT(widget));
	gchar msg[256];
	if (name && name[0] != '\0') {
		gchar *type;
		set_cursor_waiting(TRUE);
		const char *ext = get_filename_ext(name);
		if (!strcmp(ext, "ser")) {
			name[strlen(name) - 1] = 'q';
			type = " SER";
#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
		} else if (!check_for_film_extensions(ext)) {
			int len = strlen(ext);
			strncpy(name + strlen(name) - len - 1, "seq", len + 1);
			type = " AVI";
#endif
		} else
			type = "";
		g_snprintf(msg, sizeof(msg), _("Selected %s sequence %s..."), type, name);
		progress_bar_set_text(msg);
		set_seq(name);
		set_cursor_waiting(FALSE);
		progress_bar_reset_ready();
	}
	g_free(name);
}

/* signal handler for the gray window layer change */
void on_notebook1_switch_page(GtkNotebook *notebook, GtkWidget *page,
		guint page_num, gpointer user_data) {
	com.cvport = page_num;
	set_cutoff_sliders_values();// load the previous known values for sliders
	set_display_mode();		// change the mode in the combo box if needed
	redraw(com.cvport, REMAP_ONLY);
	calculate_fwhm(com.vport[com.cvport]);
	fill_sequence_list(&com.seq, com.cvport);
}

void on_checkseqbutton_clicked(GtkButton *button, gpointer user_data) {
	static GtkToggleButton *forceButton = NULL;
	int force;

	if (forceButton == NULL) {
		forceButton = GTK_TOGGLE_BUTTON(lookup_widget("checkforceseq"));
	}
	force = gtk_toggle_button_get_active(forceButton);

	set_cursor_waiting(TRUE);
	progress_bar_set_text(
			_("Searching for sequences in the current working directory..."));
	if (!check_seq(force))
		update_sequences_list(NULL);

	/* it's better to uncheck the force button each time it is used */
	if (force)
		gtk_toggle_button_set_active(forceButton, FALSE);
	progress_bar_reset_ready();
	set_cursor_waiting(FALSE);
}

void on_confirmok_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("confirm_dialog"));
	switch (confirm) {
	case CD_INCALL:
		sequence_setselect_all(TRUE);
		break;

	case CD_EXCALL:
		sequence_setselect_all(FALSE);
		break;

	case CD_NULL:
		break;

	case CD_QUIT:
		undo_flush();
		exit(EXIT_SUCCESS);
		break;
	}

	confirm = CD_NULL;
}

void on_confirmDontShowButton_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {

	com.dontShowConfirm = gtk_toggle_button_get_active(togglebutton);
	set_GUI_misc();
	writeinitfile();
}

void on_confirmcancel_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("confirm_dialog"));
	confirm = CD_NULL;
}

gboolean on_drawingarea_key_press_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {

	/*** ZOOM shortcuts ***/
	double oldzoom;
	//fprintf(stdout, "drawing area key event\n");
	oldzoom = com.zoom_value;
	is_shift_on = FALSE;

	switch (event->keyval) {
	case GDK_KEY_plus:
	case GDK_KEY_KP_Add:
		if (oldzoom < 0)
			com.zoom_value = 1.0;
		else
			com.zoom_value = min(ZOOM_MAX, oldzoom * 2.0);
		break;
	case GDK_KEY_minus:
	case GDK_KEY_KP_Subtract:
		if (oldzoom < 0)
			com.zoom_value = 1.0;
		else
			com.zoom_value = max(ZOOM_MIN, oldzoom / 2.0);
		break;
	case GDK_KEY_equal:
	case GDK_KEY_KP_Multiply:
		com.zoom_value = 1.0;
		break;
	case GDK_KEY_KP_0:
	case GDK_KEY_0:
		com.zoom_value = -1.0;
		break;
	case GDK_KEY_Shift_L:
	case GDK_KEY_Shift_R:
		is_shift_on = TRUE;
		break;
	default:
		//~ fprintf(stdout, "No bind found for key '%x'.\n", event->keyval);
		break;
	}
	if (com.zoom_value != oldzoom) {
		fprintf(stdout, _("new zoom value: %f\n"), com.zoom_value);
		zoomcombo_update_display_for_zoom();
		adjust_vport_size_to_image();
		redraw(com.cvport, REMAP_NONE);
	}
	return FALSE;
}

void on_dialog1_OK(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("dialog1"));
}

void on_button_data_ok_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("data_dialog"));
}

void on_menuitemgray_toggled(GtkCheckMenuItem *checkmenuitem,
		gpointer user_data) {
	if (gtk_check_menu_item_get_active(checkmenuitem))
		gtk_widget_show_all(lookup_widget("main_window"));
	else
		gtk_widget_hide(lookup_widget("main_window"));
}

void on_menuitemcolor_toggled(GtkCheckMenuItem *checkmenuitem,
		gpointer user_data) {
	if (gtk_check_menu_item_get_active(checkmenuitem))
		gtk_widget_show_all(lookup_widget("rgb_window"));
	else
		gtk_widget_hide(lookup_widget("rgb_window"));
}

gboolean rgb_area_popup_menu_handler(GtkWidget *widget) {
	do_popup_rgbmenu(widget, NULL);
	return TRUE;
}

void on_rgb_window_hide(GtkWidget *object, gpointer user_data) {
	GtkCheckMenuItem *rgbcheck = GTK_CHECK_MENU_ITEM(
			gtk_builder_get_object(builder, "menuitemcolor"));
	gtk_check_menu_item_set_active(rgbcheck, FALSE);
}

void on_gray_window_hide(GtkWidget *object, gpointer user_data) {
	GtkCheckMenuItem *graycheck = GTK_CHECK_MENU_ITEM(
			gtk_builder_get_object(builder, "menuitemgray"));
	gtk_check_menu_item_set_active(graycheck, FALSE);
}

void toggle_histogram_window_visibility(GtkToolButton *button, gpointer user_data) {
	GtkWidget *window = lookup_widget("histogram_window");
	set_cursor_waiting(TRUE);
	compute_histo_for_gfit(1);// it needs to be forced in the case where operation like background extraction have been done
	if (gtk_widget_get_visible(window))
		gtk_widget_hide(window);
	else
		gtk_widget_show(window);
	set_cursor_waiting(FALSE);
}

void on_combozoom_changed(GtkComboBox *widget, gpointer user_data) {
	gint active = gtk_combo_box_get_active(widget);
	switch (active) {
	case 0: /* 16:1 */
		com.zoom_value = 16.;
		break;
	case 1: /* 8:1 */
		com.zoom_value = 8.;
		break;
	case 2: /* 4:1 */
		com.zoom_value = 4.;
		break;
	case 3: /* 2:1 */
		com.zoom_value = 2.;
		break;
	case -1:
	case 4: /* 1:1 */
		com.zoom_value = 1.;
		break;
	case 5: /* 1:2 */
		com.zoom_value = .5;
		break;
	case 6: /* 1:4 */
		com.zoom_value = .25;
		break;
	case 7: /* 1:8 */
		com.zoom_value = .125;
		break;
	case 8: /* fit to window */
		com.zoom_value = -1.;
		break;
	}
	fprintf(stdout, "zoom is now %f\n", com.zoom_value);
	adjust_vport_size_to_image();
	redraw(com.cvport, REMAP_NONE);
}

void on_comboboxreglayer_changed(GtkComboBox *widget, gpointer user_data) {
	free_reference_image();
	update_stack_interface();
}

void scrollbars_hadjustment_changed_handler(GtkAdjustment *adjustment,
		gpointer user_data) {
	int i;
	double value = gtk_adjustment_get_value(adjustment);

	for (i = 0; i < MAXVPORT; i++) {
		if (com.hadj[i] != adjustment) {
			gtk_adjustment_set_value(com.hadj[i], value);
		}
	}
}

void scrollbars_vadjustment_changed_handler(GtkAdjustment *adjustment,
		gpointer user_data) {
	int i;
	double value = gtk_adjustment_get_value(adjustment);
	for (i = 0; i < MAXVPORT; i++) {
		if (com.vadj[i] != adjustment) {
			gtk_adjustment_set_value(com.vadj[i], value);
		}
	}
}

void on_menu_rgb_savefits_activate(GtkMenuItem *menuitem, gpointer user_data) {
	static GtkNotebook* notebookFormat = NULL;
	static GtkWidget *savepopup = NULL;
	GtkToggleButton *b8bit = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_save_fit8"));
	GtkToggleButton *b16bitu = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_save_fit16"));
	GtkToggleButton *b16bits = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_save_fit16s"));

	if (notebookFormat == NULL) {
		notebookFormat = GTK_NOTEBOOK(
				gtk_builder_get_object(builder, "notebookFormat"));
		savepopup = lookup_widget("savepopup");
	}
	if (single_image_is_loaded() || sequence_is_loaded()) {
		gtk_window_set_title(GTK_WINDOW(savepopup), "Saving FITS");
		gtk_widget_show_all(savepopup);
		gtk_notebook_set_current_page(notebookFormat, 2);
		whichminisave = TYPEFITS;
	}
	switch(gfit.bitpix) {
	case BYTE_IMG:
		gtk_toggle_button_set_active(b8bit, TRUE);
		break;
	case SHORT_IMG:
		gtk_toggle_button_set_active(b16bits, TRUE);
		break;
	default:
		gtk_toggle_button_set_active(b16bitu, TRUE);
	}
}

void on_menu_rgb_savetiff_activate(GtkMenuItem *menuitem, gpointer user_data) {
	static GtkNotebook* notebookFormat = NULL;
	static GtkWidget *savepopup = NULL;


	if (notebookFormat == NULL) {
		notebookFormat = GTK_NOTEBOOK(
				gtk_builder_get_object(builder, "notebookFormat"));
		savepopup = lookup_widget("savepopup");
	}

	if (single_image_is_loaded() || sequence_is_loaded()) {
		Set_Programm_name_in_TIFF(); //Write "Siril Version X.Y in Copyright_Txt
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving TIFF"));
		gtk_widget_show_all(savepopup);
		gtk_notebook_set_current_page(notebookFormat, 0);
		whichminisave = TYPETIFF;
	}
}

void on_menu_rgb_save8ppm_activate(GtkMenuItem *menuitem, gpointer user_data) {
	static GtkNotebook* notebookFormat = NULL;
	static GtkWidget *savepopup = NULL;

	if (notebookFormat == NULL) {
		notebookFormat = GTK_NOTEBOOK(
				gtk_builder_get_object(builder, "notebookFormat"));
		savepopup = lookup_widget("savepopup");
	}

	if (single_image_is_loaded() || sequence_is_loaded()) {
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving Netpbm"));
		gtk_widget_show_all(savepopup);
		gtk_notebook_set_current_page(notebookFormat, 3);
		whichminisave = TYPEPNM;
	}
}

void on_menu_rgb_savebmp_activate(GtkMenuItem *menuitem, gpointer user_data) {
	static GtkNotebook* notebookFormat = NULL;
	static GtkWidget *savepopup = NULL;

	if (notebookFormat == NULL) {
		notebookFormat = GTK_NOTEBOOK(
				gtk_builder_get_object(builder, "notebookFormat"));
		savepopup = lookup_widget("savepopup");
	}

	if (single_image_is_loaded() || sequence_is_loaded()) {
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving BMP"));
		gtk_widget_show_all(savepopup);
		gtk_notebook_set_current_page(notebookFormat, 3);
		whichminisave = TYPEBMP;
	}
}

void on_menu_rgb_savejpg_activate(GtkMenuItem *menuitem, gpointer user_data) {
	static GtkNotebook* notebookFormat = NULL;
	static GtkWidget *savepopup = NULL;

	if (notebookFormat == NULL) {
		notebookFormat = GTK_NOTEBOOK(
				gtk_builder_get_object(builder, "notebookFormat"));
		savepopup = lookup_widget("savepopup");
	}

	if (single_image_is_loaded() || sequence_is_loaded()) {
		gtk_window_set_title(GTK_WINDOW(savepopup), _("Saving JPG"));
		if (sequence_is_loaded() && !single_image_is_loaded()) {
			char filename[256];
			/* set the output file name default as the current image.jpg */
			GtkEntry *entry = GTK_ENTRY(lookup_widget("savetxt"));
			seq_get_image_filename(&com.seq, com.seq.current, filename);
			gtk_entry_set_text(entry, filename);
		}
		gtk_widget_show_all(savepopup);
		gtk_notebook_set_current_page(notebookFormat, 1);
		whichminisave = TYPEJPG;
	}
}

void on_savetxt_changed(GtkEditable *editable, gpointer user_data) {
	GtkEntry *entry = GTK_ENTRY(editable);
	GtkWidget *button = lookup_widget("button_savepopup");

	const gchar *name = gtk_entry_get_text(entry);
	gtk_widget_set_sensitive(button, (name[0] != '\0'));
}

void on_button_savepopup_clicked(GtkButton *button, gpointer user_data) {
	set_cursor_waiting(TRUE);
	minisavedial();
	set_cursor_waiting(FALSE);
}

void on_button_cancelpopup_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("savepopup"));
}

void on_removegreen_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded() && isrgb(&gfit))
		gtk_widget_show_all(lookup_widget("SCNR_dialog"));
}

void on_menuitem_satu_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded() && isrgb(&gfit))
		gtk_widget_show_all(lookup_widget("satu_dialog"));
}

void on_satu_cancel_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("satu_dialog"));
}

void on_satu_Apply_clicked(GtkButton *button, gpointer user_data) {
	static GtkComboBox *combo_satu = NULL;
	int combo;

	if (get_thread_run()) {
		siril_log_message(
				"Another task is already in progress, ignoring new request.\n");
		return;
	}

	struct enhance_saturation_data *args = malloc(sizeof(struct enhance_saturation_data));


	args->coeff = gtk_range_get_value(
			GTK_RANGE(gtk_builder_get_object(builder, "scale_satu")));
	gboolean preserve = gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("preserve_bg")));

	if (args->coeff == 0.0) {
		free(args);
		return;
	}
	undo_save_state("Processing: Saturation enhancement (%lf)", args->coeff);

	if (!combo_satu)
		combo_satu = GTK_COMBO_BOX(lookup_widget("combo_saturation"));
	combo = gtk_combo_box_get_active(combo_satu);

	set_cursor_waiting(TRUE);

	switch (combo) {
	case 0:		// Pink-Red to Red-Orange
		args->h_min = 346.0;
		args->h_max = 20.0;
		break;
	case 1:		// Orange-Brown to Yellow
		args->h_min = 21.0;
		args->h_max = 60.0;
		break;
	case 2:		// Yellow-Green to Green-Cyan
		args->h_min = 61.0;
		args->h_max = 200.0;
		break;
	case 3:		// Cyan
		args->h_min = 170.0;
		args->h_max = 200.0;
		break;
	case 4:		// Cyan-Blue to Blue-Magenta
		args->h_min = 201.0;
		args->h_max = 280.0;
		break;
	case 5:		// Magenta to Pink
		args->h_min = 281.0;
		args->h_max = 345.0;
		break;
	default:
	case 6:		// Global
		args->h_min = 0.0;
		args->h_max = 360.0;
	}

	args->fit = &gfit;
	args->preserve = preserve;
	set_cursor_waiting(TRUE);
	start_in_new_thread(enhance_saturation, args);
}

void on_SCNR_dialog_show(GtkWidget *widget, gpointer user_data) {
	GtkComboBox *comboscnr = GTK_COMBO_BOX(
			gtk_builder_get_object(builder, "combo_scnr"));
	int type = gtk_combo_box_get_active(comboscnr);

	if (type == -1)
		gtk_combo_box_set_active(comboscnr, 0);
}

void on_SCNR_Apply_clicked(GtkButton *button, gpointer user_data) {
	/* Type 0: Average Neutral protection
	 * Type 1: Maximum Neutral protection
	 */
	int type = gtk_combo_box_get_active(
			GTK_COMBO_BOX(gtk_builder_get_object(builder, "combo_scnr")));
	GtkToggleButton *light_button = GTK_TOGGLE_BUTTON(
			gtk_builder_get_object(builder, "preserve_light"));
	gboolean preserve = gtk_toggle_button_get_active(light_button);
	double amount = gtk_range_get_value(
			GTK_RANGE(gtk_builder_get_object(builder, "scale_scnr")));

	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return;
	}

	struct scnr_data *args = malloc(sizeof(struct scnr_data));
	undo_save_state("Processing: SCNR (type: %d, amount: %0.2lf, preserve: %s)",
			type, amount, preserve ? "TRUE" : "FALSE");

	args->fit = &gfit;
	args->type = type;
	args->amount = amount;
	args->preserve = preserve;
	set_cursor_waiting(TRUE);
	start_in_new_thread(scnr, args);
}

void on_SCNR_cancel_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("SCNR_dialog"));
}

void on_combo_scnr_changed(GtkComboBoxText *box, gpointer user_data) {
	int type = gtk_combo_box_get_active(
			GTK_COMBO_BOX(gtk_builder_get_object(builder, "combo_scnr")));
	GtkScale *scale = GTK_SCALE(lookup_widget("scale_scnr"));
	GtkLabel *label = GTK_LABEL(lookup_widget("label56"));

	gtk_widget_set_sensitive(GTK_WIDGET(scale), type > 1);
	gtk_widget_set_sensitive(GTK_WIDGET(label), type > 1);
}

#ifdef HAVE_OPENCV
/* Resample */
void on_menuitem_resample_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded())
		gtk_widget_show_all(lookup_widget("resample_dialog"));
}

void on_button_resample_ok_clicked(GtkButton *button, gpointer user_data) {
	double sample[2];
	sample[0] = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_X")));
	sample[1] = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_Y")));
	int interpolation = gtk_combo_box_get_active(
			GTK_COMBO_BOX(lookup_widget("combo_interpolation")));

	set_cursor_waiting(TRUE);
	int toX = round_to_int((sample[0] / 100.0) * gfit.rx);
	int toY = round_to_int((sample[1] / 100.0) * gfit.ry);
	undo_save_state("Processing: Resample (%g - %g)",
			sample[0]/100.0, sample[1]/100.0);
	verbose_resize_gaussian(&gfit, toX, toY, interpolation);
	update_used_memory();
	adjust_vport_size_to_image();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void on_button_resample_close_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("resample_dialog"));
}

void on_spinbutton_resample_X_value_changed(GtkSpinButton *spinbutton,
		gpointer user_data) {
	GtkToggleButton *ratio = GTK_TOGGLE_BUTTON(
			lookup_widget("button_sample_ratio"));
	double xvalue = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_X")));

	if (gtk_toggle_button_get_active(ratio))
		gtk_spin_button_set_value(
				GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_Y")),
				xvalue);
}

void on_spinbutton_resample_Y_value_changed(GtkSpinButton *spinbutton,
		gpointer user_data) {
	GtkToggleButton *ratio = GTK_TOGGLE_BUTTON(
			lookup_widget("button_sample_ratio"));
	double yvalue = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_Y")));

	if (gtk_toggle_button_get_active(ratio))
		gtk_spin_button_set_value(
				GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_X")),
				yvalue);
}

void on_button_sample_ratio_toggled(GtkToggleButton *button, gpointer user_data) {
	double xvalue = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_X")));

	if (gtk_toggle_button_get_active(button))
		gtk_spin_button_set_value(
				GTK_SPIN_BUTTON(lookup_widget("spinbutton_resample_Y")),
				xvalue);
}

/* Rotation */
void on_menuitem_rotation90_activate(GtkMenuItem *menuitem, gpointer user_data) {
	static GtkToggleButton *crop_rotation = NULL;
	int cropped;

	if (crop_rotation == NULL) {
		crop_rotation = GTK_TOGGLE_BUTTON(
				lookup_widget("checkbutton_rotation_crop"));
	}
	cropped = gtk_toggle_button_get_active(crop_rotation);

	set_cursor_waiting(TRUE);
	undo_save_state("Processing: Rotation (90.0deg)");
	verbose_rotate_image(&gfit, 90.0, -1, cropped);	// fast rotation, no interpolation, no crop
	adjust_vport_size_to_image();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void on_menuitem_rotation270_activate(GtkMenuItem *menuitem, gpointer user_data) {
	static GtkToggleButton *crop_rotation = NULL;
	int cropped;

	if (crop_rotation == NULL) {
		crop_rotation = GTK_TOGGLE_BUTTON(
				lookup_widget("checkbutton_rotation_crop"));
	}
	cropped = gtk_toggle_button_get_active(crop_rotation);

	set_cursor_waiting(TRUE);
	undo_save_state("Processing: Rotation (-90.0deg)");
	verbose_rotate_image(&gfit, 270.0, -1, cropped);// fast rotation, no interpolation, no crop
	adjust_vport_size_to_image();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void on_menuitem_rotation_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded())
		gtk_widget_show_all(lookup_widget("rotation_dialog"));
}

void on_button_rotation_close_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("rotation_dialog"));
}

void on_button_rotation_ok_clicked(GtkButton *button, gpointer user_data) {
	static GtkToggleButton *crop_rotation = NULL;
	double angle = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_rotation")));
	int interpolation = gtk_combo_box_get_active(
			GTK_COMBO_BOX(lookup_widget("combo_interpolation_rotation")));
	int cropped;

	if (crop_rotation == NULL) {
		crop_rotation = GTK_TOGGLE_BUTTON(
				lookup_widget("checkbutton_rotation_crop"));
	}
	cropped = gtk_toggle_button_get_active(crop_rotation);

	set_cursor_waiting(TRUE);
	undo_save_state("Processing: Rotation (%.1lfdeg, cropped=%s)",
			angle, cropped ? "TRUE" : "FALSE");
	verbose_rotate_image(&gfit, angle, interpolation, cropped);
	update_used_memory();
	adjust_vport_size_to_image();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

#endif

void on_menuitem_mirrorx_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded()) {
		set_cursor_waiting(TRUE);
		undo_save_state("Processing: Mirror X");
		mirrorx(&gfit, TRUE);
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
		set_cursor_waiting(FALSE);
	}
}

void on_menuitem_mirrory_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded()) {
		set_cursor_waiting(TRUE);
		undo_save_state("Processing: Mirror Y");
		mirrory(&gfit, TRUE);
		redraw(com.cvport, REMAP_ALL);
		redraw_previews();
		set_cursor_waiting(FALSE);
	}
}

void on_menuitem_noise_activate(GtkMenuItem *menuitem, gpointer user_data) {

	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return;
	}

	struct noise_data *args = malloc(sizeof(struct noise_data));

	set_cursor_waiting(TRUE);
	control_window_switch_to_tab(OUTPUT_LOGS);

	args->fit = &gfit;
	args->verbose = TRUE;
	memset(args->bgnoise, 0.0, sizeof(double[3]));
	start_in_new_thread(noise, args);
}

void on_menuitem_stat_activate(GtkMenuItem *menuitem, gpointer user_data) {
	set_cursor_waiting(TRUE);
	computeStat();
	gtk_widget_show_all(lookup_widget("StatWindow"));
	set_cursor_waiting(FALSE);
}
/**************** GUI for Background extraction *******************/

void on_menuitem_bkg_extraction_activate(GtkMenuItem *menuitem,
		gpointer user_data) {
	if (single_image_is_loaded()) {
		update_bkg_interface();
		gtk_widget_show(lookup_widget("Bkg_extract_window"));
	}
}

void on_bkgButtonManual_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {

	update_bkg_interface();
	redraw(com.cvport, REMAP_NONE);
	redraw_previews();
}

void on_bkgCompute_clicked(GtkButton *button, gpointer user_data) {
	static GtkToggleButton *imgbutton = NULL, *bgkAutoButton = NULL;
	gboolean automatic;

	if (imgbutton == NULL) {
		imgbutton = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_bkg_img"));
		bgkAutoButton = GTK_TOGGLE_BUTTON(lookup_widget("bkgButtonAuto"));
	}
	automatic = gtk_toggle_button_get_active(bgkAutoButton);

	if (!gtk_toggle_button_get_active(imgbutton)) {
		char *msg =	siril_log_message(_("Background cannot be extracted"
				" from itself. Please, click on Show Image\n"));
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		return;
	}

	set_cursor_waiting(TRUE);
	bkgExtractBackground(&wfit[0], automatic);
	redraw(com.cvport, REMAP_NONE);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void on_button_bkg_correct_clicked(GtkButton *button, gpointer user_data) {
	static GtkToggleButton *imgbutton = NULL;

	if (imgbutton == NULL)
		imgbutton = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_bkg_img"));

	if (!gtk_toggle_button_get_active(imgbutton)) {
		char *msg =
				siril_log_message(
						_("Please, apply correction on the image by clicking on Show Image\n"));
		show_dialog(msg, _("Error"), "gtk-dialog-error");
		return;
	}

	int layer;
	GtkComboBox *operation = GTK_COMBO_BOX(lookup_widget("combo_correction"));
	int correction = gtk_combo_box_get_active(operation);

	set_cursor_waiting(TRUE);
	undo_save_state("Processing: Background extraction (Correction: %s)",
			correction ? "Division" : "Subtraction");

	switch (correction) {
	default:
	case 0:
		for (layer = 0; layer < com.uniq->nb_layers; layer++) {
			if (sub_background(&gfit, &wfit[0], layer)) {
				set_cursor_waiting(FALSE);
				return;
			}
		}
		siril_log_message(_("Subtraction done ...\n"));
		break;
	case 1:
		if (ndiv(&gfit, &wfit[0])) {
			set_cursor_waiting(FALSE);
			return;
		}
		siril_log_message(_("Division done ...\n"));
		break;
	}

	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	set_cursor_waiting(FALSE);
}

void on_checkbutton_bkg_boxes_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	com.grad_boxes_drawn = gtk_toggle_button_get_active(togglebutton);
	redraw(com.cvport, REMAP_NONE);
	redraw_previews();
}

void on_radiobutton_bkg_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	memcpy(&wfit[1], &gfit, sizeof(fits));
	memcpy(&gfit, &wfit[0], sizeof(fits));
	memcpy(&wfit[0], &wfit[1], sizeof(fits));
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
}

void on_combobox_gradient_inter_changed(GtkComboBox* box, gpointer user_data) {
	GtkNotebook* notebook = GTK_NOTEBOOK(
			gtk_builder_get_object(builder, "notebook_bkg"));

	gtk_notebook_set_current_page(notebook, gtk_combo_box_get_active(box));
}

void on_bkgClearSamples_clicked(GtkButton *button, gpointer user_data) {
	static GtkToggleButton *imgbutton = NULL, *bkgbutton = NULL;
	int remap_option;

	if (imgbutton == NULL) {
		imgbutton = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_bkg_img"));
		bkgbutton = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_bkg_bkg"));
	}
	gtk_widget_set_sensitive(lookup_widget("frame_bkg_tools"), FALSE);
	gtk_widget_set_sensitive(lookup_widget("button_bkg_correct"), FALSE);
	remap_option = REMAP_NONE;

	set_cursor_waiting(TRUE);
	if (gtk_toggle_button_get_active(bkgbutton)) {
		gtk_toggle_button_set_active(imgbutton, TRUE);
		remap_option = REMAP_ALL;
	}
	clearSamples();
	redraw(com.cvport, remap_option);
	redraw_previews();
	clearfits(&wfit[0]);
	set_cursor_waiting(FALSE);
}

void on_button_bkg_extract_close_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("Bkg_extract_window"));
}

void on_Bkg_extract_window_hide(GtkWidget *widget, gpointer user_data) {
	static GtkToggleButton *imgbutton = NULL, *bkgbutton = NULL, *bgkManButton = NULL;
	int remap_option;

	if (imgbutton == NULL) {
		imgbutton = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_bkg_img"));
		bkgbutton = GTK_TOGGLE_BUTTON(lookup_widget("radiobutton_bkg_bkg"));
		bgkManButton = GTK_TOGGLE_BUTTON(lookup_widget("bkgButtonManual"));
	}
	gtk_widget_set_sensitive(lookup_widget("frame_bkg_tools"), FALSE);
	gtk_widget_set_sensitive(lookup_widget("button_bkg_correct"), FALSE);
	gtk_toggle_button_set_active(bgkManButton, TRUE);
	remap_option = REMAP_NONE;

	set_cursor_waiting(TRUE);
	if (gtk_toggle_button_get_active(bkgbutton)) {
		gtk_toggle_button_set_active(imgbutton, TRUE);
		remap_option = REMAP_ALL;
	}
	clearSamples();
	mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
	redraw(com.cvport, remap_option);
	redraw_previews();
	clearfits(&wfit[0]);
	set_cursor_waiting(FALSE);
}

/**********************************************************************/

void on_menu_channel_separation_activate(GtkMenuItem *menuitem,
		gpointer user_data) {
	if (single_image_is_loaded() && isrgb(&gfit))
		gtk_widget_show_all(lookup_widget("extract_channel_dialog"));
}

void on_menuitem_histo_activate(GtkMenuItem *menuitem, gpointer user_data) {
	GtkWidget *window = lookup_widget("histogram_window");
	if (single_image_is_loaded() || sequence_is_loaded()) {
		set_cursor_waiting(TRUE);
		compute_histo_for_gfit(1);
		gtk_widget_show(window);
		set_cursor_waiting(FALSE);
	}
}

void on_menuitemcalibration_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded() && isrgb(&gfit)) {
		initialize_calibration_interface();
		gtk_widget_show(lookup_widget("color_calibration"));
	}
}

void on_menuitemPSF_toggled(GtkCheckMenuItem *checkmenuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(checkmenuitem))
		gtk_widget_show_all(lookup_widget("stars_list_window"));
	else
		gtk_widget_hide(lookup_widget("stars_list_window"));
}

void on_stars_list_window_hide(GtkWidget *object, gpointer user_data) {
	GtkCheckMenuItem *PSFcheck = GTK_CHECK_MENU_ITEM(
			gtk_builder_get_object(builder, "menuitemPSF"));
	gtk_check_menu_item_set_active(PSFcheck, FALSE);
	com.selected_star = -1;
}

void on_sum_button_clicked(GtkButton *button, gpointer user_data) {
	display_PSF(com.stars);
}

gboolean on_Stars_stored_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	GtkTreeView *tree = GTK_TREE_VIEW(
			gtk_builder_get_object(builder, "Stars_stored"));
	GtkTreeSelection *selection = GTK_TREE_SELECTION(
			gtk_builder_get_object(builder, "treeview-selection"));
	GtkTreeModel *tree_stars = gtk_tree_view_get_model(tree);
	GtkTreeIter iter;

	if (event->button == 1) {
		if (com.stars) {
			if (gtk_tree_model_get_iter_first(tree_stars, &iter) == FALSE)
				return FALSE;	//The tree is empty
			if (gtk_tree_selection_get_selected(selection, &tree_stars,
					&iter)) {	//get selected item
				GtkTreePath *path = gtk_tree_model_get_path(tree_stars, &iter);
				int *index = gtk_tree_path_get_indices(path);
				if (!index)
					return FALSE;
				com.selected_star = index[0];
				display_status(com.selected_star);
				gtk_tree_path_free(path);
				redraw(com.cvport, REMAP_NONE);
				redraw_previews();
			}
		}
	}
	return TRUE;
}

void on_Stars_stored_key_release_event(GtkWidget *widget, GdkEventKey *event,
		gpointer user_data) {
	if (event->keyval == GDK_KEY_Delete || event->keyval == GDK_KEY_KP_Delete
			|| event->keyval == GDK_KEY_BackSpace) {
		remove_selected_line();
	}
	move_selected_line();
}

void on_remove_button_clicked(GtkButton *button, gpointer user_data) {
	remove_selected_line();
}

void on_remove_all_button_clicked(GtkButton *button, gpointer user_data) {
	remove_all_lines();
}

void on_process_starfinder_button_clicked(GtkButton *button, gpointer user_data) {
	int layer = RLAYER;
	starFinder sf;

	memset(&sf, 0, sizeof(starFinder));

	if (!single_image_is_loaded() && !sequence_is_loaded()) {
		siril_log_color_message(_("Load an image first, aborted.\n"), "red");
		return;
	}
	set_cursor_waiting(TRUE);
	if (gfit.naxes[2] == 3)
		layer = GLAYER;
	delete_selected_area();
	com.stars = peaker(&gfit, layer, &sf, NULL);
	refresh_stars_list(com.stars);
	set_cursor_waiting(FALSE);
}

void on_export_button_clicked(GtkButton *button, gpointer user_data) {
	int i = 0;
	if (!com.stars)
		return;
	FILE *f = fopen("stars.lst", "w");

	if (f == NULL)
		return;
	while (com.stars[i]) {
		fprintf(f,
				"%d\t%d\t%10.6f %10.6f %10.2f %10.2f %10.2f %10.2f %3.2f %10.3e\n",
				i + 1, com.stars[i]->layer, com.stars[i]->B, com.stars[i]->A,
				com.stars[i]->xpos, com.stars[i]->ypos, com.stars[i]->fwhmx,
				com.stars[i]->fwhmy, com.stars[i]->angle, com.stars[i]->rmse);
		i++;
	}
	fclose(f);
	siril_log_message(_("The file stars.lst has been created.\n"));
}

void on_stars_list_window_show(GtkWidget *widget, gpointer user_data) {
	fill_stars_list(com.stars);
}

void on_button_stars_list_ok_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("stars_list_window"));
}

void on_extract_channel_button_close_clicked(GtkButton *button,
		gpointer user_data) {
	gtk_widget_hide(lookup_widget("extract_channel_dialog"));
}

void on_extract_channel_button_ok_clicked(GtkButton *button, gpointer user_data) {
	static GtkEntry *channel_extract_entry[3] = { NULL, NULL, NULL };
	static GtkComboBox *combo_extract_channel = NULL;

	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return;
	}

	struct extract_channels_data *args = malloc(
			sizeof(struct extract_channels_data));

	if (combo_extract_channel == NULL) {
		combo_extract_channel = GTK_COMBO_BOX(
				lookup_widget("combo_extract_colors"));
		channel_extract_entry[0] = GTK_ENTRY(
				lookup_widget("Ch1_extract_channel_entry"));
		channel_extract_entry[1] = GTK_ENTRY(
				lookup_widget("Ch2_extract_channel_entry"));
		channel_extract_entry[2] = GTK_ENTRY(
				lookup_widget("Ch3_extract_channel_entry"));
	}

	args->type = gtk_combo_box_get_active(combo_extract_channel);
	args->str_type = gtk_combo_box_get_active_id(combo_extract_channel);

	args->channel[0] = gtk_entry_get_text(channel_extract_entry[0]);
	args->channel[1] = gtk_entry_get_text(channel_extract_entry[1]);
	args->channel[2] = gtk_entry_get_text(channel_extract_entry[2]);

	if ((args->channel[0][0] != '\0') && (args->channel[1][0] != '\0')
			&& (args->channel[2][0] != '\0')) {
		args->fit = calloc(1, sizeof(fits));
		set_cursor_waiting(TRUE);
		copyfits(&gfit, args->fit, CP_ALLOC | CP_COPYA | CP_FORMAT, 0);
		start_in_new_thread(extract_channels, args);
	}
}

/******************* POPUP GRAY MENU *******************************/

void on_menu_gray_psf_activate(GtkMenuItem *menuitem, gpointer user_data) {
	char msg[512];
	fitted_PSF *result = NULL;
	int layer = match_drawing_area_widget(com.vport[com.cvport], FALSE);
	char *str;

	if (layer == -1)
		return;
	if (!(com.selection.h && com.selection.w))
		return;
	if (com.selection.w > 300 || com.selection.h > 300) {
		show_dialog(
				_("Current selection is too large.\n"
						"To determine the PSF, please make a selection around a star.\n"),
						_("Warning"), "gtk-dialog-warning");
		return;
	}
	result = psf_get_minimisation(&gfit, layer, &com.selection);
	if (!result)
		return;

	if (com.magOffset > 0.0)
		str = "true reduced";
	else
		str = "relative";
	g_snprintf(msg, sizeof(msg),
			_("Centroid Coordinates:\n\t\tx0=%.2fpx\n\t\ty0=%.2fpx\n\n"
					"Full Width Half Maximum:\n\t\tFWHMx=%.2f%s\n\t\tFWHMy=%.2f%s\n\n"
					"Angle:\n\t\t%0.2fdeg\n\n"
					"Background Value:\n\t\tB=%.6f\n\n"
					"Maximal Intensity:\n\t\tA=%.6f\n\n"
					"Magnitude (%s):\n\t\tm=%.2f\n\n"
					"RMSE:\n\t\tRMSE=%.3e"), result->x0 + com.selection.x,
			com.selection.y + com.selection.h - result->y0, result->fwhmx,
			result->units, result->fwhmy, result->units, result->angle, result->B,
			result->A, str, result->mag + com.magOffset, result->rmse);
	show_data_dialog(msg, "PSF Results");
	free(result);
}

void on_menu_gray_seqpsf_activate(GtkMenuItem *menuitem, gpointer user_data) {
	process_seq_psf(0);
}

void on_menu_gray_pick_star_activate(GtkMenuItem *menuitem, gpointer user_data) {
	int layer = match_drawing_area_widget(com.vport[com.cvport], FALSE);
	int new_index;
	GtkCheckMenuItem *PSFcheck = GTK_CHECK_MENU_ITEM(
			gtk_builder_get_object(builder, "menuitemPSF"));
	GtkWidget *window = lookup_widget("stars_list_window");

	if (layer != -1) {
		if (!(com.selection.h && com.selection.w))
			return;
		if (com.selection.w > 300 || com.selection.h > 300) {
			char *msg =
					siril_log_message(
							_("Current selection is too large.\nTo determine the PSF, please make a selection around a star.\n"));
			show_dialog(msg, _("Warning"), "gtk-dialog-warning");
			return;
		}
		fitted_PSF *new_star = add_star(&gfit, layer, &new_index);
		if (new_star) {
			add_star_to_list(new_star);
			if (!(gtk_widget_get_visible(window)))//We open the stars_list_window
				gtk_widget_show_all(window);
			gtk_check_menu_item_set_active(PSFcheck, TRUE);
		} else
			return;
	}
	redraw(com.cvport, REMAP_NONE);
}

void on_menu_gray_crop_activate(GtkMenuItem *menuitem, gpointer user_data) {
	undo_save_state("Processing: Crop (x=%d, y=%d, w=%d, h=%d)",
			com.selection.x, com.selection.y, com.selection.w, com.selection.h);
	crop(&gfit, &com.selection);
	delete_selected_area();
	adjust_cutoff_from_updated_gfit();
	redraw(com.cvport, REMAP_ALL);
	redraw_previews();
	update_used_memory();
}

void on_menu_gray_crop_seq_activate(GtkMenuItem *menuitem, gpointer user_data) {
	gtk_widget_show(lookup_widget("crop_dialog"));
}

void on_menu_gray_stat_activate(GtkMenuItem *menuitem, gpointer user_data) {
	computeStat();
	gtk_widget_show_all(lookup_widget("StatWindow"));
}

/************************* GUI for FFT ********************************/

void on_button_fft_apply_clicked(GtkButton *button, gpointer user_data) {
	const char *mag, *phase;
	char *type, page;
	int type_order = -1;
	static GtkToggleButton *order = NULL;
	static GtkNotebook* notebookFFT = NULL;

	if (get_thread_run()) {
		siril_log_message(
				_("Another task is already in progress, ignoring new request.\n"));
		return;
	}

	if (notebookFFT == NULL) {
		notebookFFT = GTK_NOTEBOOK(
				gtk_builder_get_object(builder, "notebook_fft"));
		order = GTK_TOGGLE_BUTTON(
				gtk_builder_get_object(builder, "fft_centered"));
	}

	page = gtk_notebook_get_current_page(notebookFFT);

	if (page == 0) {
		if (sequence_is_loaded()) {
			char *msg = siril_log_message(_("FFT does not work with sequences !\n"));
			show_dialog(msg, _("Error"), "gtk-dialog-error");
			set_cursor_waiting(FALSE);
			return;
		}
		if (!single_image_is_loaded()) {
			char *msg = siril_log_message(_("Open an image first !\n"));
			show_dialog(msg, _("Error"), "gtk-dialog-error");
			set_cursor_waiting(FALSE);
			return;

		}

		GtkEntry *entry_mag = GTK_ENTRY(lookup_widget("fftd_mag_entry"));
		GtkEntry *entry_phase = GTK_ENTRY(lookup_widget("fftd_phase_entry"));

		if (gtk_toggle_button_get_active(order) == TRUE)
			type_order = 0;		// centered fftd
		else
			type_order = 1;	// regular fftd
		type = strdup("fftd");
		mag = gtk_entry_get_text(entry_mag);
		phase = gtk_entry_get_text(entry_phase);
	} else {
		type = strdup("ffti");
		mag = gtk_file_chooser_get_filename(
				GTK_FILE_CHOOSER(lookup_widget("filechooser_mag")));
		phase = gtk_file_chooser_get_filename(
				GTK_FILE_CHOOSER(lookup_widget("filechooser_phase")));

		if (mag == NULL || phase == NULL) {
			char *msg = siril_log_message(_("Select magnitude and phase before !\n"));
			show_dialog(msg, _("Error"), "gtk-dialog-error");
			set_cursor_waiting(FALSE);
			free(type);
			return;
		}
		close_single_image();
		open_single_image(mag);
	}

	if ((mag != NULL) && (phase != NULL)) {
		set_cursor_waiting(TRUE);
		struct fft_data *args = malloc(sizeof(struct fft_data));
		args->fit = &gfit;
		args->type = type;
		args->modulus = mag;
		args->phase = phase;
		args->type_order = type_order;
		set_cursor_waiting(TRUE);
		start_in_new_thread(fourier_transform, args);
	}
}

void on_button_fft_close_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("dialog_FFT"));
}

void on_menuitem_fft_activate(GtkMenuItem *menuitem, gpointer user_data) {
	GtkFileChooserButton *magbutton, *phasebutton;

	magbutton = GTK_FILE_CHOOSER_BUTTON(lookup_widget("filechooser_mag"));
	phasebutton = GTK_FILE_CHOOSER_BUTTON(lookup_widget("filechooser_phase"));
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(magbutton), com.wd);
	gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(phasebutton), com.wd);
	gtk_widget_show_all(lookup_widget("dialog_FFT"));
}

void on_menuitem_medianfilter_activate(GtkMenuItem *menuitem,
		gpointer user_data) {
	if (single_image_is_loaded())
		gtk_widget_show(lookup_widget("Median_dialog"));
}

void on_spin_w_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	gtk_widget_set_sensitive(lookup_widget("button_apply_w"), TRUE);
}

void on_menuitem_wavelets_activate(GtkMenuItem *menuitem, gpointer user_data) {

	if (single_image_is_loaded()) {
		reset_scale_w();
		gtk_widget_show_all(lookup_widget("wavelets_dialog"));
	}
}

void on_wavelets_dialog_hide(GtkWidget *widget, gpointer user_data) {
	gtk_widget_set_sensitive(lookup_widget("grid_w"), FALSE);
}

void on_button_apply_w_clicked(GtkButton *button, gpointer user_data) {
	update_wavelets();
	gtk_widget_set_sensitive(lookup_widget("button_apply_w"), FALSE);
}

void on_button_reset_w_clicked(GtkButton *button, gpointer user_data) {
	float scale[6];
	static GtkSpinButton *spin_w[6] = { NULL, NULL, NULL, NULL, NULL, NULL };
	int i;

	if (spin_w[0] == NULL) {
		spin_w[0] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w0"));
		spin_w[1] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w1"));
		spin_w[2] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w2"));
		spin_w[3] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w3"));
		spin_w[4] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w4"));
		spin_w[5] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w5"));
	}

	for (i = 0; i < 6; i++)
		scale[i] = (float) gtk_spin_button_get_value(spin_w[i]);

	if (scale[0] == 1.f && scale[1] == 1.f && scale[2] == 1.f && scale[3] == 1.f
			&& scale[4] == 1.f && scale[5] == 1.f)
		return;
	reset_scale_w();
	update_wavelets();
}

void on_button_ok_w_clicked(GtkButton *button, gpointer user_data) {
	int need_to_be_updated = gtk_widget_get_sensitive(
			lookup_widget("button_apply_w"));
	if (need_to_be_updated) {
		update_wavelets();
		gtk_widget_set_sensitive(lookup_widget("button_apply_w"), FALSE);
	}
	gtk_widget_hide(lookup_widget("wavelets_dialog"));
}

void on_button_cancel_w_clicked(GtkButton *button, gpointer user_data) {
	float scale[6];
	static GtkSpinButton *spin_w[6] = { NULL, NULL, NULL, NULL, NULL, NULL };
	int i;

	if (spin_w[0] == NULL) {
		spin_w[0] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w0"));
		spin_w[1] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w1"));
		spin_w[2] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w2"));
		spin_w[3] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w3"));
		spin_w[4] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w4"));
		spin_w[5] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spin_w5"));
	}

	for (i = 0; i < 6; i++)
		scale[i] = (float) gtk_spin_button_get_value(spin_w[i]);

	if (!(scale[0] == 1.f && scale[1] == 1.f && scale[2] == 1.f
			&& scale[3] == 1.f && scale[4] == 1.f && scale[5] == 1.f)) {
		if (gtk_widget_get_sensitive(lookup_widget("grid_w")) == TRUE) {
			reset_scale_w();
			update_wavelets();
		}
	}
	gtk_widget_hide(lookup_widget("wavelets_dialog"));
}

void on_button_compute_w_clicked(GtkButton *button, gpointer user_data) {
	int Type_Transform, Nbr_Plan, maxplan, mins, i;
	int nb_chan = gfit.naxes[2];
	float *Imag;
	char *File_Name_Transform[3] = { "r_rawdata.wave", "g_rawdata.wave",
			"b_rawdata.wave" }, *dir[3];
	const char *tmpdir;

	assert(nb_chan == 1 || nb_chan == 3);

	tmpdir = g_get_tmp_dir();

	Nbr_Plan = gtk_spin_button_get_value(
			GTK_SPIN_BUTTON(lookup_widget("spinbutton_plans_w")));
	Type_Transform = gtk_combo_box_get_active(
			GTK_COMBO_BOX(lookup_widget("combobox_type_w"))) + 1;

	mins = min(gfit.rx, gfit.ry);
	maxplan = log(mins) / log(2) - 2;

	if (Nbr_Plan > maxplan) {
		char *msg = siril_log_message(
				_("Wavelet: maximum number of plans for this image size is %d\n"),
				maxplan);
		show_dialog(msg, _("Warning"), "gtk-dialog-warning");
		Nbr_Plan = maxplan;
		gtk_spin_button_set_value(
				GTK_SPIN_BUTTON(lookup_widget("spinbutton_plans_w")), Nbr_Plan);
	}

	if (Type_Transform != TO_PAVE_LINEAR && Type_Transform != TO_PAVE_BSPLINE) {
		char *msg = siril_log_message(_("Wavelet: type must be %d or %d\n"),
		TO_PAVE_LINEAR, TO_PAVE_BSPLINE);
		show_dialog(msg, _("Warning"), "gtk-dialog-warning");
	}

	set_cursor_waiting(TRUE);

	Imag = (float *) malloc(gfit.rx * gfit.ry * sizeof(float));

	for (i = 0; i < nb_chan; i++) {
		dir[i] = malloc(strlen(tmpdir) + strlen(File_Name_Transform[i]) + 2);
		strcpy(dir[i], tmpdir);
		strcat(dir[i], "/");
		strcat(dir[i], File_Name_Transform[i]);
		wavelet_transform_file(Imag, gfit.ry, gfit.rx, dir[i], Type_Transform, Nbr_Plan,
				gfit.pdata[i]);
		free(dir[i]);
	}

	free(Imag);
	Imag = NULL;
	gtk_widget_set_sensitive(lookup_widget("grid_w"), TRUE);
	set_cursor_waiting(FALSE);
	return;
}

/****************** GUI for Wavelet Layers Extraction *****************/

void on_menu_wavelet_separation_activate(GtkMenuItem *menuitem,
		gpointer user_data) {

	if (single_image_is_loaded()) {
		reset_scale_w();
		gtk_widget_show_all(lookup_widget("extract_wavelets_layers_dialog"));
	}
}

void on_button_extract_w_ok_clicked(GtkButton *button, gpointer user_data) {
	fits *fit;
	int Nbr_Plan, Type, maxplan, mins, i;
	static GtkSpinButton *Spin_Nbr_Plan = NULL;
	static GtkComboBox *Combo_Wavelets_Type = NULL;

	if (Spin_Nbr_Plan == NULL) {
		Spin_Nbr_Plan = GTK_SPIN_BUTTON(lookup_widget("spinbutton_extract_w"));
		Combo_Wavelets_Type = GTK_COMBO_BOX(
				lookup_widget("combo_interpolation_extract_w"));
	}

	Nbr_Plan = gtk_spin_button_get_value(Spin_Nbr_Plan);
	Type = gtk_combo_box_get_active(Combo_Wavelets_Type) + 1;// 1: linear, 2: bspline

	set_cursor_waiting(TRUE);
	mins = min(gfit.rx, gfit.ry);
	maxplan = log(mins) / log(2) - 2;

	if (Nbr_Plan > maxplan) {
		char *msg = siril_log_message(
				_("Wavelet: maximum number of plans for this image size is %d\n"),
				maxplan);
		show_dialog(msg, _("Warning"), "gtk-dialog-warning");
		set_cursor_waiting(FALSE);
		return;
	}
	fit = calloc(1, sizeof(fits));
	copyfits(&gfit, fit, CP_ALLOC | CP_COPYA | CP_FORMAT, 0);

	for (i = 0; i < Nbr_Plan; i++) {
		char filename[256];

		g_snprintf(filename, sizeof(filename), "layer%02d", i);
		get_wavelet_layers(fit, Nbr_Plan, i, Type, -1);
		savefits(filename, fit);
	}
	clearfits(fit);
	update_used_memory();
	set_cursor_waiting(FALSE);
}

void on_button_extract_w_close_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("extract_wavelets_layers_dialog"));
}

/********************** GUI for Median Filter *************************/

void on_Median_cancel_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("Median_dialog"));
}

void on_Median_Apply_clicked(GtkButton *button, gpointer user_data) {
	int combo_size = gtk_combo_box_get_active(
			GTK_COMBO_BOX(
					gtk_builder_get_object(builder, "combo_ksize_median")));
	double amount = gtk_range_get_value(
			GTK_RANGE(gtk_builder_get_object(builder, "scale_median")));
	int iterations = round_to_int(
			gtk_spin_button_get_value(
					GTK_SPIN_BUTTON(
							gtk_builder_get_object(builder,
									"median_button_iterations"))));

	if (get_thread_run()) {
		siril_log_message(
				_(	"Another task is already in progress, ignoring new request.\n"));
		return;
	}

	struct median_filter_data *args = malloc(sizeof(struct median_filter_data));

	switch (combo_size) {
	default:
	case 0:
		args->ksize = 3;
		break;
	case 1:
		args->ksize = 5;
		break;
	case 2:
		args->ksize = 7;
		break;
	case 3:
		args->ksize = 9;
		break;
	case 4:
		args->ksize = 11;
		break;
	case 5:
		args->ksize = 13;
		break;
	case 6:
		args->ksize = 15;
		break;
	}
	undo_save_state("Processing: Median Filter (filter=%dx%d px)",
			args->ksize, args->ksize);

	args->fit = &gfit;
	args->amount = amount;
	args->iterations = iterations;
	set_cursor_waiting(TRUE);
	start_in_new_thread(median_filter, args);

}

/****************** GUI for Cosmetic Correction *********************/

void on_menuitem_cosmetic_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (sequence_is_loaded()) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkCosmeticSeq")), TRUE);
	}
	else if (single_image_is_loaded()) {
		// not a processing result
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkCosmeticSeq")), FALSE);
	}
	else
		return;
	gtk_widget_show(lookup_widget("cosmetic_dialog"));
}

void on_button_cosmetic_close_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("cosmetic_dialog"));
}

void on_checkSigCosme_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	static GtkWidget *cosmeticApply = NULL;
	static GtkToggleButton *checkCosmeSigCold = NULL;
	static GtkToggleButton *checkCosmeSigHot = NULL;
	gboolean checkCold, checkHot;

	if (cosmeticApply == NULL) {
		cosmeticApply = lookup_widget("button_cosmetic_ok");
		checkCosmeSigCold = GTK_TOGGLE_BUTTON(lookup_widget("checkSigColdBox"));
		checkCosmeSigHot = GTK_TOGGLE_BUTTON(lookup_widget("checkSigHotBox"));
	}
	checkCold = gtk_toggle_button_get_active(checkCosmeSigCold);
	checkHot = gtk_toggle_button_get_active(checkCosmeSigHot);
	gtk_widget_set_sensitive(cosmeticApply, checkCold || checkHot);
}

void on_button_cosmetic_ok_clicked(GtkButton *button, gpointer user_data) {
	GtkEntry *cosmeticSeqEntry;
	GtkToggleButton *CFA, *seq;
	GtkSpinButton *sigma[2];
	GtkAdjustment *adjCosmeAmount;

	CFA = GTK_TOGGLE_BUTTON(gtk_builder_get_object(builder,"cosmCFACheckBox"));
	sigma[0] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder,"spinSigCosmeColdBox"));
	sigma[1] = GTK_SPIN_BUTTON(gtk_builder_get_object(builder,"spinSigCosmeHotBox"));
	seq = GTK_TOGGLE_BUTTON(lookup_widget("checkCosmeticSeq"));
	cosmeticSeqEntry = GTK_ENTRY(lookup_widget("entryCosmeticSeq"));
	adjCosmeAmount = GTK_ADJUSTMENT(gtk_builder_get_object(builder, "adjCosmeAmount"));

	struct cosmetic_data *args = malloc(sizeof(struct cosmetic_data));

	if (gtk_toggle_button_get_active(
			GTK_TOGGLE_BUTTON(lookup_widget("checkSigColdBox"))))
		args->sigma[0] = gtk_spin_button_get_value(sigma[0]);
	else
		args->sigma[0] = -1.0;

	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget("checkSigHotBox"))))
		args->sigma[1] = gtk_spin_button_get_value(sigma[1]);
	else
		args->sigma[1] = -1.0;

	args->is_cfa = gtk_toggle_button_get_active(CFA);
	args->amount = gtk_adjustment_get_value(adjCosmeAmount);

	args->fit = &gfit;
	args->seqEntry = gtk_entry_get_text(cosmeticSeqEntry);
	set_cursor_waiting(TRUE);

	if (gtk_toggle_button_get_active(seq) && sequence_is_loaded()) {
		if (args->seqEntry && args->seqEntry[0] == '\0')
			args->seqEntry = "cc_";
		apply_cosmetic_to_sequence(args);
	} else {
		undo_save_state("Processing: Cosmetic Correction");
		start_in_new_thread(autoDetectThreaded, args);
	}
}

/***************** GUI for Canon Banding Reduction ********************/

void on_menuitem_fixbanding_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (sequence_is_loaded()) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkBandingSeq")), TRUE);
	}
	else if (single_image_is_loaded()) {
		// not a processing result
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(lookup_widget("checkBandingSeq")), FALSE);
	}
	else
		return;
	gtk_widget_show(lookup_widget("canon_fixbanding_dialog"));
}

void on_button_ok_fixbanding_clicked(GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("canon_fixbanding_dialog"));
}

void on_button_apply_fixbanding_clicked(GtkButton *button, gpointer user_data) {
	static GtkRange *range_amount = NULL;
	static GtkRange *range_invsigma = NULL;
	static GtkToggleButton *toggle_protect_highlights_banding = NULL,
		*vertical = NULL, *seq = NULL;
	static GtkEntry *bandingSeqEntry = NULL;
	double amount, invsigma;
	gboolean protect_highlights;

	if (get_thread_run()) {
		siril_log_message(
				_(	"Another task is already in progress, ignoring new request.\n"));
		return;
	}

	struct banding_data *args = malloc(sizeof(struct banding_data));

	if (range_amount == NULL) {
		range_amount = GTK_RANGE(lookup_widget("scale_fixbanding_amount"));
		range_invsigma = GTK_RANGE(lookup_widget("scale_fixbanding_invsigma"));
		toggle_protect_highlights_banding = GTK_TOGGLE_BUTTON(
				lookup_widget("checkbutton_fixbanding"));
		vertical = GTK_TOGGLE_BUTTON(lookup_widget("checkBandingVertical"));
		seq = GTK_TOGGLE_BUTTON(lookup_widget("checkBandingSeq"));
		bandingSeqEntry = GTK_ENTRY(lookup_widget("entryBandingSeq"));
	}
	amount = gtk_range_get_value(range_amount);
	invsigma = gtk_range_get_value(range_invsigma);
	protect_highlights = gtk_toggle_button_get_active(
			toggle_protect_highlights_banding);

	if (!protect_highlights)
		undo_save_state("Processing: Canon Banding Reduction (amount=%.2lf)", amount);
	else
		undo_save_state("Processing: Canon Banding Reduction (amount=%.2lf, Protect=TRUE, invsigma=%.2lf)",
				amount, invsigma);

	args->fit = &gfit;
	args->protect_highlights = protect_highlights;
	args->amount = amount;
	args->sigma = invsigma;
	args->applyRotation = gtk_toggle_button_get_active(vertical);
	args->seqEntry = gtk_entry_get_text(bandingSeqEntry);
	set_cursor_waiting(TRUE);

	if (gtk_toggle_button_get_active(seq) && sequence_is_loaded()) {
		if (args->seqEntry && args->seqEntry[0] == '\0')
			args->seqEntry = "unband_";
		apply_banding_to_sequence(args);
	} else {
		start_in_new_thread(BandingEngineThreaded, args);
	}
}

void on_checkbutton_fixbanding_toggled(GtkToggleButton *togglebutton,
		gpointer user_data) {
	static GtkWidget *bandingHighlightBox = NULL;
	gboolean is_active;

	if (bandingHighlightBox == NULL) {
		bandingHighlightBox = lookup_widget("bandingHighlightBox");
	}

	is_active = gtk_toggle_button_get_active(togglebutton);
	gtk_widget_set_sensitive(bandingHighlightBox, is_active);
}

void on_select_convert_button_clicked(GtkButton *button, gpointer user_data) {
	whichdial = OD_CONVERT;
	opendial();
}

void on_clear_convert_button_clicked(GtkButton *button, gpointer user_data) {
	get_convert_list_store();
	gtk_list_store_clear(liststore_convert);
	check_for_conversion_form_completeness();
}

void on_remove_convert_button_clicked(GtkWidget *button, gpointer user_data) {
	GtkTreeSelection *selection = GTK_TREE_SELECTION(
			gtk_builder_get_object(builder, "treeview-selection5"));
	GtkTreeIter iter;
	GtkTreeModel *model;
	gchararray string;

	if (gtk_tree_selection_get_selected(selection, &model, &iter)) {//get selected item
		gtk_tree_model_get(model, &iter, 0, &string, -1);
		if (string) {
			gtk_list_store_remove(GTK_LIST_STORE(model), &iter);
			gtk_tree_selection_unselect_all(selection);
		}
	}
	check_for_conversion_form_completeness();
}

void on_spinCPU_value_changed (GtkSpinButton *spinbutton, gpointer user_data) {
	com.max_thread = (int)gtk_spin_button_get_value(spinbutton);
}

/*** GUI for crop sequence */
void on_crop_Apply_clicked (GtkButton *button, gpointer user_data) {
	if (get_thread_run()) {
		siril_log_message(
				"Another task is already in progress, ignoring new request.\n");
		return;
	}

#if defined(HAVE_FFMS2_1) || defined(HAVE_FFMS2_2)
	if (com.seq.type == SEQ_AVI) {
		siril_log_message("Crop does not work with avi film. Please, convert your file to SER first.\n");
		return;
	}
#endif
	if (com.seq.type == SEQ_INTERNAL) {
		siril_log_message("Not a valid sequence for cropping.\n");
	}

	struct crop_sequence_data *args = malloc(sizeof(struct crop_sequence_data));

	GtkEntry *cropped_entry = GTK_ENTRY(lookup_widget("cropped_entry"));

	args->seq = &com.seq;
	args->area = &com.selection;
	args->prefix = gtk_entry_get_text(cropped_entry);

	set_cursor_waiting(TRUE);
	start_in_new_thread(crop_sequence, args);
}

void on_crop_close_clicked (GtkButton *button, gpointer user_data) {
	gtk_widget_hide(lookup_widget("crop_dialog"));
}

void on_undo_item_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded() && is_undo_available()) {
		set_cursor_waiting(TRUE);
		undo_display_data(UNDO);
		set_cursor_waiting(FALSE);
	}

	/* update menus */
	update_MenuItem();
}

void on_redo_item_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded() && is_redo_available()) {
		set_cursor_waiting(TRUE);
		undo_display_data(REDO);
		set_cursor_waiting(FALSE);
	}

	/* update menus */
	update_MenuItem();
}

void on_undo_item1_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded() && is_undo_available()) {
		set_cursor_waiting(TRUE);
		undo_display_data(UNDO);
		set_cursor_waiting(FALSE);
	}
}

void on_redo_item1_activate(GtkMenuItem *menuitem, gpointer user_data) {
	if (single_image_is_loaded() && is_redo_available()) {
		set_cursor_waiting(TRUE);
		undo_display_data(REDO);
		set_cursor_waiting(FALSE);
	}
}

void on_darkThemeCheck_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	com.have_dark_theme = gtk_toggle_button_get_active(togglebutton);
}

void on_entryAviWidth_changed(GtkEditable *editable, gpointer user_data) {
	double ratio, width, height;
	char c_height[6];
	GtkEntry *heightEntry = GTK_ENTRY(lookup_widget("entryAviHeight"));

	if (com.selection.w && com.selection.h) return;
	ratio = (double) com.seq.ry / (double) com.seq.rx;
	width = atof(gtk_entry_get_text(GTK_ENTRY(editable)));
	height = ratio * width;
	g_snprintf(c_height, sizeof(c_height), "%d", (int)(height));

	g_signal_handlers_block_by_func(heightEntry, on_entryAviHeight_changed, NULL);
	gtk_entry_set_text(heightEntry, c_height);
	g_signal_handlers_unblock_by_func(heightEntry, on_entryAviHeight_changed, NULL);
}

void on_entryAviHeight_changed(GtkEditable *editable, gpointer user_data) {
	double ratio, width, height;
	char c_width[6];
	GtkEntry *widthEntry = GTK_ENTRY(lookup_widget("entryAviWidth"));

	if (com.selection.w && com.selection.h) return;
	ratio = (double) com.seq.rx / (double) com.seq.ry;
	height = atof(gtk_entry_get_text(GTK_ENTRY(editable)));
	width = ratio * height;
	g_snprintf(c_width, sizeof(c_width), "%d", (int)(width));

	g_signal_handlers_block_by_func(widthEntry, on_entryAviWidth_changed, NULL);
	gtk_entry_set_text(widthEntry, c_width);
	g_signal_handlers_unblock_by_func(widthEntry, on_entryAviWidth_changed, NULL);

}


void on_menu_rgb_align_select(GtkMenuItem *menuitem, gpointer user_data) {
	gboolean sel_is_drawn = ((com.selection.w > 0.0) && (com.selection.h > 0.0));

	gtk_widget_set_sensitive(lookup_widget("rgb_align_dft"), sel_is_drawn);
	gtk_widget_set_sensitive(lookup_widget("rgb_align_psf"), sel_is_drawn);
}

void on_rgb_align_dft_activate(GtkMenuItem *menuitem, gpointer user_data) {
	undo_save_state("Processing: RGB alignment (DFT)");
	rgb_align(1);
}

void on_rgb_align_psf_activate(GtkMenuItem *menuitem, gpointer user_data) {
	undo_save_state("Processing: RGB alignment (PSF)");
	rgb_align(0);
}
