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

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/undo.h"
#include "algos/background_extraction.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "image_interactions.h"
#include "image_display.h"
#include "callbacks.h"
#include "progress_and_log.h"
#include "message_dialog.h"

mouse_status_enum mouse_status;

/* mouse callbacks */
static double margin_size = 10;
static gboolean is_shift_on = FALSE;

static gboolean is_over_the_left_side_of_sel(double zoomedX, double zoomedY, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = margin_size / zoom;
	if ((zoomedX > com.selection.x - s && zoomedX < com.selection.x + s)) {
		if (zoomedY > com.selection.y - s
				&& zoomedY < com.selection.y + com.selection.h + s)
			return TRUE;
	}

	return FALSE;
}

static gboolean is_over_the_right_side_of_sel(double zoomedX, double zoomedY, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = margin_size / zoom;
	if ((zoomedX > com.selection.x + com.selection.w - s
				&& zoomedX < com.selection.x + com.selection.w + s)) {
		if (zoomedY > com.selection.y - s
				&& zoomedY < com.selection.y + com.selection.h + s)
			return TRUE;
	}
	return FALSE;
}

static gboolean is_over_the_bottom_of_sel(double zoomedX, double zoomedY, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = margin_size / zoom;
	if ((zoomedY > com.selection.y + com.selection.h - s
				&& zoomedY < com.selection.y + com.selection.h + s)) {
		if (zoomedX > com.selection.x - s
				&& zoomedX < com.selection.x + com.selection.w + s)
			return TRUE;
	}
	return FALSE;
}

static gboolean is_over_the_top_of_sel(double zoomedX, double zoomedY, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = margin_size / zoom;
	if ((zoomedY > com.selection.y - s && zoomedY < com.selection.y + s)) {
		if (zoomedX > com.selection.x - s
				&& zoomedX < com.selection.x + com.selection.w + s)
			return TRUE;
	}
	return FALSE;
}


static gboolean inimage(GdkEvent *event) {
	double zoom = get_zoom_val();
	return ((GdkEventButton*) event)->x > 0
		&& ((GdkEventButton*) event)->x < gfit.rx * zoom
		&& ((GdkEventButton*) event)->y > 0
		&& ((GdkEventButton*) event)->y < gfit.ry * zoom;
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
	siril_debug_print("selection: %d,%d,\t%dx%d\n", com.selection.x, com.selection.y,
			com.selection.w, com.selection.h);
	for (i = 0; i < _nb_registered_callbacks; ++i) {
		_registered_callbacks[i]();
	}
	redraw(com.cvport, REMAP_NONE);
}

void delete_selected_area() {
	memset(&com.selection, 0, sizeof(rectangle));
	new_selection_zone();
}

/*
 * Button events
 */

static void do_popup_rgbmenu(GtkWidget *my_widget, GdkEventButton *event) {
	static GtkMenu *menu = NULL;

	if (!menu) {
		menu = GTK_MENU(gtk_builder_get_object(builder, "menurgb"));
		gtk_menu_attach_to_widget(GTK_MENU(menu), my_widget, NULL);
	}

#if GTK_CHECK_VERSION(3, 22, 0)
	gtk_menu_popup_at_pointer(GTK_MENU(menu), NULL);
#else
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
#endif
}

/* Gray popup menu */

static void do_popup_graymenu(GtkWidget *my_widget, GdkEventButton *event) {
	static GtkMenu *menu = NULL;
	gboolean selected;

	gboolean is_a_single_image_loaded = single_image_is_loaded() && (!sequence_is_loaded()
			|| (sequence_is_loaded() && (com.seq.current == RESULT_IMAGE
					|| com.seq.current == SCALED_IMAGE)));

	if (!menu) {
		menu = GTK_MENU(gtk_builder_get_object(builder, "menugray"));
		gtk_menu_attach_to_widget(GTK_MENU(menu), my_widget, NULL);
	}

	selected = com.selection.w && com.selection.h;
	gtk_widget_set_sensitive(lookup_widget("undo_item1"), is_undo_available());
	gtk_widget_set_sensitive(lookup_widget("redo_item1"), is_redo_available());
	gtk_widget_set_sensitive(lookup_widget("menu_gray_psf"), selected);
	gtk_widget_set_sensitive(lookup_widget("menu_gray_stat"), is_a_single_image_loaded || sequence_is_loaded());
	gtk_widget_set_sensitive(lookup_widget("menu_gray_seqpsf"), selected);
	gtk_widget_set_sensitive(lookup_widget("menu_gray_pick_star"), selected);
	gtk_widget_set_sensitive(lookup_widget("menu_gray_crop"), selected && is_a_single_image_loaded);
	gtk_widget_set_sensitive(lookup_widget("menu_gray_crop_seq"), selected && sequence_is_loaded());

#if GTK_CHECK_VERSION(3, 22, 0)
	gtk_menu_popup_at_pointer(GTK_MENU(menu), NULL);
#else
	int button, event_time;

	if (event) {
		button = event->button;
		event_time = event->time;
	} else {
		button = 0;
		event_time = gtk_get_current_event_time();
	}

	gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL, button, event_time);
#endif
}


gboolean rgb_area_popup_menu_handler(GtkWidget *widget) {
	do_popup_rgbmenu(widget, NULL);
	return TRUE;
}

gboolean on_drawingarea_button_press_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	if (inimage((GdkEvent *) event)) {
		/* click on RGB image */
		if (widget == com.vport[RGB_VPORT]) {
			if (event->button == 1) {	// left click
				siril_message_dialog(GTK_MESSAGE_INFO, _("Only for visualization"),
						_("The RGB tab is only for visualization. Operations must be done on R, G, and B channels"));
			} else if (event->button == 3) {	// right click
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
					if (is_over_the_left_side_of_sel(event->x / zoom,
								event->y / zoom, zoom)) {
						com.drawing = TRUE;
						com.startX = com.selection.x + com.selection.w;
						com.freezeY = TRUE;
						com.freezeX = FALSE;
					} else if (is_over_the_right_side_of_sel(event->x / zoom,
								event->y / zoom, zoom)) {
						com.drawing = TRUE;
						com.startX = com.selection.x;
						com.freezeY = TRUE;
						com.freezeX = FALSE;
					} else if (is_over_the_bottom_of_sel(event->x / zoom,
								event->y / zoom, zoom)) {
						com.drawing = TRUE;
						com.startY = com.selection.y;
						com.freezeY = FALSE;
						com.freezeX = TRUE;
					} else if (is_over_the_top_of_sel(event->x / zoom,
								event->y / zoom, zoom)) {
						com.drawing = TRUE;
						com.startY = com.selection.y + com.selection.h;
						com.freezeY = FALSE;
						com.freezeX = TRUE;
					} else {
						com.drawing = TRUE;
						com.startX = event->x / zoom;
						com.startY = event->y / zoom;
						com.selection.h = 0;
						com.selection.w = 0;
						com.freezeX = com.freezeY = FALSE;
					}
				}
				gtk_widget_queue_draw(widget);
			} else if (mouse_status == MOUSE_ACTION_DRAW_SAMPLES) {
				double zoom = get_zoom_val();
				point pt;
				int radius = get_sample_radius();

				pt.x = (event->x / zoom) - radius;
				pt.y = (event->y / zoom) - radius;

				if (pt.x + radius <= gfit.rx && pt.y + radius <= gfit.ry
						&& pt.x - radius >= 0 && pt.y - radius >= 0) {
					com.grad_samples = add_background_sample(com.grad_samples, &gfit, pt);

					redraw(com.cvport, REMAP_NONE);
					redraw_previews();
				}
			}
		} else if (event->button == 3) {	// right click
			if (mouse_status == MOUSE_ACTION_DRAW_SAMPLES) {
				double zoom = get_zoom_val();
				point pt;
				int radius = (int) (25 / 2);

				pt.x = (event->x / zoom) - radius;
				pt.y = (event->y / zoom) - radius;

				if (pt.x + radius <= gfit.rx && pt.y + radius <= gfit.ry
						&& pt.x - radius >= 0 && pt.y - radius >= 0) {
					com.grad_samples = remove_background_sample(com.grad_samples, &gfit, pt);

					redraw(com.cvport, REMAP_NONE);
					redraw_previews();
				}
			}
		}
	}
	return FALSE;
}

gboolean on_drawingarea_button_release_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	double zoom = get_zoom_val();
	gdouble zoomedX, zoomedY;

	if (inimage((GdkEvent *) event)) {
		zoomedX = event->x / zoom;
		zoomedY = event->y / zoom;
	} else {
		if (event->x < 0)
			zoomedX = 0.0;
		else if (event->x > gfit.rx * zoom)
			zoomedX = gfit.rx;
		else
			zoomedX = event->x / zoom;
		if (event->y < 0)
			zoomedY = 0.0;
		else if (event->y > gfit.ry * zoom)
			zoomedY = gfit.ry;
		else
			zoomedY = event->y / zoom;
	}
	if (event->button == 1) {	// left click
		if (com.drawing && mouse_status == MOUSE_ACTION_SELECT_REG_AREA) {
			com.drawing = FALSE;
			/* finalize selection rectangle coordinates */
			if (!com.freezeX) {
				if (zoomedX > com.startX) {
					com.selection.x = com.startX;
					com.selection.w = zoomedX - com.selection.x;
				} else {
					com.selection.x = zoomedX;
					com.selection.w = com.startX - zoomedX;
				}
			}
			if (!com.freezeY) {
				if (zoomedY > com.startY) {
					com.selection.y = com.startY;
					com.selection.h = zoomedY - com.selection.y;
				} else {
					com.selection.y = zoomedY;
					com.selection.h = com.startY - zoomedY;
				}
			}
			/* we have a new rectangular selection zone,
			 * or an unselection (empty zone) */
			new_selection_zone();

			/* calculate and display FWHM - not in event
			 * callbacks because it's in the same file and
			 * requires a special argument */
			calculate_fwhm(widget);
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
	} else if (event->button == 2) {	// middle click
		if (inimage((GdkEvent *) event)) {
			double dX, dY, w, h;

			dX = 1.5 * com.pref.phot_set.outer;
			dY = dX;
			w = 3 * com.pref.phot_set.outer;
			h = w;

			if ((dX <= zoomedX) && (dY <= zoomedY)
					&& (zoomedX - dX + w < gfit.rx)
					&& (zoomedY - dY + h < gfit.ry)) {

				com.selection.x = zoomedX - dX;
				com.selection.y = zoomedY - dY;
				com.selection.w = w;
				com.selection.h = h;

				new_selection_zone();
				calculate_fwhm(widget);
			}
		}

	} else if (event->button == 3) {	// right click
		if (mouse_status != MOUSE_ACTION_DRAW_SAMPLES) {
			do_popup_graymenu(widget, NULL);
		}
	}
	return FALSE;
}

gboolean on_drawingarea_motion_notify_event(GtkWidget *widget,
		GdkEventMotion *event, gpointer user_data) {
	fits *fit = &(gfit);
	double zoom = get_zoom_val();
	gint zoomedX = 0, zoomedY = 0;
	const char *suffix = untranslated_vport_number_to_name(com.cvport);
	gchar *label = g_strdup_printf("labeldensity_%s", suffix);

	if (fit->type == DATA_UNSUPPORTED) return FALSE;

	if (inimage((GdkEvent *) event)) {
		char *buffer;
		char *format;
		int coords_width = 3;
		zoomedX = (gint) (event->x / zoom);
		zoomedY = (gint) (event->y / zoom);

		if (fit->rx >= 1000 || fit->ry >= 1000)
			coords_width = 4;
		if (fit->type == DATA_USHORT) {
			int val_width = 3;
			char *format_base_ushort = "x: %%.%dd y: %%.%dd = %%.%dd";
			if (fit->hi >= 1000)
				val_width = 4;
			if (fit->hi >= 10000)
				val_width = 5;
			format = g_strdup_printf(format_base_ushort,
					coords_width, coords_width, val_width);
			buffer = g_strdup_printf(format, zoomedX, zoomedY,
					fit->pdata[com.cvport][fit->rx * (fit->ry - zoomedY - 1)
					+ zoomedX]);
		} else if (fit->type == DATA_FLOAT) {
			char *format_base_float = "x: %%.%dd y: %%.%dd = %%f";
			format = g_strdup_printf(format_base_float,
					coords_width, coords_width);
			buffer = g_strdup_printf(format, zoomedX, zoomedY,
					fit->fpdata[com.cvport][fit->rx * (fit->ry - zoomedY - 1)
					+ zoomedX]);
		} else return FALSE;

		gtk_label_set_text(GTK_LABEL(lookup_widget(label)), buffer);
		g_free(buffer);
		g_free(format);
	} else {
		gtk_label_set_text(GTK_LABEL(lookup_widget(label)), "");
	}

	g_free(label);

	if (com.drawing) {	// with button 1 down
		if (!inimage((GdkEvent *) event)) {
			set_cursor("crosshair");
			if (event->x < 0)
				zoomedX = 0;
			else if (event->x > gfit.rx * zoom)
				zoomedX = gfit.rx;
			else
				zoomedX = round_to_int(event->x / zoom);
			if (event->y < 0)
				zoomedY = 0;
			else if (event->y > gfit.ry * zoom)
				zoomedY = gfit.ry;
			else
				zoomedY = round_to_int(event->y / zoom);
		}
		if (!com.freezeX) {
			if (zoomedX > com.startX) {
				com.selection.x = com.startX;
				com.selection.w = zoomedX - com.selection.x;
			} else {
				com.selection.x = zoomedX;
				com.selection.w = com.startX - zoomedX;
			}
		}
		if (!com.freezeY) {
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
		}
		gtk_widget_queue_draw(widget);
	}
	if (inimage((GdkEvent *) event)) {
		if (mouse_status == MOUSE_ACTION_DRAW_SAMPLES) {
			set_cursor("cell");
		} else {
			if (!com.drawing) {
				if (is_over_the_left_side_of_sel(zoomedX, zoomedY, zoom)) {
					set_cursor("w-resize");
				} else if (is_over_the_right_side_of_sel(zoomedX, zoomedY, zoom)) {
					set_cursor("e-resize");
				} else if (is_over_the_bottom_of_sel(zoomedX, zoomedY, zoom)) {
					set_cursor("s-resize");
				} else if (is_over_the_top_of_sel(zoomedX, zoomedY, zoom)) {
					set_cursor("n-resize");
				} else {
					set_cursor("crosshair");
				}
			}
		}
	}

	return FALSE;
}

void on_drawingarea_leave_notify_event(GtkWidget *widget, GdkEvent *event,
		gpointer user_data) {

	if (get_thread_run()) {
		set_cursor_waiting(TRUE);
	} else {
		/* trick to get default cursor */
		set_cursor_waiting(FALSE);
	}
}

static void get_scroll_position(GtkWidget *widget, int *x, int *y, int *width, int *height) {
	// the GtkScrolledWindow are the grand parents of the drawing areas
	GtkWidget *scrwindow = gtk_widget_get_parent(gtk_widget_get_parent(widget));
	GtkAdjustment *hadj = gtk_scrolled_window_get_hadjustment(GTK_SCROLLED_WINDOW(scrwindow));
	GtkAdjustment *vadj = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(scrwindow));

	*x = gtk_adjustment_get_value(hadj);
	*y = gtk_adjustment_get_value(vadj);
	*width = gtk_adjustment_get_page_size(hadj);
	*height = gtk_adjustment_get_page_size(vadj);

	siril_debug_print("get_scroll_position: %d %d %d %d\n", *x, *y, *width, *height);
}

static void set_scroll_position(GtkWidget *widget, int x, int y) {
	GtkWidget *scrwindow = gtk_widget_get_parent(gtk_widget_get_parent(widget));
	GtkAdjustment *hadj = gtk_scrolled_window_get_hadjustment(GTK_SCROLLED_WINDOW(scrwindow));
	GtkAdjustment *vadj = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(scrwindow));

	siril_debug_print("set_scroll_position: %d %d (for adj of size %g and %g)\n", x,
			y, gtk_adjustment_get_upper(hadj), gtk_adjustment_get_upper(vadj));

	x = CLAMP(x, 0, gtk_adjustment_get_upper(hadj) - gtk_adjustment_get_page_size(hadj));
	y = CLAMP(y, 0, gtk_adjustment_get_upper(vadj) - gtk_adjustment_get_page_size(vadj));

	gtk_adjustment_set_value(hadj, x);
	gtk_adjustment_set_value(vadj, y);
}

static GdkModifierType get_primary() {
	return gdk_keymap_get_modifier_mask(
			gdk_keymap_get_for_display(gdk_display_get_default()),
			GDK_MODIFIER_INTENT_PRIMARY_ACCELERATOR);
}

gboolean on_drawingarea_scroll_event(GtkWidget *widget, GdkEventScroll *event, gpointer user_data) {
	GtkToggleToolButton *button;
	gdouble delta_x, delta_y, evpos_x, evpos_y;
	gboolean handled = FALSE;
	int pix_x, pix_y, pix_width, pix_height;

	if (!single_image_is_loaded() && !sequence_is_loaded())
		return FALSE;

	if (event->state & get_primary()) {
		button = (GtkToggleToolButton *)user_data;
		if (gtk_toggle_tool_button_get_active(button))
			gtk_toggle_tool_button_set_active(button, FALSE);

		get_scroll_position(widget, &pix_x, &pix_y, &pix_width, &pix_height);
		// event position in image coordinates before changing the zoom value
		evpos_x = (event->x) / com.zoom_value;
		evpos_y = (event->y) / com.zoom_value;

		switch (event->direction) {
		case GDK_SCROLL_SMOOTH:	// what's that?
			handled = TRUE;
			gdk_event_get_scroll_deltas((GdkEvent*) event, &delta_x, &delta_y);
			if (delta_y < 0) {
				if (com.zoom_value * 1.5 > ZOOM_MAX) {
					return handled;
				}
				com.zoom_value *= 1.5;
			}
			if (delta_y > 0) {
				if (com.zoom_value / 1.5 < ZOOM_MIN) {
					return handled;
				}
				com.zoom_value /= 1.5 ;
			}
			adjust_vport_size_to_image();
			set_scroll_position(widget, evpos_x * com.zoom_value - pix_width / 2,
					evpos_y * com.zoom_value - pix_height / 2);
			redraw(com.cvport, REMAP_NONE);
			break;
		case GDK_SCROLL_DOWN:
			handled = TRUE;
			if (com.zoom_value / 1.5 < ZOOM_MIN) {
				return handled;
			}
			com.zoom_value /= 1.5 ;
			adjust_vport_size_to_image();
			// event->[xy] - pix_[xy] are the coordinates of the click in the widget
			siril_debug_print("zoom out (%f) at %f,%f in image, %f,%f on area %d,%d,%d,%d\n",
					com.zoom_value, evpos_x, evpos_y, event->x - pix_x,
					event->y - pix_y, pix_x, pix_y, pix_width, pix_height);
			// evpos_[xy] * zoom_value are the coordinates of the event in the new zoom value
			set_scroll_position(widget, evpos_x * com.zoom_value - (event->x - pix_x),
					 evpos_y * com.zoom_value - (event->y - pix_y));
			redraw(com.cvport, REMAP_NONE);
			break;
		case GDK_SCROLL_UP:
			handled = TRUE;
			if (com.zoom_value * 1.5 > ZOOM_MAX) {
				return handled;
			}
			com.zoom_value *= 1.5;
			adjust_vport_size_to_image();
			siril_debug_print("zoom in (%f) at %f,%f in image, %f,%f on area %d,%d,%d,%d\n",
					com.zoom_value, evpos_x, evpos_y, event->x - pix_x,
					event->y - pix_y, pix_x, pix_y, pix_width, pix_height);
			set_scroll_position(widget, evpos_x * com.zoom_value - (event->x - pix_x),
					 evpos_y * com.zoom_value - (event->y - pix_y));
			redraw(com.cvport, REMAP_NONE);
			break;
		default:
			handled = FALSE;
		}
	}
	return handled;
}

void on_zoom_to_fit_check_button_toggled(GtkToggleToolButton *button, gpointer data) {
	if (gtk_toggle_tool_button_get_active(button)) {
		com.zoom_value = -1;
		adjust_vport_size_to_image();
		redraw(com.cvport, REMAP_NONE);
	} else {
		com.zoom_value = get_zoom_val();
	}
}

void on_zoom_to_one_button_clicked(GtkToolButton *button, gpointer user_data) {
	GtkToggleToolButton *tbutton = (GtkToggleToolButton*) user_data;
	if (gtk_toggle_tool_button_get_active(tbutton))
		gtk_toggle_tool_button_set_active(tbutton, FALSE);
	com.zoom_value = 1;
	adjust_vport_size_to_image();
	redraw(com.cvport, REMAP_NONE);
}
