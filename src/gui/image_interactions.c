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

static gboolean is_over_the_left_side_of_sel(double zoomedX, double zoomedY, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = margin_size / zoom;
	if (zoomedX > com.selection.x - s && zoomedX < com.selection.x + s) {
		if (zoomedY > com.selection.y - s
				&& zoomedY < com.selection.y + com.selection.h + s)
			return TRUE;
	}

	return FALSE;
}

static gboolean is_over_the_right_side_of_sel(double zoomedX, double zoomedY, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = margin_size / zoom;
	if (zoomedX > com.selection.x + com.selection.w - s
				&& zoomedX < com.selection.x + com.selection.w + s) {
		if (zoomedY > com.selection.y - s
				&& zoomedY < com.selection.y + com.selection.h + s)
			return TRUE;
	}
	return FALSE;
}

static gboolean is_over_the_bottom_of_sel(double zoomedX, double zoomedY, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = margin_size / zoom;
	if (zoomedY > com.selection.y + com.selection.h - s
				&& zoomedY < com.selection.y + com.selection.h + s) {
		if (zoomedX > com.selection.x - s
				&& zoomedX < com.selection.x + com.selection.w + s)
			return TRUE;
	}
	return FALSE;
}

static gboolean is_over_the_top_of_sel(double zoomedX, double zoomedY, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = margin_size / zoom;
	if (zoomedY > com.selection.y - s && zoomedY < com.selection.y + s) {
		if (zoomedX > com.selection.x - s
				&& zoomedX < com.selection.x + com.selection.w + s)
			return TRUE;
	}
	return FALSE;
}

static gboolean is_inside_of_sel(double zoomedX, double zoomedY, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = margin_size / zoom;
	if (zoomedX >= com.selection.x + s && zoomedX <= com.selection.x + com.selection.w - s) {
		if (zoomedY >= com.selection.y + s && zoomedY <= com.selection.y + com.selection.h - s)
			return TRUE;
	}
	return FALSE;
}

/* Clamp given coordinates to image boundaries.
   Returns true if point was inside, false otherwise.
*/
static gboolean clamp2image(gint* x, gint* y) {
	gboolean x_inside = 0;
	gboolean y_inside = 0;
	if (*x < 0) {
		*x = 0;
	} else if (*x >= gfit.rx) {
		*x = gfit.rx - 1;
	} else {
		x_inside = 1;
	}

	if (*y < 0) {
		*y = 0;
	} else if (*y >= gfit.ry) {
		*y = gfit.ry - 1;
	} else {
		y_inside = 1;
	}
	return x_inside && y_inside;
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
	gtk_widget_set_sensitive(lookup_widget("menu_gray_seqpsf"), selected && sequence_is_loaded());
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
	double zoom = get_zoom_val();

	// evpos_x/evpos_y = cursor position in image coordinate
	double evpos_x = event->x;
	double evpos_y = event->y;
	cairo_matrix_transform_point(&com.image_matrix, &evpos_x, &evpos_y);

	// same as evpos but rounded to integer and clamped to image bounds
	gint zoomedX = round_to_int(evpos_x);
	gint zoomedY = round_to_int(evpos_y);
	gboolean inside = clamp2image(&zoomedX, &zoomedY);

	if (inside) {
		if ((event->state & GDK_CONTROL_MASK) == GDK_CONTROL_MASK) {
			if (event->button == GDK_BUTTON_PRIMARY) {
				// viewport translation
				com.translating = TRUE;
				com.start.x = round_to_int(event->x);
				com.start.y = round_to_int(event->y);
				return TRUE;
			}
		}

		/* click on RGB image */
		if (widget == com.vport[RGB_VPORT]) {
			if (event->button == GDK_BUTTON_PRIMARY) {	// left click
				siril_message_dialog(GTK_MESSAGE_INFO, _("Only for visualization"),
						_("The RGB tab is only for visualization. Operations must be done on R, G, and B channels"));
			} else if (event->button == GDK_BUTTON_SECONDARY) {	// right click
				do_popup_rgbmenu(widget, event);
				return TRUE;
			}
			return FALSE;
		}

		/* else, click on gray image */
		if (event->button == GDK_BUTTON_PRIMARY) {	// left click
			if (mouse_status == MOUSE_ACTION_SELECT_REG_AREA) {
				if (com.drawing) {
					com.drawing = FALSE;
				} else {
					com.drawing = TRUE;
					if (is_inside_of_sel(zoomedX, zoomedY, zoom)) {
						// Move selection
						com.freezeX = com.freezeY = TRUE;
						com.start.x = zoomedX;
						com.start.y = zoomedY;
						com.origin.x = com.selection.x;
						com.origin.y = com.selection.y;
					} else {
						// Default values
						com.freezeX = com.freezeY = FALSE;
						// The order matters if the selection is so small that edge detection overlaps
						// and need to be the same as in the on_drawingarea_motion_notify_event()
						gboolean right = is_over_the_right_side_of_sel(zoomedX, zoomedY, zoom);
						gboolean left = is_over_the_left_side_of_sel(zoomedX, zoomedY, zoom);
						gboolean bottom = is_over_the_bottom_of_sel(zoomedX, zoomedY, zoom);
						gboolean top = is_over_the_top_of_sel(zoomedX, zoomedY, zoom);
						if (right || left || bottom || top) {
							// Freeze one axis when grabbing an edges far enough from a corner
							if (right) {
								com.start.x = com.selection.x;
								if (!bottom && !top)
									com.freezeY = TRUE;
							} else if (left) {
								com.start.x = com.selection.x + com.selection.w;
								if (!bottom && !top)
									com.freezeY = TRUE;
							}
							if (bottom) {
								com.start.y = com.selection.y;
								if (!left && !right)
									com.freezeX = TRUE;
							} else if (top) {
								com.start.y = com.selection.y + com.selection.h;
								if (!left && !right)
									com.freezeX = TRUE;
							}
						} else {
							com.start.x = zoomedX;
							com.start.y = zoomedY;
							com.selection.h = 0;
							com.selection.w = 0;
						}
					}
				}
				gtk_widget_queue_draw(widget);
			} else if (mouse_status == MOUSE_ACTION_DRAW_SAMPLES) {
				point pt;
				int radius = get_sample_radius();

				pt.x = zoomedX - radius;
				pt.y = zoomedY - radius;

				if (pt.x + radius <= gfit.rx && pt.y + radius <= gfit.ry
						&& pt.x - radius >= 0 && pt.y - radius >= 0) {
					com.grad_samples = add_background_sample(com.grad_samples, &gfit, pt);

					redraw(com.cvport, REMAP_NONE);
					redraw_previews();
				}
			}
		} else if (event->button == GDK_BUTTON_SECONDARY) {	// right click
			if (mouse_status == MOUSE_ACTION_DRAW_SAMPLES) {
				point pt;
				int radius = (int) (25 / 2);

				pt.x = zoomedX - radius;
				pt.y = zoomedY - radius;

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
	// evpos_x/evpos_y = cursor position in image coordinate
	double evpos_x = event->x;
	double evpos_y = event->y;
	cairo_matrix_transform_point(&com.image_matrix, &evpos_x, &evpos_y);

	// same as evpos but rounded to integer and clamped to image bounds
	gint zoomedX = round_to_int(evpos_x);
	gint zoomedY = round_to_int(evpos_y);
	gboolean inside = clamp2image(&zoomedX, &zoomedY);

	if (event->button == GDK_BUTTON_PRIMARY) {	// left click
		if (com.translating) {
			com.translating = FALSE;
		} else if (com.drawing && mouse_status == MOUSE_ACTION_SELECT_REG_AREA) {
			com.drawing = FALSE;
			/* finalize selection rectangle coordinates */
			if (!com.freezeX) {
				if (zoomedX > com.start.x) {
					com.selection.x = com.start.x;
					com.selection.w = zoomedX - com.selection.x;
				} else {
					com.selection.x = zoomedX;
					com.selection.w = com.start.x - zoomedX;
				}
			}
			if (!com.freezeY) {
				if (zoomedY > com.start.y) {
					com.selection.y = com.start.y;
					com.selection.h = zoomedY - com.selection.y;
				} else {
					com.selection.y = zoomedY;
					com.selection.h = com.start.y - zoomedY;
				}
			}
			if (com.freezeX && com.freezeY) {
				// Move selection
				com.selection.x = (zoomedX - com.start.x) + com.origin.x;
				com.selection.y = (zoomedY - com.start.y) + com.origin.y;
				// clamp it inside the image
				if (com.selection.x < 0)
					com.selection.x = 0;
				else if (com.selection.x + com.selection.w > gfit.rx)
					com.selection.x = gfit.rx - com.selection.w;
				if (com.selection.y < 0)
					com.selection.y = 0;
				else if (com.selection.y + com.selection.h > gfit.ry)
					com.selection.y = gfit.ry - com.selection.h;
			}
			/* we have a new rectangular selection zone,
			 * or an unselection (empty zone) */
			new_selection_zone();
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
	} else if (event->button == GDK_BUTTON_MIDDLE) {	// middle click
		if (inside) {
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
			}
		}

	} else if (event->button == GDK_BUTTON_SECONDARY) {	// right click
		if (mouse_status != MOUSE_ACTION_DRAW_SAMPLES) {
			do_popup_graymenu(widget, NULL);
		}
	}
	return FALSE;
}

static GdkModifierType get_primary() {
	return gdk_keymap_get_modifier_mask(
			gdk_keymap_get_for_display(gdk_display_get_default()),
			GDK_MODIFIER_INTENT_PRIMARY_ACCELERATOR);
}

gboolean on_drawingarea_motion_notify_event(GtkWidget *widget,
		GdkEventMotion *event, gpointer user_data) {
	if (gfit.type == DATA_UNSUPPORTED) return FALSE;

	double zoom = get_zoom_val();

	// evpos_x/evpos_y = cursor position in image coordinate
	double evpos_x = event->x;
	double evpos_y = event->y;
	cairo_matrix_transform_point(&com.image_matrix, &evpos_x, &evpos_y);

	// same as evpos but rounded to integer and clamped to image bounds
	gint zoomedX = round_to_int(evpos_x);
	gint zoomedY = round_to_int(evpos_y);
	gboolean inside = clamp2image(&zoomedX, &zoomedY);

	const char *suffix = untranslated_vport_number_to_name(com.cvport);
	gchar *label = g_strdup_printf("labeldensity_%s", suffix);

	if (inside && com.cvport < RGB_VPORT) {
		char *buffer = NULL;
		char *format = NULL;
		int coords_width = 3;

		if (gfit.rx >= 1000 || gfit.ry >= 1000)
			coords_width = 4;
		if (gfit.type == DATA_USHORT && gfit.pdata[com.cvport] != NULL) {
			int val_width = 3;
			char *format_base_ushort = "x: %%.%dd y: %%.%dd = %%.%dd";
			if (gfit.hi >= 1000)
				val_width = 4;
			if (gfit.hi >= 10000)
				val_width = 5;
			format = g_strdup_printf(format_base_ushort,
					coords_width, coords_width, val_width);
			buffer = g_strdup_printf(format, zoomedX, zoomedY,
					gfit.pdata[com.cvport][gfit.rx * (gfit.ry - zoomedY - 1)
					+ zoomedX]);
		} else if (gfit.type == DATA_FLOAT  && gfit.fpdata[com.cvport] != NULL) {
			char *format_base_float = "x: %%.%dd y: %%.%dd = %%f";
			format = g_strdup_printf(format_base_float,
					coords_width, coords_width);
			buffer = g_strdup_printf(format, zoomedX, zoomedY,
					gfit.fpdata[com.cvport][gfit.rx * (gfit.ry - zoomedY - 1)
					+ zoomedX]);
		}

		if (buffer) {
			gtk_label_set_text(GTK_LABEL(lookup_widget(label)), buffer);
			g_free(buffer);
			g_free(format);
		}
	} else if (widget != com.vport[RGB_VPORT]) {
		gtk_label_set_text(GTK_LABEL(lookup_widget(label)), "");
	}

	g_free(label);

	if (com.translating) {
		gint ev_x = round_to_int(event->x);
		gint ev_y = round_to_int(event->y);
		double delta_x = ev_x - com.start.x; com.start.x = ev_x;
		double delta_y = ev_y - com.start.y; com.start.y = ev_y;
		com.display_offset_x += delta_x;
		com.display_offset_y += delta_y;
		adjust_vport_size_to_image();
		gtk_widget_queue_draw(widget);
	} else if (com.drawing) {	// with button 1 down
		if (!com.freezeX) {
			if (zoomedX > com.start.x) {
				com.selection.x = com.start.x;
				com.selection.w = zoomedX - com.selection.x;
			} else {
				com.selection.x = zoomedX;
				com.selection.w = com.start.x - zoomedX;
			}
		}
		if (!com.freezeY) {
			if (zoomedY > com.start.y) {
				com.selection.y = com.start.y;
				com.selection.h = zoomedY - com.selection.y;
			} else {
				com.selection.y = zoomedY;
				com.selection.h = com.start.y - zoomedY;
			}
		}
		if (com.freezeX && com.freezeY) {
			// Move selection
			com.selection.x = (zoomedX - com.start.x) + com.origin.x;
			com.selection.y = (zoomedY - com.start.y) + com.origin.y;
			// clamp it inside the image
			if (com.selection.x < 0)
				com.selection.x = 0;
			else if (com.selection.x + com.selection.w > gfit.rx)
				com.selection.x = gfit.rx - com.selection.w;
			if (com.selection.y < 0)
				com.selection.y = 0;
			else if (com.selection.y + com.selection.h > gfit.ry)
				com.selection.y = gfit.ry - com.selection.h;
		}
		// Display the dimensions of the selection while drawing it
		update_display_selection();
		gtk_widget_queue_draw(widget);
	}

	if (inside && com.cvport < RGB_VPORT) {
		if (mouse_status == MOUSE_ACTION_DRAW_SAMPLES) {
			set_cursor("cell");
		} else {
			if (!com.drawing && !com.translating) {
				// The order matters if the selection is so small that edge detection overlaps
				// and need to be the same as in the on_drawingarea_button_press_event()
				gboolean right = is_over_the_right_side_of_sel(zoomedX, zoomedY, zoom);
				gboolean left = is_over_the_left_side_of_sel(zoomedX, zoomedY, zoom);
				gboolean bottom = is_over_the_bottom_of_sel(zoomedX, zoomedY, zoom);
				gboolean top = is_over_the_top_of_sel(zoomedX, zoomedY, zoom);
				if (bottom && right) {
					set_cursor("se-resize");
				} else if (top && right) {
					set_cursor("ne-resize");
				} else if (right) {
					set_cursor("e-resize");
				} else if (bottom && left) {
					set_cursor("sw-resize");
				} else if (top && left) {
					set_cursor("nw-resize");
				} else if (left) {
					set_cursor("w-resize");
				} else if (bottom) {
					set_cursor("s-resize");
				} else if (top) {
					set_cursor("n-resize");
				} else if (is_inside_of_sel(zoomedX, zoomedY, zoom)) {
					set_cursor("all-scroll");
				} else {
					if (event->state & get_primary()) {
						set_cursor("all-scroll");
					} else {
						set_cursor("crosshair");
					}
				}
			} else {
				if (event->state & get_primary()) {
					set_cursor("all-scroll");
				} else {
					set_cursor("crosshair");
				}
			}
		}
	} else {
		set_cursor("default");
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

gboolean on_drawingarea_scroll_event(GtkWidget *widget, GdkEventScroll *event, gpointer user_data) {
	gboolean handled = FALSE;

	if (!single_image_is_loaded() && !sequence_is_loaded())
		return FALSE;

	if (event->state & get_primary()) {
		GtkToggleToolButton *button = (GtkToggleToolButton *)user_data;
		if (gtk_toggle_tool_button_get_active(button))
			gtk_toggle_tool_button_set_active(button, FALSE);

		// event position in image coordinates before changing the zoom value
		gdouble evpos_x = event->x;
		gdouble evpos_y = event->y;
		cairo_matrix_transform_point(&com.image_matrix, &evpos_x, &evpos_y);

		gdouble delta_x, delta_y, factor;

		switch (event->direction) {
		case GDK_SCROLL_SMOOTH:	// what's that?
			handled = TRUE;
			gdk_event_get_scroll_deltas((GdkEvent*) event, &delta_x, &delta_y);
			if (delta_y < 0) {
				if (com.zoom_value * 1.5 > ZOOM_MAX) {
					return handled;
				}
				factor = 1.5;
			}
			if (delta_y > 0) {
				if (com.zoom_value / 1.5 < ZOOM_MIN) {
					return handled;
				}
				factor = 1. / 1.5;
			}
			break;
		case GDK_SCROLL_DOWN:
			handled = TRUE;
			if (com.zoom_value / 1.5 < ZOOM_MIN) {
				return handled;
			}
			factor = 1. / 1.5;
			break;
		case GDK_SCROLL_UP:
			handled = TRUE;
			if (com.zoom_value * 1.5 > ZOOM_MAX) {
				return handled;
			}
			factor = 1.5;
			break;
		default:
			handled = FALSE;
		}

		if (handled) {
			com.zoom_value *= factor;
			adjust_vport_size_to_image();
			cairo_matrix_transform_point(&com.display_matrix, &evpos_x, &evpos_y);
			com.display_offset_x += event->x - evpos_x;
			com.display_offset_y += event->y - evpos_y;
			adjust_vport_size_to_image();
			redraw(com.cvport, REMAP_NONE);
		}
	}
	return handled;
}

void on_zoom_to_fit_check_button_toggled(GtkToggleToolButton *button, gpointer data) {
	if (gtk_toggle_tool_button_get_active(button)) {
		com.zoom_value = -1;
		com.display_offset_x = 0;
		com.display_offset_y = 0;
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
