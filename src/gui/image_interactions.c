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

static gboolean is_over_the_left_side_of_sel(pointi zoomed, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	int s = round_to_int(margin_size / zoom);
	if (zoomed.x > com.selection.x - s && zoomed.x < com.selection.x + s) {
		if (zoomed.y > com.selection.y - s
				&& zoomed.y < com.selection.y + com.selection.h + s)
			return TRUE;
	}

	return FALSE;
}

static gboolean is_over_the_right_side_of_sel(pointi zoomed, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = margin_size / zoom;
	if (zoomed.x > com.selection.x + com.selection.w - s
				&& zoomed.x < com.selection.x + com.selection.w + s) {
		if (zoomed.y > com.selection.y - s
				&& zoomed.y < com.selection.y + com.selection.h + s)
			return TRUE;
	}
	return FALSE;
}

static gboolean is_over_the_bottom_of_sel(pointi zoomed, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = margin_size / zoom;
	if (zoomed.y > com.selection.y + com.selection.h - s
				&& zoomed.y < com.selection.y + com.selection.h + s) {
		if (zoomed.x > com.selection.x - s
				&& zoomed.x < com.selection.x + com.selection.w + s)
			return TRUE;
	}
	return FALSE;
}

static gboolean is_over_the_top_of_sel(pointi zoomed, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = margin_size / zoom;
	if (zoomed.y > com.selection.y - s && zoomed.y < com.selection.y + s) {
		if (zoomed.x > com.selection.x - s
				&& zoomed.x < com.selection.x + com.selection.w + s)
			return TRUE;
	}
	return FALSE;
}

static gboolean is_inside_of_sel(pointi zoomed, double zoom) {
	if (com.selection.w == 0 && com.selection.h == 0) return FALSE;
	double s = margin_size / zoom;
	if (zoomed.x >= com.selection.x + s && zoomed.x <= com.selection.x + com.selection.w - s) {
		if (zoomed.y >= com.selection.y + s && zoomed.y <= com.selection.y + com.selection.h - s)
			return TRUE;
	}
	return FALSE;
}

/* Clamp given coordinates to image boundaries.
   Returns true if point was inside, false otherwise.
*/
static gboolean clamp2image(pointi* pt) {
	gboolean x_inside = 0;
	gboolean y_inside = 0;
	if (pt->x < 0) {
		pt->x = 0;
	} else if (pt->x >= gfit.rx) {
		pt->x = gfit.rx - 1;
	} else {
		x_inside = 1;
	}

	if (pt->y < 0) {
		pt->y = 0;
	} else if (pt->y >= gfit.ry) {
		pt->y = gfit.ry - 1;
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
	siril_debug_print("selection: %d,%d,\t%dx%d (%lf)\n", com.selection.x, com.selection.y,
			com.selection.w, com.selection.h, com.ratio);
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

	// selection submenu
	double original_ratio = (double)gfit.rx / (double)gfit.ry;
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_free")), com.ratio == 0.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_preserve")), com.ratio == original_ratio);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_16_9")), com.ratio == 16.0 / 9.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_4_3")), com.ratio == 4.0 / 3.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_1_1")), com.ratio == 1.0 / 1.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_3_4")), com.ratio == 3.0 / 4.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_9_16")), com.ratio == 9.0 / 16.0);
	gtk_widget_set_sensitive(lookup_widget("menuitem_selection_preserve"), is_a_single_image_loaded || sequence_is_loaded());
	gtk_widget_set_sensitive(lookup_widget("menuitem_selection_all"), is_a_single_image_loaded || sequence_is_loaded());

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


static void enforce_ratio_and_clamp() {
	if (com.ratio > 0.0
		&& !(com.freezeX && com.freezeY)) {
		// Enforce a ratio for the selection
		if (com.freezeY) {
			com.selection.h = round_to_int(com.selection.w / com.ratio);
		} else if (com.freezeX) {
			com.selection.w = round_to_int(com.selection.h * com.ratio);
		} else {
			gint delta_w = round_to_int(com.selection.h * com.ratio) - com.selection.w;
			com.selection.w += delta_w;

			if (com.selection.x < com.start.x) { // Changing selection from the left
				com.selection.x -= delta_w;
			}
		}

		// clamp the selection dimensions
		if (com.selection.w > gfit.rx) {
			com.selection.w = gfit.rx;
			com.selection.h = round_to_int(com.selection.w / com.ratio);
		}
		else if (com.selection.h > gfit.ry) {
			com.selection.h = gfit.ry;
			com.selection.w = round_to_int(com.selection.h * com.ratio);
		}
	}

	// clamp the selection inside the image (needed when enforcing a ratio or moving)
	com.selection.x = set_int_in_interval(com.selection.x, 0, gfit.rx - com.selection.w);
	com.selection.y = set_int_in_interval(com.selection.y, 0, gfit.ry - com.selection.h);
}

static void set_selection_ratio(double ratio) {
	com.ratio = ratio;
	enforce_ratio_and_clamp();
	update_display_selection();
	new_selection_zone();
	gtk_widget_queue_draw(com.vport[com.cvport]);
}

void on_menuitem_selection_free_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		com.ratio = 0.0;
	}
}

void on_menuitem_selection_preserve_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		set_selection_ratio((double)gfit.rx / (double)gfit.ry);
	}
}

void on_menuitem_selection_16_9_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		set_selection_ratio(16.0 / 9.0);
	}
}

void on_menuitem_selection_4_3_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		set_selection_ratio(4.0 / 3.0);
	}
}

void on_menuitem_selection_1_1_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		set_selection_ratio(1.0 / 1.0);
	}
}

void on_menuitem_selection_3_4_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		set_selection_ratio(3.0 / 4.0);
	}
}

void on_menuitem_selection_9_16_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		set_selection_ratio(9.0 / 16.0);
	}
}

void on_menuitem_selection_all_activate(GtkMenuItem *menuitem, gpointer user_data) {
	com.selection.x = 0;
	com.selection.y = 0;
	com.selection.w = gfit.rx;
	com.selection.h = gfit.ry;
	// "Select All" need to reset any enforced ratio that would not match the ratio of the image
	// 1. it's nice to NOT enforce a ratio when the user just want to select the whole image
	// 2. it's nice to keep the enforced ratio if it does match the image
	if (com.ratio != ((double)gfit.rx / (double)gfit.ry)) {
		set_selection_ratio(0.0);
	} else {
		set_selection_ratio((double)gfit.rx / (double)gfit.ry); // triggers the new_selection() callbacks etc.
	}
}

gboolean rgb_area_popup_menu_handler(GtkWidget *widget) {
	do_popup_rgbmenu(widget, NULL);
	return TRUE;
}

gboolean on_drawingarea_button_press_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {
	double zoom = get_zoom_val();

	// evpos.x/evpos.y = cursor position in image coordinate
	point evpos = { event->x, event->y };
	cairo_matrix_transform_point(&com.image_matrix, &evpos.x, &evpos.y);

	// same as evpos but rounded to integer and clamped to image bounds
	pointi zoomed = { round_to_int(evpos.x), round_to_int(evpos.y) };
	gboolean inside = clamp2image(&zoomed);

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
					if (is_inside_of_sel(zoomed, zoom)) {
						// Move selection
						com.freezeX = com.freezeY = TRUE;
						com.start = zoomed;
						com.origin.x = com.selection.x;
						com.origin.y = com.selection.y;
					} else {
						// Default values
						com.freezeX = com.freezeY = FALSE;
						// The order matters if the selection is so small that edge detection overlaps
						// and need to be the same as in the on_drawingarea_motion_notify_event()
						gboolean right = is_over_the_right_side_of_sel(zoomed, zoom);
						gboolean left = is_over_the_left_side_of_sel(zoomed, zoom);
						gboolean bottom = is_over_the_bottom_of_sel(zoomed, zoom);
						gboolean top = is_over_the_top_of_sel(zoomed, zoom);
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
							com.start = zoomed;
							com.selection.h = 0;
							com.selection.w = 0;
						}
					}
				}
				gtk_widget_queue_draw(widget);
			} else if (mouse_status == MOUSE_ACTION_DRAW_SAMPLES) {
				point pt;
				int radius = get_sample_radius();

				pt.x = zoomed.x - radius;
				pt.y = zoomed.y - radius;

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

				pt.x = zoomed.x - radius;
				pt.y = zoomed.y - radius;

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
	// evpos.x/evpos.y = cursor position in image coordinate
	point evpos = { event->x, event->y };
	cairo_matrix_transform_point(&com.image_matrix, &evpos.x, &evpos.y);

	// same as evpos but rounded to integer and clamped to image bounds
	pointi zoomed = { round_to_int(evpos.x), round_to_int(evpos.y) };
	gboolean inside = clamp2image(&zoomed);

	if (event->button == GDK_BUTTON_PRIMARY) {	// left click
		if (com.translating) {
			com.translating = FALSE;
		} else if (com.drawing && mouse_status == MOUSE_ACTION_SELECT_REG_AREA) {
			com.drawing = FALSE;
			/* finalize selection rectangle coordinates */
			if (!com.freezeX) {
				if (zoomed.x > com.start.x) {
					com.selection.x = com.start.x;
					com.selection.w = zoomed.x - com.selection.x;
				} else {
					com.selection.x = zoomed.x;
					com.selection.w = com.start.x - zoomed.x;
				}
			}
			if (!com.freezeY) {
				if (zoomed.y > com.start.y) {
					com.selection.y = com.start.y;
					com.selection.h = zoomed.y - com.selection.y;
				} else {
					com.selection.y = zoomed.y;
					com.selection.h = com.start.y - zoomed.y;
				}
			}

			if (com.freezeX && com.freezeY) { // Move selection
				com.selection.x = (zoomed.x - com.start.x) + com.origin.x;
				com.selection.y = (zoomed.y - com.start.y) + com.origin.y;
			}

			// Enforce a ratio and clamp selection to the image
			enforce_ratio_and_clamp();

			/* we have a new rectangular selection zone,
			 * or an unselection (empty zone) */
			new_selection_zone();

			// Terminate any specific selection modification mode
			com.freezeX = com.freezeY = FALSE;
		} else if (mouse_status == MOUSE_ACTION_SELECT_PREVIEW1) {
			set_preview_area(0, zoomed.x, zoomed.y);
			mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
			// redraw to get the position of the new preview area
			gtk_widget_queue_draw(widget);
		} else if (mouse_status == MOUSE_ACTION_SELECT_PREVIEW2) {
			set_preview_area(1, zoomed.x, zoomed.y);
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

			if ((dX <= zoomed.x) && (dY <= zoomed.y)
					&& (zoomed.x - dX + w < gfit.rx)
					&& (zoomed.y - dY + h < gfit.ry)) {

				com.selection.x = zoomed.x - dX;
				com.selection.y = zoomed.y - dY;
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

	// evpos.x/evpos.y = cursor position in image coordinate
	point evpos = { event->x, event->y };
	cairo_matrix_transform_point(&com.image_matrix, &evpos.x, &evpos.y);

	// same as evpos but rounded to integer and clamped to image bounds
	pointi zoomed = { round_to_int(evpos.x), round_to_int(evpos.y) };
	gboolean inside = clamp2image(&zoomed);

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
			buffer = g_strdup_printf(format, zoomed.x, zoomed.y,
					gfit.pdata[com.cvport][gfit.rx * (gfit.ry - zoomed.y - 1)
					+ zoomed.x]);
		} else if (gfit.type == DATA_FLOAT  && gfit.fpdata[com.cvport] != NULL) {
			char *format_base_float = "x: %%.%dd y: %%.%dd = %%f";
			format = g_strdup_printf(format_base_float,
					coords_width, coords_width);
			buffer = g_strdup_printf(format, zoomed.x, zoomed.y,
					gfit.fpdata[com.cvport][gfit.rx * (gfit.ry - zoomed.y - 1)
					+ zoomed.x]);
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
		pointi ev = { round_to_int(event->x), round_to_int(event->y) };
		point delta = { ev.x - com.start.x , ev.y - com.start.y };
		com.start = ev;
		com.display_offset.x += delta.x;
		com.display_offset.y += delta.y;
		adjust_vport_size_to_image();
		gtk_widget_queue_draw(widget);
	} else if (com.drawing) {	// with button 1 down
		if (!com.freezeX) {
			if (zoomed.x > com.start.x) {
				com.selection.x = com.start.x;
				com.selection.w = zoomed.x - com.selection.x;
			} else {
				com.selection.x = zoomed.x;
				com.selection.w = com.start.x - zoomed.x;
			}
		}
		if (!com.freezeY) {
			if (zoomed.y > com.start.y) {
				com.selection.y = com.start.y;
				com.selection.h = zoomed.y - com.selection.y;
			} else {
				com.selection.y = zoomed.y;
				com.selection.h = com.start.y - zoomed.y;
			}
		}

		if (com.freezeX && com.freezeY) { // Move selection
			com.selection.x = (zoomed.x - com.start.x) + com.origin.x;
			com.selection.y = (zoomed.y - com.start.y) + com.origin.y;
		}

		// Enforce a ratio and clamp selection to the image
		enforce_ratio_and_clamp();

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
				gboolean right = is_over_the_right_side_of_sel(zoomed, zoom);
				gboolean left = is_over_the_left_side_of_sel(zoomed, zoom);
				gboolean bottom = is_over_the_bottom_of_sel(zoomed, zoom);
				gboolean top = is_over_the_top_of_sel(zoomed, zoom);
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
				} else if (is_inside_of_sel(zoomed, zoom)) {
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
		point evpos = { event->x, event->y };
		cairo_matrix_transform_point(&com.image_matrix, &evpos.x, &evpos.y);

		point delta;
		gdouble factor;

		switch (event->direction) {
		case GDK_SCROLL_SMOOTH:	// what's that?
			handled = TRUE;
			gdk_event_get_scroll_deltas((GdkEvent*) event, &delta.x, &delta.y);
			if (delta.y < 0) {
				if (com.zoom_value * 1.5 > ZOOM_MAX) {
					return handled;
				}
				factor = 1.5;
			}
			if (delta.y > 0) {
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
			cairo_matrix_transform_point(&com.display_matrix, &evpos.x, &evpos.y);
			com.display_offset.x += event->x - evpos.x;
			com.display_offset.y += event->y - evpos.y;
			adjust_vport_size_to_image();
			redraw(com.cvport, REMAP_NONE);
		}
	}
	return handled;
}

void on_zoom_to_fit_check_button_toggled(GtkToggleToolButton *button, gpointer data) {
	if (gtk_toggle_tool_button_get_active(button)) {
		com.zoom_value = -1;
		com.display_offset.x = 0;
		com.display_offset.y = 0;
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
	com.display_offset.x = 0;
	com.display_offset.y = 0;
	adjust_vport_size_to_image();
	redraw(com.cvport, REMAP_NONE);
}
