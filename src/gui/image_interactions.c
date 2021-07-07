/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2021 team free-astro (see more in AUTHORS file)
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

#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/undo.h"
#include "core/siril_world_cs.h"
#include "algos/background_extraction.h"
#include "algos/siril_wcs.h"
#include "io/single_image.h"
#include "io/sequence.h"
#include "gui/open_dialog.h"
#include "gui/PSF_list.h"
#include "image_interactions.h"
#include "image_display.h"
#include "gui/callbacks.h"
#include "gui/utils.h"
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

#define MAX_CALLBACKS_PER_EVENT 10
/* selection zone event management */
static selection_update_callback _registered_selection_callbacks[MAX_CALLBACKS_PER_EVENT];
static int _nb_selection_callbacks = 0;

void register_selection_update_callback(selection_update_callback f) {
	if (_nb_selection_callbacks < MAX_CALLBACKS_PER_EVENT) {
		_registered_selection_callbacks[_nb_selection_callbacks] = f;
		_nb_selection_callbacks++;
	}
}

void unregister_selection_update_callback(selection_update_callback f) {
	int i;
	for (i = 0; i < _nb_selection_callbacks; ++i) {
		if (_registered_selection_callbacks[i] == f) {
			_registered_selection_callbacks[i] =
				_registered_selection_callbacks[_nb_selection_callbacks];
			_registered_selection_callbacks[_nb_selection_callbacks] = NULL;
			_nb_selection_callbacks--;
			return;
		}
	}
}

// send the events
static void new_selection_zone() {
	int i;
	siril_debug_print("selection: %d,%d,\t%dx%d (%lf)\n", com.selection.x, com.selection.y,
			com.selection.w, com.selection.h, com.ratio);
	for (i = 0; i < _nb_selection_callbacks; ++i) {
		_registered_selection_callbacks[i]();
	}
	redraw(com.cvport, REMAP_NONE);
}


void delete_selected_area() {
	memset(&com.selection, 0, sizeof(rectangle));
	new_selection_zone();
}

void reset_display_offset() {
	com.display_offset.x = 0;
	com.display_offset.y = 0;
}

void reset_zoom_default() {
	com.zoom_value = ZOOM_DEFAULT;
	if (com.zoom_value == ZOOM_FIT && !com.script) {
		gtk_toggle_tool_button_set_active(GTK_TOGGLE_TOOL_BUTTON(lookup_widget("zoom_to_fit_check_button")), TRUE);
	}
}

void update_zoom_fit_button() {
	GtkToggleToolButton *button = GTK_TOGGLE_TOOL_BUTTON(lookup_widget("zoom_to_fit_check_button"));
	if (gtk_toggle_tool_button_get_active(button)) {
		gtk_toggle_tool_button_set_active(button, FALSE);
	}
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
#ifdef HAVE_WCSLIB
	gtk_widget_set_sensitive(lookup_widget("menu_gray_search"), has_wcs(&gfit));
#else
	gtk_widget_set_sensitive(lookup_widget("menu_gray_search"), FALSE);
#endif

	// selection submenu
	double original_ratio = (double)gfit.rx / (double)gfit.ry;
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_free")), com.ratio == 0.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_preserve")), com.ratio == original_ratio);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_16_9")), com.ratio == 16.0 / 9.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_4_3")), com.ratio == 4.0 / 3.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_3_2")), com.ratio == 3.0 / 2.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_1_1")), com.ratio == 1.0 / 1.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_3_4")), com.ratio == 3.0 / 4.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_2_3")), com.ratio == 2.0 / 3.0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_9_16")), com.ratio == 9.0 / 16.0);
	gtk_widget_set_sensitive(lookup_widget("menuitem_selection_preserve"), is_a_single_image_loaded || sequence_is_loaded());
	gtk_widget_set_sensitive(lookup_widget("menuitem_selection_all"), is_a_single_image_loaded || sequence_is_loaded());
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_guides_0")), com.pref.selection_guides == 0);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_guides_2")), com.pref.selection_guides == 2);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_guides_3")), com.pref.selection_guides == 3);
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(lookup_widget("menuitem_selection_guides_5")), com.pref.selection_guides == 5);

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

void on_menuitem_selection_3_2_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		set_selection_ratio(3.0 / 2.0);
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

void on_menuitem_selection_2_3_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		set_selection_ratio(2.0 / 3.0);
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

void menuitem_selection_guides_0_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		com.pref.selection_guides = 0;
	}
}

void menuitem_selection_guides_2_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		com.pref.selection_guides = 2;
	}
}

void menuitem_selection_guides_3_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		com.pref.selection_guides = 3;
	}
}

void menuitem_selection_guides_5_toggled(GtkCheckMenuItem *menuitem, gpointer user_data) {
	if (gtk_check_menu_item_get_active(menuitem)) {
		com.pref.selection_guides = 5;
	}
}

gboolean rgb_area_popup_menu_handler(GtkWidget *widget) {
	do_popup_rgbmenu(widget, NULL);
	return TRUE;
}

static GdkModifierType get_primary() {
	return gdk_keymap_get_modifier_mask(
			gdk_keymap_get_for_display(gdk_display_get_default()),
			GDK_MODIFIER_INTENT_PRIMARY_ACCELERATOR);
}

gboolean on_drawingarea_button_press_event(GtkWidget *widget,
		GdkEventButton *event, gpointer user_data) {

	/* when double clicking on drawing area  (if no images loaded)
	 * you can load an image This feature is in GIMP and I really
	 * love it: lazy world :).
	 */
	if (!single_image_is_loaded() && !sequence_is_loaded()) {
		if (event->button == GDK_BUTTON_PRIMARY
				&& event->type == GDK_DOUBLE_BUTTON_PRESS) {
			header_open_button_clicked();
		}
		return FALSE;
	}

	double zoom = get_zoom_val();

	// evpos.x/evpos.y = cursor position in image coordinate
	point evpos = { event->x, event->y };
	cairo_matrix_transform_point(&com.image_matrix, &evpos.x, &evpos.y);

	// same as evpos but rounded to integer and clamped to image bounds
	pointi zoomed = { (int)(evpos.x), (int)(evpos.y) };
	gboolean inside = clamp2image(&zoomed);

	if (inside) {
		if (event->state & get_primary()) {
			if (event->button == GDK_BUTTON_PRIMARY) {
				// viewport translation
				com.translating = TRUE;
				com.start.x = (int)(event->x);
				com.start.y = (int)(event->y);
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

				pt.x = (gdouble) zoomed.x;
				pt.y = (gdouble) zoomed.y;

				if (pt.x + radius <= gfit.rx && pt.y + radius <= gfit.ry
						&& pt.x - radius >= 0 && pt.y - radius >= 0) {
					com.grad_samples = add_background_sample(com.grad_samples, &gfit, pt);

					redraw(com.cvport, REMAP_NONE);
					redraw_previews();
				}
			} else if (mouse_status == MOUSE_ACTION_PHOTOMETRY) {
				int s = com.pref.phot_set.outer;
				rectangle area = { zoomed.x - s, zoomed.y - s, s * 2, s * 2 };
				if (area.x - s > 0 && area.x + s < gfit.rx
						&& area.y - s > 0 && area.y + s < gfit.ry) {
					com.qphot = psf_get_minimisation(&gfit, com.cvport, &area, TRUE, TRUE, TRUE);
					if (com.qphot) {
						com.qphot->xpos = com.qphot->x0 + area.x;
						if (gfit.top_down)
							com.qphot->ypos = com.qphot->y0 + area.y;
						else
							com.qphot->ypos = area.y + area.h - com.qphot->y0;
						redraw(com.cvport, REMAP_NONE);
						popup_psf_result(com.qphot);
					}
				}
			}
		} else if (event->button == GDK_BUTTON_SECONDARY) {	// right click
			if (mouse_status == MOUSE_ACTION_DRAW_SAMPLES) {
				point pt;
				int radius = (int) (25 / 2);

				pt.x = (gdouble) zoomed.x;
				pt.y = (gdouble) zoomed.y;

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
	pointi zoomed = { (int)(evpos.x), (int)(evpos.y) };
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

gboolean on_drawingarea_motion_notify_event(GtkWidget *widget,
		GdkEventMotion *event, gpointer user_data) {
	if (gfit.type == DATA_UNSUPPORTED) return FALSE;

	double zoom = get_zoom_val();

	// evpos.x/evpos.y = cursor position in image coordinate
	point evpos = { event->x, event->y };
	cairo_matrix_transform_point(&com.image_matrix, &evpos.x, &evpos.y);

	// same as evpos but rounded to integer and clamped to image bounds
	pointi zoomed = { (int)(evpos.x), (int)(evpos.y) };
	gboolean inside = clamp2image(&zoomed);

	static const gchar *label_density[] = { "labeldensity_red", "labeldensity_green", "labeldensity_blue", "labeldensity_rgb"};
	static const gchar *label_wcs[] = { "labelwcs_red", "labelwcs_green", "labelwcs_blue", "labelwcs_rgb" };

	if (com.cvport < RGB_VPORT) {
		gtk_label_set_text(GTK_LABEL(lookup_widget(label_density[com.cvport])), "");
		gtk_label_set_text(GTK_LABEL(lookup_widget(label_wcs[com.cvport])), "");

		if (inside) {
			static gchar buffer[256] = { 0 };
			static gchar wcs_buffer[256] = { 0 };
			static gchar format[256] = { 0 };
			int coords_width = 3;

			if (gfit.rx >= 1000 || gfit.ry >= 1000)
			coords_width = 4;
			if (gfit.type == DATA_USHORT && gfit.pdata[com.cvport] != NULL) {
				int val_width = 3;
				char *format_base_ushort = "x: %%.%dd y: %%.%dd (=%%.%dd)";
				if (gfit.hi >= 1000)
					val_width = 4;
				if (gfit.hi >= 10000)
					val_width = 5;
				g_sprintf(format, format_base_ushort,
						coords_width, coords_width, val_width);
				g_sprintf(buffer, format, zoomed.x, zoomed.y,
						gfit.pdata[com.cvport][gfit.rx * (gfit.ry - zoomed.y - 1)
											   + zoomed.x]);
			} else if (gfit.type == DATA_FLOAT  && gfit.fpdata[com.cvport] != NULL) {
				char *format_base_float = "x: %%.%dd y: %%.%dd (=%%f)";
				g_sprintf(format, format_base_float,
						coords_width, coords_width);
				g_sprintf(buffer, format, zoomed.x, zoomed.y,
						gfit.fpdata[com.cvport][gfit.rx * (gfit.ry - zoomed.y - 1) + zoomed.x]);
			}
			if (has_wcs(&gfit)) {
				double world_x, world_y;
				pix2wcs(&gfit, (double) zoomed.x, (double) (gfit.ry - zoomed.y - 1), &world_x, &world_y);
				if (world_x >= 0.0 && !isnan(world_x) && !isnan(world_y)) {
					SirilWorldCS *world_cs;

					world_cs = siril_world_cs_new_from_a_d(world_x, world_y);
					if (world_cs) {
						gchar *ra = siril_world_cs_alpha_format(world_cs, "%02dh%02dm%02ds");
						gchar *dec = siril_world_cs_delta_format(world_cs, "%c%02d°%02d\'%02d\"");
						g_sprintf(wcs_buffer, "α: %s δ: %s", ra, dec);

						gtk_label_set_text(GTK_LABEL(lookup_widget(label_wcs[com.cvport])), wcs_buffer);

						g_free(ra);
						g_free(dec);
						siril_world_cs_unref(world_cs);
					}
				}
			}

			if (buffer[0] != '\0') {
				gtk_label_set_text(GTK_LABEL(lookup_widget(label_density[com.cvport])), buffer);
			}
		}
	}

	if (com.translating) {
		update_zoom_fit_button();

		pointi ev = { (int)(event->x), (int)(event->y) };
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
					set_cursor("crosshair");
				}
			} else {
				if ((event->state & get_primary()) || is_inside_of_sel(zoomed, zoom)) {
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

gboolean update_zoom(gdouble x, gdouble y, double scale) {
	// event position in image coordinates before changing the zoom value
	point evpos = { x, y };
	cairo_matrix_transform_point(&com.image_matrix, &evpos.x, &evpos.y);
	gdouble factor;
	gboolean zoomed = FALSE;

	update_zoom_fit_button();

	com.zoom_value = get_zoom_val();

	factor = com.zoom_value * scale;

	if (factor >= ZOOM_MIN && factor <= ZOOM_MAX) {
		zoomed = TRUE;
		com.zoom_value = factor;
		adjust_vport_size_to_image();
		cairo_matrix_transform_point(&com.display_matrix, &evpos.x, &evpos.y);
		com.display_offset.x += x - evpos.x;
		com.display_offset.y += y - evpos.y;
		adjust_vport_size_to_image();
		redraw(com.cvport, REMAP_NONE);
	}
	return zoomed;
}

gboolean on_drawingarea_scroll_event(GtkWidget *widget, GdkEventScroll *event, gpointer user_data) {
	gboolean handled = FALSE;

	if (!single_image_is_loaded() && !sequence_is_loaded())
		return FALSE;

	if (event->state & get_primary()) {
		point delta;

		switch (event->direction) {
		case GDK_SCROLL_SMOOTH:	// what's that?
			gdk_event_get_scroll_deltas((GdkEvent*) event, &delta.x, &delta.y);
			if (delta.y < 0) {
				return update_zoom(event->x, event->y, ZOOM_IN);
			}
			if (delta.y > 0) {
				return update_zoom(event->x, event->y, ZOOM_OUT);
			}
			break;
		case GDK_SCROLL_DOWN:
			return update_zoom(event->x, event->y, ZOOM_OUT);
		case GDK_SCROLL_UP:
			return update_zoom(event->x, event->y, ZOOM_IN);
		default:
			handled = FALSE;
		}
	}
	return handled;
}
