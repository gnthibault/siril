/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2017 team free-astro (see more in AUTHORS file)
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
#include "gui/callbacks.h"
#include "registration/registration.h"

mouse_status_enum mouse_status;

gboolean redraw_preview(GtkWidget *widget, cairo_t *cr, gpointer data) {
	int current_preview, shiftx = 0, shifty = 0;
	static GtkToggleButton *check_display_ref = NULL;
	static GtkWidget *labelRegRef = NULL;
	guint area_width = gtk_widget_get_allocated_width (widget);
	guint area_height = gtk_widget_get_allocated_height (widget);
	gboolean display_ref_image;

	if (check_display_ref == NULL) {
		check_display_ref = GTK_TOGGLE_BUTTON(gtk_builder_get_object(builder, "checkbutton_displayref"));
		labelRegRef = lookup_widget("labelRegRef");
	}
	display_ref_image = com.refimage_regbuffer && com.refimage_surface
			&& gtk_toggle_button_get_active(check_display_ref)
							&& !gtk_widget_get_visible(labelRegRef);

	if (widget == com.preview_area[0]) current_preview = 0;
	else if (widget == com.preview_area[1]) current_preview = 1;
	else {
		fprintf(stderr, "Uninitialized com.preview_area or unknown drawing area!\n");
		return TRUE;
	}
	// update previewW and previewH in case the drawing area was resized
	com.seq.previewW[current_preview] = area_width;
	com.seq.previewH[current_preview] = area_height;

	/* display current image with shifts */
	if (!com.preview_surface[current_preview]) {
		gchar text[32];
		cairo_set_source_rgba(cr, 0.5, 0.5, 0.5, 0.5);
		cairo_rectangle(cr, 0, 0, area_width, area_height);
		cairo_fill(cr);
		cairo_select_font_face(cr, "Arial", CAIRO_FONT_SLANT_NORMAL,
				CAIRO_FONT_WEIGHT_BOLD);
		cairo_set_font_size(cr, 15);
		cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
		cairo_move_to(cr, area_width / 2.0 - 30, area_height / 2.0);
		g_snprintf(text, sizeof(text), _("Preview %d"), current_preview + 1);
		cairo_show_text(cr, text);
		return TRUE;
	}
	cairo_translate(cr, area_width / 2.0 - com.seq.previewX[current_preview],
			area_height/2.0-com.seq.previewY[current_preview]);
	if (com.seq.regparam && com.seq.regparam[com.cvport]) {
		shiftx = com.seq.regparam[com.cvport][com.seq.current].shiftx;
		shifty = com.seq.regparam[com.cvport][com.seq.current].shifty;
	}
	if (shiftx || shifty)
		cairo_translate(cr, shiftx, -shifty);
	cairo_set_source_surface(cr, com.preview_surface[current_preview], 0, 0);
	/*if (display_ref_image)
		cairo_paint_with_alpha(cr, 0.7);
	else*/ cairo_paint(cr);

	/* display reference image */
	if (display_ref_image) {
		/* already translated above but undo the shift translate */
		if (shiftx || shifty)
			cairo_translate(cr, -shiftx, shifty);
		cairo_set_source_surface(cr, com.refimage_surface, 0, 0);
		cairo_paint_with_alpha(cr, 0.5);
	}

	return FALSE;
}

void redraw_previews() {
	int i;
	for (i=0; i<PREVIEW_NB; i++)
		gtk_widget_queue_draw(com.preview_area[i]);
}

void set_preview_area(int preview_area, int centerX, int centerY) {
	/* equivalent to remap() for full image visualization, called in the
	 * mouse click callback and image change.
	 * load the preview area from reference image, from current image,
	 * sets the cairo_surface_t */
	guint area_width = gtk_widget_get_allocated_width (com.preview_area[preview_area]);
	guint area_height = gtk_widget_get_allocated_height (com.preview_area[preview_area]);

	com.seq.previewX[preview_area] = centerX;
	com.seq.previewY[preview_area] = centerY;
	com.seq.previewW[preview_area] = area_width;
	com.seq.previewH[preview_area] = area_height;

	if (cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, gfit.rx) !=
			com.surface_stride[com.cvport] ||
			gfit.ry != com.surface_height[com.cvport] ||
			!com.preview_surface[preview_area]) {
		if (com.preview_surface[preview_area])
			cairo_surface_destroy(com.preview_surface[preview_area]);
		com.preview_surface[preview_area] =
			cairo_image_surface_create_for_data(com.graybuf[com.cvport],
					CAIRO_FORMAT_RGB24,
					gfit.rx, gfit.ry,
					com.surface_stride[com.cvport]);
		if (cairo_surface_status(com.preview_surface[preview_area]) != CAIRO_STATUS_SUCCESS) {
			fprintf(stderr, "Error creating the Cairo image surface for preview %d\n", preview_area);
			cairo_surface_destroy(com.preview_surface[preview_area]);
			com.preview_surface[preview_area] = NULL;
		}
	}
	// queue for redraw
	gtk_widget_queue_draw(com.preview_area[preview_area]);
	// flush to ensure all writing to the image was done and redraw the surface
	//cairo_surface_flush (com.preview_surface[preview_area]);
	//cairo_surface_mark_dirty (com.preview_surface[preview_area]);
}

void on_toggle_preview_toggled(GtkToggleButton *toggle, gpointer user_data) {
	static GtkToggleButton *preview1 = NULL;
	if (preview1 == NULL)
		preview1 = GTK_TOGGLE_BUTTON(gtk_builder_get_object(builder, "togglebutton2"));
	if (gtk_toggle_button_get_active(toggle)) {
		if (toggle == preview1)
			mouse_status = MOUSE_ACTION_SELECT_PREVIEW1;
		else mouse_status = MOUSE_ACTION_SELECT_PREVIEW2;
	}
	else {
		/* deactivate preview */
		int preview_area;
		mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
		if (toggle == preview1)
			preview_area = 0;
		else preview_area = 1;
		cairo_surface_destroy(com.preview_surface[preview_area]);
		com.preview_surface[preview_area] = NULL;
		com.seq.previewX[preview_area] = -1;
		com.seq.previewY[preview_area] = -1;
		com.seq.previewH[preview_area] = 0;
		com.seq.previewW[preview_area] = 0;
		// queue for redraw
		gtk_widget_queue_draw(com.preview_area[preview_area]);
	}
}

void on_checkbutton_displayref_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	redraw_previews();
}

void init_mouse() {
	mouse_status = MOUSE_ACTION_SELECT_REG_AREA;
}

/* display registration data (shift{x|y} for now) in the manual adjustments */
void adjust_reginfo() {
	static GtkSpinButton *spin_shiftx = NULL, *spin_shifty = NULL;
	gboolean set_sensitive;

	if (spin_shiftx == NULL) {
		spin_shiftx = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spinbut_shiftx"));
		spin_shifty = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spinbut_shifty"));
	}
	
	g_signal_handlers_block_by_func(spin_shiftx, on_spinbut_shift_value_change, NULL);
	g_signal_handlers_block_by_func(spin_shifty, on_spinbut_shift_value_change, NULL);
	if (com.seq.regparam == NULL || com.seq.regparam[com.cvport] == NULL) {
		gtk_spin_button_set_value(spin_shiftx, 0.);
		gtk_spin_button_set_value(spin_shifty, 0.);
	} else {
		gtk_spin_button_set_value(spin_shiftx,
				com.seq.regparam[com.cvport][com.seq.current].shiftx);
		gtk_spin_button_set_value(spin_shifty,
				com.seq.regparam[com.cvport][com.seq.current].shifty);
	}
	g_signal_handlers_unblock_by_func(spin_shiftx, on_spinbut_shift_value_change, NULL);
	g_signal_handlers_unblock_by_func(spin_shifty, on_spinbut_shift_value_change, NULL);
	set_sensitive = (com.seq.current != com.seq.reference_image);
	gtk_widget_set_sensitive(GTK_WIDGET(spin_shiftx), set_sensitive);
	gtk_widget_set_sensitive(GTK_WIDGET(spin_shifty), set_sensitive);
}

void on_spinbut_shift_value_change(GtkSpinButton *spinbutton, gpointer user_data) {
	static GtkSpinButton *spin_shiftx = NULL;
	static GtkComboBox *cbbt_layers = NULL;
	int current_layer, new_value;
	if (spin_shiftx == NULL) {
		spin_shiftx = GTK_SPIN_BUTTON(gtk_builder_get_object(builder, "spinbut_shiftx"));
		cbbt_layers = GTK_COMBO_BOX(gtk_builder_get_object(builder, "comboboxreglayer"));
	}
	if (com.seq.regparam == NULL) {
		/* allocated when the number of layers is loaded from the sequence,
		 * = shouldn't happen */
		fprintf(stderr, "regparam not allocated, sequence not loaded, never displayed or malformed\n");
		return;
	}

	/* current layer: com.cvport or selected reg layer? */
	current_layer = gtk_combo_box_get_active(cbbt_layers);
	if (current_layer != com.cvport) {
		siril_log_message(_("Switching to the registration layer\n"));
		activate_tab(current_layer);
	}
	if (com.seq.regparam[current_layer] == NULL) {
		printf("Allocating registration data for this layer\n");
		com.seq.regparam[current_layer] = calloc(com.seq.number, sizeof(regdata));
		if (com.seq.regparam[current_layer] == NULL) {
			printf("FAILED: could not allocate registration data\n");
			return;
		}
	}

	new_value = gtk_spin_button_get_value_as_int(spinbutton);
	if (spinbutton == spin_shiftx)
		com.seq.regparam[current_layer][com.seq.current].shiftx = new_value;
	else com.seq.regparam[current_layer][com.seq.current].shifty = new_value;
	writeseqfile(&com.seq);
	fill_sequence_list(&com.seq, current_layer);	// update list with new regparam
	redraw_previews();
}

