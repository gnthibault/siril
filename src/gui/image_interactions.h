#ifndef _IMAGE_INTERACTIONS_H_
#define _IMAGE_INTERACTIONS_H_

#include "core/siril.h"

typedef void (*selection_update_callback)();
typedef void (*star_selection_callback)(pointi);
gboolean update_zoom(gdouble x, gdouble y, double scale);
void update_zoom_fit_button();

void register_selection_update_callback(selection_update_callback f);
void unregister_selection_update_callback(selection_update_callback f);

void new_selection_zone();
void delete_selected_area();
void reset_display_offset();
void reset_zoom_default();

void enforce_ratio_and_clamp();

gboolean display_quick_photo();

/* mouse behaviour */
typedef enum {
	MOUSE_ACTION_NONE,
	MOUSE_ACTION_SELECT_REG_AREA,
	MOUSE_ACTION_SELECT_PREVIEW1,
	MOUSE_ACTION_SELECT_PREVIEW2,
	MOUSE_ACTION_DRAW_SAMPLES,
	MOUSE_ACTION_PHOTOMETRY
} mouse_status_enum;

extern mouse_status_enum mouse_status;	// defined in registration_preview.c

#endif
