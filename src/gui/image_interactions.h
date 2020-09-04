#ifndef _IMAGE_INTERACTIONS_H_
#define _IMAGE_INTERACTIONS_H_

typedef void (*selection_update_callback)();
void register_selection_update_callback(selection_update_callback f);
void unregister_selection_update_callback(selection_update_callback f);
void delete_selected_area();
void reset_display_offset();
void reset_zoom_default();

/* mouse behaviour */
typedef enum {
	MOUSE_ACTION_NONE,
	MOUSE_ACTION_SELECT_REG_AREA,
	MOUSE_ACTION_SELECT_PREVIEW1,
	MOUSE_ACTION_SELECT_PREVIEW2,
	MOUSE_ACTION_DRAW_SAMPLES,
} mouse_status_enum;

extern mouse_status_enum mouse_status;	// defined in registration_preview.c

#endif
