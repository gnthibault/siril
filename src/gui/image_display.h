#ifndef _IMAGE_DISPLAY_H_
#define _IMAGE_DISPLAY_H_

void initialize_image_display();

void queue_redraw(int doremap);
void redraw(int vport, int remap);
double get_zoom_val();	// for image_interactions

#endif

