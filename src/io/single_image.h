#ifndef SINGLE_H_
#define SINGLE_H_

void close_single_image();
void free_image_data();
int read_single_image(const char* filename, fits *dest, char **realname_out);
int open_single_image(const char* filename);
void open_single_image_from_gfit(char *realname);

int image_find_minmax(fits *fit, int force_minmax);			// was minmax()
void init_layers_hi_and_lo_values(sliders_mode force_minmax);
void adjust_cutoff_from_updated_gfit();		// was level_adjust()
void unique_free_preprocessing_data(single *uniq);
int single_image_is_loaded();

#endif
