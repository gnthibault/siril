#ifndef CONVERSION_H
#define CONVERSION_H

typedef struct {
	char *extension;			// name of the extension of raw
	char *manufacturer;			// name of the manufacturer
	sensor_pattern suggested_pattern;// type of bayer pattern. Not used for now
} supported_raw_list;

extern supported_raw_list supported_raw[];	//supported raw extensions
int get_nb_raw_supported();

void list_format_available();
void check_for_conversion_form_completeness();
image_type get_type_for_extension(const char *extension);
void initialize_converters();
int count_selected_files();
int save_to_target_fits(fits *fit, const char *dest_filename);
int debayer_if_needed(image_type imagetype, fits *fit, gboolean compatibility);
fits *any_to_new_fits(image_type imagetype, const char *source, gboolean compatibility);
int any_to_fits(image_type imagetype, const char *source, fits *dest);

#endif
