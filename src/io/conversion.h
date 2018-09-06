#ifndef CONVERSION_H
#define CONVERSION_H

typedef struct {
	char *extension;			// name of the extension of raw
	char *manufacturer;			// name of the manufacturer
	sensor_pattern suggested_pattern;// type of bayer pattern. Not used for now
} supported_raw_list;

struct _convert_data {
	struct timeval t_start;
	GDir *dir;
	GList *list;
	int start;
	int total;
	int nb_converted;
	gboolean compatibility;
	gboolean command_line;
	gchar *destroot;
};

extern supported_raw_list supported_raw[];	//supported raw extensions
int get_nb_raw_supported();

void list_format_available();
int retrieveBayerPattern(char *bayer);
void check_for_conversion_form_completeness();
image_type get_type_for_extension(const char *extension);
void initialize_converters();
int count_selected_files();
gpointer convert_thread_worker(gpointer p);
int debayer_if_needed(image_type imagetype, fits *fit, gboolean compatibility, gboolean force_debayer);
int any_to_fits(image_type imagetype, const char *source, fits *dest);
void set_debayer_in_convflags();
void unset_debayer_in_convflags();

#endif
