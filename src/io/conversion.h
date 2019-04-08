#ifndef CONVERSION_H
#define CONVERSION_H

enum {
	COLUMN_FILENAME,		// string
	COLUMN_DATE,		// string
	N_COLUMNS_CONVERT
};

typedef struct {
	char *extension;			// name of the extension of raw
	char *manufacturer;			// name of the manufacturer
} supported_raw_list;

struct _convert_data {
	struct timeval t_start;
	GDir *dir;
	GList *list;
	int start;
	int total;
	int nb_converted;
	gboolean compatibility;
	gboolean stretch_cfa;
	gboolean command_line;
	gboolean several_type_of_files;
	gchar *destroot;
};

extern supported_raw_list supported_raw[];	//supported raw extensions
int get_nb_raw_supported();

void fill_convert_list(GSList *list);

void list_format_available();
image_type get_type_for_extension(const char *extension);
gchar *initialize_converters();
int count_selected_files();
int count_converted_files();
gpointer convert_thread_worker(gpointer p);
int debayer_if_needed(image_type imagetype, fits *fit, gboolean compatibility, gboolean force_debayer, gboolean stretch_cfa);
int any_to_fits(image_type imagetype, const char *source, fits *dest);
void set_debayer_in_convflags();

#endif
