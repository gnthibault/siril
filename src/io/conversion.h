#ifndef CONVERSION_H
#define CONVERSION_H

#include <glib.h>
#include "core/siril.h" // for image_type

#define XTRANS_1 4
#define XTRANS_2 5

typedef struct {
	char *extension;			// name of the extension of raw
	char *manufacturer;			// name of the manufacturer
} supported_raw_list;

struct _convert_data {
	struct timeval t_start;
	gchar **list;
	int start;
	int total;
	int nb_converted_files;
	gboolean input_has_a_seq;
	gboolean make_link;
	gchar *destroot;
	int retval;

	gboolean debayer;
	sequence_type output_type;
	gboolean multiple_output;	// multiple SER output
};

#define MAX_EXTENSIONS 50	// actual size of supported_extensions

extern supported_raw_list supported_raw[];	//supported raw extensions
extern char *supported_extensions[MAX_EXTENSIONS];
extern const char *filter_pattern[];

int retrieveBayerPatternFromChar(char *bayer);
const gchar *get_bayer_pattern_from_preferences();
int get_nb_raw_supported();

void list_format_available();
image_type get_type_for_extension(const char *extension);
gchar *initialize_converters();
gpointer convert_thread_worker(gpointer p);
int debayer_if_needed(image_type imagetype, fits *fit, gboolean force_debayer);
int any_to_fits(image_type imagetype, const char *source, fits *dest, gboolean interactive, gboolean force_float, gboolean debayer);

#endif
