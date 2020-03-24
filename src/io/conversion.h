#ifndef CONVERSION_H
#define CONVERSION_H

#include <glib.h>
#include "core/siril.h" // for image_type

typedef struct {
	char *extension;			// name of the extension of raw
	char *manufacturer;			// name of the manufacturer
} supported_raw_list;

struct _convert_data {
	struct timeval t_start;
	GDir *dir;
	gchar **list;
	int start;
	int total;
	int nb_converted;
	gboolean compatibility;
	gboolean command_line;
	gboolean input_has_a_seq;
	gchar *destroot;
	int retval;
};

#define MAX_EXTENSIONS 50	// actual size of supported_extensions

extern supported_raw_list supported_raw[];	//supported raw extensions
extern char *supported_extensions[MAX_EXTENSIONS];
extern char *filter_pattern[];
extern unsigned int convflags;

int get_nb_raw_supported();

void list_format_available();
image_type get_type_for_extension(const char *extension);
gchar *initialize_converters();
gpointer convert_thread_worker(gpointer p);
int debayer_if_needed(image_type imagetype, fits *fit, gboolean compatibility, gboolean force_debayer);
int any_to_fits(image_type imagetype, const char *source, fits *dest, gboolean interactive, gboolean force_float);
void set_debayer_in_convflags();
void unset_debayer_in_convflags();

#endif
