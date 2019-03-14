#ifndef CONVERSION_H
#define CONVERSION_H

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
void check_for_conversion_form_completeness();
image_type get_type_for_extension(const char *extension);
void initialize_converters();
int count_selected_files();
gpointer convert_thread_worker(gpointer p);
int debayer_if_needed(image_type imagetype, fits *fit, gboolean compatibility, gboolean force_debayer, gboolean stretch_cfa);
int any_to_fits(image_type imagetype, const char *source, fits *dest);
void set_debayer_in_convflags();
void unset_debayer_in_convflags();

void on_treeview_convert_drag_data_received(GtkWidget *widget,
		GdkDragContext *context, gint x, gint y,
		GtkSelectionData *selection_data, guint info, guint time,
		gpointer user_data);

#endif
