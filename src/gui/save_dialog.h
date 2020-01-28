#ifndef SRC_GUI_SAVE_DIALOG_H_
#define SRC_GUI_SAVE_DIALOG_H_

/* Savedialog data from GUI */
struct savedial_data {
	GtkEntry *entry;
	gint bitspersamples;
	gint quality;
	const gchar *filename;
	int bitpix;
	gboolean update_hilo;
	int retval;
};

enum {
	PAGE_TIFF, PAGE_JPG, PAGE_FITS, PAGE_MISC
};

void on_header_save_button_clicked();

#endif /* SRC_GUI_SAVE_DIALOG_H_ */
