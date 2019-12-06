#ifndef SRC_GUI_OPEN_DIALOG_H_
#define SRC_GUI_OPEN_DIALOG_H_

/* cookies for the file chooser */
#define OD_NULL 	0
#define OD_FLAT 	1
#define OD_DARK 	2
#define OD_OFFSET 	3
#define OD_CWD 		4
#define OD_OPEN 	5
#define OD_CONVERT 	6

void header_open_button_clicked();
void cwd_btton_clicked();

typedef struct _fileChooserPreview fileChooserPreview;

struct _fileChooserPreview {
	GtkWidget *image;
	GtkWidget *dim_label;
	GtkWidget *size_label;
};

#endif /* SRC_GUI_OPEN_DIALOG_H_ */
