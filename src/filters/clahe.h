#ifndef SRC_GUI_CLAHE_H_
#define SRC_GUI_CLAHE_H_

/* Lucy-Richardson data from GUI */
struct CLAHE_data {
	fits *fit;
	double clip;
	int tileSize;
};

gpointer clahe(gpointer p);
void apply_clahe_cancel();

#endif /* SRC_GUI_CLAHE_H_ */
