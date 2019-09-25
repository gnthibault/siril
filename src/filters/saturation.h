#ifndef SRC_GUI_SATURATION_H_
#define SRC_GUI_SATURATION_H_

/* color saturation data from GUI */
struct enhance_saturation_data {
	fits *input, *output;
	double coeff, h_min, h_max;
	gboolean preserve;
};

void apply_satu_cancel();
gpointer enhance_saturation(gpointer p);

#endif /* SRC_GUI_SATURATION_H_ */
