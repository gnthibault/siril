#ifndef SRC_GUI_WAVELETS_H_
#define SRC_GUI_WAVELETS_H_

#include "core/siril.h"

/* wavelets filter data from GUI */
struct wavelets_filter_data {
	fits *fit;
	int Nbr_Plan;
	int Type;
};

void apply_wavelets_cancel();
int get_wavelet_layers(fits *fit, int Nbr_Plan, int Plan, int Type, int reqlayer);
gpointer extract_plans(gpointer p);

#endif /* SRC_GUI_WAVELETS_H_ */
