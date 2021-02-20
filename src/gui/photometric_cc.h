#ifndef SRC_GUI_PHOTOMETRIC_CC_H_
#define SRC_GUI_PHOTOMETRIC_CC_H_

#include <stdio.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/PSF.h"

typedef struct struct_coeff {
	float value;
	int channel;
} coeff;

struct photometric_cc_data {
	fits *fit;
	fitted_PSF **stars;
	GInputStream *bv_stream;
	rectangle bg_area;
	gboolean bg_auto;
	int n_channel;
};

void initialize_photometric_cc_dialog();
int apply_photometric_cc();
int get_photometry_catalog();

#endif /* SRC_GUI_PHOTOMETRIC_CC_H_ */
