#ifndef SRC_FILTERS_BANDING_H_
#define SRC_FILTERS_BANDING_H_

#include "core/siril.h"

/* Banding data from GUI */
struct banding_data {
	fits *fit;
	double sigma;
	double amount;
	gboolean protect_highlights;
	gboolean applyRotation;
	const gchar *seqEntry;
};

void apply_banding_to_sequence(struct banding_data *banding_args);
gpointer BandingEngineThreaded(gpointer p);
int BandingEngine(fits *fit, double sigma, double amount, gboolean protect_highlights, gboolean applyRotation);

#endif /* SRC_FILTERS_BANDING_H_ */
