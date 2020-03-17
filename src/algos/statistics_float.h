#ifndef _SIRIL_STATS_FLOAT_H
#define _SIRIL_STATS_FLOAT_H

#include "core/siril.h"

imstats* statistics_internal_float(fits *fit, int layer, rectangle *selection,
		int option, imstats *stats, gboolean multithread);
int compute_means_from_flat_cfa_float(fits *fit, float mean[4]);
#endif
