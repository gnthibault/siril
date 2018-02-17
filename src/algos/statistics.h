#ifndef _SIRIL_STATS_H
#define _SIRIL_STATS_H

#define STATS_BASIC	(1 << 1)	// median, mean, sigma, noise, min, max
#define STATS_AVGDEV	(1 << 2)	// average absolute deviation
#define STATS_MAD	(1 << 3)	// median absolute deviation
#define STATS_BWMV	(1 << 5)	// bidweight midvariance
#define STATS_MAIN	STATS_BASIC | STATS_AVGDEV | STATS_MAD | STATS_BWMV

// Iterative K-sigma Estimator of Location and Scale. Takes time, needed only for stacking
#define STATS_IKSS	(1 << 6)	
#define STATS_EXTRA	STATS_MAIN | STATS_IKSS

#define	STATS_ZERO_NONE 0
#define	STATS_ZERO_NULLCHECK (!STATS_ZERO_NONE)

#include "core/siril.h"

imstats* statistics(sequence *seq, int image_index, fits *fit, int layer,
		rectangle *selection, int option, int nullcheck);

void allocate_stats(imstats **stat);
imstats* free_stats(imstats *stat);
void clear_stats(sequence *seq, int layer);

void add_stats_to_fit(fits *fit, int layer, imstats *stat);
void add_stats_to_seq(sequence *seq, int image_index, int layer, imstats *stat);

void copy_seq_stats_to_fit(sequence *seq, int index, fits *fit);
void save_stats_from_fit(fits *fit, sequence *seq, int index);
void invalidate_stats_from_fit(fits *fit);

#endif

