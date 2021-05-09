#ifndef _SIRIL_STATS_H
#define _SIRIL_STATS_H

struct stat_data {
	fits *fit;
	int option;
	gchar **list;
	sequence *seq;
	gchar *csv_name;
	const gchar *seqEntry;	// not used for stats

};

#define NULL_STATS -999999.0

#define STATS_MINMAX            (1 << 0)    // min, max
#define STATS_SIGMEAN           (1 << 1)    // noise, mean, sigma
#define STATS_BASIC             (1 << 2)    // median, mean, sigma, noise, min, max
#define STATS_AVGDEV            (1 << 3)    // average absolute deviation
#define STATS_MAD               (1 << 4)    // median absolute deviation
#define STATS_BWMV              (1 << 5)    // bidweight midvariance
#define STATS_IKSS              (1 << 6)    // IKSS, used for normalisation

#define STATS_MAIN              STATS_BASIC | STATS_AVGDEV | STATS_MAD | STATS_BWMV
#define STATS_EXTRA             STATS_MAIN | STATS_IKSS
#define STATS_NORM              STATS_BASIC | STATS_MAD | STATS_IKSS   // IKSSlite needs ngoodpix, median, mad and IKSS

#include "core/siril.h"

imstats* statistics(sequence *seq, int image_index, fits *fit, int layer,
		rectangle *selection, int option, gboolean multithread);

int compute_means_from_flat_cfa(fits *fit, float mean[4]);

void allocate_stats(imstats **stat);
imstats* free_stats(imstats *stat);
void clear_stats(sequence *seq, int layer);

void add_stats_to_fit(fits *fit, int layer, imstats *stat);
void add_stats_to_seq(sequence *seq, int image_index, int layer, imstats *stat);
void add_stats_to_seq_backup(sequence *seq, int image_index, int layer, imstats *stat);

void copy_seq_stats_to_fit(sequence *seq, int index, fits *fit);
void save_stats_from_fit(fits *fit, sequence *seq, int index);
void invalidate_stats_from_fit(fits *fit);
void full_stats_invalidation_from_fit(fits *fit);

void apply_stats_to_sequence(struct stat_data *stat_args);

float siril_stats_ushort_sd_64(const WORD data[], const int N);
float siril_stats_ushort_sd_32(const WORD data[], const int N);
float siril_stats_float_sd(const float data[], const int N, float *mean);

#endif

