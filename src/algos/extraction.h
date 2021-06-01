#ifndef _EXTRACTION_H
#define _EXTRACTION_H

#include <glib.h>
#include "core/siril.h"
#include "io/ser.h"
#include "io/fits_sequence.h"

struct split_cfa_data {
	fits *fit;
	sequence *seq;
	const gchar *seqEntry;	// not used for Ha-OIII split

	/* below: internal algorithm usage */
	struct ser_struct *new_ser_ha;
	fitseq *new_fitseq_ha;

	struct ser_struct *new_ser_oiii;
	fitseq *new_fitseq_oiii;

	GList *processed_images;
};

int extractHa_ushort(fits *in, fits *Ha, sensor_pattern pattern);
int extractHa_float(fits *in, fits *Ha, sensor_pattern pattern);
void apply_extractHa_to_sequence(struct split_cfa_data *split_cfa_args);

int extractGreen_ushort(fits *in, fits *green, sensor_pattern pattern);
int extractGreen_float(fits *in, fits *green, sensor_pattern pattern);
void apply_extractGreen_to_sequence(struct split_cfa_data *split_cfa_args);

int extractHaOIII_ushort(fits *in, fits *Ha, fits *OIII, sensor_pattern pattern);
int extractHaOIII_float(fits *in, fits *Ha, fits *OIII, sensor_pattern pattern);
void apply_extractHaOIII_to_sequence(struct split_cfa_data *split_cfa_args);

int split_cfa_ushort(fits *in, fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3);
int split_cfa_float(fits *in, fits *cfa0, fits *cfa1, fits *cfa2, fits *cfa3);
void apply_split_cfa_to_sequence(struct split_cfa_data *split_cfa_args);

#endif
