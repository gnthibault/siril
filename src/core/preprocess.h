#ifndef _PREPRO_H_
#define _PREPRO_H_

#include "siril.h"
#include "filters/cosmetic_correction.h"

/* preprocessing data from GUI */
struct preprocessing_data {
	struct timeval t_start;
	gboolean use_bias, use_dark, use_flat;
	fits *bias, *dark, *flat;
	gboolean use_dark_optim, use_cosmetic_correction;
	gboolean is_sequence;
	sequence *seq;
	gboolean autolevel;
	double sigma[2];
	long icold, ihot;
	deviant_pixel *dev;
	gboolean is_cfa;
	gboolean debayer;
	gboolean compatibility;
	gboolean stretch_cfa;
	gboolean equalize_cfa;
	float normalisation;
	int retval;
	const char *ppprefix;	 // prefix for output files
};

int preprocess_single_image(struct preprocessing_data *args);
void start_sequence_preprocessing(struct preprocessing_data *prepro, gboolean from_script);

#endif
