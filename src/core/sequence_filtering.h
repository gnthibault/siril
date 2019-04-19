#ifndef _SEQ_FILTERING_H_
#define _SEQ_FILTERING_H_

#include "core/siril.h"

struct stacking_args;
struct stacking_configuration;

typedef int (*seq_image_filter)(sequence *seq, int img_index, double param);

struct filtering_tuple {
	seq_image_filter filter;
	double param;
};

int seq_filter_all(sequence *seq, int nb_img, double any);
int seq_filter_included(sequence *seq, int nb_img, double any);
int seq_filter_fwhm(sequence *seq, int nb_img, double max_fwhm);
int seq_filter_quality(sequence *seq, int nb_img, double max_quality);
int seq_filter_roundness(sequence *seq, int nb_img, double min_rnd);

int compute_nb_filtered_images(sequence *seq, seq_image_filter filtering_criterion, double filtering_parameter);

seq_image_filter create_filter_prefixed_nonexisting_output(const char *prefix);

seq_image_filter create_multiple_filter(seq_image_filter filter1, double fparam1, ...);
seq_image_filter create_multiple_filter_from_list(struct filtering_tuple *filters);

int convert_stack_data_to_filter(struct stacking_configuration *arg, struct stacking_args *stackargs);
int setup_filtered_data(struct stacking_args *args);

int stack_fill_list_of_unfiltered_images(struct stacking_args *args);
double compute_highest_accepted_fwhm(sequence *seq, int layer, double percent);
double compute_lowest_accepted_quality(sequence *seq, int layer, double percent);
double compute_lowest_accepted_roundness(sequence *seq, int layer, double percent);

char *describe_filter(sequence *seq, seq_image_filter filtering_criterion, double filtering_parameter);

#endif
