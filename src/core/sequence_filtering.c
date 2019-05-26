/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2019 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdarg.h>
#include <float.h> // DBL_MIN, DBL_MAX
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include "sequence_filtering.h"
#include "proto.h" // is_readable_file()
#include "registration/registration.h"
#include "io/sequence.h"
#include "stacking/stacking.h"
#include "gui/progress_and_log.h"
#include "algos/sorting.h"

#define MAX_FILTERS 5

/******************* IMAGE FILTERING CRITERIA *******************/

/* a criterion exists for each image filtering method, and is called in a
 * processing to verify if an image should be included or not.
 * These functions have the same signature, defined in stacking.h as
 * seq_filter, and return 1 if the tested image will be included and 0 if not.
 * Several filters can be applied at the same time, using the multiple filter
 * that executes a list of filter functions (seq_filter_multiple()).
 */

int seq_filter_all(sequence *seq, int nb_img, double any) {
	return 1;
}

int seq_filter_included(sequence *seq, int nb_img, double any) {
	return (seq->imgparam[nb_img].incl);
}

/* filter for deep-sky */
int seq_filter_fwhm(sequence *seq, int nb_img, double max_fwhm) {
	int layer;
	if (!seq->regparam) return 0;
	layer = get_registration_layer(seq);
	if (layer == -1) return 0;
	if (!seq->regparam[layer]) return 0;
	if (seq->regparam[layer][nb_img].fwhm > 0.0f)
		return seq->regparam[layer][nb_img].fwhm <= max_fwhm;
	else return 0;
}

int seq_filter_roundness(sequence *seq, int nb_img, double min_rnd) {
	int layer;
	if (!seq->regparam) return 0;
	layer = get_registration_layer(seq);
	if (layer == -1) return 0;
	if (!seq->regparam[layer]) return 0;
	if (seq->regparam[layer][nb_img].roundness > 0.0f)
		return seq->regparam[layer][nb_img].roundness >= min_rnd;
	else return 0;
}

/* filter for planetary */
int seq_filter_quality(sequence *seq, int nb_img, double max_quality) {
	int layer;
	if (!seq->regparam) return 0;
	layer = get_registration_layer(seq);
	if (layer == -1) return 0;
	if (!seq->regparam[layer]) return 0;
	if (seq->regparam[layer][nb_img].quality > 0.0)
		return seq->regparam[layer][nb_img].quality >= max_quality;
	else return 0;
}

/* browse the images to konw how many fit the criterion, from global data */
// ensure that com.seq is loaded before passing it as seq: sequence_is_loaded()
int compute_nb_filtered_images(sequence *seq, seq_image_filter filtering_criterion, double filtering_parameter) {
	int i, count = 0;
	for (i = 0; i < seq->number; i++) {
		if (filtering_criterion(seq, i, filtering_parameter))
			count++;
	}
	fprintf(stdout, "number of filtered-in images: %d\n", count);
	return count;
}

/************** The funny existing file sequence filtering function *****************/

static const char *_filter_prefix;

static int seq_filter_output_doesnt_already_exists(sequence *seq, int in_index, double any) {
	if (!_filter_prefix)
		fprintf(stderr, "USING FILTER PREFIX WITHOUT INITIALIZATION\n");
	if (seq->type != SEQ_REGULAR || !_filter_prefix)
		return 0;
	char *dest = fit_sequence_get_image_filename_prefixed(seq, _filter_prefix, in_index);
	int retval = is_readable_file(dest);
	//fprintf(stdout, "file %s exists: %d\n", dest, retval);
	free(dest);
	return !retval;
}

/* in C, there's no dynamic function creation, like lambda functions, so we
 * simulate that by keeping the captured variables in a static field. The
 * returned function is simply the regular function above, which is kept
 * private in order to avoid its uninitialized use. The limitation with this
 * static field is that it can be used once at a time.
 */
seq_image_filter create_filter_prefixed_nonexisting_output(const char *prefix) {
	_filter_prefix = prefix;
	fprintf(stdout, "created filter for prefixed output `%s'\n", prefix);
	return seq_filter_output_doesnt_already_exists;
}

/****************** MULTIPLE FILTERING ********************/

static struct filtering_tuple _filters[MAX_FILTERS];

static int seq_filter_multiple(sequence *seq, int img_index, double any) {
	int f = 0;
	while (f < MAX_FILTERS && _filters[f].filter) {
		if (!_filters[f].filter(seq, img_index, _filters[f].param))
			return 0;
		f++;
	}
	return 1;
}

/* configure the multiple filter */
seq_image_filter create_multiple_filter(seq_image_filter filter1, double fparam1, ...) {
	va_list args;
	int nb_filters = 1;
	_filters[0].filter = filter1;
	_filters[0].param = fparam1;
	va_start(args, fparam1);
	seq_image_filter f;
	do {
		_filters[nb_filters].filter = NULL;
		f = va_arg(args, seq_image_filter);
		if (f) {
			assert(f == seq_filter_multiple); // no multiple of multiple
			_filters[nb_filters].filter = f;
			_filters[nb_filters].param = va_arg(args, double);
			nb_filters++;
		}
	} while (f && nb_filters < MAX_FILTERS);
	fprintf(stdout, "created multiple filter (%d filters)\n", nb_filters);
	va_end(args);
	return seq_filter_multiple;
}

seq_image_filter create_multiple_filter_from_list(struct filtering_tuple *filters) {
	int nb_filters = 0;
	while (nb_filters < MAX_FILTERS && filters[nb_filters].filter) {
		assert(filters[nb_filters].filter != seq_filter_multiple); // no multiple of multiple
		_filters[nb_filters].filter = filters[nb_filters].filter;
		_filters[nb_filters].param = filters[nb_filters].param;
		nb_filters++;
	}
	if (nb_filters < MAX_FILTERS)
		_filters[nb_filters].filter = NULL;
	fprintf(stdout, "created multiple filter (%d filters)\n", nb_filters);
	return seq_filter_multiple;
}

/******************* filtering set-up **********************/

// creates the filtering criterion from a stacking configuration
// raises an error if the configuration has duplicates
// at least two images to be selected
int convert_stack_data_to_filter(struct stacking_configuration *arg, struct stacking_args *stackargs) {
	int nb_filters = 0;
	int layer = get_registration_layer(stackargs->seq);
	struct filtering_tuple filters[5] = { { NULL, 0.0 } };

	if ((arg->f_fwhm_p > 0.0f && arg->f_fwhm > 0.0f) ||
			(arg->f_round_p > 0.0f && arg->f_round > 0.0f) ||
			(arg->f_quality_p > 0.0f && arg->f_quality > 0.0f)) {
		siril_log_message(_("Sequence filter: values can only be either literal or percent\n"));
		return 1;
	}
	if (arg->filter_included) {
		filters[nb_filters].filter = seq_filter_included;
		filters[nb_filters].param = stackargs->seq->selnum;
		siril_log_message(_("Using selected images filter (%d/%d of the sequence)\n"),
				stackargs->seq->selnum, stackargs->seq->number);
		nb_filters++;
	}
	if (arg->f_fwhm_p > 0.0f || arg->f_fwhm > 0.0f) {
		filters[nb_filters].filter = seq_filter_fwhm;
		filters[nb_filters].param = arg->f_fwhm > 0.f ? arg->f_fwhm :
			compute_highest_accepted_fwhm(stackargs->seq, layer, arg->f_fwhm_p);
		siril_log_message(_("Using star FWHM images filter (below %f)\n"),
				filters[nb_filters].param);
	       	nb_filters++;
	}
	if (arg->f_round_p > 0.0f || arg->f_round > 0.0f) {
		filters[nb_filters].filter = seq_filter_roundness;
		filters[nb_filters].param = arg->f_round > 0.f ? arg->f_round :
			compute_lowest_accepted_roundness(stackargs->seq, layer, arg->f_round_p);
		siril_log_message(_("Using star roundness images filter (above %f)\n"),
				filters[nb_filters].param);
		nb_filters++;
	}
	if (arg->f_quality_p > 0.0f || arg->f_quality > 0.0f) {
		filters[nb_filters].filter = seq_filter_quality;
		filters[nb_filters].param = arg->f_quality > 0.f ? arg->f_quality :
			compute_lowest_accepted_quality(stackargs->seq, layer, arg->f_quality_p);
		siril_log_message(_("Using image quality filter (below %f)\n"),
				filters[nb_filters].param);
		nb_filters++;
	}
	filters[nb_filters].filter = NULL;

	if (nb_filters == 0) {
		stackargs->filtering_criterion = seq_filter_all;
		stackargs->filtering_parameter = 0.0;
	}
	else if (nb_filters == 1) {
		stackargs->filtering_criterion = filters[0].filter;
		stackargs->filtering_parameter = filters[0].param;
	}
	else {
		stackargs->filtering_criterion = create_multiple_filter_from_list(filters);
		stackargs->filtering_parameter = -1.0;
	}
	return 0;
}

// have to be set or init before calling:
// seq, filtering_criterion, filtering_parameter, image_indices
int setup_filtered_data(struct stacking_args *args) {
	args->nb_images_to_stack = compute_nb_filtered_images(args->seq,
			args->filtering_criterion, args->filtering_parameter);
	if (args->nb_images_to_stack < 2) {
		siril_log_message(_("Provided filtering options do not allow at least two images to be processed.\n"));
		return 1;
	}
	if (args->image_indices)
		free(args->image_indices);
	return stack_fill_list_of_unfiltered_images(args);
}

/* fill the image_indices mapping for the args->image_indices array.
 * args->image_indices will be allocated to nb_images_to_stack. */
int stack_fill_list_of_unfiltered_images(struct stacking_args *args) {
	int i, j;
	if (args->image_indices) {
		int *newptr = realloc(args->image_indices, args->nb_images_to_stack * sizeof(int));
		if (!newptr) {
			PRINT_ALLOC_ERR;
			free(args->image_indices);
			args->image_indices = NULL;
			return 1;
		}
		args->image_indices = newptr;
	} else {
		args->image_indices = malloc(args->nb_images_to_stack * sizeof(int));
	}

	for (i = 0, j = 0; i < args->seq->number; i++) {
		if (args->filtering_criterion(
					args->seq, i,
					args->filtering_parameter)) {
			args->image_indices[j] = i;
			j++;
		}
		else if (i == args->ref_image) {
			siril_log_color_message(_("The reference image is not in the selected set of images. "
					"To avoid issues, please change it or change the filtering parameters.\n"), "red");
			args->ref_image = -1;
		}
	}
	if (args->ref_image == -1) {
		args->ref_image = args->image_indices[0];
		siril_log_message(_("Using image %d as temporary reference image\n"), args->ref_image);
	}
	return j != args->nb_images_to_stack;
}

typedef double (*regdata_selector)(regdata *reg);

static double regdata_fwhm(regdata *reg) { return reg->fwhm; }
static double regdata_roundness(regdata *reg) { return reg->roundness; }
static double regdata_quality(regdata *reg) { return reg->quality; }

/* from a percentage, find the lowest or highest accepted registration property
 * value for image filtering in sequences. */
static double generic_compute_accepted_value(sequence *seq, int layer, double percent, gboolean lower_is_better, regdata_selector datasel) {
	int i, number_images_with_data;
	double threshold, *val;
	double extreme_value = lower_is_better ? DBL_MAX : DBL_MIN;
	if (layer < 0 || !seq->regparam || !seq->regparam[layer]) {
		return 0.0;
	}
	val = malloc(seq->number * sizeof(double));
	if (!val) return 0.0;

	// copy values
	for (i = 0; i < seq->number; i++) {
		double data = datasel(&seq->regparam[layer][i]);
		val[i] = data <= 0.0f ? extreme_value : data;
	}

	//sort values
	quicksort_d(val, seq->number);

	if (val[seq->number-1] != extreme_value) {
		number_images_with_data = seq->number;
	} else {
		for (i = 0; i < seq->number; i++)
			if (val[i] == extreme_value)
				break;
		number_images_with_data = i;
		siril_log_message(_("Warning: some images don't have information available for best "
				"images selection, using only available data (%d images on %d).\n"),
				number_images_with_data, seq->number);
	}

	// get highest or lowest accepted (threshold)
	double images_number = (double)(number_images_with_data - 1);
	if (lower_is_better) {
		threshold = val[(int)(percent * images_number / 100.0)];
		if (threshold == extreme_value)
			threshold = 0.0;
	} else {
		threshold = val[(int)((100.0 - percent) * images_number / 100.0)];
	}

	free(val);
	return threshold;
}

double compute_highest_accepted_fwhm(sequence *seq, int layer, double percent) {
	return generic_compute_accepted_value(seq, layer, percent, TRUE, regdata_fwhm);
}

double compute_lowest_accepted_quality(sequence *seq, int layer, double percent) {
	return generic_compute_accepted_value(seq, layer, percent, FALSE, regdata_quality);
}

double compute_lowest_accepted_roundness(sequence *seq, int layer, double percent) {
	return generic_compute_accepted_value(seq, layer, percent, FALSE, regdata_roundness);
}

char *describe_filter(sequence *seq, seq_image_filter filtering_criterion, double filtering_parameter) {
	GString *str = g_string_sized_new(100);
	int nb_images_to_stack = compute_nb_filtered_images(seq,
			filtering_criterion, filtering_parameter);

	if (filtering_criterion == seq_filter_all) {
		g_string_printf(str, _("Processing all images in the sequence (%d)\n"), seq->number);
	}
	else if (filtering_criterion == seq_filter_included) {
		g_string_printf(str, _("Processing only selected images in the sequence (%d)\n"), seq->selnum);
	}
	else if (filtering_criterion == seq_filter_fwhm) {
		g_string_printf(str, _("Processing images of the sequence "
					"with a FWHM lower or equal than %g (%d)\n"),
				filtering_parameter, nb_images_to_stack);
	}
	else if (filtering_criterion == seq_filter_roundness) {
		g_string_printf(str, _("Processing images of the sequence "
					"with a roundness higher or equal than %g (%d)\n"),
				filtering_parameter, nb_images_to_stack);
	}
	else if (filtering_criterion == seq_filter_quality) {
		g_string_printf(str, _("Processing images of the sequence "
					"with a quality higher or equal than %g (%d)\n"),
				filtering_parameter, nb_images_to_stack);
	}
	else if (filtering_criterion == seq_filter_output_doesnt_already_exists) {
		g_string_printf(str, _("Processing images whose output don't already exist (%d)"),
				nb_images_to_stack);
	}
	else if (filtering_criterion == seq_filter_multiple) {
		int f = 0;
		while (f < MAX_FILTERS && _filters[f].filter) {
			struct filtering_tuple *filter = _filters + f;
			char *descr = describe_filter(seq, filter->filter, filter->param);
			if (descr && descr[0] != '\0') {
				descr[strlen(descr)-1] = '\0';	// remove the new line
				if (f) descr[0] = tolower(descr[0]);
			}
			g_string_append(str, descr);
			g_string_append(str, ", ");
			g_free(descr);
			f++;
		}
		g_string_append_printf(str, _("for a total of images processed of %d)\n"),
				nb_images_to_stack);
	}
	//fprintf(stdout, "FILTERING DECRIPTION: %s", str->str);
	return g_string_free(str, FALSE);
}
