#ifndef SRC_ALGOS_BACKGROUND_EXTRACTION_H_
#define SRC_ALGOS_BACKGROUND_EXTRACTION_H_

#include "core/siril.h"

/* ORDER OF POLYNOMES */
typedef enum {
	BACKGROUND_POLY_1,
	BACKGROUND_POLY_2,
	BACKGROUND_POLY_3,
	BACKGROUND_POLY_4,
} poly_order;

typedef struct sample {
	double median[3]; // median of each channel of the sample (if color)
	double mean; // mean of the 3 channel of the sample (if color)
	double min, max;
	size_t size;
	point position;
	gboolean valid;
} background_sample;

struct background_data {
	int nb_of_samples;
	double tolerance;
	int correction;
	poly_order degree;
	fits *fit;
	sequence *seq;
	const gchar *seqEntry;
};

int get_sample_radius();
void free_background_sample_list(GSList *list);
GSList* add_background_sample(GSList *list, fits *fit, point pt);
GSList* remove_background_sample(GSList *orig, fits *fit, point pt);
void generate_background_samples(int nb_of_samples, double tolerance);
gboolean remove_gradient_from_image(int correction, poly_order degree);
void apply_background_extraction_to_sequence(struct background_data *background_args);

#endif /* SRC_ALGOS_BACKGROUND_EXTRACTION_H_ */
