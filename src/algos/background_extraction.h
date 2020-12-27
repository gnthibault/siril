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

struct background_data {
	int nb_of_samples;
	double tolerance;
	int correction;
	poly_order degree;
	fits *fit;
	sequence *seq;
	const gchar *seqEntry;
};

typedef struct sample background_sample;

int get_sample_radius();
void free_background_sample_list(GSList *list);
GSList* add_background_sample(GSList *list, fits *fit, point pt);
GSList* remove_background_sample(GSList *orig, fits *fit, point pt);
void generate_background_samples(int nb_of_samples, double tolerance);
gboolean remove_gradient_from_image(int correction, poly_order degree);
void apply_background_extraction_to_sequence(struct background_data *background_args);

gboolean background_sample_is_valid(background_sample *sample);
gdouble background_sample_get_size(background_sample *sample);
point background_sample_get_position(background_sample *sample);

#endif /* SRC_ALGOS_BACKGROUND_EXTRACTION_H_ */
