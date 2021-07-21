#if !defined(MATCH_H)
#define MATCH_H

#include "core/siril.h"


#define NB_OF_MATCHING_TRY 3


int new_star_match(psf_star **s1, psf_star **s2, int n, int nobj_override,
		double s_min, double s_max,
		Homography *H, gboolean print_output, transformation_type type);

#endif   /* MATCH_H */
