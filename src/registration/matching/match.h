#if !defined(MATCH_H)
#define MATCH_H

#include "core/siril.h"

/*
 * this is default value for the base of the output file names.
 * The 4 output files will have extensions
 *
 *   ".mtA"       stars from first list which DO match
 *   ".mtB"       stars from second list which DO match
 *   ".unA"       stars from first list which DON'T match
 *   ".unB"       stars from second list which DON'T match
 */
#define AT_MATCH_OUTFILE   "matched"

int star_match(fitted_PSF **s1, fitted_PSF **s2, int n, TRANS *trans);

#endif   /* MATCH_H */
