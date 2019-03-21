#ifndef FWHM_LIST_H_
#define FWHM_LIST_H_

#include "core/siril.h"
#include "algos/PSF.h"

void add_star_to_list(fitted_PSF *);
void fill_stars_list(fits *fit, fitted_PSF **);
void refresh_stars_list(fitted_PSF **);
void clear_stars_list();

#endif
