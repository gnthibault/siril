#ifndef FWHM_LIST_H_
#define FWHM_LIST_H_

#include "core/siril.h"
#include "algos/PSF.h"

void add_star_to_list(psf_star *);
void fill_stars_list(fits *fit, psf_star **);
void refresh_star_list(psf_star **);
void clear_stars_list();
void pick_a_star();

void popup_psf_result(psf_star *result, rectangle *area);

#endif
