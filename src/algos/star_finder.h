#ifndef FINDER_H_
#define FINDER_H_

void init_peaker_GUI();
void init_peaker_default();
void update_peaker_GUI();
fitted_PSF **peaker(fits *fit, int layer, star_finder_params *sf, int *nb_stars, rectangle *area, gboolean showtime);
fitted_PSF *add_star(fits *fit, int layer, int *index);
int remove_star(int index);
void sort_stars(fitted_PSF **stars, int total);
void free_fitted_stars(fitted_PSF **stars);
int count_stars(fitted_PSF **stars);
void FWHM_average(fitted_PSF **stars, int nb, float *FWHMx, float *FWHMy, char **units);

#endif
