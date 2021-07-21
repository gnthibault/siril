#ifndef FINDER_H_
#define FINDER_H_

void init_peaker_GUI();
void init_peaker_default();
void update_peaker_GUI();
void confirm_peaker_GUI();
psf_star **peaker(fits *fit, int layer, star_finder_params *sf, int *nb_stars, rectangle *area, gboolean showtime, gboolean limit_nbstars);
psf_star *add_star(fits *fit, int layer, int *index);
int remove_star(int index);
void sort_stars(psf_star **stars, int total);
psf_star **new_fitted_stars(size_t n);
void free_fitted_stars(psf_star **stars);
int count_stars(psf_star **stars);
void FWHM_average(psf_star **stars, int nb, float *FWHMx, float *FWHMy, char **units);

typedef struct star_candidate_struct starc;

#endif
