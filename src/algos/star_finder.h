#ifndef FINDER_H_
#define FINDER_H_

typedef struct starFinder_struct starFinder;

struct starFinder_struct {
	int radius;
	double sigma;
	double roundness;
	int nb_stars;
};

fitted_PSF **peaker(fits *fit, int layer, starFinder *sf, rectangle *area);
fitted_PSF *add_star(fits *fit, int layer, int *index);
int remove_star(int index);
void sort_stars(fitted_PSF **stars, int total);
int count_stars(fitted_PSF **stars);
void FWHM_average(fitted_PSF **stars, float *FWHMx, float *FWHMy, int max);

#endif
