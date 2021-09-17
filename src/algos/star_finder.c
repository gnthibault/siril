/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2021 team free-astro (see more in AUTHORS file)
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
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <locale.h>
#include <gsl/gsl_matrix.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/Def_Wavelet.h"
#include "gui/utils.h"
#include "gui/message_dialog.h"
#include "gui/image_display.h"
#include "gui/progress_and_log.h"
#include "gui/PSF_list.h"
#include "algos/PSF.h"
#include "algos/star_finder.h"
#include "algos/statistics.h"
#include "filters/wavelets.h"
#include "io/image_format_fits.h"
#include "core/OS_utils.h"
#include "opencv/opencv.h"

#define WAVELET_SCALE 3

struct star_candidate_struct {
	int x, y;
	float mag_est;
};

static double guess_resolution(fits *fit) {
	double focal = fit->focal_length;
	double size = fit->pixel_size_x;
	double bin;

	if ((focal <= 0.0) || (size <= 0.0)) {
		focal = com.pref.focal;
		size = com.pref.pitch;
		/* if we have no way to guess, we return the
		 * flag -1
		 */
		if ((focal <= 0.0) || (size <= 0.0))
			return -1.0;
	}
	bin = ((fit->binning_x + fit->binning_y) / 2.0);
	if (bin <= 0) bin = 1.0;

	double res = RADCONV / focal * size * bin;

	/* test for high value. In this case we increase
	 * the number of detected star in is_star function
	 */
	if (res > 20.0) return -1.0;
	/* if res > 1.0 we use default radius value */
	if (res < 0.1 || res > 1.0) return 1.0;
	return res;
}

static float compute_threshold(fits *fit, double ksigma, int layer, rectangle *area, float *norm, double *bg, double *bgnoise) {
	float threshold;
	imstats *stat;

	assert(layer <= 3);

	stat = statistics(NULL, -1, fit, layer, area, STATS_BASIC, FALSE);
	if (!stat) {
		siril_log_message(_("Error: statistics computation failed.\n"));
		*norm = 0;
		*bg = 0.0;
		return 0;
	}
	threshold = (float)(stat->median + ksigma * stat->bgnoise);
	*norm = stat->normValue;
	*bg = stat->median;
	*bgnoise = stat->bgnoise;
	free_stats(stat);

	return threshold;
}

static gboolean is_star(psf_star *result, star_finder_params *sf ) {

	if (isnan(result->fwhmx) || isnan(result->fwhmy))
		return FALSE;
	if (isnan(result->x0) || isnan(result->y0))
		return FALSE;
	if (isnan(result->mag))
		return FALSE;
	if ((result->x0 <= 0.0) || (result->y0 <= 0.0))
		return FALSE;
	if (result->sx > 10 * sf->adj_radius || result->sy > 10 * sf->adj_radius)
		return FALSE;
	if (result->fwhmx <= 0.0 || result->fwhmy <= 0.0)
		return FALSE;
	if ((result->fwhmy / result->fwhmx) < sf->roundness)
		return FALSE;
	return TRUE;
}

int star_cmp(const void *a, const void *b)
{
    starc *a1 = (starc *)a;
    starc *a2 = (starc *)b;
    if ((*a1).mag_est > (*a2).mag_est)
        return -1;
    else if ((*a1).mag_est < (*a2).mag_est)
        return 1;
    else
        return 0;
}

static void get_structure(star_finder_params *sf) {
	static GtkSpinButton *spin_radius = NULL, *spin_sigma = NULL,
			*spin_roundness = NULL;
	static GtkToggleButton *toggle_adjust = NULL;

	if (spin_radius == NULL) {
		spin_radius = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_radius"));
		spin_sigma = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_threshold"));
		spin_roundness = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_round"));
		toggle_adjust = GTK_TOGGLE_BUTTON(lookup_widget("toggle_radius_adjust"));
	}
	sf->radius = (int) gtk_spin_button_get_value(spin_radius);
	sf->sigma = gtk_spin_button_get_value(spin_sigma);
	sf->adjust = gtk_toggle_button_get_active(toggle_adjust);
	sf->roundness = gtk_spin_button_get_value(spin_roundness);
}

void init_peaker_GUI() {
	/* TODO someday: read values from conf file and set them in the GUI.
	 * Until then, storing values in com.starfinder_conf instead of getting
	 * them in the GUI while running the peaker.
	 * see also init_peaker_default below */
	get_structure(&com.starfinder_conf);
}

void init_peaker_default() {
	/* values taken from siril3.glade */
	com.starfinder_conf.radius = 10;
	com.starfinder_conf.adjust = TRUE;
	com.starfinder_conf.sigma = 1.0;
	com.starfinder_conf.roundness = 0.5;
}

void on_toggle_radius_adjust_toggled(GtkToggleButton *togglebutton, gpointer user_data) {
	com.starfinder_conf.adjust = gtk_toggle_button_get_active(togglebutton);
}

void on_spin_sf_radius_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	com.starfinder_conf.radius = (int)gtk_spin_button_get_value(spinbutton);
}

void on_spin_sf_threshold_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	com.starfinder_conf.sigma = gtk_spin_button_get_value(spinbutton);
}

void on_spin_sf_roundness_changed(GtkSpinButton *spinbutton, gpointer user_data) {
	com.starfinder_conf.roundness = gtk_spin_button_get_value(spinbutton);
}

void update_peaker_GUI() {
	static GtkSpinButton *spin_radius = NULL, *spin_sigma = NULL,
			*spin_roundness = NULL;
	static GtkToggleButton *toggle_adjust = NULL;

	if (spin_radius == NULL) {
		spin_radius = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_radius"));
		spin_sigma = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_threshold"));
		spin_roundness = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_round"));
		toggle_adjust = GTK_TOGGLE_BUTTON(lookup_widget("toggle_radius_adjust"));
	}
	gtk_spin_button_set_value(spin_radius, (double) com.starfinder_conf.radius);
	gtk_toggle_button_set_active(toggle_adjust, com.starfinder_conf.adjust);
	gtk_spin_button_set_value(spin_sigma, com.starfinder_conf.sigma);
	gtk_spin_button_set_value(spin_roundness, com.starfinder_conf.roundness);
}

void confirm_peaker_GUI() {
	static GtkSpinButton *spin_radius = NULL, *spin_sigma = NULL,
			*spin_roundness = NULL;

	if (spin_radius == NULL) {
		spin_radius = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_radius"));
		spin_sigma = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_threshold"));
		spin_roundness = GTK_SPIN_BUTTON(lookup_widget("spinstarfinder_round"));
	}
	gtk_spin_button_update(spin_radius);
	gtk_spin_button_update(spin_sigma);
	gtk_spin_button_update(spin_roundness);
}


/*
 This is an implementation of a simple peak detector algorithm which
 identifies any pixel that is greater than any of its eight neighbors.

 Original algorithm come from:
 Copyleft (L) 1998 Kenneth J. Mighell (Kitt Peak National Observatory)
 */
static int minimize_candidates(fits *image, star_finder_params *sf, double bg, starc *candidates, int nb_candidates, int layer, psf_star ***retval, gboolean limit_nbstars);

psf_star **peaker(fits *fit, int layer, star_finder_params *sf, int *nb_stars, rectangle *area, gboolean showtime, gboolean limit_nbstars) {
	int nx = fit->rx;
	int ny = fit->ry;
	int areaX0 = 0;
	int areaY0 = 0;
	int areaX1 = nx;
	int areaY1 = ny;
	int nbstars = 0;
	double bg, bgnoise;
	float threshold, norm;
	float **smooth_image;
	fits smooth_fit = { 0 };
	starc *candidates;

	struct timeval t_start, t_end;

	assert(nx > 0 && ny > 0);

	siril_log_color_message(_("Findstar: processing...\n"), "green");
	gettimeofday(&t_start, NULL);

	/* running statistics on the input image is best as it caches them */
	threshold = compute_threshold(fit, sf->sigma * 5.0, layer, area, &norm, &bg, &bgnoise);
	if (norm == 0.0f)
		return NULL;

	siril_debug_print("Threshold: %f (background: %f, norm: %f)\n", threshold, bg, norm);

	/* Removing wavelets and applying a Gaussian filter to select candidates
	 */
	if (extract_fits(fit, &smooth_fit, layer, TRUE)) {
		siril_log_color_message(_("Failed to copy the image for processing\n"), "red");
		return NULL;
	}

	if (cvUnsharpFilter(&smooth_fit, 3, 0)) {
		siril_log_color_message(_("Could not apply Gaussian filter, aborting\n"), "red");
		clearfits(&smooth_fit);
		return NULL;
	}

	/* Build 2D representation of smoothed image upside-down */
	smooth_image = malloc(ny * sizeof(float *));
	if (!smooth_image) {
		PRINT_ALLOC_ERR;
		clearfits(&smooth_fit);
		return NULL;
	}
	for (int k = 0; k < ny; k++) {
		smooth_image[ny - k - 1] = smooth_fit.fdata + k * nx;
	}

	if (area && area->w != 0 && area->h != 0) {
		areaX0 = area->x;
		areaY0 = area->y;
		areaX1 = area->w + areaX0;
		areaY1 = area->h + areaY0;

		if (areaX1 > nx || areaY1 > ny) {
			siril_log_color_message(_("Selection is larger than image\n"), "red");
			return NULL;
		}
	}

	candidates = malloc(MAX_STARS * sizeof(starc));
	if (!candidates) {
		clearfits(&smooth_fit);
		free(smooth_image);
		PRINT_ALLOC_ERR;
		return NULL;
	}

	double res = guess_resolution(fit);

	if (res < 0) {
		res = 1.0;
	}

	sf->adj_radius = sf->adjust ? sf->radius / res : sf->radius;
	// siril_log_message("Adjusted radius: %d\n", sf->adj_radius);
	int r = sf->adj_radius;
	int boxsize = (2 * r + 1)*(2 * r + 1);
	double locthreshold = sf->sigma * 5.0 * bgnoise;

	/* Search for candidate stars in the filtered image */
	for (int y = r + areaY0; y < areaY1 - r; y++) {
		for (int x = r + areaX0; x < areaX1 - r; x++) {
			float pixel = smooth_image[y][x];
			if (pixel > threshold && pixel < norm) {
				gboolean bingo = TRUE;
				float neighbor;
				double mean = 0., meanhigh = 0.;
				int count = 0;
				// making sure the central pixel is a local max compared to all neighbors in the search box
				for (int yy = y - r; yy <= y + r; yy++) {
					for (int xx = x - r; xx <= x + r; xx++) {
						if (xx == x && yy == y)
							continue;
						neighbor = smooth_image[yy][xx];
						if (neighbor > pixel) {
							bingo = FALSE;
							break;
						} else if (neighbor == pixel) {
							if ((xx <= x && yy <= y) || (xx > x && yy < y)) {
								bingo = FALSE;
								break;
							}
						}
						if (!bingo) break;
						mean += neighbor;
						count ++;
					}
				}
				if (count < boxsize - 1) {
					bingo = FALSE;
					continue;
				}
				x += r; //no other local max can be found in the next r pixels band

				count = 0;
				for (int yy = y - 1; yy <= y + 1; yy++) {
					for (int xx = x - r - 1; xx <= x - r + 1; xx++) { // x corrected by the anticipated shift
						if (xx == x - r && yy == y) {
							continue;
						}
						neighbor = smooth_image[yy][xx]; 
						if (neighbor <= threshold) {
							bingo = FALSE;
							break;
						}
						if (!bingo) break;
						meanhigh += neighbor;
						count++;
					}
				}
				if (count < 8) {
					bingo = FALSE;
					continue;
				}
				mean = (mean - meanhigh) / (boxsize - 9); // (boxsize - 1) pix in mean and 8 pix in meanhigh
				/* trying to remove false positives in nebs
				mean of 9 central pixels must be above mean of the whole search box excluding them
				i.e, an approximation of the local background, by a significant amount
				*/
				meanhigh = (meanhigh + pixel) / 9;
				if (meanhigh - mean <= locthreshold) { 
					bingo = FALSE;
					continue;
				}

				if (bingo && nbstars < MAX_STARS) {
					candidates[nbstars].x = x - r; // x corrected by the anticipated shift
					candidates[nbstars].y = y;
					candidates[nbstars].mag_est = meanhigh;
					nbstars++;
					if (nbstars == MAX_STARS) break;
				}
			}
		}
		if (nbstars == MAX_STARS) break;
	}
	clearfits(&smooth_fit);
	siril_debug_print("Candidates for stars: %d\n", nbstars);

	/* Check if candidates are stars by minimizing a PSF on each */
	psf_star **results;
	nbstars = minimize_candidates(fit, sf, bg, candidates, nbstars, layer, &results, limit_nbstars);
	if (nbstars == 0)
		results = NULL;
	sort_stars(results, nbstars);
	free(smooth_image);
	free(candidates);

	gettimeofday(&t_end, NULL);
	if (showtime)
		show_time(t_start, t_end);
	if (nb_stars)
		*nb_stars = nbstars;
	return results;
}

/* returns number of stars found, result is in parameters */
static int minimize_candidates(fits *image, star_finder_params *sf, double bg, starc *candidates, int nb_candidates, int layer, psf_star ***retval, gboolean limit_nbstars) {
	int radius = sf->adj_radius;
	int nx = image->rx;
	int ny = image->ry;
	gsl_matrix *z = gsl_matrix_alloc(radius * 2, radius * 2);
	WORD **image_ushort = NULL;
	float **image_float = NULL;
	gint nbstars = 0;

	if (image->type == DATA_USHORT) {
		image_ushort = malloc(ny * sizeof(WORD *));
		for (int k = 0; k < ny; k++)
			image_ushort[ny - k - 1] = image->pdata[layer] + k * nx;
	}
	else if (image->type == DATA_FLOAT) {
		image_float = malloc(ny * sizeof(float *));
		for (int k = 0; k < ny; k++)
			image_float[ny - k - 1] = image->fpdata[layer] + k * nx;
	}
	else return 0;

	psf_star **results = new_fitted_stars(nb_candidates);
	if (!results) {
		PRINT_ALLOC_ERR;
		return 0;
	}

	//sorting candidates by starc.mean values as an estimator of mag
	qsort(candidates, nb_candidates, sizeof(starc), star_cmp);

	for (int candidate = 0; candidate < nb_candidates; candidate++) {
		int x = candidates[candidate].x, y = candidates[candidate].y;
		int ii, jj, i, j;
		/* FILL z */
		if (image->type == DATA_USHORT) {
			for (jj = 0, j = y - radius; j < y + radius; j++, jj++) {
				for (ii = 0, i = x - radius; i < x + radius;
						i++, ii++) {
					gsl_matrix_set(z, ii, jj, (double)image_ushort[j][i]);
				}
			}
		} else {
			for (jj = 0, j = y - radius; j < y + radius; j++, jj++) {
				for (ii = 0, i = x - radius; i < x + radius;
						i++, ii++) {
					gsl_matrix_set(z, ii, jj, (double)image_float[j][i]);
				}
			}
		}

		psf_star *cur_star = psf_global_minimisation(z, bg, FALSE, FALSE, FALSE);
		if (cur_star) {
			if (is_star(cur_star, sf)) {
				//fwhm_to_arcsec_if_needed(image, cur_star);	// should we do this here?
				cur_star->layer = layer;
				int result_index = g_atomic_int_add(&nbstars, 1);
				cur_star->xpos = (x - radius) + cur_star->x0 - 1.0;
				cur_star->ypos = (y - radius) + cur_star->y0 - 1.0;
				results[result_index] = cur_star;
				//fprintf(stdout, "%03d: %11f %11f %f\n",
				//		result_index, cur_star->xpos, cur_star->ypos, cur_star->mag);
			}
			else free_psf(cur_star);
			if ((nbstars >= MAX_STARS_FITTED) && limit_nbstars) break;
		}
	}

	results[nbstars] = NULL;
	if (retval)
		*retval = results;
	gsl_matrix_free(z);
	if (image_ushort) free(image_ushort);
	if (image_float) free(image_float);
	return nbstars;
}

/* Function to add star one by one, from the selection rectangle, the
 * minimization is run and the star is detected and added to the list of stars.
 *
 * IF A STAR IS FOUND and not already present in com.stars, the return value is
 * the new star and index is set to the index of the new star in com.stars.
 * IF NO NEW STAR WAS FOUND, either because it was already in the list, or a
 * star failed to be detected in the selection, or any other error, the return
 * value is NULL and index is set to -1.
 */
psf_star *add_star(fits *fit, int layer, int *index) {
	int i = 0;
	gboolean already_found = FALSE;

	*index = -1;
	psf_star *result = psf_get_minimisation(&gfit, layer, &com.selection, FALSE, TRUE, TRUE);
	if (!result)
		return NULL;
	/* We do not check if it's matching with the "is_star()" criteria.
	 * Indeed, in this case the user can add manually stars missed by star_finder */

	if (com.stars && !com.star_is_seqdata) {
		// check if the star was already detected/peaked
		while (com.stars[i]) {
			if (fabs(result->x0 + com.selection.x - com.stars[i]->xpos) < 0.9
					&& fabs(com.selection.y + com.selection.h - result->y0
									- com.stars[i]->ypos) < 0.9)
				already_found = TRUE;
			i++;
		}
	} else {
		if (com.star_is_seqdata) {
			/* com.stars was allocated with a size of 2, we need to free it before reallocating */
			clear_stars_list();
		}
		com.stars = new_fitted_stars(MAX_STARS);
		if (!com.stars) {
			PRINT_ALLOC_ERR;
			return NULL;
		}
		com.star_is_seqdata = FALSE;
	}

	if (already_found) {
		free_psf(result);
		result = NULL;
		char *msg = siril_log_message(_("This star has already been picked !\n"));
		siril_message_dialog( GTK_MESSAGE_INFO, _("Peaker"), msg);
	} else {
		if (i < MAX_STARS) {
			result->xpos = result->x0 + com.selection.x - 0.5;
			result->ypos = com.selection.y + com.selection.h - result->y0 - 0.5;
			com.stars[i] = result;
			com.stars[i + 1] = NULL;
			*index = i;
		} else {
			free_psf(result);
			result = NULL;
		}
	}
	return result;
}

int get_size_star_tab() {
	int i = 0;
	while (com.stars[i])
		i++;
	return i;
}

/* Remove a star from com.stars, at index index. The star is freed. */
int remove_star(int index) {
	if (index < 0 || !com.stars || !com.stars[index])
		return 1;

	int N = get_size_star_tab() + 1;

	free_psf(com.stars[index]);
	memmove(&com.stars[index], &com.stars[index + 1],
			(N - index - 1) * sizeof(*com.stars));
	redraw(com.cvport, REMAP_NONE);
	return 0;
}

int compare_stars(const void* star1, const void* star2) {
	psf_star *s1 = *(psf_star**) star1;
	psf_star *s2 = *(psf_star**) star2;

	if (s1->mag < s2->mag)
		return -1;
	else if (s1->mag > s2->mag)
		return 1;
	else
		return 0;
}

void sort_stars(psf_star **stars, int total) {
	if (*(&stars))
		qsort(*(&stars), total, sizeof(psf_star*), compare_stars);
}

/* allocates a new psf_star structure with a size of n + 1. First element is initialized to NULL */
psf_star **new_fitted_stars(size_t n) {
	psf_star **stars = malloc((n + 1) * sizeof(psf_star *));
	if (stars) stars[0] = NULL;

	return stars;
}

void free_fitted_stars(psf_star **stars) {
	int i = 0;
	while (stars && stars[i])
		free_psf(stars[i++]);
	free(stars);
}

void FWHM_average(psf_star **stars, int nb, float *FWHMx, float *FWHMy, char **units) {
	if (stars && stars[0]) {
		int i = 0;

		*FWHMx = 0.0f;
		*FWHMy = 0.0f;
		*units = stars[0]->units;
		while (i < nb) {
			*FWHMx += (float) stars[i]->fwhmx;
			*FWHMy += (float) stars[i]->fwhmy;
			i++;
		}
		*FWHMx = *FWHMx / (float)i;
		*FWHMy = *FWHMy / (float)i;
	}
}
