#ifndef SRC_ALGOS_PLATESOLVER_H_
#define SRC_ALGOS_PLATESOLVER_H_

#include "core/siril.h"
#include "core/siril_world_cs.h"
#include "registration/matching/degtorad.h"

#define BRIGHTEST_STARS 2500
#define AT_MATCH_CATALOG_NBRIGHT   60
#define CROP_ALLOWANCE 1.2

#define RADtoASEC (3600.0 * 180.0 / M_PI)

#define CDSSESAME "http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame"
#define VIZIERSESAME "http://vizier.cfa.harvard.edu/viz-bin/nph-sesame"

typedef struct image_solved_struct image_solved;

typedef enum {
	TYCHO2,
	NOMAD,
	GAIA,
	PPMXL,
	BRIGHT_STARS,
	APASS
} online_catalog;

struct plate_solver_data {
	online_catalog onlineCatalog;
	gchar *catalogStars;
	gboolean for_photometry_cc;
	gboolean downsample;
	gboolean autocrop;
	double scale; // scale (resolution)
	double cropfactor; // image cropping for wide fields
	fits *fit;
	fits *fit_backup;
	gchar *message; // error message
	int ret; // return value
	double pixel_size; // pixel size in µm
	gboolean manual; // Manual platesolving
	gboolean flip_image;
};

void open_astrometry_dialog();
gchar *search_in_catalogs(const gchar *object);
int fill_plate_solver_structure(struct plate_solver_data *args);
gpointer match_catalog(gpointer p);

gboolean confirm_delete_wcs_keywords(fits *fit);
void invalidate_WCS_keywords(fits *fit);

SirilWorldCS *get_image_solved_px_cat_center(image_solved *image);
SirilWorldCS *get_image_solved_image_center(image_solved *image);
double get_image_solved_x(image_solved *image);
double get_image_solved_y(image_solved *image);
void set_focal_and_pixel_pitch();
void update_image_center_coord(image_solved *image, gdouble alpha, gdouble delta);

#endif /* SRC_ALGOS_PLATESOLVER_H_ */
