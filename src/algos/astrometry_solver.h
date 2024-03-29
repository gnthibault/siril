#ifndef SRC_ALGOS_ASTROMETRY_SOLVER_H_
#define SRC_ALGOS_ASTROMETRY_SOLVER_H_

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
	GAIAEDR3,
	PPMXL,
	BRIGHT_STARS,
	APASS
} online_catalog;

enum {
	WCS_FORMALISM_2,
	WCS_FORMALISM_1
};

struct astrometry_data {
	image_solved *solution;
	online_catalog onlineCatalog;
	SirilWorldCS *cat_center;
	GFile *catalog_name;
	gchar *catalogStars;
	gboolean for_photometry_cc;
	gboolean downsample;
	gboolean use_cache;
	gboolean autocrop;
	double scale; // scale (resolution)
	double cropfactor; // image cropping for wide fields
	rectangle solvearea;
	double xoffset, yoffset; //offset of solvearea center wrt. image center
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
int fill_plate_solver_structure(struct astrometry_data *args);
void wcs_cd_to_pc(double cd[2][2], double pc[2][2], double cdelt[2]);
void wcs_pc_to_cd(double pc[2][2], double cdelt[2], double cd[2][2]);
gpointer match_catalog(gpointer p);

gboolean confirm_delete_wcs_keywords(fits *fit);
void invalidate_WCS_keywords(fits *fit);
void flip_bottom_up_astrometry_data(fits *fit);
void flip_left_right_astrometry_data(fits *fit);
void rotate_astrometry_data(fits *fit, point center, double angle, gboolean cropped);
void crop_astrometry_data(fits *fit, point shift);

SirilWorldCS *get_image_solved_px_cat_center(image_solved *image);
SirilWorldCS *get_image_solved_image_center(image_solved *image);
void set_focal_and_pixel_pitch();

#endif /* SRC_ALGOS_ASTROMETRY_SOLVER_H_ */
