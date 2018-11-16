#ifndef SRC_ALGOS_PLATESOLVER_H_
#define SRC_ALGOS_PLATESOLVER_H_

#include "core/siril.h"
#include "registration/matching/degtorad.h"

/* multiply by this to convert degrees to radians */
#ifndef PI
#define PI 3.14159265359
#endif

#define BRIGHTEST_STARS 2500
#define AT_MATCH_CATALOG_NBRIGHT   60

#define RADtoASEC (3600.0 * 180.0 / M_PI)

#define CDSSESAME "http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame"
#define VIZIERSESAME "http://vizier.cfa.harvard.edu/viz-bin/nph-sesame"

typedef enum {
	TYCHO2, NOMAD
} online_catalog;

enum resolver {
	RESOLVER_NED,
	RESOLVER_SIMBAD,
	RESOLVER_VIZIER,
	RESOLVER_NUMBER
};

/* median filter data from GUI */
struct plate_solver_data {
	online_catalog onlineCatalog;
	gchar *catalogStars;
	double scale; // scale (resolution)
	fits *fit;
	gchar *message; // error message
	int ret; // return value
	double pixel_size; // pixel size in µm
	gboolean manual; // Manual platesolving
};

struct RA_struct {
	int hour;
	int min;
	double sec;
};
typedef struct RA_struct RA;

struct Dec_struct {
	int degree;
	int min;
	double sec;
};
typedef struct Dec_struct DEC;

struct object {
	gchar *name;
	double radius;
	int maxRecords;
	RA RA;
	DEC Dec;
	point imageCenter;
	gboolean south;
};

struct image_solved_struct {
	point px_size;
	point px_cat_center;
	point fov;
	double x, y;
	double ra, dec;
	double resolution, pixel_size, focal;
};
typedef struct image_solved_struct image_solved;


gboolean confirm_delete_wcs_keywords(fits *fit);
void invalidate_WCS_keywords(fits *fit);

#endif /* SRC_ALGOS_PLATESOLVER_H_ */
