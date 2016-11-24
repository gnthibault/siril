#ifndef COSMETIC_CORRECTION_H_
#define COSMETIC_CORRECTION_H_

#include "core/siril.h"

typedef struct deviant_struct deviant_pixel;

/* Banding data from GUI */
struct cosmetic_data {
	fits *fit;
	double sigma[2];
	long icold;
	long ihot;
	double amount;
	gboolean is_cfa;
	const gchar *seqEntry;
};

typedef enum {
	COLD_PIXEL, HOT_PIXEL
} typeOfDeviant;

struct deviant_struct {
	point p;
	typeOfDeviant type;
};

long count_deviant_pixels(fits *fit, double sig[2], long *icold, long *ihot);
deviant_pixel *find_deviant_pixels(fits *fit, double sig[2], long *icold, long *ihot);
int autoDetect(fits *fit, int layer, double sig[2], long *icold, long *ihot,
		double amount, gboolean is_cfa);
void apply_cosmetic_to_sequence(struct cosmetic_data *cosme_args);
gpointer autoDetectThreaded(gpointer p);
int cosmeticCorrection(fits *fit, deviant_pixel *dev, int size, gboolean is_CFA);
int cosmeticCorrOneLine(fits *fit, deviant_pixel dev, gboolean is_cfa);
int cosmeticCorrOnePoint(fits *fit, deviant_pixel dev, gboolean is_cfa);

#endif /* COSMETIC_CORRECTION_H_ */
