#ifndef _FILTERS_H_
#define _FILTERS_H_

#include <gtk/gtk.h>

typedef struct {
	char *name;			// name of the layer (a filter name)
	double wavelength;		// the wavelength of the filter, in nanometres
} narrow_filter;

typedef struct {
	char *name;			// name of the filter (number and color)
	char *rgb;			// rgb equivalent (GdkRGBA parsable)
} broad_filter;

extern narrow_filter narrow_band_filters[];
extern int get_nb_narrow_filters();
extern broad_filter broad_band_filters[];

void wavelength_to_RGB(double wavelength, GdkRGBA *rgb);

#endif
