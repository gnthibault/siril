/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2016 team free-astro (see more in AUTHORS file)
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

#include "filters.h"
#include <math.h>

/* A common narrow-band filter list. */
narrow_filter narrow_band_filters[] = {
	/* narrow band filters */
	{ "H-alpha", 656.28 },
	{ "H-beta", 486.1 },
	{ "O III", 500.7 },
	{ "S II", 671.7 },
	{ "N II", 658.35 }
};

int get_nb_narrow_filters() {
	return sizeof(narrow_band_filters) / sizeof(narrow_filter);
}

/*
 * TODO: Add common broad band filters with their number
 * example list: http://www.myastroshop.com.au/guides/filters.asp
 */
broad_filter broad_band_filters[] = {
	/* broad band filters */
	{ "#1 (red)", "#ff0000" },
	{ "#2 (blue)", "#0000ff" },
};

/** Taken from Earl F. Glynn's web page:
 * http://www.efg2.com/Lab/ScienceAndEngineering/Spectra.htm
 * http://www.physics.sfasu.edu/astro/color/spectra.html
 *
 * Fills the rgb struct passed as arg with RGB values corresponding to the
 * wavelength.
 *
 * I couldn't find why gamma is 0.8 anywhere. */
static double Gamma = 0.80;

void wavelength_to_RGB(double wavelength, GdkRGBA *rgb){
	double factor;
	double Red,Green,Blue;

	if((wavelength >= 380) && (wavelength<440)){
		Red = -(wavelength - 440) / (440 - 380);
		Green = 0.0;
		Blue = 1.0;
	}else if((wavelength >= 440) && (wavelength<490)){
		Red = 0.0;
		Green = (wavelength - 440) / (490 - 440);
		Blue = 1.0;
	}else if((wavelength >= 490) && (wavelength<510)){
		Red = 0.0;
		Green = 1.0;
		Blue = -(wavelength - 510) / (510 - 490);
	}else if((wavelength >= 510) && (wavelength<580)){
		Red = (wavelength - 510) / (580 - 510);
		Green = 1.0;
		Blue = 0.0;
	}else if((wavelength >= 580) && (wavelength<645)){
		Red = 1.0;
		Green = -(wavelength - 645) / (645 - 580);
		Blue = 0.0;
	}else if((wavelength >= 645) && (wavelength<781)){
		Red = 1.0;
		Green = 0.0;
		Blue = 0.0;
	}else{
		Red = 0.0;
		Green = 0.0;
		Blue = 0.0;
	};

	// Let the intensity fall off near the vision limits
	if((wavelength >= 380) && (wavelength<420)){
		factor = 0.3 + 0.7*(wavelength - 380) / (420 - 380);
	}else if((wavelength >= 420) && (wavelength<701)){
		factor = 1.0;
	}else if((wavelength >= 701) && (wavelength<781)){
		factor = 0.3 + 0.7*(780 - wavelength) / (780 - 700);
	}else{
		factor = 0.0;
	};

	// Don't want 0^x = 1 for x <> 0
	rgb->red = Red==0.0 ? 0 : pow(Red * factor, Gamma);
	rgb->green = Green==0.0 ? 0 : pow(Green * factor, Gamma);
	rgb->blue = Blue==0.0 ? 0 : pow(Blue * factor, Gamma);
	rgb->alpha = 1.0;
}
