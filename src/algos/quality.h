/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2015 team free-astro (see more in AUTHORS file)
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
#ifndef SRC_QUALITY_H_
#define SRC_QUALITY_H_

// How many bright pixels to we average to get the real maximum value
#define MAXP 6
#define QMARGIN 0.1
#define QSUBSAMPLE_INC 1
#define QSUBSAMPLE_MAX 5
#define QSUBSAMPLE_MIN 3
#define THRESHOLD 40
#define QF_APERTURE_RADIUS 0

#undef DEBUG

enum {
	QUALTYPE_NORMAL, QUALTYPE_NINOX
};

double QualityEstimate(fits *fit, int layer, int qtype);
int FindCentre(fits *fit, double *x_avg, double *y_avg);

#endif /* SRC_QUALITY_H_ */
