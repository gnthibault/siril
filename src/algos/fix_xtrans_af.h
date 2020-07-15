/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
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
#ifndef SRC_ALGOS_FIX_XTRANS_AF_H_
#define SRC_ALGOS_FIX_XTRANS_AF_H_

#include "core/siril.h"

typedef struct {
	char *model;       // Model of X-Trans camera.
	rectangle af;      // Bounding box for auto focus pixels.
	rectangle sample;  // Sample size for computing adjustment.
} supported_xtrans_list;

int fix_xtrans_ac(fits *fit);

#endif /* SRC_ALGOS_FIX_XTRANS_AF_H_ */
