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


#ifndef SRC_REGISTRATION_MATCHING_APPLY_MATCH_H_
#define SRC_REGISTRATION_MATCHING_APPLY_MATCH_H_

#include "misc.h"
#include "core/siril_world_cs.h"


int apply_match(SirilWorldCS *px_cat_center, double *crpix, TRANS trans, double *alpha, double *delta);

#endif /* SRC_REGISTRATION_MATCHING_APPLY_MATCH_H_ */
