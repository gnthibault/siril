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

#ifndef SRC_CORE_INITFILE_H_
#define SRC_CORE_INITFILE_H_

enum token_index {
	WD = 0,		/* Working Directory */
	RAW = 1,	/* Raw settings */
	BAY = 2,	/* Bayer settings */
	PRE = 3,	/* Preprocessing settings */
	REG = 4,	/* Registration settings */
	STK = 5,	/* Stacking settings */
	MISC = 6,	/* Miscellaneous settings */
	NOTOK
};

int	writeinitfile();
int	checkinitfile();

#endif /* SRC_CORE_INITFILE_H_ */
