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

#include "siril_fit_linear.h"

/* code from gsl */
/* changed for siril */
int siril_fit_linear(const float *x, const float *y, const float m_x,
		const float m_dx2, const size_t n, float *c0, float *c1) {

	float m_y = y[0];
	for (size_t i = 1; i < n; i++) {
		m_y += (y[i] - m_y) * x[i];
	}

	float m_dxdy = 0.f;
	float dx = -m_x;
	for (size_t i = 0; i < n; i++, dx += 1.f) {
		const float dy = y[i] - m_y;

		m_dxdy += (dx * dy - m_dxdy) * x[i];
	}

	/* In terms of y = a + b x */

	const float b = m_dxdy * m_dx2;
	const float a = m_y - m_x * b;

	*c0 = a;
	*c1 = b;

	return 0;
}
