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
int siril_fit_linear(const float *x, const float *y, const size_t n, float *c0,
		float *c1) {
	float m_x = 0.f, m_y = 0.f, m_dx2 = 0.f, m_dxdy = 0.f;

	size_t i;

	for (i = 0; i < n; i++) {
		m_x += (x[i] - m_x) / (i + 1.f);
		m_y += (y[i] - m_y) / (i + 1.f);
	}

	for (i = 0; i < n; i++) {
		const float dx = x[i] - m_x;
		const float dy = y[i] - m_y;

		m_dx2 += (dx * dx - m_dx2) / (i + 1.f);
		m_dxdy += (dx * dy - m_dxdy) / (i + 1.f);
	}

	/* In terms of y = a + b x */

	{
		float s2 = 0.f, d2 = 0.f;
		float b = m_dxdy / m_dx2;
		float a = m_y - m_x * b;

		*c0 = a;
		*c1 = b;

		/* Compute chi^2 = \sum (y_i - (a + b * x_i))^2 */

		for (i = 0; i < n; i++) {
			const float dx = x[i] - m_x;
			const float dy = y[i] - m_y;
			const float d = dy - b * dx;
			d2 += d * d;
		}

		s2 = d2 / (n - 2.f); /* chisq per degree of freedom */

	}

	return 0;
}
