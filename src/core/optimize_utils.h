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
#ifndef SRC_CORE_OPTIMIZE_UTILS_H_
#define SRC_CORE_OPTIMIZE_UTILS_H_

#include <glib.h>
#ifdef __SSE2__
#include <x86intrin.h>
#endif

// some convenience functions from RT
inline __attribute__((always_inline)) float intpf(float a, float b, float c)
{
    // calculate a * b + (1 - a) * c
    // following is valid:
    // intp(a, b+x, c+x) = intp(a, b, c) + x
    // intp(a, b*x, c*x) = intp(a, b, c) * x
    return a * (b - c) + c;
}
#ifdef __SSE2__
inline __attribute__((always_inline)) __m128 intpsse(__m128 a, __m128 b, __m128 c)
{
    // calculate a * b + (1 - a) * c
    // following is valid:
    // intp(a, b+x, c+x) = intp(a, b, c) + x
    // intp(a, b*x, c*x) = intp(a, b, c) * x
    return a * (b - c) + c;
}
#endif

inline __attribute__((always_inline)) gboolean inInterval(float val, float low, float high)
{
    // returns TRUE if val is in [low;high]
    float maxVal = max(val, low);
    return val == min(maxVal, high);
}


#endif /* SRC_CORE_OPTIMIZE_UTILS_H_ */
