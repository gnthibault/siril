/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Ingo Weyrich (heckflosse67@gmx.de)
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <cmath>
#include <iostream>

#include "array2D.h"
#include "LUT.h"
#include "rt_math.h"
#include "rt_algo.h"
//#define BENCHMARK
#include "StopWatch.h"
#include "opthelper.h"
#include "core/sleef.h"
#include "filters/deconv.h"
#undef min
#undef max
#undef SQR

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize ("O3")
#endif

namespace {

inline float computeXYZ2LabY(float f, const LUTf& cachefy)
{
    if (f < 0.f) {
        constexpr static float kappa = 24389.f / 27.f;
        return 327.68f * (kappa * f / 65535.f);
    } else if (f > 65535.f) {
        return 327.68f * (116.f * xcbrtf(f / 65535.f) - 16.f);
    } else {
        return cachefy[f];
    }
}

void RGB2L(float *R, float *G, float *B, float *L, const float wp[3][3], int width, const LUTf& cachefy)
{

#ifdef __SSE2__
    const vfloat maxvalfv = F2V(65535.f);
    const vfloat rmv = F2V(wp[1][0]);
    const vfloat gmv = F2V(wp[1][1]);
    const vfloat bmv = F2V(wp[1][2]);
#endif
    int i = 0;
    
#ifdef __SSE2__
    for(; i < width - 3; i+=4) {
        const vfloat rv = LVFU(R[i]);
        const vfloat gv = LVFU(G[i]);
        const vfloat bv = LVFU(B[i]);
        const vfloat yv = rmv * rv + gmv * gv + bmv * bv;

        if (_mm_movemask_ps((vfloat)vorm(vmaskf_gt(yv, maxvalfv), vmaskf_lt(yv, ZEROV)))) {
            // take slower code path for all 4 pixels if one of the values is > 65535.f. Still faster than non SSE2 version
            for(int k = 0; k < 4; ++k) {
                float y = yv[k];
                L[i + k] = computeXYZ2LabY(y, cachefy);
            }
        } else {
            STVFU(L[i], cachefy[yv]);
        }
    }
#endif
    for(; i < width; ++i) {
        const float rv = R[i];
        const float gv = G[i];
        const float bv = B[i];
        float y = wp[1][0] * rv + wp[1][1] * gv + wp[1][2] * bv;

        L[i] = computeXYZ2LabY(y, cachefy);
    }
}

void RGB2Y(const float* R, const float* G, const float* B, float* Y1, float * Y2, int W) {
    int i = 0;
#ifdef __SSE2__
    const vfloat c1v = F2V(0.2627f);
    const vfloat c2v = F2V(0.6780f);
    const vfloat c3v = F2V(0.0593f);
    for (; i < W - 3; i += 4) {
        const vfloat Rv = vmaxf(LVFU(R[i]), ZEROV);
        const vfloat Gv = vmaxf(LVFU(G[i]), ZEROV);
        const vfloat Bv = vmaxf(LVFU(B[i]), ZEROV);
        vfloat yv = c1v * Rv + c2v * Gv + c3v * Bv;
        STVFU(Y1[i], yv);
        STVFU(Y2[i], yv);
    }
#endif
    for (; i < W; ++i) {
        const float r = std::max(R[i], 0.f);
        const float g = std::max(G[i], 0.f);
        const float b = std::max(B[i], 0.f);
        Y1[i] = Y2[i] = 0.2627f * r + 0.6780f * g + 0.0593f * b;
    }
}

void compute13x13kernel(float sigma, float kernel[13][13]) {

    const double temp = -2.f * rtengine::SQR(sigma);
    float sum = 0.f;
    for (int i = -6; i <= 6; ++i) {
        for (int j = -6; j <= 6; ++j) {
            if((rtengine::SQR(i) + rtengine::SQR(j)) <= rtengine::SQR(3.0 * 2.0)) {
                kernel[i + 6][j + 6] = std::exp((rtengine::SQR(i) + rtengine::SQR(j)) / temp);
                sum += kernel[i + 6][j + 6];
            } else {
                kernel[i + 6][j + 6] = 0.f;
            }
        }
    }

    for (int i = 0; i < 13; ++i) {
        for (int j = 0; j < 13; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

void compute9x9kernel(float sigma, float kernel[9][9]) {

    const double temp = -2.f * rtengine::SQR(sigma);
    float sum = 0.f;
    for (int i = -4; i <= 4; ++i) {
        for (int j = -4; j <= 4; ++j) {
            if((rtengine::SQR(i) + rtengine::SQR(j)) <= rtengine::SQR(3.0 * 1.5)) {
                kernel[i + 4][j + 4] = std::exp((rtengine::SQR(i) + rtengine::SQR(j)) / temp);
                sum += kernel[i + 4][j + 4];
            } else {
                kernel[i + 4][j + 4] = 0.f;
            }
        }
    }

    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 9; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

void compute7x7kernel(float sigma, float kernel[7][7]) {

    const double temp = -2.f * rtengine::SQR(sigma);
    float sum = 0.f;
    for (int i = -3; i <= 3; ++i) {
        for (int j = -3; j <= 3; ++j) {
            if((rtengine::SQR(i) + rtengine::SQR(j)) <= rtengine::SQR(3.0 * 1.15)) {
                kernel[i + 3][j + 3] = std::exp((rtengine::SQR(i) + rtengine::SQR(j)) / temp);
                sum += kernel[i + 3][j + 3];
            } else {
                kernel[i + 3][j + 3] = 0.f;
            }
        }
    }

    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

void compute5x5kernel(float sigma, float kernel[5][5]) {

    const double temp = -2.f * rtengine::SQR(sigma);
    float sum = 0.f;
    for (int i = -2; i <= 2; ++i) {
        for (int j = -2; j <= 2; ++j) {
            if((rtengine::SQR(i) + rtengine::SQR(j)) <= rtengine::SQR(3.0 * 0.84)) {
                kernel[i + 2][j + 2] = std::exp((rtengine::SQR(i) + rtengine::SQR(j)) / temp);
                sum += kernel[i + 2][j + 2];
            } else {
                kernel[i + 2][j + 2] = 0.f;
            }
        }
    }

    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

void compute3x3kernel(float sigma, float kernel[3][3]) {

    const double temp = -2.f * rtengine::SQR(sigma);
    float sum = 0.f;
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            if((rtengine::SQR(i) + rtengine::SQR(j)) <= rtengine::SQR(3.0 * 0.84)) {
                kernel[i + 1][j + 1] = std::exp((rtengine::SQR(i) + rtengine::SQR(j)) / temp);
                sum += kernel[i + 1][j + 1];
            } else {
                kernel[i + 1][j + 1] = 0.f;
            }
        }
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

void gauss3x3div (float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[3][3])
{

    const float c11 = kernel[0][0];
    const float c10 = kernel[0][1];
    const float c00 = kernel[1][1];

    for (int i = 1; i < tileSize - 1; i++) {
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 1; j < tileSize - 1; j++) {
            const float val = c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) + 
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) + 
                              c00 * src[i][j];
            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

void gauss5x5div (float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[5][5])
{

    const float c21 = kernel[0][1];
    const float c20 = kernel[0][2];
    const float c11 = kernel[1][1];
    const float c10 = kernel[1][2];
    const float c00 = kernel[2][2];

    for (int i = 2; i < tileSize - 2; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 2; j < tileSize - 2; ++j) {
            const float val = c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

void gauss13x13div(float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[13][13])
{
    const float c60 = kernel[0][6];
    const float c53 = kernel[1][3];
    const float c52 = kernel[1][4];
    const float c51 = kernel[1][5];
    const float c50 = kernel[1][6];
    const float c44 = kernel[2][2];
    const float c42 = kernel[2][4];
    const float c41 = kernel[2][5];
    const float c40 = kernel[2][6];
    const float c33 = kernel[3][3];
    const float c32 = kernel[3][4];
    const float c31 = kernel[3][5];
    const float c30 = kernel[3][6];
    const float c22 = kernel[4][4];
    const float c21 = kernel[4][5];
    const float c20 = kernel[4][6];
    const float c11 = kernel[5][5];
    const float c10 = kernel[5][6];
    const float c00 = kernel[6][6];

    for (int i = 6; i < tileSize - 6; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 6; j < tileSize - 6; ++j) {
            const float val = c60 * (src[i - 6][j] + src[i][j - 6] + src[i][j + 6] + src[i + 6][j]) +
                              c53 * ((src[i - 5][j - 3] + src[i - 5][j + 3]) + (src[i - 3][j - 5] + src[i - 3][j + 5]) + (src[i + 3][j - 5] + src[i + 3][j + 5]) + (src[i + 5][j - 3] + src[i + 5][j + 3])) +
                              c52 * ((src[i - 5][j - 2] + src[i - 5][j + 2]) + (src[i - 2][j - 5] + src[i - 2][j + 5]) + (src[i + 2][j - 5] + src[i + 2][j + 5]) + (src[i + 5][j - 2] + src[i + 5][j + 2])) +
                              c51 * ((src[i - 5][j - 1] + src[i - 5][j + 1]) + (src[i - 1][j - 5] + src[i - 1][j + 5]) + (src[i + 1][j - 5] + src[i + 1][j + 5]) + (src[i + 5][j - 1] + src[i + 5][j + 1])) +
                              c50 * ((src[i - 5][j] + src[i][j - 5] + src[i][j + 5] + src[i + 5][j]) + ((src[i - 4][j - 3] + src[i - 4][j + 3]) + (src[i - 3][j - 4] + src[i - 3][j + 4]) + (src[i + 3][j - 4] + src[i + 3][j + 4]) + (src[i + 4][j - 3] + src[i + 4][j + 3]))) +
                              c44 * (src[i - 4][j - 4] + src[i - 4][j + 4] + src[i + 4][j - 4] + src[i + 4][j + 4]) +
                              c42 * ((src[i - 4][j - 2] + src[i - 4][j + 2]) + (src[i - 2][j - 4] + src[i - 2][j + 4]) + (src[i + 2][j - 4] + src[i + 2][j + 4]) + (src[i + 4][j - 2] + src[i + 4][j + 2])) +
                              c41 * ((src[i - 4][j - 1] + src[i - 4][j + 1]) + (src[i - 1][j - 4] + src[i - 1][j + 4]) + (src[i + 1][j - 4] + src[i + 1][j + 4]) + (src[i + 4][j - 1] + src[i + 4][j + 1])) +
                              c40 * (src[i - 4][j] + src[i][j - 4] + src[i][j + 4] + src[i + 4][j]) +
                              c33 * (src[i - 3][j - 3] + src[i - 3][j + 3] + src[i + 3][j - 3] + src[i + 3][j + 3]) +
                              c32 * ((src[i - 3][j - 2] + src[i - 3][j + 2]) + (src[i - 2][j - 3] + src[i - 2][j + 3]) + (src[i + 2][j - 3] + src[i + 2][j + 3]) + (src[i + 3][j - 2] + src[i + 3][j + 2])) +
                              c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

void gauss9x9div(float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[9][9])
{

    const float c42 = kernel[0][2];
    const float c41 = kernel[0][3];
    const float c40 = kernel[0][4];
    const float c33 = kernel[1][1];
    const float c32 = kernel[1][2];
    const float c31 = kernel[1][3];
    const float c30 = kernel[1][4];
    const float c22 = kernel[2][2];
    const float c21 = kernel[2][3];
    const float c20 = kernel[2][4];
    const float c11 = kernel[3][3];
    const float c10 = kernel[3][4];
    const float c00 = kernel[4][4];

    for (int i = 4; i < tileSize - 4; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 4; j < tileSize - 4; ++j) {
            const float val = c42 * ((src[i - 4][j - 2] + src[i - 4][j + 2]) + (src[i - 2][j - 4] + src[i - 2][j + 4]) + (src[i + 2][j - 4] + src[i + 2][j + 4]) + (src[i + 4][j - 2] + src[i + 4][j + 2])) +
                              c41 * ((src[i - 4][j - 1] + src[i - 4][j + 1]) + (src[i - 1][j - 4] + src[i - 1][j + 4]) + (src[i + 1][j - 4] + src[i + 1][j + 4]) + (src[i + 4][j - 1] + src[i + 4][j + 1])) +
                              c40 * (src[i - 4][j] + src[i][j - 4] + src[i][j + 4] + src[i + 4][j]) +
                              c33 * (src[i - 3][j - 3] + src[i - 3][j + 3] + src[i + 3][j - 3] + src[i + 3][j + 3]) +
                              c32 * ((src[i - 3][j - 2] + src[i - 3][j + 2]) + (src[i - 2][j - 3] + src[i - 2][j + 3]) + (src[i + 2][j - 3] + src[i + 2][j + 3]) + (src[i + 3][j - 2] + src[i + 3][j + 2])) +
                              c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

void gauss7x7div(float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[7][7])
{

    const float c31 = kernel[0][2];
    const float c30 = kernel[0][3];
    const float c22 = kernel[1][1];
    const float c21 = kernel[1][2];
    const float c20 = kernel[1][3];
    const float c11 = kernel[2][2];
    const float c10 = kernel[2][3];
    const float c00 = kernel[3][3];

    for (int i = 3; i < tileSize - 3; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 3; j < tileSize - 3; ++j) {
            const float val = c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

void gauss3x3mult(float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[3][3])
{
    const float c11 = kernel[0][0];
    const float c10 = kernel[0][1];
    const float c00 = kernel[1][1];

    for (int i = 1; i < tileSize - 1; i++) {
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 1; j < tileSize - 1; j++) {
            const float val = c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) + 
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) + 
                              c00 * src[i][j];
            dst[i][j] *= val;
        }
    }

}

void gauss5x5mult (float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[5][5])
{

    const float c21 = kernel[0][1];
    const float c20 = kernel[0][2];
    const float c11 = kernel[1][1];
    const float c10 = kernel[1][2];
    const float c00 = kernel[2][2];

    for (int i = 2; i < tileSize - 2; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 2; j < tileSize - 2; ++j) {
            const float val = c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] *= val;
        }
    }
}

void gauss13x13mult(float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[13][13])
{

    const float c60 = kernel[0][6];
    const float c53 = kernel[1][3];
    const float c52 = kernel[1][4];
    const float c51 = kernel[1][5];
    const float c50 = kernel[1][6];
    const float c44 = kernel[2][2];
    const float c42 = kernel[2][4];
    const float c41 = kernel[2][5];
    const float c40 = kernel[2][6];
    const float c33 = kernel[3][3];
    const float c32 = kernel[3][4];
    const float c31 = kernel[3][5];
    const float c30 = kernel[3][6];
    const float c22 = kernel[4][4];
    const float c21 = kernel[4][5];
    const float c20 = kernel[4][6];
    const float c11 = kernel[5][5];
    const float c10 = kernel[5][6];
    const float c00 = kernel[6][6];

    for (int i = 6; i < tileSize - 6; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 6; j < tileSize - 6; ++j) {
            const float val = c60 * (src[i - 6][j] + src[i][j - 6] + src[i][j + 6] + src[i + 6][j]) +
                              c53 * ((src[i - 5][j - 3] + src[i - 5][j + 3]) + (src[i - 3][j - 5] + src[i - 3][j + 5]) + (src[i + 3][j - 5] + src[i + 3][j + 5]) + (src[i + 5][j - 3] + src[i + 5][j + 3])) +
                              c52 * ((src[i - 5][j - 2] + src[i - 5][j + 2]) + (src[i - 2][j - 5] + src[i - 2][j + 5]) + (src[i + 2][j - 5] + src[i + 2][j + 5]) + (src[i + 5][j - 2] + src[i + 5][j + 2])) +
                              c51 * ((src[i - 5][j - 1] + src[i - 5][j + 1]) + (src[i - 1][j - 5] + src[i - 1][j + 5]) + (src[i + 1][j - 5] + src[i + 1][j + 5]) + (src[i + 5][j - 1] + src[i + 5][j + 1])) +
                              c50 * ((src[i - 5][j] + src[i][j - 5] + src[i][j + 5] + src[i + 5][j]) + ((src[i - 4][j - 3] + src[i - 4][j + 3]) + (src[i - 3][j - 4] + src[i - 3][j + 4]) + (src[i + 3][j - 4] + src[i + 3][j + 4]) + (src[i + 4][j - 3] + src[i + 4][j + 3]))) +
                              c44 * (src[i - 4][j - 4] + src[i - 4][j + 4] + src[i + 4][j - 4] + src[i + 4][j + 4]) +
                              c42 * ((src[i - 4][j - 2] + src[i - 4][j + 2]) + (src[i - 2][j - 4] + src[i - 2][j + 4]) + (src[i + 2][j - 4] + src[i + 2][j + 4]) + (src[i + 4][j - 2] + src[i + 4][j + 2])) +
                              c41 * ((src[i - 4][j - 1] + src[i - 4][j + 1]) + (src[i - 1][j - 4] + src[i - 1][j + 4]) + (src[i + 1][j - 4] + src[i + 1][j + 4]) + (src[i + 4][j - 1] + src[i + 4][j + 1])) +
                              c40 * (src[i - 4][j] + src[i][j - 4] + src[i][j + 4] + src[i + 4][j]) +
                              c33 * (src[i - 3][j - 3] + src[i - 3][j + 3] + src[i + 3][j - 3] + src[i + 3][j + 3]) +
                              c32 * ((src[i - 3][j - 2] + src[i - 3][j + 2]) + (src[i - 2][j - 3] + src[i - 2][j + 3]) + (src[i + 2][j - 3] + src[i + 2][j + 3]) + (src[i + 3][j - 2] + src[i + 3][j + 2])) +
                              c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] *= val;
        }
    }
}

void gauss9x9mult(float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[9][9])
{

    const float c42 = kernel[0][2];
    const float c41 = kernel[0][3];
    const float c40 = kernel[0][4];
    const float c33 = kernel[1][1];
    const float c32 = kernel[1][2];
    const float c31 = kernel[1][3];
    const float c30 = kernel[1][4];
    const float c22 = kernel[2][2];
    const float c21 = kernel[2][3];
    const float c20 = kernel[2][4];
    const float c11 = kernel[3][3];
    const float c10 = kernel[3][4];
    const float c00 = kernel[4][4];

    for (int i = 4; i < tileSize - 4; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 4; j < tileSize - 4; ++j) {
            const float val = c42 * ((src[i - 4][j - 2] + src[i - 4][j + 2]) + (src[i - 2][j - 4] + src[i - 2][j + 4]) + (src[i + 2][j - 4] + src[i + 2][j + 4]) + (src[i + 4][j - 2] + src[i + 4][j + 2])) +
                              c41 * ((src[i - 4][j - 1] + src[i - 4][j + 1]) + (src[i - 1][j - 4] + src[i - 1][j + 4]) + (src[i + 1][j - 4] + src[i + 1][j + 4]) + (src[i + 4][j - 1] + src[i + 4][j + 1])) +
                              c40 * (src[i - 4][j] + src[i][j - 4] + src[i][j + 4] + src[i + 4][j]) +
                              c33 * (src[i - 3][j - 3] + src[i - 3][j + 3] + src[i + 3][j - 3] + src[i + 3][j + 3]) +
                              c32 * ((src[i - 3][j - 2] + src[i - 3][j + 2]) + (src[i - 2][j - 3] + src[i - 2][j + 3]) + (src[i + 2][j - 3] + src[i + 2][j + 3]) + (src[i + 3][j - 2] + src[i + 3][j + 2])) +
                              c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];
            dst[i][j] *= val;
        }
    }
}

void gauss7x7mult(float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[7][7])
{

    const float c31 = kernel[0][2];
    const float c30 = kernel[0][3];
    const float c22 = kernel[1][1];
    const float c21 = kernel[1][2];
    const float c20 = kernel[1][3];
    const float c11 = kernel[2][2];
    const float c10 = kernel[2][3];
    const float c00 = kernel[3][3];

    for (int i = 3; i < tileSize - 3; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 3; j < tileSize - 3; ++j) {
            const float val = c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] *= val;
        }
    }
}

void buildClipMaskOneChannel(const float * const *channel, int W, int H, float** clipMask, float white)
{

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int row = 0; row < H; ++row) {
        for (int col = 0; col < W; ++col) {
            clipMask[row][col] = 1.f;
        }
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int row = 2; row < H - 2; ++row) {
        for (int col = 2; col < W - 2; ++col) {
            if (channel[row][col] >= white) {
                clipMask[row - 2][col - 1] = clipMask[row - 2][col] = clipMask[row - 2][col + 1] = 0.f;
                clipMask[row - 1][col - 2] = clipMask[row - 1][col - 1] = clipMask[row - 1][col] = clipMask[row - 1][col + 1] = clipMask[row - 1][col + 2] = 0.f;
                clipMask[row][col - 2] = clipMask[row][col - 1] = clipMask[row][col] = clipMask[row][col + 1] = clipMask[row][col + 2] = 0.f;
                clipMask[row + 1][col - 2] = clipMask[row + 1][col - 1] = clipMask[row + 1][col] = clipMask[row + 1][col + 1] = clipMask[row + 1][col + 2] = 0.f;
                clipMask[row + 2][col - 1] = clipMask[row + 2][col] = clipMask[row + 2][col + 1] = 0.f;
            }
        }
    }
}

void buildClipMaskThreeChannels(const float * const *channel1, const float * const *channel2, const float * const *channel3, int W, int H, float** clipMask, float white)
{

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int row = 0; row < H; ++row) {
        for (int col = 0; col < W; ++col) {
            clipMask[row][col] = 1.f;
        }
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int row = 2; row < H - 2; ++row) {
        for (int col = 2; col < W - 2; ++col) {
            if (rtengine::min(channel1[row][col], channel2[row][col], channel3[row][col]) >= white) {
                clipMask[row - 2][col - 1] = clipMask[row - 2][col] = clipMask[row - 2][col + 1] = 0.f;
                clipMask[row - 1][col - 2] = clipMask[row - 1][col - 1] = clipMask[row - 1][col] = clipMask[row - 1][col + 1] = clipMask[row - 1][col + 2] = 0.f;
                clipMask[row][col - 2] = clipMask[row][col - 1] = clipMask[row][col] = clipMask[row][col + 1] = clipMask[row][col + 2] = 0.f;
                clipMask[row + 1][col - 2] = clipMask[row + 1][col - 1] = clipMask[row + 1][col] = clipMask[row + 1][col + 1] = clipMask[row + 1][col + 2] = 0.f;
                clipMask[row + 2][col - 1] = clipMask[row + 2][col] = clipMask[row + 2][col + 1] = 0.f;
            }
        }
    }
}

bool checkForStop(float** tmpIThr, float** iterCheck, int fullTileSize, int border)
{
    for (int ii = border; ii < fullTileSize - border; ++ii) {
#ifdef __SSE2__
        for (int jj = border; jj < fullTileSize - border; jj += 4) {
            if (UNLIKELY(_mm_movemask_ps((vfloat)vmaskf_lt(LVFU(tmpIThr[ii][jj]), LVFU(iterCheck[ii - border][jj - border]))))) {
                return true;
            }
        }
#else
        for (int jj = border; jj < fullTileSize - border; ++jj) {
            if (tmpIThr[ii][jj] < iterCheck[ii - border][jj - border]) {
                return true;
            }
        }
#endif
    }
    return false;
}

void CaptureDeconvSharpening (float** luminance, const float* const * oldLuminance, const float * const * blend, int W, int H, double sigma, double sigmaCornerOffset, int iterations, bool checkIterStop)
{
BENCHFUN
    const bool is9x9 = (sigma <= 1.50 && sigmaCornerOffset == 0.0);
    const bool is7x7 = (sigma <= 1.15 && sigmaCornerOffset == 0.0);
    const bool is5x5 = (sigma <= 0.84 && sigmaCornerOffset == 0.0);
    const bool is3x3 = (sigma < 0.6 && sigmaCornerOffset == 0.0);
    float kernel13[13][13];
    float kernel9[9][9];
    float kernel7[7][7];
    float kernel5[5][5];
    float kernel3[3][3];
    if (is3x3) {
        compute3x3kernel(sigma, kernel3);
    } else if (is5x5) {
        compute5x5kernel(sigma, kernel5);
    } else if (is7x7) {
        compute7x7kernel(sigma, kernel7);
    } else if (is9x9) {
        compute9x9kernel(sigma, kernel9);
    } else {
        compute13x13kernel(sigma, kernel13);
    }

    constexpr int tileSize = 32;
    const int border = (is3x3 || is5x5 || is7x7) ? iterations <= 30 ? 5 : 7 : 8;
    const int fullTileSize = tileSize + 2 * border;
    const float cornerRadius = std::min<float>(2.f, sigma + sigmaCornerOffset);
    const float cornerDistance = sqrt(rtengine::SQR(W * 0.5f) + rtengine::SQR(H * 0.5f));
    const float distanceFactor = (cornerRadius - sigma) / cornerDistance;

    constexpr float minBlend = 0.01f;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        array2D<float> tmpIThr(fullTileSize, fullTileSize);
        array2D<float> tmpThr(fullTileSize, fullTileSize);
        tmpThr.fill(1.f);
        array2D<float> lumThr(fullTileSize, fullTileSize);
        array2D<float> iterCheck(tileSize, tileSize);
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16) collapse(2)
#endif
        for (int i = border; i < H - border; i+= tileSize) {
            for(int j = border; j < W - border; j+= tileSize) {
                const bool endOfCol = (i + tileSize + border) >= H;
                const bool endOfRow = (j + tileSize + border) >= W;
                // fill tiles
                if (endOfRow || endOfCol) {
                    // special handling for small tiles at end of row or column
                    float maxVal = 0.f;
                    if (checkIterStop) {
                        for (int k = 0, ii = endOfCol ? H - fullTileSize + border : i; k < tileSize; ++k, ++ii) {
                            for (int l = 0, jj = endOfRow ? W - fullTileSize + border : j; l < tileSize; ++l, ++jj) {
                                iterCheck[k][l] = oldLuminance[ii][jj] * blend[ii][jj] * 0.5f;
                                maxVal = std::max(maxVal, blend[ii][jj]);
                            }
                        }
                    } else {
                        for (int k = 0, ii = endOfCol ? H - fullTileSize + border : i; k < tileSize; ++k, ++ii) {
                            for (int l = 0, jj = endOfRow ? W - fullTileSize + border : j; l < tileSize; ++l, ++jj) {
                                maxVal = std::max(maxVal, blend[ii][jj]);
                            }
                        }
                    }
                    if (maxVal < minBlend) {
                        // no pixel of the tile has a blend factor >= minBlend => skip the tile
                        continue;
                    }
                    for (int k = 0, ii = endOfCol ? H - fullTileSize : i - border; k < fullTileSize; ++k, ++ii) {
                        for (int l = 0, jj = endOfRow ? W - fullTileSize : j - border; l < fullTileSize; ++l, ++jj) {
                            tmpIThr[k][l] = oldLuminance[ii][jj];
                            lumThr[k][l] = oldLuminance[ii][jj];
                        }
                    }
                } else {
                    float maxVal = 0.f;
                    if (checkIterStop) {
                        for (int ii = 0; ii < tileSize; ++ii) {
                            for (int jj = 0; jj < tileSize; ++jj) {
                                iterCheck[ii][jj] = oldLuminance[i + ii][j + jj] * blend[i + ii][j + jj] * 0.5f;
                                maxVal = std::max(maxVal, blend[i + ii][j + jj]);
                            }
                        }
                    } else {
                        for (int ii = 0; ii < tileSize; ++ii) {
                            for (int jj = 0; jj < tileSize; ++jj) {
                                maxVal = std::max(maxVal, blend[i + ii][j + jj]);
                            }
                        }
                    }
                    if (maxVal < minBlend) {
                        // no pixel of the tile has a blend factor >= minBlend => skip the tile
                        continue;
                    }
                    for (int ii = i; ii < i + fullTileSize; ++ii) {
                        for (int jj = j; jj < j + fullTileSize; ++jj) {
                            tmpIThr[ii - i][jj - j] = oldLuminance[ii - border][jj - border];
                            lumThr[ii - i][jj - j] = oldLuminance[ii - border][jj - border];
                        }
                    }
                }
                if (is3x3) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 3x3 gaussian blur and divide luminance by result of gaussian blur
                        gauss3x3div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel3);
                        gauss3x3mult(tmpThr, tmpIThr, fullTileSize, kernel3);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else if (is5x5) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 5x5 gaussian blur and divide luminance by result of gaussian blur
                        gauss5x5div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel5);
                        gauss5x5mult(tmpThr, tmpIThr, fullTileSize, kernel5);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else if (is7x7) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 5x5 gaussian blur and divide luminance by result of gaussian blur
                        gauss7x7div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel7);
                        gauss7x7mult(tmpThr, tmpIThr, fullTileSize, kernel7);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else if (is9x9) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 5x5 gaussian blur and divide luminance by result of gaussian blur
                        gauss9x9div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel9);
                        gauss9x9mult(tmpThr, tmpIThr, fullTileSize, kernel9);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else {
                    if (sigmaCornerOffset != 0.0) {
                        const float distance = sqrt(rtengine::SQR(i + tileSize / 2 - H / 2) + rtengine::SQR(j + tileSize / 2 - W / 2));
                        const float sigmaTile = static_cast<float>(sigma) + distanceFactor * distance;
                        if (sigmaTile >= 0.4f) {
                            if (sigmaTile > 1.50) { // have to use 13x13 kernel
                                float lkernel13[13][13];
                                compute13x13kernel(static_cast<float>(sigma) + distanceFactor * distance, lkernel13);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 13x13 gaussian blur and divide luminance by result of gaussian blur
                                    gauss13x13div(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel13);
                                    gauss13x13mult(tmpThr, tmpIThr, fullTileSize, lkernel13);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            } else if (sigmaTile > 1.15) { // have to use 9x9 kernel
                                float lkernel9[9][9];
                                compute9x9kernel(static_cast<float>(sigma) + distanceFactor * distance, lkernel9);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 9x9 gaussian blur and divide luminance by result of gaussian blur
                                    gauss9x9div(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel9);
                                    gauss9x9mult(tmpThr, tmpIThr, fullTileSize, lkernel9);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            } else if (sigmaTile > 0.84) { // have to use 7x7 kernel
                                float lkernel7[7][7];
                                compute7x7kernel(static_cast<float>(sigma) + distanceFactor * distance, lkernel7);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 7x7 gaussian blur and divide luminance by result of gaussian blur
                                    gauss7x7div(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel7);
                                    gauss7x7mult(tmpThr, tmpIThr, fullTileSize, lkernel7);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            } else { // can use 5x5 kernel
                                float lkernel5[5][5];
                                compute5x5kernel(static_cast<float>(sigma) + distanceFactor * distance, lkernel5);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 7x7 gaussian blur and divide luminance by result of gaussian blur
                                    gauss5x5div(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel5);
                                    gauss5x5mult(tmpThr, tmpIThr, fullTileSize, lkernel5);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            }
                        }
                    } else {
                        for (int k = 0; k < iterations; ++k) {
                            // apply 13x13 gaussian blur and divide luminance by result of gaussian blur
                            gauss13x13div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel13);
                            gauss13x13mult(tmpThr, tmpIThr, fullTileSize, kernel13);
                            if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                break;
                            }
                        }
                    }
                }
                if (endOfRow || endOfCol) {
                    // special handling for small tiles at end of row or column
                    for (int k = border, ii = endOfCol ? H - fullTileSize : i - border; k < fullTileSize - border; ++k) {
                        for (int l = border, jj = endOfRow ? W - fullTileSize : j - border; l < fullTileSize - border; ++l) {
                            luminance[ii + k][jj + l] = rtengine::intp(blend[ii + k][jj + l], tmpIThr[k][l], luminance[ii + k][jj + l]);
                        }
                    }
                } else {
                    for (int ii = border; ii < fullTileSize - border; ++ii) {
                        for (int jj = border; jj < fullTileSize - border; ++jj) {
                            luminance[i + ii - border][j + jj - border] = rtengine::intp(blend[i + ii - border][j + jj - border], tmpIThr[ii][jj], luminance[i + ii - border][j + jj - border]);
                        }
                    }
                }
            }
        }
    }
}

}

extern "C" gpointer deconvolution(gpointer p);

gpointer deconvolution(gpointer p) {


BENCHFUN

    auto args = (struct deconv_data *) p;
    const int W = args->fit->rx;
    const int H = args->fit->ry;
    const int channels = args->fit->naxes[2];
    constexpr float xyz_rgb[3][3] = {          // XYZ from RGB
                                    { 0.412453, 0.357580, 0.180423 },
                                    { 0.212671, 0.715160, 0.072169 },
                                    { 0.019334, 0.119193, 0.950227 }
                                };

    float contrast = args->contrast_threshold / 100.f;

    array2D<float> clipMask(W, H);

    constexpr float clipLimit = 0.95f;

    LUTf cachefy(65536, LUT_CLIP_BELOW);
    {
        int i = 0;
        int epsmaxint = 65535.0 * 216.0 / 24389.0;

        for (; i <= epsmaxint; i++)
        {
            cachefy[i] = 327.68 * (24389.0 / 27.0 * i / 65535.0);
        }

        for(; i < 65536; i++)
        {
            cachefy[i] = 327.68 * (116.0 * std::cbrt((double)i / 65535.0) - 16.0);
        }
    }

    array2D<float> redVals(W,H);
    array2D<float> greenVals(W,H);
    array2D<float> blueVals(W,H);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int i = 0; i < H; ++i) {
        int fitn = i * W;
        for (int j = 0; j < W; ++j) {
        	if (args->fit->type == DATA_USHORT) {
                redVals[i][j] = args->fit->pdata[RLAYER][fitn];
                if (channels == 3) {
                    greenVals[i][j] = args->fit->pdata[GLAYER][fitn];
                    blueVals[i][j] = args->fit->pdata[BLAYER][fitn];
                } else { // for single layer we have to fill green and blue with red
                    greenVals[i][j] = args->fit->pdata[RLAYER][fitn];
                    blueVals[i][j] = args->fit->pdata[RLAYER][fitn];
                }
        	} else {
        		redVals[i][j] = args->fit->fpdata[RLAYER][fitn] > 0.f ? args->fit->fpdata[RLAYER][fitn] * USHRT_MAX_SINGLE : 0.f;
                if (channels == 3) {
                    greenVals[i][j] = args->fit->fpdata[GLAYER][fitn] > 0.f ? args->fit->fpdata[GLAYER][fitn] * USHRT_MAX_SINGLE : 0.f;
                    blueVals[i][j] = args->fit->fpdata[BLAYER][fitn] > 0.f ? args->fit->fpdata[BLAYER][fitn]* USHRT_MAX_SINGLE : 0.f;
                } else { // for single layer we have to fill green and blue with red
                    greenVals[i][j] = args->fit->fpdata[RLAYER][fitn] > 0.f ? args->fit->fpdata[RLAYER][fitn] * USHRT_MAX_SINGLE : 0.f;
                    blueVals[i][j] = args->fit->fpdata[RLAYER][fitn] > 0.f ? args->fit->fpdata[RLAYER][fitn]* USHRT_MAX_SINGLE : 0.f;
                }
        	}

            ++fitn;
        }
    }

    if (channels == 1) {
        buildClipMaskOneChannel(redVals, W, H, clipMask, args->clip * clipLimit);
    } else if (channels == 3) {
        buildClipMaskThreeChannels(redVals, greenVals, blueVals, W, H, clipMask, args->clip * clipLimit);
    }

    array2D<float> L(W, H);
    array2D<float> YOld(W, H);
    array2D<float> YNew(W, H);
    
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int i = 0; i < H; ++i) {
        RGB2L(redVals[i], greenVals[i], blueVals[i], L[i], xyz_rgb, W, cachefy);
        RGB2Y(redVals[i], greenVals[i], blueVals[i], YOld[i], YNew[i], W);
    }

    // calculate contrast based blend factors to reduce sharpening in regions with low contrast
    buildBlendMask(L, clipMask, W, H, contrast, args->auto_contrast_threshold, clipMask);
//    std::cout << "contrast : " << contrast << std::endl;
    if (args->auto_contrast_threshold) {
        args->contrast_threshold = contrast * 100.0;
    }

    CaptureDeconvSharpening(YNew, YOld, clipMask, W, H, args->sigma, args->corner_radius, args->iterations, args->auto_limit);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int i = 0; i < H; ++i) {
        int fitn = i * W;
        for (int j = 0; j < W; ++j) {
            const float factor = YNew[i][j] / std::max(YOld[i][j], 0.00001f);
            if (args->fit->type == DATA_USHORT) {
            	args->fit->pdata[RLAYER][fitn] = rtengine::CLIP<float>(redVals[i][j] * factor);
            	if (channels == 3) {
            		args->fit->pdata[GLAYER][fitn] = rtengine::CLIP<float>(greenVals[i][j] * factor);
            		args->fit->pdata[BLAYER][fitn] = rtengine::CLIP<float>(blueVals[i][j] * factor);
            	}
            } else {
            	args->fit->fpdata[RLAYER][fitn] = (redVals[i][j] * factor) / USHRT_MAX_SINGLE;
            	if (channels == 3) {
            		args->fit->fpdata[GLAYER][fitn] = (greenVals[i][j] * factor) / USHRT_MAX_SINGLE;
            		args->fit->fpdata[BLAYER][fitn] = (blueVals[i][j] * factor) / USHRT_MAX_SINGLE;
            	}
            }
            ++fitn;
        }
    }

    return GINT_TO_POINTER(0);
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif
