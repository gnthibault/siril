#pragma once

#ifdef __SSE2__
#include <x86intrin.h>
#endif

float median3x3(float* array);
float median5x5(float* array);
float median7x7(float* array);
float median9x9(float* array);

#ifdef __SSE2__
__m128 median3x3sse(__m128* array);
__m128 median5x5sse(__m128* array);
__m128 median7x7sse(__m128* array);
__m128 median9x9sse(__m128* array);
#endif
