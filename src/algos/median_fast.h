#pragma once

#ifdef __SSE2__
#include <x86intrin.h>
#endif

float median3x3(float* array);
float median24(float *array);
float median5x5(float* array);
float median7x7(float* array);
float median9x9(float* array);

inline __attribute__((always_inline)) float mymin(float a, float b) {
    return b < a ? b : a;
}

inline __attribute__((always_inline))float mymax(float a, float b) {
    return a < b ? b : a;
}

inline __attribute__((always_inline)) float median9f(float array0, float array1, float array2, float array3, float array4, float array5, float array6, float array7, float array8)
{
    float tmp = mymin(array1, array2);
    array2 = mymax(array1, array2);
    array1 = tmp;
    tmp = mymin(array4, array5);
    array5 = mymax(array4, array5);
    array4 = tmp;
    tmp = mymin(array7, array8);
    array8 = mymax(array7, array8);
    array7 = tmp;
    tmp = mymin(array0, array1);
    array1 = mymax(array0, array1);
    array0 = tmp;
    tmp = mymin(array3, array4);
    array4 = mymax(array3, array4);
    array3 = tmp;
    tmp = mymin(array6, array7);
    array7 = mymax(array6, array7);
    array6 = tmp;
    tmp = mymin(array1, array2);
    array2 = mymax(array1, array2);
    array1 = tmp;
    tmp = mymin(array4, array5);
    array5 = mymax(array4, array5);
    array4 = tmp;
    tmp = mymin(array7, array8);
    array8 = mymax(array7, array8);
    array3 = mymax(array0, array3);
    array5 = mymin(array5, array8);
    array7 = mymax(array4, tmp);
    tmp = mymin(array4, tmp);
    array6 = mymax(array3, array6);
    array4 = mymax(array1, tmp);
    array2 = mymin(array2, array5);
    array4 = mymin(array4, array7);
    tmp = mymin(array4, array2);
    array2 = mymax(array4, array2);
    array4 = mymax(array6, tmp);
    return mymin(array4, array2);
}

#ifdef __SSE2__
inline __attribute__((always_inline)) __m128 median9sse(__m128 array0, __m128 array1, __m128 array2, __m128 array3, __m128 array4, __m128 array5, __m128 array6, __m128 array7, __m128 array8)
{
    __m128 tmp = _mm_min_ps(array1, array2);
    array2 = _mm_max_ps(array1, array2);
    array1 = tmp;
    tmp = _mm_min_ps(array4, array5);
    array5 = _mm_max_ps(array4, array5);
    array4 = tmp;
    tmp = _mm_min_ps(array7, array8);
    array8 = _mm_max_ps(array7, array8);
    array7 = tmp;
    tmp = _mm_min_ps(array0, array1);
    array1 = _mm_max_ps(array0, array1);
    array0 = tmp;
    tmp = _mm_min_ps(array3, array4);
    array4 = _mm_max_ps(array3, array4);
    array3 = tmp;
    tmp = _mm_min_ps(array6, array7);
    array7 = _mm_max_ps(array6, array7);
    array6 = tmp;
    tmp = _mm_min_ps(array1, array2);
    array2 = _mm_max_ps(array1, array2);
    array1 = tmp;
    tmp = _mm_min_ps(array4, array5);
    array5 = _mm_max_ps(array4, array5);
    array4 = tmp;
    tmp = _mm_min_ps(array7, array8);
    array8 = _mm_max_ps(array7, array8);
    array3 = _mm_max_ps(array0, array3);
    array5 = _mm_min_ps(array5, array8);
    array7 = _mm_max_ps(array4, tmp);
    tmp = _mm_min_ps(array4, tmp);
    array6 = _mm_max_ps(array3, array6);
    array4 = _mm_max_ps(array1, tmp);
    array2 = _mm_min_ps(array2, array5);
    array4 = _mm_min_ps(array4, array7);
    tmp = _mm_min_ps(array4, array2);
    array2 = _mm_max_ps(array4, array2);
    array4 = _mm_max_ps(array6, tmp);
    return _mm_min_ps(array4, array2);
}

__m128 median3x3sse(__m128* array);
__m128 median5x5sse(__m128* array);
__m128 median7x7sse(__m128* array);
__m128 median9x9sse(__m128* array);
#endif
