#ifndef SORTING_H_
#define SORTING_H_

#include "core/siril.h"	// for types

/* the quicksorts */
void quicksort_d (double *a, size_t n);
void quicksort_f (float *a, size_t n);
void quicksort_s (WORD *a, size_t n);

/* Quick median based on quick select */
double quickmedian (WORD *a, size_t n);
double quickmedian_float (float *a, size_t n);
double quickmedian_double(double *a, size_t n);
double quickmedian_int (int *a, size_t n);

/* Histogram median for very large array of unsigned short */
double histogram_median (WORD *a, size_t n, gboolean multithread);
double histogram_median_float(float *a, size_t n, gboolean multithread);

/* Sorting netnork */
double sortnet_median (WORD *a, size_t n);
double sortnet_median_double (double *a, size_t n);
void sortnet (WORD *a, size_t n);

gint strcompare(gconstpointer *a, gconstpointer *b);

#endif
