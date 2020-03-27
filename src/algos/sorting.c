#include <string.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"

#include "rt/rt_algo.h"
#include "sorting.h"


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

/* This files contains all sorting algorithms of siril.
 * See src/test/sorting.c for testing and metrics.
 */

/**
 * In-place insertion sort of array of double a of size n
 * @param a array to sort
 * @param n size of the array
 */
static void insertionSort_d(double a[], size_t n) {
	for (size_t i = 1; i < n; i++) {
		const double val = a[i];
		size_t j = i - 1;

		/* Move elements of a[0..i-1], that are greater than val, to one position ahead of their current position */
		while (j >= 0 && a[j] > val) {
			a[j + 1] = a[j];
			j = j - 1;
		}
		a[j + 1] = val;
	}
}

/**
 * In-place quick sort of array of double a of size n
 * @param a array to sort
 * @param n size of the array
 */
void quicksort_d (double *a, size_t n) {
	if (n <= 32) {
		return insertionSort_d(a, n);
	}

	double pivot = a[n / 2];
	double *left = a;
	double *right = a + n - 1;
	register double t;

	while (left <= right) {
		if (*left < pivot) {
			left++;
			continue;
		}
		if (*right > pivot) {
			right--;
			continue;
		}
		t = *left;
		*left++ = *right;
		*right-- = t;
	}
	quicksort_d(a, right - a + 1);
	quicksort_d(left, a + n - left);
}

/**
 * In-place insertion sort of array of float a of size n
 * @param a array to sort
 * @param n size of the array
 */
 static void insertionSort_f(float a[], size_t n) {
	for (size_t i = 1; i < n; i++) {
		const float val = a[i];
		size_t j = i - 1;

		/* Move elements of a[0..i-1], that are greater than val, to one position ahead of their current position */
		while (j >= 0 && a[j] > val) {
			a[j + 1] = a[j];
			j = j - 1;
		}
		a[j + 1] = val;
	}
}

/**
 * In-place quick sort of array of float a of size n
 * @param a array to sort
 * @param n size of the array
 */
void quicksort_f (float *a, size_t n) {
	if (n <= 32) {
		return insertionSort_f(a, n);
	}

	float pivot = a[n / 2];
	float *left = a;
	float *right = a + n - 1;
	register float t;

	while (left <= right) {
		if (*left < pivot) {
			left++;
			continue;
		}
		if (*right > pivot) {
			right--;
			continue;
		}
		t = *left;
		*left++ = *right;
		*right-- = t;
	}
	quicksort_f(a, right - a + 1);
	quicksort_f(left, a + n - left);
}

/**
 * In-place insertion sort of array of WORD a of size n
 * @param a array to sort
 * @param n size of the array
 */
static void insertionSort_s(WORD a[], size_t n) {
	for (size_t i = 1; i < n; i++) {
		const WORD val = a[i];
		size_t j = i - 1;

		/* Move elements of a[0..i-1], that are greater than val, to one position ahead of their current position */
		while (j >= 0 && a[j] > val) {
			a[j + 1] = a[j];
			--j;
		}
		a[j + 1] = val;
	}
}

/**
 * In-place quick sort of array of WORD a of size n
 * @param a array to sort
 * @param n size of the array
 */
void quicksort_s(WORD *a, size_t n) {
	if (n <= 32) {
		return insertionSort_s(a, n);
	}
	WORD pivot = a[n / 2];
	WORD *left = a;
	WORD *right = a + n - 1;
	register WORD t;

	while (left <= right) {
		if (*left < pivot) {
			left++;
			continue;
		}
		if (*right > pivot) {
			right--;
			continue;
		}
		t = *left;
		*left++ = *right;
		*right-- = t;
	}
	quicksort_s(a, right - a + 1);
	quicksort_s(left, a + n - left);
}

/* quickmedian returns the median from array of length n
 * Derived from the original quickselect algorithm from Hoare
 * warning: data are sorted in place
 * non recurssive version modified to return median value
 * @param a array of WORD to search
 * @param n size of the array
 * @return median as double for even size average the middle two elements
*/
double quickmedian(WORD *a, size_t n) {
	// Use faster and robust sorting network for small size array
	if (n < 9)
		return sortnet_median(a, n);

	WORD pivot, tmp;
	size_t k = n / 2;	// size to sort
	size_t pindex;		// pivot index
	size_t left = 0; 	// left index
	size_t right = n - 1; 	// right index

	while (left < right) { //we stop when our indicies have crossed
		pindex = (left + right) / 2; // pivot selection, this can be whatever
		pivot = a[pindex];
		a[pindex] = a[right];
		a[right] = pivot; // SWAP(pivot,right)

		for (size_t i = pindex = left; i < right; i++) {
			if (a[i] < pivot) { // SWAP
				tmp = a[pindex];
				a[pindex] = a[i];
				a[i] = tmp;
				pindex++;
			}
		}
		a[right] = a[pindex];
		a[pindex] = pivot; // SWAP

		if (pindex < k)
			left = pindex + 1;
		else
			// pindex >= k
			right = pindex;
	}
	return (n % 2 == 0) ?
			((double) a[k - 1] + (double) a[k]) / 2.0 : (double) a[k];
}

/* quickmedian_float returns the median from array of length n
 * Derived from the original quickselect algorithm from Hoare
 * warning: data are sorted in place
 * non recurssive version modified to return median value
 * @param a array of float to search
 * @param n size of the array
 * @return median as double for even size average the middle two elements
 */
double quickmedian_float(float *a, size_t n) {
	size_t k = n / 2;	// size to sort
	size_t pindex;		// pivot index
	size_t left = 0; 	// left index
	size_t right = n - 1; 	// right index
	float pivot, tmp;

	while (left < right) { //we stop when our indicies have crossed
		pindex = (left + right) / 2; // pivot selection, this can be whatever
		pivot = a[pindex];
		a[pindex] = a[right];
		a[right] = pivot; // SWAP(pivot,right)

		for (size_t i = pindex = left; i < right; i++) {
			if (a[i] < pivot) { // SWAP
				tmp = a[pindex];
				a[pindex] = a[i];
				a[i] = tmp;
				pindex++;
			}
		}
		a[right] = a[pindex];
		a[pindex] = pivot; // SWAP(right,j)

		if (pindex < k)
			left = pindex + 1;
		else
			// pindex >= k
			right = pindex;
	}
	return (n % 2 == 0) ? ((double) a[k - 1] + a[k]) / 2.0 : (double) a[k];
}

/* quickmedian_double returns the median from array of length n
 * Derived from the original quickselect algorithm from Hoare
 * warning: data are sorted in place
 * non recurssive version modified to return median value
 * @param a array of double to search
 * @param n size of the array
 * @return median as double for even size average the middle two elements
 */
double quickmedian_double(double *a, size_t n) {
	size_t k = n / 2;	// size to sort
	size_t pindex;		// pivot index
	size_t left = 0; 	// left index
	size_t right = n - 1; 	// right index
	double pivot, tmp;

	while (left < right) { //we stop when our indicies have crossed
		pindex = (left + right) / 2; // pivot selection, this can be whatever
		pivot = a[pindex];
		a[pindex] = a[right];
		a[right] = pivot; // SWAP(pivot,right)

		for (size_t i = pindex = left; i < right; i++) {
			if (a[i] < pivot) { // SWAP
				tmp = a[pindex];
				a[pindex] = a[i];
				a[i] = tmp;
				pindex++;
			}
		}
		a[right] = a[pindex];
		a[pindex] = pivot; // SWAP(right,j)

		if (pindex < k)
			left = pindex + 1;
		else
			// pindex >= k
			right = pindex;
	}
	return (n % 2 == 0) ? (a[k - 1] + a[k]) / 2.0 : a[k];
}

/*
 * quickmedian_int returns the median from array of int of of length n
 * Derived from original quickselect algorithm from Hoare
 * @param a array to search
 * @param n size of the array
 * @return median as double
 */
double quickmedian_int(int *a, size_t n) {
	size_t k = n / 2;	// size to sort
	size_t pindex;		// pivot index
	size_t left = 0; 	// left index
	size_t right = n - 1; 	// right index
	int pivot, tmp;

	while (left < right) { //we stop when our indicies have crossed
		pindex = (left + right) / 2; // pivot selection, this can be whatever
		pivot = a[pindex];
		a[pindex] = a[right];
		a[right] = pivot; // SWAP

		for (size_t i = pindex = left; i < right; i++) {
			if (a[i] < pivot) { // SWAP
				tmp = a[pindex];
				a[pindex] = a[i];
				a[i] = tmp;
				pindex++;
			}
		}
		a[right] = a[pindex];
		a[pindex] = pivot; // SWAP

		if (pindex < k)
			left = pindex + 1;
		else
			// pindex >= k
			right = pindex;
	}
	return (n % 2 == 0) ?
			((double) a[k - 1] + (double) a[k]) / 2.0 : (double) a[k];
}

/*
 * Optimal sorting network array of size [2,9] to retrieve median value
 * (C) Emmanuel Brandt 2019-02
 * @param a array of double to sort (warning in place sorting)
 * @param n size of the array to sort [2,9]
 * warning in-place sorting
 */
#define sw(i,j) if(a[i] > a[j]) { register WORD t=a[i]; a[i]=a[j]; a[j]=t; }
double sortnet_median (WORD *a, size_t n) {
	size_t k = n / 2;

	switch (n) {
		case 1: return a[0]; break;

		case 2: sw(0,1); break;

		case 3: sw(0,1); sw(1,2); sw(0,1); break;

		case 4: sw(0,1); sw(2,3); sw(0,2); sw(1,3); sw(1,2); break;

		case 5: sw(0,1); sw(2,3); sw(1,3); sw(2,4); sw(0,2);
			sw(1,4); sw(1,2); sw(3,4); sw(2,3); break;

		case 6: sw(0,1); sw(2,3); sw(4,5); sw(0,2); sw(3,5);
			sw(1,4); sw(0,1); sw(2,3); sw(4,5); sw(1,2); sw(3,4);
			sw(2,3); break;

		case 7: sw(1,2); sw(3,4); sw(5,6); sw(0,2); sw(4,6); sw(3,5);
			sw(2,6); sw(1,5); sw(0,4); sw(2,5); sw(0,3); sw(2,4);
			sw(1,3); sw(0,1); sw(2,3); sw(4,5); break;

		case 8: sw(0,1); sw(2,3); sw(4,5); sw(6,7);
			sw(0,2); sw(1,3); sw(4,6); sw(5,7);
			sw(1,2); sw(5,6);
			sw(0,4); sw(1,5); sw(2,6); sw(3,7);
			sw(2,4); sw(3,5);
			sw(1,2); sw(3,4); sw(5,6); break;

		case 9: sw(1,8); sw(2,7); sw(3,6); sw(4,5);
			sw(1,4); sw(5,8);
			sw(0,2); sw(6,7);
			sw(2,6); sw(7,8);
			sw(0,3); sw(4,5);
			sw(0,1); sw(3,5); sw(6,7);
			sw(2,4);
			sw(1,3); sw(5,7);
			sw(4,6);
			sw(1,2); sw(3,4); sw(5,6); sw(7,8);
			sw(2,3); sw(4,5); break;

		default: return 0.0; break; // no sort
	}
	return (n % 2 == 0) ? (a[k - 1] + a[k]) / 2.0 : a[k];
}
#undef sw

/*
 * Optimal sorting network [2,9]
 * (C) Emmanuel Brandt 2019-02
 * @param a array of double to sort (warning in place sorting)
 * @param n size of the array to sort [2,9]
 * warning in-place sorting
 */
#define sw(i,j) if(a[i] > a[j]) { register WORD t=a[i]; a[i]=a[j]; a[j]=t; }
void sortnet (WORD *a, size_t n) {

	switch (n) {
		case 2: sw(0,1); break;

		case 3: sw(0,1); sw(1,2); sw(0,1); break;

		case 4: sw(0,1); sw(2,3); sw(0,2); sw(1,3); sw(1,2); break;

		case 5: sw(0,1); sw(2,3); sw(1,3); sw(2,4); sw(0,2);
			sw(1,4); sw(1,2); sw(3,4); sw(2,3); break;

		case 6: sw(0,1); sw(2,3); sw(4,5); sw(0,2); sw(3,5);
			sw(1,4); sw(0,1); sw(2,3); sw(4,5); sw(1,2); sw(3,4);
			sw(2,3); break;

		case 7: sw(1,2); sw(3,4); sw(5,6); sw(0,2); sw(4,6); sw(3,5);
			sw(2,6); sw(1,5); sw(0,4); sw(2,5); sw(0,3); sw(2,4);
			sw(1,3); sw(0,1); sw(2,3); sw(4,5); break;

		case 8: sw(0,1); sw(2,3); sw(4,5); sw(6,7);
			sw(0,2); sw(1,3); sw(4,6); sw(5,7);
			sw(1,2); sw(5,6);
			sw(0,4); sw(1,5); sw(2,6); sw(3,7);
			sw(2,4); sw(3,5);
			sw(1,2); sw(3,4); sw(5,6); break;

		case 9: sw(1,8); sw(2,7); sw(3,6); sw(4,5);
			sw(1,4); sw(5,8);
			sw(0,2); sw(6,7);
			sw(2,6); sw(7,8);
			sw(0,3); sw(4,5);
			sw(0,1); sw(3,5); sw(6,7);
			sw(2,4);
			sw(1,3); sw(5,7);
			sw(4,6);
			sw(1,2); sw(3,4); sw(5,6); sw(7,8);
			sw(2,3); sw(4,5); break;

		default: break; // no sort
	}
}
#undef sw

/*
 * Histogram median for very large array of unsigned short
 * (C) Emmanuel Brandt 2019-02
 * @param a array of unsigned short to search
 * @param n size of the array
 * @return median as a double (for n odd)
 * Use temp storage h to build the histogram. Complexity O(2*N)
 */
double histogram_median(WORD *a, size_t n, gboolean mutlithread) {
	// For arrays n < 10 histogram is use fast and simple sortnet_median
	if (n < 10)
		return sortnet_median(a, n);

	const size_t s = sizeof(unsigned int);
	unsigned int *h = (unsigned int*) calloc(USHRT_MAX + 1, s);
	if (!h) {
		PRINT_ALLOC_ERR;
		return -1.0;
	}

#ifdef _OPENMP
#pragma omp parallel num_threads(com.max_thread) if (mutlithread)
#endif
	{
		unsigned int *hthr = (unsigned int*) calloc(USHRT_MAX + 1, s);
		if (!hthr) {
			PRINT_ALLOC_ERR;
		}
		else {
#ifdef _OPENMP
#pragma omp for nowait
#endif
			for (size_t i = 0; i < n; i++) {
				hthr[a[i]]++;
			}
#ifdef _OPENMP
#pragma omp critical
#endif
			{
				// add per thread histogram to main histogram
#ifdef _OPENMP
#pragma omp simd
#endif
				for (int ii = 0; ii <= USHRT_MAX; ++ii) {
					h[ii] += hthr[ii];
				}
			}
			free(hthr);
		}
	}
	unsigned int i= 0, j = 0, k = n / 2;

	unsigned int sum = 0;
	if (n % 2 == 0) {
		for (; sum <= k - 1; j++)
			sum += h[j];
		i = j;
	}

	for (; sum <= k; i++)
		sum += h[i];

	free(h);
	return (n % 2 == 0) ? (double) (i + j) / 2.0 : (double) i;
}

double histogram_median_float(float *a, size_t n, gboolean multithread) {
	float median;
	findMinMaxPercentile(a, n, 0.5f, &median, 0.5f, &median, multithread);
	return median;
}

/**
 * Compares a and b like strcmp()
 * @param a a gconstpointer
 * @param b a gconstpointer
 * @return an integer less than, equal to, or greater than zero, if a is than b .
 */
gint strcompare(gconstpointer *a, gconstpointer *b) {
	gchar *collate_key1, *collate_key2;
	gint result;

	const gchar *s1 = (const gchar *) a;
	const gchar *s2 = (const gchar *) b;

	collate_key1 = g_utf8_collate_key_for_filename(s1, strlen(s1));
	collate_key2 = g_utf8_collate_key_for_filename(s2, strlen(s2));

	result = g_strcmp0(collate_key1, collate_key2);
	g_free(collate_key1);
	g_free(collate_key2);

	return result;
}
