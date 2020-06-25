#include <criterion/criterion.h>

#include "../core/siril.h"

#define USE_ALL_SORTING_ALGOS
#include "../algos/sorting.h"

#include <stdlib.h>
#include <time.h>

/* This program tests the new implementation of the quickselect from Emmanuel
 * and the old from statistics.
 * It compares the results with the quicksort*/

#define NBTRIES 200	// for result checking, unit test of implementations
#define USE_MULTITHREADING TRUE

cominfo com;	// the main data struct

double median_from_sorted_array(WORD *arr, int size)
{
	if (size % 2)
		return arr[(size-1)/2];
	int sum = (int)arr[(size-1)/2] + (int)arr[size/2];
	return (double)sum/2.0;
}

int compare_median_algos(int datasize)
{
	WORD *data, *data_backup;
	double result_qsel1, result_qsel2, result_qsort;
	int i, retval = 0;

	data = malloc(datasize * sizeof(WORD));
	data_backup = malloc(datasize * sizeof(WORD));
	for (i=0; i<datasize; i++) {
		int val = rand() % USHRT_MAX;
		data[i] = (WORD)val;
		data_backup[i] = (WORD)val;
	}

	quicksort_s(data, datasize);
	result_qsort = median_from_sorted_array(data, datasize);
	memcpy(data_backup, data, datasize * sizeof(WORD));

	result_qsel1 = quickmedian(data, datasize);
	memcpy(data_backup, data, datasize * sizeof(WORD));

	result_qsel2 = histogram_median(data, datasize, USE_MULTITHREADING);
	memcpy(data_backup, data, datasize * sizeof(WORD));

	if (result_qsel1 != result_qsort || result_qsel2 != result_qsort) {
		cr_log_error("got %g (quickmedian), %g (histogram_median) and %g (qsort)\n",
					 result_qsel1, result_qsel2, result_qsort);
		retval = 1;
	}

	free(data);
	free(data_backup);
	return retval;
}

void common_setup()
{
	srand(time(NULL));
	com.max_thread = g_get_num_processors();
}

TestSuite(Sorting, .init=common_setup);

Test(Sorting, Median)
{
	int size = 0;
	for (int i = 0; i < NBTRIES; i++, size++) {
		cr_assert(compare_median_algos(size) == 0, "Failed at size=%u", size);
	}
}
