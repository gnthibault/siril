#include "../core/siril.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define USE_ALL_SORTING_ALGOS
#include "../algos/sorting.h"

/* This program tests the new implementation of the quickselect from Emmanuel
 * and the old from statistics.
 * It compares the results with the quicksort*/

#define NBTRIES 200	// for result checking, unit test of implementations

double median_from_sorted_array(WORD *arr, int size) {
	if (size % 2)
		return arr[(size-1)/2];
	int sum = (int)arr[(size-1)/2] + (int)arr[size/2];
	return (double)sum/2.0;
}

int test_quickselect(int datasize) {
	WORD *data1, *data2, *data3;
	double result_qsel1, result_qsel2, result_qsort;
	int i, retval = 0;

	data1 = malloc(datasize * sizeof(WORD));
	data2 = malloc(datasize * sizeof(WORD));
	data3 = malloc(datasize * sizeof(WORD));
	for (i=0; i<datasize; i++) {
		int val = rand() % USHRT_MAX;
		data1[i] = (WORD)val;
		data2[i] = (WORD)val;
		data3[i] = (WORD)val;
	}

	quicksort_s(data1, datasize);
	result_qsort = median_from_sorted_array(data1, datasize);
	result_qsel1 = quickmedian(data2, datasize);
	result_qsel2 = histogram_median(data3, datasize);

	if (result_qsel1 != result_qsort || result_qsel2 != result_qsort) {
		for (i=0; i<datasize; i++)
			fprintf(stdout, "%hu ", data1[i]);
		fputc('\n', stdout);
		fprintf(stdout, "got %g (qsel1), %g (qsel2) and %g (qsort)\n",
				result_qsel1, result_qsel2, result_qsort);
		retval = 1;
	}
	else {
		//fprintf(stdout, "OK (size %d) got %hu (qsel1), %hu (qsel2) and %hu (qsort)\n",
		//		datasize, result_qsel1, result_qsel2, result_qsort);
	}
	free(data1);
	free(data2);
	free(data3);
	return retval;
}

void measure_big() {
	WORD *data1, *data2, *data3, result_qsel1, result_qsel2, result_qsort;
	int i, datasize = 30000000;
	clock_t t1, t2, t3, t4;

	data1 = malloc(datasize * sizeof(WORD));
	data2 = malloc(datasize * sizeof(WORD));
	data3 = malloc(datasize * sizeof(WORD));
	for (i=0; i<datasize; i++) {
		int val = rand() % USHRT_MAX;
		data1[i] = (WORD)val;
		data2[i] = (WORD)val;
		data3[i] = (WORD)val;
	}

	t1 = clock();
	quicksort_s(data1, datasize);
	result_qsort = median_from_sorted_array(data1, datasize);
	t2 = clock();
	result_qsel1 = quickmedian(data2, datasize);
	t3 = clock();
	result_qsel2 = histogram_median(data3, datasize);
	t4 = clock();
	fprintf(stdout, "siril quicksort time:\t%ld\n", t2-t1);
	fprintf(stdout, "quickmedian time:\t%ld\n", t3-t2);
	fprintf(stdout, "histogram_median time:\t%ld\n", t4-t3);

	free(data1);
	free(data2);
	free(data3);
}

void measure_small() {
	WORD *data, *data_backup, result_qsel1, result_qsel2, result_qsort;
	int i, draws, times, datasize = 8, nb_draws = 100, nb_times_each = 200000;
	clock_t t1, t2, t3, t4;

	data = malloc(datasize * sizeof(WORD));
	data_backup = malloc(datasize * sizeof(WORD));
	t1 = clock();

	for (draws = 0; draws < nb_draws; draws++) {
		for (i=0; i<datasize; i++) {
			int val = rand() % USHRT_MAX;
			data[i] = (WORD)val;
			data_backup[i] = (WORD)val;
		}

		for (times = 0; times < nb_times_each; times++) {
			quicksort_s(data, datasize);
			result_qsort = median_from_sorted_array(data, datasize);

			for (i=0; i<datasize; i++)
				data[i] = data_backup[i];
			data[times % datasize] = times % USHRT_MAX;
		}
	}
	t2 = clock();

	for (draws = 0; draws < nb_draws; draws++) {
		for (i=0; i<datasize; i++) {
			int val = rand() % USHRT_MAX;
			data[i] = (WORD)val;
			data_backup[i] = (WORD)val;
		}

		for (times = 0; times < nb_times_each; times++) {
			result_qsel1 = quickmedian(data, datasize);

			for (i=0; i<datasize; i++)
				data[i] = data_backup[i];
			data[times % datasize] = times % USHRT_MAX;
		}
	}
	t3 = clock();

	for (draws = 0; draws < nb_draws; draws++) {
		for (i=0; i<datasize; i++) {
			int val = rand() % USHRT_MAX;
			data[i] = (WORD)val;
			data_backup[i] = (WORD)val;
		}

		for (times = 0; times < nb_times_each; times++) {
			result_qsel2 = histogram_median(data, datasize);

			for (i=0; i<datasize; i++)
				data[i] = data_backup[i];
			data[times % datasize] = times % USHRT_MAX;
		}
	}
	t4 = clock();

	fprintf(stdout, "siril quicksort time:\t%ld\n", t2-t1);
	fprintf(stdout, "quickmedian time:\t%ld\n", t3-t2);
	fprintf(stdout, "histogram_median time:\t%ld\n", t4-t3);

	free(data);
	free(data_backup);
}

int main(void)
{
	int i, size = 1;
	srand(time(NULL));
	for (i = 0; i < NBTRIES; i++, size++) {
		if (test_quickselect(size)) {
			fprintf(stderr, "FAILED\n");
			exit(1);
		}
	}

	measure_big();
	measure_small();
	return 0;
}


