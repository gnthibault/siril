#include "../core/siril.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* This program tests the new implementation of the quickselect from Emmanuel
 * and the old from statistics.
 * It compares the results with the quicksort*/

#define NBTRIES 200

WORD quickselect_s(WORD *a, int n, int k);
void quicksort_s(WORD *a, int n);
WORD siril_stats_ushort_median(WORD *arr, int n);	// from algos/statistics.c

int test_quickselect(int datasize) {
	WORD *data1, *data2, *data3, result_qsel1, result_qsel2, result_qsort;
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
	result_qsort = data1[(datasize-1)/2];
	result_qsel1 = quickselect_s(data2, datasize, (datasize-1)/2);
	result_qsel2 = siril_stats_ushort_median(data3, datasize);

	if (result_qsel1 != result_qsort || result_qsel2 != result_qsort) {
		for (i=0; i<datasize; i++)
			fprintf(stdout, "%hu ", data1[i]);
		fputc('\n', stdout);
		fprintf(stdout, "got %hu (qsel1), %hu (qsel2) and %hu (qsort)\n",
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

void measure_quickselect() {
	WORD *data1, *data2, *data3, result_qsel1, result_qsel2, result_qsort;
	int i, datasize = 40000000;
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
	result_qsort = data1[(datasize-1)/2];
	t2 = clock();
	result_qsel1 = quickselect_s(data2, datasize, (datasize-1)/2);
	t3 = clock();
	result_qsel2 = siril_stats_ushort_median(data3, datasize);
	t4 = clock();
	fprintf(stdout, "siril quicksort time:\t%ld\n", t2-t1);
	fprintf(stdout, "quicksselect_s time:\t%ld\n", t3-t2);
	fprintf(stdout, "stats_median time:\t%ld\n", t4-t3);

	free(data1);
	free(data2);
	free(data3);
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

	measure_quickselect();
	return 0;
}


