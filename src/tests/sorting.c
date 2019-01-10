#include "../core/siril.h"
#include <stdio.h>
#include <stdlib.h>

#define NBTRIES 200

WORD quickselect_s(WORD *a, int n, int k);
void quicksort_s(WORD *a, int n);

int test_quickselect(int datasize) {
	WORD *data, result_qsel, result_qsort;
	int i;

	data = malloc(datasize * sizeof(WORD));
	for (i=0; i<datasize; i++) {
		int val = rand() % USHRT_MAX;
		data[i] = (WORD)val;
	}

	result_qsel = quickselect_s(data, datasize, datasize/2);
	quicksort_s(data, datasize);
	result_qsort = data[datasize/2];
	free(data);

	if (result_qsel != result_qsort) {
		for (i=0; i<datasize; i++)
			fprintf(stdout, "%hu ", data[i]);
		fputc('\n', stdout);
		fprintf(stdout, "got %hu (qsel) and %hu (qsort)\n", result_qsel, result_qsort);
		return 1;
	}
	fprintf(stdout, "OK (size %d) got %hu (qsel) and %hu (qsort)\n", datasize, result_qsel, result_qsort);
	return 0;
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
	return 0;
}


