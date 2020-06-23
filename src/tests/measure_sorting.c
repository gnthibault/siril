#include "../core/siril.h"

#define USE_ALL_SORTING_ALGOS
#include "../algos/sorting.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define USE_MULTITHREADING TRUE

cominfo com;	// the main data struct

double median_from_sorted_array(WORD *arr, int size)
{
	if (size % 2)
		return arr[(size-1)/2];
	int sum = (int)arr[(size-1)/2] + (int)arr[size/2];
	return (double)sum/2.0;
}

double _siril_qsort(WORD *data, size_t datasize)
{
	quicksort_s(data, datasize);
	return median_from_sorted_array(data, datasize);
}

double _histogram_sort(WORD *data, size_t datasize)
{
	return histogram_median(data, datasize, USE_MULTITHREADING);
}

clock_t perf_test(double (*function)(WORD *data, size_t datasize),
				  int datasize,
				  int nb_draws,
				  int nb_times_each)
{
	WORD *data = malloc(datasize * sizeof(WORD));
	WORD *data_backup = malloc(datasize * sizeof(WORD));
	int i, draws, times;

	clock_t t_start = clock();
	for (draws = 0; draws < nb_draws; draws++) {
		for (i=0; i<datasize; i++) {
			int val = rand() % USHRT_MAX;
			data[i] = (WORD)val;
			data_backup[i] = (WORD)val;
		}

		for (times = 0; times < nb_times_each; times++) {
			function(data, datasize);

			memcpy(data, data_backup, datasize * sizeof(WORD));
			data[times % datasize] = times % USHRT_MAX;
		}
	}
	clock_t t_end = clock();

	free(data);
	free(data_backup);

	return t_end - t_start;;
}

void MeasureSmall()
{
	int datasize = 8;
	int nb_draws = 100;
	int nb_times_each = 200000;

	fprintf(stdout, "== small dataset (%d elements, %d different draws run %d times)\n",
			datasize, nb_draws, nb_times_each);

	clock_t t_siril = perf_test(_siril_qsort, datasize, nb_draws, nb_times_each);
	clock_t t_quick = perf_test(quickmedian, datasize, nb_draws, nb_times_each);
	clock_t t_hist = perf_test(_histogram_sort, datasize, nb_draws, nb_times_each);

	fprintf(stdout, "siril quicksort time:\t%ld\n", t_siril);
	fprintf(stdout, "quickmedian time:\t%ld\n", t_quick);
	fprintf(stdout, "histogram_median time:\t%ld\n", t_hist);
}

void MeasureBig()
{
	int datasize = 30000000;

	fprintf(stdout, "== large dataset (%d elements, same for each)\n", datasize);

	clock_t t_siril = perf_test(_siril_qsort, datasize, 1, 1);
	clock_t t_quick = perf_test(quickmedian, datasize, 1, 1);
	clock_t t_hist = perf_test(_histogram_sort, datasize, 1, 1);

	fprintf(stdout, "siril quicksort time:\t%ld\n", t_siril);
	fprintf(stdout, "quickmedian time:\t%ld\n", t_quick);
	fprintf(stdout, "histogram_median time:\t%ld\n", t_hist);
}

int main()
{
	srand(time(NULL));
	com.max_thread = g_get_num_processors();

	fputc('\n', stdout);
	MeasureBig();
	fputc('\n', stdout);
	MeasureSmall();
	fputc('\n', stdout);
	return 0;
}
