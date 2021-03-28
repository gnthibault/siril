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

#include <criterion/criterion.h>
#include <gsl/gsl_statistics_float.h>
#include <gsl/gsl_cdf.h>

#include "core/siril.h"
#include "stacking/rejection_float.c"

cominfo com;	// the main data struct
GtkBuilder *builder = NULL;	// get widget references anywhere
fits gfit;	// currently loaded image

float y[] = { 145, 125, 190, 135, 220, 130, 210, 3, 165, 165, 150, 350, 170, 180, 195, 440, 215, 135, 410, 40, 140, 175 };
guint64 count[2] = { 0, 0 };

static void print_outliers(struct outliers *rej, int N) {
	fprintf(stderr, "outliers are: ");
	for (int i = 0; i < N; i++) {
		if (rej[i].out) printf("%.6f ", rej[i].x);
	}
	fprintf(stderr, "\n");
}

static float calculate_critical_value(int size, float alpha) {
	float t_dist = gsl_cdf_tdist_Pinv(1 - alpha / (2 * size), size - 2);
	float numerator = (size - 1) * t_dist;
	float denominator = sqrtf(size) * sqrtf(size - 2 + (t_dist * t_dist));
	return numerator / denominator;
}

static void ESD_test(float *stack, int size, float alpha, int max_outliers) {
	struct outliers *rej = malloc(max_outliers * sizeof(struct outliers));
	int *rejected = malloc(size * sizeof(int));

	quicksort_f(stack, size);
	double median = gsl_stats_float_median_from_sorted_data(stack, 1, size);

	for (int iter = 0; iter < max_outliers; iter++) {
		float Gstat, Gcritical;
		int max_index = 0;

		Gcritical = calculate_critical_value(size, alpha);
		grubbs_stat(stack, size, &Gstat, &max_index);
		rej[iter].out = check_G_values(Gstat, Gcritical);
		rej[iter].x = stack[max_index];
		rej[iter].i = max_index;
		remove_element(stack, max_index, size);
		size--;
	}
	confirm_outliers(rej, max_outliers, median, rejected, count);
	print_outliers(rej, max_outliers);

	cr_expect_eq(count[0], 2);
	cr_expect_eq(count[1], 3);

	cr_expect_float_eq(rej[0].x, 440.0f, 1e-6);
	cr_expect_float_eq(rej[1].x, 410.0f, 1e-6);
	cr_expect_float_eq(rej[2].x, 350.0f, 1e-6);
	cr_expect_float_eq(rej[3].x, 3.000f, 1e-6);
	cr_expect_float_eq(rej[4].x, 40.00f, 1e-6);

	free(rejected);
	free(rej);
}

void test_GESDT_float() {
	ESD_test(y, sizeof(y)/sizeof(float), 0.05, 7);
}

Test(science, psf_float) { test_GESDT_float(); }
