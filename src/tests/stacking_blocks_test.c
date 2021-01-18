//#define WITH_MAIN
#ifndef WITH_MAIN
#include <criterion/criterion.h>
#endif

#include "../core/siril.h"
#include "stacking/stacking.h"
#include <stdio.h>

cominfo com;	// the main data struct

/* This is the test file for the function that allocates data block to thread
 * for siril median or mean stacking.
 */

/* inputs : image 3D dimensions, number of images, channel depth, configured
 *	memory limits, actual memory available, configured thread limit
 * outputs : number of blocks and their size
 */

#ifdef WITH_MAIN
#define CHECK(cond, ...) \
	if (cond) { \
		fprintf(stderr, __VA_ARGS__); \
		return 1; \
	}
#else
#define CHECK cr_expect
#endif

int check_that_blocks_cover_the_image(long naxes[3], struct _image_block *blocks, int nb_blocks) {
	// assuming that blocks are adjacent 
	long ok_up_to_y = 0L, ok_up_to_z = 0L;
	for (int i = 0; i < nb_blocks; i++) {
		if (blocks[i].start_row != ok_up_to_y)
			return 0;
		ok_up_to_y = blocks[i].end_row + 1L;
		if (ok_up_to_y == naxes[1]) {
			ok_up_to_y = 0L;
			ok_up_to_z++;
		}
	}
	return ok_up_to_y == 0L && ok_up_to_z == naxes[2];
}

/* test  chans       memory          threads
 *   1     1      enough                1
 *   2     1      enough                8
 *   3     1      not enough            1
 *   4     1      not enough            8
 *   5     3      enough                1
 *   6     3      enough                8
 *   7     3      not enough for 3      1
 *   8     3      not enough for 1      1
 *   9     3      not enough for 3      8
 *  10     3      not enough for 2      8
 *  11     3      not enough for 1      8
 *  12     3      not enough for 2     12
 *
 */

int test1() {
	// intputs
	long naxes[] = { 1000L, 1000L, 1L };
	long max_rows = 1001L;
	int nb_threads = 1;
	// outputs
	struct _image_block *blocks = NULL;
	int retval, nb_blocks = -1;
	long largest_block = -1;

	/* case 1: single channel images, enough memory, one thread */
	retval = stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block, &nb_blocks);
	CHECK(retval, "retval indicates function failed\n");
	CHECK(nb_blocks != 1 && nb_blocks != 2, "number of blocks returned is %d (expected 1 or 2)\n", nb_blocks);
	CHECK(!blocks, "blocks is null\n");
	CHECK(!check_that_blocks_cover_the_image(naxes, blocks, nb_blocks), "blocks don't cover the whole image\n");
	CHECK(largest_block * nb_threads > max_rows, "this is solution is going out of memory\n");
	fprintf(stdout, "* test 1 passed *\n");
	return 0;
}

int test2() {
	// intputs
	long naxes[] = { 1000L, 1000L, 1L };
	long max_rows = 1001L;
	int nb_threads = 8;
	// outputs
	struct _image_block *blocks = NULL;
	int retval, nb_blocks = -1;
	long largest_block = -1;

	/* case 2: single channel images, enough memory, 8 threads */
	retval = stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block, &nb_blocks);
	CHECK(retval, "retval indicates function failed\n");
	CHECK(nb_blocks != 8, "number of blocks returned is %d (expected 8)\n", nb_blocks);
	CHECK(!blocks, "blocks is null\n");
	CHECK(!check_that_blocks_cover_the_image(naxes, blocks, nb_blocks), "blocks don't cover the whole image\n");
	CHECK(largest_block * nb_threads > max_rows, "this is solution is going out of memory\n");
	fprintf(stdout, "* test 2 passed *\n");
	return 0;
}

int test3() {
	// intputs
	long naxes[] = { 1000L, 1000L, 1L };
	long max_rows = 999L;
	int nb_threads = 1;
	// outputs
	struct _image_block *blocks = NULL;
	int retval, nb_blocks = -1;
	long largest_block = -1;


	/* case 3: single channel images, not enough memory, 1 thread */
	retval = stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block, &nb_blocks);
	CHECK(retval, "retval indicates function failed\n");
	CHECK(nb_blocks != 2, "number of blocks returned is %d (expected 2)\n", nb_blocks);
	CHECK(!blocks, "blocks is null\n");
	CHECK(!check_that_blocks_cover_the_image(naxes, blocks, nb_blocks), "blocks don't cover the whole image\n");
	CHECK(largest_block * nb_threads > max_rows, "this is solution is going out of memory\n");
	fprintf(stdout, "* test 3 passed *\n");
	return 0;
}

int test4() {
	// intputs
	long naxes[] = { 1000L, 1000L, 1L };
	long max_rows = 999L;
	int nb_threads = 8;
	// outputs
	struct _image_block *blocks = NULL;
	int retval, nb_blocks = -1;
	long largest_block = -1;

	/* case 4: single channel images, not enough memory, 8 threads */
	retval = stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block, &nb_blocks);
	CHECK(retval, "retval indicates function failed\n");
	CHECK(nb_blocks <= 8, "number of blocks returned is %d (expected more than 8, ideally 16)\n", nb_blocks);
	CHECK(!blocks, "blocks is null\n");
	CHECK(!check_that_blocks_cover_the_image(naxes, blocks, nb_blocks), "blocks don't cover the whole image\n");
	CHECK(largest_block * nb_threads > max_rows, "this is solution is going out of memory\n");
	fprintf(stdout, "* test 4 passed *\n");
	return 0;
}

int test5() {
	// intputs
	long naxes[] = { 1000L, 1000L, 3L };
	long max_rows = 3001L;
	int nb_threads = 1;
	// outputs
	struct _image_block *blocks = NULL;
	int retval, nb_blocks = -1;
	long largest_block = -1;

	/* case 5: three channel images, enough memory, one thread */
	retval = stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block, &nb_blocks);
	CHECK(retval, "retval indicates function failed\n");
	CHECK(nb_blocks != 3, "number of blocks returned is %d (expected 3)\n", nb_blocks);
	CHECK(!blocks, "blocks is null\n");
	CHECK(!check_that_blocks_cover_the_image(naxes, blocks, nb_blocks), "blocks don't cover the whole image\n");
	CHECK(largest_block * nb_threads > max_rows, "this is solution is going out of memory\n");
	fprintf(stdout, "* test 5 passed *\n");
	return 0;
}

int test6() {
	// intputs
	long naxes[] = { 1000L, 1000L, 3L };
	long max_rows = 3001L;
	int nb_threads = 8;
	// outputs
	struct _image_block *blocks = NULL;
	int retval, nb_blocks = -1;
	long largest_block = -1;

	/* case 6: three channel images, enough memory, 8 threads */
	retval = stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block, &nb_blocks);
	CHECK(retval, "retval indicates function failed\n");
	CHECK(nb_blocks < 6, "number of blocks returned is %d (expected at least 6)\n", nb_blocks);
	CHECK(!blocks, "blocks is null\n");
	CHECK(!check_that_blocks_cover_the_image(naxes, blocks, nb_blocks), "blocks don't cover the whole image\n");
	CHECK(largest_block * nb_threads > max_rows, "this is solution is going out of memory\n");
	fprintf(stdout, "* test 6 passed *\n");
	return 0;
}

int test7() {
	// intputs
	long naxes[] = { 1000L, 1000L, 3L };
	long max_rows = 2999L;
	int nb_threads = 1;
	// outputs
	struct _image_block *blocks = NULL;
	int retval, nb_blocks = -1;
	long largest_block = -1;


	/* case 7: three channel images, not enough memory for all, 1 thread */
	retval = stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block, &nb_blocks);
	CHECK(retval, "retval indicates function failed\n");
	CHECK(nb_blocks != 3, "number of blocks returned is %d (expected 3)\n", nb_blocks);
	CHECK(!blocks, "blocks is null\n");
	CHECK(!check_that_blocks_cover_the_image(naxes, blocks, nb_blocks), "blocks don't cover the whole image\n");
	CHECK(largest_block * nb_threads > max_rows, "this is solution is going out of memory\n");
	fprintf(stdout, "* test 7 passed *\n");
	return 0;
}

int test8() {
	// intputs
	long naxes[] = { 1000L, 1000L, 3L };
	long max_rows = 999L;
	int nb_threads = 1;
	// outputs
	struct _image_block *blocks = NULL;
	int retval, nb_blocks = -1;
	long largest_block = -1;

	/* case 8: three channel images, not enough memory for one, 1 thread */
	retval = stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block, &nb_blocks);
	CHECK(retval, "retval indicates function failed\n");
	CHECK(nb_blocks != 6, "number of blocks returned is %d (expected 6)\n", nb_blocks);
	CHECK(!blocks, "blocks is null\n");
	CHECK(!check_that_blocks_cover_the_image(naxes, blocks, nb_blocks), "blocks don't cover the whole image\n");
	CHECK(largest_block * nb_threads > max_rows, "this is solution is going out of memory\n");
	fprintf(stdout, "* test 8 passed *\n");
	return 0;
}

int test9() {
	// intputs
	long naxes[] = { 1000L, 1000L, 3L };
	long max_rows = 2999L;
	int nb_threads = 8;
	// outputs
	struct _image_block *blocks = NULL;
	int retval, nb_blocks = -1;
	long largest_block = -1;

	/* case 9: three channel images, not enough memory for all, 8 threads */
	/* best solutions are 15 * 200 or 24 * 125 */
	retval = stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block, &nb_blocks);
	CHECK(retval, "retval indicates function failed\n");
	CHECK(nb_blocks <= 14, "number of blocks returned is %d (expected more than 14)\n", nb_blocks);
	CHECK(!blocks, "blocks is null\n");
	CHECK(!check_that_blocks_cover_the_image(naxes, blocks, nb_blocks), "blocks don't cover the whole image\n");
	CHECK(largest_block * nb_threads > max_rows, "this is solution is going out of memory\n");
	fprintf(stdout, "* test 9 passed *\n");
	return 0;
}

int test10() {
	// intputs
	long naxes[] = { 1000L, 1000L, 3L };
	long max_rows = 1200L;
	int nb_threads = 8;
	// outputs
	struct _image_block *blocks = NULL;
	int retval, nb_blocks = -1;
	long largest_block = -1;

	/* case 10: three channel images, not enough memory for two, 8 threads */
	/* best solutions are 21 * 143 or 24 * 125 */
	retval = stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block, &nb_blocks);
	CHECK(retval, "retval indicates function failed\n");
	CHECK(nb_blocks <= 20, "number of blocks returned is %d (expected more than 20)\n", nb_blocks);
	CHECK(!blocks, "blocks is null\n");
	CHECK(!check_that_blocks_cover_the_image(naxes, blocks, nb_blocks), "blocks don't cover the whole image\n");
	CHECK(largest_block * nb_threads > max_rows, "this is solution is going out of memory\n");
	fprintf(stdout, "* test 10 passed *\n");
	return 0;
}

int test11() {
	// intputs
	long naxes[] = { 1000L, 1000L, 3L };
	long max_rows = 999L;
	int nb_threads = 8;
	// outputs
	struct _image_block *blocks = NULL;
	int retval, nb_blocks = -1;
	long largest_block = -1;

	/* case 11: three channel images, not enough memory for one, 8 threads */
	/* best solution is 30 * 100, 27 * 111 is not good for threads */
	retval = stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block, &nb_blocks);
	CHECK(retval, "retval indicates function failed\n");
	CHECK(nb_blocks <= 29, "number of blocks returned is %d (expected more than 29)\n", nb_blocks);
	CHECK(!blocks, "blocks is null\n");
	CHECK(!check_that_blocks_cover_the_image(naxes, blocks, nb_blocks), "blocks don't cover the whole image\n");
	CHECK(largest_block * nb_threads > max_rows, "this is solution is going out of memory\n");
	fprintf(stdout, "* test 11 passed *\n");
	return 0;
}

int test12() {
	// intputs
	long naxes[] = { 6024L, 4024L, 3L };
	int nb_threads = 12;
	int nb_images = 209;
	long max_rows = 27295481856 / (nb_images * naxes[0] * 4); // 27295481856 is available mem

	// outputs
	struct _image_block *blocks = NULL;
	int retval, nb_blocks = -1;
	long largest_block = -1;

	/* case 12: real data that triggered the bug: 209 big images, 32 GB of
	 * memory, 12 threads, that's enough memory to process a bit more than one
	 * channel at time, similar to test 10.
	 * Typical solution: 33 parallel blocks of size 365 (+9)
	 */
	retval = stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block, &nb_blocks);
	CHECK(retval, "retval indicates function failed\n");
	CHECK(nb_blocks < 33, "number of blocks returned is %d (expected at least 33)\n", nb_blocks);
	CHECK(!blocks, "blocks is null\n");
	CHECK(!check_that_blocks_cover_the_image(naxes, blocks, nb_blocks), "blocks don't cover the whole image\n");
	CHECK(largest_block * nb_threads > max_rows, "this is solution is going out of memory\n");
	fprintf(stdout, "* test 12 passed *\n");
	return 0;
}

#ifdef WITH_MAIN
int main() {
	int retval = 0;
	retval |= test1();
	retval |= test2();
	retval |= test3();
	retval |= test4();
	retval |= test5();
	retval |= test6();
	retval |= test7();
	retval |= test8();
	retval |= test9();
	retval |= test10();
	retval |= test11();
	retval |= test12();
	if (retval)
		fprintf(stderr, "TESTS FAILED\n");
	else fprintf(stderr, "ALL TESTS PASSED\n");
	return retval;
}
#else //with criterion

Test(stacking_blocks, test1) { cr_assert(test1()); }
Test(stacking_blocks, test2) { cr_assert(test2()); }
Test(stacking_blocks, test3) { cr_assert(test3()); }
Test(stacking_blocks, test4) { cr_assert(test4()); }
Test(stacking_blocks, test5) { cr_assert(test5()); }
Test(stacking_blocks, test6) { cr_assert(test6()); }
Test(stacking_blocks, test7) { cr_assert(test7()); }
Test(stacking_blocks, test8) { cr_assert(test8()); }
Test(stacking_blocks, test9) { cr_assert(test9()); }
Test(stacking_blocks, test10) { cr_assert(test10()); }
Test(stacking_blocks, test11) { cr_assert(test11()); }
Test(stacking_blocks, test12) { cr_assert(test12()); }

#endif
