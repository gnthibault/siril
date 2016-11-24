/*
 *  match: a package to match lists of stars (or other items)
 *  Copyright (C) 2000  Michael William Richmond
 *
 *  Contact: Michael William Richmond
 *           Physics Department
 *           Rochester Institute of Technology
 *           85 Lomb Memorial Drive
 *           Rochester, NY  14623-5603
 *           E-mail: mwrsps@rit.edu
 *
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */

/*
 * <AUTO>
 * FILE: atpmatch.c
 *
 * <HTML>
 * This file contains routines that try to match up items
 * in two different lists, which might have very different
 * coordinate systems.
 *
 * Stars must have been placed into "s_star" structures before
 * being passed to the functions in this file.  Note that
 * the "x" and "y" fields of an s_star may contain (RA, Dec),
 * or (row, col), or some other coordinates.
 * </HTML>
 *
 * </AUTO>
 *
 */

/*
 * -------------------------------------------------------------------------
 * atFindTrans             public  Find TRANS to match coord systems of
 *                                 two lists of items
 * atApplyTrans            public  Apply a TRANS to a list of items,
 *                                 modifying two of the elements in each item
 * atMatchLists            public  Find all pairs of items on two lists
 *                                 which match (and don't match)
 * atRecalcTrans           public  Calculate a TRANS, given two lists of
 *                                 stars which _already_ have been matched
 * atFindMedtf             public  Calculates MEDTF statistics, assuming
 *                                 that two lists have same scale, rotation
 *
 * All public functions appear at the start of this source-code file.
 * "Private" static functions appear following them.
 *
 * Conditional compilation is controlled by the following macros:
 *
 *    DEBUG
 *    DEBUG2
 *
 * AUTHORS:  SHIVA Creation date: Jan 22, 1996
 *           Michael Richmond
 *
 *           Modified for stand-alone Linux operation: April 26, 1996
 *           Michael Richmond
 *
 *           Added more digits to printout in "print_trans": Aug 18, 1996
 *           Michael Richmond
 *
 *           Changed 'iter_trans()' to reject 3-sigma outliers,
 *             instead of 2-sigma outliers.  Yields many more matches
 *             at end of iterating, and works better for TASS Mark IV
 *             data.  May 25, 2000
 *           Michael Richmond
 *
 *           Added new public function "atRecalcTrans", which uses
 *             existing matched lists of stars to calculate a TRANS.
 *             Allows us to use _all_ the stars in both lists,
 *             not just the brightest N.  June 1, 2000
 *           Michael Richmond
 *
 *           Added 'recalc_flag' argument to 'iter_trans()', to allow
 *             different behavior on first iteration if we call it
 *             from atRecalcTrans or not.  June 2, 2000
 *           Michael Richmond
 *
 *           Changed the "3" in "3-sigma" outliers rejected in iter_trans
 *             into a #define value AT_MATCH_NSIGMA, defined in the
 *             .h file.  June 2, 2000
 *           Michael Richmond
 *
 *           Modified so that this single file contains routines to handle
 *             the linear, quadratic, and cubic transformation cases.
 *             June 10, 2000
 *           Michael Richmond
 *
 *           Replaced old 'gauss_jordon' routine to solve matrix equation
 *             with new 'gauss_matrix' routine; uses Gaussian elimination
 *             with back-substitution instead of Gauss-Jordon technique.
 *             Not associated with "Numerical Recipes" in any way -- hah!
 *             June 19, 2000
 *           Michael Richmond
 *
 *           Added MEDTF calculations and TRANS diagnostics, as suggested
 *             by John Blakeslee.
 *             Dec 12, 2001
 *           Michael Richmond
 *
 *           Fixed off-by-one bug in "remove_elem", as suggested by
 *             Andrew Bennett.
 *           Also fixed small error in location of paranthesis in
 *             the "gauss_matrix" routine, again thanks to Andrew.
 *             Dec 28, 2001
 *           Michael Richmond
 *
 *           Fixed bug in "atCalcRMS" which caused assertion to fail
 *             when the routine was given two empty lists.
 *           Similar bugs addressed by making early checks to the
 *             number of stars in the passed lists in other funcs:
 *                atFindMedtf
 *             Nov 22, 2002
 *           Michael Richmond
 *
 *           Added the ratio of triangle sizes to a debugging printf
 *             statement in make_vote_matrix.
 *             June 26, 2010
 *           Michael Richmond
 *
 *           Fixed bug in find_quick_match() which caused premature
 *             termination if one particular match failed.
 *             Use new routine "check_trans_properties()" to see
 *             if the current match is good enough for user's
 *             needs.  This duplicates some code, but simplifies
 *             things from an overall view.
 *             Aug 7, 2012.
 *           Michael Richmond
 */

#include <stdio.h>
#include <math.h>           /* need this for 'sqrt' in calc_distances */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "core/siril.h"
#include "registration/matching/misc.h"
#include "registration/matching/atpmatch.h"

#undef DEBUG           /* get some of diagnostic output */
#undef DEBUG2          /* get LOTS more diagnostic output */
#undef DEBUG3         /* run 'test_routine' to test matrix inversion */

/*
 * used in the 'gauss_matrix' matrix inversion routine,
 * as a check for very very small numbers which might cause
 * the matrix solution to be unstable.
 */
#define MATRIX_TOL     1.0e-12

/*
 * To evaluate the quality of a match between two sets of stars,
 *   we look at the differences in their positions after transforming
 *   those in list A to the coordinate system of list B.  We sort
 *   those distances and pick the one closest to this percentile
 *   to characterize the distribution.   One stdev should include
 *   about 68% of the data.
 * This is used in routine 'iter_trans'.
 */
#define ONE_STDEV_PERCENTILE  0.683

/*
 * these values are used to tell iter_trans() whether it is being
 *   called from atRecalcTrans or not.  If yes, then it is safe
 *   to include _all_ the matched pairs in finding TRANS, in the
 *   very first iteration.  If no, then only use the best AT_MATCH_STARTN
 *   pairs in the first iteration.
 */
#define RECALC_YES       1
#define RECALC_NO        0

/*
 * the following are "private" functions, used internally only.
 */

/* this typedef is used several sorting routines */
typedef int (*PFI)();

/*
 * use this structure to create an auxiliary array of stars
 *    sorted by one of its coordinates, so that we can quickly
 *    locate a star in one of the main arrays when we are
 *    looking for a match
 */
typedef struct s_star_coord {
	int index; /* index of star in main star array */
	double x; /* "X" value of star */
	double y; /* "Y" value of star */
} s_star_coord;

static int set_star(s_star *star, double x, double y, double mag);
static void copy_star(s_star *from_ptr, s_star *to_ptr);
static void copy_star_array(s_star *from_array, s_star *to_array, int num);
#ifdef DEBUG
static void print_star_array(s_star *array, int num);
#endif
static double **calc_distances(s_star *star_array, int numstars);
static void free_distances(double **array, int num);
#ifdef DEBUG
static void print_dist_matrix(double **matrix, int num);
#endif
static void set_triangle(s_triangle *triangle, s_star *star_array, int i, int j,
		int k, double **dist_matrix);
#ifdef DEBUG2
static void print_one_triangle(s_triangle *triangle, s_star *star_array);
static void print_triangle_array(s_triangle *t_array, int numtriangles,
		s_star *star_array, int numstars);
#endif
static s_triangle *stars_to_triangles(s_star *star_array, int numstars,
		int nbright, int *numtriangles);
static void sort_star_by_mag(s_star *array, int num);
static int compare_star_by_mag(s_star *star1, s_star *star2);
static void sort_star_by_x(s_star *array, int num);
static int compare_star_by_x(s_star *star1, s_star *star2);
static void sort_star_by_match_id(s_star *array, int num);
static int compare_star_by_match_id(s_star *star1, s_star *star2);
static void sort_star_coord_by_x(s_star_coord *array, int num);
static int compare_star_coord_by_x(s_star_coord *sc1, s_star_coord *sc2);
static int fill_triangle_array(s_star *star_array, int numstars,
		double **dist_matrix, int numtriangles, s_triangle *t_array);
static void sort_triangle_array(s_triangle *array, int num);
static void sort_triangle_by_yt(s_triangle *array, int num);
static void sort_triangle_by_D(s_triangle *array, int num);
static int compare_triangle(s_triangle *triangle1, s_triangle *triangle2);
static int compare_triangle_by_yt(s_triangle *triangle1, s_triangle *triangle2);
static int compare_triangle_by_D(s_triangle *triangle1, s_triangle *triangle2);
static int find_ba_triangle(s_triangle *array, int num, double ba0);
static int find_yt_triangle(s_triangle *array, int num, double yt0);
static int find_star_coord_by_x(s_star_coord *array, int num, double x0);
static void prune_triangle_array(s_triangle *t_array, int *numtriangles);
static int **make_vote_matrix(s_star *star_array_A, int num_stars_A,
		s_star *star_array_B, int num_stars_B, s_triangle *t_array_A,
		int num_triangles_A, s_triangle *t_array_B, int num_triangles_B,
		int nbright, double radius, double min_scale, double max_scale,
		double rotation_deg, double tolerance_deg);
static int find_quick_match(s_star *star_array_A, int num_stars_A,
		s_star *star_array_B, int num_stars_B, s_triangle *t_array_A,
		int num_triangles_A, s_triangle *t_array_B, int num_triangles_B,
		int nbright, double star_match_radius, double radius, int max_iter,
		double max_sigma, double min_scale, double max_scale,
		double rotation_deg, double tolerance_deg, int min_req_pairs,
		TRANS *trans);
#ifdef DEBUG
static void print_vote_matrix(int **vote_matrix, int numcells);
#endif
static int top_vote_getters(int **vote_matrix, int num, int **winner_votes,
		int **winner_index_A, int **winner_index_B);
static int calc_trans(int nbright, s_star *star_array_A, int num_stars_A,
		s_star *star_array_B, int num_stars_B, int *winner_votes,
		int *winner_index_A, int *winner_index_B, TRANS *trans);
static s_star *list_to_array(int num_stars, struct s_star *list);
static void reset_array_ids(struct s_star *star_list, int num_stars,
		struct s_star *star_array);

/*
 * these are functions used to solve a matrix equation which
 * gives us the transformation from one coord system to the other
 */
static int gauss_matrix(double **matrix, int num, double *vector);
static int gauss_pivot(double **matrix, int num, double *vector,
		double *biggest_val, int row);

static double ** alloc_matrix(int n);
static void free_matrix(double **matrix, int n);
#ifdef DEBUG
static void print_matrix(double **matrix, int n);
#endif

#ifdef DEBUG3
static void test_routine(void);
#endif

static int iter_trans(int nbright, s_star *star_array_A, int num_stars_A,
		s_star *star_array_B, int num_stars_B, int *winner_votes,
		int *winner_index_A, int *winner_index_B, int recalc_flag, int max_iter,
		double halt_sigma, TRANS *trans);
static int compare_double(double *f1, double *f2);
static double find_percentile(double *array, int num, double perc);
static int calc_trans_coords(s_star *star, TRANS *trans, double *newx,
		double *newy);
static int apply_trans(s_star *star_array, int num_stars, TRANS *trans);
static int double_sort_by_match_id(s_star *star_array_A, int num_stars_A,
		s_star *star_array_B, int num_stars_B);
static int match_arrays_slow(s_star *star_array_A, int num_stars_A,
		s_star *star_array_B, int num_stars_B, double radius,
		s_star **star_array_J, int *num_stars_J, s_star **star_array_K,
		int *num_stars_K, s_star **star_array_L, int *num_stars_L,
		s_star **star_array_M, int *num_stars_M);
static int add_element(s_star *new_star, s_star **star_array, int *total_num,
		int *current_num);
static void remove_elem(s_star *star_array, int num, int *num_stars);
static int remove_repeated_elements(s_star *star_array_1, int *num_stars_1,
		s_star *star_array_2, int *num_stars_2);
static void remove_same_elements(s_star *star_array_1, int num_stars_1,
		s_star *star_array_2, int *num_stars_2);
static void write_array(int num_stars, struct s_star *star_array,
		char *filename);
static int is_desired_rotation(struct s_triangle *tri_1,
		struct s_triangle *tri_2, double want_angle_deg, double tolerance_deg,
		double *actual_angle_deg);
static int apply_trans_and_find_matches(s_star *star_array_A, int num_stars_A,
		s_star *star_array_B, int num_stars_B, s_star_coord *sorted_array_B,
		double star_match_radius, TRANS *trans, int *num_winners,
		int *winner_index_A, int *winner_index_B);
static int compute_match_distance_stats(s_star *star_array_A, int num_stars_A,
		s_star *star_array_B, int num_stars_B, int num_matches,
		int *match_index_A, int *match_index_B, double *mean, double *stdev);
static int prune_matched_pairs(s_star *star_array_A, int num_stars_A,
		s_star *star_array_B, int num_stars_B, int num_matches,
		int *match_index_A, int *match_index_B, double critical_distance,
		int *remaining_pairs);

static int eval_trans_quality(s_star *star_array_A, int num_stars_A,
		s_star *star_array_B, int num_stars_B, double critical_distance,
		TRANS *trans);
static int calc_trans_sig(int num_matches, s_star *star_array_A,
		int num_stars_A, s_star *star_array_B, int num_stars_B,
		int *winner_votes, int *winner_index_A, int *winner_index_B,
		TRANS *trans);
static int is_trans_good_enough(int min_matches, double max_stdev, TRANS *trans);

/*
 * we have three different versions of a routine to do the
 * dirty work of calculating the TRANS which best turns
 * coords of system A into coords of system B.
 * There is a version for linear transformations, quadratic ones,
 * and cubic ones.
 *
 * All are called by the 'calc_trans' function; it uses the
 * trans->order value to figure out which one is appropriate.
 */
static int calc_trans_linear(int nbright, s_star *star_array_A, int num_stars_A,
		s_star *star_array_B, int num_stars_B, int *winner_votes,
		int *winner_index_A, int *winner_index_B, TRANS *trans);
static int calc_trans_quadratic(int nbright, s_star *star_array_A,
		int num_stars_A, s_star *star_array_B, int num_stars_B,
		int *winner_votes, int *winner_index_A, int *winner_index_B,
		TRANS *trans);
static int calc_trans_cubic(int nbright, s_star *star_array_A, int num_stars_A,
		s_star *star_array_B, int num_stars_B, int *winner_votes,
		int *winner_index_A, int *winner_index_B, TRANS *trans);
static int
check_trans_properties(TRANS *trans, double min_scale, double max_scale,
		double rotation_deg, double tolerance_deg);

/************************************************************************
 * <AUTO EXTRACT>
 *
 * ROUTINE: atFindTrans
 *
 * DESCRIPTION:
 * This function is based on the algorithm described in Valdes et al.,
 * PASP 107, 1119 (1995).  It tries to
 *         a. match up objects in the two chains
 *         a. find a coordinate transformation that takes coords in
 *               objects in chainA and changes to those in chainB.
 *
 *
 * Actually, this is a top-level "driver" routine that calls smaller
 * functions to perform actual tasks.  It mostly creates the proper
 * inputs and outputs for the smaller routines.
 *
 * RETURN:
 *    SH_SUCCESS         if all goes well
 *    SH_GENERIC_ERROR   if an error occurs
 *
 * </AUTO>
 */

int atFindTrans(int numA, /* I: number of stars in list A */
struct s_star *listA, /* I: match this set of objects with list B */
int numB, /* I: number of stars in list B */
struct s_star *listB, /* I: match this set of objects with list A */
double star_match_radius, /* I: max radius in star-space allowed for */
/*       a pair of stars to match */
double radius, /* I: max radius in triangle-space allowed for */
/*       a pair of triangles to match */
int nobj, /* I: max number of bright stars to use in creating */
/*       triangles for matching from each list */
double min_scale, /* I: minimum permitted relative scale factor */
/*       if -1, any scale factor is allowed */
double max_scale, /* I: maximum permitted relative scale factor */
/*       if -1, any scale factor is allowed */
double rotation_deg, /* I: desired relative angle of coord systems (deg) */
/*       if AT_MATCH_NOANGLE, any orientation is allowed */
double tolerance_deg, /* I: allowed range of orientation angles (deg) */
/*       if AT_MATCH_NOANGLE, any orientation is allowed */
int max_iter, /* I: go through at most this many iterations */
/*       in the iter_trans() loop. */
double max_sigma, /* I: if the mean residual becomes this small */
/*       then the match was a success */
int min_req_pairs, /* I: must have at least this many matched pairs */
/*       of stars be count as successful match */
TRANS *trans /* O: place into this TRANS structure's fields */
/*       the coeffs which convert coords of chainA */
/*       into coords of chainB system. */
) {
	int nbright, min;
	int num_stars_A; /* number of stars in chain A */
	int num_stars_B; /* number of stars in chain B */
	int num_triangles_A; /* number of triangles formed from chain A */
	int num_triangles_B; /* number of triangles formed from chain B */
	int start_pairs = 0;
	s_star *star_array_A = NULL;
	s_star *star_array_B = NULL;
	s_triangle *triangle_array_A = NULL;
	s_triangle *triangle_array_B = NULL;

#ifdef DEBUG
	printf(" entering atFindTrans \n");
#endif

	num_stars_A = numA;
	num_stars_B = numB;
	star_array_A = list_to_array(numA, listA);
	star_array_B = list_to_array(numB, listB);

#ifdef DEBUG3
	test_routine();
#endif

	shAssert(star_array_A != NULL);
	shAssert(star_array_B != NULL);

	switch (trans->order) {
	case AT_TRANS_LINEAR:
		start_pairs = AT_MATCH_STARTN_LINEAR;
		break;
	case AT_TRANS_QUADRATIC:
		start_pairs = AT_MATCH_STARTN_QUADRATIC;
		break;
	case AT_TRANS_CUBIC:
		start_pairs = AT_MATCH_STARTN_CUBIC;
		break;
	default:
		shError("atFindTrans: invalid trans->order %d ", trans->order);
		break;
	}

	/*
	 * here we check to see if each list of stars contains a
	 * required minimum number of stars.  If not, we return with
	 * an error message, and SH_GENERIC_ERROR.
	 *
	 * In addition, we check to see that each list has at least 'nobj'
	 * items.  If not, we set 'nbright' to the minimum of the two
	 * list lengths, and print a warning message so the user knows
	 * that we're using fewer stars than he asked.
	 *
	 * On the other hand, if the user specifies a value of "nobj" which
	 * is too SMALL, then we ignore it and use the smallest valid
	 * value (which is start_pairs).
	 */
	min = (num_stars_A < num_stars_B ? num_stars_A : num_stars_B);
	if (min < start_pairs) {
		shError("atFindTrans: only %d stars in list(s), require at least %d",
				min, start_pairs);
		free_star_array(star_array_A);
		free_star_array(star_array_B);
		return (SH_GENERIC_ERROR);
	}
	if (nobj > min) {
		shDebug(AT_MATCH_ERRLEVEL,
				"atFindTrans: using only %d stars, fewer than requested %d",
				min, nobj);
		nbright = min;
	} else {
		nbright = nobj;
	}
	if (nbright < start_pairs) {
		shDebug(AT_MATCH_ERRLEVEL,
				"atFindTrans: must use %d stars, more than requested %d",
				start_pairs, nobj);
		nbright = start_pairs;
	}

	/* this is a sanity check on the above checks */
	shAssert((nbright >= start_pairs) && (nbright <= min));

#ifdef DEBUG
	printf("here comes star array A\n");
	print_star_array(star_array_A, num_stars_A);
	printf("here comes star array B\n");
	print_star_array(star_array_B, num_stars_B);
#endif

	/*
	 * we now convert each list of stars into a list of triangles,
	 * using only a subset of the "nbright" brightest items in each list.
	 */
	triangle_array_A = stars_to_triangles(star_array_A, num_stars_A, nbright,
			&num_triangles_A);
	shAssert(triangle_array_A != NULL);
	triangle_array_B = stars_to_triangles(star_array_B, num_stars_B, nbright,
			&num_triangles_B);
	shAssert(triangle_array_B != NULL);

	/*
	 * sort all triangles in list A by their D value
	 */
	sort_triangle_by_D(triangle_array_A, num_triangles_A);
#ifdef DEBUG2
	printf("after sorting by D, here comes triangle array A\n");
	print_triangle_array(triangle_array_A, num_triangles_A, star_array_A,
			num_stars_A);
#endif

	/*
	 * sort all triangles in list B by their yt values
	 */
	sort_triangle_by_yt(triangle_array_B, num_triangles_B);
#ifdef DEBUG2
	printf("after sorting by yt, here comes triangle array B\n");
	print_triangle_array(triangle_array_B, num_triangles_B, star_array_B,
			num_stars_B);
#endif

	/*
	 * The important step: walk through list A, checking
	 * as we go for matches in list B.  If we find a
	 * possible match, evaluate it using a large number
	 * of objects; if it's a good match, terminate the
	 * search.
	 */
	if (find_quick_match(star_array_A, num_stars_A, star_array_B, num_stars_B,
			triangle_array_A, num_triangles_A, triangle_array_B,
			num_triangles_B, nbright, star_match_radius, radius, max_iter,
			max_sigma, min_scale, max_scale, rotation_deg, tolerance_deg,
			min_req_pairs, trans) == SH_SUCCESS) {
		/* we found a match */
#ifdef DEBUG
		printf("find_quick_match returns with success! \n");
#endif
	} else {
#ifdef DEBUG
		printf("find_quick_match returns with failure \n");
#endif
		shError("atFindTrans: find_quick_match unable to create a valid TRANS");
		free_star_array(star_array_A);
		free_star_array(star_array_B);
		return (SH_GENERIC_ERROR);
	}

	/*
	 * clean up memory we allocated during the matching process
	 */
	free_star_array(star_array_A);
	free_star_array(star_array_B);

	return (SH_SUCCESS);
}

/************************************************************************
 * <AUTO EXTRACT>
 *
 * ROUTINE: atRecalcTrans
 *
 * DESCRIPTION:
 * Given two lists of stars which ALREADY have been matched,
 * this routine finds a coord transformation which takes coords
 * of stars in list A to those in list B.
 *
 * We can skip all the matching-triangles business, which makes this
 * _much_ faster than atFindTrans.
 *
 * RETURN:
 *    SH_SUCCESS         if all goes well
 *    SH_GENERIC_ERROR   if an error occurs
 *
 * </AUTO>
 */

int atRecalcTrans(int numA, /* I: number of stars in list A */
struct s_star *listA, /* I: match this set of objects with list B */
int numB, /* I: number of stars in list B */
struct s_star *listB, /* I: match this set of objects with list A */
int max_iter, /* I: go through at most this many iterations */
/*       in the iter_trans() loop. */
double halt_sigma, /* I: halt the fitting procedure if the mean */
/*       residual becomes this small */
TRANS *trans /* O: place into this TRANS structure's fields */
/*       the coeffs which convert coords of chainA */
/*       into coords of chainB system. */
) {
	int i, nbright, min;
	int num_stars_A; /* number of stars in chain A */
	int num_stars_B; /* number of stars in chain B */
	int *winner_votes; /* # votes gotten by top pairs of matched stars */
	int *winner_index_A; /* elem i in this array is index in star array A */
	/*    which matches ... */
	int *winner_index_B; /* elem i in this array, index in star array B */
	int start_pairs = 0;
	s_star *star_array_A = NULL;
	s_star *star_array_B = NULL;

	num_stars_A = numA;
	num_stars_B = numB;
	star_array_A = list_to_array(numA, listA);
	star_array_B = list_to_array(numB, listB);

	shAssert(star_array_A != NULL);
	shAssert(star_array_B != NULL);

	switch (trans->order) {
	case AT_TRANS_LINEAR:
		start_pairs = AT_MATCH_STARTN_LINEAR;
		break;
	case AT_TRANS_QUADRATIC:
		start_pairs = AT_MATCH_STARTN_QUADRATIC;
		break;
	case AT_TRANS_CUBIC:
		start_pairs = AT_MATCH_STARTN_CUBIC;
		break;
	default:
		shError("atRecalcTrans: invalid trans->order %d ", trans->order);
		break;
	}

	/*
	 * here we check to see if each list of stars contains a
	 * required minimum number of stars.  If not, we return with
	 * an error message, and SH_GENERIC_ERROR.
	 *
	 * We set 'nbright' to the minimum of the two list lengths
	 */
	min = (num_stars_A < num_stars_B ? num_stars_A : num_stars_B);
	if (min < start_pairs) {
		shError("atRecalcTrans: only %d stars in list(s), require at least %d",
				min, start_pairs);
		free_star_array(star_array_A);
		free_star_array(star_array_B);
		return (SH_GENERIC_ERROR);
	}
	nbright = min;

	/* this is a sanity check on the above checks */
	shAssert((nbright >= start_pairs) && (nbright <= min));

#ifdef DEBUG
	printf("here comes star array A\n");
	print_star_array(star_array_A, num_stars_A);
	printf("here comes star array B\n");
	print_star_array(star_array_B, num_stars_B);
#endif

	/*
	 * We need to create dummy arrays for 'winner_votes', and the
	 * 'winner_index' arrays.  We already know that all these stars
	 * are good matches, and so we can just create some arrays
	 * and fill them with identical numbers.  They aren't used by
	 * iter_trans(), anyway.
	 */
	winner_votes = (int *) shMalloc(nbright * sizeof(int));
	winner_index_A = (int *) shMalloc(nbright * sizeof(int));
	winner_index_B = (int *) shMalloc(nbright * sizeof(int));
	for (i = 0; i < nbright; i++) {
		winner_votes[i] = 100;
		winner_index_A[i] = i;
		winner_index_B[i] = i;
	}

	/*
	 * next, we take ALL the matched pairs of coodinates, and
	 * figure out a transformation of the form
	 *
	 *       x' = A + Bx + Cx
	 *       y' = D + Ex + Fy
	 *
	 * (i.e. a TRANS structure) which converts the coordinates
	 * of objects in list A to those in list B
	 */
	if (iter_trans(nbright, star_array_A, num_stars_A, star_array_B,
			num_stars_B, winner_votes, winner_index_A, winner_index_B,
			RECALC_YES, max_iter, halt_sigma, trans) != SH_SUCCESS) {

		shError("atRecalcTrans: iter_trans unable to create a valid TRANS");
		free_star_array(star_array_A);
		free_star_array(star_array_B);
		return (SH_GENERIC_ERROR);
	}

#ifdef DEBUG
	printf("  after calculating new TRANS structure, here it is\n");
	print_trans(trans);
#endif

	/*
	 * clean up memory we allocated during the matching process
	 */
	shFree(winner_votes);
	shFree(winner_index_A);
	shFree(winner_index_B);
	free_star_array(star_array_A);
	free_star_array(star_array_B);

	return (SH_SUCCESS);
}

/************************************************************************
 * <AUTO EXTRACT>
 *
 * ROUTINE: atApplyTrans
 *
 * DESCRIPTION:
 * Given a list of s_star structures, apply the given TRANS to each item in
 * the list, modifying the "x" and "y" values.
 *
 * The TRANS structure has 6 coefficients, which are used as follows:
 *
 *       x' = A + Bx + Cx
 *       y' = D + Ex + Fy
 *
 * Actually, this is a top-level "driver" routine that calls smaller
 * functions to perform actual tasks.  It mostly creates the proper
 * inputs and outputs for the smaller routines.
 *
 *
 * RETURN:
 *    SH_SUCCESS         if all goes well
 *    SH_GENERIC_ERROR   if an error occurs
 *
 * </AUTO>
 */

int atApplyTrans(int num, /* I: number of stars in the linked list */
s_star *star_list, /* I/O: modify x,y coords of objects in this list */
TRANS *trans /* I: use this TRANS to transform the coords of */
/*       items in the list */
) {
	int i;
	struct s_star *ptr, *star_array;

	shAssert(star_list != NULL);

	/*
	 * convert the linked list to an array
	 */
	star_array = list_to_array(num, star_list);

#ifdef DEBUG
	printf("before applying TRANS \n");
	print_star_array(star_array, num);
#endif

	/*
	 * next, apply the transformation to each element of the array
	 */
	apply_trans(star_array, num, trans);

#ifdef DEBUG
	printf("after applying TRANS \n");
	print_star_array(star_array, num);
#endif

	/*
	 * transfer the coord values from the array back into the list
	 */
	for (ptr = star_list, i = 0; i < num; i++, ptr = ptr->next) {
		shAssert(ptr != NULL);
		ptr->x = star_array[i].x;
		ptr->y = star_array[i].y;
	}

	/*
	 * delete the array
	 */
	free_star_array(star_array);

	/*
	 * all done!
	 */

	return (SH_SUCCESS);
}

/************************************************************************
 * <AUTO EXTRACT>
 *
 * ROUTINE: atMatchLists
 *
 * DESCRIPTION:
 * Given 2 lists of s_star structures,
 * which have ALREADY been transformed so that the "x"
 * and "y" coordinates of each list are close to each other
 * (i.e. matching items from each list have very similar "x" and "y")
 * this routine attempts to find all instances of matching items
 * from the 2 lists.
 *
 * We consider a "match" to be the closest coincidence of centers
 * which are within "radius" pixels of each other.
 *
 * Use a slow, but sure, algorithm.
 *
 * We will match objects from A --> B.  It is possible to have several
 * As that match to the same B:
 *
 *           A1 -> B5   and A2 -> B5
 *
 * This function finds such multiple-match items and deletes all but
 * the closest of the matches.
 *
 * place the elems of A that are matches into output list J
 *                    B that are matches into output list K
 *                    A that are not matches into output list L
 *                    B that are not matches into output list M
 *
 * Place a count of the number of matching pairs into the final
 * argument, 'num_matches'.
 *
 *
 * RETURN:
 *    SH_SUCCESS         if all goes well
 *    SH_GENERIC_ERROR   if an error occurs
 *
 * </AUTO>
 */

int atMatchLists(int numA, /* I: number of stars in list A */
s_star *listA, /* I: first input list of items to be matched */
int numB, /* I: number of stars in list B */
s_star *listB, /* I: second list of items to be matched */
double radius, /* I: maximum radius for items to be a match */
char *basename, /* I: base of filenames used to store the */
/*      output; extension indicates which */
/*      .mtA    items from A that matched */
/*      .mtB    items from B that matched */
/*      .unA    items from A that didn't match */
/*      .unB    items from A that didn't match */
int *num_matches /* O: number of matching pairs we find */
) {
	s_star *star_array_A;
	int num_stars_A;
	s_star *star_array_B;
	int num_stars_B;
	s_star *star_array_J, *star_array_K, *star_array_L, *star_array_M;
	int num_stars_J, num_stars_K, num_stars_L, num_stars_M;
	char filename[100];

	shAssert(listA != NULL);
	shAssert(listB != NULL);

	/* convert from linked lists to arrays */
	num_stars_A = numA;
	num_stars_B = numB;
	star_array_A = list_to_array(numA, listA);
	star_array_B = list_to_array(numB, listB);
	shAssert(star_array_A != NULL);
	shAssert(star_array_B != NULL);

	/* reset the 'id' fields in the arrays to match those in the lists */
	reset_array_ids(listA, numA, star_array_A);
	reset_array_ids(listB, numB, star_array_B);

	/* do the matching process */
	if (match_arrays_slow(star_array_A, num_stars_A, star_array_B, num_stars_B,
			radius, &star_array_J, &num_stars_J, &star_array_K, &num_stars_K,
			&star_array_L, &num_stars_L, &star_array_M, &num_stars_M)
			!= SH_SUCCESS) {
		shError("atMatchLists: match_arrays_slow fails");
		return (SH_GENERIC_ERROR);
	}

	/*
	 * Set the 'num_matches' value to the number of matching pairs
	 *   (we could just as well use num_stars_K)
	 */
	*num_matches = num_stars_J;

	/*
	 * now write the output into ASCII text files, each of which starts
	 * with 'basename', but has a different extension.
	 *
	 *    basename.mtA    stars from list A that did match         array J
	 *    basename.mtB    stars from list A that did match         array K
	 *    basename.unA    stars from list A that did NOT match     array L
	 *    basename.unB    stars from list A that did NOT match     array M
	 */
	sprintf(filename, "%s.mtA", basename);
	write_array(num_stars_J, star_array_J, filename);
	sprintf(filename, "%s.mtB", basename);
	write_array(num_stars_K, star_array_K, filename);
	sprintf(filename, "%s.unA", basename);
	write_array(num_stars_L, star_array_L, filename);
	sprintf(filename, "%s.unB", basename);
	write_array(num_stars_M, star_array_M, filename);

	/*
	 * all done!
	 */
	free_star_array(star_array_J);
	free_star_array(star_array_K);
	free_star_array(star_array_L);
	free_star_array(star_array_M);

	return (SH_SUCCESS);
}

/************************************************************************
 * <AUTO EXTRACT>
 *
 * ROUTINE: atBuildSmallFile
 *
 * DESCRIPTION:
 * This function is basically the first half of "atFindTrans".  It takes
 * a single list of s_star structures, performs some checks, and calculates
 * the array of triangles for the bright stars in the list.
 *
 * At that point, it creates a new file with name "outfile", and writes
 * into that file the subset of bright stars and their triangles.
 * We keep that small file for future use.
 *
 * RETURN:
 *    SH_SUCCESS         if all goes well
 *    SH_GENERIC_ERROR   if an error occurs
 *
 * </AUTO>
 */

#if 0
int atBuildSmallFile(double ra, /* I: Right Ascension of field center, in degrees */
double dec, /* I: Declination of field center, in degrees */
int numA, /* I: number of stars in list A */
struct s_star *listA, /* I: create an array of triangles for these stars */
int nobj, /* I: max number of bright stars to use in creating */
/*       triangles for matching from each list */
char *outfile /* I: create a file with this name, and place lists */
/*       of stars and triangles into it */
) {
	int nbright, min;
	int num_stars_A; /* number of stars in chain A */
	int num_triangles_A; /* number of triangles formed from chain A */
	int start_pairs;
	s_star *star_array_A = NULL;
	s_triangle *triangle_array_A = NULL;

	num_stars_A = numA;
	star_array_A = list_to_array(numA, listA);

	shAssert(star_array_A != NULL);

	start_pairs = AT_MATCH_STARTN_LINEAR;

	/*
	 * here we check to see if each list of stars contains a
	 * required minimum number of stars.  If not, we return with
	 * an error message, and SH_GENERIC_ERROR.
	 *
	 * In addition, we check to see that each list has at least 'nobj'
	 * items.  If not, we set 'nbright' to the minimum of the two
	 * list lengths, and print a warning message so the user knows
	 * that we're using fewer stars than he asked.
	 *
	 * On the other hand, if the user specifies a value of "nobj" which
	 * is too SMALL, then we ignore it and use the smallest valid
	 * value (which is start_pairs).
	 */
	min = num_stars_A;
	if (min < start_pairs) {
		shError("atBuildSmallFile: only %d stars in list, require at least %d",
				min, start_pairs);
		free_star_array(star_array_A);
		return (SH_GENERIC_ERROR);
	}
	if (nobj > min) {
		shDebug(AT_MATCH_ERRLEVEL,
				"atBuildSmallFile: using only %d stars, fewer than requested %d",
				min, nobj);
		nbright = min;
	} else {
		nbright = nobj;
	}
	if (nbright < start_pairs) {
		shDebug(AT_MATCH_ERRLEVEL,
				"atBuildSmallFile: must use %d stars, more than requested %d",
				start_pairs, nobj);
		nbright = start_pairs;
	}

	/* this is a sanity check on the above checks */
	shAssert((nbright >= start_pairs) && (nbright <= min));

#ifdef DEBUG
	printf("here comes star array A\n");
	print_star_array(star_array_A, num_stars_A);
#endif

	/*
	 * we now convert the list of stars into a list of triangles,
	 * using only a subset of the "nbright" brightest items in the list.
	 */
	triangle_array_A = stars_to_triangles(star_array_A, num_stars_A, nbright,
			&num_triangles_A);
	shAssert(triangle_array_A != NULL);

	/*
	 * Now we prune the triangle array to eliminate those with
	 * ratios (b/a) > AT_MATCH_RATIO,
	 * since Valdes et al. say that this speeds things up and eliminates
	 * lots of closely-packed triangles.
	 */
	prune_triangle_array(triangle_array_A, &num_triangles_A);
#ifdef DEBUG2
	printf("after pruning, here comes triangle array A\n");
	print_triangle_array(triangle_array_A, num_triangles_A, star_array_A,
			num_stars_A);
#endif

	if (write_small_arrays(ra, dec, num_stars_A, star_array_A, nbright,
			num_triangles_A, triangle_array_A, outfile) != SH_SUCCESS) {
		free_star_array(star_array_A);
		shFree(triangle_array_A);
		shError("atBuildSmallFile: write_small_arrays returns with error");
		return (SH_GENERIC_ERROR);
	}

	/* clean up memory */
	free_star_array(star_array_A);
	shFree(triangle_array_A);

	return (SH_SUCCESS);
}
#endif

/************************************************************************
 * <AUTO EXTRACT>
 *
 * ROUTINE: atSmallTrans
 *
 * DESCRIPTION:
 * This function is basically the second half of the "atFindTrans"
 * function.  We pass it a list of detected stars, and a pre-made
 * array of catalog stars and triangles.
 *
 * The first time we call this routine, we convert the _list_ of
 * detected stars into an array, and create an array of triangles
 * for them.  On all subsequent calls, we re-use these arrays.
 *
 *
 * RETURN:
 *    SH_SUCCESS         if all goes well
 *    SH_GENERIC_ERROR   if an error occurs
 *
 * </AUTO>
 */

int atSmallTrans(int numA, /* I: number of stars in list A */
struct s_star *listA, /* I: match this set of objects with list B */
int numB, /* I: number of stars in pre-made array B */
struct s_star *star_array_B,
/* I: match this array of stars with list A */
int num_triangles_B, /* I: number of pre-made triangles from set B */
struct s_triangle *triangle_array_B,
/* I: pre-made array of triangles */
double radius, /* I: max radius in triangle-space allowed for */
/*       a pair of triangles to match */
int nobj, /* I: max num of bright stars to use in creating */
/*       triangles for matching from each list */
double min_scale, /* I: minimum permitted relative scale factor */
/*       if -1, any scale factor is allowed */
double max_scale, /* I: maximum permitted relative scale factor */
/*       if -1, any scale factor is allowed */
double rotation_deg, /* I: desired relative angle of coord systems (deg) */
/*       if AT_MATCH_NOANGLE, any orientation is allowed */
double tolerance_deg, /* I: allowed range of orientation angles (deg) */
/*       if AT_MATCH_NOANGLE, any orientation is allowed */
int max_iter, /* I: go through at most this many iterations */
/*       in the iter_trans() loop. */
double halt_sigma, /* I: halt the fitting procedure if the mean */
/*       residual becomes this small */
TRANS *trans, /* O: place into this TRANS structure's fields */
/*       the coeffs which convert coords of "A" */
/*       into coords of "B" system. */
int *ntop, /* O: number of top "vote getters" in the */
/*       top_votes[] array */
int **top_votes /* O: array of votes gotten by the ntop stars */
/*       will help evaluate the quality of matches */
) {
	int i, nbright, min, ret;
	int num_stars_A; /* number of stars in set A  */
	int num_stars_B; /* number of stars in set B */
	static int num_triangles_A; /* number of triangles formed from set A */
	int **vote_matrix;
	int *winner_votes; /* # votes gotten by top pairs of matched stars */
	int *winner_index_A; /* elem i in this array is index in star array A */
	/*    which matches ... */
	int *winner_index_B; /* elem i in this array, index in star array B */
	int start_pairs = 0;
	static s_star *star_array_A = NULL;
	static s_triangle *triangle_array_A = NULL;
	static int first = 1;
	int first_flag = 0;

	/*
	 * the first time we call this routine, create an array of stars from
	 * the items in listA
	 */
	if (first == 1) {
		first = 0;
		first_flag = 1;
		star_array_A = list_to_array(numA, listA);
	}

	num_stars_A = numA;
	num_stars_B = numB;

	shAssert(star_array_A != NULL);
	shAssert(star_array_B != NULL);

	switch (trans->order) {
	case AT_TRANS_LINEAR:
		start_pairs = AT_MATCH_STARTN_LINEAR;
		break;
	case AT_TRANS_QUADRATIC:
		start_pairs = AT_MATCH_STARTN_QUADRATIC;
		break;
	case AT_TRANS_CUBIC:
		start_pairs = AT_MATCH_STARTN_CUBIC;
		break;
	default:
		shError("atFindTrans: invalid trans->order %d ", trans->order);
		break;
	}

	/*
	 * here we check to see if each list of stars contains a
	 * required minimum number of stars.  If not, we return with
	 * an error message, and SH_GENERIC_ERROR.
	 *
	 * In addition, we check to see that each list has at least 'nobj'
	 * items.  If not, we set 'nbright' to the minimum of the two
	 * list lengths, and print a warning message so the user knows
	 * that we're using fewer stars than he asked.
	 *
	 * On the other hand, if the user specifies a value of "nobj" which
	 * is too SMALL, then we ignore it and use the smallest valid
	 * value (which is start_pairs).
	 */
	min = (num_stars_A < num_stars_B ? num_stars_A : num_stars_B);
	if (min < start_pairs) {
		shError("atSmallTrans: only %d stars in list(s), require at least %d",
				min, start_pairs);
		return (SH_GENERIC_ERROR);
	}
	if (nobj > min) {
		shDebug(AT_MATCH_ERRLEVEL,
				"atSmallTrans: using only %d stars, fewer than requested %d",
				min, nobj);
		nbright = min;
	} else {
		nbright = nobj;
	}
	if (nbright < start_pairs) {
		shDebug(AT_MATCH_ERRLEVEL,
				"atSmallTrans: must use %d stars, more than requested %d",
				start_pairs, nobj);
		nbright = start_pairs;
	}

	/* this is a sanity check on the above checks */
	shAssert((nbright >= start_pairs) && (nbright <= min));

#ifdef DEBUG
	printf("here comes star array A\n");
	print_star_array(star_array_A, num_stars_A);
	printf("here comes star array B\n");
	print_star_array(star_array_B, num_stars_B);
	fflush(stdout);
#endif

	/*
	 * If this is the first time we've entered this function,
	 * we now convert list A of stars into a list of triangles,
	 * using only a subset of the "nbright" brightest items in the list.
	 */
	if (first_flag == 1) {
		triangle_array_A = stars_to_triangles(star_array_A, num_stars_A,
				nbright, &num_triangles_A);
	}
	shAssert(triangle_array_A != NULL);
	shAssert(triangle_array_B != NULL);

	/*
	 * If this is the first time through the function,
	 * we now prune the "A" triangle arrays to eliminate those with
	 * ratios (b/a) > AT_MATCH_RATIO,
	 * since Valdes et al. say that this speeds things up and eliminates
	 * lots of closely-packed triangles.
	 *
	 * Triangle array "B" should have been pruned when it was created,
	 * so we don't have to do it again.
	 */
	if (first_flag == 1) {
		prune_triangle_array(triangle_array_A, &num_triangles_A);
	}
#ifdef DEBUG2
	printf("after pruning, here comes triangle array A\n");
	print_triangle_array(triangle_array_A, num_triangles_A, star_array_A,
			num_stars_A);
	printf("after pruning, here comes triangle array B\n");
	print_triangle_array(triangle_array_B, num_triangles_B, star_array_B,
			num_stars_B);
	fflush(stdout);
#endif

	/*
	 * Next, we want to try to match triangles in the two arrays.
	 * What we do is to create a "vote matrix", which is a 2-D array
	 * with "nbright"-by-"nbright" cells.  The cell with
	 * coords [i][j] holds the number of matched triangles in which
	 *
	 *        item [i] in star_array_A matches item [j] in star_array_B
	 *
	 * We'll use this "vote_matrix" to figure out a first guess
	 * at the transformation between coord systems.
	 *
	 * Note that if there are fewer than "nbright" stars
	 * in either list, we'll still make the vote_matrix
	 * contain "nbright"-by-"nbright" cells ...
	 * there will just be a lot of cells filled with zero.
	 */
	vote_matrix = make_vote_matrix(star_array_A, num_stars_A, star_array_B,
			num_stars_B, triangle_array_A, num_triangles_A, triangle_array_B,
			num_triangles_B, nbright, radius, min_scale, max_scale,
			rotation_deg, tolerance_deg);

	/*
	 * having made the vote_matrix, we next need to pick the
	 * top 'nbright' vote-getters.  We call 'top_vote_getters'
	 * and are given, in its output arguments, pointers to three
	 * arrays, each of which has 'nbright' elements pertaining
	 * to a matched pair of STARS:
	 *
	 *       winner_votes[]    number of votes of winners, in descending order
	 *       winner_index_A[]  index of star in star_array_A
	 *       winner_index_B[]  index of star in star_array_B
	 *
	 * Thus, the pair of stars which matched in the largest number
	 * of triangles will be
	 *
	 *       star_array_A[winner_index_A[0]]    from array A
	 *       star_array_B[winner_index_A[0]]    from array B
	 *
	 * and the pair of stars which matched in the second-largest number
	 * of triangles will be
	 *
	 *       star_array_A[winner_index_A[1]]    from array A
	 *       star_array_B[winner_index_A[1]]    from array B
	 *
	 * and so on.
	 */
	top_vote_getters(vote_matrix, nbright, &winner_votes, &winner_index_A,
			&winner_index_B);

	/*
	 * here, we disqualify any of the top vote-getters which have
	 * fewer than AT_MATCH_MINVOTES votes.  This may decrease the
	 * number of valid matched pairs below 'nbright', so we
	 * re-set nbright if necessary.
	 */
	for (i = 0; i < nbright; i++) {
		if (winner_votes[i] < AT_MATCH_MINVOTES) {
#ifdef DEBUG
			printf(
					"disqualifying all winners after number %d, nbright now %d\n",
					i, i);
#endif
			nbright = i;
			break;
		}
	}

	/*
	 * next, we take the "top" matched pairs of coodinates, and
	 * figure out a transformation of the form
	 *
	 *       x' = A + Bx + Cx
	 *       y' = D + Ex + Fy
	 *
	 * (i.e. a TRANS structure) which converts the coordinates
	 * of objects in chainA to those in chainB.
	 */
	ret = iter_trans(nbright, star_array_A, num_stars_A, star_array_B,
			num_stars_B, winner_votes, winner_index_A, winner_index_B,
			RECALC_NO, max_iter, halt_sigma, trans);
	if (ret != SH_SUCCESS) {
		shDebug(AT_MATCH_ERRLEVEL,
				"atSmallTrans: iter_trans unable to create a valid TRANS");
		return (SH_GENERIC_ERROR);
	}

#ifdef DEBUG
	printf("  after calculating new TRANS structure, here it is\n");
	print_trans(trans);
#endif

	/*
	 * set the output args "ntop" and "winner_votes"
	 */
	*ntop = nbright;
	*top_votes = winner_votes;

	return (SH_SUCCESS);
}

/*                    end of PUBLIC information                          */
/*-----------------------------------------------------------------------*/
/*                  start of PRIVATE information                         */

/*
 * the functions listed from here on are intended to be used only
 * "internally", called by the PUBLIC functions above.  Users
 * should be discouraged from accessing them directly.
 */

/************************************************************************
 *
 *
 * ROUTINE: set_star
 *
 * DESCRIPTION:
 * Given a pointer to an EXISTING s_star, initialize its values
 * and set x, y, and mag to the given values.
 *
 * RETURN:
 *    SH_SUCCESS        if all goes well
 *    SH_GENERIC_ERROR  if not
 *
 * </AUTO>
 */

static int set_star(s_star *star, /* I: pointer to existing s_star structure */
double x, /* I: star's "X" coordinate */
double y, /* I: star's "Y" coordinate */
double mag /* I: star's "mag" coordinate */
) {
	static int id_number = 0;

	if (star == NULL) {
		shError("set_star: given a NULL star");
		return (SH_GENERIC_ERROR);
	}
	star->id = id_number++;
	star->index = -1;
	star->x = x;
	star->y = y;
	star->mag = mag;
	star->match_id = -1;
	star->next = (s_star *) NULL;

	return (SH_SUCCESS);
}

/************************************************************************
 *
 *
 * ROUTINE: copy_star
 *
 * DESCRIPTION:
 * Copy the contents of the "s_star" to which "from_ptr" points
 * to the "s_star" to which "to_ptr" points.
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

static void copy_star(s_star *from_ptr, /* I: copy contents of _this_ star ... */
s_star *to_ptr /* O: into _this_ star */
) {
	shAssert(from_ptr != NULL);
	shAssert(to_ptr != NULL);

	to_ptr->id = from_ptr->id;
	to_ptr->index = from_ptr->index;
	to_ptr->x = from_ptr->x;
	to_ptr->y = from_ptr->y;
	to_ptr->mag = from_ptr->mag;
	to_ptr->match_id = from_ptr->match_id;
	to_ptr->next = from_ptr->next;

}

/************************************************************************
 *
 *
 * ROUTINE: copy_star_array
 *
 * DESCRIPTION:
 * Given to arrays of "s_star" structures, EACH OF WHICH MUST
 * ALREADY HAVE BEEN ALLOCATED and have "num" elements,
 * copy the contents of the items in "from_array"
 * to those in "to_array".
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

static void copy_star_array(s_star *from_array, /* I: copy contents of _this_ array ... */
s_star *to_array, /* O: into _this_ array */
int num_stars /* I: each aray must have this many elements */
) {
	int i;
	s_star *from_ptr, *to_ptr;

	shAssert(from_array != NULL);
	shAssert(to_array != NULL);

	for (i = 0; i < num_stars; i++) {
		from_ptr = &(from_array[i]);
		to_ptr = &(to_array[i]);
		shAssert(from_ptr != NULL);
		shAssert(to_ptr != NULL);

		to_ptr->id = from_ptr->id;
		to_ptr->index = from_ptr->index;
		to_ptr->x = from_ptr->x;
		to_ptr->y = from_ptr->y;
		to_ptr->mag = from_ptr->mag;
		to_ptr->match_id = from_ptr->match_id;
		to_ptr->next = from_ptr->next;
	}

}

/************************************************************************
 *
 *
 * ROUTINE: free_star_array
 *
 * DESCRIPTION:
 * Delete an array of "num" s_star structures.
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

void free_star_array(s_star *first /* first star in the array to be deleted */
) {
	shFree(first);
}

/************************************************************************
 *
 *
 * ROUTINE: print_star_array
 *
 * DESCRIPTION:
 * Given an array of "num" s_star structures, print out
 * a bit of information on each in a single line.
 *
 * For debugging purposes.
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

#ifdef DEBUG

static void print_star_array(s_star *array, /* I: first star in array */
int num /* I: number of stars in the array to print */
) {
	int i;
	s_star *star;

	for (i = 0; i < num; i++) {
		star = &(array[i]);
		shAssert(star != NULL);
		printf(" %4d %4d %11.4e %11.4e %6.2f\n", i, star->id, star->x, star->y,
				star->mag);
	}
}

#endif  /* DEBUG */

/************************************************************************
 *
 *
 * ROUTINE: calc_distances
 *
 * DESCRIPTION:
 * Given an array of N='numstars' s_star structures, create a 2-D array
 * called "matrix" with NxN elements and fill it by setting
 *
 *         matrix[i][j] = distance between stars i and j
 *
 * where 'i' and 'j' are the indices of their respective stars in
 * the 1-D array.
 *
 * RETURN:
 *    double **array      pointer to array of pointers to each row of array
 *    NULL                if something goes wrong.
 *
 * </AUTO>
 */

static double **
calc_distances(s_star *star_array, /* I: array of s_stars */
int numstars /* I: with this many elements */
) {
	int i, j;
	double **matrix;
	double dx, dy, dist;

	if (numstars == 0) {
		shError("calc_distances: given an array of zero stars");
		return (NULL);
	}

	/* allocate the array, row-by-row */
	matrix = (double **) shMalloc(numstars * sizeof(double *));
	for (i = 0; i < numstars; i++) {
		matrix[i] = (double *) shMalloc(numstars * sizeof(double));
	}

	/* fill up the array */
	for (i = 0; i < numstars - 1; i++) {
		for (j = i + 1; j < numstars; j++) {
			dx = star_array[i].x - star_array[j].x;
			dy = star_array[i].y - star_array[j].y;
			dist = sqrt(dx * dx + dy * dy);
			matrix[i][j] = (double) dist;
			matrix[j][i] = (double) dist;
		}
	}
	/* for safety's sake, let's fill the diagonal elements with zeros */
	for (i = 0; i < numstars; i++) {
		matrix[i][i] = 0.0;
	}

	/* okay, we're done.  return a pointer to the array */
	return (matrix);
}

/************************************************************************
 *
 *
 * ROUTINE: free_distances
 *
 * DESCRIPTION:
 * Given a 2-D array of "num"-by-"num" double elements, free up
 * each row of the array, then free the array itself.
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

static void free_distances(double **array, /* I: square array we'll free */
int num /* I: number of elems in each row */
) {
	int i;

	for (i = 0; i < num; i++) {
		shFree(array[i]);
	}
	shFree(array);
}

/************************************************************************
 *
 *
 * ROUTINE: print_dist_matrix
 *
 * DESCRIPTION:
 * Given a 2-D array of "num"-by-"num" distances between pairs of
 * stars, print out the 2-D array in a neat fashion.
 *
 * For debugging purposes.
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

#ifdef DEBUG

static void print_dist_matrix(double **matrix, /* I: pointer to start of 2-D square array */
int num /* I: number of rows and columns in the array */
) {
	int i, j;

	for (i = 0; i < num; i++) {
		shAssert(matrix[i] != NULL);
		for (j = 0; j < num; j++) {
			printf("%11.4e ", matrix[i][j]);
		}
		printf("\n");
	}
}

#endif /* DEBUG */

/************************************************************************
 *
 *
 * ROUTINE: set_triangle
 *
 * DESCRIPTION:
 * Set the elements of some given, EXISTING instance of an "s_triangle"
 * structure, given (the indices to) three s_star structures for its vertices.
 * We check to make sure
 * that the three stars are three DIFFERENT stars, asserting
 * if not.
 *
 * The triangle's "a_index" is set to the position of the star opposite
 * its side "a" in its star array, and similarly for "b_index" and "c_index".
 *
 * Also set the Tabur-style fields
 *     xt  =  dot product of two long sides
 *     yt  =  ratio of longest to shortest side  ( = 1 / tri->ca)
 *     D   =  product of xt and yt
 *
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

static void set_triangle(s_triangle *tri, /* we set fields of this existing structure */
s_star *star_array, /* use stars in this array as vertices */
int s1, /* index in 'star_array' of one vertex */
int s2, /* index in 'star_array' of one vertex */
int s3, /* index in 'star_array' of one vertex */
double **darray /* array of distances between stars */
) {
	static int id_number = 0;
	double d12, d23, d13;
	double a = 0.0, b = 0.0, c = 0.0;
	s_star *star1, *star2, *star3;

	shAssert(tri != NULL);
	shAssert((s1 != s2) && (s1 != s3) && (s2 != s3));
	star1 = &star_array[s1];
	star2 = &star_array[s2];
	star3 = &star_array[s3];
	shAssert((star1 != NULL) && (star2 != NULL) && (star3 != NULL));

	tri->id = id_number++;
	tri->index = -1;

	/*
	 * figure out which sides is longest and shortest, and assign
	 *
	 *     "a" to the length of the longest side
	 *     "b"                      intermediate
	 *     "c"                      shortest
	 *
	 * We use temp variables   d12 = distance between stars 1 and 2
	 *                         d23 = distance between stars 2 and 3
	 *                         d13 = distance between stars 1 and 3
	 * for convenience.
	 *
	 */
	d12 = darray[s1][s2];
	d23 = darray[s2][s3];
	d13 = darray[s1][s3];

	/* sanity check */
	shAssert(d12 >= 0.0);
	shAssert(d23 >= 0.0);
	shAssert(d13 >= 0.0);

	if ((d12 >= d23) && (d12 >= d13)) {
		/* this applies if the longest side connects stars 1 and 2 */
		tri->a_index = star3->index;
		a = d12;
		if (d23 >= d13) {
			tri->b_index = star1->index;
			b = d23;
			tri->c_index = star2->index;
			c = d13;
		} else {
			tri->b_index = star2->index;
			b = d13;
			tri->c_index = star1->index;
			c = d23;
		}
	} else if ((d23 > d12) && (d23 >= d13)) {
		/* this applies if the longest side connects stars 2 and 3 */
		tri->a_index = star1->index;
		a = d23;
		if (d12 > d13) {
			tri->b_index = star3->index;
			b = d12;
			tri->c_index = star2->index;
			c = d13;
		} else {
			tri->b_index = star2->index;
			b = d13;
			tri->c_index = star3->index;
			c = d12;
		}
	} else if ((d13 > d12) && (d13 > d23)) {
		/* this applies if the longest side connects stars 1 and 3 */
		tri->a_index = star2->index;
		a = d13;
		if (d12 > d23) {
			tri->b_index = star3->index;
			b = d12;
			tri->c_index = star1->index;
			c = d23;
		} else {
			tri->b_index = star1->index;
			b = d23;
			tri->c_index = star3->index;
			c = d12;
		}
	} else {
		/* we should never get here! */
		shError("set_triangle: impossible situation?!");
		shAssert(0);
	}

	/*
	 * now that we've figured out the longest, etc., sides, we can
	 * fill in the rest of the triangle's elements
	 *
	 * We need to make a special check, in case a == 0.  In that
	 * case, we'll just set the ratios ba and ca = 1.0, and hope
	 * that these triangles are ignored.
	 *
	 * Likewise, if the length of b == 0, then we can't compute
	 * the ratio of side c to side b; in that case, set cb
	 * to 1.0 also.
	 */
	tri->a_length = a;
	if (a > 0.0) {
		tri->ba = b / a;
		tri->ca = c / a;
		if (b > 0.0) {
			tri->cb = c / b;
		} else {
			tri->cb = 1.0;
		}
	} else {
		tri->ba = 1.0;
		tri->ca = 1.0;
		tri->cb = 1.0;
	}
	tri->side_a_angle = atan2(
			star_array[tri->a_index].y - star_array[tri->b_index].y,
			star_array[tri->a_index].x - star_array[tri->b_index].x);
#ifdef DEBUG2
	printf(" triangle %5d  has side_a_angle %lf = %lf deg \n", tri->id,
			tri->side_a_angle, tri->side_a_angle * (180 / 3.14159));
#endif

	/*
	 * Now set the parameters used by Tabur algorithms.
	 *    xt   is   dot product of two longest sides
	 *    yt   is   ratio of longest to shortest side
	 *                  (easy: just the reciprocal of tri->ca)
	 *    D    is   product of xt and yt
	 */
	{
		double xdot, ydot;

		xdot = ((star_array[tri->a_index].x - star_array[tri->c_index].x)
				* (star_array[tri->b_index].x - star_array[tri->c_index].x));
		ydot = ((star_array[tri->a_index].y - star_array[tri->c_index].y)
				* (star_array[tri->b_index].y - star_array[tri->c_index].y));
		tri->xt = xdot + ydot;
		tri->yt = 1.0 / (tri->ca);
		tri->D = tri->xt * tri->yt;
	}
#ifdef DEBUG2
	printf(" triangle %5d  has xt %lf   yt %lf   D %lf  \n", tri->id, tri->xt,
			tri->yt, tri->D);
#endif

	tri->match_id = -1;
	tri->next = (s_triangle *) NULL;
}

/************************************************************************
 *
 *
 * ROUTINE: print_triangle_array
 *
 * DESCRIPTION:
 * Given an array of "numtriangle" s_triangle structures,
 * and an array of "numstars" s_star structures that make them up,
 * print out
 * a bit of information on each triangle in a single line.
 *
 * For debugging purposes.
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

#ifdef DEBUG2

static void print_one_triangle(s_triangle *triangle, /* I: print this triangle */
s_star *star_array /* I: array of stars for this triangle */
) {
	s_star *sa, *sb, *sc;

	shAssert(triangle != NULL);

	sa = &(star_array[triangle->a_index]);
	sb = &(star_array[triangle->b_index]);
	sc = &(star_array[triangle->c_index]);

	printf(
			"%4d %3d (%5.1f,%5.1f) %3d (%5.1f,%5.1f) %3d (%5.1f, %5.1f)  %5.3f %5.3f %5.3f  %8.1f  %6.1f  %7.1e \n",
			triangle->id, triangle->a_index, sa->x, sa->y, triangle->b_index,
			sb->x, sb->y, triangle->c_index, sc->x, sc->y, triangle->ba,
			triangle->ca, triangle->cb, triangle->xt, triangle->yt,
			triangle->D);
}

static void print_triangle_array(s_triangle *t_array, /* I: first triangle in array */
int numtriangles, /* I: number of triangles in the array to print */
s_star *star_array, /* I: array of stars which appear in triangles */
int numstars /* I: number of stars in star_array */
) {
	int i;
	s_triangle *triangle;
	s_star *sa, *sb, *sc;

	for (i = 0; i < numtriangles; i++) {
		triangle = &(t_array[i]);
		shAssert(triangle != NULL);

		printf("%4d ", i);
		print_one_triangle(triangle, star_array);
	}
}

#endif /* DEBUG */

/************************************************************************
 *
 *
 * ROUTINE: stars_to_triangles
 *
 * DESCRIPTION:
 * Convert an array of s_stars to an array of s_triangles.
 * We use only the brightest 'nbright' objects in the linked list.
 * The steps we need to take are:
 *
 *     1. sort the array of s_stars by magnitude, and
 *             set "index" values in the sorted list.
 *     2. calculate star-star distances in the sorted list,
 *             (for the first 'nbright' objects only)
 *             (creates a 2-D array of distances)
 *     3. create a linked list of all possible triangles
 *             (again using the first 'nbright' objects only)
 *     4. clean up -- delete the 2-D array of distances
 *
 * We place the number of triangles in the final argument, and
 * return a pointer to the new array of s_triangle structures.
 *
 * RETURN:
 *    s_triangle *             pointer to new array of triangles
 *                                  (and # of triangles put into output arg)
 *    NULL                     if error occurs
 *
 * </AUTO>
 */

static s_triangle *
stars_to_triangles(s_star *star_array, /* I: array of s_stars */
int numstars, /* I: the total number of stars in the array */
int nbright, /* I: use only the 'nbright' brightest stars */
int *numtriangles /* O: number of triangles we create */
) {
	int numt;
	double **dist_matrix;
	s_triangle *triangle_array;

	/*
	 * check to see if 'nbright' > 'numstars' ... if so, we re-set
	 *          nbright = numstars
	 *
	 * so that we don't have to try to keep track of them separately.
	 * We'll be able to use 'nbright' safely from then on in this function.
	 */
	if (numstars < nbright) {
		nbright = numstars;
	}

	/*
	 * sort the stars in the array by their 'mag' field, so that we get
	 * them in order "brightest-first".
	 */
	sort_star_by_mag(star_array, numstars);

#ifdef DEBUG
	printf("stars_to_triangles: here comes star array after sorting\n");
	print_star_array(star_array, numstars);
#endif

	/*
	 * calculate the distances between each pair of stars, placing them
	 * into the newly-created 2D array called "dist_matrix".  Note that
	 * we only need to include the first 'nbright' stars in the
	 * distance calculations.
	 */
	dist_matrix = calc_distances(star_array, nbright);
	shAssert(dist_matrix != NULL);

#ifdef DEBUG
	printf("stars_to_triangles: here comes distance matrix\n");
	print_dist_matrix(dist_matrix, nbright);
#endif

	/*
	 * create an array of the appropriate number of triangles that
	 * can be formed from the 'nbright' objects.
	 */
	numt = (nbright * (nbright - 1) * (nbright - 2)) / 6;
	*numtriangles = numt;
	triangle_array = (s_triangle *) shMalloc(numt * sizeof(s_triangle));

	/*
	 * now let's fill that array by making all the possible triangles
	 * out of the first 'nbright' objects in the array of stars.
	 */
	fill_triangle_array(star_array, nbright, dist_matrix, *numtriangles,
			triangle_array);

#ifdef DEBUG2
	printf("stars_to_triangles: here comes the triangle array\n");
	print_triangle_array(triangle_array, *numtriangles, star_array, nbright);
#endif

	/*
	 * we've successfully created the array of triangles, so we can
	 * now get rid of the "dist_matrix" array.  We won't need it
	 * any more.
	 */
	free_distances(dist_matrix, nbright);

	return (triangle_array);
}

/************************************************************************
 *
 *
 * ROUTINE: sort_star_by_mag
 *
 * DESCRIPTION:
 * Given an array of "num" s_star structures, sort it in order
 * of increasing magnitude.
 *
 * After sorting, walk through the array and set each star's
 * "index" field equal to the star's position in the array.
 * Thus, the first star will have index=0, and the second index=1,
 * and so forth.
 *
 * Calls the "compare_star_by_mag" function, below.
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

static void sort_star_by_mag(s_star *array, /* I: array of structures to be sorted */
int num /* I: number of stars in the array */
) {
	int i;

	qsort((char *) array, num, sizeof(s_star), (PFI) compare_star_by_mag);

	/* now set the "index" field for each star */
	for (i = 0; i < num; i++) {
		array[i].index = i;
	}
}

/************************************************************************
 *
 *
 * ROUTINE: compare_star_by_mag
 *
 * DESCRIPTION:
 * Given two s_star structures, compare their "mag" values.
 * Used by "sort_star_by_mag".
 *
 * RETURN:
 *    1                  if first star has larger "mag"
 *    0                  if the two have equal "mag"
 *   -1                  if the first has smaller "mag"
 *
 * </AUTO>
 */

static int compare_star_by_mag(s_star *star1, /* I: compare "mag" field of THIS star ... */
s_star *star2 /* I:  ... with THIS star  */
) {
	shAssert((star1 != NULL) && (star2 != NULL));

	if (star1->mag > star2->mag) {
		return (1);
	}
	if (star1->mag < star2->mag) {
		return (-1);
	}
	return (0);
}

/************************************************************************
 *
 *
 * ROUTINE: sort_star_by_x
 *
 * DESCRIPTION:
 * Given an array of "num" s_star structures, sort it in order
 * of increasing "x" values.
 *
 * In this case, we do NOT re-set the "index" field of each
 * s_star after sorting!
 *
 * Calls the "compare_star_by_x" function, below.
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

static void sort_star_by_x(s_star *array, /* I: array of structures to be sorted */
int num /* I: number of stars in the array */
) {
	qsort((char *) array, num, sizeof(s_star), (PFI) compare_star_by_x);
}

/************************************************************************
 *
 *
 * ROUTINE: compare_star_by_x
 *
 * DESCRIPTION:
 * Given two s_star structures, compare their "x" values.
 * Used by "sort_star_by_x".
 *
 * RETURN:
 *    1                  if first star has larger "x"
 *    0                  if the two have equal "x"
 *   -1                  if the first has smaller "x"
 *
 * </AUTO>
 */

static int compare_star_by_x(s_star *star1, /* I: compare "x" field of THIS star ... */
s_star *star2 /* I:  ... with THIS star  */
) {
	shAssert((star1 != NULL) && (star2 != NULL));

	if (star1->x > star2->x) {
		return (1);
	}
	if (star1->x < star2->x) {
		return (-1);
	}
	return (0);
}

/************************************************************************
 *
 *
 * ROUTINE: sort_star_by_match_id
 *
 * DESCRIPTION:
 * Given an array of "num" s_star structures, sort it in order
 * of increasing "match_id" values.
 *
 * In this case, we do NOT re-set the "index" field of each
 * s_star after sorting!
 *
 * Calls the "compare_star_by_match_id" function, below.
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

static void sort_star_by_match_id(s_star *array, /* I: array of structures to be sorted */
int num /* I: number of stars in the array */
) {
	qsort((char *) array, num, sizeof(s_star), (PFI) compare_star_by_match_id);
}

/************************************************************************
 *
 *
 * ROUTINE: compare_star_by_match_id
 *
 * DESCRIPTION:
 * Given two s_star structures, compare their "match_id" values.
 * Used by "sort_star_by_match_id".
 *
 * RETURN:
 *    1                  if first star has larger "match_id"
 *    0                  if the two have equal "match_id"
 *   -1                  if the first has smaller "match_id"
 *
 * </AUTO>
 */

static int compare_star_by_match_id(s_star *star1, /* I: compare "match_id" field of THIS star ... */
s_star *star2 /* I:  ... with THIS star  */
) {
	shAssert((star1 != NULL) && (star2 != NULL));

	if (star1->match_id > star2->match_id) {
		return (1);
	}
	if (star1->match_id < star2->match_id) {
		return (-1);
	}
	return (0);
}

/************************************************************************
 *
 *
 * ROUTINE: fill_triangle_array
 *
 * DESCRIPTION:
 * Given an array of stars, and a matrix of distances between them,
 * form all the triangles possible; place the properties of these
 * triangles into the "t_array" array, which must already have
 * been allocated and contain "numtriangles" elements.
 *
 * RETURN:
 *    SH_SUCCESS           if all goes well
 *    SH_GENERIC_ERROR     if error occurs
 *
 * </AUTO>
 */

static int fill_triangle_array(s_star *star_array, /* I: array of stars we use to form triangles */
int numstars, /* I: use this many stars from the array */
double **dist_matrix, /* I: numstars-by-numstars matrix of distances */
/*       between stars in the star_array */
int numtriangles, /* I: number of triangles in the t_array */
s_triangle *t_array /* O: we'll fill properties of triangles in  */
/*       this array, which must already exist */
) {
	int i, j, k, n;
	s_triangle *triangle;

	shAssert(
			(star_array != NULL) && (dist_matrix != NULL) && (t_array != NULL));

	n = 0;
	for (i = 0; i < numstars - 2; i++) {
		for (j = i + 1; j < numstars - 1; j++) {
			for (k = j + 1; k < numstars; k++) {

				triangle = &(t_array[n]);
				set_triangle(triangle, star_array, i, j, k, dist_matrix);

				n++;
			}
		}
	}
	shAssert(n == numtriangles);

	return (SH_SUCCESS);
}

/************************************************************************
 *
 *
 * ROUTINE: sort_triangle_array
 *
 * DESCRIPTION:
 * Given an array of "num" s_triangle structures, sort it in order
 * of increasing "ba" value (where "ba" is the ratio of lengths
 * of side b to side a).
 *
 * Calls the "compare_triangle" function, below.
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

static void sort_triangle_array(s_triangle *array, /* I: array of structures to be sorted */
int num /* I: number of triangles in the array */
) {
	qsort((char *) array, num, sizeof(s_triangle), (PFI) compare_triangle);
}

/************************************************************************
 *
 *
 * ROUTINE: compare_triangle
 *
 * DESCRIPTION:
 * Given two s_triangle structures, compare their "ba" values.
 * Used by "sort_triangle_array".
 *
 * RETURN:
 *    1                  if first star has larger "ba"
 *    0                  if the two have equal "ba"
 *   -1                  if the first has smaller "ba"
 *
 * </AUTO>
 */

static int compare_triangle(s_triangle *triangle1, /* I: compare "ba" field of THIS triangle ... */
s_triangle *triangle2 /* I:  ... with THIS triangle  */
) {
	shAssert((triangle1 != NULL) && (triangle2 != NULL));

	if (triangle1->ba > triangle2->ba) {
		return (1);
	}
	if (triangle1->ba < triangle2->ba) {
		return (-1);
	}
	return (0);
}

/************************************************************************
 *
 *
 * ROUTINE: find_ba_triangle
 *
 * DESCRIPTION:
 * Given an array of "num" s_triangle structures, which have already
 * been sorted in order of increasing "ba" value, and given one
 * particular "ba" value ba0, return the index of the first triangle
 * in the array which has "ba" >= ba0.
 *
 * We use a binary search, on the "ba" element of each structure.
 *
 * If there is no such triangle, just return the index of the last
 * triangle in the list.
 *
 * Calls the "compare_triangle" function, above.
 *
 * RETURN:
 *    index of closest triangle in array         if all goes well
 *    index of last triangle in array            if nothing close
 *
 * </AUTO>
 */

static int find_ba_triangle(s_triangle *array, /* I: array of structures which been sorted */
int num, /* I: number of triangles in the array */
double ba0 /* I: value of "ba" we seek */
) {
	int top, bottom, mid;

#ifdef DEBUG2
	printf("find_ba_triangle: looking for ba = %.2f\n", ba0);
#endif

	top = 0;
	if ((bottom = num - 1) < 0) {
		bottom = 0;
	}

	while (bottom - top > 2) {
		mid = (top + bottom) / 2;
#ifdef DEBUG2
		printf(" array[%4d] ba=%.2f   array[%4d] ba=%.2f  array[%4d] ba=%.2f\n",
				top, array[top].ba, mid, array[mid].ba, bottom,
				array[bottom].ba);
#endif
		if (array[mid].ba < ba0) {
			top = mid;
		} else {
			bottom = mid;
		}
	}
#ifdef DEBUG2
	printf(" array[%4d] ba=%.2f                       array[%4d] ba=%.2f\n",
			top, array[top].ba, bottom, array[bottom].ba);
#endif

	/*
	 * if we get here, then the item we seek is either "top" or "bottom"
	 * (which may point to the same item in the array).
	 */
	if (array[top].ba < ba0) {
#ifdef DEBUG2
		printf(" returning array[%4d] ba=%.2f \n", bottom, array[bottom].ba);
#endif
		return (bottom);
	} else {
#ifdef DEBUG2
		printf(" returning array[%4d] ba=%.2f \n", top, array[top].ba);
#endif
		return (top);
	}
}

/************************************************************************
 *
 *
 * ROUTINE: prune_triangle_array
 *
 * DESCRIPTION:
 * Given an array of triangles, sort them in increasing order
 * of the side ratio (b/a), and then "ignore" all triangles
 * with (b/a) > AT_MATCH_RATIO.
 *
 * We re-set the arg "numtriangles" as needed, but leave the
 * space in the array allocated (since the array was allocated
 * as a single block).
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

static void prune_triangle_array(s_triangle *t_array, /* I/O: array of triangles to sort and prune  */
int *numtriangles /* I/O: number of triangles in the t_array */
) {
	int i;

	shAssert(t_array != NULL);
	shAssert(numtriangles != NULL);

	/* first, sort the array */
	sort_triangle_array(t_array, *numtriangles);

	/*
	 * now, find the first triangle with "ba" > AT_MATCH_RATIO and
	 * re-set "numtriangles" to be just before it.
	 *
	 * if this would make "numtriangles" < 1, assert
	 */
	for (i = (*numtriangles) - 1; i >= 0; i--) {
		if (t_array[i].ba <= AT_MATCH_RATIO) {
			break;
		}
	}
	*numtriangles = i;
	shAssert(*numtriangles >= 0);
}

/************************************************************************
 *
 *
 * ROUTINE: make_vote_matrix
 *
 * DESCRIPTION:
 * Given two arrays of triangles, and the arrays of stars that make
 * up each set of triangles, try to match up triangles in the two
 * arrays.  Triangles can be considered to match only when the
 * Euclidean distance in "triangle space", created from the two
 * coordinates "ba" and "ca", is within "max_radius".  That is,
 * for two triangles to match, we must satisfy
 *
 *     sqrt[ (t1.ba - t2.ba)^2 + (t1.ca - t2.ca)^2 ] <= max_radius
 *
 * Note that there may be more than one triangle from array A which
 * matches a particular triangle from array B!  That's okay --
 * we treat any 2 which satisfy the above equation as "matched".
 * We rely upon the "vote_array" to weed out false matches.
 *
 * If "min_scale" and "max_scale" are not both -1, then disallow
 * any match for which the
 * ratio of triangles (indicated by "a_length" members)
 * is outside the given values.
 *
 * If "rotation_deg" and "tolerance_deg" are not both AT_MATCH_NOANGLE,
 * then disallow any match for which the two triangles
 * are not oriented at an angle of "rotation_deg" degrees
 * relative to each other (with a tolerance of "tolerance_deg" degrees).
 *
 * For each pair of triangles that matches, increment
 * the "vote" in each "vote cell" for each pair of matching
 * vertices.
 *
 * The "vote matrix" is a 2-D array of 'nbright'-by-'nbright'
 * integers.  We allocate the array in this function, and
 * return a pointer to the array.  Each cell in the array, vote[i][j],
 * contains the number of triangles in which
 *
 *        star_array_A[i] matched star_array_B[j]
 *
 *
 * RETURN:
 *    int **             pointer to new "vote matrix"
 *
 * </AUTO>
 */

static int **
make_vote_matrix(s_star *star_array_A, /* I: first array of stars */
int num_stars_A, /* I: number of stars in star_array_A  */
s_star *star_array_B, /* I: second array of stars */
int num_stars_B, /* I: number of stars in star_array_B  */
s_triangle *t_array_A, /* I: array of triangles from star_array_A */
int num_triangles_A, /* I: number of triangles in t_array_A */
s_triangle *t_array_B, /* I: array of triangles from star_array_B */
int num_triangles_B, /* I: number of triangles in t_array_B */
int nbright, /* I: consider at most this many stars */
/*       from each array; also the size */
/*       of the output "vote_matrix". */
double max_radius, /* I: max radius in triangle-space allowed */
/*       for 2 triangles to be considered */
/*       a matching pair. */
double min_scale, /* I: minimum permitted relative scale factor */
/*       if -1, any scale factor is allowed */
double max_scale, /* I: maximum permitted relative scale factor */
/*       if -1, any scale factor is allowed */
double rotation_deg, /* I: desired relative angle of coord systems (deg) */
/*       if AT_MATCH_NOANGLE, any orientation is allowed */
double tolerance_deg /* I: allowed range of orientation angles (deg) */
/*       if AT_MATCH_NOANGLE, any orientation is allowed */
) {
	int i, j, start_index;
	int **vote_matrix;
	double ba_A, ba_B, ca_A, ca_B, ba_min, ba_max;
	double rad2;
	double ratio;
	double actual_angle_deg;
	struct s_triangle *tri;

	shAssert(star_array_A != NULL);
	shAssert(star_array_B != NULL);
	shAssert(t_array_A != NULL);
	shAssert(t_array_B != NULL);
	shAssert(nbright > 0);
	if (min_scale != -1) {
		shAssert((max_scale != -1) && (min_scale <= max_scale));
	}
	if (max_scale != -1) {
		shAssert((min_scale != -1) && (min_scale <= max_scale));
	}

	/* allocate and initialize the "vote_matrix" */
	vote_matrix = (int **) shMalloc(nbright * sizeof(int *));
	for (i = 0; i < nbright; i++) {
		vote_matrix[i] = (int *) shMalloc(nbright * sizeof(int));
		for (j = 0; j < nbright; j++) {
			vote_matrix[i][j] = 0;
		}
	}

	/*
	 * now, the triangles in "t_array_A" have been sorted by their "ba"
	 * values.  Therefore, we walk through the OTHER array, "t_array_B",
	 * and for each triangle tri_B in it
	 *
	 *      1a. set  ba_min = tri_B.ba - max_radius
	 *      1b. set  ba_max = tri_B.ba + max_radius
	 *
	 * We'll use these values to limit our selection from t_array_A.
	 *
	 *      2. find the first triangle in t_array_A which has
	 *                     ba > ba_min
	 *      3. starting there, step through t_array_A, calculating the
	 *                     Euclidean distance between tri_B and
	 *                     the current triangle from array A.
	 *      4. stop stepping through t_array_A when we each a triangle
	 *                     with ba > ba_max
	 */
	rad2 = max_radius * max_radius;
	for (j = 0; j < num_triangles_B; j++) {

		/*
		 * make sure that this triangle doesn't have a vertex with index
		 * greater than n_bright (because, if it did, we'd overwrite memory
		 * when we tried to increment the vote_matrix array element).
		 *
		 * This is only a problem when called from "smallTrans", with
		 *             num_stars_A > nbright
		 * or
		 *             num_stars_B > nbright
		 */
		tri = &(t_array_B[j]);
		if ((tri->a_index >= nbright) || (tri->b_index >= nbright)
				|| (tri->c_index >= nbright)) {
#ifdef DEBUG2
			printf("make_vote_matrix: skipping B triangle %d\n", j);
#endif
			continue;
		}

#ifdef DEBUG2
		printf("make_vote_matrix: looking for matches to B %d\n", j);
#endif
		ba_B = t_array_B[j].ba;
		ca_B = t_array_B[j].ca;
		ba_min = ba_B - max_radius;
		ba_max = ba_B + max_radius;
#ifdef DEBUG2
		printf("   ba_min = %7.3f  ba_max = %7.3f\n", ba_min, ba_max);
#endif

		start_index = find_ba_triangle(t_array_A, num_triangles_A, ba_min);
		for (i = start_index; i < num_triangles_A; i++) {

			/*
			 * again, skip any triangle which has a vertex with ID > nbright
			 */
			tri = &(t_array_A[i]);
			if ((tri->a_index >= nbright) || (tri->b_index >= nbright)
					|| (tri->c_index >= nbright)) {
#ifdef DEBUG2
				printf("make_vote_matrix: skipping A triangle %d\n", i);
#endif
				continue;
			}

#ifdef DEBUG2
			printf("   looking at A %d\n", i);
#endif
			ba_A = t_array_A[i].ba;
			ca_A = t_array_A[i].ca;

			/* check to see if we can stop looking through A yet */
			if (ba_A > ba_max) {
				break;
			}

			if ((ba_A - ba_B) * (ba_A - ba_B) + (ca_A - ca_B) * (ca_A - ca_B)
					< rad2) {

				/*
				 * check the ratio of lengths of side "a", and discard this
				 * candidate if its outside the allowed range
				 */
				if (min_scale != -1) {
					ratio = t_array_A[i].a_length / t_array_B[j].a_length;
					if (ratio < min_scale || ratio > max_scale) {
						continue;
					}
				}

				/*
				 * check the relative orientations of the triangles.
				 * If they don't match the desired rotation angle,
				 * discard this match.
				 */
				if (rotation_deg != AT_MATCH_NOANGLE) {
					if (is_desired_rotation(&(t_array_A[i]), &(t_array_B[j]),
							rotation_deg, tolerance_deg, &actual_angle_deg)
							== 0) {
						continue;
					}
				}

				/* we have a (possible) match! */
#ifdef DEBUG2
				ratio = t_array_A[i].a_length / t_array_B[j].a_length;
				printf(
						"   match!  A: (%6.3f, %6.3f)   B: (%6.3f, %6.3f)  ratio %9.4e  angle %9.4e\n",
						ba_A, ca_A, ba_B, ca_B, ratio, actual_angle_deg);
#endif
				/*
				 * increment the vote_matrix cell for each matching pair
				 * of stars, one at each vertex
				 */
				vote_matrix[t_array_A[i].a_index][t_array_B[j].a_index]++;
				vote_matrix[t_array_A[i].b_index][t_array_B[j].b_index]++;
				vote_matrix[t_array_A[i].c_index][t_array_B[j].c_index]++;

			}
		}
	}

#ifdef DEBUG
	print_vote_matrix(vote_matrix, nbright);
#endif

	return (vote_matrix);
}

/************************************************************************
 *
 *
 * ROUTINE: print_vote_matrix
 *
 * DESCRIPTION:
 * Print out the "vote_matrix" in a nice format.
 *
 * For debugging purposes.
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

#ifdef DEBUG

static void print_vote_matrix(int **vote_matrix, /* I: the 2-D array we'll print out */
int numcells /* I: number of cells in each row and col of matrix */
) {
	int i, j;

	printf("here comes vote matrix\n");
	for (i = 0; i < numcells; i++) {
		for (j = 0; j < numcells; j++) {
			printf(" %3d", vote_matrix[i][j]);
		}
		printf("\n");
	}
}

#endif /* DEBUG */

/************************************************************************
 *
 *
 * ROUTINE: top_vote_getters
 *
 * DESCRIPTION:
 * Given a vote_matrix which has been filled in,
 * which has 'num' rows and columns, we need to pick the
 * top 'num' vote-getters.  We call 'top_vote_getters'
 * and are given, in its output arguments, pointers to three
 * arrays, each of which has 'num' elements pertaining
 * to a matched pair of STARS:
 *
 *       winner_votes[]    number of votes of winners, in descending order
 *       winner_index_A[]  index of star in star_array_A
 *       winner_index_B[]  index of star in star_array_B
 *
 * Thus, the pair of stars which matched in the largest number
 * of triangles will be
 *
 *       star_array_A[winner_index_A[0]]    from array A
 *       star_array_B[winner_index_A[0]]    from array B
 *
 * and the pair of stars which matched in the second-largest number
 * of triangles will be
 *
 *       star_array_A[winner_index_A[1]]    from array A
 *       star_array_B[winner_index_A[1]]    from array B
 *
 * and so on.
 *
 * RETURN:
 *    SH_SUCCESS         if all goes well
 *    SH_GENERIC_ERROR   if not
 *
 * </AUTO>
 */

static int top_vote_getters(int **vote_matrix, /* I: the 2-D array, already filled in */
int num, /* I: # of rows and cols in vote_matrix */
/*      also the number of elements in the next */
/*      three output arrays */
int **winner_votes, /* O: create this array of # of votes for the */
/*      'num' cells with the most votes */
int **winner_index_A, /* O: create this array of index into star array A */
/*      of the 'num' cells with most votes */
int **winner_index_B /* O: create this array of index into star array B */
/*      of the 'num' cells with most votes */
) {
	int i, j, k, l;
	int *w_votes; /* local ptr to (*winner_votes), for convenience */
	int *w_index_A; /* local ptr to (*winner_index_A), for convenience */
	int *w_index_B; /* local ptr to (*winner_index_B), for convenience */

	/* first, create the output arrays */
	*winner_votes = (int *) shMalloc(num * sizeof(int));
	*winner_index_A = (int *) shMalloc(num * sizeof(int));
	*winner_index_B = (int *) shMalloc(num * sizeof(int));

	/* this will simplify code inside this function */
	w_votes = *winner_votes;
	w_index_A = *winner_index_A;
	w_index_B = *winner_index_B;

	/*
	 * initialize all elements of the output arrays.  Use -1 as the
	 * index in "w_index" arrays, to indicate an empty place
	 * with no real winner.
	 */
	for (i = 0; i < num; i++) {
		w_votes[i] = 0;
		w_index_A[i] = -1;
		w_index_B[i] = -1;
	}

	/*
	 * now walk through the vote_matrix, using insertion sort to place
	 * a cell into the "winner" arrays if it has more votes than the
	 * least popular winner so far (i.e. w_votes[num - 1])
	 */
	for (i = 0; i < num; i++) {
		for (j = 0; j < num; j++) {
			if (vote_matrix[i][j] > w_votes[num - 1]) {

				/* have to insert this cell's values into the winner arrays */
				for (k = 0; k < num; k++) {
					if (vote_matrix[i][j] > w_votes[k]) {

						/* move all other winners down one place */
						for (l = num - 2; l >= k; l--) {
							w_votes[l + 1] = w_votes[l];
							w_index_A[l + 1] = w_index_A[l];
							w_index_B[l + 1] = w_index_B[l];
						}
						/* insert the new item in its place */
						w_votes[k] = vote_matrix[i][j];
						w_index_A[k] = i;
						w_index_B[k] = j;
						break;
					}
				}
			}
		}
	}

#ifdef DEBUG
	printf("  in top_vote_getters, we have top %d \n", num);
	for (i = 0; i < num; i++) {
		printf("   index_A %4d    index_B %4d    votes %4d\n", w_index_A[i],
				w_index_B[i], w_votes[i]);
	}
#endif

	return (SH_SUCCESS);
}

/************************************************************************
 *
 *
 * ROUTINE: calc_trans
 *
 * DESCRIPTION:
 * Given a set of "nbright" matched pairs of stars, which we can
 * extract from the "winner_index" and "star_array" arrays,
 * figure out a TRANS structure which takes coordinates of
 * objects in set A and transforms then into coords for set B.
 * A TRANS contains 6, 12, or 16 coefficients in equations like this:
 *
 *   if linear terms only:
 *
 *       x' = A + B*x + C*y
 *       y' = D + E*x + F*y
 *
 *   if linear plus quadratic terms,
 *
 *      x' =  A + Bx + Cy + Dxx + Exy + Fyy
 *      y' =  G + Hx + Iy + Jxx + Kxy + Lyy
 *
 *   if linear plus quadratic plus cubic,
 *
 *      x' =  A + Bx + Cy + Dxx + Exy + Fyy + Gx(xx+yy) + Hy(xx+yy)
 *      y' =  I + Jx + Ky + Lxx + Mxy + Nyy + Ox(xx+yy) + Py(xx+yy)
 *
 * where (x,y) are coords in set A and (x',y') are corresponding
 * coords in set B.
 *
 * This function simply checks the value of the TRANS 'order' field,
 * and calls the appropriate function to do the actual work.
 *
 *
 * RETURN:
 *    SH_SUCCESS           if all goes well
 *    SH_GENERIC_ERROR     if we can't find a solution
 *
 * </AUTO>
 */

static int calc_trans(int nbright, /* I: max number of stars we use in calculating */
/*      the transformation; we may cut down to */
/*      a more well-behaved subset. */
s_star *star_array_A, /* I: first array of s_star structure we match */
/*      the output TRANS takes their coords */
/*      into those of array B */
int num_stars_A, /* I: total number of stars in star_array_A */
s_star *star_array_B, /* I: second array of s_star structure we match */
int num_stars_B, /* I: total number of stars in star_array_B */
int *winner_votes, /* I: number of votes gotten by the top 'nbright' */
/*      matched pairs of stars */
int *winner_index_A, /* I: index into "star_array_A" of top */
/*      vote-getters */
int *winner_index_B, /* I: index into "star_array_B" of top */
/*      vote-getters */
TRANS *trans /* I/O: place solved coefficients into this */
/*      existing structure's fields */
/*      "order" field must be set before calling */
) {

	/*
	 * using the trans->order value, call the appropriate function
	 */
	switch (trans->order) {
	case AT_TRANS_LINEAR:
		if (calc_trans_linear(nbright, star_array_A, num_stars_A, star_array_B,
				num_stars_B, winner_votes, winner_index_A, winner_index_B,
				trans) != SH_SUCCESS) {
			shError("calc_trans: calc_trans_linear returns with error");
			return (SH_GENERIC_ERROR);
		}
		break;

	case AT_TRANS_QUADRATIC:
		if (calc_trans_quadratic(nbright, star_array_A, num_stars_A,
				star_array_B, num_stars_B, winner_votes, winner_index_A,
				winner_index_B, trans) != SH_SUCCESS) {
			shError("calc_trans: calc_trans_quadratic returns with error");
			return (SH_GENERIC_ERROR);
		}
		break;

	case AT_TRANS_CUBIC:
		if (calc_trans_cubic(nbright, star_array_A, num_stars_A, star_array_B,
				num_stars_B, winner_votes, winner_index_A, winner_index_B,
				trans) != SH_SUCCESS) {
			shError("calc_trans: calc_trans_cubic returns with error");
			return (SH_GENERIC_ERROR);
		}
		break;

	default:
		shFatal("calc_trans: called with invalid trans->order %d \n",
				trans->order);
		break;
	}

	/*
	 * we can set the 'nr' field to the number of matched pairs used
	 * to define the TRANS structure
	 */
	trans->nr = nbright;

	/* we can now compute the value of the TRANS "sig" field */
	if (calc_trans_sig(nbright, star_array_A, num_stars_A, star_array_B,
			num_stars_B, winner_votes, winner_index_A, winner_index_B, trans)
			!= SH_SUCCESS) {
		shError("calc_trans: calc_trans_sig returns with error");
		return (SH_GENERIC_ERROR);
	}

	return (SH_SUCCESS);
}

/************************************************************************
 *
 *
 * ROUTINE: alloc_matrix
 *
 * DESCRIPTION:
 * Allocate space for an NxN matrix of double values,
 * return a pointer to the new matrix.
 *
 * RETURNS:
 *   double **           pointer to new matrix
 *
 *
 * </AUTO>
 */

static double **
alloc_matrix(int n /* I: number of elements in each row and col */
) {
	int i;
	double **matrix;

	matrix = (double **) shMalloc(n * sizeof(double *));
	for (i = 0; i < n; i++) {
		matrix[i] = (double *) shMalloc(n * sizeof(double));
	}

	return (matrix);
}

/************************************************************************
 *
 *
 * ROUTINE: free_matrix
 *
 * DESCRIPTION:
 * Free the space allocated for the given nxn matrix.
 *
 * RETURNS:
 *   nothing
 *
 * </AUTO>
 */

static void free_matrix(double **matrix, /* I: pointer to 2-D array to be freed */
int n /* I: number of elements in each row and col */
) {
	int i;

	for (i = 0; i < n; i++) {
		shFree(matrix[i]);
	}
	shFree(matrix);
}

/************************************************************************
 *
 *
 * ROUTINE: print_matrix
 *
 * DESCRIPTION:
 * print out a nice picture of the given matrix.
 *
 * For debugging purposes.
 *
 * RETURNS:
 *   nothing
 *
 * </AUTO>
 */

#ifdef DEBUG

static void print_matrix(double **matrix, /* I: pointer to 2-D array to be printed */
int n /* I: number of elements in each row and col */
) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf(" %12.5e", matrix[i][j]);
		}
		printf("\n");
	}
}

#endif /* DEBUG */

/*
 * check to see if my versions of NR routines have bugs.
 * Try to invert a matrix.
 *
 * debugging only.
 */

#ifdef DEBUG3

static void test_routine(void) {
	int i, j, k, n;
	int *permutations;
	double **matrix1, **matrix2, **inverse;
	double *vector;
	double *col;
	double sum;

	fflush(stdout);
	fflush(stderr);
	n = 2;
	matrix1 = (double **) shMalloc(n * sizeof(double *));
	matrix2 = (double **) shMalloc(n * sizeof(double *));
	inverse = (double **) shMalloc(n * sizeof(double *));
	vector = (double *) shMalloc(n * sizeof(double));
	for (i = 0; i < n; i++) {
		matrix1[i] = (double *) shMalloc(n * sizeof(double));
		matrix2[i] = (double *) shMalloc(n * sizeof(double));
		inverse[i] = (double *) shMalloc(n * sizeof(double));
	}
	permutations = (int *) shMalloc(n * sizeof(int));
	col = (double *) shMalloc(n * sizeof(double));

	/* fill the matrix */
	matrix1[0][0] = 1.0;
	matrix1[0][1] = 2.0;
	matrix1[1][0] = 3.0;
	matrix1[1][1] = 4.0;

	/* fill the vector */
	for (i = 0; i < n; i++) {
		vector[i] = 0;
	}

	/* copy matrix1 into matrix2, so we can compare them later */
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			matrix2[i][j] = matrix1[i][j];
		}
	}

	/* now check */
	printf(" here comes original matrix \n");
	print_matrix(matrix1, n);

	/* now invert matrix1 */
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			inverse[i][j] = matrix1[i][j];
		}
	}
	gauss_matrix(inverse, n, vector);

	/* now check */
	printf(" here comes inverse matrix \n");
	print_matrix(inverse, n);

	/* find out if the product of "inverse" and "matrix2" is identity */
	sum = 0.0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++) {
				sum += inverse[i][k] * matrix2[k][j];
			}
			matrix1[i][j] = sum;
			sum = 0.0;
		}
	}

	printf(" here comes what we hope is identity matrix \n");
	print_matrix(matrix1, n);

	fflush(stdout);
	fflush(stderr);
}

#endif /* DEBUG3 */

/************************************************************************
 *
 *
 * ROUTINE: iter_trans
 *
 * DESCRIPTION:
 * We want to find a TRANS structures that takes coords of objects in
 * set A and transforms to coords of objects in set B.  We have a
 * a subset of 'nmatched' candidates for matched pairs of points.
 * However, some of these may be false matches.  Here's how we try
 * to eliminate them, and use all remaining true matches to derive
 * the transformation.
 *
 *    1. start with nbright matched candidate pairs of points
 *    2. choose N best pairs
 *          if recalc_flag == RECALC_NO,  set N = AT_MATCH_STARTN
 *          if recalc_flag == RECALC_YES, set N = nbright
 *       (assert that N >= AT_MATCH_REQUIRE)
 *    3. set Nr = N
 *    4. calculate a TRANS structure using the best Nr points
 *            (where "best" means "highest in winner_index" arrays)
 *    5.   transform all Nr points from coords in A to coords in B
 *    6.   calculate Euclidean square-of-distance between all Nr points
 *                             in coord system B
 *    7.   sort these Euclidean values
 *    8.   pick the AT_MATCH_PERCENTILE'th value from the sorted array
 *             (call it "sigma")
 *    9.   let Nb = number of candidate matched pairs which have
 *                       square-of-distance > AT_MATCH_NSIGMA*sigma
 *   10.   if Nb == 0, we're done -- quit
 *   11.   if Nb > 0,
 *                       remove all Nb candidates from matched pair arrays
 *                       set Nr = Nr - Nb
 *                       go to step 4
 *
 * Note that if we run out of candidate pairs, so that Nr < AT_MATCH_REQUIRE,
 * we print an error message and return SH_GENERIC_ERROR.
 *
 * The "recalc_flag" is used to distinguish two cases:
 *    if RECALC_NO,  then we are calling 'iter_trans()' with a bunch of
 *                   matches which probably contain some bad ones.
 *                   In order to prevent the bad ones from ruining the
 *                   initial calculation, we pick only the few best
 *                   on the first iteration.
 *    if RECALC_YES, we are calling 'iter_trans()' with a set of matches
 *                   which have already passed a test: namely, they are
 *                   based on a previously-determined TRANS and all matches
 *                   are within 'matchrad' in the coord system of list B.
 *                   In this case, we start out using all the matched
 *                   pairs in the very first iteration.
 *
 *
 * RETURNS:
 *   SH_SUCCESS          if we were able to determine a good TRANS
 *   SH_GENERIC_ERROR    if we couldn't
 *
 * </AUTO>
 */

static int iter_trans(int nbright, /* I: max number of stars we use in calculating */
/*      the transformation; we may cut down to */
/*      a more well-behaved subset. */
s_star *star_array_A, /* I: first array of s_star structure we match */
/*      the output TRANS takes their coords */
/*      into those of array B */
int num_stars_A, /* I: total number of stars in star_array_A */
s_star *star_array_B, /* I: second array of s_star structure we match */
int num_stars_B, /* I: total number of stars in star_array_B */
int *winner_votes, /* I: number of votes gotten by the top 'nbright' */
/*      matched pairs of stars */
/*      We may modify this array */
int *winner_index_A, /* I: index into "star_array_A" of top */
/*      vote-getters */
/*      We may modify this array */
int *winner_index_B, /* I: index into "star_array_B" of top */
/*      vote-getters */
/*      We may modify this array */
int recalc_flag, /* I: should we use only a few best pairs for */
/*      the first iteration, or all? */
int max_iterations, /* I: iterate at most this many times.  If we */
/*      reach this limit, stop iterating */
/*      and declare success */
double halt_sigma, /* I: if the residuals from solution drop to */
/*      this level, stop iterating and */
/*      declare success */
TRANS *trans /* O: place solved coefficients into this */
/*      existing structure's fields */
) {
	int i, j;
	int nr; /* number of matched pairs remaining in solution */
	int nb; /* number of bad pairs in any iteration */
	int initial_pairs;
	int is_ok;
	int required_pairs = 0, start_pairs = 0;
	int iters_so_far;
	double *dist2, *dist2_sorted;
	double xdiff, ydiff;
	double sigma;
	double max_dist2;
	double newx = 0.0, newy = 0.0;
	s_star *sa, *sb;
	s_star *a_prime; /* will hold transformed version of stars in set A */

	/*
	 * set some variables depending on the order of the fit to be
	 * performed.
	 */
	switch (trans->order) {
	case AT_TRANS_LINEAR:
		required_pairs = AT_MATCH_REQUIRE_LINEAR;
		start_pairs = AT_MATCH_STARTN_LINEAR;
		break;
	case AT_TRANS_QUADRATIC:
		required_pairs = AT_MATCH_REQUIRE_QUADRATIC;
		start_pairs = AT_MATCH_STARTN_QUADRATIC;
		break;
	case AT_TRANS_CUBIC:
		required_pairs = AT_MATCH_REQUIRE_CUBIC;
		start_pairs = AT_MATCH_STARTN_CUBIC;
		break;
	default:
		shFatal("iter_trans: invalid trans->order %d \n", trans->order);
		break;
	}

	if (nbright < required_pairs) {
#ifdef DEBUG
		printf("iter_trans: only %d items supplied, need %d\n", nbright,
				required_pairs);
#endif
		return (SH_GENERIC_ERROR);
	}

	shAssert(star_array_A != NULL);
	shAssert(star_array_B != NULL);
	shAssert(winner_votes != NULL);
	shAssert(winner_index_A != NULL);
	shAssert(winner_index_B != NULL);
	shAssert(trans != NULL);

	/* these should already have been checked, but it doesn't hurt */
	shAssert(num_stars_A >= nbright);
	shAssert(num_stars_A >= nbright);

	/*
	 * make a first guess at TRANS;
	 *    use all the pairs, if we are sure they are "safe",
	 *    or only the best 'start_pairs', if we're not sure
	 */
	if (recalc_flag == RECALC_YES) {
		initial_pairs = nbright;
	} else {
		initial_pairs = start_pairs;
	}
#ifdef DEBUG
	printf("   on initial calc, use %d pairs\n", initial_pairs);
#endif
	if (calc_trans(initial_pairs, star_array_A, num_stars_A, star_array_B,
			num_stars_B, winner_votes, winner_index_A, winner_index_B, trans)
			!= SH_SUCCESS) {
		shError("iter_trans: calc_trans returns with error");
		return (SH_GENERIC_ERROR);
	}
#ifdef DEBUG
	printf("   here comes initial TRANS\n");
	print_trans(trans);
#endif

	/*
	 * Now, we are going to enter the iteration with a set of the "best"
	 * matched pairs.  Recall that
	 * "winner_index" arrays are already sorted in decreasing order
	 * of goodness, so that "winner_index_A[0]" is the best.
	 * As we iterate, we may discard some matches, and then 'nr' will
	 * get smaller.  It must always be more than AT_MATCH_REQUIRE,
	 * or else 'calc_trans' will fail.
	 */
	nr = nbright;

	/*
	 * We're going to need an array of (at most) 'nbright' stars
	 * which hold the coordinates of stars in set A, after they've
	 * been transformed into coordinates system of set B.
	 */
	a_prime = (s_star *) shMalloc(nbright * sizeof(s_star));

	/*
	 * And this will be an array to hold the Euclidean square-of-distance
	 * between a transformed star from set A and its partner from set B.
	 *
	 * "dist2_sorted" is a copy of the array which we'll sort .. but we need
	 * to keep the original order, too.
	 */
	dist2 = (double *) shMalloc(nbright * sizeof(double));
	dist2_sorted = (double *) shMalloc(nbright * sizeof(double));

	/*
	 * we don't allow any candidate matches which cause the stars to
	 * differ by more than this much in the common coord system.
	 */
	max_dist2 = AT_MATCH_MAXDIST * AT_MATCH_MAXDIST;

	/*
	 * now, we enter a loop that may execute several times.
	 * We calculate the transformation for current 'nr' best points,
	 * then check to see if we should throw out any matches because
	 * the resulting transformed coordinates are too discrepant.
	 * We break out of this loop near the bottom, with a status
	 * provided by "is_ok"
	 *
	 *       is_ok = 1              all went well, can return success
	 *       is_ok = 0              we failed for some reason.
	 */
	is_ok = 1;

	iters_so_far = 0;
	while (iters_so_far < max_iterations) {

#ifdef DEBUG
		printf("iter_trans: at top of loop, nr=%4d iters_so_far=%4d\n", nr,
				iters_so_far);
#endif

		nb = 0;

		/*
		 * apply the TRANS to the A stars in all 'nr' matched pairs.
		 * we make a new set of s_stars with the transformed coordinates,
		 * called "a_prime".
		 */
		for (i = 0; i < nr; i++) {
			sa = &(star_array_A[winner_index_A[i]]);
			if (calc_trans_coords(sa, trans, &newx, &newy) != SH_SUCCESS) {
				shError("iter_trans: calc_trans_coords fails");
				return (SH_GENERIC_ERROR);
			}
			a_prime[i].x = newx;
			a_prime[i].y = newy;
		}

		/*
		 * calculate the square-of-distance between a transformed star
		 * (from set A) and its partner from set B, in the coordinate system
		 * of set B.
		 */
		for (i = 0; i < nr; i++) {
			sb = &(star_array_B[winner_index_B[i]]);
			xdiff = a_prime[i].x - sb->x;
			ydiff = a_prime[i].y - sb->y;
			dist2[i] = (xdiff * xdiff + ydiff * ydiff);
			dist2_sorted[i] = dist2[i];
#ifdef DEBUG
			printf(
					"   match %3d  (%12.5e,%12.5e) vs. (%12.5e,%12.5e)  d2=%12.6e\n",
					i, a_prime[i].x, a_prime[i].y, sb->x, sb->y, dist2[i]);
#endif
		}

		/*
		 * sort the array of square-of-distances
		 */
		qsort((char *) dist2_sorted, nr, sizeof(double), (PFI) compare_double);

		/*
		 * now, check to see if any matches have dist2 > max_dist2.
		 * If so,
		 *
		 *     - remove them from the winner_votes and winner_index arrays
		 *     - decrement 'nr'
		 *     - also decrement the loop counter 'i', because we're going
		 *            to move up all items in the "winner" and "dist2" arrays
		 *            as we discard the bad match
		 *     - increment 'nb'
		 */
		for (i = 0; i < nr; i++) {
			if (dist2[i] > max_dist2) {

				/*
				 * remove the entry for the "bad" match from the "winner" arrays
				 * and from the "dist2" array
				 */
#ifdef DEBUG
				printf("  removing old match with d2=%9.4e\n", dist2[i]);
#endif
				for (j = i + 1; j < nr; j++) {
					winner_votes[j - 1] = winner_votes[j];
					winner_index_A[j - 1] = winner_index_A[j];
					winner_index_B[j - 1] = winner_index_B[j];
					dist2[j - 1] = dist2[j];
				}

				/*
				 * and modify our counters of "remaining good matches" and
				 * "bad matches this time", too.
				 */
				nr--; /* one fewer good match remains */
				nb++; /* one more bad match during this iteration */

				/*
				 * and decrement 'i', too, since we must moved element
				 * i+1 to the place i used to be, and we must check _it_.
				 */
				i--;
			}
		}
#ifdef DEBUG
		printf("   nr now %4d, nb now %4d\n", nr, nb);
#endif

		/*
		 * pick the square-of-distance which occurs at the AT_MATCH_PERCENTILE
		 * place in the sorted array.  Call this value "sigma".  We'll clip
		 * any matches that are more than AT_MATCH_NSIGMA*"sigma".
		 *
		 * However, if we have fewer than 2 objects, don't bother with this
		 * step -- just set "sigma" equal to 0 and prepare for later
		 * failure....
		 */
		if (nr < 2) {
			sigma = 0.0;
#ifdef DEBUG
			printf("   sigma = %10.5e  (only %d matches) \n", sigma, nr);
#endif
		} else {
			sigma = find_percentile(dist2_sorted, nr,
					(double) AT_MATCH_PERCENTILE);
#ifdef DEBUG
			printf("   sigma = %10.5e\n", sigma);
#endif
		}

		/*
		 * If the current "sigma" value is less than the "halt_sigma" value,
		 * then we have succeeded.   Stop iterating.
		 */
		if (sigma <= halt_sigma) {
#ifdef DEBUG
			printf("   SUCCESS  sigma = %10.5e  <  halt_sigma %10.5e \n", sigma,
					halt_sigma);
#endif
			is_ok = 1;
			break;
		}

		/*
		 * now, check to see if any matches have dist2 > AT_MATCH_NSIGMA*sigma.
		 * If so,
		 *
		 *     - remove them from the winner_votes and winner_index arrays
		 *     - decrement 'nr'
		 *     - also decrement the loop counter 'i', because we're going
		 *            to move up all items in the "winner" and "dist2" arrays
		 *            as we discard the bad match
		 *     - increment 'nb'
		 */
		for (i = 0; i < nr; i++) {
			if (dist2[i] > AT_MATCH_NSIGMA * sigma) {

				/*
				 * remove the entry for the "bad" match from the "winner" arrays
				 * and from the "dist2" array
				 */
#ifdef DEBUG
				printf("  removing old match with d2=%9.4e\n", dist2[i]);
#endif
				for (j = i + 1; j < nr; j++) {
					winner_votes[j - 1] = winner_votes[j];
					winner_index_A[j - 1] = winner_index_A[j];
					winner_index_B[j - 1] = winner_index_B[j];
					dist2[j - 1] = dist2[j];
				}

				/*
				 * and modify our counters of "remaining good matches" and
				 * "bad matches this time", too.
				 */
				nr--; /* one fewer good match remains */
				nb++; /* one more bad match during this iteration */

				/*
				 * and decrement 'i', too, since we must moved element
				 * i+1 to the place i used to be, and we must check _it_.
				 */
				i--;
			}
		}
#ifdef DEBUG
		printf("   nr now %4d, nb now %4d\n", nr, nb);
#endif

		/*
		 * Okay, let's evaluate what has happened so far:
		 *    - if nb == 0, then all remaining matches are good
		 *    - if nb > 0, we need to iterate again
		 *    - if nr < required_pairs, we've thrown out too many points,
		 *             and must quit in shame
		 */
		if (nb == 0) {
#ifdef DEBUG
			printf("   SUCCESS  nb = 0, no more pairs to discard \n");
#endif
			is_ok = 1;
			break;
		}

		if (nr < required_pairs) {
			shDebug(AT_MATCH_ERRLEVEL,
					"iter_trans: only %d points remain, fewer than %d required",
					nr, required_pairs);
			is_ok = 0;
			break;
		}

		/*
		 * calculate the TRANS for the remaining set of matches
		 */
#ifdef DEBUG
		printf("   on this iter,    use %d pairs\n", nr);
#endif
		if (calc_trans(nr, star_array_A, num_stars_A, star_array_B, num_stars_B,
				winner_votes, winner_index_A, winner_index_B, trans)
				!= SH_SUCCESS) {
			shError("iter_trans: calc_trans returns with error");
			return (SH_GENERIC_ERROR);
		}

#ifdef DEBUG
		printf("   here comes latest TRANS\n");
		print_trans(trans);
#endif

		iters_so_far++;

	}

	if (iters_so_far == max_iterations) {
#ifdef DEBUG
		printf("   SUCCESS(?): iters_so_far %d = max_iterations\n",
				iters_so_far);
#endif
	}

	/*
	 * Here we summarize the result of our work in two of the
	 *   elements of the TRANS structure:
	 *         trans->nr   =  number of pairs used to find transformation
	 *         trans->sig  =  stdev of separation of matching pairs,
	 *                                in units of coord system B
	 */
	trans->nr = nr;
	trans->sig = find_percentile(dist2_sorted, nr, ONE_STDEV_PERCENTILE);

	/*
	 * free up the arrays we allocated
	 */
	shFree(a_prime);
	shFree(dist2);
	shFree(dist2_sorted);

	/*
	 * and decide whether we succeeded, or failed
	 */
	if (is_ok == 0) {
		return (SH_GENERIC_ERROR);
	} else {
		return (SH_SUCCESS);
	}
}

/************************************************************************
 *
 *
 * ROUTINE: compare_double
 *
 * DESCRIPTION:
 * Given pointers to two double numbers, return the comparison.
 * Used by "iter_trans"
 *
 * RETURN:
 *    1                  if first double is larger than second
 *    0                  if the two are equal
 *   -1                  if first double is smaller than second
 *
 * </AUTO>
 */

static int compare_double(double *f1, /* I: compare size of FIRST double value */
double *f2 /* I:  ... with SECOND double value  */
) {
	shAssert((f1 != NULL) && (f2 != NULL));

	if (*f1 > *f2) {
		return (1);
	}
	if (*f1 < *f2) {
		return (-1);
	}
	return (0);
}

/************************************************************************
 *
 *
 * ROUTINE: find_percentile
 *
 * DESCRIPTION:
 * Given an array of 'num' double values, which have been
 * sorted, find the value corresponding to the value which is at
 * the 'perc'th percentile in the list array.  Return this value.
 *
 * RETURN:
 *   double                value of the number at 'perc'th percentile in array
 *
 * </AUTO>
 */

static double find_percentile(double *array, /* I: look in this SORTED array */
int num, /* I: which has this many elements */
double perc /* I: for entry at this percentile */
) {
	int index;

	shAssert(array != NULL);
	shAssert(num > 0);
	shAssert((perc > 0.0) && (perc <= 1.0));

	index = (int) floor(num * perc + 0.5);
	if (index >= num) {
		index = num - 1;
	}
	return (array[index]);
}

/************************************************************************
 *
 *
 * ROUTINE: calc_trans_coords
 *
 * DESCRIPTION:
 * Given a single s_star structure, apply the
 * given TRANS structure to its coordinates.
 * Place the converted coordinates into the given output args
 * "newx" and "newy".
 *
 * We use the trans->order value to flag the type of transformation
 * to calculate.
 *
 *
 * RETURN:
 *   SH_SUCCESS             if all goes well
 *   SH_GENERIC_ERROR       if some problem occurs
 *
 * </AUTO>
 */

static int calc_trans_coords(s_star *star, /* I: use this STAR's coords as input */
TRANS *trans, /* I: contains coefficients of transformation */
double *newx, /* O: contains output x coord */
double *newy /* O: contains output y coord */
) {
	double rsquared;

	shAssert(star != NULL);
	shAssert(trans != NULL);

	switch (trans->order) {
	case AT_TRANS_LINEAR:
		*newx = trans->a + trans->b * star->x + trans->c * star->y;
		*newy = trans->d + trans->e * star->x + trans->f * star->y;
		break;

	case AT_TRANS_QUADRATIC:
		*newx = trans->a + trans->b * star->x + trans->c * star->y
				+ trans->d * star->x * star->x + trans->e * star->x * star->y
				+ trans->f * star->y * star->y;
		*newy = trans->g + trans->h * star->x + trans->i * star->y
				+ trans->j * star->x * star->x + trans->k * star->x * star->y
				+ trans->l * star->y * star->y;
		break;

	case AT_TRANS_CUBIC:
		rsquared = star->x * star->x + star->y * star->y;
		*newx = trans->a + trans->b * star->x + trans->c * star->y
				+ trans->d * star->x * star->x + trans->e * star->x * star->y
				+ trans->f * star->y * star->y + trans->g * star->x * rsquared
				+ trans->h * star->y * rsquared;

		*newy = trans->i + trans->j * star->x + trans->k * star->y
				+ trans->l * star->x * star->x + trans->m * star->x * star->y
				+ trans->n * star->y * star->y + trans->o * star->x * rsquared
				+ trans->p * star->y * rsquared;
		break;

	default:
		shFatal("calc_trans_coords: given invalid trans->order %d \n",
				trans->order);
		break;
	}

	return (SH_SUCCESS);
}

/************************************************************************
 *
 *
 * ROUTINE: apply_trans
 *
 * DESCRIPTION:
 * Given an array of 'num_stars' s_star structures, apply the
 * given TRANS structure to the coordinates of each one.
 *
 *
 * RETURN:
 *   SH_SUCCESS             if all goes well
 *   SH_GENERIC_ERROR       if some problem occurs
 *
 * </AUTO>
 */

static int apply_trans(s_star *star_array, /* I/O: array of structures to modify */
int num_stars, /* I: number of stars in the array */
TRANS *trans /* I: contains coefficients of transformation */
) {
	int i;
	double newx = 0.0, newy = 0.0;
	s_star *sp;

	if (num_stars == 0) {
		return (SH_SUCCESS);
	}
	shAssert(star_array != NULL);
	shAssert(trans != NULL);

	for (i = 0; i < num_stars; i++) {
		sp = &(star_array[i]);
		if (calc_trans_coords(sp, trans, &newx, &newy) != SH_SUCCESS) {
			shError("apply_trans: calc_trans_coords fails");
			return (SH_GENERIC_ERROR);
		}
		sp->x = newx;
		sp->y = newy;
	}

	return (SH_SUCCESS);
}

/***************************************************************************
 *
 *
 * ROUTINE: double_sort_by_match_id
 *
 * DESCRIPTION:
 * sort all the elements of the first array of "s_star" in increasing
 * order by "match_id" value.  Also, reorder the
 * elements of the _second_ array in exactly the same way, so that
 * the elements of both array which matched BEFORE the sorting
 * will match again _after_ the sorting.
 *
 * return:
 *   SH_SUCCESS                 if all goes well
 *   SH_GENERIC_ERROR           if not
 *
 * </AUTO>
 */

static int double_sort_by_match_id(s_star *star_array_A, /* I/O: array to be sorted */
int num_stars_A, /* I: number of stars in array A */
s_star *star_array_B, /* I/O: array to be re-ordered just as A */
int num_stars_B /* I: number of stars in array B */
) {
	int i;
	struct s_star *temp_array;
	struct s_star *sb, *stemp;

	shAssert(num_stars_A == num_stars_B);
	if (num_stars_A == 0) {
		return (SH_SUCCESS);
	}
	shAssert(star_array_A != NULL);
	shAssert(star_array_B != NULL);

	/*
	 * first, let's set the "index" field of each element of each
	 * star_array its position in the array.
	 */
	for (i = 0; i < num_stars_A; i++) {
		star_array_A[i].index = i;
		star_array_B[i].index = i;
	}

	/*
	 * next, we create a temporary array of the same size as A and B.
	 */
	temp_array = (s_star *) shMalloc(num_stars_A * sizeof(s_star));

	/*
	 * Now, the two arrays A and B are currently arranged so that
	 * star_array_A[i] matches star_array_B[i].  We want to sort
	 * star_array_A, and re-arrange star_array_B so that the
	 * corresponding elements still match up afterwards.
	 *
	 *    - sort star_array_A
	 *    - loop i through sorted star_array_A
	 *           copy star_array_B element matching star_array_A[i]
	 *                                                 into temp_array[i]
	 *    - loop i through star_array_B
	 *           copy temp_array[i] into star_array_B[i]
	 *
	 *    - delete temp_array
	 *
	 * We end up with star_array_A sorted by "x", and star_array_B
	 * re-arranged in exactly the same order.
	 */

	sort_star_by_match_id(star_array_A, num_stars_A);
	for (i = 0; i < num_stars_A; i++) {
		sb = &(star_array_B[star_array_A[i].index]);
		shAssert(sb != NULL);
		stemp = &(temp_array[i]);
		shAssert(stemp != NULL);
		copy_star(sb, stemp);
	}

	/*
	 * now copy the elements of the temp_array back into star_array_B
	 */
	for (i = 0; i < num_stars_A; i++) {
		sb = &(star_array_B[i]);
		shAssert(sb != NULL);
		stemp = &(temp_array[i]);
		shAssert(stemp != NULL);
		copy_star(stemp, sb);
	}

	/*
	 * and we're done!  Delete the temporary array
	 */
	free_star_array(temp_array);

	return (SH_SUCCESS);
}

/***************************************************************************
 *
 *
 * ROUTINE: match_arrays_slow
 *
 * DESCRIPTION:
 * given two arrays of s_stars [A and B], find all matching elements,
 * where a match is coincidence of centers to within "radius" pixels.
 *
 * Use a slow, but sure, algorithm (and an inefficient implementation,
 * I'm sure.  As of 1/18/96, trying for correctness, not speed).
 *
 * We will match objects from A --> B.  It is possible to have several
 * As that match to the same B:
 *
 *           A1 -> B5   and A2 -> B5
 *
 * This function finds such multiple-match items and deletes all but
 * the closest of the matches.
 *
 * This array creates 4 new arrays of s_stars, and returns a pointer
 * to each array, as well as the number of stars in each array.
 *
 * place the elems of A that are matches into output array J
 *                    B that are matches into output array K
 *                    A that are not matches into output array L
 *                    B that are not matches into output array M
 *
 * return: SH_SUCCESS          if all goes well
 *         SH_GENERIC_ERROR    if not
 *
 * </AUTO>
 */

static int match_arrays_slow(s_star *star_array_A, /* I: first array of s_stars to be matched */
int num_stars_A, /* I: number of stars in A */
s_star *star_array_B, /* I: second array of s_stars to be matched */
int num_stars_B, /* I: number of stars in B */
double radius, /* I: matching radius */
s_star **star_array_J, /* O: all stars in A which match put in here */
int *num_stars_J, /* O: number of stars in output array J */
s_star **star_array_K, /* O: all stars in B which match put in here */
int *num_stars_K, /* O: number of stars in output array K */
s_star **star_array_L, /* O: all stars in A which don't match put here */
int *num_stars_L, /* O: number of stars in output array L */
s_star **star_array_M, /* O: all stars in B which don't match put here */
int *num_stars_M /* O: number of stars in output array M */
) {
	double Ax, Ay, Bx, By;
	double dist, limit;
	int posA, posB;
	int current_num_J, current_num_K;
	double deltax, deltay;
	double Axm, Axp, Aym, Ayp;
	s_star *sa, *sb;

#ifdef DEBUG
	printf("entering match_arrays_slow ");
#endif

	/*
	 * first, we create each of the 4 output arrays.  We start with
	 * each as big as the input arrays, but we'll shrink them down
	 * to their proper sizes before we return.
	 */
	*star_array_J = (s_star *) shMalloc(num_stars_A * sizeof(s_star));
	*num_stars_J = num_stars_A;
	*star_array_K = (s_star *) shMalloc(num_stars_B * sizeof(s_star));
	*num_stars_K = num_stars_B;
	*star_array_L = (s_star *) shMalloc(num_stars_A * sizeof(s_star));
	*num_stars_L = num_stars_A;
	*star_array_M = (s_star *) shMalloc(num_stars_B * sizeof(s_star));
	*num_stars_M = num_stars_B;

	/*
	 * make some sanity checks
	 */
	shAssert(num_stars_A >= 0);
	shAssert(num_stars_B >= 0);
	if ((num_stars_A == 0) || (num_stars_B == 0)) {
		return (SH_SUCCESS);
	}
	shAssert(star_array_A != NULL);
	shAssert(star_array_B != NULL);

	/*
	 * First, we sort arrays A and B by their "x" coordinates,
	 * to facilitate matching.
	 */
	sort_star_by_x(star_array_A, num_stars_A);
	sort_star_by_x(star_array_B, num_stars_B);

	/*
	 * We copy array A into L, and array B into M.
	 * We will remove all non-matching elements from these
	 * output arrays later on in this function.
	 */

	copy_star_array(star_array_A, *star_array_L, num_stars_A);
	copy_star_array(star_array_B, *star_array_M, num_stars_B);

	/*
	 * this is the largest distance that two stars can be from
	 * each other and still be a match.
	 */
	limit = radius * radius;

	/*
	 * the first step is to go slowly through array A, checking against
	 * every object in array B.  If there's a match, we copy the matching
	 * elements onto lists J and K, respectively.  We do NOT check
	 * yet to see if there are multiply-matched elements.
	 *
	 * This implementation could be speeded up a LOT by sorting the
	 * two arrays in "x" and then making use of the information to check
	 * only stars which are close to each other in "x".  Do that
	 * some time later.... MWR 1/18/96.
	 */
#ifdef DEBUG
	printf(" size of array A is %d, array B is %d\n", num_stars_A, num_stars_B);
	printf(" about to step through array A looking for matches\n");
#endif

	current_num_J = 0;
	current_num_K = 0;

	for (posA = 0; posA < num_stars_A; posA++) {

		shAssert((sa = &(star_array_A[posA])) != NULL);
		Ax = sa->x;
		Ay = sa->y;

		Axm = Ax - radius;
		Axp = Ax + radius;
		Aym = Ay - radius;
		Ayp = Ay + radius;

		for (posB = 0; posB < num_stars_B; posB++) {

			shAssert((sb = &(star_array_B[posB])) != NULL);
			Bx = sb->x;
			By = sb->y;

			/* check quickly to see if we can avoid a multiply */
			if ((Bx < Axm) || (Bx > Axp) || (By < Aym) || (By > Ayp)) {
				continue;
			}

			/* okay, we actually have to calculate a distance here. */
			deltax = Ax - Bx;
			deltay = Ay - By;
			dist = deltax * deltax + deltay * deltay;
			if (dist < limit) {

				/*
				 * we have a match (at least, a possible match).  So, copy
				 * objA onto listJ and objB onto listK.  But do NOT remove
				 * these objects from listA and listB!  We may end up
				 * matching another objA to the same objB later on, and
				 * we will continue trying to match this same objA to other
				 * objBs.
				 */
				add_element(sa, star_array_J, num_stars_J, &current_num_J);
				add_element(sb, star_array_K, num_stars_K, &current_num_K);

			}
		}
	}

	/*
	 * at this point, let's re-set "*num_stars_J" to the proper number.
	 * Recall that the "add_element" function may increase "*num_stars_J"
	 * by factors of 2, while the variable "current_num_J" keeps track
	 * of the actual number of stars in the array.  It ought to be the
	 * case that
	 *              num_stars_J <= *num_stars_J
	 *
	 * and likewise for K.
	 */
	*num_stars_J = current_num_J;
	*num_stars_K = current_num_K;

#ifdef DEBUG
	printf(" done with stepping through array A \n");
	printf(" array J has %d, array K has %d \n", current_num_J, current_num_K);
#endif

#ifdef DEBUG
	/* for debugging only */
	for (posA = 0; posA < *num_stars_J; posA++) {
		sa = &((*star_array_J)[posA]);
		sb = &((*star_array_K)[posA]);
		printf(" %4d  J: %4d (%8.2f, %8.2f)  K: %4d (%8.2f, %8.2f) \n", posA,
				sa->match_id, sa->x, sa->y, sb->match_id, sb->x, sb->y);
	}
#endif

	/*
	 * at this point, all _possible_ matches have been placed into
	 * corresponding elements of arrays J and K.  Now, we go through
	 * array J to find elements which appear more than once.  We'll
	 * resolve them by throwing out all but the closest match.
	 */

	/*
	 * first, sort array J by the "match_id" values.  This allows us to find
	 * repeated elements easily.  Re-order array K in exactly the same
	 * way, so matching elements still match.
	 */
#ifdef DEBUG
	printf(" sorting array J by match_id\n");
#endif
	if (double_sort_by_match_id(*star_array_J, *num_stars_J, *star_array_K,
			*num_stars_K) != SH_SUCCESS) {
		shError("match_arrays_slow: can't sort array J");
		return (SH_GENERIC_ERROR);
	}
#ifdef DEBUG
	for (posA = 0; posA < *num_stars_J; posA++) {
		sa = &((*star_array_J)[posA]);
		sb = &((*star_array_K)[posA]);
		printf(" %4d  J: %4d (%8.2f, %8.2f)  K: %4d (%8.2f, %8.2f) \n", posA,
				sa->match_id, sa->x, sa->y, sb->match_id, sb->x, sb->y);
	}
#endif

	/*
	 * now remove repeated elements from array J, keeping the closest matches
	 */
#ifdef DEBUG
	printf(" before remove_repeated_elements, array J has %d\n", *num_stars_J);
#endif
	if (remove_repeated_elements(*star_array_J, num_stars_J, *star_array_K,
			num_stars_K) != SH_SUCCESS) {
		shError(
				"match_arrays_slow: remove_repeated_elements fails for array J");
		return (SH_GENERIC_ERROR);
	}
#ifdef DEBUG
	printf(" after remove_repeated_elements, array J has %d\n", *num_stars_J);
	for (posA = 0; posA < *num_stars_J; posA++) {
		sa = &((*star_array_J)[posA]);
		sb = &((*star_array_K)[posA]);
		printf(" %4d  J: %4d (%8.2f, %8.2f)  K: %4d (%8.2f, %8.2f) \n", posA,
				sa->match_id, sa->x, sa->y, sb->match_id, sb->x, sb->y);
	}
#endif
	shAssert(*num_stars_J == *num_stars_K);

	/*
	 * next, do the same for array K: sort it by "match_id"
	 * (and re-arrange array J to match),
	 * then find and remove any
	 * repeated elements, keeping only the closest matches.
	 */
#ifdef DEBUG
	printf(" sorting array K by match_id\n");
#endif
	if (double_sort_by_match_id(*star_array_K, *num_stars_K, *star_array_J,
			*num_stars_J) != SH_SUCCESS) {
		shError("match_arrays_slow: can't sort array K");
		return (SH_GENERIC_ERROR);
	}
#ifdef DEBUG
	for (posA = 0; posA < *num_stars_J; posA++) {
		sa = &((*star_array_J)[posA]);
		sb = &((*star_array_K)[posA]);
		printf(" %4d  J: %4d (%8.2f, %8.2f)  K: %4d (%8.2f, %8.2f) \n", posA,
				sa->match_id, sa->x, sa->y, sb->match_id, sb->x, sb->y);
	}
#endif

#ifdef DEBUG
	printf(" before remove_repeated_elements, array K has %d\n", *num_stars_K);
	for (posA = 0; posA < *num_stars_J; posA++) {
		sa = &((*star_array_J)[posA]);
		sb = &((*star_array_K)[posA]);
		printf(" %4d  J: %4d (%8.2f, %8.2f)  K: %4d (%8.2f, %8.2f) \n", posA,
				sa->match_id, sa->x, sa->y, sb->match_id, sb->x, sb->y);
	}
#endif
	if (remove_repeated_elements(*star_array_K, num_stars_K, *star_array_J,
			num_stars_J) != SH_SUCCESS) {
		shError(
				"match_arrays_slow: remove_repeated_elements fails for array K");
		return (SH_GENERIC_ERROR);
	}
#ifdef DEBUG
	printf(" after remove_repeated_elements, arrary K has %d\n", *num_stars_K);
	for (posA = 0; posA < *num_stars_J; posA++) {
		sa = &((*star_array_J)[posA]);
		sb = &((*star_array_K)[posA]);
		printf(" %4d  J: %4d (%8.2f, %8.2f)  K: %4d (%8.2f, %8.2f) \n", posA,
				sa->match_id, sa->x, sa->y, sb->match_id, sb->x, sb->y);
	}
#endif
	shAssert(*num_stars_J == *num_stars_K);

	/*
	 * finally, we have unique set of closest-pair matching elements
	 * in arrays J and K.  Now we can remove any element from array L
	 * which appears in array J, and remove any element from array M
	 * which appears in array K.  First, we'll sort arrays L and M
	 * to make the comparisons easier.
	 */
#ifdef DEBUG
	printf(" sorting array L \n");
#endif
	sort_star_by_match_id(*star_array_L, *num_stars_L);
#ifdef DEBUG
	printf(" sorting array M \n");
#endif
	sort_star_by_match_id(*star_array_M, *num_stars_M);

	/*
	 * Recall that array K is already sorted by "match_id", but that
	 * we may have thrown J out of order when we forced it to follow
	 * the sorting of K.  So, first we'll sort J by "match_id",
	 * (and re-order K match it), then we can remove repeated elements
	 * from L easily.
	 */
#ifdef DEBUG
	printf(" sorting array J by match_id\n");
#endif
	if (double_sort_by_match_id(*star_array_J, *num_stars_J, *star_array_K,
			*num_stars_K) != SH_SUCCESS) {
		shError("match_arrays_slow: can't sort array J");
		return (SH_GENERIC_ERROR);
	}
#ifdef DEBUG
	printf(" after double_sort_by_match_id (J, K)\n");
	for (posA = 0; posA < *num_stars_J; posA++) {
		sa = &((*star_array_J)[posA]);
		sb = &((*star_array_K)[posA]);
		printf(" %4d  J: %4d (%8.2f, %8.2f)  K: %4d (%8.2f, %8.2f) \n", posA,
				sa->match_id, sa->x, sa->y, sb->match_id, sb->x, sb->y);
	}
#endif
	/*
	 * now remove elements from array L which appear in array J
	 */
#ifdef DEBUG
	printf(" before remove_same_elements, array L has %d\n", *num_stars_L);
	for (posA = 0; posA < *num_stars_L; posA++) {
		sa = &((*star_array_L)[posA]);
		sb = &((*star_array_M)[posA]);
		printf(" %4d  L: %4d (%8.2f, %8.2f)  M: %4d (%8.2f, %8.2f) \n", posA,
				sa->match_id, sa->x, sa->y, sb->match_id, sb->x, sb->y);
	}
#endif
	remove_same_elements(*star_array_J, *num_stars_J, *star_array_L,
			num_stars_L);
#ifdef DEBUG
	printf(" after remove_same_elements, array L has %d\n", *num_stars_L);
	for (posA = 0; posA < *num_stars_L; posA++) {
		sa = &((*star_array_L)[posA]);
		sb = &((*star_array_M)[posA]);
		printf(" %4d  L: %4d (%8.2f, %8.2f)  M: %4d (%8.2f, %8.2f) \n", posA,
				sa->match_id, sa->x, sa->y, sb->match_id, sb->x, sb->y);
	}
#endif

	/*
	 * Recall that we threw K out of order when we forced it to follow
	 * the sorting of J.  So, we'll sort K by "match_id",
	 * (and re-order J match it), then we can remove repeated elements
	 * from M easily.
	 */
#ifdef DEBUG
	printf(" sorting array K by match_id\n");
#endif
	if (double_sort_by_match_id(*star_array_K, *num_stars_K, *star_array_J,
			*num_stars_J) != SH_SUCCESS) {
		shError("match_arrays_slow: can't sort array K");
		return (SH_GENERIC_ERROR);
	}
#ifdef DEBUG
	printf(" after double_sort_by_match_id (K, J)\n");
	for (posA = 0; posA < *num_stars_J; posA++) {
		sa = &((*star_array_J)[posA]);
		sb = &((*star_array_K)[posA]);
		printf(" %4d  J: %4d (%8.2f, %8.2f)  K: %4d (%8.2f, %8.2f) \n", posA,
				sa->match_id, sa->x, sa->y, sb->match_id, sb->x, sb->y);
	}
#endif
	/*
	 * and remove elements from array M which appear in array K
	 */
#ifdef DEBUG
	printf(" before remove_same_elements, array M has %d\n", *num_stars_M);
	for (posA = 0; posA < *num_stars_L; posA++) {
		sa = &((*star_array_L)[posA]);
		sb = &((*star_array_M)[posA]);
		printf(" %4d  L: %4d (%8.2f, %8.2f)  M: %4d (%8.2f, %8.2f) \n", posA,
				sa->match_id, sa->x, sa->y, sb->match_id, sb->x, sb->y);
	}
#endif
	remove_same_elements(*star_array_K, *num_stars_K, *star_array_M,
			num_stars_M);
#ifdef DEBUG
	printf(" after remove_same_elements, array M has %d\n", *num_stars_M);
	for (posA = 0; posA < *num_stars_L; posA++) {
		sa = &((*star_array_L)[posA]);
		sb = &((*star_array_M)[posA]);
		printf(" %4d  L: %4d (%8.2f, %8.2f)  M: %4d (%8.2f, %8.2f) \n", posA,
				sa->match_id, sa->x, sa->y, sb->match_id, sb->x, sb->y);
	}
#endif

	return (SH_SUCCESS);
}

/**************************************************************************
 *
 *
 * ROUTINE: add_element
 *
 * DESCRIPTION:
 * We are given a pointer to s_star, an array of "total_num" s_stars,
 * and a count of the current number of s_stars set in the array.
 *
 * We want to copy the contents of the single star into
 * the "current_num"'th element of the array.
 *
 *   If current_num < total_num,    just perform copy,
 *                                  increment current_num
 *
 *   If current_num == total_num,   we must allocate more space in array
 *                                  allocate an array 2x as big as total_num
 *                                  copy existing elements into new array
 *                                  copy new element into new array
 *                                  free old array
 *                                  make old array pointer point to new array
 *                                  increment current_num
 *
 * We could avoid all this by using linked lists, but I think
 * that we will only rarely have to increase the size of an array,
 * and never increase its size more than once.  So this isn't so bad.
 *
 * RETURN: SH_SUCCESS          if all goes well
 *         SH_GENERIC_ERROR    if not
 *
 */

static int add_element(s_star *new_star, /* I: want to copy this into next slot in array */
s_star **star_array, /* I/O: will copy into this array */
/*       if necessary, will allocate a new array, */
/*       copy entire contents into it, including */
/*       new_star, then free the old array */
int *total_num, /* I/O: total number of stars allocated in */
/*       star_array.  We may increase this if */
/*       we have to extend star_array */
int *current_num /* I/O: current number of stars in star_array */
/*       which have been set.  This number should */
/*       always increase by 1 if we succeed in */
/*       in adding the "new_star" */
) {
	int num;
	s_star *new_array;

	shAssert(new_star != NULL);
	shAssert((star_array != NULL) && (*star_array != NULL));
	shAssert((total_num != NULL) && (*total_num >= 0));
	shAssert((current_num != NULL) && (*current_num >= 0));

	/*
	 * check for the easy case: if current_num < total_num, we can
	 * just set star_array[current_num] and increment current_num.
	 */
	if (*current_num < *total_num) {
		copy_star(new_star, &((*star_array)[*current_num]));
		(*current_num)++;
	} else if (*current_num == *total_num) {

		/*
		 * this is the tricky case, in which we have to allocate space
		 * for a larger array, copy all existing elements, and then
		 * copy over the new_star.
		 */
		num = (*total_num) * 2;
		new_array = (s_star *) shMalloc(num * sizeof(s_star));
		copy_star_array((*star_array), new_array, (*total_num));
		free_star_array(*star_array);
		*star_array = new_array;
		*total_num = num;
		copy_star(new_star, &((*star_array)[*current_num]));
		(*current_num)++;
	} else {

		/*
		 * this should never occur!
		 */
		shAssert(0);
	}

	return (SH_SUCCESS);
}

/*********************************************************************
 *
 * ROUTINE: remove_repeated_elements
 *
 * DESCRIPTION:
 * step through the first array argument, star_array_1, checking for
 * successive elements which are the same. for each such pair, calculate the
 * distance between the matching elements of objects in arrays 1 and 2.
 * Throw the less-close pair out of the two array, modifying the number
 * of elements in each accordingly (and moving all other elements
 * up one place in the array).
 *
 * The two arrays must have the same number of elements,
 * and array 1 must already have been sorted by the "match_id" field.
 *
 * RETURN:
 *    SH_SUCCESS              if all goes well
 *    SH_GENERIC_ERROR        if something goes wrong
 *
 */

static int remove_repeated_elements(s_star *star_array_1, /* I/O: look in this array for repeats */
int *num_stars_1, /* I/O: number of stars in array 1 */
s_star *star_array_2, /* I/O: do to this array what we do to array 1 */
int *num_stars_2 /* I/O: number of stars in array 2 */
) {
	int pos1, pos2;
	double thisdist, lastdist;
	s_star *s1, *s2;
	s_star *last1, *last2;

	shAssert(star_array_1 != NULL);
	shAssert(star_array_2 != NULL);
	shAssert(*num_stars_1 == *num_stars_2);

	pos1 = 0;
	pos2 = 0;

	last1 = NULL;
	last2 = NULL;
	while (pos1 < *num_stars_1) {

		s1 = &(star_array_1[pos1]);
		s2 = &(star_array_2[pos2]);
		if ((s1 == NULL) || (s2 == NULL)) {
			shError("remove_repeated_elements: missing elem in array 1 or 2");
			return (SH_GENERIC_ERROR);
		}

		if (last1 == NULL) {
			last1 = s1;
			last2 = s2;
		} else if (s1->match_id == last1->match_id) {

			/* there is a repeated element.  We must find the closer match */
			thisdist = (s1->x - s2->x) * (s1->x - s2->x)
					+ (s1->y - s2->y) * (s1->y - s2->y);
			lastdist = (last1->x - last2->x) * (last1->x - last2->x)
					+ (last1->y - last2->y) * (last1->y - last2->y);

			if (thisdist < lastdist) {

				/*
				 * remove the "last" item from arrays 1 and 2.
				 * We move the "current" items up one position in the arrays,
				 * (into spaces [pos1 - 1] and [pos2 - 1]), and make
				 * them the new "last" items.
				 */
				remove_elem(star_array_1, pos1 - 1, num_stars_1);
				remove_elem(star_array_2, pos2 - 1, num_stars_2);
				last1 = &(star_array_1[pos1 - 1]);
				last2 = &(star_array_2[pos2 - 1]);
			} else {

				/*
				 * remove the current item from arrays 1 and 2.
				 * We can leave the "last" items as they are, since
				 * we haven't moved them.
				 */
				remove_elem(star_array_1, pos1, num_stars_1);
				remove_elem(star_array_2, pos2, num_stars_2);
			}
			pos1--;
			pos2--;
		} else {

			/* no repeated element.  Prepare for next step forward */
			last1 = s1;
			last2 = s2;
		}
		pos1++;
		pos2++;
	}
	return (SH_SUCCESS);
}

/*********************************************************************
 *
 * ROUTINE: remove_elem
 *
 * DESCRIPTION:
 * Remove the i'th element from the given array.
 *
 * What we do (slow as it is) is
 *
 *       1. move all elements after i up by 1
 *       2. subtract 1 from the number of elements in the array
 *
 * There's probably a better way of doing this, but let's
 * go with it for now.  1/19/96  MWR
 *
 * RETURN:
 *    nothing
 *
 */

static void remove_elem(s_star *star_array, /* I/O: we remove one element from this array */
int num, /* I: remove _this_ element */
int *num_stars /* I/O: on input: number of stars in array */
/*      on output: ditto, now smaller by one */
) {
	int i;
	s_star *s1, *s2;

	shAssert(star_array != NULL);
	shAssert(num < *num_stars);
	shAssert(num >= 0);
	shAssert(*num_stars > 0);

	s1 = &(star_array[num]);
	s2 = &(star_array[num + 1]);
	for (i = num; i < ((*num_stars) - 1); i++, s1++, s2++) {
		copy_star(s2, s1);
	}

	(*num_stars)--;
}

/*********************************************************************
 *
 * ROUTINE: remove_same_elements
 *
 * DESCRIPTION:
 * given two arrays of s_stars which have been sorted by their
 * "match_id" values, try to find s_stars which appear
 * in both arrays.  Remove any such s_stars from the second array.
 *
 * RETURN:
 *   nothing
 *
 */

static void remove_same_elements(s_star *star_array_1, /* I: look for elems which match those in array 2 */
int num_stars_1, /* I: number of elems in array 1 */
s_star *star_array_2, /* I/O: remove elems which match those in array 1 */
int *num_stars_2 /* I/O: number of elems in array 2 */
/*         will probably be smaller on output */
) {
	int pos1, pos2, pos2_top;
	s_star *s1, *s2;

	shAssert(star_array_1 != NULL);
	shAssert(star_array_2 != NULL);
	shAssert(num_stars_2 != NULL);

	pos1 = 0;
	pos2_top = 0;

	while (pos1 < num_stars_1) {

		s1 = &(star_array_1[pos1]);
		shAssert(s1 != NULL);

		for (pos2 = pos2_top; pos2 < *num_stars_2; pos2++) {
			s2 = &(star_array_2[pos2]);
			shAssert(s2 != NULL);

			if (s1->match_id == s2->match_id) {
				remove_elem(star_array_2, pos2, num_stars_2);
				if (--pos2_top < 0) {
					pos2_top = 0;
				}
			} else {
				if (s2->match_id < s1->match_id) {
					pos2_top = pos2 + 1;
				}
			}
		}
		pos1++;
	}
}

/***********************************************************************
 * ROUTINE: list_to_array
 *
 * DESCRIPTION:
 * Create an array of s_star structures, identical to the given linked
 * list.  Just make a copy of each structure.
 *
 * Sure, this is inefficient, but I'm using legacy code ...
 *
 * Return a pointer to the complete, filled array.
 *
 */

static s_star *
list_to_array(int num_stars, /* I: number of stars in the list */
struct s_star *list /* I: the linked list */
) {
	int i;
	struct s_star *array = NULL;
	struct s_star *ptr;
	struct s_star *star;

	/*
	 * okay, now we can walk down the CHAIN and create a new s_star
	 * for each item on the CHAIN.
	 */
	array = (s_star *) shMalloc(num_stars * sizeof(s_star));
	shAssert(array != NULL);
	for (i = 0, ptr = list; i < num_stars; i++, ptr = ptr->next) {
		shAssert(ptr != NULL);
		star = &(array[i]);
		shAssert(star != NULL);
		set_star(star, ptr->x, ptr->y, ptr->mag);
		star->match_id = i;
	}

	return (array);
}

/***********************************************************************
 * ROUTINE: write_array
 *
 * DESCRIPTION:
 * Given an array of s_star structures, write them to an ASCII text
 * file, with the following format:
 *
 *   ID    xvalue    yvalue    magvalue
 *
 * The 'ID' value is one assigned internally by these routines --
 * it doesn't correspond to any input ID value.
 *
 * RETURNS:
 *     nothing
 */

static void write_array(int num_stars, /* I: number of stars in the array */
struct s_star *star_array, /* I: the array of stars */
char *filename /* I: write into this file */
) {
	int i;
	FILE *fp;

	if ((fp = fopen(filename, "w")) == NULL) {
		shFatal("write_array: can't open file %s", filename);
	}

	for (i = 0; i < num_stars; i++) {
		fprintf(fp, "%6d %13.7f %13.7f %6.2f\n", star_array[i].id,
				star_array[i].x, star_array[i].y, star_array[i].mag);
	}

	fclose(fp);
}

/***********************************************************************
 * FUNCTION: reset_array_ids
 *
 * Modify the 'id' field values in the given array
 *   so that they will match the 'id' values in the corresponding
 *   stars of the given list.
 *
 * RETURNS
 *   nothing
 */

static void reset_array_ids(struct s_star *star_list, /* I: a list of stars */
int num_stars, /* I: number of stars in list and array */
struct s_star *star_array /* I/O: reset 'id' fields in this array */
) {
	int i;
	struct s_star *star_in_list, *star_in_array;

	star_in_list = star_list;
	for (i = 0; i < num_stars; i++) {

		star_in_array = &(star_array[i]);
		shAssert(star_in_list != NULL);
		shAssert(star_in_array != NULL);

		star_in_array->id = star_in_list->id;

		star_in_list = star_in_list->next;
	}
}

/************************************************************************
 *
 *
 * ROUTINE: calc_trans_linear
 *
 * DESCRIPTION:
 * Given a set of "nbright" matched pairs of stars, which we can
 * extract from the "winner_index" and "star_array" arrays,
 * figure out a TRANS structure which takes coordinates of
 * objects in set A and transforms then into coords for set B.
 *
 * In this case, a TRANS contains 6 coefficients in equations like this:
 *
 *       x' = A + B*x + C*y
 *       y' = D + E*x + F*y
 *
 * where (x,y) are coords in set A and (x',y') are corresponding
 * coords in set B.
 *
 * Internally, I'm going to solve for the very similar equations
 *
 *                x' = Ax + By + C
 *                y' = Dx + Ey + F
 *
 * and then just re-arrange the coefficients at the very end.  OK?
 *
 *
 * What we do is to treat each of the two equations above
 * separately.  We can write down 3 equations relating quantities
 * in the two sets of points (there are more than 3 such equations,
 * but we don't seek an exhaustive list).  For example,
 *
 *       a.       x'    =  Ax     + By    +  C
 *       b.       x'x   =  Ax^2   + Bxy   +  Cx      (mult both sides by x)
 *       c.       x'y   =  Axy    + By^2  +  Cy      (mult both sides by y)
 *
 * Now, since we have "nbright" matched pairs, we can take each of
 * the above 3 equations and form the sums on both sides, over
 * all "nbright" points.  So, if S(x) represents the sum of the quantity
 * "x" over all nbright points, and if we let N=nbright, then
 *
 *       a.     S(x')   =  AS(x)   + BS(y)   +  CN
 *       b.     S(x'x)  =  AS(x^2) + BS(xy)  +  CS(x)
 *       c.     S(x'y)  =  AS(xy)  + BS(y^2) +  CS(y)
 *
 * At this point, we have a set of three equations, and 3 unknowns: A, B, C.
 * We can write this set of equations as a matrix equation
 *
 *               b       = M * v
 *
 * where we KNOW the quantities
 *
 *        vector b = ( S(x'), S(x'x), S(x'y) )
 *
 *        matrix M = ( S(x)   S(y)    1      )
 *                   ( S(x^2) S(xy)   S(x)   )
 *                   ( S(xy)  S(y^2)  S(y)   )
 *
 *
 * and we want to FIND the unknown
 *
 *        vector v = ( A,     B,      C      )
 *
 * So, how to solve this matrix equation?  We use a Gaussian-elimination
 * method (see notes in 'gauss_matrix' function).   We solve
 * for A, B, C (and equivalently for D, E, F), then fill in the fields
 * of the given TRANS structure argument.
 *
 * It's possible that the matrix will be singular, and we can't find
 * a solution.  In that case, we print an error message and don't touch
 * the TRANS' fields.
 *
 *    [should explain how we make an iterative solution here,
 *     but will put in comments later.  MWR ]
 *
 * RETURN:
 *    SH_SUCCESS           if all goes well
 *    SH_GENERIC_ERROR     if we can't find a solution
 *
 * </AUTO>
 */

static int calc_trans_linear(int nbright, /* I: max number of stars we use in calculating */
/*      the transformation; we may cut down to */
/*      a more well-behaved subset. */
s_star *star_array_A, /* I: first array of s_star structure we match */
/*      the output TRANS takes their coords */
/*      into those of array B */
int num_stars_A, /* I: total number of stars in star_array_A */
s_star *star_array_B, /* I: second array of s_star structure we match */
int num_stars_B, /* I: total number of stars in star_array_B */
int *winner_votes, /* I: number of votes gotten by the top 'nbright' */
/*      matched pairs of stars */
int *winner_index_A, /* I: index into "star_array_A" of top */
/*      vote-getters */
int *winner_index_B, /* I: index into "star_array_B" of top */
/*      vote-getters */
TRANS *trans /* O: place solved coefficients into this */
/*      existing structure's fields */
) {
	int i;
	double **matrix;
	double vector[3];
	double solved_a, solved_b, solved_c, solved_d, solved_e, solved_f;
	s_star *s1, *s2;
	/* */
	double sum, sumx1, sumy1, sumx2, sumy2;
	double sumx1sq, sumy1sq;
	double sumx1y1, sumx1x2, sumx1y2;
	double sumy1x2, sumy1y2;

	shAssert(nbright >= AT_MATCH_REQUIRE_LINEAR);
	shAssert(trans->order == AT_TRANS_LINEAR);

	/*
	 * allocate a matrix we'll need for this function
	 */
	matrix = alloc_matrix(3);

	/*
	 * first, we consider the coefficients A, B, C in the trans.
	 * we form the sums that make up the elements of matrix M
	 */
	sum = 0.0;
	sumx1 = 0.0;
	sumy1 = 0.0;
	sumx2 = 0.0;
	sumy2 = 0.0;
	sumx1sq = 0.0;
	sumy1sq = 0.0;
	sumx1x2 = 0.0;
	sumx1y1 = 0.0;
	sumx1y2 = 0.0;
	sumy1x2 = 0.0;
	sumy1y2 = 0.0;

#ifdef DEBUG2
	printf(" in calc_trans_linear, num_A %5d num_B %5d here are %5d winners \n",
			num_stars_A, num_stars_B, nbright);
	for (i = 0; i < nbright; i++) {
		printf("   i %5d  winner_index_A %5d  winner_index_B %5d \n", i,
				winner_index_A[i], winner_index_B[i]);
	}
#endif

	for (i = 0; i < nbright; i++) {

		/* sanity checks */
		shAssert(winner_index_A[i] < num_stars_A);
		s1 = &(star_array_A[winner_index_A[i]]);
		shAssert(winner_index_B[i] < num_stars_B);
		s2 = &(star_array_B[winner_index_B[i]]);

		/* elements of the matrix */
		sum += 1.0;
		sumx1 += s1->x;
		sumx2 += s2->x;
		sumy1 += s1->y;
		sumy2 += s2->y;
		sumx1sq += s1->x * s1->x;
		sumy1sq += s1->y * s1->y;
		sumx1x2 += s1->x * s2->x;
		sumx1y1 += s1->x * s1->y;
		sumx1y2 += s1->x * s2->y;
		sumy1x2 += s1->y * s2->x;
		sumy1y2 += s1->y * s2->y;

	}

	/*
	 * now turn these sums into a matrix and a vector
	 */
	matrix[0][0] = sumx1sq;
	matrix[0][1] = sumx1y1;
	matrix[0][2] = sumx1;
	matrix[1][0] = sumx1y1;
	matrix[1][1] = sumy1sq;
	matrix[1][2] = sumy1;
	matrix[2][0] = sumx1;
	matrix[2][1] = sumy1;
	matrix[2][2] = sum;

	vector[0] = sumx1x2;
	vector[1] = sumy1x2;
	vector[2] = sumx2;

#ifdef DEBUG
	printf("before calling solution routines for ABC, here's matrix\n");
	print_matrix(matrix, 3);
#endif

	/*
	 * and now call the Gaussian-elimination routines to solve the matrix.
	 * The solution for TRANS coefficients A, B, C will be placed
	 * into the elements on "vector" after "gauss_matrix" finishes.
	 */
	if (gauss_matrix(matrix, 3, vector) != SH_SUCCESS) {
		shError("calc_trans_linear: can't solve for coeffs A,B,C ");
		free_matrix(matrix, 3);
		return (SH_GENERIC_ERROR);
	}

#ifdef DEBUG
	printf("after  calling solution routines, here's matrix\n");
	print_matrix(matrix, 3);
#endif

	solved_a = vector[0];
	solved_b = vector[1];
	solved_c = vector[2];

	/*
	 * Okay, now we solve for TRANS coefficients D, E, F, using the
	 * set of equations that relates y' to (x,y)
	 *
	 *       a.       y'    =  Dx     + Ey    +  F
	 *       b.       y'x   =  Dx^2   + Exy   +  Fx      (mult both sides by x)
	 *       c.       y'y   =  Dxy    + Ey^2  +  Fy      (mult both sides by y)
	 *
	 */
	matrix[0][0] = sumx1sq;
	matrix[0][1] = sumx1y1;
	matrix[0][2] = sumx1;
	matrix[1][0] = sumx1y1;
	matrix[1][1] = sumy1sq;
	matrix[1][2] = sumy1;
	matrix[2][0] = sumx1;
	matrix[2][1] = sumy1;
	matrix[2][2] = sum;

	vector[0] = sumx1y2;
	vector[1] = sumy1y2;
	vector[2] = sumy2;

#ifdef DEBUG
	printf("before calling solution routines for DEF, here's matrix\n");
	print_matrix(matrix, 3);
#endif

	/*
	 * and now call the Gaussian-elimination routines to solve the matrix.
	 * The solution for TRANS coefficients D, E, F will be placed
	 * into the elements on "vector" after "gauss_matrix" finishes.
	 */
	if (gauss_matrix(matrix, 3, vector) != SH_SUCCESS) {
		shError("calc_trans_linear: can't solve for coeffs D,E,F ");
		free_matrix(matrix, 3);
		return (SH_GENERIC_ERROR);
	}

#ifdef DEBUG
	printf("after  calling solution routines, here's matrix\n");
	print_matrix(matrix, 3);
#endif

	solved_d = vector[0];
	solved_e = vector[1];
	solved_f = vector[2];

	/*
	 * assign the coefficients we've just calculated to the output
	 * TRANS structure.  Recall that we've solved equations
	 *
	 *     x' = Ax + By + C
	 *     y' = Dx + Ey + F
	 *
	 * but that the TRANS structure assigns its coefficients assuming
	 *
	 *     x' = A + Bx + Cy
	 *     y' = D + Ex + Fy
	 *
	 * so, here, we have to re-arrange the coefficients a bit.
	 */
	trans->a = solved_c;
	trans->b = solved_a;
	trans->c = solved_b;
	trans->d = solved_f;
	trans->e = solved_d;
	trans->f = solved_e;

	/*
	 * and set the 'nr' field of the TRANS to the number of
	 * matching objects
	 */
	trans->nr = nbright;

	/*
	 * free up memory we allocated for this function
	 */
	free_matrix(matrix, 3);

	return (SH_SUCCESS);
}

/************************************************************************
 *
 *
 * ROUTINE: calc_trans_quadratic
 *
 * DESCRIPTION:
 * Given a set of "nbright" matched pairs of stars, which we can
 * extract from the "winner_index" and "star_array" arrays,
 * figure out a TRANS structure which takes coordinates of
 * objects in set A and transforms then into coords for set B.
 * In this case, a TRANS contains the twelve coefficients in the equations
 *
 *      x' =  A + Bx + Cy + Dxx + Exy + Fyy
 *      y' =  G + Hx + Iy + Jxx + Kxy + Lyy
 *
 * where (x,y) are coords in set A and (x',y') are corresponding
 * coords in set B.
 *
 *
 * What we do is to treat each of the two equations above
 * separately.  We can write down 6 equations relating quantities
 * in the two sets of points (there are more than 6 such equations,
 * but we don't seek an exhaustive list).  For example,
 *
 *  a.  x'    =  A    + Bx   + Cy    + Dxx   + Exy   +  Fyy
 *  b.  x'x   =  Ax   + Bxx  + Cxy   + Dxxx  + Exxy  +  Fxyy
 *  c.  x'y   =  Ay   + Bxy  + Cyy   + Dxxy  + Exyy  +  Fyyy
 *  d.  x'xx  =  Axx  + Bxxx + Cxxy  + Dxxxx + Exxxy +  Fxxyy
 *  e.  x'xy  =  Axy  + Bxxy + Cxyy  + Dxxxy + Exxyy +  Fxyyy
 *  f.  x'yy  =  Ayy  + Bxyy + Cyyy  + Dxxyy + Exyyy +  Fyyyy
 *
 * Now, since we have "nbright" matched pairs, we can take each of
 * the above 6 equations and form the sums on both sides, over
 * all "nbright" points.  So, if S(x) represents the sum of the quantity
 * "x" over all nbright points, and if we let N=nbright, then
 *
 *  a. S(x')   =  AN     + BS(x)   + CS(y)   + DS(xx)   + ES(xy)   +  FS(yy)
 *  b. S(x'x)  =  AS(x)  + BS(xx)  + CS(xy)  + DS(xxx)  + ES(xxy)  +  FS(xyy)
 *  c. S(x'y)  =  AS(y)  + BS(xy)  + CS(yy)  + DS(xxy)  + ES(xyy)  +  FS(yyy)
 *  d. S(x'xx) =  AS(xx) + BS(xxx) + CS(xxy) + DS(xxxx) + ES(xxxy) +  FS(xxyy)
 *  e. S(x'xy) =  AS(xy) + BS(xxy) + CS(xyy) + DS(xxxy) + ES(xxyy) +  FS(xyyy)
 *  f. S(x'yy) =  AS(yy) + BS(xyy) + CS(yyy) + DS(xxyy) + ES(xyyy) +  FS(yyyy)
 *
 * At this point, we have a set of 6 equations, and 6 unknowns:
 *        A, B, C, D, E, F
 *
 * We can write this set of equations as a matrix equation
 *
 *               b       = M * v
 *
 * where we KNOW the quantities
 *
 *        vector b = ( S(x'), S(x'x), S(x'y), S(x'xx), S(x'xy), S(x'yy) )
 *
 *        matrix M = [ N      S(x)    S(y)   S(xx)   S(xy)   S(yy)    ]
 *                   [ S(x)   S(xx)   S(xy)  S(xxx)  S(xxy)  S(xyy)   ]
 *                   [ S(y)   S(xy)   S(yy)  S(xxy)  S(xyy)  S(yyy)   ]
 *                   [ S(xx)  S(xxx)  S(xxy) S(xxxx) S(xxxy) S(xxyy)  ]
 *                   [ S(xy)  S(xxy)  S(xyy) S(xxxy) S(xxyy) S(xyyy)  ]
 *                   [ S(yy)  S(xyy)  S(yyy) S(xxyy) S(xyyy) S(yyyy)  ]
 *
 * and we want to FIND the unknown
 *
 *        vector v = ( A,     B,      C,     D,      E,      F )
 *
 * So, how to solve this matrix equation?  We use a Gaussian-elimination
 * method (see notes in 'gauss_matrix' function).   We solve
 * for A, B, C, D, E, F (and equivalently for G, H, I, J, K, L),
 * then fill in the fields
 * of the given TRANS structure argument.
 *
 * It's possible that the matrix will be singular, and we can't find
 * a solution.  In that case, we print an error message and don't touch
 * the TRANS' fields.
 *
 *    [should explain how we make an iterative solution here,
 *     but will put in comments later.  MWR ]
 *
 * RETURN:
 *    SH_SUCCESS           if all goes well
 *    SH_GENERIC_ERROR     if we can't find a solution
 *
 * </AUTO>
 */

static int calc_trans_quadratic(int nbright, /* I: max number of stars we use in calculating */
/*      the transformation; we may cut down to */
/*      a more well-behaved subset. */
s_star *star_array_A, /* I: first array of s_star structure we match */
/*      the output TRANS takes their coords */
/*      into those of array B */
int num_stars_A, /* I: total number of stars in star_array_A */
s_star *star_array_B, /* I: second array of s_star structure we match */
int num_stars_B, /* I: total number of stars in star_array_B */
int *winner_votes, /* I: number of votes gotten by the top 'nbright' */
/*      matched pairs of stars */
int *winner_index_A, /* I: index into "star_array_A" of top */
/*      vote-getters */
int *winner_index_B, /* I: index into "star_array_B" of top */
/*      vote-getters */
TRANS *trans /* O: place solved coefficients into this */
/*      existing structure's fields */
) {
	int i;
	double **matrix;
	double vector[6];
	double solved_a, solved_b, solved_c, solved_d, solved_e, solved_f;
	double solved_g, solved_h, solved_i, solved_j, solved_k, solved_l;
	s_star *s1, *s2;

	/*
	 * in variable names below, a '1' refers to coordinate of star s1
	 *   (which appear on both sides of the matrix equation)
	 *                      and a '2' refers to coordinate of star s2
	 *   (which appears only on left hand side of matrix equation)    o
	 */
	double sumx2, sumx2x1, sumx2y1, sumx2x1sq, sumx2x1y1, sumx2y1sq;
	double sumy2, sumy2x1, sumy2y1, sumy2x1sq, sumy2x1y1, sumy2y1sq;

	double sum, sumx1, sumy1, sumx1sq, sumx1y1, sumy1sq;
	double sumx1cu, sumx1sqy1, sumx1y1sq;
	double sumy1cu;
	double sumx1qu, sumx1cuy1, sumx1sqy1sq;
	double sumx1y1cu;
	double sumy1qu;

	shAssert(nbright >= AT_MATCH_REQUIRE_QUADRATIC);
	shAssert(trans->order == AT_TRANS_QUADRATIC);

	/*
	 * allocate a matrix we'll need for this function
	 */
	matrix = alloc_matrix(6);

	/*
	 * first, we consider the coefficients A, B, C, D, E, F in the trans.
	 * we form the sums that make up the elements of matrix M
	 */

	sum = 0.0;
	sumx1 = 0.0;
	sumy1 = 0.0;
	sumx1sq = 0.0;
	sumx1y1 = 0.0;
	sumy1sq = 0.0;
	sumx1cu = 0.0;
	sumx1sqy1 = 0.0;
	sumx1y1sq = 0.0;
	sumy1cu = 0.0;
	sumx1qu = 0.0;
	sumx1cuy1 = 0.0;
	sumx1sqy1sq = 0.0;
	sumx1y1cu = 0.0;
	sumy1qu = 0.0;

	sumx2 = 0.0;
	sumx2x1 = 0.0;
	sumx2y1 = 0.0;
	sumx2x1sq = 0.0;
	sumx2x1y1 = 0.0;
	sumx2y1sq = 0.0;
	sumy2 = 0.0;
	sumy2x1 = 0.0;
	sumy2y1 = 0.0;
	sumy2x1sq = 0.0;
	sumy2x1y1 = 0.0;
	sumy2y1sq = 0.0;

	for (i = 0; i < nbright; i++) {

		/* sanity checks */
		shAssert(winner_index_A[i] < num_stars_A);
		s1 = &(star_array_A[winner_index_A[i]]);
		shAssert(winner_index_B[i] < num_stars_B);
		s2 = &(star_array_B[winner_index_B[i]]);

		/* elements of the vectors */
		sumx2 += s2->x;
		sumx2x1 += s2->x * s1->x;
		sumx2y1 += s2->x * s1->y;
		sumx2x1sq += s2->x * s1->x * s1->x;
		sumx2x1y1 += s2->x * s1->x * s1->y;
		sumx2y1sq += s2->x * s1->y * s1->y;

		sumy2 += s2->y;
		sumy2x1 += s2->y * s1->x;
		sumy2y1 += s2->y * s1->y;
		sumy2x1sq += s2->y * s1->x * s1->x;
		sumy2x1y1 += s2->y * s1->x * s1->y;
		sumy2y1sq += s2->y * s1->y * s1->y;

		/* elements of the matrix */
		sum += 1.0;
		sumx1 += s1->x;
		sumy1 += s1->y;

		sumx1sq += s1->x * s1->x;
		sumx1y1 += s1->x * s1->y;
		sumy1sq += s1->y * s1->y;

		sumx1cu += s1->x * s1->x * s1->x;
		sumx1sqy1 += s1->x * s1->x * s1->y;
		sumx1y1sq += s1->x * s1->y * s1->y;
		sumy1cu += s1->y * s1->y * s1->y;

		sumx1qu += s1->x * s1->x * s1->x * s1->x;
		sumx1cuy1 += s1->x * s1->x * s1->x * s1->y;
		sumx1sqy1sq += s1->x * s1->x * s1->y * s1->y;
		sumx1y1cu += s1->x * s1->y * s1->y * s1->y;
		sumy1qu += s1->y * s1->y * s1->y * s1->y;

	}

	/*
	 * now turn these sums into a matrix and a vector
	 */
	matrix[0][0] = sum;
	matrix[0][1] = sumx1;
	matrix[0][2] = sumy1;
	matrix[0][3] = sumx1sq;
	matrix[0][4] = sumx1y1;
	matrix[0][5] = sumy1sq;

	matrix[1][0] = sumx1;
	matrix[1][1] = sumx1sq;
	matrix[1][2] = sumx1y1;
	matrix[1][3] = sumx1cu;
	matrix[1][4] = sumx1sqy1;
	matrix[1][5] = sumx1y1sq;

	matrix[2][0] = sumy1;
	matrix[2][1] = sumx1y1;
	matrix[2][2] = sumy1sq;
	matrix[2][3] = sumx1sqy1;
	matrix[2][4] = sumx1y1sq;
	matrix[2][5] = sumy1cu;

	matrix[3][0] = sumx1sq;
	matrix[3][1] = sumx1cu;
	matrix[3][2] = sumx1sqy1;
	matrix[3][3] = sumx1qu;
	matrix[3][4] = sumx1cuy1;
	matrix[3][5] = sumx1sqy1sq;

	matrix[4][0] = sumx1y1;
	matrix[4][1] = sumx1sqy1;
	matrix[4][2] = sumx1y1sq;
	matrix[4][3] = sumx1cuy1;
	matrix[4][4] = sumx1sqy1sq;
	matrix[4][5] = sumx1y1cu;

	matrix[5][0] = sumy1sq;
	matrix[5][1] = sumx1y1sq;
	matrix[5][2] = sumy1cu;
	matrix[5][3] = sumx1sqy1sq;
	matrix[5][4] = sumx1y1cu;
	matrix[5][5] = sumy1qu;

	vector[0] = sumx2;
	vector[1] = sumx2x1;
	vector[2] = sumx2y1;
	vector[3] = sumx2x1sq;
	vector[4] = sumx2x1y1;
	vector[5] = sumx2y1sq;

#ifdef DEBUG
	printf("before calling solution routines for ABCDEF, here's matrix\n");
	print_matrix(matrix, 6);
#endif

	/*
	 * and now call the Gaussian-elimination routines to solve the matrix.
	 * The solution for TRANS coefficients A, B, C, D, E, F will be placed
	 * into the elements on "vector" after "gauss_matrix" finishes.
	 */
	if (gauss_matrix(matrix, 6, vector) != SH_SUCCESS) {
		shError("calc_trans_quadratic: can't solve for coeffs A,B,C,D,E,F ");
		return (SH_GENERIC_ERROR);
	}

#ifdef DEBUG
	printf("after  calling solution routines, here's matrix\n");
	print_matrix(matrix, 6);
#endif

	solved_a = vector[0];
	solved_b = vector[1];
	solved_c = vector[2];
	solved_d = vector[3];
	solved_e = vector[4];
	solved_f = vector[5];

	/*
	 * Okay, now we solve for TRANS coefficients G, H, I, J, K, L, using the
	 * set of equations that relates y' to (x,y)
	 *
	 *      y'    =  G    + Hx   + Iy    + Jxx   + Kxy   +  Lyy
	 *      y'x   =  Gx   + Hxx  + Ixy   + Jxxx  + Kxxy  +  Lxyy
	 *      y'y   =  Gy   + Hxy  + Iyy   + Jxxy  + Kxyy  +  Lyyy
	 *      y'xx  =  Gxx  + Hxxx + Ixxy  + Jxxxx + Kxxxy +  Lxxyy
	 *      y'xy  =  Gxy  + Hxxy + Ixyy  + Jxxxy + Kxxyy +  Lxyyy
	 *      y'yy  =  Gyy  + Hxyy + Iyyy  + Jxxyy + Kxyyy +  Lyyyy
	 *
	 */
	matrix[0][0] = sum;
	matrix[0][1] = sumx1;
	matrix[0][2] = sumy1;
	matrix[0][3] = sumx1sq;
	matrix[0][4] = sumx1y1;
	matrix[0][5] = sumy1sq;

	matrix[1][0] = sumx1;
	matrix[1][1] = sumx1sq;
	matrix[1][2] = sumx1y1;
	matrix[1][3] = sumx1cu;
	matrix[1][4] = sumx1sqy1;
	matrix[1][5] = sumx1y1sq;

	matrix[2][0] = sumy1;
	matrix[2][1] = sumx1y1;
	matrix[2][2] = sumy1sq;
	matrix[2][3] = sumx1sqy1;
	matrix[2][4] = sumx1y1sq;
	matrix[2][5] = sumy1cu;

	matrix[3][0] = sumx1sq;
	matrix[3][1] = sumx1cu;
	matrix[3][2] = sumx1sqy1;
	matrix[3][3] = sumx1qu;
	matrix[3][4] = sumx1cuy1;
	matrix[3][5] = sumx1sqy1sq;

	matrix[4][0] = sumx1y1;
	matrix[4][1] = sumx1sqy1;
	matrix[4][2] = sumx1y1sq;
	matrix[4][3] = sumx1cuy1;
	matrix[4][4] = sumx1sqy1sq;
	matrix[4][5] = sumx1y1cu;

	matrix[5][0] = sumy1sq;
	matrix[5][1] = sumx1y1sq;
	matrix[5][2] = sumy1cu;
	matrix[5][3] = sumx1sqy1sq;
	matrix[5][4] = sumx1y1cu;
	matrix[5][5] = sumy1qu;

	vector[0] = sumy2;
	vector[1] = sumy2x1;
	vector[2] = sumy2y1;
	vector[3] = sumy2x1sq;
	vector[4] = sumy2x1y1;
	vector[5] = sumy2y1sq;

#ifdef DEBUG
	printf("before calling solution routines for GHIJKL, here's matrix\n");
	print_matrix(matrix, 6);
#endif

	/*
	 * and now call the Gaussian-elimination routines to solve the matrix.
	 * The solution for TRANS coefficients G, H, I, J, K, L will be placed
	 * into the elements on "vector" after "gauss_matrix" finishes.
	 */
	if (gauss_matrix(matrix, 6, vector) != SH_SUCCESS) {
		shError("calc_trans_quadratic: can't solve for coeffs G,H,I,J,K,L ");
		return (SH_GENERIC_ERROR);
	}

#ifdef DEBUG
	printf("after  calling solution routines, here's matrix\n");
	print_matrix(matrix, 6);
#endif

	solved_g = vector[0];
	solved_h = vector[1];
	solved_i = vector[2];
	solved_j = vector[3];
	solved_k = vector[4];
	solved_l = vector[5];

	/*
	 * assign the coefficients we've just calculated to the output
	 * TRANS structure.
	 */
	trans->a = solved_a;
	trans->b = solved_b;
	trans->c = solved_c;
	trans->d = solved_d;
	trans->e = solved_e;
	trans->f = solved_f;
	trans->g = solved_g;
	trans->h = solved_h;
	trans->i = solved_i;
	trans->j = solved_j;
	trans->k = solved_k;
	trans->l = solved_l;

	/*
	 * free up memory we allocated for this function
	 */
	free_matrix(matrix, 6);

	return (SH_SUCCESS);
}

/************************************************************************
 *
 *
 * ROUTINE: calc_trans_cubic
 *
 * DESCRIPTION:
 * Given a set of "nbright" matched pairs of stars, which we can
 * extract from the "winner_index" and "star_array" arrays,
 * figure out a TRANS structure which takes coordinates of
 * objects in set A and transforms then into coords for set B.
 * In this case, a TRANS contains the sixteen coefficients in the equations
 *
 *      x' =  A + Bx + Cy + Dxx + Exy + Fyy + Gx(xx+yy) + Hy(xx+yy)
 *      y' =  I + Jx + Ky + Lxx + Mxy + Nyy + Ox(xx+yy) + Py(xx+yy)
 *
 * where (x,y) are coords in set A and (x',y') are corresponding
 * coords in set B.
 *
 *
 * What we do is to treat each of the two equations above
 * separately.  We can write down 8 equations relating quantities
 * in the two sets of points (there are more than 8 such equations,
 * but we don't seek an exhaustive list).  For example,
 *
 *   x'    =  A    + Bx   + Cy    + Dxx   + Exy   +  Fyy   + GxR   + HyR
 *   x'x   =  Ax   + Bxx  + Cxy   + Dxxx  + Exxy  +  Fxyy  + GxxR  + HxyR
 *   x'y   =  Ay   + Bxy  + Cyy   + Dxxy  + Exyy  +  Fyyy  + GxyR  + HyyR
 *   x'xx  =  Axx  + Bxxx + Cxxy  + Dxxxx + Exxxy +  Fxxyy + GxxxR + HxxyR
 *   x'xy  =  Axy  + Bxxy + Cxyy  + Dxxxy + Exxyy +  Fxyyy + GxxyR + HxyyR
 *   x'yy  =  Ayy  + Bxyy + Cyyy  + Dxxyy + Exyyy +  Fyyyy + GxyyR + HyyyR
 *   x'xR  =  AxR  + BxxR + CxyR  + DxxxR + ExxyR +  FxyyR + GxxRR + HxyRR
 *   x'yR  =  AyR  + BxyR + CyyR  + DxxyR + ExyyR +  FyyyR + GxyRR + HyyRR
 *
 * (where we have used 'R' as an abbreviation for (xx + yy))
 *
 * Now, since we have "nbright" matched pairs, we can take each of
 * the above 8 equations and form the sums on both sides, over
 * all "nbright" points.  So, if S(x) represents the sum of the quantity
 * "x" over all nbright points, and if we let N=nbright, then
 *
 *  S(x')   =  AN     + BS(x)   + CS(y)   + DS(xx)   + ES(xy)   +  FS(yy)
 *                                                + GS(xR)   +  HS(yR)
 *  S(x'x)  =  AS(x)  + BS(xx)  + CS(xy)  + DS(xxx)  + ES(xxy)  +  FS(xyy)
 *                                                + GS(xxR)  +  HS(xyR)
 *  S(x'y)  =  AS(y)  + BS(xy)  + CS(yy)  + DS(xxy)  + ES(xyy)  +  FS(yyy)
 *                                                + GS(xyR)  +  HS(yyR)
 *  S(x'xx) =  AS(xx) + BS(xxx) + CS(xxy) + DS(xxxx) + ES(xxxy) +  FS(xxyy)
 *                                                + GS(xxxR) +  HS(xxyR)
 *  S(x'xy) =  AS(xy) + BS(xxy) + CS(xyy) + DS(xxxy) + ES(xxyy) +  FS(xyyy)
 *                                                + GS(xxyR) +  HS(xyyR)
 *  S(x'yy) =  AS(yy) + BS(xyy) + CS(yyy) + DS(xxyy) + ES(xyyy) +  FS(yyyy)
 *                                                + GS(xyyR) +  HS(yyyR)
 *  S(x'xR) =  AS(xR) + BS(xxR) + CS(xyR) + DS(xxxR) + ES(xxyR) +  FS(xyyR)
 *                                                + GS(xxRR) +  HS(xyRR)
 *  S(x'yR) =  AS(yR) + BS(xyR) + CS(yyR) + DS(xxyR) + ES(xyyR) +  FS(yyyR)
 *                                                + GS(xyRR) +  HS(yyRR)
 *
 * At this point, we have a set of 8 equations, and 8 unknowns:
 *        A, B, C, D, E, F, G, H
 *
 * We can write this set of equations as a matrix equation
 *
 *               b       = M * v
 *
 * where we KNOW the quantities
 *
 *  b = ( S(x'), S(x'x), S(x'y), S(x'xx), S(x'xy), S(x'yy), S(x'xR), S(x'rR) )
 *
 * matr M = [ N      S(x)    S(y)   S(xx)   S(xy)   S(yy)   S(xR)   S(yR)   ]
 *          [ S(x)   S(xx)   S(xy)  S(xxx)  S(xxy)  S(xyy)  S(xxR)  S(xyR)  ]
 *          [ S(y)   S(xy)   S(yy)  S(xxy)  S(xyy)  S(yyy)  S(xyR)  S(yyR)  ]
 *          [ S(xx)  S(xxx)  S(xxy) S(xxxx) S(xxxy) S(xxyy) S(xxxR) S(xxyR) ]
 *          [ S(xy)  S(xxy)  S(xyy) S(xxxy) S(xxyy) S(xyyy) S(xxyR) S(xyyR) ]
 *          [ S(yy)  S(xyy)  S(yyy) S(xxyy) S(xyyy) S(yyyy) S(xyyR) S(yyyR) ]
 *          [ S(xR)  S(xxR)  S(xyR) S(xxxR) S(xxyR) S(xyyR) S(xxRR) S(xyRR) ]
 *          [ S(yR)  S(xyR)  S(yyR) S(xxyR) S(xyyR) S(yyyR) S(xyRR) S(yyRR) ]
 *
 * and we want to FIND the unknown
 *
 *        vector v = ( A,     B,      C,     D,      E,      F,     G,     H )
 *
 * So, how to solve this matrix equation?  We use a Gaussian-elimination
 * method (see notes in 'gauss_matrix' function).   We solve
 * for A, B, C, D, E, F, G, H (and equivalently for I, J, K, L, M, N, O, P),
 * then fill in the fields
 * of the given TRANS structure argument.
 *
 * It's possible that the matrix will be singular, and we can't find
 * a solution.  In that case, we print an error message and don't touch
 * the TRANS' fields.
 *
 *    [should explain how we make an iterative solution here,
 *     but will put in comments later.  MWR ]
 *
 * RETURN:
 *    SH_SUCCESS           if all goes well
 *    SH_GENERIC_ERROR     if we can't find a solution
 *
 * </AUTO>
 */

static int calc_trans_cubic(int nbright, /* I: max number of stars we use in calculating */
/*      the transformation; we may cut down to */
/*      a more well-behaved subset. */
s_star *star_array_A, /* I: first array of s_star structure we match */
/*      the output TRANS takes their coords */
/*      into those of array B */
int num_stars_A, /* I: total number of stars in star_array_A */
s_star *star_array_B, /* I: second array of s_star structure we match */
int num_stars_B, /* I: total number of stars in star_array_B */
int *winner_votes, /* I: number of votes gotten by the top 'nbright' */
/*      matched pairs of stars */
int *winner_index_A, /* I: index into "star_array_A" of top */
/*      vote-getters */
int *winner_index_B, /* I: index into "star_array_B" of top */
/*      vote-getters */
TRANS *trans /* O: place solved coefficients into this */
/*      existing structure's fields */
) {
	int i;
	double **matrix;
	double vector[8];
	double solved_a, solved_b, solved_c, solved_d, solved_e, solved_f;
	double solved_g, solved_h;
	double solved_i, solved_j, solved_k, solved_l, solved_m, solved_n;
	double solved_o, solved_p;
	s_star *s1, *s2;

	/*
	 * the variable 'R' will hold the value (x1*x1 + y1*y1);
	 *   in other words, the square of the distance of (x1, y1)
	 *   from the origin.
	 */
	double R;
	/*
	 * in variable names below, a '1' refers to coordinate of star s1
	 *   (which appear on both sides of the matrix equation)
	 *                      and a '2' refers to coordinate of star s2
	 *   (which appears only on left hand side of matrix equation)    o
	 */
	double sumx2, sumx2x1, sumx2y1, sumx2x1sq, sumx2x1y1, sumx2y1sq;
	double sumx2x1R, sumx2y1R;
	double sumy2, sumy2x1, sumy2y1, sumy2x1sq, sumy2x1y1, sumy2y1sq;
	double sumy2x1R, sumy2y1R;

	double sum, sumx1, sumy1, sumx1sq, sumx1y1, sumy1sq;
	double sumx1cu, sumx1sqy1, sumx1y1sq;
	double sumy1cu;
	double sumx1R, sumy1R, sumx1sqR, sumx1y1R, sumy1sqR;
	double sumx1cuR, sumx1sqy1R, sumx1y1sqR, sumy1cuR;
	double sumx1qu, sumx1cuy1, sumx1sqy1sq;
	double sumx1y1cu;
	double sumy1qu;
	double sumx1sqRsq, sumx1y1Rsq, sumy1sqRsq;

	shAssert(nbright >= AT_MATCH_REQUIRE_CUBIC);
	shAssert(trans->order == AT_TRANS_CUBIC);

	/*
	 * allocate a matrix we'll need for this function
	 */
	matrix = alloc_matrix(8);

	/*
	 * first, we consider the coefficients A, B, C, D, E, F, G, H in the trans.
	 * we form the sums that make up the elements of matrix M
	 */

	sum = 0.0;
	sumx1 = 0.0;
	sumy1 = 0.0;
	sumx1sq = 0.0;
	sumx1y1 = 0.0;
	sumy1sq = 0.0;
	sumx1cu = 0.0;
	sumx1sqy1 = 0.0;
	sumx1y1sq = 0.0;
	sumy1cu = 0.0;
	sumx1qu = 0.0;
	sumx1cuy1 = 0.0;
	sumx1sqy1sq = 0.0;
	sumx1y1cu = 0.0;
	sumy1qu = 0.0;
	sumx1R = 0.0;
	sumy1R = 0.0;
	sumx1sqR = 0.0;
	sumx1y1R = 0.0;
	sumy1sqR = 0.0;
	sumx1cuR = 0.0;
	sumx1sqy1R = 0.0;
	sumx1y1sqR = 0.0;
	sumy1cuR = 0.0;
	sumx1sqRsq = 0.0;
	sumx1y1Rsq = 0.0;
	sumy1sqRsq = 0.0;

	sumx2 = 0.0;
	sumx2x1 = 0.0;
	sumx2y1 = 0.0;
	sumx2x1sq = 0.0;
	sumx2x1y1 = 0.0;
	sumx2y1sq = 0.0;
	sumx2x1R = 0.0;
	sumx2y1R = 0.0;
	sumy2 = 0.0;
	sumy2x1 = 0.0;
	sumy2y1 = 0.0;
	sumy2x1sq = 0.0;
	sumy2x1y1 = 0.0;
	sumy2y1sq = 0.0;
	sumy2x1R = 0.0;
	sumy2y1R = 0.0;

	for (i = 0; i < nbright; i++) {

		/* sanity checks */
		shAssert(winner_index_A[i] < num_stars_A);
		s1 = &(star_array_A[winner_index_A[i]]);
		shAssert(winner_index_B[i] < num_stars_B);
		s2 = &(star_array_B[winner_index_B[i]]);

		/* elements of the vectors */
		R = (s1->x * s1->x + s1->y * s1->y);

		sumx2 += s2->x;
		sumx2x1 += s2->x * s1->x;
		sumx2y1 += s2->x * s1->y;
		sumx2x1sq += s2->x * s1->x * s1->x;
		sumx2x1y1 += s2->x * s1->x * s1->y;
		sumx2y1sq += s2->x * s1->y * s1->y;
		sumx2x1R += s2->x * s1->x * R;
		sumx2y1R += s2->x * s1->y * R;

		sumy2 += s2->y;
		sumy2x1 += s2->y * s1->x;
		sumy2y1 += s2->y * s1->y;
		sumy2x1sq += s2->y * s1->x * s1->x;
		sumy2x1y1 += s2->y * s1->x * s1->y;
		sumy2y1sq += s2->y * s1->y * s1->y;
		sumy2x1R += s2->y * s1->x * R;
		sumy2y1R += s2->y * s1->y * R;

		/* elements of the matrix */
		sum += 1.0;
		sumx1 += s1->x;
		sumy1 += s1->y;

		sumx1sq += s1->x * s1->x;
		sumx1y1 += s1->x * s1->y;
		sumy1sq += s1->y * s1->y;

		sumx1cu += s1->x * s1->x * s1->x;
		sumx1sqy1 += s1->x * s1->x * s1->y;
		sumx1y1sq += s1->x * s1->y * s1->y;
		sumy1cu += s1->y * s1->y * s1->y;

		sumx1qu += s1->x * s1->x * s1->x * s1->x;
		sumx1cuy1 += s1->x * s1->x * s1->x * s1->y;
		sumx1sqy1sq += s1->x * s1->x * s1->y * s1->y;
		sumx1y1cu += s1->x * s1->y * s1->y * s1->y;
		sumy1qu += s1->y * s1->y * s1->y * s1->y;

		sumx1R += s1->x * R;
		sumy1R += s1->y * R;
		sumx1sqR += s1->x * s1->x * R;
		sumx1y1R += s1->x * s1->y * R;
		sumy1sqR += s1->y * s1->y * R;

		sumx1cuR += s1->x * s1->x * s1->x * R;
		sumx1sqy1R += s1->x * s1->x * s1->y * R;
		sumx1y1sqR += s1->x * s1->y * s1->y * R;
		sumy1cuR += s1->y * s1->y * s1->y * R;

		sumx1sqRsq += s1->x * s1->x * R * R;
		sumx1y1Rsq += s1->x * s1->y * R * R;
		sumy1sqRsq += s1->y * s1->y * R * R;

	}

	/*
	 * now turn these sums into a matrix and a vector
	 */
	matrix[0][0] = sum;
	matrix[0][1] = sumx1;
	matrix[0][2] = sumy1;
	matrix[0][3] = sumx1sq;
	matrix[0][4] = sumx1y1;
	matrix[0][5] = sumy1sq;
	matrix[0][6] = sumx1R;
	matrix[0][7] = sumy1R;

	matrix[1][0] = sumx1;
	matrix[1][1] = sumx1sq;
	matrix[1][2] = sumx1y1;
	matrix[1][3] = sumx1cu;
	matrix[1][4] = sumx1sqy1;
	matrix[1][5] = sumx1y1sq;
	matrix[1][6] = sumx1sqR;
	matrix[1][7] = sumx1y1R;

	matrix[2][0] = sumy1;
	matrix[2][1] = sumx1y1;
	matrix[2][2] = sumy1sq;
	matrix[2][3] = sumx1sqy1;
	matrix[2][4] = sumx1y1sq;
	matrix[2][5] = sumy1cu;
	matrix[2][6] = sumx1y1R;
	matrix[2][7] = sumy1sqR;

	matrix[3][0] = sumx1sq;
	matrix[3][1] = sumx1cu;
	matrix[3][2] = sumx1sqy1;
	matrix[3][3] = sumx1qu;
	matrix[3][4] = sumx1cuy1;
	matrix[3][5] = sumx1sqy1sq;
	matrix[3][6] = sumx1cuR;
	matrix[3][7] = sumx1sqy1R;

	matrix[4][0] = sumx1y1;
	matrix[4][1] = sumx1sqy1;
	matrix[4][2] = sumx1y1sq;
	matrix[4][3] = sumx1cuy1;
	matrix[4][4] = sumx1sqy1sq;
	matrix[4][5] = sumx1y1cu;
	matrix[4][6] = sumx1sqy1R;
	matrix[4][7] = sumx1y1sqR;

	matrix[5][0] = sumy1sq;
	matrix[5][1] = sumx1y1sq;
	matrix[5][2] = sumy1cu;
	matrix[5][3] = sumx1sqy1sq;
	matrix[5][4] = sumx1y1cu;
	matrix[5][5] = sumy1qu;
	matrix[5][6] = sumx1y1sqR;
	matrix[5][7] = sumy1cuR;

	matrix[6][0] = sumx1R;
	matrix[6][1] = sumx1sqR;
	matrix[6][2] = sumx1y1R;
	matrix[6][3] = sumx1cuR;
	matrix[6][4] = sumx1sqy1R;
	matrix[6][5] = sumx1y1sqR;
	matrix[6][6] = sumx1sqRsq;
	matrix[6][7] = sumx1y1Rsq;

	matrix[7][0] = sumy1R;
	matrix[7][1] = sumx1y1R;
	matrix[7][2] = sumy1sqR;
	matrix[7][3] = sumx1sqy1R;
	matrix[7][4] = sumx1y1sqR;
	matrix[7][5] = sumy1cuR;
	matrix[7][6] = sumx1y1Rsq;
	matrix[7][7] = sumy1sqRsq;

	vector[0] = sumx2;
	vector[1] = sumx2x1;
	vector[2] = sumx2y1;
	vector[3] = sumx2x1sq;
	vector[4] = sumx2x1y1;
	vector[5] = sumx2y1sq;
	vector[6] = sumx2x1R;
	vector[7] = sumx2y1R;

#ifdef DEBUG
	printf("before calling solution routines for ABCDEFGH, here's matrix\n");
	print_matrix(matrix, 8);
#endif

	/*
	 * and now call the Gaussian-elimination routines to solve the matrix.
	 * The solution for TRANS coefficients A, B, C, D, E, F will be placed
	 * into the elements on "vector" after "gauss_matrix" finishes.
	 */
	if (gauss_matrix(matrix, 8, vector) != SH_SUCCESS) {
		shError("calc_trans_cubic: can't solve for coeffs A,B,C,D,E,F,G,H ");
		return (SH_GENERIC_ERROR);
	}

#ifdef DEBUG
	printf("after  calling solution routines, here's matrix\n");
	print_matrix(matrix, 8);
#endif

	solved_a = vector[0];
	solved_b = vector[1];
	solved_c = vector[2];
	solved_d = vector[3];
	solved_e = vector[4];
	solved_f = vector[5];
	solved_g = vector[6];
	solved_h = vector[7];

	/*
	 * Okay, now we solve for TRANS coefficients I, J, K, L, M, N, O, P
	 * using the * set of equations that relates y' to (x,y)
	 *
	 *  y'    =  I    + Jx   + Ky    + Lxx   + Mxy   +  Nyy   + OxR   + PyR
	 *  y'x   =  Ix   + Jxx  + Kxy   + Lxxx  + Mxxy  +  Nxyy  + OxxR  + PxyR
	 *  y'y   =  Iy   + Jxy  + Kyy   + Lxxy  + Mxyy  +  Nyyy  + OxyR  + PyyR
	 *  y'xx  =  Ixx  + Jxxx + Kxxy  + Lxxxx + Mxxxy +  Nxxyy + OxxxR + PxxyR
	 *  y'xy  =  Ixy  + Jxxy + Kxyy  + Lxxxy + Mxxyy +  Nxyyy + OxxyR + PxyyR
	 *  y'yy  =  Iyy  + Jxyy + Kyyy  + Lxxyy + Mxyyy +  Nyyyy + OxyyR + PyyyR
	 *  y'xR  =  IxR  + JxxR + KxyR  + LxxxR + MxxyR +  NxyyR + OxxRR + PxyRR
	 *  y'yR  =  IyR  + JxyR + KyyR  + LxxyR + MxyyR +  NyyyR + OxyRR + PyyRR
	 *
	 */
	matrix[0][0] = sum;
	matrix[0][1] = sumx1;
	matrix[0][2] = sumy1;
	matrix[0][3] = sumx1sq;
	matrix[0][4] = sumx1y1;
	matrix[0][5] = sumy1sq;
	matrix[0][6] = sumx1R;
	matrix[0][7] = sumy1R;

	matrix[1][0] = sumx1;
	matrix[1][1] = sumx1sq;
	matrix[1][2] = sumx1y1;
	matrix[1][3] = sumx1cu;
	matrix[1][4] = sumx1sqy1;
	matrix[1][5] = sumx1y1sq;
	matrix[1][6] = sumx1sqR;
	matrix[1][7] = sumx1y1R;

	matrix[2][0] = sumy1;
	matrix[2][1] = sumx1y1;
	matrix[2][2] = sumy1sq;
	matrix[2][3] = sumx1sqy1;
	matrix[2][4] = sumx1y1sq;
	matrix[2][5] = sumy1cu;
	matrix[2][6] = sumx1y1R;
	matrix[2][7] = sumy1sqR;

	matrix[3][0] = sumx1sq;
	matrix[3][1] = sumx1cu;
	matrix[3][2] = sumx1sqy1;
	matrix[3][3] = sumx1qu;
	matrix[3][4] = sumx1cuy1;
	matrix[3][5] = sumx1sqy1sq;
	matrix[3][6] = sumx1cuR;
	matrix[3][7] = sumx1sqy1R;

	matrix[4][0] = sumx1y1;
	matrix[4][1] = sumx1sqy1;
	matrix[4][2] = sumx1y1sq;
	matrix[4][3] = sumx1cuy1;
	matrix[4][4] = sumx1sqy1sq;
	matrix[4][5] = sumx1y1cu;
	matrix[4][6] = sumx1sqy1R;
	matrix[4][7] = sumx1y1sqR;

	matrix[5][0] = sumy1sq;
	matrix[5][1] = sumx1y1sq;
	matrix[5][2] = sumy1cu;
	matrix[5][3] = sumx1sqy1sq;
	matrix[5][4] = sumx1y1cu;
	matrix[5][5] = sumy1qu;
	matrix[5][6] = sumx1y1sqR;
	matrix[5][7] = sumy1cuR;

	matrix[6][0] = sumx1R;
	matrix[6][1] = sumx1sqR;
	matrix[6][2] = sumx1y1R;
	matrix[6][3] = sumx1cuR;
	matrix[6][4] = sumx1sqy1R;
	matrix[6][5] = sumx1y1sqR;
	matrix[6][6] = sumx1sqRsq;
	matrix[6][7] = sumx1y1Rsq;

	matrix[7][0] = sumy1R;
	matrix[7][1] = sumx1y1R;
	matrix[7][2] = sumy1sqR;
	matrix[7][3] = sumx1sqy1R;
	matrix[7][4] = sumx1y1sqR;
	matrix[7][5] = sumy1cuR;
	matrix[7][6] = sumx1y1Rsq;
	matrix[7][7] = sumy1sqRsq;

	vector[0] = sumy2;
	vector[1] = sumy2x1;
	vector[2] = sumy2y1;
	vector[3] = sumy2x1sq;
	vector[4] = sumy2x1y1;
	vector[5] = sumy2y1sq;
	vector[6] = sumy2x1R;
	vector[7] = sumy2y1R;

#ifdef DEBUG
	printf("before calling solution routines for IJKLMNOP, here's matrix\n");
	print_matrix(matrix, 8);
#endif

	/*
	 * and now call the Gaussian-elimination routines to solve the matrix.
	 * The solution for TRANS coefficients I, J, K, L, M, N, O, P will be placed
	 * into the elements on "vector" after "gauss_matrix" finishes.
	 */
	if (gauss_matrix(matrix, 8, vector) != SH_SUCCESS) {
		shError("calc_trans_cubic: can't solve for coeffs I,J,K,L,M,N,O,P ");
		return (SH_GENERIC_ERROR);
	}

#ifdef DEBUG
	printf("after  calling solution routines, here's matrix\n");
	print_matrix(matrix, 8);
#endif

	solved_i = vector[0];
	solved_j = vector[1];
	solved_k = vector[2];
	solved_l = vector[3];
	solved_m = vector[4];
	solved_n = vector[5];
	solved_o = vector[6];
	solved_p = vector[7];

	/*
	 * assign the coefficients we've just calculated to the output
	 * TRANS structure.
	 */
	trans->a = solved_a;
	trans->b = solved_b;
	trans->c = solved_c;
	trans->d = solved_d;
	trans->e = solved_e;
	trans->f = solved_f;
	trans->g = solved_g;
	trans->h = solved_h;
	trans->i = solved_i;
	trans->j = solved_j;
	trans->k = solved_k;
	trans->l = solved_l;
	trans->m = solved_m;
	trans->n = solved_n;
	trans->o = solved_o;
	trans->p = solved_p;

	/*
	 * free up memory we allocated for this function
	 */
	free_matrix(matrix, 8);

	return (SH_SUCCESS);
}

/***************************************************************************
 * PROCEDURE: gauss_matrix
 *
 * DESCRIPTION:
 * Given a square 2-D 'num'-by-'num' matrix, called "matrix", and given
 * a 1-D vector "vector" of 'num' elements, find the 1-D vector
 * called "solution_vector" which satisfies the equation
 *
 *      matrix * solution_vector  =  vector
 *
 * where the * above represents matrix multiplication.
 *
 * What we do is to use Gaussian elimination (with partial pivoting)
 * and back-substitution to find the solution_vector.
 * We do not pivot in place, but physically move values -- it
 * doesn't take much time in this application.  After we have found the
 * "solution_vector", we replace the contents of "vector" with the
 * "solution_vector".
 *
 * This is a common algorithm.  See any book on linear algebra or
 * numerical solutions; for example, "Numerical Methods for Engineers,"
 * by Steven C. Chapra and Raymond P. Canale, McGraw-Hill, 1998,
 * Chapter 9.
 *
 * If an error occurs (if the matrix is singular), this prints an error
 * message and returns with error code.
 *
 * RETURN:
 *    SH_SUCCESS          if all goes well
 *    SH_GENERIC_ERROR    if not -- if matrix is singular
 *
 * </AUTO>
 */

static int gauss_matrix(double **matrix, /* I/O: the square 2-D matrix we'll invert */
/*      will hold inverse matrix on output */
int num, /* I: number of rows and cols in matrix */
double *vector /* I/O: vector which holds "b" values in input */
/*      and the solution vector "x" on output */
) {
	int i, j, k;
	double *biggest_val;
	double *solution_vector;
	double factor;
	double sum;

#ifdef DEBUG
	print_matrix(matrix, num);
#endif

	biggest_val = (double *) shMalloc(num * sizeof(double));
	solution_vector = (double *) shMalloc(num * sizeof(double));

	/*
	 * step 1: we find the largest value in each row of matrix,
	 *         and store those values in 'biggest_val' array.
	 *         We use this information to pivot the matrix.
	 */
	for (i = 0; i < num; i++) {
		biggest_val[i] = fabs(matrix[i][0]);
		for (j = 1; j < num; j++) {
			if (fabs(matrix[i][j]) > biggest_val[i]) {
				biggest_val[i] = fabs(matrix[i][j]);
			}
		}
		if (biggest_val[i] == 0.0) {
			shError("gauss_matrix: biggest val in row %d is zero", i);
			shFree(biggest_val);
			shFree(solution_vector);
			return (SH_GENERIC_ERROR);
		}
	}

	/*
	 * step 2: we use Gaussian elimination to convert the "matrix"
	 *         into a triangular matrix, in which the values of all
	 *         elements below the diagonal is zero.
	 */
	for (i = 0; i < num - 1; i++) {

		/* pivot this row (if necessary) */
		if (gauss_pivot(matrix, num, vector, biggest_val, i)
				== SH_GENERIC_ERROR) {
			shError("gauss_matrix: singular matrix");
			shFree(biggest_val);
			shFree(solution_vector);
			return (SH_GENERIC_ERROR);
		}

		if (fabs(matrix[i][i] / biggest_val[i]) < MATRIX_TOL) {
			shError("gauss_matrix: Y: row %d has tiny value %f / %f", i,
					matrix[i][i], biggest_val[i]);
			shFree(biggest_val);
			shFree(solution_vector);
			return (SH_GENERIC_ERROR);
		}

		/* we eliminate this variable in all rows below the current one */
		for (j = i + 1; j < num; j++) {
			factor = matrix[j][i] / matrix[i][i];
			for (k = i + 1; k < num; k++) {
				matrix[j][k] -= factor * matrix[i][k];
			}
			/* and in the vector, too */
			vector[j] -= factor * vector[i];
		}

	}

	/*
	 * make sure that the last row's single remaining element
	 * isn't too tiny
	 */
	if (fabs(matrix[num - 1][num - 1] / biggest_val[num - 1]) < MATRIX_TOL) {
		shError("gauss_matrix: Z: row %d has tiny value %f / %f", num,
				matrix[num - 1][num - 1], biggest_val[num - 1]);
		shFree(biggest_val);
		shFree(solution_vector);
		return (SH_GENERIC_ERROR);
	}

	/*
	 * step 3: we can now calculate the solution_vector values
	 *         via back-substitution; we start at the last value in the
	 *         vector (at the "bottom" of the vector) and work
	 *         upwards towards the top.
	 */
	solution_vector[num - 1] = vector[num - 1] / matrix[num - 1][num - 1];
	for (i = num - 2; i >= 0; i--) {
		sum = 0.0;
		for (j = i + 1; j < num; j++) {
			sum += matrix[i][j] * solution_vector[j];
		}
		solution_vector[i] = (vector[i] - sum) / matrix[i][i];
	}

	/*
	 * step 4: okay, we've found the values in the solution vector!
	 *         We now replace the input values in 'vector' with these
	 *         solution_vector values, and we're done.
	 */
	for (i = 0; i < num; i++) {
		vector[i] = solution_vector[i];
	}

	/* clean up */
	shFree(solution_vector);
	shFree(biggest_val);

	return (SH_SUCCESS);
}

/***************************************************************************
 * PROCEDURE: gauss_pivot
 *
 * DESCRIPTION:
 * This routine is called by "gauss_matrix".  Given a square "matrix"
 * of "num"-by-"num" elements, and given a "vector" of "num" elements,
 * and given a particular "row" value, this routine finds the largest
 * value in the matrix at/below the given "row" position.  If that
 * largest value isn't in the given "row", this routine switches
 * rows in the matrix (and in the vector) so that the largest value
 * will now be in "row".
 *
 * RETURN:
 *    SH_SUCCESS          if all goes well
 *    SH_GENERIC_ERROR    if not -- if matrix is singular
 *
 * </AUTO>
 */

#define SWAP(a,b)  { double temp = (a); (a) = (b); (b) = temp; }

static int gauss_pivot(double **matrix, /* I/O: a square 2-D matrix we are inverting */
int num, /* I: number of rows and cols in matrix */
double *vector, /* I/O: vector which holds "b" values in input */
double *biggest_val, /* I: largest value in each row of matrix */
int row /* I: want to pivot around this row */
) {
	int i;
	int col, pivot_row;
	double big, other_big;

	/* sanity checks */
	shAssert(matrix != NULL);
	shAssert(vector != NULL);
	shAssert(row < num);

	pivot_row = row;
	big = fabs(matrix[row][row] / biggest_val[row]);
#ifdef DEBUG
	print_matrix(matrix, num);
	printf(" biggest_val:  ");
	for (i = 0; i < num; i++) {
		printf("%9.4e ", biggest_val[i]);
	}
	printf("\n");
	printf("  gauss_pivot: row %3d  %9.4e %9.4e %12.5e ", row, matrix[row][row],
			biggest_val[row], big);
#endif

	for (i = row + 1; i < num; i++) {
		other_big = fabs(matrix[i][row] / biggest_val[i]);
		if (other_big > big) {
			big = other_big;
			pivot_row = i;
		}
	}
#ifdef DEBUG
	printf("  pivot_row %3d  %9.4e %9.4e %12.5e ", pivot_row,
			matrix[pivot_row][pivot_row], biggest_val[pivot_row], big);
#endif

	/*
	 * if another row is better for pivoting, switch it with 'row'
	 *    and switch the corresponding elements in 'vector'
	 *    and switch the corresponding elements in 'biggest_val'
	 */
	if (pivot_row != row) {
#ifdef DEBUG
		printf("   will swap \n");
#endif
		for (col = row; col < num; col++) {
			SWAP(matrix[pivot_row][col], matrix[row][col]);
		}
		SWAP(vector[pivot_row], vector[row]);
		SWAP(biggest_val[pivot_row], biggest_val[row]);
	} else {
#ifdef DEBUG
		printf("    no swap \n");
#endif
	}

	return (SH_SUCCESS);
}

/***************************************************************************
 * PROCEDURE: atFindMedtf
 *
 * DESCRIPTION:
 * Assume that the two input lists of stars were taken by the same
 * instrument, so that they have the same scale and rotation.
 * If this is true, then a simple translation ought to register
 * the two lists.  This routine tries to characterize that
 * translation: it finds both the mean shift in (x, y),
 * and the median shift in (x, y).  It also calculates the
 * (clipped) standard deviation from the mean shift.  The results of all
 * these calculations are placed into the given MEDTF structure.
 *
 * RETURN:
 *    SH_SUCCESS          if all goes well
 *    SH_GENERIC_ERROR    if not
 *
 * </AUTO>
 */

int atFindMedtf(int num_matched_A, /* I: number of matched stars in list A */
struct s_star *listA, /* I: list of matched stars from set A */
int num_matched_B, /* I: number of matched stars in list B */
struct s_star *listB, /* I: list of matched stars from set B */
double medsigclip, /* I: sigma-clipping factor */
MEDTF *medtf /* O: we place results into this structure */
) {
	int i, nstar, num_within_clip;
	double *dx, *dy, *dxclip, *dyclip;
	double xdist, ydist;
	double Dx_sum, Dx_sum2, Dx_rms, Dx_ave, Dx_med;
	double Dy_sum, Dy_sum2, Dy_rms, Dy_ave, Dy_med, clip;
	struct s_star *A_star, *B_star;

	if (num_matched_A < 3) {
		shError("atFindMedtf: fewer than 3 matched pairs; cannot find MEDTF");
		return (SH_GENERIC_ERROR);
	}
	nstar = num_matched_A;

	/* sanity checks */
	shAssert(num_matched_A == num_matched_B);
	shAssert(listA != NULL);
	shAssert(listB != NULL);
	shAssert(medsigclip >= 0);
	shAssert(medtf != NULL);

	/*
	 * allocate space for arrays we use to sort distances
	 * between matched items in the lists
	 */
	dx = (double *) shMalloc(nstar * sizeof(double));
	dy = (double *) shMalloc(nstar * sizeof(double));

	/*
	 * Step 1: calculate distances between matched stars,
	 *         and fill the "dx" and "dy" arrays for later calculation
	 *         of the median
	 */
	Dx_sum = 0.0;
	Dy_sum = 0.0;
	Dx_sum2 = 0.0;
	Dy_sum2 = 0.0;
	A_star = listA;
	B_star = listB;
	for (i = 0; i < nstar; i++) {
		shAssert(A_star != NULL);
		shAssert(B_star != NULL);

		xdist = (double) B_star->x - (double) A_star->x;
		ydist = (double) B_star->y - (double) A_star->y;
		dx[i] = xdist;
		dy[i] = ydist;
#ifdef DEBUG
		printf("  medtf:  %4d %4d  xa %7.4f  xb %7.4f   ya %7.4f  yb %7.4f\n",
				A_star->id, B_star->id, (double) A_star->x, (double) B_star->x,
				(double) A_star->y, (double) B_star->y);
		printf("  medtf:  %4d   dx %7.4f  dy %7.4f \n", i, xdist, ydist);
#endif
		Dx_sum += xdist;
		Dy_sum += ydist;
		Dx_sum2 += xdist * xdist;
		Dy_sum2 += ydist * ydist;

		A_star = A_star->next;
		B_star = B_star->next;
	}

	/*
	 * Step 2: calculate the mean distances and (unclipped) stdev
	 */
	Dx_ave = Dx_sum / nstar;
	Dy_ave = Dy_sum / nstar;
	Dx_rms = sqrt(Dx_sum2 / nstar - Dx_ave * Dx_ave);
	Dy_rms = sqrt(Dy_sum2 / nstar - Dy_ave * Dy_ave);

	/*
	 * Step 3: calculate the median distances
	 */
	qsort((char *) dx, nstar, sizeof(double), (PFI) compare_double);
	Dx_med = find_percentile(dx, nstar, (double) 0.50);
	qsort((char *) dy, nstar, sizeof(double), (PFI) compare_double);
	Dy_med = find_percentile(dy, nstar, (double) 0.50);

	/*
	 * Step 4 (if desired): recalculate statistics, using only pairs of
	 *                      stars separated by 'medsigclip' standard
	 *                      deviations from the mean.
	 */
	if (medsigclip > 0) {
		if ((Dx_rms <= 0.0) || (Dy_rms <= 0.0)) {
			shError(
					"atFindMedtf: RMS <= 0.0, so can't calculate clipped values");
		} else {

			/*
			 * we need another pair of arrays, this time containing only
			 * those distances within the clipping criterion
			 */
			dxclip = (double *) shMalloc(nstar * sizeof(double));
			dyclip = (double *) shMalloc(nstar * sizeof(double));

			/* calculate the maximum acceptable distance */
			clip = medsigclip * sqrt(Dx_rms * Dy_rms);
			if ((fabs(Dx_med - Dx_ave) > 0.5 * clip)
					|| (fabs(Dy_med - Dy_ave) > 0.5 * clip)) {
				shError("atFindMedtf: dangerous skewness in shifts");
			}

			/*
			 * recalculate statistics, discarding values outside the
			 * clipping criterion
			 */
			Dx_sum = 0.0;
			Dy_sum = 0.0;
			Dx_sum2 = 0.0;
			Dy_sum2 = 0.0;
			num_within_clip = 0;
			for (i = 0; i < nstar; i++) {
				if (fabs(dx[i] - Dx_med) > clip) {
					continue;
				}
				if (fabs(dy[i] - Dy_med) > clip) {
					continue;
				}

				xdist = dx[i];
				ydist = dy[i];
				dxclip[num_within_clip] = xdist;
				dyclip[num_within_clip] = ydist;
				Dx_sum += xdist;
				Dy_sum += ydist;
				Dx_sum2 += xdist * xdist;
				Dy_sum2 += ydist * ydist;
				num_within_clip++;
			}

			Dx_ave = Dx_sum / num_within_clip;
			Dy_ave = Dy_sum / num_within_clip;
			Dx_rms = sqrt(Dx_sum2 / num_within_clip - Dx_ave * Dx_ave);
			Dy_rms = sqrt(Dy_sum2 / num_within_clip - Dy_ave * Dy_ave);

			qsort((char *) dxclip, num_within_clip, sizeof(double),
					(PFI) compare_double);
			Dx_med = find_percentile(dxclip, num_within_clip, (double) 0.50);
			qsort((char *) dyclip, num_within_clip, sizeof(double),
					(PFI) compare_double);
			Dy_med = find_percentile(dyclip, num_within_clip, (double) 0.50);

			shFree(dxclip);
			shFree(dyclip);

			nstar = num_within_clip;
		}
	}

	/* finally, set the values in the MEDTF structure */
	medtf->mdx = Dx_med;
	medtf->mdy = Dy_med;
	medtf->adx = Dx_ave;
	medtf->ady = Dy_ave;
	medtf->sdx = Dx_rms;
	medtf->sdy = Dy_rms;
	medtf->nm = nstar;

	/* free the memory we've allocated */
	shFree(dx);
	shFree(dy);

	return (SH_SUCCESS);
}

/***************************************************************************
 * PROCEDURE: atCalcRMS
 * Added 17 Jan 2002 by JPB
 *
 * DESCRIPTION:
 * This routine takes two matched lists and calculates the
 * rms of the differences.  Just for simple diagnostic purposes.
 *
 * If the two lists are empty, set the output args to 0.00 and
 * return with code SH_SUCCESS.
 *
 * RETURN:
 *    SH_SUCCESS          if all goes well (even if lists were empty)
 *    SH_GENERIC_ERROR    if not
 *
 */
int atCalcRMS(int num_A, /* I: number of matched stars in list A */
struct s_star *mlistA, /* I: list of matched stars from set A */
int num_B, /* I: number of matched stars in list B */
struct s_star *mlistB, /* I: list of matched stars from set B */
double *Dx_rms, double *Dy_rms) {
	double Dxterm, Dyterm, Dx_sum2, Dy_sum2, xms, yms;
	int ii, Nstar, Ntoss;
	struct s_star *A_star, *B_star;

	shAssert(num_A == num_B);
	if (num_A == 0) {
		*Dx_rms = 0.0;
		*Dy_rms = 0.0;
		return (SH_SUCCESS);
	}

	shAssert(mlistA != NULL);
	shAssert(mlistB != NULL);
	Nstar = num_A;
	Dx_sum2 = Dy_sum2 = 0.0;

	for (ii = 0, A_star = mlistA, B_star = mlistB; ii < Nstar; ii++, A_star =
			A_star->next, B_star = B_star->next) {
		shAssert(A_star != NULL);
		shAssert(B_star != NULL);
		Dxterm = (double) (B_star->x - A_star->x);
		Dyterm = (double) (B_star->y - A_star->y);
		Dx_sum2 += Dxterm * Dxterm;
		Dy_sum2 += Dyterm * Dyterm;
	}
	xms = (Dx_sum2 / Nstar); /* these are mean-squared! */
	yms = (Dy_sum2 / Nstar);

	/* Ok, we just do a quick conservative 3-sigma clip here */

	Dx_sum2 = Dy_sum2 = 0.0;
	for (Ntoss = ii = 0, A_star = mlistA, B_star = mlistB; ii < Nstar;
			ii++, A_star = A_star->next, B_star = B_star->next) {
		Dxterm = (double) (B_star->x - A_star->x);
		Dyterm = (double) (B_star->y - A_star->y);

		/* Note: squaring these terms, unlike loop above! */
		Dxterm *= Dxterm;
		Dyterm *= Dyterm;

		if (Dxterm < 9 * xms && Dyterm < 9 * yms) {
			Dx_sum2 += Dxterm;
			Dy_sum2 += Dyterm;
		} else {
			Ntoss++;
		}
	}
	if (Dx_sum2 <= 0.0) {
		*Dx_rms = 0.0;
	} else {
		*Dx_rms = sqrt(Dx_sum2 / (Nstar - Ntoss));
	}
	if (Dy_sum2 <= 0.0) {
		*Dy_rms = 0.0;
	} else {
		*Dy_rms = sqrt(Dy_sum2 / (Nstar - Ntoss));
	}

	return (SH_SUCCESS);
}

/***************************************************************************
 * PROCEDURE: is_desired_rotation
 *
 * DESCRIPTION:
 * We are given two triangles and a desired relative rotation angle
 * (and tolerance).
 * Our job is to determine if the two triangles are rotated relative
 * to each other by the given angle, within the tolerance.
 *
 * This is a little tricky, because we need to allow for angles
 * wrapping around the range of (-pi, +pi) in radians, or (-180, +180)
 * in degrees.  We therefore make two checks, first without any
 * wrapping, and then with the appropriate wrapping; if either check
 * succeeds, we return with SH_SUCCESS.
 *
 * We place the actual measured rotation angle (in degrees) between
 * these two triangles into the 'actual_angle' argument.
 *
 * RETURN:
 *    1            if the two triangles are rotated by the desired angle
 *    0            if not
 *
 * </AUTO>
 */

static int is_desired_rotation(struct s_triangle *tri_A, /* I: first triangle */
struct s_triangle *tri_B, /* I: second triangle */
double want_angle_deg, /* I: desired rotation, in degrees */
double tolerance_deg, /* I: accept any rotation within this amount */
/* I:   of desired value, in degrees */
double *actual_angle_deg /* O: the measured rotation angle between */
/*      the triangles, in degrees */
) {
	int is_good_angle;
	double min_angle_deg, max_angle_deg, wrapped_delta_deg;
	double delta_angle, delta_angle_deg;

	is_good_angle = 0;
	min_angle_deg = want_angle_deg - tolerance_deg;
	max_angle_deg = want_angle_deg + tolerance_deg;

	delta_angle = tri_A->side_a_angle - tri_B->side_a_angle;
	delta_angle_deg = delta_angle * (180.0 / 3.15159);
#ifdef DEBUG2
	printf("delta_angle_deg %7.2f ", delta_angle_deg);
#endif

	if ((delta_angle_deg >= min_angle_deg)
			&& (delta_angle_deg <= max_angle_deg)) {
#ifdef DEBUG2
		printf("   is good ");
#endif
		is_good_angle = 1;
	}
#ifdef DEBUG2
	else {
		printf("   is bad  ");
	}
#endif

	/*
	 * now we check to see if the wrapped version of the angle falls
	 * into the desired range
	 */
	if (delta_angle_deg > 0) {
		wrapped_delta_deg = delta_angle_deg - 360.0;
	} else {
		wrapped_delta_deg = delta_angle_deg + 360.0;
	}
#ifdef DEBUG2
	printf("wrapped_delta_deg %7.2f ", wrapped_delta_deg);
#endif
	if ((wrapped_delta_deg >= min_angle_deg)
			&& (wrapped_delta_deg <= max_angle_deg)) {
#ifdef DEBUG2
		printf("   is good\n");
#endif
		is_good_angle = 1;
	}
#ifdef DEBUG2
	else {
		printf("   is bad \n");
	}
#endif

	*actual_angle_deg = delta_angle_deg;

	if (is_good_angle == 1) {
		return (1);
	} else {
		return (0);
	}
}

/************************************************************************
 *
 *
 * ROUTINE: sort_triangle_by_yt
 *
 * DESCRIPTION:
 * Given an array of "num" s_triangle structures, sort it in order
 * of decreasing "yt" value (where "yt" is the ratio of lengths
 * of longest to shortest side)
 *
 * Calls the "compare_triangle_by_yt" function, below.
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

static void sort_triangle_by_yt(s_triangle *array, /* I: array of structures to be sorted */
int num /* I: number of triangles in the array */
) {
	qsort((char *) array, num, sizeof(s_triangle),
			(PFI) compare_triangle_by_yt);
}

/************************************************************************
 *
 *
 * ROUTINE: compare_triangle_by_yt
 *
 * DESCRIPTION:
 * Given two s_triangle structures, compare their "yt" values.
 * Used by "sort_triangle_by_yt".
 *
 * RETURN:
 *    1                  if first star has smaller "yt"
 *    0                  if the two have equal "yt"
 *   -1                  if the first has larger "yt"
 *
 * </AUTO>
 */

static int compare_triangle_by_yt(s_triangle *triangle1, /* I: compare "yt" field of THIS triangle ... */
s_triangle *triangle2 /* I:  ... with THIS triangle  */
) {
	shAssert((triangle1 != NULL) && (triangle2 != NULL));

	if (triangle1->yt < triangle2->yt) {
		return (1);
	}
	if (triangle1->yt > triangle2->yt) {
		return (-1);
	}
	return (0);
}

/************************************************************************
 *
 *
 * ROUTINE: sort_triangle_by_D
 *
 * DESCRIPTION:
 * Given an array of "num" s_triangle structures, sort it in order
 * of decreasing "D" value (where "D" is the product of the
 * xt and yt fields)
 *
 * Calls the "compare_triangle_by_D" function, below.
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

static void sort_triangle_by_D(s_triangle *array, /* I: array of structures to be sorted */
int num /* I: number of triangles in the array */
) {
	qsort((char *) array, num, sizeof(s_triangle), (PFI) compare_triangle_by_D);
}

/************************************************************************
 *
 *
 * ROUTINE: compare_triangle_by_D
 *
 * DESCRIPTION:
 * Given two s_triangle structures, compare their "D" values.
 * Used by "sort_triangle_by_D".
 *
 * RETURN:
 *    1                  if first star has smaller "D"
 *    0                  if the two have equal "D"
 *   -1                  if the first has larger "D"
 *
 * </AUTO>
 */

static int compare_triangle_by_D(s_triangle *triangle1, /* I: compare "D" field of THIS triangle ... */
s_triangle *triangle2 /* I:  ... with THIS triangle  */
) {
	shAssert((triangle1 != NULL) && (triangle2 != NULL));

	if (triangle1->D < triangle2->D) {
		return (1);
	}
	if (triangle1->D > triangle2->D) {
		return (-1);
	}
	return (0);
}

/************************************************************************
 *
 *
 * ROUTINE: find_quick_match
 *
 * DESCRIPTION:
 * Given two arrays of triangles, and the arrays of stars that make
 * up each set of triangles, try to match up triangles in the two
 * arrays.  Triangles can be considered to match only when the
 * Euclidean distance in "triangle space", created from the two
 * coordinates "ba" and "ca", is within "max_radius".  That is,
 * for two triangles to match, we must satisfy
 *
 *     sqrt[ (t1.ba - t2.ba)^2 + (t1.ca - t2.ca)^2 ] <= max_radius
 *
 * Note that there may be more than one triangle from array A which
 * matches a particular triangle from array B!  That's okay --
 * we treat any 2 which satisfy the above equation as "matched".
 * We rely upon the "vote_array" to weed out false matches.
 *
 * If "min_scale" and "max_scale" are not both -1, then disallow
 * any match for which the
 * ratio of triangles (indicated by "a_length" members)
 * is outside the given values.
 *
 * If "rotation_deg" and "tolerance_deg" are not both AT_MATCH_NOANGLE,
 * then disallow any match for which the two triangles
 * are not oriented at an angle of "rotation_deg" degrees
 * relative to each other (with a tolerance of "tolerance_deg" degrees).
 *
 * For each pair of triangles that matches, increment
 * the "vote" in each "vote cell" for each pair of matching
 * vertices.
 *
 * RETURN:
 *           SH_SUCCESS            if a good match is found
 *           SH_GENERIC_FAILURE    if a good match is NOT found
 *
 * </AUTO>
 */

static int find_quick_match(s_star *star_array_A, /* I: first array of stars */
int num_stars_A, /* I: number of stars in star_array_A  */
s_star *star_array_B, /* I: second array of stars */
int num_stars_B, /* I: number of stars in star_array_B  */
s_triangle *t_array_A, /* I: array of triangles from star_array_A */
int num_triangles_A, /* I: number of triangles in t_array_A */
s_triangle *t_array_B, /* I: array of triangles from star_array_B */
int num_triangles_B, /* I: number of triangles in t_array_B */
int nbright, /* I: consider at most this many stars */
/*       from each array; also the size */
/*       of the output "vote_matrix". */
double star_match_radius, /* I: max radius in star-space allowed */
/*       for 2 stars to be considered */
/*       a matching pair. */
double max_radius, /* I: max radius in triangle-space allowed */
/*       for 2 triangles to be considered */
/*       a matching pair. */
int max_iterations, /* I: iterate at most this many times */
/*       while refining the TRANS */
double max_sigma, /* I: the stdev of matched pairs must be */
/*       this small for successful match */
double min_scale, /* I: minimum permitted relative scale factor */
/*       if -1, any scale factor is allowed */
double max_scale, /* I: maximum permitted relative scale factor */
/*       if -1, any scale factor is allowed */
double rotation_deg, /* I: desired relative angle of coord systems (deg) */
/*       if AT_MATCH_NOANGLE, any orientation is allowed */
double tolerance_deg, /* I: allowed range of orientation angles (deg) */
/*       if AT_MATCH_NOANGLE, any orientation is allowed */
int min_req_pairs, /* I: must have this many matched pairs */
/*       of stars to count as successful match */
TRANS *output_trans /* I/O: "order" field must be set on input */
/*       other fields filled if success for output */
) {
	int i_triA;
	int failure_flag;
	struct s_triangle *triA, *triB;
	struct s_star_coord *star_coord_array_B;

#ifdef DEBUG
	printf(" entering find_quick_match \n");
#endif

	shAssert(star_array_A != NULL);
	shAssert(star_array_B != NULL);
	shAssert(t_array_A != NULL);
	shAssert(t_array_B != NULL);
	shAssert(nbright > 0);
	if (min_scale != -1) {
		shAssert((max_scale != -1) && (min_scale <= max_scale));
	}
	if (max_scale != -1) {
		shAssert((min_scale != -1) && (min_scale <= max_scale));
	}

	/*
	 * It will increase efficiency later on
	 *    (in apply_trans_and_find_matches) to have a
	 *    way to access the elements of star_array_B
	 *    in order of their "X" coordinates.
	 */
	{
		int i;

		star_coord_array_B = shMalloc(
				sizeof(struct s_star_coord) * num_stars_B);
		for (i = 0; i < num_stars_B; i++) {
			star_coord_array_B[i].index = i;
			star_coord_array_B[i].x = star_array_B[i].x;
			star_coord_array_B[i].y = star_array_B[i].y;
		}
		sort_star_coord_by_x(star_coord_array_B, num_stars_B);
#ifdef DEBUG2
		printf("  here comes star_coord_array_B sorted by x \n");
		for (i = 0; i < num_stars_B; i++) {
			printf("   i %5d  index %5d  x %9.4e  y %9.4e \n", i,
					star_coord_array_B[i].index, star_coord_array_B[i].x,
					star_coord_array_B[i].y);
		}
#endif
	}

	/*
	 * We walk through triangles in list A, which have been sorted by
	 *    the values of D = dot product xt*yt.  In other words,
	 *    we start checking the (rare) large triangles, since they
	 *    are more likely to yield a match quickly.
	 */
#ifdef DEBUG2
	printf("  find_quick_match: about to walk through A  %6d\n",
			num_triangles_A);
#endif
	for (i_triA = 0; i_triA < num_triangles_A; i_triA++) {

		int start_index, end_index, b_index;
		double yt_eps;

#ifdef DEBUG2
		printf(" top of find_quick_match, i_triA = %6d \n", i_triA);
#endif

		triA = &(t_array_A[i_triA]);
#ifdef DEBUG2
		print_one_triangle(triA, star_array_A);
#endif

		/*
		 * We want to search for a triangle in list B which
		 *    matches triA.  We make a binary search
		 *    inside the t_array_B using the "yt" values
		 *    to choose the first and last triangle in
		 *    t_array_B which could possibly be a
		 *    match to triA.
		 *
		 * Note that we slide the index one element
		 *    toward the extremes to make our
		 *    starting and ending places definitely
		 *    contain the region of interest.
		 */
		yt_eps = triA->yt * (AT_QUICK_YT_PERCENT * 0.01);
		start_index = find_yt_triangle(t_array_B, num_triangles_B,
				triA->yt + yt_eps);
		if (start_index > 0) {
			start_index--;
		}
		end_index = find_yt_triangle(t_array_B, num_triangles_B,
				triA->yt - yt_eps);
		if (end_index < (num_triangles_B - 1)) {
			end_index++;
		}
#ifdef DEBUG2
		printf("  find_quick_match: start_index %5d end_index %5d \n",
				start_index, end_index);
#endif

		/*
		 * Okay, we have a range of triangles in list B
		 *    we can try to match to the current triA from
		 *    list A.
		 *
		 */
		for (b_index = start_index; b_index <= end_index; b_index++) {
			double ba_diff, ca_diff, cb_diff;
			double actual_angle_deg, ratio;

			triB = &(t_array_B[b_index]);
#ifdef DEBUG2
			printf("      i_triA is %6d, looking at triangle B %5d \n", i_triA,
					triB->id);
#endif

			ba_diff = triA->ba - triB->ba;
			if (fabs(ba_diff) > AT_QUICK_RATIO_DIFF) {
#ifdef DEBUG2
				printf(
						"  find_quick_match: ba values %9.4e %9.4e diff %6.3f > %6.3f \n",
						triA->ba, triB->ba, ba_diff, AT_QUICK_RATIO_DIFF);
#endif
				continue;
			}
			ca_diff = triA->ca - triB->ca;
			if (fabs(ca_diff) > AT_QUICK_RATIO_DIFF) {
#ifdef DEBUG2
				printf(
						"  find_quick_match: ca values %9.4e %9.4e diff %6.3f > %6.3f \n",
						triA->ca, triB->ca, ca_diff, AT_QUICK_RATIO_DIFF);
#endif
				continue;
			}
			cb_diff = triA->cb - triB->cb;
			if (fabs(cb_diff) > AT_QUICK_RATIO_DIFF) {
#ifdef DEBUG2
				printf(
						"  find_quick_match: cb values %9.4e %9.4e diff %6.3f > %6.3f \n",
						triA->cb, triB->cb, cb_diff, AT_QUICK_RATIO_DIFF);
#endif
				continue;
			}
#ifdef DEBUG2
			printf("  find_quick_match: good match for A %5d B %5d \n",
					triA->id, triB->id);
#endif

			/*
			 * check the relative orientations of the triangles.
			 * If they don't match the desired rotation angle,
			 * discard this match.
			 */
			if (rotation_deg != AT_MATCH_NOANGLE) {
				if (is_desired_rotation(triA, triB, rotation_deg, tolerance_deg,
						&actual_angle_deg) == 0) {
#ifdef DEBUG2
					printf(
							"  find_quick_match: bad angle %lf is not desired %lf \n",
							actual_angle_deg, rotation_deg);
#endif
					continue;
				}
			}

			/*
			 * check the ratio of lengths of side "a", and discard this
			 * candidate if it is outside the allowed range
			 */
			if (min_scale != -1) {
				ratio = triA->a_length / triB->a_length;
				if (ratio < min_scale || ratio > max_scale) {
#ifdef DEBUG2
					printf(
							"  find_quick_match: bad scale %lf is outside %lf  %lf \n",
							ratio, min_scale, max_scale);
#endif
					continue;
				}
			}

#ifdef DEBUG2
			printf("  find_quick_match: match found for A %5d = B %5d \n",
					triA->id, triB->id);
#endif

			/*
			 * now we are going to put these two triangles to
			 *    more stringent tests
			 *        triA   =   pointer to triangle in list A
			 *        triB   =   pointer to MATCHING triangle in list B
			 */

			/*
			 * okay, we have a pair of triangles which match.
			 * This might be a "true" match, which reveals the actual
			 * transformation between the two coordinate systems;
			 * or it might be an unlucky false match.
			 *
			 * To determine if this is a "true" match, we will
			 *
			 *    a) use the 3 stars in each triangle to compute
			 *         a linear TRANS between the two coord systems
			 *
			 *    b) apply the TRANS to all stars in the two lists
			 *
			 *    c) count the number of stars which end up matching
			 *         in the transformed lists as a result.
			 */

			failure_flag = 0;

			/*
			 * the "winner_index" arrays give the ID numbers of
			 *   corresponding stars from each list.
			 */
			{
				int *winner_index_A, *winner_index_B;
				int nbright;
				int *winner_votes, num_winners = 0;
				TRANS *test_trans;

				/*
				 * right now, we'll only need 3 elements for these
				 *   arrays, but later on in this routine, we'll
				 *   need enough space for every star in each list.
				 */
				winner_index_A = shMalloc(num_stars_A * sizeof(int));
				winner_index_B = shMalloc(num_stars_B * sizeof(int));
				if (num_stars_A > num_stars_B) {
					winner_votes = shMalloc(num_stars_A * sizeof(int));
				} else {
					winner_votes = shMalloc(num_stars_B * sizeof(int));
				}

				winner_index_A[0] = triA->a_index;
				winner_index_A[1] = triA->b_index;
				winner_index_A[2] = triA->c_index;
				winner_index_B[0] = triB->a_index;
				winner_index_B[1] = triB->b_index;
				winner_index_B[2] = triB->c_index;
				nbright = 3;
#ifdef DEBUG2
				printf("  list A stars %5d=%-5d %5d=%-5d %5d=%-5d \n",
						winner_index_A[0], star_array_A[winner_index_A[0]].id,
						winner_index_A[1], star_array_A[winner_index_A[1]].id,
						winner_index_A[2], star_array_A[winner_index_A[2]].id);
				printf("  list B stars %5d=%-5d %5d=%-5d %5d=%-5d \n",
						winner_index_B[0], star_array_B[winner_index_B[0]].id,
						winner_index_B[1], star_array_B[winner_index_B[1]].id,
						winner_index_B[2], star_array_B[winner_index_B[2]].id);
#endif

				/* we create a linear TRANS just for this test */
				test_trans = atTransNew();
				test_trans->order = AT_TRANS_LINEAR;

				/*
				 * we iterate the following process:
				 *
				 *    1. find a TRANS
				 *    2. apply TRANS to all stars, and look for matches
				 *          (and then prune poor matches)
				 *    3. compute statistics of the matched pairs
				 *    4. decide whether to quit or continue
				 *       4a.   quit, or
				 *       4b.   go to step 1
				 */
				{
					int iter, num_iter;

					num_iter = max_iterations;

					for (iter = 0; iter < num_iter; iter++) {

						if (calc_trans(nbright, star_array_A, num_stars_A,
								star_array_B, num_stars_B, winner_votes,
								&(winner_index_A[0]), &(winner_index_B[0]),
								test_trans) != SH_SUCCESS) {
#ifdef DEBUG
							printf(
									"find_quick_match: iter %3d, calc_trans fails",
									iter);
#endif
							failure_flag = 1;
							break;
						}
#ifdef DEBUG2
						printf("  find_quick_match: calc_trans succeeds \n");
						print_trans(test_trans);
#endif

						/*
						 * Does this TRANS satisfy the constraints on scale
						 * and rotation angle which the user specified?
						 * If no, we discard this TRANS and start all over.
						 */
						if (check_trans_properties(test_trans, min_scale,
								max_scale, rotation_deg, tolerance_deg)
								!= SH_SUCCESS) {
#ifdef DEBUG2
							printf(
									"check_trans_properties fails, so give up on this TRANS \n");
#endif
							failure_flag = 1;
							break;
						}

						/*
						 * Is this really a good TRANS?  To find out, we apply
						 * it to all the stars in list A, then check to see how
						 * many match with stars in list B.
						 *
						 * After running "apply_trans_and_find_matches()",
						 *    a) the number of matching pairs will be set
						 *          in the "num_winners" variable
						 *    b) the indices of matching pairs
						 *          from each list will be set in
						 *          the "winner_index_A[]" and
						 *          "winner_index_B[[]" arrays
						 */

#ifdef DEBUG2
						printf(
								"  before apply_trans_and_find_matches, num_stars_B %5d \n",
								num_stars_B);
#endif
						if (apply_trans_and_find_matches(star_array_A,
								num_stars_A, star_array_B, num_stars_B,
								star_coord_array_B, star_match_radius,
								test_trans, &num_winners, &(winner_index_A[0]),
								&(winner_index_B[0])) != SH_SUCCESS) {
#ifdef DEBUG2
							printf("apply_trans_and_find_matches fails \n");
#endif
						} else {
#ifdef DEBUG2
							printf("apply_trans_and_find_matches succeeds \n");
							printf("   num_winners = %5d \n", num_winners);
#endif

						}
#ifdef DEBUG2
						printf(
								"  after  apply_trans_and_find_matches, num_stars_B %5d \n",
								num_stars_B);
#endif

						nbright = num_winners;

						/* evaluate the quality of the matches given by this TRANS */
						if (eval_trans_quality(star_array_A, num_stars_A,
								star_array_B, num_stars_B, star_match_radius,
								test_trans) != SH_SUCCESS) {
							printf("eval_trans_quality fails ?!\n");
							return (SH_GENERIC_ERROR);
						}
#ifdef DEBUG
						printf(" after iter %2d, here comes TRANS \n", iter);
						print_trans(test_trans);
#endif

						/*
						 * at this point, we could could eliminate
						 * some of the matching pairs
						 */

						/*
						 * we used a linear TRANS for the initial comparison of stars,
						 * but if the user requested a higher-order TRANS, we now
						 * switch to using that higher-order for all subsequent
						 * iterations.
						 */
						if (output_trans->order != AT_TRANS_LINEAR) {
#ifdef DEBUG2
							printf(
									"   iter %d: switching from order 1 to order %d  \n",
									iter, output_trans->order);
#endif
							test_trans->order = output_trans->order;
						}

					} /* end of loop over iterations of calc_trans/apply_trans */

				}

				if (failure_flag == 1) {
					continue;
				}

				/* we should check here to see if "num_winners" is high enough ... */

				// nbright = num_winners;

				/*
				 * At this point, we have created the best TRANS we can, starting
				 * with the current pair of triangles.  The question is "is the
				 * TRANS good enough for us to declare success and return?"
				 *
				 *    If the answer is yes, we break out of loop, return
				 *    If no, we keep walking through the list of triangles,
				 *           looking for another pair of matching triangles
				 *           to serve as a starting point.
				 */
				if (is_trans_good_enough(min_req_pairs, max_sigma, test_trans)
						== SH_SUCCESS) {
#ifdef DEBUG2
					printf(" is_trans_good_enough returns yes \n");
#endif

					/* make sure the TRANS has proper scale and rotation */
					if (check_trans_properties(test_trans, min_scale, max_scale,
							rotation_deg, tolerance_deg) == SH_SUCCESS) {

						/* yes!  It's good enough! */

						/* assign the properties of the test_trans to output_trans */
						copyTrans(test_trans, output_trans);

						/* and return immediately -- no need to look for better matches */
						return (SH_SUCCESS);

					} else {
#ifdef DEBUG2
						printf(
								" but check_trans_properties fails, so give up on this TRANS \n");
#endif
					}
				} else {
#ifdef DEBUG2
					printf(" is_trans_good_enough returns no \n");
#endif
				}

			}

		} /* end of loop over all triangles in B list */

	} /* end of loop over all triangles in A list */

	/* if we reach this point, we did NOT find a good match */
	return (SH_GENERIC_ERROR);
}

/************************************************************************
 *
 *
 * ROUTINE: find_yt_triangle
 *
 * DESCRIPTION:
 * Given an array of "num" s_triangle structures, which have already
 * been sorted in order of DEcreasing "yt" value, and given one
 * particular "yt" value yt0, return the index of the first triangle
 * in the array which has "yt" <= yt0.
 *
 * We use a binary search, on the "yt" element of each structure.
 *
 * If there is no such triangle, just return the index of the last
 * triangle in the list.
 *
 * RETURN:
 *    index of closest triangle in array         if all goes well
 *    index of last triangle in array            if nothing close
 *
 * </AUTO>
 */

static int find_yt_triangle(s_triangle *array, /* I: array of structures which been sorted */
int num, /* I: number of triangles in the array */
double yt0 /* I: value of "yt" we seek */
) {
	int top, bottom, mid;

#ifdef DEBUG2
	printf("find_yt_triangle: looking for yt = %10.2f\n", yt0);
#endif

	top = 0;
	if ((bottom = num - 1) < 0) {
		bottom = 0;
	}

	while (bottom - top > 1) {
		mid = (top + bottom) / 2;
#ifdef DEBUG2
		printf(
				" array[%4d] yt=%10.2f   array[%4d] yt=%10.2f  array[%4d] yt=%10.2f\n",
				top, array[top].yt, mid, array[mid].yt, bottom,
				array[bottom].yt);
#endif
		if (array[mid].yt > yt0) {
			top = mid;
		} else {
			bottom = mid;
		}
	}
#ifdef DEBUG2
	printf(
			" array[%4d] yt=%10.2f                              array[%4d] yt=%10.2f\n",
			top, array[top].yt, bottom, array[bottom].yt);
#endif

	/*
	 * if we get here, then the item we seek is either "top" or "bottom"
	 * (which may point to the same item in the array).
	 */
	if (array[top].yt > yt0) {
#ifdef DEBUG2
		printf(" returning array[%4d] yt=%10.2f \n", bottom, array[bottom].yt);
#endif
		return (bottom);
	} else {
#ifdef DEBUG2
		printf(" returning array[%4d] yt=%10.2f \n", top, array[top].yt);
#endif
		return (top);
	}
}

/************************************************************************
 *
 *
 * ROUTINE: apply_trans_and_find_matches
 *
 * DESCRIPTION:
 * We already have a candidate TRANS structures that takes coords of
 * objects in set A and transforms to coords of objects in set B.
 * We want to find out if this TRANS is really a good one.
 * The basic idea is: apply the TRANS to all objects in list A,
 * check to see how many match an object in list B.
 *
 * In order to speed up the checks for matches, we create
 * a copy of the stars in list B and sort the copy by
 * the "x" coordinate.  As we walk through the list of
 * transformed stars from list A, we can use a binary
 * search by their "x" coordinate to pick out the small
 * subset of stars in list B which might be matches.
 *
 * When we do find a matching pair of stars, we
 *
 *   a) increment "num_winners"
 *
 *   b) set elements of the "winner_index_A[]"
 *         and "winner_index_B[]" arrays to contain
 *         the indices of the two stars.  The index
 *         for list B refers to the original list B,
 *         not the one we've sorted by "x" coord.
 *
 *
 * RETURNS:
 *   SH_SUCCESS          if no errors occur
 *   SH_GENERIC_ERROR    if something goes wrong
 *
 * </AUTO>
 */

static int apply_trans_and_find_matches(s_star *star_array_A, /* I: first array of s_star structure we match */
/*      the TRANS takes their coords */
/*      into those of array B */
int num_stars_A, /* I: total number of stars in star_array_A */
s_star *star_array_B, /* I: second array of s_star structure we match */
int num_stars_B, /* I: total number of stars in star_array_B */
s_star_coord *sorted_B, /* I: elements of star_array_B which have been */
/*      sorted by the "x" coordinate value */
double star_match_radius, /* I: max distance in star-space for two */
/*      stars to be considered a match */
TRANS *trans, /* I: use this TRANS to transform coords of */
/*      stars in list A to system of B */
int *num_winners, /* O: number of matching pairs we find */
int *winner_index_A, /* O: index into "star_array_A" of stars */
/*      which are part of a matched pair */
int *winner_index_B /* O: index into "star_array_B" of stars */
/*      which are part of a matched pair */
) {
	int i, j;
	int num_matched;
	int closest_B_index;
	double star_match_radius_sq, closest_dist_sq;
	double mean, stdev;
	s_star *star_A, *star_B;
	s_star *transformed_star_array_A;

#ifdef DEBUG2
	int start_index, end_index;
	printf("entering apply_trans_and_find_matches \n");
#endif

	num_matched = 0;
	star_match_radius_sq = star_match_radius * star_match_radius;

	/*
	 * create a copy of stars in list A
	 *    and apply the TRANS to them
	 */
	transformed_star_array_A = shMalloc(sizeof(s_star) * num_stars_A);
	copy_star_array(star_array_A, transformed_star_array_A, num_stars_A);
	if (apply_trans(transformed_star_array_A, num_stars_A, trans)
			!= SH_SUCCESS) {
		shError("apply_trans_and_find_matches: apply_trans fails on list A");
		return (SH_GENERIC_ERROR);
	}

	/* walk through the list of all stars in transformed list */
	for (i = 0; i < num_stars_A; i++) {

		int start_sc_index, end_sc_index;
		double x;

		star_A = &(transformed_star_array_A[i]);
#ifdef DEBUG2
		printf("  star A has x %9.4e \n", star_A->x);
#endif

		/*
		 * use binary search to locate starting and ending
		 *     elements we need to check for matches in list B
		 */
		x = star_A->x - star_match_radius;
		start_sc_index = find_star_coord_by_x(sorted_B, num_stars_B, x);
		if (start_sc_index > 0) {
			start_sc_index--;
		}
#ifdef DEBUG2
		start_index = sorted_B[start_sc_index].index;
#endif

		x = star_A->x + star_match_radius;
		end_sc_index = find_star_coord_by_x(sorted_B, num_stars_B, x);
		if (end_sc_index < num_stars_B - 1) {
			end_sc_index++;
		}

#ifdef DEBUG2
		end_index = sorted_B[end_sc_index].index;
		printf(
				"    start_index %5d %5d  x %9.4e   end_index %5d %5d x %9.4e \n",
				start_sc_index, start_index, star_array_B[start_index].x,
				end_sc_index, end_index, star_array_B[end_index].x);
#endif

		/*
		 * now go through all the stars in sorted list B which lie in the range
		 *   of possible matches to star A.  Find the star in sorted list B
		 *   which is closest to star A, and set the "closest_B_index"
		 *   value to its index in star_array_B[].
		 */
		closest_B_index = -1;
		closest_dist_sq = star_match_radius_sq * 2;
		for (j = start_sc_index; j <= end_sc_index; j++) {

			double dx, dy, distsq;

			star_B = &(star_array_B[sorted_B[j].index]);
#ifdef DEBUG
			printf("     star_B is index %5d \n", j);
#endif

			dx = fabs(star_A->x - star_B->x);
			if (dx > star_match_radius) {
#ifdef DEBUG
				printf("      dx %11.4f > %11.4f so go to next B \n", dx,
						star_match_radius);
#endif
				continue;
			}

			dy = star_A->y - star_B->y;
			distsq = dx * dx + dy * dy;
			if (distsq < closest_dist_sq) {
#ifdef DEBUG2
				printf("      new closest distsq %11.4e \n", distsq);
#endif
				closest_dist_sq = distsq;
				closest_B_index = sorted_B[j].index;
			}
		}

		/*
		 * now, check to see if the closest star was within the
		 *    star_match_radius.  If yes, we have a match;
		 *    If not, we have no match.
		 */
		if (closest_dist_sq < star_match_radius_sq) {
#ifdef DEBUG2
			printf("       good match: star A %5d B %5d  dist %11.4e \n", i, j,
					sqrt(closest_dist_sq));
#endif
			winner_index_A[num_matched] = i;
			winner_index_B[num_matched] = closest_B_index;
			num_matched++;
		} else {
#ifdef DEBUG2
			printf("       no   match: star A %5d  min dist %11.4e \n", i,
					sqrt(closest_dist_sq));
#endif
		}

	}

#ifdef DEBUG2
	printf("  apply_trans_and_find_matches: %5d matches  ", num_matched);
#endif

	*num_winners = num_matched;

	{

		/*
		 * compute the mean and stdev of distance between stars
		 * in the matched pairs
		 */
		if (compute_match_distance_stats(transformed_star_array_A, num_stars_A,
				star_array_B, num_stars_B, num_matched, winner_index_A,
				winner_index_B, &mean, &stdev) != 0) {
			printf(
					"apply_trans_and_find_matches: compute_match_distance_stats fails \n");
			return (SH_GENERIC_ERROR);
		}

		/*
		 * discard pairs which have separations more than
		 * 3 sigma from the mean separation
		 */
		double critical_distance = mean + 3 * stdev;
		int remaining_pairs;

		if (prune_matched_pairs(transformed_star_array_A, num_stars_A,
				star_array_B, num_stars_B, num_matched, winner_index_A,
				winner_index_B, critical_distance, &remaining_pairs) != 0) {
#ifdef DEBUG
			printf(
					"apply_trans_and_find_matches: prune_matched_pairs fails \n");
#endif
			return (SH_GENERIC_ERROR);
		} else {
#ifdef DEBUG
			printf("  after pruning, %d pairs remain \n", remaining_pairs);
#endif
		}

		*num_winners = remaining_pairs;

	}

	return (SH_SUCCESS);
}

/************************************************************************
 *
 *
 * ROUTINE: sort_star_coord_by_x
 *
 * DESCRIPTION:
 * Given an array of "num" s_star_coord structures, sort it in order
 * of increasing "x" values.
 *
 * Calls the "compare_star_coord_by_x" function, below.
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

static void sort_star_coord_by_x(s_star_coord *array, /* I: array of structures to be sorted */
int num /* I: number of stars in the array */
) {
	qsort((char *) array, num, sizeof(s_star_coord),
			(PFI) compare_star_coord_by_x);
}

/************************************************************************
 *
 *
 * ROUTINE: compare_star_coord_by_x
 *
 * DESCRIPTION:
 * Given two s_star_coord structures, compare their "x" values.
 * Used by "sort_star_coord_by_x".
 *
 * RETURN:
 *    1                  if first star has larger "x"
 *    0                  if the two have equal "x"
 *   -1                  if the first has smaller "x"
 *
 * </AUTO>
 */

static int compare_star_coord_by_x(s_star_coord *star1, /* I: compare "x" field of THIS star ... */
s_star_coord *star2 /* I:  ... with THIS star  */
) {
	shAssert((star1 != NULL) && (star2 != NULL));

	if (star1->x > star2->x) {
		return (1);
	}
	if (star1->x < star2->x) {
		return (-1);
	}
	return (0);
}

/************************************************************************
 *
 *
 * ROUTINE: find_star_coord_by_x
 *
 * DESCRIPTION:
 * Given an array of "num" s_star_coord structures, which have already
 * been sorted in order of increasing "x" value, and given one
 * particular "x" value x0, return the "index" value of the first star
 * in the array which has "x" >= x0.
 *
 * We use a binary search, on the "x" element of each structure.
 *
 * If there is no such star, just return the index of the last
 * star in the list.
 *
 * RETURN:
 *    index of closest star_coord in array         if all goes well
 *    index of last star_coord in array            if nothing close
 *
 * </AUTO>
 */

static int find_star_coord_by_x(s_star_coord *array, /* I: array of structures which been sorted */
int num, /* I: number of elements in the array */
double x0 /* I: value of "x" we seek */
) {
	int top, bottom, mid;

#ifdef DEBUG3
	printf("find_star_coord_by_x: looking for x = %11.4e\n", x0);
#endif

	top = 0;
	if ((bottom = num - 1) < 0) {
		bottom = 0;
	}

	while (bottom - top > 1) {
		mid = (top + bottom) / 2;
#ifdef DEBUG3
		printf(
				" array[%4d] x=%11.4e   array[%4d] x=%11.4e  array[%4d] x=%11.4e\n",
				top, array[top].x, mid, array[mid].x, bottom, array[bottom].x);
#endif
		if (array[mid].x < x0) {
			top = mid;
		} else {
			bottom = mid;
		}
	}
#ifdef DEBUG3
	printf(" array[%4d] x=%11.4e                       array[%4d] x=%11.4e\n",
			top, array[top].x, bottom, array[bottom].x);
#endif

	/*
	 * if we get here, then the item we seek is either "top" or "bottom"
	 * (which may point to the same item in the array).
	 */
	if (array[top].x < x0) {
#ifdef DEBUG3
		printf(" returning array[%4d] x=%11.4e \n", bottom, array[bottom].x);
#endif
		return (bottom);
	} else {
#ifdef DEBUG3
		printf(" returning array[%4d] x=%11.4e \n", top, array[top].x);
#endif
		return (top);
	}
}

/************************************************************************
 * ROUTINE: compute_match_distance_stats
 *
 * DESCRIPTION:
 * Given two arrays of "num_matches" elements, which are the indices
 * of corresponding stars in lists A and B, go through the list and
 * for each pair
 *
 *    a) compute the distance between the stars
 *    b) add to running sums of distance and distance-squared
 *
 * so that we can can in the end compute the mean and standard
 * deviation of the distances between stars.
 *
 * Note that the stars in list A must have been TRANSformed
 * into the coordinate system of list B _before_ they are passed
 * to this routine.
 *
 * Place the mean and stdev values into the output arguments.
 *
 * RETURN:
 *    SH_SUCCESS                 if all goes well
 *    SH_GENERIC_ERROR           if an error occurs
 *
 * </AUTO>
 */

static int compute_match_distance_stats(s_star *star_array_A, /* I: first array of stars */
int num_stars_A, /* I: number of stars in star_array_A  */
s_star *star_array_B, /* I: second array of stars */
int num_stars_B, /* I: number of stars in star_array_B  */
int num_matches, /* I: number of matched pairs of stars */
int *match_index_A, /* I: array with indices of stars in list A */
int *match_index_B, /* I: array with indices of stars in list B */
double *mean, /* O: mean distance between matched stars */
double *stdev /* O: stdev of distance between matched stars */
) {
	int i;
	double sum, sumsq;

#ifdef DEBUG2
	printf("entering compute_match_distance_stats \n");
#endif

	/* sanity checks */
	if (num_matches < 1) {
		shError(
				"compute_match_distance_stats: given invalid num_matches = %d \n",
				num_matches);
		return (SH_GENERIC_ERROR);
	}
	shAssert(star_array_A != NULL);
	shAssert(star_array_B != NULL);
	shAssert(match_index_A != NULL);
	shAssert(match_index_B != NULL);

	sum = 0;
	sumsq = 0;

	for (i = 0; i < num_matches; i++) {
		double dx, dy, dist, distsq;
		s_star *star_A, *star_B;

		star_A = &(star_array_A[match_index_A[i]]);
		shAssert(star_A != NULL);
		star_B = &(star_array_B[match_index_B[i]]);
		shAssert(star_B != NULL);

		dx = star_A->x - star_B->x;
		dy = star_A->y - star_B->y;
		distsq = (dx * dx + dy * dy);
		dist = sqrt(distsq);

		sum += dist;
		sumsq += distsq;
#ifdef DEBUG2
		printf("    i %5d  dx %9.4e dy %9.4e  dist %9.4e  sum %9.4e \n", i, dx,
				dy, dist, sum);
#endif
	}

	*mean = sum / num_matches;
	if (num_matches > 1) {
		*stdev = sqrt(
				(sumsq - num_matches * (*mean) * (*mean)) / (num_matches - 1));
	} else {
		*stdev = 0.0;
	}
#ifdef DEBUG
	printf("compute_match_distance_stats: num %5d mean %9.4e stdev %9.4e \n",
			num_matches, *mean, *stdev);
#endif

	return (SH_SUCCESS);
}

/************************************************************************
 * ROUTINE: prune_matched_pairs
 *
 * DESCRIPTION:
 * Given two arrays of "num_matches" elements, which are the indices
 * of corresponding stars in lists A and B, and given some critical
 * matching distance, go through the list of matches and
 * for each pair
 *
 *    a) compute the distance between the stars
 *    b) if the distance is larger than the critical distance,
 *           discard the pair from the list
 *
 * We will modify the elements of the "match_index_A" and "match_index_B"
 * arrays.
 *
 * Place the number of remaining matched pairs into the
 * "remaining_pairs" output argument.
 *
 * RETURN:
 *    SH_SUCCESS                 if all goes well
 *    SH_GENERIC_ERROR           if an error occurs
 *
 * </AUTO>
 */

static int prune_matched_pairs(s_star *star_array_A, /* I: first array of stars */
int num_stars_A, /* I: number of stars in star_array_A  */
s_star *star_array_B, /* I: second array of stars */
int num_stars_B, /* I: number of stars in star_array_B  */
int num_matches, /* I: number of matched pairs of stars */
int *match_index_A, /* I/O: array with indices of stars in list A */
int *match_index_B, /* I/O: array with indices of stars in list B */
double critical_distance, /* I: discard pairs which differ by more */
/*       than this distance */
int *remaining_pairs /* O: number of remaining matched pairs */
) {
	int i, current_num_matches;

#ifdef DEBUG2
	printf("entering prune_matched_pairs \n");
#endif

	/* sanity checks */
	if (num_matches < 1) {
		shError("prune_matched_pairs: given invalid num_matches = %d \n",
				num_matches);
		return (SH_GENERIC_ERROR);
	}
	current_num_matches = num_matches;
	shAssert(star_array_A != NULL);
	shAssert(star_array_B != NULL);
	shAssert(match_index_A != NULL);
	shAssert(match_index_B != NULL);

	for (i = 0; i < current_num_matches; i++) {
		double dx, dy, dist, distsq;
		s_star *star_A, *star_B;

		star_A = &(star_array_A[match_index_A[i]]);
		shAssert(star_A != NULL);
		star_B = &(star_array_B[match_index_B[i]]);
		shAssert(star_B != NULL);

		dx = star_A->x - star_B->x;
		dy = star_A->y - star_B->y;
		distsq = (dx * dx + dy * dy);
		dist = sqrt(distsq);

		if (dist > critical_distance) {
			int j;

#ifdef DEBUG2
			printf("  going to remove pair i=%5d  with dist %9.4e > %9.4e \n",
					i, dist, critical_distance);
#endif

			for (j = i; j < current_num_matches - 1; j++) {
#ifdef DEBUG2
				printf("   shift match_index_A[%3d] from %5d to %5d \n", j,
						match_index_A[j], match_index_A[j + 1]);
				printf("   shift match_index_B[%3d] from %5d to %5d \n", j,
						match_index_B[j], match_index_B[j + 1]);
#endif
				match_index_A[j] = match_index_A[j + 1];
				match_index_B[j] = match_index_B[j + 1];
			}
			current_num_matches--;
			i--;
		}
	}

#ifdef DEBUG
	printf("compute_match_distance_stats: %5d pairs remain from %5d \n",
			current_num_matches, num_matches);
	for (i = 0; i < current_num_matches; i++) {
		double dx, dy, dist, distsq;
		s_star *star_A, *star_B;

		star_A = &(star_array_A[match_index_A[i]]);
		shAssert(star_A != NULL);
		star_B = &(star_array_B[match_index_B[i]]);
		shAssert(star_B != NULL);

		dx = star_A->x - star_B->x;
		dy = star_A->y - star_B->y;
		distsq = (dx * dx + dy * dy);
		dist = sqrt(distsq);

		if (dist > critical_distance) {
			shFatal(
					"prune_matched_pairs: pruned pair %5d has dist %9.4e > %9.4d ?! \n",
					i, dist, critical_distance);
		}

		printf("  pair %5d has dist %9.4e  < %9.4e \n", i, dist,
				critical_distance);
	}

#endif

	*remaining_pairs = current_num_matches;

	return (SH_SUCCESS);
}

/************************************************************************
 * ROUTINE: eval_trans_quality
 *
 * DESCRIPTION:
 * Given two lists of stars, A and B, each in its own coordinate system,
 * a matching radius, and a TRANS, we
 *
 *       1) apply the TRANS to all objects in list A, creating
 *              a temporary copy of the list with new coords in system B
 *
 *       2) create an auxiliary array of objects in list B,
 *              sorted by X coord (to speed up matching)
 *
 *       3) walk through list A, searching for matches in list B;
 *              if find several, keep only the closest one
 *
 *       4) compute statistics from the matching pairs
 *
 *       5) set fields in the TRANS to hold these statistics
 *              Nm = number of items within matching radius
 *              sig = stdev of total differences in position of matched pairs
 *              sx = stdev of differences in x position of matched pairs
 *              sy = stdev of differences in y position of matched pairs
 *
 * RETURN:
 *    SH_SUCCESS                 if all goes well
 *    SH_GENERIC_ERROR           if an error occurs
 *
 * </AUTO>
 */

static int eval_trans_quality(s_star *star_array_A, /* I: first array of stars */
int num_stars_A, /* I: number of stars in star_array_A  */
s_star *star_array_B, /* I: second array of stars */
int num_stars_B, /* I: number of stars in star_array_B  */
double star_match_radius, /* I: count as "matches" stars within this */
/*       distance, in units of list B */
TRANS *trans /* I/O: apply this trans to object in a list A */
/*       and set its fields with statistics */
/*       of the matched pairs on output */
) {
	int i, num_matched, num_possible_matches;
	int *matched_index_A, *matched_index_B;
	struct s_star *transformed_star_array_A;
	struct s_star_coord *star_coord_array_B;

#ifdef DEBUG2
	printf("entering eval_trans_quality \n");
	printf("    star_match_radius is %9.4e \n", star_match_radius);
#endif

	/* sanity checks */
	if (num_stars_A < 1) {
		shError("eval_trans_quality: given invalid num_stars_A = %d \n",
				num_stars_A);
		return (SH_GENERIC_ERROR);
	}
	if (num_stars_B < 1) {
		shError("eval_trans_quality: given invalid num_stars_B = %d \n",
				num_stars_B);
		return (SH_GENERIC_ERROR);
	}
	shAssert(star_array_A != NULL);
	shAssert(star_array_B != NULL);
	shAssert(trans != NULL);

	/*
	 * step 0: allocate space for arrays which will hold indices
	 *         of the elements of each match we find
	 */
	num_possible_matches = (
			num_stars_A > num_stars_B ? num_stars_A : num_stars_B);
	matched_index_A = shMalloc(num_possible_matches * sizeof(int));
	matched_index_B = shMalloc(num_possible_matches * sizeof(int));

	/*
	 * step 1: create a copy of star list A, which we then transform
	 *         into the coords of list B using the given TRANS
	 */
#if 0
	printf(" A: trans order is %d \n", trans->order);
#endif
	transformed_star_array_A = shMalloc(sizeof(struct s_star) * num_stars_A);
#if 0
	printf(" B: trans order is %d \n", trans->order);
#endif
	copy_star_array(star_array_A, transformed_star_array_A, num_stars_A);
#if 0
	printf(" C: trans order is %d \n", trans->order);
#endif
	if (apply_trans(transformed_star_array_A, num_stars_A, trans)
			!= SH_SUCCESS) {
		shError("eval_trans_quality: apply_trans fails \n");
		return (SH_GENERIC_ERROR);
	}

	/*
	 * step 2: create an auxiliary array of star_coord structures for
	 *         the objects in list B, and then sort that array by
	 *         'x' coordinate.  This will speed up matching in next step.
	 */
	star_coord_array_B = shMalloc(sizeof(struct s_star_coord) * num_stars_B);
	for (i = 0; i < num_stars_B; i++) {
		star_coord_array_B[i].index = i;
		star_coord_array_B[i].x = star_array_B[i].x;
		star_coord_array_B[i].y = star_array_B[i].y;
	}
	sort_star_coord_by_x(star_coord_array_B, num_stars_B);

	/*
	 * step 3: walk through the transformed stars in list A,
	 *         searching for matching stars in list B.  If we find more
	 *         than one match, keep the closest match.
	 */
	num_matched = 0;
	for (i = 0; i < num_stars_A; i++) {

		int j;
		int closest_B_index;
		//int start_index, end_index;
		int start_sc_index, end_sc_index;
		double x;
		double closest_dist_sq, star_match_radius_sq;
		struct s_star *star_A, *star_B;

		star_A = &(transformed_star_array_A[i]);
		star_match_radius_sq = star_match_radius * star_match_radius;

		/*
		 * use binary search to locate starting and ending
		 *     elements we need to check for matches in list B
		 */
#ifdef DEBUG2
		printf("   looking for match to trans A %5d  x %9.2f y %9.2f \n", i,
				star_A->x, star_A->y);
#endif
		x = star_A->x - star_match_radius;
		start_sc_index = find_star_coord_by_x(star_coord_array_B, num_stars_B,
				x);
		if (start_sc_index > 0) {
			start_sc_index--;
		}
		//start_index = star_coord_array_B[start_sc_index].index;

		x = star_A->x + star_match_radius;
		end_sc_index = find_star_coord_by_x(star_coord_array_B, num_stars_B, x);
		if (end_sc_index < num_stars_B - 1) {
			end_sc_index++;
		}
		//end_index = star_coord_array_B[end_sc_index].index;

#ifdef DEBUG3
		printf(
				"    start_index %5d %5d  x %9.4e   end_index %5d %5d x %9.4e \n",
				start_sc_index, start_index, star_array_B[start_index].x,
				end_sc_index, end_index, star_array_B[end_index].x);
#endif

		/*
		 * now go through all the stars in sorted list B which lie in the range
		 *   of possible matches to star A.  Find the star in sorted list B
		 *   which is closest to star A, and set the "closest_B_index"
		 *   value to its index in star_array_B[].
		 */
		closest_B_index = -1;
		closest_dist_sq = star_match_radius_sq * 2;
		for (j = start_sc_index; j <= end_sc_index; j++) {

			double dx, dy, distsq;

			star_B = &(star_array_B[star_coord_array_B[j].index]);
#ifdef DEBUG3
			printf("     star_B is index %5d \n", j);
#endif

			dx = fabs(star_A->x - star_B->x);
			if (dx > star_match_radius) {
#ifdef DEBUG3
				printf("      dx %11.4f > %11.4f so go to next B \n", dx,
						star_match_radius);
#endif
				continue;
			}

			dy = star_A->y - star_B->y;
			distsq = dx * dx + dy * dy;
			if (distsq < closest_dist_sq) {
#ifdef DEBUG3
				printf("      new closest distsq %11.4e \n", distsq);
#endif
				closest_dist_sq = distsq;
				closest_B_index = star_coord_array_B[j].index;
			}
		}

		/*
		 * now, check to see if the closest star was within the
		 *    star_match_radius.  If yes, we have a match;
		 *    If not, we have no match.
		 */
		if (closest_dist_sq < star_match_radius_sq) {
#ifdef DEBUG2
			printf(
					"       good match: star A %5d %9.2f %9.2f  B %5d %9.2f %9.2f  dist %11.4e \n",
					i, star_A->x, star_A->y, closest_B_index,
					star_array_B[closest_B_index].x,
					star_array_B[closest_B_index].y, sqrt(closest_dist_sq));
#endif
			matched_index_A[num_matched] = i;
			matched_index_B[num_matched] = closest_B_index;
			num_matched++;
		} else {
#ifdef DEBUG2
			printf("       no   match: star A %5d  min dist %11.4e \n", i,
					sqrt(closest_dist_sq));
#endif
		}

#ifdef DEBUG2
		printf(" at bottom of loop in eval_trans_quality, i = %5d \n", i);
		printf(" num_stars_B is %d  num_matched is %d \n", num_stars_B,
				num_matched);
		fflush(NULL);
#endif

	}

	/*
	 * now we compute the statistical properties of the matches
	 */
	{
		double sumx, sumy, sumx_sq, sumy_sq, sumtot, sumtot_sq;
		double mean_x, mean_y, stdev_x, stdev_y, mean_tot, stdev_tot;

		trans->nm = num_matched;
		sumx = 0;
		sumx_sq = 0;
		sumy = 0;
		sumy_sq = 0;
		sumtot = 0;
		sumtot_sq = 0;
		for (i = 0; i < num_matched; i++) {

			struct s_star *star_A, *star_B;
			double dx, dy, tot_sq;

			star_A = &(transformed_star_array_A[matched_index_A[i]]);
			star_B = &(star_array_B[matched_index_B[i]]);
			dx = star_A->x - star_B->x;
			dy = star_A->y - star_B->y;
#ifdef DEBUG2
			printf(
					"  item %4d  A %9.2f %9.2f  B %9.2f %9.2f  dx %9.2f dy %9.2f \n",
					i, star_A->x, star_A->y, star_B->x, star_B->y, dx, dy);
#endif
			sumx += dx;
			sumy += dy;
			sumx_sq += dx * dx;
			sumy_sq += dy * dy;
			tot_sq = (dx * dx + dy * dy);
			sumtot_sq += tot_sq;
			sumtot += sqrt(tot_sq);

		}
		if (num_matched == 1) {
			stdev_x = 0;
			stdev_y = 0;
			stdev_tot = 0;
		} else {
			mean_x = sumx / (double) num_matched;
			mean_y = sumy / (double) num_matched;
			mean_tot = sumtot / (double) num_matched;
			stdev_x = sqrt(
					(sumx_sq - num_matched * mean_x * mean_x)
							/ (num_matched - 1.0));
			stdev_y = sqrt(
					(sumy_sq - num_matched * mean_y * mean_y)
							/ (num_matched - 1.0));
			stdev_tot = sqrt(
					(sumtot_sq - num_matched * mean_tot * mean_tot)
							/ (num_matched - 1.0));
		}
		trans->sx = stdev_x;
		trans->sy = stdev_y;
		trans->sig = stdev_tot;
#ifdef DEBUG2
		printf("   num_matched %4d  sx %9.4e sy = %9.4e sig = %9.4e \n",
				num_matched, trans->sx, trans->sy, trans->sig);
#endif

	}

	/* de-allocate the arrays we've created */
	shFree(star_coord_array_B);
	shFree(transformed_star_array_A);
	shFree(matched_index_A);
	shFree(matched_index_B);

	return (SH_SUCCESS);
}

/************************************************************************
 * ROUTINE: calc_trans_sig
 *
 * DESCRIPTION:
 * Given two lists of stars, A and B, each in its own coordinate system,
 * the number of matching pairs to use, and a set of arrays which
 * define the matching members of each list, we determine the standard
 * deviation of the offsets between matching stars.  Specifically,
 *
 *     for each matching pair
 *
 *        1. transform the coords of star A into coord system of star B
 *        2. compute difference in position between the two stars
 *        3. add difference to running sums
 *
 *     and, at end of loop, compute the stdev of the differences.
 *
 * We can then set the 'sig' field of the given TRANS to this value.
 *
 * RETURN:
 *    SH_SUCCESS                 if all goes well
 *    SH_GENERIC_ERROR           if an error occurs
 *
 * </AUTO>
 */

static int calc_trans_sig(int num_matches, /* I: number of matched pairs of stars */
s_star *star_array_A, /* I: first array of stars */
int num_stars_A, /* I: number of stars in star_array_A  */
s_star *star_array_B, /* I: second array of stars */
int num_stars_B, /* I: number of stars in star_array_B  */
int *winner_votes, /* I: number of votes gotten by the top 'nbright' */
/*      matched pairs of stars */
int *winner_index_A, /* I: index into "star_array_A" of top */
/*      vote-getters */
int *winner_index_B, /* I: index into "star_array_B" of top */
/*      vote-getters */
TRANS *trans /* I/O: apply this trans to object in a list A */
/*       and set its 'sig' field when done */
) {
	int i;
	double dx, dy, dist, dist_sq;
	double new_A_x = 0.0, new_A_y = 0.0;
	double sum = 0.0, sum_sq = 0.0, mean = 0.0, stdev = 0.0;
	struct s_star *star_A, *star_B;

#ifdef DEBUG2
	printf("entering calc_trans_sig \n");
#endif

	shAssert(num_matches > 0);

	for (i = 0; i < num_matches; i++) {

		star_A = &(star_array_A[winner_index_A[i]]);
		star_B = &(star_array_B[winner_index_A[i]]);

		/* transform star A into coord system of star B */
		if (calc_trans_coords(star_A, trans, &new_A_x, &new_A_y)
				!= SH_SUCCESS) {
			shError("calc_trans_sig: calc_trans_coords fails");
			return (SH_GENERIC_ERROR);
		}

		dx = new_A_x - star_B->x;
		dy = new_A_y - star_B->y;
		dist_sq = (dx * dx + dy * dy);
		dist = sqrt(dist_sq);

		sum += dist;
		sum_sq += dist_sq;
	}

	if (num_matches == 1) {
		stdev = 0;
	} else {
		mean = sum / (double) num_matches;
		stdev = sqrt(
				(sum_sq - num_matches * mean * mean) / (num_matches - 1.0));
	}
#ifdef DEBUG
	printf("  calc_trans_sig:  num %4d  mean %9.4e stdev %9.4e \n", num_matches,
			mean, stdev);
#endif

	trans->sig = stdev;

	return (SH_SUCCESS);

}

/************************************************************************
 * ROUTINE: is_trans_good_enough
 *
 * DESCRIPTION:
 * Given a TRANS, and some parameters, determine if the TRANS
 * is good enough to meet the conditions for success.
 *
 * The current conditions are:
 *
 *    1. is the number of matching stars GREATER THAN or equal to
 *                  a given value N?
 *
 *    2. is the variance of offset(*) between stars in matching pairs
 *                  LESS than a given value S?  Note that this is
 *                  the variance = square of stdev, not the stdev itself.
 *
 *           (*) as measured in units of the list B
 *
 * If both conditions are true, we declare success; otherwise,
 * we declare failure.
 *
 * RETURN:
 *    SH_SUCCESS                 if conditions are met
 *    SH_GENERIC_ERROR           if conditions are not met
 *
 * </AUTO>
 */

static int is_trans_good_enough(int min_matches, /* I: there must be this many matched stars */
double max_stdev, /* I: stdev of offset must be less than this */
TRANS *trans /* I: check properties of this TRANS */
) {
	double variance;

#ifdef DEBUG2
	printf("entering is_trans_good_enough \n");
#endif

	shAssert(trans != NULL);
	shAssert(min_matches > 0);
	shAssert(max_stdev > 0.0);

	/*
	 * Note that (for historical reasons) we look at the variance
	 * in the offset between matched pairs of stars, not the stdev.
	 */
	variance = trans->sig * trans->sig;

	if (trans->nm >= min_matches) {
		if (variance <= max_stdev) {
			return (SH_SUCCESS);
		}
	}

	return (SH_GENERIC_ERROR);
}

/************************************************************************
 * ROUTINE: check_trans_properties
 *
 * DESCRIPTION:
 * Given a TRANS, and some parameters provided by the user,
 * verify that that the scale factor (=magnification)
 * and rotation of the TRANS match the user's desired values.
 *
 * RETURN:
 *    SH_SUCCESS                 if values match user's
 *    SH_GENERIC_ERROR           otherwise
 *
 * </AUTO>
 */

static int check_trans_properties(TRANS *trans, /* I: check properties of this TRANS */
double min_scale, /* I: scale must be at least this big ... */
double max_scale, /* I: ... and less than this value */
double rotation_deg, /* I: rotation angle must be close to this .. */
double tolerance_deg /* I: ... within this many degrees */
) {
	int scale_ok = -1;
	int rot_ok = -1;
	double scale;

#ifdef DEBUG2
	printf("entering check_trans_properties \n");
#endif

	shAssert(trans != NULL);

	/*
	 * Figure out the scale factor (=magnification) of the TRANS
	 */
	scale = sqrt(trans->b * trans->b + trans->c * trans->c);

	/*
	 * If the "min_scale" and "max_scale" values are both -1,
	 * then the user did not specify any constraints.
	 */
	if ((min_scale == -1) && (max_scale == -1)) {
		/* any value is okay */
		scale_ok = 1;
	} else {
		if ((scale < min_scale) || (scale > max_scale)) {
#ifdef DEBUG2
			printf(
					" check_trans_properties: scale %.2f not within %.2f - %.2f \n",
					scale, min_scale, max_scale);
#endif
		} else {
			scale_ok = 1;
		}
	}

	/*
	 * figure out the rotation angle, in degrees.  Force it into
	 * the range -180 to +180.
	 */
	if ((rotation_deg == AT_MATCH_NOANGLE)
			&& (tolerance_deg == AT_MATCH_NOANGLE)) {
		rot_ok = 1;
	} else {
		double min_angle_deg, max_angle_deg;
		double trans_angle_rad, trans_angle_deg;

		min_angle_deg = rotation_deg - tolerance_deg;
		max_angle_deg = rotation_deg + tolerance_deg;

		trans_angle_rad = atan2(trans->c, trans->b);
		trans_angle_deg = trans_angle_rad * (180.0 / 3.14159);

		if ((trans_angle_deg >= min_angle_deg)
				&& (trans_angle_deg <= max_angle_deg)) {
			/* okay */
			rot_ok = 1;
		} else {
#ifdef DEBUG2
			printf(
					" check_trans_properties: TRANS %.1lf deg not in range %.1f - %.1f \n",
					trans_angle_deg, min_angle_deg, max_angle_deg);
#endif
			rot_ok = -1;
		}
	}

	if ((scale_ok == 1) && (rot_ok == 1)) {
		return (SH_SUCCESS);
	} else {
		return (SH_GENERIC_ERROR);
	}
}
