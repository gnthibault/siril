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
 * FILE: match.c
 *
 * <HTML>
 * This file contains the 'main' routine, which takes instructions
 * from the user, reads information from files, calls the matching
 * routines, and writes information back onto disk.
 *
 * </HTML>
 *
 *
 *  7/18/96   - Added transonly option.  MWR
 *
 *  6/1/2000  - Added 'recalc' option; this goes through the usual steps:
 *
 *                     1. read in complete lists of stars from two sources
 *                     2. using only N brightest, find a set of matching pairs
 *                     3. calc TRANS
 *                     4. apply TRANS to all of list A
 *                     5. match all of list A (with TRANS'ed coords) against
 *                              all of list B to find _many_ pairs
 *
 *              but now it also adds the following
 *
 *                     6. using these _many_ pairs, calc TRANS again
 *
 *              The point is to avoid running the similar-triangles code
 *              on the _many_ pairs, since that's the slow part;
 *              we can very quickly calc a new, more accurate TRANS from
 *              the _many_ pairs.
 *
 *              MWR
 *
 * 6/14/2000  - Added "--version" and "--help" options, plus GPL.
 *              MWR
 *
 * 7/30/2000 - Added "id1" and "id2" options.  MWR
 *
 * 1/21/2001 - Added "max_iter" and "halt_sigma" options, plus some sanity
 *             checks on values of user-supplied arguments.
 *             MWR
 *
 * 12/14/2001- Added additional output values in TRANS which describe
 *             the quality of the fit.   Added optional arguments which
 *             cause calculation of additional statistics.
 *             Thanks to John Blakeslee.
 *             MWR
 *
 * 12/31/2001- Added options:
 *                 allow user to specify an input TRANS
 *                 allow user to specify a pure offset between coords
 *             Thanks to John Blakeslee.
 *             MWR
 *
 * 11/23/2002- Fixed bug in "reset_A_coords()" which caused an assert
 *                 failure if given empty lists.  Now the routine
 *                 simply prints a warning message and returns
 *                 with normal code.
 *             MWR
 *
 * 11/4/2011 - Added new options "rotangle=" and "rottol=",
 *                 so that user may specify a range of allowed rotations
 *                 between the two lists of stars.
 *             MWR
 *
 * </AUTO>
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "core/siril.h"
#include "gui/callbacks.h"
#include "algos/PSF.h"
#include "registration/matching/misc.h"
#include "registration/matching/match.h"
#include "registration/matching/atpmatch.h"

#undef DEBUG           /* get some of diagnostic output */

static void reset_copy_ids(int numA, struct s_star *star_list_A,
		struct s_star *star_list_A_copy);

static int reset_A_coords(int numA, struct s_star *post_list_A,
		struct s_star *pre_list_A);
static int prepare_to_recalc(char *outfile, int *num_matched_A,
		struct s_star **matched_list_A, int *num_matched_B,
		struct s_star **matched_list_B, struct s_star *star_list_A_copy,
		TRANS *trans);

char progname[CMDBUFLEN + 1];

int star_match(fitted_PSF **s1, fitted_PSF **s2, int n, TRANS *t) {
	int ret;
	int numA, numB;
	int num_matched_A, num_matched_B;
	int numA_copy;
	int max_iter = AT_MATCH_MAXITER;
	int min_req_pairs = AT_MATCH_MINPAIRS;
	int trans_order = AT_TRANS_LINEAR;
	double triangle_radius = AT_TRIANGLE_RADIUS; /* in triangle-space coords */
	double match_radius = AT_MATCH_RADIUS; /* in units of list B */
	double scale = -1.0;
	double min_scale = 0.99;
	double max_scale = 1.01;
	double rot_angle = AT_MATCH_NOANGLE; /* by default, any angle is okay */
	double rot_tol = AT_MATCH_NOANGLE;
	double halt_sigma = AT_MATCH_HALTSIGMA;
	double medsigclip = 0.0; /* used in MEDTF calcs */
	double xshift = 0.0; /* guessed shift in X dir */
	double yshift = 0.0; /* guessed shift in y dir */
	int nobj = AT_MATCH_NBRIGHT;
	int transonly = 0; /* if 1, only find TRANS */
	int recalc = 1; /* if 1, calc TRANS again */
	int num_matches = 0; /* number of matching pairs */
	int medtf_flag = 0; /* calculate MEDTF stats? */
	int intrans = 0; /* use given input TRANS? */
	int identity = 0; /* use identity as TRANS? */
	char outfile[CMDBUFLEN + 1];
	char intransfile[CMDBUFLEN + 1];
	const char *tmpdir = g_get_tmp_dir();
	struct s_star *star_list_A, *star_list_B;
	struct s_star *star_list_A_copy;
	struct s_star *matched_list_A, *matched_list_B;
	MEDTF *medtf;
	TRANS *trans;

	/* buffer overflow paranoia */
	outfile[CMDBUFLEN] = '\0';
	progname[CMDBUFLEN] = '\0';
	intransfile[CMDBUFLEN] = '\0';

	strncpy(outfile, tmpdir, CMDBUFLEN);
	strcat(outfile, "/matched");

	/*
	 * We can only calculate _clipped_ MEDTF statistics if we calculate
	 *   the regular ones.  So it makes no sense to specify 'medsigclip',
	 *   but not 'medtf'.
	 */
	if ((medtf_flag == 0) && (medsigclip != 0.0)) {
		shError("WARNING: medsigclip requires the 'medtf' option");
		shError("WARNING: Setting medtf ");
		medtf_flag = 1;
	}
	/*
	 * Impossible to calculation the MEDTF statistics unless we match
	 *   up stars after finding the TRANS.  So, if the user does want
	 *   the MEDTF stats, we force the necessary step.
	 */
	if ((medtf_flag == 1) && (recalc == 0)) {
		shError("WARNING: 'medtf' requires the 'recalc' option");
		shError("WARNING: Setting recalc ");
		recalc = 1;
	}

	/*
	 * Check to make sure that the user did exactly one of the following:
	 *    a. did not specify "scale" or "min_scale" or "max_scale"
	 *    b. did specify "scale" only
	 *    c. did specify "min_scale" and "max_scale", but not "scale"
	 *
	 * If choice b., then translate the single "scale" value into
	 * a pair of "min_scale" and "max_scale" limits.  We'll always
	 * pass these limits to the matching procedures.
	 */
	if ((scale == -1) && (min_scale == -1) && (max_scale == -1)) {
		/* okay */
		;
	} else if ((scale != -1) && (min_scale == -1) && (max_scale == -1)) {
		/* okay */
		min_scale = scale - (0.01 * AT_MATCH_PERCENT * scale);
		max_scale = scale + (0.01 * AT_MATCH_PERCENT * scale);
	} else if ((scale == -1) && (min_scale != -1) && (max_scale != -1)) {
		/* okay */
		if (min_scale > max_scale) {
			shFatal("min_scale must be smaller than max_scale");
			return (SH_GENERIC_ERROR);
		}
	} else {
		/* not okay */
		shFatal("invalid combination of 'scale', 'min_scale', 'max_scale'");
		return (SH_GENERIC_ERROR);
	}
#ifdef DEBUG
	if ((scale == -1) && (min_scale == -1) && (max_scale == -1)) {
		printf("No limits set on relative scales for matching. \n");
	} else {
		printf("using min_scale %f  max_scale %f \n", min_scale, max_scale);
	}
#endif

	/*
	 * Check to make sure that the user specified
	 *
	 *       a)   neither "rotangle" nor "rottol"
	 *   or
	 *       b)   both "rotangle" and "rottol"
	 */
	if ((rot_angle == AT_MATCH_NOANGLE) && (rot_tol == AT_MATCH_NOANGLE)) {
		/* this is okay */
	} else if ((rot_angle != AT_MATCH_NOANGLE)
			&& (rot_tol != AT_MATCH_NOANGLE)) {
		/* this is okay */
	} else {
		/* this is NOT okay */
		shFatal("Must specify both 'rotangle' and 'rottol', or neither ");
		return (SH_GENERIC_ERROR);
	}

	/* Can only specify one of 'identity' and 'intrans' */
	if ((intrans != 0) && (identity != 0)) {
		shFatal("Cannot specify both 'identity' and an input TRANS file");
		return (SH_GENERIC_ERROR);
	}

	/* Makes no sense to specify 'identity' and 'quadratic' or 'cubic' */
	if ((identity != 0) && (trans_order != AT_TRANS_LINEAR)) {
		shFatal("Cannot specify both 'identity' and any order beyond linear");
		return (SH_GENERIC_ERROR);
	}
	/*
	 * Likewise, makes no sense to specify 'intrans=' and
	 * 'quadratic' or 'cubic' (or 'linear', actually, but we have
	 * no easy way to check for that mistake :-(
	 */
	if ((intrans != 0) && (trans_order != AT_TRANS_LINEAR)) {
		shFatal("Cannot specify both 'intrans' and any order ");
		return (SH_GENERIC_ERROR);
	}

	/*
	 * The minimum number of pairs required for a successful match
	 * must be at least large enough to define the given transformation.
	 * If the user has specified a value for 'min_req_pairs'
	 * which is LESS than this minimum value for the desired TRANS,
	 * print an error and terminate.
	 */
	if ((trans_order == AT_TRANS_LINEAR) && (min_req_pairs < 3)) {
		shFatal("linear trans requires at least 3 matched pairs ");
		return (SH_GENERIC_ERROR);
	}
	if ((trans_order == AT_TRANS_QUADRATIC) && (min_req_pairs < 6)) {
		shFatal("quadratic trans requires at least 6 matched pairs ");
		return (SH_GENERIC_ERROR);
	}
	if ((trans_order == AT_TRANS_CUBIC) && (min_req_pairs < 10)) {
		shFatal("cubic trans requires at least 10 matched pairs ");
		return (SH_GENERIC_ERROR);
	}

	/*
	 * There are three possibilities for the initial TRANS values:
	 *
	 *   1. The user specified "identity", so we start with a linear TRANS
	 *      which has no scale change or rotation; it has no translation,
	 *      unless the user gave 'xsh' and 'ysh'.
	 *
	 *   2. The user specified "intrans=file", in which case we read
	 *      the initial TRANS values from the given file.
	 *
	 *   3. The user didn't give us either, so we start with an empty
	 *      TRANS and will have to find one ourselves via atFindTrans().
	 *
	 * In either case 1 or 2, we will end up with a non-zero initial
	 * TRANS structure -- which means that we want to use IT to find
	 * matches, instead of atFindTrans().  We'll therefore set the 'intrans'
	 * flag to 1, to indicate "don't try to find a TRANS, just use
	 * the one the user provided".
	 */
	if (identity == 1) {
		trans = getIdentityTrans();
		trans_order = trans->order;
		if (xshift != 0.0) {
			trans->a = xshift;
		}
		if (yshift != 0.0) {
			trans->d = yshift;
		}

		intrans = 1;

	} else if (intrans == 1) {
		if ((trans = getGuessTrans(intransfile)) == NULL) {
			shFatal("GetGuessTrans fails -- must give up");
			return (SH_GENERIC_ERROR);
		}
		trans_order = trans->order;
	} else {
		/* this will be an "empty" TRANS; atFindTrans will try to fill it */
		atTransOrderSet(trans_order);
		trans = atTransNew();
		trans->order = trans_order;
	}
#ifdef DEBUG
	printf("using trans_order %d\n", trans_order);
#endif

	/* read information from the first file */
	/*	if (read_star_file(fileA, x1col, y1col, mag1col, id1col, -1, &numA,
	 &star_list_A) != SH_SUCCESS) {
	 shError("can't read data from file %s", fileA);
	 exit(1);
	 }*/
	if (get_stars(s1, n, &numA, &star_list_A)) {
		printf("can't read data\n");
		return (SH_GENERIC_ERROR);
	}

	/*
	 * We always (whether given initial TRANS or not) will
	 *   make a second calculation of the TRANS, using only objects in
	 *   the set of matched pairs.
	 *   As input to this second calculation, we will need items
	 *   from list A with their original coordinates.
	 *   Therefore, we now create a copy of the stars in set A,
	 *   so that we can restore the output matched coords
	 *   (which have been converted to those in set B) with the original coords.
	 */
	/*	if (read_star_file(fileA, x1col, y1col, mag1col, id1col, -1, &numA_copy,
	 &star_list_A_copy) != SH_SUCCESS) {
	 shFatal("can't read data from file %s", fileA);
	 }*/
	if (get_stars(s1, n, &numA_copy, &star_list_A_copy)) {
		printf("can't read data\n");
		return (SH_GENERIC_ERROR);
	}
	/* sanity check */
	shAssert(numA_copy == numA);

	/*
	 * reset the 'id' field values in the star_list_A_copy so that they
	 * match the 'id' field values in their counterparts in star_list_A
	 */
	reset_copy_ids(numA, star_list_A, star_list_A_copy);

	/* read information from the second file */
	/*	if (read_star_file(fileB, x2col, y2col, mag2col, id2col, -1, &numB,
	 &star_list_B) != SH_SUCCESS) {
	 shError("can't read data from file %s", fileB);
	 exit(1);
	 }*/
	if (get_stars(s2, n, &numB, &star_list_B)) {
		printf("can't read data\n");
		return (SH_GENERIC_ERROR);
	}

	/*
	 * Now, if the has has not given us an initial TRANS structure, we need
	 * to find one ourselves.
	 */
	if (intrans == 0) {
		ret = atFindTrans(numA, star_list_A, numB, star_list_B, match_radius,
				triangle_radius, nobj, min_scale, max_scale, rot_angle, rot_tol,
				max_iter, halt_sigma, min_req_pairs, trans);
		if (ret != SH_SUCCESS) {
			shFatal("initial call to atFindTrans fails");
			/** */
			atTransDel(trans);
			free_star_array(star_list_A);
			free_star_array(star_list_A_copy);
			free_star_array(star_list_B);
			/** */
			return (SH_GENERIC_ERROR);
		}
	}

#ifdef DEBUG
	printf("using trans_order %d.\n", trans_order);
	printf("Initial trans structure:\n");
	print_trans(trans);
#endif

	/*
	 * if the "transonly" flag is set, stop at this point
	 */
	if (transonly == 1) {
		print_trans(trans);
		return (0);
	}

	/*
	 * having found (or been given) the TRANS that takes A -> B, let us apply
	 * it to all the elements in A; thus, we'll have two sets of
	 * of stars, each in the same coordinate system
	 */
	atApplyTrans(numA, star_list_A, trans);

	/*
	 * now match up the two sets of items, and find four subsets:
	 *
	 *     those from list A that do     have matches in list B
	 *     those from list B that do     have matches in list A
	 *     those from list A that do NOT have matches in list B
	 *     those from list B that do NOT have matches in list A
	 *
	 * We may use the two sets of matched objects for further processing,
	 * so we put the names of the files containing those matched objects
	 * into 'matched_file_A' and 'matched_file_B' for easy reference.
	 */
	atMatchLists(numA, star_list_A, numB, star_list_B, match_radius, outfile,
			&num_matches);
	trans->nm = num_matches;

	/*
	 * Okay, there are two possibilities at this point:
	 *
	 *    1. The user gave us information about the initial TRANS,
	 *       either via 'intrans=' or 'identity'.  We have applied
	 *       his initial TRANS to the input list A, and looked for
	 *       matched pairs.
	 *
	 *    2. The user didn't give us any information about an initial
	 *       TRANS, so we called 'atFindTrans()' to find one.
	 *       We have applied this TRANS to input list A, and
	 *       looked for matched pairs.
	 *
	 * Now, we want to improve the initial TRANS, whether it was supplied
	 * by user or determined by 'atFindTrans()'.  We do so by applying
	 * the initial TRANS to only the matched objects in list A, and then
	 * calling 'atRecalcTrans()' on the matched objects only.  This should
	 * give an improved TRANS, since it very likely won't be contaminated
	 * by any spurious matches.
	 */

	/* need to send trans to prepare_to_recalc because it adds sdx,sdy to it */
	if (prepare_to_recalc(outfile, &num_matched_A, &matched_list_A,
			&num_matched_B, &matched_list_B, star_list_A_copy, trans) != 0) {
		shFatal("prepare_to_recalc fails");
		/** */
		atTransDel(trans);
		free_star_array(matched_list_A);
		free_star_array(matched_list_B);
		free_star_array(star_list_A_copy);
		/** */
		return (SH_GENERIC_ERROR);
	}
	/* okay, now we're ready to call atRecalcTrans, on matched items only */
	if (atRecalcTrans(num_matched_A, matched_list_A, num_matched_B,
			matched_list_B, max_iter, halt_sigma, trans) != SH_SUCCESS) {
		shFatal("atRecalcTrans fails on matched pairs only");
		/** */
		atTransDel(trans);
		free_star_array(matched_list_A);
		free_star_array(matched_list_B);
		free_star_array(star_list_A_copy);
		/** */
		return (SH_GENERIC_ERROR);
	}
#ifdef DEBUG
	printf("TRANS based on matches only :\n");
	print_trans(trans);
#endif

	/*
	 * At this point, we have a TRANS which is based solely on those items
	 * which matched.  If the user wishes, we can improve the TRANS
	 * even more by applying the current transformation to ALL items
	 * in list A, making a second round of matching pairs, and then
	 * using these pairs to calculate a new and better TRANS.
	 *
	 * The point is that we'll probably end up with more matched pairs
	 * if we start with the current TRANS, instead of the initial TRANS.
	 */
	if (recalc == 1) {

		/* re-set coords of all items in star A */
		if (reset_A_coords(numA, star_list_A, star_list_A_copy) != 0) {
			shFatal("reset_A_coords returns with error before recalc");
			/** */
			atTransDel(trans);
			free_star_array(matched_list_A);
			free_star_array(matched_list_B);
			free_star_array(star_list_A_copy);
			/** */
			return (SH_GENERIC_ERROR);
		}

		/*
		 * apply the current TRANS (which is probably much better than
		 * the initial TRANS) to all items in list A
		 */
		atApplyTrans(numA, star_list_A, trans);

		/*
		 * Match items in list A to those in list B
		 */
		atMatchLists(numA, star_list_A, numB, star_list_B, match_radius,
				outfile, &num_matches);
		trans->nm = num_matches;
#ifdef DEBUG
		printf("After tuning with recalc, num matches is %d\n", num_matches);
		print_trans(trans);
#endif

		/* prepare to call atRecalcTrans one last time */
		/* need to send trans to prepare_to_recalc because it adds sdx,sdy */
		if (prepare_to_recalc(outfile, &num_matched_A, &matched_list_A,
				&num_matched_B, &matched_list_B, star_list_A_copy, trans)
				!= 0) {
			shFatal("prepare_to_recalc fails");
			/** */
			atTransDel(trans);
			free_star_array(matched_list_A);
			free_star_array(matched_list_B);
			free_star_array(star_list_A_copy);
			/** */
			return (SH_GENERIC_ERROR);
		}

		/* final call atRecalcTrans, on matched items only */
		if (atRecalcTrans(num_matched_A, matched_list_A, num_matched_B,
				matched_list_B, max_iter, halt_sigma, trans) != SH_SUCCESS) {
			shFatal("atRecalcTrans fails on matched pairs only");
			/** */
			atTransDel(trans);
			free_star_array(matched_list_A);
			free_star_array(matched_list_B);
			free_star_array(star_list_A_copy);
			/** */
			return (SH_GENERIC_ERROR);
		}

#ifdef DEBUG
		printf("TRANS based on recalculated matches is \n");
		print_trans(trans);
#endif
	}

	/*
	 * If the user wants summary statistics on the straight translation
	 * differences between the lists (without any scale or rotation),
	 * we calculate them now.  We take all the pairs of stars which
	 * match -- based on their TRANSFORMED coordinates -- and go back
	 * to look at the differences between their ORIGINAL coordinates.
	 *
	 * This only makes sense if a pure translation connects the lists.
	 */
	if (medtf_flag) {
		if (reset_A_coords(num_matched_A, matched_list_A, star_list_A_copy)
				!= 0) {
			shFatal("second call to reset_A_coords returns with error");
			/** */
			atTransDel(trans);
			free_star_array(matched_list_A);
			free_star_array(matched_list_B);
			free_star_array(star_list_A_copy);
			/** */
			return (SH_GENERIC_ERROR);
		}

		medtf = atMedtfNew();
		if (atFindMedtf(num_matched_A, matched_list_A, num_matched_B,
				matched_list_B, medsigclip, medtf) != SH_SUCCESS) {
			shError("atFindMedtf fails");
		} else {
			print_medtf(medtf);
		}
		atMedtfDel(medtf);
	}

	print_trans(trans);
	*t = *trans;

	return (0);
}

/***********************************************************************
 * FUNCTION: reset_copy_ids
 *
 * Modify the 'id' field values in the given list (a copy of list A)
 *   so that they will match the 'id' values in the corresponding
 *   stars of list A.
 *
 * We have to do this because the routine which creates new s_star
 *   structs keeps incrementing the 'id' values, and so the copies
 *   have a different value.
 *
 * RETURNS
 *   nothing
 */

static void reset_copy_ids(int numA, /* I: number of stars in list A */
struct s_star *star_list_A, /* I: original star A list */
struct s_star *star_list_A_copy /* I/O: copy of star A list */
) {
	int i;
	struct s_star *star, *star_copy;

	star = star_list_A;
	star_copy = star_list_A_copy;
	for (i = 0; i < numA; i++) {
		shAssert(star != NULL);
		shAssert(star_copy != NULL);
		star_copy->id = star->id;

		star = star->next;
		star_copy = star_copy->next;
	}
}

/***********************************************************************
 * FUNCTION: reset_A_coords
 *
 * Given the number of elements in list 'post_list_A',
 *   and two versions of the
 *   stars in list A (after conversion to coord system of list B,
 *   and before conversion), restore the original coords of stars
 *   in the the matched list.
 *
 * RETURNS
 *   0            if all goes well
 *   1            otherwise
 */

static int reset_A_coords(int numA, /* I: number of stars in the list */
struct s_star *post_list_A, /* I/O: stars in A, after they have */
/*         been matched to stars in  */
/*         list B; coords have been */
/*         converted.  We'll reset  */
/*         coords in this list */
struct s_star *pre_list_A /* I: stars in A, with original coords */
) {
	int post_index;
	int found_it;
	struct s_star *pre_star, *post_star;

	/* if handed empty list, do nothing. */
	if (numA == 0) {
		shError("reset_A_coords: handed empty list, will do nothing");
		return (0);
	}

	/* sanity checks */
	shAssert(post_list_A != NULL);
	shAssert(pre_list_A != NULL);

	for (post_index = 0, post_star = post_list_A; post_index < numA;
			post_index++, post_star = post_star->next) {

		shAssert(post_star != NULL);

		found_it = 0;
		pre_star = pre_list_A;
		while (pre_star != NULL) {

			if (pre_star->id == post_star->id) {
				post_star->x = pre_star->x;
				post_star->y = pre_star->y;
				found_it = 1;
				break;
			}
			pre_star = pre_star->next;
		}
		if (found_it == 0) {
			printf("reset_A_coords: no match for post_star %d?\n", post_index);
			shAssert(0);
		}
	}

	return (0);
}

/**********************************************************************
 * PROCEDURE: prepare_to_recalc
 *
 * DESCRIPTION: This function sets us up to call "atRecalcTrans".
 *              We have already found (or been given) a TRANS, and
 *              used it to match up items from list A and list B.
 *              Those matched items are in a pair of files
 *              with names based on 'outfile', but with
 *              extensions "mtA" and "mtB".
 *              Ex: if 'outfile' is "matched",
 *              then the two sets of matched items are in
 *
 *                    matched.mtA      matched.mtB
 *
 *              The format of each file is one star per line,
 *              with 4 fields:
 *
 *                    internal_ID     xval   yval     mag
 *
 *              where the coords (xval, yval) are in system of list B.
 *
 *              We are about to use these good, matched items to
 *              find an improved TRANS -- which should take objects
 *              from coord system A to coord system B.
 *
 *              In order to do that, we need to
 *
 *                  a. Read in the good, matched items.  The coordinates
 *                     of objects in list A will have been transformed
 *                     to their corresponding values in coord system
 *                     of list B, so ...
 *
 *                  b. We must then re-set the coords of the items in
 *                     list A to their original values, so that we can
 *                     re-calculate a TRANS which takes the coords
 *                     from system A to system B.
 *
 *              We also take this opportunity to compare the transformed
 *              positions of items in list A against the positions of
 *              the matching objects in list B.  We calculate the
 *              RMS of the differences in both "x" and "y" directions,
 *              and place them into the "sx" and "sy" members of
 *              the current TRANS.
 *
 *
 * RETURNS:
 *    0             if all goes well
 *    1             if there's an error
 */

static int prepare_to_recalc(char *outfile, /* I: stem of files with matched items */
int *num_matched_A, /* O: number of stars in matched set */
/*      from list A */
struct s_star **matched_list_A, /* O: fill this with matched items from */
/*      list A, in coord system B */
int *num_matched_B, /* O: number of stars in matched set */
/*      from list B */
struct s_star **matched_list_B, /* O: fill this with matched items from */
/*      list B, in coord system B */
struct s_star *star_list_A_copy, /* O: fill this with matched items from */
/*      list A, with their orig coords  */
TRANS *trans /* O: we calc herein the sx, sy fields  */
/*      so put them into this TRANS */
) {
	char matched_file_A[CMDBUFLEN + 4];
	char matched_file_B[CMDBUFLEN + 4];
	double Xrms, Yrms;

	shAssert(outfile != NULL);

	sprintf(matched_file_A, "%s.mtA", outfile);
	if (read_matched_file(matched_file_A, num_matched_A,
			matched_list_A) != SH_SUCCESS) {
		shError("read_matched_file can't read data from file %s",
				matched_file_A);
		return (1);
	}
	sprintf(matched_file_B, "%s.mtB", outfile);
	if (read_matched_file(matched_file_B, num_matched_B,
			matched_list_B) != SH_SUCCESS) {
		shError("read_matched_file can't read data from file %s",
				matched_file_B);
		return (1);
	}

	/* here we find the rms of those stars we read in -- JPB 17/Jan/02 */
	if (atCalcRMS(*num_matched_A, *matched_list_A, *num_matched_B,
			*matched_list_B, &Xrms, &Yrms) != SH_SUCCESS) {
		shFatal("atCalcRMS fails on matched pairs");
	}
	trans->sx = Xrms;
	trans->sy = Yrms;
	/************************************************/

	if (reset_A_coords(*num_matched_A, *matched_list_A, star_list_A_copy)
			!= 0) {
		shError("prepare_to_recalc: reset_A_coords returns with error");
		return (1);
	}

	return (0);
}
