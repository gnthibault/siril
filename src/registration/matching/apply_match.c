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
 * FILE: apply_match.c
 *
 * <HTML>
 * Given 
 *    - an ASCII file consisting of a list of stars, with one line per
 *        star and multiple columns of information separated by white space
 *    - the numbers of the columns containing "X" and "Y" coords of stars
 *    - a central RA and Dec, each in decimal degrees
 *    - coefficients of a TRANS structure which convert "X" and "Y"
 *        coordinates from the ASCII file's system to plate coordinates
 *        (xi, eta), which are projections onto the tangent plane 
 *        centered on the central RA and Dec.
 *
 * run through the data file.  For each entry, calculate the (RA, Dec) of
 * the star, then replace the "X" and "Y" values with (RA, Dec).  Leave all
 * other information in the ASCII file as-is.
 *
 * The TRANS structure converts "X" and "Y" to (xi, eta) in one of three
 * ways.  If the user specifies a 'linear' transformation (which is the
 * default), then
 *
 *     xi = A + Bx + Cy
 *    eta = D + Ex + Fy
 *
 * In the case of 'quadratic', 
 *
 *     xi =  A + Bx + Cy + Dxx + Exy + Fyy
 *    eta =  G + Hx + Iy + Jxx + Kxy + Lyy
 *
 * In the case of 'cubic', 
 *
 *     xi =  A + Bx + Cy + Dxx + Exy + Fyy + Gx(xx+yy) + Hy(xx+yy)
 *    eta =  I + Jx + Ky + Lxx + Mxy + Nyy + Ox(xx+yy) + Py(xx+yy)
 *
 * where "xi" and "eta" are in radians, and measure the distance
 * of the star from the center of field from which the TRANS was calculated.
 * We assume that the given "ra" and "dec" values are the same as this
 * central position.
 *
 * We force all values of RA to lie between 0 < RA < 360
 *
 * Print the results to stdout, or place them into the file given
 * by the optional "outfile" command-line argument.
 *
 * Usage: apply_match starfile1 xcol ycol ra dec linear|quadratic|cubic
 *                    a b c d e f [g h i j k [l m n o ]] [outfile=] 
 *
 * </HTML>
 * </AUTO>
 *
 * modified to assume that the TRANS structure takes (x, y) coords and
 *   turns them into (xi, eta), in radians, rather than in arcseconds.
 *   MWR 5/24/2000
 *
 * modified to handle the three cases of linear, quadratic, or cubic
 *   TRANSformations. 
 *   MWR 6/11/2000
 *
 * fixed equations in proc_star_file() so that they handle properly
 *   the coordinate transformations near the celestial poles.
 *   MWR 5/19/2003
 *
 * added 10 more %s in the "sscanf" statement in proc_star_file()
 *   MWR 6/11/2008
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_LIBCURL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "core/siril.h"
#include "algos/PSF.h"
#include "algos/plateSolver.h"
#include "apply_match.h"
#include "misc.h"
#include "degtorad.h"

#undef DEBUG           /* get some of diagnostic output */

static int
proc_star_file(fitted_PSF **s, image_solved *solved, TRANS *trans);

int apply_match(image_solved *solved, TRANS trans) {

	/* now walk through the file and do the dirty work */
	if (proc_star_file(com.stars, solved, &trans) != SH_SUCCESS) {
		shError("can't process data for platesolving");
		return (1);
	}

	return (0);
}

/****************************************************************************
 * ROUTINE: proc_star_file
 *
 * walk through the given file, one line at a time.  
 *
 * If the line starts with COMMENT_CHAR, place it into the output stream.
 * If the line is completely blank, place it into the output stream.
 *
 * Otherwise, 
 *   - read in the entire line, 
 *   - figure out the "X" and "Y" coords
 *   - transform the "X" and "Y" coords to be (RA, Dec) from the central
 *         "ra" and "dec", in units of arcseconds
 *   - transform from the tangent plane back to the spherical sky, so that
 *         we have genuine (RA, Dec) for each star
 *   - print out the line, replacing the "X" with RA, the "Y" with Dec
 *
 * RETURNS:
 *   SH_SUCCESS            if all goes well
 *   SH_GENERIC_ERROR      if not
 */

static int proc_star_file(fitted_PSF **s, /* I: name of input data with star list */
image_solved *solved, TRANS *trans /* I: TRANS taking (x,y) -> (ra, dec) */
) {
//	int i = 0;
	double xval, yval;
	double r_dec;
	double z, alpha, delta;
	double delta_ra, delta_dec;
	double rsquared;
	double ra = solved->px_cat_center.x;
	double dec = solved->px_cat_center.y;

	r_dec = dec * DEGTORAD;

//	while (s && s[i]) {
	xval = solved->x; //s[i]->xpos;
	yval = solved->y; //s[i]->ypos;
	/*
	 * let's transform from (x,y) to (delta_ra, delta_dec),
	 * using either a linear, quadratic, or cubic transformation
	 * (signalled by the 'order' field of the TRANS)
	 */
	switch (trans->order) {
	default:
	case AT_TRANS_LINEAR:
		delta_ra = trans->a + trans->b * xval + trans->c * yval;
		delta_dec = trans->d + trans->e * xval + trans->f * yval;
		break;
	case AT_TRANS_QUADRATIC:
		delta_ra = trans->a + trans->b * xval + trans->c * yval
				+ trans->d * xval * xval + trans->e * xval * yval
				+ trans->f * yval * yval;
		delta_dec = trans->g + trans->h * xval + trans->i * yval
				+ trans->j * xval * xval + trans->k * xval * yval
				+ trans->l * yval * yval;
		break;
	case AT_TRANS_CUBIC:
		rsquared = xval * xval + yval * yval;
		delta_ra = trans->a + trans->b * xval + trans->c * yval
				+ trans->d * xval * xval + trans->e * xval * yval
				+ trans->f * yval * yval + trans->g * xval * rsquared
				+ trans->h * yval * rsquared;
		delta_dec = trans->i + trans->j * xval + trans->k * yval
				+ trans->l * xval * xval + trans->m * xval * yval
				+ trans->n * yval * yval + trans->o * xval * rsquared
				+ trans->p * yval * rsquared;
		break;
	}

#if 1
	/*
	 * and now convert from arcseconds to radians
	 * (convenient for calculations)
	 */
	delta_ra = (delta_ra / 3600.0) * DEGTORAD;
	delta_dec = (delta_dec / 3600.0) * DEGTORAD;
#endif

	/*
	 * we have (delta_ra, delta_dec), in radians; these give the distance
	 * of this star from the central (RA, Dec).  Now we can de-project from
	 * the tangent plane (centered on RA,Dec) and calculate the actual
	 * RA, Dec of the star (in degrees)
	 */
	{
		double zz;

		z = cos(r_dec) - delta_dec * sin(r_dec);
		zz = atan2(delta_ra, z) / DEGTORAD;
		alpha = zz + ra;

		zz = cos((alpha - ra) * DEGTORAD)
				* (sin(r_dec) + delta_dec * cos(r_dec));
		delta = atan2(zz, z) / DEGTORAD;
	}

	/*
	 * make sure new RA lies in range  0 < RA < 360
	 */
	if (alpha < 0) {
		alpha += 360.0;
	}
	if (alpha >= 360.0) {
		alpha -= 360.0;
	}

	/*
	 * make sure Dec lies in range  -90 < Dec < +90
	 */
	if (delta < -90) {
		delta += 180;
	}
	if (delta > 90) {
		delta -= 180;
	}
#ifdef DEBUG
	fprintf(stdout, "new RA = %10.5f, new dec = %10.5f\n", alpha, delta);
#endif
	solved->ra = alpha;
	solved->dec = delta;

//	i++;
//	}

	return (SH_SUCCESS);
}
#endif
