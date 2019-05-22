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
 *
 *  Changed the format of the "print_trans()" output so that TRANS elements
 *    are now printed as %15.9e instead of %15.9f.  This will make a big
 *    difference in accuracy when the elements have small values -- as they
 *    will when one of the coordinate systems is in radians and the field
 *    size is that of a typical CCD (10 arminutes).
 *    June 26, 2010
 *  Michael Richmond
 */

/*
 * little support functions for the matching code
 *
 */

#include "core/siril.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include "gui/callbacks.h"
#include "algos/PSF.h"
#include "registration/matching/misc.h"
#include "registration/matching/atpmatch.h"

#undef DEBUG

/*
 * this variable holds the value of the TRANS order used in this
 * instance of the program.  It signals whether we're using
 * linear, quadratic, or cubic terms in the transformation.
 *
 * See
 *       atTransOrderSet
 *       atTransOrderGet
 */
static int at_trans_order = -1;

/*********************************************************************
 * ROUTINE: shMalloc
 *
 * Attempt to allocate the given number of bytes.  Return the
 * memory, if we succeeded, or print an error message and
 * exit with error code if we failed.
 *
 * RETURNS:
 *      void *             to new memory, if we got it
 */

void *
shMalloc(int nbytes /* I: allocate a chunk of this many bytes */
) {
	void *vptr;

	if ((vptr = (void *) malloc(nbytes)) == NULL) {
		shError("shMalloc: failed to allocate for %d bytes", nbytes);
		exit(1);
	}
	return (vptr);
}

/*********************************************************************
 * ROUTINE: shFree
 *
 * Attempt to free the given piece of memory.
 *
 * RETURNS:
 *      nothing
 */

void shFree(void *vptr /* I: free this chunk of memory */
) {
	free(vptr);
}


/*********************************************************************
 * ROUTINE: shError
 *
 * Print the given error message to stderr, but continue to execute.
 *
 * RETURNS:
 *      nothing
 */

void shError(char *format, /* I: format part of printf statement */
... /* I: optional arguments to printf */
) {
	va_list ap;

	va_start(ap, format);
	(void) vfprintf(stderr, (const char *) format, ap);
	fputc('\n', stderr);
	fflush(stdout);
	fflush(stderr);
	va_end(ap);
}

/*********************************************************************
 * ROUTINE: shFatal
 *
 * Print the given error message to stderr, and halt program execution.
 *
 * RETURNS:
 *      nothing
 */

void shFatal(char *format, /* I: format part of printf statement */
... /* I: optional arguments to printf */
) {
	va_list ap;

	va_start(ap, format);
	(void) vfprintf(stderr, (const char *) format, ap);
	fputc('\n', stderr);
	fflush(stdout);
	fflush(stderr);
	va_end(ap);
	//exit(1);
}

/*********************************************************************
 * ROUTINE: shDebugSet, shDebug
 *
 * shDebugSet: sets the static variable 'debug_level' to the
 * given integer; it starts at 0, but the user may set it to
 * higher levels.  The higher the level, the more messages
 * may be printed during execution.
 *
 * shDebug: If the current 'debug level' is >= the passed 'level',
 * then print the given message to stdout, and continue execution.
 * Otherwise, just continue.
 *
 * RETURNS:
 *      nothing
 */

static int debug_level = 0;

void shDebug(int level, /* I: debug level at which we print this */
char *format, /* I: format part of printf statement */
... /* I: optional arguments to printf */
) {
	va_list ap;

	if (level > debug_level) {
		return;
	}

	va_start(ap, format);
	(void) vfprintf(stdout, (const char *) format, ap);
	fputc('\n', stdout);
	fflush(stdout);
	va_end(ap);
}

/************************************************************************
 * ROUTINE: atTransOrderSet
 *
 * DESCRIPTION:
 * Set the value of the order we'll use for TRANS structures.
 * Possibilities are:
 *
 *      AT_TRANS_LINEAR      linear transformation
 *      AT_TRANS_QUADRATIC   linear plus quadratic terms
 *      AT_TRANS_CUBIC       linear plus quadratic plus cubic terms
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

void atTransOrderSet(int order /* I: order for all TRANS structures */
) {
	at_trans_order = order;
}

/************************************************************************
 * ROUTINE: atTransOrderGet
 *
 * DESCRIPTION:
 * Get the value of the order we're using in this instance of the program.
 * Possibilities are:
 *
 *      AT_TRANS_LINEAR      linear transformation
 *      AT_TRANS_QUADRATIC   linear plus quadratic terms
 *      AT_TRANS_CUBIC       linear plus quadratic plus cubic terms
 *
 * RETURN:
 *    the order value
 *
 * </AUTO>
 */

int atTransOrderGet() {
	/* sanity check -- make sure it's been set */
	if (at_trans_order == -1) {
		shFatal("atTransOrderGet: at_trans_order not set yet");
	}

	return (at_trans_order);
}

/************************************************************************
 *
 *
 * ROUTINE: atTransNew
 *
 * DESCRIPTION:
 * Create a new TRANS structure, and return a pointer to it.
 *
 * RETURN:
 *    TRANS *           if all goes well
 *    NULL              if not
 *
 * </AUTO>
 */

TRANS *
atTransNew() {
	TRANS *new;
	static int id_number = 0;

	new = shMalloc(sizeof(TRANS));
	new->id = id_number++;
	new->order = atTransOrderGet();
	new->nr = 0;
	new->nm = 0;
	new->sig = 0.0;
	new->sx = 0.0;
	new->sy = 0.0;
	return (new);
}

Homography *
atHNew() {
	Homography *new;
	new = shMalloc(sizeof(Homography));
	new->pair_matched = 0;
	new->Inliers = 0;
	return (new);
}

/************************************************************************
 *
 *
 * ROUTINE: atTransDel
 *
 * DESCRIPTION:
 * Delete the given TRANS structure
 *
 * RETURN:
 *    nothing
 *
 * </AUTO>
 */

void atTransDel(TRANS *trans /* I: structure to delete */
) {
	shFree(trans);
}

void atHDel(Homography *H /* I: structure to delete */
) {
	shFree(H);
}

/************************************************************************
 *
 * ROUTINE: atStarNew
 *
 * DESCRIPTION:
 * Create a new s_star structure, and return a pointer to it.  Fill
 * it with values 'x', 'y' and 'mag'.
 *
 * RETURN:
 *    s_star *          if all goes well
 *    NULL              if not
 *
 * </AUTO>
 */

struct s_star *
atStarNew(double x, /* I: x value for new star */
double y, /* I: y value for new star */
double mag, /* I: mag value for new star */
double BV /* I: BV value for new star */
) {
	struct s_star *new;
	static int id_number = 0;

	new = (struct s_star *) shMalloc(sizeof(struct s_star));
	new->id = id_number++;
	new->index = -1;
	new->x = x;
	new->y = y;
	new->mag = mag;
	new->BV = BV;
	new->match_id = -1;
	new->next = (struct s_star *) NULL;
	return (new);
}

/**********************************************************************
 * ROUTINE: is_blank
 *
 * If the given string consists only of whitespace, return 1.
 * Otherwise, return 0.
 */

int
is_blank
   (
   char *line                /* I: string to be checked */
   )
{
   char *p;

   for (p = line; (*p != '\0') && (isspace(*p)); p++) {
      ;
   }
   /* 17/Dec/01, jpb: fixed a bug that said "if (p=='\0')" */
   if (*p == '\0') {
      return(1);
   }
   else {
      return(0);
   }
}


/**********************************************************************
 * ROUTINE: get_value
 *
 * Given a string containing a numerical value, read the numerical
 * value and place it into the given double argument.
 *
 * Return 0 if all goes well.
 * Return 1 if there is an error.
 */

int
get_value
   (
   char *str,                /* I: string to be converted to double */
   double *val               /* O: place value here */
   )
{
   if (sscanf(str, "%lf", val) != 1) {
      return(1);
   }
   else {
      return(0);
   }
}

/**********************************************************************
 * ROUTINE: get_stars
 *
 *
 * Return 0 if all goes well.
 * Return 1 if there is an error.
 */

int get_stars(fitted_PSF **s, int n, int *num_stars, struct s_star **list) {
	int i = 0;
	struct s_star *head, *last, *new;

	head = (struct s_star *) NULL;
	last = head;

	while (i < n) {
		new = atStarNew(s[i]->xpos, s[i]->ypos, s[i]->mag, s[i]->BV);
		new->id = i;

		if (head == NULL) {
			head = new;
			last = new;
		} else {
			last->next = new;
			last = new;
		}
		i++;
	}

	*num_stars = i;
	*list = head;

	return (SH_SUCCESS);
}

void free_stars(struct s_star *head) {
	struct s_star* tmp;

	while (head != NULL) {
		tmp = head;
		head = head->next;
		shFree(tmp);
	}
}

void print_H(Homography *H) {
	printf("Transformation Matrix:\n");
	printf("%+*.5f %+*.5f %+*.5f\n", 11, H->h00, 11, H->h01, 11, H->h02);
	printf("%+*.5f %+*.5f %+*.5f\n", 11, H->h10, 11, H->h11, 11, H->h12);
	printf("%+*.5f %+*.5f %+*.5f\n", 11, H->h20, 11, H->h21, 11, H->h22);
}

/************************************************************************
 *
 *
 * ROUTINE: print_trans
 *
 * DESCRIPTION:
 * Print the elements of a TRANS structure.
 *
 * RETURNS:
 *   nothing
 *
 * </AUTO>
 */


void
print_trans
   (
   TRANS *trans       /* I: TRANS to print out */
   )
{
   switch (trans->order) {

   case 1:  /* linear transformation */
      printf("TRANS: a=%-15.9e b=%-15.9e c=%-15.9e d=%-15.9e e=%-15.9e f=%-15.9e",
            trans->a, trans->b, trans->c, trans->d, trans->e, trans->f);
      break;

   case 2:  /* quadratic terms */
      printf("TRANS: a=%-15.9e b=%-15.9e c=%-15.9e d=%-15.9e e=%-15.9e f=%-15.9e ",
          trans->a, trans->b, trans->c, trans->d, trans->e, trans->f);
      printf("       g=%-15.9e h=%-15.9e i=%-15.9e j=%-15.9e k=%-15.9e l=%-15.9e",
          trans->g, trans->h, trans->i, trans->j, trans->k, trans->l);
      break;

   case 3:  /* cubic terms */
      printf("TRANS: a=%-15.9e b=%-15.9e c=%-15.9e d=%-15.9e e=%-15.9e f=%-15.9e g=%-15.9e h=%-15.9e",
         trans->a, trans->b, trans->c, trans->d, trans->e, trans->f,
         trans->g, trans->h);
      printf("       i=%-15.9e j=%-15.9e k=%-15.9e l=%-15.9e m=%-15.9e n=%-15.9e o=%-15.9e p=%-15.9e",
         trans->i, trans->j, trans->k, trans->l, trans->m, trans->n,
         trans->o, trans->p);
      break;

   default:
      shFatal("print_trans: invalid trans->order %d \n", trans->order);
      exit(1);
   }

	/*
	 * we always print this information about the match at the end
	 * of the line ]
	 */
	printf(" sig=%-.4e Nr=%d Nm=%d sx=%-.4e sy=%-.4e",
	       trans->sig, trans->nr, trans->nm, trans->sx, trans->sy);
	printf(" \n");

}


/************************************************************************/
