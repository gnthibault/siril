#if !defined(MISC_H)
#define MISC_H

#include "core/siril.h"

/*
 *
 * FILE: misc.h
 *
 * DESCRIPTION:
 * Definitions for structures and functions used in support
 * of the matching code.  Much of this is material that appears
 * in SHIVA, the SDSS environment.
 *
 */


#define SH_SUCCESS        0          /* indicates that all went well */
#define SH_GENERIC_ERROR  1          /* indicates that error occurred */


	/* a buffer used for parsing command-line arguments */
#define CMDBUFLEN       500

   /* max length of lines in input files */
#define LINELEN         300

   /* ignore any lines in input files that start with this */
#define COMMENT_CHAR   '#'

   /* data files can have this many data columns, at most */
#define MAX_DATA_COL    30

   /* each column in the data file can have at most this many characters */
#define MAX_COL_LENGTH 50


   /* possibilities for the "order" field of TRANS structures */
#define AT_TRANS_LINEAR      1           /* linear terms only */
#define AT_TRANS_QUADRATIC   2           /* linear plus quadratic */
#define AT_TRANS_CUBIC       3           /* linear plus quadratic plus cubic */

	/* maximum possible number of coefficients in a TRANS */
#define AT_TRANS_MAXCOEFF   16           /* for cubic case */


   /*
    * little wrappers around 'malloc' and 'free' functions.
    */

void *
shMalloc(int nbytes);

void
shFree(void *vptr);

void
shError(char *format, ...);

void
shFatal(char *format, ...);

void
shDebug(int level, char *format, ...);

   /*
    * This is a preprocessor macro, which acts like a function.  The user
    * calls it with one argument, a condition which should evaluate to
    * 0 or 1.  If it evaluates to 1, then nothing happens; but if it
    * evaluates to 0, then the program prints an error message, giving
    * location of the error, and halts execution with an error code.
    *
    * Thus, one typically uses it to make a 'sanity check', such as
    * making sure that a pointer (which really, really should have
    * a valid value) isn't NULL:
    *
    *     fp = fopen("file", "r");
    *        ...
    *     shAssert(fp != NULL);
    *     fgets(line, LINELEN, fp);
    *        ...
    */
#define shAssert(x) if((x)!=1){fprintf(stderr,"assertion fails in file %s, line %d\n",__FILE__,__LINE__);}


   /*
    * This is a generic transformation from one coordinate system to another.
    * Given the measured (x, y), the transformed coords (x', y') are
    * calculated like this:
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
    *
    *  The 'order' field of the TRANS structure signals which
    *  of the above cases is to be used.
	 *
	 *  The 'nr' field contains the number of pairs ultimately used to
	 *  determine the transform between the coordinate system.
	 *
	 *  The 'nm' field contains the number of pairs which -- after
	 *  the transform has been determined and applied -- match up,
	 *  and can be used to determine the quality of the fit.
	 *
	 *  The 'sig' field holds the standard deviation of the separation
	 *  between matched stars in the two sets, after they have been
	 *  transformed into the coordinate system of the second list.
	 *  Only items which were used to derive the TRANS are included
	 *  in this calculation.
	 *
	 *  The 'sx' and 'sy' fields hold the standard deviation of the
	 *  1-D separations between corresponding items in the two sets,
	 *  after they have been transformed into the coordinate system
	 *  of the second list.  A single iteration of 3-sigma clipping
	 *  is used in the calculation.
	 *  If the user specifies the "recalc" option, then the "sx" and "sy"
	 *  values will be based on all items which can be matched in
	 *  the two frames -- not just those few which were used to derive
	 *  the TRANS.
    */

typedef struct Trans {
  int id;
  int order;
  double a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p;
  int nr;
  int nm;
  double sig;
  double sx;
  double sy;
} TRANS;

void atTransOrderSet(int order);
int atTransOrderGet(void);
TRANS *atTransNew(void);
Homography *atHNew(void);
void print_trans(TRANS *trans);
void atTransDel(TRANS *trans);
void atHDel(Homography *H);
void print_H(Homography *H);

   /*
    * create a new s_star structure
    */

struct s_star *
atStarNew(double x, double y, double mag, double BV);

   /*
    * read an ASCII file with a catalog of stars, and create a list
    * of s_star structures
    */
int
read_star_file(char *filename, int xcolumn, int ycolumn, int magcolumn,
               int idcolumn, int ra_hours_col, int *num_stars,
               struct s_star **list);

   /*
    * read an ASCII file with a list of matched stars, created by
    * the atMatchLists function, and create a list of s_star structures.
    */

int
get_stars(fitted_PSF **s, int n, int *num_stars, struct s_star **list);

void
free_stars(struct s_star *head);

int
is_blank(char *line);

int
get_value(char *str, double *val);


#endif    /* MISC_H */
