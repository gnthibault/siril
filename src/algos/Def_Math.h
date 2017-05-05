/****************************************************************************
**           Copyright (C) 1993 by European Southern Observatory
*****************************************************************************
**
**    UNIT
**
**    Version: 19.1
**
**    Author: Jean-Luc Starck
**
**    Date:  03/02/25
**    
**    File:  Def_Math.h
**
*****************************************************************************
**
**    DECRIPTION
**    ----------
**
**   Mathematical Definitions
**
****************************************************************************/



#define OK 1
#define KO 0
#define STRING_SIZE 80

#ifndef PI
#define PI 3.1415926536 
#endif

#define	ZERO	1.0e-20
#define	WAVELET_INFINITY 1.0e+20

	/* The smallest value such that 1.0 + EPSILON != 1.0 */
#define FLOAT_EPSILON 5.96047e-08
#define DOUBLE_EPSILON 1.11077e-16
	/* DOUBLE_EPSILON software = 1.11022e-16  f68881 = 1.11077e-16 */

#define	GREATEST(x,y)  ( ((x) > (y)) ? (x) : (y) )
	/* Returns the largest of its arguments	 */

#define INT_POW(x,y,z) { int l,xx,yy; xx = (x) ; yy = (y);  for (l=0,(z)=1;l<yy;++ l,z *= xx); } 
	/*Finds the value of x raised to the power y, x and y are integer
 	expressions, z MUST be an integer VARIABLE 
	x or y are permitted to equal z  */

#define MAX_ABS(a,b) ((abs(a) > abs(b)) ? a : b)
#define MAX_FABS(a,b) ((fabs(a) > fabs(b)) ? a : b)

#define MIN_ABS(a,b) ((abs(a) > abs(b)) ? b : a)
#define MIN_FABS(a,b) ((fabs(a) > fabs(b)) ? b : a)

#define COPY_IMA(Ima1, Ima2, N) \
   { \
       int i; \
       for (i = 0; i < N; i++) \
           Ima1[i] = Ima2[i]; \
   }

#define MAX_IMA(Ima, N, Max) \
   { \
       int i; \
       Max = 0.; \
       for (i = 0; i < N; i++)  Max = MAX_FABS (Max, Ima[i]); \
   }

#define MIN_MAX_IMA(Ima, N, Min, Max) \
   { \
       int i; \
       Max = Min = Ima[0]; \
       for (i = 1; i < N; i++) \
       { \
            Max = MAX (Max, Ima[i]); \
            Min = MIN (Min, Ima[i]); \
       }\
   }

#define PRINT_MIN_MAX(Ima, N,Mes) \
   { \
       int i; \
       float Min,Max, Energ; \
       Energ = Max = Min = Ima[0]; \
       for (i = 1; i < N; i++) \
       { \
            Max = MAX (Max, Ima[i]); \
            Min = MIN (Min, Ima[i]); \
            Energ += Ima[i];\
       }\
       Energ /= (float) N;\
       printf ("%s: Min = %f, Max = %f, Ener_Moy = %f\n", Mes, Min,Max, Energ);\
   }

#define PRINT_MIN_MAX_CF(Ima, N,Mes) \
   { \
       int i; \
       float Min,Max, Energ; \
       Energ = Max = Min = Ima[0].re; \
       for (i = 1; i < N; i++) \
       { \
            Max = MAX (Max, Ima[i].re); \
            Min = MIN (Min, Ima[i].re); \
            Energ += Ima[i].re;\
       }\
       Energ /= (float) N;\
       printf ("%s:Min = %f, Max = %f, Ener_Moy = %f\n", Mes,Min, Max, Energ);\
   }



#define NORM_ENERG(Ima, N) \
   { \
       int i; \
       float Energ = 0.;\
       for (i = 0; i < N; i++)  Energ += Ima[i]; \
        for (i = 0; i < N; i++)  Ima[i] /= Energ; \
   }

#define DIV_IMA(Ima, N, Val) \
   { \
       int i; \
       for (i = 0; i < N; i++)  Ima[i] /= Val; \
   }

#define NORM_TO_1(Ima, N) \
   { \
       int i; \
       float Max=0.;\
       for (i = 0; i < N; i++)  Max = MAX_FABS (Max, Ima[i]); \
       for (i = 0; i < N; i++)  Ima[i] /= Max; \
   }

#define RAZ_IMA(Ima, N) \
   { \
       int i; \
       for (i = 0; i < N; i++) Ima[i] = 0.; \
   }

#define CF_MLT(x,y,z) \
   { \
	z.re = x.re * y.re - x.im * y.im ; \
	z.im = x.re * y.im + x.im * y.re ; \
   } 

#undef CF_DIF
#define CF_DIF(x,y,z) \
   { \
	z.re = x.re - y.re; \
	z.im = x.im - y.im ; \
   } 

#define CF_ADD(x,y,z) \
   { \
	z.re = x.re + y.re; \
	z.im = x.im + y.im ; \
   } 

#define CF_ASS(x,y) \
   { \
	y.re = x.re; \
	y.im = x.im; \
   } 


#define CF_DIV(x,y,ret) \
  { \
        double mod; \
	mod = (y.re*y.re + y.im*y.im); \
	if (mod < ZERO) ret.re = ret.im = 0.0; \
        else \
	{ \
	    ret.re = (x.re*y.re +x.im*y.im) / mod;  \
	    ret.im = (x.im*y.re -x.re*y.im) / mod; \
	} \
  } 

/* Compute the module from the real and imainary parts */
#define MOD(a,b) (sqrt((float)(a*a+b*b)))

/* Compute the module of a complex number */
#define MOD_CF(a) (sqrt((float)(a.re*a.re+a.im*a.im)))

/*  Compute the argument from the real and imainary parts */
#define ARG(a,b,Arg) \
   { \
      float Val,Va,Vb;\
      Va = (float) a; Vb = (float) b;\
      if (fabs(Va) < FLOAT_EPSILON) \
      {\
          if (fabs(Vb) < FLOAT_EPSILON)  Arg = 0.; \
          else if (Vb < 0.) Arg = PI / 2.; \
               else Arg =  - PI / 2.; \
      }\
      else \
      {\
          Val = Vb / Va; \
          Arg = atan(Val);\
      }\
   }

/* Compute the argument of a complex number */
#define ARG_CF(a,Arg) \
   { \
      double Val;\
      if (fabs((float) a.re) < FLOAT_EPSILON) \
                if (fabs((float) a.im) < FLOAT_EPSILON)  Arg = 0.; \
                else if (a.im < 0.) Arg = PI / 2.; \
                     {\
                        Val = (float) a.im / (float) a.re; \
                        Arg = atan(Val);\
                     }\
   }

typedef char string[80];

enum    direction       {FORWARD, REVERSE};


extern int ft_cf_any_power_of_2();
/*
extern int ft_cf_any_power_of_2(dat,direction,length)
        complex_float   *dat;
        int             direction;
        int             length;
 
 Takes the array pointed to by dat which is assumed to be of side-length
length and performs a fourier transform of its contents

External function calls
        None
 
Return codes
0 :     No problems encountered.
1 :     Length was not a power of 2 FATAL_ERROR
 
*/


double	gaussian_num();
/*
double	gaussian_num(mean,std,finesse)
		double	mean;
		double	std;
		int	finesse;

Generates a gaussian number of the requested mean and standard deviation by
repeated summation of linear random variables.
The parameter finesse indicates the number of linear random numbers to be summed,
a value of 9 will give excellent gaussian statistics. 
*/

