/******************************************************************************
**           Copyright (C) 1993 by European Southern Observatory
*******************************************************************************
**
**    UNIT
**
**    Version: 19.1
**
**    Author: Jean-Luc Starck
**
**    Date:  03/02/25
**    
**    File:  Def_Mem.h
**
*******************************************************************************
**
**    DECRIPTION  Memory and errors definition
**    ---------- 
**
******************************************************************************/
/*		Commented out by dulle who will alloc himself
extern  float   *f_vector_alloc();            
extern  double  *d_vector_alloc(); 
extern  int  *i_vector_alloc();           
extern  int     **i_matrix_alloc();           
extern  void    i_matrix_free();               
extern  float   **f_matrix_alloc();            
extern  void    f_matrix_free();              
extern  complex_float  *cf_vector_alloc();           
extern  complex_float  **cf_matrix_alloc();        
extern  void    cf_matrix_free();     
*/

/* extern  char    *malloc(); */


/* error definitions */ 

#define NUMBER_ERROR 13

#define ERR_TRANSF     1
#define ERR_POWER_OF_2 2
#define ERR_READ_DATA  3
#define ERR_WRITE_DATA 4
#define ERR_OPEN_FILE  5
#define ERR_CLOSE_FILE 6
#define ERR_ALLOC_MEMO 7
#define ERR_PLANE_NUMBER 8
#define ERR_IMAGE_SQUARE 9
#define ERR_IMAGE_SIZE   10
#define ERR_NOT_IMPLEMENTED  11
#define ERR_NUMBER_OF_PLANES  12
#define ERR_INTERP_PLANE_NUMBER 13

/*  io_err_message_exit (Num_Err, Mes);*/

/* error definitions */
#define MAX_TAB NUMBER_ERROR+1

