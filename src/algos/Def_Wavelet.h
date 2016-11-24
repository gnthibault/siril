
/*************************************************************************
**           Copyright (C) 1993 by European Southern Observatory
**************************************************************************
**
**    UNIT
**
**    Version: 19.1
**
**    Author: Jean-Luc Starck
**
**    Date:  03/02/25
**    
**    File:  Def_Wavelet.h
**
*******************************************************************************
**
**    DESCRIPTION  Data structures definitions for the wavelet package
**    -----------   
******************************************************************************/ 

/* Wavelet transform algorithm number  */
#define TO_PAVE_LINEAR   1
#define TO_PAVE_BSPLINE  2
#define TO_PAVE_BSPLINE_FFT 3
#define TO_PYR_LINEAR 4
#define TO_PYR_BSPLINE 5
#define TO_PYR_FFT_DIFF_RESOL 6
#define TO_PYR_FFT_DIFF_SQUARE_RESOL 7
#define TO_MALLAT_BARLAUD 8

#define MAX_SIZE_NAME_IMAG 100
#define MAX_PLAN_WAVELET 20

/* Function definitions for the algorithms using FFT  */
#define SCALING_FUNCTION 1
#define FILTER_H 2
#define FILTER_H_TILDE 3
#define FILTER_G 4
#define FILTER_G_TILDE 5
#define WAVELET 6

/* Pyramid data structure */
typedef struct 	{
	int Tab_Nl[MAX_PLAN_WAVELET];
	int Tab_Col[MAX_PLAN_WAVELET];
	int Tab_Pos[MAX_PLAN_WAVELET];
	int Size;
	float Freq_Coup;
	float *Data;
	        } pyramid_f_des;

/* Complex pyramid data structure */
typedef struct 	{
	int Tab_Nl[MAX_PLAN_WAVELET];
	int Tab_Col[MAX_PLAN_WAVELET];
	int Tab_Pos[MAX_PLAN_WAVELET];
	int Size;
	float Freq_Coup; /* Frequency cutt-off */
	float *Data;
	        } pyramid_cf_des;

/* Data  structure for an algorithm without reduction of sampling */
typedef struct 	{
	float *Data;
	float Freq_Coup;  /* Frequency cutt-off */
	        } pave_f_des;


/* Data structure for Mallat's algorithm */
struct mallat_plan_des	{
	int Nl,Nc;
        float *Coef_Horiz;
        float *Coef_Diag;
        float *Coef_Vert;
        float *Low_Resol;
        struct mallat_plan_des *Smooth_Imag;
        } mallat_plan_des;

/* Data structure for a wavelet transform */
typedef struct 	{
	/* Image name */
	char Name_Imag [MAX_SIZE_NAME_IMAG];
	/* Line and column number */
	int Nbr_Ligne, Nbr_Col;
	/* Scale Number */
	int Nbr_Plan;
	/* Transform algorithm choosen */
	int Type_Wave_Transform;
	/* Buffer for the data */
	pyramid_f_des Pyramid;
	pave_f_des Pave;
	struct mallat_plan_des Mallat;
	         } wave_transf_des;

/* Data structure for image information */
typedef struct {
        float Sigma;
        float Mean;
        float Min, Max;
        float Energ, Entrop;
        float Correl_Plan[MAX_PLAN_WAVELET];
	         } plan_info_des;

/* Filtering */

    /* Thresholding */
#define FILTER_TRESHOLD 1     

    /* Adaptative thresholding */
#define FILTER_HIERARCHICAL_TRESHOLD 2

    /* Hierarchical Wiener filtering */
#define FILTER_HIERARCHICAL 3

    /* Multiresolution Wiener filtering */
#define FILTER_MULTI_RES_WIENER 4

#define TO1_FRENCH 1
#define TO1_MEX 2
#define TO1_LINEAR 3
#define TO1_B1SPLINE 4
#define TO1_B3SPLINE 5
#define TO1_MORLET 6
#define TO1_ROBUST 7
#define TO1_D1GAUS 8

int wave_io_size_data (int Nl, int Nc, int Nbr_Plan, int Type_Wave_Transform);
int wave_io_read (char *File_Name_In, wave_transf_des *Wave_Trans);
int wave_io_write (char *File_Name_In, wave_transf_des *Wave_Trans);
int wave_io_free (wave_transf_des *Wave_Trans);
int wave_io_alloc (wave_transf_des *Wave_Trans, int Type_Transform, int Nbr_Plan, int Nl, int Nc);
float *f_vector_alloc(int Nbr_Elem);
int wavelet_transform_file (float *Imag, int Nl, int Nc, char *File_Name_Transform, int Type_Transform, int Nbr_Plan, WORD *data);
int wavelet_transform_data (float *Imag, int Nl, int Nc, wave_transf_des *Wavelet, int Type_Transform, int Nbr_Plan);
int pave_2d_linear_smooth (float *Imag, float *Smooth, int Nl, int Nc, int Num_Plan);
int pave_2d_tfo (float *Pict, float *Pave, int Nl, int Nc, int Nbr_Plan, int Type_To);
int pave_2d_build (float *Pave, float *Imag, int Nl, int Nc, int Nbr_Plan, float *coef);
int pave_2d_extract_plan (float *Pave, float *Imag, int Nl, int Nc, int Num_Plan);
int pave_2d_bspline_smooth (float *Imag, float *Smooth, int Nl, int Nc, int Num_Plan);
int prepare_rawdata(float *Imag, int Nl, int Nc, WORD *data);
int wavelet_reconstruct_data (wave_transf_des *Wavelet, float *Imag, float *coef);
int wavelet_reconstruct_file (char *File_Name_Transform, float *coef, WORD *data);
int reget_rawdata(float *Imag, int Nl, int Nc, WORD *buf);
