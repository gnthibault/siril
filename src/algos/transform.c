/* @(#)transform.c	19.1 (ES0-DMD) 02/25/03 13:34:40 */
/*===========================================================================
  Copyright (C) 1995 European Southern Observatory (ESO)
 
  This program is free software; you can redistribute it and/or 
  modify it under the terms of the GNU General Public License as 
  published by the Free Software Foundation; either version 2 of 
  the License, or (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public 
  License along with this program;
  If not, see <http://www.gnu.org/licenses/>.
 
  Corresponding concerning ESO-MIDAS should be addressed as follows:
	Internet e-mail: midas@eso.org
	Postal address: European Southern Observatory
			Data Management Division 
			Karl-Schwarzschild-Strasse 2
			D 85748 Garching bei Muenchen 
			GERMANY
===========================================================================*/

/******************************************************************************
**              Copyright (C) 1993 by European Southern Observatory
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
**    File:  transform.c
**
*******************************************************************************
**
**    DESCRIPTION  This module contains wavelet transform routines
**    -----------                 
******************************************************************************
**
** wavelet_transform_file (File_Name_Imag, File_Name_Transform, 
**                         Type_Transform, Fc, Nbr_Plan)
** char *File_Name_Imag, *File_Name_Transform;
** int Type_Transform;
** float Fc;
** int Nbr_Plan;
** 
** Computes the wavelet transform of the image with a file name File_Name_Imag
** and write the result in the file File_Name_Transform
** 
** File_Name_Imag = File name of the imput image
** File_Name_Transform = File name of the output wavelet transform
** Type_Transform = wavelet transform algorithm number 
** Fc = cut-off frequency if the algorithm use the FFT
** Nbr_Plan = number of scales
**
******************************************************************************
**
** wavelet_transform_data (Imag, Nl, Nc, Wavelet, Type_Transform, Fc, Nbr_Plan)
** float *Imag;
** int Nl, Nc;
** wave_transf_des *Wavelet;
** int Type_Transform;
** float Fc;
** int Nbr_Plan;
** 
** Computes the wavelet transform of the image Imag
** and write the result in the structure Wavelet
**
** Imag = INPUT:image
** Wavelet = OUTPUT:wavelet 
** Nl,Nc = INPUT: number of lines and columns 
** Type_Transform = INPUT:wavelet transform algorithm number 
**           Which_Algorithm = 1...7
**           1 ==> a trous algorithm with a linear scaling function
**                 the wavelet function is the difference betwwen two resolutions
**           2 ==> a trous algorithm with a B3-spline scaling function
**                 the wavelet function is the difference betwwen two resolutions
**           3 ==> algorithm using the Fourier transform:
**                 without any reduction of the samples between two scales
**                 the Fourier transform of the scaling function is a b3-spline
**                 the wavelet function is the difference between two resolutions
**           4 ==> pyramidal algorithm in the direct space, with a linear 
**                 scaling function
**           5 ==> pyramidal algorithm in the direct space, with a b3-spline 
**                 scaling function
**           6 ==> algorithm using the Fourier transform:
**                 with a reduction of the samples between two scales
**                 the Fourier transform of the scaling function is a b3-spline
**                 the wavelet function is the difference between two resolutions
**           7 ==> algorithm using the Fourier transform:
**                 with a reduction of the samples between two scales
**                 the Fourier transform of the scaling function is a b3-spline
**                 the wavelet function is the difference between the square of
**                 two resolutions
**           8 ==> Mallat's Algorithm with biorthogonal filters.
** Fc = INPUT:cut-off frequency if the algorithm use the FFT
** Nbr_Plan = INPUT:number of scales
**
******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "core/siril.h"
#include "gui/callbacks.h"
#include "algos/Def_Math.h"
#include "algos/Def_Mem.h"
#include "algos/Def_Wavelet.h"

int prepare_rawdata(float *Imag, int Nl, int Nc, WORD *buf){
	float *im=Imag;
	int i;

	for (i=0;i<(Nl)*(Nc);++i){
		(im[i])=(float)buf[i];
	}
	return 0;
}


/****************************************************************************/

float *f_vector_alloc(Nbr_Elem)
	/* allocates a vector of float */
	int  Nbr_Elem;
{
	float *Vector;

	Vector = (float*) calloc ((unsigned)Nbr_Elem * sizeof (float),1);
	if (Vector == NULL) {
		printf("wavelet: memory error\n");
	}	
	return(Vector);
}

/*****************************************************************************/

int wavelet_transform_file (Imag, Nl, Nc, File_Name_Transform, Type_Transform, Nbr_Plan, data)
float *Imag;
int Nl;
int Nc;
char *File_Name_Transform;
int Type_Transform;
int Nbr_Plan;
WORD *data;
{
	wave_transf_des Wavelet;
    
    memset(&Wavelet, 0, sizeof(wave_transf_des));

	/* read the input image */
	prepare_rawdata (Imag, Nl, Nc, data);
	snprintf (Wavelet.Name_Imag, MAX_SIZE_NAME_IMAG-1, "%s", File_Name_Transform);
	Wavelet.Name_Imag[MAX_SIZE_NAME_IMAG-1] = '\0';

	if (wavelet_transform_data (Imag, Nl, Nc, &Wavelet, Type_Transform, Nbr_Plan)) {
		return 1;
	}
	if (wave_io_write (File_Name_Transform, &Wavelet)) {
		return 1;
	}

	wave_io_free (&Wavelet);
	return 0;
}

/*****************************************************************************/

int wavelet_transform_data (Imag, Nl, Nc, Wavelet, Type_Transform, Nbr_Plan)
float *Imag;
int Nl, Nc;
wave_transf_des *Wavelet;
int Type_Transform;
int Nbr_Plan;
{
	float *Pave;
	double Exp;
	int Size,Min,temp;
	
	Wavelet->Nbr_Ligne = Nl;
	Wavelet->Nbr_Col = Nc;
	Wavelet->Nbr_Plan = Nbr_Plan;
	Wavelet->Type_Wave_Transform = Type_Transform;
	
	/* test if the number of planes is not too high */
	Min = (Nl < Nc) ? Nl : Nc;
	Exp = (double) Nbr_Plan + 2.;
	temp = pow(2., Exp) + 0.5;
	if (Min < temp) {
		siril_log_message(_("wavelet_transform_data: bad plane number\n"));
		return 1;
	}
	
	switch (Type_Transform)	{
		case TO_PAVE_LINEAR:
		case TO_PAVE_BSPLINE:
			Size = Nl * Nc * Nbr_Plan;
			Wavelet->Pave.Data = f_vector_alloc (Size);
			if (Wavelet->Pave.Data == NULL){
				return 1;
			}
			Pave = Wavelet->Pave.Data;
			pave_2d_tfo (Imag, Pave, Nl, Nc, Nbr_Plan, Type_Transform);
			break;
		default:
			printf("wavelet_transform_data: wrong transform type\n");
			return 1;
			break;
	} 
    return 0;
}

/*****************************************************************************/

