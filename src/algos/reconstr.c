/* @(#)reconstr.c	19.1 (ES0-DMD) 02/25/03 13:34:40 */
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
**    File:  reconstr.c
**
*******************************************************************************
**
**    DESCRIPTION  reconstruction routines
**    -----------
**
******************************************************************************
**
** wavelet_reconstruct_file (File_Name_Imag, 
**                           File_Name_Transform, Build_Direct_Ok)
** char *File_Name_Imag, *File_Name_Transform;
** Build_Direct_Ok;
**
** Reconstructs an image from its the wavelet transform of file name 
** File_Name_Transform and write the result in the file File_Name_Imag
**
** File_Name_Imag = File name of the output image
** File_Name_Transform = File name of the input wavelet transform
** Build_Direct_Ok = input paramater (TRUE=1 or FALSE=0)
**    if the wavelet transform algorithm used the FFT and is pyramidal
**        (Type_Wave_Transform =  TO_PYR_FFT_DIFF_RESOL
**                             or TO_PYR_FFT_DIFF_SQUARE_RESOL)
**        then if Build_Direct_Ok = 1 (TRUE) then
**                       the reconstruction is done by addition
**                       of the wavelet coefficient in the Fourier space
**             else (Build_Direct_Ok = 0 (FALSE))  
**                       the reconstruction is done from a least mean square
**                       estimation
**
******************************************************************************
**
** wavelet_reconstruct_data (Wavelet, Imag, Build_Direct_Ok)
** float *Imag;
** wave_transf_des *Wavelet;
** int Build_Direct_Ok;
**
** Reconstructs an image from its the wavelet transform 
**
** Imag = OUTPUT:image
** Wavelet = INPUT:wavelet
** Build_Direct_Ok = input paramater (TRUE=1 or FALSE=0)
**    if Wavelet->Type_Wave_Transform =  TO_PYR_FFT_DIFF_RESOL
**                                    or TO_PYR_FFT_DIFF_SQUARE_RESOL
**        then if Build_Direct_Ok = 1 (TRUE) then
**                       the reconstruction is done by addition
**                       of the wavelet coefficient in the Fourier space
**             else (Build_Direct_Ok = 0 (FALSE))  
**                       the reconstruction is done from a least mean square
**                       estimation

******************************************************************************
**
** int W_Pyr_Rec_Iter_Number = 1; 
**
** Global variable used to define the number of iterations in
** Van Ciitert's iterative reconstruction. This parameter is used if
**  Wavelet->Type_Wave_Transform = TO_PYR_BSPLINE
**                              or TO_PYR_LINEAR
** 
** This parameter can be modified by an external program by:
**
** extern int W_Pyr_Rec_Iter_Number;
** W_Pyr_Rec_Iter_Number = Value;
**
******************************************************************************/


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gtk/gtk.h>

#include "gui/callbacks.h"
#include "core/siril.h"
#include "core/proto.h"
#include "algos/Def_Math.h"
#include "algos/Def_Mem.h"
#include "algos/Def_Wavelet.h"

int reget_rawdata(float *Imag, int Nl, int Nc, WORD *buf){
	float *im = Imag;
	float maximum = 0.f;
	double ratio;
	int i;

	gfit.ry = Nl;
	gfit.rx = Nc;

	for (i=0;i<Nl*Nc;++i){
		maximum = max(maximum, im[i]);
	}
	if (maximum > USHRT_MAX)
		ratio = USHRT_MAX_DOUBLE / (double)maximum;
	else	ratio = 1.0;
	
	for (i=0;i<Nl*Nc;++i){
		buf[i]=round_to_WORD(im[i] * ratio);
	}
	return 0;
}

/*****************************************************************************/

int wavelet_reconstruct_file (File_Name_Transform, coef, data)
	char *File_Name_Transform;
	float *coef;
	WORD *data;
{
	float *Imag;
	wave_transf_des Wavelet;
	int Nl, Nc;

	/* read the wavelet file */
	if (wave_io_read (File_Name_Transform, &Wavelet))
		return 1;

	Nl =  Wavelet.Nbr_Ligne;
	Nc = Wavelet.Nbr_Col;
	Imag = f_vector_alloc (Nl*Nc);
	if (Imag == NULL) return 1;
	wavelet_reconstruct_data (&Wavelet, Imag, coef);

	/* get and view result */
	reget_rawdata (Imag, Nl, Nc, data);

	wave_io_free (&Wavelet);
	free ((char *) Imag);
	return 0;
}

/*****************************************************************************/

int wavelet_reconstruct_data (Wavelet, Imag, coef)
	float *Imag, *coef;
	wave_transf_des *Wavelet;
{
	float *Pave;
	int Nl, Nc, Nbr_Plan;

	Nl = Wavelet->Nbr_Ligne;
	Nc = Wavelet->Nbr_Col;
	Nbr_Plan = Wavelet->Nbr_Plan;
	switch (Wavelet->Type_Wave_Transform)
	{
		case TO_PAVE_LINEAR:
		case TO_PAVE_BSPLINE:
			Pave = Wavelet->Pave.Data;
			pave_2d_build (Pave, Imag, Nl, Nc, Nbr_Plan, coef);
			break;
		default:
			siril_log_message ("Unknown transform\n");
			return 1;
			break;
	}   
	return 0;
}

/*****************************************************************************/
