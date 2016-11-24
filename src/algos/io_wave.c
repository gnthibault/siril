/* @(#)io_wave.c	19.1 (ES0-DMD) 02/25/03 13:34:39 */
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
**                   Copyright (C) 1993 by European Southern Observatory
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
**    File:  io_wave.c
**
*******************************************************************************
**
**    DESCRIPTION  inpout output routines for the wavelet files
**    -----------  
**
*******************************************************************************
**
** wave_io_read (File_Name, Wave_Trans)
** char *File_Name;
** wave_transf_des *Wave_Trans;
**
** read a wavelet transform file and stores the results in a structure
**
** File_Name = INPUT: file name
** Wave_Trans = OUTPUT: wavelet
**
*******************************************************************************
**
** wave_io_write (File_Name, Wave_Trans)
** char *File_Name;
** wave_transf_des *Wave_Trans;
**
** writes a wavelet transform in a file
**
** File_Name = INPUT: file name
** Wave_Trans = INPUT: wavelet
**
*******************************************************************************
**
** wave_io_free (Wave_Trans)
** wave_transf_des *Wave_Trans;
**
** deallocates the memory of a wavelet
**
*******************************************************************************
**
** int wave_io_size_data (Nl, Nc, Nbr_Plan, Type_Wave_Transform)
** int Nl, Nc, Type_Wave_Transform
**
** computes the necessary memory size for a wavelet transform algorithm
** defined by Type_Wave_Transform and for a number of scales Nbr_Plan
**
** Nl, Nc = original image size
** Nbr_Plan = number of scales = number of planes
** Type_Wave_Transform = choosen wavelet transform algorithm
** 
*******************************************************************************
**
** wave_io_alloc (Wave_Trans, Type_Transform, Nbr_Plan, Nl, Nc)
** wave_transf_des *Wave_Trans;
** int Type_Transform;
** int Nbr_Plan, Nl, Nc;
** 
** allocates the memory for a wavelet transform algorithm
** defined by Type_Wave_Transform and for a number of scales Nbr_Plan
** 
** Nl, Nc = original image size
** Nbr_Plan = number of scales = number of planes
** Type_Wave_Transform = choosen wavelet transform algorithm
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

/****************************************************************************/

int wave_io_size_data (Nl, Nc, Nbr_Plan, Type_Wave_Transform)
int Nl, Nc, Nbr_Plan, Type_Wave_Transform;
{
    int Size;

    switch (Type_Wave_Transform)
    {
        case TO_PAVE_LINEAR:
        case TO_PAVE_BSPLINE:
                Size = Nbr_Plan * Nl * Nc;
                break;
        default : 
		printf("wave_io_read: wrong transform type\n");
                return(1);
                break;
    }
    return (Size);

}
/***************************************************************************/

int static wave_io_name (File_Name_In, File_Name_Out)
char *File_Name_In, *File_Name_Out;
{
    int L;

    strcpy (File_Name_Out, File_Name_In);

    L = strlen (File_Name_In);
    if ((L < 5) || (File_Name_In[L-1] != 'e')
                || (File_Name_In[L-2] != 'v')
                || (File_Name_In[L-3] != 'a')
                || (File_Name_In[L-4] != 'w')
                || (File_Name_In[L-5] != '.'))
    {
        strcat (File_Name_Out, ".wave");
    }
return 0;
}

/****************************************************************************/

int wave_io_read (File_Name_In, Wave_Trans)
	char *File_Name_In;
	wave_transf_des *Wave_Trans;
{
	FILE *File_Des;
	int Nbr,Nl,Nc,Nbr_Plan,Size;
	char *Ptr;
	string File_Name;

	/* find the exact file name which must finish by .wave */
	wave_io_name (File_Name_In, File_Name);

	/* open the file */
	File_Des = fopen (File_Name, "r"); 
	if (File_Des == NULL){
		printf("wave_io_read: error opening file: %s\n", File_Name_In);
		return 1;
	}

	/* read the descriptor */
	Nbr = fread ((char *) Wave_Trans, sizeof(wave_transf_des), 1, File_Des);

	if (Nbr <= 0) {
		printf("wave_io_read: error reading data\n");
		fclose(File_Des);
		return 1;
	}

	Nl = Wave_Trans -> Nbr_Ligne;
	Nc = Wave_Trans -> Nbr_Col;
	Nbr_Plan = Wave_Trans -> Nbr_Plan;

	switch (Wave_Trans -> Type_Wave_Transform)
	{
		case TO_PAVE_LINEAR:
		case TO_PAVE_BSPLINE:
			Size = Nbr_Plan * Nl * Nc;
			(Wave_Trans -> Pave).Data = f_vector_alloc(Size);
			Ptr = (char *) ((Wave_Trans -> Pave).Data);
			Nbr = fread (Ptr, sizeof(float), Size, File_Des);
			if (Nbr <= 0){
				printf("wave_io_read: error reading data\n");
				fclose(File_Des);
				return 1;
			}
			break;
		default : 
			printf("wave_io_read: wrong transform type\n");
			fclose(File_Des);
			return 1;
	}
	if (fclose (File_Des) != 0) {
		printf("wave_io_read: error closing file\n");
		return 1;
	}
	return 0;
}

/****************************************************************************/

int wave_io_write (File_Name_In, Wave_Trans)
	char *File_Name_In;
	wave_transf_des *Wave_Trans;
{
	FILE *File_Des;
	int Nbr,Nl,Nc,Nbr_Plan,Size;
	char *Ptr;
	string File_Name;

	/* find the exact file name */
	wave_io_name (File_Name_In, File_Name);

	/* open the file */
	File_Des = fopen (File_Name, "w");
	if (File_Des == NULL) {
		printf("wave_io_write: error opening file\n");
		return 1;
	}

	/* write the descriptor */
	Nbr = fwrite ((char *) Wave_Trans, sizeof(wave_transf_des), 1, File_Des);
	if (Nbr <= 0) {
		printf("wave_io_write: error writing data\n");
		fclose(File_Des);
		return 1;
	}

	Nl = Wave_Trans -> Nbr_Ligne;
	Nc = Wave_Trans -> Nbr_Col;
	Nbr_Plan = Wave_Trans -> Nbr_Plan;

	switch (Wave_Trans -> Type_Wave_Transform)
	{
		case TO_PAVE_LINEAR:
		case TO_PAVE_BSPLINE:
			Ptr = (char *) ((Wave_Trans -> Pave).Data);
			Size = Nbr_Plan * Nl * Nc;
			Nbr = fwrite (Ptr, sizeof(float), Size, File_Des);
			if (Nbr <= 0) {
				printf("wave_io_write: error writing data\n");
				fclose(File_Des);
				return 1;
			}
			break;
		default : 
			printf("wave_io_write: wrong transform type\n");
			fclose(File_Des);
			return 1;
	}
	if (fclose (File_Des) != 0){
		printf("wave_io_write: error closing file\n");
		return 1;
	}
	return 0;
}

/****************************************************************************/

int wave_io_free (Wave_Trans)
wave_transf_des *Wave_Trans;
{
	switch (Wave_Trans -> Type_Wave_Transform)
    {
		case TO_PAVE_LINEAR:
		case TO_PAVE_BSPLINE:
        	if ((Wave_Trans -> Pave).Data) free ((char *) ((Wave_Trans -> Pave).Data));
        	break;
        default : 
		printf("wave_io_free: wrong transform type\n");
			return 1;
			break;
    }
	return 0;
}



/****************************************************************************/

int wave_io_alloc (Wave_Trans, Type_Transform, Nbr_Plan, Nl, Nc)
wave_transf_des *Wave_Trans;
int Type_Transform;
int Nbr_Plan, Nl, Nc;
{
    int Size;
  
    Wave_Trans -> Nbr_Ligne = Nl;
    Wave_Trans -> Nbr_Col = Nc;
    Wave_Trans -> Nbr_Plan = Nbr_Plan;
    Wave_Trans -> Type_Wave_Transform = Type_Transform;
    Wave_Trans->Pyramid.Freq_Coup = (float)1;

    switch (Wave_Trans -> Type_Wave_Transform)
    {
		case TO_PAVE_LINEAR:
		case TO_PAVE_BSPLINE:
        	Size = Nbr_Plan * Nl * Nc;
        	(Wave_Trans -> Pave).Data = f_vector_alloc(Size);
        	break;
		default :
			printf("wave_io_alloc: wrong transform type\n");
			return 1;
			break;
	}
    return 0;
}

/****************************************************************************/
