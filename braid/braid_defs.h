/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory.
 * 
 * This file is part of XBraid. For support, post issues to the XBraid Github page.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (as published by the Free Software
 * Foundation) version 2.1 dated February 1999.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc., 59
 * Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 ***********************************************************************EHEADER*/
 

/** \file braid_defs.h
 * \brief Definitions of braid types, error flags, etc... 
 *
 */

#ifndef braid_defs_HEADER
#define braid_defs_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------------
 * Define basic types
 *--------------------------------------------------------------------------*/

/**
 * Defines integer type
 **/
typedef int    braid_Int;
#define braid_Int_Max INT_MAX;
#define braid_Int_Min INT_MIN;

/**
 * Defines floating point type
 * Switch beween single and double precision by un-/commenting lines. 
 **/
// typedef float braid_Real;
typedef double braid_Real;

/*--------------------------------------------------------------------------
 * MPI stuff
 *--------------------------------------------------------------------------*/

 /* Switch beween single and double precision by un-/commenting lines. */
// #define braid_MPI_REAL  MPI_FLOAT
#define braid_MPI_REAL  MPI_DOUBLE
#define braid_MPI_INT   MPI_INT
#define braid_MPI_Comm  MPI_Comm


struct _braid_Vector_struct;
/**
 * This defines (roughly) a state vector at a certain time value.  
 * It could also contain any other information related to this vector which is 
 * needed to evolve the vector to the next time value, like mesh information.
 **/
typedef struct _braid_Vector_struct *braid_Vector;

#ifdef __cplusplus
}
#endif

#endif

