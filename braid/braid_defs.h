/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
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
 * \brief Definitions of types, error flags, etc... 
 *
 */

#ifndef braiddefs_HEADER
#define braiddefs_HEADER

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

/*--------------------------------------------------------------------------
 * Error handling
 *--------------------------------------------------------------------------*/

/** 
 * This is the global XBraid error flag.  If it is ever nonzero, an error has 
 * occurred. 
 **/
extern braid_Int _braid_error_flag;

void _braid_ErrorHandler(const char *filename, braid_Int line, braid_Int ierr, const char *msg);
#define _braid_Error(IERR, msg)       _braid_ErrorHandler(__FILE__, __LINE__, IERR, msg)
#define _braid_ErrorInArg(IARG, msg)  _braid_Error(HYPRE_ERROR_ARG | IARG<<3, msg)

/*--------------------------------------------------------------------------
 * Memory allocation macros
 *--------------------------------------------------------------------------*/

/** 
 * Allocation macro 
 **/
#define _braid_TAlloc(type, count) \
( (type *)malloc((size_t)(sizeof(type) * (count))) )

/** 
 * Allocation macro 
 **/
#define _braid_CTAlloc(type, count) \
( (type *)calloc((size_t)(count), (size_t)sizeof(type)) )

/** 
 * Re-allocation macro 
 **/
#define _braid_TReAlloc(ptr, type, count) \
( (type *)realloc((char *)ptr, (size_t)(sizeof(type) * (count))) )

/** 
 * Free memory macro 
 **/
#define _braid_TFree(ptr) \
( free((char *)ptr), ptr = NULL )

/*--------------------------------------------------------------------------
 * Miscellaneous macros and functions 
 *--------------------------------------------------------------------------*/

#ifndef _braid_max
#define _braid_max(a,b)  (((a)<(b)) ? (b) : (a))
#endif
#ifndef _braid_min
#define _braid_min(a,b)  (((a)<(b)) ? (a) : (b))
#endif
#ifndef braid_isnan
#define braid_isnan(a) (a != a)
#endif

/**
 * Machine independent pseudo-random number generator is defined in Braid.c
 */
#ifndef braid_RAND_MAX
#define braid_RAND_MAX 32768
#endif

#ifdef __cplusplus
}
#endif

#endif

