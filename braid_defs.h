/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. Email xbraid-support@llnl.gov for support.
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

/**
 * Defines floating point type
 **/
typedef double braid_Real;

/*--------------------------------------------------------------------------
 * MPI stuff
 *--------------------------------------------------------------------------*/

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
 * Miscellaneous macros (RDF - Is this the right file location?)
 *--------------------------------------------------------------------------*/

#ifndef _braid_max
#define _braid_max(a,b)  (((a)<(b)) ? (b) : (a))
#endif
#ifndef _braid_min
#define _braid_min(a,b)  (((a)<(b)) ? (a) : (b))
#endif


#ifdef __cplusplus
}
#endif

#endif

