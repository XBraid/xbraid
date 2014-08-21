/*BHEADER**********************************************************************
 * Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of XBraid.  See file COPYRIGHT for details.
 *
 * XBraid is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 ***********************************************************************EHEADER*/

/** \file braid_defs.h
 * \brief Definitions of types, error flags, etc... 
 *
 */

#ifndef braiddefs_HEADER
#define braiddefs_HEADER

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
 * Error handling
 *--------------------------------------------------------------------------*/

/** 
 * This is the global Braid error flag.  If it is ever nonzero, an error has 
 * occurred. 
 **/
extern braid_Int _braid_error_flag;


#ifdef __cplusplus
}
#endif

#endif

