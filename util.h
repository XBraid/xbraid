/*BHEADER**********************************************************************
 * Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of WARP.  See file COPYRIGHT for details.
 *
 * WARP is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 ***********************************************************************EHEADER*/

/** \file util.h
 * \brief Define headers for utility routines.
 *
 * This file contains the headers for utility routines. Essentially,
 * if a routine does not take warp_Core (or other warp specific structs) 
 * as an argument, then it's a utility routine.
 */

#ifndef warp_util_HEADER
#define warp_util_HEADER

#include "_warp.h"

/**
 * Project an interval onto a strided index space that contains the index
 * 'index' and has stride 'stride'.  An empty projection is represented by
 * ilower > iupper.
 **/
warp_Int
_warp_ProjectInterval( warp_Int   ilower,
                       warp_Int   iupper,
                       warp_Int   index,
                       warp_Int   stride,
                       warp_Int  *pilower,
                       warp_Int  *piupper );

/**
 * Determine the accuracy used for the spatial solves based on the ratio of
 * the current residual norm and the stopping tolerance. 
 **/
warp_Int
_warp_SetAccuracy( warp_Real   rnorm,
                   warp_Real   loose_tol,
                   warp_Real   tight_tol,
                   warp_Real   oldAccuracy,
                   warp_Real   tol,
                   warp_Real  *paccuracy );

/**
 * If set, print to @ref _warp_printfile and then flush.  
 * Else print to standard out.\n
 *
 * The string *format* can be multiple parameters
 * in the standard * C-format, like\n
 * `` format = '%1.2e is a format string', 1.24 ``
 **/
warp_Int
_warp_printf( const char *format, ...);

/**
 * This is a function that allows for "sane" printing
 * of information in parallel.  Currently, only 
 * myid = 0 prints, but this can be updated as needs change.
 *
 * The string *message* is printed and can be multiple parameters
 * in the standard * C-format, like\n
 * `` message = '%1.2e is a format string', 1.24 ``
 **/
warp_Int
_warp_ParFprintfFlush(FILE * file, 
                      warp_Int myid,
                      char * message, 
                      ...);

#endif

