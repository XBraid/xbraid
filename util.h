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

/** \file util.h
 * \brief Define headers for utility routines.
 *
 * This file contains the headers for utility routines. Essentially,
 * if a routine does not take braid_Core (or other Braid specific structs) 
 * as an argument, then it's a utility routine.
 */

#ifndef braid_util_HEADER
#define braid_util_HEADER

#include "_braid.h"

/**
 * Project an interval onto a strided index space that contains the index
 * 'index' and has stride 'stride'.  An empty projection is represented by
 * ilower > iupper.
 **/
braid_Int
_braid_ProjectInterval( braid_Int   ilower,
                        braid_Int   iupper,
                        braid_Int   index,
                        braid_Int   stride,
                        braid_Int  *pilower,
                        braid_Int  *piupper );

/**
 * Determine the accuracy used for the spatial solves based on the ratio of
 * the current residual norm and the stopping tolerance. 
 **/
braid_Int
_braid_SetAccuracy( braid_Real   rnorm,
                    braid_Real   loose_tol,
                    braid_Real   tight_tol,
                    braid_Real   oldAccuracy,
                    braid_Real   tol,
                    braid_Real  *paccuracy );

/**
 * If set, print to @ref _braid_printfile and then flush.  
 * Else print to standard out.\n
 *
 * The string *format* can be multiple parameters
 * in the standard * C-format, like\n
 * `` format = '%1.2e is a format string', 1.24 ``
 **/
braid_Int
_braid_printf( const char *format, ...);

/**
 * This is a function that allows for "sane" printing
 * of information in parallel.  Currently, only 
 * myid = 0 prints, but this can be updated as needs change.
 *
 * The string *message* is printed and can be multiple parameters
 * in the standard * C-format, like\n
 * `` message = '%1.2e is a format string', 1.24 ``
 **/
braid_Int
_braid_ParFprintfFlush(FILE * file, 
                       braid_Int myid,
                       char * message, 
                       ...);

#endif

