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

warp_Int
_warp_ProjectInterval( warp_Int   ilower,
                       warp_Int   iupper,
                       warp_Int   index,
                       warp_Int   stride,
                       warp_Int  *pilower,
                       warp_Int  *piupper );

warp_Int
_warp_SetAccuracy( warp_Float  rnorm,
                   warp_Float  loose_tol,
                   warp_Float  tight_tol,
                   warp_Float  oldAccuracy,
                   warp_Float  tol,
                   warp_Float *paccuracy );


#endif
