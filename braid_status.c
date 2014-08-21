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

/** \file braid_status.c
 * \brief Source code for status interface routines.  See braid_status.h 
 * for more information.
 *
 */

#include "braid_status.h"

#ifndef DEBUG
#define DEBUG 0
#endif

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_AccessStatusGetResidual(braid_AccessStatus  status,
                              braid_Real         *rnorm_ptr)
{
   *rnorm_ptr = _braid_StatusElt(status, rnorm);
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_AccessStatusGetIter(braid_AccessStatus  status,
                          braid_Int          *iter_ptr)
{
   *iter_ptr = _braid_StatusElt(status, iter);
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_AccessStatusGetLevel(braid_AccessStatus  status,
                           braid_Int          *level_ptr)
{
   *level_ptr = _braid_StatusElt(status, level);
   return _braid_error_flag;
}


/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_AccessStatusGetDone(braid_AccessStatus  status,
                          braid_Int          *done_ptr)
{
   *done_ptr = _braid_StatusElt(status, done);
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Destroy a braid_AccessStatus structure
 *--------------------------------------------------------------------------*/
braid_Int
_braid_AccessStatusDestroy(braid_AccessStatus  status)
{
   if (status)
   {
      _braid_TFree(status);
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Initialize members of a braid_AccessStatus structure 
 *--------------------------------------------------------------------------*/
braid_Int
_braid_AccessStatusInit(braid_Real           rnorm,
                        braid_Int            iter,
                        braid_Int            level,
                        braid_Int            done,
                        braid_AccessStatus   status)
{
   _braid_StatusElt(status, level) = level;
   _braid_StatusElt(status, rnorm) = rnorm;
   _braid_StatusElt(status, done)  = done;
   _braid_StatusElt(status, iter)  = iter; 

   return _braid_error_flag;
}

