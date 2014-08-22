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
 * AccessStatus Routines
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


/*--------------------------------------------------------------------------
 * CoarsenRefStatus Routines
 *--------------------------------------------------------------------------*/

braid_Int
_braid_CoarsenRefStatusInit(braid_Real              tstart,  
                            braid_Real              f_tminus,
                            braid_Real              f_tplus, 
                            braid_Real              c_tminus,
                            braid_Real              c_tplus,
                            braid_CoarsenRefStatus  status)
{
   _braid_StatusElt(status, tstart) = tstart;
   _braid_StatusElt(status, f_tminus) = f_tminus;
   _braid_StatusElt(status, f_tplus)  = f_tplus;
   _braid_StatusElt(status, c_tminus) = c_tminus;
   _braid_StatusElt(status, c_tplus)  = c_tplus;

   return _braid_error_flag;
}

braid_Int
_braid_CoarsenRefStatusDestroy(braid_CoarsenRefStatus  status)
{
   if (status)
   {
      _braid_TFree(status);
   }

   return _braid_error_flag;
}

braid_Int
braid_CoarsenRefStatusGetTstart(braid_CoarsenRefStatus  status,
                                braid_Real             *tstart_ptr
                                )
{
   *tstart_ptr = _braid_StatusElt(status, tstart);
   return _braid_error_flag;
}

braid_Int
braid_CoarsenRefStatusGetFTplus(braid_CoarsenRefStatus  status,
                                braid_Real             *f_tplus_ptr
                                )
{
   *f_tplus_ptr = _braid_StatusElt(status, f_tplus);
   return _braid_error_flag;
}

braid_Int
braid_CoarsenRefStatusGetFTminus(braid_CoarsenRefStatus  status,
                                 braid_Real             *f_tminus_ptr
                                 )
{
   *f_tminus_ptr = _braid_StatusElt(status, f_tminus);
   return _braid_error_flag;
}

braid_Int
braid_CoarsenRefStatusGetCTplus(braid_CoarsenRefStatus  status,
                                braid_Real             *c_tplus_ptr
                                )
{
   *c_tplus_ptr = _braid_StatusElt(status, c_tplus);
   return _braid_error_flag;
}

braid_Int
braid_CoarsenRefStatusGetCTminus(braid_CoarsenRefStatus  status,
                                 braid_Real             *c_tminus_ptr
                                 )
{
   *c_tminus_ptr = _braid_StatusElt(status, c_tminus);
   return _braid_error_flag;
}


/*--------------------------------------------------------------------------
 * CoarsenRefStatus Routines
 *--------------------------------------------------------------------------*/
braid_Int
_braid_PhiStatusInit(braid_Real       tstart,
                     braid_Real       tplus,
                     braid_Real       accuracy,
                     braid_PhiStatus  status)
{
   _braid_StatusElt(status, tstart)   = tstart;
   _braid_StatusElt(status, tplus)    = tplus;
   _braid_StatusElt(status, accuracy) = accuracy;

   return _braid_error_flag;
}

braid_Int
_braid_PhiStatusDestroy(braid_PhiStatus  status)
{
   if (status)
   {
      _braid_TFree(status);
   }

   return _braid_error_flag;
}


braid_Int
braid_PhiStatusGetTstart(braid_PhiStatus  status,
                         braid_Real      *tstart_ptr)
{
   *tstart_ptr = _braid_StatusElt(status, tstart);
   return _braid_error_flag;
}

braid_Int
braid_PhiStatusGetTplus(braid_PhiStatus  status,
                        braid_Real      *tplus_ptr)
{
   *tplus_ptr = _braid_StatusElt(status, tplus);
   return _braid_error_flag;
}

braid_Int
braid_PhiStatusGetAccuracy(braid_PhiStatus  status,
                           braid_Real      *accuracy_ptr)
{
   *accuracy_ptr = _braid_StatusElt(status, accuracy);
   return _braid_error_flag;
}

braid_Int
braid_PhiStatusSetRFactor(braid_PhiStatus  status,
                          braid_Real       rfactor)
{
   _braid_StatusElt(status, rfactor) = rfactor;
   return _braid_error_flag;
}

