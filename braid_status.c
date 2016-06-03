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
 
/** \file braid_status.c
 * \brief Source code for status interface routines.  See braid_status.h 
 * for more information.
 *
 */

#include "braid_status.h"
#include "_braid.h"
#include "braid_defs.h"
#include "_util.h"

#ifndef DEBUG
#define DEBUG 0
#endif

/*--------------------------------------------------------------------------
 * AccessStatus Routines
 *--------------------------------------------------------------------------*/

braid_Int
braid_AccessStatusGetT(braid_AccessStatus  status,
                       braid_Real         *t_ptr)
{
   *t_ptr = _braid_StatusElt(status, t);
   return _braid_error_flag;
}

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
braid_AccessStatusGetNRefine(braid_AccessStatus  status,
                             braid_Int          *nrefine_ptr)
{
   *nrefine_ptr = _braid_StatusElt(status, nrefine);
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_AccessStatusGetNTPoints(braid_AccessStatus  status,
                              braid_Int          *ntpoints_ptr)
{
   *ntpoints_ptr = _braid_StatusElt(status, gupper);
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
braid_AccessStatusGetTILD(braid_AccessStatus  status,
                          braid_Real          *t_ptr,
                          braid_Int           *iter_ptr,
                          braid_Int           *level_ptr,
                          braid_Int           *done_ptr)
{
   *t_ptr     =_braid_StatusElt(status, t);
   *level_ptr = _braid_StatusElt(status, level);
   *done_ptr  = _braid_StatusElt(status, done);
   *iter_ptr  = _braid_StatusElt(status, iter); 

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
braid_AccessStatusGetWrapperTest(braid_AccessStatus  status,
                                 braid_Int          *wtest_ptr
                                 )
{
   *wtest_ptr = _braid_StatusElt(status, wrapper_test);
   return _braid_error_flag;
}


/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/
braid_Int
_braid_AccessStatusInit(braid_Real           t,
                        braid_Real           rnorm,
                        braid_Int            iter,
                        braid_Int            level,
                        braid_Int            nrefine,
                        braid_Int            gupper,
                        braid_Int            done,
                        braid_Int            wrapper_test,
                        braid_AccessStatus   status)
{
   _braid_StatusElt(status, t)            = t;
   _braid_StatusElt(status, level)        = level;
   _braid_StatusElt(status, nrefine)      = nrefine;
   _braid_StatusElt(status, gupper)       = gupper;
   _braid_StatusElt(status, rnorm)        = rnorm;
   _braid_StatusElt(status, done)         = done;
   _braid_StatusElt(status, iter)         = iter; 
   _braid_StatusElt(status, wrapper_test) = wrapper_test;

   return _braid_error_flag;
}


/*--------------------------------------------------------------------------
 * CoarsenRefStatus Routines
 *--------------------------------------------------------------------------*/

braid_Int
_braid_CoarsenRefStatusInit(braid_Real              tstart,  
                            braid_Real              f_tprior,
                            braid_Real              f_tstop, 
                            braid_Real              c_tprior,
                            braid_Real              c_tstop,
                            braid_Int               level,
                            braid_Int               nrefine,
                            braid_Int               gupper,
                            braid_CoarsenRefStatus  status)
{
   _braid_StatusElt(status, tstart)   = tstart;
   _braid_StatusElt(status, f_tprior) = f_tprior;
   _braid_StatusElt(status, f_tstop)  = f_tstop;
   _braid_StatusElt(status, c_tprior) = c_tprior;
   _braid_StatusElt(status, c_tstop)  = c_tstop;
   _braid_StatusElt(status, level)    = level;
   _braid_StatusElt(status, nrefine)  = nrefine;
   _braid_StatusElt(status, gupper)   = gupper;

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
braid_CoarsenRefStatusGetFTstop(braid_CoarsenRefStatus  status,
                                braid_Real             *f_tstop_ptr
                                )
{
   *f_tstop_ptr = _braid_StatusElt(status, f_tstop);
   return _braid_error_flag;
}

braid_Int
braid_CoarsenRefStatusGetFTprior(braid_CoarsenRefStatus  status,
                                 braid_Real             *f_tprior_ptr
                                 )
{
   *f_tprior_ptr = _braid_StatusElt(status, f_tprior);
   return _braid_error_flag;
}

braid_Int
braid_CoarsenRefStatusGetCTstop(braid_CoarsenRefStatus  status,
                                braid_Real             *c_tstop_ptr
                                )
{
   *c_tstop_ptr = _braid_StatusElt(status, c_tstop);
   return _braid_error_flag;
}

braid_Int
braid_CoarsenRefStatusGetCTprior(braid_CoarsenRefStatus  status,
                                 braid_Real             *c_tprior_ptr
                                 )
{
   *c_tprior_ptr = _braid_StatusElt(status, c_tprior);
   return _braid_error_flag;
}

braid_Int
braid_CoarsenRefStatusGetTpriorTstop(braid_CoarsenRefStatus  status,
                                     braid_Real              *tstart_ptr,
                                     braid_Real              *f_tprior_ptr,
                                     braid_Real              *f_tstop_ptr, 
                                     braid_Real              *c_tprior_ptr,
                                     braid_Real              *c_tstop_ptr)
{
   *tstart_ptr = _braid_StatusElt(status, tstart);
   *f_tprior_ptr = _braid_StatusElt(status, f_tprior);
   *f_tstop_ptr = _braid_StatusElt(status, f_tstop);
   *c_tprior_ptr = _braid_StatusElt(status, c_tprior);
   *c_tstop_ptr = _braid_StatusElt(status, c_tstop);

   return _braid_error_flag;
}

braid_Int
braid_CoarsenRefStatusGetLevel(braid_CoarsenRefStatus  status,
                               braid_Int              *level_ptr
                               )
{
   *level_ptr = _braid_StatusElt(status, level);
   return _braid_error_flag;
}

braid_Int
braid_CoarsenRefStatusGetNRefine(braid_CoarsenRefStatus  status,
                                 braid_Int              *nrefine_ptr
                                )
{
   *nrefine_ptr = _braid_StatusElt(status, nrefine);
   return _braid_error_flag;
}

braid_Int
braid_CoarsenRefStatusGetNTPoints(braid_CoarsenRefStatus  status,
                                  braid_Int              *ntpoints_ptr)
{
   *ntpoints_ptr = _braid_StatusElt(status, gupper);
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * StepStatus Routines
 *--------------------------------------------------------------------------*/
braid_Int
_braid_StepStatusInit(braid_Real       tstart,
                      braid_Real       tstop,
                      braid_Real       tol,
                      braid_Int        iter,
                      braid_Int        level,
                      braid_Int        nrefine,
                      braid_Int        gupper,
                      braid_StepStatus status)
{
   _braid_StatusElt(status, tstart)    = tstart;
   _braid_StatusElt(status, tstop)     = tstop;
   _braid_StatusElt(status, tol)       = tol;
   _braid_StatusElt(status, iter)      = iter;
   _braid_StatusElt(status, level)     = level;
   _braid_StatusElt(status, nrefine)   = nrefine;
   _braid_StatusElt(status, gupper)    = gupper;
   _braid_StatusElt(status, rfactor)   = 1;
   _braid_StatusElt(status, r_space)   = 0;

   return _braid_error_flag;
}

braid_Int
_braid_StepStatusDestroy(braid_StepStatus  status)
{
   if (status)
   {
      _braid_TFree(status);
   }

   return _braid_error_flag;
}


braid_Int
braid_StepStatusGetTstart(braid_StepStatus  status,
                          braid_Real       *tstart_ptr)
{
   *tstart_ptr = _braid_StatusElt(status, tstart);
   return _braid_error_flag;
}

braid_Int
braid_StepStatusGetTstop(braid_StepStatus  status,
                         braid_Real       *tstop_ptr)
{
   *tstop_ptr = _braid_StatusElt(status, tstop);
   return _braid_error_flag;
}

braid_Int
braid_StepStatusGetLevel(braid_StepStatus  status,
                         braid_Int        *level_ptr
                        )
{
   *level_ptr = _braid_StatusElt(status, level);
   return _braid_error_flag;
}

braid_Int
braid_StepStatusGetNRefine(braid_StepStatus  status,
                           braid_Int        *nrefine_ptr
                          )
{
   *nrefine_ptr = _braid_StatusElt(status, nrefine);
   return _braid_error_flag;
}

braid_Int
braid_StepStatusGetNTPoints(braid_StepStatus  status,
                            braid_Int        *ntpoints_ptr)
{
   *ntpoints_ptr = _braid_StatusElt(status, gupper);
   return _braid_error_flag;
}

braid_Int
braid_StepStatusSetRFactor(braid_StepStatus  status,
                           braid_Real        rfactor)
{
   _braid_StatusElt(status, rfactor) = rfactor;
   return _braid_error_flag;
}

braid_Int
braid_StepStatusSetRSpace(braid_StepStatus  status,
                          braid_Int         r_space)
{
   _braid_StatusElt(status, r_space) = r_space;
   return _braid_error_flag;
}   

braid_Int
braid_StepStatusGetTstartTstop(braid_StepStatus  status,
                               braid_Real       *tstart_ptr,
                               braid_Real       *tstop_ptr)
{
   *tstart_ptr = _braid_StatusElt(status, tstart);
   *tstop_ptr = _braid_StatusElt(status, tstop);

   return _braid_error_flag;
}

braid_Int
braid_StepStatusGetTol(braid_StepStatus  status,
                       braid_Real       *tol_ptr
                       )
{
   *tol_ptr = _braid_StatusElt(status, tol);
   return _braid_error_flag;
}

braid_Int
braid_StepStatusGetIter(braid_StepStatus  status,
                        braid_Int        *iter_ptr
                        )
{
   *iter_ptr = _braid_StatusElt(status, iter);
   return _braid_error_flag;
}

braid_Int
braid_StepStatusGetRNorms(braid_StepStatus  status,
                          braid_Int        *nrequest_ptr,
                          braid_Real       *rnorms 
                          )
{
   braid_Real     *_rnorms   = _braid_StatusElt(status, rnorms);
   braid_Int      rnorms_len = _braid_StatusElt(status, iter) + 1;
   
   _braid_GetNEntries(_rnorms, rnorms_len, nrequest_ptr, rnorms);
   return _braid_error_flag;
}

braid_Int
braid_StepStatusGetOldFineTolx(braid_StepStatus  status,
                               braid_Real       *old_fine_tolx_ptr
                               )
{
   *old_fine_tolx_ptr = _braid_StatusElt(status, old_fine_tolx);
   return _braid_error_flag;
}

braid_Int
braid_StepStatusSetOldFineTolx(braid_StepStatus  status,
                               braid_Real        old_fine_tolx_ptr
                               )
{
   _braid_StatusElt(status, old_fine_tolx) = old_fine_tolx_ptr;
   return _braid_error_flag;
}

braid_Int
braid_StepStatusSetTightFineTolx(braid_StepStatus  status,
                                 braid_Int         tight_fine_tolx
                                 )
{
   _braid_StatusElt(status, tight_fine_tolx) = tight_fine_tolx;
   return _braid_error_flag;
}


/*--------------------------------------------------------------------------
 * BufferStatus Routines
 *--------------------------------------------------------------------------*/


braid_Int
_braid_BufferStatusInit(braid_Int        messagetype,
                        braid_Int        size,
                        braid_BufferStatus status)
{
   _braid_StatusElt(status, messagetype)    = messagetype;
   _braid_StatusElt(status, size)           = size;
   return _braid_error_flag;
}

braid_Int
_braid_BufferStatusDestroy(braid_BufferStatus  status)
{
   if (status)
   {
      _braid_TFree(status);
   }

   return _braid_error_flag;
}

braid_Int
braid_BufferStatusGetMessageType(braid_BufferStatus  status,
                                 braid_Int           *messagetype_ptr)
{
   *messagetype_ptr = _braid_StatusElt(status, messagetype);
   return _braid_error_flag;
}

braid_Int
braid_BufferStatusSetSize(braid_BufferStatus  status,
                          braid_Int           size)
{
   _braid_StatusElt(status, size ) = size; 
   return _braid_error_flag;
}

