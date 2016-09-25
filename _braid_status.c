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
 * \brief Source code for status interface routines.  See braid_status.h and
 * _braid_status.h for more information.
 *
 */

#include "_braid_status.h"
#include "_braid.h"
#include "braid_defs.h"
#include "_util.h"

#ifndef DEBUG
#define DEBUG 0
#endif

#define ACCESSOR_FUNCTION_GET1(stype,param,vtype1) \
   braid_Int braid_##stype##StatusGet##param(braid_##stype##Status s, braid_##vtype1 *v1) \
   {return braid_StatusGet##param((braid_Status)s, v1);}
#define ACCESSOR_FUNCTION_GET2(stype,param,vtype1,vtype2) \
   braid_Int braid_##stype##StatusGet##param(braid_##stype##Status s, braid_##vtype1 *v1, braid_##vtype2 *v2) \
   {return braid_StatusGet##param((braid_Status)s, v1, v2);}
#define ACCESSOR_FUNCTION_GET4(stype,param,vtype1,vtype2,vtype3,vtype4) \
   braid_Int braid_##stype##StatusGet##param(braid_##stype##Status s, braid_##vtype1 *v1, braid_##vtype2 *v2, braid_##vtype3 *v3, braid_##vtype4 *v4) \
   {return braid_StatusGet##param((braid_Status)s, v1, v2, v3, v4);}
#define ACCESSOR_FUNCTION_GET5(stype,param,vtype1,vtype2,vtype3,vtype4,vtype5) \
   braid_Int braid_##stype##StatusGet##param(braid_##stype##Status s, braid_##vtype1 *v1, braid_##vtype2 *v2, braid_##vtype3 *v3, braid_##vtype4 *v4, braid_##vtype5 *v5) \
   {return braid_StatusGet##param((braid_Status)s, v1, v2, v3, v4, v5);}
#define ACCESSOR_FUNCTION_SET1(stype,param,vtype1) \
   braid_Int braid_##stype##StatusSet##param(braid_##stype##Status s, braid_##vtype1 v1) \
   {return braid_StatusSet##param((braid_Status)s, v1);}


braid_Int
_braid_StatusDestroy(braid_Status status)
{
   if (status)
   {
      _braid_TFree(status);
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Status Routines
 *--------------------------------------------------------------------------*/

braid_Int
braid_StatusGetT(braid_Status status,
                 braid_Real  *t_ptr
                 )
{
   *t_ptr = _braid_StatusElt(status, t);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetTIndex(braid_Status status,
                     braid_Int    *idx_ptr
                     )
{
   *idx_ptr = _braid_StatusElt(status, idx);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetIter(braid_Status status,
                    braid_Int   *iter_ptr
                    )
{
   *iter_ptr = _braid_StatusElt(status, niter);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetLevel(braid_Status status,
                     braid_Int   *level_ptr
                     )
{
   *level_ptr = _braid_StatusElt(status, level);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetNRefine(braid_Status status,
                       braid_Int   *nrefine_ptr
                       )
{
   *nrefine_ptr = _braid_StatusElt(status, nrefine);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetNTPoints(braid_Status status,
                        braid_Int   *ntpoints_ptr
                        )
{
   *ntpoints_ptr = _braid_StatusElt(status, gupper);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetResidual(braid_Status status,
                        braid_Real  *rnorm_ptr
                        )
{
   *rnorm_ptr = _braid_StatusElt(status, rnorm);
   return _braid_error_flag;
}
braid_Int
braid_StatusGetDone(braid_Status status,
                    braid_Int   *done_ptr
                    )
{
   *done_ptr = _braid_StatusElt(status, done);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetTILD(braid_Status status,
                    braid_Real  *t_ptr,
                    braid_Int   *iter_ptr,
                    braid_Int   *level_ptr,
                    braid_Int   *done_ptr
                    )
{
   *t_ptr = _braid_StatusElt(status, t);
   *iter_ptr = _braid_StatusElt(status, niter);
   *level_ptr = _braid_StatusElt(status, level);
   *done_ptr = _braid_StatusElt(status, done);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetWrapperTest(braid_Status status,
                           braid_Int   *wtest_ptr
                           )
{
   *wtest_ptr = _braid_StatusElt(status, wrapper_test);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetCallingFunction(braid_Status status,
                               braid_Int   *cfunction_ptr
                               )
{
   *cfunction_ptr = _braid_StatusElt(status, calling_function);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetCTprior(braid_Status status,
                       braid_Real  *ctprior_ptr
                       )
{
   *ctprior_ptr = _braid_StatusElt(status, c_tprior);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetCTstop(braid_Status status,
                      braid_Real  *ctstop_ptr
                      )
{
   *ctstop_ptr = _braid_StatusElt(status, c_tstop);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetFTprior(braid_Status status,
                       braid_Real  *ftprior_ptr
                       )
{
   *ftprior_ptr = _braid_StatusElt(status, f_tprior);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetFTstop(braid_Status status,
                      braid_Real  *ftstop_ptr
                      )
{
   *ftstop_ptr = _braid_StatusElt(status, f_tstop);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetTpriorTstop(braid_Status status,
                           braid_Real  *t_ptr,
                           braid_Real  *ftprior_ptr,
                           braid_Real  *ftstop_ptr,
                           braid_Real  *ctprior_ptr,
                           braid_Real  *ctstop_ptr
                           )
{
   *t_ptr = _braid_StatusElt(status, t);
   *ctprior_ptr = _braid_StatusElt(status, c_tprior);
   *ctstop_ptr = _braid_StatusElt(status, c_tstop);
   *ftprior_ptr = _braid_StatusElt(status, f_tprior);
   *ftstop_ptr = _braid_StatusElt(status, f_tstop);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetTstop(braid_Status status,
                     braid_Real  *tstop_ptr
                     )
{
   *tstop_ptr = _braid_StatusElt(status, tnext);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetTstartTstop(braid_Status status,
                           braid_Real  *tstart_ptr,
                           braid_Real  *tstop_ptr
                           )
{
   *tstart_ptr = _braid_StatusElt(status, t);
   *tstop_ptr = _braid_StatusElt(status, tnext);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetTol(braid_Status status,
                   braid_Real  *tol_ptr
                   )
{
   *tol_ptr = _braid_StatusElt(status, tol);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetRNorms(braid_Status status,
                      braid_Int   *nrequest_ptr,
                      braid_Real  *rnorms_ptr
                      )
{
   braid_Real *_rnorms    = _braid_StatusElt(status, rnorms);
   braid_Int   rnorms_len = _braid_StatusElt(status, niter) + 1;

   _braid_GetNEntries(_rnorms, rnorms_len, nrequest_ptr, rnorms_ptr);
   return _braid_error_flag;
}

braid_Int
braid_StatusGetOldFineTolx(braid_Status status,
                           braid_Real  *old_fine_tolx_ptr
                           )
{
   *old_fine_tolx_ptr = _braid_StatusElt(status, old_fine_tolx);
   return _braid_error_flag;
}

braid_Int
braid_StatusSetOldFineTolx(braid_Status status,
                           braid_Real   old_fine_tolx
                           )
{
   _braid_StatusElt(status, old_fine_tolx) = old_fine_tolx;
   return _braid_error_flag;
}

braid_Int
braid_StatusSetTightFineTolx(braid_Status status,
                             braid_Real   tight_fine_tolx
                             )
{
   _braid_StatusElt(status, tight_fine_tolx) = tight_fine_tolx;
   return _braid_error_flag;
}

braid_Int
braid_StatusSetRFactor(braid_Status status,
                       braid_Real   rfactor
                       )
{
   _braid_StatusElt(status, rfactor) = rfactor;
   return _braid_error_flag;
}

braid_Int
braid_StatusSetRSpace(braid_Status status,
                      braid_Real   r_space
                      )
{
   if ( !_braid_StatusElt(status, r_space) && r_space )
      _braid_StatusElt(status, r_space) = 1;
   return _braid_error_flag;
}

braid_Int
braid_StatusGetMessageType(braid_Status status,
                           braid_Int   *messagetype_ptr
                           )
{
   *messagetype_ptr = _braid_StatusElt(status, messagetype);
   return _braid_error_flag;
}

braid_Int
braid_StatusSetSize(braid_Status status,
                    braid_Real   size
                    )
{
   _braid_StatusElt(status, size_buffer ) = size;
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * AccessStatus Routines
 *--------------------------------------------------------------------------*/

braid_Int
_braid_AccessStatusInit(braid_Real           t,
                        braid_Int            idx,
                        braid_Real           rnorm,
                        braid_Int            iter,
                        braid_Int            level,
                        braid_Int            nrefine,
                        braid_Int            gupper,
                        braid_Int            done,
                        braid_Int            wrapper_test,
                        braid_Int            calling_function,
                        braid_AccessStatus   status)
{
   _braid_StatusElt(status, t)            = t;
   _braid_StatusElt(status, idx)          = idx;
   _braid_StatusElt(status, level)        = level;
   _braid_StatusElt(status, nrefine)      = nrefine;
   _braid_StatusElt(status, gupper)       = gupper;
   _braid_StatusElt(status, rnorm)        = rnorm;
   _braid_StatusElt(status, done)         = done;
   _braid_StatusElt(status, niter)        = iter;
   _braid_StatusElt(status, wrapper_test) = wrapper_test;
   _braid_StatusElt(status, calling_function) = calling_function;
   return _braid_error_flag;
}
ACCESSOR_FUNCTION_GET1(Access, T,               Real)
ACCESSOR_FUNCTION_GET1(Access, TIndex,          Int)
ACCESSOR_FUNCTION_GET1(Access, Iter,            Int)
ACCESSOR_FUNCTION_GET1(Access, Level,           Int)
ACCESSOR_FUNCTION_GET1(Access, NRefine,         Int)
ACCESSOR_FUNCTION_GET1(Access, NTPoints,        Int)
ACCESSOR_FUNCTION_GET1(Access, Residual,        Real)
ACCESSOR_FUNCTION_GET1(Access, Done,            Int)
ACCESSOR_FUNCTION_GET4(Access, TILD,            Real, Int, Int, Int)
ACCESSOR_FUNCTION_GET1(Access, WrapperTest,     Int)
ACCESSOR_FUNCTION_GET1(Access, CallingFunction, Int)

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
   _braid_StatusElt(status, t)        = tstart;
   _braid_StatusElt(status, f_tprior) = f_tprior;
   _braid_StatusElt(status, f_tstop)  = f_tstop;
   _braid_StatusElt(status, c_tprior) = c_tprior;
   _braid_StatusElt(status, c_tstop)  = c_tstop;
   _braid_StatusElt(status, level)    = level;
   _braid_StatusElt(status, nrefine)  = nrefine;
   _braid_StatusElt(status, gupper)   = gupper;

   return _braid_error_flag;
}
ACCESSOR_FUNCTION_GET1(CoarsenRef, T,           Real)
ACCESSOR_FUNCTION_GET1(CoarsenRef, TIndex,      Int)
ACCESSOR_FUNCTION_GET1(CoarsenRef, Iter,        Int)
ACCESSOR_FUNCTION_GET1(CoarsenRef, Level,       Int)
ACCESSOR_FUNCTION_GET1(CoarsenRef, NRefine,     Int)
ACCESSOR_FUNCTION_GET1(CoarsenRef, NTPoints,    Int)
ACCESSOR_FUNCTION_GET1(CoarsenRef, CTprior,     Real)
ACCESSOR_FUNCTION_GET1(CoarsenRef, CTstop,      Real)
ACCESSOR_FUNCTION_GET1(CoarsenRef, FTprior,     Real)
ACCESSOR_FUNCTION_GET1(CoarsenRef, FTstop,      Real)
ACCESSOR_FUNCTION_GET5(CoarsenRef, TpriorTstop, Real, Real, Real, Real, Real)

/*--------------------------------------------------------------------------
 * StepStatus Routines
 *--------------------------------------------------------------------------*/

braid_Int
_braid_StepStatusInit(braid_Real       tstart,
                      braid_Real       tstop,
                      braid_Int        idx,
                      braid_Real       tol,
                      braid_Int        iter,
                      braid_Int        level,
                      braid_Int        nrefine,
                      braid_Int        gupper,
                      braid_StepStatus status)
{
   _braid_StatusElt(status, t)         = tstart;
   _braid_StatusElt(status, tnext)     = tstop;
   _braid_StatusElt(status, idx)       = idx;
   _braid_StatusElt(status, tol)       = tol;
   _braid_StatusElt(status, niter)     = iter;
   _braid_StatusElt(status, level)     = level;
   _braid_StatusElt(status, nrefine)   = nrefine;
   _braid_StatusElt(status, gupper)    = gupper;
   _braid_StatusElt(status, rfactor)   = 1;
   _braid_StatusElt(status, r_space)   = 0;

   return _braid_error_flag;
}
ACCESSOR_FUNCTION_GET1(Step, T,             Real)
ACCESSOR_FUNCTION_GET1(Step, TIndex,        Int)
ACCESSOR_FUNCTION_GET1(Step, Iter,          Int)
ACCESSOR_FUNCTION_GET1(Step, Level,         Int)
ACCESSOR_FUNCTION_GET1(Step, NRefine,       Int)
ACCESSOR_FUNCTION_GET1(Step, NTPoints,      Int)
ACCESSOR_FUNCTION_GET1(Step, Tstop,         Real)
ACCESSOR_FUNCTION_GET2(Step, TstartTstop,   Real, Real)
ACCESSOR_FUNCTION_GET1(Step, Tol,           Real)
ACCESSOR_FUNCTION_GET2(Step, RNorms,        Int,  Real)
ACCESSOR_FUNCTION_GET1(Step, OldFineTolx,   Real)
ACCESSOR_FUNCTION_SET1(Step, OldFineTolx,   Real)
ACCESSOR_FUNCTION_SET1(Step, TightFineTolx, Real)
ACCESSOR_FUNCTION_SET1(Step, RFactor,       Real)
ACCESSOR_FUNCTION_SET1(Step, RSpace,        Real)

/*--------------------------------------------------------------------------
 * BufferStatus Routines
 *--------------------------------------------------------------------------*/

braid_Int
_braid_BufferStatusInit(braid_Int        messagetype,
                        braid_Int        size,
                        braid_BufferStatus status)
{
   _braid_StatusElt(status, messagetype)    = messagetype;
   _braid_StatusElt(status, size_buffer)    = size;
   return _braid_error_flag;
}
ACCESSOR_FUNCTION_GET1(Buffer, MessageType, Int)
ACCESSOR_FUNCTION_SET1(Buffer, Size,        Real)
