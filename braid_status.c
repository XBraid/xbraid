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
 * GlobalStatus Routines
 *--------------------------------------------------------------------------*/

braid_Int
braid_GlobalStatusGetT(braid_GlobalStatus status,                      /**< structure containing current simulation info */
                       braid_Real        *t_ptr                        /**< output, current time */
                       )
{
   *t_ptr = _braid_StatusElt(status->gs, t);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetTIndex(braid_GlobalStatus status,                  /**< structure containing current simulation info */
                           braid_Int         *idx_ptr                  /**< output, global index value corresponding to current time value */
                           )
{
   *idx_ptr = _braid_StatusElt(status->gs, idx);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetIter(braid_GlobalStatus status,                   /**< structure containing current simulation info */
                          braid_Int         *iter_ptr                  /**< output, current XBraid iteration number*/
                          )
{
   *iter_ptr = _braid_StatusElt(status->gs, iter);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetLevel(braid_GlobalStatus status,                  /**< structure containing current simulation info */
                           braid_Int         *level_ptr                /**< output, current level in XBraid */
                           )
{
   *level_ptr = _braid_StatusElt(status->gs, level);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetNRefine(braid_GlobalStatus status,                /**< structure containing current simulation info */
                             braid_Int         *nrefine_ptr            /**< output, number of refinements done */
                             )
{
   *nrefine_ptr = _braid_StatusElt(status->gs, nrefine);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetNTPoints(braid_GlobalStatus status,               /**< structure containing current simulation info */
                              braid_Int         *ntpoints_ptr          /**< output, number of time points on the fine grid */
                              )
{
   *ntpoints_ptr = _braid_StatusElt(status->gs, gupper);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetResidual(braid_GlobalStatus status,               /**< structure containing current simulation info */
                              braid_Real        *rnorm_ptr             /**< output, current residual norm */
                              )
{
   *rnorm_ptr = _braid_StatusElt(status->gs, rnorm);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetDone(braid_GlobalStatus status,                   /**< structure containing current simulation info */
                          braid_Int         *done_ptr                  /**< output,  =1 if XBraid has finished, else =0 */
                          )
{
   *done_ptr = _braid_StatusElt(status->gs, done);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetTILD(braid_GlobalStatus status,                   /**< structure containing current simulation info */
                          braid_Real        *t_ptr,                    /**< output, current time */
                          braid_Int         *iter_ptr,                 /**< output, current XBraid iteration number*/
                          braid_Int         *level_ptr,                /**< output, current level in XBraid */
                          braid_Int         *done_ptr                  /**< output,  =1 if XBraid has finished, else =0 */
                          )
{
   *t_ptr = _braid_StatusElt(status->gs, t);
   *iter_ptr = _braid_StatusElt(status->gs, iter);
   *level_ptr = _braid_StatusElt(status->gs, level);
   *done_ptr = _braid_StatusElt(status->gs, done);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetWrapperTest(braid_GlobalStatus status,            /**< structure containing current simulation info */
                                 braid_Int         *wtest_ptr          /**< output, =1 if this is a wrapper test, =0 if XBraid run */
                                 )
{
   *wtest_ptr = _braid_StatusElt(status->gs, wrapper_test);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetCallingFunction(braid_GlobalStatus status,        /**< structure containing current simulation info */
                                     braid_Int         *cfunction_ptr  /**< output, function number (0=FInterp, 1=FRestrict, 2=FRefine, 3=FAccess) */
                                     )
{
   *cfunction_ptr = _braid_StatusElt(status->gs, calling_function);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetCTprior(braid_GlobalStatus status,                /**< structure containing current simulation info */
                             braid_Real        *ctprior_ptr            /**< output, time value to the left of current time value on coarse grid */
                             )
{
   *ctprior_ptr = _braid_StatusElt(status->gs, c_tprior);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetCTstop(braid_GlobalStatus status,                 /**< structure containing current simulation info */
                            braid_Real        *ctstop_ptr              /**< output, time value to the right of current time value on coarse grid */
                            )
{
   *ctstop_ptr = _braid_StatusElt(status->gs, c_tstop);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetFTprior(braid_GlobalStatus status,                /**< structure containing current simulation info */
                             braid_Real        *ftprior_ptr            /**< output, time value to the left of current time value on fine grid */
                             )
{
   *ftprior_ptr = _braid_StatusElt(status->gs, f_tprior);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetFTstop(braid_GlobalStatus status,                 /**< structure containing current simulation info */
                            braid_Real        *ftstop_ptr              /**< output, time value to the right of current time value on fine grid */
                            )
{
   *ftstop_ptr = _braid_StatusElt(status->gs, f_tstop);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetTpriorTstop(braid_GlobalStatus status,            /**< structure containing current simulation info */
                                 braid_Real        *t_ptr,             /**< output, current time */
                                 braid_Real        *ftprior_ptr,       /**< output, time value to the left of current time value on fine grid */
                                 braid_Real        *ftstop_ptr,        /**< output, time value to the right of current time value on fine grid */
                                 braid_Real        *ctprior_ptr,       /**< output, time value to the left of current time value on coarse grid */
                                 braid_Real        *ctstop_ptr         /**< output, time value to the right of current time value on coarse grid */
                                 )
{
   *ctprior_ptr = _braid_StatusElt(status->gs, c_tprior);
   *ctstop_ptr = _braid_StatusElt(status->gs, c_tstop);
   *ftprior_ptr = _braid_StatusElt(status->gs, f_tprior);
   *ftstop_ptr = _braid_StatusElt(status->gs, f_tstop);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetTstop(braid_GlobalStatus status,                  /**< structure containing current simulation info */
                           braid_Real        *tstop_ptr                /**< output, next time value to evolve towards */
                           )
{
   *tstop_ptr = _braid_StatusElt(status->gs, tstop);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetTstartTstop(braid_GlobalStatus status,            /**< structure containing current simulation info */
                                 braid_Real        *tstart_ptr,        /**< output, current time */
                                 braid_Real        *tstop_ptr          /**< output, next time value to evolve towards */
                                 )
{
   *tstart_ptr = _braid_StatusElt(status->gs, t);
   *tstop_ptr = _braid_StatusElt(status->gs, tstop);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetTol(braid_GlobalStatus status,                    /**< structure containing current simulation info */
                         braid_Real        *tol_ptr                    /**< output, current XBraid stopping tolerance */
                         )
{
   *tol_ptr = _braid_StatusElt(status->gs, tol);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetRNorms(braid_GlobalStatus status,                 /**< structure containing current simulation info */
                            braid_Int         *nrequest_ptr,           /**< input/output, input: number of requested residual norms, output: number actually copied */
                            braid_Real        *rnorms_ptr              /**< output, XBraid residual norm history, of length *nrequest_ptr* */
                            )
{
   braid_Real *_rnorms    = _braid_StatusElt(status->gs, rnorms);
   braid_Int   rnorms_len = _braid_StatusElt(status->gs, iter) + 1;

   _braid_GetNEntries(_rnorms, rnorms_len, nrequest_ptr, rnorms_ptr);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetOldFineTolx(braid_GlobalStatus status,            /**< structure containing current simulation info */
                                 braid_Real        *old_fine_tolx_ptr  /**< output, previous *old_fine_tolx*, set through *braid_StepStatusSetOldFineTolx* */
                                 )
{
   *old_fine_tolx_ptr = _braid_StatusElt(status->gs, old_fine_tolx);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusSetOldFineTolx(braid_GlobalStatus status,            /**< structure containing current simulation info */
                                 braid_Real         old_fine_tolx      /**< input, the last used fine_tolx */
                                 )
{
   _braid_StatusElt(status->gs, old_fine_tolx) = old_fine_tolx;
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusSetTightFineTolx(braid_GlobalStatus status,          /**< structure containing current simulation info */
                                   braid_Real         tight_fine_tolx  /**< input, boolean indicating whether the tight tolx has been used */
                                   )
{
   _braid_StatusElt(status->gs, tight_fine_tolx) = tight_fine_tolx;
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusSetRFactor(braid_GlobalStatus status,                /**< structure containing current simulation info */
                             braid_Real         rfactor                /**< input, user-determined desired rfactor */
                             )
{
   _braid_StatusElt(status->gs, rfactor) = rfactor;
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusSetRSpace(braid_GlobalStatus status,                 /**< structure containing current simulation info */
                            braid_Real         r_space
                            )
{
   _braid_StatusElt(status->gs, r_space) = r_space;
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusGetMessageType(braid_GlobalStatus status,            /**< structure containing current simulation info */
                                 braid_Int         *messagetype_ptr    /**< output, type of message, 0: for Step(), 1: for load balancing */
                                 )
{
   *messagetype_ptr = _braid_StatusElt(status->gs, messagetype);
   return _braid_error_flag;
}
braid_Int
braid_GlobalStatusSetSize(braid_GlobalStatus status,                   /**< structure containing current simulation info */
                          braid_Real         size                      /**< input, size of the send buffer */
                          )
{
   _braid_StatusElt(status->gs, size ) = size;
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
   _braid_StatusElt(status->gs, t)            = t;
   _braid_StatusElt(status->gs, idx)          = idx;
   _braid_StatusElt(status->gs, level)        = level;
   _braid_StatusElt(status->gs, nrefine)      = nrefine;
   _braid_StatusElt(status->gs, gupper)       = gupper;
   _braid_StatusElt(status->gs, rnorm)        = rnorm;
   _braid_StatusElt(status->gs, done)         = done;
   _braid_StatusElt(status->gs, iter)         = iter;
   _braid_StatusElt(status->gs, wrapper_test) = wrapper_test;
   _braid_StatusElt(status->gs, calling_function) = calling_function;
   return _braid_error_flag;
}
braid_Int
braid_AccessStatusGetT(braid_AccessStatus status,                      /**< structure containing current simulation info */
                       braid_Real        *t_ptr                        /**< output, current time */
                       )
{
   return braid_GlobalStatusGetT((braid_GlobalStatus)status, t_ptr);
}
braid_Int
braid_AccessStatusGetTIndex(braid_AccessStatus status,                  /**< structure containing current simulation info */
                           braid_Int         *idx_ptr                  /**< output, global index value corresponding to current time value */
                           )
{
   return braid_GlobalStatusGetTIndex((braid_GlobalStatus)status, idx_ptr);
}
braid_Int
braid_AccessStatusGetIter(braid_AccessStatus status,                   /**< structure containing current simulation info */
                          braid_Int         *iter_ptr                  /**< output, current XBraid iteration number*/
                          )
{
   return braid_GlobalStatusGetIter((braid_GlobalStatus)status, iter_ptr);
}
braid_Int
braid_AccessStatusGetLevel(braid_AccessStatus status,                  /**< structure containing current simulation info */
                           braid_Int         *level_ptr                /**< output, current level in XBraid */
                           )
{
   return braid_GlobalStatusGetLevel((braid_GlobalStatus)status, level_ptr);
}
braid_Int
braid_AccessStatusGetNRefine(braid_AccessStatus status,                /**< structure containing current simulation info */
                             braid_Int         *nrefine_ptr            /**< output, number of refinements done */
                             )
{
   return braid_GlobalStatusGetNRefine((braid_GlobalStatus)status, nrefine_ptr);
}
braid_Int
braid_AccessStatusGetNTPoints(braid_AccessStatus status,               /**< structure containing current simulation info */
                              braid_Int         *ntpoints_ptr          /**< output, number of time points on the fine grid */
                              )
{
   return braid_GlobalStatusGetNTPoints((braid_GlobalStatus)status, ntpoints_ptr);
}
braid_Int
braid_AccessStatusGetResidual(braid_AccessStatus status,               /**< structure containing current simulation info */
                              braid_Real        *rnorm_ptr             /**< output, current residual norm */
                              )
{
   return braid_GlobalStatusGetResidual((braid_GlobalStatus)status, rnorm_ptr);
}
braid_Int
braid_AccessStatusGetDone(braid_AccessStatus status,                   /**< structure containing current simulation info */
                          braid_Int         *done_ptr                  /**< output,  =1 if XBraid has finished, else =0 */
                          )
{
   return braid_GlobalStatusGetDone((braid_GlobalStatus)status, done_ptr);
}
braid_Int
braid_AccessStatusGetTILD(braid_AccessStatus status,                   /**< structure containing current simulation info */
                          braid_Real        *t_ptr,                    /**< output, current time */
                          braid_Int         *iter_ptr,                 /**< output, current XBraid iteration number*/
                          braid_Int         *level_ptr,                /**< output, current level in XBraid */
                          braid_Int         *done_ptr                  /**< output,  =1 if XBraid has finished, else =0 */
                          )
{
   return braid_GlobalStatusGetTILD((braid_GlobalStatus)status, t_ptr, iter_ptr, level_ptr, done_ptr);
}
braid_Int
braid_AccessStatusGetWrapperTest(braid_AccessStatus status,            /**< structure containing current simulation info */
                                 braid_Int         *wtest_ptr          /**< output, =1 if this is a wrapper test, =0 if XBraid run */
                                 )
{
   return braid_GlobalStatusGetWrapperTest((braid_GlobalStatus)status, wtest_ptr);
}
braid_Int
braid_AccessStatusGetCallingFunction(braid_AccessStatus status,        /**< structure containing current simulation info */
                                     braid_Int         *cfunction_ptr  /**< output, function number (0=FInterp, 1=FRestrict, 2=FRefine, 3=FAccess) */
                                     )
{
   return braid_GlobalStatusGetCallingFunction((braid_GlobalStatus)status, cfunction_ptr);
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
   _braid_StatusElt(status->gs, t)        = tstart;
   _braid_StatusElt(status->gs, f_tprior) = f_tprior;
   _braid_StatusElt(status->gs, f_tstop)  = f_tstop;
   _braid_StatusElt(status->gs, c_tprior) = c_tprior;
   _braid_StatusElt(status->gs, c_tstop)  = c_tstop;
   _braid_StatusElt(status->gs, level)    = level;
   _braid_StatusElt(status->gs, nrefine)  = nrefine;
   _braid_StatusElt(status->gs, gupper)   = gupper;

   return _braid_error_flag;
}
braid_Int
braid_CoarsenRefStatusGetT(braid_CoarsenRefStatus status,                      /**< structure containing current simulation info */
                           braid_Real            *t_ptr                        /**< output, current time */
                           )
{
   return braid_GlobalStatusGetT((braid_GlobalStatus)status, t_ptr);
}
braid_Int
braid_CoarsenRefStatusGetTIndex(braid_CoarsenRefStatus status,                  /**< structure containing current simulation info */
                               braid_Int             *idx_ptr                  /**< output, global index value corresponding to current time value */
                               )
{
   return braid_GlobalStatusGetTIndex((braid_GlobalStatus)status, idx_ptr);
}
braid_Int
braid_CoarsenRefStatusGetIter(braid_CoarsenRefStatus status,                   /**< structure containing current simulation info */
                              braid_Int             *iter_ptr                  /**< output, current XBraid iteration number*/
                              )
{
   return braid_GlobalStatusGetIter((braid_GlobalStatus)status, iter_ptr);
}
braid_Int
braid_CoarsenRefStatusGetLevel(braid_CoarsenRefStatus status,                  /**< structure containing current simulation info */
                               braid_Int             *level_ptr                /**< output, current level in XBraid */
                               )
{
   return braid_GlobalStatusGetLevel((braid_GlobalStatus)status, level_ptr);
}
braid_Int
braid_CoarsenRefStatusGetNRefine(braid_CoarsenRefStatus status,                /**< structure containing current simulation info */
                                 braid_Int             *nrefine_ptr            /**< output, number of refinements done */
                                 )
{
   return braid_GlobalStatusGetNRefine((braid_GlobalStatus)status, nrefine_ptr);
}
braid_Int
braid_CoarsenRefStatusGetNTPoints(braid_CoarsenRefStatus status,               /**< structure containing current simulation info */
                                  braid_Int             *ntpoints_ptr          /**< output, number of time points on the fine grid */
                                  )
{
   return braid_GlobalStatusGetNTPoints((braid_GlobalStatus)status, ntpoints_ptr);
}
braid_Int
braid_CoarsenRefStatusGetCTprior(braid_CoarsenRefStatus status,                /**< structure containing current simulation info */
                                 braid_Real            *ctprior_ptr            /**< output, time value to the left of current time value on coarse grid */
                                 )
{
   return braid_GlobalStatusGetCTprior((braid_GlobalStatus)status, ctprior_ptr);
}
braid_Int
braid_CoarsenRefStatusGetCTstop(braid_CoarsenRefStatus status,                 /**< structure containing current simulation info */
                                braid_Real            *ctstop_ptr              /**< output, time value to the right of current time value on coarse grid */
                                )
{
   return braid_GlobalStatusGetCTstop((braid_GlobalStatus)status, ctstop_ptr);
}
braid_Int
braid_CoarsenRefStatusGetFTprior(braid_CoarsenRefStatus status,                /**< structure containing current simulation info */
                                 braid_Real            *ftprior_ptr            /**< output, time value to the left of current time value on fine grid */
                                 )
{
   return braid_GlobalStatusGetFTprior((braid_GlobalStatus)status, ftprior_ptr);
}
braid_Int
braid_CoarsenRefStatusGetFTstop(braid_CoarsenRefStatus status,                 /**< structure containing current simulation info */
                                braid_Real            *ftstop_ptr              /**< output, time value to the right of current time value on fine grid */
                                )
{
   return braid_GlobalStatusGetFTstop((braid_GlobalStatus)status, ftstop_ptr);
}
braid_Int
braid_CoarsenRefStatusGetTpriorTstop(braid_CoarsenRefStatus status,            /**< structure containing current simulation info */
                                     braid_Real            *t_ptr,             /**< output, current time */
                                     braid_Real            *ftprior_ptr,       /**< output, time value to the left of current time value on fine grid */
                                     braid_Real            *ftstop_ptr,        /**< output, time value to the right of current time value on fine grid */
                                     braid_Real            *ctprior_ptr,       /**< output, time value to the left of current time value on coarse grid */
                                     braid_Real            *ctstop_ptr         /**< output, time value to the right of current time value on coarse grid */
                                     )
{
   return braid_GlobalStatusGetTpriorTstop((braid_GlobalStatus)status, t_ptr, ftprior_ptr, ftstop_ptr, ctprior_ptr, ctstop_ptr);
}

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
   _braid_StatusElt(status->gs, t)         = tstart;
   _braid_StatusElt(status->gs, tstop)     = tstop;
   _braid_StatusElt(status->gs, idx)       = idx;
   _braid_StatusElt(status->gs, tol)       = tol;
   _braid_StatusElt(status->gs, iter)      = iter;
   _braid_StatusElt(status->gs, level)     = level;
   _braid_StatusElt(status->gs, nrefine)   = nrefine;
   _braid_StatusElt(status->gs, gupper)    = gupper;
   _braid_StatusElt(status->gs, rfactor)   = 1;
   _braid_StatusElt(status->gs, r_space)   = 0;

   return _braid_error_flag;
}
braid_Int
braid_StepStatusGetT(braid_StepStatus status,                      /**< structure containing current simulation info */
                     braid_Real      *t_ptr                        /**< output, current time */
                     )
{
   return braid_GlobalStatusGetT((braid_GlobalStatus)status, t_ptr);
}
braid_Int
braid_StepStatusGetTIndex(braid_StepStatus status,                  /**< structure containing current simulation info */
                         braid_Int       *idx_ptr                  /**< output, global index value corresponding to current time value */
                         )
{
   return braid_GlobalStatusGetTIndex((braid_GlobalStatus)status, idx_ptr);
}
braid_Int
braid_StepStatusGetIter(braid_StepStatus status,                   /**< structure containing current simulation info */
                        braid_Int       *iter_ptr                  /**< output, current XBraid iteration number*/
                        )
{
   return braid_GlobalStatusGetIter((braid_GlobalStatus)status, iter_ptr);
}
braid_Int
braid_StepStatusGetLevel(braid_StepStatus status,                  /**< structure containing current simulation info */
                         braid_Int       *level_ptr                /**< output, current level in XBraid */
                         )
{
   return braid_GlobalStatusGetLevel((braid_GlobalStatus)status, level_ptr);
}
braid_Int
braid_StepStatusGetNRefine(braid_StepStatus status,                /**< structure containing current simulation info */
                           braid_Int       *nrefine_ptr            /**< output, number of refinements done */
                           )
{
   return braid_GlobalStatusGetNRefine((braid_GlobalStatus)status, nrefine_ptr);
}
braid_Int
braid_StepStatusGetNTPoints(braid_StepStatus status,               /**< structure containing current simulation info */
                            braid_Int       *ntpoints_ptr          /**< output, number of time points on the fine grid */
                            )
{
   return braid_GlobalStatusGetNTPoints((braid_GlobalStatus)status, ntpoints_ptr);
}
braid_Int
braid_StepStatusGetTstop(braid_StepStatus status,                  /**< structure containing current simulation info */
                         braid_Real      *tstop_ptr                /**< output, next time value to evolve towards */
                         )
{
   return braid_GlobalStatusGetTstop((braid_GlobalStatus)status, tstop_ptr);
}
braid_Int
braid_StepStatusGetTstartTstop(braid_StepStatus status,            /**< structure containing current simulation info */
                               braid_Real      *tstart_ptr,        /**< output, current time */
                               braid_Real      *tstop_ptr          /**< output, next time value to evolve towards */
                               )
{
   return braid_GlobalStatusGetTstartTstop((braid_GlobalStatus)status, tstart_ptr, tstop_ptr);
}
braid_Int
braid_StepStatusGetTol(braid_StepStatus status,                    /**< structure containing current simulation info */
                       braid_Real      *tol_ptr                    /**< output, current XBraid stopping tolerance */
                       )
{
   return braid_GlobalStatusGetTol((braid_GlobalStatus)status, tol_ptr);
}
braid_Int
braid_StepStatusGetRNorms(braid_StepStatus status,                 /**< structure containing current simulation info */
                          braid_Int       *nrequest_ptr,           /**< input/output, input: number of requested residual norms, output: number actually copied */
                          braid_Real      *rnorms_ptr              /**< output, XBraid residual norm history, of length *nrequest_ptr* */
                          )
{
   return braid_GlobalStatusGetRNorms((braid_GlobalStatus)status, nrequest_ptr, rnorms_ptr);
}
braid_Int
braid_StepStatusGetOldFineTolx(braid_StepStatus status,            /**< structure containing current simulation info */
                               braid_Real      *old_fine_tolx_ptr  /**< output, previous *old_fine_tolx*, set through *braid_StepStatusSetOldFineTolx* */
                               )
{
   return braid_GlobalStatusGetOldFineTolx((braid_GlobalStatus)status, old_fine_tolx_ptr);
}
braid_Int
braid_StepStatusSetOldFineTolx(braid_StepStatus status,            /**< structure containing current simulation info */
                               braid_Real       old_fine_tolx      /**< input, the last used fine_tolx */
                               )
{
   return braid_GlobalStatusSetOldFineTolx((braid_GlobalStatus)status, old_fine_tolx);
}
braid_Int
braid_StepStatusSetTightFineTolx(braid_StepStatus status,          /**< structure containing current simulation info */
                                 braid_Real       tight_fine_tolx  /**< input, boolean indicating whether the tight tolx has been used */
                                 )
{
   return braid_GlobalStatusSetTightFineTolx((braid_GlobalStatus)status, tight_fine_tolx);
}
braid_Int
braid_StepStatusSetRFactor(braid_StepStatus status,                /**< structure containing current simulation info */
                           braid_Real       rfactor                /**< input, user-determined desired rfactor */
                           )
{
   return braid_GlobalStatusSetRFactor((braid_GlobalStatus)status, rfactor);
}
braid_Int
braid_StepStatusSetRSpace(braid_StepStatus status,                 /**< structure containing current simulation info */
                          braid_Real       r_space
                          )
{
   return braid_GlobalStatusSetRSpace((braid_GlobalStatus)status, r_space);
}

/*--------------------------------------------------------------------------
 * BufferStatus Routines
 *--------------------------------------------------------------------------*/

braid_Int
_braid_BufferStatusInit(braid_Int        messagetype,
                        braid_Int        size,
                        braid_BufferStatus status)
{
   _braid_StatusElt(status->gs, messagetype)    = messagetype;
   _braid_StatusElt(status->gs, size)           = size;
   return _braid_error_flag;
}
braid_Int
braid_BufferStatusGetMessageType(braid_BufferStatus status,            /**< structure containing current simulation info */
                                 braid_Int         *messagetype_ptr    /**< output, type of message, 0: for Step(), 1: for load balancing */
                                 )
{
   return braid_GlobalStatusGetMessageType((braid_GlobalStatus)status, messagetype_ptr);
}
braid_Int
braid_BufferStatusSetSize(braid_BufferStatus status,                   /**< structure containing current simulation info */
                          braid_Real         size                      /**< input, size of the send buffer */
                          )
{
   return braid_GlobalStatusSetSize((braid_GlobalStatus)status, size);
}
