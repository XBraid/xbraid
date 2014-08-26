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
_braid_AccessStatusInit(braid_Real           t,
                        braid_Real           rnorm,
                        braid_Int            iter,
                        braid_Int            level,
                        braid_Int            done,
                        braid_AccessStatus   status)
{
   _braid_StatusElt(status, t)     = t;
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
                            braid_Real              f_tprior,
                            braid_Real              f_tstop, 
                            braid_Real              c_tprior,
                            braid_Real              c_tstop,
                            braid_CoarsenRefStatus  status)
{
   _braid_StatusElt(status, tstart) = tstart;
   _braid_StatusElt(status, f_tprior) = f_tprior;
   _braid_StatusElt(status, f_tstop)  = f_tstop;
   _braid_StatusElt(status, c_tprior) = c_tprior;
   _braid_StatusElt(status, c_tstop)  = c_tstop;

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


/*--------------------------------------------------------------------------
 * CoarsenRefStatus Routines
 *--------------------------------------------------------------------------*/
braid_Int
_braid_PhiStatusInit(braid_Real       tstart,
                     braid_Real       tstop,
                     braid_Real       accuracy,
                     braid_PhiStatus  status)
{
   _braid_StatusElt(status, tstart)   = tstart;
   _braid_StatusElt(status, tstop)    = tstop;
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
braid_PhiStatusGetTstop(braid_PhiStatus  status,
                        braid_Real      *tstop_ptr)
{
   *tstop_ptr = _braid_StatusElt(status, tstop);
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

braid_Int
braid_PhiStatusGetTstartTstop(braid_PhiStatus  status,
                              braid_Real       *tstart_ptr,
                              braid_Real       *tstop_ptr)
{
   *tstart_ptr = _braid_StatusElt(status, tstart);
   *tstop_ptr = _braid_StatusElt(status, tstop);

   return _braid_error_flag;
}


