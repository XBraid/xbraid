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
 

/** \file braid_status.h
 * \brief Define headers for XBraid status structures and headers for the user 
 * functions allowing the user to get/set status structure values. 
 *
 */

#ifndef braid_status_HEADER
#define braid_status_HEADER

#include "braid_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------------
 * These accessor macros allow for a "generic" braid_StatusGet function to be
 * used by all the various Status structures, e.g., StepStatus and AccessStatus
 *--------------------------------------------------------------------------*/
/** Macros allowing for auto-generation of `inherited' StatusGet functions */
#define ACCESSOR_HEADER_GET1(stype,param,vtype1) \
  braid_Int braid_##stype##StatusGet##param(braid_##stype##Status s, braid_##vtype1 *v1);
#define ACCESSOR_HEADER_GET2(stype,param,vtype1,vtype2) \
  braid_Int braid_##stype##StatusGet##param(braid_##stype##Status s, braid_##vtype1 *v1, braid_##vtype2 *v2);
#define ACCESSOR_HEADER_GET4(stype,param,vtype1,vtype2,vtype3,vtype4) \
  braid_Int braid_##stype##StatusGet##param(braid_##stype##Status s, braid_##vtype1 *v1, braid_##vtype2 *v2, braid_##vtype3 *v3, braid_##vtype4 *v4);
#define ACCESSOR_HEADER_GET5(stype,param,vtype1,vtype2,vtype3,vtype4,vtype5) \
  braid_Int braid_##stype##StatusGet##param(braid_##stype##Status s, braid_##vtype1 *v1, braid_##vtype2 *v2, braid_##vtype3 *v3, braid_##vtype4 *v4, braid_##vtype5 *v5);
#define ACCESSOR_HEADER_SET1(stype,param,vtype1) \
  braid_Int braid_##stype##StatusSet##param(braid_##stype##Status s, braid_##vtype1 v1);

/*----------------------------------------------------------------------------------
 * Define Status Structure. The `base' class is braid_Status, and all the other
 * status structures derive from this class.
 *----------------------------------------------------------------------------------*/
/** \defgroup braidstatusstruct XBraid status structures 
 *  \ingroup userinterface
 *  
 *  Define the different status types.
 *
 *  @{
 */

struct _braid_Status_struct;
/**
 * This is the main Status structure, that contains the properties of all the status.
 * The user does not have access to this structure, but only to the derived Status
 * structures. This class is accessed only inside XBraid code.
 */
typedef struct _braid_Status_struct *braid_Status;


/**
 * AccessStatus structure which defines the status of XBraid at a given instant
 * on some level during a run.  The user accesses it through
 * _braid_AccessStatusGet**()_ functions. This is just a pointer to the braid_Status.
 */
typedef struct _braid_AccessStatus_struct *braid_AccessStatus;

/**
 * The user's step routine routine will receive a StepStatus structure, which
 * defines the status of XBraid at the given instant for step evaluation on some level
 * during a run.  The user accesses it through _braid_StepStatusGet**()_ functions.
 * This is just a pointer to the braid_Status.
 */
typedef struct _braid_StepStatus_struct *braid_StepStatus;

/**
 * The user coarsen and refine routines will receive a CoarsenRefStatus structure, which
 * defines the status of XBraid at a given instant of coarsening or refinement on some level
 * during a run.  The user accesses it through _braid_CoarsenRefStatusGet**()_ functions.
 * This is just a pointer to the braid_Status.
 */
typedef struct _braid_CoarsenRefStatus_struct *braid_CoarsenRefStatus;

/**
 * The user's bufpack, bufunpack and bufsize routines will receive a BufferStatus structure, which
 * defines the status of XBraid at a given buff (un)pack instance.  The user accesses it
 * through _braid_BufferStatusGet**()_ functions. This is just a pointer to the braid_Status.
 */
typedef struct _braid_BufferStatus_struct *braid_BufferStatus;

/** @}*/


/*--------------------------------------------------------------------------
 * Routines for user to access XBraid status structures
 *--------------------------------------------------------------------------*/

/** \defgroup braidstatusroutines XBraid status routines
 *  \ingroup userinterface
 *  
 *  XBraid status structures and associated Get/Set routines are what tell 
 *  the user the status of the simulation  when their routines (step, 
 *  coarsen/refine, access) are called. 
 *
 *  @{
 */

/*--------------------------------------------------------------------------
 * Global Status Prototypes
 *--------------------------------------------------------------------------*/

/**
 * Return the current time from the Status structure.
 **/
braid_Int
braid_StatusGetT(braid_Status status,                      /**< structure containing current simulation info */
                 braid_Real  *t_ptr                        /**< output, current time */
                 );

/**
 * Return the index value corresponding to the current time value
 * from the Status structure.
 **/
braid_Int
braid_StatusGetTIndex(braid_Status status,                  /**< structure containing current simulation info */
                      braid_Int   *idx_ptr                  /**< output, global index value corresponding to current time value */
                      );

/**
 * Return the current iteration from the Status structure.
 **/
braid_Int
braid_StatusGetIter(braid_Status status,                   /**< structure containing current simulation info */
                    braid_Int   *iter_ptr                  /**< output, current XBraid iteration number*/
                    );

/**
 * Return the current XBraid level from the Status structure.
 **/
braid_Int
braid_StatusGetLevel(braid_Status status,                  /**< structure containing current simulation info */
                     braid_Int   *level_ptr                /**< output, current level in XBraid */
                     );

/**
 * Return the number of refinements done.
 **/
braid_Int
braid_StatusGetNRefine(braid_Status status,                /**< structure containing current simulation info */
                       braid_Int   *nrefine_ptr            /**< output, number of refinements done */
                       );

/**
 * Return the global number of time points on the fine grid.
 **/
braid_Int
braid_StatusGetNTPoints(braid_Status status,               /**< structure containing current simulation info */
                        braid_Int   *ntpoints_ptr          /**< output, number of time points on the fine grid */
                        );

/**
 * Return the current residual norm from the Status structure.
 **/
braid_Int
braid_StatusGetResidual(braid_Status status,               /**< structure containing current simulation info */
                        braid_Real  *rnorm_ptr             /**< output, current residual norm */
                        );

/**
 * Return whether XBraid is done for the current simulation.
 *
 * *done_ptr = 1* indicates that XBraid has finished iterating, 
 * (either maxiter has been reached, or the tolerance has been met).
 **/
braid_Int
braid_StatusGetDone(braid_Status status,                   /**< structure containing current simulation info */
                    braid_Int   *done_ptr                  /**< output,  =1 if XBraid has finished, else =0 */
                    );

/**
 * Return XBraid status for the current simulation. Four values are 
 * returned.
 *
 * TILD : time, iteration, level, done
 *
 * These values are also available through individual Get routines. 
 * These individual routines are the location of detailed documentation on 
 * each parameter, e.g., see *braid_StatusGetDone* for more information
 * on the *done* value.
 **/
braid_Int
braid_StatusGetTILD(braid_Status status,                   /**< structure containing current simulation info */
                    braid_Real  *t_ptr,                    /**< output, current time */
                    braid_Int   *iter_ptr,                 /**< output, current XBraid iteration number*/
                    braid_Int   *level_ptr,                /**< output, current level in XBraid */
                    braid_Int   *done_ptr                  /**< output,  =1 if XBraid has finished, else =0 */
                    );

/**
 * Return whether this is a wrapper test or an XBraid run
 **/
braid_Int
braid_StatusGetWrapperTest(braid_Status status,            /**< structure containing current simulation info */
                           braid_Int   *wtest_ptr          /**< output, =1 if this is a wrapper test, =0 if XBraid run */
                           );

/**
 * Return flag indicating from which function the vector is accessed
 **/
braid_Int
braid_StatusGetCallingFunction(braid_Status status,        /**< structure containing current simulation info */
                               braid_Int   *cfunction_ptr  /**< output, function number (0=FInterp, 1=FRestrict, 2=FRefine, 3=FAccess) */
                               );

/**
 * Return the **coarse grid** time value to the left of the current time value from
 * the Status structure.
 **/
braid_Int
braid_StatusGetCTprior(braid_Status status,                /**< structure containing current simulation info */
                       braid_Real  *ctprior_ptr            /**< output, time value to the left of current time value on coarse grid */
                       );

/**
 * Return the **coarse grid** time value to the right of the current time value from
 * the Status structure.
 **/
braid_Int
braid_StatusGetCTstop(braid_Status status,                 /**< structure containing current simulation info */
                      braid_Real  *ctstop_ptr              /**< output, time value to the right of current time value on coarse grid */
                      );

/**
 * Return the **fine grid** time value to the left of the current time value from
 * the Status structure.
 **/
braid_Int
braid_StatusGetFTprior(braid_Status status,                /**< structure containing current simulation info */
                       braid_Real  *ftprior_ptr            /**< output, time value to the left of current time value on fine grid */
                       );

/**
 * Return the **fine grid** time value to the right of the current time value from
 * the Status structure.
 **/
braid_Int
braid_StatusGetFTstop(braid_Status status,                 /**< structure containing current simulation info */
                      braid_Real  *ftstop_ptr              /**< output, time value to the right of current time value on fine grid */
                      );

/**
 * Return XBraid status for the current simulation. Five values are 
 * returned, tstart, f_tprior, f_tstop, c_tprior,  c_tstop. 
 *
 * These values are also available through individual Get routines. 
 * These individual routines are the location of detailed documentation on 
 * each parameter, e.g., see *braid_StatusGetCTprior* for more
 * information on the *c_tprior* value.
 **/
braid_Int
braid_StatusGetTpriorTstop(braid_Status status,            /**< structure containing current simulation info */
                           braid_Real  *t_ptr,             /**< output, current time */
                           braid_Real  *ftprior_ptr,       /**< output, time value to the left of current time value on fine grid */
                           braid_Real  *ftstop_ptr,        /**< output, time value to the right of current time value on fine grid */
                           braid_Real  *ctprior_ptr,       /**< output, time value to the left of current time value on coarse grid */
                           braid_Real  *ctstop_ptr         /**< output, time value to the right of current time value on coarse grid */
                           );

/**
 * Return the time value to the right of the current time value from
 * the Status structure.
 **/
braid_Int
braid_StatusGetTstop(braid_Status status,                  /**< structure containing current simulation info */
                     braid_Real  *tstop_ptr                /**< output, next time value to evolve towards */
                     );

/**
 * Return XBraid status for the current simulation. Two values are
 * returned, tstart and tstop.
 *
 * These values are also available through individual Get routines.
 * These individual routines are the location of detailed documentation on
 * each parameter, e.g., see *braid_StatusGetTstart* for more information
 * on the *tstart* value.
 **/
braid_Int
braid_StatusGetTstartTstop(braid_Status status,            /**< structure containing current simulation info */
                           braid_Real  *tstart_ptr,        /**< output, current time */
                           braid_Real  *tstop_ptr          /**< output, next time value to evolve towards */
                           );

/**
 * Return the current XBraid stopping tolerance
 **/
braid_Int
braid_StatusGetTol(braid_Status status,                    /**< structure containing current simulation info */
                   braid_Real  *tol_ptr                    /**< output, current XBraid stopping tolerance */
                   );

/**
 * Return the current XBraid residual history.  If *nrequest_ptr*
 * is negative, return the last *nrequest_ptr* residual norms.  If
 * positive, return the first *nrequest_ptr* residual norms.  Upon
 * exit, *nrequest_ptr* holds the number of residuals actually
 * returned.
 **/
braid_Int
braid_StatusGetRNorms(braid_Status status,                 /**< structure containing current simulation info */
                      braid_Int   *nrequest_ptr,           /**< input/output, input: number of requested residual norms, output: number actually copied */
                      braid_Real  *rnorms_ptr              /**< output, XBraid residual norm history, of length *nrequest_ptr* */
                      );

/**
 * Return the previous *old_fine_tolx* set through *braid_StatusSetOldFineTolx*
 * This is used especially by *braid_GetSpatialAccuracy
 **/
braid_Int
braid_StatusGetOldFineTolx(braid_Status status,            /**< structure containing current simulation info */
                           braid_Real  *old_fine_tolx_ptr  /**< output, previous *old_fine_tolx*, set through *braid_StepStatusSetOldFineTolx* */
                           );

/**
 * Set *old_fine_tolx*, available for retrieval through *braid_StatusGetOldFineTolx*
 * This is used especially by *braid_GetSpatialAccuracy
 **/
braid_Int
braid_StatusSetOldFineTolx(braid_Status status,            /**< structure containing current simulation info */
                           braid_Real   old_fine_tolx      /**< input, the last used fine_tolx */
                           );

/**
 * Set *tight_fine_tolx*, boolean variable indicating whether the tightest
 * tolerance has been used for spatial solves (implicit schemes).  This value
 * must be 1 in order for XBraid to halt (unless maxiter is reached)
 **/
braid_Int
braid_StatusSetTightFineTolx(braid_Status status,          /**< structure containing current simulation info */
                             braid_Real   tight_fine_tolx  /**< input, boolean indicating whether the tight tolx has been used */
                             );

/**
 * Set the rfactor, a desired refinement factor for this interval.  rfactor=1
 * indicates no refinement, otherwise, this inteval is subdivided rfactor
 * times. 
 **/
braid_Int
braid_StatusSetRFactor(braid_Status status,                /**< structure containing current simulation info */
                       braid_Real   rfactor                /**< input, user-determined desired rfactor */
                       );

/**
 * Set the r_space flag. When set = 1, spatial coarsening will be called,
 * for all local time points, following the  completion of the current
 * iteration, provided rfactors are not set at any global time point. This
 * allows for spatial refinment without temporal refinment
 **/
braid_Int
braid_StatusSetRSpace(braid_Status status,                 /**< structure containing current simulation info */
                      braid_Real   r_space                 /**< input, if 1, call spatial refinement on finest grid after this iter */
                      );

/**
 * Return the current message type from the Status structure.
 **/
braid_Int
braid_StatusGetMessageType(braid_Status status,            /**< structure containing current simulation info */
                           braid_Int   *messagetype_ptr    /**< output, type of message, 0: for Step(), 1: for load balancing */
                           );

/**
 * Set the size of the buffer. If set by user, the send buffer will
   be "size" bytes in length. If not, BufSize is used.
 **/
braid_Int
braid_StatusSetSize(braid_Status status,                   /**< structure containing current simulation info */
                    braid_Real   size                      /**< input, size of the send buffer */
                    );
/** @}*/


/*--------------------------------------------------------------------------
 * Begin definition of `inherited' Get/Set functions from the base class, 
 * use macros to accomplish this
 *--------------------------------------------------------------------------*/
/** \defgroup braidstatusroutinesinherited Inherited XBraid status routines
 *  \ingroup userinterface
 *  
 *  These are the `inherited' Status Get/Set functions.  See the 
 *  *XBraid status routines* section for the description of each function.
 *  For example, for braid_StepStatusGetT(...), you would look up 
 *  braid_StatusGetT(...)
 *
 *  @{
 */

/*--------------------------------------------------------------------------
 * AccessStatus Prototypes: They just wrap the corresponding Status accessors
 *--------------------------------------------------------------------------*/

ACCESSOR_HEADER_GET1(Access, T,               Real)
ACCESSOR_HEADER_GET1(Access, TIndex,          Int)
ACCESSOR_HEADER_GET1(Access, Iter,            Int)
ACCESSOR_HEADER_GET1(Access, Level,           Int)
ACCESSOR_HEADER_GET1(Access, NRefine,         Int)
ACCESSOR_HEADER_GET1(Access, NTPoints,        Int)
ACCESSOR_HEADER_GET1(Access, Residual,        Real)
ACCESSOR_HEADER_GET1(Access, Done,            Int)
ACCESSOR_HEADER_GET4(Access, TILD,            Real, Int, Int, Int)
ACCESSOR_HEADER_GET1(Access, WrapperTest,     Int)
ACCESSOR_HEADER_GET1(Access, CallingFunction, Int)

/*--------------------------------------------------------------------------
 * CoarsenRefStatus Prototypes: They just wrap the corresponding Status accessors
 *--------------------------------------------------------------------------*/

ACCESSOR_HEADER_GET1(CoarsenRef, T,           Real)
ACCESSOR_HEADER_GET1(CoarsenRef, TIndex,      Int)
ACCESSOR_HEADER_GET1(CoarsenRef, Iter,        Int)
ACCESSOR_HEADER_GET1(CoarsenRef, Level,       Int)
ACCESSOR_HEADER_GET1(CoarsenRef, NRefine,     Int)
ACCESSOR_HEADER_GET1(CoarsenRef, NTPoints,    Int)
ACCESSOR_HEADER_GET1(CoarsenRef, CTprior,     Real)
ACCESSOR_HEADER_GET1(CoarsenRef, CTstop,      Real)
ACCESSOR_HEADER_GET1(CoarsenRef, FTprior,     Real)
ACCESSOR_HEADER_GET1(CoarsenRef, FTstop,      Real)
ACCESSOR_HEADER_GET5(CoarsenRef, TpriorTstop, Real, Real, Real, Real, Real)

/*--------------------------------------------------------------------------
 * StepStatus Prototypes: They just wrap the corresponding Status accessors
 *--------------------------------------------------------------------------*/

ACCESSOR_HEADER_GET1(Step, T,             Real)
ACCESSOR_HEADER_GET1(Step, TIndex,        Int)
ACCESSOR_HEADER_GET1(Step, Iter,          Int)
ACCESSOR_HEADER_GET1(Step, Level,         Int)
ACCESSOR_HEADER_GET1(Step, NRefine,       Int)
ACCESSOR_HEADER_GET1(Step, NTPoints,      Int)
ACCESSOR_HEADER_GET1(Step, Tstop,         Real)
ACCESSOR_HEADER_GET2(Step, TstartTstop,   Real, Real)
ACCESSOR_HEADER_GET1(Step, Tol,           Real)
ACCESSOR_HEADER_GET2(Step, RNorms,        Int,  Real)
ACCESSOR_HEADER_GET1(Step, OldFineTolx,   Real)
ACCESSOR_HEADER_SET1(Step, OldFineTolx,   Real)
ACCESSOR_HEADER_SET1(Step, TightFineTolx, Real)
ACCESSOR_HEADER_SET1(Step, RFactor,       Real)
ACCESSOR_HEADER_SET1(Step, RSpace,        Real)

/*--------------------------------------------------------------------------
 * BufferStatus Prototypes: They just wrap the corresponding Status accessors
 *--------------------------------------------------------------------------*/

ACCESSOR_HEADER_GET1(Buffer, MessageType, Int)
ACCESSOR_HEADER_SET1(Buffer, Size,        Real)

/** @}*/


/*--------------------------------------------------------------------------
 * Macros 
 *--------------------------------------------------------------------------*/
/** \defgroup braidstatusmacros XBraid status macros
 *  \ingroup userinterface
 * Macros defining Status values that the user can obtain during runtime, which will
 * tell the user where in Braid the current cycle is, e.g. in the FInterp function.
 *  @{
 */

/** When CallingFunction equals 0, Braid is in FInterp */
#define braid_ASCaller_FInterp   0
/** When CallingFunction equals 0, Braid is in FRestrict */
#define braid_ASCaller_FRestrict 1
/** When CallingFunction equals 0, Braid is in FRefine */
#define braid_ASCaller_FRefine   2
/** When CallingFunction equals 0, Braid is in FAccess */
#define braid_ASCaller_FAccess   3

/** @}*/


#ifdef __cplusplus
}
#endif

#endif



