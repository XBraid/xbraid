/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory.
 * 
 * This file is part of XBraid. For support, post issues to the XBraid Github page.
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
 * \brief Define headers for the user-interface with the XBraid status structures, 
 * allowing the user to get/set status structure values. 
 *
 */

#ifndef braid_status_HEADER
#define braid_status_HEADER

#include "_braid.h"
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
#define ACCESSOR_HEADER_GET1_IN1(stype,param,vtype1,vtype2) \
   braid_Int braid_##stype##StatusGet##param(braid_##stype##Status s, braid_##vtype1 *v1, braid_##vtype2 v2);
#define ACCESSOR_HEADER_GET1_IN2(stype,param,vtype1,vtype2,vtype3) \
   braid_Int braid_##stype##StatusGet##param(braid_##stype##Status s, braid_##vtype1 *v1, braid_##vtype2 v2, braid_##vtype3 v3);
#define ACCESSOR_HEADER_GET1_IN3(stype,param,vtype1,vtype2,vtype3,vtype4) \
   braid_Int braid_##stype##StatusGet##param(braid_##stype##Status s, braid_##vtype1 *v1, braid_##vtype2 v2, braid_##vtype3 v3, braid_##vtype4 v4);
#define ACCESSOR_HEADER_GET2(stype,param,vtype1,vtype2) \
  braid_Int braid_##stype##StatusGet##param(braid_##stype##Status s, braid_##vtype1 *v1, braid_##vtype2 *v2);
#define ACCESSOR_HEADER_GET2_IN1(stype,param,vtype1,vtype2,vtype3) \
   braid_Int braid_##stype##StatusGet##param(braid_##stype##Status s, braid_##vtype1 *v1, braid_##vtype2 *v2, braid_##vtype3 v3);
#define ACCESSOR_HEADER_GET3(stype,param,vtype1,vtype2,vtype3) \
  braid_Int braid_##stype##StatusGet##param(braid_##stype##Status s, braid_##vtype1 *v1, braid_##vtype2 *v2, braid_##vtype3 *v3);
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
 * SyncStatus structure which provides the status of XBraid at a given instant
 * on some level during a run. This is vector independent and called once per
 * processor. The user accesses it through _braid_SyncStatusGet**()_ functions.
 * This is just a pointer to the braid_Status.
 */
typedef struct _braid_SyncStatus_struct *braid_SyncStatus;

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

/**
 * The user's objectiveT and PostprocessObjective will receive an ObjectiveStatus structure, which
 * defines the status of XBraid at a given instance of evaluating the objective function.  The user accesses it
 * through _braid_ObjectiveStatusGet**()_ functions. This is just a pointer to the braid_Status.
 */
typedef struct _braid_ObjectiveStatus_struct *braid_ObjectiveStatus;


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

/** Return the index value corresponding to the current time value from the
 * Status structure.  
 *
 * For Step(), this corresponds to the time-index of "tstart", as this is the
 * time-index of the input vector.  That is, NOT the time-index of "tstop".
 * For Access, this corresponds just simply to the time-index of the input
 * vector.
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
 * Return the total number of XBraid levels from the Status structure.
 **/
braid_Int
braid_StatusGetNLevels(braid_Status status,                  /**< structure containing current simulation info */
                       braid_Int   *nlevels_ptr              /**< output, number of levels in XBraid */
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
 * Returns upper and lower time point indices on this processor. Two
 * values are returned. Requires the user to specify which level they
 * want the time point indices from.
 **/
braid_Int
braid_StatusGetTIUL(braid_Status status,                   /**< structure containing current simulation info */
                    braid_Int   *iloc_upper,               /**< output, the upper time point index on this processor */
                    braid_Int   *iloc_lower,               /**< output, the lower time point index on this processor */
                    braid_Int    level                     /**< input, level for the desired indices */
                    );

/**
 * Returns an array of time values corresponding to the given inputs.
 * The inputs are the level you want the time values from, the upper
 * time point index you want the value of, and the lower time point
 * index you want the time value of. The output is then filled with
 * all time values from the upper index to the lower index, inclusive.
 *
 * The caller is responsible for allocating and managing the memory
 * for the array. Time values are filled in so that tvalues_ptr[0]
 * corresponds to the lower time index.
 **/
braid_Int
braid_StatusGetTimeValues(braid_Status status,             /**< structure containing current simulation info */
                          braid_Real **tvalues_ptr,        /**< output, time point values for the requested range of indices */
                          braid_Int    i_upper,            /**< input, upper index of the desired time value range (inclusive) */
                          braid_Int    i_lower,            /**< input, lower index of the desired time value range (inclusive) */
                          braid_Int    level               /**< input, level for the desired time values */
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
                               braid_Int   *cfunction_ptr  /**< output, function number (0=FInterp, 1=FRestrict, 2=FRefine, 3=FAccess, 4=FRefine after refinement, 5=Drive Top of Cycle) */
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
 * Returns the processor number in *proc_ptr* on which the time step *index*
 * lives for the given *level*.  Returns -1 if *index* is out of range.
 * This is used especially by the _braid_SyncStatus functionality
 **/
braid_Int
braid_StatusGetProc(braid_Status  status,                  /**< structure containing current simulation info */
                    braid_Int    *proc_ptr,                /**< output, the processor number corresponding to the level and time point index inputs */
                    braid_Int     level,                   /**< input, level for the desired processor */
                    braid_Int     index                    /**< input, the global time point index for the desired processor */
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
 * times (uniform refinement).
 **/
braid_Int
braid_StatusSetRFactor(braid_Status status,                /**< structure containing current simulation info */
                       braid_Real   rfactor                /**< input, user-determined desired rfactor */
                       );


/**
 * Set time step sizes for refining the time interval non-uniformly.
 **/
braid_Int
braid_StatusSetRefinementDtValues(braid_Status status,      /**< structure containing current simulation info */
                                  braid_Real   rfactor,     /**< input, number of subintervals */
                                  braid_Real*   dtarray     /**< input, array of dt values for non-uniform refinement */
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

/**
 * Set the size of the buffer for basis vectors. 
 * If set by user, the send buffer will
 * be "size" bytes in length. If not, BufSize is used.
 **/
braid_Int
braid_StatusSetBasisSize(braid_Status status,                   /**< structure containing current simulation info */
                         braid_Real   size                      /**< input, size of the send buffer */
                         );

/** 
 * Get the Richardson based error estimate at the single time point currently
 * being "Stepped", i.e., return the current error estimate for the time point
 * at "tstart".
 *
 * Note that Step needs specific logic distinct from Access, hence please use
 * [braid_AccessStatusGetSingleErrorEstAccess](@ref braid_AccessStatusGetSingleErrorEstAccess)
 * for the user Access() function.
 */

braid_Int
braid_StatusGetSingleErrorEstStep(braid_Status   status,           /**< structure containing current simulation info */
                                  braid_Real    *estimate          /**< output, error estimate, equals -1 if not available yet (e.g., before iteration 1, or after refinement) */
                                  );

/** 
 * Get the Richardson based error estimate at the single time point currently 
 * accessible from Access.
 * 
 * Note that Access needs specific logic distinct from Step, hence please use
 * [braid_StepStatusGetSingleErrorEstStep](@ref braid_StepStatusGetSingleErrorEstStep)
 * for the user Step() function. 
 */
braid_Int
braid_StatusGetSingleErrorEstAccess(braid_Status   status,           /**< structure containing current simulation info */
                                    braid_Real    *estimate          /**< output, error estimate, equals -1 if not available yet (e.g., before iteration 1, or after refinement) */
                                    );



/** 
 * Get the number of local Richardson-based error estimates stored on this
 * processor.  Use this function in conjuction with GetAllErrorEst().
 * Workflow: use this function to get the size of the needed user-array that
 * will hold the error estimates, then pre-allocate array, then call
 * GetAllErrorEst() to write error estimates to the user-array, then
 * post-process array in user-code.  This post-processing will often occur in
 * the Sync function.  See examples/ex-06.c.
 */
braid_Int
braid_StatusGetNumErrorEst(braid_Status   status,        /**< structure containing current simulation info */
                           braid_Int     *npoints        /**< output, number of locally stored Richardson error estimates */
                           );

/** Get All the Richardson based error estimates, e.g. from inside Sync.  Use
 * this function in conjuction with GetNumErrorEst().  Workflow: use
 * GetNumErrorEst() to get the size of the needed user-array that will hold the
 * error estimates, then pre-allocate array, then call this function to write
 * error estimates to the user-array, then post-process array in user-code.
 * This post-processing will often occur in the Sync function.  See
 * examples/ex-06.c.
 *
 * The error_est array must be user-allocated.
 */
braid_Int
braid_StatusGetAllErrorEst(braid_Status    status,       /**< structure containing current simulation info */
                           braid_Real     *error_est     /**< output, user-allocated error estimate array, written by Braid, equals -1 if not available yet (e.g., before iteration 1, or after refinement) */
                           );

/**
 * Gets accces to the temporal communicator. Allows this processor
 * to access other temporal processors.
 * This is used especially by Sync.
 **/
braid_Int
braid_StatusGetTComm(braid_Status  status,            /**< structure containing current simulation info */
                     MPI_Comm     *comm_ptr           /**< output, temporal communicator */
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
ACCESSOR_HEADER_GET1(Access, NLevels,         Int)
ACCESSOR_HEADER_GET1(Access, NRefine,         Int)
ACCESSOR_HEADER_GET1(Access, NTPoints,        Int)
ACCESSOR_HEADER_GET1(Access, Residual,        Real)
ACCESSOR_HEADER_GET1(Access, Done,            Int)
ACCESSOR_HEADER_GET4(Access, TILD,            Real, Int, Int, Int)
ACCESSOR_HEADER_GET1(Access, WrapperTest,     Int)
ACCESSOR_HEADER_GET1(Access, CallingFunction, Int)
ACCESSOR_HEADER_GET1(Access, SingleErrorEstAccess, Real)
ACCESSOR_HEADER_GET1(Access, DeltaRank, Int)
ACCESSOR_HEADER_GET1_IN1(Access, BasisVec, Vector, Int)

/*--------------------------------------------------------------------------
 * SyncStatus Prototypes: They just wrap the corresponding Status accessors
 *--------------------------------------------------------------------------*/

ACCESSOR_HEADER_GET2_IN1(Sync, TIUL,         Int, Int, Int)
ACCESSOR_HEADER_GET1_IN3(Sync, TimeValues,   Real*, Int, Int, Int)
ACCESSOR_HEADER_GET1_IN2(Sync, Proc,         Int, Int, Int)
ACCESSOR_HEADER_GET1(Sync, Iter,             Int)
ACCESSOR_HEADER_GET1(Sync, Level,            Int)
ACCESSOR_HEADER_GET1(Sync, NLevels,          Int)
ACCESSOR_HEADER_GET1(Sync, NRefine,          Int)
ACCESSOR_HEADER_GET1(Sync, NTPoints,         Int)
ACCESSOR_HEADER_GET1(Sync, Done,             Int)
ACCESSOR_HEADER_GET1(Sync, CallingFunction,  Int)
ACCESSOR_HEADER_GET1(Sync, NumErrorEst,      Int)
ACCESSOR_HEADER_GET1(Sync, AllErrorEst,      Real)
ACCESSOR_HEADER_GET1(Sync, TComm,            MPI_Comm)

/*--------------------------------------------------------------------------
 * CoarsenRefStatus Prototypes: They just wrap the corresponding Status accessors
 *--------------------------------------------------------------------------*/

ACCESSOR_HEADER_GET1(CoarsenRef, T,           Real)
ACCESSOR_HEADER_GET1(CoarsenRef, TIndex,      Int)
ACCESSOR_HEADER_GET1(CoarsenRef, Iter,        Int)
ACCESSOR_HEADER_GET1(CoarsenRef, Level,       Int)
ACCESSOR_HEADER_GET1(CoarsenRef, NLevels,     Int)
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

ACCESSOR_HEADER_GET2_IN1(Step, TIUL,      Int, Int, Int)
ACCESSOR_HEADER_GET1(Step, T,             Real)
ACCESSOR_HEADER_GET1(Step, TIndex,        Int)
ACCESSOR_HEADER_GET1(Step, Iter,          Int)
ACCESSOR_HEADER_GET1(Step, Level,         Int)
ACCESSOR_HEADER_GET1(Step, NLevels,       Int)
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
ACCESSOR_HEADER_GET1(Step, Done,          Int)
ACCESSOR_HEADER_GET1(Step, SingleErrorEstStep, Real)
ACCESSOR_HEADER_GET1(Step, CallingFunction,    Int)
ACCESSOR_HEADER_GET1(Step, DeltaRank, Int)
ACCESSOR_HEADER_GET1_IN1(Step, BasisVec, Vector, Int)

/*--------------------------------------------------------------------------
 * BufferStatus Prototypes: They just wrap the corresponding Status accessors
 *--------------------------------------------------------------------------*/

ACCESSOR_HEADER_GET1(Buffer, MessageType, Int)
ACCESSOR_HEADER_SET1(Buffer, Size,        Real)
ACCESSOR_HEADER_SET1(Buffer, BasisSize,   Real)

/*--------------------------------------------------------------------------
 * ObjectiveStatus Prototypes: They just wrap the corresponding Status accessors
 *--------------------------------------------------------------------------*/

ACCESSOR_HEADER_GET1(Objective, T,             Real)
ACCESSOR_HEADER_GET1(Objective, TIndex,        Int)
ACCESSOR_HEADER_GET1(Objective, Iter,          Int)
ACCESSOR_HEADER_GET1(Objective, Level,         Int)
ACCESSOR_HEADER_GET1(Objective, NLevels,       Int)
ACCESSOR_HEADER_GET1(Objective, NRefine,       Int)
ACCESSOR_HEADER_GET1(Objective, NTPoints,      Int)
ACCESSOR_HEADER_GET1(Objective, Tol,           Real)


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
#define braid_ASCaller_FInterp 0
/** When CallingFunction equals 1, Braid is in FRestrict */
#define braid_ASCaller_FRestrict 1
/** When CallingFunction equals 2, Braid is in FRefine */
#define braid_ASCaller_FRefine 2
/** When CallingFunction equals 3, Braid is in FAccess */
#define braid_ASCaller_FAccess 3
/** When CallingFunction equals 4, Braid is inside FRefine after the new finest
 * level has been initialized */
#define braid_ASCaller_FRefine_AfterInitHier 4
/** When CallingFunction equals 5, Braid is at the top of the cycle */
#define braid_ASCaller_Drive_TopCycle 5
/** When CallingFunction equals 6, Braid is in FCrelax */
#define braid_ASCaller_FCRelax 6
/** When CallingFunction equals 7, Braid just finished initialization */
#define braid_ASCaller_Drive_AfterInit 7
/** When CallingFunction equals 8, Braid is in BaseStep_diff */
#define braid_ASCaller_BaseStep_diff 8
/** When CallingFunction equals 9, Braid is in ComputeFullRNorm */
#define braid_ASCaller_ComputeFullRNorm 9
/** When CallingFunction equals 10, Braid is in FASResidual */
#define braid_ASCaller_FASResidual 10
/** When CallingFunction equals 11, Braid is in Residual, immediately after restriction */
#define braid_ASCaller_Residual 11
/** When CallingFunction equals 12, Braid is in InitGuess */
#define braid_ASCaller_InitGuess 12

/** @}*/


#ifdef __cplusplus
}
#endif

#endif



