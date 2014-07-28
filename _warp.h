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

/** \file _warp.h
 * \brief Define headers for developer routines.
 *
 * This file contains the headers for developer routines.
 */

#ifndef _warp_HEADER
#define _warp_HEADER

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <math.h>

#include "warp.h"

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------------
 * Main data structures and accessor macros
 *--------------------------------------------------------------------------*/

/**
 * Points to the status structure which defines the status of Warp 
 * at a given instant on a some level during a run.  The user accesses
 * it through warp_Get**Status() functions.
 **/
typedef struct _warp_Status_struct
{
   warp_Int     iter;         /**< warp iteration number */
   warp_Int     level;        /**< current level in warp */
   warp_Real    rnorm;        /**< residual norm */
   warp_Int     done;         /**< boolean describing whether warp has finished */
   
} _warp_Status;

/**
 * Warp comm handle structure\n
 * Used for initiating and completing nonblocking communication to pass
 * warp_Vectors between processors.
 **/
typedef struct
{
   warp_Int     request_type;    /**< two values: recv type = 1, and send type = 0 */
   warp_Int     num_requests;    /**< number of active requests for this handle, usually 1 */
   MPI_Request *requests;        /**< MPI request structure */
   MPI_Status  *status;          /**< MPI status */
   void        *buffer;          /**< Buffer for message */
   warp_Vector *vector_ptr;      /**< warp_vector being sent/received */
   
} _warp_CommHandle;

/**
 * warp Accuracy Handle, used for controlling the accuracy of solves during
 * implicit time stepping.  For instance, to do less accurate solves on coarse
 * time grids
 **/
typedef struct
{
   warp_Int   matchF;
   warp_Real  value;      /**< accuracy value */
   warp_Real  old_value;  /**< old accuracy value used in FRestrict */
   warp_Real  loose;      /**< loose accuracy for spatial solves */
   warp_Real  tight;      /**< tight accuracy for spatial solves */
   warp_Int   tight_used; /**< tight accuracy used (1) or not (0) */
} _warp_AccuracyHandle;

/**
 * Warp Grid structure for a certain time level\n
 * Holds all the information for a processor related to the temporal
 * grid at this level.
 **/
typedef struct
{
   warp_Int     level;              /**< Level that grid is on */
   warp_Int     ilower, iupper;     /**< smallest and largest time indices at this level*/
   warp_Int     clower, cupper;     /**< smallest C point index, largest C point index */
   warp_Int     cfactor, ncpoints;  /**< coarsening factor and number of C points */

   warp_Vector *ua;                 /**< unknown vectors            (C-points only)*/
   warp_Real   *ta;                 /**< time values                (all points) */
   warp_Vector *va;                 /**< restricted unknown vectors (all points, NULL on level 0) */
   warp_Vector *wa;                 /**< rhs vectors f-v            (all points, NULL on level 0) */

   warp_Int          recv_index;    /**<  -1 means no receive */
   warp_Int          send_index;    /**<  -1 means no send */
   _warp_CommHandle *recv_handle;   /**<  Handle for nonblocking receives of warp_Vectors */
   _warp_CommHandle *send_handle;   /**<  Handle for nonblocking sends of warp_Vectors */

   /* pointers to the original memory allocation for ua, ta, va, and wa */
   warp_Vector *ua_alloc;
   warp_Real   *ta_alloc;
   warp_Vector *va_alloc;
   warp_Vector *wa_alloc;

} _warp_Grid;

/**
 * The typedef _warp_Core struct is a **critical** part of warp and 
 * is passed to *each* routine in warp.  It thus allows each routine access 
 * to Warp attributes.
 **/
typedef struct _warp_Core_struct
{
   MPI_Comm              comm_world;
   MPI_Comm              comm;         /**< communicator for the time dimension */
   warp_Real             tstart;       /**< start time */
   warp_Real             tstop;        /**< stop time */
   warp_Int              ntime;        /**< initial number of time intervals */
   warp_App              app;          /**< application data for the user */
   
   warp_PtFcnPhi         phi;          /**< apply phi function */
   warp_PtFcnInit        init;         /**< return an initial solution vector */
   warp_PtFcnClone       clone;        /**< clone a vector */
   warp_PtFcnFree        free;         /**< free up a vector */
   warp_PtFcnSum         sum;          /**< vector sum */
   warp_PtFcnDot         dot;          /**< dot product */
   warp_PtFcnWrite       write;        /**< write the vector */
   warp_PtFcnBufSize     bufsize;      /**< return buffer size */
   warp_PtFcnBufPack     bufpack;      /**< pack a buffer */
   warp_PtFcnBufUnpack   bufunpack;    /**< unpack a buffer */
   warp_PtFcnCoarsen     coarsen;      /**< (optional) return a coarsened vector */
   warp_PtFcnRefine      refine;       /**< (optional) return a refined vector */

   warp_Int              write_level;  /**< determines how often to call the user's write routine */ 
   warp_Int              print_level;  /**< determines amount of output printed to screem (0,1,2) */ 
   warp_Int              max_levels;   /**< maximum number of temporal grid levels */
   warp_Int              max_coarse;   /**< maximum allowed coarse grid size  (in terms of C-points) */
   warp_Real             tol;          /**< stopping tolerance */
   warp_Int              rtol;         /**< use relative tolerance */
   warp_Int             *nrels;        /**< number of pre-relaxations on each level */
   warp_Int              nrdefault;    /**< default number of pre-relaxations */
   warp_Int             *cfactors;     /**< coarsening factors */
   warp_Int              cfdefault;    /**< default coarsening factor */
   warp_Int              max_iter;     /**< maximum number of multigrid in time iterations */
   warp_Int              niter;        /**< number of iterations */
   warp_Real             rnorm;        /**< residual norm */
   warp_Int              fmg;          /**< use FMG cycle */
   warp_Int              nfmg_Vcyc;    /**< number of V-cycle calls at each level in FMG */
   _warp_AccuracyHandle *accuracy;     /**< accuracy of spatial solves on different levels */

   warp_Int              gupper;       /**< global upper index on the fine grid */

   warp_Int             *rfactors;     /**< refinement factors for finest grid (if any) */

   warp_Int              nlevels;      /**< number of temporal grid levels */
   _warp_Grid          **grids;        /**< pointer to temporal grid structures for each level*/

   warp_Real             localtime;    /**< local wall time for warp_Drive() */
   warp_Real             globaltime;   /**< global wall time for warp_Drive() */

} _warp_Core;

/* Accessor macros */

#define _warp_CommHandleElt(handle, elt)  ((handle) -> elt)

#define _warp_GridElt(grid, elt)  ((grid) -> elt)

#define _warp_StatusElt(status, elt) ( (status) -> elt )
#define _warp_CoreElt(core, elt)     (  (core)  -> elt )
#define _warp_CoreFcn(core, fcn)     (*((core)  -> fcn))

/*--------------------------------------------------------------------------
 * Memory allocation macros
 *--------------------------------------------------------------------------*/

#define _warp_TAlloc(type, count) \
( (type *)malloc((size_t)(sizeof(type) * (count))) )

#define _warp_CTAlloc(type, count) \
( (type *)calloc((size_t)(count), (size_t)sizeof(type)) )

#define _warp_TReAlloc(ptr, type, count) \
( (type *)realloc((char *)ptr, (size_t)(sizeof(type) * (count))) )

#define _warp_TFree(ptr) \
( free((char *)ptr), ptr = NULL )

/*--------------------------------------------------------------------------
 * Error handling
 *--------------------------------------------------------------------------*/

extern warp_Int _warp_error_flag;

/*--------------------------------------------------------------------------
 * Print file for redirecting stdout when needed
 *--------------------------------------------------------------------------*/

extern FILE *_warp_printfile;

/*--------------------------------------------------------------------------
 * Coarsening macros
 *--------------------------------------------------------------------------*/

#define _warp_MapFineToCoarse(findex, cfactor, cindex) \
( cindex = (findex)/(cfactor) )
#define _warp_MapCoarseToFine(cindex, cfactor, findex) \
( findex = (cindex)*(cfactor) )
#define _warp_IsFPoint(index, cfactor) \
( (index)%(cfactor) )
#define _warp_IsCPoint(index, cfactor) \
( !_warp_IsFPoint(index, cfactor) )

/*--------------------------------------------------------------------------
 * Prototypes
 *--------------------------------------------------------------------------*/

/**
 * Determine processor distribution.  This must agree with GetProc().
 * For the processor rank calling this function, it returns the smallest
 * and largest time indices ( *ilower_ptr* and *iupper_ptr*) that belong to 
 * that processor (the indices may * be F or C points).
 */
warp_Int
_warp_GetDistribution(warp_Core   core,
                      warp_Int   *ilower_ptr,
                      warp_Int   *iupper_ptr);

/**
 * Return the processor number in *proc_ptr* on which the time step *index* 
 * lives for the given *level*.  * Returns -1 if *index* is out of range
 */
warp_Int
_warp_GetProc(warp_Core   core,
              warp_Int    level,
              warp_Int    index,
              warp_Int   *proc_ptr);

/**
 * Initialize a receive to go into *vector_ptr* for the given time *index* on *level*.  
 * Also return a comm handle *handle_ptr* for querying later, to see if the receive has 
 * occurred.
 */
warp_Int
_warp_CommRecvInit(warp_Core           core,
                   warp_Int            level,
                   warp_Int            index,
                   warp_Vector        *vector_ptr,
                   _warp_CommHandle  **handle_ptr);

/**
 * Initialize a send of *vector* for the given time *index* on *level*.  
 * Also return a comm handle *handle_ptr* for querying later, to see if the 
 * send has occurred.
 */
warp_Int
_warp_CommSendInit(warp_Core           core,
                   warp_Int            level,
                   warp_Int            index,
                   warp_Vector         vector,
                   _warp_CommHandle  **handle_ptr);

/**
 * Block on the comm handle *handle_ptr* until the MPI operation (send or recv)
 * has completed
 */
warp_Int
_warp_CommWait(warp_Core          core,
               _warp_CommHandle **handle_ptr);

/**
 * Working on all intervals\n
 * At *level*, post a receive for the point to the left of ilower (regardless 
 * whether ilower is F or C).  Then, post a send of iupper if iupper is a C 
 * point.
 */
warp_Int
_warp_UCommInit(warp_Core  core,
                warp_Int   level);

/**
 * Working only on F-pt intervals\n
 * At *level*, **only** post a receive for the point to the left of ilower 
 * if ilower is an F point.  Then, post a send of iupper if iupper is a C point.
 */
warp_Int
_warp_UCommInitF(warp_Core  core,
                 warp_Int   level);

/**
 * Finish up communication\n
 * On *level*, wait on both the recv and send handles at this level.
 */
warp_Int
_warp_UCommWait(warp_Core  core,
                warp_Int   level);

/**
 * Retrieve the time step indices at this *level* which correspond to the FC interval
 * given by *interval_index*.  *ci_ptr* is the time step index for the C point
 * and *flo_ptr* and *fhi_ptr* are the smallest and largest F point indices in this
 * interval.  *flo* = *ci* +1, and *fhi* = *ci* + coarsening_factor - 1
 */
warp_Int
_warp_UGetInterval(warp_Core   core,
                   warp_Int    level,
                   warp_Int    interval_index,
                   warp_Int   *flo_ptr,
                   warp_Int   *fhi_ptr,
                   warp_Int   *ci_ptr);

/**
 * Returns a reference to the local u-vector in *u_ptr* for the grid *level* at 
 * point *index*.  Caveat: if *index* is not a C-point and within my index range, 
 * NULL is returned.
 */
warp_Int
_warp_UGetVectorRef(warp_Core     core,
                    warp_Int      level,
                    warp_Int      index,
                    warp_Vector  *u_ptr);

/**
 * Stores a reference to the vector *u* on grid *level* at point *index*.
 * If *index* is not a C-point and within this processor's range of time points, 
 * then nothing is done.
 */
warp_Int
_warp_USetVectorRef(warp_Core    core,
                    warp_Int     level,
                    warp_Int     index,
                    warp_Vector  u);

/**
 * Returns the u-vector in *u_ptr* on grid *level* at point *index*.  If *index* is my
 * "receive index" (as set by UCommInit(), for example), the u-vector  will be
 * received from a neighbor processor.  If *index* is within my index range and
 * is also a C-point, the saved value of u will be used.  A NULL value is
 * returned otherwise.
 */
warp_Int
_warp_UGetVector(warp_Core     core,
                 warp_Int      level,
                 warp_Int      index,
                 warp_Vector  *u_ptr);

/**
 * Sets the u-vector on grid *level* at point *index*.  If *index* is my "send
 * index" (as set by UCommInit(), for example), a send is initiated to a neighbor 
 * processor.  If *index* is within my index range and is also a C-point, the 
 * value is saved locally.
 */
warp_Int
_warp_USetVector(warp_Core    core,
                 warp_Int     level,
                 warp_Int     index,
                 warp_Vector  u);

/**
 * Call the user's write function to write *u* which is the vector corresponding to 
 * time step *index* on *level*.  *status* holds state information about the current 
 * Warp iteration, time value, etc...
 */
warp_Int
_warp_UWriteVector(warp_Core    core,
                   warp_Int     level,
                   warp_Int     index,
                   warp_Status  status,
                   warp_Vector  u);

/**
 * Apply Phi to the vector *u*\n
 * This is the vector corresponding to the time step *index* on *level*.
 * *accuracy* is a user set variable to allow for tuning the accuracy of 
 * spatial solvesfor implicit stepping. And, *rfactor* allows the user to
 * subdivide time intervals for accuracy purposes.
 */
warp_Int
_warp_Phi(warp_Core     core,
          warp_Int      level,
          warp_Int      index,
          warp_Real     accuracy,
          warp_Vector   u,
          warp_Int     *rfactor);


/**
 * Integrate one time step at time step *index* to time step *index*+1\n
 */
warp_Int
_warp_Step(warp_Core     core,
           warp_Int      level,
           warp_Int      index,
           warp_Real     accuracy,
           warp_Vector   u);

/**
 * Coarsen in space on *level* by calling the user's coarsen function.
 * The vector corresponding to the time step index *f_index* on the fine 
 * grid is coarsened to the time step index *c_index* on the coarse grid.
 * The output goes in *cvector* and the input vector is *fvector*.
 */
warp_Int
_warp_Coarsen(warp_Core     core,
              warp_Int      level,    /* coarse level */
              warp_Int      f_index,  /* fine index */
              warp_Int      c_index,  /* coarse index */
              warp_Vector   fvector,
              warp_Vector  *cvector);

/**
 * Refine in space on *level* by calling the user's refine function.
 * The vector corresponding to the time step index *c_index* on the coarse
 * grid is refined to the time step index *f_index* on the fine grid.
 * The output goes in *fvector* and the input vector is *cvector*.

 */
warp_Int
_warp_Refine(warp_Core     core,
             warp_Int      level,    /* fine level */
             warp_Int      f_index,  /* fine index */
             warp_Int      c_index,  /* coarse index */
             warp_Vector   cvector,
             warp_Vector  *fvector);

/**
 * Create a new grid object *grid_ptr* in core at *level*\n
 * *ilower* and *iupper* correspond to the lower and upper time index values
 * for this processor on this grid.
 */
warp_Int
_warp_GridInit(warp_Core     core,
               warp_Int      level,
               warp_Int      ilower,
               warp_Int      iupper,
               _warp_Grid  **grid_ptr);

/**
 * Destroy a Warp *grid*
 */
warp_Int
_warp_GridDestroy(warp_Core    core,
                  _warp_Grid  *grid);

/**
 * Set initial guess at C-points on *level*
 */
warp_Int
_warp_InitGuess(warp_Core  core,
                warp_Int   level);

/**
 * Do nu sweeps of F-then-C relaxation on *level*
 */
warp_Int
_warp_CFRelax(warp_Core  core,
              warp_Int   level);

/**
 * F-Relax on *level* and then restrict to *level+1*\n
 * Output:
 * - The restricted vectors *va* and *wa* will be created, representing 
 *    *level+1* versions of the unknown and rhs vectors.
 *    
 *    If set, the user-defined coarsening routine is called.
 *
 * - If *level==0*, then *rnorm_ptr* will contain the
 *   residual norm.
 */
warp_Int
_warp_FRestrict(warp_Core   core,       /**< warp_Core (_warp_Core) struct */   
                warp_Int    level,      /**< restrict from level to level+1 */
                warp_Int    iter,       /**< current iteration number (for user info) */
                warp_Real  *rnorm_ptr   /**< pointer to residual norm (if level 0) */
                );

/**
 * F-Relax on *level* and interpolate to *level-1*\n
 * Output:
 * - The unknown vector *u* on *level* is created by interpolating
 *   from *level+1*.
 *    
 *   If set, the user-defined refinement routine is called.
 */
warp_Int
_warp_FInterp(warp_Core  core,           /**< warp_Core (_warp_Core) struct */  
              warp_Int   level,          /**< interp from level to level+1 */
              warp_Int   iter,           /**< current iteration number (for user info) */
              warp_Real  rnorm           /**< residual norm (if level 0) */
              );

/**
 * Create a new fine grid based on user refinement factor information, then
 * F-relax and interpolate to the new fine grid and create a new multigrid
 * hierarchy.  In general, this will require load re-balancing as well.
 *
 * RDF: Todo, routine is unwritten
 */
warp_Int
_warp_FRefine(warp_Core   core,
              warp_Int   *refined_ptr);

/**
 * Write out the solution on grid *level* at Warp iteration *iter*.\n
 * *rnorm* denotes the last computed residual norm, and *done* is a boolean
 * indicating whether Warp has finished iterating and this is the last Write
 * call.
 */
warp_Int
_warp_FWrite(warp_Core     core,
             warp_Real     rnorm,
             warp_Int      iter,
             warp_Int      level,
             warp_Int      done);

/**
 * Initialize (and re-initialize) hierarchy
 */
warp_Int
_warp_InitHierarchy(warp_Core    core,
                    _warp_Grid  *fine_grid);

/**
 * Initialize a warp_Status structure in *status_ptr*, setting
 * the status values *rnorm*, *iter*, *level*, *done*
 **/
warp_Int
_warp_InitStatus(warp_Real        rnorm,
                 warp_Int         iter,
                 warp_Int         level,
                 warp_Int         done,
                 warp_Status     *status_ptr);

/**
 * Destroy a warp_Status structure
 **/
warp_Int
_warp_DestroyStatus(warp_Status  status);


#ifdef __cplusplus
}
#endif

#endif

