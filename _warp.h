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
#include <stdio.h>
#include <math.h>

#include "warp.h"

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------------
 * Main data structures and accessor macros
 *--------------------------------------------------------------------------*/

typedef struct
{
   warp_Int     request_type; /* recv type = 1 */
   warp_Int     num_requests;
   MPI_Request *requests;
   MPI_Status  *status;
   void        *buffer;
   warp_Vector *vector_ptr;
   
} _warp_CommHandle;

typedef struct
{
   warp_Int   matchF;
   warp_Float value;      /* accuracy value */
   warp_Float old_value;  /* old accuracy value used in FRestrict */
   warp_Float loose;      /* loose accuracy for spatial solves */
   warp_Float tight;      /* tight accuracy for spatial solves */
   warp_Int   tight_used; /* tight accuracy used (1) or not (0) */
} _warp_AccuracyHandle;

typedef struct
{
   warp_Int     level;
   warp_Int     ilower, iupper;
   warp_Int     clower, cupper, cfactor, ncpoints;

   warp_Vector *ua;  /* unknown vectors            (C-points only)*/
   warp_Float  *ta;  /* time values                (all points) */
   warp_Vector *va;  /* restricted unknown vectors (all points, NULL on level 0) */
   warp_Vector *wa;  /* rhs vectors f-v            (all points, NULL on level 0) */

   warp_Int          recv_index; /* -1 means no receive */
   warp_Int          send_index; /* -1 means no send */
   _warp_CommHandle *recv_handle;
   _warp_CommHandle *send_handle;

   /* pointers to the original memory allocation for ua, ta, va, and wa */
   warp_Vector *ua_alloc;
   warp_Float  *ta_alloc;
   warp_Vector *va_alloc;
   warp_Vector *wa_alloc;

} _warp_Grid;

/**
 * The typedef _warp_Core struct is a **critical** part of warp and 
 * is passed to *each* routine in warp.  It thus allows each routine access 
 * to the basic warp attributes such as the user's
 * - phi (\f$\phi \f$) function where
 *   \f{equation*}{
 *   \phi(\mathbf{u}_i) = \mathbf{u}_{i+1} \f} 
 *   evolves the solution
 * - warp_App structure
 *
 **/
typedef struct _warp_Core_struct
{
   MPI_Comm              comm_world;
   MPI_Comm              comm;      /**< communicator for the time dimension */
   warp_Float            tstart;    /**< start time */
   warp_Float            tstop;     /**< stop time */
   warp_Int              ntime;     /**< initial number of time intervals */
   warp_App              app;       /**< application data for the user */
   
   warp_PtFcnPhi         phi;       /**< apply phi function */
   warp_PtFcnInit        init;      /**< return an initial solution vector */
   warp_PtFcnClone       clone;     /**< clone a vector */
   warp_PtFcnFree        free;      /**< free up a vector */
   warp_PtFcnSum         sum;       /**< vector sum */
   warp_PtFcnDot         dot;       /**< dot product */
   warp_PtFcnWrite       write;     /**< write the vector */
   warp_PtFcnBufSize     bufsize;   /**< return buffer size */
   warp_PtFcnBufPack     bufpack;   /**< pack a buffer */
   warp_PtFcnBufUnpack   bufunpack; /**< unpack a buffer */
   warp_PtFcnCoarsen     coarsen;   /**< (optional) return a coarsened vector */
   warp_PtFcnRefine      refine;    /**< (optional) return a refined vector */

   warp_Int              max_levels;/**< maximum number of temporal grid levels */
   warp_Float            tol;       /**< stopping tolerance */
   warp_Int              rtol;      /**< use relative tolerance */
   warp_Int             *nrels;     /**< number of pre-relaxations on each level */
   warp_Int              nrdefault; /**< default number of pre-relaxations */
   warp_Int             *cfactors;  /**< coarsening factors */
   warp_Int              cfdefault; /**< default coarsening factor */
   warp_Int              max_iter;  /**< maximum number of multigrid in time iterations */
   warp_Int              niter;     /**< number of iterations */
   warp_Float            rnorm;     /**< residual norm */
   warp_Int              fmg;       /**< use FMG cycle */
   _warp_AccuracyHandle *accuracy;  /**< accuracy of spatial solves on different levels */

   warp_Int              gupper;    /**< global upper index on the fine grid */

   warp_Int             *rfactors;  /**< refinement factors for finest grid (if any) */

   warp_Int              nlevels;   /**< number of temporal grid levels */
   _warp_Grid          **grids;     /**< pointer to temporal grid structures for each level*/

} _warp_Core;

/* Accessor macros */

#define _warp_CommHandleElt(handle, elt)  ((handle) -> elt)

#define _warp_GridElt(grid, elt)  ((grid) -> elt)

#define _warp_CoreElt(core, elt)  (  (core) -> elt )
#define _warp_CoreFcn(core, fcn)  (*((core) -> fcn))

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
 */
warp_Int
_warp_GetDistribution(warp_Core   core,
                      warp_Int   *ilower_ptr,
                      warp_Int   *iupper_ptr);

/**
 * Returns -1 if index is out of range
 */
warp_Int
_warp_GetProc(warp_Core   core,
              warp_Int    level,
              warp_Int    index,
              warp_Int   *proc_ptr);

/**
 * Blah..
 */
warp_Int
_warp_CommRecvInit(warp_Core           core,
                   warp_Int            level,
                   warp_Int            index,
                   warp_Vector        *vector_ptr,
                   _warp_CommHandle  **handle_ptr);

/**
 * Blah..
 */
warp_Int
_warp_CommSendInit(warp_Core           core,
                   warp_Int            level,
                   warp_Int            index,
                   warp_Vector         vector,
                   _warp_CommHandle  **handle_ptr);

/**
 * Blah..
 */
warp_Int
_warp_CommWait(warp_Core          core,
               _warp_CommHandle **handle_ptr);

/**
 * Working on all intervals
 */
warp_Int
_warp_UCommInit(warp_Core  core,
                warp_Int   level);

/**
 * Working only on F-pt intervals
 */
warp_Int
_warp_UCommInitF(warp_Core  core,
                 warp_Int   level);

/**
 * Finish up communication
 */
warp_Int
_warp_UCommWait(warp_Core  core,
                warp_Int   level);

/**
 * Blah..
 */
warp_Int
_warp_UGetInterval(warp_Core   core,
                   warp_Int    level,
                   warp_Int    interval_index,
                   warp_Int   *flo_ptr,
                   warp_Int   *fhi_ptr,
                   warp_Int   *ci_ptr);

/**
 * Returns a reference to the local u-vector on grid 'level' at point 'index'.
 * If 'index' is not a C-point and within my index range, NULL is returned.
 */
warp_Int
_warp_UGetVectorRef(warp_Core     core,
                    warp_Int      level,
                    warp_Int      index,
                    warp_Vector  *u_ptr);

/**
 * Stores a reference to the u-vector on grid 'level' at point 'index'.
 * If 'index' is not a C-point and within my index range, nothing is done.
 */
warp_Int
_warp_USetVectorRef(warp_Core    core,
                    warp_Int     level,
                    warp_Int     index,
                    warp_Vector  u);

/**
 * Returns the u-vector on grid 'level' at point 'index'.  If 'index' is my
 * "receive index" (as set by UCommInit(), for example), the u-vector will be
 * received from a neighbor processor.  If 'index' is within my index range and
 * is also a C-point, the saved value of u will be used.  A NULL value is
 * returned otherwise.
 */
warp_Int
_warp_UGetVector(warp_Core     core,
                 warp_Int      level,
                 warp_Int      index,
                 warp_Vector  *u_ptr);

/**
 * Sets the u-vector on grid 'level' at point 'index'.  If 'index' is my "send
 * index", a send is initiated to a neighbor processor.  If 'index' is within my
 * index range and is also a C-point, the value is saved locally.
 */
warp_Int
_warp_USetVector(warp_Core    core,
                 warp_Int     level,
                 warp_Int     index,
                 warp_Vector  u);

/**
 * Blah..
 */
warp_Int
_warp_UWriteVector(warp_Core    core,
                   warp_Int     level,
                   warp_Int     index,
                   warp_Vector  u);

/**
* Apply Phi
 */
warp_Int
_warp_Phi(warp_Core     core,
          warp_Int      level,
          warp_Int      index,
          warp_Float    accuracy,
          warp_Int      gzero,
          warp_Vector   u,
          warp_Int     *rfactor);


/**
 * Integrate one time step
 */
warp_Int
_warp_Step(warp_Core     core,
           warp_Int      level,
           warp_Int      index,
           warp_Float    accuracy,
           warp_Vector   u);

/**
 * Coarsen in space
 */
warp_Int
_warp_Coarsen(warp_Core     core,
              warp_Int      level,    /* coarse level */
              warp_Int      f_index,  /* fine index */
              warp_Int      c_index,  /* coarse index */
              warp_Vector   fvector,
              warp_Vector  *cvector);

/**
 * Refine in space
 */
warp_Int
_warp_Refine(warp_Core     core,
             warp_Int      level,    /* fine level */
             warp_Int      f_index,  /* fine index */
             warp_Int      c_index,  /* coarse index */
             warp_Vector   cvector,
             warp_Vector  *fvector);

/**
 * Create a new grid object
 */
warp_Int
_warp_GridInit(warp_Core     core,
               warp_Int      level,
               warp_Int      ilower,
               warp_Int      iupper,
               _warp_Grid  **grid_ptr);

/**
 * Blah..
 */
warp_Int
_warp_GridDestroy(warp_Core    core,
                  _warp_Grid  *grid);

/**
 * Set initial guess at C-points
 */
warp_Int
_warp_InitGuess(warp_Core  core,
                warp_Int   level);

/**
 * Do nu sweeps of F-then-C relaxation
 */
warp_Int
_warp_CFRelax(warp_Core  core,
              warp_Int   level);

/**
 * F-Relax on level and restrict to level+1
 */
warp_Int
_warp_FRestrict(warp_Core    core,
                warp_Int     level,
                warp_Float  *rnorm_ptr);


/**
 * F-Relax on level and interpolate to level-1
 */
warp_Int
_warp_FInterp(warp_Core  core,
              warp_Int   level);

/**
 * Create a new fine grid based on user refinement factor information, then
 * F-relax and interpolate to the new fine grid and create a new multigrid
 * hierarchy.  In general, this will require load re-balancing as well.
 *
 * RDF: Todo
 */
warp_Int
_warp_FRefine(warp_Core   core,
              warp_Int   *refined_ptr);

/**
 * Write out the solution on grid level
 */
warp_Int
_warp_FWrite(warp_Core  core,
             warp_Int   level);

/**
 * Initialize (and re-initialize) hierarchy
 */
warp_Int
_warp_InitHierarchy(warp_Core    core,
                    _warp_Grid  *fine_grid);


#ifdef __cplusplus
}
#endif

#endif

