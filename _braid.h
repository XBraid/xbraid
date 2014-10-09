/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder schroder2@llnl.gov, Rob Falgout falgout2@llnl.gov,
 * Tzanio Kolev kolev1@llnl.gov, Ulrike Yang yang11@llnl.gov, 
 * Veselin Dobrev dobrev1@llnl.gov, et al. 
 * LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. Email schroder2@llnl.gov on how to download. 
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
 

/** \file _braid.h
 * \brief Define headers for developer routines.
 *
 * This file contains the headers for developer routines.
 */

#ifndef _braid_HEADER
#define _braid_HEADER

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <math.h>

#include "braid.h"

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------------
 * Main data structures and accessor macros
 *--------------------------------------------------------------------------*/

/**
 * XBraid comm handle structure
 *
 * Used for initiating and completing nonblocking communication to pass
 * braid_Vectors between processors.
 **/
typedef struct
{
   braid_Int       request_type;    /**< two values: recv type = 1, and send type = 0 */
   braid_Int       num_requests;    /**< number of active requests for this handle, usually 1 */
   MPI_Request    *requests;        /**< MPI request structure */
   MPI_Status     *status;          /**< MPI status */
   void           *buffer;          /**< Buffer for message */
   braid_Vector   *vector_ptr;      /**< braid_vector being sent/received */
   
} _braid_CommHandle;

/**
 * XBraid Accuracy Handle, used for controlling the accuracy of solves during
 * implicit time stepping.  For instance, to do less accurate solves on coarse
 * time grids
 **/
typedef struct
{
   braid_Int   matchF;
   braid_Real  value;      /**< accuracy value */
   braid_Real  old_value;  /**< old accuracy value used in FRestrict */
   braid_Real  loose;      /**< loose accuracy for spatial solves */
   braid_Real  tight;      /**< tight accuracy for spatial solves */
   braid_Int   tight_used; /**< tight accuracy used (1) or not (0) */
} _braid_AccuracyHandle;

/**
 * XBraid Grid structure for a certain time level
 *
 * Holds all the information for a processor related to the temporal
 * grid at this level.
 **/
typedef struct
{
   braid_Int          level;         /**< Level that grid is on */
   braid_Int          ilower;        /**< smallest time index at this level*/
   braid_Int          iupper;        /**< largest time index at this level*/
   braid_Int          clower;        /**< smallest C point index */
   braid_Int          cupper;        /**< largest C point index */
   braid_Int          gupper;        /**< global size of the grid */
   braid_Int          cfactor;       /**< coarsening factor */
   braid_Int          ncpoints;      /**< number of C points */
                  
   braid_Vector      *ua;            /**< unknown vectors            (C-points only)*/
   braid_Real        *ta;            /**< time values                (all points) */
   braid_Vector      *va;            /**< restricted unknown vectors (all points, NULL on level 0) */
   braid_Vector      *wa;            /**< rhs vectors f-v            (all points, NULL on level 0) */

   braid_Int          recv_index;    /**<  -1 means no receive */
   braid_Int          send_index;    /**<  -1 means no send */
   _braid_CommHandle *recv_handle;   /**<  Handle for nonblocking receives of braid_Vectors */
   _braid_CommHandle *send_handle;   /**<  Handle for nonblocking sends of braid_Vectors */

   braid_Vector      *ua_alloc;      /**< original memory allocation for ua */
   braid_Real        *ta_alloc;      /**< original memory allocation for ta */
   braid_Vector      *va_alloc;      /**< original memory allocation for va */
   braid_Vector      *wa_alloc;      /**< original memory allocation for wa */

} _braid_Grid;

/**
 * The typedef _braid_Core struct is a **critical** part of XBraid and 
 * is passed to *each* routine in XBraid.  It thus allows each routine access 
 * to XBraid attributes.
 **/
typedef struct _braid_Core_struct
{
   MPI_Comm               comm_world;
   MPI_Comm               comm;         /**< communicator for the time dimension */
   braid_Real             tstart;       /**< start time */
   braid_Real             tstop;        /**< stop time */
   braid_Int              ntime;        /**< initial number of time intervals */
   braid_App              app;          /**< application data for the user */
   
   braid_PtFcnPhi         phi;          /**< apply phi function */
   braid_PtFcnInit        init;         /**< return an initialized braid_Vector */
   braid_PtFcnClone       clone;        /**< clone a vector */
   braid_PtFcnFree        free;         /**< free up a vector */
   braid_PtFcnSum         sum;          /**< vector sum */
   braid_PtFcnSpatialNorm spatialnorm;  /**< Compute norm of a braid_Vector, this is a norm only over space */
   braid_PtFcnAccess      access;       /**< user access function to XBraid and current vector */
   braid_PtFcnBufSize     bufsize;      /**< return buffer size */
   braid_PtFcnBufPack     bufpack;      /**< pack a buffer */
   braid_PtFcnBufUnpack   bufunpack;    /**< unpack a buffer */
   braid_PtFcnCoarsen     coarsen;      /**< (optional) return a coarsened vector */
   braid_PtFcnRefine      refine;       /**< (optional) return a refined vector */

   braid_Int              access_level; /**< determines how often to call the user's access routine */ 
   braid_Int              print_level;  /**< determines amount of output printed to screem (0,1,2) */ 
   braid_Int              max_levels;   /**< maximum number of temporal grid levels */
   braid_Int              min_coarse;   /**< minimum possible coarse grid size */
   braid_Real             tol;          /**< stopping tolerance */
   braid_Int              rtol;         /**< use relative tolerance */
   braid_Int             *nrels;        /**< number of pre-relaxations on each level */
   braid_Int              nrdefault;    /**< default number of pre-relaxations */
   braid_Int             *cfactors;     /**< coarsening factors */
   braid_Int              cfdefault;    /**< default coarsening factor */
   braid_Int              max_iter;     /**< maximum number of multigrid in time iterations */
   braid_Int              niter;        /**< number of iterations */
   braid_Real             rnorm;        /**< residual norm */
   braid_Int              fmg;          /**< use FMG cycle */
   braid_Int              nfmg_Vcyc;    /**< number of V-cycle calls at each level in FMG */
   braid_Int              tnorm;        /**< choice of temporal norm */
   braid_Real            *tnorm_a;      /**< local array of residual norms on a proc's interval, used for inf-norm */
   _braid_AccuracyHandle *accuracy;     /**< accuracy of spatial solves on different levels */

   braid_AccessStatus     astatus;      /**< status structure passed to user-written Access routine */
   braid_CoarsenRefStatus cstatus;      /**< status structure passed to user-written coarsen/refine routines */
   braid_PhiStatus        pstatus;      /**< status structure passed to user-written phi routines */

   braid_Int              gupper;       /**< global upper index on the fine grid */

   braid_Int             *rfactors;     /**< refinement factors for finest grid (if any) */

   braid_Int              nlevels;      /**< number of temporal grid levels */
   _braid_Grid          **grids;        /**< pointer to temporal grid structures for each level*/

   braid_Real             localtime;    /**< local wall time for braid_Drive() */
   braid_Real             globaltime;   /**< global wall time for braid_Drive() */

} _braid_Core;

/*--------------------------------------------------------------------------
 * Accessor macros 
 *--------------------------------------------------------------------------*/

/** 
 * Accessor for _braid_CommHandle attributes 
 **/
#define _braid_CommHandleElt(handle, elt)  ((handle) -> elt)

/**
 * Accessor for _braid_Grid attributes 
 **/
#define _braid_GridElt(grid, elt)  ((grid) -> elt)

/** 
 * Accessor for _braid_Core attributes 
 **/
#define _braid_CoreElt(core, elt)     (  (core)  -> elt )

/** 
 * Accessor for _braid_Core functions 
 **/
#define _braid_CoreFcn(core, fcn)     (*((core)  -> fcn))

/*--------------------------------------------------------------------------
 * Print file for redirecting stdout when needed
 *--------------------------------------------------------------------------*/

/** 
 * This is the print file for redirecting stdout for all XBraid screen output
 **/
extern FILE *_braid_printfile;

/*--------------------------------------------------------------------------
 * Coarsening macros
 *--------------------------------------------------------------------------*/

/** 
 * Map a fine time index to a coarse time index, assumes a uniform coarsening
 * factor.
 **/
#define _braid_MapFineToCoarse(findex, cfactor, cindex) \
( cindex = (findex)/(cfactor) )

/** 
 * Map a coarse time index to a fine time index, assumes a uniform coarsening
 * factor.
 **/
#define _braid_MapCoarseToFine(cindex, cfactor, findex) \
( findex = (cindex)*(cfactor) )

/** 
 * Boolean, returns whether a time index is an F-point
 **/
#define _braid_IsFPoint(index, cfactor) \
( (index)%(cfactor) )

/** 
 * Boolean, returns whether a time index is an C-point
 **/
#define _braid_IsCPoint(index, cfactor) \
( !_braid_IsFPoint(index, cfactor) )

/*--------------------------------------------------------------------------
 * Prototypes
 *--------------------------------------------------------------------------*/

/**
 * Determine processor distribution.  This must agree with GetProc().
 * For the processor rank calling this function, it returns the smallest
 * and largest time indices ( *ilower_ptr* and *iupper_ptr*) that belong to 
 * that processor (the indices may * be F or C points).
 */
braid_Int
_braid_GetDistribution(braid_Core   core,
                       braid_Int   *ilower_ptr,
                       braid_Int   *iupper_ptr);

/**
 * Return the processor number in *proc_ptr* on which the time step *index* 
 * lives for the given *level*.  * Returns -1 if *index* is out of range
 */
braid_Int
_braid_GetProc(braid_Core   core,
               braid_Int    level,
               braid_Int    index,
               braid_Int   *proc_ptr);

/**
 * Initialize a receive to go into *vector_ptr* for the given time *index* on *level*.  
 * Also return a comm handle *handle_ptr* for querying later, to see if the receive has 
 * occurred.
 */
braid_Int
_braid_CommRecvInit(braid_Core           core,
                    braid_Int            level,
                    braid_Int            index,
                    braid_Vector        *vector_ptr,
                    _braid_CommHandle  **handle_ptr);

/**
 * Initialize a send of *vector* for the given time *index* on *level*.  
 * Also return a comm handle *handle_ptr* for querying later, to see if the 
 * send has occurred.
 */
braid_Int
_braid_CommSendInit(braid_Core           core,
                    braid_Int            level,
                    braid_Int            index,
                    braid_Vector         vector,
                    _braid_CommHandle  **handle_ptr);

/**
 * Block on the comm handle *handle_ptr* until the MPI operation (send or recv)
 * has completed
 */
braid_Int
_braid_CommWait(braid_Core          core,
               _braid_CommHandle **handle_ptr);

/**
 * Working on all intervals
 *
 * At *level*, post a receive for the point to the left of ilower (regardless 
 * whether ilower is F or C).  Then, post a send of iupper if iupper is a C 
 * point.
 */
braid_Int
_braid_UCommInit(braid_Core  core,
                 braid_Int   level);

/**
 * Working only on F-pt intervals
 *
 * At *level*, **only** post a receive for the point to the left of ilower 
 * if ilower is an F point.  Then, post a send of iupper if iupper is a C point.
 */
braid_Int
_braid_UCommInitF(braid_Core  core,
                  braid_Int   level);

/**
 * Finish up communication
 *
 * On *level*, wait on both the recv and send handles at this level.
 */
braid_Int
_braid_UCommWait(braid_Core  core,
                 braid_Int   level);

/**
 * Retrieve the time step indices at this *level* which correspond to the FC interval
 * given by *interval_index*.  *ci_ptr* is the time step index for the C point
 * and *flo_ptr* and *fhi_ptr* are the smallest and largest F point indices in this
 * interval.  *flo* = *ci* +1, and *fhi* = *ci* + coarsening_factor - 1
 */
braid_Int
_braid_UGetInterval(braid_Core   core,
                    braid_Int    level,
                    braid_Int    interval_index,
                    braid_Int   *flo_ptr,
                    braid_Int   *fhi_ptr,
                    braid_Int   *ci_ptr);

/**
 * Returns a reference to the local u-vector in *u_ptr* for the grid *level* at 
 * point *index*.  Caveat: if *index* is not a C-point and within my index range, 
 * NULL is returned.
 */
braid_Int
_braid_UGetVectorRef(braid_Core     core,
                     braid_Int      level,
                     braid_Int      index,
                     braid_Vector  *u_ptr);

/**
 * Stores a reference to the vector *u* on grid *level* at point *index*.
 * If *index* is not a C-point and within this processor's range of time points, 
 * then nothing is done.
 */
braid_Int
_braid_USetVectorRef(braid_Core    core,
                     braid_Int     level,
                     braid_Int     index,
                     braid_Vector  u);

/**
 * Returns the u-vector in *u_ptr* on grid *level* at point *index*.  If *index* is my
 * "receive index" (as set by UCommInit(), for example), the u-vector  will be
 * received from a neighbor processor.  If *index* is within my index range and
 * is also a C-point, the saved value of u will be used.  A NULL value is
 * returned otherwise.
 */
braid_Int
_braid_UGetVector(braid_Core     core,
                  braid_Int      level,
                  braid_Int      index,
                  braid_Vector  *u_ptr);

/**
 * Sets the u-vector on grid *level* at point *index*.  If *index* is my "send
 * index" (as set by UCommInit(), for example), a send is initiated to a neighbor 
 * processor.  If *index* is within my index range and is also a C-point, the 
 * value is saved locally.
 */
braid_Int
_braid_USetVector(braid_Core    core,
                  braid_Int     level,
                  braid_Int     index,
                  braid_Vector  u);

/** 
 * Call user's access function in order to give access to XBraid and the
 * current vector.  Most commonly, this lets the user write *u* to screen,
 * disk, etc...  The vector *u* corresponds to time step *index* on *level*.
 * *status* holds state information about the current XBraid iteration, time
 * value, etc...
 */
braid_Int
_braid_UAccessVector(braid_Core         core,
                     braid_AccessStatus status,
                     braid_Vector       u);

/**
 * Apply Phi to the vector *u*
 *
 * This is the vector corresponding to the time step *index* on *level*.
 * *accuracy* is a user set variable to allow for tuning the accuracy of 
 * spatial solvesfor implicit stepping. And, *rfactor* allows the user to
 * subdivide time intervals for accuracy purposes.
 */
braid_Int
_braid_Phi(braid_Core     core,
           braid_Int      level,
           braid_Int      index,
           braid_Real     accuracy,
           braid_Vector   u,
           braid_Int     *rfactor);


/**
 * Integrate one time step at time step *index* to time step *index*+1
 *
 */
braid_Int
_braid_Step(braid_Core     core,
            braid_Int      level,
            braid_Int      index,
            braid_Real     accuracy,
            braid_Vector   u);

/**
 * Coarsen in space on *level* by calling the user's coarsen function.
 * The vector corresponding to the time step index *f_index* on the fine 
 * grid is coarsened to the time step index *c_index* on the coarse grid.
 * The output goes in *cvector* and the input vector is *fvector*.
 */
braid_Int
_braid_Coarsen(braid_Core     core,
               braid_Int      level,    /* coarse level */
               braid_Int      f_index,  /* fine index */
               braid_Int      c_index,  /* coarse index */
               braid_Vector   fvector,
               braid_Vector  *cvector);

/**
 * Refine in space on *level* by calling the user's refine function.
 * The vector corresponding to the time step index *c_index* on the coarse
 * grid is refined to the time step index *f_index* on the fine grid.
 * The output goes in *fvector* and the input vector is *cvector*.

 */
braid_Int
_braid_Refine(braid_Core     core,
              braid_Int      level,    /* fine level */
              braid_Int      f_index,  /* fine index */
              braid_Int      c_index,  /* coarse index */
              braid_Vector   cvector,
              braid_Vector  *fvector);

/**
 * Create a new grid object *grid_ptr* in core at *level*
 *
 * *ilower* and *iupper* correspond to the lower and upper time index values
 * for this processor on this grid.
 */
braid_Int
_braid_GridInit(braid_Core     core,
                braid_Int      level,
                braid_Int      ilower,
                braid_Int      iupper,
                _braid_Grid  **grid_ptr);

/**
 * Destroy an XBraid *grid*
 */
braid_Int
_braid_GridDestroy(braid_Core    core,
                   _braid_Grid  *grid);

/**
 * Set initial guess at C-points on *level*
 */
braid_Int
_braid_InitGuess(braid_Core  core,
                 braid_Int   level);

/**
 * Do nu sweeps of F-then-C relaxation on *level*
 */
braid_Int
_braid_FCRelax(braid_Core  core,
               braid_Int   level);

/**
 * F-Relax on *level* and then restrict to *level+1*
 * 
 * The output is set in the braid_Grid in core, so that the restricted 
 * vectors *va* and *wa* will be created, representing 
 * *level+1* versions of the unknown and rhs vectors.
 * 
 * If the user has set spatial coarsening, then this user-defined routine is 
 * also called.
 *
 * If *level==0*, then *rnorm_ptr* will contain the residual norm.
 */
braid_Int
_braid_FRestrict(braid_Core   core,       /**< braid_Core (_braid_Core) struct */   
                 braid_Int    level,      /**< restrict from level to level+1 */
                 braid_Int    iter,       /**< current iteration number (for user info) */
                 braid_Real  *rnorm_ptr   /**< pointer to residual norm (if level 0) */
                 );

/** F-Relax on *level* and interpolate to *level-1*
 *
 * The output is set in the braid_Grid in core, so that the 
 * vector *u* on *level* is created by interpolating from *level+1*.
 *
 * If the user has set spatial refinement, then this user-defined routine is 
 * also called.
 */
braid_Int
_braid_FInterp(braid_Core  core,           /**< braid_Core (_braid_Core) struct */  
               braid_Int   level,          /**< interp from level to level+1 */
               braid_Int   iter,           /**< current iteration number (for user info) */
               braid_Real  rnorm           /**< residual norm (if level 0) */
               );

/**
 * Create a new fine grid based on user refinement factor information, then
 * F-relax and interpolate to the new fine grid and create a new multigrid
 * hierarchy.  In general, this will require load re-balancing as well.
 *
 * RDF: Todo, routine is unwritten
 */
braid_Int
_braid_FRefine(braid_Core   core,
               braid_Int   *refined_ptr);

/** 
 * Call the user's access function in order to give access to XBraid and 
 * the current vector at grid *level and iteration *iter*.  Most commonly, 
 * this lets the user write solutions to screen, disk, etc... 
 * The quantity *rnorm* denotes the last computed residual
 * norm, and *done* is a boolean indicating whether XBraid has finished
 * iterating and this is the last Access call. 
 */
braid_Int
_braid_FAccess(braid_Core     core,
               braid_Real     rnorm,
               braid_Int      iter,
               braid_Int      level,
               braid_Int      done);

/**
 * Initialize (and re-initialize) hierarchy
 */
braid_Int
_braid_InitHierarchy(braid_Core    core,
                     _braid_Grid  *fine_grid);


#ifdef __cplusplus
}
#endif

#endif

