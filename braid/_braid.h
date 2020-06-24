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

/** \file _braid.h
 * \brief Define headers for XBraid internal (developer) routines and XBraid internal structure declarations.
 *
 * This file contains the headers for XBraid internal (developer) routines and
 * structure declarations.
 */

#ifndef _braid_HEADER
#define _braid_HEADER

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <math.h>

#include "braid.h"
#include "tape.h"

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------------
 * Error handling
 *--------------------------------------------------------------------------*/

/** 
 * This is the global XBraid error flag.  If it is ever nonzero, an error has 
 * occurred. 
 **/
extern braid_Int _braid_error_flag;

void _braid_ErrorHandler(const char *filename, braid_Int line, braid_Int ierr, const char *msg);
#define _braid_Error(IERR, msg)       _braid_ErrorHandler(__FILE__, __LINE__, IERR, msg)
#define _braid_ErrorInArg(IARG, msg)  _braid_Error(HYPRE_ERROR_ARG | IARG<<3, msg)

/*--------------------------------------------------------------------------
 * Memory allocation macros
 *--------------------------------------------------------------------------*/

/** 
 * Allocation macro 
 **/
#define _braid_TAlloc(type, count) \
( (type *)malloc((size_t)(sizeof(type) * (count))) )

/** 
 * Allocation macro 
 **/
#define _braid_CTAlloc(type, count) \
( (type *)calloc((size_t)(count), (size_t)sizeof(type)) )

/** 
 * Re-allocation macro 
 **/
#define _braid_TReAlloc(ptr, type, count) \
( (type *)realloc((char *)ptr, (size_t)(sizeof(type) * (count))) )

/** 
 * Free memory macro 
 **/
#define _braid_TFree(ptr) \
( free((char *)ptr), ptr = NULL )

/*--------------------------------------------------------------------------
 * Miscellaneous macros and functions 
 *--------------------------------------------------------------------------*/

#ifndef _braid_max
#define _braid_max(a,b)  (((a)<(b)) ? (b) : (a))
#endif
#ifndef _braid_min
#define _braid_min(a,b)  (((a)<(b)) ? (a) : (b))
#endif
#ifndef _braid_isnan
#define _braid_isnan(a) (a != a)
#endif

/* Used to implement periodic feature */
#define _braid_SendIndexNull -2
#define _braid_RecvIndexNull -2
#define _braid_MapPeriodic(index, npoints) \
( index = ((index)+(npoints)) % (npoints) )  /* this also handles negative indexes */

/** 
 * Braid Vector Structures:
 *
 * There are three vector structures
 *   _braid_VectorBar      Defined below
 *   braid_Vector          Defined in braid.h
 *   braid_BaseVector      Defined below
 *
 * The braid_BaseVector is the main internal Vector class, which is 
 * stored at each time point.  It basically wraps the Vector and 
 * braid_VectorBar (see below).  The braid_VectorBar is only used if the
 * adjoint capability is used, when it stores adjoint variables.  It's 
 * basically a smart pointer wrapper around a braid_Vector. Note that it
 * is always the braid_Vector that's passed to user-routines.   
 *
 */

/**
 * Shared pointer implementation for storing the intermediat AD-bar variables while taping.
 * This is essentially the same as a userVector, except we need shared pointer 
 * capabilities to know when to delete.
 */
struct _braid_VectorBar_struct
{
   braid_Vector userVector;         /**< holds the u_bar data */
   braid_Int    useCount;           /**< counts the number of pointers to this struct */
}; 
typedef struct _braid_VectorBar_struct *braid_VectorBar;

/** 
 * Braid vector used for storage of all state and (if needed) adjoint
 * information.  Stores both the user's primal vector (braid_Vector type) and
 * the associated bar vector (braid_VectorBar type) if the adjoint
 * functionality is being used.  If adjoint is not being used, bar==NULL. 
 */
struct _braid_BaseVector_struct
{
   braid_Vector    userVector;      /**< holds the users primal vector */
   braid_VectorBar bar;             /**< holds the bar vector (shared pointer implementation) */
};
typedef struct _braid_BaseVector_struct *braid_BaseVector;

/** 
 * Data structure for storing the optimization variables
 */
struct _braid_Optimization_struct
{
   braid_Real       sum_user_obj;     /**< sum of user's objective function values over time */
   braid_Real       objective;        /**< global objective function value */
   braid_Real       tstart_obj;       /**< starting time for evaluating the user's local objective */
   braid_Real       tstop_obj;        /**< stopping time for evaluating the user's local objective */
   braid_Real       f_bar;            /**< contains the seed for tape evaluation  */
   braid_Real       rnorm_adj;        /**< norm of the adjoint residual */
   braid_Real       rnorm0_adj;       /**< initial norm of the adjoint residual */
   braid_Real       rnorm;            /**< norm of the state residual */
   braid_Real       rnorm0;           /**< initial norm of the state residual */
   braid_Real       tol_adj;          /**< tolerance of adjoint residual */
   braid_Int        rtol_adj;         /**< flag: use relative tolerance for adjoint */
   braid_Vector    *adjoints;         /**< vector for the adjoint variables */
   braid_VectorBar *tapeinput;        /**< helper: store pointer to input of one braid iteration */
   void            *sendbuffer;       /**< helper: Memory for BufUnPackDiff communication */
   MPI_Request     *request;          /**< helper: Storing the MPI request of BufUnPackDiff */
};
typedef struct _braid_Optimization_struct *braid_Optim;

/*--------------------------------------------------------------------------
 * Main data structures and accessor macros
 *--------------------------------------------------------------------------*/

/**
 * XBraid comm handle structure
 *
 * Used for initiating and completing nonblocking communication to pass
 * braid_BaseVectors between processors.
 **/
typedef struct
{
   braid_Int         request_type;    /**< two values: recv type = 1, and send type = 0 */
   braid_Int         num_requests;    /**< number of active requests for this handle, usually 1 */
   MPI_Request      *requests;        /**< MPI request structure */
   MPI_Status       *status;          /**< MPI status */
   void             *buffer;          /**< Buffer for message */
   braid_BaseVector *vector_ptr;      /**< braid_vector being sent/received */
   
} _braid_CommHandle;

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

   braid_Int          nupoints;      /**< number of unknown vector points */
   braid_BaseVector  *ua;            /**< unknown vectors            (C-points at least)*/
   braid_Real        *ta;            /**< time values                (all points) */
   braid_BaseVector  *va;            /**< restricted unknown vectors (all points, NULL on level 0) */
   braid_BaseVector  *fa;            /**< rhs vectors f              (all points, NULL on level 0) */

   braid_Int          recv_index;    /**<  -1 means no receive */
   braid_Int          send_index;    /**<  -1 means no send */
   _braid_CommHandle *recv_handle;   /**<  Handle for nonblocking receives of braid_BaseVectors */
   _braid_CommHandle *send_handle;   /**<  Handle for nonblocking sends of braid_BaseVectors */

   braid_BaseVector  *ua_alloc;      /**< original memory allocation for ua */
   braid_Real        *ta_alloc;      /**< original memory allocation for ta */
   braid_BaseVector  *va_alloc;      /**< original memory allocation for va */
   braid_BaseVector  *fa_alloc;      /**< original memory allocation for fa */

} _braid_Grid;

/**
 * The typedef _braid_Core struct is a **critical** part of XBraid and 
 * is passed to *each* routine in XBraid.  It thus allows each routine access 
 * to XBraid attributes.
 **/
typedef struct _braid_Core_struct
{
   MPI_Comm               comm_world;
   MPI_Comm               comm;             /**< communicator for the time dimension */
   braid_Int              myid_world;       /**< my rank in the world communicator */
   braid_Int              myid;             /**< my rank in the time communicator */
   braid_Real             tstart;           /**< start time */
   braid_Real             tstop;            /**< stop time */
   braid_Int              ntime;            /**< initial number of time intervals */
   braid_App              app;              /**< application data for the user */

   braid_PtFcnStep        step;             /**< apply step function */
   braid_PtFcnInit        init;             /**< return an initialized braid_BaseVector */
   braid_PtFcnSInit       sinit;            /**< (optional) return an initialized shell of braid_BaseVector */
   braid_PtFcnClone       clone;            /**< clone a vector */
   braid_PtFcnSClone      sclone;           /**< (optional) clone the shell of a vector */
   braid_PtFcnFree        free;             /**< free up a vector */
   braid_PtFcnSFree       sfree;            /**< (optional) free up the data of a vector, keep the shell */
   braid_PtFcnSum         sum;              /**< vector sum */
   braid_PtFcnSpatialNorm spatialnorm;      /**< Compute norm of a braid_BaseVector, this is a norm only over space */
   braid_PtFcnAccess      access;           /**< user access function to XBraid and current vector */
   braid_PtFcnBufSize     bufsize;          /**< return buffer size */
   braid_PtFcnBufPack     bufpack;          /**< pack a buffer */
   braid_PtFcnBufUnpack   bufunpack;        /**< unpack a buffer */
   braid_PtFcnResidual    residual;         /**< (optional) compute residual */
   braid_PtFcnSCoarsen    scoarsen;         /**< (optional) return a spatially coarsened vector */
   braid_PtFcnSRefine     srefine;          /**< (optional) return a spatially refined vector */
   braid_PtFcnSync        sync;             /**< (optional) user access to app once-per-processor */
   braid_PtFcnTimeGrid    tgrid;            /**< (optional) return time point values on level 0 */
   braid_Int              periodic;         /**< determines if periodic */
   braid_Int              initiali;         /**< initial condition grid index (0: default; -1: periodic ) */

   braid_Int              access_level;     /**< determines how often to call the user's access routine */ 
   braid_Int              print_level;      /**< determines amount of output printed to screen (0,1,2,3) */
   braid_Int              io_level;         /**< determines amount of output printed to files (0,1) */
   braid_Int              seq_soln;         /**< boolean, controls if the initial guess is from sequential time stepping*/
   braid_Int              max_levels;       /**< maximum number of temporal grid levels */
   braid_Int              incr_max_levels;  /**< After doing refinement, increase the max number of levels by 1 (0=false, 1=true)*/
   braid_Int              min_coarse;       /**< minimum possible coarse grid size */
   braid_Real             tol;              /**< stopping tolerance */
   braid_Int              rtol;             /**< use relative tolerance */
   braid_Int             *nrels;            /**< number of pre-relaxations on each level */
   braid_Int              nrdefault;        /**< default number of pre-relaxations */
   braid_Real            *CWts;             /**< C-relaxation weight for each level */
   braid_Real             CWt_default;      /**< default C-relaxtion weight */
   braid_Int             *cfactors;         /**< coarsening factors */
   braid_Int              cfdefault;        /**< default coarsening factor */
   braid_Int              max_iter;         /**< maximum number of multigrid in time iterations */
   braid_Int              niter;            /**< number of iterations */
   braid_Int              fmg;              /**< use FMG cycle */
   braid_Int              nfmg;             /**< number of fmg cycles to do initially before switching to V-cycles */
   braid_Int              nfmg_Vcyc;        /**< number of V-cycle calls at each level in FMG */
   braid_Int              warm_restart;     /**< boolean, indicates whether this is a warm restart of an existing braid_Core */
   braid_Int              tnorm;            /**< choice of temporal norm */
   braid_Real            *tnorm_a;          /**< local array of residual norms on a proc's interval, used for inf-norm */
   braid_Real             rnorm0;           /**< initial residual norm */
   braid_Real            *rnorms;           /**< residual norm history */
   braid_PtFcnResidual    full_rnorm_res;   /**< (optional) used to compute full residual norm */
   braid_Real             full_rnorm0;      /**< (optional) initial full residual norm */
   braid_Real            *full_rnorms;      /**< (optional) full residual norm history */

   braid_Int              storage;          /**< storage = 0 (C-points), = 1 (all) */
   braid_Int              useshell;         /**< activate the shell structure of vectors */

   braid_Int              gupper;           /**< global size of the fine grid */

   braid_Int              refine;           /**< refine in time (refine = 1) */
   braid_Int             *rfactors;         /**< refinement factors for finest grid (if any) */
   braid_Real           **rdtvalues;        /**< Array of pointers to arrays of dt values for non-uniform refinement  */
   braid_Int              r_space;          /**< spatial refinement flag */
   braid_Int              rstopped;         /**< refinement stopped at iteration rstopped */
   braid_Int              nrefine;          /**< number of refinements done */
   braid_Int              max_refinements;  /**< maximum number of refinements */
   braid_Int              tpoints_cutoff;   /**< refinements halt after the number of time steps exceed this value */

   braid_Int              skip;             /**< boolean, controls skipping any work on first down-cycle */

   braid_Int              nlevels;          /**< number of temporal grid levels */
   _braid_Grid          **grids;            /**< pointer to temporal grid structures for each level*/

   braid_Real             localtime;        /**< local wall time for braid_Drive() */
   braid_Real             globaltime;       /**< global wall time for braid_Drive() */

   /* Data for adjoint and optimization */
   braid_Optim            optim;             /**< structure that stores optimization variables (objective function, etc.) */ 
   braid_Int              adjoint;           /**< determines if adjoint run is performed (1) or not (0) */
   braid_Int              record;            /**< determines if actions are recorded to the tape or not.  This separate 
                                                  flag from adjoint is needed, because the final FAccess call should 
                                                  not be recorded unless nlevels==1, but the adjoint flag must be true 
                                                  even if nlevels==1. */
   braid_Int              obj_only;          /**< determines if adjoint code computes ONLY objective, no gradients. */
   braid_Int              verbose_adj;       /**< verbosity of the adjoint tape, displays the actions that are pushed / popped to the tape*/

   _braid_Tape*          actionTape;         /**< tape storing the actions while recording */
   _braid_Tape*          userVectorTape;     /**< tape storing primal braid_vectors while recording */
   _braid_Tape*          barTape;            /**< tape storing intermediate AD-bar variables while recording */
      
   braid_PtFcnObjectiveT                objectiveT;           /**< User function: evaluate objective function at time t */
   braid_PtFcnStepDiff                  step_diff;            /**< User function: apply differentiated step function */
   braid_PtFcnObjectiveTDiff            objT_diff;            /**< User function: apply differentiated objective function */
   braid_PtFcnResetGradient             reset_gradient;       /**< User function: Set the gradient to zero. Is called before each iteration */
   braid_PtFcnPostprocessObjective      postprocess_obj;      /**< Optional user function: Modify the time-averaged objective function, e.g. for inverse design problems, adding relaxation term etc. */
   braid_PtFcnPostprocessObjective_diff postprocess_obj_diff; /**< Optional user function: Derivative of postprocessing function  */

   /** Data elements required for the Status structures */
   /** Common Status properties */
   braid_Real    t;                /**< current time */
   braid_Int     idx;              /**< time point index value corresponding to t on the global time grid */
   braid_Int     level;            /**< current level in XBraid*/
   /** AccessStatus properties */
   braid_Real    rnorm;            /**< residual norm */
   braid_Int     done;             /**< boolean describing whether XBraid has finished */
   braid_Int     wrapper_test;     /**< boolean describing whether this call is only a wrapper test */
   braid_Int     calling_function; /**< from which function are we accessing the vector */
   /** CoarsenRefStatus properties*/
   braid_Real    f_tprior;         /**< time value to the left of tstart on fine grid */
   braid_Real    f_tstop;          /**< time value to the right of tstart  on fine grid */
   braid_Real    c_tprior;         /**< time value to the left of tstart on coarse grid */
   braid_Real    c_tstop;          /**< time value to the right of tstart on coarse grid */
   /** StepStatus properties */
   braid_Real    tnext;            /**< time value to evolve towards, time value to the right of tstart */
   braid_Real    old_fine_tolx;    /**< Allows for storing the previously used fine tolerance from GetSpatialAccuracy */
   braid_Int     tight_fine_tolx;  /**< Boolean, indicating whether the tightest fine tolx has been used, condition for halting */
   braid_Int     rfactor;          /**< if set by user, allows for subdivision of this interval for better time accuracy */
   /** BufferStatus properties */
   braid_Int    messagetype;       /**< message type, 0: for Step(), 1: for load balancing */
   braid_Int    size_buffer;       /**< if set by user, send buffer will be "size" bytes in length */
   braid_Int    send_recv_rank;    /***< holds the rank of the source / receiver from MPI_Send / MPI_Recv calls. */
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

/** 
 * Returns the index for the next C-point to the right of index (inclusive)
 **/
#define _braid_NextCPoint(index, cfactor) \
( ((braid_Int)((index)+(cfactor)-1)/(cfactor))*(cfactor) )

/**
 * Returns the index for the previous C-point to the left of index (inclusive)
 **/
#define _braid_PriorCPoint(index, cfactor) \
( ((braid_Int)(index)/(cfactor))*(cfactor) )

/*--------------------------------------------------------------------------
 * Prototypes
 *--------------------------------------------------------------------------*/

/* distribution.c */

/**
 * Returns the index interval for *proc* in a blocked data distribution.
 */
braid_Int
_braid_GetBlockDistInterval(braid_Int   npoints,
                            braid_Int   nprocs,
                            braid_Int   proc,
                            braid_Int  *ilower_ptr,
                            braid_Int  *iupper_ptr);

/**
 * Returns the processor that owns *index* in a blocked data distribution
 * (returns -1 if *index* is out of range).
 */
braid_Int
_braid_GetBlockDistProc(braid_Int   npoints,
                        braid_Int   nprocs,
                        braid_Int   index,
                        braid_Int   periodic,
                        braid_Int  *proc_ptr);

/**
 * Returns the index interval for my processor on the finest grid level.
 * For the processor rank calling this function, it returns the smallest
 * and largest time indices ( *ilower_ptr* and *iupper_ptr*) that belong to 
 * that processor (the indices may be F or C points).
 */
braid_Int
_braid_GetDistribution(braid_Core   core,
                       braid_Int   *ilower_ptr,
                       braid_Int   *iupper_ptr);

/**
 * Returns the processor number in *proc_ptr* on which the time step *index*
 * lives for the given *level*.  Returns -1 if *index* is out of range.
 */
braid_Int
_braid_GetProc(braid_Core   core,
               braid_Int    level,
               braid_Int    index,
               braid_Int   *proc_ptr);

/* communication.c */

/**
 * Initialize a receive to go into *vector_ptr* for the given time *index* on
 * *level*.  Also return a comm handle *handle_ptr* for querying later, to see
 * if the receive has occurred.
 */
braid_Int
_braid_CommRecvInit(braid_Core           core,
                    braid_Int            level,
                    braid_Int            index,
                    braid_BaseVector    *vector_ptr,
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
                    braid_BaseVector     vector,
                    _braid_CommHandle  **handle_ptr);

/**
 * Block on the comm handle *handle_ptr* until the MPI operation (send or recv)
 * has completed
 */
braid_Int
_braid_CommWait(braid_Core         core,
               _braid_CommHandle **handle_ptr);

/* uvector.c */

/**
 * Returns an index into the local u-vector for grid *level* at point *index*, 
 * and information on the storage status of the point. If nothing is stored at
 * that point, uindex = -1 and store_flag = -2. If only the shell is stored
 * store_flag = -1, and if the whole u-vector is stored, store_flag = 0.
 */
braid_Int
_braid_UGetIndex(braid_Core   core,
                 braid_Int    level,
                 braid_Int    index,
                 braid_Int   *uindex_ptr,
                 braid_Int   *store_flag_ptr);

/**
 * Returns a reference to the local u-vector on grid *level* at point *index*.
 * If the u-vector is not stored, returns NULL.
 */
braid_Int
_braid_UGetVectorRef(braid_Core        core,
                     braid_Int         level,
                     braid_Int         index,
                     braid_BaseVector *u_ptr);

/**
 * Stores a reference to the local u-vector on grid *level* at point *index*.
 * If the u-vector is not stored, nothing is done.
 */
braid_Int
_braid_USetVectorRef(braid_Core       core,
                     braid_Int        level,
                     braid_Int        index,
                     braid_BaseVector u);

/**
 * Returns a copy of the u-vector on grid *level* at point *index*.  If *index*
 * is my "receive index" (as set by UCommInit(), for example), the u-vector will
 * be received from a neighbor processor.  If the u-vector is not stored, NULL
 * is returned.
 */
braid_Int
_braid_UGetVector(braid_Core        core,
                  braid_Int         level,
                  braid_Int         index,
                  braid_BaseVector *u_ptr);

/**
 * Stores the u-vector on grid *level* at point *index*.  If *index* is my "send
 * index", a send is initiated to a neighbor processor.  If *move* is true, the
 * u-vector is moved into core storage instead of copied.  If the u-vector is
 * not stored, nothing is done.
 */
braid_Int
_braid_USetVector(braid_Core        core,
                  braid_Int         level,
                  braid_Int         index,
                  braid_BaseVector  u,
                  braid_Int         move);

/**
 * Basic communication (from the left, to the right).  Arguments *recv_msg* and
 * *send_msg* are booleans that indicate whether or not to initiate a receive
 * from the left and a send to the right respectively.  Argument *send_now*
 * indicates that the send should be initiated immediately.
 */
braid_Int
_braid_UCommInitBasic(braid_Core  core,
                      braid_Int   level,
                      braid_Int   recv_msg,
                      braid_Int   send_msg,
                      braid_Int   send_now);

/**
 * This routine initiates communication under the assumption that work will be
 * done on all intervals (F or C) on *level*.  It posts a receive for the point
 * to the left of ilower (regardless whether ilower is F or C), and it posts a
 * send of iupper if iupper is a C point.
 */
braid_Int
_braid_UCommInit(braid_Core  core,
                 braid_Int   level);

/**
 * This routine initiates communication under the assumption that work will be
 * done on only F-pt intervals on *level*.  It only posts a receive for the
 * point to the left of ilower if ilower is an F point, and it posts a send of
 * iupper if iupper is a C point.
 */
braid_Int
_braid_UCommInitF(braid_Core  core,
                  braid_Int   level);

/**
 * Finish up communication.  On *level*, wait on both the recv and send handles
 * at this level.
 */
braid_Int
_braid_UCommWait(braid_Core  core,
                 braid_Int   level);

/* step.c */

/**
 * Integrate one time step at time step *index* to time step *index*+1.
 */
braid_Int
_braid_Step(braid_Core        core,
            braid_Int         level,
            braid_Int         index,
            braid_BaseVector  ustop,
            braid_BaseVector  u);

/**
 * Return an initial guess in *ustop_ptr* to use in the step routine for
 * implicit schemes.  The value returned depends on the storage options used.
 * If the return value is NULL, no initial guess is available.
 */
braid_Int
_braid_GetUInit(braid_Core         core,
                braid_Int          level,
                braid_Int          index,
                braid_BaseVector   u,
                braid_BaseVector  *ustop_ptr);

/* residual.c */

/**
 * Compute residual *r*
 */
braid_Int
_braid_Residual(braid_Core       core,
                braid_Int        level,
                braid_Int        index,
                braid_BaseVector ustop,
                braid_BaseVector r);

/**
 * Compute FAS residual = f - residual
 */
braid_Int
_braid_FASResidual(braid_Core       core,
                   braid_Int        level,
                   braid_Int        index,
                   braid_BaseVector ustop,
                   braid_BaseVector r);

/* space.c */

/**
 * Coarsen in space on *level* by calling the user's coarsen function.
 * The vector corresponding to the time step index *f_index* on the fine 
 * grid is coarsened to the time step index *c_index* on the coarse grid.
 * The output goes in *cvector* and the input vector is *fvector*.
 */
braid_Int
_braid_Coarsen(braid_Core        core,
               braid_Int         level,    /* coarse level */
               braid_Int         f_index,  /* fine index */
               braid_Int         c_index,  /* coarse index */
               braid_BaseVector  fvector,
               braid_BaseVector *cvector);

/**
 * Refine in space (basic routine)
 */
braid_Int
_braid_RefineBasic(braid_Core        core,
                   braid_Int         level,    /* fine level */
                   braid_Int         c_index,  /* coarse time index */
                   braid_Real       *f_ta,     /* pointer into fine time array */
                   braid_Real       *c_ta,     /* pointer into coarse time array */
                   braid_BaseVector  cvector,
                   braid_BaseVector *fvector);

/**
 * Refine in space on *level* by calling the user's refine function.
 * The vector corresponding to the time step index *c_index* on the coarse
 * grid is refined to the time step index *f_index* on the fine grid.
 * The output goes in *fvector* and the input vector is *cvector*.

 */
braid_Int
_braid_Refine(braid_Core        core,
              braid_Int         level,    /* fine level */
              braid_Int         f_index,  /* fine index */
              braid_Int         c_index,  /* coarse index */
              braid_BaseVector  cvector,
              braid_BaseVector *fvector);

/** 
 * Call spatial refinement on all local time steps if r_space has been set on
 * the local processor.  Returns refined_ptr == 2 if refinment was completed at
 * any point globally, otherwise returns 0.  This is a helper function for
 * _braid_FRefine().
 */
braid_Int
_braid_FRefineSpace(braid_Core   core,
                    braid_Int   *refined_ptr);

/* grid.c */

/**
 * Create a new grid object *grid_ptr* with level indicator *level*.  Arguments
 * *ilower* and *iupper* correspond to the lower and upper time index values for
 * this processor on this grid.
 */
braid_Int
_braid_GridInit(braid_Core     core,
                braid_Int      level,
                braid_Int      ilower,
                braid_Int      iupper,
                _braid_Grid  **grid_ptr);

/**
 * Destroy the vectors on *grid*
 */
braid_Int
_braid_GridClean(braid_Core    core,
                 _braid_Grid  *grid);

/**
 * Destroy *grid*
 */
braid_Int
_braid_GridDestroy(braid_Core    core,
                   _braid_Grid  *grid);

/* norm.c */

/**
 * Set the residual norm for iteration iter.  If iter < 0, set the rnorm for the
 * last iteration minus |iter|-1.  Also set the initial residual norm.
 */
braid_Int
_braid_SetRNorm(braid_Core  core,
                braid_Int   iter,
                braid_Real  rnorm);

/**
 * Get the residual norm for iteration iter.  If iter < 0, get the rnorm for the
 * last iteration minus |iter|-1.
 */
braid_Int
_braid_GetRNorm(braid_Core  core,
                braid_Int   iter,
                braid_Real *rnorm_ptr);

/**
 * Same as SetRNorm, but sets full residual norm.
 */
braid_Int
_braid_SetFullRNorm(braid_Core  core,
                    braid_Int   iter,
                    braid_Real  rnorm);

/**
 * Same as GetRNorm, but gets full residual norm.
 */
braid_Int
_braid_GetFullRNorm(braid_Core  core,
                    braid_Int   iter,
                    braid_Real *rnorm_ptr);

/**
 * Compute full temporal residual norm with user-provided residual routine. 
 * Output goes in *return_rnorm. 
 */
braid_Int
_braid_ComputeFullRNorm(braid_Core  core,
                        braid_Int   level,
                        braid_Real *return_rnorm);

/**
 * Print out the residual norm for every C-point.
 * Processor 0 gathers all the rnorms and prints them
 * in order through a gatherv operator
 */
braid_Int
_braid_PrintSpatialNorms(braid_Core    core,
                         braid_Real   *rnorms,
                         braid_Int     n);

/* relax.c */

/**
 * Do nu sweeps of F-then-C relaxation on *level*
 */
braid_Int
_braid_FCRelax(braid_Core  core,
               braid_Int   level);

/* restrict.c */

/**
 * F-Relax on *level* and then restrict to *level+1*
 * 
 * The output is set in the braid_Grid in core, so that the restricted vectors
 * *va* and *fa* will be created, representing *level+1* versions of the unknown
 * and rhs vectors.
 * 
 * If the user has set spatial coarsening, then this user-defined routine is
 * also called.
 *
 * If *level==0*, then *rnorm_ptr* will contain the residual norm.
 */
braid_Int
_braid_FRestrict(braid_Core   core,       /**< braid_Core (_braid_Core) struct */   
                 braid_Int    level       /**< restrict from level to level+1 */
                 );

/* interp.c */

/**
 * F-Relax on *level* and interpolate to *level-1*
 *
 * The output is set in the braid_Grid in core, so that the vector *u* on
 * *level* is created by interpolating from *level+1*.
 *
 * If the user has set spatial refinement, then this user-defined routine is
 * also called.
 */
braid_Int
_braid_FInterp(braid_Core  core,   /**< braid_Core (_braid_Core) struct */  
               braid_Int   level   /**< interp from level to level+1 */
               );

/* refine.c */

/**
 * Create a new fine grid (level 0) and corresponding grid hierarchy by refining
 * the current fine grid based on user-provided refinement factors.  Return the
 * boolean *refined_ptr* to indicate whether grid refinement was actually done.
 * To simplify the algorithm, refinement factors are automatically capped to be
 * no greater than the coarsening factor (for level 0).  The grid data is also
 * redistributed to achieve good load balance in the temporal dimension.  If the
 * refinement factor is 1 in each time interval, no refinement is done.
 */
braid_Int
_braid_FRefine(braid_Core   core,
               braid_Int   *refined_ptr);

/* access.c */

/** 
 * Call the user's access function in order to give access to XBraid and the
 * current vector at grid *level and iteration *iter*.  Most commonly, this lets
 * the user write solutions to screen, disk, etc...  The quantity *rnorm*
 * denotes the last computed residual norm, and *done* is a boolean indicating
 * whether XBraid has finished iterating and this is the last Access call.
 */
braid_Int
_braid_FAccess(braid_Core     core,
               braid_Int      level,
               braid_Int      done);

/** 
 * Call user's access function in order to give access to XBraid and the current
 * vector.  Most commonly, this lets the user write *u* to screen, disk, etc...
 * The vector *u* corresponds to time step *index* on *level*.  *status* holds
 * state information about the current XBraid iteration, time value, etc...
 */
braid_Int
_braid_AccessVector(braid_Core         core,
                    braid_AccessStatus status,
                    braid_BaseVector   u);

/**
 * Call user's sync function in order to give access to XBraid and the user's
 * app. This is called once-per-processor at various points in XBraid in
 * order to allow the user to perform any book-keeping operations. *status*
 * provides state information about the current XBraid status and processor.
 */
braid_Int
_braid_Sync(braid_Core       core,
            braid_SyncStatus status);

/* hierarchy.c */

/**
 * Initialize grid hierarchy with *fine_grid* serving as the finest grid.
 * Boolean *refined* indicates whether *fine_grid* was created by refining a
 * coarser grid (in the FRefine() routine), which has implications on how to
 * define CF-intervals.
 */
braid_Int
_braid_InitHierarchy(braid_Core    core,
                     _braid_Grid  *fine_grid,
                     braid_Int     refined);

/**
 * Returns the coarsening factor to use on grid *level*.
 */
braid_Int
_braid_GetCFactor(braid_Core   core,
                  braid_Int    level,
                  braid_Int   *cfactor_ptr);

/**
 * Set initial guess on *level*.  Only C-pts are initialized on level 0,
 * otherwise stored values are initialized based on restricted fine-grid values.
 */
braid_Int
_braid_InitGuess(braid_Core  core,
                 braid_Int   level);

/** 
 * Copy the initialized C-points on the fine grid, to all coarse levels.  For
 * instance, if a point k on level m corresponds to point p on level 0, then
 * they are equivalent after this function.  The only exception is any spatial
 * coarsening the user decides to do.  This function allows XBraid to skip all
 * work on the first down cycle and start in FMG style on the coarsest level.
 * Assumes level 0 C-points are initialized.
 */
braid_Int
_braid_CopyFineToCoarse(braid_Core  core);

/* drive.c */

/**
 * Main loop for MGRIT
 */
braid_Int
_braid_Drive(braid_Core core, 
             braid_Real localtime);

#ifdef __cplusplus
}
#endif

#include "status.h"
#include "base.h"
#include "adjoint.h"

#endif

