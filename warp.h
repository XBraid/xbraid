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

/** \file warp.h
 * \brief Define headers for user interface routines.
 *
 * This file contains routines used to allow the user to initialize, run
 * and get and set warp. 
 */

#ifndef warp_HEADER
#define warp_HEADER

#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

typedef int    warp_Int;
typedef double warp_Real;

/*--------------------------------------------------------------------------
 * User-written routines
 *--------------------------------------------------------------------------*/

struct _warp_App_struct;
/**
 * Blah...
 **/
typedef struct _warp_App_struct *warp_App;

struct _warp_Vector_struct;
/**
 * Blah...
 **/
typedef struct _warp_Vector_struct *warp_Vector;

struct _warp_Status_struct;
/**
 * Points to the status structure defined in _warp.h 
 **/
typedef struct _warp_Status_struct *warp_Status;


/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnPhi)(warp_App      app,
                 warp_Real     tstart,
                 warp_Real     tstop,
                 warp_Real     accuracy,
                 warp_Vector   u,
                 warp_Int     *rfactor_ptr);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnInit)(warp_App      app,
                  warp_Real     t,
                  warp_Vector  *u_ptr);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnClone)(warp_App      app,
                   warp_Vector   u,
                   warp_Vector  *v_ptr);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnFree)(warp_App     app,
                  warp_Vector  u);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnSum)(warp_App     app,
                 warp_Real    alpha,
                 warp_Vector  x,
                 warp_Real    beta,
                 warp_Vector  y);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnDot)(warp_App      app,
                 warp_Vector   u,
                 warp_Vector   v,
                 warp_Real    *dot_ptr);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnWrite)(warp_App      app,
                   warp_Real     t,
                   warp_Status   status,
                   warp_Vector   u);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnBufSize)(warp_App   app,
                     warp_Int  *size_ptr);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnBufPack)(warp_App      app,
                     warp_Vector   u,
                     void         *buffer);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnBufUnpack)(warp_App      app,
                       void         *buffer,
                       warp_Vector  *u_ptr);

/**
 * (Optional) Blah...
 **/
typedef warp_Int
(*warp_PtFcnCoarsen)(warp_App      app,
                     warp_Real     tstart,
                     warp_Real     f_tminus,
                     warp_Real     f_tplus,
                     warp_Real     c_tminus,
                     warp_Real     c_tplus,
                     warp_Vector   fu,
                     warp_Vector  *cu_ptr);

/**
 * (Optional) Blah...
 **/
typedef warp_Int
(*warp_PtFcnRefine)(warp_App      app,
                    warp_Real     tstart,
                    warp_Real     f_tminus,
                    warp_Real     f_tplus,
                    warp_Real     c_tminus,
                    warp_Real     c_tplus,
                    warp_Vector   cu,
                    warp_Vector  *fu_ptr);


/*--------------------------------------------------------------------------
 * User interface routines
 *--------------------------------------------------------------------------*/

struct _warp_Core_struct;
/**
 * Points to the core structure defined in _warp.h 
 **/
typedef struct _warp_Core_struct *warp_Core;


/**
 * Create a core object with the required initial data.\n
 * Output:
 * - *core_ptr* will point to the newly created warp_Core structure 
 **/
warp_Int
warp_Init(MPI_Comm              comm_world,  /**< Global communicator for space and time */
          MPI_Comm              comm,        /**< Communicator for temporal dimension*/
          warp_Real             tstart,      /**< start time */
          warp_Real             tstop,       /**< End time*/
          warp_Int              ntime,       /**< Initial number of temporal grid values*/
          warp_App              app,         /**< User defined structure to hold *state* information */
          warp_PtFcnPhi         phi,         /**< User time stepping routine to advance state one time value*/
          warp_PtFcnInit        init,        /**< Initialize a warp_Vector function on finest temporal grid*/
          warp_PtFcnClone       clone,       /**< Clone a warp_Vector*/
          warp_PtFcnFree        free,        /**< Free a temporal state warp_Vector*/
          warp_PtFcnSum         sum,         /**< Compute vector sum of two temporal states*/
          warp_PtFcnDot         dot,         /**< Compute dot product between two temporal states*/
          warp_PtFcnWrite       write,       /**< *Writes* (file, screen..) upon completion. */
          warp_PtFcnBufSize     bufsize,     /**< Computes size for MPI buffer for one */
          warp_PtFcnBufPack     bufpack,     /**< Packs MPI buffer to contain one temporal state*/
          warp_PtFcnBufUnpack   bufunpack,   /**< Unpacks MPI buffer containing one temporal state*/
          warp_Core            *core_ptr     /**< Pointer to warp_Core (_warp_Core) struct*/   
          );

/**
 * Integrate in time.
 **/
warp_Int
warp_Drive(warp_Core  core);

/**
 * Destroy core.
 **/
warp_Int
warp_Destroy(warp_Core  core);

/**
 * Print statistics.
 **/
warp_Int
warp_PrintStats(warp_Core  core);

/**
 * Set loose stopping tolerance for spatial solves on grid level
 * *level* (level 0 is the finest grid).
 **/
warp_Int
warp_SetLoosexTol(warp_Core  core,
                  warp_Int   level,
                  warp_Real  loose_tol);

/**
 * Set tight stopping tolerance for spatial solves on grid level
 * *level* (level 0 is the finest grid).
 **/
warp_Int
warp_SetTightxTol(warp_Core  core,
                  warp_Int   level,
                  warp_Real  tight_tol);

/**
 * Set max number of multigrid levels.
 **/
warp_Int
warp_SetMaxLevels(warp_Core  core,
                  warp_Int   max_levels);

/**
 * Set absolute stopping tolerance.
 **/
warp_Int
warp_SetAbsTol(warp_Core  core,
               warp_Real  atol);

/**
 * Set absolute stopping tolerance.
 **/
warp_Int
warp_SetRelTol(warp_Core  core,
               warp_Real  rtol);

/**
 * Set the number of relaxation sweeps *nrelax* on grid level *level*
 * (level 0 is the finest grid).  The default is 1 on all levels.  To change the
 * default factor, use *level* = -1.
 **/
warp_Int
warp_SetNRelax(warp_Core  core,
               warp_Int   level,
               warp_Int   nrelax);

/**
 * Set the coarsening factor *cfactor* on grid level *level* (level 0 is
 * the finest grid).  The default factor is 2 on all levels.  To change the
 * default factor, use *level* = -1.
 **/
warp_Int
warp_SetCFactor(warp_Core  core,
                warp_Int   level,
                warp_Int   cfactor);

/**
 * Set max number of multigrid iterations.
 **/
warp_Int
warp_SetMaxIter(warp_Core  core,
                warp_Int   max_iter);

/**
 * Use FMG cycling.
 **/
warp_Int
warp_SetFMG(warp_Core  core);

/**
 * Set spatial coarsening routine with user-defined routine.
 * Default is no spatial refinment or coarsening.
 **/
warp_Int
warp_SetSpatialCoarsen(warp_Core  core, 
                       warp_PtFcnCoarsen coarsen);

/**
 * Set spatial refinement routine with user-defined routine.
 * Default is no spatial refinment or coarsening.
 **/
warp_Int
warp_SetSpatialRefine(warp_Core  core,
                      warp_PtFcnRefine refine);

/**
 * Set print level for warp.  This controls how much information is 
 * printed to standard out.
 * 
 * Level 0: no output
 * Level 1: print typical information like a residual history, 
 *    number of levels in the Warp hierarchy, and so on.
 * Level 2: level 1 output, plus debug level output.
 * 
 * Default is level 1.
 **/
warp_Int
warp_SetPrintLevel(warp_Core  core,
                   warp_Int   print_level);
/**
 * Set write level for warp.  This controls how often the user's
 * write routine is called.
 * 
 * Level 0:  Never call the user's write routine
 * Level 1:  Only call the user's write routine after Warp is finished
 * Level 2:  Call the user's write routine every iteration in _warp_FRestrict(),
 *   which is during the down-cycle part of a Warp iteration 
 * 
 * Default is level 1.
 **/
warp_Int
warp_SetWriteLevel(warp_Core  core,
                   warp_Int   write_level);

/**
 * Split MPI commworld into comm_x and comm_t, the 
 * spatial and temporal communicators
 **/
warp_Int
warp_SplitCommworld(warp_Core core,    /**< Warp Core*/
                    warp_Int  px,      /**< Number of processors parallelizing space for a single time step*/
                    MPI_Comm  comm_x,  /**< Spatial communicator (written as output) */
                    MPI_Comm  comm_t   /**< Temporal communicator (written as output) */
                    );

/**
 * Return the residual for the current status object
 **/
warp_Int
warp_GetStatusResidual(warp_Status  status,
                       warp_Real   *rnorm_ptr);

/**
 * Return the iteration for the current status object
 **/
warp_Int
warp_GetStatusIter(warp_Status  status,
                   warp_Int    *iter_ptr);

/**
 * Return the warp level for the current status object
 **/
warp_Int
warp_GetStatusLevel(warp_Status  status,
                    warp_Int    *level_ptr);

/**
 * Return whether warp is done for the current status object
 **/
warp_Int
warp_GetStatusDone(warp_Status  status,
                   warp_Int    *done_ptr);


#ifdef __cplusplus
}
#endif

#endif

