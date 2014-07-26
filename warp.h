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
 * Define basic types
 *--------------------------------------------------------------------------*/
/**
 * Defines integer type
 **/
typedef int    warp_Int;

/**
 * Defines floating point type
 **/
typedef double warp_Real;

/*--------------------------------------------------------------------------
 * User-written routines
 *--------------------------------------------------------------------------*/
/** \defgroup userwritten User-written routines
 *  
 *  These are all user-written data structures and routines 
 *
 *  @{
 */
struct _warp_App_struct;
/**
 * This holds a wide variety of information and is ``global`` in that it
 * is passed to every function.  This structure holds everything that the user 
 * will need to carry out a simulation.  For a simple example, this could just
 * hold the global MPI communicator and a few values describing the temporal domain.
 **/
typedef struct _warp_App_struct *warp_App;

struct _warp_Vector_struct;
/**
 * This defines something like a state vector at a certain time value.  It 
 * could contain any other information related to this vector needed to evolve
 * the vector to the next time value, like mesh information.
 **/
typedef struct _warp_Vector_struct *warp_Vector;

struct _warp_Status_struct;
/**
 * Points to the status structure defined in _warp.h 
 * This is NOT a user-defined structure.
 **/
typedef struct _warp_Status_struct *warp_Status;


/**
 * Defines the central time stepping function that the user must write.
 * The user must advance the vector u from time tstart to time tstop.
 * Here advancing the solution just involves the scalar \f$ \lambda \f$.  
 * The rfactor_ptr and accuracy inputs are advanced topics.  rfactor_ptr
 * allows the user to tell Warp to refine this time interval.
 **/
typedef warp_Int
(*warp_PtFcnPhi)(warp_App      app,                /**< user-defined _warp_App structure */
                 warp_Real     tstart,             /**< time value for u */
                 warp_Real     tstop,              /**< time value to evolve u towards */
                 warp_Real     accuracy,           /**< advanced option */
                 warp_Vector   u,                  /**< output, vector to evolve */
                 warp_Int     *rfactor_ptr         /**< output, allows user to subdivide this interval for accuracy */
                 );

/**
 * Initializes a vector at time t
 **/
typedef warp_Int
(*warp_PtFcnInit)(warp_App      app,               /**< user-defined _warp_App structure */
                  warp_Real     t,                 /**< time value for u */
                  warp_Vector  *u_ptr              /**< output, newly allocated and initialized vector */
                  );

/**
 * Clone u into v_ptr
 **/
typedef warp_Int
(*warp_PtFcnClone)(warp_App      app,              /**< user-defined _warp_App structure */
                   warp_Vector   u,                /**< vector to clone */ 
                   warp_Vector  *v_ptr             /**< output, newly allocated and cloned vector */
                   );

/**
 * Free and deallocate u
 **/
typedef warp_Int
(*warp_PtFcnFree)(warp_App     app,               /**< user-defined _warp_App structure */
                  warp_Vector  u                  /**< vector to free */
                  );

/**
 * AXPY, alpha*x + beta*y --> y
 **/
typedef warp_Int
(*warp_PtFcnSum)(warp_App     app,                /**< user-defined _warp_App structure */
                 warp_Real    alpha,              /**< scalar for AXPY */
                 warp_Vector  x,                  /**< vector for AXPY */
                 warp_Real    beta,               /**< scalar for AXPY */
                 warp_Vector  y                   /**< output and vector for AXPY */
                 );

/**
 *  Carry out a dot product
 *  *dot_ptr = <u, v>
 **/
typedef warp_Int
(*warp_PtFcnDot)(warp_App      app,                /**< user-defined _warp_App structure */
                 warp_Vector   u,                  /**< vector to dot */
                 warp_Vector   v,                  /**< vector to dot */
                 warp_Real    *dot_ptr             /**< output, scalar dot product value */
                 );

/**
 * Write the vector u at time t.  The user decides whether to write to
 * file, screen or do nothing.  Only solution values at certain times or
 * Warp iterations may be desired.
 **/
typedef warp_Int
(*warp_PtFcnWrite)(warp_App      app,              /**< user-defined _warp_App structure */
                   warp_Real     t,                /**< time value for u */
                   warp_Status   status,           /**< can be querried for info like Warp Iteration */
                   warp_Vector   u                 /**< vector to write */
                   );

/**
 * Computes an upper bound for the MPI Buffer size in bytes for an arbitrary warp_Vector
 **/
typedef warp_Int
(*warp_PtFcnBufSize)(warp_App   app,               /**< user-defined _warp_App structure */
                     warp_Int  *size_ptr           /**< upper bound on vector size in bytes */
                     );      

/**
 * Packs a vector u into a void * buffer for MPI
 **/
typedef warp_Int
(*warp_PtFcnBufPack)(warp_App      app,            /**< user-defined _warp_App structure */
                     warp_Vector   u,              /**< vector to back into buffer */
                     void         *buffer          /**< output, MPI buffer containing u */
                     );
/**
 * Unpack that buffer from void * to warp_Vector
 **/
typedef warp_Int
(*warp_PtFcnBufUnpack)(warp_App      app,          /**< user-defined _warp_App structure */
                       void         *buffer,       /**< MPI Buffer to unpack and place in u_ptr */
                       warp_Vector  *u_ptr         /**< output, warp_Vector containing buffer's data */
                       );
/**
 * Spatial coarsening (optional).  Allows the user to coarsen
 * when going from a fine time grid to a coarse time grid.
 * This function is called on every vector at each level, thus
 * you can coarsem the entire space time domain.
 **/
typedef warp_Int
(*warp_PtFcnCoarsen)(warp_App      app,         /**< user-defined _warp_App structure */
                     warp_Real     tstart,      /**< time value for cu */                          
                     warp_Real     f_tminus,    /**< time value for cu to the left on fine grid */ 
                     warp_Real     f_tplus,     /**< time value for cu to the right on fine grid */
                     warp_Real     c_tminus,    /**< time value for cu to the left on coarse grid */
                     warp_Real     c_tplus,     /**< time value for cu to the right on coarse grid */
                     warp_Vector   fu,          /**< warp_Vector to refine*/                       
                     warp_Vector  *cu_ptr       /**< output, refined vector */    
                     );

/**
 * Spatial refinement (optional). Allows the user to refine 
 * when going from a coarse time grid to a fine time grid.  
 * This function is called on every vector at each level, thus
 * you can refine the entire space time domain.
 **/
typedef warp_Int
(*warp_PtFcnRefine)(warp_App      app,          /**< user-defined _warp_App structure */
                    warp_Real     tstart,       /**< time value for cu */                          
                    warp_Real     f_tminus,     /**< time value for cu to the left on fine grid */ 
                    warp_Real     f_tplus,      /**< time value for cu to the right on fine grid */
                    warp_Real     c_tminus,     /**< time value for cu to the left on coarse grid */
                    warp_Real     c_tplus,      /**< time value for cu to the right on coarse grid */
                    warp_Vector   cu,           /**< warp_Vector to refine*/                       
                    warp_Vector  *fu_ptr        /**< output, refined vector */       
                    );
/** @}*/

/*--------------------------------------------------------------------------
 * User interface routines
 *--------------------------------------------------------------------------*/
/** \defgroup userinterface User interface routines
 *  
 *  These are interface routines to initialize and run Warp
 *
 *  @{
 */

struct _warp_Core_struct;
/**
 * Points to the core structure defined in _warp.h 
 **/
typedef struct _warp_Core_struct *warp_Core;


/**
 * Create a core object with the required initial data.\n
 * The output is the core_ptr pointing to the newly created warp_Core structure. 
 **/
warp_Int
warp_Init(MPI_Comm              comm_world,  /**< Global communicator for space and time */
          MPI_Comm              comm,        /**< Communicator for temporal dimension*/
          warp_Real             tstart,      /**< start time */
          warp_Real             tstop,       /**< End time*/
          warp_Int              ntime,       /**< Initial number of temporal grid values*/
          warp_App              app,         /**< User-defined _warp_App structure */
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
 * Carry out a simulation with Warp.  Integrate in time.
 **/
warp_Int
warp_Drive(warp_Core  core                /**< warp_Core (_warp_Core) struct*/
          );

/**
 * Clean up and destroy core.
 **/
warp_Int
warp_Destroy(warp_Core  core              /**< warp_Core (_warp_Core) struct*/
            );

/**
 * Print statistics after a Warp run.
 **/
warp_Int
warp_PrintStats(warp_Core  core           /**< warp_Core (_warp_Core) struct*/
               );

/**
 * Set loose stopping tolerance for spatial solves on grid level
 * *level* (level 0 is the finest grid).
 **/
warp_Int
warp_SetLoosexTol(warp_Core  core,        /**< warp_Core (_warp_Core) struct*/
                  warp_Int   level,       /**< level to set loose_tol */
                  warp_Real  loose_tol    /**< tolerance to set */
                  );

/**
 * Set tight stopping tolerance for spatial solves on grid level
 * *level* (level 0 is the finest grid).
 **/
warp_Int
warp_SetTightxTol(warp_Core  core,        /**< warp_Core (_warp_Core) struct*/
                  warp_Int   level,       /**< level to set tight_tol */
                  warp_Real  tight_tol    /**< tolerance to set */
                  );

/**
 * Set max number of multigrid levels.
 **/
warp_Int
warp_SetMaxLevels(warp_Core  core,        /**< warp_Core (_warp_Core) struct*/
                  warp_Int   max_levels   /**< maximum levels */
                  );

/**
 * Set max allowed coarse grid size (in terms of C-points) 
 **/
warp_Int
warp_SetMaxCoarse(warp_Core  core,        /**< warp_Core (_warp_Core) struct*/
                  warp_Int   max_coarse   /** maximum coarse grid size */
                  );

/**
 * Set absolute stopping tolerance.\n
 * **Recommended Option**
 **/
warp_Int
warp_SetAbsTol(warp_Core  core,           /**< warp_Core (_warp_Core) struct*/
               warp_Real  atol            /**< absolute stopping tolerance */
               );

/**
 * Set relative stopping tolerance, relative to the initial residual.  Be 
 * careful.  If your initial guess is all zero, then the initial residual
 * may only be nonzero over one or two time values, and this will skew the
 * relative tolerance.  Absolute tolerances are recommended.
 **/
warp_Int
warp_SetRelTol(warp_Core  core,           /**< warp_Core (_warp_Core) struct*/
               warp_Real  rtol            /**< relative stopping tolerance */
               );

/**
 * Set the number of relaxation sweeps *nrelax* on grid level *level*
 * (level 0 is the finest grid).  The default is 1 on all levels.  To change the
 * default factor, use *level* = -1.  One sweep is a CF relaxation sweep.
 **/
warp_Int
warp_SetNRelax(warp_Core  core,           /**< warp_Core (_warp_Core) struct*/           
               warp_Int   level,          /**< level to set nrelax on */
               warp_Int   nrelax          /**< number of relaxations to do on level */
               );

/**
 * Set the coarsening factor *cfactor* on grid level *level* (level 0 is
 * the finest grid).  The default factor is 2 on all levels.  To change the
 * default factor, use *level* = -1.
 **/
warp_Int
warp_SetCFactor(warp_Core  core,          /**< warp_Core (_warp_Core) struct*/
                warp_Int   level,         /**< level to set coarsening factor on */
                warp_Int   cfactor        /**< desired coarsening factor */
                );

/**
 * Set max number of multigrid iterations.
 **/
warp_Int
warp_SetMaxIter(warp_Core  core,          /**< warp_Core (_warp_Core) struct*/
                warp_Int   max_iter       /**< maximum iterations to allow */
                );

/**
 * Use FMG cycling.
 **/
warp_Int
warp_SetFMG(warp_Core  core               /**< warp_Core (_warp_Core) struct*/
            );

/**
 * Set number of V cycles to use at each FMG level (standard is 1)
 **/
warp_Int
warp_SetNFMGVcyc(warp_Core  core,         /**< warp_Core (_warp_Core) struct*/
                 warp_Int   nfmg_Vcyc     /**< number of V-cycles to do each FMG level */
                 );


/**
 * Set spatial coarsening routine with user-defined routine.
 * Default is no spatial refinment or coarsening.
 **/
warp_Int
warp_SetSpatialCoarsen(warp_Core  core,            /**< warp_Core (_warp_Core) struct*/ 
                       warp_PtFcnCoarsen coarsen   /**< function pointer to spatial coarsening routine */
                       );

/**
 * Set spatial refinement routine with user-defined routine.
 * Default is no spatial refinment or coarsening.
 **/
warp_Int
warp_SetSpatialRefine(warp_Core  core,             /**< warp_Core (_warp_Core) struct*/
                      warp_PtFcnRefine refine      /**< function pointer to spatial refinement routine */
                      );

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
warp_SetPrintLevel(warp_Core  core,          /**< warp_Core (_warp_Core) struct*/
                   warp_Int   print_level    /**< desired print level */
                   );

/**
 * Set output file for runtime print messages.  Level of printing is 
 * controlled by warp_SetPrintLevel.  Default is stdout.
 **/
warp_Int
warp_SetPrintFile(warp_Core   core,             /**< warp_Core (_warp_Core) struct*/
                  const char *printfile_name    /**< output file for Warp runtime output */
                  );

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
warp_SetWriteLevel(warp_Core  core,          /**< warp_Core (_warp_Core) struct*/
                   warp_Int   write_level    /**< desired write_level */
                   );

/**
 * Split MPI commworld into comm_x and comm_t, the 
 * spatial and temporal communicators
 **/
warp_Int
warp_SplitCommworld(const MPI_Comm  *comm_world,  /**< Global communicator to split */
                          warp_Int  px,           /**< Number of processors parallelizing space for a single time step*/
                          MPI_Comm  *comm_x,      /**< Spatial communicator (written as output) */
                          MPI_Comm  *comm_t       /**< Temporal communicator (written as output) */
                          );

/**
 * Return the residual for the current status object
 **/
warp_Int
warp_GetStatusResidual(warp_Status  status,     /**< structure containing current simulation info */
                       warp_Real   *rnorm_ptr   /**< output, current residual norm */
                       );

/**
 * Return the iteration for the current status object
 **/
warp_Int
warp_GetStatusIter(warp_Status  status,         /**< structure containing current simulation info */
                   warp_Int    *iter_ptr        /**< output, current iteration number*/
                   );

/**
 * Return the warp level for the current status object
 **/
warp_Int
warp_GetStatusLevel(warp_Status  status,        /**< structure containing current simulation info */
                    warp_Int    *level_ptr      /**< output, current level in Warp */
                    );

/**
 * Return whether warp is done for the current status object
 **/
warp_Int
warp_GetStatusDone(warp_Status  status,         /**< structure containing current simulation info */
                   warp_Int    *done_ptr        /**< output,  =1 if warp has finished and this is the final Write, else =0 */
                   );

/**
 * After Drive() finishes, this returns the number of iterations taken
 **/
warp_Int
warp_GetNumIter(warp_Core  core,          /**< warp_Core (_warp_Core) struct*/
                warp_Int   *niter_ptr     /**< output, holds number of iterations taken */
                );


/**
 * After Drive() finishes, this returns the last measured residual norm
 **/
warp_Int
warp_GetRNorm(warp_Core  core,            /**< warp_Core (_warp_Core) struct*/
              warp_Real  *rnorm_ptr       /**< output, holds final residual norm */
              );

/** @}*/

#ifdef __cplusplus
}
#endif

#endif

