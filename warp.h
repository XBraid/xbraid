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
 *  These are all user-written data structures and routines.  There are two
 *  data structures (@ref warp_App and @ref warp_Vector) for the user to define.
 *  And, there are a variety of function interfaces (defined through function pointer
 *  declarations) that the user must implement.
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
 * This defines (roughly) a state vector at a certain time value.  
 * It could also contain any other information related to this vector which is 
 * needed to evolve the vector to the next time value, like mesh information.
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
 * The user must advance the vector *u* from time *tstart* to time *tstop*.
 * The *rfactor_ptr* and *accuracy* inputs are advanced topics.  *rfactor_ptr*
 * allows the user to tell Warp to refine this time interval.
 **/
typedef warp_Int
(*warp_PtFcnPhi)(warp_App      app,                /**< user-defined _warp_App structure */
                 warp_Real     tstart,             /**< time value for *u* */
                 warp_Real     tstop,              /**< time value to evolve *u* towards */
                 warp_Real     accuracy,           /**< advanced option */
                 warp_Vector   u,                  /**< output, vector to evolve */
                 warp_Int     *rfactor_ptr         /**< output, allows user to subdivide this interval for accuracy */
                 );

/**
 * Initializes a vector *u_ptr* at time *t*
 **/
typedef warp_Int
(*warp_PtFcnInit)(warp_App      app,               /**< user-defined _warp_App structure */
                  warp_Real     t,                 /**< time value for *u_ptr* */
                  warp_Vector  *u_ptr              /**< output, newly allocated and initialized vector */
                  );

/**
 * Clone *u* into *v_ptr*
 **/
typedef warp_Int
(*warp_PtFcnClone)(warp_App      app,              /**< user-defined _warp_App structure */
                   warp_Vector   u,                /**< vector to clone */ 
                   warp_Vector  *v_ptr             /**< output, newly allocated and cloned vector */
                   );

/**
 * Free and deallocate *u*
 **/
typedef warp_Int
(*warp_PtFcnFree)(warp_App     app,               /**< user-defined _warp_App structure */
                  warp_Vector  u                  /**< vector to free */
                  );

/**
 * AXPY, *alpha* *x* + *beta* *y* --> *y*
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
 *  *dot_ptr* = <*u*, *v*>
 **/
typedef warp_Int
(*warp_PtFcnDot)(warp_App      app,                /**< user-defined _warp_App structure */
                 warp_Vector   u,                  /**< vector to dot */
                 warp_Vector   v,                  /**< vector to dot */
                 warp_Real    *dot_ptr             /**< output, scalar dot product value */
                 );

/**
 * Write the vector *u* at time *t* to screen, file, etc...  The user decides 
 * what is appropriate.  Notice how you are told the time value of the 
 * vector *u* and even more information in *status*.  This lets you tailor the 
 * output to only certain time values. \n 
 * 
 * If write_level is 2 (see [warp_SetWriteLevel](@ref warp_SetWriteLevel) ), then 
 * *write* is called every Warp iteration and on every Warp level.  In this case, 
 * *status* can be querried using the warp_Get**Status() functions, to determine the 
 * current warp level and iteration.  This allows for even more detailed tracking of the
 * simulation. 
 **/
typedef warp_int
(*warp_PtFcnWrite)(warp_App      app,              /**< user-defined _warp_App structure */
                   warp_Real     t,                /**< time value for *u* */
                   warp_Status   status,           /**< can be querried for info like Warp Iteration */
                   warp_Vector   u                 /**< vector to write */
                   );

/**
 * This routine tells Warp message sizes by computing an upper bound in bytes for 
 * an arbitrary warp_Vector.  This size must be an upper bound for what BufPack and BufUnPack 
 * will assume.
 **/
typedef warp_int
(*warp_PtFcnBufSize)(warp_App   app,               /**< user-defined _warp_App structure */
                     warp_Int  *size_ptr           /**< upper bound on vector size in bytes */
                     );      

/**
 * This allows warp to send messages containing warp_Vectors.  This routine
 * packs a vector _u_ into a _void \*  buffer_ for MPI.
 **/
typedef warp_int
(*warp_PtFcnBufPack)(warp_App      app,            /**< user-defined _warp_App structure */
                     warp_Vector   u,              /**< vector to back into buffer */
                     void         *buffer          /**< output, MPI buffer containing u */
                     );
/**
 * This allows warp to receive messages containing warp_Vectors.  This routine
 * unpacks a _void * buffer_ from MPI into a warp_Vector.
 **/
typedef warp_int
(*warp_PtFcnBufUnpack)(warp_App      app,          /**< user-defined _warp_App structure */
                       void         *buffer,       /**< MPI Buffer to unpack and place in u_ptr */
                       warp_Vector  *u_ptr         /**< output, warp_Vector containing buffer's data */
                       );
/**
 * spatial coarsening (optional).  Allows the user to coarsen
 * when going from a fine time grid to a coarse time grid.
 * This function is called on every vector at each level, thus
 * you can coarsem the entire space time domain.  The action of 
 * this function should match the @ref warp_PtFcnRefine function.
 **/
typedef warp_int
(*warp_PtFcnCoarsen)(warp_App      app,         /**< user-defined _warp_App structure */
                     warp_Real     tstart,      /**< time value for *cu* */                          
                     warp_Real     f_tminus,    /**< time value for *cu* to the left on fine grid */ 
                     warp_Real     f_tplus,     /**< time value for *cu* to the right on fine grid */
                     warp_Real     c_tminus,    /**< time value for *cu* to the left on coarse grid */
                     warp_Real     c_tplus,     /**< time value for *cu* to the right on coarse grid */
                     warp_Vector   fu,          /**< warp_Vector to refine*/                       
                     warp_Vector  *cu_ptr       /**< output, refined vector */    
                     );

/**
 * spatial refinement (optional). Allows the user to refine 
 * when going from a coarse time grid to a fine time grid.  
 * This function is called on every vector at each level, thus
 * you can refine the entire space time domain. The action of 
 * this function should match the @ref warp_PtFcnCoarsen function.
 **/
typedef warp_int
(*warp_PtFcnRefine)(warp_App      app,          /**< user-defined _warp_App structure */
                    warp_Real     tstart,       /**< time value for *cu* */                          
                    warp_Real     f_tminus,     /**< time value for *cu* to the left on fine grid */ 
                    warp_Real     f_tplus,      /**< time value for *cu* to the right on fine grid */
                    warp_Real     c_tminus,     /**< time value for *cu* to the left on coarse grid */
                    warp_Real     c_tplus,      /**< time value for *cu* to the right on coarse grid */
                    warp_Vector   cu,           /**< warp_Vector to refine*/                       
                    warp_Vector  *fu_ptr        /**< output, refined vector */       
                    );
/** @}*/

/*--------------------------------------------------------------------------
 * user interface routines
 *--------------------------------------------------------------------------*/
/** \defgroup userinterface User interface routines
 *  
 *  these are interface routines to initialize and run Warp
 *
 *  @{
 */

struct _warp_core_struct;
/**
 * points to the core structure defined in _warp.h 
 **/
typedef struct _warp_Core_struct *warp_Core;


/**
 * Create a core object with the required initial data.\n
 * This core is used by Warp for internal data structures. 
 * The output is *core_ptr* which points to the newly created 
 * warp_Core structure. 
 **/
warp_int
warp_init(mpi_comm              comm_world,  /**< Global communicator for space and time */
          mpi_comm              comm,        /**< Communicator for temporal dimension*/
          warp_real             tstart,      /**< start time */
          warp_real             tstop,       /**< End time*/
          warp_int              ntime,       /**< Initial number of temporal grid values*/
          warp_app              app,         /**< User-defined _warp_App structure */
          warp_ptFcnPhi         phi,         /**< User time stepping routine to advance a warp_Vector forward one step */
          warp_ptFcnInit        init,        /**< Initialize a warp_Vector on the finest temporal grid*/
          warp_ptFcnClone       clone,       /**< Clone a warp_Vector*/
          warp_ptFcnFree        free,        /**< Free a warp_Vector*/
          warp_ptFcnSum         sum,         /**< Compute vector sum of two warp_Vectors*/
          warp_ptFcnDot         dot,         /**< Compute dot product between two warp_Vectors*/
          warp_ptFcnWrite       write,       /**< Writes a warp_Vector to file, screen */
          warp_ptFcnBufSize     bufsize,     /**< Computes size for MPI buffer for one warp_Vector */
          warp_ptFcnBufPack     bufpack,     /**< Packs MPI buffer to contain one warp_Vector*/
          warp_ptFcnBufUnpack   bufunpack,   /**< Unpacks MPI buffer into a warp_Vector */
          warp_core            *core_ptr     /**< Pointer to warp_Core (_warp_Core) struct*/   
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
 * Set loose stopping tolerance *loose_tol* for spatial solves on grid
 * *level* (level 0 is the finest grid).
 **/
warp_Int
warp_SetLoosexTol(warp_Core  core,        /**< warp_Core (_warp_Core) struct*/
                  warp_Int   level,       /**< level to set *loose_tol* */
                  warp_Real  loose_tol    /**< tolerance to set */
                  );

/**
 * Set tight stopping tolerance *tight_tol* for spatial solves on grid 
 * *level* (level 0 is the finest grid).
 **/
warp_Int
warp_SetTightxTol(warp_Core  core,        /**< warp_Core (_warp_Core) struct*/
                  warp_Int   level,       /**< level to set *tight_tol* */
                  warp_Real  tight_tol    /**< tolerance to set */
                  );

/**
 * Set max number of multigrid levels.
 **/
warp_Int
warp_SetMaxLevels(warp_Core  core,        /**< warp_Core (_warp_Core) struct*/
                  warp_Int   max_levels   /**< maximum levels to allow */
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
 * **Recommended option over relative tolerance**
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
 * Set the number of relaxation sweeps *nrelax* on grid *level*
 * (level 0 is the finest grid).  The default is 1 on all levels.  To change the
 * default factor, use *level = -1*.  One sweep is a CF relaxation sweep.
 **/
warp_Int
warp_SetNRelax(warp_Core  core,           /**< warp_Core (_warp_Core) struct*/           
               warp_Int   level,          /**< *level* to set *nrelax* on */
               warp_Int   nrelax          /**< number of relaxations to do on *level* */
               );

/**
 * Set the coarsening factor *cfactor* on grid *level* (level 0 is
 * the finest grid).  The default factor is 2 on all levels.  To change the
 * default factor, use *level = -1*.
 **/
warp_Int
warp_SetCFactor(warp_Core  core,          /**< warp_Core (_warp_Core) struct*/
                warp_Int   level,         /**< *level* to set coarsening factor on */
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
 * Once called, Warp will use FMG (i.e., F-cycles.
 **/
warp_Int
warp_SetFMG(warp_Core  core               /**< warp_Core (_warp_Core) struct*/
            );

/**
 * Set number of V-cycles to use at each FMG level (standard is 1)
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
 * printed to the Warp print file (@ref warp_SetPrintFile).
 * 
 * - Level 0: no output
 * - Level 1: print typical information like a residual history, 
 *    number of levels in the Warp hierarchy, and so on.
 * - Level 2: level 1 output, plus debug level output.
 * 
 * Default is level 1.
 **/
warp_Int
warp_SetPrintLevel(warp_Core  core,          /**< warp_Core (_warp_Core) struct*/
                   warp_Int   print_level    /**< desired print level */
                   );

/**
 * Set output file for runtime print messages.  Level of printing is 
 * controlled by @ref warp_SetPrintLevel.  Default is stdout.
 **/
warp_Int
warp_SetPrintFile(warp_Core   core,             /**< warp_Core (_warp_Core) struct*/
                  const char *printfile_name    /**< output file for Warp runtime output */
                  );

/**
 * Set write level for warp.  This controls how often the user's
 * write routine is called.
 * 
 * - Level 0:  Never call the user's write routine
 * - Level 1:  Only call the user's write routine after Warp is finished
 * - Level 2:  Call the user's write routine every iteration and on every level.
 *             This is during @ref _warp_FRestrict, during the down-cycle part of a Warp iteration. 
 * 
 * Default is level 1.
 **/
warp_Int
warp_SetWriteLevel(warp_Core  core,          /**< warp_Core (_warp_Core) struct*/
                   warp_Int   write_level    /**< desired write_level */
                   );

/**
 * Split MPI commworld into *comm_x* and *comm_t*, the 
 * spatial and temporal communicators.  The total number of processors
 * will equal Px*Pt, there Px is the number of procs in space, and Pt is 
 * the number of procs in time.
 **/
warp_Int
warp_SplitCommworld(const MPI_Comm  *comm_world,  /**< Global communicator to split */
                          warp_Int  px,           /**< Number of processors parallelizing space for a single time step*/
                          MPI_Comm  *comm_x,      /**< Spatial communicator (written as output) */
                          MPI_Comm  *comm_t       /**< Temporal communicator (written as output) */
                          );

/**
 * Return the residual for the current status object.
 **/
warp_Int
warp_GetStatusResidual(warp_Status  status,     /**< structure containing current simulation info */
                       warp_Real   *rnorm_ptr   /**< output, current residual norm */
                       );

/**
 * Return the iteration for the current status object.
 **/
warp_Int
warp_GetStatusIter(warp_Status  status,         /**< structure containing current simulation info */
                   warp_Int    *iter_ptr        /**< output, current iteration number*/
                   );

/**
 * Return the warp level for the current status object.
 **/
warp_Int
warp_GetStatusLevel(warp_Status  status,        /**< structure containing current simulation info */
                    warp_Int    *level_ptr      /**< output, current level in Warp */
                    );

/**
 * Return whether warp is done for the current status object\n
 * *done_ptr = 1* indicates that this is the last call to Write, because
 * Warp has stopped iterating (either maxiter has been reached, or the tolerance
 * has been met).
 **/
warp_Int
warp_GetStatusDone(warp_Status  status,         /**< structure containing current simulation info */
                   warp_Int    *done_ptr        /**< output,  =1 if warp has finished and this is the final Write, else =0 */
                   );

/**
 * After Drive() finishes, this returns the number of iterations taken.
 **/
warp_Int
warp_GetNumIter(warp_Core  core,          /**< warp_Core (_warp_Core) struct*/
                warp_Int   *niter_ptr     /**< output, holds number of iterations taken */
                );


/**
 * After Drive() finishes, this returns the last measured residual norm.
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

