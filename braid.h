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
 
/** \file braid.h
 * \brief Define headers for user interface routines.
 *
 * This file contains routines used to allow the user to initialize, run
 * and get and set a XBraid solver. 
 */

#ifndef braid_HEADER
#define braid_HEADER

#include "mpi.h"
#include "braid_defs.h"
#include "braid_status.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Turn on Fortran 90 interface options manually */
#define braid_FMANGLE 1
#define braid_Fortran_SpatialCoarsen 0
#define braid_Fortran_Residual 0

/** Value used to represent an invalid residual norm */
#define braid_INVALID_RNORM -1

/*--------------------------------------------------------------------------
 * Error codes
 *--------------------------------------------------------------------------*/

#define braid_ERROR_GENERIC         1   /* generic error */
#define braid_ERROR_MEMORY          2   /* unable to allocate memory */
#define braid_ERROR_ARG             4   /* argument error */
/* bits 4-8 are reserved for the index of the argument error */

/*--------------------------------------------------------------------------
 * User-written routines
 *--------------------------------------------------------------------------*/
/** \defgroup userwritten User-written routines
 *  
 *  These are all user-written data structures and routines.  There are two
 *  data structures (@ref braid_App and @ref braid_Vector) for the user to define.
 *  And, there are a variety of function interfaces (defined through function pointer
 *  declarations) that the user must implement.
 *
 *  @{
 */
struct _braid_App_struct;
/**
 * This holds a wide variety of information and is ``global`` in that it
 * is passed to every function.  This structure holds everything that the user 
 * will need to carry out a simulation.  For a simple example, this could just
 * hold the global MPI communicator and a few values describing the temporal domain.
 **/
typedef struct _braid_App_struct *braid_App;

struct _braid_Vector_struct;
/**
 * This defines (roughly) a state vector at a certain time value.  
 * It could also contain any other information related to this vector which is 
 * needed to evolve the vector to the next time value, like mesh information.
 **/
typedef struct _braid_Vector_struct *braid_Vector;

/**
 * Defines the central time stepping function that the user must write.
 *
 * The user must advance the vector *u* from time *tstart* to *tstop*.  The time
 * step is taken assuming the right-hand-side vector *fstop* at time *tstop*.
 * The vector *ustop* may be the same vector as *u* (in the case where not all
 * unknowns are stored).  The vector *fstop* is set to NULL to indicate a zero
 * right-hand-side.
 *
 * Query the status structure with *braid_StepStatusGetTstart(status, &tstart)*
 * and *braid_StepStatusGetTstop(status, &tstop)* to get *tstart* and *tstop*.
 * The status structure also allows for steering.  For example,
 * *braid_StepStatusSetRFactor(...)* allows for setting rfactor, which tells
 * XBraid to refine this time interval.
 **/
typedef braid_Int
(*braid_PtFcnStep)(braid_App        app,    /**< user-defined _braid_App structure */
                   braid_Vector     ustop,  /**< input, u vector at *tstop* */
                   braid_Vector     fstop,  /**< input, right-hand-side at *tstop* */
                   braid_Vector     u     , /**< input/output, initially u vector at *tstart*, upon exit, u vector at *tstop* */
                   braid_StepStatus status  /**< query this struct for info about u (e.g., tstart and tstop), allows for steering (e.g., set rfactor) */ 
                   );

/**
 * Initializes a vector *u_ptr* at time *t*
 **/
typedef braid_Int
(*braid_PtFcnInit)(braid_App      app,           /**< user-defined _braid_App structure */
                   braid_Real     t,             /**< time value for *u_ptr* */
                   braid_Vector  *u_ptr          /**< output, newly allocated and initialized vector */
                   );

/**
 * Clone *u* into *v_ptr*
 **/
typedef braid_Int
(*braid_PtFcnClone)(braid_App      app,          /**< user-defined _braid_App structure */
                    braid_Vector   u,            /**< vector to clone */ 
                    braid_Vector  *v_ptr         /**< output, newly allocated and cloned vector */
                    );

/**
 * Free and deallocate *u*
 **/
typedef braid_Int
(*braid_PtFcnFree)(braid_App     app,            /**< user-defined _braid_App structure */
                   braid_Vector  u               /**< vector to free */
                   );

/**
 * AXPY, *alpha* *x* + *beta* *y* --> *y*
 **/ 
typedef braid_Int
(*braid_PtFcnSum)(braid_App     app,             /**< user-defined _braid_App structure */
                  braid_Real    alpha,           /**< scalar for AXPY */
                  braid_Vector  x,               /**< vector for AXPY */
                  braid_Real    beta,            /**< scalar for AXPY */
                  braid_Vector  y                /**< output and vector for AXPY */
                  );

/** 
 * Carry out a spatial norm by taking the norm of a braid_Vector 
 * *norm_ptr* = || *u* || 
 * A common choice is the standard Eucliden norm, but many other choices are
 * possible, such as an L2-norm based on a finite element space.  See
 * [braid_SetTemporalNorm](@ref braid_SetTemporalNorm) for information on how the
 * spatial norm is combined over time for a global space-time residual norm.
 * This global norm then controls halting. 
 **/
typedef braid_Int
(*braid_PtFcnSpatialNorm)(braid_App      app,      /**< user-defined _braid_App structure */
                          braid_Vector   u,        /**< vector to norm */
                          braid_Real    *norm_ptr  /**< output, norm of braid_Vector (this is a spatial norm) */ 
                          );

/**
 * Gives user access to XBraid and to the current vector *u* at time *t*.  Most
 * commonly, this lets the user write the vector to screen, file, etc...  The
 * user decides what is appropriate.  Note how you are told the time value *t*
 * of the vector *u* and other information in *status*.  This lets you tailor
 * the output, e.g., for only certain time values at certain XBraid iterations.
 * Querrying status for such information is done through
 * _braid_AccessStatusGet**(..)_ routines.
 * 
 * The frequency of XBraid's calls to *access* is controlled through
 * [braid_SetAccessLevel](@ref braid_SetAccessLevel).  For instance, if
 * access_level is set to 2, then *access* is called every XBraid iteration and
 * on every XBraid level.  In this case, querrying *status* to determine the
 * current XBraid level and iteration will be useful. This scenario allows for
 * even more detailed tracking of the simulation.
 *
 * Eventually, access will be broadened to allow the user to steer XBraid.
 **/
typedef braid_Int
(*braid_PtFcnAccess)(braid_App           app,              /**< user-defined _braid_App structure */
                     braid_Vector        u,                /**< vector to be accessed */
                     braid_AccessStatus  status            /**< can be querried for info like the current XBraid Iteration */
                     );

/**
 * This routine tells XBraid message sizes by computing an upper bound in bytes
 * for an arbitrary braid_Vector.  This size must be an upper bound for what
 * BufPack and BufUnPack will assume.
 **/
typedef braid_Int
(*braid_PtFcnBufSize)(braid_App   app,               /**< user-defined _braid_App structure */
                      braid_Int  *size_ptr           /**< upper bound on vector size in bytes */
                      );      

/**
 * This allows XBraid to send messages containing braid_Vectors.  This routine
 * packs a vector _u_ into a _void \*  buffer_ for MPI.
 **/
typedef braid_Int
(*braid_PtFcnBufPack)(braid_App      app,            /**< user-defined _braid_App structure */
                      braid_Vector   u,              /**< vector to back into buffer */
                      void          *buffer,         /**< output, MPI buffer containing u */
                      braid_Int     *size_ptr        /**< output, number of bytes packed, must be less than or equal to value returned by BufSize */
                      );

/**
 * This allows XBraid to receive messages containing braid_Vectors.  This routine
 * unpacks a _void * buffer_ from MPI into a braid_Vector.
 **/
typedef braid_Int
(*braid_PtFcnBufUnpack)(braid_App      app,          /**< user-defined _braid_App structure */
                        void          *buffer,       /**< MPI Buffer to unpack and place in u_ptr */
                        braid_Vector  *u_ptr         /**< output, braid_Vector containing buffer's data */
                        );

/**
 * This function (optional) computes the residual *r* at time *tstop*.  On
 * input, *r* holds the value of *u* at *tstart*, and *ustop* is the value of
 * *u* at *tstop*.  If used, set with @ref braid_SetResidual.
 *
 * Query the status structure with *braid_StepStatusGetTstart(status, &tstart)*
 * and *braid_StepStatusGetTstop(status, &tstop)* to get *tstart* and *tstop*.
 **/
typedef braid_Int
(*braid_PtFcnResidual)(braid_App        app,    /**< user-defined _braid_App structure */
                       braid_Vector     ustop,  /**< input, u vector at *tstop* */
                       braid_Vector     r     , /**< output, residual at *tstop* */
                       braid_StepStatus status  /**< query this struct for info about u (e.g., tstart and tstop) */ 
                       );

/**
 * Spatial coarsening (optional).  Allows the user to coarsen
 * when going from a fine time grid to a coarse time grid.
 * This function is called on every vector at each level, thus
 * you can coarsen the entire space time domain.  The action of 
 * this function should match the @ref braid_PtFcnSRefine function.
 *
 * The user should query the status structure at run time with 
 * _braid_CoarsenRefGet**()_ calls in order to determine how to coarsen.  
 * For instance, status tells you what the current time value is, and 
 * what the time step sizes on the fine and coarse levels are.
 **/
typedef braid_Int
(*braid_PtFcnSCoarsen)(braid_App               app,    /**< user-defined _braid_App structure */
                       braid_Vector            fu,     /**< braid_Vector to refine*/
                       braid_Vector           *cu_ptr, /**< output, refined vector */   
                       braid_CoarsenRefStatus  status  /**< query this struct for info about fu and cu (e.g., where in time fu and cu are)  */ 
                       );

/**
 * Spatial refinement (optional). Allows the user to refine 
 * when going from a coarse time grid to a fine time grid.  
 * This function is called on every vector at each level, thus
 * you can refine the entire space time domain. The action of 
 * this function should match the @ref braid_PtFcnSCoarsen function.
 *
 * The user should query the status structure at run time with 
 * _braid_CoarsenRefGet**()_ calls in order to determine how to coarsen.  
 * For instance, status tells you what the current time value is, and 
 * what the time step sizes on the fine and coarse levels are.
 **/
typedef braid_Int
(*braid_PtFcnSRefine)(braid_App               app,    /**< user-defined _braid_App structure */
                      braid_Vector            cu,     /**< braid_Vector to refine*/
                      braid_Vector           *fu_ptr, /**< output, refined vector */       
                      braid_CoarsenRefStatus  status  /**< query this struct for info about fu and cu (e.g., where in time fu and cu are)  */ 
                      );
/** @}*/

/*--------------------------------------------------------------------------
 * User Interface Routines
 *--------------------------------------------------------------------------*/
/** \defgroup userinterface User interface routines
 *  
 *  These are all the user interface routines.
 *  @{
 */
/** @}*/

/** \defgroup generalinterface General Interface routines
 *  \ingroup userinterface
 *
 *  These are general interface routines, e.g., routines to initialize and 
 *  run a XBraid solver, or to split a communicator into spatial and temporal
 *  components.
 *
 *  @{
 */

struct _braid_Core_struct;
/**
 * points to the core structure defined in _braid.h 
 **/
typedef struct _braid_Core_struct *braid_Core;


/**
 * Create a core object with the required initial data.
 *
 * This core is used by XBraid for internal data structures. 
 * The output is *core_ptr* which points to the newly created 
 * braid_Core structure. 
 **/
braid_Int
braid_Init(MPI_Comm               comm_world,  /**< Global communicator for space and time */
           MPI_Comm               comm,        /**< Communicator for temporal dimension*/
           braid_Real             tstart,      /**< start time */
           braid_Real             tstop,       /**< End time*/
           braid_Int              ntime,       /**< Initial number of temporal grid values*/
           braid_App              app,         /**< User-defined _braid_App structure */
           braid_PtFcnStep        step,        /**< User time stepping routine to advance a braid_Vector forward one step */
           braid_PtFcnInit        init,        /**< Initialize a braid_Vector on the finest temporal grid*/
           braid_PtFcnClone       clone,       /**< Clone a braid_Vector*/
           braid_PtFcnFree        free,        /**< Free a braid_Vector*/
           braid_PtFcnSum         sum,         /**< Compute vector sum of two braid_Vectors*/
           braid_PtFcnSpatialNorm spatialnorm, /**< Compute norm of a braid_Vector, this is a norm only over space */
           braid_PtFcnAccess      access,      /**< Allows access to XBraid and current braid_Vector */
           braid_PtFcnBufSize     bufsize,     /**< Computes size for MPI buffer for one braid_Vector */
           braid_PtFcnBufPack     bufpack,     /**< Packs MPI buffer to contain one braid_Vector*/
           braid_PtFcnBufUnpack   bufunpack,   /**< Unpacks MPI buffer into a braid_Vector */
           braid_Core            *core_ptr     /**< Pointer to braid_Core (_braid_Core) struct*/   
           );

/**
 * Carry out a simulation with XBraid.  Integrate in time.
 **/
braid_Int
braid_Drive(braid_Core  core                /**< braid_Core (_braid_Core) struct*/
            );

/**
 * Clean up and destroy core.
 **/
braid_Int
braid_Destroy(braid_Core  core              /**< braid_Core (_braid_Core) struct*/
              );

/**
 * Print statistics after a XBraid run.
 **/
braid_Int
braid_PrintStats(braid_Core  core           /**< braid_Core (_braid_Core) struct*/
                 );

/**
 * Set max number of multigrid levels.
 **/
braid_Int
braid_SetMaxLevels(braid_Core  core,        /**< braid_Core (_braid_Core) struct*/
                   braid_Int   max_levels   /**< maximum levels to allow */
                   );

/**
 * Set max number of multigrid levels.
 **/
braid_Int
braid_SetSkip(braid_Core  core,        /**< braid_Core (_braid_Core) struct*/
              braid_Int   skip         /**< boolean, whether to skip all work on first down-cycle */
              );

/**
 * Turn time refinement on (refine = 1) or off (refine = 0).
 **/
braid_Int
braid_SetRefine(braid_Core  core,    /**< braid_Core (_braid_Core) struct*/
                braid_Int   refine   /**< boolean, refine in time or not */
                );
/**
 * Set the max number of refinements.
 **/
braid_Int
braid_SetMaxRefinements(braid_Core  core,             /**< braid_Core (_braid_Core) struct*/
                        braid_Int   max_refinements   /**< maximum levels to allow */
                       );
/**
 * Set the number of time steps, beyond with refinements stop.
 * If num(tpoints) > tpoints_cutoff, then stop doing refinements.
 **/
braid_Int
braid_SetTPointsCutoff(braid_Core  core,                /**< braid_Core (_braid_Core) struct*/
                       braid_Int   tpoints_cutoff       /**< cutoff for stopping refinements */
                      );

/**
 * Set minimum allowed coarse grid size. XBraid stops coarsening whenever
 * creating the next coarser grid will result in a grid smaller than
 * min_coarse.  The maximum possible coarse grid size will be
 * min_coarse*coarsening_factor.
 **/
braid_Int
braid_SetMinCoarse(braid_Core  core,        /**< braid_Core (_braid_Core) struct*/
                   braid_Int   min_coarse   /** minimum coarse grid size */
                   );

/**
 * Set absolute stopping tolerance.
 *
 * **Recommended option over relative tolerance**
 **/
braid_Int
braid_SetAbsTol(braid_Core  core,           /**< braid_Core (_braid_Core) struct*/
                braid_Real  atol            /**< absolute stopping tolerance */
                );

/**
 * Set relative stopping tolerance, relative to the initial residual.  Be 
 * careful.  If your initial guess is all zero, then the initial residual
 * may only be nonzero over one or two time values, and this will skew the
 * relative tolerance.  Absolute tolerances are recommended.
 **/
braid_Int
braid_SetRelTol(braid_Core  core,           /**< braid_Core (_braid_Core) struct*/
                braid_Real  rtol            /**< relative stopping tolerance */
                );

/**
 * Set the number of relaxation sweeps *nrelax* on grid *level*
 * (level 0 is the finest grid).  The default is 1 on all levels.  To change the
 * default factor, use *level = -1*.  One sweep is a CF relaxation sweep.
 **/
braid_Int
braid_SetNRelax(braid_Core  core,           /**< braid_Core (_braid_Core) struct*/           
                braid_Int   level,          /**< *level* to set *nrelax* on */
                braid_Int   nrelax          /**< number of relaxations to do on *level* */
                );

/**
 * Set the coarsening factor *cfactor* on grid *level* (level 0 is
 * the finest grid).  The default factor is 2 on all levels.  To change the
 * default factor, use *level = -1*.
 **/
braid_Int
braid_SetCFactor(braid_Core  core,          /**< braid_Core (_braid_Core) struct*/
                 braid_Int   level,         /**< *level* to set coarsening factor on */
                 braid_Int   cfactor        /**< desired coarsening factor */
                 );

/**
 * Set max number of multigrid iterations.
 **/
braid_Int
braid_SetMaxIter(braid_Core  core,          /**< braid_Core (_braid_Core) struct*/
                 braid_Int   max_iter       /**< maximum iterations to allow */
                 );

/**
 * Once called, XBraid will use FMG (i.e., F-cycles.
 **/
braid_Int
braid_SetFMG(braid_Core  core               /**< braid_Core (_braid_Core) struct*/
             );

/**
 * Once called, XBraid will use FMG (i.e., F-cycles.
 **/
braid_Int
braid_SetNFMG(braid_Core  core,              /**< braid_Core (_braid_Core) struct*/
              braid_Int   k                  /**< number of initial F-cycles to do before switching to V-cycles */
             );

/**
 * Set number of V-cycles to use at each FMG level (standard is 1)
 **/
braid_Int
braid_SetNFMGVcyc(braid_Core  core,         /**< braid_Core (_braid_Core) struct*/
                  braid_Int   nfmg_Vcyc     /**< number of V-cycles to do each FMG level */
                  );


/**
 * Sets the storage properties of the code.
 **/
braid_Int
braid_SetStorage(braid_Core  core,          /**< braid_Core (_braid_Core) struct*/
                 braid_Int   storage        /**< store C-points (0), all points (1) */
                );

/** 
 * Sets XBraid temporal norm.
 *
 * This option determines how to obtain a global space-time residual norm.
 * That is, this decides how to combine the spatial norms returned by
 * [braid_PtFcnSpatialNorm](@ref braid_PtFcnSpatialNorm) at each time step to
 * obtain a global norm over space and time.  It is this global norm  that 
 * then controls halting.
 *
 * There are three options for setting *tnorm*.  See section @ref halting 
 * for a more detailed discussion (in Introduction.md).
 *
 * - *tnorm=1*: One-norm summation of spatial norms
 * - *tnorm=2*: Two-norm summation of spatial norms 
 * - *tnorm=3*: Infinity-norm combination of spatial norms
 *
 * **The default choice is _tnorm=2_**
 *
 **/
braid_Int
braid_SetTemporalNorm(braid_Core  core,   /**< braid_Core (_braid_Core) struct*/
                      braid_Int   tnorm   /**< choice of temporal norm*/
                      );

/**
 * Set user-defined residual routine.
 **/
braid_Int
braid_SetResidual(braid_Core          core,     /**< braid_Core (_braid_Core) struct*/ 
                  braid_PtFcnResidual residual  /**< function pointer to residual routine */
                  );

/**
 * Set user-defined residual routine for computing full residual norm (all C/F points).
 **/
braid_Int
braid_SetFullRNormRes(braid_Core          core,     /**< braid_Core (_braid_Core) struct*/ 
                      braid_PtFcnResidual residual  /**< function pointer to residual routine */
                      );

/**
 * Set spatial coarsening routine with user-defined routine.
 * Default is no spatial refinment or coarsening.
 **/
braid_Int
braid_SetSpatialCoarsen(braid_Core          core,    /**< braid_Core (_braid_Core) struct*/ 
                        braid_PtFcnSCoarsen scoarsen /**< function pointer to spatial coarsening routine */
                        );

/**
 * Set spatial refinement routine with user-defined routine.
 * Default is no spatial refinment or coarsening.
 **/
braid_Int
braid_SetSpatialRefine(braid_Core         core,   /**< braid_Core (_braid_Core) struct*/
                       braid_PtFcnSRefine srefine /**< function pointer to spatial refinement routine */
                       );

/**
 * Set print level for XBraid.  This controls how much information is 
 * printed to the XBraid print file (@ref braid_SetPrintFile).
 * 
 * - Level 0: no output
 * - Level 1: print typical information like a residual history, 
 *    number of levels in the XBraid hierarchy, and so on.
 * - Level 2: level 1 output, plus debug level output.
 * 
 * Default is level 1.
 **/
braid_Int
braid_SetPrintLevel(braid_Core  core,          /**< braid_Core (_braid_Core) struct*/
                    braid_Int   print_level    /**< desired print level */
                    );

/**
 * Set output file for runtime print messages.  Level of printing is 
 * controlled by @ref braid_SetPrintLevel.  Default is stdout.
 **/
braid_Int
braid_SetPrintFile(braid_Core     core,             /**< braid_Core (_braid_Core) struct*/
                   const char    *printfile_name    /**< output file for XBraid runtime output */
                   );

/** Use default filename, *braid_runtime.out* for runtime print messages.  This
 * function is particularly useful for Fortran codes, where passing filename
 * strings between C and Fortran is troublesome.  Level of printing is
 * controlled by @ref braid_SetPrintLevel.
 **/
braid_Int
braid_SetDefaultPrintFile(braid_Core     core       /**< braid_Core (_braid_Core) struct*/
                          );

/**
 * Set access level for XBraid.  This controls how often the user's
 * access routine is called.
 * 
 * - Level 0:  Never call the user's access routine
 * - Level 1:  Only call the user's access routine after XBraid is finished
 * - Level 2:  Call the user's access routine every iteration and on every level.
 *             This is during _braid_FRestrict, during the down-cycle part of a XBraid iteration. 
 * 
 * Default is level 1.
 **/
braid_Int
braid_SetAccessLevel(braid_Core  core,          /**< braid_Core (_braid_Core) struct*/
                     braid_Int   access_level   /**< desired access_level */
                     );

/**
 * Split MPI commworld into *comm_x* and *comm_t*, the 
 * spatial and temporal communicators.  The total number of processors
 * will equal Px*Pt, there Px is the number of procs in space, and Pt is 
 * the number of procs in time.
 **/
braid_Int
braid_SplitCommworld(const MPI_Comm  *comm_world,  /**< Global communicator to split */
                     braid_Int        px,          /**< Number of processors parallelizing space for a single time step*/
                     MPI_Comm        *comm_x,      /**< Spatial communicator (written as output) */
                     MPI_Comm        *comm_t       /**< Temporal communicator (written as output) */
                     );

/**
 * After Drive() finishes, this returns the number of iterations taken.
 **/
braid_Int
braid_GetNumIter(braid_Core  core,          /**< braid_Core (_braid_Core) struct*/
                 braid_Int   *niter_ptr     /**< output, holds number of iterations taken */
                 );


/** 
 * After Drive() finishes, this returns  XBraid residual history.  If
 * *nrequest_ptr* is negative, return the last *nrequest_ptr* residual norms.
 * If positive, return the first *nrequest_ptr* residual norms.  Upon exit,
 * *nrequest_ptr* holds the number of residuals actually returned.
 **/
braid_Int
braid_GetRNorms(braid_Core  core,           /**< braid_Core (_braid_Core) struct */
                braid_Int   *nrequest_ptr,  /**< input/output, input: num requested resid norms, output: num actually returned */
                braid_Real  *rnorms         /**< output, holds residual norm history array */
                );

/**
 * After Drive() finishes, this returns the number of XBraid levels
 **/
braid_Int
braid_GetNLevels(braid_Core  core,          /**< braid_Core (_braid_Core) struct*/
                 braid_Int  *nlevels_ptr    /**< output, holds the number of XBraid levels */
                 );

/** Example function to compute a tapered stopping tolerance for implicit time
 * stepping routines, i.e., a tolerance *tol_ptr* for the spatial solves.  This
 * tapering only occurs on the fine grid.
 *
 * This rule must be followed.  The same tolerance must be returned over all
 * processors, for a given XBraid and XBraid level.  Different levels may have
 * different tolerances and the same level may vary its tolerance from 
 * iteration to iteration, but for the same iteration and level, the tolerance
 * must be constant.
 *
 * This additional rule must be followed.  The fine grid tolerance is never 
 * reduced (this is important for convergence)
 *
 * On the fine level,the spatial stopping tolerance *tol_ptr* is interpolated
 * from *loose_tol* to *tight_tol* based on the relationship between 
 *     *rnorm* / *rnorm0* and *tol*.  
 * Remember when 
 *     *rnorm* / *rnorm0* < *tol*, 
 * XBraid halts.  Thus, this function lets us have a loose stopping tolerance
 * while the Braid residual is still relatively large, and then we transition
 * to a tight stopping tolerance as the Braid residual is reduced.
 *
 * If the user has not defined a residual function, tight_tol is always returned.
 *
 * The loose_tol is always used on coarse grids, excepting the above mentioned 
 * residual computations.
 *
 * This function will normally be called from the user's step routine.
 *
 * This function is also meant as a guide for users to develop their own 
 * routine.
 *
 **/
braid_Int
braid_GetSpatialAccuracy( braid_StepStatus  status,         /**< Current XBraid step status */
                          braid_Real        loose_tol,      /**< Loosest allowed spatial solve stopping tol on fine grid*/
                          braid_Real        tight_tol,      /**< Tightest allowed spatial solve stopping tol on fine grid*/
                          braid_Real       *tol_ptr         /**< output, holds the computed spatial solve stopping tol */
                         );

/**
 * Set the initial guess to XBraid as the sequential time stepping solution.
 * This is primarily for debugging.  When used with storage=-2, the initial
 * residual should evaluate to exactly 0.  The residual can also be 0 for other
 * storage options if the time stepping is *exact*, e.g., the implicit solve in
 * Step is done to full precision.
 *
 * The value *seq_soln* is a Boolean
 * - 0: The user's Init() function initializes the state vector (default)
 * - 1: Sequential time stepping, with the user's initial condition from
 *      Init(t=0) initializes the state vector
 *
 * Default is 0.
 **/
braid_Int
braid_SetSeqSoln(braid_Core  core,          /**< braid_Core (_braid_Core) struct*/
                 braid_Int   seq_soln       /**< 1: Init with sequential time stepping soln, 0: Use user's Init()*/
                 );
/** @}*/

#ifdef __cplusplus
}
#endif

#endif

