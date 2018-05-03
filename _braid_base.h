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


 /** \file _braid_base.h
 *  \brief Define headers for wrapper routines of user-defined functions. 
 * 
 *  Routines defined here wrap the user-defined routines. If this is a primal 
 * XBraid run, they just call the user's routines. If this is an XBraid_Adjoint
 * run, they also record themselves to the action tape and push state
 * and bar vectors to the primal and the bar tape, respectively. *_diff routines 
 * then pop these vectors from the tape and perform the differentiated actions.
 **/


#ifndef _braid_base_HEADER
#define _braid_base_HEADER

#include "_braid.h"


/**
 * This calls the user's step routine.
 * If (adjoint): record the action, and push state and bar vector to primal and bar tapes. 
 */
braid_Int 
_braid_BaseStep(braid_Core       core,       /**< braid_Core structure */                 
                braid_App        app,        /**< user-defined _braid_App structure */  
                braid_BaseVector ustop,      /**< input, u vector at *tstop* */
                braid_BaseVector fstop,      /**< input, right-hand-side at *tstop* */
                braid_BaseVector u,          /**< input/output, initially u vector at *tstart*, upon exit, u vector at *tstop* */   
                braid_Int        level,      /**< current time grid level */ 
                braid_StepStatus status );   /**< braid_Status structure (pointer to the core) */    


/**
 * This initializes a braid_BaseVector and calls the user's 
 * init routine. 
 * If (adjoint): record the action, initialize barVector with zero and push to the bar tape. 
 */
braid_Int
_braid_BaseInit(braid_Core         core,     /**< braid_Core structure */
                braid_App          app,      /**< user-defined _braid_App structure */ 
                braid_Real         t,        /**< current time value for (*u_ptr) */
                braid_BaseVector  *u_ptr     /**< output, newly allocated and initialized vector */
                );

/**
 * This initializes a braid_BaseVector and calls the user's 
 * clone routine.
 * If (adjoint): record the action, initialize a barVector with 
 * zero and push to the bar tape
 */ 
braid_Int
_braid_BaseClone(braid_Core         core,     /**< braid_Core structure */
                 braid_App          app,      /**< user-defined _braid_App structure */ 
                 braid_BaseVector   u,        /**< vector to clone */ 
                 braid_BaseVector  *v_ptr     /**< output, newly allocated and cloned vector */ 
                 );
/**
 * This calls the user's free routine.
 * If (adjoint): record the action, and free the bar vector. 
 */ 
braid_Int
_braid_BaseFree(braid_Core        core,      /**< braid_Core structure */
                braid_App         app,       /**< user-defined _braid_App structure */
                braid_BaseVector  u          /**< vector to free */ 
                );

/**
 * This calls the user's sum routine.
 * If (adjoint): record the action, and push to the bar tape. 
 */ 
braid_Int
_braid_BaseSum(braid_Core        core,       /**< braid_Core structure */
               braid_App         app,        /**< user-defined _braid_App structure */ 
               braid_Real        alpha,      /**< scalar for AXPY */ 
               braid_BaseVector  x,          /**< vector for AXPY */
               braid_Real        beta,       /**< scalar for AXPY */ 
               braid_BaseVector  y           /**< output and vector for AXPY */ 
               );

/**
 * This calls the user's SpatialNorm routine. 
 */ 
braid_Int
_braid_BaseSpatialNorm(braid_Core        core,      /**< braid_Core structure */
                       braid_App         app,       /**< user-defined _braid_App structure */     
                       braid_BaseVector  u,         /**< vector to norm */
                       braid_Real       *norm_ptr   /**< output, norm of braid_Vector (this is a spatial norm) */
                       );

/**
 * This calls the user's Access routine. 
 * If (adjoint): record the action
 */
braid_Int
_braid_BaseAccess(braid_Core          core,        /**< braid_Core structure */
                  braid_App           app,         /**< user-defined _braid_App structure */
                  braid_BaseVector    u,           /**< vector to be accessed */
                  braid_AccessStatus  status       /**< can be querried for info like the current XBraid Iteration */
                  );

/** 
 * This calls the user's BufSize routine.
 */
braid_Int
_braid_BaseBufSize(braid_Core          core,          /**< braid_Core structure */
                   braid_App           app,           /**< user-defined _braid_App structure */
                   braid_Int          *size_ptr,      /**< upper bound on vector size in bytes */
                   braid_BufferStatus  status         /**< can be querried for info on the message type */
                   );      


/** 
 * This calls the user's BufPack routine.
 * If (adjoint): record the action, and push to the bar tape. 
 */
braid_Int
_braid_BaseBufPack(braid_Core          core,      /**< braid_Core structure */
                   braid_App           app,       /**< user-defined _braid_App structure */
                   braid_BaseVector    u,         /**< vector to back into buffer */
                   void               *buffer,    /**< output, MPI buffer containing u */
                   braid_BufferStatus  status     /**< can be queeried for info on the message type required */
                   );

/** 
 * This calls the user's BufUnPack routine.
 * If (adjoint): record the action, initialize the bar vector 
 * with zero and push it to the bar tape. 
 */
braid_Int
_braid_BaseBufUnpack(braid_Core           core,      /**< braid_Core structure */   
                     braid_App            app,       /**< user-defined _braid_App structure */    
                     void                *buffer,    /**< MPI Buffer to unpack and place in u_ptr */ 
                     braid_BaseVector    *u_ptr,     /**< output, braid_Vector containing buffer's data */
                     braid_BufferStatus   status     /**< braid_Status structure (pointer to the core) */
                     );

/** 
 * If (adjoint): This calls the user's ObjectiveT routine.
 * Further: record the action, push to the state and bar tapes. 
 */
braid_Int
_braid_BaseObjectiveT(braid_Core             core,     /**< braid_Core structure */
                      braid_App              app,      /**< user-defined _braid_App structure */
                      braid_BaseVector       u,        /**< input: state vector at current time */
                      braid_ObjectiveStatus  ostatus,  /**< status structure for querying time, index, etc. */
                      braid_Real            *objT_ptr  /**< output: objective function value at current time */
                      );


/** 
 * This calls the user's Residual routine.
 */
braid_Int
_braid_BaseResidual(braid_Core       core,             /**< braid_Core structure */
                    braid_App        app,              /**< user-defined _braid_App structure */
                    braid_BaseVector ustop,            /**< input, u vector at *tstop* */
                    braid_BaseVector r,                /**< output, residual at *tstop* (at input, equals *u* at *tstart*) */
                    braid_StepStatus status            /**< braid_Status structure (pointer to the core) */ 
                    );

/** 
 * This calls the user's FullResidual routine (full_rnorm_res).
 */
braid_Int
_braid_BaseFullResidual(braid_Core        core,        /**< braid_Core structure */
                        braid_App         app,         /**< user-defined _braid_App structure */
                        braid_BaseVector  r,           /**< output, residual at *tstop* */
                        braid_BaseVector  u,           /**< input, u vector at *tstop* */
                        braid_StepStatus  status       /**< braid_Status structure (pointer to the core) */ 
                       );    

/**
 * This initializes a baseVector and calls the user's SCoarsen routine. 
 */
braid_Int
_braid_BaseSCoarsen(braid_Core              core,      /**< braid_Core structure */
                    braid_App               app,       /**< user-defined _braid_App structure */
                    braid_BaseVector        fu,        /**< braid_BaseVector to refine*/
                    braid_BaseVector       *cu_ptr,    /**< output, refined vector */   
                    braid_CoarsenRefStatus  status     /**< braid_Status structure (pointer to the core) */ 
                    );

/**
 * This initializes a baseVector and calls the user's SRefine routine. 
 */
braid_Int
_braid_BaseSRefine(braid_Core              core,      /**< braid_Core structure */
                   braid_App               app,       /**< user-defined _braid_App structure */
                   braid_BaseVector        cu,        /**< braid_BaseVector to refine*/
                   braid_BaseVector       *fu_ptr,    /**< output, refined vector */       
                   braid_CoarsenRefStatus  status     /**< braid_Status structure (pointer to the core) */ 
                   );

/**
 * This initializes a baseVector and call's the user's SInit routine.
 */
braid_Int
_braid_BaseSInit(braid_Core        core,        /**< braid_Core structure */
                 braid_App         app,         /**< user-defined _braid_App structure */
                 braid_Real        t,           /**< time value for *u_ptr* */
                 braid_BaseVector *u_ptr        /**< output, newly allocated and initialized vector shell */
                 );

/**
 * This initializes a baseVector and call's the user's SClone routine.
 */
braid_Int
_braid_BaseSClone(braid_Core         core,      /**< braid_Core structure */ 
                  braid_App          app,       /**< user-defined _braid_App structure */
                  braid_BaseVector   u,         /**< vector to clone */ 
                  braid_BaseVector  *v_ptr      /**< output, newly allocated and cloned vector shell */
                  );

/**
 * This call's the user's SFree routine.
 */
braid_Int
_braid_BaseSFree(braid_Core        core,        /**< braid_Core structure */
                 braid_App         app,         /**< user-defined _braid_App structure */
                 braid_BaseVector  u            /**< vector to free (keeping the shell) */
                 );


/**
 * This calls the user's TimeGrid routine. 
 */
braid_Int
_braid_BaseTimeGrid(braid_Core   core,      /**< braid_Core structure */
                    braid_App    app,       /**< user-defined _braid_App structure */
                    braid_Real  *ta,        /**< temporal grid on level 0 (slice per processor) */
                    braid_Int   *ilower,    /**< lower time index value for this processor */
                    braid_Int   *iupper     /**< upper time index value for this processor */
                    );


/*-------------------------------------------------------------
 * XBraid_Adjoint: Differentiated routines 
 *-----------------------------------------------------------*/

/**
 * This pops state and bar vectors from the tape and calls the user's
 * differentiated step routine (step_diff):
 * ubar = ( d(step)/d(u) )^T * ubar
 */
braid_Int
_braid_BaseStep_diff(_braid_Action *action      /**< _braid_Action structure, holds information about the primal XBraid action */
                     );


/**
 * This pops bar vectors from the tape and performs the 
 * differentiated clone action:
 *  ubar += vbar
 *  vbar  = 0.0
 */ 
braid_Int
_braid_BaseClone_diff(_braid_Action *action      /**< _braid_Action structure, holds information about the primal XBraid action */
                      );

/**
 * This pops bar vectors from the tape and performs the 
 * differentiated sum action:
 *  xbar += alpha * ybar
 *  ybar  = beta  * ybar
 */ 
braid_Int
_braid_BaseSum_diff(_braid_Action *action      /**< _braid_Action structure, holds information about the primal XBraid action */
                    );

/**
 * This pops state and bar vectors from the tape and calls the user's
 * differentiated ObjectiveT routine (objT_diff):
 * ubar = ( d(objectiveT)/d(u) )^T * f_bar
 */ 
braid_Int
_braid_BaseObjectiveT_diff(_braid_Action *action      /**< _braid_Action structure, holds information about the primal XBraid action */
                           );


/**
 * This pops the bar vector from the tape and performs the differentiated BufPack action:
 * MPI_Recv(utmp)
 * ubar += utmp
 */ 
braid_Int
_braid_BaseBufPack_diff(_braid_Action *action      /**< _braid_Action structure, holds information about the primal XBraid action */
                        );

/**
 * This pops the bar vector from the tape and performs the differentiated BufUnPack action:
 * MPI_Send(ubar);
 * ubar = 0.0
 */ 
braid_Int
_braid_BaseBufUnpack_diff(_braid_Action *action      /**< _braid_Action structure, holds information about the primal XBraid action */
                          );


/**
 * This pops the bar vector from the tape and call's the user's
 * differentiated init routine (init_diff), if set.
 */ 
braid_Int
_braid_BaseInit_diff(_braid_Action *action      /**< _braid_Action structure, holds information about the primal XBraid action */
                     );

#endif