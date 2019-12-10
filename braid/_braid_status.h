/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
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
 

/** \file _braid_status.h
 * \brief Define the internals of the XBraid status structures and internal 
 * status structure functions, like destroy. 
 */

#ifndef _braid_status_HEADER
#define _braid_status_HEADER

#include "braid_status.h"
#include "_braid.h"

#ifdef __cplusplus
extern "C" {
#endif


/*--------------------------------------------------------------------------
 * Define base Status structure as a pointer to core, and all other derived
 * Status structures as pointers to the base class.  
 *
 * See braid_status.h for a description of each Status structure.
 *--------------------------------------------------------------------------*/

struct _braid_Status_struct
{
   _braid_Core core;
};
typedef struct _braid_Status_struct _braid_Status;

struct _braid_AccessStatus_struct
{
   _braid_Status status;
};

struct _braid_SyncStatus_struct
{
   _braid_Status status;
};

struct _braid_StepStatus_struct
{
   _braid_Status status;
};

struct _braid_CoarsenRefStatus_struct
{
   _braid_Status status;
};

struct _braid_BufferStatus_struct
{
   _braid_Status status;
};

struct _braid_ObjectiveStatus_struct
{
   _braid_Status status;
};


/*--------------------------------------------------------------------------
 * Begin headers for internal Braid Status functions, like Destroy, and StatusInit
 *--------------------------------------------------------------------------*/

#define _braid_StatusElt(status, elt) ( ((braid_Core)status) -> elt )

braid_Int
_braid_StatusDestroy(braid_Status status);

/**
 * Initialize a braid_AccessStatus structure
 */
braid_Int
_braid_AccessStatusInit(braid_Real          t,                /**< current time */
                        braid_Int           idx,              /**< time point index value corresponding to t on the global time grid */
                        braid_Real          rnorm,            /**< current residual norm in XBraid */
                        braid_Int           iter,             /**< current iteration in XBraid*/
                        braid_Int           level,            /**< current level in XBraid */
                        braid_Int           nrefine,          /**< number of refinements done */
                        braid_Int           gupper,           /**< global size of the fine grid */
                        braid_Int           done,             /**< boolean describing whether XBraid has finished */
                        braid_Int           wrapper_test,     /**< boolean describing whether this call is only a wrapper test */
                        braid_Int           calling_function, /**< from which function are we accessing the vector */
                        braid_AccessStatus  status            /**< structure to initialize */
                        );

/**
 * Initialize a braid_SyncStatus structure
 */
braid_Int
_braid_SyncStatusInit(braid_Int           iter,             /**< current iteration in XBraid*/
                      braid_Int           level,            /**< current level in XBraid */
                      braid_Int           nrefine,          /**< number of refinements done */
                      braid_Int           gupper,           /**< global size of the fine grid */
                      braid_Int           done,             /**< boolean describing whether XBraid has finished */
                      braid_Int           calling_function, /**< from which function are we accessing braid */
                      braid_SyncStatus    status            /**< structure to initialize */
                      );

/**
 * Initialize a braid_CoarsenRefStatus structure
 */
braid_Int
_braid_CoarsenRefStatusInit(braid_Real              tstart,      /**< time value for current vector */
                            braid_Real              f_tprior,    /**< time value to the left of tstart on fine grid */
                            braid_Real              f_tstop,     /**< time value to the right of tstart on fine grid */
                            braid_Real              c_tprior,    /**< time value to the left of tstart on coarse grid */
                            braid_Real              c_tstop,     /**< time value to the right of tstart on coarse grid */
                            braid_Int               level,       /**< current fine level in XBraid */
                            braid_Int               nrefine,     /**< number of refinements done */
                            braid_Int               gupper,      /**< global size of the fine grid */
                            braid_Int               c_index,     /**< coarse time index refining from */
                            braid_CoarsenRefStatus  status       /**< structure to initialize */
                            );

/**
 * Initialize a braid_StepStatus structure
 */
braid_Int
_braid_StepStatusInit(braid_Real        tstart,      /**< current time value  */
                      braid_Real        tstop,       /**< time value to evolve towards, time value to the right of tstart */
                      braid_Int         idx,         /**< time point index value corresponding to tstart on the global time grid */
                      braid_Real        tol,         /**< Current XBraid stopping tolerance */
                      braid_Int         iter,        /**< Current XBraid iteration (also equal to length of rnorms) */
                      braid_Int         level,       /**< current level in XBraid */
                      braid_Int         nrefine,     /**< number of refinements done */
                      braid_Int         gupper,      /**< global size of the fine grid */
                      braid_StepStatus  status       /**< structure to initialize */
                      );

/**
 * Initialize a braid_BufferStatus structure 
 */
braid_Int
_braid_BufferStatusInit(braid_Int           messagetype,  /**< message type, 0: for Step(), 1: for load balancing */
                        braid_Int           size,         /**< if set by user, size of send buffer is "size" bytes */
                        braid_BufferStatus  status        /**< structure to initialize */
                        );


/**
 * Initialize a braid_ObjectiveStatus structure */
braid_Int
_braid_ObjectiveStatusInit(braid_Real            tstart,
                           braid_Int             idx,
                           braid_Int             iter,
                           braid_Int             level,
                           braid_Int             nrefine,
                           braid_Int             gupper,
                           braid_ObjectiveStatus status);                     

#ifdef __cplusplus
}
#endif

#endif
