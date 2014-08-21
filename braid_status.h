/*BHEADER**********************************************************************
 * Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of XBraid.  See file COPYRIGHT for details.
 *
 * XBraid is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 ***********************************************************************EHEADER*/

/** \file braid_status.h
 * \brief Define headers for Braid status structures, status get/set routines
 * and status create/destroy routines.
 *
 */

#ifndef braid_status_HEADER
#define braid_status_HEADER

#include "braid_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------------
 * Routines for user to access Braid status structures
 *--------------------------------------------------------------------------*/
/** \defgroup braidstatus Braid status routines
 *  
 *  Braid status structures are what the user accesses to determine the status
 *  of the simulation when their routines (phi, coarsen/refine, access) are 
 *  called. 
 *
 *  @{
 */

struct _braid_Status_struct;
/**
 * The user routines will receive braid_Status structures, which will be 
 * pointers to the actual _braid_Status_struct
 **/
typedef struct _braid_Status_struct *braid_Status;


/**
 * Points to the status structure which defines the status of Braid
 * at a given instant on a some level during a run.  The user accesses
 * it through braid_Get**Status() functions.
 **/
typedef struct _braid_Status_struct
{
   braid_Int     iter;         /**< Braid iteration number */
   braid_Int     level;        /**< current level in Braid*/
   braid_Real    rnorm;        /**< residual norm */
   braid_Int     done;         /**< boolean describing whether Braid has finished */
   
} _braid_Status;


/*--------------------------------------------------------------------------
 * Accessor macros 
 *--------------------------------------------------------------------------*/

/**
 * Accessor for _braid_Status attributes 
 **/
#define _braid_StatusElt(status, elt) ( (status) -> elt )


/*--------------------------------------------------------------------------
 * Prototypes
 *--------------------------------------------------------------------------*/

/**
 * Initialize a braid_Status structure in *status_ptr*, setting
 * the status values *rnorm*, *iter*, *level*, *done*
 **/
braid_Int
_braid_InitStatus(braid_Real        rnorm,
                  braid_Int         iter,
                  braid_Int         level,
                  braid_Int         done,
                  braid_Status     *status_ptr);

/**
 * Destroy a braid_Status structure
 **/
braid_Int
_braid_DestroyStatus(braid_Status  status);

/**
 * Return the residual for the current status object.
 **/
braid_Int
braid_GetStatusResidual(braid_Status  status,     /**< structure containing current simulation info */
                        braid_Real   *rnorm_ptr   /**< output, current residual norm */
                        );

/**
 * Return the iteration for the current status object.
 **/
braid_Int
braid_GetStatusIter(braid_Status  status,         /**< structure containing current simulation info */
                    braid_Int    *iter_ptr        /**< output, current iteration number*/
                    );

/**
 * Return the Braid level for the current status object.
 **/
braid_Int
braid_GetStatusLevel(braid_Status  status,        /**< structure containing current simulation info */
                     braid_Int    *level_ptr      /**< output, current level in Braid */
                     );

/**
 * Return whether Braid is done for the current status object\n
 * *done_ptr = 1* indicates that Braid has finished iterating, 
 * (either maxiter has been reached, or the tolerance has been met).
 **/
braid_Int
braid_GetStatusDone(braid_Status  status,         /**< structure containing current simulation info */
                    braid_Int    *done_ptr        /**< output,  =1 if Braid has finished, else =0 */
                    );


/** @}*/

#ifdef __cplusplus
}
#endif

#endif


