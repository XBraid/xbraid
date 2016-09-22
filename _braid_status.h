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
 

/** \file _braid_status.h
 * \brief Define headers for XBraid status structures, status get/set routines
 * and status create/destroy routines.
 *
 */

#ifndef _braid_status_HEADER
#define _braid_status_HEADER

#include "braid_status.h"

#ifdef __cplusplus
extern "C" {
#endif

struct _braid_Status_struct
{
   /** Common properties */
   braid_Real    t;                /**< current time */
   braid_Int     idx;              /**< time point index value corresponding to t on the global time grid */
   braid_Int     iter;             /**< XBraid iteration number */
   braid_Int     level;            /**< current level in XBraid*/
   braid_Int     nrefine;          /**< number of refinements done */
   braid_Int     gupper;           /**< global size of the fine grid */
   /** Access properties */
   braid_Real    rnorm;            /**< residual norm */
   braid_Int     done;             /**< boolean describing whether XBraid has finished */
   braid_Int     wrapper_test;     /**< boolean describing whether this call is only a wrapper test */
   braid_Int     calling_function; /**< from which function are we accessing the vector */
   /** CoarsenRef properties*/
   braid_Real    f_tprior;         /**< time value to the left of tstart on fine grid */
   braid_Real    f_tstop;          /**< time value to the right of tstart  on fine grid */
   braid_Real    c_tprior;         /**< time value to the left of tstart on coarse grid */
   braid_Real    c_tstop;          /**< time value to the right of tstart on coarse grid */
   /** Step properties */
   braid_Real    tstop;            /**< time value to evolve towards, time value to the right of tstart */
   braid_Real    tol;              /**< Current stopping tolerance */
   braid_Real*   rnorms;           /**< residual norm history, (points to Core->rnorms object) */
   braid_Real    old_fine_tolx;    /**< Allows for storing the previously used fine tolerance from GetSpatialAccuracy */
   braid_Int     tight_fine_tolx;  /**< Boolean, indicating whether the tightest fine tolx has been used, condition for halting */
   braid_Int     rfactor;          /**< if set by user, allows for subdivision of this interval for better time accuracy */
   braid_Int     r_space;          /**< if set by the user, spatial coarsening function will be called following the vcycle */
   /** Buffer properties */
   braid_Int    messagetype;       /**< message type, 0: for Step(), 1: for load balancing */
   braid_Int    size_buffer;       /**< if set by user, send buffer will be "size" bytes in length */
};

/**
 * AccessStatus structure which defines the status of XBraid at a given instant
 * on some level during a run.  The user accesses it through
 * _braid_AccessStatusGet**()_ functions. This is just a pointer to the braid_Status
 **/
struct _braid_AccessStatus_struct
{
   _braid_Status status;
};

/**
 * The user's step routine routine will receive a StepStatus structure, which
 * defines the status of XBraid at the given instant for step evaluation on some level
 * during a run.  The user accesses it through _braid_StepStatusGet**()_ functions.
 * This is just a pointer to the braid_Status.
 **/
struct _braid_StepStatus_struct
{
   _braid_Status status;
};

/**
 * The user coarsen and refine routines will receive a CoarsenRefStatus structure, which
 * defines the status of XBraid at a given instant of coarsening or refinement on some level
 * during a run.  The user accesses it through _braid_CoarsenRefStatusGet**()_ functions.
 * This is just a pointer to the braid_Status.
 **/
struct _braid_CoarsenRefStatus_struct
{
   _braid_Status status;
};

/**
 * The user's bufpack, bufunpack and bufsize routines will receive a BufferStatus structure, which
 * defines the status of XBraid at a given buff (un)pack instance.  The user accesses it
 * through _braid_BufferStatusGet**()_ functions. This is just a pointer to the braid_Status.
 **/
struct _braid_BufferStatus_struct
{
   _braid_Status status;
};


#ifdef __cplusplus
}
#endif

#endif
