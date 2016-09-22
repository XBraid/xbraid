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
#include "_braid.h"

#ifdef __cplusplus
extern "C" {
#endif

struct _braid_Status_struct
{
   _braid_Core core;
};
typedef struct _braid_Status_struct _braid_Status;

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
