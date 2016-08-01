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

/** \file blockdist.h
 * \brief Header files for the data distribution functions.
 */
#ifndef _braiddist_HEADER
#define _braiddist_HEADER
#include "_braid.h"

/**************************************************************************
 * Block Data distribution functions. These functions use an O(1) global
 * function to distribute time points evenly across the available time
 * processors. This is the defualt distribution used by xbraid, and is also
 * the distribution used to calculate the assumed partition, used to create
 * communication maps when a non global distribution is used (such as the
 * weighted distribution below). To save on communication, the data distribution
 * functions must also calculate the new fine grid intervals if time refinment
 * is to be completed.
 **************************************************************************/

/**
 * Returns a bstruct containing all the necessary information to
 * distribute time steps evenly across temporal processors, on a new
 * (possibly refined ) grid, using the block data distribution.
 */

braid_Int
_braid_BlockDist(braid_Core            core,
                 braid_Int             *done,
                 braid_Int              npoints,
                 _braid_BalanceStruct  *bstruct );


/**
 * Helper function for _braid_BlockDist. Calculates the new fine
 * grid interval, and stores it inside the _braid_BalanceStruct.
 */

braid_Int
_braid_GetBlockDistInterval(braid_Core core,
                             _braid_BalanceStruct *bstruct );

/**
 * Returns the index interval for *proc* in a blocked data distribution.
 */

braid_Int
_braid_GetBlockDistInterval_basic(braid_Int   npoints,
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
                        braid_Int  *proc_ptr);

/*******************************************************************************
 * Weighted Data distribution: The following functions impliment the weigthed data
 * distribution. These functions calculate the new data distribution when using
 * load balancing based on the user defined wfactors. This distribution balances the
 * total weight across all processors.
 ********************************************************************************/

/* _braid_WeightedStruct is used to pass information about the weights to the corresponding
 * functions. This is used for simplicicty, as it allows all communication needed during
 * the load balancing and temporal refinment steps to be completed at the same time.
 */

typedef struct
{
    braid_Real local_max;   /* local max of the weights */
    braid_Real local_min;   /* local min of the weights */
    braid_Real local_sum;   /* local sum of the weights */
    braid_Real global_max;  /* global max of the weights */
    braid_Real global_min;  /* global min of the weights */
    braid_Real global_sum;  /* global sum of the weights */
    braid_Real local_start; /* cumulative sum of the weights exclusive */
    braid_Real local_stop;  /* cumulative sum of the weights inclusive */

} _braid_WeightedStruct;

/**
 * Accessor for _braid_WeightedStruct attibutes
 **/

#define _braid_WeightedElt( weighted, elt)  ( (weighted) -> elt )

/**
 * Init the weighted struct class
 */
braid_Int
_braid_WeightedStructInit( _braid_WeightedStruct *wstruct );


/**
 * Destroy the _braid_WeightedStruct structure
 */
braid_Int
_braid_WeightedStructDestroy( _braid_WeightedStruct *wstruct );

/**
 * Build the Balance structure required to load balance on a (possibly refined )
 * new grid, using a weighted distribution.
 */

braid_Int
_braid_WeightedDist(braid_Core            core,
                    braid_Int             *done,
                    braid_Real             *wfactors,
                    braid_Int              npoints,
                    _braid_BalanceStruct  *bstruct );

/**
 * Helper function for _braid_WeightedDist. Gets the new fine grid interval
 * based on the user defined weights, and stores that information in the
 * Balance structure.
 */

braid_Int
_braid_GetWeightedInterval(braid_Core core,
                           braid_Real  *wfactors,
                           _braid_WeightedStruct *wstruct,
                           _braid_BalanceStruct *bstruct );

/**
 * MPI User function, used to combine 4 all reduce commands into a single call.
 * In this case the function accepts an array of length 4. It returns the sum,
 * of the first elements, the sum of the second elements, the max of the third elements
 * and the min of the fourth elements.
 */

void SumSumMaxMin(braid_Real *in, braid_Real *inout, braid_Int *len, MPI_Datatype *datatype);

#endif
