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
 * \brief Header files for Weighted Distrabution functions. 
 */

#ifndef weighteddist_HEADER
#define weighteddist_HEADER
#include "_braid.h"
#include "blockdist.h"

/**
 * Accessor for _braid_WeightedStruct attibutes
 **/
#define _braid_WeightedElt( weighted, elt)  ( (weighted) -> elt )

typedef struct
{
    braid_Int local_max;  
    braid_Int local_min;
    braid_Int local_sum;
    braid_Int global_max;
    braid_Int global_min;
    braid_Int global_sum;
    braid_Int local_start; /* cumulative sum of the weights exclusive */
    braid_Int local_stop;  /* cumulative sum of the weights inclusive */

} _braid_WeightedStruct; 

/** 
 * Init the weighted struct class );
 */ 
braid_Int
_braid_WeightedStructInit( _braid_WeightedStruct *wstruct );


/**
 * Destroy the _braid_WeightedStruct structure
 */
braid_Int
_braid_WeightedStructDestroy( _braid_WeightedStruct *wstruct );

/** 
 * Build the Balance structure using a weighted distrobution 
 */
braid_Int
_braid_WeightedDist(braid_Core            core,
                           braid_Int             *done,
                           braid_Int             *wfactors,
                           braid_Int              npoints,
                           _braid_WeightedStruct *wstruct,
                           _braid_BalanceStruct  *bstruct );

/**
 * Get the new fine grid interval based on the weighted distrobution
 */ 
braid_Int
_braid_GetWeightedInterval(braid_Core core,
                           braid_Int  *wfactors,
                           _braid_WeightedStruct *wstruct, 
                           _braid_BalanceStruct *bstruct );

void SumSumMaxMin(int *in, int *inout, int *len, MPI_Datatype *datatype);

#endif 
