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
 * \brief Header files for Block Distrabution functions. 
 */
#ifndef _blockdist_HEADER
#define _blockdist_HEADER
#include "_braid.h"

/** 
 * Returns the Balance struct neccesary to create a new grid
 * with load balancing and/or time refinment 
 */

braid_Int
_braid_BlockDist(braid_Core            core,
                 braid_Int             *done,
                 braid_Int              npoints,
                 _braid_BalanceStruct  *bstruct );


/**
 * Used by BlockDist, fills in the new fine time interval 
 */

braid_Int
_braid_GetBlockDistInterval1(braid_Core core,
                             _braid_BalanceStruct *bstruct );

/**
 * Returns the index interval for *proc* in a blocked data distribution.
 */
braid_Int
_braid_GetBlockDistInterval(braid_Int   npoints,
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

#endif 
