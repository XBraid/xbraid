/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory.
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
 
#ifndef _braid_adjoint_HEADER
#define _braid_adjoint_HEADER

#include "_braid.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Shallow copy a braid_VectorBar shared pointer, bar_ptr is set to bar
 * and the useCount is incremented by one.
 */
braid_Int
_braid_VectorBarCopy(braid_VectorBar  bar,       
                     braid_VectorBar *bar_ptr);  

/**
 * Reduce the useCount of a braid_VectorBar shared pointer 
 * Free the pointer memory if useCount is zero.  
 */ 
braid_Int
_braid_VectorBarDelete(braid_Core      core, 
                       braid_VectorBar bar);

/**
 * Free memory of the optimization structure 
 */
braid_Int
_braid_OptimDestroy( braid_Core core);

/**
 * Update the adjoint variables and compute adjoint residual norm
 * Returns the tnorm of adjoint residual
 */
braid_Int
_braid_UpdateAdjoint(braid_Core  core,
                     braid_Real *rnorm_adj_ptr);

/**
 * Set adjoint residual norm 
 */
braid_Int 
_braid_SetRNormAdjoint(braid_Core  core, 
                       braid_Int   iter, 
                       braid_Real  rnorm_adj);

/** 
 * Evaluate the user's local objective function at time *t* and add it to the
 * time-averaged objective function
 */
braid_Int
_braid_AddToObjective(braid_Core             core, 
                      braid_BaseVector       u, 
                      braid_ObjectiveStatus  ostatus);

/**
 * Evaluate the objective function:
 * MPI_Allreduce the time average and postprocess the objective 
 */
braid_Int
_braid_EvalObjective(braid_Core core);

/** 
 * Differentiated objective function 
 */
braid_Int
_braid_EvalObjective_diff(braid_Core core);

/**
 * Allocate and initialize the adjoint variables 
 */
braid_Int
_braid_InitAdjointVars(braid_Core   core, 
                       _braid_Grid *fine_grid);

/**
 * Sanity check for non-supported adjoint features
 */
braid_Int
_braid_AdjointFeatureCheck(braid_Core core);

#ifdef __cplusplus
}
#endif

#endif
