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

/** \file delta.h
 * \brief Define internal XBraid headers for Delta correction.
 *
 * This file contains the internal XBraid headers for Delta correction,
 */

#ifndef _braid_DELTA_HEADER
#define _braid_DELTA_HEADER

#include "_braid.h"

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * macro for determining when to compute Delta correction
 */
#define _braid_DoDeltaCorrect(core, level, niter) \
( _braid_CoreElt(core, delta_correct) && niter >= _braid_CoreElt(core, delta_defer_iter) && level >= _braid_CoreElt(core, delta_defer_lvl) )

/**
 * macro for determining when we can use computed Delta correction
 */
#define _braid_UseDeltaCorrect(core, level, niter) \
( _braid_CoreElt(core, delta_correct) && niter >= _braid_CoreElt(core, delta_defer_iter) && level > _braid_CoreElt(core, delta_defer_lvl) )

/**
 * Compute the action of the low-rank approximation to Delta on a vector
 */

braid_Int
_braid_LRDeltaDot(braid_Core core,
                  braid_App app,
                  braid_Vector u,
                  braid_Basis delta,
                  braid_Basis basis);

/**
 * Compute the action of the low-rank approximation to Delta on a basis
 */

braid_Int
_braid_LRDeltaDotMat(braid_Core core,
                     braid_App app,
                     braid_Basis psi,
                     braid_Basis delta,
                     braid_Basis basis);

                  
/**
 * Perform modified Gram-Schmidt orthonormalization on a basis, while also
 * computing local Lyapunov exponents
 */
braid_Int
_braid_GramSchmidt(braid_Core   core,
                     braid_App    app,
                     braid_Basis  basis,
                     braid_Real  *exps);


#ifdef __cplusplus
}
#endif

#endif