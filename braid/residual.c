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

#include "_braid.h"
#include "util.h"

/*----------------------------------------------------------------------------
 * Compute residual
 *----------------------------------------------------------------------------*/

braid_Int
_braid_Residual(braid_Core        core,
                braid_Int         level,
                braid_Int         index,
                braid_Int         calling_function,
                braid_BaseVector  ustop,
                braid_BaseVector  r)
{
   braid_App        app      = _braid_CoreElt(core, app);
   braid_Real       tol      = _braid_CoreElt(core, tol);
   braid_Int        iter     = _braid_CoreElt(core, niter);
   _braid_Grid    **grids    = _braid_CoreElt(core, grids);
   braid_StepStatus status   = (braid_StepStatus)core;
   braid_Int        nrefine  = _braid_CoreElt(core, nrefine);
   braid_Int        gupper   = _braid_CoreElt(core, gupper);
   braid_Int        ilower   = _braid_GridElt(grids[level], ilower);
   braid_Real      *ta       = _braid_GridElt(grids[level], ta);

   braid_Int  delta_correct = _braid_CoreElt(core, delta_correct);
   braid_Int  delta_rank    = _braid_CoreElt(core, delta_rank);

   braid_BaseVector rstop;
   braid_Int        ii;

   ii = index-ilower;
   _braid_StepStatusInit(ta[ii-1], ta[ii], index-1, tol, iter, level, nrefine, gupper, calling_function, r->basis, status);
   if ( _braid_CoreElt(core, residual) == NULL )
   {
      /* By default: r = ustop - \Phi(ustart)*/
      _braid_GetUInit(core, level, index, r, &rstop);
      _braid_BaseStep(core, app,  rstop, NULL, r, level, status);
      _braid_BaseSum(core, app,  1.0, ustop, -1.0, r);
   }
   else
   {
      /* Call the user's residual routine */
      _braid_BaseResidual(core, app, ustop, r, status);
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Compute FAS residual = f - residual
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FASResidual(braid_Core        core,
                   braid_Int         level,
                   braid_Int         index,
                   braid_BaseVector  ustop,
                   braid_BaseVector  r)
{
   braid_App          app    = _braid_CoreElt(core, app);
   _braid_Grid      **grids  = _braid_CoreElt(core, grids);
   braid_Int          ilower = _braid_GridElt(grids[level], ilower);
   braid_BaseVector  *fa     = _braid_GridElt(grids[level], fa);

   braid_Int ii;

   _braid_Residual(core, level, index, braid_ASCaller_FASResidual, ustop, r);

   if ( (level == 0) || ((ii = index-ilower) && (fa[ii] == NULL)) )
   {
      _braid_BaseSum(core, app,  0.0, r, -1.0, r);
   }
   else
   {
      _braid_BaseSum(core, app,  1.0, fa[ii], -1.0, r);
   }

   return _braid_error_flag;
}

