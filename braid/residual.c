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

   /* this decides whether Delta correction needs to be computed on this level,
   *  which depends on from where the residual is called
   */
   braid_Int delta_correct;
   delta_correct = _braid_CoreElt(core, delta_correct)
                   &&  ( iter  >= _braid_CoreElt(core, delta_defer_iter) )
                   &&  ( level >= _braid_CoreElt(core, delta_defer_lvl) )
                   && !( calling_function == braid_ASCaller_Residual
                      && level == _braid_CoreElt(core, delta_defer_lvl) );

   braid_BaseVector rstop;
   braid_Int        ii;

   ii = index-ilower;
   
   if ( delta_correct )
   {  /* Give the user access to the basis vectors through StepStatusGetBasisVec */
      _braid_StepStatusInit(ta[ii-1], ta[ii], index-1, tol, iter, level, nrefine, gupper, calling_function, r->basis, status);

      /* By default: r = ustop - \Phi(ustart)*/
      _braid_GetUInit(core, level, index, r, &rstop);
      _braid_BaseStep(core, app, rstop, NULL, r, level, status);
      _braid_CoreFcn(core, sum)(app, 1.0, ustop->userVector, -1.0, r->userVector);
      _braid_BaseSumBasis(core, app, 0., r->basis, -1.0, r->basis);

      return _braid_error_flag;
   }
   /* else, default behavior */

   /* initialize status struct */
   _braid_StepStatusInit(ta[ii-1], ta[ii], index-1, tol, iter, level, nrefine, gupper, calling_function, NULL, status);
   if ( _braid_CoreElt(core, residual) == NULL )
   {
      /* By default: r = ustop - \Phi(ustart)*/
      _braid_GetUInit(core, level, index, r, &rstop);
      _braid_BaseStep(core, app, rstop, NULL, r, level, status);
      _braid_BaseSum(core, app, 1.0, ustop, -1.0, r);
   }
   else   /* can we make residual option compatible with Delta correction? */
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
   braid_BaseVector  *va     = _braid_GridElt(grids[level], va);
   // braid_Basis       *ba     = _braid_GridElt(grids[level], ba);


   braid_Int ii = index-ilower;

   /* && short circuits, so (fa[ii] != NULL) is never evaluated if (level == 0) */
   braid_Int delta_correct = _braid_CoreElt(core, delta_correct)
                             && (_braid_CoreElt(core, niter) >= _braid_CoreElt(core, delta_defer_iter))
                             && (level > _braid_CoreElt(core, delta_defer_lvl))
                             && (fa[ii] != NULL);

   braid_BaseVector delta;    /* temporary storage for delta correction */

   if ( delta_correct )
   {  /* compute the Delta correction for the state and Lyapunov vectors */
      _braid_BaseClone(core, app, r, &delta);
      // _braid_LRDeltaDot(core, app, delta->userVector, fa[ii]->basis, ba[ii]);
      // _braid_LRDeltaDotMat(core, app, delta->basis, fa[ii]->basis, ba[ii]);
      _braid_LRDeltaDot(core, app, delta->userVector, fa[ii]->basis, va[ii-1]->basis);
      _braid_LRDeltaDotMat(core, app, delta->basis, fa[ii]->basis, va[ii-1]->basis);
   }

   _braid_Residual(core, level, index, braid_ASCaller_FASResidual, ustop, r);

   /* || short circuits, so (fa[ii] == NULL) is never evaluated if level == 0 is true */
   if ( (level == 0) || (fa[ii] == NULL) )
   {
      _braid_BaseSum(core, app, 0.0, r, -1.0, r);
   }
   else
   {
      if ( delta_correct )
      {
         /* delta correction */
         _braid_CoreFcn(core, sum)(app, 1.0, fa[ii]->userVector, 1.0, delta->userVector);

         /* tau correction */
         _braid_CoreFcn(core, sum)(app, 1.0, delta->userVector, -1.0, r->userVector);
         _braid_BaseSumBasis(core, app, 1.0, delta->basis, -1.0, r->basis);

         _braid_BaseFree(core, app, delta);
      }
      else
      {
         /* tau correction */
         _braid_BaseSum(core, app, 1.0, fa[ii], -1.0, r);
      }
   }

   return _braid_error_flag;
}

