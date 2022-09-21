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
 * F-Relax on level and interpolate to level-1
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FInterp(braid_Core  core,
               braid_Int   level)
{
   braid_App            app          = _braid_CoreElt(core, app);
   _braid_Grid        **grids        = _braid_CoreElt(core, grids);
   braid_AccessStatus   astatus      = (braid_AccessStatus)core;
   braid_Int            iter         = _braid_CoreElt(core, niter);
   braid_Int            access_level = _braid_CoreElt(core, access_level);
   braid_Int            nrefine      = _braid_CoreElt(core, nrefine);
   braid_Int            gupper       = _braid_CoreElt(core, gupper);
   braid_Int            ilower       = _braid_GridElt(grids[level], ilower);
   braid_Int            ncpoints     = _braid_GridElt(grids[level], ncpoints);
   braid_BaseVector    *va           = _braid_GridElt(grids[level], va);
   braid_Real          *ta           = _braid_GridElt(grids[level], ta);
   
   braid_Real         rnorm;
   braid_Int          f_level, f_cfactor, f_index;
   braid_BaseVector       f_u, f_e;

   braid_BaseVector       u, e;
   braid_Int          flo, fhi, fi, ci;
   braid_Int          interval;

   f_level   = level-1;
   f_cfactor = _braid_GridElt(grids[f_level], cfactor);

   _braid_GetRNorm(core, -1, &rnorm);
   
   _braid_UCommInitF(core, level);

   /**
    * Start from the right-most interval 
    *
    * First, generate the coarse-grid F-points through F-relaxation and
    * interpolate them to the fine grid, where they are C-points.  Second,
    * interpolate the coarse-grid C-points to the fine-grid.  The user-defined
    * spatial refinement (if set) is also called.  
    **/
   for (interval = ncpoints; interval > -1; interval--)
   {
      _braid_GetInterval(core, level, interval, &flo, &fhi, &ci);

      /* Relax and interpolate F-points, refining in space if needed */
      if (flo <= fhi)
      {
         _braid_UGetVector(core, level, flo-1, &u);
      }
      for (fi = flo; fi <= fhi; fi++)
      {
         _braid_Step(core, level, fi, braid_ASCaller_FInterp, NULL, u);
         _braid_USetVector(core, level, fi, u, 0);
         /* Allow user to process current vector */
         if( (access_level >= 3) )
         {
            _braid_AccessStatusInit(ta[fi-ilower], fi, rnorm, iter, level, nrefine, gupper,
                                    0, 0, braid_ASCaller_FInterp, u->basis, astatus);
            _braid_AccessVector(core, astatus, u);
         }
         e = va[fi-ilower];
         _braid_BaseSum(core, app,  1.0, u, -1.0, e);
         _braid_MapCoarseToFine(fi, f_cfactor, f_index);
         _braid_Refine(core, f_level, f_index, fi, e, &f_e);
         _braid_UGetVectorRef(core, f_level, f_index, &f_u);
         _braid_BaseSum(core, app,  1.0, f_e, 1.0, f_u);
         _braid_USetVectorRef(core, f_level, f_index, f_u);
         _braid_BaseFree(core, app,  f_e);
         /* Allow user to process current vector on the FINEST level*/
         if( (access_level >= 3) && (f_level == 0) )
         {
            _braid_AccessStatusInit(ta[fi-ilower], f_index, rnorm, iter, f_level, nrefine, gupper,
                                    0, 0, braid_ASCaller_FInterp, u->basis, astatus);
            _braid_AccessVector(core, astatus, f_u);
         }

      }
      if (flo <= fhi)
      {
         _braid_BaseFree(core, app,  u);
      }

      /* Interpolate C-points, refining in space if needed */
      if (ci > _braid_CoreElt(core, initiali))
      {
         _braid_UGetVectorRef(core, level, ci, &u);
         /* Allow user to process current C-point */
         if( (access_level >= 3) )
         {
            _braid_AccessStatusInit(ta[ci-ilower], ci, rnorm, iter, level, nrefine, gupper,
                                    0, 0, braid_ASCaller_FInterp, u->basis, astatus);
            _braid_AccessVector(core, astatus, u);
         }
         e = va[ci-ilower];
         _braid_BaseSum(core, app,  1.0, u, -1.0, e);
         _braid_MapCoarseToFine(ci, f_cfactor, f_index);
         _braid_Refine(core, f_level, f_index, ci, e, &f_e);
         _braid_UGetVectorRef(core, f_level, f_index, &f_u);
         _braid_BaseSum(core, app,  1.0, f_e, 1.0, f_u);
         _braid_USetVectorRef(core, f_level, f_index, f_u);
         _braid_BaseFree(core, app,  f_e);
         /* Allow user to process current C-point on the FINEST level*/
         if( (access_level >= 3) && (f_level == 0) )
         {
            _braid_AccessStatusInit(ta[ci-ilower], f_index, rnorm, iter, f_level, nrefine, gupper,
                                    0, 0, braid_ASCaller_FInterp, u->basis, astatus);
            _braid_AccessVector(core, astatus, f_u);
         }

      }
   }

   _braid_UCommWait(core, level);

   /* Clean up */
   _braid_GridClean(core, grids[level]);

   return _braid_error_flag;
}

