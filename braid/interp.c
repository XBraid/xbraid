/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
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
#include "_util.h"

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
   braid_Int            ichunk       = _braid_CoreElt(core, ichunk);
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
         _braid_Step(core, level, fi, NULL, u);
         _braid_USetVector(core, level, fi, u, 0);
         /* Allow user to process current vector */
         if( (access_level >= 3) )
         {
            _braid_AccessStatusInit(ta[fi-ilower], fi, ichunk,  rnorm, iter, level, nrefine, gupper,
                                    0, 0, braid_ASCaller_FInterp, astatus);
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
      }
      if (flo <= fhi)
      {
         _braid_BaseFree(core, app,  u);
      }

      /* Interpolate C-points, refining in space if needed */
      if (ci > 0)
      {
         _braid_UGetVectorRef(core, level, ci, &u);
         /* Allow user to process current C-point */
         if( (access_level >= 3) )
         {
            _braid_AccessStatusInit(ta[ci-ilower], ci, ichunk, rnorm, iter, level, nrefine, gupper,
                                    0, 0, braid_ASCaller_FInterp, astatus);
            _braid_AccessVector(core, astatus, u);
         }
         e = va[ci-ilower];
         _braid_BaseSum(core, app, 1.0, u, -1.0, e);
         _braid_MapCoarseToFine(ci, f_cfactor, f_index);
         _braid_Refine(core, f_level, f_index, ci, e, &f_e);
         _braid_UGetVectorRef(core, f_level, f_index, &f_u);
         _braid_BaseSum(core, app, 1.0, f_e, 1.0, f_u);
         _braid_USetVectorRef(core, f_level, f_index, f_u);
         _braid_BaseFree(core, app, f_e);
      }
   }

   _braid_UCommWait(core, level);

   /* Clean up */
   _braid_GridClean(core, grids[level]);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_TriInterp(braid_Core   core,
                 braid_Int    level)
{
   braid_App            app     = _braid_CoreElt(core, app);
   _braid_Grid        **grids   = _braid_CoreElt(core, grids);
   braid_TriStatus      status  = (braid_TriStatus)core;
   braid_Int            ilower  = _braid_GridElt(grids[level], ilower);
   braid_Int            iupper  = _braid_GridElt(grids[level], iupper);
   braid_BaseVector    *va      = _braid_GridElt(grids[level], va);

   braid_Real           rnorm;
   braid_Int            f_level, f_cfactor, f_i;
   braid_BaseVector     f_u, f_e;

   braid_Int            i;
   braid_BaseVector     u, e;

   /* Update status (core) */
   _braid_GetRNorm(core, -1, &rnorm);
   _braid_StatusElt(status, rnorm) = rnorm;
   _braid_StatusElt(status, level) = level;

   f_level   = level-1;
   f_cfactor = _braid_GridElt(grids[f_level], cfactor);

#if 1
{
   /* Update u at C-points on the fine grid */
   for (i = ilower; i <= iupper; i++)
   {
      e = va[i-ilower];
      _braid_UGetVectorRef(core, level, i, &u);
      _braid_BaseSum(core, app, 1.0, u, -1.0, e);          // e = u - u0
      _braid_MapCoarseToFine(i, f_cfactor, f_i);
      _braid_Refine(core, f_level, f_i, i, e, &f_e);       // f_e = P_space e
      _braid_UGetVectorRef(core, f_level, f_i, &f_u);
      _braid_BaseSum(core, app, 1.0, f_e, 1.0, f_u);       // f_u = f_u + f_e
      _braid_USetVectorRef(core, f_level, f_i, f_u);
      _braid_BaseFree(core, app, f_e);
   }

   /* Update u at F-points with F-relaxation */
   _braid_TriFCFRelax(core, f_level, 0);
}
#else
{
   braid_BaseVector    *fa      = _braid_GridElt(grids[level], fa);
   braid_Int            f_ilower, f_iupper;

   f_ilower  = _braid_GridElt(grids[f_level], ilower);
   f_iupper  = _braid_GridElt(grids[f_level], iupper);

   /* Update u at C-points on the fine grid, store error corrections in fa temporarily */
   for (i = ilower; i <= iupper; i++)
   {
      braid_Real  *ta = _braid_GridElt(grids[level], ta);
      braid_Int    ii = i-ilower;

      /* Update status (core) */
      _braid_StatusElt(status, t)     = ta[ii];
      _braid_StatusElt(status, tprev) = ta[ii-1];
      _braid_StatusElt(status, tnext) = ta[ii+1];
      _braid_StatusElt(status, idx)   = i;

      e = va[i-ilower];
      _braid_UGetVectorRef(core, level, i, &u);
      _braid_BaseSum(core, app, 1.0, u, -1.0, e);          // e = u - u0
      _braid_MapCoarseToFine(i, f_cfactor, f_i);
      _braid_Refine(core, f_level, f_i, i, e, &f_e);       // f_e = P_space e
      _braid_UGetVectorRef(core, f_level, f_i, &f_u);
      _braid_BaseSum(core, app, 1.0, f_e, 1.0, f_u);       // f_u = f_u + f_e
      _braid_USetVectorRef(core, f_level, f_i, f_u);

      /* Store error correction in fa */
      _braid_BaseFree(core, app, fa[i-ilower]);
      fa[i-ilower] = f_e;
   }

   /* ZTODO: Communicate fa boundary values here */
   /* ZTODO: This only works for cfactor 2 at the moment.  It would be nice to
    * be able to just do F-relaxation with _braid_TriFCFRelax() on a grid with
    * the error corrections at C-points, zero at F-points initially, and using a
    * homogeneous problem with no FAS rhs. */

   /* Update u at F-points on the fine grid */
   _braid_StatusElt(status, level) = f_level;
   for (f_i = f_ilower; f_i <= f_iupper; f_i++)
   {
      if ( _braid_IsFPoint(f_i, f_cfactor) )
      {
         braid_Real  *f_ta = _braid_GridElt(grids[f_level], ta);
         braid_Int    f_ii = f_i-f_ilower;

         /* Update status (core) */
         _braid_StatusElt(status, t)     = f_ta[f_ii];
         _braid_StatusElt(status, tprev) = f_ta[f_ii-1];
         _braid_StatusElt(status, tnext) = f_ta[f_ii+1];
         _braid_StatusElt(status, idx)   = f_i;

         /* Create zero initial correction vector */
         _braid_UGetVectorRef(core, f_level, f_i, &f_u);
         _braid_BaseClone(core, app,  f_u, &f_e);
         _braid_BaseSum(core, app, 0.0, f_u, 0.0, f_e);    // f_e = 0

         /* Homogeneous solve at F-point with initial guess of zero */
         _braid_MapFineToCoarse(f_i-1, f_cfactor, i);
         _braid_BaseTriSolve(core, app, fa[i-ilower], fa[i+1-ilower], NULL, f_e, 1, status);

         _braid_BaseSum(core, app, 1.0, f_e, 1.0, f_u);    // f_u = f_u + f_e
         _braid_USetVectorRef(core, f_level, f_i, f_u);
         _braid_BaseFree(core, app, f_e);
      }
   }
}
#endif

   /* Clean up */
   _braid_GridClean(core, grids[level]);

   return _braid_error_flag;
}

