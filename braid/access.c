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
 * Access to XBraid on grid level
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FAccess(braid_Core     core,
               braid_Int      level,
               braid_Int      done)
{
   braid_App              app          = _braid_CoreElt(core, app);
   _braid_Grid          **grids        = _braid_CoreElt(core, grids);
   braid_AccessStatus     astatus      = (braid_AccessStatus)core;
   braid_ObjectiveStatus  ostatus      = (braid_ObjectiveStatus)core;
   braid_Int              iter         = _braid_CoreElt(core, niter);
   braid_Int              ichunk       = _braid_CoreElt(core, ichunk);
   braid_Int              nrefine      = _braid_CoreElt(core, nrefine);
   braid_Int              gupper       = _braid_CoreElt(core, gupper);
   braid_Int              access_level = _braid_CoreElt(core, access_level);
   braid_Int              ncpoints     = _braid_GridElt(grids[level], ncpoints);
   braid_Real             *ta          = _braid_GridElt(grids[level], ta);
   braid_Int              ilower       = _braid_GridElt(grids[level], ilower);

   braid_Real        rnorm;
   braid_BaseVector  u;
   braid_Int         interval, flo, fhi, fi, ci;

   _braid_UCommInitF(core, level);
   
   _braid_GetRNorm(core, -1, &rnorm);

   /* Start from the right-most interval */
   for (interval = ncpoints; interval > -1; interval--)
   {
      _braid_GetInterval(core, level, interval, &flo, &fhi, &ci);

      /* Give access at F-points */
      if (flo <= fhi)
      {
         _braid_UGetVector(core, level, flo-1, &u);
      }
      for (fi = flo; fi <= fhi; fi++)
      {
         _braid_Step(core, level, fi, NULL, u);
         _braid_USetVector(core, level, fi, u, 0);

         if (access_level >= 1)
         {
            _braid_AccessStatusInit( ta[fi-ilower], fi, ichunk,  rnorm, iter, level, nrefine, gupper,
                                     done, 0, braid_ASCaller_FAccess, astatus);
            _braid_AccessVector(core, astatus, u);
         }

         /* If time-serial run: Evaluate the user's local objective function at F-points on finest grid */
         if ( _braid_CoreElt(core, adjoint) && 
              _braid_CoreElt(core, max_levels <=1) ) 
         {
            _braid_ObjectiveStatusInit(ta[fi-ilower], fi, ichunk, iter, level, nrefine, gupper, ostatus);
            _braid_AddToObjective(core, u, ostatus);
         }

         /* store last time step */
         if (_braid_CoreElt(core, storage) < 0 && level == 0 )
         {
            if (fi == _braid_CoreElt(core, ntime))
            {
              if (_braid_GridElt(grids[level], ulast) != NULL)
              {
                _braid_BaseFree(core, app, _braid_GridElt(grids[level], ulast));
                _braid_GridElt(grids[level], ulast) = NULL;
              }
              _braid_BaseClone(core, app,  u, &(_braid_GridElt(grids[level], ulast)));
            }
         }
      }

      if (flo <= fhi)
      {
         /* Free u */
         _braid_BaseFree(core, app,  u);
      }

      /* Give access at C-points */
      if ( ci > -1 )
      {
         _braid_UGetVectorRef(core, level, ci, &u);

         if ( access_level >= 1 )
         {
            _braid_AccessStatusInit( ta[ci-ilower], ci, ichunk, rnorm, iter, level, nrefine, gupper,
                                     done, 0, braid_ASCaller_FAccess, astatus);
            _braid_AccessVector(core, astatus, u);
         }

         /* If time-serial: Evaluate the user's local objective function at CPoints on finest grid */
         if ( _braid_CoreElt(core, adjoint)   && 
              _braid_CoreElt(core, max_levels <=1) ) 
         {
            _braid_ObjectiveStatusInit(ta[ci-ilower], ci, ichunk, iter, level, nrefine, gupper, ostatus);
            _braid_AddToObjective(core, u, ostatus);
         }

         /* store last time step */
         if (_braid_CoreElt(core, storage) < 0 && level == 0 )
         {
            if (ci == _braid_CoreElt(core, ntime))
            {
              if (_braid_GridElt(grids[level], ulast) != NULL)
              {
                _braid_BaseFree(core, app, _braid_GridElt(grids[level], ulast));
                _braid_GridElt(grids[level], ulast) = NULL;
              }
              _braid_BaseClone(core, app,  u, &(_braid_GridElt(grids[level], ulast)));
            }
         }
       }
   }
   _braid_UCommWait(core, level);

   return _braid_error_flag;
}

