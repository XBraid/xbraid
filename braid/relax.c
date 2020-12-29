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
 * Do nu sweeps of F-then-C relaxation
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FCRelax(braid_Core  core,
               braid_Int   level)
{
   braid_App       app          = _braid_CoreElt(core, app);
   braid_Int      *nrels        = _braid_CoreElt(core, nrels);
   _braid_Grid   **grids        = _braid_CoreElt(core, grids);
   braid_Int       ncpoints     = _braid_GridElt(grids[level], ncpoints);
   braid_Int       done         = _braid_CoreElt(core, done);
   braid_Int       ntime        = _braid_CoreElt(core, ntime);
   braid_Int       storage      = _braid_CoreElt(core, storage);
   
   braid_AccessStatus   astatus      = (braid_AccessStatus)core;
   braid_Real          *ta           = _braid_GridElt(grids[level], ta);
   braid_Int            access_level = _braid_CoreElt(core, access_level);
   braid_Int            f_ilower     = _braid_GridElt(grids[level], ilower);
   braid_Int            iter         = _braid_CoreElt(core, niter);
   braid_Int            nrefine      = _braid_CoreElt(core, nrefine);
   braid_Int            gupper       = _braid_CoreElt(core, gupper);

   braid_BaseVector  u;
   braid_Int         flo, fhi, fi, ci;
   braid_Int         nu, nrelax, interval;
   braid_Real        rnm;

   nrelax  = nrels[level];

   for (nu = 0; nu < nrelax; nu++)
   {
      _braid_UCommInit(core, level);

      /* Start from the right-most interval */
      for (interval = ncpoints; interval > -1; interval--)
      {
         _braid_GetInterval(core, level, interval, &flo, &fhi, &ci);

         if (flo <= fhi)
         {
            _braid_UGetVector(core, level, flo-1, &u);
         }
         else if (ci > _braid_CoreElt(core, initiali))
         {
            _braid_UGetVector(core, level, ci-1, &u);
         }

         /* F-relaxation */
         _braid_GetRNorm(core, -1, &rnm);
         for (fi = flo; fi <= fhi; fi++)
         {
            _braid_Step(core, level, fi, NULL, u);
            _braid_USetVector(core, level, fi, u, 0);

            /* Allow user to process current vector */
            if( (access_level >= 3) || (done == 1) )
            {
               _braid_AccessStatusInit(ta[fi-f_ilower], fi, rnm, iter, level, nrefine, gupper,
                                       done, 0, braid_ASCaller_FCRelax, astatus);
               _braid_AccessVector(core, astatus, u);
            }

            /* If braid is finished, store last time step */
            if (done && fi == ntime && level == 0 && storage < 0)
            {
               if (_braid_GridElt(grids[level], ulast) != NULL)
               {
                  _braid_BaseFree(core, app, _braid_GridElt(grids[level], ulast));
                  _braid_GridElt(grids[level], ulast) = NULL;
               }
               _braid_BaseClone(core, app,  u, &(_braid_GridElt(grids[level], ulast)));
            }

         }

         /* C-relaxation */
         if (ci > _braid_CoreElt(core, initiali))
         {
            _braid_Step(core, level, ci, NULL, u);
            _braid_USetVector(core, level, ci, u, 1);
            
            /* Allow user to process current vector, noting that if done */
            if( (access_level >= 3) || (done == 1) )
            {
               _braid_AccessStatusInit(ta[ci-f_ilower], ci, rnm, iter, level, nrefine, gupper,
                                       done, 0, braid_ASCaller_FCRelax, astatus);
               _braid_AccessVector(core, astatus, u);
            }
         
            /* If braid is finished, store last time step */
            if (done && ci == ntime && level == 0 && storage < 0)
            {
               if (_braid_GridElt(grids[level], ulast) != NULL)
               {
                  _braid_BaseFree(core, app, _braid_GridElt(grids[level], ulast));
                  _braid_GridElt(grids[level], ulast) = NULL;
               }
               _braid_BaseClone(core, app,  u, &(_braid_GridElt(grids[level], ulast)));
            }
         }
         else if( (level == 0) && (ci == _braid_CoreElt(core, initiali)) && done)
         {
            printf("Tlingit\n");
            /* If final opportunity for user access, provide access to initial condition */
            _braid_UGetVector(core, level, ci, &u);
            _braid_AccessStatusInit(ta[ci-f_ilower], ci, rnm, iter, level, nrefine, gupper,
                                    done, 0, braid_ASCaller_FCRelax, astatus);
            _braid_AccessVector(core, astatus, u);
         }

         /* if ((flo <= fhi) && (interval == ncpoints)) */
         if ((flo <= fhi) && !(ci > _braid_CoreElt(core, initiali)))
         {
            _braid_BaseFree(core, app,  u);
         }
      }
      _braid_UCommWait(core, level);
   }

   return _braid_error_flag;
}

