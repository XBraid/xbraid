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

braid_Int
_braid_LRDeltaDot(braid_Core core,
                  braid_App app,
                  braid_Vector u,
                  braid_Basis delta,
                  braid_Basis basis)
{
   braid_Int rank = _braid_CoreElt(core, delta_rank);
   braid_Real *coords = _braid_TAlloc(braid_Real, rank);

   if (rank == 0)
   {
      return _braid_error_flag;
   }

   /* project u onto the basis */
   for (braid_Int i = 0; i < rank; i++)
   {
      _braid_CoreFcn(core, inner_prod)(app, basis->userVecs[i], u, &coords[i]);
   }

   /* compute the dot product of action with the new coordinate vector */
   _braid_CoreFcn(core, sum)(app, coords[0], delta->userVecs[0], 0., u);
   for (braid_Int i = 1; i < rank; i++)
   {
      _braid_CoreFcn(core, sum)(app, coords[i], delta->userVecs[i], 1., u);
   }

   return _braid_error_flag;
}

braid_Int
_braid_LRDeltaDotMat(braid_Core core,
                     braid_App app,
                     braid_Basis psi,
                     braid_Basis delta,
                     braid_Basis basis)
{
   braid_Int rank = _braid_CoreElt(core, delta_rank);
   for (braid_Int i = 0; i < rank; i++)
   {
      _braid_LRDeltaDot(core, app, psi->userVecs[i], delta, basis);
   }

   return _braid_error_flag;
}

braid_Int
_braid_Normalize(braid_Core    core,
                 braid_App     app,
                 braid_Vector  u,
                 braid_Real   *norm_ptr)
{
   braid_Real norm;
   _braid_CoreFcn(core, inner_prod)(app, u, u, &norm);
   _braid_CoreFcn(core, sum)(app, 0., u, 1/sqrt(norm), u);

   if (norm_ptr)
   {
      *norm_ptr = norm;
   } 

   return _braid_error_flag;
}

braid_Int
_braid_GramSchmidt(braid_Core   core,
                   braid_App    app,
                   braid_Basis  basis,
                   braid_Real  *exps)
{
   for (braid_Int i = 0; i < basis->rank; i++)
   {
      braid_Real prod;
      /* subtract projections of the columns to the left */
      for (braid_Int j = 0; j < i; j++)
      {
         _braid_CoreFcn(core, inner_prod)(app, basis->userVecs[i], basis->userVecs[j], &prod);
         _braid_CoreFcn(core, sum)(app, -prod, basis->userVecs[j], 1., basis->userVecs[i]);
      }

      /* normalize this column */
      _braid_CoreFcn(core, inner_prod)(app, basis->userVecs[i], basis->userVecs[i], &prod);
      _braid_CoreFcn(core, sum)(app, 0., basis->userVecs[i], 1/sqrt(prod), basis->userVecs[i]);

      if (exps)
      {
         /* save local exponents (diagonals of R) */
         exps[i] = log(fabs(prod))/2;
      }
   }

   return _braid_error_flag;
}