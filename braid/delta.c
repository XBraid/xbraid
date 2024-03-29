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
   braid_Real time = 0.;
   braid_Int rank = _braid_CoreElt(core, delta_rank);
   braid_Real *coords = _braid_TAlloc(braid_Real, rank);

   if (rank == 0)
   {
      return _braid_error_flag;
   }

   /* project u onto the basis */
   time = _braid_MPI_Wtime(core, 2);
   for (braid_Int i = 0; i < rank; i++)
   {
      _braid_CoreFcn(core, inner_prod)(app, basis->userVecs[i], u, &coords[i]);
   }
   _braid_CoreElt(core, timer_user_inner_prod) += _braid_MPI_Wtime(core, 2) - time;

   /* compute the dot product of action with the new coordinate vector */
   time = _braid_MPI_Wtime(core, 2);
   _braid_CoreFcn(core, sum)(app, coords[0], delta->userVecs[0], 0., u);
   for (braid_Int i = 1; i < rank; i++)
   {
      _braid_CoreFcn(core, sum)(app, coords[i], delta->userVecs[i], 1., u);
   }
   _braid_CoreElt(core, timer_user_sum) += _braid_MPI_Wtime(core, 2) - time;

   _braid_TFree(coords);

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
_braid_GramSchmidt(braid_Core   core,
                   braid_App    app,
                   braid_Basis  basis,
                   braid_Real  *exps)
{
   braid_Real innerprod_time = 0.; 
   braid_Real sum_time       = 0.;
   braid_Real time;

   for (braid_Int i = 0; i < basis->rank; i++)
   {
      braid_Real prod;
      /* subtract projections of the columns to the left */
      for (braid_Int j = 0; j < i; j++)
      {
         time = _braid_MPI_Wtime(core, 2);
         _braid_CoreFcn(core, inner_prod)(app, basis->userVecs[i], basis->userVecs[j], &prod);
         innerprod_time += _braid_MPI_Wtime(core, 2) - time;

         time = _braid_MPI_Wtime(core, 2);
         _braid_CoreFcn(core, sum)(app, -prod, basis->userVecs[j], 1., basis->userVecs[i]);
         sum_time += _braid_MPI_Wtime(core, 2) - time;
      }

      /* normalize this column */
      time = _braid_MPI_Wtime(core, 2);
      _braid_CoreFcn(core, inner_prod)(app, basis->userVecs[i], basis->userVecs[i], &prod);
      innerprod_time += _braid_MPI_Wtime(core, 2) - time;

      if (prod == 0)
      {
         _braid_Error(braid_ERROR_GENERIC, "Division by zero in _braid_GramSchmidt. Check that your basis vectors are not linearly dependent\n");
      }
      time = _braid_MPI_Wtime(core, 2);
      _braid_CoreFcn(core, sum)(app, 0., basis->userVecs[i], 1/sqrt(prod), basis->userVecs[i]);
      sum_time += _braid_MPI_Wtime(core, 2) - time;

      if (exps)
      {
         /* save local exponents (diagonals of R) */
         exps[i] = log(fabs(prod))/2;
      }
   }

   _braid_CoreElt(core, timer_user_inner_prod) += innerprod_time;
   _braid_CoreElt(core, timer_user_sum)       += sum_time;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_DeltaFeatureCheck(braid_Core core)
{

   braid_PtFcnResidual  residual  = _braid_CoreElt(core, residual);
   braid_PtFcnResidual  fullres   = _braid_CoreElt(core, full_rnorm_res);
   braid_PtFcnSCoarsen  scoarsen  = _braid_CoreElt(core, scoarsen);
   braid_PtFcnSRefine   srefine   = _braid_CoreElt(core, srefine);
   braid_Int            useshell  = _braid_CoreElt(core, useshell);
   braid_Int            trefine   = _braid_CoreElt(core, refine);
   braid_Int            adjoint   = _braid_CoreElt(core, adjoint);
   braid_Int            richardson= _braid_CoreElt(core, richardson);
   braid_Int err;
   char* err_char;

   err = 0;
   if ( residual != NULL ) 
   {
      err_char = "User-defined residual" ;
      err = 1;
   }
   if ( fullres  != NULL ) 
   {
      err_char = "User-defined full residual";
      err = 1;
   }
   if ( scoarsen != NULL ) 
   {
      err_char = "Spatial coarsening";
      err = 1;
   }
   if ( srefine  != NULL ) 
   {
      err_char = "Spatial refinement";
      err = 1;
   }
   if ( useshell ) 
   {
      err_char = "Shell-vector feature";
      err = 1;
   }
   if ( trefine )  
   {
      err_char = "Time refinement";
      err = 1;
   }
   if ( adjoint )  
   {
      err_char = "Adjoint feature";
      err = 1;
   }
   if ( richardson )  
   {
      err_char = "Richardson feature";
      err = 1;
   }

   // Print error message if needed 
   if ( err )
   {
      _braid_printf(" \n\n WARNING! %s is not yet supported for Delta correction (or at least not tested).\n", err_char); 
      _braid_printf("          If the code still runs, check solution carefully!\n\n\n", err_char); 
   }

   return _braid_error_flag;
}
 
