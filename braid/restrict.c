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
 * F-Relax on level and restrict to level+1
 *
 * At the bottom of FRestrict(), the coarse-grid right-hand side is computed and
 * stored in c_fa[].  The result stored in c_fa[] must be considered carefully
 * to maintain a fixed point nature of the algorithm.  The FAS coarse-grid
 * equation is
 *
 *   A_c(u_c) = R(f - A(u)) + A_c(R(u))
 *
 * where R() is restriction, and A() and A_c() are the fine and coarse
 * operators.  We desire our algorithm to be a fixed-point method, in that the
 * exact solution (from sequential time stepping) should yield a zero initial
 * residual and that the residual should stay zero as Braid iterates.  This is
 * actually a subtle point.  In this case, the solution u_c should be equivalent
 * to the solution u.  Since R(f - A(u)) = 0, this implies that
 *
 *   A_c(u_c) = A_c(R(u))
 *
 * This further implies that the initial guess provided to the implicit solve
 * routines in A_c() must be consistent on the right and left hand sides.
 *
 * The function (@ref _braid_GetUInit()) provides this functionality by
 * encapsulating all decisions about how to set this initial guess (called
 * ustop) The ustop values for A_c(R(u)) are set in (@ _braid_Residual) and for
 * A_c(u_c), they are set in (@ _braid_step).
 *
 * storage == -2 implies that (@ref _braid_GetUInit()) will always use the
 * previous time step value as the initial guess, even if you have better
 * information, which is the case at C-points on the fine grid, and at all
 * points on coarse-levels.
 *
 * storage == -1 implies that (@ref _braid_GetUInit()) will always use the best
 * information.  The value at the previous time step is used as the initial
 * guess only at F-points on the fine-grid.  The va[] values provide all initial
 * guess at F- and C-points on coarse-grids.
 *
 * storage == 0 is the same as -1, except that storage also exists at fine-grid
 * F-points.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FRestrict(braid_Core   core,
                 braid_Int    level)
{
   MPI_Comm              comm         = _braid_CoreElt(core, comm);
   braid_App             app          = _braid_CoreElt(core, app);
   _braid_Grid         **grids        = _braid_CoreElt(core, grids);
   braid_AccessStatus    astatus      = (braid_AccessStatus)core;
   braid_ObjectiveStatus ostatus      = (braid_ObjectiveStatus)core;
   braid_Int             iter         = _braid_CoreElt(core, niter);
   braid_Int             ichunk       = _braid_CoreElt(core, ichunk);
   braid_Int             print_level  = _braid_CoreElt(core, print_level);
   braid_Int             access_level = _braid_CoreElt(core, access_level);
   braid_Int             tnorm        = _braid_CoreElt(core, tnorm);
   braid_Real           *tnorm_a      = _braid_CoreElt(core, tnorm_a);
   braid_Int             nrefine      = _braid_CoreElt(core, nrefine);
   braid_Int             gupper       = _braid_CoreElt(core, gupper);
   braid_Int             cfactor      = _braid_GridElt(grids[level], cfactor);
   braid_Int             ncpoints     = _braid_GridElt(grids[level], ncpoints);
   braid_Real           *ta           = _braid_GridElt(grids[level], ta);
   braid_Int             f_ilower     = _braid_GridElt(grids[level], ilower);
   _braid_CommHandle    *recv_handle  = NULL;
   _braid_CommHandle    *send_handle  = NULL;

   braid_Int            c_level, c_ilower, c_iupper, c_index, c_i, c_ii;
   braid_BaseVector     c_u, *c_va, *c_fa;

   braid_BaseVector     u, r;
   braid_Int            interval, flo, fhi, fi, ci;
   braid_Real           rnorm, grnorm, rnorm_temp, rnm;

   c_level  = level+1;
   c_ilower = _braid_GridElt(grids[c_level], ilower);
   c_iupper = _braid_GridElt(grids[c_level], iupper);
   c_va     = _braid_GridElt(grids[c_level], va);
   c_fa     = _braid_GridElt(grids[c_level], fa);

   rnorm = 0.0;

   _braid_UCommInit(core, level);

   /* Start from the right-most interval.
    * 
    * Do an F-relax and then a C-relax.  These relaxations are needed to compute
    * the residual, which is needed for the coarse-grid right-hand-side and for
    * convergence checking on the finest grid.  This loop updates va and fa. */
   for (interval = ncpoints; interval > -1; interval--)
   {
      _braid_GetInterval(core, level, interval, &flo, &fhi, &ci);

      if (flo <= fhi)
      {
         _braid_UGetVector(core, level, flo-1, &r);
      }
      else if (ci > 0)
      {
         _braid_UGetVector(core, level, ci-1, &r);
      }

      /* F-relaxation */
      _braid_GetRNorm(core, -1, &rnm);
      for (fi = flo; fi <= fhi; fi++)
      {
         _braid_Step(core, level, fi, NULL, r);
         _braid_USetVector(core, level, fi, r, 0);
         
         /* Allow user to process current vector, note that r here is
          * temporarily holding the state vector */
         if( (access_level >= 3) )
         {
            _braid_AccessStatusInit(ta[fi-f_ilower], fi, ichunk, rnm, iter, level, nrefine, gupper,
                                    0, 0, braid_ASCaller_FRestrict, astatus);
            _braid_AccessVector(core, astatus, r);
         }

         /* Evaluate the user's local objective function at F-points on finest grid */
         if ( _braid_CoreElt(core, adjoint) && level == 0)
         {
            _braid_ObjectiveStatusInit(ta[fi-f_ilower], fi, ichunk, iter, level, nrefine, gupper, ostatus);
            _braid_AddToObjective(core, r, ostatus);
         }

      }

      /* Allow user to process current C-point */
      if( (access_level>= 3) && (ci > -1) )
      {
         _braid_AccessStatusInit(ta[ci-f_ilower], ci, ichunk, rnm, iter, level, nrefine, gupper,
                                 0, 0, braid_ASCaller_FRestrict, astatus);
         _braid_UGetVectorRef(core, level, ci, &u);
         _braid_AccessVector(core, astatus, u);
      }

      /* Evaluate the user's local objective function at CPoints on finest grid */
      if (_braid_CoreElt(core, adjoint) && level == 0 && (ci > -1) )
      {
         _braid_ObjectiveStatusInit(ta[ci-f_ilower], ci, ichunk, iter, level, nrefine, gupper, ostatus);
         _braid_UGetVectorRef(core, 0, ci, &u);
         _braid_AddToObjective(core, u, ostatus);
      }
         
      
      /* Compute residual and restrict */
      if (ci > 0)
      {
         /* Compute FAS residual */
         _braid_UGetVectorRef(core, level, ci, &u);
         _braid_FASResidual(core, level, ci, u, r);

         /* Compute rnorm (only on level 0) */
         if (level == 0)
         {
            _braid_BaseSpatialNorm(core, app,  r, &rnorm_temp);
            tnorm_a[interval] = rnorm_temp;       /* inf-norm uses tnorm_a */
            if(tnorm == 1) 
            {  
               rnorm += rnorm_temp;               /* one-norm combination */ 
            }
            else if(tnorm == 2)
            {  
               rnorm += (rnorm_temp*rnorm_temp);  /* two-norm combination */
            }
         }

         /* Restrict u and residual, coarsening in space if needed */
         _braid_MapFineToCoarse(ci, cfactor, c_index);
         _braid_Coarsen(core, c_level, ci, c_index, u, &c_va[c_index-c_ilower]);
         _braid_Coarsen(core, c_level, ci, c_index, r, &c_fa[c_index-c_ilower]);
      }
      else if (ci == 0)
      {
         /* Restrict initial condition, coarsening in space if needed */
         _braid_UGetVectorRef(core, level, 0, &u);
         _braid_Coarsen(core, c_level, 0, 0, u, &c_va[0]);
      }

      if ((flo <= fhi) || (ci > 0))
      {
         _braid_BaseFree(core, app,  r);
      }
   }
   _braid_UCommWait(core, level);

   /* If debug printing, print out tnorm_a for this interval. This
    * should show the serial propagation of the exact solution */
   if ((print_level > 2) && (level == 0) )
   {
      _braid_PrintSpatialNorms(core, tnorm_a, ncpoints);
   }

   /* Compute rnorm (only on level 0) */
   if (level == 0)
   {
      if(tnorm == 1)          /* one-norm reduction */
      {  
         MPI_Allreduce(&rnorm, &grnorm, 1, braid_MPI_REAL, MPI_SUM, comm);
      }
      else if(tnorm == 3)     /* inf-norm reduction */
      {  
         _braid_Max(tnorm_a, ncpoints, &rnorm); 
         MPI_Allreduce(&rnorm, &grnorm, 1, braid_MPI_REAL, MPI_MAX, comm);
      }
      else                    /* default two-norm reduction */
      {  
         MPI_Allreduce(&rnorm, &grnorm, 1, braid_MPI_REAL, MPI_SUM, comm);
         grnorm = sqrt(grnorm);
      }

      /* Store new rnorm */
      _braid_SetRNorm(core, -1, grnorm);
   }
   
   /* Now apply coarse residual to update fa values */

   /* Set initial guess on coarse level */
   _braid_InitGuess(core, c_level);

   /* Initialize update of c_va[-1] boundary */
   if (c_ilower <= c_iupper)
   {
      _braid_CommRecvInit(core, c_level, c_ilower-1, &c_va[-1],
                          &recv_handle);
      _braid_CommSendInit(core, c_level, c_iupper, c_va[c_iupper-c_ilower],
                          &send_handle);
   }

   /* Start with rightmost point */
   for (c_i = c_iupper; c_i >= c_ilower; c_i--)
   {
      if (c_i > 0)
      {
         c_ii = c_i - c_ilower;
         if (c_ii == 0)
         {
            /* Finalize update of c_va[-1] */
            _braid_CommWait(core, &recv_handle);
         }
         _braid_BaseClone(core, app,  c_va[c_ii-1], &c_u);
         _braid_Residual(core, c_level, c_i, c_va[c_ii], c_u);
         _braid_BaseSum(core, app,  1.0, c_u, 1.0, c_fa[c_ii]);
         _braid_BaseFree(core, app,  c_u);
      }
   }
   _braid_CommWait(core, &send_handle);
  
   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * ZTODO: Update access status
 *----------------------------------------------------------------------------*/

braid_Int
_braid_TriRestrict(braid_Core   core,
                   braid_Int    level)
{
   MPI_Comm             comm    = _braid_CoreElt(core, comm);
   braid_App            app     = _braid_CoreElt(core, app);
   _braid_Grid        **grids   = _braid_CoreElt(core, grids);
   braid_Int            clower  = _braid_GridElt(grids[level], clower);
   braid_Int            cupper  = _braid_GridElt(grids[level], cupper);
   braid_Int            cfactor = _braid_GridElt(grids[level], cfactor);

   braid_Int            nrequests;
   MPI_Request         *requests;
   MPI_Status          *statuses;
   void               **buffers;

   braid_Int            c_level, c_ilower, c_iupper, c_i;
   braid_BaseVector     c_u, *c_va, *c_fa;

   braid_Int            ci;
   braid_BaseVector     u, r, c_r;
   braid_Real           rnorm;

   c_level  = level+1;
   c_ilower = _braid_GridElt(grids[c_level], ilower);
   c_iupper = _braid_GridElt(grids[c_level], iupper);
   c_va     = _braid_GridElt(grids[c_level], va);
   c_fa     = _braid_GridElt(grids[c_level], fa);

   rnorm = 0.0;

   /* Compute residuals at C-points and restrict */

   /* Communicate u boundaries (don't worry about communication overlap yet) */
   _braid_TriCommInit(core, level, &nrequests, &requests, &statuses, &buffers);
   _braid_TriCommWait(core, level,  nrequests, &requests, &statuses, &buffers);

   for (ci = clower; ci <= cupper; ci += cfactor)
   {
      _braid_UGetVectorRef(core, level, ci, &u);
      _braid_TriResidual(core, level, ci, &r);

      /* Compute rnorm (only on level 0) */
      if (level == 0)
      {
         braid_Real  srnorm;

         _braid_BaseSpatialNorm(core, app, r, &srnorm);
         rnorm += (srnorm*srnorm);  /* two-norm */
      }

      /* Restrict u to coarse va and coarse u (this initializes coarse u)
       * Restrict residual to coarse fa
       * Coarsen in space if needed */
      _braid_MapFineToCoarse(ci, cfactor, c_i);
      _braid_Coarsen(core, c_level, ci, c_i, u, &c_va[c_i - c_ilower]);
      _braid_Coarsen(core, c_level, ci, c_i, r, &c_fa[c_i - c_ilower]);
      _braid_BaseFree(core, app, r);
      _braid_BaseClone(core, app, c_va[c_i - c_ilower], &c_u);
      _braid_USetVectorRef(core, c_level, c_i, c_u);
   }

   /* Compute rnorm (only on level 0) */
   if (level == 0)
   {
      braid_Real  grnorm;

      MPI_Allreduce(&rnorm, &grnorm, 1, braid_MPI_REAL, MPI_SUM, comm);
      grnorm = sqrt(grnorm);

      /* Store new rnorm */
      _braid_SetRNorm(core, -1, grnorm);
   }

   /* Finish FAS right-hand-side: A_c(u_c) = R(f - A(u)) + A_c(R(u))
    * Currently, the rhs holds R(A(u) - f) */

   /* Communicate c_u boundaries (don't worry about communication overlap yet) */
   _braid_TriCommInit(core, c_level, &nrequests, &requests, &statuses, &buffers);
   _braid_TriCommWait(core, c_level,  nrequests, &requests, &statuses, &buffers);

   for (c_i = c_ilower; c_i <= c_iupper; c_i++)
   {
      _braid_TriResidual(core, c_level, c_i, &c_r);
      _braid_BaseSum(core, app, 1.0, c_r, -1.0, c_fa[c_i - c_ilower]);
      _braid_BaseFree(core, app, c_r);
   }
  
   return _braid_error_flag;
}

