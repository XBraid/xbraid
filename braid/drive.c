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

/*--------------------------------------------------------------------------
 * Cycle state structure
 *--------------------------------------------------------------------------*/

typedef struct
{
   braid_Int  down;
   braid_Int  try_refine;
   braid_Int  fmglevel;
   braid_Int  fmg_Vcyc;
   FILE      *outfile;

} _braid_CycleState;

/*--------------------------------------------------------------------------
 * This is a locally scoped helper function for braid_Drive(), not a user
 * function.
 *--------------------------------------------------------------------------*/

braid_Int
_braid_DriveInitCycle(braid_Core          core,
                      _braid_CycleState  *cycle_ptr)
{
   braid_Int  myid       = _braid_CoreElt(core, myid_world);
   braid_Int  fmg        = _braid_CoreElt(core, fmg);
   braid_Int  nfmg       = _braid_CoreElt(core, nfmg);
   braid_Int  nlevels    = _braid_CoreElt(core, nlevels);
   braid_Int  io_level   = _braid_CoreElt(core, io_level);

   _braid_CycleState  cycle;

   cycle.down       = 1;
   cycle.try_refine = 0;
   cycle.fmglevel   = 0;
   cycle.fmg_Vcyc   = 0;

   if (fmg && (nfmg != 0))
   {
      cycle.fmglevel = nlevels-1;
   }

   /* Open cycle output file */
   if (myid == 0 && io_level>=1)
   {
      cycle.outfile = fopen("braid.out.cycle", "w");
   }

   *cycle_ptr = cycle;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * This is a locally scoped helper function for braid_Drive(), not a user
 * function.
 *
 * This routine determines the cycle direction (down or up) based on the current
 * grid level, iteration number, and cycle state.  The resulting cycle direction
 * is expected to produce three basic actions as follows:
 *
 *   Direction   Level                 Expected Action
 *   down        0...(nlevels-2)       relaxation/restriction
 *   up          1...(nlevels-1)       interpolation
 *   up          0                     refine or check convergence
 *--------------------------------------------------------------------------*/

braid_Int
_braid_DriveUpdateCycle(braid_Core          core,
                        braid_Int           level,
                        braid_Int           iter,
                        _braid_CycleState  *cycle_ptr)
{
   braid_Int      myid      = _braid_CoreElt(core, myid_world);
   braid_Int      fmg       = _braid_CoreElt(core, fmg);
   braid_Int      nfmg      = _braid_CoreElt(core, nfmg);
   braid_Int      nfmg_Vcyc = _braid_CoreElt(core, nfmg_Vcyc);
   braid_Int      nlevels   = _braid_CoreElt(core, nlevels);
   braid_Int      io_level  = _braid_CoreElt(core, io_level);
   _braid_CycleState  cycle = *cycle_ptr;
   braid_Real     rnorm;

   _braid_GetRNorm(core, -1, &rnorm);

   if (cycle.down)
   {
      /* Down cycle */

      /* If we are on the coarsest grid, go up */
      if (level == (nlevels-1))
      {
         cycle.down = 0;
      }
   }

   if (!cycle.down)
   {
      /* Up cycle */

      if (level > 0)
      {
         /* If we are not on the finest grid and we are doing F-cycles, make
          * sure we do nfmg_Vcyc V-cycles at each grid level */
         if (level < cycle.fmglevel)
         {
            cycle.fmg_Vcyc++;
            if ( cycle.fmg_Vcyc == nfmg_Vcyc )
            {
               cycle.fmg_Vcyc = 0;
               cycle.fmglevel--;
            }

            cycle.down = 1;
         }
      }
      else
      {
         /* If we are on the finest grid, first try to refine, then go down */
         if (!cycle.try_refine || nlevels == 1)
         {
            cycle.try_refine = 1;
         }
         else
         {
            /* Reset fmglevel if not done doing fmg */
            if (fmg && (nfmg != iter))
            {
               cycle.fmglevel = nlevels-1;
            }

            cycle.try_refine = 0;
            cycle.down = 1;
         }
      }
   }

   /* Print to cycle output file */
   if (myid == 0 && io_level>=1)
   {
      braid_Int            nrefine = _braid_CoreElt(core, nrefine);
      braid_Int            gupper  = _braid_CoreElt(core, gupper);
      braid_Real           tol     = _braid_CoreElt(core, tol);
      braid_Int            rtol    = _braid_CoreElt(core, rtol);
      braid_PtFcnResidual  fullres = _braid_CoreElt(core, full_rnorm_res);
      braid_Real           rnorm0;

      /* If using a relative tolerance, adjust tol */
      if (rtol)
      {
         if (fullres != NULL) {
            rnorm0 = _braid_CoreElt(core, full_rnorm0);
         }
         else {
            rnorm0 = _braid_CoreElt(core, rnorm0);
         }
         tol *= rnorm0;
      }

      _braid_ParFprintfFlush(cycle.outfile, myid, "%d %d %d %d %1.15e %1.15e\n",
                             level, nrefine, iter, gupper, rnorm, tol);
   }

   *cycle_ptr = cycle;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * This is a locally scoped helper function for braid_Drive(), not a user
 * function.
 *--------------------------------------------------------------------------*/

braid_Int
_braid_DriveEndCycle(braid_Core          core,
                     _braid_CycleState  *cycle_ptr)
{
   braid_Int  myid     = _braid_CoreElt(core, myid_world);
   braid_Int  io_level = _braid_CoreElt(core, io_level);

   _braid_CycleState  cycle = *cycle_ptr;

   /* Close cycle output file */
   if (myid == 0 && io_level>=1)
   {
      fclose(cycle.outfile);
   }

   *cycle_ptr = cycle;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * This is a locally scoped helper function for braid_Drive(), not a user
 * function.
 *--------------------------------------------------------------------------*/

braid_Int
_braid_DriveCheckConvergence(braid_Core  core,
                             braid_Int   iter,
                             braid_Int  *done_ptr)
{
   braid_Int            myid            = _braid_CoreElt(core, myid_world);
   braid_Real           tol             = _braid_CoreElt(core, tol);
   braid_Int            rtol            = _braid_CoreElt(core, rtol);
   braid_Int            max_iter        = _braid_CoreElt(core, max_iter);
   braid_PtFcnResidual  fullres         = _braid_CoreElt(core, full_rnorm_res);
   braid_Int            tight_fine_tolx = _braid_CoreElt(core, tight_fine_tolx);
   braid_Optim          optim           = _braid_CoreElt(core, optim);
   braid_Int            adjoint         = _braid_CoreElt(core, adjoint);
   braid_Int            obj_only        = _braid_CoreElt(core, obj_only);
   braid_Real           rnorm, rnorm0;
   braid_Real           rnorm_adj, rnorm0_adj;
   braid_Real           tol_adj, rtol_adj;

   braid_Int            done = *done_ptr;

   /* Use the full rnorm, if provided */
   if (fullres != NULL)
   {
      _braid_GetFullRNorm(core, -1, &rnorm);
      rnorm0 = _braid_CoreElt(core, full_rnorm0);
   }
   else
   {
      _braid_GetRNorm(core, -1, &rnorm);
      rnorm0 = _braid_CoreElt(core, rnorm0);
   }

   if ( adjoint )
   {
      /* Store state norm in the optimization structure */
      optim->rnorm  = rnorm;
      optim->rnorm0 = rnorm0;

      /* Get information from the optim structure */
      rnorm_adj  = optim->rnorm_adj;
      rnorm0_adj = optim->rnorm0_adj;
      tol_adj    = optim->tol_adj;
      rtol_adj   = optim->rtol_adj;
   }

   /* If using a relative tolerance, adjust tol */
   if (rtol)
   {
      tol *= rnorm0;
   }
   if ( adjoint )
   {
      if (rtol_adj)
      {
         tol_adj  *= rnorm0_adj;
      }
   }

   if ( (rnorm != braid_INVALID_RNORM) && (rnorm < tol) && (tight_fine_tolx == 1) )
   {
      done = 1;

      if ( adjoint && !obj_only )
      {
         /* Keep iterating, if adjoint not converged yet. */
         if ( ! (rnorm_adj < tol_adj) )
         {
            done = 0;
         }
      }
   }
   else if ( _braid_isnan(rnorm) )
   {
      if (myid == 0)
      {
         _braid_printf("  Braid: Iterations diverged.\n");
      }
      done = 1;
   }
   else if ( adjoint && _braid_isnan(rnorm_adj) )
   {
      if (myid == 0)
      {
         _braid_printf("  Braid: Adjoint iterations diverged.\n");
      }
      done = 1;
   }

   if (iter == max_iter-1 )
   {
      if (myid == 0)
      {
         _braid_printf("  Braid: Max. iterations reached.\n\n");
      }
      done = 1;
   }

   *done_ptr = done;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * This is a locally scoped helper function for braid_Drive(), not a user
 * function.
 *--------------------------------------------------------------------------*/

braid_Int
_braid_DrivePrintStatus(braid_Core  core,
                        braid_Int   level,
                        braid_Int   iter,
                        braid_Int   refined,
                        braid_Real  localtime)
{
   braid_Int            myid            = _braid_CoreElt(core, myid_world);
   braid_PtFcnResidual  fullres         = _braid_CoreElt(core, full_rnorm_res);
   braid_Int            rstopped        = _braid_CoreElt(core, rstopped);
   braid_Int            print_level     = _braid_CoreElt(core, print_level);
   braid_Optim          optim;
   braid_Real           rnorm, rnorm_prev, cfactor, wtime;
   braid_Real           rnorm_adj, objective;

   /* If my processor is not 0, or if print_level is not set high enough, return */
   if ((myid != 0) || (print_level < 1))
   {
      return _braid_error_flag;
   }

   wtime = MPI_Wtime() - localtime;

   if (_braid_CoreElt(core, adjoint))
   {
      optim     = _braid_CoreElt(core, optim);
      rnorm_adj = optim->rnorm_adj;
      objective = optim->objective;
      if (_braid_CoreElt(core, obj_only))
      {
         rnorm_adj = -1.0;
      }
   }

   _braid_GetRNorm(core, -1, &rnorm);
   _braid_GetRNorm(core, -2, &rnorm_prev);
   cfactor = 1.0;
   if (rnorm_prev != braid_INVALID_RNORM)
   {
      cfactor = rnorm / rnorm_prev;
   }
   if (rnorm != braid_INVALID_RNORM)
   {
      if (!_braid_CoreElt(core, adjoint))
      {
         _braid_printf("  Braid: || r_%d || = %1.6e, conv factor = %1.2e, wall time = %1.2e\n",
                       iter, rnorm, cfactor, wtime);
      }
      else
      {
         _braid_printf("  Braid: %3d  %1.6e  %1.6e  %1.8e\n", iter, rnorm, rnorm_adj, objective);
      }
   }
   else
   {
      _braid_printf("  Braid: || r_%d || not available, wall time = %1.2e\n", iter, wtime);
   }

   if (fullres != NULL)
   {
      _braid_GetFullRNorm(core, -1, &rnorm);
      _braid_GetFullRNorm(core, -2, &rnorm_prev);
      cfactor = 1.0;
      if (rnorm_prev != braid_INVALID_RNORM)
      {
         cfactor = rnorm / rnorm_prev;
      }
      if (rnorm != braid_INVALID_RNORM)
      {
         _braid_printf("  Braid: Full || r_%d || = %1.6e, conv factor = %1.2e\n",
                       iter, rnorm, cfactor);
      }
   }

   if (refined == 1)
   {
      _braid_printf("  Braid: Temporal refinement occurred, %d time steps\n",
                    _braid_CoreElt(core, gupper));
   }
   else if (refined == 2)
   {
      _braid_printf(" Braid: Spatial refinement occured, %d time steps\n",
                    _braid_CoreElt(core, gupper));
   }

   if ((rstopped > -1) && (rstopped == iter))
   {
      _braid_printf("  Braid: Temporal refinement stopped, %d time steps, %d refinements done\n",
                    _braid_CoreElt(core, gupper),
                    _braid_CoreElt(core, nrefine));
      _braid_printf("  Braid: Temporal refinement limits: time steps = %d, refinements = %d\n",
                    _braid_CoreElt(core, tpoints_cutoff),
                    _braid_CoreElt(core, max_refinements));
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
_braid_Drive(braid_Core  core, 
             braid_Real  localtime)
{
   braid_Int            skip            = _braid_CoreElt(core, skip);
   braid_Int            max_levels      = _braid_CoreElt(core, max_levels);
   braid_Int            incr_max_levels = _braid_CoreElt(core, incr_max_levels);
   braid_Int            access_level    = _braid_CoreElt(core, access_level);
   braid_PtFcnResidual  fullres         = _braid_CoreElt(core, full_rnorm_res);
   braid_Int            obj_only        = _braid_CoreElt(core, obj_only);
   braid_Int            adjoint         = _braid_CoreElt(core, adjoint);
   braid_Int            seq_soln        = _braid_CoreElt(core, seq_soln);
   braid_SyncStatus     sstatus         = (braid_SyncStatus)core;

   braid_Int     *nrels;
   braid_Int      nlevels;
   braid_Int      ilower, iupper;
   braid_Real     rnorm_adj;

   /* Cycle state variables */
   _braid_CycleState  cycle;
   braid_Int          iter, level, done, refined;

   /* Initialize cycle state */
   _braid_DriveInitCycle(core, &cycle);

   nlevels = _braid_CoreElt(core, nlevels);
   done  = 0;
   if (max_levels <= 1 && incr_max_levels == 0)
   {
      /* Just do sequential time marching */
      done = 1;
   }

   level = 0;
   if (skip)
   {
      /* Skip work on first down cycle */
      level = nlevels-1;
      _braid_CopyFineToCoarse(core);
   }

   iter = 0;
   _braid_CoreElt(core, niter) = iter;
   while (!done)
   {
      /* When there is just one grid level, do sequential time marching */
      /* RDF: Can this be done more efficiently?  It looks like this is needed
       * only to get the 'refine' factors, then the solve is recomputed in
       * either FRefine() or FAccess(). */
      if (nlevels == 1)
      {
         braid_Int  nrel0;
         nrels = _braid_CoreElt(core, nrels);
         nrel0 = nrels[0];
         nrels[0] = 1;
         _braid_FCRelax(core, 0);
         nrels[0] = nrel0;
         _braid_SetRNorm(core, -1, 0.0);
      }

      /* Update cycle state and direction based on level and iter */
      _braid_DriveUpdateCycle(core, level, iter, &cycle);

      if (cycle.down)
      {
         /* Down cycle */

         /* CF-relaxation */
         _braid_FCRelax(core, level);

         /* F-relax then restrict (note that FRestrict computes a new rnorm) */
         /* if adjoint: This computes the local objective function at each step on finest grid. */
         _braid_FRestrict(core, level);

         /* Compute full residual norm if requested */
         if ( (level == 0) &&  (fullres != NULL) )
         {
            braid_Real  full_rnorm;
            _braid_ComputeFullRNorm(core, level, &full_rnorm);
            _braid_SetFullRNorm(core, -1, full_rnorm);
         }

         level++;
      }
      else
      {
         /* Up cycle */

         if (level > 0)
         {
            /* Periodic coarsest grid solve */
            if ( (level == (nlevels-1)) && _braid_CoreElt(core, periodic) )
            {
               braid_Int  nrel;
               nrels = _braid_CoreElt(core, nrels);
               nrel  = nrels[level];
               nrels[level] = nrels[max_levels-1];
               _braid_FCRelax(core, level);
               nrels[level] = nrel;
            }

            /* F-relax then interpolate */
            _braid_FInterp(core, level);

            level--;
         }
         else
         {
            _braid_SyncStatusInit(iter, level, _braid_CoreElt(core, nrefine),
                                  _braid_CoreElt(core, gupper), done,
                                  braid_ASCaller_Drive_TopCycle, sstatus);
            _braid_Sync(core, sstatus);

            // Output the solution at the end of each cycle
            // Copy the rfactors because the call to FAccess will modify them
            if (access_level >= 2)
            {
               _braid_Grid **grids = _braid_CoreElt(core, grids);
               ilower = _braid_GridElt(grids[0], ilower);
               iupper = _braid_GridElt(grids[0], iupper);
               braid_Int *saved_rfactors = _braid_CTAlloc(braid_Int,iupper-ilower+2);
               braid_Int *rfactors       = _braid_CoreElt(core, rfactors);
               int ii,i;
               for (i=ilower; i<=iupper+1; i++)
               {
                  ii=i-ilower;
                  saved_rfactors[ii]=rfactors[ii];
               }
               _braid_FAccess(core, 0, 0);
               for (i=ilower; i<=iupper+1; i++)
               {
                  ii=i-ilower;
                  rfactors[ii]=saved_rfactors[ii];
               }
               _braid_TFree(saved_rfactors);
            }

            /* Finest grid - refine grid if desired */
            _braid_FRefine(core, &refined);
            nlevels = _braid_CoreElt(core, nlevels);

            // If we are done refining and doing a fixed point test,
            // then compute the sequential solution on the finest grid
            braid_Int rstopped = _braid_CoreElt(core, rstopped);
            if( (rstopped > -1) && (rstopped == iter) && (seq_soln == 1) )
            {
               _braid_InitGuess(core, 0);
            }

            if ( adjoint )
            {
               /* Compute the objective function */
               _braid_EvalObjective(core);

               /* Compute differentiated objective function */
               _braid_EvalObjective_diff(core);

               /* Set the adjoint seed at coarse points on level 0 */
               _braid_TapeSetSeed(core);

               /* Evaluate (and clear) the action tape */
               _braid_TapeEvaluate(core);

               /* Update adjoints and compute residual norm */
               _braid_UpdateAdjoint(core, &rnorm_adj);
               _braid_SetRNormAdjoint(core, iter, rnorm_adj);
            }

            /* Print current status */
            _braid_DrivePrintStatus(core, level, iter, refined, localtime);

            /* If no refinement was done, check for convergence */
            if (!refined)
            {
               /* Check convergence */
               _braid_DriveCheckConvergence(core, iter, &done);
            }

            if ( adjoint)
            {
               /* Prepare for the next iteration */
               _braid_CoreElt(core, optim)->sum_user_obj = 0.0;
               _braid_CoreElt(core, optim)->f_bar        = 0.0;

               if (!done && !obj_only)
               {
                  _braid_CoreFcn(core, reset_gradient)(_braid_CoreElt(core, app));
               }

               /* Reset the pointer to input variables */
               _braid_TapeResetInput(core);
            }

            /* Increase MGRIT iteration counter */
            if (!refined)
            {
               iter++;
               _braid_CoreElt(core, niter) = iter;
            }
         }
      }
   }

   /* By default, set the final residual norm to be the same as the previous */
   {
      braid_Real  rnorm;
      _braid_GetRNorm(core, -2, &rnorm);
      _braid_SetRNorm(core, -1, rnorm);
   }

   /* Compute final full residual norms if requested */
   if (fullres != NULL)
   {
      braid_Real  full_rnorm;
      _braid_ComputeFullRNorm(core, level, &full_rnorm);
      _braid_SetFullRNorm(core, -1, full_rnorm);
      /* JBS: Ben S wanted a final rnorm, we should move this final residual
       * computation to FAccess to save work */
      _braid_FRestrict(core, level);
   }

   /* Allow final access to Braid by carrying out an F-relax to generate points */
   /* Record it only if sequential time stepping */
   if (max_levels > 1)
   {
      _braid_CoreElt(core, record) = 0;
   }
   _braid_FAccess(core, 0, 1);

   /* If sequential time-marching, evaluate the tape */
   if ( adjoint && max_levels <= 1 )
   {
      /* Compute the objective function */
      _braid_EvalObjective(core);

      /* Compute differentiated objective function */
      _braid_EvalObjective_diff(core);

      /* Evaluate (and clear) the action tape */
      _braid_TapeEvaluate(core);
   }

   /* End cycle */
   _braid_DriveEndCycle(core, &cycle);

   return _braid_error_flag;
}

