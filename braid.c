/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. Email xbraid-support@llnl.gov for support.
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

/** \file braid.c
 * \brief Source code for user interface routines.  See braid.h for more information.
 *
 */

#include "_braid.h"
#include "braid_defs.h"
#include "_braid_tape.h"
#include "_util.h"
#include "_braid_base.h"

#ifndef DEBUG
#define DEBUG 0
#endif

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
   braid_Real           rnorm, rnorm0;

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

   /* If using a relative tolerance, adjust tol */
   if (rtol)
   {
      tol *= rnorm0;
   }

   if ( (rnorm != braid_INVALID_RNORM) && (rnorm < tol) && (tight_fine_tolx == 1) )
   {
      done = 1;
   } 
   else if (iter == max_iter-1)
   {
      done = 1;
   } 
   else if ( braid_isnan(rnorm) )
   {
      if (myid == 0)
      {
         _braid_printf("  Iterations diverged.\n");
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
   braid_Real           rnorm, rnorm_prev, cfactor, wtime;

   /* If my processor is not 0, or if print_level is not set high enough, return */
   if ((myid != 0) || (print_level < 1))
   {
      return _braid_error_flag;
   }

   wtime = MPI_Wtime() - localtime;

   _braid_GetRNorm(core, -1, &rnorm);
   _braid_GetRNorm(core, -2, &rnorm_prev);
   cfactor = 1.0;
   if (rnorm_prev != braid_INVALID_RNORM)
   {
      cfactor = rnorm / rnorm_prev;
   }
   if (rnorm != braid_INVALID_RNORM)
   {
       _braid_printf("  Braid: || r_%d || = %1.6e, conv factor = %1.2e, wall time = %1.2e\n",
                       iter, rnorm, cfactor, wtime);
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
braid_Drive(braid_Core  core)
{
   MPI_Comm             comm_world      = _braid_CoreElt(core, comm_world);
   braid_Int            myid            = _braid_CoreElt(core, myid_world);
   braid_Real           tstart          = _braid_CoreElt(core, tstart);
   braid_Real           tstop           = _braid_CoreElt(core, tstop);
   braid_Int            ntime           = _braid_CoreElt(core, ntime);
   braid_Int            skip            = _braid_CoreElt(core, skip);
   braid_Int            max_levels      = _braid_CoreElt(core, max_levels);
   braid_Int            print_level     = _braid_CoreElt(core, print_level);
   braid_Int            access_level    = _braid_CoreElt(core, access_level);
   braid_App            app             = _braid_CoreElt(core, app);
   braid_PtFcnResidual  fullres         = _braid_CoreElt(core, full_rnorm_res);

   braid_Int     *nrels, nrel0;
   braid_Int      nlevels;
   braid_Int      ilower, iupper, i;
   braid_Real    *ta;
   _braid_Grid   *grid;
   braid_Real     localtime, globaltime;
   braid_Real     localobjective, globalobjective;

   /* Cycle state variables */
   _braid_CycleState  cycle;
   braid_Int          iter, level, done, refined;

   if (myid == 0)
   { 
      _braid_printf("  Braid: Begin simulation, %d time steps\n\n",
                    _braid_CoreElt(core, gupper));
   }

   /* Start timer */
   localtime = MPI_Wtime();

   /* Create fine grid */
   _braid_GetDistribution(core, &ilower, &iupper);
   _braid_GridInit(core, 0, ilower, iupper, &grid);

   /* Set t values */
   ta = _braid_GridElt(grid, ta);
   if ( _braid_CoreElt(core, tgrid) != NULL )
   {
      /* Call the user's time grid routine */
      _braid_BaseTimeGrid(core, app, ta, &ilower, &iupper);
   }
   else
   {
      for (i = ilower; i <= iupper; i++)
      {
         ta[i-ilower] = tstart + (((braid_Real)i)/ntime)*(tstop-tstart);
      }
   }

   /* Create a grid hierarchy */
   _braid_InitHierarchy(core, grid, 0);
   nlevels = _braid_CoreElt(core, nlevels);

   /* Set initial values */
   _braid_InitGuess(core, 0);

   /* Initialize cycle state */
   _braid_DriveInitCycle(core, &cycle);
   
   done  = 0;
   if (max_levels <= 1)
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
            /* F-relax then interpolate */
            _braid_FInterp(core, level);

            level--;
         }
         else
         {

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

            /* Finest grid - refine grid if desired, else check convergence */
            _braid_FRefine(core, &refined);
            nlevels = _braid_CoreElt(core, nlevels);

            /* Print current status */
            _braid_DrivePrintStatus(core, level, iter, refined, localtime);

            /* If no refinement was done, check for convergence */
            if (!refined)
            {
               /* Check convergence */
               _braid_DriveCheckConvergence(core, iter, &done);

               iter++;
               _braid_CoreElt(core, niter) = iter;
            }

            /*--- Evaluate the adjoint sensitivities after one cycle ---*/

            // printf("Eval Adjoint. iter %d level %d\n", iter, level);
            if (_braid_CoreElt(core,adjoint))
            {
               /* Compute the time-averaged objective function. */
               localobjective = _braid_CoreElt(core, optim)->objective;
               MPI_Allreduce(&localobjective, &globalobjective, 1, braid_MPI_REAL, MPI_SUM, comm_world);
               _braid_CoreElt(core, optim)->objective = globalobjective / ( ntime + 1 );
               printf("  Objective = %1.14e\n", _braid_CoreElt(core, optim)->objective);

               /* Evaluate (and clear) the action tape */
               _braid_TapeEvaluate(core);

               /* Reset the objective function for the next iteration */
               _braid_CoreElt(core, optim)->objective = 0.0;

               /* Reset the gradient */
      
               /* Stop iterating */
               done = 1;
               _braid_CoreElt(core, done) = 1;
      
            } /* End of Adjoint */

         }
      }
   }

   /* Stop adjoint recording */
   printf("\nStop recording\n\n");
   _braid_CoreElt(core, record) = 0;

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

   /* If sequential time-marching, record it! */
   if (max_levels <= 1)
   {
      _braid_CoreElt(core, record) = 1; 
   }

   /* Allow final access to Braid by carrying out an F-relax to generate points */
   _braid_FAccess(core, 0, 1);


   /* If sequential time-marching, evaluate the tape */
   if (max_levels <= 1 && _braid_CoreElt(core, adjoint))
   {
      /* Compute the time-averaged objective function. */
      localobjective = _braid_CoreElt(core, optim)->objective;
      MPI_Allreduce(&localobjective, &globalobjective, 1, braid_MPI_REAL, MPI_SUM, comm_world);
      _braid_CoreElt(core, optim)->objective = globalobjective / ( ntime + 1 );
      printf("  Objective = %1.14e\n", _braid_CoreElt(core, optim)->objective);
      /* Evaluate (and clear) the action tape */
      _braid_TapeEvaluate(core);
   }


   /* End cycle */
   _braid_DriveEndCycle(core, &cycle);

   /* Stop timer */
   localtime = MPI_Wtime() - localtime;
   MPI_Allreduce(&localtime, &globaltime, 1, braid_MPI_REAL, MPI_MAX, comm_world);
   _braid_CoreElt(core, localtime)  = localtime;
   _braid_CoreElt(core, globaltime) = globaltime;

   /* Print statistics for this run */
   if ( (print_level > 0) && (myid == 0) )
   {
      braid_PrintStats(core);
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_Init(MPI_Comm               comm_world,
           MPI_Comm               comm,
           braid_Real             tstart,
           braid_Real             tstop,
           braid_Int              ntime,
           braid_App              app,
           braid_PtFcnStep        step,
           braid_PtFcnInit        init,
           braid_PtFcnClone       clone,
           braid_PtFcnFree        free,
           braid_PtFcnSum         sum,
           braid_PtFcnSpatialNorm spatialnorm,
           braid_PtFcnAccess      access,
           braid_PtFcnBufSize     bufsize,
           braid_PtFcnBufPack     bufpack,
           braid_PtFcnBufUnpack   bufunpack,
           braid_Core            *core_ptr)
{
   _braid_Core           *core;

   /* Braid default values */
   braid_Int              cfdefault       = 2;              /* Default coarsening factor */
   braid_Int              nrdefault       = 1;              /* Default number of FC sweeps on each level */
   braid_Int              fmg             = 0;              /* Default fmg (0 is off) */
   braid_Int              nfmg            = -1;             /* Default fmg cycles is -1, indicating all fmg-cycles (if fmg=1) */
   braid_Int              nfmg_Vcyc       = 1;              /* Default num V-cycles at each fmg level is 1 */
   braid_Int              max_iter        = 100;            /* Default max_iter */
   braid_Int              max_levels      = 30;             /* Default max_levels */
   braid_Int              min_coarse      = 2;              /* Default min_coarse */
   braid_Int              seq_soln        = 0;              /* Default initial guess is from user's Init() function */
   braid_Int              print_level     = 1;              /* Default print level */
   braid_Int              io_level        = 1;              /* Default output-to-file level */
   braid_Int              access_level    = 1;              /* Default access level */
   braid_Int              tnorm           = 2;              /* Default temporal norm */
   braid_Real             tol             = 1.0e-09;        /* Default absolute tolerance */
   braid_Int              rtol            = 1;              /* Use relative tolerance */
   braid_Int              skip            = 1;              /* Default skip value, skips all work on first down-cycle */
   braid_Int              max_refinements = 200;            /* Maximum number of F-refinements */
   braid_Int              tpoints_cutoff  = braid_Int_Max;  /* Maximum number of time steps, controls FRefine()*/ 
   braid_Int              adjoint         = 0;              /* Default adjoint run: Turned off */
   braid_Int              record          = 0;              /* Default action recording: Turned off */
   braid_Int              verbose         = 0;              /* Default verbosity Turned off */

   braid_Int              myid_world,  myid;

   MPI_Comm_rank(comm_world, &myid_world);
   MPI_Comm_rank(comm, &myid);

   core = _braid_CTAlloc(_braid_Core, 1);

   _braid_CoreElt(core, comm_world)      = comm_world;
   _braid_CoreElt(core, comm)            = comm;
   _braid_CoreElt(core, myid_world)      = myid_world;
   _braid_CoreElt(core, myid)            = myid;
   _braid_CoreElt(core, tstart)          = tstart;
   _braid_CoreElt(core, tstop)           = tstop;
   _braid_CoreElt(core, ntime)           = ntime;
   _braid_CoreElt(core, app)             = app;

   _braid_CoreElt(core, step)            = step;
   _braid_CoreElt(core, init)            = init;
   _braid_CoreElt(core, sinit)           = NULL;
   _braid_CoreElt(core, clone)           = clone;
   _braid_CoreElt(core, sclone)          = NULL;
   _braid_CoreElt(core, free)            = free;
   _braid_CoreElt(core, sfree)           = NULL;
   _braid_CoreElt(core, sum)             = sum;
   _braid_CoreElt(core, spatialnorm)     = spatialnorm;
   _braid_CoreElt(core, access)          = access;
   _braid_CoreElt(core, bufsize)         = bufsize;
   _braid_CoreElt(core, bufpack)         = bufpack;
   _braid_CoreElt(core, bufunpack)       = bufunpack;
   _braid_CoreElt(core, residual)        = NULL;
   _braid_CoreElt(core, scoarsen)        = NULL;
   _braid_CoreElt(core, srefine)         = NULL;
   _braid_CoreElt(core, tgrid)           = NULL;

   _braid_CoreElt(core, access_level)    = access_level;
   _braid_CoreElt(core, tnorm)           = tnorm;
   _braid_CoreElt(core, print_level)     = print_level;
   _braid_CoreElt(core, io_level)        = io_level;
   _braid_CoreElt(core, max_levels)      = 0; /* Set with SetMaxLevels() below */
   _braid_CoreElt(core, min_coarse)      = min_coarse;
   _braid_CoreElt(core, seq_soln)        = seq_soln;
   _braid_CoreElt(core, tol)             = tol;
   _braid_CoreElt(core, rtol)            = rtol;

   _braid_CoreElt(core, nrels)           = NULL; /* Set with SetMaxLevels() below */
   _braid_CoreElt(core, nrdefault)       = nrdefault;

   _braid_CoreElt(core, cfactors)        = NULL; /* Set with SetMaxLevels() below */
   _braid_CoreElt(core, cfdefault)       = cfdefault;

   _braid_CoreElt(core, max_iter)        = 0; /* Set with SetMaxIter() below */
   _braid_CoreElt(core, niter)           = 0;
   _braid_CoreElt(core, fmg)             = fmg;
   _braid_CoreElt(core, nfmg)            = nfmg;
   _braid_CoreElt(core, nfmg_Vcyc)       = nfmg_Vcyc;

   _braid_CoreElt(core, storage)         = -1;            /* only store C-points */
   _braid_CoreElt(core, useshell)         = 0;

   _braid_CoreElt(core, gupper)          = ntime;

   _braid_CoreElt(core, refine)          = 0;  /* Time refinement off by default */
   _braid_CoreElt(core, rfactors)        = NULL;
   _braid_CoreElt(core, r_space)         = 0;
   _braid_CoreElt(core, rstopped)        = -1;
   _braid_CoreElt(core, nrefine)         = 0;
   _braid_CoreElt(core, max_refinements) = max_refinements;
   _braid_CoreElt(core, tpoints_cutoff)  = tpoints_cutoff;

   _braid_CoreElt(core, nlevels)         = 0;
   _braid_CoreElt(core, grids)           = NULL; /* Set with SetMaxLevels() below */

   _braid_CoreElt(core, skip)            = skip;

   _braid_CoreElt(core, adjoint)         = adjoint;
   _braid_CoreElt(core, record)          = record;
   _braid_CoreElt(core, verbose)         = verbose;
   _braid_CoreElt(core, actiontape)      = NULL;
   _braid_CoreElt(core, primaltape)      = NULL;
   _braid_CoreElt(core, adjointtape)     = NULL;
   _braid_CoreElt(core, optim)           = NULL;
   _braid_CoreElt(core, step_diff)       = NULL;
   _braid_CoreElt(core, objT_diff)       = NULL;
   _braid_CoreElt(core, objectiveT)      = NULL;

   
   /* Residual history and accuracy tracking for StepStatus*/
   _braid_CoreElt(core, rnorm0)              = braid_INVALID_RNORM;
   _braid_CoreElt(core, rnorms)              = NULL; /* Set with SetMaxIter() below */
   _braid_CoreElt(core, full_rnorm_res)      = NULL;
   _braid_CoreElt(core, full_rnorm0)         = braid_INVALID_RNORM;
   _braid_CoreElt(core, full_rnorms)         = NULL; /* Set with SetMaxIter() below */
   _braid_CoreElt(core, old_fine_tolx)       = -1.0;
   _braid_CoreElt(core, tight_fine_tolx)     = 1;

   braid_SetMaxLevels(core, max_levels);
   braid_SetMaxIter(core, max_iter);

   *core_ptr = core;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/
braid_Int
braid_Init_Adjoint(braid_PtFcnObjectiveT      objectiveT,
                   braid_PtFcnStepDiff        step_diff, 
                   braid_PtFcnObjectiveTDiff  objT_diff,
                   braid_Core                 *core_ptr)
{
  if( _braid_CoreElt(*core_ptr, refine) )
  {
    printf("ERROR: Adjoint of FRefine not supported yet!\n");
    exit(1);
  }

   /* Set the adjoint flag */ 
   _braid_CoreElt(*core_ptr, adjoint) = 1;
   _braid_CoreElt(*core_ptr, record) = 1;

   /* Turn on adjoint verbosity */
   _braid_CoreElt(*core_ptr, verbose) = 0;

  /* Turn off skip on downcycle */
  _braid_CoreElt(*core_ptr, skip) = 0;

   /* Initialize the tapes */
   _braid_TapeInit( _braid_CoreElt(*core_ptr,actiontape) );
   _braid_TapeInit( _braid_CoreElt(*core_ptr,primaltape) );
   _braid_TapeInit( _braid_CoreElt(*core_ptr,adjointtape) );

   /* Initialize the optimization structure */
   braid_Optim optim;
   _braid_OptimInit( core_ptr, &optim);
   _braid_CoreElt(*core_ptr, optim) = optim; 

   /* Additional user functions */
   _braid_CoreElt(*core_ptr, objectiveT)     = objectiveT;

   /* Differentiated functions */
   _braid_CoreElt(*core_ptr, step_diff)   = step_diff;
   _braid_CoreElt(*core_ptr, objT_diff)   = objT_diff;
  

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_Destroy(braid_Core  core)
{
   if (core)
   {
      braid_Int               nlevels    = _braid_CoreElt(core, nlevels);
      _braid_Grid           **grids      = _braid_CoreElt(core, grids);
      braid_Int               level;

      _braid_TFree(_braid_CoreElt(core, nrels));
      _braid_TFree(_braid_CoreElt(core, rnorms));
      _braid_TFree(_braid_CoreElt(core, full_rnorms));
      _braid_TFree(_braid_CoreElt(core, cfactors));
      _braid_TFree(_braid_CoreElt(core, rfactors));
      _braid_TFree(_braid_CoreElt(core, tnorm_a));
      _braid_TFree(_braid_CoreElt(core, optim));
      
      for (level = 0; level < nlevels; level++)
      {
         _braid_GridDestroy(core, grids[level]);
      }
      _braid_TFree(grids);

      _braid_TFree(core);
   }

   if (_braid_printfile != NULL)
   {
      fclose(_braid_printfile);
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_PrintStats(braid_Core  core)
{
   braid_Int     myid          = _braid_CoreElt(core, myid_world);
   braid_Real    tstart        = _braid_CoreElt(core, tstart);
   braid_Real    tstop         = _braid_CoreElt(core, tstop);
   braid_Int     gupper        = _braid_CoreElt(core, gupper);
   braid_Int     max_levels    = _braid_CoreElt(core, max_levels);
   braid_Int     min_coarse    = _braid_CoreElt(core, min_coarse);
   braid_Int     seq_soln      = _braid_CoreElt(core, seq_soln);
   braid_Int     storage       = _braid_CoreElt(core, storage);
   braid_Real    tol           = _braid_CoreElt(core, tol);
   braid_Int     rtol          = _braid_CoreElt(core, rtol);
   braid_Int    *nrels         = _braid_CoreElt(core, nrels);
   /*braid_Int    *cfactors     = _braid_CoreElt(core, cfactors);*/
   braid_Int     max_iter      = _braid_CoreElt(core, max_iter);
   braid_Int     nrefine       = _braid_CoreElt(core, nrefine);
   braid_Int     niter         = _braid_CoreElt(core, niter);
   braid_Int     nlevels       = _braid_CoreElt(core, nlevels);
   braid_Int     tnorm         = _braid_CoreElt(core, tnorm); 
   braid_Int     fmg           = _braid_CoreElt(core, fmg); 
   braid_Int     nfmg          = _braid_CoreElt(core, nfmg); 
   braid_Int     nfmg_Vcyc     = _braid_CoreElt(core, nfmg_Vcyc); 
   braid_Int     access_level  = _braid_CoreElt(core, access_level); 
   braid_Int     print_level   = _braid_CoreElt(core, print_level); 
   braid_Int     skip          = _braid_CoreElt(core, skip); 
   braid_Real    globaltime    = _braid_CoreElt(core, globaltime);
   braid_PtFcnResidual fullres = _braid_CoreElt(core, full_rnorm_res);
   _braid_Grid **grids         = _braid_CoreElt(core, grids);

   braid_Real    rnorm;
   braid_Int     level;

   _braid_GetRNorm(core, -1, &rnorm);
   
   if ( myid == 0 )
   {
      _braid_printf("\n");
      _braid_printf("  start time = %e\n", tstart);
      _braid_printf("  stop time  = %e\n", tstop);
      _braid_printf("  time steps = %d\n", gupper);
      _braid_printf("\n");
      _braid_printf("  use seq soln?         = %d\n", seq_soln);
      _braid_printf("  storage               = %d\n", storage);
      _braid_printf("\n");
      _braid_printf("  stopping tolerance    = %e\n", tol);
      _braid_printf("  use relative tol?     = %d\n", rtol);
      _braid_printf("  max iterations        = %d\n", max_iter);
      _braid_printf("  iterations            = %d\n", niter);
      _braid_printf("  residual norm         = %e\n", rnorm);
      if (tnorm == 1)
      {
         _braid_printf("                         --> 1-norm TemporalNorm \n");
      }
      else if (tnorm == 2)
      {
         _braid_printf("                         --> 2-norm TemporalNorm \n");
      }
      else if (tnorm == 3)
      {
         _braid_printf("                         --> Inf-norm TemporalNorm \n");
      }
      if (fullres != NULL)
      {
         _braid_GetFullRNorm(core, -1, &rnorm);
         _braid_printf("  Global res 2-norm     = %e\n", rnorm);
      }
      _braid_printf("\n");

      _braid_printf("  use fmg?              = %d\n", fmg);
      if ( fmg )
      {
         _braid_printf("  V-cycles / FMG level  = %d\n", nfmg_Vcyc);
         if ( nfmg != -1 )
         {
            _braid_printf("  number fmg cycles     = %d\n", nfmg);
         }
         else
         {
            _braid_printf("  fmg-cycles for all iteratons\n");         
         }
      }
      _braid_printf("  access_level          = %d\n", access_level);
      _braid_printf("  print_level           = %d\n\n", print_level);
      _braid_printf("  max number of levels  = %d\n", max_levels);
      _braid_printf("  min coarse            = %d\n", min_coarse);
      _braid_printf("  number of levels      = %d\n", nlevels);
      _braid_printf("  skip down cycle       = %d\n", skip);
      _braid_printf("  number of refinements = %d\n", nrefine);
      _braid_printf("\n");
      _braid_printf("  level   time-pts   cfactor   nrelax\n");
      for (level = 0; level < nlevels-1; level++)
      {
         _braid_printf("  % 5d  % 8d  % 7d   % 6d\n",
                      level, _braid_GridElt(grids[level], gupper), 
                      _braid_GridElt(grids[level], cfactor), nrels[level]);
      }
      /* Print out coarsest level information */
      _braid_printf("  % 5d  % 8d  \n",
                      level, _braid_GridElt(grids[level], gupper) );
      _braid_printf("\n");
      _braid_printf("  wall time = %f\n", globaltime);
      _braid_printf("\n");
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetMaxLevels(braid_Core  core,
                   braid_Int   max_levels)
{
   braid_Int              old_max_levels = _braid_CoreElt(core, max_levels);
   braid_Int             *nrels          = _braid_CoreElt(core, nrels);
   braid_Int             *cfactors       = _braid_CoreElt(core, cfactors);
   _braid_Grid          **grids          = _braid_CoreElt(core, grids);
   braid_Int              level;

   _braid_CoreElt(core, max_levels) = max_levels;

   nrels = _braid_TReAlloc(nrels, braid_Int, max_levels);
   cfactors = _braid_TReAlloc(cfactors, braid_Int, max_levels);
   grids    = _braid_TReAlloc(grids, _braid_Grid *, max_levels);
   for (level = old_max_levels; level < max_levels; level++)
   {
      nrels[level]    = -1;
      cfactors[level] = 0;
      grids[level]    = NULL;
   }
   _braid_CoreElt(core, nrels)    = nrels;
   _braid_CoreElt(core, cfactors) = cfactors;
   _braid_CoreElt(core, grids)    = grids;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetSkip(braid_Core  core,
              braid_Int   skip)
{
   /* Do not set skip=1 if we do sequential integration first */
   if (_braid_CoreElt(core, seq_soln) == 0)
      _braid_CoreElt(core, skip) = skip;
   else if (skip == 1)
      _braid_printf("  Braid: The skip option is not compatible with SeqSoln option\n");

   return _braid_error_flag;
}


/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetMinCoarse(braid_Core  core,
                   braid_Int   min_coarse)
{
   _braid_CoreElt(core, min_coarse) = min_coarse;

   return _braid_error_flag;
}


/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetPrintLevel(braid_Core  core,
                    braid_Int   print_level)
{
   _braid_CoreElt(core, print_level) = print_level;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetFileIOLevel(braid_Core  core,
                     braid_Int   io_level)
{
   _braid_CoreElt(core, io_level) = io_level;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetPrintFile(braid_Core     core,
                   const char    *printfile_name)
{
   braid_Int  myid = _braid_CoreElt(core, myid_world);
   
   if (myid == 0)
   {
      if ((_braid_printfile = fopen(printfile_name, "w")) == NULL)
      {
         printf("Error: can't open output file %s\n", printfile_name);
         exit(1);
      }
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetDefaultPrintFile(braid_Core     core)
{
   const char fname[] = "braid_runtime.out";
   braid_SetPrintFile(core, fname); 
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetAccessLevel(braid_Core  core,
                     braid_Int   access_level)
{
   _braid_CoreElt(core, access_level) = access_level;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SplitCommworld(const MPI_Comm  *comm_world,
                     braid_Int       px,
                     MPI_Comm        *comm_x,
                     MPI_Comm        *comm_t)
{
   braid_Int myid, xcolor, tcolor;

   /* Create communicators for the time and space dimensions */
   /* The communicators are based on colors and keys (= myid) */
   MPI_Comm_rank( *comm_world, &myid );
   xcolor = myid / px;
   tcolor = myid % px;

   MPI_Comm_split( *comm_world, xcolor, myid, comm_x );
   MPI_Comm_split( *comm_world, tcolor, myid, comm_t );

   return _braid_error_flag;
}


/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetAbsTol(braid_Core  core,
                braid_Real  tol)
{
   _braid_CoreElt(core, tol)  = tol;
   _braid_CoreElt(core, rtol) = 0;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetRelTol(braid_Core  core,
                braid_Real  tol)
{
   _braid_CoreElt(core, tol)  = tol;
   _braid_CoreElt(core, rtol) = 1;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetNRelax(braid_Core  core,
                braid_Int   level,
                braid_Int   nrelax)
{
   braid_Int  *nrels = _braid_CoreElt(core, nrels);

   if (level < 0)
   {
      /* Set default value */
      _braid_CoreElt(core, nrdefault) = nrelax;
   }
   else
   {
      /* Set factor on specified level */
      nrels[level] = nrelax;
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetCFactor(braid_Core  core,
                 braid_Int   level,
                 braid_Int   cfactor)
{
   braid_Int  *cfactors = _braid_CoreElt(core, cfactors);

   if (level < 0)
   {
      /* Set default value */
      _braid_CoreElt(core, cfdefault) = cfactor;
   }
   else
   {
      /* Set factor on specified level */
      cfactors[level] = cfactor;
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetMaxIter(braid_Core  core,
                 braid_Int   max_iter)
{
   braid_Real  *rnorms      = _braid_CoreElt(core, rnorms);
   braid_Real  *full_rnorms = _braid_CoreElt(core, full_rnorms);
   braid_Int    next_iter   = _braid_CoreElt(core, niter) + 1;
   braid_Int    i;

   /* If rnorms has never been allocated, make sure all entries are initialized */
   if (rnorms == NULL)
   {
      next_iter = 0;
   }

   _braid_CoreElt(core, max_iter) = max_iter;

   rnorms = _braid_TReAlloc(rnorms, braid_Real, max_iter+1);
   for (i = next_iter; i <= max_iter; i++)
   {
      rnorms[i] = braid_INVALID_RNORM;
   }

   /* Allocate even if not using full rnorms (simplifies intialization) */
   full_rnorms = _braid_TReAlloc(full_rnorms, braid_Real, max_iter+1);
   for (i = next_iter; i <= max_iter; i++)
   {
      full_rnorms[i] = braid_INVALID_RNORM;
   }

   _braid_CoreElt(core, rnorms) = rnorms;
   _braid_CoreElt(core, full_rnorms) = full_rnorms;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetRefine(braid_Core  core,
                braid_Int   refine)
{
   _braid_CoreElt(core, refine) = refine;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetMaxRefinements(braid_Core  core,
                        braid_Int   max_refinements)
{
   _braid_CoreElt(core, max_refinements) = max_refinements;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetTPointsCutoff(braid_Core  core,
                    braid_Int   tpoints_cutoff)
{
   _braid_CoreElt(core, tpoints_cutoff) = tpoints_cutoff;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetFMG(braid_Core  core)
{
   _braid_CoreElt(core, fmg) = 1;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetNFMGVcyc(braid_Core  core,
                  braid_Int   nfmg_Vcyc)
{
   _braid_CoreElt(core, nfmg_Vcyc) = nfmg_Vcyc;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetNFMG(braid_Core  core,
              braid_Int   k)
{
   _braid_CoreElt(core, nfmg) = k;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetStorage(braid_Core  core,
                 braid_Int   storage)
{
   _braid_CoreElt(core, storage) = storage;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetTemporalNorm(braid_Core  core,
                      braid_Int   tnorm)
{
   _braid_CoreElt(core, tnorm) = tnorm;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetResidual(braid_Core          core, 
                  braid_PtFcnResidual residual)
{
   _braid_CoreElt(core, residual) = residual;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetFullRNormRes(braid_Core          core, 
                      braid_PtFcnResidual residual)
{
   _braid_CoreElt(core, full_rnorm_res) = residual;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetTimeGrid(braid_Core          core,
                  braid_PtFcnTimeGrid tgrid
                  )
{
   _braid_CoreElt(core, tgrid) = tgrid;
   
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetSpatialCoarsen(braid_Core          core, 
                        braid_PtFcnSCoarsen scoarsen)
{
   _braid_CoreElt(core, scoarsen) = scoarsen;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetSpatialRefine(braid_Core         core,
                       braid_PtFcnSRefine srefine)
{
   _braid_CoreElt(core, srefine) = srefine;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetShell(braid_Core          core,
               braid_PtFcnSInit    sinit,
               braid_PtFcnSClone   sclone,
               braid_PtFcnSFree    sfree)
{
   _braid_CoreElt(core, sinit) = sinit;
   _braid_CoreElt(core, sclone) = sclone;
   _braid_CoreElt(core, sfree) = sfree;
   _braid_CoreElt(core, useshell) = 1;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_GetNumIter(braid_Core   core,
                 braid_Int   *niter_ptr)
{
   *niter_ptr =  _braid_CoreElt(core, niter);
   return _braid_error_flag;
} 

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_GetRNorms(braid_Core   core,
                braid_Int   *nrequest_ptr,
                braid_Real  *rnorms)
{
   braid_Real  *rnorms_all = _braid_CoreElt(core, rnorms);
   braid_Int    rnorms_len = _braid_CoreElt(core, niter) + 1;
   
   _braid_GetNEntries(rnorms_all, rnorms_len, nrequest_ptr, rnorms);
   return _braid_error_flag;
} 

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_GetNLevels(braid_Core  core,
                 braid_Int  *nlevels_ptr)

{
   *nlevels_ptr = _braid_CoreElt(core, nlevels);
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_GetSpatialAccuracy( braid_StepStatus  status,
                          braid_Real        loose_tol,
                          braid_Real        tight_tol,
                          braid_Real       *tol_ptr )
{
   braid_Int nrequest   = 2;
   braid_Real stol, tol, rnorm, rnorm0, old_fine_tolx;
   braid_Int level;
   braid_Real l_rnorm, l_ltol, l_ttol, l_tol;
   braid_Real *rnorms = (braid_Real *) malloc( 2*sizeof(braid_Real) ); 
   
   braid_StepStatusGetTol(status, &tol);
   braid_StepStatusGetLevel(status, &level);
   braid_StepStatusGetOldFineTolx(status, &old_fine_tolx);

   /* Get the first and then the current residual norms */
   rnorms[0] = -1.0; rnorms[1] = -1.0;
   braid_StepStatusGetRNorms(status, &nrequest, rnorms);
   if((rnorms[0] == -1.0) && (rnorms[1] != -1.0)){
      rnorm0 = rnorms[1];
   }
   else{
      rnorm0 = rnorms[0];
   }
   nrequest = -2;
   braid_StepStatusGetRNorms(status, &nrequest, rnorms);
   if((rnorms[1] == -1.0) && (rnorms[0] != -1.0)){
      rnorm = rnorms[0];
   }
   else{
      rnorm = rnorms[1];
   }
 

   if ( (level > 0) || (nrequest == 0) || (rnorm0 == -1.0) )
   {
      /* Always return the loose tolerance, if
       * (1) On a coarse grid computation
       * (2) There is no residual history yet (this is the first Braid iteration with skip turned on) */
      *tol_ptr = loose_tol;
   }
   else
   {
      /* Else, do a variable tolerance for the fine grid */
      l_rnorm = -log10(rnorm / rnorm0);
      l_tol   = -log10(tol / rnorm0);
      l_ltol  = -log10(loose_tol);
      l_ttol  = -log10(tight_tol);
         
      if ( l_rnorm >= (7.0/8.0)*l_tol )
      {
         /* Close to convergence, return tight_tol */
         *tol_ptr = tight_tol;
      }
      else
      {
         /* Do linear interpolation between loose_tol and tight_tol (but with respect to log10) */
         stol = (l_rnorm / l_tol) * (l_ttol - l_ltol) + l_ltol;
         *tol_ptr = pow(10, -stol);

         /* The fine grid tolerance MUST never decrease */
         if ( ((*tol_ptr) > old_fine_tolx) && (old_fine_tolx > 0) )
         {
            *tol_ptr = old_fine_tolx;
         }
      }
   }
   
   if (level == 0)
   {
      /* Store this fine grid tolerance */
      braid_StepStatusSetOldFineTolx(status, (*tol_ptr));
      
      /* If we've reached the "tight tolerance", then indicate to Braid that we can halt */
      if ( *tol_ptr == tight_tol )
      {
        braid_StepStatusSetTightFineTolx(status, 1); 
      }
      else
      {
        braid_StepStatusSetTightFineTolx(status, 0); 
      }
   }

   free(rnorms);
    /* printf( "lev: %d, accuracy: %1.2e, nreq: %d, rnorm: %1.2e, rnorm0: %1.2e, loose: %1.2e, tight: %1.2e, old: %1.2e, braid_tol: %1.2e \n", level, *tol_ptr, nrequest, rnorm, rnorm0, loose_tol, tight_tol, old_fine_tolx, tol); */
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetSeqSoln(braid_Core  core,
                 braid_Int   seq_soln)
{
   /* Skip needs to be 0 if we do a sequential integration first */
   _braid_CoreElt(core, seq_soln) = seq_soln;
   if (seq_soln == 1)
      _braid_CoreElt(core, skip) = 0;

   return _braid_error_flag;
}


/*----------------------------------------------------------------------------
* Optimization 
*-----------------------------------------------------------------------------*/
braid_Int
braid_SetTStartTimeaverage(braid_Core core, 
                           braid_Real tstart_obj)
{
  _braid_CoreElt(core, optim->tstart_obj) = tstart_obj;
  
  /* Sanity check */
  if ( tstart_obj < _braid_CoreElt(core, tstart) )
  {
    _braid_printf("\n WARNING: tstart_objective < tstart ! Using default tstart now.\n\n");
    _braid_CoreElt(core, optim->tstart_obj) = _braid_CoreElt(core, tstart);
  }

  return _braid_error_flag;
}

braid_Int
braid_GetTStartTimeaverage(braid_Core core,
                           braid_Real *tstart_obj_ptr)
{
  *tstart_obj_ptr = _braid_CoreElt(core, optim->tstart_obj);

  return _braid_error_flag;
}

braid_Int
braid_SetTStopTimeaverage(braid_Core core, 
                          braid_Real tstop_obj)
{
  _braid_CoreElt(core, optim->tstop_obj) = tstop_obj;
  
  /* Sanity check */
  if ( tstop_obj > _braid_CoreElt(core, tstop) )
  {
    _braid_printf("\n WARNING: tstop_objective > tstop ! Using default tstop now.\n\n");
    _braid_CoreElt(core, optim->tstop_obj) = _braid_CoreElt(core, tstop);
  }


  return _braid_error_flag;
}

braid_Int
braid_GetTStopTimeaverage(braid_Core core,
                          braid_Real *tstop_obj_ptr)
{
  *tstop_obj_ptr = _braid_CoreElt(core, optim->tstop_obj);

  return _braid_error_flag;
}