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
#include "util.h"

#ifndef DEBUG
#define DEBUG 0
#endif

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
   braid_Int              cfdefault    = 2;         /* Default coarsening factor */
   braid_Int              nrdefault    = 1;         /* Default number of FC sweeps on each level */
   braid_Int              fmg          = 0;         /* Default fmg (0 is off) */
   braid_Int              nfmg         = -1;        /* Default fmg cycles is -1, indicating all fmg-cycles (if fmg=1) */
   braid_Int              nfmg_Vcyc    = 1;         /* Default num V-cycles at each fmg level is 1 */
   braid_Int              max_iter     = 100;       /* Default max_iter */
   braid_Int              max_levels   = 30;        /* Default max_levels */
   braid_Int              min_coarse   = 2;         /* Default min_coarse */
   braid_Int              print_level  = 1;         /* Default print level */
   braid_Int              access_level = 1;         /* Default access level */
   braid_Int              tnorm        = 2;         /* Default temporal norm */
   braid_Real             tol          = 1.0e-09;   /* Default absolute tolerance */
   braid_Real             rtol         = 1;         /* Use relative tolerance */

   core = _braid_CTAlloc(_braid_Core, 1);

   _braid_CoreElt(core, comm_world)    = comm_world;
   _braid_CoreElt(core, comm)          = comm;
   _braid_CoreElt(core, tstart)        = tstart;
   _braid_CoreElt(core, tstop)         = tstop;
   _braid_CoreElt(core, ntime)         = ntime;
   _braid_CoreElt(core, app)           = app;

   _braid_CoreElt(core, step)           = step;
   _braid_CoreElt(core, init)          = init;
   _braid_CoreElt(core, clone)         = clone;
   _braid_CoreElt(core, free)          = free;
   _braid_CoreElt(core, sum)           = sum;
   _braid_CoreElt(core, spatialnorm)   = spatialnorm;
   _braid_CoreElt(core, access)        = access;
   _braid_CoreElt(core, bufsize)       = bufsize;
   _braid_CoreElt(core, bufpack)       = bufpack;
   _braid_CoreElt(core, bufunpack)     = bufunpack;
   _braid_CoreElt(core, residual)      = NULL;
   _braid_CoreElt(core, globresidual)  = NULL;
   _braid_CoreElt(core, coarsen)       = NULL;
   _braid_CoreElt(core, refine)        = NULL;

   _braid_CoreElt(core, access_level)  = access_level;
   _braid_CoreElt(core, tnorm)         = tnorm;
   _braid_CoreElt(core, print_level)   = print_level;
   _braid_CoreElt(core, max_levels)    = 0; /* Set with SetMaxLevels() below */
   _braid_CoreElt(core, min_coarse)    = min_coarse;
   _braid_CoreElt(core, tol)           = tol;
   _braid_CoreElt(core, rtol)          = rtol;

   _braid_CoreElt(core, nrels)      = NULL; /* Set with SetMaxLevels() below */
   _braid_CoreElt(core, nrdefault)  = nrdefault;

   _braid_CoreElt(core, cfactors)   = NULL; /* Set with SetMaxLevels() below */
   _braid_CoreElt(core, cfdefault)  = cfdefault;

   _braid_CoreElt(core, max_iter)   = 0; /* Set with SetMaxIter() below */
   _braid_CoreElt(core, niter)      = 0;
   _braid_CoreElt(core, fmg)        = fmg;
   _braid_CoreElt(core, nfmg)       = nfmg;
   _braid_CoreElt(core, nfmg_Vcyc)  = nfmg_Vcyc;

   _braid_CoreElt(core, astatus)    = _braid_CTAlloc(_braid_AccessStatus, 1);
   _braid_CoreElt(core, sstatus)    = _braid_CTAlloc(_braid_StepStatus, 1);
   _braid_CoreElt(core, cstatus)    = _braid_CTAlloc(_braid_CoarsenRefStatus, 1);

   _braid_CoreElt(core, storage)    = -1;            /* only store C-points */

   _braid_CoreElt(core, gupper)     = ntime;

   _braid_CoreElt(core, rfactors)   = NULL;
   _braid_CoreElt(core, nrefine)    = 0;

   _braid_CoreElt(core, nlevels)    = 0;
   _braid_CoreElt(core, grids)      = NULL; /* Set with SetMaxLevels() below */

   /* Residual history and accuracy tracking for StepStatus*/
   _braid_CoreElt(core, rnorms)                       = NULL; /* Set with SetMaxIter() below */
   _braid_StatusElt(
      _braid_CoreElt(core, sstatus), rnorms)          = NULL; /* Set with SetMaxIter() below */
   _braid_StatusElt(
      _braid_CoreElt(core, sstatus), rnorms_len_ptr)  = &(_braid_CoreElt(core, rnorms_len));
   _braid_CoreElt(core, rnorms_len)                                  = 0;
   _braid_StatusElt( _braid_CoreElt(core, sstatus), old_fine_tolx)   = -1.0;
   _braid_StatusElt( _braid_CoreElt(core, sstatus), tight_fine_tolx) = 1;

   braid_SetMaxLevels(core, max_levels);
   braid_SetMaxIter(core, max_iter);

   *core_ptr = core;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_Drive(braid_Core  core)
{
   braid_Real     tstart       = _braid_CoreElt(core, tstart);
   braid_Real     tstop        = _braid_CoreElt(core, tstop);
   braid_Int      ntime        = _braid_CoreElt(core, ntime);
   braid_Real     tol          = _braid_CoreElt(core, tol);
   braid_Int      rtol         = _braid_CoreElt(core, rtol);
   braid_Int      fmg          = _braid_CoreElt(core, fmg);
   braid_Int      nfmg         = _braid_CoreElt(core, nfmg);
   braid_Int      max_levels   = _braid_CoreElt(core, max_levels);
   braid_Int      max_iter     = _braid_CoreElt(core, max_iter);
   braid_Int      print_level  = _braid_CoreElt(core, print_level);
   braid_Int      access_level = _braid_CoreElt(core, access_level);
   braid_Int      nfmg_Vcyc    = _braid_CoreElt(core, nfmg_Vcyc); 
   braid_PtFcnResidual globres = _braid_CoreElt(core, globresidual);
   braid_Real*    rnorms       = _braid_CoreElt(core, rnorms);
   braid_StepStatus    sstatus = _braid_CoreElt(core, sstatus);


   braid_Int     *nrels, nrel0;
   braid_Int      nlevels, iter, nprocs, tight_fine_tolx;
   braid_Real     global_rnorm, old_globalrnorm, braid_rnorm;
   braid_Int      ilower, iupper;
   braid_Real    *ta;
   braid_Int      level, fmglevel, fmg_Vcyc, down, done, i, refined;
   _braid_Grid   *grid;
   braid_Real     localtime, globaltime;

   MPI_Comm       comm       = _braid_CoreElt(core, comm);
   MPI_Comm       comm_world = _braid_CoreElt(core, comm_world);
   braid_Int      myid;

   /* Check that nprocs <= npoints */
   MPI_Comm_size(comm, &nprocs);

   /* Start timer */
   localtime = MPI_Wtime();

   MPI_Comm_rank(comm_world, &myid);

   level = 0;
   global_rnorm = -1.0;
   rnorms[0] = -1.0;

   /* Create fine grid */
   _braid_GetDistribution(core, &ilower, &iupper);
   _braid_GridInit(core, 0, ilower, iupper, &grid);

   /* Set t values */
   ta = _braid_GridElt(grid, ta);
   for (i = ilower; i <= iupper; i++)
   {
      ta[i-ilower] = tstart + (((braid_Real)i)/ntime)*(tstop-tstart);
   }

   /* Create a grid hierarchy */
   _braid_InitHierarchy(core, grid, 0);
   nlevels = _braid_CoreElt(core, nlevels);

   /* Set initial values */
   _braid_InitGuess(core, 0);

   /* Set cycling variables */
   fmglevel = 0;
   fmg_Vcyc = 0;
   if (nfmg == 0)
   {
      fmg = 0;
   }
   if (fmg)
   {
      fmglevel = nlevels-1;
   }
   down = 1;
   done = 0;
   if (max_levels <= 1)
   {
      /* Just do sequential time marching */
      done = 1;
   }

   iter = 0;
   while (!done)
   {
      /* Down cycle */
      if (down)
      {
         if (level < (nlevels-1))
         {
            /* CF-relaxation */
            _braid_FCRelax(core, level);

            /* Check to use user provided global residual */
            if( (level == 0) &&  (globres != NULL) )
            {
               old_globalrnorm = global_rnorm;
               _braid_GetFullResidual(core, level, &global_rnorm);
            }
            
            /* F-relax then restrict
             * Note that FRestrict computes a new rnorm */
            _braid_FRestrict(core, level);
            braid_rnorm = rnorms[ _braid_CoreElt(core, rnorms_len)-1 ];
            
            /* Adjust tolerance */
            if ((level == 0) && (iter == 0) && (rtol))
            {
               /* Check to use user provided global residual */
               if(globres != NULL)
               {
                  tol *= global_rnorm; 
               }
               else
               {
                  tol *= braid_rnorm;
               }
               _braid_CoreElt(core, tol) = tol;
            }

            level++;
         }
         else
         {
            if (nlevels == 1)
            {
               /* Do sequential time marching - refine on the up-cycle */
               nrels = _braid_CoreElt(core, nrels);
               nrel0 = nrels[0];
               nrels[0] = 1;
               _braid_FCRelax(core, 0);
               nrels[0] = nrel0;
               down = 0;
            }
            else
            {
               /* Coarsest grid - solve on the up-cycle */
               down = 0;
            }
         }
      }

      /* Up cycle */

      if (!down)
      {
         if (level > 0)
         {
            if (level >= fmglevel)
            {
               /* F-relax then interpolate */
               _braid_FInterp(core, level);
               level--;
            }
            else
            {  
               /* Do nfmg_Vcyc number of V-cycles at each level in FMG */
               fmg_Vcyc += 1;
               if ( fmg_Vcyc == nfmg_Vcyc )
               {
                  fmg_Vcyc = 0;
                  fmglevel--;
               }
               
               down = 1;
            }
         }
         else
         {
            /* Finest grid - refine grid if desired, else check convergence */
            _braid_FRefine(core, &refined);

            if (refined)
            {
               nlevels = _braid_CoreElt(core, nlevels);
            }
            else
            {
               if (nlevels == 1)
               {
                  rnorms[0] = 0.0;
                  _braid_CoreElt(core, rnorms_len) = 1;
               }

               /* Note that this residual is based on an earlier iterate */
               if( (print_level >= 1) && (myid == 0) )
               {
                  if (iter == 0)
                  {
                     if (globres != NULL)
                     {
                        _braid_printf("  Braid: Global || r_%d || = %1.6e\n", iter, global_rnorm);
                     }
                     _braid_printf("  Braid: || r_%d || = %1.6e,  wall time = %1.2e\n", 
                        iter,  braid_rnorm, (MPI_Wtime() - localtime));
                  }
                  else
                  {
                     if (globres != NULL)
                     {
                        _braid_printf("  Braid: Global || r_%d || = %1.6e,  conv. factor = %1.2e \n",
                           iter, global_rnorm, global_rnorm/old_globalrnorm);
                     }
                     _braid_printf("  Braid: || r_%d || = %1.6e,  conv. factor = %1.2e, wall time = %1.2e\n",  
                        iter, braid_rnorm, braid_rnorm/rnorms[ _braid_CoreElt(core, rnorms_len)-2 ], 
                        (MPI_Wtime() - localtime));
                  }
               }
               
               /* Use user provided global residual as stopping criterion? */
               tight_fine_tolx = _braid_StatusElt( sstatus, tight_fine_tolx);
               if (globres != NULL)
               {
                  if ( ((global_rnorm < tol) && (tight_fine_tolx == 1))  ||  (iter == max_iter-1) )
                  {
                     done = 1;
                  }
                  else if ( isnan(global_rnorm) )
                  {
                     if (myid == 0)
                     {
                        _braid_printf("  Iterations diverged.\n");
                     }
                     done = 1; 
                  }
               }
               else 
               {
                  if ( ((braid_rnorm < tol) && (tight_fine_tolx == 1))  ||  (iter == max_iter-1) )
                  {
                     done = 1;
                  } 
                  else if ( isnan(braid_rnorm) )
                  {
                     if (myid == 0)
                     {
                        _braid_printf("  Iterations diverged.\n");
                     }
                     done = 1; 
                  }
               }

               iter++;
              _braid_CoreElt(core, niter) = iter;
               
               /* Check if done doing fmg, if not reset fmglevel. */
               if (fmg)
               {
                  if ( (nfmg-iter) == 0 )
                  {
                     fmg = 0;
                     fmglevel = 0;
                  }
                  else {
                     fmglevel = nlevels-1;                     
                  }
               }
            }

            down = 1;
         }
      }
   }

   /* Get final residual */
   if (globres != NULL)
   {
      _braid_GetFullResidual(core, level, &global_rnorm);
      /* RDF - Why are these next two lines needed? 
       * JBS - Ben S wanted a final rnorm, we should move this final residual computation to FAccess to save work */
      _braid_FRestrict(core, level);
   }

   /* Allow final access to Braid by carrying out an F-relax to generate points */
   if( access_level >= 1)
   {
      _braid_FAccess(core, 0, 1);
   }

   /* Stop timer */
   localtime = MPI_Wtime() - localtime;
   MPI_Allreduce(&localtime, &globaltime, 1, braid_MPI_REAL, MPI_MAX, comm_world);
   _braid_CoreElt(core, localtime)  = localtime;
   _braid_CoreElt(core, globaltime) = globaltime;
   _braid_CoreElt(core, global_rnorm) = global_rnorm;

   /* Print statistics for this run */
   if( (print_level > 0) && (myid == 0) )
   {
      braid_PrintStats(core);
   }

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
      braid_AccessStatus      astatus    = _braid_CoreElt(core, astatus);
      braid_CoarsenRefStatus  cstatus    = _braid_CoreElt(core, cstatus);
      braid_StepStatus        sstatus    = _braid_CoreElt(core, sstatus);
      braid_Int               level;

      _braid_TFree(_braid_CoreElt(core, nrels));
      _braid_TFree(_braid_CoreElt(core, rnorms));
      _braid_TFree(_braid_CoreElt(core, cfactors));
      _braid_TFree(_braid_CoreElt(core, rfactors));
      _braid_TFree(_braid_CoreElt(core, tnorm_a));
      _braid_AccessStatusDestroy(astatus);
      _braid_StepStatusDestroy(sstatus);
      _braid_CoarsenRefStatusDestroy(cstatus);
      
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
   MPI_Comm      comm_world    = _braid_CoreElt(core, comm_world);
   braid_Real    tstart        = _braid_CoreElt(core, tstart);
   braid_Real    tstop         = _braid_CoreElt(core, tstop);
   braid_Int     ntime         = _braid_CoreElt(core, ntime);
   braid_Int     max_levels    = _braid_CoreElt(core, max_levels);
   braid_Int     min_coarse    = _braid_CoreElt(core, min_coarse);
   braid_Real    tol           = _braid_CoreElt(core, tol);
   braid_Int     rtol          = _braid_CoreElt(core, rtol);
   braid_Int    *nrels         = _braid_CoreElt(core, nrels);
   /*braid_Int    *cfactors     = _braid_CoreElt(core, cfactors);*/
   braid_Int     max_iter      = _braid_CoreElt(core, max_iter);
   braid_Int     niter         = _braid_CoreElt(core, niter);
   braid_Real*   rnorms        = _braid_CoreElt(core, rnorms);
   braid_Int     rnorms_len    = _braid_CoreElt(core, rnorms_len);
   braid_Real    global_rnorm  = _braid_CoreElt(core, global_rnorm);
   braid_Int     nlevels       = _braid_CoreElt(core, nlevels);
   braid_Int     tnorm         = _braid_CoreElt(core, tnorm); 
   braid_Int     fmg           = _braid_CoreElt(core, fmg); 
   braid_Int     nfmg          = _braid_CoreElt(core, nfmg); 
   braid_Int     nfmg_Vcyc     = _braid_CoreElt(core, nfmg_Vcyc); 
   braid_Int     access_level  = _braid_CoreElt(core, access_level); 
   braid_Int     print_level   = _braid_CoreElt(core, print_level); 
   braid_Real    globaltime    = _braid_CoreElt(core, globaltime);
   braid_PtFcnResidual globres = _braid_CoreElt(core, globresidual);
   _braid_Grid **grids         = _braid_CoreElt(core, grids);

   braid_Int     myid, level;

   MPI_Comm_rank(comm_world, &myid);
   if ( myid == 0 )
   {
      _braid_printf("\n");
      _braid_printf("  start time = %e\n", tstart);
      _braid_printf("  stop time  = %e\n", tstop);
      _braid_printf("  time steps = %d\n", ntime);
      _braid_printf("\n");
      _braid_printf("  stopping tolerance   = %e\n", tol);
      _braid_printf("  use relative tol?    = %d\n", rtol);
      _braid_printf("  max iterations       = %d\n", max_iter);
      _braid_printf("  iterations           = %d\n", niter);
      _braid_printf("  residual norm        = %e\n", rnorms[rnorms_len-1]);
      if(tnorm == 1){
         _braid_printf("                        --> 1-norm TemporalNorm \n"); 
      }
      else if(tnorm == 2){
         _braid_printf("                        --> 2-norm TemporalNorm \n"); 
      }
      else if(tnorm == 3){
         _braid_printf("                        --> Inf-norm TemporalNorm \n"); 
      }
      if (globres != NULL) {
         _braid_printf("  Global res 2-norm    = %e\n", global_rnorm);
      }
      _braid_printf("\n");

      _braid_printf("  use fmg?             = %d\n", fmg);
      if ( fmg ) {
         _braid_printf("  V-cycles / FMG level = %d\n", nfmg_Vcyc);
         if ( nfmg != -1 ) {
            _braid_printf("  number fmg cycles    = %d\n", nfmg);                     
         }
         else {
            _braid_printf("  fmg-cycles for all iteratons\n");         
         }
      }
      _braid_printf("  access_level         = %d\n", access_level);
      _braid_printf("  print_level          = %d\n\n", print_level);
      _braid_printf("  max number of levels = %d\n", max_levels);
      _braid_printf("  min coarse           = %d\n", min_coarse);
      _braid_printf("  number of levels     = %d\n", nlevels);
      _braid_printf("\n");
      _braid_printf("  level   time-pts   cfactor   nrelax\n", globaltime);
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
braid_SetPrintFile(braid_Core     core,
                   const char    *printfile_name)
{
   MPI_Comm   comm_world = _braid_CoreElt(core, comm_world);
   braid_Int     myid;
   
   MPI_Comm_rank(comm_world, &myid);
   
   if(myid == 0)
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
   braid_Real*  rnorms = _braid_CoreElt(core, rnorms);

   _braid_CoreElt(core, max_iter) = max_iter;

   rnorms = _braid_TReAlloc(rnorms, braid_Real, max_iter);

   _braid_CoreElt(core, rnorms) = rnorms;
   _braid_StatusElt( _braid_CoreElt(core, sstatus), rnorms) = rnorms;

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
braid_SetResidual(braid_Core  core, 
                  braid_PtFcnResidual residual)
{
   _braid_CoreElt(core, residual) = residual;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetGlobalResidual(braid_Core  core, 
                        braid_PtFcnResidual residual)
{
   _braid_CoreElt(core, globresidual) = residual;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetSpatialCoarsen(braid_Core  core, 
                        braid_PtFcnCoarsen coarsen)
{
   _braid_CoreElt(core, coarsen) = coarsen;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_SetSpatialRefine(braid_Core  core,
                    braid_PtFcnRefine refine)
{
   _braid_CoreElt(core, refine) = refine;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_GetNumIter(braid_Core  core,
              braid_Int   *niter_ptr)
{
   *niter_ptr =  _braid_CoreElt(core, niter);
   return _braid_error_flag;
} 

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_GetRNorms(braid_Core  core,
                braid_Int   *nrequest_ptr,
                braid_Real  *rnorms)

{
   braid_Real     *_rnorms   = _braid_CoreElt(core, rnorms);
   braid_Int      rnorms_len = _braid_CoreElt(core, rnorms_len);
   
   _braid_GetNEntries(_rnorms, rnorms_len, nrequest_ptr, rnorms);
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
   braid_Int nrequest   = 1;
   braid_Real stol, tol, rnorm, rnorm0, old_fine_tolx;
   braid_Int level;
   braid_Real l_rnorm, l_ltol, l_ttol, l_tol;
   
   braid_StepStatusGetTol(status, &tol);
   braid_StepStatusGetLevel(status, &level);
   braid_StepStatusGetOldFineTolx(status, &old_fine_tolx);

   /* Get the first and then the current residual norms */
   braid_StepStatusGetRNorms(status, &nrequest, &rnorm0);
   nrequest = -1;
   braid_StepStatusGetRNorms(status, &nrequest, &rnorm);

   if( (level > 0) || (nrequest == 0) )
   {
      /* Always return the loose tolerance, if
       * (1) On a coarse grid computation
       * (2) There is no residual history yet (this is the first Braid iteration) */
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
         if( ((*tol_ptr) > old_fine_tolx) && (old_fine_tolx > 0) )
         {
            *tol_ptr = old_fine_tolx;
         }
      }
   }
   
   if(level == 0)
   {
      /* Store this fine grid tolerance */
      braid_StepStatusSetOldFineTolx(status, (*tol_ptr));
      
      /* If we've reached the "tight tolerance", then indicate to Braid that we can halt */
      if( *tol_ptr == tight_tol )
      {
        braid_StepStatusSetTightFineTolx(status, 1); 
      }
      else
      {
        braid_StepStatusSetTightFineTolx(status, 0); 
      }
   }

    /* printf( "lev: %d, accuracy: %1.2e, nreq: %d, rnorm: %1.2e, rnorm0: %1.2e, loose: %1.2e, tight: %1.2e, old: %1.2e, braid_tol: %1.2e \n", level, *tol_ptr, nrequest, rnorm, rnorm0, loose_tol, tight_tol, old_fine_tolx, tol); */
   return _braid_error_flag;
}
