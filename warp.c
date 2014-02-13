/*BHEADER**********************************************************************
 * Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of WARP.  See file COPYRIGHT for details.
 *
 * WARP is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 ***********************************************************************EHEADER*/

/** \file warp.c
 * \brief Source code for user interface routines.  See warp.h for more information.
 *
 */

#include "_warp.h"
#include "util.h"

#ifndef DEBUG
#define DEBUG 0
#endif

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_Init(MPI_Comm              comm_world,
          MPI_Comm              comm,
          warp_Real             tstart,
          warp_Real             tstop,
          warp_Int              ntime,
          warp_App              app,
          warp_PtFcnPhi         phi,
          warp_PtFcnInit        init,
          warp_PtFcnClone       clone,
          warp_PtFcnFree        free,
          warp_PtFcnSum         sum,
          warp_PtFcnDot         dot,
          warp_PtFcnWrite       write,
          warp_PtFcnBufSize     bufsize,
          warp_PtFcnBufPack     bufpack,
          warp_PtFcnBufUnpack   bufunpack,
          warp_Core            *core_ptr)
{
   _warp_Core           *core;
   warp_Int             *nrels;
   warp_Int              level;
   _warp_AccuracyHandle *accuracy;

   /* Warp default values */
   warp_Int              cfdefault = 2;         /* Default coarsening factor */
   warp_Int              nrdefault = 1;         /* Default number of FC sweeps on each level */
   warp_Int              fmg = 0;               /* Default fmg (0 is off) */
   warp_Int              max_iter = 100;        /* Default max_iter */
   warp_Int              max_levels = 30;       /* Default max_levels */
   warp_Int              print_level = 1;       /* Default print level */
   warp_Int              write_level = 0;       /* Default write level */
   warp_Real             tol = 1.0e-09;         /* Default absolute tolerance */
   warp_Real             rtol = 1.0e-09;        /* Default relative tolerance */

   core = _warp_CTAlloc(_warp_Core, 1);

   _warp_CoreElt(core, comm_world) = comm_world;
   _warp_CoreElt(core, comm)       = comm;
   _warp_CoreElt(core, tstart)     = tstart;
   _warp_CoreElt(core, tstop)      = tstop;
   _warp_CoreElt(core, ntime)      = ntime;
   _warp_CoreElt(core, app)        = app;

   _warp_CoreElt(core, phi)        = phi;
   _warp_CoreElt(core, init)       = init;
   _warp_CoreElt(core, clone)      = clone;
   _warp_CoreElt(core, free)       = free;
   _warp_CoreElt(core, sum)        = sum;
   _warp_CoreElt(core, dot)        = dot;
   _warp_CoreElt(core, write)      = write;
   _warp_CoreElt(core, bufsize)    = bufsize;
   _warp_CoreElt(core, bufpack)    = bufpack;
   _warp_CoreElt(core, bufunpack)  = bufunpack;
   _warp_CoreElt(core, coarsen)    = NULL;
   _warp_CoreElt(core, refine)     = NULL;

   _warp_CoreElt(core, write_level) = write_level;
   _warp_CoreElt(core, print_level) = print_level;
   _warp_CoreElt(core, max_levels) = max_levels;
   _warp_CoreElt(core, tol)        = tol;
   _warp_CoreElt(core, rtol)       = rtol;

   nrels = _warp_TAlloc(warp_Int, max_levels);
   for (level = 0; level < max_levels; level++)
   {
      nrels[level] = -1;
   }
   _warp_CoreElt(core, nrels)      = nrels;
   _warp_CoreElt(core, nrdefault)  = nrdefault;

   _warp_CoreElt(core, cfactors)   = _warp_CTAlloc(warp_Int, max_levels);
   _warp_CoreElt(core, cfdefault)  = cfdefault;

   _warp_CoreElt(core, max_iter)   = max_iter;
   _warp_CoreElt(core, niter)      = 0;
   _warp_CoreElt(core, rnorm)      = 0.0;
   _warp_CoreElt(core, fmg)        = fmg;

   /* Accuracy for spatial solves for using implicit schemes
    *  - accuracy[0] refers to accuracy on level 0
    *  - accuracy[1] refers to accuracy on all levels > 0 */
   accuracy                        = _warp_TAlloc(_warp_AccuracyHandle, 2);
   accuracy[0].matchF              = 0;
   accuracy[0].value               = 1.0e-02;
   accuracy[0].old_value           = 1.0e-02;
   accuracy[0].loose               = 1.0e-02;
   accuracy[0].tight               = 1.0e-02;
   accuracy[0].tight_used          = 0;

   accuracy[1].matchF              = 0;
   accuracy[1].value               = 1.0e-02;
   accuracy[1].old_value           = 1.0e-02;
   accuracy[1].loose               = 1.0e-02;
   accuracy[1].tight               = 1.0e-02;
   accuracy[1].tight_used          = 0;
   _warp_CoreElt(core, accuracy)   = accuracy;

   _warp_CoreElt(core, gupper)     = ntime;

   _warp_CoreElt(core, rfactors)   = NULL;

   _warp_CoreElt(core, nlevels)    = 0;
   _warp_CoreElt(core, grids)      = _warp_CTAlloc(_warp_Grid *, max_levels);

   *core_ptr = core;

   return _warp_error_flag;
}


/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_Drive(warp_Core  core)
{
   warp_Real     tstart      = _warp_CoreElt(core, tstart);
   warp_Real     tstop       = _warp_CoreElt(core, tstop);
   warp_Int      ntime       = _warp_CoreElt(core, ntime);
   warp_Real     tol         = _warp_CoreElt(core, tol);
   warp_Int      rtol        = _warp_CoreElt(core, rtol);
   warp_Int      fmg         = _warp_CoreElt(core, fmg);
   warp_Int      max_iter    = _warp_CoreElt(core, max_iter);
   warp_Int      print_level = _warp_CoreElt(core, print_level);

   warp_Int      nlevels, iter;
   warp_Real     rnorm;
   warp_Real     accuracy;
   warp_Int      ilower, iupper;
   warp_Real    *ta;
   warp_Int      level, fmglevel, down, done, i, refined;
   _warp_Grid   *grid;

   MPI_Comm     comm_world = _warp_CoreElt(core, comm_world);
   warp_Int     myid;
   MPI_Comm_rank(comm_world, &myid);

   level = 0;
   rnorm = -1.0;

   /* Create fine grid */
   _warp_GetDistribution(core, &ilower, &iupper);
   _warp_GridInit(core, 0, ilower, iupper, &grid);

   /* Set t values */
   ta = _warp_GridElt(grid, ta);
   for (i = ilower; i <= iupper; i++)
   {
      ta[i-ilower] = tstart + (((warp_Real)i)/ntime)*(tstop-tstart);
   }

   /* Create a grid hierarchy */
   _warp_InitHierarchy(core, grid);
   nlevels = _warp_CoreElt(core, nlevels);

   /* Set initial values at C-points */
   _warp_InitGuess(core, 0);

   /* Set cycling variables */
   fmglevel = 0;
   if (fmg)
   {
      fmglevel = nlevels-1;
   }
   down = 1;
   done = 0;
   if ((nlevels <= 1) || (tol <= 0.0))
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
            _warp_CFRelax(core, level);

            /* F-relax then restrict */
            _warp_FRestrict(core, level, rnorm, iter, &rnorm);
            /* Set initial guess on next coarser level */
            _warp_InitGuess(core, level+1);

            /* Adjust tolerance */
            if ((level == 0) && (iter == 0))
            {
               if (rtol){
                  tol *= rnorm;
               }
            }

            if (level == 0)
            {
               /* Adjust accuracy of spatial solves for level 0 */
               _warp_SetAccuracy(rnorm, _warp_CoreElt(core, accuracy[0].loose), 
                                 _warp_CoreElt(core, accuracy[0].tight),
                                 _warp_CoreElt(core, accuracy[0].value), tol, &accuracy);
               _warp_CoreElt(core, accuracy[0].old_value) = _warp_CoreElt(core, accuracy[0].value);
               _warp_CoreElt(core, accuracy[0].value)     = accuracy;
               _warp_CoreElt(core, accuracy[0].matchF)    = 1;
               
               if( (print_level >= 2) && (myid == 0) )
               {
                  printf("  **** Accuracy changed to %.2e ****\n", accuracy);
               }
            }

            level++;
         }
         else
         {
            /* Coarsest grid - solve on the up-cycle */

            down = 0;
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
               _warp_FInterp(core, level);
               
               level--;
            }
            else
            {
               fmglevel--;
               down = 1;
            }
         }
         else
         {
            /* Finest grid - refine grid if desired, else check convergence */
            _warp_FRefine(core, &refined);

            if (refined)
            {
               nlevels = _warp_CoreElt(core, nlevels);
            }
            else
            {
               /* Note that this residual is based on an earlier iterate */
               if( (print_level >= 1) && (myid == 0) )
               {
                  printf("  Warp || r_%d || = %e\n", iter, rnorm);
               }

               if (((rnorm < tol) && (_warp_CoreElt(core, accuracy[0].tight_used) == 1)) || (iter == max_iter-1))
               {
                  done = 1;
               }
            }

            iter++;

            if (fmg)
            {
               fmglevel = nlevels-1;
            }
            down = 1;
         }
      }

   }

   /* F-relax and write solution to file */
   _warp_FWrite(core, rnorm, iter, 0, 1);

   _warp_CoreElt(core, niter) = iter;
   _warp_CoreElt(core, rnorm) = rnorm;

   /* Print statistics for this run */
   if( (print_level >= 1) && (myid == 0) )
   {
      warp_PrintStats(core);
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_Destroy(warp_Core  core)
{
   if (core)
   {
      warp_Int      nlevels    = _warp_CoreElt(core, nlevels);
      _warp_Grid  **grids      = _warp_CoreElt(core, grids);
      warp_Int      level;

      _warp_TFree(_warp_CoreElt(core, nrels));
      _warp_TFree(_warp_CoreElt(core, cfactors));
      _warp_TFree(_warp_CoreElt(core, accuracy));
      _warp_TFree(_warp_CoreElt(core, rfactors));
      for (level = 0; level < nlevels; level++)
      {
         _warp_GridDestroy(core, grids[level]);
      }
      _warp_TFree(grids);

      _warp_TFree(core);
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_PrintStats(warp_Core  core)
{
   MPI_Comm     comm_world = _warp_CoreElt(core, comm_world);
   warp_Real    tstart     = _warp_CoreElt(core, tstart);
   warp_Real    tstop      = _warp_CoreElt(core, tstop);
   warp_Int     ntime      = _warp_CoreElt(core, ntime);
   warp_Int     max_levels = _warp_CoreElt(core, max_levels);
   warp_Int    *nrels      = _warp_CoreElt(core, nrels);
   warp_Int     nrdefault  = _warp_CoreElt(core, nrdefault);
   warp_Real    tol        = _warp_CoreElt(core, tol);
   warp_Int     rtol       = _warp_CoreElt(core, rtol);
   /*warp_Int    *cfactors   = _warp_CoreElt(core, cfactors);*/
   warp_Int     cfdefault  = _warp_CoreElt(core, cfdefault);
   warp_Int     max_iter   = _warp_CoreElt(core, max_iter);
   warp_Int     niter      = _warp_CoreElt(core, niter);
   warp_Real    rnorm      = _warp_CoreElt(core, rnorm);
   warp_Int     nlevels    = _warp_CoreElt(core, nlevels);

   warp_Int     myid;

   MPI_Comm_rank(comm_world, &myid);
   if ( myid == 0 )
   {
      printf("\n");
      printf("  start time = %e\n", tstart);
      printf("  stop time  = %e\n", tstop);
      printf("  time steps = %d\n", ntime);
      printf("\n");
      printf("  max number of levels = %d\n", max_levels);
      printf("  number of levels     = %d\n", nlevels);
      printf("  coarsening factor    = %d\n", cfdefault);
      printf("  num F-C relaxations  = %d\n", nrdefault);
      printf("  num rels on level 0  = %d\n", nrels[0]);
      printf("  stopping tolerance   = %e\n", tol);
      printf("  relative tolerance?  = %d\n", rtol);
      printf("  max iterations       = %d\n", max_iter);
      printf("  iterations           = %d\n", niter);
      printf("  residual norm        = %e\n", rnorm);
      printf("\n");
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetLoosexTol(warp_Core  core,
                  warp_Int   level,
                  warp_Real  loose_tol)
{
   if (level < 0)
   {
      /* Set the loose tolerance on all levels. 
       * Index 0 corresponds to level 0, index 1 to all levels > 0. */
      _warp_CoreElt(core, accuracy[0].loose)     = loose_tol;
      _warp_CoreElt(core, accuracy[1].loose)     = loose_tol;

      /* Initialize the current and old value with loose_tol. */
      _warp_CoreElt(core, accuracy[0].value)     = loose_tol;
      _warp_CoreElt(core, accuracy[0].old_value) = loose_tol;
      
      _warp_CoreElt(core, accuracy[1].value)     = loose_tol;
      _warp_CoreElt(core, accuracy[1].old_value) = loose_tol;
   }
   else
   {  
      /* Set the loose tolerance on level level and initialize
       * the current and old value for that level with loose_tol. */
      _warp_CoreElt(core, accuracy[level].loose)     = loose_tol;
      _warp_CoreElt(core, accuracy[level].value)     = loose_tol;
      _warp_CoreElt(core, accuracy[level].old_value) = loose_tol;
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetTightxTol(warp_Core  core,
                  warp_Int   level,
                  warp_Real  tight_tol)
{
   if (level < 0)
   {
      /* Set tight tolerance on all levels. */
      _warp_CoreElt(core, accuracy[0].tight)     = tight_tol;
      _warp_CoreElt(core, accuracy[1].tight)     = tight_tol;
   }
   else
   {  
      /* Set tight tolerance on level level. */
      _warp_CoreElt(core, accuracy[level].tight) = tight_tol;
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetMaxLevels(warp_Core  core,
                  warp_Int   max_levels)
{
   _warp_CoreElt(core, max_levels) = max_levels;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetPrintLevel(warp_Core  core,
                   warp_Int   print_level)
{
   _warp_CoreElt(core, print_level) = print_level;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetWriteLevel(warp_Core  core,
                   warp_Int   write_level)
{
   _warp_CoreElt(core, write_level) = write_level;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SplitCommworld(warp_Core core,
                    warp_Int  px,
                    MPI_Comm  comm_x,
                    MPI_Comm  comm_t)
{
   MPI_Comm comm_world =_warp_CoreElt(core, comm_world);
   warp_Int myid, xcolor, tcolor;

   /* Create communicators for the time and space dimensions */
   /* The communicators are based on colors and keys (= myid) */
   MPI_Comm_rank( comm_world, &myid );
   xcolor = myid / px;
   tcolor = myid % px;

   MPI_Comm_split( comm_world, xcolor, myid, &comm_x );
   MPI_Comm_split( comm_world, tcolor, myid, &comm_t );

   return _warp_error_flag;
}


/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetAbsTol(warp_Core  core,
               warp_Real  tol)
{
   _warp_CoreElt(core, tol) = tol;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetRelTol(warp_Core  core,
               warp_Real  tol)
{
   _warp_CoreElt(core, tol)  = tol;
   _warp_CoreElt(core, rtol) = 1;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetNRelax(warp_Core  core,
               warp_Int   level,
               warp_Int   nrelax)
{
   warp_Int  *nrels = _warp_CoreElt(core, nrels);

   if (level < 0)
   {
      /* Set default value */
      _warp_CoreElt(core, nrdefault) = nrelax;
   }
   else
   {
      /* Set factor on specified level */
      nrels[level] = nrelax;
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetCFactor(warp_Core  core,
                warp_Int   level,
                warp_Int   cfactor)
{
   warp_Int  *cfactors = _warp_CoreElt(core, cfactors);

   if (level < 0)
   {
      /* Set default value */
      _warp_CoreElt(core, cfdefault) = cfactor;
   }
   else
   {
      /* Set factor on specified level */
      cfactors[level] = cfactor;
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetMaxIter(warp_Core  core,
                warp_Int   max_iter)
{
   _warp_CoreElt(core, max_iter) = max_iter;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetFMG(warp_Core  core)
{
   _warp_CoreElt(core, fmg) = 1;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetSpatialCoarsen(warp_Core  core, 
                       warp_PtFcnCoarsen coarsen)
{
   _warp_CoreElt(core, coarsen) = coarsen;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetSpatialRefine(warp_Core  core,
                      warp_PtFcnRefine refine)
{
   _warp_CoreElt(core, refine) = refine;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_GetStatusResidual(warp_Status  status,
                       warp_Real   *rnorm_ptr)
{
   *rnorm_ptr = _warp_StatusElt(status, rnorm);
   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_GetStatusIter(warp_Status  status,
                   warp_Int    *iter_ptr)
{
   *iter_ptr = _warp_StatusElt(status, iter);
   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_GetStatusLevel(warp_Status  status,
                    warp_Int    *level_ptr)
{
   *level_ptr = _warp_StatusElt(status, level);
   return _warp_error_flag;
}


/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_GetStatusDone(warp_Status  status,
                   warp_Int    *done_ptr)
{
   *done_ptr = _warp_StatusElt(status, done);
   return _warp_error_flag;
}


/*--------------------------------------------------------------------------
 * Run some simple tests on the init, write and free routines 
 *--------------------------------------------------------------------------*/
warp_Int
warp_TestInitWrite( warp_App              app, 
                    MPI_Comm              comm_x,
                    warp_Real             t,
                    warp_PtFcnInit        init, 
                    warp_PtFcnWrite       write,
                    warp_PtFcnFree        free)
{
   
   warp_Vector    u ;
   warp_Status    status;
   warp_Int       myid_x;
   char           header[255];
   char           message[255];
   
   _warp_InitStatus( 0.0, 0, 0, 0, &status);
   MPI_Comm_rank( comm_x, &myid_x );
   sprintf(header,  "   warp_TestInitWrite:   ");

   /* Print intro */
   sprintf(message, "\nStarting warp_TestInitWrite\n\n");
   _warp_ParFprintfFlush(stdout, "", message, myid_x);

   /* Test */
   sprintf(message, "Starting Test 1\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "u = init(t=%1.2e)\n", t);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   init(app, t, &u);

   sprintf(message, "write(u) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   write(app, t, status, u);

   sprintf(message, "check output: wrote u for initial condition at t=%1.2e. \n\n",t);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   /* Free variables */
   sprintf(message, "free(u) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, u);
   
   sprintf(message, "Finished\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   return 0;
}

warp_Int
warp_TestClone( warp_App              app,  
                MPI_Comm              comm_x,
                warp_Real             t,
                warp_PtFcnInit        init, 
                warp_PtFcnWrite       write,
                warp_PtFcnFree        free,
                warp_PtFcnClone       clone)
{
   
   warp_Vector  u, v;
   warp_Status  status;
   warp_Int     myid_x;
   char         header[255];
   char         message[255];
   
   _warp_InitStatus( 0.0, 0, 0, 0, &status);
   MPI_Comm_rank( comm_x, &myid_x );
   sprintf(header,  "   warp_TestClone:   ");

   /* Print intro */
   sprintf(message, "\nStarting warp_TestClone\n\n");
   _warp_ParFprintfFlush(stdout, "", message, myid_x);

   /* Test 1 */
   sprintf(message, "Starting Test 1\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "u = init(t=%1.2e)\n", t);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   init(app, t, &u);

   sprintf(message, "v = clone(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, u, &v);
   
   sprintf(message, "write(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   write(app, t, status, u);

   sprintf(message, "write(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   write(app, t, status, v);
   
   sprintf(message, "check output:  wrote u and v for initial condition at t=%1.2e.\n\n", t);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   
   /* Free variables */
   sprintf(message, "free(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, u);

   sprintf(message, "free(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, v);

   sprintf(message, "Finished\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   
   return 0;
}



warp_Int
warp_TestSum( warp_App              app, 
              MPI_Comm              comm_x,
              warp_Real             t,
              warp_PtFcnInit        init, 
              warp_PtFcnWrite       write,
              warp_PtFcnFree        free, 
              warp_PtFcnClone       clone,
              warp_PtFcnSum         sum )  
{
   
   warp_Vector  u, v;
   warp_Status  status;
   warp_Int     myid_x;
   char         header[255];
   char         message[255];
   
   _warp_InitStatus( 0.0, 0, 0, 0, &status);
   MPI_Comm_rank( comm_x, &myid_x );
   sprintf(header,  "   warp_TestSum:   ");

   /* Print intro */
   sprintf(message, "\nStarting warp_TestSum\n\n");
   _warp_ParFprintfFlush(stdout, "", message, myid_x);
   
   /* Test 1 */
   sprintf(message, "Starting Test 1\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "u = init(t=%1.2e)\n", t);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   init(app, t, &u);

   sprintf(message, "v = clone(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, u, &v);

   sprintf(message, "v = u - v\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, u, -1.0, v); 

   sprintf(message, "write(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   write(app, t, status, v);
   
   sprintf(message, "check output:  v should equal the zero vector\n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   
   /* Test 2 */
   sprintf(message, "Starting Test 2\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "v = 2*u + v\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 2.0, u, 1.0, v); 

   sprintf(message, "write(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   write(app, t, status, v);
   
   sprintf(message, "write(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   write(app, t, status, u);
   
   sprintf(message, "check output:  v should equal 2*u \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   /* Free variables */
   sprintf(message, "free(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, u);

   sprintf(message, "free(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, v);

   sprintf(message, "Finished\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   return 0;
}

warp_Int
warp_TestDot( warp_App              app, 
              MPI_Comm              comm_x,
              warp_Real             t,
              warp_PtFcnInit        init, 
              warp_PtFcnFree        free, 
              warp_PtFcnClone       clone,
              warp_PtFcnSum         sum,  
              warp_PtFcnDot         dot,
              warp_Int             *correct)  
{   
   warp_Vector  u, v, w;
   warp_Real    result1, result2, result3;
   warp_Int     myid_x;
   char         header[255];
   char         message[255];
   double       wiggle = 1e-12;
   
   MPI_Comm_rank( comm_x, &myid_x );
   sprintf(header,  "   warp_TestDot:   ");

   /* Initialize the correct flag */
   *correct = 1;

   /* Print intro */
   sprintf(message, "\nStarting warp_TestDot\n\n");
   _warp_ParFprintfFlush(stdout, "", message, myid_x);
   
   /* Test 1 */
   sprintf(message, "Starting Test 1\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "u = init(t=%1.2e)\n", t);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   init(app, t, &u);

   sprintf(message, "dot(u,u) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, u, u, &result1);
   if( fabs(result1) == 0.0)
   {
      sprintf(message, "Warning:  dot(u,u) = 0.0\n"); 
      _warp_ParFprintfFlush(stdout, header, message, myid_x);
   }
   else if( isnan(result1) )
   {
      sprintf(message, "Warning:  dot(u,u) = nan\n"); 
      _warp_ParFprintfFlush(stdout, header, message, myid_x);
      *correct = 0;
   }

   sprintf(message, "v = clone(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, u, &v);

   sprintf(message, "v = u - v \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, u, -1.0, v); 

   sprintf(message, "dot(v,v) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, v, v, &result1);
   if( (fabs(result1) > wiggle) || isnan(result1) )
   {
      *correct = 0;
      _warp_ParFprintfFlush(stdout, header, "Test 1 Failed\n", myid_x);
   }
   else
   {
      _warp_ParFprintfFlush(stdout, header, "Test 1 Passed\n", myid_x);
   }
   sprintf(message, "actual output:    dot(v,v) = %1.2e  \n", result1);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "expected output:  dot(v,v) = 0.0 \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   

   /* Test 2 */
   sprintf(message, "Starting Test 2\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "w = clone(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, u, &w);

   sprintf(message, "w = u + w \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, u, 1.0, w); 

   sprintf(message, "dot(u,u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, u, u, &result1);
   
   sprintf(message, "dot(w,w)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, w, w, &result2);
   if( (fabs(result2/result1 - 4.0) > wiggle) || isnan(result2/result1) )
   {
      *correct = 0;
      _warp_ParFprintfFlush(stdout, header, "Test 2 Failed\n", myid_x);
   }
   else
   {
      _warp_ParFprintfFlush(stdout, header, "Test 2 Passed\n", myid_x);
   }
   sprintf(message, "actual output:    dot(w,w) / dot(u,u) = %1.2e / %1.2e = %1.2e \n",
         result2, result1, result2/result1);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "expected output:  dot(w,w) / dot(u,u) = 4.0 \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   /* Test 3 */
   sprintf(message, "Starting Test 3\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "free(w)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, w);

   sprintf(message, "w = clone(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, u, &w);

   sprintf(message, "w = 0.0*u + 0.5*w \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 0.0, u, 0.5, w); 

   sprintf(message, "dot(u,u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, u, u, &result1);
   
   sprintf(message, "dot(w,w)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, w, w, &result2);
   /* Check Result */
   if( (fabs(result2/result1 - 0.25) > wiggle) || isnan(result2/result1) )

   {
      *correct = 0;
      _warp_ParFprintfFlush(stdout, header, "Test 3 Failed\n", myid_x);
   }
   else
   {
      _warp_ParFprintfFlush(stdout, header, "Test 3 Passed\n", myid_x);
   }
   sprintf(message, "actual output:    dot(w,w) / dot(u,u) = %1.2e / %1.2e = %1.2e \n",
         result2, result1, result2/result1);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "expected output:  dot(w,w) / dot(u,u) = 0.25 \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   /* Test 4 */
   sprintf(message, "Starting Test 4\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   sprintf(message, "free(w)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, w);
   
   sprintf(message, "w = clone(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, u, &w);
   
   sprintf(message, "w = u + 0.5*w \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, u, 0.5, w); 

   sprintf(message, "dot(u,u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, u, u, &result1);
   
   sprintf(message, "dot(w,u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, w, u, &result2);

   sprintf(message, "actual output:    dot(w,u) + dot(u,u) = %1.2e + %1.2e = %1.2e\n", result2, result1, result2+result1);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   
   sprintf(message, "v = u + w \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, u, 1.0, w);   
   
   sprintf(message, "dot(v,u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, w, u, &result3);

   /* Check Result */
   if( (fabs(result2 + result1 - result3)/fabs(result3) > wiggle) || 
       isnan(result2) || isnan(result1) || isnan(result3) )
   {
      *correct = 0;
      _warp_ParFprintfFlush(stdout, header, "Test 4 Failed\n", myid_x);
   }
   else
   {
      _warp_ParFprintfFlush(stdout, header, "Test 4 Passed\n", myid_x);
   }
   sprintf(message, "actual output:    dot(v,u) = %1.2e  \n", result3);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "expected output:  dot(v,u) equals dot(w,u) + dot(u,u) \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
 

   /* Free variables */
   sprintf(message, "free(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, u);

   sprintf(message, "free(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, v);
   
   sprintf(message, "free(w)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, w);

   if(*correct == 1) 
      _warp_ParFprintfFlush(stdout, header, "Finished, all tests passed successfully\n", myid_x);
   else      
      _warp_ParFprintfFlush(stdout, header, "Finished, some tests failed\n", myid_x);

   return 0;
}

warp_Int
warp_TestBuf( warp_App              app,
              MPI_Comm              comm_x,
              warp_Real             t,
              warp_PtFcnInit        init,
              warp_PtFcnFree        free,
              warp_PtFcnSum         sum,  
              warp_PtFcnDot         dot,
              warp_PtFcnBufSize     bufsize,
              warp_PtFcnBufPack     bufpack,
              warp_PtFcnBufUnpack   bufunpack,
              warp_Int             *correct)
{   
   warp_Vector  u, v;
   warp_Real    result1;
   warp_Int     myid_x, size;
   void        *buffer;
   char         header[255];
   char         message[255];
   double       wiggle = 1e-12;
   
   MPI_Comm_rank( comm_x, &myid_x );
   sprintf(header,  "   warp_TestBuf:   ");

   /* Initialize the correct flag */
   *correct = 1;

   /* Print intro */
   sprintf(message, "\nStarting warp_TestBuf\n\n");
   _warp_ParFprintfFlush(stdout, "", message, myid_x);
   
   /* Test 1 */
   sprintf(message, "Starting Test 1\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "u = init(t=%1.2e)\n", t);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   init(app, t, &u);
   
   sprintf(message, "dot(u,u) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, u, u, &result1);
   if( fabs(result1) == 0.0)
   {
      sprintf(message, "Warning:  dot(u,u) = 0.0\n"); 
      _warp_ParFprintfFlush(stdout, header, message, myid_x);
   }
   else if( isnan(result1) )
   {
      sprintf(message, "Warning:  dot(u,u) = nan\n"); 
      _warp_ParFprintfFlush(stdout, header, message, myid_x);
      *correct = 0;
   }

   sprintf(message, "size = bufsize()\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   bufsize(app, &size);
   
   sprintf(message, "buffer = malloc(size)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   buffer = malloc(size);

   sprintf(message, "buffer = bufpack(u, buffer))\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   bufpack(app, u, buffer);

   sprintf(message, "v = bufunpack(buffer)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   bufunpack(app, buffer, &v);

   sprintf(message, "v = u - v \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, u, -1.0, v); 

   sprintf(message, "dot(v,v) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, v, v, &result1);
   if( (fabs(result1) > wiggle) || isnan(result1) )

   {
      *correct = 0;
      _warp_ParFprintfFlush(stdout, header, "Test 1 Failed\n", myid_x);
   }
   else
   {
      _warp_ParFprintfFlush(stdout, header, "Test 1 Passed\n", myid_x);
   }
   sprintf(message, "actual output:    dot(v,v) = %1.2e  \n", result1);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "expected output:  dot(v,v) = 0.0 \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   
   /* Free variables */
   sprintf(message, "free(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, u);

   sprintf(message, "free(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, v);
   
   if(*correct == 1) 
      _warp_ParFprintfFlush(stdout, header, "Finished, all tests passed successfully\n", myid_x);
   else      
      _warp_ParFprintfFlush(stdout, header, "Finished, some tests failed\n", myid_x);

   return 0;
}

warp_Int
warp_TestCoarsenRefine( warp_App          app,
                        MPI_Comm          comm_x,
                        warp_Real         t,
                        warp_Real         f_tminus,
                        warp_Real         f_tplus,
                        warp_Real         c_tminus,
                        warp_Real         c_tplus,
                        warp_PtFcnInit    init,
                        warp_PtFcnWrite   write,
                        warp_PtFcnFree    free,
                        warp_PtFcnClone   clone,
                        warp_PtFcnSum     sum,
                        warp_PtFcnDot     dot,
                        warp_PtFcnCoarsen coarsen,
                        warp_PtFcnRefine  refine,
                        warp_Int         *correct)
 {   
   warp_Vector  u, v, w, uc, vc, wc;
   warp_Real    result1;
   warp_Int     myid_x, level;
   char         header[255];
   char         message[255];
   warp_Status  status;
   
   MPI_Comm_rank( comm_x, &myid_x );
   sprintf(header,  "   warp_TestCoarsenRefine:   ");

   /* Initialize the correct flag */
   *correct = 1;

   /* Print intro */
   sprintf(message, "\nStarting warp_TestCoarsenRefine\n\n");
   _warp_ParFprintfFlush(stdout, "", message, myid_x);
   
   /* Test 1 */
   sprintf(message, "Starting Test 1\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "u = init(t=%1.2e)\n", t);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   init(app, t, &u);

   sprintf(message, "dot(u,u) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, u, u, &result1);
   if( fabs(result1) == 0.0)
   {
      sprintf(message, "Warning:  dot(u,u) = 0.0\n"); 
      _warp_ParFprintfFlush(stdout, header, message, myid_x);
   }
   else if( isnan(result1) )
   {
      sprintf(message, "Warning:  dot(u,u) = nan\n"); 
      _warp_ParFprintfFlush(stdout, header, message, myid_x);
      *correct = 0;
   }

   sprintf(message, "uc = coarsen(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   coarsen(app, t, f_tminus, f_tplus, c_tminus, c_tplus, u, &uc); 

   sprintf(message, "write(uc) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   level = 1;
   _warp_InitStatus( 0.0, 0, level, 0, &status);
   write(app, t, status, uc);

   sprintf(message, "write(u) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   level = 0;
   _warp_InitStatus( 0.0, 0, level, 0, &status);
   write(app, t, status, u);

   sprintf(message, "actual output:   wrote u and spatially coarsened u \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   /* Test 2 */
   sprintf(message, "Starting Test 2\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   
   sprintf(message, "v = clone(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, u, &v);

   sprintf(message, "vc = coarsen(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   coarsen(app, t, f_tminus, f_tplus, c_tminus, c_tplus, v, &vc); 
   
   sprintf(message, "wc = clone(vc)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, vc, &wc);

   sprintf(message, "wc = uc - wc \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, uc, -1.0, wc); 

   sprintf(message, "dot(wc,wc)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, wc, wc, &result1);
   
   /* We expect exact equality between uc and vc */
   if( (fabs(result1) != 0.0) || isnan(result1) )
   {
      *correct = 0;
      _warp_ParFprintfFlush(stdout, header, "Test 2 Failed\n", myid_x);
   }
   else
   {
      _warp_ParFprintfFlush(stdout, header, "Test 2 Passed\n", myid_x);
   }
   sprintf(message, "actual output:    dot(wc,wc) = %1.2e \n", result1);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "expected output:  dot(wc,wc) = 0.0 \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   /* Test 3 */
   sprintf(message, "Starting Test 3\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   
   sprintf(message, "w = clone(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, u, &w);

   sprintf(message, "free(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, u);

   sprintf(message, "free(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, v);

   sprintf(message, "v = refine(vc)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   refine(app, t, f_tminus, f_tplus, c_tminus, c_tplus, vc, &v); 

   sprintf(message, "u = refine(uc)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   refine(app, t, f_tminus, f_tplus, c_tminus, c_tplus, uc, &u); 

   sprintf(message, "v = u - v \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, u, -1.0, v); 

   sprintf(message, "dot(v,v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, v, v, &result1);
   
   /* We expect exact equality between u and v */
   if( (fabs(result1) != 0.0) || isnan(result1) )
   {
      *correct = 0;
      _warp_ParFprintfFlush(stdout, header, "Test 3 Failed\n", myid_x);
   }
   else
   {
      _warp_ParFprintfFlush(stdout, header, "Test 3 Passed\n", myid_x);
   }
   sprintf(message, "actual output:    dot(v,v) = %1.2e \n", result1);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "expected output:  dot(v,v) = 0.0 \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   /* Test 4 */
   sprintf(message, "Starting Test 4\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   
   sprintf(message, "w = w - u \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, u, -1.0, w); 

   sprintf(message, "dot(w,w)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, w, w, &result1);
   
   sprintf(message, "actual output:    dot(w,w) = %1.2e \n", result1);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "%s%s%s", "expected output:  For simple interpolation formulas\n",
                   "                             (e.g., bilinear) and a known function\n",
                   "                             (e.g., constant), dot(w,w) should = 0\n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   /* Free variables */
   sprintf(message, "free(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, u);

   sprintf(message, "free(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, v);

   sprintf(message, "free(w)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, w);
   
   sprintf(message, "free(uc)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, uc);

   sprintf(message, "free(vc)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, vc);

   sprintf(message, "free(wc)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, wc);


   if(*correct == 1) 
      _warp_ParFprintfFlush(stdout, header, "Finished, all tests passed successfully\n", myid_x);
   else      
      _warp_ParFprintfFlush(stdout, header, "Finished, some tests failed\n", myid_x);
   
   return 0;
}


