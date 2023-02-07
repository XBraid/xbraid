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

/**
 * Example:       ex-01-expanded.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-01-expanded
 *
 * Help with:     ex-01-expanded -help
 *
 * Sample run:    mpirun -np 2 ex-01-expanded
 *
 * Description:   solve the scalar ODE 
 *                   u' = lambda u, 
 *                   with lambda=-1 and y(0) = 1
 *
 *                Same as ex-01, only implements more advanced XBraid features.
 *                
 *                When run with the default 10 time steps, the solution is:
 *                $ ./ex-01-expanded
 *                $ cat ex-01-expanded.out.00*
 *                  1.00000000000000e+00
 *                  6.66666666666667e-01
 *                  4.44444444444444e-01
 *                  2.96296296296296e-01
 *                  1.97530864197531e-01
 *                  1.31687242798354e-01
 *                  8.77914951989026e-02
 *                  5.85276634659351e-02
 *                  3.90184423106234e-02
 *                  2.60122948737489e-02
 *                  1.73415299158326e-02
 **/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.h"

/*--------------------------------------------------------------------------
 * User-defined routines and structures
 *--------------------------------------------------------------------------*/

/* App structure can contain anything, and be named anything as well */
typedef struct _braid_App_struct
{
   MPI_Comm  comm;
   double    tstart;
   double    tstop;
   double   *dt;
   int       mydt;
   int       ntime;
   int       rank;
   int       num_syncs;

} my_App;

/* Vector structure can contain anything, and be name anything as well */
typedef struct _braid_Vector_struct
{
   double value;

} my_Vector;

/* helper function for my_timegrid, that generates the time-step sizes. Called
 * before braid_init.  Alternatively, this function could read time step sizes
 * from a file */
int
init_TimeSteps(braid_App  app)
{
   
   int ntime = app->ntime;
   double dt = (app->tstop - app->tstart) / app->ntime;
   int i;
   (app->dt) = (double*) malloc(app->ntime*sizeof(double));
   
   /* example of varying time step size */
   if (app->mydt == 2)
   {
      for (i = 0; i < ntime*0.5; i++)
      {
         app->dt[i] = dt * 0.5;
      }
      for (i = ntime*0.5; i < ntime; i++)
      {
         app->dt[i] = dt * 1.5;
      }
   }
   /* default to constant time step size */
   else
   {
      for (i = 0; i < ntime; i++)
      {
         app->dt[i] = dt;
      }
   }
   
   return 0;
}

int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   int    istart;             /* time point index value corresponding to tstart on (global) fine grid */
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   braid_StepStatusGetTIndex(status, &istart);

   /* Account for XBraid right-hand-side */
   if (fstop != NULL)
   {
      (u->value) += (fstop->value);
   }

   /* Use backward Euler to propagate solution */
   (u->value) = 1./(1. + tstop-tstart)*(u->value);

   /* no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}

int
my_Residual(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     r,
            braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);

   /* compute A(u_i, u_{i-1}) */
   (r->value) = (1. + tstop-tstart)*(ustop->value) - (r->value);

   return 0;
}

int
print_my_timegrid(braid_App        app,
                  braid_Real      *ta,
                  braid_Int       *ilower,
                  braid_Int       *iupper)
{
   int   i, lower, upper;
   char  filename[255];
   FILE *file;
   
   lower = *ilower;
   upper = *iupper;

   /* filename could be anything that helps you track the current time grid */
   sprintf(filename, "timegrid.%d.info", app->rank);
   file = fopen(filename, "w");
   if (file != NULL) {
      for (i = lower; i <= upper; i++)
      {
         fprintf(file, "%d %f\n", i, ta[i-lower]);
      }
   }
   else
   {
      return 1;
   }
   fclose(file);

   return 0;
}

int
my_timegrid(braid_App        app,
            braid_Real      *ta,
            braid_Int       *ilower,
            braid_Int       *iupper)
{
   double tstart;             /* time corresponding to ilower, i.e. lower time index value for this processor */
   int i, lower, upper;
   
   lower = *ilower;
   upper = *iupper;
   
   /* Start from the global tstart to compute the local tstart */
   tstart = app->tstart;
   for (i = 0; i < lower; i++)
   {
      tstart += app->dt[i];
   }
   /* Assign time point values for local time point index values lower:upper */
   for (i = lower; i <= upper; i++)
   {
      ta[i-lower]  = tstart;
      tstart       += app->dt[i];
   }
   print_my_timegrid(app, ta, &lower, &upper);

   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   if (t == app->tstart) /* Initial condition */
   {
      (u->value) = 1.0;
   }
   else /* All other time points set to arbitrary value */
   {
      (u->value) = 0.456;
   }
   *u_ptr = u;

   return 0;
}

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   my_Vector *v;

   v = (my_Vector *) malloc(sizeof(my_Vector));
   (v->value) = (u->value);
   *v_ptr = v;

   return 0;
}

int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u);

   return 0;
}

int
my_Sum(braid_App     app,
       double        alpha,
       braid_Vector  x,
       double        beta,
       braid_Vector  y)
{
   (y->value) = alpha*(x->value) + beta*(y->value);

   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   double dot;

   dot = (u->value)*(u->value);
   *norm_ptr = sqrt(dot);

   return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int        index;
   char       filename[255];
   FILE      *file;
   
   braid_AccessStatusGetTIndex(astatus, &index);
   sprintf(filename, "%s.%04d.%03d", "ex-01-expanded.out", index, app->rank);
   file = fopen(filename, "w");
   fprintf(file, "%.14e\n", (u->value));
   fflush(file);
   fclose(file);

   return 0;
}

int
my_BufSize(braid_App          app,
           int                *size_ptr,
           braid_BufferStatus bstatus)
{
   *size_ptr = sizeof(double);
   return 0;
}

int
my_BufPack(braid_App          app,
           braid_Vector       u,
           void               *buffer,
           braid_BufferStatus bstatus)
{
   double *dbuffer = buffer;

   dbuffer[0] = (u->value);
   braid_BufferStatusSetSize( bstatus, sizeof(double) );

   return 0;
}

int
my_BufUnpack(braid_App          app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus bstatus)
{
   double    *dbuffer = buffer;
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   (u->value) = dbuffer[0];
   *u_ptr = u;

   return 0;
}

int
my_BufAlloc(braid_App          app,
            void               **buffer,
            braid_Int          nbytes,
            braid_BufferStatus bstatus)
{
   *buffer = malloc(nbytes);
   return 0;
}

int
my_BufFree(braid_App          app,
           void               **buffer)
{
   free((char *) *buffer);
   *buffer = NULL;
   return 0;
}

int my_Sync(braid_App        app,
            braid_SyncStatus status)
{
   braid_Int calling_fcn;
   braid_SyncStatusGetCallingFunction(status, &calling_fcn);
   if(calling_fcn == braid_ASCaller_Drive_TopCycle)
      app->num_syncs += 1;

   return 0;
}
/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
   braid_Core    core;
   MPI_Comm      comm;
   my_App       *app;
   double        tstart, tstop;
   int           ntime;

   int           max_levels    = 2;
   int           nrelax        = 1;
   int           nrelax0       = -1;
   int           nrelaxc       = 5;
   double        tol           = 1.0e-06;
   int           cfactor       = 2;
   int           max_iter      = 100;
   int           timings       = 0;
   int           fmg           = 0;
   int           skip          = 1;
   int           res           = 0;
   int           mydt          = 0;
   int           sync          = 0;
   int           periodic      = 0;
   int           relax_only_cg = 0;
   int           bufalloc      = 0;

   int           arg_index;
   int           rank;

   /* Initialize MPI */
   comm   = MPI_COMM_WORLD;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(comm, &rank);

   /* Define time domain: ntime intervals */
   ntime  = 10;
   tstart = 0.0;
   tstop  = tstart + ntime/2.;
   
   /* Parse command line */
   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         if ( rank == 0 )
         {
            printf("\nExample 1: Solve a scalar ODE \n\n");
            printf("  -ntime <ntime>    : set num time points\n");
            printf("  -ml  <max_levels> : set max levels\n");
            printf("  -nu  <nrelax>     : set num F-C relaxations\n");
            printf("  -nu0 <nrelax>     : set num F-C relaxations on level 0\n");
            printf("  -nuc <nrelax>     : set num F-C relaxations on coarsest grid\n");
            printf("  -tol <tol>        : set stopping tolerance\n");
            printf("  -cf  <cfactor>    : set coarsening factor\n");
            printf("  -mi  <max_iter>   : set max iterations\n");
            printf("  -skip <set_skip>  : set skip relaxations on first down-cycle; 0: no skip;  1: skip\n");
            printf("  -timings <bool>   : turn XBraid internal timings on/off\n");
            printf("  -fmg              : use FMG cycling\n");
            printf("  -res              : use my residual\n");
            printf("  -sync             : enable calls to the sync function\n");
            printf("  -periodic         : solve a periodic problem\n");
            printf("  -bufalloc         : user-defined MPI buffer allocation\n");
            printf("  -tg <mydt>        : use user-specified time grid as global fine time grid, options are\n");
            printf("                      1 - uniform time grid\n");
            printf("                      2 - nonuniform time grid, where dt*0.5 for n = 1, ..., nt/2; dt*1.5 for n = nt/2+1, ..., nt\n\n");
         }
         exit(1);
      } 
      else if ( strcmp(argv[arg_index], "-ntime") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
         tstop  = tstart + ntime/2.;
      }
      else if ( strcmp(argv[arg_index], "-ml") == 0 )
      {
         arg_index++;
         max_levels = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu") == 0 )
      {
         arg_index++;
         nrelax = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu0") == 0 )
      {
         arg_index++;
         nrelax0 = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nuc") == 0 )
      {
         arg_index++;
         nrelaxc = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tol") == 0 ) 
      {
         arg_index++;
         tol = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-cf") == 0 )
      {
         arg_index++;
         cfactor = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mi") == 0 )
      {
         arg_index++;
         max_iter = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-timings") == 0 )
      {
         arg_index++;
         timings = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-fmg") == 0 )
      {
         arg_index++;
         fmg = 1;
      }
      else if ( strcmp(argv[arg_index], "-res") == 0 )
      {
         arg_index++;
         res = 1;
      }
      else if ( strcmp(argv[arg_index], "-tg") == 0 )
      {
         arg_index++;
         mydt = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-sync") == 0 )
      {
         arg_index++;
         sync = 1;
      }
      else if( strcmp(argv[arg_index], "-periodic") == 0 )
      {
         arg_index++;
         periodic = 1;
      }
      else if( strcmp(argv[arg_index], "-bufalloc") == 0 )
      {
         arg_index++;
         bufalloc = 1;
      }
      else if( strcmp(argv[arg_index], "-relax_only_cg") == 0 )
      {
         arg_index++;
         relax_only_cg = 1;
      }
      else if ( strcmp(argv[arg_index], "-skip") == 0 )
      {
         arg_index++;
         skip = atoi(argv[arg_index++]);
      }
      else
      {
         arg_index++;
      }
   }

   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->dt)     = NULL;
   (app->mydt)   = mydt;
   (app->comm)   = comm;
   (app->tstart) = tstart;
   (app->tstop)  = tstop;
   (app->ntime)  = ntime;
   (app->rank)   = rank;
   (app->num_syncs) = 0;

   /* initialize XBraid and set options */
   braid_Init(comm, comm, tstart, tstop, ntime, app,
             my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
             my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

   braid_SetPrintLevel( core, 2);
   braid_SetMaxLevels(core, max_levels);
   braid_SetNRelax(core, -1, nrelax);
   if (nrelax0 > -1)
   {
      braid_SetNRelax(core,  0, nrelax0);
   }
   braid_SetNRelax(core, max_levels-1, nrelaxc);
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, cfactor);
   braid_SetMaxIter(core, max_iter);
   braid_SetTimings(core, timings);
   braid_SetSkip(core, skip);
   if (fmg)
   {
      braid_SetFMG(core);
   }
   if (res)
   {
      braid_SetResidual(core, my_Residual);
   }
   if (app->mydt > 0)
   {
      init_TimeSteps(app);
      braid_SetTimeGrid(core, my_timegrid);
   }
   if (sync)
   {
      braid_SetSync(core, my_Sync);
   }
   if (periodic)
   {
      braid_SetPeriodic(core, periodic);
   }
   if (relax_only_cg)
   {
      braid_SetRelaxOnlyCG(core, relax_only_cg);
   }
   if (bufalloc)
   {
      braid_SetBufAllocFree(core, my_BufAlloc, my_BufFree);
   }

   /* Run simulation, and then clean up */
   braid_Drive(core);

   if (sync && rank == 0)
      printf("  num_syncs             = %d\n\n", (app->num_syncs));

   braid_Destroy(core);
   free(app);
   MPI_Finalize();

   return (0);
}
