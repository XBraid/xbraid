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


/*
   Example 01

   Compile with: make ex-01

   Sample run:   mpirun -np 2 ex-01

   Description:

   Blah...

   When run with the default 10 time steps, the solution is as follows:

      1.00000000000000e+00
      5.00000000000000e-01
      2.50000000000000e-01
      1.25000000000000e-01
      6.25000000000000e-02
      3.12500000000000e-02
      1.56250000000000e-02
      7.81250000000000e-03
      3.90625000000000e-03
      1.95312500000000e-03
      9.76562500000000e-04

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.h"

/*--------------------------------------------------------------------------
 * My integration routines
 *--------------------------------------------------------------------------*/

/* can put anything in my app and name it anything as well */
typedef struct _braid_App_struct
{
   MPI_Comm  comm, comm_t, comm_x;
   double    tstart;
   double    tstop;
   double   *dt;
   int       mydt;
   int       ntime;
   int       ptx, pt, px;
   int       rank, rank_t, rank_x;

} my_App;

/* can put anything in my vector and name it anything as well */
typedef struct _braid_Vector_struct
{
   double value;

} my_Vector;

int
init_TimeSteps(braid_App  app)
{
   
   int ntime = app->ntime;
   double dt = (app->tstop - app->tstart) / app->ntime;
   int i;
   (app->dt) = (double*) malloc(app->ntime*sizeof(double));
   
   /* example on varying time step size */
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
   int    istop;              /* time point index value corresponding to tstop on (global) fine grid */
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   braid_StepStatusGetIstop(status, &istop);

   /* On the finest grid, each value is half the previous value */
   (u->value) = pow(0.5, tstop-tstart)*(u->value);

   if (fstop != NULL)
   {
      /* Nonzero rhs */
      (u->value) += (fstop->value);
   }

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

   /* On the finest grid, each value is half the previous value */
   (r->value) = (ustop->value) - pow(0.5, tstop-tstart)*(r->value);

   return 0;
}

int
print_my_timegrid(braid_App        app,
                  braid_Real      *ta,
                  braid_Int        ilower,
                  braid_Int        iupper)
{
   int   i;
   char  filename[255];
   FILE *file;

   /* filename could be anything that helps you track the current time grid */
   sprintf(filename, "timegrid.%d.%d.info", app->rank_t, app->rank_x);
   file = fopen(filename, "w");
   if (file != NULL) {
      for (i = ilower; i <= iupper; i++)
      {
         fprintf(file, "%d %f\n", i, ta[i-ilower]);
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
            braid_Int        ilower,
            braid_Int        iupper)
{
   double tstart;             /* time corresponding to ilower, i.e. lower time index value for this processor */
   int i;
   
   /* Start from the global tstart to compute the local tstart */
   tstart = app->tstart;
//   printf("tstart = %f\n", tstart);
   for (i = 0; i < ilower; i++)
   {
      tstart += app->dt[i];
//      printf("tstart = %f\n", tstart);
   }
   /* Assign time point values for local time point index values ilower:iupper */
   for (i = ilower; i <= iupper; i++)
   {
      ta[i-ilower]  = tstart;
//      printf("ta[%d] = %f\n", i-ilower, tstart);
      tstart       += app->dt[i];
   }
   print_my_timegrid(app, ta, ilower, iupper);

   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   if (t == 0.0)
   {
      /* Initial condition */
      (u->value) = 1.0;
   }
   else
   {
      /* Initialize all other time points */
      (u->value) = 0.456;//((double)rand()) / RAND_MAX;
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
   MPI_Comm   comm   = (app->comm);
   int        index, rank;
   char       filename[255];
   FILE      *file;
   double     t;
   
   braid_AccessStatusGetT(astatus,     &t);
   braid_AccessStatusGetIstop(astatus, &index);
   MPI_Comm_rank(comm, &rank);

   sprintf(filename, "%s.%07d.%05d", "ex-01.out", index, rank);
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

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
   braid_Core    core;
   my_App       *app;
   MPI_Comm      comm, comm_t, comm_x;
   double        tstart, tstop;
   int           ntime;

   int           max_levels = 1;
   int           nrelax     = 1;
   int           nrelax0    = -1;
   double        tol        = 1.0e-06;
   int           cfactor    = 2;
   int           max_iter   = 100;
   int           fmg        = 0;
   int           res        = 0;
   int           mydt       = 0;

   int           arg_index;
   int           ptx, pt, px;
   int           rank, rank_t, rank_x;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);

   /* default to 1 processor in space and time */
   ptx = pt = px = 1;

   /* ntime time intervals with spacing 1 */
   comm   = MPI_COMM_WORLD;
   MPI_Comm_rank(comm, &rank);
   ntime  = 32;
   tstart = 0.0;
   tstop  = tstart + ntime;
   
   /* Parse command line */

   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         if ( rank == 0 )
         {
            printf("\n");
            printf("  -pgrid  <pt px>   : set number of processors to be used in time and space\n");
            printf("  -ml  <max_levels> : set max levels\n");
            printf("  -nu  <nrelax>     : set num F-C relaxations\n");
            printf("  -nu0 <nrelax>     : set num F-C relaxations on level 0\n");
            printf("  -tol <tol>        : set stopping tolerance\n");
            printf("  -cf  <cfactor>    : set coarsening factor\n");
            printf("  -mi  <max_iter>   : set max iterations\n");
            printf("  -fmg              : use FMG cycling\n");
            printf("  -res              : use my residual\n");
            printf("  -dt <mydt>        : if specified, user supplied time grid is used as global fine grid, options in this example are\n");
            printf("                      1 - constant time step size dt\n");
            printf("                      2 - dt*0.5 for n = 1, ..., nt/2; dt*1.5 for n = nt/2+1, ..., nt\n");
            printf("\n");
         }
         exit(1);
      }
      else if ( strcmp(argv[arg_index], "-pgrid") == 0 ){
         arg_index++;
         pt = atoi(argv[arg_index++]);
         px = atoi(argv[arg_index++]);
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
      else if ( strcmp(argv[arg_index], "-dt") == 0 )
      {
         arg_index++;
         mydt = atoi(argv[arg_index++]);
      }
      else
      {
         arg_index++;
         /*break;*/
      }
   }

   /* Check the processor grid (px x pt = ptx?). */
   MPI_Comm_size( comm, &ptx );
   if( (px*pt) != ptx )
   {
       if( rank == 0 )
           printf("Error: px(%d) x pt(%d) does not equal the number of processors ptx(%d)!\n", px, pt, ptx);
       MPI_Finalize();
       return (0);
   }

   /* Create communicators for the time and space dimensions */
   braid_SplitCommworld(&comm, px, &comm_x, &comm_t);
   MPI_Comm_rank(comm_t, &rank_t);
   MPI_Comm_rank(comm_x, &rank_x);
   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->dt)     = NULL;
   (app->mydt)   = mydt;
   (app->comm)   = comm;
   (app->comm_t) = comm_t;
   (app->comm_x) = comm_x;
   (app->ptx)    = ptx;
   (app->pt)     = pt;
   (app->px)     = px;
   (app->tstart) = tstart;
   (app->tstop)  = tstop;
   (app->ntime)  = ntime;
   (app->rank)   = rank;
   (app->rank_t) = rank_t;
   (app->rank_x) = rank_x;

   braid_Init(comm, comm_t, tstart, tstop, ntime, app,
             my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
             my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

   braid_SetPrintLevel( core, 1);
   braid_SetMaxLevels(core, max_levels);
   braid_SetNRelax(core, -1, nrelax);
   if (nrelax0 > -1)
   {
      braid_SetNRelax(core,  0, nrelax0);
   }
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, cfactor);
   /*braid_SetCFactor(core,  0, 10);*/
   braid_SetMaxIter(core, max_iter);
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

   braid_Drive(core);

   braid_Destroy(core);
   MPI_Comm_free( &comm_x );
   MPI_Comm_free( &comm_t );

   /* Finalize MPI */
   MPI_Finalize();

   return (0);
}
