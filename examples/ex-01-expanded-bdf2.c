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
 * Example:       ex-01-expanded-bdf2.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-01-expanded-bdf2
 *
 * Help with:     ex-01-expanded-bdf2 -help
 *
 * Sample run:    mpirun -np 2 ex-01-expanded-bdf2
 *
 * Description:   solve the scalar ODE 
 *                   u' = lambda u, 
 *                   with lambda=-1 and y(0) = 1
 *
 *                Similar to ex-01-expanded, but demonstrates a multistep integration scheme
 *                using the shellvector feature to store time value information at all points
 *                in time. However for simplicity, some of the extras from ex-01-expanded 
 *                (like my_residual and my_timegrid) are not presented here.  The purpose of
 *                this example is to show a way to solve with BDF2.
 *                
 *                For the implicit BDF2 scheme used, the solve is hard-coded for this
 *                specific equation in the helper function BDF2_Expo1D.
 *                
 *                When run with the default 10 time steps, the solution is:
 *                $ ./ex-01-expanded-bdf2
 *                $ cat ex-01-expanded-bdf2.out.00*
 *                  0.00000000000000e+00 1.00000000000000e+00
 *                  5.00000000000000e-01 6.06530659712633e-01
 *                  1.00000000000000e+00 3.56530659712633e-01
 *                  1.50000000000000e+00 2.04897994784475e-01
 *                  2.00000000000000e+00 1.15765329856317e-01
 *                  2.50000000000000e+00 6.45408311601979e-02
 *                  3.00000000000000e+00 3.55994986961188e-02
 *                  3.50000000000000e+00 1.94642909060693e-02
 *                  4.00000000000000e+00 1.05644162320396e-02
 *                  4.50000000000000e+00 5.69834350552227e-03
 *                  5.00000000000000e+00 3.05723944751237e-03
 *                  5.50000000000000e+00 1.63265357113180e-03
 **/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.h"

/*--------------------------------------------------------------------------
 * User-defined routines and structures
 *--------------------------------------------------------------------------*/

/* Define ODE coefficient */
double lambda = -1.0;

/* App structure can contain anything, and be named anything as well */
typedef struct _braid_App_struct
{
   MPI_Comm  comm;
   double    tstart;
   double    tstop;
   int       ntime;
   int       rank;
} my_App;

/* Vector structure can contain anything, and be name anything as well */
typedef struct _braid_Vector_struct
{
   // The Shell part consistes of the time-step values, tprev[0]=t_{n-1}, tprev[1]=t_{n-2}
   double* tprev;

   // The Vector part consistes of the solution values, yprev[0]=y_{n-1}, yprev[1]=y_{n-2}
   double* yprev;

} my_Vector;

/* Helper function for my_Step */
int
BDF2_Expo1D(double  t0,
            double  t1,
            double  t2,
            double  y1,
            double  y2,
            double* y0)
{
   double h0 = t0-t1;
   double rho = h0/(t1-t2);
   double c0 = (2*rho+1)/(rho+1)-lambda*h0;
   
   *y0 = (rho+1)/c0*y1-(rho*rho/(rho+1))/c0*y2;
   return 0;
}

int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   double new_y2, new_y1;

   /* Solve for the new vector, placing it in ustop.
    * We use the fact that each braid_Vector always contains at least the
    * Shell information containing time values.  Thus, there is no need to query 
    * StepStatus for tstart and tstop. */
   BDF2_Expo1D(ustop->tprev[1], u->tprev[0], u->tprev[1], u->yprev[0], u->yprev[1], &new_y2);
   BDF2_Expo1D(ustop->tprev[0], ustop->tprev[1], u->tprev[0], new_y2, u->yprev[0], &new_y1);

   /* Write the updated solution into the output vector, u*/
   u->tprev[0] = ustop->tprev[0];
   u->tprev[1] = ustop->tprev[1];
   u->yprev[0] = new_y1;
   u->yprev[1] = new_y2;

   /* no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;
   
   u = (my_Vector *) malloc(sizeof(my_Vector));
   u->tprev = (double*)malloc(2*sizeof(double));
   u->yprev = (double*)malloc(2*sizeof(double));

   double dt = (app->tstop-app->tstart)/app->ntime;
   
   /* Initialize time values */
   u->tprev[1] = t;
   u->tprev[0] = t+0.5*dt;

   /* Initialize solution values */
   if (t == 0.0) /* Initial condition, note how it now contains two time values */
   {
      (u->yprev[1]) = 1.0;
      (u->yprev[0]) = exp(lambda*0.5*dt);
   }
   else /* All other time points set to arbitrary value */
   {
      (u->yprev[1]) = 0.456;
      (u->yprev[0]) = 0.456;
   }

   *u_ptr = u;
   return 0;
}

/* Initialization of the shell part only (i.e., initialize the time-values in
 * tprev, not solution values in yprev */
int
my_InitShell(braid_App     app,
             double        t,
             braid_Vector *u_ptr)
{
   my_Vector *u;
   
   u = (my_Vector *) malloc(sizeof(my_Vector));
   u->tprev = (double*)malloc(2*sizeof(double));
   u->yprev = NULL;
 
   double dt = (app->tstop-app->tstart)/app->ntime;
   u->tprev[1] = t;
   u->tprev[0] = t + 0.5*dt;

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
   v->tprev = (double*)malloc(2*sizeof(double));
   v->yprev = (double*)malloc(2*sizeof(double));

   v->tprev[0] = u->tprev[0];
   v->tprev[1] = u->tprev[1];
   v->yprev[0] = u->yprev[0];
   v->yprev[1] = u->yprev[1];

   *v_ptr = v;
   return 0;
}

/* Cloning only the shell part of the vector*/
int
my_CloneShell(braid_App     app,
              braid_Vector  u,
              braid_Vector *v_ptr)
{
   my_Vector *v;

   v = (my_Vector *) malloc(sizeof(my_Vector));
   v->tprev = (double*)malloc(2*sizeof(double));
   v->yprev = NULL;

   v->tprev[0] = u->tprev[0];
   v->tprev[1] = u->tprev[1];

   *v_ptr = v;
   return 0;
}

int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u->tprev);
   free(u->yprev);
   u->tprev = NULL;
   u->yprev = NULL;
   free(u);

   return 0;
}

/* Free only the non-shell part of the vector, but keep the shell part*/
int
my_FreeShell(braid_App    app,
             braid_Vector u)
{
   if (u->yprev != NULL)
   {
      free(u->yprev);
   }
   u->yprev = NULL;

   return 0;
}

int
my_Sum(braid_App     app,
       double        alpha,
       braid_Vector  x,
       double        beta,
       braid_Vector  y)
{
   y->yprev[0] = alpha*x->yprev[0] + beta*y->yprev[0];
   y->yprev[1] = alpha*x->yprev[1] + beta*y->yprev[1];

   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   double sqdot = pow(u->yprev[0],2) + pow(u->yprev[1],2);
   *norm_ptr = sqrt(sqdot);
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
   sprintf(filename, "%s.%04d.%03d", "ex-01-expanded-bdf2.out", index, app->rank);
   file = fopen(filename, "w");
   fprintf(file,"%.14e %.14e\n%.14e %.14e\n",u->tprev[1],u->yprev[1],u->tprev[0],u->yprev[0]);
   fflush(file);
   fclose(file);
   
   return 0;
}

int
my_BufSize(braid_App  app,
           int       *size_ptr,
           braid_BufferStatus  bstatus)
{
   *size_ptr = 4*sizeof(double);
   return 0;
}

int
my_BufPack(braid_App     app,
           braid_Vector  u,
           void         *buffer,
           braid_BufferStatus  bstatus)
{
   double* dbuffer=(double*)buffer;

   dbuffer[0] = u->tprev[1];
   dbuffer[1] = u->yprev[1];
   dbuffer[2] = u->tprev[0];
   dbuffer[3] = u->yprev[0];

   braid_BufferStatusSetSize(bstatus, 4*sizeof(double));
   return 0;
}

int
my_BufUnpack(braid_App           app,
             void                *buffer,
             braid_Vector        *u_ptr,
             braid_BufferStatus  status)
{
   double* dbuffer = (double*)buffer;

   my_Vector *u;
   u = (my_Vector *) malloc(sizeof(my_Vector));
   u->tprev = (double*)malloc(2*sizeof(double));
   u->yprev = (double*)malloc(2*sizeof(double));

   u->tprev[1] = dbuffer[0];
   u->yprev[1] = dbuffer[1];
   u->tprev[0] = dbuffer[2];
   u->yprev[0] = dbuffer[3];

   *u_ptr = u;
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

   int           max_levels = 2;
   int           nrelax     = 1;
   int           nrelax0    = -1;
   double        tol        = 1.0e-06;
   int           cfactor    = 2;
   int           max_iter   = 100;
   int           fmg        = 0;

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
            printf("\nExample 1: Solve a scalar ODE, using BDF2 \n\n");
            printf("  -ntime <ntime>    : set num time points\n");
            printf("  -ml  <max_levels> : set max levels\n");
            printf("  -nu  <nrelax>     : set num F-C relaxations\n");
            printf("  -nu0 <nrelax>     : set num F-C relaxations on level 0\n");
            printf("  -tol <tol>        : set stopping tolerance\n");
            printf("  -cf  <cfactor>    : set coarsening factor\n");
            printf("  -mi  <max_iter>   : set max iterations\n");
            printf("  -fmg              : use FMG cycling\n\n");
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
      else
      {
         arg_index++;
         /*break;*/
      }
   }

   // Time points are lumped into pairs of "2" for BDF 2
   ntime = (ntime+1)/2;

   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->tstart) = tstart;
   (app->tstop)  = tstop;
   (app->ntime)  = ntime;
   (app->rank)   = rank;
   (app->comm)   = comm;

   /* initialize XBraid and set options */
   braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
             my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
             my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
   
   braid_SetPrintLevel(core, 2);
   braid_SetMaxLevels(core, max_levels);
   braid_SetNRelax(core, -1, nrelax);
   if (nrelax0 > -1)
   {
      braid_SetNRelax(core,  0, nrelax0);
   }
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, cfactor);
   braid_SetMaxIter(core, max_iter);
   if (fmg)
   {
      braid_SetFMG(core);
   }

   braid_SetShell(core, my_InitShell, my_CloneShell, my_FreeShell);
   braid_SetStorage(core, -1);

   braid_Drive(core);
   braid_Destroy(core);
   // check TWO regression tests, harmonize the F90 file and check its parameters
   free(app);
   
   MPI_Finalize();
   return 0;
}
