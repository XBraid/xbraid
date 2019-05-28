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
 * Example:       ex-04.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-04-adjoint
 *
 * Description:  Solves a simple optimal control problem in time-parallel:
 * 
 *                 min   \int_0^1 u_1(t)^2 + u_2(t)^2 + gamma c(t)^2  dt
 * 
 *                  s.t.  d/dt u_1(t) = u_2(t)
 *                        d/dt u_2(t) = -u_2(t) + c(t)
 * 
 *               with initial condition u_1(0) = 0, u_2(0) = -1
 *               and piecewise constant control c(t).  
 *
 *               Implements a steepest-descent optimization iteration
 *               using fixed step size for design updates.   
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "braid.h"
#include "braid_test.h"
#include "ex-04-lib.c"

/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/

typedef struct _braid_App_struct
{
   int     myid;        /* Rank of the processor */
   double  gamma;       /* Relaxation parameter for objective function */
   int     ntime;       /* Total number of time-steps */
} my_App;


/* Define the state vector at one time-step */
typedef struct _braid_Vector_struct
{
   double *values;     /* Holds the R^2 state vector (u_1, u_2) */

} my_Vector;

/*--------------------------------------------------------------------------
 * Vector utility routines
 *--------------------------------------------------------------------------*/

void
vec_create(int size, double **vec_ptr)
{
   *vec_ptr = (double*) malloc( size*sizeof(double) );
}

void
vec_destroy(double *vec)
{
   free(vec);
}

/*------------------------------------*/

void
vec_copy(int size, double *invec, double *outvec)
{
   int i;
   for (i = 0; i < size; i++)
   {
      outvec[i] = invec[i];
   }
}

/*------------------------------------*/

void
vec_axpy(int size, double alpha, double *x, double *y)
{
   int i;
   for (i = 0; i < size; i++)
   {
      y[i] = y[i] + alpha*x[i];
   }
}

/*--------------------------------------------------------------------------
 * KKT component routines
 *--------------------------------------------------------------------------*/

void
apply_Phi(double dt, double *u)
{
   u[0] = u[0] + dt*u[1];
   u[1] = u[1] - dt*u[1];
}

/*------------------------------------*/

void
apply_PhiAdjoint(double dt, double *w)
{
   w[1] = w[1] - dt*w[1] + dt*w[0];
}

/*------------------------------------*/

void
apply_Uinv(double dt, double *u)
{
   u[0] /= 2*dt;
   u[1] /= 2*dt;
}

/*------------------------------------*/

void
apply_Vinv(double dt, double gamma, double *v)
{
   v[0] /= 2*gamma*dt;
}

/*------------------------------------*/

/* Maps v to w (which may be the same pointer) */
void
apply_D(double dt, double *v, double *w)
{
   double design = v[0];
   w[0] = 0.0;
   w[1] = dt*design;
}

/*------------------------------------*/

/* Maps w to v (which may be the same pointer) */
void
apply_DAdjoint(double dt, double *w, double *v)
{
   double design = dt*w[1];
   v[0] = design;
}

/*--------------------------------------------------------------------------
 * TriMGRIT wrapper routines
 *--------------------------------------------------------------------------*/

/* Compute A(u) - f */
int
my_TriResidual(braid_App       app,
               braid_Vector    uleft,
               braid_Vector    uright,
               braid_Vector    f,
               braid_Vector    r,
               braid_TriStatus status)
{
   double  t, tprev, tnext;
   double  dtprev, dtnext;
   double  gamma = (app->gamma);
   double *rtmp, *utmp;
   
   /* Get the time-step size */
   braid_TriStatusGetTriT(status, &t, &tprev, &tnext);
   dtprev = t - tprev;
   dtnext = tnext - t;

   /* Create temporary vectors */
   vec_create(2, &rtmp);
   vec_create(2, &utmp);

   /* Initialize temporary residual vector */
   vec_copy(2, (r->values), rtmp);

   /* Compute action of center block */
   vec_copy(2, (r->values), utmp);
   apply_DAdjoint(dtnext, utmp, utmp);
   apply_Vinv(dtnext, gamma, utmp);
   apply_D(dtnext, utmp, utmp);
   vec_axpy(2, 1.0, utmp, rtmp);

   vec_copy(2, (r->values), utmp);
   apply_PhiAdjoint(dtnext, utmp);
   apply_Uinv(dtprev, utmp);
   apply_Phi(dtprev, utmp);
   vec_axpy(2, 1.0, utmp, rtmp);

   vec_copy(2, (r->values), utmp);
   apply_Uinv(dtnext, utmp);
   vec_axpy(2, 1.0, utmp, rtmp);

   /* Compute action of west block */
   vec_copy(2, (uleft->values), utmp);
   apply_Uinv(dtprev, utmp);
   apply_Phi(dtprev, utmp);
   vec_axpy(2, 1.0, utmp, rtmp);
   
   /* Compute action of east block */
   vec_copy(2, (uright->values), utmp);
   apply_PhiAdjoint(dtnext, utmp);
   apply_Uinv(dtnext, utmp);
   vec_axpy(2, 1.0, utmp, rtmp);

   /* Subtract rhs f */
   if (f != NULL)
   {
      vec_axpy(2, -1.0, (f->values), rtmp);
   }
   
   /* Copy temporary residual vector into residual */
   vec_copy(2, rtmp, (r->values));
   
   /* Destroy temporary vectors */
   vec_destroy(rtmp);
   vec_destroy(utmp);
   
   return 0;
}   

/*------------------------------------*/

/* Solve A(u) = f */
int
my_TriSolve(braid_App       app,
            braid_Vector    uleft,
            braid_Vector    uright,
            braid_Vector    f,
            braid_Vector    u,
            braid_TriStatus status)
{
   double  t, tprev, tnext;
   double  dtnext;
   double  gamma = (app->gamma);
   double *utmp, *rtmp;
   
   /* Get the time-step size */
   braid_TriStatusGetTriT(status, &t, &tprev, &tnext);
   dtnext = tnext - t;

   /* Create temporary vector */
   vec_create(2, &utmp);

   /* Initialize temporary solution vector */
   vec_copy(2, (u->values), utmp);
   
   /* Compute residual */
   my_TriResidual(app, uleft, uright, f, u, status);

   /* Apply center block preconditioner */
   rtmp = (u->values);
   rtmp[0] = -rtmp[0]*dtnext;
   rtmp[1] = -rtmp[1]/(1/dtnext + dtnext/(2*gamma));

   /* Complete residual update */
   vec_axpy(2, 1.0, utmp, (u->values));
   
   /* no refinement */
   braid_TriStatusSetRFactor(status, 1);

   return 0;
}   

/*------------------------------------*/

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{

   my_Vector *u;

   /* Allocate the vector */
   u = (my_Vector *) malloc(sizeof(my_Vector));
   vec_create(2, &(u->values));

   /* Initialize the vector */
   if (t == 0.0)
   {
      /* This is -g because we are solving the Schur-complement system */
      u->values[0] = 0.0;
      u->values[1] = 1.0;
   }
   else
   {
      u->values[0] = 0.0;
      u->values[1] = 0.0;
   }

   *u_ptr = u;

   return 0;
}

/*------------------------------------*/

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   my_Vector *v;

   /* Allocate the vector */
   v = (my_Vector *) malloc(sizeof(my_Vector));
   vec_create(2, &(v->values));

   /* Clone the values */
   v->values[0] = u->values[0];
   v->values[1] = u->values[1];

   *v_ptr = v;

   return 0;
}

/*------------------------------------*/

int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u->values);
   free(u);

   return 0;
}

/*------------------------------------*/

int
my_Sum(braid_App     app,
       double        alpha,
       braid_Vector  x,
       double        beta,
       braid_Vector  y)
{

   (y->values)[0] = alpha*(x->values)[0] + beta*(y->values)[0];
   (y->values)[1] = alpha*(x->values)[1] + beta*(y->values)[1];

   return 0;
}

/*------------------------------------*/

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   int i;
   double dot = 0.0;

   for (i = 0; i < 2; i++)
   {
      dot += (u->values)[i]*(u->values)[i];
   }
   *norm_ptr = sqrt(dot);

   return 0;
}

/*------------------------------------*/

// ZTODO: Need to compute u from adjoint and it reqires communication

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int   done, index;
   char  filename[255];
   FILE * file;

   /* Print solution to file if simulation is over */
   braid_AccessStatusGetDone(astatus, &done);

   if (done)
   {
      braid_AccessStatusGetTIndex(astatus, &index);
      sprintf(filename, "%s.%04d.%03d", "ex-04.out", index, app->myid);
      file = fopen(filename, "w");
      fprintf(file, "%1.14e, %1.14e\n", (u->values)[0], (u->values)[1]);
      fflush(file);
      fclose(file);
   }


   return 0;
}

/*------------------------------------*/

int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
   *size_ptr = 2*sizeof(double);
   return 0;
}

/*------------------------------------*/

int
my_BufPack(braid_App           app,
           braid_Vector        u,
           void               *buffer,
           braid_BufferStatus  bstatus)
{
   double *dbuffer = buffer;
   int i;

   for(i = 0; i < 2; i++)
   {
      dbuffer[i] = (u->values)[i];
   }

   braid_BufferStatusSetSize( bstatus,  2*sizeof(double));

   return 0;
}

/*------------------------------------*/

int
my_BufUnpack(braid_App           app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus  bstatus)
{
   my_Vector *u = NULL;
   double    *dbuffer = buffer;
   int i;

   /* Allocate memory */
   u = (my_Vector *) malloc(sizeof(my_Vector));
   u->values = (double*) malloc( 2*sizeof(double) );

   /* Unpack the buffer */
   for(i = 0; i < 2; i++)
   {
      (u->values)[i] = dbuffer[i];
   }

   *u_ptr = u;
   return 0;
}

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
   braid_Core  core;
   my_App     *app;
         
   double      tstart, tstop; 
   int         rank, ntime, arg_index;
   double      gamma;
   int         max_levels, cfactor, access_level, print_level, braid_maxiter;
   double      braid_tol;

   /* Define time domain */
   ntime  = 20;              /* Total number of time-steps */
   tstart = 0.0;             /* Beginning of time domain */
   tstop  = 1.0;             /* End of time domain*/

   /* Define some optimization parameters */
   gamma = 0.005;            /* Relaxation parameter in the objective function */

   /* Define some Braid parameters */
   max_levels     = 4;
   braid_maxiter  = 10;
   cfactor        = 2;
   braid_tol      = 1.0e-6;
   access_level   = 0;
   print_level    = 0;

   /* Parse command line */
   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         printf("\n");
         printf(" Solves a simple optimal control problem in time-serial on [0, 1] \n\n");
         printf("  min   \\int_0^1 u_1(t)^2 + u_2(t)^2 + gamma c(t)^2  dt \n\n");
         printf("  s.t.  d/dt u_1(t) = u_2(t) \n");
         printf("        d/dt u_2(t) = -u_2(t) + c(t) \n\n");
         printf("  -ntime <ntime>          : set num points in time\n");
         printf("  -gamma <gamma>          : Relaxation parameter in the objective function \n");
         printf("  -ml <max_levels>        : Max number of braid levels \n");
         printf("  -bmi <braid_maxiter>    : Braid max_iter \n");
         printf("  -cf <cfactor>           : Coarsening factor \n");
         printf("  -btol <braid_tol>       : Braid halting tolerance \n");
         printf("  -access <access_level>  : Braid access level \n");
         printf("  -print <print_level>    : Braid print level \n");
         exit(1);
      }
      else if ( strcmp(argv[arg_index], "-ntime") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-gamma") == 0 )
      {
         arg_index++;
         gamma = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ml") == 0 )
      {
         arg_index++;
         max_levels = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-bmi") == 0 )
      {
         arg_index++;
         braid_maxiter = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-cf") == 0 )
      {
         arg_index++;
         cfactor = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-btol") == 0 )
      {
         arg_index++;
         braid_tol = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-access") == 0 )
      {
         arg_index++;
         access_level = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-print") == 0 )
      {
         arg_index++;
         print_level = atoi(argv[arg_index++]);
      }
      else
      {
         printf("ABORTING: incorrect command line parameter %s\n", argv[arg_index]);
         return (0);
      }
   }

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   /* Set up the app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->myid     = rank;
   app->ntime    = ntime;
   app->gamma    = gamma;

   /* Initialize XBraid */
   braid_InitTriMGRIT(MPI_COMM_WORLD, MPI_COMM_WORLD, tstart, tstop, ntime, app,
                      my_TriResidual, my_TriSolve, my_Init, my_Clone, my_Free,
                      my_Sum, my_SpatialNorm, my_Access,
                      my_BufSize, my_BufPack, my_BufUnpack, &core);

   /* Set some XBraid(_Adjoint) parameters */
   braid_SetMaxLevels(core, max_levels);
   braid_SetCFactor(core, -1, cfactor);
   braid_SetAccessLevel(core, access_level);
   braid_SetPrintLevel( core, print_level);       
   braid_SetMaxIter(core, braid_maxiter);
   braid_SetAbsTol(core, braid_tol);

   /* Parallel-in-time TriMGRIT simulation */
   braid_Drive(core);




//    /* Get the state and adjoint residual norms */
//    braid_GetRNorms(core, &nreq, &rnorm);
//    braid_GetRNormAdjoint(core, &rnorm_adj);
   
//   /* Output */
//   if (rank == 0)
//   {
//      /* Print some statistics about the optimization run */
//      printf("\n");
//      printf("  Optimization has converged.\n");
//      printf("\n"); 
//      printf("  Objective function value = %1.8e\n", objective);
//      printf("\n");
//   }
   braid_PrintStats(core);

   free(app);
   
   braid_Destroy(core);
   MPI_Finalize();

   return (0);
}
