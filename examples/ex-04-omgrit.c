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
   double *design;      /* Holds time-dependent design (i.e. control) vector */
   double *gradient;    /* Holds the gradient vector */
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
   double *vec = (double*) malloc( size*sizeof(double) );
}

void
vec_destroy(double *vec)
{
   free(vec);
}

void
vec_copy(int size, double *invec, double *outvec)
{
   int i;
   for (i = 0; i < size; i++)
   {
      outvec[i] = invec[i];
   }
}

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

void
apply_PhiAdjoint(double dt, double *w)
{
   w[1] = w[1] - dt*w[1] + dt*w[0];
}

void
apply_Uinv(int size, double dt, double *u)
{
   u[0] /= 2*dt;
   u[1] /= 2*dt;
}

void
apply_Vinv(int size, double dt, double design, double *v)
{
   v[0] /= 2*dt;
}

void
apply_D(int size, double dt, double design, double *v)
{
}

void
apply_DAdjoint(int size, double dt, double design, double *w)
{
}

/*--------------------------------------------------------------------------
 * TriMGRIT wrapper routines
 *--------------------------------------------------------------------------*/

/* Compute A(u) - f */
my_TriResidual(braid_App       app,
               braid_Vector    wleft,
               braid_Vector    wright,
               braid_Vector    f,
               braid_Vector    r,
               braid_TriStatus status)
{
   int     index;
   double  t, tprev, tnext;
   double  dtprev, dtnext;
   double  design, gamma = (app->gamma);
   double *rtmp, *wtmp;
   
   /* Get the time-step size */
   braid_StepStatusGetTriT(status, &t, &tprev, &tnext);
   dtprev = t - tprev;
   dtnext = tnext - t;

   /* Get the current design from the app */
   braid_StepStatusGetTIndex(status, &index);
   design = (app->design)[index];

   /* Create temporary vectors */
   vec_create(2, &rtmp);
   vec_create(2, &wtmp);

   /* Initialize temporary residual vector */
   vec_copy(2, (r->values), rtmp);

   /* Compute action of center block */
   vec_copy(2, (r->values), wtmp);
   apply_DAdjoint(2, dtnext, design, wtmp);
   apply_Vinv(2, dtnext, gamma, wtmp);
   apply_D(2, dtnext, design, wtmp);
   vec_axpy(2, 1.0, wtmp, rtmp);

   vec_copy(2, (r->values), wtmp);
   apply_PhiAdjoint(2, dtnext, wtmp);
   apply_Uinv(2, dtprev, wtmp);
   apply_Phi(2, dtprev, wtmp);
   vec_axpy(2, 1.0, wtmp, rtmp);

   vec_copy(2, (r->values), wtmp);
   apply_Uinv(2, dtnext, wtmp);
   vec_axpy(2, 1.0, wtmp, rtmp);

   /* Compute action of west block */
   vec_copy(2, (wleft->values), wtmp);
   apply_Uinv(2, dtprev, wtmp);
   apply_Phi(2, dtprev, wtmp);
   vec_axpy(2, 1.0, wtmp, rtmp);
   
   /* Compute action of east block */
   vec_copy(2, (wright->values), wtmp);
   apply_PhiAdjoint(2, dtnext, wtmp);
   apply_Uinv(2, dtnext, wtmp);
   vec_axpy(2, 1.0, wtmp, rtmp);

   /* Subtract rhs f */
   if (f != NULL)
   {
      vec_axpy(2, -1.0, (f->values), rtmp);
   }
   
   /* Copy temporary residual vector into residual */
   vec_copy(2, rtmp, (r->values));
   
   /* Destroy temporary vectors */
   vec_destroy(rtmp);
   vec_destroy(wtmp);
   
   return 0;
}   


/* Solve A(u) = f */
my_TriSolve(braid_App       app,
            braid_Vector    wleft,
            braid_Vector    wright,
            braid_Vector    f,
            braid_Vector    w,
            braid_TriStatus status)
{
   int     index;
   double  t, tprev, tnext;
   double  dtprev, dtnext;
   double  design, gamma = (app->gamma);
   double *wtmp, *rtmp;
   
   /* Get the time-step size */
   braid_StepStatusGetTriT(status, &t, &tprev, &tnext);
   dtprev = t - tprev;
   dtnext = tnext - t;

   /* Get the current design from the app */
   braid_StepStatusGetTIndex(status, &index);
   design = (app->design)[index];

   /* Create temporary vector */
   vec_create(2, &wtmp);

   /* Initialize temporary solution vector */
   vec_copy(2, (w->values), wtmp);
   
   /* Compute residual */
   my_TriResidual(app, wleft, wright, f, w, status);

   /* Apply center block preconditioner */
   rtmp = (w->values);
   rtmp[0] = -rtmp[0]*dtnext;
   rtmp[1] = -rtmp[1]/(1/dtnext + dt/(2*gamma));

   /* Complete residual update */
   vec_axpy(2, 1.0, wtmp, (w->values));
   
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

   /* Allocate the vector */
   u = (my_Vector *) malloc(sizeof(my_Vector));
   u->values = (double*) malloc( 2*sizeof(double) );

   /* Initialize the vector */
   if (t == 0.0)
   {
      u->values[0] = 0.0;
      u->values[1] = -1.0;
   }
   else
   {
      u->values[0] = 0.0;
      u->values[1] = 0.0;
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

   /* Allocate the vector */
   v = (my_Vector *) malloc(sizeof(my_Vector));
   v->values = (double*) malloc( 2*sizeof(double) );

   /* Clone the values */
   v->values[0] = u->values[0];
   v->values[1] = u->values[1];

   *v_ptr = v;

   return 0;
}


int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u->values);
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

   (y->values)[0] = alpha*(x->values)[0] + beta*(y->values)[0];
   (y->values)[1] = alpha*(x->values)[1] + beta*(y->values)[1];

   return 0;
}


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


int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
   *size_ptr = 2*sizeof(double);
   return 0;
}


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

/* Evaluate one term of the time-dependent discretized objective function for
 * vector u.  The result over all time values will be summed by Braid. */
int 
my_ObjectiveT(braid_App              app,
              braid_Vector           u,
              braid_ObjectiveStatus  ostatus,
              double                *objectiveT_ptr)
{
   double objT;
   double design;
   int    index;
   double deltaT = 1./app->ntime;
   double gamma  = app->gamma;

   /* Get the time index*/
   braid_ObjectiveStatusGetTIndex(ostatus, &index);

   /* Evaluate the objective function after the first step */
   if ( index > 0)
   {
      /* Get the design from the app */
      design = app->design[index-1];

      /* Evaluate objective */
      objT = evalObjectiveT( u->values, design, deltaT, gamma);
   }
   else
   {
      objT = 0.0;
   }

   *objectiveT_ptr = objT;
   
   return 0;
}


/* Transposed partial derivatives of objectiveT */ 
int
my_ObjectiveT_diff(braid_App            app,
                  braid_Vector          u,
                  braid_Vector          u_bar,
                  braid_Real            F_bar,
                  braid_ObjectiveStatus ostatus)
{
   int     index;
   double  design;
   double  gamma   = app->gamma;
   double  deltaT  = 1. / app->ntime;
   
   /* Get the design from the app */
   braid_ObjectiveStatusGetTIndex(ostatus, &index);

   if ( index > 0 )
   {
      /* Get the design from the app */
      design = app->design[index-1];

      /* Partial derivatives of objective */
      app->gradient[index-1] += evalObjectiveT_diff(u_bar->values, u->values, design, gamma, deltaT);
   }
   
   return 0;
}

/* Transposed partial derivatives of step times u_bar */
int
my_Step_diff(braid_App              app,
                braid_Vector        ustop,
                braid_Vector        u,
                braid_Vector        ustop_bar,
                braid_Vector        u_bar,
                braid_StepStatus    status)
{

   double  tstop, tstart;
   int     tidx;

   /* Get time and time index  */
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   braid_StepStatusGetTIndex(status, &tidx);
   double deltaT = tstop - tstart;

   /* transposed derivative of take_step times u_bar */
   app->gradient[tidx] += take_step_diff(u_bar->values, deltaT);

   return 0;
}

/* Set the gradient to zero */
int 
my_ResetGradient(braid_App app)
{
   int ts;

   for(ts = 0; ts < app->ntime; ts++) 
   {
      app->gradient[ts] = 0.0;
   }

   return 0;
}



/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/
int main (int argc, char *argv[])
{

   braid_Core core;
   my_App     *app;
         
   double   tstart, tstop; 
   int      rank, ntime, ts, iter, maxiter, nreq, arg_index;
   double  *design; 
   double  *gradient; 
   double   objective, gamma, stepsize, mygnorm, gnorm, gtol, rnorm, rnorm_adj;
   int      max_levels, cfactor, access_level, print_level, braid_maxiter;
   double   braid_tol, braid_adjtol;

   /* Define time domain */
   ntime  = 20;              /* Total number of time-steps */
   tstart = 0.0;             /* Beginning of time domain */
   tstop  = 1.0;             /* End of time domain*/

   /* Define some optimization parameters */
   gamma    = 0.005;         /* Relaxation parameter in the objective function */
   stepsize = 50.0;          /* Step size for design updates */
   maxiter  = 500;           /* Maximum number of optimization iterations */
   gtol     = 1e-6;          /* Stopping criterion on the gradient norm */

   /* Define some Braid parameters */
   max_levels     = 4;
   braid_maxiter  = 10;
   cfactor        = 2;
   braid_tol      = 1.0e-6;
   braid_adjtol   = 1.0e-6;
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
         printf("  -stepsize <stepsize>    : Step size for design updates \n");
         printf("  -mi <maxiter>           : Maximum number of optimization iterations \n");
         printf("  -gtol <gtol>            : Stopping criterion on the gradient norm \n");
         printf("  -ml <max_levels>        : Max number of braid levels \n");
         printf("  -bmi <braid_maxiter>    : Braid max_iter \n");
         printf("  -cf <cfactor>           : Coarsening factor \n");
         printf("  -btol <braid_tol>       : Braid halting tolerance \n");
         printf("  -batol <braid_adjtol>   : Braid adjoint halting tolerance \n");
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
      else if ( strcmp(argv[arg_index], "-stepsize") == 0 )
      {
         arg_index++;
         stepsize = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mi") == 0 )
      {
         arg_index++;
         maxiter = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-gtol") == 0 )
      {
         arg_index++;
         gtol = atof(argv[arg_index++]);
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
      else if ( strcmp(argv[arg_index], "-batol") == 0 )
      {
         arg_index++;
         braid_adjtol = atof(argv[arg_index++]);
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

   /* Initialize optimization */
   design   = (double*) malloc( ntime*sizeof(double) );    /* design vector (control c) */
   gradient = (double*) malloc( ntime*sizeof(double) );    /* gradient vector */
   for (ts = 0; ts < ntime; ts++)
   {
      design[ts]   = 0.;
      gradient[ts] = 0.;
   }

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   /* Set up the app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->myid     = rank;
   app->ntime    = ntime;
   app->design   = design;
   app->gradient = gradient;
   app->gamma    = gamma;

   /* Initialize XBraid */
   braid_InitTriMGRIT(MPI_COMM_WORLD, MPI_COMM_WORLD, tstart, tstop, ntime, app,
                      my_TriSolve, my_TriResidual, my_Clone, my_Free, my_Sum, my_SpatialNorm,
                      my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

   /* Set some XBraid(_Adjoint) parameters */
   braid_SetMaxLevels(core, max_levels);
   braid_SetCFactor(core, -1, cfactor);
   braid_SetAccessLevel(core, access_level);
   braid_SetPrintLevel( core, print_level);       
   braid_SetMaxIter(core, braid_maxiter);
   braid_SetAbsTol(core, braid_tol);
   braid_SetAbsTolAdjoint(core, braid_adjtol);

   /* Prepare optimization output */
   if (rank == 0)
   {
      printf("\nOptimization:         || r ||        || r_adj ||        Objective           || Gradient ||\n");
   }


   /* Optimization iteration */
   for (iter = 0; iter < maxiter; iter++)
   {

      /* Parallel-in-time simulation and gradient computation */
      braid_Drive(core);

      /* Get objective function value */
      nreq = -1;
      braid_GetObjective(core, &objective);

      /* Get the state and adjoint residual norms */
      braid_GetRNorms(core, &nreq, &rnorm);
      braid_GetRNormAdjoint(core, &rnorm_adj);

      /* Compute the norm of the gradient */
      mygnorm = compute_sqnorm(app->gradient, ntime);
      MPI_Allreduce(&mygnorm, &gnorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      gnorm = sqrt(gnorm);

      /* Output */
      if (rank == 0)
      {
         printf("Optimization: %3d  %1.8e  %1.8e  %1.14e  %1.14e\n", iter, rnorm, rnorm_adj, objective, gnorm);
      }

      /* Check optimization convergence */
      if (gnorm < gtol)
      {
         break;
      }

      /* Design update */
      for(ts = 0; ts < ntime; ts++) 
      {
         app->design[ts] -= stepsize * app->gradient[ts];
      }

   }

   
   /* Output */
   if (rank == 0)
   {
      if (iter == maxiter)
      {
         printf("\n Max. number of iterations reached! \n\n"); 
      }
      else
      {
         /* Print some statistics about the optimization run */
         printf("\n");
         printf("  Optimization has converged.\n");
         printf("\n"); 
         printf("  Objective function value = %1.8e\n", objective);
         printf("  Gradient norm            = %1.8e\n", gnorm);
         printf("\n");
         printf("  optimization iterations  = %d\n", iter);
         printf("  max optim iterations     = %d\n", maxiter);
         printf("  gradient norm tolerance  = %1.1e\n", gtol);
         printf("\n");
      }
   }
   braid_PrintStats(core);



   /* Write final design to file */
   write_design_vec("design", design, ntime);

   /* Clean up */
   free(design);
   free(gradient);
   free(app);
   
   braid_Destroy(core);
   MPI_Finalize();

   return (0);
}
