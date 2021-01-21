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

 /** Driver:       drive-simple-optimization.c
 *
 * Interface:    C
 * 
 * Requires:     only C-language support     
 *
 * Compile with: make drive-simple-optimization
 *
 * Description:  Uses XBraid to solve the adjoint equation of the following 
 *               simple optimal control problem
 * 
 *                 min   \int_0^1 u_1(t)^2 + u_2(t)^2 dt 
 * 
 *                  s.t.  d/dt u_1(t) = u_2(t), 
 *                        d/dt u_2(t) = -u_2(t) + c(t)^2
 * 
 *               with initial condition u_1(0) = 0, u_2(0) = -1 and piecewise
 *               constant control c(t).  
 *
 *               This example is similar to examples/ex-04.c, but instead of
 *               using the XBraid_adjoint feature to compute the gradient, 
 *               it uses XBraid to solve the adjoint equation directly. 
 *               It builds two separate XBraid hierarchies -- one for the 
 *               primal forward pass and another for the adjoint backwards 
 *               pass.  This example additionally allows for a distributed 
 *               storage of the control variable across processors in time.  
 **/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "braid.h"
#include "_braid.h"
#include "braid_test.h"
#include "drive-simple-optimization-lib.c"


/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/

typedef struct _braid_App_struct
{
   int     myid;        /* Rank of the processor */
   int     ntime;       /* Total number of time-steps */
   int     tstop;        /* Final time */
   double* design;       /* Local design vector (local on this processor ) */
   int     ilower;       /* index of first time point on this processor */

   braid_Core primalcore;  /* XBraid core for solving the forward pass (ODE) */
} my_App;


/* Define the state vector at one time-step */
typedef struct _braid_Vector_struct
{
   /* The Vector part holds the R^2 state vector (u_1, u_2) at one time step */
   double *values;    

   /* The Shell part consists of the design variable c at one time step */
   double  design;
   /* Store the local gradient in u. */
   double gradient;

} my_Vector;


/* Get the index of the local design on this processor */
int GetLocalDesignIndex(braid_App app,
                        int       ts)
{
        return ts - app->ilower;
}

/* Get the time step index */
int GetTimeStepIndex(braid_App app,
                     double    t)
{
   int ts = t * app->ntime / app->tstop; 
   return ts;
}  

int GetPrimalTimeIndex(braid_App app,
                       int       ts)
{
   int idx = app->ntime - ts;
   return idx;
}
/*--------------------------------------------------------------------------
 * Integration routines
 *--------------------------------------------------------------------------*/

/* Take one step forward from tstart to tstop */
int 
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   int    ts_start, ts_stop;
   double tstart, tstop;
   double design;
   double deltaT;
   
   /* Get the time-steps */
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   braid_StepStatusGetTIndex(status, &ts_start);
   ts_start = GetTimeStepIndex(app, tstart); 
   ts_stop  = GetTimeStepIndex(app, tstop); 
   deltaT   = tstop - tstart;

   /* Get the current design from the shell part of the vector */
   design = u->design;

   /* Take one step forward */
   take_step(u->values, design, deltaT);

   /* Get new design from the app and store in shell part of the vector */
   int localindex = GetLocalDesignIndex(app, ts_stop);
   u->design = app->design[localindex];

   /* no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}   


/* Take one adjoint step backwards from tstop to tstart */
int 
my_Step_Adj(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     fstop,
            braid_Vector     u,
            braid_StepStatus status)
{
   int    ts_start, ts_stop;
   double tstart, tstop;
   double design;
   double deltaT;
   double localgrad = 0.0;
   
   /* Get the time-steps */
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   braid_StepStatusGetTIndex(status, &ts_start);
   ts_start = GetTimeStepIndex(app, tstart); 
   ts_stop  = GetTimeStepIndex(app, tstop); 
   deltaT   = tstop - tstart;

   /* Get the current design from the primal core */
   braid_BaseVector uprimal;
   int primaltimestep = GetPrimalTimeIndex(app, ts_stop);
   _braid_UGetVectorRef(app->primalcore, 0, primaltimestep, &uprimal);
   design = uprimal->userVector->design;

  /* Take one step backwards */
   localgrad = take_step_diff(u->values, deltaT, design);
   
   /* Add derivative of objective function */
   u->values[0] += 2 * uprimal->userVector->values[0];
   u->values[1] += 2 * uprimal->userVector->values[1];

   /* Store local gradient */
   u->gradient = localgrad;

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
   u->values     = (double*) malloc( 2*sizeof(double) );

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

   /* Initialize the design */
   int ts  = GetTimeStepIndex(app, t); 
   int idx = GetLocalDesignIndex(app, ts);
   u->design = app->design[idx];

   *u_ptr = u;
   return 0;
}


int
my_Init_Adj(braid_App     app,
            double        t,
            braid_Vector *u_ptr)
{

   int ts  = GetTimeStepIndex(app, t); 
   my_Vector *u;

   /* Allocate the vector */
   u = (my_Vector *) malloc(sizeof(my_Vector));
   u->values     = (double*) malloc( 2*sizeof(double) );

   /* Initialize the adjoint vector with derivative of objective function at last time step */
   if (t == 0)
   {
         int primaltimestep = GetPrimalTimeIndex(app, ts) ;
         braid_BaseVector uprimal;
         _braid_UGetVectorRef(app->primalcore, 0, primaltimestep, &uprimal);
         u->values[0] = 2 * uprimal->userVector->values[0];
         u->values[1] = 2 * uprimal->userVector->values[1];
   }
   else 
   {
         u->values[0] = 0.0;
         u->values[1] = 0.0;
   }

   /* Reset the rest */
   u->design   = 0.0;
   u->gradient = 0.0;

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
   v->values[0]  = u->values[0];
   v->values[1]  = u->values[1];
   v->design     = u->design;
   v->gradient   = u->gradient;

   *v_ptr = v;
   return 0;
}


int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u->values);
   u->values     = NULL;
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
      fprintf(file, "%1.14e, %1.14e, %1.14e\n", (u->values)[0], (u->values)[1], u->design);
      fflush(file);
      fclose(file);
   }


   return 0;
}

int
my_Access_Adj(braid_App          app,
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
      sprintf(filename, "%s.%04d.%03d", "ex-04.out.adjoint", index, app->myid);
      file = fopen(filename, "w");
      fprintf(file, "%1.14e, %1.14e, %1.14e\n", (u->values)[0], (u->values)[1], u->design);
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
   *size_ptr = 3*sizeof(double); // (design + two values)
   return 0;
}


int
my_BufPack(braid_App           app,
           braid_Vector        u,
           void               *buffer,
           braid_BufferStatus  bstatus)
{
   double *dbuffer = buffer;

   dbuffer[0] = u->values[0];
   dbuffer[1] = u->values[1];
   dbuffer[2] = u->design;

   braid_BufferStatusSetSize( bstatus,  3*sizeof(double));

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

   /* Allocate memory */
   u = (my_Vector *) malloc(sizeof(my_Vector));
   u->values     = (double*) malloc( 2*sizeof(double) );

   /* Unpack primal */
   u->values[0] = dbuffer[0];
   u->values[1] = dbuffer[1];
   u->design    = dbuffer[2];

   *u_ptr = u;
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
   int      rank, size, arg_index;
   int      ntime, maxiter;
   double   objective, stepsize, mygnorm, gnorm, gtol;
   int      max_levels, cfactor, access_level, print_level, braid_maxiter;
   double   braid_tol, braid_adjtol;

   /* Define time domain */
   ntime  = 20;              /* Total number of time-steps */
   tstart = 0.0;             /* Beginning of time domain */
   tstop  = 1.0;             /* End of time domain*/

   /* Define some optimization parameters */
   stepsize = 5.0;            /* Step size for design updates */
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
         printf("  min   \\int_0^1 u_1(t)^2 + u_2(t)^2 dt \n\n");
         printf("  s.t.  d/dt u_1(t) = u_2(t) \n");
         printf("        d/dt u_2(t) = -u_2(t) + c(t) \n\n");
         printf("  -ntime <ntime>          : set num points in time\n");
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

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   /* Set up the app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->myid     = rank;
   app->ntime    = ntime;
   app->tstop    = tstop;

   /* Initialize XBraid core for the forward pass */
   braid_Init(MPI_COMM_WORLD, MPI_COMM_WORLD, tstart, tstop, ntime, app, my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

   /* Get xbraid's time-grid distribution and store it in the app */
   int ilower, iupper;
   _braid_GetDistribution(core, &ilower, &iupper);
   app->ilower = ilower;

   /* Initialize local design vector and store it in the app */
   int ndesign = iupper - ilower + 1;
   double* design   = (double*) malloc( ndesign*sizeof(double) );    /* local design vector */
   for (int i = 0; i < ndesign; i++)
   {
        int ts = i + ilower;
        design[i] = ts * 0.1 + 1.0;
      //   printf("%d: design[%d,%d] = %1.14e\n", rank, i, ts, esign[i]);
   }
   app->design = design;

   printf("%d: Grid distribution: [%d, %d], total %d\n", rank, ilower, iupper, ndesign);

   /* Store all points (needed for adjoint) */
   braid_SetStorage(core, 0);

   /* Set some XBraid parameters */
   braid_SetMaxLevels(core, max_levels);
   braid_SetSkip(core, 0);
   braid_SetCFactor(core, -1, cfactor);
   braid_SetAccessLevel(core, access_level);
   braid_SetPrintLevel( core, print_level);       
   braid_SetMaxIter(core, braid_maxiter);
   braid_SetAbsTol(core, braid_tol);


   /* Store the primal core in the app so that the adjoint core can access primal values from xbraid */
   app->primalcore = core;

   /* Solve the forward ODE with XBraid */
   if (rank==0) printf("Solving ODE with XBraid...\n");
   braid_Drive(core);

   /* Evaluate the objective function */
   braid_BaseVector u;
   double obj;
   objective = 0.0;
   for (int n = 0; n <= ntime; n++)
   {
      /* Get braid vector at this time step */
      _braid_UGetVectorRef(core, 0, n, &u);
      /* Compute objective */
      if (u != NULL) // this is only true on one processor (the one that stores u)
      {
         obj = pow(u->userVector->values[0],2) + pow(u->userVector->values[1],2);
         objective += obj;
         // printf("%d: localobj %d, using %1.14e\n",rank, n, u->userVector->values[0]);
      }
   }
   MPI_Allreduce(&objective, &objective, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if (rank == 0) printf("Objective = %1.14e\n", objective);


   /* Set up XBraid core to solve the adjoint equation */
   braid_Core core_adj;
   braid_Init(MPI_COMM_WORLD, MPI_COMM_WORLD, tstart, tstop, ntime, app, my_Step_Adj, my_Init_Adj, my_Clone, my_Free, my_Sum, my_SpatialNorm, my_Access_Adj, my_BufSize, my_BufPack, my_BufUnpack, &core_adj);

   /* Set some XBraid parameters */
   braid_SetMaxLevels(core_adj, max_levels);
   braid_SetSkip(core_adj, 0);
   braid_SetCFactor(core_adj, -1, cfactor);
   braid_SetAccessLevel(core_adj, access_level);
   braid_SetPrintLevel( core_adj, print_level);       
   braid_SetMaxIter(core_adj, braid_maxiter);
   braid_SetAbsTol(core_adj, braid_tol);

   /* Tell XBraid to use reverted processor ranks when solving adjoint equation */
   braid_SetRevertedRanks(core_adj, 1);

   /* Store all points */
   braid_SetStorage(core_adj, 0);

   /* Solve the adjoint equation with XBraid */
   if (rank==0) printf("Solving adjoint equation with XBraid...\n");
   braid_Drive(core_adj);

   /* Compute gradient norm */
   gnorm = 0.0;
   for (int n = 0; n <= ntime; n++)
   {
      /* Get braid vector at this time step */
      _braid_UGetVectorRef(core_adj, 0, n, &u);
      /* Compute objective */
      if (u != NULL) // this is only true on one processor (the one that stores u)
      {
         double grad_ts = u->userVector->gradient;
         gnorm += grad_ts;
         // printf("%d: (%d) u->gradient = %1.14e\n", rank, ntime - n, grad_ts);
      }
   }
   mygnorm = gnorm;
   MPI_Allreduce(&mygnorm, &gnorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if (rank == 0) printf("Gradient norm: %1.14e\n", gnorm);


   /* Clean up */
   free(design);
   braid_Destroy(core);
   braid_Destroy(core_adj);
   free(app);
   
   MPI_Finalize();

   return (0);
}
