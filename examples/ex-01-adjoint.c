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

/**
 * Example:       ex-01.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-01
 *
 * Help with:     this is the simplest example available, read the source
 *
 * Sample run:    mpirun -np 2 ex-01
 *
 * Description:   solve the scalar ODE 
 *                   u' = lambda u, 
 *                   with lambda=-1 and y(0) = 1
 *                in a very simplified XBraid setting.
 *                
 *                When run with the default 10 time steps, the solution is:
 *                $ ./ex-01
 *                $ cat ex-01.out.00*
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
 * 
 *                 The time-averaged objective after convergence is 1.63614492461104e-01. 
 *                 The time-averaged objective after one iteration  1.77166086544457e-01. 
 * 
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.h"

#include "braid_test.h"

/*--------------------------------------------------------------------------
 * User-defined routines and structures
 *--------------------------------------------------------------------------*/

/* App structure can contain anything, and be named anything as well */
typedef struct _braid_App_struct
{
   int       rank;
   double    ntime;
   int       myid;
   double    design;          /* Store the design variables in the app */
   double    gradient;        /* Store the gradient in the app - should be of same size as design! */
   double    stepsize;        /* Stepsize for design updates */
   double    target;          /* Target for tracking type objective function (inverse design) */
   double    relax;           /* Relaxation parameter for the tracking type objective function */
} my_App;

/* Vector structure can contain anything, and be name anything as well */
typedef struct _braid_Vector_struct
{
   double value;
} my_Vector;

int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);

   /* Grab the design from the app */
   double lambda = app->design;

   /* Use backward Euler to propagate solution */
   (u->value) = 1./(1. - lambda* (tstop-tstart))*(u->value);
   
   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   if (t == 0.0) /* Initial condition */
   {
      (u->value) = 1.0;
   }
   else /* All other time points set to arbitrary value */
   {
      // (u->value) = 0.456;
      (u->value) = t+1;
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
   sprintf(filename, "%s.%04d.%03d", "ex-01.out", index, app->rank);
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
my_ObjectiveT(braid_App          app,
              braid_Vector       u,
              braid_AccessStatus astatus,
              double            *objectiveT_ptr)
{
   /* 1/N * f(u(t),lambda) = 1/N * u(t)**2 */
   double objT = 1. / (app->ntime) * (u->value) * (u->value);

   *objectiveT_ptr = objT;
   
   return 0;
}

int
my_PostprocessObjective(braid_App   app,
                        braid_Real  timeavg,
                        double     *postprocess
                        )
{
   double J;

   /* Tracking-type function */
   J  = 1./2. * pow(timeavg - app->target,2);
   /* Regularization term */
   J += + 1./2. * (app->relax) * pow(app->design,2);

   *postprocess = J;

   return 0;
}

int
my_PostprocessObjective_diff(braid_App   app,
                             braid_Real  timeavg,
                             double     *timeavg_bar
                             )
{
   double J_bar;

   /* Derivative of tracking type function */
   J_bar = timeavg - app->target;

   /* Derivative of regularization term */
   app->gradient = (app->relax) * (app->design);

   *timeavg_bar= J_bar;

   return 0;
}


int
my_Step_diff(braid_App              app,
                // braid_Vector     ustop,
                // braid_Vector     fstop,
                braid_Vector        u,
                braid_Vector        u_bar,
                braid_StepStatus    status)
{

   /* Get time that has been used in primal step evaluation  */
   double tstop, tstart, deltat;
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   deltat = tstop - tstart;

   double du, ddesign;

   /* Grab the design from the app */
   double lambda = app->design;

   /* Transposed derivative of step wrt u times u_bar */
   du = 1./(1. - lambda * deltat) * (u_bar->value);
 
   /* Transposed derivative of step wrt design times u_bar */
   ddesign = (deltat * (u->value)) / pow(1. - deltat*lambda,2) * (u_bar->value);

   /* Update u_bar and gradient */
   u_bar->value      = du;              // Make sure to do "=" here, not "+="! 
   app->gradient    += ddesign;


   return 0;
}


int
my_ObjectiveT_diff(braid_App         app,
                  braid_Vector       u,
                  braid_Vector       u_bar,
                  braid_Real         f_bar,
                  braid_Real         t, 
                  braid_Int          tidx)
{
   double du, ddesign; 

   /* Partial derivative with respect to u times f_bar */
   du = 2. / (app->ntime) * u->value * f_bar;

   /* Partial derivative with respect to design times f_bar*/
   ddesign = 0.0 * f_bar;

   /* Update u_bar and gradient */
   u_bar->value  += du;
   app->gradient += ddesign;

   return 0;
}

int 
my_AccessGradient(braid_App app)
{
   /* Print the gradient */
   if (app->myid==0)
   {
     printf("Gradient: %1.14e\n", app->gradient);
   } 

   return 0;
}

int
my_AllreduceGradient(braid_App app, 
                     MPI_Comm comm)
{
   double localgradient;
   double globalgradient;

   /* Collect sensitivities from all time-processors */
   localgradient = app->gradient;
   MPI_Allreduce(&localgradient, &globalgradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   /* Broadcast the global gradient to all processor */
   app->gradient = globalgradient; 

   return 0;
}                  


int 
my_ResetGradient(braid_App app)
{
   /* Set the gradient to zero */
   app->gradient = 0.0;

   return 0;
}

int
my_GradientNorm(braid_App app,
                double   *gradient_norm_prt)

{
   double gnorm;

   /* Norm of gradient */
   gnorm = sqrt( (app->gradient)*(app->gradient) );

   /* Return the gradient norm */
   *gradient_norm_prt = gnorm;

   return 0;
}
             

int 
my_UpdateDesign(braid_App app, 
                double    objective,
                double    rnorm,
                double    rnorm_adj)
{
   double Hinv;

   /* Hessian approximation */
   Hinv = 1. ;

   /* Design update */
   app->design = app->design - (app->stepsize) * Hinv * (app->gradient);
   printf("Design: %1.14e\n", app->design);

   return 0;
}             


/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
   braid_Core    core;
   my_App       *app;
   double        tstart, tstop;
   int           ntime, rank;
   double        design; 
   double        gradient;
   double        target;
   double        relax;
   double        stepsize;

   /* Define time domain: ntime intervals */
   ntime  = 50;
   tstart = 0.0;
   tstop  = tstart + ntime/2.;

   /* Initialize optimization */
   design           = -1.0;                 // Initial design 
   gradient         = 0.0;                  // Initial gradient
   target           = 1.15231184218078e-01; // Inverse Design Target: precomputed from design = -0.2, ntime=50
   relax            = 0.0005;               // Relaxation parameter
   stepsize         = 6.0;                  // Stepsize for design updates 

   /* DEBUG gradient */
   // design += 1e-8;
   
   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->rank             = rank;
   app->design           = design;
   app->gradient         = gradient;
   app->ntime            = ntime;
   app->myid             = rank;
   app->stepsize         = stepsize;
   app->target           = target;
   app->relax            = relax;

   /* initialize XBraid */
   braid_Init(MPI_COMM_WORLD, MPI_COMM_WORLD, tstart, tstop, ntime, app,
             my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
             my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);


   /* Initialize adjoint XBraid */
   braid_Init_Adjoint( my_ObjectiveT, my_Step_diff, my_ObjectiveT_diff, my_AllreduceGradient, my_ResetGradient, my_AccessGradient, my_GradientNorm, my_UpdateDesign, &core);
   
   /* Set some typical Braid parameters */
   braid_SetPrintLevel( core, 1);
   braid_SetMaxLevels(core, 2);
   braid_SetAbsTol(core, 1.0e-06);
   braid_SetCFactor(core, -1, 2);
   braid_SetAccessLevel(core, 1);
   braid_SetVerbosity(core, 0);
   braid_SetMaxIter(core, 10);
   braid_SetMaxOptimIter(core, 100);
   // braid_SetGradientAccessLevel(core, 1);   
   // braid_SetTolAdjoint(core, 1e-9);        
   // braid_SetTolGradient(core, 1e-9);
   // braid_SetTolDesignUpdate(core, 1e-8);

   /* debug: don't skip work on downcycle for comparison with adjoint run.*/
   braid_SetSkip(core, 0);

   /* Optional: Set the time for starting time-average */
   // braid_SetTStartTimeaverage( core, 1.0);
   // braid_SetTStopTimeaverage( core, tstop);

   /* Optional: Set the tracking type objective function and derivative */
   braid_SetPostprocessObjective(core, my_PostprocessObjective);
   braid_SetPostprocessObjective_diff(core, my_PostprocessObjective_diff);


   /* Run simulation */
   braid_Drive(core);


   /* Clean up */
   free(app);
   braid_Destroy(core);
   MPI_Finalize();

   return (0);
}
