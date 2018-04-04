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
 * Example:       ex-04-adjoint.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-04-adjoint
 *
 * Description:   TODO: Describe
 *
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
   int     myid;      /* Rank of the processor */
   double *design;    /* Holds the time-dependent design vector */
   double *gradient;  /* Holds the gradient vector */
   double  gamma;     /* Relaxation parameter for objective function */
   int     ntime;     /* number of time-steps */

} my_App;


/* Define the state vector at one time-step */
typedef struct _braid_Vector_struct
{
   double *values;     /* Holds an R^2 vector */

} my_Vector;


/*--------------------------------------------------------------------------
 * Integration routines
 *--------------------------------------------------------------------------*/

int 
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   int    index;
   double tstart, tstop;
   double design;
   double deltaT;
   
   /* Get the time-step size */
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   deltaT = tstop - tstart;

   /* Get the current design */
   braid_StepStatusGetTIndex(status, &index);
   design = app->design[index];

   /* Take one step */
   take_step(u->values, design, deltaT);

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
   
   for (int i = 0; i < 2; i++)
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

   /* Allocate memory */
   u = (my_Vector *) malloc(sizeof(my_Vector));
   u->values = (double*) malloc( 2*sizeof(double) );

   /* Unpack the buffer */
   for (int i = 0; i < 2; i++)
   {
      (u->values)[i] = dbuffer[i];
   }

   *u_ptr = u;
   return 0;
}

int 
my_ObjectiveT(braid_App          app,
              braid_Vector       u,
              braid_AccessStatus astatus,
              double            *objectiveT_ptr)
{
   double objT;
   double design;
   int    index;
   double deltaT = 1./app->ntime;

   /* Get the time index*/
   braid_AccessStatusGetTIndex(astatus, &index);

   /* Evaluate the objective function after first step */
   if ( index > 0)
   {
      /* Get the current design */
      design = app->design[index-1];

      /* Evaluate objective */
      objT = evalObjectiveT( u->values, design, deltaT, app->gamma);
   }

   *objectiveT_ptr = objT;
   
   return 0;
}


int
my_ObjectiveT_diff(braid_App         app,
                  braid_Vector       u,
                  braid_Vector       u_bar,
                  braid_Real         f_bar,
                  braid_Real         t, 
                  braid_Int          index)
{
   double deltaT = 1. / app->ntime;
   double  ddesign; 
   double *du = (double*) malloc(2*sizeof(double));

   /* Compute derivative if time index > 0 */
   if ( index > 0 )
   {
      /* Transposed derivative of objectiveT wrt u times f_bar */
      du[0] = -2. * deltaT * u->values[0] * f_bar;
      du[1] = -2. * deltaT * u->values[1] * f_bar;

      /* Transposed derivative of objective wrt design times f_bar */
      ddesign = 2. * deltaT * app->gamma * app->design[index-1] * f_bar;

      /* Update u_bar and gradient */
      u_bar->values[0]       += du[0];
      u_bar->values[1]       += du[1];
      app->gradient[index-1] += ddesign;
   }
   
   free(du);

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

   double tstop, tstart, deltaT;
   int tidx;
   double *du;
   double ddesign;

   /* Get time and index that have been */
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   braid_StepStatusGetTIndex(status, &tidx);
   deltaT = tstop - tstart;

   /* Transposed derivative of step wrt u times u_bar */
   du = (double*) malloc(2*sizeof(double));
   du[0] = u_bar->values[0];
   du[1] = u_bar->values[1] + deltaT * u_bar->values[0] - deltaT * u_bar->values[1];

   /* Transposed derivative of step wrt design times u_bar */
   ddesign = - deltaT * u_bar->values[1];

   /* Update u_bar and gradient */
   u_bar->values[0]     = du[0];
   u_bar->values[1]     = du[1];
   app->gradient[tidx] += ddesign;

   // printf("%d: mystep_diff %d: %1.14e\n", app->myid, tidx, ddesign);

   free(du);
   return 0;
}


int
my_AllreduceGradient(braid_App app, MPI_Comm comm)
{
   double localgradient;
   double globalgradient;

   /* Collect sensitivities from all time-processors */
   for (int ts = 0; ts < app->ntime; ts++)
   {
       localgradient = app->gradient[ts];
      MPI_Allreduce(&localgradient, &globalgradient, 1, MPI_DOUBLE, MPI_SUM, comm);
      /* Broadcast the global gradient to all processor */
      app->gradient[ts] = globalgradient;
   }

   return 0;
}


int 
my_ResetGradient(braid_App app)
{
   int ts;
   /* Set the gradient to zero */
   for(ts = 0; ts < app->ntime; ts++) 
   {
      app->gradient[ts] = 0.0;
   }

   return 0;
}

int 
my_UpdateDesign(braid_App app, 
                double    objective,
                double    rnorm,
                double    rnorm_adj,
                double   *gradient_norm_ptr)
{
   
   *gradient_norm_ptr = 0.0;
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
   int      ntime, ts;
   double   gamma;
   double  *design; 
   double  *gradient; 
   int      rank;

   /* Define time domain */
   ntime  = 20;
   tstart = 0.0;
   tstop  = 1.0;

   /* Define some optimization parameters */
   gamma = 0.005;          /* Relaxation parameter */

   /* Initialize optimization */
   design   = (double*) malloc( ntime*sizeof(double) ); /* design vector */
   gradient = (double*) malloc( ntime*sizeof(double) ); /* gradient vector */

   /* Initialize the optimization variables */
   for (ts = 0; ts < ntime; ts++)
   {
      design[ts]   = 0.;
      gradient[ts] = 0.;
   }

   /* DEBUG gradient */
   // design[9] += 1e-8;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);


   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->myid     = rank;
   app->ntime    = ntime;
   app->design   = design;
   app->gradient = gradient;
   app->gamma    = gamma;

   /* Initialize XBraid */
   braid_Init(MPI_COMM_WORLD, MPI_COMM_WORLD, tstart, tstop, ntime, app, my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

   /* Initialize adjoint XBraid */
   braid_Init_Adjoint( my_ObjectiveT, my_Step_diff, my_ObjectiveT_diff, my_AllreduceGradient, my_ResetGradient, my_UpdateDesign, &core);

   /* Set some Braid parameters */
   braid_SetPrintLevel( core, 1);
   braid_SetMaxLevels(core, 2);
   braid_SetAbsTol(core, 1.0e-06);
   braid_SetCFactor(core, -1, 2);
   braid_SetAccessLevel(core, 1);
   braid_SetMaxIter(core, 10);
   // braid_SetGradientAccessLevel(core, 2);

   // /* debug: never skip work on downcycle for comparing primal and adjoint run.*/
   braid_SetSkip(core, 0);


   /* Run a Braid simulation */
   braid_Drive(core);


   /* Clean up */
   free(design);
   free(gradient);
   free(app);

   braid_Destroy(core);
   MPI_Finalize();

   return (0);
}
