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
 * Example:       ex-04-serial.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-04-serial
 *
 * 
 * Description:  Solves a simple optimal control problem in time-serial:
 * 
 *                 min   \int_0^1 u_1(t)^2 + u_2(t)^2 + gamma c(t)^2  dt
 * 
 *                  s.t.  d/dt u_1(t) = u_2(t)
 *                        d/dt u_2(t) = -u_2(t) + c(t)
 * 
 *                 with initial condition u_1(0) = 0, u_2(0) = -1
 *                 and piecewise constant control c(t).  
 * 
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ex-04-lib.c"



int main (int argc, char *argv[])
{
   double   tstart, tstop, deltaT;
   int      ntime, ts;
   int      maxiter,iter;
   double   objective;
   double   gamma;
   double  *u, *w;
   double  *u0, *u1;
   double  *design; 
   double  *gradient;
   double   gnorm;
   double   gtol;
   double   stepsize;
   
   /* Define time domain */
   ntime  = 20;                   /* Total number of time-steps */
   tstart = 0.0;                  /* Beginning of time-domain */
   tstop  = 1.0;                  /* End of time-domain */
   deltaT = ( tstop - tstart ) / ntime;  /* Time-step size */

   /* Define some optimization parameters */
   gamma    = 0.005;          /* Relaxation parameter in the objective function */
   stepsize = 50.0;           /* Step size for design updates */
   maxiter  = 300;            /* Maximum number of optimization iterations */
   gtol     = 1e-6;           /* Stopping criterion on the gradient norm */
   

   /* Initialize the optimization variables */
   u        = (double*) malloc( 2*sizeof(double) );      /* State at on time-step */
   w        = (double*) malloc( 2*sizeof(double) );      /* Adjoint at on time-step */
   u0       = (double*) malloc( ntime*sizeof(double) );  /* Stores the space-time state, 1st component */
   u1       = (double*) malloc( ntime*sizeof(double) );  /* Stores the space-time state, 2nd component */
   design   = (double*) malloc( ntime*sizeof(double) );  /* design vector (control c) */
   gradient = (double*) malloc( ntime*sizeof(double) );  /* gradient vector */

   for (ts = 0; ts < ntime; ts++)
   {
      u0[ts]       = 0.;
      u1[ts]       = 0.;
      design[ts]   = 0.;
      gradient[ts] = 0.;
   }

   /* Prepare optimization output */
   printf("\n#    Objective             || Gradient ||\n");

   /* Optimization iteration */
   for (iter = 0; iter < maxiter; iter++)
   {

      /* Set initial condition for the state */
      u[0] =  0.0;
      u[1] = -1.0;

      /* Main forward time-stepping loop */
      objective = 0.0;
      for (ts = 1; ts <= ntime; ts++)
      {
         /* Take the step */
         take_step(u, design[ts-1], deltaT );

         /* Add to the objective function */
         objective += evalObjectiveT(u, design[ts-1], deltaT, gamma);

         /* Store the state variable */
         u0[ts-1] = u[0];
         u1[ts-1] = u[1];
      }

      /* Set terminal condition for the adjoint */
      w[0] = 0.0;
      w[1] = 0.0;

      /* Main backward adjoint time-stepping loop */
      for (ts = ntime; ts >= 1; ts--)
      {
         /* Take adjoint step */
         take_adjoint_step(w, u0[ts-1], u1[ts-1], deltaT);

         /* Evaluate the gradient */
         gradient[ts-1] = evalGradientT(w, design[ts-1], deltaT, gamma);
      }

      /* Compute norm of gradient */
      gnorm = compute_sqnorm(gradient, ntime);
      gnorm = sqrt(gnorm);


      /* Output */
      printf("%3d  %1.14e  %1.14e\n", iter, objective, gnorm);


      /* Check optimization convergence */
      if (gnorm < gtol)
      {
         break;
      }

      /* Design update */
      for(ts = 0; ts < ntime; ts++) 
      {
         design[ts] -= stepsize * gradient[ts];
      }
   }

   if (iter == maxiter)
   {
      printf("\n Max. number of iterations reached.\n\n");
   }
   else
   {
      printf("\n Optimization has converged.\n\n");
   }

   /* Write final state to file */
   write_vec("state", u0, u1, ntime);

   /* Write final design to file */
   FILE *designfile;
   designfile = fopen("ex-04.out.design","w");
   for (ts = 0; ts < ntime; ts++)
   {
      fprintf(designfile, "%1.14e\n", design[ts]);
   }
   fflush(designfile);
   fclose(designfile);

   /* Free memory */
   free(design);
   free(u0);
   free(u1);
   free(u);
   free(gradient);
   free(w);

   return (0);
}
