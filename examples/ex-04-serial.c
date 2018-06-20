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
 *               Implements a steepest-descent optimization iteration
 *               using fixed step size for design updates.   
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ex-04-lib.c"



int main (int argc, char *argv[])
{
   double   tstart, tstop, deltaT;
   int      ntime, ts, arg_index;
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
   ntime  = 20;                         /* Total number of time-steps */
   tstart = 0.0;                        /* Beginning of time-domain */
   tstop  = 1.0;                        /* End of time-domain */
   deltaT = (tstop - tstart) / ntime;   /* Time-step size */

   /* Define optimization parameters */
   gamma    = 0.005;                    /* Relaxation parameter in the objective function */
   stepsize = 50.0;                     /* Step size for design updates */
   maxiter  = 500;                      /* Maximum number of optimization iterations */
   gtol     = 1e-6;                     /* Stopping criterion on the gradient norm */
   
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
         printf("  -ntime <ntime>       : set num points in time\n");
         printf("  -gamma <gamma>       : Relaxation parameter in the objective function \n");
         printf("  -stepsize <stepsize> : Step size for design updates \n");
         printf("  -maxiter <maxiter>   : Maximum number of optimization iterations \n");
         printf("  -gtol <gtol>         : Stopping criterion on the gradient norm \n");
         exit(1);
      }
      else if ( strcmp(argv[arg_index], "-ntime") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
         deltaT = (tstop - tstart) / ntime;  /* recompute */
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
      else if ( strcmp(argv[arg_index], "-maxiter") == 0 )
      {
         arg_index++;
         maxiter = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-gtol") == 0 )
      {
         arg_index++;
         gtol = atof(argv[arg_index++]);
      }
      else
      {
         printf("ABORTING: incorrect command line parameter %s\n", argv[arg_index]);
         return (0);
      }
   }

   /* Initialize the optimization variables */
   u        = (double*) malloc( 2*sizeof(double) );      /* State at a time-step */
   w        = (double*) malloc( 2*sizeof(double) );      /* Adjoint at a time-step */
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
         /* Restore the state variable */
         u[0] = u0[ts-1];
         u[1] = u1[ts-1];

         /* Take adjoint step and evaluate gradient */
         gradient[ts-1] = take_adjoint_step(w, u, design[ts-1], gamma, deltaT);
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

   /* Write final state and design to file */
   write_vec("state", u0, u1, ntime);
   write_design_vec("design", design, ntime);

   /* Free memory */
   free(design);
   free(u0);
   free(u1);
   free(u);
   free(gradient);
   free(w);

   return (0);
}
