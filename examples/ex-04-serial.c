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
 * Example:       ex-04-adjoint-serial.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-04-adjoint-serial
 *
 * Description:   TODO: Describtion
 * 
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ex-04-lib.c"



int main (int argc, char *argv[])
{
   double   tstart, tstop, deltaT;
   int      ntime;
   int      ts;
   double  *u, *w;
   double  *design; 
   double  *u0, *u1;
   double  *w0, *w1;
   double   objective;
   double   gamma;
   double  *gradient;
   
   /* Define time domain */
   ntime  = 20;
   tstart = 0.0;
   tstop  = 1.0;
   deltaT = ( tstop - tstart ) / ntime;

   /* Define some optimization parameters */
   gamma = 0.005;          /* Relaxation parameter */

   /* Allocate memory */
   u0 = (double*) malloc( ntime*sizeof(double) );   /* state, 1st component */
   u1 = (double*) malloc( ntime*sizeof(double) );   /* state, 2nd component */
   w0 = (double*) malloc( ntime*sizeof(double) );   /* adjoint, 1st component */
   w1 = (double*) malloc( ntime*sizeof(double) );   /* adjoint, 2nd component */
   design   = (double*) malloc( ntime*sizeof(double) ); /* design vector */
   gradient = (double*) malloc( ntime*sizeof(double) ); /* gradient vector */
   u  = (double*) malloc( 2*sizeof(double) );       /* Temporary state */
   w  = (double*) malloc( 2*sizeof(double) );       /* Temporary adjoint */

   /* Initialize the optimization variables */
   for (ts = 0; ts < ntime; ts++)
   {
      u0[ts] = 0.;
      u1[ts] = 0.;
      w0[ts] = 0.;
      w1[ts] = 0.;
      design[ts]  = 0.;
      gradient[ts] = 0.;
   }

   /* Set the state initial condition */
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

   /* Output */
   write_vec("state", u0, u1, ntime);
   printf("\n Objective: %1.14e\n\n", objective);

   /* Set the adjoint terminal condition */
   w[0] = 0.0;
   w[1] = 0.0;

   /* Main backward adjoint time-stepping loop */
   for (ts = ntime; ts >= 1; ts--)
   {
      /* Take adjoint step */
      take_adjoint_step(w, u0[ts-1], u1[ts-1], deltaT);

      /* Evaluate the gradient */
      gradient[ts-1] = evalGradientT(w, design[ts-1], deltaT, gamma);

      /* Store the adjoint variable */
      w0[ts-1] = w[0];
      w1[ts-1] = w[1];
   }

   /* Output */
   write_vec("adjoint", w0, w1, ntime);
   for (ts = 0; ts < ntime; ts++)
   {
      printf(" Gradient %d %1.14e\n", ts, gradient[ts]);
   }



   /* Free memory */
   free(design);
   free(u0);
   free(u1);
   free(u);
   free(gradient);
   free(w0);
   free(w1);
   free(w);

   return (0);
}
