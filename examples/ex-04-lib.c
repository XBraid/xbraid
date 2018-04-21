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
 *
 * This file contains library functions for ex-04.c and ex-04.serial.c. 
 * Together, these files solve a simple optimal control problem. 
 * 
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/*--------------------------------------------------------------------------
 * Routines for solving the ODE constraint and evaluating the objective
 *--------------------------------------------------------------------------*/


/* Advance a state vector u forward in time 
 * Implements simple explicite Euler integration scheme */
void 
take_step(double* u,         /* state at current time */
          double  design,    /* design at current time */
          double  deltaT )   /* time step size */
{
   double du0, du1;

   /* Compute the update */
   du0 =  deltaT * u[1];
   du1 =  deltaT * ( -u[1] + design);
   
   /* Update the state */
   u[0] = u[0] +  du0;
   u[1] = u[1] +  du1;

}

/* Evaluate the time-dependent objective function at the current time */
double 
evalObjectiveT(double* u,        /* state at current time */
               double  design,   /* design at current time */
               double  deltaT,   /* time step size */
               double  gamma)    /* relaxation parameter */
{
   double objectiveT;

   objectiveT  =  deltaT * ( u[0]*u[0] + u[1]*u[1] )
                + deltaT * gamma * ( design*design );

   return objectiveT;
}


/*--------------------------------------------------------------------------
 * Routines for solving the adjoint equation 
 *--------------------------------------------------------------------------*/

/* Compute the transposed partial derivative of take_step multiplied with current adjoint w */
double
take_step_diff(double *w, 
               double  deltaT)
{
   double *ddu = (double*)malloc(2*sizeof(double));
   double gradientT;

   /* derivative with respect to u  */
   ddu[0] = w[0];
   ddu[1] = w[1] + deltaT * w[0] - deltaT * w[1];

   /* derivative with respect to c  */
   gradientT = -deltaT * w[1];

   /* Update adjoint */
   w[0] = ddu[0];
   w[1] = ddu[1];

   free(ddu);

   return gradientT;
}            

/* Partial derivative of the objective */
double
evalObjectiveT_diff(double *w, 
                    double *u, 
                    double  design,
                    double  gamma,
                    double  deltaT)
{
   double gradientT;

   /* derivative with respect to u */
   w[0] -= 2. * deltaT * u[0];
   w[1] -= 2. * deltaT * u[1];

   /* derivative with respect to c */
   gradientT = 2.* deltaT * gamma * design;
   
   return gradientT;
}

/* Advance an adjoint variable backwards in time 
 * and evaluate local gradient at that time step */
double
take_adjoint_step(double *w,         /* adjoint variable that gets propagated backwards */
                  double *u,         /* state variable at the current time */
                  double  design,    /* design variable at the current time */
                  double  gamma,     /* relaxation parameter in the objective function */
                  double  deltaT)    /* time step size */
{
   double gradientT = 0.0;

   /* transposed derivatives of take_step times w */
   gradientT += take_step_diff(w, deltaT);

   /* transposed derivatives of evalObjectiveT */
   gradientT += evalObjectiveT_diff(w, u, design, gamma, deltaT);

   return gradientT;
}                 


/*--------------------------------------------------------------------------
 * Utility routines
 *--------------------------------------------------------------------------*/

/* Write state or adjoint vector to file */
void 
write_vec(char*   name,     /* Filename extension (ex-04.out.name) */
          double* var0,     /* first components of the vector */
          double* var1,     /* second components of the vector */
          int     ntime )   /* total number of time steps */
{
   char  filename[255];
   FILE *file;

   /* Open file for output  */
   sprintf(filename, "ex-04.out.%s", name);
   file = fopen(filename, "w");
   for (int ts = 0; ts < ntime; ts++)
   {
      fprintf(file, "%1.14e  %1.14e\n", var0[ts], var1[ts]);
   }
   fflush(file);
   fclose(file);
}              

/* Compute the squared norm of a vector */
double 
compute_sqnorm(double *vector, 
               int     size)

{
   double norm = 0.0;

   for(int i = 0; i < size; i++) 
   {
      norm += (vector[i]) * (vector[i]) ;
   }

   return norm;
}

