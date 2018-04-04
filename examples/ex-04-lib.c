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
 * This file contains library functions for ex-04-adjoint.c and ex-04-adjoint.serial.c. 
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

/* Advance a state vector u forward in time */
void 
take_step(double* u,        /* state at current time */
          double design,    /* design at current time */
          double deltaT )   /* time step size */
{
   double du0, du1;

   /* Compute the update */
   du0 =  deltaT * u[1];
   du1 =  deltaT * ( -u[1] + design);
   
   /* Update the state */
   u[0] = u[0] +  du0;
   u[1] = u[1] +  du1;

}

/* Evaluate the local objective function at the current time */
double 
evalObjectiveT(double* u,       /* state at current time */
               double design,   /* design at current time */
               double deltaT,   /* time step size */
               double gamma)    /* relaxation parameter */
{
   double objectiveT;

   objectiveT  =  deltaT * ( u[0]*u[0] + u[1]*u[1] )
                + deltaT * gamma * ( design*design );

   return objectiveT;
}


/*--------------------------------------------------------------------------
 * Routines for solving the adjoint equation 
 *--------------------------------------------------------------------------*/

/* Advance an adjoint variable backwards in time */
void 
take_adjoint_step(double *w,      /* adjoint at current time */
                  double u0,      /* 1st component of state */
                  double u1,      /* 2nd omponent of state */
                  double deltaT ) /* time step size */
{
   double dw0, dw1;
   
   /* Transposed derivative of step wrt u times w */
   dw0  = w[0];
   dw1  = w[1] + deltaT * w[0] - deltaT * w[1];

   /* Transposed derivative of the objective wrt u */
   dw0 -= 2. * deltaT * u0;
   dw1 -= 2. * deltaT * u1;

   /* Update */
   w[0] = dw0;
   w[1] = dw1;
}


/* Evaluate the local gradient at the current time */
double 
evalGradientT(double* w,       /* adjoint at current time */
              double  design,  /* design  at current time */
              double  deltaT,  /* time-step size */
              double  gamma )  /* relaxation parameter */
{
   double gradient;

   /* Derivative of the step wrt design times w*/
   gradient = - deltaT * w[1];

   /* Derivative of objective wrt design */
   gradient += 2.* deltaT * gamma * design;

   return gradient;
}


/*--------------------------------------------------------------------------
 * Utility routines
 *--------------------------------------------------------------------------*/

/* Write state or adjoint vector to file */
void 
write_vec(char*   varname,  /* Name of the vector that is to be written */
      double* var0,     /* first components of the vector */
      double* var1,     /* second components of the vector */
      int     ntime )   /* total number of time steps */
{
   char  filename[255];
   FILE *file;

   /* Open file for output  */
   sprintf(filename, "ex-04.out.%s", varname);
   file = fopen(filename, "w");
   for (int ts = 0; ts < ntime; ts++)
   {
      fprintf(file, "%1.14e  %1.14e\n", var0[ts], var1[ts]);
   }
   fflush(file);
   fclose(file);
}              
