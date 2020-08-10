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
 *  This file contains library functions for ex-02.c and ex-02-serial.c. 
 *  Together, these files solve the 1D heat equation. 
 *                                   
 *  For more details on the discretization, see the header comment in ex-02.c.
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifdef M_PI
   #define PI M_PI
#else
   #define PI 3.14159265358979
#endif

/*--------------------------------------------------------------------------
 * Routines defining solution, forcing term and how to compute the error 
 *--------------------------------------------------------------------------*/

/* Exact solution */ 
double exact(double t, double x)
{
    return sin(x)*cos(t);
}

/* Initialize array of values to solution at time t */
void      
get_solution(double * values,
             int      size,
             double   t, 
             double   xstart, 
             double   deltaX)
{
   int i;

   for(i=0; i < size; i++)
   {
      values[i] = exact(t, xstart);
      xstart += deltaX;
   }
}

/* Forcing term F(t, x) for PDE,  u_t = u_xx + F(t,x) */
double forcing(double t, double x)
{
   return (-1.0)*sin(x)*sin(t) + sin(x)*cos(t);
}

/* Compute L2-norm of the error at a point in time */
double compute_error_norm(double * values, 
                          double   xstart, 
                          double   xstop, 
                          int      nspace, 
                          double   t)
{
   int i;
   double deltaX = (xstop - xstart) / (nspace - 1.0);
   double x = xstart;
   double error = 0.0;

   for(i = 0; i < nspace; i++)
   {
      error = error + pow( values[i] - exact(t, x), 2);
      x += deltaX;
   }

   return sqrt( error*deltaX ); 
}


/*--------------------------------------------------------------------------
 * Routines for computing time-steps and residuals 
 *--------------------------------------------------------------------------*/

/* Helper function for tridiagonal solver */
double dabs(double x)
{
   if (x < 0.0)
   {
      return -1.0;
   }
   else if (x > 0.0)
   {
      return 1.0;
   }
   else
   {
      return 0.0;
   }
}


/* Helper function for Step: Tridiagonal system solver (Thomas algorithm) */
void 
solve_tridiag(double *x, double *g, int N, double* matrix)
{
   /**
    * solves Ax = v where A is a tridiagonal matrix with stencil [ a, b, c].
    * 
    * There is a built in assumption that the first and last rows are the 
    * identity (boundary conditions)
    *
    * Input
    * -----
    * x - initially contains v
    * g - temp data array for the algorithm
    * N - length of vectors x and g
    * matrix - length three array representing the tridiagonal stencil
    * 
    * Output
    * ------
    * x - contains solution upon output (over-written)
    * g - is a working array, and will be modified as such
    **/

   int i;
   double m;
       
   g[0] = 0.0;  /* Assume the first row is the identity */
   
   /* loop from 1 to N - 2 inclusive, performing the forward sweep */
   for (i = 1; i < N - 1; i++) 
   {
       m = 1.0 / (matrix[1] - matrix[0]*g[i-1]);
       g[i] = m*matrix[2];
       x[i] = m*(x[i] - matrix[0]*x[i-1]);
   }
   
   /* Do nothing for x[N-1], assume last row is the identity */
   
   /* loop from N - 2 to 1 inclusive to perform the back substitution */
   for (i = N - 2; i >= 1; i--)
   {
       x[i] = x[i] - g[i] * x[i+1];
   }
}

/* Helper function for Residual: Carry out mat-vec with tridiagonal stencil */
void 
matvec_tridiag(double *x, double *g, int N, double* matrix)
{
   /**
    * Matvec solves g <-- Ax, where A is a tridiagonal matrix with stencil [ a, b, c].
    * 
    * There is a built in assumption that the first and last rows are the 
    * identity (boundary conditions)
    *
    * Input
    * -----
    * x - input vector 
    * N - length of vectors x and g
    * matrix - length three array representing the tridiagonal stencil
    * 
    * Output
    * ------
    * g - Equals A*x
    **/

   int i;
       
   /* loop from 1 to N - 2 inclusive, performing the matvec */
   for (i = 1; i < N - 1; i++) 
   {
       g[i] = matrix[0]*x[i-1] + matrix[1]*x[i] + matrix[2]*x[i+1];
   }
   
   /* boundary points */
   g[0] = x[0];
   g[N-1] = x[N-1];

}

/* Compute three point backward Euler stencil */
void
compute_stencil(double   deltaX, 
                double   deltaT,
                double * matrix)
{
   double cfl;

   cfl = (deltaT/(deltaX*deltaX));
   matrix[0] = -cfl;
   matrix[1] = 1.0 + 2*cfl;
   matrix[2] = -cfl;
}

/* Take a backward Euler step */
void
take_step(double * values,    /* state vector to evolve */
          int      size,      /* size of vector */
          double   t,         /* time value to evolve to */
          double   xstart,    /* left-most spatial grid point */
          double   deltaX,    /* spatial mesh size */
          double   deltaT,    /* temporal mesh size */
          double * matrix,    /* three-point matrix stencil */
          double * temp)      /* temporary array of size 'size' */
{
   int i; 
   double x;
   
   /* Set up matrix stencil for 1D heat equation*/
   compute_stencil(deltaX, deltaT, matrix);

   /* Apply boundary conditions */
   values[0] = exact(t, 0.0);
   values[size-1] = exact(t, PI);

   /* PDE forcing */
   x = xstart;
   for(i = 0; i < size; i++)
   {
      values[i] = values[i] + deltaT*forcing(t, x);
      x = x + deltaX;
   }

   /* Backward Euler step */
   solve_tridiag(values, temp, size, matrix);
}



/*--------------------------------------------------------------------------
 * Routines for coarsening and interpolating in space 
 *--------------------------------------------------------------------------*/

/* Bilinear interpolation from grid size 2^{k-1} + 1 
 * to grid size 2^{k} + 1 */
void
interpolate_1D(double * cvalues,
               double * fvalues, 
               int      csize,
               int      fsize)
{
   int i;
   
   for (i = 1; i < fsize-1; i++)
   {
      if(i%2 == 1)
         (fvalues)[i] = 0.5*cvalues[i/2] + 0.5*cvalues[(i+1)/2];
      else
         (fvalues)[i] = cvalues[i/2];
   }
   
   /* Boundary Conditions */
   (fvalues)[0] = cvalues[0];
   (fvalues)[fsize-1] = cvalues[csize-1];
}

/* Bilinear coarsening from grid size 2^{k1} + 1 
 * to grid size 2^{k-1} + 1 */
void
coarsen_1D(double * cvalues,
           double * fvalues, 
           int      csize,
           int      fsize)
{
   int i, fidx;
   
   for (i = 1; i < csize-1; i++)
   {
      fidx = 2*i;
      cvalues[i] = 0.5*fvalues[fidx] + 0.25*fvalues[fidx+1] + 0.25*fvalues[fidx-1];
   }
   
   /* Boundary Conditions */
   cvalues[0] = fvalues[0];
   cvalues[csize-1] = fvalues[fsize-1];
}


/*--------------------------------------------------------------------------
 * Helper routines for output
 *--------------------------------------------------------------------------*/

/* Print accumulated info on space-time grids visited during the simulation */
void
print_sc_info(double *sc_info, 
              int     max_levels)
{
   int i;
   double dx, dt; 

   printf("\n-----------------------------------------------------------------\n"); 
   printf("-----------------------------------------------------------------\n\n"); 
   printf( " Per level diagnostic information \n\n");
   printf("level       dx          dt        dt/dx^2\n"); 
   printf("-----------------------------------------------------------------\n"); 
   for( i = 0; i < max_levels; i++)
   {
      dx = sc_info[i*2];
      dt = sc_info[i*2+1];
      if (dx == -1){
         break;
      }
      printf(" %2d   |   %1.2e    %1.2e    %1.2e\n", i, dx, dt, dt/(dx*dx) ); 
   }
   printf( "\n" );
}

/* Print a solution vector to file with the format
 *
 *    ntime_steps
 *    tstart
 *    tstop
 *    nspace_points
 *    xstart
 *    xstop
 *    x[0]
 *    x[1]
 *      .
 *      .
 *      .
 *    x[k]
 *
 **/
void
save_solution(char   *filename, 
              double *values, 
              int     size, 
              double  xstart,
              double  xstop,
              int     ntime,
              double  tstart,
              double  tstop,
              double  current_t)
{
   FILE      *file;
   int i;
   
   file = fopen(filename, "w");
   fprintf(file, "%d\n",    ntime +1);
   fprintf(file, "%.14e\n", tstart );
   fprintf(file, "%.14e\n", tstop );
   fprintf(file, "%.14e\n", current_t );
   fprintf(file, "%d\n",    size );
   fprintf(file, "%.14e\n", xstart );
   fprintf(file, "%.14e\n", xstop );
   
   for (i = 0; i < size; i++)
   {
      fprintf(file, "%.14e\n", values[i]);
   }
   
   fflush(file);
   fclose(file);
}
