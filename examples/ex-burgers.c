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


/*
   Example 01

   Compile with: make ex-burgers

   Sample run:   mpirun -np 2 ex-burgers

   Description:

   Burger's equation

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.h"
#include "braid_test.h"

/*--------------------------------------------------------------------------
 * My integration routines
 *--------------------------------------------------------------------------*/

/* can put anything in my app and name it anything as well */
typedef struct _braid_App_struct
{
   MPI_Comm  comm;
   double    tstart;
   double    tstop;
   int       ntime;
   double    xstart;
   double    xstop;
   double    xLeft;         /* this is the value of x on the left part of the domain for the initial condition, 
                               this is also the Dirichlet boundary condition on the left */
   double    xRight;        /* likewise, this is the value of x on the right */
   double    epsilon;       /* Diffusion coefficient for viscous Burgers*/
   int       max_iter_x[2]; /* length 2 array of max Newton iterations on fine and coarse grids */
   int       nspace;
   int       problem;       /* test problem, 0: linear advection, 1: Burger's equation */
   double    a;             /* if is_linear == 1, then use this as the linear advection constant */

} my_App;

/* can put anything in my vector and name it anything as well */
typedef struct _braid_Vector_struct
{
   int     size;
   double *values;

} my_Vector;

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


/* Helper function for Step: Tridiagonal system solver */
void 
solve_tridiag(double * x,
	      int N, 
	      double * a,
	      double * b, 
	      double * c) 
{
    /*
     solves Ax = v where A is a tridiagonal matrix consisting of vectors a, b, c
     x - initially contains v, returns solution x. indexed from 0 to N - 1
     N - length of vector x
     a - subdiagonal, indexed from 1 to N - 1
     b - main diagonal, indexed from 0 to N - 1
     c - superdiagonal, indexed from 0 to N - 2
     
     Note: contents of c will be modified
     */

    int ix;
    double m;
        
    c[0] = c[0] / b[0];
    x[0] = x[0] / b[0];
    
    /* loop from 1 to N - 2 inclusive, performing the forward sweep */
    for (ix = 1; ix < N - 1; ix++) 
    {
        m = 1.0 / (b[ix] - a[ix] * c[ix - 1]);
        c[ix] = c[ix] * m;
        x[ix] = (x[ix] - a[ix] * x[ix - 1]) * m;
    }

    x[N - 1] = (x[N - 1] - a[N - 1] * x[N - 2]) / (b[N - 1] - a[N - 1] * c[N - 2]);
    
    /* loop from N - 2 to 0 inclusive to perform the back substitution */
    for (ix = N - 2; ix >= 0; ix--)
        x[ix] = x[ix] - c[ix] * x[ix + 1];
}

/* helper function for my_Step*/
double compute_fstar(double uk_plus, double uk)
{
   return 0.5*( 0.5*uk_plus*uk_plus +  0.5*uk*uk) - 0.5*(fabs(0.5*uk) + fabs(0.5*uk_plus))  *(uk_plus - uk);
}

/* helper function for my_Step*/
double compute_fstar_linear(double uk_plus, double uk, double a)
{
   return 0.5*( a*uk_plus +  a*uk) - 0.5*fabs(a)*(uk_plus - uk);
}



int my_StepBE(braid_App        app,
              braid_Vector     ustop,
              braid_Vector     fstop,
              braid_Vector     u,
              braid_StepStatus status)
{
   int k, j, num_iter, level;
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   double deltaT, fstar_plus = 0.0, fstar_minus = 0.0, uk, uk_minus, uk_plus;
   double *u_old, *v, *a, *b, *c;
   double epsilon = app->epsilon;
   double deltaX = (app->xstop - app->xstart) / (u->size - 1.0);

   braid_StepStatusGetLevel(status, &level);
   if(level == 0) {
      num_iter = app->max_iter_x[0];
   }
   else{
      num_iter = app->max_iter_x[1];
   }
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   deltaT = tstop - tstart;

   /* allocate memory for vectors */
   v = (double *) malloc(sizeof(double)*u->size);
   a = (double *) malloc(sizeof(double)*u->size);
   b = (double *) malloc(sizeof(double)*u->size);
   c = (double *) malloc(sizeof(double)*u->size);
   
   /* copy u */
   u_old = (double *) malloc(sizeof(double)*u->size);
   for(k = 0; k < u->size; k++)
      u_old[k] = u->values[k];

   /* Newton Solver Loop */
   for(j = 0; j < num_iter; j++)
   {
      /* Create RHS vector */

      uk_minus = app->xLeft;
      uk = u->values[0];
      uk_plus = u->values[1];

      if(app->problem == 0)
      {
         fstar_plus = compute_fstar_linear(uk_plus, uk, app->a);
         fstar_minus = compute_fstar_linear(uk, uk_minus, app->a);
      }
      else if(app->problem == 1)
      {
         fstar_plus = compute_fstar(uk_plus, uk);
         fstar_minus = compute_fstar(uk, uk_minus); 
      }

      v[0] = uk - u_old[0] + (deltaT/deltaX)*(fstar_plus - fstar_minus) - 
            epsilon*(deltaT/(deltaX*deltaX))*(uk_plus - 2*uk + uk_minus);

      for(k = 1; k < u->size-1; k++)
      {
          uk_minus = u->values[k-1];
          uk = u->values[k];
          uk_plus = u->values[k+1];

          if(app->problem == 0)
          {
             fstar_plus = compute_fstar_linear(uk_plus, uk, app->a);
             fstar_minus = compute_fstar_linear(uk, uk_minus, app->a);
          }
          else if(app->problem == 1)
          {
             fstar_plus = compute_fstar(uk_plus, uk);
             fstar_minus = compute_fstar(uk, uk_minus); 
          }

          v[k] = uk - u_old[k] + (deltaT/deltaX)*(fstar_plus - fstar_minus) - 
               epsilon*(deltaT/(deltaX*deltaX))*(uk_plus - 2*uk + uk_minus);
      }

      uk_minus = u->values[u->size-2];
      uk = u->values[u->size-1];
      uk_plus = app->xRight;

      if(app->problem == 0)
      {
         fstar_plus = compute_fstar_linear(uk_plus, uk, app->a);
         fstar_minus = compute_fstar_linear(uk, uk_minus, app->a);
      }
      else if(app->problem == 1)
      {
         fstar_plus = compute_fstar(uk_plus, uk);
         fstar_minus = compute_fstar(uk, uk_minus); 
      }

      v[u->size-1] = uk - u_old[u->size-1] + (deltaT/deltaX)*(fstar_plus - fstar_minus) - 
            epsilon*(deltaT/(deltaX*deltaX))*(uk_plus - 2*uk + uk_minus);

      /* Create Jacobian - tridiag(a,b,c) */

      if(app->problem == 0)
      {
          uk_minus = app->xLeft;
          uk = u->values[0];
          uk_plus = u->values[1];

          b[0] = 1 + fabs(app->a)*(deltaT/deltaX) + 2*epsilon*(deltaT/(deltaX*deltaX));
          c[0] = 0.5*(app->a-fabs(app->a))*(deltaT/deltaX) - epsilon*(deltaT/(deltaX*deltaX));

    
          for(k = 1; k < u->size-1; k++)
          {
              uk_minus = u->values[k-1];
              uk = u->values[k];
              uk_plus = u->values[k+1];

              a[k] = -0.5*(app->a+fabs(app->a))*(deltaT/deltaX) - epsilon*(deltaT/(deltaX*deltaX));
              b[k] = 1 + fabs(app->a)*(deltaT/deltaX) + 2*epsilon*(deltaT/(deltaX*deltaX));
              c[k] =  0.5*(app->a-fabs(app->a))*(deltaT/deltaX) - epsilon*(deltaT/(deltaX*deltaX));
          }

          uk_minus = u->values[u->size-2];
          uk = u->values[u->size-1];
          uk_plus = app->xRight;

          a[u->size-1] = -0.5*(app->a+fabs(app->a))*(deltaT/deltaX) - epsilon*(deltaT/(deltaX*deltaX));
          b[u->size-1] = 1 + fabs(app->a)*(deltaT/deltaX) + 2*epsilon*(deltaT/(deltaX*deltaX));
      }
      else if(app->problem == 1)
      {
          uk_minus = app->xLeft;
          uk = u->values[0];
          uk_plus = u->values[1];

          b[0] = 1 + 0.25*(deltaT/deltaX)*(-dabs(uk)*(uk_plus - 2*uk + uk_minus) + fabs(uk_plus) + 2*fabs(uk) + fabs(uk_minus)) + 2*epsilon*(deltaT/(deltaX*deltaX));
          c[0] = (deltaT/deltaX)*(0.5*uk_plus - 0.25*dabs(uk_plus)*(uk_plus - uk) - 0.25*(fabs(uk_plus) + fabs(uk))) - epsilon*(deltaT/(deltaX*deltaX));
    
          for(k = 1; k < u->size-1; k++)
          {
              uk_minus = u->values[k-1];
              uk = u->values[k];
              uk_plus = u->values[k+1];

              a[k] = (deltaT/deltaX)*(-0.5*uk_minus + 0.25*dabs(uk_minus)*(uk - uk_minus) - 0.25*(fabs(uk) + fabs(uk_minus))) - epsilon*(deltaT/(deltaX*deltaX));
              b[k] = 1 + 0.25*(deltaT/deltaX)*(-dabs(uk)*(uk_plus - 2*uk + uk_minus) + fabs(uk_plus) + 2*fabs(uk) + fabs(uk_minus)) + 2*epsilon*(deltaT/(deltaX*deltaX));
              c[k] = (deltaT/deltaX)*(0.5*uk_plus - 0.25*dabs(uk_plus)*(uk_plus - uk) - 0.25*(fabs(uk_plus) + fabs(uk))) - epsilon*(deltaT/(deltaX*deltaX));
          }

          uk_minus = u->values[u->size-2];
          uk = u->values[u->size-1];
          uk_plus = app->xRight;

          a[u->size-1] = (deltaT/deltaX)*(-0.5*uk_minus + 0.25*dabs(uk_minus)*(uk - uk_minus) - 0.25*(fabs(uk) + fabs(uk_minus))) - epsilon*(deltaT/(deltaX*deltaX));
          b[u->size-1] = 1 + 0.25*(deltaT/deltaX)*(-dabs(uk)*(uk_plus - 2*uk + uk_minus) + fabs(uk_plus) + 2*fabs(uk) + fabs(uk_minus)) + 2*epsilon*(deltaT/(deltaX*deltaX)); 
      }

      solve_tridiag(v, u->size, a, b, c);
      for(k = 0; k < u->size; k++)
      {
         u->values[k] = u->values[k] - v[k];
      }
      
   }

   /* Free up */
   free(u_old);
   free(v);
   free(a);
   free(b);
   free(c);

   /* no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}


int
my_StepFE(braid_App        app,
          braid_Vector     ustop,
          braid_Vector     fstop,
          braid_Vector     u,
          braid_StepStatus status)
{
   int k;
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   double deltaT, fstar_plus = 0.0, fstar_minus = 0.0, uk, uk_minus, uk_plus;
   double *u_old;
   double epsilon = app->epsilon;
   double deltaX = (app->xstop - app->xstart) / (u->size - 1.0);
   
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   deltaT = tstop - tstart;
   
   /* copy u */
   u_old = (double *) malloc(sizeof(double)*u->size);
   for(k = 0; k < u->size; k++)
      u_old[k] = u->values[k];

   /* update interior of domain */
   for(k = 1; k < u->size-1; k++)
   {
      uk = u_old[k];
      uk_minus = u_old[k-1];
      uk_plus = u_old[k+1];
      
      if(app->problem == 0)
      {
         fstar_plus = compute_fstar_linear(uk_plus, uk, app->a);
         fstar_minus = compute_fstar_linear(uk, uk_minus, app->a);
      }
      else if(app->problem == 1)
      {
         fstar_plus = compute_fstar(uk_plus, uk);
         fstar_minus = compute_fstar(uk, uk_minus); 
      }

      u->values[k] = uk - (deltaT/deltaX)*(fstar_plus - fstar_minus) + 
               epsilon*(deltaT/(deltaX*deltaX))*(uk_plus - 2*uk + uk_minus);
   }

   /* update left boundary point (Dirichlet 0.0) */
   uk = u_old[0];
   uk_minus = app->xLeft;
   uk_plus = u_old[1];

   if(app->problem == 0)
   {
      fstar_plus = compute_fstar_linear(uk_plus, uk, app->a);
      fstar_minus = compute_fstar_linear(uk, uk_minus, app->a);
   }
   else if(app->problem == 1)
   {
      fstar_plus = compute_fstar(uk_plus, uk);
      fstar_minus = compute_fstar(uk, uk_minus); 
   }
   
   u->values[0] = uk - (deltaT/deltaX)*(fstar_plus - fstar_minus) + 
               epsilon*(deltaT/(deltaX*deltaX))*(uk_plus - 2*uk + uk_minus);


   
   /* update right boundary point (Dirichlet 1.0) */ 
   uk = u_old[u->size-1];
   uk_minus = u_old[u->size-2];
   uk_plus = app->xRight;

   if(app->problem == 0)
   {
      fstar_plus = compute_fstar_linear(uk_plus, uk, app->a);
      fstar_minus = compute_fstar_linear(uk, uk_minus, app->a);
   }
   else if(app->problem == 1)
   {
      fstar_plus = compute_fstar(uk_plus, uk);
      fstar_minus = compute_fstar(uk, uk_minus); 
   }
   
   u->values[u->size-1] = uk - (deltaT/deltaX)*(fstar_plus - fstar_minus) + 
               epsilon*(deltaT/(deltaX*deltaX))*(uk_plus - 2*uk + uk_minus);


   /* Free up */
   free(u_old);

   /* no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}

int
my_Residual(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     r,
            braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);


   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;
   int    i, nspace = (app->nspace);
   double xstart, xstop, x1, x2, x;
   
   u = (my_Vector *) malloc(sizeof(my_Vector));
   (u->size)   = nspace+1;
   (u->values) = (double *) malloc((u->size)*sizeof(double));

   xstart = (app->xstart);
   xstop  = (app->xstop);
   x1     = -1.0;
   x2     =  0.0;

   /* Initial guess */
   for (i = 0; i <= nspace; i++)
   {
      x = xstart + ((double)i/nspace)*(xstop - xstart);
      if (x < x1)
      {
         (u->values)[i] = app->xLeft;
      }
      else if (x < x2)
      {
         (u->values)[i] = app->xLeft - (x - x1)*(app->xLeft - app->xRight);
      }
      else
      {
         (u->values)[i] = app->xRight;
      }
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
   int i, size = (u->size);

   v = (my_Vector *) malloc(sizeof(my_Vector));
   (v->size)   = size;
   (v->values) = (double *) malloc(size*sizeof(double));
   for (i = 0; i < size; i++)
   {
      (v->values)[i] = (u->values)[i];
   }
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
   int i, size = (y->size);

   for (i = 0; i < size; i++)
   {
      (y->values)[i] = alpha*(x->values)[i] + beta*(y->values)[i];
   }

   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   int    i, size = (u->size);
   double dot = 0.0;

   for (i = 0; i < size; i++)
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
   MPI_Comm   comm   = (app->comm);
   double     tstart = (app->tstart);
   double     tstop  = (app->tstop);
   int        ntime  = (app->ntime);
   int        size = (u->size);
   int        i, index, myid;
   char       filename[255];
   FILE      *file;
   double     t;
   
   braid_AccessStatusGetT(astatus, &t);
   index = ((t-tstart) / ((tstop-tstart)/ntime) + 0.1);

   MPI_Comm_rank(comm, &myid);

   sprintf(filename, "%s.%07d.%05d", "ex-burgers.out", index, myid);
   file = fopen(filename, "w");
   fprintf(file, "%d\n", ntime+1);
   fprintf(file, "%.14e\n", app->tstart);
   fprintf(file, "%.14e\n", app->tstop);
   fprintf(file, "%d\n", app->nspace+1);
   fprintf(file, "%.14e\n", app->xstart);
   fprintf(file, "%.14e\n", app->xstop);
   for (i = 0; i < size; i++)
   {
      fprintf(file, "%.14e\n", (u->values)[i]);
   }
   fflush(file);
   fclose(file);

   return 0;
}

int
my_BufSize(braid_App  app,
           int       *size_ptr)
{
   int size = (app->nspace)+1;
   *size_ptr = (size+1)*sizeof(double);
   return 0;
}

int
my_BufPack(braid_App     app,
           braid_Vector  u,
           void         *buffer,
           braid_Int    *size_ptr)
{
   double *dbuffer = buffer;
   int i, size = (u->size);
   
   dbuffer[0] = size;
   for (i = 0; i < size; i++)
   {
      dbuffer[i+1] = (u->values)[i];
   }

   *size_ptr = (size+1)*sizeof(double);

   return 0;
}

int
my_BufUnpack(braid_App     app,
             void         *buffer,
             braid_Vector *u_ptr)
{
   my_Vector *u;
   double    *dbuffer = buffer;
   int        i, size;

   size = dbuffer[0];
   
   u = (my_Vector *) malloc(sizeof(my_Vector));
   (u->size)   = size;
   (u->values) = (double *) malloc(size*sizeof(double));
   for (i = 0; i < size; i++)
   {
      (u->values)[i] = dbuffer[i+1];
   }
   *u_ptr = u;

   return 0;
}


int
my_CoarsenLinear(braid_App              app,           
                   braid_Vector           fu,
                   braid_Vector          *cu_ptr,
                   braid_CoarsenRefStatus status)
{

   int i, csize, fidx, level;
   double *fvals = fu->values;
   my_Vector *v;
   
   csize = (fu->size - 1)/2 + 1;
   braid_CoarsenRefStatusGetLevel(status, &level);
   if(fu->size == 3)
      printf("Warning, coarsening in space down to 1 point on level %d!!\n", level);

   v = (my_Vector *) malloc(sizeof(my_Vector));
   (v->size)   = csize;
   (v->values) = (double *) malloc(csize*sizeof(double));
   for (i = 1; i < csize-1; i++)
   {
      fidx = 2*i;
      (v->values)[i] = 0.5*fvals[fidx] + 0.25*fvals[fidx+1] + 0.25*fvals[fidx-1];
   }

   /* Boundary Conditions */
   (v->values)[0] = (4./3.)*0.5*fvals[0] + (4./3.)*0.25*fvals[1];
   (v->values)[csize-1] = (4./3.)*0.5*fvals[fu->size-1] + (4./3.)*0.25*fvals[fu->size-2];
   
   *cu_ptr = v;
   
   return 0;
}

int
my_CoarsenOneSided(braid_App              app,           
                   braid_Vector           fu,
                   braid_Vector          *cu_ptr,
                   braid_CoarsenRefStatus status)
{

   int i, csize, fidx, level;
   double *fvals = fu->values;
   my_Vector *v;
   
   csize = (fu->size - 1)/2 + 1;
   braid_CoarsenRefStatusGetLevel(status, &level);
   if(fu->size == 3)
      printf("Warning, coarsening in space down to 1 point on level %d!!\n", level);
   
   /* Do averaging (assuming positive wave-speed) */
   v = (my_Vector *) malloc(sizeof(my_Vector));
   (v->size)   = csize;
   (v->values) = (double *) malloc(csize*sizeof(double));
   for (i = 1; i < csize-1; i++)
   {
      fidx = 2*i;
      (v->values)[i] = 0.5*fvals[fidx] + 0.5*fvals[fidx-1];
   }

   /* Boundary Conditions */
   (v->values)[0] = 0.5*fvals[0];
   (v->values)[csize-1] = 0.5*fvals[fu->size-1] + 0.5*fvals[fu->size-2];

   *cu_ptr = v;
   
   return 0;
}

int
my_InterpLinear(braid_App              app,           
                  braid_Vector           cu,
                  braid_Vector          *fu_ptr,
                  braid_CoarsenRefStatus status)
{

   int i, fsize;
   double *cvals = cu->values;
   my_Vector *v;
   
   fsize = (cu->size - 1)*2 + 1;

   v = (my_Vector *) malloc(sizeof(my_Vector));
   (v->size)   = fsize;
   (v->values) = (double *) malloc(fsize*sizeof(double));
   for (i = 1; i < fsize-1; i++)
   {
      if(i%2 == 1)
         (v->values)[i] = 0.5*cvals[i/2] + 0.5*cvals[(i+1)/2];
      else
         (v->values)[i] = cvals[i/2];
   }

   /* Boundary Conditions */
   (v->values)[0] = cvals[0];
   (v->values)[fsize-1] = cvals[cu->size-1];

   *fu_ptr = v;
   
   return 0;
}

int
my_InterpOneSided(braid_App              app,           
                  braid_Vector           cu,
                  braid_Vector          *fu_ptr,
                  braid_CoarsenRefStatus status)
{
   int i, fsize;
   double *cvals = cu->values;
   my_Vector *v;
   
   fsize = (cu->size - 1)*2 + 1;

   /* Inject values to fine grid */
   v = (my_Vector *) malloc(sizeof(my_Vector));
   (v->size)   = fsize;
   (v->values) = (double *) malloc(fsize*sizeof(double));
   for (i = 1; i < fsize-1; i++)
       (v->values)[i] = cvals[i/2];

   /* Boundary Conditions */
   (v->values)[0] = cvals[0];
   (v->values)[fsize-1] = cvals[cu->size-1];

   *fu_ptr = v;
   
   return 0;
}


/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
   braid_Core    core;
   my_App       *app;
   MPI_Comm      comm;
   double        tstart, tstop, epsilon, a;
   int           ntime;
   double        xstart, xstop, xLeft, xRight;
   int           nspace, problem;

   int           max_levels = 1;
   int           nrelax     = 1;
   int           nrelax0    = -1;
   double        tol        = 1.0e-07;
   int           cfactor    = 2;
   int           max_iter   = 30;
   int           fmg        = 0;
   int           scoarsen   = 0;
   int           res        = 0;
   int           stepper    = 0;
   int           max_iter_x[2];

   int           arg_index;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);

   /* dt = dx = 0.1 by default */
   comm   = MPI_COMM_WORLD;
   tstart =  0.0;
   tstop  =  2.0;
   ntime  =  60;
   xstart = -2.0;
   xstop  =  3.0;
   nspace =  16;
   problem = 1;
   a = 1.0;
   xLeft = 1.0;
   xRight = 0.0;
   epsilon = 0.0;
   max_iter_x[0] = 10;
   max_iter_x[1] = 10;
   
   /* Parse command line */

   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         int  myid;
         MPI_Comm_rank(comm, &myid);
         if ( myid == 0 )
         {
            printf("\n");
            printf("  -ml   <max_levels> : set max levels\n");
            printf("  -nu   <nrelax>     : set num F-C relaxations\n");
            printf("  -nx   <nspace>     : set num points in space\n");
            printf("  -nt   <ntime>      : set num points in time\n");
            printf("  -xL   <xLeft>      : set the left x-value (both as boundary condition and as initial condition)\n");
            printf("  -xR   <xRight>     : set the right x-value (both as boundary condition and as initial condition)\n");
            printf("  -eps  <epsilon>    : set the diffusion coefficient for a viscous problem \n");
            printf("  -prob <problem>    : set problem, problem 0 is linear advection, problem 1 is Burger's equation\n");
            printf("  -a    <a>          : set the speed of the linear advection problem (not used if problem is 1)\n");
            printf("  -st   <stepper>    : set the time stepper, 0: forward Euler, 1: backward Euler\n");
            printf("  -nu0  <nrelax>     : set num F-C relaxations on level 0\n");
            printf("  -tol  <tol>        : set stopping tolerance (scaled by sqrt(dt) sqrt(dx))\n");
            printf("  -cf   <cfactor>    : set coarsening factor\n");
            printf("  -mi   <max_iter>   : set max iterations\n");
            printf("  -mix  <mif  mic>   : set max Newton iterations on fine (mif) and all coarse levels (mic)\n");
            printf("  -sc   <scoarsen>   : use spatial coarsening by factor of 2 each level; must use 2^k sized grids, 1: bilinear, 2: one-sided, 3: ...\n");
            printf("  -fmg              : use FMG cycling\n");
            printf("  -res              : use my residual\n");
            printf("\n");
         }
         exit(1);
      }
      else if ( strcmp(argv[arg_index], "-ml") == 0 )
      {
         arg_index++;
         max_levels = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu") == 0 )
      {
         arg_index++;
         nrelax = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nt") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nx") == 0 )
      {
         arg_index++;
         nspace = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-prob") == 0 )
      {
         arg_index++;
         problem = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-a") == 0 )
      {
         arg_index++;
         a = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-xL") == 0 )
      {
         arg_index++;
         xLeft = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-xR") == 0 )
      {
         arg_index++;
         xRight = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-eps") == 0 )
      {
         arg_index++;
         epsilon = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-st") == 0 )
      {
         arg_index++;
         stepper = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu0") == 0 )
      {
         arg_index++;
         nrelax0 = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tol") == 0 ) 
      {
         arg_index++;
         tol = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-cf") == 0 )
      {
         arg_index++;
         cfactor = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mi") == 0 )
      {
         arg_index++;
         max_iter = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mix") == 0 )
      {
         arg_index++;
         max_iter_x[0] = atoi(argv[arg_index++]);
         max_iter_x[1] = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-fmg") == 0 )
      {
         arg_index++;
         fmg = 1;
      }
      else if ( strcmp(argv[arg_index], "-sc") == 0 )
      {
         arg_index++;
         scoarsen  = atoi(argv[arg_index++]);
      }

      else if ( strcmp(argv[arg_index], "-res") == 0 )
      {
         arg_index++;
         res = 1;
      }
      else
      {
         printf("ABORTING: incorrect command line parameter %s\n", argv[arg_index]);
         MPI_Finalize();
         return (0);

      }
   }

   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->comm)          = comm;
   (app->tstart)        = tstart;
   (app->tstop)         = tstop;
   (app->ntime)         = ntime;
   (app->xstart)        = xstart;
   (app->xstop)         = xstop;
   (app->nspace)        = nspace;
   (app->xLeft)         = xLeft;
   (app->xRight)        = xRight;
   (app->epsilon)       = epsilon;
   (app->max_iter_x[0]) = max_iter_x[0];
   (app->max_iter_x[1]) = max_iter_x[1];
   (app->problem)       = problem;
   (app->a)             = a;

   if(stepper == 0)
   {
      braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
             my_StepFE, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
             my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
   }
   else if(stepper == 1)
   {
      braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
             my_StepBE, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
             my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
   }
   else
   {
         printf("ABORTING: incorrect choice of stepper parameter 'st'\n");
         MPI_Finalize();
         return (0);
   }

   /* Scale tol by domain */
   tol = tol/( sqrt((tstop - tstart)/(ntime-1))*sqrt((xstop - xstart)/(nspace-1)) );

   braid_SetPrintLevel( core, 1);
   braid_SetMaxLevels(core, max_levels);
   braid_SetNRelax(core, -1, nrelax);
   if (nrelax0 > -1)
   {
      braid_SetNRelax(core,  0, nrelax0);
   }
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, cfactor);
   /*braid_SetCFactor(core,  0, 10);*/
   braid_SetMaxIter(core, max_iter);
   if (fmg)
   {
      braid_SetFMG(core);
   }
   if (res)
   {
      braid_SetResidual(core, my_Residual);
   }
   if (scoarsen)
   {
      if( fabs(pow(2.0, round(log2((double) nspace))) - nspace) > 0.5 )
      {
         printf("ABORTING: when using spatial coarsening, grid size must be a power of two.  Grid size = %d.\n",nspace);
         MPI_Finalize();
         return (0);
      }
      
      if (scoarsen == 1)
      {
         braid_SetSpatialCoarsen(core, my_CoarsenLinear);
         braid_SetSpatialRefine(core,  my_InterpLinear);
      }
      if (scoarsen == 2)
      {
         braid_SetSpatialCoarsen(core, my_CoarsenOneSided);
         braid_SetSpatialRefine(core,  my_InterpOneSided);
      }
      else
      {
         printf("Invalid scoarsen choice.  Ignoring this parameter\n");
      }

   }

   //braid_TestCoarsenRefine(app, 0, stdout, 0.0, 0.1, 0.2, my_Init,
   //      my_Access, my_Free, my_Clone, my_Sum, my_SpatialNorm, my_CoarsenBilinear, 
   //      my_InterpBilinear);

   braid_Drive(core);

   braid_Destroy(core);

   /* Finalize MPI */
   MPI_Finalize();

   return (0);
}
