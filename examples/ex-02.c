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
 * Example:       ex-02.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-02
 *
 * Help with:     ex-02 -help
 *
 * Sample run:    mpirun -np 2 ex-02
 *
 * Description:   Solves the 1D heat equation, with only parallelism in time
 *                Space-time domain:  [0, PI] x [0, 2*PI]
 *                Exact solution:     u(t,x) = sin(x)*cos(t)
 *
 **/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.h"
#include "braid_test.h"

#ifdef M_PI
   #define PI M_PI
#else
   #define PI 3.14159265358979
#endif

/*--------------------------------------------------------------------------
 * My integration routines
 *--------------------------------------------------------------------------*/

// have same interface and bells and whistles as f90 and expanded
//
// if you do spatial coarsening, then you'll have to re-compute the matrix stencils based on x.
//
// mem-check
//
// don't forget that you want a serial time evolution code as well...
//
// Get the residual option to work
//
// Put in Test All

/* can put anything in my app and name it anything as well */
typedef struct _braid_App_struct
{
   MPI_Comm  comm;
   double    tstart;       /* Define the temporal domain */
   double    tstop;
   int       ntime;
   double    xstart;       /* Define the spatial domain */
   double    xstop;
   int       nspace;
   double    xLeft;        /* the left spatial boundary condition (for all t-values) */
   double    xRight;       /* the right spatial boundary condition (for all t-values) */
   double    matrix[3];    /* the three point spatial discretization stencil */
   double *  g;            /* temporary vector for inversions and mat-vecs */
   double *  sc_info;      /* Runtime information on CFL's encountered and spatial discretizations used */
} my_App;

/* Can put anything in my vector and name it anything as well */
typedef struct _braid_Vector_struct
{
   int     size;
   double *values;

} my_Vector;


/* Exact solution */ 
double exact(double t, double x)
{
    return sin(x)*cos(t);
}

/* Forcing term F(t, x) for PDE,  u_t = u_xx + F(t,x) */
double forcing(double t, double x)
{
   return (-1.0)*sin(x)*sin(t) + sin(x)*cos(t);
}

/* Compute L2-norm of the error at a point in time */
double compute_error_norm(double * values, double xstart, double xstop, int nspace, double t)
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

int my_Step(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     fstop,
            braid_Vector     u,
            braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   int level, i;
   double cfl, x;
   double deltaT;
   double deltaX = (app->xstop - app->xstart) / (u->size - 1.0);

   braid_StepStatusGetLevel(status, &level);
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   deltaT = tstop - tstart;
   
   /* Store information on CFL and spatial coarsening for later output */
   (app->sc_info)[ (2*level) ] = deltaX;
   (app->sc_info)[ (2*level) + 1] = deltaT;
   
   /* Set up matrix stencil for 1D heat equation*/
   cfl = (deltaT/(deltaX*deltaX));
   app->matrix[0] = -cfl;
   app->matrix[1] = 1.0 + 2*cfl;
   app->matrix[2] = -cfl;

   /* Braid forcing */
   if(fstop != NULL)
   {
      for(i = 0; i < u->size; i++)
      {
         u->values[i] = u->values[i] + fstop->values[i];
      }
   }

   /* Apply boundary conditions */
   (u->values[0]) = (app->xLeft);
   (u->values[u->size-1]) = (app->xRight);

   /* PDE forcing */
   x = app->xstart;
   for(i = 0; i < u->size; i++)
   {
      u->values[i] = u->values[i] + deltaT*forcing(tstop, x);
      x = x + deltaX;
   }

   /* Backward Euler step */
   solve_tridiag(u->values, app->g, u->size, app->matrix);
   
   /* no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
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


int
my_Residual(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     r,
            braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   int i;
   double cfl, x;
   double deltaT;
   double deltaX = (app->xstop - app->xstart) / (ustop->size - 1.0);

   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   deltaT = tstop - tstart;
   
   /* Set up matrix stencil for 1D heat equation*/
   cfl = (deltaT/(deltaX*deltaX));
   app->matrix[0] = -cfl;
   app->matrix[1] = 1.0 + 2*cfl;
   app->matrix[2] = -cfl;

   /* Residual r = A*ustop - r - forcing */
   matvec_tridiag(ustop->values, app->g, ustop->size, app->matrix);
   x = app->xstart;
   for(i = 0; i < r->size; i++)
   {
      r->values[i] = app->g[i] - r->values[i] - deltaT*forcing(tstop, x);
      x = x + deltaX;
   }

   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;
   int    i; 
   int nspace = (app->nspace);
   double deltaX = (app->xstop - app->xstart) / (nspace - 1.0);
   double x = app->xstart;

   /* Allocate vector */
   u = (my_Vector *) malloc(sizeof(my_Vector));
   (u->size)   = nspace;
   (u->values) = (double *) malloc((u->size)*sizeof(double));

   
   /* Initialize vector (with correct boundary conditions) */
   if(t == 0.0)
   {
      for(i=0; i < nspace; i++)
      {
         (u->values)[i] = exact(t, x);
         x += deltaX;
      }
   }
   else
   {
      /* Use random values for u(t>0), this measures asymptotic convergence rate */
      srand(0);
      for(i=0; i < nspace; i++)
      {
         (u->values)[i] = ((double)rand())/RAND_MAX;
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
   int size = (u->size);
   int i;

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
   int i;
   int size = (y->size);

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
   int    i;
   int size   = (u->size);
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
   int        i, index, rank;
   char       filename[255];
   FILE      *file;
   double     t, error;
   
   braid_AccessStatusGetT(astatus, &t);
   braid_AccessStatusGetTIndex(astatus, &index);

   MPI_Comm_rank( (app->comm), &rank);

   sprintf(filename, "%s.%07d.%05d", "ex-02.out", index, rank);
   file = fopen(filename, "w");
   fprintf(file, "%d\n",    (app->ntime) +1);
   fprintf(file, "%.14e\n", (app->tstart) );
   fprintf(file, "%.14e\n", (app->tstop) );
   fprintf(file, "%d\n",    (u->size) );
   fprintf(file, "%.14e\n", (app->xstart) );
   fprintf(file, "%.14e\n", (app->xstop) );
   for (i = 0; i < (u->size); i++)
   {
      fprintf(file, "%.14e\n", (u->values)[i]);
   }
   fflush(file);
   fclose(file);

   if(index == app->ntime)
   {
      error = compute_error_norm(u->values, app->xstart, app->xstop, u->size, t);
      printf("\n  Discretization error at final time:  %1.4e\n", error);
   }

   return 0;
}

int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
   int size = (app->nspace);
   *size_ptr = (size+1)*sizeof(double);
   return 0;
}

int
my_BufPack(braid_App           app,
           braid_Vector        u,
           void               *buffer,
           braid_BufferStatus  bstatus)
{
   double *dbuffer = buffer;
   int i, size = (u->size);
   
   dbuffer[0] = size;
   for (i = 0; i < size; i++)
   {
      dbuffer[i+1] = (u->values)[i];
   }

   braid_BufferStatusSetSize( bstatus,  (size+1)*sizeof(double));

   return 0;
}

int
my_BufUnpack(braid_App           app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus  bstatus)
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


/* Bilinear Coarsening */
int
my_CoarsenLinear(braid_App              app,           
                 braid_Vector           fu,
                 braid_Vector          *cu_ptr,
                 braid_CoarsenRefStatus status)
{

   int i, csize, fidx, level;
   double *fvals = fu->values;
   my_Vector *v;
   
   /* This returns the level for fu */
   braid_CoarsenRefStatusGetLevel(status, &level);
   
   /* There are only enough points in space to coarsen on the first
    * floor(log2(app->nspace)) levels.  The finest level is 0 */
   if( level < floor(log2(app->nspace)) )
   {
      /* Interpolate to interior points */
      csize = (fu->size - 1)/2 + 1;
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
   }
   else
   {
      /* No coarsening, clone the vector */
      my_Clone(app, fu, &v);
   }   


   *cu_ptr = v;

   return 0;
}

/* Bilinear interpolation */
int
my_InterpLinear(braid_App              app,           
                braid_Vector           cu,
                braid_Vector          *fu_ptr,
                braid_CoarsenRefStatus status)
{

   int i, fsize, level;
   double *cvals = cu->values;
   my_Vector *v;

   /* This returns the level for fu_ptr */
   braid_CoarsenRefStatusGetLevel(status, &level);
   
   /* There are only enough points in space to coarsen on the first
    * floor(log2(app->nspace)) levels.  The finest level is 0 */
   if( level < floor(log2(app->nspace)) )
   {
      /* Interpolate to interior points */
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
   }
   else
   {
      /* No refinement, clone the vector */
      my_Clone(app, cu, &v); 
   }

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
   double        tstart, tstop, dx, dt;
   int           i, ntime, rank;
   double        xstart, xstop, xLeft, xRight;
   int           nspace;

   int           max_levels    = 2;
   int           nrelax        = 1;
   int           nrelax0       = -1;
   int           skip          = 0;
   double        tol           = 1.0e-07;
   int           cfactor       = 2;
   int           max_iter      = 30;
   int           fmg           = 0;
   int           nfmg_Vcyc     = 1;
   int           scoarsen      = 0;
   int           res           = 0;

   int           arg_index;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   comm   = MPI_COMM_WORLD;
   MPI_Comm_rank(comm, &rank);

   tstart =  0.0;
   tstop  =  2*PI;
   ntime  =  64;
   xstart =  0.0;
   xstop  =  PI;
   nspace =  32;
   xLeft  =  exact(0, 0.0);
   xRight =  exact(0, PI);
   
   /* Parse command line */

   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         if ( rank == 0 )
         {
            printf("\n");
            printf(" Solve the 1D heat equation on space-time domain:  [0, PI] x [0, 2*PI]\n");
            printf(" with exact solution u(t,x) = sin(x)*cos(t) \n\n");
            printf("  -ntime <ntime>       : set num points in time\n");
            printf("  -nspace <nspace>     : set num points in space\n");
            printf("  -ml   <max_levels>   : set max levels\n");
            printf("  -nu   <nrelax>       : set num F-C relaxations\n");
            printf("  -skip <set_skip>     : set skip relaxations on first down-cycle; 0: no skip;  1: skip\n");
            printf("  -nu0  <nrelax>       : set num F-C relaxations on level 0\n");
            printf("  -tol  <tol>          : set stopping tolerance (scaled later by sqrt(dt) sqrt(dx))\n");
            printf("  -cf   <cfactor>      : set coarsening factor\n");
            printf("  -mi   <max_iter>     : set max iterations\n");
            printf("  -fmg  <nfmg_Vcyc>    : use FMG cycling, nfmg_Vcyc V-cycles at each fmg level\n");
            printf("  -sc                  : use spatial coarsening by factor of 2 each level\n");
            printf("  -res                 : use my residual\n");
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
      else if ( strcmp(argv[arg_index], "-skip") == 0 )
      {
         arg_index++;
         skip = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ntime") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nspace") == 0 )
      {
         arg_index++;
         nspace = atoi(argv[arg_index++]);
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
      else if ( strcmp(argv[arg_index], "-fmg") == 0 )
      {
         arg_index++;
         fmg = 1;
         nfmg_Vcyc = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-sc") == 0 )
      {
         arg_index++;
         scoarsen  = 1;
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
   (app->g)             = (double*) malloc( nspace*sizeof(double) );
   (app->comm)          = comm;
   (app->tstart)        = tstart;
   (app->tstop)         = tstop;
   (app->ntime)         = ntime;
   (app->xstart)        = xstart;
   (app->xstop)         = xstop;
   (app->nspace)        = nspace;
   (app->xLeft)         = xLeft;
   (app->xRight)        = xRight;

   /* Initialize the storage structure for recording spatial coarsening information */ 
   app->sc_info = (double*) malloc( 2*max_levels*sizeof(double) );
   for( i = 0; i < 2*max_levels; i++) {
      app->sc_info[i] = -1.0;
   }
   
   /* Initialize Braid */
   braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
          my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
          my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

   /* Scale tol by domain */
   tol = tol/( sqrt((tstop - tstart)/(ntime-1))*sqrt((xstop - xstart)/(nspace-1)) );

   /* Set Braid options */
   braid_SetPrintLevel( core, 1);
   braid_SetMaxLevels(core, max_levels);
   braid_SetSkip(core, skip);
   braid_SetNRelax(core, -1, nrelax);
   if (nrelax0 > -1)
   {
      braid_SetNRelax(core,  0, nrelax0);
   }
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, cfactor);
   braid_SetMaxIter(core, max_iter);
   if (fmg)
   {
      braid_SetFMG(core);
      braid_SetNFMGVcyc(core, nfmg_Vcyc);
   }
   if (res)
   {
      braid_SetResidual(core, my_Residual);
   }
   if (scoarsen)
   {
      if( fabs(log2(nspace - 1.0) - round(log2(nspace - 1.0))) > 1e-10 )
      {
         if(rank == 0)
         {
            fprintf(stderr, "\nWarning!\nWhen using spatial coarsening, spatial grids "
                            "must be power of 2 + 1, \ni.e., nspace = 2^k + 1.  Your "
                            "spatial grid is of size %d.  Spatial \ncoarsening is "
                            "therefore being ignored.\n\n", nspace);
         }
      }
      else
      {
         braid_SetSpatialCoarsen(core, my_CoarsenLinear);
         braid_SetSpatialRefine(core,  my_InterpLinear);
      }
   }

   braid_Drive(core);

   if(rank == 0)
   {
      printf("\n-----------------------------------------------------------------\n"); 
      printf("-----------------------------------------------------------------\n\n"); 
      printf( " Per level diagnostic information \n\n");
      printf("level       dx          dt        dt/dx^2\n"); 
      printf("-----------------------------------------------------------------\n"); 
      for( i = 0; i < max_levels; i++)
      {
         dx = app->sc_info[i*2];
         dt = app->sc_info[i*2+1];
         if (dx == -1){
            break;
         }
         printf(" %2d   |   %1.2e    %1.2e    %1.2e\n", i, dx, dt, dt/(dx*dx) ); 
      }
      printf( "\n" );
   }

   braid_Destroy(core);
   free( app->sc_info);
   free( app->g);
   free( app );

   /* Finalize MPI */
   MPI_Finalize();

   return (0);
}
