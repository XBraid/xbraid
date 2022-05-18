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
 * Driver:        drive-lorenz.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make drive-lorenz
 *
 * Help with:     drive-lorenz -help
 *
 * Sample run:    mpirun -np 2 drive-lorenz
 *
 * Description:   Solves the Lorenz problem 
 *
 *                Use vis-burgers-1D.py to visualize.
 **/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.h"

/*--------------------------------------------------------------------------
 * My integration routines
 *--------------------------------------------------------------------------*/

#define VecSize 3
typedef double Vec[VecSize];

/* can put anything in my app and name it anything as well */
typedef struct _braid_App_struct
{
   MPI_Comm  comm;
   double    tstart;
   double    tstop;
   int       ntime;
   FILE     *file;

} my_App;

/* can put anything in my vector and name it anything as well */
typedef struct _braid_Vector_struct
{
   Vec values;

} my_Vector;

double sigma = 10.0;
double beta  = 8.0/3.0;
double rho   = 28.0;

int
VecSet(Vec     x,
       double  alpha)
{
   int  i;

   for (i = 0; i < VecSize; i++)
   {
      x[i] = alpha;
   }

   return 0;
}

int
VecAxpy(double  alpha,
        Vec     x,
        double  beta,
        Vec     y)
{
   int  i;

   for (i = 0; i < VecSize; i++)
   {
      y[i] = alpha*x[i] + beta*y[i];
   }

   return 0;
}

int
RK4(double  tstart,
    Vec     ystart,
    double  h,
    int     (*func)(double t, Vec y, Vec ynew))
{
   int     s = 4;
   double  b[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
   double  c[3] = {1.0/2.0, 1.0/2.0,   1.0};
   double  a[3][3] = {{1.0/2.0,     0.0,   0.0},
                      {    0.0, 1.0/2.0,   0.0},
                      {    0.0,     0.0,   1.0}};
   Vec     k[4], y;
   int     i, j;

   VecSet(y, 0);
   (*func)(tstart, ystart, k[0]);
   for (i = 1; i < s; i++)
   {
      double  t;
      
      t = tstart + h*c[i-1];
      VecAxpy(a[i-1][0], k[0], 0, y);
      for (j = 1; j < i; j++)
      {
         VecAxpy(a[i-1][j], k[j], 1, y);
      }
      VecAxpy(1, ystart, h, y); // y = ystart + h*y
      (*func)(t, y, k[i]);
   }
   VecAxpy(b[0], k[0], 0, y); // y = b[0]*k[0]
   for (i = 1; i < s; i++)
   {
      VecAxpy(b[i], k[i], 1, y); // y += b[i]*k[i]
   }
   VecAxpy(h, y, 1, ystart); // ystop = ystart + h*y
   
   return 0;
}

int
Lorenz(double t,
       Vec    y,
       Vec    ynew)
{
   ynew[0] = sigma*(y[1] - y[0]);
   ynew[1] = y[0]*(rho - y[2]) - y[1];
   ynew[2] = y[0]*y[1] - beta*y[2];

   return 0;
}

int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time */
   double h;                  /* dt value */
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);

   h = tstop - tstart;
   /* On the finest grid, each value is half the previous value */
   RK4(tstart, (u->values), h, Lorenz);

   if (fstop != NULL)
   {
      /* Nonzero rhs */
      VecAxpy(1, (fstop->values), 1, (u->values));
   }

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

   /* On the finest grid, each value is half the previous value */
//   (r->value) = (ustop->value) - pow(0.5, tstop-tstart)*(r->value);

   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   if (t == 0.0)
   {
      /* Initial guess */
      VecSet((u->values), 0.1);
   }
   else
   {
      /* Initialize all other time points */
      VecSet((u->values), 0.1);
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

   v = (my_Vector *) calloc(1, sizeof(my_Vector));
   VecAxpy(1, (u->values), 0, (v->values));
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
   VecAxpy(alpha, (x->values), beta, (y->values));

   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   int    i;
   double dot = 0;

   for (i = 0; i < VecSize; i++)
   {
      dot += (u->values[i])*(u->values[i]);
   }
   *norm_ptr = sqrt(dot);

   return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   double     tstart = (app->tstart);
   double     tstop  = (app->tstop);
   int        ntime  = (app->ntime);
   FILE      *file   = (app->file);
   int        index, i;
   double     t;
   
   braid_AccessStatusGetT(astatus, &t);
   index = ((t-tstart) / ((tstop-tstart)/ntime) + 0.1);

   fprintf(file, "%d", index);
   for (i = 0; i < VecSize; i++)
   {
      fprintf(file, " %.14e", (u->values[i]));
   }
   fprintf(file, "\n");
   fflush(file);

   return 0;
}

int
my_BufSize(braid_App  app,
           int       *size_ptr,
           braid_BufferStatus  bstatus)
{
   *size_ptr = VecSize*sizeof(double);
   return 0;
}

int
my_BufPack(braid_App     app,
           braid_Vector  u,
           void         *buffer,
           braid_BufferStatus  bstatus)
{
   double *dbuffer = buffer;
   int     i;

   for (i = 0; i < VecSize; i++)
   {
      dbuffer[i] = (u->values[i]);
   }
   braid_BufferStatusSetSize( bstatus, VecSize*sizeof(double));

   return 0;
}

int
my_BufUnpack(braid_App           app,
             void                *buffer,
             braid_Vector        *u_ptr,
             braid_BufferStatus  status)
{
   double    *dbuffer = buffer;
   my_Vector *u;
   int        i;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   for (i = 0; i < VecSize; i++)
   {
      (u->values[i]) = dbuffer[i];
   }
   *u_ptr = u;

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
   double        tstart, tstop;
   int           ntime;

   int           max_levels = 1;
   int           nrelax     = 1;
   int           nrelax0    = -1;
   double        tol        = 1.0e-06;
   int           cfactor    = 2;
   int           max_iter   = 100;
   int           fmg        = 0;
   int           res        = 0;

   int           arg_index, myid, nprocs;
   char          filename[255];
   FILE         *file;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);

   /* ntime time intervals with spacing 1 */
   comm   = MPI_COMM_WORLD;
   ntime  = 10000;
   tstart =   0.0;
   tstop  = 100.0;

   MPI_Comm_rank(comm, &myid);
   MPI_Comm_size(comm, &nprocs);
   
   /* Parse command line */

   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         if ( myid == 0 )
         {
            printf("\n");
            printf("  -ntime <ntime>    : set num time points (default %d)\n", ntime);
            printf("  -tstop <tstop>    : set end time (default %lf)\n", tstop);
            printf("  -ml  <max_levels> : set max levels\n");
            printf("  -nu  <nrelax>     : set num F-C relaxations\n");
            printf("  -nu0 <nrelax>     : set num F-C relaxations on level 0\n");
            printf("  -tol <tol>        : set stopping tolerance\n");
            printf("  -cf  <cfactor>    : set coarsening factor\n");
            printf("  -mi  <max_iter>   : set max iterations\n");
            printf("  -fmg              : use FMG cycling\n");
            printf("  -res              : use my residual\n");
            printf("\n");
         }
         exit(1);
      }
      else if ( strcmp(argv[arg_index], "-ntime") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tstop") == 0 ) 
      {
         arg_index++;
         tstop = atof(argv[arg_index++]);
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
      }
      else if ( strcmp(argv[arg_index], "-res") == 0 )
      {
         arg_index++;
         res = 1;
      }
      else
      {
         arg_index++;
         /*break;*/
      }
   }

   /* initialize the output file */
   sprintf(filename, "%s.%05d", "ex-lorenz.out", myid);
   file = fopen(filename, "w");

   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->comm)   = comm;
   (app->tstart) = tstart;
   (app->tstop)  = tstop;
   (app->ntime)  = ntime;
   (app->file)   = file;

   braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
             my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
             my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

   braid_SetPrintLevel( core, 2);
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

   braid_Drive(core);

   braid_Destroy(core);

   /* reorder the output file */
   if (myid == 0)
   {
      int   ti, index, i, npoints = (ntime+1);
      Vec  *solution;
      fclose(file);
      solution = (Vec *) malloc(npoints*sizeof(Vec));
      file = fopen(filename, "r");
      for (ti = 0; ti < npoints; ti++)
      {
         fscanf(file, "%d", &index);
         for (i = 0; i < VecSize; i++)
         {
            fscanf(file, "%le", &solution[index][i]);
         }
      }
      fclose(file);
      file = fopen("ex-lorenz.out", "w");
      for (ti = 0; ti < ntime; ti++)
      {
         for (i = 0; i < VecSize; i++)
         {
            fprintf(file, " %.14e", solution[ti][i]);
         }
         fprintf(file, "\n");
      }
      free(solution);
      fclose(file);
   }

   /* Finalize MPI */
   MPI_Finalize();
   
   return (0);
}
