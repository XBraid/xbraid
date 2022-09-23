/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory. Written by
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
 *
 * This file is part of XBraid. For support, post issues to the XBraid Github
 *page.
 *
 * This program is free software; you can redistribute it and/or modify it
 *under the terms of the GNU General Public License (as published by the Free
 *Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 *ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
 *FOR A PARTICULAR PURPOSE. See the terms and conditions of the GNU General
 *Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 *along with this program; if not, write to the Free Software Foundation, Inc.,
 *59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 ***********************************************************************EHEADER*/

/**
 * Example:       ex-07.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-07
 *
 * Help with:     ./ex-07 -help
 *
 * Sample run:    mpirun -np 2 ex-07 -rank 3
 *
 * Description:   Solve the chaotic Lorenz system,
 * 
 *                  {x' = sigma (y - x),    x(0) = 0.1,
 *                  {y' = x (rho - z) - y,  y(0) = 0.1,
 *                  {z' = xy - beta z,      z(0) = 0.1.
 * 
 *                with sigma=10, beta=8/3, and rho=28,
 *                and estimate the backward Lyapunov vectors along the trajectory
 *                parallel-in-time using MGRIT equipped with Delta correction.
 *                
 *                TODO: Give some high-level info about the algorithm.
 * 
 *                See (https://arxiv.org/abs/2208.12629) for a detailed introduction
 *                to the technique used here.
 * 
 * Visualize with: python viz-ex-07.py (requires NumPy and MatPlotLib)
 *
 **/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "braid.h"
#include "braid_test.h"

/*-------------------------*
 | linear algebra routines |
 *-------------------------*/

#define VecSize 3
typedef double VEC[VecSize];
typedef double MAT[VecSize][VecSize];

void
VecCopy(const VEC x, VEC y)
{
   for (int i = 0; i < VecSize; i++)
   {
      y[i] = x[i];
   }
}

void
VecSet(VEC x, double alpha)
{
   for (int i = 0; i < VecSize; i++)
   {
      x[i] = alpha;
   }
}

void
MatSetEye(MAT I)
{
   for (size_t i = 0; i < VecSize; i++)
   {
      VecSet(I[i], 0);
      I[i][i] = 1.;
   }
}

void
VecAxpy(const double alpha, const VEC x, VEC y)
{
   for (int i = 0; i < VecSize; i++)
   {
      y[i] = alpha * x[i] + y[i];
   }
}

void
MatAxpy(const double alpha, const MAT A, MAT B)
{
   for (int i = 0; i < VecSize; i++)
   {
      VecAxpy(alpha, A[i], B[i]);
   }
}

void
MatVec(const MAT A, VEC x)
{
   VEC tmp;
   VecCopy(x, tmp);

   for (int i = 0; i < VecSize; i++)
   {
      x[i] = 0.;
      for (int j = 0; j < VecSize; j++)
      {
         x[i] += A[i][j] * tmp[j];
      }
   }
}

/*-------------------------*
 | my integration routines |
 *-------------------------*/

/* can put anything in my app and name it anything as well */
typedef struct _braid_App_struct
{
   MPI_Comm comm;
   double tstart;
   double tstop;
   int ntime;
   FILE *file;    // saves the state vectors
   FILE *file_lv; // saves the Lyapunov vectors 

} my_App;

int
my_App_Init(my_App *app,
         MPI_Comm comm_,
         double tstart_,
         double tstop_,
         int ntime_,
         FILE *file)
{
   (app->comm) = comm_;
   (app->tstart) = tstart_;
   (app->tstop) = tstop_;
   (app->ntime) = ntime_;
   (app->file) = file;
   return 0;
}

/* can put anything in my vector and name it anything as well */
typedef struct _braid_Vector_struct
{
   VEC values;
} my_Vector;

const double sigma = 10.0;
const double beta = 8.0 / 3.0;
const double rho = 28.0;

// returns the derivative at point *u*
void
Lorenz(const VEC u, VEC u_out)
{
   u_out[0] = sigma * (u[1] - u[0]);
   u_out[1] = u[0] * (rho - u[2]) - u[1];
   u_out[2] = u[0] * u[1] - beta * u[2];
}

// returns the time-step scaled Jacobian at point *u*
void
Lorenz_du(const VEC u, MAT J)
{
   J[0][0] = -sigma;
   J[0][1] = sigma;
   J[0][2] = 0.;
   J[1][0] = rho - u[2];
   J[1][1] = -1.;
   J[1][2] = -u[0];
   J[2][0] = u[1];
   J[2][1] = u[0];
   J[2][2] = -beta;
}

/*
 * Computes a time-step with euler's method
 * Optionally computes the jacobian of the time-step and
 * writes the result in Psi (if not NULL)
 */
void
Euler(VEC u0, double h, MAT *F)
{
   VEC k1;
   Lorenz(u0, k1);

   if (F) // compute Jacobian of the time-step
   {
      MAT J;
      MatSetEye(*F);

      // F = I + hJ(u0)
      Lorenz_du(u0, J);
      MatAxpy(h, J, *F);
   }

   // u1 = u0 + hf(u0)
   VecAxpy(h, k1, u0);
}

int
my_Step(braid_App app,
        braid_Vector ustop,
        braid_Vector fstop,
        braid_Vector u,
        braid_StepStatus status)
{
   double tstart; /* current time */
   double tstop;  /* evolve to this time */
   double h;      /* dt value */
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);

   // get the number of Lyapunov vectors we need to propagate
   int rank; /* rank of Delta correction */
   braid_StepStatusGetDeltaRank(status, &rank);
   MAT *LinProp;

   if (rank == 0) // we are not propagating Lyapunov vectors
   {
      LinProp = NULL;
   }

   h = tstop - tstart;
   Euler((u->values), h, LinProp);

   for (int i = 0; i < rank; i++)
   {
      my_Vector *psi; // will point to the i-th basis vector

      braid_StepStatusGetBasisVec(status, &psi, i);

      // propagate the basis vector from tstart to tstop
      // Here, we are using the full Jacobian
      MatVec(*LinProp, psi->values);
   }

   /* no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}

int
my_Init(braid_App app, double t, braid_Vector *u_ptr)
{
   my_Vector *u;

   u = (my_Vector *)malloc(sizeof(my_Vector));
   if (t == 0.0)
   {
      /* Initial condition */
      VecSet(u->values, 0.1);
   }
   else
   {
      /* Initial guess at all other time points */
      VecSet(u->values, 0.1);
      u->values[0] = 0.2;
   }
   *u_ptr = u;

   return 0;
}

int
my_InitBasis(braid_App app, double t, int index, braid_Vector *u_ptr)
{
   my_Vector *u;

   u = (my_Vector *)malloc(sizeof(my_Vector));

   // initialize with the columns of the identity matrix
   // TODO: should we require these to be ortho-normal, or should we
   //       do Gram-schmidt on these after initializing?
   VecSet(u->values, 0.);
   u->values[index] = 1.;

   return 0;
}

int
my_Clone(braid_App app, braid_Vector u, braid_Vector *v_ptr)
{
   my_Vector *v;

   v = (my_Vector *)calloc(1, sizeof(my_Vector));
   VecCopy(u->values, v->values);

   *v_ptr = v;

   return 0;
}

int
my_Free(braid_App app, braid_Vector u)
{
   free(u);

   return 0;
}

int
my_Sum(braid_App app,
       double alpha,
       braid_Vector x,
       double beta,
       braid_Vector y)
{
   for (int i = 0; i < VecSize; i++)
   {
      y->values[i] = alpha * x->values[i] + beta * y->values[i];
   }

   return 0;
}

int
my_InnerProd(braid_App app, braid_Vector u, braid_Vector v, double *prod_ptr)
{
   double dot = 0.;

   for (int i = 0; i < VecSize; i++)
   {
      dot += (u->values[i]) * (v->values[i]);
   }
   *prod_ptr = dot;
   return 0;
}

int
my_SpatialNorm(braid_App app, braid_Vector u, double *norm_ptr)
{
   /* our inner product induces the Euclidean norm */
   my_InnerProd(app, u, u, norm_ptr);
   *norm_ptr = sqrt(*norm_ptr);
   return 0;
}

int
my_Access(braid_App app, braid_Vector u, braid_AccessStatus astatus)
{
   double tstart = (app->tstart);
   double tstop = (app->tstop);
   int ntime = (app->ntime);
   FILE *file = (app->file);
   int index, i;
   double t;

   braid_AccessStatusGetT(astatus, &t);
   index = ((t - tstart) / ((tstop - tstart) / ntime) + 0.1);

   fprintf(file, "%d", index);
   for (i = 0; i < VecSize; i++)
   {
      fprintf(file, " %.14e", (u->values[i]));
   }
   fprintf(file, "\n");
   fflush(file);

   /* write the lyapunov vectors to file */
   file = app->file_lv;
   int rank;
   braid_AccessStatusGetDeltaRank(astatus, &rank);

   fprintf(file, "%d", index);
   for (int j = 0; j < rank; j++)
   {
      my_Vector *psi;
      braid_AccessStatusGetBasisVec(astatus, &psi, j);
      for (i = 0; i < VecSize; i++)
      {
         fprintf(file, " %.14e", (psi->values[i]));
      }
      fprintf(file, "\n");
      fflush(file);
   }
   

   return 0;
}

int
my_BufSize(braid_App app, int *size_ptr, braid_BufferStatus bstatus)
{
   *size_ptr = VecSize * sizeof(double);
   
   /* tell braid the size of the basis vectors */
   braid_BufferStatusSetBasisSize(bstatus, VecSize * sizeof(double));
   /* Note: this isn't necessary here, but for more complicated examples the buffer-size may be different */
   return 0;
}

int
my_BufPack(braid_App app,
           braid_Vector u,
           void *buffer,
           braid_BufferStatus bstatus)
{
   double *dbuffer = buffer;
   int i;

   for (i = 0; i < VecSize; i++)
   {
      dbuffer[i] = (u->values[i]);
   }
   braid_BufferStatusSetSize(bstatus, VecSize * sizeof(double));

   return 0;
}

int
my_BufUnpack(braid_App app,
             void *buffer,
             braid_Vector *u_ptr,
             braid_BufferStatus status)
{
   double *dbuffer = buffer;
   my_Vector *u;
   int i;

   u = (my_Vector *)malloc(sizeof(my_Vector));
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

int
main(int argc, char *argv[])
{
   braid_Core core;
   my_App *app;
   MPI_Comm comm;
   double tstart, tstop;
   int ntime;

   int max_levels = 2;
   int nrelax = 1;
   int nrelax0 = -1;
   double tol = 1.0e-06;
   int cfactor = 2;
   int max_iter = 100;
   int fmg = 0;
   int test = 1;

   int arg_index, myid, nprocs;
   char filename[255];
   FILE *file;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);

   /* ntime time intervals */
   comm = MPI_COMM_WORLD;
   ntime = 1024;
   tstart = 0.0;
   tstop = 20.0;

   MPI_Comm_rank(comm, &myid);
   MPI_Comm_size(comm, &nprocs);

   /* Parse command line */

   arg_index = 1;
   while (arg_index < argc)
   {
      if (strcmp(argv[arg_index], "-help") == 0)
      {
         if (myid == 0)
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
            printf("  -test             : run wrapper tests\n");
            printf("\n");
         }
         exit(1);
      }
      else if (strcmp(argv[arg_index], "-ntime") == 0)
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-tstop") == 0)
      {
         arg_index++;
         tstop = atof(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-ml") == 0)
      {
         arg_index++;
         max_levels = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-nu") == 0)
      {
         arg_index++;
         nrelax = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-nu0") == 0)
      {
         arg_index++;
         nrelax0 = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-tol") == 0)
      {
         arg_index++;
         tol = atof(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-cf") == 0)
      {
         arg_index++;
         cfactor = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-mi") == 0)
      {
         arg_index++;
         max_iter = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-fmg") == 0)
      {
         arg_index++;
         fmg = 1;
      }
      else if (strcmp(argv[arg_index], "-test") == 0)
      {
         arg_index++;
         test = 1;
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
   app = (my_App *)malloc(sizeof(my_App));
   my_App_Init(app, MPI_COMM_WORLD, tstart, tstop, ntime, file);

   braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app, my_Step,
              my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, my_Access,
              my_BufSize, my_BufPack, my_BufUnpack, &core);

   /**
    * this allows braid to perform Gram-Schmidt orthonormalization on
    * a set of my_Vector objects, which is required for the estimation
    * of the Lyapunov vectors
    */
   //   braid_SetInnerProd(core, my_InnerProd);

   braid_SetDeltaCorrection(core, VecSize, my_InitBasis, my_InnerProd);

   braid_SetPrintLevel(core, 2);
   braid_SetMaxLevels(core, max_levels);
   braid_SetNRelax(core, -1, nrelax);
   if (nrelax0 > -1)
   {
      braid_SetNRelax(core, 0, nrelax0);
   }
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, cfactor);
   braid_SetMaxIter(core, max_iter);
   if (fmg)
   {
      braid_SetFMG(core);
   }

   if (test)
   {
      braid_TestInitAccess(app, MPI_COMM_WORLD, stdout, 0., my_Init, my_Access, my_Free);
      braid_TestInitAccess(app, MPI_COMM_WORLD, stdout, 1., my_Init, my_Access, my_Free);
      braid_TestClone(app, MPI_COMM_WORLD, stdout, 0., my_Init, my_Access, my_Free, my_Clone);
      braid_TestSum(app, MPI_COMM_WORLD, stdout, 0., my_Init, my_Access, my_Free, my_Clone, my_Sum);
      braid_TestSpatialNorm(app, MPI_COMM_WORLD, stdout, 0., my_Init, my_Free, my_Clone, my_Sum, my_SpatialNorm);
      braid_TestBuf(app, MPI_COMM_WORLD, stdout, 0., my_Init, my_Free, my_Sum, my_SpatialNorm, my_BufSize, my_BufPack, my_BufUnpack);
      braid_TestInnerProd(app, MPI_COMM_WORLD, stdout, 0., 1., my_Init, my_Free, my_Sum, my_InnerProd);

      braid_Destroy(core);
      MPI_Finalize();

      return 0;
   }

   braid_Drive(core);

   braid_Destroy(core);

   /* reorder the output file */
   if (myid == 0)
   {
      int ti, index, i, npoints = (ntime + 1);
      VEC *solution;
      fclose(file);
      solution = (VEC *)malloc(npoints * sizeof(VEC));
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

   return 0;
}
