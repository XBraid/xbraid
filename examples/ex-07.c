/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory. Written by
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
 *
 * This file is part of XBraid. For support, post issues to the XBraid Github
 * page.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the terms and conditions of the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
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
 *                  {x' = σ (y - x),        x(0) = 9.614521163788712,
 *                  {y' = x (ρ - z) - y,    y(0) = 17.519127090401067,
 *                  {z' = xy - β z,         z(0) = 13.94230712158098,
 *
 *                with sigma=10, beta=8/3, and rho=28,
 *                and estimate the backward Lyapunov vectors along the trajectory
 *                parallel-in-time using MGRIT equipped with Delta correction.
 *                The initial condition should be a point close to the strange attractor.
 *
 *                See the user manual for more information about this example.
 *                See (https://arxiv.org/abs/2208.12629) for a more detailed introduction
 *                to the Delta correction technique used here.
 *
 * Visualize with: python viz-ex-07.py (requires NumPy and MatPlotLib)
 *
 *                 This will plot the trajectory, along with any estimate Lyapunov vectors
 *                 computed by XBraid.
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

const double sigma = 10.0;
const double beta = 8.0 / 3.0;
const double rho = 28.0;

void VecCopy(const VEC x, VEC y)
{
   for (int i = 0; i < VecSize; i++)
   {
      y[i] = x[i];
   }
}

void VecSet(VEC x, double alpha)
{
   for (int i = 0; i < VecSize; i++)
   {
      x[i] = alpha;
   }
}

void MatSetEye(MAT I)
{
   for (size_t i = 0; i < VecSize; i++)
   {
      VecSet(I[i], 0);
      I[i][i] = 1.;
   }
}

void VecAxpy(const double alpha, const VEC x, VEC y)
{
   for (int i = 0; i < VecSize; i++)
   {
      y[i] = alpha * x[i] + y[i];
   }
}

void MatAxpy(const double alpha, const MAT A, MAT B)
{
   for (int i = 0; i < VecSize; i++)
   {
      VecAxpy(alpha, A[i], B[i]);
   }
}

void MatVec(const MAT A, VEC x)
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

// returns the derivative at point *u*
void Lorenz(const VEC u, VEC u_out)
{
   u_out[0] = sigma * (u[1] - u[0]);
   u_out[1] = u[0] * (rho - u[2]) - u[1];
   u_out[2] = u[0] * u[1] - beta * u[2];
}

// returns the Jacobian at point *u*
void Lorenz_du(const VEC u, MAT J)
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
 * writes the result in F (if not NULL)
 */
void Euler(VEC u0, double h, MAT *F)
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

/*--------------------------------------------------------------------------
 * User-defined routines and structures required by Braid
 *--------------------------------------------------------------------------*/

/* can put anything in my app and name it anything as well */
typedef struct _braid_App_struct
{
   MPI_Comm comm;
   double tstart;
   double tstop;
   int ntime;
   int rank;
   double lyap_exps[VecSize]; // temporary storage for lyapunov exponents
   FILE *file;                // saves the state vectors
   FILE *file_lv;             // saves the Lyapunov vectors and exponents

} my_App;

int my_App_Init(my_App *app,
                MPI_Comm comm_,
                double tstart_,
                double tstop_,
                int ntime_,
                int rank_,
                FILE *file,
                FILE *file_lv)
{
   (app->comm) = comm_;
   (app->tstart) = tstart_;
   (app->tstop) = tstop_;
   (app->ntime) = ntime_;
   (app->rank) = rank_;
   (app->file) = file;
   (app->file_lv) = file_lv;
   for (int i = 0; i < rank_; i++)
   {
      (app->lyap_exps)[i] = 0.;
   }

   return 0;
}

/* can put anything in my vector and name it anything as well */
typedef struct _braid_Vector_struct
{
   VEC values;
} my_Vector;

int my_Step(braid_App app,
            braid_Vector ustop,
            braid_Vector fstop,
            braid_Vector u,
            braid_StepStatus status)
{
   /* for Delta correction, the user must propagate the solution vector (as in a traditional Braid code)
    * as well as the Lyapunov vectors. The Lyapunov vectors are available through the StepStatus structure,
    * and are propagated by the Jacobian of the time-step function. (see below)
    */
   double tstart; /* current time */
   double tstop;  /* evolve to this time */
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);

   double h; /* dt value */
   h = tstop - tstart;

   // get the number of Lyapunov vectors we need to propagate
   int rank; /* rank of Delta correction */
   braid_StepStatusGetDeltaRank(status, &rank);
   MAT Jacobian = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

   if (rank > 0) // we are propagating Lyapunov vectors
   {
      Euler((u->values), h, &Jacobian);
   }
   else
   {
      Euler((u->values), h, NULL);
   }

   for (int i = 0; i < rank; i++)
   {
      // get a reference to the ith Lyapunov vector
      my_Vector *psi;
      braid_StepStatusGetBasisVec(status, &psi, i);

      // propagate the vector from tstart to tstop
      if (psi)
      {
         MatVec(Jacobian, psi->values);
      }
   }

   /* no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}

int my_Init(braid_App app, double t, braid_Vector *u_ptr)
{
   my_Vector *u;

   u = (my_Vector *)malloc(sizeof(my_Vector));

   /* Initial condition */
   u->values[0] = 9.614521163788712;
   u->values[1] = 17.519127090401067;
   u->values[2] = 13.94230712158098;

   *u_ptr = u;

   return 0;
}

int my_InitBasis(braid_App app, double t, int index, braid_Vector *u_ptr)
{
   /*
    *  For Delta correction, an initial guess is needed for the Lyapunov basis vectors.
    *  This function initializes the basis vector with spatial index *index* at time *t*.
    *  Note that the vectors at each index *index* must be linearly independent.
    */
   my_Vector *u;

   u = (my_Vector *)malloc(sizeof(my_Vector));

   // initialize with the columns of the identity matrix
   VecSet(u->values, 0.);
   u->values[index] = 1.;

   *u_ptr = u;

   return 0;
}

int my_Clone(braid_App app, braid_Vector u, braid_Vector *v_ptr)
{
   my_Vector *v;

   v = (my_Vector *)calloc(1, sizeof(my_Vector));
   VecCopy(u->values, v->values);

   *v_ptr = v;

   return 0;
}

int my_Free(braid_App app, braid_Vector u)
{
   free(u);

   return 0;
}

int my_Sum(braid_App app,
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

int my_InnerProd(braid_App app, braid_Vector u, braid_Vector v, double *prod_ptr)
{
   /*
    *  For Delta correction, braid needs to be able to compute an inner product between two user vectors,
    *  which is used to project the user's vector onto the Lyapunov basis for low-rank Delta correction.
    *  This function should define a valid inner product between the vectors *u* and *v*.
    */
   double dot = 0.;

   for (int i = 0; i < VecSize; i++)
   {
      dot += (u->values[i]) * (v->values[i]);
   }
   *prod_ptr = dot;
   return 0;
}

int my_SpatialNorm(braid_App app, braid_Vector u, double *norm_ptr)
{
   /* our inner product induces the Euclidean norm */
   my_InnerProd(app, u, u, norm_ptr);
   *norm_ptr = sqrt(*norm_ptr);
   return 0;
}

int my_Access(braid_App app, braid_Vector u, braid_AccessStatus astatus)
{
   FILE *file = (app->file);
   int index, i;
   double t;

   braid_AccessStatusGetT(astatus, &t);
   braid_AccessStatusGetTIndex(astatus, &index);

   fprintf(file, "%d", index);
   for (i = 0; i < VecSize; i++)
   {
      fprintf(file, " %.14e", (u->values[i]));
   }
   fprintf(file, "\n");
   fflush(file);

   /* write the lyapunov vectors to file */
   file = app->file_lv;
   int local_rank, num_exp;
   braid_AccessStatusGetDeltaRank(astatus, &local_rank);
   num_exp = local_rank;
   double *exponents = malloc(local_rank * sizeof(double));
   if (num_exp > 0)
   {
      braid_AccessStatusGetLocalLyapExponents(astatus, exponents, &num_exp);
   }

   fprintf(file, "%d", index);
   for (int j = 0; j < local_rank; j++)
   {
      my_Vector *psi;
      braid_AccessStatusGetBasisVec(astatus, &psi, j);
      if (psi)
      {
         if (j < num_exp)
         {
            (app->lyap_exps)[j] += exponents[j];
            fprintf(file, " %.14e", exponents[j]);
         }
         else
         {
            fprintf(file, " %.14e", 0.);
         }
         for (i = 0; i < VecSize; i++)
         {
            fprintf(file, " %.14e", (psi->values[i]));
         }
      }
   }
   fprintf(file, "\n ");
   fflush(file);
   free(exponents);

   return 0;
}

int my_BufSize(braid_App app, int *size_ptr, braid_BufferStatus bstatus)
{
   /* Tell Braid the size of a state vector */
   *size_ptr = VecSize * sizeof(double);

   /* 
    * In contrast with traditional Braid, you may also specify the size of a single Lyapunov basis vector, 
    * in case it is different from the size of a state vector.
    * Note: this isn't necessary here, but for more complicated applications this size may be different.
    */
   braid_BufferStatusSetBasisSize(bstatus, VecSize * sizeof(double));
   return 0;
}

int my_BufPack(braid_App app,
               braid_Vector u,
               void *buffer,
               braid_BufferStatus bstatus)
{
   /* This function is used to pack both the state vector and Lyapunov basis vectors */
   double *dbuffer = buffer;

   for (int i = 0; i < VecSize; i++)
   {
      dbuffer[i] = (u->values[i]);
   }
   braid_BufferStatusSetSize(bstatus, VecSize * sizeof(double));

   return 0;
}

int my_BufUnpack(braid_App app,
                 void *buffer,
                 braid_Vector *u_ptr,
                 braid_BufferStatus status)
{
   /* This function is used to unpack both the state vector and Lyapunov basis vectors */
   double *dbuffer = buffer;
   my_Vector *u;

   u = (my_Vector *)malloc(sizeof(my_Vector));
   for (int i = 0; i < VecSize; i++)
   {
      (u->values[i]) = dbuffer[i];
   }
   *u_ptr = u;

   return 0;
}

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
   braid_Core core;
   my_App *app;
   MPI_Comm comm;
   double tstart, tstop;
   int ntime;

   double tlyap = 4.;
   int max_levels = 3;
   int nrelax = 0;
   int nrelax0 = -1;
   double tol = 1.0e-06;
   int cfactor = 2;
   int max_iter = 100;
   int fmg = 0;
   int test = 0;
   int lyap = 1;
   int relax_lyap = 0;
   int delta_rank = 3;
   int defer_lvl = 0;
   int defer_iter = 0;

   int arg_index, myid, nprocs;
   char filename[255], filename_lv[255];
   FILE *file, *file_lv;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);

   /* ntime time intervals */
   comm = MPI_COMM_WORLD;
   ntime = 2048;
   tstart = 0.0;
   tstop = log(10.) / 0.9 * tlyap;

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
            printf("Solve the Lorenz system with the following parameters\n");
            printf("  -ntime <ntime>    : set num time points (default %d)\n", ntime);
            printf("  -tstop <tstop>    : set end time in Lyapunov time (default %lf)\n", tlyap);
            printf("  -ml  <max_levels> : set max levels\n");
            printf("  -nu  <nrelax>     : set num F-C relaxations\n");
            printf("  -nu0 <nrelax>     : set num F-C relaxations on level 0\n");
            printf("  -tol <tol>        : set stopping tolerance\n");
            printf("  -cf  <cfactor>    : set coarsening factor\n");
            printf("  -mi  <max_iter>   : set max iterations\n");
            printf("  -fmg              : use FMG cycling\n");
            printf("  -rank             : rank of Delta correction (integer values in [0, 3])\n");
            printf("  -noLyap           : turn off estimation of Lyapunov vectors (static basis)\n");
            printf("  -fcfLyap          : turn on relaxation of Lyapunov vectors\n");
            printf("  -defer-lvl        : defer Delta correction to given level (default 0)\n");
            printf("  -defer-iter       : defer Delta correction until given iteration (default 0)\n");
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
         tlyap = atof(argv[arg_index++]);
         tstop = log(10.) / 0.9 * tlyap; /* the leading Lyapunov exponent is ~0.9 */
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
      else if (strcmp(argv[arg_index], "-rank") == 0)
      {
         arg_index++;
         delta_rank = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-noLyap") == 0)
      {
         arg_index++;
         lyap = 0;
      }
      else if (strcmp(argv[arg_index], "-fcfLyap") == 0)
      {
         arg_index++;
         relax_lyap = 1;
      }
      else if (strcmp(argv[arg_index], "-defer-lvl") == 0)
      {
         arg_index++;
         defer_lvl = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-defer-iter") == 0)
      {
         arg_index++;
         defer_iter = atoi(argv[arg_index++]);
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
   sprintf(filename, "%s.%05d", "ex-07.out", myid);
   sprintf(filename_lv, "%s.%05d", "ex-07-lv.out", myid);
   if (test)
   {
      file = stdout;
      file_lv = stdout;
   }
   else
   {
      file = fopen(filename, "w");
      file_lv = fopen(filename_lv, "w");
   }

   /* set up app structure */
   app = (my_App *)malloc(sizeof(my_App));
   my_App_Init(app, MPI_COMM_WORLD, tstart, tstop, ntime, delta_rank, file, file_lv);

   if (test)
   {
      braid_TestInitAccess(app, MPI_COMM_WORLD, stdout, 0., my_Init, my_Access, my_Free);
      braid_TestInitAccess(app, MPI_COMM_WORLD, stdout, 1., my_Init, my_Access, my_Free);
      braid_TestClone(app, MPI_COMM_WORLD, stdout, 0., my_Init, my_Access, my_Free, my_Clone);
      braid_TestSum(app, MPI_COMM_WORLD, stdout, 0., my_Init, my_Access, my_Free, my_Clone, my_Sum);
      braid_TestSpatialNorm(app, MPI_COMM_WORLD, stdout, 0., my_Init, my_Free, my_Clone, my_Sum, my_SpatialNorm);
      braid_TestBuf(app, MPI_COMM_WORLD, stdout, 0., my_Init, my_Free, my_Sum, my_SpatialNorm, my_BufSize, my_BufPack, my_BufUnpack);
      /* if the basic wrapper tests pass, test the routines specific to Delta correction */
      braid_TestDelta(app, MPI_COMM_WORLD, stdout, 0., 0.01, delta_rank, my_Init, my_InitBasis, my_Access, my_Free, my_Clone, my_Sum, my_BufSize, my_BufPack, my_BufUnpack, my_InnerProd, my_Step);

      MPI_Finalize();
      free(app);

      return 0;
   }

   braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app, my_Step,
              my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, my_Access,
              my_BufSize, my_BufPack, my_BufUnpack, &core);

   if (delta_rank > 0)
   {
      braid_SetDeltaCorrection(core, delta_rank, my_InitBasis, my_InnerProd);
      braid_SetDeferDelta(core, defer_lvl, defer_iter);
      braid_SetLyapunovEstimation(core, relax_lyap, lyap, relax_lyap || lyap);
   }

   braid_SetSkip(core, 0);
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

   braid_Drive(core);

   /* get average Lyapunov exponents */
   if (app->rank > 0)
   {
      double global_exps[VecSize];
      MPI_Reduce(&(app->lyap_exps), &global_exps, VecSize, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      // if (1)
      if (myid == 0)
      {  
         /* print the average Lyapunov exponents */
         printf("  Average Lyapunov Exponents: [");
         for (int i = 0; i < app->rank; i++)
         {
            global_exps[i] /= tstop;
            printf("%4f", global_exps[i]);
            if (i != app->rank - 1)
            {
               printf(", ");
            }
         }
         printf("]\n");
      }
   }

   braid_Destroy(core);

   /* reorder the output file */
   if (myid == 0)
   {
      int findex = 0;
      int npoints = (ntime + 1);
      VEC *solution;
      VEC *lyapunov;
      double *exponents;
      fclose(file);
      fclose(file_lv);
      solution = (VEC *)malloc(npoints * sizeof(VEC));
      lyapunov = (VEC *)malloc(npoints * delta_rank * sizeof(VEC));
      exponents = malloc(npoints * delta_rank * sizeof(double));
      file = fopen(filename, "r");
      file_lv = fopen(filename_lv, "r");
      for (size_t ti = 0; ti < npoints; ti++)
      {
         fscanf(file, "%d", &findex);
         for (size_t i = 0; i < VecSize; i++)
         {
            fscanf(file, "%le", &solution[findex][i]);
         }

         fscanf(file_lv, "%d", &findex);
         for (size_t j = 0; j < delta_rank; j++)
         {
            fscanf(file_lv, "%le", &exponents[findex * delta_rank + j]);
            for (size_t i = 0; i < VecSize; i++)
            {
               fscanf(file_lv, "%le", &lyapunov[findex * delta_rank + j][i]);
            }
         }
      }
      fclose(file);
      fclose(file_lv);
      file = fopen("ex-07.out", "w");
      file_lv = fopen("ex-07-lv.out", "w");
      for (size_t ti = 0; ti < npoints; ti++)
      {
         if (file)
         {
            for (size_t i = 0; i < VecSize; i++)
            {
               fprintf(file, " %.14e", solution[ti][i]);
            }
            fprintf(file, "\n");
         }

         if (file_lv)
         {
            for (size_t j = 0; j < delta_rank; j++)
            {
               fprintf(file_lv, " %.14e", exponents[ti * delta_rank + j]);
               for (size_t i = 0; i < VecSize; i++)
               {
                  fprintf(file_lv, " %.14e", lyapunov[ti * delta_rank + j][i]);
               }
            }
            fprintf(file_lv, "\n");
         }
      }
      free(solution);
      free(lyapunov);
      free(exponents);
      fclose(file);
      fclose(file_lv);
   }

   /* Finalize MPI */
   MPI_Finalize();

   return 0;
}
