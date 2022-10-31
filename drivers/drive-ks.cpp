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

//
// Driver:        drive-ks.cpp
//
// Interface:     C++
//
// Requires:      C-language and C++ support, Eigen linear algebra library:   https://gitlab.com/libeigen/eigen
//                umfpack
//
// Compile with:  make drive-ks
//
// Help with:     ./drive-ks -help
//
// Sample run:    mpirun -np 2 drive-ks
//
// Description:   solve the Kuramoto-Shivasinsky equation using 4th order finite differencing and Lobatto IIIC with optional low rank Delta correction
//
//

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string.h>
#include <vector>

#include "braid.hpp"
#include "drive-ks-lib/ks_lib.hpp"

// --------------------------------------------------------------------------
// User-defined routines and objects
// --------------------------------------------------------------------------

// Define BraidVector, can contain anything, and be named anything
// --> Put all time-dependent information here
class BraidVector
{
   public:
   VEC state;
   VEC guess;
   int stages;

   // constructors
   BraidVector(const BraidVector &other) : state(other.state), guess(other.guess), stages(other.stages) {}
   BraidVector(VEC state_, VEC guess_, int stages_) : state(state_), guess(guess_), stages(stages_) {}
   BraidVector(int nx, int stages) : state(VEC::Zero(nx)), guess(VEC::Zero(stages * nx)), stages(stages) {}

   // Deconstructor
   virtual ~BraidVector(){};
};

// Wrapper for BRAID's App object
// --> Put all time INDEPENDENT information here
class MyBraidApp : public BraidApp
{
   protected:
   // BraidApp defines tstart, tstop, ntime and comm_t

   public:
   int cfactor;
   int cfactor0;
   int newton_iters;
   int DeltaRank;
   int DeltaLevel;
   int delayDelta;
   int max_levels;
   bool cglv;
   bool useTheta;
   bool doRefine;
   bool getInit;
   bool output;
   int stages;
   std::vector<int> sclevels;
   std::vector<double> thetas;
   std::vector<double> err_ests;
   double est_norm;

   int nx;
   KSDiscretization disc;
   std::vector<KSDiscretization> DiscTable;
   VEC initial_data;

   // Constructor
   MyBraidApp(MPI_Comm comm_t_,
              int rank_,
              double tstart_,
              double tstop_,
              int ntime_,
              int target_nt_,
              int cfactor_,
              int cfactor0_,
              int DeltaLevel_,
              int DeltaRank_,
              bool cglv_,
              bool useTheta_,
              bool doRefine_,
              bool getInit_,
              bool output_,
              int newton_iters_,
              int max_levels_,
              int nx,
              double length,
              int order,
              const VEC &initial_data_);

   // We will need the MPI Rank
   int rank;

   // Deconstructor
   virtual ~MyBraidApp(){};

   bool IsCPoint(int i, int level);

   bool
   IsFPoint(int i, int level)
   {
      return !IsCPoint(i, level);
   }

   double getTheta(int level);

   int spatialLevel(int level);

   // Define all the Braid Wrapper routines
   // Note: braid_Vector == BraidVector*
   virtual int Step(braid_Vector u_,
                    braid_Vector ustop_,
                    braid_Vector fstop_,
                    BraidStepStatus &pstatus);

   virtual int Clone(braid_Vector u_,
                     braid_Vector *v_ptr);

   virtual int Init(double t,
                    braid_Vector *u_ptr);

   virtual int InitBasis(double t,
                         int index,
                         braid_Vector *u_ptr);

   virtual int Free(braid_Vector u_);

   virtual int Sum(double alpha,
                   braid_Vector x_,
                   double beta,
                   braid_Vector y_);

   virtual int SpatialNorm(braid_Vector u_,
                           double *norm_ptr);

   virtual int InnerProd(braid_Vector u_,
                         braid_Vector v_,
                         double *prod_ptr);

   virtual int BufSize(int *size_ptr,
                       BraidBufferStatus &status);

   virtual int BufPack(braid_Vector u_,
                       void *buffer,
                       BraidBufferStatus &status);

   virtual int BufUnpack(void *buffer,
                         braid_Vector *u_ptr,
                         BraidBufferStatus &status);

   virtual int Access(braid_Vector u_,
                      BraidAccessStatus &astatus);

   virtual int Sync(BraidSyncStatus &status);

   virtual int Coarsen(braid_Vector fu_,
                       braid_Vector *cu_ptr,
                       BraidCoarsenRefStatus &status);

   virtual int Refine(braid_Vector cu_,
                      braid_Vector *fu_ptr,
                      BraidCoarsenRefStatus &status);

   // not needed in this driver
   virtual int
   Residual(braid_Vector u_,
            braid_Vector r_,
            BraidStepStatus &pstatus)
   {
      return 0;
   }
};

// Braid App Constructor
MyBraidApp::MyBraidApp(MPI_Comm comm_t_,
                       int rank_,
                       double tstart_,
                       double tstop_,
                       int ntime_,
                       int target_nt,
                       int cfactor_,
                       int cfactor0_,
                       int DeltaLevel_,
                       int DeltaRank_,
                       bool cglv_,
                       bool useTheta_,
                       bool doRefine_,
                       bool getInit_,
                       bool output_,
                       int newton_iters_,
                       int max_levels_,
                       int nx,
                       double length,
                       int order,
                       const VEC &initial_data_) : BraidApp(comm_t_, tstart_, tstop_, ntime_)
{
   rank = rank_;
   cfactor = cfactor_;
   cfactor0 = cfactor0_;
   newton_iters = newton_iters_;
   DeltaRank = DeltaRank_;
   DeltaLevel = DeltaLevel_;
   delayDelta = 0;
   cglv = cglv_;
   useTheta = useTheta_;
   doRefine = doRefine_;
   getInit = getInit_;
   output = output_;
   max_levels = max_levels_;
   est_norm = -1.;
   int worldsize, time_chunk;
   MPI_Comm_size(comm_t, &worldsize);
   time_chunk = target_nt / worldsize + 1;
   err_ests.assign(time_chunk, 0.);

   switch (order)
   {
   case 2:
      stages = 2;
      break;

      // todo: implement 4th order method
      // case 4:
      //    stages = 3;
      //    break;

   default:
      throw std::invalid_argument("Order not implemented");
      break;
   }

   // fourth order
   Stencil d1 = Stencil({ 1. / 12, -2. / 3, 0., 2. / 3, -1. / 12 });
   Stencil d2 = Stencil({ -1. / 12, 4. / 3, -5. / 2, 4. / 3, -1. / 12 });
   Stencil d4 = Stencil({ -1. / 6, 2., -13. / 2, 28. / 3, -13. / 2, 2., -1. / 6 });
   disc = KSDiscretization(nx, length, d1, d2, d4);
   initial_data = initial_data_;

   // compute theta based on effective grid level
   thetas.assign(max_levels, 1.);
   if (useTheta)
   {
      double total_cf;
      for (int level = 1; level < max_levels; level++)
      {
         total_cf = cfactor0 * intpow(cfactor, level - 1);

         if (order == 2)
         {
            thetas[level] = (total_cf - sqrt(3 * total_cf * total_cf + 6) / 3) / total_cf; // theta2
         }

         if (order == 4)
         {
            double cf_4 = std::pow(total_cf, 4);
            if (cf_4 == 0.) // overflow
            {
               thetas[level] = 1 - std::sqrt(10) / 5.;
            }
            else
            {
               thetas[level] = 1. - std::sqrt(10 * cf_4 + 15) / (5 * total_cf * total_cf);
            }
         }
      }
   }
}

// Helper function to check if current point is a C point for this level
bool
MyBraidApp::IsCPoint(int i, int level)
{
   if (level == 0)
   {
      return ((i % cfactor0) == 0);
   }
   return ((i % cfactor) == 0);
}

// Helper function to compute theta for a given level
double
MyBraidApp::getTheta(int level)
{
   return thetas[level];
}

void
getGuessTheta2(VEC &guess, const VEC &u, const VEC &ustop, const KSDiscretization &disc, double dt)
{
   const int stages = 2;
   int nx = u.size();
   assert(guess.size() == stages * nx);

   guess.tail(nx) = dt * disc.f_ks(ustop);
   guess.head(nx) = 2 * (ustop - u) - guess.tail(nx);
}

VEC
theta2(const VEC &u, VEC &guess, const KSDiscretization &disc, BraidStepStatus &pstatus, double dt, double th_A, double th_B, double th_C, int newton_iters, double tol, double *err_est)
{
   const int stages = 2;
   int nx = u.size();
   assert(guess.size() == stages * nx);

   double th_Cs, a11, a12, a21, a22;
   th_Cs = 1. - th_A - th_B - th_C;
   a11 = (th_B + th_C) / 2;
   a12 = -(th_C) / 2;
   a21 = (th_A + th_B + th_C) / 2 + th_Cs;
   a22 = (th_A + th_C) / 2;

   VEC k(guess);
   VEC p(guess);
   VEC u1(u), u2(u);

   SPMAT A(stages * nx, stages * nx);
   Eigen::UmfPackLU<SPMAT> solver; // best option I have tried so far!
   SPMAT eye(nx, nx);
   eye.setIdentity();
   VEC rhs(stages * nx);
   for (int i = 0; i < newton_iters; i++)
   {
      u1 = u + a11 * k.head(nx) + a12 * k.tail(nx);
      u2 = u + a21 * k.head(nx) + a22 * k.tail(nx);

      rhs << k.head(nx) - dt * disc.f_ks(u1),
          k.tail(nx) - dt * disc.f_ks(u2);

      if (i > 0 && inf_norm(p) <= tol)
      {
         break;
      }

      // solve block system
      SPMAT B1 = dt * disc.f_ks_du(u1);
      SPMAT B2 = dt * disc.f_ks_du(u2);
      setup_nxnbmat(A,
                    { { eye - a11 * B1, -a12 * B1 },
                      { -a21 * B2, eye - a22 * B2 } },
                    nx);

      solver.compute(A);
      p = solver.solve(rhs);
      k -= p;
   }
   // propagate tangent vectors
   int delta_rank;
   pstatus.GetDeltaRank(&delta_rank);
   if (delta_rank > 0)
   {
      // here we propagate tangent vectors implicitly without forming the full lin. of Phi.
      BraidVector* col;
      VEC tmp(stages * nx);
      SPMAT K1 = dt * disc.f_ks_du(u1);
      SPMAT K2 = dt * disc.f_ks_du(u2);

      setup_nxnbmat(A, { { eye - a11 * K1, -a12 * K1 }, { -a21 * K2, eye - a22 * K2 } }, nx);
      solver.compute(A); // does the factorization

      // solve one column at a time
      for (Index j = 0; j < delta_rank; j++)
      {
         pstatus.GetBasisVec(((braid_Vector*)&col), j);

         rhs << K1 * col->state, K2 * col->state;
         tmp = solver.solve(rhs);
         col->state += (tmp.head(nx) + tmp.tail(nx)) / 2.;
      }
   }

   // std::cout << "guess accuracy: " << (k - guess).norm() << '\n';
   guess = k;
   u2 = u + (k.head(nx) + k.tail(nx)) / 2;
   if (err_est)
   {
      // use embedded 1st order method to estimate relative local discretization error
      u1 = u + k.head(nx);
      *err_est = (u1 - u2).lpNorm<Eigen::Infinity>() / u2.lpNorm<Eigen::Infinity>();
   }
   return u2;
}

int
MyBraidApp::Step(braid_Vector u_,
                 braid_Vector ustop_,
                 braid_Vector fstop_,
                 BraidStepStatus &pstatus)
{
   BraidVector *u = (BraidVector *)u_;
   BraidVector *ustop = (BraidVector *)ustop_;

   int level, lvl_eff, nlevels, nrefine, T_index;
   double tstart, tstop, err_est;
   int done, calling_func, mgrit_iter;

   pstatus.GetLevel(&level);
   pstatus.GetTstartTstop(&tstart, &tstop);
   pstatus.GetNLevels(&nlevels);
   pstatus.GetNRefine(&nrefine);
   pstatus.GetDone(&done);
   pstatus.GetCallingFunction(&calling_func);
   pstatus.GetIter(&mgrit_iter);
   pstatus.GetTIndex(&T_index);

   // get effective level:
   lvl_eff = level;
   if (doRefine)
   {
      lvl_eff = (level - nrefine) + (max_levels - 2);
   }
   double dt = tstop - tstart;

   bool useGuess, toCPoint, coarseGrid;
   toCPoint = IsCPoint(T_index + 1, lvl_eff);
   coarseGrid = (level == nlevels - 1);

   // the initial guess is only valid for f-points on fine-grids
   useGuess = (coarseGrid || !toCPoint);

   // We want the coarse grid to be cheap, so we save the best
   // initial guess for the coarse grid at C-points
   if (useGuess)
   {
      u->guess = ustop->guess;
   }
   else if (stages == 2)
   {
      // this is exact on the fine-grid, a lot closer on coarse grids
      getGuessTheta2(u->guess, u->state, ustop->state, disc, dt);
   }
   else
   {
      u->guess = VEC::Zero(stages * nx);
   }

   // tolerance for Newton's method
   double tol = std::sqrt(disc.nx) * 1e-10;
   // dynamic solver tolerance
   // double tol, tight{1e-13/std::sqrt(disc.nx)}, loose{1e-11/std::sqrt(disc.nx)};
   // if (level > 0)
   // {
   //    tol = tight;
   // }
   // else
   // {
   //    pstatus.GetSpatialAccuracy(loose, tight, &tol);
   // }
   // if (nlevels == 1)
   // {
   //    tol = tight;
   // }

   int iters = newton_iters;

   // second order theta method
   if (stages == 2)
   {
      if (lvl_eff == 0 || !useTheta)
      {
         u->state = theta2(u->state, u->guess, disc, pstatus, dt, 0., 0., 1., iters, tol, &err_est);
      }
      else
      {
         double theta = getTheta(lvl_eff);
         double cf = cfactor0 * intpow(cfactor, lvl_eff - 1);
         u->state = theta2(u->state, u->guess, disc, pstatus, dt, theta, theta, 2. / 3 + 1 / (3 * cf * cf) - theta, iters, tol, &err_est);
      }
   }
   else
   {
      // fourth order theta method **not implemented**
      // TODO: Get this working
      if (lvl_eff == 0 || !useTheta)
      {
         u->state = theta4(u->state, u->guess, disc, dt, 0., 0., 1., NULL, iters, tol);
      }
      else
      {
         double cf = cfactor0 * intpow(cfactor, lvl_eff - 1);
         double theta = getTheta(lvl_eff);
         double cf_4 = cf * cf * cf * cf;
         bool overflow = (cf_4 == 0.);
         double th_C = 7. / 10 - theta;
         if (!overflow)
         {
            th_C += 3. / (10 * cf_4);
         }
         u->state = theta4(u->state, u->guess, disc, dt, theta, theta, th_C, NULL, iters, tol);
      }
   }

   // never overwrite initial guess at C-points
   if (!useGuess)
   {
      u->guess = ustop->guess; // the old value
   }

   // refinement/error estimate control
   if (level > 0 || !doRefine)
   {
      // no refinement, or error estimate
      pstatus.SetRFactor(1);
      return 0;
   }

   int tu, tl, T_irel;
   pstatus.GetTIUL(&tu, &tl, level);
   T_irel = T_index;
   if (rank > 0)
   {
      T_irel = T_index - tl + 1;
   }

   err_ests[T_irel] = err_est;

   // refine based on comparison of the residual to the local discretization error
   int nrequest = -2;
   double *rnorm = (double *)(malloc(4 * sizeof(double)));
   pstatus.GetRNorms(&nrequest, rnorm);
   double last_norm = rnorm[1];
   if (rnorm[1] < 0)
   {
      last_norm = rnorm[0];
   }

   // if rnorm is -1, we don't have residual yet; if max estimate is -1, we haven't called sync yet
   if (last_norm > 0. && est_norm > 0. && last_norm <= est_norm)
   {
      if (lvl_eff == 1)
      {
         pstatus.SetRFactor(cfactor0);
      }
      else
      {
         pstatus.SetRFactor(cfactor);
      }
   }
   free(rnorm);
   return 0;
}

int
MyBraidApp::Init(double t,
                 braid_Vector *u_ptr)
{
   // initial the state and intial guess vectors for Newton's method
   BraidVector *u = new BraidVector(disc.nx, stages);
   u->state = initial_data;

   *u_ptr = (braid_Vector)u;
   return 0;
}

int
MyBraidApp::InitBasis(double t,
                      int index,
                      braid_Vector *u_ptr)
{
   // No need to store initial guesses for the lyapunov vectors
   BraidVector *u = new BraidVector(disc.nx, 0);
   u->state = FourierMode(index, disc.nx, disc.len);

   *u_ptr = (braid_Vector)u;
   return 0;
}

int
MyBraidApp::Clone(braid_Vector u_,
                  braid_Vector *v_ptr)
{
   BraidVector *u = (BraidVector *)u_;
   BraidVector *v = new BraidVector(*u);
   *v_ptr = (braid_Vector)v;

   return 0;
}

int
MyBraidApp::Free(braid_Vector u_)
{
   BraidVector *u = (BraidVector *)u_;
   delete u;
   return 0;
}

int
MyBraidApp::Sum(double alpha,
                braid_Vector x_,
                double beta,
                braid_Vector y_)
{
   BraidVector *x = (BraidVector *)x_;
   BraidVector *y = (BraidVector *)y_;
   (y->state) = alpha * (x->state) + beta * (y->state);
   // (y->guess) = alpha * (x->guess) + beta * (y->guess);
   return 0;
}

int
MyBraidApp::SpatialNorm(braid_Vector u_,
                        double *norm_ptr)
{
   BraidVector *u = (BraidVector *)u_;
   *norm_ptr = u->state.norm() / std::sqrt(disc.nx); // normalized like l2 norm
   return 0;
}

int
MyBraidApp::InnerProd(braid_Vector u_,
                      braid_Vector v_,
                      double *prod_ptr)
{
   BraidVector *u = (BraidVector *)u_;
   BraidVector *v = (BraidVector *)v_;
   *prod_ptr = (u->state).dot(v->state);
   return 0;
}

int
MyBraidApp::BufSize(int *size_ptr,
                    BraidBufferStatus &status)
{
   int headersize = sizeof(int);
   *size_ptr = headersize + (1 + stages) * disc.nx * sizeof(double);

   // tell Braid the size of the basis vectors (which don't have initial guesses)
   status.SetBasisSize(headersize + disc.nx * sizeof(double));

   return 0;
}

int
MyBraidApp::BufPack(braid_Vector u_,
                    void *buffer,
                    BraidBufferStatus &status)
{
   BraidVector *u = (BraidVector *)u_;

   // need to pack the number of guess vectors, so we can know how many to unpack later
   int *header = (int *)buffer;
   int headersize = sizeof(int);

   *header = u->stages;
   header += 1;

   double *dbuffer = (double *)header;
   size_t bf_size = 0;
   bf_pack_help(dbuffer, u->state, disc.nx, bf_size);
   if (u->stages > 0)
   {
      bf_pack_help(dbuffer, u->guess, u->stages * disc.nx, bf_size);
   }
   status.SetSize(bf_size * sizeof(double) + headersize);
   return 0;
}

int
MyBraidApp::BufUnpack(void *buffer,
                      braid_Vector *u_ptr,
                      BraidBufferStatus &status)
{


   int *header = (int *)buffer;

   BraidVector *u = new BraidVector(disc.nx, *header);
   header += 1;
   double *dbuffer = (double *)header;

   size_t bf_size = 0;
   bf_unpack_help(dbuffer, u->state, disc.nx, bf_size);
   if (u->stages > 0)
   {
      bf_unpack_help(dbuffer, u->guess, u->stages * disc.nx, bf_size);
   }

   *u_ptr = (braid_Vector)u;

   return 0;
}

int
MyBraidApp::Coarsen(braid_Vector fu_,
                    braid_Vector *cu_ptr,
                    BraidCoarsenRefStatus &status)
{
   Clone(fu_, cu_ptr);
   return 0;
}

int
MyBraidApp::Refine(braid_Vector cu_,
                   braid_Vector *fu_ptr,
                   BraidCoarsenRefStatus &status)
{
   Clone(cu_, fu_ptr);
   return 0;
}

int
MyBraidApp::Access(braid_Vector u_,
                   BraidAccessStatus &astatus)
{
   char filename[255];
   char lv_fname[255];
   std::ofstream file;
   BraidVector *u = (BraidVector *)u_;

   // Extract information from astatus
   int done, level, iter, index, nt, delta_rank, testing;
   double t;
   astatus.GetTILD(&t, &iter, &level, &done);
   astatus.GetNTPoints(&nt);
   astatus.GetTIndex(&index);
   astatus.GetDeltaRank(&delta_rank);
   astatus.GetWrapperTest(&testing);

   // Print information to file
   if (getInit && done && level == 0 && index == nt)
   {
      sprintf(filename, "drive-ks-init-nx%d.out", disc.nx);
      file.open(filename);
      pack_array(file, u->state);
      file.close();
   }

   if (testing || (output && done && level == 0 && IsCPoint(index, level)))
   {
      sprintf(filename, "%s.%04d", "drive-ks.out", index);
      file.open(filename);
      pack_array(file, u->state);
      file.close();

      if (delta_rank > 0)
      {
         MAT Psi = MAT(disc.nx, delta_rank);
         BraidVector *col;

         for (int i = 0; i < delta_rank; i++)
         {
            astatus.GetBasisVec((braid_Vector*)&col, i);
            Psi.col(i) = col->state;
         }

         sprintf(lv_fname, "%s.%04d", "drive-ks-lv.out", index);
         file.open(lv_fname);
         pack_darray(file, Psi);
         file.close();
      }
   }

   return 0;
}

int
MyBraidApp::Sync(BraidSyncStatus &status)
{
   /* Sync can be called from two places, at the top of each Braid Cycle or,
    * from inside of Refine.  Here, we don't care, but show you this for your
    * reference. This is how you detect your calling function */
   braid_Int calling_fcn;
   status.GetCallingFunction(&calling_fcn);
   if (calling_fcn == braid_ASCaller_FRefine_AfterInitHier)
   {
      est_norm = -1;
      return 0;
   }

   if (calling_fcn == braid_ASCaller_Drive_TopCycle)
   {
      /*
       * Find the global max of the error estimate
       *
       * Note: if there is spatial parallelism, then you would have to take a
       * max over the spatial communicator, and then over the temporal
       * communicator.
       *
       * */
      double my_sum = 0.;
      for (auto &&est : err_ests)
      {
         my_sum += est * est;
      }

      MPI_Allreduce(&my_sum, &(est_norm), 1, MPI_DOUBLE, MPI_SUM, comm_t);
      est_norm = std::sqrt(est_norm);

      // if (rank == 0)
      // {
      //    printf("  Braid: || est || =  %1.5e\n", est_norm);
      // }
   }
   return 0;
}

// --------------------------------------------------------------------------
// Main driver
// --------------------------------------------------------------------------

int
read_init(std::string fname, VEC &out, int nx)
{
   out.resize(nx);
   std::ifstream inf;
   std::string entry;

   inf.open(fname);
   if (!inf)
   {
      std::cout << "!!!init file not found!!!\n";
      return 0;
   }

   for (int i = 0; i < nx; i++)
   {
      getline(inf, entry, ',');
      out[i] = std::stof(entry);
   }
   return 1;
}

int
del_extra(int nt, int cfactor, std::string fname)
{
   std::ifstream inf;
   char ifname[255];
   // loop over all possible files with that name
   for (int i = 0; i <= nt; i += cfactor)
   {
      // check if they exist
      sprintf(ifname, "%s.%s.%04d", fname.c_str(), "out", i);
      inf.open(ifname);
      if (inf)
      {
         // delete if they do
         inf.close();
         remove(ifname);
      }
   }
   return 0;
}

int
collate_files(int nt, int cfactor, std::string fname)
{
   // collate files:
   std::ifstream inf;
   std::ofstream of;
   char ifname[255];
   char ofname[255];

   // state vectors:
   sprintf(ofname, "%s.out", fname.c_str());
   of.open(ofname);
   for (int i = 0; i <= nt; i += cfactor)
   {
      sprintf(ifname, "%s.%s.%04d", fname.c_str(), "out", i);
      inf.open(ifname);
      // if we can't open one file, delete all extra files
      if (!inf)
      {
         std::cout << "!!!error collating files- deleting extra!!!\n";
         del_extra(nt, cfactor, fname);
         return 1;
      }
      of << inf.rdbuf();
      inf.close();
      // delete the extra files
      remove(ifname);
   }
   of.close();
   return 0;
}

int
main(int argc, char *argv[])
{
   double tstart, tstop;
   int rank;

   int nt = 2048;
   double Tf_lyap = 4;
   int max_levels = 4;
   int nrelax = 1;
   int nrelax0 = 1;
   double tol = 1e-8;
   int cfactor = 4;
   int cf0 = 4;
   int max_iter = 25;
   int newton_iters = 3;
   bool useFMG = false;
   int DeltaLvl = 0;
   int DeltaRank = 0;
   bool useTheta = false;
   bool doRefine = false;
   bool wrapperTests = false;
   bool getInit = false;
   bool output = false;
   std::vector<int> sclevels;
   sclevels.reserve(100);
   int ord = 2;
   bool cglv = false;
   double weight = 1.;
   bool rcg = false;

   // KS parameters
   double len = 64;
   int nx = 64;

   int arg_index;

   // Initialize MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if (rank == 0)
   {
      int world_size;
      MPI_Comm_size(MPI_COMM_WORLD, &world_size);
      std::cout << "nprocs: " << world_size << '\n';
   }

   // Parse command line
   arg_index = 1;
   while (arg_index < argc)
   {
      if (strcmp(argv[arg_index], "-help") == 0)
      {
         if (rank == 0)
         {
            printf("\n");
            printf("  -nt         : set num time points (default %d)\n", nt);
            printf("  -nx         : set num space points (default %d)\n", nx);
            printf("  -len        : set length of spatial domain (default %f)\n", len);
            printf("  -ord        : set order of implicit solver (2 or 4, default %d)\n", ord);
            printf("  -tf         : set end time, in Lyapunov time (default %lf)\n", Tf_lyap);
            printf("  -ml         : set max levels\n");
            printf("  -nu         : set num F-C relaxations\n");
            printf("  -nu0        : set num F-C relaxations on level 0\n");
            printf("  -tol        : set stopping tolerance\n");
            printf("  -cf         : set coarsening factor\n");
            printf("  -cf0        : set coarsening factor on level 0\n");
            printf("  -mi         : set max iterations\n");
            printf("  -niters     : set number of newton iters for spatial solver\n");
            printf("  -fmg        : use F cycles\n");
            printf("  -Delta      : use Delta correction\n");
            printf("  -Deltalvl   : Delta correction is deferred until this level\n");
            printf("  -cglv       : Propagate Lyapunov Vectors on the coarsest grid\n");
            printf("  -rank       : set rank of delta correction\n");
            printf("  -theta      : use theta method\n");
            printf("  -refine     : use global refinement to approximate true FMG cycle\n");
            printf("  -sclevels   : levels on which to coarsen in space (comma separated, e.g. '1,3,5', no spaces\n");
            printf("  -out        : write output to file (for visualization)\n");
            printf("  -test       : run wrapper tests\n");
            printf("  -getInit    : write solution at final time to file\n");
            printf("  -weight     : weight for weighted relaxation\n");
            printf("  -rcg        : use relaxation only on coarsest grid\n");
            printf("\n");
         }
         exit(0);
      }
      else if (strcmp(argv[arg_index], "-nt") == 0)
      {
         arg_index++;
         nt = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-nx") == 0)
      {
         arg_index++;
         nx = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-len") == 0)
      {
         arg_index++;
         len = atof(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-ord") == 0)
      {
         arg_index++;
         ord = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-tf") == 0)
      {
         arg_index++;
         Tf_lyap = atof(argv[arg_index++]);
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
      else if (strcmp(argv[arg_index], "-cf0") == 0)
      {
         arg_index++;
         cf0 = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-mi") == 0)
      {
         arg_index++;
         max_iter = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-niters") == 0)
      {
         arg_index++;
         newton_iters = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-fmg") == 0)
      {
         arg_index++;
         useFMG = true;
      }
      else if (strcmp(argv[arg_index], "-Deltalvl") == 0)
      {
         arg_index++;
         DeltaLvl = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-cglv") == 0)
      {
         arg_index++;
         cglv = true;
      }
      else if (strcmp(argv[arg_index], "-rank") == 0)
      {
         arg_index++;
         DeltaRank = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-theta") == 0)
      {
         arg_index++;
         useTheta = true;
      }
      else if (strcmp(argv[arg_index], "-refine") == 0)
      {
         arg_index++;
         doRefine = true;
      }
      else if (strcmp(argv[arg_index], "-sclevels") == 0)
      {
         arg_index++;
         // create a string stream and read the values from it
         std::stringstream s_stream(argv[arg_index++]);
         while (s_stream.good())
         {
            std::string substr;
            std::getline(s_stream, substr, ',');
            sclevels.push_back(std::stoi(substr));
         }
      }
      else if (strcmp(argv[arg_index], "-out") == 0)
      {
         arg_index++;
         output = true;
      }
      else if (strcmp(argv[arg_index], "-getInit") == 0)
      {
         arg_index++;
         getInit = true;
      }
      else if (strcmp(argv[arg_index], "-test") == 0)
      {
         arg_index++;
         wrapperTests = true;
      }
      else if (strcmp(argv[arg_index], "-weight") == 0)
      {
         arg_index++;
         weight = atof(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-rcg") == 0)
      {
         arg_index++;
         rcg = true;
      }
      else
      {
         arg_index++;
         /*break;*/
      }
   }

   tstart = 0.0;
   tstop = Tf_lyap * T_lyap;

   if (!sclevels.empty())
   {
      std::cout << "using spatial coarsening on levels: ";
      for (auto &&i : sclevels)
      {
         std::cout << i << ", ";
      }
      std::cout << '\n';
   }
   if (rank == 0 && DeltaLvl > 0)
   {
      std::cout << "Using Delta correction\n";
   }
   if (rank == 0 && useTheta)
   {
      std::cout << "Using theta method\n";
   }

   // make sure we do enough iterations when sequential time-stepping
   if (max_levels == 1)
   {
      newton_iters = 10;
   }

   // calculate min coarse grid size:
   double max_dt = 1.5;
   if (useTheta)
   {
      max_dt = 3.;
   }

   int min_cg = std::ceil(tstop / max_dt);
   int cnt = nt / cf0;
   int actual_lvls = max_levels;

   // figure out how many levels we will actually have
   for (int i = 1; i < max_levels; i++)
   {
      cnt /= cfactor;
      if (cnt < min_cg)
      {
         actual_lvls = i + 1;
         break;
      }
   }
   cnt *= cfactor * cfactor;

   int nt_init = nt;
   if (doRefine)
   {
      nt_init = cnt;
   }

   // get initial data
   VEC u0;
   u0 = FourierMode(1, nx, len) + smoothed_noise(nx, 4);
   if (!getInit)
   {
      using namespace Eigen;
      int init_nx = 512;
      int cfx = init_nx / nx;
      VEC u0fine;
      if (read_init("drive-ks-init-nx512.out", u0fine, init_nx))
      {
         u0 = u0fine(seq(0, last, cfx));
      }
   }

   // set up app structure
   MyBraidApp app(MPI_COMM_WORLD,
                  rank, tstart, tstop, nt_init, nt, cfactor, cf0,
                  DeltaLvl, DeltaRank, cglv,
                  useTheta, doRefine, getInit, output,
                  newton_iters, actual_lvls,
                  nx, len, ord, u0);

   // wrapper tests
   if (wrapperTests)
   {
      if (rank != 0)
      {
         // Clean up
         MPI_Finalize();
         return 0;
      }
      BraidUtil Util = BraidUtil();
      // FILE *ftest = fopen("KSwrapperTests.txt", "w");
      Util.TestInitAccess(&app, MPI_COMM_WORLD, stdout, 0.);
      Util.TestInitAccess(&app, MPI_COMM_WORLD, stdout, 1.);
      Util.TestClone(&app, MPI_COMM_WORLD, stdout, 0.);
      Util.TestSpatialNorm(&app, MPI_COMM_WORLD, stdout, 0.);
      Util.TestBuf(&app, MPI_COMM_WORLD, stdout, 0.);

      if (DeltaRank > 0)
      {
         Util.TestDelta(&app, MPI_COMM_WORLD, stdout, 0., 0.01, DeltaRank);
      }
      else
      {
         Util.TestDelta(&app, MPI_COMM_WORLD, stdout, 0., 0.01, 9);
      }

      // Clean up
      MPI_Finalize();
      return 0;
   }

   // Initialize Braid Core Object and set some solver options
   BraidCore core(MPI_COMM_WORLD, &app);
   if (DeltaRank > 0)
   {
      core.SetDeltaCorrection(DeltaRank);
      if (cglv)
      {
         core.SetLyapunovEstimation();
      }
   }
   if (useFMG)
   {
      core.SetFMG();
      core.SetNFMG(1);
   }
   core.SetPrintLevel(2);
   core.SetMinCoarse(min_cg);
   core.SetMaxLevels(actual_lvls);
   if (doRefine)
   {
      core.SetMaxLevels(2);
      core.SetIncrMaxLevels();
   }
   core.SetMaxIter(max_iter);
   core.SetAbsTol(tol);
   core.SetCFactor(-1, cfactor);
   core.SetCFactor(0, cf0);
   core.SetNRelax(-1, nrelax);
   core.SetNRelax(0, nrelax0);
   core.SetRefine(doRefine);
   core.SetMaxRefinements(max_levels - 2);
   core.SetTPointsCutoff(nt);
   core.SetSkip(1);
   core.SetStorage(0);
   core.SetTemporalNorm(2);
   core.SetAccessLevel(output || getInit);
   core.SetSync();
   core.SetCRelaxWt(-1, weight);
   core.SetRelaxOnlyCG(rcg);
   core.SetNRelax(actual_lvls - 1, 1);

   // Run Simulation
   core.Drive();

   if (rank == 0 && (output))
   {
      collate_files(nt, cf0, "drive-ks");
      collate_files(nt, cf0, "drive-ks-lv");
   }

   // Clean up
   MPI_Finalize();

   return 0;
}
