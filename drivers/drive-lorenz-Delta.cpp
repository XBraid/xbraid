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
// Driver:       drive-lorenz-Delta.cpp
//
// Interface:     C++
//
// Requires:      C-language and C++ support, Eigen linear algebra library:
//                https://gitlab.com/libeigen/eigen
//                Get the library file "Eigen" from the link above and copy it into
//                drivers/drive-lorenz-Delta-lib/
//
// Compile with:  make drive-lorenz-Delta
//
// Help with:     ./drive-lorenz-Delta -help
//
// Sample run:    mpirun -np 2 drive-lorenz-Delta
//
// Description:   solve the Lorenz system using forward Euler and Delta correction
//
//

#include <iostream>
#include <cmath>
#include <vector>
#include <string.h>
#include <fstream>

#include "braid.hpp"
#include "drive-lorenz-Delta-lib/lorenz_lib.hpp"

// --------------------------------------------------------------------------
// User-defined routines and objects
// --------------------------------------------------------------------------

// Define BraidVector, can contain anything, and be named anything
// --> Put all time-dependent information here
class BraidVector
{
public:
   // Each vector holds the state vector at a particular time
   VEC state;
   VEC guess;
   VEC action;
   MAT Delta;

   // Construct a BraidVector for a given vector of doubles
   BraidVector(VEC state_, VEC guess_, VEC prev_c_point_, MAT Delta_) : state(state_), guess(guess_), action(prev_c_point_), Delta(Delta_) {}
   BraidVector() : state(VEC::Zero()), guess(VEC::Zero()), action(VEC::Zero()), Delta(MAT::Identity()) {}

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
   int newton_iters;  // only used if useTheta
   double newton_tol; // only used if useTheta
   int DeltaLevel;
   bool useDelta;
   bool useTheta;
   std::vector<double> thetas;

   // Constructor
   MyBraidApp(MPI_Comm comm_t_, int rank_, double tstart_, double tstop_, int ntime_, int cfactor_, int cfactor0_, bool useDelta_, int DeltaLevel_, bool useTheta_, int newton_iters_, double newton_tol_, int max_levels);

   // We will need the MPI Rank
   int rank;

   // Deconstructor
   virtual ~MyBraidApp(){};

   int IsCPoint(int i, int level);

   double getTheta(int level);

   VEC baseStep(VEC u,
                VEC &uguess,
                double dt,
                int level,
                int nlevels,
                MAT *P_tan_ptr = nullptr);

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

   virtual int Free(braid_Vector u_);

   virtual int Sum(double alpha,
                   braid_Vector x_,
                   double beta,
                   braid_Vector y_);

   virtual int SpatialNorm(braid_Vector u_,
                           double *norm_ptr);

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

   virtual int Residual(braid_Vector u_,
                        braid_Vector f_,
                        braid_Vector r_,
                        BraidStepStatus &pstatus);

   // Not needed in this driver
   virtual int Coarsen(braid_Vector fu_,
                       braid_Vector *cu_ptr,
                       BraidCoarsenRefStatus &status) { return 0; }

   // Not needed in this driver
   virtual int Refine(braid_Vector cu_,
                      braid_Vector *fu_ptr,
                      BraidCoarsenRefStatus &status) { return 0; }
};

// Braid App Constructor
MyBraidApp::MyBraidApp(MPI_Comm comm_t_, int rank_, double tstart_, double tstop_, int ntime_, int cfactor_, int cfactor0_, bool useDelta_, int DeltaLevel_, bool useTheta_, int newton_iters_, double newton_tol_, int max_levels)
    : BraidApp(comm_t_, tstart_, tstop_, ntime_)
{
   rank = rank_;
   cfactor = cfactor_;
   cfactor0 = cfactor0_;
   newton_iters = newton_iters_;
   newton_tol = newton_tol_;
   DeltaLevel = DeltaLevel_;
   useDelta = useDelta_;
   useTheta = useTheta_;
   thetas.assign(max_levels, 1.);
   if (useTheta)
   {
      double total_cf;
      double &cf_pow = total_cf;
      for (int level = 1; level < max_levels; level++)
      {
         total_cf = cfactor0 * intpow(cfactor, level - 1);
         // thetas[level] = (1 + total_cf) / (2 * total_cf); // asymptotic values for forward Euler
         cf_pow = pow(total_cf, 4);
         if (cf_pow == 0) // overflow
         {
            thetas[level] = 1. / 5.;
         }
         else
         {
            thetas[level] = (4 + cf_pow) / (double)(5 * cf_pow); // asymptotic values for rk4
         }
         // thetas[level] = (total_cf - sqrt(3*total_cf*total_cf - 3)/3)/total_cf; // theta2
      }
   }
}

// Helper function to check if current point is a C point for this level
int MyBraidApp::IsCPoint(int i, int level)
{
   if (level == 0)
   {
      return ((i % cfactor0) == 0);
   }
   // else:
   return ((i % cfactor) == 0);
}

// Helper function to compute theta for a given level
double MyBraidApp::getTheta(int level)
{
   return thetas[level];
}

VEC MyBraidApp::baseStep(const VEC u, VEC &uguess, double dt, int level, int nlevels, MAT *P_tan_ptr)
{
   // this should return Phi(u) and set u->guess

   // first order theta method
   // if (level == 0 || !useTheta)
   // {
   //    // forward euler
   //    if (P_tan_ptr)
   //    {
   //       *P_tan_ptr = euler_du(u, dt);
   //    }
   //    return euler(u, dt);
   // }
   // double theta = getTheta(level);
   // if (uguess.isZero())
   // {
   //    uguess = u; // much safer initial guess
   // }
   // VEC ustop = theta1(u, uguess, dt, theta, P_tan_ptr, newton_iters, newton_tol);
   // uguess = ustop;
   // return ustop;
   // no initial guess:
   // return theta1(u, u, dt, theta, P_tan_ptr, newton_iters, 1e-14);

   // second order theta method
   // double theta = getTheta(level);
   // if (level == 0 || !useTheta)
   // {
   //    return crank_nicolson(u, uguess, dt, P_tan_ptr, newton_iters);
   //    // return theta2(u, uguess, dt, 1., 0., 0., P_tan_ptr, newton_iters);
   // }
   // double cf = intpow(cfactor, level);
   // // A-stable
   // return theta2(u, uguess, dt, theta, theta, 0.5 - theta, P_tan_ptr, newton_iters);
   // L-stable
   // return theta2(u, uguess, dt, theta, theta, 2./3. - 1/(6*cf*cf) - theta, P_tan_ptr, newton_iters);

   // fourth order theta method
   if (level == 0 || !useTheta)
   {
      return rk4(u, dt, P_tan_ptr);
   }
   double theta = getTheta(level);
   // else:
   int iters = newton_iters;
   if (level == nlevels - 1)
   {
      // limit newton iterations on coarsest grid!
      iters = std::min(newton_iters, 2);
   }

   // no initial guess
   // uguess = u;
   return theta4(u, uguess, dt, theta, P_tan_ptr, iters, 1e-10);
}

//
int MyBraidApp::Step(braid_Vector u_,
                     braid_Vector ustop_,
                     braid_Vector fstop_,
                     BraidStepStatus &pstatus)
{

   BraidVector *u = (BraidVector *)u_;
   BraidVector *f = (BraidVector *)fstop_;
   BraidVector *ustop = (BraidVector *)ustop_;

   double tstart; // current time
   double tstop;  // evolve to this time
   int level, nlevels, T_index, calling_fnc;

   pstatus.GetTstartTstop(&tstart, &tstop);
   pstatus.GetLevel(&level);
   pstatus.GetNLevels(&nlevels);
   pstatus.GetTIndex(&T_index); // this is the index of tstart
   pstatus.GetCallingFunction(&calling_fnc);

   double dt = tstop - tstart;

   // no refinement
   pstatus.SetRFactor(1);

   bool DeltaCorrect = (useDelta) && (level > DeltaLevel);
   bool computeDeltas = (useDelta) && (level >= DeltaLevel);
   computeDeltas = (computeDeltas) && (calling_fnc == braid_ASCaller_FRestrict); // only compute Deltas when in FRestrict

   // std::cout << "Stp called @ " << T_index << " on level " << level << '\n';
   // std::cout << "f is null: " << f_is_null << '\n';
   // if (level > 0)
   // {
   //    std::cout << level << ": ";
   //    if (calling_fnc == braid_ASCaller_FRestrict || calling_fnc == braid_ASCaller_FInterp || calling_fnc == braid_ASCaller_FCRelax)
   //    {
   //       std::cout << "Stp: ";
   //    }
   //    else if (calling_fnc == braid_ASCaller_FASResidual || calling_fnc == braid_ASCaller_Residual)
   //    {
   //       std::cout << "Res: ";
   //    }
   //    else
   //    {
   //       std::cout << calling_fnc << ": ";
   //    }
   // }

   VEC guess;
   guess << ustop->guess;

   VEC utmp;
   MAT Ptmp;

   // take step and compute linear tangent propagators
   if (computeDeltas)
   {
      // store state and compute linear tangent propagator
      utmp = baseStep(u->state, guess, dt, level, nlevels, &Ptmp);
      if (f && DeltaCorrect) // use Delta correction
      {
         Ptmp += f->Delta;
      }

      if (IsCPoint(T_index, level))
      {
         // Need to store the value at the previous c-point for tau correction later
         u->action = u->state;

         // we implicitly set u->Delta = I at each C-point
         u->Delta = Ptmp;
      }
      else
      {
         // u->action already stores the correct value
         // u->Delta stores the previous tangent matrix
         u->Delta = Ptmp * u->Delta;
      }
   }
   else
   {
      utmp = baseStep(u->state, guess, dt, level, nlevels);
   }

   // on fine grids, only overwrite u->guess if ustop is an f-point:
   if (!IsCPoint(T_index + 1, level) || level == nlevels - 1)
   {
      u->guess = guess;          // the new value
   }
   else
   {
      u->guess = ustop->guess;   // the old value
   }

   if (f)
   {
      utmp += f->state;
      if (DeltaCorrect)
      {
         utmp += f->Delta * u->state - f->action;
      }
   }

   u->state = utmp;
   return 0;
}

int MyBraidApp::Residual(braid_Vector u_,
                         braid_Vector f_,
                         braid_Vector r_,
                         BraidStepStatus &pstatus)
{
   BraidVector *u = (BraidVector *)u_;
   BraidVector *f = (BraidVector *)f_;
   BraidVector *r = (BraidVector *)r_;

   double tstart; // current time
   double tstop;  // evolve to this time
   int level, nlevels, T_index, calling_fnc;

   pstatus.GetTstartTstop(&tstart, &tstop);
   pstatus.GetLevel(&level);
   pstatus.GetNLevels(&nlevels);
   pstatus.GetTIndex(&T_index);
   pstatus.GetCallingFunction(&calling_fnc);

   bool up = (calling_fnc == braid_ASCaller_Residual);
   double dt = tstop - tstart;

   // get initial guess
   VEC guess = u->guess;
   VEC utmp;

   if (!useDelta || level < DeltaLevel || (level == DeltaLevel && up))
   {
      // default residual (no linearizations computed)
      utmp = baseStep(r->state, guess, dt, level, nlevels);
      if (f)
      {
         utmp += f->state;
      }
      r->state = u->state - utmp;
      return 0;
   }
   // else we need to treat the Delta correction specially:

   MAT Ptmp;
   utmp = baseStep(r->state, guess, dt, level, nlevels, &Ptmp);

   if (up)
   { // this is called on the coarse grid right after restriction (no tau correction available)
      r->Delta = -Ptmp;
      r->action = r->Delta * r->state;
      r->state = u->state - utmp;
      return 0;
   }
   // else this is called on the fine grid right after F-relax

   if (f)
   { // do delta correction and tau correction
      Ptmp += f->Delta;
      utmp += f->Delta * r->state + (f->state - f->action);
   }

   r->Delta = -Ptmp * r->Delta;      // -D \Phi^m
   r->action = r->Delta * r->action; // -[D Phi^m]u_{i-m}
   r->state = u->state - utmp;       // u_i - Phi^m(u_{i-m})
   return 0;
}

int MyBraidApp::Init(double t,
                     braid_Vector *u_ptr)
{
   // this empty constructor should take care of everything
   BraidVector *u = new BraidVector();
   // set initial condition
   if (t == tstart)
   {
      // u->state << -1.8430428, -0.07036326, 23.15614636; // some point near the attractor
      u->state << -5.789, -9.434, 15.962; // some other point near the attractor
      u->guess << u->state;
      u->action << u->state;
   }

   *u_ptr = (braid_Vector)u;
   return 0;
}

int MyBraidApp::Clone(braid_Vector u_,
                      braid_Vector *v_ptr)
{
   BraidVector *u = (BraidVector *)u_;
   BraidVector *v = new BraidVector(u->state, u->guess, u->action, u->Delta);
   *v_ptr = (braid_Vector)v;

   return 0;
}

int MyBraidApp::Free(braid_Vector u_)
{
   BraidVector *u = (BraidVector *)u_;
   delete u;
   return 0;
}

int MyBraidApp::Sum(double alpha,
                    braid_Vector x_,
                    double beta,
                    braid_Vector y_)
{
   BraidVector *x = (BraidVector *)x_;
   BraidVector *y = (BraidVector *)y_;
   (y->state) = alpha * (x->state) + beta * (y->state);
   (y->guess) = alpha * (x->guess) + beta * (y->guess);
   (y->action) = alpha * (x->action) + beta * (y->action);
   (y->Delta) = alpha * (x->Delta) + beta * (y->Delta);
   return 0;
}

int MyBraidApp::SpatialNorm(braid_Vector u_,
                            double *norm_ptr)
{
   BraidVector *u = (BraidVector *)u_;
   *norm_ptr = u->state.norm();
   return 0;
}

int MyBraidApp::BufSize(int *size_ptr,
                        BraidBufferStatus &status)
{
   *size_ptr = (3 * DIM + DIM * DIM) * sizeof(double);
   return 0;
}

int MyBraidApp::BufPack(braid_Vector u_,
                        void *buffer,
                        BraidBufferStatus &status)
{
   BraidVector *u = (BraidVector *)u_;
   double *dbuffer = (double *)buffer;

   size_t bf_size = 0;

   bf_pack_help(dbuffer, u->state, DIM, bf_size);
   bf_pack_help(dbuffer, u->guess, DIM, bf_size);

   if (useDelta)
   {
      bf_pack_help(dbuffer, u->action, DIM, bf_size);
      bf_pack_help(dbuffer, u->Delta, DIM * DIM, bf_size);
   }
   status.SetSize(bf_size * sizeof(double));

   return 0;
}

int MyBraidApp::BufUnpack(void *buffer,
                          braid_Vector *u_ptr,
                          BraidBufferStatus &status)
{
   double *dbuffer = (double *)buffer;
   size_t bf_size = 0;

   BraidVector *u = new BraidVector();

   bf_unpack_help(dbuffer, u->state, DIM, bf_size);
   bf_unpack_help(dbuffer, u->guess, DIM, bf_size);

   if (useDelta)
   {
      bf_unpack_help(dbuffer, u->action, DIM, bf_size);
      bf_unpack_help(dbuffer, u->Delta, DIM * DIM, bf_size);
   }
   *u_ptr = (braid_Vector)u;

   return 0;
}

int MyBraidApp::Access(braid_Vector u_,
                       BraidAccessStatus &astatus)
{
   char filename[255];
   std::ofstream file;
   BraidVector *u = (BraidVector *)u_;

   // Extract information from astatus
   int done, level, iter, index;
   double t;
   astatus.GetTILD(&t, &iter, &level, &done);
   astatus.GetTIndex(&index);
   int nt = MyBraidApp::ntime;

   // Print information to file
   if ((index == nt - 1) && done)
   {
      sprintf(filename, "%s.%04d.%03d", "lorenz-Delta.out", index, rank);
      file.open(filename);
      pack_array(file, u->state);
      file.close();
   }

   return 0;
}

// --------------------------------------------------------------------------
// Main driver
// --------------------------------------------------------------------------

int main(int argc, char *argv[])
{
   double tstart, tstop;
   int rank;

   int nt = 4096;
   double Tf_lyap = 4;
   int max_levels = 4;
   int nrelax = 1;
   int nrelax0 = 0;
   double tol = 1e-5;
   int cfactor = 4;
   int cf0 = 4;
   int max_iter = 25;
   int newton_iters = 2;
   int DeltaLevel = 1;
   bool useFMG = false;
   bool useDelta = false;
   bool useTheta = true;

   int arg_index;

   // Initialize MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
            printf("  -tf         : set end time, in Lyapunov time (default %lf)\n", Tf_lyap);
            printf("  -ml         : set max levels\n");
            printf("  -nu         : set num F-C relaxations\n");
            printf("  -nu0        : set num F-C relaxations on level 0\n");
            printf("  -tol        : set stopping tolerance\n");
            printf("  -cf         : set coarsening factor\n");
            printf("  -cf0        : set coarsening factor for level 0\n");
            printf("  -mi         : set max iterations\n");
            printf("  -niters     : set number of newton iters for theta method\n");
            printf("  -fmg        : use FMG cycling\n");
            printf("  -Delta      : use delta correction\n");
            printf("  -Deltalvl   : defer delta correction until this level\n");
            printf("  -theta      : use first order theta method\n");
            printf("\n");
         }
         exit(0);
      }
      else if (strcmp(argv[arg_index], "-nt") == 0)
      {
         arg_index++;
         nt = atoi(argv[arg_index++]);
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
      else if (strcmp(argv[arg_index], "-Deltalvl") == 0)
      {
         arg_index++;
         DeltaLevel = atoi(argv[arg_index++]);
      }
      else if (strcmp(argv[arg_index], "-fmg") == 0)
      {
         arg_index++;
         useFMG = true;
      }
      else if (strcmp(argv[arg_index], "-Delta") == 0)
      {
         arg_index++;
         useDelta = true;
      }
      else if (strcmp(argv[arg_index], "-theta") == 0)
      {
         arg_index++;
         useTheta = true;
      }
      else
      {
         arg_index++;
         /*break;*/
      }
   }

   tstart = 0.0;
   tstop = Tf_lyap * T_lyap;

   if (useDelta && rank == 0)
   {
      std::cout << "Using Delta correction\n";
   }
   if (useTheta && rank == 0)
   {
      std::cout << "Using theta method\n";
   }

   double newton_tol = 1e-2 * tol;

   // set up app structure
   MyBraidApp app(MPI_COMM_WORLD, rank, tstart, tstop, nt, cfactor, cf0, useDelta, DeltaLevel, useTheta, newton_iters, newton_tol, max_levels);

   // Initialize Braid Core Object and set some solver options
   BraidCore core(MPI_COMM_WORLD, &app);
   // if (useDelta)
   // {
   //    core.SetResidual();
   // }
   core.SetResidual();
   if (useFMG)
   {
      core.SetFMG();
      core.SetNFMG(1);
   }
   core.SetPrintLevel(2);
   core.SetAccessLevel(1);
   core.SetMaxLevels(max_levels);
   core.SetMaxIter(max_iter);
   core.SetAbsTol(tol);
   core.SetCFactor(-1, cfactor);
   core.SetCFactor(0, cf0);
   core.SetNRelax(-1, nrelax);
   core.SetNRelax(0, nrelax0);
   core.SetSkip(0);
   core.SetStorage(1);
   core.SetTemporalNorm(3);

   // Run Simulation
   core.Drive();

   // Clean up
   MPI_Finalize();

   return (0);
}
