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
// Example:       ex-01-pp.cpp
//
// Interface:     C++
// 
// Requires:      C-language and C++ support     
//
// Compile with:  make ex-01-pp
//
// Help with:     ex-01-pp -help
//
// Sample run:    mpirun -np 2 ex-01-pp
//
// Description:   solve the scalar ODE 
//                   u' = lambda u, 
//                   with lambda=-1 and y(0) = 1
//
//                Same as ex-01, only implements more advanced XBraid features.
//                
//                When run with the default 10 time steps, the solution is:
//                $ ./ex-01-pp
//                $ cat ex-01.out.00*
//                  1.00000000000000e+00
//                  6.66666666666667e-01
//                  4.44444444444444e-01
//                  2.96296296296296e-01
//                  1.97530864197531e-01
//                  1.31687242798354e-01
//                  8.77914951989026e-02
//                  5.85276634659351e-02
//                  3.90184423106234e-02
//                  2.60122948737489e-02
//                  1.73415299158326e-02
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.hpp"

// --------------------------------------------------------------------------
// User-defined routines and objects 
// --------------------------------------------------------------------------

// Define BraidVector, can contain anything, and be named anything
// --> Put all time-dependent information here
class BraidVector
{
public:
   // Each vector holds the scalar solution value at a particular time 
   double value;

   // Construct a BraidVector for a given double 
   BraidVector(double value_) : value(value_) { }

   // Deconstructor
   virtual ~BraidVector() {};

};


// Wrapper for BRAID's App object 
// --> Put all time INDEPENDENT information here
class MyBraidApp : public BraidApp
{
protected:
   // BraidApp defines tstart, tstop, ntime and comm_t

public:
   // Constructor 
   MyBraidApp(MPI_Comm comm_t_, int rank_, double tstart_ = 0.0, double tstop_ = 1.0, int ntime_ = 100);
   
   // We will need the MPI Rank
   int rank;

   // Deconstructor
   virtual ~MyBraidApp() {};

   // Define all the Braid Wrapper routines
   // Note: braid_Vector == BraidVector*
   virtual int Step(braid_Vector    u_,
                    braid_Vector    ustop_,
                    braid_Vector    fstop_,
                    BraidStepStatus &pstatus);

   virtual int Clone(braid_Vector  u_,
                     braid_Vector *v_ptr);

   virtual int Init(double        t,
                    braid_Vector *u_ptr);

   virtual int Free(braid_Vector u_);

   virtual int Sum(double       alpha,
                   braid_Vector x_,
                   double       beta,
                   braid_Vector y_);

   virtual int SpatialNorm(braid_Vector  u_,
                           double       *norm_ptr);

   virtual int BufSize(int *size_ptr,
                       BraidBufferStatus  &status);

   virtual int BufPack(braid_Vector  u_,
                       void         *buffer,
                       BraidBufferStatus  &status);

   virtual int BufUnpack(void         *buffer,
                         braid_Vector *u_ptr,
                         BraidBufferStatus  &status);

   virtual int Access(braid_Vector       u_,
                      BraidAccessStatus &astatus);

   // Not needed in this example
   virtual int Residual(braid_Vector     u_,
                        braid_Vector     r_,
                        BraidStepStatus &pstatus) { return 0; }

   // Not needed in this example
   virtual int BufAlloc(void              **buffer,
                        int               nbytes,
                        BraidBufferStatus &bstatus) { return 0; }

   // Not needed in this example
   virtual braid_Int BufFree(void          **buffer) { return 0; }

   // Not needed in this example
   virtual int Coarsen(braid_Vector   fu_,
                       braid_Vector  *cu_ptr,
                       BraidCoarsenRefStatus &status) { return 0; }

   // Not needed in this example
   virtual int Refine(braid_Vector   cu_,
                      braid_Vector  *fu_ptr,
                      BraidCoarsenRefStatus &status)  { return 0; }

};

// Braid App Constructor
MyBraidApp::MyBraidApp(MPI_Comm comm_t_, int rank_, double tstart_, double tstop_, int ntime_)
   : BraidApp(comm_t_, tstart_, tstop_, ntime_)
{
   rank = rank_;
}

// 
int MyBraidApp::Step(braid_Vector    u_,
                     braid_Vector    ustop_,
                     braid_Vector    fstop_,
                     BraidStepStatus &pstatus)
{
   
   BraidVector *u = (BraidVector*) u_;
   double tstart;             // current time
   double tstop;              // evolve to this time

   // Get time step information
   pstatus.GetTstartTstop(&tstart, &tstop);

   // Use backward Euler to propagate solution
   (u->value) = 1./(1. + tstop-tstart)*(u->value);
   
   // no refinement
   pstatus.SetRFactor(1);
   
   return 0;

}

int MyBraidApp::Init(double        t,
                       braid_Vector *u_ptr)
{
   BraidVector *u = new BraidVector(0.0);
   if (t != tstart)
   {
      u->value = 0.456;
   }
   else
   {
      u->value = 1.0;
   }

   *u_ptr = (braid_Vector) u;
   return 0;

}

int MyBraidApp::Clone(braid_Vector  u_,
                        braid_Vector *v_ptr)
{
   BraidVector *u = (BraidVector*) u_;
   BraidVector *v = new BraidVector(u->value); 
   *v_ptr = (braid_Vector) v;

   return 0;
}


int MyBraidApp::Free(braid_Vector u_)
{
   BraidVector *u = (BraidVector*) u_;
   delete u;
   return 0;
}

int MyBraidApp::Sum(double       alpha,
                      braid_Vector x_,
                      double       beta,
                      braid_Vector y_)
{
   BraidVector *x = (BraidVector*) x_;
   BraidVector *y = (BraidVector*) y_;
   (y->value) = alpha*(x->value) + beta*(y->value);
   return 0;
}

int MyBraidApp::SpatialNorm(braid_Vector  u_,
                              double       *norm_ptr)
{
   double dot;
   BraidVector *u = (BraidVector*) u_;
   dot = (u->value)*(u->value);
   *norm_ptr = sqrt(dot);
   return 0;
}

int MyBraidApp::BufSize(int                *size_ptr,
                          BraidBufferStatus  &status)                           
{
   *size_ptr = sizeof(double);
   return 0;
}

int MyBraidApp::BufPack(braid_Vector       u_,
                          void               *buffer,
                          BraidBufferStatus  &status)
{
   BraidVector *u = (BraidVector*) u_;
   double *dbuffer = (double *) buffer;

   dbuffer[0] = (u->value);
   status.SetSize(sizeof(double));

   return 0;
}

int MyBraidApp::BufUnpack(void              *buffer,
                            braid_Vector      *u_ptr,
                            BraidBufferStatus &status)
{
   double *dbuffer = (double *) buffer;
   
   BraidVector *u = new BraidVector(dbuffer[0]); 
   *u_ptr = (braid_Vector) u;

   return 0;
}

int MyBraidApp::Access(braid_Vector       u_,
                         BraidAccessStatus &astatus)
{
   char       filename[255];
   FILE      *file;
   BraidVector *u = (BraidVector*) u_;

   // Extract information from astatus
   int done, level, iter, index;
   double t;
   astatus.GetTILD(&t, &iter, &level, &done);
   astatus.GetTIndex(&index);
   
   // Print information to file
   sprintf(filename, "%s.%04d.%03d", "ex-01.out", index, rank);
   file = fopen(filename, "w");
   fprintf(file, "%.14e\n", (u->value));
   fflush(file);
   fclose(file);

   return 0;
}


// --------------------------------------------------------------------------
// Main driver
// --------------------------------------------------------------------------

int main (int argc, char *argv[])
{

   double        tstart, tstop;
   int           ntime, rank;

   // Define time domain: ntime intervals
   ntime  = 10;
   tstart = 0.0;
   tstop  = tstart + ntime/2.;
  
   // Initialize MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
   // set up app structure
   MyBraidApp app(MPI_COMM_WORLD, rank, tstart, tstop, ntime);

   // Initialize Braid Core Object and set some solver options
   BraidCore core(MPI_COMM_WORLD, &app);
   core.SetPrintLevel(2);
   core.SetMaxLevels(2);
   core.SetAbsTol(1.0e-6);
   core.SetCFactor(-1, 2);
   
   // Run Simulation
   core.Drive();

   // Clean up
   MPI_Finalize();

   return (0);
}



