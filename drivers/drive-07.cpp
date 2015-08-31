// Copyright (c) 2013, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. Written by
// Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin
// Dobrev, et al. LLNL-CODE-660355. All rights reserved.
//
// This file is part of XBraid. Email xbraid-support@llnl.gov for support.
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License (as published by the Free Software
// Foundation) version 2.1 dated February 1999.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 59
// Temple Place, Suite 330, Boston, MA 02111-1307 USA
//

/* Drive 07 2D/3D Diffusion equation . Requires MFEM, Hpyre, Metis and GlVis

         Interface: MFEM 

         Compile with: make drive-07 -- Modify Makefile to point to metis, mfem and hypre libraries

         Help with:    drive-06 -help

         Sample run:   mpirun -n12 ./drive-07 -prob -1

         Description:  This code solves the 2D/3D diffusion equaiton with
                        time dependent coefficients 

                                       u_t - div( A(t) grad(u) ) = b,
                                             u(x,0) = U0(x)                                  
                                    A(t)grad(u) . n = g(x,t) on boundary

                        on a given mesh using MFEM where A(t) is a d dimentional tensor

                                       A(t) = [a(t) , 0; 0 , b(t)] in 2D
                                 A(t) = [ a(t), 0 , 0; 0, c(t), 0; 0, 0, b(t)] in 3D

                     a(t) and b(t) are specidied in the problem class. Examples include
                                 prob < 0: a(t) = b(t) = c(t) = -prob
                                 prob = 10: a(t) = 1, c(t) = 1, b(t) = cosine
                                 prob = 11: a(t) = 1, c(t) = 1, b(t) ~ t
                                 prob = 12: a(t) = 1, c(t) = 1, b(t) jumps between 0.001 and 1*/

#include "braid_mfem.hpp"
#include <fstream>
#include <iostream>

//Extra Hypre Functions
#include "hypre_extra.hpp"
using namespace hypre;

using namespace std;

/*************************************************************************************
  Class that contains all coefficients and variables relating to the manufactured solution */
class Problem
{
   protected:
      static double kap;  //Spatial frequency of manufactured solution
      static double tau;  //Temporal frequency of manufactured solution
      static double problem; // Choose the coefficient
      static double tstop; // final stopping time
   public:

      ConstantCoefficient one;          //One 

      MatrixFunctionCoefficient *dc;     //Diffusion Coefficient
      FunctionCoefficient *source;      //Source Term
      VectorFunctionCoefficient *nbc; //Boundary Conditions
      FunctionCoefficient *exact_sol;    //Exact Solution
      FunctionCoefficient *ic;			 // Initial Condition 

      Problem(int dim, double _kap, double _tau, double _problem, double _tstop);

      //Define the time dependent coefficients in the Tensor
      static double a(double t);
      static double b(double t);
      static double c(double t);

      //Neuman boundary Condition
      static void NBC(const Vector &p, double t, Vector &v);

      //Exact Solution
      static double ExSol(Vector &p, double t);

      //Source term
      static double ST(Vector &p, double t);

      //Initial Condition
      static double IC(Vector &x);

      //Diffusion Coefficient
      static void DC(const Vector &p, double t, DenseMatrix &m);

      virtual ~Problem();
};

/*************************************************************************************
  Struct to hold all command line and defualt options */
struct DiffusionOptions : public BraidOptions
{  

   int ode_solver_type; 
   double dt; 
   int sc_print_level; //Spatial coarsening print level

   double kap, tau;    
   double problem;

   //Visualization Options
   const char *vishost;  
   int         visport;
   int vis_time_steps;   // Visualize every x time steps
   int vis_braid_steps;  // Visualize every x braid iterations 
   bool vis_screenshots; // Save Screenshots t/f

   Problem *prob;
   HypreParVector *U;

   DiffusionOptions(int argc, char *argv[]);

   void InitializeProblem(int dim);
   
   virtual ~DiffusionOptions();
};

/*************************************************************************************
  A time dependent Operator that solves the nonlinear system ode using newtons method */
class DiffusionODE : public TimeDependentOperator
{
   private:
      DiffusionOptions *opts;

      ParBilinearForm *a;	            
      ParLinearForm *b_form;
      Coefficient *b_coeff;
      VectorCoefficient *b_vcoeff;
      MatrixFunctionCoefficient *d_coeff;

      mutable HypreParVector *b;
      mutable HypreParMatrix *A;	
      HypreParMatrix *B, *M;
      HypreParVector *X, *Y, *Z;
      HypreBoomerAMG *amg;
      HyprePCG *pcg;
      HypreBoomerAMG *B_amg;
      HyprePCG *B_pcg;

      double current_dt;
   public:
      DiffusionODE(ParBilinearForm *_a,
            HypreParMatrix  *_M,
            ParLinearForm   *_b_form,
            DiffusionOptions *_opts);

      //Assemble diffusion matrix
      void AssembleDiffusionMatrix() const;	 

      // Assemble the rhs
      void AssembleBVector() const;

      //Diffusion Multiplication -- Not implimented yet -- Not needed??
      virtual void Mult(const Vector &x, Vector &y) const;

      //Diffusion Solve using Newtons Method 
      virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k);

      virtual ~DiffusionODE(); 
};

/*************************************************************************************
  Class defining the DiffusionApp */ 
class DiffusionApp : public MFEMBraidApp
{
   protected:

      DiffusionOptions *opts;
      FiniteElementCollection *fe_coll;
      Problem *prob;

      //Allocate Arrays needed by the DiffusionApp
      virtual void AllocLevels(int num_levels);

      //Construct the Meshes on each level
      virtual ParFiniteElementSpace *ConstructFESpace(ParMesh *pmesh);

      //Init the ode/solver/max_dt for each level
      virtual void InitLevel(int l);

   public:
      SpaceTimeMeshInfo MeshInfo;	

      DiffusionApp(DiffusionOptions *opts, MPI_Comm comm_t_, ParMesh *pmesh);

      //Modified MFEMBraidApp::Step to also record space-time mesh info
       virtual int Step(braid_Vector    u_,
                   braid_Vector    ustop_,
                   braid_Vector    fstop_,
                   BraidStepStatus &pstatus);  

      virtual ~DiffusionApp();

};

/*************************************************************************************
  Main */
int main (int argc, char *argv[])
{
   //Init MPI
   MPI_Init(&argc, &argv);

   MPI_Comm comm = MPI_COMM_WORLD;
   int myid;
   MPI_Comm_rank(comm, &myid);

   //Parse Command Line
   DiffusionOptions opts(argc, argv);

   if (!opts.Good())
   {
      if (myid == 0)
      {
         opts.PrintUsage(cout);
      }
      MPI_Finalize();
      return 1;
   }
   // Print the used options
   if (myid == 0)
   {
      opts.PrintOptions(cout);
   }

   //Load the Mesh and refine in serial
   Mesh *mesh = opts.LoadMeshAndSerialRefine();
   
   if (!mesh)
   {
      if (myid == 0)
      {
         cerr << "\nError loading mesh file: " << opts.mesh_file
            << '\n' << endl;
      }
      MPI_Finalize();
      return 2;
   }

   //Init the problem class -- needs dimension of the mesh
   int dim = mesh->Dimension();
   opts.InitializeProblem(dim);
   
   // Split comm (MPI_COMM_WORLD) into spatial and temporal communicators
   MPI_Comm comm_x, comm_t;
   BraidUtil util;
   util.SplitCommworld(&comm, opts.num_procs_x, &comm_x, &comm_t);

   // Partition the mesh accross the spatial communicator
   ParMesh *pmesh = new ParMesh(comm_x, *mesh);
   delete mesh;

   //Init braid App
   DiffusionApp app(&opts, comm_t, pmesh);

   //Run braid
   BraidCore core(comm, &app);
   opts.SetBraidCoreOptions(core);
   core.Drive();   

   if (opts.sc_print_level)
      app.MeshInfo.Print(comm);

   MPI_Finalize();
   return 0;
}	

/*************************************************************************************
  Implimentaion::Problem */

//Initailize static data members
double Problem::kap = 0;
double Problem::tau = 0;
double Problem::problem = -1;
double Problem::tstop = 1;

Problem::Problem(int dim, double _kap, double _tau, double _problem, double _tstop) :one(1.0)
{
   kap = _kap; tau = _tau; problem = _problem; tstop = _tstop;
   source = new FunctionCoefficient(ST);
   exact_sol = new FunctionCoefficient(ExSol);
   ic = new FunctionCoefficient(IC);
   nbc = new VectorFunctionCoefficient(dim, NBC);
   dc = new MatrixFunctionCoefficient(dim, DC);
}

//problem < 0 are all time independent problems
//problem > 0 are all time dependent problems
double Problem::a(double t)
{  
   if (problem < 0)
      return -problem;  // run with both = -problem
   else 
      return 1;
}  
double Problem::b(double t)
{
   if (problem < 0)
      return -problem;
   else if (problem == 10) 
      return cos(M_PI/2*t/(tstop));
   else if (problem == 11)
      return t/tstop + 0.001;
   else if (problem == 12)
   {
      if (t <= tstop/4)
         return 0.001;
      else if (t > tstop/4 && t <= tstop/2)
         return 1;
      else if (t > tstop/2 && t <= 3*tstop/4)
         return 0.001;
      else
         return 1;
   }
   else
      return 1;
}
double Problem::c(double t)
{ 
   if (problem < 0)
      return -problem;
   else
      return 1;
}   

void Problem::DC(const Vector &p, double t, DenseMatrix &m)
{   
   int dim = p.Size();
   if (dim == 2)
   {
      double aa[4] = {-Problem::a(t), 0,
         0, -Problem::b(t)};
      m = aa;
   }
   else
   {
      double aa[9] = {-Problem::a(t), 0, 0 
         , 0 , -Problem::c(t), 0,
            0, 0, -Problem::b(t)};
      m = aa;
   }   
}
void Problem::NBC(const Vector &p, double t, Vector &v)
{
   int dim = p.Size();
   if (dim == 2)
   {
      v(0) = a(t)*kap*cos(kap*p(0))*sin(kap*p(1))*sin(tau*t);
      v(1) = b(t)*kap*sin(kap*p(0))*cos(kap*p(1))*sin(tau*t);
   }
   else
   {
      v(0) = a(t)*kap*cos(kap*p(0))*sin(kap*p(1))*sin(kap*p(2))*sin(tau*t);
      v(1) = c(t)*kap*sin(kap*p(0))*cos(kap*p(1))*sin(kap*p(2))*sin(tau*t);
      v(2) = b(t)*kap*sin(kap*p(0))*sin(kap*p(1))*cos(kap*p(2))*sin(tau*t);
   }
}

double Problem::ExSol(Vector &p, double t)
{
   int dim = p.Size();
   if (dim == 2)
      return sin(kap*p(0))*sin(kap*p(1))*sin(tau*t);
   else
      return sin(kap*p(0))*sin(kap*p(1))*sin(kap*p(2))*sin(tau*t);
}

double Problem::ST(Vector &p, double t)
{  
   int dim = p.Size();
   if (dim == 2)
      return (tau*cos(tau*t) + kap*kap*(a(t)+b(t))*sin(tau*t))*sin(kap*p(0))*sin(kap*p(1));
   else
      return (tau*cos(tau*t) + kap*kap*(a(t)+b(t)+c(t))*sin(tau*t))
         *sin(kap*p(0))*sin(kap*p(1))*sin(kap*p(2));
}

double Problem::IC(Vector &x)
{
   return ExSol(x,0.0);
}

Problem::~Problem() 
{
   delete source;
   delete exact_sol;
   delete ic;
   delete nbc;
   delete dc;
}
/*************************************************************************************
  Implimentatio::DiffusionOptions */
   DiffusionOptions::DiffusionOptions(int argc, char *argv[])
: BraidOptions(argc, argv)
{
   //Set inherited default options
   num_time_steps = 32;

   //Set defualts for the mesh/refinement inherited options
   mesh_file			= "../../mfem/data/star.mesh";
   ser_ref_levels = 1;
   par_ref_levels = 1;
   AddMeshOptions();

   //Set defualts for the DiffusionOptions specific options
   kap = M_PI;
   tau = (2. + 1./6.)*M_PI;
   problem = -1;
 
   
   sc_print_level = 0;
   ode_solver_type = -1;

   vishost         = "localhost";
   visport         = 19916;
   vis_time_steps  = 0;
   vis_braid_steps = 1;
   vis_screenshots = false;

   //Parse Command Line
   AddOption(&sc_print_level, "-scprint", "--spatial-coarsening-print_level",
         "Print level for spatial coarseing"); 
   AddOption(&ode_solver_type, "-s", "--ode-solver",
         "ODE solver: 1 - Forward Euler, 2 - RK2 SSP, 3 - RK3 SSP,"
         " 4 - RK4, 6 - RK6, -1 - Backward Euler, -2 - SDIRK22 "
         " -3 - SDIRK33 -4 -SDIRK23, -5 - SDIRK34.");
   AddOption(&problem, "-prob", "-diffusion-coefficient",
         "Pick which diffusion coefficient to use -- see Problem::a(t) and b(t)");
   AddOption(&vishost, "-vh", "--visualization-host",
         "Set the GLVis host.");
   AddOption(&visport, "-vp", "--visualization-port",
         "Set the GLVis port.");
   AddOption(&vis_time_steps, "-vts", "--visualize-time-steps",
         "Visualize every n-th time step (0:final only).");
   AddOption(&vis_braid_steps, "-vbs", "--visualize-braid-steps",
         "Visualize every n-th Braid step (0:final only).");
   AddOption(&vis_screenshots, "-vss", "--vis-screenshots",
         "-no-vss", "--vis-no-screenshots", "Enable/disable the saving of"
         " screenshots of the visualized data.");

   Parse();

   dt = (t_final - t_start) / num_time_steps;
}

void DiffusionOptions::InitializeProblem(int dim)
{
   prob = new Problem(dim, kap, tau, problem, t_final);
}  

DiffusionOptions::~DiffusionOptions()
{
   delete prob;
}

 

/*************************************************************************************
  Implimentation::DiffusionODE */
DiffusionODE::DiffusionODE(ParBilinearForm *_a,
      HypreParMatrix  *_M,
      ParLinearForm   *_b_form,
      DiffusionOptions *_opts) :
   a(_a), M(_M), b_form(_b_form), b_coeff(_opts->prob->source), b_vcoeff(_opts->prob->nbc),
   d_coeff(_opts->prob->dc), opts(_opts)
{        
   d_coeff->SetTime(GetTime());
   a->Assemble();
   a->Finalize();
   A = a->ParallelAssemble();
   width = A->Width();
   height = A->Height();

   double tmp; // workaround to avoid memory leak
   X = new HypreParVector(A->GetComm(), A->GetGlobalNumCols(), &tmp,
         A->GetColStarts());
   Y = new HypreParVector(A->GetComm(), A->GetGlobalNumCols(), &tmp,
         A->GetColStarts());
   //Set Up BoomerAMG for taking inverse of mass matrix 
   amg = new HypreBoomerAMG(*M);
   amg->SetPrintLevel(-1);
   pcg = new HyprePCG(*M);
   pcg->SetTol(1e-15);
   pcg->SetMaxIter(2000);
   pcg->SetPrintLevel(0);
   pcg->SetPreconditioner(*amg);
   pcg->SetZeroInintialIterate();

   b = NULL;
   B = NULL; 
   B_amg = NULL;
   B_pcg = NULL;

   Z = new HypreParVector(*A);
   current_dt = -1;
}		 

void DiffusionODE::AssembleDiffusionMatrix() const
{
   //If diffusion coefficients are time dependent then update matrix
   if (opts->problem > 0)
   {
      delete A;
      d_coeff->SetTime(GetTime());
      a->Update();
      a->Assemble();
      a->Finalize();
      A = a->ParallelAssemble();
   }   
}

void DiffusionODE::AssembleBVector() const
{

   b_coeff->SetTime(GetTime());
   b_vcoeff->SetTime(GetTime());
   b_form->Assemble();

   if (!b)
      b = b_form->ParallelAssemble();
   else
      b_form->ParallelAssemble(*b);
}

void DiffusionODE::Mult(const Vector &x, Vector &y) const
{
   X->SetData(x.GetData());
   Y->SetData(y.GetData());

   AssembleBVector();

   AssembleDiffusionMatrix();

   A->Mult(*X, *Z); // M^{-1} A X
   if (b)
      *Z += *b; // M^{-1} (A X + b)
   pcg->Mult(*Z, *Y);
}

//Iterate over the system M*k^(j+1) = A(x + dt*k(j)) + b(x,t)
void DiffusionODE::ImplicitSolve(const double dt, const Vector &x, Vector &k)
{
   //Wrap the data (X == x, Y == k)
   X->SetData(x.GetData());
   Y->SetData(k.GetData());

   //Assemble the source term vector 
   AssembleBVector();
   AssembleDiffusionMatrix();  

   //Form the matrix: B = M - dt A
   if (!B)
      B = new HypreParMatrix(hypre_ParCSRMatrixAdd(*M, *A));
   //If A is time dependent or if the time step size has changed from B
   if (opts->problem > 0 || current_dt != dt)
   {
      hypre_ParCSRMatrixSetConstantValues(*B, 0.0);
      hypre_ParCSRMatrixSum(*B, 1.0, *M);
      hypre_ParCSRMatrixSum(*B, -dt, *A);

      //Set up BOOMERamg for inverting B

      delete B_amg;
      B_amg = new HypreBoomerAMG(*B);
      B_amg->SetPrintLevel(-1);

      delete B_pcg;
      B_pcg = new HyprePCG(*B);
      B_pcg->SetTol(1e-15);
      B_pcg->SetMaxIter(2000);
      B_pcg->SetPrintLevel(0);
      B_pcg->SetPreconditioner(*B_amg);
   }
   //Form the rhs: A x + b
   A->Mult(*X, *Z);
   if (b)
      *Z += *b;

   //Solve the system B Y = Z
   *Y = *opts->U;  //Set the initial guess
   B_pcg->Mult(*Z, *Y);
}

DiffusionODE::~DiffusionODE() 
{
   delete X;
   delete Y;
   delete Z;
   delete A;   
   delete B;
   delete b;

}


/*************************************************************************************
  Implimentation::DiffusionApp */
DiffusionApp::DiffusionApp(
      DiffusionOptions *_opts, MPI_Comm comm_t_, ParMesh *pmesh)
: MFEMBraidApp(comm_t_, _opts->t_start, _opts->t_final, _opts->num_time_steps),
   opts(_opts), MeshInfo(opts->max_levels)

{
   //Init all meshes/odes/solvers/etc for all space levels
   fe_coll = new LinearFECollection;
   InitMultilevelApp(pmesh, opts->par_ref_levels, opts->spatial_coarsen);

   //Set the initial/exact solutions
   x[0]->ProjectCoefficient(*(opts->prob->ic));
   HypreParVector *U = x[0]->GetTrueDofs();
   SetInitialCondition(U);
   SetExactSolution(opts->prob->exact_sol);

   //Init glvis visualization
   SetVisHostAndPort(opts->vishost, opts->visport);
   SetVisSampling(opts->vis_time_steps, opts->vis_braid_steps);
   SetVisScreenshots(opts->vis_screenshots);  
}

int DiffusionApp::Step(braid_Vector    u_,
                   braid_Vector    ustop_,
                   braid_Vector    fstop_,
                   BraidStepStatus &pstatus)
{
   // This contains one small change over the default Step, we store the Space-Time mesh info
   
   // Extract info from braid_vector
   BraidVector *u = (BraidVector*) u_;
   BraidVector *ustop = (BraidVector*) ustop_;
   
   int braid_level;
   int braid_iter;

   pstatus.GetLevel(&braid_level);
   pstatus.GetIter(&braid_iter);
   
   double tstart, tstop, dt;
   pstatus.GetTstartTstop(&tstart, &tstop);
   dt = tstop - tstart;

   //Store run-time mesh info
   if (opts->sc_print_level)
        MeshInfo.SetRow(braid_level, u->level, x[u->level]->ParFESpace()->GetParMesh(), dt);
   opts->U = ustop;
   MFEMBraidApp::Step(u_,ustop_, fstop_, pstatus);       

   return 0;
}

void DiffusionApp::AllocLevels(int num_levels)
{
}

ParFiniteElementSpace* DiffusionApp::ConstructFESpace(ParMesh *pmesh)
{
   return new ParFiniteElementSpace(pmesh, fe_coll);
}

void DiffusionApp::InitLevel(int l)
{ 
   //Create Bilinear and linear forms
   ParBilinearForm *a = new ParBilinearForm(fe_space[l]);
   ParBilinearForm *m = new ParBilinearForm(fe_space[l]);
   ParLinearForm   *b = new ParLinearForm(fe_space[l]);

   //Add Boundary and domain integrators 
   m->AddDomainIntegrator(new MassIntegrator(opts->prob->one));
   b->AddDomainIntegrator(new DomainLFIntegrator(*opts->prob->source));
   b->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(*opts->prob->nbc));
   a->AddDomainIntegrator(new DiffusionIntegrator(*opts->prob->dc));

   // Local assembly. Use precomputed sparsity to ensure that 'a' and 'm' have
   // the same sparsity patterns.
   a->UsePrecomputedSparsity();
   m->UsePrecomputedSparsity();

   //Finalize the constant Mass Matrix 
   m->Assemble();
   m->Finalize();
   HypreParMatrix *M = m->ParallelAssemble();
   delete m;

   //Create a nonlinearsystem to be solved on this level
   ode[l] = new DiffusionODE(a, M, b, opts );

   //Set ODESolver and inititialize with ode[l]
   ODESolver *ode_solver = NULL;
   switch (opts->ode_solver_type)
   {
      case 1: ode_solver = new ForwardEulerSolver; break;
      case 2: ode_solver = new RK2Solver(1.0); break;
      case 3: ode_solver = new RK3SSPSolver; break;
      case 4: ode_solver = new RK4Solver; break;
      case 6: ode_solver = new RK6Solver; break;
      case -1: ode_solver = new BackwardEulerSolver; break;
      case -2: ode_solver = new SDIRK23Solver(2); break;
      case -3: ode_solver = new SDIRK33Solver; break;
      case -4: ode_solver = new SDIRK23Solver; break;
      case -5: ode_solver = new SDIRK34Solver; break;

      default: ode_solver = new BackwardEulerSolver; break;
   }
   solver[l] = ode_solver;
   solver[l]->Init(*ode[l]);
	
   // Set max_dt[l] -- the maximum time step allowed on this level
   if ( l == 0 && opts->spatial_coarsen)
   {        
      for (int i = opts->par_ref_levels; i >= 0; i--)
      {
         max_dt[i] = (opts->t_final-opts->t_start)/(opts->min_coarse + 1 ) / pow(opts->cfactor, opts->par_ref_levels - i);
      }
   }  


}

DiffusionApp:: ~DiffusionApp()
{
   delete fe_coll;
}			  
