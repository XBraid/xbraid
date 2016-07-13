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



// Example 04
//
// Compile with: make drive-04
//
// Sample runs:  mpirun -np 4 drive-04
//               mpirun -np 4 drive-04 -pref 1 -px 2 -nt 100 -odesolver -2
//               mpirun -np 4 drive-04 -mesh ../../mfem/data/beam-quad.mesh -sref 2 -nobraid
//               mpirun -np 16 drive-04 -mesh ../../mfem/data/escher-p2.mesh -nobraid
//
// Description:  Unstructured 2D/3D ODE problems in MFEM.

#include <fstream>
#include "braid_mfem.hpp"
#include "hypre_extra.hpp"

using namespace std;
using namespace mfem;
using namespace hypre;


// Right-hand side for the set of scalar ODEs du_i/dt = f(u_i,t)
class ScalarODE : public TimeDependentOperator
{
private:
   int option;

   inline double f(double u, double t) const
   {
      if (option == 0)
         return 1.0; // u(t) = t+C
      else if (option == 1)
         return 4*pow(t,3)-3*pow(t,2)+2*t-1; // u(t) = t^4-t^3+t^2-t+C
      else if (option == 2)
         return u*u+1; // u(t) = tan(t+C)
      else if (option == 3)
         return 2*t*u+t; // u(t) = -1/2 + exp(t^2)*C
      else if (option == 4)
         return -2000.*(u - cos(t)); // stiff
      else
         return u*u+1; // option 2
   }

   inline double ik(const double &u, const double &t, const double &dt)
   {
      if (option == 0 || option == 1)
         return f(0, t); // f(.,.) is independent of u
      else if (option == 2)
         return ((1-2*dt*u) - sqrt(1-4*dt*(u+dt)))/(2*dt*dt); // f = u^2+1
      else if (option == 3)
         return (2*t*u+t)/(1-2*t*dt); // f = 2*t*u+t
      else if (option == 4)
         return -2000*(u - cos(t))/(1. + 2000.*dt); // stiff
      else
         return ((1-2*dt*u) - sqrt(1-4*dt*(u+dt)))/(2*dt*dt); // option 2
   }

public:
   ScalarODE(int _size, int _option=2) : TimeDependentOperator(_size)
   { option=_option; }

   virtual void Mult(const Vector &x, Vector &y) const
   {
      double t = GetTime();
      for (int i = 0; i < x.Size(); i++)
         y(i) = f(x(i),t);
   }

   /// Solve the equation: k = f(x + dt*k, t), for the unknown k.
   virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k)
   {
      double t = GetTime();
      for (int i = 0; i < x.Size(); i++)
         k(i) = ik(x(i), t, dt);

      if (k.CheckFinite() > 0)
         cerr << "Detected nan/inf!" << endl;
   }
};


// Right-hand side for the linear system of ODEs
//               M du/dt = A u + b(t),
// where A and M are time-dependent.
class LinearSystemODE : public TimeDependentOperator
{
private:
   HypreParMatrix *M;
   mutable HypreParMatrix *A;
   HypreParVector *X, *Y, *Z;
   mutable HypreParVector *b;
   ParLinearForm *b_form;
   mutable ParBilinearForm *aform;
   Coefficient *b_coeff;
   VectorCoefficient *b_vcoeff;
   HypreBoomerAMG *amg;
   HyprePCG *pcg;
   mutable Array<double> dts;
   mutable Array<HypreParMatrix*> B; // B = M - dt A, for ImplicitSolve
   mutable Array<HypreBoomerAMG*> B_amg;
   mutable Array<CGSolver*> B_pcg;
   int timedep_A; // Is 'a' time-dependent?

   int GetDtIndex(double dt) const
   {
      for (int i = 0; i < dts.Size(); i++)
      {
         if (std::abs(dts[i]-dt) < 1e-10*dt)
         {
            if (timedep_A)
            {
               MFEM_ABORT("TODO");
            }
            return i;
         }
      }
      dts.Append(dt);
      B.Append(new HypreParMatrix(hypre_ParCSRMatrixAdd(*M, *A)));
      hypre_ParCSRMatrixSetConstantValues(*B.Last(), 0.0);
      hypre_ParCSRMatrixSum(*B.Last(), 1.0, *M);
      hypre_ParCSRMatrixSum(*B.Last(), -dt, *A);

      B_amg.Append(new HypreBoomerAMG);
      B_amg.Last()->SetPrintLevel(-1);

      B_pcg.Append(new CGSolver(B.Last()->GetComm()));
      B_pcg.Last()->SetRelTol(1e-12);
      B_pcg.Last()->SetMaxIter(200);
      B_pcg.Last()->SetPrintLevel(0);
      B_pcg.Last()->SetPreconditioner(*B_amg.Last());
      B_pcg.Last()->SetOperator(*B.Last());
      B_pcg.Last()->iterative_mode = false;
      return dts.Size()-1;
   }

public:
   LinearSystemODE(ParBilinearForm *_aform,
                   int _timedep_A = 1,
                   HypreParMatrix *_M = NULL,
                   ParLinearForm *_b_form = NULL,
                   Coefficient *_b_coeff = NULL,
                   VectorCoefficient *_b_vcoeff = NULL)
      : M(_M), b_form(_b_form), aform(_aform), b_coeff(_b_coeff),
        b_vcoeff(_b_vcoeff), timedep_A(_timedep_A)
   {
      // Inital Assembly of Stiffness Matrix
      aform->Assemble();
      aform->Finalize();
      A = aform->ParallelAssemble();
      height = A->Height();
      width = A->Width();

      double tmp; // workaround to avoid memory leak
      X = new HypreParVector(A->GetComm(), A->GetGlobalNumCols(), &tmp,
                             A->GetColStarts());
      Y = new HypreParVector(A->GetComm(), A->GetGlobalNumCols(), &tmp,
                             A->GetColStarts());

      if (!M)
      {
         amg = NULL;
         pcg = NULL;
         Z = NULL;
      }
      else
      {
         amg = new HypreBoomerAMG(*M);
         amg->SetPrintLevel(-1);

         pcg = new HyprePCG(*M);
         pcg->SetTol(1e-12);
         pcg->SetMaxIter(200);
         pcg->SetPrintLevel(0);
         pcg->SetPreconditioner(*amg);
         pcg->SetZeroInintialIterate();

         Z = new HypreParVector(*A);
      }

      b = NULL;
   }

   void AssembleDiffMatrix() const
   {
      if (timedep_A)
      {
         delete A;

         aform->Update();
         aform->Assemble();
         aform->Finalize();
         A = aform->ParallelAssemble();
      }
   }

   // Assemble the Source Term Vector at each time t
   void AssembleBVector() const
   {
      if (b_form)
      {
         // set the time for the b-coefficients
         b_coeff->SetTime(GetTime());
         b_vcoeff->SetTime(GetTime());

         // assemble b
         b_form->Assemble();
         if (!b)
            b = b_form->ParallelAssemble();
         else
            b_form->ParallelAssemble(*b);
      }
   }

   virtual void Mult(const Vector &x, Vector &y) const
   {
      // Wrap the Data
      X->SetData(x.GetData());
      Y->SetData(y.GetData());

      // Assemble
      AssembleBVector();
      AssembleDiffMatrix();

      // Solve
      if (!M)
      {
         A->Mult(*X, *Y); // A X
         if (b)
            *Y += *b; // A X + b
      }
      else
      {
         A->Mult(*X, *Z); // Z = A X
         if (b)
            *Z += *b; // Z = (A X + b)

         pcg->Mult(*Z, *Y); // Y = M^{-1} Z
      }
   }

   /** Solve for k in the equation: M k = A(x + dt k, t) + b(t).
       This method is needed for backward Euler and other DIRK methods. */
   virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k)
   {
      // Wrap the Data
      X->SetData(x.GetData());
      Y->SetData(k.GetData());

      // Assemble
      AssembleBVector();
      AssembleDiffMatrix();

      // 1) Form the matrix: B = M - dt A
      int i = GetDtIndex(dt);

      // 2) Form the rhs: A x + b
      A->Mult(*X, *Z);
      if (b)
         *Z += *b;

      // 3) Solve the system B Y = Z
      B_pcg[i]->Mult(*Z, *Y);
   }

   double InnerProduct(const Vector &x, const Vector &y) const
   {
         // Wrap the Data
         X->SetData(x.GetData());
         Y->SetData(y.GetData());

         M->Mult(*X, *Z, 1.0, 0.0);
         return mfem::InnerProduct(*Z, *Y);
   }

   virtual ~LinearSystemODE()
   {
      for (int i = dts.Size()-1; i >= 0; i--)
      {
         delete B_pcg[i];
         delete B_amg[i];
         delete B[i];
      }
      delete b;
      delete A;
      delete Z;
      delete pcg;
      delete amg;
      delete Y;
      delete X;
   }
};


class DiffBraidApp : public MFEMBraidApp
{
protected:
   LinearSystemODE *ode;

public:
   DiffBraidApp(MPI_Comm comm_t_, LinearSystemODE *ode_,
                HypreParVector *X0_, ParGridFunction *x_, ODESolver *solver_,
                double tstart_, double tstop_, int ntime_)
      : MFEMBraidApp(comm_t_, ode_, X0_, x_, solver_,
                     tstart_, tstop_, ntime_),
        ode(ode_)
   { }

   virtual int SpatialNorm(braid_Vector  u_,
                           double       *norm_ptr)
   {
      BraidVector *u = (BraidVector*) u_;
      *norm_ptr = sqrt(ode->InnerProduct(*u, *u));
      return 0;
   }
};


// Solve the given ODE with the given initial condtion using uniform
// time-stepping and (optionally) visualize the solution.
void SolveODE(TimeDependentOperator *ode, HypreParVector *X0,
              ODESolver *solver, double tstart = 0.0, double tstop = 1.0,
              int ntime = 100, ParGridFunction *x = NULL,
              Coefficient *exsol = NULL)
{
   solver->Init(*ode);

   double t = tstart;
   double dt = (tstop-tstart)/ntime;

   socketstream *sol_sock = NULL;
   ParMesh *pmesh = NULL;
   if (x)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      pmesh = x->ParFESpace()->GetParMesh();
      sol_sock = new socketstream(vishost, visport);

      *x = *X0;
      *sol_sock << "parallel " << pmesh->GetNRanks()
                << " " << pmesh->GetMyRank() << "\n";
      *sol_sock << "solution\n";
      sol_sock->precision(8);
      pmesh->Print(*sol_sock);
      x->Save(*sol_sock);
      *sol_sock << "valuerange -1 1\n";
      *sol_sock << "autoscale off\n";
      //*sol_sock << "autoscale value\n";
      *sol_sock << "keys cmAaa\n";
      *sol_sock << flush;
   }

   int good, all_good = 1;
   for (int i = 0; i < ntime; i++)
   {
      solver->Step(*X0, t, dt);
      if (x && all_good)
      {
         good = sol_sock->good();
         MPI_Allreduce(&good, &all_good, 1, MPI_INT, MPI_LAND,
                       pmesh->GetComm());
         if (!all_good)
         {
            if (pmesh->GetMyRank() == 0)
               cout << "\nVisualization closed." << endl;
            continue;
         }

         *x = *X0;
         *sol_sock << "parallel " << pmesh->GetNRanks()
                   << " " << pmesh->GetMyRank() << "\n";
         *sol_sock << "solution\n";
         pmesh->Print(*sol_sock);
         x->Save(*sol_sock);
         *sol_sock << flush;
      }
   }

   if (x && exsol)
   {
      *x = *X0;

      exsol->SetTime(t);

      double err = x->ComputeL2Error(*exsol);

      if (pmesh->GetMyRank() == 0)
         cout << "\nL2 norm of the error = " << err << endl;
   }

   delete sol_sock;
}


// Exact solution parameters
const double tau = (2.+1./6.)*M_PI;
double kap;

// Exact solution
double ExSol(Vector &p, double t)
{
   int dim = p.Size();

   if (dim == 2)
      return sin(kap*p(0))*sin(kap*p(1))*sin(tau*t);
   else
      return sin(kap*p(0))*sin(kap*p(1))*sin(kap*p(2))*sin(tau*t);
}

// Initial condition
double IC(Vector &x)
{
   return ExSol(x, 0.0);
}

// Neumann boundary conditions.
// Defined as a vector field, g, so that exact solutions will work with
// arbitrary domains: du/dn|_{\partial\Omega} = g.n.
void NBC(const Vector &p, double t, Vector &v)
{
   int dim = p.Size();
   if (dim == 2)
   {
      v(0) = kap*cos(kap*p(0))*sin(kap*p(1))*sin(tau*t);
      v(1) = kap*sin(kap*p(0))*cos(kap*p(1))*sin(tau*t);
   }
   else
   {
      v(0) = kap*cos(kap*p(0))*sin(kap*p(1))*sin(kap*p(2))*sin(tau*t);
      v(1) = kap*sin(kap*p(0))*cos(kap*p(1))*sin(kap*p(2))*sin(tau*t);
      v(2) = kap*sin(kap*p(0))*sin(kap*p(1))*cos(kap*p(2))*sin(tau*t);
   }
}

// Source term
double ST(Vector &p, double t)
{
   int dim = p.Size();
   if (dim == 2)
   {
      return ((tau*cos(tau*t) + 2*kap*kap*sin(tau*t))*
              sin(kap*p(0))*sin(kap*p(1)));
   }
   else
   {
      return ((tau*cos(tau*t) + 3*kap*kap*sin(tau*t))*
              sin(kap*p(0))*sin(kap*p(1))*sin(kap*p(2)));
   }
}

int main(int argc, char *argv[])
{
   // Initialize MPI and split the space and time communicators
   MPI_Comm comm, comm_x, comm_t;
   int num_procs, myid, myid_x, myid_t;
   comm = MPI_COMM_WORLD;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(comm, &num_procs);
   MPI_Comm_rank(comm, &myid);

   // Variables used by BraidUtil
   BraidUtil util;
   int correct = 1;
   double test_t, test_dt;

   // Default parameters:
   const char *meshfile = "../../mfem/data/star.mesh";
   int    sref          = 1;
   int    pref          = 1;
   int    use_braid     = 1;
   int    num_procs_x   = 1;
   double tstart        = 0.0;
   double tstop         = 1.0;
   int    ntime         = 32;

   kap = M_PI;

   // double cfl         = 1.0;
   int    fe_degree   = 1;
   int    ode_solver  = -1;

   int heat_equation = 1;
   int scalar_ode_option = 2;

   // Vis parameters
   const char * vishost = "localhost";
   int visport = 19916;


   // BRAID default parameters:
   int    max_levels  = 30;
   int    min_coarse  = 2;
   int    nrelax      = 1;
   int    nrelax0     = -1;
   double tol         = 1e-9;
   int    rtol        = 1;
   int    tnorm       = 2;
   int    skip        = 1;
   int    cfactor     = 2;
   int    cfactor0    = -1;
   int    max_iter    = 100;
   int    fmg         = 0;
   int    nfmg_Vcyc   = 1;
   int    access_level= 1;
   bool   wrapper_tests = false;
   bool   one_wrapper_test = false;

   // Use random initial vectors.
   int    init_rand   = 0;

   /* Parse command line */
   int print_usage = 0;
   int arg_index   = 1;
   while (arg_index < argc)
   {
      if (strcmp(argv[arg_index], "-mesh") == 0)
      {
         meshfile = argv[++arg_index];
      }
      else if (strcmp(argv[arg_index], "-sref") == 0)
      {
         sref = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-pref") == 0)
      {
         pref = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-nobraid") == 0)
      {
         use_braid = 0;
      }
      else if (strcmp(argv[arg_index], "-px") == 0)
      {
         num_procs_x = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-tstart") == 0)
      {
         tstart = atof(argv[++arg_index]);
      }
         else if (strcmp(argv[arg_index], "-kap") == 0)
      {
         kap = atof(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-tstop") == 0)
      {
         tstop = atof(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-nt") == 0)
      {
         ntime = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-ode") == 0)
      {
         arg_index++;
         if (strcmp(argv[arg_index], "heat") == 0)
            heat_equation = 1;
         else if (strcmp(argv[arg_index], "scalar") == 0)
         {
            heat_equation = 0;
            arg_index++;
            if (arg_index >= argc)
            {
               if (myid == 0)
                  cerr << "Missing option id after -ode scalar" << endl;
            }
            else
               scalar_ode_option = atoi(argv[arg_index]);
         }
         else
            if (myid == 0)
               cerr << "Unknown -ode parameter: " << argv[arg_index] << endl;
      }
      else if (strcmp(argv[arg_index], "-fe") == 0)
      {
         fe_degree = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-odesolver") == 0)
      {
         ode_solver = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-ml") == 0)
      {
         max_levels = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-mc") == 0)
      {
         min_coarse = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-nu") == 0)
      {
         nrelax = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-nu0") == 0)
      {
         nrelax0 = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-tol") == 0)
      {
         tol = atof(argv[++arg_index]);
      }
      else if( strcmp(argv[arg_index], "-rtol") == 0 )
      {
          rtol = atoi(argv[++arg_index]);
      }
      else if( strcmp(argv[arg_index], "-tnorm") == 0 )
      {
          tnorm = atoi(argv[++arg_index]);
      }
      else if( strcmp(argv[arg_index], "-skip") == 0 )
      {
          skip = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-cf") == 0)
      {
         cfactor = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-cf0") == 0)
      {
         cfactor0 = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-mi") == 0)
      {
         max_iter = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-fmg") == 0)
      {
         fmg = 1;
         nfmg_Vcyc = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-rand") == 0)
      {
         init_rand = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-wrapper_tests") == 0)
      {
         wrapper_tests = true;
      }
      else if (strcmp(argv[arg_index], "-one_wrapper_test") == 0)
      {
         one_wrapper_test = true;
      }
      else if (strcmp(argv[arg_index], "-access") == 0)
      {
         access_level = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-visport") == 0)
      {
         visport = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-vishost") == 0)
      {
         vishost = argv[++arg_index];
      }
      else if (strcmp(argv[arg_index], "-help") == 0)
      {
         print_usage = 1;
         break;
      }
      else
      {
         if (myid == 0)
            cerr << "Unknown parameter: " << argv[arg_index] << endl;
         MPI_Finalize();
         return 1;
      }
      arg_index++;
   }

   if (num_procs % num_procs_x)
   {
      if (myid == 0)
         cerr << "Total number of processors must be a multiple of the number"
            " of spatial processors!" << endl;
      MPI_Finalize();
      return 1;
   }

   if ((print_usage) && (myid == 0))
   {
      cout <<
         "\n"
         "Usage: " << argv[0] << " [<options>]\n"
         "\n"
         "  -mesh <file>      : spatial mesh (default: " << meshfile << ")\n"
         "  -wrapper_tests:   : quick run of the Braid wrapper tests\n"
         "                      (do not combine with temporal parallelism)\n"
         "  -one_wrapper_test : run only one wrapper test. can comment out/in\n"
         "                      the wrapper test that you want to focus on.\n"
         "                      (do not combine with temporal parallelism)\n"
         "  -sref <num>       : levels of serial refinements (default: 1)\n"
         "  -pref <num>       : levels of parallel refinements (default: 1)\n"
         "  -nobraid          : sequential time integration (default: off)\n"
         "  -px  <px>         : processors in space (default: 1)\n"
         "  -tstart <tstart>  : start time (default: 0.0)\n"
         "  -tstop  <tstop>   : stop time (default: 1.0)\n"
         "  -nt  <nt>         : number of time steps (default: 32)\n"
         "  -ode <name>       : ODE to solve, the options are:\n"
         "                         heat - Heat equation (default)\n"
         "                         scalar <id> - predefined scalar ODE\n"
         "  -fe <degree>      : FE degree/order (default: 1)\n"
         "  -odesolver <id>   : ODE solver id\n"
         "                      explicit methods:\n"
         "                         1 - Forward Euler\n"
         "                         2 - RK2 SSP\n"
         "                         3 - RK3 SSP\n"
         "                         4 - RK4\n"
         "                      implicit methods:\n"
         "                        -1  - Backward Euler (default)\n"
         "                        -2  - SDIRK(2,2)\n"
         "                        -3  - SDIRK(3,3)\n"
         "                        -21 - Implicit midpoint\n"
         "                        -31 - SDIRK(2,3)\n"
         "                        -41 - SDIRK(3,4)\n"
         "  -kap <kap>        : set the frequency of the exact solution sin wave (default: pi)\n"
         "  -ml  <max_levels> : set max number of time levels (default: 10)\n"
         "  -mc  <min_coarse> : set minimum possible coarse level size (default: 3)\n"
         "  -nu  <nrelax>     : set num F-C relaxations (default: 1)\n"
         "  -nu0 <nrelax>     : set num F-C relaxations on level 0\n"
         "  -tol <tol>        : set stopping tolerance (default: 1e-9)\n"
         "  -rtol <0/1>       : use relative or absolute stopping tolerance (default: 1)\n"
         "  -tnorm <tnorm>    : set temporal norm \n"
         "                      1 - One-norm \n"
         "                      2 - Two-norm (default) \n"
         "                      3 - Infinity-norm \n"
         "  -skip <skip>      : set skip, if true do no work on the first down cycle (default 1)\n"
         "  -cf  <cfactor>    : set coarsening factor (default: 2)\n"
         "  -cf0 <cfactor0>   : set aggressive coarsening (default: off)\n"
         "  -mi  <max_iter>   : set max iterations (default: 100)\n"
         "  -fmg <nfmg_Vcyc>  : use FMG cycling with nfmg_Vcyc V-cycles at each fmg level\n"
         "  -rand <0/1>       : use random initial vectors (default: 0)\n"
         "  -access  <a>      : set access_level (default: 1) \n"
         "  -vishost <vh>     : set glvis visualisation host (default: 'localhost') \n"
         "  -visport <vp>     : set glvis visualisation port (default: 19916) \n"
         "\n";
   }


   if (print_usage)
   {
      MPI_Finalize();
      return 0;
   }

   if (!use_braid)
      num_procs_x = num_procs;

   // Split comm into spatial and temporal communicators
   util.SplitCommworld(&comm, num_procs_x, &comm_x, &comm_t);
   MPI_Comm_rank(comm_x, &myid_x);
   MPI_Comm_rank(comm_t, &myid_t);

   Mesh *mesh;

   // Read the (serial) mesh from the given mesh file on all processors.
   ifstream imesh(meshfile);
   if (!imesh)
   {
      if (myid == 0)
         cerr << "\nCan not open mesh file: " << meshfile << '\n' << endl;
      MPI_Finalize();
      return 2;
   }
   mesh = new Mesh(imesh, 1, 1);
   imesh.close();

   // Refine the serial mesh on all processors to increase the resolution.
   for (int l = 0; l < sref; l++)
      mesh->UniformRefinement();

   if (mesh->NURBSext)
   {
      mesh->SetCurvature(fe_degree);
   }

   // Define a parallel mesh by a partitioning of the serial mesh. Refine this
   // mesh further in parallel to increase the resolution.
   ParMesh *pmesh = new ParMesh(comm_x, *mesh);
   delete mesh;
   for (int l = 0; l < pref; l++)
      pmesh->UniformRefinement();

   if (myid == 0)
      cout << endl;
   if (myid_t == 0)
      pmesh->PrintInfo();

   int dim = pmesh->Dimension();

   // Define a parallel finite element space on the parallel mesh. We use the
   // finite elements coming from the mesh nodes (linear by default).
   FiniteElementCollection *fec;
   fec = new H1_FECollection(fe_degree, dim);
   ParFiniteElementSpace *fespace = new ParFiniteElementSpace(pmesh, fec);
   int size = fespace->GlobalTrueVSize();
   if (myid == 0)
      cout << endl << "Number of unknowns: " << size << endl;

   // Define the initial condition vector x0 as a parallel finite element grid
   // function corresponding to fespace.
   ParGridFunction x0(fespace);
   FunctionCoefficient ic(&IC);
   x0.ProjectCoefficient(ic);
   HypreParVector *X0 = x0.ParallelAverage();

   // Set up the stiffness and mass parallel bilinear forms a(.,.) and m(.,.),
   // assemble and extract the corresponding parallel matrices.
   ConstantCoefficient one(1.0);
   ConstantCoefficient minus_one(-1.0);
   FunctionCoefficient source(&ST);
   VectorFunctionCoefficient nbc(dim, &NBC);
   FunctionCoefficient exact_sol(&ExSol);
   ParBilinearForm *a = new ParBilinearForm(fespace);
   ParBilinearForm *m = new ParBilinearForm(fespace);
   ParLinearForm   *b = new ParLinearForm(fespace);

   m->AddDomainIntegrator(new MassIntegrator(one));
   b->AddDomainIntegrator(new DomainLFIntegrator(source));
   b->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(nbc));

   // Local assembly. Use precomputed sparsity to ensure that 'a' and 'm' have
   // the same sparsity patterns.
   a->UsePrecomputedSparsity();
   m->UsePrecomputedSparsity(); m->Assemble(); m->Finalize();
   HypreParMatrix *M = m->ParallelAssemble();

   delete m;
   // Define the righ-hand side of the ODE and call MFEM or BRAID to solve it.
   TimeDependentOperator *ode;
   LinearSystemODE *diff_ode = NULL;
   if (heat_equation)
   {
      a->AddDomainIntegrator(new DiffusionIntegrator(minus_one));
      const int timedep_a = 0;
      diff_ode = new LinearSystemODE(a, timedep_a, M, b, &source, &nbc);
      ode = diff_ode;
   }
   else
   {
      ode = new ScalarODE(M->Height(), scalar_ode_option);
   }
   ODESolver *solver;

   if (ode_solver == 4)
      solver = new RK4Solver;
   else if (ode_solver == 3)
      solver = new RK3SSPSolver;
   else if (ode_solver == 2)
      solver = new RK2Solver(1.0);
   else if (ode_solver == 1)
      solver = new ForwardEulerSolver;
   else if (ode_solver == -1)
      solver = new BackwardEulerSolver;
   else if (ode_solver == -2)
      solver = new SDIRK23Solver(2);
   else if (ode_solver == -3)
      solver = new SDIRK33Solver;
   else if (ode_solver == -21)
      solver = new ImplicitMidpointSolver;
   else if (ode_solver == -31)
      solver = new SDIRK23Solver;
   else if (ode_solver == -41)
      solver = new SDIRK34Solver;
   else
   {
      if (myid == 0)
      {
         cerr << "Unknown ODE solver option: " << ode_solver
              << ", using default." << endl;
      }
      solver = new BackwardEulerSolver;
   }

   if (!use_braid)
   {
      if (heat_equation)
         SolveODE(ode, X0, solver, tstart, tstop, ntime, &x0, &exact_sol);
      else
         SolveODE(ode, X0, solver, tstart, tstop, ntime, &x0);
   }
   else
   {
      DiffBraidApp app(comm_t, diff_ode, X0, &x0, solver, tstart, tstop, ntime);
      if (init_rand)
         app.SetRandomInitVectors(285136749 + myid);
      app.SetVisHostAndPort(vishost, visport);

      if (wrapper_tests)
      {
         test_t = (app.tstop - app.tstart)/ (double) app.ntime;
         correct = util.TestAll(&app, comm_x, stdout, 0.0, test_t, 2*test_t);

         if(correct == 0)
         {
           cout << "Drive-04 Failed: at least one of the tests failed\n";
         }
      }
      else if(one_wrapper_test)
      {
         // Comment in/out the wrapper test that you want to focus on

         // Change the time value used by the test routines
         //test_t = app.tstart;
         test_t = app.tstop;
         test_dt = (app.tstop - app.tstart)/ (double) app.ntime;

         // These tests are cumulative.  You should test InitAccess before
         // Clone, Clone before Sum and so on.
         util.TestInitAccess( &app, comm_x, stdout, test_t);
         util.TestClone( &app, comm_x, stdout, test_t);
         util.TestSum( &app, comm_x, stdout, test_t);
         correct = util.TestSpatialNorm( &app, comm_x, stdout, test_t);
         correct = util.TestBuf( &app, comm_x, stdout, test_t);
         correct = util.TestCoarsenRefine(&app, comm_x, stdout,
                                     2*test_dt, test_dt/2, test_dt);

         if(correct == 0)
         {
           cout << "Drive-04 Failed: at least one of the tests failed\n";
         }
      }
      else
      {
         // Run a Braid simulation
         BraidCore core(comm, &app);
         if (heat_equation)
            app.SetExactSolution(&exact_sol);

         core.SetAccessLevel(access_level);
         core.SetPrintLevel(1);
         core.SetMaxLevels(max_levels);
         core.SetMinCoarse(min_coarse);
         core.SetNRelax(-1, nrelax);
         if (nrelax0 > -1)
            core.SetNRelax(0, nrelax0);
         rtol ? core.SetRelTol(tol) : core.SetAbsTol(tol);
         core.SetCFactor(-1, cfactor);
         core.SetAggCFactor(cfactor0);
         core.SetMaxIter(max_iter);
         core.SetTemporalNorm(tnorm);
         core.SetSkip(skip);
         if (fmg)
         {
            core.SetFMG();
            core.SetNFMGVcyc(nfmg_Vcyc);
         }
         core.Drive();
      }
   }

   // Free memory.
   delete solver;
   delete ode;

   delete X0;
   delete M;
   delete b;
   delete a;

   delete fespace;
   delete fec;
   delete pmesh;

   MPI_Comm_free(&comm_x);
   MPI_Comm_free(&comm_t);

   MPI_Finalize();

   return 0;
}
