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

#include "braid_mfem.hpp"
#include "mfem_arnoldi.hpp"
#include <fstream>
#include <iostream>

using namespace std;

// Choice for the problem setup. The fluid velocity, initial condition and
// inflow boundary condition are chosen based on this parameter.
int problem;

// Velocity coefficient
void velocity_function(const Vector &x, Vector &v);

// Initial condition
double u0_function(Vector &x);

// Inflow boundary condition
double inflow_function(Vector &x);


/** A time-dependent operator for the right-hand side of the ODE. The DG weak
    form of du/dt = v.grad(u) is M du/dt = K u + b, where M and K are the mass
    and advection matrices, and b describes the flow on the boundary. This can
    be written as a general ODE, du/dt = M^{-1} (K u + b), and this class is
    used to evaluate the right-hand side. */
class FE_Evolution : public TimeDependentOperator
{
private:
   HypreParMatrix &M, &K;
   const Vector &b;
   HypreSmoother M_prec;
   CGSolver M_solver;

   mutable Vector z;

public:
   FE_Evolution(HypreParMatrix &_M, HypreParMatrix &_K, const Vector &_b);

   virtual void Mult(const Vector &x, Vector &y) const;

   virtual ~FE_Evolution() { }
};


struct DGAdvectionOptions : public BraidOptions
{
   int    problem;
   double diffusion;
   int    order;
   int    ode_solver_type;
   bool   lump_mass;
   double dt; // derived from t_start, t_final, and num_time_steps
   bool   krylov_coarse;
   int    krylov_size;

   const char *vishost;
   int         visport;

   int vis_time_steps;
   int vis_braid_steps;
   bool vis_screenshots;

   DGAdvectionOptions(int argc, char *argv[]);
};


class DGAdvectionApp : public MFEMBraidApp
{
protected:
   DGAdvectionOptions &options;

   DG_FECollection fe_coll;

   VectorFunctionCoefficient velocity;
   FunctionCoefficient inflow;
   FunctionCoefficient u0;
   ConstantCoefficient diff;

   // Data for all discretization levels
   Array<HypreParMatrix *> M; // mass matrices
   Array<HypreParMatrix *> K; // advection matrices
   Array<HypreParVector *> B; // r.h.s. vectors

   // Allocate data structures for the given number of spatial levels. Used by
   // InitMultilevelApp.
   virtual void AllocLevels(int num_levels);

   // Construct the ParFiniteElmentSpace for the given mesh. Used by
   // InitMultilevelApp.
   virtual ParFiniteElementSpace *ConstructFESpace(ParMesh *pmesh);

   // Assuming mesh[l] and fe_space[l] are set, initialize, ode[l], solver[l],
   // and max_dt[l]. Used by InitMultilevelApp.
   virtual void InitLevel(int l);

public:
   DGAdvectionApp(DGAdvectionOptions &opts, MPI_Comm comm_t_, ParMesh *pmesh);
   
   SpaceTimeMeshInfo MeshInfo;

   // braid_Vector == BraidVector*
   virtual int Step(braid_Vector    u_,
                    braid_Vector    u_stop_,
                    braid_Vector    fstop_,
                    BraidStepStatus &pstatus);

   virtual ~DGAdvectionApp();
};


int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);

   MPI_Comm comm = MPI_COMM_WORLD;
   int myid;
   MPI_Comm_rank(comm, &myid);

   int precision = 8;
   cout.precision(precision);

   // Parse command line options
   DGAdvectionOptions opts(argc, argv);
   // Check for errors
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

   // Load the mesh and refine it in each processor (serial refinement)
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

   // If the mesh is NURBS, convert it to curved mesh
   if (mesh->NURBSext)
   {
      int dim = mesh->Dimension();
      int mesh_order = std::max(opts.order, 1);
      FiniteElementCollection *mfec = new H1_FECollection(mesh_order, dim);
      FiniteElementSpace *mfes = new FiniteElementSpace(mesh, mfec, dim);
      mesh->SetNodalFESpace(mfes);
      mesh->GetNodes()->MakeOwner(mfec);
   }

   // Split comm (MPI_COMM_WORLD) into spatial and temporal communicators
   MPI_Comm comm_x, comm_t;
   BraidUtil util;
   util.SplitCommworld(&comm, opts.num_procs_x, &comm_x, &comm_t);

   // Partition the mesh accross the spatial communicator
   ParMesh *pmesh = new ParMesh(comm_x, *mesh);
   delete mesh;

   // Create and initialize a DGAdvectionApp (including parallel refinement)
   DGAdvectionApp app(opts, comm_t, pmesh);

   // Run a Braid simulation
   BraidCore core(comm, &app);
   opts.SetBraidCoreOptions(core);

   core.Drive();
   app.MeshInfo.Print(comm);
   
   MPI_Finalize();
   return 0;
}


// Implementation of class FE_Evolution
FE_Evolution::FE_Evolution(HypreParMatrix &_M, HypreParMatrix &_K,
                           const Vector &_b)
   : TimeDependentOperator(_M.Height()),
     M(_M), K(_K), b(_b), M_solver(M.GetComm()), z(_M.Height())
{
   M_prec.SetType(HypreSmoother::Jacobi);
   M_solver.SetPreconditioner(M_prec);
   M_solver.SetOperator(M);

   M_solver.iterative_mode = false;
   M_solver.SetRelTol(1e-9);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(100);
   M_solver.SetPrintLevel(0);
}

void FE_Evolution::Mult(const Vector &x, Vector &y) const
{
   // y = M^{-1} (K x + b)
   K.Mult(x, z);
   z += b;
   M_solver.Mult(z, y);
}


DGAdvectionOptions::DGAdvectionOptions(int argc, char *argv[])
   : BraidOptions(argc, argv)
{
   // set defaults for inherited options
   t_start        = 0.0;
   t_final        = 10.0;
   num_time_steps = 1000;
   num_procs_x    = 1;

   // set defaults for the (optional) mesh/refinement inherited options
   mesh_file       = "../../mfem/data/periodic-hexagon.mesh";
   ser_ref_levels  = 2;
   par_ref_levels  = 0;
   AddMeshOptions();

   // set defaults for the DGAdvection-specific options
   problem         = 0;
   diffusion       = 0.0;
   order           = 3;
   ode_solver_type = 4;
   lump_mass       = false;
   krylov_coarse   = false;
   krylov_size     = 4;
   vishost         = "localhost";
   visport         = 19916;
   vis_time_steps  = 0;
   vis_braid_steps = 0;
   vis_screenshots = false;

   AddOption(&problem, "-p", "--problem",
             "Problem setup to use. See options in velocity_function().");
   AddOption(&diffusion, "-diff", "--diffusion", "Diffusion coefficient.");
   AddOption(&order, "-o", "--order",
             "Order (degree) of the finite elements.");
   AddOption(&ode_solver_type, "-s", "--ode-solver",
             "ODE solver: 1 - Forward Euler, 2 - RK2 SSP, 3 - RK3 SSP,"
             " 4 - RK4, 6 - RK6.");
   AddOption(&lump_mass, "-lump", "--lump-mass-matrix", "-dont-lump",
             "--dont-lump-mass-matrix", "Enable/disable lumping of the mass"
             " matrix.");
   AddOption(&krylov_coarse, "-kc", "--krylov-coarse", "-dont-kc",
             "--dont-krylov-coarse", "Enable/disable use of Arnoldi-based coarse-grid"
             " time-stepping.");
   AddOption(&krylov_size, "-ks", "--krylov-size", "Set size of the Krylov space for the"
                           " Arnoldi-based coarse-grid time-stepping.");
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


DGAdvectionApp::DGAdvectionApp(
      DGAdvectionOptions &opts, MPI_Comm comm_t_, ParMesh *pmesh)

   : MFEMBraidApp(comm_t_, opts.t_start, opts.t_final, opts.num_time_steps),
     options(opts),
     fe_coll(opts.order, pmesh->Dimension()),
     velocity(pmesh->Dimension(), velocity_function),
     inflow(inflow_function),
     u0(u0_function),
     diff(-opts.diffusion), // diffusion goes in the r.h.s. with a minus
     MeshInfo(opts.max_levels)
{
   problem = opts.problem;

   InitMultilevelApp(pmesh, opts.par_ref_levels, opts.spatial_coarsen);

   // Set the initial condition
   x[0]->ProjectCoefficient(u0);
   HypreParVector *U = x[0]->GetTrueDofs();
   SetInitialCondition(U);

   SetVisHostAndPort(opts.vishost, opts.visport);
   SetVisSampling(opts.vis_time_steps, opts.vis_braid_steps);
   SetVisScreenshots(opts.vis_screenshots);
}

DGAdvectionApp::~DGAdvectionApp()
{
   for (int l = 0; l < GetNumSpaceLevels(); l++)
   {
      delete B[l];
      delete K[l];
      delete M[l];
   }
}


int DGAdvectionApp::Step(braid_Vector    u_,
                         braid_Vector    ustop_,
                         braid_Vector    fstop_,
                         BraidStepStatus &pstatus)
{
   // This contains one small change over the default Step, we store the Space-Time mesh info
   
   // Store Space-Time Mesh Info
   BraidVector *u     = (BraidVector*) u_;
   
 //BraidVector *start  = new BraidVector(0, *X0);
 //braid_Vector start2 = (braid_Vector) start;

   int level = u->level;
   double tstart, tstop, dt;
   int braid_level, iter;
   
   pstatus.GetTstartTstop(&tstart, &tstop);
   pstatus.GetLevel(&braid_level);
   pstatus.GetIter(&iter);
   dt = tstop - tstart;
 
 // This causes some kind of bug sometimes...
 //if(iter > 2)
 //{
 //   *start = 0.0;
 //   Sum(1.0, ustop_, 0.0, start2);
 //   Sum(iter*(1 - 0.1), u_, iter*0.1, start2);
 //}
 //else
 //{
 //   *start = 0.0;
 //   Sum(1.0, u_, 0.0, start2);
 //}

   MeshInfo.SetRow(braid_level, level, x[u->level]->ParFESpace()->GetParMesh(), dt);
   
   if((braid_level == 0) || (options.krylov_coarse == false) )
   {
       // Call the default Step
       MFEMBraidApp::Step(u_,ustop_,fstop_, pstatus);
   }
   else
   {
       double norm_u_sq = InnerProduct(u, u);
       
       if(norm_u_sq > 0.0)
       {
           // Generate Krylov Space
           Arnoldi arn(options.krylov_size*braid_level, mesh[0]->GetComm());
           //Arnoldi arn(options.krylov_size, mesh[0]->GetComm());
           arn.SetOperator(*ode[0]);
           arn.GenKrylovSpace(*u);
           Vector ubar;
           arn.ApplyVT(*u, ubar);
           
           // Setup the ODE solver
           ODESolver *ode_solver = NULL;
           switch (options.ode_solver_type)
           {
              case 1: ode_solver = new ForwardEulerSolver; break;
              case 2: ode_solver = new RK2Solver(1.0); break;
              case 3: ode_solver = new RK3SSPSolver; break;
              case 4: ode_solver = new RK4Solver; break;
              case 6: ode_solver = new RK6Solver; break;
           
              default: ode_solver = new RK4Solver; break;
           }
           ode_solver->Init(arn.GetH());
           
           // Do time-stepping
           int m;
           if( options.cfactor0 == -1)
           {
               m = pow(options.cfactor, braid_level);
           }
           else
           {
               m = options.cfactor0*pow(options.cfactor, braid_level-1);
           }
           
           double dt_fine = dt / ( (double) m);
           for(int k = 0; k < m; k++)
           {
               // Step() advances tstart by dt_fine
               ode_solver->Step(ubar, tstart, dt_fine);
           }
           
           // Convert ubar back to full space
           arn.ApplyV(ubar, *u);
           
           delete ode_solver;
       }

       // no refinement
       pstatus.SetRFactor(1);
       
   }

   return 0;
}

// Allocate data structures for the given number of spatial levels. Used by
// InitMultilevelApp.
void DGAdvectionApp::AllocLevels(int num_levels)
{
   M.SetSize(num_levels, NULL);
   K.SetSize(num_levels, NULL);
   B.SetSize(num_levels, NULL);
}

// Construct the ParFiniteElmentSpace for the given mesh. Used by
// InitMultilevelApp.
ParFiniteElementSpace *DGAdvectionApp::ConstructFESpace(ParMesh *pmesh)
{
   return new ParFiniteElementSpace(pmesh, &fe_coll);
}

// Assuming mesh[l] and fe_space[l] are set, initialize, ode[l], solver[l],
// and max_dt[l]. Used by InitMultilevelApp.
void DGAdvectionApp::InitLevel(int l)
{
   // Create the components needed for the time-dependent operator, ode[l]:
   // the mass matrix, M; the advection (+diffusion) matrix, K; and the inflow
   // b.c. vector, B.
   ParBilinearForm *m = new ParBilinearForm(fe_space[l]);
   if (!options.lump_mass)
      m->AddDomainIntegrator(new MassIntegrator);
   else
      m->AddDomainIntegrator(new LumpedIntegrator(new MassIntegrator));
   ParBilinearForm *k = new ParBilinearForm(fe_space[l]);
   k->AddDomainIntegrator(new ConvectionIntegrator(velocity, -1.0));
   k->AddInteriorFaceIntegrator(
      new TransposeIntegrator(new DGTraceIntegrator(velocity, 1.0, -0.5)));
   k->AddBdrFaceIntegrator(
      new TransposeIntegrator(new DGTraceIntegrator(velocity, 1.0, -0.5)));
   if (options.diffusion > 0.0)
   {
      double sigma, kappa;
      if (1)
      {
         // IP method
         sigma = -1.0;
         kappa = (options.order+1)*(options.order+1);
      }
      else
      {
         // NIPG method
         sigma = 1.0;
         kappa = 1.0;
      }
      k->AddDomainIntegrator(new DiffusionIntegrator(diff));
      k->AddInteriorFaceIntegrator(
         new DGDiffusionIntegrator(diff, sigma, kappa));
   }

   ParLinearForm *b = new ParLinearForm(fe_space[l]);
   b->AddBdrFaceIntegrator(
      new BoundaryFlowIntegrator(inflow, velocity, -1.0, -0.5));

   m->Assemble();
   m->Finalize();
   int skip_zeros = 0;
   k->Assemble(skip_zeros);
   k->Finalize(skip_zeros);
   b->Assemble();

   M[l] = m->ParallelAssemble();
   K[l] = k->ParallelAssemble();
   B[l] = b->ParallelAssemble();

   delete b;
   delete k;
   delete m;

   // Create the time-dependent operator, ode[l]
   ode[l] = new FE_Evolution(*M[l], *K[l], *B[l]);

   // Setup the ODE solver, solver[l]
   ODESolver *ode_solver = NULL;
   switch (options.ode_solver_type)
   {
      case 1: ode_solver = new ForwardEulerSolver; break;
      case 2: ode_solver = new RK2Solver(1.0); break;
      case 3: ode_solver = new RK3SSPSolver; break;
      case 4: ode_solver = new RK4Solver; break;
      case 6: ode_solver = new RK6Solver; break;

      default: ode_solver = new RK4Solver; break;
   }
   solver[l] = ode_solver;

   // Associate the ODE operator, ode[l], with the ODE solver, solver[l]
   solver[l]->Init(*ode[l]);

   // Set max_dt[l] = 1.01*dt[0]*(2^l)
   max_dt[l] = 1.01 * options.dt * (1 << (l+1));
  //max_dt[l] = 1.01 * options.dt * (1 << l);
}


// Velocity coefficient
void velocity_function(const Vector &x, Vector &v)
{
   int dim = x.Size();

   switch (problem)
   {
      case 0:
      {
         // Translations in 1D, 2D, and 3D
         switch (dim)
         {
            case 1: v(0) = 1.0; break;
            case 2: v(0) = sqrt(2./3.); v(1) = sqrt(1./3.); break;
            case 3: v(0) = sqrt(3./6.); v(1) = sqrt(2./6.); v(2) = sqrt(1./6.);
               break;
         }
         break;
      }
      case 1:
      case 2:
      {
         // Clockwise rotation in 2D around the origin
         const double w = M_PI/2;
         switch (dim)
         {
            case 1: v(0) = 1.0; break;
            case 2: v(0) = w*x(1); v(1) = -w*x(0); break;
            case 3: v(0) = w*x(1); v(1) = -w*x(0); v(2) = 0.0; break;
         }
         break;
      }
      case 3:
      {
         // Clockwise twisting rotation in 2D around the origin
         const double w = M_PI/2;
         double d = max((x(0)+1.)*(1.-x(0)),0.) * max((x(1)+1.)*(1.-x(1)),0.);
         d = d*d;
         switch (dim)
         {
            case 1: v(0) = 1.0; break;
            case 2: v(0) = d*w*x(1); v(1) = -d*w*x(0); break;
            case 3: v(0) = d*w*x(1); v(1) = -d*w*x(0); v(2) = 0.0; break;
         }
         break;
      }
   }
}

// Initial condition
double u0_function(Vector &x)
{
   int dim = x.Size();

   switch (problem)
   {
      case 0:
      case 1:
      {
         switch (dim)
         {
            case 1:
               return exp(-40.*pow(x(0)-0.5,2));
            case 2:
            case 3:
            {
               double rx = 0.45, ry = 0.25, cx = 0., cy = -0.2, w = 10.;
               if (dim == 3)
               {
                  const double s = (1. + 0.25*cos(2*M_PI*x(2)));
                  rx *= s;
                  ry *= s;
               }
               return ( erfc(w*(x(0)-cx-rx))*erfc(-w*(x(0)-cx+rx)) *
                        erfc(w*(x(1)-cy-ry))*erfc(-w*(x(1)-cy+ry)) )/16;
            }
         }
      }
      case 2:
      {
         const double r = sqrt(8.);
         double x_ = x(0), y_ = x(1), rho, phi;
         rho = hypot(x_, y_) / r;
         phi = atan2(y_, x_);
         return pow(sin(M_PI*rho),2)*sin(3*phi);
      }
      case 3:
      {
         const double f = M_PI;
         return sin(f*x(0))*sin(f*x(1));
      }
   }
   return 0.0;
}

// Inflow boundary condition (zero for the problems considered in this example)
double inflow_function(Vector &x)
{
   switch (problem)
   {
      case 0:
      case 1:
      case 2:
      case 3: return 0.0;
   }
   return 0.0;
}
