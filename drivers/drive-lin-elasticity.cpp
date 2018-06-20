// Copyright (c) 2013, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. Written by
// Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin
// Dobrev, et al. LLNL-CODE-660355. All rights reserved.
//
// This file is part of XBraid. For support, post issues to the XBraid Github page.
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


// Driver:        drive-lin-elasticity.cpp
//
// Interface:     C++, through MFEM 
// 
// Requires:      MFEM, Hypre, Metis and GlVis
//                Modify Makefile to point to metis, mfem and hypre libraries
//
// Compile with:  make drive-lin-elasticity
//
// Sample run:    mpirun -np 4 ./drive-lin-elasticity -cf 4 -nu 2 -skip 0 -no-vis
//
// Help with:     drive-lin-elasticity -help
//
// Description:   Solves time-dependent linearized elasticity 
//                
//                After spatial discretization, the elasticity model can be written 
//                as a system of ODEs:
//               
//                    dv/dt = -M^{-1}*(H*x + S*v)
//                    dx/dt = v,
//                
//                where x is the vector representing the displacement, v is the
//                velocity field, M is the mass matrix, S is the viscosity matrix, 
//                and H is the elasticity matrix.


#include "braid_mfem.hpp"
#include <memory>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;
using namespace mfem;

const double ref_density = 1.0; // density in the reference configuration
void InitialDeformation(const Vector &x, Vector &y);
void InitialVelocity(const Vector &x, Vector &v);


//  Class LinearElasticOperator represents the right-hand side of the above
//  system of ODEs in terms of dv/dt and dx/dt 
class LinearElasticOperator : public TimeDependentOperator
{
protected:
   ParFiniteElementSpace &fespace;

   mutable ParBilinearForm M;
   ParBilinearForm S, H;
   double viscosity;

   mutable HypreParMatrix *Mmat; // Mass matrix from ParallelAssemble()
   mutable CGSolver M_solver;    // CG solver for inverting the mass matrix M
   HypreSmoother M_prec;         // Preconditioner for the mass matrix M

   mutable Vector z; // auxiliary vector

   // Data used by ImplicitSolve.
   double current_dt;
   Array<double> J_dt;
   mutable Array<HypreParMatrix*> J;
   Array<int> ess_tdof_list;
   mutable Vector w;

   // Preconditioners for ImplicitSolve.
   Array<HypreBoomerAMG*> J_amg;

public:
   LinearElasticOperator(ParFiniteElementSpace &f, Array<int> &ess_bdr,
                         double visc, double mu, double K);

   virtual void Mult(const Vector &vx, Vector &dvx_dt) const;
   /** Solve the Backward-Euler equation: k = f(x + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k);

   double ElasticEnergy(ParGridFunction &x) const;
   double KineticEnergy(ParGridFunction &v) const;

   double Norm(const Vector &vx) const;

   virtual ~LinearElasticOperator();
};


struct LinearElasticityOptions : public BraidOptions
{
   int    order;
   int    ode_solver_type;
   double visc;
   double mu;
   double lambda;
   bool   visualization;
   int    vis_steps;
   bool   screenshots;

   // derived options
   double dt;

   LinearElasticityOptions(int argc, char *argv[]);

   void Parse();
};


class LinearElasticityApp : public BraidApp
{
public:
   enum ErrorType
   {
      E_NONE      = 0,
      E_ODESolver = 1
   };

protected:
   ErrorType                error;
   LinearElasticityOptions &opts;
   ParMesh                 &pmesh;
   H1_FECollection          fe_coll;
   ParFiniteElementSpace    fespace;
   ODESolver               *ode_solver;
   Array<int>               true_offsets;
   BlockVector             *vx;
   ParGridFunction          v_gf;
   ParGridFunction          x_gf;
   ParGridFunction          x_ref;
   Array<int>               ess_bdr;
   LinearElasticOperator   *ode_oper;
   socketstream             vis;

   double ee0, ke0;

public:
   LinearElasticityApp(LinearElasticityOptions &opts_,
                       MPI_Comm                 comm_t_,
                       ParMesh                 &pmesh_);

   ErrorType GetError() { return error; }

   void PrintInfo(bool root);

   virtual ~LinearElasticityApp();

   // ---------- Virtual Braid methods ----------
   // braid_Vector = BlockVector *

   virtual int Step(
         braid_Vector     u_,
         braid_Vector     ustop_,
         braid_Vector     fstop_,
         BraidStepStatus &pstatus);

   virtual int Residual(
         braid_Vector     u_,
         braid_Vector     r_,
         BraidStepStatus &pstatus);

   virtual int Clone(
         braid_Vector  u_,
         braid_Vector *v_ptr);

   virtual int Init(
         double        t,
         braid_Vector *u_ptr);

   virtual int Free(
         braid_Vector u_);

   virtual int Sum(
         double       alpha,
         braid_Vector x_,
         double       beta,
         braid_Vector y_);

   virtual int SpatialNorm(
         braid_Vector  u_,
         double       *norm_ptr);

   virtual int BufSize(
         int               *size_ptr,
         BraidBufferStatus &status);

   virtual int BufPack(
         braid_Vector       u_,
         void              *buffer,
         BraidBufferStatus &status);

   virtual int BufUnpack(
         void              *buffer,
         braid_Vector      *u_ptr,
         BraidBufferStatus &status);

   // optional
   virtual int Coarsen(
         braid_Vector           fu_,
         braid_Vector          *cu_ptr,
         BraidCoarsenRefStatus &status);

   // optional
   virtual int Refine(
         braid_Vector           cu_,
         braid_Vector          *fu_ptr,
         BraidCoarsenRefStatus &status);

   virtual int Access(
         braid_Vector       u_,
         BraidAccessStatus &astatus);
};


int main(int argc, char *argv[])
{
   // Initialize MPI. Automatically finalizes MPI at destruction.
   MPI_Session mpi(argc, argv);

   // Command-line options parser, extended with Braid and linear elasticity
   // options.
   LinearElasticityOptions opts(argc, argv);

   // Parse command-line options.
   opts.Parse();
   if (!opts.Good())
   {
      if (mpi.Root())
      {
         opts.PrintUsage(cout);
      }
      return 1;
   }
   if (mpi.Root())
   {
      opts.PrintOptions(cout);
   }

   // Split MPI_COMM_WORLD into spatial, 'x', and time, 't', components.
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Comm comm_x, comm_t;
   BraidUtil util;
   util.SplitCommworld(&comm, opts.num_procs_x, &comm_x, &comm_t);

   // Read the mesh and refine it as specified on the command line.
   ParMesh *pmesh = opts.LoadMeshAndRefine(comm_x);
   if (!pmesh)
   {
      if (mpi.Root())
      {
         cout << "Error loading mesh file: " << opts.mesh_file << endl;
      }
      return 2;
   }

   // Create the main linear elasticity 'app'.
   LinearElasticityApp app(opts, comm_t, *pmesh);

   // Check for errors, currently just one error is possible E_ODESolver.
   if (app.GetError())
   {
      if (mpi.Root())
      {
         cout << "Unknown ODE solver type: " << opts.ode_solver_type << '\n';
      }
      delete pmesh;
      return 3;
   }

   app.PrintInfo(mpi.Root());

   // Run a Braid simulation
   BraidCore core(comm, &app);
   opts.SetBraidCoreOptions(core);

   core.Drive();

   // Free the used memory.
   delete pmesh;

   return 0;
}


void visualize(ostream &out, ParMesh *mesh, ParGridFunction *displacement,
               ParGridFunction *field, const char *field_name = NULL,
               bool init_vis = false)
{
   if (!out) { return; }

   ParGridFunction x_new(displacement->ParFESpace());
   mesh->GetNodes(x_new);
   x_new += *displacement;
   GridFunction *nodes = &x_new;
   int owns_nodes = 0;

   mesh->SwapNodes(nodes, owns_nodes);

   out << "parallel " << mesh->GetNRanks() << " " << mesh->GetMyRank() << "\n";
   out << "solution\n" << *mesh << *field;

   mesh->SwapNodes(nodes, owns_nodes);

   if (init_vis)
   {
      out << "window_size 800 800\n";
      out << "window_title '" << field_name << "'\n";
      if (mesh->SpaceDimension() == 2)
      {
         out << "view 0 0\n"; // view from top
         out << "keys jl\n";  // turn off perspective and light
      }
      out << "keys cmA\n";         // show colorbar and mesh; multisampling on
      out << "autoscale value\n"; // update value-range; keep mesh-extents fixed
      // out << "pause\n";
   }
   out << flush;
}


LinearElasticOperator::LinearElasticOperator(ParFiniteElementSpace &f,
                                             Array<int> &ess_bdr, double visc,
                                             double lambda, double mu)
   : TimeDependentOperator(2*f.TrueVSize(), 0.0), fespace(f),
     M(&fespace), S(&fespace), H(&fespace), Mmat(NULL), M_solver(f.GetComm()),
     z(height/2), current_dt(0.0)
{
   const double rel_tol = 1e-8;
   const int skip_zero_entries = 0;

   ConstantCoefficient rho0(ref_density);
   M.AddDomainIntegrator(new VectorMassIntegrator(rho0));
   M.Assemble(skip_zero_entries);
   M.EliminateEssentialBC(ess_bdr);
   M.Finalize(skip_zero_entries);

   M_solver.iterative_mode = false;
   M_solver.SetRelTol(rel_tol);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(30);
   M_solver.SetPrintLevel(0);
   M_prec.SetType(HypreSmoother::Jacobi);
   M_solver.SetPreconditioner(M_prec);

   ConstantCoefficient lambda_c(lambda);
   ConstantCoefficient mu_c(mu);
   H.AddDomainIntegrator(new ElasticityIntegrator(lambda_c, mu_c));
   H.Assemble(skip_zero_entries);
   H.EliminateEssentialBC(ess_bdr);
   H.Finalize(skip_zero_entries);

   viscosity = visc;
   ConstantCoefficient visc_coeff(viscosity);
   S.AddDomainIntegrator(new VectorDiffusionIntegrator(visc_coeff));
   S.Assemble(skip_zero_entries);
   S.EliminateEssentialBC(ess_bdr);
   S.Finalize(skip_zero_entries);

   // ess_tdof_list is needed for ImplicitSolve
   fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
}

void LinearElasticOperator::Mult(const Vector &vx, Vector &dvx_dt) const
{
   // Create views to the sub-vectors v, x of vx, and dv_dt, dx_dt of dvx_dt
   // Here, v is velocaity, and x is displacement.
   int sc = height/2;
   Vector v(vx.GetData() +  0, sc);
   Vector x(vx.GetData() + sc, sc);
   Vector dv_dt(dvx_dt.GetData() +  0, sc);
   Vector dx_dt(dvx_dt.GetData() + sc, sc);

   if (Mmat == NULL)
   {
      Mmat = M.ParallelAssemble();
      M_solver.SetOperator(*Mmat);
   }

   H.TrueAddMult(x, z);
   if (viscosity != 0.0)
   {
      S.TrueAddMult(v, z);
   }
   z.Neg(); // z = -z
   M_solver.Mult(z, dv_dt);

   dx_dt = v;
}

void LinearElasticOperator::ImplicitSolve(const double dt,
                                          const Vector &vx, Vector &dvx_dt)
{
   // Create views to the sub-vectors v, x of vx, and dv_dt, dx_dt of dvx_dt
   // Here, v is velocaity, and x is displacement.
   int sc = height/2;
   Vector v(vx.GetData() +  0, sc);
   Vector x(vx.GetData() + sc, sc);
   Vector dv_dt(dvx_dt.GetData() +  0, sc);
   Vector dx_dt(dvx_dt.GetData() + sc, sc);
   // Allocate data in 'w', if not already allocated.
   w.SetSize(sc);

   // Setup the matrix of the system we need to solve:
   //    J = M + dt*S + dt*dt*H.
   int j = -1;
   for (int i = 0; i < J_dt.Size(); i++)
   {
      if (J_dt[i] == dt) { j = i; break; }
   }
   if (j == -1)
   {
      // cout << "LinearElasticOperator: add dt = " << dt << endl;
      J_dt.Append(dt);
      SparseMatrix *localJ = Add(1.0, M.SpMat(), dt, S.SpMat());
      localJ->Add(dt*dt, H.SpMat());
      J.Append(M.ParallelAssemble(localJ));
      delete localJ;
      HypreParMatrix *Je = J.Last()->EliminateRowsCols(ess_tdof_list);
      delete Je;

      J_amg.Append(new HypreBoomerAMG(*J.Last()));
      J_amg.Last()->SetPrintLevel(0);
      J_amg.Last()->SetSystemsOptions(fespace.GetMesh()->Dimension());
      j = J.Size() - 1;
   }
   if (current_dt != dt)
   {
      if (fespace.GetMyRank() == 0)
      {
         cout << "ImplicitSolve: " << current_dt << " --> " << dt << endl;
      }
      current_dt = dt;
   }

   // By eliminating kx from the coupled system:
   //    kv = -M^{-1}*[H*(x + dt*kx) + S*(v + dt*kv)]
   //    kx = v + dt*kv
   // we reduce it to a linear equation for kv, i.e. dv_dt.
   add(x, dt, v, z);    // z = x + dt*v
   for (int i = 0; i < ess_tdof_list.Size(); i++)
   {
      z(ess_tdof_list[i]) = 0.0;
   }
   w = 0.;
   H.TrueAddMult(z, w);  // w = H*(x + dt*v)
   S.TrueAddMult(v, w);  // w = H*(x + dt*v) + S*v
   w.Neg();              // w = -[H*(x + dt*v) + S*v]
   for (int i = 0; i < ess_tdof_list.Size(); i++)
   {
      w(ess_tdof_list[i]) = 0.0;
   }

   CGSolver J_pcg(fespace.GetComm());
   J_pcg.iterative_mode = false;
   J_pcg.SetRelTol(1e-8);
   J_pcg.SetMaxIter(500);
   J_pcg.SetOperator(*J[j]);
   J_pcg.SetPreconditioner(*J_amg[j]);
   J_pcg.Mult(w, dv_dt);

   add(v, dt, dv_dt, dx_dt);
}

double LinearElasticOperator::ElasticEnergy(ParGridFunction &x) const
{
   double loc_energy = 0.5*H.InnerProduct(x, x);
   double energy;
   MPI_Allreduce(&loc_energy, &energy, 1, MPI_DOUBLE, MPI_SUM,
                 fespace.GetComm());
   return energy;
}

double LinearElasticOperator::KineticEnergy(ParGridFunction &v) const
{
   double loc_energy = 0.5*M.InnerProduct(v, v);
   double energy;
   MPI_Allreduce(&loc_energy, &energy, 1, MPI_DOUBLE, MPI_SUM,
                 fespace.GetComm());
   return energy;
}

double LinearElasticOperator::Norm(const Vector &vx) const
{
   int sc = height/2;
   Vector v(vx.GetData() +  0, sc);
   Vector x(vx.GetData() + sc, sc);

   for (int i = 0; i < ess_tdof_list.Size(); i++)
   {
      MFEM_VERIFY(v(ess_tdof_list[i]) == 0.0,
                  "expect zero velocity at essential boundary");
      MFEM_VERIFY(x(ess_tdof_list[i]) == 0.0,
                  "expect zero displacement at essential boundary");
   }

   double local = M.InnerProduct(v, v) + M.InnerProduct(x, x);
   double global;
   MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, fespace.GetComm());
   return sqrt(global);
}

LinearElasticOperator::~LinearElasticOperator()
{
   if (fespace.GetMyRank() == 0)
   {
      cout << "LinearElasticOperator: num dts = " << J.Size() << endl;
   }
   for (int i = J.Size()-1; i >= 0; i--)
   {
      delete J[i];
      delete J_amg[i];
   }
   delete Mmat;
}


LinearElasticityOptions::LinearElasticityOptions(int argc, char *argv[])
   : BraidOptions(argc, argv)
{
   // change derived default values
   t_final         = 300.0;

   // add mesh + refinement options
   AddMeshOptions();
   mesh_file      = "../../mfem/data/beam-quad.mesh";
   ser_ref_levels = 2;
   par_ref_levels = 0;

   // linear elasticity options
   order           = 2;
   ode_solver_type = 3;
   visc            = 1e-2;
   mu              = 0.25;
   lambda          = 5.-1./6; // K-(2/3)*mu, with K=5
   visualization   = true;
   vis_steps       = 1;
   screenshots     = false;

   AddOption(&order, "-o", "--order",
             "Order (degree) of the finite elements.");
   AddOption(&ode_solver_type, "-s", "--ode-solver",
             "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3,\n\t"
             "\t   11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.");
   AddOption(&visc, "-v", "--viscosity", "Viscosity coefficient.");
   AddOption(&mu, "-mu", "--shear-modulus",
             "Shear modulus in the elastic model.");
   AddOption(&lambda, "-l", "--lambda",
             "Lame's first parameter in the elastic model.");
   AddOption(&visualization, "-vis", "--visualization", "-no-vis",
             "--no-visualization",
             "Enable or disable GLVis visualization.");
   AddOption(&vis_steps, "-vs", "--visualization-steps",
             "Visualize every n-th timestep.");
   AddOption(&screenshots, "-ss", "--take-screenshots", "-no-ss",
             "--dont-take-screenshots",
             "Take or not screenshots of the GLVis windows.");

   // derived options
   // dt = (t_final - t_start)/num_time_steps;
}

void LinearElasticityOptions::Parse()
{
   BraidOptions::Parse();
   dt = (t_final - t_start)/num_time_steps;
}


LinearElasticityApp::LinearElasticityApp(
      LinearElasticityOptions &opts_,
      MPI_Comm                 comm_t_,
      ParMesh                 &pmesh_)

   : BraidApp(comm_t_,
              opts_.t_start,
              opts_.t_final,
              opts_.num_time_steps),
     error(E_NONE),
     opts(opts_),
     pmesh(pmesh_),
     fe_coll(opts.order,
             pmesh.Dimension()),
     fespace(&pmesh,
             &fe_coll,
             pmesh.Dimension()),
     ode_solver(NULL),
     true_offsets(3),
     vx(NULL),
     v_gf(&fespace),
     x_gf(&fespace),
     x_ref(&fespace),
     ode_oper(NULL),
     ee0(0.0),
     ke0(0.0)
{
   // Define the ODE solver used for time integration.
   switch (opts.ode_solver_type)
   {
      // Implicit L-stable methods
      case 1:  ode_solver = new BackwardEulerSolver; break;
      case 2:  ode_solver = new SDIRK23Solver(2); break;
      case 3:  ode_solver = new SDIRK33Solver; break;
      // Explicit methods
      case 11: ode_solver = new ForwardEulerSolver; break;
      case 12: ode_solver = new RK2Solver(0.5); break; // midpoint method
      case 13: ode_solver = new RK3SSPSolver; break;
      case 14: ode_solver = new RK4Solver; break;
      // Implicit A-stable methods (not L-stable)
      case 22: ode_solver = new ImplicitMidpointSolver; break;
      case 23: ode_solver = new SDIRK23Solver; break;
      case 24: ode_solver = new SDIRK34Solver; break;
      default: ode_solver = NULL; error = E_ODESolver; return;
   }

   int true_size = fespace.TrueVSize();
   true_offsets[0] = 0;
   true_offsets[1] = true_size;
   true_offsets[2] = 2*true_size;

   vx = new BlockVector(true_offsets);

   // Get the initial mesh configuration.
   pmesh.GetNodes(x_ref);

   // Set the initial conditions for v_gf, x_gf and vx, and define the
   // boundary conditions on a beam-like mesh (see MFEM's example 10).
   VectorFunctionCoefficient velo(pmesh.Dimension(), InitialVelocity);
   v_gf.ProjectCoefficient(velo);
   VectorFunctionCoefficient deform(pmesh.Dimension(), InitialDeformation);
   x_gf.ProjectCoefficient(deform);
   x_gf -= x_ref;

   v_gf.GetTrueDofs(vx->GetBlock(0));
   x_gf.GetTrueDofs(vx->GetBlock(1));

   ess_bdr.SetSize(fespace.GetMesh()->bdr_attributes.Max(), 0);
   ess_bdr[0] = 1; // boundary attribute 1 (index 0) is fixed

   ode_oper = new LinearElasticOperator(fespace, ess_bdr, opts.visc,
                                        opts.lambda, opts.mu);
   ode_solver->Init(*ode_oper);

   v_gf.Distribute(vx->GetBlock(0));
   x_gf.Distribute(vx->GetBlock(1));
   ee0 = ode_oper->ElasticEnergy(x_gf);
   ke0 = ode_oper->KineticEnergy(v_gf);

   if (opts.visualization)
   {
      vis.open("localhost", 19916);
      visualize(vis, &pmesh, &x_gf, &v_gf, "initial velocity", true);
   }
}

void LinearElasticityApp::PrintInfo(bool root)
{
   // P-wave speed
   double vp = sqrt((opts.lambda + 2*opts.mu)/ref_density);
   HYPRE_Int glob_size = fespace.GlobalTrueVSize();
   if (root)
   {
      cout <<
         "Number of velocity/displacement unknowns: " << glob_size << "\n"
         "P-wave speed = " << vp << "\n\n"
         "initial elastic energy (EE) = " << ee0 << "\n"
         "initial kinetic energy (KE) = " << ke0 << "\n"
         "initial   total energy (TE) = " << (ee0 + ke0) << endl;
   }
}

LinearElasticityApp::~LinearElasticityApp()
{
   delete ode_oper;
   delete vx;
   delete ode_solver;
}

int LinearElasticityApp::Step(
      braid_Vector     u_,
      braid_Vector     ustop_,
      braid_Vector     fstop_,
      BraidStepStatus &pstatus)
{
   BlockVector *u     = (BlockVector *) u_;
   // BlockVector *ustop = (BlockVector *) ustop_;
   BlockVector *fstop = (BlockVector *) fstop_;

   MFEM_VERIFY(fstop == NULL, "TODO");

   double tstart, tstop, dt = opts.dt;
   pstatus.GetTstartTstop(&tstart, &tstop);
   double mult = (tstop - tstart)/dt;
   MFEM_VERIFY(abs(round(mult)-mult) < 1e-12, "mult = " << mult);
   dt *= round(mult);

   ode_solver->Step(*u, tstart, dt);

   pstatus.SetRFactor(1); // no refinement. optional?
   return 0;
}

int LinearElasticityApp::Residual(
      braid_Vector     u_,
      braid_Vector     r_,
      BraidStepStatus &pstatus)
{
   // BlockVector *u = (BlockVector *) u_; // input
   // BlockVector *r = (BlockVector *) r_; // input,output
   // u - input: approximate solution at tstop
   // r - input: approximate solution at tstart
   //   - output: residual at tstop
   // tstart, tstop - from pstatus (not needed here)
   cout << _MFEM_FUNC_NAME << ": not implemented" << endl;
   return 0;
}

int LinearElasticityApp::Clone(
      braid_Vector  u_,
      braid_Vector *v_ptr)
{
   BlockVector *u = (BlockVector *) u_;
   BlockVector *v = new BlockVector(*u); // copy data
   *v_ptr = (braid_Vector) v;
   return 0;
}

int LinearElasticityApp::Init(
      double        t,
      braid_Vector *u_ptr)
{
   BlockVector *u = new BlockVector(true_offsets);
   if (t == tstart)
   {
      *u = *vx; // initial condition
   }
   else
   {
      *u = 0.0;
   }
   *u_ptr = (braid_Vector) u;
   return 0;
}

int LinearElasticityApp::Free(
      braid_Vector u_)
{
   BlockVector *u = (BlockVector *) u_;
   delete u;
   return 0;
}

int LinearElasticityApp::Sum(
      double       alpha,
      braid_Vector x_,
      double       beta,
      braid_Vector y_)
{
   BlockVector *x = (BlockVector *) x_;
   BlockVector *y = (BlockVector *) y_;
   add(alpha, *x, beta, *y, *y);
   return 0;
}

int LinearElasticityApp::SpatialNorm(
      braid_Vector  u_,
      double       *norm_ptr)
{
   BlockVector *u = (BlockVector *) u_;
   *norm_ptr = ode_oper->Norm(*u);
   return 0;
}

int LinearElasticityApp::BufSize(
      int               *size_ptr,
      BraidBufferStatus &status)
{
   *size_ptr = sizeof(double)*true_offsets.Last();
   return 0;
}

int LinearElasticityApp::BufPack(
      braid_Vector       u_,
      void              *buffer,
      BraidBufferStatus &status)
{
   BlockVector *u = (BlockVector *) u_;
   std::copy(u->GetData(), u->GetData() + u->Size(), (double *)buffer);
   status.SetSize(sizeof(double)*u->Size()); // optional?
   return 0;
}

int LinearElasticityApp::BufUnpack(
      void              *buffer,
      braid_Vector      *u_ptr,
      BraidBufferStatus &status)
{
   BlockVector *u = new BlockVector(true_offsets);
   std::copy((double *)buffer, (double *)buffer + u->Size(), u->GetData());
   *u_ptr = (braid_Vector) u;
   return 0;
}

// optional
int LinearElasticityApp::Coarsen(
      braid_Vector           fu_,
      braid_Vector          *cu_ptr,
      BraidCoarsenRefStatus &status)
{ return 0; }

// optional
int LinearElasticityApp::Refine(
      braid_Vector           cu_,
      braid_Vector          *fu_ptr,
      BraidCoarsenRefStatus &status)
{ return 0; }

int LinearElasticityApp::Access(
      braid_Vector       u_,
      BraidAccessStatus &astatus)
{
   BlockVector *u = (BlockVector *) u_;
   double t;
   int iter, level, done;
   astatus.GetTILD(&t, &iter, &level, &done);
   int cycle = (int) round((t-tstart)/opts.dt);
   v_gf.Distribute(u->GetBlock(0));
   x_gf.Distribute(u->GetBlock(1));
   double ee = ode_oper->ElasticEnergy(x_gf);
   double ke = ode_oper->KineticEnergy(v_gf);
   bool root = (fespace.GetMyRank() == 0);
   if (root)
   {
      cout << "Access: lev: " << setw(2) << level
           << ", iter: " << setw(2) << iter
           << ", cycle: " << setw(3) << cycle
           << ", t: " << setw(10) << t
           << ", EE: " << ee << ", KE: " << ke << ", ΔTE: "
           << (ee+ke)-(ee0+ke0) << endl;
   }
   if (opts.visualization && (cycle % opts.vis_steps == 0))
   {
      if (root)
      {
         cout << "Visualizing ..." << endl;
      }
      visualize(vis, &pmesh, &x_gf, &v_gf);
      vis << "window_title 'velocity, cycle " << cycle
          << ", iter " << iter << "'" << endl;
      if (opts.screenshots)
      {
         if (root)
         {
            cout << "Screenshot ..." << endl;
         }
         vis << "screenshot braid_" << setw(2) << setfill('0') << iter << "_"
             << setw(6) << setfill('0') << cycle << ".png" << endl;
      }
   }
   return 0;
}


void InitialDeformation(const Vector &x, Vector &y)
{
   // set the initial configuration to be the same as the reference, stress
   // free, configuration
   y = x;
}

void InitialVelocity(const Vector &x, Vector &v)
{
   const int dim = x.Size();
   const double s = 0.1/64.;

   v = 0.0;
   v(dim-1) = s*x(0)*x(0)*(8.0-x(0));
   v(0) = -s*x(0)*x(0);
}
