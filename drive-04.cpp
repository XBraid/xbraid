// Example 04
//
// Compile with: make drive-04
//
// Sample run:   mpirun -np 2 drive-04 ../mfem/data/star.mesh
//
// Description:  Unstructured 2D/3D parabolic problem in MFEM.

#include <fstream>
#include "mfem.hpp"
#include "warp.h"


// Right-hand side for the set of scalar ODEs du_i/dt = f(u_i,t)
class ScalarODE : public TimeDependentOperator
{
private:
   inline double f(double u, double t) const
   {
      // return 1.0;                         // u(t) = t+u(0)
      // return 4*pow(t,3)-3*pow(t,2)+2*t-1; // u(t) = t^4-t^3+t^2-t+u(0)
      // return pow(u,2)+1;                  // u(t) = tan(t) for u(0) = 0
      return 2*t*u+t;                        // u(t) = (1/2+u(0))*e^(t^2)-1/2
   }

public:
   ScalarODE(int _size) { size = _size; }

   void Mult(const Vector &x, Vector &y) const
   {
      double t = GetTime();
      for (int i = 0; i < x.Size(); i++)
         y(i) = f(x(i),t);
   }
};


// Right-hand side for the linear system of ODEs
//               M du/dt = A u + b,
// where A, M and b are time-independent
class LinearSystemODE : public TimeDependentOperator
{
private:
   HypreParMatrix *A;
   HypreParVector *b;
   HypreParMatrix *M;
   HypreBoomerAMG *amg;
   HyprePCG *pcg;
   HypreParVector *X, *Y, *Z;

public:
   LinearSystemODE(HypreParMatrix *_A,
                   HypreParMatrix *_M = NULL,
                   HypreParVector *_b = NULL)
      : A(_A), b(_b), M(_M)
   {
      size = A->Size();
      X = new HypreParVector(*A);
      Y = new HypreParVector(*A);

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

         Z = new HypreParVector(*A);
      }
   }

   void Mult(const Vector &x, Vector &y) const
   {
      X -> SetData(x.GetData());
      Y -> SetData(y.GetData());

      if (!M)
      {
         if (!b)  // A X
            A->Mult(*X, *Y);
         else     // A X + b
         {
            *Y = *b;
            A->Mult(*X, *Y, 1.0, 1.0);
         }
      }
      else
      {
         if (!b)  // M^{-1} A X
            A->Mult(*X, *Z);
         else     // M^{-1} (A X + b)
         {
            *Z = *b;
            A->Mult(*X, *Z, 1.0, 1.0);
         }
         *Y = 0.0;
         pcg->Mult(*Z, *Y);
      }
   }

   ~LinearSystemODE()
   {
      delete X;
      delete Y;
      if (amg)
         delete amg;
      if (pcg)
         delete pcg;
      if (Z)
         delete Z;
   }
};


// Solve the given ODE with the given initial condtion using uniform
// time-stepping and (optionally) visualize the solution.
void SolveODE(TimeDependentOperator *ode, HypreParVector *X0,
              int ode_solver = 1, double tstart = 0.0, double tstop = 1.0,
              int ntime = 100, ParGridFunction *x = NULL)
{
   ODESolver *solver;

   if (ode_solver == 4)
      solver = new RK4Solver;
   else if (ode_solver == 3)
      solver = new RK3SSPSolver;
   else if (ode_solver == 2)
      solver = new RK2Solver(1.0);
   else if (ode_solver == 1)
      solver = new ForwardEulerSolver;

   solver->Init(*ode);

   double t = tstart;
   double dt = (tstop-tstart)/ntime;

   socketstream *sol_sock = NULL;
   ParMesh *pmesh = NULL;
   if (x)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      sol_sock = new socketstream(vishost, visport);
      pmesh = x->ParFESpace()->GetParMesh();

      *x = *X0;
      *sol_sock << "parallel " << pmesh->GetNRanks()
                << " " << pmesh->GetMyRank() << "\n";
      *sol_sock << "solution\n";
      sol_sock->precision(8);
      pmesh->Print(*sol_sock);
      x->Save(*sol_sock);
      *sol_sock << "autoscale value\n";
      *sol_sock << "keys cmA\n";
      *sol_sock << flush;
   }

   for (int i = 0; i < ntime; i++)
   {
      solver->Step(*X0, t, dt);
      if (x)
      {
         *x = *X0;
         *sol_sock << "parallel " << pmesh->GetNRanks()
                   << " " << pmesh->GetMyRank() << "\n";
         *sol_sock << "solution\n";
         pmesh->Print(*sol_sock);
         x->Save(*sol_sock);
         *sol_sock << flush;
      }
   }

   if (x && !pmesh->GetMyRank())
      cout << endl;

   delete sol_sock;
}


// Wrapper for WARP's App object
class WarpApp
{
public:
   TimeDependentOperator *ode;
   ODESolver *solver;

   MPI_Comm comm_t;
   double   tstart;
   double   tstop;
   int      ntime;

   HypreParVector *X0;
   ParGridFunction *x;
   int buff_size;
   socketstream *sol_sock;
   ParMesh *pmesh;

   WarpApp(MPI_Comm _comm_t, TimeDependentOperator *_ode,
           HypreParVector *_X0, ParGridFunction *_x, int ode_solver = 1,
           double _tstart = 0.0, double _tstop = 1.0, int _ntime = 100)
      : ode(_ode), comm_t(_comm_t), tstart(_tstart), tstop(_tstop), ntime(_ntime),
        X0(_X0), x(_x)
   {
      buff_size = sizeof(double) * X0->Size();

      if (ode_solver == 4)
         solver = new RK4Solver;
      else if (ode_solver == 3)
         solver = new RK3SSPSolver;
      else if (ode_solver == 2)
         solver = new RK2Solver(1.0);
      else if (ode_solver == 1)
         solver = new ForwardEulerSolver;

      solver->Init(*ode);

      char vishost[] = "localhost";
      int  visport   = 19916;
      sol_sock = new socketstream(vishost, visport);
      pmesh = x->ParFESpace()->GetParMesh();

      *x = *X0;
      *sol_sock << "parallel " << pmesh->GetNRanks()
                << " " << pmesh->GetMyRank() << "\n";
      *sol_sock << "solution\n";
      sol_sock->precision(8);
      pmesh->Print(*sol_sock);
      x->Save(*sol_sock);
      *sol_sock << "autoscale value\n";
      *sol_sock << "keys cmA\n";
      *sol_sock << flush;
   }

   // Below warp_Vector == HypreParVector* and  warp_App == WarpApp*

   static int Phi(warp_App     _app,
                  double       tstart,
                  double       tstop,
                  double       accuracy,
                  int          gzero,
                  warp_Vector  _u,
                  int         *rfactor_ptr)
   {
      WarpApp *app = (WarpApp*) _app;
      HypreParVector *u = (HypreParVector*) _u;

      double t = tstart;
      double dt = tstop-tstart;

      app->solver->Step(*u, t, dt);

      if (!gzero)
      {
         // this could be handled in the LinearSystemODE class
      }

      // no refinement
      *rfactor_ptr = 1;

      return 0;
   }

   static int Clone(warp_App     _app,
                    warp_Vector  _u,
                    warp_Vector *v_ptr)
   {
      HypreParVector *u = (HypreParVector*) _u;
      HypreParVector *v = new HypreParVector((const HypreParVector &)*u);
      *v = *u;
      *v_ptr = (warp_Vector) v;
      return 0;
   }

   static int Init(warp_App    _app,
                   double       t,
                   warp_Vector *u_ptr)
   {
      WarpApp *app = (WarpApp*) _app;
      Clone(_app, (warp_Vector)app->X0, u_ptr);
      if (t != app->tstart)
      {
         HypreParVector *u = (HypreParVector*) *u_ptr;
         u->Randomize(2142);
      }
      return 0;
   }

   static int Free(warp_App    _app,
                   warp_Vector _u)
   {
      HypreParVector *u = (HypreParVector*) _u;
      delete u;
      return 0;
   }

   static int Sum(warp_App    _app,
                  double      alpha,
                  warp_Vector _x,
                  double      beta,
                  warp_Vector _y)
   {
      HypreParVector *x = (HypreParVector*) _x;
      HypreParVector *y = (HypreParVector*) _y;
      add(alpha, *x, beta, *y, *y);
      return 0;
   }

   static int Dot(warp_App     _app,
                  warp_Vector  _u,
                  warp_Vector  _v,
                  double      *dot_ptr)
   {
      HypreParVector *u = (HypreParVector*) _u;
      HypreParVector *v = (HypreParVector*) _v;
      *dot_ptr = InnerProduct(u,v);
      return 0;
   }

   static int Write(warp_App     _app,
                    double       t,
                    warp_Vector  _u)
   {
      WarpApp *app = (WarpApp*) _app;
      HypreParVector *u = (HypreParVector*) _u;
      (*app->x) = *u;

      if (t != app->tstart)
      {
         *app->sol_sock << "parallel " << app->pmesh->GetNRanks()
                        << " " << app->pmesh->GetMyRank() << "\n";
         *app->sol_sock << "solution\n";
         app->pmesh->Print(*app->sol_sock);
         app->x->Save(*app->sol_sock);
         *app->sol_sock << flush;
      }

      return 0;
   }

   static int BufSize(warp_App  _app,
                      int      *size_ptr)
   {
      WarpApp *app = (WarpApp*) _app;
      *size_ptr = app->buff_size;
      return 0;
   }

   static int BufPack(warp_App     _app,
                      warp_Vector  _u,
                      void        *buffer)
   {
      WarpApp *app = (WarpApp*) _app;
      HypreParVector *u = (HypreParVector*) _u;
      memcpy(buffer, u->GetData(), app->buff_size);
      return 0;
   }

   static int BufUnpack(warp_App     _app,
                        void        *buffer,
                        warp_Vector *u_ptr)
   {
      WarpApp *app = (WarpApp*) _app;
      Clone(_app, (warp_Vector)app->X0, u_ptr);
      HypreParVector *u = (HypreParVector*) *u_ptr;
      memcpy(u->GetData(), buffer, app->buff_size);
      return 0;
   }
};


// Wrapper for WARP's core object
class WarpCore
{
private:
   warp_Core core;

public:
   WarpCore(MPI_Comm comm_world, WarpApp *app)
   {
      warp_Init(comm_world,
                app->comm_t, app->tstart, app->tstop, app->ntime, (warp_App)app,
                WarpApp::Phi, WarpApp::Init, WarpApp::Clone, WarpApp::Free,
                WarpApp::Sum, WarpApp::Dot, WarpApp::Write,
                WarpApp::BufSize, WarpApp::BufPack, WarpApp::BufUnpack, &core);
   }

   void SetMaxLevels(int max_levels) { warp_SetMaxLevels(core, max_levels); }
   void SetNRelax(int level, int nrelax) { warp_SetNRelax(core, level, nrelax); }
   void SetAbsTol(double tol) { warp_SetAbsTol(core, tol); }
   void SetCFactor(int level, int cfactor) { warp_SetCFactor(core, level, cfactor); }
   void SetMaxIter(int max_iter) { warp_SetMaxIter(core, max_iter); }
   void SetFMG() { warp_SetFMG(core); }

   void SetDefaults()
   {
      int    max_levels = 1;
      int    nrelax     = 1;
      int    nrelax0    = -1;
      double tol        = 1.0e-06;
      int    cfactor    = 2;
      int    max_iter   = 100;
      int    fmg        = 0;

      SetMaxLevels(max_levels);
      SetNRelax(-1, nrelax);
      if (nrelax0 > -1)
      {
         SetNRelax(0, nrelax0);
      }
      SetAbsTol(tol);
      SetCFactor(-1, cfactor);
      SetMaxIter(max_iter);
      if (fmg)
         SetFMG();
   }

   void Drive() { warp_Drive(core); }

   void PrintStats() { warp_PrintStats(core); }

   ~WarpCore() { warp_Destroy(core); }
};


// Boundary conditions
double BC(Vector &x)
{
   return 0;
}


int main (int argc, char *argv[])
{
   // Initialize MPI and split the space and time communicators
   MPI_Comm comm, comm_x, comm_t;
   int num_procs, myid, myid_x;
   comm = MPI_COMM_WORLD;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(comm, &num_procs);
   MPI_Comm_rank(comm, &myid);

   int num_procs_x = 1;
   int xcolor = myid / num_procs_x;
   int tcolor = myid % num_procs_x;
   MPI_Comm_split(comm, xcolor, myid, &comm_x);
   MPI_Comm_split(comm, tcolor, myid, &comm_t);

   MPI_Comm_rank(comm_x, &myid_x);

   Mesh *mesh;

   if (argc == 1)
   {
      if (myid == 0)
         cout << "\nUsage: mpirun -np <np> ex1p <mesh_file>\n" << endl;
      MPI_Finalize();
      return 1;
   }

   // Read the (serial) mesh from the given mesh file on all processors.
   ifstream imesh(argv[1]);
   if (!imesh)
   {
      if (myid == 0)
         cerr << "\nCan not open mesh file: " << argv[1] << '\n' << endl;
      MPI_Finalize();
      return 2;
   }
   mesh = new Mesh(imesh, 1, 1);
   imesh.close();

   // Refine the serial mesh on all processors to increase the resolution.
   {
      int ref_levels = 1;
      for (int l = 0; l < ref_levels; l++)
         mesh->UniformRefinement();
   }

   // Define a parallel mesh by a partitioning of the serial mesh. Refine this
   // mesh further in parallel to increase the resolution.
   ParMesh *pmesh = new ParMesh(comm_x, *mesh);
   delete mesh;
   {
      int par_ref_levels = 0;
      for (int l = 0; l < par_ref_levels; l++)
         pmesh->UniformRefinement();
   }

   // Define a parallel finite element space on the parallel mesh. We use the
   // finite elements coming from the mesh nodes (linear by default).
   FiniteElementCollection *fec;
   if (pmesh->GetNodes())
      fec = pmesh->GetNodes()->OwnFEC();
   else
      fec = new LinearFECollection;
   ParFiniteElementSpace *fespace = new ParFiniteElementSpace(pmesh, fec);
   int size = fespace->GlobalTrueVSize();
   if (myid == 0)
      cout << "Number of unknowns: " << size << endl;

   // Define the initial guess vector x0 as a parallel finite element grid
   // function corresponding to fespace.
   ParGridFunction x0(fespace);
   // x0.Randomize(1323);
   x0 = 1.0;
   Array<int> ess_bdr(pmesh->bdr_attributes.Max());
   ess_bdr = 1;
   FunctionCoefficient bc(&BC);
   x0.ProjectBdrCoefficient(bc, ess_bdr);
   HypreParVector *X0 = x0.ParallelAverage();

   // Set up the stiffness and mass parallel bilinear forms a(.,.) and m(.,.),
   // assemble and extract the corresponding parallel matrices.
   ConstantCoefficient minus_one(-1.0);
   ConstantCoefficient one(1.0);

   ParBilinearForm *a = new ParBilinearForm(fespace);
   ParBilinearForm *m = new ParBilinearForm(fespace);

   a->AddDomainIntegrator(new DiffusionIntegrator(minus_one));
   m->AddDomainIntegrator(new MassIntegrator(one));

   a->Assemble(); a->Finalize();
   m->Assemble(); m->Finalize();

   HypreParMatrix *A = a->ParallelAssemble();
   HypreParMatrix *M = m->ParallelAssemble();

   delete a;
   delete m;

   // Define the righ-hand side of the ODE and call MFEM or WARP to solve it.
   TimeDependentOperator *ode;
   // ode = new ScalarODE(A->Size());
   ode = new LinearSystemODE(A, M);

   int    ode_solver = 4; // RK4
   double tstart = 0.0;
   double tstop  = 1.0;
   int    ntime  = 200;

#if 1
   SolveODE(ode, X0, ode_solver, tstart, tstop, ntime, &x0);
#else
   WarpApp app(comm_t, ode, X0, &x0, ode_solver, tstart, tstop, ntime);
   WarpCore core(comm, &app);

   core.SetDefaults();
   core.Drive();
   core.PrintStats();
#endif

   // Free memory.
   delete ode;

   delete X0;
   delete M;
   delete A;

   delete fespace;
   if (!pmesh->GetNodes())
      delete fec;
   delete pmesh;

   MPI_Comm_free(&comm_x);
   MPI_Comm_free(&comm_t);

   MPI_Finalize();

   return 0;
}
