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

using namespace std;
using namespace mfem;

// Additional "HYPRE" functions
namespace hypre
{
   /**------------------------------------------------------------------------
    * Perform the operation A += beta*B, assuming that the sparsity pattern of
    * A contains that of B.
    * -----------------------------------------------------------------------*/
   HYPRE_Int
   hypre_CSRMatrixSum( hypre_CSRMatrix *A,
                       HYPRE_Complex    beta,
                       hypre_CSRMatrix *B )
   {
      HYPRE_Complex    *A_data   = hypre_CSRMatrixData(A);
      HYPRE_Int        *A_i      = hypre_CSRMatrixI(A);
      HYPRE_Int        *A_j      = hypre_CSRMatrixJ(A);
      HYPRE_Int         nrows_A  = hypre_CSRMatrixNumRows(A);
      HYPRE_Int         ncols_A  = hypre_CSRMatrixNumCols(A);
      HYPRE_Complex    *B_data   = hypre_CSRMatrixData(B);
      HYPRE_Int        *B_i      = hypre_CSRMatrixI(B);
      HYPRE_Int        *B_j      = hypre_CSRMatrixJ(B);
      HYPRE_Int         nrows_B  = hypre_CSRMatrixNumRows(B);
      HYPRE_Int         ncols_B  = hypre_CSRMatrixNumCols(B);

      HYPRE_Int         ia, j, pos;
      HYPRE_Int        *marker;

      if (nrows_A != nrows_B || ncols_A != ncols_B)
      {
         hypre_printf("Warning! incompatible matrix dimensions!\n");
         return -1;
      }

      marker = hypre_CTAlloc(HYPRE_Int, ncols_A);
      for (ia = 0; ia < ncols_A; ia++)
         marker[ia] = -1;

      for (ia = 0; ia < nrows_A; ia++)
      {
         for (j = A_i[ia]; j < A_i[ia+1]; j++)
            marker[A_j[j]] = j;

         for (j = B_i[ia]; j < B_i[ia+1]; j++)
         {
            pos = marker[B_j[j]];
            if (pos < A_i[ia])
            {
               hypre_printf("Warning! hypre_CSRMatrixSum: Incorrect input!\n");
               return -2;
            }
            A_data[pos] += beta * B_data[j];
         }
      }

      hypre_TFree(marker);
      return 0;
   }

   /**------------------------------------------------------------------------
    * Return a new matrix containing the sum of A and B, assuming that both
    * matrices use the same row and column partitions and the same col_map_offd
    * arrays.
    * -----------------------------------------------------------------------*/
   hypre_ParCSRMatrix *
   hypre_ParCSRMatrixAdd( hypre_ParCSRMatrix *A,
                          hypre_ParCSRMatrix *B )
   {
      MPI_Comm            comm   = hypre_ParCSRMatrixComm(A);
      hypre_CSRMatrix    *A_diag = hypre_ParCSRMatrixDiag(A);
      hypre_CSRMatrix    *A_offd = hypre_ParCSRMatrixOffd(A);
      HYPRE_Int          *A_cmap = hypre_ParCSRMatrixColMapOffd(A);
      HYPRE_Int           A_cmap_size = hypre_CSRMatrixNumCols(A_offd);
      hypre_CSRMatrix    *B_diag = hypre_ParCSRMatrixDiag(B);
      hypre_CSRMatrix    *B_offd = hypre_ParCSRMatrixOffd(B);
      HYPRE_Int          *B_cmap = hypre_ParCSRMatrixColMapOffd(B);
      HYPRE_Int           B_cmap_size = hypre_CSRMatrixNumCols(B_offd);
      hypre_ParCSRMatrix *C;
      hypre_CSRMatrix    *C_diag;
      hypre_CSRMatrix    *C_offd;
      HYPRE_Int          *C_cmap;
      HYPRE_Int           im;

      /* Make sure A_cmap and B_cmap are the same. */
      if (A_cmap_size != B_cmap_size)
      {
         hypre_printf("Warning! hypre_ParCSRMatrixAdd: cmap_size!\n");
         return NULL;
      }
      for (im = 0; im < A_cmap_size; im++)
      {
         if (A_cmap[im] != B_cmap[im])
         {
            hypre_printf("Warning! hypre_ParCSRMatrixAdd: cmap!\n");
            return NULL;
         }
      }

      /* Add diagonals, off-diagonals, copy cmap. */
      C_diag = hypre_CSRMatrixAdd(A_diag, B_diag);
      C_offd = hypre_CSRMatrixAdd(A_offd, B_offd);
      C_cmap = hypre_TAlloc(HYPRE_Int, A_cmap_size);
      for (im = 0; im < A_cmap_size; im++)
         C_cmap[im] = A_cmap[im];

      C = hypre_ParCSRMatrixCreate( comm,
                                    hypre_ParCSRMatrixGlobalNumRows(A),
                                    hypre_ParCSRMatrixGlobalNumCols(A),
                                    hypre_ParCSRMatrixRowStarts(A),
                                    hypre_ParCSRMatrixColStarts(A),
                                    hypre_CSRMatrixNumCols(C_offd),
                                    hypre_CSRMatrixNumNonzeros(C_diag),
                                    hypre_CSRMatrixNumNonzeros(C_offd) );

      /* In C, destroy diag/offd (allocated by Create) and replace them with
         C_diag/C_offd. */
      hypre_CSRMatrixDestroy(hypre_ParCSRMatrixDiag(C));
      hypre_CSRMatrixDestroy(hypre_ParCSRMatrixOffd(C));
      hypre_ParCSRMatrixDiag(C) = C_diag;
      hypre_ParCSRMatrixOffd(C) = C_offd;

      hypre_ParCSRMatrixColMapOffd(C) = C_cmap;

      /* hypre_ParCSRMatrixSetNumNonzeros(A); */

      /* Make sure that the first entry in each row is the diagonal one. */
      hypre_CSRMatrixReorder(hypre_ParCSRMatrixDiag(C));

      /* C owns diag, offd, and cmap. */
      hypre_ParCSRMatrixSetDataOwner(C, 1);
      /* C does not own row and column starts. */
      hypre_ParCSRMatrixSetRowStartsOwner(C, 0);
      hypre_ParCSRMatrixSetColStartsOwner(C, 0);

      return C;
   }

   /**------------------------------------------------------------------------
    * Perform the operation A += beta*B, assuming that both matrices use the
    * same row and column partitions and the same col_map_offd arrays. Also,
    * it is assumed that the sparsity pattern of A contains that of B.
    * -----------------------------------------------------------------------*/
   HYPRE_Int
   hypre_ParCSRMatrixSum( hypre_ParCSRMatrix *A,
                          HYPRE_Complex       beta,
                          hypre_ParCSRMatrix *B )
   {
      hypre_CSRMatrix *A_diag = hypre_ParCSRMatrixDiag(A);
      hypre_CSRMatrix *A_offd = hypre_ParCSRMatrixOffd(A);
      hypre_CSRMatrix *B_diag = hypre_ParCSRMatrixDiag(B);
      hypre_CSRMatrix *B_offd = hypre_ParCSRMatrixOffd(B);
      HYPRE_Int error;

      error = hypre_CSRMatrixSum(A_diag, beta, B_diag);
      if (!error)
         error = hypre_CSRMatrixSum(A_offd, beta, B_offd);

      return error;
   }

   HYPRE_Int
   hypre_CSRMatrixSetConstantValues( hypre_CSRMatrix *A,
                                     HYPRE_Complex    value )
   {
      HYPRE_Complex *A_data = hypre_CSRMatrixData(A);
      HYPRE_Int      A_nnz  = hypre_CSRMatrixNumNonzeros(A);
      HYPRE_Int      ia;

      for (ia = 0; ia < A_nnz; ia++)
         A_data[ia] = value;

      return 0;
   }

   HYPRE_Int
   hypre_ParCSRMatrixSetConstantValues( hypre_ParCSRMatrix *A,
                                        HYPRE_Complex       value )
   {
      hypre_CSRMatrixSetConstantValues(hypre_ParCSRMatrixDiag(A), value);
      hypre_CSRMatrixSetConstantValues(hypre_ParCSRMatrixOffd(A), value);

      return 0;
   }
}

using namespace hypre;


class FuncDepGridFuncCoefficient : public Coefficient
{
protected:
   ParGridFunction *gfunc;
   double power;
   Vector grad;
public:

   FuncDepGridFuncCoefficient(ParGridFunction *vd, double _power)
   {
      gfunc = vd;
      power = _power;
   }

   void SetGridFunction(HypreParVector *X)
   {
      *gfunc = *X; // distribute (communication)
   }

   void SetPower(double _power) { power = _power; }

   double Eval(ElementTransformation &T, const IntegrationPoint &ip)
   {
      gfunc->GetGradient(T, grad);
      return -pow(grad * grad, power/2);
   }

   virtual void Read(istream &in) { }
};


// All variables required for a nonlinear solve with variable tolerance.
class NonlinearOptions
{
   // timeint = tstart-tstop, tolf and tolc are fine and coarse grid nonlinear
   // solve tolerances and ntime is the number of time steps
protected:
   int    max_it, ntime;
   double timeint, tolf, tolc;

public:
   NonlinearOptions(int _maxit=100, int _ntime=32, double _timeint=1,
                    double _tolf=0.001, double _tolc=0.001)
   {
      tolf = _tolf;
      tolc = _tolc;
      max_it = _maxit;
      ntime = _ntime;
      timeint = _timeint;
   }

   double GetMaxIteration() {return max_it;}

   double GetTolerance(double dt = 0)
   {
      if (dt - timeint/ntime < 0.000001)
         return tolf;
      else
         return tolc;
   }
};


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
   HypreParMatrix *M, *B; // B = M - dt A, for ImplicitSolve
   mutable HypreParMatrix *A;
   HypreParVector *X, *Y, *Z;
   HypreParVector *W, *H, *D;
   mutable HypreParVector *b;
   ParLinearForm *b_form;
   mutable ParBilinearForm *aform;
   Coefficient *b_coeff;
   VectorCoefficient *b_vcoeff;
   MatrixCoefficient *mcoeff;
   FuncDepGridFuncCoefficient *gcoeff;
   HypreBoomerAMG *amg;
   HyprePCG *pcg;
   double current_dt; // dt used in B
   HypreBoomerAMG *B_amg;
   HyprePCG *B_pcg;
   NonlinearOptions *NLO;
   int timedep_A; // Is 'a' time-dependent or nonlinear?

public:
   LinearSystemODE(ParBilinearForm *_aform,
                   int _timedep_A = 1,
                   HypreParMatrix *_M = NULL,
                   ParLinearForm *_b_form = NULL,
                   Coefficient *_b_coeff = NULL,
                   VectorCoefficient *_b_vcoeff = NULL,
                   MatrixCoefficient *_mcoeff = NULL,
                   FuncDepGridFuncCoefficient *_gcoeff = NULL,
                   NonlinearOptions *_NLO = NULL)
      : M(_M), b_form(_b_form), aform(_aform), b_coeff(_b_coeff),
        b_vcoeff(_b_vcoeff), mcoeff(_mcoeff), gcoeff(_gcoeff), NLO(_NLO),
        timedep_A(_timedep_A)
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
         Z = W = H = D = NULL;
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
         W = new HypreParVector(*A);
         H = new HypreParVector(*A);
         D = new HypreParVector(*A);
      }

      b = NULL;
      B = NULL; // Allocated in ImplicitSolve if needed
      current_dt = -1.;
      B_amg = NULL;
      B_pcg = NULL;
   }

   void AssembleDiffMatrix() const
   {
      if (timedep_A)
      {
         delete A;

         if (mcoeff)
            mcoeff->SetTime(GetTime());

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
      if (gcoeff)
         gcoeff->SetGridFunction(X);
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
         A->Mult(*X, *Z); // M^{-1} A X
         if (b)
            *Z += *b; // M^{-1} (A X + b)

         pcg->Mult(*Z, *Y);
      }
   }

   /** Solve for k in the equation: M k = A(x + dt k, t) + b(t).
       In the nonlinear case we solve using a fixed point iteration for the
       equation written as:
       *     k = [M - dt*AA(x + dt*k, t)]^{-1} [AA(x + dt*k, t)*x + b(t)],
       where AA(z,t) is the matrix assembled using the grid function z in the
       coefficient of the bilinear form, so that the nonlinear operator A(z,t)
       can be written as A(z,t) = AA(z,t)*z.
       This method is needed for backward Euler and other DIRK methods. */
   virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k)
   {
      double resid = 100;
      int i = 0;
      int iteration = 0;
      double tol = 1;

      if (NLO)
      {
         iteration = NLO->GetMaxIteration();
         tol = NLO->GetTolerance(dt);
      }

      // Wrap the Data
      X->SetData(x.GetData());
      Y->SetData(k.GetData());

      // Assemble
      AssembleBVector();
      if (gcoeff)
         gcoeff->SetGridFunction(X);
      AssembleDiffMatrix();

      while (resid >= tol)
      {
         if (!M)
         {
            cerr << "Not implemented!" << endl;
            abort();
         }

         // 1) Form the matrix: B = M - dt A
         if (!B)
            B = new HypreParMatrix(hypre_ParCSRMatrixAdd(*M, *A));
         if (timedep_A || current_dt != dt)
         {
            hypre_ParCSRMatrixSetConstantValues(*B, 0.0);
            hypre_ParCSRMatrixSum(*B, 1.0, *M);
            hypre_ParCSRMatrixSum(*B, -dt, *A);

            delete B_amg;
            B_amg = new HypreBoomerAMG(*B);
            B_amg->SetPrintLevel(-1);

            delete B_pcg;
            B_pcg = new HyprePCG(*B);
            B_pcg->SetTol(1e-12);
            B_pcg->SetMaxIter(200);
            B_pcg->SetPrintLevel(0);
            B_pcg->SetPreconditioner(*B_amg);
            B_pcg->SetZeroInintialIterate();

            current_dt = dt;
         }

         // 2) Form the rhs: A x + b
         A->Mult(*X, *Z);
         if (b)
            *Z += *b;

         // 3) Solve the system B Y = Z
         B_pcg->Mult(*Z, *Y);

         if (++i >= iteration)
            break;

         // Update the matrix for the next iteration and calculate the residual;
         add(*X, dt, *Y, *D); // D = X + k*dt
         if (gcoeff)
         {
            gcoeff->SetGridFunction(D);
            AssembleDiffMatrix();
         }

         A->Mult(*D, *H);
         *H += *b;   // H = A(x + k dt) + b(t)
         M->Mult(*Y, *W);
         *W -= *H;   // W = M k - A(x + k*dt) - b(t)
         resid = sqrt(InnerProduct(W, W));

         cout << "LinearSystemODE::ImplicitSolve : resid = "
              << resid << endl;
      }

      cout << "LinearSystemODE::ImplicitSolve : number of iterations = "
           << i << endl;
   }

   virtual ~LinearSystemODE()
   {
      delete B_pcg;
      delete B_amg;
      delete B;
      delete b;
      delete A;
      delete Z;
      delete pcg;
      delete amg;
      delete Y;
      delete X;
      delete W;
      delete H;
      delete D;
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

//Flagterms
double diff_term;
double tfinal;
double power;


// Exact solution
double ExSol(Vector &p, double t)
{
   int dim = p.Size();

   if (dim == 2)
      return sin(kap*p(0))*sin(kap*p(1))*sin(tau*t);
   else
      return sin(kap*p(0))*sin(kap*p(1))*sin(kap*p(2))*sin(tau*t);
}

double DiffFunction(double t)
{
   if (diff_term < 10)
      return diff_term;
   else if (diff_term == 10)
      return cos(M_PI/2*t/(tfinal+ 0.01));
   else if (diff_term == 11)
      return t/tfinal + 0.001;
   else if (diff_term == 12)
   {
      if (t <= tfinal/4)
         return 0.001;
      else if (t > tfinal/4 && t <= tfinal/2)
         return 1;
      else if (t > tfinal/2 && t <= 3*tfinal/4)
         return 0.001;
      else
         return 1;
   }
   else
      return 1;
}

// matrix diffusion coefficient term
void Diff(const Vector &p, double t, DenseMatrix &m)
{
   int dim = p.Size();
   if (dim == 2)
   {
      double aa[4]={-1,0,0,-DiffFunction(t)};
      m = aa;
   }
   else
   {
      double aa[9] = {-1,0,0,0,-1,0,0,0,-DiffFunction(t)};
      m = aa;
   }
}

int timedep()
{
   if (diff_term < 10)
      return 0;
   else
      return 1;
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
   if (diff_term > 20)
   {
      double Ux, Uy;
      Ux = kap*cos(kap*p(0))*sin(kap*p(1))*sin(tau*t);
      Uy = kap*sin(kap*p(0))*cos(kap*p(1))*sin(tau*t);
      v(0) = pow(Ux*Ux + Uy*Uy,power/2)*Ux;
      v(1) = pow(Ux*Ux + Uy*Uy,power/2)*Uy;
      return;
   }
   else if (diff_term <= 0)
   {
      v(0) = -diff_term*kap*cos(kap*p(0))*sin(kap*p(1))*sin(tau*t);
      v(1) = -diff_term*kap*sin(kap*p(0))*cos(kap*p(1))*sin(tau*t);
      return;
   }
   else if (dim == 2)
   {
      v(0) = kap*cos(kap*p(0))*sin(kap*p(1))*sin(tau*t);
      v(1) = DiffFunction(t)*kap*sin(kap*p(0))*cos(kap*p(1))*sin(tau*t);
      return;
   }
   // 3D if only implimented for the matrix case and not the nonlinear case.
   v(0) = kap*cos(kap*p(0))*sin(kap*p(1))*sin(kap*p(2))*sin(tau*t);
   v(1) = kap*sin(kap*p(0))*cos(kap*p(1))*sin(kap*p(2))*sin(tau*t);
   v(2) = DiffFunction(t)*kap*sin(kap*p(0))*sin(kap*p(1))*cos(kap*p(2))*sin(tau*t);
}

// Source term
double ST(Vector &p, double t)
{
   int dim = p.Size();
   if (diff_term > 20)
   {
      double Ut,Ux,Uxx,Uy,Uxy;
      Ut = tau*sin(kap*p(0))*sin(kap*p(1))*cos(tau*t);
      Ux = kap*cos(kap*p(0))*sin(kap*p(1))*sin(tau*t);
      Uy = kap*sin(kap*p(0))*cos(kap*p(1))*sin(tau*t);
      Uxx = -kap*kap*sin(kap*p(0))*sin(kap*p(1))*sin(tau*t);
      Uxy =  kap*kap*cos(kap*p(0))*cos(kap*p(1))*sin(tau*t);
      if (power == 0)
         return  Ut-(power+2)*pow((Ux*Ux +Uy*Uy),power/2)*Uxx; //
      else
         return (Ut-(power+2)*pow((Ux*Ux +Uy*Uy),power/2)*Uxx -
                 2*power*(Ux*Uy*Uxy)*pow(Ux*Ux+Uy*Uy,power/2-1));
   }
   else if (diff_term <= 0)
      return ((tau*cos(tau*t) -(2*diff_term)*kap*kap*sin(tau*t))*
              sin(kap*p(0))*sin(kap*p(1)));
   else if (dim == 2)
      return ((tau*cos(tau*t) +(1+DiffFunction(t))*kap*kap*sin(tau*t))*
              sin(kap*p(0))*sin(kap*p(1)));
   else
      return ((tau*cos(tau*t) + (2+DiffFunction(t))*kap*kap*sin(tau*t))*
              sin(kap*p(0))*sin(kap*p(1))*sin(kap*p(2)));
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
   int    picard        = 5;
   double tolf          = 0.0001;
   double tolc          = tolf;


   extern double diff_term;
   extern double tfinal;
   extern double power;
   extern double kap;

   kap = M_PI;
   power = 0;
   diff_term = -1;

   // double cfl         = 1.0;
   int    ode_solver  = -1;

   int heat_equation = 1;
   int scalar_ode_option = 2;

   // Vis parameters
   const char * vishost = "localhost";
   int visport = 19916;


   // BRAID default parameters:
   int    max_levels  = 10;
   int    min_coarse  = 3;
   int    nrelax      = 1;
   int    nrelax0     = -1;
   double tol         = 1e-9;
   int    tnorm       = 2;
   int    cfactor     = 2;
   int    cfactor0    = -1;
   int    max_iter    = 100;
   int    fmg         = 0;
   int    nfmg_Vcyc   = 1;
   int    access_level= 1;
   bool   wrapper_tests = false;
   bool   one_wrapper_test = false;

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
      else if (strcmp(argv[arg_index], "-picard") == 0)
      {
         picard = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-tolf") == 0)
      {
         tolf = atof(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-tolc") == 0)
      {
         tolc = atof(argv[++arg_index]);
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
      else if (strcmp(argv[arg_index], "-dt") == 0)
      {
         diff_term = atof(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-pow") == 0)
      {
         power = atof(argv[++arg_index]);
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
      else if( strcmp(argv[arg_index], "-tnorm") == 0 ){
          arg_index++;
          tnorm = atoi(argv[arg_index++]);
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
         "  -dt <diff_term>           :   -dt<0 - constant coefficient with value = dt \n "
         "                    : set D(t) = [1 0 ; 0 f(t)], the diffusion tensor \n"
         "                         0<dt<10 f(t) = dt\n"
         "                         10 - f(t) = cos(pi*t/(2*tstop))+0.001\n"
         "                         11 - f(t) = t/tstop + 0.001 \n"
         "                         12 - f(t) = jump function \n"
         "                         21 - nonlinear coefficient a = |grad(u)|^p\n"
         "  -pow <power>      : set the value of the power in nonlinear equation (dt = 21)\n"
         "  -kap <kap>        : set the frequency of the exact solution sin wave\n"
         "  -picard <picard>      : set maximum number of picard iterations for each nonlinear solve\n"
         "  -tolf <tolf>      : set the nonlinear solve tolerance for the fine grid (default 0.001)\n"
         "  -tolc <tolc>      : set the nonlinear solve tolerance for the coarse grids (default 0.001)\n"
         "  -ml  <max_levels> : set max number of time levels (default: 10)\n"
         "  -mc  <min_coarse> : set minimum possible coarse level size (default: 3)\n"
         "  -nu  <nrelax>     : set num F-C relaxations (default: 1)\n"
         "  -nu0 <nrelax>     : set num F-C relaxations on level 0\n"
         "  -tol <tol>        : set stopping tolerance (default: 1e-9)\n"
         "  -tnorm <tnorm>    : set temporal norm \n"
         "                      1 - One-norm \n"
         "                      2 - Two-norm (default) \n"
         "                      3 - Infinity-norm \n"
         "  -cf  <cfactor>    : set coarsening factor (default: 2)\n"
         "  -cf0 <cfactor0>   : set aggressive coarsening (default: off)\n"
         "  -mi  <max_iter>   : set max iterations (default: 100)\n"
         "  -fmg <nfmg_Vcyc>  : use FMG cycling with nfmg_Vcyc V-cycles at each fmg level\n"
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
   if (pmesh->GetNodes())
      fec = pmesh->GetNodes()->OwnFEC();
   else
      fec = new LinearFECollection;
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
   ConstantCoefficient minus(diff_term);
   FunctionCoefficient source(&ST);
   VectorFunctionCoefficient nbc(dim, &NBC);
   FunctionCoefficient exact_sol(&ExSol);
   MatrixFunctionCoefficient dc(dim, &Diff);
   FuncDepGridFuncCoefficient gfc(&x0, power);
   ParBilinearForm *a = new ParBilinearForm(fespace);
   ParBilinearForm *m = new ParBilinearForm(fespace);
   ParLinearForm   *b = new ParLinearForm(fespace);
   NonlinearOptions *NLO = new NonlinearOptions(picard, ntime, tstop - tstart,
                                                tolf, tolc);

   tfinal = tstop;

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
   if (heat_equation)
   {
      if (diff_term <= 0)
      {
         // constant diffusion coefficient = (-diff_term)
         a->AddDomainIntegrator(new DiffusionIntegrator(minus));
         ode = new LinearSystemODE(a, 0, M, b, &source, &nbc);
      }
      else if (diff_term < 20)
      {
         // matrix diffusion coefficient (possibly time dependent)
         a->AddDomainIntegrator(new DiffusionIntegrator(dc));
         ode = new LinearSystemODE(a, timedep(), M, b, &source, &nbc, &dc);
      }
      else
      {
         // nonlinear diffusion coefficient = |grad(u)|^p
         a->AddDomainIntegrator(new DiffusionIntegrator(gfc));
         ode = new LinearSystemODE(a, timedep(), M, b, &source, &nbc, NULL,
                                   &gfc, NLO);
      }
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
      MFEMBraidApp app(comm_t, ode, X0, &x0, solver, tstart, tstop, ntime);
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
         core.SetAbsTol(tol);
         core.SetCFactor(-1, cfactor);
         core.SetAggCFactor(cfactor0);
         core.SetMaxIter(max_iter);
         core.SetTemporalNorm(tnorm);
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
   if (!pmesh->GetNodes())
      delete fec;
   delete pmesh;

   MPI_Comm_free(&comm_x);
   MPI_Comm_free(&comm_t);

   MPI_Finalize();

   return 0;
}
