// Example 04
//
// Compile with: make drive-04
//
// Sample runs:  mpirun -np 4 drive-04
//               mpirun -np 4 drive-04 -pref 1 -px 2 -nt 100 -odesolver -2
//               mpirun -np 4 drive-04 -mesh ../../mfem/data/beam-quad.mesh -sref 2 -nowarp
//               mpirun -np 16 drive-04 -mesh ../../mfem/data/escher-p2.mesh -nowarp
//
// Description:  Unstructured 2D/3D ODE problems in MFEM.

#include <fstream>
#include "mfem.hpp"
#include "_warp.h"
#include "warp.h"
#include "warp_test.h"

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


// Features that will become part of MFEM
namespace mfem
{
   /** A time-dependent operator, f(x,t), that defines a method to solve the
       following equation: k = f(x + dt*k, t), for the unknown k.

       This method allows for the abstract implementation of some time
       integration methods, including diagonal implicit Runge-Kutta (DIRK)
       methods and backward Euler method in particular.

       Another (simpler) option is to make the solve method part of the
       TimeDependentOperator class. */
   class Diagonal_Implicit_Evolution_Operator : public TimeDependentOperator
   {
   public:
      /// Solve the equation: k = f(x + dt*k, t), for the unknown k.
      virtual void ImplicitSolve(const double dt, const Vector &x,
                                 Vector &k) = 0;
   };


   /// Backward Euler ODE solver. L-stable.
   class BackwardEulerSolver : public ODESolver
   {
   protected:
      Vector k;

   public:
      virtual void Init(TimeDependentOperator &_f)
      {
         ODESolver::Init(_f);
         k.SetSize(f->Size());
      }

      virtual void Step(Vector &x, double &t, double &dt)
      {
         Diagonal_Implicit_Evolution_Operator *g =
            dynamic_cast<Diagonal_Implicit_Evolution_Operator *>(f);

         g->SetTime(t);
         g->ImplicitSolve(dt, x, k); // solve for k: k = f(x + dt*k, t)
         x.Add(dt, k);
         t += dt;
      }
   };

   /// Implicit midpoint method. A-stable, not L-stable.
   class ImplicitMidpointSolver : public ODESolver
   {
   protected:
      Vector k;

   public:
      virtual void Init(TimeDependentOperator &_f)
      {
         ODESolver::Init(_f);
         k.SetSize(f->Size());
      }

      virtual void Step(Vector &x, double &t, double &dt)
      {
         Diagonal_Implicit_Evolution_Operator *g =
            dynamic_cast<Diagonal_Implicit_Evolution_Operator *>(f);

         g->SetTime(t + dt/2);
         g->ImplicitSolve(dt/2, x, k);
         x.Add(dt, k);
         t += dt;
      }
   };

   /** Two stage, singly diagonal implicit Runge-Kutta (SDIRK) methods;
       the choices for gamma_opt are:
       0 - 3rd order method, not A-stable
       1 - 3rd order method, A-stable, not L-stable (default)
       2 - 2nd order method, L-stable
       3 - 2nd order method, L-stable (has solves ouside [t,t+dt]). */
   class SDIRK23Solver : public ODESolver
   {
   protected:
      double gamma;
      Vector k, y;

   public:
      SDIRK23Solver(int gamma_opt = 1)
      {
         if (gamma_opt == 0)
            gamma = (3. - sqrt(3.))/6.; // not A-stable, order 3
         else if (gamma_opt == 2)
            gamma = (2. - sqrt(2.))/2.; // L-stable, order 2
         else if (gamma_opt == 3)
            gamma = (2. + sqrt(2.))/2.; // L-stable, order 2
         else
            gamma = (3. + sqrt(3.))/6.; // A-stable, order 3
      }

      virtual void Init(TimeDependentOperator &_f)
      {
         ODESolver::Init(_f);
         k.SetSize(f->Size());
         y.SetSize(f->Size());
      }

      virtual void Step(Vector &x, double &t, double &dt)
      {
         Diagonal_Implicit_Evolution_Operator *g =
            dynamic_cast<Diagonal_Implicit_Evolution_Operator *>(f);

         // with a = gamma:
         //   a   |   a
         //  1-a  |  1-2a  a
         // ------+-----------
         //       |  1/2  1/2
         // note: with gamma_opt=3, both solve are outside [t,t+dt] since a>1
         g->SetTime(t + gamma*dt);
         g->ImplicitSolve(gamma*dt, x, k);
         add(x, (1.-2.*gamma)*dt, k, y); // y = x + (1-2*gamma)*dt*k
         x.Add(dt/2, k);

         g->SetTime(t + (1.-gamma)*dt);
         g->ImplicitSolve(gamma*dt, y, k);
         x.Add(dt/2, k);
         t += dt;
      }
   };

   /** Three stage, singly diagonal implicit Runge-Kutta (SDIRK) method of
       order 4. A-stable, not L-stable. */
   class SDIRK34Solver : public ODESolver
   {
   protected:
      Vector k, y, z;

   public:
      virtual void Init(TimeDependentOperator &_f)
      {
         ODESolver::Init(_f);
         k.SetSize(f->Size());
         y.SetSize(f->Size());
         z.SetSize(f->Size());
      }

      virtual void Step(Vector &x, double &t, double &dt)
      {
         Diagonal_Implicit_Evolution_Operator *g =
            dynamic_cast<Diagonal_Implicit_Evolution_Operator *>(f);

         //   a   |    a
         //  1/2  |  1/2-a    a
         //  1-a  |   2a    1-4a   a
         // ------+--------------------
         //       |    b    1-2b   b
         // note: two solves are outside [t,t+dt] since c1=a>1, c3=1-a<0
         const double a = 1./sqrt(3.)*cos(M_PI/18.) + 0.5;
         const double b = 1./(6.*(2.*a-1.)*(2.*a-1.));

         g->SetTime(t + a*dt);
         g->ImplicitSolve(a*dt, x, k);
         add(x, (0.5-a)*dt, k, y);
         add(x,  (2.*a)*dt, k, z);
         x.Add(b*dt, k);

         g->SetTime(t + dt/2);
         g->ImplicitSolve(a*dt, y, k);
         z.Add((1.-4.*a)*dt, k);
         x.Add((1.-2.*b)*dt, k);

         g->SetTime(t + (1.-a)*dt);
         g->ImplicitSolve(a*dt, z, k);
         x.Add(b*dt, k);
         t += dt;
      }
   };

   /** Three stage, singly diagonal implicit Runge-Kutta (SDIRK) method of
       order 3. L-stable. */
   class SDIRK33Solver : public ODESolver
   {
   protected:
      Vector k, y;

   public:
      virtual void Init(TimeDependentOperator &_f)
      {
         ODESolver::Init(_f);
         k.SetSize(f->Size());
         y.SetSize(f->Size());
      }

      virtual void Step(Vector &x, double &t, double &dt)
      {
         Diagonal_Implicit_Evolution_Operator *g =
            dynamic_cast<Diagonal_Implicit_Evolution_Operator *>(f);

         //   a  |   a
         //   c  |  c-a    a
         //   1  |   b   1-a-b  a
         // -----+----------------
         //      |   b   1-a-b  a
         const double a = 0.435866521508458999416019;
         const double b = 1.20849664917601007033648;
         const double c = 0.717933260754229499708010;

         g->SetTime(t + a*dt);
         g->ImplicitSolve(a*dt, x, k);
         add(x, (c-a)*dt, k, y);
         x.Add(b*dt, k);

         g->SetTime(t + c*dt);
         g->ImplicitSolve(a*dt, y, k);
         x.Add((1.-a-b)*dt, k);

         g->SetTime(t + dt);
         g->ImplicitSolve(a*dt, x, k);
         x.Add(a*dt, k);
         t += dt;
      }
   };

   /// Abstract class for time dependent coefficient
   class TimeDependentCoefficient : public Coefficient
   {
   protected:
      double time;

   public:
      TimeDependentCoefficient() { time = 0.; }
      void SetTime(double t) { time = t; }
      double GetTime() { return time; }
   };

   /// Abstract class for time dependent vector coefficient
   class TimeDependentVectorCoefficient : public VectorCoefficient
   {
   protected:
      double time;

   public:
      TimeDependentVectorCoefficient(int vd)
         : VectorCoefficient(vd) { time = 0.; }
      void SetTime(double t) { time = t; }
      double GetTime() { return time; }
   };

   /// Class for time dependent C-function coefficient
   class TimeDepFunctionCoefficient : public TimeDependentCoefficient
   {
   protected:
      double (*F)(Vector &, double);

   public:
      TimeDepFunctionCoefficient( double (*f)(Vector &, double) ) : F(f) { }

      virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip)
      {
         double x[3];
         Vector p(x, 3);

         T.Transform(ip, p);

         return (*F)(p, GetTime());
      }

      virtual void Read(istream &in) { }
   };

   class TimeDepVectorFunctionCoefficient
      : public TimeDependentVectorCoefficient
   {
   protected:
      void (*F)(const Vector &, double, Vector &);

   public:
      TimeDepVectorFunctionCoefficient(
         int vd, void (*f)(const Vector &, double, Vector &))
         : TimeDependentVectorCoefficient(vd), F(f) { }

      using VectorCoefficient::Eval;
      virtual void Eval(Vector &V, ElementTransformation &T,
                        const IntegrationPoint &ip)
      {
         double x[3];
         Vector p(x, 3);

         T.Transform(ip, p);

         V.SetSize(vdim);
         (*F)(p, GetTime(), V);
      }
   };
}

using namespace mfem;


// Right-hand side for the set of scalar ODEs du_i/dt = f(u_i,t)
class ScalarODE : public Diagonal_Implicit_Evolution_Operator
// class ScalarODE : public TimeDependentOperator
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
         return 4*pow(t,3)-3*pow(t,2)+2*t-1; // f = 4t^3-3t^2+2t-1
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
   ScalarODE(int _size, int _option=2) { size = _size; option=_option; }

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
// where A and M are time-independent.
class LinearSystemODE : public Diagonal_Implicit_Evolution_Operator
//class LinearSystemODE : public TimeDependentOperator
{
private:
   HypreParMatrix *A, *M, *B; // B = M - dt A, for ImplicitSolve
   ParLinearForm *b_form;
   TimeDependentCoefficient *b_coeff;
   TimeDependentVectorCoefficient *b_vcoeff;
   mutable HypreParVector *b;
   HypreBoomerAMG *amg;
   HyprePCG *pcg;
   HypreParVector *X, *Y, *Z;
   double current_dt; // dt used in B
   HypreBoomerAMG *B_amg;
   HyprePCG *B_pcg;

public:
   LinearSystemODE(HypreParMatrix *_A,
                   HypreParMatrix *_M = NULL,
                   ParLinearForm *_b_form = NULL,
                   TimeDependentCoefficient *_b_coeff = NULL,
                   TimeDependentVectorCoefficient *_b_vcoeff = NULL)
      : A(_A), M(_M), b_form(_b_form), b_coeff(_b_coeff), b_vcoeff(_b_vcoeff)
   {
      size = A->Size();
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

      B = NULL; // Allocated in ImplicitSolve if needed
      current_dt = -1.;
      B_amg = NULL;
      B_pcg = NULL;
   }

   void AssembleBVector() const
   {
      if (b_form)
      {
         // set the time for the linear form ?
         // b_form->SetTime(GetTime());     // ?

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
      X->SetData(x.GetData());
      Y->SetData(y.GetData());

      AssembleBVector();

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

   /** Solve for k in the equation: M k = A x + dt A k + b.
       This method is needed for backward Euler and other DIRK methods. */
   virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k)
   {
      if (!M)
      {
         cerr << "Not implemented!" << endl;
         abort();
      }

      // Make X and Y wrappers for x and k, respectively.
      X->SetData(x.GetData());
      Y->SetData(k.GetData());

      AssembleBVector();

      // 1) Form the matrix: B = M - dt A
      if (!B)
         B = new HypreParMatrix(hypre_ParCSRMatrixAdd(*M, *A));
      if (dt != current_dt)
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

      // 3) Solve the system: B Y = Z
      B_pcg->Mult(*Z, *Y);
   }

   virtual ~LinearSystemODE()
   {
      delete B_pcg;
      delete B_amg;
      delete B;
      delete b;
      delete Z;
      delete pcg;
      delete amg;
      delete Y;
      delete X;
   }
};


// Solve the given ODE with the given initial condtion using uniform
// time-stepping and (optionally) visualize the solution.
void SolveODE(TimeDependentOperator *ode, HypreParVector *X0,
              ODESolver *solver, double tstart = 0.0, double tstop = 1.0,
              int ntime = 100, ParGridFunction *x = NULL,
              TimeDependentCoefficient *exsol = NULL)
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
      *sol_sock << "autoscale value\n";
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

// Wrapper for WARP's status object
class WarpStatus
{
   private:
      warp_Status status;
   
   public:
      WarpStatus(warp_Status _status)
      {
         status = _status;
      }

      void GetDone(warp_Int *done_ptr) { warp_GetStatusDone(status, done_ptr); }
      void GetLevel(warp_Int *level_ptr) { warp_GetStatusLevel(status, level_ptr); }
      void GetIter(warp_Int *iter_ptr) { warp_GetStatusIter(status, iter_ptr); }
      void GetResidual(warp_Real *rnorm_ptr) { warp_GetStatusResidual(status, rnorm_ptr); }
      
      // The warp_Status structure is deallocated inside of Warp
      // This class is just to make code consistently look object oriented
      ~WarpStatus() { }
};


// Wrapper for WARP's App object
class WarpApp
{
public:
   TimeDependentOperator *ode;
   ODESolver *solver;

   TimeDependentCoefficient *exact_sol;

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
           HypreParVector *_X0, ParGridFunction *_x, ODESolver *_solver,
           double _tstart = 0.0, double _tstop = 1.0, int _ntime = 100)
      : ode(_ode), solver(_solver), comm_t(_comm_t), tstart(_tstart),
        tstop(_tstop), ntime(_ntime), X0(_X0), x(_x)
   {
      buff_size = sizeof(double) * X0->Size();

      solver->Init(*ode);
      *x = *X0;
      pmesh = x->ParFESpace()->GetParMesh();
      sol_sock = NULL;
      exact_sol = NULL;
   }

   void SetExactSolution(TimeDependentCoefficient *exsol)
   { exact_sol = exsol; }

   ~WarpApp()
   {
      delete sol_sock;
   }

   // Below warp_Vector == HypreParVector* and warp_App == WarpApp*

   static int Phi(warp_App     _app,
                  double       tstart,
                  double       tstop,
                  double       accuracy,
                  warp_Vector  _u,
                  int         *rfactor_ptr)
   {
      WarpApp *app = (WarpApp*) _app;
      HypreParVector *u = (HypreParVector*) _u;

      double t = tstart;
      double dt = tstop-tstart;

      app->solver->Step(*u, t, dt);

      // no refinement
      *rfactor_ptr = 1;

      return 0;
   }

   static int Clone(warp_App     _app,
                    warp_Vector  _u,
                    warp_Vector *v_ptr)
   {
      HypreParVector *u = (HypreParVector*) _u;
      HypreParVector *v = new HypreParVector(*u);
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
         // u->Randomize(2142);
         *u = 0.;
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
                    warp_Status  _status,
                    warp_Vector  _u)
   {
      WarpApp *app = (WarpApp*) _app;
      WarpStatus status = WarpStatus(_status);
      HypreParVector *u = (HypreParVector*) _u;
      
      // Extract information from status
      int done, level, iter;
      double rnorm;
      status.GetDone(&done);
      status.GetLevel(&level);
      status.GetIter(&iter);
      status.GetResidual(&rnorm); 

      // if (t == app->tstart || t == app->tstop)
      if ( (t == app->tstop) && (level == 0) )
      {
         (*app->x) = *u; // Distribute

         // Opening multiple 'parallel' connections to GLVis simultaneously
         // (from different time intervals in this case) may cause incorrect
         // behavior.

         ParMesh *pmesh = app->pmesh;

         char vishost[] = "localhost";
         int  visport   = 19916;
         int  init_sock = 0;
         if (!app->sol_sock)
         {
            app->sol_sock = new socketstream(vishost, visport);
            init_sock = 1;
         }

         socketstream &sol_sock = *app->sol_sock;

         int good, all_good;
         good = sol_sock.good();
         MPI_Allreduce(&good, &all_good, 1, MPI_INT, MPI_LAND,
                       pmesh->GetComm());

         if (all_good)
         {
            sol_sock << "parallel " << pmesh->GetNRanks()
                     << " " << pmesh->GetMyRank() << "\n";
            sol_sock << "solution\n";
            sol_sock.precision(8);
            pmesh->Print(sol_sock);
            app->x->Save(sol_sock);

            if (init_sock)
            {
               int comm_t_rank;
               MPI_Comm_rank(app->comm_t, &comm_t_rank);

               // sol_sock << "valuerange 0 1\n";
               // sol_sock << "autoscale off\n";
               sol_sock << "keys cmAaa\n";
               sol_sock << "window_title 'comm_t rank: " << comm_t_rank
                        << ", t = " << t << "'\n";
            }
            sol_sock << flush;
            if (pmesh->GetMyRank() == 0)
               cout << "Solution updated." << flush;
         }

         if (app->exact_sol)
         {
            app->exact_sol->SetTime(t);

            double err = app->x->ComputeL2Error(*app->exact_sol);

            if (pmesh->GetMyRank() == 0)
               cout << " L2 norm of the error = " << err << endl;
         }
         else if (pmesh->GetMyRank() == 0)
            cout << endl;
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

// Wrapper for WARP's testing routines
class WarpTest
{
private:

public:
   
   // Empty constructor
   WarpTest( ){ }
   
   // Test Function for Init and Write function
   void TestInitWrite( WarpApp              *app,
                       MPI_Comm              comm_x,
                       double                t,
                       warp_PtFcnInit        init,
                       warp_PtFcnWrite       write,
                       warp_PtFcnFree        free)
   { warp_TestInitWrite((warp_App) app, comm_x, t, init, write, free); }

   // Test Function for Clone 
   void TestClone( WarpApp              *app,
                   MPI_Comm              comm_x,
                   double                t,
                   warp_PtFcnInit        init,
                   warp_PtFcnWrite       write,
                   warp_PtFcnFree        free,
                   warp_PtFcnClone       clone)
   { warp_TestClone((warp_App) app, comm_x, t, init, write, free, clone); }
   
   // Test Function for Sum 
   void TestSum( WarpApp              *app,
                 MPI_Comm              comm_x,
                 double                t,
                 warp_PtFcnInit        init,
                 warp_PtFcnWrite       write,
                 warp_PtFcnFree        free,
                 warp_PtFcnClone       clone,
                 warp_PtFcnSum         sum)
   { warp_TestSum((warp_App) app, comm_x, t, init, write, free, clone, sum); }
   
   // Test Function for Dot 
   void TestDot( WarpApp              *app,
                 MPI_Comm              comm_x,
                 double                t,
                 warp_PtFcnInit        init,
                 warp_PtFcnFree        free,
                 warp_PtFcnClone       clone,
                 warp_PtFcnSum         sum,
                 warp_PtFcnDot         dot,
                 int                  *correct)
   { warp_TestDot((warp_App) app, comm_x, t, init, free, clone, sum, dot, correct); }

   // Test Functions BufSize, BufPack, BufUnpack
   void TestBuf( WarpApp              *app,
                 MPI_Comm              comm_x,
                 double                t,
                 warp_PtFcnInit        init,
                 warp_PtFcnFree        free,
                 warp_PtFcnSum         sum,  
                 warp_PtFcnDot         dot,
                 warp_PtFcnBufSize     bufsize,
                 warp_PtFcnBufPack     bufpack,
                 warp_PtFcnBufUnpack   bufunpack,
                 int                  *correct)
   { warp_TestBuf((warp_App) app, comm_x, t, init, free, sum, dot, bufsize, bufpack, bufunpack, correct); }

   // Test Functions Coarsen and Refine
   void TestCoarsenRefine(WarpApp          *app,
                          MPI_Comm          comm_x,
                          double            t,
                          double            f_tminus,
                          double            f_tplus,
                          double            c_tminus,
                          double            c_tplus,
                          warp_PtFcnInit    init,
                          warp_PtFcnWrite   write,
                          warp_PtFcnFree    free,
                          warp_PtFcnClone   clone,
                          warp_PtFcnSum     sum,
                          warp_PtFcnDot     dot,
                          warp_PtFcnCoarsen coarsen,
                          warp_PtFcnRefine  refine,
                          warp_Int         *correct)
   { warp_TestCoarsenRefine( (warp_App) app, comm_x, t, f_tminus, f_tplus, c_tminus, c_tplus, init,
                            write, free, clone, sum, dot, coarsen, refine, correct); }


   ~WarpTest() { }

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

   void SetNRelax(int level, int nrelax)
   { warp_SetNRelax(core, level, nrelax); }

   void SetAbsTol(double tol) { warp_SetAbsTol(core, tol); }

   void SetCFactor(int level, int cfactor)
   { warp_SetCFactor(core, level, cfactor); }

   /** Use cfactor0 on all levels until there are < cfactor0 points
       on each processor. */
   void SetAggCFactor(int cfactor0)
   {
      WarpApp *app = (WarpApp *) core->app;
      int nt = app->ntime, pt;
      MPI_Comm_size(app->comm_t, &pt);
      if (cfactor0 > -1)
      {
         int level = (int) (log10((nt + 1) / pt) / log10(cfactor0));
         for (int i = 0; i < level; i++)
            warp_SetCFactor(core, i, cfactor0);
      }
   }

   void SetMaxIter(int max_iter) { warp_SetMaxIter(core, max_iter); }
   
   void SetPrintLevel(int print_level) { warp_SetPrintLevel(core, print_level); }
   
   void SetWriteLevel(int write_level) { warp_SetWriteLevel(core, write_level); }

   void SetFMG() { warp_SetFMG(core); }

   void Drive() { warp_Drive(core); }

   ~WarpCore() { warp_Destroy(core); }
};


// Exact solution parameters
const double tau = (2.+1./6.)*M_PI;
const double kap = M_PI;

// Exact solution
double ExSol(Vector &p, double t)
{
   int dim = p.Size();

   if (dim == 2)
      return sin(kap*p(0))*sin(kap*p(1))*sin(tau*t);

   return sin(kap*p(0))*sin(kap*p(1))*sin(kap*p(2))*sin(tau*t);
}

// Initial condition
double IC(Vector &x)
{
   // return exp(-5.*(x(0)*x(0)+x(1)*x(1)));
   return ExSol(x, 0.);
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
      return;
   }

   v(0) = kap*cos(kap*p(0))*sin(kap*p(1))*sin(kap*p(2))*sin(tau*t);
   v(1) = kap*sin(kap*p(0))*cos(kap*p(1))*sin(kap*p(2))*sin(tau*t);
   v(2) = kap*sin(kap*p(0))*sin(kap*p(1))*cos(kap*p(2))*sin(tau*t);
}

// Source term
double ST(Vector &p, double t)
{
   int dim = p.Size();

   if (dim == 2)
      return (tau*cos(tau*t) + 2*kap*kap*sin(tau*t))*
         sin(kap*p(0))*sin(kap*p(1));

   return (tau*cos(tau*t) + 3*kap*kap*sin(tau*t))*
      sin(kap*p(0))*sin(kap*p(1))*sin(kap*p(2));
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

   // Variables used by WarpTest
   int correct;
   double dt;

   // Default parameters:
   const char *meshfile = "../../mfem/data/star.mesh";
   int    sref        = 1;
   int    pref        = 1;
   int    use_warp    = 1;
   int    num_procs_x = 1;
   double tstart      = 0.0;
   double tstop       = 1.0;
   int    ntime       = 32;
   // double cfl         = 1.0;
   int    ode_solver  = -1;

   int heat_equation = 1;
   int scalar_ode_option = 2;

   // WARP default parameters:
   int    max_levels  = 10;
   int    nrelax      = 1;
   int    nrelax0     = -1;
   double tol         = 1e-9;
   int    cfactor     = 2;
   int    cfactor0    = -1;
   int    max_iter    = 100;
   int    fmg         = 0;
   int    write_level = 0;
   bool   run_wrapper_tests = false;

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
      else if (strcmp(argv[arg_index], "-nowarp") == 0)
      {
         use_warp = 0;
      }
      else if (strcmp(argv[arg_index], "-px") == 0)
      {
         num_procs_x = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-tstart") == 0)
      {
         tstart = atof(argv[++arg_index]);
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
      else if (strcmp(argv[arg_index], "-odesolver") == 0)
      {
         ode_solver = atoi(argv[++arg_index]);
      }
      else if (strcmp(argv[arg_index], "-ml") == 0)
      {
         max_levels = atoi(argv[++arg_index]);
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
      }
      else if (strcmp(argv[arg_index], "-run_wrapper_tests") == 0)
      {
         run_wrapper_tests = true;
      }
      else if (strcmp(argv[arg_index], "-write") == 0)
      {
         write_level = 1;
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
         "  -run_wrapper_tests: Only run the Warp wrapper tests\n"
         "                      (do not combine with temporal parallelism)\n"
         "  -sref <num>       : levels of serial refinements (default: 1)\n"
         "  -pref <num>       : levels of parallel refinements (default: 1)\n"
         "  -nowarp           : sequential time integration (default: off)\n"
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
         "  -ml  <max_levels> : set max number of time levels (default: 10)\n"
         "  -nu  <nrelax>     : set num F-C relaxations (default: 1)\n"
         "  -nu0 <nrelax>     : set num F-C relaxations on level 0\n"
         "  -tol <tol>        : set stopping tolerance (default: 1e-9)\n"
         "  -cf  <cfactor>    : set coarsening factor (default: 2)\n"
         "  -cf0 <cfactor0>   : set aggressive coarsening (default: off)\n"
         "  -mi  <max_iter>   : set max iterations (default: 100)\n"
         "  -fmg              : use FMG cycling\n"
         "  -write            : set write_level (default: 0) \n"
         "\n";
   }

   if (print_usage)
   {
      MPI_Finalize();
      return 0;
   }

   if (!use_warp)
      num_procs_x = num_procs;

   int xcolor = myid / num_procs_x;
   int tcolor = myid % num_procs_x;
   MPI_Comm_split(comm, xcolor, myid, &comm_x);
   MPI_Comm_split(comm, tcolor, myid, &comm_t);

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
   ConstantCoefficient minus_one(-1.0);
   ConstantCoefficient one(1.0);
   TimeDepFunctionCoefficient source(&ST);
   TimeDepVectorFunctionCoefficient nbc(dim, &NBC);
   TimeDepFunctionCoefficient exact_sol(&ExSol);

   ParBilinearForm *a = new ParBilinearForm(fespace);
   ParBilinearForm *m = new ParBilinearForm(fespace);
   ParLinearForm   *b = new ParLinearForm(fespace);

   a->AddDomainIntegrator(new DiffusionIntegrator(minus_one));
   m->AddDomainIntegrator(new MassIntegrator(one));
   b->AddDomainIntegrator(new DomainLFIntegrator(source));
   // b->AddBoundaryIntegrator(new BoundaryLFIntegrator(one));
   b->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(nbc));

   // Local assembly. Use precomputed sparsity to ensure that 'a' and 'm' have
   // the same sparsity patterns.
   a->UsePrecomputedSparsity(); a->Assemble(); a->Finalize();
   m->UsePrecomputedSparsity(); m->Assemble(); m->Finalize();

   HypreParMatrix *A = a->ParallelAssemble();
   HypreParMatrix *M = m->ParallelAssemble();

   delete a;
   delete m;

   // Define the righ-hand side of the ODE and call MFEM or WARP to solve it.
   TimeDependentOperator *ode;
   if (heat_equation)
      ode = new LinearSystemODE(A, M, b, &source, &nbc);
   else
      ode = new ScalarODE(A->Size(), scalar_ode_option);

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

   if (!use_warp)
   {
      if (heat_equation)
         SolveODE(ode, X0, solver, tstart, tstop, ntime, &x0, &exact_sol);
      else
         SolveODE(ode, X0, solver, tstart, tstop, ntime, &x0);
   }
   else
   {
      WarpApp app(comm_t, ode, X0, &x0, solver, tstart, tstop, ntime);

      if(run_wrapper_tests)
      {
         // Simple tests for the wrappers 
         WarpTest test;
         dt = (app.tstop - app.tstart)/ (double) app.ntime;

         // Test init(), write(), free()
         test.TestInitWrite( &app, comm_x, 0.0, WarpApp::Init, 
                             WarpApp::Write, WarpApp::Free);
         cout << endl << "Press Enter to continue " << endl;
         cin.get();
         test.TestInitWrite( &app, comm_x, dt, WarpApp::Init, 
                             WarpApp::Write, WarpApp::Free);
         cout << endl << "Press Enter to continue " << endl;
         cin.get();

         // Test clone()
         test.TestClone( &app, comm_x, 0.0, WarpApp::Init, 
                         WarpApp::Write, WarpApp::Free, 
                         WarpApp::Clone);
         cout << endl << "Press Enter to continue " << endl;
         cin.get();
         test.TestClone( &app, comm_x, dt, WarpApp::Init, 
                         WarpApp::Write, WarpApp::Free, 
                         WarpApp::Clone);
         cout << endl << "Press Enter to continue " << endl;
         cin.get();

         // Test sum() 
         test.TestSum( &app, comm_x, 0.0, WarpApp::Init, 
                       WarpApp::Write, WarpApp::Free, 
                       WarpApp::Clone, WarpApp::Sum);
         cout << endl << "Press Enter to continue " << endl;
         cin.get();
         test.TestSum( &app, comm_x, dt, WarpApp::Init, 
                       WarpApp::Write, WarpApp::Free, 
                       WarpApp::Clone, WarpApp::Sum);
         cout << endl << "Press Enter to continue " << endl;
         cin.get();

         // Test dot()
         test.TestDot( &app, comm_x, 0.0, WarpApp::Init, WarpApp::Free, 
                       WarpApp::Clone, WarpApp::Sum, WarpApp::Dot, 
                       &correct);
         cout << endl << "Press Enter to continue " << endl;
         cin.get();
         test.TestDot( &app, comm_x, dt, WarpApp::Init, WarpApp::Free,
                       WarpApp::Clone, WarpApp::Sum, WarpApp::Dot, 
                       &correct);
         cout << endl << "Press Enter to continue " << endl;
         cin.get();

         // Test bufsize(), bufpack(), bufunpack()
         test.TestBuf( &app, comm_x, 0.0, WarpApp::Init, WarpApp::Free, 
                       WarpApp::Sum, WarpApp::Dot, WarpApp::BufSize, 
                       WarpApp::BufPack, WarpApp::BufUnpack, &correct);
         cout << endl << "Press Enter to continue " << endl;
         cin.get();
         test.TestBuf( &app, comm_x, dt, WarpApp::Init, WarpApp::Free,
                        WarpApp::Sum, WarpApp::Dot, WarpApp::BufSize, 
                        WarpApp::BufPack, WarpApp::BufUnpack, &correct);
         cout << endl << "Press Enter to continue " << endl;
         cin.get();
      
      }
      else
      {
         // Run a warp simulation
         WarpCore core(comm, &app);

         if (heat_equation)
            app.SetExactSolution(&exact_sol);

         core.SetWriteLevel(write_level);
         core.SetPrintLevel(1);
         core.SetMaxLevels(max_levels);
         core.SetNRelax(-1, nrelax);
         if (nrelax0 > -1)
            core.SetNRelax(0, nrelax0);
         core.SetAbsTol(tol);
         core.SetCFactor(-1, cfactor);
         core.SetAggCFactor(cfactor0);
         core.SetMaxIter(max_iter);
         if (fmg)
            core.SetFMG();

         core.Drive();
      }
   }

   // Free memory.
   delete solver;
   delete ode;

   delete X0;
   delete M;
   delete A;
   delete b;

   delete fespace;
   if (!pmesh->GetNodes())
      delete fec;
   delete pmesh;

   MPI_Comm_free(&comm_x);
   MPI_Comm_free(&comm_t);

   MPI_Finalize();

   return 0;
}
