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

/* Drive 06: 2D nonlinear p-laplacian equation. Requires MFEM, Hpyre, Metis and GlVis
   
   Interface: MFEM 

   Compile with: make drive-06 -- Modify Makefile to point to metis, mfem and hypre libraries

   Help with:    drive-06 -help

   Sample run:   mpirun -n12 ./drive-06

   Description:  This code solves the 2D p -laplacian problem 
      
                                     u_t - div( |grad(u)|^2 grad(u) ) = b,
                                                 u(x,0) = U0(x)                                  
                                 |grad(u)|^2 grad(u) . n = g(x,t) on boundary
            
                 on a given mesh using MFEM. The coressponding weak form is
 
                     <u_t,v> + <|grad(u)|^2 grad(u), grad(v) = <b,v> + <g,v>_b = <f,v>
                 
                 Implicit methods use newtons method. Consider BDF1. Let u_k = u(t_k) and define

                               L(u)(v) = <u,v> + dt*<|grad(u)|^2 grad(u),grad(v)>
                                       f(v) = <u_k + dt*b, v> + dt*<g,v>_b
                L'(u)(v)[w] = <w,v> + dt*<(|grad(u)|^2 + 2 grad(u)grad(u)^t ) grad(w), grad(v)>  
                                       
                where L' is the Frechet derivitive of L. Each implicit time step solves

                                            L(u_{k+1})(v) - f(v) = 0 
                
                Letting superscript denote iteration number, Newtons method is

                            u^{j+1} = u^j -  L'(u^j)(v)[L(u^j)(v) - f(v)] = u^j - m^j
                
                where the previous time step value is used for u^0 and m^j solves
                the linear equaiton 
                              
                              L'(u^j)(v)[m^j] = L(u_j)(v) - f(v)  
                               
                Hypre BOOMERamg is used as the linear solver */

#include "braid_mfem.hpp"
#include <fstream>
#include <iostream>

//Extra Hypre Functions

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

using namespace std;
   
/**************************************************************************************
Return the coefficient in the diffusion term of the weak form
               
               <|nl_coeff(u,t)*grad(u),grad(v)> 
       nl_coeff(u,t) = -a*grad(u)|^2*I - 2*b*grad(u)*grad(u)^t

where:    1. a = 1, b = 0 is is the coeff for the nolinear diffusion term 
          2. a = b = dt is the coeff for the Newton Linearization diffusion term
*/
class NonlinearCoefficient : public MatrixCoefficient
{
protected:
    DenseMatrix fd;
    int dim;
    double power, a, b, *data;
    ParGridFunction *gfunc;
    Vector grad;
   
public:
    NonlinearCoefficient(ParGridFunction *_gfunc, double _power, int _dim);
 
    void SetGridFunction( HypreParVector &X, double _a, double _b);
    void Eval(DenseMatrix &K, ElementTransformation &T, const IntegrationPoint &ip);

    virtual void Read(istream &in) { }

    virtual ~NonlinearCoefficient();
};

/*************************************************************************************
Class that contains all coefficients and variables relating to the manufactured solution */
class Problem
{
protected:
   static double kap;  //Spatial frequency of manufactured solution
   static double tau;  //Temporal frequency of manufactured solution
   static double power; // p in p-laplacian equation
   static double t_start;
 
public:
	
   ConstantCoefficient one;          //One 
   FunctionCoefficient *source;      //Source Term
   VectorFunctionCoefficient *nbc; //Boundary Conditions
   FunctionCoefficient *exact_sol;    //Exact Solution
   FunctionCoefficient *ic;			 // Initial Condition 
   Array<NonlinearCoefficient *> nlc;   // Nonlinear Coefficient
 

   Problem(double _t_start, double _kap, double _tau, double _power, int num_space_levels);
   
   //Neuman boundary Condition
	static void NBC(const Vector &p, double t, Vector &v);

	//Exact Solution
	static double ExSol(Vector &p, double t);

	//Source term
	static double ST(Vector &p, double t);

	//Initial Condition
	static double IC(Vector &x);

   virtual ~Problem();
};

/*************************************************************************************
Struct to hold all command line and defualt options */
struct NonlinearOptions : public BraidOptions
{  

	int mi_fine, mi_coarse;    //Max nonlinear iterations on f/c grids
	double tol_fine, tol_coarse; //Nonlinear tolerance on f/c grids 
   int ode_solver_type;         //What type of ode solver to use
	
   int nl_print_level; //Nonlinear print level
   int sc_print_level; //Spatial coarsening print level
   
   int step_type;   //Type of step being completed
   int braid_level; //Current temporal level of braid
   int braid_iter;  //Current braid Iteration
   
   double kap, tau; //Exact Solution parameters    
   double power;    //P in p-laplacian 
   double dt;       //time step on the finest grid.
   double dx;       //Space step on fine grid.
   double cg_dt;    //Coarse grid time step;
   //Visualization Options
   const char *vishost;  
   int         visport;
   int vis_time_steps;   // Visualize every x time steps
   int vis_braid_steps;  // Visualize every x braid iterations 
   bool vis_screenshots; // Save Screenshots t
   Problem *prob;       // All of the problem specific info
   HypreParVector *u;  //Stores a pointer to the initial guess for nonlinear solver
   BraidUtil util;     // Braid utility for tolerance scaling. 

   NonlinearOptions(int argc, char *argv[]);

   virtual ~NonlinearOptions();
};

/*************************************************************************************
A time dependent Operator that solves the nonlinear system ode using newtons method */
class NonlinearSystemODE : public TimeDependentOperator
{
private:


	ParBilinearForm *a;	            
	ParLinearForm *b_form;

  
   mutable HypreParVector *b;
	mutable HypreParMatrix *A;	
	HypreParMatrix *B, *M;
   Coefficient *b_coeff;
   VectorCoefficient *b_vcoeff;
   NonlinearOptions *opts;
   NonlinearCoefficient *nl_coeff;
   HypreParVector *X, *Y, *Z;
	HypreParVector *U, *V, *W;
	HypreBoomerAMG *amg;
   HyprePCG *pcg;
   HyprePCG *B_pcg;
	HypreBoomerAMG *B_amg;	
   int braid_iter, step_type;
   int v,newt;
public:
	NonlinearSystemODE(int              l,
                      ParBilinearForm *_a,
							 HypreParMatrix  *_M,
							 ParLinearForm   *_b_form,
							 NonlinearOptions *_opts);
	 

	 void AssembleDiffusionMatrix() const;	 
	 
    // Assemble the rhs
    void AssembleBVector() const;

    //Update the RHS and calculate residual
    void UpdateResidual(double dt, double &resid);

    //Do a line search for the newton solver
    int LineSearch(double dt, double &resid, double a);

    //Nonlinear Multiplication -- Not implimented yet -- Not needed??
    virtual void Mult(const Vector &x, Vector &y) const;
    
    //Nonlinear Solve using Newtons Method 
    virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k);
	 
    //Set the hypre parameters
    void SetHypre(const double dt,const double tol);

    virtual ~NonlinearSystemODE(); 
};

/*************************************************************************************
Class defining the NonlinearApp */ 
class NonlinearApp : public MFEMBraidApp
{
protected:

	NonlinearOptions *opts;
	FiniteElementCollection *fe_coll;
   Problem *prob;
   int braid_iter; 
 
   
   //Allocate Arrays needed by the NonlinearApp
   virtual void AllocLevels(int num_levels);
	
   //Construct the Meshes on each level
   virtual ParFiniteElementSpace *ConstructFESpace(ParMesh *pmesh);

   //Init the ode/solver/max_dt for each level
	virtual void InitLevel(int l);            

public:
   SpaceTimeMeshInfo MeshInfo;	
   
   NonlinearApp(NonlinearOptions *opts, MPI_Comm comm_t_, ParMesh *pmesh);
	
   //Modified MFEMBraidApp::Phi to also record space-time mesh info
   virtual int Step(braid_Vector    u_,
                   braid_Vector    ustop_,
                   braid_Vector    fstop_,
                   BraidStepStatus &pstatus);  

   virtual ~NonlinearApp();

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
	NonlinearOptions opts(argc, argv);
   
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

   // Split comm (MPI_COMM_WORLD) into spatial and temporal communicators
   MPI_Comm comm_x, comm_t;
   BraidUtil util;
   util.SplitCommworld(&comm, opts.num_procs_x, &comm_x, &comm_t);

   // Partition the mesh accross the spatial communicator
   ParMesh *pmesh = new ParMesh(comm_x, *mesh);
   delete mesh;

   //Init braid App
	NonlinearApp app(&opts, comm_t, pmesh);
   //Run braid
   //
   double dx = 1/pow(2,opts.ser_ref_levels + opts.par_ref_levels);
   
   opts.tol = opts.tol/(opts.dt*dx*dx);
   opts.dx = dx;
   pmesh->PrintInfo();
   BraidCore core(comm, &app);
   opts.SetBraidCoreOptions(core);
   core.Drive();   
   
   if (opts.sc_print_level && myid == 0)
      app.MeshInfo.Print(comm);
 
   MPI_Finalize();
   return 0;
}	

//*************************************************************************************//
/* Implimentation::NonlinearCoefficient */
NonlinearCoefficient::NonlinearCoefficient(ParGridFunction *_gfunc, double _power, int _dim)
         :MatrixCoefficient(_dim), fd(_dim), dim(_dim), power(_power), gfunc(_gfunc)
{
   a = 1;
   b = 0;
   *gfunc = 1.0;
}

void NonlinearCoefficient::SetGridFunction(HypreParVector &X, double _a, double _b)
{
        *gfunc = X; // distribute (communication)
        a = _a;
        b = _b;
}

void NonlinearCoefficient::Eval(DenseMatrix &K, ElementTransformation &T, const IntegrationPoint &ip)
{
        //-a*<|grad(u)|^2*I*grad(u),grad(v)> - b*<2*grad(u)grad(u)^t grad(u),grad(v)>
        gfunc->GetGradient(T, grad);
        data = grad.GetData();
        double coeff = -a*pow(grad * grad, power/2.) ;
        double coeff1;
        if (power == 0)
           coeff1 = 0;
        else
           coeff1 = -b*power*pow(grad * grad, (power-2.)/2.);
        double aa[4] = {data[0]*data[0],data[0]*data[1],data[0]*data[1],data[1]*data[1]};
        fd = aa;
       
        fd *= coeff1;
        K.Diag(coeff, dim);
        K+=fd;
             
}

NonlinearCoefficient::~NonlinearCoefficient()
{
}

/*************************************************************************************
Implimentaion::Problem */

//Initailize static data members
double Problem::kap = 0;
double Problem::power = 0;
double Problem::tau = 0;
double Problem::t_start = 0;

Problem::Problem(double _t_start, double _kap, double _tau, double _power, int num_space_levels)
                  :  one(1.0), nlc(num_space_levels)
{
      t_start = _t_start;
      kap = _kap; tau = _tau; power = _power; 
      source = new FunctionCoefficient(ST);
      exact_sol = new FunctionCoefficient(ExSol);
      ic = new FunctionCoefficient(IC);
      nbc = new VectorFunctionCoefficient(2,NBC);
}

void Problem::NBC(const Vector &p, double t, Vector &v)
{
      double Ux, Uy;
      Ux = kap*cos(kap*p(0))*sin(kap*p(1))*cos(tau*t);
      Uy = kap*sin(kap*p(0))*cos(kap*p(1))*cos(tau*t);
      v(0) = pow(Ux*Ux + Uy*Uy,power/2)*Ux;
      v(1) = pow(Ux*Ux + Uy*Uy,power/2)*Uy;
      return;
}

double Problem::ExSol(Vector &p, double t)
{
  return sin(kap*p(0))*sin(kap*p(1))*cos(tau*t);
}

double Problem::ST(Vector &p, double t)
{
      double Ut,Ux,Uxx,Uy,Uxy;
      Ut = -tau*sin(kap*p(0))*sin(kap*p(1))*sin(tau*t);
      Ux = kap*cos(kap*p(0))*sin(kap*p(1))*cos(tau*t);
      Uy = kap*sin(kap*p(0))*cos(kap*p(1))*cos(tau*t);
      Uxx = -kap*kap*sin(kap*p(0))*sin(kap*p(1))*cos(tau*t);
      Uxy =  kap*kap*cos(kap*p(0))*cos(kap*p(1))*cos(tau*t);
    
      if (power == 0)
         return  Ut-(power+2.)*pow((Ux*Ux +Uy*Uy),power/2)*Uxx; 
      else
         return (Ut-(power+2.)*pow((Ux*Ux +Uy*Uy),power/2)*Uxx -
                 2*power*(Ux*Uy*Uxy)*pow(Ux*Ux+Uy*Uy,power/2.-1.));
}

double Problem::IC(Vector &x)
{
	return ExSol(x,Problem::t_start);
}

Problem::~Problem() 
{
   delete source;
   delete exact_sol;
   delete ic;
   delete nbc;

   for (int l = 0; l < nlc.Size() ; l++)
   { 
      delete nlc[l];
   }
}
/*************************************************************************************
Implimentation::NonlinearOptions */
NonlinearOptions::NonlinearOptions(int argc, char *argv[])
            : BraidOptions(argc, argv)
{
   //Set inherited default options
	num_time_steps = 32;
   t_start = 0.0;
   t_final = 4.0;
   double a = 1;
	//Set defualts for the mesh/refinement inherited options
	mesh_file			= "/g/g90/oneill10/Programs/mfem/data/inline-quad.mesh";
	ser_ref_levels = 1;
	par_ref_levels = 1;
	AddMeshOptions();

	//Set defualts for the NonlinearOptions specific options
   mi_fine = 100;
	mi_coarse = 100;
	tol_fine = 0.0000001;
	tol_coarse = tol_fine; 
   kap = M_PI;
   tau = (2 + 1./6.)*M_PI; 
   power = 2;
   
   nl_print_level = 0;
   sc_print_level = 0;
	ode_solver_type = -1;
   no_down = 1;
   vishost         = "localhost";
   visport         = 19916;
   vis_time_steps  = 0;
   vis_braid_steps = 1;
   vis_screenshots = false;

   //Parse Command Line
   AddOption(&mi_fine, "-mif", "--max-iter-fine ,",
             "Max Newton iterations on the fine grid");
   AddOption(&mi_coarse, "-mic", "--max-iter-coarse ,",
             "Max Newton iterations on the coarse grids ");
   AddOption(&tol_fine, "-rtf", "--resid-tol-fine ,",
             "Nonlinear Residual tolerance on the fine grid ");
   AddOption(&tol_coarse, "-rtc", "--resid-tol-coarse ,",
             "Nonlinear Residual tolerance on the coarse grid");
   AddOption(&sc_print_level, "-scprint", "--spatial-coarsening-print_level ,",
             "Print level for spatial coarseing"); 
   AddOption(&nl_print_level, "-nlprint", "--nonlinear-print_level ,",
             "Print level for nonlinear sovler");    
   AddOption(&ode_solver_type, "-s", "--ode-solver ,",
             "ODE solver: 1 - Forward Euler, 2 - RK2 SSP, 3 - RK3 SSP,"
             " 4 - RK4, 6 - RK6, -1 - Backward Euler, -2 - SDIRK22 "
             " -3 - SDIRK33 -4 -SDIRK23, -5 - SDIRK34.");
   AddOption(&power, "-pow", "--p-laplacian ,",
             "Value of p in the p-laplacian");
   AddOption(&vishost, "-vh", "--visualization-host ,",
             "Set the GLVis host."); 
   AddOption(&visport, "-vp", "--visualization-port ,",
             "Set the GLVis port.");
   AddOption(&vis_time_steps, "-vts", "--visualize-time-steps ,",
             "Visualize every n-th time step (0:final only).");
   AddOption(&vis_braid_steps, "-vbs", "--visualize-braid-steps ,",
             "Visualize every n-th Braid step (0:final only).");
   AddOption(&vis_screenshots, "-vss", "--vis-screenshots , 1 ,",
             "-no-vss", "--vis-no-screenshots, 0", "Enable/disable the saving of"
             " screenshots of the visualized data.");
          
   Parse();
   
   //Init Problem class 

   if (spatial_coarsen)
      prob = new Problem(t_start, kap, tau, power, par_ref_levels + 1);
   else
      prob = new Problem(t_start, kap, tau, power, 1);
   
   dt = (t_final - t_start) / num_time_steps;
   double nt = num_time_steps;
   double l = 0;
   while (nt/cfactor > min_coarse && l < max_levels)
   { 
      nt = nt/cfactor;
      l = l+1;
   }
   cg_dt = (t_final-t_start)/nt;
}  


NonlinearOptions::~NonlinearOptions()
{
   delete prob;
}

/*************************************************************************************
Implimentation::NonlinearSystemODE */
NonlinearSystemODE::NonlinearSystemODE(int             l,
                                       ParBilinearForm *_a,
							                  HypreParMatrix  *_M,
							                  ParLinearForm   *_b_form,
							                  NonlinearOptions *_opts) :
        a(_a), b_form(_b_form), M(_M), b_coeff(_opts->prob->source), b_vcoeff(_opts->prob->nbc),
		  opts(_opts), nl_coeff(_opts->prob->nlc[l])
{
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
		  
		  U = new HypreParVector(*A);
        V = new HypreParVector(*A);
        W = new HypreParVector(*A);
        Z = new HypreParVector(*A);
		  *V = 0.0;
        newt = 0; v = 0;
}

void NonlinearSystemODE::AssembleDiffusionMatrix() const
{
	delete A;
   a->Update();
	a->Assemble();
	a->Finalize();
	A = a->ParallelAssemble();
}

void NonlinearSystemODE::AssembleBVector() const
{
	b_coeff->SetTime(GetTime());
	b_vcoeff->SetTime(GetTime());
	b_form->Assemble();
	if (!b)
		b = b_form->ParallelAssemble();
	else
		b_form->ParallelAssemble(*b);
}

void NonlinearSystemODE::Mult(const Vector &x, Vector &y) const
{
        X->SetData(x.GetData());
        Y->SetData(y.GetData());
        
        nl_coeff->SetGridFunction(*X, 1.0, 0.0);
        AssembleBVector() ;
        AssembleDiffusionMatrix();

        A->Mult(*X, *Z); // M^{-1} A X
        if (b)
            *Z += *b; // M^{-1} (A X + b)
        pcg->Mult(*Z, *Y);
}

//Iterate over the system M*k^(j+1) = A(x + dt*k(j)) + b(x,t)
void NonlinearSystemODE::ImplicitSolve(const double dt, const Vector &x, Vector &k)
{

        //Set the nonlinear options
        int iter = 1;
        double resid = 1;
        double tol;
        int max_it;
        
        if (opts->braid_level == 0)
        {
            tol = opts->tol_fine;
            max_it = opts->mi_fine;
        }
        else
        {
            tol = opts->tol_coarse;
            max_it = opts->mi_coarse;   
        }  

        //Wrap the data (X == x, Y == k)
        X->SetData(x.GetData());
        Y->SetData(k.GetData());
        
        //Assemble the source term vector 
        AssembleBVector();

        //Set Inital guess and residual
       *Y = *opts->u; //Set Y to be the initial guess -- used by hypre
       nl_coeff->SetGridFunction(*Y, 1.0,0.0);
       AssembleDiffusionMatrix(); //A = A[Y] 
       
       if (!B)
         B = new HypreParMatrix(hypre_ParCSRMatrixAdd(*M, *A));
       
       hypre_ParCSRMatrixSetConstantValues(*B, 0.0);
       hypre_ParCSRMatrixSum(*B, 1.0, *M);
       hypre_ParCSRMatrixSum(*B, -dt, *A);
       
       
       B_amg = new HypreBoomerAMG(*B);
       B_amg->SetPrintLevel(-1); 
       B_pcg = new HyprePCG(*B);
       B_pcg->SetTol(tol/100);
       B_pcg->SetMaxIter(2000);
       B_pcg->SetPrintLevel(0);
       B_pcg->SetPreconditioner(*B_amg); 

       //Form the rhs: A x + bi
       A->Mult(*X, *Z);
       *Z += *b;
       //Solve the system B Y = Z 
       B_pcg->Mult(*Z, *Y);

       UpdateResidual(dt, resid); 
       double resid_start = resid;
       int temp;
       B_pcg->GetNumIterations(temp);
       delete B_pcg;
       delete B_amg;
       //Newton Iteration -- Makes sure we always reduce the residual
		 while ((iter < max_it && resid > tol) || resid >= resid_start )
		 {  
            iter+=1;
			   nl_coeff->SetGridFunction(*Z, 1.0,1.0);
			   AssembleDiffusionMatrix();  //A = Newton Linearization Matrix
            
            hypre_ParCSRMatrixSetConstantValues(*B, 0.0);
            hypre_ParCSRMatrixSum(*B, 1.0, *M);
            hypre_ParCSRMatrixSum(*B, -dt, *A);

            B_amg = new HypreBoomerAMG(*B);
            B_amg->SetPrintLevel(-1); 
            B_pcg = new HyprePCG(*B);
            B_pcg->SetTol(tol);
            B_pcg->SetMaxIter(2000);
            B_pcg->SetPrintLevel(0);
            B_pcg->SetPreconditioner(*B_amg); 

            *U = 0.0; // Set initial condition for the linear solver 
            B_pcg->Mult(*W, *U); //W is the RHS calculated in Update residual
            B_pcg->GetNumIterations(temp);
            v+= temp;
            *Y -= *U;         
            delete B_pcg;
            delete B_amg;
            UpdateResidual(dt, resid); //Updates W=RHS and Z = u_k+1 to match y==k         
       }  

       braid_iter = opts->braid_iter;
       newt += iter;
        if (opts->nl_print_level)
           cout << "NonLinear Solver:" << dt << " , " << setw(10) << t << " , "  << setw(10)
           << newt << " , " << setw(10) << braid_iter << " , " 
           << setw(10) << v << endl;
}

void NonlinearSystemODE::UpdateResidual(double dt, double &resid)
{
      add(*X, dt, *Y, *Z);  
      nl_coeff->SetGridFunction(*Z, 1.0,0.0); 
      AssembleDiffusionMatrix(); //A = A[x + k(j)*dt] 
      A->Mult(*Z,*U);			// U = A[x+k(j)*dt](x+k(j)*dt) 
      M->Mult(*Y,*W);			 //b = Mk(j)
      *W -= *U;
      *W -= *b;  // = Mk(j) - A(x+dt*k(j)) - Source vector
      resid = opts->dx*sqrt(InnerProduct(W,W));
    
}

// Solve the linear system Mk = A[x + k*dt] + b, where A is linear, for k
void NonlinearSystemODE::SetHypre(double dt, double tol) 
{
   //Set up BOOMERamg for inverting B
   B_amg = new HypreBoomerAMG(*B);
   B_amg->SetPrintLevel(-1); 
   B_pcg = new HyprePCG(*B);
   B_pcg->SetTol(tol);
   B_pcg->SetMaxIter(2000);
   B_pcg->SetPrintLevel(0);
   B_pcg->SetPreconditioner(*B_amg); 
  
}

NonlinearSystemODE::~NonlinearSystemODE() 
{  
   delete B;
   delete b;
   delete A;
   delete pcg;
   delete amg;
   delete b_form;
   delete a;
   delete M;

   delete X;
   delete Y;
   delete Z;
   delete U;
   delete V;
   delete W;
   
}



/*************************************************************************************
Implimentation::NonlinearApp */
NonlinearApp::NonlinearApp(
			NonlinearOptions *_opts, MPI_Comm comm_t_, ParMesh *pmesh)
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

int NonlinearApp::Step(braid_Vector    u_,
                   braid_Vector    ustop_,
                   braid_Vector    fstop_,
                   BraidStepStatus &pstatus)
{
   // This contains one small change over the default Phi, we store the Space-Time mesh info
   
   // Extract info from braid_vector
   BraidVector *u = (BraidVector*) u_;
   BraidVector *ustop = (BraidVector*) ustop_;
 
   int braid_level; 
   int braid_iter;
   pstatus.GetLevel(&braid_level);
   pstatus.GetIter(&braid_iter);

   opts->braid_iter = braid_iter;
   opts->braid_level =  braid_level;
   
   double tstart, tstop, dt;
   pstatus.GetTstartTstop(&tstart, &tstop);
   dt = tstop - tstart;
   
   //Store run-time mesh info
   if (opts->sc_print_level)
        MeshInfo.SetRow(braid_level, u->level, x[u->level]->ParFESpace()->GetParMesh(), dt);

   opts->u = ustop; 
   MFEMBraidApp::Step(u_,ustop_, fstop_, pstatus);       

   // no refinement
   pstatus.SetRFactor(1);
   return 0;
}

void NonlinearApp::AllocLevels(int num_levels){}

ParFiniteElementSpace* NonlinearApp::ConstructFESpace(ParMesh *pmesh)
{
	 return new ParFiniteElementSpace(pmesh, fe_coll);
}

void NonlinearApp::InitLevel(int l)
{ 
    //Create Bilinear and linear forms
    ParBilinearForm *a = new ParBilinearForm(fe_space[l]);
    ParBilinearForm *m = new ParBilinearForm(fe_space[l]);
    ParLinearForm   *b = new ParLinearForm(fe_space[l]);
    
    //Create a nonlinear coefficient for this level  
    opts->prob->nlc[l] = new NonlinearCoefficient(x[l], opts->power, 2);

    //Add Boundary and domain integrators 
    m->AddDomainIntegrator(new MassIntegrator(opts->prob->one));
    b->AddDomainIntegrator(new DomainLFIntegrator(*opts->prob->source));
    b->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(*opts->prob->nbc));
    a->AddDomainIntegrator(new DiffusionIntegrator(*opts->prob->nlc[l]));
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
    ode[l] = new NonlinearSystemODE(l, a, M, b, opts );

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

NonlinearApp:: ~NonlinearApp()
{
   delete fe_coll;
}			  
