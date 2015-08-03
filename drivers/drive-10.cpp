#include <fstream>
#include <iostream>
#include "braid_mfem_block.hpp"

using namespace std;


/* -------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------- */


void InitialDeformation(const Vector &x, Vector &y);


void InitialVelocity(const Vector &x, Vector &v);


class BackwardEulerOperator;


/** After spatial discretization, the hyperelastic model can be written as a
 *  system of ODEs:
 *	  dv/dt = -M^{-1}*(H(x) + S*v)
 *	  dx/dt = v,
 *  where x is the vector representing the deformation, v is the velocity field,
 *  M is the mass matrix, S is the viscosity matrix, and H(x) is the nonlinear
 *  hyperelastic operator.
 *
 *  Class HyperelasticOperator represents the right-hand side of the above
 *  system of ODEs. */
class HyperelasticOperator : public TimeDependentOperator
{
protected:
	ParFiniteElementSpace &fespace;

	ParBilinearForm M, S;
	ParNonlinearForm H;
	double viscosity;
	HyperelasticModel *model;

	HypreParMatrix *Mmat; // Mass matrix from ParallelAssemble()
	CGSolver M_solver;	 // Krylov solver for inverting the mass matrix M
	HypreSmoother M_prec; // Preconditioner for the mass matrix M

	/** Nonlinear operator defining the reduced backward Euler equation for the
		 velocity. Used in the implementation of method ImplicitSolve. */
	BackwardEulerOperator *backward_euler_oper;
	// Newton solver for the backward Euler equation
	NewtonSolver newton_solver;
	// Solver for the Jacobian solve in the Newton method
	Solver *J_solver;
	// Preconditioner for the Jacobian
	Solver *J_prec;
	// auxiliary vector
	mutable Vector z; 
	// initial guess for Newton Solver
	BraidVector *guess;

public:

	HyperelasticOperator(ParFiniteElementSpace &f, Array<int> &ess_bdr,
						 double visc);

	virtual void Mult(const Vector &vx, Vector &dvx_dt) const;
	/** Solve the Backward-Euler equation: k = f(x + dt*k, t), for the unknown k.
		 This is the only requirement for high-order SDIRK implicit integration.*/
	virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k);

	double ElasticEnergy(ParGridFunction &x) const;
	double KineticEnergy(ParGridFunction &v) const;
	void GetElasticEnergyDensity(ParGridFunction &x, ParGridFunction &w) const;
	void SetInitialGuess(BraidVector *guess_);

	virtual ~HyperelasticOperator();
};


// Nonlinear operator of the form:
//		 k --> (M + dt*S)*k + H(x + dt*v + dt^2*k) + S*v,
// where M and S are given BilinearForms, H is a given NonlinearForm, v and x
// are given vectors, and dt is a scalar.
class BackwardEulerOperator : public Operator
{
private:
	ParBilinearForm *M, *S;
	ParNonlinearForm *H;
	mutable HypreParMatrix *Jacobian;
	const Vector *v, *x;
	double dt;
	mutable Vector w, z;

public:
	BackwardEulerOperator(ParBilinearForm *M_, ParBilinearForm *S_,
						  ParNonlinearForm *H_);
	void SetParameters(double dt_, const Vector *v_, const Vector *x_);
	virtual void Mult(const Vector &k, Vector &y) const;
	virtual Operator &GetGradient(const Vector &k) const;
	virtual ~BackwardEulerOperator();
};


/** Function representing the elastic energy density for the given hyperelastic
 *  model+deformation. Used in HyperelasticOperator::GetElasticEnergyDensity. */
class ElasticEnergyCoefficient : public Coefficient
{
private:
	HyperelasticModel &model;
	ParGridFunction	&x;
	DenseMatrix J;

public:
	ElasticEnergyCoefficient(HyperelasticModel &m, ParGridFunction &x_)
		: model(m), x(x_) { }
	virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
	virtual ~ElasticEnergyCoefficient() { }
};


/** Structure including options parsed from command line. Inherited from
 *  OptionsParser -> BraidOptions -> BeamOptions */
struct BeamOptions : public BraidOptions 
{
	int	order,
		ode_solver_type,
		visport,
		vis_time_steps,
		vis_braid_steps;

	bool vis_screenshots,
		 visualization;

	double dt, // derived from t_start, t_final, and num_time_steps
		   visc;

	const char *vishost;

	BeamOptions(int argc, char *argv[]);

};


/** App to be passed to XBraid containing all necessary variables and functions.
 *  Inherited from BraidApp -> MFEMBraidApp -> BeamApp */
class BeamApp : public MFEMBraidApp
{
protected:

	BeamOptions &opts; 
	H1_FECollection fe_coll;

	virtual void AllocLevels(int num_space_levels);
	virtual ParFiniteElementSpace *ConstructFESpace(ParMesh *pmesh);
	virtual void InitLevel(int l);


public:

	BeamApp(BeamOptions &options, MPI_Comm comm_t, ParMesh *pmesh);

	virtual int Step(braid_Vector	 u_,
					 braid_Vector	 ustop_,
					 braid_Vector	 fstop_,
					 BraidStepStatus &pstatus);

	virtual ~BeamApp() { };

};


/* -------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------- */


int main(int argc, char *argv[])
{
	// Initialize MPI.
	int myid, num_procs;
	MPI_Init(&argc, &argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	// Parse command line for Braid, MFEM and this problem. 
	BeamOptions opts(argc,argv);
	if (!opts.Good()) // check valid input, from MFEM OptionsParser class. 
	{
		if (myid == 0)
		{
			opts.PrintUsage(cout);
		}
		MPI_Finalize();
		return 1;
	}
	if (myid == 0) 
	{
		opts.PrintOptions(cout);
	}

	// Read serial mesh from given mesh file on all processors;
	// refine uniformly in serial to increase resolution.  
	Mesh *mesh = opts.LoadMeshAndSerialRefine();
	if (!mesh)
	{
	  if (myid == 0)
	  {
		  cerr << "\nError loading mesh file: " << opts.mesh_file << '\n' << endl;
	  }
	  MPI_Finalize();
	  return 2;
	}

	// Split global MPI communicator into spatial and temporal communicators.
	// Define parallel mesh by partitioning serial mesh. Parallel refinement
	// is done in MFEMBraidApp::InitMultilevelApp(). Once parallel mesh is
	// formed we can delete serial mesh.
	BraidUtil util;
	MPI_Comm comm_x, comm_t;
	util.SplitCommworld(&comm, opts.num_procs_x, &comm_x, &comm_t);
	ParMesh *pmesh = new ParMesh(comm_x, *mesh);
	delete mesh;

	// Create and initialize BeamApp. 
	BeamApp app(opts, comm_t, pmesh);

	// Run Braid simulation.
	BraidCore core(comm, &app);
	opts.SetBraidCoreOptions(core);
	core.Drive();

	MPI_Finalize();
	return 0;
}


/* -------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------- */


BeamApp::BeamApp(BeamOptions &options, MPI_Comm comm_t, ParMesh *pmesh) 
	
	: MFEMBraidApp(comm_t, options.t_start, options.t_final, options.num_time_steps),
	  fe_coll(options.order, pmesh->Dimension()),
	  opts(options)
{
	// Initialize multilevel structures.
	InitMultilevelApp(pmesh, opts.par_ref_levels, opts.spatial_coarsen);

	// Define parallel vector finite element spaces representing mesh deformation
	// x_gf, velocity v_gf, and initial configuration x_ref. Define also the elastic
	// energy density, w_gf, which is in a discontinuous higher-order space. Since
	// x and v are integrated in time as a system, we group them together in block
	// vector vx, on the unique parallel degrees of freedom, with offsets given by
	// array true_offset.

	int glob_size = (fe_space[0])->GlobalTrueVSize(),
		true_size = (fe_space[0])->TrueVSize();
	// if (myid == 0)
	// 	cout << "Number of velocity/deformation unknowns: " << glob_size << endl;

	Array<int> true_offset(3);
	true_offset[0] = 0;
	true_offset[1] = true_size;
	true_offset[2] = 2*true_size;
	X0 = new BlockVector(true_offset);
	// Lose ownership of data so not deleted when true_offset is out of scope
	true_offset.LoseData();

	ParGridFunction v_gf(fe_space[0]),
					x_gf(fe_space[0]),
					x_ref(fe_space[0]);
	pmesh->GetNodes(x_ref);

	L2_FECollection w_fec(opts.order + 1, pmesh->Dimension());
	ParFiniteElementSpace w_fespace(pmesh, &w_fec);
	ParGridFunction w_gf(&w_fespace);

	// Set initial conditions for v_gf, x_gf, vx, and define boundary conditions
	// on a beam-like mesh.
	VectorFunctionCoefficient velo(pmesh->Dimension(), InitialVelocity);
	v_gf.ProjectCoefficient(velo);
	VectorFunctionCoefficient deform(pmesh->Dimension(), InitialDeformation);
	x_gf.ProjectCoefficient(deform);

	v_gf.GetTrueDofs(X0->GetBlock(0));
	x_gf.GetTrueDofs(X0->GetBlock(1));
}


void BeamApp::AllocLevels(int num_space_levels)
{

}


ParFiniteElementSpace* BeamApp::ConstructFESpace(ParMesh *pmesh)
{
	return new ParFiniteElementSpace(pmesh, &fe_coll, pmesh->Dimension());
}


void BeamApp::InitLevel(int l)
{
	Array<int> ess_bdr(fe_space[l]->GetMesh()->bdr_attributes.Max());
	ess_bdr = 0;
	ess_bdr[0] = 1; // boundary attribute 1 (index 0) is fixed

	// Initialize hyperelastic operator for level l.
	ode[l] = new HyperelasticOperator(*fe_space[l], ess_bdr, opts.visc);

	// Define and initialize ODE solver used for time integration on level l.
	ODESolver *ode_solver = NULL;
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
		default: ode_solver = new BackwardEulerSolver; break;
	}
	solver[l] = ode_solver;
	solver[l]->Init(*ode[l]);

	// Set maximum dt to keep ratio between dt, dx appropriate. 
	// For this problem we want dt ~ dx, where we coarsen dx by 2 each level. 
	max_dt[l] = 1.01 * opts.dt * pow(2,l);
}


int BeamApp::Step(braid_Vector	 u_,
				  braid_Vector	 ustop_,
				  braid_Vector	 fstop_,
				  BraidStepStatus &pstatus)
{
	BraidVector *u	  = (BraidVector*) u_;
	BraidVector *ustop = (BraidVector*) ustop_;
	BraidVector *fstop = (BraidVector*) fstop_;
	int spatial_level = u->spatial_level;
	double tstart, tstop, accuracy, t, dt;
	int braid_level;

	// Get time step information
	pstatus.GetTstartTstop(&tstart, &tstop);
	// pstatus.GetAccuracy(&accuracy);
	pstatus.GetLevel(&braid_level);

	t = tstart;
	dt = tstop - tstart;

	// Set initial guess for Newton solver and solve. 
	HyperelasticOperator *he_op;
	he_op = dynamic_cast<HyperelasticOperator *>(ode[spatial_level]);
	he_op->SetInitialGuess(ustop);
	solver[spatial_level]->Step(*u, t, dt); // Note, this is MFEM Step function
	he_op->SetInitialGuess(NULL);

	// no refinement
	pstatus.SetRFactor(1);

	return 0;
}


BeamOptions::BeamOptions(int argc, char *argv[])
	: BraidOptions(argc, argv)
{
	// Set default values for time.
	t_start		  = 0.0;
	t_final		  = 300.0;
	num_time_steps = 100;

	// Set default values for mesh.
	mesh_file = "../data/beam-quad.mesh";
	ser_ref_levels = 2;
	par_ref_levels = 0;

	// Set default values for beam problem. 
	order 			= 2;
	ode_solver_type = 1;
	visc			= 1e-2;
	
	AddOption(&order, "-o", "--order",
			  "Order (degree) of the finite elements.");
	AddOption(&ode_solver_type, "-s", "--ode-solver",
			  "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3,\n\t"
			  "\t	11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.");
	AddOption(&visc, "-v", "--viscosity",
			  "Viscosity coefficient.");
	Parse();

	dt = (t_final - t_start) / num_time_steps;
}


/* -------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------- */


BackwardEulerOperator::BackwardEulerOperator(
	ParBilinearForm *M_, ParBilinearForm *S_, ParNonlinearForm *H_)
	: Operator(M_->ParFESpace()->TrueVSize()), M(M_), S(S_), H(H_),
	  Jacobian(NULL), v(NULL), x(NULL), dt(0.0), w(height), z(height)
{ 

}


void BackwardEulerOperator::SetParameters(double dt_, const Vector *v_,
										  const Vector *x_)
{
	dt = dt_;  v = v_;  x = x_;
}


void BackwardEulerOperator::Mult(const Vector &k, Vector &y) const
{
	// compute: y = H(x + dt*(v + dt*k)) + M*k + S*(v + dt*k)
	// w = (*v) + dt*k and z = (*x) + dt*w
	add(*v, dt, k, w);
	add(*x, dt, w, z);
	H->Mult(z, y);
	M->TrueAddMult(k, y);
	S->TrueAddMult(w, y);
}


Operator &BackwardEulerOperator::GetGradient(const Vector &k) const
{
	delete Jacobian;
	SparseMatrix *localJ = Add(1.0, M->SpMat(), dt, S->SpMat());
	add(*v, dt, k, w);
	add(*x, dt, w, z);
	localJ->Add(dt*dt, H->GetLocalGradient(z));
	Jacobian = M->ParallelAssemble(localJ);
	delete localJ;
	return *Jacobian;
}


BackwardEulerOperator::~BackwardEulerOperator()
{
	delete Jacobian;
}


/* -------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------- */


HyperelasticOperator::HyperelasticOperator(ParFiniteElementSpace &f,
										   Array<int> &ess_bdr, double visc)
	: TimeDependentOperator(2*f.TrueVSize(), 0.0), fespace(f),
	  M(&fespace), S(&fespace), H(&fespace), M_solver(f.GetComm()),
	  newton_solver(f.GetComm()), z(height/2)
{
	const double rel_tol = 1e-10;
	const int skip_zero_entries = 0;

	const double ref_density = 1.0; // density in the reference configuration
	ConstantCoefficient rho0(ref_density);
	M.AddDomainIntegrator(new VectorMassIntegrator(rho0));
	M.Assemble(skip_zero_entries);
	M.EliminateEssentialBC(ess_bdr);
	M.Finalize(skip_zero_entries);
	Mmat = M.ParallelAssemble();

	M_solver.iterative_mode = false;
	M_solver.SetRelTol(rel_tol);
	M_solver.SetAbsTol(0.0);
	M_solver.SetMaxIter(30);
	M_solver.SetPrintLevel(0);
	M_prec.SetType(HypreSmoother::Jacobi);
	M_solver.SetPreconditioner(M_prec);
	M_solver.SetOperator(*Mmat);

	double mu = 0.25; // shear modulus
	double K  = 5.0;  // bulk modulus
	model = new NeoHookeanModel(mu, K);
	H.AddDomainIntegrator(new HyperelasticNLFIntegrator(model));
	H.SetEssentialBC(ess_bdr);

	viscosity = visc;
	ConstantCoefficient visc_coeff(viscosity);
	S.AddDomainIntegrator(new VectorDiffusionIntegrator(visc_coeff));
	S.Assemble(skip_zero_entries);
	S.EliminateEssentialBC(ess_bdr);
	S.Finalize(skip_zero_entries);

	backward_euler_oper = new BackwardEulerOperator(&M, &S, &H);

	HypreSmoother *J_hypreSmoother = new HypreSmoother;
	J_hypreSmoother->SetType(HypreSmoother::l1Jacobi);
	J_prec = J_hypreSmoother;

	MINRESSolver *J_minres = new MINRESSolver(f.GetComm());
	J_minres->SetRelTol(rel_tol);
	J_minres->SetAbsTol(0.0);
	J_minres->SetMaxIter(300);
	J_minres->SetPrintLevel(-1);
	J_minres->SetPreconditioner(*J_prec);
	J_solver = J_minres;

	newton_solver.iterative_mode = true; // provide inital guess for Newton solves
	newton_solver.SetSolver(*J_solver);
	newton_solver.SetOperator(*backward_euler_oper);
	newton_solver.SetPrintLevel(-1); // don't print Newton iterations
	newton_solver.SetRelTol(rel_tol);
	newton_solver.SetAbsTol(1e-12);
	newton_solver.SetMaxIter(20);
}


void HyperelasticOperator::Mult(const Vector &vx, Vector &dvx_dt) const
{
	// Create views to the sub-vectors v, x of vx, and dv_dt, dx_dt of dvx_dt
	int sc = height/2;
	Vector v(vx.GetData() +  0, sc);
	Vector x(vx.GetData() + sc, sc);
	Vector dv_dt(dvx_dt.GetData() +  0, sc);
	Vector dx_dt(dvx_dt.GetData() + sc, sc);

	H.Mult(x, z);
	if (viscosity != 0.0)
		S.TrueAddMult(v, z);
	z.Neg(); // z = -z
	M_solver.Mult(z, dv_dt);

	dx_dt = v;
}


void HyperelasticOperator::ImplicitSolve(const double dt,
										 const Vector &vx, Vector &dvx_dt)
{
	int sc = height/2;
	Vector v(vx.GetData() +  0, sc);
	Vector x(vx.GetData() + sc, sc);
	Vector dv_dt(dvx_dt.GetData() +  0, sc);
	Vector dx_dt(dvx_dt.GetData() + sc, sc);

	// Set initial guess using slope between previous time step at current iteration and
	// current time step at previous iteration --> dv_dt = (y_{n+1,old} - y_{n,new}) / dt. 
	Vector g( guess->GetBlock(0) );
	dv_dt  = g;
	dv_dt -= v;
	dv_dt /= dt;

	// By eliminating kx from the coupled system:
	//	 kv = -M^{-1}*[H(x + dt*kx) + S*(v + dt*kv)]
	//	 kx = v + dt*kv
	// we reduce it to a nonlinear equation for kv, represented by the
	// backward_euler_oper. This equation is solved with the newton_solver
	// object (using J_solver and J_prec internally).
	backward_euler_oper->SetParameters(dt, &v, &x);
	Vector zero; // empty vector is interpreted as zero r.h.s. by NewtonSolver
	newton_solver.Mult(zero, dv_dt);
	add(v, dt, dv_dt, dx_dt);

	// MFEM_VERIFY(newton_solver.GetConverged(), "Newton Solver did not converge.");
}


double HyperelasticOperator::ElasticEnergy(ParGridFunction &x) const
{
	return H.GetEnergy(x);
}


double HyperelasticOperator::KineticEnergy(ParGridFunction &v) const
{
	double loc_energy = 0.5*M.InnerProduct(v, v);
	double energy;
	MPI_Allreduce(&loc_energy, &energy, 1, MPI_DOUBLE, MPI_SUM,
				   fespace.GetComm());
	return energy;
}


void HyperelasticOperator::GetElasticEnergyDensity(
	ParGridFunction &x, ParGridFunction &w) const
{
	ElasticEnergyCoefficient w_coeff(*model, x);
	w.ProjectCoefficient(w_coeff);
}


void HyperelasticOperator::SetInitialGuess(BraidVector *guess_)
{
	guess = guess_;
}


HyperelasticOperator::~HyperelasticOperator()
{
	delete model;
	delete backward_euler_oper;
	delete J_solver;
	delete J_prec;
	delete Mmat;
}


/* -------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------- */


double ElasticEnergyCoefficient::Eval(ElementTransformation &T,
									  const IntegrationPoint &ip)
{
	model.SetTransformation(T);
	x.GetVectorGradient(T, J);
	// return model.EvalW(J);  // in reference configuration
	return model.EvalW(J)/J.Det(); // in deformed configuration
}


/* -------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------- */


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