#include "mfem.hpp"
#include "braid_mfem_block.hpp"
//Extra Hypre Functions
#include "hypre_extra.hpp"

#include <fstream>
#include <iostream>

using namespace hypre;
using namespace std;
using namespace mfem;



/* Load vector from file. */
void LoadVector(Vector &data, double t, int type)
{
	ostringstream filestream;
	if (type == 1) {
		filestream << "./sol/Vector_" << t << ".txt";
	}
	else {
		filestream << "./sol/Mesh_" << t << ".txt";
	}

	string filename = filestream.str();
    ifstream input_file;
    input_file.open(filename.c_str());

    int i=0;
    if (input_file) {
    	while(!input_file.eof()) {
    		input_file >> data(i);
    		i++;
    	}
	}

	if (type == 1) {
		cout << "Vector at t = " << t << endl;
	}
	else {
		cout << "Mesh at t = " << t << endl;
	}
	data.Print();
}


/* Initial conditions of vector. */
double InitialConditions(Vector &x0);


/* Gaussian bump in (x,t) with given center and width in each dimensions. */
double GuassianBlip(const double &x, const double &c, const double &w, const double &scale);


/* Forcing function f(x,t). */
void Forcing(const Vector &x, double t, Vector &f);


class DiffusionOperator : public TimeDependentOperator
{
private:

	ConstantCoefficient			*alpha_;
	VectorFunctionCoefficient 	*f_;
	ParLinearForm 				*b_form_;
	ParBilinearForm 			*a_form_,
								*m_form_,
								*mesh_form_;
	HypreParMatrix 				*A_,
								*M_,
								*B_,
								*Amesh_;
	HypreParVector 				*z_,
								*k_,
								*x_;
	mutable HypreParVector 		*b_;
	HypreBoomerAMG 				*B_amg_,
								*M_amg_;
	HyprePCG 					*B_pcg_,
								*M_pcg_;
	Vector 						mesh_monitor_,
								mesh_displace_, // displace = new_nodes - old_nodes.
								mesh_nodes_;
	double 						current_dt_;
	const double				dzeta_,
								mesh_tau_;
	int 						numNodes_;
	Array<int> 					ess_bdr_;

	void UpdateSpatialSolvers(const double & dt = -1.0);
	void Assemble_b_Vector(const double &t = -1) const;
	// Arclength monitor function.
	void UpdateMeshMonitor(const Vector &u_curr);
	void UpdateMeshSolver();
	void InterpolateSolution(Vector &u_curr);

public:

	DiffusionOperator(ParMesh *pmesh, ParFiniteElementSpace *fespace, BlockVector *&X0,
					  const double &dzeta, const double &alpha = 1.0, const double &tau = 1.0);
	void GetResidual(Vector &u_stop, Vector &u_start, Vector &temp, const double &dt,
					 const double &t, const int &update_ops);
	void UpdateMesh(ParMesh *pmesh, Vector &u_curr, int &update_ops, int &move_mesh, const double &dt);
	virtual void Mult(const Vector &x, Vector &y) const;
	virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k);
	virtual ~DiffusionOperator();
};



struct MovingOptions: public BraidOptions
{

   int    num_intervals,
   		  order,
	      ode_solver_type,
	      x_dim,
	      forcing_eqn,
	      visport,
	   	  vis_time_steps,
	   	  vis_braid_steps;

   bool	  vis_screenshots,
		  visualization;

   double dt, // derived from t_start, t_final, and num_time_steps
   		  dzeta,
   		  max_x,
   		  alpha, 
   		  tau;

   const char *vishost;

   MovingOptions(int argc, char *argv[]);

};


class MovingApp : public MFEMBraidApp
{
protected:

	MovingOptions &opts;
	H1_FECollection fe_coll; 

   // Allocate data structures for the given number of spatial levels. Used by
   // InitMultilevelApp.
   virtual void AllocLevels(int num_space_levels);

   // Construct ParFiniteElmentSpace the given mesh. Used by InitMultilevelApp.
   virtual ParFiniteElementSpace *ConstructFESpace(ParMesh *pmesh);

   // Assuming mesh[l] and fe_space[l] are set, initializes ode[l], solver[l],
   // and max_dt[l]. Used by InitMultilevelApp.
   virtual void InitLevel(int l);


public:

	MovingApp(MovingOptions &options, MPI_Comm comm_t, ParMesh *pmesh);

	void InterpolateSolution(Vector &u_curr, Vector &old_node, Vector &displace);
	virtual int Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_,
					 BraidStepStatus &pstatus);
	virtual int Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus &pstatus);
	virtual int Init(double t, braid_Vector *u_ptr);
	virtual int Sum(double alpha, braid_Vector a_, double beta, braid_Vector b_);
	virtual int SpatialNorm(braid_Vector u_, double *norm_ptr);
	virtual ~MovingApp();

};


/* -------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------- */

// Global variable on which forcing funciton to use. 
int forcing_eqn;


int main(int argc, char *argv[]) 
{

	// Initialize MPI.
	int myid, num_procs;
	MPI_Init(&argc, &argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &myid);
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	// Parse command line for Braid, MFEM and this problem. 
	MovingOptions opts(argc, argv);
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
	forcing_eqn = opts.forcing_eqn;

	// Block scope so that objects are destroyed before MPI_Finalize();
	{
		// Construct serial mesh on all processors. 
		Mesh mesh(opts.num_intervals,opts.max_x);

		// Split global MPI communicator into spatial and temporal communicators.
		// Define parallel mesh by partitioning serial mesh. Parallel refinement
		// is done in MFEMBraidApp::InitMultilevelApp(). Once parallel mesh is
		// formed we can delete serial mesh.
		BraidUtil util;
		MPI_Comm comm_x, comm_t;
		util.SplitCommworld(&comm, opts.num_procs_x, &comm_x, &comm_t);
		ParMesh *pmesh = new ParMesh(comm_x,mesh);

		// Create and initialize MovingApp. 
		MovingApp app(opts, comm_t, pmesh);

		// Run Braid simulation.
		BraidCore core(comm, &app);
		opts.SetBraidCoreOptions(core);
		core.Drive();

                // MPI_Comm_free( &comm );
                MPI_Comm_free( &comm_x );
                MPI_Comm_free( &comm_t );
	}

	MPI_Finalize();
	return 0;
}


/* -------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------- */


MovingOptions::MovingOptions(int argc, char *argv[])
	: BraidOptions(argc, argv)
{
	// Set default values for time.
	t_start        = 0.0;
	t_final        = 10.0;
	num_time_steps = 50;

	// Set default values for mesh.
	ser_ref_levels = 0;
	par_ref_levels = 0;
	forcing_eqn = 1;
	// AddMeshOptions();

	// Set default values for diffusion problem. 
	num_intervals	= 10;
	order 			= 1;
	ode_solver_type = 1;
	alpha			= 1.0;
	x_dim 			= 1;
	max_x			= 1.0;
	tau 			= 1.0;

	// Set default values for visualization.
	vis_time_steps 	= 10;
	vis_braid_steps = 1;
	visualization 	= false;
	vis_screenshots = false;
	visport 		= 19916;
	vishost 		= "localhost";

	AddOption(&num_intervals, "-n","--nintervals",
			  "Numer of spatial intervals.");
	AddOption(&alpha, "-alpha", "--diffusivity",
	          "Thermal diffusivity.");
	AddOption(&tau, "-tau", "--tau",
	          "Moving mesh speed parameter.");
	AddOption(&forcing_eqn, "-forcing", "--forcing-equation",
			  "Choice of forcing funciton; 0 = none, 1 = moving bump, 2 = five fixed bumps.");
	AddOption(&order, "-o", "--order",
	          "Order (degree) of the finite elements.");
	AddOption(&ode_solver_type, "-s", "--ode-solver",
	          "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3,\n\t"
	          "\t   11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.");
	AddOption(&visualization, "-vis", "--visualization", "-no-vis", "--no-visualization",
	          "Enable or disable GLVis visualization.");
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

	dt    = (t_final - t_start) / num_time_steps;
	dzeta = max_x / num_intervals;

	// scale stopping tolerance by grid element size.
	tol *= (num_intervals*num_time_steps);
}


/* -------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------- */



MovingApp::MovingApp(MovingOptions &options, MPI_Comm comm_t, ParMesh *pmesh) 
	
	: MFEMBraidApp(comm_t, options.t_start, options.t_final, options.num_time_steps),
	  fe_coll(options.order, pmesh->Dimension()),
	  opts(options)
{
	// Initialize multilevel structures.
	InitMultilevelApp(pmesh, opts.par_ref_levels, opts.spatial_coarsen);
}



void MovingApp::AllocLevels(int num_space_levels)
{

}


ParFiniteElementSpace* MovingApp::ConstructFESpace(ParMesh *pmesh)
{
	return new ParFiniteElementSpace(pmesh, &fe_coll, pmesh->Dimension());
}


void MovingApp::InitLevel(int level)
{
	// Initialize diffusion operator for level l.
	ode[level] = new DiffusionOperator(mesh[level], fe_space[level], this->X0, 
									   opts.dzeta, opts.alpha, opts.tau);

	// Define and initialize ODE solver used for time integration on given level.
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
	solver[level] = ode_solver;
	solver[level]->Init(*ode[level]);

	// Set maximum dt to keep ratio between dt, dx appropriate. 
	// For this problem we want dt ~ dx, where we coarsen dx by 2 each level. 
	max_dt[level] = 1.01 * opts.dt * pow(2,level);
}


void MovingApp::InterpolateSolution(Vector &u_curr, Vector &old_nodes, Vector &displace)
{
	int    lower,
	   	   upper,
	   	   size = u_curr.Size();
	double temp[size-2];
	double dx,
		   new_node;

	// Find j s.t. new_node(i) in [old_node(j),old_node(j+1)]
	for (int i=1; i<(size-1); i++) {
		if (abs(displace(i)) < 1e-6) {
			temp[i-1] = u_curr(i);
		}
		else {
			new_node = old_nodes(i) + displace(i);
			if (displace(i) < 0) {
				upper = i;
				while (old_nodes(upper) > new_node) {
					upper -= 1;
				}
				lower  = upper; 
				upper += 1;
			}
			else {
				lower = i;
				while (old_nodes(lower) < new_node) {
					lower += 1;
				}
				upper  = lower;
				lower -= 1;
			}
			// Interpolate solution linearly to new mesh. 
			dx		  = old_nodes(upper) - old_nodes(lower);
			temp[i-1] = ( (new_node - old_nodes(lower))*u_curr(upper) +
						  (old_nodes(upper) - new_node)*u_curr(lower) ) / dx;
		}
	}
	for (int i=1; i<(size-1); i++) {
		u_curr(i) = temp[i-1];
	}
}


int MovingApp::Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_,
                    BraidStepStatus &pstatus)
{
	BraidVector *u  = (BraidVector*) u_;
	Vector &u_grid	= u->GetBlock(0);
	Vector &u_value = u->GetBlock(1);

	int    spatial_level = u->spatial_level,
		   braid_level,
		   braid_iter,
		   update_ops = 0,
		   move_mesh = 1;
	double tstart,
		   tstop,
		   t,
		   dt;

	DiffusionOperator *diff_op;
	diff_op = dynamic_cast<DiffusionOperator *>(ode[spatial_level]);
	ParMesh *mesh_temp = this->mesh[spatial_level];

	// Get time step information.
	pstatus.GetTstartTstop(&tstart, &tstop);
	pstatus.GetLevel(&braid_level);
	pstatus.GetIter(&braid_iter);
	dt = tstop - tstart;
	t  = tstart;

	// Check if mesh in class object is different than mesh with current vector.
	Vector op_mesh_displace;
	mesh_temp->GetNodes(op_mesh_displace);
	op_mesh_displace -= u_grid;
	if ( abs(op_mesh_displace.Max()) > 1e-6 || 
		 abs(op_mesh_displace.Min()) > 1e-6 ) {	
		update_ops = 1;
		mesh_temp->SetNodes(u_grid);
	}


/* Save grid and vector at time t=0 to .csv file w.r.t. braid level and iteration. */
#if 0
if (t == 0) {
	ofstream output_file;
	ostringstream grid_stream;
	grid_stream << "./mesh_data/Grid_lev" << braid_level << "_iter" << braid_iter << ".csv";
	string grid_name = grid_stream.str();
	output_file.open(grid_name.c_str(), ios::app);
	output_file << setprecision(6) << t << ",";
	for (int i=0; i<u_grid.Size(); i++) {
		output_file << setprecision(14) << u_grid(i) << ",";
	}
	output_file << "\n";
	output_file.close();

	ostringstream vec_stream;
	vec_stream << "./mesh_data/Vec_lev" << braid_level << "_iter" << braid_iter << ".csv";
	string vec_name = vec_stream.str();
	output_file.open(vec_name.c_str(), ios::app);
	output_file << setprecision(6) << t << ",";
	for (int i=0; i<u_value.Size(); i++) {
		output_file << setprecision(14) << u_value(i) << ",";
	}
	output_file << "\n";
	output_file.close();
}
#endif


	// Move and upate mesh.
	diff_op->UpdateMesh(mesh_temp, u_value, update_ops, move_mesh, dt);
	if (move_mesh == 1) {
		mesh_temp->GetNodes(u_grid);
	}

	// Take time step.
	solver[spatial_level]->Step(u_value, t, dt);


/* Save grid and vector to .csv file w.r.t. braid level and iteration. */
#if 0
ofstream output_file;
ostringstream grid_stream;
grid_stream << "./mesh_data/Grid_lev" << braid_level << "_iter" << braid_iter << ".csv";
string grid_name = grid_stream.str();
output_file.open(grid_name.c_str(), ios::app);
output_file << setprecision(6) << t << ",";
for (int i=0; i<u_grid.Size(); i++) {
	output_file << setprecision(14) << u_grid(i) << ",";
}
output_file << "\n";
output_file.close();

ostringstream vec_stream;
vec_stream << "./mesh_data/Vec_lev" << braid_level << "_iter" << braid_iter << ".csv";
string vec_name = vec_stream.str();
output_file.open(vec_name.c_str(), ios::app);
output_file << setprecision(6) << t << ",";
for (int i=0; i<u_value.Size(); i++) {
	output_file << setprecision(14) << u_value(i) << ",";
}
output_file << "\n";
output_file.close();
#endif



	// no refinement
	pstatus.SetRFactor(1);

	mesh_temp    = NULL;
	diff_op 	 = NULL;
	u 			 = NULL;
	return 0;
}

/* Don't think this is correct. */
int MovingApp::Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus &pstatus)
{
	double tstart,
		   tstop,
		   dt;
	pstatus.GetTstartTstop(&tstart, &tstop);
	dt = tstop - tstart;

	BraidVector *u_stop  = (BraidVector*) u_;
	BraidVector *u_start = (BraidVector*) r_;
	Vector &grid_start   = u_start->GetBlock(0);
	Vector &grid_stop    = u_stop->GetBlock(0);
	Vector &val_start    = u_start->GetBlock(1);
	Vector &val_stop     = u_stop->GetBlock(1);
	
	Vector temp = grid_stop;
		   temp -= grid_start;

	// Interpolate u_start to grid associated with u_stop
	if ( abs(temp.Max()) > 1e-6 || 
		 abs(temp.Min()) > 1e-6 ) {	   
		InterpolateSolution(val_start, grid_start, temp);
		grid_start = grid_stop;
	}

	int spatial_level = u_stop->spatial_level;
	DiffusionOperator *diff_op;
	diff_op = dynamic_cast<DiffusionOperator *>(ode[spatial_level]);
	ParMesh *mesh_temp = this->mesh[spatial_level];

	// Update mesh in ParMesh to be same as u_stop
	int update_ops = 0;
	mesh_temp->GetNodes(temp);
	temp -= grid_stop;
	if ( abs(temp.Max()) > 1e-6 || 
		 abs(temp.Min()) > 1e-6 ) {	
		update_ops = 1;
		mesh_temp->SetNodes(grid_stop);
	}

	// Get residual
	diff_op->GetResidual(val_stop, val_start, temp, dt, tstop, update_ops);

	mesh_temp    = NULL;
	diff_op 	 = NULL;
	u_stop		 = NULL;
	u_start		 = NULL;
	return 0;
}


int MovingApp::Init(double        t,
                    braid_Vector *u_ptr)
{
	int spatial_level = 0;

	// Set initial mesh points equal for all time t.
	BraidVector *u = new BraidVector(spatial_level, *X0);

	// Set initial values equal to zero except on first block.
	if (t > 0) {
		u->GetBlock(1) = 0.;
	}


// Load C-point vectors from file
// Vector &values = u->GetBlock(1);
// Vector &mesh = u->GetBlock(0);
// LoadVector(values, t, 1);
// LoadVector(mesh, t, 0);



	*u_ptr = (braid_Vector) u;
	u = NULL;
	return 0;
}


int MovingApp::SpatialNorm(braid_Vector u_,
                           double      *norm_ptr)
{
	// double dot;
	BraidVector *u = (BraidVector*) u_;
	Vector &values = u->GetBlock(1);
	// dot = InnerProduct(u, u);
	// *norm_ptr = sqrt(dot);
	*norm_ptr = values.Norml2(); // Note this may be a serial implementation?
	return 0;
}


int MovingApp::Sum(double       alpha,
                   braid_Vector a_,
                   double       beta,
                   braid_Vector b_)
{
	BraidVector *a  = (BraidVector*) a_;
	BraidVector *b  = (BraidVector*) b_;
	Vector &a_grid  = a->GetBlock(0);
	Vector &b_grid  = b->GetBlock(0);
	Vector &a_value = a->GetBlock(1);
	Vector &b_value = b->GetBlock(1);

	// Add vectors. 
	add(alpha, a_value, beta, b_value, b_value);
	add(alpha, a_grid, beta, b_grid, b_grid);

	a = NULL;
	b = NULL;
	return 0;
}


MovingApp::~MovingApp()
{

}



/* -------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------- */



DiffusionOperator::DiffusionOperator(ParMesh *pmesh, ParFiniteElementSpace *fespace,
									 BlockVector *&X0, const double &dzeta,
									 const double &alpha, const double &tau)
	: dzeta_(dzeta), mesh_tau_(tau)
{
	// Assemble parallel linear form for right hand side, with vector Forcing
	// function, (f,phi_i) where phi_i are the basis functions in fespace.
	int dim = pmesh->Dimension();
	f_ 		= new VectorFunctionCoefficient(dim, &Forcing);
	b_form_ = new ParLinearForm(fespace);
	b_form_->AddDomainIntegrator(new VectorDomainLFIntegrator(*f_));
	b_form_->Assemble();

	// Assemble parallel bilinear forms and initial conditions.
	FunctionCoefficient u0(InitialConditions);
	ParGridFunction *u = new ParGridFunction(fespace);
	alpha_  = new ConstantCoefficient(alpha);
	a_form_ = new ParBilinearForm(fespace);
	m_form_ = new ParBilinearForm(fespace);
	a_form_->AddDomainIntegrator(new DiffusionIntegrator(*alpha_));
	m_form_->AddDomainIntegrator(new MassIntegrator);
	a_form_->Assemble(0);
	m_form_->Assemble(0);

	// Fix Dirichlet boundary conditions.
	ess_bdr_.SetSize(pmesh->bdr_attributes.Max());
	ess_bdr_ = 1;
	a_form_->EliminateEssentialBC(ess_bdr_, *u, *b_form_);
	m_form_->EliminateEssentialBC(ess_bdr_, *u, *b_form_);

	// Construct bilinear form for mesh equations. Add 1d diffusion to make
	// a tridiagonal operator.
	mesh_form_ = new ParBilinearForm(fespace);
	mesh_form_->AddDomainIntegrator(new DiffusionIntegrator(*alpha_));
	mesh_form_->Assemble();
	mesh_form_->Finalize(0);

	// Construct initial condition block vector, with mesh in first block
	// and solution values in second block. 
	int true_size = fespace->TrueVSize();
	Array<int> true_offset(3);
	true_offset[0] = 0;
	true_offset[1] = true_size;
	true_offset[2] = 2*true_size;
	X0 = new BlockVector(true_offset);
	true_offset.LoseData();

	u->ProjectCoefficient(u0);
	pmesh->GetNodes(X0->GetBlock(0));
	u->GetTrueDofs(X0->GetBlock(1));
	delete u;

	// Finalize bilinear forms and form parallel operators 
	a_form_->Finalize(0);
	m_form_->Finalize(0);
	A_ = a_form_->ParallelAssemble();
	M_ = m_form_->ParallelAssemble();

	// Define other class variables.
	width 		= A_->Width();
	height 		= A_->Height();
	numNodes_ 	= height;
	b_ 			= NULL;
	B_amg_ 		= NULL;
	B_pcg_ 		= NULL;
	M_amg_ 		= NULL;
	M_pcg_ 		= NULL;
	Amesh_ 		= NULL;
	current_dt_ = -1.0;
	mesh_monitor_.SetSize(numNodes_);
	mesh_displace_.SetSize(numNodes_);
	mesh_nodes_.SetSize(numNodes_);

	double tmp; // workaround to avoid memory leak
	x_ = new HypreParVector(A_->GetComm(), A_->GetGlobalNumCols(), &tmp,
						    A_->GetColStarts());
	k_ = new HypreParVector(A_->GetComm(), A_->GetGlobalNumCols(), &tmp,
						    A_->GetColStarts());
	z_ = new HypreParVector(*A_);

	// Workaround to create HypreParMatrix w/ sparsity pattern of A and M
	B_ = new HypreParMatrix(hypre_ParCSRMatrixAdd(*M_, *A_));
}


void DiffusionOperator::Assemble_b_Vector(const double &t) const
{
	if (t > 0) {
		f_->SetTime(t);
	}
	else {
		f_->SetTime(this->GetTime());
	}
	b_form_->Assemble();

	if (b_ == NULL) {
		b_ = b_form_->ParallelAssemble();
	}
	else {
		b_form_->ParallelAssemble(*b_);
	}
}


void DiffusionOperator::Mult(const Vector &x, Vector &y) const
{

}


void DiffusionOperator::ImplicitSolve(const double dt, const Vector &x, Vector &k)
{
	// Set HypreParVectors x_,k_ to point to data in Vectors x,k
	x_->SetData(x.GetData());
	k_->SetData(k.GetData());

	// Update b vector and operators to be at the current time 
	Assemble_b_Vector();
	if (dt != current_dt_) {
		UpdateSpatialSolvers(dt);
	}

	// Solve system (M + dtA)k = -Ax + b.
	A_->Mult(*x_, *z_);
	(*z_) *= -1.0;
	if (b_ != NULL) {
		(*z_) += (*b_);	     // z = -Ax + b
	}
	B_pcg_->Mult(*z_, *k_);  // k = (M + dtA)^(-1)z
}


void DiffusionOperator::GetResidual(Vector &u_stop, Vector &u_start, Vector &temp,
									const double &dt, const double &t, const int &update_ops)
{

	// Update operators to correct mesh if needed, update RHS
	if (update_ops == 1) {	
		delete A_;
		delete M_;
		delete B_;

		a_form_->BilinearForm::operator=(0.0);  // not directly inherited
		m_form_->BilinearForm::operator=(0.0);  // from BilinearForm
		a_form_->Assemble(0);
		m_form_->Assemble(0);
		a_form_->EliminateEssentialBC(ess_bdr_);
		m_form_->EliminateEssentialBC(ess_bdr_);
		a_form_->Finalize(0);
		m_form_->Finalize(0);
		A_ = a_form_->ParallelAssemble();
		M_ = m_form_->ParallelAssemble();
		B_ = new HypreParMatrix(hypre_ParCSRMatrixAdd(*M_, *A_));
		UpdateSpatialSolvers(dt);
	}

	Assemble_b_Vector(t);

	// Compute residual r = b_(i+1) - (M+dtA)u_(i+1) + Mu_i
	B_->Mult(u_stop, temp);
	M_->Mult(u_start, u_start);
	u_start -= temp; 
	u_start += *b_;
}


void DiffusionOperator::UpdateSpatialSolvers(const double & dt)
{
	// If new dt was provided, update class dt. Default dt = -1,
	// and class dt will not change if a valid dt is not provided.
	if (dt > 0){
		current_dt_ = dt;
	}

	hypre_ParCSRMatrixSetConstantValues(*B_, 0.0);
	hypre_ParCSRMatrixSum(*B_, 1.0, *M_);
	hypre_ParCSRMatrixSum(*B_, current_dt_, *A_);

	delete B_amg_;
	B_amg_ = new HypreBoomerAMG(*B_);
	B_amg_->SetPrintLevel(-1);

	delete B_pcg_;
	B_pcg_ = new HyprePCG(*B_);
	B_pcg_->SetTol(1e-14);
	B_pcg_->SetMaxIter(2000);
	B_pcg_->SetPrintLevel(0);
	B_pcg_->SetPreconditioner(*B_amg_);
	B_pcg_->SetZeroInintialIterate();
}


void DiffusionOperator::UpdateMeshMonitor(const Vector &u_curr)
{
	// Use central difference to approximate derivative of current solution.
	double u_x;
	for (int i=1; i<(numNodes_-1); i++) {
		u_x = (u_curr(i+1) - u_curr(i-1)) / (mesh_nodes_(i+1) - mesh_nodes_(i-1));
		mesh_monitor_(i) = sqrt(1.0 + u_x*u_x);
	}
	// Use second order forward/backward difference on first and last node, respectively. 
	u_x = (-1.5*u_curr(0) + 2.0*u_curr(1) - 0.5*u_curr(2)) / (mesh_nodes_(1) - mesh_nodes_(0));
	mesh_monitor_(0) = sqrt(1.0 + u_x*u_x);
	u_x = (0.5*u_curr(numNodes_-3) - 2.0*u_curr(numNodes_-2) + 1.5*u_curr(numNodes_-1)) / 
		  (mesh_nodes_(numNodes_-1) - mesh_nodes_(numNodes_-2));
	mesh_monitor_(numNodes_-1) = sqrt(1.0 + u_x*u_x);
}


void DiffusionOperator::UpdateMeshSolver()
{
	delete Amesh_;
	delete M_amg_;
	delete M_pcg_;

	double a0 = 0.0,
		   a1 = 0.0,
		   a2 = 0.0;
	SparseMatrix &A = mesh_form_->SpMat();

	// This is specific to 1-dimension in regards to boundary nodes being
	// first and last in vector, and in mesh movement equation as well. 
	for (int i=1; i<(numNodes_-1); i++) {

		a1 = 1 + current_dt_*(mesh_monitor_(i-1) + 2.0*mesh_monitor_(i) + mesh_monitor_(i+1)) /
			  (2.0*mesh_tau_*dzeta_*dzeta_);
		A.Set(i,i,a1);

		a0 = -current_dt_*(mesh_monitor_(i-1) + mesh_monitor_(i)) /
		  (2.0*mesh_tau_*dzeta_*dzeta_);
		A.Set(i,i-1,a0);

		a2 = -current_dt_*(mesh_monitor_(i) + mesh_monitor_(i+1)) /
		  (2.0*mesh_tau_*dzeta_*dzeta_);
		A.Set(i,i+1,a2);

	}
	// Fix boundaries.
	A.Set(0,0,1.0);
	A.Set(0,1,0.0);
	A.Set(numNodes_-1,numNodes_-1,1.0);
	A.Set(numNodes_-1,numNodes_-2,-0.0);
	Amesh_ = mesh_form_->ParallelAssemble();

	// Construct Hypre operators to solve moving mesh equation.
	M_amg_ = new HypreBoomerAMG(*Amesh_);
	M_amg_->SetPrintLevel(-1);
        HYPRE_BoomerAMGSetCycleRelaxType((HYPRE_Solver) *M_amg_, 3, 3);
        HYPRE_BoomerAMGSetCycleNumSweeps((HYPRE_Solver) *M_amg_, 1, 3);
	M_pcg_ = new HyprePCG(*Amesh_);
	M_pcg_->SetTol(1e-14);
	M_pcg_->SetMaxIter(2000);
	M_pcg_->SetPrintLevel(-1);
	M_pcg_->SetPreconditioner(*M_amg_);
}


void DiffusionOperator::InterpolateSolution(Vector &u_curr)
{
	double temp[numNodes_-2];
	double dx,
		   new_node;
	int    lower,
		   upper;

	// Find j s.t. new_node(i) in [old_node(j),old_node(j+1)]
	for (int i=1; i<(numNodes_-1); i++) {
		if (abs(mesh_displace_(i)) < 1e-6) {
			temp[i-1] = u_curr(i);
		}
		else {
			new_node = mesh_nodes_(i) + mesh_displace_(i);
			if (mesh_displace_(i) < 0) {
				upper = i;
				while (mesh_nodes_(upper) > new_node) {
					upper -= 1;
				}
				lower  = upper; 
				upper += 1;
			}
			else {
				lower = i;
				while (mesh_nodes_(lower) < new_node) {
					lower += 1;
				}
				upper  = lower;
				lower -= 1;
			}
			// Interpolate solution linearly to new mesh. 
			dx		  = mesh_nodes_(upper) - mesh_nodes_(lower);
			temp[i-1] = ( (new_node - mesh_nodes_(lower))*u_curr(upper) +
						  (mesh_nodes_(upper) - new_node)*u_curr(lower) ) / dx;
		}
	}
	for (int i=1; i<(numNodes_-1); i++) {
		u_curr(i) = temp[i-1];
	}
}


void DiffusionOperator::UpdateMesh(ParMesh *pmesh, Vector &u_curr, int &update_ops,
								   int &move_mesh, const double &dt)
{
	if( dt != current_dt_) {
		update_ops = 1;
		current_dt_ = dt;
	}

	if (move_mesh == 1) {
		// Get nodes, update MMPDE linear solver.
		pmesh->GetNodes(mesh_nodes_);
		for (int i=0; i<numNodes_; i++) {
			mesh_displace_(i) = mesh_nodes_(i);
		} 
		UpdateMeshMonitor(u_curr);
		UpdateMeshSolver();

		// Solve MMPDE to get new mesh location, and solve for displacement.
		M_pcg_->Mult(mesh_nodes_,mesh_displace_);
		mesh_displace_ -= mesh_nodes_; 

		if ( abs(mesh_displace_.Max()) > 1e-6 || 
			 abs(mesh_displace_.Min()) > 1e-6 ) {
			update_ops = 1;
			pmesh->MoveNodes(mesh_displace_);
			InterpolateSolution(u_curr);
		}
		else {
			move_mesh = 0;
		}
	}

	// Update mesh and interpolate solution if significant displacement has been
	// made, or operators in class are not up to date with current vector mesh. 
	if (update_ops == 1) {	

		delete A_;
		delete M_;
		delete B_;

		a_form_->BilinearForm::operator=(0.0);  // not directly inherited
		m_form_->BilinearForm::operator=(0.0);  // from BilinearForm
		a_form_->Assemble(0);
		m_form_->Assemble(0);
		a_form_->EliminateEssentialBC(ess_bdr_);
		m_form_->EliminateEssentialBC(ess_bdr_);
		a_form_->Finalize(0);
		m_form_->Finalize(0);
		if (b_ == NULL) {
			b_ = b_form_->ParallelAssemble();
		}
		else {
			b_form_->ParallelAssemble(*b_);
		}
		A_ = a_form_->ParallelAssemble();
		M_ = m_form_->ParallelAssemble();
		B_ = new HypreParMatrix(hypre_ParCSRMatrixAdd(*M_, *A_));

		UpdateSpatialSolvers();
	}
}


DiffusionOperator::~DiffusionOperator()
{
	delete Amesh_;
	delete f_;
	delete alpha_;
	delete A_;
	delete M_;
	delete B_;
	delete b_;
	delete z_;
	delete x_;
	delete k_;
	delete B_amg_;
	delete B_pcg_;
	delete M_amg_;
	delete M_pcg_;
	delete a_form_;
	delete m_form_;
	delete b_form_;
	delete mesh_form_;
};



/* -------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------- */


double InitialConditions(Vector &x0)
{
	int dim = x0.Size();
	double scale = 1.0;
	switch (dim) {
		case 1:
			// if ( x0(0) == 0 || x0(0) == 1) {
			// 	return 0.;
			// }
			// else {
			// 	return scale*sin(3.14159265359*x0(0));				
			// }
			return 0.;
		case 2:
			// return sin(x(0)) + sin(x(1));
			return 0.;
		default:
			return 0.;
	}
}


/* Gaussian bump in (x,t) with given center and width in each dimensions. */
double GuassianBlip(const double &x, const double &c, const double &w, const double &scale)
{
	double e = 2.718281828;
	if (abs(x-c) < w) {
		double d = (x-c)/w;
		double k = -1./(1 - d*d);
		return scale*pow(e,k);		
	}
	else {
		return 0.;
	}
}


/* Forcing function f(x,t). */
void Forcing(const Vector &x, double t, Vector &f) 
{
	double t_cent  = 0.2,
		   t_width = 0.1,
		   t_scale = 1.0,
		   x_cent  = (t+0.25)/2,
		   x_width = 0.05,
		   x_scale = 1.0;

	int dim = x.Size();

	for (int i=0; i<dim; i++) {
		f(i) = 0.0;

	// Gaussian sources moving across spatial domain over time.
	if (forcing_eqn == 1) {
		for (int i=0; i<dim; i++) {
			f(i) = GuassianBlip(x(i), x_cent=(t+0.25)/2, x_width=0.05, x_scale=50.0);
		}
	}
	// Five time-dependent Gaussian sources.
	else if (forcing_eqn == 2) {
		for (int i=0; i<dim; i++) {
			f(i) = 0.0;
			f(i) += GuassianBlip(t,    t_cent=0.1,  t_width=0.05, t_scale=50)*
					GuassianBlip(x(i), x_cent=0.9,  x_width=0.05, x_scale=30);
			f(i) += GuassianBlip(t,    t_cent=0.25, t_width=0.2,  t_scale=30)*
					GuassianBlip(x(i), x_cent=0.3,  x_width=0.15, x_scale=30);
			f(i) += GuassianBlip(t,    t_cent=0.6,  t_width=0.1,  t_scale=10)*
					GuassianBlip(x(i), x_cent=0.5,  x_width=0.3,  x_scale=20);
			f(i) += GuassianBlip(t,    t_cent=0.8,  t_width=0.3,  t_scale=40)*
					GuassianBlip(x(i), x_cent=0.8,  x_width=0.1,  x_scale=30);
			f(i) += GuassianBlip(t,    t_cent=0.9, t_width=0.1,  t_scale=30)*
					GuassianBlip(x(i), x_cent=0.3,  x_width=0.2,  x_scale=30);
		}
	}
	else {
		f = 0.0;
	}

	// Single Gaussian source in space time. 
	#if 0
		f(i) += GuassianBlip(t, t_cent, t_width)*GuassianBlip(x(i),x_cent=(t+0.25)/2,x_width=0.05,);
	#endif			


	}


}



/* -------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------- */




