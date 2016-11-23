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


// Driver:        drive-diffusion-1D-moving-mesh-serial.cpp
//
// Interface:     C++, through MFEM 
//
// Requires:      MFEM, Hypre, Metis and GlVis
//                Modify Makefile to point to metis, mfem and hypre libraries
//
// Compile with:  make drive-diffusion-1D-moving-mesh-serial
//
// Help with:     drive-diffusion-1D-moving-mesh-serial -help
//
// Sample run:    ./drive-diffusion-1D-moving-mesh-serial
//
// Description:   This code runs a serial simulation of the 1D heat equation, with a
//                moving mesh that adapts to the forcing function so that the mesh
//                equidistributes the arc-length of the solution 
//
//                See the write-up in 
//                braid_notes/summer_students/2015/Southworth/MovingMesh_Diffusion
//                for more details


#include "mfem.hpp"

//Extra Hypre Functions
#include "hypre_extra.hpp"
#include <fstream>

using namespace hypre;
using namespace std;
using namespace mfem;


void WriteVector(Vector &data, double t, int type)
{
	ostringstream filestream;
	if (type == 1) {
		filestream << "./sol/Vector_" << t << ".txt";
	}
	else {
		filestream << "./sol/Mesh_" << t << ".txt";
	}

	string filename = filestream.str();
    ofstream output_file(filename.c_str());

    if (output_file.is_open()) {
    	for (int i=0; i<data.Size(); i++) {
    		output_file << setprecision(18) << data[i] << " ";
    	}
	}
}



/* Class to take time steps and move mesh. */
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
								mesh_displace_,
								mesh_nodes_;
	double 						current_dt_;
	const double				dzeta_,
								mesh_tau_;
	int 						numNodes_;
	Array<int> 					ess_bdr_;

	void UpdateSpatialSolvers(const double & dt = -1.0);
	void Assemble_b_Vector() const;
	// Arclength monitor function.
	void UpdateMeshMonitor(const Vector &u_curr);
	void UpdateMeshSolver();
	void InterpolateSolution(Vector &u_curr);

public:

	DiffusionOperator(ParMesh *pmesh, ParFiniteElementSpace *fespace, HypreParVector *U,
					  const double &dzeta, const double &alpha = 1.0, const double &tau = 1.0);
	void UpdateMesh(ParMesh *pmesh, ParFiniteElementSpace *fespace, Vector &u_curr, const double &dt);
	virtual void Mult(const Vector &x, Vector &y) const;
	virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k);
	virtual ~DiffusionOperator();
};


/* Initial conditions of vector. */
double InitialConditions(Vector &x0);


/* Gaussian bump in (x,t) with given center and width in each dimensions. */
double GuassianBlip(const double &x, const double &c, const double &w, const double &scale);


/* Forcing function f(x,t). */
void Forcing(const Vector &x, double t, Vector &f);


// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //


int main(int argc, char* argv[]) 
{

	// Initialize MPI.
	int num_procs, myid;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	int    arg_index		= 0,
		   num_interval		= 10,
		   order			= 1,
		   dim				= 1,
		   ode_solver_type  = 1,
		   vis_steps 		= 1,
		   precision		= 8,
		   num_time			= 10;

	bool visualization = true,
		 visit 		   = false,
		 move_mesh	   = true;

	double max_x 		= 1.0,
		   dt 			= 0.1,
		   t_final		= 1.0,
		   alpha	 	= 1.0,
		   tau 			= 1.0,
		   dzeta;

	// Parse command line.
	while( arg_index < argc ){
		if( strcmp(argv[arg_index], "-n") == 0 ){
			arg_index++;
			num_interval = atoi(argv[arg_index++]);
		}
		else if( strcmp(argv[arg_index], "-tf") == 0 ){
			arg_index++;
			t_final = atof(argv[arg_index++]);
		}
		else if( strcmp(argv[arg_index], "-nt") == 0 ){
			arg_index++;
			num_time = atof(argv[arg_index++]);
		}
		else if( strcmp(argv[arg_index], "-alpha") == 0 ){
			arg_index++;
			alpha = atof(argv[arg_index++]);
		}
		else if( strcmp(argv[arg_index], "-tau") == 0 ){
			arg_index++;
			tau = atof(argv[arg_index++]);
		}
		else if( strcmp(argv[arg_index], "-mesh") == 0 ){
			arg_index++;
			move_mesh = atoi(argv[arg_index++]);
		}
		else{
			if(arg_index > 1){
				printf("UNUSED command line paramter %s\n", argv[arg_index]);
			}
			arg_index++;
		}
	}
	dt 	  = double(t_final) / num_time;
	dzeta = max_x / num_interval;

	// Create mesh.
	Mesh 	mesh(num_interval,max_x);
	ParMesh pmesh(MPI_COMM_WORLD,mesh);

	// Create finite element collection and space.
	H1_FECollection fec(order, dim);
	ParFiniteElementSpace fespace(&pmesh, &fec);
	cout << "Number of unknowns: " << fespace.GetVSize() << endl;

	// Assemble initial conditions.
	FunctionCoefficient u0(InitialConditions);
	ParGridFunction *u = new ParGridFunction(&fespace);
	u->ProjectCoefficient(u0);
	HypreParVector *U = u->GetTrueDofs();

// -------------------------------------------------------------------------- //

	VisItDataCollection visit_dc("Example9-Parallel", &pmesh);
	if (visit) {
		ostringstream mesh_name, sol_name;
		mesh_name << "ex9-mesh." << setfill('0') << setw(6) << myid;
		sol_name << "ex9-init." << setfill('0') << setw(6) << myid;
		ofstream omesh(mesh_name.str().c_str());
		omesh.precision(precision);
		pmesh.Print(omesh);
		ofstream osol(sol_name.str().c_str());
		osol.precision(precision);
		u->Save(osol);

		visit_dc.RegisterField("solution", u);
		visit_dc.SetCycle(0);
		visit_dc.SetTime(0.0);
		visit_dc.Save();
	}

	socketstream sout;
	if (visualization) {
		char vishost[] = "localhost";
		int  visport   = 19916;
		sout.open(vishost, visport);
		if (!sout) {
			if (myid == 0) {
				cout << "Unable to connect to GLVis server at "
				<< vishost << ':' << visport << endl
				<< "GLVis visualization disabled.\n";  	
			}
			visualization = false;
		}
		else {
			sout << "parallel " << num_procs << " " << myid << "\n";
			sout.precision(precision);
			sout << "solution\n" << pmesh << *u;
			sout << "pause\n";
			sout << flush;
			if (myid == 0) {
				cout << "GLVis visualization paused."
				<< " Press space (in the GLVis window) to resume it.\n";         	
			}
		}
	}


// -------------------------------------------------------------------------- //

	// Define and initialize ODE solver used for time integration on level l.
	ODESolver *ode_solver = NULL;
	switch (ode_solver_type)
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
	DiffusionOperator diff(&pmesh, &fespace, U, dzeta, alpha, tau);
	ode_solver->Init(diff);

	double t 	 = 0.0;
	int timestep = 0;

	// Vector type pointing in memory to current step. Used in moving mesh.
	Vector *u_vec = (Vector*) U;
	Vector nodes; 
	pmesh.GetNodes(nodes);
	while (t < t_final - dt/2.)
	{

#if 0
ofstream output_file;
ostringstream grid_stream;
grid_stream << "./mesh_data/Serial_grid.csv";
string grid_name = grid_stream.str();
output_file.open(grid_name.c_str(), ios::app);
output_file << setprecision(6) << t << ",";
for (int i=0; i<nodes.Size(); i++) {
	output_file << setprecision(14) << nodes(i) << ",";
}
output_file << "\n";
output_file.close();

ostringstream vec_stream;
vec_stream << "./mesh_data/Serial_vec.csv";
string vec_name = vec_stream.str();
output_file.open(vec_name.c_str(), ios::app);
output_file << setprecision(6) << t << ",";
for (int i=0; i<u_vec->Size(); i++) {
	output_file << setprecision(14) << (*u_vec)(i) << ",";
}
output_file << "\n";
output_file.close();
#endif

		// Move mesh.
		if (move_mesh == true) {
			WriteVector(*u_vec, t, 1);
			WriteVector(nodes, t, 0);
			diff.UpdateMesh(&pmesh, &fespace, *u_vec, dt);
			// u_vec->Print();
			pmesh.GetNodes(nodes);
			// nodes.Print();
			// cout << endl;		
		}

		// Update u,t one time step of size dt.
		ode_solver->Step(*U, t, dt);
		timestep += 1;

		if (timestep % vis_steps == 0) {

			// 11. Extract the parallel grid function corresponding to the finite
			//     element approximation U (the local solution on each processor).
			*u = *U;

			if (visualization) {

				sout << "parallel " << num_procs << " " << myid << "\n";
				sout << "solution\n" << pmesh << *u << "pause\n" << flush;
			}

			if (visit) {
				visit_dc.SetCycle(timestep);
				visit_dc.SetTime(t);
				visit_dc.Save();
			}
		}
	}


#if 0
ofstream output_file;
ostringstream grid_stream;
grid_stream << "./mesh_data/Serial_grid.csv";
string grid_name = grid_stream.str();
output_file.open(grid_name.c_str(), ios::app);
output_file << setprecision(6) << t << ",";
for (int i=0; i<nodes.Size(); i++) {
	output_file << setprecision(14) << nodes(i) << ",";
}
output_file << "\n";
output_file.close();

ostringstream vec_stream;
vec_stream << "./mesh_data/Serial_vec.csv";
string vec_name = vec_stream.str();
output_file.open(vec_name.c_str(), ios::app);
output_file << setprecision(6) << t << ",";
for (int i=0; i<u_vec->Size(); i++) {
	output_file << setprecision(14) << (*u_vec)(i) << ",";
}
output_file << "\n";
output_file.close();
#endif



	pmesh.GetNodes(nodes);
	WriteVector(*u_vec, t, 1);
	WriteVector(nodes, t, 0);

	cout << "\nMesh at t = " << t << endl;
	nodes.Print();
	cout << "\nSolution at t = " << t << endl;
	u_vec->Print();


	delete U;
	delete u;
	delete ode_solver;
	return 0;
}


// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //


DiffusionOperator::DiffusionOperator(ParMesh *pmesh, ParFiniteElementSpace *fespace,
									 HypreParVector *U, const double &dzeta,
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
	a_form_->EliminateEssentialBC(ess_bdr_, *U, *b_form_);
	m_form_->EliminateEssentialBC(ess_bdr_, *U, *b_form_);

	// Construct bilinear form for mesh equations. Add 1d diffusion to make
	// a tridiagonal operator.
	mesh_form_ = new ParBilinearForm(fespace);
	mesh_form_->AddDomainIntegrator(new DiffusionIntegrator(*alpha_));
	mesh_form_->Assemble();
	mesh_form_->Finalize(0);

	// Finalize bilinear forms and form parallel operators. 
	a_form_->Finalize();
	m_form_->Finalize();
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


void DiffusionOperator::Assemble_b_Vector() const
{
	// Check that this runs properly w/ this->, previously did not have it.
	f_->SetTime(this->GetTime());
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


void DiffusionOperator::UpdateMesh(ParMesh *pmesh, ParFiniteElementSpace *fespace,
								   Vector &u_curr, const double &dt)
{
	int update_ops = 0;

	if( dt != current_dt_) {
		update_ops = 1;
		current_dt_ = dt;
	}

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
	}

	// Update mesh and interpolate solution if significant displacement has been
	// made. Update operators correspondingly. 
	if (update_ops == 1) {
		InterpolateSolution(u_curr);
		pmesh->MoveNodes(mesh_displace_);

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



// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //


double InitialConditions(Vector &x0)
{
	int dim = x0.Size();
	switch (dim) {
		case 1:
	      // double scale = 1.0;
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

