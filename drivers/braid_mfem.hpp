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

#ifndef braid_mfem_HEADER
#define braid_mfem_HEADER

#include "braid.hpp"
#include "mfem.hpp"
#include <fstream>
#include <cstdlib> // rand, srand

using namespace mfem;


class BraidVector : public HypreParVector
{
public:
   // The spatial coarsening level; level 0 is the finest level.
   int level;

   /// Construct a BraidVector given a level and a ParFiniteElementSpace.
   BraidVector(int source_level, ParFiniteElementSpace *pfes)
      : HypreParVector(pfes), level(source_level) { }

   /** "Copy" constructor: create a BraidVector compatible with source_vector
        without copying the data. */
   BraidVector(const BraidVector &source_vector)
      : HypreParVector(source_vector), level(source_vector.level) { }

   /** Construct a BraidVector compatible with source_vector (a HypreParVector)
        and the level set to source_level. */
   BraidVector(int source_level, const HypreParVector &source_vector)
      : HypreParVector(source_vector), level(source_level) { }

   /// Copy the data from source_vector (a Vector).
   BraidVector &operator=(const Vector &source_vector)
   { Vector::operator=(source_vector); return *this; }

   /// Initialize the Braidvector with a double value.
   BraidVector &operator=(double value)
   { Vector::operator=(value); return *this; }

   /// Clone function (copy the data as well into the new BraidVector)
   BraidVector *clone()
   {
      BraidVector *y = new BraidVector(*this);
      *y = *this;
      return y;
   }
};


// Wrapper for BRAID's App object
class MFEMBraidApp : public BraidApp
{
protected:
   // BraidApp defines tstart, tstop, ntime and comm_t

   // Data for multiple spatial levels, level 0 is the finest level
   Array<ParMesh *>               mesh;
   Array<ParFiniteElementSpace *> fe_space;
   Array<ParGridFunction *>       x;  // auxiliary ParGridFunctions
   Array<const SparseMatrix *>    P;  // local prolongation matrices, l+1 --> l
   Array<TimeDependentOperator *> ode;
   Array<ODESolver *>        solver;
   Array<int>                buff_size;
   Array<double>             max_dt; // maximal safe dt on each spatial mesh

   HypreParVector *X0;    // Initial condition (at the finest level 0)
   bool init_rand;     /* If true, use std::rand() to initialize BraidVectors in
                          the Init() method. The default value is false:
                          initialize with 0.0. */

   // ownership of mesh, fe_space, x, R, ode, solver, and X0
   bool own_data;

   // Exact solution at the final time (optional); used in the Access method
   Coefficient *exact_sol;

   // Used for output to GLVIS in the Access method
   socketstream sol_sock;
   const char  *vishost;
   int          visport;

   static const char *vishost_default;
   static const int   visport_default;

   int vis_time_steps;
   int vis_braid_steps;
   bool vis_screenshots;

   // Return the level index l such that max_dt[l-1] < dt <= max_dt[l],
   // This index corresponds to the finest compatible mesh for given time step
   int ComputeSpaceLevel(double tstart, double tprior, double tstop);

   static int EvalBufSize(int vector_size)
   {
      // Add 1 to the buffer size to account for BraidVector::level
      return sizeof(double) * (vector_size + 1);
   }

   // Construct MFEMBraidApp with empty (spatial) multilevel structures. Can be
   // used in derived classes for basic initialization.
   MFEMBraidApp(MPI_Comm comm_t_,
                double tstart_ = 0.0, double tstop_ = 1.0, int ntime_ = 100);

   // This method can be used by derived classes to initialize the multilevel
   // structures: mesh, fe_space, etc. The MFEMBraidApp assumes ownership of the
   // ParMesh object and all other constructed multilevel structures.
   //
   // It uses the virtual methods AllocLevels, ConstructFESpace and InitLevel
   // which must be defined by the derived class.
   void InitMultilevelApp(ParMesh *pmesh, int pref, bool scoarsen);

   // Allocate data structures for the given number of spatial levels. Used by
   // InitMultilevelApp.
   //
   // This virtual method can be re-defined by derived classes that use the
   // InitMultilevelApp method to allocate structures that depend on the actual
   // number of spatial levels as determined by InitMultilevelApp.
   virtual void AllocLevels(int num_levels) { }

   // Construct the ParFiniteElmentSpace for the given mesh. Used by
   // InitMultilevelApp.
   //
   // This virtual method must be defined by derived classes that use the
   // InitMultilevelApp method.
   virtual ParFiniteElementSpace *ConstructFESpace(ParMesh *pmesh);

   // Assuming mesh[l] and fe_space[l] are set, initialize, ode[l], solver[l],
   // and max_dt[l]. Used by InitMultilevelApp.
   //
   // This virtual method must be defined by derived classes that use the
   // InitMultilevelApp method.
   virtual void InitLevel(int l);

public:
   /// Construct an MFEMBraidApp with one spatial level.
   MFEMBraidApp(MPI_Comm comm_t_, TimeDependentOperator *ode_,
                HypreParVector *X0_, ParGridFunction *x_, ODESolver *solver_,
                double tstart_ = 0.0, double tstop_ = 1.0, int ntime_ = 100);

   /** Begin the construction of an MFEMBraidApp with multiple spatial levels.
       To complete the construction, call SetSpaceLevel (for all spatial levels)
       and SetInitialCondition. */
   MFEMBraidApp(MPI_Comm comm_t_, const int num_space_levels,
                double tstart_ = 0.0, double tstop_ = 1.0, int ntime_ = 100);

   virtual ~MFEMBraidApp();

   int GetNumSpaceLevels() const { return ode.Size(); }

   void SetSpaceLevel(int l, TimeDependentOperator *ode_l,
                      ODESolver *solver_l, ParGridFunction *x_l,
                      SparseMatrix *P_l, double max_dt_l);

   void SetInitialCondition(HypreParVector *_X0) { X0 = _X0; }

   void SetRandomInitVectors(unsigned seed)
   {
      init_rand = true;
      std::srand(seed);
   }

   void SetExactSolution(Coefficient *exsol) { exact_sol = exsol; }

   void SetVisHostAndPort(const char *vh, int vp);

   void SetVisSampling(int time_steps, int braid_steps)
   {
      vis_time_steps = time_steps;
      vis_braid_steps = braid_steps;
   }

   void SetVisScreenshots(bool ss) { vis_screenshots = ss; }

   // Below braid_Vector == BraidVector*

   virtual int Step(braid_Vector    u_,
                    braid_Vector    ustop_,
                    braid_Vector    fstop_,
                    BraidStepStatus &pstatus);

   virtual int Clone(braid_Vector  u_,
                     braid_Vector *v_ptr);

   virtual int Init(double        t,
                    braid_Vector *u_ptr);

   virtual int Free(braid_Vector u_);

   virtual int Sum(double       alpha,
                   braid_Vector x_,
                   double       beta,
                   braid_Vector y_);

   virtual int SpatialNorm(braid_Vector  u_,
                           double       *norm_ptr);

   virtual int BufSize(int *size_ptr,
                       BraidBufferStatus  &status);

   virtual int BufPack(braid_Vector  u_,
                       void         *buffer,
                       BraidBufferStatus  &status);

   virtual int BufUnpack(void         *buffer,
                         braid_Vector *u_ptr,
                         BraidBufferStatus  &status);

   virtual int Coarsen(braid_Vector   fu_,
                       braid_Vector  *cu_ptr,
                       BraidCoarsenRefStatus &status);

   virtual int Refine(braid_Vector   cu_,
                      braid_Vector  *fu_ptr,
                      BraidCoarsenRefStatus &status);

   virtual int Access(braid_Vector       u_,
                      BraidAccessStatus &astatus);

   virtual int Residual(braid_Vector     u_,
                        braid_Vector     r_,
                        BraidStepStatus &pstatus) { return 0; }
};


struct BraidOptions : public OptionsParser
{
   double t_start;
   double t_final;
   int    num_time_steps;
   int    num_procs_x;
   int    skip;
   int    max_levels;
   int    min_coarse;
   int    nrelax;
   int    nrelax0;
   double tol;
   bool   rtol;
   int    tnorm;
   int    storage;
   int    cfactor;
   int    cfactor0;
   int    max_iter;
   int    nfmg_Vcyc;
   bool   spatial_coarsen;
   int    access_level;
   int    print_level;
   int    use_seq_soln;

   // (optional) mesh and refinement options
   const char *mesh_file;
   int         ser_ref_levels;
   int         par_ref_levels;

   BraidOptions(int argc, char *argv[]);

   // Add (enable) the mesh and refinement options
   void AddMeshOptions();

   Mesh *LoadMeshAndSerialRefine();

   ParMesh *LoadMeshAndRefine(MPI_Comm comm_x);

   virtual void SetBraidCoreOptions(BraidCore &core);
};

// Stores and prints runtime information regarding spatial coarsening.
class SpaceTimeMeshInfo
{
public:

   // This table is (num Braid levels) x 5 table, where each row is a Braid level
   // and the columns store, respectively, the spatial level used at that Braid level,
   // the max mesh width h_max, the min mesh width h_min and the time step size dt,
   Array<double> mesh_table;
   Array<double> mesh_table_global;
   int           max_levels;

   // Unused rows in the table are -1.0
   SpaceTimeMeshInfo(int _max_levels) :
      mesh_table(4*_max_levels),
      mesh_table_global(4*_max_levels),
      max_levels(_max_levels)
   {
      mesh_table = -1.0;
      mesh_table_global = -1.0;
   }

   // Reinitialize the arrays
   void Reinitialize(int _max_levels);

   // Helper function to compute mesh size
   void ComputeMeshSize( ParMesh *pmesh, double * h_min_ptr, double * h_max_ptr);

   // Fill in a row of the table
   void SetRow(int braid_level, int vec_level, ParMesh * pmesh, double dt);

   // Print the table to screen, reducing the entries over comm
   void Print(MPI_Comm comm);
};



// Implementations

// Initialize static data members
const char *MFEMBraidApp::vishost_default = "localhost";
const int   MFEMBraidApp::visport_default = 19916;

// Construct MFEMBraidApp with empty (spatial) multilevel structures
MFEMBraidApp::MFEMBraidApp(
      MPI_Comm comm_t_, double tstart_, double tstop_, int ntime_)

   : BraidApp(comm_t_, tstart_, tstop_, ntime_)
{
   X0 = NULL;
   init_rand = false;
   own_data = false;

   exact_sol = NULL;

   vishost = vishost_default;
   visport = visport_default;
}

// Construct an MFEMBriadApp with one spatial level
MFEMBraidApp::MFEMBraidApp(
      MPI_Comm comm_t_, TimeDependentOperator *ode_,
      HypreParVector *X0_, ParGridFunction *x_, ODESolver *solver_,
      double tstart_, double tstop_, int ntime_)

   : BraidApp(comm_t_, tstart_, tstop_, ntime_),
     mesh(1), fe_space(1), x(1), ode(1), solver(1), buff_size(1), max_dt(1),
     X0(X0_), init_rand(false)
{
   x[0] = x_;
   fe_space[0] = x[0]->ParFESpace();
   mesh[0] = fe_space[0]->GetParMesh();
   ode[0] = ode_;
   solver[0] = solver_;
   solver[0]->Init(*ode[0]);
   max_dt[0] = 2 * (tstop_ - tstart_);
   buff_size[0] = EvalBufSize(fe_space[0]->TrueVSize());
   own_data = false;

   exact_sol = NULL;

   vishost = vishost_default;
   visport = visport_default;
}

// Begin the construction of an MFEMBriadApp with multiple spatial levels
MFEMBraidApp::MFEMBraidApp(
      MPI_Comm comm_t_, const int num_space_levels,
      double tstart_, double tstop_, int ntime_)

   : BraidApp(comm_t_, tstart_, tstop_, ntime_),
     mesh(num_space_levels), fe_space(num_space_levels), x(num_space_levels),
     P(num_space_levels - 1), ode(num_space_levels), solver(num_space_levels),
     buff_size(num_space_levels), max_dt(num_space_levels)
{
   // initialize the arrays element-wise
   mesh = NULL;
   fe_space = NULL;
   x = NULL;
   P = NULL;
   ode = NULL;
   solver = NULL;
   buff_size = 0;
   max_dt = 0.0;

   X0 = NULL;
   init_rand = false;
   own_data = false;
   exact_sol = NULL;

   vishost = vishost_default;
   visport = visport_default;
}

MFEMBraidApp::~MFEMBraidApp()
{
   if (own_data)
   {
      delete X0;
      for (int i = 0; i < mesh.Size(); i++)
      {
         delete solver[i];
         delete ode[i];
         if (i < P.Size()) delete P[i];
         delete x[i];
         delete fe_space[i];
         delete mesh[i];
      }
   }
}

void MFEMBraidApp::InitMultilevelApp(ParMesh *pmesh, int pref, bool scoarsen)
{
   int num_levels = pref+1;

   if (!scoarsen)
   {
      // Parallel Refinement
      for (int l = 0; l < pref; l++)
         pmesh->UniformRefinement();
      num_levels = 1;
   }

   mesh.SetSize(num_levels);
   fe_space.SetSize(num_levels);
   x.SetSize(num_levels);
   P.SetSize(num_levels-1);
   ode.SetSize(num_levels);
   solver.SetSize(num_levels);
   buff_size.SetSize(num_levels);
   max_dt.SetSize(num_levels);
   AllocLevels(num_levels);

   ParFiniteElementSpace *pfes = ConstructFESpace(pmesh);
   for (int l = num_levels-1; l > 0; l--)
   {
      mesh[l] = new ParMesh(*pmesh);
      fe_space[l] = ConstructFESpace(mesh[l]);
      x[l] = new ParGridFunction(fe_space[l]);
      InitLevel(l); // initialize ode[l], solver[l], and max_dt[l]
      buff_size[l] = EvalBufSize(fe_space[l]->TrueVSize());

#if 0 // pre-rebalance-dev version
      pmesh->UseTwoLevelState(1);
      pmesh->UniformRefinement();
      pfes->Update();
      R[l-1] = pfes->GlobalRestrictionMatrix(fe_space[l], 0);
      pmesh->SetState(Mesh::NORMAL);
#else
      pmesh->UniformRefinement();
      pfes->Update();
      P[l-1] = dynamic_cast<const SparseMatrix*>(pfes->GetUpdateOperator());
      MFEM_VERIFY(P[l-1] != NULL, "update operator is not a SparseMatrix");
      pfes->SetUpdateOperatorOwner(false);
      pfes->UpdatesFinished();
#endif
   }
   mesh[0] = pmesh;
   fe_space[0] = pfes;
   x[0] = new ParGridFunction(fe_space[0]);
   InitLevel(0); // initialize ode[0], solver[0], and max_dt[0]
   buff_size[0] = EvalBufSize(fe_space[0]->TrueVSize());

   // Print mesh info
   int myid_t;
   MPI_Comm_rank(comm_t, &myid_t);
   if(myid_t == 0)
      pmesh->PrintInfo();

   // Print size of finest-grid spatial matrix
   int myid, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   size = fe_space[0]->GlobalTrueVSize();
   if (myid == 0)
      std::cout << "\nNumber of spatial unknowns on finest-grid: "
                << size << "\n\n";

   own_data = true;
}

ParFiniteElementSpace *MFEMBraidApp::ConstructFESpace(ParMesh *pmesh)
{
   MFEM_ABORT("this virtual method must be defined by derived classes that use"
              " the InitMultilevelApp method!");
   return NULL;
}

void MFEMBraidApp::InitLevel(int l)
{
   MFEM_ABORT("this virtual method must be defined by derived classes that use"
              " the InitMultilevelApp method!");
}

void MFEMBraidApp::SetSpaceLevel(int l, TimeDependentOperator *ode_l,
                                 ODESolver *solver_l, ParGridFunction *x_l,
                                 SparseMatrix *P_l, double max_dt_l)
{
   x[l] = x_l;
   fe_space[l] = x[l]->ParFESpace();
   mesh[l] = fe_space[l]->GetParMesh();
   ode[l] = ode_l;
   solver[l] = solver_l;
   max_dt[l] = max_dt_l;
   buff_size[l] = EvalBufSize(fe_space[l]->TrueVSize());
   solver[l]->Init(*ode[l]);
   if (l < GetNumSpaceLevels() - 1)
      P[l] = P_l;
}

void MFEMBraidApp::SetVisHostAndPort(const char *vh, int vp)
{
   vishost = vh;
   visport = vp;
}

// Return the level index l such that max_dt[l-1] < dt <= max_dt[l],
// This index corresponds to the finest compatible mesh for given time step
int MFEMBraidApp::ComputeSpaceLevel(double tstart, double tprior, double tstop)
{
   double dt = fmax(tstop - tstart, tstart - tprior);
   int level = GetNumSpaceLevels()-1;
   if (dt > max_dt[level])
   {
      std::cerr << "Warning: incompatible dt = " << dt
                << " for the coarsest spatial level of " << level << std::endl;
   }
   else
   {
      for ( ; level > 0; level--)
      {
         if (dt > max_dt[level-1])
            break;
      }
   }

   return level;
}

int MFEMBraidApp::Step(braid_Vector    u_,
                       braid_Vector    ustop_,
                       braid_Vector    fstop_,
                       BraidStepStatus &pstatus)
{
   BraidVector *u = (BraidVector*) u_;
   int level = u->level;
   double tstart, tstop, t, dt;
   int braid_level;

   // Get time step information
   pstatus.GetTstartTstop(&tstart, &tstop);
   pstatus.GetLevel(&braid_level);

   t = tstart;
   dt = tstop - tstart;
   solver[level]->Step(*u, t, dt);

   // no refinement
   pstatus.SetRFactor(1);

   return 0;
}

int MFEMBraidApp::Clone(braid_Vector  u_,
                        braid_Vector *v_ptr)
{
   BraidVector *u = (BraidVector*) u_;
   BraidVector *v = u->clone();
   *v_ptr = (braid_Vector) v;
   return 0;
}

int MFEMBraidApp::Init(double        t,
                       braid_Vector *u_ptr)
{
   int level = 0;
   BraidVector *u = new BraidVector(level, *X0);
   if (t != tstart)
   {
      if (init_rand)
      {
         const int s = u->Size();
         double *data = u->GetData();
         for (int i = 0; i < s; i++)
         {
            data[i] = double(std::rand())/(RAND_MAX);
         }
      }
      else
         *u = 0.0;
   }
   else
      *u = *X0;

   *u_ptr = (braid_Vector) u;
   return 0;
}

int MFEMBraidApp::Free(braid_Vector u_)
{
   BraidVector *u = (BraidVector*) u_;
   delete u;
   return 0;
}

int MFEMBraidApp::Sum(double       alpha,
                      braid_Vector x_,
                      double       beta,
                      braid_Vector y_)
{
   BraidVector *x = (BraidVector*) x_;
   BraidVector *y = (BraidVector*) y_;
   add(alpha, *x, beta, *y, *y);
   return 0;
}

int MFEMBraidApp::SpatialNorm(braid_Vector  u_,
                              double       *norm_ptr)
{
   double dot;
   BraidVector *u = (BraidVector*) u_;
   dot = InnerProduct(u, u);
   *norm_ptr = sqrt(dot);
   return 0;
}

int MFEMBraidApp::BufSize(int                *size_ptr,
                          BraidBufferStatus  &status)                           
{
   *size_ptr = buff_size[0];
   return 0;
}

int MFEMBraidApp::BufPack(braid_Vector       u_,
                          void               *buffer,
                          BraidBufferStatus  &status)
{
   BraidVector *u = (BraidVector*) u_;
   double *dbuf = (double *) buffer;

   dbuf[0] = (double) u->level;
   memcpy(dbuf + 1, u->GetData(), buff_size[u->level] - sizeof(double));
   status.SetSize( buff_size[u->level] );
   return 0;
}

int MFEMBraidApp::BufUnpack(void              *buffer,
                            braid_Vector      *u_ptr,
                            BraidBufferStatus &status)
{
   double *dbuf = (double *) buffer;
   int level = (int) dbuf[0];
   BraidVector *u = new BraidVector(level, fe_space[level]);
   memcpy(u->GetData(), dbuf + 1, buff_size[u->level] - sizeof(double));
   *u_ptr = (braid_Vector) u;
   return 0;
}

int MFEMBraidApp::Coarsen(braid_Vector   fu_,
                          braid_Vector  *cu_ptr,
                          BraidCoarsenRefStatus &status)
{
   double tstart, c_tprior, c_tstop, f_tprior, f_tstop;
   status.GetTpriorTstop(&tstart, &f_tprior, &f_tstop, &c_tprior, &c_tstop);

   BraidVector *fu = (BraidVector *) fu_;
   int flevel      = fu->level;
   int clevel      = ComputeSpaceLevel(tstart, c_tprior, c_tstop);
   BraidVector *cu = NULL;

   MFEM_VERIFY(flevel == ComputeSpaceLevel(tstart, f_tprior, f_tstop),
               "ComputeSpaceLevel returned incorrect level for fine vector");

   if (clevel < GetNumSpaceLevels() && clevel > flevel)
   {
      for (int lev = flevel+1; lev <= clevel; lev++)
      {
         cu = new BraidVector(lev, fe_space[lev]);

         // Maps true dofs in 'fu' to local dofs in ParGridFunction x,
         // with the caveat that local dofs that fu does not own are set to 0
         // in x.
         fe_space[lev-1]->GetRestrictionMatrix()->MultTranspose(*fu, *x[lev-1]);

         // Rescale before restriction
          *x[lev-1] *= 1./4;

         // Apply local restriction, P^T (no communication)
         P[lev-1]->MultTranspose(*x[lev-1], *x[lev]);

         // Assemble x by condensing it into the HypreParVector
         x[lev]->ParallelAssemble(*cu);

         // Delete fu, unless it is the input vector
         if (lev > flevel+1)
         {
            delete fu;
         }
         // Change fu to point to the "new" fine vector
         fu = cu;
      }
   }
   else if (clevel == flevel)
   {
      cu = fu->clone();
   }
   else
   {
      MFEM_ABORT("ComputeSpaceLevel returned incorrect level for"
                 " coarse vector");
   }

   // check cu for correct size
   MFEM_VERIFY(cu->Size() == ode[cu->level]->Height(),
         "incorrect coarse vector size!");

   *cu_ptr = (braid_Vector) cu;
   return 0;
}

int MFEMBraidApp::Refine(braid_Vector   cu_,
                         braid_Vector  *fu_ptr,
                         BraidCoarsenRefStatus &status)
{
   double tstart, c_tprior, c_tstop, f_tprior, f_tstop;
   status.GetTpriorTstop(&tstart, &f_tprior, &f_tstop, &c_tprior, &c_tstop);

   BraidVector *cu = (BraidVector *) cu_;
   int clevel      = cu->level;
   int flevel      = ComputeSpaceLevel(tstart, f_tprior, f_tstop);
   BraidVector *fu = NULL;

   MFEM_VERIFY(clevel == ComputeSpaceLevel(tstart, c_tprior, c_tstop),
               "ComputeSpaceLevel returned incorrect level for coarse vector");

   if (flevel >= 0 && flevel < clevel)
   {
      for (int lev = clevel-1; lev >=flevel; lev--)
      {
         fu = new BraidVector(lev, fe_space[lev]);

         // Distribute cu into x, including shared dofs requiring communication
         x[lev+1]->Distribute(cu);

         // Apply local interpolation (no communication)
         P[lev]->Mult(*x[lev+1], *x[lev]);

         // Map only the true dofs from x into fu
         x[lev]->GetTrueDofs(*fu);

         // Delete cu, unless it is the input vector
         if (lev < clevel-1)
         {
            delete cu;
         }
         // Change cu to point to the "new" fine vector
         cu = fu;
      }
   }
   else if (clevel == flevel)
   {
      fu = cu->clone();
   }
   else
   {
      MFEM_ABORT("ComputeSpaceLevel returned incorrect level for fine vector");
   }

   // check fu for correct size
   MFEM_VERIFY(fu->Size() == ode[fu->level]->Height(),
         "incorrect fine vector size!");

   *fu_ptr = (braid_Vector) fu;
   return 0;
}

int MFEMBraidApp::Access(braid_Vector       u_,
                         BraidAccessStatus &astatus)
{
   BraidVector *u = (BraidVector*) u_;

   // Extract information from astatus
   int done, level, iter;
   double rnorm, t;
   astatus.GetTILD(&t, &iter, &level, &done);
   astatus.GetResidual(&rnorm);
   int cycle = (int)std::floor((t - tstart) / (tstop - tstart) * ntime + 0.5);

   if ( (level == 0) &&

        ( (vis_time_steps > 0 && cycle % vis_time_steps == 0) ||
          (cycle == ntime) ) &&

        ( (vis_braid_steps > 0 && iter % vis_braid_steps == 0) || done ) )
   {
      (*x[level]) = *u; // Distribute

      // Opening multiple 'parallel' connections to GLVis simultaneously
      // (from different time intervals in this case) may cause incorrect
      // behavior.

      int init_sock = 0;
      if (!sol_sock.is_open())
      {
         sol_sock.open(vishost, visport);
         init_sock = 1;
      }

      int good, all_good;
      good = sol_sock.good();
      MPI_Allreduce(&good, &all_good, 1, MPI_INT, MPI_LAND,
                    mesh[level]->GetComm());

      if (all_good)
      {
         sol_sock << "parallel " << mesh[level]->GetNRanks()
                  << " " << mesh[level]->GetMyRank() << "\n";
         sol_sock << "solution\n";
         sol_sock.precision(8);
         mesh[level]->Print(sol_sock);
         x[level]->Save(sol_sock);

         if (init_sock)
         {
            int comm_t_rank;
            MPI_Comm_rank(comm_t, &comm_t_rank);

            // sol_sock << "valuerange 0 1\n";
            // sol_sock << "autoscale off\n";
            sol_sock << "keys cmAaa\n";
            sol_sock << "window_title 'comm_t rank: " << comm_t_rank
                     << ", t = " << t << "'\n";
            if (vis_screenshots)
            {
               sol_sock << "pause\n";
               if (mesh[level]->GetMyRank() == 0)
                  std::cout << "Visualization paused. Press space (in GLVis) to"
                               " continue." << std::endl;
            }
         }
         sol_sock << std::flush;
         if (vis_screenshots)
         {
            sol_sock << "screenshot braid_" << std::setfill('0') << std::setw(2)
                     << iter << "_" << std::setw(6) << cycle << ".png"
                     << std::endl;
         }
         if (mesh[level]->GetMyRank() == 0)
            std::cout << "Visualization updated." << std::flush;
      }

      if (exact_sol)
      {
         exact_sol->SetTime(t);

         double err = x[level]->ComputeL2Error(*exact_sol);

         if (mesh[level]->GetMyRank() == 0)
            std::cout << " L2 norm of the error = " << err << std::endl;
      }
      else if (mesh[level]->GetMyRank() == 0)
         std::cout << std::endl;
   }

   return 0;
}

BraidOptions::BraidOptions(int argc, char *argv[])
   : OptionsParser(argc, argv)
{
   t_start          = 0.0;
   t_final          = 1.0;
   num_time_steps   = 100;
   num_procs_x      = 1;

   // For the braid defaults, see braid_Init() in braid.c
   // The values below should match the braid defaults.
   skip             = 1;
   max_levels       = 30;
   min_coarse       = 2;
   nrelax           = 1;
   nrelax0          = -1;    // if > -1, use as nrelax for level 0
   tol              = 1e-9;
   rtol             = true;
   tnorm            = 2;
   storage          = -1;
   cfactor          = 2;
   cfactor0         = -1;    // if > -1, use as cfactor for level 0
   max_iter         = 100;
   nfmg_Vcyc        = 0;     // if > 0, enable FMG and use as nfmg_Vcyc
   spatial_coarsen  = false; // if true, enable spatial coarsening
   access_level     = 1;
   print_level      = 2;
   use_seq_soln     = 0;

   AddOption(&t_start, "-ts", "--t-start", "Start time.");
   AddOption(&t_final, "-tf", "--t-final", "Final time.");
   AddOption(&num_time_steps, "-nt", "--num-time-steps",
             "Number of time steps.");
   AddOption(&num_procs_x, "-px", "--num-spatial-procs",
             "Number of processors to use in spatial direction.");
   AddOption(&max_levels, "-ml", "--max-levels",
             "Maximum number of time levels.");
   AddOption(&min_coarse, "-mc", "--min-coarse",
             "Minimum possible coarse level size.");
   AddOption(&nrelax, "-nu", "--num-fc-relax",
             "Number of F-C relaxations.");
   AddOption(&nrelax0, "-nu0", "--num-fc-relax-level-0",
             "Number of F-C relaxations on level 0.");
   AddOption(&tol, "-tol", "--tolerance", "Stopping tolerance.");
   AddOption(&rtol, "-reltol", "--relative-tolerance", "-abstol",
             "--absolute-tolerance",
             "Use relative or absolute stopping tolerance.");
   AddOption(&tnorm, "-tnorm", "--temporal-norm",
             "Temporal norm to use: 1:one-norm, 2:two-norm, or "
             "3:max-norm.");
   AddOption(&cfactor, "-cf", "--coarsen-factor",
             "Coarsening factor.");
   AddOption(&storage, "-store", "--storage-option",
             "Storage to use: 0:store C points, 1:store all points.");
   AddOption(&cfactor0, "-cf0", "--agg-coarsen-factor",
             "Aggressive coarsening factor, -1:off.");
   AddOption(&max_iter, "-mi", "--max-iter",
             "Maximum number of iterations.");
   AddOption(&nfmg_Vcyc, "-fmg", "--fmg-v-cycles",
             "Number of V-cycles to use at each FMG level (0:off).");
   AddOption(&spatial_coarsen, "-sc", "--spatial-coarsen", "-no-sc",
             "--no-spatial-coarsen", "Enable/disable spatial coarsening.");
   AddOption(&access_level, "-access", "--access-level",
             "Set the access level."); // TODO: what are the options?
   AddOption(&print_level, "-print", "--print-level",
             "Set the print level."); // TODO: what are the options?
   AddOption(&use_seq_soln, "-seq_soln", "--use-sequential-solution",
             "If 1, use the sequential time stepping solution as the initial guess.");
   AddOption(&skip, "-skip", "--skip-work-on-first-down-cycle",
             "If True, skip all work on the first down cycle (0:False, 1:True)");

}

void BraidOptions::AddMeshOptions()
{
   AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   AddOption(&ser_ref_levels, "-rs", "--refine-serial",
             "Number of times to refine the mesh uniformly in serial.");
   AddOption(&par_ref_levels, "-rp", "--refine-parallel",
             "Number of times to refine the mesh uniformly in parallel.");
}

Mesh *BraidOptions::LoadMeshAndSerialRefine()
{
   Mesh *mesh = NULL;
   std::ifstream imesh(mesh_file);
   if (imesh)
   {
      mesh = new Mesh(imesh, 1, 1);

      for (int lev = 0; lev < ser_ref_levels; lev++)
      {
         mesh->UniformRefinement();
      }
   }
   return mesh;
}

ParMesh *BraidOptions::LoadMeshAndRefine(MPI_Comm comm_x)
{
   ParMesh *pmesh = NULL;
   Mesh *mesh = LoadMeshAndSerialRefine();
   if (mesh)
   {
      pmesh = new ParMesh(comm_x, *mesh);
      delete mesh;
      for (int lev = 0; lev < par_ref_levels; lev++)
      {
         pmesh->UniformRefinement();
      }
   }
   return pmesh;
}

void BraidOptions::SetBraidCoreOptions(BraidCore &core)
{
   core.SetAccessLevel(access_level);
   core.SetPrintLevel(print_level);
   core.SetSeqSoln(use_seq_soln);
   core.SetMaxLevels(max_levels);
   core.SetMinCoarse(min_coarse);
   core.SetNRelax(-1, nrelax);
   if (nrelax0 > -1)
      core.SetNRelax(0, nrelax0);
   rtol ? core.SetRelTol(tol) : core.SetAbsTol(tol);
   core.SetCFactor(-1, cfactor);
   core.SetAggCFactor(cfactor0);
   core.SetMaxIter(max_iter);
   core.SetStorage(storage);
   core.SetTemporalNorm(tnorm);
   core.SetSkip(skip);
   if (spatial_coarsen)
   {
      core.SetSpatialCoarsenAndRefine();
   }
   if (nfmg_Vcyc > 0)
   {
      core.SetFMG();
      core.SetNFMGVcyc(nfmg_Vcyc);
   }
}

void SpaceTimeMeshInfo::Reinitialize(int _max_levels)
{
   max_levels = _max_levels;
   mesh_table.SetSize(4*max_levels);
   mesh_table_global.SetSize(4*max_levels);
   mesh_table = -1.0;
   mesh_table_global = -1.0;
}

void SpaceTimeMeshInfo::SetRow(int braid_level, int vec_level, ParMesh * pmesh, double dt)
{
   // Fill in the level-th row with runtime spatial coarsening information
   double h_min, h_max;

   ComputeMeshSize(pmesh, &h_min, &h_max);
   mesh_table[ braid_level*4 ]     = std::max(mesh_table[ braid_level*4 ], (double) vec_level );
   mesh_table[ braid_level*4 + 1 ] = std::max(mesh_table[ braid_level*4 + 1 ], h_min);
   mesh_table[ braid_level*4 + 2 ] = std::max(mesh_table[ braid_level*4 + 2 ], h_max);
   mesh_table[ braid_level*4 + 3 ] = std::max(mesh_table[ braid_level*4 + 3 ], dt);
}

void SpaceTimeMeshInfo::ComputeMeshSize( ParMesh *pmesh, double * h_min_ptr, double * h_max_ptr)
{
    // Compute the maximum and minimum mesh sizes for a pmesh
    int i;
    double h, h_min, h_max;
    int NumOfElements = pmesh->GetNE();
    MPI_Comm MyComm = pmesh->GetComm();

    for (i = 0; i < NumOfElements; i++)
    {
        h = pmesh->GetElementSize(i);
        if (i == 0)
        {
            h_min = h_max = h;
        }
        else
        {
            if (h < h_min)  h_min = h;
            if (h > h_max)  h_max = h;
        }
    }

    double gh_min, gh_max;
    MPI_Reduce(&h_min, &gh_min, 1, MPI_DOUBLE, MPI_MIN, 0, MyComm);
    MPI_Reduce(&h_max, &gh_max, 1, MPI_DOUBLE, MPI_MAX, 0, MyComm);

    *h_max_ptr = gh_max;
    *h_min_ptr = gh_min;
}

// Print the all the rows that are not -1, i.e., that aren't empy
void SpaceTimeMeshInfo::Print(MPI_Comm comm)
{
   int myid, level;
   double h_min, h_max, dt;

   MPI_Comm_rank(comm, &myid);

   // Reduce all the table entries over processors
   for(int i = 0; i < max_levels; i++)
   {
      for(int j = 0; j < 4; j++){
         MPI_Allreduce( &(mesh_table[i*4 + j]), &(mesh_table_global[i*4 + j]), 1, MPI_DOUBLE, MPI_MAX, comm );
      }
   }

   if(myid == 0)
   {
      std::cout << std::endl;
      std::cout << " Space-Time Mesh Information Gathered at Run-Time" << std::endl;
      std::cout << " XBraid level | Spatial level |    h_min    |     h_max   |      dt     |  dt/h_max   | dt/h_max^2 " << std::endl;
      std::cout << " --------------------------------------------------------------------------------------------------" << std::endl;
      for(int i = 0; i < max_levels; i++)
      {
         level = (int) mesh_table_global[i*4];
         h_min =  mesh_table_global[i*4 + 1];
         h_max =  mesh_table_global[i*4 + 2];
         dt =  mesh_table_global[i*4 + 3];
         if (dt < 0)
            break;

         std::cout.precision(4);
         std::cout << std::scientific;
         std::cout << "    " << std::setw(10) << std::left << i << "|"
                   << "    " << std::setw(11) << std::left << level << "|"
                   << " "    << std::setw(12) << std::left << h_min << "|"
                   << " "    << std::setw(12) << std::left << h_max << "|"
                   << " "    << std::setw(12) << std::left << dt << "|"
                   << " "    << std::setw(12) << std::left << dt/h_max << "|"
                   << " "    << std::setw(12) << std::left << dt/(h_max*h_max) << std::endl;
      }
      std::cout << std::endl;
   }
}

#endif // braid_mfem_HEADER
