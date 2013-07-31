/*
   Example 02
   2D diffusion problem

   Interface:    SStructured interface (SStruct)

   Compile with: make drive-02.c

   Sample run:   srun -n 16 diff2d -pgrid 4 4 -nl 32 32 -nt 16

   Description:  This code solves the diffusion problem 

                                  u_t - K * div u = 0

                 in the unit square subject to 
                    zero Dirichlet boundary conditions,u_b = 0, and 
                    initial condition U0 at time t = 0.
	         The domain is split into an p_x x p_y processor grid.
		 Each processor has a n_x x n_y grid, with nodes
                 connected by a 5-point stencil. More precisely, we use
                 central FD in space and backward Euler in time,
                 definining the stencil

                                    -K*dt/(dy^2)
                   -K*dt/(dx^2)  1+2*K*(dt/(dx^2)+dt/(dy^2)) -K*dt/(dx^2)
                                    -K*dt/(dy^2)   
                     
                  We use cell-centered variables, and, therefore, the 
                  nodes are not shared.

                  To incorporate the boundary conditions, we do the 
                  following: Let x_i and x_b be the interior and boundary
                  parts of the solution vector x, respectively, and let 
                  u_b the boundary condition. If we split the matrix A as

                             A = [A_ii A_ib; 
                                  A_bi A_bb],

                  then we solve

                             [A_ii 0; [x_i;    [b_i - A_ib u_b;
                              0    I]  x_b] =        u_b       ].

                  For zero boundary conditions, u_b = 0, we are just 
                  solving

                          A_ii x_i = b_i.

                  Note that this approach is useful for more general types 
                  of boundary conditions.

                  We use a structured solver for the spatial solve at each 
                  time step.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "_hypre_utilities.h"
#include "HYPRE_sstruct_ls.h"

#include "warp.h"

#ifdef M_PI
   #define PI M_PI
#else
   #define PI 3.14159265358979
#endif

/* --------------------------------------------------------------------
 * Structs for my integration routines
 * -------------------------------------------------------------------- */
/* struct my_App contains general information about the problem, its
 * discretizaion, spatial distribution, and solver used at each time point
 *
 *   comm_t      communicator for parallelizing in time
 *   comm_x      communicator for parallelizing in space 
 *   dim_x       spatial dimension
 *   K           diffusion coefficient
 *   nlx         local problem size in x-dimension
 *   nly         local problem size in y-dimension
 *   tstart      initial time
 *   tstop       integration stop time
 *   nt          number of time steps
 *   dx          spatial step size in x-direction
 *   dy          spatial step size in y-direction
 *   dt          time step on finest grid
 *   max_levels  maximum number of time levels
 *   nparts      number of structured parts of the spatial grid 
 *   nvars       number of variable types on structured parts of
 *               the spatial grid
 *   vartypes    (nvars)-dimensional array of variable types on 
 *               structured parts of the spatial grid
 *   grid_x      spatial grid
 *   stencil     discretization stencil object
 *   graph       graph object that determine the non-zero structure
 *               of the discretization matrices
 *   A           array of discretization matrices (one per time level)
 *   dt_A        array of time steps for which discretization matrix
 *               has been created
 *   nA          number of discretization matrices that have been
 *               created
 *   sym         symmetric storage (1) or not (0)
 *   px          number of processors in x-dimension
 *   py          number of processors in y-dimension
 *   pi          x-coordinate of position in processor grid
 *   pj          y-coordinate of position in processor grid
 *   ilower_x    (dim_x)-dimensional array with integer indices of 
 *               local space interval lower bounds
 *   iupper_x    (dim_x)-dimensional array with integer indices of
 *               local space interval upper bounds
 *   object_type object type of vector to access different hypre solvers 
 *   solver      array of solvers used at each time step on different
 *               time levels
 *   n_pre       number of pre-relaxation sweeps in spatial solve
 *   n_post      number of post-relaxation sweeps in spatial solve
 *   rap         coarse-grid operator type for spatial solver
 *                  rap = 0: Galerkin (default)
 *                  rap = 1: non-Galerkin ParFlow operators
 *                  rap = 2: Galerkin, general operators 
 *   relax       relaxation type for spatial solver
 *                  relax = 0: Jacobi
 *                  relax = 1: Weighted Jacobi 
 *                  relax = 2: R/B Gauss-Seidel (default)
 *                  relax = 3: R/B Gauss-Seidel (non-symmetric)
 *   skip        skip levels in spatial PFMG solver (0 or 1)
 */
typedef struct _warp_App_struct {
   MPI_Comm                comm_t;
   MPI_Comm                comm_x;
   int                     dim_x;
   double                  K;
   int                     nlx, nly;
   double                  tstart;
   double                  tstop;
   int                     nt;
   double                  dx, dy;
   double                  dt;
   int                     max_levels;
   int                     nparts; 
   int                     nvars; 
   HYPRE_SStructVariable  *vartypes;
   HYPRE_SStructGrid       grid_x;
   HYPRE_SStructStencil    stencil;
   HYPRE_SStructGraph      graph;
   HYPRE_SStructMatrix    *A;
   double                 *dt_A;
   int                     nA;
   int                     sym;
   int                     px, py;
   int                     pi, pj;
   int                     ilower_x[2], iupper_x[2];
   int                     object_type;
   HYPRE_StructSolver     *solver;
   int                    *max_num_iterations;
   int                     n_pre, n_post;
   int                     rap, relax, skip;
} my_App;

/* struct my_Vector contains local information specific to one time point
 *   x            spatial vector 
 */
typedef struct _warp_Vector_struct
{
   HYPRE_SStructVector   x;
} my_Vector;



int max( int a, int b ){
  return (a >= b ? a : b );
}

/* --------------------------------------------------------------------
 * Initial condition 
 * -------------------------------------------------------------------- */
double U0(double x, double y){
    return sin(x)*sin(y);
}

/* --------------------------------------------------------------------
 * Boundary condition: zero Dirichlet condition for now
 * -------------------------------------------------------------------- */
double B0(double x, double y){
    return 0.0;
}

/* --------------------------------------------------------------------
 * Apply boundary conditions:
 * To incorporate the boundary conditions, we removed the connections 
 * between the interior and boundary nodes in the discretization matrix. 
 * We adjust for removing these connections by appropriately modifying
 * the corresponding RHS entries. 
 * -------------------------------------------------------------------- */
void 
addBoundaryToRHS( HYPRE_SStructVector  b, 
                  double               K, 
                  double               dx, 
                  double               dy, 
                  double               dt,
                  int                  nlx, 
                  int                  nly,
                  int                  px, 
                  int                  py, 
		  int                  pi, 
                  int                  pj, 
                  int                  zeroDirichlet )
{
   int i, j, m;
   
   int bc_ilower[2];
   int bc_iupper[2];
   int istart, iend, jstart, jend;
   double *bvalues;

   /* We have one part and one variable. */
   int part = 0;
   int var = 0;
       
   /* Allocate vector for values on boundary planes */
   bvalues = (double *) malloc( (max(nlx, nly)+1)*sizeof(double) );     
  
       
   /* a) boundaries y = 0 or y = PI */
   /* Processors at y = 0 */
   if (pj == 0){
      bc_ilower[0] = pi*nlx;
      bc_ilower[1] = pj*nly;
           
      bc_iupper[0] = bc_ilower[0] + nlx-1;
      bc_iupper[1] = bc_ilower[1];
           
      /* Put the boundary conditions in b */
      for( i = 0; i < nlx; i++ )
          bvalues[i] = B0( (bc_ilower[0]+i)*dx,
                            bc_ilower[1]*dy );
           
      HYPRE_SStructVectorSetBoxValues(b, part, bc_ilower,
                                      bc_iupper, var, bvalues);
   }
       
   /* Processors at y = PI */
   if (pj == py-1){
      bc_ilower[0] = pi*nlx;
      bc_ilower[1] = pj*nly + nly-1;
           
      bc_iupper[0] = bc_ilower[0] + nlx-1;
      bc_iupper[1] = bc_ilower[1];
           
      /* Put the boundary conditions in b */
      for( i = 0; i < nlx; i++ )
          bvalues[i] = B0( (bc_ilower[0]+i)*dx,
                            bc_ilower[1]*dy );
           
      HYPRE_SStructVectorSetBoxValues(b, part, bc_ilower,
                                      bc_iupper, var, bvalues);
   }
       
   /* b) boundaries x = 0 or x = PI */
   /* Processors at x = 0 */
   if (pi == 0){
      bc_ilower[0] = pi*nlx;
      bc_ilower[1] = pj*nly;
           
      bc_iupper[0] = bc_ilower[0];
      bc_iupper[1] = bc_ilower[1] + nly-1;
          
      /* Put the boundary conditions in b */
      for( j = 0; j < nly; j++ )
          bvalues[j] = B0(  bc_ilower[0]*dx,
                           (bc_ilower[1]+j)*dy );
           
      HYPRE_SStructVectorSetBoxValues(b, part, bc_ilower,
                                      bc_iupper, var, bvalues);
   }
       
   /* Processors at x = PI */
   if (pi == px-1){
      bc_ilower[0] = pi*nlx + nlx-1;
      bc_ilower[1] = pj*nly;
           
      bc_iupper[0] = bc_ilower[0];
      bc_iupper[1] = bc_ilower[1] + nly-1;
          
      /* Put the boundary conditions in b */
      for( j = 0; j < nly; j++ )
          bvalues[j] = B0(  bc_ilower[0]*dx,
                           (bc_ilower[1]+j)*dy );
           
      HYPRE_SStructVectorSetBoxValues(b, part, bc_ilower,
                                      bc_iupper, var, bvalues);
   }

   
   /* For zero Dirichlet boundary conditions, no further changes
    * are necessary. Otherwise, we have to modify the right-hand
    * side vector to adjust for removing connections between
    * the interior and boundary nodes in the discretization
    * matrix. */
   if( zeroDirichlet == 0 ){ 
      /* Neighbors of boundary nodes of boundary y = 0 
       * Neighbors are either
       *   i) on same processor as boundary nodes (pj = 0)
       * or
       *   ii) on neighboring processor (pj = 1) 
       * Case ii) only applies if nly = 1 */
       
      /* Neighbors of boundary on same processor */
      if( (nly > 1) && (pj == 0) )
      {
         bc_ilower[0] = pi*nlx;
         bc_ilower[1] = pj*nly + 1;
           
         bc_iupper[0] = bc_ilower[0] + nlx-1;
         bc_iupper[1] = bc_ilower[1];
           
         istart = 0; iend = nlx;
           
         /* Adjust box to not include boundary nodes */
	 if( pi == 0 ){
	    bc_ilower[0] += 1;
	    istart += 1;
	 }

	 if( pi == px-1 ){
	    bc_iupper[0] -= 1;
	    iend -= 1;
	 }
        
         /* Adjust for removing connections between the boundary
          * and interior nodes in the discretization matrix. */
         for( m = 0, i = istart; i < iend; i++, m++ )
            bvalues[m] = K*(dt/(dy*dy))*
                           B0( (bc_ilower[0]+i-istart)*dx,
                               (bc_ilower[1]-1)*dy );

         HYPRE_SStructVectorAddToBoxValues(b, part, bc_ilower,
                                           bc_iupper, var, bvalues);
      }
      /* Neighbors of boundary on neighboring processor */
      if( (nly == 1) && (pj == 1) )
      {
         bc_ilower[0] = pi*nlx;
         bc_ilower[1] = pj*nly;
           
         bc_iupper[0] = bc_ilower[0] + nlx-1;
         bc_iupper[1] = bc_ilower[1];
           
         istart = 0; iend = nlx;
           
         /* Adjust box to not include boundary nodes */
	 if( pi == 0 ){
	    bc_ilower[0] += 1;
	    istart += 1;
	 }

	 if( pi == px-1 ){
	    bc_iupper[0] -= 1;
	    iend -= 1;
	 }

         /* Adjust for removing connections between the boundary
          * and interior nodes in the discretization matrix. */
         for( m = 0, i = istart; i < iend; i++, m++ )
             bvalues[m] = K*(dt/(dy*dy))*
                           B0( (bc_ilower[0]+i-istart)*dx,
                               (bc_ilower[1]-1)*dy );

         HYPRE_SStructVectorAddToBoxValues(b, part, bc_ilower,
                                           bc_iupper, var, bvalues);
      }

      /* Neighbors of boundary nodes of boundary y = PI
       * Neighbors are either
       *   i) on same processor as boundary nodes (pj = py-1)
       * or
       *   ii) on neighboring processor (pj = py-2) 
       * Case ii) only applies if nly = 1 */
       
      /* Neighbors of boundary on same processor */
      if( (nly > 1) && (pj == py-1) )
      {
         bc_ilower[0] = pi*nlx;
         bc_ilower[1] = pj*nly + nly-1 - 1;
           
         bc_iupper[0] = bc_ilower[0] + nlx-1;
         bc_iupper[1] = bc_ilower[1];
           
         istart = 0; iend = nlx;
           
         /* Adjust box to not include boundary nodes */
         if( pi == 0 ){
   	    bc_ilower[0] += 1;
	    istart += 1;
	 }

	 if( pi == px-1 ){
	    bc_iupper[0] -= 1;
	    iend -= 1;
	 }
        
         /* Adjust for removing connections between the boundary
          * and interior nodes in the discretization matrix. */ 
         for( m = 0, i = istart; i < iend; i++, m++ )
            bvalues[m] = K*(dt/(dy*dy))*
                           B0( (bc_ilower[0]+i-istart)*dx,
                               (bc_ilower[1]+1)*dy );

         HYPRE_SStructVectorAddToBoxValues(b, part, bc_ilower,
                                           bc_iupper, var, bvalues);
      }
      /* Neighbors of boundary on neighboring processor */
      if( (nly == 1) && (pj == py-2) )
      {
         bc_ilower[0] = pi*nlx;
         bc_ilower[1] = pj*nly + nly-1;
           
         bc_iupper[0] = bc_ilower[0] + nlx-1;
         bc_iupper[1] = bc_ilower[1];
           
         istart = 0; iend = nlx;
           
         /* Adjust box to not include boundary nodes */
         if( pi == 0 ){
	    bc_ilower[0] += 1;
	    istart += 1;
	 }

	 if( pi == px-1 ){
	    bc_iupper[0] -= 1;
	    iend -= 1;
	 }
           
         /* Adjust for removing connections between the boundary
          * and interior nodes in the discretization matrix. */
         for( m = 0, i = istart; i < iend; i++, m++ )
            bvalues[m] = K*(dt/(dy*dy))*
                           B0( (bc_ilower[0]+i-istart)*dx,
                               (bc_ilower[1]+1)*dy );

         HYPRE_SStructVectorAddToBoxValues(b, part, bc_ilower,
                                           bc_iupper, var, bvalues);
      }

      /* Neighbors of boundary nodes of boundary x = 0 
       * Neighbors are either
       *   i) on same processor as boundary nodes (pi = 0)
       * or
       *   ii) on neighboring processor (pi = 1) 
       * Case ii) only applies if nlx = 1 */
       
      /* Neighbors of boundary on same processor */
      if( (nlx > 1) && (pi == 0) )
      {
         bc_ilower[0] = pi*nlx + 1;
         bc_ilower[1] = pj*nly;
           
         bc_iupper[0] = bc_ilower[0];
         bc_iupper[1] = bc_ilower[1] + nly-1;
           
         jstart = 0; jend = nly;
           
         /* Adjust box to not include boundary nodes */
         if( pj == 0 ){
	    bc_ilower[1] += 1;
	    jstart += 1;
	 }

	 if( pj == py-1 ){
	    bc_iupper[1] -= 1;
	    jend -= 1;
	 }          

         /* Adjust for removing connections between the boundary
          * and interior nodes in the discretization matrix. */
         for( m = 0, j = jstart; j < jend; j++, m++ )
            bvalues[m] = K*(dt/(dx*dx))*
                           B0( (bc_ilower[0]-1)*dx,
                               (bc_ilower[1]+j-jstart)*dy );

         HYPRE_SStructVectorAddToBoxValues(b, part, bc_ilower,
                                           bc_iupper, var, bvalues);
      }
      /* Neighbors of boundary on neighboring processor */
      if( (nlx == 1) && (pi == 1) )
      {
         bc_ilower[0] = pi*nlx;
         bc_ilower[1] = pj*nly;
          
         bc_iupper[0] = bc_ilower[0];
         bc_iupper[1] = bc_ilower[1] + nly-1;
           
         jstart = 0; jend = nly;
           
         /* Adjust box to not include boundary nodes */
         if( pj == 0 ){
	    bc_ilower[1] += 1;
	    jstart += 1;
	 }

	 if( pj == py-1 ){
	    bc_iupper[1] -= 1;
	    jend -= 1;
	 }       

         /* Adjust for removing connections between the boundary
          * and interior nodes in the discretization matrix. */
         for( m = 0, j = jstart; j < jend; j++, m++ )
            bvalues[m] = K*(dt/(dx*dx))*
                           B0( (bc_ilower[0]-1)*dx,
                               (bc_ilower[1]+j-jstart)*dy );

         HYPRE_SStructVectorAddToBoxValues(b, part, bc_ilower,
                                           bc_iupper, var, bvalues);
      }

      /* Neighbors of boundary nodes of boundary x = PI
       * Neighbors are either
       *   i) on same processor as boundary nodes (pi = px-1)
       * or
       *   ii) on neighboring processor (pi = px-2) 
       * Case ii) only applies if nlx = 1 */
       
      /* Neighbors of boundary on same processor */
      if( (nlx > 1) && (pi == px-1) )
      {
         bc_ilower[0] = pi*nlx + nlx-1 - 1;
         bc_ilower[1] = pj*nly;
           
         bc_iupper[0] = bc_ilower[0];
         bc_iupper[1] = bc_ilower[1] + nly-1;

         jstart = 0; jend = nly;
           
         /* Adjust box to not include boundary nodes */
         if( pj == 0 ){
	    bc_ilower[1] += 1;
	    jstart += 1;
	 }

	 if( pj == py-1 ){
	    bc_iupper[1] -= 1;
	    jend -= 1;
	 }

         /* Adjust for removing connections between the boundary
          * and interior nodes in the discretization matrix. */
         for( m = 0, j = jstart; j < jend; j++, m++ )
            bvalues[m] = K*(dt/(dx*dx))*
                           B0( (bc_ilower[0]+1)*dx,
                               (bc_ilower[1]+j-jstart)*dy );

         HYPRE_SStructVectorAddToBoxValues(b, part, bc_ilower,
                                           bc_iupper, var, bvalues);
      }
      /* Neighbors of boundary on neighboring processor */
      if( (nlx == 1) && (pi == px-2) )
      {
         bc_ilower[0] = pi*nlx + nlx-1;
         bc_ilower[1] = pj*nly;
           
         bc_iupper[0] = bc_ilower[0];
         bc_iupper[1] = bc_ilower[1] + nly-1;

         jstart = 0; jend = nly;
           
         /* Adjust box to not include boundary nodes */
         if( pj == 0 ){
            bc_ilower[1] += 1;
	    jstart += 1;
	 }

	 if( pj == py-1 ){
	    bc_iupper[1] -= 1;
	    jend -= 1;
	 }
        
         /* Adjust for removing connections between the boundary
          * and interior nodes in the discretization matrix. */
         for( m = 0, j = jstart; j < jend; j++, m++ )
            bvalues[m] = K*(dt/(dx*dx))*
                           B0( (bc_ilower[0]+1)*dx,
                               (bc_ilower[1]+j-jstart)*dy );

         HYPRE_SStructVectorAddToBoxValues(b, part, bc_ilower,
                                           bc_iupper, var, bvalues);
      }
   } /* end if zeroDirichlet == 0 */
   free(bvalues);

   /* Finalize the vector assembly. */
   HYPRE_SStructVectorAssemble(b);
}


/* --------------------------------------------------------------------
 * Set up a 2D grid.
 * -------------------------------------------------------------------- */
void
setUp2Dgrid( MPI_Comm               comm,
             HYPRE_SStructGrid     *grid_ptr,
             int                    ndim,
             int                    nparts,
             int                   *ilower, 
             int                   *iupper,
             HYPRE_SStructVariable *vartypes )
{
   /* We have one part and one variable. */
   int nvars = 1;
   int part = 0;

   HYPRE_SStructGrid grid;

   /* Create an empty 2D grid object. */
   HYPRE_SStructGridCreate( comm, ndim, nparts, &grid );

   /* Add a new box to the grid. */
   HYPRE_SStructGridSetExtents( grid, part, ilower, iupper );

   /* Set the variable type for each part */
   for( part = 0; part < nparts; part++ )
      HYPRE_SStructGridSetVariables( grid, part, nvars, vartypes );

   /* This is a collective call finalizing the grid assembly.
    * The grid is now ``ready to be used''. */
   HYPRE_SStructGridAssemble( grid );

   *grid_ptr = grid;
}


/* --------------------------------------------------------------------
 * Define the discretization stencil.
 * -------------------------------------------------------------------- */
void
set5ptStencil( HYPRE_SStructStencil *stencil_ptr,
               int                   sym,
               int                   ndim )
{
   /* We have one variable. */
   int var = 0;

   int i;

   HYPRE_SStructStencil stencil;

   if (sym == 0)
   {
      /* Create an empty 2D, 5-pt stencil object. */
      HYPRE_SStructStencilCreate( ndim, 5, &stencil );

      /* Define the geometry of the stencil. */
      int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};

      /* Assign stencil entries. */
      for( i = 0; i < 5; i++ )
         HYPRE_SStructStencilSetEntry( stencil, i, offsets[i], var );
   }
   else /* Symmetric storage */
   {
      /* Create an empty 2D, 3-pt stencil object. */
      HYPRE_SStructStencilCreate( ndim, 3, &stencil );
 
      /* Define the geometry of the stencil. */
      int offsets[3][2] = {{0,0}, {1,0}, {0,1}};

      /* Assign stencil entries */
      for( i = 0; i < 3; i++ )
         HYPRE_SStructStencilSetEntry( stencil, i, offsets[i], var );
   }

   *stencil_ptr = stencil;
}


/* --------------------------------------------------------------------
 * Set up the graph -- this defines the non-zero structure of the
 * discretization matrix.
 * -------------------------------------------------------------------- */
void
setUpGraph( MPI_Comm               comm,
            HYPRE_SStructGraph    *graph_ptr,
            HYPRE_SStructGrid      grid,
            int                    object_type,
            HYPRE_SStructStencil   stencil )
{
   int var = 0;
   int part = 0;

   HYPRE_SStructGraph graph;

   /* Create the graph object. */
   HYPRE_SStructGraphCreate( comm, grid, &graph );

   /* Set the object type. */
   HYPRE_SStructGraphSetObjectType( graph, object_type );

   /* Now we need to tell the graph which stencil to use for each
    * variable on each part (we only have one variable and one part) */
   HYPRE_SStructGraphSetStencil( graph, part, var, stencil );

   /* Here we could establish connections between parts if we
    * had more than one part. */

   /* Assemble the graph. */
   HYPRE_SStructGraphAssemble( graph );

   *graph_ptr = graph;
}


/* --------------------------------------------------------------------
 * Set up the discretization matrix. 
 * First, we set the stencil values at every node neglecting the 
 * boundary. Then, we correct the matrix stencil at boundary nodes.
 * We have to eliminate the coefficients reaching outside of the domain
 * boundary. Furthermore, to incorporate boundary conditions, we remove
 * the connections between the interior nodes and boundary nodes.
 * -------------------------------------------------------------------- */
void
setUpMatrix( MPI_Comm             comm,
             HYPRE_SStructMatrix *A_ptr,
             HYPRE_SStructGraph   graph, 
             int                  sym,
             int                  object_type,
             double               K, 
             double               dx, 
             double               dy, 
             double               dt,
             int                 *ilower, 
             int                 *iupper,
             int                  nlx, 
             int                  nly,
             int                  px, 
             int                  py,
             int                  pi, 
             int                  pj )   
{
   int i, j, m, idx;

   /* We have one part and one variable. */
   int part = 0;
   int var = 0;

   int stencil_indices[5] = {0, 1, 2, 3, 4};
   int nentries;

   double *values, *bvalues;
   int bc_ilower[2];
   int bc_iupper[2];

   HYPRE_SStructMatrix A;

   /* Create an empty matrix object. */
   HYPRE_SStructMatrixCreate( comm, graph, &A);

   /* Use symmetric storage? The function below is for symmetric stencil
    * entries (use HYPRE_SStructMatrixSetNSSymmetric for non-stencil 
    * entries). */
   HYPRE_SStructMatrixSetSymmetric( A, part, var, var, sym );

   /* Set the object type. */
   HYPRE_SStructMatrixSetObjectType( A, object_type );

   /* Indicate that the matrix coefficients are ready to be set. */
   HYPRE_SStructMatrixInitialize( A );

   /* 1. neglect boundaries */
   /* Set the stencil values. */
   if (sym == 0)
   {
      nentries = 5;
      values = (double *) malloc( nentries*(nlx*nly)*sizeof(double) );

      /* The order is left-to-right, bottom-to-top */
      for( m = 0, j = 0; j < nly; j++ )
         for( i = 0; i < nlx; i++, m+=nentries ){
            for( idx = 1; idx <= 2; idx++ )
               values[m+idx] = -K*(dt/(dx*dx));
            for( idx = 3; idx <= 4; idx++ )
               values[m+idx] = -K*(dt/(dy*dy));
                     
            values[m] = 1.0 + 2*K*( (dt/(dx*dx)) + (dt/(dy*dy)) );
         }

      HYPRE_SStructMatrixSetBoxValues( A, part, ilower, iupper, var, 
                                       nentries, stencil_indices, 
                                       values );

      free(values);
   }
   else /* Symmetric storage */
   {
      nentries = 3;
      values = (double *) malloc( nentries*(nlx*nly)*sizeof(double) );

      /* The order is left-to-right, bottom-to-top */
      for( m = 0, j = 0; j < nly; j++ )
         for( i = 0; i < nlx; i++, m += nentries ){
            values[m+1] = -K*(dt/(dx*dx));
            values[m+2] = -K*(dt/(dy*dy));  
            values[m]   = 1.0 + 2*K*( (dt/(dx*dx)) + (dt/(dy*dy)) );
         }

      HYPRE_SStructMatrixSetBoxValues( A, part, ilower, iupper, var,
                                       nentries, stencil_indices, 
                                       values );

      free(values);
   }

   /* 2. correct stencils at boundary nodes */     
   /* Allocate vectors for values on boundary planes */
   values  = (double *) malloc( nentries*( max(nlx,nly)+1 )*
                                         sizeof(double) );
   bvalues = (double *) malloc( (max(nlx,nly)+1)*sizeof(double) );     
       
   /* a) boundaries y = 0 or y = PI */
   /* The stencil at the boundary nodes is 1-0-0-0-0. Because
    * we have I x_b = u_b. */
   for( i = 0; i < nentries*nlx; i+= nentries ){
      values[i] = 1.0;
      for( idx = 1; idx < nentries; idx++ )
         values[i+idx] = 0.0;
   }
       
   /* Processors at y = 0 */
   if( pj == 0 ){
      bc_ilower[0] = pi*nlx;
      bc_ilower[1] = pj*nly;
           
      bc_iupper[0] = bc_ilower[0] + nlx-1;
      bc_iupper[1] = bc_ilower[1];
           
      /* Modify the matrix. */
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, nentries,
                                      stencil_indices, values);
   }
       
   /* Processors at y = PI */
   if( pj == py-1 ){
      bc_ilower[0] = pi*nlx;
      bc_ilower[1] = pj*nly + nly-1;
           
      bc_iupper[0] = bc_ilower[0] + nlx-1;
      bc_iupper[1] = bc_ilower[1];
           
      /* Modify the matrix */
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, nentries,
                                      stencil_indices, values);
   }
       
   /* b) boundaries x = 0 or x = PI */    
   /* The stencil at the boundary nodes is 1-0-0-0-0. Because
    * we have I x_b = u_b. */
   for( j = 0; j < nentries*nly; j+= nentries ){
      values[j] = 1.0;
      for( idx = 1; idx < nentries; idx++ )
         values[j+idx] = 0.0;
   }
       
   /* Processors at x = 0 */
   if( pi == 0 ){
      bc_ilower[0] = pi*nlx;
      bc_ilower[1] = pj*nly;
           
      bc_iupper[0] = bc_ilower[0];
      bc_iupper[1] = bc_ilower[1] + nly-1;
           
      /* Modify the matrix */
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, nentries,
                                      stencil_indices, values);
   }
       
   /* Processors at x = PI */
   if( pi == px-1 ){
      bc_ilower[0] = pi*nlx + nlx-1;
      bc_ilower[1] = pj*nly;
           
      bc_iupper[0] = bc_ilower[0];
      bc_iupper[1] = bc_ilower[1] + nly-1;
           
      /* Modify the matrix */
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, nentries,
                                      stencil_indices, values);
   }

   free(values);    
       
   /* Recall that the system we are solving is:
    *
    *   [A_ii 0; [x_i;    [b_i - A_ib*u_b;
    *      0  I]  x_b ] =        u_b       ].
    * 
    * This requires removing the connections between the interior
    * and boundary nodes that we have set up when we set the
    * 5pt stencil at each node. For the symmetric ordering 
    * scheme, we just do the top and right boundary. */
       
   /* a) Neighbors of boundary nodes of boundary y = 0.
    *    Neighbors are either
    *      i) on same processor as boundary nodes (pj = 0)
    *    or
    *      ii) on neighboring processor (pj = 1) 
    *    Case ii) only applies if nly = 1 */
       
   /* Neighbors of boundary on same processor. */
   if( (nly > 1) && (pj == 0) && (sym == 0) )
   {
      bc_ilower[0] = pi*nlx;
      bc_ilower[1] = pj*nly + 1;
           
      bc_iupper[0] = bc_ilower[0] + nlx-1;
      bc_iupper[1] = bc_ilower[1];
           
      /* Adjust box to not include boundary nodes */
      if( pi == 0 )
         bc_ilower[0] += 1;

      if( pi == px-1 )
         bc_iupper[0] -= 1;
           
      stencil_indices[0] = 3;
           
      /* Modify the matrix */
      for( m = 0; m < nlx; m++ )
         bvalues[m] = 0.0;
           
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, 1,
                                      stencil_indices, bvalues);
   }
   /* Neighbors of boundary on neighboring processor. */
   if( (nly == 1) && (pj == 1) && (sym == 0) )
   {
      bc_ilower[0] = pi*nlx;
      bc_ilower[1] = pj*nly;
           
      bc_iupper[0] = bc_ilower[0] + nlx-1;
      bc_iupper[1] = bc_ilower[1];
           
      /* Adjust box to not include boundary nodes */
      if( pi == 0 )
         bc_ilower[0] += 1;

      if( pi == px-1 )
         bc_iupper[0] -= 1;
           
      stencil_indices[0] = 3;
           
      /* Modify the matrix */
      for( m = 0; m < nlx; m++ )
          bvalues[m] = 0.0;
           
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, 1,
                                      stencil_indices, bvalues);
   }

   /* b) Neighbors of boundary nodes of boundary y = PI.
    *    Neighbors are either
    *      i) on same processor as boundary nodes (pj = py-1)
    *    or
    *      ii) on neighboring processor (pj = py-2) 
    *    Case ii) only applies if nly = 1 */  
   /* Neighbors of boundary on same processor */
   if( (nly > 1) && (pj == py-1) )
   {
      bc_ilower[0] = pi*nlx;
      bc_ilower[1] = pj*nly + nly-1 - 1;
           
      bc_iupper[0] = bc_ilower[0] + nlx-1;
      bc_iupper[1] = bc_ilower[1];
           
           
      /* Adjust box to not include boundary nodes */
      if( pi == 0 )
         bc_ilower[0] += 1;

      if( pi == px-1 )
         bc_iupper[0] -= 1;
           
      if( sym == 0 )
         stencil_indices[0] = 4;
      else
         stencil_indices[0] = 2;
           
      /* Modify the matrix */
      for( m = 0; m < nlx; m++ )
         bvalues[m] = 0.0;
           
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, 1, stencil_indices, bvalues);
   }
   /* Neighbors of boundary on neighboring processor */
   if( (nly == 1) && (pj == py-2) )
   {
      bc_ilower[0] = pi*nlx;
      bc_ilower[1] = pj*nly + nly-1;
           
      bc_iupper[0] = bc_ilower[0] + nlx-1;
      bc_iupper[1] = bc_ilower[1];
           
      /* Adjust box to not include boundary nodes */
      if( pi == 0 )
         bc_ilower[0] += 1;

      if( pi == px-1 )
         bc_iupper[0] -= 1;
           
      if( sym == 0 )
         stencil_indices[0] = 4;
      else
         stencil_indices[0] = 2;
           
      /* Modify the matrix */
      for( m = 0; m < nlx; m++ )
          bvalues[m] = 0.0;
           
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, 1, stencil_indices, bvalues);
   }

   /* c) Neighbors of boundary nodes of boundary x = 0. 
    *    Neighbors are either
    *      i) on same processor as boundary nodes (pi = 0)
    *    or
    *      ii) on neighboring processor (pi = 1) 
    *    Case ii) only applies if nlx = 1 */  
   /* Neighbors of boundary on same processor */
   if( (nlx > 1) && (pi == 0) && (sym == 0) )
   {
      bc_ilower[0] = pi*nlx + 1;
      bc_ilower[1] = pj*nly;
           
      bc_iupper[0] = bc_ilower[0];
      bc_iupper[1] = bc_ilower[1] + nly-1;
           
      /* Adjust box to not include boundary nodes */
      if( pj == 0 )
         bc_ilower[1] += 1;

      if( pj == py-1 )
         bc_iupper[1] -= 1;

      stencil_indices[0] = 1;
           
      /* Modify the matrix */
      for( m = 0; m < nly; m++ )
         bvalues[m] = 0.0;
           
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, 1,
                                      stencil_indices, bvalues);
   }
   /* Neighbors of boundary on neighboring processor */
   if( (nlx == 1) && (pi == 1) && (sym == 0) )
   {
      bc_ilower[0] = pi*nlx;
      bc_ilower[1] = pj*nly;
           
      bc_iupper[0] = bc_ilower[0];
      bc_iupper[1] = bc_ilower[1] + nly-1;
           
      /* Adjust box to not include boundary nodes */
      if( pj == 0 )
         bc_ilower[1] += 1;

      if( pj == py-1 )
         bc_iupper[1] -= 1;

      stencil_indices[0] = 1;
           
      /* Modify the matrix */
      for( m = 0; m < nly; m++ )
         bvalues[m] = 0.0;
           
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, 1,
                                      stencil_indices, bvalues);
   }

   /* d) Neighbors of boundary nodes of boundary x = PI.
    *    Neighbors are either
    *      i) on same processor as boundary nodes (pi = px-1)
    *    or
    *      ii) on neighboring processor (pi = px-2) 
    *    Case ii) only applies if nlx = 1 */  
   /* Neighbors of boundary on same processor */
   if( (nlx > 1) && (pi == px-1) )
   {
      bc_ilower[0] = pi*nlx + nlx-1 - 1;
      bc_ilower[1] = pj*nly;
           
      bc_iupper[0] = bc_ilower[0];
      bc_iupper[1] = bc_ilower[1] + nly-1;
           
      /* Adjust box to not include boundary nodes */
      if( pj == 0 )
         bc_ilower[1] += 1;

      if( pj == py-1 )
         bc_iupper[1] -= 1;
           
      if( sym == 0 )
         stencil_indices[0] = 2;
      else
         stencil_indices[0] = 1;
           
      /* Modify the matrix */
      for( m = 0; m < nly; m++ )
         bvalues[m] = 0.0;
           
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, 1, stencil_indices, bvalues);
   }
   /* Neighbors of boundary on neighboring processor */
   if( (nlx == 1) && (pi == px-2) )
   {
      bc_ilower[0] = pi*nlx + nlx-1;
      bc_ilower[1] = pj*nly;
           
      bc_iupper[0] = bc_ilower[0];
      bc_iupper[1] = bc_ilower[1] + nly-1;
           
      /* Adjust box to not include boundary nodes */
      if( pj == 0 )
         bc_ilower[1] += 1;

      if( pj == py-1 )
         bc_iupper[1] -= 1;
           
      if( sym == 0 )
         stencil_indices[0] = 2;
      else
         stencil_indices[0] = 1;
           
      /* Modify the matrix */
      for( m = 0; m < nly; m++ )
         bvalues[m] = 0.0;
           
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, 1, stencil_indices, bvalues);
   }

   free(bvalues);

   /* Finalize the matrix assembly. */
   HYPRE_SStructMatrixAssemble( A );

   *A_ptr = A;
}


/* --------------------------------------------------------------------
 * Set up PFMG solver.
 * -------------------------------------------------------------------- */
void
setUpStructSolver( MPI_Comm             comm,
                   HYPRE_StructSolver  *solver_ptr,
                   HYPRE_StructMatrix  *sA_ptr,
                   HYPRE_StructVector  *sb_ptr,
                   HYPRE_StructVector  *sx_ptr,
                   HYPRE_SStructMatrix  A,
                   HYPRE_SStructVector  b,
                   HYPRE_SStructVector  x,
                   int                  rap,
                   int                  relax,
                   int                  n_pre,
                   int                  n_post,
                   int                  skip )  
{
   HYPRE_StructSolver  solver;
   HYPRE_StructMatrix  sA;
   HYPRE_StructVector  sb;
   HYPRE_StructVector  sx;

   HYPRE_SStructMatrixGetObject( A, (void **) &sA );
   HYPRE_SStructVectorGetObject( b, (void **) &sb );
   HYPRE_SStructVectorGetObject( x, (void **) &sx );

   /* Set PFMG options. */
   HYPRE_StructPFMGCreate( comm, &solver );
   HYPRE_StructPFMGSetMaxIter( solver, 50 );
   HYPRE_StructPFMGSetTol( solver, 1.0e-09 );
   HYPRE_StructPFMGSetRelChange( solver, 0 );
   HYPRE_StructPFMGSetRAPType( solver, rap );
   HYPRE_StructPFMGSetRelaxType( solver, relax );
   HYPRE_StructPFMGSetNumPreRelax( solver, n_pre );
   HYPRE_StructPFMGSetNumPostRelax( solver, n_post );
   HYPRE_StructPFMGSetSkipRelax( solver, skip );
   HYPRE_StructPFMGSetPrintLevel( solver, 1 );
   HYPRE_StructPFMGSetLogging( solver, 1 );

   /* Set up PFMG solver. */
   HYPRE_StructPFMGSetup( solver, sA, sb, sx );

   *solver_ptr = solver;
   *sA_ptr     = sA;
   *sb_ptr     = sb;
   *sx_ptr     = sx;
}


/* ====================================================================
 * My integration routines.
 * ==================================================================== */

/* --------------------------------------------------------------------
 * Time integrator routine.
 * This routine performs the update
 *   u_i = Phi_i(u_{i-1})        (if gzero = 0)
 *   u_i = Phi_i(u_{i-1}) + g_i  (if gzero = 1)
 * Note that the first case corresponds to assuming zero Dirichlet BCs
 * and a zero RHS of the PDE.
 * When Phi is called, u is u_{i-1}. At the end of the routine, u is 
 * set to u_i.
 * -------------------------------------------------------------------- */
int
my_Phi(warp_App     app,
       double       tstart,
       double       tstop,
       warp_Int     gzero,
       warp_Vector  u,
       int         *rfactor_ptr)
{
   int i, A_idx;
   double *values;

   HYPRE_SStructVector b;
   HYPRE_StructMatrix  sA;
   HYPRE_StructVector  sb;
   HYPRE_StructVector  sx;
   
   /* We have one part and one variable. */
   int part = 0;
   int var = 0;

   int num_iterations;

   /* -----------------------------------------------------------------
    * Set up the discretization matrix.
    * If no variable coefficients, check matrix lookup table if matrix 
    * has already been created for time step size tstop-tstart.
    * ----------------------------------------------------------------- */
   A_idx = -1.0;
   for( i = 0; i < app->max_levels; i++ ){
      if( app->dt_A[i] == -1.0 )
         break;
      if( fabs( app->dt_A[i] - (tstop-tstart) ) <= (app->dt)/10 ){
         A_idx = i;
         break;
      }
   }

   if( A_idx == -1.0 ){
      A_idx = i;
      app->nA++;
/*   printf( "Create new matrix %d\n", A_idx );*/
      /* No matrix for time step tstop-tstart exists. 
       * We can add if statement here to compare to spatial step size
       * and use explicit or implicit scheme. */   
      /* Add entry to matrix lookup table. */
      app->dt_A[A_idx] = tstop-tstart;

      /* Set up the discretization matrix. */
      setUpMatrix( app->comm_x, &(app->A[A_idx]), app->graph, 
                   app->sym, app->object_type, app->K, app->dx, app->dy,
                   app->dt_A[A_idx], app->ilower_x, app->iupper_x, 
                   app->nlx, app->nly, app->px, app->py, app->pi,
                   app->pj );

      
      /* Set up the PFMG solver using u->x as dummy vectors. */
      setUpStructSolver( app->comm_x, &(app->solver[A_idx]), &sA, &sb, 
                         &sx, app->A[A_idx], u->x, u->x, app->rap,
                         app->relax, app->n_pre, app->n_post, 
                         app->skip );
   } 

   /* -----------------------------------------------------------------
    * Set up the right-hand side vector. 
    * The right-hand side is the solution from the previous time step
    * modified to incorporate the boundary conditions and the right-
    * hand side of the PDE (here, we assume a zero RHS of the PDE).
    * ----------------------------------------------------------------- */
   values = (double *) malloc( (app->nlx)*(app->nly)*sizeof(double) );
   HYPRE_SStructVectorGather( u->x );
   HYPRE_SStructVectorGetBoxValues( u->x, part, app->ilower_x,
                                    app->iupper_x, var, values );
   HYPRE_SStructVectorCreate( app->comm_x, app->grid_x, &b );
   HYPRE_SStructVectorSetObjectType( b, app->object_type );
   HYPRE_SStructVectorInitialize( b );
   HYPRE_SStructVectorSetBoxValues( b, part, app->ilower_x,
                                    app->iupper_x, var, values );
   free( values );
   addBoundaryToRHS( b, app->K, app->dx, app->dy, tstop-tstart,
                     app->nlx, app->nly, app->px, app->py, 
                     app->pi, app->pj, gzero );
   /* add infos from RHS of PDE here */ 

   /* -----------------------------------------------------------------
    * Time integration to next time point: Solve the system Ax = b.
    * ----------------------------------------------------------------- */
   HYPRE_SStructMatrixGetObject( app->A[A_idx], (void **) &sA );
   HYPRE_SStructVectorGetObject( b, (void **) &sb );
   HYPRE_SStructVectorGetObject( u->x, (void **) &sx );
   
   HYPRE_StructPFMGSolve( app->solver[A_idx], sA, sb, sx );
   HYPRE_StructPFMGGetNumIterations( app->solver[A_idx], &num_iterations );
   app->max_num_iterations[A_idx] = max((app->max_num_iterations[A_idx]),
                                        num_iterations);

   /* free memory */
   HYPRE_SStructVectorDestroy( b );

   /* no refinement */
   *rfactor_ptr = 1;

   return 0;
}


/* --------------------------------------------------------------------
 * Create a vector object for a given time point.
 * -------------------------------------------------------------------- */
int
my_Init(warp_App     app,
        double       t,
        warp_Vector *u_ptr)
{
   my_Vector *u;
   double    *values;
   int        i, j, m;
   int        part;
   int        var;

   u = (my_Vector *) malloc( sizeof(my_Vector) );

   /* Create an empty vector object. */
   HYPRE_SStructVectorCreate( app->comm_x, app->grid_x, &(u->x) );
   
   /* Set the object type (by default HYPRE_SSTRUCT). */
   HYPRE_SStructVectorSetObjectType( u->x, app->object_type );

   /* Indicate that the vector coefficients are ready to be set. */
   HYPRE_SStructVectorInitialize( u->x );

   /* Set the values in left-to-right, bottom-to-top order. 
    * Use initial condition at start time and a random initial guess
    * for all other times. */
   values = (double *) malloc( (app->nlx)*(app->nly)*sizeof(double) );
   if( t == app->tstart )
   {
      /* Initial guess */
      for( m = 0, j = 0; j < app->nly; j++ )
         for (i = 0; i < app->nlx; i++, m++)
            values[m] = U0( (app->ilower_x[0]+i)*(app->dx), 
                            (app->ilower_x[1]+j)*(app->dy) );
   }
   else
   {
      /* Random between 0 and 1 */
      for( m = 0; m < ((app->nlx)*(app->nly)); m++ )
         /* Random between 0 and 1 */
         values[m] = rand()/RAND_MAX;
   }
   for( part = 0; part < app->nparts; part++ )
      for( var = 0; var < app->nvars; var++ )
         HYPRE_SStructVectorSetBoxValues( u->x, part, app->ilower_x,
                                          app->iupper_x, var, values ); 

   HYPRE_SStructVectorAssemble( u->x );
   free(values);

   *u_ptr = u;

   return 0;
}


/* --------------------------------------------------------------------
 * Create a a copy of a vector object.
 * -------------------------------------------------------------------- */
int
my_Clone(warp_App     app,
         warp_Vector  u,
         warp_Vector *v_ptr)
{
   my_Vector *v;
   double    *values;
   int        part;
   int        var;

   v = (my_Vector *) malloc(sizeof(my_Vector));

   /* Create an empty vector object. */
   HYPRE_SStructVectorCreate( app->comm_x, app->grid_x, &(v->x) );
   
   /* Set the object type (by default HYPRE_SSTRUCT) and initialize. */
   HYPRE_SStructVectorSetObjectType( v->x, app->object_type );
   HYPRE_SStructVectorInitialize( v->x );

   /* Set the values. */
   values = (double *) malloc( (app->nlx)*(app->nly)*sizeof(double) );
   HYPRE_SStructVectorGather( u->x );
   for( part = 0; part < app->nparts; part++ )
      for( var = 0; var < app->nvars; var++ ){
         HYPRE_SStructVectorGetBoxValues( u->x, part, app->ilower_x,
                                          app->iupper_x, var, values );
         HYPRE_SStructVectorSetBoxValues( v->x, part, app->ilower_x,
                                          app->iupper_x, var, values );
      }
   free( values );
   HYPRE_SStructVectorAssemble( v->x );

   *v_ptr = v;

   return 0;
}


/* --------------------------------------------------------------------
 * Destroy vector object.
 * -------------------------------------------------------------------- */
int
my_Free(warp_App    app,
        warp_Vector u)
{
   HYPRE_SStructVectorDestroy( u->x );
   free( u );

   return 0;
}


/* --------------------------------------------------------------------
 * Compute vector sum y = alpha*x + beta*y.
 * -------------------------------------------------------------------- */
int
my_Sum(warp_App    app,
       double      alpha,
       warp_Vector x,
       double      beta,
       warp_Vector y)
{
   int i;
   double *values_x, *values_y;
   int part;
   int var;
   
   values_x = (double *) malloc( (app->nlx)*(app->nly)*sizeof(double) );
   values_y = (double *) malloc( (app->nlx)*(app->nly)*sizeof(double) );

   for( part = 0; part < app->nparts; part++ )
      for( var = 0; var < app->nvars; var++ ){
         HYPRE_SStructVectorGather( x->x );
         HYPRE_SStructVectorGetBoxValues( x->x, part, app->ilower_x,
                                          app->iupper_x, var, values_x );

         HYPRE_SStructVectorGather( y->x );
         HYPRE_SStructVectorGetBoxValues( y->x, part, app->ilower_x,
                                          app->iupper_x, var, values_y );

         for( i = 0; i < (app->nlx)*(app->nly); i++ )
            values_y[i] = alpha*values_x[i] + beta*values_y[i];

         HYPRE_SStructVectorSetBoxValues( y->x, part, app->ilower_x,
                                          app->iupper_x, var, values_y );
      }

   free( values_x );
   free( values_y );

   return 0;
}


/* --------------------------------------------------------------------
 * Compute dot product.
 * -------------------------------------------------------------------- */
int
my_Dot(warp_App     app,
       warp_Vector  u,
       warp_Vector  v,
       double      *dot_ptr)
{
   double dot;
   int i;
   double *values_u, *values_v;
   int part;
   int var;
   
   values_u = (double *) malloc( (app->nlx)*(app->nly)*sizeof(double) );
   values_v = (double *) malloc( (app->nlx)*(app->nly)*sizeof(double) );

   for( part = 0; part < app->nparts; part++ )
      for( var = 0; var < app->nvars; var++ ){
         HYPRE_SStructVectorGather( u->x );
         HYPRE_SStructVectorGetBoxValues( u->x, part, app->ilower_x,
                                          app->iupper_x, var, values_u );

         HYPRE_SStructVectorGather( v->x );
         HYPRE_SStructVectorGetBoxValues( v->x, part, app->ilower_x,
                                          app->iupper_x, var, values_v );

         dot = 0;
         for( i = 0; i < (app->nlx)*(app->nly); i++ )
            dot += values_u[i]*values_v[i];
      }
   
   *dot_ptr = dot;

   free( values_u );
   free( values_v );

   return 0;
}


/* --------------------------------------------------------------------
 * Write the vector.
 * -------------------------------------------------------------------- */
int
my_Write(warp_App     app,
         warp_Float   t,
         warp_Vector  u)
{
   MPI_Comm   comm   = MPI_COMM_WORLD;
   double     tstart = (app->tstart);
   double     tstop  = (app->tstop);
   int        ntime  = (app->nt);
   int        index, myid;
   char       filename[255];
   FILE      *file;

   double    *values;
   /* We have one part and one variable. */
   int part = 0;
   int var = 0;

   /*  damping factor */
   double damping = 1.0 / ( 1 + ((2*(app->K)*(app->dt))/
                                 ((app->dx)*(app->dx)))*(1-cos(app->dx)) 
                              + ((2*(app->K)*(app->dt))/
                                 ((app->dy)*(app->dy)))*(1-cos(app->dy)));
   double damping_nt = 1.0; 

   int i, j, m;
   int  pi  = (app->pi);
   int  pj  = (app->pj);
   int  px  = (app->px);
   int  nlx = (app->nlx);
   int  nly = (app->nly);

   index = ((t-tstart) / ((tstop-tstart)/ntime) + 0.1);

   /* damping factor after index time steps */
   for( i = 1; i <= index; i++ )
      damping_nt *= damping;

   MPI_Comm_rank(comm, &myid);

   sprintf(filename, "%s.%07d.%05d", "drive-02.out", index, myid);
   file = fopen(filename, "w");

   values = (double *) malloc( nlx*nly*sizeof(double) );
   HYPRE_SStructVectorGetBoxValues( u->x, part, app->ilower_x, 
                                    app->iupper_x, var, values );
   m = 0;
   for( j = 0; j < nly; j++ )
      for( i = 0; i < nlx; i++ )
         fprintf(file, "%06d %.14e --- discrete solution %.14e\n", 
                 pj*px*nlx*nly+pi*nlx+j*px*nlx+i,
                 values[m++], 
                 damping_nt*sin((app->ilower_x[0]+i)*(app->dx))
                           *sin((app->ilower_x[1]+j)*(app->dy)));
   fflush(file);
   fclose(file);
   free( values );

   return 0;
}


/* --------------------------------------------------------------------
 * Return buffer size for vector object buffer. Vector object contains
 * values at every grid point and thus, the buffer size is the number
 * of grid points.
 * -------------------------------------------------------------------- */
int
my_BufSize(warp_App  app,
           int      *size_ptr)
{
   *size_ptr = (app->nlx)*(app->nly)*sizeof(double);
   return 0;
}


/* --------------------------------------------------------------------
 * Pack a vector object in a buffer.
 * -------------------------------------------------------------------- */
int
my_BufPack(warp_App     app,
           warp_Vector  u,
           void        *buffer)
{
   double *dbuffer = buffer;

   /* We have one variable and one part. */
   int        part = 0;
   int        var  = 0;

   HYPRE_SStructVectorGather( u->x );
   HYPRE_SStructVectorGetBoxValues( u->x, part, app->ilower_x,
                                    app->iupper_x, var, dbuffer );

   return 0;
}


/* --------------------------------------------------------------------
 * Unpack a vector object from a buffer.
 * -------------------------------------------------------------------- */
int
my_BufUnpack(warp_App     app,
             void        *buffer,
             warp_Vector *u_ptr)
{
   double    *dbuffer = buffer;
   my_Vector *u;

   /* We have one variable and one part. */
   int        part = 0;
   int        var  = 0;

   u = (my_Vector *) malloc( sizeof(my_Vector) );

   /* Create an empty vector object. */
   HYPRE_SStructVectorCreate( app->comm_x, app->grid_x, &(u->x) );
   
   /* Set the object type (by default HYPRE_SSTRUCT). */
   HYPRE_SStructVectorSetObjectType( u->x, app->object_type );

   /* Indicate that the vector coefficients are ready to be set. */
   HYPRE_SStructVectorInitialize( u->x );

   /* Set the values. */
   HYPRE_SStructVectorSetBoxValues( u->x, part, app->ilower_x,
                                    app->iupper_x, var, dbuffer );

   HYPRE_SStructVectorAssemble( u->x );

   *u_ptr = u;

   return 0;
}


/* --------------------------------------------------------------------
 * Main driver
 * -------------------------------------------------------------------- */
int main (int argc, char *argv[])
{
   int i;
   int arg_index;
   int print_usage = 0;

   warp_Core  core;
   my_App    *app;
   int        max_levels;
   int        nrelax;
   double     tol;
   int        cfactor;
   int        max_iter;

   MPI_Comm    comm, comm_x, comm_t;
   int         myid, num_procs;
   int         xcolor, tcolor;

   /* We consider a 2D problem. */
   int ndim = 2;

   /* diffusion coefficient */
   double K; 

   int nlx, nly;
   double tstart, tstop;
   int nt;
   double dx, dy, dt;

   /* We have one part and one variable. */
   int nparts = 1;
   int nvars = 1;
   int var;

   /* We want to use a struct solver. */
   int object_type = HYPRE_STRUCT;
   
   int sym; 

   int px, py, pt;
   int pi, pj;

   int n_pre, n_post;
   int rap, relax, skip;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);

   /* Default parameters. */
   comm        = MPI_COMM_WORLD;
   max_levels  = 1;
   nrelax      = 1;
   tol         = 1.0e-09;
   cfactor     = 2;
   max_iter    = 100;
   K           = 1.0;
   nlx         = 16;
   nly         = 16;
   tstart      = 0.0;
   nt          = 32;
   sym         = 0;
   px          = 1;
   py          = 1;
   pt          = 1;  
   n_pre       = 1;
   n_post      = 1;
   rap         = 0;
   relax       = 2;
   skip        = 0;

   MPI_Comm_rank( comm, &myid );
   MPI_Comm_size( comm, &num_procs );

   /* Parse command line */
   arg_index = 0;
   while( arg_index < argc ){
      if( strcmp(argv[arg_index], "-pgrid") == 0 ){
         arg_index++;
         px = atoi(argv[arg_index++]);
         py = atoi(argv[arg_index++]);
         pt = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-nl") == 0 ){
         arg_index++;
         nlx = atoi(argv[arg_index++]);
	 nly = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-nt") == 0 ){
          arg_index++;
          nt = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-ml") == 0 ){
          arg_index++;
          max_levels = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-nu") == 0 ){
          arg_index++;
          nrelax = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-tol") == 0 ){
          arg_index++;
          tol = atof(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-cf") == 0 ){
          arg_index++;
          cfactor = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-mi") == 0 ){
          arg_index++;
          max_iter = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-v") == 0 ){
         arg_index++;
         n_pre = atoi(argv[arg_index++]);
         n_post = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-rap") == 0 ){
         arg_index++;
         rap = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-relax") == 0 ){
         arg_index++;
         relax = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-skip") == 0 ){
         arg_index++;
         skip = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-sym") == 0 ){
         arg_index++;
         sym = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-help") == 0 )
      {
         print_usage = 1;
         break;
      }
      else
      {
         arg_index++;
      }
   }

   if((print_usage) && (myid == 0)){
      printf("\n");
      printf("Usage: %s [<options>]\n", argv[0]);
      printf("\n");
      printf("  -pgrid  <px py pt>  : processors in each dimension (default: 1 1 1)\n");
      printf("  -nl  <nlx nly>      : problem size per processor (default: 16 16)\n");
      printf("  -nt  <n>            : number of time steps (default: 32)\n"); 
      printf("  -ml  <max_levels>   : set max number of time levels (default: 1)\n");
      printf("  -nu  <nrelax>       : set num F-C relaxations (default: 1)\n");
      printf("  -tol <tol>          : set stopping tolerance (default: 1e-09)\n");
      printf("  -cf  <cfactor>      : set coarsening factor (default: 2)\n");
      printf("  -mi  <max_iter>     : set max iterations (default: 100)\n");
      printf("  -v <n_pre> <n_post> : number of pre and post relaxations in PFMG\n");
      printf("  -rap <r>            : coarse grid operator type in PFMG\n");
      printf("                        0 - Galerkin (default)\n");
      printf("                        1 - non-Galerkin ParFlow operators\n");
      printf("                        2 - Galerkin, general operators\n");
      printf("  -relax <r>          : relaxation type in PFMG\n");
      printf("                        0 - Jacobi\n");
      printf("                        1 - Weighted Jacobi\n");
      printf("                        2 - R/B Gauss-Seidel (default)\n");
      printf("                        3 - R/B Gauss-Seidel (nonsymmetric)\n");
      printf("  -skip <s>           : skip levels in PFMG (0 or 1)\n");      
      printf("  -sym <s>            : symmetric storage (1) or not (0)\n");
      printf("\n");
   }

   if( print_usage ){
      MPI_Finalize();
      return (0);
   }

   /* Check the processor grid (px x py x pt = num_procs?). */
   if( (px*py*pt) > num_procs)
   {
       if( myid == 0 )
           printf("Error: px x py x pt is greater than the number of processors!\n");
       MPI_Finalize();
       return (0);
   }
   if( (px*py*pt) < num_procs)
   {
       if( myid == 0 )
            printf("Error: px x py x pt is less than the number of processors!\n");
       MPI_Finalize();
       return (0);
   }

   /* Create communicators for the time and space dimensions */
   /* The communicators are based on colors and keys (= myid) */
   xcolor = myid / ( px*py );
   tcolor = myid % ( px*py );

   MPI_Comm_split( comm, xcolor, myid, &comm_x );
   MPI_Comm_split( comm, tcolor, myid, &comm_t );

   MPI_Comm_size( comm_t, &num_procs );
/*   printf( "number of processors in time:  %d\n", num_procs );*/
   MPI_Comm_size( comm_x, &num_procs );
/*   printf( "number of processors in space: %d\n", num_procs );*/
   
   /* Determine position in the processor grid. */
   MPI_Comm_rank( comm_x, &myid );
   pi = myid % px;
   pj = ( (myid - pi)/px ) % py;

   /* Compute grid spacing. */
   dx = PI / (px*nlx - 1);
   dy = PI / (py*nly - 1);

   /* Set time-step size. */
   dt = 1.0*(dx*dx);
   /* Determine tstop. */
   tstop =  tstart + nt*dt;

/*   printf("Time interval [%.4lf, %.4lf]\n", tstart, tstop);*/

   /* -----------------------------------------------------------------
    * Set up App structure.
    * ----------------------------------------------------------------- */
   app = (my_App *) malloc(sizeof(my_App));
   (app->comm_t)      = comm_t;
   (app->comm_x)      = comm_x;
   (app->dim_x)       = ndim;
   (app->K)           = K;
   (app->nlx)         = nlx;
   (app->nly)         = nly;
   (app->tstart)      = tstart;
   (app->tstop)       = tstop;
   (app->nt)          = nt;
   (app->dx)          = dx;
   (app->dy)          = dy;
   (app->dt)          = dt;
   (app->max_levels)  = max_levels;
   (app->nparts)      = nparts;
   (app->nvars)       = nvars;
   (app->sym)         = sym;
   (app->px)          = px;
   (app->py)          = py;
   (app->pi)          = pi;
   (app->pj)          = pj;
   (app->object_type) = object_type;
   (app->n_pre)       = n_pre;
   (app->n_post)      = n_post;
   (app->rap)         = rap;
   (app->relax)       = relax;
   (app->skip)        = skip;

   /* Set the variable types. */
   (app->vartypes)   = (HYPRE_SStructVariable*) malloc( nvars* 
                                         sizeof(HYPRE_SStructVariable) );
   for( var = 0; var < nvars; var++ )
      app->vartypes[var] = HYPRE_SSTRUCT_VARIABLE_CELL;
 
   /* Define the nodes owned by the current processor (each processor's
    * piece of the global grid) */
   app->ilower_x[0] = pi*nlx;
   app->ilower_x[1] = pj*nly;
   app->iupper_x[0] = app->ilower_x[0] + nlx-1;
   app->iupper_x[1] = app->ilower_x[1] + nly-1;

   /* Set up a 2D grid. */
   setUp2Dgrid( app->comm_x, &(app->grid_x), app->dim_x, app->nparts, 
                app->ilower_x, app->iupper_x, app->vartypes );

   /* Define the discretization stencil. */
   set5ptStencil( &(app->stencil), app->sym, app->dim_x );

   /* Set up the graph - this determines the non-zero structure of the
    * matrix. */ 
   setUpGraph( app->comm_x, &(app->graph), app->grid_x, app->object_type, 
               app->stencil );

   /* Allocate memory for array of discretization matrices. */
   app->A = (HYPRE_SStructMatrix*) malloc( (app->max_levels)*
                                           sizeof(HYPRE_SStructMatrix));
   /* Create empty matrix lookup table. */
   app->dt_A = (double*) malloc( (app->max_levels)*sizeof(double) );
   for( i = 0; i < app->max_levels; i++ )
      app->dt_A[i] = -1.0;
   app->nA = 0;

   /* Allocate memory for array of solvers. */
   app->solver = (HYPRE_StructSolver*) malloc( (app->max_levels)*
                                               sizeof(HYPRE_StructSolver));

   /* Allocate memory for array of iteration counts. */
   app->max_num_iterations = (int*) calloc( (app->max_levels),  sizeof(int) );

   warp_Init(comm, comm_t, tstart, tstop, nt, app,
             my_Phi, my_Init, my_Clone, my_Free, my_Sum, my_Dot, 
             my_Write, my_BufSize, my_BufPack, my_BufUnpack,
             &core);

   warp_SetMaxLevels( core, max_levels );

   warp_SetNRelax(core, nrelax);

   warp_SetTol(core, tol);

   warp_SetCFactor(core, -1, cfactor);
   /*warp_SetCFactor(core,  0, 10);*/
   
   warp_SetMaxIter(core, max_iter);

   warp_Drive(core);

   warp_PrintStats(core);

   /* Print some additional statistics */
   MPI_Comm_rank( comm, &myid );
   if( myid == 0 )
   {
      for( i = 0; i < app->nA; i++ ){
         printf( "max iterations on level %d: %d\n",
                 i, app->max_num_iterations[i] );
      }
      printf( "\n");
   }

   /* Free memory */
   HYPRE_SStructGridDestroy( app->grid_x );
   HYPRE_SStructStencilDestroy( app->stencil );
   HYPRE_SStructGraphDestroy( app->graph );
   free( app->vartypes );
   free( app->dt_A );
   for( i = 0; i < app->nA; i++ ){
      HYPRE_SStructMatrixDestroy( app->A[i] );
      HYPRE_StructPFMGDestroy( app->solver[i] );
   }
   free( app->A );
   free( app->solver );
   free( app->max_num_iterations );
   free( app );
   warp_Destroy(core);
   MPI_Comm_free( &comm_x );
   MPI_Comm_free( &comm_t );

   /* Finalize MPI */
   MPI_Finalize();

   return 0;
}

