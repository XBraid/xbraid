/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. Email xbraid-support@llnl.gov for support.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (as published by the Free Software
 * Foundation) version 2.1 dated February 1999.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc., 59
 * Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 ***********************************************************************EHEADER*/


/*
   Example 05
   2D diffusion problem

   Interface:    SStructured interface (SStruct)

   Compile with: make drive-05

   Help with:    drive-05 -help

   Sample run:   mpirun -np 8 drive-05 -pgrid 1 1 8 -ml 15 -nt 128 -nx 33 33 -mi 100 -expl -scoarsen  

   Notes:        The " -expl and -scoarsen " options should be used together
                 to coarsen spatially and do explicit time stepping.  

   Description:  This code solves the diffusion problem 

                                  u_t - K * div u = 0

                 in the unit square subject to zero Dirichlet boundary
                 conditions,u_b = 0, and initial condition U0 at time t = 0.
                 
                 The domain is split into an p_x x p_y processor grid.  Each
                 processor has a n_x x n_y grid, with nodes connected by a
                 5-point stencil. More precisely, we use central FD in space
                 and forward or backward Euler in time, definining the stencil

                                    -K*dt/(dy^2)
                   -K*dt/(dx^2)  1+2*K*(dt/(dx^2)+dt/(dy^2)) -K*dt/(dx^2)
                                    -K*dt/(dy^2)   
                     
                  We use cell-centered variables, and, therefore, the nodes are
                  not shared.

                  To incorporate the boundary conditions, we do the following:
                  Let x_i and x_b be the interior and boundary parts of the
                  solution vector x, respectively, and let u_b the boundary
                  condition. If we split the matrix A as

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
#include "_hypre_sstruct_mv.h"

#include "braid.h"
#include "braid_test.h"

#include "vis.c"

#define DEBUG 0

#ifdef M_PI
   #define PI M_PI
#else
   #define PI 3.14159265358979
#endif


/* 
 * This will store the spatial discretization information for interpolating to and 
 * from a pair of levels denoted by the tuple of time step sizes, (cdt, fdt).
 *
 * The app structure will contain a table of spatial_discretization structures and 
 * each braid_Vector will contain an integer index, spatial_disc_idx, which will 
 * denote with entry in the table describes that vector's spatial grid.
 *
 * This table will essentially be a database with the tuples (cdt, fdt) as the keys.
 *
 * fdt                  fine grid time step size         
 * cdt                  coarse grid time step size       
 * dx                   coarse grid dx
 * dy                   coarse grid dy
 * nx                   coarse grid nx (num global x-points)
 * ny                   coarse grid ny (num global y-points)
 * nlx                  coarse grid nlx (num local x-points)
 * nly                  coarse grid nly (num local y-points)
 * ilower               coarse grid ilower
 * iupper               coarse grid iupper
 * grid_x_vector        SStructGrid for coarse grid vectors (has expanded ghost layer for interp)
 * grid_x_matrix        SStructGrid for coarse grid matrices
 * graph_matrix         SStructGraph for coarse grid matrices
 * graph_vector         SStructGraph for coarse grid vectors (has expanded ghost layer for interp)
 * fspatial_disc_idx    integer index into app->spatial_disc_table for the 
 *                      fine spatial mesh used to generate this mesh
 * ncoarsen             the number of coarsenings used to generate this grid from the grid specified
 *                      by fspatial_disc_idx
 * */
typedef struct _spatial_discretization
{
   double fdt, cdt, dx, dy; 
   int ncoarsen, nx, ny, nlx, nly, fspatial_disc_idx;
   int ilower[2], iupper[2];
   HYPRE_SStructGrid       grid_x_vector;
   HYPRE_SStructGrid       grid_x_matrix;
   HYPRE_SStructGraph      graph_vector;
   HYPRE_SStructGraph      graph_matrix;
} spatial_discretization;


/* --------------------------------------------------------------------
 * Structs for my integration routines
 * -------------------------------------------------------------------- */
/* struct my_App contains general information about the problem, its
 * discretizaion, spatial distribution, and solver used at each time point
 *
 *   comm_t              communicator for parallelizing in time
 *   comm_x              communicator for parallelizing in space 
 *   dim_x               spatial dimension
 *   K                   diffusion coefficient
 *   nlx                 local problem size in x-dimension
 *   nly                 local problem size in y-dimension
 *   nx                  global problem size in x-dimension
 *   ny                  global problem size in y-dimension
 *   tstart              initial time
 *   tstop               integration stop time
 *   nt                  number of time steps
 *   dx                  spatial step size in x-direction
 *   dy                  spatial step size in y-direction
 *   dt                  time step on finest grid
 *   max_levels          maximum number of time levels
 *   vartypes            (nvars)-dimensional array of variable types on 
 *                       structured parts of the spatial grid
 *   grid_x              spatial grid
 *   stencil             discretization stencil object
 *   graph               graph object that determine the non-zero structure
 *                       of the discretization matrices
 *   A                   array of discretization matrices (one per time level)
 *   dt_A                array of time steps for which discretization matrix
 *                       has been created
 *   dx_A                array of x-direction mesh sizes for which a discretization 
 *                       matrix has been created
 *   dy_A                array of y-direction mesh sizes for which a discretization 
 *                       matrix has been created
 *   nA                  number of discretization matrices that have been
 *                       created
 *   px                  number of processors in x-dimension
 *   py                  number of processors in y-dimension
 *   pi                  x-coordinate of position in processor grid
 *   pj                  y-coordinate of position in processor grid
 *   ilower              (dim_x)-dimensional array with integer indices of 
 *                       local space interval lower bounds
 *   iupper              (dim_x)-dimensional array with integer indices of
 *                       local space interval upper bounds
 *   object_type         object type of vector to access different hypre solvers 
 *   solver              array of solvers used at each time step on different
 *                       time levels
 *   n_pre               number of pre-relaxation sweeps in spatial solve
 *   n_post              number of post-relaxation sweeps in spatial solve
 *   rap                 coarse-grid operator type for spatial solver
 *                          rap = 0: Galerkin (default)
 *                          rap = 1: non-Galerkin ParFlow operators
 *                          rap = 2: Galerkin, general operators 
 *   relax               relaxation type for spatial solver
 *                          relax = 0: Jacobi
 *                          relax = 1: Weighted Jacobi 
 *                          relax = 2: R/B Gauss-Seidel (default)
 *                          relax = 3: R/B Gauss-Seidel (non-symmetric)
 *   skip                skip levels in spatial PFMG solver (0 or 1)
 *   max_iter_x          expensive and cheap maximum number of spatial MG iterations
 *   tol_x               loose and tight stopping tolerance for spatial MG
 *   explicit            use explicit discretization (1) or not (0)
 *   scoarsen            use spatial refinement and coarsening
 *   num_scoarsen_CFL    the number of entries in scoarsen_CFL
 *   scoarsenCFL         array of CFL numbers such that if the actual CFL at level k is 
 *                       greater that scoarsenCFL[k] then do spatial coarsening until the
 *                       CFL < scoarsenCFL[k]
 *                       Note: scoarsen must also > 0 to turn on spatial coarsening
 *   spatial_disc_table  Lookup table recording dx, dy when coarsening spatially.
 *                       tuples (cdt, fdt) are the keys used to access the spatial
 *                       discretizations stored in the table.  The accessor function
 *                       is get_coarse_spatial_disc( ).  See it for more documentation.
 *   last_tsize          output related flag that lets my_Phi know when the 
 *                       time step size has changed
 *   scheme              int array of integration scheme used: explicit or
 *                       implicit 
 *   output_files        save the solution/error/error norm to files
 *   output_vis          save the error for GLVis visualization
 */
typedef struct _braid_App_struct {
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
   HYPRE_SStructVariable  *vartypes;
   HYPRE_SStructGrid       grid_x;
   HYPRE_SStructStencil    stencil;
   HYPRE_SStructGraph      graph;
   HYPRE_SStructMatrix    *A;
   double                 *dt_A;
   double                 *dx_A;
   double                 *dy_A;
   int                     nA;
   int                     px, py;
   int                     pi, pj;
   int                     ilower[2], iupper[2];
   int                     object_type;
   HYPRE_StructSolver     *solver;
   int                    *max_num_iterations;
   int                     n_pre, n_post;
   int                     rap, relax, skip;
   int                    *max_iter_x;
   double                 *tol_x;
   double                 *scoarsen_table;
   int                     explicit;
   int                     use_rand;
   int                     scoarsen;
   int                     num_scoarsenCFL;
   double                 *scoarsenCFL;
   spatial_discretization *spatial_disc_table;
   double                  last_tsize;
   int                     output_files;
   int                     output_vis;
} my_App;

/* struct my_Vector contains local information specific to one time point
 *   x            spatial vector 
 */
typedef struct _braid_Vector_struct
{
   int                   spatial_disc_idx;
   HYPRE_SStructVector   x;
} my_Vector;

int max( int a, int b ){
  return (a >= b ? a : b );
}

int min( int a, int b ){
  return (a <= b ? a : b );
}


int GetDistribution_x( int    npoints,
                       int    nprocs,
                       int    proc,
                       int   *ilower_ptr,
                       int   *iupper_ptr ){
   int  ilower, iupper;
   int  quo, rem, p;

   quo = npoints/nprocs;
   rem = npoints%nprocs;

   p = proc;
   ilower = p*quo + (p < rem ? p : rem);
   p = proc+1;
   iupper = p*quo + (p < rem ? p : rem) - 1;

   *ilower_ptr = ilower;
   *iupper_ptr = iupper;

   return 0;
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
 * Apply boundary conditions for explicit scheme:
 * Put the boundary conditions in the vector. 
 * -------------------------------------------------------------------- */
void 
addBoundary( HYPRE_SStructVector  b, 
             double               K, 
             double               dx, 
             double               dy, 
             double               dt,
             int                 *ilower,
             int                 *iupper,
             int                  nlx, 
             int                  nly,
             int                  nx, 
             int                  ny,
             int                  px, 
             int                  py, 
             int                  pi, 
             int                  pj )
{
   int i, j;
   
   int bc_ilower[2];
   int bc_iupper[2];
   double *bvalues;

   /* We have one part and one variable. */
   int part = 0;
   int var = 0;
       
   /* Allocate vector for values on boundary planes */
   bvalues = (double *) malloc( (max(nlx, nly)+1)*sizeof(double) );     
  
       
   /* a) boundaries y = 0 or y = PI */
   /* Processors at y = 0 */
   if( ilower[1] == 0 ){
      /* All of proc's x-extents */
      bc_ilower[0] = ilower[0];
      bc_iupper[0] = iupper[0];
      
      /* Only first row of y-extents */
      bc_ilower[1] = ilower[1];
      bc_iupper[1] = ilower[1];
           
      /* Only do work if your box is nonzero in size */
      if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
      {
         /* Put the boundary conditions in b */
         for( i = 0; i < nlx; i++ )
             bvalues[i] = B0( (bc_ilower[0]+i)*dx,
                               bc_ilower[1]*dy );
              
         HYPRE_SStructVectorSetBoxValues(b, part, bc_ilower,
                                         bc_iupper, var, bvalues);
      }
   }
       
   /* Processors at y = PI */
   if( iupper[1] == ny-1 ){
      /* All of proc's x-extents */
      bc_ilower[0] = ilower[0];
      bc_iupper[0] = iupper[0];
      
      /* Only the last row of the y-extents */
      bc_ilower[1] = iupper[1]; 
      bc_iupper[1] = iupper[1];


      /* Only do work if your box is nonzero in size */
      if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
      {
         /* Put the boundary conditions in b */
         for( i = 0; i < nlx; i++ )
             bvalues[i] = B0( (bc_ilower[0]+i)*dx,
                               bc_ilower[1]*dy );
              
         HYPRE_SStructVectorSetBoxValues(b, part, bc_ilower,
                                         bc_iupper, var, bvalues);
      }
   }
       
   /* b) boundaries x = 0 or x = PI */
   /* Processors at x = 0 */
   if( ilower[0] == 0 ){
      /* Only the first column of x-extents */
      bc_ilower[0] = ilower[0];
      bc_iupper[0] = ilower[0];
      
      /* All of the proc's y-extents */
      bc_ilower[1] = ilower[1];
      bc_iupper[1] = iupper[1];

      /* Only do work if your box is nonzero in size */
      if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
      {
         /* Put the boundary conditions in b */
         for( j = 0; j < nly; j++ )
             bvalues[j] = B0(  bc_ilower[0]*dx,
                              (bc_ilower[1]+j)*dy );
              
         HYPRE_SStructVectorSetBoxValues(b, part, bc_ilower,
                                         bc_iupper, var, bvalues);
      }
   }
       
   /* Processors at x = PI */
   if( iupper[0] == nx-1 ){
      
      /* Only the last column of x-extents */
      bc_ilower[0] = iupper[0];  
      bc_iupper[0] = iupper[0];
      
      /* All of the proc's y-extents */
      bc_ilower[1] = ilower[1];
      bc_iupper[1] = iupper[1]; 

      /* Only do work if your box is nonzero in size */
      if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
      {
         /* Put the boundary conditions in b */
         for( j = 0; j < nly; j++ )
             bvalues[j] = B0(  bc_ilower[0]*dx,
                              (bc_ilower[1]+j)*dy );
              
         HYPRE_SStructVectorSetBoxValues(b, part, bc_ilower,
                                         bc_iupper, var, bvalues);
      }
   }

   free(bvalues);

   /* Finalize the vector assembly. */
   HYPRE_SStructVectorAssemble(b);
}


/* --------------------------------------------------------------------
 * Apply boundary conditions for implicit scheme:
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
                  int                 *ilower,
                  int                 *iupper,
                  int                  nlx, 
                  int                  nly,
                  int                  nx, 
                  int                  ny,
                  int                  px, 
                  int                  py, 
                  int                  pi, 
                  int                  pj) 
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
   if( ilower[1] == 0 ){
      /* All of proc's x-extents */
      bc_ilower[0] = ilower[0];
      bc_iupper[0] = iupper[0];
      
      /* Only first row of y-extents */
      bc_ilower[1] = ilower[1];
      bc_iupper[1] = ilower[1];
           
      /* Only do work if your box is nonzero in size */
      if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
      {
         /* Put the boundary conditions in b */
         for( i = 0; i < nlx; i++ )
             bvalues[i] = B0( (bc_ilower[0]+i)*dx,
                               bc_ilower[1]*dy );
              
         HYPRE_SStructVectorSetBoxValues(b, part, bc_ilower,
                                         bc_iupper, var, bvalues);
      }
   }
       
   /* Processors at y = PI */
   if( iupper[1] == ny-1 ){
      /* All of proc's x-extents */
      bc_ilower[0] = ilower[0];
      bc_iupper[0] = iupper[0];
      
      /* Only the last row of the y-extents */
      bc_ilower[1] = iupper[1]; 
      bc_iupper[1] = iupper[1];


      /* Only do work if your box is nonzero in size */
      if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
      {
         /* Put the boundary conditions in b */
         for( i = 0; i < nlx; i++ )
             bvalues[i] = B0( (bc_ilower[0]+i)*dx,
                               bc_ilower[1]*dy );
              
         HYPRE_SStructVectorSetBoxValues(b, part, bc_ilower,
                                         bc_iupper, var, bvalues);
      }
   }
       
   /* b) boundaries x = 0 or x = PI */
   /* Processors at x = 0 */
   if( ilower[0] == 0 ){
      /* Only the first column of x-extents */
      bc_ilower[0] = ilower[0];
      bc_iupper[0] = ilower[0];
      
      /* All of the proc's y-extents */
      bc_ilower[1] = ilower[1];
      bc_iupper[1] = iupper[1];

      /* Only do work if your box is nonzero in size */
      if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
      {
         /* Put the boundary conditions in b */
         for( j = 0; j < nly; j++ )
             bvalues[j] = B0(  bc_ilower[0]*dx,
                              (bc_ilower[1]+j)*dy );
              
         HYPRE_SStructVectorSetBoxValues(b, part, bc_ilower,
                                         bc_iupper, var, bvalues);
      }
   }
       
   /* Processors at x = PI */
   if( iupper[0] == nx-1 ){
      
      /* Only the last column of x-extents */
      bc_ilower[0] = iupper[0];  
      bc_iupper[0] = iupper[0];
      
      /* All of the proc's y-extents */
      bc_ilower[1] = ilower[1];
      bc_iupper[1] = iupper[1]; 

      /* Only do work if your box is nonzero in size */
      if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
      {
         /* Put the boundary conditions in b */
         for( j = 0; j < nly; j++ )
             bvalues[j] = B0(  bc_ilower[0]*dx,
                              (bc_ilower[1]+j)*dy );
              
         HYPRE_SStructVectorSetBoxValues(b, part, bc_ilower,
                                         bc_iupper, var, bvalues);
      }
   }

   
   /* 
    * Now, account for the boundary conditions contributions to the
    * domain interior from A_ib u_b
    */ 
   
   /* a) Neighbors of boundary nodes of boundary y = 0.
    *    These neighbors are in row 1 
    * */
   if( (ilower[1] <=1) && (iupper[1] >= 1) )
   {
      /* All of proc's x-extents */
      bc_ilower[0] = ilower[0];
      bc_iupper[0] = iupper[0];
      
      /* Only the second row of the y-extents */
      bc_ilower[1] = 1;
      bc_iupper[1] = 1;
        
      istart = 0; iend = bc_iupper[0] - bc_ilower[0] + 1;
        
      /* Adjust box to not include boundary nodes */
      if( bc_ilower[0] == 0 ){
         bc_ilower[0] += 1;
         istart += 1;
      }

      if( bc_iupper[0] == nx-1 ){
         bc_iupper[0] -= 1;
         iend -= 1;
      }
     
      /* Only do work if your box is nonzero in size */
      if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
      {
         /* Adjust for removing connections between the boundary
          * and interior nodes in the discretization matrix. */
         for( m = 0, i = istart; i < iend; i++, m++ )
            bvalues[m] = K*(dt/(dy*dy))*
                           B0( (bc_ilower[0]+i-istart)*dx,
                               (bc_ilower[1]-1)*dy );

         HYPRE_SStructVectorAddToBoxValues(b, part, bc_ilower,
                                           bc_iupper, var, bvalues);
      }
   }

   /* b) Neighbors of boundary nodes of boundary y = PI. 
    *    These neighbors are in row ny-2 
    * */
   if( (ilower[1] <= (ny-2) ) && (iupper[1] >= (ny-2)) )
   {
      /* All of proc's x-extents */
      bc_ilower[0] = ilower[0];
      bc_iupper[0] = iupper[0];
      
      /* Only the second to last row of the y-extents */
      bc_ilower[1] = ny-2;
      bc_iupper[1] = ny-2;
        
      istart = 0; iend = bc_iupper[0] - bc_ilower[0] + 1;
        
      /* Adjust box to not include boundary nodes */
      if( bc_ilower[0] == 0 ){
         bc_ilower[0] += 1;
         istart += 1;
         }

      if( bc_iupper[0] == nx-1 ){
         bc_iupper[0] -= 1;
         iend -= 1;
      }
     
      /* Only do work if your box is nonzero in size */
      if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
      {
         /* Adjust for removing connections between the boundary
          * and interior nodes in the discretization matrix. */ 
         for( m = 0, i = istart; i < iend; i++, m++ )
            bvalues[m] = K*(dt/(dy*dy))*
                           B0( (bc_ilower[0]+i-istart)*dx,
                               (bc_ilower[1]+1)*dy );
         
         HYPRE_SStructVectorAddToBoxValues(b, part, bc_ilower,
                                           bc_iupper, var, bvalues);
      }
   }

   /* c) Neighbors of boundary nodes of boundary x = 0. 
    * These neighbors or in column 1 
    * */
   if( (ilower[0] <= 1) && (iupper[0] >= 1) )
   {
      /* Only the first column */
      bc_ilower[0] = 1;  
      bc_iupper[0] = 1;
      
      /* All of the proc's y-extents */
      bc_ilower[1] = ilower[1];
      bc_iupper[1] = iupper[1];
        
      jstart = 0; jend = bc_iupper[1] - bc_ilower[1] + 1;
        
      /* Adjust box to not include boundary nodes */
      if( bc_ilower[1] == 0 ){
         bc_ilower[1] += 1;
         jstart += 1;
      }
      
      if( bc_iupper[1] == ny-1 ){
         bc_iupper[1] -= 1;
         jend -= 1;
      }          
      
      /* Only do work if your box is nonzero in size */
      if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
      {
         /* Adjust for removing connections between the boundary
          * and interior nodes in the discretization matrix. */
         for( m = 0, j = jstart; j < jend; j++, m++ )
            bvalues[m] = K*(dt/(dx*dx))*
                           B0( (bc_ilower[0]-1)*dx,
                               (bc_ilower[1]+j-jstart)*dy );
         
         HYPRE_SStructVectorAddToBoxValues(b, part, bc_ilower,
                                           bc_iupper, var, bvalues);
      }
   }

   /* d) Neighbors of boundary nodes of boundary x = PI.
    * These neighbors or in column nx-2 
    * */
   if( (ilower[0] <= (nx-2) ) && (iupper[0] >= (nx-2)) )
   {
      /* Only the second to last column */
      bc_ilower[0] = nx-2;  
      bc_iupper[0] = nx-2;
      
      /* All of the proc's y-extents */
      bc_ilower[1] = ilower[1];
      bc_iupper[1] = iupper[1];

      jstart = 0; jend = bc_iupper[1] - bc_ilower[1] + 1;
        
      /* Adjust box to not include boundary nodes */
      if( bc_ilower[1] == 0 ){
         bc_ilower[1] += 1;
         jstart += 1;
      }

      if( bc_iupper[1] == ny-1 ){
         bc_iupper[1] -= 1;
         jend -= 1;
      }

      /* Only do work if your box is nonzero in size */
      if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
      {
         /* Adjust for removing connections between the boundary
          * and interior nodes in the discretization matrix. */
         for( m = 0, j = jstart; j < jend; j++, m++ )
            bvalues[m] = K*(dt/(dx*dx))*
                           B0( (bc_ilower[0]+1)*dx,
                               (bc_ilower[1]+j-jstart)*dy );
         
         HYPRE_SStructVectorAddToBoxValues(b, part, bc_ilower,
                                           bc_iupper, var, bvalues);
      }
   }
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
             int                   *ilower, 
             int                   *iupper,
             HYPRE_SStructVariable *vartypes,
             int                    nghost)
{
   /* We have one part and one variable. */
   int nvars = 1;
   int nparts = 1;
   int part = 0;
   int num_ghost[4];
   
   num_ghost[0] = nghost;
   num_ghost[1] = nghost;
   num_ghost[2] = nghost;
   num_ghost[3] = nghost;

   HYPRE_SStructGrid grid;

   /* Create an empty 2D grid object. */
   HYPRE_SStructGridCreate( comm, ndim, nparts, &grid );
   
   /* Set the variable type for each part 
    * This call MUST go before setting the number of ghost, i.e.,
    * nvars must be > 0 for that SetNumGhost to have any effect */
   HYPRE_SStructGridSetVariables( grid, part, nvars, vartypes );

   /* Add a new box to the grid. */
   if ((ilower[0] <= iupper[0]) && (ilower[1] <= iupper[1]))
   {
      HYPRE_SStructGridSetExtents( grid, part, ilower, iupper );
      HYPRE_SStructGridSetNumGhost(grid, num_ghost);
   }

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
               int                   ndim )
{
   /* We have one variable. */
   int var = 0;

   int i;

   HYPRE_SStructStencil stencil;

   /* Create an empty 2D, 5-pt stencil object. */
   HYPRE_SStructStencilCreate( ndim, 5, &stencil );

   /* Define the geometry of the stencil. */
   int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};

   /* Assign stencil entries. */
   for( i = 0; i < 5; i++ )
      HYPRE_SStructStencilSetEntry( stencil, i, offsets[i], var );

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
 * Set up the implicit discretization matrix. 
 * First, we set the stencil values at every node neglecting the 
 * boundary. Then, we correct the matrix stencil at boundary nodes.
 * We have to eliminate the coefficients reaching outside of the domain
 * boundary. Furthermore, to incorporate boundary conditions, we remove
 * the connections between the interior nodes and boundary nodes.
 * -------------------------------------------------------------------- */
void
setUpImplicitMatrix( MPI_Comm             comm,
                     HYPRE_SStructMatrix *A_ptr,
                     HYPRE_SStructGraph   graph, 
                     int                  object_type,
                     double               K, 
                     double               dx, 
                     double               dy, 
                     double               dt,
                     int                 *ilower, 
                     int                 *iupper,
                     int                  nlx, 
                     int                  nly,
                     int                  nx, 
                     int                  ny,
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
   HYPRE_SStructMatrixSetSymmetric( A, part, var, var, 0);

   /* Set the object type. */
   HYPRE_SStructMatrixSetObjectType( A, object_type );

   /* Indicate that the matrix coefficients are ready to be set. */
   HYPRE_SStructMatrixInitialize( A );

   /* Only do work if your box is nonzero in size */
   if( (ilower[0] <= iupper[0]) && (ilower[1] <= iupper[1]) )
   {
      /* 1. neglect boundaries */
      /* Set the stencil values. */
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

#if DEBUG
printf( "My box: [%d %d] x [%d %d]\n", ilower[0], iupper[0], ilower[1], iupper[1] );
#endif

      /* 2. correct stencils at boundary nodes */     
      /* Allocate vectors for values on boundary planes */
      values  = (double *) malloc( nentries*( max(nlx,nly)+1 )*sizeof(double) );
      bvalues = (double *) malloc( (max(nlx,nly)+1)*sizeof(double) );     
      for( i = 0; i < nentries*max(nlx,nly); i+= nentries ){
         values[i] = 1.0;
         for( idx = 1; idx < nentries; idx++ )
            values[i+idx] = 0.0;
      }
          
      /* a) boundaries y = 0 or y = PI */
      /* The stencil at the boundary nodes is 1-0-0-0-0. Because
       * we have I x_b = u_b. */
          
      /* Processors at y = 0 */
      if( ilower[1] == 0 ){
         /* All of proc's x-extents */
         bc_ilower[0] = ilower[0];
         bc_iupper[0] = iupper[0];
         
         /* Only the first row of the y-extents */
         bc_ilower[1] = ilower[1];
         bc_iupper[1] = ilower[1];
              
         /* Only do work if your box is nonzero in size */
         if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
         {
            /* Modify the matrix. */
            HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                            var, nentries,
                                            stencil_indices, values);
         }
      }
          
      /* Processors at y = PI */
      if( iupper[1] == ny-1 ){
         /* All of proc's x-extents */
         bc_ilower[0] = ilower[0];
         bc_iupper[0] = iupper[0];
         
         /* Only the last row of the y-extents */
         bc_ilower[1] = iupper[1]; 
         bc_iupper[1] = iupper[1];
              
         /* Only do work if your box is nonzero in size */
         if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
         {
            /* Modify the matrix */
            HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                            var, nentries,
                                            stencil_indices, values);
         }
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
      if( ilower[0] == 0 ){
         /* Only the first column of x-extents */
         bc_ilower[0] = ilower[0];
         bc_iupper[0] = ilower[0];
         
         /* All of the proc's y-extents */
         bc_ilower[1] = ilower[1];
         bc_iupper[1] = iupper[1];
              
         /* Only do work if your box is nonzero in size */
         if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
         {
            /* Modify the matrix */
            HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                            var, nentries,
                                            stencil_indices, values);
         }
      }
          
      /* Processors at x = PI */
      if( iupper[0] == nx-1 ){
         /* Only the last column of x-extents */
         bc_ilower[0] = iupper[0];  
         bc_iupper[0] = iupper[0];
         
         /* All of the proc's y-extents */
         bc_ilower[1] = ilower[1];
         bc_iupper[1] = iupper[1];
              
         /* Only do work if your box is nonzero in size */
         if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
         {
            /* Modify the matrix */
            HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                            var, nentries,
                                            stencil_indices, values);
         }
      }

      free(values);    
          
      /* Recall that the system we are solving is:
       *
       *   [A_ii 0; [x_i;    [b_i - A_ib*u_b;
       *      0  I]  x_b ] =        u_b       ].
       * 
       * This requires removing the connections between the interior
       * and boundary nodes that we have set up when we set the
       * 5pt stencil at each node. */
          
      /* a) Neighbors of boundary nodes of boundary y = 0.
       *    These neighbors are in row 1 
       * */
      if( (ilower[1] <=1) && (iupper[1] >= 1) )
      {
         /* All of proc's x-extents */
         bc_ilower[0] = ilower[0];
         bc_iupper[0] = iupper[0];
         
         /* Only the second row of the y-extents */
         bc_ilower[1] = 1;
         bc_iupper[1] = 1;
              
         stencil_indices[0] = 3;
              
         /* Modify the matrix */
         for( m = 0; m < (iupper[0] - ilower[0] + 1); m++ )
            bvalues[m] = 0.0;
              
         /* Only do work if your box is nonzero in size */
         if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
         {
            HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                            var, 1,
                                            stencil_indices, bvalues);
         }
      }

      /* b) Neighbors of boundary nodes of boundary y = PI. 
       *    These neighbors are in row ny-2 
       * */
      if( (ilower[1] <= (ny-2) ) && (iupper[1] >= (ny-2)) )
      {
         /* All of proc's x-extents */
         bc_ilower[0] = ilower[0];
         bc_iupper[0] = iupper[0];
         
         /* Only the second to last row of the y-extents */
         bc_ilower[1] = ny-2;
         bc_iupper[1] = ny-2;
              
         stencil_indices[0] = 4;
              
         /* Modify the matrix */
         for( m = 0; m < (iupper[0] - ilower[0] + 1); m++ )
            bvalues[m] = 0.0;
              
         /* Only do work if your box is nonzero in size */
         if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
         {
            HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                            var, 1, stencil_indices, bvalues);
         }
      }

      /* c) Neighbors of boundary nodes of boundary x = 0. 
       * These neighbors or in column 1 
       * */
      if( (ilower[0] <= 1) && (iupper[0] >= 1) )
      {
         /* Only the first column */
         bc_ilower[0] = 1;  
         bc_iupper[0] = 1;
         
         /* All of the proc's y-extents */
         bc_ilower[1] = ilower[1];
         bc_iupper[1] = iupper[1];

         stencil_indices[0] = 1;
              
         /* Modify the matrix */
         for( m = 0; m < (iupper[1] - ilower[1] + 1); m++ )
            bvalues[m] = 0.0;
              
         /* Only do work if your box is nonzero in size */
         if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
         {
            HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                            var, 1,
                                            stencil_indices, bvalues);
         }
      }

      /* d) Neighbors of boundary nodes of boundary x = PI.
       * These neighbors or in column nx-2 
       * */
      if( (ilower[0] <= (nx-2) ) && (iupper[0] >= (nx-2)) )
      {
         /* Only the second to last column */
         bc_ilower[0] = nx-2;  
         bc_iupper[0] = nx-2;
         
         /* All of the proc's y-extents */
         bc_ilower[1] = ilower[1];
         bc_iupper[1] = iupper[1];
              
         stencil_indices[0] = 2;
              
         /* Modify the matrix */
         for( m = 0; m < (iupper[1] - ilower[1] + 1); m++ )
            bvalues[m] = 0.0;
              
         /* Only do work if your box is nonzero in size */
         if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
         {
            HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                            var, 1, stencil_indices, bvalues);
         }
      }

      free(bvalues);
   }
   /* Finalize the matrix assembly. */
   HYPRE_SStructMatrixAssemble( A );

   *A_ptr = A;
}



/* --------------------------------------------------------------------
 * Set up the explicit discretization matrix. 
 * First, we set the stencil values at every node neglecting the 
 * boundary. Then, we correct the matrix stencil at boundary nodes.
 * We have to eliminate the coefficients reaching outside of the domain
 * boundary. 
 * -------------------------------------------------------------------- */
void
setUpExplicitMatrix( MPI_Comm             comm,
                     HYPRE_SStructMatrix *A_ptr,
                     HYPRE_SStructGraph   graph, 
                     int                  object_type,
                     double               K, 
                     double               dx, 
                     double               dy, 
                     double               dt,
                     int                 *ilower, 
                     int                 *iupper,
                     int                  nlx, 
                     int                  nly,
                     int                  nx, 
                     int                  ny,
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

   double *values;
   int bc_ilower[2];
   int bc_iupper[2];

   HYPRE_SStructMatrix A;

   /* Create an empty matrix object. */
   HYPRE_SStructMatrixCreate( comm, graph, &A);

   HYPRE_SStructMatrixSetSymmetric( A, part, var, var, 0);

   /* Set the object type. */
   HYPRE_SStructMatrixSetObjectType( A, object_type );

   /* Indicate that the matrix coefficients are ready to be set. */
   HYPRE_SStructMatrixInitialize( A );
   
   /* Only do work if your box is nonzero in size */
   if( (ilower[0] <= iupper[0]) && (ilower[1] <= iupper[1]) )
   {
      /* 1. neglect boundaries */
      /* Set the stencil values. */
      nentries = 5;
      values = (double *) malloc( nentries*(nlx*nly)*sizeof(double) );

      /* The order is left-to-right, bottom-to-top */
      for( m = 0, j = 0; j < nly; j++ )
         for( i = 0; i < nlx; i++, m+=nentries ){
            for( idx = 1; idx <= 2; idx++ )
               values[m+idx] = K*(dt/(dx*dx));
            for( idx = 3; idx <= 4; idx++ )
               values[m+idx] = K*(dt/(dy*dy));
                     
            values[m] = 1.0 - 2*K*( (dt/(dx*dx)) + (dt/(dy*dy)) );
         }

      HYPRE_SStructMatrixSetBoxValues( A, part, ilower, iupper, var, 
                                       nentries, stencil_indices, 
                                       values );

      free(values);

      /* 2. correct stencils at boundary nodes */     
      /* Allocate vectors for values on boundary planes */
      values  = (double *) malloc( nentries*( max(nlx,nly)+1 )*
                                            sizeof(double) );   
          
      /* a) boundaries y = 0 or y = PI */
      /* The stencil at the boundary nodes is 1-0-0-0-0. Because
       * we have I x_b = u_b. */
      for( i = 0; i < nentries*nlx; i+= nentries ){
         values[i] = 1.0;
         for( idx = 1; idx < nentries; idx++ )
            values[i+idx] = 0.0;
      }
          
      /* Processors at y = 0 */
      if( ilower[1] == 0 ){
         /* All of proc's x-extents */
         bc_ilower[0] = ilower[0];
         bc_iupper[0] = iupper[0];
         
         /* Only the first row of the y-extents */
         bc_ilower[1] = ilower[1];
         bc_iupper[1] = ilower[1];
              
         /* Only do work if your box is nonzero in size */
         if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
         {
            /* Modify the matrix. */
            HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                            var, nentries,
                                         stencil_indices, values);
         }
      }
          
      /* Processors at y = PI */
      if( iupper[1] == ny-1 ){
         /* All of proc's x-extents */
         bc_ilower[0] = ilower[0];
         bc_iupper[0] = iupper[0];
         
         /* Only the last row of the y-extents */
         bc_ilower[1] = iupper[1]; 
         bc_iupper[1] = iupper[1];
              
         /* Only do work if your box is nonzero in size */
         if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
         {
            /* Modify the matrix */
            HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                            var, nentries,
                                            stencil_indices, values);
         }
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
      if( ilower[0] == 0 ){
         /* Only the first column of x-extents */
         bc_ilower[0] = ilower[0];
         bc_iupper[0] = ilower[0];
         
         /* All of the proc's y-extents */
         bc_ilower[1] = ilower[1];
         bc_iupper[1] = iupper[1];
              
         /* Only do work if your box is nonzero in size */
         if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
         {
            /* Modify the matrix */
            HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                            var, nentries,
                                            stencil_indices, values);
         }
      }
          
      /* Processors at x = PI */
      if( iupper[0] == nx-1 ){
         /* Only the last column of x-extents */
         bc_ilower[0] = iupper[0];  
         bc_iupper[0] = iupper[0];
         
         /* All of the proc's y-extents */
         bc_ilower[1] = ilower[1];
         bc_iupper[1] = iupper[1];
              
         /* Only do work if your box is nonzero in size */
         if( (bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]) )
         {
            /* Modify the matrix */
            HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                            var, nentries,
                                            stencil_indices, values);
         }
      }

      free(values);   
   }

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
                   int                  max_iter,
                   double               tol,
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
   HYPRE_StructPFMGSetTol( solver, tol );
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

   HYPRE_StructPFMGSetMaxIter( solver, max_iter );

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
 *   u_i = Phi_i(u_{i-1}) + g_i 
 * Note that the first case corresponds to assuming zero Dirichlet BCs
 * and a zero RHS of the PDE.
 * When Phi is called, u is u_{i-1}. At the end of the routine, u is 
 * set to u_i.
 * -------------------------------------------------------------------- */
int
my_Phi(braid_App       app,
       braid_Vector    u,
       braid_PhiStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   double accuracy;
   int i, A_idx;
   double *values;
   double cfl_value;
   int ilower[2], iupper[2];
   double dx      = (app->spatial_disc_table[u->spatial_disc_idx]).dx;
   double dy      = (app->spatial_disc_table[u->spatial_disc_idx]).dy;
   int    nlx     = (app->spatial_disc_table[u->spatial_disc_idx]).nlx;
   int    nly     = (app->spatial_disc_table[u->spatial_disc_idx]).nly;
   int    nx      = (app->spatial_disc_table[u->spatial_disc_idx]).nx;
   int    ny      = (app->spatial_disc_table[u->spatial_disc_idx]).ny;

   ilower[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower[0];
   ilower[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower[1];
   iupper[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper[0];
   iupper[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper[1];

   HYPRE_SStructVector b;
   HYPRE_StructMatrix  sA;
   HYPRE_StructVector  sb;
   HYPRE_StructVector  sx;
   
   braid_PhiStatusGetTstartTstop(status, &tstart, &tstop);
   braid_PhiStatusGetAccuracy(status, &accuracy);
   
   /* We have one part and one variable. */
   int part = 0;
   int var = 0;

   int num_iterations, cfl = 0;

   /* -----------------------------------------------------------------
    * Set up the discretization matrix.
    * If no variable coefficients, check matrix lookup table if matrix 
    * has already been created for time step size tstop-tstart and the 
    * mesh size dx and dy.
    * ----------------------------------------------------------------- */
   A_idx = -1.0;
   for( i = 0; i < app->max_levels; i++ ){
      if( app->dt_A[i] == -1.0 )
      {
         break;
      }
      if( (fabs( app->dt_A[i] - (tstop-tstart) )/(tstop-tstart) < 1e-10) &&
          (fabs( app->dx_A[i] - dx )/dx < 1e-10) &&
          (fabs( app->dy_A[i] - dy )/dy < 1e-10) )
      { 
         A_idx = i;
         break;
      }
   }

   /* Check CFL condition, always switch to implicit time stepping if you violate the CFL */
   cfl_value = (app->K)*( (tstop-tstart)/((dx)*(dx)) + (tstop-tstart)/((dy)*(dy)) );
   if( cfl_value < 0.5 )
   {
      cfl = 1;
   } 
   
   /* Store information on CFL and spatial coarsening for user output */
   if( A_idx == -1.0)
   {
      (app->scoarsen_table)[ (5*i) + 1] = dx;
      (app->scoarsen_table)[ (5*i) + 2] = dy;
      (app->scoarsen_table)[ (5*i) + 3] = (tstop-tstart);
      (app->scoarsen_table)[ (5*i) + 4] = cfl_value;
   }

   if( A_idx == -1.0 ){
      A_idx = i;
      app->nA++;
#if DEBUG
      printf( "Create new matrix %d\n", A_idx );
#endif
      /* No matrix for time step tstop-tstart exists. 
       * Add entry to matrix lookup table. */   
      
      app->dt_A[A_idx] = tstop-tstart;
      app->dy_A[A_idx] = dx;
      app->dx_A[A_idx] = dy;

      /* If we want to use an explicit scheme, check CFL condition 
       * to determine whether we can still use explicit scheme for
       * this time step size or if we have to switch to implicit 
       * scheme. */
      if( app->explicit && cfl ){
         /* Set up the explicit discretization matrix. */
         setUpExplicitMatrix( app->comm_x, &(app->A[A_idx]), 
                              (app->spatial_disc_table[u->spatial_disc_idx]).graph_matrix, 
                              app->object_type, app->K, dx, 
                              dy, app->dt_A[A_idx], ilower, 
                              iupper, nlx, nly, nx, ny,
                              app->px, app->py, app->pi, app->pj );
         (app->scoarsen_table)[ (5*i) ]    = 1;
      }
      else{
         /* Set up the implicit discretization matrix. */
         setUpImplicitMatrix( app->comm_x, &(app->A[A_idx]), 
                              (app->spatial_disc_table[u->spatial_disc_idx]).graph_matrix, 
                              app->object_type, app->K, dx, 
                              dy, app->dt_A[A_idx], ilower, 
                              iupper, nlx, nly, nx, ny,
                              app->px, app->py, app->pi, app->pj );
      
         /*HYPRE_SStructMatrixPrint("sstruct.out.A",  app->A[A_idx], 0);*/
         /* Set up the PFMG solver using u->x as dummy vectors. */
         setUpStructSolver( app->comm_x, &(app->solver[A_idx]), &sA, &sb, 
                            &sx, app->A[A_idx], u->x, u->x, 
                            app->max_iter_x[0], app->tol_x[0], app->rap,
                            app->relax, app->n_pre, app->n_post, 
                            app->skip );
         (app->scoarsen_table)[ (5*i) ]    = 0;
      }
   } 

   if( app->explicit && cfl ){
      /* Set up a vector for the MatVec result and modify the solution
       * from the previous time step to incorporate the boundary 
       * conditions. */    
      values = (double *) malloc( (nlx)*(nly)*sizeof(double) );
      HYPRE_SStructVectorGather( u->x );
      HYPRE_SStructVectorGetBoxValues( u->x, part, ilower,
                                       iupper, var, values );
      HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[u->spatial_disc_idx]).grid_x_matrix, &b );
      HYPRE_SStructVectorSetObjectType( b, app->object_type );
      HYPRE_SStructVectorInitialize( b );
      HYPRE_SStructVectorSetBoxValues( b, part, ilower,
                                       iupper, var, values );
      free( values );

      addBoundary( b, app->K, dx, dy, tstop-tstart,
                   ilower, iupper, nlx, nly, nx, ny,
                   app->px, app->py, app->pi, app->pj );

      /* --------------------------------------------------------------
       * Time integration to next time point: Perform MatVec x = Ab.
       * -------------------------------------------------------------- */
      HYPRE_SStructMatrixGetObject( app->A[A_idx], (void **) &sA );
      HYPRE_SStructVectorGetObject( u->x, (void **) &sx );
      HYPRE_SStructVectorGetObject( b, (void **) &sb );
      HYPRE_StructMatrixMatvec( 1, sA, sb, 0, sx );

      /* free memory */
      HYPRE_SStructVectorDestroy( b );
   }
   else{
      /* --------------------------------------------------------------
       * Set up the right-hand side vector. 
       * The right-hand side is the solution from the previous time 
       * step modified to incorporate the boundary conditions and the 
       * right-hand side of the PDE (here, we assume a zero RHS of the
       * PDE).
       * -------------------------------------------------------------- */
      values = (double *) malloc( (nlx)*(nly)*sizeof(double) );
      HYPRE_SStructVectorGather( u->x );
      HYPRE_SStructVectorGetBoxValues( u->x, part, ilower,
                                       iupper, var, values );
      HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[u->spatial_disc_idx]).grid_x_matrix, &b );
      HYPRE_SStructVectorSetObjectType( b, app->object_type );
      HYPRE_SStructVectorInitialize( b );
      HYPRE_SStructVectorSetBoxValues( b, part, ilower,
                                       iupper, var, values );
      free( values );
      addBoundaryToRHS( b, app->K, dx, dy, tstop-tstart,
                        ilower, iupper, nlx, nly, nx, ny,
                        app->px, app->py, app->pi, app->pj );
      /* add infos from RHS of PDE here */ 

      /* --------------------------------------------------------------
       * Time integration to next time point: Solve the system Ax = b.
       * -------------------------------------------------------------- */
      HYPRE_SStructMatrixGetObject( app->A[A_idx], (void **) &sA );
      HYPRE_SStructVectorGetObject( b, (void **) &sb );
      HYPRE_SStructVectorGetObject( u->x, (void **) &sx );
      HYPRE_StructPFMGSetTol( app->solver[A_idx], accuracy );

      if( A_idx == 0 ){
         /* Fine level: use expensive max iters. */
         HYPRE_StructPFMGSetMaxIter( app->solver[A_idx], 
                                     app->max_iter_x[0] );
      }
      else{
         /* Coarse level: use cheap max iters. */
         HYPRE_StructPFMGSetMaxIter( app->solver[A_idx], 
                                     app->max_iter_x[1] );
      }

      HYPRE_StructPFMGSolve( app->solver[A_idx], sA, sb, sx );
      HYPRE_StructPFMGGetNumIterations( app->solver[A_idx], 
                                        &num_iterations );

#if DEBUG
      {
         int myid;
         MPI_Comm_rank(MPI_COMM_WORLD, &myid);
         if( myid == 0 )
            if( A_idx == 0 )
               printf( "tstart = %f: number of PFMG iterations: %d\n", 
                       tstart, num_iterations );

      }
#endif

      app->max_num_iterations[A_idx] = max((app->max_num_iterations[A_idx]),
                                           num_iterations);

      /* free memory */
      HYPRE_SStructVectorDestroy( b );
   }
   /* no refinement */
   braid_PhiStatusSetRFactor(status, 1);

   return 0;
}


/* --------------------------------------------------------------------
 * Create a vector object for a given time point.
 * This function is only called on the finest level
 * -------------------------------------------------------------------- */
int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;
   double    *values;
   int        i, j, m;
   
   /* We have one part */
   int        part = 0;
   int        var = 0;

   /* Guarantee reproducibility by seeding rand every call */
   srand(0);

   u = (my_Vector *) malloc( sizeof(my_Vector) );

   /* Create an empty vector object. */
   HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[0]).grid_x_matrix, &(u->x) );
   
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
            values[m] = U0( (app->ilower[0]+i)*(app->dx), 
                            (app->ilower[1]+j)*(app->dy) );
   }
   else
   {
      /* Random between 0 and 1 */
      for( m = 0; m < ((app->nlx)*(app->nly)); m++ )
         /* Random between 0 and 1 */
         if(app->use_rand)
         {  values[m] = ((double)rand())/RAND_MAX;}
         else
         {  values[m] = 1.42; }
   }
   HYPRE_SStructVectorSetBoxValues( u->x, part, app->ilower,
                                    app->iupper, var, values ); 

   HYPRE_SStructVectorAssemble( u->x );
   free(values);
   
   /* Store the coarsening rule (here on level 0, it's just 0) */
   u->spatial_disc_idx = 0; 

   *u_ptr = u;

   return 0;
}


/* --------------------------------------------------------------------
 * Create a a copy of a vector object.
 * -------------------------------------------------------------------- */
int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   /* We have one part and one var*/
   int        part = 0;
   int        var = 0;
   my_Vector *v;
   double    *values;
   
   int ilower[2], iupper[2];
   int    nlx     = (app->spatial_disc_table[u->spatial_disc_idx]).nlx;
   int    nly     = (app->spatial_disc_table[u->spatial_disc_idx]).nly;

   ilower[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower[0];
   ilower[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower[1];
   iupper[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper[0];
   iupper[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper[1];


   v = (my_Vector *) malloc(sizeof(my_Vector));

   /* Create an empty vector object. */
   HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[u->spatial_disc_idx]).grid_x_matrix, &(v->x) );
   
   /* Set the object type (by default HYPRE_SSTRUCT) and initialize. */
   HYPRE_SStructVectorSetObjectType( v->x, app->object_type );
   HYPRE_SStructVectorInitialize( v->x );

   /* Set the values. */
   values = (double *) malloc( (nlx)*(nly)*sizeof(double) );
   HYPRE_SStructVectorGather( u->x );
   HYPRE_SStructVectorGetBoxValues( u->x, part, ilower,
                                    iupper, var, values );
   HYPRE_SStructVectorSetBoxValues( v->x, part, ilower,
                                    iupper, var, values );
   free( values );
   HYPRE_SStructVectorAssemble( v->x );

   /* Store the spatial_disc_idx used to generate u */
   v->spatial_disc_idx = u->spatial_disc_idx; 

   *v_ptr = v;

   return 0;
}



/* --------------------------------------------------------------------
 * Destroy vector object.
 * -------------------------------------------------------------------- */
int
my_Free(braid_App    app,
        braid_Vector u)
{
   HYPRE_SStructVectorDestroy( u->x );
   free( u );

   return 0;
}


/* --------------------------------------------------------------------
 * Compute vector sum y = alpha*x + beta*y.
 * -------------------------------------------------------------------- */
int
my_Sum(braid_App    app,
       double       alpha,
       braid_Vector x,
       double       beta,
       braid_Vector y)
{
   int i;
   double *values_x, *values_y;
   /* We have one part and one variable */
   int part = 0;
   int var = 0;
   
   int ilower[2], iupper[2];
   int    nlx     = (app->spatial_disc_table[x->spatial_disc_idx]).nlx;
   int    nly     = (app->spatial_disc_table[x->spatial_disc_idx]).nly;

   ilower[0]    = (app->spatial_disc_table[x->spatial_disc_idx]).ilower[0];
   ilower[1]    = (app->spatial_disc_table[x->spatial_disc_idx]).ilower[1];
   iupper[0]    = (app->spatial_disc_table[x->spatial_disc_idx]).iupper[0];
   iupper[1]    = (app->spatial_disc_table[x->spatial_disc_idx]).iupper[1];

   values_x = (double *) malloc( nlx*nly*sizeof(double) );
   values_y = (double *) malloc( nlx*nly*sizeof(double) );

   HYPRE_SStructVectorGather( x->x );
   HYPRE_SStructVectorGetBoxValues( x->x, part, ilower,
                                    iupper, var, values_x );

   HYPRE_SStructVectorGather( y->x );
   HYPRE_SStructVectorGetBoxValues( y->x, part, ilower,
                                    iupper, var, values_y );

   for( i = 0; i < nlx*nly; i++ )
   {
      values_y[i] = alpha*values_x[i] + beta*values_y[i];
   }

   HYPRE_SStructVectorSetBoxValues( y->x, part, ilower,
                                    iupper, var, values_y );

   free( values_x );
   free( values_y );

   return 0;
}



/* --------------------------------------------------------------------
 * Compute spatial norm 
 * -------------------------------------------------------------------- */
int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   double dot;

   hypre_SStructInnerProd( u->x, u->x, &dot );

   *norm_ptr = sqrt(dot);

   return 0;
}


/* --------------------------------------------------------------------
 * Access the vector.
 * -------------------------------------------------------------------- */
int
my_Access(braid_App           app,
          braid_Vector        u,
          braid_AccessStatus  astatus)
{
   MPI_Comm   comm   = MPI_COMM_WORLD;
   double     tstart = (app->tstart);
   double     tstop  = (app->tstop);
   int        ntime  = (app->nt);
   double     rnorm;
   int        iter, level, done;
   int        index, myid, myid_x;
   static int previous_level = -5;
   char       filename[255];
   FILE      *file;
   double     t;

   double    *values;
   /* We have one part and one variable. */
   int part = 0;
   int var = 0;

   /* error norm */
   double enorm = 0.0;
   double enorm_sum = 0.0;

   /* error vector */
   HYPRE_SStructVector e;

   int i, j, m;
   
   int ilower[2], iupper[2];
   int    nlx     = (app->spatial_disc_table[u->spatial_disc_idx]).nlx;
   int    nly     = (app->spatial_disc_table[u->spatial_disc_idx]).nly;
   double dx      = (app->spatial_disc_table[u->spatial_disc_idx]).dx;
   double dy      = (app->spatial_disc_table[u->spatial_disc_idx]).dy;

   ilower[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower[0];
   ilower[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower[1];
   iupper[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper[0];
   iupper[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper[1];
   
   /* Retrieve current time from Status Object */
   braid_AccessStatusGetT(astatus, &t);

   /* Retrieve Braid State Information from Status Object */
   MPI_Comm_rank(comm, &myid);
   braid_AccessStatusGetTILD(astatus, &t, &iter, &level, &done);
   braid_AccessStatusGetResidual(astatus, &rnorm);
   if( (myid == 0) && (level != previous_level) )
   {
      previous_level = level;
      printf("  my_Access() called, iter= %d, level= %d\n", iter, level);
   }

   /* Write to files:
    *   - save computed solution and true discrete solution
    *     if we want to visualize with GLVis, we also save the solution
    *     at the initial time, middle time, and end time
    *   - save error norm at each time point
    *     if we want to visualize with GLVis, we also save the error
    *     at the initial time, middle time, and end time */
   if( app->output_files && (level == 0) )
   {

      index = ((t-tstart) / ((tstop-tstart)/ntime) + 0.1);

      values = (double *) malloc( nlx*nly*sizeof(double) );
      HYPRE_SStructVectorGetBoxValues( u->x, part, ilower, 
                                       iupper, var, values );
      
      /* Compute L2-norm of error at each time step and print to file */
      m = 0;
      for( j = 0; j < nly; j++ ){
         for( i = 0; i < nlx; i++ ){
            enorm += pow(values[m++] - 
                     exp(-2.*t)*sin((ilower[0]+i)*dx)*sin((ilower[1]+j)*dy),2);
         }
      }
      free( values );

      MPI_Comm_rank(app->comm_x, &myid_x);
      MPI_Reduce( &enorm, &enorm_sum, 1, MPI_DOUBLE, MPI_SUM, 0, app->comm_x );
      if(myid_x == 0)
      {  
         /* To complete the norm computation, multiply by dx and dy, and then take the sqrt */
         enorm_sum = sqrt(dx*dy*enorm_sum);
         
         /* Print max error at each time step */
         sprintf(filename, "%s.iter%03d.time%07d", "drive-05.error_norm", iter, index);
         file = fopen(filename, "w");
         fprintf(file, "%.14e\n", enorm_sum);
         fflush(file);
         fclose(file);
      }

   }

   /* Save the error and solution for GLVis visualization */
   if( app->output_vis && (level == 0) ){ 
      MPI_Comm_rank(app->comm_x, &myid_x);

      if( t == app->tstop ){
         values = (double *) malloc( nlx*nly*sizeof(double) );
         HYPRE_SStructVectorGetBoxValues( u->x, part, ilower, 
                                          iupper, var, values );

         HYPRE_SStructVectorCreate( app->comm_x, app->grid_x, &e );
         HYPRE_SStructVectorSetObjectType( e, app->object_type );
         HYPRE_SStructVectorInitialize( e );
      
         m = 0;
         for( j = 0; j < nly; j++ )
            for( i = 0; i < nlx; i++ )
            {
               values[m] = values[m] - exp(-2.*t)*sin((ilower[0]+i)*dx)*sin((ilower[1]+j)*dy);
               m++;
            }

         HYPRE_SStructVectorSetBoxValues( e, part, ilower,
                                          iupper, var, values );
         free( values );
         HYPRE_SStructVectorAssemble( e );

         sprintf(filename, "%s.iter%03d", "drive-05_mesh", iter);
         GLVis_PrintSStructGrid( app->grid_x, filename, 
                                 myid_x, NULL, NULL );
         sprintf(filename, "%s.iter%03d", "drive-05_err_tstop", iter);
         GLVis_PrintSStructVector( e, 0, filename, myid_x );
         sprintf(filename, "%s.iter%03d", "drive-05_sol_tstop", iter);
         GLVis_PrintSStructVector( u->x, 0, filename, myid_x );

         HYPRE_SStructVectorDestroy( e );
      }
   }

   return 0;
}


/* --------------------------------------------------------------------
 * Return buffer size for vector object buffer. Vector object contains
 * values at every grid point and thus, the buffer size is the number
 * of grid points.
 * -------------------------------------------------------------------- */
int
my_BufSize(braid_App  app,
           int       *size_ptr)
{
    /* A vector needs to contain 1 extra doubles for the coarsening rule */ 
    *size_ptr = (1 + (app->nlx)*(app->nly))*sizeof(double);
    return 0;
}


/* --------------------------------------------------------------------
 * Pack a vector object in a buffer.
 * -------------------------------------------------------------------- */
int
my_BufPack(braid_App     app,
           braid_Vector  u,
           void         *buffer,
           braid_Int    *size_ptr)
{
   double *dbuffer = buffer;
   int ilower[2], iupper[2], nlx, nly;
   
   /* Retrieve ilower and iupper */
   ilower[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower[0];
   ilower[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower[1];
   iupper[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper[0];
   iupper[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper[1];

   /* We have one variable and one part. */
   int        part = 0;
   int        var  = 0;
   
   /* Pack the spatial coarsening rule */ 
   dbuffer[0] = (double) u->spatial_disc_idx; 

   /* Pack the vector */
   HYPRE_SStructVectorGather( u->x );
   HYPRE_SStructVectorGetBoxValues( u->x, part, ilower,
                                    iupper, var, &(dbuffer[1]) );

   /* Determine number of bytes actually packed */
   nlx = iupper[0] - ilower[0] + 1;
   nly = iupper[1] - ilower[1] + 1;
   *size_ptr = (nlx*nly + 1)*sizeof(double);
   
   return 0;
}


/* --------------------------------------------------------------------
 * Unpack a vector object from a buffer.
 * -------------------------------------------------------------------- */
int
my_BufUnpack(braid_App     app,
             void         *buffer,
             braid_Vector *u_ptr)
{
   int ilower[2], iupper[2];
   double    *dbuffer = buffer;
   my_Vector *u;

   /* We have one variable and one part. */
   int        part = 0;
   int        var  = 0;

   u = (my_Vector *) malloc( sizeof(my_Vector) );

   /* Unpack the spatial coarsening rule */
   u->spatial_disc_idx = (int) dbuffer[0];

   /* Retrieve ilower and iupper */
   ilower[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower[0];
   ilower[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower[1];
   iupper[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper[0];
   iupper[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper[1];

   /* 
    * Unpack the vector 
    */

   /* Create an empty vector object. */
   HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[u->spatial_disc_idx]).grid_x_matrix, &(u->x) );
   
   /* Set the object type (by default HYPRE_SSTRUCT). */
   HYPRE_SStructVectorSetObjectType( u->x, app->object_type );

   /* Indicate that the vector coefficients are ready to be set. */
   HYPRE_SStructVectorInitialize( u->x );

   /* Set the values. */
   HYPRE_SStructVectorSetBoxValues( u->x, part, ilower,
                                    iupper, var, &(dbuffer[1]) );

   HYPRE_SStructVectorAssemble( u->x );

   *u_ptr = u;

   return 0;
}

/* Helper function for my_ComputeNumCoarsenings */
double log2( double n )  
{  
   /* log(n)/log(2) is log2. */  
   return log(n)/log(2.0);  
}

/* --------------------------------------------------------------------
 * Return the number of times to uniformly coarsen, in order
 * to minimally satisfy the CFL.  This function returns an integer beta, 
 * such that dx and dy must be multipled by 2^beta.   
 *
 * The CFL conditioner checked is:   
 *
 *       K*( dt/dx^2 + dt/dy^2 ) <  scoarsenCFL
 *
 * nx, ny are the global number of points in the x and y dimensions.  If the 
 * required number of coarsenings is too many for this grid, then zero 
 * is returned.
 * -------------------------------------------------------------------- */
int my_ComputeNumCoarsenings(double dt,
                             double dx,
                             double dy,
                             double K,
                             int nx,
                             int ny,
                             double scoarsenCFL)
{
   /* compute required coarsening to satisfy CFL */
   int ncoarsen1 = ceil( (log2( K*dt*(dx*dx + dy*dy)/(dx*dx*dy*dy) ) - log2(scoarsenCFL))/2.0 );
   int ncoarsen2 = ceil( (log2( K*0.999999999*dt*(dx*dx + dy*dy)/(dx*dx*dy*dy) ) - log2(scoarsenCFL))/2.0 );
   int ncoarsen3 = ceil( (log2( K*1.000000001*dt*(dx*dx + dy*dy)/(dx*dx*dy*dy) ) - log2(scoarsenCFL))/2.0 );
   
   /* Due to floating point arithmetic, this is our hack to make this the same
    * across processors.  dt is the only value that will vary across processors, so we try a couple
    * different tweaks up there on the order of the 10th digit. */
   int ncoarsen = max(ncoarsen1, ncoarsen2);
   ncoarsen = max(ncoarsen, ncoarsen3);

   int coarsen_factor =  (int) pow(2.0, ncoarsen);
   if( coarsen_factor == 0)
   {
      return 0;
   }

   /* Check if the existing grid can be coarsened that much */
   int cnx = (nx-1) / coarsen_factor  + 1; 
   int cny = (ny-1) / coarsen_factor  + 1; 
   if( (cnx < 2) || (cny < 2) )
   {
      return 0;
   }
   else
   {
      return ncoarsen;
   }
}

/* --------------------------------------------------------------------
 * Return the spatial discretization stored in  app->spatial_disc_table  
 * that corresponds to the (cdt, fdt) tuple.  Remember that this 
 * table is essentially a database with the tuples (cdt, fdt) as the keys.
 *
 * If spatial coarsening is required, then it is done by generating new SStruct
 * grids and graphs.
 *
 * The grids are numbered by their integer index into the spatial_disc_table
 * such that grid k is as fine or finer than grid k+1.  
 *
 * If no spatial coarsening is required to satisfy the CFL constraints when 
 * moving from grid k to grid k+1, then these two grids are the same.  
 * Otherwise they are different.  If multiple coarsenings are required to 
 * satisfy the CFl constraint when coarsening grid k, then these intermediate
 * grids are created as well (for use by the spatial interpolation/restriction
 * functions).  If two intermediate grids are needed, then grid k will spawn
 * the creation of grids k+1, k+2 and k+3 in the spatial_disc_table.  That is
 * spatial_disc_table[k+1] and spatial_disc_table[k+2] will hold grid information
 * for intermediate grids only used for spatial interpolation/restriction
 * and spatial_disc_table[k+3] will the grid information for the next XBraid level.
 *
 * If the entries cdt and fdt are both -1.0, then that denotes an empty table
 * entry.
 * 
 * Return value is 
 *    spatial_disc_idx
 * -------------------------------------------------------------------- */
void get_coarse_spatial_disc( braid_App app, 
                              double    cdt, 
                              double    fdt, 
                              double    fdx,
                              double    fdy,
                              double    scoarsenCFL,
                              int       fnx,
                              int       fny,
                              int       fspatial_disc_idx,
                              int*      filower,
                              int*      fiupper,
                              int*      spatial_disc_idx)
{
   int i, k, ncoarsen, ncoarsen_loc, cnlx, cnly, cnx, cny;
   int cilower[2], ciupper[2];
   int max_levels = app->max_levels;
   double coarsen_factor, cdx, cdy, cdt_loc;

   /* Initial spatial_disc_idx to -1, indicating that no match has yet 
    * been found */
   (*spatial_disc_idx) = -1;

   /* Search for this cdt, fdt combo and see if a spatial 
    * discretization already exists for it.  
    * Note: we start at i=1.  The first entry for the table is a dummy entry 
    *       for cloning vectors. 
    * Note: We assume that no more than 3*(app->max_levels) instances
    *       of spatial coarsening will ever occur.
    * */
   for(i = 1; i < 3*max_levels; i++)
   {
      if( fabs((app->spatial_disc_table[i]).cdt - cdt)/cdt < 1e-10 )
      {
         /* A rule for this cdt has already been stored */
         if( fabs((app->spatial_disc_table[i]).fdt - fdt)/fdt < 1e-10 )
         {
            /* We've found a match, both cdt and fdt match * 
             * This is a return value detailing which rule this is*/
            (*spatial_disc_idx) = i;
            return;
         }
      }
   }
   
   /* 
    * A spatial discretization was not found, so generate a new one, 
    * and store in the next open spot.  We loop over ncoarsen and generate
    * the intermediate grids as well (they are used for spatial 
    * interpolation/restriction).  Note that even if we don't coarsen, we 
    * still create a grid and graph for this cdt/fdt combo -- that is, the 
    * loop over k always iterates once.
    */
   
   /* Determine New sizes */
   ncoarsen = my_ComputeNumCoarsenings(cdt, fdx, fdy, app->K, fnx, fny, scoarsenCFL); 
   for(k = 1; k <= max(ncoarsen,1); k++)
   {
      if(ncoarsen == 0)
      {  
         coarsen_factor = 1.0;
         ncoarsen_loc = 0;
         cdt_loc = cdt;
      }
      else
      {  
         coarsen_factor = pow(2.0, k); 
         ncoarsen_loc = k;
         if(k == max(ncoarsen,1) )
         {  cdt_loc = cdt;}
         else
         {  cdt_loc = -1.0;}
      }

      cdx = coarsen_factor * fdx;
      cdy = coarsen_factor * fdy;
      cnx = (fnx-1) / ((int) coarsen_factor)  + 1; 
      cny = (fny-1) / ((int) coarsen_factor)  + 1;

      /* Define the nodes owned by the current processor (each processor's
       * piece of the global grid) */
      cilower[0] = ceil( (double) filower[0] / (double) coarsen_factor );
      cilower[1] = ceil( (double) filower[1] / (double) coarsen_factor );
      ciupper[0] = floor( (double) fiupper[0] / (double) coarsen_factor );
      ciupper[1] = floor( (double) fiupper[1] / (double) coarsen_factor );

      /* Determine local problem size. */
      cnlx = max(ciupper[0] - cilower[0] + 1, 0);
      cnly = max(ciupper[1] - cilower[1] + 1, 0);


      /* Note: We start at i=1.  The first entry for the table is a dummy
       *       entry for cloning vectors. 
       * Note: We assume that no more than 3*(app->max_levels) instances
       *       of spatial coarsening will ever occur.
       * */
      for(i = 1; i < 3*max_levels; i++)
      {
         if( ( (app->spatial_disc_table[i]).cdt == -1.0 ) && ( (app->spatial_disc_table[i]).fdt == -1.0 ) )
         {
            (app->spatial_disc_table[i]).cdt = cdt_loc;
            (app->spatial_disc_table[i]).fdt = fdt;
            (app->spatial_disc_table[i]).ncoarsen = ncoarsen_loc;
            (app->spatial_disc_table[i]).dx = cdx;
            (app->spatial_disc_table[i]).dy = cdy;
            (app->spatial_disc_table[i]).nx = cnx;
            (app->spatial_disc_table[i]).ny = cny;
            (app->spatial_disc_table[i]).nlx = cnlx;
            (app->spatial_disc_table[i]).nly = cnly;
            (app->spatial_disc_table[i]).ilower[0] = cilower[0];
            (app->spatial_disc_table[i]).ilower[1] = cilower[1];
            (app->spatial_disc_table[i]).iupper[0] = ciupper[0];
            (app->spatial_disc_table[i]).iupper[1] = ciupper[1];
            (app->spatial_disc_table[i]).fspatial_disc_idx = fspatial_disc_idx; 
            
            /* We also have to set up a new 2D grid and graph for the new COARSE level*/
            setUp2Dgrid( app->comm_x, &((app->spatial_disc_table[i]).grid_x_matrix), app->dim_x,  
                         cilower, ciupper, app->vartypes, 1 );
            setUpGraph( app->comm_x, &((app->spatial_disc_table[i]).graph_matrix), 
                       (app->spatial_disc_table[i]).grid_x_matrix, app->object_type, 
                        app->stencil );
            /* We temporarily let the vector grid and graph be the matrix grid and graph. Once we figure
             * out the correct ghost layer, we will replace these two values */
            (app->spatial_disc_table[i]).grid_x_vector = (app->spatial_disc_table[i]).grid_x_matrix;
            (app->spatial_disc_table[i]).graph_vector = (app->spatial_disc_table[i]).graph_matrix;
            
            /* We also have to set up a new 2D grid and graph for the fine level vectors that has the
             * correct number of ghost layers */
            setUp2Dgrid( app->comm_x, &((app->spatial_disc_table[fspatial_disc_idx]).grid_x_vector), 
                         app->dim_x, (app->spatial_disc_table[fspatial_disc_idx]).ilower, 
                         (app->spatial_disc_table[fspatial_disc_idx]).iupper, app->vartypes, 2);
            setUpGraph(app->comm_x, &((app->spatial_disc_table[fspatial_disc_idx]).graph_vector), 
                       (app->spatial_disc_table[fspatial_disc_idx]).grid_x_vector, 
                       app->object_type, app->stencil );

            /* This is a return value detailing which rule this is*/
            (*spatial_disc_idx) = i;
            break;
         }
      }
   
   } /* End ncoarsen loop */
   
   /* If no match was found, then (*spatial_disc_idx) will equal its initial value of -1 */
   return;
}

/* --------------------------------------------------------------------
 * Print a spatial vector (for diagnostics
 * -------------------------------------------------------------------- */
void print_spatial_vector(double * v,
                          int nx,
                          int ny,
                          MPI_Comm comm)
{
   int i,j,myid;
   MPI_Comm_rank(comm, &myid);
   fprintf(stderr, "\n\n");
   for(i = ny-1; i > -1; i--)
   {
      for(j = 0; j < nx; j++)
      {
         fprintf(stderr, "%11.2e", v[i*nx + j]);
      }
      fprintf(stderr, "\n");
   }
   fprintf(stderr, "\n\n");
}

/* --------------------------------------------------------------------
 * Retrieve a spatial discretization stored in app->spatial_disc_table  
 * that corresponds to the (cdt, fdt) tuple.  Remember that this 
 * table is essentially a database with the tuples (cdt, fdt) as the keys.*
 * -------------------------------------------------------------------- */
void retrieve_spatial_discretization( braid_App app, 
                                      double    cdt, 
                                      double    fdt,
                                      int       *spatial_disc_idx)
{
   int max_levels = app->max_levels;
   int i;
   
   /* Note: We start at i=1.  The first entry for the table is a dummy
    *       entry for cloning vectors.
    * Note: We assume that no more than 3*(app->max_levels) instances
    *       of spatial coarsening will ever occur.
    * */
   for(i = 1; i < 3*max_levels; i++)
   {
      if( fabs((app->spatial_disc_table[i]).cdt - cdt)/cdt < 1e-10 )
      {
         /* A rule for this cdt has already been stored */
         if( fabs((app->spatial_disc_table[i]).fdt - fdt)/fdt < 1e-10 )
         {
            /* We've found a match, both cdt and fdt match */
            (*spatial_disc_idx) = i;
            return;
         }
      }
   }
   (*spatial_disc_idx) = -1;
}


/* --------------------------------------------------------------------
 * Do a single uniform refinement of a braid_Vector.
 * Use bilinear interpolation.
 * Assume a regular grid of size 2^k + 1 in each dimension.
 *
 * The basic strategy is to (possible multiple) uniform refinments with a
 * helper function that does only one uniform refinement while an outter 
 * driver calls the helper function ncoarsen number of times.
 * -------------------------------------------------------------------- */
int
my_RefineHelper(braid_App              app,           
                braid_Vector           cu,
                braid_Vector          *fu_ptr,
                braid_CoarsenRefStatus status)
{
   my_Vector  *fu;
   double     *cvalues, *fvalues, *fvaluesplus;
   int        cidx, fidx, i, j, cstride, fstride, istart, istop, jstart, jstop;
   int        filower[2], fiupper[2];
   int        spatial_disc_idx, fspatial_disc_idx;
   int        fnlx, fnly;

   int        cilower[2], ciupper[2];
   int        cnlx, cnly;
   
   hypre_Box             *box;
   hypre_Box             *value_box;
   hypre_CommInfo        *comm_info;
   hypre_CommHandle      *comm_handle;
   hypre_CommPkg         *comm_pkg;
   int                   *num_ghost;
   
   /*double     tstart;
   braid_CoarsenRefStatusGetTstart(status, &tstart); */
   
   /* Determine current local spatial grid size.  Note that subtracting one from a 
    * spatial discretization index will give you the index for the next finer grid. */
   spatial_disc_idx  = cu->spatial_disc_idx; 
   fspatial_disc_idx = spatial_disc_idx - 1;
   cnlx              = (app->spatial_disc_table[spatial_disc_idx]).nlx;
   cnly              = (app->spatial_disc_table[spatial_disc_idx]).nly;
   cilower[0]        = (app->spatial_disc_table[spatial_disc_idx]).ilower[0];
   cilower[1]        = (app->spatial_disc_table[spatial_disc_idx]).ilower[1];
   ciupper[0]        = (app->spatial_disc_table[spatial_disc_idx]).iupper[0];
   ciupper[1]        = (app->spatial_disc_table[spatial_disc_idx]).iupper[1];

   /* Find this processor's part of the fine vector.  If there is no
    * coarsening, then fu will be the same size as cu. */    
   fnlx              = (app->spatial_disc_table[fspatial_disc_idx]).nlx;
   fnly              = (app->spatial_disc_table[fspatial_disc_idx]).nly;
   filower[0]        = (app->spatial_disc_table[fspatial_disc_idx]).ilower[0];
   filower[1]        = (app->spatial_disc_table[fspatial_disc_idx]).ilower[1];
   fiupper[0]        = (app->spatial_disc_table[fspatial_disc_idx]).iupper[0];
   fiupper[1]        = (app->spatial_disc_table[fspatial_disc_idx]).iupper[1];


   /* Now, begin the long road to interpolating cu to fu */ 
   if( (fiupper[0] >= filower[0]) && (fiupper[1] >= filower[1])) 
   {
       fvalues = (double *) malloc( (fnlx)*(fnly)*sizeof(double) );
       for(i=0; i < fnlx*fnly; i++)
       {  fvalues[i] = 0.0;}
   }

   /* Grab values from cu and inject into fvalues, if no spatial coarsening occurs, this
    * will just copy cu into fu */
   if( (ciupper[0] >= cilower[0]) && (ciupper[1] >= cilower[1]) )
   {
      HYPRE_SStructVectorGather( cu->x );
      cvalues = (double *) malloc( (cnlx)*(cnly)*sizeof(double) );
      HYPRE_SStructVectorGetBoxValues( cu->x, 0, cilower,
                                       ciupper, 0, cvalues );
      
      cstride = 2;
      cidx = 0;
      for(i = 0; i < cnly; i++)
      {   
         for(j = 0; j < cnlx; j++)
         {
            /* this maps i,j (the local c row and column numbers) to 
             * local f row and column numbers */
            fidx = ((i + cilower[1])*cstride - filower[1])*fnlx +
                   ((j + cilower[0])*cstride - filower[0]);
            fvalues[fidx] = cvalues[cidx];
            cidx ++;
         }
      }
      
      free(cvalues); 
   }
   
   /* Setup fvector -- Be sure to use the grid with a larger ghost zone !!*/
   fu = (my_Vector *) malloc(sizeof(my_Vector));
   HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[fspatial_disc_idx]).grid_x_vector, &(fu->x) );
   HYPRE_SStructVectorSetObjectType( fu->x, app->object_type );
   HYPRE_SStructVectorInitialize( fu->x );
   if( (fiupper[0] >= filower[0]) && (fiupper[1] >= filower[1])) 
   {  HYPRE_SStructVectorSetBoxValues( fu->x, 0, filower,
                                    fiupper, 0, fvalues ); }
   HYPRE_SStructVectorAssemble( fu->x );
   fu->spatial_disc_idx = fspatial_disc_idx; 
   if( (fiupper[0] >= filower[0]) && (fiupper[1] >= filower[1])) 
   { free(fvalues); } 

   /* Communicate ghost layer on fine level with fu and interpolate. 
    * Only enter this code block if ncoarsen < 0, that is if you're actually going 
    * to do interpolation.  Otherwise, you're done -- you've cloned cu into fu, 
    * which are of the same size, and done a SetBoxValues with fvalues */
   if( (fiupper[0] >= filower[0]) && (fiupper[1] >= filower[1])) 
   {
       HYPRE_SStructVectorGather( fu->x );
       hypre_StructVector *fu_struct = hypre_SStructPVectorSVector( (hypre_SStructVectorPVector((fu->x),0)), 0);
       num_ghost = hypre_StructVectorNumGhost(fu_struct);
       hypre_CreateCommInfoFromNumGhost(hypre_StructVectorGrid(fu_struct),
                                        num_ghost, &comm_info);
       hypre_CommPkgCreate(comm_info,
                           hypre_StructVectorDataSpace(fu_struct),
                           hypre_StructVectorDataSpace(fu_struct),
                           1, NULL, 0,
                           hypre_StructVectorComm(fu_struct), &comm_pkg);
       hypre_CommInfoDestroy(comm_info);
       hypre_InitializeCommunication( comm_pkg,
                                      hypre_StructVectorData(fu_struct),
                                      hypre_StructVectorData(fu_struct),
                                      0, 0,
                                      &comm_handle );
       hypre_FinalizeCommunication( comm_handle );
       box = hypre_BoxCreate(2);
       filower[0] -= num_ghost[0];
       filower[1] -= num_ghost[0];
       fiupper[0] += num_ghost[0];
       fiupper[1] += num_ghost[0];
       fnlx += 2*num_ghost[0];
       fnly += 2*num_ghost[0];
       hypre_CopyIndex(filower, hypre_BoxIMin(box));
       hypre_CopyIndex(fiupper, hypre_BoxIMax(box));
       value_box = box;
       fvaluesplus = (double *) malloc( (fnlx)*(fnly)*sizeof(double) );
       hypre_StructVectorSetBoxValues( fu_struct,
                                       box,
                                       value_box,
                                       fvaluesplus,    /* HYPRE_Complex      *values, */
                                       -1,             /* HYPRE_Int           action, */
                                       0,              /* HYPRE_Int           boxnum, */
                                       1);             /* HYPRE_Int           outside */ 
       hypre_BoxDestroy(box);
       hypre_CommPkgDestroy(comm_pkg);
       
       /* print the vectors for some time step (comment in matching print below)
        * $$ srun -N 1 -n 1 -p pdebug ./drive-05 -pgrid 1 1 1 -nt 256 -mi 2 -ml 15 -nx 9 9 -scoarsen 2 
        * 
       if( (tstart > (0.00578197133*1000.)) &&  (tstart < (0.00578397133*1001.)) )
       {    print_spatial_vector(fvaluesplus, fnlx, fnly, app->comm_x); } */ 

       /* So, now we have our ghost layer information in fvaluesplus.  Time to interpolate.
        * fvaluesplus will be nonzero at C-pts and zero at F-pts */
       cstride = 2;
       fstride = 1;
          
       /* Begin Interpolation 
        * (1) While this may be a patch of the global grid, locally we still have a 
        *     numbering that starts at (0,0) and ends at (fnlx, fnly). 
        * (2) Skip interpolation for the ghost layer (not needed), that is start and 
        *     end loops at +1, -1
        * (3) i, j are global indices below, so we can easily check if they are C/F points
        * */
       
       /* Compute the loop starts and stops
        * istart, istop, jstart and jstop are the global indices AT THIS LEVEL OF COARSENING
        * for the first and last rows and columns (i.e. box coordinates) local to this proc. */
       
       /*                     Note this integer div is a floor!           */
       jstart = cilower[0]*cstride - (abs(cilower[0]*cstride - (filower[0] + num_ghost[0]))/cstride) * cstride;
       if( (jstart - fstride) >= (filower[0] + num_ghost[0]) )
       { jstart -= fstride; }
       /* */
       jstop  = ciupper[0]*cstride + (abs(ciupper[0]*cstride - (fiupper[0] - num_ghost[0]))/cstride) * cstride;
       if( (jstop + fstride) <= (fiupper[0] - num_ghost[0]) )
       { jstop += fstride; }
       /* */
       istart = cilower[1]*cstride - (abs(cilower[1]*cstride - (filower[1] + num_ghost[0]))/cstride) * cstride;
       if( (istart - fstride) >= (filower[1] + num_ghost[0] ) )
       { istart -= fstride; }
       /* */
       istop  = ciupper[1]*cstride + (abs(ciupper[1]*cstride - (fiupper[1] - num_ghost[0]))/cstride) * cstride;
       if( (istop + fstride) <= (fiupper[1] - num_ghost[0]) )
       { istop += fstride; }

       /* Loop over rows */
       for(i = istart; i <= istop; i+= fstride)
       {
          /* fidx is the linear index into fvaluesplus for row i
           * --> have to make sure to skip the ghost layer of points 
           *     by accounting for filower*/
          fidx = (i - filower[1])*fnlx + (jstart - filower[0]);

          /* Loop over columns */
          for(j = jstart; j <= jstop; j+= fstride)
          {
             /* No rescaling needed at boundary points 
              * We either inject to C-points at the boundary or 
              * interpolate with weights = 1/2 from two neighboring 
              * C-points on the boundary -- that is, there are no 
              * interpolation stencils that don't preserve the constant 
              * */

             /* if( (i%cstride== 0) && (j%cstride== 0) )
             {
                 This is a C-point and has already been injected
             } */
             if( (i%cstride == 0) && (j%cstride == fstride) )
             {
                /* Interpolate: this is an F-point horizontally between two C-points */
                fvaluesplus[fidx] = 0.5*fvaluesplus[ fidx - fstride ] + 
                                    0.5*fvaluesplus[ fidx + fstride ];
             }
             else if( (i%cstride == fstride ) && (j%cstride == 0) )
             {
                /* Interpolate: this is an F-point vertically between two C-points */
                fvaluesplus[fidx] = 0.5*fvaluesplus[ fidx + fnlx*fstride] + 
                                    0.5*fvaluesplus[ fidx - fnlx*fstride];
             }
             else if( (i%cstride == fstride ) && (j%cstride == fstride ) )
             {
                /* Interpolate: this is an F-point in the center of a grid cell */
                fvaluesplus[fidx] = 0.25*fvaluesplus[ fidx + fnlx*fstride + fstride ] + 
                                    0.25*fvaluesplus[ fidx + fnlx*fstride - fstride ] +
                                    0.25*fvaluesplus[ fidx - fnlx*fstride + fstride ] + 
                                    0.25*fvaluesplus[ fidx - fnlx*fstride - fstride ];
             }
             fidx+=fstride;
          }
       }
       
       /*if( (tstart > (0.00578197133*1000.)) &&  (tstart < (0.00578397133*1001.)) )
       {    print_spatial_vector(fvaluesplus, fnlx, fnly, app->comm_x); } */ 
   }

   /* Do a new SetBoxValues with interpolated F-point values,
    * Be sure to use the grid with the reduced ghost zone !! */
   HYPRE_SStructVectorDestroy( fu->x );
   free( fu );
   fu = (my_Vector *) malloc(sizeof(my_Vector));
   HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[fspatial_disc_idx]).grid_x_matrix, &(fu->x) );
   HYPRE_SStructVectorSetObjectType( fu->x, app->object_type );
   HYPRE_SStructVectorInitialize( fu->x );
   if( (fiupper[0] >= filower[0]) && (fiupper[1] >= filower[1])) 
   {   HYPRE_SStructVectorSetBoxValues( fu->x, 0, filower,
                                    fiupper, 0, fvaluesplus ); }
   HYPRE_SStructVectorAssemble( fu->x );
   fu->spatial_disc_idx = fspatial_disc_idx; 
       
   if( (fiupper[0] >= filower[0]) && (fiupper[1] >= filower[1])) 
   {   free(fvaluesplus); } 

   /* Set return value */
   *fu_ptr = fu;
   
   return 0;
}

/* Driver funtion that calls my_RefineHelper to carry out possibly multiple
 * uniform refinments */
int
my_Refine(braid_App              app,           
          braid_Vector           cu,
          braid_Vector          *fu_ptr,
          braid_CoarsenRefStatus status)
{
   int        k, ncoarsen, spatial_disc_idx;
   double     cdt, fdt;
   double     tstart, f_tstop, f_tprior, c_tstop, c_tprior;

   
   /* Get coarse and fine time step sizes */
   braid_CoarsenRefStatusGetTstart(status, &tstart);
   braid_CoarsenRefStatusGetCTstop(status, &c_tstop);
   braid_CoarsenRefStatusGetCTprior(status, &c_tprior);
   braid_CoarsenRefStatusGetFTstop(status, &f_tstop);
   braid_CoarsenRefStatusGetFTprior(status, &f_tprior);
   cdt = c_tstop - tstart;
   fdt = f_tstop - tstart;

   /* If fdt or cdt is 0.0,  then this is a final time interval.  We then use a
    * simple rule to use the length of the previous time interval to represent
    * an appropriate dt. */
   if (fdt == 0.0)
   {
      fdt = tstart - f_tprior;
   }
   if (cdt == 0.0)
   {
      cdt = tstart - c_tprior;
   }

   /* Determine current local spatial grid size and negate ncoarsen because we
    * are refining, i.e., multiplying dx and dy by 2^(-ncoarsen) */
   retrieve_spatial_discretization(app, cdt, fdt, &spatial_disc_idx );
   ncoarsen          = (app->spatial_disc_table[spatial_disc_idx]).ncoarsen;
   
   /* If no coarsening, then just clone.
    * Otherwise, repeatedly coarsen by a factor of two until the 
    * desired grid is reached */
   if(ncoarsen == 0)
   {
     /* Unless multiple coarsenings/refinements happen between levels, then 
      * the spatial_disc_idx just progresses by +/- 1*/   
      my_Clone(app, cu, fu_ptr);
      (*fu_ptr)->spatial_disc_idx = spatial_disc_idx-1;
   }
   else
   {
      for(k = 0; k < ncoarsen; k++)
      {     
         my_RefineHelper(app, cu, fu_ptr, status);
         if((k > 0) && (k < (ncoarsen-1)) )
         {  my_Free(app, cu); }
         cu = *fu_ptr;
      }
   }

   return 0;
}

/* --------------------------------------------------------------------
 * Do a single uniform coarsening of a braid_Vector.
 * Use (1/4) times the transpose of bilinear interpolation.
 * Assume a regular grid of size 2^k + 1 in each dimension.
 *
 * The basic strategy is to (possible multiple) uniform coaresnings with a
 * helper function that does only one uniform coarsening while an outter 
 * driver calls the helper function ncoarsen number of times.
 * -------------------------------------------------------------------- */

int
my_CoarsenBilinearHelper(braid_App              app,           
                         braid_Vector           fu,
                         braid_Vector          *cu_ptr,
                         braid_CoarsenRefStatus status)
{
   my_Vector  *cu;
   double     *cvalues, *fvalues;
   int        counter, i, j, fidx, cstride, fstride;
   int        filower[2], fiupper[2];
   int        fspatial_disc_idx = fu->spatial_disc_idx;
   int        fnlx = (app->spatial_disc_table[fspatial_disc_idx]).nlx;
   int        fnly = (app->spatial_disc_table[fspatial_disc_idx]).nly;
   int        fnx  = (app->spatial_disc_table[fspatial_disc_idx]).nx;
   int        fny  = (app->spatial_disc_table[fspatial_disc_idx]).ny;

   int        cilower[2], ciupper[2];
   int        cnlx, cnly;
   double     one_eigth = 1.0 / 8.0;
   double     one_sixteenth = 1.0 / 16.0;
   double     scale, tstart;
   
   hypre_Box             *box;
   hypre_Box             *value_box;
   hypre_CommInfo        *comm_info;
   hypre_CommHandle      *comm_handle;
   hypre_CommPkg         *comm_pkg;
   int                   *num_ghost;

   braid_CoarsenRefStatusGetTstart(status, &tstart);
   cu = (my_Vector *) malloc(sizeof(my_Vector));
   
   filower[0]    = (app->spatial_disc_table[fspatial_disc_idx]).ilower[0];
   filower[1]    = (app->spatial_disc_table[fspatial_disc_idx]).ilower[1];
   fiupper[0]    = (app->spatial_disc_table[fspatial_disc_idx]).iupper[0];
   fiupper[1]    = (app->spatial_disc_table[fspatial_disc_idx]).iupper[1];
   
   /* Get the next coarser spatial discretization */ 
   cnlx     = (app->spatial_disc_table[fspatial_disc_idx+1]).nlx;
   cnly     = (app->spatial_disc_table[fspatial_disc_idx+1]).nly;
   cilower[0] = (app->spatial_disc_table[fspatial_disc_idx+1]).ilower[0];
   cilower[1] = (app->spatial_disc_table[fspatial_disc_idx+1]).ilower[1];
   ciupper[0] = (app->spatial_disc_table[fspatial_disc_idx+1]).iupper[0];
   ciupper[1] = (app->spatial_disc_table[fspatial_disc_idx+1]).iupper[1];

   /* Create an empty vector object. */
   HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[fspatial_disc_idx+1]).grid_x_matrix, &(cu->x) );
   
   /* Set the object type (by default HYPRE_SSTRUCT) and initialize. */
   HYPRE_SStructVectorSetObjectType( cu->x, app->object_type );
   HYPRE_SStructVectorInitialize( cu->x );

   /* Only do work if you have a non-empty fine grid */
   if( (filower[0] <= fiupper[0]) && (filower[1] <= fiupper[1]) )
   {

      /* Grab the fvalues plus the default single ghost layer  */
      HYPRE_SStructVectorGather( fu->x );
      hypre_StructVector *fu_struct = hypre_SStructPVectorSVector( (hypre_SStructVectorPVector((fu->x),0)), 0);
      num_ghost = hypre_StructVectorNumGhost(fu_struct);
      hypre_CreateCommInfoFromNumGhost(hypre_StructVectorGrid(fu_struct),
                                       num_ghost, &comm_info);
      hypre_CommPkgCreate(comm_info,
                          hypre_StructVectorDataSpace(fu_struct),
                          hypre_StructVectorDataSpace(fu_struct),
                          1, NULL, 0,
                          hypre_StructVectorComm(fu_struct), &comm_pkg);
      hypre_CommInfoDestroy(comm_info);
      hypre_InitializeCommunication( comm_pkg,
                                     hypre_StructVectorData(fu_struct),
                                     hypre_StructVectorData(fu_struct),
                                     0, 0,
                                     &comm_handle );
      hypre_FinalizeCommunication( comm_handle );
      box = hypre_BoxCreate(2);
      filower[0] -= num_ghost[0];
      filower[1] -= num_ghost[0];
      fiupper[0] += num_ghost[0];
      fiupper[1] += num_ghost[0];
      fnlx += 2*num_ghost[0];
      fnly += 2*num_ghost[0];
      hypre_CopyIndex(filower, hypre_BoxIMin(box));
      hypre_CopyIndex(fiupper, hypre_BoxIMax(box));
      value_box = box;
      fvalues = (double *) malloc( (fnlx)*(fnly)*sizeof(double) );
      hypre_StructVectorSetBoxValues( fu_struct,
                                      box,
                                      value_box,
                                      fvalues,    /* HYPRE_Complex      *values, */
                                      -1,             /* HYPRE_Int           action, */
                                      0,              /* HYPRE_Int           boxnum, */
                                      1);             /* HYPRE_Int           outside */ 
      hypre_BoxDestroy(box);
      hypre_CommPkgDestroy(comm_pkg);
     
      /* print the vectors for some time step (comment in matching print below)
       * $$ srun -N 1 -n 1 -p pdebug ./drive-05 -pgrid 1 1 1 -nt 256 -mi 2 -ml 15 -nx 9 9 -scoarsen 2 
       * 
      if( (tstart > (0.00578197133*1000.)) &&  (tstart < (0.00578397133*1001.)) )
      {    print_spatial_vector(fvalues, fnlx, fnly, app->comm_x); } */

      /* Only do if you have a nonempty c-vector to fill in */
      if( (cilower[0] <= ciupper[0]) && (cilower[1] <= ciupper[1]) )
      {
         /* Begin Restriction Loop (1/4 times bilinear interp)^T
          * Note that (i,j) are global indices below, so we can easily check if they are C/F points
          **/ 
         cstride = 2;
         fstride = 1;
         
         /* Loop over rows, skipping the ghost layer
          * Here, we just initialize the C-points for this coarsening level */
         for(i = filower[1]+num_ghost[0]; i <= fiupper[1]-num_ghost[0]; i++)
         {
            /* fidx is the linear index into fvalues for this (i,j) point.  It is 
             * the offset for skipping i-rows into fvalues.  You have to make sure 
             * to skip the ghost layer of points */
            fidx = (i - filower[1])*fnlx + num_ghost[0];

            /* Loop over columns */
            for(j = filower[0]+num_ghost[0]; j <= fiupper[0]-num_ghost[0]; j++)
            {
               if( (i%cstride == 0) && (j%cstride == 0) )
               {
                   /* Have to rescale at boundary points */
                   scale = 1.0;
                   if( (i==0) || (i==fnx-1) )
                   {
                      if( (j==0) || (j==fny-1) )
                      {  scale *= (16./9.); }
                      else
                      {  scale *= (4./3.); }
                   }
                   else if( (j==0) || (j==fny-1) )
                   {  scale *= (4./3.); }

                   /* This is a C-point */
                   fvalues[fidx] = scale*0.25*fvalues[ fidx ];
               
                   /* Interpolate from F-points horizontally to the left and right */
                   fvalues[fidx] += scale*one_eigth*(fvalues[ fidx - fstride ] + fvalues[ fidx + fstride ]); 
                   
                   /* Interpolate from F-points vertically above and below */
                   fvalues[fidx] += scale*one_eigth*(fvalues[ fidx + fnlx] + fvalues[ fidx - fnlx]);
                   
                   /* Interpolate from diagonally offset points */
                   fvalues[fidx] += scale*one_sixteenth*(fvalues[ fidx  + fnlx + fstride] +
                     fvalues[ fidx  + fnlx - fstride] + fvalues[ fidx  - fnlx + fstride] +fvalues[ fidx  - fnlx - fstride]);
               }
               fidx++;
            }
         }
            

         /* Copy over the coarse values from fvalues into cvalues, update cu and exit */
         counter = 0;
         cvalues = (double *) malloc( (cnlx)*(cnly)*sizeof(double) );
         for(i = filower[1]+num_ghost[0]; i <= fiupper[1]-num_ghost[0]; i++)
         {
            /* fidx is the linear index into fvalues for this (i,j) point.  It is 
             * the offset for skipping i-rows into fvalues.  You have to make sure 
             * to skip the ghost layer of points */
            fidx = (i - filower[1])*fnlx + num_ghost[0];

            /* Loop over columns */
            for(j = filower[0]+num_ghost[0]; j <= fiupper[0]-num_ghost[0]; j++)
            {
               
               if( (i%cstride == 0) && (j%cstride == 0) )
               {
                   /* This is a C-point */
                   cvalues[counter] = fvalues[ fidx ];
                   counter++;
               }
               fidx++;
            }
         }
   
         HYPRE_SStructVectorSetBoxValues( cu->x, 0, cilower,
                                          ciupper, 0, cvalues );
      
         /* Make sure to comment in matching print above 
         if( (tstart > (0.00578197133*1000.)) &&  (tstart < (0.00578397133*1001.)) )
         {    print_spatial_vector(cvalues, cnlx, cnly, app->comm_x); } */

         free( cvalues );
      } /* End-if for the non-empty f-level box */

      free( fvalues );
   } /* End-if for the non-empty c-level box */
   
   HYPRE_SStructVectorAssemble( cu->x );
     

   /* Unless multiple coarsenings/refinements happen between levels, then 
    * the spatial_disc_idx just progresses by +/- 1*/   
   cu->spatial_disc_idx = fspatial_disc_idx+1;
   *cu_ptr = cu;

   return 0;
}


/* Driver funtion that calls my_CoarsenBilinearHelper to carry out possibly 
 * multiple uniform coarsenings*/
int
my_CoarsenBilinear(braid_App              app,           
                   braid_Vector           fu,
                   braid_Vector          *cu_ptr,
                   braid_CoarsenRefStatus status)
{
   int        spatial_disc_idx;
   int        filower[2], fiupper[2];
   int        fspatial_disc_idx = fu->spatial_disc_idx;
   int        fnx = (app->spatial_disc_table[fspatial_disc_idx]).nx;
   int        fny = (app->spatial_disc_table[fspatial_disc_idx]).ny;
   double     fdy = (app->spatial_disc_table[fspatial_disc_idx]).dy;
   double     fdx = (app->spatial_disc_table[fspatial_disc_idx]).dx;

   int        k, ncoarsen, level;
   double     cdt, fdt, scoarsenCFL;
   double     tstart, f_tstop, f_tprior, c_tstop, c_tprior;
   
   /* Get Coarse and fine time step sizes */
   braid_CoarsenRefStatusGetTstart(status, &tstart);
   braid_CoarsenRefStatusGetCTstop(status, &c_tstop);
   braid_CoarsenRefStatusGetCTprior(status, &c_tprior);
   braid_CoarsenRefStatusGetFTstop(status, &f_tstop);
   braid_CoarsenRefStatusGetFTprior(status, &f_tprior);
   cdt = c_tstop - tstart;
   fdt = f_tstop - tstart;
   
   filower[0]    = (app->spatial_disc_table[fspatial_disc_idx]).ilower[0];
   filower[1]    = (app->spatial_disc_table[fspatial_disc_idx]).ilower[1];
   fiupper[0]    = (app->spatial_disc_table[fspatial_disc_idx]).iupper[0];
   fiupper[1]    = (app->spatial_disc_table[fspatial_disc_idx]).iupper[1];
   
   /* If fdt or cdt is 0.0,  then this is a final time interval.  We then use a
    * simple rule to use the length of the previous time interval to represent
    * an appropriate dt. */
   if (fdt == 0.0)
   {
      fdt = tstart - f_tprior;
   }
   if (cdt == 0.0)
   {
      cdt = tstart - c_tprior;
   }
   
   /* Generate the next spatial discretization, which is stored in app->spatial_disc_table[i]
    * This could be the same as the fine spatial discretization */ 
   braid_CoarsenRefStatusGetLevel(status, &level);
   scoarsenCFL = app->scoarsenCFL[min(app->num_scoarsenCFL-1, level)];
   get_coarse_spatial_disc( app, cdt, fdt, fdx, fdy, scoarsenCFL, fnx, fny, 
                            fspatial_disc_idx, filower, fiupper, &spatial_disc_idx);
   ncoarsen = (app->spatial_disc_table[spatial_disc_idx]).ncoarsen;
         
   /* If no coarsening, then just clone.
    * Otherwise, repeatedly coarsen by a factor of two until the 
    * desired grid is reached */
   if(ncoarsen == 0)
   {
      /* Unless multiple coarsenings/refinements happen between levels, then 
       * the spatial_disc_idx just progresses by +/- 1*/
      my_Clone(app, fu, cu_ptr);
      (*cu_ptr)->spatial_disc_idx = fspatial_disc_idx+1;
   }
   else
   {
      for(k = 0; k < ncoarsen; k++)
      {     
         my_CoarsenBilinearHelper(app, fu, cu_ptr, status);
         if((k > 0) && (k < (ncoarsen-1)) )
         {  my_Free(app, fu); }
         fu = *cu_ptr;
         
      
      }
   }
   

   return 0;
}


/* --------------------------------------------------------------------
 * Create a coarsened copy of a vector object.
 * Only uses injection
 * Assume a regular grid of size 2^k + 1 in each dimension
 * -------------------------------------------------------------------- */
int
my_CoarsenInjection(braid_App              app,           
                    braid_Vector           fu,
                    braid_Vector          *cu_ptr,
                    braid_CoarsenRefStatus status)
{
   my_Vector  *cu;
   double     *cvalues, *fvalues;
   int        counter, flrow, flcol, i, j, spatial_disc_idx;
   int        filower[2], fiupper[2];
   int        fspatial_disc_idx = fu->spatial_disc_idx;
   int        fnx = (app->spatial_disc_table[fspatial_disc_idx]).nx;
   int        fny = (app->spatial_disc_table[fspatial_disc_idx]).ny;
   int        fnlx = (app->spatial_disc_table[fspatial_disc_idx]).nlx;
   int        fnly = (app->spatial_disc_table[fspatial_disc_idx]).nly;
   double     fdy = (app->spatial_disc_table[fspatial_disc_idx]).dy;
   double     fdx = (app->spatial_disc_table[fspatial_disc_idx]).dx;

   int        cilower[2], ciupper[2];
   double     coarsen_factor;
   int        cnlx, cnly, ncoarsen, coarsen_factor_int, level;
   double     cdt, fdt, scoarsenCFL;
   double     tstart, f_tstop, f_tprior, c_tstop, c_tprior;

   /* Get Coarse and fine time step sizes */
   braid_CoarsenRefStatusGetTstart(status, &tstart);
   braid_CoarsenRefStatusGetCTstop(status, &c_tstop);
   braid_CoarsenRefStatusGetCTprior(status, &c_tprior);
   braid_CoarsenRefStatusGetFTstop(status, &f_tstop);
   braid_CoarsenRefStatusGetFTprior(status, &f_tprior);
   cdt = c_tstop - tstart;
   fdt = f_tstop - tstart;
   
   cu = (my_Vector *) malloc(sizeof(my_Vector));
   
   filower[0]    = (app->spatial_disc_table[fspatial_disc_idx]).ilower[0];
   filower[1]    = (app->spatial_disc_table[fspatial_disc_idx]).ilower[1];
   fiupper[0]    = (app->spatial_disc_table[fspatial_disc_idx]).iupper[0];
   fiupper[1]    = (app->spatial_disc_table[fspatial_disc_idx]).iupper[1];
   
   /* If fdt or cdt is 0.0,  then this is a final time interval.  We then use a
    * simple rule to use the length of the previous time interval to represent
    * an appropriate dt. */
   if (fdt == 0.0)
   {
      fdt = tstart - f_tprior;
   }
   if (cdt == 0.0)
   {
      cdt = tstart - c_tprior;
   }
   
   /* Generate the next spatial discretization, which is stored in app->spatial_disc_table[i]
    * This could be the same as the fine spatial discretization */ 
   braid_CoarsenRefStatusGetLevel(status, &level);
   scoarsenCFL = app->scoarsenCFL[min(app->num_scoarsenCFL-1, level)];
   get_coarse_spatial_disc( app, cdt, fdt, fdx, fdy, scoarsenCFL, fnx, fny, 
                            fspatial_disc_idx, filower, fiupper, &spatial_disc_idx);
   ncoarsen = (app->spatial_disc_table[spatial_disc_idx]).ncoarsen;
   cnlx     = (app->spatial_disc_table[spatial_disc_idx]).nlx;
   cnly     = (app->spatial_disc_table[spatial_disc_idx]).nly;
   cilower[0] = (app->spatial_disc_table[spatial_disc_idx]).ilower[0];
   cilower[1] = (app->spatial_disc_table[spatial_disc_idx]).ilower[1];
   ciupper[0] = (app->spatial_disc_table[spatial_disc_idx]).iupper[0];
   ciupper[1] = (app->spatial_disc_table[spatial_disc_idx]).iupper[1];
   coarsen_factor = pow(2.0, ncoarsen);
   coarsen_factor_int = (int) coarsen_factor;

   /* Create an empty vector object. */
   HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[spatial_disc_idx]).grid_x_matrix, &(cu->x) );
   
   /* Set the object type (by default HYPRE_SSTRUCT) and initialize. */
   HYPRE_SStructVectorSetObjectType( cu->x, app->object_type );
   HYPRE_SStructVectorInitialize( cu->x );

   /* Set the coarse values, assuming a regular grid and injection.  Just 
    * loop over the f points, and when the fpoint index tuple is (even,even)
    * then you know that you have a C-point.  Do this lexicographically,
    * and that's it.  It's really simple. */
   cvalues = (double *) malloc( (cnlx)*(cnly)*sizeof(double) );
   HYPRE_SStructVectorGather( fu->x );
   fvalues = (double *) malloc( fnlx*fnly*sizeof(double) );
   HYPRE_SStructVectorGetBoxValues( fu->x, 0, filower, fiupper, 0, fvalues );
   counter = 0;
   flrow = 0;
   for(i = filower[1]; i <= fiupper[1]; i++)
   {
      /* flrow and flcol are named for f-local-row and f-local-column, where
       * the numbering starts at 0. */
      flcol = 0;
      for(j = filower[0]; j <= fiupper[0]; j++)
      {
         /* Even indexed points are always C-points, because the grids are 
          * for standard coarsening and sizes of 2^k+1 */
         if( (i%coarsen_factor_int == 0) && (j%coarsen_factor_int == 0) )
         {
            cvalues[counter] = fvalues[flrow*fnlx + flcol];
            counter += 1;
         }
         flcol += 1;
      }
      
      flrow += 1;
   }
   
   HYPRE_SStructVectorSetBoxValues( cu->x, 0, cilower,
                                    ciupper, 0, cvalues );
   HYPRE_SStructVectorAssemble( cu->x );
   
   /* print the vectors for some time step */
   /*if( (tstart > (0.00578197133*1000.)) &&  (tstart < (0.00578397133*1001.)) )
   {    print_spatial_vector(fvalues, fnlx, fnly, app->comm_x); }
   if( (tstart > (0.00578197133*1000.)) &&  (tstart < (0.00578397133*1001.)) )
   {    print_spatial_vector(cvalues, cnlx, cnly, app->comm_x); } */
   
   free( fvalues );
   free( cvalues );

   /* Store spatial_disc_idx */
   cu->spatial_disc_idx = spatial_disc_idx;

   *cu_ptr = cu;

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

   int correct, fspatial_disc_idx;

   braid_Core    core;
   my_App       *app;
   int           max_levels;
   int           min_coarse;
   int           nrelax, nrelax0;
   double        tol;
   int           cfactor, cfactor0;
   int           max_iter;
   int           fmg;
   int           tnorm;
   int           nfmg_Vcyc;
   int           scoarsen;
   int           num_scoarsenCFL;
   int           use_rand;
   double       *scoarsenCFL;

   MPI_Comm    comm, comm_x, comm_t;
   int         myid, num_procs;
   double      mystarttime, myendtime, mytime, maxtime;

   /* We consider a 2D problem. */
   int ndim = 2;

   /* diffusion coefficient */
   double K; 

   int nx, ny, nlx, nly;
   double tstart, tstop;
   int nt;
   double dx, dy, dt, cfl;
   int ilower[2], iupper[2];

   /* We have one part and one variable.
   int nparts = 1; */
   int nvars = 1;
   int var;

   /* We want to use a struct solver. */
   int object_type = HYPRE_STRUCT;


   int px, py, pt;
   int pi, pj;

   int n_pre, n_post;
   int rap, relax, skip, max_iter_x[2];
   double tol_x[2], tol_x_coarse;

   int nA_max;
   int *max_num_iterations_global = NULL;
   double *scoarsen_table_global = NULL;

   int output_files, explicit, output_vis, print_level, access_level;
   int run_wrapper_tests;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);

   /* Default parameters. */
   comm                = MPI_COMM_WORLD;
   max_levels          = 2;  /* Must be two, in order for coarsen and refine wrapper tests to run */
   min_coarse          = 3;  
   nrelax              = 1;
   nrelax0             = -1;
   tol                 = 1.0e-09;
   tnorm               = 2;
   cfactor             = 2;
   cfactor0            = -1;
   max_iter            = 100;
   fmg                 = 0;
   nfmg_Vcyc           = 1;
   use_rand            = 1;
   scoarsen            = 0;
   num_scoarsenCFL     = 1;
   scoarsenCFL         = (double*) malloc( 1*sizeof(double) );
   scoarsenCFL[0]      = 0.5;
   K                   = 1.0;
   nx                  = 17;
   ny                  = 17;
   nlx                 = 17;
   nly                 = 17;
   tstart              = 0.0;
   nt                  = 32;
   cfl                 = 0.30;
   px                  = 1;
   py                  = 1;
   pt                  = 1;  
   n_pre               = 1;
   n_post              = 1;
   rap                 = 1;
   relax               = 3;
   skip                = 1;
   max_iter_x[0]       = 50;
   max_iter_x[1]       = 50;
   tol_x[0]            = 1.0e-09;
   tol_x[1]            = 1.0e-09;
   tol_x_coarse        = 1.0e-09;
   explicit            = 0;
   output_files        = 0;
   output_vis          = 0;
   print_level         = 1;
   access_level        = 1;
   run_wrapper_tests   = 0;

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
      else if( strcmp(argv[arg_index], "-nx") == 0 ){
         arg_index++;
         nx = atoi(argv[arg_index++]);
         ny = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-nt") == 0 ){
          arg_index++;
          nt = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-cfl") == 0 ){
          arg_index++;
          cfl = atof(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-ml") == 0 ){
          arg_index++;
          max_levels = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-mc") == 0 ){
          arg_index++;
          min_coarse = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-nu") == 0 ){
          arg_index++;
          nrelax = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-nu0") == 0 ){
          arg_index++;
          nrelax0 = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-tol") == 0 ){
          arg_index++;
          tol = atof(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-tnorm") == 0 ){
          arg_index++;
          tnorm = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-cf") == 0 ){
          arg_index++;
          cfactor = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-cf0") == 0 ){
          arg_index++;
          cfactor0 = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-mi") == 0 ){
          arg_index++;
          max_iter = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-fmg") == 0 ){
         arg_index++;
         fmg = 1;
         nfmg_Vcyc = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-use_rand") == 0 ){
         arg_index++;
         use_rand = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-scoarsen") == 0 ){
         arg_index++;
         scoarsen = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-scoarsenCFL") == 0 ){
         arg_index++;
         num_scoarsenCFL = atoi(argv[arg_index++]);
         free(scoarsenCFL);
         scoarsenCFL = (double*) malloc( num_scoarsenCFL*sizeof(double) );
         for (i = 0; i < num_scoarsenCFL; i++)
            scoarsenCFL[i] = atof(argv[arg_index++]);
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
      else if( strcmp(argv[arg_index], "-iter") == 0 ){
         arg_index++;
         max_iter_x[0] = atoi(argv[arg_index++]);
         max_iter_x[1] = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-tolx") == 0 ){
          arg_index++;
          tol_x[0] = atof(argv[arg_index++]);
          tol_x[1] = atof(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-tolxc") == 0 ){
          arg_index++;
          tol_x_coarse = atof(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-expl") == 0 ){
         arg_index++;
         explicit = 1;
      }
      else if( strcmp(argv[arg_index], "-print_level") == 0 ){
         arg_index++;
         print_level = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-access_level") == 0 ){
         arg_index++;
         access_level = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-run_wrapper_tests") == 0 ){
         arg_index++;
         run_wrapper_tests = 1;
      }
      else if( strcmp(argv[arg_index], "-output_files") == 0 ){
         arg_index++;
         output_files = 1;
      }
      else if( strcmp(argv[arg_index], "-output_vis") == 0 ){
         arg_index++;
         output_vis = 1;
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
      printf(" General XBraid configuration parameters\n");
      printf(" ---------------------------------------\n");
      printf("  -run_wrapper_tests               : Only run the Braid wrapper tests\n");
      printf("                                     (do not combine with temporal parallelism)\n");
      printf("  -pgrid  <px py pt>               : processors in each dimension (default: 1 1 1)\n");
      printf("  -nx  <nlx nly>                   : 2D spatial problem size of form 2^k+1, 2^k+1 (default: 17 17)\n");
      printf("  -nt  <n>                         : number of time steps (default: 32)\n"); 
      printf("  -cfl <cfl>                       : CFL number to run, note that 2*CFL = dt/(dx^2) (default: 0.30)\n"); 
      printf("                                     CFL < 0.5 for explicit forward Euler\n");
      printf("  -ml  <max_levels>                : set max number of time levels (default: 1)\n");
      printf("  -mc  <min_coarse>                : set min possible coarse level size (default: 3)\n");
      printf("  -nu  <nrelax>                    : set num F-C relaxations (default: 1)\n");
      printf("  -nu0 <nrelax>                    : set num F-C relaxations on level 0\n");
      printf("  -tol <tol>                       : set stopping tolerance (default: 1e-09)\n");
      printf("  -tnorm <tnorm>                   : set temporal norm \n");
      printf("                                     1 - One-norm \n");
      printf("                                     2 - Two-norm (default) \n");
      printf("                                     3 - Infinity-norm \n");
      printf("  -cf  <cfactor>                   : set coarsening factor (default: 2)\n");   
      printf("  -cf0  <cfactor>                  : set coarsening factor for level 0 \n");
      printf("  -mi  <max_iter>                  : set max iterations (default: 100)\n");
      printf("  -iter <max_iter max_iter_cheap>  : maximum number of PFMG iterations (default: 50 50)\n"); 
      printf("  -tolx <loose_tol tight_tol>      : loose and tight stopping tolerance for PFMG (default: 1e-09 1e-09)\n"); 
      printf("  -tolxc <tol_x>                   : stopping tolerance for PFMG on coarse grids (default: 1e-09)\n");
      printf("  -fmg <nfmg_Vcyc>                 : use FMG cycling, nfmg_Vcyc V-cycles at each fmg level\n");
      printf("  -use_rand <bool>                 : if nonzero, then use a uniformly random value to initialize each\n");
      printf("                                     time step.  if zero, then use the constant value 1.42 to initialize.\n");
      printf("                                     \n");
      printf("                                     \n");
      printf(" Spatial Coarsening related parameters  \n");
      printf(" -------------------------------------  \n");
      printf("                                     \n");
      printf(" In general, spatial coarsening should be used with explicit time stepping \n");
      printf(" (the \"-expl\" option) in order to maintain stability on the coarse grid.  But it \n");
      printf(" can also be used with implicit time stepping for efficiency reasons.  To use spatial \n");
      printf(" coarsening, the \"-scoarsen <s>\" option must be used.  To tailor when and when not to \n");
      printf(" use spatial coarsening, use the \"-scoarsen_CFL <n> <cfls> \" option. \n");
      printf("                                     \n");
      printf(" As a general rule, the code will switch to implicit time stepping on a coarse grid if the \n");
      printf(" CFL number exceeds the threshhold for explicit time stepping (0.5).\n");
      printf("                                     \n");
      printf(" All spatial coarsening and refinements are done on nice structured grids of size 2^k+1 x 2^k + 1.\n");
      printf(" This is just a demonstration driver, nothing fancy! \n");
      printf("                                     \n");
      printf("  -expl <e>                        : use explicit scheme\n");
      printf("  -scoarsen <s>                    : use spatial coarsening when needed to satisfy the CFL constraints \n");
      printf("                                     specified by the scoarsenCFL option \n");
      printf("                                     0 - No spatial coarsening (default) \n");
      printf("                                     1 - Use injection for spatial restriction \n");
      printf("                                     2 - Use transpose of bilinear interpolation for spatial restriction\n");
      printf("  -scoarsenCFL <n> <cfls>          : array of CFL numbers such that if the actual CFL at level k is\n");
      printf("                                     greater that scoarsenCFL[k] then do spatial coarsening until the\n");
      printf("                                     CFL < scoarsenCFL[k].  <n> specifies the length of the array and <cfls>\n");
      printf("                                     specifies the entries of the array.\n");
      printf("                                     \n");
      printf("                                     For example, -scoarsenCFl 1 0.4 would tell Braid to spatially coarsen \n");
      printf("                                     whenever the CFL number exceeds 0.4.\n");
      printf("                                     As another example, -scoarsenCFl 2 10.0 1.0 would tell Braid to \n");
      printf("                                     spatially coarsen on the first coarse grid only if the CFL number\n");
      printf("                                     is greater than 10.0.  Then on all subsequent coarse grids, spatial\n");
      printf("                                     coarsening would occur only if the CFL number is greater than 1.0.\n");
      printf("                                     \n");
      printf("                                     Notes: valid CFLs for explicit time stepping are in [0, 0.5],\n");
      printf("                                     \n");
      printf("                                     Default is n=1 and cfls=0.5.\n");
      printf("                                     and this option must be used with scoarsen > 0.\n");
      printf("                                     \n");
      printf("                                     \n");
      printf(" PFMG related parameters (only used if implicit time stepping is chosen) \n");
      printf(" ----------------------------------------------------------------------- \n");
      printf("  -v <n_pre> <n_post>              : number of pre and post relaxations in PFMG\n");
      printf("  -rap <r>                         : coarse grid operator type in PFMG\n");
      printf("                                     0 - Galerkin (default)\n");
      printf("                                     1 - non-Galerkin ParFlow operators\n");
      printf("                                     2 - Galerkin, general operators\n");
      printf("  -relax <r>                       : relaxation type in PFMG\n");
      printf("                                     0 - Jacobi\n");
      printf("                                     1 - Weighted Jacobi\n");
      printf("                                     2 - R/B Gauss-Seidel (default)\n");
      printf("                                     3 - R/B Gauss-Seidel (nonsymmetric)\n");
      printf("  -skip <s>                       : skip levels in PFMG (0 or 1)\n");      
      printf("                                     \n");
      printf("                                     \n");
      printf(" Output related parameters\n");
      printf(" -------------------------\n");
      printf("  -print_level <l>                : sets the print_level (default: 1) \n");
      printf("                                    0 - no output to standard out \n");
      printf("                                    1 - Basic convergence information and hierarchy statistics\n");
      printf("                                    2 - Debug level output \n");
      printf("  -access_level <l>               : sets the access_level (default: 1) \n");
      printf("                                    0 - never call access \n");
      printf("                                    1 - call access only after completion \n");
      printf("                                    2 - call access every iteration and level\n");
      printf("  -output_files                   : save the solution/error/error norms to files\n");
      printf("                                    frequency of file accesss is set by access_level\n");
      printf("  -output_vis                     : save the error for GLVis visualization\n");
      printf("                                    frequency of file accesss is set by access_level\n");
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

   /* Check the spatial grid size, must be power of 2 + 1*/
   if( ( (pow(2.0, floor(log2((double) nx-1))) + 1 - nx) != 0) || 
       ( (pow(2.0, floor(log2((double) ny-1))) + 1 - ny) != 0)    )
   {
      if( myid == 0 )
         printf("Error: nx and ny must be powers of 2 plus 1, e.g., 17, 33, 65 and so on\n");
      MPI_Finalize();
      return (0);
   }
   
   /* Check that domain is square */
   if( nx != ny )
   {
      if( myid == 0 )
         printf("Error: nx must equal ny.  This is a simple test driver :)  \n");
      MPI_Finalize();
      return (0);
   }

   /* Create communicators for the time and space dimensions */
   braid_SplitCommworld(&comm, px*py, &comm_x, &comm_t);

#if DEBUG
   MPI_Comm_size( comm_t, &num_procs );
   printf( "number of processors in time:  %d\n", num_procs );
   MPI_Comm_size( comm_x, &num_procs );
   printf( "number of processors in space: %d\n", num_procs );
#endif
  
   /* Determine position in the processor grid. */
   MPI_Comm_rank( comm_x, &myid );
   pi = myid % px;
   pj = ( (myid - pi)/px ) % py;

   /* Define the nodes owned by the current processor (each processor's
    * piece of the global grid) */
   GetDistribution_x( nx, px, pi, &ilower[0], &iupper[0] );
   GetDistribution_x( ny, py, pj, &ilower[1], &iupper[1] );

   /* Determine local problem size. */
   nlx = iupper[0] - ilower[0] + 1;
   nly = iupper[1] - ilower[1] + 1;

#if DEBUG
   printf( "%d = (%d, %d): nlx = %d, nly = %d\n", 
           myid, pi, pj, nlx, nly ); 
#endif

   /* Compute grid spacing. */
   dx = PI / (nx - 1);
   dy = PI / (ny - 1);

   /* Set time-step size. */
   dt = K*cfl*( (dx*dx)*(dy*dy) / ((dx*dx)+(dy*dy)) );
   /* Determine tstop. */
   tstop =  tstart + nt*dt;

   /* -----------------------------------------------------------------
    * Set up App structure.
    * ----------------------------------------------------------------- */
   app = (my_App *) malloc(sizeof(my_App));
   (app->comm_t)          = comm_t;
   (app->comm_x)          = comm_x;
   (app->dim_x)           = ndim;
   (app->K)               = K;
   (app->nlx)             = nlx;
   (app->nly)             = nly;
   (app->tstart)          = tstart;
   (app->tstop)           = tstop;
   (app->nt)              = nt;
   (app->dx)              = dx;
   (app->dy)              = dy;
   (app->dt)              = dt;
   (app->max_levels)      = max_levels;
   (app->px)              = px;
   (app->py)              = py;
   (app->pi)              = pi;
   (app->pj)              = pj;
   (app->ilower[0])       = ilower[0];
   (app->ilower[1])       = ilower[1];
   (app->iupper[0])       = iupper[0];
   (app->iupper[1])       = iupper[1];
   (app->object_type)     = object_type;
   (app->n_pre)           = n_pre;
   (app->n_post)          = n_post;
   (app->rap)             = rap;
   (app->relax)           = relax;
   (app->skip)            = skip;
   (app->explicit)        = explicit;
   (app->last_tsize)      = -1.0;
   (app->output_files)    = output_files;
   (app->output_vis)      = output_vis;

   /* Set the maximum number of PFMG iterations for expensive (index 0)
    * and cheap (index 1) solves. */
   (app->max_iter_x)      = (int*) malloc( 2*sizeof(int) );
   (app->max_iter_x[0])   = max_iter_x[0];
   (app->max_iter_x[1])   = max_iter_x[1];

   /* Set the loose (index 0) and tight (index 1) stopping tolerance
    * for PFMG. */
   (app->tol_x)           = (double*) malloc( 2*sizeof(double) );
   (app->tol_x[0])        = tol_x[0];
   (app->tol_x[1])        = tol_x[1];

   /* Set the variable types. */
   (app->vartypes)   = (HYPRE_SStructVariable*) malloc( nvars* 
                                         sizeof(HYPRE_SStructVariable) );
   for( var = 0; var < nvars; var++ )
      app->vartypes[var] = HYPRE_SSTRUCT_VARIABLE_CELL;

   /* Set up a 2D grid. */
   setUp2Dgrid( app->comm_x, &(app->grid_x), app->dim_x,
                app->ilower, app->iupper, app->vartypes, 1 );

   /* Define the discretization stencil. */
   set5ptStencil( &(app->stencil), app->dim_x );

   /* Set up the graph - this determines the non-zero structure of the
    * matrix. */ 
   setUpGraph( app->comm_x, &(app->graph), app->grid_x, app->object_type, 
               app->stencil );
   
   /* Initialize the storage structure for recording spatial coarsening information */ 
   app->scoarsen_table = (double*) malloc( 5*(app->max_levels)*sizeof(double) );
   for( i = 0; i < 5*(app->max_levels); i++)
   {
      app->scoarsen_table[i] = -1.0;
   }
   
   /* Store whether we want a random or constant initial guess */
   app->use_rand = use_rand;
   
   /* Store CFl that controls spatial coarsening */
   app->num_scoarsenCFL = num_scoarsenCFL;
   app->scoarsenCFL = scoarsenCFL;

   /* Allocate memory for array of discretization matrices. */
   app->A = (HYPRE_SStructMatrix*) malloc( (app->max_levels)*
                                           sizeof(HYPRE_SStructMatrix));
   /* Create empty matrix lookup tables for dt, dx and dy. */
   app->dt_A = (double*) malloc( (app->max_levels)*sizeof(double) );
   app->dy_A = (double*) malloc( (app->max_levels)*sizeof(double) );
   app->dx_A = (double*) malloc( (app->max_levels)*sizeof(double) );
   for( i = 0; i < app->max_levels; i++ )
   {
      app->dt_A[i] = -1.0;
      app->dx_A[i] = -1.0;
      app->dy_A[i] = -1.0;
   }
   app->nA = 0;

   /* Allocate memory for array of solvers. */
   app->solver = (HYPRE_StructSolver*) malloc( (app->max_levels)*
                                               sizeof(HYPRE_StructSolver));

   /* Allocate memory for array of iteration counts. */
   app->max_num_iterations = (int*) calloc( (app->max_levels),  sizeof(int) );
   for( i = 0; i < app->max_levels; i++ )
      app->max_num_iterations[i] = 0;

   /* Setup the lookup table that records how grids are coarsened (refined)
    * spatially.  
    * Note: that the first entry for the table is a dummy entry 
    *       for cloning vectors.  
    * Note: We assume that no more than 3*(app->max_levels) instances
    *       of spatial coarsening will ever occur.*/
   (app->spatial_disc_table) = (spatial_discretization*) malloc( 3*max_levels*sizeof(spatial_discretization) );
   for( i = 1; i < 3*(app->max_levels); i++ )
   {
      app->spatial_disc_table[i].fdt = -1.0;
      app->spatial_disc_table[i].cdt = -1.0;
   }
   (app->spatial_disc_table[0]).grid_x_matrix = (app->grid_x);
   (app->spatial_disc_table[0]).grid_x_vector = (app->grid_x);
   (app->spatial_disc_table[0]).graph_matrix = (app->graph);
   (app->spatial_disc_table[0]).graph_vector = (app->graph);
   (app->spatial_disc_table[0]).fdt = (app->dt);
   (app->spatial_disc_table[0]).cdt = (app->dt);
   (app->spatial_disc_table[0]).ncoarsen = 0;
   (app->spatial_disc_table[0]).dx = (app->dx);
   (app->spatial_disc_table[0]).dy = (app->dy);
   (app->spatial_disc_table[0]).nx = nx;
   (app->spatial_disc_table[0]).ny = ny;
   (app->spatial_disc_table[0]).nlx = (app->nlx);
   (app->spatial_disc_table[0]).nly = (app->nly);
   (app->spatial_disc_table[0]).ilower[0] = (app->ilower[0]);
   (app->spatial_disc_table[0]).ilower[1] = (app->ilower[1]);
   (app->spatial_disc_table[0]).iupper[0] = (app->iupper[0]);
   (app->spatial_disc_table[0]).iupper[1] = (app->iupper[1]);
   (app->spatial_disc_table[0]).fspatial_disc_idx = 0;
   
   
   if( run_wrapper_tests)
   {
      /* Run only the Braid wrapper tests */
      MPI_Comm_rank( comm_x, &myid );

      /* Test init(), access(), free() */
      braid_TestInitAccess( app, comm_x, stdout, 0.0, my_Init, my_Access, my_Free);
      braid_TestInitAccess( app, comm_x, stdout, dt, my_Init, my_Access, my_Free);

      /* Test clone() */
      braid_TestClone( app, comm_x, stdout, 0.0, my_Init, my_Access, my_Free, my_Clone);
      braid_TestClone( app, comm_x, stdout, dt, my_Init, my_Access, my_Free, my_Clone);

      /* Test sum() */
      braid_TestSum( app, comm_x, stdout, 0.0, my_Init, my_Access, my_Free, my_Clone, my_Sum);
      braid_TestSum( app, comm_x, stdout, dt, my_Init, my_Access, my_Free, my_Clone, my_Sum);

      /* Test spatialnorm() */
      correct = braid_TestSpatialNorm( app, comm_x, stdout, 0.0, my_Init, my_Free, my_Clone, my_Sum, my_SpatialNorm);
      correct = braid_TestSpatialNorm( app, comm_x, stdout, dt, my_Init, my_Free, my_Clone, my_Sum, my_SpatialNorm);

      /* Test bufsize(), bufpack(), bufunpack() */
      correct = braid_TestBuf( app, comm_x, stdout, 0.0, my_Init, my_Free, my_Sum, my_SpatialNorm, my_BufSize, my_BufPack, my_BufUnpack);
      correct = braid_TestBuf( app, comm_x, stdout, dt, my_Init, my_Free, my_Sum, my_SpatialNorm, my_BufSize, my_BufPack, my_BufUnpack);
       
      /* Test coarsen and refine */
      correct = braid_TestCoarsenRefine(app, comm_x, stdout, 0.0, dt, 2*dt, my_Init,
                             my_Access, my_Free, my_Clone, my_Sum, my_SpatialNorm, my_CoarsenInjection, 
                             my_Refine);

      correct = braid_TestCoarsenRefine(app, comm_x, stdout, 0.0, dt, 2*dt, my_Init,
                             my_Access, my_Free, my_Clone, my_Sum, my_SpatialNorm, my_CoarsenBilinear, 
                             my_Refine);

      if(correct == 0)
      {
        printf("Failed: at least one of the tests failed\n");
      }
      else
      {
        printf("Passed: all tests passed\n");
      }

   }
   else
   {
      /* Run a Braid simulation */
      /* Start timer. */
      mystarttime = MPI_Wtime();

      braid_Init(comm, comm_t, tstart, tstop, nt, app, my_Phi, my_Init,
            my_Clone, my_Free, my_Sum, my_SpatialNorm, my_Access, my_BufSize,
            my_BufPack, my_BufUnpack, &core);

      braid_SetLoosexTol( core, 0, tol_x[0] );
      braid_SetLoosexTol( core, 1, tol_x_coarse );

      braid_SetTightxTol( core, 0, tol_x[1] );

      braid_SetMaxLevels( core, max_levels );
      braid_SetMinCoarse( core, min_coarse );

      braid_SetPrintLevel( core, print_level);
      braid_SetAccessLevel( core, access_level);

      braid_SetNRelax(core, -1, nrelax);
      if (nrelax0 > -1)
      {
         braid_SetNRelax(core,  0, nrelax0);
      }

      /*braid_SetRelTol(core, tol);*/
      /*braid_SetAbsTol(core, tol*sqrt(px*nlx*py*nly*(nt+1)) );*/
      braid_SetAbsTol(core, tol/sqrt(dx*dy*dt));
      braid_SetTemporalNorm(core, tnorm);

      braid_SetCFactor(core, -1, cfactor);
      if( cfactor0 > 0 )
      {
         braid_SetCFactor(core,  0, cfactor0);
      }
      
      braid_SetMaxIter(core, max_iter);
      if (fmg)
      {
         braid_SetFMG(core);
         braid_SetNFMGVcyc(core, nfmg_Vcyc);
      }
      
      if (scoarsen)
      {
         app->scoarsen=1;
         if (scoarsen == 1)
         {
            braid_SetSpatialCoarsen(core, my_CoarsenInjection);
         }
         else if (scoarsen == 2)
         {
            braid_SetSpatialCoarsen(core, my_CoarsenBilinear);
         }
         else
         {
            printf("Invalid scoarsen choice.  Ignoring this parameter\n");
         }
         braid_SetSpatialRefine(core, my_Refine);
      }

      /* Print some additional statistics */
      MPI_Comm_rank( comm, &myid );

      if( myid == 0 )
      {
         printf("\n-----------------------------------------------------------------\n"); 
         printf("-----------------------------------------------------------------\n\n"); 
         
         printf("  Begin simulation \n\n");
      }
      
      braid_Drive(core);

      /* Stop timer. */
      myendtime = MPI_Wtime();
      mytime    = myendtime - mystarttime;

      /* Compute maximum time */
      MPI_Reduce( &mytime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, comm );

      /* Determine some solver information about topics like spatial coarsening
       * and the maximum number of iterations in spatial solves on each time level
       * (if implicit used) */
      MPI_Allreduce( &(app->nA), &nA_max, 1, MPI_INT, MPI_MAX, comm ); 
      max_num_iterations_global = (int*) malloc( nA_max*sizeof(int) );
      scoarsen_table_global = (double*) malloc( 5*nA_max*sizeof(double) ); 
      for( i = 0; i < nA_max; i++ ){
         /* Grab max num interations information */
         MPI_Allreduce( &(app->max_num_iterations[i]), 
                        &max_num_iterations_global[i], 1, MPI_INT, MPI_MAX, comm );

         /* Grab spatial coarsening information */
         MPI_Allreduce( &(app->scoarsen_table)[i*5],
               &scoarsen_table_global[i*5], 1, MPI_DOUBLE, MPI_MAX, comm ); 
         MPI_Allreduce( &(app->scoarsen_table)[i*5+1],
               &scoarsen_table_global[i*5+1], 1, MPI_DOUBLE, MPI_MAX, comm ); 
         MPI_Allreduce( &(app->scoarsen_table)[i*5+2],
               &scoarsen_table_global[i*5+2], 1, MPI_DOUBLE, MPI_MAX, comm ); 
         MPI_Allreduce( &(app->scoarsen_table)[i*5+3],
               &scoarsen_table_global[i*5+3], 1, MPI_DOUBLE, MPI_MAX, comm ); 
         MPI_Allreduce( &(app->scoarsen_table)[i*5+4],
               &scoarsen_table_global[i*5+4], 1, MPI_DOUBLE, MPI_MAX, comm ); 
      
      }

      if( myid == 0 )
      {
         printf( "  runtime: %.5lfs\n", maxtime );

         printf("\n\n-----------------------------------------------------------------\n"); 
         printf("-----------------------------------------------------------------\n\n"); 
         printf(" Implicit time stepping solve parameters\n\n");
         printf("   Fine-level loose stopping tol  :  %1.2e    (while ||r|| is large)\n", app->tol_x[0]);
         printf("   Fine-level tight stopping tol  :  %1.2e    (while ||r|| is small)\n", app->tol_x[0]);
         printf("   Coarse-level stopping tol      :  %1.2e    (for all ||r||) \n", tol_x_coarse);
         printf(" \n"); 
         printf("   Fine-level max iter            :  %d\n", app->max_iter_x[0]);
         printf("   Coarse-level max iter          :  %d\n", app->max_iter_x[1]);

         printf("\n-----------------------------------------------------------------\n"); 
         printf("-----------------------------------------------------------------\n\n"); 
         
         printf( " Per level diagnostic information \n\n");
         printf(" Scheme is either explicit or implicit.  If implicit then the table\n");
         printf(" entry is 'impl, %%d' where the integer represents the maximum number\n");
         printf(" of AMG iterations used by the implicit solver on that level.  'expl'\n");
         printf(" stands for explicit. \n\n"); 
         printf(" Fine level spatial problem size  : %d x %d\n\n", nx, ny );
         printf("level   scheme         dx          dy          dt          cfl\n"); 
         printf("-----------------------------------------------------------------\n"); 
         for( i = 0; i < nA_max; i++)
         {
            if( scoarsen_table_global[i*5] == 1)
            {
               printf(" %2d   |  expl       %1.2e    %1.2e    %1.2e    %1.2e\n", i,
                  scoarsen_table_global[i*5+1], scoarsen_table_global[i*5+2],
                  scoarsen_table_global[i*5+3], scoarsen_table_global[i*5+4]);
            }
            else
            {
               printf(" %2d   |  impl, %2d   %1.2e    %1.2e    %1.2e    %1.2e\n", 
                  i, max_num_iterations_global[i],
                  scoarsen_table_global[i*5+1], scoarsen_table_global[i*5+2],
                  scoarsen_table_global[i*5+3], scoarsen_table_global[i*5+4]);
            }
         }

         printf( "\n" );
      }

      braid_Destroy(core);
   }
   
   /* Free memory */
   if( max_num_iterations_global != NULL)
   {
      free( max_num_iterations_global );
   }
   HYPRE_SStructGridDestroy( app->grid_x );
   HYPRE_SStructStencilDestroy( app->stencil );
   HYPRE_SStructGraphDestroy( app->graph );
   free( app->vartypes );
   free( app->dt_A );
   free( app->dx_A );
   free( app->dy_A );
   free( app->scoarsenCFL);
   for( i = 0; i < app->nA; i++ )
   {
      HYPRE_SStructMatrixDestroy( app->A[i] );
      if( scoarsen_table_global != NULL)
      {
         if( scoarsen_table_global[i*5] == 0 )
         {
            /* If an implicit solver was used */
            HYPRE_StructPFMGDestroy( app->solver[i] );
         }
      }
   }
   free( scoarsen_table_global );
   free( app->A );
   free( app->solver );
   free( app->max_num_iterations );
   free( app->scoarsen_table );
   free( app->max_iter_x );
   
   /* Destroy all coarse grid_x and graphs
    * The app->grid_x and app->graph are destroyed above, and correspond to 
    * the i=0 matrix graph and grid entries below.  So, we still need to destory the
    * i=0 vector graph and grid entries, all middle level vector AND matrix graph 
    * and grid entries, and only the coarsest level matrix graph and grid entries.
    * At the coarsest level, the vector and matrix entries are the same, since no
    * coarsening is needed there.
    *
    * Note: that the vector grid and graph use the same grid and graph as the
    *       _matrix versions until coarsen is called and a new grid and graph are
    *       generated with an expanded ghost layer.  So, only deallocate if not
    *       at coarsest level
    * Note: We assume that no more than 3*(app->max_levels) instances
    *       of spatial coarsening will ever occur.
    */
   for( i = 1; i < 3*(app->max_levels); i++ )
   {
      /* Note that if cdt and fdt == -1.0, then this level has never been used */
      if( ( (app->spatial_disc_table[i]).cdt != -1.0 ) && ( (app->spatial_disc_table[i]).fdt != -1.0 ) )
      {
         HYPRE_SStructGridDestroy( (app->spatial_disc_table[i]).grid_x_matrix );
         HYPRE_SStructGraphDestroy( (app->spatial_disc_table[i]).graph_matrix );
         fspatial_disc_idx = (app->spatial_disc_table[i]).fspatial_disc_idx;

         /* Note that coarse spatial grid creation also touches its relative fine grid with a
          * _vector creation of grid and graph with large ghost layers.  Destroy that */
         if(fspatial_disc_idx != -1)
         {
            HYPRE_SStructGridDestroy( (app->spatial_disc_table[fspatial_disc_idx]).grid_x_vector );
            HYPRE_SStructGraphDestroy( (app->spatial_disc_table[fspatial_disc_idx]).graph_vector );
         }
      }
   }
   /* Note that this always allocated */
   free( app->spatial_disc_table );
   
   free( app->tol_x );
   free( app );
   MPI_Comm_free( &comm_x );
   MPI_Comm_free( &comm_t );

   /* Finalize MPI */
   MPI_Finalize();

   return 0;
}





