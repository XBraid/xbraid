/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder schroder2@llnl.gov, Rob Falgout falgout2@llnl.gov,
 * Tzanio Kolev kolev1@llnl.gov, Ulrike Yang yang11@llnl.gov, 
 * Veselin Dobrev dobrev1@llnl.gov, et al. 
 * LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. Email schroder2@llnl.gov on how to download. 
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

   Sample run:   mpirun -np 8 drive-05 -pgrid 1 1 8 -ml 15 -nt 128 -nx 33 33 -mi 100 -expl -scoarsen  

   Notes:        The " -expl and -scoarsen " options should be used together
                 to coarsen spatially and do explicit time stepping.  No spatial
                 parallel distribution is supported.

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
 * This will store the spatial discretization information for moving from 
 * one level with fdt to a coarser level with cdt 
 *
 * fdt                  fine grid time step size         
 * cdt                  coarse grid time step size       
 * dx                   coarse grid dx
 * dy                   coarse grid dy
 * nlx                  coarse grid nlx
 * nly                  coarse grid nly
 * ilower               coarse grid ilower
 * iupper               coarse grid iupper
 * grid_x               SStructGrid for coarse grid
 * graph                SStructGraph for coarse grid
 * fspatial_disc_idx    integer index into app->spatial_disc_table for the 
 *                      fine spatial mesh used to generate this mesh
 * */
typedef struct _spatial_discretization
{
   double fdt, cdt, dx, dy; 
   int ncoarsen, nlx, nly, fspatial_disc_idx;
   int ilower_x[2], iupper_x[2];
   HYPRE_SStructGrid       grid_x;
   HYPRE_SStructGraph      graph;
} spatial_discretization;


/* --------------------------------------------------------------------
 * Structs for my integration routines
 * -------------------------------------------------------------------- */
/* struct my_App contains general information about the problem, its
 * discretizaion, spatial distribution, and solver used at each time point
 *
 *   comm_t          communicator for parallelizing in time
 *   comm_x          communicator for parallelizing in space 
 *   dim_x           spatial dimension
 *   K               diffusion coefficient
 *   nlx             local problem size in x-dimension
 *   nly             local problem size in y-dimension
 *   tstart          initial time
 *   tstop           integration stop time
 *   nt              number of time steps
 *   dx              spatial step size in x-direction
 *   dy              spatial step size in y-direction
 *   dt              time step on finest grid
 *   max_levels      maximum number of time levels
 *   vartypes        (nvars)-dimensional array of variable types on 
 *                   structured parts of the spatial grid
 *   grid_x          spatial grid
 *   stencil         discretization stencil object
 *   graph           graph object that determine the non-zero structure
 *                   of the discretization matrices
 *   A               array of discretization matrices (one per time level)
 *   dt_A            array of time steps for which discretization matrix
 *                   has been created
 *   nA              number of discretization matrices that have been
 *                   created
 *   sym             symmetric storage (1) or not (0)
 *   px              number of processors in x-dimension
 *   py              number of processors in y-dimension
 *   pi              x-coordinate of position in processor grid
 *   pj              y-coordinate of position in processor grid
 *   ilower_x        (dim_x)-dimensional array with integer indices of 
 *                   local space interval lower bounds
 *   iupper_x        (dim_x)-dimensional array with integer indices of
 *                   local space interval upper bounds
 *   object_type     object type of vector to access different hypre solvers 
 *   solver          array of solvers used at each time step on different
 *                   time levels
 *   n_pre           number of pre-relaxation sweeps in spatial solve
 *   n_post          number of post-relaxation sweeps in spatial solve
 *   rap             coarse-grid operator type for spatial solver
 *                      rap = 0: Galerkin (default)
 *                      rap = 1: non-Galerkin ParFlow operators
 *                      rap = 2: Galerkin, general operators 
 *   relax           relaxation type for spatial solver
 *                      relax = 0: Jacobi
 *                      relax = 1: Weighted Jacobi 
 *                      relax = 2: R/B Gauss-Seidel (default)
 *                      relax = 3: R/B Gauss-Seidel (non-symmetric)
 *   skip            skip levels in spatial PFMG solver (0 or 1)
 *   max_iter_x      expensive and cheap maximum number of spatial MG iterations
 *   tol_x           loose and tight stopping tolerance for spatial MG
 *   explicit        use explicit discretization (1) or not (0)
 *   scoarsen        use spatial refinement and coarsening
 *   spatial_disc_table    Lookup table recording dx, dy when coarsening spatially
 *   last_tsize      output related flag that lets my_Phi know when the 
 *                   time step size has changed
 *   scheme          int array of integration scheme used: explicit or
 *                   implicit 
 *   output_files    save the solution/error/error norm to files
 *   output_vis      save the error for GLVis visualization
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
   int                    *max_iter_x;
   double                 *tol_x;
   double                 *scoarsen_table;
   int                     explicit;
   int                     scoarsen;
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
             int                  nlx, 
             int                  nly,
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
   if (pj == 0){
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];
           
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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1] + nly-1;
           
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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];
           
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
      bc_ilower[0] = ilower[0] + nlx-1;
      bc_ilower[1] = ilower[1];
           
      bc_iupper[0] = bc_ilower[0];
      bc_iupper[1] = bc_ilower[1] + nly-1;
          
      /* Put the boundary conditions in b */
      for( j = 0; j < nly; j++ )
          bvalues[j] = B0(  bc_ilower[0]*dx,
                           (bc_ilower[1]+j)*dy );
           
      HYPRE_SStructVectorSetBoxValues(b, part, bc_ilower,
                                      bc_iupper, var, bvalues);
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
                  int                  nlx, 
                  int                  nly,
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
   if (pj == 0){
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];
           
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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1] + nly-1;
           
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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];
           
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
      bc_ilower[0] = ilower[0] + nlx-1;
      bc_ilower[1] = ilower[1];
           
      bc_iupper[0] = bc_ilower[0];
      bc_iupper[1] = bc_ilower[1] + nly-1;
          
      /* Put the boundary conditions in b */
      for( j = 0; j < nly; j++ )
          bvalues[j] = B0(  bc_ilower[0]*dx,
                           (bc_ilower[1]+j)*dy );
           
      HYPRE_SStructVectorSetBoxValues(b, part, bc_ilower,
                                      bc_iupper, var, bvalues);
   }

   
   /* 
    * Now, account for the boundary conditions contributions to the
    * domain interior from A_ib u_b
    */ 
   
   /* Neighbors of boundary nodes of boundary y = 0 
    * Neighbors are either
    *   i) on same processor as boundary nodes (pj = 0)
    * or
    *   ii) on neighboring processor (pj = 1) 
    * Case ii) only applies if nly = 1 */
    
   /* Neighbors of boundary on same processor */
   if( (nly > 1) && (pj == 0) )
   {
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1] + 1;
        
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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];
        
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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1] + nly-1 - 1;
        
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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1] + nly-1;
        
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
      bc_ilower[0] = ilower[0] + 1;
      bc_ilower[1] = ilower[1];
        
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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];
       
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
      bc_ilower[0] = ilower[0] + nlx-1 - 1;
      bc_ilower[1] = ilower[1];
        
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
      bc_ilower[0] = ilower[0] + nlx-1;
      bc_ilower[1] = ilower[1];
        
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
             HYPRE_SStructVariable *vartypes )
{
   /* We have one part and one variable. */
   int nvars = 1;
   int nparts = 1;
   int part = 0;

   HYPRE_SStructGrid grid;

   /* Create an empty 2D grid object. */
   HYPRE_SStructGridCreate( comm, ndim, nparts, &grid );

   /* Add a new box to the grid. */
   HYPRE_SStructGridSetExtents( grid, part, ilower, iupper );

   /* Set the variable type for each part */
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

#if DEBUG
   printf( "My box: [%d %d] x [%d %d]\n", 
           ilower[0], iupper[0], ilower[1], iupper[1] );
#endif

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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];
           
      bc_iupper[0] = bc_ilower[0] + nlx-1;
      bc_iupper[1] = bc_ilower[1];
           
      /* Modify the matrix. */
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, nentries,
                                      stencil_indices, values);
   }
       
   /* Processors at y = PI */
   if( pj == py-1 ){
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1] + nly-1;
           
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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];
           
      bc_iupper[0] = bc_ilower[0];
      bc_iupper[1] = bc_ilower[1] + nly-1;
           
      /* Modify the matrix */
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, nentries,
                                      stencil_indices, values);
   }
       
   /* Processors at x = PI */
   if( pi == px-1 ){
      bc_ilower[0] = ilower[0] + nlx-1;
      bc_ilower[1] = ilower[1];
           
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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1] + 1;
           
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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];
           
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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1] + nly-1 - 1;
           
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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1] + nly-1;
           
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
      bc_ilower[0] = ilower[0] + 1;
      bc_ilower[1] = ilower[1];
           
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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];
           
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
      bc_ilower[0] = ilower[0] + nlx-1 - 1;
      bc_ilower[1] = ilower[1];
           
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
      bc_ilower[0] = ilower[0] + nlx-1;
      bc_ilower[1] = ilower[1];
           
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

   double *values;
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
               values[m+idx] = K*(dt/(dx*dx));
            for( idx = 3; idx <= 4; idx++ )
               values[m+idx] = K*(dt/(dy*dy));
                     
            values[m] = 1.0 - 2*K*( (dt/(dx*dx)) + (dt/(dy*dy)) );
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
            values[m+1] = K*(dt/(dx*dx));
            values[m+2] = K*(dt/(dy*dy));  
            values[m]   = 1.0 - 2*K*( (dt/(dx*dx)) + (dt/(dy*dy)) );
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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];
           
      bc_iupper[0] = bc_ilower[0] + nlx-1;
      bc_iupper[1] = bc_ilower[1];
           
      /* Modify the matrix. */
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, nentries,
                                      stencil_indices, values);
   }
       
   /* Processors at y = PI */
   if( pj == py-1 ){
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1] + nly-1;
           
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
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];
           
      bc_iupper[0] = bc_ilower[0];
      bc_iupper[1] = bc_ilower[1] + nly-1;
           
      /* Modify the matrix */
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, nentries,
                                      stencil_indices, values);
   }
       
   /* Processors at x = PI */
   if( pi == px-1 ){
      bc_ilower[0] = ilower[0] + nlx-1;
      bc_ilower[1] = ilower[1];
           
      bc_iupper[0] = bc_ilower[0];
      bc_iupper[1] = bc_ilower[1] + nly-1;
           
      /* Modify the matrix */
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, nentries,
                                      stencil_indices, values);
   }

   free(values);   

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
   int ilower_x[2], iupper_x[2];
   double dx      = (app->spatial_disc_table[u->spatial_disc_idx]).dx;
   double dy      = (app->spatial_disc_table[u->spatial_disc_idx]).dy;
   int    nlx     = (app->spatial_disc_table[u->spatial_disc_idx]).nlx;
   int    nly     = (app->spatial_disc_table[u->spatial_disc_idx]).nly;

   ilower_x[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower_x[0];
   ilower_x[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower_x[1];
   iupper_x[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper_x[0];
   iupper_x[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper_x[1];

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
    * has already been created for time step size tstop-tstart.
    * ----------------------------------------------------------------- */
   A_idx = -1.0;
   for( i = 0; i < app->max_levels; i++ ){
      if( app->dt_A[i] == -1.0 )
      {
         break;
      }
      if( fabs( app->dt_A[i] - (tstop-tstart) )/(tstop-tstart) < 1e-10)
      { 
         A_idx = i;
         break;
      }
   }

   /* Check CFL condition */
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

      /* If we want to use an explicit scheme, check CFL condition 
       * to determine whether we can still use explicit scheme for
       * this time step size or if we have to switch to implicit 
       * scheme. */
      if( app->explicit && cfl ){
         /* Set up the explicit discretization matrix. */
         setUpExplicitMatrix( app->comm_x, &(app->A[A_idx]), 
                              (app->spatial_disc_table[u->spatial_disc_idx]).graph, 
                              app->sym, app->object_type, app->K, dx, 
                              dy, app->dt_A[A_idx], ilower_x, 
                              iupper_x, nlx, nly, app->px,
                              app->py, app->pi, app->pj );
         (app->scoarsen_table)[ (5*i) ]    = 1;
      }
      else{
         /* Set up the implicit discretization matrix. */
         setUpImplicitMatrix( app->comm_x, &(app->A[A_idx]), 
                              (app->spatial_disc_table[u->spatial_disc_idx]).graph, 
                              app->sym, app->object_type, app->K, dx, 
                              dy, app->dt_A[A_idx], ilower_x, 
                              iupper_x, nlx, nly, app->px,
                              app->py, app->pi, app->pj );

      
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
      HYPRE_SStructVectorGetBoxValues( u->x, part, ilower_x,
                                       iupper_x, var, values );
      HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[u->spatial_disc_idx]).grid_x, &b );
      HYPRE_SStructVectorSetObjectType( b, app->object_type );
      HYPRE_SStructVectorInitialize( b );
      HYPRE_SStructVectorSetBoxValues( b, part, ilower_x,
                                       iupper_x, var, values );
      free( values );

      addBoundary( b, app->K, dx, dy, tstop-tstart,
                   ilower_x, nlx, nly, app->px, app->py, 
                   app->pi, app->pj );

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
      HYPRE_SStructVectorGetBoxValues( u->x, part, ilower_x,
                                       iupper_x, var, values );
      HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[u->spatial_disc_idx]).grid_x, &b );
      HYPRE_SStructVectorSetObjectType( b, app->object_type );
      HYPRE_SStructVectorInitialize( b );
      HYPRE_SStructVectorSetBoxValues( b, part, ilower_x,
                                       iupper_x, var, values );
      free( values );
      addBoundaryToRHS( b, app->K, dx, dy, tstop-tstart,
                        ilower_x, nlx, nly, app->px, 
                        app->py, app->pi, app->pj );
      /* add infos from RHS of PDE here */ 

      /* --------------------------------------------------------------
       * Time integration to next time point: Solve the system Ax = b.
       * -------------------------------------------------------------- */
      HYPRE_SStructMatrixGetObject( app->A[A_idx], (void **) &sA );
      HYPRE_SStructVectorGetObject( b, (void **) &sb );
      HYPRE_SStructVectorGetObject( u->x, (void **) &sx );

#if 0
      /* most accurate spatial solves */
      if( accuracy == 0.0 ){   
         /* tolerance-based accuracy */
         HYPRE_StructPFMGSetTol( app->solver[A_idx], app->tol_x[0] );

         /* make sure to change the maximum number of iterations */
         HYPRE_StructPFMGSetMaxIter( app->solver[A_idx], 
                                     app->max_iter_x[0] );
      }

      /* somewhat accurate spatial solves */
      if( accuracy == 0.5 ){    
         /* tolerance-based accuracy */
         HYPRE_StructPFMGSetTol( app->solver[A_idx], app->tol_x[1] );

         /* fixed-iteration based accuracy */  
         HYPRE_StructPFMGSetMaxIter( app->solver[A_idx], 
                                     app->max_iter_x[1] );
      }

      /* inaccurate tolerance-based spatial solves */
      if( (accuracy == 1.0) && (A_idx == 0) ){        
         /* tolerance-based accuracy */
         HYPRE_StructPFMGSetTol( app->solver[A_idx], app->tol_x[2] );

         /* fixed-iteration based accuracy */
         HYPRE_StructPFMGSetMaxIter( app->solver[A_idx], 
                                     app->max_iter_x[2] );
      }

      /* coarse levels */
      if( A_idx > 0 ){     
         /* tolerance-based accuracy */
         HYPRE_StructPFMGSetTol( app->solver[A_idx], 1.0e-02 );

         /* fixed-iteration based accuracy 
         HYPRE_StructPFMGSetMaxIter( app->solver[A_idx], 
                                     3 );*/
      }
      else
         HYPRE_StructPFMGSetTol( app->solver[A_idx], accuracy );
#endif
 
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
   HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[0]).grid_x, &(u->x) );
   
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
         values[m] = ((double)rand())/RAND_MAX;
   }
   HYPRE_SStructVectorSetBoxValues( u->x, part, app->ilower_x,
                                    app->iupper_x, var, values ); 

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
   
   int ilower_x[2], iupper_x[2];
   int    nlx     = (app->spatial_disc_table[u->spatial_disc_idx]).nlx;
   int    nly     = (app->spatial_disc_table[u->spatial_disc_idx]).nly;

   ilower_x[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower_x[0];
   ilower_x[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower_x[1];
   iupper_x[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper_x[0];
   iupper_x[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper_x[1];


   v = (my_Vector *) malloc(sizeof(my_Vector));

   /* Create an empty vector object. */
   HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[u->spatial_disc_idx]).grid_x, &(v->x) );
   
   /* Set the object type (by default HYPRE_SSTRUCT) and initialize. */
   HYPRE_SStructVectorSetObjectType( v->x, app->object_type );
   HYPRE_SStructVectorInitialize( v->x );

   /* Set the values. */
   values = (double *) malloc( (nlx)*(nly)*sizeof(double) );
   HYPRE_SStructVectorGather( u->x );
   HYPRE_SStructVectorGetBoxValues( u->x, part, ilower_x,
                                    iupper_x, var, values );
   HYPRE_SStructVectorSetBoxValues( v->x, part, ilower_x,
                                    iupper_x, var, values );
   free( values );
   HYPRE_SStructVectorAssemble( v->x );

   /* Store the spatial_disc_idx used to generate u */
   v->spatial_disc_idx = u->spatial_disc_idx; 

   *v_ptr = v;

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
 * CFL:   K*( dt/dx^2 + dt/dy^2 ) <  alpha
 *   where alpha < 0.5, and alpha is a constant up to the coder below
 *
 * nlx, nly are the number of points in the x and y dimensions.  If the 
 * required number of coarsenings is too many for this grid, then zero 
 * is returned.
 * -------------------------------------------------------------------- */
int my_ComputeNumCoarsenings(double dt,
                             double dx,
                             double dy,
                             double K,
                             int nlx,
                             int nly)
{
   /* compute required coarsening to satisfy CFL */
   float alpha = 0.49;
   int ncoarsen = ceil( (log2( K*dt*(dx*dx + dy*dy)/(dx*dx*dy*dy) ) - log2(alpha))/2.0 );
   int coarsen_factor =  (int) pow(2.0, ncoarsen);
   if( coarsen_factor == 0)
   {
      return 0;
   }

   /* Check if the existing grid can be coarsened that much */
   int cnlx = (nlx-1) / coarsen_factor  + 1; 
   int cnly = (nly-1) / coarsen_factor  + 1; 
   if( (cnlx < 2) || (cnly < 2) )
   {
      return 0;
   }
   else
   {
      return ncoarsen;
   }
}

/* --------------------------------------------------------------------
 * get the coarse spatial discretization (stored in app->spatial_disc_table)
 * A new spatial discretization is generated, if required by the CFL.
 * Return value is 
 *    spatial_disc_idx
 * -------------------------------------------------------------------- */
void get_coarse_spatial_disc( braid_App app, 
                              double    cdt, 
                              double    fdt, 
                              double    fdx,
                              double    fdy,
                              int       fnlx,
                              int       fnly,
                              int       fspatial_disc_idx,
                              int*      filower_x,
                              int*      fiupper_x,
                              int*      spatial_disc_idx)
{
   int i, ncoarsen, cnlx, cnly;
   int cilower_x[2], ciupper_x[2];
   int max_levels = app->max_levels;
   double coarsen_factor, cdx, cdy;

   /* Search for this cdt, fdt combo and see if a spatial 
    * discretization already exists for it */
   for(i = 0; i < max_levels; i++)
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
    * and store in the next open spot 
    */
   
   /* Determine New sizes */
   ncoarsen = my_ComputeNumCoarsenings(cdt, fdx, fdy, app->K, fnlx, fnly); 
   coarsen_factor = pow(2.0, ncoarsen);
   cdx = coarsen_factor * fdx;
   cdy = coarsen_factor * fdy;
   cnlx = (fnlx-1) / ((int) coarsen_factor)  + 1; 
   cnly = (fnly-1) / ((int) coarsen_factor)  + 1;
   cilower_x[0] = 0;       /* Assumes no spatial parallelism  */
   cilower_x[1] = 0;       /*        ...  */
   ciupper_x[0] = cnlx-1;  /*        ...  */
   ciupper_x[1] = cnly-1;  /*        ...  */


   for(i = 0; i < max_levels; i++)
   {
      if( (app->spatial_disc_table[i]).cdt == -1.0 )
      {
         (app->spatial_disc_table[i]).cdt = cdt;
         (app->spatial_disc_table[i]).fdt = fdt;
         (app->spatial_disc_table[i]).ncoarsen = ncoarsen;
         (app->spatial_disc_table[i]).dx = cdx;
         (app->spatial_disc_table[i]).dy = cdy;
         (app->spatial_disc_table[i]).nlx = cnlx;
         (app->spatial_disc_table[i]).nly = cnly;
         (app->spatial_disc_table[i]).ilower_x[0] = cilower_x[0];
         (app->spatial_disc_table[i]).ilower_x[1] = cilower_x[1];
         (app->spatial_disc_table[i]).iupper_x[0] = ciupper_x[0];
         (app->spatial_disc_table[i]).iupper_x[1] = ciupper_x[1];
         (app->spatial_disc_table[i]).fspatial_disc_idx = fspatial_disc_idx; 
         
         /* We also have to set up a new 2D grid and graph  */
         setUp2Dgrid( app->comm_x, &((app->spatial_disc_table[i]).grid_x), app->dim_x,  
                      cilower_x, ciupper_x, app->vartypes );
         setUpGraph( app->comm_x, &((app->spatial_disc_table[i]).graph), 
                    (app->spatial_disc_table[i]).grid_x, app->object_type, 
                     app->stencil );

         /* This is a return value detailing which rule this is*/
         (*spatial_disc_idx) = i;

         return;
      }
   }
   
   /* No match found, return error value */
   (*spatial_disc_idx) = -1;
}

/* --------------------------------------------------------------------
 * Retrieve a spatial discretization 
 * -------------------------------------------------------------------- */
void retrieve_spatial_discretization( braid_App app, 
                                      double    cdt, 
                                      double    fdt,
                                      int       *spatial_disc_idx)
{
   int max_levels = app->max_levels;
   int i;
   for(i = 0; i < max_levels; i++)
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
 * Create a refined copy of a vector object.
 * Assume a regular grid of size 2^k + 1 in each dimension
 * Assuming no spatial parallelism
 * -------------------------------------------------------------------- */
int
my_Refine(braid_App              app,           
          braid_Vector           cu,
          braid_Vector          *fu_ptr,
          braid_CoarsenRefStatus status)
{
   my_Vector *fu;
   double     *cvalues, *fvalues;
   int        counter, i, j, k;
   int        filower_x[2], fiupper_x[2];
   double     refine_factor;
   int        fnlx, fnly, ncoarsen, spatial_disc_idx, fspatial_disc_idx;
   int        fnlx_temp, fnly_temp, cnlx_temp, cnly_temp;
   double     cdt, fdt;
   double     tstart, f_tstop, f_tprior, c_tstop, c_tprior;

   int        cilower_x[2], ciupper_x[2];
   int        cnlx, cnly;
   
   /* Get Coarse and fine time step sizes */
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

   /* Determine New sizes, negate ncoarsen because we are refining, i.e., multiplying dx 
    * and dy by 2^(-ncoarsen) */
   retrieve_spatial_discretization(app, cdt, fdt, &spatial_disc_idx );
   ncoarsen       = -(app->spatial_disc_table[spatial_disc_idx]).ncoarsen;
   cnlx           = (app->spatial_disc_table[spatial_disc_idx]).nlx;
   cnly           = (app->spatial_disc_table[spatial_disc_idx]).nly;
   cilower_x[0]   = (app->spatial_disc_table[spatial_disc_idx]).ilower_x[0];
   cilower_x[1]   = (app->spatial_disc_table[spatial_disc_idx]).ilower_x[1];
   ciupper_x[0]   = (app->spatial_disc_table[spatial_disc_idx]).iupper_x[0];
   ciupper_x[1]   = (app->spatial_disc_table[spatial_disc_idx]).iupper_x[1];

   refine_factor = pow(2.0, ncoarsen);
   fnlx = (int) ((cnlx-1) / refine_factor) + 1; 
   fnly = (int) ((cnly-1) / refine_factor) + 1; 
   filower_x[0] = 0;          /* Assume no spatial parallelism */
   filower_x[1] = 0;          /*       ...  */
   fiupper_x[0] = fnlx-1;     /*       ...  */
   fiupper_x[1] = fnly-1;     /*       ...  */
    
   /* Grab values from vector cu*/
   HYPRE_SStructVectorGather( cu->x );
   cvalues = (double *) malloc( (cnlx)*(cnly)*sizeof(double) );
   HYPRE_SStructVectorGetBoxValues( cu->x, 0, cilower_x,
                                    ciupper_x, 0, cvalues );
   
   if (ncoarsen < 0)
   {   
       /* The refinement algorithm just refines by a factor of two, so we 
        * loop over it ncoarsen number of times */
       cnlx_temp = cnlx;
       cnly_temp = cnly;
       fnlx_temp = 2*(cnlx_temp-1) + 1;
       fnly_temp = 2*(cnly_temp-1) + 1;
       
       for(k = 0; k > ncoarsen; k--)
       {
          /* Set the fine values using simple bilinear interpolation */
          fvalues = (double *) malloc( (fnlx_temp)*(fnly_temp)*sizeof(double) );
          counter = 0;
          
          for(i = 0; i < fnlx_temp; i++)
          {
             for(j = 0; j < fnly_temp; j++)
             {
                if( (i%2 == 0) && (j%2 == 0) )
                {
                   /* Injection: this is the F-point analogue to a C-point */
                   fvalues[counter] = cvalues[ (i/2)*cnlx_temp + j/2 ];
                }
                else if( (i%2 == 0) && (j%2 == 1) )
                {
                   /* Interpolate: this is an F-point horizontally between two C-points */
                   fvalues[counter] = 0.5*cvalues[ (i/2)*cnlx_temp + j/2 ] + 
                                      0.5*cvalues[ (i/2)*cnlx_temp + j/2 +1 ];
                }
                else if( (i%2 == 1) && (j%2 == 0) )
                {
                   /* Interpolate: this is an F-point vertically between two C-points */
                   fvalues[counter] = 0.5*cvalues[ (i/2)*cnlx_temp + j/2 ] + 
                                      0.5*cvalues[ ((i/2)+1)*cnlx_temp + j/2 ];
                }
                else if( (i%2 == 1) && (j%2 == 1) )
                {
                   /* Interpolate: this is an F-point in the center of a grid cell */
                   fvalues[counter] = 0.25*cvalues[ (i/2)*cnlx_temp + j/2        ] + 
                                      0.25*cvalues[ (i/2)*cnlx_temp + j/2 +1     ] +
                                      0.25*cvalues[ ((i/2)+1)*cnlx_temp + j/2    ] + 
                                      0.25*cvalues[ ((i/2)+1)*cnlx_temp + j/2 + 1];
                }
                counter += 1;
             }
          }
          
          cnlx_temp = fnlx_temp;
          cnly_temp = fnly_temp;
          fnlx_temp = 2*(cnlx_temp-1) + 1;
          fnly_temp = 2*(cnly_temp-1) + 1;
          if( k-1 > ncoarsen )
          {
             free(cvalues); 
             cvalues = fvalues;
          }
       }

   }
   else
   {
      /* No refinement, just copy the vector */
      fvalues = (double *) malloc( (fnlx)*(fnly)*sizeof(double) );
      for(i = 0; i < fnlx*fnly; i++)
      {
          fvalues[i] = cvalues[i];
      }
   }
   
   
   /* Create an empty vector object. */
   fspatial_disc_idx = (app->spatial_disc_table[spatial_disc_idx]).fspatial_disc_idx;
   fu = (my_Vector *) malloc(sizeof(my_Vector));
   HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[fspatial_disc_idx]).grid_x, &(fu->x) );
   
   /* Set the object type (by default HYPRE_SSTRUCT) and initialize. */
   HYPRE_SStructVectorSetObjectType( fu->x, app->object_type );
   HYPRE_SStructVectorInitialize( fu->x );


   /* Finalize fu vector */
   HYPRE_SStructVectorSetBoxValues( fu->x, 0, filower_x,
                                    fiupper_x, 0, fvalues );
   free( cvalues );
   free( fvalues );
   HYPRE_SStructVectorAssemble( fu->x );


   /* Store the spatial_disc_idx used to generate fu */
   fu->spatial_disc_idx = fspatial_disc_idx; 

   *fu_ptr = fu;
    
   return 0;
}

/* --------------------------------------------------------------------
 * Create a coarsened copy of a vector object.
 * Only uses injection
 * Assume a regular grid of size 2^k + 1 in each dimension
 * Assuming no spatial parallelism
 * -------------------------------------------------------------------- */
int
my_CoarsenInjection(braid_App              app,           
                    braid_Vector           fu,
                    braid_Vector          *cu_ptr,
                    braid_CoarsenRefStatus status)
{
   my_Vector *cu;
   double     *cvalues, *fvalues;
   int        counter, frow, fcol, findex, i, j, spatial_disc_idx;
   int        filower_x[2], fiupper_x[2];
   int        fspatial_disc_idx = fu->spatial_disc_idx;
   int        fnlx = (app->spatial_disc_table[fspatial_disc_idx]).nlx;
   int        fnly = (app->spatial_disc_table[fspatial_disc_idx]).nly;
   double     fdy = (app->spatial_disc_table[fspatial_disc_idx]).dy;
   double     fdx = (app->spatial_disc_table[fspatial_disc_idx]).dx;

   int        cilower_x[2], ciupper_x[2];
   double     coarsen_factor;
   int        cnlx, cnly, ncoarsen;
   double     cdt, fdt;
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
   
   filower_x[0]    = (app->spatial_disc_table[fspatial_disc_idx]).ilower_x[0];
   filower_x[1]    = (app->spatial_disc_table[fspatial_disc_idx]).ilower_x[1];
   fiupper_x[0]    = (app->spatial_disc_table[fspatial_disc_idx]).iupper_x[0];
   fiupper_x[1]    = (app->spatial_disc_table[fspatial_disc_idx]).iupper_x[1];
   
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
   get_coarse_spatial_disc( app, cdt, fdt, fdx, fdy, fnlx, fnly, 
                            fspatial_disc_idx, filower_x, fiupper_x, &spatial_disc_idx);
   ncoarsen = (app->spatial_disc_table[spatial_disc_idx]).ncoarsen;
   cnlx     = (app->spatial_disc_table[spatial_disc_idx]).nlx;
   cnly     = (app->spatial_disc_table[spatial_disc_idx]).nly;
   cilower_x[0] = (app->spatial_disc_table[spatial_disc_idx]).ilower_x[0];
   cilower_x[1] = (app->spatial_disc_table[spatial_disc_idx]).ilower_x[1];
   ciupper_x[0] = (app->spatial_disc_table[spatial_disc_idx]).iupper_x[0];
   ciupper_x[1] = (app->spatial_disc_table[spatial_disc_idx]).iupper_x[1];
   coarsen_factor = pow(2.0, ncoarsen);

   /* Create an empty vector object. */
   HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[spatial_disc_idx]).grid_x, &(cu->x) );
   
   /* Set the object type (by default HYPRE_SSTRUCT) and initialize. */
   HYPRE_SStructVectorSetObjectType( cu->x, app->object_type );
   HYPRE_SStructVectorInitialize( cu->x );

   /* Set the coarse values, assuming a regular grid and injection  */
   cvalues = (double *) malloc( (cnlx)*(cnly)*sizeof(double) );
   HYPRE_SStructVectorGather( fu->x );
   fvalues = (double *) malloc( fnlx*fnly*sizeof(double) );
   HYPRE_SStructVectorGetBoxValues( fu->x, 0, filower_x, fiupper_x, 0, fvalues );
   counter = 0;
   for(i = 0; i < cnlx; i++)
   {
      frow = coarsen_factor*i;
      for(j = 0; j < cnly; j++)
      {
         fcol = coarsen_factor*j;
         findex = frow*fnlx + fcol;
         cvalues[counter] = fvalues[findex];
         counter += 1;
      }
   }
   
   HYPRE_SStructVectorSetBoxValues( cu->x, 0, cilower_x,
                                    ciupper_x, 0, cvalues );
   HYPRE_SStructVectorAssemble( cu->x );
   free( fvalues );
   free( cvalues );

   /* Store spatial_disc_idx */
   cu->spatial_disc_idx = spatial_disc_idx;

   *cu_ptr = cu;

   return 0;
}



/* --------------------------------------------------------------------
 * Create a coarsened copy of a vector object.
 * Uses the transpose of bilinear interpolation scaled by 1/4
 * Assume a regular grid of size 2^k + 1 in each dimension
 * Assuming no spatial parallelism
 * -------------------------------------------------------------------- */
int
my_CoarsenBilinear(braid_App              app,           
                   braid_Vector           fu,
                   braid_Vector          *cu_ptr,
                   braid_CoarsenRefStatus status)
{
   my_Vector *cu;
   double     *cvalues, *fvalues;
   int        counter, i, j, k, spatial_disc_idx;
   int        filower_x[2], fiupper_x[2];
   int        fspatial_disc_idx = fu->spatial_disc_idx;
   int        fnlx = (app->spatial_disc_table[fspatial_disc_idx]).nlx;
   int        fnly = (app->spatial_disc_table[fspatial_disc_idx]).nly;
   double     fdy = (app->spatial_disc_table[fspatial_disc_idx]).dy;
   double     fdx = (app->spatial_disc_table[fspatial_disc_idx]).dx;
   int        cnlx_temp, cnly_temp, fnlx_temp, fnly_temp; 

   int        cilower_x[2], ciupper_x[2];
   int        ncoarsen;
   double     scale = 0.25;
   double     cdt, fdt;
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
   
   filower_x[0]    = (app->spatial_disc_table[fspatial_disc_idx]).ilower_x[0];
   filower_x[1]    = (app->spatial_disc_table[fspatial_disc_idx]).ilower_x[1];
   fiupper_x[0]    = (app->spatial_disc_table[fspatial_disc_idx]).iupper_x[0];
   fiupper_x[1]    = (app->spatial_disc_table[fspatial_disc_idx]).iupper_x[1];
   
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
   get_coarse_spatial_disc( app, cdt, fdt, fdx, fdy, fnlx, fnly, 
                            fspatial_disc_idx, filower_x, fiupper_x, &spatial_disc_idx);
   ncoarsen = (app->spatial_disc_table[spatial_disc_idx]).ncoarsen;
   cilower_x[0] = (app->spatial_disc_table[spatial_disc_idx]).ilower_x[0];
   cilower_x[1] = (app->spatial_disc_table[spatial_disc_idx]).ilower_x[1];
   ciupper_x[0] = (app->spatial_disc_table[spatial_disc_idx]).iupper_x[0];
   ciupper_x[1] = (app->spatial_disc_table[spatial_disc_idx]).iupper_x[1];

   /* Create an empty vector object. */
   HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[spatial_disc_idx]).grid_x, &(cu->x) );
   
   /* Set the object type (by default HYPRE_SSTRUCT) and initialize. */
   HYPRE_SStructVectorSetObjectType( cu->x, app->object_type );
   HYPRE_SStructVectorInitialize( cu->x );

   /* Set the coarse values, assuming a regular grid and injection  */
   fvalues = (double *) malloc( fnlx*fnly*sizeof(double) );
   HYPRE_SStructVectorGetBoxValues( fu->x, 0, filower_x, fiupper_x, 0, fvalues );
   if (ncoarsen > 0)
   {  
      /* The coarsening algorithm just coarsens by a factor of two, so we 
       * loop over it ncoarsen number of times */
      cnlx_temp = (fnlx-1)/2+1;
      cnly_temp = (fnly-1)/2+1;
      fnlx_temp = fnlx; 
      fnly_temp = fnly; 
      
      for(k = 0; k < ncoarsen; k++)
      {
         /* Initialize cvalues */
         cvalues = (double *) malloc( (cnlx_temp)*(cnly_temp)*sizeof(double) );
         for(i=0; i < cnlx_temp*cnly_temp; i++)
         {
            cvalues[i] = 0.0;
         }

         /* Set the coarse values using transpose of simple bilinear interpolation */
         counter = 0;
         for(i = 0; i < fnlx_temp; i++)
         {
            for(j = 0; j < fnly_temp; j++)
            {
               if( (i%2 == 0) && (j%2 == 0) )
               {
                  /* Injection: this is the F-point analogue to a C-point */
                  cvalues[ (i/2)*cnlx_temp + j/2 ] += scale*fvalues[counter];
               }
               else if( (i%2 == 0) && (j%2 == 1) )
               {
                  /* Interpolate: this is an F-point horizontally between two C-points */
                  cvalues[ (i/2)*cnlx_temp + j/2 ]    += scale*0.5*fvalues[counter];
                  cvalues[ (i/2)*cnlx_temp + j/2 +1 ] += scale*0.5*fvalues[counter];                   
               }
               else if( (i%2 == 1) && (j%2 == 0) )
               {
                  /* Interpolate: this is an F-point vertically between two C-points */
                  cvalues[ (i/2)*cnlx_temp + j/2 ]     += scale*0.5*fvalues[counter];
                  cvalues[ ((i/2)+1)*cnlx_temp + j/2 ] += scale*0.5*fvalues[counter];
               }
               else if( (i%2 == 1) && (j%2 == 1) )
               {
                  /* Interpolate: this is an F-point in the center of a grid cell */
                  cvalues[ (i/2)*cnlx_temp + j/2        ] += scale*0.25*fvalues[counter];
                  cvalues[ (i/2)*cnlx_temp + j/2 +1     ] += scale*0.25*fvalues[counter];
                  cvalues[ ((i/2)+1)*cnlx_temp + j/2    ] += scale*0.25*fvalues[counter];
                  cvalues[ ((i/2)+1)*cnlx_temp + j/2 + 1] += scale*0.25*fvalues[counter];
               }
               counter += 1;
            }
         }
         
         fnlx_temp = cnlx_temp;
         fnly_temp = cnly_temp;
         cnlx_temp = (fnlx_temp-1)/2+1;
         cnly_temp = (fnly_temp-1)/2+1;
         if( k+1 < ncoarsen )
         {
            free(fvalues); 
            fvalues = cvalues;
         }
      }

   }
   else
   {
      /* No refinement, just copy the vector */
      cvalues = (double *) malloc( (fnlx)*(fnly)*sizeof(double) );
      for(i = 0; i < fnlx*fnly; i++)
      {
          cvalues[i] = fvalues[i];
      }
   }

   
   HYPRE_SStructVectorSetBoxValues( cu->x, 0, cilower_x,
                                    ciupper_x, 0, cvalues );
   HYPRE_SStructVectorAssemble( cu->x );
   free( fvalues );
   free( cvalues );

   /* Store spatial_disc_idx */
   cu->spatial_disc_idx = spatial_disc_idx;

   *cu_ptr = cu;

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
   
   int ilower_x[2], iupper_x[2];
   int    nlx     = (app->spatial_disc_table[x->spatial_disc_idx]).nlx;
   int    nly     = (app->spatial_disc_table[x->spatial_disc_idx]).nly;

   ilower_x[0]    = (app->spatial_disc_table[x->spatial_disc_idx]).ilower_x[0];
   ilower_x[1]    = (app->spatial_disc_table[x->spatial_disc_idx]).ilower_x[1];
   iupper_x[0]    = (app->spatial_disc_table[x->spatial_disc_idx]).iupper_x[0];
   iupper_x[1]    = (app->spatial_disc_table[x->spatial_disc_idx]).iupper_x[1];

   values_x = (double *) malloc( nlx*nly*sizeof(double) );
   values_y = (double *) malloc( nlx*nly*sizeof(double) );

   HYPRE_SStructVectorGather( x->x );
   HYPRE_SStructVectorGetBoxValues( x->x, part, ilower_x,
                                    iupper_x, var, values_x );

   HYPRE_SStructVectorGather( y->x );
   HYPRE_SStructVectorGetBoxValues( y->x, part, ilower_x,
                                    iupper_x, var, values_y );

   for( i = 0; i < nlx*nly; i++ )
   {
      values_y[i] = alpha*values_x[i] + beta*values_y[i];
   }

   HYPRE_SStructVectorSetBoxValues( y->x, part, ilower_x,
                                    iupper_x, var, values_y );

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
   int        index, myid;
   static int previous_level = -5;
   char       filename[255];
   FILE      *file;
   double     t;

   double    *values;
   /* We have one part and one variable. */
   int part = 0;
   int var = 0;

   /*  damping factor */
   double damping, damping_nt = 1.0; 

   /* error norm */
   double enorm = 0.0;


   /* error vector */
   HYPRE_SStructVector e;

   int i, j, m;
   
   int ilower_x[2], iupper_x[2];
   int    nlx     = (app->spatial_disc_table[u->spatial_disc_idx]).nlx;
   int    nly     = (app->spatial_disc_table[u->spatial_disc_idx]).nly;
   double dx      = (app->spatial_disc_table[u->spatial_disc_idx]).dx;
   double dy      = (app->spatial_disc_table[u->spatial_disc_idx]).dy;

   ilower_x[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower_x[0];
   ilower_x[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower_x[1];
   iupper_x[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper_x[0];
   iupper_x[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper_x[1];
   
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
   if( app->output_files && (level == 0) ){
      if( app->explicit )
         /* forward (explicit) Euler */
         damping = 1 + ((2*(app->K)*(app->dt))/
                        (dx*dx))*(cos(dx)-1) 
                  + ((2*(app->K)*(app->dt))/
                        (dy*dy))*(cos(dy)-1);
      else
         /* backward (implicit) Euler */
         damping = 1.0 / ( 1 + ((2*(app->K)*(app->dt))/
                                (dx*dx))*(1-cos(dx)) 
                             + ((2*(app->K)*(app->dt))/
                                (dy*dy))*(1-cos(dy)));

      

      index = ((t-tstart) / ((tstop-tstart)/ntime) + 0.1);

      /* damping factor after index time steps */
      for( i = 1; i <= index; i++ )
         damping_nt *= damping;

      MPI_Comm_rank(comm, &myid);

      sprintf(filename, "%s.iter%03d.time%07d.proc%05d", "drive-05.out", iter, index, myid);
      file = fopen(filename, "w");

      values = (double *) malloc( nlx*nly*sizeof(double) );
      HYPRE_SStructVectorGetBoxValues( u->x, part, ilower_x, 
                                       iupper_x, var, values );
      m = 0;
      for( j = 0; j < nly; j++ )
         for( i = 0; i < nlx; i++ ){
            enorm += pow(values[m++] - 
                         damping_nt*sin((ilower_x[0]+i)*dx)
                                   *sin((ilower_x[1]+j)*dy),2);

            /*fprintf(file, "%06d %.14e --- discrete solution %.14e\n", 
                    pj*px*nlx*nly+pi*nlx+j*px*nlx+i,
                    values[m++], 
                    damping_nt*sin((ilower_x[0]+i)*dx)
                              *sin((ilower_x[1]+j)*dy));*/
         }
      fprintf(file, "%.14e\n", sqrt(enorm));

      fflush(file);
      fclose(file);
      free( values );
   }

   /* Save the error and solution for GLVis visualization */
   if( app->output_vis && (level == 0) ){ 
      MPI_Comm_rank(app->comm_x, &myid);

      if( t == app->tstop ){
         values = (double *) malloc( nlx*nly*sizeof(double) );
         HYPRE_SStructVectorGetBoxValues( u->x, part, ilower_x, 
                                          iupper_x, var, values );

         HYPRE_SStructVectorCreate( app->comm_x, app->grid_x, &e );
         HYPRE_SStructVectorSetObjectType( e, app->object_type );
         HYPRE_SStructVectorInitialize( e );
      
         m = 0;
         for( j = 0; j < nly; j++ )
            for( i = 0; i < nlx; i++ )
            {
               values[m] = values[m] - 
                           damping_nt*sin((ilower_x[0]+i)*dx)
                                     *sin((ilower_x[1]+j)*dy);
               m++;
            }

         HYPRE_SStructVectorSetBoxValues( e, part, ilower_x,
                                          iupper_x, var, values );
         free( values );
         HYPRE_SStructVectorAssemble( e );

         sprintf(filename, "%s.iter%03d", "drive-05_mesh", iter);
         GLVis_PrintSStructGrid( app->grid_x, filename, 
                                 myid, NULL, NULL );
         sprintf(filename, "%s.iter%03d", "drive-05_err_tstop", iter);
         GLVis_PrintSStructVector( e, 0, filename, myid );
         sprintf(filename, "%s.iter%03d", "drive-05_sol_tstop", iter);
         GLVis_PrintSStructVector( u->x, 0, filename, myid );

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
   int ilower_x[2], iupper_x[2], nlx, nly;
   
   /* Retrieve ilower and iupper */
   ilower_x[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower_x[0];
   ilower_x[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower_x[1];
   iupper_x[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper_x[0];
   iupper_x[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper_x[1];

   /* We have one variable and one part. */
   int        part = 0;
   int        var  = 0;
   
   /* Pack the spatial coarsening rule */ 
   dbuffer[0] = (double) u->spatial_disc_idx; 

   /* Pack the vector */
   HYPRE_SStructVectorGather( u->x );
   HYPRE_SStructVectorGetBoxValues( u->x, part, ilower_x,
                                    iupper_x, var, &(dbuffer[1]) );

   /* Determine number of bytes actually packed */
   nlx = iupper_x[0] - ilower_x[0] + 1;
   nly = iupper_x[1] - ilower_x[1] + 1;
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
   int ilower_x[2], iupper_x[2];
   double    *dbuffer = buffer;
   my_Vector *u;

   /* We have one variable and one part. */
   int        part = 0;
   int        var  = 0;

   u = (my_Vector *) malloc( sizeof(my_Vector) );

   /* Unpack the spatial coarsening rule */
   u->spatial_disc_idx = (int) dbuffer[0];

   /* Retrieve ilower and iupper */
   ilower_x[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower_x[0];
   ilower_x[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).ilower_x[1];
   iupper_x[0]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper_x[0];
   iupper_x[1]    = (app->spatial_disc_table[u->spatial_disc_idx]).iupper_x[1];

   /* 
    * Unpack the vector 
    */

   /* Create an empty vector object. */
   HYPRE_SStructVectorCreate( app->comm_x, (app->spatial_disc_table[u->spatial_disc_idx]).grid_x, &(u->x) );
   
   /* Set the object type (by default HYPRE_SSTRUCT). */
   HYPRE_SStructVectorSetObjectType( u->x, app->object_type );

   /* Indicate that the vector coefficients are ready to be set. */
   HYPRE_SStructVectorInitialize( u->x );

   /* Set the values. */
   HYPRE_SStructVectorSetBoxValues( u->x, part, ilower_x,
                                    iupper_x, var, &(dbuffer[1]) );

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

   int correct;

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
   double dx, dy, dt, c;
   int ilower_x[2], iupper_x[2];

   /* We have one part and one variable.
   int nparts = 1; */
   int nvars = 1;
   int var;

   /* We want to use a struct solver. */
   int object_type = HYPRE_STRUCT;
   
   int sym; 

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
   scoarsen            = 0;
   K                   = 1.0;
   nx                  = 17;
   ny                  = 17;
   nlx                 = 17;
   nly                 = 17;
   tstart              = 0.0;
   nt                  = 32;
   c                   = 0.15;
   sym                 = 0;
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
      else if( strcmp(argv[arg_index], "-c") == 0 ){
          arg_index++;
          c = atof(argv[arg_index++]);
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
      else if ( strcmp(argv[arg_index], "-scoarsen") == 0 ){
         arg_index++;
         scoarsen = atoi(argv[arg_index++]);
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
      printf("  -run_wrapper_tests               : Only run the Braid wrapper tests\n");
      printf("                                     (do not combine with temporal parallelism)\n");
      printf("  -pgrid  <px py pt>               : processors in each dimension (default: 1 1 1)\n");
      printf("  -nx  <nlx nly>                   : 2D spatial problem size of form 2^k+1, 2^k+1 (default: 17 17)\n");
      printf("  -nt  <n>                         : number of time steps (default: 32)\n"); 
      printf("  -c  <c>                          : ratio dt/(dx^2) (default: 1.0)\n"); 
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
      printf("  -scoarsen                        : use spatial coarsening when needed to satisfy CFL\n");
      printf("                                     0 - No spatial coarsening (default) \n");
      printf("                                     1 - Use injection for spatial restriction \n");
      printf("                                     2 - Use transpose of bilinear interpolation for spatial restriction\n");
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
      printf("  -sym <s>                        : symmetric storage (1) or not (0)\n");  
      printf("  -expl <e>                       : use explicit scheme\n");
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
   if( (px > 1) || (py > 1) )
   {
       if( myid == 0 )
            printf("Error: px and py must be 1, i.e. no spatial parallelism.  Sorry!\n");
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
   GetDistribution_x( nx, px, pi, &ilower_x[0], &iupper_x[0] );
   GetDistribution_x( ny, py, pj, &ilower_x[1], &iupper_x[1] );

   /* Determine local problem size. */
   nlx = iupper_x[0] - ilower_x[0] + 1;
   nly = iupper_x[1] - ilower_x[1] + 1;

#if DEBUG
   printf( "%d = (%d, %d): nlx = %d, nly = %d\n", 
           myid, pi, pj, nlx, nly ); 
#endif

   /* Compute grid spacing. */
   dx = PI / (nx - 1);
   dy = PI / (ny - 1);

   /* Set time-step size. */
   dt = K*c*(dx*dx);
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
   (app->sym)             = sym;
   (app->px)              = px;
   (app->py)              = py;
   (app->pi)              = pi;
   (app->pj)              = pj;
   (app->ilower_x[0])     = ilower_x[0];
   (app->ilower_x[1])     = ilower_x[1];
   (app->iupper_x[0])     = iupper_x[0];
   (app->iupper_x[1])     = iupper_x[1];
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
                app->ilower_x, app->iupper_x, app->vartypes );

   /* Define the discretization stencil. */
   set5ptStencil( &(app->stencil), app->sym, app->dim_x );

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

   /* Allocate memory for array of discretization matrices. */
   app->A = (HYPRE_SStructMatrix*) malloc( (app->max_levels)*
                                           sizeof(HYPRE_SStructMatrix));
   /* Create empty matrix lookup table. */
   app->dt_A = (double*) malloc( (app->max_levels)*sizeof(double) );
   for( i = 0; i < app->max_levels; i++ )
   {
      app->dt_A[i] = -1.0;
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
    * spatially */
   (app->spatial_disc_table) = (spatial_discretization*) malloc( max_levels*sizeof(spatial_discretization) );
   for( i = 1; i < app->max_levels; i++ )
   {
      app->spatial_disc_table[i].fdt = -1.0;
      app->spatial_disc_table[i].cdt = -1.0;
   }
   (app->spatial_disc_table[0]).grid_x = (app->grid_x);
   (app->spatial_disc_table[0]).graph = (app->graph);
   (app->spatial_disc_table[0]).fdt = (app->dt);
   (app->spatial_disc_table[0]).cdt = (app->dt);
   (app->spatial_disc_table[0]).ncoarsen = 0;
   (app->spatial_disc_table[0]).dx = (app->dx);
   (app->spatial_disc_table[0]).dy = (app->dy);
   (app->spatial_disc_table[0]).nlx = (app->nlx);
   (app->spatial_disc_table[0]).nly = (app->nly);
   (app->spatial_disc_table[0]).ilower_x[0] = (app->ilower_x[0]);
   (app->spatial_disc_table[0]).ilower_x[1] = (app->ilower_x[1]);
   (app->spatial_disc_table[0]).iupper_x[0] = (app->iupper_x[0]);
   (app->spatial_disc_table[0]).iupper_x[1] = (app->iupper_x[1]);
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

      /* Set cfactor */
      braid_SetCFactor(core, -1, cfactor);
      if( cfactor0 > -1 ){
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

      braid_Drive(core);

      /* Stop timer. */
      myendtime = MPI_Wtime();
      mytime    = myendtime - mystarttime;

      /* Print some additional statistics */
      MPI_Comm_rank( comm, &myid );

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
         printf( "  runtime: %.5lfs\n\n", maxtime );
         printf( "spatial problem size       : %d x %d\n", nx, ny );
         printf( "spatial stopping tolerance : %.2e (loose), %.2e (tight), %.2e (coarse levels)\n", 
                 app->tol_x[0], app->tol_x[1], tol_x_coarse );
         printf( "max spatial iterations     : %d (expensive), %d (cheap)\n", 
                 app->max_iter_x[0], app->max_iter_x[1] );
         printf( "\n" );
         /*printf(" Scheme is either explicit or implicit.  If implicit then\n the table entry is 'impl, %%d' where the integer represents\n the maximum number of AMG iterations used by the implicit\n solver on that level.  'expl' stands for explicit. \n\n"); */
         printf("level   scheme         dx          dy          dt          cfl\n"); 
         printf("-----------------------------------------------------------------\n"); 
         for( i = 0; i < nA_max; i++)
         {
            if( scoarsen_table_global[i*5] == 1)
            {
               printf(" %2d   |  expl        %1.2e    %1.2e    %1.2e    %1.2e\n", i,
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
    * The app->grid_x and app->graph are destroyed above, and are 
    * the i=0 entries below
    */
   for( i = 1; i < app->max_levels; i++ )
   {
      if( (app->spatial_disc_table[i]).fdt != -1.0)
      {
         HYPRE_SStructGridDestroy( (app->spatial_disc_table[i]).grid_x );
         HYPRE_SStructGraphDestroy( (app->spatial_disc_table[i]).graph );
      }
   }
   free( app->spatial_disc_table );
   
   free( app->tol_x );
   free( app );
   MPI_Comm_free( &comm_x );
   MPI_Comm_free( &comm_t );

   /* Finalize MPI */
   MPI_Finalize();

   return 0;
}


