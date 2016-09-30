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


/**
 *  This file contains library functions and basic data structures for ex-03.c 
 *  and ex-03-serial.c.  This file is also used by drivers/driver-diffusion-2D.c.
 *  In general, these files solve the 2D heat equation with the spatial stencil,
 *                                   
 *  For more details on the discretization, see the header comment in ex-03.c.
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "_hypre_utilities.h"
#include "HYPRE_sstruct_ls.h"
#include "_hypre_sstruct_mv.h"

#include "vis.c"

#ifdef M_PI
   #define PI M_PI
#else
   #define PI 3.14159265358979
#endif


/* --------------------------------------------------------------------
 * Simulation manager structure.  
 * Holds the needed simulation data structures, e.g., discretizaion 'e' 
 * or 'i', spatial distribution, and solver used at each time point
 *
 *   comm                communicator 
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
 *   forcing             consider non-zero forcing term
 *   vartype             variable type of of the structured spatial grid 
 *   grid_x              spatial grid
 *   stencil             discretization stencil object
 *   graph               graph object that determine the non-zero structure
 *                       of the discretization matrices
 *   A                   discretization matrix
 *   px                  number of processors in x-dimension
 *   py                  number of processors in y-dimension
 *   pi                  x-coordinate of position in processor grid
 *   pj                  y-coordinate of position in processor grid
 *   ilower              (dim_x)-dimensional array with integer indices of 
 *                       local space interval lower bounds
 *   iupper              (dim_x)-dimensional array with integer indices of
 *                       local space interval upper bounds
 *   object_type         object type of vector to access different hypre solvers 
 *   solver              hypre solver (for implicit time stepping)
 *   max_iter            maximum number of spatial MG iterations
 *   tol                 stopping tolerance for spatial MG
 *   explicit            use explicit discretization (1) or not (0)
 *   output_vis            save the error for GLVis visualization
 *   output_files        save the solution/error/error norm to files
 */
typedef struct _simulation_manager_struct {
   MPI_Comm                comm;
   double                  K;
   int                     dim_x;
   int                     nlx, nly;
   int                     nx, ny;
   double                  tstart;
   double                  tstop;
   int                     nt;
   double                  dx, dy;
   double                  dt;
   int                     forcing;
   HYPRE_SStructVariable   vartype;
   HYPRE_SStructGrid       grid_x;
   HYPRE_SStructStencil    stencil;
   HYPRE_SStructGraph      graph;
   HYPRE_SStructMatrix     A;
   int                     px, py;
   int                     pi, pj;
   int                     ilower[2], iupper[2];
   int                     object_type;
   HYPRE_StructSolver      solver;
   int                     max_iter;
   double                  tol;
   int                     explicit;
   int                     output_vis;
   int                     output_files;
} simulation_manager;

int 
print_simulation_manager(simulation_manager *man)
{
   int myid;
   MPI_Comm_rank( man->comm, &myid );
   
   printf("\n\nmyid:  %d,  Simulation manager contents:\n", myid);
   printf("myid:  %d,  K:            %1.2e\n", myid, man->K);
   printf("myid:  %d,  dim_x         %d\n", myid, man->dim_x);
   printf("myid:  %d,  nlx:          %d\n", myid, man->nlx);
   printf("myid:  %d,  nly:          %d\n", myid, man->nly);
   printf("myid:  %d,  nx:           %d\n", myid, man->ny);
   printf("myid:  %d,  ny:           %d\n", myid, man->ny);
   printf("myid:  %d,  tstart:       %1.2e\n", myid, man->tstart);
   printf("myid:  %d,  tstop:        %1.2e\n", myid, man->tstop);
   printf("myid:  %d,  nt:           %d\n", myid, man->nt);
   printf("myid:  %d,  dx:           %1.2e\n", myid, man->dx);
   printf("myid:  %d,  dy:           %1.2e\n", myid, man->dy);
   printf("myid:  %d,  dt:           %1.2e\n", myid, man->dt);
   printf("myid:  %d,  forcing:      %d\n", myid, man->forcing);
   printf("myid:  %d,  px:           %d\n", myid, man->px);
   printf("myid:  %d,  py:           %d\n", myid, man->py);
   printf("myid:  %d,  pi:           %d\n", myid, man->pi);
   printf("myid:  %d,  pj:           %d\n", myid, man->pj);
   printf("myid:  %d,  ilower[0]     %d\n", myid, man->ilower[0]);
   printf("myid:  %d,  ilower[1]     %d\n", myid, man->ilower[1]);
   printf("myid:  %d,  iupper[0]     %d\n", myid, man->iupper[0]);
   printf("myid:  %d,  iupper[1]     %d\n", myid, man->iupper[1]);
   printf("myid:  %d,  max_iter:     %d\n", myid, man->max_iter);
   printf("myid:  %d,  object_type:  %d\n", myid, man->object_type);
   printf("myid:  %d,  tol:          %1.2e\n", myid, man->tol);
   printf("myid:  %d,  explicit:     %d\n", myid, man->explicit);
   printf("myid:  %d,  output_vis:   %d\n", myid, man->output_vis);
   printf("myid:  %d,  output_files: %d\n", myid, man->output_files);

   printf("\nmyid:  %d,  Note that some object members like vartype, grid_x, stencil, graph, A and solver cannot be printed\n\n", myid);
   return 0;
}

/* --------------------------------------------------------------------
 * Define max, min functions
 * -------------------------------------------------------------------- */
int max_i( int a, int b ){
  return (a >= b ? a : b );
}
int min_i( int a, int b ){
  return (a <= b ? a : b );
}
double max_d( double a, double b ){
  return (a >= b ? a : b );
}
double min_d( double a, double b ){
  return (a <= b ? a : b );
}

/* --------------------------------------------------------------------
 * Grab spatial discretization information from a manager
 * -------------------------------------------------------------------- */
int grab_manager_spatial_info(simulation_manager     *man,
                              int                     ilower[2],
                              int                     iupper[2],
                              int                    *nlx,
                              int                    *nly,
                              int                    *nx,
                              int                    *ny,
                              int                    *px,
                              int                    *py,
                              int                    *pi,
                              int                    *pj,
                              double                 *dx,
                              double                 *dy)
{
   (*nlx) = man->nlx;
   (*nly) = man->nly;
   (*nx)  = man->nx;
   (*ny)  = man->ny;
   (*dx)  = man->dx;
   (*dy)  = man->dy;

   (*px)  = man->px;
   (*py)  = man->py;
   (*pi)  = man->pi;
   (*pj)  = man->pj;

   ilower[0] = man->ilower[0];
   ilower[1] = man->ilower[1];
   iupper[0] = man->iupper[0];
   iupper[1] = man->iupper[1];

   return 0;
}

/* --------------------------------------------------------------------
 * Compute a processor's 2D block of the global 2D grid
 * -------------------------------------------------------------------- */
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
 * Exact solution 
 * -------------------------------------------------------------------- */
double U_exact(simulation_manager *man, double x, double y, double t){
   if(man->forcing)
      return sin(x)*sin(y)*cos(t);
   else
      return exp(-2.*t)*sin(x)*sin(y);
}

/* --------------------------------------------------------------------
 * Boundary condition: zero Dirichlet condition for now
 * -------------------------------------------------------------------- */
double B0(double x, double y){
    return 0.0;
}

/* --------------------------------------------------------------------
 * Forcing term: b(x,y,t) = -sin(x)*sin(y)*( sin(t) - 2*K*cos(t) )
 * -------------------------------------------------------------------- */
double Ft(double x, double y, double t, double K){
   return (-1.0)*sin(x)*sin(y)*( sin(t) - 2.0*K*cos(t) );
}

/* --------------------------------------------------------------------
 * Bundle together three common calls that initialize a vector
 * -------------------------------------------------------------------- */
int initialize_vector(simulation_manager   *man,
                      HYPRE_SStructVector  *u)
{
   HYPRE_SStructVectorCreate( man->comm, man->grid_x, u );
   HYPRE_SStructVectorSetObjectType( *u, man->object_type );
   HYPRE_SStructVectorInitialize( *u );
   return 0;
}

/* --------------------------------------------------------------------
 * Set initial condition in a vector for a time t
 *  t == 0 : results in the initial condition
 *  t > 0  : results in a zero vector
 *  t < 0  : results in a uniformly random vector
 * -------------------------------------------------------------------- */
int set_initial_condition(simulation_manager   *man,
                          HYPRE_SStructVector  *u,
                          double                t)
{
   double *values; 
   int i, j, m;

   initialize_vector(man, u);
   values = (double *) malloc( (man->nlx)*(man->nly)*sizeof(double) );
   
   if( t == 0.0){
      /* Set the values in left-to-right, bottom-to-top order. */ 
      for( m = 0, j = 0; j < man->nly; j++ )
         for (i = 0; i < man->nlx; i++, m++)
            values[m] = U0( ((man->ilower[0])+i)*(man->dx), 
                            ((man->ilower[1])+j)*(man->dy) );
   }
   else if (t < 0.0){
      for( m = 0; m < (man->nlx)*(man->nly); m++ )
            values[m] = ((double)rand())/RAND_MAX;
   }
   else{
      for( m = 0; m < (man->nlx)*(man->nly); m++ )
            values[m] = 0.0;
   }

   HYPRE_SStructVectorSetBoxValues( *u, 0, man->ilower, man->iupper, 0, values ); 
   HYPRE_SStructVectorAssemble( *u );
   free(values);
   return 0;
}

/* --------------------------------------------------------------------
 * Apply boundary conditions for explicit scheme:
 * Put the boundary conditions in the vector. 
 * -------------------------------------------------------------------- */
void 
addBoundary( simulation_manager  *man,
             HYPRE_SStructVector  b) 
{
   double     dx, dy;
   int        ilower[2], iupper[2], nlx, nly, nx, ny, px, py, pi, pj;    

   int i, j;
   
   int bc_ilower[2];
   int bc_iupper[2];
   double *bvalues;

   /* Grab info from manager */
   grab_manager_spatial_info(man, ilower, iupper, &nlx, &nly, &nx, 
                             &ny, &px, &py, &pi, &pj, &dx, &dy);
       
   /* Allocate vector for values on boundary planes */
   bvalues = (double *) malloc( (max_i(nlx, nly)+1)*sizeof(double) );     
  
       
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
              
         HYPRE_SStructVectorSetBoxValues(b, 0, bc_ilower,
                                         bc_iupper, 0, bvalues);
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
              
         HYPRE_SStructVectorSetBoxValues(b, 0, bc_ilower,
                                         bc_iupper, 0, bvalues);
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
              
         HYPRE_SStructVectorSetBoxValues(b, 0, bc_ilower,
                                         bc_iupper, 0, bvalues);
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
              
         HYPRE_SStructVectorSetBoxValues(b, 0, bc_ilower,
                                         bc_iupper, 0, bvalues);
      }
   }

   free(bvalues);

   /* Finalize the vector assembly. */
/*   HYPRE_SStructVectorAssemble(b);*/
}


/* --------------------------------------------------------------------
 * Apply boundary conditions for implicit scheme:
 * To incorporate the boundary conditions, we removed the connections 
 * between the interior and boundary nodes in the discretization matrix. 
 * We adjust for removing these connections by appropriately modifying
 * the corresponding RHS entries. 
 * -------------------------------------------------------------------- */
void 
addBoundaryToRHS( simulation_manager  *man,
                  HYPRE_SStructVector  b )
{
   double     dt          = man->dt;
   int        K           = man->K;
   double     dx, dy;
   int        ilower[2], iupper[2], nlx, nly, nx, ny, px, py, pi, pj;    
   int i, j, m;
   
   int bc_ilower[2];
   int bc_iupper[2];
   int istart, iend, jstart, jend;
   double *bvalues;

   /* Grab info from manager */
   grab_manager_spatial_info(man, ilower, iupper, &nlx, &nly, &nx, 
                             &ny, &px, &py, &pi, &pj, &dx, &dy);
   
   /* Allocate vector for values on boundary planes */
   bvalues = (double *) malloc( (max_i(nlx, nly)+1)*sizeof(double) );     
       
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
              
         HYPRE_SStructVectorSetBoxValues(b, 0, bc_ilower,
                                         bc_iupper, 0, bvalues);
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
              
         HYPRE_SStructVectorSetBoxValues(b, 0, bc_ilower,
                                         bc_iupper, 0, bvalues);
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
              
         HYPRE_SStructVectorSetBoxValues(b, 0, bc_ilower,
                                         bc_iupper, 0, bvalues);
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
              
         HYPRE_SStructVectorSetBoxValues(b, 0, bc_ilower,
                                         bc_iupper, 0, bvalues);
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

         HYPRE_SStructVectorAddToBoxValues(b, 0, bc_ilower,
                                           bc_iupper, 0, bvalues);
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
         
         HYPRE_SStructVectorAddToBoxValues(b, 0, bc_ilower,
                                           bc_iupper, 0, bvalues);
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
         
         HYPRE_SStructVectorAddToBoxValues(b, 0, bc_ilower,
                                           bc_iupper, 0, bvalues);
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
         
         HYPRE_SStructVectorAddToBoxValues(b, 0, bc_ilower,
                                           bc_iupper, 0, bvalues);
      }
   }
   free(bvalues);

   /* Finalize the vector assembly. */
/*   HYPRE_SStructVectorAssemble(b);*/
}

/* --------------------------------------------------------------------
 * Add forcing term:
 * We have to multiply the RHS of the PDE by dt.
 * -------------------------------------------------------------------- */
void
addForcingToRHS( simulation_manager *man, 
                 double               t,
                 HYPRE_SStructVector  b )                
{
   double     dt          = man->dt;
   int        K           = man->K;
   double     dx, dy;
   int        ilower[2], iupper[2], nlx, nly, nx, ny, px, py, pi, pj;  
   double    *values;
   int        i, j, m;
   int        rhs_ilower[2];
   int        rhs_iupper[2];
   int        istart, iend, jstart, jend;
   
   /* Grab info from manager */
   grab_manager_spatial_info(man, ilower, iupper, &nlx, &nly, &nx, 
                             &ny, &px, &py, &pi, &pj, &dx, &dy);
   
   /* Add the values from the RHS of the PDE in left-to-right,
    * bottom-to-top order. Entries associated with DoFs on the boundaries
    * are not considered so that boundary conditions are not messed up. */
   
   values = (double *) malloc( nlx*nly*sizeof(double) );
   
   rhs_ilower[0] = ilower[0];
   rhs_ilower[1] = ilower[1];
   istart        = 0;
   iend          = nlx;
   
   rhs_iupper[0] = iupper[0];
   rhs_iupper[1] = iupper[1];
   jstart        = 0;
   jend          = nly;
   
   /* Adjust box to not include boundary nodes. */
   
   /* a) Boundaries y = 0 or y = PI */
   /*    i) Processors at y = 0 */
   if( pj == 0 )
   {
      rhs_ilower[1] += 1;
      jstart        += 1;
   }
   /*    ii) Processors at y = PI */
   if( pj == py-1 )
   {
      rhs_iupper[1] -= 1;
      jend          -= 1;
   }
   
   /* b) Boundaries x = 0 or x = PI */
   /*    i) Processors at x = 0 */
   if( pi == 0 )
   {
      rhs_ilower[0] += 1;
      istart        += 1;
   }
   /*    ii) Processors at x = PI */
   if( pi == px-1 )
   {
      rhs_iupper[0] -= 1;
      iend          -= 1;
   }
   
   for( m = 0, j = jstart; j < jend; j++ )
      for( i = istart; i < iend; i++, m++ )
         values[m] = dt*Ft( (rhs_ilower[0]+i-istart)*dx,
                            (rhs_ilower[1]+j-jstart)*dy, t, K );
   
   HYPRE_SStructVectorAddToBoxValues(b, 0, rhs_ilower,
                                     rhs_iupper, 0, values);
   
/*   HYPRE_SStructVectorAssemble( b );*/
   free(values);
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
             HYPRE_SStructVariable  vartype,
             int                    nghost)
{
   /* We have one part and one variable. */
   int num_ghost[4];
   
   num_ghost[0] = nghost;
   num_ghost[1] = nghost;
   num_ghost[2] = nghost;
   num_ghost[3] = nghost;

   HYPRE_SStructGrid grid;

   /* Create an empty 2D grid object. */
   HYPRE_SStructGridCreate( comm, ndim, 1, &grid );
   
   /* Set the variable type for each part 
    * This call MUST go before setting the number of ghost */
   HYPRE_SStructGridSetVariables( grid, 0, 1, &vartype );

   /* Add a new box to the grid. */
   if ((ilower[0] <= iupper[0]) && (ilower[1] <= iupper[1]))
   {
      HYPRE_SStructGridSetExtents( grid, 0, ilower, iupper );
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
   int i;

   HYPRE_SStructStencil stencil;

   /* Create an empty 2D, 5-pt stencil object. */
   HYPRE_SStructStencilCreate( ndim, 5, &stencil );

   /* Define the geometry of the stencil. */
   int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};

   /* Assign stencil entries. */
   for( i = 0; i < 5; i++ )
      HYPRE_SStructStencilSetEntry( stencil, i, offsets[i], 0 );

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

   HYPRE_SStructGraph graph;

   /* Create the graph object. */
   HYPRE_SStructGraphCreate( comm, grid, &graph );

   /* Set the object type. */
   HYPRE_SStructGraphSetObjectType( graph, object_type );

   /* Now we need to tell the graph which stencil to use for each
    * variable on each part (we only have one variable and one part) */
   HYPRE_SStructGraphSetStencil( graph, 0, 0, stencil );

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
setUpImplicitMatrix( simulation_manager *man)
{
   MPI_Comm             comm        = man->comm;
   HYPRE_SStructGraph   graph       = man->graph;
   int                  object_type = man->object_type;
   double               dt          = man->dt;
   int                  K           = man->K;
   
   double               dx, dy;
   int                  ilower[2], iupper[2], nlx, nly, nx, ny, px, py, pi, pj;  

   int i, j, m, idx;

   int stencil_indices[5] = {0, 1, 2, 3, 4};
   int nentries;

   double *values, *bvalues;
   int bc_ilower[2];
   int bc_iupper[2];

   HYPRE_SStructMatrix A;

   /* Grab info from manager */
   grab_manager_spatial_info(man, ilower, iupper, &nlx, &nly, &nx, 
                             &ny, &px, &py, &pi, &pj, &dx, &dy);

   /* Create an empty matrix object. */
   HYPRE_SStructMatrixCreate( comm, graph, &A);

   /* Use symmetric storage? The function below is for symmetric stencil
    * entries (use HYPRE_SStructMatrixSetNSSymmetric for non-stencil 
    * entries). */
   HYPRE_SStructMatrixSetSymmetric( A, 0, 0, 0, 0);

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

      HYPRE_SStructMatrixSetBoxValues( A, 0, ilower, iupper, 0, 
                                       nentries, stencil_indices, 
                                       values );

      free(values);

      /* 2. correct stencils at boundary nodes */     
      /* Allocate vectors for values on boundary planes */
      values  = (double *) malloc( nentries*( max_i(nlx,nly)+1 )*sizeof(double) );
      bvalues = (double *) malloc( (max_i(nlx,nly)+1)*sizeof(double) );     
      for( i = 0; i < nentries*max_i(nlx,nly); i+= nentries ){
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
            HYPRE_SStructMatrixSetBoxValues(A, 0, bc_ilower, bc_iupper,
                                            0, nentries,
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
            HYPRE_SStructMatrixSetBoxValues(A, 0, bc_ilower, bc_iupper,
                                            0, nentries,
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
            HYPRE_SStructMatrixSetBoxValues(A, 0, bc_ilower, bc_iupper,
                                            0, nentries,
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
            HYPRE_SStructMatrixSetBoxValues(A, 0, bc_ilower, bc_iupper,
                                            0, nentries,
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
            HYPRE_SStructMatrixSetBoxValues(A, 0, bc_ilower, bc_iupper,
                                            0, 1, stencil_indices, bvalues);
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
            HYPRE_SStructMatrixSetBoxValues(A, 0, bc_ilower, bc_iupper,
                                            0, 1, stencil_indices, bvalues);
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
            HYPRE_SStructMatrixSetBoxValues(A, 0, bc_ilower, bc_iupper,
                                            0, 1, stencil_indices, bvalues);
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
            HYPRE_SStructMatrixSetBoxValues(A, 0, bc_ilower, bc_iupper,
                                            0, 1, stencil_indices, bvalues);
         }
      }

      free(bvalues);
   }
   /* Finalize the matrix assembly. */
   HYPRE_SStructMatrixAssemble( A );

   man->A = A;
}



/* --------------------------------------------------------------------
 * Set up the explicit discretization matrix. 
 * First, we set the stencil values at every node neglecting the 
 * boundary. Then, we correct the matrix stencil at boundary nodes.
 * We have to eliminate the coefficients reaching outside of the domain
 * boundary. 
 * -------------------------------------------------------------------- */
void
setUpExplicitMatrix( simulation_manager   *man )
{

   MPI_Comm             comm        = man->comm;
   HYPRE_SStructGraph   graph       = man->graph;
   int                  object_type = man->object_type;
   double               dt          = man->dt;
   int                  K           = man->K;
   
   double               dx, dy;
   int                  ilower[2], iupper[2], nlx, nly, nx, ny, px, py, pi, pj;  

   int i, j, m, idx;
   int stencil_indices[5] = {0, 1, 2, 3, 4};
   int nentries;

   double *values;
   int bc_ilower[2];
   int bc_iupper[2];

   HYPRE_SStructMatrix A;
   
   /* Grab info from manager */
   grab_manager_spatial_info(man, ilower, iupper, &nlx, &nly, &nx, 
                             &ny, &px, &py, &pi, &pj, &dx, &dy);

   /* Create an empty matrix object. */
   HYPRE_SStructMatrixCreate( comm, graph, &A);

   HYPRE_SStructMatrixSetSymmetric( A, 0, 0, 0, 0);

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

      HYPRE_SStructMatrixSetBoxValues( A, 0, ilower, iupper, 0, 
                                       nentries, stencil_indices, 
                                       values );

      free(values);

      /* 2. correct stencils at boundary nodes */     
      /* Allocate vectors for values on boundary planes */
      values  = (double *) malloc( nentries*( max_i(nlx,nly)+1 )*
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
            HYPRE_SStructMatrixSetBoxValues(A, 0, bc_ilower, bc_iupper,
                                            0, nentries, stencil_indices, values);
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
            HYPRE_SStructMatrixSetBoxValues(A, 0, bc_ilower, bc_iupper,
                                            0, nentries, stencil_indices, values);
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
            HYPRE_SStructMatrixSetBoxValues(A, 0, bc_ilower, bc_iupper,
                                            0, nentries, stencil_indices, values);
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
            HYPRE_SStructMatrixSetBoxValues(A, 0, bc_ilower, bc_iupper,
                                            0, nentries, stencil_indices, values);
         }
      }

      free(values);   
   }

   /* Finalize the matrix assembly. */
   HYPRE_SStructMatrixAssemble( A );
   man->A = A;
}


/* --------------------------------------------------------------------
 * Set up PFMG solver.
 * -------------------------------------------------------------------- */
void
setUpStructSolver( simulation_manager  *man,
                   HYPRE_SStructVector  b,
                   HYPRE_SStructVector  x )                
{
   MPI_Comm            comm     = man->comm;
   HYPRE_SStructMatrix A        = man->A;
   int                 max_iter = man->max_iter;
   double              tol      = man->tol;
   
   /* hard coded PFMG parameters */
   int n_pre               = 1;       /* number of PFMG presmoothing steps */
   int n_post              = 1;       /* number of PFMG postmoothing steps */
   int rap                 = 1;       /* type of RAP coarse grid to use in PFMG
                                         0 - Galerkin (default)
                                         1 - non-Galerkin ParFlow operators
                                         2 - Galerkin, general operators */
   int relax               = 3;       /* type of relaxation to use in PFMG 
                                         0 - Jacobi
                                         1 - Weighted Jacobi
                                         2 - R/B Gauss-Seidel (default)
                                         3 - R/B Gauss-Seidel (nonsymmetric) */
   int skip                = 1;       /* PFMG solver option (see print usage block) */



   HYPRE_StructSolver  solver;
   HYPRE_StructMatrix  sA;
   HYPRE_StructVector  sb;
   HYPRE_StructVector  sx;

   HYPRE_SStructMatrixGetObject( A, (void **) &sA );
   HYPRE_SStructVectorGetObject( b, (void **) &sb );
   HYPRE_SStructVectorGetObject( x, (void **) &sx );

   /* Set PFMG options. */
   HYPRE_StructPFMGCreate( comm, &solver );
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

   man->solver = solver;
}


/* --------------------------------------------------------------------
 * Time integration routine
 * -------------------------------------------------------------------- */
int take_step(simulation_manager * man,         /* manager holding basic sim info */
              HYPRE_SStructVector  xstop,       /* approximation at tstop */
              HYPRE_SStructVector  bstop,       /* additional RHS forcing */
              HYPRE_SStructVector  x,           /* vector to evolve */
              double               tstart,      /* evolve x from tstart to tstop */
              double               tstop, 
              int                 *iters_taken) /* if implicit, returns the number of iters taken */
{
   int      iters      = man->max_iter; /* if implicit, max iters for solve */
   double   tol        = man->tol;      /* if implicit, solver tolerance */
   int      explicit   = man->explicit; /* if true, use explicit, else implicit */
   int      forcing    = man->forcing;  /* if true, use the nonzero forcing term */
   HYPRE_Int num_iters = 0;
   
   HYPRE_SStructVector b;
   HYPRE_StructMatrix  sA;
   HYPRE_StructVector  sxstop, sbstop, sx, sb;
   
   /* Grab these object pointers for use below */
   HYPRE_SStructMatrixGetObject( man->A, (void **) &sA );
   HYPRE_SStructVectorGetObject( xstop, (void **) &sxstop );
   HYPRE_SStructVectorGetObject( x, (void **) &sx );

   /* Create temporary vector */
   initialize_vector(man, &b);
   HYPRE_SStructVectorAssemble(b);
   HYPRE_SStructVectorGetObject( b, (void **) &sb );
   HYPRE_StructVectorCopy(sx, sb);

   if( explicit )
   {
      /* Incorporate the boundary conditions. */    
      addBoundary( man, b );

      /* Time integration to next time point: Perform MatVec x = Ab. */
      HYPRE_StructMatrixMatvec( 1, sA, sb, 0, sx );

      if (forcing) {
         /* add RHS of PDE: g_i = dt*b_{i-1}, i > 0 */
         addForcingToRHS( man, tstop, x );
      }
      if (bstop != NULL)
      {
         /* Add extra forcing from braid */
         HYPRE_SStructVectorGetObject( bstop, (void **) &sbstop );
         hypre_StructAxpy(1.0, sbstop, sx);
      }
   }
   else
   {
      /* Set up the right-hand side vector, which is the solution from 
       * the previous time step modified to incorporate the boundary 
       * conditions and the right-hand side of the PDE */ 
      addBoundaryToRHS( man, b );
      
      if (forcing) {
         /* add RHS of PDE: g_i = Phi*dt*b_i, i > 0 */
         addForcingToRHS( man, tstop, b );
      }
      if (bstop != NULL)
      {
         /* Add extra forcing from braid */
         HYPRE_SStructVectorGetObject( bstop, (void **) &sbstop );
         hypre_StructAxpy(1.0, sbstop, sb);
      }

      /* Solve system */
      if (xstop != x)
      {
         /* Set initial guess */
         HYPRE_StructVectorCopy(sxstop, sx);
      }
      HYPRE_StructPFMGSetTol( man->solver, tol );
      HYPRE_StructPFMGSetMaxIter( man->solver, iters);
      HYPRE_StructPFMGSolve( man->solver, sA, sb, sx );
      HYPRE_StructPFMGGetNumIterations( man->solver, &num_iters);
      (*iters_taken) = num_iters;

   }

   /* free memory */
   HYPRE_SStructVectorDestroy( b );

   return 0;
}

/* --------------------------------------------------------------------
 * Residual routine
 * -------------------------------------------------------------------- */
int comp_res(simulation_manager * man,         /* manager holding basic sim info */
             HYPRE_SStructVector  xstop,       /* approximation at tstop */
             HYPRE_SStructVector  r,           /* approximation at tstart */
             double               tstart,
             double               tstop)
{
   int      explicit   = man->explicit; /* if true, use explicit, else implicit */
   int      forcing    = man->forcing;  /* if true, use the nonzero forcing term */
   
   HYPRE_SStructVector b;
   HYPRE_StructMatrix  sA;
   HYPRE_StructVector  sxstop, sr, sb; 
   
   /* Grab these object pointers for use below */
   HYPRE_SStructMatrixGetObject( man->A, (void **) &sA );
   HYPRE_SStructVectorGetObject( xstop, (void **) &sxstop );
   HYPRE_SStructVectorGetObject( r, (void **) &sr );

   if( explicit )
   {
   
      /* Create temporary vector */
      initialize_vector(man, &b);
      HYPRE_SStructVectorAssemble(b);
      HYPRE_SStructVectorGetObject( b, (void **) &sb );
      HYPRE_StructVectorCopy(sr, sb);

      /* Incorporate the boundary conditions. */    
      addBoundary( man, b );

      /* Residual r = xstop - A*r - forcing */
      HYPRE_StructMatrixMatvec( 1, sA, sb, 0, sr );
      if (forcing) {
         /* add RHS of PDE: g_i = dt*b_{i-1}, i > 0 */
         addForcingToRHS( man, tstop, r );
      }
      hypre_StructAxpy (-1, sxstop, sr);
      hypre_StructScale(-1, sr);

      /* free memory */
      HYPRE_SStructVectorDestroy( b );

   }
   else
   {
      /* Set up the right-hand side vector, which is the solution from 
       * the previous time step modified to incorporate the boundary 
       * conditions and the right-hand side of the PDE */ 
      addBoundaryToRHS( man, r );

      /* Residual r = A*xstop - r - forcing */
      if (forcing) {
         /* add RHS of PDE: g_i = Phi*dt*b_i, i > 0 */
         addForcingToRHS( man, tstop, r );
      }
      HYPRE_StructMatrixMatvec( -1, sA, sxstop, 1, sr );
      hypre_StructScale(-1, sr);
   }

   return 0;
}


/* --------------------------------------------------------------------
 * Compute the current error vector, relative to the true continuous 
 * solution.  Assumes that e has already been allocated and setup
 * -------------------------------------------------------------------- */
int compute_error(simulation_manager  *man, 
                  HYPRE_SStructVector  x, 
                  double               t,
                  HYPRE_SStructVector  e) 

{
   double  *values;
   int i, j, m;
   

   values = (double *) malloc( (man->nlx) * (man->nly)*sizeof(double) );
   HYPRE_SStructVectorGetBoxValues( x, 0, man->ilower, man->iupper, 0, values );
   
   /* Compute error */
   m = 0;
   for( j = 0; j < man->nly; j++ ){
      for( i = 0; i < man->nlx; i++ ){
         values[m] =  U_exact(man, 
                              ((man->ilower[0])+i)*(man->dx), 
                              ((man->ilower[1])+j)*(man->dy), t) 
                      - values[m];
         m++;
      }
   }
   HYPRE_SStructVectorSetBoxValues( e, 0, man->ilower, man->iupper, 0, values );
   free( values );
   HYPRE_SStructVectorAssemble( e );

   return 0;
}


/* --------------------------------------------------------------------
 * Compute || x ||_2 and return in norm_ptr
 * The 2-norm (Euclidean norm) is used.
 * -------------------------------------------------------------------- */
int norm(HYPRE_SStructVector  x, 
         double              *norm_ptr)
{
   double dot;
   hypre_SStructInnerProd( x, x, &dot );
   *norm_ptr = sqrt(dot);
   return 0;
}


/* --------------------------------------------------------------------
 * Compute the little l2 norm of the discretization error
 * -------------------------------------------------------------------- */
int compute_disc_err(simulation_manager  *man, 
                     HYPRE_SStructVector  u, 
                     double               tstop, 
                     HYPRE_SStructVector  e, 
                     double              *disc_err)
{
   /* be sure to scale by mesh size so that this is the little l2 norm */
   compute_error(man, u, tstop, e);
   norm(e, disc_err);
   *disc_err = sqrt( (*disc_err)*(*disc_err)*(man->dx)*(man->dy) );
   return 0;
}

/* --------------------------------------------------------------------
 * If Proc 0, output the l2 norm of the discretization error at 
 * time t for vector x.  
 * -------------------------------------------------------------------- */
int output_error_file(simulation_manager * man, 
                      double               t,
                      double               err_norm,
                      char                *filename)
{
   FILE  *file;
   int   myid;

   MPI_Comm_rank(man->comm, &myid);
   if(myid == 0)
   {  
      /* Print max error at each time step */
      file = fopen(filename, "w");
      fprintf(file, "%.14e\n", err_norm);
      fflush(file);
      fclose(file);
   }
   
   return 0;
}


/* --------------------------------------------------------------------
 * User a visualization file of the discretization error at time t 
 * for vector x
 * -------------------------------------------------------------------- */
int output_vis(simulation_manager * man, 
               HYPRE_SStructVector  x, 
               double               t,
               char *               filename_mesh, 
               char *               filename_err, 
               char *               filename_sol) 
{   
   int myid;
   HYPRE_SStructVector e;

   MPI_Comm_rank(man->comm, &myid);

   /* compute error */
   initialize_vector(man, &e);
   compute_error(man, x, t, e);
   HYPRE_SStructVectorAssemble( e );
   
   GLVis_PrintSStructGrid( man->grid_x, filename_mesh, myid, NULL, NULL );
   GLVis_PrintSStructVector( e, 0, filename_err, myid );
   GLVis_PrintSStructVector( x, 0, filename_sol, myid );

   HYPRE_SStructVectorDestroy( e );

   return 0;
}


