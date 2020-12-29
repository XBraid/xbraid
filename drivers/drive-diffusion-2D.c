/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. For support, post issues to the XBraid Github page.
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
 * Driver:        drive-diffusion-2D.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make drive-diffusion-2D
 *
 * Help with:     drive-diffusion-2D -help
 *
 * Sample run:    drive-diffusion-2D -nt 100
 *
 * Description:   Solves the 2D heat equation on a regular grid in space and time
 *                Uses forward or backward Euler in time and classic second order 
 *                finite-differencing in space.  Spatial coarsening is also 
 *                supported with linear interpolation and restriction in space. 
 *
 *                More specifically, the diffusion problem is 
 *
 *                                u_t - K * div u = b,
 *                      b(x,y,t) = -sin(x)*sin(y)*(sin(t)-2*cos(t))
 *                                        or
 *                                   b(x,y,t) = 0
 *
 *                in the unit square subject to zero Dirichlet boundary
 *                conditions,u_b = 0, and initial condition U0 at time t = 0.
 *                
 *                The domain is split into an p_x x p_y processor grid.  Each
 *                processor has a n_x x n_y grid, with nodes connected by a
 *                5-point stencil. More precisely, we use central FD in space
 *                and forward or backward Euler in time, definining the stencil
 *
 *                                  -K*dt/(dy^2)
 *                 -K*dt/(dx^2)  1+2*K*(dt/(dx^2)+dt/(dy^2)) -K*dt/(dx^2)
 *                                  -K*dt/(dy^2)   
 *                   
 *                We use cell-centered variables, and, therefore, the nodes are
 *                not shared.
 *
 *                To incorporate the boundary conditions, we do the following:
 *                Let x_i and x_b be the interior and boundary parts of the
 *                solution vector x, respectively, and let u_b the boundary
 *                condition. If we split the matrix A as
 *
 *                           A = [A_ii A_ib; 
 *                                A_bi A_bb],
 *
 *                then we solve
 *
 *                           [A_ii 0; [x_i;    [b_i - A_ib u_b;
 *                            0    I]  x_b] =        u_b       ].
 *
 *                For zero boundary conditions, u_b = 0, we are just 
 *                solving
 *
 *                        A_ii x_i = b_i.
 *
 *                Note that this approach is useful for more general types 
 *                of boundary conditions.
 *
 *                We use a structured solver for the spatial solve at each 
 *                time step.
 *
 *
 * Notes:         Uses the SStructured interface (SStruct) in hypre to build the spatial
 *                discretization matrices, and the hypre PFMG solver for implicit time 
 *                stepping
 **/


/*
   Example 05
   2D diffusion problem

   Interface:    SStructured interface (SStruct)

   Compile with: make drive-02

   Help with:    drive-02 -help

   Sample run:   mpirun -np 8 drive-02 -pgrid 1 1 8 -ml 15 -nt 128 -nx 33 33 -mi 100 -expl -scoarsen  

   Notes:        The " -expl and -scoarsen " options should be used together
                 to coarsen spatially and do explicit time stepping.  

   Description:  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "_hypre_utilities.h"
#include "HYPRE_sstruct_ls.h"
#include "_hypre_sstruct_mv.h"

#include "braid.h"
#include "braid_test.h"

#include "../examples/ex-03-lib.c"

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
 * fspatial_disc_idx    integer index into app->spatial_lookup_table for the 
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
/* struct my_App contains general information about the problem specific
 * to running the problem with XBraid.
 *
 *   comm_t                communicator for parallelizing in time
 *   comm_x                communicator for parallelizing in space 
 *   max_levels            maximum number of time levels
 *   e                     vector used for error computations
 *   A                     array of discretization matrices (one per time level)
 *   dt_A                  array of time steps for which discretization matrix
 *                         has been created
 *   dx_A                  array of x-direction mesh sizes for which a discretization 
 *                         matrix has been created
 *   dy_A                  array of y-direction mesh sizes for which a discretization 
 *                         matrix has been created
 *   nA                    number of discretization matrices that have been
 *                         created
 *   solver                array of solvers used at each time step on different
 *                         time levels
 *   runtime_max_iter      runtime information on the number of iterations taken for 
 *                         the implicit solves
 *   max_iter_x            expensive and cheap maximum number of spatial MG iterations
 *   scoarsen              use spatial refinement and coarsening
 *   num_scoarsen_CFL      the number of entries in scoarsen_CFL
 *   scoarsenCFL           array of CFL numbers such that if the actual CFL at level k is 
 *                         greater that scoarsenCFL[k] then do spatial coarsening until the
 *                         CFL < scoarsenCFL[k]
 *                         Note: scoarsen must also > 0 to turn on spatial coarsening
 *   spatial_lookup_table  Lookup table recording dx, dy when coarsening spatially.
 *                         tuples (cdt, fdt) are the keys used to access the spatial
 *                         discretizations stored in the table.  The accessor function
 *                         is get_coarse_spatial_disc( ).  See it for more documentation.
 *   runtime_scoarsen_info Runtime information on CFL's encountered and spatial 
 *                         discretizations used.
 *   tol_x[2]              Spatial stopping tolerance limits, (loose, tight), for the fine level
 *   use_rand              binary, use random initial guess (1) or zero initial guess (0)
 *   buffer_size           integer containing the largest possible MPI buffer 
 *   print_level           user-desired print level, print level 2 does lots of debug output in Step()
 */
typedef struct _braid_App_struct {
   MPI_Comm                comm_t;
   MPI_Comm                comm_x;
   simulation_manager     *man;
   int                     max_levels;
   HYPRE_SStructVector     e;
   HYPRE_SStructMatrix    *A;
   double                 *dt_A;
   double                 *dx_A;
   double                 *dy_A;
   int                     nA;
   HYPRE_StructSolver     *solver;
   int                    *runtime_max_iter;
   int                    *max_iter_x;
   int                     scoarsen;
   int                     num_scoarsenCFL;
   double                 *scoarsenCFL;
   spatial_discretization *spatial_lookup_table;
   double                 *runtime_scoarsen_info;
   double                  tol_x[2];
   int                     use_rand;
   int                     buffer_size;
   int                     print_level;
} my_App;

/* struct my_Vector contains local information specific to one time point
 *   x            spatial vector 
 */
typedef struct _braid_Vector_struct
{
   int                   spatial_disc_idx;
   HYPRE_SStructVector   x;
} my_Vector;

/* --------------------------------------------------------------------
 * Inititialization routine for the user's manager structure, so that it
 * aligns with a given vector and spatial discretization level
 * -------------------------------------------------------------------- */
   
int update_manager_from_vector(simulation_manager     *man, 
                               braid_Vector            u,
                               spatial_discretization *spatial_lookup_table)
{
   man->nlx = (spatial_lookup_table[u->spatial_disc_idx]).nlx;
   man->nly = (spatial_lookup_table[u->spatial_disc_idx]).nly;
   man->nx = (spatial_lookup_table[u->spatial_disc_idx]).nx;
   man->ny = (spatial_lookup_table[u->spatial_disc_idx]).ny;
   man->dx = (spatial_lookup_table[u->spatial_disc_idx]).dx;
   man->dy = (spatial_lookup_table[u->spatial_disc_idx]).dy;

   man->ilower[0] = (spatial_lookup_table[u->spatial_disc_idx]).ilower[0];
   man->ilower[1] = (spatial_lookup_table[u->spatial_disc_idx]).ilower[1];
   man->iupper[0] = (spatial_lookup_table[u->spatial_disc_idx]).iupper[0];
   man->iupper[1] = (spatial_lookup_table[u->spatial_disc_idx]).iupper[1];
   man->grid_x    = (spatial_lookup_table[u->spatial_disc_idx]).grid_x_matrix;

   return 0;
}

/* --------------------------------------------------------------------
 * Helper routine to grab the spatial discretization infor about 
 * a braid_Vector u
 * -------------------------------------------------------------------- */
   
int grab_vec_spatial_info(braid_Vector            u,
                          spatial_discretization *spatial_lookup_table,
                          int                     ilower[2],
                          int                     iupper[2],
                          int                    *nlx,
                          int                    *nly,
                          int                    *nx,
                          int                    *ny,
                          double                 *dx,
                          double                 *dy)
{
   (*nlx) = (spatial_lookup_table[u->spatial_disc_idx]).nlx;
   (*nly) = (spatial_lookup_table[u->spatial_disc_idx]).nly;
   (*nx)  = (spatial_lookup_table[u->spatial_disc_idx]).nx;
   (*ny)  = (spatial_lookup_table[u->spatial_disc_idx]).ny;
   (*dx)  = (spatial_lookup_table[u->spatial_disc_idx]).dx;
   (*dy)  = (spatial_lookup_table[u->spatial_disc_idx]).dy;

   ilower[0] = (spatial_lookup_table[u->spatial_disc_idx]).ilower[0];
   ilower[1] = (spatial_lookup_table[u->spatial_disc_idx]).ilower[1];
   iupper[0] = (spatial_lookup_table[u->spatial_disc_idx]).iupper[0];
   iupper[1] = (spatial_lookup_table[u->spatial_disc_idx]).iupper[1];

   return 0;
}


/* --------------------------------------------------------------------
 * Time integrator routine.
 * This routine performs the update
 *   u_i = Phi_i(u_{i-1}) + g_i 
 * Note that the first case corresponds to assuming zero Dirichlet BCs
 * and a zero RHS of the PDE.
 * When Step is called, u is u_{i-1}. At the end of the routine, u is 
 * set to u_i.
 * -------------------------------------------------------------------- */
int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   HYPRE_SStructVector  bstop;
   double accuracy, cfl_value;
   int user_explicit, level;
   int ilower[2], iupper[2], nlx, nly, nx, ny;
   double dx, dy;
   int print_level      = app->print_level;
   double temp_double;
   int temp_int;
   int cfl = 0;
   int iters_taken = -1;
   
   /* This debug output is mostly for regression testing */
   if(print_level == 3)
   {
      braid_GetSpatialAccuracy( status, 1e-2, 1e-9, &temp_double);
      printf("  braid_GetSpatialAccuracy:  %1.2e\n", temp_double);

      braid_StepStatusGetTol(status, &temp_double);
      printf("  braid_StepStatusGetTol:  %1.2e\n", temp_double);

      braid_StepStatusGetIter(status, &temp_int);
      printf("  braid_StepStatusGetIter:  %d\n", temp_int);
      
      temp_int = -1;
      braid_StepStatusGetRNorms(status, &temp_int, &temp_double);
      printf("  braid_StepStatusGetRNorms[-1]:  %1.2e\n", temp_double);

      temp_int = 1;
      braid_StepStatusGetRNorms(status, &temp_int, &temp_double);
      printf("  braid_StepStatusGetRNorms[0]:  %1.2e\n", temp_double);
   }

   /* Grab level */
   braid_StepStatusGetLevel(status, &level);

   /* Grab spatial info about u */
   grab_vec_spatial_info(u, app->spatial_lookup_table, ilower, iupper, 
                         &nlx, &nly, &nx, &ny, &dx, &dy);
   
   /* Grab status of current time step */
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   
   /* Compute the desired spatial solve accuracy */
   if(app->tol_x[0] != app->tol_x[1])
   {
      braid_GetSpatialAccuracy( status, app->tol_x[0], app->tol_x[1], &accuracy);
   }
   else
   {
      accuracy = app->tol_x[0];
   }

   /* 
    * Now, set up the discretization matrix.  Use the XBraid level to index
    * into the matrix lookup table 
    * */

   /* First: check CFL condition, always switch to implicit time stepping if
    * you violate the CFL Note: CFL = K*(dt/dx^2 + dt/dy^2) */
   cfl_value = (app->man->K)*( (tstop-tstart)/((dx)*(dx)) + (tstop-tstart)/((dy)*(dy)) );
   if( cfl_value < 0.5 )
   {
      cfl = 1;
   } 
   
   /* Store information on CFL and spatial coarsening for user output */
   (app->runtime_scoarsen_info)[ (5*level) + 1] = dx;
   (app->runtime_scoarsen_info)[ (5*level) + 2] = dy;
   (app->runtime_scoarsen_info)[ (5*level) + 3] = (tstop-tstart);
   (app->runtime_scoarsen_info)[ (5*level) + 4] = cfl_value;

   /* Update manager relative to this vector */
   update_manager_from_vector(app->man, u, app->spatial_lookup_table);
   
   /* We need to "trick" the user's manager with the new dt */
   app->man->dt = tstop - tstart;

   if( app->dt_A[level] == -1.0 ){
      app->nA++;
      
      /* No matrix for time step tstop-tstart exists. 
       * Add entry to matrix lookup table. */   
      
      app->dt_A[level] = tstop-tstart;
      app->dy_A[level] = dx;
      app->dx_A[level] = dy;

     /* We need to "trick" the user's data structure into mimicking this
      * discretization level */
      app->man->graph = (app->spatial_lookup_table[u->spatial_disc_idx]).graph_matrix;

      /* Set up the implicit or explicit discretization matrix.  If CFL is violated, 
       * automatically use implicit */
     if( app->man->explicit && cfl ){
         setUpExplicitMatrix( app->man );
         app->A[level] = app->man->A;
         
         /* Store that we used explicit on this level */
         (app->runtime_scoarsen_info)[ (5*level) ]    = 1;
      }
      else{
         setUpImplicitMatrix( app->man );
         app->A[level] = app->man->A;
      
         /* Set up the PFMG solver using u->x as dummy vectors. */
         setUpStructSolver( app->man, u->x, u->x );
         app->solver[level] = app->man->solver;

         /* Store that we used implicit on this level */
         (app->runtime_scoarsen_info)[ (5*level) ]    = 0;   
      }
   } 


   /* --------------------------------------------------------------
    * Time integration to next time point: Solve the system Ax = b.
    * 
    * Overwrite the necessary info in the manager, and then call the user's
    * time stepping routine
    * -------------------------------------------------------------- */
   
   /* update manager with a few other important data items */
   app->man->A = app->A[level];
   app->man->solver = app->solver[level];
   
   /* Use level specific max_iter by "tricking" the user's data structure*/
   if( level == 0 )
      app->man->max_iter = app->max_iter_x[0];
   else
      app->man->max_iter = app->max_iter_x[1];
   
   /* Use level specific time stepper by "tricking" the user's data structure*/
   user_explicit = app->man->explicit;
   if( app->man->explicit && cfl )
      app->man->explicit = 1;
   else
      app->man->explicit = 0;

   /* Use level specific tol by "tricking" the user's data structure.  Note that
    * this accuracy value is controlled, in part, by the -tol_x and -tol_xc command
    * line parameters. */
   app->man->tol = accuracy;
   
   /* Take step */
   if (fstop == NULL)
   {
      bstop = NULL;
   }
   else
   {
      bstop = fstop->x;
   }
   take_step(app->man, ustop->x, bstop, u->x, tstart, tstop, &iters_taken);

   /* Go back to the user's original choice of explicit */
   app->man->explicit = user_explicit;
   
   /* Store iterations taken */
   app->runtime_max_iter[level] = max_i( (app->runtime_max_iter[level]),
                                            iters_taken);
   /* Tell XBraid no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}

/* --------------------------------------------------------------------
 * -------------------------------------------------------------------- */
int
my_Residual(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     r,
            braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve u to this time*/
   int ilower[2], iupper[2], nlx, nly, nx, ny, level;
   double dx, dy;
   
   /* Grab level */
   braid_StepStatusGetLevel(status, &level);

   /* Grab spatial info about u */
   grab_vec_spatial_info(ustop, app->spatial_lookup_table, ilower, iupper, 
                         &nlx, &nly, &nx, &ny, &dx, &dy);
 
   /* Grab status of current time step */
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);

   /* Update manager relative to this vector */
   update_manager_from_vector(app->man, ustop, app->spatial_lookup_table);
   
   /* We need to "trick" the user's manager with the new dt */
   app->man->dt = tstop - tstart;

   /* 
    * Now, set up the discretization matrix.  Use the XBraid level to index
    * into the matrix lookup table 
    * */
   if( app->dt_A[level] == -1.0 ){
      app->nA++;
      app->dt_A[level] = tstop-tstart;
      app->dx_A[level] = dx;
      app->dy_A[level] = dy;
      
      /* We need to "trick" the user's data structure into mimicking this
      * discretization level */
      app->man->graph = (app->spatial_lookup_table[ustop->spatial_disc_idx]).graph_matrix;


      setUpImplicitMatrix( app->man );
      app->A[level] = app->man->A;
      
      /* Set up the PFMG solver using r->x as dummy vectors. */
      setUpStructSolver( app->man, r->x, r->x );
      app->solver[level] = app->man->solver;

     /* Store that we used implicit on this level */
     (app->runtime_scoarsen_info)[ (5*level) ]    = 0;   
   } 

   /* Compute residual Ax */
   app->man->A = app->A[level];
   comp_res(app->man, ustop->x, r->x, tstart, tstop);

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
   
   /* Initialize vector and store the spatial index, here it's level 0 */
   my_Vector * u = (my_Vector *) malloc( sizeof(my_Vector) );
   u->spatial_disc_idx = 0; 

   /* Make sure that the manager information on spatial level agrees with u */
   update_manager_from_vector(app->man, u, app->spatial_lookup_table);
   
   if( t == app->man->tstart ){
      /* Sets u_ptr as the initial condition */
      t = 0.0;
   }
   else if (app->use_rand){
      /* Sets u_ptr as uniformly random, for reproducibility use a seed */
      t = -1.0;
   }
   else{
      /* Sets u_ptr as an all zero vector*/
      t = 1.0;
   }

   set_initial_condition(app->man, &(u->x), t);

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
   my_Vector *v;
   double    *values;
   
   /* Grab spatial info about vector */
   int ilower[2], iupper[2], nlx, nly, nx, ny;
   double dx, dy;
   grab_vec_spatial_info(u, app->spatial_lookup_table, ilower, iupper, 
                     &nlx, &nly, &nx, &ny, &dx, &dy);

   v = (my_Vector *) malloc(sizeof(my_Vector));

   /* Create an empty vector object. */
   app->man->grid_x = (app->spatial_lookup_table[u->spatial_disc_idx]).grid_x_matrix;
   initialize_vector(app->man, &(v->x));

   /* Set the values. */
   values = (double *) malloc( (nlx)*(nly)*sizeof(double) );
   HYPRE_SStructVectorGather( u->x );
   HYPRE_SStructVectorGetBoxValues( u->x, 0, ilower, iupper, 0, values );
   HYPRE_SStructVectorSetBoxValues( v->x, 0, ilower, iupper, 0, values );
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
   
   /* Grab spatial info about vector */
   int ilower[2], iupper[2], nlx, nly, nx, ny;
   double dx, dy;
   grab_vec_spatial_info(x, app->spatial_lookup_table, ilower, iupper, 
                     &nlx, &nly, &nx, &ny, &dx, &dy);

   values_x = (double *) malloc( nlx*nly*sizeof(double) );
   values_y = (double *) malloc( nlx*nly*sizeof(double) );

   HYPRE_SStructVectorGather( x->x );
   HYPRE_SStructVectorGetBoxValues( x->x, 0, ilower, iupper, 0, values_x );

   HYPRE_SStructVectorGather( y->x );
   HYPRE_SStructVectorGetBoxValues( y->x, 0, ilower, iupper, 0, values_y );

   for( i = 0; i < nlx*nly; i++ ){
      values_y[i] = alpha*values_x[i] + beta*values_y[i];
   }

   HYPRE_SStructVectorSetBoxValues( y->x, 0, ilower, iupper, 0, values_y );

   free( values_x );
   free( values_y );

   return 0;
}


/* --------------------------------------------------------------------
 * Access the vector -- used for user output
 * -------------------------------------------------------------------- */
int
my_Access(braid_App           app,
          braid_Vector        u,
          braid_AccessStatus  astatus)
{
   MPI_Comm   comm   = MPI_COMM_WORLD;
   double     tstart = (app->man->tstart);
   double     tstop  = (app->man->tstop);
   double     rnorm, disc_err, t;
   int        nt = (app->man->nt);
   int        iter, level, done, calling_fcn;
   int        index, myid;
   static int previous_level = -5;
   static int previous_iter = -5;
   char       filename[255], filename_mesh[255], filename_err[255], filename_sol[255];
  
   /* Retrieve Braid State Information from Status Object */
   MPI_Comm_rank(comm, &myid);
   braid_AccessStatusGetTILD(astatus, &t, &iter, &level, &done);
   braid_AccessStatusGetResidual(astatus, &rnorm);
   braid_AccessStatusGetCallingFunction(astatus, &calling_fcn);
   if( (myid == 0) && ((level != previous_level) || (iter != previous_iter)) )
   {
      previous_level = level;
      previous_iter = iter;
      printf("  my_Access() called, iter= %d, level= %d, calling function= %d\n", iter, level, calling_fcn);
   }

   /* Write a file for each time step that contains a single scalar, the l2
    * norm of the discretization error */
   MPI_Comm_rank(app->comm_x, &myid);
   if( level == 0 )
   {
      /* update manager to correspond to this vector and its time level */
      update_manager_from_vector(app->man, u, app->spatial_lookup_table);
      
      /* Print discretization error to screen for only final time */
      index = ((t - tstart) / ((tstop - tstart)/nt) + 0.1);
      compute_disc_err(app->man, u->x, t, app->e, &disc_err);
      if( (t == app->man->tstop) && myid == 0 ) {
         printf("\n  my_Access():  Braid iter %d,  discr. error at final time:  %1.4e\n", iter, disc_err);
      }
      
      /* Write the norm of the discretization error to a separate file for each time step */
      if( app->man->output_files ){
         sprintf(filename, "%s.iter%03d.time%07d", "drive-02.error_norm", iter, index);
         output_error_file(app->man, t, disc_err, filename); 
      }
   }

  
   /* Write THREE GLVIS visualization files for the final time step:
    * (1) the discretization error (2) the true solution (3) the discrete solution */ 
   if( app->man->output_vis && (level == 0) && (t == app->man->tstop) )
   {
      /* update manager to correspond to this vector and its time level */
      update_manager_from_vector(app->man, u, app->spatial_lookup_table);
      
      /* Create file name and do output */
      sprintf(filename_mesh, "%s.iter%03d", "drive-02_mesh", iter);
      sprintf(filename_err, "%s.iter%03d", "drive-02_err_tstop", iter);
      sprintf(filename_sol, "%s.iter%03d", "drive-02_sol_tstop", iter);
      output_vis(app->man, u->x, t, filename_mesh, filename_err, filename_sol);
   }
   
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
   norm(u->x, norm_ptr);
   return 0;
}


/* --------------------------------------------------------------------
 * Return buffer size for vector object buffer. Vector object contains
 * values at every grid point and thus, the buffer size is the number
 * of grid points.
 * -------------------------------------------------------------------- */
int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
    /* A vector needs to contain 1 extra doubles for the coarsening rule */ 
    *size_ptr = (app->buffer_size);
    return 0;
}


/* --------------------------------------------------------------------
 * Pack a vector object in a buffer.
 * -------------------------------------------------------------------- */
int
my_BufPack(braid_App           app,
           braid_Vector        u,
           void               *buffer,
           braid_BufferStatus  bstatus)
{
   double *dbuffer = buffer;
   
   /* Grab spatial info about vector */
   int ilower[2], iupper[2], nlx, nly, nx, ny;
   double dx, dy;
   grab_vec_spatial_info(u, app->spatial_lookup_table, ilower, iupper, 
                     &nlx, &nly, &nx, &ny, &dx, &dy);

   /* Pack the spatial coarsening rule */ 
   dbuffer[0] = (double) u->spatial_disc_idx; 

   /* Pack the vector */
   HYPRE_SStructVectorGather( u->x );
   HYPRE_SStructVectorGetBoxValues( u->x, 0, ilower,
                                    iupper, 0, &(dbuffer[1]) );

   /* Determine number of bytes actually packed */
   braid_BufferStatusSetSize( bstatus, (nlx*nly + 1)*sizeof(double));
   
   return 0;
}


/* --------------------------------------------------------------------
 * Unpack a vector object from a buffer.
 * -------------------------------------------------------------------- */
int
my_BufUnpack(braid_App           app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus  bstatus)
{
   int ilower[2], iupper[2], nlx, nly, nx, ny;
   double dx, dy;
   double    *dbuffer = buffer;
   my_Vector *u;

   u = (my_Vector *) malloc( sizeof(my_Vector) );

   /* Unpack the spatial coarsening rule */
   u->spatial_disc_idx = (int) dbuffer[0];

   /* Grab spatial info about vector */
   grab_vec_spatial_info(u, app->spatial_lookup_table, ilower, iupper, 
                     &nlx, &nly, &nx, &ny, &dx, &dy);

   /* 
    * Unpack the vector 
    */

   /* Create an empty vector object. */
   app->man->grid_x = (app->spatial_lookup_table[u->spatial_disc_idx]).grid_x_matrix;
   initialize_vector(app->man, &(u->x));
   /* Set the values. */
   HYPRE_SStructVectorSetBoxValues( u->x, 0, ilower, iupper, 0, &(dbuffer[1]) );
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
 * The CFL condition checked is:   
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
   int ncoarsen = max_i(ncoarsen1, ncoarsen2);
   ncoarsen = max_i(ncoarsen, ncoarsen3);

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
 * Return the spatial discretization stored in  app->spatial_lookup_table  
 * that corresponds to the (cdt, fdt) tuple.  Remember that this 
 * table is essentially a database with the tuples (cdt, fdt) as the keys.
 *
 * If spatial coarsening is required, then it is done by generating new SStruct
 * grids and graphs.
 *
 * The grids are numbered by their integer index into the spatial_lookup_table
 * such that grid k is as fine or finer than grid k+1.  
 *
 * If no spatial coarsening is required to satisfy the CFL constraints when 
 * moving from grid k to grid k+1, then these two grids are the same.  
 * Otherwise they are different.  If multiple coarsenings are required to 
 * satisfy the CFl constraint when coarsening grid k, then these intermediate
 * grids are created as well (for use by the spatial interpolation/restriction
 * functions).  If two intermediate grids are needed, then grid k will spawn
 * the creation of grids k+1, k+2 and k+3 in the spatial_lookup_table.  That is
 * spatial_lookup_table[k+1] and spatial_lookup_table[k+2] will hold grid information
 * for intermediate grids only used for spatial interpolation/restriction
 * and spatial_lookup_table[k+3] will the grid information for the next XBraid level.
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
      if( fabs((app->spatial_lookup_table[i]).cdt - cdt)/cdt < 1e-10 )
      {
         /* A rule for this cdt has already been stored */
         if( fabs((app->spatial_lookup_table[i]).fdt - fdt)/fdt < 1e-10 )
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
   ncoarsen = my_ComputeNumCoarsenings(cdt, fdx, fdy, app->man->K, fnx, fny, scoarsenCFL); 
   for(k = 1; k <= max_i(ncoarsen,1); k++)
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
         if(k == max_i(ncoarsen,1) )
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
      cnlx = max_i(ciupper[0] - cilower[0] + 1, 0);
      cnly = max_i(ciupper[1] - cilower[1] + 1, 0);


      /* Note: We start at i=1.  The first entry for the table is a dummy
       *       entry for cloning vectors. 
       * Note: We assume that no more than 3*(app->max_levels) instances
       *       of spatial coarsening will ever occur.
       * */
      for(i = 1; i < 3*max_levels; i++)
      {
         if( ( (app->spatial_lookup_table[i]).cdt == -1.0 ) && ( (app->spatial_lookup_table[i]).fdt == -1.0 ) )
         {
            (app->spatial_lookup_table[i]).cdt = cdt_loc;
            (app->spatial_lookup_table[i]).fdt = fdt;
            (app->spatial_lookup_table[i]).ncoarsen = ncoarsen_loc;
            (app->spatial_lookup_table[i]).dx = cdx;
            (app->spatial_lookup_table[i]).dy = cdy;
            (app->spatial_lookup_table[i]).nx = cnx;
            (app->spatial_lookup_table[i]).ny = cny;
            (app->spatial_lookup_table[i]).nlx = cnlx;
            (app->spatial_lookup_table[i]).nly = cnly;
            (app->spatial_lookup_table[i]).ilower[0] = cilower[0];
            (app->spatial_lookup_table[i]).ilower[1] = cilower[1];
            (app->spatial_lookup_table[i]).iupper[0] = ciupper[0];
            (app->spatial_lookup_table[i]).iupper[1] = ciupper[1];
            (app->spatial_lookup_table[i]).fspatial_disc_idx = fspatial_disc_idx; 
            
            /* We also have to set up a new 2D grid and graph for the new COARSE level*/
            setUp2Dgrid( app->comm_x, &((app->spatial_lookup_table[i]).grid_x_matrix), app->man->dim_x,  
                         cilower, ciupper, app->man->vartype, 1 );
            setUpGraph( app->comm_x, &((app->spatial_lookup_table[i]).graph_matrix), 
                       (app->spatial_lookup_table[i]).grid_x_matrix, app->man->object_type, 
                        app->man->stencil );
            /* We temporarily let the vector grid and graph be the matrix grid and graph. Once we figure
             * out the correct ghost layer, we will replace these two values */
            (app->spatial_lookup_table[i]).grid_x_vector = (app->spatial_lookup_table[i]).grid_x_matrix;
            (app->spatial_lookup_table[i]).graph_vector = (app->spatial_lookup_table[i]).graph_matrix;
            
            /* We also have to set up a new 2D grid and graph for the fine level vectors that has the
             * correct number of ghost layers */
            setUp2Dgrid( app->comm_x, &((app->spatial_lookup_table[fspatial_disc_idx]).grid_x_vector), 
                         app->man->dim_x, (app->spatial_lookup_table[fspatial_disc_idx]).ilower, 
                         (app->spatial_lookup_table[fspatial_disc_idx]).iupper, app->man->vartype, 2);
            setUpGraph(app->comm_x, &((app->spatial_lookup_table[fspatial_disc_idx]).graph_vector), 
                       (app->spatial_lookup_table[fspatial_disc_idx]).grid_x_vector, 
                       app->man->object_type, app->man->stencil );

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
 * Retrieve a spatial discretization stored in app->spatial_lookup_table  
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
      if( fabs((app->spatial_lookup_table[i]).cdt - cdt)/cdt < 1e-10 )
      {
         /* A rule for this cdt has already been stored */
         if( fabs((app->spatial_lookup_table[i]).fdt - fdt)/fdt < 1e-10 )
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
   braid_CoarsenRefStatusGetT(status, &tstart); */
   
   /* Determine current local spatial grid size.  Note that subtracting one from a 
    * spatial discretization index will give you the index for the next finer grid. */
   spatial_disc_idx  = cu->spatial_disc_idx; 
   fspatial_disc_idx = spatial_disc_idx - 1;
   cnlx              = (app->spatial_lookup_table[spatial_disc_idx]).nlx;
   cnly              = (app->spatial_lookup_table[spatial_disc_idx]).nly;
   cilower[0]        = (app->spatial_lookup_table[spatial_disc_idx]).ilower[0];
   cilower[1]        = (app->spatial_lookup_table[spatial_disc_idx]).ilower[1];
   ciupper[0]        = (app->spatial_lookup_table[spatial_disc_idx]).iupper[0];
   ciupper[1]        = (app->spatial_lookup_table[spatial_disc_idx]).iupper[1];

   /* Find this processor's part of the fine vector.  If there is no
    * coarsening, then fu will be the same size as cu. */    
   fnlx              = (app->spatial_lookup_table[fspatial_disc_idx]).nlx;
   fnly              = (app->spatial_lookup_table[fspatial_disc_idx]).nly;
   filower[0]        = (app->spatial_lookup_table[fspatial_disc_idx]).ilower[0];
   filower[1]        = (app->spatial_lookup_table[fspatial_disc_idx]).ilower[1];
   fiupper[0]        = (app->spatial_lookup_table[fspatial_disc_idx]).iupper[0];
   fiupper[1]        = (app->spatial_lookup_table[fspatial_disc_idx]).iupper[1];


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
   app->man->grid_x = (app->spatial_lookup_table[fspatial_disc_idx]).grid_x_vector;
   initialize_vector(app->man, &(fu->x));
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
        * $$ srun -N 1 -n 1 -p pdebug ./drive-02 -pgrid 1 1 1 -nt 256 -mi 2 -ml 15 -nx 9 9 -scoarsen 2 
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
   app->man->grid_x = (app->spatial_lookup_table[fspatial_disc_idx]).grid_x_matrix;
   initialize_vector(app->man, &(fu->x));
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
   braid_CoarsenRefStatusGetT(status, &tstart);
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
   ncoarsen          = (app->spatial_lookup_table[spatial_disc_idx]).ncoarsen;
   
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
   int        fnlx = (app->spatial_lookup_table[fspatial_disc_idx]).nlx;
   int        fnly = (app->spatial_lookup_table[fspatial_disc_idx]).nly;
   int        fnx  = (app->spatial_lookup_table[fspatial_disc_idx]).nx;
   int        fny  = (app->spatial_lookup_table[fspatial_disc_idx]).ny;

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

   braid_CoarsenRefStatusGetT(status, &tstart);
   cu = (my_Vector *) malloc(sizeof(my_Vector));
   
   filower[0]    = (app->spatial_lookup_table[fspatial_disc_idx]).ilower[0];
   filower[1]    = (app->spatial_lookup_table[fspatial_disc_idx]).ilower[1];
   fiupper[0]    = (app->spatial_lookup_table[fspatial_disc_idx]).iupper[0];
   fiupper[1]    = (app->spatial_lookup_table[fspatial_disc_idx]).iupper[1];
   
   /* Get the next coarser spatial discretization */ 
   cnlx     = (app->spatial_lookup_table[fspatial_disc_idx+1]).nlx;
   cnly     = (app->spatial_lookup_table[fspatial_disc_idx+1]).nly;
   cilower[0] = (app->spatial_lookup_table[fspatial_disc_idx+1]).ilower[0];
   cilower[1] = (app->spatial_lookup_table[fspatial_disc_idx+1]).ilower[1];
   ciupper[0] = (app->spatial_lookup_table[fspatial_disc_idx+1]).iupper[0];
   ciupper[1] = (app->spatial_lookup_table[fspatial_disc_idx+1]).iupper[1];

   /* Create an empty vector object. */
   app->man->grid_x = (app->spatial_lookup_table[fspatial_disc_idx+1]).grid_x_matrix;
   initialize_vector(app->man, &(cu->x));

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
       * $$ srun -N 1 -n 1 -p pdebug ./drive-02 -pgrid 1 1 1 -nt 256 -mi 2 -ml 15 -nx 9 9 -scoarsen 2 
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
   int        fnx = (app->spatial_lookup_table[fspatial_disc_idx]).nx;
   int        fny = (app->spatial_lookup_table[fspatial_disc_idx]).ny;
   double     fdy = (app->spatial_lookup_table[fspatial_disc_idx]).dy;
   double     fdx = (app->spatial_lookup_table[fspatial_disc_idx]).dx;

   int        k, ncoarsen, level;
   double     cdt, fdt, scoarsenCFL;
   double     tstart, f_tstop, f_tprior, c_tstop, c_tprior;
   
   /* Get Coarse and fine time step sizes */
   braid_CoarsenRefStatusGetT(status, &tstart);
   braid_CoarsenRefStatusGetCTstop(status, &c_tstop);
   braid_CoarsenRefStatusGetCTprior(status, &c_tprior);
   braid_CoarsenRefStatusGetFTstop(status, &f_tstop);
   braid_CoarsenRefStatusGetFTprior(status, &f_tprior);
   cdt = c_tstop - tstart;
   fdt = f_tstop - tstart;
   
   filower[0]    = (app->spatial_lookup_table[fspatial_disc_idx]).ilower[0];
   filower[1]    = (app->spatial_lookup_table[fspatial_disc_idx]).ilower[1];
   fiupper[0]    = (app->spatial_lookup_table[fspatial_disc_idx]).iupper[0];
   fiupper[1]    = (app->spatial_lookup_table[fspatial_disc_idx]).iupper[1];
   
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
   
   /* Generate the next spatial discretization, which is stored in app->spatial_lookup_table[i]
    * This could be the same as the fine spatial discretization */ 
   braid_CoarsenRefStatusGetLevel(status, &level);
   scoarsenCFL = app->scoarsenCFL[min_i(app->num_scoarsenCFL-1, level)];
   get_coarse_spatial_disc( app, cdt, fdt, fdx, fdy, scoarsenCFL, fnx, fny, 
                            fspatial_disc_idx, filower, fiupper, &spatial_disc_idx);
   ncoarsen = (app->spatial_lookup_table[spatial_disc_idx]).ncoarsen;
         
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
   int        fnx = (app->spatial_lookup_table[fspatial_disc_idx]).nx;
   int        fny = (app->spatial_lookup_table[fspatial_disc_idx]).ny;
   int        fnlx = (app->spatial_lookup_table[fspatial_disc_idx]).nlx;
   int        fnly = (app->spatial_lookup_table[fspatial_disc_idx]).nly;
   double     fdy = (app->spatial_lookup_table[fspatial_disc_idx]).dy;
   double     fdx = (app->spatial_lookup_table[fspatial_disc_idx]).dx;

   int        cilower[2], ciupper[2];
   double     coarsen_factor;
   int        cnlx, cnly, ncoarsen, coarsen_factor_int, level;
   double     cdt, fdt, scoarsenCFL;
   double     tstart, f_tstop, f_tprior, c_tstop, c_tprior;

   /* Get Coarse and fine time step sizes */
   braid_CoarsenRefStatusGetT(status, &tstart);
   braid_CoarsenRefStatusGetCTstop(status, &c_tstop);
   braid_CoarsenRefStatusGetCTprior(status, &c_tprior);
   braid_CoarsenRefStatusGetFTstop(status, &f_tstop);
   braid_CoarsenRefStatusGetFTprior(status, &f_tprior);
   cdt = c_tstop - tstart;
   fdt = f_tstop - tstart;
   
   cu = (my_Vector *) malloc(sizeof(my_Vector));
   
   filower[0]    = (app->spatial_lookup_table[fspatial_disc_idx]).ilower[0];
   filower[1]    = (app->spatial_lookup_table[fspatial_disc_idx]).ilower[1];
   fiupper[0]    = (app->spatial_lookup_table[fspatial_disc_idx]).iupper[0];
   fiupper[1]    = (app->spatial_lookup_table[fspatial_disc_idx]).iupper[1];
   
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
   
   /* Generate the next spatial discretization, which is stored in app->spatial_lookup_table[i]
    * This could be the same as the fine spatial discretization */ 
   braid_CoarsenRefStatusGetLevel(status, &level);
   scoarsenCFL = app->scoarsenCFL[min_i(app->num_scoarsenCFL-1, level)];
   get_coarse_spatial_disc( app, cdt, fdt, fdx, fdy, scoarsenCFL, fnx, fny, 
                            fspatial_disc_idx, filower, fiupper, &spatial_disc_idx);
   ncoarsen = (app->spatial_lookup_table[spatial_disc_idx]).ncoarsen;
   cnlx     = (app->spatial_lookup_table[spatial_disc_idx]).nlx;
   cnly     = (app->spatial_lookup_table[spatial_disc_idx]).nly;
   cilower[0] = (app->spatial_lookup_table[spatial_disc_idx]).ilower[0];
   cilower[1] = (app->spatial_lookup_table[spatial_disc_idx]).ilower[1];
   ciupper[0] = (app->spatial_lookup_table[spatial_disc_idx]).iupper[0];
   ciupper[1] = (app->spatial_lookup_table[spatial_disc_idx]).iupper[1];
   coarsen_factor = pow(2.0, ncoarsen);
   coarsen_factor_int = (int) coarsen_factor;

   /* Create an empty vector object. */
   app->man->grid_x = (app->spatial_lookup_table[spatial_disc_idx]).grid_x_matrix;
   initialize_vector(app->man, &(cu->x));

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
   /* Declare variables -- variables explained when they are set below */
   int print_usage                       = 0;
   int object_type                       = HYPRE_STRUCT;
   int *runtime_max_iter_global    = NULL;
   double *runtime_scoarsen_info_global  = NULL;
   
   int i, arg_index, myid, num_procs;
   MPI_Comm comm, comm_x, comm_t;
   int ndim, nx, ny, nlx, nly, nt, forcing, ilower[2], iupper[2];
   double K, tstart, tstop, dx, dy, dt, cfl;
   int px, py, pt, pi, pj;
   int output_files, explicit, output_vis;

   /* Declare Braid variables -- variables explained when they are set below */
   braid_Core    core;
   my_App       *app;
   double tol_x[2], tol, *scoarsenCFL;
   double mystarttime, myendtime, mytime, maxtime;
   int run_wrapper_tests, correct, fspatial_disc_idx, max_iter_x[2];
   int print_level, access_level, nA_max, max_levels, min_coarse, skip;
   int nrelax, nrelax0, finalFC, cfactor, cfactor0, max_iter, storage, res, new_res;
   int fmg, tnorm, nfmg_Vcyc, scoarsen, num_scoarsenCFL, use_rand;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);

   /* Default parameters from user code */
   comm                = MPI_COMM_WORLD;
   ndim                = 2;       /* Two dimensional problem */
   forcing             = 0;       /* Boolean, if 1 use a nonzero forcing term, if 0 use a 0 forcing term */
   K                   = 1.0;     /* Diffusion coefficient */
   nx                  = 17;      /* number of points in the x-dim */
   ny                  = 17;      /* number of points in the y-dim */
   nlx                 = 17;      /* number of point ~local~ to this processor in the x-dim */
   nly                 = 17;      /* number of point ~local~ to this processor in the y-dim */
   tstart              = 0.0;     /* global start time */
   nt                  = 32;      /* number of time steps */
   cfl                 = 0.30;    /* CFL ratio of dt/dx^2 to use */
   px                  = 1;       /* my processor number in the x-direction, px*py=num procs in space */
   py                  = 1;       /* my processor number in the y-direction, px*py=num procs in space */
   explicit            = 0;       /* Boolean, if 1 use explicit forward Euler, else use backward Euler */
   output_files        = 0;       /* Boolean, if 1 output the norm of the discretization error to a file for each time step */
   output_vis          = 0;       /* Boolean, if 1 output GLVIS files of error, solution and true solution */
   skip                = 1;       /* Boolean, if 1 do no work on the first down cycle, if 0 do work on the first down cycle */
   
   /* Default XBraid parameters */
   max_levels          = 15;  /* Must be two, in order for coarsen and refine wrapper tests to run */
   min_coarse          = 3;  
   nrelax              = 1;
   nrelax0             = -1;
   finalFC             = 0;
   tol                 = 1.0e-09;
   tnorm               = 2;
   cfactor             = 2;
   cfactor0            = -1;
   max_iter            = 100;
   fmg                 = 0;
   nfmg_Vcyc           = 1;
   res                 = 0;
   new_res             = 0;
   storage             = -1;
   use_rand            = 1;
   scoarsen            = 0;
   num_scoarsenCFL     = 1;
   scoarsenCFL         = (double*) malloc( 1*sizeof(double) );
   scoarsenCFL[0]      = 0.5;
   pt                  = 1;  
   max_iter_x[0]       = 50;
   max_iter_x[1]       = 50;
   tol_x[0]            = 1.0e-09;
   tol_x[1]            = 1.0e-09;
   print_level         = 2;
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
      else if( strcmp(argv[arg_index], "-skip") == 0 ){
          arg_index++;
          skip = atoi(argv[arg_index++]);
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
      else if ( strcmp(argv[arg_index], "-finalFC") == 0 ){
         arg_index++;
         finalFC = 1;
      }
      else if ( strcmp(argv[arg_index], "-storage") == 0 ){
         arg_index++;
         storage = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-new_res") == 0 ){
         arg_index++;
         new_res = 1;
      }
      else if ( strcmp(argv[arg_index], "-res") == 0 ){
         arg_index++;
         res = 1;
      }
      else if ( strcmp(argv[arg_index], "-use_rand") == 0 ){
         arg_index++;
         use_rand = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-forcing") == 0 ){
         arg_index++;
         forcing = 1;
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
      else if( strcmp(argv[arg_index], "-pfmg_mi") == 0 ){
         arg_index++;
         max_iter_x[0] = atoi(argv[arg_index++]);
         max_iter_x[1] = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-pfmg_tolx") == 0 ){
          arg_index++;
          tol_x[0] = atof(argv[arg_index++]);
          tol_x[1] = atof(argv[arg_index++]);
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
      printf("  -run_wrapper_tests                 : Only run the Braid wrapper tests\n");
      printf("                                       (do not combine with temporal parallelism)\n");
      printf("  -pgrid  <px py pt>                 : processors in each dimension (default: 1 1 1)\n");
      printf("  -nx  <nx ny>                       : 2D spatial problem size of form 2^k+1, 2^k+1 (default: 17 17)\n");
      printf("  -nt  <n>                           : number of time steps (default: 32)\n"); 
      printf("  -cfl <cfl>                         : CFL number to run (default: 0.30)\n"); 
      printf("                                       Note: CFL = K*(dt/dx^2 + dt/dy^2) is used to define dt and t-final\n");
      printf("                                       Note: CFL < 0.5 for explicit forward Euler\n");
      printf("  -ml  <max_levels>                  : set max number of time levels (default: 15)\n");
      printf("  -skip <skip>                       : boolean, whether to skip all work on first down cycle (default: 1\n");
      printf("  -mc  <min_coarse>                  : set min possible coarse level size (default: 3)\n");
      printf("  -nu  <nrelax>                      : set num F-C relaxations (default: 1)\n");
      printf("  -nu0 <nrelax>                      : set num F-C relaxations on level 0\n");
      printf("  -finalFC                           : set to do a final FC-relaxation after the cycle ends\n");
      printf("  -tol <tol>                         : set stopping tolerance (default: 1e-09)\n");
      printf("  -tnorm <tnorm>                     : set temporal norm \n");
      printf("                                       1 - One-norm \n");
      printf("                                       2 - Two-norm (default) \n");
      printf("                                       3 - Infinity-norm \n");
      printf("  -cf  <cfactor>                     : set coarsening factor (default: 2)\n");   
      printf("  -cf0  <cfactor>                    : set coarsening factor for level 0 \n");
      printf("  -mi  <max_iter>                    : set max iterations (default: 100)\n");
      printf("  -pfmg_mi <mi_fine mi_coarse>       : max number of PFMG iters for fine and coarse levels (default: 50 50)\n"); 
      printf("  -pfmg_tolx <loose_tol tight_tol>   : loose and tight PFMG stopping tol (default: 1e-09 1e-09)\n"); 
      printf("  -fmg <nfmg_Vcyc>                   : use FMG cycling, nfmg_Vcyc V-cycles at each fmg level\n");
      printf("  -res                               : use my residual\n");
      printf("  -new_res                           : use user residual routine to compute full residual each iteration\n");
      printf("                                       on all grid points for stopping criterion.\n");
      printf("  -storage <level>                   : full storage on levels >= level\n");
      printf("  -forcing                           : consider non-zero RHS b(x,y,t) = -sin(x)*sin(y)*(sin(t)-2*cos(t))\n");
      printf("  -use_rand <bool>                   : if nonzero, then use a uniformly random value to initialize each\n");
      printf("                                       time step for t>0.  if zero, then use a zero initial guess.\n");
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
      printf("                                     Default is n=1 and cfls=0.5, for stable explicit time-stepping.\n");
      printf("                                     This option must be used with scoarsen > 0.\n");
      printf("                                     \n");
      printf("                                     \n");
      printf(" Output related parameters\n");
      printf(" -------------------------\n");
      printf("  -print_level <l>                : sets the print_level (default: 1) \n");
      printf("                                    0 - no output to standard out \n");
      printf("                                    2 - Basic convergence information and hierarchy statistics\n");
      printf("                                    2 - Debug level output \n");
      printf("  -access_level <l>               : sets the access_level (default: 1) \n");
      printf("                                    0 - never call access \n");
      printf("                                    1 - call access only after completion \n");
      printf("                                    2 - call access at the end of every cycle \n");
      printf("                                    3 - call access every iteration and level\n");
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
   if( (px*py*pt) != num_procs)
   {
       if( myid == 0 )
           printf("Error: px x py x pt does not equal the number of processors!\n");
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

   /* Determine position (pi, pj)  in the 2D processor grid, 
    * 0 <= pi < px,   0 <= pj < py */
   MPI_Comm_rank( comm_x, &myid );
   pi = myid % px;
   pj = ( (myid - pi)/px ) % py;

   /* Define the 2D block of the global grid owned by this processor, that is
    * (ilower[0], iupper[0]) x (ilower[1], iupper[1])
    * defines the piece of the global grid owned by this processor. */ 
   GetDistribution_x( nx, px, pi, &ilower[0], &iupper[0] );
   GetDistribution_x( ny, py, pj, &ilower[1], &iupper[1] );

   /* Determine local problem size. */
   nlx = iupper[0] - ilower[0] + 1;
   nly = iupper[1] - ilower[1] + 1;

   /* Compute grid spacing. */
   dx = PI / (nx - 1);
   dy = PI / (ny - 1);

   /* Set time-step size, noting that the CFL number definition 
    *     K*(dt/dx^2 + dt/dy^2) = CFL
    * implies that dt is equal to 
    *     dt = ( CFL dx^2 dy^2) / ( K(dx^2 + dy^2)) */
   dt = (cfl*(dx*dx)*(dy*dy)) / (K*((dx*dx)+(dy*dy)));
   /* Now using dt, compute the final time, tstop value */
   tstop =  tstart + nt*dt;

   /* -----------------------------------------------------------------
    * Set up App structure.
    * ----------------------------------------------------------------- */
   app = (my_App *) malloc(sizeof(my_App));

   /* Set up the serial simulation manager values */
   (app->man)                  = (simulation_manager *) malloc(sizeof(simulation_manager));
   (app->man->comm)            = comm_x;
   (app->man->forcing)         = forcing;
   (app->man->K)               = K;
   (app->man->dim_x)           = ndim;
   (app->man->nlx)             = nlx;
   (app->man->nly)             = nly;
   (app->man->nx)              = nx;
   (app->man->ny)              = ny;
   (app->man->tstart)          = tstart;
   (app->man->tstop)           = tstop;
   (app->man->nt)              = nt;
   (app->man->dx)              = dx;
   (app->man->dy)              = dy;
   (app->man->dt)              = dt;
   (app->man->px)              = px;
   (app->man->py)              = py;
   (app->man->pi)              = pi;
   (app->man->pj)              = pj;
   (app->man->ilower[0])       = ilower[0];
   (app->man->ilower[1])       = ilower[1];
   (app->man->iupper[0])       = iupper[0];
   (app->man->iupper[1])       = iupper[1];
   (app->man->object_type)     = object_type;
   (app->man->max_iter)        = max_iter;
   (app->man->tol)             = tol;
   (app->man->explicit)        = explicit;
   (app->man->output_vis)      = output_vis;
   (app->man->output_files)    = output_files;
   
   /* Set up the variable type, grid, stencil and matrix graph. */
   app->man->vartype           = HYPRE_SSTRUCT_VARIABLE_CELL;
   setUp2Dgrid( comm_x, &(app->man->grid_x), app->man->dim_x,
                app->man->ilower, app->man->iupper, app->man->vartype, 1 );
   set5ptStencil( &(app->man->stencil), app->man->dim_x );
   setUpGraph( comm_x, &(app->man->graph), app->man->grid_x, 
               app->man->object_type, app->man->stencil );

   /* 
    * Set up the XBraid app structure 
    */
   (app->comm_t)          = comm_t;
   (app->comm_x)          = comm_x;
   (app->max_levels)      = max_levels;
   (app->buffer_size)     = (1 + nlx*nly)*sizeof(double);


   /* Set the maximum number of PFMG iterations for expensive (index 0)
    * and cheap (index 1) solves. */
   (app->max_iter_x)      = (int*) malloc( 2*sizeof(int) );
   (app->max_iter_x[0])   = max_iter_x[0];
   (app->max_iter_x[1])   = max_iter_x[1];

   /* Initialize the storage structure for recording spatial coarsening information */ 
   app->runtime_scoarsen_info = (double*) malloc( 5*(app->max_levels)*sizeof(double) );
   for( i = 0; i < 5*(app->max_levels); i++)
   {
      app->runtime_scoarsen_info[i] = -1.0;
   }
   
   /* Store whether we want a random or constant initial guess */
   app->use_rand = use_rand;
   
   /* Store print-level (controls debug output) */
   app->print_level = print_level;


   /* Store CFl that controls spatial coarsening */
   app->num_scoarsenCFL = num_scoarsenCFL;
   app->scoarsenCFL = scoarsenCFL;
   
   /* Allocate error vector */
   initialize_vector(app->man, &(app->e));

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
   app->runtime_max_iter = (int*) calloc( (app->max_levels),  sizeof(int) );
   for( i = 0; i < app->max_levels; i++ )
      app->runtime_max_iter[i] = 0;

   /* Setup the lookup table that records how grids are coarsened (refined)
    * spatially.  
    * Note: that the first entry for the table is a dummy entry 
    *       for cloning vectors.  
    * Note: We assume that no more than 3*(app->max_levels) instances
    *       of spatial coarsening will ever occur.*/
   (app->spatial_lookup_table) = (spatial_discretization*) malloc( 3*max_levels*sizeof(spatial_discretization) );
   for( i = 1; i < 3*(app->max_levels); i++ )
   {
      app->spatial_lookup_table[i].fdt = -1.0;
      app->spatial_lookup_table[i].cdt = -1.0;
   }
   (app->spatial_lookup_table[0]).grid_x_matrix = (app->man->grid_x);
   (app->spatial_lookup_table[0]).grid_x_vector = (app->man->grid_x);
   (app->spatial_lookup_table[0]).graph_matrix = (app->man->graph);
   (app->spatial_lookup_table[0]).graph_vector = (app->man->graph);
   (app->spatial_lookup_table[0]).fdt = (app->man->dt);
   (app->spatial_lookup_table[0]).cdt = (app->man->dt);
   (app->spatial_lookup_table[0]).ncoarsen = 0;
   (app->spatial_lookup_table[0]).dx = (app->man->dx);
   (app->spatial_lookup_table[0]).dy = (app->man->dy);
   (app->spatial_lookup_table[0]).nx = nx;
   (app->spatial_lookup_table[0]).ny = ny;
   (app->spatial_lookup_table[0]).nlx = (app->man->nlx);
   (app->spatial_lookup_table[0]).nly = (app->man->nly);
   (app->spatial_lookup_table[0]).ilower[0] = (app->man->ilower[0]);
   (app->spatial_lookup_table[0]).ilower[1] = (app->man->ilower[1]);
   (app->spatial_lookup_table[0]).iupper[0] = (app->man->iupper[0]);
   (app->spatial_lookup_table[0]).iupper[1] = (app->man->iupper[1]);
   (app->spatial_lookup_table[0]).fspatial_disc_idx = 0;
   
   
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
      if( max_levels == 2){
         correct = braid_TestCoarsenRefine(app, comm_x, stdout, 0.0, dt, 2*dt, my_Init,
                             my_Access, my_Free, my_Clone, my_Sum, my_SpatialNorm, my_CoarsenInjection, 
                             my_Refine);

         correct = braid_TestCoarsenRefine(app, comm_x, stdout, 0.0, dt, 2*dt, my_Init,
                             my_Access, my_Free, my_Clone, my_Sum, my_SpatialNorm, my_CoarsenBilinear, 
                             my_Refine);
      }
      else {
        printf("\nCoarsen and refine wrapper tests not run.  They require '-ml 2'. \n\n");
      }


      if(correct == 0) {
        printf("Failed: at least one of the tests failed\n");
      }
      else {
        printf("Passed: all tests passed\n");
      }

   }
   else
   {
      /* Run a Braid simulation */
      /* Start timer. */
      mystarttime = MPI_Wtime();

      braid_Init(comm, comm_t, tstart, tstop, nt, app, my_Step, my_Init,
            my_Clone, my_Free, my_Sum, my_SpatialNorm, my_Access, my_BufSize,
            my_BufPack, my_BufUnpack, &core);

      app->tol_x[0] = tol_x[0];
      app->tol_x[1] = tol_x[1];

      braid_SetSkip( core, skip );
      braid_SetMaxLevels( core, max_levels );
      braid_SetMinCoarse( core, min_coarse );

      braid_SetPrintLevel( core, print_level);
      braid_SetAccessLevel( core, access_level);

      braid_SetNRelax(core, -1, nrelax);
      if (nrelax0 > -1)
      {
         braid_SetNRelax(core,  0, nrelax0);
      }
      if (finalFC)
      {
         braid_SetFinalFCRelax(core);
      }

      /*braid_SetRelTol(core, tol);*/
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
      if (res)
      {
         
         if(app->man->explicit) {
            printf("\nCannot mix -res and -expl.  Option -res designed to make the implicit time stepping cheaper. \nIgnoring -res\n\n");
         }
         else{
            braid_SetResidual(core, my_Residual);
         }
      }
      if (new_res) 
      {
         if(app->man->explicit) {
            printf("\nCannot mix -new_res and -expl.  This uses the residual function designed to make the implicit\ntime stepping cheaper. \nIgnoring -new_res\n\n");
         }
         else{
            braid_SetFullRNormRes(core, my_Residual);        
         }
      
      }
      if (storage >= -2)
      {
         braid_SetStorage(core, storage);
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
      
      /* Debug level printing for regression testing */
      if(print_level == 3)
      {
         int nrequest, num_iter, num_iter2;
         braid_GetNumIter(core, &num_iter);
         double *rnorms_array = (double*) malloc( (num_iter+1)*sizeof(double) ); 
         
         nrequest = 1;
         braid_GetRNorms(core, &nrequest, rnorms_array);
         printf("  braid_GetRNorms[0] = %1.2e\n", rnorms_array[0]);
         nrequest = -1;
         braid_GetRNorms(core, &nrequest, rnorms_array);
         printf("  braid_GetRNorms[-1] = %1.2e\n", rnorms_array[0]);
         
         braid_GetRNorms(core, &num_iter, rnorms_array);
         for(i=0; i < num_iter; i++) {
            printf("  braid_GetRNorms[%d] = %1.2e\n", i, rnorms_array[i]);
         }
         
         num_iter2 = -num_iter - 1;
         braid_GetRNorms(core, &num_iter2, rnorms_array);
         for(i=0; i <= num_iter; i++) {
            printf("  braid_GetRNorms[%d] = %1.2e\n", i, rnorms_array[i]);
         }

         printf("\n");
         free(rnorms_array);  
      }

      /* Stop timer. */
      myendtime = MPI_Wtime();
      mytime    = myendtime - mystarttime;

      /* Compute maximum time */
      MPI_Reduce( &mytime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, comm );

      /* Determine some solver information about topics like spatial coarsening
       * and the maximum number of iterations in spatial solves on each time level
       * (if implicit used) */
      MPI_Allreduce( &(app->nA), &nA_max, 1, MPI_INT, MPI_MAX, comm ); 
      runtime_max_iter_global = (int*) malloc( nA_max*sizeof(int) );
      runtime_scoarsen_info_global = (double*) malloc( 5*nA_max*sizeof(double) ); 
      for( i = 0; i < nA_max; i++ ){
         /* Grab max num interations information */
         MPI_Allreduce( &(app->runtime_max_iter[i]), 
                        &runtime_max_iter_global[i], 1, MPI_INT, MPI_MAX, comm );

         /* Grab spatial coarsening information */
         MPI_Allreduce( &(app->runtime_scoarsen_info)[i*5],
               &runtime_scoarsen_info_global[i*5], 1, MPI_DOUBLE, MPI_MAX, comm ); 
         MPI_Allreduce( &(app->runtime_scoarsen_info)[i*5+1],
               &runtime_scoarsen_info_global[i*5+1], 1, MPI_DOUBLE, MPI_MAX, comm ); 
         MPI_Allreduce( &(app->runtime_scoarsen_info)[i*5+2],
               &runtime_scoarsen_info_global[i*5+2], 1, MPI_DOUBLE, MPI_MAX, comm ); 
         MPI_Allreduce( &(app->runtime_scoarsen_info)[i*5+3],
               &runtime_scoarsen_info_global[i*5+3], 1, MPI_DOUBLE, MPI_MAX, comm ); 
         MPI_Allreduce( &(app->runtime_scoarsen_info)[i*5+4],
               &runtime_scoarsen_info_global[i*5+4], 1, MPI_DOUBLE, MPI_MAX, comm ); 
      
      }

      if( myid == 0 )
      {
         printf( "  runtime: %.5lfs\n", maxtime );

         printf("\n\n-----------------------------------------------------------------\n"); 
         printf("-----------------------------------------------------------------\n\n"); 
         printf(" Implicit time stepping solve parameters\n\n");
         printf("   Fine-level loose stopping tol  :  %1.2e    (while ||r|| is large)\n", tol_x[0]);
         printf("   Fine-level tight stopping tol  :  %1.2e    (while ||r|| is small)\n", tol_x[1]);
         printf("   Coarse-level stopping tol      :  %1.2e    (for all ||r||) \n", tol_x[0]);
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
            if( runtime_scoarsen_info_global[i*5] == 1)
            {
               printf(" %2d   |  expl       %1.2e    %1.2e    %1.2e    %1.2e\n", i,
                  runtime_scoarsen_info_global[i*5+1], runtime_scoarsen_info_global[i*5+2],
                  runtime_scoarsen_info_global[i*5+3], runtime_scoarsen_info_global[i*5+4]);
            }
            else
            {
               printf(" %2d   |  impl, %2d   %1.2e    %1.2e    %1.2e    %1.2e\n", 
                  i, runtime_max_iter_global[i],
                  runtime_scoarsen_info_global[i*5+1], runtime_scoarsen_info_global[i*5+2],
                  runtime_scoarsen_info_global[i*5+3], runtime_scoarsen_info_global[i*5+4]);
            }
         }

         printf( "\n" );
      }

      braid_Destroy(core);
   }
   
   /* Free memory */
   if( runtime_max_iter_global != NULL)
   {
      free( runtime_max_iter_global );
   }
   HYPRE_SStructVectorDestroy( app->e );
   HYPRE_SStructGridDestroy( app->man->grid_x );
   HYPRE_SStructStencilDestroy( app->man->stencil );
   HYPRE_SStructGraphDestroy( app->man->graph );
   free( app->dt_A );
   free( app->dx_A );
   free( app->dy_A );
   free( app->scoarsenCFL);
   for( i = 0; i < app->nA; i++ )
   {
      HYPRE_SStructMatrixDestroy( app->A[i] );
      if( runtime_scoarsen_info_global != NULL)
      {
         if( runtime_scoarsen_info_global[i*5] == 0 )
         {
            /* If an implicit solver was used */
            HYPRE_StructPFMGDestroy( app->solver[i] );
         }
      }
   }
   free( runtime_scoarsen_info_global );
   free( app->A );
   free( app->solver );
   free( app->runtime_max_iter );
   free( app->runtime_scoarsen_info );
   free( app->max_iter_x );
   
   /* Destroy all coarse grid_x and graphs
    * The app->man->grid_x and app->man->graph are destroyed above, and correspond to 
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
      if( ( (app->spatial_lookup_table[i]).cdt != -1.0 ) && ( (app->spatial_lookup_table[i]).fdt != -1.0 ) )
      {
         HYPRE_SStructGridDestroy( (app->spatial_lookup_table[i]).grid_x_matrix );
         HYPRE_SStructGraphDestroy( (app->spatial_lookup_table[i]).graph_matrix );
         fspatial_disc_idx = (app->spatial_lookup_table[i]).fspatial_disc_idx;

         /* Note that coarse spatial grid creation also touches its relative fine grid with a
          * _vector creation of grid and graph with large ghost layers.  Destroy that */
         if(fspatial_disc_idx != -1)
         {
            HYPRE_SStructGridDestroy( (app->spatial_lookup_table[fspatial_disc_idx]).grid_x_vector );
            HYPRE_SStructGraphDestroy( (app->spatial_lookup_table[fspatial_disc_idx]).graph_vector );
         }
      }
   }
   /* Note that this always allocated */
   free( app->spatial_lookup_table );
   
   free( app );
   MPI_Comm_free( &comm_x );
   MPI_Comm_free( &comm_t );

   /* Finalize MPI */
   MPI_Finalize();

   return 0;
}





