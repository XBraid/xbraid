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
 * Example:       ex-03-serial.c
 *
 * Interface:     C
 * 
 * Requires:      hypre 
 *
 * Compile with:  make ex-03-serial
 *
 * Help with:     ex-03-serial -help
 *
 * Sample run:    mpirun -np 9 ex-03-serial -pgrid 3 3 -nx 33 33
 *
 * Description:   Solves the 2D heat equation on a regular grid in space and time,
 *                using backward Euler in time and classic second order 
 *                finite-differencing in space.
 *
 *                The key difference withe ex-03.c is that this is a strictly 
 *                sequential time integration code, created for readers to 
 *                compare with ex-03.c. 
 *
 *                For more details on the discretization, see the header comment in ex-03.c.
 **/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "_hypre_utilities.h"
#include "HYPRE_sstruct_ls.h"
#include "_hypre_sstruct_mv.h"

#include "ex-03-lib.c"

int main (int argc, char *argv[])
{
   /* Declare variables -- variables explained when they are set below */
   MPI_Comm    comm;
   HYPRE_SStructVector u, e;
   simulation_manager * man;

   int print_usage = 0;
   int object_type = HYPRE_STRUCT;
   int i, arg_index, myid, num_procs, iters_taken, max_iters_taken;
   int ndim, nx, ny, nlx, nly, nt, forcing, ilower[2], iupper[2];
   double K, tstart, tstop, dx, dy, dt, cfl, tol;
   double myendtime, mystarttime, mytime, maxtime;
   double disc_err, max_disc_err, max_disc_err_time;
   int max_iter, px, py, pi, pj, output_files, vis, max_disc_err_iter;
   char filename[255], filename_mesh[255], filename_err[255], filename_sol[255];

   /* Initialize MPI */
   MPI_Init(&argc, &argv);

   /* Default parameters. */
   comm                = MPI_COMM_WORLD;
   ndim                = 2;       /* Two dimensional problem */
   forcing             = 0;       /* Boolean, if 1 use a nonzero forcing term, if 0 use a 0 forcing term */
   K                   = 1.0;     /* Diffusion coefficient */
   nx                  = 17;      /* number of points in the x-dim */
   ny                  = 17;      /* number of points in the y-dim */
   nlx                 = 17;      /* number of point ~local~ to this processor in the x-dim */
   nly                 = 17;      /* number of point ~local~ to this processor in the y-dim */
   tstart              = 0.0;     /* global start time */
   nt                  = 64;      /* number of time steps */
   cfl                 = 0.30;    /* CFL = K*(dt/dx^2 + dt/dy^2) is used to define dt and t-final */
   px                  = 1;       /* my processor number in the x-direction, px*py=num procs in space */
   py                  = 1;       /* my processor number in the y-direction, px*py=num procs in space */
   max_iter            = 50;      /* Maximum number of iterations to use inside of PFMG */
   tol                 = 1.0e-09; /* PFMG halting tolerance */
   output_files        = 0;       /* Boolean, if 1 output the norm of the discretization error to a file for each time step */
   vis                 = 0;       /* Boolean, if 1 output GLVIS files of error, solution and true solution */

   MPI_Comm_rank( comm, &myid );
   MPI_Comm_size( comm, &num_procs );

   /* Parse command line */
   arg_index = 0;
   while( arg_index < argc ){
      if( strcmp(argv[arg_index], "-pgrid") == 0 ){
         arg_index++;
         px = atoi(argv[arg_index++]);
         py = atoi(argv[arg_index++]);
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
      else if( strcmp(argv[arg_index], "-forcing") == 0 ){
         arg_index++;
         forcing = 1;
      }
      else if( strcmp(argv[arg_index], "-pfmg_mi") == 0 ){
         arg_index++;
         max_iter = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-pfmg_tol") == 0 ){
          arg_index++;
          tol = atof(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-output_files") == 0 ){
         arg_index++;
         output_files = 1;
      }
      else if( strcmp(argv[arg_index], "-output_vis") == 0 ){
         arg_index++;
         vis = 1;
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
      printf(" Usage parameters\n");
      printf(" ---------------------------------------\n");
      printf("  -pgrid  <px py>                  : processors in each dimension (default: 1 1)\n");
      printf("  -nx  <nx ny>                     : 2D spatial problem size of form 2^k+1, 2^k+1 (default: 17 17)\n");
      printf("  -nt  <n>                         : number of time steps (default: 32)\n"); 
      printf("  -cfl <cfl>                       : CFL number to run (default: 0.30)\n"); 
      printf("                                     Note: CFL = K*(dt/dx^2 + dt/dy^2) is used to define dt and t-final\n");
      printf("  -forcing                         : consider non-zero RHS b(x,y,t) = -sin(x)*sin(y)*(sin(t)-2*cos(t))\n");
      printf("  -pfmg_mi <max_iter>              : maximum number of PFMG iterations (default: 50)\n"); 
      printf("  -pfmg_tol <tol>                  : loose and tight stopping tolerance for PFMG (default: 1e-09 1e-09)\n"); 
      printf("  -output_files                    : save the solution/error/error norms to files\n");
      printf("                                     frequency of file accesss is set by access_level\n");
      printf("  -output_vis                      : save the error for GLVis visualization\n");
      printf("                                     frequency of file accesss is set by access_level\n");
      printf("\n");
   }

   if( print_usage ){
      MPI_Finalize();
      return (0);
   }

   /* Check the processor grid (px x py = num_procs?). */
   if( (px*py) != num_procs)
   {
       if( myid == 0 )
           printf("Error: px x py is not equal to the number of processors!\n");
       MPI_Finalize();
       return (0);
   }

   /* Determine position (pi, pj)  in the 2D processor grid, 
    * 0 <= pi < px,   0 <= pj < py */
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
    * Set up the manager 
    * ----------------------------------------------------------------- */
   man               = (simulation_manager *) malloc(sizeof(simulation_manager));
   man->comm         = comm;
   man->forcing      = forcing;
   man->K            = K;
   man->dim_x        = ndim;
   man->nlx          = nlx;
   man->nly          = nly;
   man->nx           = nx;
   man->ny           = ny;
   man->tstart       = tstart;
   man->tstop        = tstop;
   man->nt           = nt;
   man->dx           = dx;
   man->dy           = dy;
   man->dt           = dt;
   man->px           = px;
   man->py           = py;
   man->pi           = pi;
   man->pj           = pj;
   man->ilower[0]    = ilower[0];
   man->ilower[1]    = ilower[1];
   man->iupper[0]    = iupper[0];
   man->iupper[1]    = iupper[1];
   man->object_type  = object_type;
   man->max_iter     = max_iter;
   man->tol          = tol;
   man->output_vis   = vis;
   man->output_files = output_files;
   man->explicit     = 0;

   /* Set up the variable type, grid, stencil and matrix graph. */
   man->vartype           = HYPRE_SSTRUCT_VARIABLE_CELL;
   setUp2Dgrid( comm, &(man->grid_x), man->dim_x,
                man->ilower, man->iupper, man->vartype, 1 );
   set5ptStencil( &(man->stencil), man->dim_x );
   setUpGraph( comm, &(man->graph), man->grid_x, man->object_type, 
               man->stencil );

   /* Set up initial state vector */
   set_initial_condition(man, &u, 0.0);

   /* Set up error vector */
   initialize_vector(man, &e);
 
   /* Set up the matrix */
   setUpImplicitMatrix( man );
   setUpStructSolver( man, u, u );

   if( myid == 0 ) {
      printf("\n  --------------------- \n");
      printf("  Begin simulation \n");
      printf("  --------------------- \n\n");
   }
   
   /* Run a simulation */
   mystarttime = MPI_Wtime();
   max_iters_taken = -1;
   max_disc_err = -1.;
   tstart = 0.0;
   tstop = tstart + dt;
   for(i = 0; i < man->nt; i++)
   {
      if( myid == 0 ) {
         if( i % 50 == 0)  printf("  Taking iteration %d...\n", i);
      }
      
      /* Take Step */
      take_step(man, u, NULL, u, tstart, tstop, &iters_taken);
      if( i < man->nt-1){
         tstart = tstop;
         tstop = tstop + dt;
      }

      /* Output */
      compute_disc_err(man, u, tstop, e, &disc_err);
      if( man->output_files ){
         sprintf(filename, "%s.timestep%05d", "ex-03-serial.error_norm", i);
         output_error_file(man, tstop, disc_err, filename); 
      }

      /* Store PFMG iters taken and maximum discretization error */
      max_iters_taken = max_i(max_iters_taken, iters_taken);
      max_disc_err = max_d(max_disc_err, disc_err);
      if( max_disc_err == disc_err) {
         max_disc_err_iter = i;
         max_disc_err_time = tstop;
      }

   }
   myendtime = MPI_Wtime();

   /* Print some additional statistics */
   mytime    = myendtime - mystarttime;
   MPI_Reduce( &mytime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, comm );
   if( myid == 0 )
   {
      printf("\n  --------------------- \n");
      printf("  End simulation \n");
      printf("  --------------------- \n\n");
      
      printf("  Start time                    %1.5e\n", 0.0);
      printf("  Stop time                     %1.5e\n", tstop);
      printf("  Time step size                %1.5e\n", man->dt);
      printf("  Time steps taken:             %d\n\n", i);
      printf("  Spatial grid size:            %d,%d\n", man->nx, man->ny);
      printf("  Spatial mesh width (dx,dy):  (%1.2e, %1.2e)\n", man->dx, man->dy);           
      printf("  CFL ratio 2dt/(dx^2 + dy^2):  %1.2e\n\n",
            man->K*((man->dt)/((man->dx)*(man->dx)) + (man->dt)/((man->dy)*(man->dy))));
      printf("  Run time:                     %1.2e\n", maxtime);
      printf("  Max PFMG Iterations:          %d\n", max_iters_taken);
      printf("  Discr. error at final time:   %1.4e\n", disc_err);
      printf("  Max discr. error:             %1.2e\n",max_disc_err);
      printf("     found at iteration:        %d\n", max_disc_err_iter);
      printf("     found at time:             %1.2e\n\n", max_disc_err_time);
   }      

   /* Visualize final error */
   if( man->output_vis){
      sprintf(filename_mesh, "%s", "ex-03-serial_mesh");
      sprintf(filename_err, "%s", "ex-03-serial_err_tstop");
      sprintf(filename_sol, "%s", "ex-03-serial_sol_tstop");
      output_vis(man, u, tstop, filename_mesh, filename_err, filename_sol);
   }

   /* Free memory */
   HYPRE_SStructVectorDestroy( u );
   HYPRE_SStructVectorDestroy( e );
   HYPRE_SStructGridDestroy( man->grid_x );
   HYPRE_SStructStencilDestroy( man->stencil );
   HYPRE_SStructGraphDestroy( man->graph );
   HYPRE_SStructMatrixDestroy( man->A );
   HYPRE_StructPFMGDestroy( man->solver );
   free(man);
   
   MPI_Finalize();
   return 0;
}
