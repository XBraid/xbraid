/*BHEADER**********************************************************************
 * Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of XBraid.  See file COPYRIGHT for details.
 *
 * XBraid is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 ***********************************************************************EHEADER*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "c_array.h"
#include "advect_data.h"

#include "braid.h"

int main(int argc, char ** argv)
{
   int nsteps, pnr, taylorbc;
   /* enum bcType bcLeft=Dirichlet, bcRight=Extrapolation; */
   enum bcType bcLeft=Periodic, bcRight=Periodic;
   
   double h, cfl;
   double L, l2, li, tfinal;
   double amp, ph, om;
   double wave_speed, viscosity, restr_coeff=0.0, ad_coeff=0.0;
   
   int nstepsset, tfinalset, arg_index, print_usage=0, myid=0;
   int wave_no = 1, spatial_order=6;

   FILE *eun, *fp;
   
   advection_setup *kd_ = NULL;
   grid_fcn *gf_ = NULL, *exact_=NULL;

/* from drive-05.c */
   int i, level;

   braidCore  core;
/* my_App is called advection_setup, app = kd_ */
/*   my_App    *app; */
   int        max_levels;
   int        scoarsen;
   int        nrelax, nrelax0;
   double     tol;
   int        cfactor, cfactor0;
   int        max_iter;
   int        fmg;

   MPI_Comm    comm, comm_t; /* no spatial MPI decomposition */
   int         num_procs;
   /* int         xcolor, tcolor; */
   double      mystarttime, myendtime, mytime;

   /* We consider a 2D problem. */
   /* int ndim = 2; */

   /* /\* diffusion coefficient *\/ */
   /* double K;  */

   /* int nx, ny, nlx, nly; */
   /* int nt;  nt = nsteps  */
   double c;

   int pt; /* AP: what is this? */

   int n_pre, n_post;
   int rap, skip;
   
   double tol_x[2], tol_x_coarse;

   int write, vis;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
      
   /* Default parameters. */
   comm                = MPI_COMM_WORLD;
   comm_t              = comm;
   max_levels          = 2; /* This is where you specify the number of levels in MG */
   scoarsen            = 0; /* coarsen in space? */
   nrelax              = 1;
   nrelax0             = -1;
   tol                 = 1.0e-09; /* relative (or absolute) residual stopping criterion */
   cfactor             = 2;
   cfactor0            = -1;
   max_iter            = 50;
   fmg                 = 0;
   /* K                   = 1.0; */
   /* nt                  = 32; */
   c                   = 0.15;
   /* sym                 = 0; */
   /* px                  = 1; */
   /* py                  = 1; */
   pt                  = 1;
   n_pre               = 1;
   n_post              = 1;
   rap                 = 1;
   /* relax               = 3; */
   skip                = 1;
   tol_x[0]            = 1.0e-09;
   tol_x[1]            = 1.0e-09;
   tol_x_coarse        = 1.0e-09;
   write               = 0;
   vis                 = 0;

   MPI_Comm_rank( comm, &myid );
   MPI_Comm_size( comm, &num_procs );

/* Default problem parameters */
   L = 1.0; /* Domain length*/
   wave_speed = 1.0; /* wave speed */
   viscosity = 0.0;  /* viscosity */
   
/*!**  Twilight testing parameters*/
   amp  = 0.8;
   ph   = 0.17;
   wave_no = 1;
   
   om   = wave_no*2.0*M_PI; /* wave_no is the number of wave lengths in the domain */
   
/*!** exact solution
! pnr == 1:
!!$          x = (i-1)*h
!!$          w(i) = sin(pi*x*x+ph)*cos(t)
! pnr ==2:
!!$          x = (i-1)*h
!!$          w(i) = sin(om*(x-t)+ph)
*/
/* BC: pnr = 1 or 2, Dirichlet in u on the left, extrapolate u on the right */
   pnr = 2;
   
   cfl    = 0.5;
   h      = 0.01;
   tfinal = 1.0;

   nstepsset = 0;
   tfinalset = 0;

/*! time-dependent boundary data
! taylorbc = 0: assign exact boundary data at intermediate stages
! taylorbc = 1: Approximate boundary data at stages based on Gottlieb's recipe
! taylorbc = 3: Solve ODE for boundary data*/
   taylorbc = 1;

   /* Parse command line */
   if (argc>1)
   {
      arg_index = 1;
      while( arg_index < argc ){
         if( strcmp(argv[arg_index], "-dx") == 0 ){
            arg_index++;
            h = atof(argv[arg_index++]);
         }
         else if( strcmp(argv[arg_index], "-cfl") == 0 ){
            arg_index++;
            cfl = atof(argv[arg_index++]);
         }
         else if( strcmp(argv[arg_index], "-nu") == 0 ){
            arg_index++;
            viscosity = atof(argv[arg_index++]);
         }
         else if( strcmp(argv[arg_index], "-wn") == 0 ){
            arg_index++;
            wave_no = atoi(argv[arg_index++]);
            om = wave_no*2.0*M_PI;
         }
         else if( strcmp(argv[arg_index], "-nsteps") == 0 ){
            arg_index++;
            nsteps = atoi(argv[arg_index++]);
            nstepsset = 1;
         }
         else if( strcmp(argv[arg_index], "-order") == 0 ){
            arg_index++;
            spatial_order = atoi(argv[arg_index++]);
         }
         else if( strcmp(argv[arg_index], "-tfinal") == 0 ){
            arg_index++;
            tfinal = atof(argv[arg_index++]);
         }
         else if( strcmp(argv[arg_index], "-ad") == 0 ){
            arg_index++;
            ad_coeff = atof(argv[arg_index++]);
         }
         /* else if( strcmp(argv[arg_index], "-rc") == 0 ){ */
         /*    arg_index++; */
         /*    restr_coeff = atof(argv[arg_index++]); */
         /* } */
         else if( strcmp(argv[arg_index], "-tbc") == 0 ){
            arg_index++;
            taylorbc = atoi(argv[arg_index++]);
         }
         else if( strcmp(argv[arg_index], "-help") == 0 || strcmp(argv[arg_index], "--help") == 0 || 
                  strcmp(argv[arg_index], "-h") == 0)
         {
            print_usage = 1;
            break;
         }
         else
         {
            printf("unknown argument: %s\n", argv[arg_index]);
            print_usage = 1;
            break;
            /* arg_index++; */
         }
      }
   }
   
   if((print_usage) && (myid == 0)){
      printf("\n");
      printf("Solve the 1-D advection-diffusion equation with a SBP finite difference method and 4th order explicit RK:\n");
      printf(" du/dt + du/dx = nu d^2u/dx^2 + f(x,t), 0<x<1, t>0,\n u(0,t)=g(t),\n u(x,0)=h(x).\n");
      printf("\nUsage: %s [<options>]\n", argv[0]);
      printf("\n");
      printf("  -dx  <float>    : grid size (default 0.01)\n");
      printf("  -cfl <float>    : cfl-number (default 0.5)\n");
      printf("  -nu  <float>    : viscosity (>=0, default 0.0)\n");
      printf("  -nsteps <int>   : number of time steps (positive) (default tfinal/dt)\n");
      printf("  -order <int>    : spatial order of accuracy (positive, even, <=6) (default 6)\n");
      printf("  -tfinal <float> : end time (default 1.0)\n");
      printf("  -wn <int>       : wave number in exact solution (default 1)\n");
      printf("  -ad <float>     : artificial dissipation coefficient for coarse grids (default 0.0)\n");
      printf("  -tbc <int>      : treatment of bndry forcing at intermediate stages (0,1, or 3) (default 1)\n");
      printf("\n");
/* MPI_Finalize(); */
      return(0);
   }

   if (!(taylorbc==0 || taylorbc==1 || taylorbc==3))
   {
      printf("ERROR unknown taylorbc = %i\n", taylorbc);
      exit(-1);
   }
      
/* open file for saving solution error data */
   eun = fopen("err.dat","w");

/* setup solver meta-data */
   kd_ = malloc(sizeof(advection_setup));
   init_advection_solver(h, amp, ph, om, pnr, taylorbc, L, cfl, nstepsset, nsteps, tfinal, wave_speed,
                         viscosity, bcLeft, bcRight, max_iter, tol, restr_coeff, ad_coeff, spatial_order,
                         kd_);
   
#define bcnr(i) compute_index_1d(kd_->bcnr_, i)    

   printf("------------------------------\n");
   printf("Viscosity (nu): %e\n", kd_->nu_coeff);
   printf("Problem number (pnr): %i\n", kd_->pnr);
   printf("Boundary treatment: bcnr(left, right): %i, %i\n", bcnr(1), bcnr(2));
   printf("Wave number in exact solution: %i\n", wave_no);
   printf("Treatment of time-dependent bndry data: %i\n", kd_->taylorbc);
   printf("Solving to time %e using %i steps\n",kd_->tstop, kd_->nsteps);
   printf("Time step on finest grid is %e\n",kd_->dt_fine);
   printf("Finest grid has spacing h=%e with n=%i grid points\n", kd_->h_fine, kd_->n_fine);
   printf("Artificial damping coefficient for coarse grids %e\n", kd_->ad_coeff);
   printf("Spatial order of accuracy: %i\n", kd_->spatial_order);
   /* printf("Restriction coefficient for 2nd undivided difference %e\n", kd_->restr_coeff); */
   printf("------------------------------\n");

/* Start timer. */
   mystarttime = MPI_Wtime();

/* seed the random number generator */
   srand(1);
   
/* nt = nsteps : number of time steps */
   braidInit(comm, comm_t, kd_->tstart, kd_->tstop, kd_->nsteps, kd_,
             explicit_rk4_stepper, init_grid_fcn, copy_grid_fcn, free_grid_fcn, sum_grid_fcn, dot_grid_fcn, 
             save_grid_fcn, gridfcn_BufSize, gridfcn_BufPack, gridfcn_BufUnpack,
             &core);

   braidSetLoosexTol( core, 0, tol_x[0] );
   braidSetLoosexTol( core, 1, tol_x_coarse );

   braidSetTightxTol( core, 0, tol_x[1] );

/* set max number of MG levels */
   braidSetMaxLevels( core, max_levels );

   braidSetNRelax(core, -1, nrelax);
   if (nrelax0 > -1)
   {
      braidSetNRelax(core,  0, nrelax0);
   }

   /*braidSetAbsTol(core, tol*sqrt(px*nlx*py*nly*(nt+1)) );*/
   /* braidSetAbsTol(core, tol/sqrt(dx*dy*dt)); */

   /* braidSetAbsTol(core, tol); */

   braidSetRelTol(core, tol);

/* AP: this is probably related to grid coarsening in time */
   braidSetCFactor(core, -1, cfactor);
   if( cfactor0 > -1 ){
      /* Use cfactor0 on all levels until there are < cfactor0 points
       * on each processor. */
      level = (int) (log10((nsteps + 1) / pt) / log10(cfactor0));
      for( i = 0; i < level; i++ )
         braidSetCFactor(core,  i, cfactor0);
   }
   
   braidSetMaxIter(core, max_iter);
   if (fmg)
   {
      braidSetFMG(core);
   }
   
/* this is where the coarsen and refine routines are defined */
   if (scoarsen)
   {
      braidSetSpatialCoarsen(core, gridfcn_Coarsen);
      braidSetSpatialRefine(core, gridfcn_Refine);
   }
   
   /* control how often my save_grid_fcn routine is called. */
/* 0 is never, 1 is at convergence for the finest level, 2 is after every iteration on every level */
   braidSetAccessLevel(core, 2);

   braidDrive(core);

   /* Stop timer. */
   myendtime = MPI_Wtime();
   mytime    = myendtime - mystarttime;

   braidPrintStats(core);

/* my stuff... */
   printf("------------------------------\n");
   printf("Time-stepping completed. Solved to time t: %e using %i time steps\n", kd_->tstop, kd_->nsteps);
   printf("Fine grid h: %e, cfl=dt/h: %e, T*dt/h^2: %e\n", kd_->h_fine, kd_->dt_fine/kd_->h_fine, 
          kd_->tstop*kd_->dt_fine/kd_->h_fine/kd_->h_fine );

   if (kd_->sol_copy)
   {
/*  get a pointer to the final solution from the advection_setup structure */
      gf_ = kd_->sol_copy;

/* allocate storage for exact solution of the same size as gf_ */
      copy_grid_fcn(kd_, gf_, &exact_);

/* ! evaluate solution error */
      exact1( exact_, kd_->tstop, kd_ );
/* get exact bndry data */
      bdata( exact_, kd_->tstop, kd_);

      evaldiff( gf_, exact_, &l2, &li );

      printf("------------------------------\n");
   
      printf("Solution error in maximum norm, bndry error\n");
   
#define vsol(i) compute_index_1d(gf_->vsol_, i)   
#define ex_vsol(i) compute_index_1d(exact_->vsol_, i)   
      printf("time: %e, sol-err: %e, bndry-err: %e\n", kd_->tstop, li, fabs(ex_vsol(1)-vsol(1)));
      printf("nu/h_fine=%e\n", kd_->nu_coeff/kd_->h_fine);
      
      printf("------------------------------\n");

      printf("Saving ...\n");

      fp=fopen("sol.dat","w");
      for ( i=1; i<=gf_->n; i++)
      {
         fprintf(fp, "%e %e %e\n", (i-1)*gf_->h, gf_->sol[i], gf_->sol[i]-exact_->sol[i]);
      }
      fclose(fp);
      
      printf("done.\n");
      printf("------------------------------\n");

/* free tmp storage */
      free_grid_fcn( kd_, exact_ );
   }

}


             
