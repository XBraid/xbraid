#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "c_array.h"
#include "advect_data.h"

#include "warp.h"

int main(int argc, char ** argv)
{
   int nsteps, pnr, taylorbc;
   /* const int bdatareset=1000; */
   
   double h, cfl, bdataL, bdataR;
   double L, l2, li, tfinal;
   double t0, dtf, dtc;
   
   double amp, ph, om;
   
   int nstepsset, tfinalset, arg_index, print_usage=0, myid=0;

   
   advection_setup *kd_ = NULL;
   grid_fcn *gf_ = NULL, *coarse_gf_ = NULL, *fine_gf_ = NULL;

/* from drive-05.c */
   int i;

/*   warp_Core  core;*/
/* my_App is called advection_setup, app = kd_ */
/*   my_App    *app; */
   int        max_levels;
   int        nrelax, nrelax0;
   double     tol;
   int        cfactor, cfactor0;
   int        max_iter;
   int        fmg;

   /* MPI_Comm    comm, comm_t; /\* no spatial MPI decomposition *\/ */
   /* int         num_procs; */
   /* int         xcolor, tcolor; */

   /* We consider a 2D problem. */
   /* int ndim = 2; */

   /* /\* diffusion coefficient *\/ */
   /* double K;  */

   /* int nx, ny, nlx, nly; */
   /* int nt;  nt = nsteps  */
   double c;

   int pt; /* AP: what is this? */

   int n_pre, n_post;
   int rap, relax, skip;
   
   double tol_x[2], tol_x_coarse;

   int write, vis;

   /* Initialize MPI */
   /* MPI_Init(&argc, &argv); */
      
   /* Default parameters. */
   /* comm                = MPI_COMM_WORLD; */
   /* comm_t              = comm; */
   max_levels          = 2; /* AP changed from 1 */
   nrelax              = 1;
   nrelax0             = -1;
   tol                 = 1.0e-09;
   cfactor             = 2;
   cfactor0            = -1;
   max_iter            = 30;
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
   relax               = 3;
   skip                = 1;
   tol_x[0]            = 1.0e-09;
   tol_x[1]            = 1.0e-09;
   tol_x_coarse        = 1.0e-09;
   write               = 0;
   vis                 = 0;

   /* MPI_Comm_rank( comm, &myid ); */
   /* MPI_Comm_size( comm, &num_procs ); */

/* Default problem parameters */
/*!**  Domain length*/
   L = 1.0;
   
/*!**  Twilight testing parameters*/
   amp  = 0.8;
   ph   = 0.17;
   om = 5.0;
   
/*!** exact solution
! pnr == 1:
!!$          x = (i-1)*h
!!$          w(1,i) = sin(pi*x*x+ph)*cos(t)
!!$          w(2,i) = amp*cos(pi*x)**2*sin(t+ph)
! pnr ==2:
!!$          x = (i-1)*h
!!$          w(1,i) = sin(om*(x-t)+ph)
!!$          w(2,i) = cos(om*(x+t))
*/
/* BC: pnr = 1 or 2, Dirichlet in u on the left, extrapolate u on the right */
   pnr = 2;
   
   cfl    = 1.0;
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
   arg_index = 0;
   while( arg_index < argc ){
      if( strcmp(argv[arg_index], "-dx") == 0 ){
         arg_index++;
         h = atof(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-cfl") == 0 ){
         arg_index++;
         cfl = atof(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-nsteps") == 0 ){
          arg_index++;
          nsteps = atoi(argv[arg_index++]);
          nstepsset = 1;
      }
      else if( strcmp(argv[arg_index], "-tfinal") == 0 ){
          arg_index++;
          tfinal = atof(argv[arg_index++]);
      }
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
         arg_index++;
      }
   }

   if((print_usage) && (myid == 0)){
      printf("\n");
      printf("Solve the 1-D advection equation with a SBP finite difference method and 4th order explicit RK:\n");
      printf(" du/dt + du/dx = f(x,t), 0<x<1, t>0,\n u(0,t)=g(t),\n u(x,0)=h(x).\n");
      printf("\nUsage: %s [<options>]\n", argv[0]);
      printf("\n");
      printf("  -dx  <float>    : grid size (default 0.01)\n");
      printf("  -cfl <float>    : cfl-number (default 1.0)\n");
      printf("  -nsteps <int>   : number of time steps (positive) (default tfinal/dt)\n");
      printf("  -tfinal <float> : end time (default 1.0)\n");
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
      
/* setup solver meta-data */
   kd_ = malloc(sizeof(advection_setup));
   init_advection_solver(h, amp, ph, om, pnr, taylorbc, L, cfl, nstepsset, nsteps, tfinal, kd_);
   
#define bcnr(i) compute_index_1d(kd_->bcnr_, i)    

   printf("------------------------------\n");
   printf("Testing init, coarsen and refine routines...\n");
   printf("Problem number (pnr): %i\n", kd_->pnr);
   printf("Boundary treatment: bcnr(left, right): %i, %i\n", bcnr(1), bcnr(2));
   printf("Treatment of time-dependent bndry data: %i\n", kd_->taylorbc);
   printf("Solving to time %e using %i steps\n",kd_->tstop, kd_->nsteps);
   printf("Time step is %e\n",kd_->dt);
   printf("Finest grid has spacing h=%e with n=%i grid points\n", kd_->h_fine, kd_->n_fine);
   printf("------------------------------\n");

/* Plan: */
/* 1) call init_grid_fcn() for time t=0 to define a grid function on the finest grid */
   t0 = kd_->tstart;
   init_grid_fcn(kd_, t0, &gf_);

/* 2) call gridfcn_Coarsen() to coarsen to define a coarse version of that grid function */
   dtf = kd_->dt;
   dtc = 2*dtf;
   gridfcn_Coarsen(kd_, t0, t0, t0+dtf, t0, t0+dtc, gf_, &coarse_gf_);
   
/* 3) call gridfcn_Refine() to refine the coarsened grid function and compare to the original function */
   gridfcn_Refine(kd_, t0, t0, t0+dtf, t0, t0+dtc, coarse_gf_, &fine_gf_);

   printf("------------------------------\n");
/* evaluate interpolation error */
#define SQR(x) ((x)*(x))
   li=0;
   l2=0;
   for (i=1; i<=gf_->n; i++)
   {
      l2 += SQR(gf_->sol[i] - fine_gf_->sol[i]);
      if (fabs(gf_->sol[i] - fine_gf_->sol[i]) > li)
         li = fabs(gf_->sol[i] - fine_gf_->sol[i]);
   }
   l2 = sqrt(l2/gf_->n);
   
   printf("Difference (I - PR)u: l2 = %e, linf = %e\n", l2, li);

/* ! evaluate solution error, stick exact solution in kd_ workspace array */
   exact1( kd_->n_fine, kd_->current, kd_->h_fine, kd_->amp, kd_->ph, kd_->om, kd_->tstop, kd_->pnr );
/* get exact bndry data (bdataL) */
   twbndry1( 0.0, &bdataL, kd_->L, &bdataR, 1, kd_->tstop, kd_->dt, kd_->amp, kd_->ph, kd_->om, kd_->pnr);

/*   printf("Saving ...\n");*/
    /* open(21,file='sol.bin',form='unformatted') */
    /* fprintf(21) n, h, dt, t, (sol(i),i=1,n), (current(i),i=1,n) */
    /* close(21) */
/*   printf("done.\n");*/
/*   printf("------------------------------\n");*/
   
}


             
