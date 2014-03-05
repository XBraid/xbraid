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
   /* enum bcType bcLeft=Dirichlet, bcRight=Extrapolation; */
   enum bcType bcLeft=Periodic, bcRight=Periodic;
   
   double h, cfl, bdata[2];
   double L, l2, li, lib, tfinal, linf;
   double t0, dtf, dtc;
   
   double amp, ph, om;
   double wave_speed, viscosity;
   
   int nstepsset, tfinalset, arg_index, print_usage=0, myid=0;

   FILE *fp=NULL;
   
   advection_setup *kd_ = NULL;
   grid_fcn *gf_ = NULL, *coarse_gf_ = NULL, *fine_gf_ = NULL, *wx_gf_ = NULL, *wxx_ = NULL;
   grid_fcn *exact_=NULL, *wt_gf_=NULL;
   
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
   L = 1.0; /*  Domain length*/
   wave_speed = 1.0;
   viscosity = 0.0;

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
      else if( strcmp(argv[arg_index], "-nu") == 0 ){
         arg_index++;
         viscosity = atof(argv[arg_index++]);
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
      printf("  -nu  <float>    : viscosity (>=0, default 0.0)\n");
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
   init_advection_solver(h, amp, ph, om, pnr, taylorbc, L, cfl, nstepsset, nsteps, tfinal, 
                         wave_speed, viscosity, bcLeft, bcRight, 0, 0.0, kd_); /* stopping criteria not used */
   
/* for testing */
   kd_->tstart = 1.0;

#define bcnr(i) compute_index_1d(kd_->bcnr_, i)    

   printf("------------------------------\n");
   printf("Testing init, coarsen and refine routines at t=%e\n", kd_->tstart);
   printf("Viscosity (nu): %e\n", kd_->nu_coeff);
   printf("Problem number (pnr): %i\n", kd_->pnr);
   printf("Boundary treatment: bcnr(left, right): %i, %i\n", bcnr(1), bcnr(2));
   printf("Treatment of time-dependent bndry data: %i\n", kd_->taylorbc);
   printf("Solving to time %e using %i steps\n",kd_->tstop, kd_->nsteps);
   printf("Time step on fine grid is %e\n",kd_->dt_fine);
   printf("Finest grid has spacing h=%e with n=%i grid points\n", kd_->h_fine, kd_->n_fine);
   printf("------------------------------\n");

#define SQR(x) ((x)*(x))

/* Test truncation error of spatial operator */
/* 1) call init_grid_fcn() for time t=0 to define a grid function on the finest grid */

   t0 = kd_->tstart;
   init_grid_fcn(kd_, t0, &gf_);
/* make grid function, t=1 gives dw/dt=0 */
   init_grid_fcn(kd_, 1.0, &wx_gf_);
   init_grid_fcn(kd_, 1.0, &wt_gf_);
   init_grid_fcn(kd_, 1.0, &wxx_);

/* evaluate dw/dx */
   dwdx( gf_, wx_gf_, kd_ );

/* tmp storage */
   copy_grid_fcn( kd_, gf_, &exact_ );

/* calculate exact dudx, save in exact_ grid function */
   exact_x( exact_, t0, kd_ );

/* evaluate error in dw/dx */
   evalerr1(wx_gf_, exact_, &l2, &linf);

/* interior points only */
   li = 0.0;
   for (i=1+kd_->nb; i<=gf_->n-kd_->nb; i++)
   {
      if (fabs(wx_gf_->sol[i] - exact_->sol[i]) > li)
         li = fabs(wx_gf_->sol[i] - exact_->sol[i]);
   }
   
   lib=0; /* near bndry */
/* left bndry */
   for (i=1; i<=kd_->nb; i++)
   {
      if (fabs(wx_gf_->sol[i] - exact_->sol[i]) > lib)
         lib = fabs(wx_gf_->sol[i] - exact_->sol[i]);
   }
/* right bndry */
   for (i=gf_->n - kd_->nb + 1; i<=gf_->n; i++)
   {
      if (fabs(wx_gf_->sol[i] - exact_->sol[i]) > lib)
         lib = fabs(wx_gf_->sol[i] - exact_->sol[i]);
   }
   printf("Difference (Dw - dw/dx): all-i = %e, all-2 = %e, interior-i = %e, bndry-i = %e\n", linf, l2, li, lib);

/* save Dw and dw/dx */
   fp = fopen("wx.dat","w");
   for (i=1; i<=gf_->n; i++)
   {
      fprintf(fp, "%e  %e\n", exact_->sol[i], wx_gf_->sol[i]);
   }
   fclose(fp);

/* evaluate d2w/dx2 */
   d2wdx2( gf_, wxx_, kd_ );
/* ! evaluate solution error, save exact solution in kd_ workspace array */
   exact_xx( exact_, t0, kd_ );
   
/* evaluate error in dw/dx */
   evalerr1(wxx_, exact_, &l2, &linf);

/* interior points only */
   li=0;  /* interior */
   lib=0; /* near bndry */
   
   for (i=1+kd_->nb2; i<=gf_->n - kd_->nb2; i++)
   {
      if (fabs(wxx_->sol[i] - exact_->sol[i]) > li)
         li = fabs(wxx_->sol[i] - exact_->sol[i]);
   }

/* left bndry */
   for (i=1; i<=kd_->nb; i++)
   {
      if (fabs(wxx_->sol[i] - exact_->sol[i]) > lib)
         lib = fabs(wxx_->sol[i] - exact_->sol[i]);
   }
/* right bndry */
   for (i=gf_->n - kd_->nb + 1; i<=gf_->n; i++)
   {
      if (fabs(wxx_->sol[i] - exact_->sol[i]) > lib)
         lib = fabs(wxx_->sol[i] - exact_->sol[i]);
   }
   printf("Difference (Gw - d2w/dx2): all-i = %e, all-2 = %e, interior-i = %e, bndry-i = %e\n", linf, l2, li, lib);

/* calculate exact dudt, save in exact_ grid function */
   exact_t( exact_, t0, kd_ );

/* get boundary data */
   twbndry1( &(bdata[0]), &(bdata[1]), 1, t0, 0.0, kd_ );

/* evaluate dwdt, tstage is for evaluating forcing */
   dwdt( gf_, wt_gf_, t0, bdata, kd_ );

/* evaluate error in dw/dx */
   evalerr1(wt_gf_, exact_, &l2, &linf);

/* interior points only */
   li=0;  /* interior */
   lib=0; /* near bndry */
   
   for (i=1+kd_->nb2; i<=gf_->n - kd_->nb2; i++)
   {
      if (fabs(wt_gf_->sol[i] - exact_->sol[i]) > li)
         li = fabs(wt_gf_->sol[i] - exact_->sol[i]);
   }

/* left bndry */
   for (i=1; i<=kd_->nb; i++)
   {
      if (fabs(wt_gf_->sol[i] - exact_->sol[i]) > lib)
         lib = fabs(wt_gf_->sol[i] - exact_->sol[i]);
   }
/* right bndry */
   for (i=gf_->n - kd_->nb + 1; i<=gf_->n; i++)
   {
      if (fabs(wt_gf_->sol[i] - exact_->sol[i]) > lib)
         lib = fabs(wt_gf_->sol[i] - exact_->sol[i]);
   }

   printf("Difference (Tw - dw/dt): all-i = %e, all-2 = %e, interior-i = %e, bndry-i = %e\n", linf, l2, li, lib);

/* Test init-coarsen-restrict */

/* 1) call init_grid_fcn() for time t=0 to define a grid function on the finest grid */
   t0 = kd_->tstart;
   init_grid_fcn(kd_, t0, &gf_);

/* 2) call gridfcn_Coarsen() to coarsen to define a coarse version of that grid function */
   dtf = kd_->dt_fine;
   dtc = 2*dtf;
   gridfcn_Coarsen(kd_, t0, t0, t0+dtf, t0, t0+dtc, gf_, &coarse_gf_);
   
/* 3) call gridfcn_Refine() to refine the coarsened grid function and compare to the original function */
   gridfcn_Refine(kd_, t0, t0, t0+dtf, t0, t0+dtc, coarse_gf_, &fine_gf_);

   printf("------------------------------\n");
/* evaluate interpolation error */
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

/* get exact bndry data (bdataL). Note: when s=1, dt is not used */
   twbndry1( &(bdata[0]), &(bdata[1]), 1, kd_->tstop, 0.0, kd_);

/*   printf("Saving ...\n");*/
    /* open(21,file='sol.bin',form='unformatted') */
    /* fprintf(21) n, h, dt, t, (sol(i),i=1,n), (current(i),i=1,n) */
    /* close(21) */
/*   printf("done.\n");*/
/*   printf("------------------------------\n");*/
   
   free_grid_fcn( kd_, exact_ );

}


             
