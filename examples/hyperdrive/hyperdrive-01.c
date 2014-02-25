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
   
   double h, cfl, bdataL, bdataR;
   double L, l2, li, tfinal;
   double amp, ph, om;
   double wave_speed, viscosity;
   
   int nstepsset, tfinalset, arg_index, print_usage=0, myid=0;

   FILE *eun;
   
   advection_setup *kd_ = NULL;
   grid_fcn *gf_ = NULL, *exact_=NULL;

/* from drive-05.c */
   int i, level;

   warp_Core  core;
/* my_App is called advection_setup, app = kd_ */
/*   my_App    *app; */
   int        max_levels;
   int        scoarsen=1;
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
   int rap, relax, skip;
   
   double tol_x[2], tol_x_coarse;

   int write, vis;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
      
   /* Default parameters. */
   comm                = MPI_COMM_WORLD;
   comm_t              = comm;
   max_levels          = 2; /* This is where you specify the number of levels in MG */
   nrelax              = 1;
   nrelax0             = -1;
   tol                 = 1.0e-09;
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
   relax               = 3;
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
   om = 5.0;
   
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
      
/* open file for saving solution error data */
   eun = fopen("err.dat","w");

/* setup solver meta-data */
   kd_ = malloc(sizeof(advection_setup));
   init_advection_solver(h, amp, ph, om, pnr, taylorbc, L, cfl, nstepsset, nsteps, tfinal, 
                         wave_speed, viscosity, bcLeft, bcRight, kd_);
   
#define bcnr(i) compute_index_1d(kd_->bcnr_, i)    

   printf("------------------------------\n");
   printf("Viscosity (nu): %e\n", kd_->nu_coeff);
   printf("Problem number (pnr): %i\n", kd_->pnr);
   printf("Boundary treatment: bcnr(left, right): %i, %i\n", bcnr(1), bcnr(2));
   printf("Treatment of time-dependent bndry data: %i\n", kd_->taylorbc);
   printf("Solving to time %e using %i steps\n",kd_->tstop, kd_->nsteps);
   printf("Time step on finest grid is %e\n",kd_->dt_fine);
   printf("Finest grid has spacing h=%e with n=%i grid points\n", kd_->h_fine, kd_->n_fine);
   printf("------------------------------\n");

/* Start timer. */
   mystarttime = MPI_Wtime();

/* nt = nsteps : number of time steps */
   warp_Init(comm, comm_t, kd_->tstart, kd_->tstop, kd_->nsteps, kd_,
             explicit_rk4_stepper, init_grid_fcn, copy_grid_fcn, free_grid_fcn, sum_grid_fcn, dot_grid_fcn, 
             save_grid_fcn, gridfcn_BufSize, gridfcn_BufPack, gridfcn_BufUnpack,
             &core);

   warp_SetLoosexTol( core, 0, tol_x[0] );
   warp_SetLoosexTol( core, 1, tol_x_coarse );

   warp_SetTightxTol( core, 0, tol_x[1] );

/* set max number of MG levels */
   warp_SetMaxLevels( core, max_levels );

   warp_SetNRelax(core, -1, nrelax);
   if (nrelax0 > -1)
   {
      warp_SetNRelax(core,  0, nrelax0);
   }

   /* warp_SetRelTol(core, tol); */
   /*warp_SetAbsTol(core, tol*sqrt(px*nlx*py*nly*(nt+1)) );*/
   /* warp_SetAbsTol(core, tol/sqrt(dx*dy*dt)); */
   warp_SetAbsTol(core, tol);

/* AP: this is probably related to grid coarsening in time */
   warp_SetCFactor(core, -1, cfactor);
   if( cfactor0 > -1 ){
      /* Use cfactor0 on all levels until there are < cfactor0 points
       * on each processor. */
      level = (int) (log10((nsteps + 1) / pt) / log10(cfactor0));
      for( i = 0; i < level; i++ )
         warp_SetCFactor(core,  i, cfactor0);
   }
   
   warp_SetMaxIter(core, max_iter);
   if (fmg)
   {
      warp_SetFMG(core);
   }
   
/* this is where the coarsen and refine routines are defined */
   if (scoarsen)
   {
      warp_SetSpatialCoarsen(core, gridfcn_Coarsen);
      warp_SetSpatialRefine(core, gridfcn_Refine);
   }
   
/* control how often my write routine is called. How is this supposed to work??? */
   warp_SetWriteLevel(core, 1);
   

   warp_Drive(core);

   /* Stop timer. */
   myendtime = MPI_Wtime();
   mytime    = myendtime - mystarttime;

   warp_PrintStats(core);

/* my stuff... */
   printf("------------------------------\n");
   printf("Time-stepping completed. Solved to time t: %e\n", kd_->tstop);

/* tmp storage */
   init_grid_fcn(kd_, 0.0, &exact_);

/* ! evaluate solution error */
   exact1( exact_, kd_->tstop, kd_ );
/* get exact bndry data (bdataL). Note: dt not used when s=1 */
   twbndry1( &bdataL, &bdataR, 1, kd_->tstop, 0.0, kd_ );

   if (kd_->sol_copy)
   {
      
/*  get a pointer to the final solution from the advection_setup structure */
      gf_ = kd_->sol_copy;
      evalerr1( gf_, exact_, &l2, &li );

      printf("------------------------------\n");
   
      printf("Solution error in maximum norm, bndry error\n");
   
#define vsol(i) compute_index_1d(gf_->vsol_, i)   
      printf("time: %e, sol-err: %e, bndry-err: %e\n", kd_->tstop, li, fabs(bdataL-vsol(1)));
      printf("4*nu/h_fine=%e\n", 4.0*kd_->nu_coeff/kd_->h_fine);
      
      printf("------------------------------\n");

/*   printf("Saving ...\n");*/
    /* open(21,file='sol.bin',form='unformatted') */
    /* fprintf(21) n, h, dt, t, (sol(i),i=1,n), (current(i),i=1,n) */
    /* close(21) */
/*   printf("done.\n");*/
/*   printf("------------------------------\n");*/
   }

/* free tmp storage */
   free_grid_fcn( kd_, exact_ );
   
}


             
