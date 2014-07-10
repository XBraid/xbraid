#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "c_array.h"
#include "advect_data.h"

int main(int argc, char ** argv)
{
   int step, nsteps, pnr, taylorbc;
   /* enum bcType bcLeft=Dirichlet, bcRight=Extrapolation; */
   enum bcType bcLeft=Periodic, bcRight=Periodic;
   
   const int bdatareset=1000;
   
   double h, cfl;
   double L, t, l2, li, tfinal;
   double amp, ph, om;
   double wave_speed, viscosity, ad_coeff=0.0;
   int wave_no = 1, spatial_order=6;
   
   double dummy=0;
   int rfact_dummy=0;
   
   int nstepsset, tfinalset, arg_index, print_usage=0, myid=0;

   FILE *eun;
   
   advection_setup *kd_ = NULL;
   grid_fcn *gf_ = NULL, *exact_ = NULL;
      
/* Default problem parameters */

   L = 1.0; /* Domain length*/
   wave_speed = 1.0; /* wave speed */
   viscosity = 0.0;  /* viscosity */

/*!**  Twilight testing parameters*/
   amp  = 0.8;
   ph   = 0.17;
   om   = wave_no*2.0*M_PI;
   
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
      printf("Solve the 1-D advection-diffusion equation with a SBP finite difference method and 4th order explicit RK:\n");
      printf(" du/dt + du/dx = nu d^2u/dx^2 + f(x,t), 0<x<1, t>0,\n u(0,t)=g(t),\n u(x,0)=h(x).\n");
      printf("\nUsage: %s [<options>]\n", argv[0]);
      printf("\n");
      printf("  -dx  <float>    : grid size (default 0.01)\n");
      printf("  -cfl <float>    : cfl-number (default 1.0)\n");
      printf("  -nu  <float>    : viscosity (>=0, default 0.0)\n");
      printf("  -nsteps <int>   : number of time steps (positive) (default tfinal/dt)\n");
      printf("  -order <int>    : spatial order of accuracy (positive, even, <=6) (default 6)\n");
      printf("  -tfinal <float> : end time (default 1.0)\n");
      printf("  -wn <int>       : wave number in exact solution (default 1)\n");
      printf("  -ad <float>     : artificial dissipation coefficient (default 0.0)\n");
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

#define bcnr(i)    compute_index_1d(kd_->bcnr_, i)    
#define vsol(i)    compute_index_1d(gf_->vsol_, i)   
#define ex_vsol(i) compute_index_1d(exact_->vsol_, i)   

/* solver meta-data */
   kd_ = malloc(sizeof(advection_setup));
/* Note: stopping criteria not used by this solver */
   init_advection_solver(h, amp, ph, om, pnr, taylorbc, L, cfl, nstepsset, nsteps, tfinal, 
                         wave_speed, viscosity, bcLeft, bcRight, 0, 0.0, 0.0, ad_coeff, spatial_order, 
                         kd_); 
   
/* create solution vector */
   init_grid_fcn(kd_, 0.0, &gf_);

/* tmp storage */
   copy_grid_fcn( kd_, gf_, &exact_ );

   printf("------------------------------\n");
   printf("Viscosity (nu): %e\n", kd_->nu_coeff);
   printf("Problem number (pnr): %i\n", kd_->pnr);
   printf("Boundary treatment: bcnr(left, right): %i, %i\n", bcnr(1), bcnr(2));
   printf("Wave number in exact solution: %i\n", wave_no);
   printf("Treatment of time-dependent bndry data: %i\n", kd_->taylorbc);
   printf("Solving to time %e using %i steps\n",kd_->tstop, kd_->nsteps);
   printf("Time step on fine grid is %e\n",kd_->dt_fine);
   printf("Grid spacing is %e with %i grid points\n", kd_->h_fine, kd_->n_fine);
   printf("Artificial damping coefficient: %e\n", kd_->ad_coeff);
   printf("Spatial order of accuracy: %i\n", kd_->spatial_order);
   printf("------------------------------\n");

   t = 0.0;

/* ! time stepping... */
   for (step=1; step<=kd_->nsteps; step++)
   {

/* ! reset initial conditions for the boundary data ode */
      if (step%bdatareset == 0)
      {
         printf("Resetting bdata initial cond. step = %i\n", step);
         bdata(gf_, t, kd_);
      }

/* ! evaluate solution error (stage=1 evaluates the plain bndry data at time t)*/
/* when stage=1, dt is not used */
      bdata( exact_, t, kd_);
      exact1( exact_, t, kd_ );

      evaldiff( gf_, exact_, &l2, &li );
      
/* ! save errors on file... */
      fprintf(eun,"%e %e %e %e\n", t, li, l2, fabs(ex_vsol(1) - vsol(1)) ); 

/* ! time-stepper from t_n to t_{n+1} starts here */
      explicit_rk4_stepper(kd_, t, t+kd_->dt_fine, dummy, gf_, &rfact_dummy); /* this is my_Phi() */

      t = t + kd_->dt_fine;
/* ! time-stepper from t_n to t_{n+1} ends here */
   }   

   printf("------------------------------\n");
   printf("Time-stepping completed. Solved to time t: %e\n", t);

/* ! evaluate solution error, stick exact solution in kd_ workspace array */
   exact1( exact_, t, kd_ );
/* get exact bndry data */
   bdata( exact_, t, kd_);

   evaldiff( gf_, exact_, &l2, &li );
/* ! save errors on file... */
   fprintf(eun,"%e %e %e %e\n", t, li, l2, fabs(ex_vsol(1)-vsol(1)));

/*! close error file*/
   fclose(eun);
   
   printf("------------------------------\n");
   
   printf("Solution error in maximum norm, bndry error\n");
   
   printf("time: %e, sol-err: %e, bndry-err: %e\n", t, li, fabs(ex_vsol(1)-vsol(1)));
   printf("------------------------------\n");

/* free tmp storage */
   free_grid_fcn(kd_, exact_);

}


             
