#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "c_array.h"
#include "kreiss_data.h"

/* fcn prototypes */
void 
bdata(double_array_1d *vsol_, double amp, double ph, double om, double t, int pnr);
void
dvdtbndry(double_array_1d *vsol_, double_array_1d *dvdt_, double amp, double ph, double om, double t, int pnr);
void
twbndry1( double x0, double *bdata0, double x1, double *bdata1, int s, double t, double dt, 
          double amp, double ph, double om, int pnr );
void
exact1( int n, double *w, double h, double amp, double ph, double om, double t, int pnr);
void
evalerr1( int n, double *w, double *we, double *l2, double*li, double h );
void
bckreiss1( int n, double *w, double bdataL, double bdataR, double betapcoeff, double h, int_array_1d *bcnr_ );
void
dwdtkreiss1( int n, double *w, double *dwdt, double h, int nb, int wb, double_array_2d *bop_, 
             double_array_2d *bope_, double gh);
void
twforce1( int n, double *f, double t, double h, double amp, double ph, double om, int pnr, double Lx );
/* end propotypes */

int main(int argc, char ** argv)
{
   int step, nsteps, pnr, taylorbc;
   const int bdatareset=1000;
   
   double h, cfl, bdataL, bdataR;
   double L, t, l2, li, tfinal;
   double amp, ph, om;
   double dummy=0;
   int rfact_dummy=0;
   
   int nstepsset, tfinalset, arg_index, print_usage=0, myid=0;

   FILE *eun;
   
   kreiss_solver *kd_ = NULL;
   grid_fcn *gf_ = NULL;
      
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
      printf("Solve Kreiss equation with a SBP finite difference method and 4th order explicit RK:\n");
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
      
/* open file for saving solution error data */
   eun = fopen("err.dat","w");

#define bcnr(i) compute_index_1d(kd_->bcnr_, i)    
#define vsol(i) compute_index_1d(gf_->vsol_, i)   

/* solver meta-data */
   kd_ = malloc(sizeof(kreiss_solver));
   init_kreiss_solver(h, amp, ph, om, pnr, taylorbc, L, cfl, nstepsset, nsteps, tfinal, kd_);
   
/* create solution vector */
   init_grid_fcn(kd_, 0.0, &gf_);

   printf("------------------------------\n");
   printf("Problem number (pnr): %i\n", kd_->pnr);
   printf("Boundary treatment: bcnr(left, right): %i, %i\n", bcnr(1), bcnr(2));
   printf("Treatment of time-dependent bndry data: %i\n", kd_->taylorbc);
   printf("Solving to time %e using %i steps\n",kd_->tstop, kd_->nsteps);
   printf("Time step is %e\n",kd_->dt);
   printf("Grid spacing is %e with %i grid points\n", kd_->h, kd_->n);

   t = 0.0;

/* ! time stepping... */
   for (step=1; step<=kd_->nsteps; step++)
   {

/* ! reset initial conditions for the boudary data ode */
      if (step%bdatareset == 0)
      {
         printf("Resetting bdata initial cond. step = %i\n", step);
         bdata(gf_->vsol_, kd_->amp, kd_->ph, kd_->om, t, kd_->pnr);
      }

/* ! evaluate solution error (stage=1 evaluates the plain bndry data at time t)*/
      twbndry1( 0.0, &bdataL, kd_->L, &bdataR, 1, t, kd_->dt, kd_->amp, kd_->ph, kd_->om, kd_->pnr);
      
      exact1( kd_->n, kd_->current, kd_->h, kd_->amp, kd_->ph, kd_->om, t, kd_->pnr );
      evalerr1( kd_->n, gf_->sol, kd_->current, &l2, &li, kd_->h );
      
/* ! save errors on file... */
      fprintf(eun,"%e %e %e %e\n", t, li, l2, fabs(bdataL - vsol(1)) ); 

/* ! time-stepper from t_n to t_{n+1} starts here */
      explicit_rk4_stepper(kd_, t, t+kd_->dt, dummy, gf_, &rfact_dummy); /* this is my_Phi() */

      t = t + kd_->dt;
/* ! time-stepper from t_n to t_{n+1} ends here */
   }   

   printf("------------------------------\n");
   printf("Time-stepping completed. Solved to time t: %e\n", t);

/* ! evaluate solution error, stick exact solution in kd_ workspace array */
   exact1( kd_->n, kd_->current, kd_->h, kd_->amp, kd_->ph, kd_->om, t, kd_->pnr );
/* get exact bndry data (bdataL) */
   twbndry1( 0.0, &bdataL, kd_->L, &bdataR, 1, t, kd_->dt, kd_->amp, kd_->ph, kd_->om, kd_->pnr);

   evalerr1( kd_->n, gf_->sol, kd_->current, &l2, &li, kd_->h );
/* ! save errors on file... */
   fprintf(eun,"%e %e %e %e\n", t, li, l2, fabs(bdataL-vsol(1)));

/*! close error file*/
   fclose(eun);
   

   printf("------------------------------\n");
   
   printf("Solution error in maximum norm, bndry error\n");
   
   printf("time: %e, sol-err: %e, bndry-err: %e\n", t, li, fabs(bdataL-vsol(1)));
   printf("------------------------------\n");
/*   printf("Saving ...\n");*/
    /* open(21,file='sol.bin',form='unformatted') */
    /* fprintf(21) n, h, dt, t, (sol(i),i=1,n), (current(i),i=1,n) */
    /* close(21) */
/*   printf("done.\n");*/
/*   printf("------------------------------\n");*/
}


             
