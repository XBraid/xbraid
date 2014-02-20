#ifndef ADVECT_DATA
#define ADVECT_DATA

#include "c_array.h"
#include "warp.h"

#define MY_EPS 1e-12

/* define HD_DEBUG to get printouts from various user defined functions*/
/* #define HD_DEBUG */

typedef struct _warp_Vector_struct
{
   double *sol;
   int n;
   double h;     /* grid size */
   double_array_1d *vsol_;
} grid_fcn;

typedef struct _warp_App_struct
{
   int n_fine;     /* number of grid points in the finest grid */
   double h_fine;  /* finest grid size */
   double dt_fine; /* time step satisfying cfl constraint on finest grid */
   double amp;     /* testing parameter (amplitude) */
   double ph;      /* testing parameter (phase) */
   double om;      /* testing parameter (omega) */
   int pnr;        /* problem number (1 or 2) */
   int taylorbc;   /* treatment of bndry data at intermediate RK stages */
   int nb;         /* first dimension of bop and bope arrays */
   int wb;         /* second dimension of bop and bope arrays */
   double_array_2d *bop_;  /* regular SBP coefficients */
   double_array_2d *bope_; /* extended SBP coefficients */
   double gh;              /* ghost point coefficient */
   int nb2;         /* first dimension of bop2 array */
   int wb2;         /* second dimension of bop2 array */
   double_array_2d *bop2_;  /* SBP coefficients for 2nd derivative */
   double_array_1d *iop2_;  /* interior coefficients for 2nd derivative */
   double gh2;              /* ghost point coefficient for 2nd derivative */
   double bder[5];          /* coefficients for boundary derivative */
   double L;               /* length of 1-d domain */
   double c_coeff;         /* wave speed */
   double nu_coeff;        /* viscosity */
   double betapcoeff;      /* coefficient for Dirichlet data with ghost point */
   int_array_1d *bcnr_;    /* boundary coefficient number */
   double_array_1d *alpha_, *beta_; /* RK-4 coefficients */

/* warp specific stuff */
   grid_fcn *sol_copy; /* assigned by the call-back routine save_grid_fcn() */
   double t_copy;
   
   int write; /* flag to tell warp if it should save grid functions to file */
   double tstart, tstop;
   int nsteps;
   
} advection_setup;

/* fcn prototypes */
int
init_grid_fcn(advection_setup *kd_, double t, grid_fcn **u_handle);

void
init_advection_solver(double h, double amp, double ph, double om, int pnr, int taylorbc, 
                      double L, double cfl, int nstepsset, int nsteps, double tfinal, 
                      double wave_speed, double viscosity, advection_setup *kd_);

int
explicit_rk4_stepper(advection_setup *kd_, double t, double tend, double accuracy, grid_fcn *gf_, 
                     int *rfact_);
int
copy_grid_fcn(advection_setup    *kd_,
              grid_fcn  *u_,
              grid_fcn **v_handle);
int
free_grid_fcn(advection_setup    *kd_,
              grid_fcn  *u_);
int
sum_grid_fcn(advection_setup *kd_,
             double      alpha,
             grid_fcn *x_,
             double      beta,
             grid_fcn *y_);
int
dot_grid_fcn(advection_setup *kd_,
             grid_fcn *u_,
             grid_fcn *v_,
             double      *dot_ptr);
int
save_grid_fcn(advection_setup *kd_,
              double t,
              warp_Status   status,
              grid_fcn *u_);
int
gridfcn_BufSize(advection_setup *kd_,
                int *size_ptr);
int
gridfcn_BufPack(advection_setup *kd_,
                grid_fcn *u_,
                void *buffer);
int
gridfcn_BufUnpack(advection_setup *kd_,
                  void *buffer,
                  warp_Vector *u_handle);
int
gridfcn_Refine(advection_setup * kd_,
               double     tstart,
               double     f_tminus,
               double     f_tplus,
               double     c_tminus,
               double     c_tplus,
               grid_fcn *cu_,
               grid_fcn **fu_handle);
int
gridfcn_Coarsen(advection_setup *kd_,
                double tstart,
                double f_tminus,
                double f_tplus,
                double c_tminus,
                double c_tplus,
                grid_fcn *fu_,
                grid_fcn **cu_handle);

void
exact1( double *w, double t, advection_setup *kd_);
void
exact_x( double *w, double t, advection_setup *kd_);
void
exact_xx( double *w, double t, advection_setup *kd_);

void 
bdata(double_array_1d *vsol_, double amp, double ph, double om, double t, int pnr);
void 
bop6g(double t, double_array_2d *q06_ );
void 
diffusion_coeff_4( double_array_1d *iop2_, double_array_2d *bop2_, double *gh2, double bder[5] );
void 
sbpghost( int nb, int wb, double_array_2d * bop_, double *gh, double betapcoeff);
void
twbndry1( double *bdataL, double *bdataR, int stage, double t, double dt, advection_setup *kd_ );
void
assign_gp( int n, double *w, double bdataL, double bdataR, double betapcoeff, double h, int_array_1d *bcnr_ );
void
dwdt( int n, double *w, double *dwdt, double h, advection_setup *kd_ );
void
dwdx( int n, double *w, double *dwdt, double h, advection_setup *kd_ );
void
d2wdx2( int n, double *w, double *wxx, double h, advection_setup *kd_ );
void
twforce1( int n, double *f, double t, double h, double amp, double ph, double om, int pnr, double Lx );
void
dvdtbndry(double_array_1d *vsol_, double_array_1d *dvdt_, double amp, double ph, double om, double t, int pnr);
void
evalerr1( int n, double *w, double *we, double *l2, double*li, double h );

/* end propotypes */
   
#endif
