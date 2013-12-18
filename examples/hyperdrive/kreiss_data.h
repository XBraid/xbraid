#ifndef KREISS_DATA
#define KREISS_DATA

#include "c_array.h"

typedef struct _warp_Vector_struct
{
   double *sol;
   int n;
   double_array_1d *vsol_;
} kreiss_grid_fcn;

typedef struct _warp_App_struct
{
   int n;        /* number of grid points */
   double h;     /* grid size */
   double amp;   /* testing parameter (amplitude) */
   double ph;    /* testing parameter (phase) */
   double om;    /* testing parameter (omega) */
   int pnr;      /* problem number (1 or 2) */
   int taylorbc; /* treatment of bndry data at intermediate RK stages */
   int nb; /* first dimension of bop and bope arrays */
   int wb; /* second dimension of bop and bope arrays */
   double_array_2d *bop_;  /* regular SBP coefficients */
   double_array_2d *bope_; /* extended SBP coefficients */
   double gh;              /* ghost point coefficient */
   double L;               /* length of 1-d domain */
   double dt;              /* time step satisfying cfl constraint */
   double betapcoeff;      /* coefficient for Dirichlet data with ghost point */
   int_array_1d *bcnr_;    /* boundary coefficient number */
   double_array_1d *alpha_, *beta_; /* RK-4 coefficients */
   double_array_1d *vcur_, *veval_, *dvdt_; /* boundary ode workspace variables */
   double *eval, *current, *rhs, *force;   /* grid function workspace variables */
} kreiss_solver;

/* propotypes */
void
init_kreiss_grid_fcn(kreiss_solver *kd_, double t, kreiss_grid_fcn *gf_);

void
init_kreiss_solver(double h, double amp, double ph, double om, int pnr, int taylorbc, 
                   double L, double cfl, kreiss_solver *kd_);

void
explicit_rk4_stepper(kreiss_solver *kd_, double t, double tend, double accuracy, kreiss_grid_fcn *gf_, 
                     int *rfact_);
   
#endif
