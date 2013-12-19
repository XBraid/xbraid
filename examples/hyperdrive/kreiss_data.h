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

/* fcn prototypes */
void
init_kreiss_grid_fcn(kreiss_solver *kd_, double t, kreiss_grid_fcn *gf_);

void
init_kreiss_solver(double h, double amp, double ph, double om, int pnr, int taylorbc, 
                   double L, double cfl, kreiss_solver *kd_);

void
explicit_rk4_stepper(kreiss_solver *kd_, double t, double tend, double accuracy, kreiss_grid_fcn *gf_, 
                     int *rfact_);
void
exact1( int n, double *w, double h, double amp, double ph, double om, double t, int pnr);
void 
bdata(double_array_1d *vsol_, double amp, double ph, double om, double t, int pnr);
void 
bop6g(double t, double_array_2d *q06_ );
void 
sbpghost( int nb, int wb, double_array_2d * bop_, double *gh, double betapcoeff);
void
twbndry1( double x0, double *bdata0, double x1, double *bdata1, int s, double t, double dt, 
          double amp, double ph, double om, int pnr );
void
bckreiss1( int n, double *w, double bdataL, double bdataR, double betapcoeff, double h, int_array_1d *bcnr_ );
void
dwdtkreiss1( int n, double *w, double *dwdt, double h, int nb, int wb, double_array_2d *bop_, 
             double_array_2d *bope_, double gh);
void
twforce1( int n, double *f, double t, double h, double amp, double ph, double om, int pnr, double Lx );
void
dvdtbndry(double_array_1d *vsol_, double_array_1d *dvdt_, double amp, double ph, double om, double t, int pnr);
void
evalerr1( int n, double *w, double *we, double *l2, double*li, double h );
void
bckreiss1( int n, double *w, double bdataL, double bdataR, double betapcoeff, double h, int_array_1d *bcnr_ );
/* end propotypes */
   
#endif
