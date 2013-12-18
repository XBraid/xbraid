#include <stdlib.h>
#include <math.h>
#include "kreiss_data.h"

void
init_kreiss_grid_fcn(kreiss_solver *kd_, double t, kreiss_grid_fcn *gf_)
{
/* should make sure the kd_ and gf_ are not NULL */
   gf_->n     = kd_->n;
   gf_->sol   = malloc((kd_->n+2)*sizeof(double));
   gf_->vsol_ = create_double_array_1d(3);

/* initial conditions */
   exact1( kd_->n, gf_->sol, kd_->h, kd_->amp, kd_->ph, kd_->om, t, kd_->pnr );
   bdata(gf_->vsol_, kd_->amp, kd_->ph, kd_->om, t, kd_->pnr);

}

void
init_kreiss_solver(double h, double amp, double ph, double om, int pnr, int taylorbc, 
                   double L, double cfl, kreiss_solver *kd_)
{
   double mxeg;
   int n;
/* # ghost points */
   const int o = 1;

   n = (int) L/h+o;
   
   if( fabs( L/(n-o) - h ) > 1e-10 )
   {
      h  = L/(n-o);
      printf("NOTE: dx changed to %e\n", h);
   }

   kd_->n = n;
   kd_->h = h;
   kd_->amp = amp;
   kd_->ph = ph;
   kd_->om = om;
   kd_->pnr = pnr;
   kd_->taylorbc = taylorbc;
   kd_->L = L;
/* ! compute time step */
   mxeg = 1.0;
   kd_->dt = cfl*h/mxeg;
   
/* difference operators are currently hard-wired for 6/3 order SBP */
   kd_->nb = 6;
   kd_->wb = 9;
   kd_->bop_ = create_double_array_2d(kd_->nb,  kd_->wb);
   kd_->bope_ = create_double_array_2d(kd_->nb, kd_->wb);
/* sbp coefficients */
   bop6g( 0.7037, kd_->bop_ );
   bop6g( 0.7037, kd_->bope_ );
/* !** Extra part of SBP operator, for ghost point SBP. */
   kd_->betapcoeff = 0.03;
/* ! bope and gh contain weights for the modified SBP-operator \tilde{D} */
   sbpghost(kd_->nb, kd_->wb, kd_->bope_, &(kd_->gh), kd_->betapcoeff);
/*!** Boundary conditions on the two boundaries */
   kd_->bcnr_ = create_int_array_1d(2);
/*   code means: 1- Dirichlet in u, extrapolate v */
/*               2- Dirichlet in v, extrapolate u*/
#define bcnr(i) compute_index_1d(kd_->bcnr_, i)   
   if( kd_->pnr == 1 || kd_->pnr == 2 )
   {
      bcnr(1) = 1;
      bcnr(2) = 2;
   }

/* bndry ode workspace variables */
   kd_->vcur_  = create_double_array_1d(3);
   kd_->veval_ = create_double_array_1d(3);
   kd_->dvdt_  = create_double_array_1d(3);

/*!** RK coefficients */
   kd_->alpha_ = create_double_array_1d(4);
   kd_->beta_  = create_double_array_1d(4);

#define alpha(i) compute_index_1d(kd_->alpha_, i)   
#define beta(i) compute_index_1d(kd_->beta_, i)   

   alpha(1) = 0;
   alpha(2) = 0.5;
   alpha(3) = 0.5;
   alpha(4) = 1.0;
   
   beta(1) = 1.0/6.0;
   beta(2) = 1.0/3.0;
   beta(3) = 1.0/3.0;
   beta(4) = 1.0/6.0;

/* grid function workspace variables */
   kd_->eval = malloc((n+2)*sizeof(double));
   kd_->current = malloc((n+2)*sizeof(double));
   kd_->rhs = malloc((n+2)*sizeof(double));
   kd_->force = malloc((n+2)*sizeof(double));
   
}
