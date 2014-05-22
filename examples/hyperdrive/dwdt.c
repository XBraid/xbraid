#include <math.h>
#include "c_array.h"
#include "advect_data.h"

void
dwdt( grid_fcn *w, grid_fcn *dwdt, double t, double bdata[2], advection_setup *kd_ )
{
   int i;
   int n = w->n;
   double h = w->h, ih;
   double x;
   double pi = M_PI;
   grid_fcn *wxx = NULL;

   ih = 1.0/h;
   
   assign_gp( w, bdata[0], bdata[1], kd_ );

/* first calculate dwdx and store it in dwdt grid_fcn */
   dwdx( w, dwdt, kd_ );
/* switch the sign (solving ut+ux=f) */
   for (i=1; i<=n; i++)
   {
      dwdt->sol[i] = - dwdt->sol[i];
   }

/* tmp storage */
   copy_grid_fcn( kd_, w, &wxx );

/* evaluate and add in the viscous term, but only if the viscocity is positive */
   if (kd_->nu_coeff > 0)
   {
/* evaluate d2w/dx2 */
      d2wdx2( w, wxx, kd_ );
      for (i=1; i<=n; i++)
      {
         dwdt->sol[i] += kd_->nu_coeff*wxx->sol[i];
      }
   }
   
/* evaluate and add in artificial dissipation term, but only if the coefficient is positive */
   if (kd_->ad_coeff > 0)
   {
/* evaluate undivided spatial difference */
      adterm( w, wxx, kd_ );
      for (i=1; i<=n; i++)
      {
         dwdt->sol[i] += kd_->ad_coeff*ih*wxx->sol[i];
      }
   }
   

/* add in forcing */
   if( kd_->pnr == 1 )
   {
      for (i=0; i<=n+1; i++)
      {
         x = (i-1)*h;
/* u-eqn */
         dwdt->sol[i] += -sin(pi*x*x+kd_->ph)*sin(t) + 2.0*pi*x*cos(pi*x*x+kd_->ph)*cos(t);
      }
   }

/* free tmp storage */
   free_grid_fcn( kd_, wxx );

}

