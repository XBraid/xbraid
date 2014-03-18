#include <math.h>
#include "c_array.h"
#include "advect_data.h"

void
dwdt( grid_fcn *w, grid_fcn *dwdt, double t, double bdata[2], advection_setup *kd_ )
{
   int i;
   int n = w->n;
   double h = w->h;
   double x;
   double pi = M_PI;
   grid_fcn *wxx = NULL;
/* artificial damping coefficient */
   double ad = 0.0;
   
/* add artificial damping if this is a coarse grid */
   if (h > 1.5*kd_->h_fine)
      ad += kd_->ad_coeff;

/* new stuff */
   assign_gp( w, bdata[0], bdata[1], kd_ );
/* end new */


/* first calculate dwdx and store it in dwdt grid_fcn */
   dwdx( w, dwdt, kd_ );

/* tmp storage */
   copy_grid_fcn( kd_, w, &wxx );
/* evaluate d2w/dx2 */
   d2wdx2( w, wxx, kd_ );

   for (i=1; i<=n; i++)
   {
      dwdt->sol[i] = - dwdt->sol[i] + (kd_->nu_coeff + ad*h)*wxx->sol[i];
         /* + ad/h*(w->sol[i-1] - 2.0*w->sol[i] + w->sol[i+1]); */
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

