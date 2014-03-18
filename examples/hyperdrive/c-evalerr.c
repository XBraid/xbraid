#include <math.h>
#include "advect_data.h"

/* evaluate the L2 and Linf norm of the difference between grid functions w and we */
void
evaldiff( grid_fcn *w, grid_fcn *we, double *l2, double *li )
{
   int i;
   double locerr;
   int n = w->n;
   double h = w->h;
   
   *l2 = 0;
   *li = 0;
   
   for (i=1; i<n; i++)
   {
      locerr = fabs(w->sol[i]-we->sol[i]);
      *l2 += h*locerr*locerr;
      if( locerr > *li )
      {
         *li = locerr;
      }
   }
   
   *l2 = sqrt(*l2);
}

/* evaluate the L2 and Linf norm of grid function w */
void
evalnorm( grid_fcn *w, double *l2, double *li )
{
   int i;
   double locerr;
   int n = w->n;
   double h = w->h;
   
   *l2 = 0;
   *li = 0;
   
/* for periodic grid functions, only the first n-1 values are unique */
   for (i=1; i<n; i++)
   {
      locerr = fabs(w->sol[i]);
      *l2 += h*locerr*locerr;
      if( locerr > *li )
      {
         *li = locerr;
      }
   }
/* for periodic grid functions, only the first n-1 values are unique */
   *l2 = sqrt(*l2);
}

