#include <math.h>
#include "advect_data.h"

void
evalerr1( grid_fcn *w, grid_fcn *we, double *l2, double *li )
{
   int i;
   double locerr;
   int n = w->n;
   double h = w->h;
   
   *l2 = 0;
   *li = 0;
   
   for (i=1; i<=n; i++)
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

