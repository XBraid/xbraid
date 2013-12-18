#include <math.h>

/*    subroutine evalerr1( n, w, we, l2, li, h )*/
void
evalerr1( int n, double *w, double *we, double *l2, double*li, double h )
{
   int i;
   double locerr;
   
   *l2 = 0;
   *li = 0;
   
   for (i=1; i<=n; i++)
   {
      locerr = fabs(w[i]-we[i]);
      *l2 += h*locerr*locerr;
      if( locerr > *li )
      {
         *li = locerr;
      }
   }
   
   *l2 = sqrt(*l2);
}

