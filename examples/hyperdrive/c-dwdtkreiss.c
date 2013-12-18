#include "c_array.h"

/*    subroutine dwdtkreiss1( n, w, dwdt, h, nb, wb, bop, bope, gh)*/
void
dwdtkreiss1( int n, double *w, double *dwdt, double h, int nb, int wb, double_array_2d *bop_, 
             double_array_2d *bope_, double gh)
{
   const double d6a= 0.75, d6b=-0.15, d6c=1.0/60;
   
   int i, k;
   
/* bop: coefficients for standard SBP */
/* bope, gh: coefficients for SBP with ghost point */
/* bope(i,k) = bop(i,k) for i>=2 */
#define bop(i,j) compute_index_2d(bop_, i, j)
#define bope(i,j) compute_index_2d(bope_, i, j)
   double ih, i2h, du, dp, du1, dp1;
      
   i2h = 1.0/(2*h);
   ih = 1.0/(h);

/*! on the boundary */
   i=1;
   
   du  = 0;
   du1 = 0;
   
   for (k=1; k<=wb; k++)
   {
      du  = du  + bope(i,k) * w[k];
      du1 = du1 + bop(i,k)  * w[k];
   }
   
   du = du + gh*w[0];
   dwdt[i] = i2h*(-du - du1 );

/* near the left bndry */
   for (i=2; i<=nb; i++)
   {
      du  = 0;
      for (k=1; k<=wb; k++)
      {
         du  = du  + bope(i,k) * w[k];
      }
      dwdt[i] = ih*(-du );
   }

/* centered formula in the interior */
   for (i=nb+1; i<=n-nb; i++)
   {
      du  = d6c*(w[i+3]-w[i-3]) + d6b*(w[i+2]-w[i-2]) + d6a*(w[i+1]-w[i-1]);
      dwdt[i] = ih*(-du );
   }

/* near the right bndry */
   for (i=n-nb+1; i<=n-1; i++)
   {
      du   = 0;
      for ( k=n-wb+1; k<=n; k++)
      {
         du = du  - bope(n-i+1,n-k+1) * w[k];
      }
      dwdt[i] = ih*(-du );
   }

/*! on the right bndry */
   i=n;
   du = 0;
   du1 = 0;
   
   for (k=n-wb+1; k<=n; k++)
   {
      du  = du  - bope(n-i+1,n-k+1) * w[k];
      du1 = du1 - bop(n-i+1,n-k+1)  * w[k];
   }
   du = du - gh * w[n+1];
   dwdt[i] = i2h* (-du - du1 );
}

