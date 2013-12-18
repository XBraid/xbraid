/*    subroutine exact1( n, w, h, amp, ph, om, t, pnr ) */
#include <math.h>
void
exact1( int n, double *w, double h, double amp, double ph, double om, double t, int pnr)
{
   int i;
   double x;
   const double pi=M_PI;

   if (pnr == 1){
      for (i=0; i<=n+1; i++){
         x = (i-1)*h;
         w[i] = sin(pi*x*x+ph)*cos(t);
      }
   }
   else if (pnr == 2)
   {
      for (i=0; i<=n+1; i++){
         x = (i-1)*h;
         w[i] = sin(om*(x-t)+ph);
      }
   }
}

