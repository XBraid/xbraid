#include <stdio.h>
#include <math.h>

/*    subroutine twforce1( n, f, t, h, amp, ph, om, pnr, taylorbc )*/
void
twforce1( int n, double *f, double t, double h, double amp, double ph, double om, int pnr, double Lx )
{
   int i;
   const double pi=M_PI;
   double x;
   
   if( pnr == 1 )
   {
      for (i=0; i<=n+1; i++)
      {
         x = (i-1)*h;
/* u-eqn */
         f[i] = -sin(pi*x*x+ph)*sin(t) + 2.0*pi*x*cos(pi*x*x+ph)*cos(t);
      }
   }
   else if (pnr == 2)
   {
      for (i=0; i<=n+1; i++)
      {
         f[i] = 0;
      }
   }
   else
   {
      printf("ERROR: twforce, unknown pnr = %i\n", pnr);
   }
}

