#include <math.h>
#include "c_array.h"

/*subroutine bdata(vsol, amp, ph, om, t, pnr )*/
void 
bdata(double_array_1d *vsol_, double amp, double ph, double om, double t, int pnr)
{
   double pi, x0;
   
#define vsol(i) compute_index_1d(vsol_, i)   

   pi=M_PI;
   
   x0=0.0;
   
   if( pnr == 1 )
   {
/*! w(x,t) = sin(pi*x0**2+ph)*cos(t)*/
      vsol(1) = sin(pi*x0*x0+ph)*cos(t);
      vsol(2) = -sin(pi*x0*x0+ph)*sin(t);
      vsol(3) = -sin(pi*x0*x0+ph)*cos(t);
   }
   else if ( pnr == 2 )
   {
/*! w(x,t) = sin(om*(x0-t)+ph) */
      vsol(1) = sin(om*(x0-t)+ph);
      vsol(2) = -om*cos(om*(x0-t)+ph);
      vsol(3) = -om*om*sin(om*(x0-t)+ph);
   }
#undef vsol
}


/*subroutine dvdtbndry(vsol, dvdt, amp, ph, om, t, pnr )*/
void
dvdtbndry(double_array_1d *vsol_, double_array_1d *dvdt_, double amp, double ph, double om, double t, int pnr)
{
   double pi, x0, g;
   const double eps=0.0;
#define vsol(i) compute_index_1d(vsol_, i)   
#define dvdt(i) compute_index_1d(dvdt_, i)
   pi=M_PI;
   x0=0.0;
   
/*! forcing*/
   if( pnr == 1 )
   {
      g = sin(pi*x0*x0+ph)*sin(t);
   }
   else if ( pnr == 2 )
   {
      g = om*om*om*cos(om*(x0-t)+ph);
   }
   dvdt(1) = (1.0+eps)*vsol(2);
   dvdt(2) = (1.0+eps)*vsol(3) + eps*vsol(1);
   dvdt(3) = eps*vsol(2) + g;
}

