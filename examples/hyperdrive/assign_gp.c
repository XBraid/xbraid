#include "c_array.h"

void
assign_gp( int n, double *w, double bdataL, double bdataR, double betapcoeff, double h, int_array_1d *bcnr_ )
{
#define bcnr(i) compute_index_1d(bcnr_, i)   

/*! bcnr(1) = 1: Dirichlet in u, extrapolate v */
/*! bcnr(2) = 2: Dirichlet in v, extrapolate u */

/*!** Left side */
   if( bcnr(1) == 1 )
   {
/*!** Dirichlet in u */
      w[0] = (w[1]-bdataL)/(0.5*betapcoeff) + 4*w[1]-6*w[2]+4*w[3]-w[4];
   }
   else
   {
      printf("ERROR: assign_gp, unknown bcnr(1)= %i\n", bcnr(1));
   }

/*!** Right side */
   if( bcnr(2) == 2 )
   {
/*! extrapolate u */
      w[n+1] = 4*w[n]-6*w[n-1]+4*w[n-2]-w[n-3];
   }
   else
   {
      printf("ERROR: assign_gp, unknown bcnr(2)= %i\n", bcnr(2));
   }
}

