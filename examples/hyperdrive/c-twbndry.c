#include <math.h>
#include "advect_data.h"
void
twbndry1( double *bdata0, double *bdata1, int s, double t, double dt, advection_setup *kd_ )
{
   const double pi=M_PI;
   double dt2 = dt*dt, dt3 = dt*dt*dt;
   double x0 = 0;
   double x1 = kd_->L;
#define bcnr(i) compute_index_1d(kd_->bcnr_, i)   
   
   if (bcnr(1)==0)
   { /* periodic case does not need any boundary data */
      *bdata0 = 0.0;
      *bdata1 = 0.0;      
   }
   else
   { /* non-periodic case */
      if( kd_->pnr == 1 )
      {
         if (s==1)
         {/*  g */
            *bdata0 = sin(pi*x0*x0+kd_->ph)*cos(t);
            *bdata1 = sin(pi*x1*x1+kd_->ph)*cos(t);
         }
         else if( s == 2 )
         { /* g + dt/2*g' */
            *bdata0 = sin(pi*x0*x0+kd_->ph)*cos(t)-dt*sin(pi*x0*x0+kd_->ph)*sin(t)/2;
            *bdata1 = sin(pi*x1*x1+kd_->ph)*cos(t)-dt*sin(pi*x1*x1+kd_->ph)*sin(t)/2;
         }
         else if( s == 3 )
         { /* g + dt/2*g' + dt*dt/4*g'' */
            *bdata0 = sin(pi*x0*x0+kd_->ph)*cos(t)-dt*sin(pi*x0*x0+kd_->ph)*sin(t)/2 - dt2*sin(pi*x0*x0+kd_->ph)*cos(t)/4;
            *bdata1 = sin(pi*x1*x1+kd_->ph)*cos(t)-dt*sin(pi*x1*x1+kd_->ph)*sin(t)/2 - dt2*sin(pi*x1*x1+kd_->ph)*cos(t)/4;
         }
         else if( s == 4 )
         { /* g + dt*g' + dt*dt/2*g'' + dt*dt*dt/4*g''' */
            *bdata0 = sin(pi*x0*x0+kd_->ph)*cos(t)-dt*sin(pi*x0*x0+kd_->ph)*sin(t) - dt2*sin(pi*x0*x0+kd_->ph)*cos(t)/2 +
               dt3*sin(pi*x0*x0+kd_->ph)*sin(t)/4;
            *bdata1 = sin(pi*x1*x1+kd_->ph)*cos(t)-dt*sin(pi*x1*x1+kd_->ph)*sin(t) - dt2*sin(pi*x1*x1+kd_->ph)*cos(t)/2 +
               dt3*sin(pi*x1*x1+kd_->ph)*sin(t)/4;
         }
      }
      else if ( kd_->pnr == 2 )
      {/*  w(1,i) = sin(om*(x0-t)+ph), w(2,i) = cos(om*(x0+t)) */
         if (s==1)
         {/*  g */
            *bdata0 = sin(kd_->om*(x0-t)+kd_->ph);
            *bdata1 = sin(kd_->om*(x1-t)+kd_->ph);
         }
         else if (s == 2)
         {/* g + dt/2*g' */
            *bdata0 = sin(kd_->om*(x0-t)+kd_->ph) - 0.5*dt*kd_->om*cos(kd_->om*(x0-t)+kd_->ph);
            *bdata1 = sin(kd_->om*(x1-t)+kd_->ph) - 0.5*dt*kd_->om*cos(kd_->om*(x1-t)+kd_->ph);
         }
         else if (s == 3)
         { /* g + dt/2*g' + dt*dt/4*g'' */
            *bdata0 = sin(kd_->om*(x0-t)+kd_->ph) - 0.5*dt*kd_->om*cos(kd_->om*(x0-t)+kd_->ph) - 0.25*dt*dt*kd_->om*kd_->om*sin(kd_->om*(x0-t)+kd_->ph);
            *bdata1 = sin(kd_->om*(x1-t)+kd_->ph) - 0.5*dt*kd_->om*cos(kd_->om*(x1-t)+kd_->ph) - 0.25*dt*dt*kd_->om*kd_->om*sin(kd_->om*(x1-t)+kd_->ph);
         }
         else if (s == 4)
         { /* g + dt*g' + dt*dt/2*g'' + dt*dt*dt/4*g''' */
            *bdata0 = sin(kd_->om*(x0-t)+kd_->ph) - dt*kd_->om*cos(kd_->om*(x0-t)+kd_->ph) - 0.5*dt*dt*kd_->om*kd_->om*sin(kd_->om*(x0-t)+kd_->ph) + 
               0.25*dt*dt*dt*kd_->om*kd_->om*kd_->om*cos(kd_->om*(x0-t)+kd_->ph);
            *bdata1 = sin(kd_->om*(x1-t)+kd_->ph) - dt*kd_->om*cos(kd_->om*(x1-t)+kd_->ph) - 0.5*dt*dt*kd_->om*kd_->om*sin(kd_->om*(x1-t)+kd_->ph) + 
               0.25*dt*dt*dt*kd_->om*kd_->om*kd_->om*cos(kd_->om*(x1-t)+kd_->ph);
         }
      }
   }
#undef bcnr
}
