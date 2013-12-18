#include <math.h>

/*    subroutine twbndry1( x0, bdata0, x1, bdata1, s, t, dt, amp, ph, om, pnr, taylorbc )*/
void
twbndry1( double x0, double *bdata0, double x1, double *bdata1, int s, double t, double dt, 
          double amp, double ph, double om, int pnr )
{
   const double pi=M_PI;
   double dt2 = dt*dt, dt3 = dt*dt*dt;
   
   if( pnr == 1 )
   {
      if (s==1)
      {/*  g */
         *bdata0 = sin(pi*x0*x0+ph)*cos(t);
         *bdata1 = sin(pi*x1*x1+ph)*cos(t);
      }
      else if( s == 2 )
      { /* g + dt/2*g' */
         *bdata0 = sin(pi*x0*x0+ph)*cos(t)-dt*sin(pi*x0*x0+ph)*sin(t)/2;
         *bdata1 = sin(pi*x1*x1+ph)*cos(t)-dt*sin(pi*x1*x1+ph)*sin(t)/2;
      }
      else if( s == 3 )
      { /* g + dt/2*g' + dt*dt/4*g'' */
         *bdata0 = sin(pi*x0*x0+ph)*cos(t)-dt*sin(pi*x0*x0+ph)*sin(t)/2 - dt2*sin(pi*x0*x0+ph)*cos(t)/4;
         *bdata1 = sin(pi*x1*x1+ph)*cos(t)-dt*sin(pi*x1*x1+ph)*sin(t)/2 - dt2*sin(pi*x1*x1+ph)*cos(t)/4;
      }
      else if( s == 4 )
      { /* g + dt*g' + dt*dt/2*g'' + dt*dt*dt/4*g''' */
         *bdata0 = sin(pi*x0*x0+ph)*cos(t)-dt*sin(pi*x0*x0+ph)*sin(t) - dt2*sin(pi*x0*x0+ph)*cos(t)/2 +
            dt3*sin(pi*x0*x0+ph)*sin(t)/4;
         *bdata1 = sin(pi*x1*x1+ph)*cos(t)-dt*sin(pi*x1*x1+ph)*sin(t) - dt2*sin(pi*x1*x1+ph)*cos(t)/2 +
            dt3*sin(pi*x1*x1+ph)*sin(t)/4;
      }
   }
   else if ( pnr == 2 )
   {/*  w(1,i) = sin(om*(x0-t)+ph), w(2,i) = cos(om*(x0+t)) */
      if (s==1)
      {/*  g */
         *bdata0 = sin(om*(x0-t)+ph);
         *bdata1 = sin(om*(x1-t)+ph);
      }
      else if (s == 2)
      {/* g + dt/2*g' */
         *bdata0 = sin(om*(x0-t)+ph) - 0.5*dt*om*cos(om*(x0-t)+ph);
         *bdata1 = sin(om*(x1-t)+ph) - 0.5*dt*om*cos(om*(x1-t)+ph);
      }
      else if (s == 3)
      { /* g + dt/2*g' + dt*dt/4*g'' */
         *bdata0 = sin(om*(x0-t)+ph) - 0.5*dt*om*cos(om*(x0-t)+ph) - 0.25*dt*dt*om*om*sin(om*(x0-t)+ph);
         *bdata1 = sin(om*(x1-t)+ph) - 0.5*dt*om*cos(om*(x1-t)+ph) - 0.25*dt*dt*om*om*sin(om*(x1-t)+ph);
      }
      else if (s == 4)
      { /* g + dt*g' + dt*dt/2*g'' + dt*dt*dt/4*g''' */
         *bdata0 = sin(om*(x0-t)+ph) - dt*om*cos(om*(x0-t)+ph) - 0.5*dt*dt*om*om*sin(om*(x0-t)+ph) + 
            0.25*dt*dt*dt*om*om*om*cos(om*(x0-t)+ph);
         *bdata1 = sin(om*(x1-t)+ph) - dt*om*cos(om*(x1-t)+ph) - 0.5*dt*dt*om*om*sin(om*(x1-t)+ph) + 
            0.25*dt*dt*dt*om*om*om*cos(om*(x1-t)+ph);
      }
   }
}
