#include <math.h>
#include "advect_data.h"

#define SQR(x) ((x)*(x))
void
exact1( double *w, double t, advection_setup *kd_)
{
   int i;
   double x;
   const double pi=M_PI;

   if (kd_->pnr == 1){
      for (i=0; i<=kd_->n_fine+1; i++){
         x = (i-1)*kd_->h_fine;
         w[i] = sin(pi*x*x+kd_->ph)*cos(t);
      }
   }
   else if (kd_->pnr == 2)
   {
      for (i=0; i<=kd_->n_fine+1; i++){
         x = (i-1)*kd_->h_fine;
         w[i] = sin(kd_->om*(x-t)+kd_->ph);
      }
   }
}

void
exact_x( double *w, double t, advection_setup *kd_)
{
   int i;
   double x;
   const double pi=M_PI;

   if (kd_->pnr == 1){
      for (i=0; i<=kd_->n_fine+1; i++){
         x = (i-1)*kd_->h_fine;
         w[i] = 2*x*pi*cos(pi*x*x+kd_->ph)*cos(t);
      }
   }
   else if (kd_->pnr == 2)
   {
      for (i=0; i<=kd_->n_fine+1; i++){
         x = (i-1)*kd_->h_fine;
         w[i] = kd_->om*cos(kd_->om*(x-t)+kd_->ph);
      }
   }
}

void
exact_xx( double *w, double t, advection_setup *kd_)
{
   int i;
   double x;
   const double pi=M_PI;

   if (kd_->pnr == 1){
      for (i=0; i<=kd_->n_fine+1; i++){
         x = (i-1)*kd_->h_fine;
         w[i] = 2*x*pi*cos(pi*x*x+kd_->ph)*cos(t); /* update! */
      }
   }
   else if (kd_->pnr == 2)
   {
      for (i=0; i<=kd_->n_fine+1; i++){
         x = (i-1)*kd_->h_fine;
         w[i] = -SQR(kd_->om)*sin(kd_->om*(x-t)+kd_->ph);
      }
   }
}

