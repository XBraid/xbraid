/*BHEADER**********************************************************************
 * Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of XBraid.  See file COPYRIGHT for details.
 *
 * XBraid is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 ***********************************************************************EHEADER*/

#include <math.h>
#include "advect_data.h"

#define SQR(x) ((x)*(x))
void
exact1( grid_fcn *w, double t, advection_setup *kd_)
{
   int i;
   double x;
   int n=w->n;
   double h=w->h;
   const double pi=M_PI;

   if (kd_->pnr == 1){
      for (i=0; i<=n+1; i++){
         x = (i-1)*h;
         w->sol[i] = sin(pi*x*x+kd_->ph)*cos(t);
      }
   }
   else if (kd_->pnr == 2)
   {
      for (i=0; i<=n+1; i++){
         x = (i-1)*h;
         w->sol[i] = exp(-(kd_->nu_coeff)*SQR(kd_->om)*t)*sin(kd_->om*(x-t)+kd_->ph);
      }
   }
}

void
exact_t( grid_fcn *w, double t, advection_setup *kd_)
{
   int i;
   double x;
   int n=w->n;
   double h=w->h;

   const double pi=M_PI;

   if (kd_->pnr == 1){ /* needs to be updated */
      for (i=0; i<=n+1; i++){
         x = (i-1)*h;
         w->sol[i] = 2*x*pi*cos(pi*x*x+kd_->ph)*cos(t);
      }
   }
   else if (kd_->pnr == 2)
   {
      for (i=0; i<=n+1; i++){
         x = (i-1)*h;
         w->sol[i] = -kd_->om*exp(-(kd_->nu_coeff)*SQR(kd_->om)*t)*cos(kd_->om*(x-t)+kd_->ph) -
            kd_->nu_coeff*SQR(kd_->om)*exp(-(kd_->nu_coeff)*SQR(kd_->om)*t)*sin(kd_->om*(x-t)+kd_->ph);
      }
   }
}

void
exact_x( grid_fcn *w, double t, advection_setup *kd_)
{
   int i;
   double x;
   int n=w->n;
   double h=w->h;

   const double pi=M_PI;

   if (kd_->pnr == 1){
      for (i=0; i<=n+1; i++){
         x = (i-1)*h;
         w->sol[i] = 2*x*pi*cos(pi*x*x+kd_->ph)*cos(t);
      }
   }
   else if (kd_->pnr == 2)
   {
      for (i=0; i<=n+1; i++){
         x = (i-1)*h;
         w->sol[i] = kd_->om*exp(-(kd_->nu_coeff)*SQR(kd_->om)*t)*cos(kd_->om*(x-t)+kd_->ph);
      }
   }
}

void
exact_xx( grid_fcn *w, double t, advection_setup *kd_)
{
   int i;
   double x;
   int n=w->n;
   double h=w->h;
   
   const double pi=M_PI;

   if (kd_->pnr == 1){
      for (i=0; i<=n+1; i++){
         x = (i-1)*h;
         w->sol[i] = 2*x*pi*cos(pi*x*x+kd_->ph)*cos(t); /* update! */
      }
   }
   else if (kd_->pnr == 2)
   {
      for (i=0; i<=n+1; i++){
         x = (i-1)*h;
         w->sol[i] = -SQR(kd_->om)*exp(-(kd_->nu_coeff)*SQR(kd_->om)*t)*sin(kd_->om*(x-t)+kd_->ph);
      }
   }
}

