/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder schroder2@llnl.gov, Rob Falgout falgout2@llnl.gov,
 * Tzanio Kolev kolev1@llnl.gov, Ulrike Yang yang11@llnl.gov, 
 * Veselin Dobrev dobrev1@llnl.gov, et al. 
 * LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. Email schroder2@llnl.gov on how to download. 
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (as published by the Free Software
 * Foundation) version 2.1 dated February 1999.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc., 59
 * Temple Place, Suite 330, Boston, MA 02111-1307 USA
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

