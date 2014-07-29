/*BHEADER**********************************************************************
 * Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of WARP.  See file COPYRIGHT for details.
 *
 * WARP is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 ***********************************************************************EHEADER*/

#include <math.h>
#include "c_array.h"
#include "advect_data.h"

void 
bdata( grid_fcn *w, double t, advection_setup *kd_)
{
   double pi, x0;
   
#define vsol(i) compute_index_1d(w->vsol_, i)   

   pi=M_PI;
   
   x0=0.0;
   
   if( kd_->pnr == 1 )
   {
/*! w(x,t) = sin(pi*x0**2+ph)*cos(t)*/
      vsol(1) = sin(pi*x0*x0+kd_->ph)*cos(t);
      vsol(2) = -sin(pi*x0*x0+kd_->ph)*sin(t);
      vsol(3) = -sin(pi*x0*x0+kd_->ph)*cos(t);
   }
   else if ( kd_->pnr == 2 ) /* AP: Need to update exact solution wrt viscosity */
   {
/*! w(x,t) = sin(om*(x0-t)+ph) */
      vsol(1) = sin(kd_->om*(x0-t)+kd_->ph);
      vsol(2) = -kd_->om*cos(kd_->om*(x0-t)+kd_->ph);
      vsol(3) = -kd_->om*kd_->om*sin(kd_->om*(x0-t)+kd_->ph);
   }
#undef vsol
}


void
dvdtbndry(grid_fcn *w, grid_fcn *dwdt, double t, advection_setup *kd_)
{
   double pi, x0, g;
#define vsol(i) compute_index_1d(w->vsol_, i)   
#define dvdt(i) compute_index_1d(dwdt->vsol_, i)
   pi=M_PI;
   x0=0.0;
   
/*! forcing*/
   if( kd_->pnr == 1 ) 
   {
      g = sin(pi*x0*x0+kd_->ph)*sin(t);
   }
   else if ( kd_->pnr == 2 )
   {
      g = kd_->om*kd_->om*kd_->om*cos(kd_->om*(x0-t)+kd_->ph);
   }
   dvdt(1) = vsol(2);
   dvdt(2) = vsol(3);
   dvdt(3) = g;
}

