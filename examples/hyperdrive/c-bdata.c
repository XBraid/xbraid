/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. Email xbraid-support@llnl.gov for support.
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

