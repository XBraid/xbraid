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
dwdt( grid_fcn *w, grid_fcn *dwdt, double t, double bdata[2], advection_setup *kd_ )
{
   int i;
   int n = w->n;
   double h = w->h, ih;
   double x;
   double pi = M_PI;
   grid_fcn *wxx = NULL;

   ih = 1.0/h;
   
   assign_gp( w, bdata[0], bdata[1], kd_ );

/* first calculate dwdx and store it in dwdt grid_fcn */
   dwdx( w, dwdt, kd_ );
/* switch the sign (solving ut+ux=f) */
   for (i=1; i<=n; i++)
   {
      dwdt->sol[i] = - dwdt->sol[i];
   }

/* tmp storage */
   copy_grid_fcn( kd_, w, &wxx );

/* evaluate and add in the viscous term, but only if the viscocity is positive */
   if (kd_->nu_coeff > 0)
   {
/* evaluate d2w/dx2 */
      d2wdx2( w, wxx, kd_ );
      for (i=1; i<=n; i++)
      {
         dwdt->sol[i] += kd_->nu_coeff*wxx->sol[i];
      }
   }
   
/* evaluate and add in artificial dissipation term, but only if the coefficient is positive */
   if (kd_->ad_coeff > 0)
   {
/* evaluate undivided spatial difference */
      adterm( w, wxx, kd_ );
      for (i=1; i<=n; i++)
      {
         dwdt->sol[i] += kd_->ad_coeff*ih*wxx->sol[i];
      }
   }
   

/* add in forcing */
   if( kd_->pnr == 1 )
   {
      for (i=0; i<=n+1; i++)
      {
         x = (i-1)*h;
/* u-eqn */
         dwdt->sol[i] += -sin(pi*x*x+kd_->ph)*sin(t) + 2.0*pi*x*cos(pi*x*x+kd_->ph)*cos(t);
      }
   }

/* free tmp storage */
   free_grid_fcn( kd_, wxx );

}

