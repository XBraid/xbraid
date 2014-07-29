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

#include <stdio.h>
#include <math.h>
#include "advect_data.h"

void
twforce1( grid_fcn *f, double t, advection_setup *kd_ )
{
   int i;
   const double pi=M_PI;
   double x;
   int n = f->n;
   double h = f->h;
   
   if( kd_->pnr == 1 )
   {
      for (i=0; i<=n+1; i++)
      {
         x = (i-1)*h;
/* u-eqn */
         f->sol[i] = -sin(pi*x*x+kd_->ph)*sin(t) + 2.0*pi*x*cos(pi*x*x+kd_->ph)*cos(t);
      }
   }
   else if (kd_->pnr == 2)
   {
      for (i=0; i<=n+1; i++)
      {
         f->sol[i] = 0;
      }
   }
   else
   {
      printf("ERROR: twforce, unknown pnr = %i\n", kd_->pnr);
   }
}

