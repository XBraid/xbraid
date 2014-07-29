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

#include "c_array.h"
#include "advect_data.h"

void
assign_gp( grid_fcn *w, double bdataL, double bdataR, advection_setup *kd_ )
{
   int n = w->n;
   
#define bcnr(i) compute_index_1d(kd_->bcnr_, i)   

/*! bcnr(1) = 1: Dirichlet in u, extrapolate v */
/*! bcnr(2) = 2: Dirichlet in v, extrapolate u */

/*!** Left side */
   if( bcnr(1) == 1 )
   {
/*!** Dirichlet in u */
      w->sol[0] = (w->sol[1]-bdataL)/(0.5*kd_->betapcoeff) + 4*w->sol[1]-6*w->sol[2]+4*w->sol[3]-w->sol[4];
   }
   else if (bcnr(1) != 0)
   {
      printf("ERROR: assign_gp, unknown bcnr(1)= %i\n", bcnr(1));
   }

/*!** Right side */
   if( bcnr(2) == 2 )
   {
/*! extrapolate u */
      w->sol[n+1] = 4*w->sol[n]-6*w->sol[n-1]+4*w->sol[n-2]-w->sol[n-3];
   }
   else if (bcnr(2) != 0)
   {
      printf("ERROR: assign_gp, unknown bcnr(2)= %i\n", bcnr(2));
   }
}

