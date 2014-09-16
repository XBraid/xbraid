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

/* evaluate the L2 and Linf norm of the difference between grid functions w and we */
void
evaldiff( grid_fcn *w, grid_fcn *we, double *l2, double *li )
{
   int i;
   double locerr;
   int n = w->n;
   double h = w->h;
   
   *l2 = 0;
   *li = 0;
   
   for (i=1; i<n; i++)
   {
      locerr = fabs(w->sol[i]-we->sol[i]);
      *l2 += h*locerr*locerr;
      if( locerr > *li )
      {
         *li = locerr;
      }
   }
   
   *l2 = sqrt(*l2);
}

/* evaluate the L2 and Linf norm of grid function w */
void
evalnorm( grid_fcn *w, double *l2, double *li )
{
   int i;
   double locerr;
   int n = w->n;
   double h = w->h;
   
   *l2 = 0;
   *li = 0;
   
/* for periodic grid functions, only the first n-1 values are unique */
   for (i=1; i<n; i++)
   {
      locerr = fabs(w->sol[i]);
      *l2 += h*locerr*locerr;
      if( locerr > *li )
      {
         *li = locerr;
      }
   }
/* for periodic grid functions, only the first n-1 values are unique */
   *l2 = sqrt(*l2);
}

