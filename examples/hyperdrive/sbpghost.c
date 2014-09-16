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


#include "c_array.h"
   
/*    subroutine OPWGHOST( nb, wb, bop, gh, betapcoeff )*/
void sbpghost( int nb, int wb, double_array_2d * bop_, double *gh, double betapcoeff)
{
   double beta, p1;
#define bop(i,j) compute_index_2d(bop_, i, j)
   p1   = 13649.0/43200.0; /* this is the weight in the 6th order scalar product for point 1 (on the bndry)*/
   beta = betapcoeff/p1;
   
   *gh      =            beta*( 1);
   bop(1,1) = bop(1,1) + beta*(-4);
   bop(1,2) = bop(1,2) + beta*( 6);
   bop(1,3) = bop(1,3) + beta*(-4);
   bop(1,4) = bop(1,4) + beta*( 1);
}


