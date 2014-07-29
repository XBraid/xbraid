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


