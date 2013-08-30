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

#include "_warp.h"

/*--------------------------------------------------------------------------
 * Project an interval onto a strided index space that contains the index
 * 'index' and has stride 'stride'.  An empty projection is represented by
 * ilower > iupper.
 *--------------------------------------------------------------------------*/

warp_Int
_warp_ProjectInterval( warp_Int   ilower,
                       warp_Int   iupper,
                       warp_Int   index,
                       warp_Int   stride,
                       warp_Int  *pilower,
                       warp_Int  *piupper )
{
   warp_Int  il, iu;

   il = ilower - index;
   iu = iupper - index;

   if ( il <= 0 )
      il = (warp_Int) (il / stride);
   else
      il = (warp_Int) ((il + (stride-1)) / stride);

   if ( iu >= 0 )
      iu = (warp_Int) (iu / stride);
   else
      iu = (warp_Int) ((iu - (stride-1)) / stride);

   *pilower = index + il * stride;
   *piupper = index + iu * stride;

   return _warp_error_flag;
}


/*--------------------------------------------------------------------------
 * Determine the accuracy used for the spatial solves based on the ratio of
 * the current residual norm and the stopping tolerance. 
 *--------------------------------------------------------------------------*/

warp_Int
_warp_SetAccuracy( warp_Float  rnorm,
                   warp_Float  loose_tol,
                   warp_Float  tight_tol,
                   warp_Float  oldAccuracy,
                   warp_Float  tol,
                   warp_Float *paccuracy )
{
   warp_Float accuracy, ttol, stol;

   warp_Float loose_stol = -log10(loose_tol);
   warp_Float tight_stol = -log10(tight_tol);
   warp_Float tight_ttol = (8.0 / 9)*(-log10(tol));

   /*ttol = -log10(rnorm / rnorm0);*/
   ttol = -log10(rnorm);
   if ( ttol >= tight_ttol )
      stol = tight_stol;
   else
   {
      stol = (ttol / tight_ttol) * (tight_stol - loose_stol) + loose_stol;
   }

   accuracy = pow(10, -stol);
   if (accuracy > oldAccuracy)
      accuracy = oldAccuracy;

   *paccuracy = accuracy;

   return _warp_error_flag;
}


