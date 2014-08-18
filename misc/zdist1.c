/*BHEADER**********************************************************************
 * Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of XBraid.  See file COPYRIGHT for details.
 *
 * XBraid is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 ***********************************************************************EHEADER*/

#include <stdlib.h>
#include <stdio.h>

int
main( int argc, char *argv[] )
{
   int  arg_index;
   int  npoints, nprocs;
   int  proc, quo,rem, ilower, iupper, i, p, q, last_proc;

   if (argc > 2)
   {
      arg_index = 1;
      npoints = atoi(argv[arg_index++]);
      nprocs  = atoi(argv[arg_index++]);
   }
   else
   {
      printf("Problem with input\n");
      exit(1);
   }

   printf("\n npoints = %d, nprocs = %d\n", npoints, nprocs);

   quo = npoints/nprocs;
   rem = npoints%nprocs;

   /* GetDistribution */

   printf("\n GetDistribution:\n\n");
   for (proc = 0; proc < nprocs; proc++)
   {
      p = proc;
      ilower = p*quo + (p < rem ? p : rem);
      p = proc+1;
      iupper = p*quo + (p < rem ? p : rem) - 1;

      printf(" %2d:", proc);
      for (i = ilower; i <= iupper; i++)
      {
         printf(" %2d", i);
      }
      printf("\n");
   }

   /* GetProc */

   printf("\n GetProc:\n");
   last_proc = -1;
   for (i = 0; i < npoints; i++)
   {
      p = i/(quo+1);
      q = (i - rem*(quo+1))/quo;
      proc = (p < rem ? p : rem+q);
      
      if (proc > last_proc)
      {
         printf("\n");
         printf(" %2d:", proc);
      }
      printf(" %2d", i);

      last_proc = proc;
   }
   printf("\n\n");

   return (0);
}
