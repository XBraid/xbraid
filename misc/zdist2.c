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

#include <stdlib.h>
#include <stdio.h>

int
main( int argc, char *argv[] )
{
   int  arg_index;
   int  npoints, nprocs;
   int  proc, ilower, iupper, i, p, q, last_proc;

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

   /* GetDistribution */

   printf("\n GetDistribution:\n\n");
   for (proc = 0; proc < nprocs; proc++)
   {
      p = proc;
      ilower = npoints - ((nprocs-p)*npoints)/nprocs;
      p = proc+1;
      iupper = npoints - ((nprocs-p)*npoints)/nprocs - 1;

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
      proc = (i*nprocs)/npoints;
      
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
