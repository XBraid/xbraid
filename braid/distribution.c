/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. For support, post issues to the XBraid Github page.
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

#include "_braid.h"
#include "_util.h"

/*----------------------------------------------------------------------------
 * Returns the index interval for 'proc' in a blocked data distribution
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetBlockDistInterval(braid_Int   npoints,
                            braid_Int   nprocs,
                            braid_Int   proc,
                            braid_Int  *ilower_ptr,
                            braid_Int  *iupper_ptr)
{
   braid_Int  ilower, iupper, quo, rem, p;

   quo = npoints/nprocs;
   rem = npoints%nprocs;

   p = proc;
   ilower = p*quo + (p < rem ? p : rem);
   p = proc+1;
   iupper = p*quo + (p < rem ? p : rem) - 1;

   *ilower_ptr = ilower;
   *iupper_ptr = iupper;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Returns the processor that owns 'index' in a blocked data distribution
 * (returns -1 if index is out of range)
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetBlockDistProc(braid_Int   npoints,
                        braid_Int   nprocs,
                        braid_Int   index,
                        braid_Int  *proc_ptr)
{
   braid_Int      proc, quo, rem, p, q;

   /* Compute processor number */
   if ((index < 0) || (index > (npoints-1)))
   {
      proc = -1;
   }
   else
   {
      quo = npoints/nprocs;
      rem = npoints%nprocs;

      if (quo > 0)
      {
         p = index/(quo+1);
         q = (index - rem*(quo+1))/quo;
         proc = (p < rem ? p : rem+q);
      }
      else
      {
         proc = index;
      }
   }

   *proc_ptr = proc;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Returns the index interval for my processor on the finest grid level
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetDistribution(braid_Core   core,
                       braid_Int   *ilower_ptr,
                       braid_Int   *iupper_ptr)
{
   MPI_Comm   comm    = _braid_CoreElt(core, comm);
   braid_Int  gupper = _braid_CoreElt(core, gupper);
   braid_Int  npoints, nprocs, proc;

   npoints = gupper + 1;
   MPI_Comm_size(comm, &nprocs);
   MPI_Comm_rank(comm, &proc);

   _braid_GetBlockDistInterval(npoints, nprocs, proc, ilower_ptr, iupper_ptr);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Returns the processor that owns 'index' on the given grid 'level'
 * (returns -1 if index is out of range)
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetProc(braid_Core   core,
               braid_Int    level,
               braid_Int    index,
               braid_Int   *proc_ptr)
{
   MPI_Comm       comm   = _braid_CoreElt(core, comm);
   _braid_Grid  **grids  = _braid_CoreElt(core, grids);
   braid_Int      gupper = _braid_CoreElt(core, gupper);
   braid_Int      npoints, nprocs;
   braid_Int      l, cfactor;

   npoints = gupper + 1;
   MPI_Comm_size(comm, &nprocs);
   /* Map index to the finest grid */
   for (l = level-1; l > -1; l--)
   {
      cfactor = _braid_GridElt(grids[l], cfactor);
      _braid_MapCoarseToFine(index, cfactor, index);
   }

   _braid_GetBlockDistProc(npoints, nprocs, index, proc_ptr);

   return _braid_error_flag;
}

