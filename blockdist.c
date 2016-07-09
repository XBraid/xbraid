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

/** \file bockdist.c
 * \brief Source code for block distrobution routines.  See blockdist.h for more information.
 *
 */
#include "blockdist.h"

braid_Int
_braid_BlockDist(braid_Core            core,
                 braid_Int             *done,
                 braid_Int              npoints,
                 _braid_BalanceStruct  *bstruct )
{

    MPI_Comm comm = _braid_CoreElt( core, comm );
    braid_Int refine = _braid_BalanceElt( bstruct, refine );
    braid_Int coarse_gupper = _braid_BalanceElt( bstruct, coarse_gupper );
    braid_Int refined_gupper, refined_iupper;

    //Check if we even need to refine anywhere
    refined_gupper = coarse_gupper;
    if ( refine )
    {
        MPI_Allreduce(&npoints, &refined_gupper, 1, braid_MPI_INT, MPI_SUM, comm);
        (refined_gupper)--;
    }

    //If no refinment occured -- set refine = 0, and return
    if ( refined_gupper == coarse_gupper )
    {
        _braid_BalanceElt( bstruct, refine ) = 0;
        *done = 1;
    }
    else
    {
        //Do a scan to get refined_ilower and refined_iupper
        _braid_CoreElt(core, r_space) = 0;
        MPI_Scan(&npoints, &refined_iupper, 1, braid_MPI_INT, MPI_SUM, comm);
        _braid_BalanceElt( bstruct, refined_ilower ) = refined_iupper - npoints;
        _braid_BalanceElt( bstruct, refined_iupper ) = refined_iupper - 1;
        _braid_BalanceElt( bstruct, refined_gupper ) = refined_gupper;
        _braid_BalanceElt( bstruct, fine_gupper ) = refined_gupper;
        _braid_GetBlockDistInterval1( core, bstruct );
        *done = 0;
    }

    return _braid_error_flag;
}

braid_Int
_braid_GetBlockDistInterval1(braid_Core core,
                             _braid_BalanceStruct *bstruct )
{
    MPI_Comm comm = _braid_CoreElt( core, comm );

    braid_Int  npoints, nprocs, proc, new_ilower, new_iupper;
    npoints = _braid_BalanceElt( bstruct, refined_gupper ) + 1;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &proc);

    _braid_GetBlockDistInterval( npoints, nprocs, proc, &new_ilower, &new_iupper );
    _braid_BalanceElt( bstruct, fine_ilower ) = new_ilower;
    _braid_BalanceElt( bstruct, fine_iupper )  = new_iupper;
    _braid_BalanceElt( bstruct, fine_gupper ) = npoints - 1;

    return _braid_error_flag;
}

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

