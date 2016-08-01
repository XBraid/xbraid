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

/** \file braid_dist.c
 * \brief Source code for data distribution routines.  See braid_dist.h for more information.
 *
 */
#include "braid_dist.h"

/*----------------------------------------------------------------------------
 * Constructs the Balance struct needed by FRefine to complete load balancing and
 * temporal refinment.
 *----------------------------------------------------------------------------*/

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

    //Check if we need to refine in time.
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
        _braid_GetBlockDistInterval( core, bstruct );
        *done = 0;
    }

    return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Returns the index interval for 'proc' in a blocked data distribution
 *----------------------------------------------------------------------------*/
braid_Int
_braid_GetBlockDistInterval(braid_Core core,
                             _braid_BalanceStruct *bstruct )
{
    MPI_Comm comm = _braid_CoreElt( core, comm );

    braid_Int  npoints, nprocs, proc, new_ilower, new_iupper;
    npoints = _braid_BalanceElt( bstruct, refined_gupper ) + 1;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &proc);

    _braid_GetBlockDistInterval_basic( npoints, nprocs, proc, &new_ilower, &new_iupper );
    _braid_BalanceElt( bstruct, fine_ilower ) = new_ilower;
    _braid_BalanceElt( bstruct, fine_iupper )  = new_iupper;
    _braid_BalanceElt( bstruct, fine_gupper ) = npoints - 1;

    return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Returns the index interval for 'proc' in a blocked data distribution
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetBlockDistInterval_basic(braid_Int   npoints,
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
 * Initialize the data structure used to store the global wfactor information
 *----------------------------------------------------------------------------*/

braid_Int
_braid_WeightedStructInit( _braid_WeightedStruct *wstruct )
{
    _braid_WeightedElt( wstruct, local_sum ) = 0;
    _braid_WeightedElt( wstruct, local_min ) = braid_Int_Max;
    _braid_WeightedElt( wstruct, local_max ) = -1;
    _braid_WeightedElt( wstruct, global_max ) = -1;
    _braid_WeightedElt( wstruct, global_min) = braid_Int_Max;
    _braid_WeightedElt( wstruct, global_sum ) = -1;
    _braid_WeightedElt( wstruct, local_start ) = -1;
    _braid_WeightedElt( wstruct, local_stop ) = -1;

    return _braid_error_flag;
}


/*----------------------------------------------------------------------------
 * Destroy the weighted data structure
 *----------------------------------------------------------------------------*/

braid_Int
_braid_WeightedStructDestroy( _braid_WeightedStruct *wstruct )
{
    _braid_TFree( wstruct );
    return _braid_error_flag;
}


/*----------------------------------------------------------------------------
 * Completes the Balance structure needed by FRefine to complete load balancing
 * and temporal refinement, using the weighted data distribution.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_WeightedDist(braid_Core            core,
                    braid_Int             *done,
                    braid_Real             *wfactors,
                    braid_Int              npoints,
                    _braid_BalanceStruct  *bstruct )
{
    MPI_Comm comm = _braid_CoreElt( core, comm );

    _braid_WeightedStruct wstruct;
    _braid_WeightedStructInit( &wstruct );

    /* Get the local sum, min, and max; */
    braid_Int i;
    braid_Real local_max = 0;
    braid_Real local_sum = 0;
    braid_Real local_min = 1e100;
    braid_Int refined_gupper = 0;
    for ( i = 0; i < npoints; i++ )
    {
        if ( local_max < wfactors[i] )
            local_max = wfactors[i];
        if ( local_min > wfactors[i] && wfactors[i] > 0  )
            local_min = wfactors[i];
        local_sum += wfactors[i];
    }

    /* Do the allreduce */
    braid_Real *mpi_send = _braid_CTAlloc( braid_Real, 4 );
    braid_Real *mpi_recv = _braid_CTAlloc( braid_Real, 4 );
    mpi_send[0] = (braid_Real) npoints;
    mpi_send[1] = local_sum;
    mpi_send[2] = local_max;
    mpi_send[3] = local_min;

    MPI_Op sum_sum_max_min;
    MPI_Op_create( (MPI_User_function*) SumSumMaxMin , 1 , &sum_sum_max_min );
    MPI_Allreduce(mpi_send, mpi_recv, 4, braid_MPI_REAL, sum_sum_max_min, comm);
    MPI_Op_free( &sum_sum_max_min );
    refined_gupper = (braid_Int) mpi_recv[0];
    (refined_gupper)--;

    /* Save the information */
    _braid_WeightedElt( &wstruct, local_sum ) = local_sum;
    _braid_WeightedElt( &wstruct, local_max ) = local_max;
    _braid_WeightedElt( &wstruct, local_min ) = local_min;
    _braid_WeightedElt( &wstruct, global_sum ) = mpi_recv[1];
    _braid_WeightedElt( &wstruct, global_max ) = mpi_recv[2];
    _braid_WeightedElt( &wstruct, global_min ) = mpi_recv[3];
    _braid_BalanceElt( bstruct, refined_gupper ) = refined_gupper;

    /* Do the scan */
    braid_Real lbal = mpi_recv[3] - mpi_recv[2];
    braid_Int refine = _braid_BalanceElt( bstruct, refine );
    braid_Int coarse_gupper = _braid_BalanceElt( bstruct, coarse_gupper );
    braid_Int ref = refined_gupper - coarse_gupper ;
    braid_Int refined_iupper;

    /* Refine in time and do a weigted load balance */
    if ( refine && lbal != 0 && ref != 0 )
    {
        MPI_Scan(mpi_send, mpi_recv, 2, braid_MPI_REAL, MPI_SUM, comm);
        refined_iupper = (braid_Int) mpi_recv[0];

        _braid_BalanceElt( bstruct, refined_gupper ) = refined_gupper;
        _braid_BalanceElt( bstruct, refined_ilower ) = refined_iupper - npoints;
        _braid_BalanceElt( bstruct, refined_iupper ) = refined_iupper - 1;
        _braid_WeightedElt( &wstruct, local_start )    = mpi_recv[1] - local_sum;
        _braid_WeightedElt( &wstruct, local_stop ) = mpi_recv[1];
        *done = 0;
    }

    /* 2 Refine in time and dont use a weighted load balence because all weights are the same */
    else if ( refine && ref != 0 && lbal == 0 )
    {
        MPI_Scan(&npoints, &refined_iupper, 1, braid_MPI_INT, MPI_SUM, comm);
        _braid_BalanceElt( bstruct, refined_gupper ) = refined_gupper;
        _braid_BalanceElt( bstruct, refined_ilower ) = refined_iupper - npoints;
        _braid_BalanceElt( bstruct, refined_iupper ) = refined_iupper - 1;
        _braid_BalanceElt( bstruct, lbalance ) = 0;
        *done = 0;
    }

    //No refine in time
    else if ( !refine || ref == 0 )
    {
        /*No load balance and no refine so just return */
        if ( lbal == 0  )
        {
            _braid_BalanceElt( bstruct, refine ) = 0;
            _braid_BalanceElt( bstruct, lbalance ) = 0;
            *done = 1;
        }
        else
        {
            braid_Real c_sum_begin, c_sum_end;
            MPI_Scan(&local_sum, &c_sum_end, 1, braid_MPI_REAL, MPI_SUM, comm );
            c_sum_begin = c_sum_end - local_sum;
            _braid_BalanceElt( bstruct, refined_gupper ) = _braid_BalanceElt( bstruct, coarse_gupper );
            _braid_BalanceElt( bstruct, refined_ilower ) = _braid_BalanceElt( bstruct, coarse_ilower );
            _braid_BalanceElt( bstruct, refined_iupper ) = _braid_BalanceElt( bstruct, coarse_iupper );
            _braid_WeightedElt( &wstruct, local_start )    = c_sum_begin;
            _braid_WeightedElt( &wstruct, local_stop )     = c_sum_end;
            _braid_BalanceElt( bstruct, refine ) = 0;
            *done = 0;
        }

    }

    _braid_TFree( mpi_recv );
    _braid_TFree( mpi_send );

    if ( !(*done) )
    {
        /* Get the interval using either the weighted or block distribution */
        if ( _braid_BalanceElt( bstruct , lbalance ) )
            _braid_GetWeightedInterval( core, wfactors, &wstruct, bstruct );
        else
            _braid_GetBlockDistInterval( core, bstruct );
    }

    return _braid_error_flag;
}

/**************************************************************************************
 * Get the new fine grid interval using a weighted data distribution.
 *************************************************************************************/

braid_Int
_braid_GetWeightedInterval(braid_Core core,
                           braid_Real  *wfactors,
                           _braid_WeightedStruct *wstruct,
                           _braid_BalanceStruct *bstruct )
{

    braid_Int rank, comm_size;
    braid_Int i ,j;
    braid_Real goal_load;
    braid_Int *send_buffer_ilower, *send_buffer_iupper, num_recvs, num_sends;
    braid_Int partition_low, partition_high, num_partitions;

    MPI_Comm comm = _braid_CoreElt( core, comm );
    MPI_Request *send_requests, *recv_requests;
    MPI_Comm_rank( comm, &rank );
    MPI_Comm_size( comm, &comm_size );

    /* Extract information needed from the structures */
    braid_Int fine_ilower, fine_iupper;
    braid_Int refined_gupper = _braid_BalanceElt( bstruct, refined_gupper );
    braid_Int refined_ilower = _braid_BalanceElt( bstruct, refined_ilower );

    braid_Real global_sum =   _braid_WeightedElt( wstruct, global_sum );
    braid_Real global_max =   _braid_WeightedElt( wstruct, global_max );
    braid_Real global_min =   _braid_WeightedElt( wstruct, global_min );
    braid_Real local_sum =    _braid_WeightedElt( wstruct, local_sum );
    braid_Real c_sum_begin =  _braid_WeightedElt( wstruct, local_start );
    braid_Real c_sum_end =    _braid_WeightedElt( wstruct, local_stop );

    /* calculate the goal load on each processor and the number of processors to use */
    goal_load = global_sum / _braid_min( comm_size, refined_gupper + 1 );
    if ( goal_load >= global_max)
        num_partitions = (braid_Int) (floor( global_sum/goal_load + global_min/10000. ) - 1);
    else
    {
        goal_load  = global_max;
        num_partitions = (braid_Int) ( floor( global_sum/global_max + global_min/10000. ) - 1) ;
    }

    //Post recvs for ilower and iupper
    num_recvs = 0;
    recv_requests = _braid_CTAlloc( MPI_Request, 2 );
    if ( rank <= num_partitions )
    {
        if ( rank < num_partitions )
            MPI_Irecv( &fine_iupper, 1, MPI_INT, MPI_ANY_SOURCE, 18, comm, &recv_requests[num_recvs++] );
        else
            fine_iupper = refined_gupper;
        if ( rank > 0 )
            MPI_Irecv( &fine_ilower, 1, MPI_INT, MPI_ANY_SOURCE, 17, comm, &recv_requests[num_recvs++] );
        else
            fine_ilower = 0;
    }
    else
    {
        fine_ilower = refined_gupper + 1;
        fine_iupper = refined_gupper;
    }

    /* Check if any new partition boundaries lie between time points on the current grid. If so, send
       the ilower and iupper information to the corresponding processors.  */
    num_sends = 0;
    j = 0;

    if ( local_sum > 0 )
    {
        partition_low = (braid_Int) ceil(c_sum_begin/goal_load );
        /* Account for round off errors */
        if ( partition_low == 0 || fabs( ceil(c_sum_begin/goal_load) - floor(c_sum_begin/goal_load)) < 1e-13 )
            partition_low++;

        partition_high = (braid_Int) floor(c_sum_end/goal_load );
        partition_high = _braid_min( partition_high, num_partitions );

        send_buffer_ilower  = _braid_CTAlloc( braid_Int , (partition_high - partition_low + 1) );
        send_buffer_iupper  = _braid_CTAlloc( braid_Int , (partition_high - partition_low + 1) );
        send_requests = _braid_CTAlloc( MPI_Request, 2*(partition_high-partition_low + 1) ) ;

        for ( i = partition_low; i <= partition_high; i++ )
        {
            while ( c_sum_begin < goal_load*i )
                c_sum_begin += wfactors[j++];
            send_buffer_ilower[(i-partition_low)] = j + refined_ilower;
            send_buffer_iupper[(i-partition_low)] = j-1 + refined_ilower;

            MPI_Isend( &send_buffer_ilower[i-partition_low], 1, braid_MPI_INT, i, 17, comm, &send_requests[num_sends++] );
            MPI_Isend( &send_buffer_iupper[i-partition_low], 1, braid_MPI_INT, i - 1, 18, comm, &send_requests[num_sends++] );
        }
    }

    /* Complete communication */
    if ( num_recvs > 0 )
        MPI_Waitall( num_recvs , recv_requests , MPI_STATUS_IGNORE);
    if ( num_sends > 0 )
        MPI_Waitall( num_sends , send_requests , MPI_STATUS_IGNORE);

    if ( local_sum > 0 )
    {
        _braid_TFree( send_requests );
        _braid_TFree( send_buffer_ilower );
        _braid_TFree( send_buffer_iupper );
    }
    _braid_TFree( recv_requests );

    _braid_BalanceElt( bstruct, fine_ilower ) = fine_ilower;
    _braid_BalanceElt( bstruct, fine_iupper ) = fine_iupper;
    _braid_BalanceElt( bstruct, fine_gupper ) = refined_gupper;

    MPI_Barrier ( comm );
    return _braid_error_flag;
}

/*********************************************************************************************
 * Used in the weighted data distribution all reduce call, this function combines all neccesary
 * reductions into a single call. Note that length of the in and out buffers must be 4 */

void SumSumMaxMin(braid_Real *in, braid_Real *inout, braid_Int *len, MPI_Datatype *datatype)
{
    inout[0] = in[0] + inout[0];
    inout[1] = in[1] + inout[1];
    inout[2] = _braid_max( in[2] , inout[2] );
    inout[3] = _braid_min( in[3] , inout[3] );
}
