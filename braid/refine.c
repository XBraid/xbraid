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

#define DEBUG 0

#if DEBUG
braid_Int  FRefine_count = 0;
#endif

/*----------------------------------------------------------------------------
 * Macros used below
 *----------------------------------------------------------------------------*/

/* Compute number of reals given some number of bytes (use ceiling) */
#define _braid_NBytesToNReals(nbytes, nreals) \
nreals = nbytes / sizeof(braid_Real) + ((nbytes % sizeof(braid_Real)) != 0)

/*----------------------------------------------------------------------------
 * Create a new fine grid (level 0) and corresponding grid hierarchy by refining
 * the current fine grid based on user-provided refinement factors.  Return the
 * boolean 'refined_ptr' to indicate whether grid refinement was actually done.
 * To simplify the algorithm, refinement factors are automatically capped to be
 * no greater than the coarsening factor (for level 0).  The grid data is also
 * redistributed to achieve good load balance in the temporal dimension.  If the
 * refinement factor is 1 in each time interval, no refinement is done.
 *
 * This routine is somewhat complex, but an attempt was made to use consistent
 * terminology throughout.  We refer to the initial level 0 grid as the "coarse"
 * grid and the new level 0 grid as the "fine" grid.  The routine starts with
 * the coarse grid and creates an intermediate "refined" grid that is co-located
 * with the parent coarse grid data distribution.  On this refined grid, a
 * number of things are constructed: the mapping between the coarse and fine
 * indexes; the new fine time values; and the injected and (possibly) spatially
 * refined coarse u-vectors.  This data is then redistributed by first sending
 * the mapping and time value information to the appropriate processors (the
 * receiving side polls for arbitrary messages until each fine index is
 * accounted for).  The u-vector values are communicated in a second phase
 * (without polling).  Finally, a new hierarchy is created and the fine grid
 * values are initialized by integrating the communicated u-vector values to the
 * next C-point to the right.  Note that in the case of C-point storage, some
 * u-vector values may not need to be communicated.
 *
 * The variable names use certain conventions as well.  No prefix (usually)
 * indicates a coarse variable, and the prefixes 'r_' and 'f_' indicate data on
 * the "refined" and "fine" grids, respectively.  Characters 'c', 'f', and 'a'
 * usually mean "coarse", "fine", and "array".  Indexes 'i', 'r_i', and 'f_i'
 * refer to global indexes, while 'ii', 'r_ii', and 'f_ii' are local indexes.
 *
 * Here are some important variables along with a brief description (they are
 * illustrated further in the example below):
 *
 *   gupper, f_gupper - global upper index for coarse and fine grids
 *   ilower,   iupper,   npoints   - extents and size of local coarse interval
 *   r_ilower, r_iupper, r_npoints - extents and size of local refined interval
 *   f_ilower, f_iupper, f_npoints - extents and size of local fine interval
 *
 *   r_ca - index map from fine to coarse on the refined grid  (size 'r_npoints')
 *   r_ta - time values on the refined grid                    (size 'r_npoints+2')
 *   r_fa - index map from coarse to fine on the refined grid  (size 'npoints+1')
 *          (note the extra value)
 *   f_ca - index map from fine to coarse on the fine grid     (size 'f_npoints')
 *
 *   send_ua - array of u-vectors to send to new processors    (size 'npoints')
 *   recv_ua - array of u-vectors received from old processors (size 'f_npoints')
 *
 *   f_first - fine grid index for the first coarse point in my fine interval
 *   f_next  - fine grid index for the first coarse point in the next fine interval
 *
 * Example: Some processor p owns the coarse interval, ilower = 29, iupper = 33.
 * The coarsening factor is 3 and 'rfactors' indicates the refinement factor for
 * the coarse interval to the left.  From this, an intermediate refined grid and
 * a final fine grid are formed as described above.
 *
 *   Coarse            |-----------c-----------|-----------|-----------c
 *     Grid           29          30          31          32          33
 * 
 * rfactors            3           1           2           3           2
 * 
 *  Refined ---c---|---|-----------c-----|-----|---c---|---|-----c-----|
 *     Grid   57  58  59          60    61    62  63  64  65    66    67
 * 
 * r_ilower   57
 *     r_ca   -1  -1  29          30    -1    31  -1  -1  32    -1    33
 *     r_ta    *   *   *           *     *     *   *   *   *     *     *
 *     r_fa           59          60          62          65          67          70
 *  send_ua            *           *           *           *           *
 * 
 *     Fine      |-----|---c---|---|-----c-----|---|---c---|-----|-----c
 *     Grid     61    62  63  64  65    66    67  68  69  70    71    72
 *              
 * f_ilower     61
 *     f_ca     -1    31  -1  -1  32    -1    33  -1  -1  34    -1    35
 *  f_first           62
 *   f_next                                                                       74
 *  recv_ua      0     *   0   0   *     0     *   0   0   *     0     *
 * 
 * When storing C-pts only, data from coarse indices 29 and 34 are not needed,
 * so we have the following differences from above:
 *     r_ca   -1  -1  -1          30    -1    31  -1  -1  32    -1    33
 *  send_ua            0           *           *           *           *
 *     f_ca     -1    31  -1  -1  32    -1    33  -1  -1  -1    -1    35
 *  recv_ua      0     *   0   0   *     0     *   0   0   0     0     *
 * 
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FRefine(braid_Core   core,
               braid_Int   *refined_ptr)
{
   MPI_Comm           comm            = _braid_CoreElt(core, comm);
   braid_App          app             = _braid_CoreElt(core, app);
   braid_Int          iter            = _braid_CoreElt(core, niter);
   braid_Int          ichunk          = _braid_CoreElt(core, ichunk);
   braid_Int          refine          = _braid_CoreElt(core, refine);
   braid_Int         *rfactors        = _braid_CoreElt(core, rfactors);
   braid_Int          nrefine         = _braid_CoreElt(core, nrefine);
   braid_Int          max_refinements = _braid_CoreElt(core, max_refinements);
   braid_Int          tpoints_cutoff  = _braid_CoreElt(core, tpoints_cutoff);
   braid_AccessStatus astatus         = (braid_AccessStatus)core;
   braid_BufferStatus bstatus         = (braid_BufferStatus)core;
   braid_Int          access_level    = _braid_CoreElt(core, access_level);
   _braid_Grid      **grids           = _braid_CoreElt(core, grids);
   braid_Int          ncpoints        = _braid_GridElt(grids[0], ncpoints);

   braid_Real         rnorm;

   /* Use prefix 'r_' for the refined version of the current grid level 0, and
    * use prefix 'f_' for the fine grid in the new distribution. */

   braid_Int         npoints, ilower, iupper, gupper, i, j, ii;
   braid_Int         r_npoints, r_ilower, r_iupper, r_i, r_ii;
   braid_Int         f_npoints, f_ilower, f_iupper, f_gupper, f_i, f_j, f_ii;
   braid_Int        *r_ca, *r_fa, *f_ca, f_first, f_next, next;
   braid_Real       *ta, *r_ta_alloc, *r_ta, *f_ta;

   braid_BaseVector *send_ua, *recv_ua, u;
   braid_Int        *send_procs, *recv_procs, *send_unums, *recv_unums, *iptr;
   braid_Real       *send_buffer, *recv_buffer, **send_buffers, **recv_buffers, *bptr;
   void             *buffer;
   braid_Int         send_size, recv_size, *send_sizes, size, isize, max_usize;
   braid_Int         ncomms, nsends, nrecvs, nreceived, nprocs, proc, prevproc;
   braid_Int         unum, send_msg, recv_msg;
   MPI_Request      *requests, request;
   MPI_Status       *statuses, status;
                   
   _braid_Grid      *f_grid;
   braid_Int         cfactor, rfactor, m, interval, flo, fhi, fi, ci, f_hi, f_ci;

#if DEBUG
   braid_Int  myproc;
   /*cfactor = 6;*/ /* RDF HACKED TEST */
   MPI_Comm_rank(comm, &myproc);
#endif

   /* Only refine if refinement is turned on */
   if(refine == 0)
   {
      *refined_ptr = 0;
      return _braid_error_flag;
   }

   gupper  = _braid_CoreElt(core, gupper);
   ilower  = _braid_GridElt(grids[0], ilower);
   iupper  = _braid_GridElt(grids[0], iupper);
   npoints = iupper - ilower + 1;
   
   /* If reached max refinements or have too many time points, stop refining */
   if( !((nrefine < max_refinements) && (gupper < tpoints_cutoff)) )
   {
      _braid_CoreElt(core, refine)   = 0;
      _braid_CoreElt(core, rstopped) = iter;
      _braid_FRefineSpace(core, refined_ptr);
      return _braid_error_flag;
   }

   /*-----------------------------------------------------------------------*/
   /* 1. Compute f_gupper and the local interval extents for both the refined
    * and fine grids.  The local refined interval contains the fine grid points
    * underlying the coarse interval (ilower-1, iupper]. */

   /* Compute f_gupper and r_npoints */
   _braid_GetCFactor(core, 0, &cfactor);
   r_npoints = 0;
   for (i = ilower; i <= iupper; i++)
   {
      /* First modify rfactor to be no greater than cfactor (required) */
      ii = i - ilower;
      if (rfactors[ii] < 1)
      {
         _braid_Error(braid_ERROR_GENERIC, "Refinement factor smaller than one");
         rfactors[ii] = 1;
      }
      r_npoints += rfactors[i-ilower];
   }
   MPI_Allreduce(&r_npoints, &f_gupper, 1, braid_MPI_INT, MPI_SUM, comm);
   f_gupper--;

#if DEBUG
   for (i = ilower; i <= iupper; i++)
   {
      ii = i - ilower;
      printf("%d %d: 0 rfactor %d = %d\n", FRefine_count, myproc, i, rfactors[ii]);
   }
#endif

   /* Check to see if we need to refine, and return if not */
   if (f_gupper == gupper)
   {
      _braid_FRefineSpace(core, refined_ptr);
      return _braid_error_flag;
   }
   else
   {
      _braid_CoreElt(core, r_space) = 0;
   }
      
   /* Compute r_ilower and r_iupper */
   MPI_Scan(&r_npoints, &r_iupper, 1, braid_MPI_INT, MPI_SUM, comm);
   r_ilower = r_iupper - r_npoints;
   r_iupper = r_iupper - 1;

   /* Compute f_ilower, f_iupper, and f_npoints for the final distribution */
   MPI_Comm_size(comm, &nprocs);
   MPI_Comm_rank(comm, &proc);
   _braid_GetBlockDistInterval((f_gupper+1), nprocs, proc, &f_ilower, &f_iupper);
   f_npoints = f_iupper - f_ilower + 1;

   /* Initialize the new fine grid */
   _braid_GridInit(core, 0, f_ilower, f_iupper, &f_grid);

   /*-----------------------------------------------------------------------*/
   /* 2. On the refined grid, compute the mapping between coarse and fine
    * indexes (r_ca, r_fa) and the fine time values (r_ta). */

   r_ca = _braid_CTAlloc(braid_Int,  r_npoints);
   r_ta_alloc = _braid_CTAlloc(braid_Real, r_npoints+2);
   r_ta = &r_ta_alloc[1];
   r_fa = _braid_CTAlloc(braid_Int,  npoints+1);
   ta = _braid_GridElt(grids[0], ta);

   r_ta[-1]=ta[-1];
   r_ii = 0;
   for (i = (ilower-1); i < iupper; i++)
   {
      ii = i-ilower;
      rfactor = rfactors[ii+1];

      for (j = 1; j <= rfactor; j++)
      {
         if (j < rfactor)
         {
            r_ca[r_ii] = -1;
            /* This works because we have ta[-1] */
            r_ta[r_ii] = ta[ii] + (((braid_Real)j)/rfactor)*(ta[ii+1]-ta[ii]);
         }
         else
         {
            r_ca[r_ii] = i+1;
            r_ta[r_ii] = ta[ii+1];
            r_fa[ii+1] = r_ilower + r_ii;
         }

         r_ii++;
      }
   }

   /* Get the next r_fa value to my right */
   ncomms = 2; /* Upper bound */
   requests = _braid_CTAlloc(MPI_Request, ncomms);
   statuses = _braid_CTAlloc(MPI_Status,  ncomms);
   ncomms = 0;
   if (npoints > 0)
   {
      braid_Real send_buf[2], recv_buf[2];
      send_buf[0]=r_fa[0];
      send_buf[1]=r_ta[0];     
      /* Post r_fa receive */
      r_fa[npoints] = f_gupper+1;
      if (iupper < gupper)
      {
         MPI_Irecv(recv_buf, 2, braid_MPI_REAL, MPI_ANY_SOURCE, 2, comm,
                   &requests[ncomms++]);
      }

      /* Post r_fa send */
      if (ilower > 0)
      {
         _braid_GetBlockDistProc((gupper+1), nprocs, (ilower-1), &prevproc);
         MPI_Isend(send_buf, 2, braid_MPI_REAL, prevproc, 2, comm,
                   &requests[ncomms++]);
      }
      MPI_Waitall(ncomms, requests, statuses);
      if ( iupper < gupper )
      {
         r_fa[npoints]=recv_buf[0];
         r_ta[r_npoints]=recv_buf[1];
      }
   }
   _braid_TFree(requests);
   _braid_TFree(statuses);

   /* If storing only C-points on the fine grid (and NOT using shell vectors),
    *  modify r_ca to mark only those coarse points that need to be sent to
    * initialize the C-points */
   if (_braid_CoreElt(core, storage) != 0 && _braid_CoreElt(core, useshell) != 1)
   {
      for (ii = 0; ii < npoints; ii++)
      {
         /* If the index for the next coarse point is not larger than the index
          * for the next C-point, then the coarse point is not needed */
         if ( !(r_fa[ii+1] > _braid_NextCPoint(r_fa[ii], cfactor)) )
         {
            r_ii = r_fa[ii] - r_ilower;
            r_ca[r_ii] = -1;
         }
      }
   }

   /*-----------------------------------------------------------------------*/
   /* 3. Send the index mapping and time value information (r_ca, r_ta) to the
    * appropriate processors to build index mapping and time value information
    * for the fine grid (f_ca, f_ta).  Also compute f_first and f_next. */

   /* Post f_next receive */
   f_next = -1;
   if (f_npoints > 0)
   {
      f_next = f_gupper+1;
      if (f_iupper < f_gupper)
      {
         MPI_Irecv(&f_next, 1, braid_MPI_INT, MPI_ANY_SOURCE, 3, comm, &request);
      }
   }

   /* Compute send information and send f_next info */
   size = 2*sizeof(braid_Int);         /* size of two integers */
   _braid_NBytesToNReals(size, isize); /* convert to units of braid_Real */
   send_procs  = _braid_CTAlloc(braid_Int,  r_npoints);
   send_sizes  = _braid_CTAlloc(braid_Int,  r_npoints);
   send_buffer = _braid_CTAlloc(braid_Real, r_npoints*(1+isize+1));
   bptr = send_buffer;
   nsends = -1;
   _braid_GetBlockDistProc((f_gupper+1), nprocs, (r_ilower-1), &prevproc);
   ii = 0;
   for (r_ii = 0; r_ii < r_npoints; r_ii++)
   {
      r_i = r_ilower + r_ii;
      _braid_GetBlockDistProc((f_gupper+1), nprocs, r_i, &proc);
      if ((proc != prevproc) || (nsends < 0))
      {
         nsends++;
         send_procs[nsends] = proc;
         bptr++; /* leave room for size value */

         if ((proc != prevproc) && (prevproc > -1))
         {
            /* Send f_next info */
            MPI_Send(&r_fa[ii], 1, braid_MPI_INT, prevproc, 3, comm);
         }
         prevproc = proc;
      }
      send_sizes[nsends] += (isize+1);

      iptr = (braid_Int *) bptr;
      iptr[0] = r_i;
      iptr[1] = r_ca[r_ii];
      bptr += isize;
      bptr[0] = r_ta[r_ii];
      bptr++;

      /* Update f_next info */
      if (r_fa[ii] == r_i)
      {
         ii++;
      }
   }
   nsends++;

#if DEBUG
   for (m = 0; m < nsends; m++)
   {
      size = send_sizes[m];
      proc = send_procs[m];
      printf("%d %d: 1 send %d, proc = %d, send size = %d\n",
             FRefine_count, myproc, m, proc, size);
   }
#endif

   requests = _braid_CTAlloc(MPI_Request, nsends);
   statuses = _braid_CTAlloc(MPI_Status,  nsends);

   /* Post sends (do this first, since we will poll on receives) */
   bptr = send_buffer;
   for (m = 0; m < nsends; m++)
   {
      size = send_sizes[m];
      proc = send_procs[m];
      bptr[0] = (braid_Real) size; /* insert size at the beginning */
      MPI_Isend(bptr, (1+size), braid_MPI_REAL, proc, 4, comm, &requests[m]);
      bptr += (1+size);
   }

   /* Post receives */
   recv_size = f_npoints*(1+isize+1); /* max receive size */
   recv_buffer = _braid_CTAlloc(braid_Real, recv_size);
   f_ca = _braid_CTAlloc(braid_Int,  f_npoints);
   f_ta = _braid_GridElt(f_grid, ta);
   nreceived = 0;
   while (nreceived < f_npoints)
   {
      /* post receive from arbitrary process (should always get at least one) */
      bptr = recv_buffer;
      size = recv_size;
      MPI_Recv(bptr, (1+size), braid_MPI_REAL, MPI_ANY_SOURCE, 4, comm, &status);

      size = (braid_Int) bptr[0];
      bptr++;
      for (j = 0; j < size; j += (isize+1))
      {
         iptr = (braid_Int *) bptr;
         f_i = iptr[0];
         f_ii = f_i - f_ilower;
         f_ca[f_ii] = iptr[1];
         bptr += isize;
         f_ta[f_ii] = bptr[0];
         bptr++;
         {
#if DEBUG
            printf("%d %d: 1 f_i = %02d, f_ca = %2d, f_ta = %f  (recv %2d) \n",
                   FRefine_count, myproc, f_i, f_ca[f_ii], f_ta[f_ii], nreceived);
#endif
         }
         nreceived++;
      }

   }

   /* Finish sends and f_next receive */
   MPI_Waitall(nsends, requests, statuses);
   if (f_npoints > 0)
   {
      if (f_iupper < f_gupper)
      {
         MPI_Wait(&request, &status);
      }
   }

   /* Compute f_first */
   f_first = f_next;
   for (f_ii = 0; f_ii < f_npoints; f_ii++)
   {
      if (f_ca[f_ii] > -1)
      {
         f_first = f_ilower + f_ii;
         break;
      }
   }

#if DEBUG
   printf("%d %d: 2 f_first = %d, f_next = %d\n",
          FRefine_count, myproc, f_first, f_next);
#endif

   /* Free up some memory */
   _braid_TFree(requests);
   _braid_TFree(statuses);
   _braid_TFree(send_procs);
   _braid_TFree(send_sizes);
   _braid_TFree(send_buffer);
   _braid_TFree(recv_buffer);

   /*-----------------------------------------------------------------------*/
   /* 4. Build u-vectors on the fine grid (send_ua) by first integrating on the
    * coarse grid, then injecting and refining spatially.  Redistribute these
    * u-vectors to the fine grid (recv_ua). */

   send_ua = _braid_CTAlloc(braid_BaseVector, npoints);
   send_procs = _braid_CTAlloc(braid_Int, npoints);
   send_unums = _braid_CTAlloc(braid_Int, npoints);
   send_buffers = _braid_CTAlloc(braid_Real *, npoints);

   recv_ua = _braid_CTAlloc(braid_BaseVector, f_npoints);
   recv_procs = _braid_CTAlloc(braid_Int, f_npoints);
   recv_unums = _braid_CTAlloc(braid_Int, f_npoints);
   recv_buffers = _braid_CTAlloc(braid_Real *, f_npoints);

   _braid_GetRNorm(core, -1, &rnorm);

   _braid_UCommInitF(core, 0);

   /* Start from the right-most interval */
   for (interval = ncpoints; interval > -1; interval--)
   {
      _braid_GetInterval(core, 0, interval, &flo, &fhi, &ci);

      /* Integrate F-points and refine in space */
      if (flo <= fhi)
      {
         _braid_UGetVector(core, 0, flo-1, &u);
         for (fi = flo; fi <= fhi; fi++)
         {
            _braid_Step(core, 0, fi, NULL, u);
            _braid_USetVector(core, 0, fi, u, 0); /* needed for communication */

            /* Set send_ua */
            ii = fi - ilower;
            r_ii = r_fa[ii] - r_ilower;
            if (r_ca[r_ii] > -1)
            {
               _braid_RefineBasic(core, -1, fi, &r_ta[r_ii], &ta[ii], u, &send_ua[ii]);
            }

            /* Allow user to process current vector */
            if( (access_level >= 3) )
            {
               _braid_AccessStatusInit(ta[ii], fi, ichunk, rnorm, iter, 0, nrefine, gupper,
                                       0, 0, braid_ASCaller_FRefine, astatus);
               _braid_AccessVector(core, astatus, u);
            }
         }
         _braid_BaseFree(core, app,  u);
      }

      /* Refine C-points in space */
      if (ci > -1)
      {
         _braid_UGetVectorRef(core, 0, ci, &u);

         /* Set send_ua */
         ii = ci - ilower;
         r_ii = r_fa[ii] - r_ilower;
         if (r_ca[r_ii] > -1)
         {
            _braid_RefineBasic(core, -1, ci, &r_ta[r_ii], &ta[ii], u, &send_ua[ii]);
         }

         /* Allow user to process current vector */
         if( (access_level >= 3) )
         {
            _braid_AccessStatusInit(ta[ii], ci, ichunk, rnorm, iter, 0, nrefine, gupper,
                                    0, 0, braid_ASCaller_FRefine, astatus);
            _braid_AccessVector(core, astatus, u);
         }
      }
   }

   _braid_UCommWait(core, 0);

   /* Compute nsends, send_procs, and send_unums from send_ua array */
   nsends = -1;
   prevproc = -1;
   for (ii = 0; ii < npoints; ii++)
   {
      if (send_ua[ii] != NULL)
      {
         r_i = r_fa[ii];
         _braid_GetBlockDistProc((f_gupper+1), nprocs, r_i, &proc);
         if (proc != prevproc)
         {
            nsends++;
            send_procs[nsends] = proc;
            send_unums[nsends] = 0;
            prevproc = proc;
         }
         send_unums[nsends]++;
      }
   }
   nsends++;

   /* Compute nrecvs, recv_procs, and recv_unums from f_ca array */
   nrecvs = -1;
   prevproc = -1;
   for (f_ii = 0; f_ii < f_npoints; f_ii++)
   {
      if (f_ca[f_ii] > -1)
      {
         i = f_ca[f_ii];
         _braid_GetBlockDistProc((gupper+1), nprocs, i, &proc);
         if (proc != prevproc)
         {
            nrecvs++;
            recv_procs[nrecvs] = proc;
            recv_unums[nrecvs] = 0;
            prevproc = proc;
         }
         recv_unums[nrecvs]++;
      }
   }
   nrecvs++;

   requests = _braid_CTAlloc(MPI_Request, (nsends+nrecvs));
   statuses = _braid_CTAlloc(MPI_Status,  (nsends+nrecvs));

   _braid_BufferStatusInit( 1, 0, bstatus );
   _braid_BaseBufSize(core, app,  &max_usize, bstatus); /* max buffer size */
   _braid_NBytesToNReals(max_usize, max_usize);

   /* Post u-vector receives */
   for (m = 0; m < nrecvs; m++)
   {
      unum = recv_unums[m]; /* Number of u-vectors being received */
      recv_size = unum*(1 + max_usize);
      recv_buffers[m] = _braid_CTAlloc(braid_Real, recv_size);
      MPI_Irecv(recv_buffers[m], recv_size, braid_MPI_REAL, recv_procs[m], 5, comm,
                &requests[m]);

#if DEBUG
      proc = recv_procs[m];
      printf("%d %d: 2 recv %d, proc = %d, unum = %d, size = %d\n",
             FRefine_count, myproc, m, proc, unum, recv_size);
#endif
   }
  
   /* Post u-vector sends */
   ii = 0;
   for (m = 0; m < nsends; m++)
   {
      unum = send_unums[m]; /* Number of u-vectors being sent */
      send_size = unum*(1 + max_usize);
      send_buffers[m] = _braid_CTAlloc(braid_Real, send_size);
      send_size = 0; /* Recompute send_size and realloc buffer */
      bptr = send_buffers[m];
      while (unum > 0)
      {
         if (send_ua[ii] != NULL)
         {
            /* Pack u into buffer, adjust size, and put size into buffer */
            buffer = &bptr[1];
            _braid_BaseBufSize(core, app,  &size, bstatus);
            _braid_StatusElt( bstatus, size_buffer ) = size;
            _braid_BaseBufPack(core, app,  send_ua[ii], buffer, bstatus);
            size = _braid_StatusElt(bstatus, size_buffer);
            _braid_BaseFree(core, app,  send_ua[ii]);
            _braid_NBytesToNReals(size, size);
            bptr[0] = (braid_Int) size; /* insert size at the beginning */
            bptr += (1+size);
            send_size += (1+size);
            unum--;
         }
         ii++;
      }
      send_buffers[m] = _braid_TReAlloc(send_buffers[m], braid_Real, send_size);
      MPI_Isend(send_buffers[m], send_size, braid_MPI_REAL, send_procs[m], 5, comm,
                &requests[m + nrecvs]);

#if DEBUG
      unum = send_unums[m];
      proc = send_procs[m];
      printf("%d %d: 2 send %d, proc = %d, unum = %d, size = %d\n",
             FRefine_count, myproc, m, proc, unum, send_size);
#endif
   }

#if DEBUG
   printf("%d %d: 3\n", FRefine_count, myproc);
#endif

   /* Finish communication */
   MPI_Waitall((nsends+nrecvs), requests, statuses);

#if DEBUG
   printf("%d %d: 4\n", FRefine_count, myproc);
#endif

   /* Unpack u-vectors */
   f_ii = 0;
   for (m = 0; m < nrecvs; m++)
   {
      unum = recv_unums[m];
      bptr = recv_buffers[m];
      while (unum > 0)
      {
         if (f_ca[f_ii] > -1)
         {
            /* Unpack buffer into u-vector */
            buffer = &bptr[1];
            _braid_BaseBufUnpack(core, app, buffer, &recv_ua[f_ii], bstatus);
            size = (braid_Int) bptr[0];
            bptr += (1+size);
            unum--;
         }
         f_ii++;
      }
   }

   /* Free up some memory */
   _braid_TFree(send_ua);
   _braid_TFree(send_procs);
   _braid_TFree(send_unums);
   for (m = 0; m < nsends; m++)
   {
      _braid_TFree(send_buffers[m]);
   }
   _braid_TFree(send_buffers);
   _braid_TFree(recv_procs);
   _braid_TFree(recv_unums);
   for (m = 0; m < nrecvs; m++)
   {
      _braid_TFree(recv_buffers[m]);
   }
   _braid_TFree(recv_buffers);
   _braid_TFree(requests);
   _braid_TFree(statuses);
   _braid_TFree(r_ca);
   _braid_TFree(r_ta_alloc);
   _braid_TFree(r_fa);
   _braid_TFree(f_ca);
   {
      braid_Int  level, nlevels = _braid_CoreElt(core, nlevels);
      _braid_TFree(_braid_CoreElt(core, rfactors));
      _braid_TFree(_braid_CoreElt(core, tnorm_a));

      for (level = 0; level < nlevels; level++)
      {
         _braid_GridDestroy(core, grids[level]);
      }
   }

   /*-----------------------------------------------------------------------*/
   /* 5. Build the new fine grid hierarchy, then use recv_ua to populate the
    * initial values on grid level 0.  This is done by integrating values to the
    * next C-point the right.  Because we require that rfactor <= cfactor, each
    * C-point either has a corresponding value or has a value in the F-interval
    * immediately to the left that can be integrated.  Communication from the
    * left processor may still be needed. */

#if DEBUG
   printf("%d %d: 5\n", FRefine_count, myproc);
#endif

   /* Initialize new hierarchy */
   _braid_CoreElt(core, gupper)  = f_gupper;
   _braid_CoreElt(core, nrefine) += 1;
   /*braid_SetCFactor(core,  0, cfactor);*/ /* RDF HACKED TEST */
   _braid_InitHierarchy(core, f_grid, 1);

   /* Initialize communication */
   recv_msg = 0;
   send_msg = 0;
   if (f_first > _braid_NextCPoint(f_ilower, cfactor))
   {
      recv_msg = 1;
   }
   if (f_next > _braid_NextCPoint(f_iupper+1, cfactor))
   {
      send_msg = 1;
   }
   _braid_UCommInitBasic(core, 0, recv_msg, send_msg, 0);

#if DEBUG
   printf("%d %d: 6 recv_msg = %d, send_msg = %d\n",
          FRefine_count, myproc, recv_msg, send_msg);
#endif

   /* Start from the right-most point */
   f_i = f_iupper;
   next = f_next;
   while (f_i >= f_ilower)
   {
      /* Find the next value to the left */
      u = NULL;
      for ( ; f_i >= f_ilower; f_i--)
      {
         f_ii = f_i - f_ilower;
         if (recv_ua[f_ii] != NULL)
         {
            u = recv_ua[f_ii];
            break;
         }
      }
      if ((f_i < f_ilower) && (recv_msg))
      {
         f_i = f_ilower-1; /* receive value from left processor */
         _braid_UGetVector(core, 0, f_i, &u);
      }

      /* Integrate the value if needed and set */
      if (u != NULL)
      {
         f_j = f_i;
         f_ci = _braid_PriorCPoint(next-1, cfactor);
         if (f_i < f_ci)
         {
            /* integrate */
            f_hi = _braid_min(f_ci, f_iupper);
            for ( ; f_j < f_hi; f_j++)
            {
               _braid_USetVector(core, 0, f_j, u, 0);
               _braid_Step(core, 0, f_j+1, NULL, u);
            }
         }
         _braid_USetVector(core, 0, f_j, u, 1);
         next = f_i;
         f_i--;
      }
   }

   /* Free up some memory */
   _braid_TFree(recv_ua);

#if DEBUG
   printf("%d %d: 7\n", FRefine_count, myproc);
   fflush(stdout);
#endif

   _braid_UCommWait(core, 0);

#if DEBUG
   printf("%d %d: 8\n", FRefine_count, myproc);
   fflush(stdout);
   FRefine_count++;
#endif

   *refined_ptr = 1;

   return _braid_error_flag;
}

