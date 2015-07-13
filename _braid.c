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

/** \file _braid.c
 * \brief Source code for developer routines.  See braid.h for more information.
 *
 */

#include "_braid.h"
#include "braid_defs.h"
#include "util.h"

#define DEBUG 0

#if DEBUG
braid_Int  FRefine_count = 0;
#endif

braid_Int _braid_error_flag = 0;
FILE    *_braid_printfile  = NULL;

/*--------------------------------------------------------------------------
 * Macros used below
 *--------------------------------------------------------------------------*/

/* Compute number of reals given some number of bytes (use ceiling) */
#define _braid_NBytesToNReals(nbytes, nreals) \
nreals = nbytes / sizeof(braid_Real) + ((nbytes % sizeof(braid_Real)) != 0)

/*--------------------------------------------------------------------------
 * Returns the index interval for 'proc' in a blocked data distribution
 *--------------------------------------------------------------------------*/

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

/*--------------------------------------------------------------------------
 * Returns the processor that owns 'index' in a blocked data distribution
 * (returns -1 if index is out of range)
 *--------------------------------------------------------------------------*/

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

/*--------------------------------------------------------------------------
 * Returns the index interval for my processor on the finest grid level
 *--------------------------------------------------------------------------*/

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

/*--------------------------------------------------------------------------
 * Returns the processor that owns 'index' on the given grid 'level'
 * (returns -1 if index is out of range)
 *--------------------------------------------------------------------------*/

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

/*--------------------------------------------------------------------------
 * Returns the coarsening factor to use on grid 'level'
 *--------------------------------------------------------------------------*/

braid_Int
_braid_GetCFactor(braid_Core   core,
                  braid_Int    level,
                  braid_Int   *cfactor_ptr)
{
   braid_Int     *cfactors  = _braid_CoreElt(core, cfactors);
   braid_Int      cfdefault = _braid_CoreElt(core, cfdefault);
   braid_Int      cfactor;

   if (cfactors[level] != 0)
   {
      cfactor = cfactors[level];
   }
   else
   {
      cfactor = cfdefault;
   }
   *cfactor_ptr = cfactor;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
_braid_CommRecvInit(braid_Core           core,
                    braid_Int            level,
                    braid_Int            index,
                    braid_Vector        *vector_ptr,
                    _braid_CommHandle  **handle_ptr)
{
   MPI_Comm            comm = _braid_CoreElt(core, comm);
   braid_App           app  = _braid_CoreElt(core, app);
   _braid_CommHandle  *handle = NULL;
   void               *buffer;
   MPI_Request        *requests;
   MPI_Status         *status;
   braid_Int           proc, size, num_requests;

   _braid_GetProc(core, level, index, &proc);
   if (proc > -1)
   {
      handle = _braid_TAlloc(_braid_CommHandle, 1);

      /* Allocate buffer through user routine */
      _braid_CoreFcn(core, bufsize)(app, &size);
      buffer = malloc(size);

      num_requests = 1;
      requests = _braid_CTAlloc(MPI_Request, num_requests);
      status   = _braid_CTAlloc(MPI_Status, num_requests);
      MPI_Irecv(buffer, size, MPI_BYTE, proc, 0, comm, &requests[0]);

      _braid_CommHandleElt(handle, request_type) = 1; /* recv type = 1 */
      _braid_CommHandleElt(handle, num_requests) = num_requests;
      _braid_CommHandleElt(handle, requests)     = requests;
      _braid_CommHandleElt(handle, status)       = status;
      _braid_CommHandleElt(handle, buffer)       = buffer;
      _braid_CommHandleElt(handle, vector_ptr)   = vector_ptr;
   }

   *handle_ptr = handle;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
_braid_CommSendInit(braid_Core           core,
                    braid_Int            level,
                    braid_Int            index,
                    braid_Vector         vector,
                    _braid_CommHandle  **handle_ptr)
{
   MPI_Comm            comm = _braid_CoreElt(core, comm);
   braid_App           app  = _braid_CoreElt(core, app);
   _braid_CommHandle  *handle = NULL;
   void               *buffer;
   MPI_Request        *requests;
   MPI_Status         *status;
   braid_Int           proc, size, num_requests;

   _braid_GetProc(core, level, index+1, &proc);
   if (proc > -1)
   {
      handle = _braid_TAlloc(_braid_CommHandle, 1);

      /* Allocate buffer through user routine */
      _braid_CoreFcn(core, bufsize)(app, &size);
      buffer = malloc(size);
      /* Note that bufpack may return a size smaller than bufsize */
      _braid_CoreFcn(core, bufpack)(app, vector, buffer, &size);

      num_requests = 1;
      requests = _braid_CTAlloc(MPI_Request, num_requests);
      status   = _braid_CTAlloc(MPI_Status, num_requests);
      MPI_Isend(buffer, size, MPI_BYTE, proc, 0, comm, &requests[0]);

      _braid_CommHandleElt(handle, request_type) = 0; /* send type = 0 */
      _braid_CommHandleElt(handle, num_requests) = num_requests;
      _braid_CommHandleElt(handle, requests)     = requests;
      _braid_CommHandleElt(handle, status)       = status;
      _braid_CommHandleElt(handle, buffer)       = buffer;
   }

   *handle_ptr = handle;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
_braid_CommWait(braid_Core          core,
                _braid_CommHandle **handle_ptr)
{
   braid_App           app    = _braid_CoreElt(core, app);
   _braid_CommHandle  *handle = *handle_ptr;

   if (handle != NULL)
   {
      braid_Int      request_type = _braid_CommHandleElt(handle, request_type);
      braid_Int      num_requests = _braid_CommHandleElt(handle, num_requests);
      MPI_Request   *requests     = _braid_CommHandleElt(handle, requests);
      MPI_Status    *status       = _braid_CommHandleElt(handle, status);
      void          *buffer       = _braid_CommHandleElt(handle, buffer);

      MPI_Waitall(num_requests, requests, status);
      
      if (request_type == 1) /* recv type */
      {
         braid_Vector  *vector_ptr = _braid_CommHandleElt(handle, vector_ptr);
         
         _braid_CoreFcn(core, bufunpack)(app, buffer, vector_ptr);
      }
      
      _braid_TFree(requests);
      _braid_TFree(status);
      _braid_TFree(handle);
      _braid_TFree(buffer);

      *handle_ptr = NULL;
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Returns an index into the local u-vector for grid 'level' at point 'index'.
 * If the u-vector is not stored, returns -1.
 *--------------------------------------------------------------------------*/

braid_Int
_braid_UGetIndex(braid_Core   core,
                 braid_Int    level,
                 braid_Int    index,
                 braid_Int   *uindex_ptr)
{
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_Int            ilower      = _braid_GridElt(grids[level], ilower);
   braid_Int            iupper      = _braid_GridElt(grids[level], iupper);
   braid_Int            clower      = _braid_GridElt(grids[level], clower);
   braid_Int            cfactor     = _braid_GridElt(grids[level], cfactor);
   braid_Int            uindex, ic, iclo;

   uindex = -1;
   if ((index >= ilower) && (index <= iupper))
   {
      if ( (_braid_CoreElt(core, storage) == -1) ||
           (level < _braid_CoreElt(core, storage)) )
      {
         if ( _braid_IsCPoint(index, cfactor) )
         {
            _braid_MapFineToCoarse(index, cfactor, ic);
            _braid_MapFineToCoarse(clower, cfactor, iclo);
            uindex = ic-iclo;
         }
      }
      else
      {
         uindex = index-ilower;
      }
   }

   *uindex_ptr = uindex;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Returns a reference to the local u-vector on grid 'level' at point 'index'.
 * If the u-vector is not stored, returns NULL.
 *--------------------------------------------------------------------------*/

braid_Int
_braid_UGetVectorRef(braid_Core     core,
                     braid_Int      level,
                     braid_Int      index,
                     braid_Vector  *u_ptr)
{
   _braid_Grid        **grids = _braid_CoreElt(core, grids);
   braid_Vector        *ua    = _braid_GridElt(grids[level], ua);
   braid_Int            iu;
   braid_Vector         u = NULL;

   _braid_UGetIndex(core, level, index, &iu);
   if (iu > -1)
   {
      u = ua[iu];
   }

   *u_ptr = u;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Stores a reference to the local u-vector on grid 'level' at point 'index'.
 * If the u-vector is not stored, nothing is done.
 *--------------------------------------------------------------------------*/

braid_Int
_braid_USetVectorRef(braid_Core    core,
                     braid_Int     level,
                     braid_Int     index,
                     braid_Vector  u)
{
   _braid_Grid        **grids = _braid_CoreElt(core, grids);
   braid_Vector        *ua    = _braid_GridElt(grids[level], ua);
   braid_Int            iu;

   _braid_UGetIndex(core, level, index, &iu);
   if (iu > -1)
   {
      ua[iu] = u;
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Returns a copy of the u-vector on grid 'level' at point 'index'.  If 'index'
 * is my "receive index" (as set by UCommInit(), for example), the u-vector will
 * be received from a neighbor processor.  If the u-vector is not stored, NULL
 * is returned.
 *--------------------------------------------------------------------------*/

braid_Int
_braid_UGetVector(braid_Core     core,
                  braid_Int      level,
                  braid_Int      index,
                  braid_Vector  *u_ptr)
{
   braid_App            app         = _braid_CoreElt(core, app);
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_Vector        *ua          = _braid_GridElt(grids[level], ua);
   braid_Int            recv_index  = _braid_GridElt(grids[level], recv_index);
   _braid_CommHandle   *recv_handle = _braid_GridElt(grids[level], recv_handle);
   braid_Vector         u = NULL;
   braid_Int            iu;

   if (index == recv_index)
   {
      /* If a recv was initiated, receive u value from neighbor processor */
      if (recv_index > -1)
      {
         _braid_CommWait(core, &recv_handle);
         _braid_GridElt(grids[level], recv_index)  = -1;
         _braid_GridElt(grids[level], recv_handle) = recv_handle;
         u = ua[-1];
      }
   }
   else
   {
      _braid_UGetIndex(core, level, index, &iu);
      if (iu > -1)
      {
         _braid_CoreFcn(core, clone)(app, ua[iu], &u);
      }
   }

   *u_ptr = u;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Stores the u-vector on grid 'level' at point 'index'.  If 'index' is my "send
 * index", a send is initiated to a neighbor processor.  If 'move' is true, the
 * u-vector is moved into core storage instead of copied.  If the u-vector is
 * not stored, nothing is done.
 *--------------------------------------------------------------------------*/

braid_Int
_braid_USetVector(braid_Core    core,
                  braid_Int     level,
                  braid_Int     index,
                  braid_Vector  u,
                  braid_Int     move)
{
   braid_App            app         = _braid_CoreElt(core, app);
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_Vector        *ua          = _braid_GridElt(grids[level], ua);
   braid_Int            send_index  = _braid_GridElt(grids[level], send_index);
   _braid_CommHandle   *send_handle = _braid_GridElt(grids[level], send_handle);
   braid_Int            iu;

   if (index == send_index)
   {
      /* Post send to neighbor processor */
      _braid_CommSendInit(core, level, index, u, &send_handle);
      _braid_GridElt(grids[level], send_index)  = -1;
      _braid_GridElt(grids[level], send_handle) = send_handle;
   }

   _braid_UGetIndex(core, level, index, &iu);
   if (iu > -1)
   {
      if (ua[iu] != NULL)
      {
         _braid_CoreFcn(core, free)(app, ua[iu]);
      }
      if (move)
      {
         ua[iu] = u;                                   /* move the vector */
      }
      else
      {
         _braid_CoreFcn(core, clone)(app, u, &ua[iu]); /* copy the vector */
      }
   }
   else if (move)
   {
      _braid_CoreFcn(core, free)(app, u);              /* free the vector */
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Basic communication (from the left, to the right)
 *--------------------------------------------------------------------------*/

braid_Int
_braid_UCommInitBasic(braid_Core  core,
                      braid_Int   level,
                      braid_Int   recv_msg,
                      braid_Int   send_msg,
                      braid_Int   send_now)
{
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_Int            ilower      = _braid_GridElt(grids[level], ilower);
   braid_Int            iupper      = _braid_GridElt(grids[level], iupper);
   braid_Vector        *ua          = _braid_GridElt(grids[level], ua);
   braid_Int            recv_index  = -1;
   braid_Int            send_index  = -1;
   _braid_CommHandle   *recv_handle = NULL;
   _braid_CommHandle   *send_handle = NULL;
   braid_Int            iu;

   /* Post receive */
   if (recv_msg)
   {
      recv_index = ilower-1;
      _braid_CommRecvInit(core, level, recv_index, &ua[-1], &recv_handle);
   }

   /* Post send */
   if (send_msg)
   {
      send_index = iupper;
      if (send_now)
      {
         _braid_UGetIndex(core, level, send_index, &iu);
         _braid_CommSendInit(core, level, send_index, ua[iu], &send_handle);
         send_index = -1;
      }
   }

   _braid_GridElt(grids[level], recv_index)  = recv_index ;
   _braid_GridElt(grids[level], send_index)  = send_index ;
   _braid_GridElt(grids[level], recv_handle) = recv_handle;
   _braid_GridElt(grids[level], send_handle) = send_handle;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Working on all intervals
 *--------------------------------------------------------------------------*/

braid_Int
_braid_UCommInit(braid_Core  core,
                 braid_Int   level)
{
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_Int            ilower      = _braid_GridElt(grids[level], ilower);
   braid_Int            iupper      = _braid_GridElt(grids[level], iupper);
   braid_Int            cfactor     = _braid_GridElt(grids[level], cfactor);
   braid_Vector        *ua          = _braid_GridElt(grids[level], ua);
   braid_Int            recv_index  = -1;
   braid_Int            send_index  = -1;
   _braid_CommHandle   *recv_handle = NULL;
   _braid_CommHandle   *send_handle = NULL;
   braid_Int            iu;

   if (ilower <= iupper)
   {
      /* Post receive */
      _braid_CommRecvInit(core, level, ilower-1, &ua[-1], &recv_handle);
      recv_index = ilower-1;
      
      /* Only post send if iupper is a C-point, otherwise compute and send later */
      if ( _braid_IsCPoint(iupper, cfactor) )
      {
         _braid_UGetIndex(core, level, iupper, &iu);
         _braid_CommSendInit(core, level, iupper, ua[iu], &send_handle);
         send_index = -1;
      }
      else
      {
         send_index = iupper;
      }
   }

   _braid_GridElt(grids[level], recv_index)  = recv_index ;
   _braid_GridElt(grids[level], send_index)  = send_index ;
   _braid_GridElt(grids[level], recv_handle) = recv_handle;
   _braid_GridElt(grids[level], send_handle) = send_handle;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Working only on F-pt intervals
 *--------------------------------------------------------------------------*/

braid_Int
_braid_UCommInitF(braid_Core  core,
                  braid_Int   level)
{
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_Int            ilower      = _braid_GridElt(grids[level], ilower);
   braid_Int            iupper      = _braid_GridElt(grids[level], iupper);
   braid_Int            cfactor     = _braid_GridElt(grids[level], cfactor);
   braid_Vector        *ua          = _braid_GridElt(grids[level], ua);
   braid_Int            recv_index  = -1;
   braid_Int            send_index  = -1;
   _braid_CommHandle   *recv_handle = NULL;
   _braid_CommHandle   *send_handle = NULL;
   braid_Int            iu;

   if (ilower <= iupper)
   {
      /* Only post receive if ilower is an F-point */
      if ( _braid_IsFPoint(ilower, cfactor) )
      {
         _braid_CommRecvInit(core, level, ilower-1, &ua[-1], &recv_handle);
         recv_index = ilower-1;
      }

      /* Only post send if iupper is a C-point, otherwise compute and send later */
      if ( _braid_IsCPoint(iupper, cfactor) )
      {
         _braid_UGetIndex(core, level, iupper, &iu);
         _braid_CommSendInit(core, level, iupper, ua[iu], &send_handle);
         send_index = -1;
      }
      else if ( _braid_IsFPoint(iupper+1, cfactor) )
      {
         send_index = iupper;
      }
   }

   _braid_GridElt(grids[level], recv_index)  = recv_index ;
   _braid_GridElt(grids[level], send_index)  = send_index ;
   _braid_GridElt(grids[level], recv_handle) = recv_handle;
   _braid_GridElt(grids[level], send_handle) = send_handle;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Finish up communication
 *--------------------------------------------------------------------------*/

braid_Int
_braid_UCommWait(braid_Core  core,
                 braid_Int   level)
{
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   _braid_CommHandle   *recv_handle = _braid_GridElt(grids[level], recv_handle);
   _braid_CommHandle   *send_handle = _braid_GridElt(grids[level], send_handle);

   _braid_CommWait(core, &recv_handle);
   _braid_CommWait(core, &send_handle);
   _braid_GridElt(grids[level], recv_index)  = -1;
   _braid_GridElt(grids[level], send_index)  = -1;
   _braid_GridElt(grids[level], recv_handle) = recv_handle;
   _braid_GridElt(grids[level], send_handle) = send_handle;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
_braid_GetInterval(braid_Core   core,
                   braid_Int    level,
                   braid_Int    interval_index,
                   braid_Int   *flo_ptr,
                   braid_Int   *fhi_ptr,
                   braid_Int   *ci_ptr)
{
   _braid_Grid  **grids   = _braid_CoreElt(core, grids);
   braid_Int      ilower  = _braid_GridElt(grids[level], ilower);
   braid_Int      iupper  = _braid_GridElt(grids[level], iupper);
   braid_Int      clower  = _braid_GridElt(grids[level], clower);
   braid_Int      cupper  = _braid_GridElt(grids[level], cupper);
   braid_Int      cfactor = _braid_GridElt(grids[level], cfactor);
   braid_Int      flo, fhi, ci;

   flo = ilower;
   fhi = iupper;
   ci  = -1;

   if ( _braid_IsCPoint(clower, cfactor) )
   {
      flo = clower + (interval_index-1)*cfactor + 1;
      fhi = clower + (interval_index  )*cfactor - 1;
      if (flo < ilower)
      {
         flo = ilower;
      }
      if (fhi > iupper)
      {
         fhi = iupper;
      }

      ci = clower + interval_index*cfactor;
      if (ci > cupper)
      {
         ci = -1;  /* return -1 if no C-points */
      }
   }

   *flo_ptr = flo;
   *fhi_ptr = fhi;
   *ci_ptr  = ci;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
_braid_AccessVector(braid_Core          core,
                    braid_AccessStatus  status,
                    braid_Vector        u)
{
   braid_App      app    = _braid_CoreElt(core, app);

   _braid_CoreFcn(core, access)(app, u, status);

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Get an initial guess for ustop to use in the step routine (implicit schemes)
 *--------------------------------------------------------------------------*/

braid_Int
_braid_GetUInit(braid_Core     core,
                braid_Int      level,
                braid_Int      index,
                braid_Vector   u,
                braid_Vector  *ustop_ptr)
{
   _braid_Grid    **grids    = _braid_CoreElt(core, grids);
   braid_Int        ilower   = _braid_GridElt(grids[level], ilower);
   braid_Vector    *va       = _braid_GridElt(grids[level], va);

   braid_Vector     ustop = *ustop_ptr;
   braid_Int        ii;

   ii = index-ilower;

   _braid_UGetVectorRef(core, level, index, &ustop);

   /* Run in compatibility mode, mimic the original Braid algorithm */
   if ( _braid_CoreElt(core, storage) == -2 )
   {
      /* Here we must always approximate ustop by u at tstart. */
      ustop = u;
   }

   /* User-provided residual routine, store u-vectors only at C-points */
   else if ( (_braid_CoreElt(core, storage) == -1) ||
             (level < _braid_CoreElt(core, storage)) )
   {
      if (ustop == NULL)
      {
         /* If the u-vector is not stored, use something else for ustop */
         if (level == 0)
         {
            /* On the fine grid, approximate ustop by u */
            ustop = u;
         }
         else if (va[ii] != NULL)
         {
            /* On coarse grids, approximate ustop by the restricted fine value.
             * This ensures a fixed-point iteration. */
            ustop = va[ii];
         }
      }
   }

   /* Utilize storage of all u-vectors (F and C points) at this level */
   else
   {
      if (ustop == NULL)
      {
         /* Approximate ustop by u */
         ustop = u;
      }
   }

   *ustop_ptr = ustop;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Integrate one time step
 *--------------------------------------------------------------------------*/

braid_Int
_braid_Step(braid_Core     core,
            braid_Int      level,
            braid_Int      index,
            braid_Vector   ustop,
            braid_Vector   u)
{
   braid_App        app      = _braid_CoreElt(core, app);
   braid_Real       tol      = _braid_CoreElt(core, tol);
   braid_Int        iter     = _braid_CoreElt(core, niter);
   braid_Int       *rfactors = _braid_CoreElt(core, rfactors);
   _braid_Grid    **grids    = _braid_CoreElt(core, grids);
   braid_StepStatus status   = _braid_CoreElt(core, sstatus);
   braid_Int        nrefine  = _braid_CoreElt(core, nrefine);
   braid_Int        ilower   = _braid_GridElt(grids[level], ilower);
   braid_Real      *ta       = _braid_GridElt(grids[level], ta);
   braid_Vector    *fa       = _braid_GridElt(grids[level], fa);

   braid_Int        ii;

   ii = index-ilower;
   _braid_StepStatusInit(ta[ii-1], ta[ii], tol, iter, level, nrefine, status);

   /* If ustop is set to NULL, use a default approach for setting it */
   if (ustop == NULL)
   {
      _braid_GetUInit(core, level, index, u, &ustop);
   }

   if (level == 0)
   {
      _braid_CoreFcn(core, step)(app, ustop, NULL, u, status);
      rfactors[ii] = _braid_StatusElt(status, rfactor);
   }
   else
   {
      if ( _braid_CoreElt(core, residual) == NULL )
      {
         _braid_CoreFcn(core, step)(app, ustop, NULL, u, status);
         _braid_CoreFcn(core, sum)(app, 1.0, fa[ii], 1.0, u);
      }
      else
      {
         _braid_CoreFcn(core, step)(app, ustop, fa[ii], u, status);
      }
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Compute residual
 *--------------------------------------------------------------------------*/

braid_Int
_braid_Residual(braid_Core     core,
                braid_Int      level,
                braid_Int      index,
                braid_Vector   ustop,
                braid_Vector   r)
{
   braid_App        app      = _braid_CoreElt(core, app);
   braid_Real       tol      = _braid_CoreElt(core, tol);
   braid_Int        iter     = _braid_CoreElt(core, niter);
   _braid_Grid    **grids    = _braid_CoreElt(core, grids);
   braid_StepStatus status   = _braid_CoreElt(core, sstatus);
   braid_Int        nrefine  = _braid_CoreElt(core, nrefine);
   braid_Int        ilower   = _braid_GridElt(grids[level], ilower);
   braid_Real      *ta       = _braid_GridElt(grids[level], ta);

   braid_Vector     rstop;
   braid_Int        ii;

   ii = index-ilower;
   _braid_StepStatusInit(ta[ii-1], ta[ii], tol, iter, level, nrefine, status);
   if ( _braid_CoreElt(core, residual) == NULL )
   {
      /* By default: r = ustop - \Phi(ustart)*/
      _braid_GetUInit(core, level, index, r, &rstop);
      _braid_CoreFcn(core, step)(app, rstop, NULL, r, status);
      _braid_CoreFcn(core, sum)(app, 1.0, ustop, -1.0, r);
   }
   else
   {
      /* Call the user's residual routine */
      _braid_CoreFcn(core, residual)(app, ustop, r, status);
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Compute FAS residual = f - residual
 *--------------------------------------------------------------------------*/

braid_Int
_braid_FASResidual(braid_Core     core,
                   braid_Int      level,
                   braid_Int      index,
                   braid_Vector   ustop,
                   braid_Vector   r)
{
   braid_App        app    = _braid_CoreElt(core, app);
   _braid_Grid    **grids  = _braid_CoreElt(core, grids);
   braid_Int        ilower = _braid_GridElt(grids[level], ilower);
   braid_Vector    *fa     = _braid_GridElt(grids[level], fa);

   braid_Int        ii;

   _braid_Residual(core, level, index, ustop, r);
   if (level == 0)
   {
      _braid_CoreFcn(core, sum)(app, 0.0, r, -1.0, r);
   }
   else
   {
      ii = index-ilower;
      _braid_CoreFcn(core, sum)(app, 1.0, fa[ii], -1.0, r);
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Coarsen in space
 *--------------------------------------------------------------------------*/

braid_Int
_braid_Coarsen(braid_Core     core,
               braid_Int      level,    /* coarse level */
               braid_Int      f_index,  /* fine index */
               braid_Int      c_index,  /* coarse index */
               braid_Vector   fvector,
               braid_Vector  *cvector)
{
   braid_App      app             = _braid_CoreElt(core, app);
   _braid_Grid  **grids           = _braid_CoreElt(core, grids);
   braid_CoarsenRefStatus cstatus = _braid_CoreElt(core, cstatus);
   braid_Int      nrefine         = _braid_CoreElt(core, nrefine);
   braid_Int      c_ilower        = _braid_GridElt(grids[level], ilower);
   braid_Int      f_ilower        = _braid_GridElt(grids[level-1], ilower);
   braid_Real    *c_ta            = _braid_GridElt(grids[level], ta);
   braid_Real    *f_ta            = _braid_GridElt(grids[level-1], ta);

   braid_Int      c_ii = c_index-c_ilower;
   braid_Int      f_ii = f_index-f_ilower;
   
   if ( _braid_CoreElt(core, coarsen) == NULL )
   {
      /* No spatial coarsening needed, just clone the fine vector.*/
      _braid_CoreFcn(core, clone)(app, fvector, cvector);
   }
   else
   {
      /* Call the user's coarsening routine */
      _braid_CoarsenRefStatusInit(f_ta[f_ii], f_ta[f_ii-1], f_ta[f_ii+1], 
                                  c_ta[c_ii-1], c_ta[c_ii+1],
                                  level-1, nrefine, cstatus);
      _braid_CoreFcn(core, coarsen)(app, fvector, cvector, cstatus);
   }
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Refine in space (basic routine)
 *--------------------------------------------------------------------------*/

braid_Int
_braid_RefineBasic(braid_Core     core,
                   braid_Int      level,    /* fine level */
                   braid_Real    *f_ta,     /* pointer into fine time array */
                   braid_Real    *c_ta,     /* pointer into coarse time array */
                   braid_Vector   cvector,
                   braid_Vector  *fvector)
{
   braid_App              app     = _braid_CoreElt(core, app);
   braid_CoarsenRefStatus cstatus = _braid_CoreElt(core, cstatus);
   braid_Int              nrefine = _braid_CoreElt(core, nrefine);

   if ( _braid_CoreElt(core, coarsen) == NULL )
   {
      /* No spatial refinement needed, just clone the fine vector.*/
      _braid_CoreFcn(core, clone)(app, cvector, fvector);
   }
   else
   {
      /* Call the user's refinement routine */
      _braid_CoarsenRefStatusInit(f_ta[0], f_ta[-1], f_ta[+1], c_ta[-1], c_ta[+1],
                                  level, nrefine, cstatus);
      _braid_CoreFcn(core, refine)(app, cvector, fvector, cstatus);
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Refine in space
 *--------------------------------------------------------------------------*/

braid_Int
_braid_Refine(braid_Core     core,
              braid_Int      level,    /* fine level */
              braid_Int      f_index,  /* fine index */
              braid_Int      c_index,  /* coarse index */
              braid_Vector   cvector,
              braid_Vector  *fvector)
{
   _braid_Grid  **grids           = _braid_CoreElt(core, grids);
   braid_Int      f_ilower        = _braid_GridElt(grids[level], ilower);
   braid_Int      c_ilower        = _braid_GridElt(grids[level+1], ilower);
   braid_Real    *f_ta            = _braid_GridElt(grids[level], ta);
   braid_Real    *c_ta            = _braid_GridElt(grids[level+1], ta);

   braid_Int      c_ii = c_index-c_ilower;
   braid_Int      f_ii = f_index-f_ilower;

   _braid_RefineBasic(core, level, &f_ta[f_ii], &c_ta[c_ii], cvector, fvector);

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Create a new grid object
 *--------------------------------------------------------------------------*/

braid_Int
_braid_GridInit(braid_Core     core,
                braid_Int      level,
                braid_Int      ilower,
                braid_Int      iupper,
                _braid_Grid  **grid_ptr)
{
   _braid_Grid   *grid;
   braid_Real    *ta;

   grid = _braid_CTAlloc(_braid_Grid, 1);
   
   _braid_GridElt(grid, level)  = level;
   _braid_GridElt(grid, ilower) = ilower;
   _braid_GridElt(grid, iupper) = iupper;
   _braid_GridElt(grid, recv_index) = -1;
   _braid_GridElt(grid, send_index) = -1;
   
   /* Store each processor's time slice, plus one time value to the left 
    * and to the right */
   ta = _braid_CTAlloc(braid_Real, iupper-ilower+3);
   _braid_GridElt(grid, ta_alloc) = ta;
   _braid_GridElt(grid, ta)       = ta+1;  /* shift */
   
   *grid_ptr = grid;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
_braid_GridClean(braid_Core    core,
                 _braid_Grid  *grid)
{
   braid_App      app      = _braid_CoreElt(core, app);
   braid_Int      ilower   = _braid_GridElt(grid, ilower);
   braid_Int      iupper   = _braid_GridElt(grid, iupper);
   braid_Int      nupoints = _braid_GridElt(grid, nupoints);
   braid_Vector  *ua       = _braid_GridElt(grid, ua);
   braid_Vector  *va       = _braid_GridElt(grid, va);
   braid_Vector  *fa       = _braid_GridElt(grid, fa);
   braid_Vector  *ua_alloc = _braid_GridElt(grid, ua_alloc);
   braid_Vector  *va_alloc = _braid_GridElt(grid, va_alloc);
   braid_Vector  *fa_alloc = _braid_GridElt(grid, fa_alloc);
   
   braid_Int      ii;

   if (ua_alloc)
   {
      for (ii = 0; ii < nupoints; ii++)
      {
         if (ua[ii] != NULL)
         {
            _braid_CoreFcn(core, free)(app, ua[ii]);
            ua[ii] = NULL;
         }
      }
   }
   if (va_alloc)
   {
      for (ii = -1; ii <= (iupper-ilower); ii++)
      {
         if (va[ii] != NULL)
         {
            _braid_CoreFcn(core, free)(app, va[ii]);
            va[ii] = NULL;
         }
      }
   }
   if (fa_alloc)
   {
      for (ii = -1; ii <= (iupper-ilower); ii++)
      {
         if (fa[ii] != NULL)
         {
            _braid_CoreFcn(core, free)(app, fa[ii]);
            fa[ii] = NULL;
         }
      }
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
_braid_GridDestroy(braid_Core    core,
                   _braid_Grid  *grid)
{
   if (grid)
   {
      braid_Vector  *ua_alloc = _braid_GridElt(grid, ua_alloc);
      braid_Real    *ta_alloc = _braid_GridElt(grid, ta_alloc);
      braid_Vector  *va_alloc = _braid_GridElt(grid, va_alloc);
      braid_Vector  *fa_alloc = _braid_GridElt(grid, fa_alloc);

      _braid_GridClean(core, grid);

      if (ua_alloc)
      {
         _braid_TFree(ua_alloc);
      }
      if (ta_alloc)
      {
         _braid_TFree(ta_alloc);
      }
      if (va_alloc)
      {
         _braid_TFree(va_alloc);
      }
      if (fa_alloc)
      {
         _braid_TFree(fa_alloc);
      }

      _braid_TFree(grid);
   }
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Set initial guess at C-points
 *--------------------------------------------------------------------------*/

braid_Int
_braid_InitGuess(braid_Core  core,
                 braid_Int   level)
{
   braid_App      app     = _braid_CoreElt(core, app);
   _braid_Grid  **grids   = _braid_CoreElt(core, grids);
   braid_Int      ilower  = _braid_GridElt(grids[level], ilower);
   braid_Int      iupper  = _braid_GridElt(grids[level], iupper);
   braid_Int      clower  = _braid_GridElt(grids[level], clower);
   braid_Int      cupper  = _braid_GridElt(grids[level], cupper);
   braid_Int      cfactor = _braid_GridElt(grids[level], cfactor);
   braid_Real    *ta      = _braid_GridElt(grids[level], ta);
   braid_Vector  *va      = _braid_GridElt(grids[level], va);

   braid_Vector   u;
   braid_Int      i, iu;

   if (level == 0)
   {
      /* Only need to initialize the C-points on the finest grid */
      for (i = clower; i <= cupper; i += cfactor)
      {
         _braid_CoreFcn(core, init)(app, ta[i-ilower], &u);
         _braid_USetVectorRef(core, level, i, u);
      }
   }
   else
   {
      for (i = ilower; i <= iupper; i++)
      {
         _braid_UGetIndex(core, level, i, &iu);
         if (iu > -1)
         {
            _braid_CoreFcn(core, clone)(app, va[i-ilower], &u);
            _braid_USetVectorRef(core, level, i, u);
         }
      }
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
_braid_GetFullResidual(braid_Core  core,
                       braid_Int   level,
                       braid_Real *return_norm)
{
   MPI_Comm           comm        = _braid_CoreElt(core, comm);
   braid_App          app         = _braid_CoreElt(core, app);
   braid_Real         tol         = _braid_CoreElt(core, tol);
   braid_Int          iter        = _braid_CoreElt(core, niter);
   _braid_Grid      **grids       = _braid_CoreElt(core, grids);
   braid_StepStatus   status      = _braid_CoreElt(core, sstatus);
   braid_Int          nrefine     = _braid_CoreElt(core, nrefine);
   braid_Int          ncpoints    = _braid_GridElt(grids[level], ncpoints);
   braid_Int          tnorm       = _braid_CoreElt(core, tnorm);
   braid_Real        *ta          = _braid_GridElt(grids[level], ta);
   braid_Int          ilower      = _braid_GridElt(grids[level], ilower);
   _braid_CommHandle *send_handle;
   braid_Int          send_index;

   braid_Int        flo, fhi, fi, ci, ii, interval;
   braid_Real       rnorm_temp, rnorm = 0, global_rnorm = 0;
   braid_Vector     u, r;

   _braid_UCommInit(core, level);

   /* Start from the right-most interval. */
   for (interval = ncpoints; interval > -1; interval--)
   {
      _braid_GetInterval(core, level, interval, &flo, &fhi, &ci);
      if (flo <= fhi)
      {
         _braid_UGetVector(core, level, flo-1, &u);
      }
      else if (ci > 0)
      {
         _braid_UGetVector(core, level, ci-1, &u);
      }

      /* Generate F-points and get residual. */
      for (fi = flo; fi <= fhi; fi++)
      {
         _braid_CoreFcn(core, clone)(app, u, &r);
         _braid_Step(core, level, fi, NULL, u);

         /* Update local processor norm. */
         ii = fi-ilower;
         _braid_StepStatusInit(ta[ii-1], ta[ii], tol, iter, level, nrefine, status);
         _braid_CoreFcn(core, globresidual)(app, u, r, status);
         _braid_CoreFcn(core, spatialnorm)(app, r, &rnorm_temp); 
         if(tnorm == 1)       /* one-norm */ 
         {  
            rnorm += rnorm_temp;
         }
         else if(tnorm == 3)  /* inf-norm */
         {  
            rnorm = (((rnorm_temp) > (rnorm)) ? (rnorm_temp) : (rnorm));
         }
         else                 /* default two-norm */
         {  
            rnorm += (rnorm_temp*rnorm_temp);
         }

         /* Communicate w/ neighbor nodes. */
         send_handle = _braid_GridElt(grids[level], send_handle);
         send_index  = _braid_GridElt(grids[level], send_index);
         if (fi == send_index)
         {
            /* Post send to neighbor processor */
            _braid_CommSendInit(core, level, fi, u, &send_handle);
            _braid_GridElt(grids[level], send_index)  = -1;
            _braid_GridElt(grids[level], send_handle) = send_handle;
         }
         _braid_CoreFcn(core, free)(app, r);
      }
      /* Residual from C-point. */
      if (ci > 0)
      {
         /* Update local processor norm. */
         ii = ci-ilower;
         _braid_StepStatusInit(ta[ii-1], ta[ii], tol, iter, level, nrefine, status);
         _braid_UGetVector(core, level, ci, &r);
         _braid_CoreFcn(core, globresidual)(app, r, u, status);
         _braid_CoreFcn(core, spatialnorm)(app, u, &rnorm_temp);

         if(tnorm == 1)       /* one-norm */ 
         {  
            rnorm += rnorm_temp;
         }
         else if(tnorm == 3)  /* inf-norm */
         {  
            rnorm = (((rnorm_temp) > (rnorm)) ? (rnorm_temp) : (rnorm));
         }
         else                 /* default two-norm */
         {  
            rnorm += (rnorm_temp*rnorm_temp);
         }
         _braid_CoreFcn(core, free)(app, r);
      }

      if ((flo <= fhi) || (ci > 0))
      {
         _braid_CoreFcn(core, free)(app, u);
      }
   }

   _braid_UCommWait(core, level);

   /* Compute global residual norm. */
   if(tnorm == 1)       /* one-norm reduction */
   {  
      MPI_Allreduce(&rnorm, &global_rnorm, 1, MPI_DOUBLE, MPI_SUM, comm);
   }
   else if(tnorm == 3)  /* inf-norm reduction */
   {  
      MPI_Allreduce(&rnorm, &global_rnorm, 1, MPI_DOUBLE, MPI_MAX, comm);
   }
   else                 /* default two-norm reduction */
   {  
      MPI_Allreduce(&rnorm, &global_rnorm, 1, MPI_DOUBLE, MPI_SUM, comm);
      global_rnorm = sqrt(global_rnorm);
   }

   *return_norm = global_rnorm;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Do nu sweeps of F-then-C relaxation
 *--------------------------------------------------------------------------*/

braid_Int
_braid_FCRelax(braid_Core  core,
               braid_Int   level)
{
   braid_App       app      = _braid_CoreElt(core, app);
   braid_Int      *nrels    = _braid_CoreElt(core, nrels);
   _braid_Grid   **grids    = _braid_CoreElt(core, grids);
   braid_Int       ncpoints = _braid_GridElt(grids[level], ncpoints);

   braid_Vector    u;
   braid_Int       flo, fhi, fi, ci;
   braid_Int       nu, nrelax, interval;

   nrelax  = nrels[level];

   for (nu = 0; nu < nrelax; nu++)
   {
      _braid_UCommInit(core, level);

      /* Start from the right-most interval */
      for (interval = ncpoints; interval > -1; interval--)
      {
         _braid_GetInterval(core, level, interval, &flo, &fhi, &ci);

         if (flo <= fhi)
         {
            _braid_UGetVector(core, level, flo-1, &u);
         }
         else if (ci > 0)
         {
            _braid_UGetVector(core, level, ci-1, &u);
         }

         /* F-relaxation */
         for (fi = flo; fi <= fhi; fi++)
         {
            _braid_Step(core, level, fi, NULL, u);
            _braid_USetVector(core, level, fi, u, 0);
         }

         /* C-relaxation */
         if (ci > 0)
         {
            _braid_Step(core, level, ci, NULL, u);
            _braid_USetVector(core, level, ci, u, 1);
         }

         /* if ((flo <= fhi) && (interval == ncpoints)) */
         if ((flo <= fhi) && !(ci > 0))
         {
            _braid_CoreFcn(core, free)(app, u);
         }
      }
      _braid_UCommWait(core, level);
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * F-Relax on level and restrict to level+1
 *--------------------------------------------------------------------------*/

braid_Int
_braid_FRestrict(braid_Core   core,
                 braid_Int    level,
                 braid_Real  *rnorm_ptr)
{
   MPI_Comm             comm        = _braid_CoreElt(core, comm);
   braid_App            app         = _braid_CoreElt(core, app);
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_AccessStatus   astatus     = _braid_CoreElt(core, astatus);
   braid_Int            iter        = _braid_CoreElt(core, niter);
   braid_Int            print_level = _braid_CoreElt(core, print_level);
   braid_Int            access_level= _braid_CoreElt(core, access_level);
   braid_Int            tnorm       = _braid_CoreElt(core, tnorm);
   braid_Real          *tnorm_a     = _braid_CoreElt(core, tnorm_a);
   braid_Int            nrefine     = _braid_CoreElt(core, nrefine);
   braid_Int            cfactor     = _braid_GridElt(grids[level], cfactor);
   braid_Int            ncpoints    = _braid_GridElt(grids[level], ncpoints);
   braid_Real          *ta          = _braid_GridElt(grids[level], ta);
   braid_Int            f_ilower    = _braid_GridElt(grids[level], ilower);
   _braid_CommHandle   *recv_handle = NULL;
   _braid_CommHandle   *send_handle = NULL;
   braid_Real           old_rnorm   = *rnorm_ptr;

   braid_Int            c_level, c_ilower, c_iupper, c_index, c_i, c_ii;
   braid_Vector         c_u, *c_va, *c_fa;

   braid_Vector         u, r;
   braid_Int            interval, flo, fhi, fi, ci;
   braid_Real           rnorm, grnorm, rnorm_temp;

   c_level  = level+1;
   c_ilower = _braid_GridElt(grids[c_level], ilower);
   c_iupper = _braid_GridElt(grids[c_level], iupper);
   c_va     = _braid_GridElt(grids[c_level], va);
   c_fa     = _braid_GridElt(grids[c_level], fa);

   rnorm = 0.0;

   _braid_UCommInit(core, level);

   /* Start from the right-most interval.
    * 
    * Do an F-relax and then a C-relax.  These relaxations are needed to compute
    * the residual, which is needed for the coarse-grid right-hand-side and for
    * convergence checking on the finest grid.  This loop updates va and fa. */
   for (interval = ncpoints; interval > -1; interval--)
   {
      _braid_GetInterval(core, level, interval, &flo, &fhi, &ci);

      if (flo <= fhi)
      {
         _braid_UGetVector(core, level, flo-1, &r);
      }
      else if (ci > 0)
      {
         _braid_UGetVector(core, level, ci-1, &r);
      }

      /* F-relaxation */
      for (fi = flo; fi <= fhi; fi++)
      {
         _braid_Step(core, level, fi, NULL, r);
         _braid_USetVector(core, level, fi, r, 0);
         
         /* Allow user to process current vector, note that r here is
          * temporarily holding the state vector */
         if( (access_level >= 2) && (level == 0) )
         {
            _braid_AccessStatusInit(ta[fi-f_ilower], old_rnorm, iter,
                                    level, nrefine, 0, 0, astatus);
            _braid_AccessVector(core, astatus, r);
         }
      }

      /* Allow user to process current C-point */
      if( (access_level>= 2) && (level == 0) && (ci > -1) )
      {
         _braid_UGetVectorRef(core, level, ci, &u);
         _braid_AccessStatusInit(ta[ci-f_ilower], old_rnorm, iter,
                                 level, nrefine, 0, 0, astatus);
         _braid_AccessVector(core, astatus, u);
      }
      
      /* Compute residual and restrict */
      if (ci > 0)
      {
         /* Compute FAS residual */
         _braid_UGetVectorRef(core, level, ci, &u);
         _braid_FASResidual(core, level, ci, u, r);

         /* Compute rnorm (only on level 0) */
         if (level == 0)
         {
            _braid_CoreFcn(core, spatialnorm)(app, r, &rnorm_temp);
            if(tnorm == 1) 
            {  
               rnorm += rnorm_temp;               /* one-norm combination */ 
            }
            else if(tnorm == 2)
            {  
               rnorm += (rnorm_temp*rnorm_temp);  /* two-norm combination */
            }
            else if(tnorm == 3)
            {  
               tnorm_a[interval] = rnorm_temp;    /* inf-norm combination */
            }
            
            /* If debug printing, print out rnorm_temp for this interval. rnorm_temp
             * should show the serial propagation of the exact solution */
            if (print_level >= 2)
            {
               _braid_printf("  Braid:  time step: %6d, rnorm: %1.2e\n", fhi, rnorm_temp);
            }
         }

         /* Restrict u and residual, coarsening in space if needed */
         _braid_MapFineToCoarse(ci, cfactor, c_index);
         _braid_Coarsen(core, c_level, ci, c_index, u, &c_va[c_index-c_ilower]);
         _braid_Coarsen(core, c_level, ci, c_index, r, &c_fa[c_index-c_ilower]);
      }
      else if (ci == 0)
      {
         /* Restrict initial condition, coarsening in space if needed */
         _braid_UGetVectorRef(core, level, 0, &u);
         _braid_Coarsen(core, c_level, 0, 0, u, &c_va[0]);
      }

      if ((flo <= fhi) || (ci > 0))
      {
         _braid_CoreFcn(core, free)(app, r);
      }
   }
   _braid_UCommWait(core, level);

   /* Now apply coarse residual to update fa values */

   /* Set initial guess on coarse level */
   _braid_InitGuess(core, c_level);

   /* Initialize update of c_va[-1] boundary */
   if (c_ilower <= c_iupper)
   {
      _braid_CommRecvInit(core, c_level, c_ilower-1, &c_va[-1],
                         &recv_handle);
      _braid_CommSendInit(core, c_level, c_iupper, c_va[c_iupper-c_ilower],
                         &send_handle);
   }

   /* Start with rightmost point */
   for (c_i = c_iupper; c_i >= c_ilower; c_i--)
   {
      if (c_i > 0)
      {
         c_ii = c_i - c_ilower;
         if (c_ii == 0)
         {
            /* Finalize update of c_va[-1] */
            _braid_CommWait(core, &recv_handle);
         }
         _braid_CoreFcn(core, clone)(app, c_va[c_ii-1], &c_u);
         _braid_Residual(core, c_level, c_i, c_va[c_ii], c_u);
         _braid_CoreFcn(core, sum)(app, 1.0, c_u, 1.0, c_fa[c_ii]);
         _braid_CoreFcn(core, free)(app, c_u);
      }
   }
   _braid_CommWait(core, &send_handle);

   /* Compute rnorm (only on level 0) */
   if (level == 0)
   {
      if(tnorm == 1)          /* one-norm reduction */
      {  
         MPI_Allreduce(&rnorm, &grnorm, 1, braid_MPI_REAL, MPI_SUM, comm);
      }
      else if(tnorm == 3)     /* inf-norm reduction */
      {  
         _braid_Max(tnorm_a, ncpoints, &rnorm); 
         MPI_Allreduce(&rnorm, &grnorm, 1, braid_MPI_REAL, MPI_MAX, comm);
      }
      else                    /* default two-norm reduction */
      {  
         MPI_Allreduce(&rnorm, &grnorm, 1, braid_MPI_REAL, MPI_SUM, comm);
         grnorm = sqrt(grnorm);
      }

      *rnorm_ptr = grnorm;
   }
   
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * F-Relax on level and interpolate to level-1
 *--------------------------------------------------------------------------*/

braid_Int
_braid_FInterp(braid_Core  core,
               braid_Int   level,
               braid_Real  rnorm)
{
   braid_App          app          = _braid_CoreElt(core, app);
   _braid_Grid      **grids        = _braid_CoreElt(core, grids);
   braid_AccessStatus astatus      = _braid_CoreElt(core, astatus);
   braid_Int          iter        = _braid_CoreElt(core, niter);
   braid_Int          access_level = _braid_CoreElt(core, access_level);
   braid_Int          nrefine      = _braid_CoreElt(core, nrefine);
   braid_Int          ilower       = _braid_GridElt(grids[level], ilower);
   braid_Int          ncpoints     = _braid_GridElt(grids[level], ncpoints);
   braid_Vector      *va           = _braid_GridElt(grids[level], va);
   braid_Real        *ta           = _braid_GridElt(grids[level], ta);

   braid_Int          f_level, f_cfactor, f_index;
   braid_Vector       f_u, f_e;

   braid_Vector       u, e;
   braid_Int          flo, fhi, fi, ci;
   braid_Int          interval;

   f_level   = level-1;
   f_cfactor = _braid_GridElt(grids[f_level], cfactor);

   _braid_UCommInitF(core, level);

   /**
    * Start from the right-most interval 
    *
    * First, generate the coarse-grid F-points through F-relaxation and
    * interpolate them to the fine grid, where they are C-points.  Second,
    * interpolate the coarse-grid C-points to the fine-grid.  The user-defined
    * spatial refinement (if set) is also called.  
    **/
   for (interval = ncpoints; interval > -1; interval--)
   {
      _braid_GetInterval(core, level, interval, &flo, &fhi, &ci);

      /* Relax and interpolate F-points, refining in space if needed */
      if (flo <= fhi)
      {
         _braid_UGetVector(core, level, flo-1, &u);
      }
      for (fi = flo; fi <= fhi; fi++)
      {
         _braid_Step(core, level, fi, NULL, u);
         _braid_USetVector(core, level, fi, u, 0);
         /* Allow user to process current vector 
          * We consider this iter+1, because we are on an up-cycle.*/
         if( (access_level >= 2) )
         {
            _braid_AccessStatusInit(ta[fi-ilower], rnorm, iter+1,
                                    level, nrefine, 0, 0, astatus);
            _braid_AccessVector(core, astatus, u);
         }
         e = va[fi-ilower];
         _braid_CoreFcn(core, sum)(app, 1.0, u, -1.0, e);
         _braid_MapCoarseToFine(fi, f_cfactor, f_index);
         _braid_Refine(core, f_level, f_index, fi, e, &f_e);
         _braid_UGetVectorRef(core, f_level, f_index, &f_u);
         _braid_CoreFcn(core, sum)(app, 1.0, f_e, 1.0, f_u);
         _braid_USetVectorRef(core, f_level, f_index, f_u);
         _braid_CoreFcn(core, free)(app, f_e);
      }
      if (flo <= fhi)
      {
         _braid_CoreFcn(core, free)(app, u);
      }

      /* Interpolate C-points, refining in space if needed */
      if (ci > 0)
      {
         _braid_UGetVectorRef(core, level, ci, &u);
         /* Allow user to process current C-point
          * We consider this iter+1, because we are on an up-cycle.*/
         if( (access_level >= 2) )
         {
            _braid_AccessStatusInit(ta[ci-ilower], rnorm, iter+1,
                                    level, nrefine, 0, 0, astatus);
            _braid_AccessVector(core, astatus, u);
         }
         e = va[ci-ilower];
         _braid_CoreFcn(core, sum)(app, 1.0, u, -1.0, e);
         _braid_MapCoarseToFine(ci, f_cfactor, f_index);
         _braid_Refine(core, f_level, f_index, ci, e, &f_e);
         _braid_UGetVectorRef(core, f_level, f_index, &f_u);
         _braid_CoreFcn(core, sum)(app, 1.0, f_e, 1.0, f_u);
         _braid_USetVectorRef(core, f_level, f_index, f_u);
         _braid_CoreFcn(core, free)(app, f_e);
      }
   }

   _braid_UCommWait(core, level);

   /* Clean up */
   _braid_GridClean(core, grids[level]);

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
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
 *   r_ta - time values on the refined grid                    (size 'r_npoints')
 *   r_fa - index map from coarse to fine on the refined grid  (size 'npoints+1')
 *          (note the extra value)
 *   f_ca - index map from fine to coarse on the fine grid     (size 'f_npoints')
 *
 *   send_ua - array of u-vectors to send to new processors    (size 'npoints')
 *   recv_ua - array of u-vectors received from old processors (size 'f_npoints')
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
 *--------------------------------------------------------------------------*/

braid_Int
_braid_FRefine(braid_Core   core,
               braid_Int    iter,
               braid_Real   rnorm,
               braid_Int   *refined_ptr)
{
   MPI_Comm       comm         = _braid_CoreElt(core, comm);
   braid_App      app          = _braid_CoreElt(core, app);
   braid_Int     *rfactors     = _braid_CoreElt(core, rfactors);
   _braid_Grid  **grids        = _braid_CoreElt(core, grids);
   braid_AccessStatus astatus  = _braid_CoreElt(core, astatus);
   braid_Int      access_level = _braid_CoreElt(core, access_level);
   braid_Int      nrefine      = _braid_CoreElt(core, nrefine);
   braid_Int      ncpoints     = _braid_GridElt(grids[0], ncpoints);

   /* Use prefix 'r_' for the refined version of the current grid level 0, and
    * use prefix 'f_' for the fine grid in the new distribution. */

   braid_Int      npoints, ilower, iupper, gupper, i, j, ii;
   braid_Int      r_npoints, r_ilower, r_iupper, r_i, r_ii;
   braid_Int      f_npoints, f_ilower, f_iupper, f_gupper, f_i, f_j, f_ii;
   braid_Int     *r_ca, *r_fa, *f_ca, f_first, f_next, next;
   braid_Real    *ta, *r_ta, *f_ta;

   braid_Vector  *send_ua, *recv_ua, u;
   braid_Int     *send_procs, *recv_procs, *send_unums, *recv_unums, *iptr;
   braid_Real    *send_buffer, *recv_buffer, **send_buffers, **recv_buffers, *bptr;
   void          *buffer;
   braid_Int      send_size, recv_size, *send_sizes, size, isize, max_usize;
   braid_Int      ncomms, nsends, nrecvs, nreceived, nprocs, proc, prevproc;
   braid_Int      unum, send_msg, recv_msg;
   MPI_Request   *requests, request;
   MPI_Status    *statuses, status;

   _braid_Grid   *f_grid;
   braid_Int      cfactor, rfactor, m, interval, flo, fhi, fi, ci, f_hi, f_ci;

#if DEBUG
   braid_Int  myproc;
   /*cfactor = 6;*/ /* RDF HACKED TEST */
   MPI_Comm_rank(comm, &myproc);
#endif

   gupper  = _braid_CoreElt(core, gupper);
   ilower  = _braid_GridElt(grids[0], ilower);
   iupper  = _braid_GridElt(grids[0], iupper);
   npoints = iupper - ilower + 1;

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
      rfactors[ii] = _braid_min(rfactors[ii], cfactor);
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
      *refined_ptr = 0;
      return _braid_error_flag;
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
   r_ta = _braid_CTAlloc(braid_Real, r_npoints);
   r_fa = _braid_CTAlloc(braid_Int,  npoints+1);
   ta = _braid_GridElt(grids[0], ta);
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
      /* Post r_fa receive */
      r_fa[npoints] = f_gupper+1;
      if (iupper < gupper)
      {
         MPI_Irecv(&r_fa[npoints], 1, braid_MPI_INT, MPI_ANY_SOURCE, 2, comm,
                   &requests[ncomms++]);
      }

      /* Post r_fa send */
      if (ilower > 0)
      {
         _braid_GetBlockDistProc((gupper+1), nprocs, (ilower-1), &prevproc);
         MPI_Isend(&r_fa[0], 1, braid_MPI_INT, prevproc, 2, comm,
                   &requests[ncomms++]);
      }
   }
   MPI_Waitall(ncomms, requests, statuses);
   _braid_TFree(requests);
   _braid_TFree(statuses);

   /* If storing only C-points on the fine grid, modify r_ca to mark only those
    * coarse points that need to be sent to initialize the C-points */
   if (_braid_CoreElt(core, storage) != 0)
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

   send_ua = _braid_CTAlloc(braid_Vector, npoints);
   send_procs = _braid_CTAlloc(braid_Int, npoints);
   send_unums = _braid_CTAlloc(braid_Int, npoints);
   send_buffers = _braid_CTAlloc(braid_Real *, npoints);

   recv_ua = _braid_CTAlloc(braid_Vector, f_npoints);
   recv_procs = _braid_CTAlloc(braid_Int, f_npoints);
   recv_unums = _braid_CTAlloc(braid_Int, f_npoints);
   recv_buffers = _braid_CTAlloc(braid_Real *, f_npoints);

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
               _braid_RefineBasic(core, -1, &r_ta[r_ii], &ta[ii], u, &send_ua[ii]);
            }

            /* Allow user to process current vector */
            if( (access_level >= 2) )
            {
               _braid_AccessStatusInit(ta[ii], rnorm, iter, 0, nrefine,
                                       0, 0, astatus);
               _braid_AccessVector(core, astatus, u);
            }
         }
         _braid_CoreFcn(core, free)(app, u);
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
            _braid_RefineBasic(core, -1, &r_ta[r_ii], &ta[ii], u, &send_ua[ii]);
         }

         /* Allow user to process current vector */
         if( (access_level >= 2) )
         {
            _braid_AccessStatusInit(ta[ii], rnorm, iter, 0, nrefine,
                                    0, 0, astatus);
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

   _braid_CoreFcn(core, bufsize)(app, &max_usize); /* max buffer size */
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
            _braid_CoreFcn(core, bufpack)(app, send_ua[ii], buffer, &size);
            _braid_CoreFcn(core, free)(app, send_ua[ii]);
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
            _braid_CoreFcn(core, bufunpack)(app, buffer, &recv_ua[f_ii]);
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
   _braid_TFree(r_ta);
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
   _braid_CoreElt(core, ntime)   = f_gupper;
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
         f_ci = _braid_NextCPoint(f_i, cfactor);
         if (next > f_ci)
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

/*--------------------------------------------------------------------------
 * Access to XBraid on grid level
 *--------------------------------------------------------------------------*/

braid_Int
_braid_FAccess(braid_Core     core,
               braid_Real     rnorm,
               braid_Int      level,
               braid_Int      done)
{
   braid_App           app      = _braid_CoreElt(core, app);
   _braid_Grid       **grids    = _braid_CoreElt(core, grids);
   braid_AccessStatus  astatus  = _braid_CoreElt(core, astatus);
   braid_Int           iter     = _braid_CoreElt(core, niter);
   braid_Int           nrefine  = _braid_CoreElt(core, nrefine);
   braid_Int           ncpoints = _braid_GridElt(grids[level], ncpoints);
   braid_Real          *ta      = _braid_GridElt(grids[level], ta);
   braid_Int           ilower   = _braid_GridElt(grids[level], ilower);

   braid_Vector   u;
   braid_Int      interval, flo, fhi, fi, ci;

   _braid_UCommInitF(core, level);

   /* Start from the right-most interval */
   for (interval = ncpoints; interval > -1; interval--)
   {
      _braid_GetInterval(core, level, interval, &flo, &fhi, &ci);

      /* Give access at F-points */
      if (flo <= fhi)
      {
         _braid_UGetVector(core, level, flo-1, &u);
      }
      for (fi = flo; fi <= fhi; fi++)
      {
         _braid_Step(core, level, fi, NULL, u);
         _braid_USetVector(core, level, fi, u, 0);
         _braid_AccessStatusInit( ta[fi-ilower], rnorm, iter,
                                  level, nrefine, done, 0, astatus);
         _braid_AccessVector(core, astatus, u);
      }
      if (flo <= fhi)
      {
         _braid_CoreFcn(core, free)(app, u);
      }

      /* Give access at C-points */
      if (ci > -1)
      {
         _braid_UGetVectorRef(core, level, ci, &u);
         _braid_AccessStatusInit( ta[ci-ilower], rnorm, iter,
                                  level, nrefine, done, 0, astatus);
         _braid_AccessVector(core, astatus, u);
      }
   }
   _braid_UCommWait(core, level);

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Initialize grid hierarchy
 *--------------------------------------------------------------------------*/

braid_Int
_braid_InitHierarchy(braid_Core    core,
                     _braid_Grid  *fine_grid,
                     braid_Int     refined)
{
   MPI_Comm       comm       = _braid_CoreElt(core, comm);
   braid_Int      max_levels = _braid_CoreElt(core, max_levels);
   braid_Int      min_coarse = _braid_CoreElt(core, min_coarse);
   braid_Real     tol        = _braid_CoreElt(core, tol);
   braid_Int     *nrels      = _braid_CoreElt(core, nrels);
   braid_Int      nrdefault  = _braid_CoreElt(core, nrdefault);
   braid_Int      gupper     = _braid_CoreElt(core, gupper);
   braid_Int     *rfactors   = _braid_CoreElt(core, rfactors);
   braid_Int      nlevels    = _braid_CoreElt(core, nlevels);
   _braid_Grid  **grids      = _braid_CoreElt(core, grids);

   /**
    * These are some common index names used to refer to intervals and
    * time points.  Here's what they mean.
    *
    * cfactor            - coarsening factor, fixed on each level
    * ilower, iupper     - lowest and highest time indices on a level for one processor, 
    *                      could be C or F points
    * clower, cupper     - lowest and highest C-point indices on a level for one processor
    *                      analagous to ilower, iupper being projected onto C-points
    * clo, chi           - the coarse level indices for clower and cupper
    * f_iupper, f_ilower - lowest and highest F-point indices on a level
    * c_iupper, c_ilower - lowest and highest time indices on the coarse level 
    *                      (analagous to ilower and iupper, but on the next level down)
    * flo, fhi, ci       - describes an interval [ci, flo, flo+1, ..., fhi] 
    * gclower, gcupper   - global lowest and highest C-point indices on a level
    **/

   braid_Int      level;
   braid_Int      ilower, iupper;
   braid_Int      clower, cupper, cfactor, ncpoints, nupoints;
   braid_Real    *ta;
   braid_Vector  *ua;
   braid_Vector  *va;
   braid_Vector  *fa;

   _braid_Grid   *grid;
   braid_Real    *f_ta;
   braid_Int      i, f_i, f_ilower, clo, chi, gclower, gcupper;

   MPI_Request    request1, request2;
   MPI_Status     status;
   braid_Int      left_proc, right_proc;

   grids[0] = fine_grid;

   /* Do sequential time marching if tolerance is not positive, 
    * or min_coarse is already reached */

   if ((tol <= 0.0) && (max_levels > 1) && (gupper >= min_coarse) )
   {
      max_levels = 1;
      _braid_CoreElt(core, max_levels) = max_levels;
   }

   /* Allocate space for rfactors (and initialize to zero) */
   ilower = _braid_GridElt(grids[0], ilower);
   iupper = _braid_GridElt(grids[0], iupper);
   rfactors = _braid_CTAlloc(braid_Int, iupper-ilower+2); /* Ensures non-NULL */
   rfactors[0] = 1; /* Ensures rfactor of 1 for global index 0 */
   _braid_CoreElt(core, rfactors) = rfactors;

   /* Set up nrels array */
   for (level = 0; level < max_levels; level++)
   {
      if (nrels[level] < 0)
      {
         nrels[level] = nrdefault;
      }
   }

   /* Coarsen global grid to determine nlevels */
   gclower = 0;
   gcupper = gupper;
   for (level = 0; level < max_levels; level++)
   {
      grid = grids[level];
      ilower = _braid_GridElt(grid, ilower);
      iupper = _braid_GridElt(grid, iupper);
      if (level > 0)
      {
         /* Copy ta info from level-1 grid */
         ta       = _braid_GridElt(grid, ta);
         f_ilower = _braid_GridElt(grids[level-1], ilower);
         f_ta     = _braid_GridElt(grids[level-1], ta);
         cfactor  = _braid_GridElt(grids[level-1], cfactor);
         for (i = ilower; i <= iupper; i++)
         {
            _braid_MapCoarseToFine(i, cfactor, f_i);
            ta[i-ilower] = f_ta[f_i-f_ilower];
         }
      }

      _braid_GetCFactor(core, level, &cfactor);
      
      _braid_GridElt(grid, gupper) = gcupper;

      _braid_ProjectInterval(gclower, gcupper, 0, cfactor, &gclower, &gcupper);
      _braid_MapFineToCoarse(gclower, cfactor, gclower);
      _braid_MapFineToCoarse(gcupper, cfactor, gcupper);

      /* Coarsen */
      _braid_ProjectInterval(ilower, iupper, 0, cfactor, &clower, &cupper);
      _braid_MapFineToCoarse(clower, cfactor, clo);
      _braid_MapFineToCoarse(cupper, cfactor, chi);
      ncpoints = chi-clo+1;
      if (ncpoints < 0)
      {
         ncpoints = 0;
      }
      _braid_GridElt(grid, clower)   = clower;
      _braid_GridElt(grid, cupper)   = cupper;
      _braid_GridElt(grid, cfactor)  = cfactor;
      _braid_GridElt(grid, ncpoints) = ncpoints;
      if ( (gclower < gcupper) && (max_levels > level+1) &&
           ((gcupper - gclower) >= min_coarse) )
      {
         /* Initialize the coarse grid */
         _braid_GridInit(core, level+1, clo, chi, &grids[level+1]);
      }
      else
      {
         /* This is the coarsest level */
         if ( (level > 0) || (!refined) )
         {
            /* If this is a true coarse level (it has a fine grid above it in
             * the current hierarchy) or it is a fine level that was not built
             * by refining a coarser grid, then do serial time integration by
             * setting only one C-point and the rest F-points */
            if (ilower == 0)
            {
               ncpoints = 1;
            }
            else
            {
               ncpoints = 0;
            }
            /* clower > cupper indicates empty interval */
            _braid_GridElt(grid, clower)   = ilower;
            _braid_GridElt(grid, cupper)   = 0;
            _braid_GridElt(grid, cfactor)  = gupper+1;
            _braid_GridElt(grid, ncpoints) = ncpoints;
         }

         /* Stop coarsening */
         break;
      }
      
      if(level == 0)
      {   
         /* Allocate space for storage of residual norm at each C-point */
         _braid_CoreElt(core, tnorm_a)  = _braid_CTAlloc(braid_Real, ncpoints);
      }
   }
   nlevels = level+1;
   _braid_CoreElt(core, nlevels) = nlevels;

   /* Allocate ua, va, and fa here */
   for (level = 0; level < nlevels; level++)
   {
      grid = grids[level];
      ilower = _braid_GridElt(grid, ilower);
      iupper = _braid_GridElt(grid, iupper);
      if (level > 0)
      {
         va = _braid_CTAlloc(braid_Vector, iupper-ilower+2);
         fa = _braid_CTAlloc(braid_Vector, iupper-ilower+2);
         _braid_GridElt(grid, va_alloc) = va;
         _braid_GridElt(grid, fa_alloc) = fa;
         _braid_GridElt(grid, va)       = va+1;  /* shift */
         _braid_GridElt(grid, fa)       = fa+1;  /* shift */
      }

      if ( (_braid_CoreElt(core, storage) == -1) ||
           (level < _braid_CoreElt(core, storage)) )
      {
         nupoints = _braid_GridElt(grid, ncpoints);   /* only C-points */
      }
      else
      {
         nupoints = iupper-ilower+1;                  /* all points */
      }
      ua = _braid_CTAlloc(braid_Vector, nupoints+1);
      _braid_GridElt(grid, nupoints)  = nupoints;
      _braid_GridElt(grid, ua_alloc)  = ua;
      _braid_GridElt(grid, ua)        = ua+1;  /* shift */
   }

   /* Communicate ta[-1] and ta[iupper-ilower+1] information */
   for (level = 0; level < nlevels; level++)
   {
      grid = grids[level];
      ilower = _braid_GridElt(grid, ilower);
      iupper = _braid_GridElt(grid, iupper);
      ta     = _braid_GridElt(grid, ta);

      if (ilower <= iupper)
      {
         _braid_GetProc(core, level, ilower-1, &left_proc);
         _braid_GetProc(core, level, iupper+1, &right_proc);
         
         /* Post receive to set ta[-1] on each processor*/
         if (left_proc > -1)
         {
            MPI_Irecv(&ta[-1], sizeof(braid_Real), MPI_BYTE,
                      left_proc, 1, comm, &request1);
         }
         else
         {
             /* Place a repeat value to indicate the start of the time-line for this level */
             ta[-1] = ta[0]; 
         }
         /* Post receive to set ta[iupper-ilower+1] on each processor */
         if ( _braid_CoreElt(core, coarsen) != NULL )
         {
             if (right_proc > -1)
             {
                MPI_Irecv(&ta[iupper-ilower+1], sizeof(braid_Real), MPI_BYTE,
                          right_proc, 1, comm, &request2);
             }
             else
             {
                 /* Place a repeat value to indicate the end the time-line for this level */
                 ta[iupper-ilower+1] = ta[iupper-ilower];
             }
         }

         /* Post send that sets ta[-1] on each processor */
         if (right_proc > -1)
         {
            MPI_Send(&ta[iupper-ilower], sizeof(braid_Real), MPI_BYTE,
                     right_proc, 1, comm);
         }
         /* Post send that sets ta[iupper-ilower+1] on each processor */
         if ( (left_proc > -1) && ( _braid_CoreElt(core, coarsen) != NULL ) )
         {
            MPI_Send(&ta[0], sizeof(braid_Real), MPI_BYTE, left_proc, 1, comm);
         }

         /* Finish receive */
         if (left_proc > -1)
         {
            MPI_Wait(&request1, &status);
         }
         if ( (right_proc > -1) && ( _braid_CoreElt(core, coarsen) != NULL ) )
         {
            MPI_Wait(&request2, &status);
         }
      }
   }

   return _braid_error_flag;
}
