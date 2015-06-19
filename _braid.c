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

#ifndef DEBUG
#define DEBUG 0
#endif

braid_Int _braid_error_flag = 0;
FILE    *_braid_printfile  = NULL;

/*--------------------------------------------------------------------------
 * Determine processor distribution.  This must agree with GetProc().
 *--------------------------------------------------------------------------*/

braid_Int
_braid_GetDistribution(braid_Core   core,
                       braid_Int   *ilower_ptr,
                       braid_Int   *iupper_ptr)
{
   MPI_Comm   comm    = _braid_CoreElt(core, comm);
   braid_Int  gupper = _braid_CoreElt(core, gupper);
   braid_Int  ilower, iupper;
   braid_Int  npoints, nprocs, proc, quo, rem, p;

   MPI_Comm_size(comm, &nprocs);
   MPI_Comm_rank(comm, &proc);

   npoints = gupper + 1;
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
 * Returns -1 if index is out of range
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
   braid_Int      proc;
   braid_Int      npoints, nprocs, quo, rem, p, q;
   braid_Int      l, cfactor;

   /* MPI_Comm_rank(comm, &proc); */
   MPI_Comm_size(comm, &nprocs);

   npoints = gupper + 1;
   quo = npoints/nprocs;
   rem = npoints%nprocs;

   /* Map index to the finest grid */
   for (l = level-1; l > -1; l--)
   {
      cfactor = _braid_GridElt(grids[l], cfactor);
      _braid_MapCoarseToFine(index, cfactor, index);
   }

   /* Compute processor number */
   if ((index < 0) || (index > gupper))
   {
      proc = -1;
   }
   else
   {
      p = index/(quo+1);
      q = (index - rem*(quo+1))/quo;
      proc = (p < rem ? p : rem+q);
   }

   *proc_ptr = proc;

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
         uindex = index - ilower;
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

   /* No user-provided residual routine */
   if ( _braid_CoreElt(core, residual) == NULL )
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

   /* User-provided residual routine, store all u-vectors */
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
            braid_Int      step_type,
            braid_Int      index,
            braid_Real     accuracy,
            braid_Vector   ustop,
            braid_Vector   u)
{
   braid_App        app      = _braid_CoreElt(core, app);
   braid_Int       *rfactors = _braid_CoreElt(core, rfactors);
   _braid_Grid    **grids    = _braid_CoreElt(core, grids);
   braid_StepStatus status   = _braid_CoreElt(core, sstatus);
   braid_Int        ilower   = _braid_GridElt(grids[level], ilower);
   braid_Real      *ta       = _braid_GridElt(grids[level], ta);
   braid_Vector    *fa       = _braid_GridElt(grids[level], fa);

   braid_Int        ii;

   ii = index-ilower;
   _braid_StepStatusInit(ta[ii-1], ta[ii], accuracy, level, step_type, status);

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
                braid_Real     accuracy,
                braid_Vector   ustop,
                braid_Vector   r)
{
   braid_App        app      = _braid_CoreElt(core, app);
   _braid_Grid    **grids    = _braid_CoreElt(core, grids);
   braid_StepStatus status   = _braid_CoreElt(core, sstatus);
   braid_Int        ilower   = _braid_GridElt(grids[level], ilower);
   braid_Real      *ta       = _braid_GridElt(grids[level], ta);

   braid_Vector     rstop;
   braid_Int        ii;
   braid_Int        step_type = 2;        /* A type of 2 denotes a residual computation */

   ii = index-ilower;
   _braid_StepStatusInit(ta[ii-1], ta[ii], accuracy, level, step_type, status);
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
                   braid_Real     accuracy,
                   braid_Vector   ustop,
                   braid_Vector   r)
{
   braid_App        app    = _braid_CoreElt(core, app);
   _braid_Grid    **grids  = _braid_CoreElt(core, grids);
   braid_Int        ilower = _braid_GridElt(grids[level], ilower);
   braid_Vector    *fa     = _braid_GridElt(grids[level], fa);

   braid_Int        ii;

   _braid_Residual(core, level, index, accuracy, ustop, r);
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
                                 c_ta[c_ii-1], c_ta[c_ii+1], level-1, cstatus);
      _braid_CoreFcn(core, coarsen)(app, fvector, cvector, cstatus);
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
   braid_App      app             = _braid_CoreElt(core, app);
   _braid_Grid  **grids           = _braid_CoreElt(core, grids);
   braid_CoarsenRefStatus cstatus = _braid_CoreElt(core, cstatus);
   braid_Int      f_ilower        = _braid_GridElt(grids[level], ilower);
   braid_Int      c_ilower        = _braid_GridElt(grids[level+1], ilower);
   braid_Real    *f_ta            = _braid_GridElt(grids[level], ta);
   braid_Real    *c_ta            = _braid_GridElt(grids[level+1], ta);

   braid_Int      c_ii = c_index-c_ilower;
   braid_Int      f_ii = f_index-f_ilower;

   if ( _braid_CoreElt(core, coarsen) == NULL )
   {
      /* No spatial refinement needed, just clone the fine vector.*/
      _braid_CoreFcn(core, clone)(app, cvector, fvector);
   }
   else
   {
      /* Call the user's refinement routine */
      _braid_CoarsenRefStatusInit(f_ta[f_ii], f_ta[f_ii-1], f_ta[f_ii+1], 
                                 c_ta[c_ii-1], c_ta[c_ii+1], level, cstatus);
      _braid_CoreFcn(core, refine)(app, cvector, fvector, cstatus);
   }

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
   _braid_Grid      **grids       = _braid_CoreElt(core, grids);
   braid_Int          ncpoints    = _braid_GridElt(grids[level], ncpoints);
   braid_Int          tnorm       = _braid_CoreElt(core, tnorm);
   braid_Real        *ta          = _braid_GridElt(grids[level], ta);
   braid_Int          ilower      = _braid_GridElt(grids[level], ilower);
   braid_StepStatus status   = _braid_CoreElt(core, sstatus);
   _braid_CommHandle *send_handle;
   braid_Int          send_index;

   braid_Int        flo, fhi, fi, ci, ii, interval;
   braid_Real       accuracy, rnorm_temp, rnorm = 0, global_rnorm = 0;
   braid_Vector     u, r;
   braid_Int        step_type = 2;        /* A type of 2 denotes a residual computation */


   if ( level == 0 )
   {     
      /*accuracy = _braid_CoreElt(core, accuracy[0].value);*/
      accuracy = _braid_CoreElt(core, accuracy[0].old_value);
   }
   else
   {
      accuracy = _braid_CoreElt(core, accuracy[1].value);
   }

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
         _braid_Step(core, level, step_type, fi, accuracy, NULL, u);

         /* Update local processor norm. */
         ii = fi-ilower;
         _braid_StepStatusInit(ta[ii-1], ta[ii], accuracy, level, step_type, status);
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
         _braid_StepStatusInit(ta[ii-1], ta[ii], accuracy, level, step_type, status);
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
   braid_Real      accuracy;
   braid_Int       step_type = 0;        /* A type of 0 denotes relaxation  */

   if ( level == 0 )
   {     
      /*accuracy = _braid_CoreElt(core, accuracy[0].value);*/
      accuracy = _braid_CoreElt(core, accuracy[0].old_value);
   }
   else
   {
      accuracy = _braid_CoreElt(core, accuracy[1].value);
   }

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
            _braid_Step(core, level, step_type, fi, accuracy, NULL, u);
            _braid_USetVector(core, level, fi, u, 0);
         }

         /* C-relaxation */
         if (ci > 0)
         {
            _braid_Step(core, level, step_type, ci, accuracy, NULL, u);
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
                 braid_Int    iter,
                 braid_Real  *rnorm_ptr)
{
   MPI_Comm             comm_world  = _braid_CoreElt(core, comm_world);
   MPI_Comm             comm        = _braid_CoreElt(core, comm);
   braid_App            app         = _braid_CoreElt(core, app);
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_AccessStatus   astatus     = _braid_CoreElt(core, astatus);
   braid_Int            print_level = _braid_CoreElt(core, print_level);
   braid_Int            access_level= _braid_CoreElt(core, access_level);
   braid_Int            tnorm       = _braid_CoreElt(core, tnorm);
   braid_Real          *tnorm_a     = _braid_CoreElt(core, tnorm_a);
   braid_Int            cfactor     = _braid_GridElt(grids[level], cfactor);
   braid_Int            ncpoints    = _braid_GridElt(grids[level], ncpoints);
   braid_Real          *ta          = _braid_GridElt(grids[level], ta);
   braid_Int            f_ilower    = _braid_GridElt(grids[level], ilower);
   _braid_CommHandle   *recv_handle = NULL;
   _braid_CommHandle   *send_handle = NULL;
   braid_Real           old_rnorm   = *rnorm_ptr;

   braid_Int            myid;
   braid_Int            c_level, c_ilower, c_iupper, c_index, c_i, c_ii;
   braid_Vector         c_u, *c_va, *c_fa;

   braid_Vector         u, r;
   braid_Int            interval, flo, fhi, fi, ci;
   braid_Real           rnorm, grnorm, rnorm_temp, accuracy, c_accuracy;
   braid_Int            step_type = 2;    /* A type of 2 denotes a residual computation */

   c_level  = level+1;
   c_ilower = _braid_GridElt(grids[c_level], ilower);
   c_iupper = _braid_GridElt(grids[c_level], iupper);
   c_va     = _braid_GridElt(grids[c_level], va);
   c_fa     = _braid_GridElt(grids[c_level], fa);

   rnorm = 0.0;
   MPI_Comm_rank(comm_world, &myid);

   if (level == 0)
   {
      accuracy = _braid_CoreElt(core, accuracy[0].value);
      if (accuracy == _braid_CoreElt(core, accuracy[0].tight))
      {
         _braid_CoreElt(core, accuracy[0].tight_used) = 1;
      }
   }
   else
   {
      accuracy = _braid_CoreElt(core, accuracy[1].value);
   }
   c_accuracy = _braid_CoreElt(core, accuracy[1].value);

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
         _braid_Step(core, level, step_type, fi, accuracy, NULL, r);
         _braid_USetVector(core, level, fi, r, 0);
         
         /* Allow user to process current vector, note that r here is
          * temporarily holding the state vector */
         if( (access_level >= 2) && (level == 0) )
         {
            _braid_AccessStatusInit(ta[fi-f_ilower], old_rnorm, iter, level, 0, 0, astatus);
            _braid_AccessVector(core, astatus, r);
         }
      }

      /* Allow user to process current C-point */
      if( (access_level>= 2) && (level == 0) && (ci > -1) )
      {
         _braid_UGetVectorRef(core, level, ci, &u);
         _braid_AccessStatusInit(ta[ci-f_ilower], old_rnorm, iter, level, 0, 0, astatus);
         _braid_AccessVector(core, astatus, u);
      }
      
      /* Compute residual and restrict */
      if (ci > 0)
      {
         /* Compute FAS residual */
         _braid_UGetVectorRef(core, level, ci, &u);
         _braid_FASResidual(core, level, ci, accuracy, u, r);

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
         _braid_Residual(core, c_level, c_i, c_accuracy, c_va[c_ii], c_u);
         _braid_CoreFcn(core, sum)(app, 1.0, c_u, 1.0, c_fa[c_ii]);
         _braid_CoreFcn(core, free)(app, c_u);
      }
   }
   _braid_CommWait(core, &send_handle);

   /* Set initial guess on coarse level */
   _braid_InitGuess(core, c_level);

   /* Compute rnorm (only on level 0) */
   if (level == 0)
   {
      if(tnorm == 1)          /* one-norm reduction */
      {  
         MPI_Allreduce(&rnorm, &grnorm, 1, MPI_DOUBLE, MPI_SUM, comm);
      }
      else if(tnorm == 3)     /* inf-norm reduction */
      {  
         _braid_Max(tnorm_a, ncpoints, &rnorm); 
         MPI_Allreduce(&rnorm, &grnorm, 1, MPI_DOUBLE, MPI_MAX, comm);
      }
      else                    /* default two-norm reduction */
      {  
         MPI_Allreduce(&rnorm, &grnorm, 1, MPI_DOUBLE, MPI_SUM, comm);
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
               braid_Int   iter,
               braid_Real  rnorm)
{
   braid_App          app         = _braid_CoreElt(core, app);
   _braid_Grid      **grids       = _braid_CoreElt(core, grids);
   braid_AccessStatus astatus     = _braid_CoreElt(core, astatus);
   braid_Int          access_level= _braid_CoreElt(core, access_level);
   braid_Int          ilower      = _braid_GridElt(grids[level], ilower);
   braid_Int          ncpoints    = _braid_GridElt(grids[level], ncpoints);
   braid_Vector      *va          = _braid_GridElt(grids[level], va);
   braid_Real        *ta          = _braid_GridElt(grids[level], ta);

   braid_Int          f_level, f_cfactor, f_index;
   braid_Vector       f_u, f_e;

   braid_Vector       u, e;
   braid_Int          flo, fhi, fi, ci;
   braid_Int          interval;
   braid_Int          step_type = 1;   /* A type of 1 denotes interpolation */
   braid_Real         accuracy;

   f_level   = level-1;
   f_cfactor = _braid_GridElt(grids[f_level], cfactor);

   accuracy = _braid_CoreElt(core, accuracy[1].value);

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
         _braid_Step(core, level, step_type, fi, accuracy, NULL, u);
         _braid_USetVector(core, level, fi, u, 0);
         /* Allow user to process current vector */
         if( (access_level >= 2) )
         {
            _braid_AccessStatusInit(ta[fi-ilower], rnorm, iter, level, 0, 0, astatus);
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
         /* Allow user to process current C-point */
         if( (access_level >= 2) )
         {
            _braid_AccessStatusInit(ta[ci-ilower], rnorm, iter, level, 0, 0, astatus);
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
 * Create a new fine grid based on user refinement factor information, then
 * F-relax and interpolate to the new fine grid and create a new multigrid
 * hierarchy.  In general, this will require load re-balancing as well.
 *
 * RDF: Todo
 *--------------------------------------------------------------------------*/

braid_Int
_braid_FRefine(braid_Core   core,
               braid_Int   *refined_ptr)
{
#if 0 /* Finish writing later */
   braid_Int     *rfactors = _braid_CoreElt(core, rfactors);
   braid_Int      nlevels  = _braid_CoreElt(core, nlevels);
   _braid_Grid  **grids    = _braid_CoreElt(core, grids);
   braid_Int      ilower   = _braid_GridElt(grids[0], ilower);
   braid_Int      iupper   = _braid_GridElt(grids[0], iupper);
   braid_Int      cfactor  = _braid_GridElt(grids[0], cfactor);

   braid_Int      lrefine, refine;
   braid_Int      step_type = 1;   /* A type of 1 denotes interpolation */

   lrefine = 0;
   for (i = ilower; i <= iupper; i++)
   {
      if (rfactors[i] > 1)
      {
         lrefine = 1;
         break;
      }
   }
   MPI_Allreduce(&lrefine, &refine, 1, MPI_INT, MPI_MAX, comm);

   if (refine)
   {
      /* Compute and save interpolated values */
      ua = _braid_CTAlloc(braid_Vector, (ilower-iupper+1));

      _braid_UCommInitF(core, 0);

      /* Start from the right-most interval */
      for (interval = ncpoints; interval > -1; interval--)
      {
         _braid_GetInterval(core, 0, interval, &flo, &fhi, &ci);

         /* Relax and interpolate F-points */
         if (flo <= fhi)
         {
            _braid_UGetVector(core, 0, flo-1, &u);
         }
         for (fi = flo; fi <= fhi; fi++)
         {
            _braid_Step(core, 0, step_type, fi, 1.0, NULL, u);
            _braid_USetVector(core, 0, fi, u, 0);
            _braid_CoreFcn(core, clone)(app, u, &ua[fi-ilower]);
         }
         if (flo <= fhi)
         {
            _braid_CoreFcn(core, free)(app, u);
         }

         /* Interpolate C-points */
         if (ci > 0)
         {
            _braid_UGetVectorRef(core, 0, ci, &u);
            ua[ci-ilower] = u;
         }
      }

      _braid_UCommWait(core, 0);

      /* Create new finest grid */
      f_cfactor = cfactor;
      _braid_MapCoarseToFine(ilower, f_cfactor, f_ilower);
      _braid_MapCoarseToFine(iupper, f_cfactor, f_iupper);
      f_iupper += f_cfactor-1;
      _braid_GridInit(core, 0, f_ilower, f_iupper, &f_grid);

      /* Change gupper (this currently does not match tstop) */
      gupper = _braid_CoreElt(core, gupper);
      gupper = (gupper+1)*f_cfactor - 1;
      _braid_CoreElt(core, gupper) = gupper;

      /* Set t values */
      f_ta = _braid_GridElt(f_grid, ta);
      ta   = _braid_GridElt(grids[0], ta);
      for (i = ilower; i <= iupper; i++)
      {
         _braid_MapCoarseToFine(i, f_cfactor, f_i);
         tstart = ta[i-ilower];
         /* RDF - START HERE*/
         for (j = 0; j < f_cfactor; j++)
         {
            ta[i-ilower] = f_ta[f_i-f_ilower];
            f_ta[f_i-f_ilower] = tstart + (((braid_Real)i)/ntime)*(tstop-tstart);
         }
      }

      /* Initialize new hierarchy */
      for (level = 0; level < nlevels; level++)
      {
         _braid_GridDestroy(core, grids[level]);
      }
      _braid_InitHierarchy(core, f_grid);
      nlevels = _braid_CoreElt(core, nlevels);

      /* Do the spatial refinement here and set the fine grid unknowns */
      for (i = ilower; i <= iupper; i++)
      {
         _braid_MapCoarseToFine(i, f_cfactor, f_index);
         _braid_Refine(core, 0, f_index, i, ua[i-ilower], &f_u);
         _braid_USetVectorRef(core, 0, f_index, f_u);
         _braid_CoreFcn(core, free)(app, ua[i-ilower]);
      }

      _braid_TFree(ua);
   }

   *refined_ptr = refine;
#endif

   *refined_ptr = 0;

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Access to XBraid on grid level
 *--------------------------------------------------------------------------*/

braid_Int
_braid_FAccess(braid_Core     core,
               braid_Real     rnorm,
               braid_Int      iter,
               braid_Int      level,
               braid_Int      done)
{
   braid_App           app      = _braid_CoreElt(core, app);
   _braid_Grid       **grids    = _braid_CoreElt(core, grids);
   braid_AccessStatus  astatus  = _braid_CoreElt(core, astatus);
   braid_Int           ncpoints = _braid_GridElt(grids[level], ncpoints);
   braid_Real          *ta      = _braid_GridElt(grids[level], ta);
   braid_Int           ilower   = _braid_GridElt(grids[level], ilower);
   braid_Real          accuracy;

   braid_Vector   u;
   braid_Int      interval, flo, fhi, fi, ci;
   braid_Int      step_type = 3;   /* A type of 3 denotes an access */

   _braid_UCommInitF(core, level);

   if (level == 0)
   {    
      /*accuracy = _braid_CoreElt(core, accuracy[0].value); */
      accuracy = _braid_CoreElt(core, accuracy[0].old_value);
   }
   else
   {
      accuracy = _braid_CoreElt(core, accuracy[1].value);
   }

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
         _braid_Step(core, level, step_type, fi, accuracy, NULL, u);
         _braid_USetVector(core, level, fi, u, 0);
         _braid_AccessStatusInit( ta[fi-ilower], rnorm, iter, level, done, 0, astatus);
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
         _braid_AccessStatusInit( ta[ci-ilower], rnorm, iter, level, done, 0, astatus);
         _braid_AccessVector(core, astatus, u);
      }
   }
   _braid_UCommWait(core, level);

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 * Initialize (and re-initialize) hierarchy
 *--------------------------------------------------------------------------*/

braid_Int
_braid_InitHierarchy(braid_Core    core,
                     _braid_Grid  *fine_grid)
{
   MPI_Comm       comm       = _braid_CoreElt(core, comm);
   braid_Int      max_levels = _braid_CoreElt(core, max_levels);
   braid_Int      min_coarse = _braid_CoreElt(core, min_coarse);
   braid_Real     tol        = _braid_CoreElt(core, tol);
   braid_Int     *nrels      = _braid_CoreElt(core, nrels);
   braid_Int      nrdefault  = _braid_CoreElt(core, nrdefault);
   braid_Int     *cfactors   = _braid_CoreElt(core, cfactors);
   braid_Int      cfdefault  = _braid_CoreElt(core, cfdefault);
   braid_Int      gupper     = _braid_CoreElt(core, gupper);
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

   MPI_Request   request1, request2;
   MPI_Status    status;
   braid_Int      left_proc, right_proc;

   grids[0] = fine_grid;

   /* Do sequential time marching if tolerance is not positive, 
    * or min_coarse is already reached */

   if ((tol <= 0.0) && (max_levels > 1) && (gupper >= min_coarse) )
   {
      max_levels = 1;
      _braid_CoreElt(core, max_levels) = max_levels;
   }

   /* Allocate space for rfactors */
   ilower = _braid_GridElt(grids[0], ilower);
   iupper = _braid_GridElt(grids[0], iupper);
   _braid_CoreElt(core, rfactors) = _braid_CTAlloc(braid_Int, iupper-ilower+1);

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

      if (cfactors[level] != 0)
      {
         cfactor = cfactors[level];
      }
      else
      {
         cfactor = cfdefault;
      }
      
      _braid_GridElt(grid, gupper)   = gcupper;

      _braid_ProjectInterval(gclower, gcupper, 0, cfactor, &gclower, &gcupper);
      _braid_MapFineToCoarse(gclower, cfactor, gclower);
      _braid_MapFineToCoarse(gcupper, cfactor, gcupper);
      if ( (gclower < gcupper) && (max_levels > level+1) &&
           ((gcupper - gclower) >= min_coarse) )
      {
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

         _braid_GridInit(core, level+1, clo, chi, &grids[level+1]);
      }
      else
      {
         /* Make this the coarsest level with one C-point and the rest F-points.
          * This ensures serial time integration at this level. */
         if (ilower == 0)
         {
            ncpoints = 1;
         }
         else
         {
            ncpoints = 0;
         }
         _braid_GridElt(grid, clower)   = ilower;     /* clower > cupper indicates empty interval */
         _braid_GridElt(grid, cupper)   = 0;
         _braid_GridElt(grid, cfactor)  = gupper+1;
         _braid_GridElt(grid, ncpoints) = ncpoints;

         /* Stop coarsening */
         break;
      }
      
      if( level == 0)
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

   /* Communicate ta[-1] and ta[iupper+1] information */
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
                      left_proc, 0, comm, &request1);
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
                          right_proc, 0, comm, &request2);
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
                     right_proc, 0, comm);
         }
         /* Post send that sets ta[iupper-ilower+1] on each processor */
         if ( (left_proc > -1) && ( _braid_CoreElt(core, coarsen) != NULL ) )
         {
            MPI_Send(&ta[0], sizeof(braid_Real), MPI_BYTE, left_proc, 0, comm);
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


