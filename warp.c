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

warp_Int _warp_error_flag;

/*--------------------------------------------------------------------------
 * Determine processor distribution.  This must agree with GetProc().
 *--------------------------------------------------------------------------*/

warp_Int
_warp_GetDistribution(warp_Core   core,
                      warp_Int   *ilower_ptr,
                      warp_Int   *iupper_ptr)
{
   MPI_Comm  comm   = _warp_CoreElt(core, comm);
   warp_Int  gupper = _warp_CoreElt(core, gupper);

   warp_Int  ilower, iupper;
   warp_Int  npoints, nprocs, proc, quo, rem, p;

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

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Returns -1 if index is out of range
 *--------------------------------------------------------------------------*/

warp_Int
_warp_GetProc(warp_Core   core,
              warp_Int    level,
              warp_Int    index,
              warp_Int   *proc_ptr)
{
   MPI_Comm      comm   = _warp_CoreElt(core, comm);
   _warp_Grid  **grids  = _warp_CoreElt(core, grids);
   warp_Int      gupper = _warp_CoreElt(core, gupper);

   warp_Int      proc;
   warp_Int      npoints, nprocs, quo, rem, p, q;
   warp_Int      l, cfactor;

   /* MPI_Comm_rank(comm, &proc); */
   MPI_Comm_size(comm, &nprocs);

   npoints = gupper + 1;
   quo = npoints/nprocs;
   rem = npoints%nprocs;

   /* Map index to the finest grid */
   for (l = level-1; l > -1; l--)
   {
      cfactor = _warp_GridElt(grids[l], cfactor);
      _warp_MapCoarseToFine(index, cfactor, index);
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

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
_warp_CommRecvInit(warp_Core           core,
                   warp_Int            level,
                   warp_Int            index,
                   warp_Vector        *vector_ptr,
                   _warp_CommHandle  **handle_ptr)
{
   MPI_Comm          comm = _warp_CoreElt(core, comm);
   warp_App          app  = _warp_CoreElt(core, app);

   _warp_CommHandle  *handle = NULL;
   void              *buffer;
   MPI_Request       *requests;
   MPI_Status        *status;
   warp_Int           proc, size, num_requests;

   _warp_GetProc(core, level, index, &proc);
   if (proc > -1)
   {
      handle = _warp_TAlloc(_warp_CommHandle, 1);

      /* Allocate buffer through user routine */
      _warp_CoreFcn(core, bufsize)(app, &size);
      buffer = malloc(size);

      num_requests = 1;
      requests = _warp_CTAlloc(MPI_Request, num_requests);
      status   = _warp_CTAlloc(MPI_Status, num_requests);
      MPI_Irecv(buffer, size, MPI_BYTE, proc, 0, comm, &requests[0]);

      _warp_CommHandleElt(handle, request_type) = 1; /* recv type = 1 */
      _warp_CommHandleElt(handle, num_requests) = num_requests;
      _warp_CommHandleElt(handle, requests)     = requests;
      _warp_CommHandleElt(handle, status)       = status;
      _warp_CommHandleElt(handle, buffer)       = buffer;
      _warp_CommHandleElt(handle, vector_ptr)   = vector_ptr;
   }

   *handle_ptr = handle;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
_warp_CommSendInit(warp_Core           core,
                   warp_Int            level,
                   warp_Int            index,
                   warp_Vector         vector,
                   _warp_CommHandle  **handle_ptr)
{
   MPI_Comm           comm = _warp_CoreElt(core, comm);
   warp_App           app  = _warp_CoreElt(core, app);

   _warp_CommHandle  *handle = NULL;
   void              *buffer;
   MPI_Request       *requests;
   MPI_Status        *status;
   warp_Int           proc, size, num_requests;

   _warp_GetProc(core, level, index+1, &proc);
   if (proc > -1)
   {
      handle = _warp_TAlloc(_warp_CommHandle, 1);

      /* Allocate buffer through user routine */
      _warp_CoreFcn(core, bufsize)(app, &size);
      buffer = malloc(size);
      _warp_CoreFcn(core, bufpack)(app, vector, buffer);

      num_requests = 1;
      requests = _warp_CTAlloc(MPI_Request, num_requests);
      status   = _warp_CTAlloc(MPI_Status, num_requests);
      MPI_Isend(buffer, size, MPI_BYTE, proc, 0, comm, &requests[0]);

      _warp_CommHandleElt(handle, request_type) = 0; /* send type = 0 */
      _warp_CommHandleElt(handle, num_requests) = num_requests;
      _warp_CommHandleElt(handle, requests)     = requests;
      _warp_CommHandleElt(handle, status)       = status;
      _warp_CommHandleElt(handle, buffer)       = buffer;
   }

   *handle_ptr = handle;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
_warp_CommWait(warp_Core          core,
               _warp_CommHandle **handle_ptr)
{
   warp_App           app    = _warp_CoreElt(core, app);
   _warp_CommHandle  *handle = *handle_ptr;

   if (handle != NULL)
   {
      warp_Int      request_type = _warp_CommHandleElt(handle, request_type);
      warp_Int      num_requests = _warp_CommHandleElt(handle, num_requests);
      MPI_Request  *requests     = _warp_CommHandleElt(handle, requests);
      MPI_Status   *status       = _warp_CommHandleElt(handle, status);
      void         *buffer       = _warp_CommHandleElt(handle, buffer);

      MPI_Waitall(num_requests, requests, status);
      
      if (request_type == 1) /* recv type */
      {
         warp_Vector  *vector_ptr = _warp_CommHandleElt(handle, vector_ptr);
         
         _warp_CoreFcn(core, bufunpack)(app, buffer, vector_ptr);
      }
      
      _warp_TFree(requests);
      _warp_TFree(status);
      _warp_TFree(handle);
      _warp_TFree(buffer);

      *handle_ptr = NULL;
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Working on all intervals
 *--------------------------------------------------------------------------*/

warp_Int
_warp_UCommInit(warp_Core  core,
                warp_Int   level)
{
   _warp_Grid        **grids       = _warp_CoreElt(core, grids);
   warp_Int            ilower      = _warp_GridElt(grids[level], ilower);
   warp_Int            iupper      = _warp_GridElt(grids[level], iupper);
   warp_Int            cfactor     = _warp_GridElt(grids[level], cfactor);
   warp_Int            ncpoints    = _warp_GridElt(grids[level], ncpoints);
   warp_Vector        *ua          = _warp_GridElt(grids[level], ua);
   warp_Int            recv_index  = -1;
   warp_Int            send_index  = -1;
   _warp_CommHandle   *recv_handle = NULL;
   _warp_CommHandle   *send_handle = NULL;

   if (ilower <= iupper)
   {
      /* Post receive */
      _warp_CommRecvInit(core, level, ilower-1, &ua[-1], &recv_handle);
      recv_index = ilower-1;
      
      /* Only post send if iupper is a C-point, otherwise compute and send later */
      if ( _warp_IsCPoint(iupper, cfactor) )
      {
         _warp_CommSendInit(core, level, iupper, ua[ncpoints-1], &send_handle);
         send_index = -1;
      }
      else
      {
         send_index = iupper;
      }
   }

   _warp_GridElt(grids[level], recv_index)  = recv_index ;
   _warp_GridElt(grids[level], send_index)  = send_index ;
   _warp_GridElt(grids[level], recv_handle) = recv_handle;
   _warp_GridElt(grids[level], send_handle) = send_handle;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Working only on F-pt intervals
 *--------------------------------------------------------------------------*/

warp_Int
_warp_UCommInitF(warp_Core  core,
                 warp_Int   level)
{
   _warp_Grid        **grids       = _warp_CoreElt(core, grids);
   warp_Int            ilower      = _warp_GridElt(grids[level], ilower);
   warp_Int            iupper      = _warp_GridElt(grids[level], iupper);
   warp_Int            cfactor     = _warp_GridElt(grids[level], cfactor);
   warp_Int            ncpoints    = _warp_GridElt(grids[level], ncpoints);
   warp_Vector        *ua          = _warp_GridElt(grids[level], ua);
   warp_Int            recv_index  = -1;
   warp_Int            send_index  = -1;
   _warp_CommHandle   *recv_handle = NULL;
   _warp_CommHandle   *send_handle = NULL;

   if (ilower <= iupper)
   {
      /* Only post receive if ilower is an F-point */
      if ( _warp_IsFPoint(ilower, cfactor) )
      {
         _warp_CommRecvInit(core, level, ilower-1, &ua[-1], &recv_handle);
         recv_index = ilower-1;
      }

      /* Only post send if iupper is a C-point, otherwise compute and send later */
      if ( _warp_IsCPoint(iupper, cfactor) )
      {
         _warp_CommSendInit(core, level, iupper, ua[ncpoints-1], &send_handle);
         send_index = -1;
      }
      else if ( _warp_IsFPoint(iupper+1, cfactor) )
      {
         send_index = iupper;
      }
   }

   _warp_GridElt(grids[level], recv_index)  = recv_index ;
   _warp_GridElt(grids[level], send_index)  = send_index ;
   _warp_GridElt(grids[level], recv_handle) = recv_handle;
   _warp_GridElt(grids[level], send_handle) = send_handle;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Finish up communication
 *--------------------------------------------------------------------------*/

warp_Int
_warp_UCommWait(warp_Core  core,
                warp_Int   level)
{
   _warp_Grid        **grids       = _warp_CoreElt(core, grids);
   _warp_CommHandle   *recv_handle = _warp_GridElt(grids[level], recv_handle);
   _warp_CommHandle   *send_handle = _warp_GridElt(grids[level], send_handle);

   _warp_CommWait(core, &recv_handle);
   _warp_CommWait(core, &send_handle);
   _warp_GridElt(grids[level], recv_index)  = -1;
   _warp_GridElt(grids[level], send_index)  = -1;
   _warp_GridElt(grids[level], recv_handle) = recv_handle;
   _warp_GridElt(grids[level], send_handle) = send_handle;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
_warp_UGetInterval(warp_Core   core,
                   warp_Int    level,
                   warp_Int    interval_index,
                   warp_Int   *flo_ptr,
                   warp_Int   *fhi_ptr,
                   warp_Int   *ci_ptr)
{
   _warp_Grid  **grids   = _warp_CoreElt(core, grids);
   warp_Int      ilower  = _warp_GridElt(grids[level], ilower);
   warp_Int      iupper  = _warp_GridElt(grids[level], iupper);
   warp_Int      clower  = _warp_GridElt(grids[level], clower);
   warp_Int      cupper  = _warp_GridElt(grids[level], cupper);
   warp_Int      cfactor = _warp_GridElt(grids[level], cfactor);

   warp_Int      flo, fhi, ci;

   flo = ilower;
   fhi = iupper;
   ci  = -1;

   if ( _warp_IsCPoint(clower, cfactor) )
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

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Returns a reference to the local u-vector on grid 'level' at point 'index'.
 * If 'index' is not a C-point and within my index range, NULL is returned.
 *--------------------------------------------------------------------------*/

warp_Int
_warp_UGetVectorRef(warp_Core     core,
                    warp_Int      level,
                    warp_Int      index,
                    warp_Vector  *u_ptr)
{
   _warp_Grid        **grids       = _warp_CoreElt(core, grids);
   warp_Int            ilower      = _warp_GridElt(grids[level], ilower);
   warp_Int            iupper      = _warp_GridElt(grids[level], iupper);
   warp_Int            clower      = _warp_GridElt(grids[level], clower);
   warp_Int            cfactor     = _warp_GridElt(grids[level], cfactor);
   warp_Vector        *ua          = _warp_GridElt(grids[level], ua);

   warp_Vector         u;
   warp_Int            ic, iclo;

   u = NULL;
   if ((index >= ilower) && (index <= iupper))
   {
      if ( _warp_IsCPoint(index, cfactor) )
      {
         _warp_MapFineToCoarse(index, cfactor, ic);
         _warp_MapFineToCoarse(clower, cfactor, iclo);
         u = ua[ic-iclo];
      }
   }

   *u_ptr = u;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Stores a reference to the u-vector on grid 'level' at point 'index'.
 * If 'index' is not a C-point and within my index range, nothing is done.
 *--------------------------------------------------------------------------*/

warp_Int
_warp_USetVectorRef(warp_Core    core,
                    warp_Int     level,
                    warp_Int     index,
                    warp_Vector  u)
{
   _warp_Grid        **grids       = _warp_CoreElt(core, grids);
   warp_Int            ilower      = _warp_GridElt(grids[level], ilower);
   warp_Int            iupper      = _warp_GridElt(grids[level], iupper);
   warp_Int            clower      = _warp_GridElt(grids[level], clower);
   warp_Int            cfactor     = _warp_GridElt(grids[level], cfactor);
   warp_Vector        *ua          = _warp_GridElt(grids[level], ua);

   warp_Int            ic, iclo;

   if ((index >= ilower) && (index <= iupper))
   {
      if ( _warp_IsCPoint(index, cfactor) )
      {
         _warp_MapFineToCoarse(index, cfactor, ic);
         _warp_MapFineToCoarse(clower, cfactor, iclo);
         ua[ic-iclo] = u;
      }
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Returns the u-vector on grid 'level' at point 'index'.  If 'index' is my
 * "receive index" (as set by UCommInit(), for example), the u-vector will be
 * received from a neighbor processor.  If 'index' is within my index range and
 * is also a C-point, the saved value of u will be used.  A NULL value is
 * returned otherwise.
 *--------------------------------------------------------------------------*/

warp_Int
_warp_UGetVector(warp_Core     core,
                 warp_Int      level,
                 warp_Int      index,
                 warp_Vector  *u_ptr)
{
   warp_App            app         = _warp_CoreElt(core, app);
   _warp_Grid        **grids       = _warp_CoreElt(core, grids);
   warp_Vector        *ua          = _warp_GridElt(grids[level], ua);
   warp_Int            recv_index  = _warp_GridElt(grids[level], recv_index);
   _warp_CommHandle   *recv_handle = _warp_GridElt(grids[level], recv_handle);

   warp_Vector         u, uu;

   if (index == recv_index)
   {
      /* If a recv was initiated, receive u value from neighbor processor */
      if (recv_index > -1)
      {
         _warp_CommWait(core, &recv_handle);
         _warp_GridElt(grids[level], recv_index)  = -1;
         _warp_GridElt(grids[level], recv_handle) = recv_handle;
         u = ua[-1];
      }
      else
      {
         u = NULL;
      }
   }
   else
   {
      _warp_UGetVectorRef(core, level, index, &uu);
      if (uu != NULL)
      {
         _warp_CoreFcn(core, clone)(app, uu, &u);
      }
      else
      {
         u = NULL;
      }
   }

   *u_ptr = u;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Sets the u-vector on grid 'level' at point 'index'.  If 'index' is my "send
 * index", a send is initiated to a neighbor processor.  If 'index' is within my
 * index range and is also a C-point, the value is saved locally.
 *--------------------------------------------------------------------------*/

warp_Int
_warp_USetVector(warp_Core    core,
                 warp_Int     level,
                 warp_Int     index,
                 warp_Vector  u)
{
   warp_App            app         = _warp_CoreElt(core, app);
   _warp_Grid        **grids       = _warp_CoreElt(core, grids);
   warp_Int            send_index  = _warp_GridElt(grids[level], send_index);
   _warp_CommHandle   *send_handle = _warp_GridElt(grids[level], send_handle);

   warp_Vector         uu;

   if (index == send_index)
   {
      /* Post send to neighbor processor */
      _warp_CommSendInit(core, level, index, u, &send_handle);
      _warp_GridElt(grids[level], send_index)  = -1;
      _warp_GridElt(grids[level], send_handle) = send_handle;
   }

   _warp_UGetVectorRef(core, level, index, &uu);
   if (uu != NULL)
   {
      _warp_CoreFcn(core, free)(app, uu);
   }
   _warp_USetVectorRef(core, level, index, u);

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
_warp_UWriteVector(warp_Core    core,
                   warp_Int     level,
                   warp_Int     index,
                   warp_Vector  u)
{
   warp_App      app = _warp_CoreElt(core, app);
   _warp_Grid  **grids  = _warp_CoreElt(core, grids);
   warp_Int      ilower = _warp_GridElt(grids[level], ilower);
   warp_Float   *ta     = _warp_GridElt(grids[level], ta);

   _warp_CoreFcn(core, write)(app, ta[index-ilower], u);

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Apply Phi
 *--------------------------------------------------------------------------*/

warp_Int
_warp_Phi(warp_Core     core,
          warp_Int      level,
          warp_Int      index,
          warp_Float    accuracy,
          warp_Int      gzero,
          warp_Vector   u,
          warp_Int     *rfactor)
{
   warp_App      app    = _warp_CoreElt(core, app);
   _warp_Grid  **grids  = _warp_CoreElt(core, grids);
   warp_Int      ilower = _warp_GridElt(grids[level], ilower);
   warp_Float   *ta     = _warp_GridElt(grids[level], ta);

   warp_Int      ii;

   ii = index-ilower;
   _warp_CoreFcn(core, phi)(app, ta[ii-1], ta[ii], accuracy, gzero, u, rfactor);

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Integrate one time step
 *--------------------------------------------------------------------------*/

warp_Int
_warp_Step(warp_Core     core,
           warp_Int      level,
           warp_Int      index,
           warp_Float    accuracy,
           warp_Vector   u)
{
   warp_App       app      = _warp_CoreElt(core, app);
   warp_Int      *rfactors = _warp_CoreElt(core, rfactors);
   _warp_Grid   **grids    = _warp_CoreElt(core, grids);
   warp_Int       ilower   = _warp_GridElt(grids[level], ilower);
   warp_Vector   *va       = _warp_GridElt(grids[level], va);
   warp_Vector   *wa       = _warp_GridElt(grids[level], wa);

   warp_Int       rfactor, ii;

   ii = index-ilower;
   if (level == 0)
   {
      _warp_Phi(core, level, index, accuracy, 0, u, &rfactor);
      rfactors[ii] = rfactor;
   }
   else
   {
      _warp_Phi(core, level, index, accuracy, 1, u, &rfactor);
      _warp_CoreFcn(core, sum)(app, 1.0, va[ii], 1.0, u);
      _warp_CoreFcn(core, sum)(app, 1.0, wa[ii], 1.0, u);
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Coarsen in space
 *--------------------------------------------------------------------------*/

warp_Int
_warp_Coarsen(warp_Core     core,
              warp_Int      level,  /* coarse level */
              warp_Int      index,  /* coarse index */
              warp_Vector   fvector,
              warp_Vector  *cvector)
{
   warp_App      app    = _warp_CoreElt(core, app);
   _warp_Grid  **grids  = _warp_CoreElt(core, grids);
   warp_Int      ilower = _warp_GridElt(grids[level], ilower);
   warp_Float   *ta     = _warp_GridElt(grids[level], ta);

   warp_Int      ii;

   if ( _warp_CoreElt(core, coarsen) == NULL )
   {
      /* No spatial coarsening, just clone the fine vector */
      _warp_CoreFcn(core, clone)(app, fvector, cvector);
   }
   else
   {
      /* Call the user's coarsening routine */
      ii = index-ilower;
      _warp_CoreFcn(core, coarsen)(app, ta[ii], ta[ii+1], fvector, cvector);
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Refine in space
 *--------------------------------------------------------------------------*/

warp_Int
_warp_Refine(warp_Core     core,
             warp_Int      level,  /* fine level */
             warp_Int      index,  /* fine index */
             warp_Vector   cvector,
             warp_Vector  *fvector)
{
   warp_App      app    = _warp_CoreElt(core, app);
   _warp_Grid  **grids  = _warp_CoreElt(core, grids);
   warp_Int      ilower = _warp_GridElt(grids[level], ilower);
   warp_Float   *ta     = _warp_GridElt(grids[level], ta);

   warp_Int      ii;

   if ( _warp_CoreElt(core, refine) == NULL )
   {
      /* No spatial refinement, just clone the fine vector */
      _warp_CoreFcn(core, clone)(app, cvector, fvector);
   }
   else
   {
      /* Call the user's refinement routine */
      ii = index-ilower;
      _warp_CoreFcn(core, refine)(app, ta[ii], ta[ii+1], cvector, fvector);
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Create a new grid object
 *--------------------------------------------------------------------------*/

warp_Int
_warp_GridInit(warp_Core     core,
               warp_Int      level,
               warp_Int      ilower,
               warp_Int      iupper,
               _warp_Grid  **grid_ptr)
{
   _warp_Grid   *grid;
   warp_Float   *ta;

   grid = _warp_CTAlloc(_warp_Grid, 1);
   
   _warp_GridElt(grid, level)  = level;
   _warp_GridElt(grid, ilower) = ilower;
   _warp_GridElt(grid, iupper) = iupper;
   
   ta = _warp_CTAlloc(warp_Float,  iupper-ilower+2);
   _warp_GridElt(grid, ta_alloc) = ta;
   _warp_GridElt(grid, ta)       = ta+1;  /* shift */
   
   *grid_ptr = grid;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
_warp_GridDestroy(warp_Core    core,
                  _warp_Grid  *grid)
{
   if (grid)
   {
      warp_App      app      = _warp_CoreElt(core, app);
      warp_Int      ncpoints = _warp_GridElt(grid, ncpoints);
      warp_Vector  *ua       = _warp_GridElt(grid, ua);
      warp_Vector  *ua_alloc = _warp_GridElt(grid, ua_alloc);
      warp_Float   *ta_alloc = _warp_GridElt(grid, ta_alloc);
      warp_Vector  *va_alloc = _warp_GridElt(grid, va_alloc);
      warp_Vector  *wa_alloc = _warp_GridElt(grid, wa_alloc);

      warp_Int      i;

      if (ua_alloc)
      {
         for (i = 0; i < ncpoints; i++)
         {
            if (ua[i])
            {
               _warp_CoreFcn(core, free)(app, ua[i]);
            }
         }
         _warp_TFree(ua_alloc);
      }
      if (ta_alloc)
      {
         _warp_TFree(ta_alloc);
      }
      if (va_alloc)
      {
         _warp_TFree(va_alloc);
      }
      if (wa_alloc)
      {
         _warp_TFree(wa_alloc);
      }

      _warp_TFree(grid);
   }
   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Initialize (and re-initialize) hierarchy
 *--------------------------------------------------------------------------*/

warp_Int
_warp_InitHierarchy(warp_Core    core,
                    _warp_Grid  *fine_grid)
{
   MPI_Comm      comm       = _warp_CoreElt(core, comm);
   warp_Int      max_levels = _warp_CoreElt(core, max_levels);
   warp_Float    tol        = _warp_CoreElt(core, tol);
   warp_Int     *nrels      = _warp_CoreElt(core, nrels);
   warp_Int      nrdefault  = _warp_CoreElt(core, nrdefault);
   warp_Int     *cfactors   = _warp_CoreElt(core, cfactors);
   warp_Int      cfdefault  = _warp_CoreElt(core, cfdefault);
   warp_Int      gupper     = _warp_CoreElt(core, gupper);
   warp_Int      nlevels    = _warp_CoreElt(core, nlevels);
   _warp_Grid  **grids      = _warp_CoreElt(core, grids);

   warp_Int      level;
   warp_Int      ilower, iupper;
   warp_Int      clower, cupper, cfactor, ncpoints;
   warp_Vector  *ua;
   warp_Float   *ta;
   warp_Vector  *va;
   warp_Vector  *wa;

   _warp_Grid   *grid;
   warp_Float   *f_ta;
   warp_Int      i, f_i, f_ilower, clo, chi, gclower, gcupper;

   MPI_Request   request;
   MPI_Status    status;
   warp_Int      rproc, sproc;

   grids[0] = fine_grid;

   /* Do sequential time marching if tolerance is not positive */
   if ((tol <= 0.0) && (max_levels > 1))
   {
      max_levels = 1;
      _warp_CoreElt(core, max_levels) = max_levels;
   }

   /* Allocate space for rfactors */
   ilower = _warp_GridElt(grids[0], ilower);
   iupper = _warp_GridElt(grids[0], iupper);
   _warp_CoreElt(core, rfactors) = _warp_CTAlloc(warp_Int, iupper-ilower+1);

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
      ilower = _warp_GridElt(grid, ilower);
      iupper = _warp_GridElt(grid, iupper);
      if (level > 0)
      {
         va = _warp_CTAlloc(warp_Vector, iupper-ilower+2);
         wa = _warp_CTAlloc(warp_Vector, iupper-ilower+2);
         _warp_GridElt(grid, va_alloc) = va;
         _warp_GridElt(grid, wa_alloc) = wa;
         _warp_GridElt(grid, va)       = va+1;  /* shift */
         _warp_GridElt(grid, wa)       = wa+1;  /* shift */

         /* Copy ta info from level-1 grid */
         ta       = _warp_GridElt(grid, ta);
         f_ilower = _warp_GridElt(grids[level-1], ilower);
         f_ta     = _warp_GridElt(grids[level-1], ta);
         cfactor  = _warp_GridElt(grids[level-1], cfactor);
         for (i = ilower; i <= iupper; i++)
         {
            _warp_MapCoarseToFine(i, cfactor, f_i);
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

      _warp_ProjectInterval(gclower, gcupper, 0, cfactor, &gclower, &gcupper);
      _warp_MapFineToCoarse(gclower, cfactor, gclower);
      _warp_MapFineToCoarse(gcupper, cfactor, gcupper);
      if ( (gclower < gcupper) && (max_levels > level+1) )
      {
         /* Coarsen */
         _warp_ProjectInterval(ilower, iupper, 0, cfactor, &clower, &cupper);
         _warp_MapFineToCoarse(clower, cfactor, clo);
         _warp_MapFineToCoarse(cupper, cfactor, chi);
         ncpoints = chi-clo+1;
         if (ncpoints < 0)
         {
            ncpoints = 0;
         }
         _warp_GridElt(grid, clower)   = clower;
         _warp_GridElt(grid, cupper)   = cupper;
         _warp_GridElt(grid, cfactor)  = cfactor;
         _warp_GridElt(grid, ncpoints) = ncpoints;
         ua = _warp_CTAlloc(warp_Vector, ncpoints+1);
         _warp_GridElt(grid, ua_alloc) = ua;
         _warp_GridElt(grid, ua)       = ua+1;  /* shift */

         _warp_GridInit(core, level+1, clo, chi, &grids[level+1]);
      }
      else
      {
         /* Make this the coarsest level */
         if (ilower == 0)
         {
            ncpoints = 1;
         }
         else
         {
            ncpoints = 0;
         }
         _warp_GridElt(grid, clower)   = ilower;
         _warp_GridElt(grid, cupper)   = 0;
         _warp_GridElt(grid, cfactor)  = gupper+1;
         _warp_GridElt(grid, ncpoints) = ncpoints;
         if (ilower <= iupper)
         {
            ua = _warp_CTAlloc(warp_Vector, ncpoints+1);
            _warp_GridElt(grid, ua_alloc) = ua;
            _warp_GridElt(grid, ua)       = ua+1;  /* shift */
         }

         /* Stop coarsening */
         break;
      }
   }
   nlevels = level+1;
   _warp_CoreElt(core, nlevels) = nlevels;

   /* Communicate ta[-1] information */
   for (level = 0; level < nlevels; level++)
   {
      grid = grids[level];
      ilower = _warp_GridElt(grid, ilower);
      iupper = _warp_GridElt(grid, iupper);
      ta     = _warp_GridElt(grid, ta);

      if (ilower <= iupper)
      {
         /* Post receive */
         _warp_GetProc(core, level, ilower-1, &rproc);
         if (rproc > -1)
         {
            MPI_Irecv(&ta[-1], sizeof(warp_Float), MPI_BYTE,
                      rproc, 0, comm, &request);
         }
         /* Post send */
         _warp_GetProc(core, level, iupper+1, &sproc);
         if (sproc > -1)
         {
            MPI_Send(&ta[iupper-ilower], sizeof(warp_Float), MPI_BYTE,
                     sproc, 0, comm);
         }
         /* Finish receive */
         if (rproc > -1)
         {
            MPI_Wait(&request, &status);
         }
      }
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Set initial guess at C-points
 *--------------------------------------------------------------------------*/

warp_Int
_warp_InitGuess(warp_Core  core,
                warp_Int   level)
{
   warp_App      app     = _warp_CoreElt(core, app);
   _warp_Grid  **grids   = _warp_CoreElt(core, grids);
   warp_Int      ilower  = _warp_GridElt(grids[level], ilower);
   warp_Int      clower  = _warp_GridElt(grids[level], clower);
   warp_Int      cupper  = _warp_GridElt(grids[level], cupper);
   warp_Int      cfactor = _warp_GridElt(grids[level], cfactor);
   warp_Float   *ta      = _warp_GridElt(grids[level], ta);
   warp_Vector  *va      = _warp_GridElt(grids[level], va);

   warp_Vector    u;
   warp_Int       i;

   if (level == 0)
   {
      for (i = clower; i <= cupper; i += cfactor)
      {
         _warp_CoreFcn(core, init)(app, ta[i-ilower], &u);
         _warp_USetVectorRef(core, level, i, u);
      }
   }
   else
   {
      for (i = clower; i <= cupper; i += cfactor)
      {
         _warp_CoreFcn(core, clone)(app, va[i-ilower], &u);
         _warp_USetVectorRef(core, level, i, u);
      }
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Do nu sweeps of F-then-C relaxation
 *--------------------------------------------------------------------------*/

warp_Int
_warp_CFRelax(warp_Core  core,
              warp_Int   level)
{
   warp_App       app      = _warp_CoreElt(core, app);
   warp_Int      *nrels    = _warp_CoreElt(core, nrels);
   _warp_Grid   **grids    = _warp_CoreElt(core, grids);
   warp_Int       ncpoints = _warp_GridElt(grids[level], ncpoints);

   warp_Vector    u;
   warp_Int       flo, fhi, fi, ci;
   warp_Int       nu, nrelax, interval, c_level;
   warp_Float     accuracy;

   if (level == 0)
   {
      accuracy = 0.0;
   }
   else
   {
      accuracy = 1.0;
   }

   nrelax  = nrels[level];
   c_level = level+1;

   for (nu = 0; nu < nrelax; nu++)
   {
      _warp_UCommInit(core, level);

      /* Start from the right-most interval */
      for (interval = ncpoints; interval > -1; interval--)
      {
         _warp_UGetInterval(core, level, interval, &flo, &fhi, &ci);

         if (flo <= fhi)
         {
            _warp_UGetVector(core, level, flo-1, &u);
         }
         else if (ci > 0)
         {
            _warp_UGetVector(core, level, ci-1, &u);
         }

         /* F-relaxation */
         for (fi = flo; fi <= fhi; fi++)
         {
            _warp_Step(core, level, fi, accuracy, u);
            _warp_USetVector(core, level, fi, u);
         }

         /* C-relaxation */
         if (ci > 0)
         {
            _warp_Step(core, level, ci, accuracy, u);
            _warp_USetVector(core, level, ci, u);
         }

         if ((flo <= fhi) && (interval == ncpoints))
         {
            _warp_CoreFcn(core, free)(app, u);
         }
      }
      _warp_UCommWait(core, level);
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * F-Relax on level and restrict to level+1
 *--------------------------------------------------------------------------*/

warp_Int
_warp_FRestrict(warp_Core    core,
                warp_Int     level,
                warp_Float  *rnorm_ptr)
{
   MPI_Comm            comm        = _warp_CoreElt(core, comm);
   warp_App            app         = _warp_CoreElt(core, app);
   _warp_Grid        **grids       = _warp_CoreElt(core, grids);
   warp_Int            cfactor     = _warp_GridElt(grids[level], cfactor);
   warp_Int            ncpoints    = _warp_GridElt(grids[level], ncpoints);
   _warp_CommHandle   *recv_handle = NULL;
   _warp_CommHandle   *send_handle = NULL;

   warp_Int            c_level, c_ilower, c_iupper, c_index, c_i, c_ii;
   warp_Vector         c_u, *c_va, *c_wa;

   warp_Vector         u, r;
   warp_Int            interval, flo, fhi, fi, ci, rfactor;
   warp_Float          rnorm, grnorm, rdot, accuracy;

   c_level  = level+1;
   c_ilower = _warp_GridElt(grids[c_level], ilower);
   c_iupper = _warp_GridElt(grids[c_level], iupper);
   c_va     = _warp_GridElt(grids[c_level], va);
   c_wa     = _warp_GridElt(grids[c_level], wa);

   rnorm = 0.0;

   if (level == 0)
   {
      accuracy = 0.0;
   }
   else
   {
      accuracy = 1.0;
   }

   _warp_UCommInit(core, level);

   /* Start from the right-most interval */
   for (interval = ncpoints; interval > -1; interval--)
   {
      _warp_UGetInterval(core, level, interval, &flo, &fhi, &ci);

      if (flo <= fhi)
      {
         _warp_UGetVector(core, level, flo-1, &r);
      }
      else if (ci > 0)
      {
         _warp_UGetVector(core, level, ci-1, &r);
      }

      /* F-relaxation */
      for (fi = flo; fi <= fhi; fi++)
      {
         _warp_Step(core, level, fi, accuracy, r);
         _warp_USetVector(core, level, fi, r);
      }

      /* Compute residual and restrict */
      if (ci > 0)
      {
         /* Compute residual */
         _warp_Step(core, level, ci, accuracy, r);
         _warp_UGetVectorRef(core, level, ci, &u);
         _warp_CoreFcn(core, sum)(app, -1.0, u, 1.0, r);

         /* Compute rnorm (only on level 0) */
         if (level == 0)
         {
            _warp_CoreFcn(core, dot)(app, r, r, &rdot);
            rnorm += rdot;
         }

         /* Restrict u and residual */
         _warp_MapFineToCoarse(ci, cfactor, c_index);
         _warp_Coarsen(core, c_level, c_index, u, &c_va[c_index-c_ilower]);
         _warp_Coarsen(core, c_level, c_index, r, &c_wa[c_index-c_ilower]);
      }
      else if (ci == 0)
      {
         /* Restrict initial condition to coarse grid */
         _warp_UGetVectorRef(core, level, 0, &u);
         _warp_Coarsen(core, c_level, 0, u, &c_va[0]);
      }

      if ((flo <= fhi) || (ci > 0))
      {
         _warp_CoreFcn(core, free)(app, r);
      }
   }
   _warp_UCommWait(core, level);

   /* Now apply Phi_coarse to update wa values */

   /* Initialize update of c_va[-1] boundary */
   if (c_ilower <= c_iupper)
   {
      _warp_CommRecvInit(core, c_level, c_ilower-1, &c_va[-1],
                         &recv_handle);
      _warp_CommSendInit(core, c_level, c_iupper, c_va[c_iupper-c_ilower],
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
            _warp_CommWait(core, &recv_handle);
         }
         _warp_CoreFcn(core, clone)(app, c_va[c_ii-1], &c_u);
         _warp_Phi(core, c_level, c_i, 1.0, 1, c_u, &rfactor);
         _warp_CoreFcn(core, sum)(app, -1.0, c_u, 1.0, c_wa[c_ii]);
         _warp_CoreFcn(core, free)(app, c_u);
      }
   }
   _warp_CommWait(core, &send_handle);

   /* Compute rnorm (only on level 0) */
   if (level == 0)
   {
      MPI_Allreduce(&rnorm, &grnorm, 1, MPI_DOUBLE, MPI_SUM, comm);
      grnorm = sqrt(grnorm);

      *rnorm_ptr = grnorm;
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * F-Relax on level and interpolate to level-1
 *--------------------------------------------------------------------------*/

warp_Int
_warp_FInterp(warp_Core  core,
              warp_Int   level)
{
   warp_App       app      = _warp_CoreElt(core, app);
   _warp_Grid   **grids    = _warp_CoreElt(core, grids);
   warp_Int       ilower   = _warp_GridElt(grids[level], ilower);
   warp_Int       iupper   = _warp_GridElt(grids[level], iupper);
   warp_Int       ncpoints = _warp_GridElt(grids[level], ncpoints);
   warp_Vector   *ua       = _warp_GridElt(grids[level], ua);
   warp_Vector   *va       = _warp_GridElt(grids[level], va);
   warp_Vector   *wa       = _warp_GridElt(grids[level], wa);

   warp_Int       f_level, f_ilower, f_cfactor, f_index;
   warp_Vector    f_u, f_e;

   warp_Vector    u, e;
   warp_Int       flo, fhi, fi, ci, ii;
   warp_Int       interval;

   f_level   = level-1;
   f_ilower  = _warp_GridElt(grids[f_level], ilower);
   f_cfactor = _warp_GridElt(grids[f_level], cfactor);

   _warp_UCommInitF(core, level);

   /* Start from the right-most interval */
   for (interval = ncpoints; interval > -1; interval--)
   {
      _warp_UGetInterval(core, level, interval, &flo, &fhi, &ci);

      /* Relax and interpolate F-points */
      if (flo <= fhi)
      {
         _warp_UGetVector(core, level, flo-1, &u);
      }
      for (fi = flo; fi <= fhi; fi++)
      {
         _warp_Step(core, level, fi, 1.0, u);
         _warp_USetVector(core, level, fi, u);
         e = va[fi-ilower];
         _warp_CoreFcn(core, sum)(app, 1.0, u, -1.0, e);
         _warp_MapCoarseToFine(fi, f_cfactor, f_index);
         _warp_Refine(core, f_level, f_index, e, &f_e);
         _warp_UGetVectorRef(core, f_level, f_index, &f_u);
         _warp_CoreFcn(core, sum)(app, 1.0, f_e, 1.0, f_u);
         _warp_USetVectorRef(core, f_level, f_index, f_u);
         _warp_CoreFcn(core, free)(app, f_e);
      }
      if (flo <= fhi)
      {
         _warp_CoreFcn(core, free)(app, u);
      }

      /* Interpolate C-points */
      if (ci > 0)
      {
         _warp_UGetVectorRef(core, level, ci, &u);
         e = va[ci-ilower];
         _warp_CoreFcn(core, sum)(app, 1.0, u, -1.0, e);
         _warp_MapCoarseToFine(ci, f_cfactor, f_index);
         _warp_Refine(core, f_level, f_index, e, &f_e);
         _warp_UGetVectorRef(core, f_level, f_index, &f_u);
         _warp_CoreFcn(core, sum)(app, 1.0, f_e, 1.0, f_u);
         _warp_USetVectorRef(core, f_level, f_index, f_u);
         _warp_CoreFcn(core, free)(app, f_e);
      }
   }

   _warp_UCommWait(core, level);

   /* Clean up */
   for (ii = 0; ii < ncpoints; ii++)
   {
      _warp_CoreFcn(core, free)(app, ua[ii]);
      ua[ii] = NULL;
   }
   for (ii = -1; ii <= (iupper-ilower); ii++)
   {
      if (va[ii] != NULL)
      {
         _warp_CoreFcn(core, free)(app, va[ii]);
         va[ii] = NULL;
      }
      if (wa[ii] != NULL)
      {
         _warp_CoreFcn(core, free)(app, wa[ii]);
         wa[ii] = NULL;
      }
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Create a new fine grid based on user refinement factor information, then
 * F-relax and interpolate to the new fine grid and create a new multigrid
 * hierarchy.  In general, this will require load re-balancing as well.
 *
 * RDF: Todo
 *--------------------------------------------------------------------------*/

warp_Int
_warp_FRefine(warp_Core   core,
              warp_Int   *refined_ptr)
{
#if 0 /* Finish writing later */
   warp_Int     *rfactors = _warp_CoreElt(core, rfactors);
   warp_Int      nlevels  = _warp_CoreElt(core, nlevels);
   _warp_Grid  **grids    = _warp_CoreElt(core, grids);
   warp_Int      ilower   = _warp_GridElt(grids[0], ilower);
   warp_Int      iupper   = _warp_GridElt(grids[0], iupper);
   warp_Int      cfactor  = _warp_GridElt(grids[0], cfactor);

   warp_Int      lrefine, refine;

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
      ua = _warp_CTAlloc(warp_Vector, (ilower-iupper+1));

      _warp_UCommInitF(core, 0);

      /* Start from the right-most interval */
      for (interval = ncpoints; interval > -1; interval--)
      {
         _warp_UGetInterval(core, 0, interval, &flo, &fhi, &ci);

         /* Relax and interpolate F-points */
         if (flo <= fhi)
         {
            _warp_UGetVector(core, 0, flo-1, &u);
         }
         for (fi = flo; fi <= fhi; fi++)
         {
            _warp_Step(core, 0, fi, 1.0, u);
            _warp_USetVector(core, 0, fi, u);
            _warp_CoreFcn(core, clone)(app, u, &ua[fi-ilower]);
         }
         if (flo <= fhi)
         {
            _warp_CoreFcn(core, free)(app, u);
         }

         /* Interpolate C-points */
         if (ci > 0)
         {
            _warp_UGetVectorRef(core, 0, ci, &u);
            ua[ci-ilower] = u;
         }
      }

      _warp_UCommWait(core, 0);

      /* Create new finest grid */
      f_cfactor = cfactor;
      _warp_MapCoarseToFine(ilower, f_cfactor, f_ilower);
      _warp_MapCoarseToFine(iupper, f_cfactor, f_iupper);
      f_iupper += f_cfactor-1;
      _warp_GridInit(core, 0, f_ilower, f_iupper, &f_grid);

      /* Change gupper (this currently does not match tstop) */
      gupper = _warp_CoreElt(core, gupper);
      gupper = (gupper+1)*f_cfactor - 1;
      _warp_CoreElt(core, gupper) = gupper;

      /* Set t values */
      f_ta = _warp_GridElt(f_grid, ta);
      ta   = _warp_GridElt(grids[0], ta);
      for (i = ilower; i <= iupper; i++)
      {
         _warp_MapCoarseToFine(i, f_cfactor, f_i);
         tstart = ta[i-ilower];
         /* RDF - START HERE*/
         for (j = 0; j < f_cfactor; j++)
         {
            ta[i-ilower] = f_ta[f_i-f_ilower];
            f_ta[f_i-f_ilower] = tstart + (((warp_Float)i)/ntime)*(tstop-tstart);
         }
      }

      /* Initialize new hierarchy */
      for (level = 0; level < nlevels; level++)
      {
         _warp_GridDestroy(core, grids[level]);
      }
      _warp_InitHierarchy(core, f_grid);
      nlevels = _warp_CoreElt(core, nlevels);

      /* Do the spatial refinement here and set the fine grid unknowns */
      for (i = ilower; i <= iupper; i++)
      {
         _warp_MapCoarseToFine(i, f_cfactor, f_index);
         _warp_Refine(core, 0, f_index, ua[i-ilower], &f_u);
         _warp_USetVectorRef(core, 0, f_index, f_u);
         _warp_CoreFcn(core, free)(app, ua[i-ilower]);
      }

      _warp_TFree(ua);
   }

   *refined_ptr = refine;
#endif

   *refined_ptr = 0;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 * Write out the solution on grid level
 *--------------------------------------------------------------------------*/

warp_Int
_warp_FWrite(warp_Core  core,
             warp_Int   level)
{
   warp_App       app      = _warp_CoreElt(core, app);
   _warp_Grid   **grids    = _warp_CoreElt(core, grids);
   warp_Int       ncpoints = _warp_GridElt(grids[level], ncpoints);

   warp_Vector   u;
   warp_Int      interval, flo, fhi, fi, ci;

   _warp_UCommInitF(core, level);

   /* Start from the right-most interval */
   for (interval = ncpoints; interval > -1; interval--)
   {
      _warp_UGetInterval(core, level, interval, &flo, &fhi, &ci);

      /* Write out F-points */
      if (flo <= fhi)
      {
         _warp_UGetVector(core, level, flo-1, &u);
      }
      for (fi = flo; fi <= fhi; fi++)
      {
         _warp_Step(core, level, fi, 0.0, u);
         _warp_USetVector(core, level, fi, u);
         _warp_UWriteVector(core, level, fi, u);
      }
      if (flo <= fhi)
      {
         _warp_CoreFcn(core, free)(app, u);
      }

      /* Write out C-points */
      if (ci > -1)
      {
         _warp_UGetVectorRef(core, level, ci, &u);
         _warp_UWriteVector(core, level, ci, u);
      }
   }
   _warp_UCommWait(core, level);

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_Init(MPI_Comm              comm_world,
          MPI_Comm              comm,
          warp_Float            tstart,
          warp_Float            tstop,
          warp_Int              ntime,
          warp_App              app,
          warp_PtFcnPhi         phi,
          warp_PtFcnInit        init,
          warp_PtFcnClone       clone,
          warp_PtFcnFree        free,
          warp_PtFcnSum         sum,
          warp_PtFcnDot         dot,
          warp_PtFcnWrite       write,
          warp_PtFcnBufSize     bufsize,
          warp_PtFcnBufPack     bufpack,
          warp_PtFcnBufUnpack   bufunpack,
          warp_Core            *core_ptr)
{
   _warp_Core  *core;
   warp_Int    *nrels;
   warp_Int     level, max_levels = 30;

   core = _warp_CTAlloc(_warp_Core, 1);

   _warp_CoreElt(core, comm_world) = comm_world;
   _warp_CoreElt(core, comm)       = comm;
   _warp_CoreElt(core, tstart)     = tstart;
   _warp_CoreElt(core, tstop)      = tstop;
   _warp_CoreElt(core, ntime)      = ntime;
   _warp_CoreElt(core, app)        = app;

   _warp_CoreElt(core, phi)        = phi;
   _warp_CoreElt(core, init)       = init;
   _warp_CoreElt(core, clone)      = clone;
   _warp_CoreElt(core, free)       = free;
   _warp_CoreElt(core, sum)        = sum;
   _warp_CoreElt(core, dot)        = dot;
   _warp_CoreElt(core, write)      = write;
   _warp_CoreElt(core, bufsize)    = bufsize;
   _warp_CoreElt(core, bufpack)    = bufpack;
   _warp_CoreElt(core, bufunpack)  = bufunpack;
   _warp_CoreElt(core, coarsen)    = NULL;
   _warp_CoreElt(core, refine)     = NULL;

   _warp_CoreElt(core, max_levels) = max_levels;
   _warp_CoreElt(core, tol)        = 1.0e-09;
   _warp_CoreElt(core, rtol)       = 0;

   nrels = _warp_TAlloc(warp_Int, max_levels);
   for (level = 0; level < max_levels; level++)
   {
      nrels[level] = -1;
   }
   _warp_CoreElt(core, nrels)      = nrels;
   _warp_CoreElt(core, nrdefault)  = 1;

   _warp_CoreElt(core, cfactors)   = _warp_CTAlloc(warp_Int, max_levels);
   _warp_CoreElt(core, cfdefault)  = 2;

   _warp_CoreElt(core, max_iter)   = 100;
   _warp_CoreElt(core, max_iter)   = 100;
   _warp_CoreElt(core, niter)      = 0;
   _warp_CoreElt(core, rnorm)      = 0.0;
   _warp_CoreElt(core, fmg)        = 0;

   _warp_CoreElt(core, gupper)     = ntime;

   _warp_CoreElt(core, rfactors)   = NULL;

   _warp_CoreElt(core, nlevels)    = 0;
   _warp_CoreElt(core, grids)      = _warp_CTAlloc(_warp_Grid *, max_levels);

   *core_ptr = core;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_Drive(warp_Core  core)
{
   warp_Float    tstart   = _warp_CoreElt(core, tstart);
   warp_Float    tstop    = _warp_CoreElt(core, tstop);
   warp_Int      ntime    = _warp_CoreElt(core, ntime);
   warp_Float    tol      = _warp_CoreElt(core, tol);
   warp_Int      rtol     = _warp_CoreElt(core, rtol);
   warp_Int      fmg      = _warp_CoreElt(core, fmg);
   warp_Int      max_iter = _warp_CoreElt(core, max_iter);

   warp_Int      nlevels, iter;
   warp_Float    rnorm;
   warp_Int      ilower, iupper;
   warp_Float   *ta;
   warp_Int      level, fmglevel, down, done, i, refined;
   _warp_Grid   *grid;

   MPI_Comm     comm_world = _warp_CoreElt(core, comm_world);
   warp_Int     myid;
   MPI_Comm_rank(comm_world, &myid);

   level = 0;
   rnorm = 0.0;

   /* Create fine grid */
   _warp_GetDistribution(core, &ilower, &iupper);
   _warp_GridInit(core, 0, ilower, iupper, &grid);

   /* Set t values */
   ta = _warp_GridElt(grid, ta);
   for (i = ilower; i <= iupper; i++)
   {
      ta[i-ilower] = tstart + (((warp_Float)i)/ntime)*(tstop-tstart);
   }

   /* Create a grid hierarchy */
   _warp_InitHierarchy(core, grid);
   nlevels = _warp_CoreElt(core, nlevels);

   /* Set initial values at C-points */
   _warp_InitGuess(core, 0);

   /* Set cycling variables */
   fmglevel = 0;
   if (fmg)
   {
      fmglevel = nlevels-1;
   }
   down = 1;
   done = 0;
   if ((nlevels <= 1) || (tol <= 0.0))
   {
      /* Just do sequential time marching */
      done = 1;
   }

   iter = 0;
   if ( myid == 0 )
   {
      printf("initial stopping tolerance = %e\n", tol);
   }
   while (!done)
   {
      /* Down cycle */

      if (down)
      {
         if (level < (nlevels-1))
         {
            /* CF-relaxation */
            _warp_CFRelax(core, level);

            /* F-relax then restrict */
            _warp_FRestrict(core, level, &rnorm);
            /* Set initial guess on next coarser level */
            _warp_InitGuess(core, level+1);

            /* Adjust tolerance */
            if (rtol && (level == 0) && (iter == 0))
            {
               tol *= rnorm;
            }

            level++;
         }
         else
         {
            /* Coarsest grid - solve on the up-cycle */

            down = 0;
         }
      }

      /* Up cycle */

      if (!down)
      {
         if (level > 0)
         {
            if (level >= fmglevel)
            {
               /* F-relax then interpolate */
               _warp_FInterp(core, level);
               
               level--;
            }
            else
            {
               fmglevel--;
               down = 1;
            }
         }
         else
         {
            /* Finest grid - refine grid if desired, else check convergence */
            _warp_FRefine(core, &refined);

            if (refined)
            {
               nlevels = _warp_CoreElt(core, nlevels);
            }
            else
            {
               /* Note that this residual is based on an earlier iterate */
               if ( myid == 0 )
               {
                  printf("|| r_%d || = %e\n", iter, rnorm);
               }
               if ((rnorm < tol) || (iter == max_iter-1))
               {
                  done = 1;
               }
            }

            iter++;

            if (fmg)
            {
               fmglevel = nlevels-1;
            }
            down = 1;
         }
      }
   }

   if ( myid == 0 )
   {
      printf("modified stopping tolerance = %e\n", tol);
   }

   /* F-relax and write solution to file */
   _warp_FWrite(core, 0);

   _warp_CoreElt(core, niter) = iter;
   _warp_CoreElt(core, rnorm) = rnorm;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_Destroy(warp_Core  core)
{
   if (core)
   {
      warp_Int      nlevels    = _warp_CoreElt(core, nlevels);
      _warp_Grid  **grids      = _warp_CoreElt(core, grids);
      warp_Int      level;

      _warp_TFree(_warp_CoreElt(core, nrels));
      _warp_TFree(_warp_CoreElt(core, cfactors));
      _warp_TFree(_warp_CoreElt(core, rfactors));
      for (level = 0; level < nlevels; level++)
      {
         _warp_GridDestroy(core, grids[level]);
      }
      _warp_TFree(grids);

      _warp_TFree(core);
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_PrintStats(warp_Core  core)
{
   MPI_Comm     comm_world = _warp_CoreElt(core, comm_world);
   warp_Float   tstart     = _warp_CoreElt(core, tstart);
   warp_Float   tstop      = _warp_CoreElt(core, tstop);
   warp_Int     ntime      = _warp_CoreElt(core, ntime);
   warp_Int     max_levels = _warp_CoreElt(core, max_levels);
   warp_Int    *nrels      = _warp_CoreElt(core, nrels);
   warp_Int     nrdefault  = _warp_CoreElt(core, nrdefault);
   warp_Float   tol        = _warp_CoreElt(core, tol);
   warp_Int     rtol       = _warp_CoreElt(core, rtol);
   /*warp_Int    *cfactors   = _warp_CoreElt(core, cfactors);*/
   warp_Int     cfdefault  = _warp_CoreElt(core, cfdefault);
   warp_Int     max_iter   = _warp_CoreElt(core, max_iter);
   warp_Int     niter      = _warp_CoreElt(core, niter);
   warp_Float   rnorm      = _warp_CoreElt(core, rnorm);
   warp_Int     nlevels    = _warp_CoreElt(core, nlevels);

   warp_Int     myid;

   MPI_Comm_rank(comm_world, &myid);
   if ( myid == 0 )
   {
      printf("\n");
      printf("  start time = %e\n", tstart);
      printf("  stop time  = %e\n", tstop);
      printf("  time steps = %d\n", ntime);
      printf("\n");
      printf("  max number of levels = %d\n", max_levels);
      printf("  number of levels     = %d\n", nlevels);
      printf("  coarsening factor    = %d\n", cfdefault);
      printf("  num F-C relaxations  = %d\n", nrdefault);
      printf("  num rels on level 0  = %d\n", nrels[0]);
      printf("  stopping tolerance   = %e\n", tol);
      printf("  relative tolerance?  = %d\n", rtol);
      printf("  max iterations       = %d\n", max_iter);
      printf("  iterations           = %d\n", niter);
      printf("  residual norm        = %e\n", rnorm);
      printf("\n");
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetMaxLevels(warp_Core  core,
                  warp_Int   max_levels)
{
   _warp_CoreElt(core, max_levels) = max_levels;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetAbsTol(warp_Core   core,
               warp_Float  tol)
{
   _warp_CoreElt(core, tol) = tol;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetRelTol(warp_Core   core,
               warp_Float  tol)
{
   _warp_CoreElt(core, tol)  = tol;
   _warp_CoreElt(core, rtol) = 1;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetNRelax(warp_Core  core,
               warp_Int   level,
               warp_Int   nrelax)
{
   warp_Int  *nrels = _warp_CoreElt(core, nrels);

   if (level < 0)
   {
      /* Set default value */
      _warp_CoreElt(core, nrdefault) = nrelax;
   }
   else
   {
      /* Set factor on specified level */
      nrels[level] = nrelax;
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetCFactor(warp_Core  core,
                warp_Int   level,
                warp_Int   cfactor)
{
   warp_Int  *cfactors = _warp_CoreElt(core, cfactors);

   if (level < 0)
   {
      /* Set default value */
      _warp_CoreElt(core, cfdefault) = cfactor;
   }
   else
   {
      /* Set factor on specified level */
      cfactors[level] = cfactor;
   }

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetMaxIter(warp_Core  core,
                warp_Int   max_iter)
{
   _warp_CoreElt(core, max_iter) = max_iter;

   return _warp_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

warp_Int
warp_SetFMG(warp_Core  core)
{
   _warp_CoreElt(core, fmg) = 1;

   return _warp_error_flag;
}

