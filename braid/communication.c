/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
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
#include "util.h"

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_CommRecvInit(braid_Core           core,
                    braid_Int            level,
                    braid_Int            index,
                    braid_BaseVector     *vector_ptr,
                    _braid_CommHandle  **handle_ptr)
{
   MPI_Comm            comm = _braid_CoreElt(core, comm);
   braid_App           app  = _braid_CoreElt(core, app);
   _braid_CommHandle  *handle = NULL;
   void               *buffer;
   MPI_Request        *requests;
   MPI_Status         *status;
   braid_Int           proc, size, num_requests;
   braid_BufferStatus bstatus = (braid_BufferStatus)core;
   braid_Real         timer   = 0.0;

   _braid_GetProc(core, level, index, &proc);
   if (proc > -1)
   {
      handle = _braid_TAlloc(_braid_CommHandle, 1);

      /* Allocate buffer through user routine */
      _braid_BufferStatusInit(0, index, level, 0, bstatus);
      _braid_BaseBufSize(core, app,  &size, bstatus);
      _braid_BaseBufAlloc(core, app, &buffer, size, bstatus);

      num_requests = 1;
      requests = _braid_CTAlloc(MPI_Request, num_requests);
      status   = _braid_CTAlloc(MPI_Status, num_requests);
      timer = _braid_MPI_Wtime(core);
      MPI_Irecv(buffer, size, MPI_BYTE, proc, 0, comm, &requests[0]);
      _braid_CoreElt(core, timer_MPI_recv) += _braid_MPI_Wtime(core) - timer;

      _braid_CommHandleElt(handle, request_type) = 1; /* recv type = 1 */
      _braid_CommHandleElt(handle, num_requests) = num_requests;
      _braid_CommHandleElt(handle, index)        = index;
      _braid_CommHandleElt(handle, level)        = level;
      _braid_CommHandleElt(handle, requests)     = requests;
      _braid_CommHandleElt(handle, status)       = status;
      _braid_CommHandleElt(handle, buffer)       = buffer;
      _braid_CommHandleElt(handle, vector_ptr)   = vector_ptr;
   }

   *handle_ptr = handle;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_CommSendInit(braid_Core           core,
                    braid_Int            level,
                    braid_Int            index,
                    braid_BaseVector     vector,
                    _braid_CommHandle  **handle_ptr)
{
   MPI_Comm            comm = _braid_CoreElt(core, comm);
   braid_App           app  = _braid_CoreElt(core, app);
   _braid_CommHandle  *handle = NULL;
   void               *buffer;
   MPI_Request        *requests;
   MPI_Status         *status;
   braid_Int           proc, size, num_requests;
   braid_BufferStatus  bstatus   = (braid_BufferStatus)core;
   braid_Real          timer     = 0.0;
   

   _braid_GetProc(core, level, index+1, &proc);
   if (proc > -1)
   {
      handle = _braid_TAlloc(_braid_CommHandle, 1);

      /* Allocate buffer through user routine */
      _braid_BufferStatusInit(0, index, level, 0, bstatus);
      _braid_BaseBufSize(core, app,  &size, bstatus);
      _braid_BaseBufAlloc(core, app, &buffer, size, bstatus);

      /* Store the receiver rank in the status */
      _braid_StatusElt(bstatus, send_recv_rank) = proc;

      /* Note that bufpack may return a size smaller than bufsize */
      _braid_BaseBufPack(core, app,  vector, buffer, bstatus);
      size = _braid_StatusElt( bstatus, size_buffer );

      num_requests = 1;
      requests = _braid_CTAlloc(MPI_Request, num_requests);
      status   = _braid_CTAlloc(MPI_Status, num_requests);
      timer = _braid_MPI_Wtime(core);
      MPI_Isend(buffer, size, MPI_BYTE, proc, 0, comm, &requests[0]);
      _braid_CoreElt(core, timer_MPI_send) += _braid_MPI_Wtime(core) - timer;

      _braid_CommHandleElt(handle, request_type) = 0; /* send type = 0 */
      _braid_CommHandleElt(handle, num_requests) = num_requests;
      _braid_CommHandleElt(handle, index)        = index;
      _braid_CommHandleElt(handle, level)        = level;
      _braid_CommHandleElt(handle, requests)     = requests;
      _braid_CommHandleElt(handle, status)       = status;
      _braid_CommHandleElt(handle, buffer)       = buffer;
   }

   *handle_ptr = handle;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_CommWait(braid_Core          core,
                _braid_CommHandle **handle_ptr)
{
   braid_App           app    = _braid_CoreElt(core, app);
   _braid_CommHandle  *handle = *handle_ptr;
   braid_Real          timer     = 0.0;

   if (handle != NULL)
   {
      braid_Int      request_type = _braid_CommHandleElt(handle, request_type);
      braid_Int      num_requests = _braid_CommHandleElt(handle, num_requests);
      braid_Int      level        = _braid_CommHandleElt(handle, level);
      braid_Int      index        = _braid_CommHandleElt(handle, index);
      MPI_Request   *requests     = _braid_CommHandleElt(handle, requests);
      MPI_Status    *status       = _braid_CommHandleElt(handle, status);
      void          *buffer       = _braid_CommHandleElt(handle, buffer);
      braid_BufferStatus bstatus  = (braid_BufferStatus)core;

      // Wait on all communication, timing result for coarsest-grid separately 
      timer = _braid_MPI_Wtime(core);
      MPI_Waitall(num_requests, requests, status);
      if (_braid_CoreElt(core, level) == -11){// (_braid_CoreElt(core, nlevels)-1) ) {
         _braid_CoreElt(core, timer_MPI_wait_coarse) += _braid_MPI_Wtime(core) - timer;
      }
      else {
         _braid_CoreElt(core, timer_MPI_wait) += _braid_MPI_Wtime(core) - timer;
      }
      
      if (request_type == 1) /* recv type */
      {
         _braid_BufferStatusInit(0, index, level, 0, bstatus);
         braid_BaseVector  *vector_ptr = _braid_CommHandleElt(handle, vector_ptr);
         
         /* Store the sender rank the bufferStatus */
         _braid_StatusElt(bstatus, send_recv_rank ) = status->MPI_SOURCE;
         
         _braid_BaseBufUnpack(core, app,  buffer, vector_ptr, bstatus);
      }
      
      _braid_TFree(requests);
      _braid_TFree(status);
      _braid_TFree(handle);
      _braid_BaseBufFree(core, app,  &buffer);

      *handle_ptr = NULL;
   }

   return _braid_error_flag;
}

