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
 *----------------------------------------------------------------------------*/

braid_Int
_braid_CommRecvInit(braid_Core           core,
                    braid_Int            level,
                    braid_Int            index,
                    braid_BaseVector    *vector_ptr,
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

   _braid_GetProc(core, level, index, &proc);
   if (proc > -1)
   {
      handle = _braid_TAlloc(_braid_CommHandle, 1);

      /* Allocate buffer through user routine */
      _braid_BufferStatusInit( 0, 0, bstatus );
      _braid_BaseBufSize(core, app,  &size, bstatus);
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


   _braid_GetProc(core, level, index+1, &proc);
   if (proc > -1)
   {
      handle = _braid_TAlloc(_braid_CommHandle, 1);

      /* Allocate buffer through user routine */
      _braid_BufferStatusInit( 0, 0, bstatus );
      _braid_BaseBufSize(core, app,  &size, bstatus);
      buffer = malloc(size);

      /* Store the receiver rank in the status */
      _braid_StatusElt(bstatus, send_recv_rank) = proc;

      /* Note that bufpack may return a size smaller than bufsize */
      _braid_StatusElt(bstatus, size_buffer) = size;
      _braid_BaseBufPack(core, app,  vector, buffer, bstatus);
      size = _braid_StatusElt( bstatus, size_buffer );

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

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

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
      braid_BufferStatus bstatus  = (braid_BufferStatus)core;

      MPI_Waitall(num_requests, requests, status);

      if (request_type == 1) /* recv type */
      {
         _braid_BufferStatusInit( 0, 0, bstatus );
         braid_BaseVector  *vector_ptr = _braid_CommHandleElt(handle, vector_ptr);

         /* Store the sender rank the bufferStatus */
         _braid_StatusElt(bstatus, send_recv_rank ) = status->MPI_SOURCE;

         _braid_BaseBufUnpack(core, app,  buffer, vector_ptr, bstatus);
      }

      _braid_TFree(requests);
      _braid_TFree(status);
      _braid_TFree(handle);
      _braid_TFree(buffer);

      *handle_ptr = NULL;
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_TriCommInit(braid_Core     core,
                   braid_Int      level,
                   braid_Int     *nrequests_ptr,
                   MPI_Request  **requests_ptr,
                   MPI_Status   **statuses_ptr,
                   void        ***buffers_ptr)
{
   MPI_Comm           comm    = _braid_CoreElt(core, comm);
   braid_App          app     = _braid_CoreElt(core, app);
   _braid_Grid      **grids   = _braid_CoreElt(core, grids);
   braid_Int          ilower  = _braid_GridElt(grids[level], ilower);
   braid_Int          iupper  = _braid_GridElt(grids[level], iupper);
   braid_BufferStatus bstatus = (braid_BufferStatus)core;

   braid_Int          nrequests;
   MPI_Request       *requests;
   MPI_Status        *statuses;
   void             **buffers;
   braid_BaseVector   u;
   braid_Int          k, index, proc, size;

   requests = _braid_CTAlloc(MPI_Request, 4); /* upper bound */
   statuses = _braid_CTAlloc(MPI_Status,  4); /* upper bound */
   buffers  = _braid_CTAlloc(void *,      4); /* upper bound */

   nrequests = 0;

   /* post receives */
   _braid_BufferStatusInit(1, 0, bstatus);
   for (k = 0; k < 2; k++)
   {
      switch(k)
      {
         case 0: _braid_GetProc(core, level, ilower-1, &proc); break;
         case 1: _braid_GetProc(core, level, iupper+1, &proc); break;
      }

      if (proc > -1)
      {
         _braid_BaseBufSize(core, app, &size, bstatus);
         buffers[nrequests] = malloc(size);
         MPI_Irecv(buffers[nrequests], size, MPI_BYTE, proc, 0, comm, &requests[nrequests]);
         nrequests++;
      }
   }

   /* post sends */
   _braid_BufferStatusInit(0, 0, bstatus);
   for (k = 0; k < 2; k++)
   {
      switch(k)
      {
         case 0: index = ilower; _braid_GetProc(core, level, ilower-1, &proc); break;
         case 1: index = iupper; _braid_GetProc(core, level, iupper+1, &proc); break;
      }

      if (proc > -1)
      {
         _braid_UGetVectorRef(core, level, index, &u);
         _braid_BaseBufSize(core, app, &size, bstatus);
         buffers[nrequests] = malloc(size);
         _braid_StatusElt(bstatus, size_buffer) = size;
         _braid_BaseBufPack(core, app, u, buffers[nrequests], bstatus);
         size = _braid_StatusElt(bstatus, size_buffer);
         MPI_Isend(buffers[nrequests], size, MPI_BYTE, proc, 0, comm, &requests[nrequests]);
         nrequests++;
      }
   }

   *nrequests_ptr = nrequests;
   *requests_ptr  = requests;
   *statuses_ptr  = statuses;
   *buffers_ptr   = buffers;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_TriCommWait(braid_Core     core,
                   braid_Int      level,
                   braid_Int      nrequests,
                   MPI_Request  **requests_ptr,
                   MPI_Status   **statuses_ptr,
                   void        ***buffers_ptr)
{
   braid_App          app     = _braid_CoreElt(core, app);
   _braid_Grid      **grids   = _braid_CoreElt(core, grids);
   braid_Int          ilower  = _braid_GridElt(grids[level], ilower);
   braid_Int          iupper  = _braid_GridElt(grids[level], iupper);
   braid_BufferStatus bstatus = (braid_BufferStatus)core;

   MPI_Request       *requests = *requests_ptr;
   MPI_Status        *statuses = *statuses_ptr;
   void             **buffers  = *buffers_ptr;
   braid_BaseVector   u;
   braid_Int          k, index, proc;

   MPI_Waitall(nrequests, requests, statuses);

   nrequests = 0;

   /* unpack receives */
   _braid_BufferStatusInit(1, 0, bstatus);
   for (k = 0; k < 2; k++)
   {
      switch(k)
      {
         case 0: index = ilower-1; _braid_GetProc(core, level, ilower-1, &proc); break;
         case 1: index = iupper+1; _braid_GetProc(core, level, iupper+1, &proc); break;
      }

      if (proc > -1)
      {
         _braid_UGetVectorRef(core, level, index, &u);
         if (u != NULL)
         {
            _braid_BaseFree(core, app, u);
         }
         _braid_BaseBufUnpack(core, app, buffers[nrequests], &u, bstatus);

         nrequests++;
      }
   }

   _braid_TFree(requests);
   _braid_TFree(statuses);
   for (k = 0; k < 4; k++)
   {
      _braid_TFree(buffers[k]);
   }
   _braid_TFree(buffers);

   return _braid_error_flag;
}

