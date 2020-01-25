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
#include "_util.h"

/*----------------------------------------------------------------------------
 * Set the residual norm for iteration iter.  If iter < 0, set the rnorm for the
 * last iteration minus |iter|-1.  Also set the initial residual norm.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_SetRNorm(braid_Core  core,
                braid_Int   iter,
                braid_Real  rnorm)
{
   braid_Real  *rnorms   = _braid_CoreElt(core, rnorms);
   braid_Int    max_iter = _braid_CoreElt(core, max_iter);
   braid_Int    skip     = _braid_CoreElt(core, skip);
   braid_Int    k;

   k = iter;
   if (iter < 0)
   {
      k = _braid_CoreElt(core, niter) + 1 + iter;
   }

   if ((k > -1) && (k <= max_iter)) 
   {
      rnorms[k] = rnorm;

      /* Set initial residual norm if not already set */
      if ((k == 0) || ((k==1) && (skip == 1)))
      {
         if ( _braid_CoreElt(core, rnorm0) == braid_INVALID_RNORM )
         {
            _braid_CoreElt(core, rnorm0) = rnorm;
         }
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Get the residual norm for iteration iter.  If iter < 0, get the rnorm for the
 * last iteration minus |iter|-1.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetRNorm(braid_Core  core,
                braid_Int   iter,
                braid_Real *rnorm_ptr)
{
   braid_Real  *rnorms   = _braid_CoreElt(core, rnorms);
   braid_Int    max_iter = _braid_CoreElt(core, max_iter);
   braid_Int    k;

   /* Initialize to invalid value */
   *rnorm_ptr = braid_INVALID_RNORM;

   k = iter;
   if (iter < 0)
   {
      k = _braid_CoreElt(core, niter) + 1 + iter;
   }

   if ((k > -1) && (k <= max_iter)) 
   {
      *rnorm_ptr = rnorms[k];
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Same as SetRNorm, but sets full residual norm
 *----------------------------------------------------------------------------*/

braid_Int
_braid_SetFullRNorm(braid_Core  core,
                    braid_Int   iter,
                    braid_Real  rnorm)
{
   braid_Real  *rnorms   = _braid_CoreElt(core, full_rnorms);
   braid_Int    max_iter = _braid_CoreElt(core, max_iter);
   braid_Int    skip     = _braid_CoreElt(core, skip);
   braid_Int    k;

   if (rnorms == NULL)
   {
      return _braid_error_flag;
   }

   k = iter;
   if (iter < 0)
   {
      k = _braid_CoreElt(core, niter) + 1 + iter;
   }

   if ((k > -1) && (k <= max_iter)) 
   {
      rnorms[k] = rnorm;

      /* Set initial residual norm if not already set */
      if ((k == 0) || ((k==1) && (skip == 1)))
      {
         if ( _braid_CoreElt(core, full_rnorm0) == braid_INVALID_RNORM )
         {
            _braid_CoreElt(core, full_rnorm0) = rnorm;
         }
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Same as GetRNorm, but gets full residual norm
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetFullRNorm(braid_Core  core,
                    braid_Int   iter,
                    braid_Real *rnorm_ptr)
{
   braid_Real  *rnorms   = _braid_CoreElt(core, full_rnorms);
   braid_Int    max_iter = _braid_CoreElt(core, max_iter);
   braid_Int    k;

   /* Initialize to invalid value */
   *rnorm_ptr = braid_INVALID_RNORM;

   if (rnorms == NULL)
   {
      return _braid_error_flag;
   }

   k = iter;
   if (iter < 0)
   {
      k = _braid_CoreElt(core, niter) + 1 + iter;
   }

   if ((k > -1) && (k <= max_iter)) 
   {
      *rnorm_ptr = rnorms[k];
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_ComputeFullRNorm(braid_Core  core,
                        braid_Int   level,
                        braid_Real *return_rnorm)
{
   MPI_Comm           comm        = _braid_CoreElt(core, comm);
   braid_App          app         = _braid_CoreElt(core, app);
   braid_Real         tol         = _braid_CoreElt(core, tol);
   braid_Int          iter        = _braid_CoreElt(core, niter);
   _braid_Grid      **grids       = _braid_CoreElt(core, grids);
   braid_StepStatus   status      = (braid_StepStatus)core;
   braid_Int          nrefine     = _braid_CoreElt(core, nrefine);
   braid_Int          ncpoints    = _braid_GridElt(grids[level], ncpoints);
   braid_Int          gupper      = _braid_CoreElt(core, gupper);
   braid_Int          tnorm       = _braid_CoreElt(core, tnorm);
   braid_Real        *ta          = _braid_GridElt(grids[level], ta);
   braid_Int          ilower      = _braid_GridElt(grids[level], ilower);
   _braid_CommHandle *send_handle;
   braid_Int          send_index;

   braid_Int         flo, fhi, fi, ci, ii, interval;
   braid_Real        rnorm_temp, rnorm = 0, global_rnorm = 0;
   braid_BaseVector  u, r;

   _braid_UCommInit(core, level);

   /* Start from the right-most interval. */
   for (interval = ncpoints; interval > -1; interval--)
   {
      _braid_GetInterval(core, level, interval, &flo, &fhi, &ci);
      if (flo <= fhi)
      {
         _braid_UGetVector(core, level, flo-1, &u);
      }
      else if (ci > _braid_CoreElt(core, initiali))
      {
         _braid_UGetVector(core, level, ci-1, &u);
      }

      /* Generate F-points and get residual. */
      for (fi = flo; fi <= fhi; fi++)
      {
         _braid_BaseClone(core, app,  u, &r);
         _braid_Step(core, level, fi, NULL, u);

         /* Update local processor norm. */
         ii = fi-ilower;
         _braid_StepStatusInit(ta[ii-1], ta[ii], fi-1, tol, iter, level, nrefine, gupper, status);
         _braid_BaseFullResidual(core, app, u, r, status);
         _braid_BaseSpatialNorm(core, app,  r, &rnorm_temp); 
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
            _braid_GridElt(grids[level], send_index)  = _braid_SendIndexNull;
            _braid_GridElt(grids[level], send_handle) = send_handle;
         }
         _braid_BaseFree(core, app,  r);
      }
      /* Residual from C-point. */
      if (ci > _braid_CoreElt(core, initiali))
      {
         /* Update local processor norm. */
         ii = ci-ilower;
         _braid_StepStatusInit(ta[ii-1], ta[ii], ci-1, tol, iter, level, nrefine, gupper, status);
         _braid_UGetVector(core, level, ci, &r);
         _braid_BaseFullResidual(core, app, r, u, status);
         _braid_BaseSpatialNorm(core, app,  u, &rnorm_temp);

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
         _braid_BaseFree(core, app,  r);
      }

      if ((flo <= fhi) || (ci > _braid_CoreElt(core, initiali)))
      {
         _braid_BaseFree(core, app,  u);
      }
   }

   _braid_UCommWait(core, level);

   /* Compute global residual norm. */
   if(tnorm == 1)       /* one-norm reduction */
   {  
      MPI_Allreduce(&rnorm, &global_rnorm, 1, braid_MPI_REAL, MPI_SUM, comm);
   }
   else if(tnorm == 3)  /* inf-norm reduction */
   {  
      MPI_Allreduce(&rnorm, &global_rnorm, 1, braid_MPI_REAL, MPI_MAX, comm);
   }
   else                 /* default two-norm reduction */
   {  
      MPI_Allreduce(&rnorm, &global_rnorm, 1, braid_MPI_REAL, MPI_SUM, comm);
      global_rnorm = sqrt(global_rnorm);
   }

   *return_rnorm = global_rnorm;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Print the residual norm at ever C-point for debugging purposes 
 *----------------------------------------------------------------------------*/

braid_Int
_braid_PrintSpatialNorms(braid_Core    core,
                         braid_Real   *rnorms,     /* This processor's local residual norms at C-points */
                         braid_Int     n)          /* Length of the rnorms array */
{
   MPI_Comm       comm       = _braid_CoreElt(core, comm);
   MPI_Comm       comm_world = _braid_CoreElt(core, comm_world);
   _braid_Grid  **grids      = _braid_CoreElt(core, grids);
   braid_Int      cfactor    = _braid_GridElt(grids[0], cfactor);
   braid_Int      gupper     = _braid_CoreElt(core, gupper);

   braid_Int      g_ncpoints = ceil( ((braid_Real) gupper) / ((braid_Real) cfactor )) + 1;
   
   braid_Int     *recvcounts = NULL;
   braid_Real    *recvbuf;
   braid_Int     *displs;
   braid_Int      i, p, myid_t, myid_world, my_root_global_rank;

   MPI_Comm_size(comm, &p);
   MPI_Comm_rank(comm_world, &myid_world);
   MPI_Comm_rank(comm, &myid_t);

   /* We need to know all the processor's belonging to the temporal
    * communicator of global rank 0.  It is only these processors that are
    * involved with printing output. */
   my_root_global_rank = myid_world;
   MPI_Bcast(&my_root_global_rank, 1, braid_MPI_INT, 0, comm);

   if(my_root_global_rank == 0)
   {
      if(myid_t == 0)
      {
         recvbuf = _braid_CTAlloc(braid_Real, g_ncpoints);
         recvcounts = _braid_CTAlloc(braid_Int, p);
         displs = _braid_CTAlloc(braid_Int, p);
      }

      /* Rank 0 gather's every processor's number of C-points, which forms the
       * displacements (displs) for the Gatherv call below */ 
      MPI_Gather(&n, 1, braid_MPI_INT, recvcounts, 1, braid_MPI_INT, 0, comm);

      if(myid_t == 0)
      {
         displs[0] = 0;
         for(i = 1; i < p; i++)
         {
            displs[i] = displs[i-1] + recvcounts[i-1];
         }
      }

      /* Gather over comm */
      MPI_Gatherv(rnorms, n, braid_MPI_REAL, recvbuf, recvcounts, displs,
                  braid_MPI_REAL, 0, comm);

      if(myid_t == 0)
      {
         for(i = 0; i < g_ncpoints; i++){
            _braid_printf("  Braid:  time step: %6d, rnorm: %1.2e\n", i*cfactor, recvbuf[i]);
         }
      }

      if(myid_t == 0)
      {
         _braid_TFree(recvbuf);
         _braid_TFree(recvcounts);
         _braid_TFree(displs);
      }
   }

   return _braid_error_flag;
}

