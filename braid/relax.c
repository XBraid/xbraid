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
 * Do nu sweeps of F-then-C relaxation
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FCRelax(braid_Core  core,
               braid_Int   level)
{
   braid_App       app      = _braid_CoreElt(core, app);
   braid_Int       nlevels  = _braid_CoreElt(core, nlevels);
   braid_Int      *nrels    = _braid_CoreElt(core, nrels);
   braid_Real     *CWts     = _braid_CoreElt(core, CWts);
   _braid_Grid   **grids    = _braid_CoreElt(core, grids);
   braid_Int       ncpoints = _braid_GridElt(grids[level], ncpoints);

   braid_BaseVector  u, u_old;
   braid_Real        CWt;
   braid_Int         flo, fhi, fi, ci;
   braid_Int         nu, nrelax, interval;

   /* Required for Richardson */
   MPI_Comm            comm       = _braid_CoreElt(core, comm);
   braid_Int           ilower     = _braid_GridElt(grids[level], ilower);
   braid_Int           iupper     = _braid_GridElt(grids[level], iupper);
   braid_Int           cupper     = _braid_GridElt(grids[level], cupper);
   braid_Int           clower     = _braid_GridElt(grids[level], clower);
   braid_Int           cfactor    = _braid_GridElt(grids[level], cfactor);
   braid_Int           gupper     = _braid_GridElt(grids[level], gupper);
   braid_Real         *ta         = _braid_GridElt(grids[level], ta);
   braid_Real         *dtk        = _braid_CoreElt(core, dtk);
   braid_Int           order      = _braid_CoreElt(core, order);
   braid_BaseVector   *ua         = _braid_GridElt(grids[level], ua);
   braid_Int           richardson = _braid_CoreElt(core, richardson);
   braid_Int           iter       = _braid_CoreElt(core, niter);
   braid_Int           nrefine    = _braid_CoreElt(core, nrefine);
   braid_BufferStatus  bstatus    = (braid_BufferStatus)core;
   braid_StepStatus    status     = (braid_StepStatus)core;
   braid_Real          tol        = _braid_CoreElt(core, tol );
   
   /* Required for Richardson */
   braid_BaseVector bigstep;
   braid_Real *ta_c, time_left, a, b;
   braid_Int  size, proc;
   MPI_Request send_request, recv_request;
   void *recv_buff, *send_buff;
   braid_Int send_flag;
   send_flag = 0;


   /* In this case, nothing needs to be done regarding Richardson */
   if ( level > 0 || ncpoints == 0 || nlevels <= 1 || iter + nrefine == 0 )
   {
      richardson = 0;
   }
   
   nrelax  = nrels[level];
   CWt     = CWts[level];

   for (nu = 0; nu < nrelax; nu++)
   {
      _braid_UCommInit(core, level);
   
      /* Required for Richardson */
      braid_Int dtk_index = ncpoints;
      
      /* Richardson option requires access to the C-point in the left-most interval, so send it */
      if ( richardson )
      {
         if ( ilower > _braid_CoreElt(core, initiali))
         {
            //Need to post a recv for a C-point from the left. 
            _braid_GetProc(core, level, clower-cfactor, &proc);
            _braid_BufferStatusInit( 0, 0, bstatus );
            _braid_BaseBufSize(core, app,  &size, bstatus);
             recv_buff = malloc(size);
             
             MPI_Irecv(recv_buff, size, MPI_BYTE, proc, 84, comm, &recv_request);
         }
         if( ((iupper < gupper) && (cupper + cfactor <= gupper)) || _braid_CoreElt(core, periodic) )            
         {
            //Need to post a send of ciupper         
            _braid_GetProc(core, level, cupper+cfactor, &proc); 
            _braid_BufferStatusInit( 0 ,0 ,bstatus );
            _braid_BaseBufSize(core, app,  &size, bstatus);
            send_buff = malloc(size);
           
            braid_Int iu, is_stored;
            _braid_UGetIndex(core, level, cupper, &iu, &is_stored);
            _braid_BaseBufPack(core, app,  ua[iu], send_buff, bstatus);

            size = _braid_StatusElt( bstatus, size_buffer );
            MPI_Isend(send_buff, size, MPI_BYTE, proc, 84, comm, &send_request);
            send_flag = 1;
         }
      }      

      /* Start from the right-most interval */
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

         /* For Richardson, now receive the C-point in left-most interval */
         if ( richardson && (ci > _braid_CoreElt(core, initiali)) )  
         {
            if ( ci == clower )
            {
              /* The needed C-point is coming as a message. Wait and unpack */
               MPI_Wait( &recv_request, MPI_STATUS_IGNORE);
               _braid_BufferStatusInit(0,0,bstatus);
               _braid_BaseBufUnpack(core, app, recv_buff, &bigstep, bstatus);
               _braid_TFree(recv_buff);                  
              
               ta_c = _braid_GridElt(grids[1], ta );
               time_left = ta_c[-1];
            }
            else
            {
               /* We own the solution at previous coarse grid point, so clone it */
                braid_Int iu, is_stored;
               _braid_UGetIndex(core, level, ci - cfactor , &iu, &is_stored);

               _braid_BaseClone(core, app, ua[iu], &bigstep );
               time_left = ta[ci-cfactor-ilower];
            }
         }

         /* F-relaxation */
         for (fi = flo; fi <= fhi; fi++)
         {
            _braid_Step(core, level, fi, NULL, u);
            _braid_USetVector(core, level, fi, u, 0);
         }

         /* C-relaxation */
         if (ci > _braid_CoreElt(core, initiali))
         {
            /* If weighted Jacobi, store the previous u-value,
             *   Note, do no weighting if coarsest level*/
            if( (CWt != 1.0) && (level != (nlevels-1)) )
            {
               _braid_UGetVector(core, level, ci, &u_old);
            }

            _braid_Step(core, level, ci, NULL, u);

            if( (CWt != 1.0) && (level != (nlevels-1)) )
            {
               /* Apply weighted combination for w-Jacobi
                * u <--  omega*u_new + (1-omega)*u_old */
               _braid_BaseSum(core, app, (1.0 - CWt), u_old, CWt, u);
               _braid_BaseFree(core, app, u_old);
            }

            _braid_USetVector(core, level, ci, u, 1);

            /* Compute Richardson weights a and b */
            if ( richardson )
            {
               dtk_index--;
               braid_Real DTK = pow( ta[ci-ilower] - time_left, order );
               braid_Real dtk_temp = dtk[dtk_index];
               a = DTK / ( DTK - dtk_temp );
               b = - dtk_temp / ( DTK - dtk_temp );

               /* Note that we initialize StepStatus here in a non-standard
                * way, and hence cannot use _braid_Step(...). */
               _braid_StepStatusInit(time_left, ta[ci-ilower], ci-cfactor-1, tol,
                                     iter, level, nrefine, gupper, status);
               _braid_BaseStep(core, app, u, NULL, bigstep, level, status);

               _braid_BaseSum(core, app, a, u, b, bigstep );

                braid_Int iu, is_stored;
               _braid_UGetIndex(core, level, ci, &iu, &is_stored);

               _braid_BaseFree(core, app,  ua[iu]);

               _braid_BaseClone(core, app, bigstep, &ua[iu]);
               _braid_BaseFree(core, app,  bigstep);
            }

         }

         /* if ((flo <= fhi) && (interval == ncpoints)) */
         if ((flo <= fhi) && !(ci > _braid_CoreElt(core, initiali)))
         {
            _braid_BaseFree(core, app,  u);
         }
      }
      _braid_UCommWait(core, level);
   }

   /* If Richardson, then must wait for communication, and then free buffer */
   if (send_flag) 
   {
       MPI_Wait( &send_request, MPI_STATUS_IGNORE);
       _braid_TFree( send_buff );
   }

   return _braid_error_flag;
}

