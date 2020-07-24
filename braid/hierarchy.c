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
 * Initialize grid hierarchy
 *----------------------------------------------------------------------------*/

braid_Int
_braid_InitHierarchy(braid_Core    core,
                     _braid_Grid  *fine_grid,
                     braid_Int     refined)
{
   MPI_Comm       comm       = _braid_CoreElt(core, comm);
   braid_Int      max_levels = _braid_CoreElt(core, max_levels);
   braid_Int      min_coarse = _braid_CoreElt(core, min_coarse);
   braid_Int     *nrels      = _braid_CoreElt(core, nrels);
   braid_Int      nrdefault  = _braid_CoreElt(core, nrdefault);
   braid_Real    *CWts       = _braid_CoreElt(core, CWts);
   braid_Real     CWt_default= _braid_CoreElt(core, CWt_default);
   braid_Int      gupper     = _braid_CoreElt(core, gupper);
   braid_Int     *rfactors   = _braid_CoreElt(core, rfactors);
   braid_Real   **rdtvalues  = _braid_CoreElt(core, rdtvalues);
   braid_Int      nlevels    = _braid_CoreElt(core, nlevels);
   _braid_Grid  **grids      = _braid_CoreElt(core, grids);


   /* Required for Richardson */
   braid_Int      richardson = _braid_CoreElt(core, richardson);
   braid_Int      est_error  = _braid_CoreElt(core, est_error); 
   braid_Real    *dtk, *estimate;

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

   braid_Int         level;
   braid_Int         ilower, iupper;
   braid_Int         clower, cupper, cfactor, ncpoints, nupoints;
   braid_Real       *ta;
   braid_BaseVector *ua;
   braid_BaseVector *va;
   braid_BaseVector *fa;

   _braid_Grid      *grid;
   braid_Real       *f_ta;
   braid_Int         i, j, f_i, f_ilower, clo, chi, gclower, gcupper;
                    
   MPI_Request       request1, request2;
   MPI_Status        status;
   braid_Int         left_proc, right_proc;

   grids[0] = fine_grid;

   /* Do sequential time marching if min_coarse is already reached */
   if ( gupper <= min_coarse )
   {
      max_levels = 1;
   }

   /* Allocate space for rfactors (and initialize to zero) */
   ilower = _braid_GridElt(grids[0], ilower);
   iupper = _braid_GridElt(grids[0], iupper);
   rfactors = _braid_CTAlloc(braid_Int, iupper-ilower+2); /* Ensures non-NULL */
   for(i = 0; i < iupper-ilower+2; i++)
   {
      /* We need to ensure an rfactor of 1 for global index 0, and for
       * the case of a user-defined residual function and F-relaxation,
       * a default rfactor value at F-points */
      rfactors[i] = 1; 
   }
   _braid_CoreElt(core, rfactors) = rfactors;

   /* Allocate array of refiment dt values, initialized with NULL */
   rdtvalues = _braid_CTAlloc(braid_Real*, iupper-ilower+2); /* Ensures non-NULL */
   _braid_CoreElt(core, rdtvalues) = rdtvalues;

   /* Set up nrels and Cwt array */
   for (level = 0; level < max_levels; level++)
   {
      if (nrels[level] < 0)
      {
         nrels[level] = nrdefault;
      }
      if (CWts[level] < 0)
      {
         CWts[level] = CWt_default;
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
      
      gupper = gcupper;
      _braid_GridElt(grid, gupper) = gupper;

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
         va = _braid_CTAlloc(braid_BaseVector, iupper-ilower+2);
         fa = _braid_CTAlloc(braid_BaseVector, iupper-ilower+2);
         _braid_GridElt(grid, va_alloc) = va;
         _braid_GridElt(grid, fa_alloc) = fa;
         _braid_GridElt(grid, va)       = va+1;  /* shift */
         _braid_GridElt(grid, fa)       = fa+1;  /* shift */
      }

      // If on level that only stores C-points and not using the shell vector feature
      if ( ((_braid_CoreElt(core, storage) < 0) ||
            (level < _braid_CoreElt(core, storage))) &&
           (_braid_CoreElt(core, useshell)!=1) )
      {
         nupoints = _braid_GridElt(grid, ncpoints);   /* only C-points */
      }
      else
      {
         nupoints = iupper-ilower+1;                  /* all points */
      }

      ua = _braid_CTAlloc(braid_BaseVector, nupoints+1);
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
         if ( _braid_CoreElt(core, scoarsen) != NULL )
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
         if ( (left_proc > -1) && ( _braid_CoreElt(core, scoarsen) != NULL ) )
         {
            MPI_Send(&ta[0], sizeof(braid_Real), MPI_BYTE, left_proc, 1, comm);
         }

         /* Finish receive */
         if (left_proc > -1)
         {
            MPI_Wait(&request1, &status);
         }
         if ( (right_proc > -1) && ( _braid_CoreElt(core, scoarsen) != NULL ) )
         {
            MPI_Wait(&request2, &status);
         }
      }
   }

   /* Required for Richardson.  Allocate the dtk and estimate arrays */ 
   if ( richardson || est_error )
   {
      grid = grids[0];
      ilower   = _braid_GridElt(grid, ilower);
      iupper   = _braid_GridElt(grid, iupper);
      ncpoints = _braid_GridElt(grid, ncpoints);

      dtk = _braid_CTAlloc(braid_Real, ncpoints);
      _braid_CoreElt(core, dtk) = dtk;
      if (est_error)
      {
         estimate = _braid_CTAlloc(braid_Real, iupper-ilower+1);
         for(j = 0; j < (iupper-ilower+1); j++)
         {
            estimate[j] = -1;
         }
         _braid_CoreElt(core, estimate) = estimate;
      }   
      
      _braid_GetDtk(core);
   }


   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Returns the coarsening factor to use on grid 'level'
 *----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------
 * Set initial guess at C-points, and initialize shell at F-points when using
 * shell vectors.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_InitGuess(braid_Core  core,
                 braid_Int   level)
{
   braid_App          app      = _braid_CoreElt(core, app);
   braid_Int          seq_soln = _braid_CoreElt(core, seq_soln);
   braid_Int          nrefine  = _braid_CoreElt(core,nrefine);
   _braid_Grid      **grids    = _braid_CoreElt(core, grids);
   braid_Int          ilower   = _braid_GridElt(grids[level], ilower);
   braid_Int          iupper   = _braid_GridElt(grids[level], iupper);
   braid_Int          clower   = _braid_GridElt(grids[level], clower);
   braid_Int          cupper   = _braid_GridElt(grids[level], cupper);
   braid_Int          cfactor  = _braid_GridElt(grids[level], cfactor);
   braid_Real        *ta       = _braid_GridElt(grids[level], ta);
   braid_BaseVector  *va       = _braid_GridElt(grids[level], va);

   braid_BaseVector  u;
   braid_Int         i, iu, sflag;

   if ( (level == 0) && (seq_soln == 1) )
   {
      /* If first processor, grab initial condition */
      if(ilower == 0)
      {
         /* If we have already refined, then an initial init has already been done */
         if(nrefine > 0)
         {
            _braid_UGetVector(core, 0, 0, &u);    /* Get stored vector */
         }
         else
         {
            _braid_BaseInit(core, app,  ta[0], &u);
            _braid_USetVector(core, 0, 0, u, 0);
         }
         ilower += 1;
      }
      /* Else, receive point to the left */
      else
      {
         _braid_UCommInitBasic(core, 0, 1, 0, 0);     /* Post receive to the left*/
         _braid_UGetVector(core, 0, ilower-1, &u);    /* Wait on receive */
      }

      /* Set Flag so that USetVector initiates send when iupper is available */
      _braid_GridElt(grids[level], send_index)  = iupper;

      /* Initialize all points on the finest grid with sequential time marching */
      for(i = ilower; i <= iupper; i++)
      {
         _braid_Step(core, 0, i, NULL, u);       /* Step forward */
         _braid_USetVector(core, 0, i, u, 0);    /* Store: copy u into core,
                                                    sending to left if needed */
      }

      _braid_UCommWait(core, 0);                 /* Wait on comm to finish */
   }
   else if (level == 0)
   {
      if (_braid_CoreElt(core, useshell)==1)
      {
         for (i = ilower; i <=iupper; i++)
         {
            if (_braid_IsCPoint(i,cfactor))
            {
               // We are on a C-point, init full vector
               _braid_BaseInit(core, app,  ta[i-ilower], &u);
            }
            else
            {
               // We are on a F-point, init shell only
               _braid_BaseSInit(core,  app, ta[i-ilower], &u);
            }
            _braid_USetVectorRef(core, level, i, u);
         }
      }
      else
      {
         /* Only initialize the C-points on the finest grid */
         for (i = clower; i <= cupper; i += cfactor)
         {
            _braid_BaseInit(core, app,  ta[i-ilower], &u);
            _braid_USetVectorRef(core, level, i, u);
         }
      }
   }
   else
   {
      for (i = ilower; i <= iupper; i++)
      {
         _braid_UGetIndex(core, level, i, &iu, &sflag);
         if (sflag == 0) // Full point
         {
            _braid_BaseClone(core, app,  va[i-ilower], &u);
            _braid_USetVectorRef(core, level, i, u);
         }
         else if (sflag == -1) // Shell
         {
            _braid_BaseSClone(core,  app, va[i-ilower], &u);
            _braid_USetVectorRef(core, level, i, u);
         }
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Copy the initialized C-points on the fine grid, to all coarse levels.
 * Allows first down cycle to be skipped, in FMG fashion.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_CopyFineToCoarse(braid_Core  core)
{
   braid_App      app     = _braid_CoreElt(core, app);
   _braid_Grid  **grids   = _braid_CoreElt(core, grids);
   braid_Int      nlevels = _braid_CoreElt(core, nlevels);
   
   braid_Int      f_index, index, iu, is_stored, level, f_cfactor;
   braid_Int      ilower, iupper;
   braid_BaseVector   u, *va;

   for(level = 1; level < nlevels; level++)
   {

      f_cfactor = _braid_GridElt(grids[level-1], cfactor);
      iupper    = _braid_GridElt(grids[level], iupper);
      ilower    = _braid_GridElt(grids[level], ilower);
      va        = _braid_GridElt(grids[level], va);

      /* Loop over all points belonging to this processor, and if a C-point,
       * then carry out spatial coarsening and copy to ua and va */
      for (index=ilower; index<=iupper; index++)
      {
         _braid_MapCoarseToFine(index, f_cfactor, f_index);
         _braid_UGetVector(core, level-1, f_index, &u);
         _braid_Coarsen(core, level, f_index, index, u, &va[index-ilower]);
         
         _braid_BaseFree(core, app,  u);
         _braid_BaseClone(core, app,  va[index-ilower], &u);
         _braid_USetVectorRef(core, level, index, u);
         _braid_UGetIndex(core, level, index, &iu, &is_stored);
         if (is_stored == -2) /* Case where F-points are not stored, and we are not using shell vectors */
         {
            _braid_BaseFree(core, app,  u);
         }
         else if (is_stored == -1) /* This is a shell vector */
         {
            // We free the data in u, keeping the shell
            _braid_BaseSFree(core,  app, u);
         }
 
      }
   }

   return _braid_error_flag;
}



/*----------------------------------------------------------------------------
 * Propagate time step information required to compute the Richardson error
 * estimate at each C-point. This can be done at any time, but does require
 * some communication. This fills in error_factors at the C-points. 
 *----------------------------------------------------------------------------*/
braid_Int 
_braid_GetDtk(braid_Core core )
{
  
   _braid_Grid  **grids       = _braid_CoreElt(core, grids);
   _braid_Grid   *grid        = grids[0];
   MPI_Comm      comm         = _braid_CoreElt(core, comm);
   braid_Int     myid         = _braid_CoreElt(core, myid);
   braid_Int     order        = _braid_CoreElt(core, order); 
   braid_Int     ilower       = _braid_GridElt(grid, ilower);
   braid_Int     iupper       = _braid_GridElt(grid, iupper);
   braid_Int     gupper       = _braid_GridElt(grid, gupper);
   braid_Int     ncpoints     = _braid_GridElt(grid, ncpoints);
   braid_Int     cfactor      = _braid_GridElt(grid, cfactor);
   braid_Real   *dtk_core     = _braid_CoreElt(core, dtk ); 
   braid_Real   *ta           = _braid_GridElt(grid, ta);
   braid_Int fhi,flo,ci,interval,send_flag,recv_flag,i;
   braid_Real send_val, recv_val, dtk;
   MPI_Request recv_left, send_right ;
   
   int comm_size;
   MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
   
   if ( ilower <= iupper ) 
   {  
      send_flag = recv_flag = 0;
      send_val = recv_val = 0.0;
      
      if ( !_braid_IsCPoint( ilower - 1, cfactor) && ilower > 0 )
      {
         recv_flag = 1;
         MPI_Irecv( &recv_val, 1, braid_MPI_REAL, myid -1 , 56, comm, &recv_left );
      }

      braid_Int cpoint = ncpoints-1;
      for (interval = ncpoints; interval > -1; interval--)
      {
         _braid_GetInterval(core, 0, interval, &flo, &fhi, &ci);             
        
         dtk = 0;

         /* Wait for the value from the left */
         if ( interval == 0 && recv_flag )
         {
            MPI_Wait( &recv_left, MPI_STATUS_IGNORE);
            dtk = recv_val;

            recv_flag = 0;
         }
         
         //Finish up the interval 
         for ( i = flo; i <= fhi; i++ ) 
         {
            dtk += pow( ta[i-ilower] - ta[i-ilower-1] , order );
         }   
         /* Save dtk at the C points or send it on to the right is 
          * last point is not a C point */
         if ( ci > 0 ) 
         {
            dtk += pow( ta[ci-ilower] - ta[ci-ilower-1] , order );
            dtk_core[cpoint] = dtk;
            cpoint--;
         }
         
         if ( iupper != gupper && interval == ncpoints && !_braid_IsCPoint(iupper, cfactor) )
         {
            send_val  = dtk;
            send_flag = 1;
            MPI_Isend( &send_val, 1, braid_MPI_REAL, myid + 1, 56, comm, &send_right );
         }
      }   
      /* All recvieves should be complete. Finish up sends, and exit */
      if ( send_flag )
      {
         MPI_Wait( &send_right , MPI_STATUS_IGNORE);
      }
   }

   return _braid_error_flag;
}



