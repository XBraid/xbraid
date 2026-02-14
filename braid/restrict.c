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
 * F-Relax on level and restrict to level+1
 *
 * At the bottom of FRestrict(), the coarse-grid right-hand side is computed and
 * stored in c_fa[].  The result stored in c_fa[] must be considered carefully
 * to maintain a fixed point nature of the algorithm.  The FAS coarse-grid
 * equation is
 *
 *   A_c(u_c) = R(f - A(u)) + A_c(R(u))
 *
 * where R() is restriction, and A() and A_c() are the fine and coarse
 * operators.  We desire our algorithm to be a fixed-point method, in that the
 * exact solution (from sequential time stepping) should yield a zero initial
 * residual and that the residual should stay zero as Braid iterates.  This is
 * actually a subtle point.  In this case, the solution u_c should be equivalent
 * to the solution u.  Since R(f - A(u)) = 0, this implies that
 *
 *   A_c(u_c) = A_c(R(u))
 *
 * This further implies that the initial guess provided to the implicit solve
 * routines in A_c() must be consistent on the right and left hand sides.
 *
 * The function (@ref _braid_GetUInit()) provides this functionality by
 * encapsulating all decisions about how to set this initial guess (called
 * ustop) The ustop values for A_c(R(u)) are set in (@ _braid_Residual) and for
 * A_c(u_c), they are set in (@ _braid_step).
 *
 * storage == -2 implies that (@ref _braid_GetUInit()) will always use the
 * previous time step value as the initial guess, even if you have better
 * information, which is the case at C-points on the fine grid, and at all
 * points on coarse-levels.
 *
 * storage == -1 implies that (@ref _braid_GetUInit()) will always use the best
 * information.  The value at the previous time step is used as the initial
 * guess only at F-points on the fine-grid.  The va[] values provide all initial
 * guess at F- and C-points on coarse-grids.
 *
 * storage == 0 is the same as -1, except that storage also exists at fine-grid
 * F-points.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FRestrict(braid_Core   core,
                 braid_Int    level)
{
   MPI_Comm              comm         = _braid_CoreElt(core, comm);
   braid_App             app          = _braid_CoreElt(core, app);
   _braid_Grid         **grids        = _braid_CoreElt(core, grids);
   braid_AccessStatus    astatus      = (braid_AccessStatus)core;
   braid_ObjectiveStatus ostatus      = (braid_ObjectiveStatus)core;
   braid_Int             iter         = _braid_CoreElt(core, niter);
   braid_Int             print_level  = _braid_CoreElt(core, print_level);
   braid_Int             access_level = _braid_CoreElt(core, access_level);
   braid_Int             resid_compute= _braid_CoreElt(core, resid_compute);
   braid_Int             tnorm        = _braid_CoreElt(core, tnorm);
   braid_Real           *tnorm_a      = _braid_CoreElt(core, tnorm_a);
   braid_Int             nrefine      = _braid_CoreElt(core, nrefine);
   braid_Int             gupper       = _braid_CoreElt(core, gupper);
   braid_Int             cfactor      = _braid_GridElt(grids[level], cfactor);
   braid_Int             ncpoints     = _braid_GridElt(grids[level], ncpoints);
   braid_Real           *ta           = _braid_GridElt(grids[level], ta);
   braid_Int             f_ilower     = _braid_GridElt(grids[level], ilower);
   _braid_CommHandle    *recv_handle  = NULL;
   _braid_CommHandle    *send_handle  = NULL;

   /* Required for Richardson */
   braid_Int            richardson  = _braid_CoreElt(core, richardson);
   braid_Int            est_error   = _braid_CoreElt(core, est_error);
   braid_Int            order       = _braid_CoreElt(core, order);
   braid_Real          *ta_c        = _braid_GridElt(grids[1], ta );
   braid_Real          *dtk_core    = _braid_CoreElt(core, dtk);
   braid_Real          *estimate;    
   braid_Real           factor, dtk, DTK;

   /* Required for Delta correction */
   braid_Int    delta_correct = _braid_DoDeltaCorrect(core, level, iter);
   braid_Vector delta_action;

   braid_Int            c_level, c_ilower, c_iupper, c_index, c_i, c_ii;
   braid_BaseVector     c_u, *c_va, *c_fa;

   braid_BaseVector     u, r;
   braid_Int            interval, flo, fhi, fi, ci;
   braid_Real           rnorm, grnorm, rnorm_temp, rnm;

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
      else if (ci > _braid_CoreElt(core, initiali))
      {
         _braid_UGetVector(core, level, ci-1, &r);
      }

      /* F-relaxation */
      _braid_GetRNorm(core, -1, &rnm);
      for (fi = flo; fi <= fhi; fi++)
      {
         _braid_Step(core, level, fi, braid_ASCaller_FRestrict, NULL, r);
         _braid_USetVector(core, level, fi, r, 0);
         
         /* Allow user to process current vector, note that r here is
          * temporarily holding the state vector */
         if( (access_level >= 3) )
         {
            _braid_AccessStatusInit(ta[fi-f_ilower], fi, rnm, iter, level, nrefine, gupper,
                                    0, 0, braid_ASCaller_FRestrict, r->basis, astatus);
            _braid_AccessVector(core, astatus, r);
         }

         /* Evaluate the user's local objective function at F-points on finest grid */
         if ( _braid_CoreElt(core, adjoint) && level == 0)
         {
            _braid_ObjectiveStatusInit(ta[fi-f_ilower], fi, iter, level, nrefine, gupper, ostatus);
            _braid_AddToObjective(core, r, ostatus);
         }
      }

      /* Allow user to process current C-point */
      if( (access_level>= 3) && (ci > -1) )
      {
         _braid_UGetVectorRef(core, level, ci, &u);
         _braid_AccessStatusInit(ta[ci-f_ilower], ci, rnm, iter, level, nrefine, gupper,
                                 0, 0, braid_ASCaller_FRestrict, u->basis, astatus);
         _braid_AccessVector(core, astatus, u);
      }

      /* Evaluate the user's local objective function at CPoints on finest grid */
      if (_braid_CoreElt(core, adjoint) && level == 0 && (ci > -1) )
      {
         _braid_ObjectiveStatusInit(ta[ci-f_ilower], ci, iter, level, nrefine, gupper, ostatus);
         _braid_UGetVectorRef(core, 0, ci, &u);
         _braid_AddToObjective(core, u, ostatus);
      }
         
      
      /* Compute residual and restrict */
      if (ci > _braid_CoreElt(core, initiali))
      {
         /* Compute FAS residual */
         _braid_UGetVectorRef(core, level, ci, &u);
         _braid_FASResidual(core, level, ci, u, r);

         /* Compute rnorm (only on level 0). Richardson computes the rnorm later */
         if (level == 0 && !richardson )
         {
            _braid_BaseSpatialNorm(core, app,  r, &rnorm_temp);
            tnorm_a[interval] = rnorm_temp;       /* inf-norm uses tnorm_a */
            if(tnorm == 1) 
            {  
               rnorm += rnorm_temp;               /* one-norm combination */ 
            }
            else if(tnorm == 2)
            {  
               rnorm += (rnorm_temp*rnorm_temp);  /* two-norm combination */
            }
         }

         /* Restrict u and residual, coarsening in space if needed */
         _braid_MapFineToCoarse(ci, cfactor, c_index);
         _braid_Coarsen(core, c_level, ci, c_index, u, &c_va[c_index-c_ilower]);
         _braid_Coarsen(core, c_level, ci, c_index, r, &c_fa[c_index-c_ilower]);
         if (braid_LINEAR)
         {
            /* Set initial guess to zero on coarse grids */
            _braid_BaseSum(core, app, -1.0, c_va[c_index-c_ilower], 1.0, c_va[c_index-c_ilower]);
         }
      }
      else if (ci == 0)
      {
         /* Restrict initial condition, coarsening in space if needed */
         _braid_UGetVectorRef(core, level, 0, &u);
         _braid_Coarsen(core, c_level, 0, 0, u, &c_va[0]);
         if (braid_LINEAR)
         {
            /* Set initial guess to zero on coarse grids */
            _braid_BaseSum(core, app, -1.0, c_va[0], 1.0, c_va[0]);
         }
      }

      if ((flo <= fhi) || (ci > _braid_CoreElt(core, initiali)))
      {
         _braid_BaseFree(core, app,  r);
      }
   }
   _braid_UCommWait(core, level);
  
   /* Now apply coarse residual to update fa values */

   /* Set initial guess on coarse level */
   _braid_InitGuess(core, c_level);

   if (braid_LINEAR)
   {
      /* Don't need to update the rhs in the linear case (no Richardson either) */
   }
   else
   {
   /* Initialize update of c_va[-1] boundary */
   if (c_ilower <= c_iupper)
   {
      _braid_CommRecvInit(core, c_level, c_ilower-1, &c_va[-1], &recv_handle);
      _braid_CommSendInit(core, c_level, c_iupper, c_va[c_iupper-c_ilower], &send_handle);
   }

   /* Allocate temporary error estimate array */
   if ( level == 0 && est_error )
   {
        estimate = _braid_CTAlloc(braid_Real, c_iupper-c_ilower + 1);
   }

   /* Start with rightmost point */
   for (c_i = c_iupper; c_i >= c_ilower; c_i--)
   {
      if (c_i > _braid_CoreElt(core, initiali))
      {
         c_ii = c_i - c_ilower;
         if (c_ii == 0)
         {
            /* Finalize update of c_va[-1] */
            _braid_CommWait(core, &recv_handle);
         }

         if ( delta_correct )
         {
            /* need an extra copy of c_va[c_ii-1] for tau-correction later */
            _braid_CoreFcn(core, clone)(app, c_va[c_ii-1]->userVector, &delta_action);
         }

         _braid_BaseClone(core, app,  c_va[c_ii-1], &c_u);
         _braid_Residual(core, c_level, c_i, braid_ASCaller_Residual, c_va[c_ii], c_u);
         
         /* Richardson computes norm here, and recombines solution at C-points for higher accuracy */
         if ( level == 0 && richardson ) 
         {    
               dtk = dtk_core[c_ii];
               DTK = pow( ta_c[c_ii] - ta_c[c_ii-1], order );
               /* Factor computes \bar{a} from Richardson paper, used to scale RHS term in FAS */
               factor = DTK / ( DTK - dtk );
               _braid_BaseSum(core, app, factor, c_u, factor, c_fa[c_ii]);

               /* Compute the rnorm */
               _braid_BaseSum(core, app, 1.0, c_fa[c_ii], -1.0, c_u );
               _braid_BaseSpatialNorm(core, app, c_u, &rnorm_temp);
               
               tnorm_a[c_ii] = rnorm_temp;       /* inf-norm uses tnorm_a */
               if(tnorm == 1) 
               {  
                   rnorm += rnorm_temp;               /* one-norm combination */ 
               }
               else if(tnorm == 2)
               {  
                  rnorm += (rnorm_temp*rnorm_temp);  /* two-norm combination */
               }
         }    
         else
         {
            if ( delta_correct )
            {
               /* tau correction */
               _braid_CoreFcn(core, sum)(app, 1.0, c_u->userVector, 1.0, c_fa[c_ii]->userVector);

               /* get Delta correction */
               _braid_BaseSumBasis(core, app, 1.0, c_u->basis, 1.0, c_fa[c_ii]->basis);

               /* get the action of Delta on u_{i-1} */
               _braid_LRDeltaDot(core, app, delta_action, c_fa[c_ii]->basis, c_va[c_ii-1]->basis);
               _braid_CoreFcn(core, sum)(app, -1., delta_action, 1., c_fa[c_ii]->userVector);

               _braid_CoreFcn(core, free)(app, delta_action);
            }
            else
            {
               _braid_BaseSum(core, app,  1.0, c_u, 1.0, c_fa[c_ii]);
            }
         }

         /* Compute Richardson error estimator */
         if ( level == 0 && est_error )
         {
              braid_Real est_temp;
              dtk = dtk_core[c_ii]  ;
              DTK = pow( ta_c[c_ii] - ta_c[c_ii-1], order );
             _braid_BaseSpatialNorm(core, app, c_fa[c_ii], &est_temp);
             if ( richardson )
             {
                est_temp = est_temp / DTK;
             }
             else
             {
                /* Note, the algebra:
                 *     (DTK - dtk) = (m*dt)^k - m (dt)^k = m (dt)^k (m^k_g - 1)
                 * where k is the local order, k_g = (k-1) is the global order,
                 * dt is the fine grid step size.  The m (dt)^k is canceled out
                 * in FinalizeErrorEstimates.
                 */
                est_temp = est_temp / ( DTK - dtk );
             }
             estimate[ c_ii ] = est_temp; 
         }

         _braid_BaseFree(core, app, c_u);
      }
   }
   _braid_CommWait(core, &send_handle);
   }
   
   /* Compute global rnorm (only on level 0, and only if resid_compute is turned on) */
   if (level == 0)
   {
      if(resid_compute >= 1)
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
      }
      else
      {
         grnorm = -1.0;
      }

      /* Store new rnorm */
      _braid_SetRNorm(core, -1, grnorm);
   }
   
   /* If debug printing, print out tnorm_a for this interval. This
    * should show the serial propagation of the exact solution */
   if ((print_level > 2) && (level == 0) )
   {
      _braid_PrintSpatialNorms(core, tnorm_a, ncpoints);
   }

   if (braid_LINEAR)
   {
      /* Don't need to update the rhs in the linear case (no Richardson either) */
   }
   else
   {
   /* Need to finalize the error estimates at the F-points */ 
   if ( level == 0 && est_error )
   {
      _braid_FinalizeErrorEstimates( core, estimate , c_iupper-c_ilower + 1 );
      _braid_TFree(estimate);
   } 
   }
   
   return _braid_error_flag;
}


/*----------------------------------------------------------------------------
 * Finalize Richardson error estimates at F-points.
 *----------------------------------------------------------------------------*/
braid_Int
_braid_FinalizeErrorEstimates( braid_Core   core, 
                               braid_Real  *estimate,
                               braid_Int    length)
{

   braid_Int     myid        = _braid_CoreElt(core, myid );
   MPI_Comm      comm        = _braid_CoreElt(core, comm );
  _braid_Grid  **grids       = _braid_CoreElt(core, grids);
   braid_Real   *error_est   = _braid_CoreElt(core, estimate);
   braid_Real   *dtk_core    = _braid_CoreElt(core, dtk);
   braid_Int     ilower   = _braid_GridElt(grids[0],ilower);
   braid_Int     iupper   = _braid_GridElt(grids[0],iupper);
   braid_Int     gupper   = _braid_GridElt(grids[0],gupper);
   braid_Int     cfactor  = _braid_GridElt(grids[0],cfactor);
   braid_Int     ncpoints = _braid_GridElt(grids[0], ncpoints);

   braid_Int recv_flag, send_flag_l, send_flag_r, i, last_cpoint, estimate_index ; 
   braid_Real recv_value, send_value_l, send_value_r, factor ;
   MPI_Request recv_request, send_request_l, send_request_r;

   if ( ilower <= iupper )
   {

      /* Recv situation 1: Recv from left if ilower is > last C point */
      send_flag_l = send_flag_r = recv_flag = 0; 
      last_cpoint = ( gupper / cfactor ) * cfactor ;
      if ( ilower > last_cpoint && ilower > 0 )
      { 
         recv_flag++;
         MPI_Irecv( &recv_value, 1, braid_MPI_REAL, myid -1, 98, comm, &recv_request );
      }
      
      /* Recv situation 2: recv from right if iupper is a f point */
      else if ( !_braid_IsCPoint(iupper, cfactor) && iupper != gupper )
      {
         recv_flag++;
         MPI_Irecv( &recv_value, 1, braid_MPI_REAL, myid + 1, 98, comm, &recv_request );
      }
      
      /* Send situation 1: Send to the right if proc owns a point in last interval, but not gupper */
      if ( iupper >= last_cpoint && iupper != gupper )
      {
         /* If we have no c points, wait to recv before sending */ 
         if ( ncpoints == 0 )
         {  
            recv_flag--;
            MPI_Wait( &recv_request , MPI_STATUS_IGNORE );
            send_value_r = recv_value;
         }
         else
         {
            send_value_r = estimate[length - 1] * dtk_core[ length - 1 ];
         }

         send_flag_r++;
         MPI_Isend( &send_value_r, 1, braid_MPI_REAL, myid + 1 , 98, comm, &send_request_r );
      }
      
      /* Send situation 2: Send to the left if ilower -1 is not a C point */
      if ( ilower <= last_cpoint && ilower > 0 && !_braid_IsCPoint( ilower - 1, cfactor ) )
      {
         /* If we have no c points, wait to recv before sending */
         if ( ncpoints == 0 )
         {
            recv_flag--;
            MPI_Wait( &recv_request, MPI_STATUS_IGNORE );
            send_value_l = recv_value; 
         }
         else
         {
            send_value_l = estimate[0] * dtk_core[0];
         }

         send_flag_l++;
         MPI_Isend( &send_value_l, 1, braid_MPI_REAL, myid - 1, 98, comm, &send_request_l );
      }

      /* Update the values */
      estimate_index = 0;
      factor = 0;   
      for ( i = ilower; i <= iupper; i++ )
      {
         /* If this is an F point interval at end of proc ownership */
         if (estimate_index >= length )
         {
            if (recv_flag)
            {
               recv_flag--;
               MPI_Wait( &recv_request , MPI_STATUS_IGNORE );
            }       
            factor = recv_value; 
         }
         else
         {
            factor = estimate[ estimate_index ] * dtk_core[ estimate_index ] ;
         }

         /* Set Error Estimate 
          * Note that we have scaled factor above by dtk_core, which equals
          * m*dt^k, which is the coarsening factor m, fine-grid step size dt,
          * and local order k.  This is needed to cancel out an extra factor of
          * (1/m*dt^k) in estimate.  The error estimate then equals, 
          * || U^{fine}_i - U^{coarse}_i || / (m^{k_g} - 1
          * where k_g is the global order of the time stepping method.
          * */

         error_est[ i-ilower ] = factor ;
         
         /* increase the index at the Cpoints, unless this is the
          * last C point. In that case, use the same index */ 
         if ( _braid_IsCPoint( i, cfactor ) && i != last_cpoint )
         {
            estimate_index++;
         }
      }

      /* All receives should be completed. Finish up sends */
      if ( send_flag_r )
      {
         MPI_Wait( & send_request_r, MPI_STATUS_IGNORE );
      }
      if ( send_flag_l )
      {
         MPI_Wait( & send_request_l, MPI_STATUS_IGNORE );
      }

   }
   
   return _braid_error_flag;
}


