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

/**
 *  Source file implementing the wrapper for the user routines
 **/

#include "_braid.h"
#include "util.h"

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseVectorCreate(braid_BaseVector* u_ptr)
{
   braid_BaseVector u;

   u = _braid_TAlloc(_braid_BaseVector, 1);
   u->userVector = NULL;
   u->bar        = NULL;
   u->basis      = NULL;

   *u_ptr = u;
   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int 
_braid_BaseStep(braid_Core       core,
                braid_App        app,
                braid_BaseVector ustop,
                braid_BaseVector fstop,
                braid_BaseVector u,
                braid_Int        level,
                braid_StepStatus status )
{
   _braid_Action   *action;
   braid_Vector     u_copy, ustop_copy;
   braid_VectorBar  bar_copy, ustopbar_copy;
   braid_Int        myid        = _braid_CoreElt(core, myid);
   braid_Int        verbose_adj = _braid_CoreElt(core, verbose_adj);
   braid_Int        record      = _braid_CoreElt(core, record);
   braid_Real       t           = _braid_CoreElt(core, t);
   braid_Real       tnext       = _braid_CoreElt(core, tnext);
   braid_Int        tidx        = _braid_CoreElt(core, idx);
   braid_Int        iter        = _braid_CoreElt(core, niter);
   braid_Int        nrefine     = _braid_CoreElt(core, nrefine);
   braid_Int        gupper      = _braid_CoreElt(core, gupper);
   braid_Real       tol         = _braid_CoreElt(core, tol);
   braid_Real       timer       = 0.0;

   if (verbose_adj) _braid_printf("%d: STEP %.4f to %.4f, %d\n", myid, t, tnext, tidx);

   /* Record to the tape */
   if ( record )
   {
      /* Set up the action and push it to the actiontape */
      action             = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall  = STEP;
      action->core       = core;
      action->inTime     = t;
      action->outTime    = tnext;
      action->inTimeIdx  = tidx;
      action->myid       = myid;
      action->braid_iter = iter;
      action->level      = level;
      action->nrefine    = nrefine;
      action->gupper     = gupper;
      action->tol        = tol;
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

      /* Copy & push u & ustop to primal tape */
      _braid_CoreFcn(core, clone)(app, u->userVector, &u_copy); 
      _braid_CoreFcn(core, clone)(app, ustop->userVector, &ustop_copy);  
      _braid_CoreElt(core, userVectorTape) = _braid_TapePush( _braid_CoreElt(core, userVectorTape), u_copy);
      _braid_CoreElt(core, userVectorTape) = _braid_TapePush( _braid_CoreElt(core, userVectorTape), ustop_copy);

      /* Copy & push ubar & ustopbar to bar tape */
      _braid_VectorBarCopy(u->bar, &bar_copy);
      _braid_VectorBarCopy(ustop->bar, &ustopbar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), bar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), ustopbar_copy);
  }

   /* Call the users Step function.  If periodic and integrating to the periodic
    * point, adjust tnext to be tstop. */
   if ( _braid_CoreElt(core, periodic) && (_braid_CoreElt(core, idx) < 0) )
   {
      _braid_CoreElt(core, tnext) = _braid_CoreElt(core, tstop);
   }
   timer = _braid_MPI_Wtime(core, 2);
   if ( fstop == NULL )
   {
      _braid_CoreFcn(core, step)(app, ustop->userVector, NULL, u->userVector, status);
   }
   else
   {
      /* fstop not supported by adjoint! */
      _braid_CoreFcn(core, step)(app, ustop->userVector, fstop->userVector, u->userVector, status);
   }
   _braid_CoreElt(core, timer_user_step) += _braid_MPI_Wtime(core, 2) - timer;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseInit(braid_Core        core,
                braid_App         app, 
                braid_Real        t,   
                braid_BaseVector *u_ptr )
{
   _braid_Action    *action;
   braid_BaseVector  u;
   braid_VectorBar   ubar;
   braid_Int         myid          = _braid_CoreElt(core, myid);
   braid_Int         verbose_adj   = _braid_CoreElt(core, verbose_adj);
   braid_Int         record        = _braid_CoreElt(core, record);
   braid_Int         adjoint       = _braid_CoreElt(core, adjoint);
   braid_Int         delta_correct = _braid_CoreElt(core, delta_correct);
   braid_Real        timer       = 0.0;
    
   if (verbose_adj) _braid_printf("%d: INIT\n", myid);

   /* Allocate the braid_BaseVector */
   _braid_BaseVectorCreate(&u);

   /* Allocate and initialize the userVector */
   timer = _braid_MPI_Wtime(core, 2);
   _braid_CoreFcn(core, init)(app, t, &(u->userVector));

   /* Allocate and initialize the basis vectors */
   if ( delta_correct )
   {
      _braid_BaseInitBasis(core, app, t, &(u->basis));
   }
   _braid_CoreElt(core, timer_user_init) += _braid_MPI_Wtime(core, 2) - timer;
   
   /* Allocate and initialize the bar vector */
   if ( adjoint ) 
   {
      ubar = _braid_TAlloc(_braid_VectorBar, 1);
      ubar->useCount = 1;
      _braid_CoreFcn(core, init)(app, t, &(ubar->userVector));
      _braid_CoreFcn(core, sum)(app, -1.0, ubar->userVector, 1.0, ubar->userVector);
      u->bar = ubar;
   }

   /* Record to the tape */
   if ( record )
   {
      /* Set up and push the action */
      action            = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall = INIT;
      action->core      = core;
      action->inTime    = t;
      action->myid      = myid;
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);
   }

   /* Set the return pointer */
   *u_ptr = u;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseInitBasis(braid_Core   core,
                     braid_App    app,
                     braid_Real   t,
                     braid_Basis *psi_ptr )
{
   braid_Basis u_basis;
   u_basis = _braid_TAlloc(_braid_Basis, 1);
   u_basis->rank = _braid_CoreElt(core, delta_rank);
   u_basis->userVecs = _braid_TAlloc(braid_Vector, u_basis->rank);

   for (braid_Int i = 0; i < u_basis->rank; i++)
   {
      /* Allocate and initialize userVector for each column */
      _braid_CoreFcn(core, init_basis)(app, t, i, &(u_basis->userVecs[i]));
   }
   /* orthonormalize the columns */
   _braid_GramSchmidt(core, app, u_basis, NULL);

   *psi_ptr = u_basis;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseClone(braid_Core         core,
                 braid_App          app,  
                 braid_BaseVector   u,    
                 braid_BaseVector  *v_ptr )
{
   _braid_Action    *action;
   braid_BaseVector  v;
   braid_VectorBar   ubar;
   braid_VectorBar   ubar_copy;
   braid_VectorBar   vbar_copy;
   braid_Int         myid         = _braid_CoreElt(core, myid);
   braid_Int         verbose_adj  = _braid_CoreElt(core, verbose_adj);
   braid_Int         record       = _braid_CoreElt(core, record);
   braid_Int         adjoint      = _braid_CoreElt(core, adjoint);
   braid_Real        timer       = 0.0;

   if (verbose_adj) _braid_printf("%d: CLONE\n", myid);

   /* Allocate the braid_BaseVector */
   _braid_BaseVectorCreate(&v);

   /* Allocate and copy the userVector */
   timer = _braid_MPI_Wtime(core, 2);
   _braid_CoreFcn(core, clone)(app, u->userVector, &(v->userVector) );
   _braid_CoreElt(core, timer_user_clone) += _braid_MPI_Wtime(core, 2) - timer;

   /* Allocate and clone the basis */
   if ( u->basis )
   {
      _braid_BaseCloneBasis(core, app, u->basis, &(v->basis));
   }

   /* Allocate and initialize the bar vector to zero */
   if ( adjoint )
   {
      ubar = _braid_TAlloc(_braid_VectorBar, 1);
      ubar->useCount = 1;
      _braid_CoreFcn(core, clone)(app, u->bar->userVector, &(ubar->userVector));
      _braid_CoreFcn(core, sum)(app, -1.0, ubar->userVector, 1.0, ubar->userVector);
      v->bar = ubar;
   } 

   /* Record to the tape */
   if ( record ) 
   {
      /* Set up and push the action */
      action            = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall = CLONE;
      action->core      = core;
      action->myid      = myid;
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

      /* Copy and push both bar vectors to the bartape */
      _braid_VectorBarCopy(u->bar, &ubar_copy);
      _braid_VectorBarCopy(v->bar, &vbar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), ubar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), vbar_copy);
   }

   *v_ptr = v;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseCloneBasis(braid_Core    core,
                      braid_App     app,  
                      braid_Basis   A,    
                      braid_Basis  *B_ptr)
{
   braid_Basis basis;
   basis = _braid_TAlloc(_braid_Basis, 1);
   basis->rank = A->rank;
   basis->userVecs = _braid_TAlloc(braid_Vector, basis->rank);

   braid_Real timer = _braid_MPI_Wtime(core, 2);
   for (braid_Int i = 0; i < basis->rank; i++)
   {
      /* Allocate and clone userVector for each column */
      _braid_CoreFcn(core, clone)(app, A->userVecs[i], &(basis->userVecs[i]));
   }
   _braid_CoreElt(core, timer_user_clone) += _braid_MPI_Wtime(core, 2) - timer;
   *B_ptr = basis;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseFree(braid_Core       core,
                braid_App        app,
                braid_BaseVector u )
{
   // if (u == NULL)
   // {
   //    return _braid_error_flag;
   // }

   _braid_Action *action;
   braid_Int      myid        = _braid_CoreElt(core, myid);
   braid_Int      verbose_adj = _braid_CoreElt(core, verbose_adj);
   braid_Int      adjoint     = _braid_CoreElt(core, adjoint);
   braid_Int      delta       = _braid_CoreElt(core, delta_correct);
   braid_Int      record      = _braid_CoreElt(core, record);
   braid_Real     timer       = 0.0;

   if (verbose_adj) _braid_printf("%d: FREE\n", myid);

   /* Record to the tape */
   if ( record )
   {
      /* Set up and push the action */
      action            = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall = FREE;
      action->core      = core;
      action->myid      = myid;
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);
   }
 
   /* Free the user's vector */
   timer = _braid_MPI_Wtime(core, 2);
   _braid_CoreFcn(core, free)(app, u->userVector);
   _braid_CoreElt(core, timer_user_free) += _braid_MPI_Wtime(core, 2) - timer;

   if ( u->basis && delta )
   {
      /* Free the basis/Lyapunov vectors */
      _braid_BaseFreeBasis(core, app, u->basis);
   }

   if ( adjoint )
   {
      /* Free the bar vector */
      _braid_VectorBarDelete(core, u->bar);
   }

   /* Free the braid_BaseVector */
   free(u);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseFreeBasis(braid_Core    core,
                     braid_App     app,
                     braid_Basis   b)
{
   braid_Real timer = _braid_MPI_Wtime(core, 2);
   /* Free the basis vectors */
   for (braid_Int i = 0; i < b->rank; i++)
   {
      _braid_CoreFcn(core, free)(app, b->userVecs[i]);
   }
   _braid_CoreElt(core, timer_user_free) += _braid_MPI_Wtime(core, 2) - timer;
   _braid_TFree(b->userVecs);
   free(b);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseSum(braid_Core        core,
               braid_App         app,    
               braid_Real        alpha,  
               braid_BaseVector  x,      
               braid_Real        beta,   
               braid_BaseVector  y )
{
   _braid_Action   *action;
   braid_VectorBar  xbar_copy;
   braid_VectorBar  ybar_copy;
   braid_Int        myid         =  _braid_CoreElt(core, myid);
   braid_Int        verbose_adj  =  _braid_CoreElt(core, verbose_adj);
   braid_Int        record       =  _braid_CoreElt(core, record);
   braid_Real       timer        = 0.0;

   if ( verbose_adj ) _braid_printf("%d: SUM\n", myid);

   /* Record to the tape */
   if ( record )
   {
      /* Set up and push the action */
      action             = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall  = SUM;
      action->core       = core;
      action->sum_alpha  = alpha;
      action->sum_beta   = beta;
      action->myid       = myid;
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

      /* Copy and push both bar vector to the bar tape */
      _braid_VectorBarCopy(x->bar, &xbar_copy);
      _braid_VectorBarCopy(y->bar, &ybar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), xbar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), ybar_copy);
   }

    /* Sum up the user's vector */
   timer = _braid_MPI_Wtime(core, 2);
   _braid_CoreFcn(core, sum)(app, alpha, x->userVector, beta, y->userVector);
   _braid_CoreElt(core, timer_user_sum) += _braid_MPI_Wtime(core, 2) - timer;

   /* Sum over the basis vectors */
   if ( x->basis && y->basis )
   {  /* Add together the bases if they both exist */
      _braid_BaseSumBasis(core, app, alpha, x->basis, beta, y->basis);
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseSumBasis(braid_Core  core,
                    braid_App   app,
                    braid_Real  alpha,  
                    braid_Basis A,
                    braid_Real  beta,
                    braid_Basis B )
{
   braid_Real timer = _braid_MPI_Wtime(core, 2);
   for (braid_Int i = 0; i < A->rank; i++)
   {
      _braid_CoreFcn(core, sum)(app, alpha, A->userVecs[i], beta, B->userVecs[i]);
   }
   _braid_CoreElt(core, timer_user_sum) += _braid_MPI_Wtime(core, 2) - timer;
   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseSpatialNorm(braid_Core        core,
                       braid_App         app,      
                       braid_BaseVector  u,    
                       braid_Real       *norm_ptr )
{
   braid_Real       timer       = 0.0;
   
   /* Compute the spatial norm of the user's vector */
   timer = _braid_MPI_Wtime(core, 2);
   _braid_CoreFcn(core, spatialnorm)(app, u->userVector, norm_ptr);
   _braid_CoreElt(core, timer_user_spatialnorm) += _braid_MPI_Wtime(core, 2) - timer;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseInnerProd(braid_Core        core,
                     braid_App         app,
                     braid_Vector      u,
                     braid_Vector      v,
                     braid_Real       *prod_ptr )
{
   braid_Real       timer       = 0.0;
   
   /* Compute the inner product between two user vectors */
   timer = _braid_MPI_Wtime(core, 2);
   _braid_CoreFcn(core, inner_prod)(app, u, v, prod_ptr);
   _braid_CoreElt(core, timer_user_inner_prod) += _braid_MPI_Wtime(core, 2) - timer;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseAccess(braid_Core          core,
                  braid_App           app,   
                  braid_BaseVector    u,     
                  braid_AccessStatus  status )
{
   _braid_Action   *action;
   braid_Real       t             = _braid_CoreElt(core, t);
   braid_Int        myid          = _braid_CoreElt(core, myid);
   braid_Int        verbose_adj   = _braid_CoreElt(core, verbose_adj);
   braid_Int        record        = _braid_CoreElt(core, record);
   braid_Real       timer         = 0.0;
   
   if ( verbose_adj ) _braid_printf("%d: ACCESS\n", myid);

   /* Record to the tape */
   if ( record )
   {
      /* Set up and push the action */
      action             = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall  = ACCESS;
      action->core       = core;
      action->inTime     = t;
      action->myid       = myid;
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);
   }

   /* Access the user's vector */
   timer = _braid_MPI_Wtime(core, 2);
   _braid_CoreFcn(core, access)(app, u->userVector, status);
   _braid_CoreElt(core, timer_user_access) += _braid_MPI_Wtime(core, 2) - timer;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseSync(braid_Core          core,
                braid_App           app,
                braid_SyncStatus    status)
{
   braid_Int        myid          = _braid_CoreElt(core, myid);
   braid_Int        verbose_adj   = _braid_CoreElt(core, verbose_adj);
   braid_Real       timer         = 0.0;
   if( verbose_adj ) _braid_printf("%d: SNYC\n", myid);

   /* Do adjoint stuff here */

   /* Call the user's sync function */
   timer = _braid_MPI_Wtime(core, 2);
   _braid_CoreFcn(core, sync)(app, status);
   _braid_CoreElt(core, timer_user_sync) += _braid_MPI_Wtime(core, 2) - timer;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseBufSize(braid_Core          core,
                   braid_App           app,      
                   braid_Int          *size_ptr, 
                   braid_BufferStatus  status ) 
{
   braid_Int  myid         = _braid_CoreElt(core, myid);
   braid_Int  verbose_adj  = _braid_CoreElt(core, verbose_adj);
   braid_Real timer        = 0.0;

   if ( verbose_adj ) _braid_printf("%d: BUFSIZE\n", myid);

   /* Call the users BufSize function */
   timer = _braid_MPI_Wtime(core, 2);
   _braid_CoreFcn(core, bufsize)(app, size_ptr, status);
   _braid_CoreElt(core, timer_user_bufsize) += _braid_MPI_Wtime(core, 2) - timer;

   _braid_StatusElt(status, size_buffer) = *size_ptr;

   /* extend buffer to fit basis vectors and add header information */
   if ( _braid_CoreElt(core, delta_correct) )
   {
      braid_Int rank = _braid_CoreElt(core, delta_rank);
      braid_Int size_basis = _braid_StatusElt(status, size_basis);
      if ( size_basis > 0)
      {
         *size_ptr += rank * size_basis;
      }
      else // not set, use default
      {
         _braid_StatusElt(status, size_basis) = *size_ptr;
         *size_ptr *= 1 + rank;
      }
      /* room for header (number of vectors, size uvec, size each basis vec) */
      *size_ptr += (2 + rank) * sizeof(braid_Int);
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseBufPack(braid_Core          core,
                   braid_App           app,
                   braid_BaseVector    u,
                   void               *buffer,
                   braid_BufferStatus  status )
{
   _braid_Action   *action;
   braid_VectorBar  ubar_copy;
   braid_Int        myid         = _braid_CoreElt(core, myid);
   braid_Int        verbose_adj  = _braid_CoreElt(core, verbose_adj);
   braid_Int        record       = _braid_CoreElt(core, record);
   braid_Int        sender       = _braid_CoreElt(core, send_recv_rank);

   braid_Real       timer        = 0.0;

   if ( verbose_adj ) _braid_printf("%d: BUFPACK\n",  myid );

   /* Record to the tape */
   if ( record )
   {
      /* Set up and push the action */
      action                 = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall      = BUFPACK;
      action->core           = core;
      action->send_recv_rank = sender;
      action->messagetype    = _braid_StatusElt(status, messagetype);
      action->size_buffer    = _braid_StatusElt(status, size_buffer);
      action->inTimeIdx      = _braid_StatusElt(status, idx);
      action->level          = _braid_StatusElt(status, level);
      action->myid           = myid;
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

      /* Copy and push the bar pointer to the bar tape */
      _braid_VectorBarCopy(u->bar, &ubar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), ubar_copy);
   }

   if ( _braid_CoreElt(core, delta_correct) )
   {
      /*
       * Delta correction requires sending a full set of basis vectors alongside every state vector.
       * Using a header prepended to the buffer, we can efficiently pack all of these vectors contiguously
       * in one message.
      */
      braid_Int *header, head_size, tot_size;
      braid_Byte* data_buf;

      /*
       * header pointer is offset by one, so
       * header[-1] stores number of vectors,
       * header[0] stores size of user vector 
       * header[i] stores size of ith basis vector 
       */
      header = ((braid_Int*)buffer) + 1;
      header[-1] = u->basis->rank + 1;                   /* total number of vectors */
      head_size = (header[-1] + 1) * sizeof(braid_Int);  /* header[-1] + size of each vector */

      data_buf = ((braid_Byte*)buffer) + head_size;      /* buffer must have type for legal pointer arithmetic */

      /* BufPack the user's vector */
      timer = _braid_MPI_Wtime(core, 2);
      _braid_CoreFcn(core, bufpack)(app, u->userVector, (void*)data_buf, status);

      /* get size actually written and write it to the header */
      header[0] = _braid_StatusElt(status, size_buffer);
      data_buf += header[0]; /* increment the data pointer */
      tot_size = head_size + header[0];

      for (braid_Int i = 0; i < u->basis->rank; i++)
      {
         _braid_StatusElt(status, size_buffer) = _braid_StatusElt(status, size_basis);
         _braid_CoreFcn(core, bufpack)(app, u->basis->userVecs[i], (void*)data_buf, status);
         header[i+1] = _braid_StatusElt(status, size_buffer);
         data_buf += header[i+1];
         tot_size += header[i+1];
      }
      _braid_CoreElt(core, timer_user_bufpack) += _braid_MPI_Wtime(core, 2) - timer;

      /* set actual total size of the buffer for MPI */
      _braid_StatusElt(status, size_buffer) = tot_size;

      return _braid_error_flag;
   }
   /* else default behavior */
   
   /* BufPack the user's vector */
   timer = _braid_MPI_Wtime(core, 2);
   _braid_CoreFcn(core, bufpack)(app, u->userVector, buffer, status);
   _braid_CoreElt(core, timer_user_bufpack) += _braid_MPI_Wtime(core, 2) - timer;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseBufUnpack(braid_Core          core,
                     braid_App           app,
                     void               *buffer,
                     braid_BaseVector   *u_ptr,
                     braid_BufferStatus  status )
{
   _braid_Action   *action;
   braid_BaseVector u;
   braid_VectorBar  ubar;
   braid_VectorBar  ubar_copy;
   braid_Int        myid         = _braid_CoreElt(core, myid);
   braid_Int        verbose_adj  = _braid_CoreElt(core, verbose_adj);
   braid_Int        adjoint      = _braid_CoreElt(core, adjoint);
   braid_Int        record       = _braid_CoreElt(core, record);
   braid_Int        receiver     = _braid_CoreElt(core, send_recv_rank);
   braid_Real       timer        = 0.0;

   if ( verbose_adj ) _braid_printf("%d: BUFUNPACK\n", myid);

   /* Allocate the braid_BaseVector */
   _braid_BaseVectorCreate(&u);

   if ( _braid_CoreElt(core, delta_correct) )
   {
      /* allocate the basis */
      u->basis = _braid_TAlloc(_braid_Basis, 1);

      /*  see BaseBufPack for information on how the header is written */
      braid_Int  *header, head_size;
      braid_Byte *data_buf;

      header = ((braid_Int*)buffer) + 1;
      u->basis->rank = header[-1] - 1;
      head_size = (header[-1] + 1) * sizeof(braid_Int);  /* total number of vectors + size of each vector */

      data_buf = ((braid_Byte*)buffer) + head_size;  /* buffer needs to have type to enable legal pointer arithmetic */

      /* allocate the user vectors */
      u->basis->userVecs = _braid_TAlloc(braid_Vector, u->basis->rank);

      /* BufUnpack the user's vector */
      timer = _braid_MPI_Wtime(core, 2);
      _braid_CoreFcn(core, bufunpack)(app, (void*)data_buf, &(u->userVector), status);
      data_buf += header[0];

      /* BufUnpack the basis vectors */
      for (braid_Int i = 0; i < u->basis->rank; i++)
      {
         _braid_CoreFcn(core, bufunpack)(app, (void*)data_buf, &(u->basis->userVecs[i]), status);
         data_buf += header[i+1];
      }
      _braid_CoreElt(core, timer_user_bufunpack) += _braid_MPI_Wtime(core, 2) - timer;
   }
   else
   {
      /* BufUnpack the user's vector */
      timer = _braid_MPI_Wtime(core, 2);
      _braid_CoreFcn(core, bufunpack)(app, buffer, &(u->userVector), status);
      _braid_CoreElt(core, timer_user_bufunpack) += _braid_MPI_Wtime(core, 2) - timer;
   }

   if ( adjoint )
   {
      /* Allocate and initialize the bar vector with zero*/
      ubar = _braid_TAlloc(_braid_VectorBar, 1);
      ubar->useCount = 1;
      _braid_CoreFcn(core, init)(app, -1.0, &(ubar->userVector));
      _braid_CoreFcn(core, sum)(app, -1.0, ubar->userVector, 1.0, ubar->userVector);
      u->bar = ubar;
   }

   /* Record to the tape */
   if ( record )
   {
      /* Set up and push the action */
      action                 = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall      = BUFUNPACK;
      action->core           = core;
      action->send_recv_rank = receiver;
      action->myid           = myid;
      action->messagetype    = _braid_StatusElt(status, messagetype);
      action->size_buffer    = _braid_StatusElt(status, size_buffer);
      action->inTimeIdx      = _braid_StatusElt(status, idx);
      action->level          = _braid_StatusElt(status, level);
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

      /* Copy and push the bar vector to the bar tape */
      _braid_VectorBarCopy(u->bar, &ubar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), ubar_copy);
    }
  
   *u_ptr = u;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseBufAlloc(braid_Core          core,
                    braid_App           app,    
                    void              **buffer,
                    braid_Int           nbytes,
                    braid_BufferStatus  status )
{
   /* Call the user's buffer allocate function */
   if( _braid_CoreFcn(core, bufalloc) != NULL)
   {
      _braid_CoreFcn(core, bufalloc)(app, buffer, nbytes, status);
   }
   else
   {
      *buffer = malloc(nbytes);
   }
   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseBufFree(braid_Core          core,
                   braid_App           app,    
                   void              **buffer)
{
   /* Call the user's buffer free function */
   if( _braid_CoreFcn(core, buffree) != NULL)
   {
      _braid_CoreFcn(core, buffree)(app, buffer);
   }
   else
   {
      free((char *) *buffer);
   }

   *buffer = NULL;
   return _braid_error_flag;
}


/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseObjectiveT(braid_Core             core,
                      braid_App              app,
                      braid_BaseVector       u,
                      braid_ObjectiveStatus  ostatus,
                      braid_Real            *objT_ptr )
{
   _braid_Action   *action;
   braid_Vector     u_copy;
   braid_VectorBar  ubar_copy;
   braid_Int        verbose_adj   = _braid_CoreElt(core, verbose_adj);
   braid_Int        record        = _braid_CoreElt(core, record);
   braid_Int        myid          = _braid_CoreElt(core, myid);
   braid_Real       t             = _braid_CoreElt(core, t);
   braid_Int        idx           = _braid_CoreElt(core, idx);
   braid_Int        iter          = _braid_CoreElt(core, niter);
   braid_Int        level         = _braid_CoreElt(core, level);
   braid_Int        nrefine       = _braid_CoreElt(core, nrefine);
   braid_Int        gupper        = _braid_CoreElt(core, gupper);
   
   if ( verbose_adj ) _braid_printf("%d: OBJECTIVET\n", myid);

   /* if bar: Record to the tape */
   if ( record )
   {
      /* Set up and push the action */
      action             = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall  = OBJECTIVET;
      action->core       = core;
      action->myid       = myid;
      action->inTime     = t;
      action->inTimeIdx  = idx;
      action->braid_iter = iter;
      action->level      = level;
      action->nrefine    = nrefine;
      action->gupper     = gupper;
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

      /* Push a copy of the user's vector to the userVector tape */
      _braid_CoreFcn(core, clone)(app, u->userVector, &u_copy);     // this will accolate memory for the copy!
      _braid_CoreElt(core, userVectorTape) = _braid_TapePush( _braid_CoreElt(core, userVectorTape), u_copy);

      /* Push a copy of the bar vector to the bar tape */
      _braid_VectorBarCopy(u->bar, &ubar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), ubar_copy);
   }

   /* Evaluate the objective function at time t */
   _braid_CoreFcn(core, objectiveT)(app, u->userVector, ostatus, objT_ptr);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseResidual(braid_Core        core,
                    braid_App         app,    
                    braid_BaseVector  ustop,  
                    braid_BaseVector  r,      
                    braid_StepStatus  status )
{
   braid_Int        verbose_adj  = _braid_CoreElt(core, verbose_adj);
   braid_Int        myid         = _braid_CoreElt(core, myid);
   braid_Real       timer        = 0.0;

   if ( verbose_adj ) _braid_printf("%d: RESIDUAL\n", myid);

   /* Call the users Residual function */
   timer = _braid_MPI_Wtime(core, 2);
   _braid_CoreFcn(core, residual)(app, ustop->userVector, r->userVector, status);
   _braid_CoreElt(core, timer_user_residual) += _braid_MPI_Wtime(core, 2) - timer;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseFullResidual(braid_Core        core,
                        braid_App         app,    
                        braid_BaseVector  r,      
                        braid_BaseVector  u,  
                        braid_StepStatus  status )
{
   braid_Int        verbose_adj  = _braid_CoreElt(core, verbose_adj);
   braid_Int        myid         = _braid_CoreElt(core, myid);

   if ( verbose_adj ) _braid_printf("%d: FULLRESIDUAL\n", myid);

   /* Call the users Residual function */
   _braid_CoreFcn(core, full_rnorm_res)(app, r->userVector, u->userVector, status);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseSCoarsen(braid_Core              core,
                    braid_App               app,   
                    braid_BaseVector        fu,    
                    braid_BaseVector       *cu_ptr,
                    braid_CoarsenRefStatus  status )
{
   braid_BaseVector cu;
   braid_Int        verbose_adj  = _braid_CoreElt(core, verbose_adj);
   braid_Int        myid         = _braid_CoreElt(core, myid);
   braid_Real       timer        = 0.0;

   if ( verbose_adj ) _braid_printf("%d: SCOARSEN\n", myid);

   cu = _braid_TAlloc(_braid_BaseVector, 1);
   cu->userVector = NULL;
   cu->bar        = NULL;
   cu->basis      = NULL;

   /* Call the users SCoarsen Function */
   timer = _braid_MPI_Wtime(core, 2);
   _braid_CoreFcn(core, scoarsen)(app, fu->userVector, &(cu->userVector), status);
   _braid_CoreElt(core, timer_user_scoarsen) += _braid_MPI_Wtime(core, 2) - timer;

   *cu_ptr = cu;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseSRefine(braid_Core              core,
                   braid_App               app,    
                   braid_BaseVector        cu,     
                   braid_BaseVector       *fu_ptr, 
                   braid_CoarsenRefStatus  status )
{
   braid_BaseVector fu;
   braid_Int        verbose_adj  = _braid_CoreElt(core, verbose_adj);
   braid_Int        myid         = _braid_CoreElt(core, myid);
   braid_Real       timer        = 0.0;

   if ( verbose_adj ) _braid_printf("%d: SREFINE\n", myid);

   fu = _braid_TAlloc(_braid_BaseVector, 1);
   fu->userVector = NULL;
   fu->bar        = NULL;
   fu->basis      = NULL;

   /* Call the users SRefine */
   timer = _braid_MPI_Wtime(core, 2);
   _braid_CoreFcn(core, srefine)(app, cu->userVector, &(fu->userVector), status);
   _braid_CoreElt(core, timer_user_srefine) += _braid_MPI_Wtime(core, 2) - timer;

   *fu_ptr = fu;

   return _braid_error_flag;
}                      

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseSInit(braid_Core        core,
                 braid_App         app,  
                 braid_Real        t,    
                 braid_BaseVector *u_ptr )
{
   braid_BaseVector u;
   braid_Int        verbose_adj  = _braid_CoreElt(core, verbose_adj);
   braid_Int        myid         = _braid_CoreElt(core, myid);

   if ( verbose_adj ) _braid_printf("%d: SINIT\n", myid);

   u = _braid_TAlloc(_braid_BaseVector, 1);

   /* Call the users SInit */
   _braid_CoreFcn(core, sinit)(app, t, &(u->userVector));

   *u_ptr = u;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseSClone(braid_Core        core, 
                 braid_App          app,  
                 braid_BaseVector   u,    
                 braid_BaseVector  *v_ptr )
{

   braid_BaseVector v;
   braid_Int        verbose_adj  = _braid_CoreElt(core, verbose_adj);
   braid_Int        myid         = _braid_CoreElt(core, myid);

   if ( verbose_adj ) _braid_printf("%d: SCLONE\n", myid);

   v = _braid_TAlloc(_braid_BaseVector, 1);

   /* Call the users SClone */
   _braid_CoreFcn(core, sclone)(app, u->userVector, &(v->userVector));

   *v_ptr = v;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseSFree(braid_Core      core,
                 braid_App        app,
                 braid_BaseVector u )
{
 
   braid_Int  verbose_adj  = _braid_CoreElt(core, verbose_adj);
   braid_Int  myid         = _braid_CoreElt(core, myid);
 
   if ( verbose_adj ) _braid_printf("%d: SFREE\n", myid);

   /* Call the users sfree */
   _braid_CoreFcn(core, sfree)(app, u->userVector);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseTimeGrid(braid_Core  core,
                    braid_App   app,    
                    braid_Real *ta,     
                    braid_Int  *ilower, 
                    braid_Int  *iupper )
{
   braid_Int  verbose_adj  = _braid_CoreElt(core, verbose_adj);
   braid_Int  myid         = _braid_CoreElt(core, myid);
 
   if ( verbose_adj ) _braid_printf("%d: TIMEGRID\n", myid);

   /* Call the users timegrid function */
   _braid_CoreFcn(core, tgrid)(app, ta, ilower, iupper);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseStep_diff(_braid_Action *action)
{
   braid_Vector     u, ustop;
   braid_VectorBar  ubar, ustopbar;
   braid_Core       core         = action->core;
   braid_StepStatus status       = (braid_StepStatus) action->core;
   braid_Real       inTime       = action->inTime;
   braid_Real       outTime      = action->outTime;
   braid_Int        tidx         = action->inTimeIdx;
   braid_Real       tol          = action->tol;
   braid_Int        iter         = action->braid_iter;
   braid_Int        level        = action->level;
   braid_Int        nrefine      = action->nrefine;
   braid_Int        gupper       = action->gupper;
   braid_App        app          = _braid_CoreElt(core, app);
   braid_Int        verbose_adj  = _braid_CoreElt(core, verbose_adj);
   braid_Int        myid         = _braid_CoreElt(core, myid);

   if ( verbose_adj ) _braid_printf("%d: STEP_DIFF %.4f to %.4f, %d\n", myid, inTime, outTime, tidx);

   /* Pop ustop & u from primal tape */
   ustop = (braid_Vector)    (_braid_CoreElt(core, userVectorTape)->data_ptr);
   _braid_CoreElt(core, userVectorTape) = _braid_TapePop( _braid_CoreElt(core, userVectorTape) );
   u = (braid_Vector)    (_braid_CoreElt(core, userVectorTape)->data_ptr);
   _braid_CoreElt(core, userVectorTape) = _braid_TapePop( _braid_CoreElt(core, userVectorTape) );

   /* Pop ustopbar & ubar from bar tape */
   ustopbar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
   _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );
   ubar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
   _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );


   /* Set up the status structure */
   _braid_StepStatusInit(inTime, outTime, tidx, tol, iter, level, nrefine, gupper, braid_ASCaller_BaseStep_diff, NULL, status);

   /* Call the users's differentiated step function */
   _braid_CoreFcn(core, step_diff)(app, ustop, u, ustopbar->userVector, ubar->userVector, status);

   /* Free memory of the primal and bar vectors */
   _braid_VectorBarDelete(core, ubar);
   _braid_VectorBarDelete(core, ustopbar);
   _braid_CoreFcn(core, free)(app, u);
   _braid_CoreFcn(core, free)(app, ustop);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseClone_diff(_braid_Action *action)
{
   braid_VectorBar v_bar;
   braid_VectorBar u_bar;
   braid_Core      core        = action->core;
   braid_Int       verbose_adj = _braid_CoreElt(core, verbose_adj);
   braid_Int       myid        = _braid_CoreElt(core, myid);

   if ( verbose_adj ) _braid_printf("%d: CLONE_DIFF\n", myid);

   /* Get and pop vbar from the tape */
   v_bar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
   _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );

   /* Get and pop ubar from the tape */
   u_bar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
   _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );

   /* Perform the differentiated clone action :
   *  ub += vb
   *  vb  = 0.0
   */
   _braid_CoreFcn(core, sum)(_braid_CoreElt(core,app),1.0, v_bar->userVector,  1.0, u_bar->userVector);
   _braid_CoreFcn(core, sum)(_braid_CoreElt(core,app),1.0, v_bar->userVector, -1.0, v_bar->userVector);

   /* Free the bar vectors */
   _braid_VectorBarDelete(core, u_bar);
   _braid_VectorBarDelete(core, v_bar);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseSum_diff(_braid_Action *action)
{
   braid_VectorBar y_bar;
   braid_VectorBar x_bar;
   braid_Core      core        = action->core;
   braid_Real      alpha       = action->sum_alpha;
   braid_Real      beta        = action->sum_beta;
   braid_App       app         = _braid_CoreElt(core, app);
   braid_Int       verbose_adj = _braid_CoreElt(core, verbose_adj);
   braid_Int       myid        = _braid_CoreElt(core, myid);

   if ( verbose_adj ) _braid_printf("%d: SUM_DIFF\n", myid);

   /* Get and pop ybar from the tape */
   y_bar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
   _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );

   /* Get and pop ubar from the tape */
   x_bar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
   _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );

   /* Perform the differentiated sum action: 
   *  xb += alpha * yb
   *  yb  = beta  * yb
   */
   _braid_CoreFcn(core, sum)(app, alpha, y_bar->userVector,  1.0, x_bar->userVector);
   _braid_CoreFcn(core, sum)(app,   0.0, x_bar->userVector, beta, y_bar->userVector);

   /* Free the bar vectors */
   _braid_VectorBarDelete(core, y_bar);
   _braid_VectorBarDelete(core, x_bar);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseObjectiveT_diff(_braid_Action *action)
{
   braid_Vector           u; 
   braid_VectorBar        ubar; 
   braid_Int              myid         = action->myid;
   braid_Core             core         = action->core;
   braid_Real             t            = action->inTime;
   braid_Int              idx          = action->inTimeIdx;
   braid_Int              iter         = action->braid_iter;
   braid_Int              level        = action->level;
   braid_Int              nrefine      = action->nrefine;
   braid_Int              gupper       = action->gupper;
   braid_App              app          = _braid_CoreElt(core, app);
   braid_Int              verbose_adj  = _braid_CoreElt(core, verbose_adj);
   braid_Real             f_bar        = _braid_CoreElt(core, optim)->f_bar;
   braid_ObjectiveStatus  ostatus      = (braid_ObjectiveStatus) action->core;

   if ( verbose_adj ) _braid_printf("%d: OBJT_DIFF\n", myid);

   /* Get the primal and bar vectors from the tapes */
   u    = (braid_Vector)    (_braid_CoreElt(core, userVectorTape)->data_ptr);
   ubar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);

   /* Pop the vectors from the tapes */
   _braid_CoreElt(core, userVectorTape) = _braid_TapePop( _braid_CoreElt(core, userVectorTape) );
   _braid_CoreElt(core, barTape)        = _braid_TapePop( _braid_CoreElt(core, barTape) );

   /* Store the values of the adjoint */
   braid_Vector userbarCopy;
   _braid_CoreFcn(core, clone)(app, ubar->userVector, &userbarCopy);

  /* Call the users's differentiated objective function */
   _braid_ObjectiveStatusInit(t, idx, iter, level, nrefine, gupper, ostatus);
   _braid_CoreFcn(core, objT_diff)( app, u, ubar->userVector, f_bar, ostatus);

   /* Add the stored value */
   _braid_CoreFcn(core, sum)(app, 1., userbarCopy, 1., ubar->userVector);

   /* Free the stores value */
   _braid_CoreFcn(core, free)(app, userbarCopy);

   /* Free primal and bar vectors */
   _braid_CoreFcn(core, free)(app, u);
   _braid_VectorBarDelete(core, ubar);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseBufPack_diff(_braid_Action *action )
{
   braid_Int          size;
   void              *buffer;
   braid_Vector       u;
   braid_VectorBar    ubar;
   braid_Core         core            = action->core;
   braid_Real         send_recv_rank  = action->send_recv_rank;
   braid_Int          messagetype     = action->messagetype;
   braid_Int          size_buffer     = action->size_buffer;
   braid_Int          level           = action->level;
   braid_Int          index           = action->inTimeIdx;
   braid_App          app             = _braid_CoreElt(core, app);
   braid_Int          verbose_adj     = _braid_CoreElt(core, verbose_adj);
   braid_Int          myid            = _braid_CoreElt(core, myid);
   braid_BufferStatus bstatus         = (braid_BufferStatus) core;

   if ( verbose_adj ) _braid_printf("%d: BUFPACK_DIFF\n", myid);

   /* Get the bar vector and pop it from the tape*/
   ubar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
   _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );

   /* Allocate the buffer */
   _braid_CoreFcn(core, bufsize)(app, &size, bstatus);
   buffer = malloc(size);

   /* Receive the buffer */
   MPI_Recv(buffer, size, MPI_BYTE, send_recv_rank, 0, _braid_CoreElt(core, comm), MPI_STATUS_IGNORE); 

   /* Initialize the bstatus */
   _braid_BufferStatusInit( messagetype, index, level, size_buffer, bstatus);

   /* Unpack the buffer into u */
   _braid_CoreFcn(core, bufunpack)(app, buffer, &u, bstatus);

   /* Update ubar with u */
   _braid_CoreFcn(core, sum)( app, 1., u, 1., ubar->userVector);

   /* Free the vectors */
   _braid_VectorBarDelete(core, ubar);
   _braid_CoreFcn(core, free)(app, u);
   free(buffer);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_BaseBufUnpack_diff(_braid_Action *action)
{
   braid_VectorBar     ubar;
   void               *sendbuffer;
   braid_Int           size;
   MPI_Request        *requests;
   braid_Core          core           = action->core;;
   braid_Real          send_recv_rank = action->send_recv_rank;
   braid_Int           messagetype    = action->messagetype;
   braid_Int           size_buffer    = action->size_buffer;
   braid_Int           level          = action->level;
   braid_Int           index          = action->inTimeIdx;
   braid_BufferStatus  bstatus        = (braid_BufferStatus) core;
   braid_App           app            = _braid_CoreElt(core, app);
   braid_Int           verbose_adj    = _braid_CoreElt(core, verbose_adj);
   braid_Int           myid           = _braid_CoreElt(core, myid);

   if ( verbose_adj ) _braid_printf("%d: BUFUNPACK_DIFF\n", myid);

   /* Get the bar vector and pop it from the tape*/
   ubar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
   _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );

   /* Get the buffer size */
   _braid_CoreFcn(core, bufsize)(app, &size, bstatus);

  /* Wait for previous send request to finish */
   MPI_Request *request = _braid_CoreElt(core, optim)->request;
   MPI_Status *mpistatus = malloc(sizeof(MPI_Status));
   if (request != NULL)
   {
     MPI_Wait(request, mpistatus);
   }

   /* Initialize the bufferstatus */
   _braid_BufferStatusInit( messagetype, index, level, size_buffer, bstatus);

   /* Pack the buffer */
   sendbuffer = _braid_CoreElt(core, optim)->sendbuffer;
   _braid_CoreFcn(core, bufpack)( app, ubar->userVector, sendbuffer, bstatus);

   /* Send the buffer  */
   requests = _braid_CTAlloc(MPI_Request, 1);
   MPI_Isend(sendbuffer, size, MPI_BYTE, send_recv_rank, 0, _braid_CoreElt(core, comm), &requests[0]);

   /* Store the request */
   _braid_CoreElt(core, optim)->request = requests;
   
   /* Set ubar to zero */
   _braid_CoreFcn(core, sum)(app, -1., ubar->userVector, 1., ubar->userVector );

   /* Free ubar */
   _braid_VectorBarDelete(core, ubar);

   return _braid_error_flag;
}
