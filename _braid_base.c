/**
 *  Source file implementing the wrapper for the user routines
 **/


#ifndef _braid_base_HEADER
#define _braid_base_HEADER

#include "_braid.h"
#include "_braid_status.h"
#include "_util.c"


braid_Int 
_braid_BaseStep(braid_Core       core,
                braid_App        app,    
                braid_BaseVector ustop,
                braid_BaseVector fstop, 
                braid_BaseVector u, 
                braid_StepStatus status )
{
   _braid_Action   *action;
   braid_Vector     u_copy;
   braid_VectorBar  bar_copy;
   braid_Real       t         = _braid_CoreElt(core, t);
   braid_Real       tnext     = _braid_CoreElt(core, tnext);
   braid_Int        myid      = _braid_CoreElt(core, myid);
   braid_Int        verbose   = _braid_CoreElt(core, verbose);
   braid_Int        record    = _braid_CoreElt(core, record);
   braid_Int        tidx      = _braid_CoreElt(core, idx);

   if (verbose) printf("%d: STEP %.4f to %.4f, %d\n", myid, t, tnext, tidx);

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
      action->status     = (braid_Status) status;
      action->myid       = myid;
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

      /* Push a copy of the primal vector to the primal tape */
      _braid_CoreFcn(core, clone)(app, u->userVector, &u_copy);  // this will accolate memory for the copy!
      _braid_CoreElt(core, userVectorTape) = _braid_TapePush( _braid_CoreElt(core, userVectorTape), u_copy);

      /* Push a copy of the bar vector to the bartape */
      _braid_VectorBarCopy(u->bar, &bar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), bar_copy);
  }

   /* Call the users Step function */
   if ( _braid_CoreElt(core, residual) == NULL )
   {
      _braid_CoreFcn(core, step)(app, ustop->userVector, NULL, u->userVector, status);
   }
   else
   {
      /* TODO: Check the fstop feature! */
      _braid_CoreFcn(core, step)(app, ustop->userVector, fstop->userVector, u->userVector, status);
   }
   return 0;
}


                        
braid_Int
_braid_BaseInit(braid_Core        core,
                braid_App         app, 
                braid_Real        t,   
                braid_BaseVector *u_ptr )
{
   _braid_Action    *action;
   braid_BaseVector  u;
   braid_VectorBar   ubar;
   braid_Int         myid    = _braid_CoreElt(core, myid);
   braid_Int         verbose = _braid_CoreElt(core, verbose);
   braid_Int         record  = _braid_CoreElt(core, record);

   if (verbose) printf("%d INIT\n", myid);

   /* Allocate the braid_BaseVector */
   u = (braid_BaseVector)malloc(sizeof(braid_Vector)+sizeof(braid_VectorBar));
   u->userVector = NULL;
   u->bar        = NULL;

   /* Allocate and initialize the userVector */
   _braid_CoreFcn(core, init)(app, t, &(u->userVector));
   
   /* Allocate and initialize the bar vector */
   if ( _braid_CoreElt(core, adjoint) )
   {
      ubar = (braid_VectorBar)malloc(sizeof(braid_Vector)+sizeof(int));
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

   return 0;
}

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
   braid_Int         myid       = _braid_CoreElt(core, myid);
   braid_Int         verbose    = _braid_CoreElt(core, verbose);
   braid_Int         record     = _braid_CoreElt(core, record);

   if (verbose) printf("%d: CLONE\n", myid);

   /* Allocate the braid_BaseVector */
   v = (braid_BaseVector)malloc(sizeof(braid_Vector)+sizeof(braid_VectorBar));
   v->userVector  = NULL;
   v->bar = NULL;

   /* Allocate and and copy the userVector */
   _braid_CoreFcn(core, clone)(app, u->userVector, &(v->userVector) );

   /* Allocate and initialize the bar vector to zero*/
   if (_braid_CoreElt(core, adjoint))
   {
      ubar = (braid_VectorBar)malloc(sizeof(braid_Vector)+sizeof(int));
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

   return 0;
}


braid_Int
_braid_BaseFree(braid_Core       core,
                braid_App        app,
                braid_BaseVector u )
{

   _braid_Action *action;
   braid_Int      myid    = _braid_CoreElt(core, myid);
   braid_Int      verbose = _braid_CoreElt(core, verbose);
   braid_Int      adjoint = _braid_CoreElt(core, adjoint);
   braid_Int      record  = _braid_CoreElt(core, record);

   if (verbose) printf("%d: FREE\n", myid);

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
   _braid_CoreFcn(core, free)(app, u->userVector);

   if ( adjoint )
   {
      /* Free the bar vector */
      _braid_VectorBarDelete(core, u->bar);
   }

   /* Free the braid_BaseVector */
   free(u);

   return 0;
}


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
   braid_Int        myid      =  _braid_CoreElt(core, myid);
   braid_Int        verbose   =  _braid_CoreElt(core, verbose);
   braid_Int        record    =  _braid_CoreElt(core, record);

   if ( verbose ) printf("%d: SUM\n", myid);

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
   _braid_CoreFcn(core, sum)(app, alpha, x->userVector, beta, y->userVector);

   return 0;
}


braid_Int
_braid_BaseSpatialNorm(braid_Core        core,
                       braid_App         app,      
                       braid_BaseVector  u,    
                       braid_Real       *norm_ptr )
{

   /* Compute the spatial norm of the user's vector */
   _braid_CoreFcn(core, spatialnorm)(app, u->userVector, norm_ptr);

   return 0;
}


braid_Int
_braid_BaseAccess(braid_Core          core,
                  braid_App           app,   
                  braid_BaseVector    u,     
                  braid_AccessStatus  status )
{
   _braid_Action   *action;
   braid_Real       t         = _braid_CoreElt(core, t);
   braid_Int        myid      = _braid_CoreElt(core, myid);
   braid_Int        verbose   = _braid_CoreElt(core, verbose);
   braid_Int        record    = _braid_CoreElt(core, record);
   
   if ( verbose ) printf("%d: ACCESS\n", myid);

   /* Record to the tape */
   if ( record )
   {
      /* Set up and push the action */
      action             = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall  = ACCESS;
      action->core       = core;
      action->status     = (braid_Status) status;
      action->inTime     = t;
      action->myid       = myid;
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);
   }

   /* Access the user's vector */
   _braid_CoreFcn(core, access)(app, u->userVector, status);

   return 0;
}


braid_Int
_braid_BaseBufSize(braid_Core          core,
                   braid_App           app,      
                   braid_Int          *size_ptr, 
                   braid_BufferStatus  status ) 
{
   braid_Int  myid     = _braid_CoreElt(core, myid);
   braid_Int  verbose  = _braid_CoreElt(core, verbose);

   if ( verbose ) printf("%d: BUFSIZE\n", myid);

   /* Call the users BufSize function */
   _braid_CoreFcn(core, bufsize)(app, size_ptr, status);

   return 0;
}


braid_Int
_braid_BaseBufPack(braid_Core          core,
                   braid_App           app,       
                   braid_BaseVector    u,         
                   void               *buffer,    
                   braid_BufferStatus  status )
{
   _braid_Action   *action;
   braid_VectorBar  ubar_copy;
   braid_Int        myid     = _braid_CoreElt(core, myid);
   braid_Int        verbose  = _braid_CoreElt(core, verbose);
   braid_Int        record   = _braid_CoreElt(core, record);
   braid_Int        sender   = _braid_CoreElt(core, send_recv_rank);

   if ( verbose ) printf("%d: BUFPACK\n",  myid );

   /* Record to the tape */
   if ( record )
   {
      /* Set up and push the action */
      action                 = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall      = BUFPACK;
      action->core           = core;
      action->send_recv_rank = sender; 
      action->myid           = myid;
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

      /* Copy and push the bar pointer to the bar tape */
      _braid_VectorBarCopy(u->bar, &ubar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), ubar_copy);
   }
   
   /* BufPack the user's vector */
   _braid_CoreFcn(core, bufpack)(app, u->userVector, buffer, status);

   return 0;
}


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
   braid_Int        myid     = _braid_CoreElt(core, myid);
   braid_Int        verbose  = _braid_CoreElt(core, verbose);
   braid_Int        adjoint  = _braid_CoreElt(core, adjoint);
   braid_Int        record   = _braid_CoreElt(core, record);
   braid_Int        receiver = _braid_CoreElt(core, send_recv_rank);
   braid_Real       tstart   = _braid_CoreElt(core, tstart);

   if ( verbose ) printf("%d: BUFUNPACK\n", myid);

   /* Allocate the braid_BaseVector */
   u = (braid_BaseVector)malloc(sizeof(braid_Vector)+sizeof(braid_VectorBar));
   u->userVector  = NULL;
   u->bar = NULL;

   /* BufUnpack the user's vector */
   _braid_CoreFcn(core, bufunpack)(app, buffer, &(u->userVector), status);

   if ( adjoint )
   {
      /* Allocate and initialize the bar vector with zero*/
      ubar = (braid_VectorBar)malloc(sizeof(braid_Vector)+sizeof(int));
      ubar->useCount = 1;
      _braid_CoreFcn(core, init)(app, tstart, &(ubar->userVector));
      _braid_CoreFcn(core, sum)(app, -1.0, ubar->userVector, 1.0, ubar->userVector);
      u->bar = ubar;
   }

   /* Record to the tape */
   if (  record )
   {
      /* Set up and push the action */
      action                 = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall      = BUFUNPACK;
      action->core           = core;
      action->send_recv_rank = receiver;
      action->myid           = myid;
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

      /* Copy and push the bar vector to the bar tape */
      _braid_VectorBarCopy(u->bar, &ubar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), ubar_copy);
    }
  
   *u_ptr = u;

   return 0;
}


braid_Int
_braid_BaseObjectiveT(braid_Core          core,
                      braid_App           app,
                      braid_BaseVector    u,
                      braid_AccessStatus  astatus,
                      braid_Real         *objT_ptr )
{
   _braid_Action   *action;
   braid_Vector     u_copy;
   braid_VectorBar  ubar_copy;
   braid_Int        verbose  = _braid_CoreElt(core, verbose);
   braid_Int        record   = _braid_CoreElt(core, record);
   braid_Int        myid     = _braid_CoreElt(core, myid);
   braid_Real       t;
   braid_Int        tidx;
   
   if ( verbose ) printf("%d: OBJECTIVET\n", myid);

   /* if bar: Record to the tape */
   if ( record )
   {
      braid_AccessStatusGetT(astatus, &t);
      braid_AccessStatusGetTIndex(astatus, &tidx);

      /* Set up and push the action */
      action            = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall = OBJECTIVET;
      action->inTime    = t;
      action->inTimeIdx = tidx;
      action->core      = core;
      action->myid      = myid;
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

      /* Push a copy of the user's vector to the userVector tape */
      _braid_CoreFcn(core, clone)(app, u->userVector, &u_copy);     // this will accolate memory for the copy!
      _braid_CoreElt(core, userVectorTape) = _braid_TapePush( _braid_CoreElt(core, userVectorTape), u_copy);

      /* Push a copy of the bar vector to the bar tape */
      _braid_VectorBarCopy(u->bar, &ubar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), ubar_copy);
   }

   /* Evaluate the objective function at time t */
   _braid_CoreFcn(core, objectiveT)(app, u->userVector, astatus,objT_ptr);

   return 0;
}


braid_Int
_braid_BaseResidual(braid_Core        core,
                    braid_App         app,    
                    braid_BaseVector  ustop,  
                    braid_BaseVector  r,      
                    braid_StepStatus  status )
{
   braid_Int        verbose  = _braid_CoreElt(core, verbose);
   braid_Int        myid     = _braid_CoreElt(core, myid);

   if ( verbose ) printf("%d: RESIDUAL\n", myid);

      /* Call the users Residual function */
   _braid_CoreFcn(core, residual)(app, ustop->userVector, r->userVector, status);

   return 0;
}


braid_Int
_braid_BaseSCoarsen(braid_Core              core,
                    braid_App               app,   
                    braid_BaseVector        fu,    
                    braid_BaseVector       *cu_ptr,
                    braid_CoarsenRefStatus  status )
{
   braid_BaseVector cu;
   braid_Int        verbose  = _braid_CoreElt(core, verbose);
   braid_Int        myid     = _braid_CoreElt(core, myid);

   if ( verbose ) printf("%d: SCOARSEN\n", myid);

   cu = (braid_BaseVector)malloc(sizeof(braid_BaseVector));

   /* Call the users SCoarsen Function */
   _braid_CoreFcn(core, scoarsen)(app, fu->userVector, &(cu->userVector), status);

   *cu_ptr = cu;

   return 0;
}

braid_Int
_braid_BaseSRefine(braid_Core                 core,
                   braid_App                  app,    
                      braid_BaseVector        cu,     
                      braid_BaseVector       *fu_ptr, 
                      braid_CoarsenRefStatus  status )
{
   braid_BaseVector fu;
   braid_Int        verbose  = _braid_CoreElt(core, verbose);
   braid_Int        myid     = _braid_CoreElt(core, myid);

   if ( verbose ) printf("%d: SREFINE\n", myid);

   fu = (braid_BaseVector)malloc(sizeof(braid_BaseVector));

   /* Call the users SRefine */
   _braid_CoreFcn(core, srefine)(app, cu->userVector, &(fu->userVector), status);

   *fu_ptr = fu;

   return 0;
}                      


braid_Int
_braid_BaseSInit(braid_Core        core,
                 braid_App         app,  
                 braid_Real        t,    
                 braid_BaseVector *u_ptr )
{
   braid_BaseVector u;
   braid_Int        verbose  = _braid_CoreElt(core, verbose);
   braid_Int        myid     = _braid_CoreElt(core, myid);

   if ( verbose ) printf("%d: SINIT\n", myid);

   u = (braid_BaseVector)malloc(sizeof(braid_BaseVector));

   /* Call the users SInit */
   _braid_CoreFcn(core, sinit)(app, t, &(u->userVector));

   *u_ptr = u;

   return 0;
}


braid_Int
_braid_BaseSClone(braid_Core        core, 
                 braid_App          app,  
                 braid_BaseVector   u,    
                 braid_BaseVector  *v_ptr )
{

   braid_BaseVector v;
   braid_Int        verbose  = _braid_CoreElt(core, verbose);
   braid_Int        myid     = _braid_CoreElt(core, myid);

   if ( verbose ) printf("%d: SCLONE\n", myid);

   v = (braid_BaseVector)malloc(sizeof(braid_BaseVector));

   /* Call the users SClone */
   _braid_CoreFcn(core, sclone)(app, u->userVector, &(v->userVector));

   *v_ptr = v;

   return 0;
}


braid_Int
_braid_BaseSFree(braid_Core      core,
                 braid_App        app,
                 braid_BaseVector u )
{
 
   braid_Int  verbose  = _braid_CoreElt(core, verbose);
   braid_Int  myid     = _braid_CoreElt(core, myid);
 
   if ( verbose ) printf("%d: SFREE\n", myid);

   /* Call the users sfree */
   _braid_CoreFcn(core, sfree)(app, u->userVector);

   free(u);

   return 0;
}


braid_Int
_braid_BaseTimeGrid(braid_Core  core,
                    braid_App   app,    
                    braid_Real *ta,     
                    braid_Int  *ilower, 
                    braid_Int  *iupper )
{
   braid_Int  verbose  = _braid_CoreElt(core, verbose);
   braid_Int  myid     = _braid_CoreElt(core, myid);
 
   if ( verbose ) printf("%d: TIMEGRID\n", myid);

   /* Call the users timegrid function */
   _braid_CoreFcn(core, tgrid)(app, ta, ilower, iupper);

   return 0;
}


braid_Int
_braid_BaseStep_diff(_braid_Action *action)
{
   braid_Vector     u;
   braid_VectorBar  ubar;
   braid_Core       core     = action->core;
   braid_StepStatus status   = (braid_StepStatus) action->status;
   braid_Real       inTime   = action->inTime;
   braid_Real       outTime  = action->outTime;
   braid_Int        tidx     = action->inTimeIdx;
   braid_App        app      = _braid_CoreElt(core, app);
   braid_Int        verbose  = _braid_CoreElt(core, verbose);
   braid_Int        myid     = _braid_CoreElt(core, myid);

   if ( verbose ) printf("%d: STEP_DIFF %.4f to %.4f, %d\n", myid, inTime, outTime, tidx);

   /* Get the priamal and bar vectors and pop them from the tapes */
   u    = (braid_Vector)    (_braid_CoreElt(core, userVectorTape)->data_ptr);
   ubar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
   _braid_CoreElt(core, barTape)        = _braid_TapePop( _braid_CoreElt(core, barTape) );
   _braid_CoreElt(core, userVectorTape) = _braid_TapePop( _braid_CoreElt(core, userVectorTape) );

   /* Trigger the status with the correct tstart and tstop*/
   _braid_StatusElt(status, t)     = inTime;
   _braid_StatusElt(status, tnext) = outTime;
   _braid_StatusElt(status, idx)   = tidx;

   /* Call the users's differentiated step function */
   _braid_CoreFcn(core, step_diff)(app, u, ubar->userVector, status);

   /* Free memory of the primal and bar vectors */
   _braid_VectorBarDelete(core, ubar);
   _braid_CoreFcn(core, free)(app, u);

   return 0;
}


braid_Int
_braid_BaseClone_diff(_braid_Action *action)
{
   braid_VectorBar v_bar;
   braid_VectorBar u_bar;
   braid_Core      core    = action->core;
   braid_Int       verbose = _braid_CoreElt(core, verbose);
   braid_Int       myid    = _braid_CoreElt(core, myid);

   if ( verbose ) printf("%d: CLONE_DIFF\n", myid);

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

   return 0;
}


braid_Int
_braid_BaseSum_diff(_braid_Action *action)
{
   braid_VectorBar y_bar;
   braid_VectorBar x_bar;
   braid_Core      core  = action->core;
   braid_Real      alpha = action->sum_alpha;
   braid_Real      beta  = action->sum_beta;
   braid_App       app     = _braid_CoreElt(core, app);
   braid_Int       verbose = _braid_CoreElt(core, verbose);
   braid_Int       myid    = _braid_CoreElt(core, myid);

   if ( verbose ) printf("%d: SUM_DIFF\n", myid);

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

   return 0;
}


braid_Int
_braid_BaseObjectiveT_diff(_braid_Action *action)
{
   braid_Vector    u;
   braid_VectorBar ubar;
   braid_Core      core    = action->core;
   braid_Int       t       = action->inTime;
   braid_Int       tidx    = action->inTimeIdx;
   braid_App       app     = _braid_CoreElt(core, app);
   braid_Int       verbose = _braid_CoreElt(core, verbose);
   braid_Int       myid    = _braid_CoreElt(core, myid);
   braid_Real      f_bar   = _braid_CoreElt(core, optim)->f_bar;

   if ( verbose ) printf("%d: OBJT_DIFF\n", myid);

   /* Get the primal and bar vectors from the tapes */
   u    = (braid_Vector)    (_braid_CoreElt(core, userVectorTape)->data_ptr);
   ubar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);

   /* Pop the vectors from the tapes */
   _braid_CoreElt(core, userVectorTape) = _braid_TapePop( _braid_CoreElt(core, userVectorTape) );
   _braid_CoreElt(core, barTape)        = _braid_TapePop( _braid_CoreElt(core, barTape) );

  /* Call the users's differentiated objective function */
   _braid_CoreFcn(core, objT_diff)( app, u, ubar->userVector, f_bar, t, tidx);

   /* Free primal and bar vectors */
   _braid_CoreFcn(core, free)(app, u);
   _braid_VectorBarDelete(core, ubar);

   return 0;
}

braid_Int
_braid_BaseBufPack_diff(_braid_Action *action )
{
   braid_Int          size;
   void              *buffer;
   braid_Vector       u;
   braid_VectorBar    ubar;
   braid_Core         core            = action->core;
   braid_Real         send_recv_rank  = action->send_recv_rank;
   braid_App          app             = _braid_CoreElt(core, app);
   braid_Int          verbose         = _braid_CoreElt(core, verbose);
   braid_Int          myid            = _braid_CoreElt(core, myid);
   braid_BufferStatus bstatus         = (braid_BufferStatus) core;

   if ( verbose ) printf("%d: BUFPACK_DIFF\n", myid);

   /* Get the bar vector and pop it from the tape*/
   ubar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
   _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );

   /* Allocate the buffer */
   _braid_BufferStatusInit( 0, 0, bstatus);
   _braid_CoreFcn(core, bufsize)(app, &size, bstatus);
   buffer = malloc(size);

   /* Receive the buffer 
    * TODO: Add CommHandle / status check!!
    * WHY blocking Recv ? 
    */
   MPI_Recv(buffer, size, MPI_BYTE, send_recv_rank, 0, _braid_CoreElt(core, comm), MPI_STATUS_IGNORE); 

   /* Unpack the buffer into u */
   _braid_CoreFcn(core, bufunpack)(app, buffer, &u, bstatus);

   /* Update ubar with u */
   _braid_CoreFcn(core, sum)( app, 1., u, 1., ubar->userVector);

   /* Free the vectors */
   _braid_VectorBarDelete(core, ubar);
   _braid_CoreFcn(core, free)(app, u);
   free(buffer);

   return 0;
}

braid_Int
_braid_BaseBufUnpack_diff(_braid_Action *action)
{
   braid_VectorBar     ubar;
   void               *buffer;
   braid_Int           size;
   MPI_Request        *requests;
   braid_Core          core           = action->core;;
   braid_Real          send_recv_rank = action->send_recv_rank;
   braid_BufferStatus  bstatus        = (braid_BufferStatus) core;
   braid_App           app            = _braid_CoreElt(core, app);
   braid_Int           verbose        = _braid_CoreElt(core, verbose);
   braid_Int           myid           = _braid_CoreElt(core, myid);

   if ( verbose ) printf("%d: BUFUNPACK_DIFF\n", myid);

   /* Get the bar vector and pop it from the tape*/
   ubar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
   _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );

  /* Allocate buffer */
   _braid_BufferStatusInit( 0, 0, bstatus);
   _braid_CoreFcn(core, bufsize)(app, &size, bstatus);
   buffer = malloc(size); 

   /* Pack the buffer */
   _braid_CoreFcn(core, bufpack)( app, ubar->userVector, buffer, bstatus);

   /* Send the buffer  */
   requests = _braid_CTAlloc(MPI_Request, 1);
   MPI_Isend(buffer, size, MPI_BYTE, send_recv_rank, 0, _braid_CoreElt(core, comm), &requests[0]);
   MPI_Request_free(requests);
   
   /* Set ubar to zero */
   _braid_CoreFcn(core, sum)(app, -1., ubar->userVector, 1., ubar->userVector );

   /* Free ubar */
   _braid_VectorBarDelete(core, ubar);


   return 0;
}


braid_Int
_braid_BaseUpdateDesign(braid_Core  core,
                        braid_Int   *update_flag)
{
   braid_App   app              = _braid_CoreElt(core, app);
   braid_Optim optim            = _braid_CoreElt(core, optim);
   braid_Int   myid             = _braid_CoreElt(core, myid);
   braid_Int   iter             = _braid_CoreElt(core, niter);
   braid_Real  threshold_design = optim->threshold_design;
   braid_Real  objective        = optim->objective;
   braid_Real  rnorm            = optim->rnorm;
   braid_Real  rnorm_adj        = optim->rnorm_adj;
   braid_Real  gnorm            = optim->gnorm;

   /* Return if no update_design function has been specified -> only gradient computation */
   if (_braid_CoreElt(core, update_design) == NULL)
   {
      return 0;
   }

   /* Call the users update_design function, if state and adjoint residual norms are below the tolerance */
   if ( rnorm     < threshold_design && 
        rnorm_adj < threshold_design )
   {
      /* Never update the design at first few iterations because initialization error can be big. */
      if ( !(optim->iter == 0 && iter <= 2))
      {
         /* Update design */
         _braid_CoreFcn(core, update_design)(app, objective, rnorm, rnorm_adj);

         /* Write to optimization outfile */
         _braid_ParFprintfFlush(optim->outfile, myid, "%03d  %1.14e %1.14e %1.14e %1.14e\n", optim->iter, rnorm, rnorm_adj, gnorm, objective);
         /* Increase optimization counter */
         optim->iter++;
         *update_flag = 1;
      }
   }   

   return 0;
}


braid_Int
_braid_BaseComputeGNorm(braid_Core core,
                        braid_Int  iter)
{
   braid_App   app       = _braid_CoreElt(core, app);
   braid_Real  gnorm;

   /* Call the users update_design function */
   _braid_CoreFcn(core, compute_gnorm)(app, &gnorm);

   /* If no design updates are performed, set gradient to zero
    * This ensures, that stopping criterion is based only on state and adjoint residuals. 
    */
   if (_braid_CoreElt(core, update_design) == NULL)
   {
      gnorm = 0.0;
   }

   /* Store the gradient norm in the optim structure */
   _braid_SetGradientNorm(core, iter, gnorm);

   return 0;
}

#endif