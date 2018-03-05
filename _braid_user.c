/**
 *  Source file implementing the wrapper for the user routines
 **/


#ifndef _braid_user_HEADER
#define _braid_user_HEADER

#include "_braid.h"


braid_Int 
_braid_UserStep(braid_Core       core,
                braid_App        app,    
                braid_Vector     ustop,
                braid_Vector     fstop, 
                braid_Vector     u, 
                braid_StepStatus status )
{

   /* if adjoint: Record to the tape */
   if ( _braid_CoreElt(core, adjoint) )
   {
      /* Set up the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = STEP;
      action->inTime        = _braid_CoreElt(core, t);
      action->outTime       = _braid_CoreElt(core, tnext);
      action->status        = (braid_Status) status;
      action->myid          = _braid_CoreElt(core, myid);

      
      /* Push the action to the actiontape */
      _braid_CoreElt(core, actiontape) = _braid_TapePush( _braid_CoreElt(core, actiontape) , action);

      /* delete action ? */
   
   }

   /* Call the users Step function */
   _braid_CoreFcn(core, step)(app, ustop, fstop, u, status);
   
   return 0;
}


                        
braid_Int
_braid_UserInit(braid_Core core,
                braid_App      app, 
                braid_Real     t,   
                braid_Vector  *u_ptr
                )
{
   /* if adjoint: Record to the tape */
   if ( _braid_CoreElt(core, adjoint) )
   {
      /* Set up the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = INIT;
      action->inTime        = _braid_CoreElt(core, t);
      action->myid          = _braid_CoreElt(core, myid);
      
      /* Push the action to the action tape */
      _braid_CoreElt(core, actiontape) = _braid_TapePush( _braid_CoreElt(core, actiontape) , action);

   }
   
   /* Call the users Init function */
   _braid_CoreFcn(core, init)(app, t, u_ptr);

   return 0;
}

braid_Int
_braid_UserClone(braid_Core core,
                 braid_App app,  
                 braid_Vector   u,    
                 braid_Vector  *v_ptr 
                 )
{

   /* if adjoint: Record to the tape */
   if ( _braid_CoreElt(core, adjoint) )
   {
      /* Set up the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = CLONE;
      action->myid          = _braid_CoreElt(core, myid);
      
      /* Push the action to the action tape */
      _braid_CoreElt(core, actiontape) = _braid_TapePush( _braid_CoreElt(core, actiontape) , action);
   }
   
   /* Call the users Clone function */
   _braid_CoreFcn(core, clone)(app, u, v_ptr);

   return 0;
}


braid_Int
_braid_UserFree(braid_Core core,
                braid_App     app,
                braid_Vector  u   
                )
{


   /* if adjoint: Record to the tape */
   if ( _braid_CoreElt(core, adjoint) )
   {
      /* Set up the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = FREE;
      action->myid          = _braid_CoreElt(core, myid);
      
      /* Push the action to the action tape */
      _braid_CoreElt(core, actiontape) = _braid_TapePush( _braid_CoreElt(core, actiontape) , action);
   }
 
   /* Call the users free function */
   _braid_CoreFcn(core, free)(app, u);

   return 0;
}


braid_Int
_braid_UserSum(braid_Core core,
               braid_App        app,    
               braid_Real    alpha,  
               braid_Vector  x,      
               braid_Real    beta,   
               braid_Vector  y       
               )
{

   /* if adjoint: Record to the tape */
   if ( _braid_CoreElt(core, adjoint) )
   {
      /* Set up the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = SUM;
      action->sum_alpha     = alpha;
      action->sum_beta      = beta;
      action->myid          = _braid_CoreElt(core, myid);
      
      /* Push the action to the action tape */
      _braid_CoreElt(core, actiontape) = _braid_TapePush( _braid_CoreElt(core, actiontape) , action);
   }

    /* Call the users Sum function */
   _braid_CoreFcn(core, sum)(app, alpha, x, beta, y);

   return 0;
}


braid_Int
_braid_UserSpatialNorm(braid_Core core,
                       braid_App      app,      /**< user-defined _braid_App structure */
                       braid_Vector   u,        /**< vector to norm */
                       braid_Real    *norm_ptr  /**< output, norm of braid_Vector (this is a spatial norm) */ 
                       )
{

   /* Call the users SpatialNorm function */
   _braid_CoreFcn(core, spatialnorm)(app, u, norm_ptr);

   return 0;
}


braid_Int
_braid_UserAccess(braid_Core core,
                  braid_App           app,   
                  braid_Vector        u,     
                  braid_AccessStatus  status 
                  )
{

   /* if adjoint: Record to the tape */
   if ( _braid_CoreElt(core, adjoint) )
   {
      /* Set up the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = ACCESS;
      action->inTime        = _braid_CoreElt(core, t);
      action->myid          = _braid_CoreElt(core, myid);
      
      /* Push the action to the action tape */
      _braid_CoreElt(core, actiontape) = _braid_TapePush( _braid_CoreElt(core, actiontape) , action);
   }

   /* Call the users Access function */
   _braid_CoreFcn(core, access)(app, u, status);

   return 0;
}


braid_Int
_braid_UserBufSize(braid_Core core,
                   braid_App   app,               /**< user-defined _braid_App structure */
                   braid_Int  *size_ptr,           /**< upper bound on vector size in bytes */
                   braid_BufferStatus  status     /**< can be querried for info on the message type */
                   ) 
{

   /* Call the users BufSize function */
   _braid_CoreFcn(core, bufsize)(app, size_ptr, status);

   return 0;
}


braid_Int
_braid_UserBufPack(braid_Core core,
                   braid_App           app,       
                   braid_Vector        u,         
                   void               *buffer,    
                   braid_BufferStatus  status     
                   )
{
      /* if adjoint: Record to the tape */
   if ( _braid_CoreElt(core, adjoint) )
   {
      /* Set up the action */
      _braid_Action* action  = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall      = BUFPACK;
      action->send_recv_rank = _braid_CoreElt(core, send_recv_rank);
      action->myid           = _braid_CoreElt(core, myid);
      
      /* Push the action to the action tape */
      _braid_CoreElt(core, actiontape) = _braid_TapePush( _braid_CoreElt(core, actiontape) , action);
   }
   
   /* Call the users BufPack function */
   _braid_CoreFcn(core, bufpack)(app, u, buffer, status);

   return 0;
}


braid_Int
_braid_UserBufUnpack(braid_Core core,
                     braid_App            app,    
                     void                *buffer, 
                     braid_Vector        *u_ptr,  
                     braid_BufferStatus   status  
                     )
{
      /* if adjoint: Record to the tape */
   if ( _braid_CoreElt(core, adjoint) )
   {
      /* Set up the action */
      _braid_Action* action  = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall      = BUFUNPACK;
      action->send_recv_rank = _braid_CoreElt(core, send_recv_rank);
      action->myid           = _braid_CoreElt(core, myid);
      
      /* Push the action to the action tape */
      _braid_CoreElt(core, actiontape) = _braid_TapePush( _braid_CoreElt(core, actiontape) , action);
   }
   /* Call the users BufUnPack function */
   _braid_CoreFcn(core, bufunpack)(app, buffer, u_ptr, status);

   return 0;
}


#endif