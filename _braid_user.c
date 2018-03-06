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
   if ( _braid_CoreElt(core, adjoint) && !_braid_CoreElt(core, done) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = STEP;
      action->inTime        = _braid_CoreElt(core, t);
      action->outTime       = _braid_CoreElt(core, tnext);
      action->status        = (braid_Status) status;
      action->myid          = _braid_CoreElt(core, myid);
      // printf("STEP push\n");
      _braid_CoreElt(core, actiontape) = _braid_TapePush( _braid_CoreElt(core, actiontape) , action);

      /* Push a copy of the primal vector to the tape */
      // braid_Vector u_copy;
      // _braid_CoreFcn(core, clone)(app, u, &u_copy);  // this will accolate memory for the copy!
      // _braid_CoreElt(core, primaltape) = _braid_TapePush( _braid_CoreElt(core, primaltape), u_copy);

      /* Debug info: */
      // _braid_CoreFcn(core, access)(app, u_copy, (braid_AccessStatus) status);
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
   if ( _braid_CoreElt(core, adjoint) && !_braid_CoreElt(core, done) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = INIT;
      action->inTime        = _braid_CoreElt(core, t);
      action->myid          = _braid_CoreElt(core, myid);
      // printf("INIT push: \n");
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
   if ( _braid_CoreElt(core, adjoint) && !_braid_CoreElt(core, done) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = CLONE;
      action->myid          = _braid_CoreElt(core, myid);
      // printf("CLONE push: \n");
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
   if ( _braid_CoreElt(core, adjoint) && !_braid_CoreElt(core, done) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = FREE;
      action->myid          = _braid_CoreElt(core, myid);
      // printf("FREE push\n");
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
   if ( _braid_CoreElt(core, adjoint) && !_braid_CoreElt(core, done) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = SUM;
      action->sum_alpha     = alpha;
      action->sum_beta      = beta;
      action->myid          = _braid_CoreElt(core, myid);
      // printf("SUM push\n");
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
   if ( _braid_CoreElt(core, adjoint) && !_braid_CoreElt(core, done) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = ACCESS;
      action->inTime        = _braid_CoreElt(core, t);
      action->myid          = _braid_CoreElt(core, myid);
      // printf("ACCESS push\n");
      _braid_CoreElt(core, actiontape) = _braid_TapePush( _braid_CoreElt(core, actiontape) , action);

      /* Push a copy of the primal vector to the tape */
      // braid_Vector u_copy;
      // _braid_CoreFcn(core, clone)(app, u, &u_copy);  // this will accolate memory for the copy!
      // _braid_CoreElt(core, primaltape) = _braid_TapePush( _braid_CoreElt(core, primaltape), u_copy);

      // /* Debug info: */
      // _braid_CoreFcn(core, access)(app, u_copy, (braid_AccessStatus) status);
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
   if ( _braid_CoreElt(core, adjoint) && !_braid_CoreElt(core, done) )
   {
      /* Set up and push the action */
      _braid_Action* action  = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall      = BUFPACK;
      action->send_recv_rank = _braid_CoreElt(core, send_recv_rank);
      action->myid           = _braid_CoreElt(core, myid);
      // printf("BufPack push\n");
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
   if ( _braid_CoreElt(core, adjoint) && !_braid_CoreElt(core, done) )
   {
      /* Set up and push the action */
      _braid_Action* action  = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall      = BUFUNPACK;
      action->send_recv_rank = _braid_CoreElt(core, send_recv_rank);
      action->myid           = _braid_CoreElt(core, myid);
      // printf("BufUnpack push\n");
      _braid_CoreElt(core, actiontape) = _braid_TapePush( _braid_CoreElt(core, actiontape) , action);
   }

   /* Call the users BufUnPack function */
   _braid_CoreFcn(core, bufunpack)(app, buffer, u_ptr, status);

   return 0;
}


/*----- Adjoint of user routines ------*/
braid_Int
_braid_UserStepAdjoint( _braid_Action *action, braid_App app)
{
      printf("STEP adjoint pops from the primal tape\n");

      /* Get the braid_vector */
      // braid_Vector vector = (braid_Vector) (head->data_ptr);

      /* Pop from the stack*/
      // head = pop(head);

      /* Free memory of the data */
      // free(ptr);


//    /* DEBUG: Display the vector */
//    _braid_CoreFcn(core, access)(_braid_CoreElt(core, app), vector, astatus )

      return 0;
}

braid_Int
_braid_UserCloneAdjoint(_braid_Action *action)
{
      // printf("CLONE adjoint\n");
      return 0;
}

braid_Int
_braid_UserSumAdjoint(_braid_Action *action)
{
      // printf("SUM adjoint\n");
      return 0;
}


braid_Int
_braid_UserAccessAdjoint(_braid_Action *action, braid_App app)
{
      // printf("ACCESS adjoint\n");
      return 0;
}

braid_Int
_braid_UserBufPackAdjoint(_braid_Action *action, braid_App app)
{
      // printf("BufPack adjoint\n");
      return 0;
}

braid_Int
_braid_UserBufUnpackAdjoint(_braid_Action *action, braid_App app)
{
      // printf("BufUnack adjoint\n");
      return 0;
}
#endif