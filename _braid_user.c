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

    if (_braid_CoreElt(core, verbose)) printf("STEP\n");

   /* if adjoint: Record to the tape */
   if ( _braid_CoreElt(core, record) )
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

      // /* Push a copy of the primal vector to the tape */
      // braid_Vector u_copy;
      // _braid_CoreFcn(core, clone)(app, u, &u_copy);  // this will accolate memory for the copy!
      // _braid_CoreElt(core, primaltape) = _braid_TapePush( _braid_CoreElt(core, primaltape), u_copy);

      /* Debug info: */
      // _braid_CoreFcn(core, access)(app, u_copy, (braid_AccessStatus) status);
  }

   /* Call the users Step function */
   if ( _braid_CoreElt(core, residual) == NULL )
   {
      _braid_CoreFcn(core, step)(app, ustop->primal, NULL, u->primal, status);
   }
   else
   {
      /* TODO: Check the fstop feature!!! */
      _braid_CoreFcn(core, step)(app, ustop->primal, fstop->primal, u->primal, status);
   }
   
   return 0;
}


                        
braid_Int
_braid_UserInit(braid_Core core,
                braid_App      app, 
                braid_Real     t,   
                braid_Vector  *u_ptr
                )
{
   if (_braid_CoreElt(core, verbose)) printf("INIT\n");

   /* if adjoint: Record to the tape */
   if ( _braid_CoreElt(core, record) && !_braid_CoreElt(core, done) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = INIT;
      action->inTime        = _braid_CoreElt(core, t);
      action->myid          = _braid_CoreElt(core, myid);
      // printf("INIT push: \n");
      _braid_CoreElt(core, actiontape) = _braid_TapePush( _braid_CoreElt(core, actiontape) , action);

   }

   /* Allocate memory for the braid_Vector */
   braid_Vector u = (braid_Vector)malloc(sizeof(braid_UserVector)+sizeof(braid_Adjoint));

   /* Allocate and initialize the primal */
   _braid_CoreFcn(core, init)(app, t, &(u->primal));

   /* Allocate and initialize the adjoint with zero*/
   braid_Adjoint myadjoint = (braid_Adjoint)malloc(sizeof(braid_UserVector)+sizeof(int));
   myadjoint->useCount = 1;
   _braid_CoreFcn(core, init)(app, t, &(myadjoint->userVector));
   _braid_CoreFcn(core, sum)(app, -1.0, myadjoint->userVector, 1.0, myadjoint->userVector);
   u->adjoint = myadjoint;

  /* Debug: */
//    printf("Init adjoint, useCount: %d, adjoint value: ", u->adjoint->useCount);
//   braid_AccessStatus astatus = (braid_AccessStatus) core;
//   _braid_CoreFcn(core, access)(app, u->adjoint->userVector, astatus);


   /* Set the pointer */
   *u_ptr = u;

   return 0;
}

braid_Int
_braid_UserClone(braid_Core core,
                 braid_App app,  
                 braid_Vector   u,    
                 braid_Vector  *v_ptr 
                 )
{

   if (_braid_CoreElt(core, verbose)) printf("CLONE\n");

   /* if adjoint: Record to the tape */
   if ( _braid_CoreElt(core, record) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = CLONE;
      action->myid          = _braid_CoreElt(core, myid);
      // printf("CLONE push: \n");
      _braid_CoreElt(core, actiontape) = _braid_TapePush( _braid_CoreElt(core, actiontape) , action);
   }
   
   /* Allocate memory for the braid_Vector */
   braid_Vector v = (braid_Vector)malloc(sizeof(braid_UserVector)+sizeof(braid_Adjoint));

   /* Call the users Clone function for the primal braid_UserVector */
   _braid_CoreFcn(core, clone)(app, u->primal, &(v->primal) );

   /* Allocate and initialize the adjoint to zero*/
   braid_Adjoint myadjoint = (braid_Adjoint)malloc(sizeof(braid_UserVector)+sizeof(int));
   myadjoint->useCount = 1;
   _braid_CoreFcn(core, clone)(app, u->adjoint->userVector, &(myadjoint->userVector));
   _braid_CoreFcn(core, sum)(app, -1.0, myadjoint->userVector, 1.0, myadjoint->userVector);
   v->adjoint = myadjoint;


   *v_ptr = v;

   return 0;
}


braid_Int
_braid_UserFree(braid_Core core,
                braid_App     app,
                braid_Vector  u   
                )
{

   if (_braid_CoreElt(core, verbose)) printf("FREE\n");

   /* if adjoint: Record to the tape */
   if ( _braid_CoreElt(core, record) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = FREE;
      action->myid          = _braid_CoreElt(core, myid);
      // printf("FREE push\n");
      _braid_CoreElt(core, actiontape) = _braid_TapePush( _braid_CoreElt(core, actiontape) , action);
   }
 
   /* Call the users free function for the primal */
   _braid_CoreFcn(core, free)(app, u->primal);


   /* Debug: */
//    printf("Init adjoint, useCount: %d, adjoint value: ", u->adjoint->useCount);
//    braid_AccessStatus astatus = (braid_AccessStatus) core;
//    _braid_CoreFcn(core, access)(app, u->adjoint->userVector, astatus);

   /* Decrease usecount of the adjoint and free the adjoint if useCount == 0. */
   _braid_AdjointDelete(core, u->adjoint);

   /* Free the braid_Vector */
   free(u);

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
   if (_braid_CoreElt(core, verbose)) printf("SUM\n");

   /* if adjoint: Record to the tape */
   if ( _braid_CoreElt(core, record) )
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
   _braid_CoreFcn(core, sum)(app, alpha, x->primal, beta, y->primal);

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
   _braid_CoreFcn(core, spatialnorm)(app, u->primal, norm_ptr);

   return 0;
}


braid_Int
_braid_UserAccess(braid_Core core,
                  braid_App           app,   
                  braid_Vector        u,     
                  braid_AccessStatus  status 
                  )
{
   if (_braid_CoreElt(core, verbose)) printf("ACCESS\n");

   /* if adjoint: Record to the tape */
   if ( _braid_CoreElt(core, record) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = ACCESS;
      action->inTime        = _braid_CoreElt(core, t);
      action->myid          = _braid_CoreElt(core, myid);
      // printf("ACCESS push\n");
      _braid_CoreElt(core, actiontape) = _braid_TapePush( _braid_CoreElt(core, actiontape) , action);

      // /* Push a copy of the primal vector to the tape */
      // braid_Vector u_copy;
      // _braid_CoreFcn(core, clone)(app, u, &u_copy);  // this will accolate memory for the copy!
      // _braid_CoreElt(core, primaltape) = _braid_TapePush( _braid_CoreElt(core, primaltape), u_copy);

      // /* Debug info: */
      // _braid_CoreFcn(core, access)(app, u_copy, (braid_AccessStatus) status);
   }

   /* Call the users Access function */
   _braid_CoreFcn(core, access)(app, u->primal, status);

   return 0;
}


braid_Int
_braid_UserBufSize(braid_Core core,
                   braid_App   app,               /**< user-defined _braid_App structure */
                   braid_Int  *size_ptr,           /**< upper bound on vector size in bytes */
                   braid_BufferStatus  status     /**< can be querried for info on the message type */
                   ) 
{
   if (_braid_CoreElt(core, verbose)) printf("BUFSIZE\n");

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
   if (_braid_CoreElt(core, verbose)) printf("BUFPACK\n");

      /* if adjoint: Record to the tape */
   if ( _braid_CoreElt(core, record) )
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
   _braid_CoreFcn(core, bufpack)(app, u->primal, buffer, status);

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
   if (_braid_CoreElt(core, verbose)) printf("BUFUNPACK\n");

      /* if adjoint: Record to the tape */
   if ( _braid_CoreElt(core, record) )
   {
      /* Set up and push the action */
      _braid_Action* action  = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall      = BUFUNPACK;
      action->send_recv_rank = _braid_CoreElt(core, send_recv_rank);
      action->myid           = _braid_CoreElt(core, myid);
      // printf("BufUnpack push\n");
      _braid_CoreElt(core, actiontape) = _braid_TapePush( _braid_CoreElt(core, actiontape) , action);
   }

   braid_Vector u = (braid_Vector)malloc(sizeof(braid_Vector));

   /* Call the users BufUnPack function */
   _braid_CoreFcn(core, bufunpack)(app, buffer, &(u->primal), status);

   *u_ptr = u;

   return 0;
}

braid_Int
_braid_UserResidual(braid_Core core,
                     braid_App        app,    /**< user-defined _braid_App structure */
                       braid_Vector     ustop,  /**< input, u vector at *tstop* */
                       braid_Vector     r     , /**< output, residual at *tstop* (at input, equals *u* at *tstart*) */
                       braid_StepStatus status  /**< query this struct for info about u (e.g., tstart and tstop) */ 
                       )
{
   if (_braid_CoreElt(core, verbose)) printf("RESIDUAL\n");

      /* Call the users Residual function */
   _braid_CoreFcn(core, residual)(app, ustop->primal, r->primal, status);

   return 0;
}


braid_Int
_braid_UserSCoarsen(braid_Core core,
                    braid_App               app,    /**< user-defined _braid_App structure */
                       braid_Vector            fu,     /**< braid_UserVector to refine*/
                       braid_Vector           *cu_ptr, /**< output, refined vector */   
                       braid_CoarsenRefStatus  status  /**< query this struct for info about fu and cu (e.g., where in time fu and cu are)  */ 
                       )
{
      braid_Vector cu = (braid_Vector)malloc(sizeof(braid_Vector));

      /* Call the users SCoarsen Function */
      _braid_CoreFcn(core, scoarsen)(app, fu->primal, &(cu->primal), status);

      *cu_ptr = cu;

   return 0;
}

braid_Int
_braid_UserSRefine(braid_Core core,
                   braid_App               app,    /**< user-defined _braid_App structure */
                      braid_Vector            cu,     /**< braid_UserVector to refine*/
                      braid_Vector           *fu_ptr, /**< output, refined vector */       
                      braid_CoarsenRefStatus  status  /**< query this struct for info about fu and cu (e.g., where in time fu and cu are)  */ 
                      )
{
      braid_Vector fu = (braid_Vector)malloc(sizeof(braid_Vector));

      /* Call the users SRefine */
      _braid_CoreFcn(core, srefine)(app, cu->primal, &(fu->primal), status);

      *fu_ptr = fu;

   return 0;
}                      

braid_Int
_braid_UserSInit(braid_Core core,
                 braid_App     app,           /**< user-defined _braid_App structure */
                   braid_Real     t,             /**< time value for *u_ptr* */
                   braid_Vector  *u_ptr          /**< output, newly allocated and initialized vector shell */
                   )
{
      braid_Vector u = (braid_Vector)malloc(sizeof(braid_Vector));

      /* Call the users SInit */
      _braid_CoreFcn(core, sinit)(app, t, &(u->primal));

      *u_ptr = u;

   return 0;
}

braid_Int
_braid_UserSClone(braid_Core core, 
                 braid_App      app,          /**< user-defined _braid_App structure */
                    braid_Vector   u,            /**< vector to clone */ 
                    braid_Vector  *v_ptr         /**< output, newly allocated and cloned vector shell */
                    )
{
      braid_Vector v = (braid_Vector)malloc(sizeof(braid_Vector));

      /* Call the users SClone */
      _braid_CoreFcn(core, sclone)(app, u->primal, &(v->primal));

      *v_ptr = v;

   return 0;
}


braid_Int
_braid_UserSFree(braid_Core core,
                  braid_App     app,            /**< user-defined _braid_App structure */
                    braid_Vector  u               /**< vector to free (keeping the shell) */
                    )
{
      /* Call the users sfree */
      _braid_CoreFcn(core, sfree)(app, u->primal);

      free(u);

   return 0;
}

braid_Int
_braid_UserTimeGrid(braid_Core core,
                   braid_App         app,       /**< user-defined _braid_App structure */
                       braid_Real       *ta,        /**< temporal grid on level 0 (slice per processor) */
                       braid_Int        *ilower,    /**< lower time index value for this processor */
                       braid_Int        *iupper     /**< upper time index value for this processor */
                       )
{
      /* Call the users timegrid function */
      _braid_CoreFcn(core, tgrid)(app, ta, ilower, iupper);

   return 0;
}


/*----- Adjoint of user routines ------*/

braid_Int
_braid_UserStepAdjoint(braid_Core core, _braid_Action *action)
{
      // braid_Real        inTime;
      // braid_Real        outTime;
      // braid_StepStatus  status;
      // braid_Int         myid;
      // braid_Vector      vector;

      // /* Grab information from the action */
      // inTime  = action->inTime;
      // outTime = action->outTime;
      // status  = (braid_StepStatus) action->status;
      // myid    = action->myid;

      // /* Get the braid_vector that was used in primal run */
      // vector = (braid_Vector) (_braid_CoreElt(core, primaltape)->data_ptr);

//    /* DEBUG: Display the vector */
//       printf("STEP adjoint pops from the primal tape\n");
//    _braid_CoreFcn(core, access)(_braid_CoreElt(core, app), vector, astatus )


      /* TODO: Call the users's adjoint step function */

      // /* Pop the vector from the primal tape*/
      // _braid_CoreElt(core, primaltape) = _braid_TapePop( _braid_CoreElt(core, primaltape) );

      // /* Free memory of the braid_vector */
      // _braid_CoreFcn(core, free)(_braid_CoreElt(core, app), vector);


      return 0;
}

braid_Int
_braid_UserCloneAdjoint(_braid_Action *action)
{
      // braid_Int         myid;

      // /* Grab information from the action */
      // myid = action->myid;

      // printf("CLONE adjoint\n");
      return 0;
}

braid_Int
_braid_UserSumAdjoint(_braid_Action *action)
{
      // braid_Real alpha;
      // braid_Real beta;
      // braid_Int  myid;

      // /* Grab information from the action */
      // alpha = action->sum_alpha;
      // beta  = action->sum_beta;
      // myid  = action->myid;

      // printf("SUM adjoint\n");
      return 0;
}


braid_Int
_braid_UserAccessAdjoint(braid_Core core, _braid_Action *action)
{
//       braid_Int     inTime;
//       braid_Int     myid;
//       braid_Vector  vector;

//       /* Grab information from the app */
//       inTime = action->inTime;
//       myid   = action->myid;

//       /* Get the braid_vector that was used in primal run */
//       vector = (braid_Vector) (_braid_CoreElt(core, primaltape)->data_ptr);

// //    /* DEBUG: Display the vector */
// //       printf("ACCESS adjoint pops from the primal tape\n");
// //    _braid_CoreFcn(core, access)(_braid_CoreElt(core, app), vector, astatus )


//       /* TODO: Call the users's adjoint access function */

//       /* Pop the vector from the primal tape*/
//       _braid_CoreElt(core, primaltape) = _braid_TapePop( _braid_CoreElt(core, primaltape) );

//       /* Free memory of the braid_vector */
//       _braid_CoreFcn(core, free)(_braid_CoreElt(core, app), vector);


      return 0;
}

braid_Int
_braid_UserBufPackAdjoint(_braid_Action *action, braid_App app)
{
      // braid_Real send_recv_rank;
      // braid_Int  myid;

      // /* Grab information from the action */
      // send_recv_rank = action->send_recv_rank;
      // myid           = action->myid;

      
      // printf("BufPack adjoint\n");
      return 0;
}

braid_Int
_braid_UserBufUnpackAdjoint(_braid_Action *action, braid_App app)
{
      // braid_Real send_recv_rank;
      // braid_Int  myid;

      // /* Grab information from the action */
      // send_recv_rank = action->send_recv_rank;
      // myid           = action->myid;

      // printf("BufUnack adjoint\n");
      return 0;
}
#endif