/**
 *  Source file implementing the wrapper for the user routines
 **/


#ifndef _braid_base_HEADER
#define _braid_base_HEADER

#include "_braid.h"
#include "_braid_status.h"


braid_Int 
_braid_BaseStep(braid_Core       core,
                braid_App        app,    
                braid_BaseVector     ustop,
                braid_BaseVector     fstop, 
                braid_BaseVector     u, 
                braid_StepStatus status )
{

    if (_braid_CoreElt(core, verbose)) 
    {
       printf("%d: STEP pushes %p\n",_braid_CoreElt(core, myid), u->bar);
      // _braid_CoreFcn(core, access)( app, u->bar->userVector, NULL );
    }

   /* Record to the tape */
   if ( _braid_CoreElt(core, record) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = STEP;
      action->core          = core;
      action->inTime        = _braid_CoreElt(core, t);
      action->outTime       = _braid_CoreElt(core, tnext);
      action->status        = (braid_Status) status;
      action->myid          = _braid_CoreElt(core, myid);
      // printf("STEP push\n");
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

      /* Push a copy of the primal vector to the primal tape */
      braid_Vector u_copy;
      _braid_CoreFcn(core, clone)(app, u->userVector, &u_copy);  // this will accolate memory for the copy!
      _braid_CoreElt(core, userVectorTape) = _braid_TapePush( _braid_CoreElt(core, userVectorTape), u_copy);

      /* Push a copy of the bar vector to the bar tape */
      braid_VectorBar bar_copy;
      _braid_VectorBarCopy(u->bar, &bar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), bar_copy);

      /* DEBUG: Display the vector */
      // printf("STEP pushes userVector: ");
      // braid_AccessStatus astatus = (braid_AccessStatus) status;
      // _braid_CoreFcn(core, access)( app, u->userVector, astatus);
  }

   /* Call the users Step function */
   if ( _braid_CoreElt(core, residual) == NULL )
   {
      _braid_CoreFcn(core, step)(app, ustop->userVector, NULL, u->userVector, status);
   }
   else
   {
      /* TODO: Check the fstop feature!!! */
      _braid_CoreFcn(core, step)(app, ustop->userVector, fstop->userVector, u->userVector, status);
   }
   return 0;
}


                        
braid_Int
_braid_BaseInit(braid_Core core,
                braid_App      app, 
                braid_Real     t,   
                braid_BaseVector  *u_ptr
                )
{
   if (_braid_CoreElt(core, verbose)) printf("%d INIT\n", _braid_CoreElt(core, myid));

   /* Allocate memory for the braid_BaseVector */
   braid_BaseVector u = (braid_BaseVector)malloc(sizeof(braid_Vector)+sizeof(braid_VectorBar));
   u->userVector= NULL;
   u->bar = NULL;

   /* Call the users init function. This allocates and initializes the primal */
   _braid_CoreFcn(core, init)(app, t, &(u->userVector));
   
   /* Allocate and initialize the bar */
   if ( _braid_CoreElt(core, adjoint) )
   {
      braid_VectorBar mybar = (braid_VectorBar)malloc(sizeof(braid_Vector)+sizeof(int));
      mybar->useCount = 1;
      _braid_CoreFcn(core, init)(app, t, &(mybar->userVector));
      _braid_CoreFcn(core, sum)(app, -1.0, mybar->userVector, 1.0, mybar->userVector);
      u->bar = mybar;

      // _braid_CoreFcn(core, access)( app, u->bar->userVector, NULL );
   }

   /* Record to the tape */
   if ( _braid_CoreElt(core, record) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = INIT;
      action->core          = core;
      action->inTime        = _braid_CoreElt(core, t);
      action->myid          = _braid_CoreElt(core, myid);
      // printf("INIT push: \n");
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);


      /* DEBUG: Display the vector */
      // printf("Init creates %p\n", u->bar);
      // _braid_CoreFcn(core, access)( app, u->bar->userVector, NULL );
   
   }

  /* Debug: */
//    printf("Init bar, useCount: %d, bar value: ", u->bar->useCount);
//   braid_AccessStatus astatus = (braid_AccessStatus) core;
//   _braid_CoreFcn(core, access)(app, u->bar->userVector, astatus);

   /* Set the return pointer */
   *u_ptr = u;

   return 0;
}

braid_Int
_braid_BaseClone(braid_Core core,
                 braid_App app,  
                 braid_BaseVector   u,    
                 braid_BaseVector  *v_ptr 
                 )
{

   if (_braid_CoreElt(core, verbose)) printf("%d: CLONE\n",_braid_CoreElt(core, myid));

   /* Allocate memory for the braid_BaseVector */
   braid_BaseVector v = (braid_BaseVector)malloc(sizeof(braid_Vector)+sizeof(braid_VectorBar));
   v->userVector  = NULL;
   v->bar = NULL;

   /* Call the users Clone function for the primal braid_Vector */
   _braid_CoreFcn(core, clone)(app, u->userVector, &(v->userVector) );

   if (_braid_CoreElt(core, adjoint))
   {
      /* Allocate and initialize the bar to zero*/
      braid_VectorBar mybar = (braid_VectorBar)malloc(sizeof(braid_Vector)+sizeof(int));
      mybar->useCount = 1;
      _braid_CoreFcn(core, clone)(app, u->bar->userVector, &(mybar->userVector));
      _braid_CoreFcn(core, sum)(app, -1.0, mybar->userVector, 1.0, mybar->userVector);
      v->bar = mybar;
   } 


   /* Record to the tape */
   if ( _braid_CoreElt(core, record) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = CLONE;
      action->core          = core;
      action->myid          = _braid_CoreElt(core, myid);
      // printf("CLONE push: \n");
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

      /* Copy and push both bar pointers to the bar tape */
      braid_VectorBar bar_copy_u;
      braid_VectorBar bar_copy_v;
      _braid_VectorBarCopy(u->bar, &bar_copy_u);
      _braid_VectorBarCopy(v->bar, &bar_copy_v);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), bar_copy_u);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), bar_copy_v);
      
      /*Debug: */
      // printf("Clone pushes old u %p\n", u->bar);
      // printf("Clone creates new v %p\n", v->bar);

   }


   *v_ptr = v;

   return 0;
}


braid_Int
_braid_BaseFree(braid_Core core,
                braid_App     app,
                braid_BaseVector  u   
                )
{

   if (_braid_CoreElt(core, verbose)) printf("%d: FREE\n",_braid_CoreElt(core, myid));
   /*Debug: */
//    _braid_CoreFcn(core, access)( app, u->bar->userVector, NULL );

   /* if bar: Record to the tape */
   if ( _braid_CoreElt(core, record) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = FREE;
      action->core          = core;
      action->myid          = _braid_CoreElt(core, myid);
      // printf("FREE push\n");
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);
   }
 
   /* Call the users free function for the userVector */
   _braid_CoreFcn(core, free)(app, u->userVector);


   /* Debug: */
//    printf("Init bar, useCount: %d, bar value: ", u->bar->useCount);
//    braid_AccessStatus astatus = (braid_AccessStatus) core;
//    _braid_CoreFcn(core, access)(app, u->bar->userVector, astatus);

   if ( _braid_CoreElt(core, adjoint) )
   {
      /* Decrease usecount of the bar and free the bar if useCount == 0. */
      // printf("Free %p\n", u->bar);
      // printf(", useCount: %d \n", u->bar->useCount);
      _braid_VectorBarDelete(core, u->bar);
   }

   /* Free the braid_BaseVector */
   free(u);

   return 0;
}


braid_Int
_braid_BaseSum(braid_Core core,
               braid_App        app,    
               braid_Real    alpha,  
               braid_BaseVector  x,      
               braid_Real    beta,   
               braid_BaseVector  y       
               )
{
   if (_braid_CoreElt(core, verbose)) printf("%d: SUM\n",_braid_CoreElt(core, myid));

   /* if bar: Record to the tape */
   if ( _braid_CoreElt(core, record) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = SUM;
      action->core          = core;
      action->sum_alpha     = alpha;
      action->sum_beta      = beta;
      action->myid          = _braid_CoreElt(core, myid);
      // printf("SUM push\n");
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

      /* Copy and push both bar pointers to the bar tape */
      braid_VectorBar bar_copy_x;
      braid_VectorBar bar_copy_y;
      _braid_VectorBarCopy(x->bar, &bar_copy_x);
      _braid_VectorBarCopy(y->bar, &bar_copy_y);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), bar_copy_x);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), bar_copy_y);
   }

    /* Call the users Sum function */
   _braid_CoreFcn(core, sum)(app, alpha, x->userVector, beta, y->userVector);

   
   /*Debug: */
      // printf("Sum pushes y %p\n", y->bar);
      // printf("Sum pushes x %p\n", x->bar);
   
   return 0;
}


braid_Int
_braid_BaseSpatialNorm(braid_Core core,
                       braid_App      app,      /**< user-defined _braid_App structure */
                       braid_BaseVector   u,        /**< vector to norm */
                       braid_Real    *norm_ptr  /**< output, norm of braid_BaseVector (this is a spatial norm) */ 
                       )
{

   /* Call the users SpatialNorm function */
   _braid_CoreFcn(core, spatialnorm)(app, u->userVector, norm_ptr);

   return 0;
}


braid_Int
_braid_BaseAccess(braid_Core core,
                  braid_App           app,   
                  braid_BaseVector        u,     
                  braid_AccessStatus  status 
                  )
{
   if (_braid_CoreElt(core, verbose)) printf("%d: ACCESS\n",_braid_CoreElt(core, myid));

   /* if bar: Record to the tape */
   if ( _braid_CoreElt(core, record) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = ACCESS;
      action->core          = core;
      action->status        = (braid_Status) status;
      action->inTime        = _braid_CoreElt(core, t);
      action->myid          = _braid_CoreElt(core, myid);
      // printf("ACCESS push\n");
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

   }

   /* Call the users Access function */
   _braid_CoreFcn(core, access)(app, u->userVector, status);

   /* Debug */
//    _braid_CoreFcn(core, access)( app, u->bar->userVector, NULL );

   return 0;
}


braid_Int
_braid_BaseBufSize(braid_Core core,
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
_braid_BaseBufPack(braid_Core core,
                   braid_App           app,       
                   braid_BaseVector        u,         
                   void               *buffer,    
                   braid_BufferStatus  status     
                   )
{
   if (_braid_CoreElt(core, verbose)) printf("%d: BUFPACK\n", _braid_CoreElt(core, myid) );

      /* if bar: Record to the tape */
   if ( _braid_CoreElt(core, record) )
   {
      /* Set up and push the action */
      _braid_Action* action  = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall      = BUFPACK;
      action->core          = core;
      action->send_recv_rank = _braid_CoreElt(core, send_recv_rank);
      action->myid           = _braid_CoreElt(core, myid);
      // printf("BufPack push\n");
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

      /* Copy and push the bar pointer to the bar tape */
      braid_VectorBar bar_copy;
      _braid_VectorBarCopy(u->bar, &bar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), bar_copy);

   }
   
   /* Call the users BufPack function */
   _braid_CoreFcn(core, bufpack)(app, u->userVector, buffer, status);

   return 0;
}


braid_Int
_braid_BaseBufUnpack(braid_Core core,
                     braid_App            app,    
                     void                *buffer, 
                     braid_BaseVector        *u_ptr,  
                     braid_BufferStatus   status  
                     )
{
   if (_braid_CoreElt(core, verbose)) printf("%d: BUFUNPACK\n", _braid_CoreElt(core, myid));

   /* Allocate memory for the braid_BaseVector */
   braid_BaseVector u = (braid_BaseVector)malloc(sizeof(braid_Vector)+sizeof(braid_VectorBar));
   u->userVector  = NULL;
   u->bar = NULL;
   

   /* Call the users BufUnPack function for the userVector */
   _braid_CoreFcn(core, bufunpack)(app, buffer, &(u->userVector), status);

   if ( _braid_CoreElt(core, adjoint) )
   {
      /* Allocate and initialize the bar with zero*/
      braid_VectorBar mybar = (braid_VectorBar)malloc(sizeof(braid_Vector)+sizeof(int));
      mybar->useCount = 1;
      _braid_CoreFcn(core, init)(app, _braid_CoreElt(core, tstart), &(mybar->userVector));
      _braid_CoreFcn(core, sum)(app, -1.0, mybar->userVector, 1.0, mybar->userVector);
      u->bar = mybar;
   }

   /* Record to the tape */
   if ( _braid_CoreElt(core, record) )
   {
      /* Set up and push the action */
      _braid_Action* action  = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall      = BUFUNPACK;
      action->core          = core;
      action->send_recv_rank = _braid_CoreElt(core, send_recv_rank);
      action->myid           = _braid_CoreElt(core, myid);
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

      /* Copy and push the bar pointer to the bar tape */
      braid_VectorBar bar_copy;
      _braid_VectorBarCopy(u->bar, &bar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), bar_copy);
    }

  /* Debug: */
//    printf("Init bar, useCount: %d, bar value: ", u->bar->useCount);
//   braid_AccessStatus astatus = (braid_AccessStatus) core;
//   _braid_CoreFcn(core, access)(app, u->bar->userVector, astatus);

  
   *u_ptr = u;

   return 0;
}


braid_Int
_braid_BaseObjectiveT(braid_Core core,
                      braid_App  app,
                      braid_BaseVector u,
                      braid_AccessStatus astatus,
                      braid_Real *objT_ptr
                      )
{
   if (_braid_CoreElt(core, verbose)) printf("%d: OBJECTIVET\n", _braid_CoreElt(core, myid));

   /* if bar: Record to the tape */
   if ( _braid_CoreElt(core, record) )
   {
      /* Set up and push the action */
      _braid_Action* action = _braid_CTAlloc(_braid_Action, 1);
      action->braidCall     = OBJECTIVET;
      action->core          = core;
      action->status        = (braid_Status) astatus;
      action->inTime        = _braid_CoreElt(core, t);
      action->myid          = _braid_CoreElt(core, myid);
      // printf("ACCESS push\n");
      _braid_CoreElt(core, actionTape) = _braid_TapePush( _braid_CoreElt(core, actionTape) , action);

      /* Push a copy of the primal vector to the primal tape */
      braid_Vector u_copy;
      _braid_CoreFcn(core, clone)(app, u->userVector, &u_copy);  // this will accolate memory for the copy!
      _braid_CoreElt(core, userVectorTape) = _braid_TapePush( _braid_CoreElt(core, userVectorTape), u_copy);

      /* Push a copy of the bar vector to the bar tape */
      braid_VectorBar bar_copy;
      _braid_VectorBarCopy(u->bar, &bar_copy);
      _braid_CoreElt(core, barTape) = _braid_TapePush(_braid_CoreElt(core, barTape), bar_copy);

      /* Debug info: */
      // printf("OBJT pushes userVector: ");
      // _braid_CoreFcn(core, access)(app, u->userVector, (braid_AccessStatus) status);

   }

   /* Call the users objective function */
   _braid_CoreFcn(core, objectiveT)(app, u->userVector, astatus, objT_ptr);


   /* Debug */
      // printf("ObjT pushes %p\n", u->bar);

   return 0;
}


braid_Int
_braid_BaseResidual(braid_Core core,
                     braid_App        app,    /**< user-defined _braid_App structure */
                       braid_BaseVector     ustop,  /**< input, u vector at *tstop* */
                       braid_BaseVector     r     , /**< output, residual at *tstop* (at input, equals *u* at *tstart*) */
                       braid_StepStatus status  /**< query this struct for info about u (e.g., tstart and tstop) */ 
                       )
{
   if (_braid_CoreElt(core, verbose)) printf("RESIDUAL\n");

      /* Call the users Residual function */
   _braid_CoreFcn(core, residual)(app, ustop->userVector, r->userVector, status);

   return 0;
}


braid_Int
_braid_BaseSCoarsen(braid_Core core,
                    braid_App               app,    /**< user-defined _braid_App structure */
                       braid_BaseVector            fu,     /**< braid_Vector to refine*/
                       braid_BaseVector           *cu_ptr, /**< output, refined vector */   
                       braid_CoarsenRefStatus  status  /**< query this struct for info about fu and cu (e.g., where in time fu and cu are)  */ 
                       )
{
      braid_BaseVector cu = (braid_BaseVector)malloc(sizeof(braid_BaseVector));

      /* Call the users SCoarsen Function */
      _braid_CoreFcn(core, scoarsen)(app, fu->userVector, &(cu->userVector), status);

      *cu_ptr = cu;

   return 0;
}

braid_Int
_braid_BaseSRefine(braid_Core core,
                   braid_App               app,    /**< user-defined _braid_App structure */
                      braid_BaseVector            cu,     /**< braid_Vector to refine*/
                      braid_BaseVector           *fu_ptr, /**< output, refined vector */       
                      braid_CoarsenRefStatus  status  /**< query this struct for info about fu and cu (e.g., where in time fu and cu are)  */ 
                      )
{
      braid_BaseVector fu = (braid_BaseVector)malloc(sizeof(braid_BaseVector));

      /* Call the users SRefine */
      _braid_CoreFcn(core, srefine)(app, cu->userVector, &(fu->userVector), status);

      *fu_ptr = fu;

   return 0;
}                      

braid_Int
_braid_BaseSInit(braid_Core core,
                 braid_App     app,           /**< user-defined _braid_App structure */
                   braid_Real     t,             /**< time value for *u_ptr* */
                   braid_BaseVector  *u_ptr          /**< output, newly allocated and initialized vector shell */
                   )
{
      braid_BaseVector u = (braid_BaseVector)malloc(sizeof(braid_BaseVector));

      /* Call the users SInit */
      _braid_CoreFcn(core, sinit)(app, t, &(u->userVector));

      *u_ptr = u;

   return 0;
}

braid_Int
_braid_BaseSClone(braid_Core core, 
                 braid_App      app,          /**< user-defined _braid_App structure */
                    braid_BaseVector   u,            /**< vector to clone */ 
                    braid_BaseVector  *v_ptr         /**< output, newly allocated and cloned vector shell */
                    )
{
      braid_BaseVector v = (braid_BaseVector)malloc(sizeof(braid_BaseVector));

      /* Call the users SClone */
      _braid_CoreFcn(core, sclone)(app, u->userVector, &(v->userVector));

      *v_ptr = v;

   return 0;
}


braid_Int
_braid_BaseSFree(braid_Core core,
                  braid_App     app,            /**< user-defined _braid_App structure */
                    braid_BaseVector  u               /**< vector to free (keeping the shell) */
                    )
{
      /* Call the users sfree */
      _braid_CoreFcn(core, sfree)(app, u->userVector);

      free(u);

   return 0;
}

braid_Int
_braid_BaseTimeGrid(braid_Core core,
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


/*----- Differentiated user routines ------*/

braid_Int
_braid_BaseStep_diff(_braid_Action *action)
{

      braid_Real        inTime;
      braid_Real        outTime;
      braid_StepStatus  status;
      braid_Core        core;
      // braid_Int         myid;

      /* Grab information from the action */
      core    = action->core;
      status  = (braid_StepStatus) action->status;
      inTime  = action->inTime;
      outTime = action->outTime;
      // myid    = action->myid;

      if (_braid_CoreElt(core, verbose)) printf("%d: STEP_DIFF\n", _braid_CoreElt(core, myid));

      /* Get the braid_vector that was used in primal run */
      braid_Vector userVector;
      userVector = (braid_Vector) (_braid_CoreElt(core, userVectorTape)->data_ptr);


      /* Pop the bar vector from the tape*/
      braid_VectorBar      bar;
      bar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
      _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );

      /* DEBUG: Display the vector */
      // printf("STEP bar pops from the primal tape\n");
      // _braid_CoreFcn(core, access)(_braid_CoreElt(core, app), userVector, astatus );
      // printf("Step pops %p\n", bar);

      /* Trigger the status with the correct tstart and tstop*/
      _braid_StatusElt(status, t)     = inTime;
      _braid_StatusElt(status, tnext) = outTime;

      /* Call the users's differentiated step function */
      _braid_CoreFcn(core, step_diff)(_braid_CoreElt(core, app), userVector, bar->userVector, status);


      /* Decrease the useCount of the bar and pop the bar from the bar tape */
      _braid_VectorBarDelete(core, bar);

      /* Delete memory of the primal vector and pop it from the primal tape*/
      _braid_CoreFcn(core, free)(_braid_CoreElt(core, app), userVector);
      _braid_CoreElt(core, userVectorTape) = _braid_TapePop( _braid_CoreElt(core, userVectorTape) );

      return 0;
}

braid_Int
_braid_BaseClone_diff(_braid_Action *action)
{

      /* Grab information from the action */
      braid_Core core = action->core;
      // braid_Int         myid;
      // myid = action->myid;
      
      if (_braid_CoreElt(core, verbose)) printf("%d: CLONE_DIFF\n", _braid_CoreElt(core, myid));

      /* Pop the bar vectors from the tape */
      braid_VectorBar v_bar;
      braid_VectorBar u_bar;

      v_bar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
      _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );

      u_bar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
      _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );

      // /* DEBUG: Display the vector */
      // printf("Clone pops v %p\n", v_bar);
      // printf("Clone pops u %p\n", u_bar);


      /* Perform the bar action :
      *  ub += vb
      *  vb  = 0.0
      */
      _braid_CoreFcn(core, sum)(_braid_CoreElt(core,app),1.0, v_bar->userVector,  1.0, u_bar->userVector);
      _braid_CoreFcn(core, sum)(_braid_CoreElt(core,app),1.0, v_bar->userVector, -1.0, v_bar->userVector);


      /* Decrease the useCount of the bar */
      _braid_VectorBarDelete(core, u_bar);
      _braid_VectorBarDelete(core, v_bar);

 



      return 0;
}

braid_Int
_braid_BaseSum_diff(_braid_Action *action)
{
   braid_Core core;
   braid_Real alpha;
   braid_Real beta;
   // braid_Int  myid;

   /* Grab information from the action */
   core  = action->core;
   alpha = action->sum_alpha;
   beta  = action->sum_beta;
   // myid  = action->myid;

   if (_braid_CoreElt(core, verbose)) printf("%d: SUM_DIFF\n", _braid_CoreElt(core, myid));

   /* Pop the bar vectors from the tape */
   braid_VectorBar y_bar;
   braid_VectorBar x_bar;

   y_bar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
   _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );

   x_bar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
   _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );

   /*DEBUG*/
      // printf("SUM pops y %p\n", y_bar);
      // printf("SUM pops x %p\n", x_bar);

   /* Perform the bar action: 
   *  xb += alpha * yb
   *  yb  = beta  * yb
   */
   _braid_CoreFcn(core, sum)(_braid_CoreElt(core, app), alpha, y_bar->userVector,  1.0, x_bar->userVector);
   _braid_CoreFcn(core, sum)(_braid_CoreElt(core, app),   0.0, x_bar->userVector, beta, y_bar->userVector);

   /* Decrease the useCount of the bar */
   _braid_VectorBarDelete(core, y_bar);
   _braid_VectorBarDelete(core, x_bar);


   return 0;
}


braid_Int
_braid_BaseObjectiveT_diff(_braid_Action *action)
{
//       braid_Int     inTime;
//       braid_Int     myid;

      /* Grab information from the app */
      braid_Core core            = action->core;
      braid_AccessStatus astatus = (braid_AccessStatus) action->status;
//       inTime = action->inTime;
//       myid   = action->myid;
 
      if (_braid_CoreElt(core, verbose)) printf("%d: OBJT_DIFF\n", _braid_CoreElt(core, myid));

      /* Pop the primal vector that was used in primal access function */
      braid_Vector  userVector;
      userVector = (braid_Vector) (_braid_CoreElt(core, userVectorTape)->data_ptr);
      _braid_CoreElt(core, userVectorTape) = _braid_TapePop( _braid_CoreElt(core, userVectorTape) );

      /* Pop the bar vector from the tape*/
      braid_VectorBar      bar;
      bar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
      _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );


      // printf("ObjT pops u %p\n", bar);
      
     /* Call the users's differentiated objective function */
      braid_Real f_bar = _braid_CoreElt(core, optim)->f_bar;
      _braid_CoreFcn(core, objT_diff)(_braid_CoreElt(core, app), userVector, bar->userVector, f_bar, astatus);

      /* Free memory of the primal Vector and pop it from the tape */
      _braid_CoreFcn(core, free)(_braid_CoreElt(core, app), userVector);

      /* Decrease the useCount of the bar and pop the bar from the bar tape */
      _braid_VectorBarDelete(core, bar);


      return 0;
}

braid_Int
_braid_BaseBufPack_diff(_braid_Action *action, braid_App app)
{
      braid_Core core;
      braid_Real send_recv_rank;
      braid_Int  size;
      void       *buffer;
      braid_Vector u;
      braid_BufferStatus bstatus;
      braid_VectorBar      bar;

      /* Grab information from the action */
      core           = action->core;
      send_recv_rank = action->send_recv_rank;

      if (_braid_CoreElt(core, verbose)) printf("%d: BUFPACK_DIFF\n", _braid_CoreElt(core, myid));

      /* Pop the bar vector from the tape*/
      bar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
      _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );


      /* Allocate buffer through the user routine */
      bstatus = (braid_BufferStatus) core;
      _braid_BufferStatusInit( 0, 0, bstatus);
      _braid_CoreFcn(core, bufsize)(_braid_CoreElt(core, app), &size, bstatus);
      buffer = malloc(size);

      /* Receive the buffer 
       * TODO: Add CommHandle / status check!!
       * WHY blocking Recv ? 
       */
      // MPI_Request *requests;
      // requests = _braid_CTAlloc(MPI_Request, 1);
      // MPI_Irecv(buffer, size, MPI_BYTE, send_recv_rank, 0, _braid_CoreElt(core, comm), &requests[0]); 
      // MPI_Request_free(requests);
      MPI_Recv(buffer, size, MPI_BYTE, send_recv_rank, 0, _braid_CoreElt(core, comm), MPI_STATUS_IGNORE); 

      /* Unpack the buffer (this allocates memory for a braid_Vector */
      _braid_CoreFcn(core, bufunpack)(_braid_CoreElt(core, app), buffer, &u, bstatus);

      /* Update the bar  */
      _braid_CoreFcn(core, sum)(_braid_CoreElt(core, app), 1., u, 1., bar->userVector);

      /*DEBUG */
      // printf("%d: bufunpack_diff recvs ", _braid_CoreElt(core, myid));
      // _braid_CoreFcn(core, access)(app, bar->userVector, NULL);      

      /* Decrease the useCount of the bar and pop the bar from the bar tape */
      _braid_VectorBarDelete(core, bar);

      /* Free memory */
      free(buffer);
      _braid_CoreFcn(core, free)(_braid_CoreElt(core, app), u);

      return 0;
}

braid_Int
_braid_BaseBufUnpack_diff(_braid_Action *action, braid_App app)
{
      braid_Real         send_recv_rank;
      braid_VectorBar      bar;
      braid_Core         core;
      braid_BufferStatus bstatus;
      void               *buffer;
      braid_Int          size;
      MPI_Request        *requests;

      /* Grab information from the action */
      core           = action->core;
      send_recv_rank = action->send_recv_rank;

      if (_braid_CoreElt(core, verbose)) printf("%d: BUFUNPACK_DIFF\n", _braid_CoreElt(core, myid));

      /* Pop the bar vector from the tape*/
      bar = (braid_VectorBar) (_braid_CoreElt(core, barTape)->data_ptr);
      _braid_CoreElt(core, barTape) = _braid_TapePop( _braid_CoreElt(core, barTape) );

     /* Allocate buffer through the user routine */
      bstatus = (braid_BufferStatus) core;
      _braid_BufferStatusInit( 0, 0, bstatus);
      _braid_CoreFcn(core, bufsize)(_braid_CoreElt(core, app), &size, bstatus);
      buffer = malloc(size); 

      /* Pack the buffer with the bar variable */
      _braid_CoreFcn(core, bufpack)(_braid_CoreElt(core, app), bar->userVector, buffer, bstatus);

      /*DEBUG */
      // printf("%d: bufunpack_diff isends ", _braid_CoreElt(core, myid));
      // _braid_CoreFcn(core, access)(app, bar->userVector, NULL);      

      /* Send the buffer  */
      requests = _braid_CTAlloc(MPI_Request, 1);
      MPI_Isend(buffer, size, MPI_BYTE, send_recv_rank, 0, _braid_CoreElt(core, comm), &requests[0]);
      MPI_Request_free(requests);

      
      /* Set ubar to zero */
      _braid_CoreFcn(core, sum)(app, -1., bar->userVector, 1., bar->userVector );


      /* Decrease the useCount of the bar and pop the bar from the bar tape */
      _braid_VectorBarDelete(core, bar);


      return 0;
}
#endif