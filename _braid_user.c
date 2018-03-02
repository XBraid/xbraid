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
   /* Call the users BufUnPack function */
   _braid_CoreFcn(core, bufunpack)(app, buffer, u_ptr, status);

   return 0;
}


#endif