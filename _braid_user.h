/**
 *  Header file defining the wrapper for the user routines
 **/


#ifndef _braid_user_HEADER
#define _braid_user_HEADER

#include "_braid.h"

/* TODO: Comment these routines! */


braid_Int 
_braid_UserStep(braid_Core core, 
                braid_App        app,    
                braid_Vector     ustop,
                braid_Vector     fstop, 
                braid_Vector     u, 
                braid_StepStatus status );


                        
braid_Int
_braid_UserInit(braid_Core core,
                braid_App      app, 
                braid_Real     t,   
                braid_Vector  *u_ptr
                );


braid_Int
_braid_UserClone(braid_Core core,
                 braid_App app,  
                 braid_Vector   u,    
                 braid_Vector  *v_ptr 
                 );

braid_Int
_braid_UserFree(braid_Core core,
                braid_App     app,
                braid_Vector  u   
                );


braid_Int
_braid_UserSum(braid_Core core,
                braid_App        app,    
                braid_Real    alpha,  
                braid_Vector  x,      
                braid_Real    beta,   
                braid_Vector  y       
                );


braid_Int
_braid_UserSpatialNorm(braid_Core core,
                       braid_App      app,      
                       braid_Vector   u,        
                       braid_Real    *norm_ptr  
                       );


braid_Int
_braid_UserAccess(braid_Core core,
                  braid_App           app,   
                  braid_Vector        u,     
                  braid_AccessStatus  status 
                  );

braid_Int
_braid_UserBufSize(braid_Core core,
                   braid_App   app,           
                   braid_Int  *size_ptr,      
                   braid_BufferStatus  status 
                   );      



braid_Int
_braid_UserBufPack(braid_Core core,
                   braid_App           app,       
                   braid_Vector        u,         
                   void               *buffer,    
                   braid_BufferStatus  status     
                   );


braid_Int
_braid_UserBufUnpack(braid_Core core,
                     braid_App            app,    
                     void                *buffer, 
                     braid_Vector        *u_ptr,  
                     braid_BufferStatus   status  
                     );


/*---- Adjoint routines ---- */

braid_Int
_braid_UserStepAdjoint(_braid_Action *action, braid_App app);


braid_Int
_braid_UserCloneAdjoint(_braid_Action *action);

braid_Int
_braid_UserSumAdjoint(_braid_Action *action);

braid_Int
_braid_UserAccessAdjoint(_braid_Action *action, braid_App app);

braid_Int
_braid_UserBufPackAdjoint(_braid_Action *action, braid_App app);

braid_Int
_braid_UserBufUnpackAdjoint(_braid_Action *action, braid_App app);

#endif