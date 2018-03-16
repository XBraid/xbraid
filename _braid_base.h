/**
 *  Header file defining the wrapper for the user routines
 **/


#ifndef _braid_base_HEADER
#define _braid_base_HEADER

#include "_braid.h"

/* TODO: Comment these routines! */


braid_Int 
_braid_BaseStep(braid_Core core, 
                braid_App        app,    
                braid_BaseVector     ustop,
                braid_BaseVector     fstop, 
                braid_BaseVector     u, 
                braid_StepStatus status );


                        
braid_Int
_braid_BaseInit(braid_Core core,
                braid_App      app, 
                braid_Real     t,   
                braid_BaseVector  *u_ptr
                );


braid_Int
_braid_BaseClone(braid_Core core,
                 braid_App app,  
                 braid_BaseVector   u,    
                 braid_BaseVector  *v_ptr 
                 );

braid_Int
_braid_BaseFree(braid_Core core,
                braid_App     app,
                braid_BaseVector  u   
                );


braid_Int
_braid_BaseSum(braid_Core core,
                braid_App        app,    
                braid_Real    alpha,  
                braid_BaseVector  x,      
                braid_Real    beta,   
                braid_BaseVector  y       
                );


braid_Int
_braid_BaseSpatialNorm(braid_Core core,
                       braid_App      app,      
                       braid_BaseVector   u,        
                       braid_Real    *norm_ptr  
                       );


braid_Int
_braid_BaseAccess(braid_Core core,
                  braid_App           app,   
                  braid_BaseVector        u,     
                  braid_AccessStatus  status 
                  );

braid_Int
_braid_BaseBufSize(braid_Core core,
                   braid_App   app,           
                   braid_Int  *size_ptr,      
                   braid_BufferStatus  status 
                   );      



braid_Int
_braid_BaseBufPack(braid_Core core,
                   braid_App           app,       
                   braid_BaseVector        u,         
                   void               *buffer,    
                   braid_BufferStatus  status     
                   );


braid_Int
_braid_BaseBufUnpack(braid_Core core,
                     braid_App            app,    
                     void                *buffer, 
                     braid_BaseVector        *u_ptr,  
                     braid_BufferStatus   status  
                     );


braid_Int
_braid_BaseResidual(braid_Core core,
                     braid_App        app,    /**< user-defined _braid_App structure */
                       braid_BaseVector     ustop,  /**< input, u vector at *tstop* */
                       braid_BaseVector     r     , /**< output, residual at *tstop* (at input, equals *u* at *tstart*) */
                       braid_StepStatus status  /**< query this struct for info about u (e.g., tstart and tstop) */ 
                       );

braid_Int
_braid_BaseSCoarsen(braid_Core core,
                    braid_App               app,    /**< user-defined _braid_App structure */
                       braid_BaseVector            fu,     /**< braid_BaseVector to refine*/
                       braid_BaseVector           *cu_ptr, /**< output, refined vector */   
                       braid_CoarsenRefStatus  status  /**< query this struct for info about fu and cu (e.g., where in time fu and cu are)  */ 
                       );


braid_Int
_braid_BaseSRefine(braid_Core core,
                   braid_App               app,    /**< user-defined _braid_App structure */
                      braid_BaseVector            cu,     /**< braid_BaseVector to refine*/
                      braid_BaseVector           *fu_ptr, /**< output, refined vector */       
                      braid_CoarsenRefStatus  status  /**< query this struct for info about fu and cu (e.g., where in time fu and cu are)  */ 
                      );

braid_Int
_braid_BaseSInit(braid_Core core,
                 braid_App     app,           /**< user-defined _braid_App structure */
                   braid_Real     t,             /**< time value for *u_ptr* */
                   braid_BaseVector  *u_ptr          /**< output, newly allocated and initialized vector shell */
                   );

braid_Int
_braid_BaseSClone(braid_Core core, 
                 braid_App      app,          /**< user-defined _braid_App structure */
                    braid_BaseVector   u,            /**< vector to clone */ 
                    braid_BaseVector  *v_ptr         /**< output, newly allocated and cloned vector shell */
                    );


braid_Int
_braid_BaseSFree(braid_Core core,
                  braid_App     app,            /**< user-defined _braid_App structure */
                    braid_BaseVector  u               /**< vector to free (keeping the shell) */
                    );

braid_Int
_braid_BaseTimeGrid(braid_Core core,
                   braid_App         app,       /**< user-defined _braid_App structure */
                       braid_Real       *ta,        /**< temporal grid on level 0 (slice per processor) */
                       braid_Int        *ilower,    /**< lower time index value for this processor */
                       braid_Int        *iupper     /**< upper time index value for this processor */
                       );


/*---- Differentiated user routines ---- */

braid_Int
_braid_BaseStep_diff(_braid_Action *action);


braid_Int
_braid_BaseClone_diff(_braid_Action *action);

braid_Int
_braid_BaseSum_diff(_braid_Action *action);

braid_Int
_braid_BaseAccess_diff(_braid_Action *action);

braid_Int
_braid_BaseBufPack_diff(_braid_Action *action, braid_App app);

braid_Int
_braid_BaseBufUnpack_diff(_braid_Action *action, braid_App app);

#endif