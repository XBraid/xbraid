/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
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
 

/** \file _braid_F90_iface.c
 * \brief Define F90 interface for user routines.
 *
 * This file contains definitions to handle type changes between F90 and C, name 
 * mangling and stubs that wrap the user-written F90 routines.  This
 * layer of wrapping is required so that we can handle any type changes (e.g.,
 * communicators) between F90 and C.  And, we also need to pass function pointers, 
 * which wrap the user-written F90, to XBraid.  Note F90 does not have function 
 * pointers.
 */

#ifndef braid_F90_Iface
#define braid_F90_Iface

#include "stdio.h"
#include "braid.h"
#include "_braid.h"
#include "braid_test.h"

/* Do you have any spots where you need to map types to Fortran? */

/*--------------------------------------------------------------------------
 * Define Fortran types 
 *--------------------------------------------------------------------------*/
typedef braid_Int      braid_F90_Comm;
typedef braid_Int     *braid_F90_ObjPtr;
typedef braid_Int      braid_F90_Int;
typedef braid_Int      braid_F90_Void;
typedef braid_Real     braid_F90_Real;

/*--------------------------------------------------------------------------
 * Define macros to map types going to a Fortran routine
 *--------------------------------------------------------------------------*/

#define braid_PassF90_Comm(arg)          (MPI_Comm_c2f(*arg))
#define braid_PassF90_Int(arg)           ((braid_Int *) &arg)
#define braid_PassF90_IntPtr(arg)        ((braid_Int *) arg)
#define braid_PassF90_VoidPtr(arg)       ((braid_Int *) arg)
#define braid_PassF90_Real(arg)          ((braid_Real *) &arg)
#define braid_PassF90_RealPtr(arg)       ((braid_Real *) arg)
#define braid_PassF90_Obj(arg)           ((braid_Int *) arg)
#define braid_PassF90_ObjRef(arg)        ((braid_Int *) arg)

/*--------------------------------------------------------------------------
 * Define macros to map types coming from a Fortran routine
 *--------------------------------------------------------------------------*/

#define braid_TakeF90_Comm(arg)          (MPI_Comm_f2c(*arg))
#define braid_TakeF90_Int(arg)           ((braid_Int) *arg)
#define braid_TakeF90_IntPtr(arg)        ((braid_Int *) arg)
#define braid_TakeF90_Real(arg)          ((braid_Real) *arg)
#define braid_TakeF90_RealPtr(arg)       ((braid_Real *) arg)
#define braid_TakeF90_Obj(obj, arg)      ((obj) arg)
#define braid_TakeF90_ObjDeref(obj, arg) ((obj) *arg)
#define braid_TakeF90_ObjPtr(obj, arg)   ((obj *) arg)


/*--------------------------------------------------------------------------
 * Define name mangling macro based on pound-defined braid_FMANGLE
 * braid_FMANGLE == 1 --> append one underscore
 * braid_FMANGLE == 2 --> no change 
 * braid_FMANGLE == 3 --> append two underscore
 * braid_FMANGLE == 4 --> capitalize
 *--------------------------------------------------------------------------*/

#define braid_F90_Name_1(name,Name) name##_
#define braid_F90_Name_2(name,Name) name
#define braid_F90_Name_3(name,Name) name##__
#define braid_F90_Name_4(name,Name) Name

#if   (braid_FMANGLE == 1)
#define braid_F90_Name(name,Name) braid_F90_Name_1(name,Name)
#elif (braid_FMANGLE == 2)
#define braid_F90_Name(name,Name) braid_F90_Name_2(name,Name)
#elif (braid_FMANGLE == 3)
#define braid_F90_Name(name,Name) braid_F90_Name_3(name,Name)
#elif (braid_FMANGLE == 4)
#define braid_F90_Name(name,Name) braid_F90_Name_4(name,Name)
#else
#define braid_F90_Name(name,Name) braid_F90_Name_2(name,Name)
#endif


/*--------------------------------------------------------------------------
 * Define the prototypes for the user-written writtens in F90.
 * Define the XBraid wrappers that call these prototypes. 
 * 
 * These wrappers are needed so that we can have C-style function pointers 
 * to pass to braid_Init, braid_TestInitAccess, and so on. 
 *
 * These functions are C calling Fortran.
 *--------------------------------------------------------------------------*/

/**
 * braid_Init_Vector
 *
 * Fortran interface, first we define the prototype for the user-defined function and then we
 * provide the C-wrapper around the user-written Fortran function
 * */
void braid_F90_Name(braid_init_vec_f90, BRAID_INIT_VEC_F90)(braid_F90_ObjPtr, braid_F90_Real *, braid_F90_ObjPtr);
braid_Int
braid_Init_Vec_F90_Iface(braid_App      app,           /**< user-defined _braid_App structure */
                         braid_Real     t,             /**< time value for *u_ptr* */
                         braid_Vector  *u_ptr          /**< output, newly allocated and initialized vector */
                         )
{
   /* Temporary scalar so that the calling function's alpha and beta are not overwritten */
   braid_Real t2 = t;
   braid_F90_Name(braid_init_vec_f90, BRAID_INIT_VEC_F90)( braid_PassF90_Obj(app), braid_PassF90_Real(t2), braid_PassF90_ObjRef(u_ptr) );
   return 0;
}


/**
 * braid_Access
 *
 * Fortran interface, first we define the prototype for the user-defined function and then we
 * provide the C-wrapper around the user-written Fortran function
 * */
void braid_F90_Name(braid_access_f90, BRAID_ACCESS_F90)(braid_F90_ObjPtr, braid_F90_ObjPtr, braid_F90_ObjPtr);
braid_Int
braid_Access_F90_Iface(braid_App           app,              /**< user-defined _braid_App structure */
                       braid_Vector        u,                /**< vector to be accessed */
                       braid_AccessStatus  status            /**< can be querried for info like the current XBraid Iteration */
                       )
{
   braid_F90_Name(braid_access_f90, BRAID_ACCESS_F90)( braid_PassF90_Obj(app), braid_PassF90_Obj(u), braid_PassF90_Obj(status));
   return 0;
}


/**
 * braid_Free
 *
 * Fortran interface, first we define the prototype for the user-defined function and then we
 * provide the C-wrapper around the user-written Fortran function
 * */
void braid_F90_Name(braid_free_f90, BRAID_FREE_F90)(braid_F90_ObjPtr, braid_F90_ObjPtr);
braid_Int
braid_Free_F90_Iface(braid_App     app,            /**< user-defined _braid_App structure */
                     braid_Vector  u               /**< vector to free */
                     )
{
   /* Pass the address of u.  Fortran must have a pointer here so that it can deallocate u */
   braid_F90_Name(braid_free_f90, BRAID_FREE_F90)( braid_PassF90_Obj(    app),  
                                                   braid_PassF90_ObjRef( &u) );
   return 0;
}


/**
 * braid_Clone
 *
 * Fortran interface, first we define the prototype for the user-defined function and then we
 * provide the C-wrapper around the user-written Fortran function
 * */
void braid_F90_Name(braid_clone_f90, BRAID_CLONE_F90)(braid_F90_ObjPtr, braid_F90_ObjPtr,  braid_F90_ObjPtr);
braid_Int
braid_Clone_F90_Iface(braid_App      app,          /**< user-defined _braid_App structure */
                      braid_Vector   u,            /**< vector to clone */ 
                      braid_Vector  *v_ptr         /**< output, newly allocated and cloned vector */
                      )
{
   braid_F90_Name(braid_clone_f90, BRAID_CLONE_F90)( braid_PassF90_Obj(    app),  
                                                     braid_PassF90_Obj(    u),
                                                     braid_PassF90_ObjRef( v_ptr) );
   return 0;
}


/**
 * braid_Sum
 *
 * Fortran interface, first we define the prototype for the user-defined function and then we
 * provide the C-wrapper around the user-written Fortran function
 * */
void braid_F90_Name(braid_sum_f90, BRAID_SUM_F90)(braid_F90_ObjPtr, braid_F90_Real *, braid_F90_ObjPtr,  braid_F90_Real *, braid_F90_ObjPtr);
braid_Int
braid_Sum_F90_Iface(braid_App     app,             /**< user-defined _braid_App structure */
                    braid_Real    alpha,           /**< scalar for AXPY */
                    braid_Vector  x,               /**< vector for AXPY */
                    braid_Real    beta,            /**< scalar for AXPY */
                    braid_Vector  y                /**< output and vector for AXPY */
                    )
{
   /* Temporary scalars so that the calling function's alpha and beta are not overwritten */
   braid_Real alpha2 = alpha;
   braid_Real beta2 = beta;
   braid_F90_Name(braid_sum_f90, BRAID_SUM_F90)( braid_PassF90_Obj(    app),  
                                                 braid_PassF90_Real(   alpha2),
                                                 braid_PassF90_Obj(    x),
                                                 braid_PassF90_Real(   beta2),
                                                 braid_PassF90_Obj(    y) );
   return 0;
}


/**
 * braid_SpatialNorm
 *
 * Fortran interface, first we define the prototype for the user-defined function and then we
 * provide the C-wrapper around the user-written Fortran function
 * */
void braid_F90_Name(braid_spatialnorm_f90, BRAID_SPATIALNORM_F90)(braid_F90_ObjPtr, braid_F90_ObjPtr,  braid_F90_Real *);
braid_Int
braid_SpatialNorm_F90_Iface(braid_App      app,                /**< user-defined _braid_App structure */
                            braid_Vector   u,                  /**< vector to norm */
                            braid_Real    *norm_ptr            /**< output, norm of braid_Vector (this is a spatial norm) */ 
                            )
{
   braid_F90_Name(braid_spatialnorm_f90, BRAID_SPATIALNORM_F90)( 
                                  braid_PassF90_Obj(     app),
                                  braid_PassF90_Obj(     u),
                                  braid_PassF90_RealPtr( norm_ptr) );
   return 0;
}


/**
 * braid_Step
 *
 * Fortran interface, first we define the prototype for the user-defined function and then we
 * provide the C-wrapper around the user-written Fortran function
 * */
void braid_F90_Name(braid_step_f90, BRAID_STEP_F90)(braid_F90_ObjPtr, braid_F90_ObjPtr, braid_F90_ObjPtr,  braid_F90_Int *, braid_F90_ObjPtr, braid_F90_ObjPtr);
braid_Int
braid_Step_F90_Iface(braid_App        app,    /**< user-defined _braid_App structure */
                     braid_Vector     ustop,  /**< input, u vector at *tstop* */
                     braid_Vector     fstop,  /**< input, right-hand-side at *tstop* */
                     braid_Vector     u     , /**< output, u vector at *tstop* */
                     braid_StepStatus status  /**< query this struct for info about (e.g., tstart and tstop), allows for steering (e.g., set rfactor) */ 
                     )
{
   braid_Int fnotzero = 1;
   if (fstop == NULL)
   {
      fnotzero = 0;
   }
   braid_F90_Name(braid_step_f90, BRAID_STEP_F90)( 
                            braid_PassF90_Obj(     app),
                            braid_PassF90_Obj(     ustop),
                            braid_PassF90_Obj(     fstop),
                            braid_PassF90_Int(     fnotzero),
                            braid_PassF90_Obj(     u),
                            braid_PassF90_Obj(     status) );
   
   return 0;
}


/**
 * braid_BufSize
 *
 * Fortran interface, first we define the prototype for the user-defined function and then we
 * provide the C-wrapper around the user-written Fortran function
 * */
void braid_F90_Name(braid_bufsize_f90, BRAID_BUFSIZE_F90)(braid_F90_ObjPtr, braid_F90_Int *, braid_F90_ObjPtr);
braid_Int
braid_BufSize_F90_Iface(braid_App           app,              /**< user-defined _braid_App structure */
                        braid_Int           *size_ptr,         /**< upper bound on vector size in bytes */
                        braid_BufferStatus  status            /**< querry this struct for info on message type */
                        )
{
   braid_F90_Name(braid_bufsize_f90, BRAID_BUFSIZE_F90)( 
                            braid_PassF90_Obj(     app),
                            braid_PassF90_IntPtr(  size_ptr),
                            braid_PassF90_Obj(     status) );
   return 0;
}


/**
 * braid_BufPack
 *
 * Fortran interface, first we define the prototype for the user-defined function and then we
 * provide the C-wrapper around the user-written Fortran function
 * */
void braid_F90_Name(braid_bufpack_f90, BRAID_BUFPACK_F90)(braid_F90_ObjPtr, braid_F90_ObjPtr, braid_F90_Void *, braid_F90_ObjPtr);
braid_Int
braid_BufPack_F90_Iface(braid_App      app,            /**< user-defined _braid_App structure */
                        braid_Vector   u,              /**< vector to back into buffer */
                        void          *buffer,         /**< output, MPI buffer containing u */
                        braid_BufferStatus status      /**< querry this struct for info on the message type */
                        )
{
   braid_F90_Name(braid_bufpack_f90, BRAID_BUFPACK_F90)( 
                            braid_PassF90_Obj(     app),
                            braid_PassF90_Obj(     u),
                            braid_PassF90_VoidPtr( buffer),
                            braid_PassF90_Obj(     status) );
   return 0;
}


/**
 * braid_BufUnPack
 *
 * Fortran interface, first we define the prototype for the user-defined function and then we
 * provide the C-wrapper around the user-written Fortran function
 * */
void braid_F90_Name(braid_bufunpack_f90, BRAID_BUFUNPACK_F90)(braid_F90_ObjPtr, braid_F90_Void *, braid_F90_ObjPtr, braid_F90_ObjPtr);
braid_Int
braid_BufUnpack_F90_Iface(braid_App      app,          /**< user-defined _braid_App structure */
                          void          *buffer,       /**< MPI Buffer to unpack and place in u_ptr */
                          braid_Vector  *u_ptr,         /**< output, braid_Vector containing buffer's data */
                          braid_BufferStatus  status   /**< querry this structure for info on the message type */
                          )
{
   braid_F90_Name(braid_bufunpack_f90, BRAID_BUFUNPACK_F90)( 
                            braid_PassF90_Obj(     app),
                            braid_PassF90_VoidPtr( buffer),
                            braid_PassF90_ObjRef(  u_ptr),
                            braid_PassF90_Obj(     status) );
 return 0;
}

#if (braid_Fortran_Residual == 1)

/**
 * braid_Residual
 *
 * Fortran interface, first we define the prototype for the user-defined function and then we
 * provide the C-wrapper around the user-written Fortran function
 * */
void braid_F90_Name(braid_residual_f90, BRAID_RESIDUAL_F90)(braid_F90_ObjPtr, braid_F90_ObjPtr, braid_F90_ObjPtr, braid_F90_ObjPtr);
braid_Int
braid_Residual_F90_Iface(braid_App               app,        /**< user-defined _braid_App structure */
                         braid_Vector            ustop,      /**< braid_Vector to compute residual with*/                       
                         braid_Vector            r,          /**< output, residual vector */   
                         braid_StepStatus        status      /**< query this struct for info about the current status, like tstop and tstart */
                         )
{
   braid_F90_Name(braid_residual_f90, BRAID_RESIDUAL_F90)( 
                            braid_PassF90_Obj(   app),
                            braid_PassF90_Obj(   ustop),
                            braid_PassF90_Obj(   r),
                            braid_PassF90_Obj(   status) );
   return 0;
}

#endif

#if (braid_Fortran_TimeGrid == 1)

/**
 * braid_TimeGrid
 *
 * Fortran interface, first we define the prototype for the user-defined function and then we
 * provide the C-wrapper around the user-written Fortran function
 * */
void braid_F90_Name(braid_timegrid_f90, BRAID_TIMEGRID_F90)(braid_F90_ObjPtr, braid_F90_Real *, braid_F90_Int *, braid_F90_Int *);
braid_Int
braid_TimeGrid_F90_Iface(braid_App               app,       /**< user-defined _braid_App structure */
                         braid_Real             *ta,        /**< temporal grid on level 0 (slice per processor) */
                         braid_Int              *ilower,    /**< lower time index value for this processor */
                         braid_Int              *iupper     /**< upper time index value for this processor */
                         )
{
   /* Temporary scalars so that the calling function's ilower and iupper are not overwritten */
   braid_Int ilower2 = *ilower;
   braid_Int iupper2 = *iupper;
   braid_F90_Name(braid_timegrid_f90, BRAID_TIMEGRID_F90)( 
                            braid_PassF90_Obj(      app),
                            braid_PassF90_RealPtr(  ta),
                            braid_PassF90_Int(      ilower2),
                            braid_PassF90_Int(      iupper2) );
   return 0;
}

#endif

#if (braid_Fortran_SpatialCoarsen == 1)

/**
 * braid_Coarsen
 *
 * Fortran interface, first we define the prototype for the user-defined function and then we
 * provide the C-wrapper around the user-written Fortran function
 * */
void braid_F90_Name(braid_coarsen_f90, BRAID_COARSEN_F90)(braid_F90_ObjPtr, braid_F90_ObjPtr, braid_F90_ObjPtr, braid_F90_ObjPtr);
braid_Int
braid_Coarsen_F90_Iface(braid_App               app,         /**< user-defined _braid_App structure */
                        braid_Vector            fu,          /**< braid_Vector to refine*/                       
                        braid_Vector           *cu_ptr,      /**< output, refined vector */   
                        braid_CoarsenRefStatus  status       /**< query this struct for info about fu and cu (e.g., where in time fu and cu are)  */ 
                        )
{
   braid_F90_Name(braid_coarsen_f90, BRAID_COARSEN_F90)( 
                            braid_PassF90_Obj(     app),
                            braid_PassF90_Obj(     fu),
                            braid_PassF90_ObjRef(  cu_ptr),
                            braid_PassF90_Obj(     status) );
   return 0;
}


/**
 * braid_Refine
 *
 * Fortran interface, first we define the prototype for the user-defined function and then we
 * provide the C-wrapper around the user-written Fortran function
 * */
void braid_F90_Name(braid_refine_f90, BRAID_REFINE_F90)(braid_F90_ObjPtr, braid_F90_ObjPtr, braid_F90_ObjPtr, braid_F90_ObjPtr);
braid_Int
braid_Refine_F90_Iface(braid_App               app,    /**< user-defined _braid_App structure */
                       braid_Vector            cu,     /**< braid_Vector to refine*/                       
                       braid_Vector           *fu_ptr, /**< output, refined vector */       
                       braid_CoarsenRefStatus  status  /**< query this struct for info about fu and cu (e.g., where in time fu and cu are)  */ 
                       )
{
   braid_F90_Name(braid_refine_f90, BRAID_REFINE_F90)( 
                            braid_PassF90_Obj(     app),
                            braid_PassF90_Obj(     cu),
                            braid_PassF90_ObjRef(  fu_ptr),
                            braid_PassF90_Obj(     status) );
   return 0;
}

#endif


/*--------------------------------------------------------------------------
 * Define the Fortran 90 wrappers of XBraid functions.  The user will call 
 * these functions from Fortran to initialize, test and run XBraid. 
 *
 * These functions are Fortran calling C.
 *--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------
 * Wrap XBraid Test Functions
 *--------------------------------------------------------------------------*/

/* Wrap braid_TestInitAccess( ) */
braid_Int
braid_F90_Name(braid_test_init_access_f90, BRAID_TEST_INIT_ACCESS_F90)(
                          braid_F90_ObjPtr      app,     /**< User defined App structure */
                          braid_F90_Comm       *comm_x,  /**< Spatial communicator */
                          braid_F90_Real       *t        /**< Time value to test init with (used to initialize the vectors)*/
                          )
{
   braid_TestInitAccess( braid_TakeF90_Obj(braid_App, app), 
                         braid_TakeF90_Comm(          comm_x), 
                                                      stdout, 
                         braid_TakeF90_Real(          t), 
                                                      braid_Init_Vec_F90_Iface, 
                                                      braid_Access_F90_Iface, 
                                                      braid_Free_F90_Iface);
   return 0;
}

/* Wrap braid_TestClone() */
braid_Int
braid_F90_Name(braid_test_clone_f90, BRAID_TEST_CLONE_F90)( 
                          braid_F90_ObjPtr      app,     /**< User defined App structure */
                          braid_F90_Comm       *comm_x,  /**< Spatial communicator */
                          braid_F90_Real       *t        /**< Time value to test init with (used to initialize the vectors)*/
                          )
{   
   braid_TestClone( braid_TakeF90_Obj(braid_App, app), 
                    braid_TakeF90_Comm(          comm_x), 
                                                 stdout, 
                    braid_TakeF90_Real(          t), 
                                                 braid_Init_Vec_F90_Iface, 
                                                 braid_Access_F90_Iface, 
                                                 braid_Free_F90_Iface,
                                                 braid_Clone_F90_Iface);
   return 0;
}

/* Wrap braid_TestSum() */
braid_Int
braid_F90_Name(braid_test_sum_f90, BRAID_TEST_SUM_F90)( 
                          braid_F90_ObjPtr      app,     /**< User defined App structure */
                          braid_F90_Comm       *comm_x,  /**< Spatial communicator */
                          braid_F90_Real       *t        /**< Time value to test init with (used to initialize the vectors)*/
                          )
{   
   braid_TestSum( braid_TakeF90_Obj(braid_App, app), 
                  braid_TakeF90_Comm(          comm_x), 
                                               stdout, 
                  braid_TakeF90_Real(          t), 
                                               braid_Init_Vec_F90_Iface, 
                                               braid_Access_F90_Iface, 
                                               braid_Free_F90_Iface,
                                               braid_Clone_F90_Iface,
                                               braid_Sum_F90_Iface);
   return 0;
}

/* Wrap braid_TestSpatialNorm() */
braid_Int
braid_F90_Name(braid_test_spatialnorm_f90, BRAID_TEST_SPATIALNORM_F90)( 
                          braid_F90_ObjPtr      app,     /**< User defined App structure */
                          braid_F90_Comm       *comm_x,  /**< Spatial communicator */
                          braid_F90_Real       *t        /**< Time value to test init with (used to initialize the vectors)*/
                          )
{   
   braid_TestSpatialNorm( braid_TakeF90_Obj(braid_App, app), 
                          braid_TakeF90_Comm(          comm_x), 
                                                       stdout, 
                          braid_TakeF90_Real(          t), 
                                                       braid_Init_Vec_F90_Iface, 
                                                       braid_Free_F90_Iface,
                                                       braid_Clone_F90_Iface,
                                                       braid_Sum_F90_Iface,
                                                       braid_SpatialNorm_F90_Iface);
   return 0;
}

/* Wrap braid_TestBuf() */
braid_Int
braid_F90_Name(braid_test_buf_f90, BRAID_TEST_BUF_F90)( 
                          braid_F90_ObjPtr      app,     /**< User defined App structure */
                          braid_F90_Comm       *comm_x,  /**< Spatial communicator */
                          braid_F90_Real       *t        /**< Time value to test init with (used to initialize the vectors)*/
                          )
{   
   braid_TestBuf( braid_TakeF90_Obj(braid_App, app), 
                  braid_TakeF90_Comm(          comm_x), 
                                               stdout, 
                  braid_TakeF90_Real(          t), 
                                               braid_Init_Vec_F90_Iface, 
                                               braid_Free_F90_Iface,
                                               braid_Sum_F90_Iface,
                                               braid_SpatialNorm_F90_Iface,
                                               braid_BufSize_F90_Iface,
                                               braid_BufPack_F90_Iface,
                                               braid_BufUnpack_F90_Iface);
   return 0;
}


#if (braid_Fortran_SpatialCoarsen == 1)

/* Wrap braid_TestCoarsenRefine() */
braid_Int
braid_F90_Name(braid_test_coarsen_refine_f90, BRAID_TEST_COARSEN_REFINE_F90)( 
                          braid_F90_ObjPtr      app,     /**< User defined App structure */
                          braid_F90_Comm       *comm_x,  /**< Spatial communicator */
                          braid_F90_Real       *t,       /**< Time value to initialize test vectors */
                          braid_F90_Real       *fdt,     /**< Fine time step value that you spatially coarsen from */
                          braid_F90_Real       *cdt      /**< Coarse time step value that you coarsen to */
                          )
{   
   braid_TestCoarsenRefine( braid_TakeF90_Obj(braid_App, app), 
                            braid_TakeF90_Comm(          comm_x), 
                                                         stdout, 
                            braid_TakeF90_Real(          t), 
                            braid_TakeF90_Real(          fdt), 
                            braid_TakeF90_Real(          cdt), 
                                                         braid_Init_Vec_F90_Iface, 
                                                         braid_Access_F90_Iface, 
                                                         braid_Free_F90_Iface,
                                                         braid_Clone_F90_Iface,
                                                         braid_Sum_F90_Iface,
                                                         braid_SpatialNorm_F90_Iface,
                                                         braid_Coarsen_F90_Iface,
                                                         braid_Refine_F90_Iface);
   return 0;
}

#endif

/* Wrap braid_TestAll() */
braid_Int
braid_F90_Name(braid_test_all_f90, BRAID_TEST_ALL_F90)( 
                          braid_F90_ObjPtr      app,     /**< User defined App structure */
                          braid_F90_Comm       *comm_x,  /**< Spatial communicator */
                          braid_F90_Real       *t,       /**< Time value to initialize test vectors */
                          braid_F90_Real       *fdt,     /**< Fine time step value that you spatially coarsen from */
                          braid_F90_Real       *cdt      /**< Coarse time step value that you coarsen to */
                          )
{   


#if (braid_Fortran_SpatialCoarsen != 0)
   braid_PtFcnSCoarsen coarsen_fcn = braid_Coarsen_F90_Iface;
   braid_PtFcnSRefine refine_fcn = braid_Refine_F90_Iface;
#else
   braid_PtFcnSCoarsen coarsen_fcn = NULL;
   braid_PtFcnSRefine refine_fcn = NULL;
#endif

#if (braid_Fortran_Residual != 0)
   braid_PtFcnStep step_fcn = braid_Step_F90_Iface;
   braid_PtFcnResidual residual_fcn = braid_Residual_F90_Iface;
#else
   braid_PtFcnStep step_fcn = NULL;
   braid_PtFcnResidual residual_fcn = NULL;
#endif


   braid_TestAll( braid_TakeF90_Obj(braid_App, app), 
                  braid_TakeF90_Comm(          comm_x), 
                                               stdout, 
                  braid_TakeF90_Real(          t), 
                  braid_TakeF90_Real(          fdt), 
                  braid_TakeF90_Real(          cdt), 
                                               braid_Init_Vec_F90_Iface, 
                                               braid_Free_F90_Iface,
                                               braid_Clone_F90_Iface,
                                               braid_Sum_F90_Iface,
                                               braid_SpatialNorm_F90_Iface,
                                               braid_BufSize_F90_Iface,
                                               braid_BufPack_F90_Iface,
                                               braid_BufUnpack_F90_Iface,
                                               coarsen_fcn,
                                               refine_fcn,
                                               residual_fcn,
                                               step_fcn);

   return 0;
}


/*--------------------------------------------------------------------------
 * Wrap XBraid Status Query Functions
 *--------------------------------------------------------------------------*/

/* Wrap braid_AccessStatusGetTILD( ) */
braid_Int
braid_F90_Name(braid_access_status_get_tild_f90, BRAID_ACCESS_STATUS_GET_TILD_F90)(
                          braid_F90_ObjPtr    status,       /**< structure containing current simulation info */
                          braid_F90_Real     *t_ptr,        /**< output,  current time */
                          braid_F90_Int      *iter_ptr,     /**< output,  current iteration in XBraid*/
                          braid_F90_Int      *level_ptr,    /**< output,  current level in XBraid */
                          braid_F90_Int      *done_ptr      /**< output,  boolean describing whether XBraid has finished */
                          )
{
   braid_AccessStatusGetTILD( braid_TakeF90_Obj(braid_AccessStatus,  status),
                              braid_TakeF90_RealPtr(                 t_ptr),
                              braid_TakeF90_IntPtr(                  iter_ptr),
                              braid_TakeF90_IntPtr(                  level_ptr),
                              braid_TakeF90_IntPtr(                  done_ptr) );
   return 0;
}

/* Wrap braid_AccessStatusGetTIndex( ) */
braid_Int
braid_F90_Name(braid_access_status_get_tindex_f90, BRAID_ACCESS_STATUS_GET_TINDEX_F90)(
                              braid_F90_ObjPtr    status,        /**< structure containing current simulation info */
                              braid_F90_Int      *idx_ptr      /**< output, global index value corresponding to current time value to evolve from */
                              )
{
   braid_AccessStatusGetTIndex( braid_TakeF90_Obj(braid_AccessStatus, status),
                               braid_TakeF90_IntPtr(                 idx_ptr) );
   return 0;
}

/* Wrap braid_StepStatusGetTstartTstop( ) */
braid_Int
braid_F90_Name(braid_step_status_get_tstart_tstop_f90, BRAID_STEP_STATUS_GET_TSTART_TSTOP_F90)(
                              braid_F90_ObjPtr    status,        /**< structure containing current simulation info */
                              braid_F90_Real     *tstart_ptr,    /**< output, current time */
                              braid_F90_Real     *tstop_ptr      /**< output, next time value to evolve towards */
                              )
{
   braid_StepStatusGetTstartTstop( braid_TakeF90_Obj(braid_StepStatus, status),
                                   braid_TakeF90_RealPtr(              tstart_ptr),
                                   braid_TakeF90_RealPtr(              tstop_ptr) );
   return 0;
}

/* Wrap braid_StepStatusGetTIndex( ) */
braid_Int
braid_F90_Name(braid_step_status_get_tindex_f90, BRAID_STEP_STATUS_GET_TINDEX_F90)(
                              braid_F90_ObjPtr    status,        /**< structure containing current simulation info */
                              braid_F90_Int      *idx_ptr      /**< output, global index value corresponding to next time value to evolve towards */
                              )
{
   braid_StepStatusGetTIndex( braid_TakeF90_Obj(braid_StepStatus, status),
                               braid_TakeF90_IntPtr(             idx_ptr) );
   return 0;
}

/* Wrap braid_StepStatusGetLevel( ) */
braid_Int
braid_F90_Name(braid_step_status_get_level_f90, BRAID_STEP_STATUS_GET_LEVEL_F90)(
                                  braid_F90_ObjPtr    status,        /**< structure containing current simulation info */
                                  braid_F90_Int      *level_ptr      /**< output, current level in XBraid */
                                  )
{
   braid_StepStatusGetLevel( braid_TakeF90_Obj(braid_StepStatus,  status),
                             braid_TakeF90_IntPtr(                level_ptr) );
   return 0;
}


/* Wrap braid_StepStatusSetRFactor( ) */
braid_Int
braid_F90_Name(braid_step_status_set_rfactor_f90, BRAID_STEP_STATUS_SET_RFACTOR_F90)(
                          braid_F90_ObjPtr    status,         /**< structure containing current simulation info */
                          braid_F90_Int      *rfactor         /**< user-determined desired rfactor */
                          )
{
   braid_StepStatusSetRFactor( braid_TakeF90_Obj(braid_StepStatus, status),
                               braid_TakeF90_Int(                  rfactor));
   return 0;
}

/* Wrap braid_StepStatusGetTol( ) */
braid_Int
braid_F90_Name(braid_step_status_get_tol_f90, BRAID_STEP_STATUS_GET_TOL_F90)(
                              braid_F90_ObjPtr    status,        /**< structure containing current simulation info */
                              braid_F90_Real      *tol_ptr       /**< output, current XBraid stopping tolerance */
                              )
{
   braid_StepStatusGetTol( braid_TakeF90_Obj(braid_StepStatus, status),
                           braid_TakeF90_RealPtr(              tol_ptr) );
   return 0;
}

/* Wrap braid_StepStatusGetIter( ) */
braid_Int
braid_F90_Name(braid_step_status_get_iter_f90, BRAID_STEP_STATUS_GET_ITER_F90)(
                              braid_F90_ObjPtr    status,        /**< structure containing current simulation info */
                              braid_F90_Int       *iter_ptr      /**< output, current iteration in XBraid */
                              )
{
   braid_StepStatusGetIter( braid_TakeF90_Obj(braid_StepStatus, status),
                            braid_TakeF90_IntPtr(               iter_ptr) );
   return 0;
}

/* Wrap braid_StepStatusGetRNorms( ) */
braid_Int
braid_F90_Name(braid_step_status_get_rnorms_f90, BRAID_STEP_STATUS_GET_RNORMS_F90)(
                             braid_F90_ObjPtr  status,            /**< structure containing current simulation info */
                             braid_F90_Int     *nrequest_ptr,     /**< input/output, input: num requested resid norms, output: num actually returned */
                             braid_F90_Real    *rnorm_ptr         /**< output, holds residual norm history array */
                             )
{
   braid_StepStatusGetRNorms(braid_TakeF90_Obj(braid_StepStatus, status),
                             braid_TakeF90_IntPtr(               nrequest_ptr),
                             braid_TakeF90_RealPtr(              rnorm_ptr) );
   return 0;
}


/* Wrap braid_StepStatusGetOldFineTolx( ) */
braid_Int
braid_F90_Name(braid_step_status_get_old_fine_tolx_f90, BRAID_STEP_STATUS_GET_OLD_FINE_TOLX_F90)(
                              braid_F90_ObjPtr    status,              /**< structure containing current simulation info */
                              braid_F90_Real      *old_fine_tolx_ptr   /**< output, previous *old_fine_tolx*, set through *braid_StepStatusSetOldFineTolx* */                              
                              )
{
   braid_StepStatusGetOldFineTolx( braid_TakeF90_Obj(braid_StepStatus, status),
                                   braid_TakeF90_RealPtr(              old_fine_tolx_ptr) );
   return 0;
}

/* Wrap braid_StepStatusSetOldFineTolx( ) */
braid_Int
braid_F90_Name(braid_step_status_set_old_fine_tolx_f90, BRAID_STEP_STATUS_SET_OLD_FINE_TOLX_F90)(
                              braid_F90_ObjPtr    status,              /**< structure containing current simulation info */
                              braid_F90_Real      *old_fine_tolx       /**< input, the last used fine_tolx */                              
                              )
{
   braid_StepStatusSetOldFineTolx( braid_TakeF90_Obj(braid_StepStatus, status),
                                   braid_TakeF90_Real(                 old_fine_tolx) );
   return 0;
}

/* Wrap braid_StepStatusSetTightFineTolx( ) */
braid_Int
braid_F90_Name(braid_step_status_set_tight_fine_tolx_f90, BRAID_STEP_STATUS_SET_TIGHT_FINE_TOLX_F90)(
                              braid_F90_ObjPtr    status,              /**< structure containing current simulation info */
                              braid_F90_Real      *tight_fine_tolx     /**< input, boolean indicating whether the tight tolx has been used */
                              )
{
   braid_StepStatusSetTightFineTolx( braid_TakeF90_Obj(braid_StepStatus, status),
                                     braid_TakeF90_Real(                 tight_fine_tolx) );
   return 0;
}

/* Wrap braid_BufferStatusGetMessageType( ) */
braid_Int
braid_F90_Name(braid_buffer_status_get_message_type_f90, BRAID_BUFFER_STATUS_GET_MESSAGE_TYPE_F90)(
                              braid_F90_ObjPtr     status,            /**< structure containing current simulation info */
                              braid_F90_Int        *messagetype_ptr   /**< output, type of message, 0: for Step(), 1: for load balancing */  
                              )
{
   braid_BufferStatusGetMessageType(braid_TakeF90_Obj( braid_BufferStatus, status),
                                    braid_TakeF90_IntPtr(                      messagetype_ptr) );
   return 0;
}

/* Wrap braid_BufferStatusGetSetSize( ) */
braid_Int
braid_F90_Name(braid_buffer_status_set_size_f90, BRAID_BUFFER_STATUS_SET_SIZE_F90)(
                              braid_F90_ObjPtr     status,            /**< structure containing current simulation info */
                              braid_F90_Int        *size_ptr          /**< input, current message size in bytes */
                              )
{
   braid_BufferStatusSetSize(braid_TakeF90_Obj( braid_BufferStatus, status),
                             braid_TakeF90_Int(                size_ptr) );
   return 0;
}

/*--------------------------------------------------------------------------
 * Wrap XBraid User Interface Functions
 *--------------------------------------------------------------------------*/

/*  braid_Init( ) */
braid_Int
braid_F90_Name(braid_init_f90, BRAID_INIT_F90)(
           braid_F90_Comm        *comm_world,  /**< Global communicator for space and time */
           braid_F90_Comm        *comm,        /**< Communicator for temporal dimension*/
           braid_F90_Real        *tstart,      /**< start time */
           braid_F90_Real        *tstop,       /**< End time*/
           braid_F90_Int         *ntime,       /**< Initial number of temporal grid values*/
           braid_F90_ObjPtr       app,         /**< User-defined _braid_App structure */
           braid_F90_ObjPtr      *core_ptr     /**< Pointer to braid_Core (_braid_Core) struct*/   
           )
{
   braid_Init(braid_TakeF90_Comm(              comm_world),
              braid_TakeF90_Comm(              comm),
              braid_TakeF90_Real(              tstart),
              braid_TakeF90_Real(              tstop),
              braid_TakeF90_Int(               ntime),
              braid_TakeF90_Obj(braid_App,     app),
              braid_Step_F90_Iface, 
              braid_Init_Vec_F90_Iface, 
              braid_Clone_F90_Iface,
              braid_Free_F90_Iface,
              braid_Sum_F90_Iface,
              braid_SpatialNorm_F90_Iface,
              braid_Access_F90_Iface,
              braid_BufSize_F90_Iface,
              braid_BufPack_F90_Iface,
              braid_BufUnpack_F90_Iface,
              braid_TakeF90_ObjPtr(braid_Core, core_ptr)
              );

   return 0;
}

/*  braid_Drive( ) */
braid_Int
braid_F90_Name(braid_drive_f90, BRAID_DRIVE_F90)(
                     braid_F90_ObjPtr  *core                /**< braid_Core (_braid_Core) struct*/
                     )
{
   braid_Drive( braid_TakeF90_ObjDeref(braid_Core,     core) );
   return 0;
}

/*  braid_Destroy( ) */
braid_Int
braid_F90_Name(braid_destroy_f90, BRAID_DESTROY_F90)(
                   braid_F90_ObjPtr   *core                /**< braid_Core (_braid_Core) struct*/
                   )
{
   braid_Destroy( braid_TakeF90_ObjDeref(braid_Core,     core) );
   return 0;
}

/*  braid_PrintStats( ) */
braid_Int
braid_F90_Name(braid_print_stats_f90, BRAID_PRINT_STATS_F90)(
                   braid_F90_ObjPtr   *core                /**< braid_Core (_braid_Core) struct*/
                   )
{
   braid_PrintStats( braid_TakeF90_ObjDeref(braid_Core,     core) );
   return 0;
}

/*  braid_SetMaxLevels( ) */
braid_Int
braid_F90_Name(braid_set_max_levels_f90, BRAID_SET_MAX_LEVELS_F90)(
                   braid_F90_ObjPtr  *core,        /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Int     *max_levels   /**< maximum levels to allow */
                   )
{
   braid_SetMaxLevels(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                      braid_TakeF90_Int(                  max_levels) );
   return 0;
}

/*  braid_SetSkip( ) */
braid_Int
braid_F90_Name(braid_set_skip_f90, BRAID_SET_SKIP_F90)(
                   braid_F90_ObjPtr  *core,        /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Int     *skip         /**< boolean, whether to skip work on down cycle */
                   )
{
   braid_SetSkip(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                 braid_TakeF90_Int(                  skip) );
   return 0;
}

/*  braid_SetMinCoarse( ) */
braid_Int
braid_F90_Name(braid_set_min_coarse_f90, BRAID_SET_MIN_COARSE_F90)(
                   braid_F90_ObjPtr  *core,        /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Int     *min_coarse   /**< maximum levels to allow */
                   )
{
   braid_SetMinCoarse(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                      braid_TakeF90_Int(                  min_coarse) );
   return 0;
}

/*  braid_SetAbsTol( ) */
braid_Int
braid_F90_Name(braid_set_abs_tol_f90, BRAID_SET_ABS_TOL_F90)(
                   braid_F90_ObjPtr  *core,        /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Real    *atol         /**< absolute stopping tolerance */
                   )
{
   braid_SetAbsTol(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                   braid_TakeF90_Real(                 atol) );
   return 0;
}

/* braid_SetRelTol( ) */
braid_Int
braid_F90_Name(braid_set_rel_tol_f90, BRAID_SET_REL_TOL_F90)(
                   braid_F90_ObjPtr  *core,        /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Real    *rtol         /**< relative stopping tolerance */
                   )
{
   braid_SetRelTol(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                   braid_TakeF90_Real(                 rtol) );
   return 0;
}

/* braid_SetNRelax( ) */
braid_Int
braid_F90_Name(braid_set_nrelax_f90, BRAID_SET_NRELAX_F90)(
                   braid_F90_ObjPtr  *core,        /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Int     *level,       /**< *level* to set *nrelax* on */
                   braid_F90_Int     *nrelax       /**< number of relaxations to do on *level* */
                   )
{
   braid_SetNRelax(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                   braid_TakeF90_Int(                  level),
                   braid_TakeF90_Int(                  nrelax) );
   return 0;
}

/* braid_SetCFactor( ) */
braid_Int
braid_F90_Name(braid_set_cfactor_f90, BRAID_SET_CFACTOR_F90)(
                   braid_F90_ObjPtr  *core,        /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Int     *level,       /**< *level* to set coarsening factor on */
                   braid_F90_Int     *cfactor      /**< desired coarsening factor */
                   )
{
   braid_SetCFactor(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                    braid_TakeF90_Int(                  level),
                    braid_TakeF90_Int(                  cfactor) );
   return 0;
}

/*  braid_SetMaxIter( ) */
braid_Int
braid_F90_Name(braid_set_max_iter_f90, BRAID_SET_MAX_ITER_F90)(
                   braid_F90_ObjPtr  *core,        /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Int     *max_iter     /**< maximum iterations to allow */
                   )
{
   braid_SetMaxIter(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                    braid_TakeF90_Int(                  max_iter) );
   return 0;
}

/*  braid_SetFMG( ) */
braid_Int
braid_F90_Name(braid_set_fmg_f90, BRAID_SET_FMG_F90)(
                   braid_F90_ObjPtr   *core         /**< braid_Core (_braid_Core) struct*/
                   )
{
   braid_SetFMG(braid_TakeF90_ObjDeref(braid_Core,  core) );
   return 0;
}

/*  braid_SetNFMG( ) */
braid_Int
braid_F90_Name(braid_set_nfmg_f90, BRAID_SET_NFMG_F90)(
                   braid_F90_ObjPtr   *core,         /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Int      *nfmg          /**< number of initial F-cycles to do before switching to V-cycles */
                   )
{
   braid_SetNFMG(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                     braid_TakeF90_Int(                  nfmg) );
   return 0;
}

/* braid_SetTemporalNorm( ) */
braid_Int
braid_F90_Name(braid_set_temporal_norm_f90, BRAID_SET_TEMPORAL_NORM_F90)(
                   braid_F90_ObjPtr  *core,        /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Int     *tnorm        /**< choice of temporal norm*/
                   )
{
   braid_SetTemporalNorm(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                         braid_TakeF90_Int(                  tnorm) );
   return 0;
}

/* braid_SetNFMGVcyc( ) */
braid_Int
braid_F90_Name(braid_set_nfmg_vcyc_f90, BRAID_SET_NFMG_VCYC_F90)(
                   braid_F90_ObjPtr  *core,        /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Int     *nfmg_Vcyc    /**< number of V-cycles to do each FMG level */
                   )
{
   braid_SetNFMGVcyc(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                     braid_TakeF90_Int(                  nfmg_Vcyc) );
   return 0;
}

/* braid_SetStorage( ) */
braid_Int
braid_F90_Name(braid_set_storage_f90, BRAID_SET_STORAGE_F90)(
                   braid_F90_ObjPtr  *core,        /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Int     *storage      /**< store C-points (0), all points (1) */
                   )
{
   braid_SetStorage(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                    braid_TakeF90_Int(                  storage) );
   return 0;
}

#if (braid_Fortran_Residual == 1)

/* braid_SetResidual( ) */
braid_Int
braid_F90_Name(braid_set_residual_f90, BRAID_SET_RESIDUAL_F90)(
                   braid_F90_ObjPtr   *core         /**< braid_Core (_braid_Core) struct*/
                   )
{
   braid_SetResidual(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                    (braid_PtFcnResidual) braid_F90_Name(braid_residual_f90,  BRAID_RESIDUAL_F90) );
   return 0;
}

#endif

#if (braid_Fortran_TimeGrid == 1)

/* braid_SetTimeGrid( ) */
braid_Int
braid_F90_Name(braid_set_timegrid_f90, BRAID_SET_TIMEGRID_F90)(
                   braid_F90_ObjPtr   *core         /**< braid_Core (_braid_Core) struct*/
                   )
{
   braid_SetTimeGrid(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                     (braid_PtFcnTimeGrid) braid_F90_Name(braid_timegrid_f90, BRAID_TIMEGRID_F90) );
   return 0;
}

#endif

#if (braid_Fortran_SpatialCoarsen == 1)

/* braid_SetSpatialCoarsen( ) */
braid_Int
braid_F90_Name(braid_set_spatial_coarsen_f90, BRAID_SET_SPATIAL_COARSEN_F90)(
                   braid_F90_ObjPtr   *core        /**< braid_Core (_braid_Core) struct*/
                   )
{
   braid_SetSpatialCoarsen(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                          (braid_PtFcnSCoarsen) braid_F90_Name(braid_coarsen_f90,   BRAID_COARSEN_F90) );
   return 0;
}

/* braid_SetSpatialRefine( ) */
braid_Int
braid_F90_Name(braid_set_spatial_refine_f90, BRAID_SET_SPATIAL_REFINE_F90)(
                   braid_F90_ObjPtr   *core        /**< braid_Core (_braid_Core) struct*/
                   )
{
   braid_SetSpatialRefine(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                          (braid_PtFcnSRefine) braid_F90_Name(braid_refine_f90,    BRAID_REFINE_F90) );
   return 0;
}

#endif

/* braid_SetPrintLevel( ) */
braid_Int
braid_F90_Name(braid_set_print_level_f90, BRAID_SET_PRINT_LEVEL_F90)(
                   braid_F90_ObjPtr  *core,          /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Int     *print_level    /**< desired print level */
                   )
{
   braid_SetPrintLevel(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                       braid_TakeF90_Int(                  print_level) );
   return 0;
}

/* braid_SetPrintFile( ) */
braid_Int
braid_F90_Name(braid_set_print_file_f90, BRAID_SET_PRINT_FILE_F90)(
                   braid_F90_ObjPtr  *core,          /**< braid_Core (_braid_Core) struct*/
                   char              *fname,         /**< name of string to set */
                   braid_F90_Int     *ll             /**< string length */
                   )
{
   braid_SetPrintFile(braid_TakeF90_ObjDeref(braid_Core,  core),
                                                          fname);
   return 0;
}

/* braid_SetDefaultPrintFile( ) */
braid_Int
braid_F90_Name(braid_set_default_print_file_f90, braid_set_default_print_file_F90)(
                   braid_F90_ObjPtr  *core           /**< braid_Core (_braid_Core) struct*/
                   )
{
   braid_SetDefaultPrintFile(braid_TakeF90_ObjDeref(braid_Core,  core));
   return 0;
}

/* braid_SetAccessLevel( ) */
braid_Int
braid_F90_Name(braid_set_access_level_f90, BRAID_SET_ACCESS_LEVEL_F90)(
                   braid_F90_ObjPtr  *core,          /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Int     *access_level   /**< desired access_level */
                   )
{
   braid_SetAccessLevel(braid_TakeF90_ObjDeref(braid_Core,  core),
                        braid_TakeF90_Int(                  access_level) );
   return 0;
}


/* braid_SplitCommworld( ) */
braid_Int
braid_F90_Name(braid_split_commworld_f90, BRAID_SPLIT_COMMWORLD_F90)(
                   braid_F90_Comm  *comm_world,    /**< Global communicator to split */
                   braid_F90_Int   *px,            /**< Number of processors parallelizing space for a single time step*/
                   braid_F90_Comm  *comm_x,        /**< Spatial communicator (written as output) */
                   braid_F90_Comm  *comm_t         /**< Temporal communicator (written as output) */
)
{
   MPI_Comm Fcommworld = braid_TakeF90_Comm( comm_world);
   MPI_Comm Fcomm_x, Fcomm_t;
   braid_SplitCommworld(                       &Fcommworld,
                        braid_TakeF90_Int(     px),
                                               &Fcomm_x,
                                               &Fcomm_t);
   *comm_x = braid_PassF90_Comm( &Fcomm_x );
   *comm_t = braid_PassF90_Comm( &Fcomm_t );

   return 0;
}


/* braid_GetNumIter( ) */
braid_Int
braid_F90_Name(braid_get_num_iter_f90, BRAID_GET_NUM_ITER_F90)(
                   braid_F90_ObjPtr  *core,        /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Int     *niter_ptr    /**< output, holds number of iterations taken */
                   )
{
   braid_GetNumIter(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                    braid_TakeF90_IntPtr(               niter_ptr) );
   return 0;
}

/* braid_GetRNorms( ) */
braid_Int
braid_F90_Name(braid_get_rnorms_f90, BRAID_GET_RNORMS_F90)(
                   braid_F90_ObjPtr  *core,         /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Int     *nrequest_ptr, /**< input/output, input: num requested resid norms, output: num actually returned */
                   braid_F90_Real    *rnorm_ptr     /**< output, holds residual norm history array */
                   )
{
   braid_GetRNorms(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                   braid_TakeF90_IntPtr(               nrequest_ptr),
                   braid_TakeF90_RealPtr(              rnorm_ptr) );
   return 0;
}

/* braid_GetNLevels( ) */
braid_Int
braid_F90_Name(braid_get_nlevels_f90, BRAID_GET_NLEVELS_F90)(
                   braid_F90_ObjPtr   *core,        /**< braid_Core (_braid_Core) struct*/
                   braid_F90_Int      *nlevels_ptr  /**< output, holds the number of XBraid levels */
                   )
{
   braid_GetNLevels(braid_TakeF90_ObjDeref(braid_Core,  core) ,
                    braid_TakeF90_IntPtr(               nlevels_ptr) );
   return 0;
}

/* braid_GetSpatialAccuracy( ) */
braid_Int
braid_F90_Name(braid_get_spatial_accuracy_f90, BRAID_GET_SPATIAL_ACCURACY_F90)(
               braid_F90_ObjPtr    status,        /**< structure containing current simulation info */
               braid_F90_Real      *loose_tol,    /**< Loosest allowed spatial solve stopping tol on fine grid*/
               braid_F90_Real      *tight_tol,    /**< Tightest allowed spatial solve stopping tol on fine grid*/
               braid_F90_Real      *tol_ptr       /**< output, holds the computed spatial solve stopping tol */
               )
{
   braid_GetSpatialAccuracy(braid_TakeF90_Obj(braid_StepStatus, status),   
                            braid_TakeF90_Real(                 loose_tol),
                            braid_TakeF90_Real(                 tight_tol),
                            braid_TakeF90_RealPtr(              tol_ptr));

   return 0;
}

/* braid_SetSeqSoln( ) */
braid_Int
braid_F90_Name(braid_set_seq_soln_f90, BRAID_SET_SEQ_SOLN_F90)(
                                  braid_F90_ObjPtr      *core_ptr,     /**< braid_Core (_braid_Core) struct*/
                                  braid_F90_Int         *use_ptr       /**< 1: Init with sequential time stepping soln, 0: Use user's Init()*/
                                  )
{
   braid_SetSeqSoln(braid_TakeF90_ObjDeref(braid_Core, core_ptr),
                    braid_TakeF90_Int(                 use_ptr) );
   return 0;
}



#endif

