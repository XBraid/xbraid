/*BHEADER**********************************************************************
 * Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of XBraid.  See file COPYRIGHT for details.
 *
 * XBraid is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 ***********************************************************************EHEADER*/

/** \file braid_status.h
 * \brief Define headers for XBraid status structures, status get/set routines
 * and status create/destroy routines.
 *
 */

#ifndef braid_status_HEADER
#define braid_status_HEADER

#include "braid_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------------
 * Routines for user to access XBraid status structures
 *--------------------------------------------------------------------------*/

/** \defgroup braidstatus XBraid status routines
 *  \ingroup userinterface
 *  
 *  XBraid status structures are what tell the user the status of the simulation
 *  when their routines (phi, coarsen/refine, access) are called. 
 *
 *  @{
 */


/*--------------------------------------------------------------------------
 * Define Access Status Structure
 *--------------------------------------------------------------------------*/

struct _braid_AccessStatus_struct;
/**
 * The user access routine will receive braid_AccessStatus, which will be 
 * a pointer to the actual _braid_AccessStatus_struct
 **/
typedef struct _braid_AccessStatus_struct *braid_AccessStatus;

/** 
 * AccessStatus structure which defines the status of XBraid at a given instant
 * on some level during a run.  The user accesses it through
 * _braid_AccessStatusGet**()_ functions.
 **/
typedef struct _braid_AccessStatus_struct
{
   braid_Real    t;            /**< current time */
   braid_Int     iter;         /**< XBraid iteration number */
   braid_Int     level;        /**< current level in XBraid*/
   braid_Real    rnorm;        /**< residual norm */
   braid_Int     done;         /**< boolean describing whether XBraid has finished */
   
} _braid_AccessStatus;


/*--------------------------------------------------------------------------
 * Define CoarsenRef Status Structure
 *--------------------------------------------------------------------------*/

struct _braid_CoarsenRefStatus_struct;
/**
 * The user coarsen and refine routines will receive braid_CoarsenRefStatus, which will be 
 * a pointer to the actual _braid_CoarsenRefStatus_struct
 **/
typedef struct _braid_CoarsenRefStatus_struct *braid_CoarsenRefStatus;

/** 
 * The user coarsen and refine routines will receive a CoarsenRefStatus structure, which 
 * defines the status of XBraid at a given instant of coarsening or refinement on some level 
 * during a run.  The user accesses it through _braid_CoarsenRefStatusGet**()_ functions.
 **/
typedef struct _braid_CoarsenRefStatus_struct
{
   braid_Real     tstart;      /**< current time value */                          
   braid_Real     f_tminus;    /**< time value for to the left on fine grid */ 
   braid_Real     f_tplus;     /**< time value for to the right on fine grid */
   braid_Real     c_tminus;    /**< time value for to the left on coarse grid */
   braid_Real     c_tplus;     /**< time value for to the right on coarse grid */
   
} _braid_CoarsenRefStatus;


/*--------------------------------------------------------------------------
 * Define Phi Status Structure
 *--------------------------------------------------------------------------*/

struct _braid_PhiStatus_struct;
/** 
 * The user's phi routine will receive braid_PhiStatus, which will be a
 * pointer to the actual _braid_PhiStatus_struct
 **/
typedef struct _braid_PhiStatus_struct *braid_PhiStatus;

/** 
 * The user's phi routine routine will receive a PhiStatus structure, which 
 * defines the status of XBraid at the given instant for phi evaluation on some level 
 * during a run.  The user accesses it through _braid_PhiStatusGet**()_ functions.
 **/
typedef struct _braid_PhiStatus_struct
{
   braid_Real     tstart;          /**< current time value  */
   braid_Real     tplus;           /**< time value to evolve towards, time value to the right */
   braid_Real     accuracy;        /**< advanced option allowing variable accuracy for implicit phi*/
   braid_Int      rfactor;         /**< if set by user, allows for subdivision of this interval for bettter time accuracy */
} _braid_PhiStatus;


/*--------------------------------------------------------------------------
 * Accessor macros 
 *--------------------------------------------------------------------------*/

/**
 * Accessor for all _braid_**Status attributes 
 **/
#define _braid_StatusElt(status, elt) ( (status) -> elt )


/*--------------------------------------------------------------------------
 * AccessStatus Prototypes
 *--------------------------------------------------------------------------*/

/**
 * Initialize a braid_AccessStatus structure  
 **/
braid_Int
_braid_AccessStatusInit(braid_Real          t,           /**< current time */
                        braid_Real          rnorm,       /**< current residual norm in XBraid */
                        braid_Int           iter,        /**< current iteration in XBraid*/
                        braid_Int           level,       /**< current level in XBraid */
                        braid_Int           done,        /**< boolean describing whether XBraid has finished */
                        braid_AccessStatus  status       /**< structure to initialize */
                        );

/**
 * Destroy a braid_AccessStatus structure
 **/
braid_Int
_braid_AccessStatusDestroy(braid_AccessStatus  status);        /**< structure to be destroyed */

/**
 * Return the current time from the AccessStatus structure.
 **/
braid_Int
braid_AccessStatusGetT(braid_AccessStatus  status,     /**< structure containing current simulation info */
                       braid_Real         *t_ptr       /**< output, current time */
                       );
/**
 * Return the current residual norm from the AccessStatus structure.
 **/
braid_Int
braid_AccessStatusGetResidual(braid_AccessStatus  status,     /**< structure containing current simulation info */
                              braid_Real         *rnorm_ptr   /**< output, current residual norm */
                              );

/**
 * Return the current iteration from the AccessStatus structure.
 **/
braid_Int
braid_AccessStatusGetIter(braid_AccessStatus  status,         /**< structure containing current simulation info */
                          braid_Int          *iter_ptr        /**< output, current XBraid iteration number*/
                          );

/**
 * Return the current XBraid level from the AccessStatus structure.
 **/
braid_Int
braid_AccessStatusGetLevel(braid_AccessStatus  status,        /**< structure containing current simulation info */
                           braid_Int          *level_ptr      /**< output, current level in XBraid */
                           );

/**
 * Return whether XBraid is done for the current simulation.
 *
 * *done_ptr = 1* indicates that XBraid has finished iterating, 
 * (either maxiter has been reached, or the tolerance has been met).
 **/
braid_Int
braid_AccessStatusGetDone(braid_AccessStatus  status,         /**< structure containing current simulation info */
                          braid_Int          *done_ptr        /**< output,  =1 if XBraid has finished, else =0 */
                          );

/**
 * Return XBraid status for the current simulation. Four values are 
 * returned.
 *
 * TILD : time, iteration, level, done
 *
 * These values are also available through individual Get routines. 
 * These individual routines are the location of detailed documentation on 
 * each parameter, e.g., see *braid_AccessStatusGetDone* for more information
 * on the *done* value.
 **/
braid_Int
braid_AccessStatusGetTILD(braid_AccessStatus  status,       /**< structure containing current simulation info */
                          braid_Real          *t_ptr,       /**< output,  current time */
                          braid_Int           *iter_ptr,    /**< output,  current iteration in XBraid*/
                          braid_Int           *level_ptr,   /**< output,  current level in XBraid */
                          braid_Int           *done_ptr     /**< output,  boolean describing whether XBraid has finished */
                          );

/*--------------------------------------------------------------------------
 * CoarsenRefStatus Prototypes
 *--------------------------------------------------------------------------*/

/**
 * Initialize a braid_CoarsenRefStatus structure 
 **/
braid_Int
_braid_CoarsenRefStatusInit(braid_Real              tstart,      /**< time value for current vector */             
                            braid_Real              f_tminus,    /**< time value for to the left on fine grid */ 
                            braid_Real              f_tplus,     /**< time value for to the right on fine grid */
                            braid_Real              c_tminus,    /**< time value for to the left on coarse grid */
                            braid_Real              c_tplus,     /**< time value for to the right on coarse grid */
                            braid_CoarsenRefStatus  status       /**< structure to initialize */
                            );


/**
 * Destroy a braid_CoarsenRefStatus structure
 **/
braid_Int
_braid_CoarsenRefStatusDestroy(braid_CoarsenRefStatus  status);        /**< structure to be destroyed */

/**
 * Return the current time value from the CoarsenRefStatus structure.
 **/
braid_Int
braid_CoarsenRefStatusGetTstart(braid_CoarsenRefStatus  status,         /**< structure containing current simulation info */
                                braid_Real             *tstart_ptr      /**< output, current time */
                                );

/** 
 * Return the time value to the right of the current **fine** time grid value from
 * the CoarsenRefStatus structure.
 **/
braid_Int
braid_CoarsenRefStatusGetFTplus(braid_CoarsenRefStatus  status,         /**< structure containing current simulation info */
                                braid_Real             *f_tplus_ptr     /**< output, time value to the right on fine grid */
                                );

/**
 * Return the time value to the left of the current **fine** time grid value from
 * the CoarsenRefStatus structure.
 **/
braid_Int
braid_CoarsenRefStatusGetFTminus(braid_CoarsenRefStatus  status,        /**< structure containing current simulation info */
                                 braid_Real             *f_tminus_ptr   /**< output, time value to the left on fine grid */
                                 );

/**
 * Return the time value to the right of the current **coarse** time grid value from
 * the CoarsenRefStatus structure.
 **/
braid_Int
braid_CoarsenRefStatusGetCTplus(braid_CoarsenRefStatus  status,         /**< structure containing current simulation info */
                                braid_Real             *c_tplus_ptr     /**< output, time value to the right on coarse grid */
                                );

/**
 * Return the time value to the left of the current **coarse** time grid value from
 * the CoarsenRefStatus structure.
 **/
braid_Int
braid_CoarsenRefStatusGetCTminus(braid_CoarsenRefStatus  status,        /**< structure containing current simulation info */
                                 braid_Real             *c_tminus_ptr   /**< output, time value to the left on coarse grid */
                                 );

/**
 * Return XBraid status for the current simulation. Five values are 
 * returned, tstart, f_tminus, f_tplus, c_tminus,  c_tplus. 
 *
 * These values are also available through individual Get routines. 
 * These individual routines are the location of detailed documentation on 
 * each parameter, e.g., see *braid_CoarsenRefStatusGetCTminus* for more 
 * information on the *c_tminus* value.
 **/
braid_Int
braid_CoarsenRefStatusGetTminusTplus(braid_CoarsenRefStatus  status,           /**< structure containing current simulation info */
                                     braid_Real              *tstart_ptr,      /**< output, time value for current vector */             
                                     braid_Real              *f_tminus_ptr,    /**< output, time value for to the left on fine grid */ 
                                     braid_Real              *f_tplus_ptr,     /**< output, time value for to the right on fine grid */
                                     braid_Real              *c_tminus_ptr,    /**< output, time value for to the left on coarse grid */
                                     braid_Real              *c_tplus_ptr      /**< output, time value for to the right on coarse grid */
                                     );

/*--------------------------------------------------------------------------
 * PhiStatus Prototypes
 *--------------------------------------------------------------------------*/

/**
 * Initialize a braid_PhiStatus structure 
 **/
braid_Int
_braid_PhiStatusInit(braid_Real       tstart,      /**< current time value  */
                     braid_Real       tplus,       /**< time value to evolve towards, time value to the right */
                     braid_Real       accuracy,    /**< advanced option allowing variable accuracy for implicit phi*/
                     braid_PhiStatus  status       /**< structure to initialize */
                     );

/**
 * Destroy a braid_PhiStatus structure
 **/
braid_Int
_braid_PhiStatusDestroy(braid_PhiStatus  status);        /**< structure to be destroyed */

/**
 * Return the current time value from the PhiStatus structure.
 **/
braid_Int
braid_PhiStatusGetTstart(braid_PhiStatus  status,         /**< structure containing current simulation info */
                         braid_Real      *tstart_ptr      /**< output, current time */
                         );
/**
 * Return the time value to the right of the current time value from
 * the PhiStatus structure.
 **/
braid_Int
braid_PhiStatusGetTplus(braid_PhiStatus  status,         /**< structure containing current simulation info */
                        braid_Real      *tplus_ptr       /**< output, next time value to evolve towards */
                        );

/** 
 * Return the current accuracy value, usually between 0 and 1.0, which can
 * allow for tuning of implicit solve accuracy
 **/
braid_Int
braid_PhiStatusGetAccuracy(braid_PhiStatus  status,         /**< structure containing current simulation info */
                           braid_Real      *accuracy_ptr    /**< output, current accuracy value */
                           );

/** 
 * Set the rfactor, a desired refinement factor for this interval.  rfactor=1
 * indicates no refinement, otherwise, this inteval is subdivided rfactor
 * times. 
 **/
braid_Int
braid_PhiStatusSetRFactor(braid_PhiStatus  status,         /**< structure containing current simulation info */
                          braid_Real       rfactor         /**< user-determined desired rfactor */
                          );

/**
 * Return XBraid status for the current simulation. Two values are 
 * returned, tstart and tplus. 
 *
 * These values are also available through individual Get routines. 
 * These individual routines are the location of detailed documentation on 
 * each parameter, e.g., see *braid_PhiStatusGetTstart* for more information
 * on the *tstart* value.
 **/
braid_Int
braid_PhiStatusGetTstartTplus(braid_PhiStatus  status,         /**< structure containing current simulation info */
                              braid_Real       *tstart_ptr,    /**< output, current time */
                              braid_Real       *tplus_ptr      /**< output, next time value to evolve towards */
                              );



/** @}*/

#ifdef __cplusplus
}
#endif

#endif



