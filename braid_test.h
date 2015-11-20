/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. Email xbraid-support@llnl.gov for support.
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
 

/** \file braid_test.h
 * \brief Define headers for XBraid test routines.
 *
 * This file contains routines used to test a user's XBraid wrapper routines 
 * one-by-one.
 */

#ifndef braid_test_HEADER
#define braid_test_HEADER

#include "braid.h"
#include "_braid.h"

#ifdef __cplusplus
extern "C" {
#endif


/*--------------------------------------------------------------------------
 * Routines for user to test interface routines
 *--------------------------------------------------------------------------*/
/** \defgroup braidtest XBraid test routines
 *  
 *  These are sanity check routines to help a user test their XBraid code.
 *
 *  @{
 */

/**
 * Test the init, access and free functions.\n
 * A vector is initialized at time *t*, written, and then free-d
 **/
braid_Int
braid_TestInitAccess( braid_App              app,     /**< User defined App structure */
                      MPI_Comm               comm_x,  /**< Spatial communicator */
                      FILE                  *fp,      /**< File pointer (could be stdout or stderr) for log messages*/
                      braid_Real             t,       /**< Time value to test init with (used to initialize the vectors)*/
                      braid_PtFcnInit        init,    /**< Initialize a braid_Vector on finest temporal grid*/
                      braid_PtFcnAccess      access,  /**< Allows access to XBraid and current braid_Vector (can be NULL for no writing)*/
                      braid_PtFcnFree        free     /**< Free a braid_Vector*/
                      );

 /**
  * Test the clone function.\n
  * A vector is initialized at time *t*, cloned, and both vectors are written.
  * Then both vectors are free-d.  The user is to check (via the access function) 
  * to see if it is identical.
  **/
braid_Int
braid_TestClone( braid_App              app,         /**< User defined App structure */
                 MPI_Comm               comm_x,      /**< Spatial communicator */
                 FILE                  *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
                 braid_Real             t,           /**< Time value to test clone with  (used to initialize the vectors)*/
                 braid_PtFcnInit        init,        /**< Initialize a braid_Vector on finest temporal grid*/
                 braid_PtFcnAccess      access,      /**< Allows access to XBraid and current braid_Vector (can be NULL for no writing)*/
                 braid_PtFcnFree        free,        /**< Free a braid_Vector*/
                 braid_PtFcnClone       clone        /**< Clone a braid_Vector */
                 );



 /**
  * Test the sum function.\n
  * A vector is initialized at time *t*, cloned, and then these two vectors
  * are summed a few times, with the results written.  The vectors are then free-d.
  * The user is to check (via the access function) that the output matches the 
  * sum of the two original vectors.
  **/
braid_Int
braid_TestSum( braid_App              app,         /**< User defined App structure */
               MPI_Comm               comm_x,      /**< Spatial communicator */
               FILE                  *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
               braid_Real             t,           /**< Time value to test Sum with  (used to initialize the vectors)*/
               braid_PtFcnInit        init,        /**< Initialize a braid_Vector on finest temporal grid*/
               braid_PtFcnAccess      access,      /**< Allows access to XBraid and current braid_Vector (can be NULL for no writing)*/
               braid_PtFcnFree        free,        /**< Free a braid_Vector*/
               braid_PtFcnClone       clone,       /**< Clone a braid_Vector */
               braid_PtFcnSum         sum          /**< Compute vector sum of two braid_Vectors */
               );

/**
 * Test the spatialnorm function.\n
 * A vector is initialized at time *t* and then cloned.  Various
 * norm evaluations like || 3 v || / || v || with known output are then done. 
 *
 * - Returns 0 if the tests fail
 * - Returns 1 if the tests pass
 * - Check the log messages to see details of which tests failed.
 **/
braid_Int
braid_TestSpatialNorm( braid_App              app,         /**< User defined App structure */
                       MPI_Comm               comm_x,      /**< Spatial communicator */
                       FILE                  *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
                       braid_Real             t,           /**< Time value to test SpatialNorm with  (used to initialize the vectors)*/
                       braid_PtFcnInit        init,        /**< Initialize a braid_Vector on finest temporal grid*/
                       braid_PtFcnFree        free,        /**< Free a braid_Vector*/
                       braid_PtFcnClone       clone,       /**< Clone a braid_Vector */
                       braid_PtFcnSum         sum,         /**< Compute vector sum of two braid_Vectors */
                       braid_PtFcnSpatialNorm spatialnorm  /**< Compute norm of a braid_Vector, this is a norm only over space */
                       );
              
/**
 * Test the BufPack, BufUnpack and BufSize functions.\n
 * A vector is initialized at time *t*, packed into a buffer, then unpacked from a buffer.
 * The unpacked result must equal the original vector.  
 *
 * - Returns 0 if the tests fail
 * - Returns 1 if the tests pass
 * - Check the log messages to see details of which tests failed.
 **/
braid_Int
braid_TestBuf( braid_App              app,         /**< User defined App structure */
               MPI_Comm               comm_x,      /**< Spatial communicator */
               FILE                  *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
               braid_Real             t,           /**< Time value to test Buffer routines  (used to initialize the vectors)*/
               braid_PtFcnInit        init,        /**< Initialize a braid_Vector on finest temporal grid*/
               braid_PtFcnFree        free,        /**< Free a braid_Vector*/
               braid_PtFcnSum         sum,         /**< Compute vector sum of two braid_Vectors */
               braid_PtFcnSpatialNorm spatialnorm, /**< Compute norm of a braid_Vector, this is a norm only over space */
               braid_PtFcnBufSize     bufsize,     /**< Computes size in bytes for one braid_Vector MPI buffer */
               braid_PtFcnBufPack     bufpack,     /**< Packs MPI buffer to contain one braid_Vector */
               braid_PtFcnBufUnpack   bufunpack    /**< Unpacks MPI buffer containing one braid_Vector */
               );

/**
 * Test the Coarsen and Refine functions.\n
 * A vector is initialized at time *t*, and various sanity checks on the spatial
 * coarsening and refinement routines are run.
 *
 * - Returns 0 if the tests fail
 * - Returns 1 if the tests pass
 * - Check the log messages to see details of which tests failed.
 **/
braid_Int
braid_TestCoarsenRefine( braid_App                app,         /**< User defined App structure */
                         MPI_Comm                 comm_x,      /**< Spatial communicator */
                         FILE                    *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
                         braid_Real               t,           /**< Time value to initialize test vectors */
                         braid_Real               fdt,         /**< Fine time step value that you spatially coarsen from */
                         braid_Real               cdt,         /**< Coarse time step value that you coarsen to */
                         braid_PtFcnInit          init,        /**< Initialize a braid_Vector on finest temporal grid*/
                         braid_PtFcnAccess        access,      /**< Allows access to XBraid and current braid_Vector (can be NULL for no writing)*/
                         braid_PtFcnFree          free,        /**< Free a braid_Vector*/
                         braid_PtFcnClone         clone,       /**< Clone a braid_Vector */
                         braid_PtFcnSum           sum,         /**< Compute vector sum of two braid_Vectors */
                         braid_PtFcnSpatialNorm   spatialnorm, /**< Compute norm of a braid_Vector, this is a norm only over space */
                         braid_PtFcnSCoarsen      coarsen,     /**< Spatially coarsen a vector */
                         braid_PtFcnSRefine       refine       /**< Spatially refine a vector */
                         );

/**
 * Runs all of the individual braid_Test* routines
 *
 * - Returns 0 if the tests fail
 * - Returns 1 if the tests pass
 * - Check the log messages to see details of which tests failed.
 **/
braid_Int
braid_TestAll( braid_App                app,         /**< User defined App structure */
               MPI_Comm                 comm_x,      /**< Spatial communicator */
               FILE                    *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
               braid_Real               t,           /**< Time value to initialize test vectors with*/
               braid_Real               fdt,         /**< Fine time step value that you spatially coarsen from */
               braid_Real               cdt,         /**< Coarse time step value that you coarsen to */
               braid_PtFcnInit          init,        /**< Initialize a braid_Vector on finest temporal grid*/
               braid_PtFcnFree          free,        /**< Free a braid_Vector*/
               braid_PtFcnClone         clone,       /**< Clone a braid_Vector */
               braid_PtFcnSum           sum,         /**< Compute vector sum of two braid_Vectors */
               braid_PtFcnSpatialNorm   spatialnorm, /**< Compute norm of a braid_Vector, this is a norm only over space */
               braid_PtFcnBufSize       bufsize,     /**< Computes size in bytes for one braid_Vector MPI buffer */
               braid_PtFcnBufPack       bufpack,     /**< Packs MPI buffer to contain one braid_Vector */
               braid_PtFcnBufUnpack     bufunpack,   /**< Unpacks MPI buffer into a braid_Vector */
               braid_PtFcnSCoarsen      coarsen,     /**< Spatially coarsen a vector. If NULL, test is skipped.*/
               braid_PtFcnSRefine       refine       /**< Spatially refine a vector. If NULL, test is skipped.*/
               );

/** @}*/

#ifdef __cplusplus
}
#endif

#endif

