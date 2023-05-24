/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory.
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
 

/** \file braid_test.h
 * \brief Define headers for XBraid user-test routines.
 *
 * This file contains headers for the user to test their XBraid wrapper routines 
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
 * Test the inner_prod function.\n
 * A vector is initialized at time *t1*, then the vector is normalized under the
 * norm induced by inner_prod. A second vector is initialized at time *t2*, and the
 * Gram Schmidt process removes the component of the second vector along the 
 * direction of the first. The test is inconclusive unless both vectors are nonzero 
 * and not orthogonal.
 *
 * - Returns 0 if the tests fail
 * - Returns 1 if the tests pass
 * - Check the log messages to see details of which tests failed.
 **/
braid_Int
braid_TestInnerProd( braid_App              app,        /**< User defined App structure */
                     MPI_Comm               comm_x,     /**< Spatial communicator */  
                     FILE                  *fp,         /**< File pointer (could be stdout or stderr) for log messages */
                     braid_Real             t1,         /**< Time value used to initialize the 1st vector */
                     braid_Real             t2,         /**< Time value used to initialize the 2nd vector (t1 != t2) */
                     braid_PtFcnInit        init,       /**< Initialize a braid_Vector on finest temporal grid */
                     braid_PtFcnFree        free,       /**< Free a braid_Vector */
                     braid_PtFcnSum         sum,        /**< Compute vector sum of two braid_Vectors */ 
                     braid_PtFcnInnerProd   inner_prod  /**< Compute inner product of two braid_Vectors */
                     );
              


/**
 * Test the inner_prod function.\n
 * A vector is initialized at time *t1*, then the vector is normalized under the
 * norm induced by inner_prod. A second vector is initialized at time *t2*, and the
 * Gram Schmidt process removes the component of the second vector along the 
 * direction of the first. The test is inconclusive unless both vectors are nonzero 
 * and not orthogonal.
 *
 * - Returns 0 if the tests fail
 * - Returns 1 if the tests pass
 * - Check the log messages to see details of which tests failed.
 **/
braid_Int
braid_TestInnerProd( braid_App              app,        /**< User defined App structure */
                     MPI_Comm               comm_x,     /**< Spatial communicator */  
                     FILE                  *fp,         /**< File pointer (could be stdout or stderr) for log messages */
                     braid_Real             t1,         /**< Time value used to initialize the 1st vector */
                     braid_Real             t2,         /**< Time value used to initialize the 2nd vector (t1 != t2) */
                     braid_PtFcnInit        init,       /**< Initialize a braid_Vector on finest temporal grid */
                     braid_PtFcnFree        free,       /**< Free a braid_Vector */
                     braid_PtFcnSum         sum,        /**< Compute vector sum of two braid_Vectors */ 
                     braid_PtFcnInnerProd   inner_prod  /**< Compute inner product of two braid_Vectors */
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
 * Test compatibility of the Step and Residual functions.\n
 * A vector is initialized at time *t*, step is called with *dt*,
 * followed by an evaluation of residual, to test the condition
 * fstop - residual( step(u, fstop), u)  approx.  0
 * 
 * - Check the log messages to determine if test passed.  The result 
 *   should approximately be zero.  The more accurate the solution for 
 *   *u* is computed in step, the closer the result will be to 0. 
 * - The residual is also written to file 
 **/
braid_Int
braid_TestResidual( braid_App              app,             /**< User defined App structure */
                    MPI_Comm               comm_x,          /**< Spatial communicator */
                    FILE                   *fp,             /**< File pointer (could be stdout or stderr) for log messages*/
                    braid_Real             t,               /**< Time value to initialize test vectors */
                    braid_Real             dt,              /**< Time step value to use in step */ 
                    braid_PtFcnInit        myinit,          /**< Initialize a braid_Vector on finest temporal grid*/
                    braid_PtFcnAccess      myaccess,        /**< Allows access to XBraid and current braid_Vector (can be NULL for no writing)*/ 
                    braid_PtFcnFree        myfree,          /**< Free a braid_Vector*/
                    braid_PtFcnClone       clone,           /**< Clone a braid_Vector */
                    braid_PtFcnSum         sum,             /**< Compute vector sum of two braid_Vectors */
                    braid_PtFcnSpatialNorm spatialnorm,     /**< Compute norm of a braid_Vector, this is a norm only over space */
                    braid_PtFcnResidual    residual,        /**< Compute a residual given two consectuive braid_Vectors */
                    braid_PtFcnStep        step             /**< Compute a time step with a braid_Vector */
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
               braid_PtFcnSRefine       refine,      /**< Spatially refine a vector. If NULL, test is skipped.*/
               braid_PtFcnResidual      residual,    /**< Compute a residual given two consectuive braid_Vectors */
               braid_PtFcnStep          step         /**< Compute a time step with a braid_Vector */
               );


/**
 * Test functions required for Delta correction.
 * Initializes a braid_Vector and braid_Basis at time 0, then
 * tests the inner product function with braid_TestInnerProd,
 * then checks that the basis vectors are not linearly dependent
 * with the Gram Schmidt process. Finally, compares the user
 * propagation of tangent vectors against a finite difference
 * approximation:
 * [step_du(u)] Psi_i - (step(u + eps Psi_i) - step(u)/eps ~= 0 
 * 
 */
braid_Int
braid_TestDelta(braid_App               app,          /**< User defined App structure */
                MPI_Comm                comm_x,       /**< Spatial communicator */
                FILE                   *fp,           /**< File pointer (could be stdout or stderr) for log messages*/
                braid_Real              t,            /**< Time value to initialize test vectors with*/
                braid_Real              dt,           /**< time step size */
                braid_Int               rank,         /**< rank (number of columns) of basis */
                braid_PtFcnInit         myinit,       /**< Initialize a braid_Vector*/
                braid_PtFcnInitBasis    myinit_basis, /**< Initialize the ith column of a basis set of braid_Vectors */
                braid_PtFcnAccess       myaccess,     /**< Allows access to XBraid and current braid_Vector and braid_Basis */ 
                braid_PtFcnFree         myfree,       /**< Free a braid_Vector*/
                braid_PtFcnClone        myclone,      /**< Clone a braid_Vector */
                braid_PtFcnSum          mysum,        /**< Compute vector sum of two braid_Vectors */
                braid_PtFcnBufSize      bufsize,      /**< Computes size in bytes for one braid_Vector MPI buffer */
                braid_PtFcnBufPack      bufpack,      /**< Packs MPI buffer to contain one braid_Vector */
                braid_PtFcnBufUnpack    bufunpack,    /**< Unpacks MPI buffer containing one braid_Vector */
                braid_PtFcnInnerProd    myinner_prod, /**< Compute inner product of two braid_Vectors */
                braid_PtFcnStep         mystep        /**< Compute a time step with a braid_Vector */
                );

/**
 * Warmup calls all user functions once with sensible arguments.
 * This is especially useful for the Julia wrapper, since Julia is JIT compiled,
 * and this force compiles all the functions used in braid_Drive.
 */
braid_Int
braid_Warmup( braid_App                app,         /**< User defined App structure */
              MPI_Comm                 comm_x,      /**< Spatial communicator */
              braid_Real               t,           /**< Time value to initialize test vectors with*/
              braid_Real               fdt,         /**< Fine time step value that you spatially coarsen from */
              braid_Real               cdt,         /**< Coarse time step value that you coarsen to */
              braid_PtFcnInit          init,        /**< Initialize a braid_Vector on finest temporal grid*/
              braid_PtFcnAccess        access,      /**< Allows access to XBraid and current braid_Vector (can be NULL for no writing)*/
              braid_PtFcnFree          free,        /**< Free a braid_Vector*/
              braid_PtFcnClone         clone,       /**< Clone a braid_Vector */
              braid_PtFcnSum           sum,         /**< Compute vector sum of two braid_Vectors */
              braid_PtFcnSpatialNorm   spatialnorm, /**< Compute norm of a braid_Vector, this is a norm only over space */
              braid_PtFcnBufSize       bufsize,     /**< Computes size in bytes for one braid_Vector MPI buffer */
              braid_PtFcnBufPack       bufpack,     /**< Packs MPI buffer to contain one braid_Vector */
              braid_PtFcnBufUnpack     bufunpack,   /**< Unpacks MPI buffer into a braid_Vector */
              braid_PtFcnSCoarsen      coarsen,     /**< Spatially coarsen a vector. If NULL, skipped.*/
              braid_PtFcnSRefine       refine,      /**< Spatially refine a vector. If NULL, skipped.*/
              braid_PtFcnStep          step,        /**< Compute a time step with a braid_Vector */
              braid_PtFcnInitBasis     init_basis,  /**< Initialize a basis vector at time t and spatial index i. If NULL, skipped */
              braid_PtFcnInnerProd     innerprod    /**< Compute the inner product between two braid_Vectors. If NULL, skipped. */
              );
/** @}*/

#ifdef __cplusplus
}
#endif

#endif

