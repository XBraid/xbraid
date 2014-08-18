
/*BHEADER**********************************************************************
 * Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of WARP.  See file COPYRIGHT for details.
 *
 * WARP is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 ***********************************************************************EHEADER*/

/** \file warp_test.h
 * \brief Define headers for Warp test routines.
 *
 * This file contains routines used to test a user's Warp wrapper routines 
 * one-by-one.
 */

#ifndef warp_test_HEADER
#define warp_test_HEADER

#include "warp.h"
#include "_warp.h"

#ifdef __cplusplus
extern "C" {
#endif


/*--------------------------------------------------------------------------
 * Routines for user to test interface routines
 *--------------------------------------------------------------------------*/
/** \defgroup warptest Warp test routines
 *  
 *  These are sanity check routines to help a user test their Warp code.
 *
 *  @{
 */

/**
 * Test the init, write and free functions.\n
 * A vector is initialized at time *t*, written, and then free-d
 **/
warp_Int
warp_TestInitWrite( warp_App              app,     /**< User defined App structure */
                    MPI_Comm              comm_x,  /**< Spatial communicator */
                    FILE                 *fp,      /**< File pointer (could be stdout or stderr) for log messages*/
                    warp_Real             t,       /**< Time value to test init with (used to initialize the vectors)*/
                    warp_PtFcnInit        init,    /**< Initialize a warp_Vector on finest temporal grid*/
                    warp_PtFcnWrite       write,   /**< Write a warp_Vector (can be NULL for no writing)*/
                    warp_PtFcnFree        free     /**< Free a warp_Vector*/
                    );

 /**
  * Test the clone function.\n
  * A vector is initialized at time *t*, cloned, and both vectors are written.
  * Then both vectors are free-d.  The user is to check (via the write function) 
  * to see if it is identical.
  **/
warp_Int
warp_TestClone( warp_App              app,         /**< User defined App structure */
                MPI_Comm              comm_x,      /**< Spatial communicator */
                FILE                 *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
                warp_Real             t,           /**< Time value to test clone with  (used to initialize the vectors)*/
                warp_PtFcnInit        init,        /**< Initialize a warp_Vector on finest temporal grid*/
                warp_PtFcnWrite       write,       /**< Write a warp_Vector (can be NULL for no writing)*/
                warp_PtFcnFree        free,        /**< Free a warp_Vector*/
                warp_PtFcnClone       clone        /**< Clone a warp_Vector */
                );



 /**
  * Test the sum function.\n
  * A vector is initialized at time *t*, cloned, and then these two vectors
  * are summed a few times, with the results written.  The vectors are then free-d.
  * The user is to check (via the write function) that the output matches the 
  * sum of the two original vectors.
  **/
warp_Int
warp_TestSum( warp_App              app,         /**< User defined App structure */
              MPI_Comm              comm_x,      /**< Spatial communicator */
              FILE                 *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
              warp_Real             t,           /**< Time value to test Sum with  (used to initialize the vectors)*/
              warp_PtFcnInit        init,        /**< Initialize a warp_Vector on finest temporal grid*/
              warp_PtFcnWrite       write,       /**< Write a warp_Vector (can be NULL for no writing)*/
              warp_PtFcnFree        free,        /**< Free a warp_Vector*/
              warp_PtFcnClone       clone,       /**< Clone a warp_Vector */
              warp_PtFcnSum         sum          /**< Compute vector sum of two warp_Vectors */
              );

/**
 * Test the dot function.\n
 * A vector is initialized at time *t* and then cloned.  Various
 * dot products like <3 v, v>/<v, v> are computed with known output, e.g.,
 * <3 v, v>/<v, v> must equal 3.  If all the tests pass, then 1 is returned.
 *
 * - Returns 0 if the tests fail
 * - Returns 1 if the tests pass
 * - Check the log messages to see details of which tests failed.
 **/
warp_Int
warp_TestDot( warp_App              app,         /**< User defined App structure */
              MPI_Comm              comm_x,      /**< Spatial communicator */
              FILE                 *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
              warp_Real             t,           /**< Time value to test Dot with  (used to initialize the vectors)*/
              warp_PtFcnInit        init,        /**< Initialize a warp_Vector on finest temporal grid*/
              warp_PtFcnFree        free,        /**< Free a warp_Vector*/
              warp_PtFcnClone       clone,       /**< Clone a warp_Vector */
              warp_PtFcnSum         sum,         /**< Compute vector sum of two warp_Vectors */
              warp_PtFcnDot         dot          /**< Compute dot product of two warp_Vectors */
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
warp_Int
warp_TestBuf( warp_App              app,         /**< User defined App structure */
              MPI_Comm              comm_x,      /**< Spatial communicator */
              FILE                 *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
              warp_Real             t,           /**< Time value to test Buffer routines  (used to initialize the vectors)*/
              warp_PtFcnInit        init,        /**< Initialize a warp_Vector on finest temporal grid*/
              warp_PtFcnFree        free,        /**< Free a warp_Vector*/
              warp_PtFcnSum         sum,         /**< Compute vector sum of two warp_Vectors */
              warp_PtFcnDot         dot,         /**< Compute dot product of two warp_Vectors */
              warp_PtFcnBufSize     bufsize,     /**< Computes size in bytes for one warp_Vector MPI buffer */
              warp_PtFcnBufPack     bufpack,     /**< Packs MPI buffer to contain one warp_Vector */
              warp_PtFcnBufUnpack   bufunpack    /**< Unpacks MPI buffer containing one warp_Vector */
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
warp_Int
warp_TestCoarsenRefine( warp_App          app,         /**< User defined App structure */
                        MPI_Comm          comm_x,      /**< Spatial communicator */
                        FILE             *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
                        warp_Real         t,           /**< Time value to initialize test vectors */
                        warp_Real         fdt,         /**< Fine time step value that you spatially coarsen from */
                        warp_Real         cdt,         /**< Coarse time step value that you coarsen to */
                        warp_PtFcnInit    init,        /**< Initialize a warp_Vector on finest temporal grid*/
                        warp_PtFcnWrite   write,       /**< Write a warp_Vector (can be NULL for no writing)*/
                        warp_PtFcnFree    free,        /**< Free a warp_Vector*/
                        warp_PtFcnClone   clone,       /**< Clone a warp_Vector */
                        warp_PtFcnSum     sum,         /**< Compute vector sum of two warp_Vectors */
                        warp_PtFcnDot     dot,         /**< Compute dot product of two warp_Vectors */
                        warp_PtFcnCoarsen coarsen,     /**< Spatially coarsen a vector */
                        warp_PtFcnRefine  refine       /**< Spatially refine a vector */
                        );
/**
 * Runs all of the individual warp_Test* routines
 *
 * - Returns 0 if the tests fail
 * - Returns 1 if the tests pass
 * - Check the log messages to see details of which tests failed.
 **/
warp_Int
warp_TestAll( warp_App             app,         /**< User defined App structure */
              MPI_Comm             comm_x,      /**< Spatial communicator */
              FILE                *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
              warp_Real            t,           /**< Time value to initialize test vectors with*/
              warp_Real            fdt,         /**< Fine time step value that you spatially coarsen from */
              warp_Real            cdt,         /**< Coarse time step value that you coarsen to */
              warp_PtFcnInit       init,        /**< Initialize a warp_Vector on finest temporal grid*/
              warp_PtFcnFree       free,        /**< Free a warp_Vector*/
              warp_PtFcnClone      clone,       /**< Clone a warp_Vector */
              warp_PtFcnSum        sum,         /**< Compute vector sum of two warp_Vectors */
              warp_PtFcnDot        dot,         /**< Compute dot product of two warp_Vectors */
              warp_PtFcnBufSize    bufsize,     /**< Computes size in bytes for one warp_Vector MPI buffer */
              warp_PtFcnBufPack    bufpack,     /**< Packs MPI buffer to contain one warp_Vector */
              warp_PtFcnBufUnpack  bufunpack,   /**< Unpacks MPI buffer into a warp_Vector */
              warp_PtFcnCoarsen    coarsen,     /**< Spatially coarsen a vector. If NULL, test is skipped.*/
              warp_PtFcnRefine     refine       /**< Spatially refine a vector. If NULL, test is skipped.*/
              );

/** @}*/

#ifdef __cplusplus
}
#endif

#endif

