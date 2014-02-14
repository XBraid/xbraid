
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

/** \file warp.h
 * \brief Define headers for user interface routines.
 *
 * This file contains routines used to allow the user to initialize, run
 * and get and set warp. 
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


warp_Int
warp_TestInitWrite( warp_App              app,     /**< User defined App structure */
                    MPI_Comm              comm_x,  /**< Spatial communicator */
                    FILE                 *fp,      /**< File pointer (could be stdout or stderr) for log messages*/
                    warp_Real             t,       /**< Time value to test init with */
                    warp_PtFcnInit        init,    /**< Initialize a warp_Vector function on finest temporal grid*/
                    warp_PtFcnWrite       write,   /**< Write temporal state warp_Vector (can be NULL for no writing)*/
                    warp_PtFcnFree        free     /**< Free a temporal state warp_Vector*/
                    );
warp_Int
warp_TestClone( warp_App              app,         /**< User defined App structure */
                MPI_Comm              comm_x,      /**< Spatial communicator */
                FILE                 *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
                warp_Real             t,           /**< Time value to test clone with  */
                warp_PtFcnInit        init,        /**< Initialize a warp_Vector function on finest temporal grid*/
                warp_PtFcnWrite       write,       /**< Write temporal state warp_Vector (can be NULL for no writing)*/
                warp_PtFcnFree        free,        /**< Free a temporal state warp_Vector*/
                warp_PtFcnClone       clone        /**< Clone a temporal state warp_Vector */
                );



warp_Int
warp_TestSum( warp_App              app,         /**< User defined App structure */
              MPI_Comm              comm_x,      /**< Spatial communicator */
              FILE                 *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
              warp_Real             t,           /**< Time value to test Sum with  (used to initialize the vectors*/
              warp_PtFcnInit        init,        /**< Initialize a warp_Vector function on finest temporal grid*/
              warp_PtFcnWrite       write,       /**< Write temporal state warp_Vector (can be NULL for no writing)*/
              warp_PtFcnFree        free,        /**< Free a temporal state warp_Vector*/
              warp_PtFcnClone       clone,       /**< Clone a temporal state warp_Vector */
              warp_PtFcnSum         sum          /**< Compute vector sum of two temporal states*/
              );

warp_Int
warp_TestDot( warp_App              app,         /**< User defined App structure */
              MPI_Comm              comm_x,      /**< Spatial communicator */
              FILE                 *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
              warp_Real             t,           /**< Time value to test Dot with  (used to initialize the vectors*/
              warp_PtFcnInit        init,        /**< Initialize a warp_Vector function on finest temporal grid*/
              warp_PtFcnFree        free,        /**< Free a temporal state warp_Vector*/
              warp_PtFcnClone       clone,       /**< Clone a temporal state warp_Vector */
              warp_PtFcnSum         sum,         /**< Compute vector sum of two temporal states*/
              warp_PtFcnDot         dot,         /**< Compute dot product of two temporal states*/
              warp_Int             *correct      /**< Boolean describing whether all the tests passed*/
              );
              
warp_Int
warp_TestBuf( warp_App              app,         /**< User defined App structure */
              MPI_Comm              comm_x,      /**< Spatial communicator */
              FILE                 *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
              warp_Real             t,           /**< Time value to test Buffer routines  (used to initialize the vectors*/
              warp_PtFcnInit        init,        /**< Initialize a warp_Vector function on finest temporal grid*/
              warp_PtFcnFree        free,        /**< Free a temporal state warp_Vector*/
              warp_PtFcnSum         sum,         /**< Compute vector sum of two temporal states*/
              warp_PtFcnDot         dot,         /**< Compute dot product of two temporal states*/
              warp_PtFcnBufSize     bufsize,     /**< Computes size for MPI buffer for one */
              warp_PtFcnBufPack     bufpack,     /**< Packs MPI buffer to contain one temporal state*/
              warp_PtFcnBufUnpack   bufunpack,   /**< Unpacks MPI buffer containing one temporal state*/
              warp_Int             *correct      /**< Boolean describing whether all the tests passed*/
              );

warp_Int
warp_TestCoarsenRefine( warp_App          app,         /**< User defined App structure */
                        MPI_Comm          comm_x,      /**< Spatial communicator */
                        FILE             *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
                        warp_Real         t,           /**< Time value to initialize vector */
                        warp_Real         f_tminus,    /**< Fine time value before t, used by coarsen and 
                                                            refine to determine appropriate coarsening and 
                                                            refinement algorithm */
                        warp_Real         f_tplus,     /**< Fine time value after t, used by coarsen and 
                                                            refine to determine appropriate coarsening and 
                                                            refinement algorithm, e.g. (f_tplus -t) equals 
                                                            the fine delta t*/
                        warp_Real         c_tminus,    /**< Coarse time value before t, used by coarsen 
                                                            and refine to determine appropriate coarsening 
                                                            and refinement algorithm */
                        warp_Real         c_tplus,     /**< Coarse time value after t, used by coarsen and 
                                                            refine to determine appropriate coarsening and 
                                                            refinement algorithm, e.g. (c_tplus -t) equals 
                                                            the coarse delta t */
                        warp_PtFcnInit    init,        /**< Initialize a warp_Vector function on finest temporal grid*/
                        warp_PtFcnWrite   write,       /**< Write temporal state warp_Vector (can be NULL for no writing)*/
                        warp_PtFcnFree    free,        /**< Free a temporal state warp_Vector*/
                        warp_PtFcnClone   clone,       /**< Clone a temporal state warp_Vector */
                        warp_PtFcnSum     sum,         /**< Compute vector sum of two temporal states*/
                        warp_PtFcnDot     dot,         /**< Compute dot product of two temporal states*/
                        warp_PtFcnCoarsen coarsen,     /**< Spatially coarsen a vector */
                        warp_PtFcnRefine  refine,      /**< Spatially refine a vector */
                        warp_Int         *correct      /**< Boolean describing whether all the tests passed*/
                        );
warp_Int
warp_TestAll( warp_App             app,         /**< User defined App structure */
              MPI_Comm             comm_x,      /**< Spatial communicator */
              FILE                *fp,          /**< File pointer (could be stdout or stderr) for log messages*/
              warp_Real            t,           /**< Time value to initialize vector */
              warp_Real            f_tminus,    /**< Fine time value before t.  If t==0, this may also
                                                     be 0.0 */ 
              warp_Real            f_tplus,     /**< Fine time value after t.  If t==0, this may also
                                                     be 0.0.  Also, f_tplus -t equals the fine delta t*/
              warp_Real            c_tminus,    /**< Coarse time value before t. If t==0, this may also
                                                     be 0.0 */ 
              warp_Real            c_tplus,     /**< Coarse time value after t. If t==0, this may also
                                                     be 0.0 */
              warp_PtFcnInit       init,        /**< Initialize a warp_Vector function on finest temporal grid*/
              warp_PtFcnFree       free,        /**< Free a temporal state warp_Vector*/
              warp_PtFcnClone      clone,       /**< Clone a temporal state warp_Vector */
              warp_PtFcnSum        sum,         /**< Compute vector sum of two temporal states*/
              warp_PtFcnDot        dot,         /**< Compute dot product of two temporal states*/
              warp_PtFcnBufSize    bufsize,     /**< Computes size for MPI buffer for one */
              warp_PtFcnBufPack    bufpack,     /**< Packs MPI buffer to contain one temporal state*/
              warp_PtFcnBufUnpack  bufunpack,   /**< Unpacks MPI buffer containing one temporal state*/
              warp_PtFcnCoarsen    coarsen,     /**< Spatially coarsen a vector. If null, test is skipped.*/
              warp_PtFcnRefine     refine,      /**< Spatially refine a vector. If null, test is skipped.*/
              warp_Int            *correct      /**< Boolean describing whether all the tests passed*/
              );


#ifdef __cplusplus
}
#endif

#endif

