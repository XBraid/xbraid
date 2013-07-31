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

#ifndef warp_HEADER
#define warp_HEADER

#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

typedef int    warp_Int;
typedef double warp_Float;

/*--------------------------------------------------------------------------
 * User-written routines
 *--------------------------------------------------------------------------*/

struct _warp_App_struct;
/**
 * Blah...
 **/
typedef struct _warp_App_struct *warp_App;

struct _warp_Vector_struct;
/**
 * Blah...
 **/
typedef struct _warp_Vector_struct *warp_Vector;

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnPhi)(warp_App      app,
                 warp_Float    tstart,
                 warp_Float    tstop,
                 warp_Int      gzero,
                 warp_Vector   u,
                 warp_Int     *rfactor_ptr);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnInit)(warp_App      app,
                  warp_Float    t,
                  warp_Vector  *u_ptr);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnClone)(warp_App      app,
                   warp_Vector   u,
                   warp_Vector  *v_ptr);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnFree)(warp_App     app,
                  warp_Vector  u);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnSum)(warp_App     app,
                 warp_Float   alpha,
                 warp_Vector  x,
                 warp_Float   beta,
                 warp_Vector  y);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnDot)(warp_App      app,
                 warp_Vector   u,
                 warp_Vector   v,
                 warp_Float   *dot_ptr);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnWrite)(warp_App      app,
                   warp_Float    t,
                   warp_Vector   u);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnBufSize)(warp_App   app,
                     warp_Int  *size_ptr);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnBufPack)(warp_App      app,
                     warp_Vector   u,
                     void         *buffer);

/**
 * Blah...
 **/
typedef warp_Int
(*warp_PtFcnBufUnpack)(warp_App      app,
                       void         *buffer,
                       warp_Vector  *u_ptr);

/**
 * (Optional) Blah...
 **/
typedef warp_Int
(*warp_PtFcnCoarsen)(warp_App      app,
                     warp_Float    tstart,
                     warp_Float    tstop,
                     warp_Vector   fu,
                     warp_Vector  *cu_ptr);

/**
 * (Optional) Blah...
 **/
typedef warp_Int
(*warp_PtFcnRefine)(warp_App      app,
                    warp_Float    tstart,
                    warp_Float    tstop,
                    warp_Vector   cu,
                    warp_Vector  *fu_ptr);

/*--------------------------------------------------------------------------
 * User interface routines
 *--------------------------------------------------------------------------*/

struct _warp_Core_struct;
/**
 * Blah...
 **/
typedef struct _warp_Core_struct *warp_Core;

/**
 * Create a core object with the required initial data.
 **/
warp_Int
warp_Init(MPI_Comm              comm_world,
          MPI_Comm              comm,
          warp_Float            tstart,
          warp_Float            tstop,
          warp_Int              ntime,
          warp_App              app,
          warp_PtFcnPhi         phi,
          warp_PtFcnInit        init,
          warp_PtFcnClone       clone,
          warp_PtFcnFree        free,
          warp_PtFcnSum         sum,
          warp_PtFcnDot         dot,
          warp_PtFcnWrite       write,
          warp_PtFcnBufSize     bufsize,
          warp_PtFcnBufPack     bufpack,
          warp_PtFcnBufUnpack   bufunpack,
          warp_Core            *core_ptr);

/**
 * Integrate in time.
 **/
warp_Int
warp_Drive(warp_Core  core);

/**
 * Destroy core.
 **/
warp_Int
warp_Destroy(warp_Core  core);

/**
 * Print statistics.
 **/
warp_Int
warp_PrintStats(warp_Core  core);

/**
 * Set max number of multigrid levels.
 **/
warp_Int
warp_SetMaxLevels(warp_Core  core,
                  warp_Int   max_levels);

/**
 * Set number of multigrid relaxation sweeps.
 **/
warp_Int
warp_SetNRelax(warp_Core  core,
               warp_Int   nrelax);

/**
 * Set absolute stopping tolerance.
 **/
warp_Int
warp_SetTol(warp_Core   core,
            warp_Float  tol);

/**
 * Set the coarsening factor {\tt cfactor} on grid level {\tt level} (level 0 is
 * the finest grid).  The default factor is 2 on all levels.  To change the
 * default factor, use {\tt level} = -1.
 **/
warp_Int
warp_SetCFactor(warp_Core  core,
                warp_Int   level,
                warp_Int   cfactor);

/**
 * Set max number of multigrid iterations.
 **/
warp_Int
warp_SetMaxIter(warp_Core  core,
                warp_Int   max_iter);

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif

