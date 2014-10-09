// Copyright (c) 2013, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. Written by
// Jacob Schroder schroder2@llnl.gov, Rob Falgout falgout2@llnl.gov,
// Tzanio Kolev kolev1@llnl.gov, Ulrike Yang yang11@llnl.gov,
// Veselin Dobrev dobrev1@llnl.gov, et al.
// LLNL-CODE-660355. All rights reserved.
//
// This file is part of XBraid. Email schroder2@llnl.gov on how to download.
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License (as published by the Free Software
// Foundation) version 2.1 dated February 1999.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 59
// Temple Place, Suite 330, Boston, MA 02111-1307 USA
//

#ifndef braid_hpp_HEADER
#define braid_hpp_HEADER

#include "_braid.h"
#include "braid.h"
#include "braid_test.h"

class BraidAccessStatus;
class BraidPhiStatus;
class BraidCoarsenRefStatus;

// Wrapper for BRAID's App object. Users should inherit this class and implement
// the purely virtual functions (see braid.h for descriptions).
class BraidApp
{
public:
   // User must set these four values as they are used by Braid.
   // This is done by the BraidApp constructor.
   MPI_Comm     comm_t;
   braid_Real   tstart;
   braid_Real   tstop;
   braid_Int    ntime;

   BraidApp(MPI_Comm   _comm_t,
            braid_Real _tstart = 0.0,
            braid_Real _tstop  = 1.0,
            braid_Int  _ntime  = 100)
      : comm_t(_comm_t), tstart(_tstart), tstop(_tstop), ntime(_ntime) { }

   virtual ~BraidApp() { }

   virtual braid_Int Phi(braid_Vector    _u,
                         BraidPhiStatus &pstatus) = 0;

   virtual braid_Int Clone(braid_Vector  _u,
                           braid_Vector *v_ptr) = 0;

   virtual braid_Int Init(braid_Real    t,
                          braid_Vector *u_ptr) = 0;

   virtual braid_Int Free(braid_Vector _u) = 0;

   virtual braid_Int Sum(braid_Real   alpha,
                         braid_Vector _x,
                         braid_Real   beta,
                         braid_Vector _y) = 0;

   virtual braid_Int SpatialNorm(braid_Vector  _u,
                                 braid_Real   *norm_ptr) = 0;

   virtual braid_Int Access(braid_Vector       _u,
                            BraidAccessStatus &astatus) = 0;

   virtual braid_Int BufSize(braid_Int *size_ptr) = 0;

   virtual braid_Int BufPack(braid_Vector  _u,
                             void         *buffer,
                             braid_Int    *size_ptr) = 0;

   virtual braid_Int BufUnpack(void         *buffer,
                               braid_Vector *u_ptr) = 0;

   // These two functions may be optionally defined by the user, if spatial
   // coarsening is desired (see documentation for more details).  To turn on
   // spatial coarsening, use core.SetSpatialCoarsenAndRefine()
   virtual braid_Int Coarsen(braid_Vector           _fu,
                             braid_Vector          *cu_ptr,
                             BraidCoarsenRefStatus &status)
   {
      fprintf(stderr, "Braid C++ Wrapper Warning: turn off spatial coarsening until Coarsen and Refine have been user implemented\n");
      Clone(_fu, cu_ptr);
      return 0;
   }

   virtual braid_Int Refine(braid_Vector           _cu,
                            braid_Vector          *fu_ptr,
                            BraidCoarsenRefStatus &status)
   {
      fprintf(stderr, "Braid C++ Wrapper Warning: turn off spatial coarsening until Coarsen and Refine have been user implemented\n");
      Clone(_cu, fu_ptr);
      return 0;
   }

};

// Wrapper for BRAID's AccessStatus object
class BraidAccessStatus
{
   private:
      braid_AccessStatus astatus;

   public:
      BraidAccessStatus(braid_AccessStatus _astatus)
      {
         astatus = _astatus;
      }

      void GetTILD(braid_Real *t_ptr, braid_Int *iter_ptr, braid_Int *level_ptr, braid_Int *done_ptr) { braid_AccessStatusGetTILD(astatus, t_ptr, iter_ptr, level_ptr, done_ptr); }
      void GetT(braid_Real *t_ptr)              { braid_AccessStatusGetT(astatus, t_ptr); }
      void GetDone(braid_Int *done_ptr)         { braid_AccessStatusGetDone(astatus, done_ptr); }
      void GetLevel(braid_Int *level_ptr)       { braid_AccessStatusGetLevel(astatus, level_ptr); }
      void GetIter(braid_Int *iter_ptr)         { braid_AccessStatusGetIter(astatus, iter_ptr); }
      void GetWrapperTest(braid_Int *wtest_ptr) { braid_AccessStatusGetWrapperTest(astatus, wtest_ptr); }
      void GetResidual(braid_Real *rnorm_ptr)   { braid_AccessStatusGetResidual(astatus, rnorm_ptr); }

      // The braid_AccessStatus structure is deallocated inside of Braid
      // This class is just to make code consistently look object oriented
      ~BraidAccessStatus() { }
};

// Wrapper for BRAID's PhiStatus object
class BraidPhiStatus
{
   private:
      braid_PhiStatus pstatus;

   public:
      BraidPhiStatus(braid_PhiStatus _pstatus)
      {
         pstatus = _pstatus;
      }

      void GetTstartTstop(braid_Real *tstart_ptr, braid_Real *tstop_ptr) { braid_PhiStatusGetTstartTstop(pstatus, tstart_ptr, tstop_ptr); }
      void GetTstart(braid_Real *tstart_ptr)     { braid_PhiStatusGetTstart(pstatus, tstart_ptr); }
      void GetTstop(braid_Real *tstop_ptr)       { braid_PhiStatusGetTstop(pstatus, tstop_ptr); }
      void GetAccuracy(braid_Real *accuracy_ptr)  { braid_PhiStatusGetAccuracy(pstatus, accuracy_ptr); }
      void SetRFactor(braid_Int rfactor)         { braid_PhiStatusSetRFactor(pstatus, rfactor); }

      // The braid_PhiStatus structure is deallocated inside of Braid
      // This class is just to make code consistently look object oriented
      ~BraidPhiStatus() { }
};

// Wrapper for BRAID's CoarsenRefStatus object
class BraidCoarsenRefStatus
{
   private:
      braid_CoarsenRefStatus cstatus;

   public:
      BraidCoarsenRefStatus(braid_CoarsenRefStatus  _cstatus)
      {
         cstatus = _cstatus;
      }

      void GetTpriorTstop(braid_Real *tstart_ptr, braid_Real *f_tprior_ptr, braid_Real *f_tstop_ptr, braid_Real *c_tprior_ptr, braid_Real *c_tstop_ptr) { braid_CoarsenRefStatusGetTpriorTstop(cstatus, tstart_ptr, f_tprior_ptr, f_tstop_ptr, c_tprior_ptr, c_tstop_ptr); }
      void GetTstart(braid_Real *tstart_ptr)    { braid_CoarsenRefStatusGetTstart(cstatus, tstart_ptr); }
      void GetFTstop(braid_Real *f_tstop_ptr)   { braid_CoarsenRefStatusGetFTstop(cstatus, f_tstop_ptr); }
      void GetFTprior(braid_Real *f_tprior_ptr) { braid_CoarsenRefStatusGetFTprior(cstatus, f_tprior_ptr); }
      void GetCTstop(braid_Real *c_tstop_ptr)   { braid_CoarsenRefStatusGetCTstop(cstatus, c_tstop_ptr); }
      void GetCTprior(braid_Real *c_tprior_ptr) { braid_CoarsenRefStatusGetCTprior(cstatus, c_tprior_ptr); }

      // The braid_CoarsenRefStatus structure is deallocated inside of Braid
      // This class is just to make code consistently look object oriented
      ~BraidCoarsenRefStatus() { }
};

// Static functions passed to Braid, with braid_App == BraidApp*
static braid_Int _BraidAppPhi(braid_App       _app,
                              braid_Vector    _u,
                              braid_PhiStatus _pstatus)
{
   BraidApp *app = (BraidApp*)_app;
   BraidPhiStatus pstatus = BraidPhiStatus(_pstatus);
   return app -> Phi(_u, pstatus);
}

static braid_Int _BraidAppClone(braid_App     _app,
                                braid_Vector  _u,
                                braid_Vector *v_ptr)
{
   BraidApp *app = (BraidApp*)_app;
   return app -> Clone(_u, v_ptr);
}

static braid_Int _BraidAppInit(braid_App     _app,
                               braid_Real    t,
                               braid_Vector *u_ptr)
{
   BraidApp *app = (BraidApp*)_app;
   return app -> Init(t, u_ptr);
}

static braid_Int _BraidAppFree(braid_App    _app,
                               braid_Vector _u)
{
   BraidApp *app = (BraidApp*)_app;
   return app -> Free(_u);
}

static braid_Int _BraidAppSum(braid_App    _app,
                              braid_Real   alpha,
                              braid_Vector _x,
                              braid_Real   beta,
                              braid_Vector _y)
{
   BraidApp *app = (BraidApp*)_app;
   return app -> Sum(alpha, _x, beta, _y);
}

static braid_Int _BraidAppSpatialNorm(braid_App     _app,
                                      braid_Vector  _u,
                                      braid_Real   *norm_ptr)
{
   BraidApp *app = (BraidApp*)_app;
   return app -> SpatialNorm(_u, norm_ptr);
}

static braid_Int _BraidAppAccess(braid_App          _app,
                                 braid_Vector       _u,
                                 braid_AccessStatus _astatus)
{
   BraidApp *app = (BraidApp*)_app;
   BraidAccessStatus astatus = BraidAccessStatus(_astatus);
   return app -> Access(_u, astatus);
}

static braid_Int _BraidAppBufSize(braid_App  _app,
                                  braid_Int *size_ptr)
{
   BraidApp *app = (BraidApp*)_app;
   return app -> BufSize(size_ptr);
}

static braid_Int _BraidAppBufPack(braid_App     _app,
                                  braid_Vector  _u,
                                  void         *buffer,
                                  braid_Int    *size_ptr)
{
   BraidApp *app = (BraidApp*)_app;
   return app -> BufPack(_u, buffer, size_ptr);
}

static braid_Int _BraidAppBufUnpack(braid_App     _app,
                                    void         *buffer,
                                    braid_Vector *u_ptr)
{
   BraidApp *app = (BraidApp*)_app;
   return app -> BufUnpack(buffer, u_ptr);
}

static braid_Int _BraidAppCoarsen(braid_App               _app,
                                  braid_Vector            _fu,
                                  braid_Vector           *cu_ptr,
                                  braid_CoarsenRefStatus  _cstatus)
{
   BraidApp *app = (BraidApp*)_app;
   BraidCoarsenRefStatus cstatus = BraidCoarsenRefStatus(_cstatus);
   return app -> Coarsen(_fu, cu_ptr, cstatus);
}

static braid_Int _BraidAppRefine(braid_App               _app,
                                 braid_Vector            _fu,
                                 braid_Vector           *cu_ptr,
                                 braid_CoarsenRefStatus  _cstatus)
{
   BraidApp *app = (BraidApp*)_app;
   BraidCoarsenRefStatus cstatus = BraidCoarsenRefStatus(_cstatus);
   return app -> Refine(_fu, cu_ptr, cstatus);
}

// Wrapper for BRAID's core object
class BraidCore
{
private:
   braid_Core core;

public:
   BraidCore(MPI_Comm comm_world, BraidApp *app)
   {
      braid_Init(comm_world,
                 app->comm_t, app->tstart, app->tstop, app->ntime, (braid_App)app,
                 _BraidAppPhi, _BraidAppInit, _BraidAppClone, _BraidAppFree,
                 _BraidAppSum, _BraidAppSpatialNorm, _BraidAppAccess,
                 _BraidAppBufSize, _BraidAppBufPack, _BraidAppBufUnpack, &core);
   }

   void SetMaxLevels(braid_Int max_levels) { braid_SetMaxLevels(core, max_levels); }

   void SetMinCoarse(braid_Int min_coarse) { braid_SetMinCoarse(core, min_coarse); }

   void SetNRelax(braid_Int level, braid_Int nrelax)
   { braid_SetNRelax(core, level, nrelax); }

   void SetAbsTol(braid_Real tol) { braid_SetAbsTol(core, tol); }

   void SetRelTol(braid_Real tol) { braid_SetRelTol(core, tol); }

   void SetTemporalNorm(braid_Int tnorm) { braid_SetTemporalNorm(core, tnorm); }

   void SetCFactor(braid_Int level, braid_Int cfactor)
   { braid_SetCFactor(core, level, cfactor); }

   /** Use cfactor0 on all levels until there are < cfactor0 points
       on each processor. */
   void SetAggCFactor(braid_Int cfactor0)
   {
      BraidApp *app = (BraidApp *) core->app;
      braid_Int nt = app->ntime, pt;
      MPI_Comm_size(app->comm_t, &pt);
      if (cfactor0 > -1)
      {
         braid_Int level = (braid_Int) (log10((nt + 1) / pt) / log10(cfactor0));
         for (braid_Int i = 0; i < level; i++)
            braid_SetCFactor(core, i, cfactor0);
      }
   }

   void SetSpatialCoarsenAndRefine()
   {
      braid_SetSpatialCoarsen(core, _BraidAppCoarsen);
      braid_SetSpatialRefine(core, _BraidAppRefine);
   }

   void SetMaxIter(braid_Int max_iter) { braid_SetMaxIter(core, max_iter); }

   void SetPrintLevel(braid_Int print_level) { braid_SetPrintLevel(core, print_level); }

   void SetPrintFile(const char *printfile_name) { braid_SetPrintFile(core, printfile_name); }

   void SetAccessLevel(braid_Int access_level) { braid_SetAccessLevel(core, access_level); }

   void SetFMG() { braid_SetFMG(core); }

   void SetNFMGVcyc(braid_Int nfmg_Vcyc) { braid_SetNFMGVcyc(core, nfmg_Vcyc); }

   void GetNumIter(braid_Int *niter_ptr) { braid_GetNumIter(core, niter_ptr); }

   void GetRNorm(braid_Real *rnorm_ptr) { braid_GetRNorm(core, rnorm_ptr); }

   void Drive() { braid_Drive(core); }

   ~BraidCore() { braid_Destroy(core); }
};

// Wrapper for BRAID utilities that help the user,
// includes all the braid_Test* routines for testing the
// user-written wrappers.
class BraidUtil
{
private:

public:

   // Empty constructor
   BraidUtil() { }

   // Split comm_world into comm_x and comm_t, the spatial
   // and temporal communicators
   void SplitCommworld(const MPI_Comm  *comm_world,
                       braid_Int  px,
                       MPI_Comm  *comm_x,
                       MPI_Comm  *comm_t)
   { braid_SplitCommworld(comm_world, px, comm_x, comm_t); }

   // Test Function for Init and Access function
   void TestInitAccess(BraidApp   *app,
                       MPI_Comm    comm_x,
                       FILE       *fp,
                       braid_Real  t)
   { braid_TestInitAccess((braid_App) app, comm_x, fp, t,
         _BraidAppInit, _BraidAppAccess, _BraidAppFree); }

   // Test Function for Clone
   void TestClone(BraidApp   *app,
                  MPI_Comm    comm_x,
                  FILE       *fp,
                  braid_Real  t)
   { braid_TestClone((braid_App) app, comm_x, fp, t,
         _BraidAppInit, _BraidAppAccess, _BraidAppFree,
         _BraidAppClone); }

   // Test Function for Sum
   void TestSum(BraidApp   *app,
                MPI_Comm    comm_x,
                FILE       *fp,
                braid_Real  t)
   { braid_TestSum((braid_App) app, comm_x, fp, t,
         _BraidAppInit, _BraidAppAccess, _BraidAppFree,
         _BraidAppClone, _BraidAppSum); }

   // Test Function for SpatialNorm
   braid_Int TestSpatialNorm(BraidApp   *app,
                             MPI_Comm    comm_x,
                             FILE       *fp,
                             braid_Real  t)
   { return braid_TestSpatialNorm((braid_App) app, comm_x, fp, t,
         _BraidAppInit, _BraidAppFree, _BraidAppClone,
         _BraidAppSum, _BraidAppSpatialNorm); }

   // Test Functions BufSize, BufPack, BufUnpack
   braid_Int TestBuf(BraidApp   *app,
                     MPI_Comm    comm_x,
                     FILE       *fp,
                     braid_Real  t)
   { return braid_TestBuf((braid_App) app, comm_x, fp, t,
         _BraidAppInit, _BraidAppFree, _BraidAppSum,
         _BraidAppSpatialNorm, _BraidAppBufSize, _BraidAppBufPack,
         _BraidAppBufUnpack); }

   // Test Functions Coarsen and Refine
   braid_Int TestCoarsenRefine(BraidApp   *app,
                               MPI_Comm    comm_x,
                               FILE       *fp,
                               braid_Real  t,
                               braid_Real  fdt,
                               braid_Real  cdt)
   { return braid_TestCoarsenRefine((braid_App) app, comm_x, fp, t, fdt, cdt,
         _BraidAppInit, _BraidAppAccess, _BraidAppFree,
         _BraidAppClone, _BraidAppSum, _BraidAppSpatialNorm,
         _BraidAppCoarsen, _BraidAppRefine); }

   braid_Int TestAll(BraidApp   *app,
                     MPI_Comm    comm_x,
                     FILE       *fp,
                     braid_Real  t,
                     braid_Real  fdt,
                     braid_Real  cdt)
   { return braid_TestAll((braid_App) app, comm_x, fp, t, fdt, cdt,
         _BraidAppInit, _BraidAppFree, _BraidAppClone,
         _BraidAppSum, _BraidAppSpatialNorm, _BraidAppBufSize,
         _BraidAppBufPack, _BraidAppBufUnpack, _BraidAppCoarsen,
         _BraidAppRefine); }

   ~BraidUtil() { }
};

#endif
