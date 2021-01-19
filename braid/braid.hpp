// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// Produced at the Lawrence Livermore National Laboratory.
// 
// This file is part of XBraid. For support, post issues to the XBraid Github page.
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
class BraidSyncStatus;
class BraidStepStatus;
class BraidCoarsenRefStatus;
class BraidBufferStatus;
class BraidObjectiveStatus;

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

   /** @brief Apply the time stepping routine to the input vector @a u_
       corresponding to time @a tstart, and return in the same vector @a u_ the
       computed result for time @a tstop. The values of @a tstart and @a tstop
       can be obtained from @a pstatus.

       @param[in,out] u_ Input: approximate solution at time @a tstart.
                         Output: computed solution at time @a tstop.
       @param[in] ustop_ Previous approximate solution at @a tstop?
       @param[in] fstop_ Additional source at time @a tstop. May be set to NULL,
                         indicating no additional source.

       @see braid_PtFcnStep. */
   virtual braid_Int Step(braid_Vector     u_,
                          braid_Vector     ustop_,
                          braid_Vector     fstop_,
                          BraidStepStatus &pstatus) = 0;

   /** @brief Compute the residual at time @a tstop, given the approximate
       solutions at @a tstart and @a tstop. The values of @a tstart and @a tstop
       can be obtained from @a pstatus.

       @param[in]     u_ Input: approximate solution at time @a tstop.
       @param[in,out] r_ Input: approximate solution at time @a tstart.
                         Output: residual at time @a tstop.

       @see braid_PtFcnResidual.
   */
   virtual braid_Int Residual(braid_Vector     u_,
                              braid_Vector     r_,
                              BraidStepStatus &pstatus) = 0;

   /// Allocate a new vector in @a *v_ptr, which is a deep copy of @a u_.
   /// @see braid_PtFcnClone.
   virtual braid_Int Clone(braid_Vector  u_,
                           braid_Vector *v_ptr) = 0;

   /** @brief Allocate a new vector in @a *u_ptr and initialize it with an
       initial guess appropriate for time @a t. If @a t is the starting time,
       this method should set @a *u_ptr to the initial value vector of the ODE
       problem.
       @see braid_PtFcnInit. */
   virtual braid_Int Init(braid_Real    t,
                          braid_Vector *u_ptr) = 0;

   /// De-allocate the vector @a u_.
   /// @see braid_PtFcnFree.
   virtual braid_Int Free(braid_Vector u_) = 0;

   /// Perform the operation: @a y_ = @a alpha * @a x_ + @a beta * @a y_.
   /// @see braid_PtFcnSum.
   virtual braid_Int Sum(braid_Real   alpha,
                         braid_Vector x_,
                         braid_Real   beta,
                         braid_Vector y_) = 0;

   /// Compute in @a *norm_ptr an appropriate spatial norm of @a u_.
   /// @see braid_PtFcnSpatialNorm.
   virtual braid_Int SpatialNorm(braid_Vector  u_,
                                 braid_Real   *norm_ptr) = 0;

   /// @see braid_PtFcnAccess.
   virtual braid_Int Access(braid_Vector       u_,
                            BraidAccessStatus &astatus) = 0;

   /// @see braid_PtFcnBufSize.
   virtual braid_Int BufSize(braid_Int         *size_ptr,
                             BraidBufferStatus &bstatus) = 0;

   /// @see braid_PtFcnBufPack.
   virtual braid_Int BufPack(braid_Vector       u_,
                             void              *buffer,
                             BraidBufferStatus &bstatus) = 0;

   /// @see braid_PtFcnBufUnpack.
   virtual braid_Int BufUnpack(void              *buffer,
                               braid_Vector      *u_ptr,
                               BraidBufferStatus &bstatus) = 0;

   // These two functions may be optionally defined by the user, if spatial
   // coarsening is desired (see documentation for more details).  To turn on
   // spatial coarsening, use core.SetSpatialCoarsenAndRefine()
   /// @see braid_PtFcnSCoarsen.
   virtual braid_Int Coarsen(braid_Vector           fu_,
                             braid_Vector          *cu_ptr,
                             BraidCoarsenRefStatus &status)
   {
      fprintf(stderr, "Braid C++ Wrapper Warning: turn off spatial coarsening "
                      "until Coarsen and Refine have been user implemented\n");
      Clone(fu_, cu_ptr);
      return 0;
   }

   /// @see braid_PtFcnSRefine.
   virtual braid_Int Refine(braid_Vector           cu_,
                            braid_Vector          *fu_ptr,
                            BraidCoarsenRefStatus &status)
   {
      fprintf(stderr, "Braid C++ Wrapper Warning: turn off spatial coarsening "
                      "until Coarsen and Refine have been user implemented\n");
      Clone(cu_, fu_ptr);
      return 0;
   }

   /// This function may be optionally defined by the user. This provides
   /// a once-per-processor access to the user's app function at certain
   /// points during the code (see documentation for more details).
   /// To turn on sync, use core.SetSync()
   /// @see braid_PtFcnSync.
   virtual braid_Int Sync(BraidSyncStatus &sstatus)
   {
      fprintf(stderr, "Braid C++ Wrapper Warning: turn off sync "
              "until the Sync function been user implemented\n");
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

      void GetTILD(braid_Real *t_ptr,
                   braid_Int  *iter_ptr,
                   braid_Int  *level_ptr,
                   braid_Int  *done_ptr)
      { braid_AccessStatusGetTILD(astatus, t_ptr, iter_ptr, level_ptr, done_ptr); }
      void GetT(braid_Real *t_ptr)                       { braid_AccessStatusGetT(astatus, t_ptr); }
      void GetTIndex(braid_Int *tindex_ptr)              { braid_AccessStatusGetTIndex(astatus, tindex_ptr); }
      void GetDone(braid_Int *done_ptr)                  { braid_AccessStatusGetDone(astatus, done_ptr); }
      void GetLevel(braid_Int *level_ptr)                { braid_AccessStatusGetLevel(astatus, level_ptr); }
      void GetNLevels(braid_Int *nlevels_ptr)            { braid_AccessStatusGetNLevels(astatus, nlevels_ptr);}
      void GetIter(braid_Int *iter_ptr)                  { braid_AccessStatusGetIter(astatus, iter_ptr); }
      void GetWrapperTest(braid_Int *wtest_ptr)          { braid_AccessStatusGetWrapperTest(astatus, wtest_ptr); }
      void GetResidual(braid_Real *rnorm_ptr)            { braid_AccessStatusGetResidual(astatus, rnorm_ptr); }
      void GetNRefine(braid_Int *nrefine_ptr)            { braid_AccessStatusGetNRefine(astatus, nrefine_ptr); }
      void GetNTPoints(braid_Int *ntpoints_ptr)          { braid_AccessStatusGetNTPoints(astatus, ntpoints_ptr); }
      void GetSingleErrorEstAccess(braid_Real *estimate_ptr)   { braid_AccessStatusGetSingleErrorEstAccess(astatus, estimate_ptr); }
      void GetCallingFunction(braid_Int *callingfcn_ptr)
      {
         braid_AccessStatusGetCallingFunction(astatus, callingfcn_ptr);
      }

      // The braid_AccessStatus structure is deallocated inside of Braid
      // This class is just to make code consistently look object oriented
      ~BraidAccessStatus() { }
};

// Wrapper for BRAID's SyncStatus object
class BraidSyncStatus
{
   private:
      braid_SyncStatus sstatus;

   public:
      BraidSyncStatus(braid_SyncStatus _sstatus)
      {
         sstatus = _sstatus;
      }
      void GetTIUL(braid_Int *i_upper,
                   braid_Int *i_lower,
                   braid_Int  level)
      { braid_SyncStatusGetTIUL(sstatus, i_upper, i_lower, level); }
      void GetTimeValues(braid_Real **tvalues_ptr,
                         braid_Int    i_upper,
                         braid_Int    i_lower,
                         braid_Int    level)
      { braid_SyncStatusGetTimeValues(sstatus, tvalues_ptr, i_upper, i_lower, level); }
      void GetNLevels(braid_Int *nlevels_ptr)        { braid_SyncStatusGetNLevels(sstatus, nlevels_ptr);}
      void GetIter(braid_Int *iter_ptr)              { braid_SyncStatusGetIter(sstatus, iter_ptr); }
      void GetLevel(braid_Int *level_ptr)            { braid_SyncStatusGetLevel(sstatus, level_ptr); }
      void GetNRefine(braid_Int *nrefine_ptr)        { braid_SyncStatusGetNRefine(sstatus, nrefine_ptr); }
      void GetNTPoints(braid_Int *ntpoints_ptr)      { braid_SyncStatusGetNTPoints(sstatus, ntpoints_ptr); }
      void GetDone(braid_Int *done_ptr)              { braid_SyncStatusGetDone(sstatus, done_ptr); }
      void GetNumErrorEst(braid_Int *npoints_ptr)    { braid_SyncStatusGetNumErrorEst(sstatus, npoints_ptr); }
      void GetAllErrorEst(braid_Real *error_est_ptr) { braid_SyncStatusGetAllErrorEst(sstatus, error_est_ptr); }
      void GetCallingFunction(braid_Int *callingfcn_ptr)
      {
         braid_SyncStatusGetCallingFunction(sstatus, callingfcn_ptr);
      }

      // The braid_SyncStatus structure is deallocated inside of Braid
      // This class is just to make code consistently look object oriented
      ~BraidSyncStatus() { }
};

// Wrapper for BRAID's StepStatus object
class BraidStepStatus
{
   private:
      braid_StepStatus pstatus;

   public:
      BraidStepStatus(braid_StepStatus _pstatus)
      {
         pstatus = _pstatus;
      }

      void GetRNorms(braid_Int *nrequest_ptr, braid_Real *rnorms)
      { braid_StepStatusGetRNorms(pstatus, nrequest_ptr, rnorms); }
      void GetTstartTstop(braid_Real *tstart_ptr, braid_Real *tstop_ptr)
      { braid_StepStatusGetTstartTstop(pstatus, tstart_ptr, tstop_ptr); }
      void GetT(braid_Real *tstart_ptr)                  { braid_StepStatusGetT(pstatus, tstart_ptr); }
      void GetTstop(braid_Real *tstop_ptr)               { braid_StepStatusGetTstop(pstatus, tstop_ptr); }
      void GetTIndex(braid_Int *tindex_ptr)              { braid_StepStatusGetTIndex(pstatus, tindex_ptr); }
      void GetLevel(braid_Int *level_ptr)                { braid_StepStatusGetLevel(pstatus, level_ptr); }
      void GetNLevels(braid_Int *nlevels_ptr)            { braid_StepStatusGetNLevels(pstatus, nlevels_ptr); }
      void GetNRefine(braid_Int *nrefine_ptr)            { braid_StepStatusGetNRefine(pstatus, nrefine_ptr); }
      void GetNTPoints(braid_Int *ntpoints_ptr)          { braid_StepStatusGetNTPoints(pstatus, ntpoints_ptr); }
      void SetRFactor(braid_Int rfactor)                 { braid_StepStatusSetRFactor(pstatus, rfactor); }
      void SetRSpace(braid_Int rspace)                   { braid_StepStatusSetRSpace(pstatus, rspace); }
      void GetTol(braid_Real *tol_ptr)                   { braid_StepStatusGetTol(pstatus, tol_ptr); }
      void GetIter(braid_Int *iter_ptr)                  { braid_StepStatusGetIter(pstatus, iter_ptr); }
      void GetOldFineTolx(braid_Real *old_fine_tolx_ptr) { braid_StepStatusGetOldFineTolx(pstatus, old_fine_tolx_ptr); }
      void SetOldFineTolx(braid_Real old_fine_tolx)      { braid_StepStatusSetOldFineTolx(pstatus, old_fine_tolx); }
      void SetTightFineTolx(braid_Int tight_fine_tolx)   { braid_StepStatusSetTightFineTolx(pstatus, tight_fine_tolx); }
      void GetSingleErrorEstStep(braid_Real *estimate_ptr)   { braid_StepStatusGetSingleErrorEstStep(pstatus, estimate_ptr); }

         // The braid_StepStatus structure is deallocated inside of Braid
      // This class is just to make code consistently look object oriented
      ~BraidStepStatus() { }

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

      void GetTpriorTstop(braid_Real *tstart_ptr,
                          braid_Real *f_tprior_ptr,
                          braid_Real *f_tstop_ptr,
                          braid_Real *c_tprior_ptr,
                          braid_Real *c_tstop_ptr)
      {
         braid_CoarsenRefStatusGetTpriorTstop(
            cstatus, tstart_ptr, f_tprior_ptr, f_tstop_ptr,
            c_tprior_ptr, c_tstop_ptr);
      }
      void GetT(braid_Real *tstart_ptr)         { braid_CoarsenRefStatusGetT(cstatus, tstart_ptr); }
      void GetTIndex(braid_Int *tindex_ptr)     { braid_CoarsenRefStatusGetTIndex(cstatus, tindex_ptr); }
      void GetIter(braid_Int *iter_ptr)         { braid_CoarsenRefStatusGetIter(cstatus, iter_ptr); }
      void GetFTstop(braid_Real *f_tstop_ptr)   { braid_CoarsenRefStatusGetFTstop(cstatus, f_tstop_ptr); }
      void GetFTprior(braid_Real *f_tprior_ptr) { braid_CoarsenRefStatusGetFTprior(cstatus, f_tprior_ptr); }
      void GetCTstop(braid_Real *c_tstop_ptr)   { braid_CoarsenRefStatusGetCTstop(cstatus, c_tstop_ptr); }
      void GetCTprior(braid_Real *c_tprior_ptr) { braid_CoarsenRefStatusGetCTprior(cstatus, c_tprior_ptr); }
      void GetLevel(braid_Int *level_ptr)       { braid_CoarsenRefStatusGetLevel(cstatus, level_ptr); }
      void GetNLevels(braid_Int *nlevels_ptr)   { braid_CoarsenRefStatusGetNLevels(cstatus, nlevels_ptr); }
      void GetNRefine(braid_Int *nrefine_ptr)   { braid_CoarsenRefStatusGetNRefine(cstatus, nrefine_ptr); }
      void GetNTPoints(braid_Int *ntpoints_ptr) { braid_CoarsenRefStatusGetNTPoints(cstatus, ntpoints_ptr); }

      // The braid_CoarsenRefStatus structure is deallocated inside of Braid
      // This class is just to make code consistently look object oriented
      ~BraidCoarsenRefStatus() { }
};

// Wrapper for BRAID's BufferStatus object
class BraidBufferStatus
{
   private:
      braid_BufferStatus bstatus;
   
   public:
      BraidBufferStatus( braid_BufferStatus _bstatus )
      {
         bstatus = _bstatus;
      }

      void GetMessageType( braid_Int *messagetype_ptr ) { braid_BufferStatusGetMessageType( bstatus, messagetype_ptr); }
      void SetSize( braid_Int size ) { braid_BufferStatusSetSize( bstatus, size ); }
      ~BraidBufferStatus() {} 
};

// Wrapper for BRAID's ObjectiveStatus object
class BraidObjectiveStatus
{
   private:
      braid_ObjectiveStatus ostatus;

   public:
      BraidObjectiveStatus(braid_ObjectiveStatus  _ostatus)
      {
         ostatus = _ostatus;
      }

      void GetT(braid_Real *tstart_ptr)         { braid_ObjectiveStatusGetT(ostatus, tstart_ptr); }
      void GetTIndex(braid_Int *tindex_ptr)     { braid_ObjectiveStatusGetTIndex(ostatus, tindex_ptr); }
      void GetIter(braid_Int *iter_ptr)         { braid_ObjectiveStatusGetIter(ostatus, iter_ptr); }
      void GetLevel(braid_Int *level_ptr)       { braid_ObjectiveStatusGetLevel(ostatus, level_ptr); }
      void GetNLevels(braid_Int *nlevels_ptr)   { braid_ObjectiveStatusGetNLevels(ostatus, nlevels_ptr); }
      void GetNRefine(braid_Int *nrefine_ptr)   { braid_ObjectiveStatusGetNRefine(ostatus, nrefine_ptr); }
      void GetNTPoints(braid_Int *ntpoints_ptr) { braid_ObjectiveStatusGetNTPoints(ostatus, ntpoints_ptr); }
      void GetTol(braid_Real *tol_ptr)          { braid_ObjectiveStatusGetTol(ostatus, tol_ptr); }

      // The braid_ObjectiveStatus structure is deallocated inside of Braid
      // This class is just to make code consistently look object oriented
      ~BraidObjectiveStatus() { }
};

// Static functions passed to Braid, with braid_App == BraidApp*
static braid_Int _BraidAppStep(braid_App       _app,
                              braid_Vector    _ustop,
                              braid_Vector    _fstop,
                              braid_Vector    _u,
                              braid_StepStatus _pstatus)
{
   BraidApp *app = (BraidApp*)_app;
   BraidStepStatus pstatus(_pstatus);
   return app -> Step(_u,_ustop, _fstop, pstatus);
}


static braid_Int _BraidAppResidual(braid_App     _app,
                              braid_Vector  _ustop,
                              braid_Vector  _r,
                              braid_StepStatus _pstatus)
{
   BraidApp *app = (BraidApp*)_app;
   BraidStepStatus pstatus(_pstatus);
   return app -> Residual(_ustop,_r, pstatus);
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
   BraidAccessStatus astatus(_astatus);
   return app -> Access(_u, astatus);
}

static braid_Int _BraidAppSync(braid_App        _app,
                               braid_SyncStatus _sstatus)
{
   BraidApp *app = (BraidApp*)_app;
   BraidSyncStatus sstatus(_sstatus);
   return app -> Sync(sstatus);
}

static braid_Int _BraidAppBufSize(braid_App  _app,
                                  braid_Int *size_ptr,
                                  braid_BufferStatus _bstatus)
{
   BraidApp *app = (BraidApp*)_app;
   BraidBufferStatus bstatus( _bstatus );
   return app -> BufSize(size_ptr, bstatus);
}


static braid_Int _BraidAppBufPack(braid_App     _app,
                                  braid_Vector  _u,
                                  void         *buffer,
                                  braid_BufferStatus  _bstatus)
{
   BraidApp *app = (BraidApp*)_app;
   BraidBufferStatus bstatus( _bstatus );
   return app -> BufPack(_u, buffer, bstatus);
}


static braid_Int _BraidAppBufUnpack(braid_App     _app,
                                    void         *buffer,
                                    braid_Vector *u_ptr,
                                    braid_BufferStatus _bstatus)
{
   BraidApp *app = (BraidApp*)_app;
   BraidBufferStatus bstatus( _bstatus );
   return app -> BufUnpack(buffer, u_ptr, bstatus);
}


static braid_Int _BraidAppCoarsen(braid_App               _app,
                                  braid_Vector            _fu,
                                  braid_Vector           *cu_ptr,
                                  braid_CoarsenRefStatus  _cstatus)
{
   BraidApp *app = (BraidApp*)_app;
   BraidCoarsenRefStatus cstatus(_cstatus);
   return app -> Coarsen(_fu, cu_ptr, cstatus);
}


static braid_Int _BraidAppRefine(braid_App               _app,
                                 braid_Vector            _fu,
                                 braid_Vector           *cu_ptr,
                                 braid_CoarsenRefStatus  _cstatus)
{
   BraidApp *app = (BraidApp*)_app;
   BraidCoarsenRefStatus cstatus(_cstatus);
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
                 _BraidAppStep, _BraidAppInit, _BraidAppClone, _BraidAppFree,
                 _BraidAppSum, _BraidAppSpatialNorm, _BraidAppAccess,
                 _BraidAppBufSize, _BraidAppBufPack, _BraidAppBufUnpack, &core);
   }

   void SetMaxLevels(braid_Int max_levels) { braid_SetMaxLevels(core, max_levels); }

   void SetIncrMaxLevels() { braid_SetIncrMaxLevels(core); }

   void SetSkip(braid_Int skip) { braid_SetSkip(core, skip); }

   void SetMinCoarse(braid_Int min_coarse) { braid_SetMinCoarse(core, min_coarse); }

   void SetNRelax(braid_Int level, braid_Int nrelax)
   { braid_SetNRelax(core, level, nrelax); }

   void SetAbsTol(braid_Real tol) { braid_SetAbsTol(core, tol); }

   void SetRelTol(braid_Real tol) { braid_SetRelTol(core, tol); }

   void SetTemporalNorm(braid_Int tnorm) { braid_SetTemporalNorm(core, tnorm); }

   void SetCFactor(braid_Int level, braid_Int cfactor)
   { braid_SetCFactor(core, level, cfactor); }

   /// If cfactor0 > -1, set the cfactor for level 0 to cfactor0.
   void SetAggCFactor(braid_Int cfactor0)
   {
      if (cfactor0 > -1)
         braid_SetCFactor(core, 0, cfactor0);

    /* Use cfactor0 on all levels until there are < cfactor0 points
       on each processor. */
    //BraidApp *app = (BraidApp *) core->app;
    //braid_Int nt = app->ntime, pt;
    //MPI_Comm_size(app->comm_t, &pt);
    //if (cfactor0 > -1)
    //{
    //   braid_Int level = (braid_Int) (log10((nt + 1) / pt) / log10(cfactor0));
    //   for (braid_Int i = 0; i < level; i++)
    //      braid_SetCFactor(core, i, cfactor0);
    //}
   }

   void SetSpatialCoarsenAndRefine()
   {
      braid_SetSpatialCoarsen(core, _BraidAppCoarsen);
      braid_SetSpatialRefine(core, _BraidAppRefine);
   }

   void SetPeriodic(braid_Int periodic) { braid_SetPeriodic(core, periodic); }

   void SetSync() { braid_SetSync(core, _BraidAppSync); }

   void SetResidual() { braid_SetResidual(core, _BraidAppResidual); }

   void SetMaxIter(braid_Int max_iter) { braid_SetMaxIter(core, max_iter); }

   void SetPrintLevel(braid_Int print_level) { braid_SetPrintLevel(core, print_level); }

   void SetSeqSoln(braid_Int use_seq_soln) { braid_SetSeqSoln(core, use_seq_soln); }

   void SetPrintFile(const char *printfile_name) { braid_SetPrintFile(core, printfile_name); }

   void SetAccessLevel(braid_Int access_level) { braid_SetAccessLevel(core, access_level); }

   void SetFMG() { braid_SetFMG(core); }

   void SetNFMG(braid_Int k) { braid_SetNFMG(core, k); }

   void SetNFMGVcyc(braid_Int nfmg_Vcyc) { braid_SetNFMGVcyc(core, nfmg_Vcyc); }

   void SetStorage(braid_Int storage) { braid_SetStorage(core, storage); }

   void SetRefine(braid_Int refine) {braid_SetRefine(core, refine);}

   void SetMaxRefinements(braid_Int max_refinements) {braid_SetMaxRefinements(core, max_refinements);}

   void SetRichardsonEstimation(braid_Int est_error, braid_Int richardson, braid_Int local_order) {braid_SetRichardsonEstimation(core, est_error, richardson, local_order);}

   void GetNumIter(braid_Int *niter_ptr) { braid_GetNumIter(core, niter_ptr); }

   void GetRNorms(braid_Int *nrequest_ptr, braid_Real *rnorms) { braid_GetRNorms(core, nrequest_ptr, rnorms); }
   
   void GetNLevels(braid_Int *nlevels_ptr) { braid_GetNLevels(core, nlevels_ptr); }

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

   // Return an appropriate spatial tolerance
   void GetSpatialAccuracy(braid_StepStatus  sstatus,
                           braid_Real        loose_tol,
                           braid_Real        tight_tol,
                           braid_Real        *tol_ptr)
   { braid_GetSpatialAccuracy(sstatus, loose_tol, tight_tol, tol_ptr); }

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
   
   // Test Functions Coarsen and Refine
   braid_Int TestResidual(BraidApp   *app,
                          MPI_Comm    comm_x,
                          FILE       *fp,
                          braid_Real  t,
                          braid_Real  dt)
   { return braid_TestResidual((braid_App) app, comm_x, fp, t, dt,
         _BraidAppInit, _BraidAppAccess, _BraidAppFree,
         _BraidAppClone, _BraidAppSum, _BraidAppSpatialNorm,
         _BraidAppResidual, _BraidAppStep); }

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
         _BraidAppRefine, _BraidAppResidual, _BraidAppStep); }

   ~BraidUtil() { }
};

#endif
