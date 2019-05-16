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

/** \file _braid.c
 * \brief Source code for developer routines.  See braid.h for more information.
 *
 */

#include "_braid.h"
#include "_util.h"

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int _braid_error_flag = 0;
FILE    *_braid_printfile  = NULL;

braid_Int 
_braid_VectorBarCopy(braid_VectorBar bar, braid_VectorBar *bar_ptr)
{
  bar->useCount++;
  *bar_ptr= bar;

  return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_VectorBarDelete(braid_Core core, braid_VectorBar bar)
{
   /* Decrease the useCount */
   bar->useCount--;

   /* Free memory, if no pointer is left */
   if (bar->useCount==0)
   {
       _braid_CoreFcn(core, free)(_braid_CoreElt(core, app), bar->userVector);
       free(bar);
   }
 
   /* Sanity check */
   else if (bar->useCount < 0)
   {
     printf("ERROR: useCount < 0 !\n");
     exit(0);
   }
 
   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_OptimDestroy( braid_Core core)
{
   braid_App    app        = _braid_CoreElt(core, app);
   braid_Optim  optim      = _braid_CoreElt(core, optim);

   /* Only deallocate if braid_Drive() has been called, i.e., a grid structure 
    * has been created */
   if (_braid_CoreElt(core, grids)[0] != NULL)
   {
     _braid_Grid *fine_grid  = _braid_CoreElt(core, grids)[0];
     braid_Int    storage    = _braid_CoreElt(core, storage);
     braid_Int    clower     = _braid_GridElt(fine_grid, clower);
     braid_Int    iupper     = _braid_GridElt(fine_grid, iupper);
     braid_Int    ilower     = _braid_GridElt(fine_grid, ilower);
     braid_Int    cfactor    = _braid_GridElt(fine_grid, cfactor);
     braid_Int    ic, iclocal, sflag, increment, destroy_flag;

     /* Get the number of adjoint vectors on finest level */
     if (storage < 0 ) 
     {
        /* Only C-point storage */
        ilower    = clower;
        increment = cfactor;
     }
     else
     {
        /* All points */
        increment = 1;
     }

     /* Free the adjoint variables and tapeinput */
     for (ic=ilower; ic <= iupper; ic += increment)
     {
        destroy_flag = 1;

        /* if only C-point storage, destroy only at C-points */
        if (storage < 0 &&  !(_braid_IsCPoint(ic, cfactor)) )
        {
           destroy_flag = 0;
        } 

        if (destroy_flag)
        {
           _braid_UGetIndex(core, 0, ic, &iclocal, &sflag);
           _braid_CoreFcn(core, free)( app, optim->adjoints[iclocal]);
           _braid_VectorBarDelete(core, optim->tapeinput[iclocal] );
        }
     }

   }

   free(optim->sendbuffer);
   free(optim->adjoints);
   free(optim->tapeinput);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_UpdateAdjoint(braid_Core core,
                     braid_Real *rnorm_adj_ptr)
{
   braid_App    app       = _braid_CoreElt(core, app);
   MPI_Comm     comm      = _braid_CoreElt(core, comm);
   braid_Optim  optim     = _braid_CoreElt(core, optim);
   braid_Int    tnorm     = _braid_CoreElt(core, tnorm);
   braid_Int    storage   = _braid_CoreElt(core, storage);
   braid_Int    iter      = _braid_CoreElt(core, niter);
   _braid_Grid *fine_grid = _braid_CoreElt(core, grids)[0];
   braid_Int    clower    = _braid_GridElt(fine_grid, clower);
   braid_Int    iupper    = _braid_GridElt(fine_grid, iupper);
   braid_Int    ilower    = _braid_GridElt(fine_grid, ilower);
   braid_Int    cfactor   = _braid_GridElt(fine_grid, cfactor);
   braid_Real   rnorm_adj, rnorm_temp, global_rnorm;
   braid_Vector tape_vec, adjoint_vec;
   braid_Int    ic, iclocal, sflag, increment, upd_flag;
 
   rnorm_adj    = 0.;
   global_rnorm = 0.;

   /* Get the number of adjoint vectors on finest level */
   if (storage < 0 ) 
   {
      /* Only C-point storage */
      ilower    = clower;
      increment = cfactor;
   }
   else
   {
      /* All points */
      increment = 1;
   }

   /* Loop over all adjoint vectors on the fine grid */
   for (ic=ilower; ic <= iupper; ic += increment)
   {
      upd_flag = 1;
      
      /* If first iteration, update only C-points */
      if (iter == 0)
      {
         if (!_braid_IsCPoint(ic, cfactor))
         {
            upd_flag = 0;
         }
      }

      /* If only C-point storage, update only C-points, else update all */
      if (storage < 0 && !_braid_IsCPoint(ic, cfactor) )
      {
         upd_flag = 0;
      }

      /* Compute norm and update */
      if(upd_flag)
      {
         /* Get the local index of the points */
         _braid_UGetIndex(core, 0, ic, &iclocal, &sflag);

         tape_vec    = optim->tapeinput[iclocal]->userVector;   
         adjoint_vec = optim->adjoints[iclocal];

         if (ic > 0)
         {
            /* Compute the norm of the adjoint residual */
            _braid_CoreFcn(core, sum)(app, 1., tape_vec, -1., adjoint_vec);
            _braid_CoreFcn(core, spatialnorm)(app, adjoint_vec, &rnorm_temp);
            if(tnorm == 1)       /* one-norm */ 
            {  
               rnorm_adj += rnorm_temp;
            }
            else if(tnorm == 3)  /* inf-norm */
            {  
               rnorm_adj = (((rnorm_temp) > (rnorm_adj)) ? (rnorm_temp) : (rnorm_adj));
            }
            else                 /* default two-norm */
            {  
               rnorm_adj += (rnorm_temp*rnorm_temp);
            }

            /* Update the adjoint variables */
            _braid_CoreFcn(core, sum)(app, 1., tape_vec , 0., adjoint_vec);
         }

         /* Delete the pointer */
         _braid_VectorBarDelete(core, optim->tapeinput[iclocal]);
      }
   }

   /* Compute global residual norm. */
   if(tnorm == 1)       /* one-norm reduction */
   {  
      MPI_Allreduce(&rnorm_adj, &global_rnorm, 1, braid_MPI_REAL, MPI_SUM, comm);
   }
   else if(tnorm == 3)  /* inf-norm reduction */
   {  
      MPI_Allreduce(&rnorm_adj, &global_rnorm, 1, braid_MPI_REAL, MPI_MAX, comm);
   }
   else                 /* default two-norm reduction */
   {  
      MPI_Allreduce(&rnorm_adj, &global_rnorm, 1, braid_MPI_REAL, MPI_SUM, comm);
      global_rnorm = sqrt(global_rnorm);
   }

   *rnorm_adj_ptr = global_rnorm;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int 
_braid_SetRNormAdjoint(braid_Core core, 
                       braid_Int iter, 
                       braid_Real rnorm_adj)
{
   braid_Optim optim = _braid_CoreElt(core, optim);

   /* Set the adjoint norm */
   optim->rnorm_adj = rnorm_adj;

   /* Set initial norm if not already set */
   if (iter == 0)
   {
      if (optim->rnorm0_adj == braid_INVALID_RNORM )
      {
         if (rnorm_adj > 0)
         {
            optim->rnorm0_adj = rnorm_adj;
         }
         else
         {
            optim->rnorm0_adj = 1.0;
         }
      }
   }

   return _braid_error_flag;
}                           

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_AddToObjective(braid_Core              core, 
                      braid_BaseVector      u, 
                      braid_ObjectiveStatus ostatus)
{
   braid_App   app        = _braid_CoreElt(core, app);
   braid_Optim optim      = _braid_CoreElt(core, optim);
   braid_Real  t          = _braid_CoreElt(core, t);
   braid_Real  tstart_obj = optim->tstart_obj;
   braid_Real  tstop_obj  = optim->tstop_obj;
   braid_Real  objT;

   if ( tstart_obj <= t && t <= tstop_obj )
   {
      /* Evaluate objective at time t */
     _braid_BaseObjectiveT(core, app, u, ostatus, &objT);

     /* Add to the time-averaged objective function */
     optim->sum_user_obj += objT;
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_EvalObjective(braid_Core core)
{
   MPI_Comm    comm         = _braid_CoreElt(core, comm);
   braid_App   app          = _braid_CoreElt(core, app);
   braid_Optim optim        = _braid_CoreElt(core, optim);
   braid_Real  sum_user_obj = optim->sum_user_obj;
   braid_Real localsum, globalsum, posttmp;
   braid_Real objective;
   
   /* Compute the global time average */
   localsum = sum_user_obj;
   MPI_Allreduce(&localsum, &globalsum, 1, braid_MPI_REAL, MPI_SUM, comm);

   /* Compute the global objective function, if set */
   if (_braid_CoreElt(core, postprocess_obj) != NULL)
   {
      _braid_CoreFcn(core, postprocess_obj) (app, globalsum, &posttmp);
      objective = posttmp;
   }
   else
   {
      objective = globalsum;
   }

   /* Store the timeaverage and objective function in the optim structure */
   optim->sum_user_obj = globalsum;
   optim->objective    = objective;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Differentiated objective function 
 *----------------------------------------------------------------------------*/

braid_Int
_braid_EvalObjective_diff(braid_Core core)
{
   braid_Optim optim        = _braid_CoreElt(core, optim);
   braid_App   app          = _braid_CoreElt(core, app);
   braid_Int   myid         = _braid_CoreElt(core, myid);
   braid_Real  sum_user_obj = optim->sum_user_obj;
   braid_Real  sum_user_obj_bar;

   /* Don't evaluate derivative if objective_only mode.
    * This function may not be defined by the user.*/
   if (_braid_CoreElt(core, obj_only))
   {
     return 0;
   }

   /* Differentiate the postprocessing objective, if set */
   if (_braid_CoreElt(core, postprocess_obj_diff != NULL))
   {
      _braid_CoreFcn(core, postprocess_obj_diff)(app, sum_user_obj, &sum_user_obj_bar);
   }
   else
   {
      sum_user_obj_bar = 1.0;
   }

   /* Differentiate the time average */
   optim->f_bar = sum_user_obj_bar;

   /* Reset the gradient on all but one processor (one should hold the differentiated relaxation term) */
   if ( myid != 0)
   {
      _braid_CoreFcn(core, reset_gradient)(app);
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_InitAdjointVars(braid_Core   core, 
                       _braid_Grid *fine_grid)
{
   braid_Real          tstart    = _braid_CoreElt(core, tstart);
   braid_App           app       = _braid_CoreElt(core, app);
   braid_Int           storage   = _braid_CoreElt(core, storage);
   braid_Int           ncpoints  = _braid_GridElt(fine_grid, ncpoints);
   braid_Int           clower    = _braid_GridElt(fine_grid, clower);
   braid_Int           iupper    = _braid_GridElt(fine_grid, iupper);
   braid_Int           ilower    = _braid_GridElt(fine_grid, ilower);
   braid_Int           cfactor   = _braid_GridElt(fine_grid, cfactor);
   braid_BufferStatus  bstatus   = (braid_BufferStatus) core;
   braid_Vector       *adjoints  = NULL; 
   braid_VectorBar    *tapeinput = NULL; 
   braid_BaseVector    u; 
   braid_VectorBar     bar_copy;
   braid_Vector        mybar;
   braid_Int           ic, iclocal, sflag, nupoints, increment;
   braid_Int           bufsize;
   void*               sendbuffer;
   MPI_Request*        request;
  

   /* Get the number of adjoint vectors on finest level */
   if (storage < 0 ) 
   {
      /* Only C-point storage */
      nupoints  = ncpoints;
      ilower    = clower;
      increment = cfactor;
   }
   else
   {
      /* All points */
      nupoints  = iupper-ilower+1; 
      increment = 1;
   }

   /* Allocate the adjoint variables and pointer to input variables */
   adjoints            = _braid_CTAlloc(braid_Vector, nupoints);
   tapeinput           = _braid_CTAlloc(braid_VectorBar, nupoints);

   /* Initialize adjoints and tapeinput*/
   braid_Int store = 0;
   for (ic=ilower; ic <= iupper; ic += increment)
   {
      if (storage < 0)
      {
         /* Init only C-points  */
         if( _braid_IsCPoint(ic, cfactor))
         {
            store = 1;
         }
      }
      else
      {
         /* Init all */
         store = 1;
      }

      if (store)
      {
         _braid_UGetIndex(core, 0, ic, &iclocal, &sflag);

         /* Initialize adjoint variables with zeros */
         _braid_CoreFcn(core, init)(app, tstart, &mybar);
         _braid_CoreFcn(core, sum)( app, -1.0, mybar, 1.0, mybar);
         adjoints[iclocal] = mybar;

         /* initialize the tapeinput with u_bar only at C-points */
         if( _braid_IsCPoint(ic, cfactor))
         {
            _braid_UGetVectorRef(core, 0, ic, &u);
            _braid_VectorBarCopy(u->bar, &bar_copy);
            tapeinput[iclocal] = bar_copy;
         }
      }
   }

   /* Allocate a buffer for BufUnpackDiff*/
   _braid_CoreFcn(core, bufsize)(app, &bufsize, bstatus);
   sendbuffer = malloc(bufsize); 
   request = NULL;

   /* Pass to the optimization structure */
   _braid_CoreElt(core, optim)->sendbuffer = sendbuffer;
   _braid_CoreElt(core, optim)->request    = request;
   _braid_CoreElt(core, optim)->adjoints  = adjoints;
   _braid_CoreElt(core, optim)->tapeinput = tapeinput;

   return _braid_error_flag;
}                

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_SetVerbosity(braid_Core  core,
                    braid_Int   verbose_adj)
{
   _braid_CoreElt(core, verbose_adj) = verbose_adj;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_AdjointFeatureCheck(braid_Core core)
{

   braid_PtFcnResidual  residual  = _braid_CoreElt(core, residual);
   braid_PtFcnResidual  fullres   = _braid_CoreElt(core, full_rnorm_res);
   braid_PtFcnSCoarsen  scoarsen  = _braid_CoreElt(core, scoarsen);
   braid_PtFcnSRefine   srefine   = _braid_CoreElt(core, srefine);
   braid_PtFcnTimeGrid  tgrid     = _braid_CoreElt(core, tgrid);
   braid_Int            storage   = _braid_CoreElt(core, storage);  
   braid_Int            useshell  = _braid_CoreElt(core, useshell);
   braid_Int            trefine   = _braid_CoreElt(core, refine);
   braid_Int err;
   char* err_char;

   err = 0;
   if ( residual != NULL ) 
   {
      err_char = "User-defined residual" ;
      err = 1;
   }
   if ( fullres  != NULL ) 
   {
      err_char = "User-defined full residual";
      err = 1;
   }
   if ( scoarsen != NULL ) 
   {
      err_char = "Spatial coarsening";
      err = 1;
   }
   if ( srefine  != NULL ) 
   {
      err_char = "Spatial refinement";
      err = 1;
   }
   if ( tgrid    != NULL ) 
   {
      err_char = "User-defined time grid";
      err = 1;
   }
   if ( useshell ) 
   {
      err_char = "Shell-vector feature";
      err = 1;
   }
   if ( trefine )  
   {
      err_char = "Time refinement";
      err = 1;
   }
   if (storage  >= 1 ) 
   {
      err_char = "Storage >= 1";
      err = 1;
   }

    // r_space?
   if ( err )
   {
      _braid_printf(" \n\n WARNING! %s is not yet supported for adjoint sensitivities (or at least not tested).\n", err_char); 
      _braid_printf("          If the code still runs, check the derivatives carefully!\n\n\n", err_char); 
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_ChunkFeatureCheck(braid_Core core)
{
   braid_Int            nchunks   = _braid_CoreElt(core, nchunks);
   braid_Int            ntime0    = _braid_CoreElt(core, ntime);
   braid_Int            useshell  = _braid_CoreElt(core, useshell);
   braid_PtFcnTimeGrid  tgrid     = _braid_CoreElt(core, tgrid);



   /* Sanity check */
   if (ntime0 % nchunks != 0)
   {
      _braid_printf("\n Error: ntimes must be a multiple of nchunks!\n");
      exit(1);
   }
   if ( tgrid != NULL && nchunks > 1 ) 
   {
      _braid_printf("\nUser-defined time grid not supported with more than one time chunk!\n");
      exit(1);
   }
   if ( useshell && nchunks > 1) 
   {
      _braid_printf("\nShell-vector feature not supported with more than one time chunk!\n");
      exit(1);
   }
   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * ZTODO: Should we use the error handling facility here and above?
 *----------------------------------------------------------------------------*/

braid_Int
_braid_TriMGRITFeatureCheck(braid_Core core)
{
#if 0
   braid_Int  nchunks  = _braid_CoreElt(core, nchunks );
   braid_Int  seq_soln = _braid_CoreElt(core, seq_soln);
   braid_Int  fmg      = _braid_CoreElt(core, fmg     );
   braid_Int  storage  = _braid_CoreElt(core, storage );
   braid_Int  useshell = _braid_CoreElt(core, useshell);
   braid_Int  refine   = _braid_CoreElt(core, refine  );
   braid_Int  r_space  = _braid_CoreElt(core, r_space );
   braid_Int  skip     = _braid_CoreElt(core, skip    );

   braid_Int  err = 0;
   char      *err_char;

   if ( err )
   {
      _braid_printf(" \n\n WARNING! %s is not yet supported\n", err_char); 
   }
#endif

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FeatureCheck(braid_Core core)
{
   braid_Int adjoint  = _braid_CoreElt(core, adjoint);
   braid_Int nchunks  = _braid_CoreElt(core, nchunks);
   braid_Int trimgrit = _braid_CoreElt(core, trimgrit);

   if (adjoint)
   {
      _braid_AdjointFeatureCheck(core);
   }
   if (nchunks > 1)
   {
      _braid_ChunkFeatureCheck(core);
   }
   if (trimgrit)
   {
      _braid_TriMGRITFeatureCheck(core);
   }  

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetInterval(braid_Core   core,
                   braid_Int    level,
                   braid_Int    interval_index,
                   braid_Int   *flo_ptr,
                   braid_Int   *fhi_ptr,
                   braid_Int   *ci_ptr)
{
   _braid_Grid  **grids   = _braid_CoreElt(core, grids);
   braid_Int      ilower  = _braid_GridElt(grids[level], ilower);
   braid_Int      iupper  = _braid_GridElt(grids[level], iupper);
   braid_Int      clower  = _braid_GridElt(grids[level], clower);
   braid_Int      cupper  = _braid_GridElt(grids[level], cupper);
   braid_Int      cfactor = _braid_GridElt(grids[level], cfactor);
   braid_Int      flo, fhi, ci;

   flo = ilower;
   fhi = iupper;
   ci  = -1;

   if ( _braid_IsCPoint(clower, cfactor) )
   {
      flo = clower + (interval_index-1)*cfactor + 1;
      fhi = clower + (interval_index  )*cfactor - 1;
      if (flo < ilower)
      {
         flo = ilower;
      }
      if (fhi > iupper)
      {
         fhi = iupper;
      }

      ci = clower + interval_index*cfactor;
      if (ci > cupper)
      {
         ci = -1;  /* return -1 if no C-points */
      }
   }

   *flo_ptr = flo;
   *fhi_ptr = fhi;
   *ci_ptr  = ci;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_AccessVector(braid_Core          core,
                    braid_AccessStatus  status,
                    braid_BaseVector    u)
{
   braid_App      app    = _braid_CoreElt(core, app);

   _braid_BaseAccess(core, app,  u, status);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Set initial guess at C-points, and initialize shell at F-points when using
 * shell vectors.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_InitGuess(braid_Core  core,
                 braid_Int   level)
{
   braid_App          app      = _braid_CoreElt(core, app);
   braid_Int          seq_soln = _braid_CoreElt(core, seq_soln);
   braid_Int          nrefine  = _braid_CoreElt(core,nrefine);
   _braid_Grid      **grids    = _braid_CoreElt(core, grids);
   braid_Int          ilower   = _braid_GridElt(grids[level], ilower);
   braid_Int          iupper   = _braid_GridElt(grids[level], iupper);
   braid_Int          clower   = _braid_GridElt(grids[level], clower);
   braid_Int          cupper   = _braid_GridElt(grids[level], cupper);
   braid_Int          cfactor  = _braid_GridElt(grids[level], cfactor);
   braid_Real        *ta       = _braid_GridElt(grids[level], ta);
   braid_BaseVector  *va       = _braid_GridElt(grids[level], va);

   braid_BaseVector  u;
   braid_Int         i, iu, sflag;

   if ( _braid_CoreElt(core, trimgrit) )
   {
      /* Initialize and store all points (only called from level 0) */
      for (i = ilower; i <= iupper; i++)
      {
         _braid_BaseInit(core, app, ta[i-ilower], &u);
         _braid_USetVectorRef(core, level, i, u);
      }

      return _braid_error_flag;
   }

   if ( (level == 0) && (seq_soln == 1) )
   {
      /* If first processor, grab initial condition */
      if(ilower == 0)
      {
         /* If we have already refined, then an initial init has already been done */
         if(nrefine > 0)
         {
            _braid_UGetVector(core, 0, 0, &u);    /* Get stored vector */
         }
         else
         {
            _braid_BaseInit(core, app,  ta[0], &u);
            _braid_USetVector(core, 0, 0, u, 0);
         }
         ilower += 1;
      }
      /* Else, receive point to the left */
      else
      {
         _braid_UCommInitBasic(core, 0, 1, 0, 0);     /* Post receive to the left*/
         _braid_UGetVector(core, 0, ilower-1, &u);    /* Wait on receive */
      }

      /* Set Flag so that USetVector initiates send when iupper is available */
      _braid_GridElt(grids[level], send_index)  = iupper;

      /* Initialize all points on the finest grid with sequential time marching */
      for(i = ilower; i <= iupper; i++)
      {
         _braid_Step(core, 0, i, NULL, u);       /* Step forward */
         _braid_USetVector(core, 0, i, u, 0);    /* Store: copy u into core,
                                                    sending to left if needed */
      }

      _braid_UCommWait(core, 0);                 /* Wait on comm to finish */
   }
   else if (level == 0)
   {
      if (_braid_CoreElt(core, useshell)==1)
      {
         for (i = ilower; i <=iupper; i++)
         {
            if (_braid_IsCPoint(i,cfactor))
            {
               // We are on a C-point, init full vector
               _braid_BaseInit(core, app,  ta[i-ilower], &u);
            }
            else
            {
               // We are on a F-point, init shell only
               _braid_BaseSInit(core,  app, ta[i-ilower], &u);
            }
            _braid_USetVectorRef(core, level, i, u);
         }
      }
      else
      {
         /* Only initialize the C-points on the finest grid */
         for (i = clower; i <= cupper; i += cfactor)
         {
            _braid_BaseInit(core, app,  ta[i-ilower], &u);
            _braid_USetVectorRef(core, level, i, u);
         }
      }
   }
   else
   {
      for (i = ilower; i <= iupper; i++)
      {
         _braid_UGetIndex(core, level, i, &iu, &sflag);
         if (sflag == 0) // Full point
         {
            _braid_BaseClone(core, app,  va[i-ilower], &u);
            _braid_USetVectorRef(core, level, i, u);
         }
         else if (sflag == -1) // Shell
         {
            _braid_BaseSClone(core,  app, va[i-ilower], &u);
            _braid_USetVectorRef(core, level, i, u);
         }
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Copy the initialized C-points on the fine grid, to all coarse levels.
 * Allows first down cycle to be skipped, in FMG fashion.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_CopyFineToCoarse(braid_Core  core)
{
   braid_App      app     = _braid_CoreElt(core, app);
   _braid_Grid  **grids   = _braid_CoreElt(core, grids);
   braid_Int      nlevels = _braid_CoreElt(core, nlevels);
   
   braid_Int      f_index, index, iu, is_stored, level, f_cfactor;
   braid_Int      ilower, iupper;
   braid_BaseVector   u, *va;

   for(level = 1; level < nlevels; level++)
   {

      f_cfactor = _braid_GridElt(grids[level-1], cfactor);
      iupper    = _braid_GridElt(grids[level], iupper);
      ilower    = _braid_GridElt(grids[level], ilower);
      va        = _braid_GridElt(grids[level], va);

      /* Loop over all points belonging to this processor, and if a C-point,
       * then carry out spatial coarsening and copy to ua and va */
      for (index=ilower; index<=iupper; index++)
      {
         _braid_MapCoarseToFine(index, f_cfactor, f_index);
         _braid_UGetVector(core, level-1, f_index, &u);
         _braid_Coarsen(core, level, f_index, index, u, &va[index-ilower]);
         
         _braid_BaseFree(core, app,  u);
         _braid_BaseClone(core, app,  va[index-ilower], &u);
         _braid_USetVectorRef(core, level, index, u);
         _braid_UGetIndex(core, level, index, &iu, &is_stored);
         if (is_stored == -2) /* Case where F-points are not stored, and we are not using shell vectors */
         {
            _braid_BaseFree(core, app,  u);
         }
         else if (is_stored == -1) /* This is a shell vector */
         {
            // We free the data in u, keeping the shell
            _braid_BaseSFree(core,  app, u);
         }
 
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Process the error with code ierr raised at the given line and source file
 *----------------------------------------------------------------------------*/

void
_braid_ErrorHandler(const char *filename,
                    braid_Int   line,
                    braid_Int   ierr,
                    const char *msg)
{
   _braid_error_flag |= ierr;

#ifdef braid_PRINT_ERRORS
   if (msg)
   {
      _braid_printf("braid error in file \"%s\", line %d, error code = %d - %s\n",
                    filename, line, ierr, msg);
   }
   else
   {
      _braid_printf("braid error in file \"%s\", line %d, error code = %d\n",
                    filename, line, ierr);
   }
#endif
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_ChunkSetInitialCondition(braid_Core core)
{
   
   braid_App            app       = _braid_CoreElt(core, app);
   braid_Int            ntime     = _braid_CoreElt(core, ntime);
   MPI_Comm             comm      = _braid_CoreElt(core, comm);
   braid_BufferStatus   bstatus   = (braid_BufferStatus)core;
   void           *sendbuffer;
   void           *recvbuffer;
   MPI_Request    *sendrequests;
   MPI_Request    *recvrequests;
   MPI_Status     status;
   braid_Int      size;
   braid_Int      sender, receiver;
   braid_Int      num_requests = 0;
   braid_BaseVector ulast;
   braid_BaseVector ufirst;
   braid_Int nprocs, myid;

   MPI_Comm_rank(comm, &myid);
   MPI_Comm_size(comm, &nprocs);

   /* Send last time step to first processor*/
   if (myid == nprocs - 1)     // only true on last processor 
   {
      /* Get vector at last time point of this chunk */
      _braid_UGetLast(core, &ulast);   

      /* Allocate buffer through user routine */
      _braid_BufferStatusInit( 0, 0, bstatus );
      _braid_BaseBufSize(core, app,  &size, bstatus);
      sendbuffer= malloc(size);

      /* Pack the buffer. Note that bufpack may return a size smaller than bufsize */ 
      _braid_StatusElt(bstatus, size_buffer) = size;
      _braid_BaseBufPack(core, app,  ulast, sendbuffer, bstatus);
      size = _braid_StatusElt( bstatus, size_buffer );

      /* Send the buffer to first processor*/
      num_requests = 1;
      sendrequests = _braid_CTAlloc(MPI_Request, num_requests);
      _braid_GetProc(core, 0, 0, &receiver);
      MPI_Isend(sendbuffer, size, MPI_BYTE, receiver, 0, comm, &sendrequests[0]);
   }

   /* Receive last time step from last processor */
   if (myid == 0)
   {
      /* Allocate buffer through user routine */
      _braid_BufferStatusInit( 0, 0, bstatus );
      _braid_BaseBufSize(core, app,  &size, bstatus);
      recvbuffer = malloc(size);

      /* Receive the buffer */
      num_requests = 1;
      _braid_GetProc(core, 0, ntime, &sender);
      recvrequests = _braid_CTAlloc(MPI_Request, num_requests);
      MPI_Irecv(recvbuffer, size, MPI_BYTE, sender, 0, comm, &recvrequests[0]);
   }

   /* Set the new initial condition */
   if (myid == 0)
   {
      /* Unpack the buffer into first time step */
      MPI_Wait(recvrequests, &status);
      _braid_BufferStatusInit( 0, 0, bstatus );
      _braid_BaseBufUnpack(core, app,  recvbuffer, &ufirst, bstatus);

      /* Set initial condition */
      _braid_USetVector(core, 0, 0, ufirst, 1);    // copy or move?
      // _braid_BaseFree(core, app,  ufirst);
      
      /* Free the buffer */
      _braid_TFree(recvbuffer);
      _braid_TFree(recvrequests);
   }

   /* Free the send buffer */
   if (myid == nprocs - 1) 
   {
      MPI_Wait(sendrequests, &status);
      _braid_TFree(sendbuffer);
      _braid_TFree(sendrequests);
   }

   return _braid_error_flag;
}
