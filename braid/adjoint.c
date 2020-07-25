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

#include "_braid.h"
#include "util.h"

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

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
      _braid_printf("ERROR: useCount < 0 !\n");
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
_braid_AdjointFeatureCheck(braid_Core core)
{

   braid_PtFcnResidual  residual  = _braid_CoreElt(core, residual);
   braid_PtFcnResidual  fullres   = _braid_CoreElt(core, full_rnorm_res);
   braid_PtFcnSCoarsen  scoarsen  = _braid_CoreElt(core, scoarsen);
   braid_PtFcnSRefine   srefine   = _braid_CoreElt(core, srefine);
   braid_PtFcnTimeGrid  tgrid     = _braid_CoreElt(core, tgrid);
   braid_PtFcnSync      sync      = _braid_CoreElt(core, sync);
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
   if ( sync != NULL )
   {
      err_char = "Sync";
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

