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

#include "_braid.h"
#include "_util.h"

/*----------------------------------------------------------------------------
 * Coarsen in space
 *----------------------------------------------------------------------------*/

braid_Int
_braid_Coarsen(braid_Core        core,
               braid_Int         level,    /* coarse level */
               braid_Int         f_index,  /* fine index */
               braid_Int         c_index,  /* coarse index */
               braid_BaseVector  fvector,
               braid_BaseVector *cvector)
{
   braid_App      app             = _braid_CoreElt(core, app);
   _braid_Grid  **grids           = _braid_CoreElt(core, grids);
   braid_CoarsenRefStatus cstatus = (braid_CoarsenRefStatus)core;
   braid_Int      nrefine         = _braid_CoreElt(core, nrefine);
   braid_Int      gupper          = _braid_CoreElt(core, gupper);
   braid_Int      c_ilower        = _braid_GridElt(grids[level], ilower);
   braid_Int      f_ilower        = _braid_GridElt(grids[level-1], ilower);
   braid_Real    *c_ta            = _braid_GridElt(grids[level], ta);
   braid_Real    *f_ta            = _braid_GridElt(grids[level-1], ta);

   braid_Int      c_ii = c_index-c_ilower;
   braid_Int      f_ii = f_index-f_ilower;
   
   if ( _braid_CoreElt(core, scoarsen) == NULL )
   {
      /* No spatial coarsening needed, just clone the fine vector.*/
      _braid_BaseClone(core, app,  fvector, cvector);
   }
   else
   {
      /* Call the user's coarsening routine */
      _braid_CoarsenRefStatusInit(f_ta[f_ii], f_ta[f_ii-1], f_ta[f_ii+1], 
                                  c_ta[c_ii-1], c_ta[c_ii+1],
                                  level-1, nrefine, gupper, c_index, cstatus);
      _braid_BaseSCoarsen(core, app, fvector, cvector, cstatus);
   }
   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Refine in space (basic routine)
 *----------------------------------------------------------------------------*/

braid_Int
_braid_RefineBasic(braid_Core        core,
                   braid_Int         level,    /* fine level */
                   braid_Int         c_index,  /* coarse time index */
                   braid_Real       *f_ta,     /* pointer into fine time array */
                   braid_Real       *c_ta,     /* pointer into coarse time array */
                   braid_BaseVector  cvector,
                   braid_BaseVector *fvector)
{
   braid_App              app     = _braid_CoreElt(core, app);
   braid_CoarsenRefStatus cstatus = (braid_CoarsenRefStatus)core;
   braid_Int              nrefine = _braid_CoreElt(core, nrefine);
   braid_Int              gupper  = _braid_CoreElt(core, gupper);

   if ( _braid_CoreElt(core, scoarsen) == NULL )
   {
      /* No spatial refinement needed, just clone the fine vector.*/
      _braid_BaseClone(core, app,  cvector, fvector);
   }
   else
   {
      /* Call the user's refinement routine */
      _braid_CoarsenRefStatusInit(f_ta[0], f_ta[-1], f_ta[+1], c_ta[-1], c_ta[+1],
                                  level, nrefine, gupper, c_index, cstatus);
      _braid_BaseSRefine(core,  app, cvector, fvector, cstatus);
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Refine in space
 *----------------------------------------------------------------------------*/

braid_Int
_braid_Refine(braid_Core        core,
              braid_Int         level,    /* fine level */
              braid_Int         f_index,  /* fine index */
              braid_Int         c_index,  /* coarse index */
              braid_BaseVector  cvector,
              braid_BaseVector *fvector)
{
   _braid_Grid  **grids           = _braid_CoreElt(core, grids);
   braid_Int      f_ilower        = _braid_GridElt(grids[level], ilower);
   braid_Int      c_ilower        = _braid_GridElt(grids[level+1], ilower);
   braid_Real    *f_ta            = _braid_GridElt(grids[level], ta);
   braid_Real    *c_ta            = _braid_GridElt(grids[level+1], ta);

   braid_Int      c_ii = c_index-c_ilower;
   braid_Int      f_ii = f_index-f_ilower;

   _braid_RefineBasic(core, level, c_index, &f_ta[f_ii], &c_ta[c_ii], cvector, fvector);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 * Refine in Space at every point on processor if r_space is set
 *----------------------------------------------------------------------------*/

braid_Int
_braid_FRefineSpace(braid_Core   core,
                    braid_Int   *refined_ptr)
{

   MPI_Comm      comm    = _braid_CoreElt(core, comm);
   _braid_Grid **grids   = _braid_CoreElt(core, grids);   
   
   braid_Int     r_space = _braid_CoreElt(core, r_space);
   braid_Int     ilower  = _braid_GridElt(grids[0], ilower);
   braid_Int     iupper  = _braid_GridElt(grids[0], iupper);
   braid_Real   *ta      = _braid_GridElt(grids[0], ta); 

   braid_Int     global_r_space;
   braid_Int     i, ii;
   braid_BaseVector  c_vec, f_vec;  

   if ( _braid_CoreElt(core, scoarsen) != NULL )
   {
       
      if ( r_space > 0 )  
      {         
         for ( i = ilower; i <= iupper; i++ )
         {
            ii = i-ilower;
            _braid_UGetVectorRef(core, 0, i , &c_vec);
             
            if ( c_vec != NULL )
            {
               _braid_RefineBasic(core, -1, i, &ta[ii], &ta[ii], c_vec, &f_vec);
               _braid_USetVectorRef(core, 0, i, f_vec);
            }
         }                       
      }
    
      /* Check if any refinment was completed globally. If true then refine the 
         initial time point if not done already, increase nrefine, and return 2 */
      MPI_Allreduce(&r_space, &global_r_space, 1, braid_MPI_INT, MPI_MAX, comm);
      if (global_r_space > 0)
      { 
         /* Need to make sure to refine the first point. If it is a lone C point
            on a processor r_space then can never be set for that processor */ 
         if ( ilower == 0 && r_space == 0 )
         {
            _braid_UGetVectorRef(core, 0, 0, &c_vec);
            _braid_RefineBasic(core, -1, 0, &ta[0], &ta[0], c_vec, &f_vec);
            _braid_USetVectorRef(core, 0, 0, f_vec);     
         }       
                
         *refined_ptr = 2;
         _braid_CoreElt(core, nrefine) += 1;
      }
      else
      {
         *refined_ptr = 0;
      }
   }
   else
   {
      *refined_ptr = 0;
   }
    
   /* Reset r_space */
   _braid_CoreElt(core, r_space) = 0;
   return _braid_error_flag;
}

