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
 * Create a new grid object
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GridInit(braid_Core     core,
                braid_Int      level,
                braid_Int      ilower,
                braid_Int      iupper,
                _braid_Grid  **grid_ptr)
{
   _braid_Grid   *grid;
   braid_Real    *ta;

   grid = _braid_CTAlloc(_braid_Grid, 1);
   
   _braid_GridElt(grid, level)  = level;
   _braid_GridElt(grid, ilower) = ilower;
   _braid_GridElt(grid, iupper) = iupper;
   _braid_GridElt(grid, recv_index) = _braid_RecvIndexNull;
   _braid_GridElt(grid, send_index) = _braid_SendIndexNull;
   
   /* Store each processor's time slice, plus one time value to the left 
    * and to the right */
   ta = _braid_CTAlloc(braid_Real, iupper-ilower+3);
   _braid_GridElt(grid, ta_alloc) = ta;
   _braid_GridElt(grid, ta)       = ta+1;  /* shift */
   
   *grid_ptr = grid;

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GridClean(braid_Core    core,
                 _braid_Grid  *grid)
{
   braid_App      app          = _braid_CoreElt(core, app);
   braid_Int      ilower       = _braid_GridElt(grid, ilower);
   braid_Int      iupper       = _braid_GridElt(grid, iupper);
   braid_Int      nupoints     = _braid_GridElt(grid, nupoints);
   braid_BaseVector  *ua       = _braid_GridElt(grid, ua);
   braid_BaseVector  *va       = _braid_GridElt(grid, va);
   braid_BaseVector  *fa       = _braid_GridElt(grid, fa);
   braid_BaseVector  *ua_alloc = _braid_GridElt(grid, ua_alloc);
   braid_BaseVector  *va_alloc = _braid_GridElt(grid, va_alloc);
   braid_BaseVector  *fa_alloc = _braid_GridElt(grid, fa_alloc);
   
   braid_Int      ii;

   if (ua_alloc)
   {
      if (_braid_CoreElt(core, trimgrit))
      {
         /* For TriMGRIT, free ghost layers too */
         for (ii = -1; ii < nupoints+1; ii++)
         {
            if (ua[ii] != NULL)
            {
               _braid_BaseFree(core, app,  ua[ii]);
               ua[ii] = NULL;
            }
         }
      }
      else
      {
         for (ii = 0; ii < nupoints; ii++)
         {
            if (ua[ii] != NULL)
            {
               _braid_BaseFree(core, app,  ua[ii]);
               ua[ii] = NULL;
            }
         }
      }
   }
   if (va_alloc)
   {
      for (ii = -1; ii <= (iupper-ilower); ii++)
      {
         if (va[ii] != NULL)
         {
            _braid_BaseFree(core, app,  va[ii]);
            va[ii] = NULL;
         }
      }
   }
   if (fa_alloc)
   {
      for (ii = -1; ii <= (iupper-ilower); ii++)
      {
         if (fa[ii] != NULL)
         {
            _braid_BaseFree(core, app,  fa[ii]);
            fa[ii] = NULL;
         }
      }
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GridDestroy(braid_Core    core,
                   _braid_Grid  *grid)
{
   if (grid)
   {
      braid_BaseVector  *ua_alloc = _braid_GridElt(grid, ua_alloc);
      braid_Real        *ta_alloc = _braid_GridElt(grid, ta_alloc);
      braid_BaseVector  *va_alloc = _braid_GridElt(grid, va_alloc);
      braid_BaseVector  *fa_alloc = _braid_GridElt(grid, fa_alloc);

      _braid_GridClean(core, grid);

      if (ua_alloc)
      {
         _braid_TFree(ua_alloc);
      }
      if (ta_alloc)
      {
         _braid_TFree(ta_alloc);
      }
      if (va_alloc)
      {
         _braid_TFree(va_alloc);
      }
      if (fa_alloc)
      {
         _braid_TFree(fa_alloc);
      }

      _braid_TFree(grid);
   }
   return _braid_error_flag;
}

