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

/** \file util.c
 * \brief Source code for utility routines. See util.h for more documentation.
 *
 */

#include "_braid.h"
#include "util.h"

#ifndef DEBUG
#define DEBUG 0
#endif

/*--------------------------------------------------------------------------
 * Project an interval onto a strided index space that contains the index
 * 'index' and has stride 'stride'.  An empty projection is represented by
 * ilower > iupper.
 *--------------------------------------------------------------------------*/

braid_Int
_braid_ProjectInterval( braid_Int   ilower,
                        braid_Int   iupper,
                        braid_Int   index,
                        braid_Int   stride,
                        braid_Int  *pilower,
                        braid_Int  *piupper )
{
   braid_Int  il, iu;

   il = ilower - index;
   iu = iupper - index;

   if ( il <= 0 )
      il = (braid_Int) (il / stride);
   else
      il = (braid_Int) ((il + (stride-1)) / stride);

   if ( iu >= 0 )
      iu = (braid_Int) (iu / stride);
   else
      iu = (braid_Int) ((iu - (stride-1)) / stride);

   *pilower = index + il * stride;
   *piupper = index + iu * stride;

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
_braid_SetVerbosity(braid_Core  core,
                    braid_Int   verbose_adj)
{
   _braid_CoreElt(core, verbose_adj) = verbose_adj;

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
_braid_printf( const char *format, ...)
{
   va_list   ap;

   va_start(ap, format);
   if (_braid_printfile != NULL)
   {
      vfprintf(_braid_printfile, format, ap);
      fflush(_braid_printfile);
   }
   else
   {
      vfprintf(stdout, format, ap);
   }
   va_end(ap);

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_ParFprintfFlush(FILE       *file, 
                       braid_Int   myid,
                       const char *format,
                       ...)
{
   /* Print to file */
   if (myid == 0)
   {
      va_list   ap;
      
      va_start(ap, format);
      vfprintf(file, format, ap);
      va_end(ap);
      fflush(file);
   }

   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_Max(braid_Real  *array, 
           braid_Int    size,
           braid_Real  *max_val)
{  
   braid_Real val;
   braid_Int i;
   if(size > 0)
   {
      val = array[0];
      for(i = 1; i < size; i++)
      {
         if(array[i] > val)
         {  val = array[i]; }
      }
      *max_val = val;
   }
   
   return _braid_error_flag;
}

/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Int
_braid_GetNEntries(braid_Real   *_array, 
                   braid_Int    array_len, 
                   braid_Int    *k_ptr, 
                   braid_Real   *array)
{   
   braid_Int      start;
   braid_Int      n = *k_ptr;

   if(n < 0)
   {
      /* If negative, copy the last n residual norms */
      n = -n;
      if(array_len < n)
      {
         n = array_len;
      }
      start = array_len - n;
   }
   else
   {
      /* If positive copy the first n residual norms */
      if(array_len < n)
      {
         n = array_len;
      }
      start = 0;
   }

   if(n > 0)
   {
      memcpy(array, &(_array[start]), n*sizeof(braid_Real) );
   }
   else
   {
      array[0] = -1.0;
   }
   (*k_ptr) = n;

   return _braid_error_flag;
}


/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------*/

braid_Real
_braid_MPI_Wtime(braid_Core core, braid_Int timing_level)
{
   if( _braid_CoreElt(core, timings) >= timing_level) {
      return MPI_Wtime();
   }
   else{
      return -1.0;
   }
}
