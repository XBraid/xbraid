/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. Email xbraid-support@llnl.gov for support.
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
 

/** \file _util.c
 * \brief Source code for utility routines. See _util.h for more documentation.
 *
 */

#include "_braid.h"
#include "_util.h"

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



void 
_braid_TapeInit(_braid_Tape* head)
{
   head = NULL;
}

_braid_Tape* 
_braid_TapePush(_braid_Tape* head, void* data_ptr)
{
   _braid_Tape* tmp = (_braid_Tape*)malloc(sizeof(_braid_Tape));
   if (tmp == NULL)
   {
      exit(0);
      printf("MALLOC ERROR!\n") ;
   }
   tmp->data_ptr = data_ptr;
   tmp->next = head;
   head = tmp;

   return head;
}

_braid_Tape* 
_braid_TapePop(_braid_Tape *head)
{
    _braid_Tape* tmp = head;
    head = head->next;
    free(tmp);

    return head;
}

braid_Int 
_braid_TapeIsEmpty(_braid_Tape* head)
{
    return head == NULL ? 1 : 0;
}


void _braid_TapeDisplay(_braid_Tape* head)
{
    _braid_Tape* current;
    current = head;
    if (current!=NULL)
    {
        printf("\n\nTape: \n\n");
        do
        {
            /* ONLY WORKS FOR THE ACTION TAPE !!! */
            _braid_Action* my_action = (_braid_Action*)(current->data_ptr);
            printf("Action: %d ", my_action->braidCall);
            printf("\n");
            current=current->next;
        }
        while (current!=NULL);
        printf("\n\n");
    }
    else
    {
        printf("Stack is empty\n");
    }
}