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
 

/** \file _util.h
 * \brief Define headers for utility routines.
 *
 * This file contains the headers for utility routines. Essentially,
 * if a routine does not take braid_Core (or other XBraid specific structs) 
 * as an argument, then it's a utility routine.
 */

#ifndef braid_util_HEADER
#define braid_util_HEADER

#include "_braid.h"

/**
 * Project an interval onto a strided index space that contains the index
 * 'index' and has stride 'stride'.  An empty projection is represented by
 * ilower > iupper.
 **/
braid_Int
_braid_ProjectInterval( braid_Int   ilower,
                        braid_Int   iupper,
                        braid_Int   index,
                        braid_Int   stride,
                        braid_Int  *pilower,
                        braid_Int  *piupper );

/**
 * If set, print to @ref _braid_printfile and then flush.  
 * Else print to standard out.\n
 *
 * The string *format* can be multiple parameters
 * in the standard * C-format, like\n
 * `` format = '%1.2e is a format string', 1.24 ``
 **/
braid_Int
_braid_printf( const char *format, ...);

/**
 * This is a function that allows for "sane" printing
 * of information in parallel.  Currently, only 
 * myid = 0 prints, but this can be updated as needs change.
 *
 * The string *message* is printed and can be multiple parameters
 * in the standard * C-format, like\n
 * `` message = '%1.2e is a format string', 1.24 ``
 **/
braid_Int
_braid_ParFprintfFlush(FILE * file, 
                       braid_Int myid,
                       char * message, 
                       ...);

/**
 * This function finds the maximum value in a braid_Real array 
 **/
braid_Int
_braid_Max(braid_Real * array, 
           braid_Int    size,
           braid_Real * max_val);


/**
 * Copy *k* entries from *_array* into *array*.  If *k*
 * is negative, return the last *k* entries.  If 
 * positive, return the first *k* entries.  Upon
 * exit, *k* holds the number of residuals actually 
 * returned (in the case that |k| > array_len.)
 *
 * If no entries are copied, *k=0*, *array[0] = -1.0*
 **/
braid_Int
_braid_GetNEntries(braid_Real   *_array, 
                   braid_Int    array_len, 
                   braid_Int    *k_ptr, 
                   braid_Real   *array);


#endif

