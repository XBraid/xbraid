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

/** \file tape.h
 * \brief Define headers for the tape routines (linked list for AD)
 *
 */

#ifndef _braid_tape_HEADER
#define _braid_tape_HEADER

#include "_braid.h"
#include "braid.h"

/**
 * 
 * C-Implementation of a linked list storing pointers to generic data
 * This structure represents one tape element, holding a pointer to data and a pointer to the next element
 * int size holds the number of all elements in the tape
 **/ 
typedef struct _braid_tape_struct
{
    int size;
    void *data_ptr;
    struct _braid_tape_struct *next;

} _braid_Tape;


/** 
 * Enumerator for identifying performed action 
 **/
typedef enum _braid_Call_enum
{
   STEP       = 1,
   INIT       = 2,
   CLONE      = 3,
   FREE       = 4,
   SUM        = 5,
   BUFPACK    = 6,
   BUFUNPACK  = 7,
   ACCESS     = 8,
   OBJECTIVET = 9

} _braid_Call;

/**
 * XBraid Action structure
 *
 * Holds information for the called user routines
 **/
typedef struct _braid_Action_struct
{
   _braid_Call       braidCall;        /**< type of the user routine */
   braid_Core        core;             /**< pointer to braid's core structure */
   braid_Real        inTime;           /**< time of the input vector */
   braid_Real        outTime;          /**< time of the output vector */
   braid_Int         inTimeIdx;        /**< index of time of input vector */
   braid_Real        sum_alpha;        /**< first coefficient of my_sum */
   braid_Real        sum_beta;         /**< second coefficient of my_sum */
   braid_Int         send_recv_rank;   /**< processor rank of sender / receiver in my_bufpack / my_bufunpack */
   braid_Int         braid_iter;       /**< iteration number of xBraid */
   braid_Int         myid;             /**< processors id */
   braid_Int         level;            /**< current level in Braid */
   braid_Int         nrefine;          /**< number of refinements done */
   braid_Int         gupper;           /**< global size of the fine grid */
   braid_Real        tol;              /**< primal stopping tolerance */      
   braid_Int         messagetype;      /**< message type, 0: for Step(), 1: for load balancing */
   braid_Int         size_buffer;      /**< if set by user, size of send buffer is "size" bytes */

} _braid_Action;
 

/**
 * Initialize the tape
 * Set head to NULL
 **/
braid_Int 
_braid_TapeInit(_braid_Tape* head);

/**
 * Push data on the tape 
 * Return pointer to head
 **/
_braid_Tape* 
_braid_TapePush(_braid_Tape* head, void* ptr);

/**
 * Pop an element from the tape 
 * Return pointer to head
 **/
_braid_Tape* 
_braid_TapePop(_braid_Tape* head);

/** 
 * Test if tape is empty
 * return 1 if tape is empty, otherwise returns 0
 **/
braid_Int 
_braid_TapeIsEmpty(_braid_Tape* head);

/**
 * Returns the number of elements in the tape
 */
braid_Int
_braid_TapeGetSize(_braid_Tape* head);

/** 
 * Display the tape in reverse order, calls the display function at each element
 * Input: - pointer to the braid core 
 *        - pointer to the display function
 */
braid_Int
_braid_TapeDisplayBackwards(braid_Core core, _braid_Tape* head, void (*fctptr)(braid_Core core, void* data_ptr));

/** 
 * Evaluate the action tape in reverse order. This will clear the action tape!
 * Input: - pointer to the braid core 
 *        - pointer to the head of the action tape
 */
braid_Int
_braid_TapeEvaluate(braid_Core core);

/**
 * Call differentiated action 
 */
braid_Int
_braid_DiffCall(_braid_Action* action);

/** 
 * Set the adjoint seed for tape evaluation, i.e., set u->bar at stored points
 * on level 0 to the values contained in core->optim->adjoints
 */
braid_Int
_braid_TapeSetSeed(braid_Core core);

/**
 * Set the pointers in tapeinput to the input of an xbraid iteration (ua).
 */
braid_Int 
_braid_TapeResetInput(braid_Core core);

/**
 * Return the name of a _braid_Call (action name)
 */
const char* _braid_CallGetName(_braid_Call call);

#endif
