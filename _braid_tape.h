
/** \file _braid_tape.h
 * \brief Define headers for tape routines.
 *
 */

#ifndef braid_tape_HEADER
#define braid_tape_HEADER

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
   STEP      = 1,
   INIT      = 2,
   CLONE     = 3,
   FREE      = 4,
   SUM       = 5,
   BUFPACK   = 6,
   BUFUNPACK = 7,
   ACCESS    = 8,

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
   braid_Real        inTime;           /*< time of the input vector */
   braid_Real        outTime;          /*< time of the output vector */
   braid_Status      status;           /*< braid's status used in step and access */
   braid_Real        sum_alpha;        /*< first coefficient of my_sum */
   braid_Real        sum_beta;         /*< second coefficient of my_sum */
   braid_Int         send_recv_rank;   /*< processor rank of sender / receiver in my_bufpack / my_bufunpack */
   braid_Int         braid_iter;       /*< iteration number of xBraid */
   braid_Int         myid;             /*< processors id */

} _braid_Action;
 


/**
 * Initialize the tape
 * Set head to NULL
 **/
void 
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
void
_braid_TapeDisplayBackwards(braid_Core core, _braid_Tape* head, void (*fctptr)(braid_Core core, void* data_ptr));

/** 
 * Evaluate the action tape in reverse order. This will clear the action tape!
 * Input: - pointer to the braid core 
 *        - pointer to the head of the action tape
 */
void
_braid_TapeEvaluate(braid_Core core);

/**
 * Call differentiated action 
 */
void
_braid_DiffCall(_braid_Action* action);

// /**
//  * Display function for a _braid_action
//  * Input: Pointer to an action 
//  */
// void
// _braid_TapeDisplayAction(braid_Core core,void* data_ptr);

// /**
//  * Display function for a _braid_vector
//  * Input: Pointer to a _braid_Vector
//  */
// void
// _braid_TapeDisplayPrimal(braid_Core core,void* data_ptr);


// void
// _braid_TapeDisplayInt(braid_Core core,void* data_ptr);


/**
 * Return the name of a _braid_Call (action name)
 */
const char* _braid_CallGetName(_braid_Call call);


#endif