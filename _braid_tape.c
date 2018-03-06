/** \file _braid_tape.c
 * \brief Source code for utility routines. See _util.h for more documentation.
 *
 */

#include "_braid.h"
#include "_braid_tape.h"

#ifndef DEBUG
#define DEBUG 0
#endif


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
      printf("MALLOC ERROR!\n") ;
      exit(1);
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


void 
_braid_TapeIterateBackwards(braid_Core core, _braid_Tape* head, void (*displayfct)(braid_Core core, void* data_ptr))
{
    _braid_Tape* current;
    current = head;
    if (current!=NULL)
    {
        do 
        {
            /* Call the display function */
            (*displayfct)(core, current->data_ptr);
            /* Move to next element */
            current = current->next;
        }
        while (current!=NULL);
    }
    else
    {
        printf("Tape is empty\n");
    }
}


void
_braid_TapeDisplayAction(braid_Core core,void* data_ptr){

    /* Get the action */    
    _braid_Action* my_action = (_braid_Action*)(data_ptr);
    /* Print action information */
    printf("%d: %s ", my_action->myid, _braid_CallGetName(my_action->braidCall));
    printf("\n");

}

void
_braid_TapeDisplayPrimal(braid_Core core,void* data_ptr)
{
    /* Get the braid_vector */
   braid_Vector vector = (braid_Vector) (data_ptr);

   /*--- Display the vector --*/
   braid_AccessStatus   astatus = (braid_AccessStatus)core;
   _braid_CoreFcn(core, access)(_braid_CoreElt(core, app), vector, astatus );
}

const char* _braid_CallGetName(_braid_Call call)
{
    switch (call)
    {
        case STEP:      return "STEP";
        case INIT:      return "INIT";
        case CLONE:     return "CLONE";
        case FREE:      return "FREE";
        case SUM:       return "SUM";
        case BUFPACK:   return "BUFPACK";
        case BUFUNPACK: return "BUFUNPACK";
        case ACCESS:    return "ACCESS";
    }
    return "Not a _braid_Call!";
}