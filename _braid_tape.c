/** \file _braid_tape.c
 * \brief Source code for utility routines. See _util.h for more documentation.
 *
 */

#include "_braid.h"
#include "_braid_tape.h"
#include "_braid_user.h"

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
   if ( _braid_TapeIsEmpty(head) ) tmp->size = 1;
   else tmp->size = head->size + 1;
   tmp->data_ptr = data_ptr;
   tmp->next = head;
   head = tmp;

   return head;
}

_braid_Tape* 
_braid_TapePop(_braid_Tape *head)
{
    _braid_Tape* tmp = head;
    head->size = tmp->size--;
    head = head->next;
    free(tmp);

    return head;
}

braid_Int 
_braid_TapeIsEmpty(_braid_Tape* head)
{
    return head == NULL ? 1 : 0;
}

braid_Int
_braid_TapeGetSize(_braid_Tape* head)
{
    if ( _braid_TapeIsEmpty(head) ) return 0;
    else return head->size;
}


void 
_braid_TapeDisplayBackwards(braid_Core core, _braid_Tape* head, void (*displayfct)(braid_Core core, void* data_ptr))
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
_braid_TapeEvaluateAdjoint(braid_Core core)
{
    /* Get the tape */
    _braid_Tape* actiontape = _braid_CoreElt(core, actiontape);

    while ( !_braid_TapeIsEmpty(actiontape) )
    {
       /* Get the action */
       _braid_Action* action = (_braid_Action*) actiontape->data_ptr;

       /* Call the adjoint action */
       _braid_AdjointCall(core, action);

       /* Pop the action from the tape */
       actiontape = _braid_TapePop( actiontape );

       /* Free memory of the action */
       free(action);
    }
    
    /* Update the actiontape in the core */
    _braid_CoreElt(core, actiontape) = actiontape;
}

void
_braid_AdjointCall(braid_Core core, _braid_Action* action)
{
   
   /* Call the corresponding adjoint action */
   switch (action->braidCall)
   {
      case STEP : 
      {
         _braid_UserStepAdjoint(core, action);
         break;
      }
      case INIT: 
      {
         /* Do nothing */
         break;
      }
      case CLONE: 
      {
         _braid_UserCloneAdjoint(action);
         break;
      }
      case FREE: 
      {
         /* Do nothing */
         break;
      }
      case SUM: 
      {
         _braid_UserSumAdjoint(action);
         break;
      }
      case ACCESS: 
      {
         _braid_UserAccessAdjoint(core, action);
         break;
      }
      case BUFPACK: 
      {
         _braid_UserBufPackAdjoint(action, _braid_CoreElt(core, app));
         break;
      }
      case BUFUNPACK: 
      {
         _braid_UserBufUnpackAdjoint(action, _braid_CoreElt(core, app));
         break;
      } 
   }
}


void
_braid_TapeDisplayAction(braid_Core core,void* data_ptr){

    /* Get the action */    
    _braid_Action* action = (_braid_Action*)(data_ptr);
    /* Print action information */
    printf("%d: %s ", action->myid, _braid_CallGetName(action->braidCall));
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

void
_braid_TapeDisplayInt(braid_Core core,void* data_ptr)
{
    /* Get the integer*/
    int* int_ptr = (int* ) data_ptr;

   /*--- Display the integer --*/
    printf(" Integer: %d", (*int_ptr) );
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