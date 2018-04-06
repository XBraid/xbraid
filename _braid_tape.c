/** \file _braid_tape.c
 * \brief Source code for tape routines. See _braid_tape.h for more documentation.
 *
 */

#include "_braid.h"
#include "_braid_tape.h"
#include "braid_defs.h"
#include "_braid_base.h"

#ifndef DEBUG
#define DEBUG 0
#endif


braid_Int 
_braid_TapeInit(_braid_Tape* head)
{
   head = NULL;

   return _braid_error_flag;
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
   braid_Int size = 0;

   if ( ! _braid_TapeIsEmpty(head) ) 
   {
      size = head->size;
   }

   return size;
}


braid_Int
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
  
  return _braid_error_flag;
}


braid_Int
_braid_TapeEvaluate(braid_Core core)
{
   _braid_Action *action;
   _braid_Tape   *actionTape = _braid_CoreElt(core, actionTape);

   while ( !_braid_TapeIsEmpty(actionTape) )
   {
      /* Get the action */
      action = (_braid_Action*) actionTape->data_ptr;

      /* Call the differentiated action */
      _braid_DiffCall(action);

      /* Pop the action from the tape */
      actionTape = _braid_TapePop( actionTape );

     /* Free the action */
     free(action);
   }
  
   /* Update the actionTape in the core */
   _braid_CoreElt(core, actionTape) = actionTape;

   return _braid_error_flag;
}

braid_Int
_braid_DiffCall(_braid_Action* action)
{
   
   /* Call the corresponding differentiated action */
   switch (action->braidCall)
   {
      case STEP : 
      {
         _braid_BaseStep_diff(action);
         break;
      }
      case INIT: 
      {
         /* TODO: so something, if design is part of  initial conditions!!  */
         break;
      }
      case CLONE: 
      {
         _braid_BaseClone_diff(action);
         break;
      }
      case FREE: 
      {
         /* Do nothing (memory management for bar vectors is handled by shared pointer) */
         break;
      }
      case SUM: 
      {
         _braid_BaseSum_diff(action);
         break;
      }
      case ACCESS: 
      {
         /* Do nothing (access is only for output!) */
         break;
      }
      case BUFPACK: 
      {
         _braid_BaseBufPack_diff(action);
         break;
      }
      case BUFUNPACK: 
      {
         _braid_BaseBufUnpack_diff(action);
         break;
      } 
      case OBJECTIVET:
      {
          _braid_BaseObjectiveT_diff(action);
          break;
      }
   }

   return _braid_error_flag;
}


braid_Int
_braid_TapeSetSeed(braid_Core core)
{
   braid_BaseVector u_out;
   braid_Int        iclocal, sflag; 
   _braid_Grid     *fine_grid = _braid_CoreElt(core, grids)[0];
   braid_Int        clower    = _braid_GridElt(fine_grid, clower);
   braid_Int        iupper    = _braid_GridElt(fine_grid, iupper);
   braid_Int        cfactor   = _braid_GridElt(fine_grid, cfactor);
   braid_Optim      optim     = _braid_CoreElt(core, optim);
   braid_App        app       = _braid_CoreElt(core, app);


   /* Loop over all C-points */
   for (braid_Int ic=clower; ic <= iupper; ic += cfactor)
   {
      if( _braid_IsCPoint(ic, cfactor))
      {
         /* Get the vector and its local index in ua */
         _braid_UGetIndex(core, 0, ic, &iclocal, &sflag);
         _braid_UGetVectorRef(core, 0, ic, &u_out);

         /* Set the seed using the optimization adjoints */
         _braid_CoreFcn(core, sum)(app, 1.0, optim->adjoints[iclocal], 0.0, u_out->bar->userVector);
      }
   }

   return _braid_error_flag;
}

braid_Int
_braid_TapeResetInput(braid_Core core)
{
   braid_BaseVector u;
   braid_VectorBar  ubar_copy;
   braid_Int        iclocal, sflag, ic;
   _braid_Grid     *fine_grid = _braid_CoreElt(core, grids)[0];
   braid_Int        clower    = _braid_GridElt(fine_grid, clower);
   braid_Int        iupper    = _braid_GridElt(fine_grid, iupper);
   braid_Int        cfactor   = _braid_GridElt(fine_grid, cfactor);
 
   /* Loop over all C-points */
   for (ic=clower; ic <= iupper; ic += cfactor)
   {
      if( _braid_IsCPoint(ic, cfactor))
      {
         /* Get the vector and its local index in ua */
         _braid_UGetIndex(core, 0, ic, &iclocal, &sflag);
         _braid_UGetVectorRef(core, 0, ic, &u);

         /* Copy a pointer to its bar vector to the tapeinput vector */
         _braid_VectorBarCopy(u->bar, &ubar_copy );
         _braid_CoreElt(core, optim)->tapeinput[iclocal] = ubar_copy;
      }
   }

   return _braid_error_flag;
}


const char* _braid_CallGetName(_braid_Call call)
{
    switch (call)
    {
        case STEP:       return "STEP";
        case INIT:       return "INIT";
        case CLONE:      return "CLONE";
        case FREE:       return "FREE";
        case SUM:        return "SUM";
        case BUFPACK:    return "BUFPACK";
        case BUFUNPACK:  return "BUFUNPACK";
        case ACCESS:     return "ACCESS";
        case OBJECTIVET: return "OBJECTIVET";
    }
    return "Not a _braid_Call!";
}