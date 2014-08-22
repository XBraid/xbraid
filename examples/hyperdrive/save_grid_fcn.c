/*BHEADER**********************************************************************
 * Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of XBraid.  See file COPYRIGHT for details.
 *
 * XBraid is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 ***********************************************************************EHEADER*/

#include <stdlib.h>
#include <math.h>
#include "advect_data.h"

/* --------------------------------------------------------------------
 * Save a copy of a grid function.
 * -------------------------------------------------------------------- */
int 
save_grid_fcn(advection_setup *kd_,
              grid_fcn *u_,
              braidStatus   status)
{
   MPI_Comm   comm   = MPI_COMM_WORLD;
   int        myid, level=-1, done_iterating=-1, iteration=-1;
   double residual, t;
   
   /* char       filename[255]; */
   /* FILE      *file; */

   /* copy the final solution to the solver structure */
   if( kd_->write )
   {
     MPI_Comm_rank(comm, &myid);
   
     /* printf("Inside save_grid_fcn, myRank=%i, t=%e\n", myid, t); */

     if (fabs(t-kd_->tstop)<1e-12)
     {
/* get more details about about the current level, iteration, etc */
       braidGetStatusT(status, &t);
       braidGetStatusIter(status, &iteration);
       braidGetStatusResidual(status, &residual);
       braidGetStatusDone(status, &done_iterating);
       braidGetStatusLevel(status, &level);

/* done_iterating is only well-defined on the finest level */
/*       done_iterating = (residual < kd_->braidResidualLevel || iteration >= kd_->braidMaxIter); */
       done_iterating = (residual >= 0.0 || iteration >= kd_->braidMaxIter);

       printf("Inside save_grid_fcn, final time t=%e, iter=%i, residual=%e, done iterating=%i, level=%i\n", 
              t, iteration, residual, done_iterating, level);

/* only save the final solution at the matching level */
       if (level == kd_->copy_level && done_iterating)
       {
          printf("save_grid_fcn: copying grid function, iter=%i, done iterating=%i, level=%i\n", 
                 iteration, done_iterating, level);
      
/* is there a previously saved solution that needs to be de-allocated? */
          if (kd_->sol_copy != NULL)
             free_grid_fcn(kd_, kd_->sol_copy);
          kd_->sol_copy = NULL;
      
          copy_grid_fcn(kd_, u_, &(kd_->sol_copy));
          kd_->t_copy = t;
       }
       
     }
     
   }
   return 0;
}

