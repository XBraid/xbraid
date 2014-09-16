/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder schroder2@llnl.gov, Rob Falgout falgout2@llnl.gov,
 * Tzanio Kolev kolev1@llnl.gov, Ulrike Yang yang11@llnl.gov, 
 * Veselin Dobrev dobrev1@llnl.gov, et al. 
 * LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. Email schroder2@llnl.gov on how to download. 
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

