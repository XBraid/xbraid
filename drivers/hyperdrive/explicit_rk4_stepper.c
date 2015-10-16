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


#include "advect_data.h"

int
explicit_rk4_stepper(advection_setup *kd_, grid_fcn *ustop_, grid_fcn *fstop_, grid_fcn *gf_, braid_StepStatus status) 
{
/* this is a 'my_Phi()' routine, where kd_ is of type 'my_App' and gf_ is of type 'my_Vector' */
/* take 1 time step with RK-4 */
/* currently not using rfact_ */
   
   int i, stage;
   double tstage, tbc, b_data[2], dt;   
   double t, tend; 
   grid_fcn *eval_, *current_, *dwdt_;
   int level;

#define vsol(i) compute_index_1d(gf_->vsol_, i)   
#define vcur(i) compute_index_1d(current_->vsol_, i)   
#define veval(i) compute_index_1d(eval_->vsol_, i)   
#define dvdt(i) compute_index_1d(dwdt_->vsol_, i)

#define alpha(i) compute_index_1d(kd_->alpha_, i)   
#define beta(i) compute_index_1d(kd_->beta_, i)   
#define bcnr(i) compute_index_1d(kd_->bcnr_, i)   

/* get the time-step from tend and t */
#ifndef TEST_GF
   braid_StepStatusGetTstart(status, &t);
   braid_StepStatusGetTstop(status, &tend);
   braid_StepStatusGetLevel(status, &level);
#else
   t = kd_->tstart;
   tend = t + kd_->dt_fine;
   level = 0;
#endif
   
   dt = tend - t;
   
#ifdef HD_DEBUG
   grid_fcn *exact_;
   double l2, li;
   copy_grid_fcn( kd_, gf_, &exact_ );
   
   exact1( exact_,  t, kd_ );
   evaldiff( gf_, exact_, &l2, &li );

   printf("Explict RK4 stepper: tstart=%e, level=%d, dt=%e, h=%e, n=%i, l2-sol-err=%e\n", t, level, dt, gf_->h, gf_->n,l2);
#endif

/* allocate & initialize eval, current, dwdt, and force grid functions by copying from input gf_ */
   copy_grid_fcn( kd_, gf_, &eval_ );
   copy_grid_fcn( kd_, gf_, &current_ );
   copy_grid_fcn( kd_, gf_, &dwdt_ );
      
/* ! 4 stages of RK-4 */
   for (stage=1; stage<=4; stage++)
   {
      tstage = t + alpha(stage)*dt;
/*  impose bc on 'eval' grid function  */
      if (kd_->taylorbc == 1)
      {
         tbc = t;
         twbndry1( &(b_data[0]), &(b_data[1]), stage, tbc, dt, kd_ );
      }
      else if (kd_->taylorbc == 0)
      {
/* evaluate exact boundary data for the stage */
/* set stage==1 to evaluate boundary data at t=tstage */
         twbndry1( &(b_data[0]), &(b_data[1]), 1, tstage, dt, kd_ );
      }
      else if (kd_->taylorbc == 3)
      {
/* get the boundary data from integrating an ODE system */
         b_data[0] = veval(1);
      }

/* evaluate dwdt, tstage is for evaluating forcing */
      dwdt( eval_, dwdt_, tstage, b_data, kd_ );

/* evaluate boundary ODE */
      dvdtbndry(eval_, dwdt_, tstage, kd_);

      if (stage < 4)
      {
/* compute next stage grid function 'eval'*/
         for (i=1; i<=gf_->n; i++)
         {
            eval_->sol[i] = current_->sol[i] + dt*alpha(stage+1)*( dwdt_->sol[i] );
         }
/* compute next stage boundary data */
         for (i=1; i<=3; i++)
         {
            veval(i) = vcur(i) + dt*alpha(stage+1)*dvdt(i);
         }
      }
         
/* accumulate stage contributions to the solution at next time step */
      for (i=1; i<=gf_->n; i++)
      {
         gf_->sol[i] = gf_->sol[i] + dt*beta(stage)*( dwdt_->sol[i] );
      }
         
/* accumulate bndry data */
      for (i=1; i<=3; i++)
      {
         vsol(i) = vsol(i) + dt*beta(stage)*dvdt(i);
      }

   }/* end stage */

#ifdef HD_DEBUG
   exact1( exact_,  tend, kd_ );
   evaldiff( gf_, exact_, &l2, &li );

   printf("Explict RK4 stepper: tend=%e, level=%d, dt=%e, h=%e, n=%i, l2-sol-err=%e\n", tend, level, dt, gf_->h, gf_->n, l2);
   free_grid_fcn(kd_, exact_);
   
#endif
   

/* free tmp storage */
   free_grid_fcn(kd_, current_);
   free_grid_fcn(kd_, eval_);
   free_grid_fcn(kd_, dwdt_);
   
   return 0;
}
