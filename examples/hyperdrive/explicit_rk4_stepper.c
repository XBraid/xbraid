#include "advect_data.h"

int
explicit_rk4_stepper(advection_setup *kd_, double t, double tend, double accuracy, grid_fcn *gf_, 
                     int *rfact_)
{
/* this is a 'my_Phi()' routine, where kd_ is of type 'my_App' and gf_ is of type 'my_Vector' */
/* take 1 time step with RK-4 */
/* currently not using accuracy or rfact_ */
   
   int i, stage;
   double tstage, tbc, bdataL, bdataR, dt;   
   grid_fcn *eval_, *current_, *dwdt_, *force_;

#define vsol(i) compute_index_1d(gf_->vsol_, i)   
#define vcur(i) compute_index_1d(current_->vsol_, i)   
#define veval(i) compute_index_1d(eval_->vsol_, i)   
#define dvdt(i) compute_index_1d(dwdt_->vsol_, i)

#define alpha(i) compute_index_1d(kd_->alpha_, i)   
#define beta(i) compute_index_1d(kd_->beta_, i)   
#define bcnr(i) compute_index_1d(kd_->bcnr_, i)   

/* get the time-step from tend and t */
   dt = tend - t;
   
#ifdef HD_DEBUG
   printf("Explict RK4 stepper: t=%e, dt=%e, h=%e, n=%i\n", t, dt, gf_->h, gf_->n);
#endif

/* allocate & initialize eval, current, dwdt, and force grid functions by copying from input gf_ */
   copy_grid_fcn( kd_, gf_, &eval_ );
   copy_grid_fcn( kd_, gf_, &current_ );
   copy_grid_fcn( kd_, gf_, &dwdt_ );
   copy_grid_fcn( kd_, gf_, &force_ );
      
/* ! 4 stages of RK-4 */
   for (stage=1; stage<=4; stage++)
   {
      tstage = t + alpha(stage)*dt;
/*  impose bc on 'eval' grid function  */
      if (kd_->taylorbc == 1)
      {
         tbc = t;
         twbndry1( &bdataL, &bdataR, stage, tbc, dt, kd_ );
      }
      else if (kd_->taylorbc == 0)
      {
/* evaluate exact boundary data for the stage */
/* set stage==1 to evaluate boundary data at t=tstage */
         twbndry1( &bdataL, &bdataR, 1, tstage, dt, kd_ );
      }
      else if (kd_->taylorbc == 3)
      {
/* get the boundary data from integrating an ODE system */
         bdataL = veval(1);
      }

/* enforce bc for the 'eval' grid function */
      assign_gp( gf_->n, eval_->sol, bdataL, bdataR, kd_->betapcoeff, gf_->h, kd_->bcnr_ );

/* evaluate dwdt */
      dwdt( gf_->n, eval_->sol, dwdt_->sol, gf_->h, kd_ );

/* should fold the force into dwdt() */         
      twforce1( gf_->n, force_->sol, tstage, gf_->h, kd_->amp, kd_->ph, kd_->om, kd_->pnr, kd_->L );

/* evaluate boundary ODE */
      dvdtbndry(eval_->vsol_, dwdt_->vsol_, kd_->amp, kd_->ph, kd_->om, tstage, kd_->pnr);

      if (stage < 4)
      {
/* compute next stage grid function 'eval'*/
         for (i=1; i<=gf_->n; i++)
         {
            eval_->sol[i] = current_->sol[i] + dt*alpha(stage+1)*(dwdt_->sol[i] + force_->sol[i]);
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
         gf_->sol[i] = gf_->sol[i] + dt*beta(stage)*( dwdt_->sol[i] + force_->sol[i] );
      }
         
/* accumulate bndry data */
      for (i=1; i<=3; i++)
      {
         vsol(i) = vsol(i) + dt*beta(stage)*dvdt(i);
      }

   }/* end stage */

/* free tmp storage */
   free_grid_fcn(kd_, current_);
   free_grid_fcn(kd_, eval_);
   free_grid_fcn(kd_, dwdt_);
   free_grid_fcn(kd_, force_);
   
   return 0;
}
