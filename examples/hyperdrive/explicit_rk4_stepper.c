#include "kreiss_data.h"

int
explicit_rk4_stepper(kreiss_solver *kd_, double t, double tend, double accuracy, grid_fcn *gf_, 
                     int *rfact_)
{
/* this is a 'my_Phi()' routine, where kd_ is of type 'my_App' and gf_ is of type 'my_Vector' */
/* take 1 time step with RK-4 */
/* currently ignore tend. Just use the time step in kd_->dt */
/* currently not using accuracy or rfact_ */
   
   int i, stage;
   double tstage, tbc, bdataL, bdataR, dt;   

#define vsol(i) compute_index_1d(gf_->vsol_, i)   
#define vcur(i) compute_index_1d(kd_->vcur_, i)   
#define veval(i) compute_index_1d(kd_->veval_, i)   
#define dvdt(i) compute_index_1d(kd_->dvdt_, i)

#define alpha(i) compute_index_1d(kd_->alpha_, i)   
#define beta(i) compute_index_1d(kd_->beta_, i)   
#define bcnr(i) compute_index_1d(kd_->bcnr_, i)   

/* get the time-step from tend and t */
   dt = tend - t;
   
#ifdef HD_DEBUG
   printf("Explict RK4 stepper: t=%e, dt=%e, h=%e, n=%i\n", t, dt, gf_->h, gf_->n);
#endif

/* ! boundary values */
   for (i=1; i<=3; i++)
   {
      veval(i) = vsol(i);
      vcur(i)  = vsol(i);
   }
      
/* ! grid functions */
   for (i=0; i<=gf_->n+1; i++)
   {
      kd_->eval[i] = gf_->sol[i];
      kd_->current[i] = gf_->sol[i];
   }
      
/* ! 4 stages of RK-4 */
   for (stage=1; stage<=4; stage++)
   {
      tstage = t + alpha(stage)*dt;
/*  impose bc on 'eval' grid function  */
      if (kd_->taylorbc == 1)
      {
         tbc = t;
         twbndry1( 0.0, &bdataL, kd_->L, &bdataR, stage, tbc, dt, kd_->amp, kd_->ph, kd_->om, kd_->pnr);
      }
      else if (kd_->taylorbc == 0)
      {
/* evaluate exact boundary data for the stage */
/* set stage==1 to evaluate boundary data at t=tstage */
         twbndry1( 0.0, &bdataL, kd_->L, &bdataR, 1, tstage, dt, kd_->amp, kd_->ph, kd_->om, kd_->pnr);
      }
      else if (kd_->taylorbc == 3)
      {
/* get the boundary data from integrating an ODE system */
         bdataL = veval(1);
      }

/* enforce bc for the 'eval' grid function */
      bckreiss1( gf_->n, kd_->eval, bdataL, bdataR, kd_->betapcoeff, gf_->h, kd_->bcnr_ );

/* evaluate dwdt */
      dwdtkreiss1( gf_->n, kd_->eval, kd_->rhs, gf_->h, kd_->nb, kd_->wb, kd_->bop_, kd_->bope_, kd_->gh );
         
      twforce1( gf_->n, kd_->force, tstage, gf_->h, kd_->amp, kd_->ph, kd_->om, kd_->pnr, kd_->L );

/* evaluate boundary ODE */
      dvdtbndry(kd_->veval_, kd_->dvdt_, kd_->amp, kd_->ph, kd_->om, tstage, kd_->pnr);

      if (stage < 4)
      {
/* compute next stage grid function 'eval'*/
         for (i=1; i<=gf_->n; i++)
         {
            kd_->eval[i] = kd_->current[i] + dt*alpha(stage+1)*(kd_->rhs[i] + kd_->force[i]);
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
         gf_->sol[i] = gf_->sol[i] + dt*beta(stage)*( kd_->rhs[i] + kd_->force[i] );
      }
         
/* accumulate bndry data */
      for (i=1; i<=3; i++)
      {
         vsol(i) = vsol(i) + dt*beta(stage)*dvdt(i);
      }

   }/* end stage */
   return 0;
}
