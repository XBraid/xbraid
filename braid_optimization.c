/**
 * Optimization template - Implements a simple parallel-in-time optimization loop.
 * 
 * Calls the user's design_update routine in each outer optimization iteration. 
 * A line-search procedure should be added by the user. 
 *  
 **/

#include "braid.h"


braid_Int
braid_DriveOptimization(braid_Core                   core,    /**< braid_Core (_braid_Core) struct*/
                        braid_App                    app,
                        braid_PtFcnDesignUpdate      design_update,
                        braid_PtFcnGradientNorm      gradient_norm,
                        braid_PtFcnGradientAllreduce gradient_allreduce
                        )
{
   braid_Int  iter, maxiter;
   braid_Real objective;
   braid_Real gnorm, gnorm0;
   braid_Real gnorm_tol, tol;
   braid_Int  rtol;
   braid_Int  myid;
   braid_Real stepsize;

   /* Get processor id */
   braid_GetMyID(core, &myid);

   /* Get maximum number of iterations */
   braid_GetMaxOptimIter(core,&maxiter);

   /* Get stopping tolerance for optimization */
   braid_GetTolOptim(core, &gnorm_tol);

   /* Check if this is relative or absolute criterion */
   braid_IsRelTolOptim(core, &rtol);

   /* Get the initial stepsize */
   braid_GetStepsize(core, &stepsize);

   /* Optimization iteration */
   for (iter = 0; iter < maxiter; iter++)
   {
      /* Run adjoint XBraid to compute objective function and gradient */
      braid_Drive(core);

      /* Get the objective function value from XBraid */
      braid_GetObjective(core, &objective);

      /* Collect sensitivities from all processors */
      gradient_allreduce(app);

      /* Compute norm of the gradient */
      gradient_norm(app, &gnorm);
      if (iter == 0)
      {
         gnorm0 = gnorm;
      }

      /* Output */
      if (myid == 0) 
      {
         printf("\n %3d: Objective = %1.14e,  || Gradient || = %1.14e\n", iter, objective, gnorm);
      }
      
      /* Adjust the tolerance, if relative stopping criterion is used */
      if (rtol)
      {
         tol = gnorm_tol * gnorm0; 
      }
      else
      {
         tol = gnorm_tol;
      }

      /* Check optimization convergence */
      if (gnorm < tol)
      {
         break;
      }

      /*** Implement your favorite line-search procedure here! ***/

      /* Design update */
      design_update(app, stepsize);
   }
  

   if (iter == maxiter)
   {
      if (myid == 0) printf("\n Max. number of iterations reached.\n\n"); }
   else
   {
      if (myid == 0) printf("\n Optimization has converged.\n\n");
   }



   return 0;
}                    



