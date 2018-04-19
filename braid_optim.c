/**
 *  Source file implementing the optimization routines. 
 **/


#include "braid.h"

braid_Int
braid_DriveOptimization(braid_Core                   core,    /**< braid_Core (_braid_Core) struct*/
                        braid_App                    app,
                        MPI_Comm                     comm,
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

   /* Get processor id */
   MPI_Comm_rank(comm, &myid);

   /* Get maximum number of iterations */
   braid_GetMaxOptimIter(core,&maxiter);

   /* Get stopping criterion for optimization */
   braid_GetTolOptim(core, &gnorm_tol);
   braid_IsRelTolOptim(core, &rtol);

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
      if (rtol && iter > 0)
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

      /* Design update */
      design_update(app);
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


