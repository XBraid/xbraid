#include <stdlib.h>
#include <math.h>
#include "advect_data.h"

/* --------------------------------------------------------------------
 * Save a copy of a grid function.
 * -------------------------------------------------------------------- */
int 
save_grid_fcn(advection_setup *kd_,
              warp_Real t,
              warp_Status   status,
              grid_fcn *u_)
{
   MPI_Comm   comm   = MPI_COMM_WORLD;
   int        myid, level=-1;
   int doneflag=-1;
   /* char       filename[255]; */
   /* FILE      *file; */

   /* copy the final solution to the solver structure */
   if( kd_->write )
   {
     MPI_Comm_rank(comm, &myid);
   
     /* printf("Inside save_grid_fcn, myRank=%i, t=%e\n", myid, t); */

     if (fabs(t-kd_->tstop)<1e-12)
     {
#ifdef HD_DEBUG
       printf("...copying the final solution at t=%e\n", t);
#endif

       warp_GetStatusDone(status, &doneflag);
       warp_GetStatusLevel(status, &level);

       printf("Inside save_grid_fcn, final time t=%e, done=%i, level=%i\n", t, doneflag, level);
      
/* is there a previously saved solution that needs to be de-allocated? */
       if (kd_->sol_copy != NULL)
         free_grid_fcn(kd_, kd_->sol_copy);
       kd_->sol_copy = NULL;
      
       copy_grid_fcn(kd_, u_, &(kd_->sol_copy));
       kd_->t_copy = t;
     }
     
   }
   return 0;
}

