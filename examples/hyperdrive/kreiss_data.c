#include <stdlib.h>
#include <math.h>
#include "kreiss_data.h"

int
init_grid_fcn(kreiss_solver *kd_, double t, grid_fcn **u_handle)
{
   grid_fcn * u_;
   int i, offset;
   const double eps=1e-12;
   
/* allocate memory */
   u_ = (grid_fcn *) malloc(sizeof(grid_fcn));
   
/* should make sure the kd_ and gf_ are not NULL */
   if (kd_ == NULL)
   {
      printf("ERROR init_grid_fcn: kd_ == NULL\n");
      return 1;
   }
      
   u_->n     = kd_->n;
   u_->sol   = malloc((kd_->n+2)*sizeof(double));
   u_->vsol_ = create_double_array_1d(3);

/* initial conditions */
   if (fabs(t - kd_->tstart) < eps)
   {
      exact1( kd_->n, u_->sol, kd_->h, kd_->amp, kd_->ph, kd_->om, t, kd_->pnr );
      bdata(u_->vsol_, kd_->amp, kd_->ph, kd_->om, t, kd_->pnr);
   }
   else /* set the grid function to zero, but could assign random values instead */
   {
      for (i=0; i< kd_->n+2; i++)
         u_->sol[i] = 0;

/* then the 3 values in the bndry ode */
#define uvsol(i) compute_index_1d(u_->vsol_, i)
      offset = kd_->n+1;
      for (i=1; i<=3; i++)
         uvsol(i) = 0;
#undef uvsol
   }
   
/* make the grid function useful outside this routine */   
   *u_handle = u_;
   
   return 0;
}

void
init_kreiss_solver(double h, double amp, double ph, double om, int pnr, int taylorbc, 
                   double L, double cfl, int nstepsset, int nsteps, double tfinal, kreiss_solver *kd_)
{
   double mxeg;
   int n;
/* # ghost points */
   const int o = 1;

   n = (int) L/h+o;
   
   if( fabs( L/(n-o) - h ) > 1e-10 )
   {
      h  = L/(n-o);
      printf("NOTE: dx changed to %e\n", h);
   }

   kd_->n = n;
   kd_->h = h;
   kd_->amp = amp;
   kd_->ph = ph;
   kd_->om = om;
   kd_->pnr = pnr;
   kd_->taylorbc = taylorbc;
   kd_->L = L;
/* ! compute time step */
   mxeg = 1.0;
   kd_->dt = cfl*h/mxeg;
   
/* difference operators are currently hard-wired for 6/3 order SBP */
   kd_->nb = 6;
   kd_->wb = 9;
   kd_->bop_ = create_double_array_2d(kd_->nb,  kd_->wb);
   kd_->bope_ = create_double_array_2d(kd_->nb, kd_->wb);
/* sbp coefficients */
   bop6g( 0.7037, kd_->bop_ );
   bop6g( 0.7037, kd_->bope_ );
/* !** Extra part of SBP operator, for ghost point SBP. */
   kd_->betapcoeff = 0.03;
/* ! bope and gh contain weights for the modified SBP-operator \tilde{D} */
   sbpghost(kd_->nb, kd_->wb, kd_->bope_, &(kd_->gh), kd_->betapcoeff);
/*!** Boundary conditions on the two boundaries */
   kd_->bcnr_ = create_int_array_1d(2);
/*   code means: 1- Dirichlet in u, extrapolate v */
/*               2- Dirichlet in v, extrapolate u*/
#define bcnr(i) compute_index_1d(kd_->bcnr_, i)   
   if( kd_->pnr == 1 || kd_->pnr == 2 )
   {
      bcnr(1) = 1;
      bcnr(2) = 2;
   }

/* bndry ode workspace variables */
   kd_->vcur_  = create_double_array_1d(3);
   kd_->veval_ = create_double_array_1d(3);
   kd_->dvdt_  = create_double_array_1d(3);

/*!** RK coefficients */
   kd_->alpha_ = create_double_array_1d(4);
   kd_->beta_  = create_double_array_1d(4);

#define alpha(i) compute_index_1d(kd_->alpha_, i)   
#define beta(i) compute_index_1d(kd_->beta_, i)   

   alpha(1) = 0;
   alpha(2) = 0.5;
   alpha(3) = 0.5;
   alpha(4) = 1.0;
   
   beta(1) = 1.0/6.0;
   beta(2) = 1.0/3.0;
   beta(3) = 1.0/3.0;
   beta(4) = 1.0/6.0;

/* grid function workspace variables */
   kd_->eval = malloc((n+2)*sizeof(double));
   kd_->current = malloc((n+2)*sizeof(double));
   kd_->rhs = malloc((n+2)*sizeof(double));
   kd_->force = malloc((n+2)*sizeof(double));

/* ! compute final time or number of time steps */
   kd_->tstart = 0;
   
   if( nstepsset )
   {
      kd_->tstop = nsteps*kd_->dt;
   }
   else
   {
      nsteps = tfinal/kd_->dt;
      kd_->dt = tfinal/nsteps;
      kd_->tstop = tfinal;
   }
   kd_->nsteps = nsteps;
   
/* always save solution for now...*/
   kd_->write=1;
/* initialize solution pointer and time */
   kd_->sol_copy = NULL;
   kd_->t_copy   = -1;
}

/* --------------------------------------------------------------------
 * Create a a copy of a vector object.
 * -------------------------------------------------------------------- */
int
copy_grid_fcn(kreiss_solver    *kd_,
              grid_fcn  *u_,
              grid_fcn **v_handle)
{
/* create new grid_fcn, copy all fields from u_ */
   grid_fcn *v_;
   int i;
   
/* allocate memory */
   v_ = (grid_fcn *) malloc(sizeof(grid_fcn));

   v_->n = u_->n;
   v_->sol   = malloc((v_->n+2)*sizeof(double));
   v_->vsol_ = create_double_array_1d(3);

/* copy the array data */
   for (i=0; i<v_->n+2; i++)
      v_->sol[i] = u_->sol[i];

#define vvsol(i) compute_index_1d(v_->vsol_, i)
#define uvsol(i) compute_index_1d(u_->vsol_, i)   
   for (i=1; i<=3; i++)
      vvsol(i) = uvsol(i);
#undef vvsol
#undef uvsol

/* assign the handle to make it useful outside this routine (call by value stuff...)*/
   *v_handle = v_;

   return 0;
}

/* --------------------------------------------------------------------
 * Destroy vector object.
 * -------------------------------------------------------------------- */
int
free_grid_fcn(kreiss_solver    *kd_,
              grid_fcn  *u_)
{
/* de-allocate everything inside u_ */
   if (u_->sol) free(u_->sol);
   delete_double_array_1d( u_->vsol_);
/* de-allocate u_ itself */
   free( u_ );

   return 0;
}

/* --------------------------------------------------------------------
 * Compute vector sum y = alpha*x + beta*y.
 * NOTE: over-writes values of y
 * -------------------------------------------------------------------- */
int
sum_grid_fcn(kreiss_solver *kd_,
             double      alpha,
             grid_fcn *x_,
             double      beta,
             grid_fcn *y_)
{
   int i;
/* make sure the grid functions have the same number of grid points */
   if (x_->n != y_->n)
   {
      printf("ERROR: incompatible sizes in sum_grid_fcn: x_->n=%i but y_->n=%i\n", x_->n, y_->n);
      return 1;
   }
   for (i=0; i<x_->n+2; i++) /* not sure it makes sense to sum the ghost point values... */
      y_->sol[i] = alpha*x_->sol[i] + beta*y_->sol[i];

#define xvsol(i) compute_index_1d(x_->vsol_, i)
#define yvsol(i) compute_index_1d(y_->vsol_, i)   
   for (i=1; i<=3; i++)
      yvsol(i) = alpha*xvsol(i) + beta*yvsol(i);
   
#undef xvsol
#undef yvsol
   
   return 0;
}

/* --------------------------------------------------------------------
 * Compute dot product of 2 grid functions
 * -------------------------------------------------------------------- */
int
dot_grid_fcn(kreiss_solver *kd_,
             grid_fcn *u_,
             grid_fcn *v_,
             double      *dot_ptr)
{
   double dot=0;
   int i;

/* make sure the grid functions have the same number of grid points */
   if (u_->n != v_->n)
   {
      printf("ERROR: incompatible sizes in dot_grid_fcn: u_->n=%i but v_->n=%i\n", u_->n, v_->n);
      return 1;
   }
/* I am NOT including the ghost point values... */
   for (i=1; i<=u_->n; i++) 
      dot += u_->sol[i] * v_->sol[i];

#define uvsol(i) compute_index_1d(u_->vsol_, i)
#define vvsol(i) compute_index_1d(v_->vsol_, i)   
/* it is not clear if the boundary variable should have the same weights as the interior solution??? */
   for (i=1; i<=3; i++)
      dot += uvsol(i)*vvsol(i);
#undef uvsol
#undef xvsol

/* make the result useful outside this routine */
   *dot_ptr = dot;

   return 0;
}


/* --------------------------------------------------------------------
 * Save a grid function to file.
 * -------------------------------------------------------------------- */


int 
save_grid_fcn(kreiss_solver *kd_,
              warp_Real t,
              grid_fcn *u_)
{
   MPI_Comm   comm   = MPI_COMM_WORLD;
   int        myid;
   /* char       filename[255]; */
   /* FILE      *file; */

   /* copy the final solution to the solver structure */
   if( kd_->write )
   {
     MPI_Comm_rank(comm, &myid);
   
     printf("Inside save_grid_fcn, myRank=%i, t=%e\n", myid, t);
     if (fabs(t-kd_->tstop)<1e-12)
     {
       printf("...copying the final solution at t=%e\n", t);
      
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

/* --------------------------------------------------------------------
 * Return buffer size for vector object buffer. Vector objects contains
 * values at every grid point plus boundary ODE data
 * -------------------------------------------------------------------- */
int
gridfcn_BufSize(kreiss_solver *kd_,
                int *size_ptr)
{
/* a grid function currently hold sol[n+2] and vsol(3) values */
   *size_ptr = kd_->n+2+3;
   
   return 0;
}

/* --------------------------------------------------------------------
 * Pack a vector object in a buffer.
 * -------------------------------------------------------------------- */
int
gridfcn_BufPack(kreiss_solver *kd_,
                grid_fcn *u_,
                void *buffer)
{
   double *dbuff = buffer;
   int i, offset;
/* first n+2 values from sol */
   for (i=0; i<kd_->n+2; i++)
      dbuff[i] = u_->sol[i];
   
/* then the 3 values in the bndry ode */
#define uvsol(i) compute_index_1d(u_->vsol_, i)
   offset = kd_->n+1;
   for (i=1; i<=3; i++)
      dbuff[offset+i] = uvsol(i);
#undef uvsol

   return 0;
}

/* --------------------------------------------------------------------
 * Allocate a grid function object and copy values from a buffer.
 * -------------------------------------------------------------------- */
int
gridfcn_BufUnpack(kreiss_solver *kd_,
             void *buffer,
             warp_Vector *u_handle)
{
   int i, offset;
   grid_fcn *u_;
   double *dbuff=buffer;
   
   u_ = (grid_fcn *) malloc(sizeof(grid_fcn));
   u_->n     = kd_->n;
   u_->sol   = malloc((kd_->n+2)*sizeof(double));
   u_->vsol_ = create_double_array_1d(3);

/* copy values from buffer */
   for (i=0; i< kd_->n+2; i++)
      u_->sol[i] = dbuff[i];

/* then the 3 values in the bndry ode */
#define uvsol(i) compute_index_1d(u_->vsol_, i)
   offset = kd_->n+1;
   for (i=1; i<=3; i++)
      uvsol(i) = dbuff[offset+i];
#undef uvsol
   
/* make the new grid function useful outside this routine */
   *u_handle = u_;
   return 0;
   
}
