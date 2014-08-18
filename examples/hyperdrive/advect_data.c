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

#define SQR(x) ((x)*(x))
/**< Initialize a tw_Vector function on finest temporal grid*/
int
init_grid_fcn(advection_setup *kd_, double t, grid_fcn **u_handle)
{
   grid_fcn * u_;
   int i;
   
/* allocate memory */
   u_ = (grid_fcn *) malloc(sizeof(grid_fcn));
   
/* make sure the kd_ pointer is not NULL */
   if (kd_ == NULL)
   {
      printf("ERROR init_grid_fcn: kd_ == NULL\n");
      return 1;
   }
      
   u_->n     = kd_->n_fine;
   u_->h     = kd_->h_fine;
   u_->sol   = malloc((u_->n+2)*sizeof(double));
   u_->vsol_ = create_double_array_1d(3);

/* initial conditions */
   if (fabs(t - kd_->tstart) < MY_EPS)
   {
#ifdef HD_DEBUG
      printf("Init: assigning initial data at t=tstart=%e\n", t);
#endif
      exact1( u_, t, kd_ );
      bdata( u_, t, kd_);
   }
   else /* set the grid function to zero, but could assign random values instead */
   {
#ifdef HD_DEBUG
      printf("Init: assigning zero grid function at t=%e\n", t);
#endif
      for (i=0; i< u_->n+2; i++)
         u_->sol[i] = ((double)rand())/RAND_MAX;

/* then the 3 values in the bndry ode */
#define uvsol(i) compute_index_1d(u_->vsol_, i)
      for (i=1; i<=3; i++)
         uvsol(i) = ((double)rand())/RAND_MAX;
#undef uvsol
   }
   
/* make the grid function useful outside this routine */   
   *u_handle = u_;
   
   return 0;
}

void
init_advection_solver(double h, double amp, double ph, double om, int pnr, int taylorbc, 
                      double L, double cfl, int nstepsset, int nsteps, double tfinal, 
                      double wave_speed, double viscosity, int bcLeft, int bcRight, int braidMaxIter, 
                      double braidResidualLevel, double restr_coeff, double ad_coeff, int spatial_order,
                      advection_setup *kd_)
{
   double mxeg, p1, beta, pi=M_PI, omL2pi;
   int n;
/* # ghost points */
   const int o = 1;

   n = (int) rint(L/h) + o;
   
   if( fabs( L/(n-o) - h ) > 1e-10 )
   {
      h  = L/(n-o);
      printf("NOTE: dx changed to %e\n", h);
   }

   kd_->n_fine = n;
   kd_->h_fine = h;
   kd_->amp = amp;
   kd_->ph = ph;
   kd_->pnr = pnr;
/* if periodic bc, adjust om such that om*L = 2*pi */
   if (bcLeft==Periodic || bcRight== Periodic)
   {
      omL2pi = om*L/2.0/pi;
//      om = 2*pi/L;
      printf("NOTE: Periodic bc, om*L/(2*pi)=%e, must be integer\n", omL2pi);      
   }
   kd_->om = om;
// scaled wave number
   printf("Scaled wave number om*h/pi = %e\n", om*h/pi);

   kd_->taylorbc = taylorbc;
   kd_->L = L;
   if (wave_speed < 0.0)
   {
      wave_speed *= -1;
      printf("NOTE: wave_speed changed to %e\n", wave_speed);
   }
   
   kd_->c_coeff  = wave_speed;
   kd_->nu_coeff = viscosity;
   
/* restriction stuff */
   kd_->restr_coeff = restr_coeff;

/* artificial dissipation coefficient */
   kd_->ad_coeff = ad_coeff;

/* spatial order of accuracy */
   kd_->spatial_order = spatial_order;
   
/* ! compute time step */
   if (kd_->c_coeff <= 4.0*kd_->nu_coeff / h)
   {
      mxeg = 4.0*kd_->nu_coeff / h ;
      printf("Viscous time step calc with nu/(c*h) = %e\n", kd_->nu_coeff/h/kd_->c_coeff);
   }
   else
   {
      mxeg = kd_->c_coeff;
      printf("Inviscid time step calc with nu/(c*h) = %e\n", kd_->nu_coeff/h/kd_->c_coeff);
      /* if (4.0*(kd_->nu_coeff) / h > 1.0) */
      /* { */
      /*    printf("WARNING: viscos term likely to casuse instability because nu/h = %e > 0.25\n", kd_->nu_coeff/h); */
      /* } */
   }
      
   kd_->dt_fine = cfl*h/mxeg;
   
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

   p1   = 13649.0/43200.0; /* this is the weight in the 6th order scalar product for point 1 (on the bndry)*/
   beta = kd_->betapcoeff/p1;
   
/* ! bope and gh contain weights for the modified SBP-operator \tilde{D} */
#define bope(i,j) compute_index_2d(kd_->bope_, i, j)
   kd_->gh   =             beta*( 1);
   bope(1,1) = bope(1,1) + beta*(-4);
   bope(1,2) = bope(1,2) + beta*( 6);
   bope(1,3) = bope(1,3) + beta*(-4);
   bope(1,4) = bope(1,4) + beta*( 1);
#undef bope

/* compute coefficients for diffusion term */
   kd_->nb2 = 6;
   kd_->wb2 = 9;
   kd_->bop2_ = create_double_array_2d(kd_->nb2,  kd_->wb2);
   kd_->iop2_ = create_double_array_1d(7);   

   diffusion_coeff_6( kd_->iop2_, kd_->bop2_, &(kd_->gh2), kd_->bder );

/*!** Boundary conditions on the two boundaries */
   kd_->bcnr_ = create_int_array_1d(2);
/*   code means: 1- Dirichlet in u, extrapolate v */
/*               2- Dirichlet in v, extrapolate u*/
#define bcnr(i) compute_index_1d(kd_->bcnr_, i)   
   bcnr(1) = bcLeft;
   bcnr(2) = bcRight;

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

/* ! compute final time or number of time steps */
   kd_->tstart = 0;
   
   if( nstepsset )
   {
      kd_->tstop = nsteps*kd_->dt_fine;
   }
   else
   {
      nsteps       = tfinal/kd_->dt_fine;
      if (nsteps%2 == 1)
      {
         nsteps++;
         printf("After adding an extra time step to make nsteps=%i even\n", nsteps);
      }
      
      kd_->dt_fine = tfinal/nsteps;
      kd_->tstop   = tfinal;
   }
   kd_->nsteps = nsteps;
   
/* always save solution for now...*/
   kd_->write=1;
/* save final solution at this level */
   kd_->copy_level=1; /* level=0 is the finest, level=1 is coarser, and so on */
/* keep track of braid convergence criteria */
   kd_->braidMaxIter = braidMaxIter;
   kd_->braidResidualLevel = braidResidualLevel;
   
/* initialize solution pointer and time */
   kd_->sol_copy = NULL;
   kd_->t_copy   = -1;
#undef bcnr
}

/* --------------------------------------------------------------------
 * Create a a copy of a vector object.
 * -------------------------------------------------------------------- */
int
copy_grid_fcn(advection_setup *kd_,
              grid_fcn *u_,
              grid_fcn **v_handle)
{
/* create new grid_fcn, copy all fields from u_ */
   grid_fcn *v_;
   int i;
   
/* allocate memory */
   v_ = (grid_fcn *) malloc(sizeof(grid_fcn));

   v_->n = u_->n;
   v_->h = u_->h;
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
free_grid_fcn(advection_setup    *kd_,
              grid_fcn  *u_)
{
   if (u_ == NULL) return 0;
   
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
sum_grid_fcn(advection_setup *kd_,
             double alpha,
             grid_fcn *x_,
             double beta,
             grid_fcn *y_)
{
   int i;
/* make sure the grid functions have the same number of grid points */
   if (x_->n != y_->n)
   {
      printf("ERROR: incompatible sizes in sum_grid_fcn: x_->n=%i but y_->n=%i\n", x_->n, y_->n);
      return 1;
   }
/* not sure it makes sense to sum the ghost point values, but that is the prescription */
   for (i=0; i<x_->n+2; i++) 
      y_->sol[i] = alpha*x_->sol[i] + beta*y_->sol[i];

#define xvsol(i) compute_index_1d(x_->vsol_, i)
#define yvsol(i) compute_index_1d(y_->vsol_, i)   
/* sum up the boundary ode variables too */
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
dot_grid_fcn(advection_setup *kd_,
             grid_fcn *u_,
             grid_fcn *v_,
             double *dot_ptr)
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

/* I am not including the boundary ODE variables either */
#define uvsol(i) compute_index_1d(u_->vsol_, i)
#define vvsol(i) compute_index_1d(v_->vsol_, i)   
/* it is not clear if the boundary variable should have the same weights as the interior solution??? */
   /* for (i=1; i<=3; i++) */
   /*    dot += uvsol(i)*vvsol(i); */
#undef uvsol
#undef xvsol

/* make the result useful outside this routine */
   *dot_ptr = dot;

   return 0;
}


/* --------------------------------------------------------------------
 * Return buffer size for vector object buffer. Vector objects contains
 * values at every grid point plus boundary ODE data
 * -------------------------------------------------------------------- */
int
gridfcn_BufSize(advection_setup *kd_,
                int *size_ptr)
{
/* allocate enough storage for the finest grid function */
/* a grid function currently holds sol[n+2] and vsol(3) values */

/* Add space for one double and one int (sized as a double) to hold the grid size and number of grid points? */   
   *size_ptr = (kd_->n_fine+2+3+1+1)*sizeof(double);

   return 0;
}

/* --------------------------------------------------------------------
 * Pack a vector object in a buffer.
 * -------------------------------------------------------------------- */
int
gridfcn_BufPack(advection_setup *kd_,
                grid_fcn *u_,
                void *buffer)
{
/* not all elements are doubles, but for simplicity we give all elements the size of a double */
   double *dbuff = buffer;
   int *ibuff = buffer;
   int i, offset;
/* first the size (n) */
   ibuff[0] = u_->n;
/* then the grid size */
   dbuff[1] = u_->h;
/* first n+2 values from sol */
   offset = 2;
   for (i=0; i<u_->n+2; i++)
      dbuff[i+offset] = u_->sol[i];
   
/* then the 3 values in the bndry ode */
#define uvsol(i) compute_index_1d(u_->vsol_, i)
   offset = u_->n+4;
   for (i=0; i<3; i++)
      dbuff[offset+i] = uvsol(i+1); /* uvsol(i) is base 1 */ 
#undef uvsol

   return 0;
}

/* --------------------------------------------------------------------
 * Allocate a grid function object and copy values from a buffer.
 * -------------------------------------------------------------------- */
int
gridfcn_BufUnpack(advection_setup *kd_,
             void *buffer,
             tw_Vector *u_handle)
{
   int i, offset, n;
   double h;
   grid_fcn *u_;
   double *dbuff=buffer;
   double *ibuff=buffer;
   
/* start by reading the number of grid points */
   n = ibuff[0];
/* then the grid size */
   h = dbuff[1];
   
/* now allocate space for the grid function */
   u_ = (grid_fcn *) malloc(sizeof(grid_fcn));
   u_->n     = n;
   u_->h     = h;
   
   u_->sol   = malloc((u_->n+2)*sizeof(double));
   u_->vsol_ = create_double_array_1d(3);

/* copy values from buffer */
   offset = 2;
   for (i=0; i< u_->n+2; i++)
      u_->sol[i] = dbuff[i+offset];

/* then the 3 values in the bndry ode */
#define uvsol(i) compute_index_1d(u_->vsol_, i)
   offset = u_->n+4;
   for (i=0; i<3; i++)
      uvsol(i+1) = dbuff[offset+i];
#undef uvsol
   
/* make the new grid function useful outside this routine */
   *u_handle = u_;
   return 0;
}

#define MAX(a,b) (a<b? b:a)
int
gridfcn_Coarsen(advection_setup *kd_,
                double tstart,
                double f_tminus,
                double f_tplus,
                double c_tminus,
                double c_tplus,
                grid_fcn *gf_, /* pointer to the fine grid function */
                grid_fcn **cu_handle) /* handle to the coarse grid function */
{
   grid_fcn * u_;
   int i, nf, nc, ifine;
   double dt_f, dt_c;
#define bcnr(i) compute_index_1d(kd_->bcnr_, i)   

   dt_f = MAX(f_tplus - tstart, tstart - f_tminus);
   dt_c = MAX(c_tplus - tstart, tstart - c_tminus);
   
#ifdef HD_DEBUG
   printf("Coarsen: tstart=%e, dt_f = %e, dt_c=%e, grid pts (fine)=%i\n", 
          tstart, dt_f, dt_c, gf_->n);
#endif
   
/* are the time steps the same??? */
   if (fabs(dt_f-dt_c)<MY_EPS)
   {
/* just clone the grid_fcn */
      printf("Same coarse and fine time steps=%e. Copying fine to coarse grid function\n", dt_f);
      copy_grid_fcn(kd_, gf_, cu_handle);
   }
   
/* can we deal with the requested coarsening? */
   if (fabs(dt_c/dt_f - 2.0) > MY_EPS)
   {
      printf("The requested coarsening factor dt_c/dt_f=%e is not implemented\n", dt_c/dt_f);
      return 1;
   }

/* from here on the coarsening factor equals 2 */
   
/* allocate memory for the coarse grid function */
   u_ = (grid_fcn *) malloc(sizeof(grid_fcn));
   
/* should make sure the kd_ and gf_ are not NULL */
   nf = gf_->n;
   nc = u_->n = (nf-1)/2 + 1;
   u_->h     = gf_->h * 2.0;
   u_->sol   = malloc((u_->n+2)*sizeof(double));
   u_->vsol_ = create_double_array_1d(3);

/* assign the coarse grid function */
/* u_ has nc interior points */
/* gf_ has nf interior grid points */

/* impose bc on fine grid function */
   if (bcnr(1) == Periodic)
   {
      gf_->sol[0]    = gf_->sol[nf-1];
      gf_->sol[nf+1] = gf_->sol[2];
   }

/* inject the fine grid solution into coarse grid */
   for (i=1; i<=nc; i++)
   {
      ifine = 2*i-1;
/* add 2nd undivided difference term here... */
      u_->sol[i] = gf_->sol[2*i-1] + kd_->restr_coeff*(gf_->sol[ifine-1] - 2.0*gf_->sol[ifine] + gf_->sol[ifine+1]);
   }

/* enforce boundary conditions */
/* is this really necessary? */

#define fvsol(i) compute_index_1d(gf_->vsol_, i)

   if (bcnr(1) == Periodic)
   {
      u_->sol[0] = u_->sol[nc-1];
      u_->sol[nc+1] = u_->sol[2];
   }
   else
   {
/* set 0 ghost point values */
      u_->sol[0] = 0;
      u_->sol[nc+1] = 0;
   }
   
/* should check u_->sol, esp ghost point values */
/* copy the 3 values in the bndry ode */
#define uvsol(i) compute_index_1d(u_->vsol_, i)
   for (i=1; i<=3; i++)
      uvsol(i) = fvsol(i);
#undef uvsol
#undef fvsol
   
/* make the grid function useful outside this routine */   
   *cu_handle = u_;

#undef bcnr
   return 0;   
}

int
gridfcn_Refine(advection_setup * kd_,
               double tstart,
               double f_tminus,
               double f_tplus,
               double c_tminus,
               double c_tplus,
               grid_fcn *gf_, /* pointer to the coarse grid function */
               grid_fcn **fu_handle) /* handle to the fine grid function */
{
   grid_fcn *u_;
   int i, nc, nf, ifine, ig;   
   double dt_f, dt_c;
#define bcnr(i) compute_index_1d(kd_->bcnr_, i)   

   dt_f = MAX(f_tplus - tstart, tstart - f_tminus);
   dt_c = MAX(c_tplus - tstart, tstart - c_tminus);
   
#ifdef HD_DEBUG
   printf("Refine: tstart=%e, dt_f = %e, dt_c=%e, grid pts (coarse)=%i\n", 
          tstart, dt_f, dt_c, gf_->n);
#endif

/* are the time steps the same??? */
   if (fabs(dt_f-dt_c)<MY_EPS)
   {
/* just clone the vector */
      printf("Same coarse and fine time steps=%e. Copying coarse to fine grid function\n", dt_f);
      copy_grid_fcn(kd_, gf_, fu_handle);
   }
   
/* can we deal with the requested refinement? */
   if (fabs(dt_c/dt_f - 2.0) > MY_EPS)
   {
      printf("The requested refinement factor dt_c/dt_f=%e is not implemented\n", dt_c/dt_f);
      return 1;
   }

/* from here on the refinement factor equals 2 */
   
/* allocate memory for the fine grid function */
   u_ = (grid_fcn *) malloc(sizeof(grid_fcn));
   
/* should make sure the kd_ and gf_ are not NULL */
   nc = gf_->n;
   nf = u_->n = (nc-1)*2 + 1;
   u_->h      = gf_->h * 0.5;
   u_->sol    = malloc((u_->n+2)*sizeof(double));
   u_->vsol_  = create_double_array_1d(3);

/* assign the fine grid function */

/* gf_ has nc interior grid points */
/* u_ has nf interior points */

/* inject the coinciding coarse grid solution into the fine grid */
   for (i=1; i<=nc; i++)
   {
      u_->sol[2*i-1] = gf_->sol[i];
   }
/* for now, do linear interpolation to define the intermediate fine grid points */
   /* for (i=2; i<=nf-1; i+=2) */
   /* { */
   /*    u_->sol[i] = 0.5*(u_->sol[i-1] + u_->sol[i+1]); */
   /* } */

   if (bcnr(1) == Periodic)
   {
      gf_->sol[0] = gf_->sol[nc-1];
      gf_->sol[nc+1] = gf_->sol[2];

/* fourth order interpolation */
      for (ig=1; ig<=nc-1; ig++)
      {
         ifine = 2*ig; /* this is the index on the fine mesh between ig and ig+1 */
         u_->sol[ifine] = ( -gf_->sol[ig-1] - gf_->sol[ig+2] + 9.*gf_->sol[ig] + 9.*gf_->sol[ig+1] )/16.0;
      }

/* periodic fine grid function */
      u_->sol[0] = u_->sol[nf-1];
      u_->sol[nf+1] = u_->sol[2];
   }
   else
   {
/* fourth order interpolation */
      for (ig=2; ig<=nc-2; ig++)
      {
         ifine = 2*ig; /* this is the index on the fine mesh between ig and ig+1 */
         u_->sol[ifine] = ( -gf_->sol[ig-1] - gf_->sol[ig+2] + 9.*gf_->sol[ig] + 9.*gf_->sol[ig+1] )/16.0;
      }

/*! left bndry */
      ig = 1;
      ifine = 2*ig; /* ! this is the index on the fine mesh between the coarse points ig and ig+1*/
      u_->sol[ifine] = ( 5.*gf_->sol[ig] + 15.*gf_->sol[ig+1] - 5.*gf_->sol[ig+2] + gf_->sol[ig+3] )/16.0;
   
/* ! right bndry */
      ig = nc;
      ifine = nf-1; /* ! this is the index on the fine mesh between the coarse points nxG and nxG-1*/
      u_->sol[ifine] = ( 5.*gf_->sol[ig] + 15.*gf_->sol[ig-1] - 5.*gf_->sol[ig-2] + gf_->sol[ig-3] )/16.0;
   
/* enforce boundary conditions */

#define fvsol(i) compute_index_1d(gf_->vsol_, i)

/* set 0 ghost point values */
      u_->sol[0] = 0;
      u_->sol[nf+1] = 0;
   } /* end non-periodic case */
   
/* copy the 3 values in the bndry ode */
#define uvsol(i) compute_index_1d(u_->vsol_, i)
#define cvsol(i) compute_index_1d(gf_->vsol_, i)
   for (i=1; i<=3; i++)
      uvsol(i) = cvsol(i);
#undef uvsol
#undef cvsol
   
/* make the fine grid function useful outside this routine */   
   *fu_handle = u_;

   return 0;
#undef bcnr
}


