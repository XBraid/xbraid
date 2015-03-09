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


#ifndef ADVECT_DATA
#define ADVECT_DATA

#include "c_array.h"
#include "braid.h"

#define MY_EPS 1e-12

/* define HD_DEBUG to get printouts from various user defined functions*/
/* #define HD_DEBUG */

enum bcType {Periodic, Dirichlet, Extrapolation};

typedef struct _braid_Vector_struct
{
   double *sol;
   int n;
   double h;     /* grid size */
   double_array_1d *vsol_;
} grid_fcn;

typedef struct _braid_App_struct
{
   int n_fine;     /* number of grid points in the finest grid */
   double h_fine;  /* finest grid size */
   double dt_fine; /* time step satisfying cfl constraint on finest grid */
   double amp;     /* testing parameter (amplitude) */
   double ph;      /* testing parameter (phase) */
   double om;      /* testing parameter (omega) */
   int pnr;        /* problem number (1 or 2) */
   int taylorbc;   /* treatment of bndry data at intermediate RK stages */
   int nb;         /* first dimension of bop and bope arrays */
   int wb;         /* second dimension of bop and bope arrays */
   double_array_2d *bop_;  /* regular SBP coefficients */
   double_array_2d *bope_; /* extended SBP coefficients */
   double gh;              /* ghost point coefficient */
   int nb2;         /* first dimension of bop2 array */
   int wb2;         /* second dimension of bop2 array */
   double_array_2d *bop2_;  /* SBP coefficients for 2nd derivative */
   double_array_1d *iop2_;  /* interior coefficients for 2nd derivative */
   double gh2;              /* ghost point coefficient for 2nd derivative */
   double bder[7];          /* coefficients for 6th order boundary derivative */
   double L;               /* length of 1-d domain */
   double c_coeff;         /* wave speed */
   double nu_coeff;        /* viscosity */
   double betapcoeff;      /* coefficient for Dirichlet data with ghost point */
   double restr_coeff;     /* coeff for undivided 2nd difference in restriction operator */
   double ad_coeff;        /* coeff for artificial damping */
   int_array_1d *bcnr_;    /* boundary coefficient number */
   double_array_1d *alpha_, *beta_; /* RK-4 coefficients */
   int spatial_order;      /* spatial order of accuracy */

/* braid specific stuff */
   int braidMaxIter;
   double braidResidualLevel;
   int copy_level; /* copy the solution at this level */
   grid_fcn *sol_copy; /* assigned by the call-back routine save_grid_fcn() */
   double t_copy;
   
   int write; /* flag to tell braid if it should save grid functions to file */
   double tstart, tstop;
   int nsteps;
   
} advection_setup;

/* fcn prototypes */
int
init_grid_fcn(advection_setup *kd_, double t, grid_fcn **u_handle);

void
init_advection_solver(double h, double amp, double ph, double om, int pnr, int taylorbc, 
                      double L, double cfl, int nstepsset, int nsteps, double tfinal, 
                      double wave_speed, double viscosity, int bcLeft, int bcRight, int braidMaxIter, 
                      double braidResidualLevel, double restr_coeff, double ad_coeff, int spatial_order,
                      advection_setup *kd_);
int
explicit_rk4_stepper(advection_setup *kd_, grid_fcn *gf_, braid_PhiStatus status);

int
copy_grid_fcn(advection_setup    *kd_,
              grid_fcn  *u_,
              grid_fcn **v_handle);
int
free_grid_fcn(advection_setup    *kd_,
              grid_fcn  *u_);
int
sum_grid_fcn(advection_setup *kd_,
             double      alpha,
             grid_fcn *x_,
             double      beta,
             grid_fcn *y_);
int
norm_grid_fcn(advection_setup *kd_,
              grid_fcn *u_,
              grid_fcn *v_,
              double      *norm_ptr);
int
save_grid_fcn(advection_setup *kd_,
              grid_fcn *u_,
              braid_AccessStatus   status);
int
gridfcn_BufSize(advection_setup *kd_,
                int *size_ptr);
int
gridfcn_BufPack(advection_setup *kd_,
                grid_fcn *u_,
                void *buffer,
                int *size_ptr);
int
gridfcn_BufUnpack(advection_setup *kd_,
                  void *buffer,
                  braid_Vector *u_handle);
int
gridfcn_Refine(advection_setup * kd_,
               grid_fcn *cu_,
               grid_fcn **fu_handle,
               braid_CoarsenRefStatus status);

int
gridfcn_Coarsen(advection_setup *kd_,
                grid_fcn *fu_,
                grid_fcn **cu_handle,
                braid_CoarsenRefStatus status);

void
exact1( grid_fcn *w, double t, advection_setup *kd_);
void
exact_t( grid_fcn *w, double t, advection_setup *kd_);
void
exact_x( grid_fcn *w, double t, advection_setup *kd_);
void
exact_xx( grid_fcn *w, double t, advection_setup *kd_);
void 
bdata( grid_fcn *w, double t, advection_setup *kd_);
void 
bop6g(double t, double_array_2d *q06_ );
void 
diffusion_coeff_4( double_array_1d *iop2_, double_array_2d *bop2_, double *gh2, double bder[5] );
void 
diffusion_coeff_6( double_array_1d *iop2_, double_array_2d *bop2_, double *gh2, double bder[7] );
void
twbndry1( double *bdataL, double *bdataR, int stage, double t, double dt, advection_setup *kd_ );
void
assign_gp( grid_fcn *w, double bdataL, double bdataR, advection_setup *kd_ );
void
dwdt( grid_fcn *w, grid_fcn *dwdt, double t, double bdata[2], advection_setup *kd_ );
void
dwdx( grid_fcn *w, grid_fcn *dwdt, advection_setup *kd_ );
void
d2wdx2( grid_fcn *w, grid_fcn *wxx, advection_setup *kd_ );
void
dvdtbndry(grid_fcn *w, grid_fcn *dwdt, double t, advection_setup *kd_);
void
evaldiff( grid_fcn *w, grid_fcn *we, double *l2, double*li );
void
evalnorm( grid_fcn *w, double *l2, double *li );
void
adterm( grid_fcn *w, grid_fcn *dpwdxp, advection_setup *kd_ );
/* end propotypes */
   
#endif
