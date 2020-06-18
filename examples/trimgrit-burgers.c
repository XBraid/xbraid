/*BHEADER**********************************************************************
 * Written by Isaiah Meyers, Joseph Munar, Eric Neville, Tom Overman
 * 
 * This file is part of XBraid. For support, post issues to the XBraid Github page.
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

 /**
 * Example:       trimgrit-burgers.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make trimgrit-burgers
 *
 * Description:   Solves a nonlinear optimal control problem in time-parallel:
 * 
 *                min   0.5\int_0^T \int_0^1 (u(x,t)-u0(x))^2+alpha v(x,t)^2 dxdt
 * 
 *                s.t.  du/dt + u*du/dx - nu d^2u/dx^2 = v(x,t)
 *                      u(0,t)=u(1,t)=0
 *                      u(x,0)=u0(x)
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "braid.h"
#include "braid_test.h"
#define PI 3.14159265
#define g(dt,dx) dt/(2*dx)
//#define g(dt,dx) 0.0
#define b(dt,dx,nu) nu*dt/(dx*dx)

// RDF HACK
#include "_braid.h"
#include "util.h"

/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/

typedef struct _braid_App_struct
{
   int      myid;        /* Rank of the processor */
   double   alpha;       /* Relaxation parameter for objective function, v(x,t) */
   double   nu;          /* Diffusion coefficent, which we take to be large */
   int      ntime;       /* Total number of time-steps (starting at time 0) */
   int      mspace;      /* Number of space points included in our state vector */
                         /* So including boundaries we have M+2 space points */
   double   dx;          /* Spatial mesh spacing */

   int      ilower;      /* Lower index for my proc */
   int      iupper;      /* Upper index for my proc */
   int      npoints;     /* Number of time points on my proc */

   double **w;           /* Adjoint vectors at each time point on my proc */
   double  *u0;          /* Initial condition */
   double  *scr;         /* Scratch space (enough for 3 spatial vectors) */

   FILE    *ufile, *vfile, *wfile;

} my_App;

/* Define the state vector at one time-step */
typedef struct _braid_Vector_struct
{
   double *uvals;     /* Holds the R^M state vector */
   double *vvals;     /* Holds the R^M control vector */
   double *wvals;     /* Holds the R^M adjoint vector */

} my_Vector;

/*--------------------------------------------------------------------------
 * Vector utility routines
 *--------------------------------------------------------------------------*/

void
vec_create(int size, double **vec_ptr)
{
   *vec_ptr = (double*) malloc( size*sizeof(double) );
}

void
vec_destroy(double *vec)
{
   free(vec);
}

/*------------------------------------*/

void
vec_copy(int size, double *invec, double *outvec)
{
   int i;
   for (i = 0; i < size; i++)
   {
      outvec[i] = invec[i];
   }
}

/*------------------------------------*/

void
vec_axpy(int size, double alpha, double *x, double beta, double *y)
{
   int i;
   for (i = 0; i < size; i++)
   {
      y[i] = alpha*x[i] + beta*y[i];
   }
}

/*------------------------------------*/

void
vec_scale(int size, double alpha, double *x)
{
   int i;
   for (i = 0; i < size; i++)
   {
      x[i] = alpha*x[i];
   }
}

/*--------------------------------------------------------------------------
 * KKT component routines
 *--------------------------------------------------------------------------*/

/* This computes u <-- A^{-1} u for tridiagonal matrix A = [al_{i-1} ac_i au_i]
 * by solving the system Ax = u.  Note: the diagonals 'al' and 'ac' are modified
 * during the solve. */

void
apply_AInverse(int M, double *ac, double *al, double *au, double *u)
{   
   double *v = al;
   double  li, lisave;
   int     i;

   /* Compute A = LU and solve Lv = u */
   lisave = al[0]/ac[0];
   v[0] = u[0];                             // Note: this modifies al[0]
   for(i = 1; i < M; i++)
   {
      li = lisave;
      ac[i] = ac[i] - li*au[i-1];
      lisave = al[i]/ac[i];
      v[i] = u[i] - li*v[i-1];              // Note: this modifies al[i]
   }

   /* Now solve Uu = v */ 
   u[M-1] = v[M-1]/ac[M-1];
   for (i = M-2; i >= 0; i--)
   {
      u[i] = (v[i] - au[i]*u[i+1])/ac[i];      
   }
}

/*------------------------------------*/

/* This computes v <-- Au for tridiagonal matrix A = [al_{i-1} ac_i au_i] */

void
apply_A(int M, double *ac, double *al, double *au, double *u, double *v)
{   
   int  i;

   v[0] = ac[0]*u[0] + au[0]*u[1];
   for(i = 1; i < M-1; i++)
   {
      v[i] = ac[i]*u[i] + al[i-1]*u[i-1] + au[i]*u[i+1];
   }
   v[M-1] = ac[M-1]*u[M-1] + al[M-2]*u[M-2];
}

/*------------------------------------*/

/* Compute tridiagonal Jacobian matrix A'(u) = [al_{i-1} ac_i au_i] */

void
compute_A(double dt, double dx, double nu, int M, double *u,
          double *ac, double *al, double *au)
{
   double  beta  = b(dt,dx,nu);
   double  gamma = g(dt,dx);
   int     i;

#if 0

   /* This gives advection-diffusion (and matches advec-diff-rms2.c) */
   ac[0] = (1 + 2*beta);
   al[0] = -beta - gamma;
   au[0] = -beta + gamma;
   for (i = 1; i < M-1; i++)
   {
      ac[i] = (1 + 2*beta);
      al[i] = -beta - gamma;
      au[i] = -beta + gamma;
   }
   ac[M-1] = (1 + 2*beta);
   al[M-2] = -beta - gamma;

//   /* This sets A to the identity */
//   for (i = 0; i < M; i++)
//   {
//      ac[i] = 1.0;
//      al[i] = 0.0;
//      au[i] = 0.0;
//   }

#else

   /* Lax-Friedrichs flux approximation (see adaptive spatial coarsening paper) */
   {
      double gamma2, absu, absum1, absup1;
      gamma2 = gamma/2;
      absu   = fabs(u[0]);
      absup1 = fabs(u[1]);
      ac[0] = (1 + 2*beta) + gamma2*(2*absu + absup1);
      al[0] = -beta + gamma2*(-u[0] - absup1 - absu);
      au[0] = -beta + gamma2*( u[1] - absup1 - absu);
      for (i = 1; i < M-1; i++)
      {
         absum1 = absu;
         absu   = absup1;
         absup1 = fabs(u[i+1]);
         ac[i] = (1 + 2*beta) + gamma2*(2*absu + absum1 + absup1);
         al[i] = -beta + gamma2*(-u[i]   - absup1 - absu);
         au[i] = -beta + gamma2*( u[i+1] - absup1 - absu);
      }
      absum1 = absu;
      absu   = absup1;
      ac[M-1] = (1 + 2*beta) + gamma2*(2*absu + absum1);
   }

#endif

//   /* This BTCS scheme is unstable due to centered difference on dx term.
//    *
//    * With u = (-a, a, -a, a, ...), the M x M tri-diagonal matrix (M even)
//    * approaches the following singular matrix as a -> infinity:
//    *
//    *    | 1 -1          |
//    *    |-1  0  1       |
//    *    |    1  0 -1    |
//    *    |         ...   |
//    *    |         -1  1 |
//    *
//    * */
//
//   ac[0] = (1 + 2*beta) + gamma*u[1];
//   al[0] = -beta - gamma*u[1];
//   au[0] = -beta + gamma*u[0];
//   for (i = 1; i < M-1; i++)
//   {
//      ac[i] = (1 + 2*beta) + gamma*(u[i+1] - u[i-1]);
//      al[i] = -beta - gamma*u[i+1];
//      au[i] = -beta + gamma*u[i];
//   }
//   ac[M-1] = (1 + 2*beta) - gamma*u[M-2];
}

/*------------------------------------*/

/* Apply Phi(u,g) by solving the nonlinear system A(u) = g with Newton to some
 * tolerance or with a fixed number of iterations. */

void
apply_Phi(double dt, double dx, double nu, int M, double *g, double *u, double *scr)
{
   double *ac    = scr;
   double *al    = scr + M;
   double *au    = scr + 2*M;
   double *du    = scr + 3*M;
//   double  tol   = 0.0;          // TODO: Pass this in
   int     maxit = 1;            // TODO: Pass this in
   int     iter;

   for (iter = 0; iter < maxit; iter++)
   {
      /* Compute Jacobian A'(u) */
      compute_A( dt, dx, nu, M, u, ac, al, au );

      /* Compute F(u) = A(u) - g (store in 'du') */
      apply_A( M, ac, al, au, u, du );                  // note: A(u) = A'(u)u for this problem
      vec_axpy(M, -1.0, g, 1.0, du);

      /* Apply F'(u)^{-1} F(u) :  note that F'(u) = A'(u) */
      apply_AInverse( M, ac, al, au, du );

      /* Update u = u - du */
      vec_axpy(M, -1.0, du, 1.0, u);
   }
}

/*------------------------------------*/

/* Apply Du(Phi)(u) (or its adjoint) by computing w <-- A'(u)^{-1} w */

void
apply_DuPhi(double dt, double dx, double nu, int M, double *u, double *w,
            int adjoint, double *scr)
{

   double *ac    = scr;
   double *al    = scr + M;
   double *au    = scr + 2*M;

   /* Compute Jacobian A'(u) */
   compute_A( dt, dx, nu, M, u, ac, al, au );

   if (adjoint)
   {
      /* Apply A'(u)^{-T} w */
      apply_AInverse( M, ac, au, al, w );   // note order of 'au' and 'al'
   }
   else
   {
      /* Apply A'(u)^{-1} w */
      apply_AInverse( M, ac, al, au, w );
   }
}

/*------------------------------------*/

/* Apply Dv(Phi)(u) (or its adjoint) by computing w <-- dt A'(u)^{-1} w */

void
apply_DvPhi(double dt, double dx, double nu, int M, double *u, double *w,
            int adjoint, double *scr)
{

   double *ac    = scr;
   double *al    = scr + M;
   double *au    = scr + 2*M;

   /* Compute Jacobian A'(u) */
   compute_A( dt, dx, nu, M, u, ac, al, au );

   if (adjoint)
   {
      /* Apply A'(u)^{-T} w */
      apply_AInverse( M, ac, au, al, w );   // note order of 'au' and 'al'
   }
   else
   {
      /* Apply A'(u)^{-1} w */
      apply_AInverse( M, ac, al, au, w );
   }
   vec_scale(M, dt, w);
}

/*------------------------------------*/

/* Apply Gu(J)(u) by computing u <-- dx dt (u - u0) */

void
apply_GuObjective(double dt, double dx, double *u0, int M, double *u)
{
   vec_axpy(M, -1.0, u0, 1.0, u);
   vec_scale(M, dx*dt, u);
}

/*------------------------------------*/

/* Apply Gv(J)(v) by computing v <-- alpha dx dt v */

void
apply_GvObjective(double dt, double dx, double alpha, int M, double *v)
{
   vec_scale(M, alpha*dx*dt, v);
}

/*------------------------------------*/

int
apply_TriResidual(my_App     *app,
                  my_Vector  *uleft,
                  my_Vector  *uright,
                  my_Vector  *f,
                  my_Vector  *r,
                  int         homogeneous,
                  double      dt)
{
   double  nu     = (app->nu);
   double  alpha  = (app->alpha);
   int     mspace = (app->mspace);
   double  dx     = (app->dx);
   double *u0     = (app->u0);
   double *scr    = (app->scr);

   double *rutmp, *rvtmp, *rwtmp, *vectmp;

   /* Create temporary vectors */
   vec_create(mspace, &rutmp);
   vec_create(mspace, &rvtmp);
   vec_create(mspace, &rwtmp);
   vec_create(mspace, &vectmp);

   /* Grad-u of Lagrangian residual equation */
   vec_copy(mspace, (r->uvals), rutmp);
   apply_GuObjective(dt, dx, u0, mspace, rutmp);
   vec_axpy(mspace, 1.0, (r->wvals), 1.0, rutmp);
   if (uright != NULL)
   {
      vec_copy(mspace, (uright->wvals), vectmp);
      apply_DuPhi(dt, dx, nu, mspace, (uright->uvals), vectmp, 1, scr);
      vec_axpy(mspace, -1.0, vectmp, 1.0, rutmp);
   }

   /* Grad-v of Lagrangian residual equation */
   vec_copy(mspace, (r->vvals), rvtmp);
   apply_GvObjective(dt, dx, alpha, mspace, rvtmp);
   vec_copy(mspace, (r->wvals), vectmp);
   apply_DvPhi(dt, dx, nu, mspace, (r->uvals), vectmp, 1, scr);  // note: u needed for implicit
   vec_axpy(mspace, -1.0, vectmp, 1.0, rvtmp);

   /* Grad-w of Lagrangian residual equation */
   vec_copy(mspace, (r->vvals), vectmp);
   vec_scale(mspace, dt, vectmp);
   if (uleft != NULL)
   {
      vec_axpy(mspace, 1.0, (uleft->uvals), 1.0, vectmp);
   }
   else
   {
      vec_axpy(mspace, 1.0, u0, 1.0, vectmp);
   }
   vec_copy(mspace, (r->uvals), rwtmp);
   apply_Phi(dt, dx, nu, mspace, vectmp, rwtmp, scr);
   vec_axpy(mspace, 1.0, (r->uvals), -1.0, rwtmp);

   /* Subtract rhs f on coarse grids */
   if (f != NULL)
   {
      vec_axpy(mspace, -1.0, (f->uvals), 1.0, rutmp);
      vec_axpy(mspace, -1.0, (f->vvals), 1.0, rvtmp);
      vec_axpy(mspace, -1.0, (f->wvals), 1.0, rwtmp);
   }

   /* Update r */
   vec_copy(mspace, rutmp, (r->uvals));
   vec_copy(mspace, rvtmp, (r->vvals));
   vec_copy(mspace, rwtmp, (r->wvals));

   /* Destroy temporary vectors */
   vec_destroy(rutmp);
   vec_destroy(rvtmp);
   vec_destroy(rwtmp);
   vec_destroy(vectmp);

   return 0;
}   

/*--------------------------------------------------------------------------
 * TriMGRIT wrapper routines
 *--------------------------------------------------------------------------*/

/* Compute A(u) - f */

int
my_TriResidual(braid_App       app,
               braid_Vector    uleft,
               braid_Vector    uright,
               braid_Vector    f,
               braid_Vector    r,
               braid_Int       homogeneous,
               braid_TriStatus status)
{
   double  t, tprev, tnext, dt;
   int     level, index;

   braid_TriStatusGetTriT(status, &t, &tprev, &tnext);
   braid_TriStatusGetLevel(status, &level);
   braid_TriStatusGetTIndex(status, &index);

   /* Get the time-step size */
   if (t < tnext)
   {
      dt = tnext - t;
   }
   else
   {
      dt = t - tprev;
   }

   /* Compute residual */
   apply_TriResidual(app, uleft, uright, f, r, homogeneous, dt);

   return 0;
}   

/*------------------------------------*/

/* Solve A(u) = f */

int
my_TriSolve(braid_App       app,
            braid_Vector    uleft,
            braid_Vector    uright,
            braid_Vector    f,
            braid_Vector    u,
            braid_Int       homogeneous,
            braid_TriStatus status)
{
   double  nu     = (app->nu);
   double  alpha  = (app->alpha);
   int     mspace = (app->mspace);
   double  dx     = (app->dx);
   double *u0     = (app->u0);
   double *scr    = (app->scr);

   double  t, tprev, tnext, dt;
   double *utmp, *vtmp, *wtmp, *rtmp, *vectmp, scale;
   int     xrelax, iter;

   // RDF HACK BEGIN
   int     level, index;
   braid_TriStatusGetLevel(status, &level);
   braid_TriStatusGetTIndex(status, &index);
   // RDF HACK END

   braid_TriStatusGetXRelax(status, &xrelax);

   /* Get the time-step size */
   braid_TriStatusGetTriT(status, &t, &tprev, &tnext);
   if (t < tnext)
   {
      dt = tnext - t;
   }
   else
   {
      dt = t - tprev;
   }

   /* Create temporary vectors */
   vec_create(mspace, &utmp);
   vec_create(mspace, &vtmp);
   vec_create(mspace, &wtmp);
   vec_create(mspace, &rtmp);
   vec_create(mspace, &vectmp);

   /* Initialize temporary solution vectors */
   vec_copy(mspace, (u->uvals), utmp);
   vec_copy(mspace, (u->vvals), vtmp);
   vec_copy(mspace, (u->wvals), wtmp);
   
   for (iter = 0; iter < 1; iter++)
   {
   /* BEGIN residual calc for Grad-u of Lagrangian at time (i-1) */
   if (!xrelax)
   {
      if (uleft != NULL)
      {
         vec_copy(mspace, (uleft->uvals), rtmp);
         apply_GuObjective(dt, dx, u0, mspace, rtmp);
         vec_axpy(mspace, 1.0, (uleft->wvals), 1.0, rtmp);
         vec_copy(mspace, (u->wvals), vectmp);
         apply_DuPhi(dt, dx, nu, mspace, (u->uvals), vectmp, 1, scr);
         vec_axpy(mspace, -1.0, vectmp, 1.0, rtmp);
         /* RDF HACK BEGIN subtract the tau-correction RHS here */
         if (level > 0)
         {
            _braid_Grid      **grids = _braid_CoreElt((braid_Core)status, grids);
            braid_BaseVector  *fa    = _braid_GridElt(grids[level], fa);
            braid_Vector       fleft = fa[index-1]->userVector;
            vec_axpy(mspace, -1.0, (fleft->uvals), 1.0, rtmp);
         }
         /* RDF HACK END   subtract the tau-correction RHS here */
         vec_scale(mspace, 1/(dx*dt), rtmp);                            // see GuObjective
         apply_DuPhi(dt, dx, nu, mspace, (u->uvals), rtmp, 0, scr);
         vec_scale(mspace, -1.0, rtmp);
      }
      else
      {
         vec_copy(mspace, (u->uvals), rtmp);
         vec_scale(mspace, 0.0, rtmp);
      }
   }
   /* END   residual calc at time (i-1) */

   /* Compute residual */
   my_TriResidual(app, uleft, uright, f, u, homogeneous, status);

   if (!xrelax)
   {
      /* Compute w_i using a Schur-complement approach on the 4x4 block system.
       * Assume the residual for the first equation is zero.  Use diagonal of
       * Schur-complement to scale residual.
       *
       *   Schur_residual = U_i^{-1}*ures - D_i*V_i^{-1}*vres - wres
       */
      vec_copy(mspace, (u->uvals), vectmp);
      vec_scale(mspace, 1/(dx*dt), vectmp);                        // see GuObjective
      vec_axpy(mspace, 1.0, vectmp, 1.0, rtmp);
      vec_copy(mspace, (u->vvals), vectmp);
      vec_scale(mspace, 1/(alpha*dx*dt), vectmp);                  // see GvObjective
      apply_DvPhi(dt, dx, nu, mspace, (u->uvals), vectmp, 0, scr); // note: u needed for implicit
      vec_axpy(mspace, -1.0, vectmp, 1.0, rtmp);
      vec_axpy(mspace, -1.0, (u->wvals), 1.0, rtmp);
      if (uleft != NULL)
      {
         scale = (1/(dx*dt))*(2 + dt*dt/alpha);
      }
      else
      {
         scale = (1/(dx*dt))*(1 + dt*dt/alpha);
      }
      vec_axpy(mspace, -1.0/scale, rtmp, 1.0, wtmp);
   }
   else
   {
      /* Compute v_i from w_i */
      vec_copy(mspace, (u->vvals), rtmp);
      vec_scale(mspace, 1/(alpha*dx*dt), rtmp);                    // see GvObjective
      vec_axpy(mspace, -1.0, rtmp, 1.0, vtmp);

      /* Compute u_i from w_i */
      vec_copy(mspace, (u->uvals), rtmp);
      vec_scale(mspace, 1/(dx*dt), rtmp);                          // see GuObjective
      vec_axpy(mspace, -1.0, rtmp, 1.0, utmp);
   }

   /* Update u */
   vec_copy(mspace, wtmp, (u->wvals));
   vec_copy(mspace, vtmp, (u->vvals));
   vec_copy(mspace, utmp, (u->uvals));
   }

   /* no refinement */
   braid_TriStatusSetRFactor(status, 1);

   /* Destroy temporary vectors */
   vec_destroy(utmp);
   vec_destroy(vtmp);
   vec_destroy(wtmp);
   vec_destroy(rtmp);
   vec_destroy(vectmp);

   return 0;
}   

/*------------------------------------*/

/* This is only called from level 0 */

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   int i;
   my_Vector *u;
   int mspace = (app->mspace);

   /* Allocate the vector */
   u = (my_Vector *) malloc(sizeof(my_Vector));
   vec_create(mspace, &(u->uvals));
   vec_create(mspace, &(u->vvals));
   vec_create(mspace, &(u->wvals));

   for (i = 0; i <= mspace-1; i++)
   {
      /* In reverse order to match against Schur approach */
      u->wvals[i] = ((double)braid_Rand())/braid_RAND_MAX;
      u->vvals[i] = ((double)braid_Rand())/braid_RAND_MAX;
      u->uvals[i] = ((double)braid_Rand())/braid_RAND_MAX;
   }

   *u_ptr = u;

   return 0;
}

/*------------------------------------*/

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   int mspace = (app->mspace);
   my_Vector *v;

   /* Allocate the vector */
   v = (my_Vector *) malloc(sizeof(my_Vector));
   vec_create(mspace, &(v->uvals));
   vec_create(mspace, &(v->vvals));
   vec_create(mspace, &(v->wvals));
   vec_copy(mspace, (u->uvals), (v->uvals));
   vec_copy(mspace, (u->vvals), (v->vvals));
   vec_copy(mspace, (u->wvals), (v->wvals));

   *v_ptr = v;

   return 0;
}

/*------------------------------------*/

int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u->uvals);
   free(u->vvals);
   free(u->wvals);
   free(u);

   return 0;
}

/*------------------------------------*/

int
my_Sum(braid_App     app,
       double        alpha,
       braid_Vector  x,
       double        beta,
       braid_Vector  y)
{
   vec_axpy((app->mspace), alpha, (x->uvals), beta, (y->uvals));
   vec_axpy((app->mspace), alpha, (x->vvals), beta, (y->vvals));
   vec_axpy((app->mspace), alpha, (x->wvals), beta, (y->wvals));

   return 0;
}

/*------------------------------------*/

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   int i;
   double dot = 0.0;
   int mspace = (app->mspace);
   for (i = 0; i < mspace; i++)
   {
      dot += (u->uvals)[i]*(u->uvals)[i];
      dot += (u->vvals)[i]*(u->vvals)[i];
      dot += (u->wvals)[i]*(u->wvals)[i];
   }
   *norm_ptr = sqrt(dot);

   return 0;
}

/*------------------------------------*/

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int   done, index;
   int   mspace  = (app->mspace);

   /* Print solution to file if simulation is over */
   braid_AccessStatusGetDone(astatus, &done);
   if (done)
   {
      int  j;

      braid_AccessStatusGetTIndex(astatus, &index);

      fprintf((app->ufile), "%05d: ", index);
      fprintf((app->vfile), "%05d: ", index);
      fprintf((app->wfile), "%05d: ", index);
      for(j = 0; j < (mspace-1); j++)
      {
         fprintf((app->ufile), "% 1.14e, ", (u->uvals[j]));
         fprintf((app->vfile), "% 1.14e, ", (u->vvals[j]));
         fprintf((app->wfile), "% 1.14e, ", (u->wvals[j]));
      }
      fprintf((app->ufile), "% 1.14e\n", (u->uvals[j]));
      fprintf((app->vfile), "% 1.14e\n", (u->vvals[j]));
      fprintf((app->wfile), "% 1.14e\n", (u->wvals[j]));
   }

   return 0;
}

/*------------------------------------*/

int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
   *size_ptr = 3*(app->mspace)*sizeof(double);
   return 0;
}

/*------------------------------------*/

int
my_BufPack(braid_App           app,
           braid_Vector        u,
           void               *buffer,
           braid_BufferStatus  bstatus)
{
   double *dbuffer = buffer;

   vec_copy((app->mspace), (u->uvals), dbuffer);
   dbuffer += (app->mspace);
   vec_copy((app->mspace), (u->vvals), dbuffer);
   dbuffer += (app->mspace);
   vec_copy((app->mspace), (u->wvals), dbuffer);
   braid_BufferStatusSetSize(bstatus, 3*(app->mspace)*sizeof(double));

   return 0;
}

/*------------------------------------*/

int
my_BufUnpack(braid_App           app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus  bstatus)
{
   my_Vector *u = NULL;
   double    *dbuffer = buffer;

   /* Allocate memory */
   u = (my_Vector *) malloc(sizeof(my_Vector));
   vec_create((app->mspace), &(u->uvals));
   vec_create((app->mspace), &(u->vvals));
   vec_create((app->mspace), &(u->wvals));

   /* Unpack the buffer */
   vec_copy((app->mspace), dbuffer, (u->uvals));
   dbuffer += (app->mspace);
   vec_copy((app->mspace), dbuffer, (u->vvals));
   dbuffer += (app->mspace);
   vec_copy((app->mspace), dbuffer, (u->wvals));

   *u_ptr = u;
   return 0;
}

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int
main(int argc, char *argv[])
{
   braid_Core  core;
   my_App     *app;
         
   char        filename[255];
   double      tstart, tstop, dt, dx, start, end; 
   int         rank, ntime, mspace, arg_index;
   double      alpha, nu;
   int         max_levels, min_coarse, nrelax, nrelaxc, cfactor, maxiter, fmg;
   int         access_level, print_level, forward;
   double      tol;
   double      time;
   double     *u0, *scr;
   int         i;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   /* Define space domain. Space domain is between 0 and 1, mspace defines the number of steps */
   mspace = 8;
   ntime  = 256;

   /* Define some optimization parameters */
   alpha = .005;            /* parameter in the objective function */
   nu    = .03;             /* parameter in PDE (used 0.3 in RIPS) */

   /* Define some Braid parameters */
   max_levels     = 30;
   min_coarse     = 1;
   nrelax         = 1;
   nrelaxc        = 10;
   maxiter        = 50;
   cfactor        = 2;
   tol            = 1.0e-6;
   fmg            = 0;
   access_level   = 2;
   print_level    = 2;
   forward        = 0;

   /* Parse command line */
   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         printf("\n");
         printf(" Solves the advection-diffusion model problem \n\n");
         printf("  min  1/2 \\int_0^T\\int_0^1 (u(x,t)-ubar(x))^2 + alpha*v(x,t)^2  dxdt \n\n");
         printf("  s.t.  u_t + u_x - nu*u_xx = v(x,t) \n");
         printf("        u(0,t) = u(1,t) = 0 \n\n");
         printf("        u(x,0) = u0(x) \n");
         printf("  -tstop <tstop>          : Upper integration limit for time\n");
         printf("  -ntime <ntime>          : Num points in time\n");
         printf("  -mspace <mspace>        : Num points in space\n");
         printf("  -nu <nu>                : Constant Parameter in PDE  \n");
         printf("  -alpha <alpha>          : Constant Parameter in Objective Function  \n");
         printf("  -ml <max_levels>        : Max number of braid levels \n");
         printf("  -num  <nrelax>          : Num F-C relaxations\n");
         printf("  -nuc <nrelaxc>          : Num F-C relaxations on coarsest grid\n");
         printf("  -mi <maxiter>           : Max iterations \n");
         printf("  -cf <cfactor>           : Coarsening factor \n");
         printf("  -tol <tol>              : Stopping tolerance \n");
         printf("  -fmg                    : use FMG cycling\n");
         printf("  -access <access_level>  : Braid access level \n");
         printf("  -print <print_level>    : Braid print level \n");
         exit(1);
      }
      else if ( strcmp(argv[arg_index], "-ntime") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tstop") == 0 )
      {
         arg_index++;
         tstop = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mspace") == 0 )
      {
         arg_index++;
         mspace = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ml") == 0 )
      {
         arg_index++;
         max_levels = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu") == 0 )
      {
         arg_index++;
         nu = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-alpha") == 0 )
      {
         arg_index++;
         alpha = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-num") == 0 )
      {
         arg_index++;
         nrelax = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nuc") == 0 )
      {
         arg_index++;
         nrelaxc = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mi") == 0 )
      {
         arg_index++;
         maxiter = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-cf") == 0 )
      {
         arg_index++;
         cfactor = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tol") == 0 )
      {
         arg_index++;
         tol = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-fmg") == 0 )
      {
         arg_index++;
         fmg = 1;
      }
      else if ( strcmp(argv[arg_index], "-access") == 0 )
      {
         arg_index++;
         access_level = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-print") == 0 )
      {
         arg_index++;
         print_level = atoi(argv[arg_index++]);
      }
      else
      {
         printf("ABORTING: incorrect command line parameter %s\n", argv[arg_index]);
         return (0);
      }
   }

   /* Define the space step */
   dx = (double)1/(mspace+1);

   /* Define time domain and step */
   tstart = 0.0;             /* Beginning of time domain */
   tstop  = 1.0;             /* End of time domain*/
   dt     = (tstop-tstart)/ntime; 
    
   /* Set up initial condition */
   vec_create(mspace, &u0);
   for(i=0; i<mspace/2; i++)
   {
      u0[i] = 1;
//      u0[i] = 0; // RDF TMP
   }
   for(i=mspace/2; i<mspace; i++)
   {
      u0[i] = 0;
//      u0[i] = 1; // RDF TMP
   }

   /* Set up scratch space */
   vec_create(4*mspace, &scr);

   /* Run uncontrolled forward problem in serial */
   if (forward)
   {
      FILE    *file;
      double  *g, *u;
      int      j;

      file = fopen("trimgrit-burgers.out.forward", "w");

      vec_create(mspace, &g);
      vec_create(mspace, &u);

      for (i = 0; i < ntime; i++)
      {
         if (i == 0)
         {
            vec_copy(mspace, u0, g);
         }
         else
         {
            vec_copy(mspace, u, g);
         }

         apply_Phi(dt, dx, nu, mspace, g, u, scr);

         for(j = 0; j < (mspace-1); j++)
         {
            fprintf(file, "% 1.14e, ", u[j]);
         }
         fprintf(file, "\n");
      }

      fclose(file);
      return (0);
   }

   /* Set up the app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->myid     = rank;
   app->ntime    = ntime;
   app->mspace   = mspace;
   app->dx       = dx;
   app->nu       = nu;
   app->alpha    = alpha;
   app->w        = NULL;
   app->u0       = u0;
   app->scr      = scr;

   sprintf(filename, "%s.%03d", "trimgrit-burgers.out.u", (app->myid));
   (app->ufile) = fopen(filename, "w");
   sprintf(filename, "%s.%03d", "trimgrit-burgers.out.v", (app->myid));
   (app->vfile) = fopen(filename, "w");
   sprintf(filename, "%s.%03d", "trimgrit-burgers.out.w", (app->myid));
   (app->wfile) = fopen(filename, "w");

   /* Initialize XBraid */
   braid_InitTriMGRIT(MPI_COMM_WORLD, MPI_COMM_WORLD, dt, tstop, ntime-1, app,
                      my_TriResidual, my_TriSolve, my_Init, my_Clone, my_Free,
                      my_Sum, my_SpatialNorm, my_Access,
                      my_BufSize, my_BufPack, my_BufUnpack, &core);

   /* Set some XBraid(_Adjoint) parameters */
   braid_SetMaxLevels(core, max_levels);
   braid_SetMinCoarse(core, min_coarse);
   braid_SetNRelax(core, -1, nrelax);
   if (max_levels > 1)
   {
      braid_SetNRelax(core, max_levels-1, nrelaxc); /* nrelax on coarsest level */
   }
   braid_SetCFactor(core, -1, cfactor);
   braid_SetAccessLevel(core, access_level);
   braid_SetPrintLevel( core, print_level);       
   braid_SetMaxIter(core, maxiter);
   braid_SetAbsTol(core, tol);
   if (fmg)
   {
      braid_SetFMG(core);
   }

   /* Parallel-in-time TriMGRIT simulation */
   start=clock();
   braid_Drive(core);
   end=clock();

   /* Print runtime to file (for runtime comparisons)*/
   time = (double)(end-start)/CLOCKS_PER_SEC;
   printf("Total Run Time: %f s \n", time);

//   {
//      char    filename[255];
//      FILE   *file;
//      
//      //Note that this out file appends the number of time steps
//      sprintf(filename, "%s.%d", "trimgrit-burgers.time", ntime);
//
//      file = fopen(filename, "w");
//      fprintf(file, "%f", time);
//      fflush(file);
//      fclose(file);
//   }

//   /* RDF Testing */
//   {
//      double     cscale = (1/(dx*dt))*(2 + dt*dt/alpha);
//      double     lscale = (1/(dx*dt))*(1 + dt*dt/alpha);
//      my_Vector *e = (my_Vector *) malloc(sizeof(my_Vector));
//      my_Vector *z = (my_Vector *) malloc(sizeof(my_Vector));
//
//      (e->values) = (double*) malloc( mspace*sizeof(double) );
//      (z->values) = (double*) malloc( mspace*sizeof(double) );
//      for(i = 0; i < mspace; i++)
//      {
//         (e->values[i]) = 1.0;
//         (z->values[i]) = 0.0;
//      }
//      apply_TriResidual(app, z, z, NULL, e, 1, dt);
//
//      for(i = 0; i < mspace; i++)
//      {
//         (e->values[i]) = 1.0;
//         (z->values[i]) = 0.0;
//      }
//      apply_TriResidual(app, NULL, z, NULL, e, 1, dt);
//
//      for(i = 0; i < mspace; i++)
//      {
//         (e->values[i]) = 1.0;
//         (z->values[i]) = 0.0;
//      }
//      apply_TriResidual(app, z, NULL, NULL, e, 1, dt);
//   }

   vec_destroy(app->u0);
   vec_destroy(app->scr);
   fflush(app->ufile);
   fflush(app->vfile);
   fflush(app->wfile);
   fclose(app->ufile);
   fclose(app->vfile);
   fclose(app->wfile);
   free(app);

   braid_Destroy(core);
   MPI_Finalize();

   return (0);
}
