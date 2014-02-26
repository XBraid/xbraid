#include <stdlib.h>
#include "c_array.h"
#include "advect_data.h"
void
dwdx( grid_fcn *w, grid_fcn *wx, advection_setup *kd_ )
{
   const double d6a= 0.75, d6b=-0.15, d6c=1.0/60;
   
   int i, k;
/* w: input grid function with ghost points */
/* wx: output grid function with ghost points. Only 1 <= i <= n are assigned */
   
/* bop: coefficients for standard SBP */
/* bope, gh: coefficients for SBP with ghost point */
/* bope(i,k) = bop(i,k) for i>=2 */
/* nb: 1st dim of bop & bope arrays */
/* wb: 2nd dim of -- " " --         */
#define bcnr(i)  compute_index_1d(kd_->bcnr_, i)
#define bop(i,j)  compute_index_2d(kd_->bop_, i, j)
#define bope(i,j) compute_index_2d(kd_->bope_, i, j)
   double ih, i2h, du, du1;
   double h = w->h;
   int n = w->n;
   const int off=2;
   double *wpad;
      
   i2h = 1.0/(2*h);
   ih = 1.0/(h);

   if (bcnr(1) == Periodic)
   {
/* allocate tmp space */
      wpad=(double*) malloc((n+6)*sizeof(double));
/* copy data (assume w->sol[1] = w->sol[n] */
      for (i=1; i<=n; i++)
      {
         wpad[i+off] = w->sol[i];
      }
      wpad[off-0] = w->sol[n-1];
      wpad[off-1] = w->sol[n-2];
      wpad[off-2] = w->sol[n-3];

      wpad[n+1+off] = w->sol[2];
      wpad[n+2+off] = w->sol[3];
      wpad[n+3+off] = w->sol[4];
/* centered formula for all points */
      for (i=1; i<=n; i++)
      {
         du  = d6c*(wpad[i+3+off]-wpad[i-3+off]) + d6b*(wpad[i+2+off]-wpad[i-2+off]) + d6a*(wpad[i+1+off]-wpad[i-1+off]);
         wx->sol[i] = ih*(du);
      }
/* free tmp space */
      free(wpad);
   }
   else
   { /* non-periodic case... */
      
/*! on the boundary */
      i=1;
   
      du  = 0;
      du1 = 0;
   
      for (k=1; k<=kd_->wb; k++)
      {
         du  = du  + bope(i,k) * w->sol[k];
         du1 = du1 + bop(i,k)  * w->sol[k];
      }
   
      du = du + kd_->gh*w->sol[0];
      wx->sol[i] = i2h*(du + du1 );

/* near the left bndry */
      for (i=2; i<=kd_->nb; i++)
      {
         du  = 0;
         for (k=1; k<=kd_->wb; k++)
         {
            du  = du  + bope(i,k) * w->sol[k];
         }
         wx->sol[i] = ih*(du );
      }

/* centered formula in the interior */
      for (i=kd_->nb+1; i<=n-kd_->nb; i++)
      {
         du  = d6c*(w->sol[i+3]-w->sol[i-3]) + d6b*(w->sol[i+2]-w->sol[i-2]) + d6a*(w->sol[i+1]-w->sol[i-1]);
         wx->sol[i] = ih*(du );
      }

/* near the right bndry */
      for (i=n-kd_->nb+1; i<=n-1; i++)
      {
         du   = 0;
         for ( k=n-kd_->wb+1; k<=n; k++)
         {
            du = du  - bope(n-i+1,n-k+1) * w->sol[k];
         }
         wx->sol[i] = ih*(du );
      }

/*! on the right bndry */
      i=n;
      du = 0;
      du1 = 0;
   
      for (k=n-kd_->wb+1; k<=n; k++)
      {
         du  = du  - bope(n-i+1,n-k+1) * w->sol[k];
         du1 = du1 - bop(n-i+1,n-k+1)  * w->sol[k];
      }
      du = du - kd_->gh * w->sol[n+1];
      wx->sol[i] = i2h* (du + du1 );
   }
   
#undef bop
#undef bope
}

/* 2nd derivative */
void
d2wdx2( grid_fcn *w, grid_fcn *wxx, advection_setup *kd_ )
/*d2wdx2( int n, double *w, double *wxx, double h, advection_setup *kd_ )*/
/* change the double pointers to *grid_fcn */
{
   int i, k;

/* w: input grid function with ghost points */
/* wx: output grid function with ghost points. Only 1 <= i <= n are assigned */
   
/* bop2: boundary coefficients for standard SBP */
/* gh2:  ghost point coefficient */
/* iop2(i,k): interior coefficients for i>= nb2 + 1 */
/* nb2: 1st dim of bop2 array */
/* wb2: 2nd dim of -- " " --         */
#define bcnr(i)  compute_index_1d(kd_->bcnr_, i)
#define bop2(i,j)  compute_index_2d(kd_->bop2_, i, j)
#define iop2(i)  compute_index_1d(kd_->iop2_, i)
   double ih2, du, d60, d61, d62, d63;
   int n=w->n;
   const int off=2;
   double *wpad;
      
   ih2 = 1.0/(w->h*w->h);

/* coefficients for centered formula */
   d60 = iop2(4);
   d61 = iop2(3);
   d62 = iop2(2);
   d63 = iop2(1);

   if (bcnr(1) == Periodic)
   {
/* allocate tmp space */
      wpad=(double*) malloc((n+6)*sizeof(double));
/* copy data (assume w->sol[1] = w->sol[n] */
      for (i=1; i<=n; i++)
      {
         wpad[i+off] = w->sol[i];
      }
      wpad[off-0] = w->sol[n-1];
      wpad[off-1] = w->sol[n-2];
      wpad[off-2] = w->sol[n-3];

      wpad[n+1+off] = w->sol[2];
      wpad[n+2+off] = w->sol[3];
      wpad[n+3+off] = w->sol[4];
/* centered formula for all points */
      for (i=1; i<=n; i++)
      {
         du = d60*wpad[i+off] + d61*(wpad[i-1+off] + wpad[i+1+off]) + d62*(wpad[i-2+off] + wpad[i+2+off])
            + d63*(wpad[i-3+off] + wpad[i+3+off]);
         wxx->sol[i] = ih2*(du );
      }
/* free tmp space */
      free(wpad);
   }
   else
   { /* non-periodic case... */
/*! on the boundary */
      i=1;
   
      du  = 0;
      for (k=1; k<=kd_->wb2; k++)
      {
         du  = du  + bop2(i,k) * w->sol[k];
      }
   
      du = du + kd_->gh2*w->sol[0];
      wxx->sol[i] = ih2*(du);

/* near the left bndry */
      for (i=2; i<=kd_->nb2; i++)
      {
         du  = 0;
         for (k=1; k<=kd_->wb2; k++)
         {
            du  = du  + bop2(i,k) * w->sol[k];
         }
         wxx->sol[i] = ih2*(du );
      }

/* centered formula in the interior */
      for (i=kd_->nb2 + 1; i<=n - kd_->nb2; i++)
      {
         du = iop2(4)*w->sol[i] + iop2(1)*(w->sol[i-3] + w->sol[i+3])  + iop2(2)*(w->sol[i-2] + w->sol[i+2])  + iop2(3)*(w->sol[i-1] + w->sol[i+1]);
         wxx->sol[i] = ih2*(du );
      }

/* near the right bndry */
      for (i=n - kd_->nb2 + 1; i<=n-1; i++)
      {
         du   = 0;
         for (k=n - kd_->wb2 + 1; k<= n; k++)
         {
            du = du  + bop2(n-i+1,n-k+1) * w->sol[k];
         }
         wxx->sol[i] = ih2*(du );
      }

/*! on the right bndry */
      i=n;
      du = 0;
   
      for (k=n - kd_->wb2+1; k<=n; k++)
      {
         du  = du  + bop2(n-i+1, n-k+1) * w->sol[k];
      }
      du = du + kd_->gh2 * w->sol[n+1];
      wxx->sol[i] = ih2* (du );
   } /* end non-periodic case */
#undef bop
#undef bope
}

