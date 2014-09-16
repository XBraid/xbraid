/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder schroder2@llnl.gov, Rob Falgout falgout2@llnl.gov,
 * Tzanio Kolev kolev1@llnl.gov, Ulrike Yang yang11@llnl.gov, 
 * Veselin Dobrev dobrev1@llnl.gov, et al. 
 * LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. Email schroder2@llnl.gov on how to download. 
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
 *
 * This is an example of a C code illustrating the indentation used
 * for Braid.  This code does not illustrate issues related to
 * error handling, efficiency of implementation or naming conventions.
 * 
 * The most important item here is consistent indentation of the following
 * structures:
 *    - for loops
 *    - if statements
 *
 * We use the ellemtel style with tab and shift widths of 3.  See 
 * sample.emacs and sample.vimrc for settings you can use. 
 *
 * This example also shows a sample comment block for a function and 
 * parameters that works well with doxygen. 
 *
 */ 

#include "headers.h"

/**
 * Matvec - matrix-vector function.
 * 
 * Calculates y = alpha * A * x + beta * y
 * where A is a matrix stored in compressed sparse row format, x and y
 * are n vectors, and alpha and beta are scalars.
 * 
 * Output:
 * - *y* will point to the vector
 **/
void
Matvec( double  alpha,  /**< Describe parameters like this for doxygen*/
        Matrix *A,      /**< Matrix for mat-vec operation */
        Vector *x,      /**< Vector to be multiplied */
        double  beta,   /**< scalar */
        Vector *y       /**< Vector for addition */    
        )
{
   double           *a  = MatrixData(A);  /* element values for matrix A */
   HYPRE_Int        *ia = MatrixIA(A);    /* pointer to start of each row */
   HYPRE_Int        *ja = MatrixJA(A);    /* column values for matrix elements */
   HYPRE_Int         n  = MatrixSize(A);  /* size of matrix */

   double     *xp = VectorData(x);
   double     *yp = VectorData(y);

   double      temp;

   HYPRE_Int         i, j, jj;

   /*
    * Do (alpha == 0.0) computation 
    */

   if (alpha == 0.0)
   {
      for (i = 0; i < n; i++)
         yp[i] *= beta;

      return;
   }

   /* 
    * y = (beta/alpha)*y
    */
   
   temp = beta / alpha;
   
   if (temp != 1.0)
   {
      if (temp == 0.0)
      {
         for (i = 0; i < n; i++)
            yp[i] = 0.0;
      }
      else
      {
         for (i = 0; i < n; i++)
            yp[i] *= temp;
      }
   }

   /*
    * y += A*x
    */

   for (i = 0; i < n; i++)
   {
      for (jj = ia[i]-1; jj < ia[i+1]-1; jj++)
      {
         j = ja[jj]-1;
         yp[i] += a[jj] * xp[j];
      }
   }

   /*
    * y = alpha*y
    */

   if (alpha != 1.0)
   {
      for (i = 0; i < n; i++)
         yp[i] *= alpha;
   }
}
