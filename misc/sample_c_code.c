/*BHEADER**********************************************************************
 * Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of WARP.  See file COPYRIGHT for details.
 *
 * WARP is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 ***********************************************************************EHEADER*/


/******************************************************************************
 *
 * This is an example of a C code illustrating the indentation used
 * for Warp.  This code does not illustrate issues related to
 * error handling, efficiency of implementation or naming conventions.
 * 
 * The most important item here is consistent indentation of the following
 * structures:
 *    - for loops
 *    - if statements
 *
 * We use the ellemtel style with tab and shift widths of 3.  See 
 * sample.emacs and sample.vimrc for settings you should use. 
 *
 * This example also shows a sample comment block for a function and it's 
 * parameters that will work well with doxygen. 
 *
 *****************************************************************************/

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
   double           *a  = MatrixData(A);        /* element values for matrix A */
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
