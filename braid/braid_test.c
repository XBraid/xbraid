/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory.
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
 

/** \file braid_test.c
 * \brief Define XBraid test routines.
 */

#include "_braid.h"
#include "util.h"

/*--------------------------------------------------------------------------
 * Some simple tests on the myinit, myaccess and myfree routines
 *--------------------------------------------------------------------------*/
braid_Int
braid_TestInitAccess( braid_App           app, 
                     MPI_Comm            comm_x,
                     FILE               *fp, 
                     braid_Real          t,
                     braid_PtFcnInit     myinit,
                     braid_PtFcnAccess   myaccess,
                     braid_PtFcnFree     myfree)
{
   
   braid_Vector          u ;
   braid_Status          status = _braid_CTAlloc(_braid_Status, 1);
   braid_AccessStatus    astatus = (braid_AccessStatus)status;
   braid_Int             myid_x = 0;
   
   _braid_AccessStatusInit(t, 0, 0.0, 0, 0, 0, 0, 0, 1, -1, NULL, astatus);
   MPI_Comm_rank( comm_x, &myid_x );

   /* Print intro */
   _braid_ParFprintfFlush(fp, myid_x, "\nStarting braid_TestInitAccess\n\n");

   /* Test */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInitAccess:   Starting Test 1\n");
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInitAccess:   u = init(t=%1.2e)\n", t);
   myinit(app, t, &u);
   
   if(myaccess != NULL)
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInitAccess:   access(u) \n");
      myaccess(app, u, astatus);

      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInitAccess:   check output: wrote u for initial condition at t=%1.2e. \n\n",t);
   }

   /* Free variables */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInitAccess:   free(u) \n");
   myfree(app, u);
   _braid_StatusDestroy(status);
   
   _braid_ParFprintfFlush(fp, myid_x, "Finished braid_TestInitAccess\n");

   return 0;
}

braid_Int
braid_TestClone( braid_App        app,  
              MPI_Comm            comm_x,
              FILE               *fp, 
              braid_Real          t,
              braid_PtFcnInit     myinit,
              braid_PtFcnAccess   myaccess,
              braid_PtFcnFree     myfree,
              braid_PtFcnClone    clone)
{
   
   braid_Vector        u, v;
   braid_Status        status = _braid_CTAlloc(_braid_Status, 1);
   braid_AccessStatus  astatus = (braid_AccessStatus)status;
   braid_Int           myid_x;
   
   _braid_AccessStatusInit(t, 0, 0.0, 0, 0, 0, 0, 0, 1, -1, NULL, astatus);
   MPI_Comm_rank( comm_x, &myid_x );

   /* Print intro */
   _braid_ParFprintfFlush(fp, myid_x, "\nStarting braid_TestClone\n\n");

   /* Test 1 */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestClone:   Starting Test 1\n");
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestClone:   u = init(t=%1.2e)\n", t);
   myinit(app, t, &u);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestClone:   v = clone(u)\n");
   clone(app, u, &v);
   
   if(myaccess != NULL)
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestClone:   access(u)\n");
      myaccess(app, u, astatus);

      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestClone:   access(v)\n");
      myaccess(app, v, astatus);
      
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestClone:   check output:  wrote u and v for initial condition at t=%1.2e.\n\n", t);

   }

   /* Free variables */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestClone:   free(u)\n");
   myfree(app, u);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestClone:   free(v)\n");
   myfree(app, v);

   _braid_StatusDestroy(status);
   
   _braid_ParFprintfFlush(fp, myid_x, "Finished braid_TestClone\n");
   
   return 0;
}



braid_Int
braid_TestSum( braid_App        app, 
            MPI_Comm            comm_x,
            FILE               *fp,
            braid_Real          t,
            braid_PtFcnInit     myinit,
            braid_PtFcnAccess   myaccess,
            braid_PtFcnFree     myfree,
            braid_PtFcnClone    clone,
            braid_PtFcnSum      sum )
{
   
   braid_Vector        u, v;
   braid_Status        status  = _braid_CTAlloc(_braid_Status, 1);
   braid_AccessStatus  astatus = (braid_AccessStatus)status;
   braid_Int           myid_x;
   
   _braid_AccessStatusInit(t, 0, 0.0, 0, 0, 0, 0, 0, 1, -1, NULL, astatus);
   MPI_Comm_rank( comm_x, &myid_x );

   /* Print intro */
   _braid_ParFprintfFlush(fp, myid_x, "\nStarting braid_TestSum\n\n");
   
   /* Test 1 */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSum:   Starting Test 1\n");
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSum:   u = init(t=%1.2e)\n", t);
   myinit(app, t, &u);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSum:   v = clone(u)\n");
   clone(app, u, &v);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSum:   v = u - v\n");
   sum(app, 1.0, u, -1.0, v); 

   if(myaccess != NULL)
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSum:   access(v)\n");
      myaccess(app, v, astatus);
      
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSum:   check output:  v should equal the zero vector\n\n");
   }

   /* Test 2 */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSum:   Starting Test 2\n");
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSum:   v = 2*u + v\n");
   sum(app, 2.0, u, 1.0, v); 

   if(myaccess != NULL)
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSum:   access(v)\n");
      myaccess(app, v, astatus);
      
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSum:   access(u)\n");
      myaccess(app, u, astatus);
   }

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSum:   check output:  v should equal 2*u \n\n");

   /* Free variables */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSum:   free(u)\n");
   myfree(app, u);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSum:   free(v)\n");
   myfree(app, v);

   _braid_StatusDestroy(status);
   
   _braid_ParFprintfFlush(fp, myid_x, "Finished braid_TestSum\n");

   return 0;
}

braid_Int
braid_TestSpatialNorm( braid_App              app, 
                       MPI_Comm               comm_x,
                       FILE                  *fp, 
                       braid_Real             t,
                       braid_PtFcnInit        myinit,
                       braid_PtFcnFree        myfree,
                       braid_PtFcnClone       clone,
                       braid_PtFcnSum         sum,  
                       braid_PtFcnSpatialNorm spatialnorm) 
{   
   braid_Vector  u, v, w;
   braid_Real    result1, result2;
   braid_Int     myid_x, correct;
   double     wiggle = 1e-12;
   
   MPI_Comm_rank( comm_x, &myid_x );

   /* Initialize the flags */
   correct = 1;
   braid_Int zero_flag = 0;

   /* Print intro */
   _braid_ParFprintfFlush(fp, myid_x, "\nStarting braid_TestSpatialNorm\n\n");
   
   /* Test 1 */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   Starting Test 1\n");
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   u = init(t=%1.2e)\n", t);
   myinit(app, t, &u);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   spatialnorm(u) \n");
   spatialnorm(app, u, &result1);
   if( fabs(result1) == 0.0)
   {
      zero_flag = 1;
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   Warning:  spatialnorm(u) = 0.0\n"); 
   }
   else if( _braid_isnan(result1) )
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   Warning:  spatialnorm(u) = nan\n"); 
      correct = 0;
   }

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   v = clone(u)\n");
   clone(app, u, &v);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   v = u - v \n");
   sum(app, 1.0, u, -1.0, v); 

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   spatialnorm(v) \n");
   spatialnorm(app, v, &result1);
   if( (fabs(result1) > wiggle) || _braid_isnan(result1) )
   {
      correct = 0;
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   Test 1 Failed\n");
   }
   else
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   Test 1 Passed\n");
   }
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   actual output:    spatialnorm(v) = %1.2e  \n", result1);
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   expected output:  spatialnorm(v) = 0.0 \n\n");
   

   /* Test 2 */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   Starting Test 2\n");
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   w = clone(u)\n");
   clone(app, u, &w);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   w = u + w \n");
   sum(app, 1.0, u, 1.0, w); 

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   spatialnorm(u)\n");
   spatialnorm(app, u, &result1);
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   spatialnorm(w)\n");
   spatialnorm(app, w, &result2);
   if( (fabs(result2/result1 - 2.0) > wiggle) || _braid_isnan(result2/result1) )
   {
      correct = 0;
      if(zero_flag)
      {
         _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   Test 2 Failed, Likely due to u = 0\n");
      }
      else
      {
         _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   Test 2 Failed\n");
      }
   }
   else
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   Test 2 Passed\n");
   }
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   actual output:    spatialnorm(w) / spatialnorm(u) = %1.2e / %1.2e = %1.2e \n",
         result2, result1, result2/result1);
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   expected output:  spatialnorm(w) / spatialnorm(u) = 2.0 \n\n");

   /* Test 3 */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   Starting Test 3\n");
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   free(w)\n");
   myfree(app, w);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   w = clone(u)\n");
   clone(app, u, &w);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   w = 0.0*u + 0.5*w \n");
   sum(app, 0.0, u, 0.5, w); 

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   spatialnorm(u)\n");
   spatialnorm(app, u, &result1);
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   spatialnorm(w)\n");
   spatialnorm(app, w, &result2);
   /* Check Result */
   if( (fabs(result2/result1 - 0.5) > wiggle) || _braid_isnan(result2/result1) )

   {
      correct = 0;
      if(zero_flag)
      {
         _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   Test 3 Failed, Likely due to u = 0\n");
      }
      else
      {
         _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   Test 3 Failed\n");
      }
   }
   else
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   Test 3 Passed\n");
   }
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   actual output:    spatialnorm(w) / spatialnorm(u) = %1.2e / %1.2e = %1.2e \n",
         result2, result1, result2/result1);
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   expected output:  spatialnorm(w) / spatialnorm(u) = 0.5 \n\n");

   /* Free variables */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   free(u)\n");
   myfree(app, u);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   free(v)\n");
   myfree(app, v);
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestSpatialNorm:   free(w)\n");
   myfree(app, w);

   if(correct == 1) 
      _braid_ParFprintfFlush(fp, myid_x, "Finished braid_TestSpatialNorm: all tests passed successfully\n");
   else
   {
      if(zero_flag)
      {
         _braid_ParFprintfFlush(fp, myid_x, "Finished braid_TestSpatialNorm: some tests failed, possibly due to u = 0\n");
      }
      else
      {
         _braid_ParFprintfFlush(fp, myid_x, "Finished braid_TestSpatialNorm: some tests failed\n");
      }
   }

   return correct;
}

braid_Int
braid_TestInnerProd( braid_App              app, 
                     MPI_Comm               comm_x,
                     FILE                  *fp,
                     braid_Real             t1,
                     braid_Real             t2,
                     braid_PtFcnInit        myinit,
                     braid_PtFcnFree        myfree,
                     braid_PtFcnSum         sum,
                     braid_PtFcnInnerProd   inner_prod)
{
   braid_Vector u, v;
   braid_Real   result1, result2;
   braid_Int    myid_x, zero_flag1, zero_flag2;

   double wiggle = 1e-12;

   MPI_Comm_rank(comm_x, &myid_x);

   /* Initialize the flags */
   braid_Int correct = 1;
   zero_flag1 = 0;
   zero_flag2 = 0;

   /* Print intro */
   _braid_ParFprintfFlush(fp, myid_x, "\nStarting braid_TestInnerProd\n\n");
   
   /* Test 1 */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   Starting Test 1\n");
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   u = init(t=%1.2e)\n", t1);
   myinit(app, t1, &u);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   inner_prod(u, u)\n");
   inner_prod(app, u, u, &result1);

   if( fabs(result1) == 0.0)
   {
      zero_flag1 = 1;
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   Warning:  inner_prod(u, u) = 0.0\n"); 
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   Test 1 Inconclusive (make sure u is not the zero vector)\n"); 
   }
   else if( _braid_isnan(result1) )
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   Warning:  inner_prod(u, u) = nan\n");
      correct = 0;
   }

   if (!zero_flag1)
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   u = u / sqrt(inner_prod(u, u))\n");
      sum(app, 0., u, 1/sqrt(result1), u);

      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   inner_prod(u, u)\n");
      inner_prod(app, u, u, &result2);
      
      if( (fabs(1. - result2) > wiggle) || _braid_isnan(result2) )
      {
         correct = 0;
         _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   Test 1 Failed\n");
      }
      else
      {
         _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   Test 1 Passed\n");
      }
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   actual output:    inner_prod(u, u) = %1.2e  \n", result2);
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   expected output:  inner_prod(u, u) = 1.0 \n\n");
   }

   /* Test 2 */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   Starting Test 2\n");
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   v = init(t=%1.2e)\n", t2);
   myinit(app, t2, &v);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   inner_prod(u, v)\n");
   inner_prod(app, u, v, &result1);

   if( fabs(result1) == 0.0 )
   {
      zero_flag2 = 1;
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   Warning:  inner_prod(u, v) = 0.0\n"); 
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   Test 2 Inconclusive (make sure v is not the zero vector and is not orthogonal to u)\n"); 
   }
   else if( _braid_isnan(result1) )
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   Warning:  inner_prod(u, v) = nan\n");
      correct = 0;
   }

   if ( !zero_flag1 && !zero_flag2 )
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   v = v - inner_prod(u, v)*u\n");
      sum(app, -result1, u, 1., v);

      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   inner_prod(u, v)\n");
      inner_prod(app, u, v, &result2);
      
      if( (fabs(result2) > wiggle) || _braid_isnan(result2) )
      {
         correct = 0;
         _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   Test 2 Failed\n");
      }
      else
      {
         _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   Test 2 Passed\n");
      }
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   actual output:    inner_prod(u, v) = %1.2e  \n", result2);
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestInnerProd:   expected output:  inner_prod(u, v) = 0. \n\n");
   }

   return correct;
} 

braid_Int
braid_TestBuf( braid_App              app,
               MPI_Comm               comm_x,
               FILE                   *fp, 
               braid_Real             t,
               braid_PtFcnInit        myinit,
               braid_PtFcnFree        myfree,
               braid_PtFcnSum         sum,  
               braid_PtFcnSpatialNorm spatialnorm, 
               braid_PtFcnBufSize     bufsize,
               braid_PtFcnBufPack     bufpack,
               braid_PtFcnBufUnpack   bufunpack)
{   
   braid_Vector  u, v;
   braid_Real    result1;
   braid_Int     myid_x, size, correct;
   void      *buffer;
   double     wiggle = 1e-12;
   
   MPI_Comm_rank( comm_x, &myid_x );
   
   braid_Status            status = _braid_CTAlloc(_braid_Status, 1);
   braid_BufferStatus      bstatus = (braid_BufferStatus)status;
   _braid_BufferStatusInit( 0, 0, bstatus );
   /* Initialize the correct flag */
   correct = 1;

   /* Print intro */
   _braid_ParFprintfFlush(fp, myid_x, "\nStarting braid_TestBuf\n\n");
   
   /* Test 1 */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   Starting Test 1\n");
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   u = init(t=%1.2e)\n", t);
   myinit(app, t, &u);
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   spatialnorm(u) \n");
   spatialnorm(app, u, &result1);
   if( fabs(result1) == 0.0)
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   Warning:  spatialnorm(u) = 0.0\n"); 
   }
   else if( _braid_isnan(result1) )
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   Warning:  spatialnorm(u) = nan\n"); 
      correct = 0;
   }

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   size = bufsize()\n");
   bufsize(app, &size, bstatus);
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   buffer = malloc(size)\n");
   buffer = malloc(size);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   buffer = bufpack(u, buffer))\n");
   
   _braid_StatusElt( bstatus, size_buffer ) = size;
   bufpack(app, u, buffer, bstatus);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   v = bufunpack(buffer)\n");
   bufunpack(app, buffer, &v, bstatus);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   v = u - v \n");
   sum(app, 1.0, u, -1.0, v); 

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   spatialnorm(v) \n");
   spatialnorm(app, v, &result1);
   if( (fabs(result1) > wiggle) || _braid_isnan(result1) )

   {
      correct = 0;
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   Test 1 Failed\n");
   }
   else
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   Test 1 Passed\n");
   }
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   actual output:    spatialnorm(v) = %1.2e  \n", result1);
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   expected output:  spatialnorm(v) = 0.0 \n\n");
   
   /* Free variables */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   free(u)\n");
   myfree(app, u);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestBuf:   free(v)\n");
   myfree(app, v);
   
   if(correct == 1) 
      _braid_ParFprintfFlush(fp, myid_x, "Finished braid_TestBuf: all tests passed successfully\n");
   else      
      _braid_ParFprintfFlush(fp, myid_x, "Finished braid_TestBuf: some tests failed\n");

   _braid_StatusDestroy(status);
   return correct;
}

braid_Int
braid_TestCoarsenRefine( braid_App           app,
                      MPI_Comm               comm_x,
                      FILE                   *fp, 
                      braid_Real             t,
                      braid_Real             fdt,
                      braid_Real             cdt,
                      braid_PtFcnInit        myinit,
                      braid_PtFcnAccess      myaccess,
                      braid_PtFcnFree        myfree,
                      braid_PtFcnClone       clone,
                      braid_PtFcnSum         sum,
                      braid_PtFcnSpatialNorm spatialnorm, 
                      braid_PtFcnSCoarsen    coarsen,
                      braid_PtFcnSRefine     refine)
 {   
   braid_Vector            u, v, w, uc, vc, wc;
   braid_Real              result1;
   braid_Int               myid_x, level, correct;
   braid_Status            status  = _braid_CTAlloc(_braid_Status, 1);
   braid_AccessStatus      astatus = (braid_AccessStatus)status;
   braid_CoarsenRefStatus  cstatus = (braid_CoarsenRefStatus)status;
   
   _braid_CoarsenRefStatusInit(t, t-fdt, t+fdt, t-cdt, t+cdt, 0, 0, 0, 0, cstatus);
   MPI_Comm_rank( comm_x, &myid_x );

   /* Initialize the correct flag */
   correct = 1;

   /* Print intro */
   _braid_ParFprintfFlush(fp, myid_x, "\nStarting braid_TestCoarsenRefine\n\n");
   
   /* Test 1 */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   Starting Test 1\n");
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   u = init(t=%1.2e)\n", t);
   myinit(app, t, &u);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   spatialnorm(u) \n");
   spatialnorm(app, u, &result1);
   if( fabs(result1) == 0.0)
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   Warning:  spatialnorm(u) = 0.0\n"); 
   }
   else if( _braid_isnan(result1) )
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   Warning:  spatialnorm(u) = nan\n"); 
      correct = 0;
   }

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   uc = coarsen(u)\n");
   coarsen(app, u, &uc, cstatus); 

   if(myaccess != NULL)
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   access(uc) \n");
      level = 1;
      _braid_AccessStatusInit(t, 0, 0.0, 0, level, 0, 0, 0, 1, -1, NULL, astatus);
      myaccess(app, uc, astatus);

      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   access(u) \n");
      level = 0;
      _braid_AccessStatusInit(t, 0, 0.0, 0, level, 0, 0, 0, 1, -1, NULL, astatus);
      myaccess(app, u, astatus);
   }

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   actual output:   wrote u and spatially coarsened u \n\n");

   /* Test 2 */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   Starting Test 2\n");
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   v = clone(u)\n");
   clone(app, u, &v);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   vc = coarsen(v)\n");
   coarsen(app, v, &vc, cstatus); 
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   wc = clone(vc)\n");
   clone(app, vc, &wc);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   wc = uc - wc \n");
   sum(app, 1.0, uc, -1.0, wc); 

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   spatialnorm(wc)\n");
   spatialnorm(app, wc, &result1);
   
   /* We expect exact equality between uc and vc */
   if( (fabs(result1) != 0.0) || _braid_isnan(result1) )
   {
      correct = 0;
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   Test 2 Failed\n");
   }
   else
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   Test 2 Passed\n");
   }
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   actual output:    spatialnorm(wc) = %1.2e \n", result1);
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   expected output:  spatialnorm(wc) = 0.0 \n\n");

   /* Test 3 */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   Starting Test 3\n");
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   w = clone(u)\n");
   clone(app, u, &w);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   free(u)\n");
   myfree(app, u);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   free(v)\n");
   myfree(app, v);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   v = refine(vc)\n");
   refine(app, vc, &v, cstatus); 

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   u = refine(uc)\n");
   refine(app, uc, &u, cstatus); 

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   v = u - v \n");
   sum(app, 1.0, u, -1.0, v); 

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   spatialnorm(v)\n");
   spatialnorm(app, v, &result1);
   
   /* We expect exact equality between u and v */
   if( (fabs(result1) != 0.0) || _braid_isnan(result1) )
   {
      correct = 0;
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   Test 3 Failed\n");
   }
   else
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   Test 3 Passed\n");
   }
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   actual output:    spatialnorm(v) = %1.2e \n", result1);
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   expected output:  spatialnorm(v) = 0.0 \n\n");

   /* Test 4 */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   Starting Test 4\n");
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   w = w - u \n");
   sum(app, 1.0, u, -1.0, w); 

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   spatialnorm(w)\n");
   spatialnorm(app, w, &result1);
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   actual output:    spatialnorm(w) = %1.2e \n", result1);
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   %s%s%s", 
                   "expected output:  For simple interpolation formulas\n",
                   "                             (e.g., bilinear) and a known function\n",
                   "                             (e.g., constant), spatialnorm(w) should = 0\n\n");

   /* Free variables */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   free(u)\n");
   myfree(app, u);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   free(v)\n");
   myfree(app, v);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   free(w)\n");
   myfree(app, w);
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   free(uc)\n");
   myfree(app, uc);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   free(vc)\n");
   myfree(app, vc);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestCoarsenRefine:   free(wc)\n");
   myfree(app, wc);

   _braid_StatusDestroy(status);

   if(correct == 1) 
      _braid_ParFprintfFlush(fp, myid_x, "Finished braid_TestCoarsenRefine: all tests passed successfully\n");
   else      
      _braid_ParFprintfFlush(fp, myid_x, "Finished braid_TestCoarsenRefine: some tests failed\n");
   
   return correct;
}

braid_Int
braid_TestResidual( braid_App              app,
                    MPI_Comm               comm_x,
                    FILE                   *fp, 
                    braid_Real             t,
                    braid_Real             dt,
                    braid_PtFcnInit        myinit,
                    braid_PtFcnAccess      myaccess,
                    braid_PtFcnFree        myfree,
                    braid_PtFcnClone       clone,
                    braid_PtFcnSum         sum,
                    braid_PtFcnSpatialNorm spatialnorm, 
                    braid_PtFcnResidual    residual,
                    braid_PtFcnStep        step)
 {   
   braid_Vector            u, unext, ustop, fstop;
   braid_Real              result1;
   braid_Int               myid_x, result_int;
   braid_Status            status  = _braid_CTAlloc(_braid_Status, 1);
   braid_AccessStatus      astatus = (braid_AccessStatus) status;
   braid_StepStatus        sstatus = (braid_StepStatus) status;
 

   /* 
    * Next, we must initialize status so that the user may call
    * braid_StepStatusSetRFactor() from inside of step(), which many users will
    * do.  Calling braid_StepStatusSetRFactor() requires that status (which is
    * really a braid_Core) to have allocated two pieces of data, 
    * (1) core->rfactors and (2) core->grids[0]->ilower
    */
   braid_Core              core    = (braid_Core) status;
   _braid_Grid          **grids    = _braid_CoreElt(core, grids);
   braid_Int              *rfactors;
   _braid_Grid            *fine_grid;
   /* 1) Initialize array of grids and fine-grid, then set ilower */ 
   grids    = _braid_TReAlloc(grids, _braid_Grid *, 1);
   fine_grid = _braid_CTAlloc(_braid_Grid, 1);
   _braid_GridElt(fine_grid, ilower) = 0;
   grids[0] = fine_grid;
   _braid_CoreElt(core, grids) = grids;
   /* 2) Initialize rfactors */
   rfactors = _braid_CTAlloc(braid_Int, 4); 
   _braid_CoreElt(core, rfactors) = rfactors;

   _braid_StepStatusInit(t, t+dt, 0, 1e-16, 0, 0, 0, 2, 0, NULL, sstatus);
   _braid_AccessStatusInit(t, 0, 0.0, 0, 0, 0, 2, 0, 1, -1, NULL, astatus);

   MPI_Comm_rank( comm_x, &myid_x );

   /* Print intro */
   _braid_ParFprintfFlush(fp, myid_x, "\nStarting braid_TestResidual\n\n");
   
   /* Test 1 */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   Starting Test 1\n");
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   u = init(t=%1.2e)\n", t);
   myinit(app, t, &u);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   spatialnorm(u) \n");
   spatialnorm(app, u, &result1);
   if( fabs(result1) == 0.0)
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   Warning:  spatialnorm(u) = 0.0\n"); 
   }
   else if( _braid_isnan(result1) )
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   Warning:  spatialnorm(u) = nan\n"); 
   }
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   unext = clone(u)\n");
   clone(app, u, &unext);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   ustop = clone(u)\n");
   clone(app, u, &ustop);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   fstop = clone(u)\n");
   clone(app, u, &fstop);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   fstop = ustop + fstop \n");
   sum(app, 1.0, ustop, 1.0, fstop); 

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   unext = step(ustop, fstop, unext) \n");
   step(app, ustop, fstop, unext, sstatus); 
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   ustop = clone(u)\n");
   clone(app, u, &ustop);
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   fstop = clone(u)\n");
   clone(app, u, &fstop);
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   fstop = ustop + fstop \n");
   sum(app, 1.0, ustop, 1.0, fstop); 

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   r = residual(unext, u) \n");
   residual(app, unext, u, sstatus); 
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   r = fstop - r \n");
   sum(app, 1.0, fstop, -1.0, u); 
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   spatialnorm(r)\n");
   spatialnorm(app, u, &result1);

   if(myaccess != NULL)
   {
      _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   access(r) \n");
      myaccess(app, u, astatus);
   }

   /* We expect the result to be close to zero */
   result_int = (braid_Int) round(-log10(result1));
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   actual output:    spatialnorm(r) approx. 1.0e-%d \n", result_int);
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   expected output:  spatialnorm(r) approx. 0.0 \n\n");

   /* Free variables */
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   free(u)\n");
   myfree(app, u);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   free(unext)\n");
   myfree(app, unext);

   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   free(ustop)\n");
   myfree(app, ustop);
   
   _braid_ParFprintfFlush(fp, myid_x, "   braid_TestResidual:   free(fstop)\n");
   myfree(app, fstop);

   _braid_StatusDestroy(status);

   _braid_TFree(rfactors);
   _braid_TFree(grids);
   _braid_TFree(fine_grid);

   _braid_ParFprintfFlush(fp, myid_x, "Finished braid_TestResidual \n");
   
   return 0;
}

braid_Int
braid_TestAll( braid_App            app,
            MPI_Comm                comm_x,
            FILE                    *fp, 
            braid_Real              t,
            braid_Real              fdt,
            braid_Real              cdt,
            braid_PtFcnInit         myinit,
            braid_PtFcnFree         myfree,
            braid_PtFcnClone        clone,
            braid_PtFcnSum          sum,
            braid_PtFcnSpatialNorm  spatialnorm, 
            braid_PtFcnBufSize      bufsize,
            braid_PtFcnBufPack      bufpack,
            braid_PtFcnBufUnpack    bufunpack,
            braid_PtFcnSCoarsen     coarsen,
            braid_PtFcnSRefine      refine,
            braid_PtFcnResidual     residual,
            braid_PtFcnStep         step)

{
   braid_Int    myid_x, flag = 0, correct = 1;
   MPI_Comm_rank( comm_x, &myid_x );
   
   /** 
    * We set the myaccess parameter to NULL below, because this function is
    * is designed to return only one value, the boolean correct
    **/

   /* Test myinit(), myfree() */
   braid_TestInitAccess( app, comm_x, fp, t, myinit, NULL, myfree);
   braid_TestInitAccess( app, comm_x, fp, fdt, myinit, NULL, myfree);

   /* Test clone() */
   braid_TestClone( app, comm_x, fp, t, myinit, NULL, myfree, clone);
   braid_TestClone( app, comm_x, fp, fdt, myinit, NULL, myfree, clone);

   /* Test sum() */
   braid_TestSum( app, comm_x, fp, t, myinit, NULL, myfree, clone, sum);
   braid_TestSum( app, comm_x, fp, fdt, myinit, NULL, myfree, clone, sum);

   /* Test spatialnorm() */
   flag = braid_TestSpatialNorm( app, comm_x, fp, t, myinit, myfree, clone, sum, spatialnorm);
   if(flag == 0)
   {
      _braid_ParFprintfFlush(fp, myid_x, "-> braid_TestAll:   TestSpatialNorm 1 Failed\n");
      correct = 0;
   }
   flag = braid_TestSpatialNorm( app, comm_x, fp, fdt, myinit, myfree, clone, sum, spatialnorm);
   if(flag == 0)
   {
      _braid_ParFprintfFlush(fp, myid_x, "-> braid_TestAll:   TestSpatialNorm 2 Failed\n");
      correct = 0;
   }

   /* Test bufsize(), bufpack(), bufunpack() */
   flag = braid_TestBuf( app, comm_x, fp, t, myinit, myfree, sum, spatialnorm, bufsize, bufpack, bufunpack);
   if(flag == 0)
   {
      _braid_ParFprintfFlush(fp, myid_x, "-> braid_TestAll:   TestBuf 1 Failed\n");
      correct = 0;
   }
   flag = braid_TestBuf( app, comm_x, fp, fdt, myinit, myfree, sum, spatialnorm, bufsize, bufpack, bufunpack);
   if(flag == 0)
   {
      _braid_ParFprintfFlush(fp, myid_x, "-> braid_TestAll:   TestBuf 2 Failed\n");
      correct = 0;
   }
 
   /* Test coarsen and refine */
   if( (coarsen != NULL) && (refine != NULL) )
   {
      flag = braid_TestCoarsenRefine(app, comm_x, fp, t, fdt, cdt, myinit,
                          NULL, myfree, clone, sum, spatialnorm, coarsen, refine);
      if(flag == 0)
      {
         _braid_ParFprintfFlush(fp, myid_x, "-> braid_TestAll:   TestCoarsenRefine 1 Failed\n");
         correct = 0;
      }
   }

   /* Test step and residual */
   if( (step != NULL) && (residual != NULL) )
   {
      braid_TestResidual(app, comm_x, fp, t, fdt, myinit, NULL, myfree, clone, sum, spatialnorm, residual, step);
   }

   if(correct == 1)
   {
      _braid_ParFprintfFlush(fp, myid_x, "\n\n");
      _braid_ParFprintfFlush(fp, myid_x, "Finished braid_TestAll: no fails detected, however some results must be\nuser-interpreted in the log messages\n\n");
   }
   else
   {
      _braid_ParFprintfFlush(fp, myid_x, "\n\n");
      _braid_ParFprintfFlush(fp, myid_x, "Finished braid_TestAll: some tests failed\n\n");
   }
   
   return correct;
}
