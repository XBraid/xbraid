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

#include "warp.h"
#include "_warp.h"
#include "util.h"

/*--------------------------------------------------------------------------
 * Some simple tests on the init, write and free routines 
 *--------------------------------------------------------------------------*/
warp_Int
warp_TestInitWrite( warp_App              app, 
                    MPI_Comm              comm_x,
                    FILE                 *fp, 
                    warp_Real             t,
                    warp_PtFcnInit        init, 
                    warp_PtFcnWrite       write,
                    warp_PtFcnFree        free)
{
   
   warp_Vector    u ;
   warp_Status    status;
   warp_Int       myid_x;
   
   _warp_InitStatus( 0.0, 0, 0, 0, &status);
   MPI_Comm_rank( comm_x, &myid_x );

   /* Print intro */
   _warp_ParFprintfFlush(fp, myid_x, "\nStarting warp_TestInitWrite\n\n");

   /* Test */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestInitWrite:   Starting Test 1\n");
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestInitWrite:   u = init(t=%1.2e)\n", t);
   init(app, t, &u);
   
   if(write != NULL)
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestInitWrite:   write(u) \n");
      write(app, t, status, u);

      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestInitWrite:   check output: wrote u for initial condition at t=%1.2e. \n\n",t);
   }

   /* Free variables */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestInitWrite:   free(u) \n");
   free(app, u);
   
   _warp_ParFprintfFlush(fp, myid_x, "Finished warp_TestInitWrite\n");

   return 0;
}

warp_Int
warp_TestClone( warp_App              app,  
                MPI_Comm              comm_x,
                FILE                 *fp, 
                warp_Real             t,
                warp_PtFcnInit        init, 
                warp_PtFcnWrite       write,
                warp_PtFcnFree        free,
                warp_PtFcnClone       clone)
{
   
   warp_Vector  u, v;
   warp_Status  status;
   warp_Int     myid_x;
   
   _warp_InitStatus( 0.0, 0, 0, 0, &status);
   MPI_Comm_rank( comm_x, &myid_x );

   /* Print intro */
   _warp_ParFprintfFlush(fp, myid_x, "\nStarting warp_TestClone\n\n");

   /* Test 1 */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestClone:   Starting Test 1\n");
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestClone:   u = init(t=%1.2e)\n", t);
   init(app, t, &u);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestClone:   v = clone(u)\n");
   clone(app, u, &v);
   
   if(write != NULL)
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestClone:   write(u)\n");
      write(app, t, status, u);

      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestClone:   write(v)\n");
      write(app, t, status, v);
      
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestClone:   check output:  wrote u and v for initial condition at t=%1.2e.\n\n", t);

   }

   /* Free variables */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestClone:   free(u)\n");
   free(app, u);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestClone:   free(v)\n");
   free(app, v);

   _warp_ParFprintfFlush(fp, myid_x, "Finished warp_TestClone\n");
   
   return 0;
}



warp_Int
warp_TestSum( warp_App              app, 
              MPI_Comm              comm_x,
              FILE                 *fp, 
              warp_Real             t,
              warp_PtFcnInit        init, 
              warp_PtFcnWrite       write,
              warp_PtFcnFree        free, 
              warp_PtFcnClone       clone,
              warp_PtFcnSum         sum )  
{
   
   warp_Vector  u, v;
   warp_Status  status;
   warp_Int     myid_x;
   
   _warp_InitStatus( 0.0, 0, 0, 0, &status);
   MPI_Comm_rank( comm_x, &myid_x );

   /* Print intro */
   _warp_ParFprintfFlush(fp, myid_x, "\nStarting warp_TestSum\n\n");
   
   /* Test 1 */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestSum:   Starting Test 1\n");
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestSum:   u = init(t=%1.2e)\n", t);
   init(app, t, &u);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestSum:   v = clone(u)\n");
   clone(app, u, &v);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestSum:   v = u - v\n");
   sum(app, 1.0, u, -1.0, v); 

   if(write != NULL)
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestSum:   write(v)\n");
      write(app, t, status, v);
      
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestSum:   check output:  v should equal the zero vector\n\n");
   }

   /* Test 2 */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestSum:   Starting Test 2\n");
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestSum:   v = 2*u + v\n");
   sum(app, 2.0, u, 1.0, v); 

   if(write != NULL)
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestSum:   write(v)\n");
      write(app, t, status, v);
      
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestSum:   write(u)\n");
      write(app, t, status, u);
   }

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestSum:   check output:  v should equal 2*u \n\n");

   /* Free variables */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestSum:   free(u)\n");
   free(app, u);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestSum:   free(v)\n");
   free(app, v);

   _warp_ParFprintfFlush(fp, myid_x, "Finished warp_TestSum\n");

   return 0;
}

warp_Int
warp_TestDot( warp_App              app, 
              MPI_Comm              comm_x,
              FILE                 *fp, 
              warp_Real             t,
              warp_PtFcnInit        init, 
              warp_PtFcnFree        free, 
              warp_PtFcnClone       clone,
              warp_PtFcnSum         sum,  
              warp_PtFcnDot         dot)
{   
   warp_Vector  u, v, w;
   warp_Real    result1, result2, result3;
   warp_Int     myid_x, correct;
   double       wiggle = 1e-12;
   
   MPI_Comm_rank( comm_x, &myid_x );

   /* Initialize the flags */
   correct = 1;
   warp_Int zero_flag = 0;

   /* Print intro */
   _warp_ParFprintfFlush(fp, myid_x, "\nStarting warp_TestDot\n\n");
   
   /* Test 1 */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   Starting Test 1\n");
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   u = init(t=%1.2e)\n", t);
   init(app, t, &u);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   dot(u,u) \n");
   dot(app, u, u, &result1);
   if( fabs(result1) == 0.0)
   {
      zero_flag = 1;
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   Warning:  dot(u,u) = 0.0\n"); 
   }
   else if( isnan(result1) )
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   Warning:  dot(u,u) = nan\n"); 
      correct = 0;
   }

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   v = clone(u)\n");
   clone(app, u, &v);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   v = u - v \n");
   sum(app, 1.0, u, -1.0, v); 

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   dot(v,v) \n");
   dot(app, v, v, &result1);
   if( (fabs(result1) > wiggle) || isnan(result1) )
   {
      correct = 0;
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   Test 1 Failed\n");
   }
   else
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   Test 1 Passed\n");
   }
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   actual output:    dot(v,v) = %1.2e  \n", result1);
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   expected output:  dot(v,v) = 0.0 \n\n");
   

   /* Test 2 */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   Starting Test 2\n");
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   w = clone(u)\n");
   clone(app, u, &w);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   w = u + w \n");
   sum(app, 1.0, u, 1.0, w); 

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   dot(u,u)\n");
   dot(app, u, u, &result1);
   
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   dot(w,w)\n");
   dot(app, w, w, &result2);
   if( (fabs(result2/result1 - 4.0) > wiggle) || isnan(result2/result1) )
   {
      correct = 0;
      if(zero_flag)
      {
         _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   Test 2 Failed, Likely due to u = 0\n");
      }
      else
      {
         _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   Test 2 Failed\n");
      }
   }
   else
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   Test 2 Passed\n");
   }
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   actual output:    dot(w,w) / dot(u,u) = %1.2e / %1.2e = %1.2e \n",
         result2, result1, result2/result1);
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   expected output:  dot(w,w) / dot(u,u) = 4.0 \n\n");

   /* Test 3 */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   Starting Test 3\n");
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   free(w)\n");
   free(app, w);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   w = clone(u)\n");
   clone(app, u, &w);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   w = 0.0*u + 0.5*w \n");
   sum(app, 0.0, u, 0.5, w); 

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   dot(u,u)\n");
   dot(app, u, u, &result1);
   
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   dot(w,w)\n");
   dot(app, w, w, &result2);
   /* Check Result */
   if( (fabs(result2/result1 - 0.25) > wiggle) || isnan(result2/result1) )

   {
      correct = 0;
      if(zero_flag)
      {
         _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   Test 3 Failed, Likely due to u = 0\n");
      }
      else
      {
         _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   Test 3 Failed\n");
      }
   }
   else
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   Test 3 Passed\n");
   }
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   actual output:    dot(w,w) / dot(u,u) = %1.2e / %1.2e = %1.2e \n",
         result2, result1, result2/result1);
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   expected output:  dot(w,w) / dot(u,u) = 0.25 \n\n");

   /* Test 4 */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   Starting Test 4\n");

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   free(w)\n");
   free(app, w);
   
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   w = clone(u)\n");
   clone(app, u, &w);
   
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   w = u + 0.5*w \n");
   sum(app, 1.0, u, 0.5, w); 

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   dot(u,u)\n");
   dot(app, u, u, &result1);
   
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   dot(w,u)\n");
   dot(app, w, u, &result2);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   actual output:    dot(w,u) + dot(u,u) = %1.2e + %1.2e = %1.2e\n", 
      result2, result1, result2+result1);
   
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   v = u + w \n");
   sum(app, 1.0, u, 1.0, w);   
   
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   dot(v,u)\n");
   dot(app, w, u, &result3);

   /* Check Result */
   if( (fabs(result2 + result1 - result3)/fabs(result3) > wiggle) || 
       isnan(result2) || isnan(result1) || isnan(result3) )
   {
      correct = 0;
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   Test 4 Failed\n");
   }
   else
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   Test 4 Passed\n");
   }
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   actual output:    dot(v,u) = %1.2e  \n", result3);
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   expected output:  dot(v,u) = dot(w,u) + dot(u,u) \n\n");
 

   /* Free variables */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   free(u)\n");
   free(app, u);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   free(v)\n");
   free(app, v);
   
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestDot:   free(w)\n");
   free(app, w);

   if(correct == 1) 
      _warp_ParFprintfFlush(fp, myid_x, "Finished warp_TestDot: all tests passed successfully\n");
   else
   {
      if(zero_flag)
      {
         _warp_ParFprintfFlush(fp, myid_x, "Finished warp_TestDot: some tests failed, possibly due to u = 0\n");
      }
      else
      {
         _warp_ParFprintfFlush(fp, myid_x, "Finished warp_TestDot: some tests failed\n");
      }
   }

   return correct;
}

warp_Int
warp_TestBuf( warp_App              app,
              MPI_Comm              comm_x,
              FILE                 *fp, 
              warp_Real             t,
              warp_PtFcnInit        init,
              warp_PtFcnFree        free,
              warp_PtFcnSum         sum,  
              warp_PtFcnDot         dot,
              warp_PtFcnBufSize     bufsize,
              warp_PtFcnBufPack     bufpack,
              warp_PtFcnBufUnpack   bufunpack)
{   
   warp_Vector  u, v;
   warp_Real    result1;
   warp_Int     myid_x, size, correct;
   void        *buffer;
   double       wiggle = 1e-12;
   
   MPI_Comm_rank( comm_x, &myid_x );

   /* Initialize the correct flag */
   correct = 1;

   /* Print intro */
   _warp_ParFprintfFlush(fp, myid_x, "\nStarting warp_TestBuf\n\n");
   
   /* Test 1 */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   Starting Test 1\n");
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   u = init(t=%1.2e)\n", t);
   init(app, t, &u);
   
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   dot(u,u) \n");
   dot(app, u, u, &result1);
   if( fabs(result1) == 0.0)
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   Warning:  dot(u,u) = 0.0\n"); 
   }
   else if( isnan(result1) )
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   Warning:  dot(u,u) = nan\n"); 
      correct = 0;
   }

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   size = bufsize()\n");
   bufsize(app, &size);
   
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   buffer = malloc(size)\n");
   buffer = malloc(size);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   buffer = bufpack(u, buffer))\n");
   bufpack(app, u, buffer);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   v = bufunpack(buffer)\n");
   bufunpack(app, buffer, &v);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   v = u - v \n");
   sum(app, 1.0, u, -1.0, v); 

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   dot(v,v) \n");
   dot(app, v, v, &result1);
   if( (fabs(result1) > wiggle) || isnan(result1) )

   {
      correct = 0;
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   Test 1 Failed\n");
   }
   else
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   Test 1 Passed\n");
   }
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   actual output:    dot(v,v) = %1.2e  \n", result1);
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   expected output:  dot(v,v) = 0.0 \n\n");
   
   /* Free variables */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   free(u)\n");
   free(app, u);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestBuf:   free(v)\n");
   free(app, v);
   
   if(correct == 1) 
      _warp_ParFprintfFlush(fp, myid_x, "Finished warp_TestBuf: all tests passed successfully\n");
   else      
      _warp_ParFprintfFlush(fp, myid_x, "Finished warp_TestBuf: some tests failed\n");

   return correct;
}

warp_Int
warp_TestCoarsenRefine( warp_App          app,
                        MPI_Comm          comm_x,
                        FILE             *fp, 
                        warp_Real         t,
                        warp_Real         fdt,
                        warp_Real         cdt,
                        warp_PtFcnInit    init,
                        warp_PtFcnWrite   write,
                        warp_PtFcnFree    free,
                        warp_PtFcnClone   clone,
                        warp_PtFcnSum     sum,
                        warp_PtFcnDot     dot,
                        warp_PtFcnCoarsen coarsen,
                        warp_PtFcnRefine  refine)
 {   
   warp_Vector  u, v, w, uc, vc, wc;
   warp_Real    result1;
   warp_Int     myid_x, level, correct;
   warp_Status  status;
   
   MPI_Comm_rank( comm_x, &myid_x );

   /* Initialize the correct flag */
   correct = 1;

   /* Print intro */
   _warp_ParFprintfFlush(fp, myid_x, "\nStarting warp_TestCoarsenRefine\n\n");
   
   /* Test 1 */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   Starting Test 1\n");
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   u = init(t=%1.2e)\n", t);
   init(app, t, &u);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   dot(u,u) \n");
   dot(app, u, u, &result1);
   if( fabs(result1) == 0.0)
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   Warning:  dot(u,u) = 0.0\n"); 
   }
   else if( isnan(result1) )
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   Warning:  dot(u,u) = nan\n"); 
      correct = 0;
   }

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   uc = coarsen(u)\n");
   coarsen(app, t, t-fdt, t+fdt, t-cdt, t+cdt, u, &uc); 

   if(write != NULL)
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   write(uc) \n");
      level = 1;
      _warp_InitStatus( 0.0, 0, level, 0, &status);
      write(app, t, status, uc);

      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   write(u) \n");
      level = 0;
      _warp_InitStatus( 0.0, 0, level, 0, &status);
      write(app, t, status, u);
   }

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   actual output:   wrote u and spatially coarsened u \n\n");

   /* Test 2 */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   Starting Test 2\n");
   
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   v = clone(u)\n");
   clone(app, u, &v);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   vc = coarsen(v)\n");
   coarsen(app, t, t-fdt, t+fdt, t-cdt, t+cdt, v, &vc); 
   
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   wc = clone(vc)\n");
   clone(app, vc, &wc);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   wc = uc - wc \n");
   sum(app, 1.0, uc, -1.0, wc); 

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   dot(wc,wc)\n");
   dot(app, wc, wc, &result1);
   
   /* We expect exact equality between uc and vc */
   if( (fabs(result1) != 0.0) || isnan(result1) )
   {
      correct = 0;
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   Test 2 Failed\n");
   }
   else
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   Test 2 Passed\n");
   }
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   actual output:    dot(wc,wc) = %1.2e \n", result1);
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   expected output:  dot(wc,wc) = 0.0 \n\n");

   /* Test 3 */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   Starting Test 3\n");
   
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   w = clone(u)\n");
   clone(app, u, &w);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   free(u)\n");
   free(app, u);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   free(v)\n");
   free(app, v);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   v = refine(vc)\n");
   refine(app, t, t-fdt, t+fdt, t-cdt, t+cdt, vc, &v); 

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   u = refine(uc)\n");
   refine(app, t, t-fdt, t+fdt, t-cdt, t+cdt, uc, &u); 

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   v = u - v \n");
   sum(app, 1.0, u, -1.0, v); 

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   dot(v,v)\n");
   dot(app, v, v, &result1);
   
   /* We expect exact equality between u and v */
   if( (fabs(result1) != 0.0) || isnan(result1) )
   {
      correct = 0;
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   Test 3 Failed\n");
   }
   else
   {
      _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   Test 3 Passed\n");
   }
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   actual output:    dot(v,v) = %1.2e \n", result1);
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   expected output:  dot(v,v) = 0.0 \n\n");

   /* Test 4 */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   Starting Test 4\n");
   
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   w = w - u \n");
   sum(app, 1.0, u, -1.0, w); 

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   dot(w,w)\n");
   dot(app, w, w, &result1);
   
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   actual output:    dot(w,w) = %1.2e \n", result1);
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   %s%s%s", 
                   "expected output:  For simple interpolation formulas\n",
                   "                             (e.g., bilinear) and a known function\n",
                   "                             (e.g., constant), dot(w,w) should = 0\n\n");

   /* Free variables */
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   free(u)\n");
   free(app, u);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   free(v)\n");
   free(app, v);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   free(w)\n");
   free(app, w);
   
   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   free(uc)\n");
   free(app, uc);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   free(vc)\n");
   free(app, vc);

   _warp_ParFprintfFlush(fp, myid_x, "   warp_TestCoarsenRefine:   free(wc)\n");
   free(app, wc);


   if(correct == 1) 
      _warp_ParFprintfFlush(fp, myid_x, "Finished warp_TestCoarsenRefine: all tests passed successfully\n");
   else      
      _warp_ParFprintfFlush(fp, myid_x, "Finished warp_TestCoarsenRefine: some tests failed\n");
   
   return correct;
}

warp_Int
warp_TestAll( warp_App             app,
              MPI_Comm             comm_x,
              FILE                *fp, 
              warp_Real            t,
              warp_Real            fdt,
              warp_Real            cdt,
              warp_PtFcnInit       init,
              warp_PtFcnFree       free,
              warp_PtFcnClone      clone,
              warp_PtFcnSum        sum,
              warp_PtFcnDot        dot,
              warp_PtFcnBufSize    bufsize,
              warp_PtFcnBufPack    bufpack,
              warp_PtFcnBufUnpack  bufunpack,
              warp_PtFcnCoarsen    coarsen,
              warp_PtFcnRefine     refine)
{
   warp_Int    myid_x, flag = 0, correct = 1;
   
   MPI_Comm_rank( comm_x, &myid_x );
   
   /** 
    * We set the write parameter to NULL below, because this function is
    * is designed to return only one value, the boolean correct
    **/

   /* Test init(), free() */
   warp_TestInitWrite( app, comm_x, fp, t, init, NULL, free);
   warp_TestInitWrite( app, comm_x, fp, fdt, init, NULL, free);

   /* Test clone() */
   warp_TestClone( app, comm_x, fp, t, init, NULL, free, clone);
   warp_TestClone( app, comm_x, fp, fdt, init, NULL, free, clone);

   /* Test sum() */
   warp_TestSum( app, comm_x, fp, t, init, NULL, free, clone, sum);
   warp_TestSum( app, comm_x, fp, fdt, init, NULL, free, clone, sum);

   /* Test dot() */
   flag = warp_TestDot( app, comm_x, fp, t, init, free, clone, sum, dot);
   if(flag == 0)
   {
      _warp_ParFprintfFlush(fp, myid_x, "-> warp_TestAll:   TestDot 1 Failed\n");
      correct = 0;
   }
   flag = warp_TestDot( app, comm_x, fp, fdt, init, free, clone, sum, dot);
   if(flag == 0)
   {
      _warp_ParFprintfFlush(fp, myid_x, "-> warp_TestAll:   TestDot 2 Failed\n");
      correct = 0;
   }

   /* Test bufsize(), bufpack(), bufunpack() */
   flag = warp_TestBuf( app, comm_x, fp, t, init, free, sum, dot, bufsize, bufpack, bufunpack);
   if(flag == 0)
   {
      _warp_ParFprintfFlush(fp, myid_x, "-> warp_TestAll:   TestBuf 1 Failed\n");
      correct = 0;
   }
   flag = warp_TestBuf( app, comm_x, fp, fdt, init, free, sum, dot, bufsize, bufpack, bufunpack);
   if(flag == 0)
   {
      _warp_ParFprintfFlush(fp, myid_x, "-> warp_TestAll:   TestBuf 2 Failed\n");
      correct = 0;
   }
 
   /* Test coarsen and refine */
   if( (coarsen != NULL) && (refine != NULL) )
   {
      flag = warp_TestCoarsenRefine(app, comm_x, fp, t, fdt, cdt, init,
                          NULL, free, clone, sum, dot, coarsen, refine);
      if(flag == 0)
      {
         _warp_ParFprintfFlush(fp, myid_x, "-> warp_TestAll:   TestCoarsenRefine 1 Failed\n");
         correct = 0;
      }
   }

   if(correct == 1)
   {
      _warp_ParFprintfFlush(fp, myid_x, "\n\n");
      _warp_ParFprintfFlush(fp, myid_x, "Finished warp_TestAll: all tests passed successfully\n\n");
   }
   else
   {
      _warp_ParFprintfFlush(fp, myid_x, "\n\n");
      _warp_ParFprintfFlush(fp, myid_x, "Finished warp_TestAll: some tests failed\n\n");
   }
   
   return correct;
}
