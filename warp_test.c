#include "warp.h"
#include "_warp.h"
#include "util.h"

/*--------------------------------------------------------------------------
 * Some simple tests on the init, write and free routines 
 *--------------------------------------------------------------------------*/
warp_Int
warp_TestInitWrite( warp_App              app, 
                    MPI_Comm              comm_x,
                    warp_Real             t,
                    warp_PtFcnInit        init, 
                    warp_PtFcnWrite       write,
                    warp_PtFcnFree        free)
{
   
   warp_Vector    u ;
   warp_Status    status;
   warp_Int       myid_x;
   char           header[255];
   char           message[255];
   
   _warp_InitStatus( 0.0, 0, 0, 0, &status);
   MPI_Comm_rank( comm_x, &myid_x );
   sprintf(header,  "   warp_TestInitWrite:   ");

   /* Print intro */
   sprintf(message, "\nStarting warp_TestInitWrite\n\n");
   _warp_ParFprintfFlush(stdout, "", message, myid_x);

   /* Test */
   sprintf(message, "Starting Test 1\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "u = init(t=%1.2e)\n", t);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   init(app, t, &u);

   sprintf(message, "write(u) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   write(app, t, status, u);

   sprintf(message, "check output: wrote u for initial condition at t=%1.2e. \n\n",t);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   /* Free variables */
   sprintf(message, "free(u) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, u);
   
   sprintf(message, "Finished\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   return 0;
}

warp_Int
warp_TestClone( warp_App              app,  
                MPI_Comm              comm_x,
                warp_Real             t,
                warp_PtFcnInit        init, 
                warp_PtFcnWrite       write,
                warp_PtFcnFree        free,
                warp_PtFcnClone       clone)
{
   
   warp_Vector  u, v;
   warp_Status  status;
   warp_Int     myid_x;
   char         header[255];
   char         message[255];
   
   _warp_InitStatus( 0.0, 0, 0, 0, &status);
   MPI_Comm_rank( comm_x, &myid_x );
   sprintf(header,  "   warp_TestClone:   ");

   /* Print intro */
   sprintf(message, "\nStarting warp_TestClone\n\n");
   _warp_ParFprintfFlush(stdout, "", message, myid_x);

   /* Test 1 */
   sprintf(message, "Starting Test 1\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "u = init(t=%1.2e)\n", t);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   init(app, t, &u);

   sprintf(message, "v = clone(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, u, &v);
   
   sprintf(message, "write(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   write(app, t, status, u);

   sprintf(message, "write(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   write(app, t, status, v);
   
   sprintf(message, "check output:  wrote u and v for initial condition at t=%1.2e.\n\n", t);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   
   /* Free variables */
   sprintf(message, "free(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, u);

   sprintf(message, "free(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, v);

   sprintf(message, "Finished\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   
   return 0;
}



warp_Int
warp_TestSum( warp_App              app, 
              MPI_Comm              comm_x,
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
   char         header[255];
   char         message[255];
   
   _warp_InitStatus( 0.0, 0, 0, 0, &status);
   MPI_Comm_rank( comm_x, &myid_x );
   sprintf(header,  "   warp_TestSum:   ");

   /* Print intro */
   sprintf(message, "\nStarting warp_TestSum\n\n");
   _warp_ParFprintfFlush(stdout, "", message, myid_x);
   
   /* Test 1 */
   sprintf(message, "Starting Test 1\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "u = init(t=%1.2e)\n", t);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   init(app, t, &u);

   sprintf(message, "v = clone(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, u, &v);

   sprintf(message, "v = u - v\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, u, -1.0, v); 

   sprintf(message, "write(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   write(app, t, status, v);
   
   sprintf(message, "check output:  v should equal the zero vector\n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   
   /* Test 2 */
   sprintf(message, "Starting Test 2\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "v = 2*u + v\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 2.0, u, 1.0, v); 

   sprintf(message, "write(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   write(app, t, status, v);
   
   sprintf(message, "write(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   write(app, t, status, u);
   
   sprintf(message, "check output:  v should equal 2*u \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   /* Free variables */
   sprintf(message, "free(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, u);

   sprintf(message, "free(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, v);

   sprintf(message, "Finished\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   return 0;
}

warp_Int
warp_TestDot( warp_App              app, 
              MPI_Comm              comm_x,
              warp_Real             t,
              warp_PtFcnInit        init, 
              warp_PtFcnFree        free, 
              warp_PtFcnClone       clone,
              warp_PtFcnSum         sum,  
              warp_PtFcnDot         dot,
              warp_Int             *correct)  
{   
   warp_Vector  u, v, w;
   warp_Real    result1, result2, result3;
   warp_Int     myid_x;
   char         header[255];
   char         message[255];
   double       wiggle = 1e-12;
   
   MPI_Comm_rank( comm_x, &myid_x );
   sprintf(header,  "   warp_TestDot:   ");

   /* Initialize the correct flag */
   *correct = 1;

   /* Print intro */
   sprintf(message, "\nStarting warp_TestDot\n\n");
   _warp_ParFprintfFlush(stdout, "", message, myid_x);
   
   /* Test 1 */
   sprintf(message, "Starting Test 1\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "u = init(t=%1.2e)\n", t);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   init(app, t, &u);

   sprintf(message, "dot(u,u) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, u, u, &result1);
   if( fabs(result1) == 0.0)
   {
      sprintf(message, "Warning:  dot(u,u) = 0.0\n"); 
      _warp_ParFprintfFlush(stdout, header, message, myid_x);
   }
   else if( isnan(result1) )
   {
      sprintf(message, "Warning:  dot(u,u) = nan\n"); 
      _warp_ParFprintfFlush(stdout, header, message, myid_x);
      *correct = 0;
   }

   sprintf(message, "v = clone(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, u, &v);

   sprintf(message, "v = u - v \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, u, -1.0, v); 

   sprintf(message, "dot(v,v) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, v, v, &result1);
   if( (fabs(result1) > wiggle) || isnan(result1) )
   {
      *correct = 0;
      _warp_ParFprintfFlush(stdout, header, "Test 1 Failed\n", myid_x);
   }
   else
   {
      _warp_ParFprintfFlush(stdout, header, "Test 1 Passed\n", myid_x);
   }
   sprintf(message, "actual output:    dot(v,v) = %1.2e  \n", result1);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "expected output:  dot(v,v) = 0.0 \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   

   /* Test 2 */
   sprintf(message, "Starting Test 2\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "w = clone(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, u, &w);

   sprintf(message, "w = u + w \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, u, 1.0, w); 

   sprintf(message, "dot(u,u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, u, u, &result1);
   
   sprintf(message, "dot(w,w)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, w, w, &result2);
   if( (fabs(result2/result1 - 4.0) > wiggle) || isnan(result2/result1) )
   {
      *correct = 0;
      _warp_ParFprintfFlush(stdout, header, "Test 2 Failed\n", myid_x);
   }
   else
   {
      _warp_ParFprintfFlush(stdout, header, "Test 2 Passed\n", myid_x);
   }
   sprintf(message, "actual output:    dot(w,w) / dot(u,u) = %1.2e / %1.2e = %1.2e \n",
         result2, result1, result2/result1);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "expected output:  dot(w,w) / dot(u,u) = 4.0 \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   /* Test 3 */
   sprintf(message, "Starting Test 3\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "free(w)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, w);

   sprintf(message, "w = clone(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, u, &w);

   sprintf(message, "w = 0.0*u + 0.5*w \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 0.0, u, 0.5, w); 

   sprintf(message, "dot(u,u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, u, u, &result1);
   
   sprintf(message, "dot(w,w)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, w, w, &result2);
   /* Check Result */
   if( (fabs(result2/result1 - 0.25) > wiggle) || isnan(result2/result1) )

   {
      *correct = 0;
      _warp_ParFprintfFlush(stdout, header, "Test 3 Failed\n", myid_x);
   }
   else
   {
      _warp_ParFprintfFlush(stdout, header, "Test 3 Passed\n", myid_x);
   }
   sprintf(message, "actual output:    dot(w,w) / dot(u,u) = %1.2e / %1.2e = %1.2e \n",
         result2, result1, result2/result1);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "expected output:  dot(w,w) / dot(u,u) = 0.25 \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   /* Test 4 */
   sprintf(message, "Starting Test 4\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   sprintf(message, "free(w)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, w);
   
   sprintf(message, "w = clone(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, u, &w);
   
   sprintf(message, "w = u + 0.5*w \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, u, 0.5, w); 

   sprintf(message, "dot(u,u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, u, u, &result1);
   
   sprintf(message, "dot(w,u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, w, u, &result2);

   sprintf(message, "actual output:    dot(w,u) + dot(u,u) = %1.2e + %1.2e = %1.2e\n", result2, result1, result2+result1);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   
   sprintf(message, "v = u + w \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, u, 1.0, w);   
   
   sprintf(message, "dot(v,u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, w, u, &result3);

   /* Check Result */
   if( (fabs(result2 + result1 - result3)/fabs(result3) > wiggle) || 
       isnan(result2) || isnan(result1) || isnan(result3) )
   {
      *correct = 0;
      _warp_ParFprintfFlush(stdout, header, "Test 4 Failed\n", myid_x);
   }
   else
   {
      _warp_ParFprintfFlush(stdout, header, "Test 4 Passed\n", myid_x);
   }
   sprintf(message, "actual output:    dot(v,u) = %1.2e  \n", result3);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "expected output:  dot(v,u) equals dot(w,u) + dot(u,u) \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
 

   /* Free variables */
   sprintf(message, "free(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, u);

   sprintf(message, "free(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, v);
   
   sprintf(message, "free(w)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, w);

   if(*correct == 1) 
      _warp_ParFprintfFlush(stdout, header, "Finished, all tests passed successfully\n", myid_x);
   else      
      _warp_ParFprintfFlush(stdout, header, "Finished, some tests failed\n", myid_x);

   return 0;
}

warp_Int
warp_TestBuf( warp_App              app,
              MPI_Comm              comm_x,
              warp_Real             t,
              warp_PtFcnInit        init,
              warp_PtFcnFree        free,
              warp_PtFcnSum         sum,  
              warp_PtFcnDot         dot,
              warp_PtFcnBufSize     bufsize,
              warp_PtFcnBufPack     bufpack,
              warp_PtFcnBufUnpack   bufunpack,
              warp_Int             *correct)
{   
   warp_Vector  u, v;
   warp_Real    result1;
   warp_Int     myid_x, size;
   void        *buffer;
   char         header[255];
   char         message[255];
   double       wiggle = 1e-12;
   
   MPI_Comm_rank( comm_x, &myid_x );
   sprintf(header,  "   warp_TestBuf:   ");

   /* Initialize the correct flag */
   *correct = 1;

   /* Print intro */
   sprintf(message, "\nStarting warp_TestBuf\n\n");
   _warp_ParFprintfFlush(stdout, "", message, myid_x);
   
   /* Test 1 */
   sprintf(message, "Starting Test 1\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "u = init(t=%1.2e)\n", t);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   init(app, t, &u);
   
   sprintf(message, "dot(u,u) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, u, u, &result1);
   if( fabs(result1) == 0.0)
   {
      sprintf(message, "Warning:  dot(u,u) = 0.0\n"); 
      _warp_ParFprintfFlush(stdout, header, message, myid_x);
   }
   else if( isnan(result1) )
   {
      sprintf(message, "Warning:  dot(u,u) = nan\n"); 
      _warp_ParFprintfFlush(stdout, header, message, myid_x);
      *correct = 0;
   }

   sprintf(message, "size = bufsize()\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   bufsize(app, &size);
   
   sprintf(message, "buffer = malloc(size)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   buffer = malloc(size);

   sprintf(message, "buffer = bufpack(u, buffer))\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   bufpack(app, u, buffer);

   sprintf(message, "v = bufunpack(buffer)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   bufunpack(app, buffer, &v);

   sprintf(message, "v = u - v \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, u, -1.0, v); 

   sprintf(message, "dot(v,v) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, v, v, &result1);
   if( (fabs(result1) > wiggle) || isnan(result1) )

   {
      *correct = 0;
      _warp_ParFprintfFlush(stdout, header, "Test 1 Failed\n", myid_x);
   }
   else
   {
      _warp_ParFprintfFlush(stdout, header, "Test 1 Passed\n", myid_x);
   }
   sprintf(message, "actual output:    dot(v,v) = %1.2e  \n", result1);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "expected output:  dot(v,v) = 0.0 \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   
   /* Free variables */
   sprintf(message, "free(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, u);

   sprintf(message, "free(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, v);
   
   if(*correct == 1) 
      _warp_ParFprintfFlush(stdout, header, "Finished, all tests passed successfully\n", myid_x);
   else      
      _warp_ParFprintfFlush(stdout, header, "Finished, some tests failed\n", myid_x);

   return 0;
}

warp_Int
warp_TestCoarsenRefine( warp_App          app,
                        MPI_Comm          comm_x,
                        warp_Real         t,
                        warp_Real         f_tminus,
                        warp_Real         f_tplus,
                        warp_Real         c_tminus,
                        warp_Real         c_tplus,
                        warp_PtFcnInit    init,
                        warp_PtFcnWrite   write,
                        warp_PtFcnFree    free,
                        warp_PtFcnClone   clone,
                        warp_PtFcnSum     sum,
                        warp_PtFcnDot     dot,
                        warp_PtFcnCoarsen coarsen,
                        warp_PtFcnRefine  refine,
                        warp_Int         *correct)
 {   
   warp_Vector  u, v, w, uc, vc, wc;
   warp_Real    result1;
   warp_Int     myid_x, level;
   char         header[255];
   char         message[255];
   warp_Status  status;
   
   MPI_Comm_rank( comm_x, &myid_x );
   sprintf(header,  "   warp_TestCoarsenRefine:   ");

   /* Initialize the correct flag */
   *correct = 1;

   /* Print intro */
   sprintf(message, "\nStarting warp_TestCoarsenRefine\n\n");
   _warp_ParFprintfFlush(stdout, "", message, myid_x);
   
   /* Test 1 */
   sprintf(message, "Starting Test 1\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "u = init(t=%1.2e)\n", t);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   init(app, t, &u);

   sprintf(message, "dot(u,u) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, u, u, &result1);
   if( fabs(result1) == 0.0)
   {
      sprintf(message, "Warning:  dot(u,u) = 0.0\n"); 
      _warp_ParFprintfFlush(stdout, header, message, myid_x);
   }
   else if( isnan(result1) )
   {
      sprintf(message, "Warning:  dot(u,u) = nan\n"); 
      _warp_ParFprintfFlush(stdout, header, message, myid_x);
      *correct = 0;
   }

   sprintf(message, "uc = coarsen(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   coarsen(app, t, f_tminus, f_tplus, c_tminus, c_tplus, u, &uc); 

   sprintf(message, "write(uc) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   level = 1;
   _warp_InitStatus( 0.0, 0, level, 0, &status);
   write(app, t, status, uc);

   sprintf(message, "write(u) \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   level = 0;
   _warp_InitStatus( 0.0, 0, level, 0, &status);
   write(app, t, status, u);

   sprintf(message, "actual output:   wrote u and spatially coarsened u \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   /* Test 2 */
   sprintf(message, "Starting Test 2\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   
   sprintf(message, "v = clone(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, u, &v);

   sprintf(message, "vc = coarsen(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   coarsen(app, t, f_tminus, f_tplus, c_tminus, c_tplus, v, &vc); 
   
   sprintf(message, "wc = clone(vc)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, vc, &wc);

   sprintf(message, "wc = uc - wc \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, uc, -1.0, wc); 

   sprintf(message, "dot(wc,wc)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, wc, wc, &result1);
   
   /* We expect exact equality between uc and vc */
   if( (fabs(result1) != 0.0) || isnan(result1) )
   {
      *correct = 0;
      _warp_ParFprintfFlush(stdout, header, "Test 2 Failed\n", myid_x);
   }
   else
   {
      _warp_ParFprintfFlush(stdout, header, "Test 2 Passed\n", myid_x);
   }
   sprintf(message, "actual output:    dot(wc,wc) = %1.2e \n", result1);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "expected output:  dot(wc,wc) = 0.0 \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   /* Test 3 */
   sprintf(message, "Starting Test 3\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   
   sprintf(message, "w = clone(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   clone(app, u, &w);

   sprintf(message, "free(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, u);

   sprintf(message, "free(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, v);

   sprintf(message, "v = refine(vc)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   refine(app, t, f_tminus, f_tplus, c_tminus, c_tplus, vc, &v); 

   sprintf(message, "u = refine(uc)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   refine(app, t, f_tminus, f_tplus, c_tminus, c_tplus, uc, &u); 

   sprintf(message, "v = u - v \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, u, -1.0, v); 

   sprintf(message, "dot(v,v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, v, v, &result1);
   
   /* We expect exact equality between u and v */
   if( (fabs(result1) != 0.0) || isnan(result1) )
   {
      *correct = 0;
      _warp_ParFprintfFlush(stdout, header, "Test 3 Failed\n", myid_x);
   }
   else
   {
      _warp_ParFprintfFlush(stdout, header, "Test 3 Passed\n", myid_x);
   }
   sprintf(message, "actual output:    dot(v,v) = %1.2e \n", result1);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "expected output:  dot(v,v) = 0.0 \n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   /* Test 4 */
   sprintf(message, "Starting Test 4\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   
   sprintf(message, "w = w - u \n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sum(app, 1.0, u, -1.0, w); 

   sprintf(message, "dot(w,w)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   dot(app, w, w, &result1);
   
   sprintf(message, "actual output:    dot(w,w) = %1.2e \n", result1);
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   sprintf(message, "%s%s%s", "expected output:  For simple interpolation formulas\n",
                   "                             (e.g., bilinear) and a known function\n",
                   "                             (e.g., constant), dot(w,w) should = 0\n\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);

   /* Free variables */
   sprintf(message, "free(u)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, u);

   sprintf(message, "free(v)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, v);

   sprintf(message, "free(w)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, w);
   
   sprintf(message, "free(uc)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, uc);

   sprintf(message, "free(vc)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, vc);

   sprintf(message, "free(wc)\n");
   _warp_ParFprintfFlush(stdout, header, message, myid_x);
   free(app, wc);


   if(*correct == 1) 
      _warp_ParFprintfFlush(stdout, header, "Finished, all tests passed successfully\n", myid_x);
   else      
      _warp_ParFprintfFlush(stdout, header, "Finished, some tests failed\n", myid_x);
   
   return 0;
}


