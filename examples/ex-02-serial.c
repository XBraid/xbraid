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

/**
 * Example:       ex-02-serial.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-02-serial
 *
 * Help with:     ex-02-serial -help
 *
 * Sample run:    ex-02-serial -ntime 64 -nspace 17
 *
 * Description:   Solves the 1D heat equation on a regular grid in space and time,
 *                using backward Euler in time and classic second order 
 *                finite-differencing in space.
 *
 *                The difference with ex-02.c is that this is a strictly sequential 
 *                time integration code, created for comparison with ex-02.c.  Also,
 *                there is no parallelism in this example.
 *
 *                For details on the discretization, see the header in ex-02.c.
 **/

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "ex-02-lib.c"

int main (int argc, char *argv[])
{
   int             step, arg_index;
   double          t, error, matrix[3], elapsed;
   double         *values, *temp;
   char            filename[255];
   struct timeval  begin, end;

   /* Define space-time domain */
   double    deltaX, deltaT;
   double    tstart        =  0.0;
   double    tstop         =  2*PI;
   int       ntime         =  64;
   double    xstart        =  0.0;
   double    xstop         =  PI;
   int       nspace        =  33;

   /* Parse command line */
   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         printf("\n");
         printf(" Solve the 1D heat equation on space-time domain:  [0, PI] x [0, 2*PI]\n");
         printf(" with exact solution u(t,x) = sin(x)*cos(t) \n\n");
         printf("  -ntime <ntime>       : set num points in time\n");
         printf("  -nspace <nspace>     : set num points in space\n\n");
         exit(1);
      }
      else if ( strcmp(argv[arg_index], "-ntime") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nspace") == 0 )
      {
         arg_index++;
         nspace = atoi(argv[arg_index++]);
      }
      else
      {
         printf("ABORTING: incorrect command line parameter %s\n", argv[arg_index]);
         return (0);
      }
   }

   printf("\n  --------------------- \n");
   printf("  Begin simulation \n");
   printf("  --------------------- \n\n");
   gettimeofday(&begin, NULL);

   /* Setup workspace */
   values = (double*) malloc( nspace*sizeof(double) );
   temp   = (double*) malloc( nspace*sizeof(double) );
   deltaX =  (xstop-xstart)/(nspace-1.0);
   deltaT =  (tstop-tstart)/ntime;

   /* Initialize solution at time t=0.0 */
   get_solution(values, nspace, 0.0, xstart, deltaX);
   sprintf(filename, "%s.%07d.%05d", "ex-02-serial.out", 0, 0);
   save_solution(filename, values, nspace, xstart, xstop, ntime, tstart, tstop);

   /* Main time step loop */
   t = tstart;
   for(step=1; step <= ntime; step++)
   {
      if((step % 25) == 0) printf("  Computing step %d\n",step);
      
      /* Take Step */
      t = t + deltaT;
      take_step(values, nspace, t, xstart, deltaX, deltaT, matrix, temp);
      
      /* Output Solution */
      sprintf(filename, "%s.%07d.%05d", "ex-02-serial.out", step, 0);
      save_solution(filename, values, nspace, xstart, xstop, ntime, tstart, tstop);
   }
   error = compute_error_norm(values, xstart, xstop, nspace, tstop);
   gettimeofday(&end, NULL);
   elapsed = (end.tv_sec - begin.tv_sec) + (end.tv_usec - begin.tv_usec)/1000000.0;

   printf("\n  --------------------- \n");
   printf("  End simulation \n");
   printf("  --------------------- \n\n");
   printf("  Start time                           %1.5e\n", tstart);
   printf("  Stop time                            %1.5e\n", tstop);
   printf("  Time step size                       %1.5e\n", deltaT);
   printf("  Time steps taken:                    %d\n\n", ntime);
   printf("  Spatial points:                      %d\n", nspace);
   printf("  Spatial mesh size:                   %1.2e\n", deltaX);           
   printf("  CFL ratio dt/dx^2:                   %1.2e\n\n", deltaT/(deltaX*deltaX)); 
   printf("  Wall-clock run-time:                 %.2f\n", elapsed);
   printf("  Discretization error at final time:  %1.4e\n", error);
   
   free(values);
   free(temp);
   return 0;
}
