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


/*
   Example 01

   Compile with: make ex-01

   Sample run:   mpirun -np 2 ex-01

   Description:

   Solve the scalar ODE: u' = lambda u
   in a very simplified XBraid setting.


   When run with the default 10 time steps, the solution is:
   $ ./ex-01
   $ cat ex-01.out.00*
     1.00000000000000e+00
     5.00000000000000e-01
     2.50000000000000e-01
     1.25000000000000e-01
     6.25000000000000e-02
     3.12500000000000e-02
     1.56250000000000e-02
     7.81250000000000e-03
     3.90625000000000e-03
     1.95312500000000e-03
     9.76562500000000e-04
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.h"

/*--------------------------------------------------------------------------
 * User-defined routines and structures
 *--------------------------------------------------------------------------*/

/* App structure can contain anything, and be named anything as well */
typedef struct _braid_App_struct
{
   int       rank;
} my_App;

/* Vector structure can contain anything, and be name anything as well */
typedef struct _braid_Vector_struct
{
   double value;
} my_Vector;

int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);

   /* On the finest grid, each value is half the previous value */
   (u->value) = pow(0.5, tstop-tstart)*(u->value);
   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   if (t == 0.0) /* Initial condition */
   {
      (u->value) = 1.0;
   }
   else /* All other time points set to arbitrary value */
   {
      (u->value) = 0.456;
   }
   *u_ptr = u;

   return 0;
}

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   my_Vector *v;

   v = (my_Vector *) malloc(sizeof(my_Vector));
   (v->value) = (u->value);
   *v_ptr = v;

   return 0;
}

int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u);
   return 0;
}

int
my_Sum(braid_App     app,
       double        alpha,
       braid_Vector  x,
       double        beta,
       braid_Vector  y)
{
   (y->value) = alpha*(x->value) + beta*(y->value);
   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   double dot;

   dot = (u->value)*(u->value);
   *norm_ptr = sqrt(dot);
   return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int        index;
   char       filename[255];
   FILE      *file;
   
   braid_AccessStatusGetTIndex(astatus, &index);

   sprintf(filename, "%s.%04d.%03d", "ex-01.out", index, app->rank);
   file = fopen(filename, "w");
   fprintf(file, "%.14e\n", (u->value));
   fflush(file);
   fclose(file);

   return 0;
}

int
my_BufSize(braid_App          app,
           int                *size_ptr,
           braid_BufferStatus bstatus)
{
   *size_ptr = sizeof(double);
   return 0;
}

int
my_BufPack(braid_App          app,
           braid_Vector       u,
           void               *buffer,
           braid_BufferStatus bstatus)
{
   double *dbuffer = buffer;

   dbuffer[0] = (u->value);
   braid_BufferStatusSetSize( bstatus, sizeof(double) );

   return 0;
}

int
my_BufUnpack(braid_App          app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus bstatus)
{
   double    *dbuffer = buffer;
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   (u->value) = dbuffer[0];
   *u_ptr = u;

   return 0;
}

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
   braid_Core    core;
   my_App       *app;
   double        tstart, tstop;
   int           ntime, rank;

   /* Define time domain: ntime intervals with dt=1 */
   ntime  = 10;
   tstart = 0.0;
   tstop  = tstart + ntime;
   
   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->rank)   = rank;
   
   /* initialize XBraid and set options */
   braid_Init(MPI_COMM_WORLD, MPI_COMM_WORLD, tstart, tstop, ntime, app,
             my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
             my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
   
   /* Set some typical Braid parameters */
   braid_SetPrintLevel( core, 1);
   braid_SetMaxLevels(core, 2);
   braid_SetAbsTol(core, 1.0e-06);
   braid_SetCFactor(core, -1, 2);
   
   /* Run simulation, and then clean up */
   braid_Drive(core);

   braid_Destroy(core);
   MPI_Finalize();

   return (0);
}
