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
   Example 01_a

   Compile with: make ex-03

   Sample run:   mpirun -np 2 ex-03

   Description:

   Simple linear example. Solution is u = 3.0*t 

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "braid.h"

/*--------------------------------------------------------------------------
 * My integration routines
 *--------------------------------------------------------------------------*/

/* can put anything in my app and name it anything as well */
typedef struct _braid_App_struct
{
   MPI_Comm  comm;
   double    tstart;
   double    tstop;
   int       ntime;

} my_App;

/* can put anything in my vector and name it anything as well */
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

   /* Solution is a linear line with slope 3 */
   (u->value) = (u->value) + 3.0*(tstop - tstart) ;
  
   double wfactor = ( rand()/( (double) RAND_MAX) )*100 + 1;
   int rfactor = 2;// ( rand() % 10 ) + 1; 
   braid_StepStatusSetWFactor(status, wfactor );
   braid_StepStatusSetRFactor(status, rfactor );
   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   if (t == 0.0)
   {
      /* Initial condition */
      (u->value) = 0.0;
   }
   else
   {
      /* Initialize all other time points */
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
   double     tstart = (app->tstart);
   double     tstop  = (app->tstop);
   int        ntime  = (app->ntime);
   int        index;
   double     t;
   
   braid_AccessStatusGetT(astatus, &t);
   index = ((t-tstart) / ((tstop-tstart)/ntime) + 0.1);

   MPI_Comm_rank( MPI_COMM_WORLD, &index );
  // if ( fabs( 3.0*t - (u->value) ) > 1e-10 )
    //  printf( " %d :  %.14f  \n ",index ,  3.0*t - (u->value));

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
#include <execinfo.h>

int main (int argc, char *argv[])
{
   braid_Core    core;
   my_App       *app;
   MPI_Comm      comm;
   double        tstart, tstop;
   int           ntime;
   int           tpts       = 10000;
   int           max_levels = 1;
   int           nrelax     = 1;
   int           nrelax0    = -1;
   double        tol        = 1.0e-06;
   int           cfactor    = 2;
   int           max_iter   = 100;
   int           lbalance   = 0;
   int           refine_time = 0; 
   int           arg_index;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   /* ntime time intervals with spacing 1 */
   comm   = MPI_COMM_WORLD;
   ntime  = 16;
   tstart = 0.0;
   tstop  = tstart + ntime;
   
   /* Parse command line */

   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         int  myid;
         MPI_Comm_rank(comm, &myid);
         if ( myid == 0 )
         {
            printf("\n");
            printf("  -ml  <max_levels> : set max levels\n");
            printf("  -nu  <nrelax>     : set num F-C relaxations\n");
            printf("  -nu0 <nrelax>     : set num F-C relaxations on level 0\n");
            printf("  -tol <tol>        : set stopping tolerance\n");
            printf("  -cf  <cfactor>    : set coarsening factor\n");
            printf("  -mi  <max_iter>   : set max iterations\n");
            printf("  -lb               : use load balencing\n");
            printf("  -rt               : use temporal refinment\n");
            printf("  -tpts <tpts>      : cutoff for time refinement\n");
            printf("  -nt <ntime>       : number of time points\n");

            printf("\n");
         }
         exit(1);
      }
      else if ( strcmp(argv[arg_index], "-ml") == 0 )
      {
         arg_index++;
         max_levels = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu") == 0 )
      {
         arg_index++;
         nrelax = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu0") == 0 )
      {
         arg_index++;
         nrelax0 = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tol") == 0 ) 
      {
         arg_index++;
         tol = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-cf") == 0 )
      {
         arg_index++;
         cfactor = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tpts") == 0 )
      {
         arg_index++;
         tpts = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nt") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mi") == 0 )
      {
         arg_index++;
         max_iter = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-lb") == 0 )
      {
         arg_index++;
         lbalance = 1;
      }
      else if ( strcmp(argv[arg_index], "-rt") == 0 )
      {
         arg_index++;
         refine_time = 1;
      }
      else
      {
         arg_index++;
         /*break;*/
      }
   }
   
   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->comm)   = comm;
   (app->tstart) = tstart;
   (app->tstop)  = tstop;
   (app->ntime)  = ntime;

   braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
             my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
             my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

   braid_SetPrintLevel( core, 1);
   braid_SetMaxLevels(core, max_levels);
   braid_SetNRelax(core, -1, nrelax);
   if (nrelax0 > -1)
   {
      braid_SetNRelax(core,  0, nrelax0);
   }
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, cfactor);
   /*braid_SetCFactor(core,  0, 10);*/
   braid_SetMaxIter(core, max_iter);
   
   if (lbalance)
   {
      braid_SetLoadBalance( core, lbalance);
   }  
   if (refine_time)
   {
      braid_SetRefine( core, 1 );
      braid_SetTPointsCutoff( core, tpts); 
   }
   
   int myid;
   MPI_Comm_rank( comm, &myid );
   srand( (unsigned int) time(NULL)*( myid + 1 ) );
   braid_Drive(core);

   braid_Destroy(core);
   free( app );
   /* Finalize MPI */
   MPI_Finalize();

   return (0);
}
