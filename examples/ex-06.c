
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
   Example 07_a

   Compile with: make ex-07

   Sample run:   mpirun -np 2 ex-07

   Description:


*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "braid.h"

/*--------------------------------------------------------------------------
 * Time integration routines
 *--------------------------------------------------------------------------*/

/* Can put anything in my app and name it anything as well */
typedef struct _braid_App_struct
{
   MPI_Comm  comm;
   double    tstart;
   double    tstop;
   int       ntime;
   int       time_discr;    /* Time-discretization options: 1 is Backward Euler, 2 is SDIRK-23 */         
} my_App;

/* Can put anything in my vector and name it anything as well */
typedef struct _braid_Vector_struct
{
   double value0,error;

} my_Vector;

/* Solve k = f(x + dt k,t+dt) */ 
braid_Real 
implicit_solve(braid_Real x,
               braid_Real tstop,
               braid_Real dt)
{
   return ( -4.0*x + 1 - tstop )/( 1.0 + 4.0*dt );
}

/* Take backward Euler time-step */
braid_Real 
myBE(braid_Real x,
     braid_Real tstart,
     braid_Real dt)
{
   return x + dt*implicit_solve(x, tstart + dt , dt);
}

/* Take SDIRK-23 time-step */
braid_Real 
mysdirk23(braid_Real x,
          braid_Real tstart,
          braid_Real dt)
{
    braid_Real gamma, k, y;
    gamma = (2. - sqrt(2.))/2.;
    k = implicit_solve(x, tstart+gamma*dt, gamma*dt);
    y = x + (1.0-2.0*gamma)*dt*k;
    x = x + k*dt/2.0;
    k = implicit_solve(y,  tstart + (1.0-gamma)*dt, gamma*dt);
    x = x + k*(dt/2.0);
    return x; 
}

/* Braid time-stepping routine */
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

   double dt = tstop - tstart;
   
   /* Choose backward Euler, or SDIRK23 */
   if (app->time_discr == 1)
   {
      u->value0 = myBE(u->value0, tstart,dt);
   }
   else if (app->time_discr == 2)
   {
      u->value0 = mysdirk23(u->value0, tstart,dt);
   }
   else
   {
      u->value0 = myBE(u->value0, tstart,dt);
   }

   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   u->error = 0;
   if (t == 0.0)
   {
      /* Initial condition */
      (u->value0) =  1.0;
   }
   else
   {
      /* Initialize all other time points */
      (u->value0) = 0.456;
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
   (v->value0) = (u->value0);
   (v->error) = (u->error);
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
   (y->value0) = alpha*(x->value0) + beta*(y->value0);

   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   double dot;

   dot = (u->value0)*(u->value0);
   *norm_ptr = sqrt(dot);

   return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int        index, iter, nrefine,level,done,caller,nt;
   double     t;
   braid_AccessStatusGetLevel(astatus, &level);
   braid_AccessStatusGetCallingFunction(astatus, &caller);
   if (level > 0 )
   {
      return 0;
   }
   braid_AccessStatusGetDone(astatus,&done);
   braid_AccessStatusGetT(astatus, &t);
   
   braid_Real exact_solution;
   exact_solution = (1.0/16.0)*( -4.0*t + 11.0*exp(-4.0*t) + 5.0 );
   braid_AccessStatusGetNTPoints(astatus,&nt);
   braid_AccessStatusGetIter(astatus, &iter);
   braid_AccessStatusGetNRefine(astatus, &nrefine );
   MPI_Comm_rank( MPI_COMM_WORLD, &index );
   
   if (done && t == app->tstop ) 
      printf("Discretization error (t,sol,err): %.20f %.20f %.20f \n", t, u->value0, fabs( (u->value0)- exact_solution )  );

   return 0;
}

int
my_BufSize(braid_App          app,
           int                *size_ptr,
           braid_BufferStatus bstatus)
{
   *size_ptr = 3*sizeof(double);
   return 0;
}

int
my_BufPack(braid_App          app,
           braid_Vector       u,
           void               *buffer,
           braid_BufferStatus bstatus)
{
   double *dbuffer = buffer;

   dbuffer[0] = (u->value0);
   dbuffer[1] = (u->error);
   braid_BufferStatusSetSize( bstatus, 2*sizeof(double) );

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
   (u->value0) = dbuffer[0];
   (u->error) = dbuffer[1];
   *u_ptr = u;

   return 0;
}
/*
int my_GlobalAccess(braid_App     app,
                    braid_Vector         *u_ptr,
                    braid_GAccessStatus  _gstatus)
{

   double *wfactors, *error_est;
   int i, npoints, *rfactors, iter, nrefine;
   braid_GAccessStatusGetFactors( _gstatus, &rfactors, &wfactors, &error_est, &npoints );
   braid_GAccessStatusGetIter( _gstatus, &iter, &nrefine );
   for ( i = 0; i < npoints; i++)
   {   
      rfactors[i] = 1;
      (u_ptr[i])->error = error_est[i];
   }
   return 0 ;
}
*/


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
   int           max_levels = 100;
   int           nrelax     = 1;
   int           nrelax0    = -1;
   double        tol        = 1.0e-12;
   int           cfactor    = 2;
   int           max_iter   = 1500;
   int           refine_time = 0;
   int           arg_index;
   int           min_coarse = 8;
   int           order = -1;
   int           seq=0;
   int           time_discr = 1;
   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   /* ntime time intervals with spacing 1 */
   comm   = MPI_COMM_WORLD;
   ntime  = 256;
   tstart = 0.0;
   tstop  = 0.5;

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
            printf("  -ml  <max_levels>  : set max levels\n");
            printf("  -nu  <nrelax>      : set num F-C relaxations\n");
            printf("  -nu0 <nrelax>      : set num F-C relaxations on level 0\n");
            printf("  -tol <tol>         : set stopping tolerance\n");
            printf("  -cf  <cfactor>     : set coarsening factor\n");
            printf("  -mi  <max_iter>    : set max iterations\n");
            printf("  -rt                : use temporal refinment\n");
            printf("  -tpts <tpts>       : cutoff for time refinement\n");
            printf("  -nt <ntime>        : number of time points\n");
            printf("  -time_discr <int>  : Time-discretization: 1 is Backward Euler, 2 is SDIRK-23 \n");
  
            printf("e.g. ./ex-07                            --> BE with no Richardson \n");
            printf("e.g. ./ex-07 -time_discr 1              --> BE with no Richardson \n\n");
            printf("e.g. ./ex-07 -time_discr 2              --> SDIRK23 with no Richardson \n");
            printf("e.g. ./ex-07 -time_discr 1 -richardson  --> BE with Richardson \n");
            printf("e.g. ./ex-07 -time_discr 2 -richardson  --> SDIRK23 with Richardson \n");
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
      else if ( strcmp(argv[arg_index], "-seq") == 0 )
      {
         arg_index++;
         seq = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu0") == 0 )
      {
         arg_index++;
         nrelax0 = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-time_discr") == 0 )
      {
         arg_index++;
         time_discr = atoi(argv[arg_index++]);; 
      }
      else if ( strcmp(argv[arg_index], "-richardson") == 0 )
      {
         arg_index++;
         order = 1;
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
   (app->time_discr) = time_discr;

   braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
              my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm,
              my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
   braid_SetPrintLevel( core, 2);
   braid_SetMaxLevels(core, max_levels);
   braid_SetNRelax(core, -1, nrelax);
   braid_SetAccessLevel(core,2);

   int richardson = 0;
   if ( order > 0 )
   {
      richardson++;
      if (time_discr == 1) 
      {
      	order = 2;
      } 
      else if (time_discr == 2)
      {
        order = 3;
      }
      else
      {
         printf("Unrecognized time_discr option:  %d\n", time_discr);
         return -1;
      }
   }
   
   braid_SetErrorEstimation(core, 0, richardson, order);
   if (nrelax0 > -1)
   {
      braid_SetNRelax(core,  -1, nrelax0);
   }
   braid_SetRelTol(core, tol);
   braid_SetCFactor(core, -1, cfactor);
   braid_SetMaxIter(core, max_iter);
   braid_SetMinCoarse(core, min_coarse);
   braid_SetStorage(core, 0);
   braid_SetSkip(core,0); 
   if (refine_time)
   {
      braid_SetRefine( core, 1 );
      braid_SetTPointsCutoff( core, tpts);
   }
   int myid;
   MPI_Comm_rank( comm, &myid );
   braid_Drive(core);

   braid_Destroy(core);
   free( app );
   /* Finalize MPI */
   MPI_Finalize();

   return (0);
}
