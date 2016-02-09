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

   Compile with: make ex-burgers

   Sample run:   mpirun -np 2 ex-burgers

   Description:

   Burger's equation

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

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
   double    xstart;
   double    xstop;
   int       nspace;

} my_App;

/* can put anything in my vector and name it anything as well */
typedef struct _braid_Vector_struct
{
   int     size;
   double *values;

} my_Vector;

int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   int k;
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   double deltaT, fstar_plus, fstar_minus, uk, uk_minus, uk_plus;
   double *u_old;
   double deltaX = (app->xstop - app->xstart) / (u->size - 1.0);
   
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   deltaT = tstop - tstart;
   
   /* copy u */
   u_old = (double *) malloc(sizeof(double)*u->size);
   for(k = 0; k < u->size; k++)
      u_old[k] = u->values[k];

   /* update interior of domain */
   for(k = 1; k < u->size; k++)
   {
      uk = u_old[k];
      uk_minus = u_old[k-1];
      uk_plus = u_old[k+1];

      fstar_plus = 0.5*(uk_plus*uk_plus + uk*uk)    - 0.5*fabs(0.5*(uk + uk_plus))*(uk_plus - uk);
      fstar_minus = 0.5*(uk*uk + uk_minus*uk_minus) - 0.5*fabs(0.5*(uk_minus + uk))*(uk - uk_minus);
      
      u->values[k] = uk - (deltaT/deltaX)*(fstar_plus - fstar_minus);
   }

   /* update left boundary point (Dirichlet 0.0) */
   uk = u_old[0];
   uk_minus = 1.0;
   uk_plus = u_old[1];

   fstar_plus = 0.5*(uk_plus*uk_plus + uk*uk)    - 0.5*fabs(0.5*(2*uk + 2*uk_plus))* (uk_plus - uk);
   fstar_minus = 0.5*(uk*uk + uk_minus*uk_minus) - 0.5*fabs(0.5*(2*uk_minus + 2*uk))*(uk - uk_minus);
   
   u->values[0] = uk - (deltaT/deltaX)*(fstar_plus - fstar_minus);

   
   /* update right boundary point (Dirichlet 1.0) */ 
   uk = u_old[u->size-1];
   uk_minus = u_old[u->size-2];
   uk_plus = 0.0;

   fstar_plus = 0.5*(uk_plus*uk_plus + uk*uk)    - 0.5*fabs(0.5*(uk + uk_plus))*(uk_plus - uk);
   fstar_minus = 0.5*(uk*uk + uk_minus*uk_minus) - 0.5*fabs(0.5*(uk_minus + uk))*(uk - uk_minus);
   
   u->values[u->size] = uk - (deltaT/deltaX)*(fstar_plus - fstar_minus);


   /* Free up */
   free(u_old);


   /* no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}

int
my_Residual(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     r,
            braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);


   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;
   int    i, nspace = (app->nspace);
   double xstart, xstop, x1, x2, x;
   
   u = (my_Vector *) malloc(sizeof(my_Vector));
   (u->size)   = nspace+1;
   (u->values) = (double *) malloc((u->size)*sizeof(double));

   xstart = (app->xstart);
   xstop  = (app->xstop);
   x1     = -1.0;
   x2     =  0.0;

   /* Initial guess */
   for (i = 0; i <= nspace; i++)
   {
      x = xstart + ((double)i/nspace)*(xstop - xstart);
      if (x < x1)
      {
         (u->values)[i] = 1.0;
      }
      else if (x < x2)
      {
         (u->values)[i] = 1.0 - (x - x1);
      }
      else
      {
         (u->values)[i] = 0.0;
      }
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
   int i, size = (u->size);

   v = (my_Vector *) malloc(sizeof(my_Vector));
   (v->size)   = size;
   (v->values) = (double *) malloc(size*sizeof(double));
   for (i = 0; i < size; i++)
   {
      (v->values)[i] = (u->values)[i];
   }
   *v_ptr = v;

   return 0;
}

int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u->values);
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
   int i, size = (y->size);

   for (i = 0; i < size; i++)
   {
      (y->values)[i] = alpha*(x->values)[i] + beta*(y->values)[i];
   }

   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   int    i, size = (u->size);
   double dot = 0.0;

   for (i = 0; i < size; i++)
   {
      dot += (u->values)[i]*(u->values)[i];
   }
   *norm_ptr = sqrt(dot);

   return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   MPI_Comm   comm   = (app->comm);
   double     tstart = (app->tstart);
   double     tstop  = (app->tstop);
   int        ntime  = (app->ntime);
   int        size = (u->size);
   int        i, index, myid;
   char       filename[255];
   FILE      *file;
   double     t;
   
   braid_AccessStatusGetT(astatus, &t);
   index = ((t-tstart) / ((tstop-tstart)/ntime) + 0.1);

   MPI_Comm_rank(comm, &myid);

   sprintf(filename, "%s.%07d.%05d", "ex-burgers.out", index, myid);
   file = fopen(filename, "w");
   for (i = 0; i < size; i++)
   {
      fprintf(file, "%.14e\n", (u->values)[i]);
   }
   fflush(file);
   fclose(file);

   return 0;
}

int
my_BufSize(braid_App  app,
           int       *size_ptr)
{
   int size = (app->nspace)+1;
   *size_ptr = size*sizeof(double) + 1;
   return 0;
}

int
my_BufPack(braid_App     app,
           braid_Vector  u,
           void         *buffer,
           braid_Int    *size_ptr)
{
   double *dbuffer = buffer;
   int i, size = (u->size);
   
   dbuffer[0] = size;
   for (i = 0; i < size; i++)
   {
      dbuffer[i+1] = (u->values)[i];
   }

   *size_ptr = size*sizeof(double);

   return 0;
}

int
my_BufUnpack(braid_App     app,
             void         *buffer,
             braid_Vector *u_ptr)
{
   my_Vector *u;
   double    *dbuffer = buffer;
   int        i, size;

   size = dbuffer[0];
   
   u = (my_Vector *) malloc(sizeof(my_Vector));
   (u->size)   = size;
   (u->values) = (double *) malloc(size*sizeof(double));
   for (i = 0; i < size; i++)
   {
      (u->values)[i] = dbuffer[i+1];
   }
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
   MPI_Comm      comm;
   double        tstart, tstop;
   int           ntime;
   double        xstart, xstop;
   int           nspace;

   int           max_levels = 1;
   int           nrelax     = 1;
   int           nrelax0    = -1;
   double        tol        = 1.0e-06;
   int           cfactor    = 2;
   int           max_iter   = 100;
   int           fmg        = 0;
   int           res        = 0;

   int           arg_index;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);

   /* dt = dx = 0.1 by default */
   comm   = MPI_COMM_WORLD;
   tstart =  0.0;
   tstop  =  2.0;
   ntime  =  20;
   xstart = -2.0;
   xstop  =  1.0;
   nspace =  30;
   
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
            printf("  -fmg              : use FMG cycling\n");
            printf("  -res              : use my residual\n");
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
      else if ( strcmp(argv[arg_index], "-mi") == 0 )
      {
         arg_index++;
         max_iter = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-fmg") == 0 )
      {
         arg_index++;
         fmg = 1;
      }
      else if ( strcmp(argv[arg_index], "-res") == 0 )
      {
         arg_index++;
         res = 1;
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
   (app->xstart) = xstart;
   (app->xstop)  = xstop;
   (app->nspace) = nspace;

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
   if (fmg)
   {
      braid_SetFMG(core);
   }
   if (res)
   {
      braid_SetResidual(core, my_Residual);
   }

   braid_Drive(core);

   braid_Destroy(core);

   /* Finalize MPI */
   MPI_Finalize();

   return (0);
}
