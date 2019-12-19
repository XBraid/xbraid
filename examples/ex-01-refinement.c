/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory. Written by
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
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

/**
 * Example:       ex-01-refinement.c
 *
 * Interface:     C
 *
 * Requires:      only C-language support
 *
 * Compile with:  make ex-01
 *
 * Help with:     this is a simple example showing usage of refinement
 *
 * Sample run:    mpirun -np 2 ex-01
 *
 * Description:   solve the scalar ODE
 *                   u' = lambda u,
 *                   with lambda=-1 and y(0) = 1
 *                in a very simplified XBraid setting.
 *
 *                When run with the default 10 time steps and no refinement, the solution is:
 *                $ ./ex-01
 *                $ cat ex-01.out.00*
 *                  1.00000000000000e+00
 *                  6.66666666666667e-01
 *                  4.44444444444444e-01
 *                  2.96296296296296e-01
 *                  1.97530864197531e-01
 *                  1.31687242798354e-01
 *                  8.77914951989026e-02
 *                  5.85276634659351e-02
 *                  3.90184423106234e-02
 *                  2.60122948737489e-02
 *                  1.73415299158326e-02
 *
 *                The refine option enables the refinement. If set to a value of 1,
 *                an arbitrary refinement is used. If set to a value of 2, the refinement
 *                is based on an estimation of the error. On this simplified problem,
 *                the option -refine 2 leads to a targeted step size of sqrt(2*tol).
 **/

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
   int       limit_rfactor;
   int       refine;
   double    tol;
   int       num_syncs;
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
   double dt;                 /* step size */
   double v;                  /* current value of the solution */
   double LTE;

   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   dt = tstop-tstart;
   v = (u->value);

   /* Use backward Euler to propagate solution */
   /* Use forward Euler to estimate the local trucation error */
   (u->value) = 1./(1. + dt)*v;
   LTE = (u->value) - (1. - dt)*v;
   LTE  = (LTE < 0) ? -LTE : LTE;

   int level, nrefine;
   braid_StepStatusGetLevel(status, &level);
   braid_StepStatusGetNRefine(status, &nrefine);

   /* XBraid only accepts refinements on level 0, and it's also a good idea to
    * cap the number of possible refinements (here capped at 8) */
   if ((level == 0) && (nrefine < 8))
   {
      int rf = 1;
      if (app->refine == 1)
      {
         if (dt>0.001)
         {
            if ( (tstart<=2.5+0.00001)&&(2.5-0.00001<=tstop) )
            {
               rf = 100;
            }
         }
      }
      else if (app->refine == 2)
      {
         rf = (int)(ceil(sqrt(0.5*LTE/(u->value)/(app->tol))));
      }
      else if (app->refine == 3)
      {
         if (dt>0.001)
         {
            if ( (tstart<=2.5+0.00001)&&(2.5-0.00001<=tstop) )
            {
               double newdt, *dtvalues;
               int    i;

               rf = 2;
               newdt = dt / (double) rf;
               dtvalues = (double*) malloc((rf-1)*sizeof(double));
               for (i=0; i<rf-1; i++)
               {
                  dtvalues[i] = newdt;
               }
               braid_StatusSetRefinementDtValues((braid_Status)status, rf, dtvalues);
               free(dtvalues);
            }
         }
      }

      rf = (rf < 1) ? 1 : rf;
      if (app->limit_rfactor > 0)
      {
         rf = (rf < app->limit_rfactor) ? rf : app->limit_rfactor;
      }
      braid_StepStatusSetRFactor(status, rf);
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
   double     t;

   braid_AccessStatusGetT(astatus, &t);
   braid_AccessStatusGetTIndex(astatus, &index);
   sprintf(filename, "%s.%04d.%03d", "ex-01-refinement.out", index, app->rank);
   file = fopen(filename, "w");
   fprintf(file, "%.14e %.14e\n", t, (u->value));
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

/* Sync status */
int
my_Sync(braid_App app,
        braid_SyncStatus status)
{
   braid_Int calling_fcn;
   braid_SyncStatusGetCallingFunction(status, &calling_fcn);

   if(calling_fcn == braid_ASCaller_FRefine_AfterInitHier)
   {
      app->num_syncs += 1;
   }
   return 0;
}

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
   braid_Core    core;
   my_App       *app;
   double        tstart, tstop, tol;
   int           ntime, rank, limit_rfactor, arg_index, print_usage;
   int           refine, output, storage, fmg, sync, incMaxLvl ;

   /* Define time domain: ntime intervals */
   ntime  = 100;
   tstart = 0.0;
   tstop  = 5.0;
   limit_rfactor = -1;
   print_usage = 0;
   tol = 1.0e-6;
   refine = 0;
   output = 1;
   storage = -1;
   fmg = 0;
   sync = 0;
   incMaxLvl = 0;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   /* Parse command line */
   arg_index = 0;
   while( arg_index < argc )
   {
      if( strcmp(argv[arg_index], "-nt") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-max_rfactor") == 0 )
      {
         arg_index++;
         limit_rfactor = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-refine") == 0 )
      {
         arg_index++;
         refine = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-tol") == 0 )
      {
         arg_index++;
         tol = atof(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-no_output") == 0 )
      {
         arg_index++;
         output = 0;
      }
      else if( strcmp(argv[arg_index], "-help") == 0 )
      {
         print_usage = 1;
         break;
      }
      else if ( strcmp(argv[arg_index], "-storage") == 0 )
      {
         arg_index++;
         storage = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-sync") == 0 )
      {
         arg_index++;
         sync = 1;
      }
      else if( strcmp(argv[arg_index], "-incMaxLvl") == 0 )
      {
         arg_index++;
         incMaxLvl = 1;
      }
      else if ( strcmp(argv[arg_index], "-fmg") == 0 )
      {
         arg_index++;
         fmg = 1;
      }
      else
      {
         if(arg_index > 1)
         {
            printf("UNUSED command line paramter %s\n", argv[arg_index]);
         }
         arg_index++;
      }
   }

   if((print_usage) && (rank == 0))
   {
      printf("\n");
      printf("Usage: %s [<options>]\n", argv[0]);
      printf("\n");
      printf(" General XBraid configuration parameters\n");
      printf(" ---------------------------------------\n");
      printf("  -nt  <n>                      : number of time steps (default: 100)\n");
      printf("  -tol <tol>                    : set the stopping tolerance (default: 1e-6)\n");
      printf("  -refine <n>                   : set the type of temporal refinement (default: 0)\n");
      printf("                                : 0 - no refinement\n");
      printf("                                : 1 - arbitrary refinement around t=2.5\n");
      printf("                                : 2 - refinement based on local truncation error\n");
      printf("                                : 3 - arbitrary refinement around t=2.5, specifying the new time-step sizes\n");
      printf("  -max_rfactor <lim>            : limit the refinement factor (default: -1)\n");
      printf("  -fmg                          : use FMG cycling\n");
      printf("  -storage <level>              : full storage on levels >= level\n");
      printf("  -sync                         : enable calls to the sync function\n");
      printf("  -incMaxLvl                    : increase max number of Braid levels after each FRefine\n");
      printf("  -no_output                    : do not save the solution in output files\n");
      printf("  -help                         : print this help and exit\n");
      printf("\n");
   }

   if( print_usage )
   {
      MPI_Finalize();
      return (0);
   }

   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->rank)   = rank;
   (app->limit_rfactor)   = limit_rfactor;
   (app->refine) = refine;
   (app->tol) = tol;
   (app->num_syncs) = 0;

   /* initialize XBraid and set options */
   braid_Init(MPI_COMM_WORLD, MPI_COMM_WORLD, tstart, tstop, ntime, app,
             my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm,
             my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

   /* Set some typical Braid parameters */
   braid_SetPrintLevel( core, 2);
   braid_SetMaxLevels(core, 15);
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, 2);
   braid_SetRefine(core, 1);
   if (fmg)
   {
      braid_SetFMG(core);
   }
   if (storage >= -2)
   {
      braid_SetStorage(core, storage);
   }
   if (!output)
   {
      braid_SetAccessLevel(core, 0);
   }
   if (sync)
   {
      braid_SetSync(core, my_Sync);
   }
   if (incMaxLvl)
   {
      braid_SetIncrMaxLevels(core);
   }

   /* Run simulation, and then clean up */
   braid_Drive(core);

   if (sync && rank == 0)
   {
      printf("  num_syncs             = %d\n\n", (app->num_syncs));
   }

   braid_Destroy(core);
   free(app);
   MPI_Finalize();

   return (0);
}
