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
 * Example:       ex-02.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-02
 *
 * Help with:     ex-02 -help
 *
 * Sample run:    mpirun -np 2 ex-02
 *
 * Description:   Solves the 1D heat equation, with only parallelism in time
 *                Space-time domain:  [0, PI] x [0, 2*PI]
 *                Exact solution:     u(t,x) = sin(x)*cos(t)
 *
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "braid.h"
#include "braid_test.h"
#include "ex-02-lib.c"

/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/

/* can put anything in my app and name it anything as well */
typedef struct _braid_App_struct
{
   MPI_Comm  comm;
   double    tstart;       /* Define the temporal domain */
   double    tstop;
   int       ntime;
   double    xstart;       /* Define the spatial domain */
   double    xstop;
   int       nspace;
   double    matrix[3];    /* the three point spatial discretization stencil */
   double *  g;            /* temporary vector for inversions and mat-vecs */
   double *  sc_info;      /* Runtime information that tracks the space-time grids visited */
   int       print_level;  /* Level of output desired by user (see the -help message below) */
   int       refine_time;  /* Boolean, controls whether adaptive refinement in time is turned on */
   double    error_tol;    /* If doing adaptive refinement, use this absolute step-wise error tol */ 
} my_App;

/* Can put anything in my vector and name it anything as well */
typedef struct _braid_Vector_struct
{
   int     size;
   double *values;

} my_Vector;

/* create and allocate a vector */
void
create_vector(my_Vector **u, 
              int size)
{
   (*u) = (my_Vector *) malloc(sizeof(my_Vector));
   ((*u)->size)   = size;
   ((*u)->values) = (double *) malloc(size*sizeof(double));
}


/*--------------------------------------------------------------------------
 * My integration routines
 *--------------------------------------------------------------------------*/

int my_Step(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     fstop,
            braid_Vector     u,
            braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   int level, i;
   double deltaX, deltaT, error_est;

   braid_StepStatusGetLevel(status, &level);
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   deltaT = tstop - tstart;
   deltaX = (app->xstop - app->xstart) / (ustop->size - 1.0);
   
   /* XBraid forcing */
   if(fstop != NULL)
   {
      for(i = 0; i < u->size; i++)
      {
         u->values[i] = u->values[i] + fstop->values[i];
      }
   }
   
   /* Take backward Euler step 
    * Note: if an iterative solver were used, ustop->values would 
    *       contain the XBraid's best initial guess. */
   take_step(u->values, u->size, tstop, app->xstart, deltaX, deltaT,
         app->matrix, app->g); 

   /* Store info on space-time grids visited during the simulation */
   (app->sc_info)[ (2*level) ] = deltaX;
   (app->sc_info)[ (2*level) + 1] = deltaT;
   

   /* Carry out any adaptive refinement in time 
    *
    * Note that the Richardson estimates are only computed after about one
    * iteration.  If the estimate is not available, error_est == -1.0 
    */
   if( (app->refine_time) && (level == 0) )
   {
      braid_StepStatusGetSingleErrorEstStep(status, &error_est); 
      if(error_est != -1.0)
      {
         int rfactor = (int) ceil( sqrt( error_est / app->error_tol) );
         braid_StepStatusSetRFactor(status, rfactor);
      }
      else
      {
         // Do no refinement in time
         braid_StepStatusSetRFactor(status, 1);
      }
   }


   return 0;
}


int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;
   int    i; 
   int nspace = (app->nspace);
   double deltaX = (app->xstop - app->xstart) / (nspace - 1.0);

   /* Allocate vector */
   create_vector(&u, nspace);
   
   /* Initialize vector (with correct boundary conditions) */
   if(t == 0.0)
   {
      /* Get the solution at time t=0 */
      get_solution(u->values, u->size, 0.0, app->xstart, deltaX);
   }
   else
   {
      /* Use random values for u(t>0), this measures asymptotic convergence rate */
      for(i=0; i < nspace; i++)
      {
         (u->values)[i] = ((double)braid_Rand())/braid_RAND_MAX;
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
   int size = (u->size);
   int i;

   create_vector(&v, size);
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
   int i;
   int size = (y->size);

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
   int    i;
   int size   = (u->size);
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
   int        index, rank, level, done, ntime;
   char       filename[255];
   double     t, error;
   
   braid_AccessStatusGetT(astatus, &t);
   braid_AccessStatusGetTIndex(astatus, &index);
   braid_AccessStatusGetLevel(astatus, &level);
   braid_AccessStatusGetDone(astatus, &done);
   braid_AccessStatusGetNTPoints(astatus, &ntime);
   
   /* Print solution to file if simulation is over */
   if(done)
   {
      MPI_Comm_rank( (app->comm), &rank);
      sprintf(filename, "%s.%07d.%05d", "ex-02.out", index, rank);
      save_solution(filename, u->values, u->size, app->xstart, 
            app->xstop, ntime, app->tstart, app->tstop, t);
   }

   /* IF on the finest level AND print_level is high enough AND at the final time,
    * THEN print out the discretization error */
   if( (level == 0) && ((app->print_level) > 0) && (index == ntime) )
   {
      error = compute_error_norm(u->values, app->xstart, app->xstop, u->size, t);
      printf("  Discretization error at final time:  %1.4e\n", error);
      fflush(stdout);
   }

   return 0;
}

int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
   int size = (app->nspace);
   *size_ptr = (size+1)*sizeof(double);
   return 0;
}

int
my_BufPack(braid_App           app,
           braid_Vector        u,
           void               *buffer,
           braid_BufferStatus  bstatus)
{
   double *dbuffer = buffer;
   int i, size = (u->size);
   
   dbuffer[0] = size;
   for (i = 0; i < size; i++)
   {
      dbuffer[i+1] = (u->values)[i];
   }

   braid_BufferStatusSetSize( bstatus,  (size+1)*sizeof(double));

   return 0;
}

int
my_BufUnpack(braid_App           app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus  bstatus)
{
   my_Vector *u = NULL;
   double    *dbuffer = buffer;
   int        i, size;

   size = dbuffer[0];
   create_vector(&u, size);

   for (i = 0; i < size; i++)
   {
      (u->values)[i] = dbuffer[i+1];
   }
   *u_ptr = u;

   return 0;
}


/*--------------------------------------------------------------------------
 * my_Residual, my_Coarsen and my_Refine are advanced XBraid options, ignore
 * them until you understand the rest of the driver.
 *--------------------------------------------------------------------------*/

int
my_Residual(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     r,
            braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   int i;
   double x, deltaX, deltaT;

   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   deltaT = tstop - tstart;
   deltaX = (app->xstop - app->xstart) / (ustop->size - 1.0);
   
   /* Set up matrix stencil for 1D heat equation*/
   compute_stencil(deltaX, deltaT, app->matrix);
   
   /* Residual r = A*xstop - r - forcing - boundary 
    *   note: there are no boundary terms here */
   matvec_tridiag(ustop->values, app->g, ustop->size, app->matrix);
   x = app->xstart;
   for(i = 0; i < r->size; i++)
   {
      r->values[i] = app->g[i] - r->values[i] - deltaT*forcing(tstop, x);
      x = x + deltaX;
   }

   return 0;
}

/* Bilinear Coarsening */
int
my_Coarsen(braid_App              app,           
           braid_Vector           fu,
           braid_Vector          *cu_ptr,
           braid_CoarsenRefStatus status)
{

   int csize, level;
   my_Vector *v;
   
   /* This returns the level for fu */
   braid_CoarsenRefStatusGetLevel(status, &level);
   
   /* Stop spatial coarsening/interpolation after reaching a grid of 3 points.
    * This is the smallest spatial grid (2 boundary points, one true DOF). */
   if( level < floor(log2(app->nspace)) - 1 )
   {
      csize = ((fu->size) - 1)/2 + 1;
      create_vector(&v, csize);
      coarsen_1D(v->values, fu->values, csize, fu->size);
   }
   else
   {
      /* No coarsening, clone the vector */
      my_Clone(app, fu, &v);
   }   

   *cu_ptr = v;

   return 0;
}

/* Bilinear interpolation */
int
my_Interp(braid_App              app,           
          braid_Vector           cu,
          braid_Vector          *fu_ptr,
          braid_CoarsenRefStatus status)
{

   int fsize, level;
   my_Vector *v;

   /* This returns the level for fu_ptr */
   braid_CoarsenRefStatusGetLevel(status, &level);
   
   /* Stop spatial coarsening/interpolation after reaching a grid of 3 points.
    * This is the smallest spatial grid (2 boundary points, one true DOF). */
   if( level < floor(log2(app->nspace)) - 1 )
   {
      fsize = (cu->size - 1)*2 + 1;
      create_vector(&v, fsize);
      interpolate_1D(cu->values, v->values, cu->size, fsize);     
   }
   else
   {
      /* No refinement, clone the vector */
      my_Clone(app, cu, &v); 
   }

   *fu_ptr = v;
   
   return 0;
}


/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
   braid_Core    core;
   my_App       *app;
   MPI_Comm      comm, comm_x, comm_t;
   int           i, rank, arg_index;
   double        loglevels;
   
   /* Define space-time domain */
   double    tstart        =  0.0;
   double    tstop         =  2*PI;
   int       ntime         =  64;
   double    xstart        =  0.0;
   double    xstop         =  PI;
   int       nspace        =  33;

   /* Define XBraid parameters 
    * See -help message for descriptions */
   int       max_levels    = 2;
   int       nrelax        = 1;
   int       skip          = 0;
   double    CWt           = 1.0;
   double    tol           = 1.0e-07;
   int       cfactor       = 2;
   int       max_iter      = 30;
   int       min_coarse    = 3;
   int       fmg           = 0;
   int       scoarsen      = 0;
   int       res           = 0;
   int       wrapper_tests = 0;
   int       print_level   = 2;
   int       access_level  = 1;
   int       use_sequential= 0;
   int       richardson    = 0;
   int       refine_time   = 0;
   double    error_tol     = 1e-3;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   comm   = MPI_COMM_WORLD;
   MPI_Comm_rank(comm, &rank);
   
   /* Parse command line */
   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         if ( rank == 0 )
         {
            printf("\n");
            printf(" Solve the 1D heat equation on space-time domain:  [0, PI] x [0, 2*PI]\n");
            printf(" with exact solution u(t,x) = sin(x)*cos(t) \n\n");
            printf("   -wrapper_tests       : run the user code + XBraid wrapper tests (no simulation)\n\n");
            printf("   -use_seq             : use the solution from sequential time stepping as the initial guess\n");
            printf("                          for XBraid. All zero residuals should be produced.\n");
            printf("   -ntime <ntime>       : set num points in time\n");
            printf("   -nspace <nspace>     : set num points in space\n\n");
            printf("   -ml   <max_levels>   : set max levels\n");
            printf("   -nu   <nrelax>       : set num F-C relaxations\n");
            printf("   -CWt  <CWt>          : set C-relaxation weight on all levels\n");
            printf("   -skip <set_skip>     : set skip relaxations on first down-cycle; 0: no skip;  1: skip\n");
            printf("   -tol  <tol>          : set stopping tolerance (scaled later by sqrt(dt) sqrt(dx))\n");
            printf("   -cf   <cfactor>      : set coarsening factor on all levels \n");
            printf("   -mc   <min_coarse>   : set min possible coarse level size (default: 3)\n");
            printf("   -mi   <max_iter>     : set max iterations\n");
            printf("   -fmg                 : use FMG cycling\n");
            printf("   -sc                  : use spatial coarsening by factor of 2 each level\n");
            printf("   -richardson          : use built-in Richardson extrapolation to improve finest grid accuracy\n");
            printf("   -refinet <tol>       : use temporal refinment with a Richardson-based error estimate and <tol> absolute tolerance\n");
            printf("                          Try: ./ex-02 -ntime 8 -refinet 1e-2\n");
            printf("                               ./python viz-ex-02.py\n");
            printf("   -res                 : use my residual\n\n");
            printf("   -print_level <l>     : sets the print_level (default: 1) \n");
            printf("                          0 - no output to standard out \n");
            printf("                          2 - Basic convergence information and hierarchy statistics\n");
            printf("                          3 - Debug level output \n");
            printf("   -access_level <l>    : sets the access_level (default: 1) \n");
            printf("                          0 - never call access \n");
            printf("                          1 - call access only after completion \n");
            printf("                          2 - call access at the end of every cycle \n");
            printf("                          3 - call access every iteration and level\n\n");
         }
         exit(1);
      }
      else if ( strcmp(argv[arg_index], "-wrapper_tests") == 0 )
      {
         arg_index++;
         wrapper_tests = 1;
      }
      else if ( strcmp(argv[arg_index], "-use_seq") == 0 )
      {
         arg_index++;
         use_sequential = 1;
         tol = -1.0;
         max_iter = 5;
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
      else if ( strcmp(argv[arg_index], "-CWt") == 0 )
      {
         arg_index++;
         CWt = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-skip") == 0 )
      {
         arg_index++;
         skip = atoi(argv[arg_index++]);
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
      else if( strcmp(argv[arg_index], "-mc") == 0 ){
          arg_index++;
          min_coarse = atoi(argv[arg_index++]);
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
      else if ( strcmp(argv[arg_index], "-sc") == 0 )
      {
         arg_index++;
         scoarsen  = 1;
      }
      else if ( strcmp(argv[arg_index], "-res") == 0 )
      {
         arg_index++;
         res = 1;
      }
      else if ( strcmp(argv[arg_index], "-richardson") == 0 )
      {
         arg_index++;
         richardson = 1;
      }
      else if ( strcmp(argv[arg_index], "-refinet") == 0 )
      {
         arg_index++;
         refine_time = 1;
         error_tol = atof(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-print_level") == 0 ){
         arg_index++;
         print_level = atoi(argv[arg_index++]);
      }
      else if( strcmp(argv[arg_index], "-access_level") == 0 ){
         arg_index++;
         access_level = atoi(argv[arg_index++]);
      }
      else
      {
         printf("ABORTING: incorrect command line parameter %s\n", argv[arg_index]);
         MPI_Finalize();
         return (0);

      }
   }

   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->g)             = (double*) malloc( nspace*sizeof(double) );
   (app->comm)          = comm;
   (app->tstart)        = tstart;
   (app->tstop)         = tstop;
   (app->ntime)         = ntime;
   (app->xstart)        = xstart;
   (app->xstop)         = xstop;
   (app->nspace)        = nspace;
   (app->print_level)   = print_level;
   (app->refine_time)   = refine_time;
   (app->error_tol)     = error_tol;

   /* Initialize storage for sc_info, for tracking space-time grids visited during the simulation */
   app->sc_info = (double*) malloc( 2*max_levels*sizeof(double) );
   for( i = 0; i < 2*max_levels; i++) {
      app->sc_info[i] = -1.0;
   }
   
   /* Initialize Braid */
   braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
          my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
          my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
   
   /* The first step before running simulations, is always to verify the wrapper tests */
   if(wrapper_tests)
   {
      /* Create spatial communicator for wrapper-tests */
      braid_SplitCommworld(&comm, 1, &comm_x, &comm_t);
      
      braid_TestAll(app, comm_x, stdout, 0.0, (tstop-tstart)/ntime, 
                    2*(tstop-tstart)/ntime, my_Init, my_Free, my_Clone, 
                    my_Sum, my_SpatialNorm, my_BufSize, my_BufPack, 
                    my_BufUnpack, my_Coarsen, my_Interp, my_Residual, my_Step);
   }
   else
   {
      /* Scale tol by domain */
      tol = tol/( sqrt((tstop - tstart)/(ntime-1))*sqrt((xstop - xstart)/(nspace-1)) );
      
      /* Set Braid options */
      braid_SetPrintLevel( core, print_level);
      braid_SetAccessLevel( core, access_level);
      braid_SetMaxLevels(core, max_levels);
      braid_SetMinCoarse( core, min_coarse );
      braid_SetSkip(core, skip);
      braid_SetNRelax(core, -1, nrelax);
      braid_SetCRelaxWt(core, -1, CWt);
      braid_SetAbsTol(core, tol);
      braid_SetCFactor(core, -1, cfactor);
      braid_SetMaxIter(core, max_iter);
      braid_SetSeqSoln(core, use_sequential);
      if (fmg)
      {
         braid_SetFMG(core);
      }
      if (res)
      {
         braid_SetResidual(core, my_Residual);
      }
      loglevels = log2(nspace - 1.0);
      if ( scoarsen && ( fabs(loglevels - round(loglevels)) > 1e-10 ))
      {
         if(rank == 0)
         {
            fprintf(stderr, "\nWarning!\nFor spatial coarsening, spatial grids must be a "
                    "power of 2 + 1, \ni.e., nspace = 2^k + 1.  Your spatial grid is of size"
                    " %d.  Spatial \ncoarsening is therefore being ignored.\n\n", nspace);
         }
      }
      else if (scoarsen)
      {
         braid_SetSpatialCoarsen(core, my_Coarsen);
         braid_SetSpatialRefine(core,  my_Interp);
      }
      
      /* Set Richardson and temporal refinement options */
      if (app->refine_time)
      {
         braid_SetRefine( core, 1 );
         braid_SetTPointsCutoff( core, 10000);
      }
      braid_SetRichardsonEstimation(core, app->refine_time, richardson, 2);

      braid_Drive(core);
      
      /* Print accumulated info on space-time grids visited during the simulation */
      if( (print_level > 0) && (rank == 0))
      {
         print_sc_info(app->sc_info, max_levels);
      }
   }

   /* Clean up */
   braid_Destroy(core);
   free( app->sc_info);
   free( app->g);
   free( app );

   /* Finalize MPI */
   MPI_Finalize();

   return (0);
}
