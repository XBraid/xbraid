/*
   Example 01

   Compile with: make drive-01

   Sample run:   mpirun -np 2 drive-01

   Description:

   Blah...

   When run with the default 10 time steps, the solution is as follows:

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

#include "warp.h"

/*--------------------------------------------------------------------------
 * My integration routines
 *--------------------------------------------------------------------------*/

/* can put anything in my app and name it anything as well */
typedef struct _warp_App_struct
{
   MPI_Comm  comm;
   double    tstart;
   double    tstop;
   int       ntime;

} my_App;

/* can put anything in my vector and name it anything as well */
typedef struct _warp_Vector_struct
{
   double value;

} my_Vector;

int
my_Phi(warp_App     app,
       double       tstart,
       double       tstop,
       int          gzero,
       warp_Vector  u,
       int         *rfactor_ptr)
{
   /* On the finest grid, each value is half the previous value */
   (u->value) = pow(0.5, tstop-tstart)*(u->value);

   if (!gzero)
   {
      /* Zero rhs for now */
      (u->value) += 0.0;
   }

   /* no refinement */
   *rfactor_ptr = 1;

   return 0;
}

int
my_Init(warp_App     app,
        double       t,
        warp_Vector *u_ptr)
{
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   if (t == 0.0)
   {
      /* Initial guess */
      (u->value) = 1.0;
   }
   else
   {
      /* Random between 0 and 1 */
      (u->value) = ((double)rand()) / RAND_MAX;
   }
   *u_ptr = u;

   return 0;
}

int
my_Clone(warp_App     app,
         warp_Vector  u,
         warp_Vector *v_ptr)
{
   my_Vector *v;

   v = (my_Vector *) malloc(sizeof(my_Vector));
   (v->value) = (u->value);
   *v_ptr = v;

   return 0;
}

int
my_Free(warp_App    app,
        warp_Vector u)
{
   free(u);

   return 0;
}

int
my_Sum(warp_App    app,
       double      alpha,
       warp_Vector x,
       double      beta,
       warp_Vector y)
{
   (y->value) = alpha*(x->value) + beta*(y->value);

   return 0;
}

int
my_Dot(warp_App     app,
       warp_Vector  u,
       warp_Vector  v,
       double      *dot_ptr)
{
   double dot;

   dot = (u->value)*(v->value);
   *dot_ptr = dot;

   return 0;
}

int
my_Write(warp_App     app,
         double       t,
         warp_Vector  u)
{
   MPI_Comm   comm   = (app->comm);
   double     tstart = (app->tstart);
   double     tstop  = (app->tstop);
   int        ntime  = (app->ntime);
   int        index, myid;
   char       filename[255];
   FILE      *file;

   index = ((t-tstart) / ((tstop-tstart)/ntime) + 0.1);

   MPI_Comm_rank(comm, &myid);

   sprintf(filename, "%s.%07d.%05d", "drive-01.out", index, myid);
   file = fopen(filename, "w");
   fprintf(file, "%.14e\n", (u->value));
   fflush(file);
   fclose(file);

   return 0;
}

int
my_BufSize(warp_App  app,
           int      *size_ptr)
{
   *size_ptr = sizeof(double);
   return 0;
}

int
my_BufPack(warp_App     app,
           warp_Vector  u,
           void        *buffer)
{
   double *dbuffer = buffer;

   dbuffer[0] = (u->value);

   return 0;
}

int
my_BufUnpack(warp_App     app,
             void        *buffer,
             warp_Vector *u_ptr)
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
   warp_Core  core;
   my_App    *app;
   MPI_Comm   comm;
   double     tstart, tstop;
   int        ntime;

   int        max_levels = 1;
   int        nrelax     = 1;
   double     tol        = 1.0e-06;
   int        cfactor    = 2;
   int        max_iter   = 100;

   int        arg_index;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);

   /* ntime time intervals with spacing 1 */
   comm   = MPI_COMM_WORLD;
   ntime  = 10;
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
            printf("  -tol <tol>        : set stopping tolerance\n");
            printf("  -cf  <cfactor>    : set coarsening factor\n");
            printf("  -mi  <max_iter>   : set max iterations\n");
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

   warp_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
             my_Phi, my_Init, my_Clone, my_Free, my_Sum, my_Dot, my_Write,
             my_BufSize, my_BufPack, my_BufUnpack,
             &core);

   warp_SetMaxLevels(core, max_levels);
   warp_SetNRelax(core, nrelax);
   warp_SetTol(core, tol);
   warp_SetCFactor(core, -1, cfactor);
   /*warp_SetCFactor(core,  0, 10);*/
   warp_SetMaxIter(core, max_iter);

   warp_Drive(core);

   warp_PrintStats(core);

   warp_Destroy(core);

   /* Finalize MPI */
   MPI_Finalize();

   return (0);
}
