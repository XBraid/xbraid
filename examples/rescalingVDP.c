#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.h"
#include "_braid.h"


void vdp_eval(double y, double z, double beta, double* rhs)
{
  rhs[0] = z;
  rhs[1] = - y + beta * (1.0 - y*y) * z ;
}

/*--------------------------------------------------------------------------
 * User-defined routines and structures
 *--------------------------------------------------------------------------*/

/* App structure can contain anything, and be named anything as well */
typedef struct _braid_App_struct
{
   int       rank;
   int       iter;
   int       ntimes;
   double    beta;
   double*   time;  

   double*   statesy;  
   double*   statesz;  
} my_App;

/* Vector structure can contain anything, and be name anything as well */
typedef struct _braid_Vector_struct
{
   double y;
   double z;

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
   double dt;
   int ts;
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   braid_StepStatusGetTIndex(status, &ts);
   dt = tstop - tstart;

   double rhs[2], rhso[2];
   double res_vdp[2];
   double yo, zo;
   yo = u->y;
   zo = u->z;

  //functional iteration for implicit equation
  int maxiter = 500;
  int k;
  for(k = 0; k < maxiter; k++)
  {
     // right hand side of vdp at actual and previous time step
     vdp_eval(u->y, u->z, app->beta, rhs);
     vdp_eval(yo, zo, app->beta, rhso);
     // compute residuum
     res_vdp[0] = u->y - yo - 0.5 * dt * (rhs[0] + rhso[0]);
     res_vdp[1] = u->z - zo - 0.5 * dt * (rhs[1] + rhso[1]);

     /* Fixed-point update */
     u->y = u->y - 0.95 * res_vdp[0];
     u->z = u->z - 0.95 * res_vdp[1];

    //Stopping criterion
    double normres = sqrt(res_vdp[0]*res_vdp[0] + res_vdp[1]*res_vdp[1]);
//     printf("%1.14e\n", normres);
    if (normres < 1e-8)
    {
      break;
    }
  }

  // convergence control
  if (k >= maxiter - 1){
          printf("NO CONVERGENCE! %d  %1.14e  %1.14e\n", ts, tstop, tstart);
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
      (u->y) = 1.30865;
      (u->z) = 0.201522;
   }
   else /* All other time points set to arbitrary value */
   {
      (u->y) = 0.0;
      (u->z) = 0.0;
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
   v->y = u->y;
   v->z = u->z;
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
   (y->y) = alpha*(x->y) + beta*(y->y);
   (y->z) = alpha*(x->z) + beta*(y->z);
   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   double dot;

   dot  = (u->y)*(u->y);
   dot += (u->z)*(u->z);
   *norm_ptr = sqrt(dot);
   return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int        index, level;
   double     t;
   char       filename[255];
   FILE      *file;
   
   braid_AccessStatusGetT(astatus, &t);
   braid_AccessStatusGetTIndex(astatus, &index);
   braid_AccessStatusGetLevel(astatus, &level);

   if( level == 0 ) 
   {
      sprintf(filename, "%s.%05d", "states.out", app->iter);
      file = fopen(filename, "a");
      fprintf(file, "%d %1.14e    %1.14e %1.14e", index, t, u->y, u->z);
      fprintf(file, "\n");
      fflush(file);
      fclose(file);
   }

   /* store in app */
   (app->statesy)[index] = u->y;
   (app->statesz)[index] = u->z;


   return 0;
}

int
my_BufSize(braid_App          app,
           int                *size_ptr,
           braid_BufferStatus bstatus)
{
   *size_ptr = 2*sizeof(double);
   return 0;
}

int
my_BufPack(braid_App          app,
           braid_Vector       u,
           void               *buffer,
           braid_BufferStatus bstatus)
{
   double *dbuffer = buffer;

   dbuffer[0] = (u->y);
   dbuffer[1] = (u->z);
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
   (u->y) = dbuffer[0];
   (u->z) = dbuffer[1];
   *u_ptr = u;

   return 0;
}

int
my_TimeGrid(braid_App         app,       /**< user-defined _braid_App structure */
            braid_Real       *ta,        /**< temporal grid on level 0 (slice per processor) */
            braid_Int        *ilower,    /**< lower time index value for this processor */
            braid_Int        *iupper)    /**< upper time index value for this processor */
{
        for (int ts = 0; ts <= app->ntimes; ts++)
        {
                ta[ts]  = app->time[ts];
        }
        // *ilower = 0;
        // *iupper = app->ntimes-1;
        printf("hi\n");

        return 0;
}
            

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
   braid_Core    core;
   my_App       *app;
   double        tstart, tstop, dt;
   double*      states[2];
   double*      states0[2];
   double*      time;
   double*      time0;
   double       dot, norm, deltat;
   double       rhs_vdp[2];
   int           ntime, rank;
   FILE*        timefile;
   char         filename[255];

   /* Define time domain: ntime intervals */
   ntime  = 40000;
   dt     = 0.005;
   tstart = 0.0;
   tstop  = ntime * dt;
   
   /* Define VDP */
   double beta = 2.0;
   states[0] = (double*) malloc(ntime*sizeof(double));
   states[1] = (double*) malloc(ntime*sizeof(double));
   states0[0] = (double*) malloc(ntime*sizeof(double));
   states0[1] = (double*) malloc(ntime*sizeof(double));
   time0     = (double*) malloc((ntime+1)*sizeof(double));
   time      = (double*) malloc((ntime+1)*sizeof(double));
   for (int ts = 0; ts <= ntime; ts++)
   {
           time[ts]  = ts*dt;
           time0[ts] = ts*dt;
           states0[0][ts] = 0.0;
           states0[1][ts] = 0.0;
   }

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->rank)   = rank;
   (app->beta)   = beta;
   app->statesy  = states[0];
   app->statesz  = states[1];
   app->ntimes   = ntime;
   app->time     = time;
   
   /* initialize XBraid and set options */
   braid_Init(MPI_COMM_WORLD, MPI_COMM_WORLD, tstart, tstop, ntime, app,
             my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
             my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
//    braid_SetTimeGrid(core, my_TimeGrid);
   
   /* Set some typical Braid parameters */
   braid_SetPrintLevel( core, 1);
   braid_SetMaxLevels(core, 5);
   braid_SetAbsTol(core, 1.0e-06);
   braid_SetCFactor(core, -1, 2);
   braid_SetAccessLevel(core, 2);
   braid_SetSkip(core, 0);
   braid_SetMaxIter(core, 1);
   
   
   // /* Run simulation, and then clean up */
   for (int i = 0; i < 13; i++)
   {
        braid_Drive(core);
        // int ts = 40000-1;
        // printf("%d %1.14e\n", ts, app->statesy[ts]);
        
        /* rescaling */
        sprintf(filename, "%s.%05d", "rescaled.out", i);
        timefile = fopen(filename, "a");
        time[0] = tstart;

        double tnorm = 0.0;

        for (int ts = 0; ts < ntime; ts++)
        {
                time0[ts+1] = time[ts+1]; 
        }
        for (int ts = 0; ts < ntime; ts++)
        {
                app->iter = i;
                /* get right hand side */
                vdp_eval(app->statesy[ts], app->statesz[ts], beta, rhs_vdp);

                /* Dotproduct */
                dot  = (app->statesy[ts+1] - app->statesy[ts]) * rhs_vdp[0];
                dot += (app->statesz[ts+1] - app->statesz[ts]) * rhs_vdp[1];

                /* norm */
                norm  = (app->statesz[ts+1] - app->statesz[ts]) * (app->statesz[ts+1] - app->statesz[ts]);
                norm += (app->statesy[ts+1] - app->statesy[ts]) * (app->statesy[ts+1] - app->statesy[ts]);

                /* update */
                deltat = norm / dot;
                time[ts+1] = time[ts] + deltat; 

                fprintf(timefile, "%d %1.14e %1.14e    %1.14e %1.14e", ts, ts*dt, time[ts], 
                        app->statesy[ts], app->statesz[ts]);
                fprintf(timefile, "\n");
                fflush(timefile);

                /* norm */
                
                
                        //double dydt = (app->statesy[ts+1] - app->statesy[ts]) / (time0[ts+1]-time0[ts]) ;
                        //double dzdt = (app->statesz[ts+1] - app->statesz[ts]) / (time0[ts+1]-time0[ts]) ;
                        double dydt = rhs_vdp[0] ;
                        double dzdt = rhs_vdp[1] ;
                        double scaledy = app->statesy[ts+1] - dydt * (deltat);
                        double scaledz = app->statesz[ts+1] - dzdt * (deltat);
                        tnorm += (scaledy - states0[0][ts+1]) * (scaledy - states0[0][ts+1]);
                        tnorm += (scaledz - states0[1][ts+1]) * (scaledz - states0[1][ts+1]);
                
                        // tnorm += (app->statesy[ts] - states0[0][ts]) * (app->statesy[ts] - states0[0][ts]);
                        // tnorm += (app->statesz[ts] - states0[1][ts]) * (app->statesz[ts] - states0[1][ts]);

                        states0[0][ts+1] = app->statesy[ts+1];
                        states0[1][ts+1] = app->statesz[ts+1];
                
        //       printf("%1.14e %1.14e %1.14e\n", scaledy - states0[0][ts], dydt, deltat);
        }
        fclose(timefile);
        tnorm = sqrt(tnorm)/ntime;
        printf("tnorm %1.14e\n", tnorm);

   }

   braid_Destroy(core);
   free(app);
   free(time);     
   free(states[0]);
   free(states[1]);
   MPI_Finalize();

   return (0);
}
