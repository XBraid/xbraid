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

/* Example BDF2: 1D linear problem with constant coefficient using BDF 2 time integrator

   Compile with: make ex-bdf2

   Help with:    ex-bdf2 -help

   Sample run:   mpirun -np 4 -ntime 10000 -cf 10 -tol 1e-10

   Uses: shellvector features

   Description: This code solves the simple 1D problem

                               y_t=lambda*y, y(0)=1

                where lambda=-1. It serves as a demonstration code
                of a multistep integration in XBraid using the
                shellvector feature to store time values information
                at all points in time. The BDF 2 time stepping is an
                implicit scheme so the solve is hard-coded for this
                specific equation.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.h"

typedef struct _braid_App_struct
{
   double  tstart;
   double  tstop;
   int     ntime;
   FILE   *file;
} my_App;

typedef struct _braid_Vector_struct
{
   // Shell part tprev[0]=t_{n-1}, tprev[1]=t_{n-2}
   double* tprev;

   // Vector part yprev[0]=y_{n-1}, yprev[1]=y_{n-2}
   double* yprev;

} my_Vector;

double lambda=-1.0;
int number_of_access=0;

int
BDF2_Expo1D(double  t0,
            double  t1,
            double  t2,
            double  y1,
            double  y2,
            double* y0)
{
   double h0=t0-t1;
   double rho=h0/(t1-t2);
   double c0=(2*rho+1)/(rho+1)-lambda*h0;
   *y0=(rho+1)/c0*y1-(rho*rho/(rho+1))/c0*y2;
   return 0;
}

int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   double new_y2, new_y1;
   /* Here, we use the fact that ustop is at least a shell containing
      the time values of the targeted time step. */
   BDF2_Expo1D(ustop->tprev[1],u->tprev[0],u->tprev[1],u->yprev[0],u->yprev[1],&new_y2);
   BDF2_Expo1D(ustop->tprev[0],ustop->tprev[1],u->tprev[0],new_y2,u->yprev[0],&new_y1);
   u->tprev[0]=ustop->tprev[0];
   u->tprev[1]=ustop->tprev[1];
   u->yprev[0]=new_y1;
   u->yprev[1]=new_y2;
   
   /* Nonzero rhs */
   if (fstop != NULL)
   {
      u->yprev[0]=u->yprev[0]+fstop->yprev[0];
      u->yprev[1]=u->yprev[1]+fstop->yprev[1];
   }

   /* no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;
   u=(my_Vector *) malloc(sizeof(my_Vector));
   u->tprev=(double*)malloc(2*sizeof(double));
   u->yprev=(double*)malloc(2*sizeof(double));

   double dt = (app->tstop-app->tstart)/app->ntime;
   u->tprev[1]=t;
   u->tprev[0]=t+0.5*dt;
   u->yprev[1]=1.0;
   u->yprev[0]=exp(lambda*0.5*dt);

   *u_ptr=u;
   return 0;
}

/* Initialization of the shell part only (=> tprev, not yprev) */
int
my_InitShell(braid_App     app,
             double        t,
             braid_Vector *u_ptr)
{
   my_Vector *u;
   u=(my_Vector *) malloc(sizeof(my_Vector));
   u->tprev=(double*)malloc(2*sizeof(double));
   u->yprev=NULL;
 
   double dt = (app->tstop-app->tstart)/app->ntime;
   u->tprev[1]=t;
   u->tprev[0]=t+0.5*dt;

   *u_ptr=u;
   return 0;
}

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   my_Vector *v;
   v=(my_Vector *) malloc(sizeof(my_Vector));
   v->tprev=(double*)malloc(2*sizeof(double));
   v->yprev=(double*)malloc(2*sizeof(double));

   v->tprev[0]=u->tprev[0];
   v->tprev[1]=u->tprev[1];
   v->yprev[0]=u->yprev[0];
   v->yprev[1]=u->yprev[1];

   *v_ptr=v;
   return 0;
}

/* Cloning only the shell part of the vector*/
int
my_CloneShell(braid_App     app,
              braid_Vector  u,
              braid_Vector *v_ptr)
{
   my_Vector *v;
   v=(my_Vector *) malloc(sizeof(my_Vector));
   v->tprev=(double*)malloc(2*sizeof(double));
   v->yprev=NULL;

   v->tprev[0]=u->tprev[0];
   v->tprev[1]=u->tprev[1];

   *v_ptr=v;
   return 0;
}

int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u->tprev);
   free(u->yprev);
   u->tprev=NULL;
   u->yprev=NULL;
   free(u);

   return 0;
}

/* Free only the non shell part of the vector, to keep the shell*/
int
my_FreeShell(braid_App    app,
             braid_Vector u)
{
   if (u->yprev != NULL)
      free(u->yprev);
   u->yprev=NULL;

   return 0;
}

int
my_Sum(braid_App     app,
       double        alpha,
       braid_Vector  x,
       double        beta,
       braid_Vector  y)
{
   y->yprev[0]=alpha*x->yprev[0]+beta*y->yprev[0];
   y->yprev[1]=alpha*x->yprev[1]+beta*y->yprev[1];

   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   double sqdot=pow(u->yprev[0],2)+pow(u->yprev[1],2);
   *norm_ptr=sqrt(sqdot);
   return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int done;
   braid_AccessStatusGetDone(astatus,&done);
   if (done==1)
   {
      int idx;
      idx=((u->tprev[1]-app->tstart)/((app->tstop-app->tstart)/(2*app->ntime))+0.1);
      fprintf(app->file,"%d %.14e %.14e\n%d %.14e %.14e\n",idx,u->tprev[1],u->yprev[1],idx+1,u->tprev[0],u->yprev[0]);
      fflush(app->file);
   }
   number_of_access++;
   return 0;
}

int
my_BufSize(braid_App  app,
           int       *size_ptr,
           braid_BufferStatus  bstatus)
{
   *size_ptr=4*sizeof(double);
   return 0;
}

int
my_BufPack(braid_App     app,
           braid_Vector  u,
           void         *buffer,
           braid_BufferStatus  bstatus)
{
   double* dbuffer=(double*)buffer;
   dbuffer[0]=u->tprev[1];
   dbuffer[1]=u->yprev[1];
   dbuffer[2]=u->tprev[0];
   dbuffer[3]=u->yprev[0];
   braid_BufferStatusSetSize(bstatus,4*sizeof(double));
   return 0;
}

int
my_BufUnpack(braid_App           app,
             void                *buffer,
             braid_Vector        *u_ptr,
             braid_BufferStatus  status)
{
   double* dbuffer=(double*)buffer;
   my_Vector *u;
   u=(my_Vector *) malloc(sizeof(my_Vector));
   u->tprev=(double*)malloc(2*sizeof(double));
   u->yprev=(double*)malloc(2*sizeof(double));

   u->tprev[1]=dbuffer[0];
   u->yprev[1]=dbuffer[1];
   u->tprev[0]=dbuffer[2];
   u->yprev[0]=dbuffer[3];

   *u_ptr=u;
   return 0;
}

int main (int argc, char *argv[])
{
   int mpi_rank, mpi_size;
   MPI_Comm comm=MPI_COMM_WORLD;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(comm,&mpi_rank);
   MPI_Comm_size(comm,&mpi_size);

   /* Default values of parameters */
   double tstart,tstop,tol;
   int ntime,max_levels,max_iter,cfactor,access_level,fmg;
   tstart=0;
   tstop=5;
   ntime=1000;
   max_levels=100;
   max_iter=100;
   cfactor=2;
   tol=1e-10;
   access_level=1;
   fmg=1;

   int arg_index=1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         if ( mpi_rank == 0 )
         {
            printf("\n");
            printf("  -ntime <ntime>    : set num time points (default %d)\n", ntime);
            printf("  -tstop <tstop>    : set end time (default %g)\n", tstop);
            printf("  -ml  <max_levels> : set max levels (default %d)\n",max_levels);
            printf("  -tol <tol>        : set stopping tolerance (default %g)\n",tol);
            printf("  -cf  <cfactor>    : set coarsening factor (default %d)\n",cfactor);
            printf("  -mi  <max_iter>   : set max iterations (default %d)\n",max_iter);
            printf("  -ac  <access_lvl> : set the access level (default %d)\n",access_level);
            printf("  -fmg              : use FMG cycling (default %d)\n",fmg);
            printf("\n");
         }
         exit(1);
      }
      else if ( strcmp(argv[arg_index], "-ntime") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tstop") == 0 ) 
      {
         arg_index++;
         tstop = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ml") == 0 )
      {
         arg_index++;
         max_levels = atoi(argv[arg_index++]);
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
      else if ( strcmp(argv[arg_index], "-ac") == 0 )
      {
         arg_index++;
         access_level = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-fmg") == 0 )
      {
         arg_index++;
         fmg = 1;
      }
      else
      {
         arg_index++;
         /*break;*/
      }
   }

   //Because of BDF 2
   ntime=(ntime+1)/2;

   my_App *app=(my_App *) malloc(sizeof(my_App));
   app->tstart=tstart;
   app->tstop=tstop;
   app->ntime=ntime;

   char filename[255];
   sprintf(filename,"ex-bdf2_%03d.out",mpi_rank);
   app->file=fopen(filename,"w");

   braid_Core core;
   braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
             my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
             my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
   braid_SetPrintLevel(core, 1);

   braid_SetMaxLevels(core, max_levels);
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, cfactor);
   braid_SetMaxIter(core, max_iter);
   braid_SetAccessLevel(core,access_level);
   if (fmg)
      braid_SetFMG(core);

   braid_SetShell(core, my_InitShell, my_CloneShell, my_FreeShell);
   braid_SetStorage(core, -1);

   braid_Drive(core);
   braid_Destroy(core);

   int total_number_of_access;
   MPI_Reduce(&number_of_access,&total_number_of_access,1,MPI_INT,MPI_SUM,0,comm);
   if (mpi_rank==0)
   {
      printf("Number of calls to my_Access: %d\n",total_number_of_access);
   }
   fclose(app->file);
   free(app);

   /* Sort output */
   ntime=ntime+1;
   double* time=(double*)malloc(2*ntime*sizeof(double));
   double* solution=(double*)malloc(2*ntime*sizeof(double));
   int nf,idx;
   FILE* file;
   for (nf=0;nf<mpi_size;nf++)
   {
      sprintf(filename,"ex-bdf2_%03d.out",nf);
      file=fopen(filename,"r");
      while (fscanf(file,"%d",&idx)>0)
      {
         fscanf(file,"%le",&time[idx]);
         fscanf(file,"%le",&solution[idx]);
      }
      fclose(file);
      remove(filename);
   }
   file=fopen("ex-bdf2.out","w");
   for (idx=0;idx<2*ntime;idx++)
      fprintf(file,"%.14e %.14e\n",time[idx],solution[idx]);
   fclose(file);
   free(solution);
   free(time);

   MPI_Finalize();
   return 0;
}
