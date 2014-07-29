<!--
 - Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 - Produced at the Lawrence Livermore National Laboratory.
 - This file is part of WARP.  See file COPYRIGHT for details.
 -
 - WARP is free software; you can redistribute it and/or modify it under the
 - terms of the GNU Lesser General Public License (as published by the Free
 - Software Foundation) version 2.1 dated February 1999.
 -->

## A Simple Example 

### User Defined Structures and Wrappers

As mentioned, the user must wrap their existing time stepping routine per the Warp interface. 
To do this, the user must define two data structures and some wrapper routines.  To make the 
idea more concrete, we now give these function definitions from examples/drive-01, which 
implements a scalar ODE, \f$ u_t = \lambda u \f$.

The two data structures are:
1. **App**: This holds a wide variety of information and is *global* in that it
  is passed to every function.  This structure holds everything that the user will need 
  to carry out a simulation.  Here, this is just the global MPI communicator and few values
  describing the temporal domain.
      
      typedef struct _warp_App_struct
      {
         MPI_Comm  comm;
         double    tstart;
         double    tstop;
         int       ntime;
       
      } my_App;


2. **Vector**: this defines (roughly) a state vector at a certain time value.  
  It could also contain any other information related to this vector which is 
  needed to evolve the vector to the next time value, like mesh information.
  Here, the vector is just a scalar double.    

      typedef struct _warp_Vector_struct
      {
         double value;
      
      } my_Vector;


The user must also define a few wrapper routines.  Note, that the app structure is the 
first argument to every function.
1. **Phi**: This function tells Warp how to take a time step, and is the core user routine. 
   The user must advance the vector *u* 
   from time *tstart* to time *tstop*.  Here advancing the solution just involves the scalar \f$ \lambda \f$.  
   The *rfactor_ptr* and *accuracy* parameters are advanced topics not used here.

   **Importantly,** the \f$ g_i \f$ function (from @ref warpoverview) must be 
   incorporated into Phi, so that \f$\Phi(u_i) \rightarrow u_{i+1} \f$

         int
         my_Phi(warp_App     app,
                double       tstart,
                double       tstop,
                double       accuracy,
                warp_Vector  u,
                int         *rfactor_ptr)
         {
            /* On the finest grid, each value is half the previous value */
            (u->value) = pow(0.5, tstop-tstart)*(u->value);

            /* Zero rhs for now */
            (u->value) += 0.0;

            /* no refinement */
            *rfactor_ptr = 1;

            return 0;
         }

2. **Init**: This function tells Warp how to initialize a vector at time *t*.  
   Here that is just allocating and setting a scalar on the heap.

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


3. **Clone**: This function tells Warp how to clone a vector 
   into a new vector.

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

4. **Free**: This function tells Warp how to free 
   a vector.

         int
         my_Free(warp_App    app,
                 warp_Vector u)
         {
            free(u);

            return 0;
         }

5. **Sum**: This function tells Warp how to sum two 
   vectors (AXPY operation).

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

6. **Dot**: This function tells Warp how to take the dot 
   product of two vectors.

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

7. **Write**: This function tells Warp how to write a vector at time *t* to screen, file, etc... 
   The user defines what is appropriate output.  Notice how you are told the time value of the 
   vector *u* and even more information in *status*.  This lets you tailor the output to only 
   certain time values.  
   
   If write_level is 2 (see [warp_SetWriteLevel](@ref warp_SetWriteLevel) ), then 
   *Write* is called every Warp iteration and on every Warp level.  In this case, 
   *status* can be querried using the warp_Get**Status() functions, to determine the 
   current Warp level and iteration.  This allows for even more detailed tracking of the
   simulation. 
   
   See examples/drive-02 and examples/drive-04 for more advanced uses of the Write function.  
   Drive-04 writes to a GLVIS visualization port, and examples/drive-02 writes to .vtu files.

         int
         my_Write(warp_App     app,
                  double       t,
                  warp_Status  status,
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


8. **BufSize**, **BufPack**, **BufUnpack**: These three routines tell Warp how to 
   communicate vectors between processors.  *BufPack* packs a vector 
   into a ``void *`` buffer for MPI and then *BufUnPack* unpacks it from ``void *`` 
   to vector.  Here doing that for a scalar is trivial.  *BufSize* computes the 
   upper bound for the size of an arbitrary vector.

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



9. **Coarsen**, **Restrict** (optional): These are advanced options that allow for coarsening
  in space while you coarsen in time.  This is useful for maintaining stable
  explicit schemes on coarse time scales and is not needed here.  See for instance
  examples/drive-04 and examples/drive-05 which use these routines.

  These functions allow you vary the spatial mesh size on Warp levels as depicted here
  where the spatial and temporal grid sizes are halved every level.
  \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.25\textwidth]{../img/spatial_coarsening.pdf}}
       \label{img:heat_results}
   \end{figure}
   \endlatexonly

10. Adaptive and variable time stepping is in the works to be implemented.  The *rfactor* parameter
in *Phi* will allow this.


### Running Warp

A typical flow of events in the ``main`` function is to first initialize the ``app``
structure.

    /* set up app structure */
    app = (my_App *) malloc(sizeof(my_App));
    (app->comm)   = comm;
    (app->tstart) = tstart;
    (app->tstop)  = tstop;
    (app->ntime)  = ntime;

Then, the data structure definitions and wrapper routines are passed to Warp.
The core structure is used by Warp for internal data structures. 

    warp_Core  core;
    warp_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
              my_Phi, my_Init, my_Clone, my_Free, my_Sum, my_Dot, my_Write,
              my_BufSize, my_BufPack, my_BufUnpack,
              &core);
    
Then, Warp options are set.

    warp_SetPrintLevel( core, 1);
    warp_SetMaxLevels(core, max_levels);
    warp_SetNRelax(core, -1, nrelax);
    warp_SetAbsTol(core, tol);
    warp_SetCFactor(core, -1, cfactor);
    warp_SetMaxIter(core, max_iter);
   
Then, the simulation is run.

    warp_Drive(core);

Then, we clean up.

    warp_Destroy(core);

Finally, to run drive-01, type

    drive-01 -ml 5

This will run drive-01. See examples/drive-0* for more extensive examples.

### Testing Warp

The best overall test for Warp, is to set the maximum number of levels to 1 
(see [warp_SetMaxLevels](@ref warp_SetMaxLevels) ) which will carry out a
sequential time stepping test.  Take the output given to you by your *Write*
function and compare it to output from a non-Warp run.  Is everything OK?
Once this is complete, repeat for multilevel Warp, and check that the solution
is correct (that is, it matches a serial run to within tolerance).

At a lower level, to do sanity checks of your data structures and wrapper routines, there are
also Warp test functions, which can be easily run.  The test routines also 
take as arguments the app structure, spatial communicator comm_x, a stream
like stdout for test output and a time step size to test dt.  After these arguments,
function pointers to wrapper routines are the rest of the arguments. Some of the tests
can return a boolean variable to indicate correctness.
      
    /* Test init(), write(), free() */
    warp_TestInitWrite( app, comm_x, stdout, dt, my_Init, my_Write, my_Free);

    /* Test clone() */
    warp_TestClone( app, comm_x, stdout, dt, my_Init, my_Write, my_Free, my_Clone);

    /* Test sum() */
    warp_TestSum( app, comm_x, stdout, dt, my_Init, my_Write, my_Free, my_Clone, my_Sum);

    /* Test dot() */
    correct = warp_TestDot( app, comm_x, stdout, dt, my_Init, my_Free, my_Clone, my_Sum, my_Dot);

    /* Test bufsize(), bufpack(), bufunpack() */
    correct = warp_TestBuf( app, comm_x, stdout, dt, my_Init, my_Free, my_Sum, my_Dot, my_BufSize, my_BufPack, my_BufUnpack);
     
    /* Test coarsen and refine */
    correct = warp_TestCoarsenRefine(app, comm_x, stdout, 0.0, dt, 2*dt, my_Init,
                           my_Write, my_Free, my_Clone, my_Sum, my_Dot, my_CoarsenInjection, 
                           my_Refine);
    correct = warp_TestCoarsenRefine(app, comm_x, stdout, 0.0, dt, 2*dt, my_Init,
                          my_Write, my_Free, my_Clone, my_Sum, my_Dot, my_CoarsenBilinear, 
                          my_Refine);

