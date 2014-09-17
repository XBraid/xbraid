<!--
  - Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
  - Produced at the Lawrence Livermore National Laboratory. Written by 
  - Jacob Schroder schroder2@llnl.gov, Rob Falgout falgout2@llnl.gov,
  - Tzanio Kolev kolev1@llnl.gov, Ulrike Yang yang11@llnl.gov, 
  - Veselin Dobrev dobrev1@llnl.gov, et al. 
  - LLNL-CODE-660355. All rights reserved.
  - 
  - This file is part of XBraid. Email schroder2@llnl.gov on how to download. 
  - 
  - This program is free software; you can redistribute it and/or modify it under
  - the terms of the GNU General Public License (as published by the Free Software
  - Foundation) version 2.1 dated February 1999.
  - 
  - This program is distributed in the hope that it will be useful, but WITHOUT ANY
  - WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
  - PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
  - License for more details.
  - 
  - You should have received a copy of the GNU Lesser General Public License along
  - with this program; if not, write to the Free Software Foundation, Inc., 59
  - Temple Place, Suite 330, Boston, MA 02111-1307 USA
 -->

## A Simple Example 

### User Defined Structures and Wrappers

As mentioned, the user must wrap their existing time stepping routine per the XBraid interface. 
To do this, the user must define two data structures and some wrapper routines.  To make the 
idea more concrete, we now give these function definitions from examples/drive-01, which 
implements a scalar ODE, \f$ u_t = \lambda u \f$.

The two data structures are:
1. **App**: This holds a wide variety of information and is *global* in that it
  is passed to every function.  This structure holds everything that the user will need 
  to carry out a simulation.  Here, this is just the global MPI communicator and few values
  describing the temporal domain.
      
      typedef struct _braid_App_struct
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

      typedef struct _braid_Vector_struct
      {
         double value;
      
      } my_Vector;


The user must also define a few wrapper routines.  Note, that the app structure is the 
first argument to every function.
1. **Phi**: This function tells XBraid how to take a time step, and is the core user routine. 
   The user must advance the vector *u* from time *tstart* to time *tstop*.
   Note how the time values are given to the user through the *status*
   structure and associated Get routines.  The *rfactor_ptr* parameter is an
   advanced topic not used here.
   
   Here advancing the solution just involves the scalar \f$ \lambda \f$.  

   **Importantly,** the \f$ g_i \f$ function (from @ref braidoverview) must be 
   incorporated into Phi, so that \f$\Phi(u_i) \rightarrow u_{i+1} \f$

         int
         my_Phi(braid_App       app,
                braid_Vector    u,
                braid_PhiStatus status)
         {
            double tstart;             /* current time */
            double tstop;              /* evolve to this time*/
            braid_PhiStatusGetTstartTstop(status, &tstart, &tstop);
      
            /* On the finest grid, each value is half the previous value */
            (u->value) = pow(0.5, tstop-tstart)*(u->value);
      
            /* Zero rhs for now */
            (u->value) += 0.0;
      
            /* no refinement */
            braid_PhiStatusSetRFactor(status, 1);
      
            return 0;
         }

2. **Init**: This function tells XBraid how to initialize a vector at time *t*.  
   Here that is just allocating and setting a scalar on the heap.

         int
         my_Init(braid_App     app,
                 double     t,
                 braid_Vector *u_ptr)
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


3. **Clone**: This function tells XBraid how to clone a vector 
   into a new vector.

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

4. **Free**: This function tells XBraid how to free 
   a vector.

         int
         my_Free(braid_App    app,
                 braid_Vector u)
         {
            free(u);

            return 0;
         }

5. **Sum**: This function tells XBraid how to sum two 
   vectors (AXPY operation).

         int
         my_Sum(braid_App    app,
                double    alpha,
                braid_Vector x,
                double    beta,
                braid_Vector y)
         {
            (y->value) = alpha*(x->value) + beta*(y->value);

            return 0;
         }

6. **SpatialNorm**: This function tells XBraid how to take the norm 
   of a braid_Vector and is used for halting. This norm is only over 
   space.  A common norm choice is the standard Eucliden norm, but 
   many other choices are possible, such as an L2-norm based on a 
   finite element space.  The norm choice should be based on what 
   makes sense for you problem.  How to accumulate spatial norm values
   to obtain a global space-time residual norm for halting decsions is 
   controlled by [braid_SetTemporalNorm](@ref braid_SetTemporalNorm).

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

7. **Access**: This function allows the user access to XBraid and the current solution vector 
   at time *t*.  This is most commonly used to print solution(s) to screen, file, etc... 
   The user defines what is appropriate output.  Notice how you are told the time value *t* of the 
   vector *u* and even more information in *astatus*.  This lets you tailor the output to only 
   certain time values at certain XBraid iterations.  Querrying *astatus* for such information 
   is done through _braid_AccessStatusGet**(..)_ routines.
   \latexonly \\ \endlatexonly

   The frequency of the calls to *access* is controlled through 
   [braid_SetAccessLevel](@ref braid_SetAccessLevel).  For instance, if access_level is 
   set to 2, then *access* is called every XBraid iteration and on every XBraid level.  In 
   this case, querrying *astatus* to determine the current XBraid level and iteration will 
   be useful. This scenario allows for even more detailed tracking of the simulation.
   \latexonly \\ \endlatexonly

   Eventually, this routine will allow for broader access to XBraid and computational steering.
   \latexonly \\ \endlatexonly
   
   See examples/drive-02 and examples/drive-04 for more advanced uses of the *access* function.  
   Drive-04 uses *access* to write solution vectors to a GLVIS visualization port, and 
   examples/drive-02 uses *access* to write to .vtu files.
   \latexonly \\ \endlatexonly

         int
         my_Access(braid_App          app,
                   braid_Vector       u,
                   braid_AccessStatus astatus)
         {
            MPI_Comm   comm   = (app->comm);
            double     tstart = (app->tstart);
            double     tstop  = (app->tstop);
            int        ntime  = (app->ntime);
            int        index, myid;
            char       filename[255];
            FILE      *file;
            double     t;
            
            braid_AccessStatusGetT(astatus, &t);
            index = ((t-tstart) / ((tstop-tstart)/ntime) + 0.1);
        
            MPI_Comm_rank(comm, &myid);
        
            sprintf(filename, "%s.%07d.%05d", "drive-01.out", index, myid);
            file = fopen(filename, "w");
            fprintf(file, "%.14e\n", (u->value));
            fflush(file);
            fclose(file);
        
            return 0;
         }


8. **BufSize**, **BufPack**, **BufUnpack**: These three routines tell XBraid how to 
   communicate vectors between processors.  *BufPack* packs a vector 
   into a ``void *`` buffer for MPI and then *BufUnPack* unpacks it from ``void *`` 
   to vector.  Here doing that for a scalar is trivial.  *BufSize* computes the 
   upper bound for the size of an arbitrary vector.
   \latexonly \\ \endlatexonly

   Note how *BufPack* also returns a size pointer.  This size pointer should be
   the exact number of bytes packed, while *BufSize* should provide only an
   upper-bound on a possible buffer size.  This flexibility allows for variable
   spatial grid sizes to result in smaller messages sent when appropriate. **To
   avoid MPI issues, it is very important that BufSize be pessimistic, provide 
   an upper bound, and return the same value across processors.**
   \latexonly \\ \endlatexonly

   In general, the buffer should be self-contained.  The receiving processor
   should be able to pull all necessary information from the buffer in order to
   properly interpret and unpack the buffer.
   \latexonly \\ \endlatexonly

         int
         my_BufSize(braid_App  app,
                    int    *size_ptr)
         {
            *size_ptr = sizeof(double);
            return 0;
         }

         int
         my_BufPack(braid_App     app,
                    braid_Vector  u,
                    void         *buffer,
                    braid_Int    *size_ptr)
         {
            double *dbuffer = buffer;
      
            dbuffer[0] = (u->value);
            *size_ptr = sizeof(double);
      
            return 0;
         }

         int
         my_BufUnpack(braid_App     app,
                      void      *buffer,
                      braid_Vector *u_ptr)
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

  These functions allow you vary the spatial mesh size on XBraid levels as depicted here
  where the spatial and temporal grid sizes are halved every level.
  \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.25\textwidth]{../img/spatial_coarsening.pdf}}
       \label{img:heat_results}
   \end{figure}
   \endlatexonly

10. Adaptive and variable time stepping is in the works to be implemented.  The *rfactor* parameter
in *Phi* will allow this.


### Running XBraid

A typical flow of events in the ``main`` function is to first initialize the ``app``
structure.

    /* set up app structure */
    app = (my_App *) malloc(sizeof(my_App));
    (app->comm)   = comm;
    (app->tstart) = tstart;
    (app->tstop)  = tstop;
    (app->ntime)  = ntime;

Then, the data structure definitions and wrapper routines are passed to XBraid.
The core structure is used by XBraid for internal data structures. 

    braid_Core  core;
    braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
            my_Phi, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
            my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

Then, XBraid options are set.

    braid_SetPrintLevel( core, 1);
    braid_SetMaxLevels(core, max_levels);
    braid_SetNRelax(core, -1, nrelax);
    braid_SetAbsTol(core, tol);
    braid_SetCFactor(core, -1, cfactor);
    braid_SetMaxIter(core, max_iter);
   
Then, the simulation is run.

    braid_Drive(core);

Then, we clean up.

    braid_Destroy(core);

Finally, to run drive-01, type

    drive-01 -ml 5

This will run drive-01. See examples/drive-0* for more extensive examples.

### Testing XBraid

The best overall test for XBraid, is to set the maximum number of levels to 1 
(see [braid_SetMaxLevels](@ref braid_SetMaxLevels) ) which will carry out a
sequential time stepping test.  Take the output given to you by your *Access*
function and compare it to output from a non-XBraid run.  Is everything OK?
Once this is complete, repeat for multilevel XBraid, and check that the solution
is correct (that is, it matches a serial run to within tolerance).

At a lower level, to do sanity checks of your data structures and wrapper routines, there are
also XBraid test functions, which can be easily run.  The test routines also 
take as arguments the app structure, spatial communicator comm_x, a stream
like stdout for test output and a time step size to test dt.  After these arguments,
function pointers to wrapper routines are the rest of the arguments. Some of the tests
can return a boolean variable to indicate correctness.
      
    /* Test init(), access(), free() */
    braid_TestInitAccess( app, comm_x, stdout, dt, my_Init, my_Access, my_Free);

    /* Test clone() */
    braid_TestClone( app, comm_x, stdout, dt, my_Init, my_Access, my_Free, my_Clone);

    /* Test sum() */
    braid_TestSum( app, comm_x, stdout, dt, my_Init, my_Access, my_Free, my_Clone, my_Sum);

    /* Test spatialnorm() */
    correct = braid_TestSpatialNorm( app, comm_x, stdout, dt, my_Init, my_Free, my_Clone, 
                              my_Sum, my_SpatialNorm);

    /* Test bufsize(), bufpack(), bufunpack() */
    correct = braid_TestBuf( app, comm_x, stdout, dt, my_Init, my_Free, my_Sum, my_SpatialNorm, 
                        my_BufSize, my_BufPack, my_BufUnpack);
     
    /* Test coarsen and refine */
    correct = braid_TestCoarsenRefine(app, comm_x, stdout, 0.0, dt, 2*dt, my_Init,
                           my_Access, my_Free, my_Clone, my_Sum, my_SpatialNorm, 
                           my_CoarsenInjection, my_Refine);
    correct = braid_TestCoarsenRefine(app, comm_x, stdout, 0.0, dt, 2*dt, my_Init,
                          my_Access, my_Free, my_Clone, my_Sum, my_SpatialNorm, 
                          my_CoarsenBilinear, my_Refine);
