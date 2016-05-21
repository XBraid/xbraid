<!--
  - Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
  - Produced at the Lawrence Livermore National Laboratory. Written by 
  - Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
  - Dobrev, et al. LLNL-CODE-660355. All rights reserved.
  - 
  - This file is part of XBraid. Email xbraid-support@llnl.gov for support.
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

# The Simplest Example {#exampleone}

## User Defined Structures and Wrappers

The user must wrap their existing time stepping routine per the XBraid interface. 
To do this, the user must define two data structures and some wrapper routines.  To make the 
idea more concrete, we now give these function definitions from ``examples/ex-01``, which 
implements a scalar ODE, \f[ u_t = \lambda u. \f]

The two data structures are:

1. **App**: This holds a wide variety of information and is *global* in that it
   is passed to every function.  This structure holds everything that the user
   will need to carry out a simulation.  Here, this is just the global MPI
   communicator and few values describing the temporal domain.
   
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

The user must also define a few wrapper routines.  Note, that the *app* structure is the 
first argument to every function.

1. **Step**: This function tells XBraid how to take a time step, and is the core
   user routine.  The user must advance the vector *u* from time *tstart* to
   time *tstop*.  Note how the time values are given to the user through the
   *status* structure and associated *Get* routine.  The *ustop* parameter
   serves as an approximation to the solution at time *tstop* and is not needed
   here.  It can be useful for implicit schemes that require an initial guess
   for a linear or nonlinear solver.  The *fstop* parameter is the right-hand
   side of the nonlinear problem on the given time grid.  It is only needed when
   providing a residual with [braid_SetResidual](@ref braid_SetResidual).  More
   information on how to use this optional feature is given below.  Also see the
   full code listing for this example in ``examples/ex-01``.
   \latexonly \\ \endlatexonly
   
   Here advancing the solution just involves the scalar \f$ \lambda \f$.
   **Important note:** the \f$ g_i \f$ function from @ref braidoverview must be
   incorporated into *Step*.
         
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
         
            /* On the finest grid, each value is half the previous value */
            (u->value) = pow(0.5, tstop-tstart)*(u->value);
         
            if (fstop != NULL)
            {
               /* Nonzero rhs */
               (u->value) += (fstop->value);
            }
         
            /* no refinement */
            braid_StepStatusSetRFactor(status, 1);
         
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
               /* Initial condition */
               (u->value) = 1.0;
            }
            else
            {
               /* Initialize all other time points */
               (u->value) = 0.456;
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

4. **Free**: This function tells XBraid how to free a vector.

         int
         my_Free(braid_App    app,
                 braid_Vector u)
         {
            free(u);

            return 0;
         }

5. **Sum**: This function tells XBraid how to sum two vectors (AXPY operation).

         int
         my_Sum(braid_App    app,
                double       alpha,
                braid_Vector x,
                double       beta,
                braid_Vector y)
         {
            (y->value) = alpha*(x->value) + beta*(y->value);

            return 0;
         }

6. **SpatialNorm**: This function tells XBraid how to take the norm 
   of a *braid_Vector* and is used for halting. This norm is only over 
   space.  A common norm choice is the standard Euclidean norm, but 
   many other choices are possible, such as an L2-norm based on a 
   finite element space.  The norm choice should be based on what 
   makes sense for your problem.  How to accumulate spatial norm values
   to obtain a global space-time residual norm for halting decisions is 
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
   certain time values at certain XBraid iterations.  Querying *astatus* for such information 
   is done through _braid_AccessStatusGet**(..)_ routines.
   \latexonly \\ \endlatexonly

   The frequency of the calls to *access* is controlled through 
   [braid_SetAccessLevel](@ref braid_SetAccessLevel).  For instance, if *access_level* is 
   set to 2, then *access* is called every XBraid iteration and on every XBraid level.  In 
   this case, querying *astatus* to determine the current XBraid level and iteration will 
   be useful. This scenario allows for even more detailed tracking of the simulation.
   \latexonly \\ \endlatexonly

   Eventually, this routine will allow for broader access to XBraid and computational steering.
   \latexonly \\ \endlatexonly
   
   See ``examples/ex-02`` and ``drivers/drive-04`` for more advanced uses of the
   *access* function.  In ``drive-04``, *access* is used to write solution
   vectors to a GLVIS visualization port, and ``ex-02`` uses *access* to write
   to .vtu files.

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
        
            sprintf(filename, "%s.%07d.%05d", "ex-01.out", index, myid);
            file = fopen(filename, "w");
            fprintf(file, "%.14e\n", (u->value));
            fflush(file);
            fclose(file);
        
            return 0;
         }

8. **BufSize**, **BufPack**, **BufUnpack**: These three routines tell XBraid how
   to communicate vectors between processors.  *BufPack* packs a vector into a
   ``void *`` buffer for MPI and then *BufUnPack* unpacks it from ``void *`` to
   vector.  Here doing that for a scalar is trivial.  *BufSize* computes the
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

9. **SCoarsen**, **SRestrict** (optional): These are advanced options that allow
   for coarsening in space while you coarsen in time.  This is useful for
   maintaining stable explicit schemes on coarse time scales and is not needed
   here.  See for instance ``drivers/drive-04`` and ``drivers/drive-02`` which
   use these routines.
   \latexonly \\ \endlatexonly

   These functions allow you vary the spatial mesh size on XBraid levels as depicted here
   where the spatial and temporal grid sizes are halved every level.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.25\textwidth]{../img/spatial_coarsening.pdf}}
       \label{img:heat_results}
   \end{figure}
   \endlatexonly

10. **Residual** (optional): A user-defined residual can be provided with the
    function [braid_SetResidual](@ref braid_SetResidual).  *Residual* defines
    the nonlinear equation to solve at each time step with the *Step* routine.
    Because XBraid assumes a one-step method, the equation to solve on each grid
    level has the form

    \f[ A(u_i, u_{i-1}) = f_i, \f]

    where \f$ A() \f$ is the *Residual* function and \f$ f_i = 0 \f$ on the
    finest grid level.  The nonzero right-hand-side on each coarse grid is
    needed to implement the FAS algorithm.  The *Step* function provides an
    approximation to the solution \f$ u_i \f$ of this equation.  That is,

    \f[ u_i \approx \Phi(u_{i-1}, f_i), \f]

    where \f$ \Phi() \f$ is the *Step* function.  When *Residual* is provided,
    *Step* does not need to produce an accurate time step, especially when all
    time points are stored (see [braid_SetStorage](@ref braid_SetStorage)).
    Users should write the *Residual* and *Step* routines such that the
    following holds:
    
    \f[ A( \Phi(u_{i-1}, f_i), u_{i-1} ) \approx f_i. \f]

    Note that when *Residual* is not provided, XBraid defines the residual in
    terms of the *Step* function as follows:

    \f[ A(u_i, u_{i-1}) = u_i - \Phi(u_{i-1}) = 0. \f]

    In this case, we have exact equality above (i.e., \f$ \approx \f$ becomes
    \f$ = \f$).  In addition, the nonzero right-hand-side needed on coarse grids
    to do FAS does not need to be provided to the user because the solution to
    the non-homogeneous equation is simply \f$ u_i = \Phi(u_{i-1}) + f_i \f$,
    and this can easily be computed internally in XBraid.  Also note that here
    *Step* must always be an accurate time step on the finest grid level.

11. Adaptive and variable time stepping are available by first calling the
    function [braid_SetRefine](@ref braid_SetRefine) in the main driver and then
    using [braid_StepStatusSetRFactor](@ref braid_StepStatusSetRFactor) in the
    *Step* routine to set a refinement factor for interval [*tstart*, *tstop*].
    In this example, no refinement is being done (factor = 1).  Currently, each
    refinement factor is constrained to be no larger than the coarsening factor.
    The final time grid is constructed adaptively in an FMG-like cycle by
    refining the initial grid according to the requested refinement factors.
    Refinement stops when the requested factors are all one or when various
    upper bounds are reached such as the max number of time points or max number
    of time grid refinement levels allowed.

## Running XBraid for this Example 

A typical flow of events in the *main* function is to first initialize the *app*
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
               my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
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

Finally, to run ex-01, type

    ex-01 -ml 5

This will run ``ex-01``. See ``examples/ex-0*`` for more extensive examples.


# Two-Dimensional Heat Equation {#exampletwo}

In this example, we assume familiarity with @ref exampleone and describe the
major ways in which this example differs.  This example is a full space-time
parallel example, as opposed to @ref exampleone, which implements only a scalar
ode for one degree-of-freedom in space.  We solve the heat equation in 2D,

\f[ \delta/\delta_t \; u(x,y,t) = \Delta\, u(x,y,t) + g(x,y,t). \f]

For spatial parallelism, we rely on the
[hypre](https://computation.llnl.gov/project/linear_solvers/software.php)
package where the SemiStruct interface is used to define our spatial
discretization stencil and form our time stepping scheme, the backward Euler
method.  The spatial discretization is just the standard 5-point finite
difference stencil (\f$[-1;-1,4,-1;-1]\f$), scaled by mesh widths, and the PFMG
solver is used for the solves required by backward Euler.  Please see the hypre
manual and examples for more information on the SemiStruct interface and PFMG.
Although, the hypre specific calls have mostly been abstracted away for this
example, and so it is not necessary to be familiar with the SemiStruct
interface for this example.

This example consists of three files and two executables.

   - examples/ex-02-serial.c:  This file compiles into its own executable
     ``ex-02-serial`` and represents a simple example user application.  This file
     supports only parallelism in space and represents a basic approach to
     doing efficient sequential time stepping with the backward Euler scheme.
     Note that the hypre solver used (PFMG) to carry out the time stepping is
     highly efficient.

   - examples/ex-02.c:  This file compiles into its own executable
     ``ex-02`` and represents a basic example of wrapping the user application 
     ``ex-02-serial``.  We will go over the wrappers below.

   - ex-02-lib.c:  This file contains shared functions used by the time-serial
     version and the time-parallel version.  This is where most of the hypre
     specific calls reside.  This file provides the basic functionality of this
     problem.  For instance, *take_step(u, tstart, tstop, ...)* carries out 
     a step, moving the vector *u* from time *tstart* to time *tstop* 
     and *setUpImplicitMatrix(...)* constructs the matrix to be inverted by PFMG
     for the backward Euler method.

## User Defined Structures and Wrappers

We now discuss in more detail the important data structures and wrapper
routines in ``examples/ex-02.c``.  The actual code for this example is quite simple
and it is recommended to read through it after this overview.

The two data structures are:
1. **App**: This holds a wide variety of information and is *global* in that it
   is passed to every user function.  This structure holds everything that the
   user will need to carry out a simulation.  One important structure contained in
   the *app* is the *simulation_manager*.  This is a structure native to the user
   code ``ex-02-lib.c``.  This structure conveniently holds the information needed by
   the user code to carry out a time step.  For instance, 

        app->man->A

   is the time stepping matrix,
      
        app->man->solver

   is the hypre PFMG solver object,
      
        app->man->dt

   is the current time step size.  The app is defined as 
      
  
      typedef struct _braid_App_struct {
         MPI_Comm             comm;             /* global communicator */
         MPI_Comm             comm_t;           /* communicator for parallelizing in time  */
         MPI_Comm             comm_x;           /* communicator for parallelizing in space  */
         int                  pt;               /* number of processors in time  */
         simulation_manager  *man;              /* user's simulation manager structure */
         HYPRE_SStructVector  e;                /* temporary vector used for error computations */
         int                  nA;               /* number of spatial matrices created */
         HYPRE_SStructMatrix *A;                /* array of spatial matrices, size nA, one per level*/
         double              *dt_A;             /* array of time step sizes, size nA, one per level*/
         HYPRE_StructSolver  *solver;           /* array of PFMG solvers, size nA, one per level*/
         int                  use_rand;         /* binary value, use random or zero initial guess */
         int                 *runtime_max_iter; /* runtime info for number of PFMG iterations*/
         int                 *max_iter_x;       /* maximum iteration limits for PFMG */
      } my_App;

   The app contains all the information needed to take a time step with the
   user code for an arbitrary time step size.  See the *Step* function below
   for more detail.

2. **Vector**: this defines a state vector at a certain time value.  
  Here, the vector is a structure containing a native hypre data-type, the
  *SStructVector*, which describes a vector over the spatial grid. Note that
  *my_Vector* is used to define *braid_Vector*.

       typedef struct _braid_Vector_struct {
          HYPRE_SStructVector   x;
       } my_Vector;


The user must also define a few wrapper routines.  Note, that the ``app`` structure
is the first argument to every function.

1. **Step**: This function tells XBraid how to take a time step, and is the core
   user routine.  This function advances the vector *u* from time *tstart* to time
   *tstop*.  A few important things to note are as follows.

  + The time values are given to the user through the *status* structure
    and associated *Get* routines.  
  
  + The basic strategy is to see if a matrix and solver already exist for 
  this \f$dt\f$ value.  If not, generate a new matrix and solver and store them in
  the *app* structure.  If they do already exist, then re-use the data.  
  
  + To carry out a step, the user routines from ``ex-02-lib.c`` rely on a few
  crucial data members *man->dt*, *man->A* and *man-solver*.  We overwrite
  these members with the correct information for the time step size in question.
  Then, we pass *man* and *u* to the user function *take_step(...)* which
  evolves *u*.

  + The forcing term \f$ g_i \f$ is wrapped into the *take_step(...)*
    function.  Thus, \f$\Phi(u_i) \rightarrow u_{i+1} \f$.  

        int my_Step(braid_App        app,
                    braid_Vector     u,
                    braid_StepStatus status)
        {
           double tstart;             /* current time */
           double tstop;              /* evolve u to this time*/
           int i, A_idx;
           int iters_taken = -1;
           
           /* Grab status of current time step */
           braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
     
           /* Check matrix lookup table to see if this matrix already exists*/
           A_idx = -1.0;
           for( i = 0; i < app->nA; i++ ){
              if( fabs( app->dt_A[i] - (tstop-tstart) )/(tstop-tstart) < 1e-10) { 
                 A_idx = i;
                 break;
              }
           }
     
           /* We need to "trick" the user's manager with the new dt */
           app->man->dt = tstop - tstart;
     
           /* Set up a new matrix and solver and store in app */
           if( A_idx == -1.0 ){
              A_idx = i;
              app->nA++;
              app->dt_A[A_idx] = tstop-tstart;
     
              setUpImplicitMatrix( app->man );
              app->A[A_idx] = app->man->A;
              
              setUpStructSolver( app->man, u->x, u->x );
              app->solver[A_idx] = app->man->solver;
           } 
     
           /* Time integration to next time point: Solve the system Ax = b.
            * First, "trick" the user's manager with the right matrix and solver */ 
           app->man->A = app->A[A_idx];
           app->man->solver = app->solver[A_idx];
           ...
           /* Take step */
           take_step(app->man, u->x, tstart, tstop);
           ...
           return 0;
        }

2. There are other functions, **Init**, **Clone**, **Free**, **Sum**,
**SpatialNorm**, **Access**, **BufSize**, **BufPack** and **BufUnpack**, which
also must be written.  These functions are all simple for this example, as for
the case of @ref exampleone. All we do here is standard operations on a spatial
vector such as initialize, clone, take an inner-product, pack, etc... We refer
the reader to ``ex-02.c``.

## Running XBraid for this Example 

To initialize and run XBraid, the procedure is similar to  @ref exampleone.
Only here, we have to both initialize the user code and XBraid.  The code that
is specific to the user's application comes directly from the existing serial
simulation code.  If you compare ``ex-02-serial.c`` and ``ex-02.c``, you will
see that most of the code setting up the user's data structures and defining
the wrapper functions are simply lifted from the serial simulation.

Taking excerpts from the function *main()* in ex-02.c,
we first initialize the user's simulation manager with code like
      
      ...
      app->man->px    = 1;   /* my processor number in the x-direction */
      app->man->py    = 1;   /* my processor number in the y-direction */
                             /* px*py=num procs in space */
      app->man->nx    = 17;  /* number of points in the x-dim */
      app->man->ny    = 17;  /* number of points in the y-dim */
      app->man->nt    = 32;  /* number of time steps */
      ...

We also define default XBraid parameters with code like
      
      ...
      max_levels      = 15; /* Max levels for XBraid solver */
      min_coarse      = 3;  /* Minimum possible coarse grid size */
      nrelax          = 1;  /* Number of CF relaxation sweeps on all levels */
      ...

The XBraid app must also be initialized with code like 
    
    ...
    app->comm   = comm;
    app->tstart = tstart;
    app->tstop  = tstop;
    app->ntime  = ntime;

Then, the data structure definitions and wrapper routines are passed to XBraid.

    braid_Core  core;
    braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
            my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
            my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

Then, XBraid options are set with calls like

    ...
    braid_SetPrintLevel( core, 1);
    braid_SetMaxLevels(core, max_levels);
    braid_SetNRelax(core, -1, nrelax);
    ...

Then, the simulation is run.

    braid_Drive(core);

Then, we clean up.

    braid_Destroy(core);

Finally, to run ex-02, type

    ex-02 -help

As a simple example, try the following.

    mpirun -np 8 ex-02 -pgrid 2 2 2 -nt 256

## Scaling Study with this Example {#twodheat_scaling}

Here, we carry out a simple strong scaling study for this example.  The "time
stepping" data set represents sequential time stepping and was generated using
``examples/ex-02-serial``.  The time-parallel data set was generated using
``examples/ex-02``.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.45\textwidth]{../img/heat_results.pdf}}
       \label{img:heat_results}
   \end{figure}
   \endlatexonly
The problem setup is as follows.
- Backwards Euler is used as the time stepper.  This is the only time stepper
supported by ``ex-02``.
- We used a Linux cluster with 4 cores per node, a Sandybridge Intel chipset, 
  and a fast Infiniband interconnect.  
- The space-time problem size was \f$ 129^2 \times 16,192 \f$ over the unit cube 
  \f$ [0,1] \times [0,1] \times [0,1] \f$ .  
- The coarsening factor was \f$ m = 16 \f$ on the finest level and \f$ m=2 \f$ 
  on coarser levels.  
- Since 16 processors optimized the serial time stepping approach, 16 processors
  in space are also used for the XBraid experiments.  So for instance 512 processrs in 
  the plot corresponds to 16 processors in space and 32 processors in time, 
  \f$ 16*32 = 512 \f$.  Thus, each processor owns a space-time hypercube of
  \f$ (129^2 / 16) \times (16,192 / 32) \f$.  See @ref decomposition for a depiction
  of how XBraid breaks the problem up.
- Various relaxation and V and F cycling strategies are experimented with.  
  + *V-cycle, FCF* denotes V-cycles and FCF-relaxation on each level.
  + *V-cycle, F-FCF* denotes V-cycles and F-relaxation on the finest level and
     FCF-relaxation on all coarser levels.
  + *F-cycle, F* denotes F-cycles and F-relaxation on each level.
- The initial guess at time values for \f$ t > 0 \f$ is zero, which is typical.
- The halting tolerance corresponds to a discrete L2-norm and was 
  \f[\mbox{tol} = \frac{10^{-8}}{\sqrt{ (h_x)^2 h_t} }, \f]
  where \f$ h_x \f$ and \f$ h_t \f$ are the spatial and temporal grid spacings, 
  respectively.

  This corresponds to passing *tol* to [braid_SetAbsTol](@ref braid_SetAbsTol),
  passing *2* to [braid_SetTemporalNorm](@ref braid_SetTemporalNorm) and
  defining [braid_PtFcnSpatialNorm](@ref braid_PtFcnSpatialNorm) to be the
  standard Euclidean 2-norm.  All together, this appropriately scales the
  space-time residual in way that is relative to the number of space-time grid
  points (i.e., it approximates the L2-norm).


To re-run this scaling study, a sample run string for ex-02 is

      mpirun -np 64 ex-02 -pgrid 4 4 4 -nx 129 129 -nt 16129 -cf0 16 -cf 2 -nu 1 -use_rand 0

To re-run the baseline sequential time stepper, ex-02-serial, try

      mpirun -np 64 ex-02-serial -pgrid 8 8 -nx 129 129 -nt 16129

For explanations of the command line parameters, type
     
      ex-02-serial -help
      ex-02 -help


Regarding the performance, we can say
- The best speedup is 10x and this would grow if more processors were available.
- Although not shown, the iteration counts here are about 10-15 XBraid iterations.
  See [Parallel Time Integration with Multigrid](https://computation.llnl.gov/project/linear_solvers/pubs/mgritPaper-2014.pdf)
  for the exact iteration counts.
- At smaller core counts, serial time stepping is faster.  But at about 256 processors,
  there is a crossover and XBraid is faster.
- You can see the impact of the cycling and relaxation strategies discussed in
  @ref cyclingrelaxation.  For instance, even though *V-cycle, F-FCF* is a 
  weaker relaxation strategy than *V-cycle, FCF* (i.e., the XBraid convergence 
  is slower), *V-cycle, F-FCF* has a faster time to solution than *V-cycle, FCF* 
  because each cycle is cheaper.
- In general, one level of aggressive coarsening (here by a factor 16) followed by
  slower coarsening was found to be best on this machine.

Achieving the best speedup can require some tuning, and it is recommended to read
[Parallel Time Integration with Multigrid](https://computation.llnl.gov/project/linear_solvers/pubs/mgritPaper-2014.pdf)
where this 2D heat equation example is explored in much more detail.


# Running and Testing XBraid

The best overall test for XBraid, is to set the maximum number of levels to 1 
(see [braid_SetMaxLevels](@ref braid_SetMaxLevels) ) which will carry out a
sequential time stepping test.  Take the output given to you by your *Access*
function and compare it to output from a non-XBraid run.  Is everything OK?
Once this is complete, repeat for multilevel XBraid, and check that the solution
is correct (that is, it matches a serial run to within tolerance).

At a lower level, to do sanity checks of your data structures and wrapper routines, there are
also XBraid test functions, which can be easily run.  The test routines also 
take as arguments the *app* structure, spatial communicator *comm_x*, a stream
like *stdout* for test output and a time step size *dt* to test.  After these arguments,
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

