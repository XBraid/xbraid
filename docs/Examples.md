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

This section is the chief *tutorial* of XBraid, illustrating how to use it through a sequence
of progressively more sophisticated examples.

# The Simplest Example {#exampleone}

## User Defined Structures and Wrappers

The user must wrap their existing time stepping routine per the XBraid interface. 
To do this, the user must define two data structures and some wrapper routines.  To make the 
idea more concrete, we now give these function definitions from ``examples/ex-01``, which 
implements a scalar ODE, \f[ u_t = \lambda u. \f]

The two data structures are:

1. **App**: This holds a wide variety of information and is *global* in that it
   is passed to every function.  This structure holds everything that the user
   will need to carry out a simulation.  Here for illustration, this is just an 
   integer storing a processor's rank. 
   
        typedef struct _braid_App_struct
        {
           int       rank;
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
   *status* structure and associated *Get* routine.  **Important note:** 
   the \f$ g_i \f$ function from @ref braidoverview must be incorporated into 
   *Step*, so that the following equation is solved by default. 
   \f[ \Phi(u_i) = 0. \f]
   The *ustop* parameter
   serves as an approximation to the solution at time *tstop* and is not needed
   here.  It can be useful for implicit schemes that require an initial guess
   for a linear or nonlinear solver.  The use of *fstop* is an advanced parameter 
   (not required) and forms the the right-hand
   side of the nonlinear problem on the given time grid.  This value is only nonzero when
   providing a residual with [braid_SetResidual](@ref braid_SetResidual).  More
   information on how to use this optional feature is given below.  
   \latexonly \\ \endlatexonly
   
   Here advancing the solution just involves the scalar \f$ \lambda \f$.
         
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
      
            /* Use backward Euler to propagate solution */
            (u->value) = 1./(1. + tstop-tstart)*(u->value);
            
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
   The default *access_level* is 1 and gives the user access only after the simulation ends
   and only on the finest time-grid.
   \latexonly \\ \endlatexonly

   Eventually, this routine will allow for broader access to XBraid and computational steering.
   \latexonly \\ \endlatexonly
   
   See ``examples/ex-03`` and ``drivers/drive-diffusion`` for more advanced
   uses of the *access* function.  In ``drive-diffusion``, *access* is used to write
   solution vectors to a GLVIS visualization port, and ``ex-03`` uses *access*
   to write to .vtu files.

         int
         my_Access(braid_App          app,
                   braid_Vector       u,
                   braid_AccessStatus astatus)
         {
            int        index;
            char       filename[255];
            FILE      *file;
            
            braid_AccessStatusGetTIndex(astatus, &index);
            sprintf(filename, "%s.%04d.%03d", "ex-01.out", index, app->rank);
            file = fopen(filename, "w");
            fprintf(file, "%.14e\n", (u->value));
            fflush(file);
            fclose(file);
      
            return 0;
         }

8. **BufSize**, **BufPack**, **BufUnpack**: These three routines tell XBraid how to 
   communicate vectors between processors.  *BufPack* packs a vector 
   into a ``void *`` buffer for MPI and then *BufUnPack* unpacks the ``void *`` buffer
   into a vector.  Here doing that for a scalar is trivial.  *BufSize* computes the 
   upper bound for the size of an arbitrary vector.
   \latexonly \\ \endlatexonly

   Note how *BufPack* also sets the size in *bstatus*.  This value is optional,
   but if set it should be the exact number of bytes packed, while *BufSize*
   should provide only an upper-bound on a possible buffer size.  This
   flexibility allows for the buffer to be allocated the fewest possible times,
   but smaller messages to be sent when needed.  For instance, this occurs when
   using variable spatial grid sizes.  **To avoid MPI issues, it is very
   important that BufSize be pessimistic, provide an upper bound, and return
   the same value across processors.**
   \latexonly \\ \endlatexonly

   In general, the buffer should be self-contained.  The receiving processor
   should be able to pull all necessary information from the buffer in order to
   properly interpret and unpack the buffer.

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

## Running XBraid for this Example 

A typical flow of events in the *main* function is to first initialize the *app*
structure.

    /* set up app structure */
    app = (my_App *) malloc(sizeof(my_App));
    (app->rank)   = rank;

Then, the data structure definitions and wrapper routines are passed to XBraid.
The core structure is used by XBraid for internal data structures. 

    braid_Core  core;
    braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
               my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
               my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

Then, XBraid options are set.

    braid_SetPrintLevel( core, 1);
    braid_SetMaxLevels(core, max_levels);
    braid_SetAbsTol(core, tol);
    braid_SetCFactor(core, -1, cfactor);
   
Then, the simulation is run.

    braid_Drive(core);

Then, we clean up.

    braid_Destroy(core);

Finally, to run ex-01, type

    ex-01



# Some Advanced Features {#advancedfeatures}

We now give an overview of some *optional* advanced features, together in one place, that will 
be implemented in some of the following examples. 

9. **SCoarsen**, **SRestrict**: These are advanced options that allow
   for coarsening in space while you coarsen in time.  This is useful for
   maintaining stable explicit schemes on coarse time scales and is not needed
   here.  See for instance ``drivers/drive-diffusion`` and ``drivers/drive-diffusion-2D`` 
   which use these routines.
   \latexonly \\ \endlatexonly

   These functions allow you to vary the spatial mesh size on XBraid levels as depicted here
   where the spatial and temporal grid sizes are halved every level.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.25\textwidth]{../img/spatial_coarsening.pdf}}
       \label{img:heat_results}
   \end{figure}
   \endlatexonly

10. **Residual**: A user-defined residual can be provided with the
    function [braid_SetResidual](@ref braid_SetResidual) and can result
    in substantial computational savings, as explained below.  
    
    However to use this advanced feature, one must first understand how XBraid 
    measures the residual.  XBraid computes residuals of this equation,
    
    \f[ A_i(u_i, u_{i-1}) = f_i, \f]

    where \f$ A_i(,) \f$ evaluates one block-row of the the global space-time
    operator \f$A\f$.  The forcing \f$f_i\f$ is the XBraid forcing, the FAS
    right-hand-side term on coarse-grids and 0 on the finest grid.  The PDE forcing
    goes inside of \f$A_i\f$. 
    
    Since XBraid assumes one-step methods, define \f$A_i()\f$ to be

    \f[ A_i(u_i, u_{i-1}) = - \Phi(u_{i-1} + \Psi(u_i), \f]
    
    i.e., the subdiagonal and diagonal blocks of \f$A\f$.
    
    **Default setting** In the default XBraid setting (no residual option
    used), the user only defines *Step()* and here that means defining only the
    \f$\Phi()\f$ function because \f$\Psi()\f$ is assumed to be the identity.
    Thus, XBraid can compute the residual using only a user-defined
    \f$\Phi()\f$ in *Step()* by combining *Step()* with the *Sum()* function, i.e.

    \f[ r_i = f_i + \Phi(u_{i-1}) -u_i. \f]
    
    The *fstop* parameter in *Step()* corresponds to \f$f_i\f$ and need not
    be operated on my user in order to compute the residual \f$r_i\f$.  Thus
    *fstop* is NULL in the default setting.  

    An implication of this is the evaluation of \f$\Phi()\f$ on the finest grid
    must be very accurate, or the residual will not be accurate.  This leads to a
    nonintrusive, but expensive algorithm.  The accuracy of \f$\Phi()\f$ can be relaxed
    on coarser grids.  
    
    **Residual setting** 
    The alternative to the above default non-intrusive strategy is to have the user
    define 

    \f[ A_i(u_i, u_{i-1}) = - \Phi(u_{i-1} + \Psi(u_i), \f]

    directly, which is the case when the *Residual* function is defined with
    [braid_PtFcnResidual](@ref braid_PtFcnResidual).  In other words, the user
    now defines each block-row of the space-time operator, rather than only
    defining \f$\Phi()\f$.  The user function computes \f$A_i(u_i, u_{i-1})\f$
    and XBraid then subtracts this from \f$f_i\f$ to compute \f$r_i\f$.

    However, more care must now be taken when defining the *Step()* function.
    In particular, the *fstop* value (\f$f_i\f$) value must be taken into account.
    Overall, the definition of *Step()* now changes.  It no longer defines  \f$\Phi()\f$,
    it instead defines a (possibly inexact) solve of the equation defined by 
    
    \f[ A_i(u_i, u_{i-1}) = f_i. \f]

    Thus, *Step()* must be compatible with *Residual()*. Essentially, *Step()*
    must now compute 
    
    \f[ u_i = \Psi^{-1}(f_i - \Phi(u_{i-1})). \f]
    
    It is clear that the *fstop* (\f$f_i\f$) value must be given to the
    *Step()* function so that this equation must be solved.  Essentially, one
    can think of *Residual()* as defining the equation, and *Step()* defining a
    preconditioner for that row of the equation, an inexact solve for
    \f$u_i\f$. 
   
    As an example, let \f$\Psi = (I + \Delta t L)\f$, where \f$L\f$ is a
    Laplacian and \f$\Phi = I\f$.  The application of the residual function
    will only be a sparse matrix-vector application, as opposed to the default
    case where \f$\Phi = (I + \Delta t L)^{-1}\f$ and  \f$\Psi = I\f$.  This
    results in considerable computational savings.  The application of *Step()*
    then involves an inexact inversion of  \f$\Psi\f$, e.g., with just 1
    spatial multigrid V-cycle, which again results in substantial computation
    savings over a full inversion.

    Another way to think about the compatibility between \f$\Psi\f$ and \f$\Phi\f$
    is that 
    
    \f[ f_i - A_i( u_i, u_{i-1} ) = 0 \f]

    must hold exactly if \f$u_i\f$ is an exact propagation of \f$u_{i-1}\f$,
    or 
    \f[ f_i - A_i( Step(u_{i-1},f_i), u_{i-1} ) = 0 \f]
    must hold.  When the accuracy of the *Step()* function is reduced (as mentioned
    above), this exact equality with 0 is lost, but this should evaluate to something
    ``small``.  There is an XBraid test function [braid_TestResidual](@ref braid_TestResidual)
    that tests for this compatibility.

    The residual feature is implemented in the examples
    ``examples/ex-01-expanded.c``, ``examples/ex-02.c``, and
    ``examples/ex-03.c``.


11. **Adaptive and variable time stepping**:  This feature is available by first calling the
   function [braid_SetRefine](@ref braid_SetRefine) in the main driver and then
   using [braid_StepStatusSetRFactor](@ref braid_StepStatusSetRFactor) in the
   *Step* routine to set a refinement factor for interval [*tstart*, *tstop*].
   In this way, user-defined criteria can subdivide intervals on the fly and adaptively
   refine in time.  For instance, returning a refinement factor of 4 in *Step* 
   will tell XBraid to subdivide that interval into 4 evenly spaced smaller intervals
   for the next iteration.  Refinement can only be done on the finest XBraid level.
   \latexonly \\ \endlatexonly

   In this example, no refinement is being done (factor = 1).  Currently, each
   refinement factor is constrained to be no larger than the coarsening factor.
   The final time grid is constructed adaptively in an FMG-like cycle by
   refining the initial grid according to the requested refinement factors.
   Refinement stops when the requested factors are all one or when various
   upper bounds are reached such as the max number of time points or max number
   of time grid refinement levels allowed.

12. **Shell-vector**: This feature supports the use of multi-step methods.
   It is used in example ``examples/ex-01-expanded-bdf2.c``.

   The strategy for BDF-K methods is to allow for the lumping of ``k`` time
   points into a single XBraid vector.  So, if the problem had 100 time points
   and the time-stepper was BDF-2, then XBraid would only ``see`` 50 time
   points but each XBraid vector would contain two separate time points.  By
   lumping 2 time points into one vector, the BDF-2 scheme remains one-step and
   compatible with XBraid.  

   However, the time-point spacing between the two points internal to the
   vector stays the same on all time grids, while the spacing between vectors
   grows on coarse time grids.  This creates an irregular spacing which is
   problematic for BDF-k methods and lends itself to the shell-vector strategy
   that lets meta-data be stored at all time points, even for F-points which are
   usually not stored.

   There are many strategies for handling the coarse time-grids (dropping the BDF
   order, adjusting time-point spacings inside the vector, etc...).  Prospective
   users are encouraged to contact xbraid-support@llnl.gov because this is active
   research.

13. **Storage**:  

   This option (see [braid_SetStorage](@ref braid_SetStorage)) allows the user
   to specify storage at all time points (C and F) or only at C-points.  This
   extra storage is useful for implicit methods, where the solution value from
   the *previous XBraid iteration* for a time step can be used as the initial
   guess to the implicit solver.  This is often a better initial guess than
   using the solution value from the previous time step.  The default is to
   store only C-point values, thus the better initial guess is only available
   there.
   
   In general, the user should always use the *ustop* parameter in
   *Step()* as the initial guess for an implicit solve.  If storage
   is turned on (set to 0), then this value will always be the improved
   initial guess for C- and F-points.  If storage is not turned on, then
   this will be the improved guess only for C-points.  For F-points,
   it will equal the solution from the previous time step.

   See ``examples/ex-03`` for an example which uses this feature.

# Simplest example expanded {#exampleoneexapanded}

These examples build on @ref exampleone, but still solve the 
scalar ODE, 

\f[ u_t = \lambda u. \f]

The goal here is to show more advanced features of XBraid.

   - ``examples/ex-01-expanded.c``:  same as ``ex-01.c`` but adds more XBraid features
     such as the residual feature, the user defined initial time-grid and full multigrid
     cycling.

   - ``examples/ex-01-expanded-bdf2.c``:  same as ex-01-expanded.c, but uses
     BDF2 instead of backward Euler.  This example makes use of the advanced 
     shell-vector feature in order to implement BDF2.
    
   - ``examples/ex-01-expanded-f.f90``:  same as ex-01-expanded.c, but
     implemented in f90


# One-Dimensional Heat Equation {#exampletwo}

In this example, we assume familiarity with @ref exampleone.  This example is a
time-only parallel example that implements the 1D heat equation, 

\f[ \delta/\delta_t \; u(x,t) = \Delta\, u(x,t) + g(x,t), \f]

as opposed to @ref exampleone, which implements only a scalar ODE for one
degree-of-freedom in space.  There is not spatial parallelism, as a serial
cyclic reduction algorithm is used to invert the tri-diagonal spatial operators.
The space-time discretization is the standard 3-point finite
difference stencil (\f$[-1,2,-1]\f$), scaled by mesh widths coupled with 
backward Euler.

This example consists of three files and two executables.

   - examples/ex-02-serial.c:  This file compiles into its own executable
     ``ex-02-serial`` and represents a simple example user application that
     does sequential time-stepping.  This file represents where a new 
     XBraid user would start, in terms of converting a sequential time-stepping
     code to XBraid. 

   - examples/ex-02.c:  This file compiles into its own executable
     ``ex-02`` and represents an XBraid wrapping of the user application 
     ``ex-02-serial``.  

   - ex-02-lib.c:  This file contains shared functions used by the time-serial
     version and the time-parallel version.   This file provides the basic functionality of this
     problem.  For instance, *take_step(u, tstart, tstop, ...)* carries out 
     a step, moving the vector *u* from time *tstart* to time *tstop*. 

# Two-Dimensional Heat Equation {#examplethree}

In this example, we assume familiarity with @ref exampleone and describe the
major ways in which this example differs.  This example is a full space-time
parallel example, as opposed to @ref exampleone, which implements only a scalar
ODE for one degree-of-freedom in space.  We solve the heat equation in 2D,

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

   - examples/ex-03-serial.c:  This file compiles into its own executable
     ``ex-03-serial`` and represents a simple example user application.  This file
     supports only parallelism in space and represents a basic approach to
     doing efficient sequential time stepping with the backward Euler scheme.
     Note that the hypre solver used (PFMG) to carry out the time stepping is
     highly efficient.

   - examples/ex-03.c:  This file compiles into its own executable
     ``ex-03`` and represents a basic example of wrapping the user application 
     ``ex-03-serial``.  We will go over the wrappers below.

   - ex-03-lib.c:  This file contains shared functions used by the time-serial
     version and the time-parallel version.  This is where most of the hypre
     specific calls reside.  This file provides the basic functionality of this
     problem.  For instance, *take_step(u, tstart, tstop, ...)* carries out 
     a step, moving the vector *u* from time *tstart* to time *tstop* 
     and *setUpImplicitMatrix(...)* constructs the matrix to be inverted by PFMG
     for the backward Euler method.

## User Defined Structures and Wrappers

We now discuss in more detail the important data structures and wrapper
routines in ``examples/ex-03.c``.  The actual code for this example is quite simple
and it is recommended to read through it after this overview.

The two data structures are:
1. **App**: This holds a wide variety of information and is *global* in that it
   is passed to every user function.  This structure holds everything that the
   user will need to carry out a simulation.  One important structure contained in
   the *app* is the *simulation_manager*.  This is a structure native to the user
   code ``ex-03-lib.c``.  This structure conveniently holds the information needed by
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
  
  + To carry out a step, the user routines from ``ex-03-lib.c`` rely on a few
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
the reader to ``ex-03.c``.

## Running XBraid for this Example 

To initialize and run XBraid, the procedure is similar to  @ref exampleone.
Only here, we have to both initialize the user code and XBraid.  The code that
is specific to the user's application comes directly from the existing serial
simulation code.  If you compare ``ex-03-serial.c`` and ``ex-03.c``, you will
see that most of the code setting up the user's data structures and defining
the wrapper functions are simply lifted from the serial simulation.

Taking excerpts from the function *main()* in ex-03.c,
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

Finally, to run ex-03, type

    ex-03 -help

As a simple example, try the following.

    mpirun -np 8 ex-03 -pgrid 2 2 2 -nt 256

## Scaling Study with this Example {#twodheat_scaling}

Here, we carry out a simple strong scaling study for this example.  The "time
stepping" data set represents sequential time stepping and was generated using
``examples/ex-03-serial``.  The time-parallel data set was generated using
``examples/ex-03``.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.45\textwidth]{../img/heat_results.pdf}}
       \label{img:heat_results}
   \end{figure}
   \endlatexonly
The problem setup is as follows.
- Backwards Euler is used as the time stepper.  This is the only time stepper
supported by ``ex-03``.
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


To re-run this scaling study, a sample run string for ex-03 is

      mpirun -np 64 ex-03 -pgrid 4 4 4 -nx 129 129 -nt 16129 -cf0 16 -cf 2 -nu 1 -use_rand 0

To re-run the baseline sequential time stepper, ex-03-serial, try

      mpirun -np 64 ex-03-serial -pgrid 8 8 -nx 129 129 -nt 16129

For explanations of the command line parameters, type
     
      ex-03-serial -help
      ex-03 -help


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

# More Complicated Examples

We have Fortran90 and C++ interfaces.  See ``examples/ex-01f.f90``, ``braid.hpp``
and the various C++ examples in ``drivers/drive-**.cpp``.  For discussion of
more complex problems please see our project
[publications website](http://computation.llnl.gov/projects/parallel-time-integration-multigrid/publications)
for our recent publications concerning some of these varied applications. 


