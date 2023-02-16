<!--
  - Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
  - Produced at the Lawrence Livermore National Laboratory. Written by 
  - Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
  - Dobrev, et al. LLNL-CODE-660355. All rights reserved.
  - 
  - This file is part of XBraid. For support, post issues to the XBraid Github page.
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

## Running XBraid for the Simplest Example {#running_simplestexample}

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

We now give an overview of some *optional* advanced features that will 
be implemented in some of the following examples. 

9. **SCoarsen**, **SRestrict**: These are advanced options that allow
   for coarsening in space while you coarsen in time.  This is useful for
   maintaining stable explicit schemes on coarse time scales and is not needed
   here.  See ``examples/ex-02`` for a simple example of this feature, and then 
   ``drivers/drive-diffusion`` and ``drivers/drive-diffusion-2D`` 
   for more advanced examples of this feature. 
   \latexonly \\ \endlatexonly

   These functions allow you to vary the spatial mesh size on XBraid levels as depicted here
   where the spatial and temporal grid sizes are halved every level.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.25\textwidth]{../img/spatial_coarsening.pdf}}
   \end{figure}
   \endlatexonly

10. **Residual**: A user-defined residual can be provided with the
   function [braid_SetResidual](@ref braid_SetResidual) and can result
   in substantial computational savings, as explained below.  
   However to use this advanced feature, one must first understand how XBraid 
   measures the residual.  XBraid computes residuals of this equation,
   
   \f[ A_i(u_i, u_{i-1}) = f_i, \f]

   where \f$ A_i(,) \f$ evaluates one block-row of the the global space-time
   operator \f$A\f$.  The forcing \f$f_i\f$ is the XBraid forcing, which is the FAS
   right-hand-side term on coarse grids and 0 on the finest grid.  The PDE forcing
   goes inside of \f$A_i\f$. 
   
   Since XBraid assumes one-step methods, \f$A_i()\f$ is defined to be

   \f[ A_i(u_i, u_{i-1}) = - \Phi(u_{i-1}) + \Psi(u_i), \f]
   
   i.e., the subdiagonal and diagonal blocks of \f$A\f$.
   \latexonly \\ \endlatexonly
   
   **Default setting**: In the default XBraid setting (no residual option
   used), the user only implements *Step()* and *Step()* will simply apply 
   \f$\Phi()\f$, because \f$\Psi()\f$ is assumed to be the identity.
   Thus, XBraid can compute the residual using only the user-defined
   *Step()* function by combining *Step()* with the *Sum()* function, i.e.

   \f[ r_i = f_i + \Phi(u_{i-1}) -u_i. \f]
   
   The *fstop* parameter in *Step()* corresponds to \f$f_i\f$, but is always
   passed in as NULL to the user in this setting and should be ignored.  This
   is because XBraid can compute the contribution of \f$f_i\f$ to the residual 
   on its own using the *Sum()* function.

   An implication of this is that the evaluation of \f$\Phi()\f$ on the finest grid
   must be very accurate, or the residual will not be accurate.  This leads to a
   nonintrusive, but expensive algorithm.  The accuracy of \f$\Phi()\f$ can be relaxed
   on coarser grids to save computations.
   \latexonly \\ \endlatexonly
    
   **Residual setting**: 
   The alternative to the above default least-intrusive strategy is to have the user
   define 

   \f[ A_i(u_i, u_{i-1}) = - \Phi(u_{i-1}) + \Psi(u_i), \f]

   directly, which is what the *Residual* function implements (set with 
   [braid_PtFcnResidual](@ref braid_PtFcnResidual)).  In other words, the user
   now defines each block-row of the space-time operator, rather than only
   defining \f$\Phi()\f$.  The user *Residual()* function computes \f$A_i(u_i, u_{i-1})\f$
   and XBraid then subtracts this from \f$f_i\f$ to compute \f$r_i\f$.

   However, more care must now be taken when defining the *Step()* function.
   In particular, the *fstop* value (i.e., the \f$f_i\f$ value) must be taken into account.
   Essentially, the definition of *Step()* changes so that it no longer defines \f$\Phi()\f$,
   but instead defines a (possibly inexact) solve of the equation defined by 
   
   \f[ A_i(u_i, u_{i-1}) = f_i. \f]

   Thus, *Step()* must be compatible with *Residual()*.  Expanding the previous
   equation, we say that *Step()* must now compute 
   
   \f[ u_i = \Psi^{-1}(f_i + \Phi(u_{i-1})). \f]
   
   It is clear that the *fstop* value (i.e., the \f$f_i\f$ value) must now be
   given to the *Step()* function so that this equation can be solved by the
   user.  In other words, *fstop* is now no longer NULL.  
   
   Essentially, one can
   think of *Residual()* as defining the equation, and *Step()* defining a
   preconditioner for that row of the equation, or an inexact solve for
   \f$u_i\f$. 
   
   As an example, let \f$\Psi = (I + \Delta t L)\f$, where \f$L\f$ is a
   Laplacian and \f$\Phi = I\f$.  The application of the residual function will
   only be a sparse matrix-vector multiply, as opposed to the default case
   where an inversion is required for \f$\Phi = (I + \Delta t L)^{-1}\f$ and
   \f$\Psi = I\f$.  This results in considerable computational savings.
   Moreover, the application of *Step()* now involves an inexact inversion of
   \f$\Psi\f$, e.g., by using just one spatial multigrid V-cycle. This again results
   in substantial computation savings when compared with the naive approach of
   a full matrix inversion.  \latexonly \\ \endlatexonly

   Another way to think about the compatibility between \f$\Psi\f$ and \f$\Phi\f$
   is that 
   
   \f[ f_i - A_i( u_i, u_{i-1} ) = 0 \f]

   must hold exactly if \f$u_i\f$ is an exact propagation of \f$u_{i-1}\f$,
   that is,
   \f[ f_i - A_i( Step(u_{i-1},f_i), u_{i-1} ) = 0 \f]
   must hold.  When the accuracy of the *Step()* function is reduced (as mentioned
   above), this exact equality with 0 is lost, but this should evaluate to something
   ``small``.  There is an XBraid test function [braid_TestResidual](@ref braid_TestResidual)
   that tests for this compatibility.
   \latexonly \\ \endlatexonly

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

   The final time grid is constructed adaptively in an
   FMG-like cycle by refining the initial grid according to the requested
   refinement factors.  Refinement stops when the requested factors are all
   one or when various upper bounds are reached such as the max number of time
   points or max number of time grid refinement levels allowed. No restriction on
   the refinement factors is applied within XBraid, so the user may want to apply his
   own upper bound on the refinement factors to avoid over-refinement.
   See ``examples/ex-01-refinement.c`` and ``examples/ex-03.c`` for an
   implementation of this.

12. **Richardson-based Error Estimation and Extrapolation**: This feature
    allows the user to access built-in Richardson-based error estimates and
    accuracy improving extrapolation.  The error estimates and/or extrapolation
    can be turned on by using 
    [braid_SetRichardsonEstimation](@ref braid_SetRichardsonEstimation) .
    Moreover, this feature can be used in conjunction with the above discussed function,
    [braid_StepStatusSetRFactor](@ref braid_StepStatusSetRFactor), to achieve easy-to-use
    adaptive refinement in time.

    Essentially, Richardson extrapolation (RE) is used to improve the accuracy
    of the solution at the C-points on the finest level.  When the built-in
    error estimate option is turned on, RE is used to estimate the local
    truncation error at each point. These estimates can be accessed through
    StepStatus and AccessStatus functions. 
   
    The Richardson-based error estimates and extrapolation are only available
    after the first Braid iteration, in that the coarse level solution must be
    available to compute the error estimate and/or extrapolation.  Thus, after
    an adaptive refinement (and new hierarchy is constructed), another
    iteration is again required for the error estimates to be available.  If
    the error estimate isn't available, Braid returns a value of -1.  See this
    example for more details 

               examples/ex-06.c
   
13. **Shell-vector**: This feature supports the use of multi-step methods.
   The strategy for BDF-K methods is to allow for the lumping of ``k`` time
   points into a single XBraid vector.  So, if the problem had 100 time points
   and the time-stepper was BDF-2, then XBraid would only ``see`` 50 time
   points but each XBraid vector would contain two separate time points.  By
   lumping 2 time points into one vector, the BDF-2 scheme remains one-step and
   compatible with XBraid.  

   However, the time-point spacing between the two points internal to the
   vector stays the same on all time grids, while the spacing between vectors
   grows on coarse time grids.  This creates an irregular spacing which is
   problematic for BDF-k methods.  Thus the shell-vector strategy
   lets meta-data be stored at all time points, even for F-points which are
   usually not stored, so that the irregular spacings can be tracked and 
   accounted for with the BDF method.  (Note, there are other possible uses
   for shell-vectors.)

   There are many strategies for handling the coarse time-grids with BDF
   methods (dropping the BDF order, adjusting time-point spacings inside the
   lumped vectors, etc...).  Prospective users are encouraged to contact
   the devlopers through the XBraid Github page and issue tracker.  This area
   is active research.

   See ``examples/ex-01-expanded-bdf2.c``.

14. **Storage**: This option (see [braid_SetStorage](@ref braid_SetStorage))
    allows the user to specify storage at all time points (C and F) or only at
    C-points.  This extra storage is useful for implicit methods, where the
    solution value from the *previous XBraid iteration* for time step \f$i\f$
    can be used as the initial guess when computing step \f$i\f$ with the
    implicit solver.  This is often a better initial guess than using the
    solution value from the previous time step \f$i-1\f$.  The default is to
    store only C-point values, thus the better initial guess is only available
    at C-points in the default setting.  When storage is turned on at F-points,
    the better initial guess becomes available everywhere.
   
   In general, the user should always use the *ustop* parameter in
   *Step()* as the initial guess for an implicit solve.  If storage
   is turned on (i.e., set to 0), then this value will always be the improved
   initial guess for C- and F-points.  If storage is not turned on, then
   this will be the improved guess only for C-points.  For F-points,
   it will equal the solution from the previous time step.

   See ``examples/ex-03`` for an example which uses this feature.

15. **Delta Correction and Lyapunov Vector Estimation**: These options (see
    [braid_SetDeltaCorrection](@ref braid_SetDeltaCorrection) and
    [braid_SetLyapunovEstimation](@ref braid_SetLyapunovEstimation)) allow
    XBraid to accelerate convergence by using Delta correction, i.e. 
    low rank approximations to the Jacobian of the fine grid time-stepper 
    as a linear correction to the coarse grid time-stepper. This can 
    converge quadratically in some cases. LyapunovEstimation is not required
    for Delta correction, but for chaotic systems, the unstable modes of 
    error, corresponding with the first few Lyapunov vectors, are often the 
    slowest to converge, so Lyapunov estimation targets these modes by
    computing estimates to the backward Lyapunov vectors of the system, then
    computing the Delta correction using these vectors as a basis.

    See ``examples/ex-07`` for an example which uses these features.


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
     implemented in f90.

   - ``examples/ex-01-refinement.c``: same as ex-01.c, but adds the refinement feature
     of XBraid. The refinement can be arbitrary or based on error estimate.

# One-Dimensional Heat Equation {#exampletwo}

In this example, we assume familiarity with @ref exampleone.  This example is a
time-only parallel example that implements the 1D heat equation, 

\f[ \delta/\delta_t \; u(x,t) = \Delta\, u(x,t) + g(x,t), \f]

as opposed to @ref exampleone, which implements only a scalar ODE for one
degree-of-freedom in space.  There is no spatial parallelism, as a serial
cyclic reduction algorithm is used to invert the tri-diagonal spatial operators.
The space-time discretization is the standard 3-point finite
difference stencil (\f$[-1,2,-1]\f$), scaled by mesh widths.  Backward Euler
is used in time. 

This example consists of three files and two executables.

   - ``examples/ex-02-serial.c``:  This file compiles into its own executable
     ``ex-02-serial`` and represents a simple example user application that
     does sequential time-stepping.  This file represents where a new 
     XBraid user would start, in terms of converting a sequential time-stepping
     code to XBraid. 

   - ``examples/ex-02.c``:  This file compiles into its own executable
     ``ex-02`` and represents a time-parallel XBraid wrapping of the user
     application ``ex-02-serial``.  

   - ``ex-02-lib.c``:  This file contains shared functions used by the time-serial
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



# Simplest XBraid_Adjoint example {#exampleoneadjoint}

The file ``examples/ex-01-adjoint.c`` extends the simple scalar ODE example in ``ex-01.c`` for computing adjoint-based sensitivities.  See @ref exampleone.  The scalar ODE is 
\f[
   u_t(t) = \lambda u(t) \quad \forall t \in (0,T),
\f]
where \f$\lambda\f$ is considered the design variable. We consider an objective function of the form
\f[ 
   J(u,\lambda) = \int_0^T \frac{1}{T}  \| u(t) \|^2 dt.
\f] 


## User Defined Structures and Wrappers

The two user-defined data structures are:

1. **Vector:** This structure is unchanged from @ref exampleone, and contains
   a single scalar representing the state at a given time.

         typedef struct _braid_Vector_struct
         {
            double value;
         } my_Vector;


2. **App:**  This structure holds two additional elements when compared to @ref exampleone : 
   the *design* and the *reduced gradient*. This ensures that both are accessible in all user routines. 
         
         typedef struct _braid_App_struct
         {
            int       rank;
            double    design;
            double    gradient;
         } my_App;

The user must also define a few *additional* wrapper routines. Note, that the
app structure continues to be the first argument to every function.

1. All user-defined routines from `examples/ex-01.c` stay the same, except
   `Step()`, which must be changed to account for the new design parameter in `app`.

2. The user's **Step** routine queries the `app` to get the design and propagates
  the `braid_Vector u` forward in time for one time step:

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

            /* Get the design variable from the app */
            double lambda = app->design;

            /* Use backward Euler to propagate the solution */
            (u->value) = 1./(1. - lambda * (tstop-tstart))*(u->value);
            
            return 0;
         }  

3. **ObjectiveT**: This new routine evaluates the time-dependent part of the
   objective function at a local time \f$t_i\f$, i.e. it returns the integrand
   \f$f(u_i, \lambda) = \frac{1}{T} \|u_i \|_2^2 \f$. 

         int 
         my_ObjectiveT(braid_App              app,
                       braid_Vector           u,
                       braid_ObjectiveStatus  ostatus,
                       double                *objectiveT_ptr)
         {
            /* Get the total number of time steps */
            braid_ObjectiveStatusGetNTPoints(ostatus, &ntime);

            /* Evaluate the local objective: 1/N u(t)^2 */
            objT = 1. / ntime * (u->value) * (u->value);

            *objectiveT_ptr = objT;
            return 0;
         }

   The `ObjectiveStatus` can be queried for information about the current
   status of XBraid (e.g., what is the current time value, time-index, number
   of time steps, current iteration number, etc...).
   
   XBraid_Adjoint calls the `ObjectiveT` function on the finest time-grid level
   during the down-cycle of the multigrid algorithm and adds the value to a
   global objective function value with a simple summation.  Thus, any
   user-specific integration formula of the objective function must be here.

4. **ObjectiveT_diff**: This new routine updates the adjoint variable `u_bar`
   and the reduced gradient with the transposed partial derivatives of
   `ObjectiveT` multiplied by the scalar input \f$\bar F\f$, i.e.,
   \f[  
      \bar u_i = \frac{\partial f(u_i, \lambda)}{\partial u_i}^T \bar F \quad \text{and} \quad
      \bar \rho += \frac{\partial f(u_i, \lambda)}{\partial \rho}^T \bar F .
   \f]
   Note that \f$\bar u_i\f$ gets overwritten (''\f$=\f$''), whereas \f$\rho\f$ is updated (''\f$+=\f$'').

         int
         my_ObjectiveT_diff(braid_App            app,
                           braid_Vector          u,
                           braid_Vector          u_bar,
                           braid_Real            F_bar,
                           braid_ObjectiveStatus ostatus)
         {
            int    ntime;
            double ddu;      /* Derivative wrt u */
            double ddesign;  /* Derivative wrt design */

            /* Get the total number of time steps */
            braid_ObjectiveStatusGetNTPoints(ostatus, &ntime);

            /* Partial derivative with respect to u times F_bar */
            ddu = 2. / ntime * u->value * F_bar;

            /* Partial derivative with respect to design times F_bar*/
            ddesign = 0.0 * F_bar;

            /* Update u_bar and gradient */
            u_bar->value   = ddu;
            app->gradient += ddesign;

            return 0;
         }

5. **Step_diff**: This new routine computes transposed partial derivatives of the
   `Step` routine multiplied with the adjoint vector `u_bar` (\f$\bar{u}_i\f$), i.e., 
   \f[
      \bar u_{i} = \left(\frac{\partial \Phi_{i+1}(u_i, \rho)}{\partial u_i}\right)^T\bar u_i \quad \text{and} \quad
      \bar \rho += \left(\frac{\partial \Phi_{i+1}(u_i, \rho)}{\partial \rho}\right)^T\bar u_i .
   \f]

         int
         my_Step_diff(braid_App           app,
                      braid_Vector        ustop,
                      braid_Vector        u,
                      braid_Vector        ustop_bar,
                      braid_Vector        u_bar,
                      braid_StepStatus    status)
         {
            double ddu;      /* Derivative wrt u */
            double ddesign;  /* Derivative wrt design */

            /* Get the time step size */
            double tstop, tstart, deltat;
            braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
            deltat = tstop - tstart;

            /* Get the design from the app */
            double lambda = app->design;

            /* Transposed derivative of step wrt u times u_bar */
            ddu = 1./(1. - lambda * deltat) * (u_bar->value);
          
            /* Transposed derivative of step wrt design times u_bar */
            ddesign = (deltat * (u->value)) / pow(1. - deltat*lambda,2) * (u_bar->value);

            /* Update u_bar and gradient */
            u_bar->value      = ddu;     
            app->gradient    += ddesign;

            return 0;
         }

   **Important note on the usage of ustop**: If the `Step` routine uses the input vector `ustop` instead of `u` (typically for initializing a (non-)linear solve within \f$\Phi\f$), then `Step_diff` must update `ustop_bar` instead of `u_bar` and set `u_bar` to zero:
   \f[
      \overline{ustop}\,  += \left(\frac{\partial \Phi_{i+1}(ustop,\rho)}{\partial\, ustop}\right) ^T \bar u_i
      \quad \text{and} \quad \bar u_i = 0.0 .
   \f]

6. **ResetGradient**: This new routine sets the gradient to zero. 

         int 
         my_ResetGradient(braid_App app)
         {
            app->gradient = 0.0;
            return 0;
         }

   XBraid_Adjoint calls this routine before each iteration such that old gradient information is removed properly.


## Running XBraid_Adjoint for this example

The workflow for computing adjoint sensitivities with XBraid_Adjoint alongside
the primal state computation closely follows XBraid's workflow.  The user's
*main* file will first set up the `app` structure, holding the additional
information on an initial design and zero gradient.  Then, all the setup calls done in
@ref running_simplestexample will also be done. 

The XBraid_Adjoint specific calls are as follows.  After `braid_Init(...)` is
called, the user initializes XBraid_Adjoint by calling

         /* Initialize XBraid_Adjoint */
         braid_InitAdjoint( my_ObjectiveT, my_ObjectiveT_diff, my_Step_diff, my_ResetGradient, &core);

Next, in addition to the usual XBraid options for controlling the multigrid iterations, the adjoint solver's accuracy is set by calling 

      braid_SetAbsTolAdjoint(core, 1e-6);

After that, one call to 

      /* Run simulation and adjoint-based gradient computation */
      braid_Drive(core);

runs the multigrid iterations with additional adjoint sensitivity computations (i.e. the piggy-back iterations). After it finishes, the objective function value can be accessed by calling

      /* Get the objective function value from XBraid */
      braid_GetObjective(core, &objective);

Further, the reduced gradient, which is stored in the user's `App` structure, holds the sensitivity information \f$ dJ/d\rho \f$. As this information is local to all the time-processors, the user is responsible for summing up the gradients from all time-processors, if necessary. This usually involves an `MPI_Allreduce` call as in 

      /* Collect sensitivities from all processors */
      double mygradient = app->gradient;
      MPI_Allreduce(&mygradient, &(app->gradient), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 

Lastly, the gradient computed with XBraid_Adjoint is verified using Finite Differences.  See the source code ``examples/ex-01-adjoint.c`` for details. 


# Optimization with the Simplest Example {#exampleoneoptimization}

The file ``examples/ex-01-optimization.c`` implements a simple optimization iteration by extending ``examples/ex-01-adjoint.c``, described in @ref exampleoneadjoint.
This example solves an inverse design problem for the simple scalar ODE example:

\f[
  \begin{matrix} \min ~ \frac 1 2 \left( \int_0^T \frac{1}{T}  \| u(t) \|^2 dt - J_{\text{Target}} \right)^2 + \frac{\gamma}{2} \| \lambda \|^2 \\ \\
  \text{s.t. } \quad \frac{\partial}{\partial t}u(t) = \lambda u(t) \quad \forall t \in (0,T) \quad
  \end{matrix}
\f]
where \f$ J_{\text{Target}} \f$ is a fixed and precomputed target value and \f$\gamma >0 \f$ is a fixed relaxation parameter.  Those fixed values are stored within the `App`. 

## User Defined Structures and Wrappers

In order to evaluate the time-independent part of the objective function (e.g. the postprocessing function \f$F\f$) and its derivative, two additional user routines are necessary. *There are no new user-defined data structures.*

1. **PostprocessObjective**: This function evaluates the tracking-type objective function and the regularization term. The input variable `integral` contains the integral-part of the objective and returns the objective that is to be minimized \f$ F(I) \f$: 


         /* Evaluate the time-independent part of the objective function */
         int
         my_PostprocessObjective(braid_App  app,
                              double        integral,
                              double       *postprocess
                              )
         {
            double F;

            /* Tracking-type functional */
            F  = 1./2. * pow(integral - app->target,2);
            
            /* Regularization term */
            F += (app->gamma) / 2. * pow(app->design,2);

            *postprocess = F;
             return 0;
         }

2. **PostprocessObjective_diff**: This provides XBraid_Adjoint with the partial derivatives of the `PostprocessObjective` routine, i.e.
\f[
   \bar F    = \frac{\partial F(I, \lambda)}{\partial I} \quad \text{and} \quad 
   \bar \rho += \frac{\partial F(I,\lambda)}{\partial \lambda} 
\f]

         int
         my_PostprocessObjective_diff(braid_App   app,
                                      double      integral,
                                      double     *F_bar
                                      )
         {

            /* Derivative of tracking type function */
            *F_bar = integral - app->target;

            /* Derivative of regularization term */
            app->gradient += (app->gamma) * (app->design);
            return 0;
         }

These routines are optional for XBraid_Adjoint. Therefore, they need to be passed to XBraid_Adjoint after the initialization with `braid_Init(...)` and `braid_InitAdjoint(...)` in the user's *main* file:

         /* Optional: Set the tracking type objective function and derivative */
         braid_SetPostprocessObjective(core, my_PostprocessObjective);
         braid_SetPostprocessObjective_diff(core, my_PostprocessObjective_diff);


## Running an Optimization Cycle with XBraid_Adjoint 

XBraid_Adjoint does not natively implement any optimization algorithms.
Instead, we provide examples showing how one can easily use XBraid_Adjoint
inside an optimization cycle.  Here, one iteration of the optimization 
cycle consists of the following steps:

1. First, we run XBraid_Adjoint to solve the primal and adjoint dynamics:

         braid_Drive(core);

2. Get the value of the objective function with 

         braid_GetObjective(core, &objective);

3. Gradient information is stored in the `app` structure. Since it is local to all temporal processors, we need to invoke an `MPI_Allreduce` call which sums up the local sensitivities:
      
         mygradient = app->gradient;
         MPI_Allreduce(&mygradient, &app->gradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   **Note**: For time-dependent design variables, summing over all processors might not be necessary, since information is needed only locally in time. See `examples/ex-04.c` for a time-dependent design example. 

4. Update the design variable using the gradient information. Here, we implement a simple steepest descent update into the direction of the negative gradient:

         app->design -= stepsize * app->gradient;

   Here, a fixed step size is used to update the design variable. Usually, a
   line-search procedure should be implemented in order to find a suitable step
   length that minimizes the objective function along the update direction.
   However to carry out a line search, we must re-evaluate the objective
   function for different design value(s).  Thus, the option
   [braid_SetObjectiveOnly(core, 1)](@ref braid_SetObjectiveOnly) can be used.
   After this option has been set, any further call to `braid_Drive(core)` will
   then only run a primal XBraid simulation and carry out an objective function
   evaluation.  No gradients will be computed, which saves computational time.
   After the line search, make sure to reset XBraid_Adjoint for gradient
   computation with `braid_SetObjectiveOnly(core, 0)`. 

5. The optimization iterations are stopped when the norm of the gradient is below a prescribed tolerance. 


# A Simple Optimal Control Problem {#optimalcontrolexample}

This example demonstrates the use of XBraid_Adjoint for solving an optimal control problem with time-dependent design variables:
\f[
   \begin{matrix}
   \min ~  \int_0^1 u_1(t)^2 + u_2(t)^2 + \gamma c(t)^2 ~ dt & \\ & \\
       \text{s.t.} \quad  \frac{\partial}{\partial t} u_1(t) = u_2(t)      &  \forall t\in(0,1) \\
              \qquad \qquad \qquad \frac{\partial}{\partial t} u_2(t) = -u_2(t) + c(t) & \forall t\in(0,1)
   \end{matrix}
\f]
with initial condition \f$u_1(0) = 0, u_2(0) = -1\f$ and piecewise constant control (design) variable \f$c(t)\f$.  

The example consists of three files, meant to indicate how one can take a
time-serial implementation for an optimal control problem and create a
corresponding XBraid_Adjoint implementation.
- `examples/ex-04-serial.c`: Compiles into its own executable
  `examples/ex-04-serial`, which solves the optimal control problem using
  time-serial forward-propagation of state variables and time-serial
  backward-propagation of the adjoint variables in each iteration of an outer
  optimization cycle.
- `examples/ex-04.c`: Compiles into `ex-04`.  This solves the same optimization
  problem in time-parallel by replacing the forward- and backward-propagation
  of state and adjoint by the time-parallel XBraid and XBraid_Adjoint solvers. 
- `examples/ex-04-lib.c`: Contains the routines that are shared by both the
  serial and the time-parallel implementation.  Study this file, and discover
  that most of the important code setting up the user-defined data structures and
  wrapper routines are simply lifted from the serial simulation.

# Chaotic Lorenz System With Delta Correction and Lyapunov Estimation {#exampleseven}

This example demonstrates acceleration of XBraid convergence and Lyapunov analysis of a system with Delta correction. Familiarity with @ref exampleone is assumed. This example solves the chaotic Lorenz system in three dimensions, defined by the system
\f[
   \begin{cases}
      x' = \sigma (y - x),   \\
      y' = x (\rho - z) - y, \\
      z' = xy - \beta z,
   \end{cases}
\f]
where \f$\sigma = 10\f$, \f$\rho = 28\f$, and \f$\beta = 8/3\f$. This system is chaotic, with the greatest Lyapunov exponent being \f$\approx 0.9\f$. Here, Delta correction is used to accelerate convergence to the solution, while Lyapunov estimation is used to simultaneously compute the Lyapunov vectors and Lyapunov exponents along the trajectory.

## User Defined Structures and Wrappers

Most of the user defined structures and wrappers are defined exactly as in previous examples, with the exception of *Step()*, *BufSize()*, and *Access()*, which are modified to accomodate the Lyapunov vectors, and *InnerProd()* and *InitBasis()*, which are new functions required by Delta correction.

1. **Step**: Here the *Step* function is required to do two things:
   1. Propagate the state vector (as in regular XBraid)
      \f[
            u \gets \Phi(u)
      \f]
   2. Propagate a number of basis vectors using the Jacobian vector product (new functionality required by Delta correction)
      \f[
            \psi_j \gets \left(\frac{d \Phi}{d u}\right) \psi_j
      \f]
   The number of basis vectors to be propagated is accessed via [braid_StepStatusGetDeltaRank](@ref braid_StatusGetDeltaRank), and references to the vectors themselves are accessed via [braid_StepStatusGetBasisVec](@ref braid_StatusGetBasisVec). In this example, the full Jacobian of *Step* is used to propagate the basis vectors, but finite differencing or even forward-mode automatic differentiation are other ways of propagating the basis vectors.

         int my_Step(braid_App app,
                     braid_Vector ustop,
                     braid_Vector fstop,
                     braid_Vector u,
                     braid_StepStatus status)
         {
            /* for Delta correction, the user must propagate the solution vector (as in a traditional Braid code)
            * as well as the Lyapunov vectors. The Lyapunov vectors are available through the StepStatus structure,
            * and are propagated by the Jacobian of the time-step function. (see below)
            */
            double tstart; /* current time */
            double tstop;  /* evolve to this time */
            braid_StepStatusGetTstartTstop(status, &tstart, &tstop);

            double h; /* dt value */
            h = tstop - tstart;

            // get the number of Lyapunov vectors we need to propagate
            int rank; /* rank of Delta correction */
            braid_StepStatusGetDeltaRank(status, &rank);
            MAT Jacobian = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

            if (rank > 0) // we are propagating Lyapunov vectors
            {
               Euler((u->values), h, &Jacobian);
            }
            else
            {
               Euler((u->values), h, NULL);
            }

            for (int i = 0; i < rank; i++)
            {
               // get a reference to the ith Lyapunov vector
               my_Vector *psi;
               braid_StepStatusGetBasisVec(status, &psi, i);

               // propagate the vector from tstart to tstop
               if (psi)
               {
                  MatVec(Jacobian, psi->values);
               }
            }

            /* no refinement */
            braid_StepStatusSetRFactor(status, 1);

            return 0;
         }

2. **BufSize()**: There is an additional option to set the size of a single basis 
   vector here, via [braid_BufferStatusSetBasisSize](@ref braid_StatusSetBasisSize).

         int my_BufSize(braid_App app, int *size_ptr, braid_BufferStatus bstatus)
         {
            /* Tell Braid the size of a state vector */
            *size_ptr = VecSize * sizeof(double);

            /* 
            * In contrast with traditional Braid, you may also specify the size of a single Lyapunov basis vector, 
            * in case it is different from the size of a state vector.
            * Note: this isn't necessary here, but for more complicated applications this size may be different.
            */
            braid_BufferStatusSetBasisSize(bstatus, VecSize * sizeof(double));
            return 0;
         }

3. **Access**: Here, the *Access* function is used to access the Lyapunov
   vector estimates via the same api as for *Step*. Also, the local Lyapunov
   exponents are accessed via [braid_AccessStatusGetLocalLyapExponents](@ref braid_StatusGetLocalLyapExponents).

         int my_Access(braid_App app, braid_Vector u, braid_AccessStatus astatus)
         {
            FILE *file = (app->file);
            int index, i;
            double t;

            braid_AccessStatusGetT(astatus, &t);
            braid_AccessStatusGetTIndex(astatus, &index);

            fprintf(file, "%d", index);
            for (i = 0; i < VecSize; i++)
            {
               fprintf(file, " %.14e", (u->values[i]));
            }
            fprintf(file, "\n");
            fflush(file);

            /* write the lyapunov vectors to file */
            file = app->file_lv;
            int local_rank, num_exp;
            braid_AccessStatusGetDeltaRank(astatus, &local_rank);
            num_exp = local_rank;
            double *exponents = malloc(local_rank * sizeof(double));
            if (num_exp > 0)
            {
               braid_AccessStatusGetLocalLyapExponents(astatus, exponents, &num_exp);
            }

            fprintf(file, "%d", index);
            for (int j = 0; j < local_rank; j++)
            {
               my_Vector *psi;
               braid_AccessStatusGetBasisVec(astatus, &psi, j);
               if (psi)
               {
                  if (j < num_exp)
                  {
                     (app->lyap_exps)[j] += exponents[j];
                     fprintf(file, " %.14e", exponents[j]);
                  }
                  else
                  {
                     fprintf(file, " %.14e", 0.);
                  }
                  for (i = 0; i < VecSize; i++)
                  {
                     fprintf(file, " %.14e", (psi->values[i]));
                  }
               }
            }
            fprintf(file, "\n ");
            fflush(file);
            free(exponents);

            return 0;
         }

4. **InnerProd**: This function tells XBraid how to compute the inner product between two *Vector*s.
   This is required by Delta correction in order to project user vectors onto the basis vectors, and
   for orthonormalization of the basis vectors. Here, the standard dot product is used.

         int my_InnerProd(braid_App app, braid_Vector u, braid_Vector v, double *prod_ptr)
         {
           /*
            *  For Delta correction, braid needs to be able to compute an inner product between two user vectors,
            *  which is used to project the user's vector onto the Lyapunov basis for low-rank Delta correction.
            *  This function should define a valid inner product between the vectors *u* and *v*.
            */
            double dot = 0.;

            for (int i = 0; i < VecSize; i++)
            {
               dot += (u->values[i]) * (v->values[i]);
            }
            *prod_ptr = dot;
            return 0;
         }

5. **InitBasis**: This function tells XBraid how to initialize a single basis vector, with spatial
   index *j* at time *t*. This initializes the *j*th column of the matrix \f$\Psi\f$ whose
   columns are the basis vectors used for Delta correction. Here, we simply use the *j*th column
   of the identity matrix. It is important that the vectors initialized by this function are
   linearly independent, or Lyapunov estimation will not work.

         int my_InitBasis(braid_App app, double t, int index, braid_Vector *u_ptr)
         {
            /*
            *  For Delta correction, an initial guess is needed for the Lyapunov basis vectors.
            *  This function initializes the basis vector with spatial index *index* at time *t*.
            *  Note that the vectors at each index *index* must be linearly independent.
            */
            my_Vector *u;

            u = (my_Vector *)malloc(sizeof(my_Vector));

            // initialize with the columns of the identity matrix
            VecSet(u->values, 0.);
            u->values[index] = 1.;

            *u_ptr = u;

            return 0;
         }

## Running XBraid with Delta correction and Lyapunov Estimation

XBraid is initialized as before, and most XBraid features are compatible, not including Richardson extrapolation, the XBraid_Adjoint feature, the *Residual* option, and spatial coarsening. Delta correction and Lyapunov estimation are turned on by calls to [braid_SetDeltaCorrection](@ref braid_SetDeltaCorrection) and [braid_SetLyapunovEstimation](@ref braid_SetDeltaCorrection), respectively, where the number of basis vectors desired (rank of low-rank Delta correction) and additional wrapper functions *InnerProd* and *InitBasis* are passed to XBraid and options regarding the estimation of Lyapunov vectors and exponents are set. Further, the function [braid_SetDeferDelta](@ref braid_SetDeferDelta) gives more options allowing Delta correction to be deferred to a later iteration, or a coarser grid. This is illustrated in the folowing exerpt from the *main()* function:

         ...
         if (delta_rank > 0)
         {
            braid_SetDeltaCorrection(core, delta_rank, my_InitBasis, my_InnerProd);
            braid_SetDeferDelta(core, defer_lvl, defer_iter);
            braid_SetLyapunovEstimation(core, relax_lyap, lyap, relax_lyap || lyap);
         }
         ...

# Running and Testing XBraid {#testingxbraid}

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

    /**
     * Test innerprod(), initbasis(), step(), bufsize(), bufpack(), bufunpack() 
     * for use with Delta correction
     */
    correct = braid_TestInnerProd(app, comm_x, stdout, 0.0, 1.0,
                                  my_Init, my_Free, my_Sum, my_InnerProd);
    correct = braid_TestDelta(app, comm_x, stdout, 0.0, dt, delta_rank, my_Init,
                              my_InitBasis, my_Access, my_Free, my_Sum, my_BufSize,
                              my_BufPack, my_BufUnpack, my_InnerProd, my_Step);

# Fortan90 Interface, C++ Interface, Python Interface, and More Complicated Examples {#complicatedexamples}

We have Fortran90, C++, and Python interfaces.  For Fortran 90, see ``examples/ex-01f.f90``.
For C++ see ``braid.hpp`` and ``examples/ex-01-pp.cpp``
For more complicated C++ examples, see the various C++ examples in
``drivers/drive-**.cpp``.  For Python, see the directories ``examples/ex-01-cython`` 
and ``examples/ex-01-cython-alt``.

For a discussion of more complex problems please see our project
[publications website](http://computation.llnl.gov/projects/parallel-time-integration-multigrid/publications)
for our recent publications concerning some of these varied applications. 
