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


# Meaning of the name {#braidname}

We chose the package name XBraid to stand for _Time-Braid_, where
X is the first letter in the Greek word for time, _Chronos_.  The
algorithm _braids_ together time-grids of different granularity in order to
create a multigrid method and achieve parallelism in the time dimension.

# Advice to users {#advice}

The field of parallel-in-time methods is in many ways under development, and
success has been shown primarily for problems with some parabolic character.
While there are ongoing projects (here and elsewhere) looking at varied
applications such as hyperbolic problems, computational fluid dynamics, power
grids, medical applications, and so on, expectations should take this fact into
account.  Please see our project
[publications website](http://computation.llnl.gov/projects/parallel-time-integration-multigrid/publications)
for our recent publications concerning some of these varied applications. 

That being said, we strongly encourage new users to try our code for
their application.  Every new application has its own issues to address and
this will help us to improve both the algorithm and the software.

For support, please email `xbraid-support@llnl.gov`.  This email address
automically interfaces with our issue tracker and notifies all developers of the
pending support request.

# Overview of the XBraid Algorithm {#braidoverview}

The goal of XBraid is to solve a problem faster than a traditional time
marching algorithm.  Instead of sequential time marching, XBraid solves the
problem iteratively by simultaneously updating a space-time solution guess over
all time values.  The initial solution guess can be anything, even a random
function over space-time.  The iterative updates to the solution guess are done
by constructing a hierarchy of temporal grids, where the finest grid contains
all of the time values for the simulation.  Each subsequent grid is a coarser
grid with fewer time values.  The coarsest grid has a trivial number of time
steps and can be quickly solved exactly.  The effect is that solutions to the
time marching problem on the coarser (i.e., cheaper) grids can be used to
correct the original finest grid solution.  Analogous to spatial multigrid, the
coarse grid correction only *corrects* and *accelerates* convergence to the
finest grid solution.  The coarse grid does not need to represent an
accurate time discretization in its own right.  Thus, a problem with many time
steps (thousands, tens of thousands or more) can be solved with 10 or 15 XBraid
iterations, and the overall time to solution can be greatly sped up.  However,
this is achieved at the cost of more computational resources.

To understand how XBraid differs from traditional time marching, consider the
simple linear advection equation, \f$ u_t = -c u_x \f$.  The next figure depicts how one
would typically evolve a solution here with sequential time stepping.  The initial condition is  
a wave, and this wave propagates sequentially across space as time increases.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.25\textwidth]{../img/sequential0.pdf}}
       \subfloat{\includegraphics[width=0.25\textwidth]{../img/sequential1.pdf}} \\
       \subfloat{\includegraphics[width=0.25\textwidth]{../img/sequential2.pdf}}
       \subfloat{\includegraphics[width=0.25\textwidth]{../img/sequential3.pdf}} 
       \caption{Sequential time stepping.}
       \label{img:seqential}
   \end{figure}
   \endlatexonly

XBraid instead begins with a solution guess over all of space-time, which for demonstration, 
we let be random. An XBraid iteration does

1. Relaxation on the fine grid, i.e., the grid that contains all of the desired time values.
   Relaxation is just a local application of the time stepping scheme, e.g., backward Euler.
2. Restriction to the first coarse grid, i.e., interpolate the problem to a grid that 
   contains fewer time values, say every second or every third time value.
3. Relaxation on the first coarse grid
4. Restriction to the second coarse grid and so on...
5. When a coarse grid of trivial size (say 2 time steps) is reached, it is solved exactly.
6. The solution is then interpolated from the coarsest grid to the finest grid

One XBraid iteration is called a *cycle* and these cycles continue until the
solution is accurate enough.  This is depicted in the next figure, where only a
few iterations are required for this simple problem.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=1.0\textwidth]{../img/3_MG_In_Time_Iterations.pdf}}
       \caption{XBraid iterations.}
       \label{img:mgrit_cycles}
   \end{figure}
   \endlatexonly

There are a few important points to make.
- The coarse time grids allow for global propagation of information across
  space-time with only one XBraid iteration.  This is visible in the above figure by observing
  how the solution is updated from iteration 0 to iteration 1.
- Using coarser (cheaper) grids to correct the fine grid is analogous to spatial multigrid.
- Only a few XBraid iterations are required to find the solution over 1024 time steps. 
  Therefore if enough processors are available to parallelize XBraid, we can see a speedup
  over traditional time stepping (more on this later).
- This is a simple example, with evenly space time steps.  XBraid is structured 
  to handle variable time step sizes and adaptive time step sizes. 

To firm up our understanding, let`s do a little math.  Assume that you
have a general system of ordinary differential equations (ODEs), 
   \f[
   \Xux^{\prime}(t) = \Xfx(t, \Xux(t)), ~~~ \Xux(0) = \Xux_0, ~~~ t \in [0,T]. 
   \f]
Next, let \f$t_i = i \delta t, i = 0, 1, ..., N\f$ be a temporal mesh with spacing 
\f$ \delta t = T/N \f$, and \f$u_i\f$ be an approximation to \f$u(t_i)\f$.
A general one-step time discretization is now given by
   \f[
   \begin{split}
      \Xux_0 =& g_0\\ 
      \Xux_i =& \Phi_i(\Xux_{i-1}) + \Xgx_i , ~~~ i = 1, 2, ... , N.
   \end{split}
   \f]

Traditional time marching would first solve for \f$ i = 1 \f$, then solve for 
\f$ i = 2 \f$, and so on.  For linear time propagators \f$\{\Phi_i\}\f$, this 
can also be expressed as applying a direct solver (a forward solve) to the following system:
   \f[
   A \Xu \equiv
   \begin{pmatrix}
    I & & & \\
   -\Phi_1 & I & & \\
    & \ddots & \ddots & \\
    & & -\Phi_N & I
   \end{pmatrix}
   \begin{pmatrix}
   \Xux_0 \\
   \Xux_1 \\
   \vdots \\
   \Xux_N
   \end{pmatrix} =
   \begin{pmatrix}
   \Xgx_0 \\
   \Xgx_1 \\
   \vdots \\
   \Xgx_N
   \end{pmatrix} \equiv \Xg 
   \f]
or 
   \f[
   A \mathbf{u} = \mathbf{g}.
   \f]
This process is optimal and O(N), but it is sequential.  XBraid achieves
parallelism in time by replacing this sequential solve with an optimal multigrid
reduction iterative method
   \latexonly
   \footnote{ Ries, Manfred, Ulrich Trottenberg, and Gerd Winter. "A note on MGR methods." Linear Algebra and its Applications 49 (1983): 1-26.}
   \endlatexonly
applied to only the time dimension. This approach is
- nonintrusive, in that it coarsens only in time and the user defines \f$ \Phi \f$. 
  Thus, users can continue using existing time stepping codes by wrapping them
  into our framework.
- optimal and O(N), but O(N) with a higher constant than time stepping.
  Thus with enough computational resources, XBraid will outperform sequential time stepping.
- highly parallel

We now describe the two-grid process in more detail, with the multilevel analogue being a recursive
application of the process.  We also assume that \f$ \Phi \f$ is constant
for notational simplicity.  XBraid coarsens in the time
dimension with factor \f$m > 1\f$ to yield a coarse time grid with 
\f$N_\Delta = N/m\f$ points and time step \f$\Delta T = m \delta t\f$.  
The corresponding coarse grid problem,
   \f[
   A_{\Delta} =
   \begin{pmatrix}
    I & & & \\
   -\Phi_{\Delta} & I & & \\
    & \ddots & \ddots & \\
    & & -\Phi_{\Delta} & I
   \end{pmatrix},
   \f]
is obtained by defining coarse grid propagators \f$\{\Phi_{\Delta}\}\f$ which
are at least as cheap to apply as the fine scale propagators \f$\{\Phi\}\f$.
The matrix \f$ A_{\Delta} \f$ has fewer rows and columns than \f$ A \f$, e.g.,
if we are coarsening in time by 2, \f$ A_{\Delta} \f$ has one half as many rows
and columns.

This coarse time grid induces a partition of the fine grid into C-points
(associated with coarse grid points) and F-points, as visualized next. 
C-points exist on both the fine and coarse time grid, but F-points exist only 
on the fine time scale.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.75\textwidth]{../img/timeline.pdf}}
       \label{img:timeline}
   \end{figure}
   \endlatexonly

Every multigrid algorithm requires a relaxation method and an approach to
transfer values between grids.  Our relaxation scheme alternates between
so-called F-relaxation and C-relaxation as illustrated next.  F-relaxation
updates the F-point values \f$\{u_j\}\f$ on interval \f$(T_i, T_{i+1})\f$ by
simply propagating the C-point value \f$u_{mi}\f$ across the interval using the
time propagator \f$\{\Phi\}\f$.  While this is a sequential
process, each F-point interval update is independent from the others and can be
computed in parallel.  Similarly, C-relaxation updates the C-point value
\f$u_{mi}\f$ based on the F-point value \f$u_{mi-1}\f$ and these updates can
also be computed in parallel.  This approach to relaxation can be thought of as
line relaxation in space in that the residual is set to 0 for an entire time
step. 

The F updates are done simultaneously in parallel,
as depicted next.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.75\textwidth]{../img/FrelaxationDetail.pdf}}
       \label{img:FrelaxDetail}
       \caption{Update all F-point intervals in parallel, using the time propagator $\Phi$.}
   \end{figure}
   \endlatexonly

Following the F sweep, the C updates are also done simultaneously in parallel, 
as depicted next.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.75\textwidth]{../img/CrelaxationDetail.pdf}}
       \label{img:CrelaxDetail}
       \caption{Update all C-points in parallel, using the time propagator $\Phi$.}
   \end{figure}
   \endlatexonly

\latexonly \newpage \endlatexonly
In general, FCF- and F-relaxation will refer to the relaxation methods used in XBraid. We can say
- FCF- or F-relaxation is highly parallel.
- But, a sequential component exists equaling the number of F-points between two C-points.
- XBraid uses regular coarsening factors, i.e., the spacing of C-points happens every \f$m\f$ points.

After relaxation, comes forming the coarse grid error correction.  To move
quantities to the coarse grid, we use the restriction operator \f$ R \f$ which
simply injects values at C-points from the fine grid to the coarse grid, 
   \f[
   R =
   \begin{pmatrix}
    I          &            &     \\
    0          &            &     \\
    \vdots     &            &     \\
    0          &            &     \\
               & I          &     \\
               & 0          &     \\
               & \vdots     &     \\
               & 0          &     \\
               &            & \ddots
   \end{pmatrix}^T.
   \f]
The spacing between each \f$ I \f$ is \f$ m-1 \f$ block rows.  While injection
is simple, XBraid always does an F-relaxation sweep before the application of \f$R\f$, 
which is equivalent to using the transpose of harmonic interpolation for restriction 
(see [Parallel Time Integration with Multigrid](https://computation.llnl.gov/project/linear_solvers/pubs/mgritPaper-2014.pdf) ).
Another interpretation is that the F-relaxation compresses the residual into the 
C-points, i.e., the residual at all F-points after an F-relaxation is 0.  Thus, 
it makes sense for restriction to be injection.

To define the coarse grid equations, we apply the Full Approximation
Scheme (FAS) method, which is a nonlinear version of multigrid.  This is to
accommodate the general case where \f$ f \f$ is a nonlinear function.  In FAS,
the solution guess and residual (i.e., \f$ \mathbf{u}, \mathbf{g} - A
\mathbf{u}\f$) are restricted.  This is in contrast to linear multigrid which
typically restricts only the residual equation to the coarse grid.  This
algorithmic change allows for the solution of general nonlinear problems.  For more
details, see this
[PDF](http://computation.llnl.gov/casc/people/henson/postscript/UCRL_JC_150259.pdf)
by Van Henson for a good introduction to FAS.  However, FAS was originally invented 
by Achi Brandt.

A central question in applying FAS is how to form the coarse grid matrix
\f$A_{\Delta}\f$, which in turn asks how to define the coarse grid time stepper \f$
\Phi_{\Delta} \f$.  One of the simplest choices (and one frequently used in
practice) is to let \f$ \Phi_{\Delta} \f$ simply be \f$\Phi\f$ but with the
coarse time step size \f$ \Delta T = m \delta t \f$.  For example, if \f$ \Phi
= (I - \delta t A)^{-1} \f$ for some backward Euler scheme, then \f$
\Phi_{\Delta} = (I - m \delta t A)^{-1} \f$ would be one choice.

With this \f$ \Phi_{\Delta} \f$ and letting \f$ \mathbf{u}_{\Delta} \f$ be the 
restricted fine grid solution and \f$ \mathbf{r}_{\Delta} \f$ be the restricted
fine grid residual, the coarse grid equation
   \f[
   A_{\Delta}(\mathbf{v}_{\Delta})  = A_{\Delta}(\mathbf{u}_{\Delta}) + \mathbf{r}_{\Delta} 
   \f]
is then solved. Finally, FAS defines a coarse grid error approximation 
\f$ \mathbf{e}_{\Delta} = \mathbf{v}_{\Delta} - \mathbf{u}_{\Delta} \f$,
which is interpolated with \f$ P_{\Phi} \f$ back to the fine grid and added to 
the current solution guess.    Interpolation is equivalent to injecting the coarse grid to 
the C-points on the fine grid, followed by an F-relaxation sweep (i.e., it is equivalent
to harmonic interpolation, as mentioned above about restriction).  That is,  
   \f[
   P_{\Phi} =
   \begin{pmatrix}
    I          &            &     \\
    \Phi       &            &     \\
    \Phi^2     &            &     \\
    \vdots     &            &     \\
    \Phi^{m-1} &            &     \\
               & I          &     \\
               & \Phi       &     \\
               & \Phi^2     &     \\
               & \vdots     &     \\
               & \Phi^{m-1} &     \\
               &            & \ddots
   \end{pmatrix},
   \f]
where \f$ m \f$ is the coarsening factor.  See @ref twogrid for a concise description of the 
FAS algorithm for MGRIT.

## Two-Grid Algorithm {#twogrid}

The two-grid FAS process is captured with this algorithm.  Using a recursive
coarse grid solve (i.e., step 3 becomes a recursive call) makes the process
multilevel.  Halting is done based on a residual tolerance.  If the operator is
linear, this FAS cycle is equivalent to standard linear multigrid.  Note that
we represent \f$ A \f$ as a function below, whereas the above notation was
simplified for the linear case.
   1. Relax on \f$A(\mathbf{u}) = \mathbf{g}\f$ using FCF-relaxation

   2. Restrict the fine grid approximation and its residual:
   \f[ \mathbf{u}_{\Delta} \leftarrow  R \mathbf{u}, \quad \mathbf{r}_{\Delta} \leftarrow  R (\mathbf{g} - A(\mathbf{u}), \f] 
   which is equivalent to updating each individual time step according to
   \f[
   u_{\Delta,i} \leftarrow u_{mi},\quad
   r_{\Delta,i} \leftarrow g_{mi} - A (\mathbf{u})_{mi}
   \quad \text{for}\quad i = 0, ..., N_\Delta.
   \f]

   3. Solve \f$A_{\Delta}(\mathbf{v}_{\Delta})  = A_{\Delta}(\mathbf{u}_{\Delta}) + \mathbf{r}_{\Delta} \f$

   4. Compute the coarse grid error approximation: \f$ \mathbf{e}_{\Delta} = \mathbf{v}_{\Delta} - \mathbf{u}_{\Delta} \f$

   5. Correct: \f$ \mathbf{u} \leftarrow \mathbf{u} + P \mathbf{e}_{\Delta} \f$

      This is equivalent to updating each individual time step by adding 
      the error to the values of \f$\mathbf{u}\f$ at the C-points:
      \f[ u_{mi} = u_{mi} + e_{\Delta,i}, \f]
      followed by an F-relaxation sweep applied to \f$\mathbf{u}\f$.


## Summary {#algorithmsummary}

In summary, a few points are
- XBraid is an iterative solver for the global space-time problem.
- The user defines the time stepping routine \f$ \Phi \f$ and can wrap 
  existing code to accomplish this. 
- XBraid convergence will depend heavily on how well \f$ \Phi_{\Delta} \f$ approximates
  \f$ \Phi^m \f$, that is how well a time step size of \f$ m \delta t = \Delta T \f$
  will approximate \f$ m \f$ applications of the same time integrator for a time step 
  size of \f$ \delta t \f$.  This is a subject of research, but this approximation
  need not capture fine scale behavior, which is instead captured by relaxation on 
  the fine grid.
- The coarsest grid is solved exactly, i.e., sequentially, which can be a bottleneck
  for two-level methods like Parareal,
  \latexonly
   \footnote{ Lions, J., Yvon Maday, and Gabriel Turinici. "A''parareal''in time discretization of PDE's." Comptes Rendus de l'Academie des Sciences Series I Mathematics 332.7 (2001): 661-668.}
  \endlatexonly
  but not for a multilevel scheme like XBraid where the coarsest grid is of trivial size.
- By forming the coarse grid to have the same sparsity structure and time stepper
  as the fine grid, the algorithm can recur easily and efficiently.
- Interpolation is ideal or exact, in that an application of 
  interpolation leaves a zero residual at all F-points.
- The process is applied recursively until a trivially sized temporal grid is
  reached, e.g., 2 or 3 time points.  Thus, the coarsening rate \f$ m \f$
  determines how many levels there are in the hierarchy.  For instance in this
  figure, a 3 level hierarchy is shown.  Three levels are chosen because there are
  six time points, \f$m = 2\f$ and \f$ m^2 < 6 \le m^3 \f$.
  If the coarsening rate had been \f$m = 4\f$ then there would only be two
  levels because there would be no more points to coarsen!
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.6\textwidth]{../img/3_levels.pdf}}
   \end{figure}
   \endlatexonly
   
   By default, XBraid will subdivide the time domain into evenly sized time
   steps.  XBraid is structured to handle variable time step sizes and adaptive
   time step sizes. 

# Overview of the XBraid Code {#codeoverview}

XBraid is designed to run in conjunction with an existing application code that
can be wrapped per our interface.  This application code will implement some
time marching simulation like fluid flow.  Essentially, the user has to
take their application code and extract a stand-alone time-stepping function
\f$ \Phi \f$ that can evolve a solution from one time value to another,
regardless of time step size.  After this is done, the XBraid code takes care of
the parallelism in the time dimension.

XBraid 
- is written in C and can easily interface with Fortran and C++
- uses MPI for parallelism
- self documents through comments in the source code and through *.md files
- functions and structures are prefixed by *braid* 
  + User routines are prefixed by `braid_`
  + Developer routines are prefixed by `_braid_` 


## Parallel decomposition and memory {#decomposition}
- XBraid decomposes the problem in parallel as depicted next. 
  \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.75\textwidth]{../img/data_layout.pdf}}
       \label{img:data_layout}
   \end{figure}
   \endlatexonly
   As you can see, traditional time stepping only stores one time step at a time, 
   but only enjoys a spatial data decomposition and spatial parallelism.  On the other
   hand, XBraid stores multiple time steps simultaneously and each processor holds a space-time
   chunk reflecting both the spatial and temporal parallelism.

- XBraid only handles temporal parallelism and is agnostic to the spatial decomposition.  
  See [braid_SplitCommworld](@ref braid_SplitCommworld).  
  Each processor owns a certain number of CF intervals of points. In the following figure, 
  processor 1 and processor 2 each own 2 CF intervals.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.75\textwidth]{../img/parallel_timeline.pdf}}
   \end{figure}
   \endlatexonly
   XBraid distributes intervals evenly on the finest grid.

- XBraid increases the parallelism significantly, but now several time steps need 
  to be stored, requiring more memory.  XBraid employs two strategies to address 
  the increased memory costs.  
  
  + First, one need not solve the whole problem at once.  Storing only one
    space-time slab is advisable.  That is, solve for as many time steps
    (say _k_ time steps) as you have available memory for.  Then move on to
    the next _k_ time steps.

  + Second, XBraid provides support for storing only C-points.  Whenever an
    F-point is needed, it is generated by F-relaxation.  More precisely, only
    the red C-point time values in the previous figure are stored.  Coarsening
    is usually aggressive with \f$ m = 8, 16, 32, ... \f$, so the storage
    requirements of XBraid are significantly reduced when compared to storing
    all of the time values.

  Overall, the memory multiplier per processor when using XBraid is \f$ O(1)
  \f$ if space-time coarsening (see @ref exampleone) is used and \f$ O(\log_m N)
  \f$ for time-only coarsening.  The time-only coarsening option is the default
  and requires no user-written spatial interpolation/restriction routines
  (which is the case for space-time coasrening).  We note that the base of the
  logarithm is \f$ m \f$, which can be quite large.

## Cycling and relaxation strategies {#cyclingrelaxation}

There are two main cycling strategies available in XBraid, F-and V-cycles.  These two
cycles differ in how often and the order in which coarse levels are visited.  A V-cycle
is depicted next, and is a simple recursive application of the @ref twogrid.
  \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.1\textwidth]{../img/VCycle.pdf}}
   \end{figure}
   \endlatexonly

An F-cycle visits coarse grids more frequently and in a different order.
Essentially, an F-cycle uses a V-cycle as the post-smoother, which is an
expensive choice for relaxation.  But, this extra work gives you a closer
approximation to a two-grid cycle, and a faster convergence rate at the extra
expense of more work.  The effectiveness of a V-cycle as a relaxation scheme
can be seen in Figure 
\latexonly \ref{img:mgrit_cycles}, \endlatexonly
where one V-cycle globally propagates and *smoothes* the error.
The cycling strategy of an F-cycle is depicted next.
  \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.27\textwidth]{../img/FCycle.pdf}}
   \end{figure}
   \endlatexonly

Next, we make a few points about F- versus V-cycles.
- One V-cycle iteration is cheaper than one F-cycle iteration.
- But, F-cycles often converge more quickly.  For some test cases, this difference can be 
  quite large.  The cycle choice for the best time to solution will be problem
  dependent.  See @ref twodheat_scaling for a case study of cycling strategies. 
- For exceptionally strong F-cycles, the option [braid_SetNFMGVcyc](@ref braid_SetNFMGVcyc)
  can be set to use multiple V-cycles as relaxation.  This has proven useful
  for some problems with a strongly advective nature.

The number of FC relaxation sweeps is another important algorithmic setting.
Note that at least one F-relaxation sweep is always 
done on a level. A few summary points about relaxation are as follows.
- Using FCF, FCFCF, or FCFCFCF relaxation corresponds to passing
  *braid_SetNRelax* a value of 1, 2 or 3 respectively, and will result in an
  XBraid cycle that converges more quickly as the number of relaxations grows.
- But as the number of relaxations grows, each XBraid cycle becomes more
  expensive.  The optimal relaxation strategy for the best time to solution will
  be problem dependent.
- However, a good first step is to try FCF on all levels (i.e., *braid_SetNRelax(core, -1, 1)* ).
- A common optimization is to first set FCF on all levels (i.e., *braid_setnrelax(core, -1, 1)* ), 
  but then overwrite the FCF option on level 0 so that only F-relaxation is done on level 0, 
  (i.e., *braid_setnrelax(core, 0, 1)* ).  Another strategy is to use F-relaxation on all levels
  together with F-cycles. 
- See @ref twodheat_scaling  for a case study of relaxation strategies. 

Last, [Parallel Time Integration with Multigrid](https://computation.llnl.gov/project/linear_solvers/pubs/mgritPaper-2014.pdf)
has a more in depth case study of cycling and relaxation strategies

## Overlapping communication and computation {#overlapping}

  XBraid effectively overlaps communication and computation.
  The main computational kernel of XBraid is one relaxation sweep touching all the
  CF intervals. At the start of a relaxation sweep, each process first posts a
  non-blocking receive at its left-most point.  It then carries out
  F-relaxation in each interval, starting with the right-most interval to send
  the data to the neighboring process as soon as possible.  If each process
  has multiple CF intervals at this XBraid level, the strategy allows for
  complete overlap.
  \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.35\textwidth]{../img/overlap.pdf}}
       \label{img:overlap}
   \end{figure}
   \endlatexonly


## Configuring the XBraid Hierarchy {#config}

Some of the more basic XBraid function calls allow you to control aspects 
discussed here.
- [braid_SetFMG](@ref braid_SetFMG): switches between using F- and V-cycles.
- [braid_SetMaxIter](@ref braid_SetMaxIter ):  sets the maximum number of XBraid iterations
- [braid_SetCFactor](@ref braid_SetCFactor): sets the coarsening factor for any (or all levels)
- [braid_SetNRelax](@ref braid_SetNRelax): sets the number of CF-relaxation sweeps for any (or all levels)
- [braid_SetRelTol](@ref braid_SetRelTol), [braid_SetAbsTol](@ref braid_SetAbsTol): sets the stopping tolerance
- [braid_SetMinCoarse](@ref braid_SetMinCoarse): sets the minimum possible coarse grid size
- [braid_SetMaxLevels](@ref braid_SetMaxLevels): sets the maximum number of levels in the XBraid hierarchy

## Halting tolerance {#halting}

Another important configuration aspect regards setting a residual halting tolerance.
Setting a tolerance involves these three XBraid options:
1. [braid_PtFcnSpatialNorm](@ref braid_PtFcnSpatialNorm)
 
    This user-defined function carries out a spatial norm by taking the norm of a
    braid_Vector.  A common choice is the standard Eucliden norm (2-norm), but
    many other choices are possible, such as an L2-norm based on a finite
    element space.  

2. [braid_SetTemporalNorm](@ref braid_SetTemporalNorm)
   
   This option determines how to obtain a global space-time residual norm.
   That is, this decides how to combine the spatial norms returned by
   [braid_PtFcnSpatialNorm](@ref braid_PtFcnSpatialNorm) at each time step to
   obtain a global norm over space and time.  It is this global norm that
   then controls halting. 

   There are three *tnorm* options supported by
   [braid_SetTemporalNorm](@ref braid_SetTemporalNorm). We let the summation
   index *i* be over all C-point values on the fine time grid, *k* refer to the
   current XBraid iteration, *r* be residual values, *space_time* norms be a
   norm over the entire space-time domain and *spatial_norm* be the
   user-defined spatial norm from [braid_PtFcnSpatialNorm](@ref braid_PtFcnSpatialNorm).
   Thus, \f$ r_i \f$ is the residual at the *ith* C-point, and \f$ r^{(k)} \f$ 
   is the residual at the *kth* XBraid iteration.  The three options are then
   defined as,
  
   - *tnorm=1*: One-norm summation of spatial norms 
        \f[ \| r^{(k)} \|_{\mbox{space\textunderscore time}} = \Sigma_i \| r^{(k)}_i \|_{\mbox{spatial\textunderscore norm}} \f]
       If [braid_PtFcnSpatialNorm](@ref braid_PtFcnSpatialNorm) is the one-norm
       over space, then this is equivalent to the one-norm of the global
       space-time residual vector.

   - *tnorm=2*: Two-norm summation of spatial norms 
         \f[ \| r^{(k)} \|_{\mbox{space\textunderscore time}} = \left( \Sigma_i \| r^{(k)}_i \|^2_{\mbox{spatial\textunderscore norm}} \right)^{1/2} \f]
       If [braid_PtFcnSpatialNorm](@ref braid_PtFcnSpatialNorm) is the
       Euclidean norm (two-norm) over space, then this is equivalent to the
       Euclidean-norm of the global space-time residual vector.
  
   - *tnorm=3*: Infinity-norm combination of spatial norms
         \f[ \| r^{(k)} \|_{\mbox{space\textunderscore time}} = \max_i \| r^{(k)}_i \|_{\mbox{spatial\textunderscore norm}} \f]
       If [braid_PtFcnSpatialNorm](@ref braid_PtFcnSpatialNorm) is the
       infinity-norm over space, then this is equivalent to the infinity-norm
       of the global space-time residual vector.
  
   **The default choice is _tnorm=2_**
  


3. [braid_SetAbsTol](@ref braid_SetAbsTol), [braid_SetRelTol](@ref braid_SetRelTol)
   
   - If an absolute tolerance is used, then 
     \f[ \| r^{(k)} \|_{\mbox{space\textunderscore time}} < \mbox{tol} \f] 
     defines when to halt.

   - If a relative tolerance is used, then 
     \f[ \frac{ \| r^{(k)} \|_{\mbox{space\textunderscore time}} } {\| r^{(0)} \|_{\mbox{space\textunderscore time}} } < \mbox{tol} \f]
     defines when to halt.  That is, the current *kth* residual is scaled by
     the initial residual before comparison to the halting tolerance.  This is 
     similar to typical relative residual halting tolerances used in spatial multigrid,
     but can be a dangerous choice in this setting.

Care should be practiced when choosing a halting tolerance.  For instance, if a relative 
tolerance is used, then issues can arise when the initial guess is zero for large numbers
of time steps.  Taking the case where the initial guess (defined by 
[braid_PtFcnInit](@ref braid_PtFcnInit)) is 0 for all time values *t > 0*,
the initial residual norm will essentially only be nonzero at the first time value,

   \f[ \| r^{(0)} \|_{\mbox{space\textunderscore time}} \approx \| r^{(k)}_1 \|_{\mbox{spatial\textunderscore norm}} \f]

This will skew the relative halting tolerance, especially if the number of time steps 
increases, but the initial residual norm does not.  

A better strategy is to choose an absolute tolerance that takes your space-time domain
size into account, as in Section @ref twodheat_scaling, or to use an infinity-norm temporal
norm option.


## Debugging XBraid {#debugging}

Wrapping and debugging a code with XBraid typically follows a few steps. 

- Test your wrapped functions with XBraid test functions, e.g.,
  [braid_TestClone](@ref braid_TestClone) or 
  [braid_TestSum](@ref braid_TestSum). 
- Set max levels to 1 ([braid_SetMaxLevels](@ref braid_SetMaxLevels)) and run
  an XBraid simulation.  You should get the exact same answer as that achieved
  with sequential time stepping.  If you make sure that the time-grids used by
  XBraid and by sequential time stepping are bit-wise the same (by using the
  user-defined time grid option [braid_SetTimeGrid](@ref braid_SetTimeGrid) ),
  then the agreement of their solutions should be bit-wise the same.
- Continue with max levels equal to 1, but switch to two processors in time.
  Check that the answer again exactly matches sequential time stepping.  This
  test checks that the information in braid_Vector is sufficient to correctly
  start the simulation on the second processor in time.
- Set max levels to 2, halting tolerance to 0.0 
  ([braid_SetAbsTol](@ref braid_SetAbsTol)), max iterations to 3 
  ([braid_SetMaxIter](@ref braid_SetMaxIter )) and turn on the option 
  [braid_SetSeqSoln](@ref braid_SetSeqSoln).  
  This will use the solution from sequential time-stepping
  as the initial guess for XBraid and then run 3 iterations.  The residual
  should be exactly 0 each iteration, verifying the fixed-point nature of
  XBraid and a (hopefully!) correct implementation.  The residual may be on
  the order of machine epsilon (or smaller).  Repeat this test for multiple
  processors in time (and space if possible).
- A similar test turns on debug level printing by passing a print level of 3 
  to [braid_SetPrintLevel](@ref braid_SetPrintLevel).  This will print out
  the residual norm at each C-point.  XBraid with FCF-relaxation has the property
  that the exact solution is propagated forward two C-points each iteration.
  Thus, this should be reflected by numerically zero residual values for the first
  so many time points. Repeat this test for multiple
  processors in time (and space if possible).
- Finally, run some multilevel tests, making sure that the XBraid results are
  within the halting tolerance of the solutions generated by sequential
  time-stepping.  Repeat this test for multiple
  processors in time (and space if possible).
- Congratulations!  Your code is now verified.

# Computing Derivatives with XBraid_Adjoint {#xbraid_adjoint}

In many application scenarios, the ODE system is driven by some independent design parameters \f$\rho\f$. Those can be any time-dependent or time-independent parameters that uniquely determine the solution of the ODE (e.g. a boundary condition, material coefficients, etc.). In a discretized ODE setting, the user's time-stepping routine might then be written as
\f[
   u_i = \Phi_i(u_{i-1}, \rho), \quad \forall i=1, \dots N
\f]
where the time-stepper \f$\Phi_i\f$, which propagates a state at a time \f$t_{i-1}\f$ to the next time step at \f$t_i\f$, now also depends on the design parameters \f$\rho\f$. In order to quantify the simulation output for the given design, a real-valued objective function can then be set up that measures the quality of the ODE solution:   
\f[
   J(\mathbf{u}, \rho) \in \mathrm{R}.
\f] 
Here, \f$\mathbf{u} =  (u_0, \dots, u_N)\f$ denotes the space-time state solution for a given design. 

XBraid_Adjoint is a consistent discrete time-parallel adjoint solver for XBraid which provides sensitivity information of the output quantity \f$J\f$ with respect to the user-defined design parameters \f$\rho\f$.
The ability to compute sensitivities can greatly improve and enhance the simulation tool, as for example for solving design optimization or optimal control problems, parameter estimation for validation and verification purposes, error estimation or uncertainty quantification techniques. XBraid_Adjoint is non-intrusive with respect to the adjoint time-stepping scheme so that existing time-serial adjoint codes can be integrated easily though an extended user-interface. 

XBraid_Adjoint has been developed in collaboration with the Scientific Computing group at TU Kaiserslautern, Germany, and in particular with Dr. Stefanie Guenther and Prof. Nicolas Gauger.



## Short Introduction to Adjoint-based Sensitivity Computation {#adjointoverview}

Adjoint-based sensitivities compute the total derivative of \f$J\f$ with respect to changes in the design parameters \f$\rho\f$ by solving additional so-called adjoint equations. We will briefly introduce the idea in the following. You can skip this section, if you are familiar with adjoint sensitivity computation in general and move to @ref xbraid_adjointalgorithm immedately. Information on the adjoint method  can e.g. be found in [Giles, Pierce, 2000]
  \latexonly
   \footnote{ 
     Giles, M.B., Pierce, N.A.: ''An introduction to the adjoint approach to design.'' Flow, Turbulence and Combustion 65(3), 393–415 (2000)
   }
  \endlatexonly
  amongst many others. 

Consider an augmented (so-called *Lagrange*) funtion  
\f[
   L(\mathbf{u}, \rho) = J(\mathbf{u}, \rho) + \mathbf{\bar u}^TA(\mathbf{u}, \rho)
\f]
where the discretized time-stepping ODE equations in 
\f[
   A(\mathbf{u},\rho) := \begin{pmatrix} \Phi_1(u_0,\rho) - u_1 \\ 
                                       \vdots \\
                                       \Phi_N(u_{N-1}, \rho) - u_{N}
                        \end{pmatrix}
\f]
have been added to the objective function, multiplied with so-called *adjoint* variables \f$\mathbf{\bar u} = ({\bar u_1 , \dots, \bar u_N}) \f$. Since the added term is zero for all design and state variables that satisfy the discrete ODE equations, the total derivative of \f$J\f$ and \f$L\f$ with respect to the design match. Using the chain rule of differentiation, this derivative can be expresses as
\f[
    \frac{\mathrm{d}J}{\mathrm{d}\rho} = \frac{\mathrm{d}L}{\mathrm{d}\rho} = \frac{\partial J}{\partial {\mathbf u} }\textcolor{red}{\frac{\mathrm{d}\mathbf u}{\mathrm{d}\rho}} + \frac{\partial J}{\partial \rho} + \bar{\mathbf u}^T \left( \frac{\partial A}{\partial \mathbf u}  \textcolor{red}{\frac{\mathrm{d}\mathbf u}{\mathrm{d}\rho}} +  \frac{ \partial A}{\partial \rho} \right)
\f]
where \f$\partial\f$ denotes partial derivatives -- in contrast to the total derivative (i.e. the sensitivity) denoted by \f$\mathrm{d}\f$.

When computing this derivative, the terms in red are the ones that are computationally most expensive. In fact, the cost for computing these sensitivities scale linearly with the number of design parameters, i.e. the dimension of \f$\rho\f$ (consider e.g. a Finite Differences setting, where a recomputation of the state would be necessary for each perturbations of the design into all unit direction of the design space). In order to avoid these costs, the adjoint method aims at setting the adjoint variable \f$\mathbf{\bar u}\f$ such that these red terms add up to zero in the above expression. Hence, when solving 
\f[
  \left(\frac{ \partial J}{\partial \mathbf u}\right)^T + \left(\frac{\partial A}{\partial \mathbf u}\right)^T \bar {\mathbf u} = 0  
\f]
for the adjoint variable \f$\mathbf{\bar u}\f$ first, then the so-called *reduced gradient* of \f$J\f$, which is the transpose of the total derivative of \f$J\f$ with respect to the design, is given by 
\f[
   \left(\frac{\mathrm{d}J}{\mathrm{d}\rho}\right)^T = \left(\frac{\partial J} {\partial \rho}\right)^T  + \left(\frac{\partial A}{\partial \rho}\right)^T\bar{\mathbf u}
\f]
The advantage of this strategy is, that in order to compute the sensitivity of \f$J\f$ with respect to \f$\rho\f$, only one additional (adjoint) equation for \f$\mathbf{\bar u}\f$ has to be solved, plus evaluating the partial derivatives. The computational cost for computing \f$\mathrm{d}J/\mathrm{d}\rho\f$ therefore does not scale anymore with the number of design parameters. 

For the time-dependent discrete ODE problem, the adjoint equation from above reads
\f[
 \text{\textcolor{red}{unsteady adjoint:}}\qquad \quad \bar u_i = \partial_{u_i} J(\mathbf{u}, \rho)^T + \left(\partial_{u_i}\Phi_{i+1}(u_i, \rho)\right)^T\bar{u}_{i+1} \qquad \forall i = N\dots, 1
\f]
using the terminal condition \f$u_{N+1} := 0 \f$. The reduced gradient is given by
\f[
   \text{\textcolor{red}{reduced gradient:}} \qquad \qquad\qquad \left(\frac{\partial J}{\partial \rho}\right)^T = \partial_{\rho} J(\mathbf{u}, \rho)^T + \sum_{i=1}^N \left(\partial_{\rho}\Phi_{i}(u_{i-1}, \rho)\right)^T\bar u_{i} 
\f] 


## Overview of the XBraid_Adjoint Algorithm {#xbraid_adjointalgorithm}

The unsteady adjoint equations can in principle be solved ``backwards in time'' in a time-serial manner, starting from the terminal condition \f$\bar u_{N+1} = 0\f$. However, the parallel-in-time XBraid_Adjoint solver offers speedup by distributing the backwards-in-time phase
onto multiple processors along the time domain. Its implementation is based on techniques of the reverse-mode of Automatic Differentiation applied to one primal XBraid iteration. To that end, each primal iteration is augmented by an objective function evaluation, followed by updates for the space-time adjoint variable \f$\mathbf{\bar u}\f$, as well as evaluation of the reduced gradient denoted by \f$\bar\rho\f$. In particular, the following so-called *piggy-back* iteration is performed:

1. **XBraid**: update the state and evaluate the objective function
   \f[
      \mathbf{u}^{(k+1)} \leftarrow \text{XBraid}(\mathbf{u}^{(k)}, \rho), \quad
      J \leftarrow J(\mathbf{u}^{(k)}, \rho) 
   \f]
2. **XBraid_Adjoint**: update the adjoint and evaluate the reduced gradient
   \f[
      \mathbf{\bar u}^{(k+1)} \leftarrow \text{XBraid\_Adjoint}(\mathbf{u}^{(k)}, \mathbf{\bar u}^{(k)} \rho),  \quad
      \bar \rho \leftarrow \left(\frac{\mathrm{d} J(\mathbf u^{(k)}, \rho)}{\mathrm{d} \rho}\right)^T
   \f]

Each XBraid_Adjoint moves backwards though the primal XBraid multigrid cycle. It collects local partial derivatives of the elemental XBraid operations in reverse order and concatenates them using the chain rule of differentiation - this is the basic idea of the reverse mode of Automatic Differentiation. 
This yields a consistent discrete time-parallel adjoint solver that inherits parallel scaling properties from the primal XBraid solver. 

Further, XBraid_Adjoint is non-intrusive to existing adjoint sequential time marching schemes. It adds additional user routines to the primal XBraid interface which provide the propagation of sensitivities of the forward time stepper backwards in time and the evaluation of partial derivatives of the local objective function at each time step. 
In cases where a time-serial unsteady adjoint solver is already available, this backwards time stepping capability can be easily wrapped into the adjoint user interface with little extra coding.

The adjoint update in the above piggy-back iteration converge at the same convergence rate as the primal state variables. However since the adjoint equations depend on the state solution, the adjoint convergence will slightly lag behind the convergence of the state.
More information on convergence results and implementational details for XBraid_Adjoint can be found in [Gunther, Gauger, Schroder, 2017].
  \latexonly
   \footnote{ 
   G{\"u}nther, S., Gauger, N.R. and Schroder, J.B. ''A Non-Intrusive Parallel-in-Time Adjoint Solver with the XBraid Library.'' Computing and Visualization in Science, Springer, (accepted), (2017)
   }
  \endlatexonly


## Overview of the XBraid_Adjoint Code {#xbraid_adjointcode}

XBraid_Adjoint offers a non-intrusive approach for time-parallelization of existing time-serial adjoint codes. To that end, an extended user-interface allows the user to wrap their existing code for evaluating the objective function and performing a backwards-in-time adjoint step into routines that are passed to XBraid_Adjoint.

### Objective function evaluation {#xbraid_adjoint_objective}


The user-interface for XBraid_Adjoint allows for objective functions of the following type:
\f[
   J =  F\left( \int_{t_0}^{t^1} f(u(t),\rho) ~ \mathrm{d}t \right)
\f]
It involves a time-integral part of some time-dependent quantity of interest \f$f\f$ as well as a *postprocessing* function \f$F\f$. 
The time-interval boundaries \f$t_0,t_1\f$ can be set using the options [braid_SetTStartObjective](@ref braid_SetTStartObjective) and [braid_SetTStopObjective](@ref braid_SetTStopObjective), otherwise the entire time domain will be considered. Note that these options can be used for objective functions that are only evaluated at one specific time instance  by setting \f$t_0 = t_1\f$ (e.g. in cases where only the last time step is of interest). 
The postprocessing function \f$F\f$ offers the possibility to further modify the time-integral, e.g. for setting up a tracking-type objective function (substract a target value and square), or adding relaxation or penalty terms.
While \f$f\f$ is mandatory for XBraid_Adjoint, the postprocessing routine \f$F\f$ is optional and is passed to XBraid_Adjoint though the optional [braid_SetPostprocessObjective](@ref braid_SetPostprocessObjective) routine. 
XBraid_Adjoint will perform the time-integration by summing up the \f$f\f$ evaluations in the given time-domain
\f[
   I \leftarrow \sum_{i = i_0}^{i_1} f(u_i, \rho) 
\f]
followed by a call to the postprocessing function \f$F\f$, if set:
\f[
   J \leftarrow F\left( I, \rho \right).
\f]


### Partial derivatives of user-routines

The user needs to provide the derivatives of the time-stepper \f$\Phi\f$ and function evaluation \f$f\f$ (and potentially \f$F\f$) for XBraid_Adjoint.
Those are provided in terms of transposed matrix-vector products in the following way:

1. **Derivatives of objective function** \f$J\f$:

   - **Postprocessing** \f$F\f$: If the postprocessing routine has been set, the user needs to provide it's transposed partial derivatives in the following way:
   \f[
      \bar F \leftarrow \frac{\partial F(I, \rho)}{\partial I}
   \f] 
   \f[
      \bar \rho \leftarrow \rho + \frac{\partial F(I,\rho)}{\partial \rho} 
   \f]
   
   - **Time-dependent part** \f$f\f$: The user provides a routine that evaluates the following transposed partial derivatives of \f$f\f$ multiplied with the scalar input \f$\bar F\f$:
   \f[  
      \bar u_i \leftarrow \left(\frac{\partial f(u_i, \rho)}{\partial u_i}\right)^T \bar F 
   \f]
   \f[
      \bar \rho \leftarrow \bar \rho + \left(\frac{\partial f(u_i, \rho)}{\partial \rho}\right)^T \bar F
   \f]
   The scalar input \f$\bar F\f$ equals \f$1.0\f$, if no postpocessing function \f$F\f$ has been set. 


2. **Derivatives of the time-stepper** \f$\Phi_i\f$:
   The user provides a routine that computes the following transposed partial derivatives of \f$\Phi_i\f$ multiplied with the adjoint input vector \f$\bar u_i\f$:
   \f[
      \bar u_{i} \leftarrow \left(\frac{\partial \Phi(u_i, \rho)}{\partial u_i}\right)^T\bar u_i
   \f]
   \f[
      \bar \rho \leftarrow \bar \rho +  \left(\frac{\partial \Phi(u_i, \rho)}{\partial \rho}\right)^T\bar u_i
   \f]

Note that the partial derivatives with respect to \f$\rho\f$ always *update* the reduced gradient  \f$\bar \rho\f$ instead of overwriting it (i.e. \f$\mathop{+}=\f$). Therefore, the gradient needs to be reset to zero before each iteration of XBraid_Adjoint, which is taken care of by XBraid_Adjoint calling an additional user-defined routine [braid_PtFcnResetGradient](@ref braid_PtFcnResetGradient). 


Depending on the nature of the design variables, it is neccessary to gather gradient information in \f$\bar \rho\f$ from all time-processors after XBraid_Adjoint has finished. It is the user's responsibility to do that, if needed, e.g. through a call to MPI_Allreduce.   


### Halting tolerance
Similar to the primal XBraid algorithm, the user can choose a halting tolerance for XBraid_Adjoint which is based on the adjoint residual norm. An absolute tolerance ([braid_SetAbsTolAdjoint](@ref braid_SetAbsTolAdjoint)) 
   \f[
      \|\bar{\mathbf{u}}^{(k)} - \bar{\mathbf{u}}^{(k-1)} \|_{\text{space\_time}} < \text{tol\_adjoint}
   \f]
   or a relative tolerance ([braid_SetRelTolAdjoint](@ref braid_SetRelTolAdjoint))
   \f[
      \frac{\|\bar{\mathbf{u}}^{(k)} - \bar{\mathbf{u}}^{(k-1)} \|_{\text{space\_time}}}{\|\bar{\mathbf{u}}^{(1)} - \bar{\mathbf{u}}^{(0)} \|_{\text{space\_time}}} < \text{tol\_adjoint}
   \f]
   can be chosen.

### Finite Difference Testing 
You can verify the gradient computed from XBraid\_Adjoint using Finite Differences. Let \f$e_i\f$ denote the \f$i\f$-th unit vector in the design space, then the i-th entry of the gradient should match with 
\f[
   i\text{-th Finite Difference: } \quad \frac{J(\bold u_{\rho + he_i}, \rho + he_i) - J(\bold u, \rho)}{h} 
\f]
for a small perturbation \f$h>0\f$. Here, \f$\bold u_{\rho + he_i}\f$ denotes the new state solution for the perturbed design variable. 
Keep in mind, that round-up errors have to be considered when computing the Finite Differences for very small perturbations \f$h \to 0\f$. Hence, you should vary the parameter to find the best fit. 

In order to save some computational work while computing the perturbed objective function value, XBraid_Adjoint can run in `ObjectiveOnly` mode, see [braid_SetObjectiveOnly](@ref braid_SetObjectiveOnly). When in this mode, XBraid_Adjoint will only solve the ODE system and evaluate the objective function, without actually computing its derivative. This option might also be useful within an optimization framework e.g. for implementing a line-search procedure.  



### Getting  started

* Look at the simple example in `examples/ex-01-adjoint.c` which implements XBraid_Adjoint sensitivity computation for a scalar ODE.

   

# Citing XBraid {#braidcite}

To cite XBraid, please state in your text the version number from the 
VERSION file, and please cite the project website in your bibliography as
   
>  [1] XBraid: Parallel multigrid in time. http://llnl.gov/casc/xbraid.

The corresponding BibTex entry is
   
     @misc{xbraid-package,
       title = {{XB}raid: Parallel multigrid in time},
       howpublished = {\url{http://llnl.gov/casc/xbraid}}
       }

# Summary {#summary}

- XBraid applies multigrid to the time dimension.
 + This exposes concurrency in the time dimension.
 + The potential for speedup is large, 10x, 100x, ...
- This is a non-intrusive approach, with an unchanged time discretization defined by user.
- Parallel time integration is only useful beyond some scale.  
  This is evidenced by the experimental results below.  For smaller numbers
  of cores sequential time stepping is faster, but at larger core counts
  XBraid is much faster.
- The more time steps that you can parallelize over, the better your speedup will be.
- XBraid is optimal for a variety of parabolic problems (see the examples directory).
- XBraid_Adjoint provides time-parallel adjoint-based sensitivities of output quantities with respect to user-defined design variables
 + It is non-intrusive with respect to existing adjoint time-marching schemes
 + It inherits parallel scaling properties from XBraid 

