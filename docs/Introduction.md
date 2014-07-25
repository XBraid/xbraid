# Overview of Warp Algorithm

The goal of Warp is to solve a problem faster than is possible with a
traditional time marching algorithm.  Instead of sequential time marching, Warp
solves the problem iteratively by simultaneously updating the current solution
guess over all time values.  The initial solution guess can be anything, even a
random function over space-time.  The iterative updates to the solution guess
are done by constructing a hierarchy of temporal grids, where the finest grid
contains all of the time values for the simulation.  Each subsequent grid is a
coarser grid with fewer time values.  The coarsest grid has a trivial number of
time steps and can be quickly solved exactly.  The overall effect is that
solutions to the time marching problem on the coarser (i.e., cheaper) grids can
be used to correct the original finest grid solution.

To understand how Warp differs from traditional time marching, consider the
simple linear advection equation, \f$ u_t = -c u_x \f$.  The next figure depicts how one
would typically evolve a solution here with sequential time stepping.  The solution propagates
sequentially across space as time increases.
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

Warp instead begins with a solution guess over all of space-time, which we let here be random for the 
sake of argument.  A Warp iterations then does

1. Relaxation on the fine grid, i.e., the grid that contains all of the desired time values
   - Relaxation is just a local application of the time stepping scheme, e.g., backward Euler
2. Restriction to the first coarse grid, i.e., a grid that contains fewer time values, 
  say every second or every third time value
3. Relaxation on the first coarse grid
4. Restriction to the second coarse grid and so on...
5. When a coarse grid of trivial size (say 2 time steps) is reached, it is solved exactly.
6. The solution is then interpolated from the coarsest grid to the finest grid

After the cycle is complete, it repeats until the solution is accurate enough. 
This is depicted in the next figure, where only a few iterations are required for this simple problem.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=1.0\textwidth]{../img/3_MG_In_Time_Iterations.pdf}}
       \caption{Warp iterations.}
       \label{img:mgrit_cycles}
   \end{figure}
   \endlatexonly

There are a few important points to make.
- The coarse time grids allow for global propagation of information across
  space-time with only one Warp iteration (c.f. how the solution is updated
  from iteration 0 to iteration 1).
- Using coarser (cheaper) grids to correct the fine grid is analagous to spatial multigrid.
- Only a few Warp iterations are required to find the solution over 1024 time steps. 
  Therefore if enough processors are available to parallelize Warp, we can see a speedup
  over traditional time stepping (more on this later).
- This is a simple example, and Warp is structured to handle variable time
  step sizes and adaptive time step sizes, and these features will be coming. 

To firm up our understanding, let`s do a little math.  Assume that you
have a general ODE, 
   \f[
   \Xux^{\prime}(t) = \Xfx(t, \Xux(t)), ~~~ \Xux(0) = \Xux_0, ~~~ t \in [0,T], 
   \f]
which you descretize with the one-step integration
   \f[
   \Xux_i = \Phi_i(\Xux_{i-1}) + \Xgx_i , ~~~ i = 1, 2, ... , N.
   \f]
Traditional time marchine would first solve for \f$ i = 1 \f$, then solve for 
\f$ i = 2 \f$, and so on.  This process is equivalent to a forward solve of this
system,
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
This process is optimal and O(N), but it is sequential.  Warp instead solves
the system iteratively, with a multigrid reduction method 
   \latexonly
   \footnote{ Ries, Manfred, Ulrich Trottenberg, and Gerd Winter. "A note on MGR methods." Linear Algebra and its Applications 49 (1983): 1-26.}
   \endlatexonly
applied in only the time dimension. This approach is
- nonintrusive, in that it coarsens only in time and asks the user to define their own 
  \f$ \Phi \f$
- optimal and O(N) with a higher constant than time stepping
- highly parallel

Warp solves this system iteratively by constructing a hierarchy of time grids.
We describe the two-grid process, with the multigrid process being a recursive
application of the process.  We also assume that \f$ \Phi \f$ is constant
for notational simplicity.  

Warp functions as follows.
The next Figure depicts a sample timeline of time values, where the time values
have been split into C and F points.  C points exist on both the fine and
coarse time grid, but F points exist only on the fine time scale.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.75\textwidth]{../img/timeline.pdf}}
       \label{img:timeline}
   \end{figure}
   \endlatexonly
Following the above steps outlining a Warp cycle, the first task is relaxation.
The simplest relaxation here alternates between C and F sweeps. An F sweep simply
updates time values by integrating with \f$ \Phi \f$ over all the F points from one C point 
to the next, as depicted here
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.15\textwidth]{../img/Frelaxation.pdf}}
       \label{img:Frelax}
   \end{figure}
   \endlatexonly

But, such an update can be done simultaneously over all F intervals in parallel, as 
depicted here
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.75\textwidth]{../img/FrelaxationDetail.pdf}}
       \label{img:FrelaxDetail}
   \end{figure}
   \endlatexonly

Following an F sweep we can also do C sweep, as depicted here
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.75\textwidth]{../img/CrelaxationDetail.pdf}}
       \label{img:CrelaxDetail}
   \end{figure}
   \endlatexonly
FCF or F relaxation will refer to how the relaxation sweeps are carried out. So, we can say
- FCF or F relaxtion is highly parallel.
- But, a sequential component exists equaling the the number of F points between two C points.
- Warp uses regular coarsening factors, i.e., the spacing of C points happens every k points.

After relaxation, then comes coarse grid correction.  The restriction operator
\f$ R \f$ maps to the coarse grid by simply injecting values at C points from
the fine grid to the coarse grid, 
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
   \end{pmatrix},
   \f]
where the spacing between each \f$ I \f$ is \f$ m-1 \f$ block rows.  This 
Warp implements an FAS (Full Approximation Scheme) multigrid cycle, and hence
the system matrix , solution guess and right-hand-side (i.e., \f$ A, \mathbf{u}, \mathbf{g}\f$)
are all restricted.  This is in contrast to linear multigrid which typically 
restricts the residual equation to the coares grid.  We choose FAS because it is
``nonlinear`` multigrid and allows us to solve nonlinear problems.

The main question here is how to form the coarse grid matrix, which in turn is 
defined by the coarse grid time stepper \f$ \Phi_{\Delta} \f$.  It is typical
to let \f$ \Phi_{\Delta} \f$ simply be \f$ \Phi\f$ but with the coarse time step
size \f$ \Delta t \f$.  Thus if 
   \f[
   A =
   \begin{pmatrix}
    I & & & \\
   -\Phi & I & & \\
    & \ddots & \ddots & \\
    & & -\Phi & I
   \end{pmatrix}
   \f]
then
   \f[
   A_{\Delta} =
   \begin{pmatrix}
    I & & & \\
   -\Phi_{\Delta} & I & & \\
    & \ddots & \ddots & \\
    & & -\Phi_{\Delta} & I
   \end{pmatrix},
   \f]
where \f$ A_{\Delta} \f$ has fewer rows and columns than \f$ A \f$, e.g., if we are coarsening
in time by 2, it has one half as many rows and columns.  This coarse grid equation
   \f[
   A_{\Delta} \mathbf{u}_{\Delta} = \mathbf{g}_{\Delta} \equiv R \mathbf{g}
   \f]
is then solved.  Finally,
\f$ \mathbf{u}_{\Delta} \f$ is interpolated with \f$ P_{\Phi} \f$ back to the
fine grid, which is equivalent to injecting the coarse grid to the C points on
the fine grid, follwed by an F relaxation sweep.  That is,  
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
where \f$ m \f$ again is the coarsening factor.  

This two-grid process is captured with this algorithm.  Using a recursive coarse
grid solves makes the process multilevel, and halting is done based on a 
residual tolerance.  If the operator is linear, this FAS cycle is equivalent
to standard linear multigrid.  Note that we represent \f$ A \f$ as function below,
whereas above the notation was simplified for the linear case.
   1. Relax on \f$A(\mathbf{u}) = \mathbf{g}\f$ using FCF-relaxation
   2. Restrict the fine grid approximation and its residual:
   \f$ \mathbf{u}_{\Delta} \leftarrow  R \mathbf{u}, \quad \mathbf{r}_{\Delta} \leftarrow  R (\mathbf{g} - A(\mathbf{u}) \f$ 
   3. Solve \f$A_{\Delta}(\mathbf{v}_{\Delta})  = A_{\Delta}(\mathbf{u}_{\Delta}) + \mathbf{r}_{\Delta} \f$
   4. Compute the coarse grid error approximation: \f$ \mathbf{e}_{\Delta} = \mathbf{v}_{\Delta} - \mathbf{u}_{\Delta} \f$
   5. Correct: \f$ \mathbf{u} \leftarrow \mathbf{u} + P \mathbf{e}_{\Delta} \f$

In summary, a few points are
- Warp is an iterative solver for the global space-time problem.
- The user defines the time stepping routine \f$ \Phi \f$.
- Warp convergence will depend heavily on how well \f$ \Phi_{\Delta} \f$ approximates
  \f$ \Phi^m \f$, that is how well a time step size of \f$ m \delta t = \Delta T \f$
  will approximate \f$ m \f$ applications of the same time integrator for a time step 
  size of \f$ \delta t \f$.  This is a subject of research, but this approximation
  need not capture fine scale behavior, which is captured by relaxation on the fine grid.
- The coarsest grid is solved exactly, i.e., sequentially, which is a bottleneck.
  This is an issue for two-level methods like Parareal,
  \latexonly
   \footnote{ Lions, J., Yvon Maday, and Gabriel Turinici. "A''parareal''in time discretization of PDE's." Comptes Rendus de l'Academie des Sciences Series I Mathematics 332.7 (2001): 661-668.}
  \endlatexonly
  but not for a multilevel scheme like Warp where the coarsest grid is of trivial size.
- By forming the coarse grid to have the same sparsity structure and time stepper
  as the fine grid, the algorithm can recur easily and efficiently.
- Interpolation and restriction are ideal or exact, in that an application of 
  interpolation leaves a zero residual at all F points.


# Overview of Warp Code

Warp is designed to run in conjunction with an existing application code that
can be wrapped per our interface.  This application code will implement some
time marching type simulation like fluid flow.  Essentially, the user has to
take their application code and extract a stand-alone time-stepping function
\f$ \Phi \f$ that can evolve a solution from one time value to another,
regardless of time step size.  After this is done, the Warp code takes care of
the parallelism in the time dimension.

Warp 
- is written in C/C++ and can easily interface with Fortran
- uses MPI for parallelism
- self documents through comments in the source code and through *.md files

## Parallel decomposition and memory
- Warp decomposes the problem in parallel as depicted in this figure.
  \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.75\textwidth]{../img/data_layout.pdf}}
       \label{img:data_layout}
   \end{figure}
   \endlatexonly
   As you can see, traditional time stepping only stores one time step as a time, 
   but only enjoys a spatial data decomposition and spatial parallelism.  On the other
   hand, Warp stores multiple time steps simultaneously and each processor holds a space-time
   chunk.

- Warp only handles temporal parallelism and is agnostic to the spatial decomposition.  
  See the [warp_SplitCommworld](@ref warp_SplitCommworld).  
  Each processor owns a certain number of CF intervals of points, as depicted here
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.75\textwidth]{../img/parallel_timeline.pdf}}
       \label{img:data_layout}
   \end{figure}
   \endlatexonly
   One interval is a C point and all following F points until the next C point.  These 
   intervals are distributed evenly on the finest grid.

- Storage is greatly minimized by only stored C points.  Whenever an F point is needed,
  it is generated by F relaxation.  That is, we only store the red C point time values.
  \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.75\textwidth]{../img/timeline.pdf}}
   \end{figure}
   \endlatexonly
   Coarsening can by aggressive with \f$ m = 8, 16, 32 \f$, so the storage requirements
   of Warp are significantly reduced when compared to storing all of the time values.

- In practice, storing only one space-time slab is advisable.  That is, solve
  for as many time steps (say k time steps) as you have available memory for.
  Then move on to the next k time steps.

## Overlapping communication and computation

  Warp effectively overlaps communication and computation.  The main computation kernel of Warp
  is an F relaxation sweep, and each processor first posts a send at it`s left-most C point, and then
  carries out F relaxation on its right-most interval.  If each processor has multiple intervals at this Warp
  level, this should allow for complete overlap.
  \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.35\textwidth]{../img/overlap.pdf}}
       \label{img:overlap}
   \end{figure}
   \endlatexonly


# Properties of Warp

- Warp applies multigrid to the time dimension
 + Exposes concurrency in the time dimension
 + Potential for speedup is large, 10x, 100x, ...
- Non-intrusive approach, with unchanged time discretization defined by user
- Parallel time integration is only useful beyond some scale.  
  This is evidenced by the experimental results below.  For smaller numbers
  of cores sequential time stepping is faster, but at larger core counts
  Warp is much faster.
- The more time steps that you can parallelize over, the better your speedup.
- Warp is optimal for a variety of parabolic problems (see the examples).

Here is some experimental data for the 2D heat equation generated by 
the drive-02 example.
  \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.45\textwidth]{../img/heat_results.pdf}}
       \label{img:heat_results}
   \end{figure}
   \endlatexonly
Here, the speedup is 10x on a Linux cluster using 4 cores per node, a
Sandybridge Intel chipset, and a fast Infiniband interconnect.  The problem
size was \f$ 129^2 \times 16,192 \f$.  The coarsening factor \f$ m = 16 \f$ 
on the finest level and 2 on coarser levels.  And various relaxation and V and F
cycling strategies are experimented with.  

