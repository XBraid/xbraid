<!--
 - Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 - Produced at the Lawrence Livermore National Laboratory.
 - This file is part of WARP.  See file COPYRIGHT for details.
 -
 - WARP is free software; you can redistribute it and/or modify it under the
 - terms of the GNU Lesser General Public License (as published by the Free
 - Software Foundation) version 2.1 dated February 1999.
 -->
 
# Introduction {#intro}

# Overview of Warp Algorithm {#warpoverview}

The goal of Warp is to solve a problem faster than a
traditional time marching algorithm.  Instead of sequential time marching, Warp
solves the problem iteratively by simultaneously updating a space-time solution
guess over all time values.  The initial solution guess can be anything, even a
random function over space-time.  The iterative updates to the solution guess
are done by constructing a hierarchy of temporal grids, where the finest grid
contains all of the time values for the simulation.  Each subsequent grid is a
coarser grid with fewer time values.  The coarsest grid has a trivial number of
time steps and can be quickly solved exactly.  The effect is that
solutions to the time marching problem on the coarser (i.e., cheaper) grids can
be used to correct the original finest grid solution.  Thus, a problem with 
many time steps (thousands, tens of thousands or more) can be solved with 10 or 15
Warp iterations, and the overall time to solution can be greatly sped up.  However, 
this is achieved at the cost of more computational resources.

To understand how Warp differs from traditional time marching, consider the
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

Warp instead begins with a solution guess over all of space-time, which for demonstration, 
we let be random. A Warp iteration then does

1. Relaxation on the fine grid, i.e., the grid that contains all of the desired time values
   - Relaxation is just a local application of the time stepping scheme, e.g., backward Euler
2. Restriction to the first coarse grid, i.e., interpolate the problem to a grid that 
   contains fewer time values, say every second or every third time value
3. Relaxation on the first coarse grid
4. Restriction to the second coarse grid and so on...
5. When a coarse grid of trivial size (say 2 time steps) is reached, it is solved exactly.
6. The solution is then interpolated from the coarsest grid to the finest grid

One Warp iteration is called a *cycle* and these cycles continue until the the
solution is accurate enough.  This is depicted in the next figure, where only a
few iterations are required for this simple problem.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=1.0\textwidth]{../img/3_MG_In_Time_Iterations.pdf}}
       \caption{Warp iterations.}
       \label{img:mgrit_cycles}
   \end{figure}
   \endlatexonly

There are a few important points to make.
- The coarse time grids allow for global propagation of information across
  space-time with only one Warp iteration.  This is visible in the above figure by observing
  how the solution is updated from iteration 0 to iteration 1.
- Using coarser (cheaper) grids to correct the fine grid is analagous to spatial multigrid.
- Only a few Warp iterations are required to find the solution over 1024 time steps. 
  Therefore if enough processors are available to parallelize Warp, we can see a speedup
  over traditional time stepping (more on this later).
- This is a simple example, with evenly space time steps.  Warp is structured 
   to handle variable time step sizes and adaptive time step sizes, and these 
   features will be coming. 

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
- nonintrusive, in that it coarsens only in time and the user defines \f$ \Phi \f$
  + Thus, users can continue using existing time stepping codes by wrapping them
    into our framework.
- optimal and O(N), but O(N) with a higher constant than time stepping
  + Thus with enough computational resources, Warp will outperform sequential time stepping.
- highly parallel

Warp solves this system iteratively by constructing a hierarchy of time grids.
We describe the two-grid process, with the multigrid process being a recursive
application of the process.  We also assume that \f$ \Phi \f$ is constant
for notational simplicity.  

Warp functions as follows.
The next figure depicts a sample timeline of time values, where the time values
have been split into C- and F-points.  C-points exist on both the fine and
coarse time grid, but F-points exist only on the fine time scale.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.75\textwidth]{../img/timeline.pdf}}
       \label{img:timeline}
   \end{figure}
   \endlatexonly
The first task is relaxation and an effective relaxation alternates between C and F 
sweeps (this is like line-relaxation in space in that the residual is set to 0 for 
an entire time step). An F sweep simply updates time values by integrating with 
\f$ \Phi \f$ over all the F-points from one C-point to the next, as depicted next.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.15\textwidth]{../img/Frelaxation.pdf}}
       \label{img:Frelax}
   \end{figure}
   \endlatexonly

But, such an update can be done simultaneously over all F intervals in parallel, as 
depicted next.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.75\textwidth]{../img/FrelaxationDetail.pdf}}
       \label{img:FrelaxDetail}
   \end{figure}
   \endlatexonly

Following an F sweep we can also do C sweep, as depicted next.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.75\textwidth]{../img/CrelaxationDetail.pdf}}
       \label{img:CrelaxDetail}
   \end{figure}
   \endlatexonly

\latexonly \newpage \endlatexonly
In general, FCF- and F-relaxation will refer to the relaxation methods used in Warp. We can say
- FCF or F-relaxtion is highly parallel.
- But, a sequential component exists equaling the the number of F-points between two C-points.
- Warp uses regular coarsening factors, i.e., the spacing of C-points happens every \f$k\f$ points.

After relaxation, comes coarse grid correction.  The restriction operator
\f$ R \f$ maps fine grid quantities to the coarse grid by simply injecting 
values at C-points from the fine grid to the coarse grid, 
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
where the spacing between each \f$ I \f$ is \f$ m-1 \f$ block rows.  
Warp implements an FAS (Full Approximation Scheme) multigrid cycle, and hence
the solution guess and residual 
(i.e., \f$ A, \mathbf{u}, \mathbf{g} - A \mathbf{u}\f$)
are restricted.  This is in contrast to linear multigrid which typically 
restricts only the residual equation to the coarse grid.  We choose FAS because it is
*nonlinear* multigrid and allows us to solve nonlinear problems.  FAS was invented
by Achi Brandt, but this 
[PDF](http://computation.llnl.gov/casc/people/henson/postscript/UCRL_JC_150259.pdf) 
by Van Henson is a good intro. 

The main question here is how to form the coarse grid matrix, which in turn asks how 
to define the coarse grid time stepper \f$ \Phi_{\Delta} \f$.  It is typical
to let \f$ \Phi_{\Delta} \f$ simply be \f$ \Phi\f$ but with the coarse time step
size \f$ \Delta T = m \delta t \f$.  Thus if 
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
in time by 2, \f$ A_{\Delta} \f$ has one half as many rows and columns.  This coarse grid equation
   \f[
   A_{\Delta} \mathbf{v}_{\Delta} = \mathbf{g}_{\Delta}
   \f]
is then solved, where the right-hand-side is defined by FAS (see @ref twogrid).  
Finally, FAS defines a coarse grid error approximation \f$ \mathbf{e}_{\Delta} \f$, which is interpolated 
with \f$ P_{\Phi} \f$ back to the fine grid and added to the current solution guess.  Interpolation 
is equivalent to injecting the coarse grid to the C-points on the fine grid, follwed by an F-relaxation 
sweep.  That is,  
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
where \f$ m \f$ is the coarsening factor.  

## Two-Grid Algorithm {#twogrid}

This two-grid process is captured with this algorithm.  Using a recursive coarse
grid solve (i.e., step 3 becomes a recursive call) makes the process multilevel. 
Halting is done based on a residual tolerance.  If the operator is linear, this FAS cycle is equivalent
to standard linear multigrid.  Note that we represent \f$ A \f$ as a function below,
whereas the above notation was simplified for the linear case.
   1. Relax on \f$A(\mathbf{u}) = \mathbf{g}\f$ using FCF-relaxation
   2. Restrict the fine grid approximation and its residual:
   \f$ \mathbf{u}_{\Delta} \leftarrow  R \mathbf{u}, \quad \mathbf{r}_{\Delta} \leftarrow  R (\mathbf{g} - A(\mathbf{u}) \f$ 
   3. Solve \f$A_{\Delta}(\mathbf{v}_{\Delta})  = A_{\Delta}(\mathbf{u}_{\Delta}) + \mathbf{r}_{\Delta} \f$
   4. Compute the coarse grid error approximation: \f$ \mathbf{e}_{\Delta} = \mathbf{v}_{\Delta} - \mathbf{u}_{\Delta} \f$
   5. Correct: \f$ \mathbf{u} \leftarrow \mathbf{u} + P \mathbf{e}_{\Delta} \f$

## Summary {#algorithmsummary}

In summary, a few points are
- Warp is an iterative solver for the global space-time problem.
- The user defines the time stepping routine \f$ \Phi \f$ and can wrap 
  existing code to accomplish this. 
- Warp convergence will depend heavily on how well \f$ \Phi_{\Delta} \f$ approximates
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
  but not for a multilevel scheme like Warp where the coarsest grid is of trivial size.
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
  levels because, there would be no more points to coarsen!
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.6\textwidth]{../img/3_levels.pdf}}
       \label{img:FrelaxDetail}
   \end{figure}
   \endlatexonly
   
   By default, Warp will subdivide the time domain into evenly sized time
   steps.  Warp is structured to handle variable time step sizes and adaptive
   time step sizes, and these features are coming. 

# Overview of Warp Code {#codeoverview}

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


## Parallel decomposition and memory {#decomposition}
- Warp decomposes the problem in parallel as depicted next. 
  \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.75\textwidth]{../img/data_layout.pdf}}
       \label{img:data_layout}
   \end{figure}
   \endlatexonly
   As you can see, traditional time stepping only stores one time step at a time, 
   but only enjoys a spatial data decomposition and spatial parallelism.  On the other
   hand, Warp stores multiple time steps simultaneously and each processor holds a space-time
   chunk reflecting both the spatial and temporal parallelism.

- Warp only handles temporal parallelism and is agnostic to the spatial decomposition.  
  See [warp_SplitCommworld](@ref warp_SplitCommworld).  
  Each processor owns a certain number of CF intervals of points, as depicted next, where
  each processor owns 2 CF intervals.
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.75\textwidth]{../img/parallel_timeline.pdf}}
       \label{img:data_layout}
   \end{figure}
   \endlatexonly
   Warp distributes Intervals evenly on the finest grid.

- Storage is greatly minimized by only storing C-points.  Whenever an F-point is needed,
  it is generated by F-relaxation.  That is, we only store the red C-point time values in the 
  previous figure.
   Coarsening can by aggressive with \f$ m = 8, 16, 32 \f$, so the storage requirements
   of Warp are significantly reduced when compared to storing all of the time values.

- In practice, storing only one space-time slab is advisable.  That is, solve
  for as many time steps (say k time steps) as you have available memory for.
  Then move on to the next k time steps.


## Cycling and relaxation strategies {#cyclingrelaxation}

There are two main cycling strategies available in Warp, F-and V-cycles.  These two
cycles differ in how often and the order in which coarse levels are visited.  A V-cycle
is depicted next, and is a simple recursive application of the @ref twogrid.
  \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.1\textwidth]{../img/VCycle.pdf}}
       \label{img:overlap}
   \end{figure}
   \endlatexonly

An F-cycle visits coarse grids more frequently and in a different order.  In essence,
this extra work gives you a closer approximation to a two-grid cycle, and a faster convergence
rate at the extra expense of more work.
  \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.27\textwidth]{../img/FCycle.pdf}}
       \label{img:overlap}
   \end{figure}
   \endlatexonly

Next, we make a few points about F- versus V-cycles.
- One V-cycle iteration is cheaper than one F-cycle iteration.
- But, F-cycles often converge more quickly.  For some test cases, this difference can be 
  quite large.  The cycle choice for the best time to solution will be problem
  dependent.  See @ref twodheat for a case study of cycling strategies. 

The number of FC relaxation sweeps is another important algorithmic setting.
Note that at least one F-relaxation sweep is always 
done on a level. A few summary points about relaxation are as follows.
- Using FCF (or even FCFCF, FCFCFCF) relaxation, corresponding to passing
  *warp_SetNRelax* a value of 1, 2 or 3 respectively, will result in a Warp cycle
  that converges more quickly as the number of relaxations grows.
- But as the number of relaxations grows, each Warp cycle becomes more expensive.  The optimal
  relaxation strategy for the best time to solution
  will be problem dependent.
- However, a good first step is to try FCF on all levels (i.e., *warp_SetNRelax(core, -1, 1)* ).
- A common optimization is to first set FCF on all levels (i.e., *warp_setnrelax(core, -1, 1)* ), 
  but then overwrite the FCF option on level 0 so that only F-relaxation is done on level 0, 
  (i.e., *warp_setnrelax(core, 0, 1)* ).  This strategy can work well with F-cycles. 
- See @ref twodheat  for a case study of relaxation strategies. 

Last, [Parallel Time Integration with Multigrid](https://computation-rnd.llnl.gov/linear_solvers/pubs/mgritPaper-2013.pdf)
has a more in depth case study of cycling and relaxation strategies

## Overlapping communication and computation {#overlapping}

  Warp effectively overlaps communication and computation.  The main
  computational kernel of Warp is relaxation (C or F).  At the start a each
  sweep, each processor first posts a send at its left-most point, and then
  carries out F-relaxation on its right-most interval in order to send the next
  processor the data that it needs.  If each processor has multiple intervals
  at this Warp level, this should allow for complete overlap.
  \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.35\textwidth]{../img/overlap.pdf}}
       \label{img:overlap}
   \end{figure}
   \endlatexonly


## Configuring the Warp Hierarchy {#config}

Some of the more basic Warp function calls allow you to control aspects 
discussed here.
- [warp_SetFMG](@ref warp_SetFMG): switches between using F- and V-cycles.
- [warp_SetMaxIter](@ref warp_SetMaxIter ):  sets the maximum number of Warp iterations
- [warp_SetCFactor](@ref warp_SetCFactor): sets the coarsening factor for any (or all levels)
- [warp_SetNRelax](@ref warp_SetNRelax): sets the number of CF-relaxation sweeps for any (or all levels)
- [warp_SetRelTol](@ref warp_SetRelTol), [warp_SetAbsTol](@ref warp_SetAbsTol): sets the stopping tolerance
- [warp_SetMaxCoarse](@ref warp_SetMaxCoarse): sets the maximum coarse grid size, in terms of C-points
- [warp_SetMaxLevels](@ref warp_SetMaxLevels): sets the maximum number of levels in the Warp hierarchy

## Heat equation example {#twodheat}

Here is some experimental data for the 2D heat equation, \f$ u_t = u_{xx} + u_{yy} \f$ 
generated by examples/drive-02.  
   \latexonly
   \begin{figure}[!ht] \centering 
       \subfloat{\includegraphics[width=0.45\textwidth]{../img/heat_results.pdf}}
       \label{img:heat_results}
   \end{figure}
   \endlatexonly
The problem setup is as follows.
- Backwards Euler is used as the time stepper.
- We used a Linux cluster with 4 cores per node, a Sandybridge Intel chipset, 
  and a fast Infiniband interconnect.  
- The space-time problem size was \f$ 129^2 \times 16,192 \f$ over the unit cube 
  \f$ [0,1] \times [0,1] \times [0,1] \f$ .  
- The coarsening factor was \f$ m = 16 \f$ on the finest level and \f$ m=2 \f$ 
  on coarser levels.  
- Since 16 processors optimized the serial time stepping approach, 16 processors
  in space are also used for the Warp experiments.  So for instance 512 processrs in 
  the plot corresponds to 16 processors in space and 32 processors in time, 
  \f$ 16*32 = 512 \f$.  Thus, each processor owns a space-time hypercube of
  \f$ (129^2 / 16) \times (16,192 / 32) \f$.  See @ref decomposition for a depiction
  of how Warp breaks the problem up.
- Various relaxation and V and F cycling strategies are experimented with.  
  + *V-cycle, FCF* denotes V-cycles and FCF-relaxation on each level.
  + *V-cycle, F-FCF* denotes V-cycles and F-relaxation on the finest level and
     FCF-relaxation on all coarser levels.
  + *F-cycle, F* denotes F-cycles and F-relaxation on each level.
- The initial guess at time values for \f$ t > 0 \f$ is zero, which is typical.

Regarding the performance, we can say
- The best speedup is 10x and this would grow if more processors were available.
- Although not shown, the iteration counts here are about 10-15 Warp iterations.
  See [Parallel Time Integration with Multigrid](https://computation-rnd.llnl.gov/linear_solvers/pubs/mgritPaper-2013.pdf)
  for the exact iteration counts.
- At smaller core counts, serial time stepping is faster.  But at about 256 processors,
  there is a crossover and Warp is faster.
- You can see the impact of the cycling and relaxation strategies discussed in
  @ref cyclingrelaxation.  For instance, even though *V-cycle, F-FCF* is a 
  weaker relaxation strategy than *V-cycle, FCF* (i.e., the Warp convergence 
  is slower), *V-cycle, F-FCF* has a faster time to solution than *V-cycle, FCF* 
  because each cycle is cheaper.
- In general, one level of aggressive coarsening (here by a factor 16) followed by
  slower coarsening was found to be best on this machine.

Achieving the best speedup can require some tuning, and it is recommended to read
[Parallel Time Integration with Multigrid](https://computation-rnd.llnl.gov/linear_solvers/pubs/mgritPaper-2013.pdf)
where this 2D heat equation example is explored in much more detail.

# Summary of Warp {#summary}

- Warp applies multigrid to the time dimension.
 + This exposes concurrency in the time dimension.
 + The potential for speedup is large, 10x, 100x, ...
- This is a non-intrusive approach, with an unchanged time discretization defined by user.
- Parallel time integration is only useful beyond some scale.  
  This is evidenced by the experimental results below.  For smaller numbers
  of cores sequential time stepping is faster, but at larger core counts
  Warp is much faster.
- The more time steps that you can parallelize over, the better your speedup will be.
- Warp is optimal for a variety of parabolic problems (see the examples directory).

