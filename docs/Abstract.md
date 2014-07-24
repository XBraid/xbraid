This package implements an optimal-scaling multigrid solver for the linear
systems that arise from the discretization of problems with evolutionary
behavior. Typically, solution algorithms for evolution equations are based on a
time-marching approach, solving sequentially for one time step after the other.
Parallelism in these traditional time-integration techniques is limited to
spatial parallelism. However, current trends in computer architectures are
leading towards systems with more, but not faster, processors. Therefore,
faster compute speeds must come from greater parallelism. One approach to
achieve parallelism in time is with multigrid, but extending classical
multigrid methods for elliptic operators to this setting is a significant
achievement. In this software, we implement a non-intrusive, optimal-scaling
time-parallel method based on multigrid reduction techniques. The examples in
the package demonstrate optimality of our multigrid-reduction-in-time algorithm
(MGRIT) for solving a variety of equations in two and three spatial
dimensions. These examples can also be used to show that MGRIT can achieve
significant speedup in comparison to sequential time marching on modern
architectures.

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

Warp is designed to run in conjunction with an existing application code that
can be wrapped per our interface.  This application code will implement some
time marching type simulation like luid flow.  Essentially, the user has to
take their application code and extract a stand-alone time-stepping function
that can evolve a solution from one time value to another.  After this is done,
the Warp code takes care of the parallelism in the time dimension.

It is **strongly recommended** that you read [Parallel Time Integration
with Multigrid](https://computation-rnd.llnl.gov/linear_solvers/pubs/mgritPaper-2013.pdf)
before proceeding. 

