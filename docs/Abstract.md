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

It is **strongly recommended** that you also read [Parallel Time Integration
with Multigrid](https://computation-rnd.llnl.gov/linear_solvers/pubs/mgritPaper-2013.pdf)
after reading the @ref braidoverview.  It is a more in depth discussion of the algorithm
and associated experiments.

